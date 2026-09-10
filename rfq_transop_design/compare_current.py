"""Generate mode-5 inputs and compare space charge on a frozen RFQ design."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
from .design_rfq import load_config, CellDesign, run_transoptr_feedback, read_fort_envelope
from .verify_phase import verify_final
from .beam_input import prepare_beam_input


def generate_input(lines, config, current_ma, emit_norm_mm_mrad):
    if not math.isfinite(current_ma) or current_ma<0 or not math.isfinite(emit_norm_mm_mrad) or emit_norm_mm_mrad<=0:
        raise ValueError('Current must be nonnegative and normalized RMS emittance positive')
    if int(lines[1].split()[1])!=5 or int(lines[5].split()[0])!=0 or int(lines[6].split()[0])!=11:
        raise ValueError('Expected uncoupled mode-5 MIRFQ 11-parameter template')
    scales=list(map(float,lines[4].split('!')[0].split()))
    beam=config['beam'];gamma=1+float(beam['energy_initial_mev'])/float(beam['mass_mev'])
    bg=math.sqrt((gamma-1)*(gamma+1));freq=float(config['rf']['frequency_hz'])
    if freq<=0:raise ValueError('Frequency must be positive')
    rms_geom=emit_norm_mm_mrad*1e-4/bg  # mm*mrad -> cm*rad
    native=5*rms_geom  # mode 5 uses sqrt(5)*RMS dimensions
    result=list(lines)
    first=result[0].split('!')[0].split()
    for index,key in [(0,'energy_initial_mev'),(3,'mass_mev'),(4,'charge_state')]:
        first[index]=str(float(beam[key]))
    charge=current_ma*1e-3/freq
    first[5]=f'{charge:.12g}'
    result[0]=' '.join(first)+' ! mode 5: charge per bunch [C], NOT current [A]'
    for index,offset in [(11,0),(14,2)]:
        values=result[index].split('!')[0].split()
        values[0]=f'{native*scales[offset]*scales[offset+1]:.12g}'
        result[index]=' '.join(values)+' ! native geometric emittance = 5 * RMS geometric'
    cfg=copy.deepcopy(config);cfg['beam']['initialization']='template_twiss'
    effective,report=prepare_beam_input(result,cfg)
    report.update({'average_current_ma':current_ma,'bunch_charge_coulomb':charge,
        'normalized_rms_emittance_mm_mrad_each_plane':emit_norm_mm_mrad,
        'rms_geometric_emittance_cm_rad':rms_geom,'native_geometric_emittance_cm_rad':native,
        'assumption':'one 90-degree native ellipsoid per RF period, all current in the bunch; not DC capture',
        'mode5_size_over_rms':math.sqrt(5)})
    return result,effective,report


def compare(config_path, output, current_ma=1., emittance=.25, prepare_only=False):
    cfg=load_config(Path(config_path));source=Path(cfg['output_dir']).resolve();out=Path(output).resolve()
    if out==source:raise ValueError('Use an output separate from the source design')
    out.mkdir(parents=True,exist_ok=True)
    lines=Path(cfg['transoptr']['template_data']).read_text().splitlines()
    override=Path(__file__).parent/'transoptr/RFQSC.f'
    cfg['transoptr']['rfq_source_override']=str(override.resolve())
    cfg['beam']['initialization']='template_twiss'
    cfg['beam']['longitudinal_full_width_deg']=90.
    cfg.setdefault('beam_diagnostics',{})['envelope_convention']='sqrt5_rms'
    history=pd.read_csv(source/'history.csv')
    cells=[CellDesign(**{k:None if pd.isna(v) else v for k,v in r.items()}) for r in history.to_dict('records')]
    records=[];profiles={};resolution_rows=[]
    for name,current in [('zero_current',0.),('with_current',current_ma)]:
        folder=out/name;folder.mkdir(exist_ok=True)
        template,effective,inputs=generate_input(lines,cfg,current,emittance)
        (folder/'template.dat').write_text('\n'.join(template)+'\n')
        (folder/'prepared_data.dat').write_text('\n'.join(effective)+'\n')
        (folder/'input_parameters.json').write_text(json.dumps(inputs,indent=2)+'\n')
        case=copy.deepcopy(cfg);case['output_dir']=str(folder)
        case['transoptr']['template_data']=str(folder/'template.dat')
        (folder/'config.yaml').write_text(yaml.safe_dump(case,sort_keys=False))
        shutil.copy2(source/'history.csv',folder/'history.csv')
        if prepare_only:continue
        print(f'Running {name}: {current:g} mA, Q={inputs["bunch_charge_coulomb"]:.6g} C',flush=True)
        validation=verify_final(cells,case,folder,run_transoptr_feedback,read_fort_envelope)
        points=validation['resolutions'][-1]['points_per_cell']
        final_folder=folder/'verification'/f'points_{points}'
        env=read_fort_envelope(final_folder/'fort.envelope').drop_duplicates('s')
        profiles[name]=env
        first=env.iloc[0];bg=math.sqrt((1+first.E/case['beam']['mass_mev'])**2-1)
        measured={axis:float(first[f'{axis}-emittance'])*bg/5*1e4 for axis in 'xy'}
        if not all(math.isclose(v,emittance,rel_tol=1e-4) for v in measured.values()):
            raise RuntimeError(f'Actual fort.envelope emittances do not match requested RMS input: {measured}')
        actual_charge=float((final_folder/'data.dat').read_text().splitlines()[0].split()[5])
        if not math.isclose(actual_charge,inputs['bunch_charge_coulomb'],rel_tol=1e-9,abs_tol=1e-25):
            raise RuntimeError('Actual TRANSOPTR bunch charge differs from request')
        optics=pd.read_csv(folder/'beam_optics.csv')
        env.to_csv(folder/'beam_envelope_full.csv',index=False)
        for r in validation['resolutions']:
            trial=read_fort_envelope(folder/'verification'/f'points_{r["points_per_cell"]}'/'fort.envelope').drop_duplicates('s')
            valid=(env.s>=trial.s.min())&(env.s<=trial.s.max())
            entry={'case':name,**r}
            for axis in 'xyz':
                ref=env.loc[valid,f'{axis}-envelope'].to_numpy()
                error=np.interp(env.loc[valid,'s'],trial.s,trial[f'{axis}-envelope'])-ref
                entry[f'{axis}_relative_envelope_error_vs_finest']=float(np.max(np.abs(error))/max(np.max(np.abs(ref)),1e-30))
            resolution_rows.append(entry)
        record={'case':name,'average_current_ma':current,'charge_coulomb':actual_charge,
            'energy_mev':float(env.E.iloc[-1]),'phase_and_refinement_passed':validation['phase_and_refinement_passed'],
            'energy_refinement_ev':validation['refinement_energy_difference_mev']*1e6,
            'max_phase_error_deg':max(r['max_phase_error_deg'] for r in validation['resolutions']),
            'exit_z_rms_mm':float(env['z-envelope'].iloc[-1])*10/math.sqrt(5),
            'exit_energy_spread_rms_kev':float(optics.energy_spread_linear_mev.iloc[-1])*1000/math.sqrt(5)}
        for axis in 'xy':
            record[f'input_{axis}_normalized_rms_mm_mrad']=measured[axis]
            record[f'max_{axis}_rms_mm']=float(env[f'{axis}-envelope'].max())*10/math.sqrt(5)
            record[f'max_{axis}_native_occupancy']=float(optics[f'{axis}_occupancy'].max())
        records.append(record)
    summary={'geometry_source':str(source),'geometry_reoptimized':False,
        'override_source':str(override.resolve()),'override_sha256':hashlib.sha256(override.read_bytes()).hexdigest(),
        'fix':'RFQSC uses CURRENT from COMMON/MOM/ instead of uninitialized BNCHARGE',
        'scope':'90-degree bunched mode-5 approximation; not continuous-beam capture or transmission prediction',
        'prepared_only':prepare_only,'cases':records}
    if not prepare_only:
        a,b=profiles['zero_current'],profiles['with_current']
        shared=(a.s>=b.s.min())&(a.s<=b.s.max())
        differences={axis:float(np.max(np.abs(np.interp(a.s[shared],b.s,b[f'{axis}-envelope'])-a.loc[shared,f'{axis}-envelope']))) for axis in 'xyz'}
        summary['max_current_effect_native_cm']=differences
        summary['nonzero_current_effect_detected']=max(differences.values())>1e-5
        pd.DataFrame(records).to_csv(out/'comparison.csv',index=False)
        pd.DataFrame(resolution_rows).to_csv(out/'refinement_comparison.csv',index=False)
        fig,axes=plt.subplots(5,1,figsize=(12,13),sharex=True)
        for name,env in profiles.items():
            optics=pd.read_csv(out/name/'beam_optics.csv')
            axes[0].plot(env.s,env.E,label=name)
            for ax,axis in zip(axes[1:4],'xyz'):ax.plot(env.s,env[f'{axis}-envelope']*10/math.sqrt(5),label=name)
            axes[4].plot(optics.z_cm,optics.energy_spread_linear_mev*1000/math.sqrt(5),label=name)
        for ax,label in zip(axes,['E [MeV]','x RMS [mm]','y RMS [mm]','z RMS [mm]','Dispersión E RMS [keV]']):
            ax.set_ylabel(label);ax.legend();ax.grid(alpha=.2)
        axes[-1].set_xlabel('z [cm]');fig.tight_layout();fig.savefig(out/'comparison.png',dpi=150);plt.close(fig)
    (out/'comparison.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('config')
    p.add_argument('--current-ma',type=float,default=1.);p.add_argument('--emittance-n-rms',type=float,default=.25)
    p.add_argument('--output',default='outputs/mirfq1_current_comparison');p.add_argument('--prepare-only',action='store_true')
    a=p.parse_args();compare(a.config,a.output,a.current_ma,a.emittance_n_rms,a.prepare_only)


if __name__=='__main__':main()
