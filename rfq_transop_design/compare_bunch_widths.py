"""Compare full longitudinal phase widths on one frozen RFQ geometry."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
from .design_rfq import CellDesign, load_config, read_fort_envelope, run_transoptr_feedback
from .verify_phase import verify_final
from .beam_input import half_length_cm


def compare(config_path, destination, widths=(180.,90.)):
    config=load_config(Path(config_path))
    source=Path(config['output_dir']).resolve()
    destination=Path(destination).resolve()
    if destination == source:
        raise ValueError('Comparison output must differ from the original design')
    destination.mkdir(parents=True,exist_ok=True)
    history=pd.read_csv(source/'history.csv')
    cells=[CellDesign(**{k:None if pd.isna(v) else v for k,v in row.items()})
           for row in history.to_dict('records')]
    # Use the input actually used in the existing verification, not a subsequently edited template.
    candidates=list((source/'verification').glob('points_*/data.dat'))
    original=max(candidates,key=lambda p:int(p.parent.name.split('_')[-1]))
    lines=original.read_text().splitlines()
    # Restore baseline RK tolerance: the finest verification template already tightens it.
    base_folder=source/'verification'/f"points_{config['phase_control']['field_points_per_cell']}"
    lines=(base_folder/'data.dat').read_text().splitlines()
    if [float(x) for x in lines[4].split('!')[0].split()][:6] != [1.]*6:
        raise ValueError('Comparison currently requires native input units cm and rad (CONX=1)')
    initial=[float(x) for x in lines[0].split('!')[0].split()]
    if int(lines[1].split()[1]) != 5:
        raise ValueError('RFQ comparison requires TRANSOPTR mode 5')
    for key,index in [('energy_initial_mev',0),('mass_mev',3),('charge_state',4)]:
        config['beam'][key]=initial[index]
    row=lines[3].split('!')[0].split()
    original_half=float(row[4])
    cm_per_degree=half_length_cm(1.,initial[0],initial[3],config['rf']['frequency_hz'])
    cases=[('baseline',original_half/cm_per_degree,original_half)]
    cases += [(f'width_{width:g}deg',width,half_length_cm(width,initial[0],initial[3],config['rf']['frequency_hz'])) for width in widths]
    summaries=[];refinements=[];profiles={}
    for name,width,half in cases:
        print(f'Running {name}: full width={width:.6g} deg, half length={half:.9g} cm',flush=True)
        folder=destination/name;folder.mkdir(exist_ok=True)
        deck=list(lines);dims=list(row);dims[4]=f'{half:.12g}'
        deck[3]=' '.join(dims)+f' ! full longitudinal width {width:.9g} deg; other dimensions unchanged'
        input_path=folder/'input.dat';input_path.write_text('\n'.join(deck)+'\n')
        cfg=copy.deepcopy(config);cfg['output_dir']=str(folder)
        cfg['beam']['initialization']='dimensions'  # replay the saved beam, not new template Twiss
        cfg['beam']['longitudinal_full_width_deg']=width if width<=360 else None
        cfg['transoptr']['template_data']=str(input_path)
        (folder/'config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
        for filename in ('history.csv','table_pipeline.txt'):
            if (source/filename).exists():shutil.copy2(source/filename,folder/filename)
        report=verify_final(cells,cfg,folder,run_transoptr_feedback,read_fort_envelope)
        resolutions=report['resolutions'];finest=resolutions[-1]['points_per_cell']
        fine=read_fort_envelope(folder/'verification'/f'points_{finest}'/'fort.envelope').drop_duplicates('s')
        profiles[name]=fine
        fine.to_csv(folder/'beam_envelope_full.csv',index=False)
        for resolution in resolutions:
            points=resolution['points_per_cell']
            env=read_fort_envelope(folder/'verification'/f'points_{points}'/'fort.envelope').drop_duplicates('s')
            shared=(fine.s>=env.s.min())&(fine.s<=env.s.max())
            ref=fine.loc[shared]
            record={'case':name,**resolution}
            for coordinate in ('x','y','z'):
                key=f'{coordinate}-envelope'
                diff=np.interp(ref.s,env.s,env[key])-ref[key].to_numpy()
                record[f'{coordinate}_envelope_max_difference_cm']=float(np.max(np.abs(diff)))
                record[f'{coordinate}_envelope_relative_max_difference']=float(np.max(np.abs(diff))/max(np.max(np.abs(ref[key])),1e-30))
            refinements.append(record)
        envelope= json.loads((folder/'beam_optics_summary.json').read_text())
        summaries.append({'case':name,'full_width_deg':width,'half_length_cm':half,
            'final_energy_mev':resolutions[-1]['energy_mev'],
            'max_phase_error_deg':max(r['max_phase_error_deg'] for r in resolutions),
            'refinement_energy_difference_ev':report['refinement_energy_difference_mev']*1e6,
            'phase_and_refinement_passed':report['phase_and_refinement_passed'],
            'target_reached':report['target_reached'],
            'exit_longitudinal_envelope_cm':float(fine['z-envelope'].iloc[-1]),
            'max_x_envelope_cm':envelope['x']['max_envelope_cm'],
            'max_y_envelope_cm':envelope['y']['max_envelope_cm']})
    table=pd.DataFrame(summaries);table.to_csv(destination/'comparison.csv',index=False)
    pd.DataFrame(refinements).to_csv(destination/'refinement_comparison.csv',index=False)
    metadata={'geometry_source':str(source),'geometry_reoptimized':False,
        'input_charge_coulomb_mode5':initial[5],
        'width_definition':'Full width of the native longitudinal ellipsoid, not RMS or FWHM',
        'fixed_inputs':'RFQ geometry, RF phase, transverse input, momentum spread and bunch charge',
        'scope':'Numerical convergence and linear envelope sensitivity; not capture of a continuous beam or a Hofmann stability validation',
        'transverse_input_warning':'Original input may have zero emittance; see per-case beam_optics_summary.json',
        'cases':summaries}
    (destination/'comparison.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    fig,axs=plt.subplots(4,1,figsize=(12,12),sharex=True)
    for name,env in profiles.items():
        for ax,column in zip(axs,['E','x-envelope','y-envelope','z-envelope']):
            ax.plot(env.s,env[column],label=name)
    for ax,label in zip(axs,['Energía [MeV]','Envolvente x [cm]','Envolvente y [cm]','Envolvente longitudinal [cm]']):
        ax.set_ylabel(label);ax.legend();ax.grid(alpha=.2)
    axs[-1].set_yscale('log');axs[-1].set_xlabel('z [cm]')
    fig.suptitle('Anchos completos de entrada; misma RFQ y carga de entrada')
    fig.tight_layout();fig.savefig(destination/'comparison.png',dpi=150);plt.close(fig)
    print(table.to_string(index=False))
    return metadata


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config');parser.add_argument('--output',default='outputs/mirfq1_width_comparison')
    args=parser.parse_args();compare(args.config,args.output)


if __name__=='__main__':main()
