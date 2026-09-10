"""Postprocess real TRANSOPTR envelope optics; never declare beam acceptance."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_twiss(path):
    """VECTIVE.f format 555: element, label, s, alpha, beta, emit, D, D', turns."""
    rows = []
    if path.exists():
        for line in path.read_text().splitlines():
            fields = line.replace('D', 'E').split()
            try:
                values = [float(v) for v in fields[-7:]]
            except ValueError:
                continue
            if len(values) == 7:
                rows.append(values)
    return pd.DataFrame(rows, columns=['z_cm','alpha','beta_cm','emit_cm_rad',
                                       'disp_cm','disp_prime','mu_turns']).drop_duplicates('z_cm')


def interpolate(z, xp, values):
    return np.interp(z, xp, values, left=np.nan, right=np.nan)


def write_beam_diagnostics(config, output_dir, source=None):
    try:
        from .design_rfq import read_fort_envelope
    except ImportError:
        from design_rfq import read_fort_envelope
    output_dir = Path(output_dir)
    if source is None:
        candidates = list((output_dir/'verification').glob('points_*/fort.envelope'))
        if not candidates:
            raise FileNotFoundError('No verified TRANSOPTR envelope; run the real design first')
        source = max(candidates, key=lambda p: int(p.parent.name.split('_')[-1])).parent
    source = Path(source)
    env = read_fort_envelope(source/'fort.envelope').drop_duplicates('s').sort_values('s')
    geo = pd.read_csv(output_dir/'vane_geometry.csv')
    sections = json.loads((output_dir/'sections.json').read_text())['sections']
    history = pd.read_csv(output_dir/'history.csv')
    settings = config.get('beam_diagnostics', {})
    multiplier = float(settings.get('envelope_multiplier', 1.))
    if not np.isfinite(multiplier) or multiplier <= 0:
        raise ValueError('beam_diagnostics.envelope_multiplier must be positive')
    z = env.s.to_numpy()
    out = pd.DataFrame({'z_cm':z})
    edges = np.array([s['z_end_cm'] for s in sections])
    out['section'] = np.array([s['section'] for s in sections])[np.minimum(np.searchsorted(edges,z),len(edges)-1)]
    warnings = []
    report = {'source':str(source.resolve()), 'envelope_multiplier':multiplier,
        'envelope_convention':settings.get('envelope_convention','unconfirmed'),
        'centroid_assumption':'centered beam; RFQ repurposes centroid slots for reference phase',
        'phase_definition':'Native VECTIVE fort.81/83 phase in turns, converted to degrees; transported optics, not periodic eigenphase or synchronous RF phase',
        'period_definition':'nonoverlapping pairs of regular cells, starting at the entrance; pairs crossing section boundaries excluded; TC/OM/UM/FF excluded',
        'mismatch_status':'not evaluated: no matched reference supplied',
        'tune_depression_status':'not evaluated: requires equivalent finite-current and zero-current optics',
        'beam_acceptance_validated':False, 'warnings':warnings}
    if report['envelope_convention'] != 'rms':
        warnings.append('Envelope RMS convention is unconfirmed; multiplier acts on native TRANSOPTR sizes, not a guaranteed particle percentile.')
    beta_rel = env['beta'].to_numpy()
    gamma = 1 + env.E.to_numpy()/float(config['beam']['mass_mev'])
    valid_phase = {}
    for axis, corr, unit in [('x','r12',81),('y','r34',83)]:
        size = env[f'{axis}-envelope'].abs().to_numpy()
        divergence = env[f"{axis}'-envelope"].abs().to_numpy()
        rho = env[corr].to_numpy()
        emit = env[f'{axis}-emittance'].to_numpy()
        out[f'{axis}_envelope_cm'] = size
        out[f'{axis}_divergence_mrad'] = 1000*divergence
        out[f'{axis}_correlation'] = rho
        out[f'{axis}_emit_normalized_cm_rad'] = emit*beta_rel*gamma
        good = np.isfinite(emit) & (emit > 0) & (np.abs(rho) < 1)
        out[f'{axis}_alpha_projected'] = np.divide(-rho*size*divergence,emit,out=np.full(len(z),np.nan),where=good)
        out[f'{axis}_beta_projected_cm'] = np.divide(size**2,emit,out=np.full(len(z),np.nan),where=good)
        # Independent of Twiss/emit validity: local envelope derivative from covariance.
        out[f'{axis}_envelope_slope'] = rho*divergence
        metal = geo.metal_profile.astype(bool)
        g = geo[metal]
        aperture = interpolate(z,g.z_cm,g[f'{axis}_plus_cm'])
        out[f'{axis}_aperture_cm'] = aperture
        out[f'{axis}_occupancy'] = multiplier*size/aperture
        twiss = read_twiss(source/f'fort.{unit}')
        valid = bool(good[0]) and len(twiss) > 1
        valid_phase[axis] = valid
        if valid:
            if twiss.z_cm.iloc[0] > 0 and z[0] == 0:
                twiss = pd.concat([pd.DataFrame([{'z_cm':0.,'mu_turns':0.}]),twiss],ignore_index=True)
            out[f'{axis}_mu_deg'] = interpolate(z,twiss.z_cm,twiss.mu_turns*360)
            if np.any(np.diff(twiss.mu_turns) < -1e-5) or np.any(np.diff(twiss.mu_turns) >= .5):
                warnings.append(f'{axis}: native phase is nonmonotonic or insufficiently sampled; period advances disabled.')
                valid_phase[axis] = False
        else:
            out[f'{axis}_mu_deg'] = np.nan
            out[f'{axis}_alpha_projected'] = np.nan
            out[f'{axis}_beta_projected_cm'] = np.nan
            warnings.append(f'{axis}: zero/invalid initial emittance or missing Twiss output; Twiss/phase unavailable. Check input beam.')
        occupancy = out[f'{axis}_occupancy'].to_numpy()
        imax = int(np.nanargmax(occupancy))
        report[axis] = {'initial_emit_cm_rad':float(emit[0]),
            'phase_usable':valid_phase[axis], 'max_envelope_cm':float(size.max()),
            'max_envelope_z_cm':float(z[np.argmax(size)]),
            'max_divergence_mrad':float(1000*divergence.max()),
            'max_occupancy':float(occupancy[imax]),'max_occupancy_z_cm':float(z[imax]),
            'envelope_exceeds_aperture':bool(np.nanmax(occupancy)>=1)}
        if report[axis]['envelope_exceeds_aperture']:
            warnings.append(f'{axis}: selected envelope exceeds the ideal vane aperture; envelope propagation does not remove lost particles.')
    # Longitudinal coordinates: z is bunch length, z-prime is dp/p despite its generic rad label.
    out['longitudinal_envelope_cm'] = env['z-envelope'].abs().to_numpy()
    out['momentum_spread_native'] = env["z'-envelope"].abs().to_numpy()
    out['phase_width_linear_deg'] = 360*float(config['rf']['frequency_hz'])*out.longitudinal_envelope_cm/(beta_rel*29979245800.)
    out['energy_spread_linear_mev'] = beta_rel**2*gamma*float(config['beam']['mass_mev'])*out.momentum_spread_native
    out['longitudinal_correlation'] = env.r56.to_numpy()
    if out.phase_width_linear_deg.max() > 180:
        warnings.append('Longitudinal envelope exceeds 180 RF degrees: local linear bunch optics cannot establish RF capture.')
    reference = settings.get('matched_reference_csv')
    for axis in ('x','y'):
        out[f'{axis}_mismatch_bmag'] = np.nan
    if reference:
        ref_path = Path(reference).expanduser().resolve()
        ref = pd.read_csv(ref_path).sort_values('z_cm')
        if len(ref) < 2 or ref.z_cm.duplicated().any():
            raise ValueError('Matched reference needs at least two distinct z_cm positions')
        for axis in ('x','y'):
            ar,br = ref[f'{axis}_alpha'],ref[f'{axis}_beta_cm']
            if not np.isfinite(ref.z_cm).all() or not np.isfinite(ar).all() or not np.isfinite(br).all() or (br<=0).any():
                raise ValueError('Matched reference Twiss must be finite with positive beta')
            alpha_ref=interpolate(z,ref.z_cm,ar)
            beta_ref=interpolate(z,ref.z_cm,br)
            alpha=out[f'{axis}_alpha_projected']
            beta=out[f'{axis}_beta_projected_cm']
            out[f'{axis}_mismatch_bmag']=.5*(beta/beta_ref+beta_ref/beta+
                (alpha*beta_ref-alpha_ref*beta)**2/(beta*beta_ref))
        report['mismatch_status']='Bmag against user supplied projected Twiss reference; 1 is matched, uncoupled interpretation only'
        report['matched_reference_csv']=str(ref_path)
    out.to_csv(output_dir/'beam_optics.csv',index=False)
    periods = []
    starts = np.r_[0.,history.Z_cm.to_numpy()[:-1]]
    for i in range(0,len(history)-1,2):
        pair = history.iloc[i:i+2]
        if pair.section.nunique()!=1 or pair.section.iloc[0] in {'TC','OM','UM','FF'}:
            continue
        a,b = float(starts[i]),float(pair.Z_cm.iloc[-1])
        row = {'first_cell':int(pair.cell.iloc[0]),'last_cell':int(pair.cell.iloc[-1]),
               'section':pair.section.iloc[0],'z_start_cm':a,'z_end_cm':b}
        for axis in ('x','y'):
            values = interpolate([a,b],z,out[f'{axis}_mu_deg'])
            row[f'{axis}_advance_deg'] = float(values[1]-values[0]) if valid_phase[axis] else np.nan
        periods.append(row)
    periods = pd.DataFrame(periods,columns=['first_cell','last_cell','section','z_start_cm','z_end_cm','x_advance_deg','y_advance_deg'])
    periods.to_csv(output_dir/'beam_phase_advance.csv',index=False)
    for axis in ('x','y'):
        advances = periods[f'{axis}_advance_deg'].dropna()
        report[axis]['max_pair_advance_deg'] = float(advances.max()) if len(advances) else None
        report[axis]['min_pair_advance_deg'] = float(advances.min()) if len(advances) else None
        report[axis]['configured_limit_exceedances'] = []
        for setting,metric in [('max_occupancy','max_occupancy'),
                               ('max_divergence_mrad','max_divergence_mrad'),
                               ('max_pair_advance_deg','max_pair_advance_deg')]:
            limit = settings.get(setting)
            if limit is not None:
                if not np.isfinite(float(limit)) or float(limit) <= 0:
                    raise ValueError(f'beam_diagnostics.{setting} must be positive')
                value = report[axis][metric]
                if value is not None and value > float(limit):
                    report[axis]['configured_limit_exceedances'].append(setting)
                    warnings.append(f'{axis}: {metric}={value:g} exceeds configured limit {limit}.')
    report['configured_limits'] = {key:settings.get(key) for key in
        ('max_occupancy','max_divergence_mrad','max_pair_advance_deg')}
    summary = []
    for section in sections:
        rows = out[out.section==section['section']]
        item = {'section':section['section'],'samples':len(rows)}
        for axis in ('x','y'):
            for quantity in ('envelope_cm','divergence_mrad','occupancy'):
                value = rows[f'{axis}_{quantity}'].max()
                item[f'max_{axis}_{quantity}'] = float(value) if pd.notna(value) else None
        summary.append(item)
    report['sections'] = summary
    (output_dir/'beam_optics_summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    fig,axes = plt.subplots(8,1,figsize=(13,20),sharex=True)
    for axis in ('x','y'):
        axes[0].plot(periods.z_end_cm,periods[f'{axis}_advance_deg'],'.-',label=axis)
        axes[1].plot(z,out[f'{axis}_occupancy'],label=axis)
        axes[2].plot(z,out[f'{axis}_divergence_mrad'],label=axis)
        axes[3].plot(z,out[f'{axis}_alpha_projected'],label=axis)
    axes[1].axhline(1,color='red',ls='--',label='Apertura')
    axes[4].plot(z,out.phase_width_linear_deg,label='Anchura de fase lineal')
    axes[5].plot(z,out.energy_spread_linear_mev*1000,label='Dispersión energética lineal')
    for axis in ('x','y'):
        axes[6].plot(z,out[f'{axis}_envelope_cm']*10,label=axis)
        axes[7].plot(z,out[f'{axis}_mismatch_bmag'],label=axis)
    if not reference:
        axes[7].text(.5,.5,'Sin referencia adaptada: no evaluado',ha='center',transform=axes[7].transAxes)
    for ax,label in zip(axes,['Avance / par [grados]','Envolvente / apertura','Divergencia [mrad]','Alpha proyectado','Anchura de fase [grados]','Dispersión E [keV]','Envolvente nativa [mm]','Mismatch Bmag']):
        ax.set_ylabel(label);ax.legend();ax.grid(alpha=.2)
        for j,s in enumerate(sections):
            ax.axvspan(s['z_start_cm'],s['z_end_cm'],color=f'C{j%10}',alpha=.07)
    for s in sections:
        if s['section'] not in {'TC','OM','UM','FF'}:
            axes[0].text((s['z_start_cm']+s['z_end_cm'])/2,1.02,s['section'],ha='center',transform=axes[0].get_xaxis_transform())
    axes[-1].set_xlabel('z [cm]')
    fig.suptitle('Óptica de envolventes — consultar advertencias en beam_optics_summary.json')
    fig.tight_layout(rect=(0,0,1,.97));fig.savefig(output_dir/'beam_optics.png',dpi=140);plt.close(fig)
    return report


def main():
    from .design_rfq import load_config
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config');args=parser.parse_args()
    config=load_config(Path(args.config))
    report=write_beam_diagnostics(config,config['output_dir'])
    print(json.dumps(report,indent=2,allow_nan=False))


if __name__ == '__main__':
    main()
