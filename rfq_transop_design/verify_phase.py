"""Verify one complete RFQ, including refinement of TRANSOPTR field samples."""
from __future__ import annotations
import copy
import json
from pathlib import Path
import numpy as np
import pandas as pd


def verify_final(cells, config, output_dir, run_feedback, read_envelope):
    output_dir = Path(output_dir)
    base_points = int(config['phase_control'].get('field_points_per_cell', 16))
    rf_phase = float(config['rf']['tracking_phase_deg'])
    comparisons = []
    factors = [f for f in (1,2,4) if len(cells)*base_points*f+1 <= 9999]
    if len(factors) < 2:
        raise ValueError('At least two field resolutions are required for verification')
    phase_mask = np.array([c.phase_match for c in cells], dtype=bool)
    for factor in factors:
        cfg = copy.deepcopy(config)
        points = base_points*factor
        if len(cells)*points+1 > 9999:
            raise ValueError('Verification exceeds TRANSOPTR 9999 field nodes')
        cfg['phase_control']['field_points_per_cell'] = points
        folder = output_dir/'verification'/f'points_{points}'
        folder.mkdir(parents=True, exist_ok=True)
        # Tighten RK tolerance at the finest resolution as an additional check.
        if factor == factors[-1]:
            template = Path(cfg['transoptr']['template_data']).read_text().splitlines()
            parts = template[1].split('!')[0].split()
            parts[3] = str(float(parts[3])*0.1)
            template[1] = ' '.join(parts)
            input_path = folder/'verification_input.dat'
            input_path.write_text('\n'.join(template)+'\n')
            cfg['transoptr']['template_data'] = str(input_path.resolve())
        energy, status = run_feedback(cells, cfg, folder)
        if status != 'optr_ok':
            raise RuntimeError(f'Final verification failed: {status}')
        env = read_envelope(folder/'fort.envelope')
        phase = np.loadtxt(folder/'fort.20', ndmin=2)
        if not np.isfinite(env[['s','E']].to_numpy()).all():
            raise ValueError('Non-finite final energy profile')
        if abs(float(env.s.iloc[-1])-cells[-1].Z_cm) > 0.02:
            raise ValueError('TRANSOPTR did not reach RFQ exit')
        z = np.array([c.Z_cm for c in cells])
        actual = np.interp(z, phase[:,0], phase[:,1])*180/np.pi + rf_phase
        error = actual-np.array([c.Phi_deg for c in cells])
        energies = np.interp(z, env.s, env.E)
        # RFQ repurposes a centroid for spatial phase, then rescales it at
        # the exit momentum change. Recover phase from ct and fort.20 instead.
        temporal = np.interp(z, env.s, env['ct'])*2*np.pi*float(config['rf']['frequency_hz'])/29979245800.
        spatial = temporal-np.interp(z, phase[:,0], phase[:,1])
        spatial_error = spatial - np.pi*np.arange(1,len(cells)+1)
        comparisons.append({'points_per_cell': points, 'energy_mev': energy,
            'max_phase_error_deg': float(np.max(abs(error[phase_mask]))),
            'max_spatial_phase_error_rad': float(np.max(abs(spatial_error))),
            'exit_position_error_cm': float(env.s.iloc[-1])-cells[-1].Z_cm})
        if factor == 1:
            pd.DataFrame({'cell': [c.cell for c in cells], 'z_cm': z,
                'section': [c.section for c in cells], 'phase_controlled': phase_mask,
                'phase_target_deg': [c.Phi_deg for c in cells],
                'phase_actual_deg': actual, 'phase_error_deg': error,
                'energy_final_run_mev': energies,
                'spatial_phase_error_rad': spatial_error,
            }).to_csv(output_dir/'phase_validation.csv', index=False)
    pipeline_path = output_dir/'table_pipeline.txt'
    if pipeline_path.exists():
        pipeline = pd.read_csv(pipeline_path, sep=r"\s+")
        profile = pd.read_csv(output_dir/'phase_validation.csv')
        if len(pipeline) != len(profile)+1:
            raise ValueError('Pipeline table / final profile length mismatch')
        pipeline.loc[pipeline.index[1:], 'Wsyn'] = profile.energy_final_run_mev.to_numpy()
        pipeline.to_csv(pipeline_path, sep=' ', index=False)
    delta = abs(comparisons[-1]['energy_mev']-comparisons[-2]['energy_mev'])
    max_phase = max(r['max_phase_error_deg'] for r in comparisons)
    report = {'rf_phase_deg': rf_phase, 'cells': len(cells),
        'length_m': cells[-1].Z_cm/100, 'resolutions': comparisons,
        'energy_target_mev': config['beam']['energy_target_mev'],
        'target_reached': abs(comparisons[-1]['energy_mev']-float(config['beam']['energy_target_mev'])) <= float(config['beam']['energy_tolerance_mev']),
        'refinement_energy_difference_mev': delta,
        'phase_and_refinement_passed': max_phase <= .5 and delta <= .001,
        'scope': 'Reference-particle synchronism and numerical convergence; not beam transmission or vane engineering validation.'}
    if config.get('sections', {}).get('enabled'):
        report['section_counts'] = {label: sum(c.section==label for c in cells) for label in dict.fromkeys(c.section for c in cells)}
        report['exit_included_in_energy_check'] = True
        report['exit_model'] = 'Two-term transition, unmodulated constant/opening vane, approximate fringe envelope; not a 3D end-field solution.'
        openings = [c for c in cells if c.section == 'OM']
        if openings:
            report['exit_opening'] = {'segments':len(openings), 'profile':'quintic smoothstep',
                'length_cm':openings[0].opening_length_cm,
                'initial_radius_cm':openings[0].opening_a_initial_cm,
                'final_radius_cm':openings[0].opening_a_final_cm}
        report['matching_validated'] = False
    (output_dir/'phase_validation.json').write_text(json.dumps(report, indent=2)+'\n')
    import matplotlib.pyplot as plt
    table = pd.read_csv(output_dir/'phase_validation.csv')
    fig, axes = plt.subplots(2,1,figsize=(9,7),sharex=True)
    selected = table.phase_controlled.astype(bool)
    axes[0].plot(table.z_cm[selected]/100,table.phase_target_deg[selected],label='Objetivo')
    axes[0].plot(table.z_cm[selected]/100,table.phase_actual_deg[selected],'--',label='TRANSOPTR, corrida completa')
    axes[0].set_ylabel('Fase relativa [grados]'); axes[0].legend()
    axes[1].plot(table.z_cm/100,table.energy_final_run_mev)
    axes[1].axhline(config['beam']['energy_target_mev'],color='gray',linestyle='--')
    axes[1].set_ylabel('Energía [MeV]'); axes[1].set_xlabel('z [m]')
    fig.tight_layout(); fig.savefig(output_dir/'phase_validation.png'); plt.close(fig)
    if config.get('sections', {}).get('enabled'):
        try:
            from .section_diagnostics import write_section_diagnostics
        except ImportError:
            from section_diagnostics import write_section_diagnostics
        write_section_diagnostics(cells, config, output_dir)
        try:
            from .beam_diagnostics import write_beam_diagnostics
        except ImportError:
            from beam_diagnostics import write_beam_diagnostics
        write_beam_diagnostics(config, output_dir, folder)
    input_report = folder/'beam_input.json'
    if input_report.exists():
        (output_dir/'beam_input.json').write_text(input_report.read_text())
    return report
