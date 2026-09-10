"""Geometric phase control with a single, fixed RF input phase."""
from __future__ import annotations
import math
import numpy as np


def field_nodes(cells, points_per_cell=16):
    """Sample a monotone Hermite spatial phase with pi radians per cell.

    k is the derivative of that phase, rather than unrelated endpoint values.
    Harmonic mean slopes avoid overshoot for gradually varying cell lengths.
    TRANSOPTR re-splines these samples, so its integrated phase must still be
    checked in the final run and by refining points_per_cell.
    """
    if points_per_cell < 4:
        raise ValueError('field_points_per_cell must be at least 4')
    lengths = np.array([c.L_cm for c in cells])
    if np.any(lengths <= 0):
        raise ValueError('Cell lengths must be positive')
    zedge = np.r_[0., np.cumsum(lengths)]
    secants = np.pi / lengths
    slopes = np.r_[secants[0], 2 / (1/secants[:-1]+1/secants[1:]), secants[-1]]
    zs, ks = [], []
    for i, length in enumerate(lengths):
        t = np.linspace(0, 1, points_per_cell, endpoint=False)
        # Derivative of cubic Hermite interpolant of phase (increments pi).
        k = 6*t*(1-t)*secants[i] + (3*t*t-4*t+1)*slopes[i] + (3*t*t-2*t)*slopes[i+1]
        zs.extend(zedge[i]+t*length)
        ks.extend(k)
    z = np.r_[zs, zedge[-1]]
    k = np.r_[ks, slopes[-1]]
    if len(z) > 9999 or np.any(k <= 0):
        raise ValueError('Invalid RFQ field sampling (k <= 0 or >9999 points)')
    entry_a = getattr(cells[0], 'entry_a_cm', None) or cells[0].a_cm
    a = np.interp(z, zedge, [entry_a]+[c.a_cm for c in cells])
    m = np.interp(z, zedge, [cells[0].m]+[c.m for c in cells])
    # Evaluate the single continuous opening law at every field node, not just
    # at three cell endpoints. This preserves zero slope/curvature at the ends.
    for cell in cells:
        if getattr(cell, 'section', '') == 'OM':
            start = cell.opening_start_cm
            end = start+cell.opening_length_cm
            mask = (z >= start) & (z <= end)
            u = np.clip((z[mask]-start)/cell.opening_length_cm,0,1)
            blend = u**3*(10+u*(-15+6*u))
            a[mask] = cell.opening_a_initial_cm + (cell.opening_a_final_cm-cell.opening_a_initial_cm)*blend
            m[mask] = 1.
            break
    a10 = (m*m-1)/(m*m*np.i0(k*a)+np.i0(k*m*a))
    a01 = (1-a10*np.i0(k*a))/(a*a)
    for i, cell in enumerate(cells):
        if getattr(cell, 'section', '') == 'FF':
            t = np.clip((z-zedge[i])/lengths[i],0,1)
            fade = 1-t*t*t*(10+t*(-15+6*t))
            a01 *= fade
            a10 *= fade
    return np.column_stack([z, a01, a10, k])


def phase_residual(path, rf_phase_deg, target_deg):
    values = np.loadtxt(path, ndmin=2)
    if values.shape[1] != 2 or not np.isfinite(values).all():
        raise ValueError('Invalid fort.20 phase output')
    # fort.20 contains spphase - phase, not spphase itself.
    actual = float(np.rad2deg(values[-1, 1]) + rf_phase_deg)
    # Do not wrap: crossing into a different RF bucket is not acceptance.
    return actual - target_deg, actual


def tune_cell(cells, config, cell_dir, run_feedback, update_geometry):
    control = config['phase_control']
    cell = cells[-1]
    original = cell.L_cm
    rf_phase = float(config['rf']['tracking_phase_deg'])
    max_trials = int(control.get('max_trials', 12))
    tolerance = float(control.get('tolerance_deg', 0.15))
    lower, upper = original*0.65, original*1.5
    last = None
    for attempt in range(max_trials):
        energy, status = run_feedback(cells, config, cell_dir)
        if status != 'optr_ok':
            raise RuntimeError(f'Phase trial failed: {status}; {cell_dir}')
        error, measured = phase_residual(cell_dir/'fort.20', rf_phase, cell.Phi_deg)
        if abs(error) <= tolerance:
            cell.phase_actual_deg = measured
            cell.phase_trials = attempt+1
            return energy, status
        length = cell.L_cm
        # Initial derivative: omega*dt/dL ~ pi/L. Secant corrects using real tracking.
        slope = 180/original
        if last is not None and abs(length-last[0]) > 1e-9:
            observed = (error-last[1])/(length-last[0])
            if observed > 0 and math.isfinite(observed):
                slope = observed
        step = np.clip(-error/slope, -0.15*original, 0.15*original)
        proposed = float(np.clip(length+step, lower, upper))
        last = (length, error)
        if abs(proposed-length) < 1e-9:
            break
        update_geometry(cell, proposed)
    raise RuntimeError(f'Phase did not converge for cell {cell.cell}: {error:.4f} deg; {cell_dir}')
