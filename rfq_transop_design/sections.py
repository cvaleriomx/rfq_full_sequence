"""Configurable RFQ section schedules, inspired by FSP/NFSP (not optimizers)."""
from __future__ import annotations
import math


def smoothstep(t):
    t = min(1., max(0., float(t)))
    return t*t*t*(10+t*(-15+6*t))


def focusing_scale(config):
    """B = |q| V[MeV] lambda[cm]^2 A01[cm^-2] (rest-mass convention)."""
    return (abs(float(config['beam']['charge_state']))
            * float(config['rf']['vane_voltage_kv'])*.001
            / float(config['beam']['mass_mev'])
            * (29979245800./float(config['rf']['frequency_hz']))**2)


def aperture_for_B(B, m, k, config, coefficients):
    lo = float(config['geometry']['min_a_cm'])
    hi = float(config['geometry']['max_a_cm'])
    scale = focusing_scale(config)
    def value(a):
        return scale*coefficients(a, m, k)[0]
    if B <= 0 or not value(hi) <= B <= value(lo):
        raise ValueError(f'B={B:g} cannot be obtained within aperture limits')
    for _ in range(45):
        mid = (lo+hi)/2
        if value(mid) > B:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def stage_schedule(config):
    s = config['sections']
    b_rm = float(s.get('B_rm', 8.))
    b_peak = float(s.get('B_peak', 12.))
    b_acc = float(s.get('B_acc', 10.))
    phi_acc = float(s.get('phi_acc_deg', -30.))
    m_sh = float(s.get('m_sh', 1.08))
    m_acc = float(s.get('m_acc', 1.8))
    scheme = s.get('scheme', 'fsp').lower()
    if scheme == 'fsp':
        return [('SH', int(s.get('sh_cells', 20)), (b_rm, m_sh, -85.)),
                ('GB', int(s.get('gb_cells', 60)), (b_rm, m_acc, phi_acc)),
                ('ACC', None, (b_rm, m_acc, phi_acc))]
    if scheme == 'nfsp':
        return [('MS', int(s.get('ms_cells_after_rm', 20)), (b_peak, m_sh, -88.)),
                ('MB', int(s.get('mb_cells', 50)), (b_peak, 1.25, -45.)),
                ('MBA', int(s.get('mba_cells', 25)), (b_acc, m_acc, phi_acc)),
                ('MA', None, (b_acc, m_acc, phi_acc))]
    raise ValueError('sections.scheme must be fsp or nfsp')


def validate_sections(config):
    s = config['sections']
    n = s.get('rm_cells', 4)
    if int(n) != n or int(n) < 1:
        raise ValueError('rm_cells must be a positive integer (4 is a starting choice, not a universal limit)')
    for _, count, _ in stage_schedule(config):
        if count is not None and count < 1:
            raise ValueError('Section lengths must be positive cell counts')
    fixed = int(n)+sum(n for _, n, _ in stage_schedule(config) if n is not None)
    if int(config['iteration']['max_cells']) <= fixed:
        raise ValueError('max_cells must leave room for the acceleration section')
    if not config.get('phase_control', {}).get('enabled'):
        raise ValueError('Sectioned design requires the dense field and real phase controller')
    e = config.get('exit', {})
    if e.get('enabled', True):
        opening = e.get('opening_cells', 0)
        if int(opening) != opening or opening < 0:
            raise ValueError('opening_cells must be a nonnegative integer')
        if opening:
            final = exit_aperture(config)
            initial = math.sqrt(focusing_scale(config)/float(s.get('B_acc', 10.) if s.get('scheme', 'fsp') == 'nfsp' else s.get('B_rm', 8.)))
            if not initial < final <= float(config['geometry']['max_a_cm']):
                raise ValueError('Exit aperture must exceed the unmodulated aperture and respect max_a_cm')
        if int(e.get('transition_cells', 1)) < 1:
            raise ValueError('At least one transition cell is required before unmodulated exit')
        if float(e.get('unmodulated_length_cm', 2.)) <= 0 or float(e.get('fringe_length_cm', 1.)) <= 0:
            raise ValueError('Exit lengths must be positive')


def section_parameters(config, index, k, coefficients):
    s = config['sections']
    nrm = int(s.get('rm_cells', 4))
    entrance_a = float(s.get('rm_entrance_aperture_cm', 1.2))
    initial_B = focusing_scale(config)/entrance_a**2
    previous = (float(s.get('B_rm', 8.)), 1., -90.)
    if index <= nrm:
        # No modulation or acceleration in the RM. B increases via aperture.
        t = smoothstep(index/nrm)
        B = initial_B + t*(previous[0]-initial_B)
        values, section = (B, 1., -90.), 'RM'
    else:
        local = index-nrm
        for name, count, end in stage_schedule(config):
            if count is None:
                values, section = end, name
                break
            if local <= count:
                t = smoothstep(local/count)
                values = tuple(a+(b-a)*t for a,b in zip(previous,end))
                section = name
                break
            local -= count
            previous = end
    B, m, phi = values
    a = aperture_for_B(B, m, k, config, coefficients)
    return a, m, B, phi, section, entrance_a


def exit_count(config):
    e = config.get('exit', {})
    return int(e.get('transition_cells', 1)) + (int(e.get('opening_cells', 0)) or 1) + 1


def exit_aperture(config):
    value = config.get('exit', {}).get('final_aperture_cm', 'rm_entrance')
    if value == 'rm_entrance':
        return float(config['sections'].get('rm_entrance_aperture_cm', 1.2))
    return float(value)


def exit_parameters(config, ordinal, last_acc, energy, kinematics, coefficients, opening_length_cm=None):
    e = config.get('exit', {})
    ntc = int(e.get('transition_cells', 1))
    opening = int(e.get('opening_cells', 0))
    if opening and ordinal > ntc:
        initial = math.sqrt(focusing_scale(config)/last_acc.B)
        final = exit_aperture(config)
        if final <= initial:
            raise ValueError('The requested exit aperture must grow beyond the TC aperture')
        local = ordinal-ntc
        if local <= opening:
            length = (opening_length_cm or opening*kinematics.cell_length_cm)/opening
            a = initial+(final-initial)*smoothstep(local/opening)
            section = 'OM'
        else:
            length = float(e.get('fringe_length_cm', 1.))
            a, section = final, 'FF'
        return a, 1., focusing_scale(config)/a**2, last_acc.Phi_deg, length, section
    if ordinal <= ntc:
        t = smoothstep(ordinal/ntc)
        m = last_acc.m+(1-last_acc.m)*t
        length = kinematics.cell_length_cm*float(e.get('transition_length_factor', .5))
        section = 'TC'
    else:
        m = 1.
        section = 'UM' if ordinal == ntc+1 else 'FF'
        length = float(e.get('unmodulated_length_cm', 2.) if section == 'UM' else e.get('fringe_length_cm', 1.))
    k = math.pi/length
    a = aperture_for_B(last_acc.B, m, k, config, coefficients)
    return a, m, last_acc.B, last_acc.Phi_deg, length, section
