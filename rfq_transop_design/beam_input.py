"""Convert the full longitudinal phase width to native TRANSOPTR input units."""
import math


def half_length_cm(full_width_deg, energy_mev, mass_mev, frequency_hz):
    values=(full_width_deg,energy_mev,mass_mev,frequency_hz)
    if not all(math.isfinite(v) and v > 0 for v in values) or full_width_deg > 360:
        raise ValueError('Width must be in (0, 360] degrees; energy, mass and frequency must be positive')
    gamma=1+energy_mev/mass_mev
    beta=math.sqrt((gamma-1)*(gamma+1))/gamma
    return beta*29979245800./frequency_hz*full_width_deg/720.


def apply_longitudinal_width(lines, config):
    width=config['beam'].get('longitudinal_full_width_deg')
    if width is None:
        return lines
    result=list(lines)
    dimensions=result[3].split('!')[0].split()
    units=result[4].split('!')[0].split()
    if len(dimensions)<6 or len(units)<6:
        raise ValueError('Invalid TRANSOPTR dimensions/units input')
    scale=float(units[4])
    if not math.isfinite(scale) or scale<=0:
        raise ValueError('Invalid longitudinal unit scale')
    half=half_length_cm(float(width),float(config['beam']['energy_initial_mev']),
                        float(config['beam']['mass_mev']),float(config['rf']['frequency_hz']))
    dimensions[4]=f'{half*scale:.12g}'
    result[3]=' '.join(dimensions)+f' ! full longitudinal width {float(width):g} deg (not RMS/FWHM)'
    return result


def prepare_beam_input(lines, config):
    """Materialize native uncoupled Twiss as CIC3 does, before TRANSOPTR starts."""
    mode=config['beam'].get('initialization','dimensions')
    if mode=='dimensions':
        result=apply_longitudinal_width(lines,config)
        return result, {'initialization':mode,'longitudinal_full_width_deg':config['beam'].get('longitudinal_full_width_deg')}
    if mode!='template_twiss':
        raise ValueError('beam.initialization must be dimensions or template_twiss')
    def numbers(line):
        return line.split('!')[0].replace('D','E').replace('d','e').split()
    units=list(map(float,numbers(lines[4])))
    if len(units)!=8 or not all(math.isfinite(u) and u>0 for u in units):
        raise ValueError('Twiss initialization requires eight positive CONX scales')
    nc=int(numbers(lines[5])[0])
    if nc!=0:
        raise ValueError('template_twiss requires an uncoupled template with zero correlation entries; it generates its own correlations')
    count=int(numbers(lines[6])[0])
    if count<11:
        raise ValueError('template_twiss requires the 11-parameter MIRFQ layout: voltage, phase, then alpha/beta/emit for x,y,z')
    parameters=[]
    for line in lines[7:7+count]:
        fields=numbers(line)
        if len(fields)<4 or int(fields[3])!=0:
            raise ValueError('template_twiss supports fixed parameters only; disable fit flags')
        parameters.append(float(fields[0]))
    width=config['beam'].get('longitudinal_full_width_deg')
    report={'initialization':mode,'native_emittance_convention':'geometric in TRANSOPTR envelope convention, not normalized RMS',
            'longitudinal_full_width_deg':width,'width_policy':'preserve alpha_z and emittance_z; derive beta_z when width is set',
            'planes':{}}
    dimensions=[];correlations=[]
    for i,axis in enumerate('xyz'):
        alpha,beta_user,emit_user=parameters[2+3*i:5+3*i]
        beta=beta_user/units[7]
        emit=emit_user/(units[2*i]*units[2*i+1])
        if not all(math.isfinite(v) for v in (alpha,beta,emit)) or beta<=0 or emit<=0:
            raise ValueError(f'{axis}: alpha must be finite, beta and emittance positive')
        declared_beta=beta
        if axis=='z' and width is not None:
            half=half_length_cm(float(width),float(config['beam']['energy_initial_mev']),
                                float(config['beam']['mass_mev']),float(config['rf']['frequency_hz']))
            beta=half**2/emit
        size=math.sqrt(beta*emit)
        spread=math.sqrt((1+alpha*alpha)*emit/beta)
        rho=-alpha/math.sqrt(1+alpha*alpha)
        dimensions.extend([size*units[2*i],spread*units[2*i+1]])
        correlations.append(f'{2*i+1} {2*i+2} {rho:.12g}')
        report['planes'][axis]={'alpha':alpha,'declared_beta_cm':declared_beta,'effective_beta_cm':beta,
            'emittance_cm_rad':emit,'envelope_cm':size,'momentum_coordinate_spread':spread,'correlation':rho}
    result=list(lines[:3])+[' '.join(f'{v:.12g}' for v in dimensions)+' ! effective beam from template Twiss',
             lines[4],'3 ! generated x-xprime, y-yprime, z-delta correlations']+correlations+list(lines[6:])
    return result,report
