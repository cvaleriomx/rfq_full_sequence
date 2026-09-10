"""Section and ideal two-term vane-tip profiles, with explicit end-model limits."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
try:
    from .phase_control import field_nodes
    from .sections import focusing_scale
except ImportError:
    from phase_control import field_nodes
    from sections import focusing_scale


def phase_samples(cells, points):
    lengths = np.array([c.L_cm for c in cells])
    secants = np.pi/lengths
    slopes = np.r_[secants[0],2/(1/secants[:-1]+1/secants[1:]),secants[-1]]
    chunks = []
    for i,L in enumerate(lengths):
        t = np.linspace(0,1,points,endpoint=False)
        chunks.extend(i*np.pi+(-2*t**3+3*t*t)*np.pi
                      +(t**3-2*t*t+t)*L*slopes[i]+(t**3-t*t)*L*slopes[i+1])
    return np.r_[chunks, len(cells)*np.pi]


def tip_root(a01, a10, k, cosine):
    # First equipotential root of A01*r^2 + A10*I0(k*r)*cos(psi)=1.
    # Scan to avoid incorrectly choosing an outer, nonphysical root.
    scale = 1/np.sqrt(a01)
    grid = np.linspace(0,2.5,161)[:,None]*scale[None,:]
    values = a01*grid**2 + a10*np.i0(k*grid)*cosine-1
    found = values >= 0
    valid = found.any(axis=0)
    first = found.argmax(axis=0)
    if not valid.all() or np.any(first == 0):
        raise ValueError('Two-term equipotential has no usable vane-tip root')
    columns = np.arange(len(a01))
    lo,hi = grid[first-1,columns],grid[first,columns]
    for _ in range(38):
        mid=(lo+hi)/2
        value=a01*mid**2+a10*np.i0(k*mid)*cosine-1
        hi=np.where(value>=0,mid,hi);lo=np.where(value<0,mid,lo)
    return (lo+hi)/2


def write_section_diagnostics(cells, config, output_dir):
    points = int(config['phase_control'].get('field_points_per_cell',16))
    nodes = field_nodes(cells,points)
    phase = phase_samples(cells,points)
    edges = np.r_[0.,[c.Z_cm for c in cells]]
    idx = np.minimum(np.searchsorted(edges[1:],nodes[:,0],side='left'),len(cells)-1)
    labels = np.array([c.section for c in cells])[idx]
    # FF is vacuum after the metal tips, represented only by a reduced field envelope.
    metal = labels != 'FF'
    metal &= nodes[:,1] > 0
    x,y = np.full(len(nodes),np.nan),np.full(len(nodes),np.nan)
    x[metal] = tip_root(nodes[metal,1],nodes[metal,2],nodes[metal,3],np.cos(phase[metal]))
    y[metal] = tip_root(nodes[metal,1],nodes[metal,2],nodes[metal,3],-np.cos(phase[metal]))
    pd.DataFrame({'z_cm':nodes[:,0],'section':labels,'metal_profile':metal,
        'x_plus_cm':x,'x_minus_cm':-x,'y_plus_cm':y,'y_minus_cm':-y,
        'A01_cm_inv2':nodes[:,1],'A10':nodes[:,2],'k_cm_inv':nodes[:,3],
        'spatial_phase_rad':phase,
        'B_two_term':nodes[:,1]*focusing_scale(config),
    }).to_csv(output_dir/'vane_geometry.csv', index=False)
    sections=[]
    for label in dict.fromkeys(c.section for c in cells):
        indices=[i for i,c in enumerate(cells) if c.section==label]
        sections.append({'section':label,'first_cell':indices[0]+1,'last_cell':indices[-1]+1,
            'count':len(indices),'z_start_cm':edges[indices[0]],'z_end_cm':edges[indices[-1]+1]})
    (output_dir/'sections.json').write_text(json.dumps({'scheme':config['sections']['scheme'],
        'section_schedule_only':True,'beam_matching_optimized':False,'sections':sections},indent=2)+'\n')
    profile=pd.read_csv(output_dir/'phase_validation.csv')
    colors=['#e8c44b','#a2d3ba','#a8c9e5','#c1b1dc','#e2a88e','#b9cdb9','#d2d2d2']
    fig,axs=plt.subplots(4,1,figsize=(12,10),sharex=True)
    controlled = profile.phase_controlled.astype(bool)
    axs[0].plot(profile.z_cm[controlled]/100,profile.phase_target_deg[controlled],label='Objetivo')
    axs[0].plot(profile.z_cm[controlled]/100,profile.phase_actual_deg[controlled],'--',label='TRANSOPTR')
    axs[0].set_ylabel('Fase [grados]');axs[0].legend(loc='lower right')
    axs[1].plot(np.array([c.Z_cm for c in cells])/100,[c.m for c in cells]);axs[1].set_ylabel('m')
    axs[2].plot(nodes[:,0]/100,nodes[:,1]*focusing_scale(config));axs[2].set_ylabel('B (dos términos)')
    axs[3].plot(profile.z_cm/100,profile.energy_final_run_mev);axs[3].set_ylabel('W [MeV]')
    axs[3].axhline(config['beam']['energy_target_mev'],color='gray',linestyle=':')
    for j,section in enumerate(sections):
        left,right=section['z_start_cm']/100,section['z_end_cm']/100
        for ax in axs:
            ax.axvspan(left,right,color=colors[j%len(colors)],alpha=.22)
        # Narrow end regions are shown clearly in the separate exit detail.
        if section['section'] not in {'TC','UM','OM','FF'}:
            axs[0].text((left+right)/2,1.03,section['section'],transform=axs[0].get_xaxis_transform(),ha='center')
    axs[-1].set_xlabel('z [m]')
    fig.suptitle('Secciones configuradas; no constituye optimización del haz',y=.995)
    fig.tight_layout();fig.savefig(output_dir/'rfq_sections.png');plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(12,7))
    for ax in axs:
        ax.plot(nodes[:,0],x,label='+x');ax.plot(nodes[:,0],y,label='+y')
        ax.plot(nodes[:,0],-x,color='C0');ax.plot(nodes[:,0],-y,color='C1')
        ax.set_ylabel('Punta de vane [cm]');ax.legend()
        for section in sections[-3:]:
            if section['section'] in {'TC','UM','OM','FF'}:
                l,r=section['z_start_cm'],section['z_end_cm']
                ax.axvspan(l,r,color=colors[{'TC':0,'UM':1,'OM':1,'FF':2}[section['section']]],alpha=.2)
    axs[1].set_xlim(max(0,sections[-3]['z_start_cm']-8),edges[-1])
    for section in sections[-3:]:
        axs[1].text((section['z_start_cm']+section['z_end_cm'])/2,1.02,section['section'],transform=axs[1].get_xaxis_transform(),ha='center')
    axs[1].set_xlabel('z [cm] — FF es campo de borde aproximado, sin perfil metálico')
    fig.tight_layout();fig.savefig(output_dir/'vane_geometry.png');plt.close(fig)
