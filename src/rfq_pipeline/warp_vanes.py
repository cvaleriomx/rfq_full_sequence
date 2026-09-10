"""Closed vane solids and DC field preparation using accelerator Warp."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd


def load_warp(serial=True):
    # Warp exposes this supported option before importing its MPI layer.
    import warpoptions
    warpoptions.ignoreUnknownArgs = True
    warpoptions.quietImport = True
    warpoptions.parser.set_defaults(serial=serial)
    import warp
    if not hasattr(warp, 'Triangles'):
        raise RuntimeError('Accelerator Warp with Triangles is required')
    return warp


def read_profiles(path):
    frame = pd.read_csv(path)
    required = ['z_cm','x_plus_cm','x_minus_cm','y_plus_cm','y_minus_cm','metal_profile']
    if not set(required).issubset(frame):
        raise ValueError('Use vane_geometry.csv exported by the sectioned designer')
    metal = frame.metal_profile.astype(str).str.lower().map({'true':True,'false':False})
    if metal.isna().any():
        raise ValueError('Invalid metal_profile flags')
    data = frame.loc[metal, required[:-1]].to_numpy(float)*.01
    if len(data)<2 or not np.isfinite(data).all() or np.any(np.diff(data[:,0])<=0):
        raise ValueError('Invalid metal profile coordinates')
    if np.any(data[:,[1,3]]<=0) or np.any(data[:,[2,4]]>=0):
        raise ValueError('Vanes must lie on their specified side of the axis')
    return data


def sweep_mesh(z, tips, radius, body_radius, arc_points=12, axis='x', sign=1):
    """Semicircular pole face with a rectangular back, swept along z.

    End caps close the metal at RM entrance and UM exit; no FF metal is added.
    Returned triangles have Warp's documented clockwise exterior winding.
    """
    if radius<=0 or body_radius<=max(tips)+radius or arc_points<4:
        raise ValueError('Body must extend beyond tip+radius, with >=4 arc samples')
    theta=np.linspace(-np.pi/2,np.pi/2,arc_points+1)
    cross=np.stack([tips[:,None]+radius*(1-np.cos(theta))[None,:],
                    np.broadcast_to(radius*np.sin(theta),(len(z),len(theta)))],axis=-1)
    backs=np.broadcast_to([[body_radius,radius],[body_radius,-radius]],(len(z),2,2))
    xy=np.concatenate([cross,backs],axis=1)
    xy[:,:,0]*=sign
    if axis=='y': xy=xy[:,:,::-1]
    rings=np.concatenate([xy,np.broadcast_to(z[:,None,None],(*xy.shape[:2],1))],axis=-1)
    nr,nv=rings.shape[:2]; vertices=rings.reshape(-1,3)
    faces=[]
    for i in range(nr-1):
        for j in range(nv):
            a=i*nv+j;b=i*nv+(j+1)%nv;c=(i+1)*nv+j;d=(i+1)*nv+(j+1)%nv
            faces.extend([[a,b,d],[a,d,c]])
    # Convex cross sections allow fan caps.
    for j in range(1,nv-1):
        faces.extend([[0,j+1,j],[(nr-1)*nv,(nr-1)*nv+j,(nr-1)*nv+j+1]])
    faces=np.array(faces,dtype=np.int64)
    tri=vertices[faces]
    volume=np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2])).sum()/6
    if volume>0: faces=faces[:,[0,2,1]]
    validate_mesh(vertices,faces)
    return vertices,faces


def validate_mesh(vertices,faces):
    tri=vertices[faces]
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    if not np.isfinite(tri).all() or np.any(area<=1e-18):
        raise ValueError('Degenerate/non-finite vane mesh')
    edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1)
    _,counts=np.unique(edges,axis=0,return_counts=True)
    if np.any(counts!=2): raise ValueError('Vane is not a closed two-manifold')


def build_geometry(config):
    s=config.section('warp_vanes'); source=config.resolve(s['profiles'])
    data=read_profiles(source)
    radius=float(s['tip_radius_m']); body=float(s['body_radius_m'])
    # Avoid intersection of adjacent pole bodies.
    if radius>=np.min(data[:,[1,3]]):
        raise ValueError('tip_radius_m must be smaller than the minimum aperture to avoid intersecting vanes')
    output=config.output_dir;output.mkdir(parents=True,exist_ok=True)
    items=[]
    for name,column,axis,sign in [('x_plus',1,'x',1),('x_minus',2,'x',-1),('y_plus',3,'y',1),('y_minus',4,'y',-1)]:
        v,f=sweep_mesh(data[:,0],abs(data[:,column]),radius,body,int(s.get('arc_points',12)),axis,sign)
        path=output/(name+'.npz');np.savez_compressed(path,vertices_m=v,faces=f)
        items.append({'name':name,'mesh':path.name,'condid':len(items)+1,
                      'voltage_v':float(s['intervane_voltage_v'])/2*(1 if axis=='x' else -1),
                      'vertices':len(v),'triangles':len(f)})
    manifest={'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
              'units':'m','winding':'clockwise viewed from exterior (Warp Triangles)',
              'tip_radius_m':radius,'body_radius_m':body,'metal_z_m':[float(data[0,0]),float(data[-1,0])],
              'geometry_model':'Swept semicircular tip and rectangular body; flat end caps, no FF metal',
              'conductors':items,'physical_validation':'not performed'}
    (output/'vanes_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    export_geometry_views(output)
    return manifest


def load_conductors(wp,manifest_path):
    path=Path(manifest_path);manifest=json.loads(path.read_text()); objects=[]
    for item in manifest['conductors']:
        with np.load(path.parent/item['mesh']) as d:
            triangles=np.asfortranarray(d['vertices_m'][d['faces']].transpose(2,1,0))
        objects.append(wp.Triangles(triangles,voltage=item['voltage_v'],condid=item['condid']))
    return objects


def prepare_warp(config,dry_run=False):
    s=config.section('warp_vanes'); data=read_profiles(config.resolve(s['profiles']))
    nx,ny,nz=[int(s[k]) for k in ('nx','ny','nz')]
    radius=float(s['domain_radius_m']); zmin=float(s.get('z_min_m',data[0,0]-.01));zmax=float(s.get('z_max_m',data[-1,0]+.02))
    if min(nx,ny,nz)<4 or nx%2 or ny%2 or zmin>=data[0,0] or zmax<=data[-1,0] or radius<=float(s['body_radius_m']):
        raise ValueError('Use even nx/ny, >=4 intervals, and a domain enclosing the complete metal')
    plan={'shape':[nx+1,ny+1,nz+1],'bounds_m':[[-radius,radius],[-radius,radius],[zmin,zmax]],
          'spacing_m':[2*radius/nx,2*radius/ny,(zmax-zmin)/nz],
          'metal_z_m':[float(data[0,0]),float(data[-1,0])],
          'field_storage_mib':4*8*(nx+1)*(ny+1)*(nz+1)/1024**2,
          'purpose':'DC geometry/field preparation; mesh convergence and beam acceptance still required'}
    if dry_run:return plan
    config.output_dir.mkdir(parents=True,exist_ok=True)
    # Never leave an old field paired with newly generated conductors after a failure.
    (config.output_dir/'fieldmap.h5').unlink(missing_ok=True)
    manifest=build_geometry(config)
    wp=load_warp(bool(s.get('serial',True)))
    wp.top.lprntpara=False;wp.top.lpsplots=False
    wp.w3d.nx,wp.w3d.ny,wp.w3d.nz=nx,ny,nz
    wp.w3d.xmmin=wp.w3d.ymmin=-radius;wp.w3d.xmmax=wp.w3d.ymmax=radius
    wp.w3d.zmmin=zmin;wp.w3d.zmmax=zmax
    wp.w3d.bound0=wp.w3d.boundnz=wp.w3d.boundxy=wp.dirichlet
    solver=wp.MultiGrid3D(mgtol=float(s.get('solver_tolerance_v',.01)),mgmaxiters=int(s.get('solver_max_iters',500)))
    wp.registersolver(solver)
    conductors=load_conductors(wp,config.output_dir/'vanes_manifest.json')
    for conductor in conductors: solver.installconductor(conductor)
    wp.package('w3d');wp.generate();wp.fieldsolve()
    residual=float(solver.mgerror)
    report={**plan,'solver_residual_v':residual,'solver_converged':bool(np.isfinite(residual) and residual<=float(s.get('solver_tolerance_v',.01))),
            'manifest':str(config.output_dir/'vanes_manifest.json'),'boundary_condition':'grounded rectangular outer boundary',
            'beam_dynamics_validated':False}
    (config.output_dir/'warp_preparation.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['solver_converged']:raise RuntimeError('Warp DC solver did not converge; inspect warp_preparation.json')
    phi=np.asarray(wp.getphi()).copy()
    if phi.shape!=(nx+1,ny+1,nz+1):raise ValueError(f'Unexpected Warp potential shape {phi.shape}')
    voltage=float(s['intervane_voltage_v']);axes=[np.linspace(-radius,radius,nx+1),np.linspace(-radius,radius,ny+1),np.linspace(zmin,zmax,nz+1)]
    fields=[-d for d in np.gradient(phi/voltage,*axes,edge_order=2)]
    if not np.isfinite(phi).all():raise ValueError('Non-finite Warp potential')
    with h5py.File(config.output_dir/'fieldmap.incomplete.h5','w') as h:
        h.attrs.update(source='Warp MultiGrid3D conductor DC solution',axis_order='xyz',normalization='per V_intervane',geometry_manifest='vanes_manifest.json')
        for name,axis in zip('xyz',axes):h.create_dataset('axes/'+name+'_m',data=axis)
        for name,field in zip(['phi_per_v','ex_vpm_per_v','ey_vpm_per_v','ez_vpm_per_v'],[phi/voltage,*fields]):
            h.create_dataset('fields/'+name,data=field,compression='gzip',shuffle=True)
    (config.output_dir/'fieldmap.incomplete.h5').replace(config.output_dir/'fieldmap.h5')
    pd.DataFrame({'z_m':axes[2],'ez_vpm_per_v':fields[2][nx//2,ny//2]}).to_csv(config.output_dir/'warp_dc_axis.csv',index=False)
    from .diagnostics import plot_field_diagnostics
    plot_field_diagnostics(config.output_dir/'fieldmap.h5',config.output_dir/'warp_dc_field.png')
    report['fieldmap']=str(config.output_dir/'fieldmap.h5')
    report['design_comparison']=compare_design_axis(config)
    (config.output_dir/'warp_preparation.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def load_scraper_conductors(wp,manifest_path):
    """Fast inside tests for the exact polygonal sweeps used in the DC mesh.

    Uses the same arc vertices and linear longitudinal interpolation. No
    nearest-triangle search is needed at every point of the scraper grid.
    """
    path=Path(manifest_path);manifest=json.loads(path.read_text());objects=[]
    radius=manifest['tip_radius_m'];body=manifest['body_radius_m']
    for item in manifest['conductors']:
        with np.load(path.parent/item['mesh']) as d: vertices=d['vertices_m']
        axis=0 if item['name'].startswith('x') else 1
        sign=1 if item['name'].endswith('plus') else -1
        z,inverse=np.unique(vertices[:,2],return_inverse=True)
        tips=np.full(len(z),np.inf);np.minimum.at(tips,inverse,vertices[:,axis]*sign)
        transverse=np.unique(vertices[:,1-axis])
        theta=np.arcsin(np.clip(transverse/radius,-1,1))
        front=radius*(1-np.cos(theta))
        class SweptScraper(wp.Assembly):
            def __init__(self,axis,sign,z,tips,transverse,front,item):
                self.axis,self.sign,self.znodes,self.tips=axis,sign,z,tips
                self.transverse,self.front=transverse,front
                wp.Assembly.__init__(self,v=item['voltage_v'],condid=item['condid'],generatord=self.signed_distance)
            def signed_distance(self,xcent,ycent,zcent,n,x,y,z,distance):
                along=(x if self.axis==0 else y)*self.sign
                across=y if self.axis==0 else x
                face=np.interp(z,self.znodes,self.tips)+np.interp(across,self.transverse,self.front)
                distance[:]=np.maximum.reduce([face-along,along-body,abs(across)-radius,self.znodes[0]-z,z-self.znodes[-1]])
            def getextent(self):
                return wp.ConductorExtent([-body,-body,float(self.znodes[0])],[body,body,float(self.znodes[-1])],[0.,0.,0.])
        objects.append(SweptScraper(axis,sign,z,tips,transverse,front,item))
    return objects


def export_geometry_views(output):
    """STL uses outward normals; NPZ retains the clockwise winding for Warp."""
    import struct
    import matplotlib.pyplot as plt
    output=Path(output);manifest=json.loads((output/'vanes_manifest.json').read_text())
    fig=plt.figure(figsize=(12,5));ax=fig.add_subplot(111,projection='3d')
    for i,item in enumerate(manifest['conductors']):
        with np.load(output/item['mesh']) as d:
            v,f=d['vertices_m'],d['faces']
        tri=v[f[:,[0,2,1]]]
        normals=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None]
        records=np.zeros(len(tri),dtype=[('normal','<f4',(3,)),('vertices','<f4',(3,3)),('attribute','<u2')])
        records['normal']=normals;records['vertices']=tri
        with (output/(item['name']+'.stl')).open('wb') as stream:
            stream.write(b'RFQ vane, units metres'.ljust(80,b' '));stream.write(struct.pack('<I',len(tri)));stream.write(records.tobytes())
        nz=len(np.unique(v[:,2]));rings=v.reshape(nz,-1,3)
        ax.plot_surface(rings[:,:,2],rings[:,:,0]*1000,rings[:,:,1]*1000,rstride=max(1,nz//120),cstride=1,color=['#dc7443','#c25230','#4285b4','#326790'][i],alpha=.85)
    ax.set_xlabel('z [m]');ax.set_ylabel('x [mm]');ax.set_zlabel('y [mm]');ax.set_box_aspect((4,1,1))
    ax.set_title('Cuatro vanes cerradas — extremos planos; escalas de ejes diferentes')
    fig.tight_layout();fig.savefig(output/'vanes_3d.png');plt.close(fig)


def compare_design_axis(config):
    """Compare the DC solution with the dense design potential on axis."""
    source=pd.read_csv(config.resolve(config.section('warp_vanes')['profiles']))
    z=source.z_cm.to_numpy()*.01
    reference_phi=.5*source.A10.to_numpy()*np.cos(source.spatial_phase_rad.to_numpy())
    with h5py.File(config.output_dir/'fieldmap.h5') as h:
        axis=h['axes/z_m'][:];shape=h['fields/ez_vpm_per_v'].shape
        ez=h['fields/ez_vpm_per_v'][shape[0]//2,shape[1]//2,:]
    reference=-np.gradient(np.interp(axis,z,reference_phi,left=0.,right=0.),axis,edge_order=2)
    selected=(axis>=z[0])&(axis<=z[-1]); residual=ez[selected]-reference[selected]
    scale=float(np.sqrt(np.mean(reference[selected]**2)))
    report={'reference':'Derivative of dense design on-axis two-term potential',
            'relative_rms_error':float(np.sqrt(np.mean(residual**2))/max(scale,1e-30)),
            'correlation':float(np.corrcoef(ez[selected],reference[selected])[0,1]),
            'physical_equivalence_validated':False,
            'note':'Finite circular pole faces and grounded boundary differ from ideal two-term electrodes; grid and boundary convergence still required.'}
    (config.output_dir/'warp_design_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(11,4));ax.plot(axis,reference,label='Diseño, potencial de dos términos');ax.plot(axis,ez,'--',label='Warp DC, conductores')
    ax.set_xlabel('z [m]');ax.set_ylabel('Ez / V_intervane [1/m]');ax.legend();fig.tight_layout();fig.savefig(config.output_dir/'warp_design_comparison.png');plt.close(fig)
    return report
