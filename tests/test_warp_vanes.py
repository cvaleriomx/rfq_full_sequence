import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from rfq_pipeline.warp_vanes import read_profiles,sweep_mesh,validate_mesh


class WarpVaneTests(unittest.TestCase):
    def test_mesh_is_closed_and_has_warp_winding(self):
        for axis in ['x','y']:
            for sign in [-1,1]:
                v,f=sweep_mesh(np.array([0.,.01,.03]),np.array([.004,.005,.004]),.002,.018,8,axis,sign)
                validate_mesh(v,f)
                t=v[f]
                volume=np.einsum('ij,ij->i',t[:,0],np.cross(t[:,1],t[:,2])).sum()/6
                self.assertLess(volume,0)
                longitudinal=0 if axis=='x' else 1
                self.assertGreaterEqual((v[:,longitudinal]*sign).min(),.004-1e-12)
                self.assertEqual(v[:,2].max(),.03)
                with self.assertRaises(ValueError):validate_mesh(v,f[:-1])

    def test_profile_excludes_fringe_and_converts_cm_to_m(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'profiles.csv'
            pd.DataFrame({'z_cm':[0,100,101],'x_plus_cm':[.4,.4,np.nan],
                'x_minus_cm':[-.4,-.4,np.nan],'y_plus_cm':[.5,.5,np.nan],
                'y_minus_cm':[-.5,-.5,np.nan],'metal_profile':[True,True,False]}).to_csv(p,index=False)
            data=read_profiles(p)
            self.assertEqual(data.shape,(2,5)); self.assertEqual(data[-1,0],1.)
            self.assertEqual(data[0,1],.004)
