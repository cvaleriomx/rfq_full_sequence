import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from rfq_transop_design.beam_diagnostics import read_twiss, write_beam_diagnostics


class BeamDiagnosticsTests(unittest.TestCase):
    def test_native_turns_and_duplicate_positions(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'fort.83'
            path.write_text('1 label 0 0 10 0.001 0 0 0\n2 label 1 0 10 0.001 0 0 0.25\n3 exit 1 0 10 0.001 0 0 0.3\n')
            table = read_twiss(path)
            self.assertEqual(len(table),2)
            self.assertEqual(table.mu_turns.iloc[-1]*360,90)

    def test_outputs_mask_degenerate_optics_and_vacuum(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            env=pd.DataFrame({'s':[0.,1.,2.,3.], 'E':[.03]*4,'beta':[.008]*4,
                'z-envelope':[.001]*4,"z'-envelope":[.001]*4,'r56':[0.]*4})
            for axis,corr,emit in [('x','r12',0.),('y','r34',.002)]:
                env[f'{axis}-envelope']=.1
                env[f"{axis}'-envelope"]=.02
                env[f'{axis}-emittance']=emit
                env[corr]=0.
                (root/f'fort.{81 if axis=="x" else 83}').write_text(
                    '\n'.join(f'{i+1} label {i} 0 5 .002 0 0 {i*.1}' for i in range(4)))
            pd.DataFrame({'z_cm':[0,1,2,3],'metal_profile':[True,True,True,False],
                'x_plus_cm':[.5,.5,.5,np.nan],'y_plus_cm':[.5,.5,.5,np.nan]}).to_csv(root/'vane_geometry.csv',index=False)
            (root/'sections.json').write_text(json.dumps({'sections':[
                {'section':'MA','z_start_cm':0,'z_end_cm':2},
                {'section':'FF','z_start_cm':2,'z_end_cm':3}]}))
            pd.DataFrame({'cell':[1,2,3],'Z_cm':[1,2,3],'section':['MA','MA','FF']}).to_csv(root/'history.csv',index=False)
            cfg={'beam':{'mass_mev':938.},'rf':{'frequency_hz':162e6},
                 'beam_diagnostics':{'envelope_multiplier':3,'max_pair_advance_deg':60}}
            reference=root/'matched.csv'
            pd.DataFrame({'z_cm':[0.,2.],'x_alpha':[0.,0.],'y_alpha':[0.,0.],
                          'x_beta_cm':[5.,5.],'y_beta_cm':[5.,5.]}).to_csv(reference,index=False)
            cfg['beam_diagnostics']['matched_reference_csv']=str(reference)
            with patch('rfq_transop_design.design_rfq.read_fort_envelope',return_value=env):
                report=write_beam_diagnostics(cfg,root,root)
            out=pd.read_csv(root/'beam_optics.csv')
            periods=pd.read_csv(root/'beam_phase_advance.csv')
            self.assertTrue(out.x_mu_deg.isna().all())
            self.assertTrue(out.x_alpha_projected.isna().all())
            self.assertAlmostEqual(out.y_occupancy.iloc[0],.6)
            self.assertTrue(np.isnan(out.y_occupancy.iloc[-1]))
            np.testing.assert_allclose(out.y_mismatch_bmag.iloc[:3],1.)
            self.assertTrue(np.isnan(out.y_mismatch_bmag.iloc[-1]))
            self.assertAlmostEqual(periods.y_advance_deg.iloc[0],72)
            self.assertFalse(report['x']['phase_usable'])
            self.assertEqual(report['y']['configured_limit_exceedances'],['max_pair_advance_deg'])
            self.assertFalse(report['beam_acceptance_validated'])
            self.assertTrue((root/'beam_optics.png').exists())


if __name__ == '__main__':
    unittest.main()
