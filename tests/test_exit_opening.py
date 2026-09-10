import unittest
from types import SimpleNamespace
import numpy as np

from rfq_transop_design.design_rfq import default_config
from rfq_transop_design.sections import exit_count, exit_aperture, exit_parameters, focusing_scale, smoothstep
from rfq_transop_design.book_pages_23_32 import kt_coefficients_from_a_m_k
from rfq_transop_design.phase_control import field_nodes


class ExitOpeningTests(unittest.TestCase):
    def config(self):
        c = default_config()
        c['sections'] = {'rm_entrance_aperture_cm':1.2}
        c['exit'] = {'transition_cells':1, 'opening_cells':3,
                     'final_aperture_cm':'rm_entrance', 'fringe_length_cm':1.}
        return c

    def test_three_equal_segments_and_final_fringe_aperture(self):
        c = self.config()
        last = SimpleNamespace(B=10., m=1.8, Phi_deg=-30.)
        kin = SimpleNamespace(cell_length_cm=4.)
        self.assertEqual(exit_count(c),5)
        self.assertEqual(exit_aperture(c),1.2)
        for ordinal in (2,3,4):
            a,m,B,phi,L,label = exit_parameters(c,ordinal,last,1.,kin,kt_coefficients_from_a_m_k,12.)
            self.assertEqual((m,L,label),(1.,4.,'OM'))
            self.assertAlmostEqual(B,focusing_scale(c)/a**2)
        self.assertAlmostEqual(a,1.2)
        end = exit_parameters(c,5,last,1.,kin,kt_coefficients_from_a_m_k,12.)
        self.assertEqual((end[0],end[4],end[5]),(1.2,1.,'FF'))
        c['sections']['rm_entrance_aperture_cm']=1.4
        self.assertEqual(exit_aperture(c),1.4)

    def test_dense_field_uses_one_global_smooth_ramp(self):
        cells=[]
        ai,af,total=.45,1.2,12.
        for i in range(3):
            cells.append(SimpleNamespace(L_cm=4.,a_cm=ai+(af-ai)*smoothstep((i+1)/3),
                m=1., section='OM', entry_a_cm=ai,
                opening_start_cm=0.,opening_length_cm=total,
                opening_a_initial_cm=ai,opening_a_final_cm=af))
        nodes=field_nodes(cells,16)
        u=nodes[:,0]/total
        expected=ai+(af-ai)*(10*u**3-15*u**4+6*u**5)
        np.testing.assert_allclose(1/np.sqrt(nodes[:,1]),expected,atol=1e-12)
        np.testing.assert_array_equal(nodes[:,2],0.)
        self.assertTrue(np.all(np.diff(expected)>=0))

    def test_legacy_exit_remains_available(self):
        c=self.config();c['exit']['opening_cells']=0
        self.assertEqual(exit_count(c),3)
