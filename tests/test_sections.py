import copy
import tempfile
import unittest
from pathlib import Path
import numpy as np
from rfq_transop_design.design_rfq import default_config, propose_cell
from rfq_transop_design.book_pages_23_32 import kt_coefficients_from_a_m_k
from rfq_transop_design.phase_control import field_nodes
from rfq_transop_design.sections import aperture_for_B, focusing_scale, stage_schedule, validate_sections
from rfq_transop_design.section_diagnostics import phase_samples, tip_root


class SectionTests(unittest.TestCase):
    def config(self):
        c=default_config()
        c['sections']={'enabled':True,'scheme':'nfsp','rm_cells':4}
        c['phase_control']={'enabled':True}
        c['iteration']['max_cells']=220
        return c

    def test_four_unmodulated_rm_cells_with_opening_aperture(self):
        c=self.config(); z=0; cells=[]
        for i in range(1,5):
            cell=propose_cell(c,i,.03,z); cells.append(cell); z=cell.Z_cm
        self.assertEqual([x.section for x in cells],['RM']*4)
        self.assertTrue(all(x.m==1 and x.A10==0 and x.Phi_deg==-90 for x in cells))
        self.assertTrue(all(x.dW_pred_MeV==0 for x in cells))
        self.assertTrue(np.all(np.diff([cells[0].entry_a_cm]+[x.a_cm for x in cells])<0))
        nodes=field_nodes(cells)
        self.assertTrue(np.all(nodes[:,2]==0))
        self.assertAlmostEqual(nodes[0,1],1/cells[0].entry_a_cm**2)
        self.assertEqual(propose_cell(c,5,.03,z).section,'MS')

    def test_B_changes_real_aperture_and_coefficients(self):
        c=self.config()
        for m in (1.,1.4,1.8):
            for k in (.7,2.,4.):
                a=aperture_for_B(10.,m,k,c,kt_coefficients_from_a_m_k)
                a01,_=kt_coefficients_from_a_m_k(a,m,k)
                self.assertAlmostEqual(focusing_scale(c)*a01,10.,places=8)
        with self.assertRaises(ValueError):
            aperture_for_B(0,1,2,c,kt_coefficients_from_a_m_k)

    def test_fsp_and_nfsp_section_order(self):
        c=self.config()
        self.assertEqual([s[0] for s in stage_schedule(c)],['MS','MB','MBA','MA'])
        c['sections']['scheme']='fsp'
        self.assertEqual([s[0] for s in stage_schedule(c)],['SH','GB','ACC'])
        c['sections']['rm_cells']=7
        validate_sections(c)  # Literature has longer matchers; no artificial five-cell cap.
        c['sections']['rm_cells']=0
        with self.assertRaises(ValueError): validate_sections(c)

    def test_fringe_ends_at_zero_without_modulation(self):
        c=self.config(); cells=[propose_cell(c,i,.03,(i-1)*.74) for i in range(1,5)]
        cells[-1].section='FF'; cells[-1].phase_match=False
        nodes=field_nodes(cells)
        self.assertEqual(nodes[-1,1],0.)
        self.assertTrue(np.all(nodes[:,2]==0))
        self.assertGreater(nodes[-17,1],0.)

    def test_tip_profile_satisfies_equipotential_and_phase_boundaries(self):
        c=self.config(); cells=[propose_cell(c,30,.05,0),propose_cell(c,31,.051,1)]
        nodes=field_nodes(cells); phase=phase_samples(cells,16)
        np.testing.assert_allclose(phase[::16],[0,np.pi,2*np.pi],atol=1e-12)
        for sign in (-1,1):
            r=tip_root(nodes[:,1],nodes[:,2],nodes[:,3],sign*np.cos(phase))
            v=nodes[:,1]*r*r+nodes[:,2]*np.i0(nodes[:,3]*r)*sign*np.cos(phase)
            np.testing.assert_allclose(v,1,atol=1e-9)
