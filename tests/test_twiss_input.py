import math
from pathlib import Path
import unittest
from rfq_transop_design.beam_input import prepare_beam_input, half_length_cm
from rfq_transop_design.design_rfq import default_config


class TwissInputTests(unittest.TestCase):
    def setUp(self):
        self.lines=Path('configs/mirfq1_transoptr.dat').read_text().splitlines()
        self.config=default_config()
        self.config['beam']['initialization']='template_twiss'

    def test_all_planes_have_requested_emittance_and_correlation(self):
        lines,report=prepare_beam_input(self.lines,self.config)
        dims=list(map(float,lines[3].split('!')[0].split()))
        self.assertEqual(lines[5].split()[0],'3')
        for i,(alpha,beta,emit) in enumerate([(2.39,31.5,.003),(2.70,33.6,.003),(2.,47.27,.00049)]):
            rho=float(lines[6+i].split()[2])
            reconstructed=dims[2*i]*dims[2*i+1]*math.sqrt(1-rho*rho)
            self.assertAlmostEqual(reconstructed,emit,places=11)
            self.assertAlmostEqual(-rho*dims[2*i]*dims[2*i+1]/emit,alpha,places=9)
            self.assertEqual(report['planes']['xyz'[i]]['declared_beta_cm'],beta)
        self.assertAlmostEqual(dims[4],half_length_cm(90.,.03,self.config['beam']['mass_mev'],162e6))
        self.assertNotEqual(report['planes']['z']['effective_beta_cm'],47.27)
        self.assertEqual(lines[9:],self.lines[6:])

    def test_null_width_uses_longitudinal_beta(self):
        self.config['beam']['longitudinal_full_width_deg']=None
        lines,report=prepare_beam_input(self.lines,self.config)
        self.assertEqual(report['planes']['z']['effective_beta_cm'],47.27)
        self.assertAlmostEqual(float(lines[3].split()[4]),math.sqrt(47.27*.00049))

    def test_nonunit_conversions_match_cic3(self):
        self.lines[4]='10 1000 10 1000 10 1000 1 10'
        _,report=prepare_beam_input(self.lines,self.config)
        self.assertAlmostEqual(report['planes']['x']['effective_beta_cm'],3.15)
        self.assertAlmostEqual(report['planes']['x']['emittance_cm_rad'],.003/10000)

    def test_reject_fitting_or_existing_coupling(self):
        self.lines[5]='1'
        with self.assertRaises(ValueError):prepare_beam_input(self.lines,self.config)
        self.lines[5]='0';self.lines[11]='.003 0 100 1'
        with self.assertRaises(ValueError):prepare_beam_input(self.lines,self.config)
