from pathlib import Path
import math
import unittest
from rfq_transop_design.compare_current import generate_input
from rfq_transop_design.design_rfq import default_config


class CurrentInputTests(unittest.TestCase):
    def test_current_and_rms_conversion(self):
        lines=Path('configs/mirfq1_transoptr.dat').read_text().splitlines()
        cfg=default_config()
        zero,ezero,rz=generate_input(lines,cfg,0.,.25)
        one,eone,r=generate_input(lines,cfg,1.,.25)
        self.assertAlmostEqual(float(one[0].split()[5]),.001/162e6,places=23)
        self.assertEqual(eone[1:],ezero[1:])
        dims=list(map(float,eone[3].split('!')[0].split()))
        bg=math.sqrt((1+.03/cfg['beam']['mass_mev'])**2-1)
        for i in (0,1):
            rho=float(eone[6+i].split()[2])
            eps=dims[2*i]*dims[2*i+1]*math.sqrt(1-rho*rho)
            self.assertAlmostEqual(eps/5*bg*1e4,.25,places=9)
        self.assertEqual(r['planes']['z']['emittance_cm_rad'],.00049)
        self.assertEqual(rz['bunch_charge_coulomb'],0.)

    def test_local_source_uses_supplied_charge(self):
        text=Path('rfq_transop_design/transoptr/RFQSC.f').read_text()
        self.assertIn('QSC=abs(CURRENT*CHARGE',text)
        self.assertNotIn('BNCHARGE*CHARGE',text)
        self.assertIn('SUBROUTINE RFQSC(',text)
