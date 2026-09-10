import unittest
from rfq_transop_design.compare_bunch_widths import half_length_cm
from rfq_transop_design.beam_input import apply_longitudinal_width
from rfq_transop_design.design_rfq import default_config


class BunchWidthTests(unittest.TestCase):
    def test_default_width_updates_only_longitudinal_size(self):
        config=default_config()
        lines=['beam','mode','field','0.1 0.02 0.1 0.02 28.0 0.01','1 1 1 1 10 1 1 1']
        result=apply_longitudinal_width(lines,config)
        expected=half_length_cm(90.,config['beam']['energy_initial_mev'],
                                config['beam']['mass_mev'],config['rf']['frequency_hz'])
        self.assertAlmostEqual(float(result[3].split()[4]),expected*10)
        self.assertEqual(result[3].split()[:4],lines[3].split()[:4])
        self.assertEqual(result[3].split()[5],'0.01')
        self.assertEqual(lines[3].split()[4],'28.0')
        config['beam']['longitudinal_full_width_deg']=None
        self.assertEqual(apply_longitudinal_width(lines,config),lines)

    def test_full_width_and_centimeters(self):
        half=half_length_cm(180.,.03,938.272089,162e6)
        self.assertAlmostEqual(half,.36995296,places=7)
        self.assertAlmostEqual(half_length_cm(90.,.03,938.272089,162e6),half/2)
        self.assertAlmostEqual(half_length_cm(180.,.03,938.272089,324e6),half/2)

    def test_reject_invalid_widths(self):
        for value in (0.,-90.,361.,float('nan')):
            with self.assertRaises(ValueError):
                half_length_cm(value,.03,938.272089,162e6)
