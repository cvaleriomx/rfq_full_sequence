from pathlib import Path
import unittest

from rfq_pipeline.design import read_design_table, validate_design


ROOT = Path(__file__).resolve().parents[1]


class DesignTests(unittest.TestCase):
    def test_canonical_isac2_table(self):
        frame = read_design_table(ROOT / "data/designs/isac2/table1.txt")
        report = validate_design(frame, expected_cells=87, expected_length_m=2.0582)
        self.assertEqual(report["last_cell"], "87T")
        self.assertAlmostEqual(report["voltage_kv_min"], 41.53)
        self.assertAlmostEqual(report["output_energy_mev"], 1.0167)


if __name__ == "__main__":
    unittest.main()

