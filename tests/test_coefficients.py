from pathlib import Path
import unittest

import numpy as np

from rfq_pipeline.coefficients import derive_two_term_coefficients
from rfq_pipeline.design import read_design_table


ROOT = Path(__file__).resolve().parents[1]


class CoefficientTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        frame = read_design_table(ROOT / "data/designs/isac2/table1.txt")
        cls.coefficients = derive_two_term_coefficients(frame)

    def test_one_row_per_physical_cell(self):
        self.assertEqual(len(self.coefficients), 87)

    def test_each_cell_adds_pi_phase(self):
        phase = self.coefficients["phase_end_rad"].to_numpy()
        self.assertTrue(np.allclose(np.diff(np.r_[0.0, phase]), np.pi))

    def test_coefficients_are_finite(self):
        columns = ["a01_per_m2", "a10", "k_rad_per_m"]
        self.assertTrue(np.isfinite(self.coefficients[columns]).all().all())


if __name__ == "__main__":
    unittest.main()

