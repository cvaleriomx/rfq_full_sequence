from pathlib import Path
import tempfile
import unittest

from rfq_pipeline.coefficients import derive_two_term_coefficients
from rfq_pipeline.design import read_design_table
from rfq_pipeline.fieldmap import generate_field_map
from rfq_pipeline.validation import validate_field_map


ROOT = Path(__file__).resolve().parents[1]


class FieldMapTests(unittest.TestCase):
    def test_small_map_has_expected_symmetry(self):
        frame = read_design_table(ROOT / "data/designs/isac2/table1.txt")
        coefficients = derive_two_term_coefficients(frame)
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "small.h5"
            metadata = generate_field_map(
                coefficients,
                target,
                x_min_m=-0.002,
                x_max_m=0.002,
                dx_m=0.001,
                y_min_m=-0.002,
                y_max_m=0.002,
                dy_m=0.001,
                z_min_m=0.0,
                z_max_m=2.0582,
                dz_m=0.1,
                chunk_z=8,
            )
            self.assertEqual(metadata["shape"][:2], [5, 5])
            report = validate_field_map(target)
            self.assertTrue(report["passed"])


if __name__ == "__main__":
    unittest.main()

