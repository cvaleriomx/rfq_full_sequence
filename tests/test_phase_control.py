import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from rfq_transop_design.phase_control import field_nodes, phase_residual, tune_cell


class PhaseControlTests(unittest.TestCase):
    def test_positive_k_and_pi_integral_per_cell(self):
        cells = [SimpleNamespace(L_cm=l, a_cm=.4, m=1.2) for l in [.7, .8, 1.1, 1.3]]
        points = 16
        nodes = field_nodes(cells, points)
        self.assertTrue(np.all(nodes[:, 3] > 0))
        for i, cell in enumerate(cells):
            k = nodes[i*points:(i+1)*points+1, 3]
            integral = cell.L_cm/points/3*(k[0]+k[-1]+4*k[1:-1:2].sum()+2*k[2:-1:2].sum())
            self.assertAlmostEqual(integral, np.pi, places=12)
        self.assertTrue(np.isfinite(nodes).all())

    def test_phase_file_requires_adding_global_phase_without_wrapping(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'fort.20'
            np.savetxt(p, [[1, np.deg2rad(420)]])
            error, actual = phase_residual(p, -90, -30)
            self.assertAlmostEqual(actual, 330)
            self.assertAlmostEqual(error, 360)

    def test_tuning_changes_geometry_and_keeps_rf_phase_fixed(self):
        with tempfile.TemporaryDirectory() as d:
            cell = SimpleNamespace(L_cm=1., Phi_deg=-30., cell=1)
            config = {'phase_control': {}, 'rf': {'tracking_phase_deg': -90.}}
            def feedback(cells, cfg, directory):
                self.assertEqual(cfg['rf']['tracking_phase_deg'], -90)
                # Root at 1.1 cm; fort.20 omits the input RF phase.
                actual = -30+180*(cells[-1].L_cm-1.1)
                np.savetxt(directory/'fort.20', [[1, np.deg2rad(actual+90)]])
                return .04, 'optr_ok'
            def update(c, length):
                c.L_cm = length
            tune_cell([cell], config, Path(d), feedback, update)
            self.assertAlmostEqual(cell.L_cm, 1.1, places=6)
            self.assertAlmostEqual(cell.phase_actual_deg, -30, places=6)

    def test_failed_tracker_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as d:
            cell = SimpleNamespace(L_cm=1., Phi_deg=-30., cell=1)
            cfg = {'phase_control': {}, 'rf': {'tracking_phase_deg': -90.}}
            with self.assertRaises(RuntimeError):
                tune_cell([cell], cfg, Path(d), lambda *args: (.03, 'optr_failed_1'), lambda *args: None)

    def test_average_gain_uses_transoptr_relative_phase(self):
        from rfq_transop_design.book_pages_23_32 import predicted_energy_gain_mev
        self.assertAlmostEqual(predicted_energy_gain_mev(50, 1, .2, -90), 0., places=12)
        self.assertAlmostEqual(predicted_energy_gain_mev(50, 1, .2, 0), np.pi/4*.05*.2)
        self.assertLess(predicted_energy_gain_mev(50, 1, .2, 120), 0.)
