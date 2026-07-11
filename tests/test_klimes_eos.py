import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts import run_literature_eos_benchmarks as benchmark


class MurnaghanTests(unittest.TestCase):
    def test_recovers_synthetic_curve(self):
        volumes = np.linspace(40.0, 60.0, 11)
        e0 = -10.0
        v0 = 50.0
        b0 = 0.02
        b1 = 4.5
        energies = e0 + b0 * benchmark.murnaghan_feature(volumes, v0, b1)

        fit = benchmark.fit_murnaghan(volumes, energies)

        self.assertAlmostEqual(float(fit["E0_hartree"]), e0, places=8)
        self.assertAlmostEqual(float(fit["V0_A3"]), v0, places=5)
        self.assertAlmostEqual(float(fit["B0_hartree_A3"]), b0, places=6)
        self.assertAlmostEqual(float(fit["B1"]), b1, places=4)

    def test_local_branch_near_reference_is_selected(self):
        factors = np.linspace(0.8, 1.2, 11)
        volumes = 100.0 * factors
        energies = -5.0 + 0.01 * benchmark.murnaghan_feature(volumes, 100.0, 4.0)
        points = [
            {"volume_factor": float(factor), "volume_A3": float(volume), "energy_hartree": float(energy)}
            for factor, volume, energy in zip(factors, volumes, energies)
        ]
        for factor in np.linspace(1.5, 1.9, 5):
            points.append(
                {
                    "volume_factor": float(factor),
                    "volume_A3": float(100.0 * factor),
                    "energy_hartree": float(-6.0 + 0.2 * (factor - 1.7) ** 2),
                }
            )

        window, fit = benchmark.select_eos_fit(points)

        self.assertLessEqual(max(float(row["volume_factor"]) for row in window), 1.2)
        self.assertAlmostEqual(float(fit["V0_A3"]), 100.0, places=4)
        self.assertEqual(benchmark.eos_fit_reason(fit), "ok")


class InputAndParserTests(unittest.TestCase):
    def test_parser_prefers_zero_temperature_energy(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "cp2k.out"
            output.write_text(
                "Total energy (extrapolated to T->0): -1.234500000000\n"
                "ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -1.235000000000\n"
                "PROGRAM ENDED AT 2026-01-01\n"
            )

            ok, energy, reason = benchmark.parse_cp2k_output(output)

        self.assertTrue(ok)
        self.assertEqual(reason, "ok")
        self.assertAlmostEqual(float(energy), -1.2345, places=12)

    def test_metal_input_uses_native_full_bloch_grid(self):
        copper = next(item for item in benchmark.load_structures("klimes") if item.system == "Cu")
        job = benchmark.Job("klimes", copper, "GFN2", 1.0, Path("cp2k.ssmp"), 1, True)

        text = benchmark.cp2k_input(job, 0.4, 500, 250)

        self.assertIn("SCHEME MACDONALD 16 16 16", text)
        self.assertIn("FULL_GRID T", text)
        self.assertIn("GFN_TYPE TBLITE", text)
        self.assertNotIn("MULTIPLE_UNIT_CELL", text)

    def test_force_option_is_available_for_reproduction(self):
        args = benchmark.parser().parse_args(["run", "--force"])
        self.assertTrue(args.force)


if __name__ == "__main__":
    unittest.main()
