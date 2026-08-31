from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def derive_two_term_coefficients(frame: pd.DataFrame) -> pd.DataFrame:
    """Derive the Kapchinsky–Teplyakov two-term coefficients in SI units."""
    cells = frame.iloc[1:].copy().reset_index(drop=True)
    if cells.empty:
        raise ValueError("No hay celdas físicas.")

    aperture = cells["a"].to_numpy(float) * 1.0e-2
    modulation = cells["m"].to_numpy(float)
    length = cells["L"].to_numpy(float) * 1.0e-2
    z_end = cells["Z"].to_numpy(float) * 1.0e-2
    k = np.pi / length
    i0_ka = np.i0(k * aperture)
    i0_kma = np.i0(k * modulation * aperture)
    a10 = (modulation**2 - 1.0) / (modulation**2 * i0_ka + i0_kma)
    a01 = (1.0 - a10 * i0_ka) / aperture**2

    result = pd.DataFrame(
        {
            "cell": cells["Cell"].astype(str),
            "z_end_m": z_end,
            "cell_length_m": length,
            "aperture_m": aperture,
            "modulation": modulation,
            "k_rad_per_m": k,
            "phase_end_rad": np.pi * np.arange(1, len(cells) + 1),
            "a01_per_m2": a01,
            "a10": a10,
            "a10_transoptr": cells["A10"].to_numpy(float),
        }
    )
    return result


def coefficient_nodes(coefficients: pd.DataFrame) -> dict[str, np.ndarray]:
    """Add the z=0 boundary used for interpolation on a regular field grid."""
    first = coefficients.iloc[0]
    return {
        "z_m": np.r_[0.0, coefficients["z_end_m"].to_numpy(float)],
        "phase_rad": np.r_[0.0, coefficients["phase_end_rad"].to_numpy(float)],
        "k_rad_per_m": np.r_[first["k_rad_per_m"], coefficients["k_rad_per_m"]],
        "a01_per_m2": np.r_[first["a01_per_m2"], coefficients["a01_per_m2"]],
        "a10": np.r_[first["a10"], coefficients["a10"]],
    }


def write_coefficients(coefficients: pd.DataFrame, csv_path: Path, fort_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(csv_path, index=False)

    # Compatibilidad con herramientas antiguas: z[cm], A01[cm^-2], A10, k[cm^-1].
    nodes = coefficient_nodes(coefficients)
    legacy = np.column_stack(
        [
            nodes["z_m"] * 100.0,
            nodes["a01_per_m2"] / 1.0e4,
            nodes["a10"],
            nodes["k_rad_per_m"] / 100.0,
        ]
    )
    np.savetxt(
        fort_path,
        legacy,
        fmt="%.10e",
        header="z_cm A01_cm^-2 A10 k_cm^-1 (DERIVADO; no es entrada canónica)",
    )

