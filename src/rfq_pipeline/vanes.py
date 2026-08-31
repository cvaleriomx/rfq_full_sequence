from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class VaneProfiles:
    z_m: np.ndarray
    x_plus_m: np.ndarray
    x_minus_m: np.ndarray
    y_plus_m: np.ndarray
    y_minus_m: np.ndarray


def build_vane_profiles(frame: pd.DataFrame, points_per_cell: int = 80) -> VaneProfiles:
    """Cosine profiles with smooth aperture matching at cell boundaries."""
    if points_per_cell < 4:
        raise ValueError("points_per_cell debe ser al menos 4.")
    cells = frame.iloc[1:].copy().reset_index(drop=True)
    numbers = cells["Cell"].astype(str).str.extract(r"(\d+)", expand=False).astype(int).to_numpy()
    a = cells["a"].to_numpy(float) * 1.0e-2
    m = cells["m"].to_numpy(float)
    z_end = cells["Z"].to_numpy(float) * 1.0e-2
    z_start = np.r_[0.0, z_end[:-1]]
    lengths = z_end - z_start
    if np.any(lengths <= 0.0):
        raise ValueError("Las fronteras longitudinales de las celdas no son válidas.")

    z_parts: list[np.ndarray] = []
    x_parts: list[np.ndarray] = []
    y_parts: list[np.ndarray] = []
    for index, (number, begin, end) in enumerate(zip(numbers, z_start, z_end)):
        include_end = index == len(cells) - 1
        z = np.linspace(begin, end, points_per_cell + int(include_end), endpoint=include_end)
        u = (z - begin) / (end - begin)
        smooth = 3.0 * u**2 - 2.0 * u**3

        previous_a = a[index - 1] if index > 0 else a[index]
        previous_ma = a[index - 1] * m[index - 1] if index > 0 else a[index] * m[index]
        next_a = a[index + 1] if index + 1 < len(a) else a[index]
        next_ma = a[index + 1] * m[index + 1] if index + 1 < len(a) else a[index] * m[index]
        aperture_factor = (1.0 - smooth) * 0.5 * (1.0 + previous_a / a[index]) + smooth * 0.5 * (
            1.0 + next_a / a[index]
        )
        ma = a[index] * m[index]
        modulation_factor = (1.0 - smooth) * 0.5 * (1.0 + previous_ma / ma) + smooth * 0.5 * (
            1.0 + next_ma / ma
        )
        aperture_profile = a[index] * aperture_factor
        ma_profile = ma * modulation_factor
        delta = ma_profile - aperture_profile
        cosine = np.cos(np.pi * u)
        if (-1.0) ** (number + 1) > 0.0:
            x = aperture_profile + 0.5 * delta * (1.0 - cosine)
            y = aperture_profile + 0.5 * delta * (1.0 + cosine)
        else:
            x = aperture_profile + 0.5 * delta * (1.0 + cosine)
            y = aperture_profile + 0.5 * delta * (1.0 - cosine)
        z_parts.append(z)
        x_parts.append(x)
        y_parts.append(y)

    z = np.concatenate(z_parts)
    x = np.concatenate(x_parts)
    y = np.concatenate(y_parts)
    return VaneProfiles(z, x, -x, y, -y)


def write_vane_profiles(profiles: VaneProfiles, csv_path: Path, plot_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "z_m": profiles.z_m,
            "x_plus_m": profiles.x_plus_m,
            "x_minus_m": profiles.x_minus_m,
            "y_plus_m": profiles.y_plus_m,
            "y_minus_m": profiles.y_minus_m,
        }
    ).to_csv(csv_path, index=False)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(12, 4))
    axis.plot(profiles.z_m, profiles.x_plus_m * 1.0e3, label="punta vane +x", lw=0.9)
    axis.plot(profiles.z_m, profiles.y_plus_m * 1.0e3, label="punta vane +y", lw=0.9)
    axis.set(xlabel="z [m]", ylabel="radio de punta [mm]", title="Perfiles de vanes derivados de table1.txt")
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)

