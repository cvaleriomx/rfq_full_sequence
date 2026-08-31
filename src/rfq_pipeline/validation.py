from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np


FIELD_NAMES = ("ex_vpm_per_v", "ey_vpm_per_v", "ez_vpm_per_v")


def validate_field_map(path: str | Path, symmetry_tolerance: float = 5.0e-3) -> dict:
    path = Path(path)
    with h5py.File(path, "r") as h5:
        x = h5["axes/x_m"][:]
        y = h5["axes/y_m"][:]
        z = h5["axes/z_m"][:]
        expected_shape = (len(x), len(y), len(z))
        for name in FIELD_NAMES:
            dataset = h5[f"fields/{name}"]
            if dataset.shape != expected_shape:
                raise ValueError(f"Forma incorrecta para {name}: {dataset.shape}")

        indices = np.unique(np.linspace(0, len(z) - 1, min(9, len(z)), dtype=int))
        ex = h5["fields/ex_vpm_per_v"][:, :, indices]
        ey = h5["fields/ey_vpm_per_v"][:, :, indices]
        ez = h5["fields/ez_vpm_per_v"][:, :, indices]

    if not all(np.isfinite(field).all() for field in (ex, ey, ez)):
        raise ValueError("El mapa contiene NaN o infinito.")

    def relative_rms(residual: np.ndarray, reference: np.ndarray) -> float:
        scale = np.sqrt(np.mean(reference.astype(float) ** 2))
        return float(np.sqrt(np.mean(residual.astype(float) ** 2)) / max(scale, 1.0e-30))

    metrics = {
        "ex_odd_x": relative_rms(ex + ex[::-1, :, :], ex),
        "ey_odd_y": relative_rms(ey + ey[:, ::-1, :], ey),
        "ez_even_x": relative_rms(ez - ez[::-1, :, :], ez),
        "ez_even_y": relative_rms(ez - ez[:, ::-1, :], ez),
    }
    worst = max(metrics.values())
    if worst > symmetry_tolerance:
        raise ValueError(
            f"Falló la simetría del campo: error relativo máximo {worst:.3e} "
            f"> {symmetry_tolerance:.3e}."
        )
    return {
        "shape": list(expected_shape),
        "sampled_z_planes": int(len(indices)),
        "symmetry_relative_rms": metrics,
        "passed": True,
    }

