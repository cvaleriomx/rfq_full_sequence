from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from .coefficients import coefficient_nodes


def regular_axis(lower: float, upper: float, requested_step: float) -> np.ndarray:
    if upper <= lower or requested_step <= 0.0:
        raise ValueError("Los límites y el paso de la malla deben ser positivos y ordenados.")
    intervals = max(1, int(round((upper - lower) / requested_step)))
    return np.linspace(lower, upper, intervals + 1, dtype=np.float64)


def _interpolated_coefficients(z: np.ndarray, nodes: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    source_z = nodes["z_m"]
    return {
        name: np.interp(z, source_z, values)
        for name, values in nodes.items()
        if name != "z_m"
    }


def _potential_per_intervane_volt(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    nodes: dict[str, np.ndarray],
) -> np.ndarray:
    """Two-term potential, normalized to a 1 V difference between opposite vanes."""
    c = _interpolated_coefficients(z, nodes)
    xx = x[:, None, None]
    yy = y[None, :, None]
    radius = np.sqrt(xx**2 + yy**2)
    k = c["k_rad_per_m"][None, None, :]
    transverse = c["a01_per_m2"][None, None, :] * (xx**2 - yy**2)
    longitudinal = c["a10"][None, None, :] * np.i0(k * radius) * np.cos(
        c["phase_rad"][None, None, :]
    )
    return 0.5 * (transverse + longitudinal)


def generate_field_map(
    coefficients: pd.DataFrame,
    output_path: str | Path,
    *,
    x_min_m: float,
    x_max_m: float,
    dx_m: float,
    y_min_m: float,
    y_max_m: float,
    dy_m: float,
    z_min_m: float,
    z_max_m: float,
    dz_m: float,
    chunk_z: int = 64,
) -> dict:
    """Write a compact, chunked HDF5 field map in Warp-compatible axis order."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    x = regular_axis(x_min_m, x_max_m, dx_m)
    y = regular_axis(y_min_m, y_max_m, dy_m)
    z = regular_axis(z_min_m, z_max_m, dz_m)
    nodes = coefficient_nodes(coefficients)
    if z[0] < nodes["z_m"][0] or z[-1] > nodes["z_m"][-1] + 1.0e-12:
        raise ValueError("La malla z sale del intervalo cubierto por los coeficientes.")
    if min(len(x), len(y), len(z)) < 3:
        raise ValueError("Cada eje necesita al menos tres puntos para calcular gradientes.")

    shape = (len(x), len(y), len(z))
    chunks = (len(x), len(y), min(chunk_z, len(z)))
    with h5py.File(output_path, "w") as h5:
        h5.attrs["format"] = "rfq-fieldmap-v1"
        h5.attrs["axis_order"] = "x,y,z"
        h5.attrs["coordinate_unit"] = "m"
        h5.attrs["field_unit"] = "V/m per V_intervane"
        h5.attrs["potential_unit"] = "V per V_intervane"
        h5.attrs["voltage_convention"] = "V_intervane = V(+x vane) - V(+y vane)"
        axes = h5.create_group("axes")
        axes.create_dataset("x_m", data=x)
        axes.create_dataset("y_m", data=y)
        axes.create_dataset("z_m", data=z)
        fields = h5.create_group("fields")
        creation = dict(shape=shape, dtype="f4", chunks=chunks, compression="gzip", shuffle=True)
        phi_ds = fields.create_dataset("phi_per_v", **creation)
        ex_ds = fields.create_dataset("ex_vpm_per_v", **creation)
        ey_ds = fields.create_dataset("ey_vpm_per_v", **creation)
        ez_ds = fields.create_dataset("ez_vpm_per_v", **creation)

        for start in range(0, len(z), chunk_z):
            stop = min(start + chunk_z, len(z))
            halo_start = max(0, start - 1)
            halo_stop = min(len(z), stop + 1)
            local_z = z[halo_start:halo_stop]
            phi = _potential_per_intervane_volt(x, y, local_z, nodes)
            ex = -np.gradient(phi, x, axis=0, edge_order=2)
            ey = -np.gradient(phi, y, axis=1, edge_order=2)
            ez = -np.gradient(phi, local_z, axis=2, edge_order=2)
            local_start = start - halo_start
            local_stop = local_start + (stop - start)
            target = np.s_[:, :, start:stop]
            source = np.s_[:, :, local_start:local_stop]
            phi_ds[target] = phi[source].astype(np.float32)
            ex_ds[target] = ex[source].astype(np.float32)
            ey_ds[target] = ey[source].astype(np.float32)
            ez_ds[target] = ez[source].astype(np.float32)

    metadata = {
        "path": str(output_path),
        "shape": list(shape),
        "points": int(np.prod(shape)),
        "x_step_m": float(x[1] - x[0]),
        "y_step_m": float(y[1] - y[0]),
        "z_step_m": float(z[1] - z[0]),
        "normalization": "V/m per V_intervane",
        "size_bytes": output_path.stat().st_size,
    }
    output_path.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata

