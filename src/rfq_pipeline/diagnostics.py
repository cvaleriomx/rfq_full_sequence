from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def plot_field_diagnostics(field_path: str | Path, output_path: str | Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with h5py.File(field_path, "r") as h5:
        x = h5["axes/x_m"][:]
        y = h5["axes/y_m"][:]
        z = h5["axes/z_m"][:]
        ix = int(np.argmin(np.abs(x)))
        iy = int(np.argmin(np.abs(y)))
        ez_axis = h5["fields/ez_vpm_per_v"][ix, iy, :]
        ex_axis = h5["fields/ex_vpm_per_v"][ix, iy, :]
        ez_yz = h5["fields/ez_vpm_per_v"][ix, :, :]

    figure, axes = plt.subplots(2, 1, figsize=(12, 7), constrained_layout=True)
    axes[0].plot(z, ez_axis, lw=0.8, label="$E_z(0,0,z)$")
    axes[0].plot(z, ex_axis, lw=0.8, label="$E_x(0,0,z)$")
    axes[0].set(xlabel="z [m]", ylabel="E [V/m/V intervane]", title="Campo sobre el eje")
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    image = axes[1].pcolormesh(z, y * 1.0e3, ez_yz, shading="auto")
    axes[1].set(xlabel="z [m]", ylabel="y [mm]", title="$E_z$ en el plano x=0")
    figure.colorbar(image, ax=axes[1], label="V/m/V intervane")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def compare_axis_fields(
    external_field_path: str | Path,
    warp_axis_csv: str | Path,
    plot_path: str | Path,
    metrics_path: str | Path,
) -> dict:
    """Compare external two-term Ez against a Warp DC axis export."""
    with h5py.File(external_field_path, "r") as h5:
        x = h5["axes/x_m"][:]
        y = h5["axes/y_m"][:]
        z = h5["axes/z_m"][:]
        external = h5["fields/ez_vpm_per_v"][np.argmin(abs(x)), np.argmin(abs(y)), :]
    warp = pd.read_csv(warp_axis_csv)
    required = {"z_m", "ez_vpm_per_v"}
    if not required.issubset(warp.columns):
        raise ValueError(f"El CSV de Warp necesita las columnas {sorted(required)}.")
    warp_interp = np.interp(z, warp["z_m"], warp["ez_vpm_per_v"])
    scale_external = max(float(np.max(np.abs(external))), 1.0e-30)
    scale_warp = max(float(np.max(np.abs(warp_interp))), 1.0e-30)
    a = external / scale_external
    b = warp_interp / scale_warp
    metrics = {
        "normalized_rmse": float(np.sqrt(np.mean((a - b) ** 2))),
        "correlation": float(np.corrcoef(a, b)[0, 1]),
        "external_peak_vpm_per_v": scale_external,
        "warp_peak_vpm_per_v": scale_warp,
    }
    Path(metrics_path).write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(12, 4))
    axis.plot(z, a, label="dos términos externo", lw=0.8)
    axis.plot(z, b, label="Warp DC", lw=0.8, alpha=0.8)
    axis.set(xlabel="z [m]", ylabel="$E_z$ normalizado", title="Comparación de campo sobre el eje")
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)
    return metrics

