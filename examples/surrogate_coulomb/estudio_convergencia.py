#!/usr/bin/env python3
"""Estudia la convergencia con N contra soluciones continuas conocidas."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np

from prueba_ml.dataset import RADIAL_FACTORS, evaluate_configuration
from prueba_ml.distributions import GEOMETRIES, GeometryParameters, sample_unit_cloud
from prueba_ml.physics import (
    gaussian_center_potential,
    gaussian_sphere_field,
    uniform_sphere_center_potential,
    uniform_sphere_field,
)


DEFAULT_PARAMETERS = {
    "sphere": GeometryParameters("sphere", 1.0, 5.0, 5.0, 5.0),
    "gaussian": GeometryParameters("gaussian", 1.0, 3.0, 3.0, 3.0),
    "kv": GeometryParameters("kv", 1.0, 4.0, 5.0, 6.0),
}


def parse_counts(value: str) -> tuple[int, ...]:
    try:
        counts = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Use una lista como 1024,2048,5000,10000.") from exc
    if len(counts) < 2 or min(counts) <= 0 or any(count % 8 for count in counts):
        raise argparse.ArgumentTypeError(
            "Use al menos dos valores positivos, todos múltiplos de 8."
        )
    return tuple(sorted(set(counts)))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--particles", type=parse_counts, default=(1024, 2048, 5000, 10000, 20000)
    )
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--softening-fraction", type=float, default=1.25)
    parser.add_argument("--output-dir", type=Path, default=Path("resultados/convergencia"))
    return parser.parse_args()


def analytic_reference(params: GeometryParameters) -> np.ndarray | None:
    radii = np.asarray(RADIAL_FACTORS) * params.axes_m[0]
    if params.geometry == "sphere":
        fields = uniform_sphere_field(radii, params.total_charge_c, params.axes_m[0])
        potential = uniform_sphere_center_potential(params.total_charge_c, params.axes_m[0])
    elif params.geometry == "gaussian":
        fields = gaussian_sphere_field(radii, params.total_charge_c, params.axes_m[0])
        potential = gaussian_center_potential(params.total_charge_c, params.axes_m[0])
    else:
        return None
    return np.concatenate((np.tile(fields, 3), (potential,)))


def run_study(args: argparse.Namespace) -> list[dict]:
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    all_outputs: dict[str, dict[int, np.ndarray]] = {}

    for geometry_index, geometry in enumerate(GEOMETRIES):
        params = DEFAULT_PARAMETERS[geometry]
        all_outputs[geometry] = {}
        for n_particles in args.particles:
            cloud = sample_unit_cloud(
                geometry,
                n_particles,
                np.random.default_rng(args.seed + 1000 * geometry_index),
                quiet_start=True,
            )
            all_outputs[geometry][n_particles] = evaluate_configuration(
                params, cloud, softening_fraction=args.softening_fraction
            )

    rows: list[dict] = []
    for geometry in GEOMETRIES:
        analytic = analytic_reference(DEFAULT_PARAMETERS[geometry])
        if analytic is None:
            reference = all_outputs[geometry][max(args.particles)]
            reference_name = f"N={max(args.particles)}"
        else:
            reference = analytic
            reference_name = "analítica continua"
        denominator = np.linalg.norm(reference)
        for n_particles in args.particles:
            value = all_outputs[geometry][n_particles]
            error = float(np.linalg.norm(value - reference) / denominator)
            rows.append(
                {
                    "geometry": geometry,
                    "n_particles": n_particles,
                    "relative_l2_error": error,
                    "reference": reference_name,
                }
            )

    csv_path = output_dir / "convergencia.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    mpl_cache = output_dir / ".mplconfig"
    mpl_cache.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(7.2, 4.8))
    for geometry in GEOMETRIES:
        selected = [row for row in rows if row["geometry"] == geometry]
        x = np.asarray([row["n_particles"] for row in selected])
        y = np.asarray([row["relative_l2_error"] for row in selected])
        positive = y > 0.0
        axis.loglog(x[positive], y[positive], "o-", label=geometry)
    axis.axhline(0.01, color="black", linestyle="--", linewidth=1, label="1 %")
    axis.set_xlabel("Número de macropartículas N")
    axis.set_ylabel("Error relativo L2")
    axis.set_title("Convergencia de los observables de Coulomb")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_dir / "convergencia.png", dpi=170)
    plt.close(figure)

    print(f"Resultados guardados en {output_dir}")
    for row in rows:
        print(
            f"  {row['geometry']:8s} N={row['n_particles']:6d}  "
            f"error={100.0 * row['relative_l2_error']:7.3f} %  "
            f"({row['reference']})"
        )
    return rows


def main() -> None:
    run_study(parse_args())


if __name__ == "__main__":
    main()
