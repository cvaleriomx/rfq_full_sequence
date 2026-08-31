#!/usr/bin/env python3
"""Genera observaciones de Coulomb para entrenar el modelo subrogado."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from prueba_ml.dataset import generate_dataset, save_dataset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=500, help="Número de configuraciones.")
    parser.add_argument(
        "--particles", type=int, default=10_000, help="Macropartículas; debe ser múltiplo de 8."
    )
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--charge-min-nc", type=float, default=0.5)
    parser.add_argument("--charge-max-nc", type=float, default=10.0)
    parser.add_argument("--axis-min-mm", type=float, default=1.0)
    parser.add_argument("--axis-max-mm", type=float, default=10.0)
    parser.add_argument("--softening-fraction", type=float, default=1.25)
    parser.add_argument(
        "--include-energy",
        action="store_true",
        help="Añade la energía de todos los pares; su costo es O(N²).",
    )
    parser.add_argument("--output", type=Path, default=Path("datos/dataset_coulomb.npz"))
    return parser.parse_args()


def run_generation(args: argparse.Namespace) -> Path:
    start = time.perf_counter()
    last_percent = -10

    def report(done: int, total: int) -> None:
        nonlocal last_percent
        percent = int(100 * done / total)
        if percent >= last_percent + 10 or done == total:
            print(f"  {done:5d}/{total} configuraciones ({percent:3d} %)")
            last_percent = percent

    print(
        f"Generando {args.samples} configuraciones con {args.particles} macropartículas..."
    )
    inputs, outputs, metadata = generate_dataset(
        n_samples=args.samples,
        n_particles=args.particles,
        seed=args.seed,
        charge_range_nc=(args.charge_min_nc, args.charge_max_nc),
        axis_range_mm=(args.axis_min_mm, args.axis_max_mm),
        softening_fraction=args.softening_fraction,
        include_energy=args.include_energy,
        progress=report,
    )
    destination = save_dataset(args.output, inputs, outputs, metadata)
    elapsed = time.perf_counter() - start
    print(f"Dataset guardado en {destination} ({elapsed:.2f} s).")
    print(f"Forma de X: {inputs.shape}; forma de y: {outputs.shape}")
    return destination


def main() -> None:
    run_generation(parse_args())


if __name__ == "__main__":
    main()
