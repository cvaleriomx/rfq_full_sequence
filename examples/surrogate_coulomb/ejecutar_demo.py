#!/usr/bin/env python3
"""Ejecuta generación, entrenamiento y una predicción de ejemplo."""

from __future__ import annotations

import argparse
from pathlib import Path

from prueba_ml.dataset import encode_parameters, generate_dataset, save_dataset
from prueba_ml.distributions import GeometryParameters
from prueba_ml.model import save_metrics, save_model_bundle, save_parity_plot, train_model


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="Usa 500 casos y 10 000 partículas; sin esta opción ejecuta una prueba rápida.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("resultados/demo"))
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.full:
        n_samples, n_particles, max_iter = 500, 20_000, 2000
    else:
        n_samples, n_particles, max_iter = 300, 2000, 1500

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"1/3 Generando {n_samples} casos con N={n_particles}...")
    inputs, outputs, metadata = generate_dataset(n_samples, n_particles, seed=args.seed)
    dataset_path = save_dataset(
        args.output_dir / "dataset_coulomb.npz", inputs, outputs, metadata
    )

    print("2/3 Entrenando el modelo...")
    bundle, metrics, y_test, predictions = train_model(
        inputs, outputs, metadata, seed=args.seed, max_iter=max_iter
    )
    model_path = save_model_bundle(args.output_dir / "modelo_subrogado.joblib", bundle)
    save_metrics(args.output_dir / "metricas.json", metrics)
    save_parity_plot(
        args.output_dir / "paridad.png", y_test, predictions, metadata["output_names"]
    )

    print("3/3 Predicción de una gaussiana nueva...")
    example = GeometryParameters("gaussian", 3.0, 2.5, 3.5, 5.0)
    predicted = bundle["model"].predict(encode_parameters(example)[None, :])[0]
    for name, value in zip(metadata["output_names"], predicted):
        print(f"  {name:24s} {value:.7e}")
    print(f"Dataset: {dataset_path}")
    print(f"Modelo:  {model_path}")
    print(f"R² global de prueba: {metrics['r2_global_variance_weighted']:.6f}")


if __name__ == "__main__":
    main()
