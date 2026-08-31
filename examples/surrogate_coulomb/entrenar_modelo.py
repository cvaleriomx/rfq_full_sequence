#!/usr/bin/env python3
"""Entrena y evalúa la red neuronal subrogada."""

from __future__ import annotations

import argparse
from pathlib import Path

from prueba_ml.dataset import load_dataset
from prueba_ml.model import save_metrics, save_model_bundle, save_parity_plot, train_model


def parse_hidden_layers(value: str) -> tuple[int, ...]:
    try:
        layers = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Use una lista como 32,32,32,32.") from exc
    if not layers or min(layers) <= 0:
        raise argparse.ArgumentTypeError("Todas las capas deben ser positivas.")
    return layers


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, default=Path("datos/dataset_coulomb.npz"))
    parser.add_argument("--output", type=Path, default=Path("resultados/modelo_subrogado.joblib"))
    parser.add_argument("--metrics", type=Path, default=Path("resultados/metricas.json"))
    parser.add_argument("--plot", type=Path, default=Path("resultados/paridad.png"))
    parser.add_argument("--hidden", type=parse_hidden_layers, default=(32, 32, 32, 32))
    parser.add_argument("--max-iter", type=int, default=2000)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def run_training(args: argparse.Namespace) -> Path:
    inputs, outputs, metadata = load_dataset(args.dataset)
    print(f"Entrenando con X={inputs.shape}, y={outputs.shape}...")
    bundle, metrics, y_test, predictions = train_model(
        inputs,
        outputs,
        metadata,
        hidden_layers=args.hidden,
        seed=args.seed,
        max_iter=args.max_iter,
        test_fraction=args.test_fraction,
    )
    model_path = save_model_bundle(args.output, bundle)
    save_metrics(args.metrics, metrics)
    save_parity_plot(args.plot, y_test, predictions, metadata["output_names"])
    print(f"Modelo guardado en {model_path}")
    print(f"R² global: {metrics['r2_global_variance_weighted']:.6f}")
    print(f"Épocas: {metrics['epochs']}")
    for name, values in metrics["outputs"].items():
        print(
            f"  {name:24s} R²={values['r2']: .5f}  "
            f"NRMSE={100.0 * values['nrmse_range']:6.2f} %"
        )
    return model_path


def main() -> None:
    run_training(parse_args())


if __name__ == "__main__":
    main()

