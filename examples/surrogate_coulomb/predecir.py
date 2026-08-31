#!/usr/bin/env python3
"""Usa el modelo entrenado para predecir una nueva configuración."""

from __future__ import annotations

import argparse

from prueba_ml.dataset import encode_parameters
from prueba_ml.distributions import GEOMETRIES, GeometryParameters
from prueba_ml.model import load_model_bundle


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", default="resultados/modelo_subrogado.joblib")
    parser.add_argument("--geometry", choices=GEOMETRIES, required=True)
    parser.add_argument("--charge-nc", type=float, required=True)
    parser.add_argument("--a-mm", type=float, required=True)
    parser.add_argument("--b-mm", type=float)
    parser.add_argument("--c-mm", type=float)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.geometry == "sphere":
        b_mm = c_mm = args.a_mm
    else:
        if args.b_mm is None or args.c_mm is None:
            raise SystemExit("gaussian y kv requieren --a-mm, --b-mm y --c-mm.")
        b_mm, c_mm = args.b_mm, args.c_mm
    params = GeometryParameters(args.geometry, args.charge_nc, args.a_mm, b_mm, c_mm)
    bundle = load_model_bundle(args.model)
    prediction = bundle["model"].predict(encode_parameters(params)[None, :])[0]
    names = bundle["dataset_metadata"]["output_names"]
    print(f"Predicción para {params}:")
    for name, value in zip(names, prediction):
        print(f"  {name:24s} {value:.8e}")


if __name__ == "__main__":
    main()

