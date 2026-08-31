#!/usr/bin/env python3
"""Grafica E_y(y) usando únicamente un modelo subrogado entrenado."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np

from prueba_ml.dataset import encode_parameters
from prueba_ml.distributions import GEOMETRIES, GeometryParameters
from prueba_ml.model import load_model_bundle


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model", type=Path, default=Path("resultados/demo/modelo_subrogado.joblib")
    )
    parser.add_argument("--geometry", choices=GEOMETRIES, required=True)
    parser.add_argument("--charge-nc", type=float, required=True)
    parser.add_argument("--a-mm", type=float, required=True)
    parser.add_argument(
        "--b-mm",
        type=float,
        help="Sigma_y o semieje y. Para sphere se usa automáticamente a_mm.",
    )
    parser.add_argument(
        "--c-mm",
        type=float,
        help="Sigma_z o semieje z. Para sphere se usa automáticamente a_mm.",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("resultados/campo_eje_y.png")
    )
    parser.add_argument(
        "--csv", type=Path, default=Path("resultados/campo_eje_y.csv")
    )
    parser.add_argument("--show", action="store_true", help="Muestra además la ventana de Matplotlib.")
    return parser.parse_args()


def build_parameters(args: argparse.Namespace) -> GeometryParameters:
    if args.geometry == "sphere":
        b_mm = c_mm = args.a_mm
    else:
        if args.b_mm is None or args.c_mm is None:
            raise SystemExit("gaussian y kv requieren --a-mm, --b-mm y --c-mm.")
        b_mm, c_mm = args.b_mm, args.c_mm
    return GeometryParameters(args.geometry, args.charge_nc, args.a_mm, b_mm, c_mm)


def main() -> None:
    args = parse_args()
    params = build_parameters(args)
    bundle = load_model_bundle(args.model)
    output_names = bundle["dataset_metadata"]["output_names"]
    required = ("E_y_0p5_V_m", "E_y_1p0_V_m", "E_y_2p0_V_m")
    missing = [name for name in required if name not in output_names]
    if missing:
        raise SystemExit(f"El modelo no contiene las salidas requeridas: {missing}")

    prediction = bundle["model"].predict(encode_parameters(params)[None, :])[0]
    positive_field = np.asarray([prediction[output_names.index(name)] for name in required])
    positive_y_mm = params.b_mm * np.asarray((0.5, 1.0, 2.0))

    # Las tres distribuciones implementadas están centradas y son simétricas:
    # E_y(-y)=-E_y(y) y E_y(0)=0.
    y_mm = np.concatenate((-positive_y_mm[::-1], (0.0,), positive_y_mm))
    field_y = np.concatenate((-positive_field[::-1], (0.0,), positive_field))

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("y_mm", "E_y_V_m", "source"))
        for y_value, field_value in zip(y_mm, field_y):
            source = "symmetry" if y_value < 0.0 else ("center" if y_value == 0.0 else "surrogate")
            writer.writerow((f"{y_value:.12g}", f"{field_value:.12g}", source))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    mpl_cache = args.output.parent / ".mplconfig"
    mpl_cache.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(7.4, 4.8))
    axis.plot(y_mm, field_y, "--", color="tab:blue", alpha=0.65, label="Guía entre puntos")
    axis.scatter(
        positive_y_mm,
        positive_field,
        s=55,
        color="tab:blue",
        zorder=3,
        label="Predicción del subrogado",
    )
    axis.scatter(
        -positive_y_mm,
        -positive_field,
        s=45,
        facecolors="white",
        edgecolors="tab:blue",
        zorder=3,
        label="Reflexión por simetría",
    )
    axis.scatter((0.0,), (0.0,), s=45, color="black", zorder=3)
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.axvline(0.0, color="black", linewidth=0.8)
    axis.set_xlabel("Posición y [mm]")
    axis.set_ylabel(r"Campo eléctrico $E_y$ [V/m]")
    axis.set_title(
        f"Campo sobre el eje y — {params.geometry}, Q={params.charge_nc:g} nC"
    )
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(args.output, dpi=180)
    if args.show:
        plt.show()
    plt.close(figure)

    print(f"Modelo:  {args.model}")
    print(f"Gráfica: {args.output}")
    print(f"Datos:   {args.csv}")
    print("Puntos calculados directamente por el subrogado:")
    for position, value in zip(positive_y_mm, positive_field):
        print(f"  y={position:10.5f} mm    E_y={value: .8e} V/m")
    print("La línea discontinua sólo conecta los puntos; no es una predicción continua.")


if __name__ == "__main__":
    main()

