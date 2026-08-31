from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = ("Cell", "V", "Wsyn", "A10", "Phi", "a", "m", "L", "Z")


def read_design_table(path: str | Path) -> pd.DataFrame:
    """Read a TRANSOPTR/PARMTEQ whitespace table without losing cell labels."""
    path = Path(path)
    frame = pd.read_csv(path, sep=r"\s+", engine="python", dtype={"Cell": str})
    frame.columns = [str(column).strip() for column in frame.columns]
    if "Cell" in frame:
        frame["Cell"] = frame["Cell"].astype(str).str.strip()
    for column in frame.columns:
        if column != "Cell":
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def validate_design(
    frame: pd.DataFrame,
    *,
    expected_cells: int | None = None,
    expected_length_m: float | None = None,
    length_tolerance_m: float = 1.0e-5,
) -> dict:
    missing = [name for name in REQUIRED_COLUMNS if name not in frame.columns]
    if missing:
        raise ValueError(f"Faltan columnas obligatorias: {missing}")
    if frame.empty:
        raise ValueError("La tabla de diseño está vacía.")
    if frame[list(REQUIRED_COLUMNS[1:])].isna().any().any():
        bad = frame[list(REQUIRED_COLUMNS[1:])].columns[
            frame[list(REQUIRED_COLUMNS[1:])].isna().any()
        ].tolist()
        raise ValueError(f"Hay valores no numéricos o vacíos en: {bad}")
    labels = frame["Cell"].str.extract(r"^(\d+)[A-Za-z]*$", expand=False)
    if labels.isna().any():
        raise ValueError("Todas las etiquetas Cell deben ser un entero con sufijo opcional.")
    cell_numbers = labels.astype(int).to_numpy()
    if np.any(np.diff(cell_numbers) <= 0):
        raise ValueError("Las celdas no están ordenadas o están repetidas.")
    if float(frame.iloc[0]["Z"]) != 0.0:
        raise ValueError("La fila Cell=0 debe comenzar en Z=0 cm.")
    physical = frame.iloc[1:]
    if (physical[["a", "m", "L"]] <= 0.0).any().any():
        raise ValueError("a, m y L deben ser positivos en las celdas físicas.")
    if np.any(np.diff(frame["Z"].to_numpy(float)) <= 0.0):
        raise ValueError("Z debe crecer estrictamente después de Cell=0.")
    if (physical["m"] < 1.0).any():
        raise ValueError("La modulación m no puede ser menor que 1.")

    physical_cells = len(physical)
    length_m = float(frame.iloc[-1]["Z"]) * 1.0e-2
    if expected_cells is not None and physical_cells != expected_cells:
        raise ValueError(
            f"Se esperaban {expected_cells} celdas físicas y se encontraron {physical_cells}."
        )
    if expected_length_m is not None and not np.isclose(
        length_m, expected_length_m, atol=length_tolerance_m, rtol=0.0
    ):
        raise ValueError(
            f"Longitud final {length_m:.8f} m; esperada {expected_length_m:.8f} m."
        )

    return {
        "rows": int(len(frame)),
        "physical_cells": int(physical_cells),
        "first_cell": str(frame.iloc[0]["Cell"]),
        "last_cell": str(frame.iloc[-1]["Cell"]),
        "length_m": length_m,
        "voltage_kv_min": float(frame["V"].min()),
        "voltage_kv_max": float(frame["V"].max()),
        "input_energy_mev": float(frame.iloc[0]["Wsyn"]),
        "output_energy_mev": float(frame.iloc[-1]["Wsyn"]),
    }


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def export_normalized_design(frame: pd.DataFrame, csv_path: Path, summary_path: Path, source: Path) -> dict:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    normalized = pd.DataFrame(
        {
            "cell": frame["Cell"],
            "voltage_kv": frame["V"],
            "synchronous_energy_mev": frame["Wsyn"],
            "synchronous_phase_deg": frame["Phi"],
            "aperture_m": frame["a"] * 1.0e-2,
            "modulation": frame["m"],
            "cell_length_m": frame["L"] * 1.0e-2,
            "z_end_m": frame["Z"] * 1.0e-2,
            "a10_transoptr": frame["A10"],
        }
    )
    normalized.to_csv(csv_path, index=False)
    summary = {
        "source": str(source),
        "source_sha256": sha256_file(source),
        "rows": int(len(frame)),
        "physical_cells": int(len(frame) - 1),
        "length_m": float(frame.iloc[-1]["Z"] * 1.0e-2),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary

