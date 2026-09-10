#!/usr/bin/env python3
"""Iterative RFQ designer driven by TRANSOPTR feedback."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import shutil
import shlex
import subprocess
import sys
from typing import Any

import pandas as pd
import numpy as np

try:
    from .phase_control import field_nodes, tune_cell
    from . import sections
except ImportError:
    from phase_control import field_nodes, tune_cell
    import sections

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover - plotting is optional
    plt = None

try:
    import yaml
except Exception:  # pragma: no cover - optional dependency
    yaml = None

try:
    from .book_pages_23_32 import (
        PROTON_MASS_MEV,
        kt_coefficients_from_a_m_k,
        linear_schedule,
        longitudinal_phase_advance_deg,
        predicted_energy_gain_mev,
        rfq_cell_kinematics,
        transverse_phase_advance_deg,
    )
except ImportError:  # Allows: python rfq_transop_design/design_rfq.py design.yaml
    from book_pages_23_32 import (  # type: ignore
        PROTON_MASS_MEV,
        kt_coefficients_from_a_m_k,
        linear_schedule,
        longitudinal_phase_advance_deg,
        predicted_energy_gain_mev,
        rfq_cell_kinematics,
        transverse_phase_advance_deg,
    )


TABLE_COLUMNS = [
    "Cell",
    "V_kV",
    "Wsyn_MeV",
    "Sig0T_deg",
    "Sig0L_deg",
    "A10",
    "Phi_deg",
    "a_cm",
    "m",
    "B",
    "L_cm",
    "Z_cm",
    "A0",
    "RFdef",
    "Oct",
    "A1",
]


@dataclass
class CellDesign:
    cell: int
    V_kV: float
    Wsyn_MeV: float
    Wout_pred_MeV: float
    Wout_transoptr_MeV: float | None
    dW_pred_MeV: float
    dW_transoptr_MeV: float | None
    Sig0T_deg: float
    Sig0L_deg: float
    A01: float
    A10: float
    Phi_deg: float
    a_cm: float
    m: float
    B: float
    L_cm: float
    Z_cm: float
    k_cm_inv: float
    status: str

    phase_actual_deg: float | None = None
    phase_trials: int = 0
    section: str = "LEGACY"
    entry_a_cm: float | None = None
    phase_match: bool = True
    opening_start_cm: float | None = None
    opening_length_cm: float | None = None
    opening_a_initial_cm: float | None = None
    opening_a_final_cm: float | None = None

    def table_row(self) -> dict[str, Any]:
        return {
            "Cell": self.cell,
            "V_kV": self.V_kV,
            "Wsyn_MeV": self.Wsyn_MeV,
            "Sig0T_deg": self.Sig0T_deg,
            "Sig0L_deg": self.Sig0L_deg,
            "A10": self.A10,
            "Phi_deg": self.Phi_deg,
            "a_cm": self.a_cm,
            "m": self.m,
            "B": self.B,
            "L_cm": self.L_cm,
            "Z_cm": self.Z_cm,
            "A0": 1.0,
            "RFdef": 0.0,
            "Oct": 0.0,
            "A1": 0.0,
        }


def deep_update(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def default_config() -> dict[str, Any]:
    return {
        "output_dir": "rfq_design_output",
        "beam": {
            "longitudinal_full_width_deg": 90.0,
            "ion": "proton",
            "charge_state": 1.0,
            "mass_mev": PROTON_MASS_MEV,
            "energy_initial_mev": 0.030,
            "energy_target_mev": 0.080,
            "energy_tolerance_mev": 0.001,
        },
        "rf": {
            "frequency_hz": 162.0e6,
            "vane_voltage_kv": 35.2,
            "phi_start_deg": -90.0,
            "phi_end_deg": -30.0,
            "energy_gain_scale": 1.0,
        },
        "geometry": {
            "a_start_cm": 0.40,
            "a_end_cm": 0.30,
            "m_start": 1.0,
            "m_end": 1.8,
            "B_start": 3.0,
            "B_end": 4.0,
            "min_a_cm": 0.1,
            "max_a_cm": 2.0,
            "min_m": 1.0,
            "max_m": 3.0,
            "min_B": 0.0,
            "max_B": 10.0,
        },
        "transverse_tuning": {
            "enabled": False,
        },
        "transoptr": {
            "runner": "mock",
            "command": "./optr",
            "timeout_sec": 120,
            "template_data": None,
            "template_sy": None,
        },
        "iteration": {
            "max_cells": 50,
            "stop_policy": "energy_target",
            "accept_first_crossing": True,
        },
        "final_run": {
            "folder": "full_rfq_run",
            "copy_optr": False,
        },
    }


def load_config(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        user_config = json.loads(text)
    elif yaml is not None:
        user_config = yaml.safe_load(text) or {}
    else:
        user_config = parse_minimal_yaml(text)

    if not isinstance(user_config, dict):
        raise ValueError("configuration root must be a mapping")
    return deep_update(default_config(), user_config)


def parse_scalar(value: str) -> Any:
    value = value.strip()
    if value in {"", "null", "None", "~"}:
        return None
    if value in {"true", "True"}:
        return True
    if value in {"false", "False"}:
        return False
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    try:
        if any(ch in value for ch in [".", "e", "E"]):
            return float(value)
        return int(value)
    except ValueError:
        return value


def parse_minimal_yaml(text: str) -> dict[str, Any]:
    """Tiny indentation-based parser for this tool's simple config files."""
    root: dict[str, Any] = {}
    stack: list[tuple[int, dict[str, Any]]] = [(-1, root)]
    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        if ":" not in line:
            raise ValueError(f"invalid config line: {raw_line!r}")
        key, value = line.strip().split(":", 1)
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]
        if value.strip() == "":
            child: dict[str, Any] = {}
            parent[key] = child
            stack.append((indent, child))
        else:
            parent[key] = parse_scalar(value)
    return root


def require_positive(name: str, value: float) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive")


def validate_config(config: dict[str, Any]) -> None:
    beam = config["beam"]
    rf = config["rf"]
    iteration = config["iteration"]
    geometry = config["geometry"]
    if config.get("sections", {}).get("enabled", False):
        sections.validate_sections(config)
    require_positive("beam.mass_mev", float(beam["mass_mev"]))
    require_positive("rf.frequency_hz", float(rf["frequency_hz"]))
    require_positive("rf.vane_voltage_kv", float(rf["vane_voltage_kv"]))
    require_positive("iteration.max_cells", int(iteration["max_cells"]))
    if float(beam["energy_target_mev"]) <= float(beam["energy_initial_mev"]):
        raise ValueError("beam.energy_target_mev must be greater than beam.energy_initial_mev")
    if geometry["min_a_cm"] > geometry["max_a_cm"]:
        raise ValueError("geometry min/max aperture limits are inconsistent")
    if geometry["min_m"] > geometry["max_m"]:
        raise ValueError("geometry min/max modulation limits are inconsistent")
    if geometry["min_B"] > geometry["max_B"]:
        raise ValueError("geometry min/max B limits are inconsistent")


def scheduled_geometry(config: dict[str, Any], cell_index: int) -> tuple[float, float, float, float]:
    geometry = config["geometry"]
    rf = config["rf"]
    max_cells = int(config["iteration"]["max_cells"])
    a_cm = linear_schedule(cell_index, max_cells, float(geometry["a_start_cm"]), float(geometry["a_end_cm"]))
    m = linear_schedule(cell_index, max_cells, float(geometry["m_start"]), float(geometry["m_end"]))
    B = linear_schedule(cell_index, max_cells, float(geometry["B_start"]), float(geometry["B_end"]))
    phi = linear_schedule(cell_index, max_cells, float(rf["phi_start_deg"]), float(rf["phi_end_deg"]))

    a_cm = min(max(a_cm, float(geometry["min_a_cm"])), float(geometry["max_a_cm"]))
    m = min(max(m, float(geometry["min_m"])), float(geometry["max_m"]))
    B = min(max(B, float(geometry["min_B"])), float(geometry["max_B"]))
    return a_cm, m, B, phi


def propose_cell(config: dict[str, Any], cell_index: int, W_current: float, z_current_cm: float) -> CellDesign:
    beam = config["beam"]
    rf = config["rf"]
    a_cm, m, B, phi = scheduled_geometry(config, cell_index)
    kin = rfq_cell_kinematics(W_current, float(rf["frequency_hz"]), float(beam["mass_mev"]))
    section, entry_a = "LEGACY", None
    if config.get("sections", {}).get("enabled", False):
        a_cm, m, B, phi, section, entry_a = sections.section_parameters(config, cell_index, kin.k_cm_inv, kt_coefficients_from_a_m_k)
    A01, A10 = kt_coefficients_from_a_m_k(a_cm, m, kin.k_cm_inv)
    dW_pred = predicted_energy_gain_mev(
        vane_voltage_kv=float(rf["vane_voltage_kv"]),
        charge_state=float(beam["charge_state"]),
        A10=A10,
        phi_deg=phi,
        gain_scale=float(rf.get("energy_gain_scale", 1.0)),
    )
    sig0t = transverse_phase_advance_deg(B, a_cm, m, float(rf["vane_voltage_kv"]))
    sig0l = longitudinal_phase_advance_deg(A10, float(rf["vane_voltage_kv"]), phi)
    return CellDesign(
        cell=cell_index,
        V_kV=float(rf["vane_voltage_kv"]),
        Wsyn_MeV=W_current,
        Wout_pred_MeV=W_current + dW_pred,
        Wout_transoptr_MeV=None,
        dW_pred_MeV=dW_pred,
        dW_transoptr_MeV=None,
        Sig0T_deg=sig0t,
        Sig0L_deg=sig0l,
        A01=A01,
        A10=A10,
        Phi_deg=phi,
        a_cm=a_cm,
        m=m,
        B=B,
        L_cm=kin.cell_length_cm,
        Z_cm=z_current_cm + kin.cell_length_cm,
        k_cm_inv=kin.k_cm_inv,
        status="proposed",
        section=section, entry_a_cm=entry_a,
    )


def write_table1(cells: list[CellDesign], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write(" ".join(TABLE_COLUMNS) + "\n")
        for cell in cells:
            row = cell.table_row()
            f.write(
                f"{row['Cell']:d} "
                f"{row['V_kV']:.10g} {row['Wsyn_MeV']:.10g} "
                f"{row['Sig0T_deg']:.10g} {row['Sig0L_deg']:.10g} "
                f"{row['A10']:.10g} {row['Phi_deg']:.10g} "
                f"{row['a_cm']:.10g} {row['m']:.10g} {row['B']:.10g} "
                f"{row['L_cm']:.10g} {row['Z_cm']:.10g} "
                f"{row['A0']:.10g} {row['RFdef']:.10g} {row['Oct']:.10g} {row['A1']:.10g}\n"
            )


def write_pipeline_table(cells: list[CellDesign], path: Path) -> None:
    """Export SI-pipeline column names; Wsyn denotes cell exit energy.

    Energies are iterative feedback values, not a matched final design.
    Keep the original table and history for provenance.
    """
    if not cells:
        return
    rows = [cell.table_row() for cell in cells]
    for row, cell in zip(rows, cells):
        row["Wsyn_MeV"] = cell.Wout_transoptr_MeV
    initial = cells[0].table_row()
    initial.update(Cell=0, L_cm=0.0, Z_cm=0.0)
    if cells[0].entry_a_cm is not None:
        initial["a_cm"] = cells[0].entry_a_cm
    frame = pd.DataFrame([initial, *rows]).rename(columns={
        "V_kV": "V", "Wsyn_MeV": "Wsyn", "Sig0T_deg": "Sig0T",
        "Sig0L_deg": "Sig0L", "Phi_deg": "Phi", "a_cm": "a",
        "L_cm": "L", "Z_cm": "Z",
    })
    frame.to_csv(path, sep=" ", index=False)


def write_fort75(cells: list[CellDesign], path: Path, config=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if config and config.get("phase_control", {}).get("enabled", False):
        nodes = field_nodes(cells, int(config["phase_control"].get("field_points_per_cell", 16)))
        np.savetxt(path, nodes, fmt="%.12e")
        return
    with path.open("w", encoding="utf-8") as f:
        if not cells:
            return
        first = cells[0]
        f.write(
            f"{0.0:18.10E} {first.A01:18.10E} {first.A10:18.10E} {first.k_cm_inv:18.10E}"
            f"   ! cell={first.cell} start, a_cm={first.a_cm:.8f}, m={first.m:.8f}\n"
        )
        for cell in cells:
            f.write(
                f"{cell.Z_cm:18.10E} {cell.A01:18.10E} {cell.A10:18.10E} {cell.k_cm_inv:18.10E}"
                f"   ! cell={cell.cell}, L_cm={cell.L_cm:.8f}, a_cm={cell.a_cm:.8f}, m={cell.m:.8f}\n"
            )


def write_sy_f(cells: list[CellDesign], config: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    n_points = len(cells) + 1 if cells else 0
    if config.get("phase_control", {}).get("enabled", False):
        n_points = len(cells)*int(config["phase_control"].get("field_points_per_cell",16))+1
    length_cm = cells[-1].Z_cm if cells else 0.0
    freq_hz = float(config["rf"]["frequency_hz"])
    content = f"""      subroutine tsystem
c
c     Auto-generated by rfq_transop_design.
c     Edit the design YAML, not this generated file.
c
      implicit none
      integer nscav
      real vane, rfp, qsc, cmps
      integer isc
      common/scparm/qsc,isc,cmps
      cmps = {float(config.get("transoptr", {}).get("print_step_cm", 0.1)):.7E}
      nscav = 0
      vane = {float(config['rf']['vane_voltage_kv']) / 1000.0:.7E}
      rfp = {float(config["rf"].get("tracking_phase_deg", 0.0)):.7E}
      call rfq(75,{n_points},vane,{length_cm:.7E},
     $ {freq_hz:.7E},rfp,nscav)
      return
      end
"""
    path.write_text(content, encoding="utf-8")


def write_history(cells: list[CellDesign], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(asdict(cells[0]).keys()) if cells else [field.name for field in CellDesign.__dataclass_fields__.values()]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for cell in cells:
            writer.writerow(asdict(cell))


def write_plots(cells: list[CellDesign], output_dir: Path) -> None:
    if plt is None:
        return
    if not cells:
        return
    z = [c.Z_cm for c in cells]
    series = [
        ("energy_vs_z.png", "Wout_transoptr_MeV", "Energy [MeV]", [c.Wout_transoptr_MeV for c in cells]),
        ("m_vs_z.png", "m", "m", [c.m for c in cells]),
        ("a_vs_z.png", "a_cm", "a [cm]", [c.a_cm for c in cells]),
        ("B_vs_z.png", "B", "B", [c.B for c in cells]),
        ("phi_vs_z.png", "Phi_deg", "Phi [deg]", [c.Phi_deg for c in cells]),
    ]
    for filename, title, ylabel, values in series:
        clean_values = [v if v is not None else math.nan for v in values]
        plt.figure(figsize=(8, 4.8))
        plt.plot(z, clean_values, marker="o")
        plt.xlabel("z [cm]")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(output_dir / filename, dpi=160)
        plt.close()


def latest_envelope_path(cells: list[CellDesign], output_dir: Path) -> Path | None:
    if not cells:
        return None
    return output_dir / "iterations" / f"cell_{cells[-1].cell:04d}" / "fort.envelope"


def read_fort_envelope(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"missing fort.envelope: {path}")
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        columns = f.readline().strip().split()
    if not columns:
        raise ValueError(f"empty fort.envelope: {path}")
    return pd.read_csv(path, sep=r"\s+", skiprows=2, names=columns, engine="python")


def write_beam_envelope_outputs(cells: list[CellDesign], output_dir: Path) -> None:
    envelope_path = latest_envelope_path(cells, output_dir)
    if envelope_path is None or not envelope_path.exists():
        return

    try:
        df = read_fort_envelope(envelope_path)
    except Exception:
        return

    required = {"s", "x-envelope", "y-envelope"}
    if not required.issubset(df.columns):
        return

    out = pd.DataFrame(
        {
            "z_cm": pd.to_numeric(df["s"], errors="coerce"),
            "x_envelope_cm": pd.to_numeric(df["x-envelope"], errors="coerce"),
            "y_envelope_cm": pd.to_numeric(df["y-envelope"], errors="coerce"),
        }
    ).dropna()
    if out.empty:
        return

    out["abs_x_envelope_cm"] = out["x_envelope_cm"].abs()
    out["abs_y_envelope_cm"] = out["y_envelope_cm"].abs()
    out.to_csv(output_dir / "beam_envelope_xy.csv", index=False)

    if plt is None:
        return

    plt.figure(figsize=(8, 4.8))
    plt.plot(out["z_cm"], out["abs_x_envelope_cm"], label="|x-envelope|", linewidth=1.8)
    plt.plot(out["z_cm"], out["abs_y_envelope_cm"], label="|y-envelope|", linewidth=1.8)
    plt.xlabel("z [cm]")
    plt.ylabel("Envelope [cm]")
    plt.title("Beam envelope evolution")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "beam_envelope_xy_vs_z.png", dpi=160)
    plt.close()


def copy_optional_template(source: str | None, destination: Path) -> None:
    if not source:
        return
    src = Path(source).expanduser().resolve()
    if src.exists():
        shutil.copy2(src, destination)


def resolve_command_executable(command: str) -> Path | None:
    command_parts = shlex.split(command)
    if not command_parts:
        return None

    command_path = Path(command_parts[0]).expanduser()
    if not command_path.is_absolute():
        command_path = (Path.cwd() / command_path).resolve()

    if command_path.exists():
        return command_path
    return None


def write_final_run_folder(cells: list[CellDesign], config: dict[str, Any], output_dir: Path) -> None:
    if not cells:
        return

    final_run = config.get("final_run", {})
    run_dir = output_dir / str(final_run.get("folder", "full_rfq_run"))
    run_dir.mkdir(parents=True, exist_ok=True)

    write_fort75(cells, run_dir / "fort.75", config)
    write_sy_f(cells, config, run_dir / "sy.f")

    transoptr = config["transoptr"]
    latest_data = output_dir / "iterations" / f"cell_{cells[-1].cell:04d}" / "data.dat"
    template_data = transoptr.get("template_data")

    if latest_data.exists():
        shutil.copy2(latest_data, run_dir / "data.dat")
        data_source = str(latest_data)
    elif template_data:
        template_path = Path(str(template_data)).expanduser()
        if template_path.exists():
            shutil.copy2(template_path, run_dir / "data.dat")
            data_source = str(template_path)
        else:
            data_source = "missing template_data"
            write_missing_data_dat(run_dir / "data.dat", template_path)
    else:
        data_source = "not configured"
        write_missing_data_dat(run_dir / "data.dat", None)

    optr_source = None
    if bool(final_run.get("copy_optr", False)):
        optr_path = resolve_command_executable(str(transoptr.get("command", "./optr")))
        if optr_path is not None:
            shutil.copy2(optr_path, run_dir / optr_path.name)
            optr_source = str(optr_path)

    write_final_run_readme(run_dir, cells, config, data_source, optr_source)


def write_missing_data_dat(path: Path, template_path: Path | None) -> None:
    requested = str(template_path) if template_path is not None else "None"
    path.write_text(
        "! WARNING: data.dat template was not available when this folder was generated.\n"
        "! Configure transoptr.template_data in the YAML and rerun the designer.\n"
        f"! Requested template_data: {requested}\n",
        encoding="utf-8",
    )


def write_final_run_readme(
    run_dir: Path,
    cells: list[CellDesign],
    config: dict[str, Any],
    data_source: str,
    optr_source: str | None,
) -> None:
    final_energy = cells[-1].Wout_transoptr_MeV
    optr_hint = "./optr" if optr_source is not None else str(config["transoptr"].get("command", "./optr"))
    run_dir.joinpath("README.txt").write_text(
        "Full RFQ TRANSOPTR run folder\n"
        "=============================\n\n"
        "Files in this folder:\n"
        "  data.dat  - TRANSOPTR input deck/template for the complete RFQ run\n"
        "  sy.f      - generated RFQ call with the final point count and length\n"
        "  fort.75   - generated RFQ coefficient table for the complete RFQ\n\n"
        f"Cells: {len(cells)}\n"
        f"Final z: {cells[-1].Z_cm:.10E} cm\n"
        f"Final energy: {final_energy if final_energy is not None else float('nan'):.10E} MeV\n"
        f"data.dat source: {data_source}\n"
        f"optr source: {optr_source if optr_source is not None else 'not copied'}\n\n"
        "Typical use:\n"
        f"  cd {run_dir}\n"
        f"  {optr_hint}\n\n"
        "If data.dat only contains WARNING comments, set transoptr.template_data\n"
        "in the design YAML and rerun the designer.\n",
        encoding="utf-8",
    )


def parse_fort_envelope_energy(path: Path) -> float:
    df = read_fort_envelope(path)
    if "E" not in df.columns:
        raise ValueError(f"fort.envelope has no E column: {path}")
    return float(pd.to_numeric(df["E"], errors="coerce").dropna().iloc[-1])


def write_mock_envelope(cell: CellDesign, path: Path) -> None:
    path.write_text(
        "s E x-envelope y-envelope z-envelope\n"
        "cm MeV cm cm deg\n"
        f"0.0 {cell.Wsyn_MeV:.10E} 0.0 0.0 0.0\n"
        f"{cell.Z_cm:.10E} {cell.Wout_pred_MeV:.10E} 0.0 0.0 0.0\n",
        encoding="utf-8",
    )


def run_transoptr_feedback(cells: list[CellDesign], config: dict[str, Any], cell_dir: Path) -> tuple[float, str]:
    transoptr = config["transoptr"]
    runner = str(transoptr.get("runner", "mock")).lower()
    write_table1(cells, cell_dir / "table1.txt")
    write_fort75(cells, cell_dir / "fort.75", config)
    write_sy_f(cells, config, cell_dir / "sy.f")
    copy_optional_template(transoptr.get("template_data"), cell_dir / "data.dat")
    copy_optional_template(transoptr.get("template_sy"), cell_dir / "sy.template.f")
    if runner == "optr":
        data_path = cell_dir / "data.dat"
        if not data_path.exists():
            raise ValueError("Real TRANSOPTR requires template_data")
        lines = data_path.read_text().splitlines()
        beam_line = lines[0].split("!")[0].split()
        if len(beam_line) < 6:
            raise ValueError("Invalid TRANSOPTR beam input line")
        for index, key in ((0, "energy_initial_mev"), (3, "mass_mev"), (4, "charge_state")):
            beam_line[index] = str(float(config["beam"][key]))
        lines[0] = " ".join(beam_line) + " ! energy/mass/charge from design config"
        try:
            from .beam_input import prepare_beam_input
        except ImportError:
            from beam_input import prepare_beam_input
        lines, input_report = prepare_beam_input(lines, config)
        import json
        (cell_dir / 'beam_input.json').write_text(json.dumps(input_report, indent=2, allow_nan=False)+'\n')
        data_path.write_text("\n".join(lines)+"\n")
        # Never accept a stale phase/energy file from an earlier trial.
        for output_name in ("fort.envelope", "fort.20"):
            (cell_dir/output_name).unlink(missing_ok=True)

    if runner == "mock":
        write_mock_envelope(cells[-1], cell_dir / "fort.envelope")
        (cell_dir / "stdout.txt").write_text("mock TRANSOPTR runner\n", encoding="utf-8")
        return cells[-1].Wout_pred_MeV, "mock_ok"

    if runner != "optr":
        raise ValueError("transoptr.runner must be 'mock' or 'optr'")

    command = str(transoptr.get("command", "./optr"))
    command_parts = shlex.split(command)
    if not command_parts:
        raise ValueError("transoptr.command is empty")

    command_path = Path(command_parts[0]).expanduser()
    if not command_path.is_absolute():
        command_path = (Path.cwd() / command_path).resolve()

    if command_path.exists():
        run_command = [str(command_path), *command_parts[1:]]
    else:
        run_command = command_parts

    result = subprocess.run(
        run_command,
        cwd=cell_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=float(transoptr.get("timeout_sec", 120)),
    )
    (cell_dir / "stdout.txt").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        return cells[-1].Wsyn_MeV, f"optr_failed_{result.returncode}"
    return parse_fort_envelope_energy(cell_dir / "fort.envelope"), "optr_ok"


def run_design(config: dict[str, Any]) -> list[CellDesign]:
    validate_config(config)
    if config.get("phase_control", {}).get("enabled", False):
        if config["transoptr"]["runner"] != "optr":
            raise ValueError("Phase control requires real TRANSOPTR")
        config["rf"].setdefault("tracking_phase_deg", float(config["rf"]["phi_start_deg"]))
    output_dir = Path(config["output_dir"]).expanduser().resolve()
    iterations_dir = output_dir / "iterations"
    output_dir.mkdir(parents=True, exist_ok=True)
    iterations_dir.mkdir(parents=True, exist_ok=True)

    target = float(config["beam"]["energy_target_mev"])
    tolerance = float(config["beam"]["energy_tolerance_mev"])
    W_current = float(config["beam"]["energy_initial_mev"])
    z_current = 0.0
    cells: list[CellDesign] = []

    max_cells = int(config["iteration"]["max_cells"])
    sectioned = config.get("sections", {}).get("enabled", False)
    with_exit = sectioned and config.get("exit", {}).get("enabled", True)
    exit_start, last_acc = None, None
    opening_start_cm, opening_length_cm = None, None
    extra = sections.exit_count(config) if with_exit else 0
    for cell_index in range(1, max_cells + extra + 1):
        if exit_start is not None:
            ordinal = cell_index-exit_start
            kin = rfq_cell_kinematics(W_current, float(config["rf"]["frequency_hz"]), float(config["beam"]["mass_mev"]))
            opening = int(config.get('exit', {}).get('opening_cells', 0))
            ntc = int(config.get('exit', {}).get('transition_cells', 1))
            if opening and ordinal == ntc+1:
                opening_start_cm = z_current
                opening_length_cm = opening*kin.cell_length_cm
            a, m, B, phi, length, label = sections.exit_parameters(config, ordinal, last_acc, W_current, kin, kt_coefficients_from_a_m_k, opening_length_cm)
            from dataclasses import replace
            cell = replace(last_acc, cell=cell_index, Wsyn_MeV=W_current, a_cm=a, m=m, B=B, Phi_deg=phi,
                L_cm=length, Z_cm=z_current+length, k_cm_inv=math.pi/length, section=label,
                phase_match=False, phase_actual_deg=None, phase_trials=0, status="proposed")
            cell.A01, cell.A10 = kt_coefficients_from_a_m_k(a,m,cell.k_cm_inv)
            if label == 'OM':
                cell.opening_start_cm = opening_start_cm
                cell.opening_length_cm = opening_length_cm
                cell.opening_a_initial_cm = math.sqrt(sections.focusing_scale(config)/last_acc.B)
                cell.opening_a_final_cm = sections.exit_aperture(config)
            cell.Wout_pred_MeV, cell.dW_pred_MeV = W_current, 0.
        else:
            cell = propose_cell(config, cell_index, W_current, z_current)
        trial_cells = cells + [cell]
        cell_dir = iterations_dir / f"cell_{cell_index:04d}"
        cell_dir.mkdir(parents=True, exist_ok=True)

        try:
            if config.get("phase_control", {}).get("enabled", False) and cell.phase_match:
                def update_geometry(c, length):
                    c.L_cm = length
                    c.Z_cm = z_current + length
                    c.k_cm_inv = math.pi/length
                    if sectioned:
                        c.a_cm = sections.aperture_for_B(c.B,c.m,c.k_cm_inv,config,kt_coefficients_from_a_m_k)
                    c.A01, c.A10 = kt_coefficients_from_a_m_k(c.a_cm, c.m, c.k_cm_inv)
                    c.dW_pred_MeV = predicted_energy_gain_mev(c.V_kV, float(config["beam"]["charge_state"]), c.A10, c.Phi_deg)
                    c.Wout_pred_MeV = c.Wsyn_MeV + c.dW_pred_MeV
                W_transoptr, run_status = tune_cell(trial_cells, config, cell_dir, run_transoptr_feedback, update_geometry)
            else:
                W_transoptr, run_status = run_transoptr_feedback(trial_cells, config, cell_dir)
        except subprocess.TimeoutExpired:
            W_transoptr, run_status = W_current, "optr_timeout"
        except Exception as exc:
            W_transoptr, run_status = W_current, f"runner_error:{exc}"

        if run_status not in {"mock_ok", "optr_ok"}:
            raise RuntimeError(f"Cell {cell_index}: {run_status}; inspect {cell_dir}")
        if not math.isfinite(W_transoptr) or W_transoptr <= 0:
            raise ValueError(f"Invalid TRANSOPTR energy: {W_transoptr}")

        cell.Wout_transoptr_MeV = W_transoptr
        cell.dW_transoptr_MeV = W_transoptr - W_current
        cell.status = run_status
        cells.append(cell)

        W_current = W_transoptr
        z_current = cell.Z_cm

        # Keep history live; render plots once, after the iterative solve.
        write_history(cells, output_dir / "history.csv")
        print(f"Cell {cell_index}: E={W_current:.7f} MeV, phase={cell.phase_actual_deg}", flush=True)

        if exit_start is not None:
            if cell_index-exit_start == extra:
                cell.status = "target_reached" if abs(W_current-target) <= tolerance else "exit_energy_outside_tolerance"
                break
            continue
        can_stop = not sectioned or cell.section in {"ACC", "MA"}
        # The final energy is always checked after the complete exit.
        reserve = float(config.get('exit', {}).get('energy_reserve_mev', 0.)) if with_exit else 0.
        if can_stop and W_current >= target - tolerance - reserve:
            if with_exit:
                exit_start, last_acc = cell_index, cell
                continue
            cell.status = "target_reached" if abs(W_current-target) <= tolerance else "crossed_target_outside_tolerance"
            break
        if cell_index >= max_cells:
            if with_exit:
                exit_start, last_acc = cell_index, cell
                continue
            cell.status = "max_cells_reached"
            break

    write_outputs(cells, config, output_dir)
    return cells


def write_outputs(cells: list[CellDesign], config: dict[str, Any], output_dir: Path) -> None:
    write_table1(cells, output_dir / "table1.txt")
    write_pipeline_table(cells, output_dir / "table_pipeline.txt")
    write_fort75(cells, output_dir / "fort.75", config)
    write_sy_f(cells, config, output_dir / "sy.f")
    write_history(cells, output_dir / "history.csv")
    write_plots(cells, output_dir)
    write_beam_envelope_outputs(cells, output_dir)
    write_final_run_folder(cells, config, output_dir)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Iterative RFQ designer with TRANSOPTR feedback.")
    parser.add_argument("config", type=Path, help="YAML or JSON design configuration")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    cells = run_design(config)
    if not cells:
        print("No cells generated.", file=sys.stderr)
        return 1

    if config.get("phase_control", {}).get("enabled", False):
        try:
            from .verify_phase import verify_final
        except ImportError:
            from verify_phase import verify_final
        report = verify_final(cells, config, Path(config["output_dir"]), run_transoptr_feedback, read_fort_envelope)
        if not report["phase_and_refinement_passed"]:
            raise RuntimeError("Final phase/refinement validation failed; inspect phase_validation.json")

    last = cells[-1]
    print(f"Generated {len(cells)} cells")
    print(f"Final energy: {last.Wout_transoptr_MeV:.8f} MeV")
    print(f"Final status: {last.status}")
    print(f"Output dir: {Path(config['output_dir']).expanduser().resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
