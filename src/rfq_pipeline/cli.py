from __future__ import annotations

import argparse
import json
from pathlib import Path

from .coefficients import derive_two_term_coefficients, write_coefficients
from .config import PipelineConfig, load_config
from .design import export_normalized_design, read_design_table, validate_design
from .diagnostics import compare_axis_fields, plot_field_diagnostics
from .fieldmap import generate_field_map
from .tracking import run_warp_tracking, tracking_plan
from .validation import validate_field_map
from .vanes import build_vane_profiles, write_vane_profiles


def _design(config: PipelineConfig):
    section = config.section("design")
    source = config.resolve(section["source"])
    frame = read_design_table(source)
    report = validate_design(
        frame,
        expected_cells=int(section["expected_physical_cells"]),
        expected_length_m=float(section["expected_length_m"]),
        length_tolerance_m=float(section.get("length_tolerance_m", 1.0e-5)),
    )
    return source, frame, report


def validate_design_stage(config: PipelineConfig) -> dict:
    source, frame, report = _design(config)
    output = config.output_dir
    output.mkdir(parents=True, exist_ok=True)
    provenance = export_normalized_design(
        frame, output / "design_si.csv", output / "design_provenance.json", source
    )
    result = {**report, **provenance, "passed": True}
    (output / "design_validation.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def coefficient_stage(config: PipelineConfig):
    _, frame, _ = _design(config)
    coefficients = derive_two_term_coefficients(frame)
    output = config.output_dir
    output.mkdir(parents=True, exist_ok=True)
    write_coefficients(coefficients, output / "coefficients.csv", output / "fort.75")
    return coefficients


def vane_stage(config: PipelineConfig) -> dict:
    _, frame, _ = _design(config)
    profiles = build_vane_profiles(frame, int(config.section("vanes")["points_per_cell"]))
    output = config.output_dir
    write_vane_profiles(profiles, output / "vane_profiles.csv", output / "vane_profiles.png")
    return {"points": int(len(profiles.z_m)), "length_m": float(profiles.z_m[-1])}


def field_stage(config: PipelineConfig) -> dict:
    coefficients = coefficient_stage(config)
    section = config.section("field")
    output_path = config.output_dir / "fieldmap.h5"
    return generate_field_map(
        coefficients,
        output_path,
        x_min_m=float(section["x_min_m"]),
        x_max_m=float(section["x_max_m"]),
        dx_m=float(section["dx_m"]),
        y_min_m=float(section["y_min_m"]),
        y_max_m=float(section["y_max_m"]),
        dy_m=float(section["dy_m"]),
        z_min_m=float(section["z_min_m"]),
        z_max_m=float(section["z_max_m"]),
        dz_m=float(section["dz_m"]),
        chunk_z=int(section.get("chunk_z", 64)),
    )


def validate_field_stage(config: PipelineConfig) -> dict:
    report = validate_field_map(
        config.output_dir / "fieldmap.h5",
        float(config.section("field").get("symmetry_tolerance", 5.0e-3)),
    )
    (config.output_dir / "field_validation.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def diagnostic_stage(config: PipelineConfig) -> dict:
    output = config.output_dir
    plot_field_diagnostics(output / "fieldmap.h5", output / "field_diagnostics.png")
    return {"plot": str(output / "field_diagnostics.png")}


def all_stage(config: PipelineConfig) -> dict:
    result = {
        "design": validate_design_stage(config),
        "coefficients": {"rows": int(len(coefficient_stage(config)))},
        "vanes": vane_stage(config),
        "field": field_stage(config),
    }
    result["field_validation"] = validate_field_stage(config)
    result["diagnostics"] = diagnostic_stage(config)
    (config.output_dir / "pipeline_summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Secuencia reproducible de diseño y simulación RFQ")
    parser.add_argument("--config", default="configs/isac2_quick.toml", help="Archivo TOML de configuración")
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("validate-design", help="Validar table1.txt y exportar unidades SI")
    commands.add_parser("coefficients", help="Derivar A01, A10, k y fase")
    commands.add_parser("vanes", help="Crear los perfiles de las cuatro vanes")
    commands.add_parser("fieldmap", help="Generar el mapa externo de dos términos")
    commands.add_parser("validate-field", help="Validar forma, finitud y simetrías")
    commands.add_parser("plot-field", help="Graficar el campo sobre el eje y el plano yz")
    track = commands.add_parser("track", help="Preparar o ejecutar tracking con Warp")
    track.add_argument("--dry-run", action="store_true", help="Validar parámetros sin importar Warp")
    compare = commands.add_parser("compare-fields", help="Comparar dos términos contra un eje exportado por Warp DC")
    compare.add_argument("warp_axis_csv", help="CSV con z_m y ez_vpm_per_v")
    commands.add_parser("all", help="Ejecutar diseño, coeficientes, vanes, campo y validaciones")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    dispatch = {
        "validate-design": validate_design_stage,
        "coefficients": coefficient_stage,
        "vanes": vane_stage,
        "fieldmap": field_stage,
        "validate-field": validate_field_stage,
        "plot-field": diagnostic_stage,
        "all": all_stage,
    }
    if args.command in dispatch:
        result = dispatch[args.command](config)
    elif args.command == "track":
        field_path = config.output_dir / "fieldmap.h5"
        if args.dry_run:
            result = tracking_plan(field_path, config.section("tracking"))
            (config.output_dir / "tracking_plan.json").write_text(
                json.dumps(result, indent=2) + "\n", encoding="utf-8"
            )
        else:
            result = run_warp_tracking(
                field_path, config.output_dir / "particles_final.npz", config.section("tracking")
            )
    else:
        result = compare_axis_fields(
            config.output_dir / "fieldmap.h5",
            args.warp_axis_csv,
            config.output_dir / "field_comparison.png",
            config.output_dir / "field_comparison.json",
        )
    if hasattr(result, "to_dict"):
        result = {"rows": int(len(result))}
    print(json.dumps(result, indent=2, default=str))
    return 0

