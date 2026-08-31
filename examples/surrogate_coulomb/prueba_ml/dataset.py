"""Muestreo de parámetros y creación de datos para el modelo subrogado."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

from .distributions import GEOMETRIES, GeometryParameters, sample_unit_cloud, scale_unit_cloud
from .physics import electric_field, electric_potential, electrostatic_energy


INPUT_NAMES = (
    "is_sphere",
    "is_gaussian",
    "is_kv",
    "charge_nC",
    "a_mm",
    "b_mm",
    "c_mm",
)
RADIAL_FACTORS = (0.5, 1.0, 2.0)
FIELD_OUTPUT_NAMES = tuple(
    f"E_{axis}_{str(factor).replace('.', 'p')}_V_m"
    for axis in "xyz"
    for factor in RADIAL_FACTORS
)
BASE_OUTPUT_NAMES = FIELD_OUTPUT_NAMES + ("V_center_V",)


def encode_parameters(params: GeometryParameters) -> np.ndarray:
    one_hot = [float(params.geometry == name) for name in GEOMETRIES]
    return np.asarray(one_hot + [params.charge_nc, params.a_mm, params.b_mm, params.c_mm])


def decode_parameters(features: Iterable[float]) -> GeometryParameters:
    vector = np.asarray(tuple(features), dtype=float)
    if vector.shape != (len(INPUT_NAMES),):
        raise ValueError(f"Se esperaban {len(INPUT_NAMES)} entradas.")
    geometry = GEOMETRIES[int(np.argmax(vector[:3]))]
    a_mm, b_mm, c_mm = vector[4:7]
    if geometry == "sphere":
        b_mm = c_mm = a_mm
    return GeometryParameters(geometry, vector[3], a_mm, b_mm, c_mm)


def latin_hypercube(n_points: int, n_dimensions: int, rng: np.random.Generator) -> np.ndarray:
    """Latin hypercube sencillo sin depender de scipy."""

    if n_points <= 0 or n_dimensions <= 0:
        raise ValueError("Las dimensiones del Latin hypercube deben ser positivas.")
    sample = np.empty((n_points, n_dimensions), dtype=float)
    for dimension in range(n_dimensions):
        sample[:, dimension] = (rng.permutation(n_points) + rng.random(n_points)) / n_points
    return sample


def sample_parameters(
    n_samples: int,
    rng: np.random.Generator,
    charge_range_nc: tuple[float, float] = (0.5, 10.0),
    axis_range_mm: tuple[float, float] = (1.0, 10.0),
) -> list[GeometryParameters]:
    """Muestra de forma estratificada las tres familias geométricas."""

    if n_samples < len(GEOMETRIES):
        raise ValueError(f"n_samples debe ser al menos {len(GEOMETRIES)}.")
    if not (0.0 < charge_range_nc[0] < charge_range_nc[1]):
        raise ValueError("charge_range_nc debe ser positivo y creciente.")
    if not (0.0 < axis_range_mm[0] < axis_range_mm[1]):
        raise ValueError("axis_range_mm debe ser positivo y creciente.")

    counts = [n_samples // len(GEOMETRIES)] * len(GEOMETRIES)
    for index in range(n_samples % len(GEOMETRIES)):
        counts[index] += 1

    result: list[GeometryParameters] = []
    q_min, q_max = charge_range_nc
    a_min, a_max = axis_range_mm
    for geometry, count in zip(GEOMETRIES, counts):
        cube = latin_hypercube(count, 4, rng)
        charge = q_min + (q_max - q_min) * cube[:, 0]
        axes = a_min + (a_max - a_min) * cube[:, 1:4]
        if geometry == "sphere":
            axes[:, 1:] = axes[:, [0]]
        for row in range(count):
            result.append(
                GeometryParameters(geometry, charge[row], axes[row, 0], axes[row, 1], axes[row, 2])
            )
    rng.shuffle(result)
    return result


def measurement_probes(params: GeometryParameters) -> np.ndarray:
    """Sondas sobre cada semieje, escaladas con el tamaño de la distribución."""

    probes = []
    for axis, length in enumerate(params.axes_m):
        for factor in RADIAL_FACTORS:
            point = np.zeros(3)
            point[axis] = factor * length
            probes.append(point)
    return np.asarray(probes)


def softening_length(params: GeometryParameters, n_particles: int, fraction: float) -> float:
    """Suavizado proporcional al espaciamiento medio característico."""

    characteristic_length = float(np.prod(params.axes_m) ** (1.0 / 3.0))
    return fraction * characteristic_length / n_particles ** (1.0 / 3.0)


def evaluate_configuration(
    params: GeometryParameters,
    unit_cloud: np.ndarray,
    softening_fraction: float = 1.25,
    block_size: int = 4096,
    include_energy: bool = False,
    energy_block_size: int = 1024,
) -> np.ndarray:
    """Calcula las etiquetas físicas de una configuración."""

    positions = scale_unit_cloud(unit_cloud, params)
    n_particles = len(positions)
    macro_charge = params.total_charge_c / n_particles
    softening = softening_length(params, n_particles, softening_fraction)
    probes = measurement_probes(params)
    field = electric_field(probes, positions, macro_charge, softening, block_size)

    axial_components = []
    for axis in range(3):
        first = axis * len(RADIAL_FACTORS)
        axial_components.extend(field[first : first + len(RADIAL_FACTORS), axis])
    potential_center = electric_potential(
        np.zeros((1, 3)), positions, macro_charge, softening, block_size
    )[0]
    outputs = axial_components + [potential_center]
    if include_energy:
        outputs.append(
            electrostatic_energy(positions, macro_charge, softening, energy_block_size)
        )
    return np.asarray(outputs, dtype=float)


def generate_dataset(
    n_samples: int,
    n_particles: int,
    seed: int = 2026,
    charge_range_nc: tuple[float, float] = (0.5, 10.0),
    axis_range_mm: tuple[float, float] = (1.0, 10.0),
    softening_fraction: float = 1.25,
    include_energy: bool = False,
    progress: Callable[[int, int], None] | None = None,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Genera un conjunto determinista entrada/salida.

    Se reutiliza una nube unidad por geometría (*common random numbers*). Así
    la red aprende la variación física con los parámetros y no una semilla
    aleatoria distinta en cada fila.
    """

    if n_particles % 8:
        raise ValueError("n_particles debe ser múltiplo de 8 para usar quiet start.")
    rng = np.random.default_rng(seed)
    parameters = sample_parameters(n_samples, rng, charge_range_nc, axis_range_mm)
    clouds = {
        geometry: sample_unit_cloud(
            geometry, n_particles, np.random.default_rng(seed + 101 * (index + 1)), True
        )
        for index, geometry in enumerate(GEOMETRIES)
    }

    output_names = BASE_OUTPUT_NAMES + (("U_electrostatic_J",) if include_energy else ())
    inputs = np.empty((n_samples, len(INPUT_NAMES)), dtype=float)
    outputs = np.empty((n_samples, len(output_names)), dtype=float)
    for index, params in enumerate(parameters):
        inputs[index] = encode_parameters(params)
        outputs[index] = evaluate_configuration(
            params,
            clouds[params.geometry],
            softening_fraction=softening_fraction,
            include_energy=include_energy,
        )
        if progress is not None:
            progress(index + 1, n_samples)

    metadata = {
        "seed": seed,
        "n_samples": n_samples,
        "n_particles": n_particles,
        "charge_range_nc": list(charge_range_nc),
        "axis_range_mm": list(axis_range_mm),
        "softening_fraction": softening_fraction,
        "include_energy": include_energy,
        "quiet_start": True,
        "kv_note": "KV representa sólo su proyección espacial uniforme en un elipsoide.",
        "input_names": list(INPUT_NAMES),
        "output_names": list(output_names),
    }
    return inputs, outputs, metadata


def save_dataset(path: str | Path, inputs: np.ndarray, outputs: np.ndarray, metadata: dict) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        destination,
        X=np.asarray(inputs),
        y=np.asarray(outputs),
        metadata_json=np.asarray(json.dumps(metadata, ensure_ascii=False)),
    )
    return destination


def load_dataset(path: str | Path) -> tuple[np.ndarray, np.ndarray, dict]:
    with np.load(path, allow_pickle=False) as data:
        inputs = np.asarray(data["X"], dtype=float)
        outputs = np.asarray(data["y"], dtype=float)
        metadata = json.loads(str(data["metadata_json"]))
    return inputs, outputs, metadata
