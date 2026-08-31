"""Generación reproducible de nubes de macropartículas en tres dimensiones."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc


GEOMETRIES = ("sphere", "gaussian", "kv")

_OCTANT_SIGNS = np.asarray(
    [
        (-1.0, -1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, -1.0),
        (1.0, 1.0, 1.0),
    ]
)


@dataclass(frozen=True)
class GeometryParameters:
    """Parámetros físicos de una distribución alineada con los ejes."""

    geometry: str
    charge_nc: float
    a_mm: float
    b_mm: float
    c_mm: float

    def __post_init__(self) -> None:
        if self.geometry not in GEOMETRIES:
            raise ValueError(f"Geometría desconocida: {self.geometry!r}")
        if self.charge_nc <= 0.0:
            raise ValueError("charge_nc debe ser positivo (magnitud de la carga).")
        if min(self.a_mm, self.b_mm, self.c_mm) <= 0.0:
            raise ValueError("Todas las longitudes deben ser positivas.")
        if self.geometry == "sphere" and not np.allclose(
            (self.a_mm, self.b_mm, self.c_mm), self.a_mm
        ):
            raise ValueError("Para sphere se requiere a_mm = b_mm = c_mm = radio.")

    @property
    def axes_m(self) -> np.ndarray:
        """Radio, sigmas o semiejes en metros, según la geometría."""

        return 1.0e-3 * np.asarray((self.a_mm, self.b_mm, self.c_mm))

    @property
    def total_charge_c(self) -> float:
        return self.charge_nc * 1.0e-9


def _uniform_unit_sphere(n: int, rng: np.random.Generator) -> np.ndarray:
    directions = rng.normal(size=(n, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    radii = rng.random(n) ** (1.0 / 3.0)
    return directions * radii[:, None]


def _quasi_random_cube(n: int, rng: np.random.Generator) -> np.ndarray:
    """Secuencia Halton aleatorizada, reproducible a partir del generador."""

    seed = int(rng.integers(0, np.iinfo(np.uint32).max))
    return qmc.Halton(d=3, scramble=True, seed=seed).random(n)


def _stratified_first_octant_sphere(n: int, rng: np.random.Generator) -> np.ndarray:
    """Esfera uniforme de baja discrepancia en r, cos(theta) y phi."""

    cube = _quasi_random_cube(n, rng)
    radius = cube[:, 0] ** (1.0 / 3.0)
    cos_theta = cube[:, 1]
    phi = 0.5 * np.pi * cube[:, 2]
    sin_theta = np.sqrt(1.0 - cos_theta**2)
    return np.column_stack(
        (
            radius * sin_theta * np.cos(phi),
            radius * sin_theta * np.sin(phi),
            radius * cos_theta,
        )
    )


def _mirror_to_octants(base: np.ndarray) -> np.ndarray:
    """Refleja puntos del primer octante para cancelar ruido dipolar."""

    return (np.abs(base)[:, None, :] * _OCTANT_SIGNS[None, :, :]).reshape(-1, 3)


def sample_unit_cloud(
    geometry: str,
    n_particles: int,
    rng: np.random.Generator,
    quiet_start: bool = True,
) -> np.ndarray:
    """Genera una nube adimensional que después puede escalarse.

    ``sphere`` y ``kv`` producen densidad espacial uniforme dentro de una
    esfera unidad. ``gaussian`` produce una normal estándar en cada eje.

    La opción ``kv`` representa sólo la proyección espacial uniforme de una
    distribución tipo KV. Un KV completo requiere también el espacio de fases
    transversal 4D.
    """

    if geometry not in GEOMETRIES:
        raise ValueError(f"Geometría desconocida: {geometry!r}")
    if n_particles <= 0:
        raise ValueError("n_particles debe ser positivo.")

    if quiet_start:
        if n_particles % 8:
            raise ValueError("Con quiet_start, n_particles debe ser múltiplo de 8.")
        n_base = n_particles // 8
        if geometry == "gaussian":
            # El valor absoluto de N(0,1) tiene CDF 2*Phi(x)-1. Muestrear
            # Phi(x) en [0.5, 1) genera directamente el primer octante.
            cube = _quasi_random_cube(n_base, rng)
            base = ndtri(0.5 + 0.5 * cube)
        else:
            base = _stratified_first_octant_sphere(n_base, rng)
        return _mirror_to_octants(base)

    if geometry == "gaussian":
        return rng.normal(size=(n_particles, 3))
    return _uniform_unit_sphere(n_particles, rng)


def scale_unit_cloud(unit_cloud: np.ndarray, params: GeometryParameters) -> np.ndarray:
    """Convierte una nube unidad a posiciones físicas en metros."""

    cloud = np.asarray(unit_cloud, dtype=float)
    if cloud.ndim != 2 or cloud.shape[1] != 3:
        raise ValueError("unit_cloud debe tener forma (N, 3).")
    return cloud * params.axes_m[None, :]


def expected_rms(params: GeometryParameters) -> np.ndarray:
    """Valores RMS continuos esperados para las posiciones."""

    if params.geometry == "gaussian":
        return params.axes_m
    return params.axes_m / np.sqrt(5.0)
