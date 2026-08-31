"""Suma de Coulomb en SI y soluciones analíticas usadas para validación."""

from __future__ import annotations

from math import erf

import numpy as np


EPSILON_0 = 8.854_187_812_8e-12
COULOMB_CONSTANT = 1.0 / (4.0 * np.pi * EPSILON_0)


def electric_field(
    probes_m: np.ndarray,
    positions_m: np.ndarray,
    macro_charge_c: float,
    softening_m: float = 0.0,
    block_size: int = 4096,
) -> np.ndarray:
    """Campo eléctrico en puntos de observación, sumando macropartículas.

    Se usa un núcleo de Plummer: ``(r² + epsilon²)^(-3/2)``. Esto evita
    singularidades artificiales cuando una macropartícula queda muy cerca de
    una sonda.
    """

    probes = np.atleast_2d(np.asarray(probes_m, dtype=float))
    positions = np.asarray(positions_m, dtype=float)
    if probes.shape[1:] != (3,) or positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("probes_m y positions_m deben tener forma (*, 3).")
    if block_size <= 0 or softening_m < 0.0:
        raise ValueError("block_size debe ser positivo y softening_m no negativo.")

    result = np.zeros((len(probes), 3), dtype=float)
    eps2 = float(softening_m) ** 2
    for start in range(0, len(positions), block_size):
        block = positions[start : start + block_size]
        delta = probes[:, None, :] - block[None, :, :]
        r2 = np.einsum("mni,mni->mn", delta, delta) + eps2
        with np.errstate(divide="ignore", invalid="ignore"):
            inv_r3 = np.where(r2 > 0.0, r2 ** -1.5, 0.0)
        result += np.einsum("mn,mni->mi", inv_r3, delta)
    return COULOMB_CONSTANT * macro_charge_c * result


def electric_potential(
    probes_m: np.ndarray,
    positions_m: np.ndarray,
    macro_charge_c: float,
    softening_m: float = 0.0,
    block_size: int = 4096,
) -> np.ndarray:
    """Potencial electrostático en puntos de observación."""

    probes = np.atleast_2d(np.asarray(probes_m, dtype=float))
    positions = np.asarray(positions_m, dtype=float)
    result = np.zeros(len(probes), dtype=float)
    eps2 = float(softening_m) ** 2
    for start in range(0, len(positions), block_size):
        block = positions[start : start + block_size]
        delta = probes[:, None, :] - block[None, :, :]
        r2 = np.einsum("mni,mni->mn", delta, delta) + eps2
        with np.errstate(divide="ignore"):
            inv_r = np.where(r2 > 0.0, r2 ** -0.5, 0.0)
        result += inv_r.sum(axis=1)
    return COULOMB_CONSTANT * macro_charge_c * result


def electrostatic_energy(
    positions_m: np.ndarray,
    macro_charge_c: float,
    softening_m: float = 0.0,
    block_size: int = 1024,
) -> float:
    """Energía de Coulomb de todos los pares, con costo O(N²).

    Esta función nunca materializa una matriz completa N×N. Aun así, debe
    activarse con cuidado al generar muchos ejemplos de entrenamiento.
    """

    positions = np.asarray(positions_m, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions_m debe tener forma (N, 3).")
    eps2 = float(softening_m) ** 2
    inverse_distance_sum = 0.0
    n_particles = len(positions)

    for i_start in range(0, n_particles, block_size):
        i_stop = min(i_start + block_size, n_particles)
        left = positions[i_start:i_stop]
        for j_start in range(i_start, n_particles, block_size):
            j_stop = min(j_start + block_size, n_particles)
            right = positions[j_start:j_stop]
            delta = left[:, None, :] - right[None, :, :]
            r2 = np.einsum("ijk,ijk->ij", delta, delta) + eps2
            if i_start == j_start:
                row, col = np.triu_indices(len(left), k=1)
                inverse_distance_sum += np.sum(r2[row, col] ** -0.5)
            else:
                inverse_distance_sum += np.sum(r2 ** -0.5)

    return float(COULOMB_CONSTANT * macro_charge_c**2 * inverse_distance_sum)


def uniform_sphere_field(radius_m: np.ndarray, total_charge_c: float, sphere_radius_m: float) -> np.ndarray:
    """Magnitud analítica de E(r) para una esfera de carga uniforme."""

    radius = np.asarray(radius_m, dtype=float)
    inside = COULOMB_CONSTANT * total_charge_c * radius / sphere_radius_m**3
    outside = np.divide(
        COULOMB_CONSTANT * total_charge_c,
        radius**2,
        out=np.zeros_like(radius),
        where=radius > 0.0,
    )
    return np.where(radius < sphere_radius_m, inside, outside)


def gaussian_sphere_field(radius_m: np.ndarray, total_charge_c: float, sigma_m: float) -> np.ndarray:
    """Magnitud analítica de E(r) para una gaussiana esférica 3D."""

    radius = np.asarray(radius_m, dtype=float)
    x = radius / sigma_m
    enclosed_fraction = np.asarray(
        [erf(value / np.sqrt(2.0)) for value in x.ravel()]
    ).reshape(x.shape)
    enclosed_fraction -= np.sqrt(2.0 / np.pi) * x * np.exp(-0.5 * x**2)
    return np.divide(
        COULOMB_CONSTANT * total_charge_c * enclosed_fraction,
        radius**2,
        out=np.zeros_like(radius),
        where=radius > 0.0,
    )


def uniform_sphere_center_potential(total_charge_c: float, sphere_radius_m: float) -> float:
    return 1.5 * COULOMB_CONSTANT * total_charge_c / sphere_radius_m


def gaussian_center_potential(total_charge_c: float, sigma_m: float) -> float:
    return COULOMB_CONSTANT * total_charge_c * np.sqrt(2.0 / np.pi) / sigma_m

