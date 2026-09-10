"""RFQ design formulas used by the iterative TRANSOPTR designer.

The functions in this module are intentionally small and unit-testable.  They
cover the kinematic and two-term potential quantities used throughout the
current workflow:

* beta and RF wavelength
* RFQ half-cell length, L = beta * lambda / 2
* spatial wave number, k = pi / L
* Kapchinsky-Teplyakov two-term coefficients A01 and A10

The user's design reference is rfq_book.pdf, pages 23-32.  The exact notation
in those pages should be mirrored here as this tool is refined; v1 keeps the
same conventions already used in the existing import_re.py/create_warp_field.py
scripts so the generated fort.75 is compatible with the rest of the repo.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


C_LIGHT_M_PER_S = 299_792_458.0
PROTON_MASS_MEV = 938.272089


@dataclass(frozen=True)
class RfqKinematics:
    beta: float
    gamma: float
    wavelength_cm: float
    cell_length_cm: float
    k_cm_inv: float


def beta_gamma_from_kinetic_energy(kinetic_energy_mev: float, mass_mev: float = PROTON_MASS_MEV) -> tuple[float, float]:
    """Return relativistic beta and gamma from kinetic energy in MeV."""
    if mass_mev <= 0:
        raise ValueError("mass_mev must be positive")
    if kinetic_energy_mev < 0:
        raise ValueError("kinetic_energy_mev must be non-negative")

    gamma = 1.0 + kinetic_energy_mev / mass_mev
    beta = math.sqrt(max(0.0, 1.0 - 1.0 / (gamma * gamma)))
    return beta, gamma


def rf_wavelength_cm(frequency_hz: float) -> float:
    """Return RF wavelength in cm."""
    if frequency_hz <= 0:
        raise ValueError("frequency_hz must be positive")
    return C_LIGHT_M_PER_S / frequency_hz * 100.0


def rfq_cell_kinematics(kinetic_energy_mev: float, frequency_hz: float, mass_mev: float = PROTON_MASS_MEV) -> RfqKinematics:
    """Return beta, gamma, lambda, half-cell length and k for one RFQ cell."""
    beta, gamma = beta_gamma_from_kinetic_energy(kinetic_energy_mev, mass_mev)
    wavelength_cm = rf_wavelength_cm(frequency_hz)
    cell_length_cm = beta * wavelength_cm / 2.0

    if cell_length_cm <= 0:
        raise ValueError("cell length is zero; check input energy and frequency")

    return RfqKinematics(
        beta=beta,
        gamma=gamma,
        wavelength_cm=wavelength_cm,
        cell_length_cm=cell_length_cm,
        k_cm_inv=math.pi / cell_length_cm,
    )


def kt_coefficients_from_a_m_k(a_cm: float, m: float, k_cm_inv: float) -> tuple[float, float]:
    """Return two-term KT coefficients A01 and A10.

    This is the convention already used in import_re.py:

        A10 = (m^2 - 1) / (m^2 I0(k a) + I0(k m a))
        A01 = (1 - A10 I0(k a)) / a^2

    a_cm and k_cm_inv must use consistent cm units.
    """
    if a_cm <= 0:
        raise ValueError("a_cm must be positive")
    if m < 1.0:
        raise ValueError("m must be >= 1")
    if k_cm_inv <= 0:
        raise ValueError("k_cm_inv must be positive")

    ka = k_cm_inv * a_cm
    kma = k_cm_inv * m * a_cm
    i0_ka = np.i0(ka)
    i0_kma = np.i0(kma)
    denom = m * m * i0_ka + i0_kma

    A10 = (m * m - 1.0) / denom
    A01 = (1.0 - A10 * i0_ka) / (a_cm * a_cm)
    return float(A01), float(A10)


def linear_schedule(cell_index: int, max_cells: int, start: float, end: float) -> float:
    """Return a linearly scheduled value from start to end over max_cells."""
    if max_cells <= 1:
        return float(end)
    t = min(max(cell_index - 1, 0), max_cells - 1) / float(max_cells - 1)
    return float(start + t * (end - start))


def predicted_energy_gain_mev(
    vane_voltage_kv: float,
    charge_state: float,
    A10: float,
    phi_deg: float,
    gain_scale: float = 1.0,
) -> float:
    """Average two-term gain for constant coefficients over one half-cell.

    phi_deg = omega*t + phi_RF - integral(k ds), as in TRANSOPTR.
    Integral sin(psi)*sin(psi+phi) dpsi from 0 to pi = pi/2*cos(phi).
    This approximation is diagnostic only; real tracking sets the energy.
    """
    return (math.pi/4 * charge_state * vane_voltage_kv * 1.0e-3
            * A10 * math.cos(math.radians(phi_deg)) * gain_scale)


def transverse_phase_advance_deg(B: float, a_cm: float, m: float, voltage_kv: float) -> float:
    """Simple report-only transverse metric for v1.

    This is not used as an optimizer yet.  It gives a monotonic diagnostic that
    rises with focusing strength and voltage, and falls with aperture.
    """
    if a_cm <= 0:
        return 0.0
    return float(max(0.0, 20.0 * math.sqrt(max(B, 0.0) * max(voltage_kv, 0.0) / a_cm) / max(m, 1.0)))


def longitudinal_phase_advance_deg(A10: float, voltage_kv: float, phi_deg: float) -> float:
    """Simple report-only longitudinal metric for v1."""
    focusing = abs(A10) * max(voltage_kv, 0.0) * max(0.0, -math.sin(math.radians(phi_deg)))
    return float(10.0 * math.sqrt(max(focusing, 0.0)))

