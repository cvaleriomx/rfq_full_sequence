from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np


PROTON_REST_ENERGY_EV = 938.27208816e6
LIGHT_SPEED_MPS = 299_792_458.0


def proton_speed(energy_ev: float) -> float:
    gamma = 1.0 + energy_ev / PROTON_REST_ENERGY_EV
    beta = np.sqrt(1.0 - gamma**-2)
    return float(beta * LIGHT_SPEED_MPS)


def tracking_plan(field_path: str | Path, tracking: dict) -> dict:
    with h5py.File(field_path, "r") as h5:
        x = h5["axes/x_m"][:]
        y = h5["axes/y_m"][:]
        z = h5["axes/z_m"][:]
        field_shape = h5["fields/ex_vpm_per_v"].shape
    speed = proton_speed(float(tracking["kinetic_energy_ev"]))
    dt = float(tracking.get("step_fraction_of_field_dz", 1.0)) * float(z[1] - z[0]) / speed
    travel = float(z[-1] - z[0])
    steps = int(np.ceil(travel / (speed * dt))) + 1
    return {
        "field_shape": list(field_shape),
        "field_bounds_m": {"x": [float(x[0]), float(x[-1])], "y": [float(y[0]), float(y[-1])], "z": [float(z[0]), float(z[-1])]},
        "particles": int(tracking["particles"]),
        "kinetic_energy_ev": float(tracking["kinetic_energy_ev"]),
        "speed_m_per_s": speed,
        "dt_s": dt,
        "steps": steps,
        "rf_frequency_hz": float(tracking["rf_frequency_hz"]),
        "intervane_voltage_v": float(tracking["intervane_voltage_v"]),
        "rf_cycles_simulated": float(steps * dt * float(tracking["rf_frequency_hz"])),
        "estimated_field_memory_mib": float(3 * np.prod(field_shape) * 4 / 1024**2),
    }


def run_warp_tracking(field_path: str | Path, output_path: str | Path, tracking: dict) -> dict:
    """Track one deterministic proton bunch through the external RFQ field."""
    try:
        import warp as wp
    except ImportError as exc:
        raise RuntimeError(
            "Warp no está disponible. Activa el ambiente warp_2 o usa track --dry-run."
        ) from exc

    plan = tracking_plan(field_path, tracking)
    rng = np.random.default_rng(int(tracking.get("random_seed", 2026)))
    with h5py.File(field_path, "r") as h5:
        x_axis = h5["axes/x_m"][:]
        y_axis = h5["axes/y_m"][:]
        z_axis = h5["axes/z_m"][:]
        ex = h5["fields/ex_vpm_per_v"][:]
        ey = h5["fields/ey_vpm_per_v"][:]
        ez = h5["fields/ez_vpm_per_v"][:]

    wp.top.lprntpara = False
    wp.top.lpsplots = False
    wp.top.nesmult = 1
    protons = wp.Species(type=wp.Proton, charge_state=1, name="protons")
    protons.ekin = float(tracking["kinetic_energy_ev"])
    protons.ibeam = float(tracking.get("beam_current_a", 0.0))
    wp.top.npmax = int(tracking["particles"])
    wp.derivqty()

    wp.w3d.nx = int(tracking.get("solver_nx", 32))
    wp.w3d.ny = int(tracking.get("solver_ny", 32))
    wp.w3d.nz = int(tracking.get("solver_nz", 64))
    wp.w3d.xmmin, wp.w3d.xmmax = float(x_axis[0]), float(x_axis[-1])
    wp.w3d.ymmin, wp.w3d.ymmax = float(y_axis[0]), float(y_axis[-1])
    wp.w3d.zmmin, wp.w3d.zmmax = float(z_axis[0]), float(z_axis[-1])
    wp.w3d.bound0 = wp.neumann
    wp.w3d.boundnz = wp.neumann
    wp.w3d.boundxy = wp.dirichlet
    wp.top.pbound0 = wp.absorb
    wp.top.pboundnz = wp.absorb
    wp.top.pboundxy = wp.absorb
    wp.top.fstype = 7 if bool(tracking.get("space_charge", False)) else -1
    wp.top.dt = plan["dt_s"]
    wp.package("w3d")
    wp.generate()

    times = np.arange(plan["steps"] + 2, dtype=float) * wp.top.dt
    signal = float(tracking["intervane_voltage_v"]) * np.cos(
        2.0 * np.pi * float(tracking["rf_frequency_hz"]) * times
        + np.deg2rad(float(tracking.get("rf_phase_deg", 0.0)))
    )
    wp.addnewegrd(
        zs=float(z_axis[0]),
        ze=float(z_axis[-1]),
        dx=float(x_axis[1] - x_axis[0]),
        dy=float(y_axis[1] - y_axis[0]),
        xs=float(x_axis[0]),
        ys=float(y_axis[0]),
        time=times,
        data=signal,
        ex=ex,
        ey=ey,
        ez=ez,
    )

    count = int(tracking["particles"])
    radius = float(tracking["beam_radius_m"])
    radial = radius * np.sqrt(rng.random(count))
    angle = 2.0 * np.pi * rng.random(count)
    x = radial * np.cos(angle)
    y = radial * np.sin(angle)
    z = np.full(count, float(z_axis[0]) + 0.25 * (z_axis[1] - z_axis[0]))
    transverse_sigma = float(tracking.get("transverse_velocity_sigma_mps", 0.0))
    vx = rng.normal(0.0, transverse_sigma, count)
    vy = rng.normal(0.0, transverse_sigma, count)
    total_speed = plan["speed_m_per_s"]
    vz = np.sqrt(np.maximum(total_speed**2 - vx**2 - vy**2, 0.0))
    protons.addparticles(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, lallindomain=True)
    wp.step(int(plan["steps"]))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    final = {
        "x_m": protons.getx(),
        "y_m": protons.gety(),
        "z_m": protons.getz(),
        "vx_mps": protons.getvx(),
        "vy_mps": protons.getvy(),
        "vz_mps": protons.getvz(),
    }
    np.savez_compressed(output_path, **final)
    summary = dict(plan)
    summary.update({"surviving_particles": int(len(final["x_m"])), "output": str(output_path)})
    output_path.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary

