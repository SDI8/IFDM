"""
Dryer and filament configuration, and simulation orchestration.

Ties together the material database, diffusion solver, and experiment
parameters into a single `simulate()` call.
"""

from dataclasses import dataclass

import numpy as np

from .materials import Material
from .model import DiffusionResult, solve_radial_diffusion


def _cross_section_area(diameter: float) -> float:
    """Cross-sectional area of a circular filament [m²]."""
    return np.pi * (diameter / 2.0) ** 2


# ---------------------------------------------------------------------------
# Air properties and psychrometrics
# ---------------------------------------------------------------------------


def saturation_vapor_pressure(T_celsius: float) -> float:
    """Water saturation vapor pressure [Pa] (Alduchov & Eskridge 1996)."""
    return 610.94 * np.exp(17.625 * T_celsius / (243.04 + T_celsius))


def compute_chamber_rh(
    ambient_rh: float,
    ambient_temp: float,
    chamber_temp: float,
) -> float:
    """RH in the dryer after heating ambient intake air.

    Heating air at constant absolute humidity lowers RH because
    saturation vapor pressure rises with temperature.
    """
    return (
        ambient_rh
        * saturation_vapor_pressure(ambient_temp)
        / saturation_vapor_pressure(chamber_temp)
    )


def air_density(T_celsius: float) -> float:
    """Dry-air density at 1 atm [kg/m³]."""
    return 101325.0 / (287.058 * (T_celsius + 273.15))


def air_dynamic_viscosity(T_celsius: float) -> float:
    """Dynamic viscosity of air [Pa·s] (Sutherland's law)."""
    T = T_celsius + 273.15
    return 1.458e-6 * T**1.5 / (T + 110.4)


def water_vapor_diffusivity_in_air(T_celsius: float) -> float:
    """Diffusivity of water vapor in air [m²/s] (Fuller et al.)."""
    T = T_celsius + 273.15
    return 2.16e-5 * (T / 273.15) ** 1.8


def compute_sherwood(Re: float, Sc: float) -> float:
    """Sherwood number for crossflow over a cylinder (Churchill-Bernstein)."""
    Sh = 0.3 + (0.62 * Re**0.5 * Sc ** (1 / 3)) / (1 + (0.4 / Sc) ** (2 / 3)) ** 0.25
    if Re > 282_000:
        Sh *= (1 + (Re / 282_000) ** (5 / 8)) ** (4 / 5)
    return Sh


def compute_mass_transfer_biot(
    D_polymer: float,
    filament_diameter: float,
    airflow_velocity: float,
    chamber_temp: float,
) -> float:
    """Biot number for convective mass transfer at the filament surface.

    Bi_m = h_m · R / D_polymer, where h_m is estimated from the
    Churchill-Bernstein (Sherwood) correlation for crossflow over a
    cylinder.
    """
    d = filament_diameter
    R = d / 2.0

    rho = air_density(chamber_temp)
    mu = air_dynamic_viscosity(chamber_temp)
    D_wa = water_vapor_diffusivity_in_air(chamber_temp)

    Re = rho * airflow_velocity * d / mu
    Sc = mu / (rho * D_wa)
    Sh = compute_sherwood(Re, Sc)

    h_m = Sh * D_wa / d  # convective mass transfer coeff [m/s]
    return h_m * R / D_polymer


@dataclass
class DryerConfig:
    """Physical parameters of the inline drying chamber.

    Attributes:
        chamber_length: Length of the heated zone [m].
        chamber_temp: Air temperature inside the chamber [°C].
        ambient_humidity: Relative humidity of the intake air at ambient_temp
            (0–1, e.g. 0.50 = 50% RH).
        ambient_temp: Temperature of the intake air [°C].
        airflow_velocity: Air velocity over the filament inside the chamber [m/s].
    """

    chamber_length: float = 0.5  # m  (500 mm default)
    chamber_temp: float = 80.0  # °C
    ambient_humidity: float = 0.50  # 50% RH at ambient temp
    ambient_temp: float = 25.0  # °C  (room temperature)
    airflow_velocity: float = 1.0  # m/s

    @property
    def chamber_humidity(self) -> float:
        """RH inside the chamber after heating ambient intake air."""
        return compute_chamber_rh(
            self.ambient_humidity,
            self.ambient_temp,
            self.chamber_temp,
        )


@dataclass
class FilamentConfig:
    """Filament properties and operating conditions.

    Attributes:
        material: Material instance from the database.
        diameter: Filament diameter [m].
        initial_moisture: Initial moisture mass fraction (e.g. 0.03 = 3 wt%).
        flow_rate: Volumetric extrusion rate [mm³/s]. This is how slicer
            software specifies speed. Linear filament speed is derived from
            flow_rate and diameter.
    """

    material: Material
    diameter: float = 1.75e-3  # m  (1.75 mm)
    initial_moisture: float = 0.03  # 3 wt% — a moderately moist filament
    flow_rate: float = 8.0  # mm³/s  (typical extrusion rate)

    @property
    def speed(self) -> float:
        """Linear filament feed speed [m/s], derived from flow_rate and diameter."""
        area_mm2 = _cross_section_area(self.diameter) * 1e6  # m² → mm²
        return (self.flow_rate / area_mm2) * 1e-3  # mm/s → m/s


@dataclass
class DryingResult:
    """Output of a drying simulation.

    Attributes:
        diffusion: Full DiffusionResult (radial profiles over time).
        transit_time: Time the filament spends in the dryer [s].
        initial_moisture: Starting moisture mass fraction.
        final_moisture: Volume-averaged moisture at exit.
        moisture_removed: Absolute moisture removed (mass fraction).
        drying_efficiency: Fraction of initial moisture removed (0–1).
        fourier_number: Dimensionless Fo = D·t/R², indicating drying depth.
        filament: FilamentConfig used.
        dryer: DryerConfig used.
    """

    diffusion: DiffusionResult
    transit_time: float
    initial_moisture: float
    final_moisture: float
    moisture_removed: float
    drying_efficiency: float
    fourier_number: float
    biot_mass: float
    filament: FilamentConfig
    dryer: DryerConfig

    def summary(self) -> str:
        """Human-readable summary of drying result."""
        lines = [
            f"Material:           {self.filament.material.name}",
            f"Filament diameter:  {self.filament.diameter * 1e3:.2f} mm",
            f"Flow rate:          {self.filament.flow_rate:.1f} mm³/s",
            f"Filament speed:     {self.filament.speed * 1e3:.2f} mm/s (linear)",
            f"Chamber length:     {self.dryer.chamber_length * 1e3:.0f} mm",
            f"Chamber temp:       {self.dryer.chamber_temp:.0f} °C",
            f"Ambient conditions: {self.dryer.ambient_temp:.0f} °C, {
                self.dryer.ambient_humidity * 100:.0f\
            }% RH",
            f"Chamber humidity:   {self.dryer.chamber_humidity * 100:.2f}% RH",
            f"Airflow velocity:   {self.dryer.airflow_velocity:.1f} m/s",
            f"Transit time:       {self.transit_time:.2f} s",
            f"Fourier number:     {self.fourier_number:.4e}",
            f"Biot (mass xfer):   {self.biot_mass:.2e}",
            f"D(T):               {
                self.filament.material.diffusivity(self.dryer.chamber_temp):.3e\
            } m²/s",
            f"Initial moisture:   {self.initial_moisture * 100:.3f} wt%",
            f"Final moisture:     {self.final_moisture * 100:.3f} wt%",
            f"Moisture removed:   {self.moisture_removed * 100:.4f} wt%",
            f"Drying efficiency:  {self.drying_efficiency * 100:.2f} %",
        ]
        return "\n".join(lines)


def simulate(
    dryer: DryerConfig,
    filament: FilamentConfig,
    N: int = 50,
    t_eval_count: int = 200,
) -> DryingResult:
    """Run an inline-drying simulation.

    Computes transit time, evaluates diffusivity at chamber temperature,
    solves the radial diffusion PDE, and returns a DryingResult.
    """
    transit_time = dryer.chamber_length / filament.speed
    R = filament.diameter / 2.0
    D = filament.material.diffusivity(dryer.chamber_temp)

    # Equilibrium moisture at the surface for the chamber air conditions.
    # Linear sorption isotherm: C_eq = M_sat × RH_chamber
    C_env = filament.material.equilibrium_moisture * dryer.chamber_humidity

    # Mass-transfer Biot number from Sherwood correlation
    biot = compute_mass_transfer_biot(
        D,
        filament.diameter,
        dryer.airflow_velocity,
        dryer.chamber_temp,
    )

    diff = solve_radial_diffusion(
        R=R,
        D=D,
        C_init=filament.initial_moisture,
        t_end=transit_time,
        C_env=C_env,
        N=N,
        t_eval_count=t_eval_count,
        biot_mass=biot,
    )

    final_moisture = diff.C_avg[-1]
    moisture_removed = filament.initial_moisture - final_moisture
    efficiency = (
        moisture_removed / filament.initial_moisture if filament.initial_moisture > 0 else 0.0
    )
    Fo = D * transit_time / R**2

    return DryingResult(
        diffusion=diff,
        transit_time=transit_time,
        initial_moisture=filament.initial_moisture,
        final_moisture=final_moisture,
        moisture_removed=moisture_removed,
        drying_efficiency=efficiency,
        fourier_number=Fo,
        biot_mass=biot,
        filament=filament,
        dryer=dryer,
    )
