"""
Dryer and filament configuration, and simulation orchestration.

Ties together the material database, diffusion solver, and experiment
parameters into a single `simulate()` call.
"""

from dataclasses import dataclass

import numpy as np

from .materials import Material


def _cross_section_area(diameter: float) -> float:
    """Cross-sectional area of a circular filament [m²]."""
    return np.pi * (diameter / 2.0) ** 2
from .model import DiffusionResult, solve_radial_diffusion, volume_average


@dataclass
class DryerConfig:
    """Physical parameters of the inline drying tube.

    Attributes:
        tube_length: Length of the heated zone [m].
        chamber_temp: Air temperature inside the tube [°C].
        chamber_humidity: Relative humidity of drying air (0–1).
            Currently unused (Dirichlet BC assumes dry air). Reserved for
            future convective BC implementation.
        airflow_velocity: Air velocity in tube [m/s].
            Currently unused. Reserved for future Sherwood-number-based
            convective mass-transfer BC.
    """

    tube_length: float = 0.5          # m  (500 mm default)
    chamber_temp: float = 80.0        # °C
    chamber_humidity: float = 0.0     # 0 = perfectly dry air
    airflow_velocity: float = 1.0     # m/s  (reserved for future use)


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
    diameter: float = 1.75e-3         # m  (1.75 mm)
    initial_moisture: float = 0.03    # 3 wt% — a moderately moist filament
    flow_rate: float = 8.0            # mm³/s  (typical extrusion rate)

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
    filament: FilamentConfig
    dryer: DryerConfig

    def summary(self) -> str:
        """Human-readable summary of drying result."""
        lines = [
            f"Material:           {self.filament.material.name}",
            f"Filament diameter:  {self.filament.diameter * 1e3:.2f} mm",
            f"Flow rate:          {self.filament.flow_rate:.1f} mm³/s",
            f"Filament speed:     {self.filament.speed * 1e3:.2f} mm/s (linear)",
            f"Tube length:        {self.dryer.tube_length * 1e3:.0f} mm",
            f"Chamber temp:       {self.dryer.chamber_temp:.0f} °C",
            f"Transit time:       {self.transit_time:.2f} s",
            f"Fourier number:     {self.fourier_number:.4e}",
            f"D(T):               {self.filament.material.diffusivity(self.dryer.chamber_temp):.3e} m²/s",
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
    transit_time = dryer.tube_length / filament.speed
    R = filament.diameter / 2.0
    D = filament.material.diffusivity(dryer.chamber_temp)

    # Environmental moisture concentration at filament surface.
    # For the Dirichlet model: C_env ≈ equilibrium moisture * RH
    C_env = filament.material.equilibrium_moisture * dryer.chamber_humidity

    diff = solve_radial_diffusion(
        R=R,
        D=D,
        C_init=filament.initial_moisture,
        t_end=transit_time,
        C_env=C_env,
        N=N,
        t_eval_count=t_eval_count,
    )

    final_moisture = diff.C_avg[-1]
    moisture_removed = filament.initial_moisture - final_moisture
    efficiency = moisture_removed / filament.initial_moisture if filament.initial_moisture > 0 else 0.0
    Fo = D * transit_time / R**2

    return DryingResult(
        diffusion=diff,
        transit_time=transit_time,
        initial_moisture=filament.initial_moisture,
        final_moisture=final_moisture,
        moisture_removed=moisture_removed,
        drying_efficiency=efficiency,
        fourier_number=Fo,
        filament=filament,
        dryer=dryer,
    )
