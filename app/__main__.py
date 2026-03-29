"""
Inline Filament Dryer — Numerical Experiments

Run with:  python -m app
"""

import numpy as np

from .dryer import DryerConfig, FilamentConfig, simulate
from .experiments import (
    compare_materials,
    sweep_chamber_length,
    sweep_length_flow_rate,
    sweep_temperature,
)
from .materials import MATERIALS, PA6
from .model import analytical_moisture_fraction


def run_default_experiment():
    """Run a baseline PA6 drying simulation and print results."""
    dryer = DryerConfig(chamber_length=0.5, chamber_temp=80.0)
    filament = FilamentConfig(material=PA6, initial_moisture=0.03, flow_rate=8.0)

    result = simulate(dryer, filament)
    print("=" * 60)
    print("  INLINE FILAMENT DRYER — Baseline Simulation")
    print("=" * 60)
    print(result.summary())
    print("=" * 60)

    # --- Validate against analytical solution ---
    R = filament.diameter / 2.0
    D = PA6.diffusivity(dryer.chamber_temp)
    C_env = PA6.equilibrium_moisture * dryer.chamber_humidity
    frac_numerical = (result.final_moisture - C_env) / (result.initial_moisture - C_env)
    frac_analytical = float(analytical_moisture_fraction(D, R, result.transit_time))
    print("\nValidation vs analytical (Dirichlet BC, constant D):")
    print(f"  (Bi_m = {result.biot_mass:.0e} → Robin BC ≈ Dirichlet)")
    print(f"  C_env = {C_env * 100:.3f} wt% (equilibrium at chamber RH)")
    print(f"  Numerical  (M-M_eq)/(M₀-M_eq) = {frac_numerical:.6f}")
    print(f"  Analytical (M-M_eq)/(M₀-M_eq) = {frac_analytical:.6f}")
    print(
        f"  Relative error                 = {abs(frac_numerical - frac_analytical) / frac_analytical * 100:.3f} %"
    )

    return result


def main():
    run_default_experiment()


if __name__ == "__main__":
    main()
