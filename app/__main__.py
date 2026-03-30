"""
Inline Filament Dryer — Numerical Experiments

Run with:  python -m app
"""

from .dryer import DryerConfig, FilamentConfig, simulate
from .materials import MATERIALS
from .model import analytical_moisture_fraction
from .optimize import OptimizationConfig, optimize


def run_default_experiment(filament: FilamentConfig):
    """Run a baseline drying simulation for the given filament and print results."""
    dryer = DryerConfig(chamber_length=0.5, chamber_temp=80.0)

    result = simulate(dryer, filament)
    print("=" * 60)
    print("  INLINE FILAMENT DRYER — Baseline Simulation")
    print("=" * 60)
    print(result.summary())
    print("=" * 60)

    # --- Validate against analytical solution ---
    R = filament.diameter / 2.0
    D = filament.material.diffusivity(dryer.chamber_temp)
    C_env = filament.material.equilibrium_moisture * dryer.chamber_humidity
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


def run_optimization_experiment(filament: FilamentConfig):
    """Find optimal dryer parameters for the given filament and compare to baseline."""

    config = OptimizationConfig(
        filament=filament,
        bounds_length=(0.1, 5.0),
        bounds_temp=(40.0, 120.0),
        bounds_airflow=(0.01, 10.0),
        weight_moisture=1.0,
        weight_length=0.1,
        weight_energy=0.1,
    )

    print("\n" + "=" * 60)
    print("  OPTIMIZATION — Finding optimal dryer parameters")
    print("  Weights: moisture=1.0, length=0.1, energy=0.1")
    print("=" * 60)

    result = optimize(config)

    print(result.summary())
    print("=" * 60)

    # Compare to baseline
    baseline = simulate(
        DryerConfig(chamber_length=0.5, chamber_temp=80.0),
        filament,
    )
    print("\nBaseline vs Optimized:")
    print(f"  {'':25s} {'Baseline':>10s}  {'Optimized':>10s}")
    print(f"  {'Chamber length [mm]':25s} {500:10.0f}  {result.optimal_length * 1e3:10.0f}")
    print(f"  {'Chamber temp [°C]':25s} {80:10.1f}  {result.optimal_temp:10.1f}")
    print(f"  {'Airflow [m/s]':25s} {1.0:10.1f}  {result.optimal_airflow:10.1f}")
    print(
        f"  {'Final moisture [wt%]':25s}"
        f" {baseline.final_moisture * 100:10.3f}"
        f"  {result.simulation.final_moisture * 100:10.3f}"
    )
    print(
        f"  {'Drying efficiency [%]':25s}"
        f" {baseline.drying_efficiency * 100:10.2f}"
        f"  {result.simulation.drying_efficiency * 100:10.2f}"
    )

    return result


def main():
    filament = FilamentConfig(material=MATERIALS["PETG"], initial_moisture=0.005, flow_rate=8.0)

    run_default_experiment(filament)
    run_optimization_experiment(filament)


if __name__ == "__main__":
    main()
