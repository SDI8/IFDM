"""
Inline Filament Dryer — Numerical Experiments

Run with:  python -m app
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
from .plotting import (
    plot_heatmap,
    plot_material_comparison,
    plot_moisture_vs_time,
    plot_radial_profile,
    plot_sweep,
)

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output"


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
    # The Crank series gives (M(t) - M_eq) / (M₀ - M_eq) for a cylinder
    # with constant surface concentration C_env.  Our Robin BC with
    # Bi >> 1 should match closely.
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


def run_visualizations(result):
    """Generate diagnostic plots."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    plot_radial_profile(result, ax=axes[0])
    plot_moisture_vs_time(result, ax=axes[1])
    fig.tight_layout()

    # --- Temperature sweep ---
    temps = np.linspace(40, 120, 30)
    dryer = DryerConfig(chamber_length=0.5, chamber_temp=80.0)
    filament = result.filament
    moist_t, _ = sweep_temperature(temps, dryer, filament)
    fig2, ax2 = plt.subplots(figsize=(7, 5))
    plot_sweep(temps, moist_t, "Chamber temperature", "°C", filament.material.name, ax=ax2)

    # --- Chamber length sweep ---
    lengths = np.linspace(0.1, 2.0, 30)
    moist_l, _ = sweep_chamber_length(lengths, dryer, filament)
    fig3, ax3 = plt.subplots(figsize=(7, 5))
    plot_sweep(lengths * 1e3, moist_l, "Chamber length", "mm", filament.material.name, ax=ax3)

    # --- 2D heatmap: length × flow rate ---
    lengths_2d = np.linspace(0.1, 2.0, 25)
    flow_rates_2d = np.linspace(1.0, 20.0, 25)  # mm³/s
    grid = sweep_length_flow_rate(lengths_2d, flow_rates_2d, dryer, filament)
    fig4, ax4 = plt.subplots(figsize=(8, 6))
    plot_heatmap(
        lengths_2d * 1e3,
        flow_rates_2d,
        grid,
        "Chamber length [mm]",
        "Flow rate [mm³/s]",
        f"{filament.material.name} — Final moisture",
        ax=ax4,
    )

    # --- Material comparison ---
    materials = [MATERIALS[k] for k in ("PA6", "PETG", "PLA")]
    dryer_cmp = DryerConfig(chamber_length=1.0, chamber_temp=80.0)
    cmp_results = compare_materials(materials, dryer_cmp)
    fig5, ax5 = plt.subplots(figsize=(8, 5))
    plot_material_comparison(cmp_results, ax=ax5)

    # Save all figures
    OUTPUT_DIR.mkdir(exist_ok=True)
    for i, f in enumerate(plt.get_fignums()):
        names = ["profiles", "temp_sweep", "length_sweep", "heatmap", "material_cmp"]
        fname = names[i] if i < len(names) else f"figure_{i}"
        plt.figure(f).savefig(OUTPUT_DIR / f"{fname}.png", dpi=150, bbox_inches="tight")
    print(f"\nPlots saved to {OUTPUT_DIR}/")


def main():
    result = run_default_experiment()
    run_visualizations(result)


if __name__ == "__main__":
    main()
