"""
Matplotlib visualizations for drying simulation results.
"""

from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from .dryer import DryingResult


def plot_radial_profile(result: DryingResult, ax: plt.Axes | None = None) -> plt.Axes:
    """Plot the radial moisture profile at end of transit."""
    if ax is None:
        _, ax = plt.subplots()

    r_mm = result.diffusion.r * 1e3  # convert to mm
    C_final = result.diffusion.C[:, -1]
    C_init = result.diffusion.C[:, 0]

    ax.plot(r_mm, C_init * 100, "--", color="0.5", label="Initial")
    ax.plot(r_mm, C_final * 100, "b-", linewidth=2, label="After drying")
    ax.set_xlabel("Radial position [mm]")
    ax.set_ylabel("Moisture content [wt%]")
    ax.set_title(
        f"{result.filament.material.name} — Radial profile\n"
        f"T={result.dryer.chamber_temp:.0f}°C, "
        f"L={result.dryer.tube_length*1e3:.0f}mm, "
        f"Q={result.filament.flow_rate:.1f}mm³/s, "
        f"t={result.transit_time:.1f}s"
    )
    ax.legend()
    ax.set_xlim(0, r_mm[-1])
    ax.set_ylim(bottom=0)
    return ax


def plot_moisture_vs_time(result: DryingResult, ax: plt.Axes | None = None) -> plt.Axes:
    """Plot volume-averaged moisture over transit time."""
    if ax is None:
        _, ax = plt.subplots()

    ax.plot(result.diffusion.t, result.diffusion.C_avg * 100, "b-", linewidth=2)
    ax.axhline(result.initial_moisture * 100, color="0.5", linestyle="--", label="Initial")
    ax.set_xlabel("Time in dryer [s]")
    ax.set_ylabel("Average moisture [wt%]")
    ax.set_title(f"{result.filament.material.name} — Moisture vs time")
    ax.legend()
    return ax


def plot_sweep(
    param_values: NDArray,
    moisture_values: NDArray,
    param_name: str,
    param_unit: str,
    material_name: str,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot a 1D parameter sweep."""
    if ax is None:
        _, ax = plt.subplots()

    ax.plot(param_values, moisture_values * 100, "o-", markersize=3)
    ax.set_xlabel(f"{param_name} [{param_unit}]")
    ax.set_ylabel("Final moisture [wt%]")
    ax.set_title(f"{material_name} — Final moisture vs {param_name}")
    return ax


def plot_heatmap(
    x_values: NDArray,
    y_values: NDArray,
    moisture_grid: NDArray,
    x_label: str,
    y_label: str,
    title: str,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot a 2D heatmap of drying efficiency over two swept parameters."""
    if ax is None:
        _, ax = plt.subplots()

    X, Y = np.meshgrid(x_values, y_values)
    pcm = ax.pcolormesh(X, Y, moisture_grid * 100, shading="auto", cmap="viridis_r")
    plt.colorbar(pcm, ax=ax, label="Final moisture [wt%]")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    return ax


def plot_material_comparison(
    results: Sequence[DryingResult],
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Compare radial profiles of multiple materials on one plot."""
    if ax is None:
        _, ax = plt.subplots()

    for res in results:
        r_mm = res.diffusion.r * 1e3
        C_final = res.diffusion.C[:, -1]
        label = f"{res.filament.material.name} ({res.drying_efficiency*100:.1f}% removed)"
        ax.plot(r_mm, C_final * 100, linewidth=2, label=label)

    ax.set_xlabel("Radial position [mm]")
    ax.set_ylabel("Moisture content [wt%]")
    ax.set_title("Material comparison — Radial profiles after drying")
    ax.legend()
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    return ax
