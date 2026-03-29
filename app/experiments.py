"""
Parameter sweep and comparison experiments.

Provides helpers to sweep one or two parameters while holding others fixed,
producing arrays of results suitable for plotting.
"""

from dataclasses import replace

import numpy as np
from numpy.typing import NDArray

from .dryer import DryerConfig, DryingResult, FilamentConfig, simulate
from .materials import Material


def sweep_chamber_length(
    lengths: NDArray,
    dryer: DryerConfig,
    filament: FilamentConfig,
) -> tuple[NDArray, list[DryingResult]]:
    """Sweep chamber length, returning final moisture for each value."""
    moistures = np.empty(len(lengths))
    results = []
    for i, L in enumerate(lengths):
        cfg = replace(dryer, chamber_length=L)
        res = simulate(cfg, filament)
        moistures[i] = res.final_moisture
        results.append(res)
    return moistures, results


def sweep_temperature(
    temperatures: NDArray,
    dryer: DryerConfig,
    filament: FilamentConfig,
) -> tuple[NDArray, list[DryingResult]]:
    """Sweep chamber temperature, returning final moisture for each value."""
    moistures = np.empty(len(temperatures))
    results = []
    for i, T in enumerate(temperatures):
        cfg = replace(dryer, chamber_temp=T)
        res = simulate(cfg, filament)
        moistures[i] = res.final_moisture
        results.append(res)
    return moistures, results


def sweep_flow_rate(
    flow_rates: NDArray,
    dryer: DryerConfig,
    filament: FilamentConfig,
) -> tuple[NDArray, list[DryingResult]]:
    """Sweep volumetric flow rate [mm³/s], returning final moisture for each value."""
    moistures = np.empty(len(flow_rates))
    results = []
    for i, q in enumerate(flow_rates):
        fcfg = FilamentConfig(
            material=filament.material,
            diameter=filament.diameter,
            initial_moisture=filament.initial_moisture,
            flow_rate=q,
        )
        res = simulate(dryer, fcfg)
        moistures[i] = res.final_moisture
        results.append(res)
    return moistures, results


def sweep_length_flow_rate(
    lengths: NDArray,
    flow_rates: NDArray,
    dryer: DryerConfig,
    filament: FilamentConfig,
) -> NDArray:
    """2D sweep over chamber length and volumetric flow rate.

    Returns:
        moisture_grid: shape (len(flow_rates), len(lengths)) — final moisture
            for each (flow_rate, length) pair. Row index = flow_rate, col = length.
    """
    grid = np.empty((len(flow_rates), len(lengths)))
    for j, q in enumerate(flow_rates):
        for i, L in enumerate(lengths):
            cfg = replace(dryer, chamber_length=L)
            fcfg = FilamentConfig(
                material=filament.material,
                diameter=filament.diameter,
                initial_moisture=filament.initial_moisture,
                flow_rate=q,
            )
            res = simulate(cfg, fcfg)
            grid[j, i] = res.final_moisture
    return grid


def compare_materials(
    materials: list[Material],
    dryer: DryerConfig,
    diameter: float = 1.75e-3,
    initial_moisture: float | None = None,
    flow_rate: float = 8.0,
) -> list[DryingResult]:
    """Run the same dryer config for multiple materials.

    If initial_moisture is None, each material uses its equilibrium_moisture
    (worst case: fully saturated filament).
    """
    results = []
    for mat in materials:
        m0 = initial_moisture if initial_moisture is not None else mat.equilibrium_moisture
        fcfg = FilamentConfig(
            material=mat,
            diameter=diameter,
            initial_moisture=m0,
            flow_rate=flow_rate,
        )
        results.append(simulate(dryer, fcfg))
    return results
