"""
Multi-objective optimization for inline filament dryer design.

Finds optimal (chamber_length, chamber_temp, airflow_velocity) by minimizing
a weighted cost of moisture remaining, chamber length, and energy usage.
Uses scipy.optimize.differential_evolution (global, gradient-free).
"""

from dataclasses import dataclass, field

from scipy.optimize import differential_evolution

from .dryer import DryerConfig, DryingResult, FilamentConfig, simulate


@dataclass
class OptimizationConfig:
    """Configuration for multi-objective dryer optimization.

    Attributes:
        filament: Fixed filament configuration (material, diameter, moisture, flow_rate).
        bounds_length: Min/max chamber length [m].
        bounds_temp: Min/max chamber temperature [°C]. Upper bound is clamped
            to material.max_temp automatically.
        bounds_airflow: Min/max airflow velocity [m/s].
        weight_moisture: Cost weight for moisture remaining fraction.
        weight_length: Cost weight for normalized chamber length.
        weight_energy: Cost weight for energy proxy.
        target_moisture: If set, final moisture above this incurs a large penalty.
        ambient_humidity: Ambient relative humidity (0–1).
        ambient_temp: Ambient temperature [°C].
    """

    filament: FilamentConfig
    bounds_length: tuple[float, float] = (0.1, 5.0)
    bounds_temp: tuple[float, float] = (40.0, 120.0)
    bounds_airflow: tuple[float, float] = (0.1, 10.0)
    weight_moisture: float = 1.0
    weight_length: float = 0.1
    weight_energy: float = 0.1
    target_moisture: float | None = None
    ambient_humidity: float = 0.50
    ambient_temp: float = 25.0


@dataclass
class OptimizationResult:
    """Result of a dryer optimization run.

    Attributes:
        optimal_length: Best chamber length [m].
        optimal_temp: Best chamber temperature [°C].
        optimal_airflow: Best airflow velocity [m/s].
        optimal_cost: Total weighted cost at the optimum.
        cost_breakdown: Individual cost components (pre-weighting).
        simulation: Full DryingResult at the optimal point.
        config: The OptimizationConfig used.
    """

    optimal_length: float
    optimal_temp: float
    optimal_airflow: float
    optimal_cost: float
    cost_breakdown: dict[str, float] = field(default_factory=dict)
    simulation: DryingResult | None = None
    config: OptimizationConfig | None = None

    def summary(self) -> str:
        lines = [
            f"Optimal chamber length:  {self.optimal_length * 1e3:.0f} mm",
            f"Optimal chamber temp:    {self.optimal_temp:.1f} °C",
            f"Optimal airflow:         {self.optimal_airflow:.2f} m/s",
            "",
            f"Total cost:              {self.optimal_cost:.6f}",
        ]
        for name, value in self.cost_breakdown.items():
            lines.append(f"  {name:25s} {value:.6f}")
        if self.simulation:
            sim = self.simulation
            lines += [
                "",
                f"Final moisture:          {sim.final_moisture * 100:.3f} wt%",
                f"Drying efficiency:       {sim.drying_efficiency * 100:.2f} %",
                f"Transit time:            {sim.transit_time:.2f} s",
                f"Fourier number:          {sim.fourier_number:.4e}",
            ]
        return "\n".join(lines)


def _objective(x: list[float], config: OptimizationConfig) -> float:
    """Scalar cost function for the optimizer.

    Cost components (each normalized to [0, 1]):
      - moisture: final_moisture / initial_moisture  (fraction remaining)
      - length:   chamber_length / max_length
      - energy:   (ΔT/ΔT_max) × (v/v_max)²
                   Heating power ∝ ΔT, fan power ∝ v².
    """
    chamber_length, chamber_temp, airflow_velocity = x

    dryer = DryerConfig(
        chamber_length=chamber_length,
        chamber_temp=chamber_temp,
        ambient_humidity=config.ambient_humidity,
        ambient_temp=config.ambient_temp,
        airflow_velocity=airflow_velocity,
    )

    result = simulate(dryer, config.filament, N=30, t_eval_count=50)

    # --- Normalized cost components ---
    cost_moisture = result.final_moisture / result.initial_moisture

    cost_length = chamber_length / config.bounds_length[1]

    delta_t = chamber_temp - config.ambient_temp
    delta_t_max = config.bounds_temp[1] - config.ambient_temp
    v_norm = airflow_velocity / config.bounds_airflow[1]
    cost_energy = (delta_t / delta_t_max) * v_norm**2

    cost = (
        config.weight_moisture * cost_moisture
        + config.weight_length * cost_length
        + config.weight_energy * cost_energy
    )

    # Soft penalty for target moisture constraint
    if config.target_moisture is not None and result.final_moisture > config.target_moisture:
        overshoot = result.final_moisture - config.target_moisture
        cost += 100.0 * overshoot

    return cost


def optimize(config: OptimizationConfig) -> OptimizationResult:
    """Run multi-objective optimization to find best dryer parameters.

    Uses differential_evolution (global, gradient-free optimizer).
    Decision variables: [chamber_length, chamber_temp, airflow_velocity].
    """
    # Clamp temperature upper bound to material safe limit
    temp_upper = min(config.bounds_temp[1], config.filament.material.max_temp)
    bounds = [
        config.bounds_length,
        (config.bounds_temp[0], temp_upper),
        config.bounds_airflow,
    ]

    result = differential_evolution(
        _objective,
        bounds=bounds,
        args=(config,),
        seed=42,
        tol=0.01,
        maxiter=1000,
        polish=True,
    )

    opt_length, opt_temp, opt_airflow = result.x

    # Re-simulate at optimum with full resolution for detailed result
    dryer = DryerConfig(
        chamber_length=opt_length,
        chamber_temp=opt_temp,
        ambient_humidity=config.ambient_humidity,
        ambient_temp=config.ambient_temp,
        airflow_velocity=opt_airflow,
    )
    sim = simulate(dryer, config.filament)

    # Recompute cost breakdown for reporting
    cost_moisture = sim.final_moisture / sim.initial_moisture
    cost_length = opt_length / config.bounds_length[1]
    delta_t = opt_temp - config.ambient_temp
    delta_t_max = config.bounds_temp[1] - config.ambient_temp
    v_norm = opt_airflow / config.bounds_airflow[1]
    cost_energy = (delta_t / delta_t_max) * v_norm**2

    breakdown = {
        f"moisture  (×{config.weight_moisture})": config.weight_moisture * cost_moisture,
        f"length   (×{config.weight_length})": config.weight_length * cost_length,
        f"energy   (×{config.weight_energy})": config.weight_energy * cost_energy,
    }

    return OptimizationResult(
        optimal_length=opt_length,
        optimal_temp=opt_temp,
        optimal_airflow=opt_airflow,
        optimal_cost=result.fun,
        cost_breakdown=breakdown,
        simulation=sim,
        config=config,
    )
