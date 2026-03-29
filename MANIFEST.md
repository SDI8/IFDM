# Inline Filament Dryer — Project Manifest

## Purpose

Numerical model for exploring the design parameters of an inline filament dryer
for 3D printing. The filament passes through a heated chamber immediately before
entering the printer; this tool simulates how much moisture is removed during
transit, helping find practical operating parameters.

## Physics Model

**1D radial Fickian diffusion** in a cylindrical filament cross-section.

Each "slice" of filament spends `t_transit = chamber_length / filament_speed` in
the dryer. During that time, moisture diffuses radially outward:

    ∂C/∂t = D(T)/r · ∂/∂r(r · ∂C/∂r)

- **Symmetry BC** at center (r = 0): ∂C/∂r = 0
- **Surface BC** at r = R: C = C_env (Dirichlet — assumes perfectly dry air)
- **Diffusivity** follows Arrhenius: D(T) = D₀ · exp(-Ea / (R·T))
- **Thermal equilibrium** assumed: heat diffuses ~100× faster than moisture,
  so filament reaches chamber temperature almost instantly.
- Axial diffusion neglected (chamber length >> filament radius).

Solved via **method of lines** (radial finite differences → ODE system) using
`scipy.integrate.solve_ivp` with BDF (implicit stiff solver).

### Key dimensionless number

**Fourier number** Fo = D·t / R² controls drying depth:
- Fo ≪ 1: only a thin surface shell dries (typical for inline drying at print speed)
- Fo ~ 1: significant core drying
- Fo > 1: nearly complete drying

### Validation

Numerical results are compared against the **Crank (1975) analytical series
solution** for constant-D, Dirichlet-BC cylinder drying (Bessel function roots).
Agreement is typically <0.5% with 50 radial grid points.

## Project Structure

```
app/
  __main__.py       Entry point — run with `python -m app`
  materials.py      Material dataclass + property database (PA6, PETG, PLA, …)
  model.py          Core 1D radial diffusion solver (method of lines + solve_ivp)
  dryer.py          DryerConfig, FilamentConfig, simulate(), DryingResult
  plotting.py       Matplotlib visualizations (profiles, sweeps, heatmaps)
  experiments.py    Parameter sweep and material comparison helpers
```

## Material Database

Properties are order-of-magnitude representative, sourced from:
- Aniskevich, Bulderberga & Stankevics (2023), Polymers 15(12):2600
- General polymer literature (Crank 1975, various PA6/PETG studies)

| Material | D₀ [m²/s] | Ea [kJ/mol] | Equil. moisture | Max dry temp |
|----------|-----------|-------------|-----------------|--------------|
| PA6      | 6.5e-5    | 38          | 7.0 wt%         | 85 °C        |
| PETG     | 1.2e-5    | 40          | 0.5 wt%         | 65 °C        |
| PLA      | 1.0e-5    | 40          | 0.8 wt%         | 55 °C        |
| TPU      | 2.0e-5    | 38          | 1.2 wt%         | 50 °C        |
| PVA      | 8.0e-5    | 36          | 10.0 wt%        | 60 °C        |
| ABS      | 0.8e-5    | 38          | 0.4 wt%         | 80 °C        |

Refine with experimental data for specific filament brands as needed.

## Dependencies

Defined in `pyproject.toml`.

**Core:** numpy, scipy, matplotlib
**Dashboard:** streamlit, plotly
**Dev:** ruff

```bash
pip install -e .               # core
pip install -e ".[dashboard]"   # + dashboard
pip install -e ".[dev]"         # + dev tools
```

## Usage

```bash
python -m app    # or: make run
```

Runs the baseline simulation (PA6, 1.75mm, 50mm/s, 80°C, 500mm chamber), prints
a summary with analytical validation, and generates diagnostic plots:
1. Radial moisture profile (before/after)
2. Volume-averaged moisture vs time
3. Temperature sweep
4. Chamber length sweep
5. 2D heatmap (chamber length × filament speed)
6. Material comparison (PA6 vs PETG vs PLA)

## Key Findings / Reality Check

For PA6 at 80°C (D ≈ 10⁻¹¹ m²/s), drying the center of a 1.75mm filament
requires Fo ~ 1, i.e., t ~ R²/D ≈ 76,000 s (~21 hours). At 50 mm/s print
speed, a 500mm chamber gives only 10s transit (Fo ≈ 10⁻⁴). Inline drying at
typical print speeds can only dry a thin surface shell. The model quantifies
exactly how thin and what parameter combinations improve it.

## Planned Enhancements

### Convective mass transfer BC (required)

The current Dirichlet BC (C_surface = 0) assumes infinite surface mass transfer.
A Robin (convective) BC is needed:

    -D · ∂C/∂r |_{r=R} = h_m · (C_s - C_env)

where h_m comes from Sherwood number correlations for the chamber geometry:

    Sh = h_m · d / D_air_water = f(Re, Sc)

This is essential for properly modeling airflow velocity and chamber humidity
effects. Infrastructure (DryerConfig fields) is already in place.

### Optimization

Use scipy.optimize to find optimal (chamber_length, temperature, speed) given
constraints (max material temp, target moisture removal, practical chamber length).

### Concentration-dependent diffusivity

Some materials show D = D(C). The method-of-lines approach handles this
naturally; needs per-material data.

## References

1. Aniskevich, A., Bulderberga, O., & Stankevics, L. (2023). "Moisture Sorption
   and Degradation of Polymer Filaments Used in 3D Printing." Polymers 15(12):2600.
2. Crank, J. (1975). The Mathematics of Diffusion. Oxford University Press.
3. pyDAMPF (github.com/govarguz/pyDAMPF) — Reviewed; not applicable. AFM simulator
   for nanoscale tip-sample interactions under variable RH, not bulk moisture diffusion.
