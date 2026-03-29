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
- **Surface BC** at r = R: Robin (convective) boundary condition:

      -D · ∂C/∂r |_{r=R} = h_m · (C_s - C_env)

  where h_m is the convective mass transfer coefficient derived from a
  Sherwood number correlation (Churchill–Bernstein for crossflow over a
  cylinder). The mass-transfer Biot number Bi_m = h_m·R/D controls how
  close the surface BC is to the ideal Dirichlet limit (Bi → ∞).
- **Surface equilibrium**: C_env = M_sat × RH_chamber, using a linear
  sorption isotherm. Chamber RH is computed from ambient conditions via
  psychrometrics (constant absolute humidity, rising saturation pressure).
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
At high Bi (typical for forced convection), the Robin BC reduces to the
Dirichlet limit and agreement with the analytical solution is <0.5% with
50 radial grid points.

## Project Structure

```
app/
  __main__.py       Entry point — run with `python -m app`
  materials.py      Material dataclass + property database (PA6, PETG, PLA, …)
  model.py          Core 1D radial diffusion solver (method of lines + solve_ivp)
  dryer.py          DryerConfig, FilamentConfig, simulate(), DryingResult
                    Air properties, psychrometrics, Sherwood/Biot calculation
  plotting.py       Matplotlib visualizations (profiles, sweeps, heatmaps)
  experiments.py    Parameter sweep and material comparison helpers
  dashboard.py      Interactive Streamlit dashboard (Plotly charts)
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

**Core:** numpy, scipy, matplotlib, streamlit, plotly
**Dev:** ruff

```bash
pip install -e .               # core + dashboard
pip install -e ".[dev]"         # + dev tools
```

## Usage

```bash
python -m app    # or: make run
```

Runs the baseline simulation (PA6, 1.75mm, 8mm³/s flow rate, 80°C, 500mm
chamber), prints a summary with analytical validation, and generates diagnostic
plots:
1. Radial moisture profile (before/after)
2. Volume-averaged moisture vs time
3. Temperature sweep
4. Chamber length sweep
5. 2D heatmap (chamber length × flow rate)
6. Material comparison (PA6 vs PETG vs PLA)

The interactive **Streamlit dashboard** (`make dashboard`) provides real-time
controls for all parameters (chamber length, temperature, airflow velocity,
material, flow rate, ambient conditions) with Plotly charts including a
filament cross-section moisture heatmap and material comparison view.

## Key Findings / Reality Check

For PA6 at 80°C (D ≈ 1.6×10⁻¹⁰ m²/s), drying the center of a 1.75mm filament
requires Fo ~ 1, i.e., t ~ R²/D ≈ 4,900 s (~1.4 hours). At 8 mm³/s volumetric
flow rate (3.3 mm/s linear speed), a 500mm chamber gives ~150s transit
(Fo ≈ 0.03). The baseline simulation removes ~1.0 wt% (from 3.0% to 2.0%,
33.5% drying efficiency), but this is predominantly surface-layer drying. The
model quantifies the radial moisture profile and what parameter combinations
improve penetration depth.

## Planned Enhancements

### Optimization

Use scipy.optimize to find optimal (chamber_length, temperature, speed) given
constraints (max material temp, target moisture removal, practical chamber length).

### Concentration-dependent diffusivity

Some materials show D = D(C). The method-of-lines approach handles this
naturally; needs per-material data.

## References

1. Aniskevich, A., Bulderberga, O., & Stankevics, L. (2023). "Moisture Sorption
   and Degradation of Polymer Filaments Used in 3D Printing." Polymers 15(12):2600.
   — Primary source for filament material properties (D₀, Ea, equilibrium
   moisture content) in the material database.
2. Crank, J. (1975). The Mathematics of Diffusion. Oxford University Press.
   — Theoretical basis for the 1D radial Fickian diffusion model and the
   analytical Bessel-function series solution used for numerical validation.
3. Churchill, S. W. & Bernstein, M. (1977). "A Correlating Equation for Forced
   Convection from Gases and Liquids to a Circular Cylinder in Crossflow."
   J. Heat Transfer 99(2):300–306.
   — Sherwood number correlation (adapted from the Nusselt form) used to
   estimate the convective mass transfer coefficient h_m at the filament surface.
4. Alduchov, O. A. & Eskridge, R. E. (1996). "Improved Magnus Form Approximation
   of Saturation Vapor Pressure." J. Appl. Meteor. 35(4):601–609.
   — Saturation vapor pressure formula used to compute chamber RH from
   ambient conditions (psychrometric heating at constant absolute humidity).
5. Fuller, E. N., Schettler, P. D. & Giddings, J. C. (1966). "A New Method for
   Prediction of Binary Gas-Phase Diffusion Coefficients." Ind. Eng. Chem.
   58(5):18–27.
   — Temperature-dependent correlation for water vapor diffusivity in air,
   needed to evaluate the Schmidt and Sherwood numbers.
