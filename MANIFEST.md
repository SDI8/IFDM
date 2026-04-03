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
dryer_core/
  __main__.py       Entry point — run with `python -m dryer_core`
  materials.py      Material dataclass + property database (PA6, PETG, PLA, …)
  model.py          Core 1D radial diffusion solver (method of lines + solve_ivp)
  dryer.py          DryerConfig, FilamentConfig, simulate(), DryingResult
                    Air properties, psychrometrics, Sherwood/Biot calculation

dashboard/
  __main__.py       Entry point — run with `python -m dashboard`
  app.py            Interactive Streamlit dashboard (Plotly charts)

convergence_study.ipynb      N and t_eval_count parameter convergence analysis
optimization_study.ipynb     Multi-objective optimization (dlib LIPO+TR global optimizer)
                             OptimizationConfig, OptimizationResult, optimize()
                             Convergence, timing, cost breakdown, multi-material comparison
```

## Material Database

Properties are order-of-magnitude representative, sourced from:
- Aniskevich, Bulderberga & Stankevics (2023), Polymers 15(12):2600
- Chaudhary, Li & Matos (2023), Results in Materials 17:100381
- Chabaud, Castro, Denoual & Le Duigou (2019), Additive Manufacturing 26:94-105
- Sayer (2014), Materials Testing 56(4):325-330
- Haghighi-Yazdi, Tang & Lee-Sullivan (2011), Polymer Degradation and Stability 96(10):1858-1865
- Manufacturer drying/processing guides (Stratasys, BASF Forward AM, Polymaker, eSUN)
- General polymer literature (Crank 1975, manufacturer TDS)

### Neat polymers

| Material | D₀ [m²/s] | Ea [kJ/mol] | Equil. moisture | Max print moisture | Max dry temp |
|----------|-----------|-------------|-----------------|--------------------|--------------| 
| PA6      | 6.5e-5    | 38          | 7.0 wt%         | 0.20 wt%           | 95 °C        |
| PA12     | 2.0e-5    | 40          | 2.0 wt%         | 0.20 wt%           | 80 °C        |
| PETG     | 1.2e-5    | 40          | 0.5 wt%         | 0.20 wt%           | 65 °C        |
| PLA      | 1.0e-5    | 40          | 0.8 wt%         | 0.30 wt%           | 65 °C        |
| TPU      | 2.0e-5    | 38          | 1.2 wt%         | 0.30 wt%           | 70 °C        |
| PVA      | 8.0e-5    | 36          | 10.0 wt%        | 0.50 wt%           | 50 °C        |
| ABS      | 0.8e-5    | 38          | 0.4 wt%         | 0.20 wt%           | 100 °C       |
| PC       | 0.5e-5    | 42          | 0.3 wt%         | 0.02 wt%           | 120 °C       |
| ASA      | 1.0e-5    | 38          | 0.6 wt%         | 0.20 wt%           | 95 °C        |
| PPA      | 3.5e-5    | 42          | 3.0 wt%         | 0.10 wt%           | 120 °C       |

**Max print moisture** is the approximate upper moisture limit for acceptable
FDM print quality (no visible bubbling/stringing). Values are from manufacturer
drying guides and processing literature.  PC is notably strict (0.02 wt%) due
to hydrolytic chain scission at its high nozzle temperature.

### Fiber-reinforced variants

Fiber-reinforced materials are **not stored as separate database entries**.
Instead, the `with_fiber(base, fiber_type, weight_fraction)` function derives
adjusted properties at runtime from any neat base polymer:

- **Equilibrium moisture** × (1 − Vf) — fibers are non-absorbing.
- **D₀** × (1 − Vf) — tortuous diffusion path around fibers.
- **Ea** unchanged — same polymer backbone.
- **Density** via rule of mixtures: ρ_matrix·(1−Vf) + ρ_fiber·Vf  
  (glass ≈ 2500 kg/m³, carbon ≈ 1800 kg/m³).
- Weight-to-volume fraction conversion is handled internally.

The dashboard exposes fiber type (None / Glass / Carbon) and fiber content
(5–40 wt%) sliders for interactive exploration.

Refine with experimental data for specific filament brands as needed.

### Arrhenius Fitting

D₀ and Ea are not directly measured — they are fit coefficients in the
Arrhenius relation  D(T) = D₀ · exp(−Ea / (R·T)).

**Procedure:**
1. Gather 1–2 literature D(T) reference points per material (e.g. D at 23 °C
   from Aniskevich 2023).
2. **Single reference point + known Ea range:** fix Ea from literature, solve
   D₀ = D_ref / exp(−Ea/(R·T_ref)).
3. **Two reference points** (D at T₁, D at T₂): solve  Ea = R·ln(D₁/D₂) /
   (1/T₂ − 1/T₁),  then back-calculate D₀.
4. **Fiber-reinforced variants** are derived at runtime via `with_fiber()`
   rather than stored as separate entries.  The corrections are:
   - Equilibrium moisture scaled by (1 − Vf) — fibers are non-absorbing.
   - D_composite ≈ D_matrix · (1 − Vf) for the tortuous diffusion path;
     Ea ≈ Ea_matrix (same polymer backbone).
   - Density via rule of mixtures.
5. All values are order-of-magnitude representative of typical commercial FFF
   filaments.  Specific brands and fiber loadings (CF ~15–20 wt%, GF ~30 wt%)
   may differ.

## Dependencies

Defined in `pyproject.toml`.

**Core:** numpy, scipy, streamlit, plotly
**Notebooks:** matplotlib, pandas, dlib-bin
**Dev:** ruff

```bash
pip install -e .               # core + dashboard
pip install -e ".[notebooks]"   # + notebook dependencies
pip install -e ".[dev]"         # + dev tools
```

## Usage

```bash
python -m app    # or: make run
```

Runs the baseline simulation (PA6, 1.75mm, 8mm³/s flow rate, 80°C, 500mm
chamber) and prints a summary with analytical validation.

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

## Optimization

The `optimization_study.ipynb` notebook finds optimal dryer parameters
by minimizing a weighted multi-objective cost using `dlib.find_min_global` with 
analysis of convergence, timing, and cost component behaviour.

**Decision variables:** chamber_length, chamber_temp, airflow_velocity.

**Cost function** — weighted sum of three normalized components:
- **Moisture** (weight 1.0): `final_moisture / initial_moisture` — fraction remaining.
- **Length** (weight 0.1): `chamber_length / max_length` — penalizes large chambers.
- **Energy** (weight 0.1): `(ΔT/ΔT_max) × (v/v_max)²` — proxy for heating
  (∝ ΔT) and fan power (∝ v²), normalized to [0, 1].

An optional `target_moisture` soft constraint adds a large penalty (100×) when
the final moisture exceeds the target.  When no explicit target is set, the
optimizer defaults to `material.max_print_moisture` — the maximum moisture
content for acceptable FDM print quality.  When no explicit target is set, the
optimizer defaults to `material.max_print_moisture` — the maximum moisture
content for acceptable FDM print quality.

Bounds are user-configurable; the temperature upper bound is automatically
clamped to `material.max_temp`.

## Planned Enhancements

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
6. Chaudhary, B., Li, H., & Matos, H. (2023). "Long-term mechanical performance
   of 3D printed thermoplastics in seawater environments." Results in Materials
   17:100381.
   — Activation energies (exponential-fit method), mass saturation, and Tg for
   Nylon, Nylon+CF, ABS, ABS+CF, PLA, PCTG, PETG, ASA.
7. Chabaud, G., Castro, M., Denoual, C., & Le Duigou, A. (2019).
   "Hygromechanical properties of 3D printed continuous carbon and glass fibre
   reinforced polyamide composite for outdoor structural applications." Additive
   Manufacturing 26:94-105.
   — Moisture uptake (5–6 wt% at 98 % RH) for continuous CF/PA and GF/PA
   composites.
8. Sayer, S. (2014). "Mechanical performance of polyamid 66 and influence of
   glass fiber content on moisture absorption." Materials Testing 56(4):325-330.
   — Glass-fiber volume-fraction proportional reduction of equilibrium moisture
   content in polyamide composites.
9. Haghighi-Yazdi, M., Tang, J.K.Y., & Lee-Sullivan, P. (2011). "Moisture uptake
   of a polycarbonate blend exposed to hygrothermal aging." Polymer Degradation
   and Stability 96(10):1858-1865.
   — Diffusion coefficients and equilibrium moisture for polycarbonate.
10. Manufacturer drying/processing guides (Stratasys, BASF Forward AM, Polymaker,
    eSUN).
    — Maximum moisture content thresholds for acceptable FDM print quality
    (bubbling, stringing, poor layer adhesion).
10. Manufacturer drying/processing guides (Stratasys, BASF Forward AM, Polymaker,
    eSUN).
    — Maximum moisture content thresholds for acceptable FDM print quality
    (bubbling, stringing, poor layer adhesion).
