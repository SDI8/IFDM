# Inline Filament Dryer Model (IFDM)

Numerical model for exploring the design space of an **inline filament dryer** for 3D printing. The filament passes through a heated tube just before the extruder; this tool simulates how much moisture is removed during transit.

## Physics

1D radial Fickian diffusion in a cylindrical cross-section, solved via method of lines (`scipy.integrate.solve_ivp`, BDF). Arrhenius-type diffusivity, validated against the Crank (1975) analytical series solution.

## Current state

- Core diffusion solver and dryer simulation working.
- Material database with order-of-magnitude properties for PA6, PETG, PLA, TPU, PVA, ABS.
- Parameter sweeps (temperature, tube length, filament speed) and heatmap generation.
- Analytical validation (< 0.5 % error with default grid).

**Key takeaway:** at typical print speeds, inline drying can only remove moisture from a thin surface shell — the model quantifies exactly how thin.

## Usage

```bash
pip install numpy scipy matplotlib
python -m app
```

See [MANIFEST.md](MANIFEST.md) for full details on the physics, material data, and project structure.
