# Inline Filament Dryer Model (IFDM)

Numerical model for exploring the design space of an **inline filament dryer** for 3D printing. The filament passes through a heated chamber just before the extruder; this tool simulates how much moisture is removed during transit.

## Physics

1D radial Fickian diffusion in a cylindrical cross-section, solved via method of lines (`scipy.integrate.solve_ivp`, BDF). Arrhenius-type diffusivity, validated against the Crank (1975) analytical series solution.

## Current state

- Core diffusion solver and dryer simulation working.
- Material database with order-of-magnitude properties for PA6, PETG, PLA, TPU, PVA, ABS.
- Parameter sweeps (temperature, chamber length, filament speed) and heatmap generation.
- Analytical validation (< 0.5 % error with default grid).

**Key takeaway:** at typical print speeds, inline drying can only remove moisture from a thin surface shell — the model quantifies exactly how thin.

## Installation

```bash
pip install -e .            # core
pip install -e ".[dev]"      # adds ruff
```

Or all at once: `pip install -e ".[dashboard,dev]"`

## Usage

### Batch / CLI (PNG export)

```bash
python -m app   # or: make run
```

### Interactive dashboard

```bash
make dashboard  # or: streamlit run app/dashboard.py
```

Adjust chamber length, temperature, material, airflow, and flow rate with sidebar controls. Charts update in real time.

### Tasks (Makefile)

| Command | Description |
|---------|-------------|
| `make install` | Install core dependencies (editable) |
| `make run` | Run CLI simulation |
| `make dashboard` | Launch Streamlit dashboard |
| `make fmt` | Format code with ruff |
| `make lint` | Lint code with ruff |
| `make fix` | Lint + auto-fix |
| `make check` | Format + lint |

See [MANIFEST.md](MANIFEST.md) for full details on the physics, material data, and project structure.
