# Inline Filament Dryer Model (IFDM)
**[Try the live dashboard →](https://inline-dryer.streamlit.app/)**


A numerical model for exploring the design space of an **inline filament dryer** for 3D printing. The filament passes through a heated chamber just before the extruder - this tool simulates how much moisture is removed during transit.

## Physics

1D radial Fickian diffusion in a cylindrical cross-section, solved via method of lines (`scipy.integrate.solve_ivp`, BDF). Arrhenius-type diffusivity with a Robin (convective) surface boundary condition — the mass-transfer coefficient is computed from a Sherwood number correlation (Churchill–Bernstein), and chamber humidity is derived from ambient conditions via psychrometrics. Validated against the Crank (1975) analytical series solution (< 0.5 % error).

## Current state

- Diffusion solver and dryer simulation working.
- Convective mass transfer boundary condition (Robin BC) with Sherwood/Biot calculation.
- Material database with order-of-magnitude properties.
- Interactive Streamlit dashboard with real-time charts.
- Analytical validation (< 0.5 % error with default grid).
- Convergence study notebook (`convergence_study.ipynb`) — N and t_eval_count analysis.
- Optimization study notebook (`optimization_study.ipynb`) — multi-objective optimizer with convergence, timing, and cost breakdown analysis.

## Installation

```bash
pip install -e .               # core + dashboard
pip install -e ".[notebooks]"   # adds matplotlib, pandas, dlib (for notebooks)
pip install -e ".[dev]"         # adds ruff
```

## Usage

### CLI

```bash
python -m dryer_core   # or: make run
```

### Interactive dashboard

```bash
make dashboard  # or: streamlit run dashboard/app.py
# or: python -m dashboard
```

### Tasks (Makefile)

| Command | Description |
|---------|-------------|
| `make install` | Install core dependencies (editable) |
| `make install-notebooks` | Install notebook dependencies (matplotlib, pandas, dlib) |
| `make run` | Run CLI simulation (prints results to stdout) |
| `make dashboard` | Launch Streamlit dashboard |
| `make fmt` | Format code with ruff |
| `make lint` | Lint code with ruff |
| `make fix` | Lint + auto-fix |
| `make check` | Format + lint |

See [MANIFEST.md](MANIFEST.md) for full details on the physics, material data, and project structure.
