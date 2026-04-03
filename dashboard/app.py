"""
Interactive Streamlit dashboard for inline filament dryer simulation.

Launch with:
    streamlit run dashboard/app.py
    python -m dashboard
"""

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from dryer_core.dryer import DryerConfig, FilamentConfig, simulate
from dryer_core.materials import MATERIALS, with_fiber

st.set_page_config(page_title="Inline Filament Dryer", layout="wide")
st.title("Inline Filament Dryer Simulator")

# ---------------------------------------------------------------------------
# Sidebar controls
# ---------------------------------------------------------------------------
st.sidebar.header("Filament parameters")

material_key = st.sidebar.selectbox(
    "Material",
    options=list(MATERIALS.keys()),
    index=0,
)
material = MATERIALS[material_key]

fiber_type = st.sidebar.selectbox(
    "Fiber reinforcement",
    options=["None", "Glass Fiber (GF)", "Carbon Fiber (CF)"],
    index=0,
)
fiber_wt_pct = 0.0
if fiber_type != "None":
    fiber_wt_pct = st.sidebar.slider(
        "Fiber content [wt%]",
        min_value=5,
        max_value=40,
        value=20,
        step=5,
    )
    ft = "glass" if "Glass" in fiber_type else "carbon"
    material = with_fiber(material, ft, fiber_wt_pct / 100.0)

initial_moisture = st.sidebar.slider(
    "Initial moisture [wt%]",
    min_value=0.1,
    max_value=material.equilibrium_moisture * 100,
    value=min(3.0, material.equilibrium_moisture * 100),
    step=0.1,
)
flow_rate = st.sidebar.slider(
    "Flow rate [mm³/s]",
    min_value=1.0,
    max_value=30.0,
    value=8.0,
    step=0.5,
)

st.sidebar.header("Dryer parameters")

chamber_length = st.sidebar.slider(
    "Chamber length [mm]",
    min_value=100,
    max_value=10000,
    value=500,
    step=100,
)
chamber_temp = st.sidebar.slider(
    "Chamber temperature [°C]",
    min_value=30.0,
    max_value=material.max_temp + 20.0,
    value=material.max_temp,
    step=5.0,
)
airflow_velocity = st.sidebar.slider(
    "Airflow velocity [m/s]",
    min_value=0.05,
    max_value=5.0,
    value=1.0,
    step=0.05,
)

st.sidebar.header("Environment")

ambient_humidity = st.sidebar.slider(
    "Ambient RH [%]",
    min_value=10,
    max_value=90,
    value=50,
    step=5,
)
ambient_temp = st.sidebar.slider(
    "Ambient temperature [°C]",
    min_value=10.0,
    max_value=40.0,
    value=25.0,
    step=1.0,
)

# ---------------------------------------------------------------------------
# Run simulation
# ---------------------------------------------------------------------------

dryer = DryerConfig(
    chamber_length=chamber_length / 1000.0,
    chamber_temp=chamber_temp,
    ambient_humidity=ambient_humidity / 100.0,
    ambient_temp=ambient_temp,
    airflow_velocity=airflow_velocity,
)
filament = FilamentConfig(
    material=material,
    initial_moisture=initial_moisture / 100.0,
    flow_rate=flow_rate,
)

result = simulate(dryer, filament)

# ---------------------------------------------------------------------------
# Key metrics
# ---------------------------------------------------------------------------

col1, col2, col3, col4, col5 = st.columns(5)
col1.metric("Final moisture", f"{result.final_moisture * 100:.3f} wt%")
col2.metric("Drying efficiency", f"{result.drying_efficiency * 100:.1f}%")
col3.metric("Transit time", f"{result.transit_time:.1f} s")
col4.metric("Fourier number", f"{result.fourier_number:.2e}")
col5.metric("Biot (mass)", f"{result.biot_mass:.2f}")

# Printability check
threshold_pct = material.max_print_moisture * 100
if result.final_moisture > material.max_print_moisture:
    st.warning(
        f"Final moisture ({result.final_moisture * 100:.2f} wt%) exceeds "
        f"max print threshold ({threshold_pct:.2f} wt%) for {material.name}."
    )
else:
    st.success(
        f"Final moisture ({result.final_moisture * 100:.2f} wt%) is below "
        f"max print threshold ({threshold_pct:.2f} wt%) for {material.name}."
    )

# ---------------------------------------------------------------------------
# Charts
# ---------------------------------------------------------------------------

left, right = st.columns(2)

# --- Radial moisture profile ---
r_mm = result.diffusion.r * 1e3
C_init = result.diffusion.C[:, 0] * 100
C_final = result.diffusion.C[:, -1] * 100

fig_radial = go.Figure()
fig_radial.add_trace(
    go.Scatter(
        x=r_mm,
        y=C_init,
        mode="lines",
        name="Initial",
        line=dict(dash="dash", color="gray"),
    )
)
fig_radial.add_trace(
    go.Scatter(
        x=r_mm,
        y=C_final,
        mode="lines",
        name="After drying",
        line=dict(width=2, color="#1f77b4"),
    )
)
fig_radial.update_layout(
    title="Radial moisture profile",
    xaxis_title="Radial position [mm]",
    yaxis_title="Moisture content [wt%]",
    yaxis_rangemode="tozero",
    height=400,
    margin=dict(t=40, b=40),
)
fig_radial.add_hline(
    y=threshold_pct,
    line_dash="dot",
    line_color="red",
    annotation_text="Max print moisture",
    annotation_position="top left",
)
left.plotly_chart(fig_radial, width="stretch")

# --- Moisture vs time ---
fig_time = go.Figure()
fig_time.add_trace(
    go.Scatter(
        x=result.diffusion.t,
        y=result.diffusion.C_avg * 100,
        mode="lines",
        name="Avg moisture",
        line=dict(width=2, color="#1f77b4"),
    )
)
fig_time.add_hline(
    y=initial_moisture,
    line_dash="dash",
    line_color="gray",
    annotation_text="Initial",
)
fig_time.update_layout(
    title="Volume-averaged moisture vs time",
    xaxis_title="Time in dryer [s]",
    yaxis_title="Average moisture [wt%]",
    height=400,
    margin=dict(t=40, b=40),
)
fig_time.add_hline(
    y=threshold_pct,
    line_dash="dot",
    line_color="red",
    annotation_text="Max print moisture",
    annotation_position="top left",
)
right.plotly_chart(fig_time, width="stretch")

# --- Cross-section moisture heatmap ---
R_mm = r_mm[-1]
n_px = 200  # resolution of the Cartesian grid
xy = np.linspace(-R_mm, R_mm, n_px)
X, Y = np.meshgrid(xy, xy)
R_grid = np.sqrt(X**2 + Y**2)

# Interpolate 1-D radial profile onto 2-D Cartesian grid
Z = np.interp(R_grid, r_mm, C_final, right=np.nan)
# Mask outside the filament circle
Z[R_grid > R_mm] = np.nan

fig_cross = go.Figure(
    data=go.Heatmap(
        x=xy,
        y=xy,
        z=Z,
        colorscale="Blues",
        colorbar_title="wt%",
        zmin=0,
        zmax=float(C_init.max()),
        hovertemplate="x: %{x:.3f} mm<br>y: %{y:.3f} mm<br>Moisture: %{z:.3f} wt%<extra></extra>",
    )
)
fig_cross.update_layout(
    title="Filament cross-section — moisture after drying",
    xaxis_title="[mm]",
    yaxis_title="[mm]",
    xaxis=dict(scaleanchor="y", constrain="domain"),
    yaxis=dict(constrain="domain"),
    height=480,
    margin=dict(t=40, b=40),
)
st.plotly_chart(fig_cross, width="stretch")


# ---------------------------------------------------------------------------
# Detailed summary
# ---------------------------------------------------------------------------

with st.expander("Simulation details"):
    st.code(result.summary(), language=None)
