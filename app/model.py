"""
1D radial Fickian diffusion solver for moisture transport in cylindrical filament.

Solves the radial diffusion PDE using the method of lines (spatial finite
differences → system of ODEs) integrated with scipy's solve_ivp.

Governing equation (Fick's second law, cylindrical symmetry):

    ∂C/∂t = D/r · ∂/∂r (r · ∂C/∂r)

Boundary conditions:
    r = 0:  ∂C/∂r = 0            (symmetry)
    r = R:  C = C_env             (Dirichlet, when biot_mass is None)
         or -D·∂C/∂r = h_m·(C_s - C_env)  (Robin / convective, when biot_mass is given)

The r = 0 singularity is handled via L'Hôpital's rule:
    lim_{r→0} (1/r)·∂C/∂r = ∂²C/∂r²
    so at center: ∂C/∂t = 2D · ∂²C/∂r²
"""

from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.special import jn_zeros


@dataclass
class DiffusionResult:
    """Result of a radial diffusion simulation.

    Attributes:
        r: Radial grid points [m], shape (N,).
        t: Time points [s], shape (M,).
        C: Concentration field C(r, t), shape (N, M).
        C_avg: Volume-averaged concentration over time, shape (M,).
    """

    r: NDArray[np.floating]
    t: NDArray[np.floating]
    C: NDArray[np.floating]
    C_avg: NDArray[np.floating]


def volume_average(r: NDArray, C: NDArray) -> float:
    """Compute volume-averaged concentration in a cylinder cross-section.

    C̄ = (2/R²) ∫₀ᴿ C(r)·r·dr   (from cylindrical area element 2πr·dr)
    """
    R = r[-1]
    return 2.0 * np.trapezoid(C * r, r) / R**2


def solve_radial_diffusion(
    R: float,
    D: float,
    C_init: float,
    t_end: float,
    C_env: float = 0.0,
    N: int = 50,
    t_eval_count: int = 200,
    biot_mass: float | None = None,
) -> DiffusionResult:
    """Solve 1D radial diffusion in a cylinder cross-section.

    Args:
        R: Filament radius [m].
        D: Diffusion coefficient [m²/s] (assumed constant during transit).
        C_init: Initial uniform moisture concentration (mass fraction).
        t_end: Total diffusion time [s] (= transit time through dryer).
        C_env: Environmental moisture concentration at surface (mass fraction).
        N: Number of radial grid points (including center and surface).
        t_eval_count: Number of time points in output.
        biot_mass: Mass-transfer Biot number Bi = h_m·R/D. When None, a
            Dirichlet BC (C_surface = C_env) is used. When provided, a
            Robin (convective) BC is applied at the surface.

    Returns:
        DiffusionResult with radial profiles over time.
    """
    # Radial grid: r[0] = 0 (center), r[N-1] = R (surface)
    r = np.linspace(0, R, N)
    dr = r[1] - r[0]

    # Initial condition: uniform moisture
    C0 = np.full(N, C_init)
    if biot_mass is None:
        C0[-1] = C_env  # Dirichlet: surface fixed immediately
    # Robin: surface starts at C_init and evolves via the BC

    def rhs(t, C):
        """Right-hand side of the ODE system dC/dt = f(C)."""
        dCdt = np.empty_like(C)

        # --- Interior nodes (i = 1 .. N-2) ---
        # Standard cylindrical Laplacian: D * (C''  + (1/r)*C')
        # Using central differences:
        #   C'' ≈ (C[i+1] - 2C[i] + C[i-1]) / dr²
        #   C'  ≈ (C[i+1] - C[i-1]) / (2·dr)
        i = np.arange(1, N - 1)
        d2C = (C[i + 1] - 2 * C[i] + C[i - 1]) / dr**2
        dC = (C[i + 1] - C[i - 1]) / (2 * dr)
        dCdt[i] = D * (d2C + dC / r[i])

        # --- Center node (r = 0): L'Hôpital ---
        # ∂C/∂t = 2D · ∂²C/∂r²
        # Use one-sided second derivative: C''(0) ≈ 2(C[1] - C[0]) / dr²
        # (from symmetry: C[-1] = C[1], so the standard central diff gives this)
        dCdt[0] = 2 * D * 2 * (C[1] - C[0]) / dr**2

        # --- Surface node (r = R) ---
        if biot_mass is None:
            # Dirichlet BC: C at surface is fixed
            dCdt[-1] = 0.0
        else:
            # Robin BC via ghost-node method:
            # -D·∂C/∂r|_{r=R} = (Bi·D/R)·(C_s - C_env)
            # β = Bi/R = h_m/D  [1/m]
            beta = biot_mass / R
            dCdt[-1] = D * (
                2 * (C[-2] - C[-1]) / dr**2
                - beta * (C[-1] - C_env) * (2 / dr + 1 / R)
            )

        return dCdt

    t_eval = np.linspace(0, t_end, t_eval_count)

    sol = solve_ivp(
        rhs,
        (0, t_end),
        C0,
        method="BDF",
        t_eval=t_eval,
        rtol=1e-8,
        atol=1e-10,
    )

    if not sol.success:
        raise RuntimeError(f"ODE integration failed: {sol.message}")

    # sol.y has shape (N, len(t_eval))
    C_field = sol.y
    C_avg = np.array([volume_average(r, C_field[:, k]) for k in range(C_field.shape[1])])

    return DiffusionResult(r=r, t=sol.t, C=C_field, C_avg=C_avg)


# ---------------------------------------------------------------------------
# Analytical solution for validation (constant D, Dirichlet BC)
# ---------------------------------------------------------------------------

def analytical_moisture_fraction(
    D: float,
    R: float,
    t: float | NDArray,
    n_terms: int = 20,
) -> float | NDArray:
    """Fractional moisture remaining M(t)/M₀ for a cylinder with Dirichlet BC.

    Uses the Crank (1975) series solution:
        M(t)/M₀ = Σ (4/γₖ²) · exp(-γₖ² · D·t/R²)
    where γₖ are the positive roots of J₀(γ) = 0.
    """
    gamma = jn_zeros(0, n_terms)  # roots of J0
    Fo = D * np.asarray(t) / R**2
    # Sum over Bessel roots
    result = np.zeros_like(np.asarray(Fo, dtype=float))
    for k in range(n_terms):
        result = result + (4.0 / gamma[k] ** 2) * np.exp(-gamma[k] ** 2 * Fo)
    return result
