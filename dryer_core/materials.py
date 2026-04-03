"""
Material property database for hygroscopic 3D printing filaments.

Provides diffusion coefficients, activation energies, and equilibrium moisture
content for common filament materials.  Fiber-reinforced variants can be
derived at runtime via ``with_fiber()``.

Sources
-------
[1] Aniskevich, Bulderberga & Stankevics (2023). "Moisture Sorption and
    Degradation of Polymer Filaments Used in 3D Printing." Polymers
    15(12):2600.  doi:10.3390/polym15122600
[2] Chaudhary, Li & Matos (2023). "Long-term mechanical performance of 3D
    printed thermoplastics in seawater environments." Results in Materials
    17:100381.  doi:10.1016/j.rinma.2023.100381
[3] Chabaud, Castro, Denoual & Le Duigou (2019). "Hygromechanical properties
    of 3D printed continuous carbon and glass fibre reinforced polyamide
    composite for outdoor structural applications." Additive Manufacturing
    26:94-105.  doi:10.1016/j.addma.2019.01.005
[4] Sayer (2014). "Mechanical performance of polyamid 66 and influence of glass
    fiber content on moisture absorption." Materials Testing 56(4):325-330.
    doi:10.3139/120.110567
[5] Haghighi-Yazdi, Tang & Lee-Sullivan (2011). "Moisture uptake of a
    polycarbonate blend exposed to hygrothermal aging." Polymer Degradation
    and Stability 96(10):1858-1865.  doi:10.1016/j.polymdegradstab.2011.07.010
[6] Crank, J. (1975). The Mathematics of Diffusion. Oxford University Press.
[7] Manufacturer TDS.
[8] Manufacturer drying guides (Stratasys, BASF Forward AM, Polymaker, eSUN)
    and processing literature — maximum moisture content thresholds above which
    visible print defects (bubbling, stringing, poor adhesion) occur.

Arrhenius fitting
-----------------
D0 and Ea are fit coefficients in D(T) = D0 * exp(-Ea / (R*T)).

1. Gather 1-2 literature D(T) reference points per material (e.g. D at 23 °C
   from [1]).
2. Single reference point + known Ea range: fix Ea from literature, solve
   D0 = D_ref / exp(-Ea / (R * T_ref)).
3. Two reference points (D at T1, D at T2): solve
   Ea = R * ln(D1/D2) / (1/T2 - 1/T1), then back-calculate D0.
4. Fiber-reinforced variants are derived at runtime via ``with_fiber()``
   rather than stored as separate database entries.  The corrections applied
   are (see [4]):
   - Equilibrium moisture scaled by (1 - Vf), where Vf is the fiber volume
     fraction, because fibers are non-absorbing.
   - D0_composite = D0_matrix * (1 - Vf) for the tortuous diffusion path,
     while Ea ~ Ea_matrix (same polymer backbone).
   - Density follows a rule-of-mixtures.
5. All values are order-of-magnitude representative of typical commercial FFF
   filaments.  Refine with experimental data for specific brands.
"""

from dataclasses import dataclass

import numpy as np

# Universal gas constant [J/(mol·K)]
R_GAS = 8.314


@dataclass
class Material:
    """Properties of a hygroscopic filament material.

    Attributes:
        name: Human-readable material name.
        D0: Pre-exponential diffusion coefficient [m²/s] in Arrhenius relation.
        Ea: Activation energy for diffusion [J/mol].
        equilibrium_moisture: Equilibrium moisture content at saturation
            (mass fraction, e.g. 0.07 = 7 wt%).
        max_temp: Maximum recommended drying temperature [°C] (below
            degradation / glass transition onset).
        density: Approximate bulk density [kg/m³].
        max_print_moisture: Maximum moisture content for acceptable FDM print
            quality (mass fraction).  Above this threshold, steam bubbling
            and surface defects become visible at typical nozzle temperatures.
    """

    name: str
    D0: float
    Ea: float
    equilibrium_moisture: float
    max_temp: float
    density: float
    max_print_moisture: float

    def diffusivity(self, T_celsius: float) -> float:
        """Moisture diffusion coefficient at temperature T via Arrhenius.

        D(T) = D0 * exp(-Ea / (R * T))

        Args:
            T_celsius: Temperature in °C.

        Returns:
            Diffusion coefficient in m²/s.
        """
        T_kelvin = T_celsius + 273.15
        return self.D0 * np.exp(-self.Ea / (R_GAS * T_kelvin))


# Typical fiber densities [kg/m³]
_FIBER_DENSITY: dict[str, float] = {
    "glass": 2500.0,
    "carbon": 1800.0,
}


def with_fiber(
    base: Material,
    fiber_type: str,
    weight_fraction: float,
) -> Material:
    """Derive a fiber-reinforced material from a neat base polymer.

    Applies corrections per [4] (Sayer 2014) and general composite theory:

    - **Equilibrium moisture** scaled by (1 - Vf): fibers are non-absorbing.
    - **D0** scaled by (1 - Vf): tortuous diffusion path around fibers.
    - **Ea** unchanged: same polymer backbone.
    - **Density** via rule of mixtures.
    - **max_temp** unchanged.

    Weight fraction is converted to volume fraction internally.

    Args:
        base: Neat (unreinforced) base material.
        fiber_type: ``"glass"`` or ``"carbon"``.
        weight_fraction: Fiber weight fraction, e.g. 0.20 for 20 wt%.

    Returns:
        A new Material with adjusted properties.
    """
    if weight_fraction <= 0.0:
        return base

    fiber_type = fiber_type.lower()
    rho_fiber = _FIBER_DENSITY[fiber_type]
    tag = "GF" if fiber_type == "glass" else "CF"

    # Convert weight fraction to volume fraction:
    #   Wf = ρf·Vf / (ρf·Vf + ρm·(1−Vf))
    #   => Vf = (Wf·ρm) / (ρf − Wf·ρf + Wf·ρm)
    rho_m = base.density
    wf = weight_fraction
    vf = (wf * rho_m) / (rho_fiber - wf * rho_fiber + wf * rho_m)

    return Material(
        name=f"{base.name} + {weight_fraction * 100:.0f}% {tag}",
        D0=base.D0 * (1.0 - vf),
        Ea=base.Ea,
        equilibrium_moisture=base.equilibrium_moisture * (1.0 - vf),
        max_temp=base.max_temp,
        density=rho_m * (1.0 - vf) + rho_fiber * vf,
        max_print_moisture=base.max_print_moisture,
    )


# ---------------------------------------------------------------------------
# Material database
# ---------------------------------------------------------------------------
# D0 and Ea are fitted such that D(T) reproduces literature values at
# reference temperatures (typically room-temp ~1e-13 to 1e-12 m²/s, and
# elevated-temp ~1e-11 m²/s for nylon).  Values are order-of-magnitude
# representative — refine with experimental data for specific filament brands.

PA6 = Material(
    name="PA6 (Nylon 6)",
    D0=6.5e-5,  # m²/s  (fitted from [1] D_ref + [2] Ea)
    Ea=38_000.0,  # J/mol (~38 kJ/mol)  Source: [2]
    equilibrium_moisture=0.07,  # ~7 wt% at high RH  Source: [1]
    max_temp=95.0,  # °C  Source: [7]
    density=1140.0,  # kg/m³  Source: [7]
    max_print_moisture=0.002,  # 0.2 wt%  Source: [8]
)

PETG = Material(
    name="PETG",
    D0=1.2e-5,  # m²/s  (fitted from [1] D_ref + [2] Ea)
    Ea=40_000.0,  # J/mol  Source: [2]
    equilibrium_moisture=0.005,  # ~0.5 wt%  Source: [1]
    max_temp=65.0,  # °C  Source: [7]
    density=1270.0,  # kg/m³  Source: [7]
    max_print_moisture=0.002,  # 0.2 wt%  Source: [8]
)

PLA = Material(
    name="PLA",
    D0=1.0e-5,  # m²/s  (fitted from [1] D_ref + [2] Ea)
    Ea=40_000.0,  # J/mol  Source: [2]
    equilibrium_moisture=0.008,  # ~0.8 wt%  Source: [1]
    max_temp=65.0,  # °C  Source: [7]
    density=1240.0,  # kg/m³  Source: [7]
    max_print_moisture=0.003,  # 0.3 wt%  Source: [8]
)

TPU = Material(
    name="TPU",
    D0=2.0e-5,  # m²/s  Source: [7]
    Ea=38_000.0,  # J/mol  Source: [7]
    equilibrium_moisture=0.012,  # ~1.2 wt%  Source: [7]
    max_temp=70.0,  # °C  Source: [7]
    density=1210.0,  # kg/m³  Source: [7]
    max_print_moisture=0.003,  # 0.3 wt%  Source: [8]
)

PVA = Material(
    name="PVA",
    D0=8.0e-5,  # m²/s  Source: [7]
    Ea=36_000.0,  # J/mol  Source: [7]
    equilibrium_moisture=0.10,  # ~10 wt%  Source: [7]
    max_temp=50.0,  # °C  Source: [7]
    density=1190.0,  # kg/m³  Source: [7]
    max_print_moisture=0.005,  # 0.5 wt%  Source: [8]
)

ABS = Material(
    name="ABS",
    D0=0.8e-5,  # m²/s  (fitted from [1] D_ref + [2] Ea)
    Ea=38_000.0,  # J/mol  Source: [2]
    equilibrium_moisture=0.004,  # ~0.4 wt%  Source: [1]
    max_temp=100.0,  # °C  Source: [7]
    density=1040.0,  # kg/m³  Source: [7]
    max_print_moisture=0.002,  # 0.2 wt%  Source: [8]
)

PA12 = Material(
    name="PA12 (Nylon 12)",
    D0=2.0e-5,  # m²/s  Source: [1]
    Ea=40_000.0,  # J/mol  Source: [1][2]
    equilibrium_moisture=0.02,  # ~2 wt%  Source: [1]
    max_temp=80.0,  # °C  Source: [7]
    density=1010.0,  # kg/m³  Source: [7]
    max_print_moisture=0.002,  # 0.2 wt%  Source: [8]
)

PC = Material(
    name="PC (Polycarbonate)",
    D0=0.5e-5,  # m²/s  Source: [5]
    Ea=42_000.0,  # J/mol  Source: [5]
    equilibrium_moisture=0.003,  # ~0.3 wt%  Source: [5]
    max_temp=120.0,  # °C  Source: [7]
    density=1200.0,  # kg/m³  Source: [7]
    max_print_moisture=0.0002,  # 0.02 wt%  Source: [8] — hydrolysis-sensitive
)

ASA = Material(
    name="ASA",
    D0=1.0e-5,  # m²/s  Source: [2]
    Ea=38_000.0,  # J/mol  Source: [2]
    equilibrium_moisture=0.006,  # ~0.6 wt%  Source: [2]
    max_temp=95.0,  # °C  Source: [7]
    density=1070.0,  # kg/m³  Source: [7]
    max_print_moisture=0.002,  # 0.2 wt%  Source: [8]
)

PPA = Material(
    name="PPA (Polyphthalamide)",
    D0=3.5e-5,  # m²/s  Source: [2]
    Ea=42_000.0,  # J/mol  Source: [2]
    equilibrium_moisture=0.03,  # ~3 wt%  Source: [2][7]
    max_temp=120.0,  # °C  Source: [7]
    density=1180.0,  # kg/m³  Source: [7]
    max_print_moisture=0.001,  # 0.1 wt%  Source: [8]
)

# Convenience lookup — neat (unreinforced) base polymers only.
# Use ``with_fiber(base, fiber_type, weight_fraction)`` to derive
# fiber-reinforced variants at runtime.
MATERIALS: dict[str, Material] = {
    "PA6": PA6,
    "PA12": PA12,
    "PETG": PETG,
    "PLA": PLA,
    "TPU": TPU,
    "PVA": PVA,
    "ABS": ABS,
    "PC": PC,
    "ASA": ASA,
    "PPA": PPA,
}
