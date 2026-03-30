"""
Material property database for hygroscopic 3D printing filaments.

Provides diffusion coefficients, activation energies, and equilibrium moisture
content for common filament materials. Properties primarily sourced from:

  Aniskevich, Bulderberga & Stankevics (2023). "Moisture Sorption and
  Degradation of Polymer Filaments Used in 3D Printing." Polymers 15(12):2600

  Crank, J. (1975). The Mathematics of Diffusion. Oxford University Press.
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
    """

    name: str
    D0: float
    Ea: float
    equilibrium_moisture: float
    max_temp: float
    density: float

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


# ---------------------------------------------------------------------------
# Material database
# ---------------------------------------------------------------------------
# D0 and Ea are fitted such that D(T) reproduces literature values at
# reference temperatures (typically room-temp ~1e-13 to 1e-12 m²/s, and
# elevated-temp ~1e-11 m²/s for nylon).  Values are order-of-magnitude
# representative — refine with experimental data for specific filament brands.

PA6 = Material(
    name="PA6 (Nylon 6)",
    D0=6.5e-5,  # m²/s  (pre-exponential, fitted)
    Ea=38_000.0,  # J/mol (~38 kJ/mol)
    equilibrium_moisture=0.07,  # ~7 wt% at high RH
    max_temp=95.0,  # °C
    density=1140.0,  # kg/m³
)

PETG = Material(
    name="PETG",
    D0=1.2e-5,  # m²/s  (pre-exponential, fitted)
    Ea=40_000.0,  # J/mol (~40 kJ/mol)
    equilibrium_moisture=0.005,  # ~0.5 wt%
    max_temp=65.0,  # °C
    density=1270.0,  # kg/m³
)

PLA = Material(
    name="PLA",
    D0=1.0e-5,
    Ea=40_000.0,
    equilibrium_moisture=0.008,  # ~0.8 wt%
    max_temp=65.0,
    density=1240.0,
)

TPU = Material(
    name="TPU",
    D0=2.0e-5,
    Ea=38_000.0,
    equilibrium_moisture=0.012,  # ~1.2 wt%
    max_temp=70.0,
    density=1210.0,
)

PVA = Material(
    name="PVA",
    D0=8.0e-5,
    Ea=36_000.0,
    equilibrium_moisture=0.10,  # ~10 wt%
    max_temp=50.0,
    density=1190.0,
)

ABS = Material(
    name="ABS",
    D0=0.8e-5,
    Ea=38_000.0,
    equilibrium_moisture=0.004,  # ~0.4 wt%
    max_temp=100.0,
    density=1040.0,
)

# Convenience lookup
MATERIALS: dict[str, Material] = {
    "PA6": PA6,
    "PETG": PETG,
    "PLA": PLA,
    "TPU": TPU,
    "PVA": PVA,
    "ABS": ABS,
}
