"""Fixed-geometry sounding-rocket aerodynamic database generation.

This package is intentionally separate from the compiled ``rsim`` extension.
It produces plain NumPy regular-grid tables that can later be consumed by the
C++ simulation without coupling database generation to runtime aerodynamics.
"""

from .database import AerodynamicOptions, AeroCoefficientTable, generate_aero_table
from .geometry import (
    ConicalBoattail,
    CylindricalBody,
    RocketGeometry,
    TangentOgiveNose,
    TrapezoidalFinSet,
)
from .grid import AeroGrid, GridAxis

__all__ = [
    "AeroCoefficientTable",
    "AeroGrid",
    "AerodynamicOptions",
    "ConicalBoattail",
    "CylindricalBody",
    "GridAxis",
    "RocketGeometry",
    "TangentOgiveNose",
    "TrapezoidalFinSet",
    "generate_aero_table",
]
