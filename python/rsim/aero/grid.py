"""Independent-variable grids for aerodynamic database generation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class GridAxis:
    """Strictly increasing finite values for one interpolation dimension."""

    values: tuple[float, ...]

    def __init__(self, values: Iterable[float]):
        array = np.asarray(tuple(values), dtype=float)
        if array.ndim != 1 or array.size == 0:
            raise ValueError("a grid axis must be a nonempty one-dimensional sequence")
        if not np.all(np.isfinite(array)):
            raise ValueError("grid-axis values must be finite")
        if np.any(np.diff(array) <= 0.0):
            raise ValueError("grid-axis values must be strictly increasing")
        object.__setattr__(self, "values", tuple(float(v) for v in array))

    @classmethod
    def evenly_spaced(
        cls, start: float, stop: float, spacing: float
    ) -> "GridAxis":
        """Create an inclusive axis with user-specified point spacing.

        ``stop-start`` must be an integer multiple of ``spacing`` so every
        interval has exactly the requested width; silent short final intervals
        are rejected because they are undesirable in interpolation tables.
        """

        if spacing <= 0.0 or stop < start:
            raise ValueError("spacing must be positive and stop must be >= start")
        intervals = (stop - start) / spacing
        rounded = round(intervals)
        if not np.isclose(intervals, rounded, rtol=1.0e-10, atol=1.0e-12):
            raise ValueError("axis range must be an integer multiple of spacing")
        return cls(np.linspace(start, stop, rounded + 1))

    def array(self) -> np.ndarray:
        """Return a new NumPy array containing the axis values."""

        return np.asarray(self.values, dtype=float)


@dataclass(frozen=True)
class AeroGrid:
    """Mach, incidence, aerodynamic-roll-angle, and Reynolds dimensions.

    The guide's low-angle baseline is approximately ``|alpha| <= 10 deg``.
    Reynolds number is based on total body length. Mach must be positive and
    Reynolds number must exceed the range where continuum/flat-plate methods
    become meaningless; this module enforces only positivity and documents
    method-specific limitations in method docstrings and table metadata.
    """

    mach: GridAxis
    alpha_deg: GridAxis
    reynolds: GridAxis
    phi_aero_deg: GridAxis = GridAxis([0.0])

    def __post_init__(self) -> None:
        if min(self.mach.values) <= 0.0:
            raise ValueError("Mach values must be greater than zero")
        if min(self.reynolds.values) <= 0.0:
            raise ValueError("Reynolds values must be greater than zero")
        if min(self.alpha_deg.values) < -10.0 or max(self.alpha_deg.values) > 10.0:
            raise ValueError("the baseline implementation is limited to |alpha| <= 10 deg")
        if min(self.phi_aero_deg.values) < -180.0 or max(self.phi_aero_deg.values) > 180.0:
            raise ValueError("phi_aero must lie within [-180, 180] degrees")

    @classmethod
    def from_spacing(
        cls,
        *,
        mach: tuple[float, float, float],
        alpha_deg: tuple[float, float, float],
        reynolds: tuple[float, float, float],
        phi_aero_deg: tuple[float, float, float] = (0.0, 0.0, 1.0),
    ) -> "AeroGrid":
        """Build all dimensions from ``(start, stop, spacing)`` tuples."""

        return cls(
            mach=GridAxis.evenly_spaced(*mach),
            alpha_deg=GridAxis.evenly_spaced(*alpha_deg),
            reynolds=GridAxis.evenly_spaced(*reynolds),
            phi_aero_deg=GridAxis.evenly_spaced(*phi_aero_deg),
        )
