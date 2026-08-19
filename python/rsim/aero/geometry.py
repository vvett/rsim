"""Fixed conventional sounding-rocket geometry definitions.

All dimensions are SI and all axial locations are measured aft from the nose.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


def _positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and greater than zero")


@dataclass(frozen=True)
class TangentOgiveNose:
    """Tangent-ogive nose with a sharp tip and circular base.

    Assumptions and range:
        The profile is axisymmetric, tangent to the following cylinder, and
        has positive length and base diameter. The aerodynamic methods are
        intended for slender sounding-rocket noses, not blunt caps.
    """

    length: float
    base_diameter: float

    def __post_init__(self) -> None:
        _positive("nose length", self.length)
        _positive("nose base diameter", self.base_diameter)

    def radius(self, x: np.ndarray | float) -> np.ndarray:
        """Return profile radius at axial location ``0 <= x <= length``.

        This is the exact tangent-ogive circle construction. It assumes a
        sharp axisymmetric ogive and is only valid within the nose interval.
        Source: standard tangent-ogive geometry, used by NWL TR-2796 Fig. 3.
        """

        locations = np.asarray(x, dtype=float)
        if np.any((locations < 0.0) | (locations > self.length)):
            raise ValueError("nose axial locations must lie within the nose")
        base_radius = 0.5 * self.base_diameter
        circle_radius = (
            self.length**2 + base_radius**2
        ) / (2.0 * base_radius)
        radicand = np.maximum(
            0.0, circle_radius**2 - (self.length - locations) ** 2
        )
        return np.sqrt(radicand) + base_radius - circle_radius

    @property
    def fineness_ratio(self) -> float:
        """Return nose length divided by base diameter."""

        return self.length / self.base_diameter


@dataclass(frozen=True)
class CylindricalBody:
    """Constant-diameter cylindrical body section."""

    length: float
    diameter: float

    def __post_init__(self) -> None:
        _positive("body length", self.length)
        _positive("body diameter", self.diameter)


@dataclass(frozen=True)
class ConicalBoattail:
    """Optional conical contraction at the aft end of the cylinder.

    Assumptions and range:
        The boattail is axisymmetric and linearly tapers from the body diameter
        to a smaller positive exit diameter. NASA TR R-100 Eq. (1), used by
        the aerodynamic module, cautions that slopes near or above 15 degrees
        can have large errors and is primarily a supersonic method.
    """

    length: float
    exit_diameter: float

    def __post_init__(self) -> None:
        _positive("boattail length", self.length)
        _positive("boattail exit diameter", self.exit_diameter)


@dataclass(frozen=True)
class TrapezoidalFinSet:
    """One fixed, evenly spaced set of thin trapezoidal aft fins.

    ``semi_span`` is exposed span from the body surface. ``tip_leading_edge_offset``
    is positive aft from the root leading edge. Thickness is the maximum fin
    thickness used by the double-wedge drag fallback.
    """

    count: int
    root_chord: float
    tip_chord: float
    semi_span: float
    tip_leading_edge_offset: float
    thickness: float
    root_leading_edge_x: float
    azimuth_offset_deg: float = 0.0

    def __post_init__(self) -> None:
        if self.count < 3:
            raise ValueError("a conventional fin set must contain at least 3 fins")
        _positive("fin root chord", self.root_chord)
        _positive("fin tip chord", self.tip_chord)
        _positive("fin semi-span", self.semi_span)
        _positive("fin thickness", self.thickness)
        if self.tip_chord > self.root_chord:
            raise ValueError("this baseline supports tip chord <= root chord")
        if self.tip_leading_edge_offset < 0.0:
            raise ValueError("forward-swept fins are outside this baseline")
        if self.root_leading_edge_x < 0.0:
            raise ValueError("fin root leading edge must be aft of the nose tip")
        if not math.isfinite(self.azimuth_offset_deg):
            raise ValueError("fin azimuth offset must be finite")

    @property
    def azimuths_rad(self) -> np.ndarray:
        """Return evenly spaced fin-plane azimuths including the fixed offset."""
        offset = math.radians(self.azimuth_offset_deg)
        return offset + 2.0 * math.pi * np.arange(self.count) / self.count

    @property
    def one_fin_area(self) -> float:
        """Return exposed planform area of one fin using guide equation FIN-01."""

        return 0.5 * (self.root_chord + self.tip_chord) * self.semi_span

    @property
    def total_area(self) -> float:
        """Return total exposed planform area for all fins."""

        return self.count * self.one_fin_area

    @property
    def taper_ratio(self) -> float:
        """Return tip chord divided by root chord using guide equation FIN-02."""

        return self.tip_chord / self.root_chord

    @property
    def mid_chord_length(self) -> float:
        """Return Barrowman mid-chord-line length using equation BARF-01.

        Assumes a straight-edged trapezoidal planform and exposed semi-span.
        """

        axial = self.tip_leading_edge_offset + 0.5 * (
            self.tip_chord - self.root_chord
        )
        return math.hypot(self.semi_span, axial)

    @property
    def aspect_ratio(self) -> float:
        """Return a fin-set aspect ratio using total tip-to-tip span squared.

        This baseline uses ``b = 2*semi_span`` and the planform area of two
        opposing exposed fins, matching a conventional planar-wing equivalent.
        Unusual fin counts or embedded/root areas require a validated convention.
        """

        full_span = 2.0 * self.semi_span
        planar_area = 2.0 * self.one_fin_area
        return full_span**2 / planar_area

    @property
    def half_chord_sweep(self) -> float:
        """Return half-chord sweep angle in radians."""

        axial = self.tip_leading_edge_offset + 0.5 * (
            self.tip_chord - self.root_chord
        )
        return math.atan2(axial, self.semi_span)

    @property
    def leading_edge_sweep(self) -> float:
        """Return leading-edge sweep angle in radians."""

        return math.atan2(self.tip_leading_edge_offset, self.semi_span)

    @property
    def thickness_ratio(self) -> float:
        """Return maximum thickness divided by root chord."""

        return self.thickness / self.root_chord

    @property
    def centroid_x(self) -> float:
        """Return the exposed trapezoid's axial area centroid."""

        cr = self.root_chord
        ct = self.tip_chord
        xt = self.tip_leading_edge_offset
        local = (cr**2 + cr * ct + ct**2 + xt * (cr + 2.0 * ct)) / (
            3.0 * (cr + ct)
        )
        return self.root_leading_edge_x + local


@dataclass(frozen=True)
class RocketGeometry:
    """Fixed conventional rocket used to generate one aerodynamic database.

    Assumptions and range:
        Axisymmetric tangent-ogive/cylinder body, optional conical boattail,
        one thin aft trapezoidal fin set, no control deflection, no
        protuberances, and a common maximum-body reference area. Geometry and
        reference quantities are SI.
    """

    nose: TangentOgiveNose
    body: CylindricalBody
    fins: TrapezoidalFinSet
    moment_reference_x: float
    boattail: ConicalBoattail | None = None

    def __post_init__(self) -> None:
        if not math.isclose(
            self.nose.base_diameter,
            self.body.diameter,
            rel_tol=1.0e-9,
            abs_tol=1.0e-12,
        ):
            raise ValueError("nose and cylindrical-body diameters must match")
        if self.boattail and self.boattail.exit_diameter >= self.body.diameter:
            raise ValueError("boattail exit diameter must be less than body diameter")
        if not 0.0 <= self.moment_reference_x <= self.length:
            raise ValueError("moment reference must lie within the rocket length")
        if self.fins.root_leading_edge_x < self.nose.length:
            raise ValueError("the baseline fin set must be mounted on the body")
        if self.fins.root_leading_edge_x + self.fins.root_chord > self.length:
            raise ValueError("fin root chord must lie within the rocket length")

    @property
    def reference_diameter(self) -> float:
        """Return maximum body diameter."""

        return self.body.diameter

    @property
    def reference_area(self) -> float:
        """Return ``pi*Dref^2/4`` using guide equation REF-01."""

        return math.pi * self.reference_diameter**2 / 4.0

    @property
    def length(self) -> float:
        """Return total nose-to-base length."""

        return self.nose.length + self.body.length + (
            self.boattail.length if self.boattail else 0.0
        )

    @property
    def base_diameter(self) -> float:
        """Return the aft base diameter."""

        return self.boattail.exit_diameter if self.boattail else self.body.diameter
