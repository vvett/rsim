"""Axisymmetric geometry and dry mass properties for liquid sounding rockets.

The body ``+X`` axis runs from aft to forward. Dimensions use meters, density
uses kg/m^3, mass uses kilograms, and inertia uses kg*m^2. Every inertia tensor
is expressed in body axes about the corresponding center of mass.
"""

from dataclasses import dataclass, field
from math import pi
from typing import Protocol

import numpy as np


def _positive(name: str, value: float) -> None:
    """Reject a nonpositive physical dimension or material property."""
    if value <= 0.0:
        raise ValueError(f"{name} must be positive")


@dataclass(frozen=True)
class Material:
    """A homogeneous material defined by its positive bulk density."""

    name: str
    density: float

    def __post_init__(self) -> None:
        _positive("material density", self.density)


@dataclass(frozen=True)
class MassProperties:
    """Mass, body-frame center of mass, and centroidal inertia tensor."""

    mass: float
    center_of_mass: np.ndarray
    moment_of_inertia: np.ndarray

    def __post_init__(self) -> None:
        if self.mass < 0.0:
            raise ValueError("mass cannot be negative")
        center = np.asarray(self.center_of_mass, dtype=float)
        inertia = np.asarray(self.moment_of_inertia, dtype=float)
        if center.shape != (3,) or inertia.shape != (3, 3):
            raise ValueError("center must be length 3 and inertia must be 3-by-3")
        object.__setattr__(self, "center_of_mass", center)
        object.__setattr__(self, "moment_of_inertia", inertia)

    @staticmethod
    def combine(items: list["MassProperties"]) -> "MassProperties":
        """Combine centroidal inertias using first moments and parallel-axis shifts."""
        total_mass = sum(item.mass for item in items)
        if total_mass <= 0.0:
            raise ValueError("at least one component must have positive mass")
        center = sum(
            (item.mass * item.center_of_mass for item in items),
            start=np.zeros(3),
        ) / total_mass
        inertia = np.zeros((3, 3))
        for item in items:
            offset = item.center_of_mass - center
            # Parallel-axis theorem: I_new = I_cm + m[(d.d)I - d*d^T].
            inertia += item.moment_of_inertia + item.mass * (
                np.dot(offset, offset) * np.eye(3) - np.outer(offset, offset)
            )
        return MassProperties(total_mass, center, inertia)


class SolidShape(Protocol):
    """Structural interface implemented by every dry solid shape."""

    def volume(self) -> float:
        """Return solid material volume in cubic meters."""

    def mass_properties(self) -> MassProperties:
        """Return the homogeneous shape's centroidal mass properties."""


@dataclass(frozen=True)
class SolidCylinder:
    """A homogeneous cylinder aligned with body X and centered at ``center_x``."""

    radius: float
    length: float
    center_x: float
    material: Material

    def __post_init__(self) -> None:
        _positive("cylinder radius", self.radius)
        _positive("cylinder length", self.length)

    def volume(self) -> float:
        """Return ``pi*r^2*L``; dimensions must be positive."""
        return pi * self.radius**2 * self.length

    def mass_properties(self) -> MassProperties:
        """Use the exact uniform solid-cylinder inertia about its center."""
        mass = self.material.density * self.volume()
        axial = 0.5 * mass * self.radius**2
        transverse = mass * (3.0 * self.radius**2 + self.length**2) / 12.0
        return MassProperties(
            mass,
            np.array([self.center_x, 0.0, 0.0]),
            np.diag([axial, transverse, transverse]),
        )


@dataclass(frozen=True)
class CylindricalShell:
    """A finite-thickness homogeneous tube aligned with body X."""

    outer_radius: float
    thickness: float
    length: float
    center_x: float
    material: Material

    def __post_init__(self) -> None:
        _positive("shell outer radius", self.outer_radius)
        _positive("shell thickness", self.thickness)
        _positive("shell length", self.length)
        if self.thickness >= self.outer_radius:
            raise ValueError("shell thickness must be smaller than outer radius")

    @property
    def inner_radius(self) -> float:
        """Return the cavity radius after subtracting wall thickness."""
        return self.outer_radius - self.thickness

    def volume(self) -> float:
        """Subtract the coaxial inner-cylinder volume from the outer volume."""
        return pi * (self.outer_radius**2 - self.inner_radius**2) * self.length

    def mass_properties(self) -> MassProperties:
        """Use the exact uniform thick-walled-cylinder inertia."""
        mass = self.material.density * self.volume()
        radial_sum = self.outer_radius**2 + self.inner_radius**2
        axial = 0.5 * mass * radial_sum
        transverse = mass * (3.0 * radial_sum + self.length**2) / 12.0
        return MassProperties(
            mass,
            np.array([self.center_x, 0.0, 0.0]),
            np.diag([axial, transverse, transverse]),
        )


@dataclass(frozen=True)
class ConicalShell:
    """A nose-cone wall modeled as equal-length outer and inner solid cones."""

    base_radius: float
    length: float
    thickness: float
    base_x: float
    material: Material

    def __post_init__(self) -> None:
        _positive("cone base radius", self.base_radius)
        _positive("cone length", self.length)
        _positive("cone thickness", self.thickness)
        if self.thickness >= self.base_radius:
            raise ValueError("cone thickness must be smaller than base radius")

    def volume(self) -> float:
        """Subtract a coaxial inner cone; this assumes an open, flat base."""
        inner_radius = self.base_radius - self.thickness
        return pi * self.length * (self.base_radius**2 - inner_radius**2) / 3.0

    def mass_properties(self) -> MassProperties:
        """Subtract exact solid-cone moments; valid for homogeneous conical walls."""
        outer_volume = pi * self.base_radius**2 * self.length / 3.0
        inner_radius = self.base_radius - self.thickness
        inner_volume = pi * inner_radius**2 * self.length / 3.0
        outer_mass = self.material.density * outer_volume
        inner_mass = self.material.density * inner_volume
        mass = outer_mass - inner_mass
        axial = 0.3 * (outer_mass * self.base_radius**2 -
                       inner_mass * inner_radius**2)
        transverse = 3.0 / 80.0 * (
            outer_mass * (4.0 * self.base_radius**2 + self.length**2) -
            inner_mass * (4.0 * inner_radius**2 + self.length**2)
        )
        # A uniform cone's center is one quarter of its length forward of its base.
        center = np.array([self.base_x + 0.25 * self.length, 0.0, 0.0])
        return MassProperties(mass, center, np.diag([axial, transverse, transverse]))


@dataclass(frozen=True)
class CylindricalTank:
    """A closed cylindrical tank with a side wall and two flat end caps."""

    outer_radius: float
    length: float
    wall_thickness: float
    center_x: float
    material: Material
    propellant_density: float
    initial_fill_fraction: float = 1.0
    mass_flow_rate: float = 0.0

    def __post_init__(self) -> None:
        _positive("tank radius", self.outer_radius)
        _positive("tank length", self.length)
        _positive("tank wall thickness", self.wall_thickness)
        _positive("propellant density", self.propellant_density)
        if 2.0 * self.wall_thickness >= self.length:
            raise ValueError("tank end caps leave no positive interior length")
        if self.wall_thickness >= self.outer_radius:
            raise ValueError("tank wall leaves no positive interior radius")
        if not 0.0 <= self.initial_fill_fraction <= 1.0:
            raise ValueError("initial fill fraction must be between zero and one")
        if self.mass_flow_rate < 0.0:
            raise ValueError("mass flow rate cannot be negative")

    @property
    def inner_radius(self) -> float:
        """Return the usable cylindrical propellant radius."""
        return self.outer_radius - self.wall_thickness

    @property
    def inner_length(self) -> float:
        """Return the usable length between flat end caps."""
        return self.length - 2.0 * self.wall_thickness

    @property
    def propellant_capacity(self) -> float:
        """Return full liquid capacity using the specified bulk density."""
        return self.propellant_density * pi * self.inner_radius**2 * self.inner_length

    def mass_properties(self) -> MassProperties:
        """Combine the exact side-wall and two solid end-cap inertias."""
        side = CylindricalShell(
            self.outer_radius, self.wall_thickness, self.inner_length,
            self.center_x, self.material,
        )
        cap_offset = 0.5 * (self.length - self.wall_thickness)
        aft_cap = SolidCylinder(
            self.outer_radius, self.wall_thickness,
            self.center_x - cap_offset, self.material,
        )
        forward_cap = SolidCylinder(
            self.outer_radius, self.wall_thickness,
            self.center_x + cap_offset, self.material,
        )
        return MassProperties.combine([
            side.mass_properties(),
            aft_cap.mass_properties(),
            forward_cap.mass_properties(),
        ])

    def cpp_definition(self) -> dict[str, float]:
        """Return arguments needed by the C++ variable-propellant tank model."""
        return {
            "radius": self.inner_radius,
            "length": self.inner_length,
            "aft_position": self.center_x - 0.5 * self.inner_length,
            "propellant_density": self.propellant_density,
            "initial_propellant_mass": (
                self.initial_fill_fraction * self.propellant_capacity
            ),
            "mass_flow_rate": self.mass_flow_rate,
        }


@dataclass
class LiquidSoundingRocketGeometry:
    """Collect dry solids and propellant tanks for one axisymmetric rocket."""

    components: list[SolidShape] = field(default_factory=list)
    tanks: list[CylindricalTank] = field(default_factory=list)

    def add_component(self, component: SolidShape) -> None:
        """Add a dry structural component such as an engine or body tube."""
        self.components.append(component)

    def add_tank(self, tank: CylindricalTank) -> None:
        """Add tank walls to dry mass and retain its liquid interior definition."""
        self.tanks.append(tank)

    def dry_mass_properties(self) -> MassProperties:
        """Combine every dry component and tank wall about the rocket dry CG."""
        items = [component.mass_properties() for component in self.components]
        items.extend(tank.mass_properties() for tank in self.tanks)
        return MassProperties.combine(items)

    def create_rigid_body(self, name: str = "rigid-body"):
        """Create the bound C++ RigidBody and install every propellant tank."""
        # The local import keeps geometry calculations usable without loading C++ first.
        from ._rsim import CylindricalPropellantTank, RigidBody

        dry = self.dry_mass_properties()
        body = RigidBody(
            dry.mass, dry.moment_of_inertia, dry.center_of_mass, name
        )
        self.configure_rigid_body(body)
        return body

    def configure_rigid_body(self, body) -> None:
        """Install this geometry on an existing RigidBody, including vehicle-owned bodies."""
        from ._rsim import CylindricalPropellantTank

        dry = self.dry_mass_properties()
        body.set_dry_mass_properties(
            dry.mass, dry.moment_of_inertia, dry.center_of_mass
        )
        for tank in self.tanks:
            body.add_propellant_tank(
                CylindricalPropellantTank(**tank.cpp_definition())
            )
