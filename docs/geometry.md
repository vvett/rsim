# Rocket geometry and variable mass properties

The `rsim.geometry` module builds dry mass properties from homogeneous solids.
It is intentionally separate from the aerodynamic geometry because structural
dimensions, material density, and tank interiors serve a different calculation.

All dimensions are in meters, densities are in kg/m³, masses are in kilograms,
and inertia tensors are in kg·m². The body `+X` axis points from aft to forward.
Shape centers lie on the X axis, and inertia tensors are expressed in body axes
about each shape's own center of mass.

## Available geometry

- `Material` assigns a name and homogeneous bulk density.
- `SolidCylinder` represents solid cylindrical equipment, engines, bulkheads,
  or other components that can reasonably be reduced to an equivalent cylinder.
- `CylindricalShell` represents a body tube or other finite-thickness tube.
- `ConicalShell` represents a nose wall as an outer solid cone minus an
  equal-length inner cone. It has an open flat base and does not model an ogive,
  rounded tip, internal ribs, or a base closure.
- `CylindricalTank` combines a thick side wall and two flat solid end caps. Its
  interior is also used to define the liquid capacity.
- `LiquidSoundingRocketGeometry` collects components and tanks and combines
  their dry mass properties.

The shape equations are the standard exact centroid and inertia equations for
uniform cylinders, thick cylindrical shells, and right circular cones. Composite
properties use mass-weighted first moments followed by the parallel-axis theorem:

```text
I_about_vehicle_cg = I_about_component_cg
                   + m ((d dot d) identity - d d_transpose)
```

where `d` runs from the vehicle center of mass to the component center of mass.

## Building a variable-mass body

```python
import rsim

aluminum = rsim.Material("aluminum", 2700.0)
rocket = rsim.LiquidSoundingRocketGeometry()

rocket.add_component(
    rsim.CylindricalShell(0.15, 0.003, 2.0, 1.0, aluminum)
)
rocket.add_component(
    rsim.SolidCylinder(0.12, 0.25, 0.0, aluminum)  # Equivalent engine solid.
)
rocket.add_component(
    rsim.ConicalShell(0.15, 0.6, 0.003, 2.0, aluminum)
)
rocket.add_tank(
    rsim.CylindricalTank(
        outer_radius=0.14,
        length=0.8,
        wall_thickness=0.004,
        center_x=1.2,
        material=aluminum,
        propellant_density=1000.0,
        initial_fill_fraction=0.9,
        mass_flow_rate=2.0,
    )
)

body = rocket.create_rigid_body("rocket-mass")
```

For a `Vehicle`, configure its already-owned body instead:

```python
rocket.configure_rigid_body(vehicle.rigid_body)
```

## Propellant assumptions

Each C++ `CylindricalPropellantTank` is axisymmetric and aligned with body X.
The liquid is treated as a uniform cylinder with the tank's full interior radius.
It occupies the tank from its aft boundary forward, so consumption shortens the
filled cylinder and moves its center of mass aft. Slosh, ullage dynamics, surface
tilt under acceleration, plumbing inventory, and separate liquid/gas phases are
not modeled.

At every scheduled `RigidBody::update()`, remaining propellant is calculated from
absolute simulator time and mass-flow rate. The body then recomputes total mass,
center of mass, and inertia from dry properties and all remaining liquids. Using
absolute time avoids drift from repeatedly subtracting `flow_rate * dt`.

Change tank state through `body.set_propellant_mass(index, mass)` or
`body.set_tank_mass_flow_rate(index, rate)`. Rate changes start at the model's
current simulation time and are not applied retroactively.
