# Gravity model

`Gravity` is a native C++ `ForceTorque` model for Newtonian central-body
gravity. Add it to a vehicle in the same list as aerodynamics, engines, or
other force models:

```python
gravity = rsim.Gravity.earth()
vehicle = rsim.Vehicle(
    "rocket",
    j2000,
    dry_mass,
    dry_inertia,
    force_torques=[gravity],
)
```

No Python function is called after setup. `Vehicle` keeps shared ownership of
the bound C++ object, and `ForceTorqueRegistry` calls `Gravity::update()` while
the simulator runs with the Python GIL released.

## Position and frame convention

The vehicle's integration frame must be inertial. Its position is the body
origin expressed in that frame, while its center of gravity is stored in body
coordinates. Gravity first calculates the center-of-gravity position

```text
r_cg = r_body_origin + R_integration_from_body r_cg_body.
```

The configured central-body center is a fixed point expressed in the same
integration frame. The relative position and acceleration are

```text
r = r_cg - r_center
a = -mu r / |r|^3,
```

where `mu` has units of m³/s². The model multiplies acceleration by the current
vehicle mass, rotates the resulting integration-frame force into body
coordinates, and publishes it through `ForceTorque`. The application point is
the current CG, so gravity creates no artificial moment about the CG.

`Gravity.earth()` uses `mu = 3.986004418e14 m³/s²`. It assumes the central-body
center is the integration-frame origin unless `center_in_integration_frame` is
provided. A geocentric J2000 position must therefore include Earth's radius;
an altitude alone is not a valid position for this model.

```python
gravity = rsim.Gravity.earth(
    name="earth-gravity",
    center_in_integration_frame=np.zeros(3),
)
```

For another spherical central body, use the direct constructor:

```python
moon_gravity = rsim.Gravity(
    "moon-gravity",
    gravitational_parameter=4.9048695e12,
    center_in_integration_frame=np.zeros(3),
)
```

The point-mass field excludes nonspherical harmonics, third bodies, tides, and
relativistic corrections. It is appropriate for sounding-rocket trajectories
when those small perturbations are outside the required fidelity. Construction
rejects nonpositive or nonfinite `mu`; updating rejects a missing vehicle
attachment, a noninertial integration frame, and zero distance from the center.
