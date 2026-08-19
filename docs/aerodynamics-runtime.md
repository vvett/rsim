# Runtime aerodynamics

The compiled runtime uses the maintained `matio` library to read MATLAB files.
`AerodynamicDatabase` still validates the specific schema written by
`AeroCoefficientTable.save_mat()` rather than accepting arbitrary MATLAB
variables. The file must use format version 2 and contain the four axes, six
body-axis load coefficient arrays, and `xCP` produced by `rsim.aero`.

## Database evaluation

`AerodynamicDatabase` interpolates in this order:

```text
Mach, alpha_deg, phi_aero_deg, Reynolds
```

Mach, alpha, and Reynolds must lie inside their axes. `phi_aero` is wrapped by
the fin-set period stored in the metadata. Scalar magnitudes repeat directly,
but transverse Y/Z force and moment components are rotated into the requested
symmetry sector; wrapping those vector components without rotating them would
reverse the load direction in some quadrants. `RegularGridInterpolator<Dimensions>`
provides the common column-major interpolation interface. It uses dedicated
linear and bilinear implementations for one and two dimensions, then dispatches
three and higher dimensions to the header-only `mlinterp` library. The ordering
matches MATLAB/NumPy column-major table storage without a runtime data transpose.

`MatFileReader` similarly presents typed `readDoubleArray()`, `readUint32()`,
and `readUtf8()` methods while keeping matio handles and cleanup out of the
aerodynamics code. Both reusable wrappers live in `src/utility`.

```python
database = rsim.AerodynamicDatabase("rocket_aero.mat")
coefficients = database.evaluate(1.5, 3.0, 45.0, 2.0e6)
```

The required values are `CA`, `CY`, `CZ`, `Cmx`, `Cmy`, `Cmz`, and `xCP`.

## Dimensional loads

`Aerodynamics` converts coefficients using database reference area `S` and
reference length `L`:

```text
F_body = q S [-CA, -CY, -CZ]
```

The generator stores axial locations from the nose toward the tail, whereas
the vehicle geometry uses X from the engine toward the nose. Runtime CP is
therefore `x_cp_body = L - x_cp_database`. The transverse force is applied at
that transformed CP, and `ForceTorque` shifts its moment to the current CG.
This avoids treating the generator's legacy transverse-moment signs as native
right-hand-rule body moments. `Cmx` remains an intrinsic roll moment.

The registered stability diagnostic is

```text
staticStabilityMargin = (x_cg_body - x_cp_body) / reference_diameter
```

It is positive when CP is aft of CG in the vehicle coordinate convention.
`centerPressureX` and `staticStabilityMargin` are available to `DataLogger`.

Conditions are currently explicit rather than inferred from atmosphere and
wind models:

```python
conditions = rsim.AerodynamicConditions()
conditions.mach = 1.5
conditions.alpha_deg = 3.0
conditions.phi_aero_deg = 45.0
conditions.reynolds = 2.0e6
conditions.dynamic_pressure = 20_000.0

aero = rsim.Aerodynamics(
    "aerodynamics",
    rsim.AerodynamicDatabase("rocket_aero.mat"),
)
aero.set_conditions(conditions)
vehicle = rsim.Vehicle(
    "rocket", j2000, mass, inertia,
    force_torques=[aero],
)
```

## Simulator connection

`Aerodynamics` is a `ForceTorque` model registered through the vehicle's force
registry exactly like propulsion, gravity, or any future native C++ force model.
The runtime update order is:

```text
RigidBody -> ForceTorqueRegistry (Aerodynamics, ...) -> Kinematics
```

Each force model stores only its own result. The registry updates models in the
order they were added and sums their loads without names, flags, or type checks.

The current runtime does not yet calculate Mach, Reynolds, dynamic pressure,
alpha, or `phi_aero` from atmosphere, wind, and kinematic state. A later flight
conditions model can populate the same `AerodynamicConditions` interface
without changing database interpolation or load dimensionalization.

## Portable dependencies

Eigen, matio, and mlinterp are pinned and downloaded by CMake `FetchContent`.
An end user clones only the rsim repository; configuration retrieves and builds
the dependency sources automatically. Matio is built statically with its test
data submodule, HDF5, and zlib disabled because rsim generates uncompressed
MAT-v5 files. Mlinterp is header-only and its upstream test dependency is also
disabled.
