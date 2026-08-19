# Aerodynamic database generation

`rsim.aero` is a pure-Python subpackage kept separate from the compiled
simulation extension. It generates fixed-geometry regular-grid aerodynamic
tables using the reproducible full-envelope route in the supplied Sounding
Rocket Aerodynamic Database Reference.

The baseline supports a tangent-ogive nose, cylindrical body, optional conical
boattail, and one fixed trapezoidal aft fin set. Dimensions and axial locations
use SI units. The independent variables are Mach, total angle of attack,
aerodynamic roll angle `phi_aero`, and body Reynolds number based on total
rocket length. Angles are in degrees, with `phi_aero = atan2(w, v)` for the
body-axis relative wind `[u, v, w]`.

```python
from rsim.aero import (
    AeroGrid,
    AerodynamicOptions,
    CylindricalBody,
    RocketGeometry,
    TangentOgiveNose,
    TrapezoidalFinSet,
    generate_aero_table,
)

geometry = RocketGeometry(
    nose=TangentOgiveNose(length=0.3, base_diameter=0.1),
    body=CylindricalBody(length=3.7, diameter=0.1),
    fins=TrapezoidalFinSet(
        count=4,
        root_chord=0.4,
        tip_chord=0.2,
        semi_span=0.15,
        tip_leading_edge_offset=0.12,
        thickness=0.008,
        root_leading_edge_x=3.4,
    ),
    moment_reference_x=2.4,
)

grid = AeroGrid.from_spacing(
    mach=(0.5, 5.0, 0.5),
    alpha_deg=(-10.0, 10.0, 1.0),
    phi_aero_deg=(-180.0, 180.0, 5.0),
    reynolds=(1.0e6, 5.0e6, 1.0e6),
)

table = generate_aero_table(
    geometry,
    grid,
    AerodynamicOptions(boundary_layer="turbulent"),
)
table.save_npz("rocket_aero.npz")
table.save_mat("rocket_aero.mat")
```

All arrays have dimension order
`(mach, alpha_deg, phi_aero_deg, reynolds)` and the mesh
arrays use `numpy.meshgrid(..., indexing="ij")`, equivalent to MATLAB
`ndgrid`. Final values are in `table.coefficients`; component contributions,
method/status flags, axes, meshes, and audit metadata are stored separately.
The NPZ file contains only named NumPy arrays plus JSON metadata—no Python
pickles—which keeps a future C++ reader straightforward.

The MAT file is an uncompressed MATLAB level-5 file and does not require
SciPy. Axes use names such as `axis__phi_aero_deg`, while coefficient arrays
use names such as `coefficient__CA`. Method flags are uint32 arrays whose
string-to-code dictionaries are stored in the UTF-8 JSON byte array
`metadata_json_utf8`.

```matlab
load rocket_aero.mat
CA = interpn(axis__mach, axis__alpha_deg, axis__phi_aero_deg, ...
             axis__reynolds, coefficient__CA, M, alpha, phi_aero, Re);
```

## Aerodynamic roll-angle synthesis

The axisymmetric body does not change coefficient magnitude with `phi_aero`;
its transverse force and moment are simply resolved into body Y/Z axes. The
discrete fin set is evaluated panel by panel. For fin-plane azimuth `phi_i`,
the implementation note's leading low-angle equivalent incidence is

```text
alpha_eq_i = alpha * sin(phi_aero - phi_i)
```

Each fin uses the existing Mach-dependent normal-force and CP methods, after
which its load and moment are rotated into body axes. A `2/N` normalization
preserves the original small-angle aggregated slope for an evenly spaced
`N`-fin set. Induced drag is scaled by the sum of individual squared fin
loads, producing discrete roll-angle dependence. Outputs include `CY`, `CZ`,
and `Cmx/Cmy/Cmz`; `CN`, `CS`, and `Cm` are the corresponding along-incidence
and cross-incidence summaries.

The note recommends runtime individual-fin synthesis when possible. This
generator applies the same synthesis while creating the grid because a full
`phi_aero` lookup dimension was requested. It is a low-angle engineering
approximation, not a direct complete-vehicle CFD or wind-tunnel roll sweep.

## Implemented routing

- Below Mach 1.2, body normal force uses slender-body/Jorgensen equations.
- Tangent-ogive transonic wave drag uses the digitized NWL TR-2796 Figure 3.
- From Mach 1.2 to the configured high-Mach switch, body pressure uses the
  guide's local shock-expansion panel fallback.
- At Mach 5 or higher, body pressure uses Modified Newtonian panels.
- Fin normal force routes through USAF DATCOM FDAT-03, the documented
  transonic smoothstep fallback, and the FSUP-02/Ackeret supersonic fallback.
- Digitized AIAA methods supply body-to-fin upwash, nonlinear fin force, and
  fin center-of-pressure corrections.
- Skin friction, base/boattail drag, fin wave drag, induced drag, final
  coefficient assembly, static slopes, pitch damping fallback, and roll
  damping are stored independently for auditing.

Every equation-level function documents its assumptions, source, and useful
range in its docstring. `method_flags` distinguishes exact equations, chart
interpolation, assembly identities, and transparent fallbacks.

## Limitations

This is an engineering database generator, not an exact Missile DATCOM clone.
It is limited to conventional symmetric rockets, approximately
`|alpha| <= 10 deg`, and `-180 <= phi_aero <= 180 deg`. It does not model
motor-on base pressure, protuberances,
control deflection, multiple fin sets, plume interaction, aeroelasticity,
unsteady transonic flow, real-gas effects, or configuration-specific
transition. Validate transonic and mission-critical results against wind-tunnel
data, CFD, Missile DATCOM, or another trusted source before flight use.
