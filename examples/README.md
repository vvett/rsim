# Liquid sounding-rocket examples

## Vespula-inspired LOX/RP-1, shared COPV, flight, and recovery

Run the longer native example from the repository root:

```bash
source .venv/bin/activate
python examples/vespula_inspired_lox_rp1.py
```

The native-runtime demonstration intentionally stops when the main parachute
is deployed. It validates ascent, burnout, apogee/drogue logic, and the
descending 2000-ft main event; it is not a simulation through touchdown.

This example is inspired by, but is not a reproduction of, the Yellow Jacket
Space Program (YJSP) Vespula. Public YJSP material describes Vespula as a
23-foot (about 7.0 m), 10-inch (0.254 m) diameter, 2,500 lbf (about 11.1 kN)
LOX/kerosene vehicle with an actively controlled bang-bang pressure-fed system.
YJSP also says that gaseous nitrogen pressurizes both propellant tanks. Its
Prometheus engine page gives a 2.0 oxidizer/fuel mixture ratio and reports a
near-full-duration 22-second LOX/RP-1 hotfire. Sources:

- [YJSP KeroLOX / Vespula overview](https://www.gtspaceprogram.com/kerolox)
- [YJSP Prometheus engine hotfire](https://www.gtspaceprogram.com/press-release/2025/10/28/ablative-engine-hotfire)

The following values are authoritative **mockup inputs supplied for this
example**, not additional public claims about the flown vehicle:

| Mockup input | Customary value | SI value used internally |
|---|---:|---:|
| Steady oxidizer/fuel ratio | 2.0 | 2.0 |
| LOX mass flow | 7.055 lb/s | 3.2 kg/s |
| RP-1 mass flow | 3.527 lb/s | 1.6 kg/s |
| LOX tank internal volume | 55 L | 0.055 m^3 |
| Initial liquid fill / ullage | 95% / 5% | 0.95 / 0.05 |
| Chamber pressure | 300 psi | 2.068 MPa |
| Injector stiffness / drop | 18% / 54 psi | 0.18 / 0.372 MPa |
| Shared nitrogen COPV | 20 L at 5.4 ksi | 0.020 m^3 at 37.23 MPa |

The fuel tank is sized for an assumed 18.0 s design burn: `18.0 s * 1.6
kg/s = 28.8 kg` of RP-1, and `28.8 / (810 kg/m^3 * 0.95) = 0.03743 m^3`
(37.43 L) internal volume. With the chosen 0.114 m inner radius this gives a
0.9167 m internal and 0.9287 m outer cylindrical length. The simulated powered
segment is 16.0 s, leaving margin before nominal depletion.

| Engineering assumption | Example value | Rationale |
|---|---:|---|
| LOX / RP-1 density | 1141 / 810 kg/m^3 | Representative sizing values |
| Tank temperature | 90 / 293.15 K | Single-temperature tank-model inputs |
| Engine expansion ratio | 6.0 | Not a public Vespula value |
| Tank operating pressure | 462 psi / 3.185 MPa | `Pc + 3 * 54 psi`, per compound-edge convention |
| LOX / RP-1 injector CdA | 1.0978e-4 / 6.5149e-5 m^2 | `mdot/sqrt(2*rho*DeltaP)` at 54 psi |
| LOX / RP-1 line diameter | 13.56 / 11.66 mm | Public `Line.mass_flow` sizing at a separate 54 psi loss |
| Valve opening | 95% | Valve CdA compensates for opening at a separate 54 psi loss |
| Tank gas closure | Isothermal (`n=1`) | N2 mass is replenished from the shared COPV |
| LOX / RP-1 pressure setpoint | 462 psi / 3.185 MPa | `300 + 54 + 54 + 54 psi` at both target flows |
| Controller thresholds | 457 / 467 psi | Open below 3.151 MPa; close above 3.220 MPa |
| Pressurant valve capacity | 4x initial makeup flow | Provides authority for hysteretic cycling as COPV pressure decays |
| Native integration / sampling | 0.01 / 0.10 s | Ten native steps per held-condition segment |

The public `CompoundEdge` currently divides its endpoint pressure difference
equally among its parts. The example therefore assigns 54 psi independently to
the line, valve, and injector, for 162 psi total tank-to-chamber drop. This is
the implemented numerical convention, not an assertion about Vespula plumbing.
Both liquid legs use public `Line + Valve + Orifice` compound edges.

One native 20 L N2 `Tank` at 5.4 ksi branches through separately sized gas
valves into both 5%-ullage propellant tanks. Two independent native
`BangBangPressureController` models command those valves before `Propulsion`
runs each step. Each opens below its 457 psi lower threshold, closes above its
467 psi upper threshold, and retains its preceding command inside the deadband.
Native transport subtracts the sum of branch flows from the COPV and deposits
each branch's mass and enthalpy into its receiver; no Python pressure-control or
pressurant bookkeeping runs during simulator stepping.

`data/lox_rp1_illustrative.csv` is a smooth illustrative software fixture, not
RocketCEA output and not design-quality data. Python loads it into an
`rsim.CombustionTable` and `.native` copies all arrays into C++ before stepping.
If RocketCEA is installed, regenerate a separate analysis table with:

```bash
python examples/generate_lox_rp1_table.py
```

Review RocketCEA inputs and generated values independently before analysis or
design. RocketCEA is optional and is not needed for the deterministic example.

### Vespula-like aerodynamic assumptions

Only overall length and diameter below are sourced YJSP facts. The remaining
external geometry is a documented best estimate for exercising rsim's public
aerodynamic database generator:

| Geometry/reference item | Value | Status |
|---|---:|---|
| Overall length / body diameter | 7.00 / 0.254 m | Rounded public dimensions |
| Nose | 0.80 m tangent ogive | Assumption |
| Cylindrical body | 6.20 m | Assumption completing overall length |
| Fins | 4 trapezoidal | Assumption |
| Root / tip chord | 0.75 / 0.32 m | Assumption |
| Exposed semi-span / thickness | 0.28 / 0.012 m | Assumption |
| Tip leading-edge offset | 0.25 m aft | Assumption |
| Root leading-edge location | 6.15 m aft of nose | Assumption |
| Reference area / length | 0.05067 m^2 / 7.00 m | Maximum-body area and overall length |
| Moment reference | 4.40 m aft of nose | Assumption |

Python generates the regular-grid aerodynamic table at setup, exports a
temporary MAT file, and constructs `rsim.AerodynamicDatabase`; its reader copies
the axes and coefficients into native C++ memory before stepping.
`Aerodynamics` is registered beside `Engine` and native Earth `Gravity`. The
native flight-state model supplies powered-ascent state tracking.

### Illustrative recovery setup

No public Vespula recovery-system dimensions are assumed here. The example
uses the requested terminal-speed targets only to exercise rsim's native
flight-state triggers and `Parachute.from_terminal_descent_speed` sizing API:

| Recovery assumption | Customary value | SI value used internally |
|---|---:|---:|
| Drogue target speed | 700 ft/s | 213.36 m/s |
| Drogue deployment | Apogee | Native downward crossing outside a 0.1 m/s deadband |
| Drogue drag coefficient | - | 0.8 |
| Main target speed | 5 ft/s | 1.524 m/s |
| Main deployment | 2,000 ft AGL descending | 609.6 m AGL descending crossing |
| Main drag coefficient | - | 1.5 |
| Recovery design mass | 275.6 lbm | 125.0 kg |
| Recovery design density | 0.0624 lbm/ft^3 | 1.0 kg/m^3 |

At those setup-time assumptions, the public terminal-speed factory gives a
0.0673 m^2 drogue area (0.293 m equivalent diameter) and a 703.7 m^2 main area
(29.93 m equivalent diameter). The requested 5 ft/s target makes the implied
main exceptionally large; these are transparent software-fixture inputs, not
a claim that such parachutes are practical or Vespula-like. Runtime parachute
area is fixed, and the example's constant density is another simplifying
assumption.

The checked run uses native state-dependent simulator timesteps: 5 ms on the
rail, 10 ms under power, 50 ms while coasting, and 1 ms for drogue/main
deployment. It physically advances powered ascent and coast through the native
apogee crossing and drogue deployment. The illustrative 5 ft/s main implies a
703.7 m² canopy and would make a full 1 ms descent needlessly expensive, so the
example then places two descending states immediately above and below 2,000 ft
AGL. Native `TruthState` crossing logic observes the main deployment, and the next
steps use the main-state timestep. This is explicitly a recovery-transition
integration check, not a continuous descent trajectory.

Until a native atmosphere/air-data coupler is available, Mach, Reynolds number,
and dynamic pressure are updated in Python at approximately 0.10 s sample
intervals using a constant sea-level atmosphere, then held during each native
`run_steps()` segment. There are no Python callbacks during stepping. Gravity
and geocentric altitude are active, but the trajectory is still illustrative:
density does not decay with altitude and the reported rail state does not
impose a mechanical rail constraint.

## Short N2O/ethanol smoke test

Run the deterministic, half-second native simulation from the repository root:

```bash
source .venv/bin/activate
python examples/liquid_sounding_rocket.py
```

The example loads `data/n2o_ethanol_illustrative.csv` in Python and copies its
NumPy arrays into C++ during setup. The complete `Simulator.run_steps()` loop is
native and does not call back into Python. The checked-in data are smooth,
plausible illustrative values for exercising the software; they are not
RocketCEA results and must not be used for engine design.

If RocketCEA is installed, generate a separate analysis table with:

```bash
python examples/generate_n2o_ethanol_table.py
```

Review generated values independently before using them for analysis or design.
