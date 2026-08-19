# Native liquid propulsion

`Propulsion` owns a finalized graph of native nodes and flow components. Build
the graph at setup time, call `finalize()`, and pass the model and its engines to
the vehicle separately:

```python
propulsion = rsim.Propulsion()
propulsion.connect(oxidizer_tank, engine, rsim.Valve(1e-5),
                   upstream_phase=rsim.FluidPhase.LIQUID,
                   downstream_role=rsim.InletRole.OXIDIZER)
propulsion.connect(fuel_tank, engine, rsim.Valve(1e-5),
                   upstream_phase=rsim.FluidPhase.LIQUID,
                   downstream_role=rsim.InletRole.FUEL)
propulsion.finalize()
vehicle = rsim.Vehicle(..., models=[propulsion], force_torques=[engine])
vehicle.rigid_body.add_variable_mass_source(oxidizer_tank)
vehicle.rigid_body.add_variable_mass_source(fuel_tank)
```

Pressure-fed branches can be commanded entirely in native code. Schedule each
controller before the propulsion network it commands:

```python
pressurant_valve = rsim.Valve(2.0e-7, opening=0.0)
pressure_control = rsim.BangBangPressureController(
    "LOX pressure control",
    lox_tank,
    pressurant_valve,
    pressure_setpoint=3.2e6,
    deadband_half_width=5.0 * 6894.757293168,
    initially_open=False,
)
propulsion.connect(copv, lox_tank, pressurant_valve,
                   upstream_phase=rsim.FluidPhase.GAS)
vehicle = rsim.Vehicle(
    ...,
    models=[pressure_control, propulsion],
)
```

`BangBangPressureController` opens its valve fully below `setpoint -
deadband_half_width`, closes it fully above `setpoint + deadband_half_width`,
and retains the preceding command at and between those thresholds. Multiple
controllers may independently command separate valves fed by one shared COPV.
The controllers are ordinary native `Model` children: insertion order makes
their commands visible to the following `Propulsion.update()` in the same
simulation step, without a Python runtime callback.

The simulator order is vehicle, additional models in insertion order, rigid
body, force registry, and kinematics. Thus an `Engine.update()` publishes the
thrust cached by the immediately preceding propulsion solve. `step()` and
`run_steps()` keep the GIL released throughout this native schedule.

## Assumptions and units

All inputs are SI. Edge flow is positive from the declared upstream node to the
downstream node. Liquid fills a cylindrical tank from its aft end in body +X.
The current implementation conserves liquid and pressurant mass and their
stored internal-energy scalars. Ullage closure is ideal-gas blowdown with
adiabatic, isentropic, or `p V^n` polytropic volume scaling. Collapse factor is
actual divided by ideal pressurant mass, so only the effective ideal fraction
contributes pressure.

`BlowdownMode.ENERGY_CONSERVING` uses the stored pressurant internal energy and
the open-system balance `dU = -h_out dm`. For a well-mixed ideal gas this gives
the usual adiabatic blowdown cooling as gas mass leaves a rigid COPV, while
pressure is recomputed from `p = m R T / V` using the current ullage volume.

Tank-to-tank edges are conservative transports. For every signed edge, the
native update identifies the actual high-to-low-pressure source, removes
`abs(mdot)*dt` from that source phase, and adds exactly the same mass to a
compatible downstream tank phase. The source's transported specific enthalpy
is removed from its stored energy and the identical energy is deposited at the
receiver. Deltas are accumulated before any tank is mutated, so one COPV can
feed arbitrary unequal branches without order dependence and its depletion is
the sum of all outgoing branches. A zero-opening valve contributes exactly zero
to both sides, and a compatible negative edge reverses both transfers.

The phase selected by `connect(..., upstream_phase=...)` applies at both tank
ports. Tank-to-tank connections require identical phase species (for example,
Nitrogen gas to Nitrogen ullage or LOX liquid to LOX inventory); incompatible
species are rejected during setup. Liquid receiver capacity and source
inventory are checked before applying any network delta. Combustors and ambient
nodes remain terminal sinks rather than finite storage volumes.

When a requested tank outflow would exceed the remaining inventory, the native
update reduces every branch drawing that phase by the same fraction and empties
the tank exactly. If that branch feeds a bipropellant combustor, both propellant
inlets are limited by the same fraction so mixture ratio is preserved during
the final partial burn. The combustor then flames out when either oxidizer or
fuel is unavailable: chamber pressure is published at its one-pascal numerical
floor, inlet flows are zero, and an `Engine` produces zero thrust. A native
`TruthEvent` can use the registered `chamberPressure` state to disable the
surrounding `Propulsion` model or command other post-burn actions.

Liquid restrictions use the signed incompressible `CdA sqrt(2 rho |dp|)` law.
Gas restrictions use the standard perfect-gas subcritical/choked relation.
Lines use Darcy--Weisbach with `64/Re`, a turbulent explicit Colebrook
approximation, and a smooth transition. A compound edge stores its intermediate
pressures and exposes one common limiting mass flow. Its current closure divides
the endpoint pressure difference into equal pressure intervals, one per
component, then selects the limiting component flow. Size each component with
that explicit convention in mind; this is not a full series pressure-allocation
solve. Reverse tank flow is accepted only when phase and species are compatible.

Engine chamber closure is `mdot = Pc At / cstar`. The cached thrust is
`(Cf_vac - epsilon Pambient/Pc) Pc At`; force direction is normalized and its
moment comes from the ordinary `ForceTorque` application point. Combustion
tables are generated or loaded in Python and immediately copied into native
vectors in MATLAB/Fortran order. Runtime never calls Python, SciPy, or RocketCEA.

`ThermalContact`, `Pump`, and `Turbine` are configuration skeletons. Nonzero
thermal contacts warn once during finalization and do not change energy yet;
finalizing a graph containing a pump or turbine fails clearly.
