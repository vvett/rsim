# Native runtime control and telemetry

`Simulator.run()` executes the complete step loop in C++ while the Python GIL
is released. Python is used to construct models, copy tables into native
memory, load a compiled flight-computer plugin, and read logged arrays after
the run. It is never called between integration steps.

## State variables and logging

A `Model` publishes a stable member from its constructor:

```cpp
addStateVariable(altitude_, "truthAltitude", "m");
addComputedStateVariable(
    [this]() -> StateValue { return calculateLoad(); }, "calculatedLoad", "N");
```

The direct form retains a read-only reference, so the member must live as long
as the non-copyable model. The computed form is native C++ and is useful when
the published quantity is derived. Aliases are registered before the run.
`DataLogger` discovers every registered model variable, stores samples in
native vectors, and creates NumPy arrays only when `to_dict()` is called.

Every logger key is canonical: `machine_namespace.alias`. The machine
namespace is deterministically derived from the display name, with spaces and
punctuation converted to underscores. Namespace collisions fail during setup.
Events always accept canonical names and accept a bare alias only when it is
unique in the recursively discovered model graph. That graph includes vehicle
children and individual force models inside `ForceTorqueRegistry`.

Floating-point, signed-integer, and Boolean channels remain typed in native
storage and become `float64`, `int64`, and `bool` NumPy arrays. The logger
records after kinematics, truth-state evaluation, and truth events, including
the final row that caused a stop.
`channel_metadata` exposes the canonical name, units, and scalar dtype without
reading the sample arrays.

## Truth events

`TruthEvent` expressions are parsed and resolved during native initialization:

```python
event = rsim.TruthEvent(
    "altitude limit",
    "truth_state.truthAltitude greater_than 20.0 and "
    "truth_state.vehicleState not_equal_to LANDED",
    stop_simulation=True,
)
```

Supported comparisons are `greater_than`, `less_than`,
`greater_than_or_equal`, `less_than_or_equal`, `equal_to`, and
`not_equal_to`. `and` binds more tightly than `or`; parentheses override that
order. Built-in aliases are `simulationTime`, `timeStep`, and `globalStep`.
An event can stop the run, enable or disable a model, or deploy a parachute.
No expression uses Python `eval`.
New event actions should use `add_actuator_command`, which preserves the same
command/actual separation as flight software. Direct model-enable and
parachute helpers remain compatibility APIs.

## Sensors, flight software, and actuators

`ScalarSensor` reads one registered truth alias at its own update rate and
holds the latest value, validity, timestamp, and sequence number. A
`FlightComputer` sees only sensor buffers and actuator outputs, not the
vehicle or truth objects.

`SpecificForceSensor` is the ideal body-axis accelerometer. It sums enabled
forces classified as applied, excludes gravity-classified forces, and divides
by current vehicle mass. A freely falling gravity-only vehicle therefore reads
zero, while thrust and aerodynamic/contact loads remain measurable.

Flight software is a separately compiled shared library using the versioned
C ABI in `FlightComputerAbi.h`. Only fixed-width/POD values, opaque state, and
caller-owned arrays cross the boundary. Each plugin descriptor declares named
sensor inputs, actuator outputs, numeric configuration entries, and telemetry
channels. The host resolves every name during setup, copies configuration into
ABI-owned structures, validates descriptor/buffer sizes, and uses positional
buffers only after resolution. The repository includes the Vespula C++ plugin,
a C compilation smoke plugin, and a Rust `cdylib` reference. Plugins are
trusted native code and must match the rsim ABI version.
Plugin loading uses `LoadLibrary`/`GetProcAddress`/`FreeLibrary` on Windows and
the corresponding `dlopen` family on POSIX systems.

The Rust source is a portability reference. Its build is verified only when a
Rust toolchain is available; the standard rsim CMake build does not require
Rust.

`ParachuteActuator`, `ValveActuator`, and `ModelEnableActuator` implement one
ideal-actuator interface and retain commanded and actual values.
`TruthState` is the authoritative class and derives recovery phases from
actual parachute deployment by default. `FlightStateModel` remains a deprecated
compatibility alias. Threshold-driven automatic recovery is available only
when `legacy_automatic_recovery` is explicitly enabled.

## Native environment and update order

`AtmosphereTable` copies altitude, density, sound speed, viscosity, and optional
three-axis wind arrays during setup. `FlightConditions` interpolates them,
transforms relative air velocity into body axes, derives alpha and aerodynamic
azimuth, and updates Mach, Reynolds number, dynamic pressure, aerodynamic
models, and parachute atmosphere natively. Values outside the table altitude
range clamp to an endpoint and set `atmosphereClamped`.

The scheduler processes environment, sensors, flight control, actuators,
plant, mass properties, forces, kinematics, and truth in that order. Insertion
order is preserved inside each phase. Events and logging follow the model
phases. State-dependent timesteps selected at the start of a step take effect
without Python involvement.

Before the first logged row, an observation-only initialization pass updates
environment, truth, and sensors with zero elapsed time. It does not integrate
dynamics or execute the flight computer.
