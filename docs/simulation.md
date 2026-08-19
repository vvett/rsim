# Model scheduling and vehicle ownership

`Model` is the small common interface for anything the simulator can schedule.
Each model has a name, an enabled flag, an exact `UpdateRate`, and the timing
context for its most recent call. A concrete model implements `update()` and can
read `simulationTime()` and `timeStep()` without receiving timing arguments
through several layers of calls.

## Exact update rates

`UpdateRate::timesPerStep(N)` calls a model `N` times during each global step.
Each call receives the timestep selected for that step divided by `N`, and
`substep()` runs from zero to `N - 1`. `UpdateRate::everyNSteps(N)` calls a
model after `N` global steps have accumulated and gives it the sum of the
actual elapsed timesteps. This remains correct when the selected timestep
changes between those steps.
`UpdateRate::everyStep()` is the common one-call-per-step case.

Only natural-number multipliers and their reciprocals are supported. Storing the
rate as integer counts avoids floating-point comparisons and ambiguous schedules.
Changing a model's rate or disabling it resets that model's accumulated schedule.
The reset is revision-based, so disabling and re-enabling a model between two
calls to `Simulator::step()` still discards its earlier partial interval.

`Simulator` stores non-owning model pointers and calls them in registration order.
The C++ caller must therefore keep registered models alive. The Python binding
enforces the same rule with `keep_alive`, so a registered Python object cannot be
garbage-collected while its simulator still exists.

## Vehicle ownership

`Vehicle` is the lifetime owner of the models that define one physical vehicle:

- An insertion-ordered shared `Model` list contains dependency models such as
  a native propulsion network.
- `RigidBody` stores mass, center of gravity, and moment of inertia.
- `ForceTorqueRegistry` updates force models in insertion order and sums them.
- `Kinematics` solves motion and owns the vehicle's body frame.

The accessors `rigidBody()`, `forces()`, and `kinematics()` return references to
those owned objects. Their Python equivalents are `vehicle.rigid_body`,
`vehicle.forces`, and `vehicle.kinematics`.

`ForceTorque` is the abstract base for one native C++ force model. Each subclass
publishes one force, intrinsic torque, and body-coordinate application point.
`ConstantForceTorque` is the concrete prescribed-load implementation used for
manual loads and tests. `Aerodynamics` is another `ForceTorque` subclass.

Force models can be supplied during Python construction without adding special
flags or simulator registration logic:

```python
aero = rsim.Aerodynamics("aero", aero_database)
engine = rsim.ConstantForceTorque("engine")
vehicle = rsim.Vehicle(
    "rocket", j2000, mass, inertia,
    models=[propulsion],
    force_torques=[aero, engine],
)
```

The vehicle retains shared C++ ownership of models passed from Python. The
registry stores their C++ pointers and calls their C++ virtual `update()` methods.
`ForceTorque` deliberately has no Python override trampoline, so Python is used
only for setup and calls such as `simulator.step()`; the model loop does not
re-enter Python while the simulation is running. The `step()` and `run_steps()`
bindings release Python's GIL for the duration of the native simulation loop.

## State-dependent global timestep

`globalTimeStep()` remains the simulator's positive finite base timestep, and
`setGlobalTimeStep()` changes that fallback for future steps. A vehicle can
request different timesteps for selected `VehicleState` values:

```python
simulator = rsim.Simulator(0.02)
simulator.set_vehicle_state_time_steps(vehicle, {
    rsim.VehicleState.ON_RAIL: 0.002,
    rsim.VehicleState.THRUSTING: 0.005,
    rsim.VehicleState.COASTING: 0.05,
    rsim.VehicleState.DROGUE_DEPLOYED: 0.02,
    rsim.VehicleState.MAIN_DEPLOYED: 0.10,
})
```

The Python dictionary is validated and copied into native C++ storage during
setup. Selecting timesteps during `step()` reads only native vehicle states and
does not call Python. An unmapped state requests the base timestep. If policies
are registered for multiple vehicles, the simulator uses the smallest request
for numerical stability.

Selection happens once at the beginning of a global step. A
`FlightStateModel` transition during that step therefore changes the timestep
on the following step, keeping one consistent interval for every model in the
current step. `lastTimeStep()` (`last_time_step` in Python) reports the timestep
selected for the most recently started step. Before the first step it reports
the configured base value.

Replace a policy by calling `set_vehicle_state_time_steps()` again. Passing an
empty dictionary or calling `clear_vehicle_state_time_steps()` removes the
policy and restores base-timestep behavior. Every configured timestep must be
finite and greater than zero.

`Simulator::addVehicle()` registers the vehicle, every additional model in
insertion order, rigid body, force registry, and kinematics in that fixed
dependency order. The registry itself calls all enabled force models before
Kinematics consumes their sum. Force models are not separately registered with
Simulator and therefore cannot be updated twice.

## Integration location

The classical RK4 calculation is the reusable `rungeKutta4` function in
`Integration.h`. It is separate from `Simulator`, because the simulator owns
scheduling but does not know how a model's state derivative is calculated.
The same header provides a table-driven additive Runge--Kutta utility and the
canonical ARS(2,2,2) IMEX tableau. Its diagonal-stage solver is supplied by the
calling model, keeping the integration utility independent of propulsion and of
any specific nonlinear-solver library.

`Kinematics::update()` advances one coupled state containing parent-frame
position and velocity, the body-to-parent quaternion, and body-frame angular
velocity. Every RK4 stage solves the rigid-body linear and angular acceleration
using that stage's angular velocity. It also uses that stage's orientation to
rotate body-frame linear acceleration into the parent frame. The quaternion
derivative uses body-frame angular velocity, and the final quaternion is
normalized before the new frame motion is published.

Forces and torques are treated as constant during one Kinematics update. A
future state-dependent force model can be evaluated inside the stage derivative
if forces must vary across RK4 stages.

The integration frame can be selected from a `FrameGraph` by its unique name:

```python
frames.add_frame(j2000)
motion = rsim.Kinematics("body", frames, "j2000", loads, rigid_body)
```

Position and translational velocity are integrated in that selected frame. The
body frame created by `Kinematics` is attached directly beneath it. Selecting a
frame marked `NON_INERTIAL` emits a C++ warning and a Python `RuntimeWarning`.
Integration is still allowed, but Coriolis, centrifugal, Euler, and moving-origin
terms are not yet added automatically, so an inertial frame should normally be
selected.
