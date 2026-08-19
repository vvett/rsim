# Flight state and parachutes

`FlightStateModel` turns native vehicle position, velocity, and selected thrust
loads into a readable flight phase. `Parachute` is an ordinary `ForceTorque`,
so any number of drogue or main parachutes can be added without changing
`Vehicle` or `Simulator`.

```python
config = rsim.FlightStateConfiguration()
config.origin_in_integration_frame = np.array([0.0, 0.0, 0.0])
config.up_direction_in_integration_frame = np.array([1.0, 0.0, 0.0])
config.rail_direction_in_integration_frame = np.array([1.0, 0.0, 0.0])
config.rail_length = 6.0
config.thrust_on_threshold = 100.0
config.thrust_off_threshold = 50.0
config.apogee_velocity_deadband = 0.1
config.drogue_trigger.deploy_at_apogee = True
config.main_trigger.maximum_altitude = 609.6  # 2000 ft AGL
config.main_trigger.require_descending_altitude_crossing = True

flight_state = rsim.FlightStateModel(
    configuration=config,
    thrust_sources=[engine],
)
drogue = rsim.Parachute.from_terminal_descent_speed(
    "drogue", drogue_cd, 213.36, apogee_mass, design_density,
    rsim.VehicleState.DROGUE_DEPLOYED,
)
main = rsim.Parachute.from_terminal_descent_speed(
    "main", main_cd, 1.524, main_deployment_mass, design_density,
    rsim.VehicleState.MAIN_DEPLOYED,
)

vehicle = rsim.Vehicle(
    "rocket", inertial_frame, mass, inertia,
    models=[propulsion, flight_state],
    force_torques=[engine, drogue, main],
)
```

All objects are connected during `Vehicle` construction. Simulation steps only
call native C++ methods, and `Simulator.step()` and `run_steps()` remain free of
Python callbacks.

## Transition order

The state machine advances in one direction:

```text
OnRail -> Thrusting -> Coasting -> Descending
                                      |
                                      v
                             DrogueDeployed -> MainDeployed -> Landed
```

- `OnRail` ends when projected displacement along the configured rail direction
  reaches the rail length. The next state is `Thrusting` when selected native
  thrust sources exceed the on threshold, otherwise it is `Coasting`.
- `Thrusting` ends below the lower off threshold. Separate on/off thresholds
  prevent noise around cutoff from changing the phase repeatedly.
- `Coasting` becomes `Descending` when projected velocity along `up_direction`
  falls below its configured threshold.
- `deploy_at_apogee` requires confirmed upward motion followed by the first
  velocity sample below the negative apogee deadband. Zero and small noisy
  samples inside the deadband do not trigger deployment. The crossing latches,
  so the drogue cannot retract or redeploy if velocity later changes sign.
- A threshold deployment trigger requires its maximum altitude, maximum
  vertical velocity, and earliest simulation time conditions. Setting
  `require_descending_altitude_crossing` additionally requires the vehicle to
  cross `maximum_altitude` from above while descending; merely beginning below
  the threshold does not count. The crossing is latched until phase ordering
  permits deployment.
- With the drogue stage enabled, main deployment is considered only after the
  drogue/apogee phase. Disabling the drogue stage explicitly permits a direct
  `Descending -> MainDeployed` transition. Either stage may be disabled.
- A descending or deployed vehicle becomes `Landed` at or below ground altitude
  when it is moving downward or no faster upward than the configured landing
  threshold.

Only one phase transition occurs per model update, and automatic transitions
never move backward. `vehicle.state` remains writable for initialization,
testing, or a setup-time manual override.

## Coordinates and scheduling

Position, wind, the origin, `up_direction`, and `rail_direction` are expressed
in the vehicle's integration frame. No Earth axis, spherical altitude, or launch
site is assumed. Altitude is the projection of displacement from the configured
origin onto the unit up direction. Rail distance uses the corresponding rail
projection.

Additional vehicle models run before `ForceTorqueRegistry`, so a deployment
transition is visible to every parachute during the same simulator step.
Kinematics then consumes the accumulated forces. Thrust phase detection reads
the most recently published native `ForceTorque` outputs, which introduces at
most one global-step delay at ignition or cutoff.

## Parachute force

After its deployment state is reached, a parachute latches deployed and applies

```text
F = -0.5 rho Cd A |v_relative| v_relative
```

where `v_relative` is vehicle velocity minus wind velocity in the integration
frame. The result is rotated into body coordinates before being published, and
the normal force-registry moment-arm calculation uses
`application_point_in_body`. Different parachutes retain independent `Cd`, area,
deployment phase, and application point; the registry sums them normally.

`from_terminal_descent_speed()` calculates a fixed setup-time area by balancing
weight and steady quadratic drag:

```text
A = 2 m g / (rho Cd v_terminal^2),  g = 9.80665 m/s^2
```

The inputs are the selected design mass and density, not live runtime values.
The returned parachute exposes both `reference_area` and the circular
`equivalent_diameter`. Its area stays fixed during simulation. Actual descent
speed therefore changes as vehicle mass, atmospheric density, wind, and flight
conditions differ from the sizing point.

The first implementation uses fixed density and fixed integration-frame wind.
`set_atmosphere()` updates those native values without installing a Python
callback. A future native atmosphere/wind model can call the same kind of native
conditions interface.

The rail state is classification only. It does not yet constrain vehicle motion
or calculate rail reaction forces.
