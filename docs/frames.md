# Frame and motion framework

The frame system separates three related ideas:

| Type | Meaning |
|---|---|
| `Transform` | The position and orientation between two frames at one instant. |
| `FrameMotion` | A `Transform` plus the linear and angular motion between those frames. |
| `MotionState` | The position, velocity, and acceleration of an object expressed in one frame. |

This separation prevents an instantaneous rotation from being mistaken for all
the information required to convert velocity or acceleration between moving
frames.

## Naming convention

Names such as `a_from_b` state the direction explicitly. They accept coordinates
expressed in B and return coordinates expressed in A:

```text
p_a = R_ab p_b + t_ab
```

`t_ab` is the location of B's origin expressed in A. Composition follows:

```cpp
a_from_c = a_from_b * b_from_c;
```

The right-hand relationship is applied first.

## Transform

`Transform` converts positions and free vectors:

```cpp
Eigen::Vector3d position_a = a_from_b.applyPosition(position_b);
Eigen::Vector3d force_a = a_from_b.applyVector(force_b);
```

A position receives rotation and translation. A direction, force, displacement,
or other free vector receives only rotation. `inverse()` reverses the conversion.

## MotionState

`MotionState` belongs to an object, not a frame relationship:

```cpp
rsim::MotionState rocket_in_ecef{
    position,
    velocity,
    acceleration
};
```

All three members must be expressed in the same frame.

## FrameMotion

`FrameMotion` describes a child frame relative to its parent. It contains:

- `transform()`: the instantaneous `parent_from_child` pose;
- `childOriginVelocity()`: velocity of the child origin relative to the parent;
- `childOriginAcceleration()`: acceleration of the child origin;
- `childAngularVelocity()`: angular velocity of the child axes;
- `childAngularAcceleration()`: angular acceleration of the child axes.

The four rate vectors are expressed in the parent frame.

For child-frame position `p`, velocity `v`, and acceleration `a`, define:

```text
r = R p
u = R v
```

`convertState()` calculates:

```text
p_parent = r + t

v_parent = R v
         + origin_velocity
         + angular_velocity × r

a_parent = R a
         + origin_acceleration
         + angular_acceleration × r
         + angular_velocity × (angular_velocity × r)
         + 2 angular_velocity × u
```

The additional acceleration terms account for angular acceleration, centripetal
acceleration, and Coriolis acceleration. `FrameMotion::fixed()` creates a
relationship whose four rate vectors are zero.

`FrameMotion::inverse()` reverses both the pose and all frame rates. Composition
with `operator*` combines the pose and rates through multiple frame edges.

## FrameDefinition

A `FrameDefinition` answers one concrete question:

> How is this child frame positioned and moving relative to its parent now?

It exposes:

```cpp
void update();
const FrameMotion& motionIntoParent() const;
```

`update()` reads the current simulation inputs and recalculates the cached
relationship. `FixedFrameDefinition` is the simple implementation for sensor
mounts and other relationships that do not change.

An Earth rotation definition, vehicle body definition, and moving sensor
definition can all implement the same interface without putting their equations
inside `Frame`.

## Frame and FrameGraph

A `Frame` is a named node with an explicit `InertialStatus`, a parent, and a
`FrameDefinition`. A root has no parent or definition. A typical tree is:

```text
J2000 (inertial root)
└── ECEF (rotating Earth)
    └── launch NED
        └── vehicle body
            └── sensor
```

The status is specified when the frame is constructed:

```cpp
auto j2000 = std::make_shared<rsim::Frame>(
    "J2000",
    rsim::InertialStatus::Inertial);

auto ecef = std::make_shared<rsim::Frame>(
    "ECEF",
    rsim::InertialStatus::NonInertial,
    j2000,
    std::move(ecef_definition));
```

`isInertial()` is a convenient Boolean query. `inertialStatus()` returns the
enum. The classification is explicit metadata: it does not add or remove motion
terms. State conversion always uses the frame's actual `FrameMotion`.

`FrameGraph::update()` updates parents before children and updates each frame at
most once per pass. It also detects cycles.

The graph exposes three levels of conversion:

```cpp
// Position and orientation only.
Transform body_from_ecef = frames.transform(body, ecef);

// Pose plus relative linear/angular rates.
FrameMotion body_from_ecef_motion = frames.motion(body, ecef);

// Complete conversion of an object's position, velocity, and acceleration.
MotionState rocket_in_body =
    frames.convertState(body, ecef, rocket_in_ecef);
```

To convert between source and destination, the graph first composes each frame
to their common root:

```text
destination_from_source =
    inverse(root_from_destination) * root_from_source
```

Disconnected frame trees cannot be converted.

## Geodetic coordinates

Latitude, longitude, and altitude are not Cartesian coordinates and do not form
a rigid frame edge. A geodetic library must first convert a geodetic position to
an ECEF Cartesian position. The resulting Cartesian state can then move
through `FrameGraph`.

See [earth-frames.md](earth-frames.md) for the simplified J2000/ECEF definition.
