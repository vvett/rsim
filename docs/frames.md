# Frame transformations

`rsim` represents Cartesian reference frames as a tree. Each child stores one
rigid transform into its parent, and the library composes those transforms to
convert a position or vector between any two frames with a common root.

The public declarations and concise API documentation are in `src/Frames.h`.
This guide develops the mathematical convention used by that API.

## Transform notation

The name `a_from_b` means “a transform that accepts coordinates expressed in B
and returns coordinates expressed in A.” The order is deliberately explicit:

```text
p_a = a_from_b.applyPosition(p_b)
```

A rigid transform contains a rotation matrix `R_ab` and translation `t_ab`:

```text
p_a = R_ab p_b + t_ab
```

`t_ab` is the position of B's origin expressed in A. It is not simply an
arbitrary displacement: its meaning depends on both the source and destination
frames.

The equivalent homogeneous-coordinate representation is:

```text
        [ R_ab  t_ab ]
T_ab =  [             ]
        [  0      1   ]
```

`Transform` stores the rotation as an Eigen quaternion instead of a 3-by-3
matrix. Quaternions are compact, compose efficiently, avoid the singularity of
Euler angles, and can be renormalized to limit numerical drift.

## Quaternion rotation

A quaternion is written as:

```text
q = w + xi + yj + zk
```

For a unit quaternion `q`, a vector `v` can be treated as the pure quaternion
`(0, v)`. Rotation is mathematically:

```text
v_rotated = q (0, v) q*
```

where `q*` is the quaternion conjugate. Eigen implements this operation through:

```cpp
Eigen::Vector3d rotated = quaternion * vector;
```

The `Transform` constructor normalizes its quaternion. A zero quaternion cannot
represent a rotation and is rejected.

## Positions and free vectors

A position identifies a point, so both rotation and translation apply:

```text
p_a = R_ab p_b + t_ab
```

This is implemented by `applyPosition()`.

A free vector represents a direction or difference between two positions. For
example, if `v_b = p2_b - p1_b`, then:

```text
v_a = (R_ab p2_b + t_ab) - (R_ab p1_b + t_ab)
    = R_ab (p2_b - p1_b)
    = R_ab v_b
```

The translations cancel. Consequently `applyVector()` performs only rotation.
Directions, displacements, forces, and angular velocities are normally free
vectors. Absolute locations must use `applyPosition()`.

## Composition

Suppose B is a child of A and C is a child of B:

```text
p_a = R_ab p_b + t_ab
p_b = R_bc p_c + t_bc
```

Substituting the second equation into the first gives:

```text
p_a = R_ab (R_bc p_c + t_bc) + t_ab
    = (R_ab R_bc) p_c + (R_ab t_bc + t_ab)
```

Therefore:

```text
R_ac = R_ab R_bc
t_ac = R_ab t_bc + t_ab
```

The API expresses this as:

```cpp
rsim::Transform a_from_c = a_from_b * b_from_c;
```

The right-hand transform is applied first, just as with ordinary matrix
multiplication. The implementation uses quaternion multiplication for the
rotation and `a_from_b.applyPosition(t_bc)` for the translation.

## Inversion

Starting with:

```text
p_a = R_ab p_b + t_ab
```

subtract the translation and premultiply by the inverse rotation:

```text
p_b = R_ab^-1 (p_a - t_ab)
    = R_ab^T p_a - R_ab^T t_ab
```

A rotation matrix is orthogonal, so its inverse equals its transpose. Thus:

```text
R_ba = R_ab^T
t_ba = -R_ab^T t_ab
```

For a unit quaternion, the conjugate is its inverse. `Transform::inverse()` uses
the quaternion conjugate and rotates the negated translation accordingly.

## Frame nodes and providers

A `Frame` is a named node. Every child owns a `TransformProvider` that supplies
`parent_from_child`. The provider makes the transform model replaceable without
putting geodesy, astronomy, or vehicle dynamics into `Frame` itself.

Examples include:

```text
ECI root
└── ECEF       EciEcefProvider, updated from time and Earth orientation
    └── NED    fixed at a particular geodetic origin
        └── Body   updated from vehicle position and attitude
```

The provider owns its cached transform. Its `update()` method reads dependencies
provided at construction—such as a simulation clock or vehicle state—and
refreshes that cache. `parentFromChild()` returns the latest cached value.

`FixedTransformProvider::update()` intentionally does nothing. It satisfies the
same interface for transformations such as a fixed sensor mounting offset.

## Graph updates

`FrameGraph::update()` recursively visits parents before children. It maintains
two sets:

- `updating` contains nodes on the active recursion path and detects cycles.
- `updated` prevents a shared ancestor from being updated more than once in a
  single pass.

Parent-first ordering matters when a child calculation depends on an ancestor's
new state. An ECEF-to-body calculation, for example, should not combine a new
body attitude with an ECEF transform from the previous timestep.

## Transforming between arbitrary frames

To find `destination_from_source`, the graph independently walks from each frame
to its root and composes the edges:

```text
root_from_destination = T_rd
root_from_source      = T_rs
```

Both paths must end at the same root. The desired transform is:

```text
destination_from_source = inverse(T_rd) * T_rs
```

This follows because:

```text
p_root        = T_rs p_source
p_destination = inverse(T_rd) p_root
```

Combining them yields:

```text
p_destination = inverse(T_rd) T_rs p_source
```

The graph currently composes paths on demand rather than caching root transforms.
That keeps invalidation simple while the frame system is small. If profiling
later shows this traversal to be significant, root transforms can be cached once
per update generation without changing the public convention.

## Example

```cpp
#include "Frames.h"

#include <Eigen/Geometry>
#include <memory>

auto eci = std::make_shared<rsim::Frame>("eci");

auto ecef = std::make_shared<rsim::Frame>(
    "ecef",
    eci,
    std::make_unique<rsim::FixedTransformProvider>(rsim::Transform{
        Eigen::Quaterniond::Identity(),
        Eigen::Vector3d::Zero()
    }));

rsim::FrameGraph graph;
graph.addFrame(eci);
graph.addFrame(ecef);
graph.update();

const rsim::Transform eci_from_ecef = graph.transform(*eci, *ecef);
const Eigen::Vector3d point_eci =
    eci_from_ecef.applyPosition(Eigen::Vector3d{1.0, 2.0, 3.0});
```

The identity provider in this small example should eventually be replaced by an
ECI/ECEF provider backed by an Earth-orientation model.

## Geodetic coordinates are not a frame edge

Latitude, longitude, and altitude are coordinates on and above an ellipsoid.
Their conversion to ECEF is nonlinear, so a single `Transform` cannot represent
it globally. Keep that operation in a separate GeographicLib adapter:

```text
geodetic position ── nonlinear conversion ──> ECEF position
                                              │
                                              └─ Cartesian FrameGraph
```

After obtaining an ECEF position, the ordinary rigid frame graph can transform
it into NED, body, sensor, or inertial coordinates. A local NED frame at a fixed
geodetic origin is a valid frame node because its relationship to ECEF is a
local rigid rotation and translation.
