# Simplified J2000 and ECEF frames

`J2000EarthRotation` defines a rotating ECEF frame beneath fixed J2000 axes:

```text
J2000
└── ECEF   J2000EarthRotation
```

This is intentionally an engineering approximation. It uses one Julian date and
one rotation about the J2000 Z axis. It does not model:

- precession or nutation;
- polar motion;
- UT1 minus UTC corrections;
- TT, TAI, or relativistic time scales;
- variations in Earth's rotation rate.

## Input time

The only input is a two-part Julian date:

```cpp
rsim::JulianDate current_date{2451545.0, elapsed_days};
```

Splitting the date is optional from the model's perspective, but retaining a
small value in `second_part` avoids losing precision as simulation time advances.
The date may be based on UTC when sub-second Earth-orientation accuracy is not
required.

J2000 is Julian date `2451545.0`, corresponding to 2000 January 1 at 12:00.

## Rotation model

The model calculates approximate Greenwich mean sidereal angle in degrees:

```text
angle = 280.46061837
      + 360.98564736629 (julian_date - 2451545.0)
```

The angle is reduced to one revolution and used to define `j2000_from_ecef`.
The ECEF origin remains at the J2000 origin.

The frame motion uses the corresponding constant angular rate:

```text
360.98564736629 degrees per day
```

Angular acceleration and relative origin velocity and acceleration are zero.
Consequently, state conversion still includes the rotational velocity,
centripetal acceleration, and Coriolis terms handled by `FrameMotion`.

## Usage

```cpp
#include "EarthFrames.h"

#include <memory>

rsim::JulianDate current_date{2451545.0, 0.0};

auto j2000 = std::make_shared<rsim::Frame>(
    "J2000",
    rsim::InertialStatus::Inertial);
auto ecef = std::make_shared<rsim::Frame>(
    "ECEF",
    rsim::InertialStatus::NonInertial,
    j2000,
    std::make_unique<rsim::J2000EarthRotation>(current_date));

rsim::FrameGraph frames;
frames.addFrame(j2000);
frames.addFrame(ecef);
frames.update();
```

Advance the referenced date before updating the graph:

```cpp
current_date.second_part += timestep_seconds / 86400.0;
frames.update();
```

Convert an Earth-fixed object state into J2000:

```cpp
rsim::MotionState vehicle_in_j2000 =
    frames.convertState(*j2000, *ecef, vehicle_in_ecef);
```

A point stationary in ECEF correctly receives nonzero velocity and centripetal
acceleration in J2000.

## Accuracy boundary

This model is appropriate when a simple rotating Earth is more useful than
arcsecond-level orientation accuracy. It should not be mixed with software that
expects GCRS/ITRS, TEME, observed polar motion, or high-precision tracking data.
Those require an explicit standards-based conversion.
