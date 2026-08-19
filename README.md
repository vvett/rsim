# rsim

C++ six-degree-of-freedom simulation components with Python bindings.

## Python development install

Use Python 3.9 or newer with `pip` and a C++17 compiler:

```bash
python -m pip install -e .
python -m unittest discover -s tests -p 'test_python_*.py'
```

Pip's isolated build installs scikit-build-core, pybind11, CMake, Ninja,
NumPy, and SciPy automatically. Do not pass
`--no-build-isolation` unless the native build tools are already installed in
the active environment.

RocketCEA is only required when regenerating combustion-performance tables.
Install that optional tool with `python -m pip install -e '.[rocketcea]'`; on
systems without a compatible RocketCEA wheel, this also requires a Fortran
compiler such as `gfortran`. Normal simulation with saved tables does not need
RocketCEA.

CMake fetches pinned Eigen, matio, and mlinterp sources automatically; users
do not need to install or clone those C++ dependencies separately.

```python
from rsim import RigidBody

inertia = ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0))
body = RigidBody(10.0, inertia)
print(body.mass)
```

## Native C++ build

```bash
cmake -S . -B build -DRSIM_BUILD_PYTHON=OFF -DRSIM_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

`main.cc` can be built by adding `-DRSIM_BUILD_DEMO=ON`.

## Adding APIs

Public C++ APIs belong in headers under `src`. Add each new implementation file to
`rsim_core`, then bind its public classes, functions, and enums in
`python/bindings.cpp` (or a focused binding source when that file grows). Add both
C++ behavior tests and Python binding tests for every new public API.

## Frame system

See [docs/frames.md](docs/frames.md) for the frame API, transform convention,
relative-motion framework, state conversion, and update behavior. See
[docs/earth-frames.md](docs/earth-frames.md) for the simplified J2000/ECEF
conversion and its single Julian-date input.

## Simulation architecture

See [docs/simulation.md](docs/simulation.md) for model update rates, simulator
scheduling, vehicle ownership, Python lifetime handling, and RK4 placement.

## Aerodynamic database generation

See [docs/aero.md](docs/aero.md) for the pure-Python `rsim.aero` package, which
generates fixed-geometry Mach/angle-of-attack/Reynolds coefficient tables.
See [docs/aerodynamics-runtime.md](docs/aerodynamics-runtime.md) for loading
those tables in C++, dimensionalizing body-axis loads, and simulator wiring.

## Geometry and variable mass

See [docs/geometry.md](docs/geometry.md) for structural solid geometry, dry mass
aggregation, liquid tank assumptions, and time-varying propellant properties.

## Liquid propulsion

See [docs/propulsion.md](docs/propulsion.md) for native network topology,
pressure-driven restrictions, combustion-table setup, engine loads, vehicle
scheduling, and generic variable-mass coupling.
