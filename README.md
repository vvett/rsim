# rsim

C++ six-degree-of-freedom simulation components with Python bindings.

## Python development install

Use Python 3.9 or newer with `pip`, a C++17 compiler, and CMake installed:

```bash
python -m pip install -e .
python -m unittest discover -s tests -p 'test_python_*.py'
```

```python
from rsim import Body

inertia = ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0))
body = Body(10.0, inertia)
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
quaternion mathematics, composition and inversion rules, and update behavior.
