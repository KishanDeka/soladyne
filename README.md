# SURYA Dynamo v3 — exact C++/Python port

An object-oriented C++17 translation of `surya_dynamo_v3.f`, wrapped with
pybind11 and presented through Python, a CLI, plots, animation, and Streamlit.
The numerical sequence and legacy text formats are retained; this is not the
earlier simplified dynamo approximation.

![Toroidal-field evolution from the supplied legacy initial state](docs/assets/dynamo-demo.gif)

## Preserved from Fortran

- fixed 257×257 field grid and 513×513 midpoint grid;
- original alpha, diffusivity, density, differential-rotation, viscosity,
  Lambda-effect, stream-function, and meridional-flow expressions;
- custom Fortran `erf` approximation and associated Legendre recurrence;
- every expanded ADI coefficient for rotation, poloidal field, and toroidal field;
- first/second half-step ordering and all Thomas tridiagonal solves;
- Lorentz source calculation before both rotation solves;
- eleven iterations of the potential-field upper boundary;
- bottom, surface, and polar boundary update order;
- `dt=0.00005`, `v0=-20`, `et0=3`, `et1=0.04`, `al0=30` defaults;
- legacy `init.dat`, `omg_lat.dat`, `diffrot.dat`, and `final.dat` layouts.

The inactive, commented magnetic-buoyancy and diagnostic blocks remain visible
in [`legacy/source/surya_dynamo_v3.f`](legacy/source/surya_dynamo_v3.f), but are
not executed—matching the supplied source.

## Structure

```text
cpp/include/solar_dynamo/model.hpp   OOP public API and owned solver state
cpp/src/model.cpp                    exact profiles, coefficients, ADI updates, I/O
cpp/src/bindings.cpp                 pybind11 NumPy interface
python/solar_dynamo/                 orchestration, CLI, plots and animation
app.py                               Streamlit experiment interface
tests/cpp/                           native solver test
tests/python/                        API, format, cadence, app and plot tests
scripts/validate_legacy_data.py      legacy text-file validator
data/README.md                       supplied-data checksums and provenance
legacy/source/                       normalized original sources
```

## Install

Requirements: Python 3.10+, CMake 3.18+, and a C++17 compiler.

```bash
python -m venv .venv
source .venv/bin/activate              # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
pip install -e ".[app,test]"
pytest
```

## Run with the supplied initial state

The input must remain the original four-column, 66,049-row format. Renaming
`init(1).dat` to `init.dat` is sufficient; no conversion is performed.

```bash
solar-dynamo --init /path/to/init.dat --steps 200 --output-dir output
```

Omit `--steps` to run all `tmax/dt = 200,000` steps:

```bash
solar-dynamo --init /path/to/init.dat --output-dir output
```

Use the original `irelax=0` analytic branch without a file:

```bash
solar-dynamo --analytic-init --steps 200 --output-dir output
```

## Python/pybind11 API

```python
from solar_dynamo import Model, Parameters

p = Parameters()                       # exact v3 defaults
model = Model(p)
model.initialize("init.dat")
model.step_with_output(200, "output") # preserves legacy cadence
state = model.snapshot()               # NumPy arrays, shape (257, 257)
model.write_final_outputs("output")
```

`step()` performs the same numerical update without disk I/O. Use
`step_with_output()` or `run()` when legacy file behavior is required. The GIL
is released during initialization and stepping.

## Exact file contracts

| File | Mode | Columns | Behavior |
|---|---|---|---|
| `init.dat` | input | `theta radius r*sin(theta)*u ub` | exactly 257×257 rows |
| `omg_lat.dat` | output | `theta time omega` | append every 200 steps |
| `diffrot.dat` during run | output | `theta radius omega time` | overwrite every 50 steps |
| `diffrot.dat` after run | output | `theta radius omega` | final overwrite |
| `final.dat` | output | `theta radius omega r*sin(theta)*u ub` | final overwrite |

Outputs occur before the Fortran `t=t+dt` statement. Thus the first
`omg_lat.dat` snapshot is labeled `0.0099500`, exactly as in the supplied file.

## Tests

```bash
pytest
cmake -S . -B build -DSOLAR_DYNAMO_BUILD_PYTHON=OFF -DSOLAR_DYNAMO_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Tests cover the exact grid, analytic initialization, boundaries, step/time
accounting, missing input handling, output row/column contracts, pre-increment
output time, 50/200-step cadence, serialization, visualization, Streamlit, and
the native C++ API.

## Streamlit

```bash
streamlit run app.py
```

The app supports the exact analytic initialization or an uploaded legacy
`init.dat`. It intentionally exposes only v3 parameters; grid size and method
cannot be changed because doing so would cease to be a faithful port.

For deployment, push the repository to GitHub and select `app.py` on Streamlit
Community Cloud. GitHub Pages cannot execute a Python server.

## Reference-output caveat

The supplied `omg_lat(1).dat` is three appended runs. The segments contain 442,
95, and 474 snapshots. The latter two have 0.008 spacing, implying a different
time step from the attached v3 source. `diffrot(1).dat` belongs to the final
segment at `t=3.7898`. These files establish formats and plausibility, but cannot
serve as a strict single-configuration parity baseline. Even the first segment,
which has v3-compatible spacing, differs from the C++ step-200 latitude field by
up to 2.2832 in Ω. It was therefore not accepted as a passing parity fixture.
A fresh output from the exact attached v3 source, compiler command, and unchanged
`init.dat` is still needed to establish numerical equality tolerances.

## Resume description

> Re-engineered a 998-line solar-dynamo Fortran program as an object-oriented
> C++17 simulation module while preserving its expanded ADI scheme, coupled
> magnetic/rotation updates, Legendre boundary treatment, and legacy I/O;
> exposed zero-copy-shaped scientific state through pybind11/NumPy and added
> automated native/Python tests, Streamlit visualization, Docker, and CI.

MIT licensed.
