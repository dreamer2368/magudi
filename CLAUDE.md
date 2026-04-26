# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

`magudi` is an MPI-parallel compressible Navier-Stokes solver in Fortran 2008 (with a small amount of C) targeted at flow control and adjoint-based optimization. It uses summation-by-parts (SBP) finite differences with simultaneous-approximation-term (SAT) boundary closures, reads/writes PLOT3D-format multi-block grids, and produces forward, adjoint, and linearized solutions.

## Container

The project ships a development container — `ghcr.io/dreamer2368/magudi/magudi_env:arm64` (image definition: [docker/Dockerfile](docker/Dockerfile)). It contains the full toolchain (`gfortran`, `mpich`, `cmake`, `valgrind`, `hdf5-tools`) and the Python deps used by the optimization framework (`numpy`, `scipy`, `tables`, `pandas`, `PyYAML`, `h5py`). All CI jobs in [.github/workflows/test.yaml](.github/workflows/test.yaml) run inside this image. **All build and test verification should run inside this container** so behavior matches CI:

```
docker pull ghcr.io/dreamer2368/magudi/magudi_env:arm64
docker run --rm -it -v "$PWD":/work -w /work \
    ghcr.io/dreamer2368/magudi/magudi_env:arm64 bash
```

The image is rebuilt and pushed by [.github/workflows/docker.yml](.github/workflows/docker.yml) whenever files in `docker/` change. Note: that workflow currently builds `linux/amd64` only; the `:arm64` tag is maintained out-of-band — when changing the Dockerfile, plan to refresh the arm64 image as well.

## Build

The project uses CMake. There is no in-tree build — create a separate build directory:

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=debug -DCMAKE_Fortran_COMPILER=mpif90 ..
make -j
```

Common CMake variables:

- `CMAKE_BUILD_TYPE` — `Debug` (default), `Release`, or `RelWithDebInfo`. Debug defines the `DEBUG` macro which enables `assert()` and disables `pure` qualifiers.
- `SCALAR_TYPE` — `real32_iso`, `real64_iso` (default), `real128_iso`, `binary128_IEEE754`, `complex`, `double_complex`. Selected via `config.h` macros (`SCALAR_KIND`, `SCALAR_TYPE`, `SCALAR_TYPE_MPI`); used everywhere via `#include "config.h"`.
- `PLOT3D_FORMAT` — `regular` (default) or `extended` (large-file 64-bit offsets).

CMake forces `CMAKE_Fortran_COMPILER=mpif90` unless the compiler ID is XL. Supported compiler families: `gfortran` (with `-fallow-argument-mismatch` for GCC ≥ 10), `ifort`, `bgxlf2008_r`. `make purge` cleans the build tree.

## Tests

Built into the same build tree. From `build/`:

```
CTEST_OUTPUT_ON_FAILURE=TRUE make test     # run all
ctest -R <name>                            # filter by regex
ctest -V -R <name>                         # verbose, single test
```

Test categories (see [test/CMakeLists.txt](test/CMakeLists.txt) and `test/{adjoint_relation,linearized_relation}/CMakeLists.txt`):

- **Serial unit tests** (`add_serial_test`): plain executables run directly.
- **MPI tests** (`add_mpi_test`): run via `mpirun -np {1,2}`.
- **Script tests** (`add_script_test`, only in `test/adjoint_relation`): the `.f90` test driver is paired with a `.sh` script that sets up grids/inputs and invokes it; CMake runs the shell script.

The `test/adjoint_relation` and `test/linearized_relation` suites verify that adjoint/linearized operators are consistent with forward operators (e.g. `<A x, y> == <x, A^T y>`).

End-to-end CI also runs two driver scripts in `.github/workflows/`:

- `optim_grad_test.sh` — gradient-accuracy check via the `OneDWave` example.
- `optim_test.sh` — full multi-point optimization loop.

Both copy `utils/optimization_ver3/*.py` and `examples/OneDWave/*` into the build tree, then drive `python3 checkGradientAccuracy.py optim.yml --mode {setup,schedule,log_result}` in a loop.

## Architecture

### Source layout

Module declarations and abstract-interface blocks live in [include/](include/) (`Foo.f90`); implementations live in [src/](src/) (`FooImpl.f90`). The two are compiled together into a single object library `magudiObj` ([CMakeLists.txt:90-231](CMakeLists.txt#L90-L231)), and every executable in [bin/](bin/), [test/](test/), and [utils/](utils/) links against `$<TARGET_OBJECTS:magudiObj>` plus a single `<name>.f90` driver.

When adding a new module, list **both** the `include/Foo.f90` and `src/FooImpl.f90` in `magudiObj_SOURCES` in the top-level [CMakeLists.txt](CMakeLists.txt). Adding a binary means a one-line `add_executable(... $<TARGET_OBJECTS:magudiObj> ...)` in the relevant `CMakeLists.txt`.

### Core types

- `t_Region` ([include/Region.f90](include/Region.f90)) — the central container. Owns the array of grids, states, patch factories, solver options, simulation flags, MPI communicators, and the optional levelset factory. All forward/adjoint/linearized operations are methods on (or take) a `t_Region`.
- `t_Solver` ([include/Solver.f90](include/Solver.f90)) — drives time integration; its `runForward`, `runAdjoint`, `runLinearized`, and `checkGradientAccuracy` methods are what the binaries call.
- `Region_enum` defines `FORWARD = +1`, `ADJOINT = -1`, `LINEARIZED = 0` — a mode selector passed throughout the code.

### Abstract types + factories

The polymorphic extension pattern is used pervasively. An abstract base in `include/` defines deferred procedures; concrete subclasses live alongside it; a factory in `include/<Base>Factory.f90` instantiates the right subclass from a string key parsed from `magudi.inp`/`bc.dat`. The base/subclass/factory triples are:

- **Patch** ([Patch.f90](include/Patch.f90)) → `SpongePatch`, `FarFieldPatch`, `ImpenetrableWall`, `IsothermalWall`, `AdiabaticWall`, `CostTargetPatch`, `ActuatorPatch`, `BlockInterfacePatch`, `ProbePatch`, `KolmogorovForcingPatch`, `ImmersedBoundaryPatch`, `SolenoidalExcitationPatch`, `JetExcitationPatch` (factory: `PatchFactory`).
- **Controller** ([Controller.f90](include/Controller.f90)) → `ThermalActuator`, `MomentumActuator`, `GenericActuator` (factory: `ControllerFactory`). Controllers compute sensitivities and update control forcing.
- **Functional** ([Functional.f90](include/Functional.f90)) → `AcousticNoise`, `PressureDrag`, `DragForce`, `ReynoldsStress`, `LighthillSource`, `LighthillTensorComponent`, `DensityGradient` (factory: `FunctionalFactory`). Functionals compute the cost J and its adjoint forcing.
- **TimeIntegrator** → `RK4Integrator`, `JamesonRK3Integrator` (factory: `TimeIntegratorFactory`).
- **ReverseMigrator** → `UniformCheckpointer` (factory: `ReverseMigratorFactory`).
- **MappingFunction** → `UniformMap` (factory: `MappingFunctionFactory`).
- **LevelsetFactory** ([LevelsetFactory.f90](include/LevelsetFactory.f90)) → `SinusoidalWallLevelset`, `StokesSecondWallLevelset` (used for immersed boundary).

When extending one of these, add the new subclass to the factory's dispatch (`connect*`) routine in `src/<Factory>Impl.f90`.

### Binaries ([bin/](bin/))

The primary executables are `forward`, `adjoint`, `linearized`, and `gradient_accuracy` (a finite-difference adjoint check). Vector arithmetic helpers (`zaxpy`, `zwxmwy`, `zxdoty`, `spatial_inner_product`, `qfile_zaxpy`, `control_space_norm`) operate on the binary `.dat` control/gradient vectors and on PLOT3D `.q` files. `slice_control_forcing` / `paste_control_forcing` / `terminal_objective` / `patchup_qfile` are utilities for multi-segment optimization.

### Utilities ([utils/](utils/))

Auxiliary executables (post-processing such as `q_criterion`, `vorticity_dilatation`, grid generators) are built into `build/utils/`. The `utils/python/` directory contains Python helpers — note that `plot3dnasa.py` and `examples/*/config.py` are **Python 2** (incompatible with Python 3); the optimization drivers are Python 3.

Three generations of Python optimization frameworks coexist: `utils/optimization/` (oldest), `utils/optimization_ver2/`, and `utils/optimization_ver3/` (current — YAML-configured, used by CI). Pull from `optimization_ver3` for new work; example config in [utils/optimization_ver3/example.yml](utils/optimization_ver3/example.yml).

### Examples ([examples/](examples/))

Each subdirectory is a self-contained simulation: a Python `config.py` that generates the PLOT3D grid/initial-condition/mollifier files (run with `python2 config.py`), a `magudi.inp` parameter file, and a `bc.dat` boundary-condition spec. [examples/magudi.inp](examples/magudi.inp) is the documented reference for all input flags.

## Conventions

- Every Fortran source begins with `#include "config.h"` to get `SCALAR_KIND`, `SCALAR_TYPE`, `STRING_LENGTH`, `SAFE_DEALLOCATE`, the debug-only `assert()`, and the `PURE_SUBROUTINE`/`PURE_FUNCTION` macros that vanish in debug mode.
- Module names use the `_mod` / `_factory` / `_enum` suffixes (e.g. `Region_mod`, `Patch_factory`, `Region_enum`).
- Type names are prefixed with `t_` (e.g. `t_Region`, `t_Patch`).
- All real/complex storage uses `SCALAR_TYPE` and `SCALAR_KIND` rather than hard-coded kinds, so the same source compiles for any of the six numeric precisions.

## Documentation

A future effort will publish API documentation to GitHub Pages using FORD (Fortran) and Sphinx (Python). Conventions for that effort are sketched in [docs/STYLE.md](docs/STYLE.md); refer to it before adding any in-source documentation so the eventual rendering pass works without rework.
