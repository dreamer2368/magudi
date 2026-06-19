# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.

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

- `run_msgrad.sh` — multi-segment gradient-accuracy check via the `OneDWave` example.
- `run_parallel.sh` — restart-equivalence test of the parallel TAO L-BFGS loop.

Both `pip install ./utils/magudi_utils` (the `magudi_utils` package) and drive its console scripts `magudi-msgrad` / `magudi-optim` against `examples/OneDWave/` staged into the build tree.

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

Auxiliary executables (post-processing such as `q_criterion`, `vorticity_dilatation`, grid generators) are built into `build/utils/`.

Python utilities live in [utils/magudi_utils/](utils/magudi_utils/) — an installable package (`pip install ./utils/magudi_utils`) that bundles PLOT3D file I/O (`plot3dnasa`, `PLOT3D`), mesh helpers (`SummationByParts`, `RoundJet`, `SingleBlockCartesian`), the FWH solver (`fwhsolver`), matplotlib helpers, and the PETSc/TAO L-BFGS optimization driver (`magudi-optim` / `magudi-msgrad` console scripts, used by CI). Examples and tests `from magudi_utils import plot3dnasa as p3d`. Architecture for the optimizer: [utils/magudi_utils/DESIGN.md](utils/magudi_utils/DESIGN.md).

Superseded code lives under [utils/legacy/](utils/legacy/) (`inexactnewton/`, `optimization/`, `optimization_ver2/`, `optimization_ver3/`, the residue `python/`); see [utils/legacy/README.md](utils/legacy/README.md). Not on the active CI path — kept for in-flight production runs only.

### Examples ([examples/](examples/))

Each subdirectory is a self-contained simulation: a Python `config.py` that generates the PLOT3D grid/initial-condition/mollifier files (run with `python2 config.py`), a `magudi.inp` parameter file, and a `bc.dat` boundary-condition spec. [examples/magudi.inp](examples/magudi.inp) is the documented reference for all input flags.

## Conventions

- Every Fortran source begins with `#include "config.h"` to get `SCALAR_KIND`, `SCALAR_TYPE`, `STRING_LENGTH`, `SAFE_DEALLOCATE`, the debug-only `assert()`, and the `PURE_SUBROUTINE`/`PURE_FUNCTION` macros that vanish in debug mode.
- Module names use the `_mod` / `_factory` / `_enum` suffixes (e.g. `Region_mod`, `Patch_factory`, `Region_enum`).
- Type names are prefixed with `t_` (e.g. `t_Region`, `t_Patch`).
- All real/complex storage uses `SCALAR_TYPE` and `SCALAR_KIND` rather than hard-coded kinds, so the same source compiles for any of the six numeric precisions.

## Documentation

A future effort will publish API documentation to GitHub Pages using FORD (Fortran) and Sphinx (Python). Conventions for that effort are sketched in [docs/STYLE.md](docs/STYLE.md); refer to it before adding any in-source documentation so the eventual rendering pass works without rework.
