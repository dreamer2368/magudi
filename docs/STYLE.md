# Documentation style guide

This document fixes the in-source documentation conventions for `magudi`. All annotation work targets two renderers:

- **FORD** for Fortran (parses `!>`, `!!`, `!<` markers; renders Markdown + MathJax)
- **Sphinx + Napoleon** for Python (parses numpy-style docstrings; renders reST)

Read this once and use the [exemplar files](#exemplars) as templates.

---

## 1. Fortran (FORD)

### Comment markers

| Marker | Position | Use |
|--------|----------|-----|
| `!>`   | preceding | Start of a doc block (above procedure / type / module). |
| `!!`   | continuation | Subsequent lines of a doc block. Markdown body. |
| `!<`   | trailing  | Inline doc on the line above (used for type fields and arguments). |

Plain `!` comments are **implementation notes** — FORD ignores them. Keep using them inside procedure bodies. The existing `! <<< Section >>>` dividers in `src/*Impl.f90` files are fine and stay as-is.

### Where docs live in this codebase

The split-file convention (`include/<Name>.f90` declarations, `src/<Name>Impl.f90` implementations) means doc blocks go **next to the abstract interface in `include/`**, not next to the implementation body. Reasons:

- `include/` is what users read to discover the API.
- Putting docs there avoids duplication between the `interface` block and the `subroutine` body.
- FORD picks them up on the interface; cross-references work either way.

The implementation file gets a single file-level header pointing back:

```fortran
!> @file
!! Implementation of [[Region_mod]]. See `include/Region.f90` for API documentation.

#include "config.h"

module RegionImpl
...
```

### Module header

Every `include/<Name>.f90` starts with a module-level block. Keep it brief — one-line summary, one-paragraph elaboration, key types/procedures inline-referenced with `[[Name]]`.

```fortran
#include "config.h"

!> @brief MPI-aware abstract base for boundary patches.
!!
!! A `t_Patch` is a contiguous index subset of a single block on which a
!! simultaneous-approximation-term (SAT) boundary closure or other localized
!! operator is applied. Concrete subclasses (e.g. [[SpongePatch_mod:t_SpongePatch]],
!! [[FarFieldPatch_mod:t_FarFieldPatch]]) implement the deferred `setup`,
!! `cleanup`, `verifyUsage`, and `updateRhs` procedures.
!!
!! See [[Patch_factory]] for instantiation and `bc.dat` for input format.
module Patch_mod
```

### Derived-type fields

Use **trailing `!<`** on each field. One line, ends in a period. Don't restate the Fortran type — say what the field *means* in the model.

```fortran
type, public :: t_Region

   type(t_Grid), allocatable :: grids(:)               !< Per-block grid container; one entry per block owned by this process.
   type(t_State), allocatable :: states(:)             !< Per-block solution state, parallel to `grids`.
   integer :: comm = MPI_COMM_NULL                     !< Communicator spanning all processes that participate in this region.
   integer :: timestep = 0                             !< Current physical time step counter; advanced by [[t_Solver:runForward]].

 contains

   procedure, pass :: setup => setupRegion             !< See [[setupRegion]].
   ...
end type t_Region
```

If a field's name is fully self-describing (e.g. `cachedValue`), a one-word doc (`!< Cached.`) is acceptable but `!<` is still required — it tells future readers "I considered this and there's nothing more to say."

### Procedures

Doc block goes **above the `subroutine`/`function` line in the interface block in `include/`**. Required tags:

- One-line summary (no leading `@brief` needed; FORD treats the first sentence as the summary).
- `@param[intent] name` for each non-`this` argument. `intent` is one of `in`, `out`, `inout`. Skip the intent if it matches the Fortran declaration (FORD reads it).
- `@return name` for functions.

Optional:

- A second paragraph of body prose for non-trivial procedures.
- `@note`, `@warning`, math (`$...$` inline, `$$...$$` displayed), Markdown lists.
- `@see [[OtherProc]]` for cross-refs (or just inline `[[OtherProc]]` in prose).

Example:

```fortran
interface

   !> Initialize the region: allocate per-block grids and states, parse the
   !! decomposition map, and split MPI communicators per block.
   !!
   !! Must be called before any other `t_Region` method except `cleanup`.
   !! Idempotent — calls `cleanup` first.
   !!
   !! @param[inout] this        The region to initialize.
   !! @param[in]    comm        MPI communicator spanning all processes that
   !!                           will participate in this region.
   !! @param[in]    globalGridSizes  Shape `(nDimensions, nBlocks)` describing
   !!                                the global extent of each block.
   !! @param[in]    simulationFlags  Optional pre-built simulation flags. When
   !!                                absent, defaults are read from the input
   !!                                file via [[t_SimulationFlags:initialize]].
   !! @param[in]    solverOptions    Optional pre-built solver options. When
   !!                                absent, defaults are constructed from
   !!                                `simulationFlags`.
   !! @param[in]    verbose     Print per-block decomposition info to stdout.
   !!                           Defaults to `.true.`.
   subroutine setupRegion(this, comm, globalGridSizes, simulationFlags,                   &
        solverOptions, verbose)
     ...
   end subroutine setupRegion

end interface
```

### Abstract / deferred procedures

Annotate the `abstract interface` block in the base type's `include/<Base>.f90`. The doc describes the **contract** — what every subclass must guarantee — not any specific implementation.

```fortran
abstract interface

   !> Apply the patch's contribution to the right-hand-side of the governing
   !! equations on a single block.
   !!
   !! Called once per block per RK stage from [[t_Region:computeRhs]].
   !! Subclasses must respect the `mode` selector (forward / adjoint /
   !! linearized) — the same patch object is used in all three passes.
   !!
   !! @param[in,out] this   The patch.
   !! @param[in]     mode   One of `FORWARD`, `ADJOINT`, `LINEARIZED` from
   !!                       [[Region_enum]].
   !! @param[in]     simulationFlags Solver-wide flags.
   !! @param[in]     solverOptions   Solver-wide options.
   !! @param[in]     grid   The block's grid.
   !! @param[in,out] state  The block's state; `state%rightHandSide` is updated.
   subroutine updateRhs(this, mode, simulationFlags, solverOptions, grid, state)
     ...
   end subroutine updateRhs

end interface
```

### Math

FORD uses MathJax. Inline: `\(...\)` or `$...$`. Displayed: `$$...$$`. Example:

```fortran
!> Compute the spatial inner product \( \langle u, v \rangle = \int_\Omega u v \, dV \)
!! using the SBP norm matrix.
```

### Cross-references

- `[[t_Region]]` → links to the type.
- `[[Region_mod:t_Region]]` → fully qualified.
- `[[setupRegion]]` → links to the procedure.
- Use these liberally in prose; they make the rendered HTML navigable.

### What NOT to do

- **Don't restate the signature.** FORD already shows `intent(in) :: comm` — your doc should explain what `comm` *is*, not that it's an integer.
- **Don't narrate the implementation.** Inline `!` comments in the body do that. Doc blocks describe contract and intent.
- **Don't write filler.** "This subroutine does X" → just "Does X". One sentence summaries are fine.
- **Don't add `@author`, `@date`, or version stamps.** No file-header convention exists in this repo and we're not adding one. Git is authoritative for who/when.
- **Don't document private helpers in factory `connect*` dispatch routines.** A module-header doc on the factory is enough.

---

## 2. Python (Sphinx + Napoleon)

### Style: numpy

Numpy docstrings render well as both plain text and HTML, and they're terser than Google style for our argument-heavy methods.

### Module docstring

First statement after the `from ... import` block. One-line summary, one paragraph body, a `Notes` section if there's a non-obvious invariant.

```python
"""Driver for the multi-point conjugate-gradient optimization loop.

Coordinates forward, adjoint, line-minimization, and bracketing stages
across multiple time segments, persisting state to an HDF5 log between
invocations so the loop can be resumed by the surrounding shell driver.

Notes
-----
This module is invoked by ``checkGradientAccuracy.py`` and the bash
optimization scripts in this directory; it does not run a complete
optimization itself but emits shell command files consumed by the driver.
"""
```

### Class docstring

Right after `class Foo:`. Sections (in order, all optional except the summary):

- One-line summary.
- Body paragraph (optional).
- `Parameters` — for `__init__` arguments.
- `Attributes` — for class- or instance-level state worth documenting.
- `Notes` — invariants, gotchas.

```python
class Optimizer:
    """Manage the multi-point optimization loop.

    Holds the running state of the conjugate-gradient line-search across
    sessions and emits shell command files that drive ``forward``,
    ``adjoint``, and the vector-arithmetic helpers in ``bin/``.

    Parameters
    ----------
    config : InputParser
        Parsed YAML config. See ``example.yml`` for the schema.

    Attributes
    ----------
    stage : Stage
        The next action scheduled for the surrounding driver to execute.
    result : Result
        Outcome of the most recently executed action.
    hyperStep, cgStep, lineStep : int
        Penalty-method, conjugate-gradient, and line-search counters.
    isInitial : bool
        ``True`` on the very first invocation; toggled on the first save.
    """
```

### Method docstring

```python
def schedule(self):
    """Advance the optimizer state machine by one step.

    Inspects ``self.stage`` and ``self.result``, emits the next batch
    of shell commands to the global command file, and updates persisted
    state. Idempotent within a single ``(stage, result)`` pair.

    Raises
    ------
    RuntimeError
        If ``self.stop`` is set, signalling external intervention.
    """
```

For methods that return a value:

```python
def nextMnbrak(self):
    """Take one step of the minimum-bracketing search.

    Returns
    -------
    bool
        ``True`` when ``a < b < c`` and ``Jb < min(Ja, Jc)`` — the bracket
        is established and line minimization can begin. ``False`` if
        another forward run is needed; the command file has been
        rewritten accordingly.
    """
```

### Type hints

Encouraged in `utils/optimization_ver3/` (Python 3 only). Don't add to legacy Py2 code. When hints exist, you can drop the type from the `Parameters` block — but keep the description.

### What NOT to do

- **Don't echo argument names without describing them.** `config : InputParser` with no description adds no signal.
- **Don't write `:param:` / `:returns:` reST style** — we standardized on numpy. Sphinx Napoleon converts numpy → reST automatically.
- **Don't document `__repr__`, `__init__` separately** — class-level `Parameters` covers `__init__`; numpy-style doesn't put a docstring on `__init__` itself.

---

## 3. Scope per item

| Scope | Required |
|-------|----------|
| Module / file | Header block (1 line + 1 paragraph) |
| Public derived type / class | One-line summary + every field/attribute |
| Public procedure / method | Summary + every argument + return |
| Abstract / deferred procedure | Summary describing the **contract** |
| Private procedure / helper | One-line `!>` summary; no `@param` needed |
| Concrete subclass that's a thin specialization | Module header only is acceptable |

Helpers (e.g. `CNSHelper`, `RhsHelper`, `MPIHelper`) that are not part of the user-facing API get a module header explaining their role and one-line `!>` summaries on each public procedure. Their argument-by-argument annotation is deferred indefinitely.

---

## 4. Exemplars

The following files are documented end-to-end and serve as the canonical templates. When in doubt, copy the structure from these:

- [include/Region.f90](../include/Region.f90) — module header, public type with field docs, public procedure interfaces with `@param`/`@return`.
- [src/RegionImpl.f90](../src/RegionImpl.f90) — file-level pointer back to `include/`; selected procedures annotated.
- [include/Patch.f90](../include/Patch.f90) — abstract base type, deferred-procedure interfaces with contract docs.
- [utils/optimization_ver3/optimizer.py](../utils/optimization_ver3/optimizer.py) — module docstring, three classes, all public methods.

Read those four files before annotating anything else.

---

## 5. Rendering and verification

All rendering happens inside the project container (see [docker/Dockerfile](../docker/Dockerfile)). From the repo root:

```bash
docker run --rm -v "$PWD":/work -w /work \
    ghcr.io/dreamer2368/magudi/magudi_env:arm64 \
    bash -c "ford docs/ford.md && sphinx-build docs/sphinx docs/sphinx/_build"
```

Then open `_ford_build/index.html` and `docs/sphinx/_build/index.html` in your browser.

The progress tracker for the full annotation pass lives in [docs/ANNOTATION_TRACKER.md](ANNOTATION_TRACKER.md).
