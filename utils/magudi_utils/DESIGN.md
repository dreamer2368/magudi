# magudi_utils optimizer — design document

This document specifies the architecture of the optimization framework now shipped inside the `magudi_utils` package, replacing the legacy [utils/legacy/optimization_ver3/](../legacy/optimization_ver3/).

## 1. Why ver4 exists

ver3 has three structural limitations that grew out of its origin as a PDE-driver harness:

1. **Control-space algebra lives in Fortran executables** (`zaxpy`, `zxdoty`, `zwxmwy`, `qfile_zaxpy`, `spatial_inner_product`, `patchup_qfile`, `slice_control_forcing`, `paste_control_forcing`) called from generated bash command files. This makes it hard to plug in alternative optimizers — every algorithmic choice has to be re-expressed in those primitives.
2. **The optimizer is hand-written**: Polak–Ribière nonlinear CG with golden-section + parabolic-interpolation line search, an explicit augmented-Lagrangian outer loop, and a five-stage state machine over `(hyperStep, cgStep, lineStep)`. Adding L-BFGS, trust region, or any other gradient method means rewriting the inner loop.
3. **Each forward/adjoint sweep is fragmented across `Nsplit` per-segment subdirectories** (`x0/0/`, `x0/1/`, …) requiring per-segment `magudi.inp` mutation via `setOption.sh`, and a multi-job dependency graph in the bash command file. At small problem scale this is manageable; at the 100 GB+ control-vector scale ver4 targets, the orchestration becomes the dominant source of fragility.

ver4 addresses all three by:

- Moving control-space algebra into Python on top of **PETSc/petsc4py distributed `Vec`s**.
- Using **PETSc/TAO's `lmvm`** (driven by `tao.solve()`) as the L-BFGS optimizer; the driver only provides the objective/gradient callback, the variable transformation, and a wall-clock signal handler. Future-extensible to a manually-owned outer loop on top of raw `MatLMVM` if production constraints demand it (see §6).
- Introducing two new Fortran executables, **`msforward`** and **`msadjoint`**, that drive a `t_Solver` over all `Nsplit` segments concatenated, eliminating the per-segment subdirectory orchestration.

## 2. Constraints and finalized choices

### Problem-scale constraints

- **Control vector size**: ≥ 100 GB at production scale (concatenated control forcing + per-segment initial conditions). Distributed parameter storage is mandatory; single-rank numpy/scipy is ruled out.
- **Async batch evaluation**: each forward/adjoint sweep is long (hours) and the parameter space is huge (TB-class allocations). Every objective / gradient / line-search-direction evaluation must persist its result to disk; a later process re-reads it. This is the same constraint that drove ver3's `--mode setup / schedule / log_result` triad.
- **One SLURM allocation may hold one or several iterations** depending on the application. The framework is built for the optimistic case ("Case A-long" in §6: one Python process spans many outer iterations within one allocation) and is documented to extend to the more constrained cases (one outer iteration per allocation; or forward and adjoint in separate allocations) if production hits those constraints.

### Finalized design choices

| # | Decision | Rationale |
|---|---|---|
| 1 | **PETSc/petsc4py + TAO `lmvm` driven by `tao.solve()`** | Native MPI-distributed Vec/Mat; no driver-side memory ceiling. Shortest implementation path — TAO owns the outer loop, line search, two-loop recursion, history, and convergence test. Drops to raw `MatLMVM` with manual outer loop only if production hits the one-iteration-per-allocation constraint (see §6). |
| 2 | **L-BFGS, history depth m = 3** | At 100 GB Vec size, m = 3 → (2m+2) · N ≈ 800 GB persistent history. Fits any cluster already running this scale. m is tunable; recommend benchmarking against m = 5 and L-SR1 before committing. |
| 3 | **Variable transformation `y = D x` with `M = DᵀD`, applied in Python** | Provably path-equivalent to true M-inner-product L-BFGS — gives mesh-independent convergence. For magudi's diagonal SBP norm, `D = diag(√m_i)`; transformation is one `VecPointwiseMult`, O(N). The Fortran exes keep ver3 conventions (no Riesz mapping inside `msadjoint`); Python wraps the I/O. |
| 4 | **Build `msforward` / `msadjoint`** | One `mpirun` per objective/gradient eval. Eliminates per-segment subdirectory orchestration. Concatenated PETSc binary I/O between Python and Fortran avoids the 100 GB "stitch Nsplit per-segment dumps" shuffle. |
| 5 | **Case A-long primary**: one Python process per SLURM allocation drives `tao.solve()` over many outer iterations | Shortest implementation path — TAO owns the line search, two-loop recursion, history, and convergence test. Signal handler at wall-clock boundary checkpoints (current `y`, MatLMVM history) and exits cleanly. Bash wrapper resubmits if not converged. Future-extensible to **Case A** (one outer iteration per allocation, raw `MatLMVM` with manual outer loop) and **Case B** (forward and adjoint on separate allocations, BB step) if production hits those constraints. |
| 6 | **Matching-condition penalty baked into J** returned by `msforward` | Same as ver3 functionally. Skip TAOALMM for now; reconsider if penalty-method convergence stalls. |

## 3. Architecture overview

```
                  ┌──────────────────────────────────────────────────┐
                  │              SLURM / Flux allocation              │
                  │                                                   │
                  │   mpirun -n N python3 optim.py optim.yml         │
                  │                                                   │
                  │   ┌──────────────────────────────────────────┐    │
                  │   │            optim.py  (driver)            │    │
                  │   │                                          │    │
                  │   │   load y.petsc, history.petsc            │    │
                  │   │   tao.setObjectiveGradient(callback)     │    │
                  │   │   tao.setLMVMMatrix(history)             │    │
                  │   │   register signal handler @ wall-clock   │    │
                  │   │                                          │    │
                  │   │   tao.solve(y)  ◄─────────────────┐     │    │
                  │   │       (drives many outer iters,    │     │    │
                  │   │        TAO owns line search, CG    │     │    │
                  │   │        history, convergence test)  │     │    │
                  │   │                                    │     │    │
                  │   │   ┌────────────────────────────┐   │     │    │
                  │   │   │  callback(y, g_y):         │───┘     │    │
                  │   │   │    x = D⁻¹ · y             │         │    │
                  │   │   │    write x.petsc           │         │    │
                  │   │   │    run msforward + msadjoint│        │    │
                  │   │   │    read J, raw_grad        │         │    │
                  │   │   │    g_y = D⁻ᵀ · raw_grad    │         │    │
                  │   │   │    return J, g_y           │         │    │
                  │   │   └────────────────────────────┘         │    │
                  │   │                                          │    │
                  │   │   on signal OR convergence:              │    │
                  │   │     extract y, MatLMVM history → disk    │    │
                  │   │     exit                                 │    │
                  │   └──────────────────────────────────────────┘    │
                  └──────────────────────────────────────────────────┘
                                         │
                                         ▼
                  ┌──────────────────────────────────────────────────┐
                  │   Surrounding bash loop (OPT.sh / OPT.flux)      │
                  │   resubmits next allocation if not converged     │
                  └──────────────────────────────────────────────────┘
```

Three things are happening:

1. **PETSc/TAO is the optimizer.** `tao.solve()` runs the L-BFGS outer loop, line search, two-loop recursion, history management, and convergence test. The Python driver only provides the objective/gradient callback, the variable transformation, and the wall-clock signal handler.
2. **The variable transformation is a pair of pointwise ops inside the callback.** TAO sees only `y`-coordinates (Euclidean inner product) and `g_y = ∇_y h`. The Fortran exes see physical `x` and write the dual derivative `dJ/dx`. The callback is the only place that knows about `D`.
3. **State persists via PETSc binary `MatView` / `Vec.view`** at allocation boundaries. Inside `tao.solve()`, all state lives in TAO's in-memory structures — no per-iteration disk I/O for optimizer state, only for `x`/`raw_grad` going to/from the Fortran exes. At allocation boundary (signal or convergence), the driver extracts the current iterate and the `MatLMVM` history matrix and writes them out.

## 4. Directory layout

For an example case with `global_prefix: OneDWave`:

```
<run-root>/
│
├── magudi.inp                          # global input — single source of truth
├── bc.dat
├── OneDWave.xyz                        # PLOT3D grid
├── OneDWave.control_mollifier.f
├── OneDWave.target_mollifier.f
├── optim.yml                           # ver4 config (YAML)
│
│   ── binary symlinks ──
├── msforward            -> ../bin/msforward
├── msadjoint            -> ../bin/msadjoint
│
│   ── canonical optimizer state (persisted at allocation boundary) ──
├── y.petsc                             # current y = D · x  — TAO's iterate
│                                       #   (saved on signal/convergence; restored on next start)
│
│   ── L-BFGS state (persisted at allocation boundary) ──
├── history.petsc                       # MatLMVM (s_k, y_k) pairs in y-space, m=3
│                                       #   Extracted via tao.getLMVMMatrix() at boundary
│
│   ── per-callback scratch (rewritten on every objective/gradient eval) ──
├── x.petsc                             # x = D⁻¹ · y written before each msforward call
│
│   ── metric (read once at setup, reused thereafter) ──
├── M_diag.petsc                        # diagonal of SBP norm matrix M
├── D_diag.petsc                        # √M elementwise (precomputed for speed)
├── D_inv_diag.petsc                    # 1/√M elementwise
│
│   ── outputs from last (msforward + msadjoint) ──
├── J.txt                               # scalar: total loss (incl. matching-condition penalty)
├── raw_grad.petsc                      # dual derivative dJ/dx, no metric applied
├── J_per_segment.txt                   # diagnostic sublosses (TSV; not consumed by optimizer)
│
│   ── persistent logs ──
├── history.h5                          # iteration history: (iter, J, ‖∇h‖, step, accepted, ...)
├── journal.txt                         # python logging output
├── line_search_state.json              # tiny: bracket points, Armijo backtracking counter
│
│   ── per-segment forward/adjoint snapshots (written by msforward/msadjoint) ──
├── snapshots/
│   ├── OneDWave-0-00030000.q
│   ├── OneDWave-0-00032400.q
│   ├── OneDWave-1-00032400.q
│   ├── OneDWave-1-00034800.q
│   └── ...
│
└── out/                                # captured stdout/stderr from msforward/msadjoint runs
    ├── msforward.iter0042.out
    └── msadjoint.iter0042.out
```

Compare against ver3's `x0/ a/ b/ c/ x/` × per-segment subdirs × auxiliary dirs (`out/ diff/ grad/ cg/ previous/ txt/ lagrangian/ baseline/ diffLog/ linminLog/`). ver4's flat layout reflects the absence of a five-tree line-search bracket (PETSc handles it differently) and the elimination of per-segment subdirectories.

## 5. Python ↔ Fortran interface

### 5.1 `msforward`

**Inputs (read by msforward)**:

| File | Contents |
|---|---|
| `magudi.inp` | Global solver settings. One copy, no per-segment mutation. |
| `bc.dat`, `OneDWave.xyz`, `*.f` | Grid + boundary conditions + mollifiers (unchanged from current). |
| `x.petsc` | PETSc binary Vec of the concatenated control vector (control forcing + per-segment ICs), length `NcontrolRegion · totalTimestep + Nsplit · NgridPoints`. Layout fixed by `optim.yml`. |
| CLI args | `--input magudi.inp --control x.petsc --output J.txt --segments-output J_per_segment.txt --snapshot-dir snapshots/` |

**Outputs (written by msforward)**:

| File | Contents |
|---|---|
| `J.txt` | Single scalar: total loss = Σ time-integral J per segment + matching-condition penalty + terminal J (if enabled). **This is the J the optimizer consumes.** |
| `J_per_segment.txt` | TSV: `segment_idx`, `time_integral_J`, `L2sq_mismatch`, `lagrangian_term`, `terminal_J`. Diagnostic only. |
| `snapshots/OneDWave-<k>-<timestep>.q` | Per-segment q-files at segment boundaries (needed by msadjoint). |

**Internals**: drives a `t_Solver` over `Nsplit` segments concatenated. Internally manages the segment-loop that ver3 expressed as a bash command file. Returns a single MPI-collective scalar to rank 0 for `J.txt`. Penalty-on-mismatch (matching-condition term) is computed inside Fortran and folded into `J`.

### 5.2 `msadjoint`

**Inputs (read by msadjoint)**:

| File | Contents |
|---|---|
| `magudi.inp`, grid, BC, mollifiers | Same as msforward. |
| `x.petsc` | Same control vector as msforward consumed. |
| `snapshots/*.q` | Forward snapshots written by the most recent msforward call. |
| CLI args | `--input magudi.inp --control x.petsc --snapshot-dir snapshots/ --output raw_grad.petsc` |

**Outputs (written by msadjoint)**:

| File | Contents |
|---|---|
| `raw_grad.petsc` | PETSc binary Vec of the **dual** derivative `dJ/dx`, same length and layout as `x.petsc`. **No `M⁻¹` or `D⁻ᵀ` applied** — Python handles the variable transformation. |

**Internals**: drives `Nsplit` adjoint sweeps backward through the snapshots, accumulating the gradient contributions from per-segment time-integral objectives, matching-condition penalty terms (with their derivative w.r.t. control and ICs), and terminal objective (if enabled). Output is one concatenated PETSc Vec.

### 5.3 Python-side variable transformation

The transformation `y = D x` with `D = diag(√m_i)` lives entirely in Python — specifically inside the TAO objective/gradient callback. TAO sees only `y` and `g_y = ∇_y h`; the Fortran exes see only physical `x` and write the raw dual derivative `dJ/dx`:

```python
# At setup (one-time): derive D and D_inv from M_diag, save for reuse
M_diag = PETSc.Vec().load(viewer('M_diag.petsc'))
D_diag = M_diag.copy(); D_diag.sqrtabs()
D_inv_diag = D_diag.copy(); D_inv_diag.reciprocal()
D_diag.view(viewer('D_diag.petsc'))
D_inv_diag.view(viewer('D_inv_diag.petsc'))

# Inside the TAO callback (one call per objective+gradient eval):
def fg(tao, y, g_y):
    # y → x for msforward
    x = y.duplicate(); x.pointwiseMult(y, D_inv_diag)   # x_i = y_i / √m_i
    x.view(viewer('x.petsc'))

    subprocess.run(['mpirun', '-n', N, './msforward', ...], check=True)
    subprocess.run(['mpirun', '-n', N, './msadjoint', ...], check=True)

    # raw dual gradient → ∇_y h
    raw_grad = PETSc.Vec().load(viewer('raw_grad.petsc'))
    g_y.pointwiseMult(raw_grad, D_inv_diag)            # (∇_y h)_i = (dJ/dx)_i / √m_i
                                                        # (note: D⁻ᵀ = D⁻¹ since D is diagonal)
    return float(open('J.txt').read())
```

The persisted optimizer state on disk (`y.petsc`, `history.petsc`) is in `y`-coordinates throughout, so resuming from a checkpoint requires no re-transformation — TAO loads `y.petsc` directly into its `Vec`.

Why this matters: M-inner-product L-BFGS gives **mesh-independent convergence** for PDE-constrained problems. Naive Riesz mapping (`grad_M = M⁻¹ · dJ/dx` with Euclidean L-BFGS on top) is *not* path-equivalent — the L-BFGS two-loop recursion's intermediate inner products would still be Euclidean, picking up unwanted `M⁻¹` factors. Variable transformation is the only first-class way to recover true M-inner-product L-BFGS using off-the-shelf Euclidean primitives.

### 5.4 The metric vector `M_diag.petsc`

Where it comes from: a one-time setup utility (call it `compute_norm`, analogous to whatever wrote `*.norm_<actuator>.dat` in ver3) walks the grid, computes the per-grid-point quadrature weights (the diagonal of the SBP norm matrix), and writes them as a PETSc Vec with the same parallel layout as the control vector.

For the control-forcing portion, the weight is `(quadrature weight at grid point) × (timestep dt)`. For the per-segment IC portion, it's just the spatial quadrature weight. The exact recipe matches ver3's `globalNormFiles` math; ver4 just packages it as a PETSc Vec instead of a `.dat` file.

This setup is run once per case, not per iteration.

## 6. Async-evaluation cases

Three physically distinct regimes, distinguished by what fits in one SLURM allocation:

| Case | Physical SLURM constraint | Python lifetime | `tao.solve()` viable? | Status in ver4 |
|---|---|---|---|---|
| **A-long** | one job ≥ many (msforward + msadjoint) iterations | one Python process spans many outer iters | ✓ Yes | **Built first (primary)** |
| **A** | one job = exactly one (msforward + msadjoint) | one outer iteration per Python invocation; bash resubmits between | No (cleanly) | **Future extension** |
| **B** | one job = one of {msforward, msadjoint}, never both | one *half*-iteration per Python invocation | No | **Future extension (degraded)** |

Production scale of `msforward`+`msadjoint` is unknown until they are implemented. Building A-long first is justified by:
- Shortest implementation path — TAO owns line search, two-loop recursion, history, and convergence test.
- Naturally collapses to Case A behavior if the SLURM allocation only fits one outer iteration before wall-clock — `tao.solve()` simply completes one iteration and the signal handler exits. The next allocation restarts `tao.solve()` from the saved state.
- The **only** thing that's lost in the "collapsed" mode is TAO's in-memory line-search bracket state mid-iteration. If one (msforward + msadjoint) eval already ≥ one full SLURM allocation, line search becomes a problem (any backtracking step would need its own allocation), and at that point Case A or Case B with a manual outer loop is the right escape hatch.

### Case A-long (primary) — TAO drives the loop inside one Python process

```python
# optim.py — single SLURM allocation; TAO does the iteration
import petsc4py; petsc4py.init(sys.argv)
from petsc4py import PETSc
import signal, subprocess

# 1. Load metric and persisted state
D_inv = PETSc.Vec().load(viewer('D_inv_diag.petsc'))
y     = PETSc.Vec().load(viewer('y.petsc'))
M_lmvm = PETSc.Mat().createLMVM(...)
if path.exists('history.petsc'):
    M_lmvm.load(viewer('history.petsc'))

# 2. Define the objective/gradient callback (lives entirely in y-space)
def fg(tao, y_in, g_y_out):
    x = y_in.duplicate(); x.pointwiseMult(y_in, D_inv)
    x.view(viewer('x.petsc'))
    subprocess.run(['mpirun', '-n', str(N), './msforward', ...], check=True)
    subprocess.run(['mpirun', '-n', str(N), './msadjoint', ...], check=True)
    raw_grad = PETSc.Vec().load(viewer('raw_grad.petsc'))
    g_y_out.pointwiseMult(raw_grad, D_inv)        # ∇_y h = D⁻ᵀ · dJ/dx
    return float(open('J.txt').read())

# 3. Configure TAO
tao = PETSc.TAO().create()
tao.setType('lmvm')
tao.setLMVMMatrix(M_lmvm)                          # seed history
tao.setObjectiveGradient(fg)
tao.setTolerances(grtol=cfg.tolerance)
tao.setMaximumIterations(cfg.max_iterations)

# 4. Wall-clock signal handler: snapshot and exit cleanly
def on_walltime(*_):
    y_now = tao.getSolution()
    y_now.view(viewer('y.petsc'))
    tao.getLMVMMatrix().view(viewer('history.petsc'))
    sys.exit(EXIT_RESUME)
signal.signal(signal.SIGUSR1, on_walltime)        # surrounding bash sends SIGUSR1 near wall-clock

# 5. Run
tao.solve(y)

# 6. Save final state and exit with convergence code
y.view(viewer('y.petsc'))
tao.getLMVMMatrix().view(viewer('history.petsc'))
sys.exit(EXIT_CONVERGED if tao.getConvergedReason() > 0 else EXIT_RESUME)
```

The surrounding bash wrapper:

```bash
#SBATCH --time=24:00:00
mpirun -n $SLURM_NTASKS python3 optim.py optim.yml &
PYPID=$!
# Send SIGUSR1 to Python `wall_clock_margin_seconds` before SLURM kills us
( sleep $((SLURM_TIME_SECONDS - WALL_MARGIN)) && kill -SIGUSR1 $PYPID ) &
wait $PYPID
status=$?
if [ $status -eq $EXIT_RESUME ]; then
    sbatch OPT.sh -l depend=$SLURM_JOBID
fi
```

### Case A (future extension) — one outer iteration per Python invocation

If `msforward` + `msadjoint` turns out to be so expensive that even one (forward + adjoint) eval consumes the entire SLURM allocation, then `tao.solve()` cannot drive the loop — TAO would die mid-iteration with no way to resume the line search. The escape hatch is to drop `tao.solve()` and use raw `MatLMVM` + a manually-owned outer loop; ~150 lines of Python in addition to the A-long driver. The architecture would mirror ver3's `--mode setup / schedule / log_result` triad: one Python invocation does one outer-iteration step (compute L-BFGS direction from history, write a trial `y`, emit a command file for the bash loop to run msforward+msadjoint, exit). State persists via `MatView` / `Vec.view`.

To support this the Python driver would gain a Case-A code path that bypasses `tao.solve()` and operates `MatLMVM` directly. The callback / variable-transformation / I/O contracts to the Fortran exes are unchanged — only the outer loop differs.

**Build only if production hits this constraint.** The decision point is Phase 2 (when msforward/msadjoint exist and we can measure per-iteration wall-clock).

### Case B (future extension, degraded) — msforward and msadjoint in separate jobs

The most-constrained regime: forward and adjoint cannot share a SLURM allocation. Line search becomes structurally hard — Wolfe conditions need both J and ∇J at trial points, and those would come from separate allocations.

Three workable patterns if this ever becomes necessary:

1. **Barzilai–Borwein step** (simplest): `α = ⟨s, s⟩ / ⟨s, y⟩` from the previous iteration's curvature. No line search; each outer iteration = exactly two allocations.
2. **Backtracking Armijo, eval-by-eval**: each backtrack step is a new forward allocation. Robust but pessimistically 3–8 forward allocations per outer iteration.
3. **Trust region** (TAO `BNTR`/`BQNK`): no line search needed, but requires Hessian-vector products → second-order adjoint in magudi → significant extra Fortran. Out of scope.

If Case B ever materializes, the Python driver gains a third code path on top of the Case A manual outer loop, with the line-search-step computation replaced by BB or Armijo backtracking-state-machine logic.

## 7. Line search

For **Case A-long (primary)**, line search is entirely TAO's responsibility. Default is `TAOLINESEARCHMT` (Moré–Thuente, satisfies strong Wolfe), with `TAOLINESEARCHARMIJO` and `TAOLINESEARCHGPCG` available via the `line_search.type` config knob. TAO's defaults are well-tuned for L-BFGS; no Python-side line-search code is needed.

For **Case A** (future, if built), the line search would move into Python:
- First iteration: `step = initial_step` (configured in `optim.yml`, default 1.0).
- Subsequent iterations with curvature info: `step = ⟨s, y⟩ / ⟨y, y⟩` (BB-style initial guess) clipped to `[1e-4, 1e4]`.
- Optional Armijo sufficient-decrease check between iterations: if `J_new > J_old + c1 · step · ⟨g_old, direction⟩`, halve the step and re-emit a new trial `y` without consuming an L-BFGS update. State persistence handles the multi-allocation backtracking.

For **Case B** (future, degraded), the same wrapper as Case A but with the Armijo branch promoted to default — or use BB step (no line search at all).

## 8. Configuration (`optim.yml`)

A first cut at the YAML schema:

```yaml
global_prefix: OneDWave

magudi:
  binary_directory: ../bin
  root_files:
    grid_file: OneDWave.xyz
    boundary_condition_file: bc.dat
    control_mollifier_file: OneDWave.control_mollifier.f
    target_mollifier_file: OneDWave.target_mollifier.f

control_space:
  number_of_actuators: 1
  actuator_list: ['controlRegion']
  number_of_segments: 6
  segment_length: 2400
  start_timestep: 30000
  initial_condition_controllability: 1.0
  matching_condition_weight: 1.6e-6
  periodic_solution: false

objective:
  include_time_integral: true
  include_terminal: false
  penalty_norm:
    type: base                      # or 'huber'

optimizer:
  algorithm: lbfgs                  # 'lbfgs' (TAO LMVM) — future: 'lsr1', 'bncg', 'bb'
  history_size: 3                   # m for L-BFGS / L-SR1
  tolerance: 1.0e-8
  max_iterations: 1000

line_search:
  type: mt                          # TAO line search type:
                                    #   'mt'     — Moré–Thuente (default, strong Wolfe)
                                    #   'armijo' — Armijo backtracking
                                    #   'gpcg'   — Gradient-projection conjugate-gradient

execution:
  case: A-long                      # 'A-long' (built first), future: 'A', 'B'
  wall_clock_margin_seconds: 600    # SIGUSR1 sent this many seconds before SLURM time limit
```

The schema mirrors ver3 where it can (so existing configs port naturally), drops the ver3-specific knobs that ver4 doesn't need (`command_scriptor.type`, `resource_distribution.jobs` — replaced by the simpler `mpirun -n $SLURM_NTASKS`), and adds the `optimizer`, `line_search`, and `execution` sections that ver3 hardcoded.

## 9. Comparison with ver3

| Aspect | ver3 | ver4 (Case A-long, primary build) |
|---|---|---|
| Optimizer | Hand-written Polak–Ribière NCG + golden-section / parabolic line search | TAO `lmvm` (`MATLMVMBFGS`, m=3) + Moré–Thuente line search; swap-in for L-SR1, BNCG via config |
| Outer loop | Hand-written state machine over `(hyperStep, cgStep, lineStep)` | `tao.solve()` |
| Algebra | Fortran exes called from bash (zaxpy, qfile_zaxpy, spatial_inner_product, …) | PETSc Vec/Mat methods inside one Python process |
| Forward eval | Bash launches `Nsplit` parallel `./forward` calls + glue (`./qfile_zaxpy`, `./spatial_inner_product`) | One `mpirun ./msforward` |
| Adjoint eval | Bash launches `Nsplit` parallel `./adjoint` calls + matching-condition setup + IC-gradient assembly | One `mpirun ./msadjoint` |
| Working dirs | `x0/ a/ b/ c/ x/` × `0/ 1/ 2/ …` + `out/ diff/ grad/ cg/ previous/ txt/ lagrangian/ baseline/ diffLog/ linminLog/` | Flat: just root files + `snapshots/` + `out/` |
| Per-segment `magudi.inp` mutation | Required (`setOption.sh` sed-edits per segment) | Eliminated — single `magudi.inp` |
| Control-space norm | Computed via `./spatial_inner_product` + `*.norm_*.dat` files | Variable transformation in Python with `D_diag.petsc` |
| Resume model | Python re-instantiates from HDF5 every `--mode schedule` call | Python stays alive across many TAO iterations; SIGUSR1 signal handler at wall-clock checkpoints `y.petsc` + `history.petsc`; bash resubmits |
| Augmented Lagrangian | Outer `hyperStep` loop, `lagrangian/*.q` files | Penalty baked into J returned by msforward; TAOALMM available as upgrade if needed |
| `command_scriptor` (bash/srun/flux) | Three subclasses generating launcher syntax | Gone — `mpirun -n $SLURM_NTASKS` in one bash file |
| Lines of Python | ~1500 | Estimated ~250 (Case A-long alone); +150 each for future Case A and Case B paths |
| New Fortran code | None | `bin/msforward.f90`, `bin/msadjoint.f90`, `bin/compute_norm.f90` |

## 10. Rollout phases

Multi-month effort. Each phase produces a working artifact.

### Phase 1 — strawman validation with `tao.solve()` *(1–2 weeks)*

- Write a minimal `optim.py` for **Case A-long**: `tao.solve()` driving `lmvm` (m=3) inside one Python process. Callback wraps the **existing** per-segment `forward` + `adjoint` exes via subprocess + bash command file. No new Fortran yet.
- Use the existing `*.norm_*.dat` files to derive `M_diag.petsc` (and `D_inv_diag.petsc`) for the variable transformation; transformation lives in the callback.
- Wire the SIGUSR1 wall-clock handler that extracts `y.petsc` + `history.petsc` from TAO at allocation boundary.
- Run on the OneDWave CI example; verify convergence vs. ver3.

**Goal**: prove that TAO + the variable transformation can drive an optimization to convergence on a real magudi case, **before** investing in the Fortran rewrite.

### Phase 2 — `msforward` / `msadjoint` prototypes *(3–4 weeks)*

- Implement `bin/msforward.f90` and `bin/msadjoint.f90`. Each drives a `t_Solver` over `Nsplit` segments concatenated.
- Match ver3's J and gradient values to round-off on OneDWave (paired test). The adjoint validation should reuse the `gradient_accuracy` finite-difference machinery currently in `bin/CheckGradientAccuracy.f90`.
- Add a CI test analogous to `.github/workflows/optim_grad_test.sh` that exercises msforward/msadjoint directly.
- **Decision point**: measure wall-clock per (msforward + msadjoint) at production scale. If many iterations comfortably fit in one SLURM allocation → stay in A-long. If exactly one fits → schedule Phase 4b (Case A code path). If only forward-or-adjoint fits → schedule Phase 4c (Case B code path).

### Phase 3 — PETSc Vec wiring *(1–2 weeks)*

- Replace `.dat` and `.q` I/O between Python and msforward with PETSc binary format throughout. msforward reads `x.petsc`; msadjoint writes `raw_grad.petsc`.
- This is the change that makes 100 GB scale workable: no driver-side numpy buffering, all data flows through distributed PETSc Vecs.
- Implement `bin/compute_norm.f90` to write `M_diag.petsc` at setup.

### Phase 4a — A-long polish *(1–2 weeks)*

- Choose / tune TAO line search (`mt` default, expose `armijo` / `gpcg` via config).
- HDF5 history log matching ver3's diagnostics (`forwardLog`, `gradientLog`); per-iteration callback into TAO via `TaoSetMonitor`.
- Journal logging, signal-handler robustness, exit-code conventions for the bash wrapper.

### Phase 4b — Case A extension *(only if Phase 2 decision triggers, ~1–2 weeks)*

- Add a Case-A code path that bypasses `tao.solve()`: raw `MatLMVM` + manually-owned outer loop, one outer iteration per Python invocation.
- Reuse the callback / variable-transformation / I/O contracts from A-long verbatim — only the outer loop differs.
- Driven by `execution.case: A` in `optim.yml`.

### Phase 4c — Case B extension *(only if Phase 2 decision triggers, ~2–3 weeks)*

- Add a Case-B code path on top of Case A's manual outer loop, with the line-search step replaced by BB or Armijo backtracking.
- Driven by `execution.case: B` in `optim.yml`.

### Phase 5 — CI port + ver3 deprecation *(1 week)*

- Port `optim_grad_test.sh` and `optim_test.sh` to ver4 (Case A-long for the CI environment, since runs are short).
- Deprecate ver3 in CLAUDE.md and `utils/legacy/optimization_ver3/`'s top-level docstrings; do **not** delete — keep as fallback for in-flight production runs.

**Total**: ~2 months focused work for the A-long-only path (Phases 1–4a + 5). Phases 4b/4c add 1–3 weeks each if production hits those constraints.

## 11. Open questions for later

These don't need resolution before Phase 1 but should be revisited during the rollout:

- **Which physical regime production hits**: A-long, A, or B. Resolved at the Phase 2 decision point once `msforward` + `msadjoint` exist and we can measure per-iteration wall-clock. Determines whether Phases 4b / 4c are needed.
- **Disk-spilling of `MatLMVM` history**: at m = 3, in-memory is fine. If a future case has even larger control vectors and m needs to grow, we may need to spill history vectors to PETSc binary files between iterations and `MatLoad` them on demand inside the two-loop. Measure first.
- **L-SR1 vs L-BFGS comparison**: same memory footprint, often converges faster on non-convex objectives. Worth running side-by-side once Phase 4a is stable.
- **Preconditioning beyond the variable transformation**: a diagonal Gauss–Newton Hessian at the initial point, or a low-frequency-mode preconditioner, can give 5–10× iteration-count reduction via `TaoLMVMSetH0`. Defer until baseline performance is characterized.
- **Augmented Lagrangian (TAOALMM)**: only revisit if penalty-method convergence stalls on production cases. The penalty-baked-into-J approach matches ver3 functionally and is simpler.
- **Container changes**: PETSc + petsc4py addition to [docker/Dockerfile](../../docker/Dockerfile). Need to coordinate the `arm64` rebuild that is currently maintained out-of-band.
