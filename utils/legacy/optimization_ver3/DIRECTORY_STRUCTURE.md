# Optimization run directory structure

This is the on-disk layout an `optimization_ver3` run uses, for an example case with `global_prefix: OneDWave`, `Nsplit: 3`, and `actuator_list: ['controlRegion']`. (The OneDWave example in CI uses `Nsplit: 6`; layout is identical, just more segment subdirectories and more `*-k.*` files.)

`Nsplit = 3` ⇒ segments `k ∈ {0, 1, 2}` ⇒ each working dir contains subdirectories `0/`, `1/`, `2/`. `NcontrolRegion = 1` ⇒ one actuator named `controlRegion`. With `start_timestep: 30000` and `segment_length (Nts): 2400`, segment-start timesteps are `30000, 32400, 34800` and segment-end timesteps are `32400, 34800, 37200`.

Code paths cited below: directory list comes from [filenames.py:108-110](filenames.py#L108-L110); per-segment file naming from [filenames.py:170-192](filenames.py#L170-L192); the working-tree creation loop from [optimizer.py:331-346](optimizer.py#L331-L346).

## Tree

```
<run-root>/                                         # cwd when you launch optimization.py
│
├── magudi.inp                                      # global input — never run from, only copied
├── bc.dat
├── OneDWave.xyz                                    # PLOT3D grid (config.getInput magudi.root_files)
├── OneDWave.control_mollifier.f
├── OneDWave.target_mollifier.f
│
├── optim.yml                                       # YAML config you passed on the CLI
│
│   ── binary symlinks (created by `--mode setup`) ──
├── forward             -> ../bin/forward
├── adjoint             -> ../bin/adjoint
├── zaxpy               -> ../bin/zaxpy
├── zxdoty              -> ../bin/zxdoty
├── zwxmwy              -> ../bin/zwxmwy
├── qfile_zaxpy         -> ../bin/qfile_zaxpy
├── spatial_inner_product -> ../bin/spatial_inner_product
├── control_space_norm  -> ../bin/control_space_norm
├── slice_control_forcing -> ../bin/slice_control_forcing
├── paste_control_forcing -> ../bin/paste_control_forcing
├── patchup_qfile       -> ../bin/patchup_qfile
│
│   ── per-iteration scratch files (regenerated each `--mode schedule`) ──
├── OneDWave.command.sh                             # the next batch of shell commands
├── setOption.sh                                    # generated helper (sed-edits magudi-<k>.inp)
│
│   ── persisted optimizer state ──
├── OneDWave.optim.h5                               # HDF5 log: hyper/cg/line counters,
│                                                   #   forward+gradient history per CG step,
│                                                   #   line-min bracket per CG step,
│                                                   #   penalty weights per hyper step
├── OneDWave.journal.txt                            # python logging output (timestamped)
├── OneDWave.line_minimization.txt                  # current line-search bracket as TSV
│                                                   #   (step, QoI, directory_index ∈ {a,b,c,x,0})
│
│   ── global-level control-space & gradient files ──
├── OneDWave.norm_controlRegion.dat                 # global norm matrix for inner-product
│                                                   #   (must exist before `--mode setup` — checked
│                                                   #    against `prerequisites` at optimizer.py:315)
│
│   ── working trees: one per line-search role ──
│   ── (all 5 created by `--mode setup`, populated during the loop) ──
├── x0/                                             # current iterate (best-known control)
│   ├── OneDWave-0.ic.q                             # initial condition for segment 0
│   ├── OneDWave-1.ic.q                             # initial condition for segment 1
│   ├── OneDWave-2.ic.q                             # initial condition for segment 2
│   ├── OneDWave.control_forcing_controlRegion.dat  # current control forcing (concatenated over time)
│   ├── 0/
│   │   ├── forward                  -> ../../../bin/forward
│   │   ├── adjoint                  -> ../../../bin/adjoint
│   │   ├── magudi-0.inp                            # per-segment input (output_prefix=OneDWave-0,
│   │   │                                           #   initial_condition_file=../../OneDWave-0.ic.q,
│   │   │                                           #   number_of_timesteps=2400, save_interval=2400,
│   │   │                                           #   controller_switch=true/false)
│   │   ├── OneDWave-0-00030000.q                   # forward solution snapshots (start of segment 0)
│   │   ├── OneDWave-0-00032400.q                   # …and end (= matchingForwardFile for k=0)
│   │   ├── OneDWave-0-00032400.adjoint.q           # adjoint snapshot (matchingAdjointFile for k=0)
│   │   ├── OneDWave-0.control_forcing.controlRegion.dat   # slice of global control forcing
│   │   ├── OneDWave-0.gradient.controlRegion.dat          # slice of gradient (moved to grad/ later)
│   │   ├── OneDWave-0.cost_functional.txt          # forward QoI history
│   │   ├── OneDWave-0.cost_sensitivity.txt         # adjoint sensitivity history
│   │   └── OneDWave-0.forward_run.txt              # final forward QoI scalar
│   ├── 1/
│   │   ├── forward                  -> ../../../bin/forward
│   │   ├── adjoint                  -> ../../../bin/adjoint
│   │   ├── magudi-1.inp
│   │   ├── OneDWave-1-00032400.q
│   │   ├── OneDWave-1-00034800.q                   # matchingForwardFile for k=1
│   │   ├── OneDWave-1-00032400.adjoint.q           # icAdjointFile for k=1
│   │   ├── OneDWave-1-00034800.adjoint.q           # matchingAdjointFile for k=1
│   │   ├── OneDWave-1.ic.adjoint.q                 # icGradientFile for k=1 (moved to grad/)
│   │   └── … (cost_functional.txt, cost_sensitivity.txt, forward_run.txt, etc.)
│   └── 2/                                          # same shape as 1/
│
├── a/                                              # left edge of line-search bracket
│   ├── OneDWave-0.ic.q                             # snapshot of x0's ICs at start of line search
│   ├── OneDWave-1.ic.q
│   ├── OneDWave-2.ic.q
│   ├── OneDWave.control_forcing_controlRegion.dat
│   ├── 0/   ← same per-segment shape as x0/0/
│   ├── 1/
│   └── 2/
├── b/                                              # interior of bracket (current best step)
│   └── … (same shape as a/)
├── c/                                              # right edge of bracket
│   └── … (same shape as a/)
├── x/                                              # next trial step (overwritten each linmin call)
│   └── … (same shape as a/, purged after each evaluation by purgeDirectoryCommand)
│
│   ── auxiliary directories (created by `--mode setup`, dirList at filenames.py:108) ──
├── out/                                            # captured stdout/stderr from each MPI job
│   ├── forward_result_0.out                        #   one per (procedure, job-index) pair
│   ├── forward_result_1.out
│   ├── adjoint_result_0.out
│   ├── zaxpy_result_0.out
│   ├── qfile-zaxpy_result_0.out
│   └── …
│
├── diff/                                           # state-mismatch q-files (segment-boundary jumps)
│   ├── OneDWave-0.diff.q
│   ├── OneDWave-1.diff.q
│   └── OneDWave-2.diff.q
│
├── grad/                                           # gathered gradient files
│   ├── OneDWave.gradient_controlRegion.dat         # global control-forcing gradient
│   ├── OneDWave-0.ic.adjoint.q                     # IC-gradient for segment 0
│   ├── OneDWave-1.ic.adjoint.q
│   └── OneDWave-2.ic.adjoint.q
│
├── cg/                                             # conjugate-gradient direction (Polak-Ribiere)
│   ├── OneDWave.conjugate_gradient_controlRegion.dat
│   ├── OneDWave-0.ic.conjugate_gradient.q
│   ├── OneDWave-1.ic.conjugate_gradient.q
│   └── OneDWave-2.ic.conjugate_gradient.q
│
├── previous/                                       # last CG step's gradient & direction
│   ├── previous.OneDWave.gradient_controlRegion.dat
│   ├── previous.OneDWave-0.ic.adjoint.q
│   ├── previous.OneDWave.conjugate_gradient_controlRegion.dat
│   └── previous.OneDWave-0.ic.conjugate_gradient.q
│
├── txt/                                            # all scalar outputs of inner-product / norm jobs
│   ├── OneDWave.inner_product_controlRegion.txt
│   ├── OneDWave.gg_controlRegion.txt               # ⟨grad, grad⟩ for control forcing
│   ├── OneDWave.dgg_controlRegion.txt              # ⟨Δgrad, grad⟩ for Polak-Ribiere
│   ├── OneDWave-0.inner_product.txt                # …same trio for IC-control space, per segment
│   ├── OneDWave-0.gg.txt
│   ├── OneDWave-0.dgg.txt
│   ├── OneDWave-0.diff.txt                         # ‖state-mismatch‖² per segment
│   ├── OneDWave-0.lagrangian.txt                   # ⟨λ, mismatch⟩ when augmented_lagrangian=true
│   ├── OneDWave.terminal_objective.txt             # terminal QoI when objective.include_terminal=true
│   ├── OneDWave.adjoint_run_controlRegion.txt
│   └── …
│
├── lagrangian/                                     # only populated when penalty_norm.augmented_lagrangian=true
│   ├── OneDWave-0.lagrangian_multiplier.q
│   ├── OneDWave-1.lagrangian_multiplier.q
│   └── OneDWave-2.lagrangian_multiplier.q
│
├── baseline/                                       # untouched baseline solutions for state_log diffs
│   ├── OneDWave-00030000.q
│   ├── OneDWave-00032400.q
│   └── OneDWave-00034800.q
│
├── diffLog/                                        # archived diff.q snapshots (saveDiffFiles=true)
│   ├── 0/                                          #   one subdir per CG step
│   │   ├── OneDWave-2.diff.0.q
│   │   └── OneDWave-1.diff.0.q
│   └── 1/
│
└── linminLog/                                      # archived line-min logs (created in dirList,
                                                    #   contents written via the lineMinLog HDF5 group
                                                    #   inside OneDWave.optim.h5; this dir is reserved
                                                    #   for any external dump/script use)
```

## How the working trees are used

```
            ┌── x0 ── current iterate (best-known control)
            │
            │   evaluate forward at x0 ─────┐
            │                               │ initialAdjoint() → grad → CG direction
            │                               ▼
            ├── a ── snapshot of x0 ICs at start of line search   (step = 0,    QoI = J(x0))
            ├── b ── interior of bracket   (step = step*,         QoI = J(x0 - step* · CG))
            ├── c ── right edge of bracket (step > b's step)
            └── x ── next trial            (overwritten each evaluation,
                                             purged by purgeDirectoryCommand after readout)

            after LINMIN converges:  b/ → x0/, a/ and c/ purged, gradient → previous/
```

## File-name rules at a glance

| Pattern                                          | Meaning                                                     |
|--------------------------------------------------|-------------------------------------------------------------|
| `<prefix>-<k>.ic.q`                              | Initial condition for segment k                             |
| `<prefix>-<k>-<timestep>.q`                      | Forward snapshot at given timestep within segment k         |
| `<prefix>-<k>-<timestep>.adjoint.q`              | Adjoint snapshot                                            |
| `<prefix>-<k>.diff.q`                            | State-mismatch (forward-end of segment k-1 minus IC of k)   |
| `<prefix>-<k>.ic.adjoint.q`                      | IC-gradient for segment k                                   |
| `<prefix>.gradient_<actuator>.dat`               | Global control-forcing gradient for an actuator             |
| `<prefix>.conjugate_gradient_<actuator>.dat`     | CG search direction for an actuator                         |
| `<prefix>-<k>.ic.conjugate_gradient.q`           | CG search direction in IC-space for segment k               |
| `previous.*`                                     | Snapshot at the previous CG step (in `previous/`)           |

## Notes

- All initial conditions and control-forcing files are kept at the **working-tree root** (`x0/`, `a/`, `b/`, `c/`, `x/`), not inside the per-segment subdirectories. The per-segment `magudi-<k>.inp` references them via relative paths like `../../OneDWave-<k>.ic.q`.
- The `_extension` (`time_splitting.use_state_mollifier=true`) variant in [base_extension.py](base_extension.py) keeps the same on-disk layout but runs segments **sequentially** (not in parallel), and inserts a `./patchup_qfile` call before each forward run to mollify the IC at the segment boundary.
- `out/` is the only directory that grows monotonically across iterations — everything else is overwritten in place. If you're tight on disk, periodically purge `out/*.out`.
- The `<run-root>` is wherever you `cd` to before invoking `python3 optimization.py optim.yml --mode setup`. The CI wires this up by copying `examples/OneDWave/*` and `utils/optimization_ver3/*.py` into a fresh build subdirectory; see [.github/workflows/optim_grad_test.sh](../../.github/workflows/optim_grad_test.sh).
