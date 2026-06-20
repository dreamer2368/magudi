# `utils/legacy/`

Quarantined Python utilities and optimization frameworks that have been superseded by the installable [`magudi_utils`](../magudi_utils/) package (which now bundles `magudi-optim` / `magudi-msgrad` and the PLOT3D / mesh / FWH helpers). Nothing here is on the active CI path; the directories are kept on disk so in-flight production runs and historical experiments can still be reproduced.

Pull from `utils/magudi_utils/` for any new work.

## Index

- [`inexactnewton/`](inexactnewton/) — original Newton-based optimization driver (`newton.py`, `inexactNewton.py`, `filenames_newton.py`). Predates the `optimization/` generations.

- [`optimization/`](optimization/) — generation 1: shell-orchestrated outer loop with a Python Fletcher-Reeves Polak-Ribiere driver (`frprmn.py`) plus per-step shell scripts (`forward.sh`, `adjoint.sh`, `ZAXPY.sh`, `OPT.sh`, …).

- [`optimization_ver2/`](optimization_ver2/) — generation 2: PBS / flux-cluster-style orchestration (`OPT.flux`, `OPT.sh`, `OPT-init.sh`). Adds the `base_extension.py` / `command_scriptor.py` abstractions that gen 3 inherited.

- [`optimization_ver3/`](optimization_ver3/) — generation 3: Python framework with `optimization.py` / `optimizer.py` driving per-segment shell calls. See [`optimization_ver3/DIRECTORY_STRUCTURE.md`](optimization_ver3/DIRECTORY_STRUCTURE.md) for the on-disk layout. Replaced by `magudi_utils`'s `magudi-optim` (multi-segment Fortran binaries `msforward`/`msadjoint` + a single TAO L-BFGS outer loop), but kept here because some production studies still drive ver3 directly.

- [`python/`](python/) — residue of the old flat `utils/python/` directory after the helpers (`plot3dnasa`, `PLOT3D`, `RoundJet`, `SingleBlockCartesian`, `SummationByParts`, `fwhsolver`, `matplotlibhelper`) moved into `magudi_utils`:
  - `checkGradientAccuracy_linearized.py` — one-off FD validator paired with the `linearized_relation` test setup.
  - `deprecate-py2/` — Python-2 originals of the helpers, kept for archeological reference; the Python-3 ports are now in `magudi_utils`.
