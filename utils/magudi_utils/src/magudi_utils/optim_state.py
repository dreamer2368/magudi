"""Scalar state machine for the manual L-BFGS driver.

The TAO driver (tao_main) keeps its outer-loop state in TAO's in-memory
structures and only writes y.petsc + history.petsc at allocation boundaries
(after tao.solve() returns). The manual driver runs one outer-loop step per
Python invocation, so it has to persist *everything* needed to resume between
allocations: the iterate (y.petsc), the L-BFGS history (history.petsc), the
accepted gradient (g.petsc), the search direction (d.petsc), and a small
scalar state machine describing what to do next.

This module owns the small scalar piece -- written to lbfgs_state.json as an
atomic rename. The PETSc Vecs are handled by ParallelIOHandler.

Format:
    {
      "schema_version": 1,
      "iter": <int>,
      "phase": "fresh" | "ls_trial" | "converged",
      "f_baseline": <float>,        # J at the current accepted iterate
      "gT_d_baseline": <float>,     # <g_accepted, d>; cached descent slope
      "converged_reason": null | <str>,
      "line_search": {              # LS-subclass-owned nested blob
        "type": "golden_section",
        ...
      }
    }

Why JSON instead of HDF5 attrs on the existing optim.h5 log: bash wrappers can
inspect "phase" with `jq -r` without spinning up Python+h5py; os.replace gives
POSIX-atomic write semantics that HDF5 cannot match; and JSON failure does not
brick the historical log.
"""
import json
import os
import tempfile
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Optional


SCHEMA_VERSION = 1


# Phase tokens. Read by the manual driver and (optionally) by the bash
# wrapper via `jq -r .phase`. Keep these as constants so a typo is caught
# at import time rather than as a silent stale-state mismatch.
PHASE_FRESH = "fresh"            # cold start: bootstrap iter-0 fwd+adj at baseline
PHASE_LS_TRIAL = "ls_trial"      # one forward at LS pending_alpha
PHASE_LS_ACCEPTED = "ls_accepted"  # one forward + one adjoint at LS-accepted alpha_star
PHASE_CONVERGED = "converged"


# Exit codes. The bash wrapper reads these to decide resubmission. Kept here
# because manual_main and the CI driver both depend on this contract.
EXIT_CONVERGED = 0
EXIT_RESUME = 42


@dataclass
class OptimState:
    """In-memory representation of lbfgs_state.json."""

    iter: int = 0
    phase: str = PHASE_FRESH
    f_baseline: float = 0.0
    gT_d_baseline: float = 0.0
    # Only meaningful when phase == ls_accepted: the LS-chosen step length
    # at which the next allocation must re-run msforward (to seed snapshots)
    # before msadjoint. The just-finished ls_trial's forward was at a
    # different alpha (the new b may not be the most recent trial).
    accepted_alpha: float = 0.0
    converged_reason: Optional[str] = None
    line_search: Dict[str, Any] = field(default_factory=dict)
    schema_version: int = SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "OptimState":
        version = d.get("schema_version", 1)
        if version != SCHEMA_VERSION:
            raise ValueError(
                f"lbfgs_state.json schema_version {version} != "
                f"current {SCHEMA_VERSION}"
            )
        return cls(
            iter=int(d["iter"]),
            phase=str(d["phase"]),
            f_baseline=float(d["f_baseline"]),
            gT_d_baseline=float(d["gT_d_baseline"]),
            accepted_alpha=float(d.get("accepted_alpha", 0.0)),
            converged_reason=d.get("converged_reason"),
            line_search=dict(d.get("line_search", {})),
            schema_version=version,
        )


def read_state(path: str) -> Optional[OptimState]:
    """Return the parsed OptimState, or None if `path` does not exist.

    A missing file is the "first invocation, fresh start" signal; the caller
    transitions through phase == "fresh".
    """
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        payload = json.load(fh)
    return OptimState.from_dict(payload)


def write_state(state: OptimState, path: str) -> None:
    """Atomically rewrite `path` with the serialized state.

    Writes to <path>.tmp on the same directory, fsyncs, then os.replace --
    POSIX guarantees that a reader observes either the prior version or the
    complete new version. Partial writes on crash leave the prior file intact.

    Must be called from rank 0 only (the JSON is rank-invariant scalars).
    """
    dirname = os.path.dirname(os.path.abspath(path)) or "."
    fd, tmp = tempfile.mkstemp(dir=dirname, prefix=".lbfgs_state.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(state.to_dict(), fh, indent=2, sort_keys=True)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise
