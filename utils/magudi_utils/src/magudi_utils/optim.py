"""magudi_utils MPI-parallel optimization entry point.

Thin dispatcher: argparse + config load + MPI init, then branch on
optimization.type (read from the YAML, default "tao") to either:

    "tao"     -> optimizers.tao_main      (tao.solve() owns the outer loop;
                                           targets DESIGN.md "Case A-long"
                                           where many outer iterations fit
                                           per SLURM allocation)

    "manual"  -> optimizers.manual_main   (manual L-BFGS state machine, one
                                           step per Python invocation;
                                           targets PLAN_case_a.md "Case A"
                                           where one allocation barely fits
                                           one msforward + one msadjoint)

Run:
    mpirun -n N_petsc magudi-optim optim.yml [--max-iter K] [--mode base|slurm]

The two paths share the same YAML schema, M_diag / D / D_inv construction,
ParallelIOHandler layout, and msforward / msadjoint subprocess interface;
the only YAML keys that change between them are `optimization.type` and
(optionally) `optimization.line_search.type`.
"""
import argparse
import sys

# mpi4py MUST import before petsc4py so PETSc binds to MPI.COMM_WORLD.
from mpi4py import MPI  # noqa: F401  (import-order side effect)
from petsc4py import PETSc

from .inputs import InputParser
from .optimizers import manual_main, tao_main


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="Path to optim.yml")
    parser.add_argument(
        "--max-iter", type=int, default=5,
        help="iterations this run (added to loaded iter if resuming)"
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="print extra progress information"
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="print debug-level diagnostics"
    )
    parser.add_argument(
        "--mode", choices=("base", "slurm"), default="base",
        help="Worker launch mode. base: mpirun (default). "
             "slurm: srun --overlap (production SLURM)."
    )
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    pcomm = PETSc.COMM_WORLD

    cfg = InputParser(args.config)
    opt_type = cfg.getInput(
        ["optimization", "type"], fallback="tao")

    if opt_type == "tao":
        return tao_main(args, cfg, comm, pcomm)
    if opt_type == "manual":
        return manual_main(args, cfg, comm, pcomm)
    raise SystemExit(
        f"optimization.type = {opt_type!r} not recognized; "
        f"expected 'tao' or 'manual'."
    )


if __name__ == "__main__":
    sys.exit(main())
