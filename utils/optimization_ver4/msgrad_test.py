"""Multi-segment gradient accuracy test for msforward / msadjoint (Phase 3).

Run as
    mpirun -n N_petsc python3 msgrad_test.py msgrad.yml [--mode full|ctrl|ic]

from the staged OneDWave_msgrad/ directory. The driver
  - parses msgrad.yml via InputParser (typed YAML reader from inputs.py);
  - spawns ./msforward, ./msadjoint, ./zaxpy, ./qfile_zaxpy via MPI_Comm_spawn
    using the resource_distribution.jobs.{forward,adjoint,petsc} block (lifted
    from optim.parallel.yml's convention);
  - runs a Taylor finite-difference test
        s(h) = ( J(x0 + h*g) - J(x0) ) / h  ->  <g,g>_M  as h -> 0
    where the reference <g,g>_M is the msadjoint .adjoint_run.txt scalar; and
  - independently cross-checks that same scalar by reading the gradient slabs
    via ParallelIOHandler and weighting them by the M_diag pieces produced by
    compute_norm. Disagreement above 1e-10 relative aborts before the FD loop
    runs (the file-aligned layout is misreading something, and the test would
    be meaningless otherwise).

state_controllability MUST be 1.0 for the <g,g>_M cross-check to hold:
compute_norm scales the IC portion of M_diag by state_controllability, but
msadjoint's per-segment IC inner product (bin/msadjoint.f90:307-318) does
not. We assert this at startup rather than silently failing the cross-check.

zaxpy and qfile_zaxpy run at n=1 -- .dat and .q files are rank-portable via
MPI-IO, and the perturbation step is O(N) memory traffic, not compute-bound.
"""
import argparse
import math
import os
import shlex
import sys

# mpi4py MUST import before petsc4py so PETSc binds to MPI.COMM_WORLD.
from mpi4py import MPI  # noqa: F401  (import-order side effect)
from petsc4py import PETSc

from inputs import InputParser
from parallel_io import ParallelIOHandler, parse_layout


GG_TOL = 1.0e-10  # relative tolerance for the <g,g>_M cross-check vs msadjoint


def _spawn_and_wait(executable, args, n):
    """Collectively spawn `n` MPI children; block until they finalize.

    Lifted verbatim from utils/optimization_ver4/optim.parallel.py. All
    N_petsc parent ranks must call this with identical arguments; only one
    set of `n` children is launched (not N_petsc * n).

    Sync is via inter-comm Barrier (true rendezvous), NOT Disconnect (which
    doesn't wait if no messages flow). The child calls a matching
    MPI_Barrier(parent) in disconnectParentIfSpawned (MPIHelperImpl.f90)
    immediately before MPI_Comm_disconnect + MPI_Finalize, so the parent's
    Barrier returns only after the child has flushed all I/O.
    """
    comm = MPI.COMM_WORLD
    log_path = os.path.abspath(
        os.path.join("out", os.path.basename(executable) + ".out")
    )
    if comm.Get_rank() == 0:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        open(log_path, "w").close()
    comm.Barrier()

    quoted_args = " ".join(shlex.quote(a) for a in args)
    redir_cmd = (
        f"exec {shlex.quote(executable)} {quoted_args} "
        # f">> {shlex.quote(log_path)} 2>&1"
    )

    inter = MPI.COMM_WORLD.Spawn(
        "/bin/bash",
        args=["-c", redir_cmd],
        maxprocs=n,
        info=MPI.INFO_NULL,
        root=0,
    )
    inter.Barrier()
    inter.Disconnect()
    MPI.COMM_WORLD.Barrier()


def _read_scalar(path, comm):
    """Rank-0 reads a single float; bcast to all ranks."""
    if comm.Get_rank() == 0:
        with open(path) as fh:
            val = float(fh.read().strip())
    else:
        val = None
    return comm.bcast(val, root=0)


def _read_sub_adjoint(path, comm):
    """Rank-0 parses .sub_adjoint_run.txt -> (sum_ctrl_IP, sum_ic_IP); bcast."""
    if comm.Get_rank() == 0:
        ctrl_total = 0.0
        ic_total = 0.0
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                ctrl_total += float(parts[1])
                ic_total += float(parts[2])
        out = (ctrl_total, ic_total)
    else:
        out = None
    return comm.bcast(out, root=0)


def _build_paths(schema, prefix, mode, ic_norm_q_path=None):
    """One path per slot for a given access mode.

    mode = 'x'      -> control vector slabs read by msforward / written by Python
                       actuator -> .control_forcing_<name>.dat
                       ic       -> -<k>.ic.q
    mode = 'grad'   -> raw_grad slabs (written by msadjoint, read by Python)
                       actuator -> .gradient_<name>.dat
                       ic       -> -<k>.ic.adjoint.q
    mode = 'metric' -> M_diag slabs (written by compute_norm)
                       actuator -> .norm_<name>.dat
                       ic       -> <ic_norm_q_path>  (shared)
    """
    paths = []
    for s in schema:
        if s.kind == "actuator":
            if mode == "x":
                paths.append(f"{prefix}.control_forcing_{s.identifier}.dat")
            elif mode == "grad":
                paths.append(f"{prefix}.gradient_{s.identifier}.dat")
            elif mode == "metric":
                paths.append(f"{prefix}.norm_{s.identifier}.dat")
            else:
                raise ValueError(f"unknown mode {mode!r}")
        else:  # ic
            k = int(s.identifier)
            if mode == "x":
                paths.append(f"{prefix}-{k}.ic.q")
            elif mode == "grad":
                paths.append(f"{prefix}-{k}.ic.adjoint.q")
            elif mode == "metric":
                if ic_norm_q_path is None:
                    raise ValueError("metric mode requires ic_norm_q_path")
                paths.append(ic_norm_q_path)
            else:
                raise ValueError(f"unknown mode {mode!r}")
    return paths


def _mask_g_by_mode(g, mode, io):
    """Zero out g-slabs of the slot kind NOT covered by `mode`.

    mode='full' -> no masking; mode='ctrl' zeros ic slabs; mode='ic' zeros
    actuator slabs. With the file-aligned layout, each rank owns at most one
    slot (or, in N_petsc==1, owns all of them contiguously) -- we can mask in
    place via the local array without any global indexing.
    """
    if mode == "full":
        return
    kind_to_zero = "ic" if mode == "ctrl" else "actuator"
    arr = g.getArray()
    if io.size == 1:
        cursor = 0
        for s in io.schema:
            if s.kind == kind_to_zero:
                arr[cursor:cursor + s.size] = 0.0
            cursor += s.size
    else:
        if kind_to_zero == "actuator" and io.is_actuator_rank:
            arr[:] = 0.0
        elif kind_to_zero == "ic" and io.is_ic_rank:
            arr[:] = 0.0
    g.assemble()


def _check_gg_against_msadjoint(M_diag, g, gg_total):
    """Independent <g,g>_M computation: g_sq.dot(M_diag) vs msadjoint's gg_total.

    Aborts via SystemExit if the two values disagree by more than GG_TOL
    relative -- a mismatch means the file-aligned layout is misreading the
    gradient files, and the downstream FD slope check would be meaningless.
    """
    g_sq = g.duplicate()
    g_sq.pointwiseMult(g, g)
    gg_py = g_sq.dot(M_diag)  # collective AllReduce
    g_sq.destroy()

    rel = abs(gg_py - gg_total) / max(abs(gg_total), 1e-30)
    PETSc.Sys.Print(
        f"<g,g>_M cross-check: Python={gg_py: .12e}  "
        f"msadjoint={gg_total: .12e}  rel_diff={rel: .3e}"
    )
    if rel > GG_TOL:
        raise SystemExit(
            f"FAIL: <g,g>_M mismatch above {GG_TOL:.0e} relative; "
            f"Python ParallelIOHandler view of the gradient does not match "
            f"msadjoint's sum. Check state_controllability, schema layout, "
            f"or .gradient/.adjoint.q write paths."
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="Path to msgrad.yml")
    parser.add_argument("--mode", choices=("full", "ctrl", "ic"), default="full",
                        help="full: perturb control+IC, reference = msadjoint <g,g> total. "
                             "ctrl: perturb only control, reference = sum control IPs. "
                             "ic:   perturb only IC (k>=1), reference = sum ic IPs.")
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    cfg = InputParser(args.config)

    prefix = cfg.getInput(["output_prefix"], datatype=str)
    h_list = cfg.getInput(["finite_difference", "step_sizes"], datatype=list)
    order_threshold = cfg.getInput(
        ["finite_difference", "order_threshold"], fallback=0.5)
    max_bad_orders = cfg.getInput(
        ["finite_difference", "max_bad_orders"], fallback=3)

    # Resource distribution: matches optim.parallel.yml's convention.
    N_forward = int(cfg.getInput(
        ["resource_distribution", "jobs", "forward"], datatype=list)[1])
    N_adjoint = int(cfg.getInput(
        ["resource_distribution", "jobs", "adjoint"], datatype=list)[1])
    N_petsc = int(cfg.getInput(
        ["resource_distribution", "jobs", "petsc"], datatype=list)[1])
    if comm.Get_size() != N_petsc:
        raise SystemExit(
            f"mpirun -n mismatch: launched with {comm.Get_size()} ranks but "
            f"resource_distribution.jobs.petsc requests {N_petsc}"
        )

    # state_controllability MUST be 1.0 -- see module docstring + compute_norm /
    # msadjoint metric mismatch note.
    state_ctrl = cfg.getInput(
        ["time_splitting", "state_controllability"], fallback=1.0)
    if abs(state_ctrl - 1.0) > 0.0:
        raise SystemExit(
            f"state_controllability must be 1.0 for the <g,g>_M cross-check "
            f"(compute_norm scales the IC norm by it, msadjoint does not); "
            f"got {state_ctrl}."
        )

    # Baseline pass: msforward then msadjoint. msforward reads the staged
    # .control_forcing.dat + .ic.q files; msadjoint writes .gradient.dat +
    # .ic.adjoint.q. After this, both x_base and g exist on disk.
    _spawn_and_wait(
        "./msforward",
        ["--input", "magudi.inp", "--output", "J0.txt"],
        N_forward,
    )
    j0 = _read_scalar("J0.txt", comm)
    _spawn_and_wait(
        "./msadjoint",
        ["--input", "magudi.inp", "--output", "gg.txt"],
        N_adjoint,
    )

    ctrl_sum, ic_sum = _read_sub_adjoint(
        f"{prefix}.sub_adjoint_run.txt", comm)
    gg_total = _read_scalar("gg.txt", comm)

    # One ParallelIOHandler drives both the cross-check and the FD perturbation
    # loop. With it we replace ./zaxpy / ./qfile_zaxpy (which don't call
    # disconnectParentIfSpawned and hang the parent's inter-comm Barrier) by
    # PETSc.Vec.axpy on the host-side and a single io.write per h.
    layout_path = f"{prefix}.layout.txt"
    ic_norm_q_path = f"{prefix}.norm_ic.q"
    schema = parse_layout(layout_path)
    io = ParallelIOHandler(schema, ic_norm_q_path, comm=comm)
    io.report_balance()

    x_paths = _build_paths(schema, prefix, "x")
    grad_paths = _build_paths(schema, prefix, "grad")
    metric_paths = _build_paths(schema, prefix, "metric", ic_norm_q_path)

    x_base = io.create_vec()
    io.read(x_base, x_paths)
    g = io.create_vec()
    io.read(g, grad_paths)
    M_diag = io.create_vec()
    io.read(M_diag, metric_paths)

    # ParallelIOHandler cross-check; aborts on mismatch.
    _check_gg_against_msadjoint(M_diag, g, gg_total)
    M_diag.destroy()

    # Build the FD-direction Vec. Mask g per --mode in place so the FD
    # perturbation only touches the requested subset of x.
    g_mode = g.duplicate()
    g.copy(g_mode)
    _mask_g_by_mode(g_mode, args.mode, io)
    g.destroy()

    if args.mode == "full":
        gg = gg_total
    elif args.mode == "ctrl":
        gg = ctrl_sum
    else:  # ic
        gg = ic_sum

    PETSc.Sys.Print("")
    PETSc.Sys.Print(f"Mode: {args.mode}")
    PETSc.Sys.Print(f"Baseline   J0    = {j0: .12e}")
    PETSc.Sys.Print(f"           <g,g> total (ctrl + ic) = {gg_total: .12e}")
    PETSc.Sys.Print(f"           sum control_forcing IP  = {ctrl_sum: .12e}")
    PETSc.Sys.Print(f"           sum ic IP               = {ic_sum: .12e}")
    PETSc.Sys.Print(f"           reference for this mode = {gg: .12e}")
    if gg <= 0.0:
        PETSc.Sys.Print(
            f"WARN: reference = {gg} is not positive; the FD test is meaningless."
        )
    PETSc.Sys.Print("")
    PETSc.Sys.Print(
        f"{'h':>14} {'Jh':>22} {'slope':>22} {'rel_err':>14} {'order':>8}"
    )

    try:
        results = []  # list of (h, Jh, slope, err)
        x_new = io.create_vec()
        for h in h_list:
            x_base.copy(x_new)
            x_new.axpy(h, g_mode)        # x_new = x_base + h * g_mode
            io.write(x_new, x_paths)

            _spawn_and_wait(
                "./msforward",
                ["--input", "magudi.inp", "--output", "Jh.txt"],
                N_forward,
            )
            jh = _read_scalar("Jh.txt", comm)
            slope = (jh - j0) / h
            err = abs(slope - gg) / abs(gg) if gg != 0.0 else abs(slope - gg)

            if results:
                order = (math.log10(err / results[-1][3])
                         / math.log10(h / results[-1][0])) if err > 0 else float("nan")
            else:
                order = float("nan")
            PETSc.Sys.Print(
                f"{h: 14.4e} {jh: 22.14e} {slope: 22.14e} {err: 14.4e} {order: 8.3f}"
            )
            results.append((h, jh, slope, err))
        x_new.destroy()

        # Find the prefix of rows up to the global error minimum -- the
        # decreasing regime before the roundoff floor.
        errs = [r[3] for r in results]
        if not errs:
            PETSc.Sys.Print("FAIL: no FD steps ran")
            return 1
        i_min = errs.index(min(errs))
        pre_floor = results[:i_min + 1]

        if len(pre_floor) < 2:
            PETSc.Sys.Print("FAIL: error did not decrease for any h")
            return 1

        orders = [
            math.log10(pre_floor[i + 1][3] / pre_floor[i][3])
            / math.log10(pre_floor[i + 1][0] / pre_floor[i][0])
            for i in range(len(pre_floor) - 1)
        ]
        n_bad = sum(o < order_threshold for o in orders)

        PETSc.Sys.Print("")
        PETSc.Sys.Print(
            f"Pre-floor orders ({len(orders)}): "
            f"{' '.join(f'{o:.3f}' for o in orders)}"
        )
        PETSc.Sys.Print(
            f"Bad orders (< {order_threshold}): {n_bad} / {len(orders)} "
            f"(allowed: {max_bad_orders})"
        )

        if n_bad <= max_bad_orders:
            PETSc.Sys.Print("PASS: gradient is discrete-exact within FD precision")
            return 0
        PETSc.Sys.Print("FAIL: gradient is not discrete-exact")
        return 1

    finally:
        # Restore the base point regardless of pass/fail. io.write preserves
        # each .ic.q's PLOT3D aux header (the existing file's header is
        # carried over -- see ParallelIOHandler._write_one_q), so msforward /
        # msadjoint see the original timestep and time after the test exits.
        io.write(x_base, x_paths)
        x_base.destroy()
        g_mode.destroy()


if __name__ == "__main__":
    sys.exit(main())
