"""Both magudi_utils optimization drivers.

`tao_main` runs the TAO-managed L-BFGS path (`tao.solve()` owns the outer
loop, two-loop recursion, line search, history). This is the original
magudi-optim behavior, intended for the "Case A-long" regime in
DESIGN.md (many outer iterations fit per SLURM allocation).

`manual_main` instantiates `ManualBQNLS` and runs its dispatch. The class
executes ONE step of the L-BFGS outer state machine per Python invocation,
persists everything to disk, and exits with code EXIT_RESUME (= 42) so a
surrounding bash wrapper resubmits. Intended for the "Case A" regime where
one SLURM allocation barely fits one (msforward + msadjoint). See
PLAN_case_a.md and DESIGN.md Section 6.

Both drivers share `setup_metric` and `_make_tao` so a YAML written for one
needs only `optimization.type` (and possibly `optimization.line_search.type`)
flipped to switch.
"""
import os
from collections import deque
from typing import Tuple

# mpi4py MUST import before petsc4py so PETSc binds to MPI.COMM_WORLD.
from mpi4py import MPI  # noqa: F401  (import-order side effect)
from petsc4py import PETSc

from .line_search import (HUGE_THRESHOLD, LSDecision, load_line_search,
                          make_line_search)
from .optim_state import (EXIT_CONVERGED, EXIT_RESUME, OptimState,
                          PHASE_CONVERGED, PHASE_FRESH, PHASE_LS_ACCEPTED,
                          PHASE_LS_TRIAL, read_state, write_state)
from .parallel_io import ParallelIOHandler, launch_and_wait


LMVM_DEPTH = 5  # matches BQNLS default Max. storage = 5


# ---------------------------------------------------------------------------
# Shared setup helpers (used by both tao_main and ManualBQNLS)
# ---------------------------------------------------------------------------


def setup_metric(io: ParallelIOHandler, comm) -> Tuple[PETSc.Vec, PETSc.Vec, PETSc.Vec]:
    """Load M_diag from disk and form D = sqrt(M), D_inv = 1/D.

    Same construction the current optim.py main() does at lines 117-139;
    factored out so tao_main and ManualBQNLS share it. Caller owns the
    returned Vecs and is responsible for .destroy()-ing them at shutdown.
    """
    M_diag = io.read_metric()
    local_min = (
        float(M_diag.getArray(readonly=True).min())
        if M_diag.getLocalSize() > 0 else float("inf")
    )
    global_min = comm.allreduce(local_min, op=MPI.MIN)
    if global_min <= 0.0:
        raise SystemExit(
            f"M_diag has non-positive entries (min={global_min}); cannot form sqrt(M). "
            f"Check controller mollifier coverage and state_controllability."
        )
    D = M_diag.duplicate(); M_diag.copy(D)
    D.sqrtabs()                 # D = sqrt(M_diag) = M^(1/2)
    D_inv = D.duplicate(); D.copy(D_inv)
    D_inv.reciprocal()          # D_inv = M^(-1/2)
    return M_diag, D, D_inv


def _make_tao(pcomm, fg, g_template, grtol, max_iter):
    """Construct a configured TAO instance of type bqnls.

    Used by tao_main as the actual outer-loop driver, and by ManualBQNLS
    purely as a container for the MatLMVM history (it calls tao.setUp() but
    never tao.solve()). The fg callback in ManualBQNLS is a stub that raises
    if ever invoked.
    """
    tao = PETSc.TAO().create(comm=pcomm)
    tao.setType("bqnls")
    tao.setObjectiveGradient(fg, g_template)
    tao.setMaximumIterations(max_iter)
    tao.setTolerances(grtol=grtol)
    return tao


def _configure_tao_line_search(cfg, verbose: bool, ls_log_file: str) -> None:
    """Push the optimization.line_search.* knobs into PETSc options for tao_main.

    No-op for ManualBQNLS (TAO's line search machinery is not active there).
    """
    ls_type = cfg.getInput(
        ["optimization", "line_search", "type"], fallback="more-thuente")
    ls_init_step = cfg.getInput(
        ["optimization", "line_search", "initial_step_size"], fallback=1.0)
    ls_step_max = cfg.getInput(
        ["optimization", "line_search", "max_step"], fallback=1.0e+15)
    ls_step_min = cfg.getInput(
        ["optimization", "line_search", "min_step"], fallback=1.0e-20)
    ls_stol = cfg.getInput(
        ["optimization", "line_search", "bracket_tol"], fallback=1.0e-4)
    ls_max_funcs = cfg.getInput(
        ["optimization", "line_search", "max_funcs"], fallback=30)
    PETSc.Options().setValue("tao_recycle_history", True)
    PETSc.Options().setValue("tao_ls_type", ls_type)
    PETSc.Options().setValue("tao_ls_stepinit", ls_init_step)
    PETSc.Options().setValue("tao_ls_stepmax", ls_step_max)
    PETSc.Options().setValue("tao_ls_stepmin", ls_step_min)
    PETSc.Options().setValue("tao_ls_rtol", ls_stol)
    PETSc.Options().setValue("tao_ls_max_funcs", ls_max_funcs)
    if verbose:
        PETSc.Options().setValue("tao_ls_monitor", ls_log_file)
    PETSc.Options().setValue("options_left", True)


def _read_resource_counts(cfg) -> Tuple[int, int, int]:
    N_forward = int(cfg.getInput(
        ["resource_distribution", "jobs", "forward"], datatype=list)[1])
    N_adjoint = int(cfg.getInput(
        ["resource_distribution", "jobs", "adjoint"], datatype=list)[1])
    N_petsc = int(cfg.getInput(
        ["resource_distribution", "jobs", "petsc"], datatype=list)[1])
    return N_forward, N_adjoint, N_petsc


# ---------------------------------------------------------------------------
# TAO driver
# ---------------------------------------------------------------------------


def tao_main(args, cfg, comm, pcomm) -> int:
    """The original magudi-optim flow: tao.solve() owns the outer loop.

    Lifted from optim.py main() at the optimization.type dispatch split.
    Logic unchanged.
    """
    N_forward, N_adjoint, N_petsc = _read_resource_counts(cfg)
    if comm.Get_size() != N_petsc:
        raise SystemExit(
            f"mpirun -n mismatch: launched with {comm.Get_size()} ranks but "
            f"resource_distribution.jobs.petsc requests {N_petsc}"
        )

    penalty_weight = cfg.getInput(
        ["magudi", "forced_inputs", "time_splitting", "penalty_weight"],
        datatype=float)
    checkpoint_path = "y.petsc"
    history_path = "history.petsc"

    grtol = cfg.getInput(["optimization", "tolerance"], fallback=1.0e-8)

    io = ParallelIOHandler(cfg, comm=comm, pcomm=pcomm)
    io.set_magudi_inp()
    io.report_balance()
    prefix = io.prefix
    ls_log_file = cfg.getInput(
        ["optimization", "line_search", "log_file"],
        fallback=f"{prefix}.line_search.h5")

    N_total = io.global_size

    M_diag, D, D_inv = setup_metric(io, comm)

    y = io.create_vec()
    g_y_template = y.duplicate(); g_y_template.set(0.0)

    iter_log = []
    snapshots = deque(maxlen=LMVM_DEPTH + 1)
    fg_count = [0]
    n_spawn = [0]
    pending = [None]
    ls_iter_seen = [-1]
    ls_bracket = [None]

    def fg(tao, y_vec, g_vec):
        io.cleanup_iteration_artifacts()

        x_dist = y_vec.duplicate()
        x_dist.pointwiseMult(y_vec, D_inv)
        io.write_x(x_dist)
        x_dist.destroy()
        comm.Barrier()

        launch_and_wait("./msforward", ["--input", "magudi.inp"], N_forward, args.mode)
        n_spawn[0] += 1

        J = io.read_scalar(io.j_path)
        if J > HUGE_THRESHOLD:
            PETSc.Sys.Print(
                f"fg eval #{fg_count[0]+1:4d}: msforward returned HUGE "
                f"(J = {J:.3e}); skipping msadjoint and returning inf to TAO."
            )
            g_vec.set(0.0)
            fg_count[0] += 1
            iter_log.append(float("inf"))
            return float("inf")

        launch_and_wait("./msadjoint", ["--input", "magudi.inp"], N_adjoint, args.mode)
        n_spawn[0] += 1

        gg_local = io.read_scalar(io.gg_path)

        raw_grad = io.read_grad()
        g_vec.pointwiseMult(raw_grad, D)
        raw_grad.destroy()

        if args.debug:
            gg_check = g_vec.dot(g_vec)
            rel = abs(gg_check - gg_local) / max(abs(gg_local), 1e-30)
            PETSc.Sys.Print(
                f"fwd eval #{fg_count[0]+1:4d}  J = {J: .5e}  "
                f"<g_y,g_y> = {gg_check: .5e}  vs gg_local = {gg_local: .5e}  "
                f"rel = {rel: .2e}"
            )

            it = tao.getIterationNumber()
            if it != ls_iter_seen[0]:
                ls_iter_seen[0] = it
                if ls_bracket[0] is not None:
                    ls_bracket[0][0].destroy()
                    ls_bracket[0][1].destroy()
                y_init = y_vec.duplicate(); y_vec.copy(y_init)
                g_init = g_vec.duplicate(); g_vec.copy(g_init)
                ls_bracket[0] = (y_init, g_init)
                PETSc.Sys.Print(f"    [LS#{it}] bracket start captured")
            else:
                y_init, g_init = ls_bracket[0]
                dy = y_vec.duplicate(); y_vec.copy(dy); dy.axpy(-1.0, y_init)
                ginit_dot_dy = g_init.dot(dy)
                gnow_dot_dy = g_vec.dot(dy)
                dy_norm = dy.norm()
                g_init_norm = g_init.norm()
                dy.destroy()
                denom = dy_norm * g_init_norm
                if denom > 0.0:
                    cos_dy_g0 = ginit_dot_dy / denom
                    cos_str = f"cos(dy,g0) = {cos_dy_g0: .6f}"
                else:
                    cos_str = "cos(dy,g0) = NaN"
                if abs(ginit_dot_dy) > 0.0:
                    ratio = gnow_dot_dy / ginit_dot_dy
                    ratio_str = f"|dg/dg0| = {abs(ratio):.3f}"
                else:
                    ratio_str = "|dg/dg0| = NaN"
                PETSc.Sys.Print(
                    f"    [LS#{it} trial] |alpha*d| = {dy_norm: .3e}  "
                    f"alpha*g0.d = {ginit_dot_dy: .3e}  "
                    f"alpha*g.d  = {gnow_dot_dy: .3e}  {ratio_str}  {cos_str}"
                )

        fg_count[0] += 1
        iter_log.append(J)
        return J

    def monitor(tao):
        it = tao.getIterationNumber()
        _, f_val, gnorm, _, _, _ = tao.getSolutionStatus()
        PETSc.Sys.Print(
            f"  iter {it:4d}  J = {f_val: .6e}  |g_y| = {gnorm: .6e}"
        )
        y_prev = pending[0][0] if pending[0] is not None else None
        if pending[0] is not None:
            snapshots.append(pending[0])
        x_now = tao.getSolution()
        g_now = tao.getGradient()[0]
        vx = x_now.duplicate(); x_now.copy(vx)
        vg = g_now.duplicate(); g_now.copy(vg)
        pending[0] = (vx, vg)

        if y_prev is not None:
            dy = vx.duplicate(); vx.copy(dy); dy.axpy(-1.0, y_prev)
            displacement = float(dy.norm())
            dy.destroy()
        else:
            displacement = 0.0
        io.log_iteration(penalty_weight=penalty_weight, line_step=displacement)

    tao = _make_tao(pcomm, fg, g_y_template, grtol, args.max_iter)
    _configure_tao_line_search(cfg, args.verbose, ls_log_file)
    try:
        tao.setMonitor(monitor)
    except (AttributeError, TypeError):
        PETSc.Options().setValue("tao_monitor", "")
    tao.setFromOptions()

    if args.verbose:
        PETSc.Sys.Print("--- TAO configuration ---")
        tao.view()
        PETSc.Sys.Print("-------------------------")

    PETSc.Sys.Print(
        f"magudi-optim (tao): prefix={prefix} N_total={N_total} "
        f"n_actuator={io.n_actuator} n_ic={io.n_ic} "
        f"N_petsc={N_petsc} N_forward={N_forward} N_adjoint={N_adjoint}"
    )

    tao.setSolution(y)
    tao.setUp()
    resumed = io.load_tao_state(tao, checkpoint_path, history_path)
    if not resumed:
        io.seed_zero_actuator_files()
        baseline_x = io.read_x()
        y.pointwiseMult(baseline_x, D)
        baseline_x.destroy()
        comm.Barrier()

    tao.solve(y)

    io.save_tao_state(tao, checkpoint_path, history_path, snapshots)

    reason = tao.getConvergedReason()
    final_J = iter_log[-1] if iter_log else float("nan")
    initial_J = iter_log[0] if iter_log else float("nan")
    PETSc.Sys.Print("")
    PETSc.Sys.Print(f"Converged reason: {reason}")
    PETSc.Sys.Print(f"Initial J = {initial_J: .6e}")
    PETSc.Sys.Print(f"Final   J = {final_J: .6e}")
    PETSc.Sys.Print(
        f"fg calls = {fg_count[0]};  spawn calls = {n_spawn[0]}"
    )
    return 0


# ---------------------------------------------------------------------------
# Manual (Case-A) driver
# ---------------------------------------------------------------------------
#
# The manual driver executes ONE step of the L-BFGS state machine per Python
# invocation. State persists via:
#
#   y.petsc           current accepted iterate
#   g.petsc           accepted y-space gradient
#   d.petsc           current search direction (-H g)
#   history.petsc     L-BFGS (x_k, g_k) snapshots + J0 diagonal
#   lbfgs_state.json  scalar state machine
#
# Each invocation runs at most one msforward (ls_trial) or one msforward + one
# msadjoint (fresh / ls_accepted). The bash wrapper resubmits while exit
# code == EXIT_RESUME (42); stops on EXIT_CONVERGED (0).

# Path constants for the manual driver. Module-level so a bash wrapper can
# grep them out of the source if needed.
_Y_PATH = "y.petsc"
_G_PATH = "g.petsc"
_D_PATH = "d.petsc"
_HIST_PATH = "history.petsc"
_STATE_PATH = "lbfgs_state.json"


def _manual_dummy_fg(tao, y, g):
    """TAO requires setObjectiveGradient before setUp(); ManualBQNLS never
    invokes tao.solve(), so this should be unreachable. Raising catches
    accidental tao.solve() calls during refactoring."""
    raise RuntimeError(
        "manual driver invoked TAO's fg callback; tao.solve() must not be "
        "called in the manual code path"
    )


class ManualBQNLS:
    """L-BFGS outer loop driven from Python; one state-machine step per
    Python invocation.

    Per-allocation budget: at most one msforward (during ls_trial) or one
    msforward + one msadjoint (during fresh / ls_accepted). The dispatch in
    `run()` reads lbfgs_state.json to pick which step to execute next.
    """

    def __init__(self, cfg, mode: str, max_iter: int, comm, pcomm,
                 verbose: bool = False):
        self.cfg = cfg
        self.mode = mode
        self.max_iter = max_iter
        self.comm = comm
        self.pcomm = pcomm
        # When True, ls_accepted loads g.petsc to compute the diagnostic
        # <s, y> curvature for the run log. When False (default), that read
        # is skipped -- the value is informational only; PETSc's updateLMVM
        # checks <s, y> > 0 internally and skips the pair if it fails.
        self.verbose = verbose

        # Resource layout (must match the launched mpirun width).
        self.N_forward, self.N_adjoint, self.N_petsc = _read_resource_counts(cfg)
        if comm.Get_size() != self.N_petsc:
            raise SystemExit(
                f"mpirun -n mismatch: launched with {comm.Get_size()} ranks "
                f"but resource_distribution.jobs.petsc requests {self.N_petsc}"
            )

        # Pull the handful of YAML knobs the dispatcher and per-phase methods
        # need, so methods don't have to re-read the cfg.
        self.penalty_weight = cfg.getInput(
            ["magudi", "forced_inputs", "time_splitting", "penalty_weight"],
            datatype=float)
        self.grtol = cfg.getInput(["optimization", "tolerance"], fallback=1.0e-8)
        self.ls_type = cfg.getInput(
            ["optimization", "line_search", "type"], fallback="golden_section")

        # I/O handler. set_magudi_inp() rewrites magudi.inp from
        # forced_inputs; report_balance() prints the rank-vs-slot mapping.
        self.io = ParallelIOHandler(cfg, comm=comm, pcomm=pcomm)
        self.io.set_magudi_inp()
        self.io.report_balance()

        # Metric: y = D * x with D = M^(1/2). D is the only piece we keep
        # around for write/read transformations; M_diag itself can be
        # released once D and D_inv have been derived.
        M_diag, self.D, self.D_inv = setup_metric(self.io, comm)
        M_diag.destroy()

        # TAO instance: container for MatLMVM (history + J0 diagonal) and
        # checkpoint-format compatibility with tao_main's save/load. We
        # never invoke tao.solve() -- the outer loop is owned by this class.
        self.y = self.io.create_vec()
        g_template = self.y.duplicate(); g_template.set(0.0)
        self.tao = _make_tao(
            pcomm, _manual_dummy_fg, g_template, self.grtol, self.max_iter)
        self.tao.setFromOptions()
        self.tao.setSolution(self.y)
        self.tao.setUp()

    # -- low-level executors ------------------------------------------------

    def msforward(self, y_new) -> float:
        """Cleanup iteration artifacts, convert y_new to x, run msforward,
        read and return J. Mutates `y_new` only via reads.
        """
        self.io.cleanup_iteration_artifacts()
        x_dist = y_new.duplicate()
        x_dist.pointwiseMult(y_new, self.D_inv)
        self.io.write_x(x_dist)
        x_dist.destroy()
        self.comm.Barrier()
        launch_and_wait(
            "./msforward", ["--input", "magudi.inp"], self.N_forward, self.mode)
        return self.io.read_scalar(self.io.j_path)

    def msadjoint(self) -> PETSc.Vec:
        """Run msadjoint at the x written by the previous msforward; read
        raw_grad; scale by D; return g_y = D * raw_grad as a fresh Vec.

        Caller owns the returned Vec and must .destroy() it.
        """
        launch_and_wait(
            "./msadjoint", ["--input", "magudi.inp"], self.N_adjoint, self.mode)
        raw_grad = self.io.read_grad()
        g_y = self.io.create_vec()
        g_y.pointwiseMult(raw_grad, self.D)
        raw_grad.destroy()
        return g_y

    def persist_petsc(self, name_to_vec) -> None:
        """Write a {path: Vec} mapping as PETSc binary files."""
        for path, vec in name_to_vec.items():
            self.io.write_petsc(vec, path)

    def direction_from_lbfgs(self, g_accepted: PETSc.Vec) -> PETSc.Vec:
        """Compute d = -H * g_accepted using the TAO LMVM matrix.

        For an empty history MatLMVM applies its J0 diagonal (default scaled
        identity), so d_0 = -g (steepest descent). After updateLMVM calls
        this becomes the proper L-BFGS two-loop recursion direction.
        """
        M = self.tao.getLMVMMat()
        d = g_accepted.duplicate()
        M.mult(g_accepted, d)
        d.scale(-1.0)
        return d

    # -- top-level dispatch -------------------------------------------------

    def run(self) -> int:
        """Read state, perform phase-specific resume, dispatch, return exit
        code.

        Resume policy:
          ls_trial   -> only y.petsc is loaded into self.y (no MatLMVM).
                        Direction is read from d.petsc inside the handler.
          ls_accepted-> full load_tao_state: y.petsc into self.y AND replay
                        history.petsc through MatLMVM (needed for the next
                        direction computation).
          fresh      -> no load; initial_step cold-starts from baseline.
          converged  -> no load; idempotent exit.

        `resumed` is determined by y.petsc existence rather than as a
        side effect of load_tao_state, so a phase that doesn't need
        MatLMVM doesn't pay for the history.petsc read.
        """
        state = read_state(_STATE_PATH)
        if state is None:
            state = OptimState(iter=0, phase=PHASE_FRESH)

        resumed = os.path.exists(_Y_PATH)
        if state.phase != PHASE_FRESH and state.phase != PHASE_CONVERGED \
                and not resumed:
            raise SystemExit(
                f"lbfgs_state.json claims phase={state.phase!r} but "
                f"y.petsc is missing. State inconsistency; "
                f"rm lbfgs_state.json to restart from scratch."
            )
        if resumed and state.phase == PHASE_FRESH:
            PETSc.Sys.Print(
                "warning: phase=fresh but y.petsc present; ignoring stale "
                "Vecs and re-bootstrapping. Delete *.petsc lbfgs_state.json "
                "to silence."
            )

        PETSc.Sys.Print(
            f"magudi-optim (manual): prefix={self.io.prefix} "
            f"N_total={self.io.global_size} "
            f"n_actuator={self.io.n_actuator} n_ic={self.io.n_ic} "
            f"N_petsc={self.N_petsc} N_forward={self.N_forward} "
            f"N_adjoint={self.N_adjoint}"
        )
        PETSc.Sys.Print(
            f"  state: iter={state.iter} phase={state.phase} "
            f"f_baseline={state.f_baseline:.6e}"
        )

        # Hard max-iter guardrail (covers stale state with iter beyond cap).
        if state.iter >= self.max_iter and state.phase != PHASE_CONVERGED:
            state.phase = PHASE_CONVERGED
            state.converged_reason = "max_iter"
            if self.comm.Get_rank() == 0:
                write_state(state, _STATE_PATH)
            self.comm.Barrier()
            PETSc.Sys.Print(
                f"iter {state.iter} >= max_iter {self.max_iter}; converged."
            )
            return EXIT_CONVERGED

        if state.phase == PHASE_FRESH:
            return self.initial_step(state)

        if state.phase == PHASE_LS_TRIAL:
            # Lightweight resume: only y.petsc into self.y. MatLMVM is not
            # needed -- the trial step only needs self.y and d.petsc.
            viewer = PETSc.Viewer().createBinary(
                _Y_PATH, mode="r", comm=self.pcomm)
            self.y.load(viewer)
            viewer.destroy()
            return self.ls_trial(state)

        if state.phase == PHASE_LS_ACCEPTED:
            # Full resume: y.petsc + history.petsc through MatLMVM. The
            # direction computation at the end of ls_accepted needs the
            # replayed history.
            self.io.load_tao_state(self.tao, _Y_PATH, _HIST_PATH)
            return self.ls_accepted(state)

        if state.phase == PHASE_CONVERGED:
            PETSc.Sys.Print(
                f"Already converged at iter {state.iter} "
                f"({state.converged_reason!r}); exiting clean."
            )
            return EXIT_CONVERGED

        raise SystemExit(f"unknown phase {state.phase!r} in {_STATE_PATH}")

    # -- phase handlers -----------------------------------------------------

    def initial_step(self, state: OptimState) -> int:
        """Iter-0 bootstrap. Pays one msforward + one msadjoint to seed
        J(y_0) and g(y_0), then sets up the first LS trial and saves state.
        """
        # Initialize accepted y from the baseline IC + zero control forcing.
        self.io.seed_zero_actuator_files()
        baseline_x = self.io.read_x()
        self.y.pointwiseMult(baseline_x, self.D)
        baseline_x.destroy()
        self.comm.Barrier()

        J = self.msforward(self.y)
        if J > HUGE_THRESHOLD:
            raise SystemExit(
                f"baseline forward returned HUGE (J={J:.3e}); cannot "
                f"bootstrap. Check IC files and resource layout."
            )
        g_accepted = self.msadjoint()

        # Prime MatLMVM with (y_0, g_0). PETSc's updateLMVM expects the
        # (X, G) pair and diffs internally against the previous call to
        # form the next (s, y) secant pair. The first call adds no pair;
        # subsequent calls in ls_accepted build the L-BFGS history.
        self.tao.getLMVMMat().updateLMVM(self.y, g_accepted)

        # Steepest-descent direction for iter 0 (no history yet).
        d = g_accepted.duplicate(); g_accepted.copy(d); d.scale(-1.0)
        gT_d = g_accepted.dot(d)

        ls = make_line_search(self.cfg, self.ls_type)
        pending_alpha = ls.first_trial(J, float(gT_d))

        # Persist y.petsc + history.petsc via save_tao_state. The snapshots
        # deque carries (y_0, g_0) so the next allocation's load_tao_state
        # re-primes MatLMVM with the same (X, G); downstream updateLMVM
        # calls in ls_accepted then form (s, y) secant pairs correctly.
        snapshots = deque(maxlen=LMVM_DEPTH + 1)
        vy = self.y.duplicate(); self.y.copy(vy)
        vg = g_accepted.duplicate(); g_accepted.copy(vg)
        snapshots.append((vy, vg))
        self.tao.setIterationNumber(0)
        self.io.save_tao_state(self.tao, _Y_PATH, _HIST_PATH, snapshots)
        vy.destroy(); vg.destroy()

        # d.petsc: required by every ls_trial. g.petsc: only if verbose
        # (used solely for the <s, y> diagnostic in subsequent ls_accepted).
        self.persist_petsc({_D_PATH: d})
        if self.verbose:
            self.persist_petsc({_G_PATH: g_accepted})

        state.iter = 0
        state.phase = PHASE_LS_TRIAL
        state.f_baseline = float(J)
        state.gT_d_baseline = float(gT_d)
        state.line_search = ls.to_dict()
        if self.comm.Get_rank() == 0:
            write_state(state, _STATE_PATH)
        self.comm.Barrier()

        # Iter-0 baseline forward+adjoint is logged with zero line step,
        # matching tao_main's first monitor row.
        self.io.log_iteration(
            penalty_weight=self.penalty_weight, line_step=0.0)

        g0_norm = float(g_accepted.norm())
        PETSc.Sys.Print(
            f"fresh: iter 0 baseline J = {J: .6e}  "
            f"|g_y^0| = {g0_norm: .6e}"
        )
        PETSc.Sys.Print(
            f"       first LS trial at alpha = {pending_alpha: .6e}; "
            f"exit RESUME"
        )

        g_accepted.destroy(); d.destroy()
        return EXIT_RESUME

    def ls_trial(self, state: OptimState) -> int:
        """One msforward at LS.pending_alpha; update LS; persist; exit.

        Decision branches:
          CONTINUE / CRASHED -> store the LS's new pending_alpha; phase
                                stays ls_trial; exit RESUME.
          ACCEPT             -> store alpha_star; phase transitions to
                                ls_accepted (next allocation runs fwd+adj
                                at alpha_star); exit RESUME.
        """
        ls = load_line_search(state.line_search, self.cfg)
        pending_alpha = ls.pending_alpha
        if pending_alpha is None:
            raise SystemExit(
                "phase=ls_trial but line_search has no pending_alpha. "
                "Stale state?"
            )

        d = self.io.read_petsc(_D_PATH)

        y_trial = self.y.duplicate(); self.y.copy(y_trial)
        y_trial.axpy(pending_alpha, d)
        J_trial = self.msforward(y_trial)
        y_trial.destroy()

        decision, next_alpha, accepted_J = ls.record_and_next(
            pending_alpha, J_trial)

        PETSc.Sys.Print(
            f"iter {state.iter} ls_trial  alpha = {pending_alpha: .6e}  "
            f"J = {J_trial: .6e}  -> {decision}  "
            f"next_alpha = {next_alpha: .6e}"
        )

        if decision == LSDecision.ACCEPT:
            # Hand off to ls_accepted. The just-finished forward's snapshots
            # are at pending_alpha; the accepted alpha_star (= LS's middle b)
            # may differ, so the next allocation re-runs msforward+msadjoint.
            state.phase = PHASE_LS_ACCEPTED
            state.accepted_alpha = float(next_alpha)
            if accepted_J is not None:
                state.f_baseline = float(accepted_J)
            state.line_search = ls.to_dict()
            if self.comm.Get_rank() == 0:
                write_state(state, _STATE_PATH)
            self.comm.Barrier()
            PETSc.Sys.Print(
                f"  LS accepted alpha* = {next_alpha: .6e}  "
                f"J* = {accepted_J: .6e}; exit RESUME "
                f"(next allocation runs fwd+adj at alpha*)"
            )
            d.destroy()
            return EXIT_RESUME

        # CONTINUE or CRASHED: stay in ls_trial with the new pending_alpha.
        state.line_search = ls.to_dict()
        if self.comm.Get_rank() == 0:
            write_state(state, _STATE_PATH)
        self.comm.Barrier()
        d.destroy()
        return EXIT_RESUME

    def ls_accepted(self, state: OptimState) -> int:
        """One msforward + one msadjoint at state.accepted_alpha.

        On success: update L-BFGS history, increment iter, check
        convergence, compute new direction, transition to ls_trial for the
        next iter (or to converged if grtol / max_iter is met).
        """
        alpha_star = state.accepted_alpha
        if alpha_star <= 0.0:
            raise SystemExit(
                f"phase=ls_accepted but accepted_alpha = {alpha_star}; "
                f"stale state?"
            )

        # Search direction at the accepted point. g_accepted (the PREVIOUS
        # iterate's gradient) is loaded only when verbose: it is used solely
        # to log the diagnostic <s, y> curvature value. PETSc's updateLMVM
        # checks <s, y> > 0 internally and skips the pair if not, so the
        # logged value is informational.
        d = self.io.read_petsc(_D_PATH)

        # Compose y_new and run fwd+adj.
        y_new = self.y.duplicate(); self.y.copy(y_new)
        y_new.axpy(alpha_star, d)

        J_new = self.msforward(y_new)
        if J_new > HUGE_THRESHOLD:
            # The LS swore J(alpha_star) was finite from its bracket cache.
            # Seeing HUGE here means state corruption between allocations;
            # don't try to recover silently.
            raise SystemExit(
                f"ls_accepted forward at alpha* = {alpha_star} returned "
                f"HUGE (J = {J_new:.3e}); state inconsistency. "
                f"Re-bootstrap (rm *.petsc lbfgs_state.json) and retry."
            )
        g_new = self.msadjoint()

        update_msg = ""
        if self.verbose:
            g_accepted = self.io.read_petsc(_G_PATH)
            s = y_new.duplicate(); y_new.copy(s); s.axpy(-1.0, self.y)
            yk = g_new.duplicate(); g_new.copy(yk); yk.axpy(-1.0, g_accepted)
            sy = float(s.dot(yk))
            s.destroy(); yk.destroy(); g_accepted.destroy()
            update_msg = (
                f"<s,y> = {sy: .3e} > 0 (updateLMVM accepts)"
                if sy > 0.0 else
                f"<s,y> = {sy: .3e} <= 0 (updateLMVM likely skips)"
            )

        # Hand the new iterate+gradient to MatLMVM; it forms the secant
        # internally against the previous (X, G) primed in initial_step.
        self.tao.getLMVMMat().updateLMVM(y_new, g_new)

        # Commit the new accepted iterate.
        y_new.copy(self.y)
        y_new.destroy()

        state.iter += 1
        state.f_baseline = float(J_new)
        gnorm = float(g_new.norm())

        converged, reason = False, None
        if state.iter >= self.max_iter:
            converged, reason = True, "max_iter"
        elif gnorm < self.grtol:
            converged, reason = True, "gatol"

        diag = f"  ({update_msg})" if update_msg else ""
        PETSc.Sys.Print(
            f"iter {state.iter} ls_accepted  alpha* = {alpha_star: .6e}  "
            f"J_new = {J_new: .6e}  |g_y| = {gnorm: .6e}{diag}"
        )

        # Re-read existing snapshots so we can append the new (y_new, g_new)
        # pair and rewrite the full history. load_snapshots returns a
        # bounded deque; if it is already at capacity we destroy the
        # popped-front pair explicitly to avoid a Vec leak.
        snapshots = self.io.load_snapshots(
            _HIST_PATH, maxlen=LMVM_DEPTH + 1)
        new_vy = self.y.duplicate(); self.y.copy(new_vy)
        new_vg = g_new.duplicate(); g_new.copy(new_vg)
        if snapshots.maxlen is not None and len(snapshots) == snapshots.maxlen:
            old_vx, old_vg = snapshots.popleft()
            old_vx.destroy(); old_vg.destroy()
        snapshots.append((new_vy, new_vg))

        # save_tao_state writes y.petsc (from tao.getSolution() = self.y)
        # and history.petsc (snapshots + J0 diag from MatLMVM). The
        # iteration counter goes into history.petsc's header, so update it
        # to match state.iter before saving.
        self.tao.setIterationNumber(state.iter)
        self.io.save_tao_state(self.tao, _Y_PATH, _HIST_PATH, snapshots)
        for vx, vg in snapshots:
            vx.destroy(); vg.destroy()

        if converged:
            state.phase = PHASE_CONVERGED
            state.converged_reason = reason
            if self.comm.Get_rank() == 0:
                write_state(state, _STATE_PATH)
            self.comm.Barrier()
            self.io.log_iteration(
                penalty_weight=self.penalty_weight, line_step=alpha_star)
            PETSc.Sys.Print(f"converged ({reason}); exit CONVERGED")
            g_new.destroy(); d.destroy()
            return EXIT_CONVERGED

        # Not converged: compute the new direction and set up the next iter.
        d_new = self.direction_from_lbfgs(g_new)
        gT_d_new = float(g_new.dot(d_new))

        # Numerical safety: if d_new is not a descent direction, the L-BFGS
        # history may be corrupted. Fall back to steepest descent.
        if gT_d_new > 0.0:
            PETSc.Sys.Print(
                f"  warning: gT_d_new = {gT_d_new: .3e} > 0 (not a descent "
                f"direction); falling back to -g"
            )
            d_new.destroy()
            d_new = g_new.duplicate(); g_new.copy(d_new); d_new.scale(-1.0)
            gT_d_new = float(g_new.dot(d_new))

        # d.petsc: needed by every ls_trial. g.petsc: only if verbose
        # (consumed only by the <s, y> diagnostic in the next ls_accepted).
        self.persist_petsc({_D_PATH: d_new})
        if self.verbose:
            self.persist_petsc({_G_PATH: g_new})

        ls = make_line_search(self.cfg, self.ls_type)
        pending_alpha = ls.first_trial(J_new, gT_d_new)

        state.phase = PHASE_LS_TRIAL
        state.gT_d_baseline = gT_d_new
        state.line_search = ls.to_dict()
        if self.comm.Get_rank() == 0:
            write_state(state, _STATE_PATH)
        self.comm.Barrier()

        # Log the iter we just completed; alpha_star is the y-space step.
        self.io.log_iteration(
            penalty_weight=self.penalty_weight, line_step=alpha_star)

        PETSc.Sys.Print(
            f"  next iter: |g_y| = {gnorm: .6e}  "
            f"pending_alpha = {pending_alpha: .6e}; exit RESUME"
        )

        g_new.destroy(); d.destroy(); d_new.destroy()
        return EXIT_RESUME


def manual_main(args, cfg, comm, pcomm) -> int:
    """Entry-point wrapper around ManualBQNLS for optim.py's dispatcher."""
    driver = ManualBQNLS(
        cfg=cfg, mode=args.mode, max_iter=args.max_iter,
        comm=comm, pcomm=pcomm, verbose=args.verbose,
    )
    return driver.run()
