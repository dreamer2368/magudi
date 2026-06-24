"""Line search algorithms for the manual L-BFGS driver.

Two ingredients:

  1. A `LineSearch` abstract base that the manual driver calls without caring
     which algorithm is wired in. Each subclass owns a small piece of state
     that round-trips through `to_dict` / `from_dict` so the LS picks up
     across allocation boundaries (one Python invocation per LS trial).

  2. `GoldenSectionLineSearch`, a port of the MNBRAK + LINMIN routine from
     legacy/optimization_ver3/optimizer.py:460-658. Forward-only at every
     trial: only J at (y_accepted + alpha * d) is consulted; gradient is
     never evaluated mid-search. Adjoint runs once per accepted step, after
     the bracket has tightened.

The manual driver invokes LS like this (illustrative; see optimizers.py):

    ls = make_line_search(cfg, ls_type)             # fresh start of an iter
    alpha = ls.first_trial(f_baseline, gT_d)        # returns the first trial
    # ... run msforward at alpha, read J_alpha ...
    decision, next_alpha, accepted_J = ls.record_and_next(alpha, J_alpha)

`record_and_next` returns one of three decisions:

  CONTINUE  : `next_alpha` is the next trial; recurse.
  ACCEPT    : LS converged. `next_alpha` is the accepted alpha*; `accepted_J`
              is J at alpha* (cached from the bracket). Driver must re-run
              msforward at alpha* to seed snapshots for msadjoint (the just-
              finished trial may have been at a different alpha).
  CRASHED   : last trial returned J >= HUGE (msforward crash sentinel from
              optim.py's HUGE_THRESHOLD). LS shrinks the step toward a known
              good interior point and returns the new trial alpha.

ArmijoLineSearch is deferred but registered in the dispatch map so future
work needs no plumbing changes -- only fill in the class.
"""
import abc
import math
from typing import Any, Dict, Optional, Tuple


# Mirrors the msforward HUGE-J sentinel in optim.py:43. Keep these two in
# sync; when msforward returns J > HUGE_THRESHOLD, both the fg callback and
# the LS treat it as a crash signal.
HUGE_THRESHOLD = 1.0e300


class LSDecision:
    """Three outcomes a line search can return after one trial.

    Plain string constants rather than enum.Enum so they JSON-serialize
    transparently and bash wrappers can grep them out of state files.
    """
    CONTINUE = "continue"
    ACCEPT = "accept"
    CRASHED = "crashed"


class LineSearch(abc.ABC):
    """Abstract one-trial-at-a-time line search.

    Subclasses are responsible for: picking the initial trial alpha, deciding
    after each trial whether to accept, continue, or shrink, and round-
    tripping their internal bracket state through `to_dict` / `from_dict`.
    """

    @abc.abstractmethod
    def first_trial(self, f_baseline: float, gT_d: float) -> float:
        """Initialize a fresh LS at y_accepted with J = f_baseline.

        gT_d is the directional derivative <g_accepted, d> at the bracket
        start. Returns the alpha to evaluate at the first trial. Drives a
        single allocation's msforward.
        """

    @abc.abstractmethod
    def record_and_next(
        self, last_alpha: float, last_J: float
    ) -> Tuple[str, float, Optional[float]]:
        """Ingest one trial result; return (decision, next_alpha, accepted_J).

        decision is one of LSDecision.{CONTINUE, ACCEPT, CRASHED}.
        For CONTINUE / CRASHED, accepted_J is None.
        For ACCEPT, accepted_J is J at the accepted alpha (cached from the
        bracket; not necessarily the just-evaluated J).
        """

    @abc.abstractmethod
    def to_dict(self) -> Dict[str, Any]:
        """Serialize internal state. Must include a "type" key so dispatch
        can pick the right subclass on resume."""

    @classmethod
    @abc.abstractmethod
    def from_dict(cls, d: Dict[str, Any], cfg) -> "LineSearch":
        """Reconstruct an LS instance from its serialized state."""


# ---------------------------------------------------------------------------
# Golden-section LS: port of ver3 MNBRAK + LINMIN
# ---------------------------------------------------------------------------


def _parabolic_interp(steps, Js):
    """Three-point parabolic minimum on (alpha, J) -- same expression as
    ver3 optimizer.py:295-301."""
    a, b, c = steps
    fa, fb, fc = Js
    num = (b - a) ** 2 * (fb - fc) - (b - c) ** 2 * (fb - fa)
    den = (b - a) * (fb - fc) - (b - c) * (fb - fa)
    if abs(den) < 1.0e-30:
        # Degenerate three points (colinear J): caller's safeguard kicks in.
        return float("nan")
    return b - 0.5 * num / den


class GoldenSectionLineSearch(LineSearch):
    """ver3 MNBRAK + LINMIN, forward-only.

    State (round-tripped through to_dict / from_dict):

      a, Ja      : lower bracket point, J at it.  Always set after first_trial.
      b, Jb      : middle bracket point, J at it. Set during MNBRAK; remains
                   the LS's "current best" through LINMIN.
      c, Jc      : upper bracket point, J at it.  Set when the bracket closes
                   (MNBRAK -> LINMIN transition) and refined during LINMIN.
      bracketed  : True after MNBRAK finds (Ja > Jb < Jc).
      pending_alpha : alpha at which the next msforward will run. Persisted
                      because the driver reads it before the trial without
                      calling first_trial / record_and_next again.
      trial_count   : guardrail against runaway expansion; aborts CRASHED
                      with min_step if trial_count > max_funcs.
      crash_alphas  : alphas where the previous trial returned J > HUGE; used
                      to clamp expansion when an old crash is the only upper
                      bound we know about.
    """

    # ver3 constants (optimizer.py:295-301, 460-560).
    GOLDEN_RATIO = (1.0 + math.sqrt(5.0)) / 2.0    # ~ 1.618
    CR = 1.0 - 1.0 / GOLDEN_RATIO                  # ~ 0.382, golden-section step
    EPS = 1.0e-15

    TYPE_TAG = "golden_section"

    def __init__(self, initial_step: float, max_step: float, min_step: float,
                 bracket_tol: float, max_funcs: int):
        self.initial_step = float(initial_step)
        self.max_step = float(max_step)
        self.min_step = float(min_step)
        self.bracket_tol = float(bracket_tol)
        self.max_funcs = int(max_funcs)

        self.a: Optional[float] = None
        self.Ja: Optional[float] = None
        self.b: Optional[float] = None
        self.Jb: Optional[float] = None
        self.c: Optional[float] = None
        self.Jc: Optional[float] = None
        self.bracketed: bool = False
        self.pending_alpha: Optional[float] = None
        self.trial_count: int = 0
        self.crash_alphas: list = []

    # -- public API ----------------------------------------------------------

    def first_trial(self, f_baseline: float, gT_d: float) -> float:
        # ver3 sets a = 0 with Ja = J(baseline); first trial alpha = initial_step.
        self.a = 0.0
        self.Ja = float(f_baseline)
        self.b = None
        self.Jb = None
        self.c = None
        self.Jc = None
        self.bracketed = False
        self.trial_count = 0
        self.crash_alphas = []
        self.pending_alpha = max(self.initial_step, self.min_step)
        return self.pending_alpha

    def record_and_next(
        self, last_alpha: float, last_J: float
    ) -> Tuple[str, float, Optional[float]]:
        self.trial_count += 1

        # Guardrail. We do not throw because in the manual driver an exception
        # here would corrupt the JSON state and force operator recovery.
        # Surfacing CRASHED with min_step lets the LS exit gracefully on the
        # next iteration's first-trial check.
        if self.trial_count > self.max_funcs:
            self.pending_alpha = self.min_step
            return (LSDecision.CRASHED, self.min_step, None)

        if last_J > HUGE_THRESHOLD:
            return self._on_crash(last_alpha)

        if not self.bracketed:
            return self._mnbrak(last_alpha, last_J)
        return self._linmin(last_alpha, last_J)

    # -- serialization -------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "type": self.TYPE_TAG,
            "initial_step": self.initial_step,
            "max_step": self.max_step,
            "min_step": self.min_step,
            "bracket_tol": self.bracket_tol,
            "max_funcs": self.max_funcs,
            "a": self.a, "Ja": self.Ja,
            "b": self.b, "Jb": self.Jb,
            "c": self.c, "Jc": self.Jc,
            "bracketed": self.bracketed,
            "pending_alpha": self.pending_alpha,
            "trial_count": self.trial_count,
            "crash_alphas": list(self.crash_alphas),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any], cfg) -> "GoldenSectionLineSearch":
        # cfg is unused on the resume path because all knobs round-tripped
        # in the state; kept in the signature to match the LineSearch ABC
        # and to allow future LS types that need cfg lookups at reload.
        del cfg
        ls = cls(
            initial_step=d["initial_step"],
            max_step=d["max_step"],
            min_step=d["min_step"],
            bracket_tol=d["bracket_tol"],
            max_funcs=d["max_funcs"],
        )
        ls.a, ls.Ja = d.get("a"), d.get("Ja")
        ls.b, ls.Jb = d.get("b"), d.get("Jb")
        ls.c, ls.Jc = d.get("c"), d.get("Jc")
        ls.bracketed = bool(d.get("bracketed", False))
        ls.pending_alpha = d.get("pending_alpha")
        ls.trial_count = int(d.get("trial_count", 0))
        ls.crash_alphas = list(d.get("crash_alphas", []))
        return ls

    # -- MNBRAK --------------------------------------------------------------

    def _on_crash(self, alpha_crashed: float) -> Tuple[str, float, Optional[float]]:
        """Shrink toward a known good interior point. Mirrors ver3
        optimizer.py:474-482."""
        self.crash_alphas.append(alpha_crashed)
        anchor = self.b if self.b is not None else self.a
        next_alpha = anchor + self.CR * (alpha_crashed - anchor)
        if next_alpha < self.min_step:
            self.pending_alpha = self.min_step
            return (LSDecision.CRASHED, self.min_step, None)
        self.pending_alpha = next_alpha
        return (LSDecision.CRASHED, next_alpha, None)

    def _mnbrak(self, alpha: float, J: float) -> Tuple[str, float, Optional[float]]:
        """Bracket-finding pass. Maintains (a, b, c) as labels for
        (lower, interior, upper) bracket points; on entry the most recent
        trial (alpha, J) needs to be absorbed.

        Transitions in plain English (ver3 optimizer.py:483-538):

          - only `a` known so far: if alpha is even better than a -> alpha
            becomes interior b, expand outward. Else alpha was too far -> alpha
            becomes upper c, contract.
          - `b` known, no `c` yet: if J > Jb the bracket is closed (J at the
            far point exceeds the interior J). Else alpha is better than b;
            slide the bracket outward (a := b, b := alpha) and try farther.
          - `c` known, no `b` yet: symmetric. If J < Ja the bracket is
            closed; else contract.
        """
        # Case 1: have only a, plus current trial.
        if self.b is None and self.c is None:
            if J > self.Ja:
                # alpha is too far out; demote alpha to c and contract.
                self.c, self.Jc = alpha, J
                next_alpha = alpha / self.GOLDEN_RATIO
            else:
                # alpha beats a; promote to interior b and expand farther.
                self.b, self.Jb = alpha, J
                next_alpha = self._mnbrak_expand_from(alpha)
            return self._maybe_continue(next_alpha)

        # Case 2: have a and b, no c yet.
        if self.b is not None and self.c is None:
            if J > self.Jb:
                # Bracket closed: a (lo) - b (interior) - alpha (hi).
                self.c, self.Jc = alpha, J
                return self._enter_linmin()
            # alpha is even better than b; slide bracket outward.
            self.a, self.Ja = self.b, self.Jb
            self.b, self.Jb = alpha, J
            next_alpha = self._mnbrak_expand_from(alpha)
            return self._maybe_continue(next_alpha)

        # Case 3: have a and c, no b yet (alpha came from a contraction step).
        if self.b is None and self.c is not None:
            if J < self.Ja:
                # Bracket closed: a (lo) - alpha (interior) - c (hi).
                self.b, self.Jb = alpha, J
                return self._enter_linmin()
            # alpha was not better than a; contract c inward.
            self.c, self.Jc = alpha, J
            next_alpha = alpha / self.GOLDEN_RATIO
            return self._maybe_continue(next_alpha)

        # Defensive: a, b, c all set but bracketed is False. Shouldn't happen
        # if state transitions follow the rules above, but treat as a re-entry
        # into LINMIN rather than crashing.
        self.bracketed = True
        return self._linmin(alpha, J)

    def _mnbrak_expand_from(self, alpha: float) -> float:
        """Pick the next outward trial. If past crashes are known and they
        bracket alpha from above, use them to cap expansion (CR of the gap);
        otherwise expand by golden ratio up to max_step."""
        upper_crashes = [c for c in self.crash_alphas if c > alpha]
        if upper_crashes:
            cap = min(upper_crashes)
            return alpha + self.CR * (cap - alpha)
        # No known upper bound; standard golden expansion clipped to max_step.
        next_alpha = alpha * self.GOLDEN_RATIO
        if next_alpha > self.max_step:
            next_alpha = self.max_step
        return next_alpha

    def _enter_linmin(self) -> Tuple[str, float, Optional[float]]:
        """Bracket is closed. Generate the first LINMIN trial and transition."""
        self.bracketed = True
        next_alpha = self._linmin_next_trial()
        # Bracket might already be tight enough at the boundary instant.
        if (self.c - self.a) < self.b * self.bracket_tol + self.EPS:
            return self._accept_b()
        return self._maybe_continue(next_alpha)

    # -- LINMIN --------------------------------------------------------------

    def _linmin(self, alpha: float, J: float) -> Tuple[str, float, Optional[float]]:
        """Refine the bracket with the new trial point, then either accept
        the new b or generate the next trial via parabolic interpolation
        (safeguarded by golden section). Mirrors ver3 optimizer.py:562-658."""
        # Merge the new trial into (a, b, c, x) and pick the new bracket as the
        # three adjacent to the J-min.
        pts = sorted(
            [(self.a, self.Ja), (self.b, self.Jb),
             (self.c, self.Jc), (alpha, J)],
            key=lambda p: p[0],
        )
        # Clamp the min index away from the boundaries so (idx-1, idx, idx+1)
        # stays in range -- matches ver3's `max(np.argmin(Js), 1)` (lower)
        # and extends to upper.
        Js = [p[1] for p in pts]
        min_idx = min(max(Js.index(min(Js)), 1), len(pts) - 2)
        (self.a, self.Ja) = pts[min_idx - 1]
        (self.b, self.Jb) = pts[min_idx]
        (self.c, self.Jc) = pts[min_idx + 1]

        if (self.c - self.a) < self.b * self.bracket_tol + self.EPS:
            return self._accept_b()

        next_alpha = self._linmin_next_trial()
        return self._maybe_continue(next_alpha)

    def _linmin_next_trial(self) -> float:
        """Parabolic interpolation through (a, b, c) with golden-section
        safeguard, per ver3 optimizer.py:631-637."""
        xs = _parabolic_interp([self.a, self.b, self.c],
                               [self.Ja, self.Jb, self.Jc])
        a, b, c = self.a, self.b, self.c
        bad = (
            math.isnan(xs)
            or xs > c
            or xs < a
            or abs(math.log10(max((c - b), self.EPS) / max((b - a), self.EPS))) > 1.0
        )
        if bad:
            if b > 0.5 * (a + c):
                xs = b - self.CR * (b - a)
            else:
                xs = b + self.CR * (c - b)
        return xs

    # -- helpers -------------------------------------------------------------

    def _accept_b(self) -> Tuple[str, float, Optional[float]]:
        self.pending_alpha = self.b
        return (LSDecision.ACCEPT, self.b, self.Jb)

    def _maybe_continue(self, next_alpha: float) -> Tuple[str, float, Optional[float]]:
        if next_alpha < self.min_step:
            self.pending_alpha = self.min_step
            return (LSDecision.CRASHED, self.min_step, None)
        self.pending_alpha = next_alpha
        return (LSDecision.CONTINUE, next_alpha, None)


# ---------------------------------------------------------------------------
# Armijo backtracking: stub, registered for the dispatch table
# ---------------------------------------------------------------------------


class ArmijoLineSearch(LineSearch):
    """Placeholder. Wire-in plan in PLAN_case_a.md ("Out of scope").

    Sufficient-decrease backtracking is straightforward to implement when
    measurement justifies it (in particular, when bracketing tends to accept
    at the first trial). Not built today because golden-section guarantees
    the L-BFGS curvature condition <s, y> > 0 by construction; Armijo does
    not, and starvation of the L-BFGS history is a real failure mode.
    """

    TYPE_TAG = "armijo"

    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "ArmijoLineSearch is not implemented yet. See PLAN_case_a.md "
            "for the design. Use line_search.type: golden_section in the YAML."
        )

    def first_trial(self, f_baseline: float, gT_d: float) -> float:  # pragma: no cover
        raise NotImplementedError

    def record_and_next(self, last_alpha, last_J):  # pragma: no cover
        raise NotImplementedError

    def to_dict(self):  # pragma: no cover
        raise NotImplementedError

    @classmethod
    def from_dict(cls, d, cfg):  # pragma: no cover
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


_REGISTRY = {
    GoldenSectionLineSearch.TYPE_TAG: GoldenSectionLineSearch,
    ArmijoLineSearch.TYPE_TAG: ArmijoLineSearch,
}


def make_line_search(cfg, ls_type: str) -> LineSearch:
    """Construct a fresh LS instance with knobs read from cfg.

    Reuses the existing optimization.line_search.* YAML keys consumed by
    tao_main (optim.py: more-thuente, etc) -- initial_step_size, max_step,
    min_step, bracket_tol, max_funcs -- so a YAML written for the tao
    driver works for the manual driver with only optimization.type and
    optimization.line_search.type changed.
    """
    initial_step = cfg.getInput(
        ["optimization", "line_search", "initial_step_size"], fallback=1.0)
    max_step = cfg.getInput(
        ["optimization", "line_search", "max_step"], fallback=1.0e+15)
    min_step = cfg.getInput(
        ["optimization", "line_search", "min_step"], fallback=1.0e-20)
    bracket_tol = cfg.getInput(
        ["optimization", "line_search", "bracket_tol"], fallback=1.0e-4)
    max_funcs = cfg.getInput(
        ["optimization", "line_search", "max_funcs"], fallback=30)
    if ls_type not in _REGISTRY:
        raise SystemExit(
            f"optimization.line_search.type = {ls_type!r} not recognized; "
            f"available: {sorted(_REGISTRY)}"
        )
    return _REGISTRY[ls_type](
        initial_step=initial_step,
        max_step=max_step,
        min_step=min_step,
        bracket_tol=bracket_tol,
        max_funcs=max_funcs,
    )


def load_line_search(d: Dict[str, Any], cfg) -> LineSearch:
    """Reconstruct an LS instance from a serialized state blob.

    Reads the "type" key to dispatch to the right subclass's `from_dict`.
    Raises if the type is unknown.
    """
    ls_type = d.get("type")
    if ls_type not in _REGISTRY:
        raise SystemExit(
            f"lbfgs_state.json line_search.type = {ls_type!r} not recognized"
        )
    return _REGISTRY[ls_type].from_dict(d, cfg)
