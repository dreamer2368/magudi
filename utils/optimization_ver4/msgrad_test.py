"""Multi-segment gradient accuracy test for msforward/msadjoint on OneDWave.

Taylor test: for x0 in the (control_forcing, ic_1, ..., ic_{Nsplit-1}) parameter
space (ic_0 is fixed, not an optimization variable) and g = grad J(x0), the
forward-difference slope

    s(h) = ( J(x0 + h*g) - J(x0) ) / h

converges to <g, g>_M at rate O(h). We get J0 from one msforward call, <g, g>_M
from one msadjoint call (it writes the value to <prefix>.adjoint_run.txt), and
the gradient pieces (control-forcing .dat per patch + ic.adjoint.q per segment
k>=1) from the same msadjoint call. Then for each h we perturb the parameter
files with zaxpy / qfile_zaxpy, re-run msforward, and check the slope.

Pinned to np=1: the actuator .dat files are MPI_File slabs concatenated;
np=1 gives a single contiguous real64 buffer, matching what msforward and
control_space_norm wrote.

Run from the staged OneDWave_msgrad/ directory:
    python3 msgrad_test.py msgrad.yml
"""
import argparse
import math
import shutil
import subprocess
import sys

import yaml


def _run(cmd):
    """Run a binary silently; on failure echo its output and re-raise."""
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(f"FAILED ({exc.returncode}): {' '.join(exc.cmd)}\n")
        if exc.stdout:
            sys.stderr.write(exc.stdout)
        if exc.stderr:
            sys.stderr.write(exc.stderr)
        raise


def _read_scalar(path):
    with open(path) as f:
        return float(f.read().strip())


def _read_sub_adjoint(path):
    """Parse <prefix>.sub_adjoint_run.txt -> (sum_ctrl_IP, sum_ic_IP) across segments."""
    ctrl_total = 0.0
    ic_total = 0.0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            ctrl_total += float(parts[1])
            ic_total += float(parts[2])
    return ctrl_total, ic_total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="Path to msgrad.yml")
    parser.add_argument("--mode", choices=("full", "ctrl", "ic"), default="full",
                        help="full: perturb control+IC, reference = msadjoint <g,g> total. "
                             "ctrl: perturb only control, reference = sum of "
                             "control_forcing_inner_product column. "
                             "ic: perturb only IC (k>=1), reference = sum of "
                             "ic_inner_product column.")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    prefix = cfg["output_prefix"]
    patches = cfg["patches"]
    nsplit = cfg["time_splitting"]["number_of_segments"]
    h_list = cfg["finite_difference"]["step_sizes"]
    order_threshold = cfg["finite_difference"].get("order_threshold", 0.5)
    max_bad_orders = cfg["finite_difference"].get("max_bad_orders", 3)

    ic_q     = lambda k: f"{prefix}-{k}.ic.q"
    ic_base  = lambda k: f"{prefix}-{k}.ic.base.q"
    ic_grad  = lambda k: f"{prefix}-{k}.ic.adjoint.q"
    cf_dat   = lambda p: f"{prefix}.control_forcing_{p}.dat"
    cf_base  = lambda p: f"{prefix}.control_forcing_{p}.base.dat"
    cf_grad  = lambda p: f"{prefix}.gradient_{p}.dat"

    # Snapshot base point (all k incl. 0 for restore; only k>=1 perturbed).
    for k in range(nsplit):
        shutil.copy(ic_q(k), ic_base(k))
    for p in patches:
        shutil.copy(cf_dat(p), cf_base(p))

    try:
        _run(["./msforward", "--input", "magudi.inp", "--output", "J0.txt"])
        j0 = _read_scalar("J0.txt")
        _run(["./msadjoint", "--input", "magudi.inp", "--output", "gg.txt"])

        # Reference inner product depends on mode.
        ctrl_sum, ic_sum = _read_sub_adjoint(f"{prefix}.sub_adjoint_run.txt")
        gg_total = _read_scalar("gg.txt")
        if args.mode == "full":
            gg = gg_total
        elif args.mode == "ctrl":
            gg = ctrl_sum
        else:  # ic
            gg = ic_sum

        print()
        print(f"Mode: {args.mode}")
        print(f"Baseline   J0    = {j0: .12e}")
        print(f"           <g,g> total (ctrl + ic) = {gg_total: .12e}")
        print(f"           sum control_forcing IP  = {ctrl_sum: .12e}")
        print(f"           sum ic IP               = {ic_sum: .12e}")
        print(f"           reference for this mode = {gg: .12e}")
        if gg <= 0.0:
            print(f"WARN: reference = {gg} is not positive; the FD test is meaningless.",
                  file=sys.stderr)
        print()
        print(f"{'h':>14} {'Jh':>22} {'slope':>22} {'rel_err':>14} {'order':>8}")

        results = []  # list of (h, Jh, slope, err)
        for h in h_list:
            if args.mode != "ic":
                for p in patches:
                    _run(["./zaxpy", cf_dat(p), f"{h:.17e}", cf_grad(p), cf_base(p)])
            if args.mode != "ctrl":
                for k in range(1, nsplit):
                    _run(["./qfile_zaxpy", ic_q(k), f"{h:.17e}", ic_grad(k), ic_base(k),
                          "--input", "magudi.inp"])

            _run(["./msforward", "--input", "magudi.inp", "--output", "Jh.txt"])
            jh = _read_scalar("Jh.txt")
            slope = (jh - j0) / h
            err = abs(slope - gg) / abs(gg) if gg != 0.0 else abs(slope - gg)

            if results:
                order = (math.log10(err / results[-1][3])
                         / math.log10(h / results[-1][0])) if err > 0 else float("nan")
            else:
                order = float("nan")
            print(f"{h: 14.4e} {jh: 22.14e} {slope: 22.14e} {err: 14.4e} {order: 8.3f}")
            results.append((h, jh, slope, err))

        # Find the prefix of rows up to the global error minimum — that is the
        # decreasing regime before the roundoff floor.
        errs = [r[3] for r in results]
        if not errs:
            print("FAIL: no FD steps ran", file=sys.stderr)
            return 1
        i_min = errs.index(min(errs))
        pre_floor = results[:i_min + 1]

        if len(pre_floor) < 2:
            print("FAIL: error did not decrease for any h", file=sys.stderr)
            return 1

        orders = [
            math.log10(pre_floor[i + 1][3] / pre_floor[i][3])
            / math.log10(pre_floor[i + 1][0] / pre_floor[i][0])
            for i in range(len(pre_floor) - 1)
        ]
        n_bad = sum(o < order_threshold for o in orders)

        print()
        print(f"Pre-floor orders ({len(orders)}): "
              f"{' '.join(f'{o:.3f}' for o in orders)}")
        print(f"Bad orders (< {order_threshold}): {n_bad} / {len(orders)} "
              f"(allowed: {max_bad_orders})")

        if n_bad <= max_bad_orders:
            print("PASS: gradient is discrete-exact within FD precision")
            return 0
        print("FAIL: gradient is not discrete-exact", file=sys.stderr)
        return 1

    finally:
        for k in range(nsplit):
            shutil.copy(ic_base(k), ic_q(k))
        for p in patches:
            shutil.copy(cf_base(p), cf_dat(p))


if __name__ == "__main__":
    sys.exit(main())
