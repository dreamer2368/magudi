"""Tests for magudi_utils.fwhsolver reference-flow generators.

Establishes the analytic (rho, u, p) acoustic field of the monopole and
dipole sources from Kim (2012) thesis section 5.3, which will be fed as
reference input into FWHSolver in a follow-on round-trip test.

Independent physical verification of the generated field (directivity /
1/r scaling / Euler residual) is deliberately postponed.
"""
import numpy as np
import pytest

from magudi_utils import fwhsolver as fwh


GAMMA = 1.4
P_INF = 1.0 / GAMMA  # ambient nondim pressure in Kim / Vishnampet units


def test_monopole_flow_shapes_and_ambient():
    """With A0=0 the analytic field is exactly the quiescent ambient state,
    and the return shapes follow the documented broadcasting contract."""
    x = np.array([0., 0., 10.])
    y = np.zeros(3)
    flow = fwh.monopole_flow(x, y, t=0., A0=0., omega=1.0)
    assert flow['rho'].shape == ()
    assert flow['u'].shape == (3,)
    assert flow['p'].shape == ()
    assert flow['rho'] == pytest.approx(1.0)
    np.testing.assert_allclose(flow['u'], 0.0)
    assert flow['p'] == pytest.approx(P_INF)


def test_dipole_flow_directivity_signs():
    """A +y-aligned dipole radiates zero p' in the equatorial plane and
    opposite-sign p' at the +y and -y poles at the same instant."""
    y = np.zeros(3)
    r = 10.
    t = 0.7
    kw = dict(A0=1e-4, omega=2.0)
    p_plus  = fwh.dipole_flow(np.array([0.,  r, 0.]), y, t, **kw)['p'] - P_INF
    p_minus = fwh.dipole_flow(np.array([0., -r, 0.]), y, t, **kw)['p'] - P_INF
    p_eq    = fwh.dipole_flow(np.array([r,  0., 0.]), y, t, **kw)['p'] - P_INF
    assert p_plus * p_minus < 0.
    assert abs(p_eq) < 1e-14


def test_get_monopole_shape_and_ambient_offset():
    """The FWH-surface adapter returns a Fortran-ordered [n1, n0, 5, size]
    chunk and, with A0=0, holds the quiescent ambient state."""
    n1, n0, size = 5, 7, 3
    theta = np.linspace(0., 2. * np.pi, n1, endpoint=False)
    xax = np.linspace(-2., 2., n0)
    R = 3.
    xyz = np.stack(np.broadcast_arrays(
        xax[None, :],
        R * np.cos(theta)[:, None],
        R * np.sin(theta)[:, None],
    ), axis=-1)  # (n1, n0, 3)
    q = fwh.get_monopole(offset=0, size=size, xyz=xyz, dt=0.1,
                         y=np.zeros(3), A0=0., omega=1.0)
    assert q.shape == (n1, n0, 5, size)
    assert q.flags['F_CONTIGUOUS']
    np.testing.assert_allclose(q[:, :, 0, :], 1.0)
    np.testing.assert_allclose(q[:, :, 1:4, :], 0.0, atol=1e-15)
    np.testing.assert_allclose(q[:, :, 4, :], P_INF / (GAMMA - 1.))
