#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
import numpy.random as random

import plot3dnasa as p3d

EPSILON = np.finfo(float).eps

class InstabilityMode:

    def __init__(self, n=0, omega=0., theta=1., mach_number=0.,
                 r_nozzle=0.5, gamma=1.4):
        self.n = n
        self.omega = omega
        self.theta = theta
        self.mach_number = mach_number
        self.r_nozzle = r_nozzle
        self.gamma = gamma

    def update(self):
        self._temperature_ratio = 1. / (1. + 0.5 * (self.gamma - 1.) *
                                        self.mach_number ** 2)
        self._u_j = self.mach_number * np.sqrt(self._temperature_ratio)
        return self

    def set_domain(self, r, r_lower, r_match, r_upper):
        self._r = r
        i_lower = np.argmin(np.abs(self._r - r_lower))
        i_match = np.argmin(np.abs(self._r - r_match))
        i_upper = np.argmin(np.abs(self._r - r_upper))
        self._num_steps_inner = i_match - i_lower
        self._num_steps_outer = i_upper - i_match
        self._r_lower = self._r[i_lower]
        self._r_match = self._r[i_match]
        self._r_upper = self._r[i_upper]
        return self

    def mean_flow(self, r):
        u = 0.5 * (1. + np.tanh(0.25 / self.theta * (
            self.r_nozzle / (r + EPSILON) - r / self.r_nozzle)))
        u_dot = -0.125 / self.theta * (1. - np.tanh(0.25 / self.theta * (
            self.r_nozzle / (r + EPSILON) - r / self.r_nozzle)) ** 2) * (
                self.r_nozzle / (r + EPSILON) ** 2 + 1. / self.r_nozzle)
        # u, \rho, a, du/dr, d\rho/dr
        mean_flow = np.array([np.empty_like(r) for i in range(5)])
        mean_flow[0] = self.mach_number * np.sqrt(self._temperature_ratio) * u
        mean_flow[1] = 1. / self._temperature_ratio / (
            0.5 * (self.gamma - 1.) * u * (1. - u) * self.mach_number ** 2 +
            u + (1. - u) / self._temperature_ratio)
        mean_flow[2] = np.sqrt(1. / mean_flow[1])
        mean_flow[3] = self.mach_number * np.sqrt(self._temperature_ratio) * \
                       u_dot
        mean_flow[4] = -self._temperature_ratio * mean_flow[1] ** 2 * u_dot * (
            0.5 * (self.gamma - 1.) * (1. - 2. * u) * self.mach_number ** 2 +
            1. - 1. / self._temperature_ratio)
        return mean_flow

    def rayleigh_rhs(self, r, p_hat, alpha):
        u, rho, a, u_r, rho_r = self.mean_flow(r)
        return np.array([p_hat[1],
                         - (1. / r - rho_r / rho + 2. * alpha /
                            (self.omega - alpha * u) * u_r) * p_hat[1] -
                         ((self.omega - alpha * u) ** 2 / a ** 2 -
                          self.n ** 2 / r ** 2 - alpha ** 2) * p_hat[0]])

    def inner_solution(self, r, alpha):
        import scipy.special
        eta = np.sqrt((self.omega - alpha * self._u_j) ** 2 /
                      self._temperature_ratio - alpha ** 2)
        return np.array([scipy.special.jn(self.n, eta * r),
                         eta * 0.5 * (scipy.special.jn(self.n - 1, eta * r) -
                                      scipy.special.jn(self.n + 1, eta * r))])

    def outer_solution(self, r, alpha):
        import scipy.special
        eta = np.sqrt(self.omega ** 2 - alpha ** 2)
        return np.array([scipy.special.hankel1(self.n,
                                               np.sign(eta.imag) * eta * r),
                         np.sign(eta.imag) * eta * 0.5 * (
                             scipy.special.hankel1(self.n - 1,
                                                   np.sign(eta.imag) *
                                                   eta * r) -
                             scipy.special.hankel1(self.n + 1,
                                                   np.sign(eta.imag) *
                                                   eta * r))])

    def newton_raphson_func(self, alpha):
        from scipy.integrate import ode
        integrator = ode(self.rayleigh_rhs)
        integrator.set_integrator('zvode', method='bdf',
                                  order=15).set_f_params(alpha)
        # Inner solution
        integrator.set_initial_value(self.inner_solution(self._r_lower, alpha),
                                     self._r_lower)
        dr = (self._r_match - self._r_lower) / self._num_steps_inner
        while integrator.successful() and \
              integrator.t < self._r_match - 0.5 * dr:
            integrator.integrate(integrator.t + dr)
        inner_solution = integrator.y
        # Outer solution
        integrator.set_initial_value(self.outer_solution(self._r_upper, alpha),
                                     self._r_upper)
        dr = (self._r_upper - self._r_match) / self._num_steps_outer
        while integrator.successful() and \
              integrator.t > self._r_match + 0.5 * dr:
            integrator.integrate(integrator.t - dr)
        outer_solution = integrator.y
        if abs(outer_solution[0] / inner_solution[0]) < 1e-6:
            return inner_solution[1] - outer_solution[1] * \
                inner_solution[0] / outer_solution[0]
        return inner_solution[1] * outer_solution[0] / inner_solution[0] - \
            outer_solution[1]

    def find_eigenvalue(self, initial_guess, max_iterations=200,
                       tolerance=1e-10):
        x = initial_guess
        y = self.newton_raphson_func(x)
        it = 1
        while abs(y) > tolerance and it <= max_iterations: # Newton-Raphson
            dx = random.uniform(1e-10, 1e-8) + 1j * \
                 random.uniform(1e-10, 1e-8) # step size for derivative
            x = x - dx * y / (self.newton_raphson_func(x + dx) - y)
            y = self.newton_raphson_func(x)
            it = it + 1
        if it >= max_iterations:
            print "Newton-Raphson failed: initial guess = ", initial_guess, \
                ", current guess = ", x, ", current function = ", abs(y)
        self.alpha = x.real - 1.j * abs(x.imag)
        return self

    def find_eigenfunction(self, tolerance=1e-8):
        from scipy.integrate import ode
        from scipy.interpolate import interp1d
        integrator = ode(self.rayleigh_rhs)
        integrator.set_integrator('zvode', method='bdf',
                                  order=15).set_f_params(self.alpha)
        p_hat = np.array([np.empty_like(self._r) for i in range(2)],
                         dtype='complex128')
        i_lower = np.argmin(np.abs(self._r - self._r_lower))
        p_hat[:,:i_lower] = self.inner_solution(self._r[:i_lower], self.alpha)
        i = i_lower
        integrator.set_initial_value(self.inner_solution(
            self._r_lower, self.alpha), self._r_lower)
        dr = (self._r_match - self._r_lower) / self._num_steps_inner
        while integrator.successful() and \
              integrator.t < self._r_match - 0.5 * dr:
            p_hat[:,i] = integrator.y
            integrator.integrate(integrator.t + dr)
            i += 1
        p_hat[:,i] = integrator.y
        i_upper = np.argmin(np.abs(self._r - self._r_upper))
        p_hat[:,i_upper:] = self.outer_solution(self._r[i_upper:], self.alpha)
        i = i_upper
        integrator.set_initial_value(self.outer_solution(
            self._r_upper, self.alpha), self._r_upper)
        dr = (self._r_upper - self._r_match) / self._num_steps_outer
        while integrator.successful() and \
              integrator.t > self._r_match + 0.5 * dr:
            p_hat[:,i] = integrator.y
            integrator.integrate(integrator.t - dr)
            i -= 1
        outer_solution = integrator.y
        scaling_factor = integrator.y[0] / p_hat[0,i]
        scaling_factor_inverse = 1. / scaling_factor
        if abs(scaling_factor) < 1e-6:
            p_hat[:,i+1:] *= scaling_factor_inverse
        else:
            p_hat[:,:i+1] *= scaling_factor
        assert abs(p_hat[1,i] - integrator.y[1]) < tolerance
        self.p_hat_dot = interp1d(self._r, p_hat[1,:], kind='linear')
        self.p_hat = interp1d(self._r, p_hat[0,:], kind='linear')
        return self

    def get_eigenmode(self, r):
        u, rho, a, u_r, rho_r = self.mean_flow(r)
        p_hat = self.p_hat(r)
        p_hat_r = self.p_hat_dot(r)
        omega = self.omega - self.alpha * u
        u_hat_r = 1. / (1.j * rho * omega) * p_hat_r
        u_hat_x = self.alpha * p_hat / (rho * omega) + \
                  u_hat_r / (1.j * omega) * u_r
        rho_hat = p_hat / a ** 2 + u_hat_r / (1.j * omega) * rho_r
        return np.array([rho_hat,
                         rho * u_hat_r,
                         self.n * p_hat / (r * omega + EPSILON),
                         rho_hat * u + rho * u_hat_x,
                         p_hat / (self.gamma - 1.) + 0.5 * rho_hat * u ** 2 +
                         rho * u * u_hat_x])

def mesh_segment(n, x_min, x_max, dx_min, dx_max=None):
    """ If dx_max is None, returns a one-sided stretched mesh from x_min
    to x_max with minimum spacing dx_min at x_min. Otherwise, returns a
    two-sided stretched mesh with spacing dx_min and dx_max at x_min and
    x_max, respectively.
    """
    from scipy.optimize import curve_fit
    if dx_max is None:
        f = lambda x, sigma: np.sinh(sigma * x) / np.sinh(sigma)
        sigma = fsolve(lambda x: ((x_max - x_min) * f(1. / (n - 1.), x) -
                                  dx_min) ** 2, 2.)
        return x_min + (x_max - x_min) * f(np.linspace(0., 1., n), sigma)
    f = lambda xi, a, b, c, d: a + b * np.tanh(c * (xi - d))
    a, b, c, d = curve_fit(f, np.array([0., 1. / (n - 1), 1. -
                                        1. / (n - 1), 1.]),
                           np.array([x_min, x_min + dx_min,
                                     x_max - dx_max, x_max]))[0]
    return f(np.linspace(0., 1., n), a, b, c, d)

def axial_coordinate(z_min=-9., z_max=34., z_physical=24., num_axial=512,
                     num_inflow_sponge=48, num_outflow_sponge=28,
                     num_potential_core=240, potential_core_length=8.,
                     dz_exit=0.04, dz_min=0.01, dz_physical_outflow=0.22):
    z = -mesh_segment(num_inflow_sponge + 1, 0., -z_min, dz_exit)[::-1]
    z = np.append(z[:-1], mesh_segment(num_potential_core + 1, 0.,
                                       potential_core_length, dz_exit, dz_min))
    z = np.append(z[:-1], mesh_segment(
        num_axial - (num_inflow_sponge + num_outflow_sponge +
                     num_potential_core), potential_core_length, z_physical,
        dz_min, dz_physical_outflow))
    z = np.append(z[:-1], mesh_segment(num_outflow_sponge + 1, z_physical,
                                       z_max, dz_physical_outflow))
    return z

def dc_func(x, n, num_const_dr, dr_min, a_inner, sigma):
    s = np.linspace(0., 1., n)
    dc = np.empty_like(s)
    dc[:num_const_dr] = x * (0.5 - 0.5 * a_inner) / (n + 1)
    dc[-num_const_dr:] = dr_min
    y = -np.tanh((s[num_const_dr:-num_const_dr] - 0.5) / sigma)
    dc[num_const_dr:-num_const_dr] = dc[-1] + (dc[0] - dc[-1]) * \
                                    (y - y.min()) / (y.max() - y.min())
    c = np.cumsum(dc) + 0.5 * a_inner - dc[0]
    return (c[-1] - 0.5) ** 2

def nozzle_quadrant(grid_size, num_const_dr=4, a_inner=0.24,
                    p_inner=1.12, dr_min=0.005):
    theta = np.linspace(0., np.pi / 2, grid_size[1])
    s = np.linspace(0., 1., grid_size[0])
    p = np.empty_like(s)
    p[:num_const_dr] = p_inner
    p[-num_const_dr:] = 2.
    p[num_const_dr:-num_const_dr] = np.linspace(p_inner, 2., grid_size[0] -
                                                2 * num_const_dr)
    sigma = 0.2
    xi = fsolve(dc_func, 1.2, args=(grid_size[0], num_const_dr,
                                    dr_min, a_inner, sigma))
    dc = np.empty_like(s)
    dc[:num_const_dr] = xi * (0.5 - 0.5 * a_inner) / (grid_size[0] + 1)
    dc[-num_const_dr:] = dr_min
    y = -np.tanh((s[num_const_dr:-num_const_dr] - 0.5) / sigma)
    dc[num_const_dr:-num_const_dr] = dc[-1] + (dc[0] - dc[-1]) * \
                                    (y - y.min()) / (y.max() - y.min())
    c = np.cumsum(dc) + 0.5 * a_inner - dc[0]
    x = np.zeros([grid_size[0], grid_size[1]])
    y = np.zeros_like(x)
    for i in range(grid_size[0]):
        r = c[i] / (np.cos(theta) ** p[i] +
                    np.sin(theta) ** p[i]) ** (1. / p[i])
        x[i,:] = r * np.cos(theta)
        y[i,:] = r * np.sin(theta)
    return x, y

def near_field_quadrant(grid_size, num_const_dr=4, dr_min=0.005, r_max=12.5):
    r = np.append(np.linspace(
        0.5, 0.5 + (num_const_dr - 1) * dr_min, num_const_dr), mesh_segment(
            grid_size[0] - num_const_dr + 2, 0.5 + (num_const_dr - 1) * dr_min,
            r_max, dr_min)[1:])
    p3d.mesh_stats(r)
    r, theta = np.meshgrid(r[1:], np.linspace(0., np.pi / 2, grid_size[1]))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x.T, y.T

def complete_rotated_blocks(g):
    n = g.get_size(0)[2]
    for i in range(2, 5):
        for k in range(n):
            g.xyz[i][:,:,k,0] = -g.xyz[i-1][:,:,k,1]
            g.xyz[i][:,:,k,1] = g.xyz[i-1][:,:,k,0]
        g.xyz[i][:,:,:,2] = g.xyz[i-1][:,:,:,2]
    return g

def elliptic_solve(x, tol=1e-4, max_iter=2000):
    y = np.zeros_like(x)
    i = 0
    while np.amax(np.abs(x - y)) > tol and i < max_iter:
        y[:,:] = x
        x[1:-1,1:-1] = (4. * x[1:-1,1:-1] + x[:-2,1:-1] + x[2:,1:-1] +
                        x[1:-1,:-2] + x[1:-1,2:]) / 8.
        i += 1
    return x

def complete_inner_block(g, smoothing_tol=1e-8):
    from scipy.ndimage.filters import uniform_filter
    n = g.get_size(0)
    g.xyz[0][0,:,0,0:2] = g.xyz[2][0,::-1,0,0:2]
    g.xyz[0][-1,:,0,0:2] = g.xyz[4][0,:,0,0:2]
    g.xyz[0][:,0,0,0:2] = g.xyz[3][0,:,0,0:2]
    g.xyz[0][:,-1,0,0:2] = g.xyz[1][0,::-1,0,0:2]
    for i in range(1, n[0] - 1):
        for j in range(1, n[1] - 1):
            g.xyz[0][i,j,0,0] = g.xyz[0][0,j,0,0] + i / (n[0] - 1.) * \
                                (g.xyz[0][-1,j,0,0] - g.xyz[0][0,j,0,0])
            g.xyz[0][i,j,0,1] = g.xyz[0][i,0,0,1] + j / (n[1] - 1.) * \
                                (g.xyz[0][i,-1,0,1] - g.xyz[0][i,0,0,1])
    g.xyz[0][:,:,0,0] = elliptic_solve(g.xyz[0][:,:,0,0], tol=smoothing_tol)
    g.xyz[0][:,:,0,1] = elliptic_solve(g.xyz[0][:,:,0,1], tol=smoothing_tol)
    for k in range(n[2]):
        g.xyz[0][:,:,k,0:2] = g.xyz[0][:,:,0,0:2]
        g.xyz[0][:,:,k,2] = g.xyz[1][0,0,k,2]
    return g

def plot_axial_spacing(z):
    import matplotlib.pyplot as plt
    ax = plt.subplot(111)
    dz = z[1:] - z[:-1]
    ax.plot(z[:-1], dz, 'ko-')
    plt.show()

def plot_radial_spacing(r):
    import matplotlib.pyplot as plt
    ax = plt.subplot(111)
    dr = r[1:] - r[:-1]
    ax.plot(r[:-1], dr, 'ko-')
    plt.show()

def plot_jacobian_continuity(g):
    import matplotlib.pyplot as plt
    from numpy.linalg import det
    ax = plt.subplot(111)
    f = p3d.compute_jacobian(g)
    for k in range(5):
        x = g.xyz[k][:,:,0,0]
        y = g.xyz[k][:,:,0,1]
        z = np.empty_like(x)
        for i in range(z.shape[0]):
            for j in range(z.shape[1]):
                z[i,j] = 1. / det(np.reshape(f.f[k][i,j,0,:], [3, 3]))
        if k == 0:
            vmin = z.min()
            vmax = z.max()
        c = ax.contour(x, y, z, levels=np.linspace(vmin, vmax, 31))
    plt.colorbar(c)
    ax.set_xlim([-0.5, 0.5])
    ax.set_ylim([-0.5, 0.5])
    ax.set_aspect('equal')
    plt.show()

def grid(num_radial_nozzle, num_radial_near_field, num_azimuthal,
         num_axial=512, a_inner=0.24, p_inner=1.12, dr_min=0.005):
    assert num_azimuthal % 4 == 0
    num_radial = num_radial_nozzle + num_radial_near_field
    g = p3d.Grid().set_size([
        [num_azimuthal / 4, num_azimuthal / 4, num_axial],
        [num_radial, num_azimuthal / 4, num_axial],
        [num_radial, num_azimuthal / 4, num_axial],
        [num_radial, num_azimuthal / 4, num_axial],
        [num_radial, num_azimuthal / 4, num_axial]], True)
    x, y = nozzle_quadrant([num_radial_nozzle, num_azimuthal / 4],
                           a_inner=a_inner, p_inner=p_inner, dr_min=dr_min)
    for k in range(num_axial):
        g.xyz[1][:num_radial_nozzle,:,k,0] = x
        g.xyz[1][:num_radial_nozzle,:,k,1] = y
    x, y = near_field_quadrant([num_radial_near_field, num_azimuthal / 4],
                               dr_min=dr_min)
    for k in range(num_axial):
        g.xyz[1][num_radial_nozzle:,:,k,0] = x
        g.xyz[1][num_radial_nozzle:,:,k,1] = y
    if num_axial > 1:
        z = axial_coordinate(num_axial=num_axial)
        p3d.mesh_stats(z)
        # plot_axial_spacing(z)
        for k in range(num_axial):
            g.xyz[1][:,:,k,2] = z[k]
    complete_rotated_blocks(g)
    for i in range(1, 5):
        x = np.copy(g.xyz[i][:,:,:,0])
        g.xyz[i][:,:,:,0] = (g.xyz[i][:,:,:,0] -
                             g.xyz[i][:,:,:,1]) / np.sqrt(2.)
        g.xyz[i][:,:,:,1] = (x + g.xyz[i][:,:,:,1]) / np.sqrt(2.)
    complete_inner_block(g)
    r = np.append(g.xyz[0][:,num_azimuthal/8,0,0],
                  g.xyz[4][1:,num_azimuthal/8,0,0])
    # plot_radial_spacing(r)
    # plot_jacobian_continuity(g)
    return g

def target_state(g, mach_number=1.3, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    temperature_ratio = 1. / (1. + 0.5 * (gamma - 1.) * mach_number ** 2)
    T_inf = 1./ (gamma - 1.)
    u_j = mach_number * np.sqrt(temperature_ratio)
    for i, xyz in enumerate(g.xyz):
        z = xyz[0,0,:,2]
        r = np.sqrt(xyz[:,:,0,0] ** 2 + xyz[:,:,0,1] ** 2)
        condlist = [z <= 0., np.logical_and(z > 0., z < 24.), z >= 24.]
        r0 = np.select(condlist, [0.5 + np.zeros_like(z), 0.04125 * z + 0.5,
                                  0.075 * z - 0.31])
        theta = np.select(condlist, [0.04 + np.zeros_like(z),
                                     0.46 * z / 24. + 0.04, 0.5])
        u = np.where(z <= 2.65, u_j + np.zeros_like(z),
                     u_j * np.exp(-(z - 2.65) / 25.))
        rho = np.where(z <= 0., 1. / temperature_ratio + np.zeros_like(z),
                       1. + (1. / temperature_ratio - 1.) * np.exp(-0.078 * z))
        for k in range(z.size):
            s.q[i][:,:,k,3] = 0.5 * u[k] * (1. + np.tanh(
                0.25 / theta[k] * (r0[k] / r - r / r0[k])))
            s.q[i][:,:,k,0] = rho[k] / (0.5 * (gamma - 1.) * s.q[i][:,:,k,3] /
                                        u[k] * (1. - s.q[i][:,:,k,3] / u[k]) *
                                        rho[k] * u[k] ** 2 + s.q[i][:,:,k,3] /
                                        u[k] + rho[k] *
                                        (1. - s.q[i][:,:,k,3] / u[k]))
    return s.fromprimitive(gamma)

def initial_condition(g, mach_number=1.3, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    temperature_ratio = 1. / (1. + 0.5 * (gamma - 1.) * mach_number ** 2)
    u_j = mach_number * np.sqrt(temperature_ratio)
    for i, xyz in enumerate(g.xyz):
        z = xyz[0,0,:,2]
        r = np.sqrt(xyz[:,:,0,0] ** 2 + xyz[:,:,0,1] ** 2) / 0.5
        condlist = [z <= 0., np.logical_and(z > 0., z < 24.), z >= 24.]
        theta = 0.04 + np.where(z <= 0., np.zeros_like(z), 0.46 * z / 34.)
        for k in range(z.size):
            s.q[i][:,:,k,3] = 0.5 * u_j * (1. + np.tanh(
                0.25 / theta[k] * (1. / r - r)))
        s.q[i][:,:,:,0] = 1. / (0.5 * (gamma - 1.) * s.q[i][:,:,:,3] / u_j *
                                (1. - s.q[i][:,:,:,3] / u_j) *
                                mach_number ** 2 + s.q[i][:,:,:,3] / u_j +
                                (1. - s.q[i][:,:,:,3] / u_j) /
                                temperature_ratio) / temperature_ratio
    return s.fromprimitive(gamma)

def plot_eigenvalues(St, alpha, u_j, theta_j):
    import matplotlib.pyplot as plt
    ax = plt.subplot(121)
    for i in range(alpha.shape[0]):
        ax.plot(St, -alpha[i,:].imag * theta_j, 'rs-', mec='r', mfc='w')
    ax.set_xlim([0.35, 0.95])
    ax.set_ylim([0.05, 0.11])
    ax = plt.subplot(122)
    omega = 2. * np.pi * u_j * St
    for i in range(alpha.shape[0]):
        ax.plot(St, omega / alpha[i,:].real / u_j, 'rs-', mec='r', mfc='w')
    ax.set_xlim([0.35, 0.95])
    ax.set_ylim([0.5, 0.75])
    plt.show()

def eigenmodes(mach_number=1.3, theta_j=0.04, show_progress=True):
    u_j = InstabilityMode(mach_number=mach_number,
                          theta=theta_j).update().mean_flow(EPSILON)[0]
    n = np.arange(1, 6)
    St = np.array([0.43, 0.51, 0.61, 0.69, 0.74, 0.88])
    alpha = np.empty([n.size, St.size], dtype='complex128')
    r_lower = 0.025
    r_match = 0.25
    r_upper = 3.
    r = np.linspace(0., r_lower, 5001)
    r = np.append(r[:-1], np.linspace(r_lower, r_match, 5001))
    r = np.append(r[:-1], np.linspace(r_match, r_upper, 5001))
    r = np.append(r[:-1], np.linspace(r_upper, 15., 5001))
    # Initialize eigenmodes.
    initial_guess = 1. - 1.j
    modes = [[InstabilityMode() for n_ in n] for St_ in St]
    if show_progress:
        from progressbar import ProgressBar, Percentage, Bar, ETA
        print 'Computing eigenmodes:'
        p = ProgressBar(widgets = [Percentage(), ' ',
                                   Bar('=', left = '[', right = ']'), ' ',
                                   ETA()], maxval = alpha.size).start()
    for j, St_ in enumerate(St):
        guess = initial_guess
        for i, n_ in enumerate(n):
            modes[j][i] = InstabilityMode(
                mach_number=mach_number, theta=theta_j, n=n_,
                omega=2. * np.pi * u_j * St_).update().set_domain(
                    r, r_lower, r_match, r_upper)
            modes[j][i].find_eigenvalue(guess).find_eigenfunction()
            guess = alpha[i,j] = modes[j][i].alpha
            if i == 0:
                initial_guess = modes[j][i].alpha
            if show_progress:
                p.update(i + n.size * j)
    if show_progress:
        p.finish()
    # plot_eigenvalues(St, alpha, u_j, theta_j)
    return modes

def extract_inflow(g, num_inflow_sponge=49):
    n = g.get_size()
    n[:,2] = num_inflow_sponge
    gi = p3d.Grid().set_size(n, True)
    for i, xyz in enumerate(g.xyz):
        gi.xyz[i] = xyz[:,:,:num_inflow_sponge,:]
    return gi

def inflow_perturbations(g, modes):
    phi_p = 2. * np.pi * random.rand(len(modes))
    phi_n = 2. * np.pi * random.rand(len(modes))
    sr = p3d.Solution().copy_from(g)
    si = p3d.Solution().copy_from(g)
    sr._format.aux_header[1] = si._format.aux_header[1] = modes[0].omega
    for i in range(g.nblocks):
        r = np.sqrt(g.xyz[i][:,:,0,0] ** 2 + g.xyz[i][:,:,0,1] ** 2)
        theta = np.arctan2(g.xyz[i][:,:,0,1], g.xyz[i][:,:,0,0])
        z = g.xyz[i][0,0,:,2]
        si.q[i].fill(0.)
        sr.q[i].fill(0.)
        for j, mode in enumerate(modes):
            q = mode.get_eigenmode(r)
            n = mode.n
            for l in range(q.shape[0]):
                if l == 2:
                    q[l] *= np.exp(1.j * n * theta + phi_p[j]) - \
                            np.exp(-1.j * n * theta + phi_n[j])
                else:
                    q[l] *= np.exp(1.j * n * theta + phi_p[j]) + \
                            np.exp(-1.j * n * theta + phi_n[j])
            u = q[1] * np.cos(theta) - q[2] * np.sin(theta)
            v = q[1] * np.sin(theta) + q[2] * np.cos(theta)
            q[1] = u
            q[2] = v
            for l in range(q.shape[0]):
                for k, z_ in enumerate(z):
                    sr.q[i][:,:,k,l] += np.real(
                        q[l] * np.exp(1.j * mode.alpha * z_))
                    si.q[i][:,:,k,l] += np.imag(
                        q[l] * np.exp(1.j * mode.alpha * z_))
    return sr, si

def target_mollifier(g):
    z_min =  0.
    z_max = 24.
    r_min =  7.
    r_max =  9.
    f = p3d.Function().copy_from(g)
    z = g.xyz[0][0,0,:,2]
    n = f.get_size()
    block_code = ['IB', 'E', 'N', 'W', 'S']
    for i, fi in enumerate(f.f):
        r = np.sqrt(g.xyz[i][:,:,0,0] ** 2 + g.xyz[i][:,:,0,1] ** 2)
        if r.max() < r_min or r.min() > r_max:
            fi.fill(0.)
            continue
        fi.fill(1.)
        for k in range(n[i][1]):
            for j in range(n[i][0]):
                fi[j,k,:,0] *= p3d.tanh_support(z, z_min, z_max, 40., 0.2)
            for j in range(n[i][2]):
                fi[:,k,j,0] *= p3d.cubic_bspline_support(r[:,k], r_min, r_max)
        kmin, kmax = p3d.find_extents(z, z_min, z_max)
        imin, imax = p3d.find_extents(np.mean(r, axis=1), r_min, r_max)
        print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
            'targetRegion.' + block_code[i], 'COST_TARGET',
            i + 1, 0, imin, imax, 1, -1, kmin, kmax)
    return f

def control_mollifier(g):
    z_min =   1.
    z_max =   3.
    r_min = 0.3
    r_max = 0.7
    f = p3d.Function().copy_from(g)
    z = g.xyz[0][0,0,:,2]
    n = f.get_size()
    block_code = ['IB', 'E', 'N', 'W', 'S']
    for i, fi in enumerate(f.f):
        r = np.sqrt(g.xyz[i][:,:,0,0] ** 2 + g.xyz[i][:,:,0,1] ** 2)
        if r.max() < r_min or r.min() > r_max:
            fi.fill(0.)
            continue
        fi.fill(1.)
        for k in range(n[i][1]):
            for j in range(n[i][0]):
                fi[j,k,:,0] *= p3d.tanh_support(z, z_min, z_max, 16., 0.2)
            for j in range(n[i][2]):
                fi[:,k,j,0] *= p3d.tanh_support(r[:,k], r_min, r_max, 20., 0.2)
        kmin, kmax = p3d.find_extents(z, z_min, z_max)
        imin = min(np.where(r[:,i] >= r_min)[0][0]
                   for i in range(n[i][1]))
        imax = max(np.where(r[:,i] <= r_max)[0][-1] + 2
                   for i in range(n[i][1]))
        print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
            'controlRegion.' + block_code[i], 'ACTUATOR',
            i + 1, 0, imin, imax, 1, -1, kmin, kmax)
    return f

def mean_pressure(s):
    s_prim = s.toprimitive()
    f = p3d.Function().copy_from(s_prim)
    print (np.shape(f.f))
    for i, fi in enumerate(f.f):
        print (i)
        print (fi.shape)
        print (np.min(s.q[i][:,:,:,4]), np.max(s.q[i][:,:,:,4]))
        #print (s.toprimitive().q[i].shape)
        fi[:,:,:,0] = s_prim.q[i][:,:,:,4]
        print (np.min(fi), np.max(fi))
    return f


if __name__ == '__main__':
    g = grid(60, 196, 132, a_inner=0.24, p_inner=1.08634735266)
    g.save('MultiblockJet.xyz')
    control_mollifier(g).save('MultiblockJet.control_mollifier.f')
    target_mollifier(g).save('MultiblockJet.target_mollifier.f')
#    target_state(g).save('MultiblockJet.target.q')
#    initial_condition(g).save('MultiblockJet.ic.q')
#    gi = extract_inflow(g)
#    gi.save('MultiblockJet.inflow.xyz')
#    modes = eigenmodes()
#    for i, mode in enumerate(modes):
#        sr, si = inflow_perturbations(gi, mode)
#        sr.save('MultiblockJet-%02d.eigenmode_real.q' % (i + 1))
#        si.save('MultiblockJet-%02d.eigenmode_imag.q' % (i + 1))
