#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d
import tempfile
import subprocess
import os
import os.path
import random

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

def draw_circular_arc(f, label, p1, p2, c):
    print >>f, 'gg::conBegin'
    print >>f, 'gg::segBegin -type CIRCULAR_ARC'
    print >>f, 'gg::segAddControlPt {%f %f %f}' % (p1[0], p1[1], p1[2])
    print >>f, 'gg::segAddControlPt {%f %f %f}' % (p2[0], p2[1], p2[2])
    print >>f, 'gg::segAddControlPt -alternate CENTER {%f %f %f}' % \
        (c[0], c[1], c[2])
    print >>f, 'gg::segEnd'
    print >>f, 'set %s [gg::conEnd]' % (label)
    return None

def draw_line3d(f, label, p1, p2):
    print >>f, 'gg::conBegin'
    print >>f, 'gg::segBegin -type 3D_LINE'
    print >>f, 'gg::segAddControlPt {%f %f %f}' % (p1[0], p1[1], p1[2])
    print >>f, 'gg::segAddControlPt {%f %f %f}' % (p2[0], p2[1], p2[2])
    print >>f, 'gg::segEnd'
    print >>f, 'set %s [gg::conEnd]' % (label)
    return None

def rotate_copy_connector(f, label, p1, p2, angle, new_label=None,
                          replace=False):
    if replace is True:
        print >>f, 'gg::conTransformBegin $%s -maintain_linkage' % (label)
    else:
        print >>f, 'gg::conCopyBegin $%s' % (label)
    print >>f, 'gg::xformRotate {%f %f %f} {%f %f %f} %f' % \
        (p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], angle)
    if replace is True:
        print >>f, 'gg::conTransformEnd'
    else:
        print >>f, 'set %s [gg::conCopyEnd]' % (new_label)
    return None

def rotate_copy_domain(f, label, p1, p2, angle, new_label=None,
                       replace=False):
    if replace is True:
        print >>f, 'gg::domTransformBegin $%s -maintain_linkage' % (label)
    else:
        print >>f, 'gg::domCopyBegin $%s' % (label)
    print >>f, 'gg::xformRotate {%f %f %f} {%f %f %f} %f' % \
        (p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], angle)
    if replace is True:
        print >>f, 'gg::domTransformEnd'
    else:
        print >>f, 'set %s [gg::domCopyEnd]' % (new_label)
    return None

def rotate_copy_block(f, label, p1, p2, angle, new_label=None,
                      replace=False):
    if replace is True:
        print >>f, 'gg::blkTransformBegin $%s -maintain_linkage' % (label)
    else:
        print >>f, 'gg::blkCopyBegin $%s' % (label)
    print >>f, 'gg::xformRotate {%f %f %f} {%f %f %f} %f' % \
        (p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], angle)
    if replace is True:
        print >>f, 'gg::blkTransformEnd'
    else:
        print >>f, 'set %s [gg::blkCopyEnd]' % (new_label)
    return None

def domain(f, edge1, edge2, edge3, edge4, label):
    print >>f, 'gg::domBegin -type STRUCTURED'
    print >>f, 'gg::edgeBegin'
    print >>f, 'gg::edgeAddCon $%s' % (edge1)
    print >>f, 'gg::edgeEnd'
    print >>f, 'gg::edgeBegin'
    print >>f, 'gg::edgeAddCon $%s' % (edge2)
    print >>f, 'gg::edgeEnd'
    print >>f, 'gg::edgeBegin'
    print >>f, 'gg::edgeAddCon $%s' % (edge3)
    print >>f, 'gg::edgeEnd'
    print >>f, 'gg::edgeBegin'
    print >>f, 'gg::edgeAddCon $%s' % (edge4)
    print >>f, 'gg::edgeEnd'
    print >>f, 'set %s [gg::domEnd]' % (label)
    return None

def write_glyph_file(filename, num_axial=1, num_radial=281, num_azimuthal=192,
                     a_inner=0.5, p_inner=1., r_nozzle=0.5, r_outer=12.5,
                     z_min=-10., z_max=34., potential_core_length=11.,
                     fraction_inside_nozzle=0.1, num_const_dr=12,
                     dr_min_fraction=0.8, dz_min_fraction=2.25,
                     interface_ds_ratio=1.):
    assert num_azimuthal % 4 == 0
    num_radial_inside_nozzle = int(np.floor(fraction_inside_nozzle *
                                            num_radial))
    dr_interface = interface_ds_ratio * a_inner / (num_azimuthal / 4)
    dr_min = dr_min_fraction * a_inner / (num_azimuthal / 4)
    dz_min = dz_min_fraction * dr_min
    # Write glyph file.
    with open(filename, 'w') as f:
        print >>f, """
        gg::memClear
        gg::defReset
        gg::tolReset
        gg::dispViewReset
        gg::dispGlide FALSE
        set cwd [file dirname [info script]]
        """
        theta = np.linspace(-np.pi / 4., np.pi / 4., 201)
        if p_inner < 2.:
            theta += np.pi / 4.
        print >>f, 'gg::dbCurveBegin -type CUBIC'
        for theta_ in theta:
            r = 0.5 * a_inner / \
                (abs(np.cos(theta_)) ** p_inner +
                 abs(np.sin(theta_)) ** p_inner) ** (1. / p_inner)
            print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % \
                (r * np.cos(theta_), r * np.sin(theta_), 0.)
        print >>f, 'set db_inner_block_E [gg::dbCurveEnd]'
        print >>f, """
        set inner_block_E [gg::conOnDBEnt $db_inner_block_E]
        """
        rotate_copy_connector(f, 'inner_block_E', [0., 0., 0.], [0., 0., 1.],
                              90., 'inner_block_N')
        rotate_copy_connector(f, 'inner_block_N', [0., 0., 0.], [0., 0., 1.],
                              90., 'inner_block_W')
        rotate_copy_connector(f, 'inner_block_W', [0., 0., 0.], [0., 0., 1.],
                              90., 'inner_block_S')
        # Dimension inner block and create domain.
        print >>f, 'gg::conDim [list $inner_block_E $inner_block_N ' \
            '$inner_block_W $inner_block_S] %i' % (num_azimuthal / 4 + 1)
        domain(f, 'inner_block_S', 'inner_block_E', 'inner_block_N',
               'inner_block_W', 'dom_inner')
        if p_inner < 2.:
            rotate_copy_domain(f, 'dom_inner', [0., 0., 0.], [0., 0., 1.],
                               -45., replace=True)
        # Create south-east connector to specify extrusion step size.
        theta_se = theta[0]
        r_se = 0.5 * a_inner / \
               (abs(np.cos(theta_se)) ** p_inner +
                abs(np.sin(theta_se)) ** p_inner) ** (1. / p_inner)
        draw_line3d(f, 'tmp_extrusion_SE',
               [r_se * np.cos(theta_se), r_se * np.sin(theta_se), 0.],
               [r_outer * np.cos(theta_se), r_outer * np.sin(theta_se), 0.])
        # Dimension south-east connector.
        print >>f, 'gg::conDim [list $tmp_extrusion_SE] %i' % (num_radial)
        print >>f, 'gg::conSetBreakPt $tmp_extrusion_SE {%f %f 0}' % \
            ((r_se + num_const_dr * dr_interface) * np.cos(theta_se),
             (r_se + num_const_dr * dr_interface) * np.sin(theta_se))
        print >>f, 'gg::conSetBreakPt $tmp_extrusion_SE {%f %f 0}' % \
            (r_nozzle * np.cos(theta_se), r_nozzle * np.sin(theta_se))
        print >>f, 'gg::conSubConDim $tmp_extrusion_SE [list %i %i %i]' % \
            (num_const_dr + 1, num_radial_inside_nozzle - num_const_dr + 1,
             num_radial - num_radial_inside_nozzle)
        print >>f, 'gg::conBeginSpacing $tmp_extrusion_SE -sub 3 %f' % (dr_min)
        print >>f, 'gg::conBeginSpacing $tmp_extrusion_SE -sub 2 %f' % \
            (dr_interface)
        print >>f, 'gg::conEndSpacing $tmp_extrusion_SE -sub 2 %f' % (dr_min)
        print >>f, 'gg::conBeginSpacing $tmp_extrusion_SE -sub 1 %f' % \
            (dr_interface)
        print >>f, 'gg::conEndSpacing $tmp_extrusion_SE -sub 1 %f' % \
            (dr_interface)
        if p_inner < 2.:
            rotate_copy_connector(f, 'tmp_extrusion_SE', [0., 0., 0.],
                                  [0., 0., 1.], -45., replace=True)
        # Extrude to create the outer blocks.
        print >>f, """
        gg::domExtrusionBegin [list $inner_block_E $inner_block_N \
        $inner_block_W $inner_block_S] -default HYPERBOLIC
        gg::domExtrusionAtt -flip
        gg::domExtrusionAtt -growth_subcon [list [list $tmp_extrusion_SE 1] \
        [list $tmp_extrusion_SE 2] [list $tmp_extrusion_SE 3]]
        """
        print >>f, 'gg::domExtrusionAtt -stop_height %f' % (r_outer)
        print >>f, 'gg::domExtrusionStep %i' % (num_radial)
        print >>f, """
        set dom_outer_blocks [gg::domExtrusionEnd]
        set dom_outer_E [lindex $dom_outer_blocks 0]
        set dom_outer_N [lindex $dom_outer_blocks 1]
        set dom_outer_W [lindex $dom_outer_blocks 2]
        set dom_outer_S [lindex $dom_outer_blocks 3]
        unset dom_outer_blocks
        """
        if num_axial > 1:
            # Move domains to axial start coordinate.
            print >>f, 'gg::domTransformBegin [list $dom_inner $dom_outer_E ' \
                '$dom_outer_N $dom_outer_W $dom_outer_S] -maintain_linkage'
            print >>f, 'gg::xformTranslate {0 0 %f}' % (z_min)
            print >>f, 'gg::domTransformEnd'
            # Create the axis connector.
            draw_line3d(f, 'axis', [0., 0., z_min], [0., 0., z_max])
            print >>f, 'gg::conDim [list $axis] %i' % (num_axial)
            print >>f, 'gg::conSetBreakPt $axis {0 0 %f}' % \
                (potential_core_length)
            print >>f, 'gg::conBreakPtSpacing $axis -breakpt 1 %f' % (dz_min)
            # Extrude blocks.
            print >>f, 'gg::blkExtrusionBegin [list $dom_inner $dom_outer_E ' \
                '$dom_outer_N $dom_outer_W $dom_outer_S] -default PATH'
            print >>f, 'gg::blkExtrusionAtt -path_connector ' \
                '[list [list $axis 1] [list $axis 2]]'
            print >>f, 'gg::blkExtrusionStep %i' % (num_axial)
            print >>f, """
            set all_blocks [gg::blkExtrusionEnd]
            set blk_inner [lindex $all_blocks 0]
            set blk_outer_E [lindex $all_blocks 1]
            set blk_outer_N [lindex $all_blocks 2]
            set blk_outer_W [lindex $all_blocks 3]
            set blk_outer_S [lindex $all_blocks 4]
            unset all_blocks
            """
            # Re-orient blocks.
            print >>f, 'gg::blkSpecifyIJK [list $blk_outer_E $blk_outer_N ' \
                '$blk_outer_W $blk_outer_S] 5 4 1'
            print >>f, 'gg::dbDelete ALL'
            # Save PLOT3D grid file.
            print >>f, 'gg::blkExport ALL \"OSUMach1.3.xyz\" -style PLOT3D ' \
                '-form UNFORMATTED -precision DOUBLE -endian NATIVE'
        else:
            print >>f, 'gg::dbDelete ALL'
            # Save PLOT3D grid file.
            print >>f, 'gg::domExport ALL \"OSUMach1.3.xyz\" -style PLOT3D ' \
                '-form UNFORMATTED -precision DOUBLE -endian NATIVE'

def grid(**kwargs):
    if not os.path.isfile('OSUMach1.3.xyz'):
        f = tempfile.NamedTemporaryFile(delete=False)
        write_glyph_file(f.name, **kwargs)
        subprocess.check_output(["gridgen", "-b", f.name])
        os.unlink(f.name)

def target_state(g, mach_number=1.3, theta_j=0.02, S=0.06,
                 potential_core_length=11., gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    temperature_ratio = 1. / (1. + 0.5 * (gamma - 1.) * mach_number ** 2)
    u_j = mach_number * np.sqrt(temperature_ratio)
    for i, xyz in enumerate(g.xyz):
        r = np.sqrt(xyz[:,:,0,0] ** 2 + xyz[:,:,0,1] ** 2) / 0.5 + EPSILON
        z = xyz[0,0,:,2]
        theta = theta_j + S * np.where(z > 0., z, np.zeros_like(z))
        for k in range(z.size):
            s.q[i][:,:,k,3] = 0.5 * u_j * (1. + np.tanh(0.25 / theta[k] *
                                                        (1. / r - r)))
            if z[k] > potential_core_length:
                s.q[i][:,:,k,3] *= (1. - np.exp(
                    1.35 / (1. - z[k] / potential_core_length)))
        s.q[i][:,:,:,0] = 1. / (0.5 * (gamma - 1.) * s.q[i][:,:,:,3] / u_j *
                                (1. - s.q[i][:,:,:,3] / u_j) *
                                mach_number ** 2 + s.q[i][:,:,:,3] / u_j +
                                (1. - s.q[i][:,:,:,3] / u_j) /
                                temperature_ratio) / temperature_ratio
    return s.fromprimitive(gamma)

def eigenmodes(mach_number=1.3, theta_j=0.02, show_progress=True):
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
    return modes

def extract_inflow(g, k=42):
    n = g.get_size()
    n[:,2] = k
    gi = p3d.Grid().set_size(n, True)
    for i, xyz in enumerate(g.xyz):
        gi.xyz[i] = xyz[:,:,:k,:]
    return gi

def inflow_perturbations(g, modes):
    phi_p = 2. * np.pi * random.rand(len(modes))
    phi_n = 2. * np.pi * random.rand(len(modes))
    sr = p3d.Solution().copy_from(g)
    si = p3d.Solution().copy_from(g)
    sr._format.aux_header[1] = si._format.aux_header[1] = mode[0].omega
    for i in range(g.nblocks):
        r = np.sqrt(g.xyz[i][:,:,0,0] ** 2 + g.xyz[i][:,:,0,1] ** 2)
        theta = np.arctan2(g.xyz[i][:,:,0,1], g.xyz[i][:,:,0,0])
        z = g.xyz[i][0,0,:,2]
        si.q[i].fill(0.)
        sr.q[i].fill(0.)
        for j, mode in enumerate(modes):
            q = mode.get_eigenmode(r)
            n = mode.n
            print q.shape
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
                fi[j,k,:,0] *= p3d.cubic_bspline_support(z, z_min, z_max)
            for j in range(n[i][2]):
                fi[:,k,j,0] *= p3d.tanh_support(r[:,k], r_min, r_max, 10., 0.2)
        kmin, kmax = p3d.find_extents(z, z_min, z_max)
        imin, imax = p3d.find_extents(np.mean(r, axis=1), r_min, r_max)
        print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
            'controlRegion.' + block_code[i], 'ACTUATOR',
            i + 1, 0, imin, imax, 1, -1, kmin, kmax)
    return f

if __name__ == '__main__':
    grid(num_axial=481, num_radial=281, num_azimuthal=192,
         p_inner=1.12148531779, interface_ds_ratio=0.999793250757)
    g = p3d.fromfile('OSUMach1.3.xyz')
    target_state(g).save('OSUMach1.3.target.q')
    target_mollifier(g).save('OSUMach1.3.target_mollifier.f')
    control_mollifier(g).save('OSUMach1.3.control_mollifier.f')
    gi = extract_inflow(g)
    gi.save('OSUMach1.3.inflow.xyz')
    modes = eigenmodes()
    for i, mode in enumerate(modes):
        sr, si = inflow_perturbations(gi, mode)
        sr.save('OSUMach1.3-%02d.eigenmode_real.q' % (i + 1))
        si.save('OSUMach1.3-%02d.eigenmode_imag.q' % (i + 1))
