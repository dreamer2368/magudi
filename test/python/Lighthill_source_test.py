#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def grid(size):
    from numpy.random import rand
    L = 10. * rand() + 5.

    x_min = 0.
    x_max = L
    y_min = 0.
    y_max = L
    g = p3d.Grid().set_size(size, True)
    x = np.linspace(x_min, x_max, g.size[0,0])
    y = np.linspace(y_min, y_max, g.size[0,1])
    g.xyz[0][:,:,0,:2] = np.transpose(np.meshgrid(x, y))
    return g

def inviscid_manufactured_solution(g, gamma=1.4):
    from numpy.random import rand

    n = g.get_size(0)[0]

    xg, yg = g.xyz[0][:,:,0,0], g.xyz[0][:,:,0,1]
    L = xg[-1,0] - xg[0,0]

    rho_coeffs = rand(1)
    rho = rho_coeffs[0]

    u_coeffs = rand(2)
    u = u_coeffs[0] * np.sin(2.*np.pi*xg/(u_coeffs[1]*L))

    v_coeffs = rand(2)
    v = v_coeffs[0] * np.cos(2.*np.pi*yg/(v_coeffs[1]*L))

    p_coeffs = rand(2)
    p = p_coeffs[0] * np.sin(2.*np.pi*yg/(p_coeffs[1]*L))

    s = p3d.Solution().copy_from(g).quiescent(gamma)
    for i, xyz in enumerate(g.xyz):
        # s.q[i][:,:,0,0] = rho0
        s.q[i][:,:,0,0] = rho
        s.q[i][:,:,0,1] = u
        s.q[i][:,:,0,2] = v
        s.q[i][:,:,0,4] = p

    f = p3d.Function().copy_from(g)
    f.f[0][:,:,0,0] = 8.*(np.pi**2)*(u_coeffs[0]**2)*rho_coeffs[0]/L/L/(u_coeffs[1]**2)*(np.cos(2.*np.pi*xg/L/u_coeffs[1]))**2          \
                        - 8.*(np.pi**2)*(v_coeffs[0]**2)*rho_coeffs[0]/L/L/(v_coeffs[1]**2)*(np.cos(2.*np.pi*yg/L/v_coeffs[1]))**2      \
                        - 8.*(np.pi**2)*(u_coeffs[0]**2)*rho_coeffs[0]/L/L/(u_coeffs[1]**2)*(np.sin(2.*np.pi*xg/L/u_coeffs[1]))**2      \
                        - 4.*(np.pi**2)*p_coeffs[0]/L/L/(p_coeffs[1]**2)*np.sin(2.*np.pi*yg/L/p_coeffs[1])                              \
                        + 8.*(np.pi**2)*(v_coeffs[0]**2)*rho_coeffs[0]/L/L/(v_coeffs[1]**2)*(np.sin(2.*np.pi*yg/L/v_coeffs[1]))**2      \
                        - 8.*(np.pi**2)*u_coeffs[0]*v_coeffs[0]*rho_coeffs[0]/L/L/u_coeffs[1]/v_coeffs[1]\
                                *np.cos(2.*np.pi*xg/L/u_coeffs[1])*np.sin(2.*np.pi*yg/L/v_coeffs[1])
    return s.fromprimitive(gamma), f

def viscous_manufactured_solution(g, gamma=1.4, Re=100.):
    from numpy.random import rand

    xg, yg = g.xyz[0][:,:,0,0], g.xyz[0][:,:,0,1]
    L = xg[-1,0] - xg[0,0]

    rho_coeffs = rand(1)
    rho = rho_coeffs[0]

    p_coeffs = rand(1)
    p = p_coeffs[0]

    u_coeffs = rand(2)
    u = u_coeffs[0] * np.sin(2.*np.pi*xg/(u_coeffs[1]*L))

    v_coeffs = rand(2)
    v = v_coeffs[0] * np.cos(2.*np.pi*yg/(v_coeffs[1]*L))

    s = p3d.Solution().copy_from(g).quiescent(gamma)
    for i, xyz in enumerate(g.xyz):
        # s.q[i][:,:,0,0] = rho0
        s.q[i][:,:,0,0] = rho
        s.q[i][:,:,0,1] = u
        s.q[i][:,:,0,2] = v
        s.q[i][:,:,0,4] = p

    T = gamma/(gamma-1.) * p/rho
    n = 0.666
    mu = ((gamma-1.)*T)**n
    lam = 0.6*mu - 2./3.*mu
    f_from_tau = - 8.*(np.pi**3)*u_coeffs[0]*lam/(L**3)/(u_coeffs[1]**3)*np.cos(2.*np.pi*xg/L/u_coeffs[1])                              \
                    - 16.*(np.pi**3)*u_coeffs[0]*mu/(L**3)/(u_coeffs[1]**3)*np.cos(2.*np.pi*xg/L/u_coeffs[1])                           \
                    + 8.*(np.pi**3)*v_coeffs[0]*lam/(L**3)/(v_coeffs[1]**3)*np.sin(2.*np.pi*yg/L/v_coeffs[1])                           \
                    + 16.*(np.pi**3)*v_coeffs[0]*mu/(L**3)/(v_coeffs[1]**3)*np.sin(2.*np.pi*yg/L/v_coeffs[1])
    f_from_tau /= -Re

    f = p3d.Function().copy_from(g)
    f.f[0][:,:,0,0] = 8.*(np.pi**2)*(u_coeffs[0]**2)*rho_coeffs[0]/L/L/(u_coeffs[1]**2)*(np.cos(2.*np.pi*xg/L/u_coeffs[1]))**2          \
                        - 8.*(np.pi**2)*(v_coeffs[0]**2)*rho_coeffs[0]/L/L/(v_coeffs[1]**2)*(np.cos(2.*np.pi*yg/L/v_coeffs[1]))**2      \
                        - 8.*(np.pi**2)*(u_coeffs[0]**2)*rho_coeffs[0]/L/L/(u_coeffs[1]**2)*(np.sin(2.*np.pi*xg/L/u_coeffs[1]))**2      \
                        + 8.*(np.pi**2)*(v_coeffs[0]**2)*rho_coeffs[0]/L/L/(v_coeffs[1]**2)*(np.sin(2.*np.pi*yg/L/v_coeffs[1]))**2      \
                        - 8.*(np.pi**2)*u_coeffs[0]*v_coeffs[0]*rho_coeffs[0]/L/L/u_coeffs[1]/v_coeffs[1]\
                                *np.cos(2.*np.pi*xg/L/u_coeffs[1])*np.sin(2.*np.pi*yg/L/v_coeffs[1])
    f.f[0][:,:,0,0] += f_from_tau
    return s.fromprimitive(gamma), f

if __name__ == '__main__':
    g = grid([201, 201])
    g.save('LighthillTest.xyz')
    s, f = inviscid_manufactured_solution(g)
    s.save('LighthillTest.inviscid.q')
    f.save('LighthillTest.inviscid.analytic.f')

    s, f = viscous_manufactured_solution(g,gamma=1.4,Re=1.7)
    s.save('LighthillTest.viscous.q')
    f.save('LighthillTest.viscous.analytic.f')
