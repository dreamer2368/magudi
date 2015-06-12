import plot3dnasa as p3d
import numpy as np

def reynolds_average(s, chunk_size=50):
    n = s.get_size(0)
    sm = p3d.Solution().set_size([n[1], 1, 1], True)
    sm.q[0].fill(0.)
    ends = [int(d) for d in np.linspace(chunk_size - 1, n[2] - 1,
                                        n[2] / chunk_size)]
    subzones = [([0, 0, ends[i-1] + 1], [-1, -1, ends[i]])
                if i > 0 else ([0, 0, 0], [-1, -1, ends[0]])
                for i in range(len(ends))]
    for subzone in subzones:
        sm.q[0][:,0,0,:] += np.sum(np.sum(s.set_subzone(
            0, *subzone).load().toprimitive().q[0], axis=2), axis=0)
    sm.q[0] /= (n[0] * n[2])
    return sm

def wall_units(y, q, Re=2266., gamma=1.4, van_driest=True):
    T = gamma * q[:,4] / ((gamma - 1.) * q[:,0])
    mu_wall = ((gamma - 1.) * T[0]) ** 0.666 / Re
    tau_wall = mu_wall * (q[1,1] - q[0,1]) / (y[1] - y[0])
    u_tau = np.sqrt(tau_wall / q[0,0])
    print np.append([0.], np.cumsum(np.sqrt(T[0] / T[:-1]) * (u[1:] - u[:-1]))).shape
    if van_driest:
        return y * q[0,0] * u_tau / mu_wall, \
            np.append([0.], np.cumsum(np.sqrt(T[0] / T[:-1]) *
                                      (u[1:] - u[:-1]))).shape / u_tau
    return y * q[0,0] * u_tau / mu_wall, q[:,1] / u_tau
