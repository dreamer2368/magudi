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

def wall_scales(y, q, Re=2266., gamma=1.4):
    T = gamma * q[:,4] / ((gamma - 1.) * q[:,0])
    mu_wall = ((gamma - 1.) * T[0]) ** 0.666 / Re
    tau_wall = mu_wall * (q[1,1] - q[0,1]) / (y[1] - y[0])
    u_tau = np.sqrt(tau_wall / q[0,0])
    return mu_wall / (q[0,0] * u_tau), u_tau

def wall_units(y, q, Re=2266., gamma=1.4, van_driest=True):
    y_tau, u_tau = wall_scales(y, q, Re, gamma)
    if not van_driest:
        return y / y_tau, q[:,1] / u_tau
    T = gamma * q[:,4] / ((gamma - 1.) * q[:,0])
    return y / y_tau, np.append([0.], np.cumsum(
        np.sqrt(T[0] / T[:-1]) * (q[1:,1] - q[:-1,1]))) / u_tau

def reynolds_stresses(s, qm, chunk_size=50):
    n = s.get_size(0)
    f = p3d.Function(ncomponents=6).set_size([n[1], 1, 1], True)
    f.f[0].fill(0.)
    ends = [int(d) for d in np.linspace(chunk_size - 1, n[2] - 1,
                                        n[2] / chunk_size)]
    subzones = [([0, 0, ends[i-1] + 1], [-1, -1, ends[i]])
                if i > 0 else ([0, 0, 0], [-1, -1, ends[0]])
                for i in range(len(ends))]
    for subzone in subzones:
        s.set_subzone(0, *subzone).load().toprimitive()
        l = 0
        for j in range(3):
            for k in range(j+1):
                for i in range(n[1]):
                    f.f[0][i,0,0,l] += np.sum((s.q[0][:,i,:,j+1] - qm[i,j+1]) *
                                              (s.q[0][:,i,:,k+1] - qm[i,k+1]))
                l += 1
    for i in range(f.ncomponents):
        f.f[0][:,0,0,i] *= qm[:,0] / (n[0] * n[2])
    return f
