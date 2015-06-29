#!/usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt

import matplotlibhelper as mh

if __name__ == '__main__':
    mh.rc_journal(width=9.0, height=5.0)
    ax = plt.subplot(111)
    a = np.loadtxt('HomogeneousMixture.gradient_error.txt', usecols=[1, 4],
                   skiprows=1)
    ax.loglog(np.abs(a[:,0]), np.abs(a[:,1]), 'rs-', mec='r', mfc='w', ms=3.5)
    mh.setup_axis(ax, [1e-7, 1e2], [1e-7, 1e1],
                  r'$\alpha$', r'$\epsilon$',
                  10 ** np.linspace(-7., 2., 10),
                  10 ** np.linspace(-7., 1., 10),
                  xscale='log', yscale='log', xlabelpad=4, ylabelpad=0)
    mh.style_tecplot(ax)
    plt.tight_layout(pad=0.5)
    filename = os.path.splitext(os.path.basename(__file__))[0] + '.pdf'
    plt.savefig(filename, dpi=300)
