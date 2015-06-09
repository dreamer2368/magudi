def rc_journal(width=None, height=None):
    from numpy import sqrt
    import matplotlib
    golden_ratio = 0.5 * (1. + sqrt(5.))
    matplotlib.rcParams['backend'] = 'pgf'
    matplotlib.rcParams['font.size'] = 8.
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['font.family'] = 'serif'
    w = width if width is not None else 4.5
    h = height if height is not None else width / golden_ratio
    matplotlib.rcParams['figure.figsize'] = (w, h)    
    matplotlib.rcParams['lines.linewidth'] = 0.5
    matplotlib.rcParams['patch.linewidth'] = 0.5
    matplotlib.rcParams['axes.linewidth'] = 0.5
    matplotlib.rcParams['grid.linewidth'] = 0.25
    matplotlib.rcParams['grid.linestyle'] = '-'
    matplotlib.rcParams['grid.color'] = '#AAAAAA'
    matplotlib.rcParams['axes.axisbelow'] = False
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    matplotlib.rcParams['legend.fontsize'] = 'medium'

def style_tecplot(ax):
    ax.minorticks_on()
    ax.tick_params(which='major', width=0.75, length=3)
    ax.tick_params(which='minor', width=0.5, length=1.5)
    return ax

def style_spine(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(which='both', axis='both', direction='out')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    return ax

def nice_labels(t, scale='linear'):
    from numpy import log10
    if scale == 'log':
        t, [r'$10^{%g}$' % t_ if t_ != 0 else r'$1$' for t_ in np.log10(t)]
    elif scale == 'linear':
        return t, [r'$%g$' % t_ for t_ in t]
    return t, [r'' for t_ in t]

def read_engauge_data(filename):
    """Reads CSV output of Engauge plot digitzer and returns a list of numpy
    ndarrays of shape (2, n) containing x-y data.
    """
    import numpy as np
    with open(filename, 'r') as f:
        raw_data = f.readlines()
    offsets = []
    for i, line in enumerate(raw_data):
        if line[0] == '#':
            offsets = offsets + [i + 1]
    curves = [ None ] * len(offsets)
    offsets = offsets + [len(raw_data)+1]
    for i in range(len(curves)):
        curves[i] = np.array([float(s) for s in ''.join(
            raw_data[offsets[i]:offsets[i+1]-1]).replace(',', ' ').split()])
        curves[i] = np.reshape(curves[i], [2, curves[i].size / 2], order='F')
    return curves    

def setup_axis(ax, xlim, ylim, xlabel, ylabel, xticks, yticks,
               xlabelpad=2., ylabelpad=2., xscale='linear', yscale='linear'):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel, labelpad=xlabelpad)
    ax.set_ylabel(ylabel, labelpad=ylabelpad)
    ticks, ticklabels = nice_labels(xticks, xscale)
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)
    ticks, ticklabels = nice_labels(yticks, yscale)
    ax.set_yticks(ticks)
    ax.set_yticklabels(ticklabels)
    return ax
