def TargetPresentation(width = None, height = None):
    import numpy as np
    import matplotlib
    GoldenRatio = 0.5 * (1 + np.sqrt(5))
    matplotlib.rcParams['backend'] = 'Agg'
    matplotlib.rcParams['font.size'] = 10.0
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}\usepackage{amsfonts}\usepackage{cmbright}"
    if width is not None:
        if height is not None:
            matplotlib.rcParams['figure.figsize'] = (width, height)
        else:
            matplotlib.rcParams['figure.figsize'] = (width, width / GoldenRatio)
    else:
        if height is not None:
            matplotlib.rcParams['figure.figsize'] = (height * GoldenRatio, height)
        else:
            matplotlib.rcParams['figure.figsize'] = (4.5, 4.5 / GoldenRatio)
    matplotlib.rcParams['lines.linewidth'] = 0.5
    matplotlib.rcParams['patch.linewidth'] = 0.5
    matplotlib.rcParams['axes.linewidth'] = 0.5
    matplotlib.rcParams['grid.linewidth'] = 0.25
    matplotlib.rcParams['grid.linestyle'] = '-'
    matplotlib.rcParams['grid.color'] = '#AAAAAA'
    matplotlib.rcParams['axes.axisbelow'] = False
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    matplotlib.rcParams['legend.fontsize'] = 'medium'

def TargetAnimation(width = None, height = None):
    import matplotlib
    TargetPresentation(width, height)
    matplotlib.rcParams['animation.bitrate'] = 18000

def TargetJournal(width = None, height = None):
    import matplotlib
    TargetPresentation(width, height)
    matplotlib.rcParams['backend'] = 'pgf'
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 8.0
    matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}\usepackage{amssymb}"

def ApplyTecplotStyle(ax):
    ax.minorticks_on()
    ax.tick_params(which = 'major', width = 0.75, length = 3)
    ax.tick_params(which = 'minor', width = 0.5, length = 1.5)
    return ax

def ApplySpineStyle(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(which = 'both', axis = 'both', direction = 'out')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    return ax

def NiceLabels(t):
    return t, [(r"$%g$" % t_) for t_ in t]

def EmptyLabels(t):
    return t, [("") for t_ in t]

def NiceLogLabels(t):
    import numpy as np
    return t, [(r"$10^{%g}$" % t_ if t_ != 0 else r"$1$") for t_ in np.log10(t)]

def ReadDigitizedLines(csvFile):
    """Reads CSV output of Engauge plot digitzer and returns a list of numpy
    ndarrays of shape (2, n) containing x-y data.
    """
    import numpy as np
    rawData = open(csvFile, 'r').readlines()
    curveOffsets = []
    for i, line in enumerate(rawData):
        if line[0] == '#':
            curveOffsets = curveOffsets + [i + 1]
    curves = [ None ] * len(curveOffsets)
    curveOffsets = curveOffsets + [len(rawData)+1]
    for i in range(len(curves)):
        curves[i] = np.array([float(s) for s in ''.join(rawData[curveOffsets[i]:curveOffsets[i+1]-1]).replace(',', ' ').split()])
        curves[i] = np.reshape(curves[i], [2, curves[i].size / 2], order = 'F')
    return curves
