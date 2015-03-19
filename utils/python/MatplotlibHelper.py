#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy import sqrt, log10
GoldenRatio = 0.5 * (1 + sqrt(5))

def TargetPresentation(width = None, height = None):
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
    TargetPresentation(width, height)
    matplotlib.rcParams['animation.bitrate'] = 18000

def TargetJournal(width = None, height = None):
    TargetPresentation(width, height)
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}\usepackage{amsfonts}"

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
    return t, [(r"$10^{%g}$" % t_) for t_ in log10(t)]

from matplotlib.collections import Collection
from matplotlib.artist import allow_rasterization

class ListCollection(Collection):
     def __init__(self, collections, **kwargs):
         Collection.__init__(self, **kwargs)
         self.set_collections(collections)
     def set_collections(self, collections):
         self._collections = collections
     def get_collections(self):
         return self._collections
     @allow_rasterization
     def draw(self, renderer):
         for _c in self._collections:
             _c.draw(renderer)

def RasterizeContour(c):
     collections = c.collections
     for _c in collections:
         _c.remove()
     cc = ListCollection(collections, rasterized = True)
     ax = plt.gca()
     ax.add_artist(cc)
     return cc
