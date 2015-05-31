import numpy as np
import PLOT3D

def getSliceCoordinates(gridFile, sliceAxis = -1, hasIBLANK = False):
    """Slices a three-dimensional single-block Cartesian grid along the
    specified axis and returns the remaining 2 coordinates as a tuple of
    2D arrays.
    """
    grid = PLOT3D.Grid()
    grid.hasIBLANK = hasIBLANK
    grid.ImportSubarray(gridFile, 0, [0, 0, 0], [0 if axis == sliceAxis else -1 for axis in range(3)])
    slice_ = [slice(0, 0) if axis == sliceAxis else slice(None) for axis in range(3)]
    axes = [0, 1, 2]
    del axes[sliceAxis]
    return np.squeeze(grid.X[0][slice_ + [axes[0]]]), np.squeeze(grid.X[0][slice_ + [axes[1]]])

def sliceAdjointContours(ax, gridFile, solutionFile, sliceAxis = -1, sliceIndex = 0, componentIndex = -1, *args, **kwargs):
    """Contours a slice of the specified adjoint component using lines (filled =
    False) or flood (filled = True) style and returns the QuadContoureSet
    object.
    """
    x, y = getSliceCoordinates(gridFile, sliceAxis, kwargs.pop('hasIBLANK', False))
    soln = PLOT3D.Solution()
    soln.ImportSubarray(solutionFile, 0, [sliceIndex if axis == sliceAxis else 0 for axis in range(3)], [sliceIndex if axis == sliceAxis else -1 for axis in range(3)])
    z = np.squeeze(soln.Q[0][[slice(0, 0) if axis == sliceAxis else slice(None) for axis in range(3)] + [componentIndex]])
    plot = ax.contourf if kwargs.pop('filled', False) is True else ax.contour
    return plot(x, y, z, *args, **kwargs)

def sliceVorticityDilatation(ax, gridFile, functionFile, sliceAxis = -1, sliceIndex = 0, **kwargs):
    """Contours vorticity and images dilatation (using a continuous flood style)
    along a slice and returns a tuple containing the QuadContourSet and
    AxesImage objects respectively.
    """
    import matplotlib.mlab as mlab
    x, y = getSliceCoordinates(gridFile, sliceAxis, kwargs.pop('hasIBLANK', False))
    func = PLOT3D.Function()
    func.ImportSubarray(functionFile, 0, [sliceIndex if axis == sliceAxis else 0 for axis in range(3)], [sliceIndex if axis == sliceAxis else -1 for axis in range(3)])
    z = np.squeeze(func.F[0][[slice(0, 0) if axis == sliceAxis else slice(None) for axis in range(3)] + [1]])
    c = ax.contour(x, y, z, **kwargs.pop('vorticity_kwargs', dict()))
    extent = [x[:,0].min(), x[:,0].max(), y[0,:].min(), y[0,:].max()]
    xi, eta = np.mgrid[extent[0]:extent[1]:512j, extent[2]:extent[3]:512j]
    im = ax.imshow(np.flipud(np.transpose(mlab.griddata(x.flatten(), y.flatten(), np.ravel(np.squeeze(func.F[0][[slice(0, 0) if axis == sliceAxis else slice(None) for axis in range(3)] + [0]])), xi, eta, interp = 'linear'))), extent = extent, **kwargs.pop('dilatation_kwargs', dict()))
    return c, im
