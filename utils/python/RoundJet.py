import numpy as np
import PLOT3D
import plot3dnasa as p3d

def mergeMeanSolutions(prefix,meanSolutionFileList, **kwargs):

    """Merge multiple successive mean solution"""

    N = len(meanSolutionFileList)

    ratioOfSpecificHeats = kwargs.pop('ratioOfSpecificHeats', 1.4)

    g = p3d.Grid('%s.xyz' % prefix)
    s = p3d.Solution().copy_from(g).quiescent(ratioOfSpecificHeats)
    for i, xyz in enumerate(g.xyz):
        s.q[i] = np.zeros_like(s.q[i])

    for k in range(0,N):
        meanSoln = p3d.fromfile(meanSolutionFileList[k])
        for i, xyz in enumerate(g.xyz):
            s.q[i] += 1./N*meanSoln.q[i]

    return s.fromprimitive(ratioOfSpecificHeats)
 
def getCenterlineSubarray(gridFile, solutionFile = None, functionFile = None, hasIBLANK = False):

    """Extracts the centerline axial coordinate, solution and function
    components and returns them as a tuple.
    """

    grid = PLOT3D.Grid()
    grid.hasIBLANK = hasIBLANK
    grid.ImportSkeleton(gridFile)
    n = grid.GetSize(0)
    grid.ImportSubarray(gridFile, 0, [n[0]/2, n[1]/2, 0], [n[0]/2, n[1]/2, -1])
    a = [grid.X[0][0,0,:,2]]

    if solutionFile is not None:
        soln = PLOT3D.Solution()
        soln.ImportSubarray(solutionFile, 0, [n[0]/2, n[1]/2, 0], [n[0]/2, n[1]/2, -1])
        a = a + [soln.Q[0][0,0,:,:]]
    if functionFile is not None:
        func = PLOT3D.Function()
        func.ImportSubarray(functionFile, 0, [n[0]/2, n[1]/2, 0], [n[0]/2, n[1]/2, -1])
        a = a + [func.F[0][0,0,:,:]]
    return tuple(a)

def getCenterlineRMSFluctuations(meanSolutionFile, solutionFileList, **kwargs):

    """Computes the RMS fluctuations of primitive variables along the centerline
    and returns them as 1D single-block solution object.
    """

    ratioOfSpecificHeats = kwargs.pop('ratioOfSpecificHeats', 1.4)
    meanSoln = PLOT3D.Solution()
    meanSoln.ImportSkeleton(meanSolutionFile)
    n = meanSoln.GetSize(0)
    meanSoln.ImportSubarray(meanSolutionFile, 0, [n[0]/2, n[1]/2, 0], [n[0]/2, n[1]/2, -1])
    meanSoln.TransformToPrimitiveVariables(ratioOfSpecificHeats)
    Q = np.zeros([meanSoln.GetSize(0)[-1], 5])

    p = None
    try:
        if kwargs.pop('showProgress', True) is True:
            from progressbar import ProgressBar, Percentage, Bar, ETA
            print 'Averaging in time:'
            p = ProgressBar(widgets = [Percentage(), ' ', Bar('=', left = '[', right = ']'), ' ', ETA()], maxval = len(solutionFileList)).start()
    except:
        pass

    for i, solutionFile in enumerate(solutionFileList):
        soln = PLOT3D.Solution()
        soln.ImportSubarray(solutionFile, 0, [n[0]/2, n[1]/2, 0], [n[0]/2, n[1]/2, -1])
        soln.TransformToPrimitiveVariables(ratioOfSpecificHeats)
        Q += (soln.Q[0][0,0,:,:] - meanSoln.Q[0][0,0,:,:]) ** 2
        if p is not None:
            p.update(i)

    if p is not None:
        p.finish()

    Q = np.sqrt(Q / len(solutionFileList))
    return Q

def interpolateSolutionAtConstantRadius(gridFile, solutionFile, radialCoordinate = 0.5, stencilSize = 5, hasIBLANK = False, showProgress = True):

    """Interpolates the grid coordinates and solution components at a
    constant-radius surface using a Barycentric polynomial interpolation scheme.
    """

    from scipy.interpolate import BarycentricInterpolator

    # Find nearest neighbors and setup Barycentric interpolators.
    grid = PLOT3D.Grid()
    grid.hasIBLANK = hasIBLANK
    grid.ImportSubarray(gridFile, 1, [0, 0, 0], [-1, -2, 0])
    interpolators = [ None ] * (grid.GetSize(0)[1])
    iRadialNearest = np.array([np.argmin(np.abs(grid.X[0][:,i,0,0] ** 2 + grid.X[0][:,i,0,1] ** 2 - radialCoordinate ** 2)) for i in range(grid.GetSize(0)[1])], dtype = int)
    iRadialOffset = iRadialNearest.min() - stencilSize / 2
    interpolators = [BarycentricInterpolator(np.sqrt(grid.X[0][iRadialNearest[i]-stencilSize/2:iRadialNearest[i]+stencilSize/2+1,i,0,0] ** 2 + grid.X[0][iRadialNearest[i]-stencilSize/2:iRadialNearest[i]+stencilSize/2+1,i,0,1] ** 2)) for i in range(iRadialNearest.size)]
    iRadialNearest -= iRadialOffset

    axialCoordinate, = getCenterlineSubarray(gridFile)
    X = np.empty([4 * len(interpolators), axialCoordinate.size, 3])
    Q = np.empty([4 * len(interpolators), axialCoordinate.size, 5])
    for k in range(axialCoordinate.size):
        X[:,k,2] = axialCoordinate[k]

    p = None
    try:
        if showProgress is True:
            from progressbar import ProgressBar, Percentage, Bar, ETA
            print 'Interpolating at constant-radius surface:'
            p = ProgressBar(widgets = [Percentage(), ' ', Bar('=', left = '[', right = ']'), ' ', ETA()], maxval = 28).start()
    except:
        pass

    soln = PLOT3D.Solution()
    for i in range(1, 5):
        grid.ImportSubarray(gridFile, i, [iRadialOffset + iRadialNearest.min() - stencilSize / 2, 0, 0], [iRadialOffset + iRadialNearest.max() + stencilSize / 2, -2, 0])
        soln.ImportSubarray(solutionFile, i, [iRadialOffset + iRadialNearest.min() - stencilSize / 2, 0, 0], [iRadialOffset + iRadialNearest.max() + stencilSize / 2, -2, -1])
        for l in range(2): # grid coordinates
            for j, interpolator in enumerate(interpolators):
                interpolator.set_yi(grid.X[0][iRadialNearest[j]-stencilSize/2:iRadialNearest[j]+stencilSize/2+1,j,0,l])
                X[(i-1)*len(interpolators)+j,:,l] = interpolator(radialCoordinate)
            if p is not None:
                p.update((i - 1) * 7 + l + 1)
        for l in range(5): # solution components
            for k in range(axialCoordinate.size):
                for j, interpolator in enumerate(interpolators):
                    interpolator.set_yi(soln.Q[0][iRadialNearest[j]-stencilSize/2:iRadialNearest[j]+stencilSize/2+1,j,k,l])
                    Q[(i-1)*len(interpolators)+j,k,l] = interpolator(radialCoordinate)
            if p is not None:
                p.update((i - 1) * 7 + l + 3)

    if p is not None:
        p.finish()

    return X, Q

def computeAzimuthalAverageFromConstantRadiusSolution(X, Q):
    """Returns the azimuthal-averaged solution from the grid coordinates and
    solution components corresponding to a constant-radius surface as returned,
    for example, using the interpolateSolutionAtConstantRadius function.
    """
    azimuthalAngle = np.arctan2(X[:,0,1], X[:,0,0])
    iRoll = -np.argmin(np.abs(azimuthalAngle))
    azimuthalAngle = np.roll(azimuthalAngle, iRoll)
    azimuthalAngle[azimuthalAngle < 0.] += 2. * np.pi
    return np.tensordot(np.roll(Q, iRoll, axis = 0), np.diff(np.append(azimuthalAngle, azimuthalAngle[0] + 2. * np.pi)), axes = ([0, 0])) / (2. * np.pi)

def getAzimuthalAverageFromConstantRadiusSolution(gridFile, solutionFile, radialCoordinate = 0.5, stencilSize = 5, hasIBLANK = False, showProgress = True):
    """Returns the azimuthal-averaged solution after interpolating it on a
    constant-radius surface.
    """
    X, Q = interpolateSolutionAtConstantRadius(gridFile, solutionFile, radialCoordinate, stencilSize, hasIBLANK, showProgress)
    return computeAzimuthalAverageFromConstantRadiusSolution(X, Q)

def getConstantRadiusRMSFluctuations(gridFile, meanSolutionFile, solutionFileList, **kwargs):

    """Computes the RMS fluctuations of primitive variables from a
    constant-radius solution obtained, for example, using the
    interpolateSolutionAtConstantRadius function.
    """

    grid = PLOT3D.Grid()
    grid.ImportSubarray(gridFile, 0, [0, 0, 0], [0, -1, 0])
    X = grid.X[0][0,:,:,:]

    ratioOfSpecificHeats = kwargs.pop('ratioOfSpecificHeats', 1.4)
    meanSoln = PLOT3D.Solution()
    meanSoln.Import(meanSolutionFile)
    meanSoln.TransformToPrimitiveVariables(ratioOfSpecificHeats)
    n = meanSoln.GetSize(0)
    Q = np.zeros([n[1], n[2], 5])

    p = None
    try:
        if kwargs.pop('showProgress', True) is True:
            from progressbar import ProgressBar, Percentage, Bar, ETA
            print 'Averaging in time:'
            p = ProgressBar(widgets = [Percentage(), ' ', Bar('=', left = '[', right = ']'), ' ', ETA()], maxval = len(solutionFileList)).start()
    except:
        pass

    for i, solutionFile in enumerate(solutionFileList):
        soln = PLOT3D.Solution()
        soln.Import(solutionFile)
        soln.TransformToPrimitiveVariables(ratioOfSpecificHeats)
        Q += (soln.Q[0][0,:,:,:] - meanSoln.Q[0][0,:,:,:]) ** 2
        if p is not None:
            p.update(i)

    if p is not None:
        p.finish()

    return computeAzimuthalAverageFromConstantRadiusSolution(X, np.sqrt(Q / len(solutionFileList)))

def getVerticalLongitudinalSubarray(gridFile, solutionFile = None, functionFile = None, hasIBLANK = False):

    """Extracts the grid coordinates, and solution and function
    components along the x = 0 plane and returns them as a tuple.
    """

    grid = PLOT3D.Grid()
    grid.hasIBLANK = False
    soln = PLOT3D.Solution()
    func = PLOT3D.Function()
    grid.ImportSkeleton(gridFile)
    n = grid.GetSize()
    assert n[0][0] % 2 == 1 and n[2][1] % 2 == 1 and n[4][1] % 2 == 1

    x = np.empty([n[4][0] + n[0][1] + n[2][0] - 2, n[0][2]])
    y = np.empty_like(x)
    if solutionFile is not None:
        Q = np.empty([x.shape[0], x.shape[1], 5])
    if functionFile is not None:
        func.ImportSkeleton(functionFile)
        F = np.empty([x.shape[0], x.shape[1], func.nComponents])

    grid.ImportSubarray(gridFile, 4, [1, n[4][1]/2, 0], [-1, n[4][1]/2, -1])
    x[:n[4][0]-1,:] = grid.X[0][::-1,0,:,1]
    y[:n[4][0]-1,:] = grid.X[0][::-1,0,:,2]
    if solutionFile is not None:
        soln.ImportSubarray(solutionFile, 4, [1, n[4][1]/2, 0], [-1, n[4][1]/2, -1])
        Q[:n[4][0]-1,:,:] = soln.Q[0][::-1,0,:,:]
    if functionFile is not None:
        func.ImportSubarray(functionFile, 4, [1, n[4][1]/2, 0], [-1, n[4][1]/2, -1])
        F[:n[4][0]-1,:,:] = func.F[0][::-1,0,:,:]

    grid.ImportSubarray(gridFile, 0, [n[0][0]/2, 0, 0], [n[0][0]/2, -2, -1])
    x[n[4][0]-1:n[4][0]+n[0][1]-2,:] = grid.X[0][0,:,:,1]
    y[n[4][0]-1:n[4][0]+n[0][1]-2,:] = grid.X[0][0,:,:,2]
    if solutionFile is not None:
        soln.ImportSubarray(solutionFile, 0, [n[0][0]/2, 0, 0], [n[0][0]/2, -2, -1])
        Q[n[4][0]-1:n[4][0]+n[0][1]-2,:,:] = soln.Q[0][0,:,:,:]
    if functionFile is not None:
        func.ImportSubarray(functionFile, 0, [n[0][0]/2, 0, 0], [n[0][0]/2, -2, -1])
        F[n[4][0]-1:n[4][0]+n[0][1]-2,:,:] = func.F[0][0,:,:,:]

    grid.ImportSubarray(gridFile, 2, [0, n[2][1]/2, 0], [-1, n[2][1]/2, -1])
    x[n[4][0]+n[0][1]-2:,:] = grid.X[0][:,0,:,1]
    y[n[4][0]+n[0][1]-2:,:] = grid.X[0][:,0,:,2]
    a = [x, y]
    if solutionFile is not None:
        soln.ImportSubarray(solutionFile, 2, [0, n[2][1]/2, 0], [-1, n[2][1]/2, -1])
        Q[n[4][0]+n[0][1]-2:,:,:] = soln.Q[0][:,0,:,:]
        a = a + [Q]
    if functionFile is not None:
        func.ImportSubarray(functionFile, 2, [0, n[2][1]/2, 0], [-1, n[2][1]/2, -1])
        F[n[4][0]+n[0][1]-2:,:,:] = func.F[0][:,0,:,:]
        a = a + [F]

    return tuple(a)

def getHorizontalLongitudinalSubarray(gridFile, solutionFile = None, functionFile = None, hasIBLANK = False):

    """Extracts the grid coordinates, and solution and function
    components along the y = 0 plane and returns them as a tuple.
    """

    grid = PLOT3D.Grid()
    grid.hasIBLANK = False
    soln = PLOT3D.Solution()
    func = PLOT3D.Function()
    grid.ImportSkeleton(gridFile)
    n = grid.GetSize()
    assert n[0][1] % 2 == 1 and n[1][1] % 2 == 1 and n[3][1] % 2 == 1

    x = np.empty([n[3][0] + n[0][0] + n[1][0] - 2, n[0][2]])
    y = np.empty_like(x)
    if solutionFile is not None:
        Q = np.empty([x.shape[0], x.shape[1], 5])
    if functionFile is not None:
        func.ImportSkeleton(functionFile)
        F = np.empty([x.shape[0], x.shape[1], func.nComponents])

    grid.ImportSubarray(gridFile, 3, [1, n[3][1]/2, 0], [-1, n[3][1]/2, -1])
    x[:n[3][0]-1,:] = grid.X[0][::-1,0,:,0]
    y[:n[3][0]-1,:] = grid.X[0][::-1,0,:,2]
    if solutionFile is not None:
        soln.ImportSubarray(solutionFile, 3, [1, n[3][1]/2, 0], [-1, n[3][1]/2, -1])
        Q[:n[3][0]-1,:,:] = soln.Q[0][::-1,0,:,:]
    if functionFile is not None:
        func.ImportSubarray(functionFile, 3, [1, n[3][1]/2, 0], [-1, n[3][1]/2, -1])
        F[:n[3][0]-1,:,:] = func.F[0][::-1,0,:,:]

    grid.ImportSubarray(gridFile, 0, [0, n[0][1]/2, 0], [-2, n[0][1]/2, -1])
    x[n[3][0]-1:n[3][0]+n[0][0]-2,:] = grid.X[0][0,:,:,0]
    y[n[3][0]-1:n[3][0]+n[0][0]-2,:] = grid.X[0][0,:,:,2]
    if solutionFile is not None:
        soln.ImportSubarray(solutionFile, 0, [0, n[0][1]/2, 0], [-2, n[0][1]/2, -1])
        Q[n[3][0]-1:n[3][0]+n[0][0]-2,:,:] = soln.Q[0][0,:,:,:]
    if functionFile is not None:
        func.ImportSubarray(functionFile, 0, [0, n[0][1]/2, 0], [-2, n[0][1]/2, -1])
        F[n[3][0]-1:n[3][0]+n[0][0]-2,:,:] = func.F[0][0,:,:,:]

    grid.ImportSubarray(gridFile, 1, [0, n[1][1]/2, 0], [-1, n[1][1]/2, -1])
    x[n[3][0]+n[0][0]-2:,:] = grid.X[0][:,0,:,0]
    y[n[3][0]+n[0][0]-2:,:] = grid.X[0][:,0,:,2]
    a = [x, y]
    if solutionFile is not None:
        soln.ImportSubarray(solutionFile, 1, [0, n[1][1]/2, 0], [-1, n[1][1]/2, -1])
        Q[n[3][0]+n[0][0]-2:,:,:] = soln.Q[0][:,0,:,:]
        a = a + [Q]
    if functionFile is not None:
        func.ImportSubarray(functionFile, 1, [0, n[1][1]/2, 0], [-1, n[1][1]/2, -1])
        F[n[3][0]+n[0][0]-2:,:,:] = func.F[0][:,0,:,:]
        a = a + [F]

    return tuple(a)

def getAzimuthalAveragedGridAndSolution(gridFile, solutionFile, hasIBLANK = False, showProgress = True):

    """Returns a tuple representing the streamwise and radial coordinates, and
    the solution on an axisymmetric grid by averaging along the azimuthal
    direction using a histogram technique.
    """

    r, z = getVerticalLongitudinalSubarray(gridFile)
    z = z[r[:,0] >= 0.,:]
    r = r[r[:,0] >= 0.,:] ** 2
    Q = np.empty([r.shape[0], r.shape[1], 5])
    bins = [r[0,0] - 0.5 * (r[1,0] - r[0,0])] + list(0.5 * (r[:-1,0] + r[1:,0])) + [r[-1,0] + 0.5 * (r[-1,0] - r[-2,0])]

    grid = PLOT3D.Grid()
    grid.hasIBLANK = hasIBLANK
    grid.ImportSkeleton(gridFile)
    soln = PLOT3D.Solution()
    x = np.empty(np.sum([np.prod(grid.GetSize(i)[0:2]) for i in range(5)]))
    y = np.empty([x.size, 5])

    j = 0
    for i in range(5):
        grid.ImportSubarray(gridFile, i, [0, 0, 0], [-1, -1, 0])
        n = np.prod(grid.GetSize(0))
        x[j:j+n] = grid.X[0][:,:,0,0].flatten() ** 2 + grid.X[0][:,:,0,1].flatten() ** 2
        j += n

    p = None
    try:
        if showProgress is True:
            from progressbar import ProgressBar, Percentage, Bar, ETA
            print 'Processing along streamwise direction:'
            p = ProgressBar(widgets = [Percentage(), ' ', Bar('=', left = '[', right = ']'), ' ', ETA()], maxval = Q.shape[1]).start()
    except:
        pass

    for k in range(Q.shape[1]):
        j = 0
        for i in range(5):
            soln.ImportSubarray(solutionFile, i, [0, 0, k], [-1, -1, k])
            n = np.prod(soln.GetSize(0))
            for l in range(5):
                y[j:j+n,l] = soln.Q[0][:,:,0,l].flatten()
            j += n
        for i in range(5):
            Q[:,k,i] = np.histogram(x, bins, weights = y[:,i])[0] / np.histogram(x, bins)[0]
        if p is not None:
            p.update(k)

    if p is not None:
        p.finish()

    return z, np.sqrt(r), Q

def interpolateOnAnnularGrid(gridFile, solutionFile, hasIBLANK = False, showProgress = True):

    import matplotlib.mlab as mlab
    import warnings

    rUniform = []
    grid = PLOT3D.Grid()
    grid.hasIBLANK = hasIBLANK
    grid.ImportSkeleton(gridFile)
    n = grid.GetSize()
    grid.ImportSubarray(gridFile, 0, [0, 0, 0], [0, 0, -1])
    axialCoordinate = grid.X[0][0,0,:,2]
    grid.ImportSubarray(gridFile, 0, [n[0][0]/2+2, n[0][1]/2, 0], [-2, n[0][1]/2, 0])
    rUniform += grid.X[0][:,0,0,0].tolist()
    grid.ImportSubarray(gridFile, 1, [0, n[1][1]/2, 0], [-2, n[1][1]/2, 0])
    rUniform = np.array(rUniform + grid.X[0][:,0,0,0].tolist())
    xi, eta = np.meshgrid(rUniform, np.linspace(0., 2. * np.pi, 1 + np.sum([n[i][1] - 1 for i in range(5)])))

    x, y = [], []
    for i in range(5):
        if i == 0:
            grid.ImportSubarray(gridFile, i, [0, 0, 0], [-1, -1, 0])
        else:
            grid.ImportSubarray(gridFile, i, [1, 1, 0], [-1, -1, 0])
        x += np.sqrt(grid.X[0][:,:,0,0] ** 2 + grid.X[0][:,:,0,1] ** 2).flatten().tolist()
        y += np.arctan2(grid.X[0][:,:,0,1], grid.X[0][:,:,0,0]).flatten().tolist()
    x, y = np.array(x), np.array(y)
    y[y < 0.] += 2. * np.pi

    overlapTolerance = 4. * (2. * np.pi / xi.shape[1])
    x = np.append(np.append(x, x[y < overlapTolerance]), x[y > 2. * np.pi - overlapTolerance])
    yOverlap = np.append(np.append(y, y[y < overlapTolerance] + 2. * np.pi), y[y > 2. * np.pi - overlapTolerance] - 2. * np.pi)
    z = np.empty([x.size, 5])

    p = None
    try:
        if showProgress is True:
            from progressbar import ProgressBar, Percentage, Bar, ETA
            print 'Processing along streamwise direction:'
            p = ProgressBar(widgets = [Percentage(), ' ', Bar('=', left = '[', right = ']'), ' ', ETA()], maxval = 5 * axialCoordinate.size).start()
    except:
        pass


    grid = PLOT3D.Grid()
    grid.nGrids = 1
    grid.SetSize(0, [xi.shape[0], xi.shape[1], axialCoordinate.size])
    for k in range(axialCoordinate.size):
        grid.X[0][:,:,k,2] = axialCoordinate[k]
        grid.X[0][:,:,k,0] = xi * np.cos(eta)
        grid.X[0][:,:,k,1] = xi * np.sin(eta)

    soln = PLOT3D.Solution()
    soln.CopyFrom(grid)

    for k in range(axialCoordinate.size):

        j = 0
        for i in range(5):
            soln_ = PLOT3D.Solution()
            if i == 0:
                soln_.ImportSubarray(solutionFile, i, [0, 0, k], [-1, -1, k])
            else:
                soln_.ImportSubarray(solutionFile, i, [1, 1, k], [-1, -1, k])
            n = np.prod(soln_.GetSize(0))
            for l in range(5):
                z[j:j+n,l] = soln_.Q[0][:,:,0,l].flatten()
            j += n
        for i in range(5):
            z[j:,i] = np.append(z[y < overlapTolerance, i], z[y > 2. * np.pi - overlapTolerance, i])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i in range(5):
                soln.Q[0][:,:,k,i] = mlab.griddata(x, yOverlap, z[:,i], xi, eta, interp = 'linear')
                if p is not None:
                    p.update(k * 5 + i)

    if p is not None:
        p.finish()

    soln.Q[0][-1,:,:,:] = soln.Q[0][0,:,:,:]
    return grid, soln


if __name__ == '__main__':
    prefix='MultiblockJet'
    meanSolutionFileList=[]
    for i in range(0,6):
        istr = str(457000 + i*100000)
        meanSolutionFileList += ['MultiblockJet-'+istr+'.100000.mean.q']
    mergeMeanSolutions(prefix,meanSolutionFileList).save('MultiblockJet.mean.q')

