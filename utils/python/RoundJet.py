import numpy as np
import PLOT3D

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

def getAzimuthalAveragedSolution(gridFile, solutionFile, radialCoordinate = 0.5, stencilSize = 5, hasIBLANK = False, showProgress = True):
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
