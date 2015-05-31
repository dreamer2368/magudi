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

def getCenterlineRMSFluctuations(gridFile, meanSolutionFile, solutionFileList, **kwargs):
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
        from progressbar import ProgressBar, SimpleProgress
        p = ProgressBar(widgets = ['Updating time-average: ', SimpleProgress()], maxval = len(solutionFileList)).start()
    except:
        pass
    for i, solutionFile in solutionFileList:
        soln.ImportSubarray(solutionFile, 0, [n[0]/2, n[1]/2, 0], [n[0]/2, n[1]/2, -1])
        soln.TransformToPrimitiveVariables(ratioOfSpecificHeats)
        Q += (soln.Q[0][0,0,:,:] - meanSoln.Q[0][0,0,:,:]) ** 2
        p.update(i) if p is not None
    p.finish() if p is not None
    Q = np.sqrt(Q / len(solutionList))
