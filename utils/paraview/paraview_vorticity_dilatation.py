# trace generated using paraview version 5.6.0-RC1

#### import the simple module from the paraview
from paraview.simple import *
from os import *
import numpy as np
chdir('/Users/Kevin/Desktop/VortexDynamics')

prefix = 'VortexDynamics'

N=421
FunctionFileNameList=[]
FigureFileNameList=[]
for i in range(0,N):
    istr = i*200
    FunctionFileNameList = FunctionFileNameList + ['%s-%08d.vorticity_dilatation.f'%(prefix,istr)]
    FigureFileNameList = FigureFileNameList + ['%s-%08d.vorticity_dilatation.png'%(prefix,istr)]

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PLOT3D Reader'
vortexDynamicsxyz = PLOT3DReader(FileName='/Users/Kevin/Desktop/VortexDynamics/VortexDynamics.xyz',
    QFileName=['/Users/Kevin/Desktop/VortexDynamics/VortexDynamics.ic.q'],
    FunctionFileName='/Users/Kevin/Desktop/VortexDynamics/VortexDynamics-00000000.vorticity_dilatation.f')
vortexDynamicsxyz.Functions = [110]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1555, 818]

# show data in view
vortexDynamicsxyzDisplay = Show(vortexDynamicsxyz, renderView1)

# get color transfer function/color map for 'Density'
densityLUT = GetColorTransferFunction('Density')
densityLUT.RGBPoints = [1.0, 0.231373, 0.298039, 0.752941, 1.0001220703125, 0.865003, 0.865003, 0.865003, 1.000244140625, 0.705882, 0.0156863, 0.14902]
densityLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Density'
densityPWF = GetOpacityTransferFunction('Density')
densityPWF.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]
densityPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
vortexDynamicsxyzDisplay.Representation = 'Surface'
vortexDynamicsxyzDisplay.ColorArrayName = ['POINTS', 'Density']
vortexDynamicsxyzDisplay.LookupTable = densityLUT
vortexDynamicsxyzDisplay.OSPRayScaleArray = 'Density'
vortexDynamicsxyzDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
vortexDynamicsxyzDisplay.SelectOrientationVectors = 'Momentum'
vortexDynamicsxyzDisplay.ScaleFactor = 153.33333333333334
vortexDynamicsxyzDisplay.SelectScaleArray = 'Density'
vortexDynamicsxyzDisplay.GlyphType = 'Arrow'
vortexDynamicsxyzDisplay.GlyphTableIndexArray = 'Density'
vortexDynamicsxyzDisplay.GaussianRadius = 7.666666666666668
vortexDynamicsxyzDisplay.SetScaleArray = ['POINTS', 'Density']
vortexDynamicsxyzDisplay.ScaleTransferFunction = 'PiecewiseFunction'
vortexDynamicsxyzDisplay.OpacityArray = ['POINTS', 'Density']
vortexDynamicsxyzDisplay.OpacityTransferFunction = 'PiecewiseFunction'
vortexDynamicsxyzDisplay.DataAxesGrid = 'GridAxesRepresentation'
vortexDynamicsxyzDisplay.SelectionCellLabelFontFile = ''
vortexDynamicsxyzDisplay.SelectionPointLabelFontFile = ''
vortexDynamicsxyzDisplay.PolarAxes = 'PolarAxesRepresentation'
vortexDynamicsxyzDisplay.ScalarOpacityFunction = densityPWF
vortexDynamicsxyzDisplay.ScalarOpacityUnitDistance = 38.18174112674134

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
vortexDynamicsxyzDisplay.OSPRayScaleFunction.Points = [-9.459567940017918e-06, 0.0, 0.5, 0.0, 1.53066275150749e-08, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
vortexDynamicsxyzDisplay.ScaleTransferFunction.Points = [-9.459567940017918e-06, 0.0, 0.5, 0.0, 1.53066275150749e-08, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
vortexDynamicsxyzDisplay.OpacityTransferFunction.Points = [-9.459567940017918e-06, 0.0, 0.5, 0.0, 1.53066275150749e-08, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
vortexDynamicsxyzDisplay.DataAxesGrid.XTitleFontFile = ''
vortexDynamicsxyzDisplay.DataAxesGrid.YTitleFontFile = ''
vortexDynamicsxyzDisplay.DataAxesGrid.ZTitleFontFile = ''
vortexDynamicsxyzDisplay.DataAxesGrid.XLabelFontFile = ''
vortexDynamicsxyzDisplay.DataAxesGrid.YLabelFontFile = ''
vortexDynamicsxyzDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
vortexDynamicsxyzDisplay.PolarAxes.PolarAxisTitleFontFile = ''
vortexDynamicsxyzDisplay.PolarAxes.PolarAxisLabelFontFile = ''
vortexDynamicsxyzDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
vortexDynamicsxyzDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 0.0, 10000.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
vortexDynamicsxyzDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(vortexDynamicsxyzDisplay, ('POINTS', 'Function0'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(densityLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
vortexDynamicsxyzDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
vortexDynamicsxyzDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Function0'
function0LUT = GetColorTransferFunction('Function0')
function0LUT.RGBPoints = [-1.5197791604976661e-05, 0.231373, 0.298039, 0.752941, -5.963111948670274e-17, 0.865003, 0.865003, 0.865003, 1.5197791604857399e-05, 0.705882, 0.0156863, 0.14902]
function0LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Function0'
function0PWF = GetOpacityTransferFunction('Function0')
function0PWF.Points = [-1.5197791604976661e-05, 0.0, 0.5, 0.0, 1.5197791604857399e-05, 1.0, 0.5, 0.0]
function0PWF.ScalarRangeInitialized = 1

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
function0LUT.ApplyPreset('X Ray', True)

# rescale color and/or opacity maps used to exactly fit the current data range
vortexDynamicsxyzDisplay.RescaleTransferFunctionToDataRange(False, True)

# create a new 'Contour'
contour1 = Contour(Input=vortexDynamicsxyz)
contour1.ContourBy = ['POINTS', 'Function1']
contour1.ComputeScalars = 1
contour1.Isosurfaces = np.linspace(5.e-19,2.0,20)
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1)

# get color transfer function/color map for 'Function1'
function1LUT = GetColorTransferFunction('Function1')
function1LUT.RGBPoints = [0.22134111111111107, 0.231373, 0.298039, 0.752941, 1.1067055555555556, 0.865003, 0.865003, 0.865003, 1.99207, 0.705882, 0.0156863, 0.14902]
function1LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'Function1']
contour1Display.LookupTable = function1LUT
contour1Display.OSPRayScaleArray = 'Function1'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'Momentum'
contour1Display.ScaleFactor = 1.5995794390945197
contour1Display.SelectScaleArray = 'Function1'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'Function1'
contour1Display.GaussianRadius = 0.07997897195472599
contour1Display.SetScaleArray = ['POINTS', 'Function1']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'Function1']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.SelectionCellLabelFontFile = ''
contour1Display.SelectionPointLabelFontFile = ''
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour1Display.OSPRayScaleFunction.Points = [-9.459567940017918e-06, 0.0, 0.5, 0.0, 1.53066275150749e-08, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [-9.459567940017918e-06, 0.0, 0.5, 0.0, 1.53066275150749e-08, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [-9.459567940017918e-06, 0.0, 0.5, 0.0, 1.53066275150749e-08, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
contour1Display.DataAxesGrid.XTitleFontFile = ''
contour1Display.DataAxesGrid.YTitleFontFile = ''
contour1Display.DataAxesGrid.ZTitleFontFile = ''
contour1Display.DataAxesGrid.XLabelFontFile = ''
contour1Display.DataAxesGrid.YLabelFontFile = ''
contour1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour1Display.PolarAxes.PolarAxisTitleFontFile = ''
contour1Display.PolarAxes.PolarAxisLabelFontFile = ''
contour1Display.PolarAxes.LastRadialAxisTextFontFile = ''
contour1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(vortexDynamicsxyz, renderView1)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'Function1'
function1PWF = GetOpacityTransferFunction('Function1')
function1PWF.Points = [0.22134111111111107, 0.0, 0.5, 0.0, 1.99207, 1.0, 0.5, 0.0]
function1PWF.ScalarRangeInitialized = 1

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
function1LUT.ApplyPreset('jet', True)

# rescale color and/or opacity maps used to exactly fit the current data range
contour1Display.RescaleTransferFunctionToDataRange(False, True)

# set active source
SetActiveSource(vortexDynamicsxyz)

# show data in view
vortexDynamicsxyzDisplay = Show(vortexDynamicsxyz, renderView1)

# show color bar/color legend
vortexDynamicsxyzDisplay.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 0.0, 10000.0]
renderView1.CameraParallelScale = 51.35179797039506

# Rescale transfer function
function0LUT.RescaleTransferFunction(-0.143467008462, 0.0874824541899)

# Rescale transfer function
function0PWF.RescaleTransferFunction(-0.143467008462, 0.0874824541899)

for k in range(N):
    # Properties modified on vortexDynamicsxyz
    vortexDynamicsxyz.FunctionFileName = FunctionFileNameList[k]

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.0, 0.0, 10000.0]
    renderView1.CameraParallelScale = 51.35179797039506

    # set active source
    SetActiveSource(contour1)

    # update the view to ensure updated data information
    renderView1.Update()

    # save screenshot
    SaveScreenshot(FigureFileNameList[k], renderView1, ImageResolution=[1555, 818])
