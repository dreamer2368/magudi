# trace generated using paraview version 5.6.0-RC1

#### import the simple module from the paraview
from paraview.simple import *
from os import *
import numpy as np
chdir('/Users/Kevin/Desktop/VortexDynamics')

prefix = 'VortexDynamics'

gamma = 1.4
p0 = 1./gamma
R = 1./0.15
wavelength = 52.5 * R
Radius = 0.5 * wavelength

N=421
QFileNameList=[]
for i in range(0,N):
	istr = i*200
	QFileNameList=QFileNameList+['%s-%08d.q'%(prefix,istr)]

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PLOT3D Reader'
vortexDynamicsxyz = PLOT3DReader(FileName='%s.xyz'%prefix,
    QFileName=QFileNameList,
    FunctionFileName='')
vortexDynamicsxyz.Functions = [110]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1555, 818]

# show data in view
vortexDynamicsxyzDisplay = Show(vortexDynamicsxyz, renderView1)

# get color transfer function/color map for 'Density'
densityLUT = GetColorTransferFunction('Density')
densityLUT.RGBPoints = [0.6184731892655716, 0.231373, 0.298039, 0.752941, 0.8105310052124406, 0.865003, 0.865003, 0.865003, 1.0025888211593097, 0.705882, 0.0156863, 0.14902]
densityLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Density'
densityPWF = GetOpacityTransferFunction('Density')
densityPWF.Points = [0.6184731892655716, 0.0, 0.5, 0.0, 1.0025888211593097, 1.0, 0.5, 0.0]
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

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
vortexDynamicsxyzDisplay.SetScalarBarVisibility(renderView1, True)

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Calculator'
calculator1 = Calculator(Input=vortexDynamicsxyz)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'p'
calculator1.Function = '(' + str(gamma) + '-1)*(StagnationEnergy-0.5/Density*( Momentum_X^2 + Momentum_Y^2 ) )'

# create a new 'Slice'
slice1 = Slice(Input=calculator1)

# Properties modified on slice1
slice1.SliceType = 'Cylinder'

# Properties modified on slice1.SliceType
slice1.SliceType.Axis = [0.0, 0.0, 1.0]
slice1.SliceType.Radius = Radius

# show data in view
slice1Display = Show(slice1, renderView1)

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'Density']
slice1Display.LookupTable = densityLUT
slice1Display.OSPRayScaleArray = 'Density'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'Momentum'
slice1Display.ScaleFactor = 59.994041067421655
slice1Display.SelectScaleArray = 'Density'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'Density'
slice1Display.GaussianRadius = 2.9997020533710828
slice1Display.SetScaleArray = ['POINTS', 'Density']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'Density']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.SelectionCellLabelFontFile = ''
slice1Display.SelectionPointLabelFontFile = ''
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [-9.459567940017918e-06, 0.0, 0.5, 0.0, 1.53066275150749e-08, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [-9.459567940017918e-06, 0.0, 0.5, 0.0, 1.53066275150749e-08, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [-9.459567940017918e-06, 0.0, 0.5, 0.0, 1.53066275150749e-08, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
slice1Display.DataAxesGrid.XTitleFontFile = ''
slice1Display.DataAxesGrid.YTitleFontFile = ''
slice1Display.DataAxesGrid.ZTitleFontFile = ''
slice1Display.DataAxesGrid.XLabelFontFile = ''
slice1Display.DataAxesGrid.YLabelFontFile = ''
slice1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
slice1Display.PolarAxes.PolarAxisTitleFontFile = ''
slice1Display.PolarAxes.PolarAxisLabelFontFile = ''
slice1Display.PolarAxes.LastRadialAxisTextFontFile = ''
slice1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(vortexDynamicsxyz, renderView1)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=slice1)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# get layout
layout1 = GetLayout()

# place view in the layout
layout1.AssignView(2, spreadSheetView1)

# show data in view
integrateVariables1Display = Show(integrateVariables1, spreadSheetView1)

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
spreadSheetView1.Update()

volume = 2.0 * np.pi * Radius
pressure = np.zeros((N,),dtype=np.double)
tsteps = vortexDynamicsxyz.TimestepValues
dt = 0.04 * 200 * 0.15
time = np.arange(N) * dt
for k in range(N-1):
    animationScene1.GoToNext()
    spreadSheetView1.Update()
    test = paraview.servermanager.Fetch(integrateVariables1)
    pressure[k+1] = test.GetPointData().GetArray('p').GetValue(0)
    print (k)
pressure /= volume
pressure -= p0

data = np.zeros([N,2],dtype=np.double)
data[:,0], data[:,1] = time, pressure
np.savetxt('R' + str(int(Radius)) + '.p.txt',data,fmt='%.16E',delimiter='\t')

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
