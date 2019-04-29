from paraview.simple import *
from os import *
#chdir('/Volumes/Seagate Backup Plus Drive/NaSt_adjoint/MultiblockJet')

prefix = 'MultiblockJet'

N=201
QFileNameList=[]
FunctionFileNameList=[]
FigureList=[]
for i in range(0,N):
	istr = 457000+i*500
	QFileNameList=QFileNameList+['%s-%08d.adjoint.q'%(prefix,istr)]
#	FunctionFileNameList=FunctionFileNameList+['MultiblockJet-00'+istr+'.vorticity_dilatation.f']
	FigureList=FigureList+['%s-%08d.adjoint.png'%(prefix,istr)]

reader = PLOT3DReader(FileName='%s.xyz'%prefix,QFileName=QFileNameList)
view = GetActiveView()
if view is None:
    view=CreateRenderView()
view.OrientationAxesVisibility=0
if (view.ViewSize[1]<800):
    view.ViewSize[1]=800
    view.ViewSize=[int(view.ViewSize[1]*43./25.),view.ViewSize[1]]

calculator = Calculator(reader)
calculator.ResultArrayName = 'SensitivityAmplitude'
calculator.Function = 'abs(StagnationEnergy)'
UpdatePipeline()

slice=Slice(calculator)
SetProperties(slice,SliceType='Plane')
slice.SliceType.Origin=[0,0,0]
slice.SliceType.Normal=[1,0,0]
dp = GetDisplayProperties()
dp.Representation = 'Surface'
dp.SetPropertyWithName('ColorArrayName',['POINTS','SensitivityAmplitude'])
#For some reason, color map for Function0 does not work.. used Density instead
Function0ColorMap = GetColorTransferFunction('SensitivityAmplitude')
Function0ColorMap.ApplyPreset('Blue to Red Rainbow',True)
Function0ColorMap.RescaleTransferFunction(1e-10,1e6)
Function0ColorMap.UseLogScale = 1
#dp.RescaleTransferFunctionToDataRange()
dp.SetScalarBarVisibility(view,False)
#scalarBar = GetScalarBar(Function0ColorMap,view)
#scalarBar.LabelColor = [0.0, 0.0, 0.0]
#scalarBar.LabelFontSize = 35
#scalarBar.LabelBold = 1
#scalarBar.TitleColor = [0.0, 0.0, 0.0]
#scalarBar.TitleFontSize = 35
#scalarBar.TitleBold = 1
#scalarBar.LabelFormat = '%-#6.3e'
#scalarBar.RangeLabelFormat = '%-#6.3e'
#scalarBar.CustomLabels = [1e-10,1e-6,1e-2,1e2,1e6]
#scalarBar.ScalarBarLength, scalarBar.ScalarBarThickness = 0.9, 32

#Time annotation
annTime = AnnotateTimeFilter(reader)
annTime.Scale = - 1.2e-3 * 500
annTime.Shift = 120.0
annTime.Format = 't = %6.2f *D/a'
dpTime = GetDisplayProperties(annTime)
dpTime.FontSize = 35
dpTime.Bold = 1
dpTime.WindowLocation = 'UpperLeftCorner'
dpTime.Color = [0.0, 0.0, 0.0]
Show(annTime)

#Recommend to adjust camera setting after Render()
Show(slice)
Render()
camera = GetActiveCamera()
camera.SetPosition(-47,0,12.5)
camera.SetFocalPoint(0,0,12.5)

dp.ColorArrayName='SensitivityAmplitude'
SetActiveSource(reader)
Render()

tsteps=reader.TimestepValues
for i in range(0,N):
    view.ViewTime=tsteps[N-1-i]
    Render()
    WriteImage(FigureList[N-1-i])


