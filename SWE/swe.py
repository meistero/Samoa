try: paraview.simple
except: from paraview.simple import *

swe = GetActiveSource()

DataRepresentation1 = GetDisplayProperties(swe)
DataRepresentation1.Visibility = 0

CleantoGrid1 = CleantoGrid()
CellDatatoPointData1 = CellDatatoPointData()

WarpByScalar1 = WarpByScalar()
WarpByScalar1.Scalars = ['POINTS', 'bathymetry']
WarpByScalar1.ScaleFactor = 5.0e-5

a1_bathymetry_PVLookupTable = GetLookupTableForArray( "bathymetry", 1, RGBPoints=[-8000, 0.23137254901960785, 0.2980392156862745, 0.7529411764705882, 0, 0.3333333333333333, 0.0, 0.0, 1300, 0.0, 0.6666666666666666, 0.0], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='RGB', HSVWrap=0, ScalarRangeInitialized=1.0 )

DataRepresentation4 = Show()
DataRepresentation4.EdgeColor = [0.0, 0.0, 0.0]
DataRepresentation4.SelectionPointFieldDataArrayName = 'bathymetry'
DataRepresentation4.DiffuseColor = [0.7568627450980392, 0.7568627450980392, 0.7568627450980392]
DataRepresentation4.SelectionCellFieldDataArrayName = 'bathymetry'
DataRepresentation4.ScalarOpacityUnitDistance = 0.05806444328815579
DataRepresentation4.BackfaceDiffuseColor = [0.7568627450980392, 0.7568627450980392, 0.7568627450980392]
DataRepresentation4.ScaleFactor = 0.1
DataRepresentation4.ColorArrayName = 'bathymetry'
DataRepresentation4.LookupTable = a1_bathymetry_PVLookupTable

SetActiveSource(CellDatatoPointData1)

WarpByScalar2 = WarpByScalar()
WarpByScalar2.Scalars = ['POINTS', 'water height']
WarpByScalar2.ScaleFactor = 5.0e-5

a1_waterheight_PVLookupTable = GetLookupTableForArray( "water height", 1, RGBPoints=[-0.01, 0.27058823529411763, 0.0, 0.6745098039215687, 0.0, 0.10196078431372549, 0.34509803921568627, 1.0, 0.01, 1.0, 1.0, 1.0], ColorSpace='RGB', LockScalarRange=1)
a1_waterheight_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

DataRepresentation5 = Show()
DataRepresentation5.EdgeColor = [0.0, 0.0, 0.0]
DataRepresentation5.SelectionPointFieldDataArrayName = 'water height'
DataRepresentation5.Opacity = 0.9
DataRepresentation5.DiffuseColor = [0.7568627450980392, 0.7568627450980392, 0.7568627450980392]
DataRepresentation5.SelectionCellFieldDataArrayName = 'water height'
DataRepresentation5.ScalarOpacityUnitDistance = 0.057107503012105605
DataRepresentation5.BackfaceDiffuseColor = [0.7568627450980392, 0.7568627450980392, 0.7568627450980392]
DataRepresentation5.ScaleFactor = 0.1
DataRepresentation5.ScalarOpacityFunction = a1_waterheight_PiecewiseFunction
DataRepresentation5.ColorArrayName = 'water height'
DataRepresentation5.LookupTable = a1_waterheight_PVLookupTable


RenderView1 = GetRenderView()
RenderView1.ViewTime = 0.0
RenderView1.CameraViewUp = [0.004157728385740771, 0.4591183137073247, 0.8883654019114056]
RenderView1.CameraPosition = [0.47392100562411804, -1.1859719410739753, 0.8161608127862509]
RenderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
RenderView1.CameraClippingRange = [0.8613699850146684, 3.2076058865387793]


Render()
