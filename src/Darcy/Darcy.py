try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

darcy = GetActiveSource()

Threshold2 = Threshold()
Threshold2.Scalars = ['CELLS', 'permeability']
Threshold2.ThresholdRange = [1e-08,1e20]

DataRepresentation0 = GetDisplayProperties(darcy)
DataRepresentation0.Visibility = 0

RenderView1 = GetRenderView()
RenderView1.InteractionMode = '3D'
RenderView1.ViewTime = 0.0
RenderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

a1_pressure_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1000000.0, 1.0, 0.5, 0.0] )

a1_saturation_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

a1_pressure_PVLookupTable = GetLookupTableForArray( "pressure", 1, NanColor=[1.0, 1.0, 0.0], RGBPoints=[5007.231157965089, 0.3176470588235294, 0.3411764705882353, 0.43137254901960786, 5645.436319215898, 0.0, 0.0, 1.0, 6283.641480466708, 0.0, 1.0, 1.0, 6884.305161643941, 0.0, 1.0, 0.0, 7522.51032289475, 1.0, 1.0, 0.0, 8160.7154841455595, 1.0, 0.0, 0.0, 8761.379165322793, 0.8784313725490196, 0.0, 1.0], ColorSpace='RGB' )

a1_saturation_PVLookupTable = GetLookupTableForArray( "saturation", 1, NanColor=[0.4980392156862745, 0.4980392156862745, 0.4980392156862745], RGBPoints=[0.1, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0], UseLogScale=0, ColorSpace='HSV', AllowDuplicateScalars=1 )

ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='pressure', Position2=[0.13, 0.5], TitleOpacity=1.0, TitleShadow=0, AutomaticLabelFormat=1, TitleFontSize=12, TitleColor=[0.0, 0.0, 0.0], AspectRatio=20.0, NumberOfLabels=5, ComponentTitle='', Resizable=1, TitleFontFamily='Arial', Visibility=0, LabelFontSize=12, LabelFontFamily='Arial', TitleItalic=0, Selectable=0, LabelItalic=0, Enabled=0, LabelColor=[0.0, 0.0, 0.0], Position=[0.87, 0.25], LabelBold=0, UseNonCompositedRenderer=1, LabelOpacity=1.0, TitleBold=0, LabelFormat='%-#6.3g', Orientation='Vertical', LabelShadow=0, LookupTable=a1_pressure_PVLookupTable, Repositionable=1 )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)

DataRepresentation1 = Show()
DataRepresentation1.CubeAxesZAxisVisibility = 1
DataRepresentation1.SelectionPointLabelColor = [0.5, 0.5, 0.5]
DataRepresentation1.SelectionPointFieldDataArrayName = 'pressure'
DataRepresentation1.SuppressLOD = 0
DataRepresentation1.CubeAxesXGridLines = 0
DataRepresentation1.CubeAxesYAxisTickVisibility = 1
DataRepresentation1.Position = [0.0, 0.0, 0.0]
DataRepresentation1.BackfaceRepresentation = 'Follow Frontface'
DataRepresentation1.SelectionOpacity = 1.0
DataRepresentation1.SelectionPointLabelShadow = 0
DataRepresentation1.CubeAxesYGridLines = 0
DataRepresentation1.CubeAxesZAxisTickVisibility = 1
DataRepresentation1.OrientationMode = 'Direction'
DataRepresentation1.Source.TipResolution = 6
DataRepresentation1.ScaleMode = 'No Data Scaling Off'
DataRepresentation1.Diffuse = 1.0
DataRepresentation1.SelectionUseOutline = 0
DataRepresentation1.SelectionPointLabelFormat = ''
DataRepresentation1.CubeAxesZTitle = 'Z-Axis'
DataRepresentation1.Specular = 0.1
DataRepresentation1.SelectionVisibility = 1
DataRepresentation1.InterpolateScalarsBeforeMapping = 1
DataRepresentation1.CustomRangeActive = [0, 0, 0]
DataRepresentation1.Origin = [0.0, 0.0, 0.0]
DataRepresentation1.Source.TipLength = 0.35
DataRepresentation1.CubeAxesVisibility = 0
DataRepresentation1.Scale = [1.0, 1.0, 1.0]
DataRepresentation1.SelectionCellLabelJustification = 'Left'
DataRepresentation1.DiffuseColor = [0.6666666666666666, 0.6666666666666666, 0.4980392156862745]
DataRepresentation1.SelectionCellLabelOpacity = 1.0
DataRepresentation1.CubeAxesInertia = 1
DataRepresentation1.Source = "Arrow"
DataRepresentation1.Source.Invert = 0
DataRepresentation1.Masking = 0
DataRepresentation1.Opacity = 1.0
DataRepresentation1.LineWidth = 1.0
DataRepresentation1.MeshVisibility = 0
DataRepresentation1.Visibility = 1
DataRepresentation1.SelectionCellLabelFontSize = 18
DataRepresentation1.CubeAxesCornerOffset = 0.0
DataRepresentation1.SelectionPointLabelJustification = 'Left'
DataRepresentation1.OriginalBoundsRangeActive = [0, 0, 0]
DataRepresentation1.SelectionPointLabelVisibility = 0
DataRepresentation1.SelectOrientationVectors = ''
DataRepresentation1.CubeAxesTickLocation = 'Inside'
DataRepresentation1.BackfaceDiffuseColor = [0.6666666666666666, 0.6666666666666666, 0.4980392156862745]
DataRepresentation1.CubeAxesYAxisVisibility = 1
DataRepresentation1.SelectionPointLabelFontFamily = 'Arial'
DataRepresentation1.Source.ShaftResolution = 6
DataRepresentation1.CubeAxesUseDefaultYTitle = 1
DataRepresentation1.SelectScaleArray = ''
DataRepresentation1.CubeAxesYTitle = 'Y-Axis'
DataRepresentation1.ColorAttributeType = 'POINT_DATA'
DataRepresentation1.AxesOrigin = [0.0, 0.0, 0.0]
DataRepresentation1.UserTransform = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
DataRepresentation1.SpecularPower = 100.0
DataRepresentation1.Texture = []
DataRepresentation1.SelectionCellLabelShadow = 0
DataRepresentation1.AmbientColor = [1.0, 1.0, 1.0]
DataRepresentation1.MapScalars = 1
DataRepresentation1.PointSize = 2.0
DataRepresentation1.CubeAxesUseDefaultXTitle = 1
DataRepresentation1.SelectionCellLabelFormat = ''
DataRepresentation1.Scaling = 0
DataRepresentation1.StaticMode = 0
DataRepresentation1.SelectionCellLabelColor = [0.0, 1.0, 0.0]
DataRepresentation1.Source.TipRadius = 0.1
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.0]
DataRepresentation1.CubeAxesXAxisTickVisibility = 1
DataRepresentation1.SelectionCellLabelVisibility = 0
DataRepresentation1.NonlinearSubdivisionLevel = 1
DataRepresentation1.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation1.Representation = 'Surface'
DataRepresentation1.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
DataRepresentation1.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
DataRepresentation1.Orientation = [0.0, 0.0, 0.0]
DataRepresentation1.CubeAxesXTitle = 'X-Axis'
DataRepresentation1.ScalarOpacityUnitDistance = 0.06435100048451363
DataRepresentation1.BackfaceOpacity = 1.0
DataRepresentation1.SelectionPointLabelFontSize = 18
DataRepresentation1.SelectionCellFieldDataArrayName = 'grid depth'
DataRepresentation1.SelectionColor = [1.0, 0.0, 1.0]
DataRepresentation1.Ambient = 0.0
DataRepresentation1.CubeAxesXAxisMinorTickVisibility = 1
DataRepresentation1.ScaleFactor = 0.1
DataRepresentation1.BackfaceAmbientColor = [1.0, 1.0, 1.0]
DataRepresentation1.Source.ShaftRadius = 0.03
DataRepresentation1.ScalarOpacityFunction = a1_saturation_PiecewiseFunction
DataRepresentation1.SelectMaskArray = ''
DataRepresentation1.SelectionLineWidth = 2.0
DataRepresentation1.CubeAxesZAxisMinorTickVisibility = 1
DataRepresentation1.CubeAxesXAxisVisibility = 1
DataRepresentation1.Interpolation = 'Gouraud'
DataRepresentation1.SelectMapper = 'Projected tetra'
DataRepresentation1.SelectionCellLabelFontFamily = 'Arial'
DataRepresentation1.SelectionCellLabelItalic = 0
DataRepresentation1.CubeAxesYAxisMinorTickVisibility = 1
DataRepresentation1.CubeAxesZGridLines = 0
DataRepresentation1.ExtractedBlockIndex = 0
DataRepresentation1.SelectionPointLabelOpacity = 1.0
DataRepresentation1.UseAxesOrigin = 0
DataRepresentation1.CubeAxesFlyMode = 'Closest Triad'
DataRepresentation1.Pickable = 1
DataRepresentation1.CustomBoundsActive = [0, 0, 0]
DataRepresentation1.CubeAxesGridLineLocation = 'All Faces'
DataRepresentation1.SelectionRepresentation = 'Wireframe'
DataRepresentation1.SelectionPointLabelBold = 0
DataRepresentation1.ColorArrayName = 'saturation'
DataRepresentation1.SelectionPointLabelItalic = 0
DataRepresentation1.AllowSpecularHighlightingWithScalarColoring = 0
DataRepresentation1.SpecularColor = [1.0, 1.0, 1.0]
DataRepresentation1.CubeAxesUseDefaultZTitle = 1
DataRepresentation1.LookupTable = a1_saturation_PVLookupTable
DataRepresentation1.SelectionPointSize = 5.0
DataRepresentation1.SelectionCellLabelBold = 0
DataRepresentation1.Orient = 0

Render()
