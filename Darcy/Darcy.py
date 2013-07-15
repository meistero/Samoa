try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

RenderView1 = CreateRenderView()
RenderView1.LightSpecularColor = [1.0, 1.0, 1.0]
RenderView1.UseOutlineForLODRendering = 0
RenderView1.KeyLightAzimuth = 10.0
RenderView1.UseTexturedBackground = 0
RenderView1.UseLight = 0
RenderView1.CameraPosition = [0.5, 0.5, 2.7320508075688776]
RenderView1.FillLightKFRatio = 3.0
RenderView1.Background2 = [0.0, 0.0, 0.165]
RenderView1.FillLightAzimuth = -10.0
RenderView1.LODResolution = 0.5
RenderView1.BackgroundTexture = []
RenderView1.InteractionMode = '2D'
RenderView1.StencilCapable = 1
RenderView1.LightIntensity = 1.0
RenderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
RenderView1.ImageReductionFactor = 2
RenderView1.CameraViewAngle = 30.0
RenderView1.CameraParallelScale = 0.571419922251881
RenderView1.EyeAngle = 2.0
RenderView1.HeadLightKHRatio = 3.0
RenderView1.StereoRender = 0
RenderView1.KeyLightIntensity = 0.75
RenderView1.BackLightAzimuth = 110.0
RenderView1.OrientationAxesInteractivity = 0
RenderView1.UseInteractiveRenderingForSceenshots = 0
RenderView1.UseOffscreenRendering = 0
RenderView1.Background = [1.0, 1.0, 1.0]
RenderView1.UseOffscreenRenderingForScreenshots = 0
RenderView1.NonInteractiveRenderDelay = 0.0
RenderView1.CenterOfRotation = [0.5, 0.5, 0.0]
RenderView1.CameraParallelProjection = 0
RenderView1.CompressorConfig = 'vtkSquirtCompressor 0 3'
RenderView1.HeadLightWarmth = 0.5
RenderView1.MaximumNumberOfPeels = 4
RenderView1.LightDiffuseColor = [1.0, 1.0, 1.0]
RenderView1.StereoType = 'Red-Blue'
RenderView1.DepthPeeling = 1
RenderView1.BackLightKBRatio = 3.5
RenderView1.StereoCapableWindow = 1
RenderView1.CameraViewUp = [0.0, 1.0, 0.0]
RenderView1.LightType = 'HeadLight'
RenderView1.LightAmbientColor = [1.0, 1.0, 1.0]
RenderView1.RemoteRenderThreshold = 3.0
RenderView1.CacheKey = 19.0
RenderView1.UseCache = 0
RenderView1.KeyLightElevation = 50.0
RenderView1.CenterAxesVisibility = 0
RenderView1.MaintainLuminance = 0
RenderView1.StillRenderImageReductionFactor = 1
RenderView1.BackLightWarmth = 0.5
RenderView1.FillLightElevation = -75.0
RenderView1.MultiSamples = 0
RenderView1.FillLightWarmth = 0.4
RenderView1.AlphaBitPlanes = 1
RenderView1.LightSwitch = 1
RenderView1.OrientationAxesVisibility = 0
RenderView1.CameraClippingRange = [2.704730299493189, 2.7730315696824106]
RenderView1.BackLightElevation = 0.0
RenderView1.ViewTime = 19.0
RenderView1.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]
RenderView1.LODThreshold = 5.0
RenderView1.CollectGeometryThreshold = 100.0
RenderView1.UseGradientBackground = 0
RenderView1.KeyLightWarmth = 0.6
RenderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]

darcy_20130127_145639_ = XMLUnstructuredGridReader( guiName="darcy_20130127_145639_*", PointArrayStatus=['pressure', 'saturation'], CellArrayStatus=['permeability', 'velocity', 'grid depth', 'refinement flag'], FileName=['/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_0.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_1.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_2.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_3.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_4.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_5.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_6.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_7.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_8.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_9.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_10.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_11.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_12.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_13.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_14.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_15.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_16.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_17.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_18.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_19.vtu', '/home/meistero/Desktop/Samoa/output/darcy_20130127_145639_20.vtu'] )

a1_pressure_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1000000.0, 1.0, 0.5, 0.0] )

a1_saturation_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

a1_pressure_PVLookupTable = GetLookupTableForArray( "pressure", 1, Discretize=1, RGBPoints=[0.0, 0.0, 0.0, 1.0, 1000000.0, 1.0, 0.0, 0.0], UseLogScale=0, VectorComponent=0, NanColor=[0.4980392156862745, 0.4980392156862745, 0.4980392156862745], NumberOfTableValues=16, EnableOpacityMapping=0, ColorSpace='HSV', IndexedLookup=0, VectorMode='Magnitude', ScalarOpacityFunction=a1_pressure_PiecewiseFunction, HSVWrap=0, ScalarRangeInitialized=1.0, AllowDuplicateScalars=1, Annotations=[], LockScalarRange=1 )

a1_saturation_PVLookupTable = GetLookupTableForArray( "saturation", 1, Discretize=1, RGBPoints=[0.0, 0.3176470588235294, 0.3411764705882353, 0.43137254901960786, 0.17, 0.0, 0.0, 1.0, 0.34, 0.0, 1.0, 1.0, 0.5, 0.0, 1.0, 0.0, 0.67, 1.0, 1.0, 0.0, 0.84, 1.0, 0.0, 0.0, 1.0, 0.8784313725490196, 0.0, 1.0], UseLogScale=0, VectorComponent=0, NanColor=[1.0, 1.0, 0.0], NumberOfTableValues=256, EnableOpacityMapping=0, ColorSpace='RGB', IndexedLookup=0, VectorMode='Magnitude', ScalarOpacityFunction=a1_saturation_PiecewiseFunction, HSVWrap=0, ScalarRangeInitialized=1.0, AllowDuplicateScalars=0, Annotations=[], LockScalarRange=0 )

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
