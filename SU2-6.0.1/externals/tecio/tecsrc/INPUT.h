/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2010 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/
#if defined EXTERN
#undef EXTERN
#endif
#if defined INITMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

/* Input Specification limits */

/* General */
EXTERN InputSpec_s          /*X*/  GridCoordInputSpec;
EXTERN InputSpec_s          /*X*/  GridCoordFloatInputSpec;
EXTERN InputSpec_s          /*X*/  XFrameCoordInputSpec;
EXTERN InputSpec_s          /*X*/  YFrameCoordInputSpec;
EXTERN InputSpec_s          /*X*/  XFrameCoordFloatInputSpec;
EXTERN InputSpec_s          /*X*/  YFrameCoordFloatInputSpec;
EXTERN InputSpec_s          /*X*/  XFrameCoordDeltaInputSpec;
EXTERN InputSpec_s          /*X*/  YFrameCoordDeltaInputSpec;
EXTERN InputSpec_s          /*X*/  XFrameCoordFloatDeltaInputSpec;
EXTERN InputSpec_s          /*X*/  YFrameCoordFloatDeltaInputSpec;
EXTERN InputSpec_s          /*X*/  FrameOffsetCoordInputSpec;
EXTERN InputSpec_s          /*X*/  XPaperCoordInputSpec;
EXTERN InputSpec_s          /*X*/  YPaperCoordInputSpec;
EXTERN InputSpec_s          /*X*/  AxisPercentageInputSpec;
EXTERN InputSpec_s          /*X*/  AngleInputSpec;
EXTERN InputSpec_s          /*X*/  AngleToApproxInputSpec;
EXTERN InputSpec_s          /*X*/  FieldOfViewInputSpec;
EXTERN InputSpec_s          /*X*/  ZeroAndAboveLgIndexInputSpec;
EXTERN InputSpec_s          /*X*/  ZeroAndAboveSmIntegerInputSpec;
EXTERN InputSpec_s          /*X*/  ZeroAndAboveDoubleInputSpec;
EXTERN InputSpec_s          /*X*/  AboveZeroLgIndexInputSpec;
EXTERN InputSpec_s          /*X*/  AboveZeroDoubleInputSpec;
EXTERN InputSpec_s          /*X*/  DoubleInputSpec;
EXTERN InputSpec_s          /*X*/  EntIndexInputSpec;
EXTERN InputSpec_s          /*X*/  EntRangeInputSpec;
EXTERN InputSpec_s          /*X*/  IndexRangeInputSpec;
EXTERN InputSpec_s          /*X*/  AboveZeroIndexRangeInputSpec;
EXTERN InputSpec_s          /*X*/  ZeroToOneInputSpec;
EXTERN InputSpec_s          /*X*/  PercentageInputSpec;
EXTERN InputSpec_s          /*X*/  AboveZeroPercentageInputSpec;
EXTERN InputSpec_s          /*X*/  SignedPercentageInputSpec;
EXTERN InputSpec_s          /*X*/  RadiansInputSpec;
EXTERN InputSpec_s          /*X*/  AboveZeroRadiansInputSpec;
EXTERN InputSpec_s          /*X*/  TimeDateDoubleInputSpec;
EXTERN InputSpec_s          /*X*/  AboveZeroTimeDateDoubleInputSpec;
EXTERN InputSpec_s          /*X*/  AboveZeroElapsedTimeInputSpec;


/* Specific */
#define MIN_VIEWPORT_SIZE 0.05
EXTERN InputSpec_s          /*X*/  SurfaceTranslucencyInputSpec;
EXTERN InputSpec_s          /*X*/  MaxDepthBufferSizeInputSpec;
EXTERN InputSpec_s          /*X*/  MaxMultiSamplesInputSpec;
EXTERN InputSpec_s          /*X*/  MaxAccumBufferSizeInputSpec;
EXTERN InputSpec_s          /*X*/  MinBitsPerRGBPlaneInputSpec;
EXTERN InputSpec_s          /*X*/  AnimationSpeedInputSpec;
EXTERN InputSpec_s          /*X*/  AnimationNumStepsInputSpec;
EXTERN InputSpec_s          /*X*/  MaxCustomColorsInInterfaceInputSpec;
EXTERN InputSpec_s          /*X*/  MaxReducedPointsInputSpec;
EXTERN InputSpec_s          /*X*/  MaxStripLengthInputSpec;
EXTERN InputSpec_s          /*X*/  MaxPrimativesPerBlockInputSpec;
EXTERN InputSpec_s          /*X*/  MaxTextureSizeInputSpec;
EXTERN InputSpec_s          /*X*/  SuperSampleFactorInputSpec;
EXTERN InputSpec_s          /*X*/  TickLengthInputSpec;
EXTERN InputSpec_s          /*X*/  BorrowLicenseInputSpec;




/* I/O Related */
EXTERN InputSpec_s          /*X*/  HardcopyPaperSizeInputSpec;
EXTERN InputSpec_s          /*X*/  HardcopyNumCopiesInputSpec;
EXTERN InputSpec_s          /*X*/  HardcopyPrecisionInputSpec;
EXTERN InputSpec_s          /*X*/  HardcopyPenSpeedInputSpec;
EXTERN InputSpec_s          /*X*/  PenPlotterPenNumberInputSpec;
EXTERN InputSpec_s          /*X*/  BitDumpDepthInputSpec;


/* Widths, physical lengths, etc. */
EXTERN InputSpec_s          /*X*/  XFrameDimensionInputSpec;
EXTERN InputSpec_s          /*X*/  YFrameDimensionInputSpec;
EXTERN InputSpec_s          /*X*/  LineThicknessInputSpec;
EXTERN InputSpec_s          /*X*/  PatternLengthInputSpec;
EXTERN InputSpec_s          /*X*/  AxisPercentageTextSizeInputSpec;
EXTERN InputSpec_s          /*X*/  FrameTextSizeInputSpec;
EXTERN InputSpec_s          /*X*/  GridTextSizeInputSpec;
EXTERN InputSpec_s          /*X*/  PointTextSizeInputSpec;
EXTERN InputSpec_s          /*X*/  TextBoxMarginInputSpec;
EXTERN InputSpec_s          /*X*/  TextLineSpacingInputSpec;
EXTERN InputSpec_s          /*X*/  ArrowheadSizeInputSpec;
EXTERN InputSpec_s          /*X*/  AxisLabelOffsetInputSpec;
EXTERN InputSpec_s          /*X*/  LegendLineSpacingInputSpec;
EXTERN InputSpec_s          /*X*/  StreamStepSizeInputSpec;
EXTERN InputSpec_s          /*X*/  StreamMaxStepsInputSpec;
EXTERN InputSpec_s          /*X*/  ArrowheadSpacingInputSpec;
EXTERN InputSpec_s          /*X*/  RulerPaddingInputSpec;
EXTERN InputSpec_s          /*X*/  RulerThicknessInputSpec;
EXTERN InputSpec_s          /*X*/  PickHandleWidthInputSpec;
EXTERN InputSpec_s          /*X*/  ImageDimensionInputSpec;
EXTERN InputSpec_s          /*X*/  ZoomScalePerFrameUnitInputSpec;
EXTERN InputSpec_s          /*X*/  RGBLegendHeightInputSpec;



/* Limit the number of objects or limit which object can be selected*/
EXTERN InputSpec_s          /*X*/  ColorMapGroupInputSpec;
EXTERN InputSpec_s          /*X*/  SliceGroupInputSpec;
EXTERN InputSpec_s          /*X*/  IsoSurfaceGroupInputSpec;
EXTERN InputSpec_s          /*X*/  ContourGroupInputSpec;
EXTERN InputSpec_s          /*X*/  ColorIndexInputSpec;
EXTERN InputSpec_s          /*X*/  NumLightSourceShadesInputSpec;
EXTERN InputSpec_s          /*X*/  NumberOfControlPointsInputSpec;
EXTERN InputSpec_s          /*X*/  CustomLabelNumberInputSpec;
EXTERN InputSpec_s          /*X*/  NumMinorTicksInputSpec;
EXTERN InputSpec_s          /*X*/  AxisEdgeNumberInputSpec;
EXTERN InputSpec_s          /*X*/  LineMapWhichXAxisInputSpec;
EXTERN InputSpec_s          /*X*/  LineMapWhichYAxisInputSpec;
EXTERN InputSpec_s          /*X*/  NumberOfCurvePointsInputSpec;
EXTERN InputSpec_s          /*X*/  NumberOfContourLevelsInputSpec;
EXTERN InputSpec_s          /*X*/  ColorMapOverrideLevelInputSpec;
EXTERN InputSpec_s          /*X*/  ColorMapOverrideNumberInputSpec;
EXTERN InputSpec_s          /*X*/  NumberOfColorMapCyclesInputSpec;
EXTERN InputSpec_s          /*X*/  NumberOfRodPointsInputSpec;
EXTERN InputSpec_s          /*X*/  NumberOfStreamtracesInputSpec;
EXTERN InputSpec_s          /*X*/  NumberOfEllipsePointsInputSpec;
EXTERN InputSpec_s          /*X*/  MaxPtsInALineInputSpec;
EXTERN InputSpec_s          /*X*/  MaxChrsTextLabelsInputSpec;
EXTERN InputSpec_s          /*X*/  MaxContourLevelsInputSpec;
EXTERN InputSpec_s          /*X*/  MaxLinkGroupsInputSpec;


/* Ratios */
EXTERN InputSpec_s          /*X*/  DataAspectRatioLimitInputSpec;
EXTERN InputSpec_s          /*X*/  DataAspectRatioResetInputSpec;
EXTERN InputSpec_s          /*X*/  AxisBoxAspectRatioLimitInputSpec;
EXTERN InputSpec_s          /*X*/  AxisBoxAspectRatioResetInputSpec;
EXTERN InputSpec_s          /*X*/  AxisRatioInputSpec;
EXTERN InputSpec_s          /*X*/  AxisBoxPaddingInputSpec;
EXTERN InputSpec_s          /*X*/  ScreenDistanceRatioInputSpec;
EXTERN InputSpec_s          /*X*/  LiftFractionInputSpec;
EXTERN InputSpec_s          /*X*/  ZClipInputSpec;
EXTERN InputSpec_s          /*X*/  VectorHeadSizeFractionInputSpec;


/* Misc */
EXTERN InputSpec_s          /*X*/  ValuePrecisionInputSpec;
EXTERN InputSpec_s          /*X*/  PolynomialOrderInputSpec;
EXTERN InputSpec_s          /*X*/  SplineSlopeInputSpec;
EXTERN InputSpec_s          /*X*/  RotationStepSizeInputSpec;
EXTERN InputSpec_s          /*X*/  SmoothRotationDegPerFrameUnitInputSpec;
EXTERN InputSpec_s          /*X*/  TranslationStepSizeInputSpec;
EXTERN InputSpec_s          /*X*/  ScaleStepSizeInputSpec;
EXTERN InputSpec_s          /*X*/  SortLevelInputSpec;
EXTERN InputSpec_s          /*X*/  AxisLabelSkipInputSpec;
EXTERN InputSpec_s          /*X*/  TextAngleInputSpec;
EXTERN InputSpec_s          /*X*/  ArrowheadAngleInputSpec;
EXTERN InputSpec_s          /*X*/  MinCreaseAngleInputSpec;
EXTERN InputSpec_s          /*X*/  ExponentInputSpec;
EXTERN InputSpec_s          /*X*/  SmoothWeightInputSpec;
EXTERN InputSpec_s          /*X*/  TriangleKeepFactorInputSpec;
EXTERN InputSpec_s          /*X*/  PlotAttrColumnWidthInputSpec;
EXTERN InputSpec_s          /*X*/  ImageQualityInputSpec;

