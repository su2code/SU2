/*
 * TECXXX.h: Copyright (C) 1988-2010 Tecplot, Inc.
 */

#if !defined TECXXX_H_
#define TECXXX_H_

#if !defined CRAY
#  define TECFOREIGN112 tecforeign112
#  define TECINI112     tecini112
#  define TECZNE112     teczne112
#  define TECDAT112     tecdat112
#  define TECNOD112     tecnod112
#  define TECNODE112     tecnode112
#  define TECGEO112     tecgeo112
#  define TECTXT112     tectxt112
#  define TECLAB112     teclab112
#  define TECFIL112     tecfil112
#  define TECEND112     tecend112
#  define TECUSR112     tecusr112
#  define TECAUXSTR112  tecauxstr112
#  define TECZAUXSTR112 teczauxstr112
#  define TECVAUXSTR112 tecvauxstr112
#  define TECFACE112    tecface112
#  define TECPOLY112    tecpoly112

#  define TECFOREIGN111 tecforeign111
#  define TECINI111     tecini111
#  define TECZNE111     teczne111
#  define TECDAT111     tecdat111
#  define TECNOD111     tecnod111
#  define TECGEO111     tecgeo111
#  define TECTXT111     tectxt111
#  define TECLAB111     teclab111
#  define TECFIL111     tecfil111
#  define TECEND111     tecend111
#  define TECUSR111     tecusr111
#  define TECAUXSTR111  tecauxstr111
#  define TECZAUXSTR111 teczauxstr111
#  define TECVAUXSTR111 tecvauxstr111
#  define TECFACE111    tecface111
#  define TECPOLY111    tecpoly111

#  define TECFOREIGN110 tecforeign110
#  define TECINI110     tecini110
#  define TECZNE110     teczne110
#  define TECDAT110     tecdat110
#  define TECNOD110     tecnod110
#  define TECGEO110     tecgeo110
#  define TECTXT110     tectxt110
#  define TECLAB110     teclab110
#  define TECFIL110     tecfil110
#  define TECEND110     tecend110
#  define TECUSR110     tecusr110
#  define TECAUXSTR110  tecauxstr110
#  define TECZAUXSTR110 teczauxstr110
#  define TECVAUXSTR110 tecvauxstr110
#  define TECFACE110    tecface110

#  define TECFOREIGN100 tecforeign100
#  define TECINI100     tecini100
#  define TECZNE100     teczne100
#  define TECDAT100     tecdat100
#  define TECNOD100     tecnod100
#  define TECGEO100     tecgeo100
#  define TECTXT100     tectxt100
#  define TECLAB100     teclab100
#  define TECFIL100     tecfil100
#  define TECEND100     tecend100
#  define TECUSR100     tecusr100
#  define TECAUXSTR100  tecauxstr100
#  define TECZAUXSTR100 teczauxstr100
#  define TECVAUXSTR100 tecvauxstr100
#  define TECFACE100    tecface100

#  define TECINI  tecini
#  define TECZNE  teczne
#  define TECDAT  tecdat
#  define TECNOD  tecnod
#  define TECGEO  tecgeo
#  define TECTXT  tectxt
#  define TECLAB  teclab
#  define TECFIL  tecfil
#  define TECEND  tecend
#  define TECUSR  tecusr
#endif


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
#define INTEGER4  int
#define INTEGER2  short
#endif

#if defined _WIN32
#if !defined MSWIN
#define MSWIN /* MSWIN */
#endif
#endif /* _WIN32 */

#if !defined (EXTERNC)
# if defined (__cplusplus)
#  define EXTERNC extern "C"
# else
#  define EXTERNC
# endif /* __cplusplus */
#endif /* EXTERN_C */

#if !defined (STDCALL)
# if defined MSWIN
#  define STDCALL __stdcall
# else /* !MSWIN */
#  define STDCALL
# endif /* MSWIN */
#endif /* STDCALL */

#if !defined (DLLEXPORT)
# if defined (MSWIN)
#  define DLLEXPORT _declspec (dllexport)
# else
#  define DLLEXPORT
# endif /* MSWIN */
#endif /* DLLEXPORT */

#if !defined (DLLIMPORT)
# if defined (MSWIN)
#  define DLLIMPORT _declspec (dllimport)
# else
#  define DLLIMPORT
# endif /* MSWIN */
#endif /* DLLIMPORT */


#if defined (TECPLOTKERNEL)
/* CORE SOURCE CODE REMOVED */
#else
    #if defined (MAKEARCHIVE)
        #define LIBCALL STDCALL
        #define LIBFUNCTION EXTERNC DLLEXPORT
    #else /* !TECPLOTKERNAL && !MAKEARCHIVE */
        #define LIBCALL STDCALL
        #define LIBFUNCTION EXTERNC DLLIMPORT
    #endif
#endif

/*
 *  V11.3 tecio functions
 */

LIBFUNCTION void LIBCALL TECFOREIGN112(INTEGER4 *OutputForeignByteOrder);

LIBFUNCTION INTEGER4 LIBCALL TECINI112(char     *Title,
                                       char     *Variables,
                                       char     *FName,
                                       char     *ScratchDir,
                                       INTEGER4 *FileType,
                                       INTEGER4 *Debug,
                                       INTEGER4 *VIsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECZNE112(char     *ZoneTitle,
                                       INTEGER4 *ZoneType,
                                       INTEGER4 *IMxOrNumPts,
                                       INTEGER4 *JMxOrNumElements,
                                       INTEGER4 *KMxOrNumFaces,
                                       INTEGER4 *ICellMx,
                                       INTEGER4 *JCellMx,
                                       INTEGER4 *KCellMx,
                                       double   *SolutionTime,
                                       INTEGER4 *StrandID,
                                       INTEGER4 *ParentZone,
                                       INTEGER4 *IsBlock,
                                       INTEGER4 *NumFaceConnections,
                                       INTEGER4 *FaceNeighborMode,
                                       INTEGER4 *TotalNumFaceNodes,
                                       INTEGER4 *NumConnectedBoundaryFaces,
                                       INTEGER4 *TotalNumBoundaryConnections,
                                       INTEGER4 *PassiveVarList,
                                       INTEGER4 *ValueLocation,
                                       INTEGER4 *ShareVarFromZone,
                                       INTEGER4 *ShareConnectivityFromZone);

LIBFUNCTION INTEGER4 LIBCALL TECDAT112(INTEGER4  *N,
                                       void      *FieldData,
                                       INTEGER4  *IsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECNOD112(INTEGER4 *NData);

LIBFUNCTION INTEGER4 LIBCALL TECNODE112(INTEGER4 *N,
                                        INTEGER4 *NData);

LIBFUNCTION INTEGER4 LIBCALL TECEND112(void);

LIBFUNCTION INTEGER4 LIBCALL TECLAB112(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECUSR112(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECGEO112(double    *XPos,
                                       double    *YPos,
                                       double    *ZPos,
                                       INTEGER4  *PosCoordMode,
                                       INTEGER4  *AttachToZone,
                                       INTEGER4  *Zone,
                                       INTEGER4  *Color,
                                       INTEGER4  *FillColor,
                                       INTEGER4  *IsFilled,
                                       INTEGER4  *GeomType,
                                       INTEGER4  *LinePattern,
                                       double    *PatternLength,
                                       double    *LineThickness,
                                       INTEGER4  *NumEllipsePts,
                                       INTEGER4  *ArrowheadStyle,
                                       INTEGER4  *ArrowheadAttachment,
                                       double    *ArrowheadSize,
                                       double    *ArrowheadAngle,
                                       INTEGER4  *Scope,
                                       INTEGER4  *Clipping,
                                       INTEGER4  *NumSegments,
                                       INTEGER4  *NumSegPts,
                                       float     *XGeomData,
                                       float     *YGeomData,
                                       float     *ZGeomData,
                                       char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECTXT112(double    *XOrThetaPos,
                                       double    *YOrRPos,
                                       double    *ZOrUnusedPos,
                                       INTEGER4  *PosCoordMode,
                                       INTEGER4  *AttachToZone,
                                       INTEGER4  *Zone,
                                       INTEGER4  *BFont,
                                       INTEGER4  *FontHeightUnits,
                                       double    *FontHeight,
                                       INTEGER4  *BoxType,
                                       double    *BoxMargin,
                                       double    *BoxLineThickness,
                                       INTEGER4  *BoxColor,
                                       INTEGER4  *BoxFillColor,
                                       double    *Angle,
                                       INTEGER4  *Anchor,
                                       double    *LineSpacing,
                                       INTEGER4  *TextColor,
                                       INTEGER4  *Scope,
                                       INTEGER4  *Clipping,
                                       char      *String,
                                       char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECFIL112(INTEGER4 *F);

LIBFUNCTION INTEGER4 LIBCALL TECAUXSTR112(char *Name,
                                          char *Value);

LIBFUNCTION INTEGER4 LIBCALL TECZAUXSTR112(char *Name,
                                           char *Value);

LIBFUNCTION INTEGER4 LIBCALL TECVAUXSTR112(INTEGER4 *Var,
                                           char     *Name,
                                           char     *Value);

LIBFUNCTION INTEGER4 LIBCALL TECFACE112(INTEGER4 *FaceConnections);

LIBFUNCTION INTEGER4 LIBCALL TECPOLY112(INTEGER4 *FaceNodeCounts,
                                        INTEGER4 *FaceNodes,
                                        INTEGER4 *FaceLeftElems,
                                        INTEGER4 *FaceRightElems,
                                        INTEGER4 *FaceBndryConnectionCounts,
                                        INTEGER4 *FaceBndryConnectionElems,
                                        INTEGER4 *FaceBndryConnectionZones);

/*
 *  V11.1 tecio functions   TODO (JN): Tecplot's version is still in flux so the .1 may change
 */

LIBFUNCTION void LIBCALL TECFOREIGN111(INTEGER4 *OutputForeignByteOrder);

LIBFUNCTION INTEGER4 LIBCALL TECINI111(char     *Title,
                                       char     *Variables,
                                       char     *FName,
                                       char     *ScratchDir,
                                       INTEGER4 *FileType,
                                       INTEGER4 *Debug,
                                       INTEGER4 *VIsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECZNE111(char     *ZoneTitle,
                                       INTEGER4 *ZoneType,
                                       INTEGER4 *IMxOrNumPts,
                                       INTEGER4 *JMxOrNumElements,
                                       INTEGER4 *KMxOrNumFaces,
                                       INTEGER4 *ICellMx,
                                       INTEGER4 *JCellMx,
                                       INTEGER4 *KCellMx,
                                       double   *SolutionTime,
                                       INTEGER4 *StrandID,
                                       INTEGER4 *ParentZone,
                                       INTEGER4 *IsBlock,
                                       INTEGER4 *NumFaceConnections,
                                       INTEGER4 *FaceNeighborMode,
                                       INTEGER4 *TotalNumFaceNodes,
                                       INTEGER4 *NumConnectedBoundaryFaces,
                                       INTEGER4 *TotalNumBoundaryConnections,
                                       INTEGER4 *PassiveVarList,
                                       INTEGER4 *ValueLocation,
                                       INTEGER4 *ShareVarFromZone,
                                       INTEGER4 *ShareConnectivityFromZone);

LIBFUNCTION INTEGER4 LIBCALL TECDAT111(INTEGER4  *N,
                                       void      *FieldData,
                                       INTEGER4  *IsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECNOD111(INTEGER4 *NData);

LIBFUNCTION INTEGER4 LIBCALL TECEND111(void);

LIBFUNCTION INTEGER4 LIBCALL TECLAB111(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECUSR111(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECGEO111(double    *XPos,
                                       double    *YPos,
                                       double    *ZPos,
                                       INTEGER4  *PosCoordMode,
                                       INTEGER4  *AttachToZone,
                                       INTEGER4  *Zone,
                                       INTEGER4  *Color,
                                       INTEGER4  *FillColor,
                                       INTEGER4  *IsFilled,
                                       INTEGER4  *GeomType,
                                       INTEGER4  *LinePattern,
                                       double    *PatternLength,
                                       double    *LineThickness,
                                       INTEGER4  *NumEllipsePts,
                                       INTEGER4  *ArrowheadStyle,
                                       INTEGER4  *ArrowheadAttachment,
                                       double    *ArrowheadSize,
                                       double    *ArrowheadAngle,
                                       INTEGER4  *Scope,
                                       INTEGER4  *Clipping,
                                       INTEGER4  *NumSegments,
                                       INTEGER4  *NumSegPts,
                                       float     *XGeomData,
                                       float     *YGeomData,
                                       float     *ZGeomData,
                                       char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECTXT111(double    *XOrThetaPos,
                                       double    *YOrRPos,
                                       double    *ZOrUnusedPos,
                                       INTEGER4  *PosCoordMode,
                                       INTEGER4  *AttachToZone,
                                       INTEGER4  *Zone,
                                       INTEGER4  *BFont,
                                       INTEGER4  *FontHeightUnits,
                                       double    *FontHeight,
                                       INTEGER4  *BoxType,
                                       double    *BoxMargin,
                                       double    *BoxLineThickness,
                                       INTEGER4  *BoxColor,
                                       INTEGER4  *BoxFillColor,
                                       double    *Angle,
                                       INTEGER4  *Anchor,
                                       double    *LineSpacing,
                                       INTEGER4  *TextColor,
                                       INTEGER4  *Scope,
                                       INTEGER4  *Clipping,
                                       char      *String,
                                       char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECFIL111(INTEGER4 *F);

LIBFUNCTION INTEGER4 LIBCALL TECAUXSTR111(char *Name,
                                          char *Value);

LIBFUNCTION INTEGER4 LIBCALL TECZAUXSTR111(char *Name,
                                           char *Value);

LIBFUNCTION INTEGER4 LIBCALL TECVAUXSTR111(INTEGER4 *Var,
                                           char     *Name,
                                           char     *Value);

LIBFUNCTION INTEGER4 LIBCALL TECFACE111(INTEGER4 *FaceConnections);

LIBFUNCTION INTEGER4 LIBCALL TECPOLY111(INTEGER4 *FaceNodeCounts,
                                        INTEGER4 *FaceNodes,
                                        INTEGER4 *FaceLeftElems,
                                        INTEGER4 *FaceRightElems,
                                        INTEGER4 *FaceBndryConnectionCounts,
                                        INTEGER4 *FaceBndryConnectionElems,
                                        INTEGER2 *FaceBndryConnectionZones);


/*
 * V11 tecio functions
 */

LIBFUNCTION void LIBCALL TECFOREIGN110(INTEGER4 *OutputForeignByteOrder);

LIBFUNCTION INTEGER4 LIBCALL TECINI110(char     *Title,
                                       char     *Variables,
                                       char     *FName,
                                       char     *ScratchDir,
                                       INTEGER4 *Debug,
                                       INTEGER4 *VIsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECZNE110(char     *ZoneTitle,
                                       INTEGER4 *ZoneType,
                                       INTEGER4 *IMxOrNumPts,
                                       INTEGER4 *JMxOrNumElements,
                                       INTEGER4 *KMxOrNumFaces,
                                       INTEGER4 *ICellMx,
                                       INTEGER4 *JCellMx,
                                       INTEGER4 *KCellMx,
                                       double   *SolutionTime,
                                       INTEGER4 *StrandID,
                                       INTEGER4 *ParentZone,
                                       INTEGER4 *IsBlock,
                                       INTEGER4 *NumFaceConnections,
                                       INTEGER4 *FaceNeighborMode,
                                       INTEGER4 *PassiveVarList,
                                       INTEGER4 *ValueLocation,
                                       INTEGER4 *ShareVarFromZone,
                                       INTEGER4 *ShareConnectivityFromZone);

LIBFUNCTION INTEGER4 LIBCALL TECDAT110(INTEGER4  *N,
                                       void      *FieldData,
                                       INTEGER4  *IsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECNOD110(INTEGER4 *NData);

LIBFUNCTION INTEGER4 LIBCALL TECEND110(void);

LIBFUNCTION INTEGER4 LIBCALL TECLAB110(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECUSR110(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECGEO110(double    *XPos,
                                       double    *YPos,
                                       double    *ZPos,
                                       INTEGER4  *PosCoordMode,
                                       INTEGER4  *AttachToZone,
                                       INTEGER4  *Zone,
                                       INTEGER4  *Color,
                                       INTEGER4  *FillColor,
                                       INTEGER4  *IsFilled,
                                       INTEGER4  *GeomType,
                                       INTEGER4  *LinePattern,
                                       double    *PatternLength,
                                       double    *LineThickness,
                                       INTEGER4  *NumEllipsePts,
                                       INTEGER4  *ArrowheadStyle,
                                       INTEGER4  *ArrowheadAttachment,
                                       double    *ArrowheadSize,
                                       double    *ArrowheadAngle,
                                       INTEGER4  *Scope,
                                       INTEGER4  *Clipping,
                                       INTEGER4  *NumSegments,
                                       INTEGER4  *NumSegPts,
                                       float     *XGeomData,
                                       float     *YGeomData,
                                       float     *ZGeomData,
                                       char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECTXT110(double    *XOrThetaPos,
                                       double    *YOrRPos,
                                       double    *ZOrUnusedPos,
                                       INTEGER4  *PosCoordMode,
                                       INTEGER4  *AttachToZone,
                                       INTEGER4  *Zone,
                                       INTEGER4  *BFont,
                                       INTEGER4  *FontHeightUnits,
                                       double    *FontHeight,
                                       INTEGER4  *BoxType,
                                       double    *BoxMargin,
                                       double    *BoxLineThickness,
                                       INTEGER4  *BoxColor,
                                       INTEGER4  *BoxFillColor,
                                       double    *Angle,
                                       INTEGER4  *Anchor,
                                       double    *LineSpacing,
                                       INTEGER4  *TextColor,
                                       INTEGER4  *Scope,
                                       INTEGER4  *Clipping,
                                       char      *String,
                                       char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECFIL110(INTEGER4 *F);

LIBFUNCTION INTEGER4 LIBCALL TECAUXSTR110(char *Name,
                                          char *Value);

LIBFUNCTION INTEGER4 LIBCALL TECZAUXSTR110(char *Name,
                                           char *Value);

LIBFUNCTION INTEGER4 LIBCALL TECVAUXSTR110(INTEGER4 *Var,
                                           char     *Name,
                                           char     *Value);

LIBFUNCTION INTEGER4 LIBCALL TECFACE110(INTEGER4 *FaceConnections);


/*
 * V10 tecio functions kept for backward compatability.
 */

LIBFUNCTION void LIBCALL TECFOREIGN100(INTEGER4 *OutputForeignByteOrder);

LIBFUNCTION INTEGER4 LIBCALL TECINI100(char     *Title,
                                       char     *Variables,
                                       char     *FName,
                                       char     *ScratchDir,
                                       INTEGER4 *Debug,
                                       INTEGER4 *VIsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECZNE100(char     *ZoneTitle,
                                       INTEGER4 *ZoneType,
                                       INTEGER4 *IMxOrNumPts,
                                       INTEGER4 *JMxOrNumElements,
                                       INTEGER4 *KMxOrNumFaces,
                                       INTEGER4 *ICellMx,
                                       INTEGER4 *JCellMx,
                                       INTEGER4 *KCellMx,
                                       INTEGER4 *IsBlock,
                                       INTEGER4 *NumFaceConnections,
                                       INTEGER4 *FaceNeighborMode,
                                       INTEGER4 *ValueLocation,
                                       INTEGER4 *ShareVarFromZone,
                                       INTEGER4 *ShareConnectivityFromZone);

LIBFUNCTION INTEGER4 LIBCALL TECDAT100(INTEGER4  *N,
                                       void      *FieldData,
                                       INTEGER4  *IsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECNOD100(INTEGER4 *NData);

LIBFUNCTION INTEGER4 LIBCALL TECEND100(void);

LIBFUNCTION INTEGER4 LIBCALL TECLAB100(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECUSR100(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECGEO100(double    *XPos,
                                       double    *YPos,
                                       double    *ZPos,
                                       INTEGER4  *PosCoordMode,
                                       INTEGER4  *AttachToZone,
                                       INTEGER4  *Zone,
                                       INTEGER4  *Color,
                                       INTEGER4  *FillColor,
                                       INTEGER4  *IsFilled,
                                       INTEGER4  *GeomType,
                                       INTEGER4  *LinePattern,
                                       double    *PatternLength,
                                       double    *LineThickness,
                                       INTEGER4  *NumEllipsePts,
                                       INTEGER4  *ArrowheadStyle,
                                       INTEGER4  *ArrowheadAttachment,
                                       double    *ArrowheadSize,
                                       double    *ArrowheadAngle,
                                       INTEGER4  *Scope,
                                       INTEGER4  *Clipping,
                                       INTEGER4  *NumSegments,
                                       INTEGER4  *NumSegPts,
                                       float     *XGeomData,
                                       float     *YGeomData,
                                       float     *ZGeomData,
                                       char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECTXT100(double    *XOrThetaPos,
                                       double    *YOrRPos,
                                       double    *ZOrUnusedPos,
                                       INTEGER4  *PosCoordMode,
                                       INTEGER4  *AttachToZone,
                                       INTEGER4  *Zone,
                                       INTEGER4  *BFont,
                                       INTEGER4  *FontHeightUnits,
                                       double    *FontHeight,
                                       INTEGER4  *BoxType,
                                       double    *BoxMargin,
                                       double    *BoxLineThickness,
                                       INTEGER4  *BoxColor,
                                       INTEGER4  *BoxFillColor,
                                       double    *Angle,
                                       INTEGER4  *Anchor,
                                       double    *LineSpacing,
                                       INTEGER4  *TextColor,
                                       INTEGER4  *Scope,
                                       INTEGER4  *Clipping,
                                       char      *String,
                                       char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECFIL100(INTEGER4 *F);

LIBFUNCTION INTEGER4 LIBCALL TECAUXSTR100(char *Name,
                                          char *Value);

LIBFUNCTION INTEGER4 LIBCALL TECZAUXSTR100(char *Name,
                                           char *Value);

LIBFUNCTION INTEGER4 LIBCALL TECVAUXSTR100(INTEGER4 *Var,
                                           char     *Name,
                                           char     *Value);

LIBFUNCTION INTEGER4 LIBCALL TECFACE100(INTEGER4 *FaceConnections);

/* Old V9 functions retained for backward compatibility */

LIBFUNCTION INTEGER4 LIBCALL TECINI(char     *Title,
                                    char     *Variables,
                                    char     *FName,
                                    char     *ScratchDir,
                                    INTEGER4 *Debug,
                                    INTEGER4 *VIsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECZNE(char     *ZoneTitle,
                                    INTEGER4 *IMx,
                                    INTEGER4 *JMx,
                                    INTEGER4 *KMx,
                                    char     *ZFormat,
                                    char     *DupList);

LIBFUNCTION INTEGER4 LIBCALL TECDAT(INTEGER4  *N,
                                    void      *FieldData,
                                    INTEGER4  *IsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECNOD(INTEGER4 *NData);

LIBFUNCTION INTEGER4 LIBCALL TECEND(void);

LIBFUNCTION INTEGER4 LIBCALL TECLAB(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECUSR(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECGEO(double    *XPos,
                                    double    *YPos,
                                    double    *ZPos,
                                    INTEGER4  *PosCoordMode,
                                    INTEGER4  *AttachToZone,
                                    INTEGER4  *Zone,
                                    INTEGER4  *Color,
                                    INTEGER4  *FillColor,
                                    INTEGER4  *IsFilled,
                                    INTEGER4  *GeomType,
                                    INTEGER4  *LinePattern,
                                    double    *PatternLength,
                                    double    *LineThickness,
                                    INTEGER4  *NumEllipsePts,
                                    INTEGER4  *ArrowheadStyle,
                                    INTEGER4  *ArrowheadAttachment,
                                    double    *ArrowheadSize,
                                    double    *ArrowheadAngle,
                                    INTEGER4  *Scope,
                                    INTEGER4  *NumSegments,
                                    INTEGER4  *NumSegPts,
                                    float     *XGeomData,
                                    float     *YGeomData,
                                    float     *ZGeomData,
                                    char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECTXT(double    *XPos,
                                    double    *YPos,
                                    INTEGER4  *PosCoordMode,
                                    INTEGER4  *AttachToZone,
                                    INTEGER4  *Zone,
                                    INTEGER4  *BFont,
                                    INTEGER4  *FontHeightUnits,
                                    double    *FontHeight,
                                    INTEGER4  *BoxType,
                                    double    *BoxMargin,
                                    double    *BoxLineThickness,
                                    INTEGER4  *BoxColor,
                                    INTEGER4  *BoxFillColor,
                                    double    *Angle,
                                    INTEGER4  *Anchor,
                                    double    *LineSpacing,
                                    INTEGER4  *TextColor,
                                    INTEGER4  *Scope,
                                    char      *Text,
                                    char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECFIL(INTEGER4 *F);

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#endif /* TECXXX_H_ */
