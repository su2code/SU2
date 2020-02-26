#pragma once

/*
 * TECIO.h: Copyright (C) 1988-2014 Tecplot, Inc.
 */

#if defined TECIOMPI
    #include "mpi.h"
#endif

#include "StandardIntegralTypes.h"

#define INTEGER8  int64_t
#define INTEGER4  int32_t
#define INTEGER2  int16_t

#if defined _WIN32
# if !defined MSWIN
#   define MSWIN
# endif
#endif

#include "tecio_Exports.h"

#if !defined STDCALL
# if defined MSWIN
#   define STDCALL __stdcall
# else
#   define STDCALL
# endif
#endif

#if !defined EXTERNC
# if defined __cplusplus
#  define EXTERNC extern "C"
# else
#  define EXTERNC
# endif
#endif

#if !defined _GLOBAL_H

    enum FieldDataType_e
    {
        FieldDataType_Float = 1,
        FieldDataType_Double,
        FieldDataType_Int32,
        FieldDataType_Int16,
        FieldDataType_Byte,
    };

    enum GeomType_e
    {
        GeomType_LineSegs = 0,
        GeomType_Rectangle = 1,
        GeomType_Square = 2,
        GeomType_Circle = 3,
        GeomType_Ellipse = 4,
        GeomType_LineSegs3D = 5
    };

    enum CoordSys_e
    {
        CoordSys_Grid = 0,
        CoordSys_Frame = 1,
        CoordSys_Grid3D = 6
    };

#endif

int32_t const FILEFORMAT_PLT = 0;
int32_t const FILEFORMAT_SZL = 1;

int32_t const FILETYPE_FULL = 0;
int32_t const FILETYPE_GRID = 1;
int32_t const FILETYPE_SOLUTION = 2;

// Use these for Tecio; they are different from the ZoneType enum for add-ons in GLOBAL.h.
int32_t const ZONETYPE_ORDERED = 0;
int32_t const ZONETYPE_FELINESEG = 1;
int32_t const ZONETYPE_FETRIANGLE = 2;
int32_t const ZONETYPE_FEQUADRILATERAL = 3;
int32_t const ZONETYPE_FETETRAHEDRON = 4;
int32_t const ZONETYPE_FEBRICK = 5;
int32_t const ZONETYPE_FEPOLYGON = 6;
int32_t const ZONETYPE_FEPOLYHEDRON = 7;

/**
 * New TecIO output routines support 64-bit output, all var types except bit, and out-of-order data output.
 * SZL output only (may support plt in the future), and no polyhedral support yet.
 */

EXTERNC tecio_API int32_t tecFileWriterOpen(
    char const* fileName,
    char const* dataSetTitle,
    char const* variableList,
    int32_t     fileFormat,
    int32_t     fileType,
    int32_t     defaultVarType,
    void*       gridFileHandle,
    void**      fileHandle);

EXTERNC tecio_API int32_t tecFileSetDiagnosticsLevel(
    void*   fileHandle,
    int32_t level);

#if defined TECIOMPI
EXTERNC tecio_API int32_t tecMPIInitialize(
    void*    fileHandle,
    MPI_Comm communicator,
    int32_t  mainRank);
#endif

EXTERNC tecio_API int32_t tecZoneCreateIJK(
    void*          fileHandle,
    char const*    zoneTitle,
    int64_t        imax,
    int64_t        jmax,
    int64_t        kmax,
    int32_t const* varTypes,
    int32_t const* shareVarFromZone,
    int32_t const* valueLocations,
    int32_t const* passiveVarList,
    int32_t        shareFaceNeighborsFromZone,
    int64_t        numFaceConnections,
    int32_t        faceNeighborMode,
    int32_t*       zone);

EXTERNC tecio_API int32_t tecZoneCreateFE(
    void*          fileHandle,
    char const*    zoneTitle,
    int32_t        zoneType,
    int64_t        numNodes,
    int64_t        numCells,
    int32_t const* varTypes,
    int32_t const* shareVarFromZone,
    int32_t const* valueLocations,
    int32_t const* passiveVarList,
    int32_t        shareConnectivityFromZone,
    int64_t        numFaceConnections,
    int32_t        faceNeighborMode,
    int32_t*       zone);

EXTERNC tecio_API int32_t tecZoneCreatePoly(
    void*          fileHandle,
    char const*    zoneTitle,
    int32_t        zoneType,
    int64_t        numNodes,
    int64_t        numFaces,
    int64_t        numCells,
    int64_t        totalNumFaceNodes,
    int32_t const* varTypes,
    int32_t const* shareVarFromZone,
    int32_t const* valueLocations,
    int32_t const* passiveVarList,
    int32_t        shareConnectivityFromZone,
    int64_t        numConnectedBoundaryFaces,
    int64_t        totalNumBoundaryConnections,
    int32_t*       zone);

EXTERNC tecio_API int32_t tecZoneSetUnsteadyOptions(
    void*   fileHandle,
    int32_t zone,
    double  solutionTime,
    int32_t strandID);

EXTERNC tecio_API int32_t tecZoneSetParentZone(
    void*   fileHandle,
    int32_t zone,
    int32_t parentZone);

#if defined TECIOMPI
EXTERNC tecio_API int32_t tecZoneMapPartitionsToMPIRanks(
    void*          fileHandle,
    int32_t&       zone, // Number may change due to communication with the main output rank
    int32_t        numPartitions,
    int32_t const* mpiRanksForPartitions);
#endif

EXTERNC tecio_API int32_t tecFEPartitionCreate32(
    void*          fileHandle,
    int32_t        zone,
    int32_t        partition,
    int64_t        numNodes,
    int64_t        numCells,
    int64_t        numGhostNodes,
    int32_t const* ghostNodes,
    int32_t const* neighborPartitions,
    int32_t const* neighborPartitionNodes,
    int64_t        numGhostCells,
    int32_t const* ghostCells);

EXTERNC tecio_API int32_t tecFEPartitionCreate64(
    void*          fileHandle,
    int32_t        zone,
    int32_t        partition,
    int64_t        numNodes,
    int64_t        numCells,
    int64_t        numGhostNodes,
    int64_t const* ghostNodes,
    int32_t const* neighborPartitions,
    int64_t const* neighborPartitionNodes,
    int64_t        numGhostCells,
    int64_t const* ghostCells);

EXTERNC tecio_API int32_t tecIJKPartitionCreate(
    void*   fileHandle,
    int32_t zone,
    int32_t partition,
    int64_t imin,
    int64_t jmin,
    int64_t kmin,
    int64_t imax,
    int64_t jmax,
    int64_t kmax);

EXTERNC tecio_API int32_t tecZoneVarWriteDoubleValues(
    void*         fileHandle,
    int32_t       zone,
    int32_t       var,
    int32_t       partition,
    int64_t       count,
    double const* values);

EXTERNC tecio_API int32_t tecZoneVarWriteFloatValues(
    void*        fileHandle,
    int32_t      zone,
    int32_t      var,
    int32_t      partition,
    int64_t      count,
    float const* values);

EXTERNC tecio_API int32_t tecZoneVarWriteInt32Values(
    void*          fileHandle,
    int32_t        zone,
    int32_t        var,
    int32_t        partition,
    int64_t        count,
    int32_t const* values);

EXTERNC tecio_API int32_t tecZoneVarWriteInt16Values(
    void*          fileHandle,
    int32_t        zone,
    int32_t        var,
    int32_t        partition,
    int64_t        count,
    int16_t const* values);

EXTERNC tecio_API int32_t tecZoneVarWriteUInt8Values(
    void*          fileHandle,
    int32_t        zone,
    int32_t        var,
    int32_t        partition,
    int64_t        count,
    uint8_t const* values);

EXTERNC tecio_API int32_t tecZoneNodeMapWrite32(
    void*          fileHandle,
    int32_t        zone,
    int32_t        partition,
    int32_t        nodesAreOneBased,
    int64_t        count,
    int32_t const* nodes);

EXTERNC tecio_API int32_t tecZoneNodeMapWrite64(
    void*          fileHandle,
    int32_t        zone,
    int32_t        partition,
    int32_t        nodeAreOneBased,
    int64_t        count,
    int64_t const* nodes);

EXTERNC tecio_API int32_t tecZoneFaceNbrWriteConnections32(
    void*          fileHandle,
    int32_t        zone,
    int32_t const* faceNeighbors);

EXTERNC tecio_API int32_t tecZoneFaceNbrWriteConnections64(
    void*          fileHandle,
    int32_t        zone,
    int64_t const* faceNeighbors);

EXTERNC tecio_API int32_t tecZoneWritePolyFaces32(
    void*          fileHandle,
    int32_t        zone,
    int32_t        partition,
    int32_t        numFaces,
    int32_t const* faceNodeCounts,
    int32_t const* faceNodes,
    int32_t const* faceLeftElems,
    int32_t const* faceRightElems,
    int32_t        isOneBased);

EXTERNC tecio_API int32_t tecZoneWritePolyFaces64(
    void*          fileHandle,
    int32_t        zone,
    int32_t        partition,
    int64_t        numFaces,
    int32_t const* faceNodeCounts,
    int64_t const* faceNodes,
    int64_t const* faceLeftElems,
    int64_t const* faceRightElems,
    int32_t        isOneBased);

EXTERNC tecio_API int32_t tecZoneWritePolyBoundaryConnections32(
    void*          fileHandle,
    int32_t        zone,
    int32_t        partition,
    int32_t        numBoundaryFaces,
    int32_t const* faceBoundaryConnectionCounts,
    int32_t const* faceBoundaryConnectionElems,
    int32_t const* faceBoundaryConnectionZones,
    int32_t        isOneBased);

EXTERNC tecio_API int32_t tecZoneWritePolyBoundaryConnections64(
    void*          fileHandle,
    int32_t        zone,
    int32_t        partition,
    int64_t        numBoundaryFaces,
    int32_t const* faceBoundaryConnectionCounts,
    int64_t const* faceBoundaryConnectionElems,
    int32_t const* faceBoundaryConnectionZones,
    int32_t        isOneBased);

EXTERNC tecio_API int32_t tecDataSetAddAuxData(
    void*       fileHandle,
    char const* name,
    char const* value);

EXTERNC tecio_API int32_t tecVarAddAuxData(
    void*       fileHandle,
    int32_t     var,
    char const* name,
    char const* value);

EXTERNC tecio_API int32_t tecZoneAddAuxData(
    void*       fileHandle,
    int32_t     zone,
    char const* name,
    char const* value);

EXTERNC tecio_API int32_t tecGeom2DLineSegmentsBegin(
    void*          fileHandle,
    double         xOrigin,
    double         yOrigin,
    int32_t        numPoints,
    double const*  relativeX,
    double const*  relativeY,
    int32_t        posCoordMode);

EXTERNC tecio_API int32_t tecGeom2DMultiLineSegmentsBegin(
    void*          fileHandle,
    double         xOrigin,
    double         yOrigin,
    int32_t        numSegments,
    int32_t const* numSegmentPoints,
    double const*  relativeX,
    double const*  relativeY,
    int32_t        posCoordMode);

EXTERNC tecio_API int32_t tecGeom3DLineSegmentsBegin(
    void*          fileHandle,
    double         xOrigin,
    double         yOrigin,
    double         zOrigin,
    int32_t        numPoints,
    double const*  relativeX,
    double const*  relativeY,
    double const*  relativeZ);

EXTERNC tecio_API int32_t tecGeom3DMultiLineSegmentsBegin(
    void*          fileHandle,
    double         xOrigin,
    double         yOrigin,
    double         zOrigin,
    int32_t        numSegments,
    int32_t const* numSegmentPoints,
    double const*  relativeX,
    double const*  relativeY,
    double const*  relativeZ);

EXTERNC tecio_API int32_t tecGeomCircleBegin(
    void*   fileHandle,
    double  xCenter,
    double  yCenter,
    double  radius,
    int32_t posCoordMode);

EXTERNC tecio_API int32_t tecGeomEllipseBegin(
    void*   fileHandle,
    double  xCenter,
    double  yCenter,
    double  width,
    double  height,
    int32_t posCoordMode);

EXTERNC tecio_API int32_t tecGeomRectangleBegin(
    void*   fileHandle,
    double  xMin,
    double  yMin,
    double  xMax,
    double  yMax,
    int32_t posCoordMode);

EXTERNC tecio_API int32_t tecGeomSquareBegin(
    void*   fileHandle,
    double  xMin,
    double  yMin,
    double  size,
    int32_t posCoordMode);

EXTERNC tecio_API int32_t tecGeomArrowheadSetInfo(
    void*   fileHandle,
    double  angle,
    int32_t attachment,
    double  size,
    int32_t style);

EXTERNC tecio_API int32_t tecGeomEllipseSetNumPoints(
    void*   fileHandle,
    int32_t numEllipsePoints);

EXTERNC tecio_API int32_t tecGeomSetClipping(
    void*   fileHandle,
    int32_t clipping);

EXTERNC tecio_API int32_t tecGeomSetLineInfo(
    void*   fileHandle,
    int32_t linePattern,
    double  patternLength,
    double  thickness,
    int32_t color);

EXTERNC tecio_API int32_t tecGeomSetMacroFunctionCmd(
    void*       fileHandle,
    char const* macroFunctionCmd);

EXTERNC tecio_API int32_t tecGeomSetScope(
    void*   fileHandle,
    int32_t scope);

EXTERNC tecio_API int32_t tecGeomAttachToZone(
    void*   fileHandle,
    int32_t zone);

EXTERNC tecio_API int32_t tecGeomFill(
    void*   fileHandle,
    int32_t fillColor);

EXTERNC tecio_API int32_t tecGeomEnd(
    void* fileHandle);

EXTERNC tecio_API int32_t tecCustomLabelsAddSet(
    void*       fileHandle,
    char const* labels);

EXTERNC tecio_API int32_t tecText2DBegin(
    void*       fileHandle,
    char const* string,
    double      x,
    double      y,
    int32_t     posCoordMode,
    double      height,
    int32_t     sizeUnits);

EXTERNC tecio_API int32_t tecText3DBegin(
    void*       fileHandle,
    char const* string,
    double      x,
    double      y,
    double      z,
    double      height,
    int32_t     sizeUnits);

EXTERNC tecio_API int32_t tecTextAttachToZone(
    void*   fileHandle,
    int32_t zone);

EXTERNC tecio_API int32_t tecTextBoxSetInfo(
    void*   fileHandle,
    int32_t boxType,
    int32_t lineColor,
    int32_t fillColor,
    double  lineThickness,
    double  margin);

EXTERNC tecio_API int32_t tecTextSetAnchor(
    void*   fileHandle,
    int32_t anchor);

EXTERNC tecio_API int32_t tecTextSetAngle(
    void*  fileHandle,
    double angle);

EXTERNC tecio_API int32_t tecTextSetClipping(
    void*   fileHandle,
    int32_t clipping);

EXTERNC tecio_API int32_t tecTextSetColor(
    void*   fileHandle,
    int32_t color);

EXTERNC tecio_API int32_t tecTextSetTypeface(
    void*       fileHandle,
    char const* family,
    int32_t     isBold,
    int32_t     isItalic);

EXTERNC tecio_API int32_t tecTextSetLineSpacing(
    void*  fileHandle,
    double lineSpacing);

EXTERNC tecio_API int32_t tecTextSetMacroFunctionCmd(
    void*       fileHandle,
    char const* macroFunctionCmd);

EXTERNC tecio_API int32_t tecTextSetScope(
    void*   fileHandle,
    int32_t scope);

EXTERNC tecio_API int32_t tecTextEnd(
    void* fileHandle);

EXTERNC tecio_API int32_t tecUserRecAdd(
    void*       fileHandle,
    char const* userRec);

EXTERNC tecio_API int32_t tecFileWriterFlush(
    void*          fileHandle,
    int32_t        numZonesToRetain,
    int32_t const* zonesToRetain);

EXTERNC tecio_API int32_t tecFileWriterClose(
    void** fileHandle);

/**
 * SZL file reading routines
 */

EXTERNC tecio_API int32_t tecCustomLabelsGetNumSets(
    void*    fileHandle,
    int32_t* numSets);

EXTERNC tecio_API int32_t
tecCustomLabelsGetSet(
    void*   fileHandle,
    int32_t whichSet,
    char**  labelSet);

EXTERNC tecio_API int32_t tecDataSetAuxDataGetItem(
    void*   fileHandle,
    int32_t whichItem,
    char**  name,
    char**  value);

EXTERNC tecio_API int32_t tecDataSetAuxDataGetNumItems(
    void*    fileHandle,
    int32_t* numItems);

EXTERNC tecio_API int32_t tecDataSetGetNumVars(
    void*    fileHandle,
    int32_t* numVars);

EXTERNC tecio_API int32_t tecDataSetGetNumZones(
    void*    fileHandle,
    int32_t* numZones);

EXTERNC tecio_API int32_t tecDataSetGetTitle(
    void*  fileHandle,
    char** title);

EXTERNC tecio_API int32_t tecFileGetType(
    void*    fileHandle,
    int32_t* fileType);

EXTERNC tecio_API int32_t tecFileReaderClose(
    void** fileHandle);

EXTERNC tecio_API int32_t tecFileReaderOpen(
    char const* fileName,
    void**      fileHandle);

EXTERNC tecio_API int32_t tecGeomArrowheadGetAngle(
    void*   fileHandle,
    int32_t geom,
    double* angle);

EXTERNC tecio_API int32_t tecGeomArrowheadGetAttach(
    void* fileHandle,
    int32_t geom,
    int32_t* attachment);

EXTERNC tecio_API int32_t tecGeomArrowheadGetSize(
    void*   fileHandle,
    int32_t geom,
    double* arrowheadSize);

EXTERNC tecio_API int32_t tecGeomArrowheadGetStyle(
    void*    fileHandle,
    int32_t  geom,
    int32_t* arrowheadStyle);

EXTERNC tecio_API int32_t tecGeomCircleGetRadius(
    void*   fileHandle,
    int32_t geom,
    double* radius);

EXTERNC tecio_API int32_t tecGeomEllipseGetNumPoints(
    void*    fileHandle,
    int32_t  geom,
    int32_t* numEllipsePoints);

EXTERNC tecio_API int32_t tecGeomEllipseGetSize(
    void*   fileHandle,
    int32_t geom,
    double* horizontalAxis,
    double* verticalAxis);

EXTERNC tecio_API int32_t tecGeomGetAnchorPos(
    void*   fileHandle,
    int32_t geom,
    double* x,
    double* y,
    double* z);

EXTERNC tecio_API int32_t tecGeomGetClipping(
    void*    fileHandle,
    int32_t  geom,
    int32_t* clipping);

EXTERNC tecio_API int32_t tecGeomGetColor(
    void*    fileHandle,
    int32_t  geom,
    int32_t* color);

EXTERNC tecio_API int32_t tecGeomGetCoordMode(
    void*    fileHandle,
    int32_t  geom,
    int32_t* coordMode);

EXTERNC tecio_API int32_t tecGeomGetFillColor(
    void*    fileHandle,
    int32_t  geom,
    int32_t* fillColor);

EXTERNC tecio_API int32_t tecGeomGetLinePattern(
    void*    fileHandle,
    int32_t  geom,
    int32_t* linePattern);

EXTERNC tecio_API int32_t tecGeomGetLineThickness(
    void*   fileHandle,
    int32_t geom,
    double* lineThickness);

EXTERNC tecio_API int32_t tecGeomGetMacroFunctionCmd(
    void*   fileHandle,
    int32_t geom,
    char**  macroFunctionCmd);

EXTERNC tecio_API int32_t tecGeomGetNumGeoms(
    void*    fileHandle,
    int32_t* numGeoms);

EXTERNC tecio_API int32_t tecGeomGetPatternLength(
    void*   fileHandle,
    int32_t geom,
    double* patternLength);

EXTERNC tecio_API int32_t tecGeomGetScope(
    void*    fileHandle,
    int32_t  geom,
    int32_t* scope);

EXTERNC tecio_API int32_t tecGeomGetType(
    void*    fileHandle,
    int32_t  geom,
    int32_t* type);

EXTERNC tecio_API int32_t tecGeomGetZone(
    void*    fileHandle,
    int32_t  geom,
    int32_t* zone);

EXTERNC tecio_API int32_t tecGeomIsAttached(
    void*    fileHandle,
    int32_t  geom,
    int32_t* isAttached);

EXTERNC tecio_API int32_t tecGeomIsFilled(
    void*    fileHandle,
    int32_t  geom,
    int32_t* isFilled);

EXTERNC tecio_API int32_t tecGeomLineGetPoint(
    void*   fileHandle,
    int32_t geom,
    int32_t segment,
    int32_t index,
    double* x,
    double* y,
    double* z);

EXTERNC tecio_API int32_t tecGeomLineGetSegmentCount(
    void*    fileHandle,
    int32_t  geom,
    int32_t* segmentCount);

EXTERNC tecio_API int32_t tecGeomLineSegmentGetPointCount(
    void*    fileHandle,
    int32_t  geom,
    int32_t  segment,
    int32_t* pointCount);

EXTERNC tecio_API int32_t tecGeomRectangleGetSize(
    void*   fileHandle,
    int32_t geom,
    double* width,
    double* height);

EXTERNC tecio_API int32_t tecGeomSquareGetSize(
    void*   fileHandle,
    int32_t geom,
    double* size);

EXTERNC tecio_API void tecStringFree(
    char** string);

EXTERNC tecio_API int32_t tecStringLength(
    char const* string);

EXTERNC tecio_API int32_t tecTextBoxGetColor(
    void*    fileHandle,
    int32_t  text,
    int32_t* boxColor);

EXTERNC tecio_API int32_t tecTextBoxGetFillColor(
    void*    fileHandle,
    int32_t  text,
    int32_t* boxFillColor);

EXTERNC tecio_API int32_t tecTextBoxGetLineThickness(
    void*   fileHandle,
    int32_t text,
    double* boxLineThickness);

EXTERNC tecio_API int32_t tecTextBoxGetMargin(
    void*   fileHandle,
    int32_t text,
    double* boxMargin);

EXTERNC tecio_API int32_t tecTextBoxGetType(
    void*    fileHandle,
    int32_t  text,
    int32_t* boxType);

EXTERNC tecio_API int32_t tecTextGetAnchor(
    void*    fileHandle,
    int32_t  text,
    int32_t* anchor);

EXTERNC tecio_API int32_t tecTextGetAnchorPos(
    void*   fileHandle,
    int32_t text,
    double* x,
    double* y,
    double* z);

EXTERNC tecio_API int32_t tecTextGetAngle(
    void*   fileHandle,
    int32_t text,
    double* angle);

EXTERNC tecio_API int32_t tecTextGetClipping(
    void*    fileHandle,
    int32_t  text,
    int32_t* clipping);

EXTERNC tecio_API int32_t tecTextGetColor(
    void*    fileHandle,
    int32_t  text,
    int32_t* color);

EXTERNC tecio_API int32_t tecTextGetCoordMode(
    void*    fileHandle,
    int32_t  text,
    int32_t* coordMode);

EXTERNC tecio_API int32_t tecTextGetHeight(
    void*   fileHandle,
    int32_t text,
    double* height);

EXTERNC tecio_API int32_t tecTextGetLineSpacing(
    void*   fileHandle,
    int32_t text,
    double* lineSpacing);

EXTERNC tecio_API int32_t tecTextGetMacroFunctionCmd(
    void*   fileHandle,
    int32_t text,
    char**  macroFunctionCmd);

EXTERNC tecio_API int32_t tecTextGetNumTexts(
    void*    fileHandle,
    int32_t* numTexts);

EXTERNC tecio_API int32_t tecTextGetScope(
    void*    fileHandle,
    int32_t  text,
    int32_t* scope);

EXTERNC tecio_API int32_t tecTextGetSizeUnits(
    void*    fileHandle,
    int32_t  text,
    int32_t* sizeUnits);

EXTERNC tecio_API int32_t tecTextGetString(
    void*   fileHandle,
    int32_t text,
    char**  string);

EXTERNC tecio_API int32_t tecTextGetTypeface(
    void*   fileHandle,
    int32_t text,
    char**  typeface);

EXTERNC tecio_API int32_t tecTextGetZone(
    void*    fileHandle,
    int32_t  text,
    int32_t* zone);

EXTERNC tecio_API int32_t tecTextIsAttached(
    void*    fileHandle,
    int32_t  text,
    int32_t* isAttached);

EXTERNC tecio_API int32_t tecTextIsBold(
    void*    fileHandle,
    int32_t  text,
    int32_t* isBold);

EXTERNC tecio_API int32_t tecTextIsItalic(
    void*    fileHandle,
    int32_t  text,
    int32_t* isItalic);

EXTERNC tecio_API int32_t tecVarAuxDataGetItem(
    void*   fileHandle,
    int32_t var,
    int32_t whichItem,
    char**  name,
    char**  value);

EXTERNC tecio_API int32_t tecVarAuxDataGetNumItems(
    void*    fileHandle,
    int32_t  var,
    int32_t* numItems);

EXTERNC tecio_API int32_t tecVarGetName(
    void*   fileHandle,
    int32_t var,
    char**  name);

EXTERNC tecio_API int32_t tecVarIsEnabled(
    void*    fileHandle,
    int32_t  var,
    int32_t* isEnabled);

EXTERNC tecio_API int32_t tecZoneAuxDataGetItem(
    void*   fileHandle,
    int32_t zone,
    int32_t whichItem,
    char**  name,
    char**  value);

EXTERNC tecio_API int32_t tecZoneAuxDataGetNumItems(
    void*    fileHandle,
    int32_t  zone,
    int32_t* numItems);

EXTERNC tecio_API int32_t tecZoneConnectivityGetSharedZone(
    void*    fileHandle,
    int32_t  zone,
    int32_t* sharedZone);

EXTERNC tecio_API int32_t tecZoneFaceNbrGetConnections(
    void*    fileHandle,
    int32_t  zone,
    int32_t* connections);

EXTERNC tecio_API int32_t tecZoneFaceNbrGetConnections64(
    void*    fileHandle,
    int32_t  zone,
    int64_t* connections);

EXTERNC tecio_API int32_t tecZoneFaceNbrGetMode(
    void*    fileHandle,
    int32_t  zone,
    int32_t* mode);

EXTERNC tecio_API int32_t tecZoneFaceNbrGetNumConnections(
    void*    fileHandle,
    int32_t  zone,
    int64_t* numConnections);

EXTERNC tecio_API int32_t tecZoneFaceNbrGetNumValues(
    void*    fileHandle,
    int32_t  zone,
    int64_t* numValues);

EXTERNC tecio_API int32_t tecZoneFaceNbrsAre64Bit(
    void*    fileHandle,
    int32_t  zone,
    int32_t* are64Bit);

EXTERNC tecio_API int32_t tecZoneGetIJK(
    void*    fileHandle,
    int32_t  zone,
    int64_t* iMax,
    int64_t* jMax,
    int64_t* kMax);

EXTERNC tecio_API int32_t tecZoneGetParentZone(
    void*    fileHandle,
    int32_t  zone,
    int32_t* parentZone);

EXTERNC tecio_API int32_t tecZoneGetSolutionTime(
    void*   fileHandle,
    int32_t zone,
    double* solutionTime);

EXTERNC tecio_API int32_t tecZoneGetStrandID(
    void*    fileHandle,
    int32_t  zone,
    int32_t* strandID);

EXTERNC tecio_API int32_t tecZoneGetTitle(
    void*   fileHandle,
    int32_t zone,
    char**  title);

EXTERNC tecio_API int32_t tecZoneGetType(
    void*    fileHandle,
    int32_t  zone,
    int32_t* type);

EXTERNC tecio_API int32_t tecZoneIsEnabled(
    void*   fileHandle,
    int32_t  zone,
    int32_t* isEnabled);

EXTERNC tecio_API int32_t tecZoneNodeMapGet(
    void*    fileHandle,
    int32_t  zone,
    int64_t  startCell,
    int64_t  numCells,
    int32_t* nodeMap);

EXTERNC tecio_API int32_t tecZoneNodeMapGet64(
    void*    fileHandle,
    int32_t  zone,
    int64_t  startCell,
    int64_t  numCells,
    int64_t* nodeMap);

EXTERNC tecio_API int32_t tecZoneNodeMapGetNumValues(
    void*    fileHandle,
    int32_t  zone,
    int64_t  numCells,
    int64_t* numValues);

EXTERNC tecio_API int32_t tecZoneNodeMapIs64Bit(
    void*    fileHandle,
    int32_t  zone,
    int32_t* is64Bit);

EXTERNC tecio_API int32_t tecZonePolyGetBoundaryConnectionCounts(
    void*    fileHandle,
    int32_t  zone,
    int64_t  startConnection,
    int64_t  numConnections,
    int32_t* connectionCounts);

EXTERNC tecio_API int32_t tecZonePolyGetBoundaryConnections(
    void*    fileHandle,
    int32_t  zone,
    int64_t  startConnection,
    int64_t  numConnections,
    int32_t* connectedElements,
    int32_t* connectedZones);

EXTERNC tecio_API int32_t tecZonePolyGetFaceElems(
    void*    fileHandle,
    int32_t  zone,
    int64_t  startFace,
    int64_t  numFaces,
    int32_t* leftElems,
    int32_t* rightElems);

EXTERNC tecio_API int32_t tecZonePolyGetFaceNodeCounts(
    void*    fileHandle,
    int32_t  zone,
    int64_t  startFace,
    int64_t  numFaces,
    int32_t* nodeCounts);

EXTERNC tecio_API int32_t tecZonePolyGetFaceNodes(
    void*    fileHandle,
    int32_t  zone,
    int64_t  startFace,
    int64_t  numFaces,
    int32_t* faceNodes);

EXTERNC tecio_API int32_t tecZonePolyGetNumConnectedBoundaryFaces(
    void*    fileHandle,
    int32_t  zone,
    int64_t* numFaces);

EXTERNC tecio_API int32_t tecZonePolyGetTotalNumFaceNodes(
    void*    fileHandle,
    int32_t  zone,
    int64_t* numNodes);

EXTERNC tecio_API int32_t tecZonePolyGetTotalNumBoundaryConnections(
    void*    fileHandle,
    int32_t  zone,
    int64_t* numConnections);

EXTERNC tecio_API int32_t tecZoneVarGetDoubleValues(
    void*   fileHandle,
    int32_t zone,
    int32_t var,
    int64_t startIndex,
    int64_t numValues,
    double* values);

EXTERNC tecio_API int32_t tecZoneVarGetFloatValues(
    void*   fileHandle,
    int32_t zone,
    int32_t var,
    int64_t startIndex,
    int64_t numValues,
    float*  values);

EXTERNC tecio_API int32_t tecZoneVarGetInt16Values(
    void*    fileHandle,
    int32_t  zone,
    int32_t  var,
    int64_t  startIndex,
    int64_t  numValues,
    int16_t* values);

EXTERNC tecio_API int32_t tecZoneVarGetInt32Values(
    void*    fileHandle,
    int32_t  zone,
    int32_t  var,
    int64_t  startIndex,
    int64_t  numValues,
    int32_t* values);

EXTERNC tecio_API int32_t tecZoneVarGetNumValues(
    void*    fileHandle,
    int32_t  zone,
    int32_t  var,
    int64_t* numValues);

EXTERNC tecio_API int32_t tecZoneVarGetSharedZone(
    void*    fileHandle,
    int32_t  zone,
    int32_t  var,
    int32_t* sharedZone);

EXTERNC tecio_API int32_t tecZoneVarGetType(
    void*    fileHandle,
    int32_t  zone,
    int32_t  var,
    int32_t* type);

EXTERNC tecio_API int32_t tecZoneVarGetUInt8Values(
    void*    fileHandle,
    int32_t  zone,
    int32_t  var,
    int64_t  startIndex,
    int64_t  numValues,
    uint8_t* values);

EXTERNC tecio_API int32_t tecZoneVarGetValueLocation(
    void*    fileHandle,
    int32_t  zone,
    int32_t  var,
    int32_t* location);

EXTERNC tecio_API int32_t tecZoneVarIsPassive(
    void*    fileHandle,
    int32_t  zone,
    int32_t  var,
    int32_t* isPassive);

/* Older routines */

#if !defined CRAY
#  define TECINI142       tecini142
#  define TECZNE142       teczne142
#  define TECDAT142       tecdat142
#  define TECDATD142      tecdatd142
#  define TECDATF142      tecdatf142
#  define TECNOD142       tecnod142
#  define TECNODE142      tecnode142
#  define TECGEO142       tecgeo142
#  define TECTXT142       tectxt142
#  define TECLAB142       teclab142
#  define TECFIL142       tecfil142
#  define TECFOREIGN142   tecforeign142
#  define TECFLUSH142     tecflush142
#  define TECEND142       tecend142
#  define TECUSR142       tecusr142
#  define TECAUXSTR142    tecauxstr142
#  define TECZAUXSTR142   teczauxstr142
#  define TECVAUXSTR142   tecvauxstr142
#  define TECFACE142      tecface142
#  define TECPOLY142      tecpoly142
#  define TECPOLYFACE142  tecpolyface142
#  define TECPOLYBCONN142 tecpolybconn142
/*
* SZL-only API:
*/
#  define TECFEPTN142   tecfeptn142
#  define TECIJKPTN142  tecijkptn142

/*
* SZL MPI-only APIs:
*/
#  define TECMPIINIT142 tecmpiinit142
#  define TECZNEMAP142  tecznemap142

/*
* Older API versions:
*/
#  define TECINI112       tecini112
#  define TECZNE112       teczne112
#  define TECDAT112       tecdat112
#  define TECNOD112       tecnod112
#  define TECNODE112      tecnode112
#  define TECGEO112       tecgeo112
#  define TECTXT112       tectxt112
#  define TECLAB112       teclab112
#  define TECFIL112       tecfil112
#  define TECFOREIGN112   tecforeign112
#  define TECEND112       tecend112
#  define TECUSR112       tecusr112
#  define TECAUXSTR112    tecauxstr112
#  define TECZAUXSTR112   teczauxstr112
#  define TECVAUXSTR112   tecvauxstr112
#  define TECFACE112      tecface112
#  define TECPOLY112      tecpoly112
#  define TECPOLYFACE112  tecpolyface112
#  define TECPOLYBCONN112 tecpolybconn112

#  define TECINI111     tecini111
#  define TECZNE111     teczne111
#  define TECDAT111     tecdat111
#  define TECNOD111     tecnod111
#  define TECGEO111     tecgeo111
#  define TECTXT111     tectxt111
#  define TECLAB111     teclab111
#  define TECFIL111     tecfil111
#  define TECFOREIGN111 tecforeign111
#  define TECEND111     tecend111
#  define TECUSR111     tecusr111
#  define TECAUXSTR111  tecauxstr111
#  define TECZAUXSTR111 teczauxstr111
#  define TECVAUXSTR111 tecvauxstr111
#  define TECFACE111    tecface111
#  define TECPOLY111    tecpoly111

#  define TECINI110     tecini110
#  define TECZNE110     teczne110
#  define TECDAT110     tecdat110
#  define TECNOD110     tecnod110
#  define TECGEO110     tecgeo110
#  define TECTXT110     tectxt110
#  define TECLAB110     teclab110
#  define TECFIL110     tecfil110
#  define TECFOREIGN110 tecforeign110
#  define TECEND110     tecend110
#  define TECUSR110     tecusr110
#  define TECAUXSTR110  tecauxstr110
#  define TECZAUXSTR110 teczauxstr110
#  define TECVAUXSTR110 tecvauxstr110
#  define TECFACE110    tecface110

#  define TECINI100     tecini100
#  define TECZNE100     teczne100
#  define TECDAT100     tecdat100
#  define TECNOD100     tecnod100
#  define TECGEO100     tecgeo100
#  define TECTXT100     tectxt100
#  define TECLAB100     teclab100
#  define TECFIL100     tecfil100
#  define TECFOREIGN100 tecforeign100
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

/*
 * TecIO API version 142 introduced the ability write both PLT and SZL formatted files.
 * No polyhedral support for SZL files.
 */

EXTERNC tecio_API INTEGER4 STDCALL TECINI142(
    char const*     Title,
    char const*     Variables,
    char const*     FName,
    char const*     ScratchDir,
    INTEGER4 const* FileFormat,
    INTEGER4 const* FileType,
    INTEGER4 const* Debug,
    INTEGER4 const* VIsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECZNE142(
    char const*     ZoneTitle,
    INTEGER4 const* ZoneType,
    INTEGER4 const* IMxOrNumPts,
    INTEGER4 const* JMxOrNumElements,
    INTEGER4 const* KMxOrNumFaces,
    INTEGER4 const* ICellMx,
    INTEGER4 const* JCellMx,
    INTEGER4 const* KCellMx,
    double const*   SolutionTime,
    INTEGER4 const* StrandID,
    INTEGER4 const* ParentZone,
    INTEGER4 const* IsBlock,
    INTEGER4 const* NumFaceConnections,
    INTEGER4 const* FaceNeighborMode,
    INTEGER4 const* TotalNumFaceNodes,
    INTEGER4 const* NumConnectedBoundaryFaces,
    INTEGER4 const* TotalNumBoundaryConnections,
    INTEGER4 const* PassiveVarList,
    INTEGER4 const* ValueLocation,
    INTEGER4 const* ShareVarFromZone,
    INTEGER4 const* ShareConnectivityFromZone);

EXTERNC tecio_API INTEGER4 STDCALL TECDAT142(
    INTEGER4 const* N,
    void const*     FieldData,
    INTEGER4 const* IsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECDATD142(
    INTEGER4 const* N,
    double const*   FieldData);

EXTERNC tecio_API INTEGER4 STDCALL TECDATF142(
    INTEGER4 const* N,
    float const*    FieldData);

EXTERNC tecio_API INTEGER4 STDCALL TECNOD142(INTEGER4 const* NData);

EXTERNC tecio_API INTEGER4 STDCALL TECNODE142(
    INTEGER4 const* N,
    INTEGER4 const* NData);

EXTERNC tecio_API INTEGER4 STDCALL TECFLUSH142(
    INTEGER4 const* NumZonesToRetain,
    INTEGER4 const* ZonesToRetain);

EXTERNC tecio_API INTEGER4 STDCALL TECEND142(void);

EXTERNC tecio_API INTEGER4 STDCALL TECLAB142(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECUSR142(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECGEO142(
    double const*   XPos,
    double const*   YPos,
    double const*   ZPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* Color,
    INTEGER4 const* FillColor,
    INTEGER4 const* IsFilled,
    INTEGER4 const* GeomType,
    INTEGER4 const* LinePattern,
    double const*   PatternLength,
    double const*   LineThickness,
    INTEGER4 const* NumEllipsePts,
    INTEGER4 const* ArrowheadStyle,
    INTEGER4 const* ArrowheadAttachment,
    double const*   ArrowheadSize,
    double const*   ArrowheadAngle,
    INTEGER4 const* Scope,
    INTEGER4 const* Clipping,
    INTEGER4 const* NumSegments,
    INTEGER4 const* NumSegPts,
    float const*    XGeomData,
    float const*    YGeomData,
    float const*    ZGeomData,
    char const*     mfc);

EXTERNC tecio_API INTEGER4 STDCALL TECTXT142(
    double const*   XOrThetaPos,
    double const*   YOrRPos,
    double const*   ZOrUnusedPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* BFont,
    INTEGER4 const* FontHeightUnits,
    double const*   FontHeight,
    INTEGER4 const* BoxType,
    double const*   BoxMargin,
    double const*   BoxLineThickness,
    INTEGER4 const* BoxColor,
    INTEGER4 const* BoxFillColor,
    double const*   Angle,
    INTEGER4 const* Anchor,
    double const*   LineSpacing,
    INTEGER4 const* TextColor,
    INTEGER4 const* Scope,
    INTEGER4 const* Clipping,
    char const*     String,
    char const*     mfc);

EXTERNC tecio_API INTEGER4 STDCALL TECFIL142(INTEGER4 const* F);

EXTERNC tecio_API void STDCALL TECFOREIGN142(INTEGER4 const* OutputForeignByteOrder);

EXTERNC tecio_API INTEGER4 STDCALL TECAUXSTR142(
    char const* Name,
    char const* Value);

EXTERNC tecio_API INTEGER4 STDCALL TECZAUXSTR142(
    char const* Name,
    char const* Value);

EXTERNC tecio_API INTEGER4 STDCALL TECVAUXSTR142(
    INTEGER4 const* Var,
    char const*     Name,
    char const*     Value);

EXTERNC tecio_API INTEGER4 STDCALL TECFACE142(INTEGER4 const* FaceConnections);

EXTERNC tecio_API INTEGER4 STDCALL TECPOLY142(
    INTEGER4 const* FaceNodeCounts,
    INTEGER4 const* FaceNodes,
    INTEGER4 const* FaceLeftElems,
    INTEGER4 const* FaceRightElems,
    INTEGER4 const* FaceBndryConnectionCounts,
    INTEGER4 const* FaceBndryConnectionElems,
    INTEGER4 const* FaceBndryConnectionZones);

EXTERNC tecio_API INTEGER4 STDCALL TECPOLYFACE142(
    INTEGER4 const* NumFaces,
    INTEGER4 const* FaceNodeCounts,
    INTEGER4 const* FaceNodes,     
    INTEGER4 const* FaceLeftElems, 
    INTEGER4 const* FaceRightElems);

EXTERNC tecio_API INTEGER4 STDCALL TECPOLYBCONN142(
    INTEGER4 const* NumBndryFaces,
    INTEGER4 const* FaceBndryConnectionCounts,
    INTEGER4 const* FaceBndryConnectionElems, 
    INTEGER4 const* FaceBndryConnectionZones);

/* SZL-only APIs: */
EXTERNC tecio_API INTEGER4 STDCALL TECFEPTN142(
    INTEGER4 const* partition,
    INTEGER4 const* numnodes,
    INTEGER4 const* numcells,
    INTEGER4 const* ngnodes,
    INTEGER4 const* gnodes,
    INTEGER4 const* gnpartitions,
    INTEGER4 const* gnpnodes,
    INTEGER4 const* ngcells,
    INTEGER4 const* gcells);

EXTERNC tecio_API INTEGER4 STDCALL TECIJKPTN142(
    INTEGER4 const* partition,
    INTEGER4 const* imin,
    INTEGER4 const* jmin,
    INTEGER4 const* kmin,
    INTEGER4 const* imax,
    INTEGER4 const* jmax,
    INTEGER4 const* kmax);

/* SZL MPI-only APIs: */
EXTERNC tecio_API INTEGER4 STDCALL TECMPIINIT142(
    void* communicator, /* MPI_Comm */
    INTEGER4 const* mainrank);

EXTERNC tecio_API INTEGER4 STDCALL TECZNEMAP142(
    INTEGER4 const* npartitions,
    INTEGER4 const* ptnranks);

/*
 *  V11.3 tecio functions
 */

EXTERNC tecio_API INTEGER4 STDCALL TECINI112(
    char const*     Title,
    char const*     Variables,
    char const*     FName,
    char const*     ScratchDir,
    INTEGER4 const* FileType,
    INTEGER4 const* Debug,
    INTEGER4 const* VIsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECZNE112(
    char const*     ZoneTitle,
    INTEGER4 const* ZoneType,
    INTEGER4 const* IMxOrNumPts,
    INTEGER4 const* JMxOrNumElements,
    INTEGER4 const* KMxOrNumFaces,
    INTEGER4 const* ICellMx,
    INTEGER4 const* JCellMx,
    INTEGER4 const* KCellMx,
    double const*   SolutionTime,
    INTEGER4 const* StrandID,
    INTEGER4 const* ParentZone,
    INTEGER4 const* IsBlock,
    INTEGER4 const* NumFaceConnections,
    INTEGER4 const* FaceNeighborMode,
    INTEGER4 const* TotalNumFaceNodes,
    INTEGER4 const* NumConnectedBoundaryFaces,
    INTEGER4 const* TotalNumBoundaryConnections,
    INTEGER4 const* PassiveVarList,
    INTEGER4 const* ValueLocation,
    INTEGER4 const* ShareVarFromZone,
    INTEGER4 const* ShareConnectivityFromZone);

EXTERNC tecio_API INTEGER4 STDCALL TECDAT112(
    INTEGER4 const* N,
    void const*     FieldData,
    INTEGER4 const* IsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECNOD112(INTEGER4 const* NData);

EXTERNC tecio_API INTEGER4 STDCALL TECNODE112(
    INTEGER4 const* N,
    INTEGER4 const* NData);

EXTERNC tecio_API INTEGER4 STDCALL TECEND112(void);

EXTERNC tecio_API INTEGER4 STDCALL TECLAB112(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECUSR112(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECGEO112(
    double const*   XPos,
    double const*   YPos,
    double const*   ZPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* Color,
    INTEGER4 const* FillColor,
    INTEGER4 const* IsFilled,
    INTEGER4 const* GeomType,
    INTEGER4 const* LinePattern,
    double const*   PatternLength,
    double const*   LineThickness,
    INTEGER4 const* NumEllipsePts,
    INTEGER4 const* ArrowheadStyle,
    INTEGER4 const* ArrowheadAttachment,
    double const*   ArrowheadSize,
    double const*   ArrowheadAngle,
    INTEGER4 const* Scope,
    INTEGER4 const* Clipping,
    INTEGER4 const* NumSegments,
    INTEGER4 const* NumSegPts,
    float const*    XGeomData,
    float const*    YGeomData,
    float const*    ZGeomData,
    char const*     mfc);

EXTERNC tecio_API INTEGER4 STDCALL TECTXT112(
    double const*   XOrThetaPos,
    double const*   YOrRPos,
    double const*   ZOrUnusedPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* BFont,
    INTEGER4 const* FontHeightUnits,
    double const*   FontHeight,
    INTEGER4 const* BoxType,
    double const*   BoxMargin,
    double const*   BoxLineThickness,
    INTEGER4 const* BoxColor,
    INTEGER4 const* BoxFillColor,
    double const*   Angle,
    INTEGER4 const* Anchor,
    double const*   LineSpacing,
    INTEGER4 const* TextColor,
    INTEGER4 const* Scope,
    INTEGER4 const* Clipping,
    char const*     String,
    char const*     mfc);

EXTERNC tecio_API INTEGER4 STDCALL TECFIL112(INTEGER4 const* F);

EXTERNC tecio_API void STDCALL TECFOREIGN112(INTEGER4 const* OutputForeignByteOrder);

EXTERNC tecio_API INTEGER4 STDCALL TECAUXSTR112(
    char const* Name,
    char const* Value);

EXTERNC tecio_API INTEGER4 STDCALL TECZAUXSTR112(
    char const* Name,
    char const* Value);

EXTERNC tecio_API INTEGER4 STDCALL TECVAUXSTR112(
    INTEGER4 const* Var,
    char const*     Name,
    char const*     Value);

EXTERNC tecio_API INTEGER4 STDCALL TECFACE112(INTEGER4 const* FaceConnections);

EXTERNC tecio_API INTEGER4 STDCALL TECPOLY112(
    INTEGER4 const* FaceNodeCounts,
    INTEGER4 const* FaceNodes,
    INTEGER4 const* FaceLeftElems,
    INTEGER4 const* FaceRightElems,
    INTEGER4 const* FaceBndryConnectionCounts,
    INTEGER4 const* FaceBndryConnectionElems,
    INTEGER4 const* FaceBndryConnectionZones);

EXTERNC tecio_API INTEGER4 STDCALL TECPOLYFACE112(
    INTEGER4 const* NumFaces,
    INTEGER4 const* FaceNodeCounts,
    INTEGER4 const* FaceNodes,     
    INTEGER4 const* FaceLeftElems, 
    INTEGER4 const* FaceRightElems);

EXTERNC tecio_API INTEGER4 STDCALL TECPOLYBCONN112(
    INTEGER4 const* NumBndryFaces,
    INTEGER4 const* FaceBndryConnectionCounts,
    INTEGER4 const* FaceBndryConnectionElems, 
    INTEGER4 const* FaceBndryConnectionZones);

/*
 *  V11.1 tecio functions   TODO (JN): Tecplot's version is still in flux so the .1 may change
 */

EXTERNC tecio_API INTEGER4 STDCALL TECINI111(
    char const*     Title,
    char const*     Variables,
    char const*     FName,
    char const*     ScratchDir,
    INTEGER4 const* FileType,
    INTEGER4 const* Debug,
    INTEGER4 const* VIsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECZNE111(
    char const*     ZoneTitle,
    INTEGER4 const* ZoneType,
    INTEGER4 const* IMxOrNumPts,
    INTEGER4 const* JMxOrNumElements,
    INTEGER4 const* KMxOrNumFaces,
    INTEGER4 const* ICellMx,
    INTEGER4 const* JCellMx,
    INTEGER4 const* KCellMx,
    double const*   SolutionTime,
    INTEGER4 const* StrandID,
    INTEGER4 const* ParentZone,
    INTEGER4 const* IsBlock,
    INTEGER4 const* NumFaceConnections,
    INTEGER4 const* FaceNeighborMode,
    INTEGER4 const* TotalNumFaceNodes,
    INTEGER4 const* NumConnectedBoundaryFaces,
    INTEGER4 const* TotalNumBoundaryConnections,
    INTEGER4 const* PassiveVarList,
    INTEGER4 const* ValueLocation,
    INTEGER4 const* ShareVarFromZone,
    INTEGER4 const* ShareConnectivityFromZone);

EXTERNC tecio_API INTEGER4 STDCALL TECDAT111(
    INTEGER4 const* N,
    void const*     FieldData,
    INTEGER4 const* IsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECNOD111(INTEGER4 const* NData);

EXTERNC tecio_API INTEGER4 STDCALL TECEND111(void);

EXTERNC tecio_API INTEGER4 STDCALL TECLAB111(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECUSR111(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECGEO111(
    double const*   XPos,
    double const*   YPos,
    double const*   ZPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* Color,
    INTEGER4 const* FillColor,
    INTEGER4 const* IsFilled,
    INTEGER4 const* GeomType,
    INTEGER4 const* LinePattern,
    double const*   PatternLength,
    double const*   LineThickness,
    INTEGER4 const* NumEllipsePts,
    INTEGER4 const* ArrowheadStyle,
    INTEGER4 const* ArrowheadAttachment,
    double const*   ArrowheadSize,
    double const*   ArrowheadAngle,
    INTEGER4 const* Scope,
    INTEGER4 const* Clipping,
    INTEGER4 const* NumSegments,
    INTEGER4 const* NumSegPts,
    float const*    XGeomData,
    float const*    YGeomData,
    float const*    ZGeomData,
    char const*     mfc);

EXTERNC tecio_API INTEGER4 STDCALL TECTXT111(
    double const*   XOrThetaPos,
    double const*   YOrRPos,
    double const*   ZOrUnusedPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* BFont,
    INTEGER4 const* FontHeightUnits,
    double const*   FontHeight,
    INTEGER4 const* BoxType,
    double const*   BoxMargin,
    double const*   BoxLineThickness,
    INTEGER4 const* BoxColor,
    INTEGER4 const* BoxFillColor,
    double const*   Angle,
    INTEGER4 const* Anchor,
    double const*   LineSpacing,
    INTEGER4 const* TextColor,
    INTEGER4 const* Scope,
    INTEGER4 const* Clipping,
    char const*     String,
    char const*     mfc);

EXTERNC tecio_API INTEGER4 STDCALL TECFIL111(INTEGER4 const* F);

EXTERNC tecio_API void STDCALL TECFOREIGN111(INTEGER4 const* OutputForeignByteOrder);

EXTERNC tecio_API INTEGER4 STDCALL TECAUXSTR111(
    char const* Name,
    char const* Value);

EXTERNC tecio_API INTEGER4 STDCALL TECZAUXSTR111(
    char const* Name,
    char const* Value);

EXTERNC tecio_API INTEGER4 STDCALL TECVAUXSTR111(
    INTEGER4 const* Var,
    char const*     Name,
    char const*     Value);

EXTERNC tecio_API INTEGER4 STDCALL TECFACE111(INTEGER4 const* FaceConnections);

EXTERNC tecio_API INTEGER4 STDCALL TECPOLY111(
    INTEGER4 const* FaceNodeCounts,
    INTEGER4 const* FaceNodes,
    INTEGER4 const* FaceLeftElems,
    INTEGER4 const* FaceRightElems,
    INTEGER4 const* FaceBndryConnectionCounts,
    INTEGER4 const* FaceBndryConnectionElems,
    INTEGER2 const* FaceBndryConnectionZones);


/*
 * V11 tecio functions
 */

EXTERNC tecio_API INTEGER4 STDCALL TECINI110(
    char const*     Title,
    char const*     Variables,
    char const*     FName,
    char const*     ScratchDir,
    INTEGER4 const* Debug,
    INTEGER4 const* VIsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECZNE110(
    char const*     ZoneTitle,
    INTEGER4 const* ZoneType,
    INTEGER4 const* IMxOrNumPts,
    INTEGER4 const* JMxOrNumElements,
    INTEGER4 const* KMxOrNumFaces,
    INTEGER4 const* ICellMx,
    INTEGER4 const* JCellMx,
    INTEGER4 const* KCellMx,
    double const*   SolutionTime,
    INTEGER4 const* StrandID,
    INTEGER4 const* ParentZone,
    INTEGER4 const* IsBlock,
    INTEGER4 const* NumFaceConnections,
    INTEGER4 const* FaceNeighborMode,
    INTEGER4 const* PassiveVarList,
    INTEGER4 const* ValueLocation,
    INTEGER4 const* ShareVarFromZone,
    INTEGER4 const* ShareConnectivityFromZone);

EXTERNC tecio_API INTEGER4 STDCALL TECDAT110(
    INTEGER4 const* N,
    void const*     FieldData,
    INTEGER4 const* IsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECNOD110(INTEGER4 const* NData);

EXTERNC tecio_API INTEGER4 STDCALL TECEND110(void);

EXTERNC tecio_API INTEGER4 STDCALL TECLAB110(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECUSR110(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECGEO110(
    double const*   XPos,
    double const*   YPos,
    double const*   ZPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* Color,
    INTEGER4 const* FillColor,
    INTEGER4 const* IsFilled,
    INTEGER4 const* GeomType,
    INTEGER4 const* LinePattern,
    double const*   PatternLength,
    double const*   LineThickness,
    INTEGER4 const* NumEllipsePts,
    INTEGER4 const* ArrowheadStyle,
    INTEGER4 const* ArrowheadAttachment,
    double const*   ArrowheadSize,
    double const*   ArrowheadAngle,
    INTEGER4 const* Scope,
    INTEGER4 const* Clipping,
    INTEGER4 const* NumSegments,
    INTEGER4 const* NumSegPts,
    float const*    XGeomData,
    float const*    YGeomData,
    float const*    ZGeomData,
    char const*     mfc);

EXTERNC tecio_API INTEGER4 STDCALL TECTXT110(
    double const*   XOrThetaPos,
    double const*   YOrRPos,
    double const*   ZOrUnusedPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* BFont,
    INTEGER4 const* FontHeightUnits,
    double const*   FontHeight,
    INTEGER4 const* BoxType,
    double const*   BoxMargin,
    double const*   BoxLineThickness,
    INTEGER4 const* BoxColor,
    INTEGER4 const* BoxFillColor,
    double const*   Angle,
    INTEGER4 const* Anchor,
    double const*   LineSpacing,
    INTEGER4 const* TextColor,
    INTEGER4 const* Scope,
    INTEGER4 const* Clipping,
    char const*     String,
    char const*     mfc);

EXTERNC tecio_API void STDCALL TECFOREIGN110(INTEGER4 const* OutputForeignByteOrder);

EXTERNC tecio_API INTEGER4 STDCALL TECFIL110(INTEGER4 const* F);

EXTERNC tecio_API INTEGER4 STDCALL TECAUXSTR110(
    char const* Name,
    char const* Value);

EXTERNC tecio_API INTEGER4 STDCALL TECZAUXSTR110(
    char const* Name,
    char const* Value);

EXTERNC tecio_API INTEGER4 STDCALL TECVAUXSTR110(
    INTEGER4 const* Var,
    char const*     Name,
    char const*     Value);

EXTERNC tecio_API INTEGER4 STDCALL TECFACE110(INTEGER4 const* FaceConnections);


/*
 * V10 tecio functions kept for backward compatability.
 */

EXTERNC tecio_API INTEGER4 STDCALL TECINI100(
    char const*     Title,
    char const*     Variables,
    char const*     FName,
    char const*     ScratchDir,
    INTEGER4 const* Debug,
    INTEGER4 const* VIsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECZNE100(
    char const*     ZoneTitle,
    INTEGER4 const* ZoneType,
    INTEGER4 const* IMxOrNumPts,
    INTEGER4 const* JMxOrNumElements,
    INTEGER4 const* KMxOrNumFaces,
    INTEGER4 const* ICellMx,
    INTEGER4 const* JCellMx,
    INTEGER4 const* KCellMx,
    INTEGER4 const* IsBlock,
    INTEGER4 const* NumFaceConnections,
    INTEGER4 const* FaceNeighborMode,
    INTEGER4 const* ValueLocation,
    INTEGER4 const* ShareVarFromZone,
    INTEGER4 const* ShareConnectivityFromZone);

EXTERNC tecio_API INTEGER4 STDCALL TECDAT100(
    INTEGER4 const* N,
    void const*     FieldData,
    INTEGER4 const* IsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECNOD100(INTEGER4 const* NData);

EXTERNC tecio_API INTEGER4 STDCALL TECEND100(void);

EXTERNC tecio_API INTEGER4 STDCALL TECLAB100(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECUSR100(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECGEO100(
    double const*   XPos,
    double const*   YPos,
    double const*   ZPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* Color,
    INTEGER4 const* FillColor,
    INTEGER4 const* IsFilled,
    INTEGER4 const* GeomType,
    INTEGER4 const* LinePattern,
    double const*   PatternLength,
    double const*   LineThickness,
    INTEGER4 const* NumEllipsePts,
    INTEGER4 const* ArrowheadStyle,
    INTEGER4 const* ArrowheadAttachment,
    double const*   ArrowheadSize,
    double const*   ArrowheadAngle,
    INTEGER4 const* Scope,
    INTEGER4 const* Clipping,
    INTEGER4 const* NumSegments,
    INTEGER4 const* NumSegPts,
    float const*    XGeomData,
    float const*    YGeomData,
    float const*    ZGeomData,
    char const*     mfc);

EXTERNC tecio_API INTEGER4 STDCALL TECTXT100(
    double const*   XOrThetaPos,
    double const*   YOrRPos,
    double const*   ZOrUnusedPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* BFont,
    INTEGER4 const* FontHeightUnits,
    double const*   FontHeight,
    INTEGER4 const* BoxType,
    double const*   BoxMargin,
    double const*   BoxLineThickness,
    INTEGER4 const* BoxColor,
    INTEGER4 const* BoxFillColor,
    double const*   Angle,
    INTEGER4 const* Anchor,
    double const*   LineSpacing,
    INTEGER4 const* TextColor,
    INTEGER4 const* Scope,
    INTEGER4 const* Clipping,
    char const*     String,
    char const*     mfc);

EXTERNC tecio_API INTEGER4 STDCALL TECFIL100(INTEGER4 const* F);

EXTERNC tecio_API void STDCALL TECFOREIGN100(INTEGER4 const* OutputForeignByteOrder);

EXTERNC tecio_API INTEGER4 STDCALL TECAUXSTR100(
    char const* Name,
    char const* Value);

EXTERNC tecio_API INTEGER4 STDCALL TECZAUXSTR100(
    char const* Name,
    char const* Value);

EXTERNC tecio_API INTEGER4 STDCALL TECVAUXSTR100(
    INTEGER4 const* Var,
    char const*     Name,
    char const*     Value);

EXTERNC tecio_API INTEGER4 STDCALL TECFACE100(INTEGER4 const* FaceConnections);

/* Old V9 functions retained for backward compatibility */

EXTERNC tecio_API INTEGER4 STDCALL TECINI(
    char const*     Title,
    char const*     Variables,
    char const*     FName,
    char const*     ScratchDir,
    INTEGER4 const* Debug,
    INTEGER4 const* VIsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECZNE(
    char const*     ZoneTitle,
    INTEGER4 const* IMx,
    INTEGER4 const* JMx,
    INTEGER4 const* KMx,
    char const*     ZFormat,
    char const*     DupList);

EXTERNC tecio_API INTEGER4 STDCALL TECDAT(
    INTEGER4 const* N,
    void const*     FieldData,
    INTEGER4 const* IsDouble);

EXTERNC tecio_API INTEGER4 STDCALL TECNOD(INTEGER4 const* NData);

EXTERNC tecio_API INTEGER4 STDCALL TECEND(void);

EXTERNC tecio_API INTEGER4 STDCALL TECLAB(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECUSR(char const* S);

EXTERNC tecio_API INTEGER4 STDCALL TECGEO(
    double const*   XPos,
    double const*   YPos,
    double const*   ZPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* Color,
    INTEGER4 const* FillColor,
    INTEGER4 const* IsFilled,
    INTEGER4 const* GeomType,
    INTEGER4 const* LinePattern,
    double const*   PatternLength,
    double const*   LineThickness,
    INTEGER4 const* NumEllipsePts,
    INTEGER4 const* ArrowheadStyle,
    INTEGER4 const* ArrowheadAttachment,
    double const*   ArrowheadSize,
    double const*   ArrowheadAngle,
    INTEGER4 const* Scope,
    INTEGER4 const* NumSegments,
    INTEGER4 const* NumSegPts,
    float const*    XGeomData,
    float const*    YGeomData,
    float const*    ZGeomData,
    char const*     mfc);

EXTERNC tecio_API INTEGER4 STDCALL TECTXT(
    double const*   XPos,
    double const*   YPos,
    INTEGER4 const* PosCoordMode,
    INTEGER4 const* AttachToZone,
    INTEGER4 const* Zone,
    INTEGER4 const* BFont,
    INTEGER4 const* FontHeightUnits,
    double const*   FontHeight,
    INTEGER4 const* BoxType,
    double const*   BoxMargin,
    double const*   BoxLineThickness,
    INTEGER4 const* BoxColor,
    INTEGER4 const* BoxFillColor,
    double const*   Angle,
    INTEGER4 const* Anchor,
    double const*   LineSpacing,
    INTEGER4 const* TextColor,
    INTEGER4 const* Scope,
    char const*     Text,
    char const*     mfc);

EXTERNC tecio_API INTEGER4 STDCALL TECFIL(INTEGER4 const* F);
