
/**
 * Read the .szplt file named in argv[1] and write it out to a .szplt file named in argv[2]
 */

#include "TECIO.h"

#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

int32_t countFaceConnections(int32_t const* faceConnections, int64_t numFaceValues, int32_t mode);

// We retrieve the return value from each tec call, but ignore it. Calls return 0 for success, non-zero for failure.
int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " infilename outfilename" << std::endl;
        return 1;
    }

    try
    {
        // Open a .szplt file for reading with TECIO
        void* inputFileHandle = NULL;
        int i = tecFileReaderOpen((char const*)(argv[1]), &inputFileHandle);

        // Read the characteristics of the data set
        char* dataSetTitle = NULL;
        i = tecDataSetGetTitle(inputFileHandle, &dataSetTitle);

        int32_t numVars;
        i = tecDataSetGetNumVars(inputFileHandle, &numVars);
        std::ostringstream outputStream;
        for (int32_t var = 1; var <= numVars; ++var)
        {
            char* name = NULL;
            i = tecVarGetName(inputFileHandle, var, &name);
            outputStream << name;
            if (var < numVars)
                outputStream << ',';
            tecStringFree(&name);
        }

        int32_t fileType;
        i = tecFileGetType(inputFileHandle, &fileType);

        int32_t numZones;
        i = tecDataSetGetNumZones(inputFileHandle, &numZones);

        int32_t isDouble = 0;
        int32_t const FieldDataType_Double = 2; // ...TecUtil types are not available to TecIO so redefine
        if (numZones > 0)
        {
            int32_t type;
            i = tecZoneVarGetType(inputFileHandle, 1, 1, &type);
            if (type == FieldDataType_Double)
                isDouble = 1;
        }

        // Begin writing a new .szplt file
        int32_t fileFormat = 1; // .szplt
        int32_t outputDebugInfo = 1;
        void* outputFileHandle = NULL;
        // Default var type is float (1), but we'll explicitly set all var types when we create each zone.
        i = tecFileWriterOpen(argv[2], dataSetTitle, outputStream.str().c_str(), fileFormat, fileType, 1, NULL, &outputFileHandle);
        i = tecFileSetDiagnosticsLevel(outputFileHandle, outputDebugInfo);
        tecStringFree(&dataSetTitle);

        for (int32_t inputZone = 1; inputZone <= numZones; ++inputZone)
        {
            // Retrieve zone characteristics
            int32_t zoneType;
            i = tecZoneGetType(inputFileHandle, inputZone, &zoneType);
            if (zoneType == 6 || zoneType == 7)
                throw std::runtime_error("Unsupported zone type.");

            char* zoneTitle = NULL;
            i = tecZoneGetTitle(inputFileHandle, inputZone, &zoneTitle);

            int64_t iMax, jMax, kMax;
            i = tecZoneGetIJK(inputFileHandle, inputZone, &iMax, &jMax, &kMax);

            std::vector<int32_t> varTypes(numVars);
            std::vector<int32_t> passiveVarList(numVars);
            std::vector<int32_t> valueLocation(numVars);
            std::vector<int32_t> shareVarFromZone(numVars);
            for (int32_t var = 1; var <= numVars; ++var)
            {
                i = tecZoneVarGetType(inputFileHandle, inputZone, var, &varTypes[var - 1]);
                i = tecZoneVarGetSharedZone(inputFileHandle, inputZone, var, &shareVarFromZone[var - 1]);
                i = tecZoneVarGetValueLocation(inputFileHandle, inputZone, var, &valueLocation[var - 1]);
                i = tecZoneVarIsPassive(inputFileHandle, inputZone, var, &passiveVarList[var - 1]);
            }

            int32_t shareConnectivityFromZone;
            i = tecZoneConnectivityGetSharedZone(inputFileHandle, inputZone, &shareConnectivityFromZone);

            int32_t faceNeighborMode;
            i = tecZoneFaceNbrGetMode(inputFileHandle, inputZone, &faceNeighborMode);

            int64_t numFaceConnections;
            i = tecZoneFaceNbrGetNumConnections(inputFileHandle, inputZone, &numFaceConnections);

            int32_t outputZone;
            if (zoneType == 0)
                i = tecZoneCreateIJK(outputFileHandle, zoneTitle, iMax, jMax, kMax, &varTypes[0],
                    &shareVarFromZone[0], &valueLocation[0], &passiveVarList[0],
                    shareConnectivityFromZone, numFaceConnections, faceNeighborMode, &outputZone);
            else
                i = tecZoneCreateFE(outputFileHandle, zoneTitle, zoneType, iMax, jMax, &varTypes[0],
                    &shareVarFromZone[0], &valueLocation[0], &passiveVarList[0],
                    shareConnectivityFromZone, numFaceConnections, faceNeighborMode, &outputZone);

            double solutionTime;
            int32_t strandID;
            i = tecZoneGetSolutionTime(inputFileHandle, inputZone, &solutionTime);
            i = tecZoneGetStrandID(inputFileHandle, inputZone, &strandID);
            if (solutionTime != 0.0 || strandID != 0)
                i = tecZoneSetUnsteadyOptions(outputFileHandle, outputZone, solutionTime, strandID);

            int32_t parentZone;
            i = tecZoneGetParentZone(inputFileHandle, inputZone, &parentZone);
            if (parentZone != 0)
                i = tecZoneSetParentZone(outputFileHandle, outputZone, parentZone);

            tecStringFree(&zoneTitle);

            // Retrieve zone data and send to the output file
            for (int32_t var = 1; var <= numVars; ++var)
            {
                if (passiveVarList[var - 1] == 0 && shareVarFromZone[var - 1] == 0)
                {
                    int64_t numValues;
                    i = tecZoneVarGetNumValues(inputFileHandle, inputZone, var, &numValues);
                    // For large zones, could "chunk" this input/output--read/write the var in pieces instead of all at once
                    switch((FieldDataType_e)varTypes[var - 1])
                    {
                    case FieldDataType_Float:
                        {
                            std::vector<float> values(numValues);
                            i = tecZoneVarGetFloatValues(inputFileHandle, inputZone, var, 1, numValues, &values[0]);
                            i = tecZoneVarWriteFloatValues(outputFileHandle, outputZone, var, 0, numValues, &values[0]);
                        }
                        break;
                    case FieldDataType_Double:
                        {
                            std::vector<double> values(numValues);
                            i = tecZoneVarGetDoubleValues(inputFileHandle, inputZone, var, 1, numValues, &values[0]);
                            i = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, var, 0, numValues, &values[0]);
                        }
                        break;
                    case FieldDataType_Int32:
                        {
                            std::vector<int32_t> values(numValues);
                            i = tecZoneVarGetInt32Values(inputFileHandle, inputZone, var, 1, numValues, &values[0]);
                            i = tecZoneVarWriteInt32Values(outputFileHandle, outputZone, var, 0, numValues, &values[0]);
                        }
                        break;
                    case FieldDataType_Int16:
                        {
                            std::vector<int16_t> values(numValues);
                            i = tecZoneVarGetInt16Values(inputFileHandle, inputZone, var, 1, numValues, &values[0]);
                            i = tecZoneVarWriteInt16Values(outputFileHandle, outputZone, var, 0, numValues, &values[0]);
                        }
                        break;
                    case FieldDataType_Byte:
                        {
                            std::vector<uint8_t> values(numValues);
                            i = tecZoneVarGetUInt8Values(inputFileHandle, inputZone, var, 1, numValues, &values[0]);
                            i = tecZoneVarWriteUInt8Values(outputFileHandle, outputZone, var, 0, numValues, &values[0]);
                        }
                        break;
                    default:
                        i = -1;
                        break;
                    }
                }
            }

            // Write zone face neighbors, if any
            if (numFaceConnections > 0)
            {
                int64_t numFaceValues;
                i = tecZoneFaceNbrGetNumValues(inputFileHandle, inputZone, &numFaceValues);
                int32_t are64Bit;
                i = tecZoneFaceNbrsAre64Bit(inputFileHandle, inputZone, &are64Bit);
                if (are64Bit)
                {
                    std::vector<int64_t> faceConnections(numFaceValues);
                    i = tecZoneFaceNbrGetConnections64(inputFileHandle, inputZone, &faceConnections[0]);
                    i = tecZoneFaceNbrWriteConnections64(outputFileHandle, outputZone, &faceConnections[0]);
                }
                else
                {
                    std::vector<int32_t> faceConnections(numFaceValues);
                    i = tecZoneFaceNbrGetConnections(inputFileHandle, inputZone, &faceConnections[0]);
                    i = tecZoneFaceNbrWriteConnections32(outputFileHandle, outputZone, &faceConnections[0]);
                }
            }

            // Retrieve zone node map, if any, and send to the output file
            if (zoneType != 0 && shareConnectivityFromZone == 0)
            {
                int64_t numValues;
                i = tecZoneNodeMapGetNumValues(inputFileHandle, inputZone, jMax, &numValues);
                int32_t is64Bit;
                i = tecZoneNodeMapIs64Bit(inputFileHandle, inputZone, &is64Bit);
                if (is64Bit)
                {
                    std::vector<int64_t> nodeMap(numValues);
                    i = tecZoneNodeMapGet64(inputFileHandle, inputZone, 1, jMax, &nodeMap[0]);
                    i = tecZoneNodeMapWrite64(outputFileHandle, outputZone, 0, 1, numValues, &nodeMap[0]);
                }
                else
                {
                    std::vector<int32_t> nodeMap(numValues);
                    i = tecZoneNodeMapGet(inputFileHandle, inputZone, 1, jMax, &nodeMap[0]);
                    i = tecZoneNodeMapWrite32(outputFileHandle, outputZone, 0, 1, numValues, &nodeMap[0]);
                }
            }

            // Retrieve and write any zone aux data
            int32_t numItems;
            i = tecZoneAuxDataGetNumItems(inputFileHandle, inputZone, &numItems);
            for (int32_t whichItem = 1; whichItem <= numItems; ++whichItem)
            {
                char* name = NULL;
                char* value = NULL;
                i = tecZoneAuxDataGetItem(inputFileHandle, inputZone, whichItem, &name, &value);
                i = tecZoneAddAuxData(outputFileHandle, outputZone, name, value);
                tecStringFree(&name);
                tecStringFree(&value);
            }
        }

        // Custom label sets
        int32_t numSets;
        i = tecCustomLabelsGetNumSets(inputFileHandle, &numSets);
        for (int32_t whichSet = 1; whichSet <= numSets; ++whichSet)
        {
            char* labelSet = NULL;
            i = tecCustomLabelsGetSet(inputFileHandle, whichSet, &labelSet); // Labels returned as "\"Mon\",\"Tues\",\"Wed\""
            i = tecCustomLabelsAddSet(outputFileHandle, labelSet);
            tecStringFree(&labelSet);
        }

        // Data set aux data
        int32_t numItems;
        i = tecDataSetAuxDataGetNumItems(inputFileHandle, &numItems);
        for (int32_t whichItem = 1; whichItem <= numItems; ++whichItem)
        {
            char* name = NULL;
            char* value = NULL;
            i = tecDataSetAuxDataGetItem(inputFileHandle, whichItem, &name, &value);
            i = tecDataSetAddAuxData(outputFileHandle, name, value);
            tecStringFree(&name);
            tecStringFree(&value);
        }

        // Var aux data
        for (int32_t var = 1; var <= numVars; ++var)
        {
            int32_t numItems;
            i = tecVarAuxDataGetNumItems(inputFileHandle, var, &numItems);
            for (int32_t whichItem = 1; whichItem <= numItems; ++whichItem)
            {
                char* name = NULL;
                char* value = NULL;
                i = tecVarAuxDataGetItem(inputFileHandle, var, whichItem, &name, &value);
                i = tecVarAddAuxData(outputFileHandle, var, name, value);
                tecStringFree(&name);
                tecStringFree(&value);
            }
        }

        // Geometries
        int32_t numGeoms;
        i = tecGeomGetNumGeoms(inputFileHandle, &numGeoms);
        for (int32_t geom = 1; geom <= numGeoms; ++geom)
        {
            // Retrieve required geom info and begin a new geom
            int32_t type;
            i = tecGeomGetType(inputFileHandle, geom, &type);

            double x, y, z;
            i = tecGeomGetAnchorPos(inputFileHandle, geom, &x, &y, &z);

            int32_t coordMode;
            i = tecGeomGetCoordMode(inputFileHandle, geom, &coordMode);

            std::vector<int32_t> numSegPts(1, 0);
            std::vector<double> xGeomData;
            std::vector<double> yGeomData;
            std::vector<double> zGeomData;
            int32_t numEllipsePoints = 2;
            int32_t numSegments;
            GeomType_e geomType = static_cast<GeomType_e>(type);
            switch (geomType)
            {
            case GeomType_LineSegs:
            case GeomType_LineSegs3D:
                i = tecGeomLineGetSegmentCount(inputFileHandle, geom, &numSegments);
                numSegPts.resize(numSegments);
                for (int32_t segment = 1; segment <= numSegments; ++segment)
                {
                    i = tecGeomLineSegmentGetPointCount(inputFileHandle, geom, segment, &numSegPts[segment - 1]);
                    for (int32_t index = 1; index <= numSegPts[segment - 1]; ++index)
                    {
                        double geomX, geomY, geomZ;
                        i = tecGeomLineGetPoint(inputFileHandle, geom, segment, index, &geomX, &geomY, &geomZ);
                        xGeomData.push_back(geomX);
                        yGeomData.push_back(geomY);
                        zGeomData.push_back(geomZ);
                    }
                }
                if (geomType == GeomType_LineSegs)
                {
                    if (numSegments == 1)
                        i = tecGeom2DLineSegmentsBegin(outputFileHandle, x, y, numSegPts[0], &xGeomData[0], &yGeomData[0], coordMode);
                    else
                        i = tecGeom2DMultiLineSegmentsBegin(outputFileHandle, x, y, numSegments, &numSegPts[0], &xGeomData[0], &yGeomData[0], coordMode);
                }
                else
                {
                    if (numSegments == 1)
                        i = tecGeom3DLineSegmentsBegin(outputFileHandle, x, y, z, numSegPts[0], &xGeomData[0], &yGeomData[0], &zGeomData[0]);
                    else
                        i = tecGeom3DMultiLineSegmentsBegin(outputFileHandle, x, y, z, numSegments, &numSegPts[0], &xGeomData[0], &yGeomData[0], &zGeomData[0]);
                }

                // The default is no arrowheads
                double arrowheadAngle;
                i = tecGeomArrowheadGetAngle(inputFileHandle, geom, &arrowheadAngle);
                int32_t arrowheadAttachment;
                i = tecGeomArrowheadGetAttach(inputFileHandle, geom, &arrowheadAttachment);
                double arrowheadSize;
                i = tecGeomArrowheadGetSize(inputFileHandle, geom, &arrowheadSize);
                int32_t arrowheadStyle;
                i = tecGeomArrowheadGetStyle(inputFileHandle, geom, &arrowheadStyle);
                if (arrowheadAttachment != 0)
                    i = tecGeomArrowheadSetInfo(outputFileHandle, arrowheadAngle, arrowheadAttachment, arrowheadSize, arrowheadStyle);
                break;
            case GeomType_Rectangle:
                double width, height;
                i = tecGeomRectangleGetSize(inputFileHandle, geom, &width, &height);
                i = tecGeomRectangleBegin(outputFileHandle, x, y, x + width, y + height, coordMode);
                break;
            case GeomType_Square:
                double squareSize;
                i = tecGeomSquareGetSize(inputFileHandle, geom, &squareSize);
                i = tecGeomSquareBegin(outputFileHandle, x, y, squareSize, coordMode);
                break;
            case GeomType_Circle:
                double radius;
                i = tecGeomEllipseGetNumPoints(inputFileHandle, geom, &numEllipsePoints);
                i = tecGeomCircleGetRadius(inputFileHandle, geom, &radius);
                i = tecGeomCircleBegin(outputFileHandle, x, y, radius, coordMode);
                i = tecGeomEllipseSetNumPoints(outputFileHandle, numEllipsePoints);
                break;
            case GeomType_Ellipse:
                i = tecGeomEllipseGetNumPoints(inputFileHandle, geom, &numEllipsePoints);
                double horizontalAxis, verticalAxis;
                i = tecGeomEllipseGetSize(inputFileHandle, geom, &horizontalAxis, &verticalAxis);
                i = tecGeomEllipseBegin(outputFileHandle, x, y, horizontalAxis, verticalAxis, coordMode);
                i = tecGeomEllipseSetNumPoints(outputFileHandle, numEllipsePoints);
                break;
            default:
                break;
            }

            // All of the below have sensible defaults, but we'll set all of it for completeness.
            int32_t clipping;
            i = tecGeomGetClipping(inputFileHandle, geom, &clipping);
            i = tecGeomSetClipping(outputFileHandle, clipping);

            int32_t isFilled;
            i = tecGeomIsFilled(inputFileHandle, geom, &isFilled);
            if (isFilled)
            {
                int32_t fillColor;
                i = tecGeomGetFillColor(inputFileHandle, geom, &fillColor);
                i = tecGeomFill(outputFileHandle, fillColor);
            }

            int32_t linePattern;
            i = tecGeomGetLinePattern(inputFileHandle, geom, &linePattern);
            double patternLength;
            i = tecGeomGetPatternLength(inputFileHandle, geom, &patternLength);
            double lineThickness;
            i = tecGeomGetLineThickness(inputFileHandle, geom, &lineThickness);
            int32_t color;
            i = tecGeomGetColor(inputFileHandle, geom, &color);
            i = tecGeomSetLineInfo(outputFileHandle, linePattern, patternLength, lineThickness, color);

            int32_t isAttached;
            i = tecGeomIsAttached(inputFileHandle, geom, &isAttached);
            if (isAttached)
            {
                int32_t attachZone;
                i = tecGeomGetZone(inputFileHandle, geom, &attachZone);
                i = tecGeomAttachToZone(outputFileHandle, attachZone);
            }

            char* macroFunctionCmd = NULL;
            i = tecGeomGetMacroFunctionCmd(inputFileHandle, geom, &macroFunctionCmd);
            if (strlen(macroFunctionCmd))
                i = tecGeomSetMacroFunctionCmd(outputFileHandle, macroFunctionCmd);
            tecStringFree(&macroFunctionCmd);

            int32_t scope;
            i = tecGeomGetScope(inputFileHandle, geom, &scope);
            i = tecGeomSetScope(outputFileHandle, scope);

            // Close the output geom
            i = tecGeomEnd(outputFileHandle);
        }

        // Texts
        int32_t numTexts;
        i = tecTextGetNumTexts(inputFileHandle, &numTexts);

        for (int32_t text = 1; text <= numTexts; ++text)
        {
            // Retrieve required information
            char* str = NULL;
            i = tecTextGetString(inputFileHandle, text, &str);

            double x, y, z;
            i = tecTextGetAnchorPos(inputFileHandle, text, &x, &y, &z);

            int32_t mode;
            i = tecTextGetCoordMode(inputFileHandle, text, &mode);
            CoordSys_e coordMode = (CoordSys_e)mode;

            double height;
            i = tecTextGetHeight(inputFileHandle, text, &height);

            int32_t sizeUnits;
            i = tecTextGetSizeUnits(inputFileHandle, text, &sizeUnits);

            // Begin a new text
            if (coordMode == CoordSys_Grid3D)
                i = tecText3DBegin(outputFileHandle, str, x, y, z, height, sizeUnits);
            else
                i = tecText2DBegin(outputFileHandle, str, x, y, mode, height, sizeUnits);
            tecStringFree(&str);

            // All of the below have sensible defaults, but we'll set everything here for completeness:
            int32_t color;
            i = tecTextGetColor(inputFileHandle, text, &color);
            i = tecTextSetColor(outputFileHandle, color);

            int32_t boxType;
            i = tecTextBoxGetType(inputFileHandle, text, &boxType);
            if (boxType != 0)
            {
                int32_t boxColor;
                i = tecTextBoxGetColor(inputFileHandle, text, &boxColor);
                int32_t boxFillColor;
                i = tecTextBoxGetFillColor(inputFileHandle, text, &boxFillColor);

                double boxLineThickness;
                i = tecTextBoxGetLineThickness(inputFileHandle, text, &boxLineThickness);
                double boxMargin;
                i = tecTextBoxGetMargin(inputFileHandle, text, &boxMargin);

                tecTextBoxSetInfo(outputFileHandle, boxType, boxColor, boxFillColor, boxLineThickness, boxMargin);
            }

            int32_t anchor;
            i = tecTextGetAnchor(inputFileHandle, text, &anchor);
            i = tecTextSetAnchor(outputFileHandle, anchor);

            int32_t isAttached;
            i = tecTextIsAttached(inputFileHandle, text, &isAttached);
            if (isAttached)
            {
                int32_t attachZone;
                i = tecTextGetZone(inputFileHandle, text, &attachZone);
                i = tecTextAttachToZone(outputFileHandle, attachZone);
            }

            double angle;
            i = tecTextGetAngle(inputFileHandle, text, &angle);
            i = tecTextSetAngle(outputFileHandle, angle);

            int32_t clipping;
            i = tecTextGetClipping(inputFileHandle, text, &clipping);
            i = tecTextSetClipping(outputFileHandle, clipping);

            int32_t scope;
            i = tecTextGetScope(inputFileHandle, text, &scope);
            i = tecTextSetScope(outputFileHandle, scope);

            char* typeface = NULL;
            i = tecTextGetTypeface(inputFileHandle, text, &typeface);
            int32_t isBold;
            i = tecTextIsBold(inputFileHandle, text, &isBold);
            int32_t isItalic;
            i = tecTextIsItalic(inputFileHandle, text, &isItalic);
            i = tecTextSetTypeface(outputFileHandle, typeface, isBold, isItalic);
            tecStringFree(&typeface);

            double lineSpacing;
            i = tecTextGetLineSpacing(inputFileHandle, text, &lineSpacing);
            i = tecTextSetLineSpacing(outputFileHandle, lineSpacing);

            char* macroFunctionCommand = NULL;
            i = tecTextGetMacroFunctionCmd(inputFileHandle, text, &macroFunctionCommand);
            if (strlen(macroFunctionCommand))
                i = tecTextSetMacroFunctionCmd(outputFileHandle, macroFunctionCommand);
            tecStringFree(&macroFunctionCommand);

            // Close the current text
            i = tecTextEnd(outputFileHandle);
        }

        // Close both old and new files
        i = tecFileReaderClose(&inputFileHandle);
        i = tecFileWriterClose(&outputFileHandle);
    }
    catch (std::runtime_error const& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

