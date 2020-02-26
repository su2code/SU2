!
!  Read a .szplt file provided as the first command-line argument,
!  and write a new .szplt file to the name provided as the second
!  command-line argument.
!
      ! Convenience to extract a returned C string to a FOrtran string
      subroutine copyCharArrayToString(charArray, length, string)
      use iso_c_binding, only : C_NULL_CHAR, c_ptr, c_f_pointer
      implicit none
      type(c_ptr) :: charArray
      integer length
      character(*) string
      
      character, pointer :: charPointer(:)
      integer i
      
      call c_f_pointer(charArray, charPointer, [length])
      
      string = ' '
      do i = 1, length
          string(i:i) = charPointer(i)
      enddo
      string(length+1:length+1) = C_NULL_CHAR
      return
      end
 
      ! Each tec function returns zero for success, non-zero for
      ! failure. We ignore these return codes for the sake of brevity.
      program main
      use iso_c_binding
      implicit none
      include "tecio.for"

      integer lastArgNum
      integer i, j
      character(256) programName
      character(256) inputFileName
      character(256) outputFileName
      character(256) dataSetTitle, zoneTitle
      character(256) name, val
      character(256) macroFunctionCmd
      character(256) textString, typeface, fontName, testFontName
      character(1024) varNames
      character(1024) labelSet
      character, pointer :: stringPtr(:)
      integer nameLen, strLen
      integer(c_int8_t), allocatable :: int8Values(:)
      integer(c_int16_t), allocatable :: int16Values(:)
      integer(c_int32_t) :: numVars, var
      integer(c_int32_t) :: fileType
      integer(c_int32_t) :: varType
      integer(c_int32_t) :: fileFormat = 1 ! .szplt
      integer(c_int32_t) :: outputDebugInfo = 1
      integer(c_int32_t) :: inputZone, numZones, outputZone
      integer(c_int32_t) :: zoneType
      integer(c_int32_t) :: strandID, parentZone, faceNeighborMode
      integer(c_int32_t) :: shareConnectivityFromZone
      integer(c_int32_t) :: zero = 0
      integer(c_int32_t) :: is64Bit
      integer(c_int32_t) :: numItems, whichItem
      integer(c_int32_t) :: numSets
      integer(c_int32_t) :: numGeoms, geom, geomType, clipping
      integer(c_int32_t) :: color, isFilled, fillColor, coordMode,
     &                      linePattern
      integer(c_int32_t) :: attachZone, isAttached, scope,
     &                      arrowheadAttachment
      integer(c_int32_t) :: arrowheadStyle, numSegments, segment, ind
      integer(c_int32_t) :: arrayLength
      integer(c_int32_t) :: numEllipsePoints = 2
      integer(c_int32_t) :: numTexts, text
      integer(c_int32_t) :: boxColor, boxFillColor, boxType, anchor
      integer(c_int32_t) :: isBold, isItalic, sizeUnits, font
      integer(c_int32_t), allocatable :: int32Values(:)
      integer(c_int32_t), allocatable :: faceConnections32(:)
      integer(c_int32_t), allocatable :: varTypes(:)
      integer(c_int32_t), allocatable :: passiveVarList(:)
      integer(c_int32_t), allocatable :: valueLocation(:)
      integer(c_int32_t), allocatable :: shareVarFromZone(:)
      integer(c_int32_t), allocatable :: nodeMap32(:)
      integer(c_int32_t), allocatable :: numSegPts(:)
      integer(c_int64_t) :: numValues, numFaceValues
      integer(c_int64_t) :: numFaceConnections
      integer(c_int64_t) :: iMax, jMax, kMax
      integer(c_int64_t), allocatable :: faceConnections64(:)
      integer(c_int64_t), allocatable :: nodeMap64(:)
      real(c_float), allocatable :: floatValues(:)
      real(c_double) :: solutionTime
      real(c_double) :: x, y, z, patternLength, lineThickness
      real(c_double) :: arrowheadAngle, arrowheadSize
      real(c_double) :: geomX, geomY, geomZ, width, height, squareSize
      real(c_double) :: radius, horizontalAxis, verticalAxis
      real(c_double) :: boxLineThickness, boxMargin, angle, lineSpacing
      real(c_double), allocatable :: doubleValues(:)
      real(c_double), allocatable :: xGeomData(:), yGeomData(:),
     &                               zGeomData(:)
      type(c_ptr) :: inputFileHandle = C_NULL_PTR
      type(c_ptr) :: outputFileHandle = C_NULL_PTR
      type(c_ptr) :: stringCPtr = C_NULL_PTR
      type(c_ptr) :: nameCPtr = C_NULL_PTR, valueCPtr = C_NULL_PTR
      
      ! Retrieve input and output file names
      lastArgNum = iargc()
      if (lastArgNum /= 2) then
          call getarg(0, programName)
        write(0,fmt='(a,a,a)')
     &    "Usage: ", trim(programName), " infilename outfilename"
        stop
      endif

      call getarg(1, inputFileName)
      inputFileName = trim(inputFileName) // C_NULL_CHAR
      call getarg(2, outputFileName)
      outputFileName = trim(outputFileName) // C_NULL_CHAR

      ! Open the input file for reading
      i = tecFileReaderOpen(inputFileName, inputFileHandle)

      ! Read info about the data set
      i = tecDataSetGetTitle(inputFileHandle, stringCPtr)
      call copyCharArrayToString(stringCPtr,
     &    tecStringLength(stringCPtr), dataSetTitle)
      call tecStringFree(stringCPtr)
      i = tecDataSetGetNumVars(inputFileHandle, numVars)

      strLen = 0
      do var = 1, numVars
          i = tecVarGetName(inputFileHandle, var, stringCPtr)
          nameLen = tecStringLength(stringCPtr)
          call c_f_pointer(stringCPtr, stringPtr, [nameLen])
          if (var .gt. 1) then
              strLen = strLen + 1
              varNames(strLen : strLen) = ','
          endif
          do j = 1, nameLen
              varNames(strLen + j : strLen + j) = stringPtr(j)
          enddo
          strLen = strLen + nameLen
          call tecStringFree(stringCPtr)
      enddo
      varNames(strLen + 1 : strlen + 1) = C_NULL_CHAR

      i = tecFileGetType(inputFileHandle, fileType)
      i = tecDataSetGetNumZones(inputFileHandle, numZones)

      ! Open the output file
      i = tecFileWriterOpen(outputFileName, dataSetTitle, varNames,
     &    fileFormat, fileType, 1, C_NULL_PTR, outputFileHandle)
      i = tecFileSetDiagnosticsLevel(outputFileHandle, outputDebugInfo)

      ! Zones
      do inputZone = 1, numZones
          i = tecZoneGetType(inputFileHandle, inputZone, zoneType)
          if (zoneType == 6 .or. zoneType == 7)
     &        stop "Unsupported inputZone type."
            
          ! Retrieve info about the inputZone
          i = tecZoneGetTitle(inputFileHandle, inputZone, stringCPtr)
          call copyCharArrayToString(stringCPtr,
     &        tecStringLength(stringCPtr), zoneTitle)
          call tecStringFree(stringCPtr)
          
          i = tecZoneGetIJK(inputFileHandle, inputZone,
     &        iMax, jMax, kMax)
     
          allocate(varTypes(numVars))
          allocate(passiveVarList(numVars))
          allocate(valueLocation(numVars))
          allocate(shareVarFromZone(numVars))
          do var = 1, numVars
              i = tecZoneVarGetType(inputFileHandle, inputZone,
     &            var, varTypes(var))
              i = tecZoneVarIsPassive(inputFileHandle, inputZone,
     &            var, passiveVarList(var))
              i = tecZoneVarGetValueLocation(inputFileHandle, inputZone,
     &            var, valueLocation(var))
              i = tecZoneVarGetSharedZone(inputFileHandle, inputZone,
     &            var, shareVarFromZone(var))
          enddo

          i = tecZoneConnectivityGetSharedZone(inputFileHandle,
     &        inputZone, shareConnectivityFromZone)
          i = tecZoneFaceNbrGetMode(inputFileHandle, inputZone,
     &        faceNeighborMode)
          if (faceNeighborMode > 4) faceNeighborMode = 1
          i = tecZoneFaceNbrGetNumConnections(inputfileHandle,
     &        inputZone, numFaceConnections)
     
          if (zoneType == 0) then
              i = tecZoneCreateIJK(outputFileHandle, zoneTitle,
     &            iMax, jMax, kMax, varTypes, shareVarFromZone,
     &            valueLocation, passiveVarList,
     &            shareConnectivityFromZone, numFaceConnections,
     &            faceNeighborMode, outputZone)
          else
              i = tecZoneCreateFE(outputFileHandle, zoneTitle, zoneType,
     &            iMax, jMax, varTypes, shareVarFromZone,
     &            valueLocation, passiveVarList,
     &            shareConnectivityFromZone, numFaceConnections,
     &            faceNeighborMode, outputZone)
          endif

          i = tecZoneGetSolutionTime(inputFileHandle, inputZone,
     &        solutionTime)
          i = tecZoneGetStrandID(inputFileHandle, inputZone, strandID)
          if (solutionTime /= 0.0 .or. strandID /= 0)
     &        i = tecZoneSetUnsteadyOptions(outputFileHandle,
     &            outputZone, solutionTime, strandID)

          i = tecZoneGetParentZone(inputFileHandle, inputZone,
     &        parentZone)
          if (parentZone /= 0)
     &        i = tecZoneSetParentZone(inputFileHandle, outputZone,
     &            parentZone)

          ! Read and write inputZone data
          do var = 1, numVars
              if (passiveVarList(var) == 0 .and.
     &              shareVarFromZone(var) == 0) then            
                  i = tecZoneVarGetNumValues(inputFileHandle,
     &                      inputZone, var, numValues)
                  select case (varTypes(var))
                  case (1) ! float
                      allocate(floatValues(numValues))
                      i = tecZoneVarGetFloatValues(inputFileHandle,
     &                    inputZone, var, 1_c_int64_t, numValues,
     &                    floatValues)
                      i = tecZoneVarWriteFloatValues(outputFileHandle,
     &                    outputZone, var, 0, numValues, floatValues)
                      deallocate(floatValues)
                  case (2) ! double
                      allocate(doubleValues(numValues))
                      i = tecZoneVarGetDoubleValues(inputFileHandle,
     &                    inputZone, var, 1_c_int64_t, numValues,
     &                    doubleValues)
                      i = tecZoneVarWriteDoubleValues(outputFileHandle,
     &                    outputZone, var, 0, numValues, doubleValues)
                      deallocate(doubleValues)
                  case (3) ! int32_t
                      allocate(int32Values(numValues))
                      i = tecZoneVarGetInt32Values(inputFileHandle,
     &                    inputZone, var, 1_c_int64_t, numValues,
     &                    int32Values)
                      i = tecZoneVarWriteInt32Values(outputFileHandle,
     &                    outputZone, var, 0, numValues, int32Values)
                      deallocate(int32Values)
                  case (4) ! int16_t
                      allocate(int16Values(numValues))
                      i = tecZoneVarGetInt16Values(inputFileHandle,
     &                    inputZone, var, 1_c_int64_t, numValues,
     &                    int16Values)
                      i = tecZoneVarWriteInt16Values(outputFileHandle,
     &                    outputZone, var, 0, numValues, int16Values)
                      deallocate(int16Values)
                  case (5) ! uint8_t
                      allocate(int8Values(numValues))
                      i = tecZoneVarGetUInt8Values(inputFileHandle,
     &                    inputZone, var, 1_c_int64_t, numValues,
     &                    int8Values)
                      i = tecZoneVarWriteUInt8Values(outputFileHandle,
     &                    outputZone, var, 0, numValues, int8Values)
                      deallocate(int8Values)
                  endselect
              endif
          enddo
          
          deallocate(varTypes)
          deallocate(passiveVarList)
          deallocate(valueLocation)
          deallocate(shareVarFromZone)

          ! Write zone face neighbors, if any
          if (numFaceConnections > 0) then
              i = tecZoneFaceNbrGetNumValues(inputFileHandle, inputZone,
     &            numFaceValues)
              i = tecZoneFaceNbrsAre64Bit(inputFileHandle, inputZone,
     &            is64Bit)
              if (is64Bit == 1) then
                  allocate(faceConnections64(numFaceValues))
                  i = tecZoneFaceNbrGetConnections64(inputFileHandle,
     &                inputZone, faceConnections64)
                  i = tecZoneFaceNbrWriteConnections64(outputFileHandle,
     &                outputZone, faceConnections64)
                  deallocate(faceConnections64)
              else
                  allocate(faceConnections32(numFaceValues))
                  i = tecZoneFaceNbrGetConnections(inputFileHandle,
     &                inputZone, faceConnections32)
                  i = tecZoneFaceNbrWriteConnections32(outputFileHandle,
     &                outputZone, faceConnections32)
                  deallocate(faceConnections32)
              endif
          endif

          ! Retrieve zone node map, if any, and send to the output file
          if (zoneType /= 0 .and.
     &        shareConnectivityFromZone == 0) then
              i = tecZoneNodeMapGetNumValues(inputFileHandle,
     &            inputZone, jMax, numValues)
              i = tecZoneNodeMapIs64Bit(inputFileHandle, inputZone,
     &            is64Bit)
              if (is64Bit == 1) then
                  allocate(nodeMap64(numValues))
                  i = tecZoneNodeMapGet64(inputFileHandle,
     &                inputZone, 1_8, jMax, nodeMap64)
                  i = tecZoneNodeMapWrite64(outputFileHandle,
     &                outputZone, 0, 1, numValues, nodeMap64)
                  deallocate(nodeMap64)
              else
                  allocate(nodeMap32(numValues))
                  i = tecZoneNodeMapGet(inputFileHandle,
     &                inputZone, 1_8, jMax, nodeMap32)
                  i = tecZoneNodeMapWrite32(outputFileHandle,
     &                outputZone, 0, 1, numValues, nodeMap32)
                  deallocate(nodeMap32)
              endif
          endif

          ! Zone aux data
          i = tecZoneAuxDataGetNumItems(inputFileHandle, inputZone,
     &        numItems)
          do whichItem = 1, numItems
              i = tecZoneAuxDataGetItem(inputFileHandle,
     &                inputZone, whichItem, nameCPtr, valueCPtr)
              call copyCharArrayToString(nameCPtr,
     &                 tecStringLength(nameCPtr), name)
              call copyCharArrayToString(valueCPtr,
     &                 tecStringLength(valueCPtr), val)
              i = tecZoneAddAuxData(outputFileHandle, outputZone,
     &            name, val) 
              call tecStringFree(nameCPtr)
              call tecStringFree(valueCPtr)
          enddo
                    
      enddo ! inputZone loop

      ! Custom label sets
      i = tecCustomLabelsGetNumSets(inputFileHandle, numSets)
      do j = 1, numSets
          i = tecCustomLabelsGetSet(inputFileHandle, j, stringCPtr) ! Labels returned as "\"Mon\",\"Tues\",\"Wed\""
          call copyCharArrayToString(stringCPtr,
     &         tecStringLength(stringCPtr), labelSet)
          i = tecCustomLabelsAddSet(outputFileHandle, labelSet)
          call tecStringFree(stringCPtr)
      enddo

      ! Data set aux data
      i = tecDataSetAuxDataGetNumItems(inputFileHandle, numItems)
      do j = 1, numItems
          i = tecDataSetAuxDataGetItem(inputFileHandle,
     &            j, nameCPtr, valueCPtr)
          call copyCharArrayToString(nameCPtr,
     &             tecStringLength(nameCPtr), name)
          call copyCharArrayToString(valueCPtr,
     &             tecStringLength(valueCPtr), val)
          i = tecDataSetAddAuxData(outputFileHandle, name, val)
          call tecStringFree(nameCPtr)
          call tecStringFree(valueCPtr)
      enddo

      ! Var aux data
      do var = 1, numVars
          i = tecVarAuxDataGetNumItems(inputFileHandle, var, numItems)
          do whichItem = 1, numItems
              i = tecVarAuxDataGetItem(inputFileHandle,
     &                var, whichItem, nameCPtr, valueCPtr)
              call copyCharArrayToString(nameCPtr,
     &                 tecStringLength(nameCPtr), name)
              call copyCharArrayToString(valueCPtr,
     &                 tecStringLength(valueCPtr), val)
              i = tecVarAddAuxData(outputFileHandle, var, name, val)
              call tecStringFree(nameCPtr);
              call tecStringFree(valueCPtr);
          enddo
      enddo

      ! Geometries
      i = tecGeomGetNumGeoms(inputFileHandle, numGeoms)
      do geom = 1, numGeoms
          i = tecGeomGetType(inputFileHandle, geom, geomType)
          i = tecGeomGetAnchorPos(inputFileHandle, geom, x, y, z)
          i = tecGeomGetCoordMode(inputFileHandle, geom, coordMode)

          select case (geomType)
          case (0, 5) ! GeomType_LineSegs, GeomType_LineSegs3D
              i = tecGeomLineGetSegmentCount(inputFileHandle,
     &            geom, numSegments)
              allocate(numSegPts(numSegments))
              arrayLength = 0;
              do segment = 1, numSegments
                  i = tecGeomLineSegmentGetPointCount(inputFileHandle,
     &                geom, segment, numSegPts(segment))
                  arrayLength = arrayLength + numSegPts(segment)
              enddo
              allocate(xGeomData(arrayLength))
              allocate(yGeomData(arrayLength))
              allocate(zGeomData(arrayLength))
              arrayLength = 0
              do segment = 1, numSegments
                  do ind = 1, numSegPts(segment)
                      arrayLength = arrayLength + 1
                      i = tecGeomLineGetPoint(inputFileHandle,
     &                    geom, segment, ind, geomX, geomY, geomZ)
                      xGeomData(arrayLength) = geomX
                      yGeomData(arrayLength) = geomY
                      zGeomData(arrayLength) = geomZ
                  enddo
              enddo
              if (geomType == 0) then
                  if (numSegments == 1) then
                      i = tecGeom2DLineSegmentsBegin(outputFileHandle,
     &                    x, y, numSegPts(1),
     &                    xGeomData, yGeomData, coordMode)
                  else
                      i = tecGeom2DMultiLineSegmentsBegin(
     &                    outputFileHandle, x, y, numSegments,
     &                    numSegPts, xGeomData, yGeomData, coordMode)
                  endif
              else
                  if (numSegments == 1) then
                      i = tecGeom3DLineSegmentsBegin(outputFileHandle,
     &                    x, y, z, numSegPts(1),
     &                    xGeomData, yGeomData, zGeomData)
                  else
                      i = tecGeom3DMultiLineSegmentsBegin(
     &                    outputFileHandle, x, y, z, numSegments,
     &                    numSegPts, xGeomData, yGeomData, zGeomData)
                  endif
              endif
              deallocate(xGeomData, yGeomData, zGeomData)
              deallocate(numSegPts)
              
              ! The default is no arrowheads
              i = tecGeomArrowheadGetAngle(inputFileHandle, geom,
     &            arrowheadAngle)
              i = tecGeomArrowheadGetAttach(inputFileHandle,
     &            geom, arrowheadAttachment)
              i = tecGeomArrowheadGetSize(inputFileHandle, geom,
     &            arrowheadSize)
              i = tecGeomArrowheadGetStyle(inputFileHandle, geom,
     &            arrowheadStyle)
              if (arrowheadAttachment /= 0)
     &            i = tecGeomArrowheadSetInfo(outputFileHandle,
     &                arrowheadAngle, arrowheadAttachment,
     &                arrowheadSize, arrowheadStyle)
          case (1) ! GeomType_Rectangle
              i = tecGeomRectangleGetSize(inputFileHandle,
     &            geom, width, height)
              i = tecGeomRectangleBegin(outputFileHandle, x, y,
     &            x + width, y + height, coordMode)
          case (2) ! GeomType_Square
              i = tecGeomSquareGetSize(inputFileHandle, geom,
     &            squareSize)
              i = tecGeomSquareBegin(outputFileHandle, x, y,
     &            squareSize, coordMode)
          case (3) ! GeomType_Circle
              i = tecGeomEllipseGetNumPoints(inputFileHandle,
     &            geom, numEllipsePoints)
              i = tecGeomCircleGetRadius(inputFileHandle, geom, radius)
              i = tecGeomCircleBegin(outputFileHandle, x, y,
     &            radius, coordMode)
              i = tecGeomEllipseSetNumPoints(outputFileHandle,
     &            numEllipsePoints)
          case (4) ! GeomType_Ellipse
              i = tecGeomEllipseGetNumPoints(inputFileHandle,
     &            geom, numEllipsePoints)
              i = tecGeomEllipseGetSize(inputFileHandle,
     &            geom, horizontalAxis, verticalAxis)
              i = tecGeomEllipseBegin(outputFileHandle, x, y,
     &            horizontalAxis, verticalAxis, coordMode)
              i = tecGeomEllipseSetNumPoints(outputFileHandle,
     &            numEllipsePoints)
          case default
              write(0,fmt='(a,a,a)') "Unsupported geometry type ",
     &            geomType, ". Skipping."
          end select
          
          i = tecGeomGetClipping(inputFileHandle, geom, clipping)
          i = tecGeomSetClipping(outputFileHandle, clipping)
          
          i = tecGeomIsFilled(inputFileHandle, geom, isFilled)
          if (isFilled == 1) then
              i = tecGeomGetFillColor(inputFileHandle, geom, fillColor)
              i = tecGeomFill(outputFileHandle, fillColor)
          endif
          
          i = tecGeomGetLinePattern(inputFileHandle, geom, linePattern)
          i = tecGeomGetPatternLength(inputFileHandle, geom,
     &        patternLength)
          i = tecGeomGetLineThickness(inputFileHandle, geom,
     &        lineThickness)
          i = tecGeomGetColor(inputFileHandle, geom, color)
          i = tecGeomSetLineInfo(outputFileHandle, linePattern,
     &        patternLength, lineThickness, color)
     
          i = tecGeomIsAttached(inputFileHandle, geom, isAttached)
          if (isAttached == 1) then
              i = tecGeomGetZone(inputFileHandle, geom, attachZone)
              i = tecGeomAttachToZone(inputFileHandle, attachZone)
          endif
          
          i = tecGeomGetMacroFunctionCmd(inputFileHandle, geom,
     &        stringCPtr)
          call copyCharArrayToString(stringCPtr,
     &         tecStringLength(stringCPtr), macroFunctionCmd)
          call tecStringFree(stringCPtr)
          i = tecGeomSetMacroFunctionCmd(outputFileHandle,
     &        macroFunctionCmd)
     
          i = tecGeomGetScope(inputFileHandle, geom, scope)
          i = tecGeomSetScope(outputFileHandle, scope)
          
          ! Close the output geom
          i = tecGeomEnd(outputFileHandle)

      enddo
      
      ! Texts
      i = tecTextGetNumTexts(inputFileHandle, numTexts)
      do text = 1, numTexts
          i = tecTextGetString(inputFileHandle, text, stringCPtr)
          call copyCharArrayToString(stringCPtr,
     &         tecStringLength(stringCPtr), textString)
          call tecStringFree(stringCPtr)
          i = tecTextGetAnchorPos(inputFileHandle, text, x, y, z)
          i = tecTextGetCoordMode(inputFileHandle, text, coordMode)
          i = tecTextGetHeight(inputFileHandle, text, height)
          i = tecTextGetSizeUnits(inputFileHandle, text, sizeUnits)
          
          ! Begin a new output text
          if (coordMode == 6) then ! 3D
              i = tecText3DBegin(outputFileHandle, textString, x, y, z,
     &            height, sizeUnits)
          else
              i = tecText2DBegin(outputFileHandle, textString, x, y,
     &            coordMode, height, sizeUnits)
          endif
          
          ! All of the below have sensible defaults, but we'll set everything here for completeness
          i = tecTextGetColor(inputFileHandle, text, color)
          i = tecTextSetColor(outputFileHandle, color)
          
          i = tecTextBoxGetType(inputFileHandle, text, boxType)
          if (boxType /= 0) then
              i = tecTextBoxGetColor(inputFileHandle, text, boxColor)
              i = tecTextBoxGetFillColor(inputFileHandle, text,
     &            boxFillColor)
              i = tecTextBoxGetLineThickness(inputFileHandle,
     &            text, boxLineThickness)
              i = tecTextBoxGetMargin(inputFileHandle, text, boxMargin)
              
              i = tecTextBoxSetInfo(outputFileHandle, boxType, boxColor,
     &            boxFillColor, boxLineThickness, boxMargin)
          endif
          
          i = tecTextGetAnchor(inputFileHandle, text, anchor)
          i = tecTextSetAnchor(outputFileHandle, anchor)
          
          i = tecTextIsAttached(inputFileHandle, text, isAttached)
          if (isAttached == 1) then
              i = tecTextGetZone(inputFileHandle, text, attachZone)
              i = tecTextAttachToZone(outputFileHandle, attachZone)
          endif
          
          i = tecTextGetAngle(inputFileHandle, text, angle)
          i = tecTextSetAngle(outputFileHandle, angle)
          
          i = tecTextGetClipping(inputFileHandle, text, clipping)
          i = tecTextSetClipping(outputFileHandle, clipping)
          
          i = tecTextGetScope(inputFileHandle, text, scope)
          i = tecTextSetScope(outputFileHandle, scope)
          
          i = tecTextGetTypeface(inputFileHandle, text, stringCPtr)
          call copyCharArrayToString(stringCPtr,
     &             tecStringLength(stringCPtr), typeface)
          call tecStringFree(stringCPtr)
          i = tecTextIsBold(inputFileHandle, text, isBold)
          i = tecTextIsItalic(inputFileHandle, text, isItalic)
          i = tecTextSetTypeface(outputFileHandle, typeface,
     &        isBold, isItalic)

          i = tecTextGetLineSpacing(inputFileHandle, text, lineSpacing)
          i = tecTextSetLineSpacing(outputFileHandle, lineSpacing)
          
          i = tecTextGetMacroFunctionCmd(inputFileHandle,
     &            text, stringCPtr)
          call copyCharArrayToString(stringCPtr,
     &             tecStringLength(stringCPtr), macroFunctionCmd)
          call tecStringFree(stringCPtr)
          i = tecTextSetMacroFunctionCmd(outputFileHandle,
     &        macroFunctionCmd)
     
          ! Close the output text
          i = tecTextEnd(outputFileHandle)

      enddo

      ! Close old and new files
      i = tecFileWriterClose(outputFileHandle)
      i = tecFileReaderClose(inputFileHandle)

      end
