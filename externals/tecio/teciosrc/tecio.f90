
interface

! Output routine interfaces
integer(c_int32_t) function tecFileWriterOpen( &
    fileName, &
    dataSetTitle, &
    variableList, &
    fileFormat, &
    fileType, &
    defaultVarType, &
    gridFileHandle, &
    fileHandle) &
    bind(c, name="tecFileWriterOpen")
    use iso_c_binding
    implicit none
    character(c_char), intent(in)         :: fileName(*)
    character(c_char), intent(in)         :: dataSetTitle(*)
    character(c_char), intent(in)         :: variableList(*)
    integer(c_int32_t), value, intent(in) :: fileFormat
    integer(c_int32_t), value, intent(in) :: fileType
    integer(c_int32_t), value, intent(in) :: defaultVarType
    type(c_ptr), value, intent(in)        :: gridFileHandle
    type(c_ptr), intent(out)              :: fileHandle
end function tecFileWriterOpen

integer(c_int32_t) function tecFileSetDiagnosticsLevel( &
    fileHandle, &
    level) &
    bind(c, name="tecFileSetDiagnosticsLevel")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: level
end function tecFileSetDiagnosticsLevel

integer(c_int32_t) function tecMPIInitialize( &
    fileHandle, &
    communicator, &
    mainRank) &
    bind(c, name="tecMPIInitializef")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: communicator
    integer(c_int32_t), value, intent(in) :: mainRank
end function tecMPIInitialize

integer(c_int32_t) function tecZoneCreateIJK( &
    fileHandle, &
    zoneTitle, &
    imax, &
    jmax, &
    kmax, &
    varTypes, &
    shareVarFromZone, &
    valueLocations, &
    passiveVarList, &
    shareFaceNeighborsFromZone, &
    numFaceConnections, &
    faceNeighborMode, &
    zone) &
    bind(c, name="tecZoneCreateIJK")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)         :: fileHandle
    character(c_char), intent(in)          :: zoneTitle(*)
    integer(c_int64_t), value, intent(in)  :: imax
    integer(c_int64_t), value, intent(in)  :: jmax
    integer(c_int64_t), value, intent(in)  :: kmax
    integer(c_int32_t), intent(in)         :: varTypes(*)
    integer(c_int32_t), intent(in)         :: shareVarFromZone(*)
    integer(c_int32_t), intent(in)         :: valueLocations(*)
    integer(c_int32_t), intent(in)         :: passiveVarList(*)
    integer(c_int32_t), value, intent(in)  :: shareFaceNeighborsFromZone
    integer(c_int64_t), value, intent(in)  :: numFaceConnections
    integer(c_int32_t), value, intent(in)  :: faceNeighborMode
    integer(c_int32_t), intent(out)        :: zone
end function tecZoneCreateIJK

integer(c_int32_t) function tecZoneCreateFE( &
    fileHandle, &
    zoneTitle, &
    zoneType, &
    numNodes, &
    numCells, &
    varTypes, &
    shareVarFromZone, &
    valueLocations, &
    passiveVarList, &
    shareConnectivityFromZone, &
    numFaceConnections, &
    faceNeighborMode, &
    zone) &
    bind(c, name="tecZoneCreateFE")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)         :: fileHandle
    character(c_char), intent(in)          :: zoneTitle(*)
    integer(c_int32_t), value, intent(in)  :: zoneType
    integer(c_int64_t), value, intent(in)  :: numNodes
    integer(c_int64_t), value, intent(in)  :: numCells
    integer(c_int32_t), intent(in)         :: varTypes(*)
    integer(c_int32_t), intent(in)         :: shareVarFromZone(*)
    integer(c_int32_t), intent(in)         :: valueLocations(*)
    integer(c_int32_t), intent(in)         :: passiveVarList(*)
    integer(c_int32_t), value, intent(in)  :: shareConnectivityFromZone
    integer(c_int64_t), value, intent(in)  :: numFaceConnections
    integer(c_int32_t), value, intent(in)  :: faceNeighborMode
    integer(c_int32_t), intent(out)        :: zone
end function tecZoneCreateFE

integer(c_int32_t) function tecZoneCreatePoly( &
    fileHandle, &
    zoneTitle, &
    zoneType, &
    numNodes, &
    numFaces, &
    numCells, &
    totalNumFaceNodes, &
    varTypes, &
    shareVarFromZone, &
    valueLocations, &
    passiveVarList, &
    shareConnectivityFromZone, &
    numConnectedBoundaryFaces, &
    totalNumBoundaryConnections, &
    zone) &
    bind(c, name="tecZoneCreatePoly")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)         :: fileHandle
    character(c_char), intent(in)          :: zoneTitle(*)
    integer(c_int32_t), value, intent(in)  :: zoneType
    integer(c_int64_t), value, intent(in)  :: numNodes
    integer(c_int64_t), value, intent(in)  :: numFaces
    integer(c_int64_t), value, intent(in)  :: numCells
    integer(c_int64_t), value, intent(in)  :: totalNumFaceNodes
    integer(c_int32_t), intent(in)         :: varTypes(*)
    integer(c_int32_t), intent(in)         :: shareVarFromZone(*)
    integer(c_int32_t), intent(in)         :: valueLocations(*)
    integer(c_int32_t), intent(in)         :: passiveVarList(*)
    integer(c_int32_t), value, intent(in)  :: shareConnectivityFromZone
    integer(c_int64_t), value, intent(in)  :: numConnectedBoundaryFaces
    integer(c_int64_t), value, intent(in)  :: totalNumBoundaryConnections
    integer(c_int32_t), intent(out)        :: zone
end function tecZoneCreatePoly

integer(c_int32_t) function tecZoneSetUnsteadyOptions( &
    fileHandle, &
    zone, &
    solutionTime, &
    strandID) &
    bind(c, name="tecZoneSetUnsteadyOptions")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    real(c_double), value, intent(in)     :: solutionTime
    integer(c_int32_t), value, intent(in) :: strandID
end function tecZoneSetUnsteadyOptions

integer(c_int32_t) function tecZoneSetParentZone( &
    fileHandle, &
    zone, &
    parentZone) &
    bind(c, name="tecZoneSetParentZone")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: parentZone
end function tecZoneSetParentZone

integer(c_int32_t) function tecZoneMapPartitionsToMPIRanks( &
    fileHandle, &
    zone, &
    numPartitions, &
    mpiRanksForPartitions) &
    bind(c, name="tecZoneMapPartitionsToMPIRanks")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), intent(inout)     :: zone
    integer(c_int32_t), value, intent(in) :: numPartitions
    integer(c_int32_t), intent(in)        :: mpiRanksForPartitions(*)
    end function tecZoneMapPartitionsToMPIRanks

integer(c_int32_t) function tecFEPartitionCreate32( &
    fileHandle, &
    zone, &
    partition, &
    numNodes, &
    numCells, &
    numGhostNodes, &
    ghostNodes, &
    neighborPartitions, &
    neighborPartitionNodes, &
    numGhostCells, &
    ghostCells) &
    bind(c, name="tecFEPartitionCreate32")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int64_t), value, intent(in) :: numNodes
    integer(c_int64_t), value, intent(in) :: numCells
    integer(c_int64_t), value, intent(in) :: numGhostNodes
    integer(c_int32_t), intent(in)        :: ghostNodes(*)
    integer(c_int32_t), intent(in)        :: neighborPartitions(*)
    integer(c_int32_t), intent(in)        :: neighborPartitionNodes(*)
    integer(c_int64_t), value, intent(in) :: numGhostCells
    integer(c_int32_t), intent(in)        :: ghostCells(*)
end function tecFEPartitionCreate32

integer(c_int32_t) function tecFEPartitionCreate64( &
    fileHandle, &
    zone, &
    partition, &
    numNodes, &
    numCells, &
    numGhostNodes, &
    ghostNodes, &
    neighborPartitions, &
    neighborPartitionNodes, &
    numGhostCells, &
    ghostCells) &
    bind(c, name="tecFEPartitionCreate64")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int64_t), value, intent(in) :: numNodes
    integer(c_int64_t), value, intent(in) :: numCells
    integer(c_int64_t), value, intent(in) :: numGhostNodes
    integer(c_int64_t), intent(in)        :: ghostNodes(*)
    integer(c_int32_t), intent(in)        :: neighborPartitions(*)
    integer(c_int64_t), intent(in)        :: neighborPartitionNodes(*)
    integer(c_int64_t), value, intent(in) :: numGhostCells
    integer(c_int64_t), intent(in)        :: ghostCells(*)
end function tecFEPartitionCreate64

integer(c_int32_t) function tecIJKPartitionCreate64( &
    fileHandle, &
    zone, &
    partition, &
    imin, &
    jmin, &
    kmin, &
    imax, &
    jmax, &
    kmax) &
    bind(c, name="tecIJKPartitionCreate64")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int64_t), value, intent(in) :: imin
    integer(c_int64_t), value, intent(in) :: jmin
    integer(c_int64_t), value, intent(in) :: kmin
    integer(c_int64_t), value, intent(in) :: imax
    integer(c_int64_t), value, intent(in) :: jmax
    integer(c_int64_t), value, intent(in) :: kmax
end function tecIJKPartitionCreate64

integer(c_int32_t) function tecZoneVarWriteDoubleValues( &
    fileHandle, &
    zone, &
    var, &
    partition, &
    count, &
    values) &
    bind(c, name="tecZoneVarWriteDoubleValues")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: var
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int64_t), value, intent(in) :: count
    real(c_double), intent(in)            :: values(*)
end function tecZoneVarWriteDoubleValues

integer(c_int32_t) function tecZoneVarWriteFloatValues( &
    fileHandle, &
    zone, &
    var, &
    partition, &
    count, &
    values) &
    bind(c, name="tecZoneVarWriteFloatValues")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: var
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int64_t), value, intent(in) :: count
    real(c_float), intent(in)             :: values(*)
end function tecZoneVarWriteFloatValues

integer(c_int32_t) function tecZoneVarWriteInt32Values( &
    fileHandle, &
    zone, &
    var, &
    partition, &
    count, &
    values) &
    bind(c, name="tecZoneVarWriteInt32Values")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: var
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int64_t), value, intent(in) :: count
    integer(c_int32_t), intent(in)        :: values(*)
end function tecZoneVarWriteInt32Values

integer(c_int32_t) function tecZoneVarWriteInt16Values( &
    fileHandle, &
    zone, &
    var, &
    partition, &
    count, &
    values) &
    bind(c, name="tecZoneVarWriteInt16Values")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: var
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int64_t), value, intent(in) :: count
    integer(c_int16_t), intent(in)        :: values(*)
end function tecZoneVarWriteInt16Values

integer(c_int32_t) function tecZoneVarWriteUInt8Values( &
    fileHandle, &
    zone, &
    var, &
    partition, &
    count, &
    values) &
    bind(c, name="tecZoneVarWriteUInt8Values")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: var
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int64_t), value, intent(in) :: count
    integer(c_int8_t), intent(in)         :: values(*)
end function tecZoneVarWriteUInt8Values

integer(c_int32_t) function tecZoneNodeMapWrite32( &
    fileHandle, &
    zone, &
    partition, &
    nodesAreOneBased, &
    count, &
    nodes) &
    bind(c, name="tecZoneNodeMapWrite32")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int32_t), value, intent(in) :: nodesAreOneBased
    integer(c_int64_t), value, intent(in) :: count
    integer(c_int32_t), intent(in)        :: nodes(*)
end function tecZoneNodeMapWrite32

integer(c_int32_t) function tecZoneNodeMapWrite64( &
    fileHandle, &
    zone, &
    partition, &
    nodeAreOneBased, &
    count, &
    nodes) &
    bind(c, name="tecZoneNodeMapWrite64")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int32_t), value, intent(in) :: nodeAreOneBased
    integer(c_int64_t), value, intent(in) :: count
    integer(c_int64_t), intent(in)        :: nodes(*)
end function tecZoneNodeMapWrite64

integer(c_int32_t) function tecZoneFaceNbrWriteConnections32( &
    fileHandle, &
    zone, &
    faceNeighbors) &
    bind(c, name="tecZoneFaceNbrWriteConnections32")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), intent(in)        :: faceNeighbors(*)
end function tecZoneFaceNbrWriteConnections32

integer(c_int32_t) function tecZoneFaceNbrWriteConnections64( &
    fileHandle, &
    zone, &
    faceNeighbors) &
    bind(c, name="tecZoneFaceNbrWriteConnections64")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int64_t), intent(in)        :: faceNeighbors(*)
end function tecZoneFaceNbrWriteConnections64

integer(c_int32_t) function tecZoneWritePolyFaces32( &
    fileHandle, &
    zone, &
    partition, &
    numFaces, &
    faceNodeCounts, &
    faceNodes, &
    faceLeftElems, &
    faceRightElems, &
    isOneBased) &
    bind(c, name="tecZoneWritePolyFaces32")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int32_t), value, intent(in) :: numFaces
    integer(c_int32_t), intent(in)        :: faceNodeCounts(*)
    integer(c_int32_t), intent(in)        :: faceNodes(*)
    integer(c_int32_t), intent(in)        :: faceLeftElems(*)
    integer(c_int32_t), intent(in)        :: faceRightElems(*)
    integer(c_int32_t), value, intent(in) :: isOneBased
end function tecZoneWritePolyFaces32

integer(c_int32_t) function tecZoneWritePolyFaces64( &
    fileHandle, &
    zone, &
    partition, &
    numFaces, &
    faceNodeCounts, &
    faceNodes, &
    faceLeftElems, &
    faceRightElems, &
    isOneBased) &
    bind(c, name="tecZoneWritePolyFaces64")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int64_t), value, intent(in) :: numFaces
    integer(c_int32_t), intent(in)        :: faceNodeCounts(*)
    integer(c_int64_t), intent(in)        :: faceNodes(*)
    integer(c_int64_t), intent(in)        :: faceLeftElems(*)
    integer(c_int64_t), intent(in)        :: faceRightElems(*)
    integer(c_int32_t), value, intent(in) :: isOneBased
end function tecZoneWritePolyFaces64

integer(c_int32_t) function tecZoneWritePolyBoundaryConnections32( &
    fileHandle, &
    zone, &
    partition, &
    numBoundaryFaces, &
    faceBoundaryConnectionCounts, &
    faceBoundaryConnectionElems, &
    faceBoundaryConnectionZones, &
    isOneBased) &
    bind(c, name="tecZoneWritePolyBoundaryConnections32")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int32_t), value, intent(in) :: numBoundaryFaces
    integer(c_int32_t), intent(in)        :: faceBoundaryConnectionCounts(*)
    integer(c_int32_t), intent(in)        :: faceBoundaryConnectionElems(*)
    integer(c_int32_t), intent(in)        :: faceBoundaryConnectionZones(*)
    integer(c_int32_t), value, intent(in) :: isOneBased
    end function tecZoneWritePolyBoundaryConnections32

integer(c_int32_t) function tecZoneWritePolyBoundaryConnections64( &
    fileHandle, &
    zone, &
    partition, &
    numBoundaryFaces, &
    faceBoundaryConnectionCounts, &
    faceBoundaryConnectionElems, &
    faceBoundaryConnectionZones, &
    isOneBased) &
    bind(c, name="tecZoneWritePolyBoundaryConnections64")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    integer(c_int32_t), value, intent(in) :: partition
    integer(c_int64_t), value, intent(in) :: numBoundaryFaces
    integer(c_int32_t), intent(in)        :: faceBoundaryConnectionCounts(*)
    integer(c_int64_t), intent(in)        :: faceBoundaryConnectionElems(*)
    integer(c_int32_t), intent(in)        :: faceBoundaryConnectionZones(*)
    integer(c_int32_t), value, intent(in) :: isOneBased
end function tecZoneWritePolyBoundaryConnections64

integer(c_int32_t) function tecDataSetAddAuxData( &
    fileHandle, &
    name, &
    value) &
    bind(c, name="tecDataSetAddAuxData")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)       :: fileHandle
    character(c_char), intent(in)        :: name(*)
    character(c_char), intent(in)        :: value(*)
end function tecDataSetAddAuxData

integer(c_int32_t) function tecVarAddAuxData( &
    fileHandle, &
    var, &
    name, &
    value) &
    bind(c, name="tecVarAddAuxData")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: var
    character(c_char), intent(in)         :: name(*)
    character(c_char), intent(in)         :: value(*)
end function tecVarAddAuxData

integer(c_int32_t) function tecZoneAddAuxData( &
    fileHandle, &
    zone, &
    name, &
    value) &
    bind(c, name="tecZoneAddAuxData")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
    character(c_char), intent(in)         :: name(*)
    character(c_char), intent(in)         :: value(*)
end function tecZoneAddAuxData

integer(c_int32_t) function tecGeom2DLineSegmentsBegin( &
    fileHandle, &
    xOrigin, &
    yOrigin, &
    numPoints, &
    relativeX, &
    relativeY, &
    posCoordMode) &
    bind(c, name="tecGeom2DLineSegmentsBegin")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    real(c_double), value, intent(in)     :: xOrigin
    real(c_double), value, intent(in)     :: yOrigin
    integer(c_int32_t), value, intent(in) :: numPoints
    real(c_double), intent(in)            :: relativeX(*)
    real(c_double), intent(in)            :: relativeY(*)
    integer(c_int32_t), value, intent(in) :: posCoordMode
end function tecGeom2DLineSegmentsBegin

integer(c_int32_t) function tecGeom2DMultiLineSegmentsBegin( &
    fileHandle, &
    xOrigin, &
    yOrigin, &
    numSegments, &
    numSegmentPoints, &
    relativeX, &
    relativeY, &
    posCoordMode) &
    bind(c, name="tecGeom2DMultiLineSegmentsBegin")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    real(c_double), value, intent(in)     :: xOrigin
    real(c_double), value, intent(in)     :: yOrigin
    integer(c_int32_t), value, intent(in) :: numSegments
    integer(c_int32_t), intent(in)        :: numSegmentPoints(*)
    real(c_double), intent(in)            :: relativeX(*)
    real(c_double), intent(in)            :: relativeY(*)
    integer(c_int32_t), value, intent(in) :: posCoordMode
end function tecGeom2DMultiLineSegmentsBegin

integer(c_int32_t) function tecGeom3DLineSegmentsBegin( &
    fileHandle, &
    xOrigin, &
    yOrigin, &
    zOrigin, &
    numPoints, &
    relativeX, &
    relativeY, &
    relativeZ) &
    bind(c, name="tecGeom3DLineSegmentsBegin")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    real(c_double), value, intent(in)     :: xOrigin
    real(c_double), value, intent(in)     :: yOrigin
    real(c_double), value, intent(in)     :: zOrigin
    integer(c_int32_t), value, intent(in) :: numPoints
    real(c_double), intent(in)            :: relativeX(*)
    real(c_double), intent(in)            :: relativeY(*)
    real(c_double), intent(in)            :: relativeZ(*)
end function tecGeom3DLineSegmentsBegin

integer(c_int32_t) function tecGeom3DMultiLineSegmentsBegin( &
    fileHandle, &
    xOrigin, &
    yOrigin, &
    zOrigin, &
    numSegments, &
    numSegmentPoints, &
    relativeX, &
    relativeY, &
    relativeZ) &
    bind(c, name="tecGeom3DMultiLineSegmentsBegin")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    real(c_double), value, intent(in)     :: xOrigin
    real(c_double), value, intent(in)     :: yOrigin
    real(c_double), value, intent(in)     :: zOrigin
    integer(c_int32_t), value, intent(in) :: numSegments
    integer(c_int32_t), intent(in)        :: numSegmentPoints(*)
    real(c_double), intent(in)            :: relativeX(*)
    real(c_double), intent(in)            :: relativeY(*)
    real(c_double), intent(in)            :: relativeZ(*)
end function tecGeom3DMultiLineSegmentsBegin

integer(c_int32_t) function tecGeomCircleBegin( &
    fileHandle, &
    xCenter, &
    yCenter, &
    radius, &
    posCoordMode) &
    bind(c, name="tecGeomCircleBegin")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    real(c_double), value, intent(in)     :: xCenter
    real(c_double), value, intent(in)     :: yCenter
    real(c_double), value, intent(in)     :: radius
    integer(c_int32_t), value, intent(in) :: posCoordMode
end function tecGeomCircleBegin

integer(c_int32_t) function tecGeomEllipseBegin( &
    fileHandle, &
    xCenter, &
    yCenter, &
    width, &
    height, &
    posCoordMode) &
    bind(c, name="tecGeomEllipseBegin")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    real(c_double), value, intent(in)     :: xCenter
    real(c_double), value, intent(in)     :: yCenter
    real(c_double), value, intent(in)     :: width
    real(c_double), value, intent(in)     :: height
    integer(c_int32_t), value, intent(in) :: posCoordMode
end function tecGeomEllipseBegin

integer(c_int32_t) function tecGeomRectangleBegin( &
    fileHandle, &
    xMin, &
    yMin, &
    xMax, &
    yMax, &
    posCoordMode) &
    bind(c, name="tecGeomRectangleBegin")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    real(c_double), value, intent(in)     :: xMin
    real(c_double), value, intent(in)     :: yMin
    real(c_double), value, intent(in)     :: xMax
    real(c_double), value, intent(in)     :: yMax
    integer(c_int32_t), value, intent(in) :: posCoordMode
end function tecGeomRectangleBegin

integer(c_int32_t) function tecGeomSquareBegin( &
    fileHandle, &
    xMin, &
    yMin, &
    size, &
    posCoordMode) &
    bind(c, name="tecGeomSquareBegin")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    real(c_double), value, intent(in)     :: xMin
    real(c_double), value, intent(in)     :: yMin
    real(c_double), value, intent(in)     :: size
    integer(c_int32_t), value, intent(in) :: posCoordMode
end function tecGeomSquareBegin

integer(c_int32_t) function tecGeomArrowheadSetInfo( &
    fileHandle, &
    angle, &
    attachment, &
    size, &
    style) &
    bind(c, name="tecGeomArrowheadSetInfo")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    real(c_double), value, intent(in)     :: angle
    integer(c_int32_t), value, intent(in) :: attachment
    real(c_double), value, intent(in)     :: size
    integer(c_int32_t), value, intent(in) :: style
end function tecGeomArrowheadSetInfo

integer(c_int32_t) function tecGeomEllipseSetNumPoints( &
    fileHandle, &
    numEllipsePoints) &
    bind(c, name="tecGeomEllipseSetNumPoints")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: numEllipsePoints
end function tecGeomEllipseSetNumPoints

integer(c_int32_t) function tecGeomSetClipping( &
    fileHandle, &
    clipping) &
    bind(c, name="tecGeomSetClipping")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: clipping
end function tecGeomSetClipping

integer(c_int32_t) function tecGeomSetLineInfo( &
    fileHandle, &
    linePattern, &
    patternLength, &
    thickness, &
    color) &
    bind(c, name="tecGeomSetLineInfo")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: linePattern
    real(c_double), value, intent(in)     :: patternLength
    real(c_double), value, intent(in)     :: thickness
    integer(c_int32_t), value, intent(in) :: color
end function tecGeomSetLineInfo

integer(c_int32_t) function tecGeomSetMacroFunctionCmd( &
    fileHandle, &
    macroFunctionCmd) &
    bind(c, name="tecGeomSetMacroFunctionCmd")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in) :: fileHandle
    character(c_char), intent(in)  :: macroFunctionCmd(*)
end function tecGeomSetMacroFunctionCmd

integer(c_int32_t) function tecGeomSetScope( &
    fileHandle, &
    scope) &
    bind(c, name="tecGeomSetScope")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: scope
end function tecGeomSetScope

integer(c_int32_t) function tecGeomAttachToZone( &
    fileHandle, &
    zone) &
    bind(c, name="tecGeomAttachToZone")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
end function tecGeomAttachToZone

integer(c_int32_t) function tecGeomFill( &
    fileHandle, &
    fillColor) &
    bind(c, name="tecGeomFill")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: fillColor
    end function tecGeomFill

integer(c_int32_t) function tecGeomEnd( &
    fileHandle) &
    bind(c, name="tecGeomEnd")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in) :: fileHandle
end function tecGeomEnd

integer(c_int32_t) function tecCustomLabelsAddSet( &
    fileHandle, &
    labels) &
    bind(c, name="tecCustomLabelsAddSet")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in) :: fileHandle
    character(c_char), intent(in)  :: labels(*)
end function tecCustomLabelsAddSet

integer(c_int32_t) function tecText2DBegin( &
    fileHandle, &
    string, &
    x, &
    y, &
    posCoordMode, &
    height, &
    sizeUnits) &
    bind(c, name="tecText2DBegin")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    character(c_char), intent(in)         :: string(*)
    real(c_double), value, intent(in)     :: x
    real(c_double), value, intent(in)     :: y
    integer(c_int32_t), value, intent(in) :: posCoordMode
    real(c_double), value, intent(in)     :: height
    integer(c_int32_t), value, intent(in) :: sizeUnits
end function tecText2DBegin

integer(c_int32_t) function tecText3DBegin( &
    fileHandle, &
    string, &
    x, &
    y, &
    z, &
    height, &
    sizeUnits) &
    bind(c, name="tecText3DBegin")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    character(c_char), intent(in)         :: string(*)
    real(c_double), value, intent(in)     :: x
    real(c_double), value, intent(in)     :: y
    real(c_double), value, intent(in)     :: z
    real(c_double), value, intent(in)     :: height
    integer(c_int32_t), value, intent(in) :: sizeUnits
end function tecText3DBegin

integer(c_int32_t) function tecTextAttachToZone( &
    fileHandle, &
    zone) &
    bind(c, name="tecTextAttachToZone")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: zone
end function tecTextAttachToZone

integer(c_int32_t) function tecTextBoxSetInfo( &
    fileHandle, &
    boxType, &
    lineColor, &
    fillColor, &
    lineThickness, &
    margin) &
    bind(c, name="tecTextBoxSetInfo")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: boxType
    integer(c_int32_t), value, intent(in) :: lineColor
    integer(c_int32_t), value, intent(in) :: fillColor
    real(c_double), value, intent(in)     :: lineThickness
    real(c_double), value, intent(in)     :: margin
end function tecTextBoxSetInfo

integer(c_int32_t) function tecTextSetAnchor( &
    fileHandle, &
    anchor) &
    bind(c, name="tecTextSetAnchor")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: anchor
end function tecTextSetAnchor

integer(c_int32_t) function tecTextSetAngle( &
    fileHandle, &
    angle) &
    bind(c, name="tecTextSetAngle")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)    :: fileHandle
    real(c_double), value, intent(in) :: angle
end function tecTextSetAngle

integer(c_int32_t) function tecTextSetClipping( &
    fileHandle, &
    clipping) &
    bind(c, name="tecTextSetClipping")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: clipping
end function tecTextSetClipping

integer(c_int32_t) function tecTextSetColor( &
    fileHandle, &
    color) &
    bind(c, name="tecTextSetColor")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: color
end function tecTextSetColor

integer(c_int32_t) function tecTextSetTypeface( &
    fileHandle, &
    family, &
    isBold, &
    isItalic) &
    bind(c, name="tecTextSetTypeface")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    character(c_char), intent(in)         :: family(*)
    integer(c_int32_t), value, intent(in) :: isBold
    integer(c_int32_t), value, intent(in) :: isItalic
end function tecTextSetTypeface

integer(c_int32_t) function tecTextSetLineSpacing( &
    fileHandle, &
    lineSpacing) &
    bind(c, name="tecTextSetLineSpacing")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)    :: fileHandle
    real(c_double), value, intent(in) :: lineSpacing
end function tecTextSetLineSpacing

integer(c_int32_t) function tecTextSetMacroFunctionCmd( &
    fileHandle, &
    macroFunctionCmd) &
    bind(c, name="tecTextSetMacroFunctionCmd")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in) :: fileHandle
    character(c_char), intent(in)  :: macroFunctionCmd(*)
end function tecTextSetMacroFunctionCmd

integer(c_int32_t) function tecTextSetScope( &
    fileHandle, &
    scope) &
    bind(c, name="tecTextSetScope")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: scope
end function tecTextSetScope

integer(c_int32_t) function tecTextEnd( &
    fileHandle) &
    bind(c, name="tecTextEnd")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in) :: fileHandle
end function tecTextEnd

integer(c_int32_t) function tecUserRecAdd( &
    fileHandle, &
    userRec) &
    bind(c, name="tecUserRecAdd")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in) :: fileHandle
    character(c_char), intent(in)  :: userRec(*)
end function tecUserRecAdd
    
integer(c_int32_t) function tecFileWriterFlush( &
    fileHandle, &
    numZonesToRetain, &
    zonesToRetain) &
    bind(c, name="tecFileWriterFlush")
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in)        :: fileHandle
    integer(c_int32_t), value, intent(in) :: numZonesToRetain
    integer(c_int32_t), intent(in)        :: zonesToRetain(*)
end function tecFileWriterFlush

integer(c_int32_t) function tecFileWriterClose( &
    fileHandle) &
    bind(c, name="tecFileWriterClose")
    use iso_c_binding
    implicit none
    type(c_ptr), intent(inout) :: fileHandle
end function tecFileWriterClose

! Input routine interfaces
integer(c_int32_t) function tecCustomLabelsGetNumSets( &
  fileHandle, &
  numSets) &
  bind (c, name="tecCustomLabelsGetNumSets")
  use iso_c_binding
  implicit none
  type(c_ptr), value , intent(in) :: fileHandle
  integer(c_int32_t), intent(out) :: numSets
end function tecCustomLabelsGetNumSets

integer(c_int32_t) function tecCustomLabelsGetSet( &
  fileHandle, &
  whichSet, &
  labelSet) &
  bind (c, name="tecCustomLabelsGetSet")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)         :: fileHandle
  integer(c_int32_t), value , intent(in) :: whichSet
  type(c_ptr), intent(out)               :: labelSet
end function tecCustomLabelsGetSet

integer(c_int32_t) function tecDataSetAuxDataGetItem( &
  fileHandle, &
  whichItem, &
  name, &
  value) &
  bind (c, name="tecDataSetAuxDataGetItem")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: whichItem
  type(c_ptr), intent(out)              :: name
  type(c_ptr), intent(out)              :: value
end function tecDataSetAuxDataGetItem

integer(c_int32_t) function tecDataSetAuxDataGetNumItems( &
  fileHandle, &
  numItems) &
  bind (c, name="tecDataSetAuxDataGetNumItems")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)  :: fileHandle
  integer(c_int32_t), intent(out) :: numItems
end function tecDataSetAuxDataGetNumItems

integer(c_int32_t) function tecDataSetGetNumVars( &
  fileHandle, &
  numVars) &
  bind (c, name="tecDataSetGetNumVars")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)  :: fileHandle
  integer(c_int32_t), intent(out) :: numVars
end function tecDataSetGetNumVars

integer(c_int32_t) function tecDataSetGetNumZones( &
  fileHandle, &
  numZones) &
  bind (c, name="tecDataSetGetNumZones")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)  :: fileHandle
  integer(c_int32_t), intent(out) :: numZones
end function tecDataSetGetNumZones

integer(c_int32_t) function tecDataSetGetTitle( &
  fileHandle, &
  title) &
  bind (c, name="tecDataSetGetTitle")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in) :: fileHandle
  type(c_ptr), intent(out)       :: title
end function tecDataSetGetTitle

integer(c_int32_t) function tecFileGetType( &
  fileHandle, &
  fileType) &
  bind (c, name="tecFileGetType")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)  :: fileHandle
  integer(c_int32_t), intent(out) :: fileType
end function tecFileGetType

integer(c_int32_t) function tecFileReaderClose( &
  fileHandle) &
  bind (c, name="tecFileReaderClose")
  use iso_c_binding
  implicit none
  type(c_ptr), intent(inout) :: fileHandle
end function tecFileReaderClose

integer(c_int32_t) function tecFileReaderOpen( &
  fileName, &
  fileHandle) &
  bind (c, name="tecFileReaderOpen")
  use iso_c_binding
  implicit none
  character(c_char), intent(in) :: fileName(*)
  type(c_ptr), intent(out)      :: fileHandle
end function tecFileReaderOpen

integer(c_int32_t) function tecGeomArrowheadGetAngle( &
  fileHandle, &
  geom, &
  angle) &
  bind (c, name="tecGeomArrowheadGetAngle")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  real(c_double), intent(out)           :: angle
end function tecGeomArrowheadGetAngle

integer(c_int32_t) function tecGeomArrowheadGetAttach( &
  fileHandle, &
  geom, &
  attachment) &
  bind (c, name="tecGeomArrowheadGetAttach")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out)       :: attachment
end function tecGeomArrowheadGetAttach

integer(c_int32_t) function tecGeomArrowheadGetSize( &
  fileHandle, &
  geom, &
  arrowheadSize) &
  bind (c, name="tecGeomArrowheadGetSize")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  real(c_double), intent(out)           :: arrowheadSize
end function tecGeomArrowheadGetSize

integer(c_int32_t) function tecGeomArrowheadGetStyle( &
  fileHandle, &
  geom, &
  arrowheadStyle) &
  bind (c, name="tecGeomArrowheadGetStyle")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out)       :: arrowheadStyle
end function tecGeomArrowheadGetStyle

integer(c_int32_t) function tecGeomCircleGetRadius( &
  fileHandle, &
  geom, &
  radius) &
  bind (c, name="tecGeomCircleGetRadius")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  real(c_double), intent(out)           :: radius
end function tecGeomCircleGetRadius

integer(c_int32_t) function tecGeomEllipseGetNumPoints( &
  fileHandle, &
  geom, &
  numEllipsePoints) &
  bind (c, name="tecGeomEllipseGetNumPoints")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out)       :: numEllipsePoints
end function tecGeomEllipseGetNumPoints

integer(c_int32_t) function tecGeomEllipseGetSize( &
  fileHandle, &
  geom, &
  horizontalAxis, &
  verticalAxis) &
  bind (c, name="tecGeomEllipseGetSize")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  real(c_double), intent(out)           :: horizontalAxis
  real(c_double), intent(out)           :: verticalAxis
end function tecGeomEllipseGetSize

integer(c_int32_t) function tecGeomGetAnchorPos( &
  fileHandle, &
  geom, &
  x, &
  y, &
  z) &
  bind (c, name="tecGeomGetAnchorPos")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  real(c_double), intent(out)           :: x
  real(c_double), intent(out)           :: y
  real(c_double), intent(out)           :: z
end function tecGeomGetAnchorPos

integer(c_int32_t) function tecGeomGetClipping( &
  fileHandle, &
  geom, &
  clipping) &
  bind (c, name="tecGeomGetClipping")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out)       :: clipping
end function tecGeomGetClipping

integer(c_int32_t) function tecGeomGetColor( &
  fileHandle, &
  geom, &
  color) &
  bind (c, name="tecGeomGetColor")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in) :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out) :: color
end function tecGeomGetColor

integer(c_int32_t) function tecGeomGetCoordMode( &
  fileHandle, &
  geom, &
  coordMode) &
  bind (c, name="tecGeomGetCoordMode")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in) :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out) :: coordMode
end function tecGeomGetCoordMode

integer(c_int32_t) function tecGeomGetFillColor( &
  fileHandle, &
  geom, &
  fillColor) &
  bind (c, name="tecGeomGetFillColor")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out)       :: fillColor
end function tecGeomGetFillColor

integer(c_int32_t) function tecGeomGetLinePattern( &
  fileHandle, &
  geom, &
  linePattern) &
  bind (c, name="tecGeomGetLinePattern")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out)       :: linePattern
end function tecGeomGetLinePattern

integer(c_int32_t) function tecGeomGetLineThickness( &
  fileHandle, &
  geom, &
  lineThickness) &
  bind (c, name="tecGeomGetLineThickness")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  real(c_double), intent(out)           :: lineThickness
end function tecGeomGetLineThickness

integer(c_int32_t) function tecGeomGetMacroFunctionCmd( &
  fileHandle, &
  geom, &
  macroFunctionCmd) &
  bind (c, name="tecGeomGetMacroFunctionCmd")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  type(c_ptr), intent(out)              :: macroFunctionCmd
end function tecGeomGetMacroFunctionCmd

integer(c_int32_t) function tecGeomGetNumGeoms( &
  fileHandle, &
  numGeoms) &
  bind (c, name="tecGeomGetNumGeoms")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)  :: fileHandle
  integer(c_int32_t), intent(out) :: numGeoms
end function tecGeomGetNumGeoms

integer(c_int32_t) function tecGeomGetPatternLength( &
  fileHandle, &
  geom, &
  patternLength) &
  bind (c, name="tecGeomGetPatternLength")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  real(c_double), intent(out)           :: patternLength
end function tecGeomGetPatternLength

integer(c_int32_t) function tecGeomGetScope( &
  fileHandle, &
  geom, &
  scope) &
  bind (c, name="tecGeomGetScope")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out)       :: scope
end function tecGeomGetScope

integer(c_int32_t) function tecGeomGetType( &
  fileHandle, &
  geom, &
  type) &
  bind (c, name="tecGeomGetType")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out)       :: type
end function tecGeomGetType

integer(c_int32_t) function tecGeomGetZone( &
  fileHandle, &
  geom, &
  zone) &
  bind (c, name="tecGeomGetZone")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out)       :: zone
end function tecGeomGetZone

integer(c_int32_t) function tecGeomIsAttached( &
  fileHandle, &
  geom, &
  isAttached) &
  bind (c, name="tecGeomIsAttached")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in) :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out) :: isAttached
end function tecGeomIsAttached

integer(c_int32_t) function tecGeomIsFilled( &
  fileHandle, &
  geom, &
  isFilled) &
  bind (c, name="tecGeomIsFilled")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out)       :: isFilled
end function tecGeomIsFilled

integer(c_int32_t) function tecGeomLineGetPoint( &
  fileHandle, &
  geom, &
  segment, &
  index, &
  x, &
  y, &
  z) &
  bind (c, name="tecGeomLineGetPoint")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), value, intent(in) :: segment
  integer(c_int32_t), value, intent(in) :: index
  real(c_double), intent(out)           :: x
  real(c_double), intent(out)           :: y
  real(c_double), intent(out)           :: z
end function tecGeomLineGetPoint

integer(c_int32_t) function tecGeomLineGetSegmentCount( &
  fileHandle, &
  geom, &
  segmentCount) &
  bind (c, name="tecGeomLineGetSegmentCount")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), intent(out)       :: segmentCount
end function tecGeomLineGetSegmentCount

integer(c_int32_t) function tecGeomLineSegmentGetPointCount( &
  fileHandle, &
  geom, &
  segment, &
  pointCount) &
  bind (c, name="tecGeomLineSegmentGetPointCount")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  integer(c_int32_t), value, intent(in) :: segment
  integer(c_int32_t), intent(out)       :: pointCount
end function tecGeomLineSegmentGetPointCount

integer(c_int32_t) function tecGeomRectangleGetSize( &
  fileHandle, &
  geom, &
  width, &
  height) &
  bind (c, name="tecGeomRectangleGetSize")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  real(c_double), intent(out)           :: width
  real(c_double), intent(out)           :: height
end function tecGeomRectangleGetSize

integer(c_int32_t) function tecGeomSquareGetSize( &
  fileHandle, &
  geom, &
  size) &
  bind (c, name="tecGeomSquareGetSize")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: geom
  real(c_double), intent(out)           :: size
end function tecGeomSquareGetSize

subroutine tecStringFree( &
  string) &
  bind (c, name="tecStringFree")
  use iso_c_binding
  implicit none
  type(c_ptr), intent(inout) :: string
end subroutine tecStringFree

integer(c_int32_t) function tecStringLength( &
  string) &
  bind(c, name = "tecStringLength")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in) :: string
  end function tecStringLength

integer(c_int32_t) function tecTextBoxGetColor( &
  fileHandle, &
  text, &
  boxColor) &
  bind (c, name="tecTextBoxGetColor")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: boxColor
end function tecTextBoxGetColor

integer(c_int32_t) function tecTextBoxGetFillColor( &
  fileHandle, &
  text, &
  boxFillColor) &
  bind (c, name="tecTextBoxGetFillColor")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: boxFillColor
end function tecTextBoxGetFillColor

integer(c_int32_t) function tecTextBoxGetLineThickness( &
  fileHandle, &
  text, &
  boxLineThickness) &
  bind (c, name="tecTextBoxGetLineThickness")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  real(c_double), intent(out)           :: boxLineThickness
end function tecTextBoxGetLineThickness

integer(c_int32_t) function tecTextBoxGetMargin( &
  fileHandle, &
  text, &
  boxMargin) &
  bind (c, name="tecTextBoxGetMargin")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  real(c_double), intent(out)           :: boxMargin
end function tecTextBoxGetMargin

integer(c_int32_t) function tecTextBoxGetType( &
  fileHandle, &
  text, &
  boxType) &
  bind (c, name="tecTextBoxGetType")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: boxType
end function tecTextBoxGetType

integer(c_int32_t) function tecTextGetAnchor( &
  fileHandle, &
  text, &
  anchor) &
  bind (c, name="tecTextGetAnchor")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: anchor
end function tecTextGetAnchor

integer(c_int32_t) function tecTextGetAnchorPos( &
  fileHandle, &
  text, &
  x, &
  y, &
  z) &
  bind (c, name="tecTextGetAnchorPos")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  real(c_double), intent(out)           :: x
  real(c_double), intent(out)           :: y
  real(c_double), intent(out)           :: z
end function tecTextGetAnchorPos

integer(c_int32_t) function tecTextGetAngle( &
  fileHandle, &
  text, &
  angle) &
  bind (c, name="tecTextGetAngle")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  real(c_double), intent(out)           :: angle
end function tecTextGetAngle

integer(c_int32_t) function tecTextGetClipping( &
  fileHandle, &
  text, &
  clipping) &
  bind (c, name="tecTextGetClipping")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: clipping
end function tecTextGetClipping

integer(c_int32_t) function tecTextGetColor( &
  fileHandle, &
  text, &
  color) &
  bind (c, name="tecTextGetColor")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: color
end function tecTextGetColor

integer(c_int32_t) function tecTextGetCoordMode( &
  fileHandle, &
  text, &
  coordMode) &
  bind (c, name="tecTextGetCoordMode")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: coordMode
end function tecTextGetCoordMode

integer(c_int32_t) function tecTextGetHeight( &
  fileHandle, &
  text, &
  height) &
  bind (c, name="tecTextGetHeight")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  real(c_double), intent(out)           :: height
end function tecTextGetHeight

integer(c_int32_t) function tecTextGetLineSpacing( &
  fileHandle, &
  text, &
  lineSpacing) &
  bind (c, name="tecTextGetLineSpacing")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  real(c_double), intent(out)           :: lineSpacing
end function tecTextGetLineSpacing

integer(c_int32_t) function tecTextGetMacroFunctionCmd( &
  fileHandle, &
  text, &
  macroFunctionCmd) &
  bind (c, name="tecTextGetMacroFunctionCmd")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  type(c_ptr), intent(out)              :: macroFunctionCmd
end function tecTextGetMacroFunctionCmd

integer(c_int32_t) function tecTextGetNumTexts( &
  fileHandle, &
  numTexts) &
  bind (c, name="tecTextGetNumTexts")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)  :: fileHandle
  integer(c_int32_t), intent(out) :: numTexts
end function tecTextGetNumTexts

integer(c_int32_t) function tecTextGetScope( &
  fileHandle, &
  text, &
  scope) &
  bind (c, name="tecTextGetScope")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: scope
end function tecTextGetScope

integer(c_int32_t) function tecTextGetSizeUnits( &
  fileHandle, &
  text, &
  sizeUnits) &
  bind (c, name="tecTextGetSizeUnits")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: sizeUnits
end function tecTextGetSizeUnits

integer(c_int32_t) function tecTextGetString( &
  fileHandle, &
  text, &
  string) &
  bind (c, name="tecTextGetString")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  type(c_ptr), intent(out)              :: string
end function tecTextGetString

integer(c_int32_t) function tecTextGetTypeface( &
  fileHandle, &
  text, &
  typeface) &
  bind (c, name="tecTextGetTypeface")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  type(c_ptr), intent(out)              :: typeface
end function tecTextGetTypeface

integer(c_int32_t) function tecTextGetZone( &
  fileHandle, &
  text, &
  zone) &
  bind (c, name="tecTextGetZone")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: zone
end function tecTextGetZone

integer(c_int32_t) function tecTextIsAttached( &
  fileHandle, &
  text, &
  isAttached) &
  bind (c, name="tecTextIsAttached")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: isAttached
end function tecTextIsAttached

integer(c_int32_t) function tecTextIsBold( &
  fileHandle, &
  text, &
  isBold) &
  bind (c, name="tecTextIsBold")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: isBold
end function tecTextIsBold

integer(c_int32_t) function tecTextIsItalic( &
  fileHandle, &
  text, &
  isItalic) &
  bind (c, name="tecTextIsItalic")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: text
  integer(c_int32_t), intent(out)       :: isItalic
end function tecTextIsItalic

integer(c_int32_t) function tecVarAuxDataGetItem( &
  fileHandle, &
  var, &
  whichItem, &
  name, &
  value) &
  bind (c, name="tecVarAuxDataGetItem")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int32_t), value, intent(in) :: whichItem
  type(c_ptr), intent(out)              :: name
  type(c_ptr), intent(out)              :: value
end function tecVarAuxDataGetItem

integer(c_int32_t) function tecVarAuxDataGetNumItems( &
  fileHandle, &
  var, &
  numItems) &
  bind (c, name="tecVarAuxDataGetNumItems")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int32_t), intent(out)       :: numItems
end function tecVarAuxDataGetNumItems

integer(c_int32_t) function tecVarGetName( &
  fileHandle, &
  var, &
  name) &
  bind (c, name="tecVarGetName")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: var
  type(c_ptr), intent(out)              :: name
end function tecVarGetName

integer(c_int32_t) function tecVarIsEnabled( &
  fileHandle, &
  var, &
  isEnabled) &
  bind (c, name="tecVarIsEnabled")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int32_t), intent(out)       :: isEnabled
end function tecVarIsEnabled

integer(c_int32_t) function tecZoneAuxDataGetItem( &
  fileHandle, &
  zone, &
  whichItem, &
  name, &
  value) &
  bind (c, name="tecZoneAuxDataGetItem")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), value, intent(in) :: whichItem
  type(c_ptr), intent(out)              :: name
  type(c_ptr), intent(out)              :: value
end function tecZoneAuxDataGetItem

integer(c_int32_t) function tecZoneAuxDataGetNumItems( &
  fileHandle, &
  zone, &
  numItems) &
  bind (c, name="tecZoneAuxDataGetNumItems")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), intent(out)       :: numItems
end function tecZoneAuxDataGetNumItems

integer(c_int32_t) function tecZoneConnectivityGetSharedZone( &
  fileHandle, &
  zone, &
  sharedZone) &
  bind (c, name="tecZoneConnectivityGetSharedZone")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in) :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), intent(out) :: sharedZone
end function tecZoneConnectivityGetSharedZone

integer(c_int32_t) function tecZoneFaceNbrGetConnections( &
  fileHandle, &
  zone, &
  connections) &
  bind (c, name="tecZoneFaceNbrGetConnections")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), intent(out)       :: connections(*)
end function tecZoneFaceNbrGetConnections

integer(c_int32_t) function tecZoneFaceNbrGetConnections64( &
  fileHandle, &
  zone, &
  connections) &
  bind (c, name="tecZoneFaceNbrGetConnections64")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), intent(out)       :: connections(*)
end function tecZoneFaceNbrGetConnections64

integer(c_int32_t) function tecZoneFaceNbrGetNumConnections( &
  fileHandle, &
  zone, &
  numConnections) &
  bind (c, name="tecZoneFaceNbrGetNumConnections")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), intent(out)       :: numConnections
end function tecZoneFaceNbrGetNumConnections

integer(c_int32_t) function tecZoneFaceNbrGetMode( &
  fileHandle, &
  zone, &
  mode) &
  bind (c, name="tecZoneFaceNbrGetMode")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), intent(out)       :: mode
end function tecZoneFaceNbrGetMode

integer(c_int32_t) function tecZoneFaceNbrGetNumValues( &
  fileHandle, &
  zone, &
  numValues) &
  bind (c, name="tecZoneFaceNbrGetNumValues")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), intent(out)       :: numValues
end function tecZoneFaceNbrGetNumValues

integer(c_int32_t) function tecZoneFaceNbrsAre64Bit( &
  fileHandle, &
  zone, &
  are64Bit) &
  bind (c, name="tecZoneFaceNbrsAre64Bit")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), intent(out)       :: are64Bit
end function tecZoneFaceNbrsAre64Bit

integer(c_int32_t) function tecZoneGetIJK( &
  fileHandle, &
  zone, &
  iMax, &
  jMax, &
  kMax) &
  bind (c, name="tecZoneGetIJK")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), intent(out)       :: iMax
  integer(c_int64_t), intent(out)       :: jMax
  integer(c_int64_t), intent(out)       :: kMax
end function tecZoneGetIJK

integer(c_int32_t) function tecZoneGetParentZone( &
  fileHandle, &
  zone, &
  parentZone) &
  bind (c, name="tecZoneGetParentZone")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), intent(out)       :: parentZone
end function tecZoneGetParentZone

integer(c_int32_t) function tecZoneGetSolutionTime( &
  fileHandle, &
  zone, &
  solutionTime) &
  bind (c, name="tecZoneGetSolutionTime")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  real(c_double), intent(out)           :: solutionTime
end function tecZoneGetSolutionTime

integer(c_int32_t) function tecZoneGetStrandID( &
  fileHandle, &
  zone, &
  strandID) &
  bind (c, name="tecZoneGetStrandID")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), intent(out)       :: strandID
end function tecZoneGetStrandID

integer(c_int32_t) function tecZoneGetTitle( &
  fileHandle, &
  zone, &
  title) &
  bind (c, name="tecZoneGetTitle")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  type(c_ptr), intent(out)              :: title
end function tecZoneGetTitle

integer(c_int32_t) function tecZoneGetType( &
  fileHandle, &
  zone, &
  type) &
  bind (c, name="tecZoneGetType")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), intent(out)       :: type
end function tecZoneGetType

integer(c_int32_t) function tecZoneIsEnabled( &
  fileHandle, &
  zone, &
  isEnabled) &
  bind (c, name="tecZoneIsEnabled")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), intent(out)       :: isEnabled
end function tecZoneIsEnabled

integer(c_int32_t) function tecZoneNodeMapGet( &
  fileHandle, &
  zone, &
  startCell, &
  numCells, &
  nodeMap) &
  bind (c, name="tecZoneNodeMapGet")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), value, intent(in) :: startCell
  integer(c_int64_t), value, intent(in) :: numCells
  integer(c_int32_t), intent(out)       :: nodeMap(*)
end function tecZoneNodeMapGet

integer(c_int32_t) function tecZoneNodeMapGet64( &
  fileHandle, &
  zone, &
  startCell, &
  numCells, &
  nodeMap) &
  bind (c, name="tecZoneNodeMapGet64")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), value, intent(in) :: startCell
  integer(c_int64_t), value, intent(in) :: numCells
  integer(c_int64_t), intent(out)       :: nodeMap(*)
end function tecZoneNodeMapGet64

integer(c_int32_t) function tecZoneNodeMapGetNumValues( &
  fileHandle, &
  zone, &
  numCells, &
  numValues) &
  bind (c, name="tecZoneNodeMapGetNumValues")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), value, intent(in) :: numCells
  integer(c_int64_t), intent(out)       :: numValues
end function tecZoneNodeMapGetNumValues

integer(c_int32_t) function tecZoneNodeMapIs64Bit( &
  fileHandle, &
  zone, &
  is64Bit) &
  bind (c, name="tecZoneNodeMapIs64Bit")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), intent(out)       :: is64Bit
end function tecZoneNodeMapIs64Bit

integer(c_int32_t) function tecZonePolyGetBoundaryConnectionCounts( &
  fileHandle, &
  zone, &
  startConnection, &
  numConnections, &
  connectionCounts) &
  bind (c, name="tecZonePolyGetBoundaryConnectionCounts")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), value, intent(in) :: startConnection
  integer(c_int64_t), value, intent(in) :: numConnections
  integer(c_int32_t), intent(out)       :: connectionCounts(*)
end function tecZonePolyGetBoundaryConnectionCounts

integer(c_int32_t) function tecZonePolyGetBoundaryConnections( &
  fileHandle, &
  zone, &
  startConnection, &
  numConnections, &
  connectedElements, &
  connectedZones) &
  bind (c, name="tecZonePolyGetBoundaryConnections")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), value, intent(in) :: startConnection
  integer(c_int64_t), value, intent(in) :: numConnections
  integer(c_int32_t), intent(out)       :: connectedElements(*)
  integer(c_int32_t), intent(out)       :: connectedZones(*)
end function tecZonePolyGetBoundaryConnections

integer(c_int32_t) function tecZonePolyGetFaceElems( &
  fileHandle, &
  zone, &
  startFace, &
  numFaces, &
  leftElems, &
  rightElems) &
  bind (c, name="tecZonePolyGetFaceElems")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), value, intent(in) :: startFace
  integer(c_int64_t), value, intent(in) :: numFaces
  integer(c_int32_t), intent(out)       :: leftElems(*)
  integer(c_int32_t), intent(out)       :: rightElems(*)
end function tecZonePolyGetFaceElems

integer(c_int32_t) function tecZonePolyGetFaceNodeCounts( &
  fileHandle, &
  zone, &
  startFace, &
  numFaces, &
  nodeCounts) &
  bind (c, name="tecZonePolyGetFaceNodeCounts")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), value, intent(in) :: startFace
  integer(c_int64_t), value, intent(in) :: numFaces
  integer(c_int32_t), intent(out)       :: nodeCounts(*)
end function tecZonePolyGetFaceNodeCounts

integer(c_int32_t) function tecZonePolyGetFaceNodes( &
  fileHandle, &
  zone, &
  startFace, &
  numFaces, &
  faceNodes) &
  bind (c, name="tecZonePolyGetFaceNodes")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), value, intent(in) :: startFace
  integer(c_int64_t), value, intent(in) :: numFaces
  integer(c_int32_t), intent(out)       :: faceNodes(*)
end function tecZonePolyGetFaceNodes

integer(c_int32_t) function tecZonePolyGetNumConnectedBoundaryFaces( &
  fileHandle, &
  zone, &
  numFaces) &
  bind (c, name="tecZonePolyGetNumConnectedBoundaryFaces")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), intent(out)       :: numFaces(*)
end function tecZonePolyGetNumConnectedBoundaryFaces

integer(c_int32_t) function tecZonePolyGetTotalNumFaceNodes( &
  fileHandle, &
  zone, &
  numNodes) &
  bind (c, name="tecZonePolyGetTotalNumFaceNodes")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), intent(out)       :: numNodes
end function tecZonePolyGetTotalNumFaceNodes

integer(c_int32_t) function tecZonePolyGetTotalNumBoundaryConnections( &
  fileHandle, &
  zone, &
  numConnections) &
  bind (c, name="tecZonePolyGetTotalNumBoundaryConnections")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int64_t), intent(out)       :: numConnections
end function tecZonePolyGetTotalNumBoundaryConnections

integer(c_int32_t) function tecZoneVarGetDoubleValues( &
  fileHandle, &
  zone, &
  var, &
  startIndex, &
  numValues, &
  values) &
  bind (c, name="tecZoneVarGetDoubleValues")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int64_t), value, intent(in) :: startIndex
  integer(c_int64_t), value, intent(in) :: numValues
  real(c_double), intent(out)           :: values(*)
end function tecZoneVarGetDoubleValues

integer(c_int32_t) function tecZoneVarGetFloatValues( &
  fileHandle, &
  zone, &
  var, &
  startIndex, &
  numValues, &
  values) &
  bind (c, name="tecZoneVarGetFloatValues")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int64_t), value, intent(in) :: startIndex
  integer(c_int64_t), value, intent(in) :: numValues
  real(c_float), intent(out)            :: values(*)
end function tecZoneVarGetFloatValues

integer(c_int32_t) function tecZoneVarGetInt16Values( &
  fileHandle, &
  zone, &
  var, &
  startIndex, &
  numValues, &
  values) &
  bind (c, name="tecZoneVarGetInt16Values")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int64_t), value, intent(in) :: startIndex
  integer(c_int64_t), value, intent(in) :: numValues
  integer(c_int16_t), intent(out)       :: values(*)
end function tecZoneVarGetInt16Values

integer(c_int32_t) function tecZoneVarGetInt32Values( &
  fileHandle, &
  zone, &
  var, &
  startIndex, &
  numValues, &
  values) &
  bind (c, name="tecZoneVarGetInt32Values")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int64_t), value, intent(in) :: startIndex
  integer(c_int64_t), value, intent(in) :: numValues
  integer(c_int32_t), intent(out)       :: values(*)
end function tecZoneVarGetInt32Values

integer(c_int32_t) function tecZoneVarGetNumValues( &
  fileHandle, &
  zone, &
  var, &
  numValues) &
  bind (c, name="tecZoneVarGetNumValues")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int64_t), intent(out)       :: numValues
end function tecZoneVarGetNumValues

integer(c_int32_t) function tecZoneVarGetSharedZone( &
  fileHandle, &
  zone, &
  var, &
  sharedZone) &
  bind (c, name="tecZoneVarGetSharedZone")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int32_t), intent(out)       :: sharedZone
end function tecZoneVarGetSharedZone

integer(c_int32_t) function tecZoneVarGetType( &
  fileHandle, &
  zone, &
  var, &
  type) &
  bind (c, name="tecZoneVarGetType")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in) :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int32_t), intent(out) :: type
end function tecZoneVarGetType

integer(c_int32_t) function tecZoneVarGetUInt8Values( &
  fileHandle, &
  zone, &
  var, &
  startIndex, &
  numValues, &
  values) &
  bind (c, name="tecZoneVarGetUInt8Values")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int64_t), value, intent(in) :: startIndex
  integer(c_int64_t), value, intent(in) :: numValues
  integer(c_int8_t), intent(out)        :: values(*)
end function tecZoneVarGetUInt8Values

integer(c_int32_t) function tecZoneVarGetValueLocation( &
  fileHandle, &
  zone, &
  var, &
  location) &
  bind (c, name="tecZoneVarGetValueLocation")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int32_t), intent(out)       :: location
end function tecZoneVarGetValueLocation

integer(c_int32_t) function tecZoneVarIsPassive( &
  fileHandle, &
  zone, &
  var, &
  isPassive) &
  bind (c, name="tecZoneVarIsPassive")
  use iso_c_binding
  implicit none
  type(c_ptr), value, intent(in)        :: fileHandle
  integer(c_int32_t), value, intent(in) :: zone
  integer(c_int32_t), value, intent(in) :: var
  integer(c_int32_t), intent(out)       :: isPassive
end function tecZoneVarIsPassive

! Older output routine interfaces
integer(c_int32_t) function tecini142( &
  Title, &
  Variables, &
  FName, &
  ScratchDir, &
  FileFormat, &
  FileType, &
  Debug, &
  VIsDouble) &
  bind (c, name="tecini142")
  use iso_c_binding
  implicit none
  character(c_char), intent(in)  :: Title(*)
  character(c_char), intent(in)  :: Variables(*)
  character(c_char), intent(in)  :: FName(*)
  character(c_char), intent(in)  :: ScratchDir(*)
  integer(c_int32_t), intent(in) :: FileFormat
  integer(c_int32_t), intent(in) :: FileType
  integer(c_int32_t), intent(in) :: Debug
  integer(c_int32_t), intent(in) :: VIsDouble
end function tecini142

integer(c_int32_t) function teczne142( &
  ZoneTitle, &
  ZoneType, &
  IMxOrNumPts, &
  JMxOrNumElements, &
  KMxOrNumFaces, &
  ICellMax, &
  JCellMax, &
  KCellMax, &
  SolutionTime, &
  StrandID, &
  ParentZone, &
  IsBlock, &
  NumFaceConnections, &
  FaceNeighborMode, &
  TotalNumFaceNodes, &
  NumConnectedBoundaryFaces, &
  TotalNumBoundaryConnections, &
  PassiveVarList, &
  ValueLocation, &
  ShareVarFromZone, &
  ShareConnectivityFromZone) &
  bind (c, name="teczne142")
  use iso_c_binding
  implicit none
  character(c_char), intent(in)  :: ZoneTitle(*)
  integer(c_int32_t), intent(in) :: ZoneType
  integer(c_int32_t), intent(in) :: IMxOrNumPts
  integer(c_int32_t), intent(in) :: JMxOrNumElements
  integer(c_int32_t), intent(in) :: KMxOrNumFaces
  integer(c_int32_t), intent(in) :: ICellMax
  integer(c_int32_t), intent(in) :: JCellMax
  integer(c_int32_t), intent(in) :: KCellMax
  real(c_double), intent(in)     :: SolutionTime
  integer(c_int32_t), intent(in) :: StrandID
  integer(c_int32_t), intent(in) :: ParentZone
  integer(c_int32_t), intent(in) :: IsBlock
  integer(c_int32_t), intent(in) :: NumFaceConnections
  integer(c_int32_t), intent(in) :: FaceNeighborMode
  integer(c_int32_t), intent(in) :: TotalNumFaceNodes
  integer(c_int32_t), intent(in) :: NumConnectedBoundaryFaces
  integer(c_int32_t), intent(in) :: TotalNumBoundaryConnections
  integer(c_int32_t), intent(in) :: PassiveVarList(*)
  integer(c_int32_t), intent(in) :: ValueLocation(*)
  integer(c_int32_t), intent(in) :: ShareVarFromZone(*)
  integer(c_int32_t), intent(in) :: ShareConnectivityFromZone
end function teczne142

integer(c_int32_t) function tecdat142( &
  N, &
  FieldData, &
  IsDouble) &
  bind (c, name="tecdat142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: N, IsDouble
  real(c_float), intent(in)   :: FieldData(*)
end function tecdat142

integer(c_int32_t) function tecdatf142( &
  N, &
  FieldData) &
  bind (c, name="tecdatf142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: N
  real(c_float), intent(in)   :: FieldData(*)
end function tecdatf142

integer(c_int32_t) function tecdatd142( &
  N, &
  FieldData) &
  bind (c, name="tecdatd142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: N
  real(c_double), intent(in)  :: FieldData(*)
end function tecdatd142

integer(c_int32_t) function tecnod142( &
  NData) &
  bind (c, name="tecnod142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: NData(*)
end function tecnod142

integer(c_int32_t) function tecnode142( &
  N, &
  NData) &
  bind (c, name="tecnode142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: N
  integer(c_int32_t), intent(in) :: NData(*)
end function tecnode142

integer(c_int32_t) function tecgeo142( &
  XPos, &
  YPos, &
  ZPos, &
  PosCoordMode, &
  AttachToZone, &
  Zone, &
  Color, &
  FillColor, &
  IsFilled, &
  GeomType, &
  LinePattern, &
  PatternLength, &
  LineThickness, &
  NumEllipsePts, &
  ArrowheadStyle, &
  ArrowheadAttachment, &
  ArrowheadSize, &
  ArrowheadAngle, &
  Scope, &
  Clipping, &
  NumSegments, &
  NumSegPts, &
  XGeomData, &
  YGeomData, &
  ZGeomData, &
  mfc) &
  bind (c, name="tecgeo142")
  use iso_c_binding
  implicit none
  real(c_double), intent(in)     :: XPos
  real(c_double), intent(in)     :: YPos
  real(c_double), intent(in)     :: ZPos
  integer(c_int32_t), intent(in) :: PosCoordMode
  integer(c_int32_t), intent(in) :: AttachToZone
  integer(c_int32_t), intent(in) :: Zone
  integer(c_int32_t), intent(in) :: Color
  integer(c_int32_t), intent(in) :: FillColor
  integer(c_int32_t), intent(in) :: IsFilled
  integer(c_int32_t), intent(in) :: GeomType
  integer(c_int32_t), intent(in) :: LinePattern
  real(c_double), intent(in)     :: PatternLength
  real(c_double), intent(in)     :: LineThickness
  integer(c_int32_t), intent(in) :: NumEllipsePts
  integer(c_int32_t), intent(in) :: ArrowheadStyle
  integer(c_int32_t), intent(in) :: ArrowheadAttachment
  real(c_double), intent(in)     :: ArrowheadSize
  real(c_double), intent(in)     :: ArrowheadAngle
  integer(c_int32_t), intent(in) :: Scope
  integer(c_int32_t), intent(in) :: Clipping
  integer(c_int32_t), intent(in) :: NumSegments
  integer(c_int32_t), intent(in) :: NumSegPts(*)
  real(c_float), intent(in)      :: XGeomData(*)
  real(c_float), intent(in)      :: YGeomData(*)
  real(c_float), intent(in)      :: ZGeomData(*)
  character(c_char), intent(in)  :: mfc(*)
end function tecgeo142

integer(c_int32_t) function tectxt142( &
  XOrThetaPos, &
  YOrRPos, &
  ZOrUnusedPos, &
  PosCoordMode, &
  AttachToZone, &
  Zone, &
  Font, &
  FontHeightUnits, &
  FontHeight, &
  BoxType, &
  BoxMargin, &
  BoxLineThickness, &
  BoxColor, &
  BoxFillColor, &
  Angle, &
  Anchor, &
  LineSpacing, &
  TextColor, &
  Scope, &
  Clipping, &
  Text, &
  mfc) &
  bind (c, name="tectxt142")
  use iso_c_binding
  implicit none
  real(c_double), intent(in)     :: XOrThetaPos
  real(c_double), intent(in)     :: YOrRPos
  real(c_double), intent(in)     :: ZOrUnusedPos
  integer(c_int32_t), intent(in) :: PosCoordMode
  integer(c_int32_t), intent(in) :: AttachToZone
  integer(c_int32_t), intent(in) :: Zone
  integer(c_int32_t), intent(in) :: Font
  integer(c_int32_t), intent(in) :: FontHeightUnits
  real(c_double), intent(in)     :: FontHeight
  integer(c_int32_t), intent(in) :: BoxType
  real(c_double), intent(in)     :: BoxMargin
  real(c_double), intent(in)     :: BoxLineThickness
  integer(c_int32_t), intent(in) :: BoxColor
  integer(c_int32_t), intent(in) :: BoxFillColor
  real(c_double), intent(in)     :: Angle
  integer(c_int32_t), intent(in) :: Anchor
  real(c_double), intent(in)     :: LineSpacing
  integer(c_int32_t), intent(in) :: TextColor
  integer(c_int32_t), intent(in) :: Scope
  integer(c_int32_t), intent(in) :: Clipping
  character(c_char), intent(in)  ::  Text(*)
  character(c_char), intent(in)  ::  mfc(*)
end function tectxt142

integer(c_int32_t) function teclab142( &
  S) &
  bind (c, name="teclab142")
  use iso_c_binding
  implicit none
  character(c_char), intent(in) :: S
end function teclab142

integer(c_int32_t) function tecfil142( &
  F) &
  bind (c, name="tecfil142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: F
end function tecfil142

subroutine tecforeign142( &
  OutputForeignByteOrder) &
  bind (c, name="tecforeign142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: OutputForeignByteOrder
end subroutine tecforeign142
    
integer(c_int32_t) function tecflush142( &
  numZonesToRetain, &
  zonesToRetain) &
  bind(c, name="tecflush142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), value, intent(in) :: numZonesToRetain
  integer(c_int32_t), intent(in)        :: zonesToRetain(*)
end function tecflush142

integer(c_int32_t) function tecend142() &
  bind (c, name="tecend142")
  use iso_c_binding
  implicit none
end function tecend142

integer(c_int32_t) function tecusr142( &
  S) &
  bind (c, name="tecusr142")
  use iso_c_binding
  implicit none
  character(c_char), intent(in) :: S(*)
end function tecusr142

integer(c_int32_t) function tecauxstr142( &
  Name, &
  Value) &
  bind (c, name="tecauxstr142")
  use iso_c_binding
  implicit none
  character(c_char), intent(in) :: Name(*)
  character(c_char), intent(in) :: Value(*)
end function tecauxstr142

integer(c_int32_t) function teczauxstr142( &
  Name, &
  Value) &
  bind (c, name="teczauxstr142")
  use iso_c_binding
  implicit none
  character(c_char), intent(in) :: Name(*)
  character(c_char), intent(in) :: Value(*)
end function teczauxstr142

integer(c_int32_t) function tecvauxstr142( &
  Var, &
  Name, &
  Value) &
  bind (c, name="tecvauxstr142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: Var
  character(c_char), intent(in)  :: Name(*)
  character(c_char), intent(in)  :: Value(*)
end function tecvauxstr142

integer(c_int32_t) function tecface142( &
  FaceConnections) &
  bind (c, name="tecface142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: FaceConnections(*)
end function tecface142

integer(c_int32_t) function tecpoly142( &
  FaceNodeCounts, &
  FaceNodes, &
  FaceLeftElems, &
  FaceRightElems, &
  FaceBndryConnectionCounts, &
  FaceBndryConnectionElems, &
  FaceBndryConnectionZones) &
  bind (c, name="tecpoly142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: FaceNodeCounts(*)
  integer(c_int32_t), intent(in) :: FaceNodes(*)
  integer(c_int32_t), intent(in) :: FaceLeftElems(*)
  integer(c_int32_t), intent(in) :: FaceRightElems(*)
  integer(c_int32_t), intent(in) :: FaceBndryConnectionCounts(*)
  integer(c_int32_t), intent(in) :: FaceBndryConnectionElems(*)
  integer(c_int32_t), intent(in) :: FaceBndryConnectionZones(*)
end function tecpoly142

integer(c_int32_t) function tecpolyface142( &
  NumFaces, &
  FaceNodeCounts, &
  FaceNodes, &
  FaceLeftElems, &
  FaceRightElems) &
  bind (c, name="tecpolyface142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: NumFaces
  integer(c_int32_t), intent(in) :: FaceNodeCounts(*)
  integer(c_int32_t), intent(in) :: FaceNodes(*)
  integer(c_int32_t), intent(in) :: FaceLeftElems(*)
  integer(c_int32_t), intent(in) :: FaceRightElems(*)
end function tecpolyface142

integer(c_int32_t) function tecpolybconn142( &
  NumBndryFaces, &
  FaceBndryConnectionCounts, &
  FaceBndryConnectionElems, &
  FaceBndryConnectionZones) &
  bind (c, name="tecpolybconn142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: NumBndryFaces
  integer(c_int32_t), intent(in) :: FaceBndryConnectionCounts(*)
  integer(c_int32_t), intent(in) :: FaceBndryConnectionElems(*) 
  integer(c_int32_t), intent(in) :: FaceBndryConnectionZones(*)
end function tecpolybconn142

integer(c_int32_t) function tecfeptn142( &
  partition, &
  numnodes, &
  numcells, &
  ngnodes, &
  gnodes, &
  gnpartitions, &
  gnpnodes, &
  ngcells, &
  gcells) &
  bind (c, name="tecfeptn142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: partition
  integer(c_int32_t), intent(in) :: numnodes
  integer(c_int32_t), intent(in) :: numcells
  integer(c_int32_t), intent(in) :: ngnodes
  integer(c_int32_t), intent(in) :: gnodes(*)
  integer(c_int32_t), intent(in) :: gnpartitions(*)
  integer(c_int32_t), intent(in) :: gnpnodes(*)
  integer(c_int32_t), intent(in) :: ngcells
  integer(c_int32_t), intent(in) :: gcells(*)
end function tecfeptn142

integer(c_int32_t) function tecijkptn142( &
  partition, &
  imin, &
  jmin, &
  kmin, &
  imax, &
  jmax, &
  kmax) &
  bind (c, name="tecijkptn142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: partition
  integer(c_int32_t), intent(in) :: imin
  integer(c_int32_t), intent(in) :: jmin
  integer(c_int32_t), intent(in) :: kmin
  integer(c_int32_t), intent(in) :: imax
  integer(c_int32_t), intent(in) :: jmax
  integer(c_int32_t), intent(in) :: kmax
end function tecijkptn142

integer(c_int32_t) function tecmpiinit142( &
  communicator, &
  mainrank) &
  bind (c, name="tecmpiinit142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: communicator
  integer(c_int32_t), intent(in) :: mainrank
end function tecmpiinit142

integer(c_int32_t) function tecznemap142( &
  npartitions, &
  ptnranks) &
  bind (c, name="tecznemap142")
  use iso_c_binding
  implicit none
  integer(c_int32_t), intent(in) :: npartitions
  integer(c_int32_t), intent(in) :: ptnranks(*)
end function tecznemap142

end interface
