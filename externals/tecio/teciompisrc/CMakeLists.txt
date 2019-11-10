project(teciompi CXX)

set(SU2_TECIOMPI_CPP_FLAGS
    "TP_PROJECT_USES_BOOST;BOOST_ALL_NO_LIB"
    CACHE STRING "Specific TecioMPI compile definitions to use")

find_package(MPI REQUIRED)

set(TECIOMPI_SOURCES
    ClassicZoneWriterAbstract.cpp
    SZLOrderedPartitionedZoneWriterMPI.cpp
    szcombine.cpp
    FileSystem.cpp
    MPIUtil.cpp
    SZLOrderedZoneHeaderWriter.cpp
    ZoneWriterFactory.cpp
    importSzPltFile.cpp
    FaceNeighborGeneratorAbstract.cpp
    PartitionTecUtilDecorator.cpp
    ZoneWriterFactoryMPI.cpp
    tecio.cpp
    TecioData.cpp
    MPIFileIOStream.cpp
    checkPercentDone.cpp
    Zone_s.cpp
    DataSetWriter.cpp
    ClassicZoneVariableWriter.cpp
    ClassicZoneFaceNeighborWriter.cpp
    UnicodeStringUtils.cpp
    ClassicZoneHeaderWriter.cpp
    SZLOrderedPartitionedZoneWriter.cpp
    SZLFEPartitionedZoneWriter.cpp
    MPICommunicator.cpp
    writeValueArray.cpp
    FECellSubzoneCompressor.cpp
    fileStuff.cpp
    MPICommunicationCache.cpp
    FileStreamReader.cpp
    SZLFEPartitionedZoneWriterMPI.cpp
    SZLOrderedPartitionWriter.cpp
    MPIFileWriter.cpp
    ZoneHeaderWriterAbstract.cpp
    TecioTecUtil.cpp
    ZoneWriterAbstract.cpp
    NodeToElemMap_s.cpp
    DataSetWriterMPI.cpp
    FileIOStream.cpp
    SZLFEPartitionedZoneHeaderWriter.cpp
    MPIFileReader.cpp
    ClassicFEZoneConnectivityWriter.cpp
    MinMaxTree.cpp
    FieldData.cpp
    SZLOrderedPartitionedZoneHeaderWriter.cpp
    AsciiOutputInfo.cpp
    FieldData_s.cpp
    IntervalTree.cpp
    ClassicFEZoneFaceNeighborGenerator.cpp
    mpiDatatype.cpp
    ClassicOrderedZoneFaceNeighborGenerator.cpp
    OrthogonalBisection.cpp
    ORBFESubzonePartitioner.cpp
    IJKZoneInfo.cpp
    IJKSubzoneInfo.cpp
    Scanner.cpp
    ClassicFEZoneWriter.cpp
    NodeMap.cpp
    SZLOrderedZoneWriter.cpp
    SZLFEPartitionWriter.cpp
    SZLFEZoneWriter.cpp
    ZoneVarMetadata.cpp
    MPINonBlockingCommunicationCollection.cpp
    SZLFEZoneHeaderWriter.cpp
    FileStreamWriter.cpp
    ClassicOrderedZoneWriter.cpp
    ZoneInfoCache.cpp
    readValueArray.cpp
    TecioSZL.cpp
    exportSubzonePlt.cpp)

add_library(teciompi STATIC ${TECIOMPI_SOURCES})
set_target_properties(teciompi PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(teciompi PRIVATE ${SU2_TECIO_CPPFLAGS} ${SU2_TECIOMPI_CPPFLAGS} ${TECIO_CPPFLAGS})
target_compile_definitions(teciompi PUBLIC TECIOMPI)

target_link_libraries(teciompi PUBLIC MPI::MPI_CXX)
target_include_directories(teciompi PUBLIC ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/../teciosrc)
add_library(Tecio::Tecio ALIAS teciompi)
