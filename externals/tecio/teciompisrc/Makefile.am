pkg_cppflags = @TECIO_CPPFLAGS@ -DUSEENUM -DTHREED -DMAKEARCHIVE -DNO_ASSERTS -DNO_THIRD_PARTY_LIBS -DTECIOMPI -DLINUX -DLINUX64 -DTP_PROJECT_USES_BOOST -DBOOST_ALL_NO_LIB
pkg_cppflags += -I$(top_srcdir)/externals/tecio
pkg_sources  = \
        AltTecUtil.h \
        AnyTypeLightweightVector.h \
        AsciiOutputInfo.h \
        AtomicMinMax.h \
        AuxData_s.h \
        basicTypes.h \
        CHARTYPE.h \
        checkPercentDone.h \
        clampToDataTypeRange.h \
        ClassicFEZoneConnectivityWriter.h \
        ClassicFEZoneFaceNeighborGenerator.h \
        ClassicFEZoneWriter.h \
        ClassicOrderedZoneFaceNeighborGenerator.h \
        ClassicOrderedZoneWriter.h \
        ClassicZoneFaceNeighborWriter.h \
        ClassicZoneFileLocations.h \
        ClassicZoneHeaderWriter.h \
        ClassicZoneVariableWriter.h \
        ClassicZoneWriterAbstract.h \
        ClassMacros.h \
        CodeContract.h \
        CszConnectivity.h \
        DataSetWriter.h \
        DataSetWriterMPI.h \
        DataWriteStatisticsInterface.h \
        exportSubzonePlt.h \
        FaceNeighborGeneratorAbstract.h \
        FECellSubzoneCompressor.h \
        FESubzonePartitionerInterface.h \
        FEZoneInfo.h \
        FieldData.h \
        FieldData_s.h \
        FileDescription.h \
        fileio.h \
        FileIOStatistics.h \
        FileIOStream.h \
        FileIOStreamInterface.h \
        FileReaderInterface.h \
        FileStreamReader.h \
        FileStreamWriter.h \
        fileStuff.h \
        FileSystem.h \
        FileWriterInterface.h \
        gatherOffsets.h \
        Geom_s.h \
        GhostInfo_s.h \
        GLOBAL.h \
        IJK.h \
        IJKPartitionTree.h \
        IJKSubzoneInfo.h \
        IJKZoneInfo.h \
        importSzPltFile.h \
        IntervalTree.h \
        ioDescription.h \
        ItemAddress.h \
        ItemSetIterator.h \
        JobControl_s.h \
        LightweightVector.h \
        MASTER.h \
        MinMax.h \
        MinMaxTree.h \
        MiscMacros.h \
        MPICommunicator.h \
        MPICommunicationCache.h \
        mpiDatatype.h \
        MPIError.h \
        MPIFileIOStream.h \
        MPIFileWriter.h \
        MPIFileReader.h \
        MPINonBlockingCommunicationCollection.h \
        MPIUtil.h \
        Mutex_s.h \
        NodeMap.h \
        NodeMap_s.h \
        NodeToElemMap_s.h \
        NoOpFESubzonePartitioner.h \
        ORBFESubzonePartitioner.h \
        OrthogonalBisection.h \
        PartitionMetadata.h \
        PartitionTecUtilDecorator.h \
        RawArray.h \
        readValueArray.h \
        Scanner.h \
        showMessage.h \
        SimpleVector.h \
        StandardIntegralTypes.h \
        stdafx.h \
        stringformat.h \
        SZLFEPartitionedZoneHeaderWriter.h \
        SZLFEPartitionedZoneWriter.h \
        SZLFEPartitionedZoneWriterMPI.h \
        SZLFEPartitionWriter.h \
        SZLFEZoneHeaderWriter.h \
        SZLFEZoneWriter.h \
        SzlFileLoader.h \
        SZLOrderedPartitionedZoneHeaderWriter.h \
        SZLOrderedPartitionedZoneWriter.h \
        SZLOrderedPartitionedZoneWriterMPI.h \
        SZLOrderedPartitionWriter.h \
        SZLOrderedZoneHeaderWriter.h \
        SZLOrderedZoneWriter.h \
        TASSERT.h \
        TECGLBL.h \
        TecioData.h \
        tecio_Exports.h \
        TECIO.h \
        TecioMPI.h \
        TecioSZL.h \
        TecioTecUtil.h \
        TecplotMinorRev.h \
        TecplotVersion.h \
        Text_s.h \
        ThirdPartyHeadersBegin.h \
        ThirdPartyHeadersEnd.h \
        TranslatedString.h \
        UnicodeStringUtils.h \
        writeValueArray.h \
        xyz.h \
        ZoneHeaderWriterAbstract.h \
        ZoneInfoCache.h \
        ZoneMetadata.h \
        Zone_s.h \
        zoneUtil.h \
        ZoneVarMetadata.h \
        ZoneWriterAbstract.h \
        ZoneWriterFactory.h \
        ZoneWriterFactoryMPI.h \
        AsciiOutputInfo.cpp \
        checkPercentDone.cpp \
        ClassicFEZoneConnectivityWriter.cpp \
        ClassicFEZoneFaceNeighborGenerator.cpp \
        ClassicFEZoneWriter.cpp \
        ClassicOrderedZoneFaceNeighborGenerator.cpp \
        ClassicOrderedZoneWriter.cpp \
        ClassicZoneFaceNeighborWriter.cpp \
        ClassicZoneHeaderWriter.cpp \
        ClassicZoneVariableWriter.cpp \
        ClassicZoneWriterAbstract.cpp \
        DataSetWriter.cpp \
        DataSetWriterMPI.cpp \
        exportSubzonePlt.cpp \
        FaceNeighborGeneratorAbstract.cpp \
        FECellSubzoneCompressor.cpp \
        FieldData.cpp \
        FieldData_s.cpp \
        FileIOStream.cpp \
        FileStreamReader.cpp \
        FileStreamWriter.cpp \
        fileStuff.cpp \
        FileSystem.cpp \
        IJKSubzoneInfo.cpp \
        IJKZoneInfo.cpp \
        importSzPltFile.cpp \
        IntervalTree.cpp \
        MinMaxTree.cpp \
        MPICommunicationCache.cpp \
        MPICommunicator.cpp \
        mpiDatatype.cpp \
        MPIFileIOStream.cpp \
        MPIFileReader.cpp \
        MPIFileWriter.cpp \
        MPIUtil.cpp \
        MPINonBlockingCommunicationCollection.cpp \
        NodeMap.cpp \
        NodeToElemMap_s.cpp \
        ORBFESubzonePartitioner.cpp \
        OrthogonalBisection.cpp \
        PartitionTecUtilDecorator.cpp \
        readValueArray.cpp \
        Scanner.cpp \
        szcombine.cpp \
        SZLFEPartitionedZoneHeaderWriter.cpp \
        SZLFEPartitionedZoneWriter.cpp \
        SZLFEPartitionedZoneWriterMPI.cpp \
        SZLFEPartitionWriter.cpp \
        SZLFEZoneHeaderWriter.cpp \
        SZLFEZoneWriter.cpp \
        SZLOrderedPartitionedZoneHeaderWriter.cpp \
        SZLOrderedPartitionedZoneWriter.cpp \
        SZLOrderedPartitionedZoneWriterMPI.cpp \
        SZLOrderedPartitionWriter.cpp \
        SZLOrderedZoneHeaderWriter.cpp \
        SZLOrderedZoneWriter.cpp \
        tecio.cpp \
        TecioData.cpp \
        TecioSZL.cpp \
        TecioTecUtil.cpp \
        UnicodeStringUtils.cpp \
        writeValueArray.cpp \
        ZoneHeaderWriterAbstract.cpp \
        ZoneInfoCache.cpp \
        Zone_s.cpp \
        ZoneVarMetadata.cpp \
        ZoneWriterAbstract.cpp \
        ZoneWriterFactory.cpp \
        ZoneWriterFactoryMPI.cpp \
        ZoneWriterFactory.cpp \
        ../boost/algorithm/string.hpp\
        ../boost/assign.hpp\
        ../boost/atomic.hpp\
        ../boost/bind.hpp\
        ../boost/foreach.hpp\
        ../boost/function.hpp\
        ../boost/make_shared.hpp\
        ../boost/ref.hpp\
        ../boost/scoped_array.hpp\
        ../boost/scoped_ptr.hpp\
        ../boost/shared_ptr.hpp\
        ../boost/static_assert.hpp\
        ../boost/tokenizer.hpp\
        ../boost/unordered_map.hpp\
        ../boost/unordered_set.hpp




######################################################################
AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS  = $(pkg_cppflags) -fPIC # any defined package-specific include paths
AM_CFLAGS    =
AM_CXXFLAGS  =
AM_LDFLAGS   = $(libmesh_LDFLAGS)

noinst_LIBRARIES = libteciompi.a

libteciompi_a_SOURCES  = $(pkg_sources)
libteciompi_a_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
libteciompi_a_CXXFLAGS = $(CXXFLAGS_DBG)
libteciompi_a_CFLAGS   = $(CFLAGS_DBG)
