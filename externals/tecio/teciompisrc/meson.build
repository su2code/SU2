tec_cxx_flags = ['-fpermissive', '-DUSEENUM','-DTECIOMPI', '-DTHREED', '-DMAKEARCHIVE', '-DNO_ASSERTS', '-DNO_THIRD_PARTY_LIBS', '-DTP_PROJECT_USES_BOOST', '-DBOOST_ALL_NO_LIB']

if (host_machine.system() == 'darwin')
  tec_cxx_flags += ['-DDARWIN', '-DMAC64']
endif

if (host_machine.system() == 'linux')
  tec_cxx_flags += ['-DLINUX', '-DLINUX64']
endif

if (host_machine.system() == 'windows')
  tec_cxx_flags += ['-DMSWIN']
endif

if build_machine.system() != 'windows' and meson.get_compiler('cpp').get_id() != 'intel'
  desired_warnings = ['-Wno-misleading-indentation',
                      '-Wno-uninitialized',
                      '-Wno-placement-new',
                      '-Wno-pedantic',
                      '-Wno-error=class-memaccess',
                      '-Wno-class-memaccess',
                      '-Wno-error=deprecated-copy',
                      '-Wno-deprecated-copy',
                      '-Wno-error=cast-function-type',
                      '-Wno-cast-function-type']
  foreach flag : desired_warnings
    if meson.get_compiler('cpp').has_argument(flag)
      tec_cxx_flags += flag
    endif
  endforeach
endif

teciompi_include = include_directories(['../', './'])

teciompi = static_library('teciompi',
           'ClassicZoneWriterAbstract.cpp',
           'SZLOrderedPartitionedZoneWriterMPI.cpp',
           'szcombine.cpp',
           'FileSystem.cpp',
           'MPIUtil.cpp',
           'SZLOrderedZoneHeaderWriter.cpp',
           'ZoneWriterFactory.cpp',
           'importSzPltFile.cpp',
           'FaceNeighborGeneratorAbstract.cpp',
           'PartitionTecUtilDecorator.cpp',
           'ZoneWriterFactoryMPI.cpp',
           'tecio.cpp',
           'TecioData.cpp',
           'MPIFileIOStream.cpp',
           'checkPercentDone.cpp',
           'Zone_s.cpp',
           'DataSetWriter.cpp',
           'ClassicZoneVariableWriter.cpp',
           'ClassicZoneFaceNeighborWriter.cpp',
           'UnicodeStringUtils.cpp',
           'ClassicZoneHeaderWriter.cpp',
           'SZLOrderedPartitionedZoneWriter.cpp',
           'SZLFEPartitionedZoneWriter.cpp',
           'MPICommunicator.cpp',
           'writeValueArray.cpp',
           'FECellSubzoneCompressor.cpp',
           'fileStuff.cpp',
           'MPICommunicationCache.cpp',
           'FileStreamReader.cpp',
           'SZLFEPartitionedZoneWriterMPI.cpp',
           'SZLOrderedPartitionWriter.cpp',
           'MPIFileWriter.cpp',
           'ZoneHeaderWriterAbstract.cpp',
           'TecioTecUtil.cpp',
           'ZoneWriterAbstract.cpp',
           'NodeToElemMap_s.cpp',
           'DataSetWriterMPI.cpp',
           'FileIOStream.cpp',
           'SZLFEPartitionedZoneHeaderWriter.cpp',
           'MPIFileReader.cpp',
           'ClassicFEZoneConnectivityWriter.cpp',
           'MinMaxTree.cpp',
           'FieldData.cpp',
           'SZLOrderedPartitionedZoneHeaderWriter.cpp',
           'AsciiOutputInfo.cpp',
           'FieldData_s.cpp',
           'IntervalTree.cpp',
           'ClassicFEZoneFaceNeighborGenerator.cpp',
           'mpiDatatype.cpp',
           'ClassicOrderedZoneFaceNeighborGenerator.cpp',
           'OrthogonalBisection.cpp',
           'ORBFESubzonePartitioner.cpp',
           'IJKZoneInfo.cpp',
           'IJKSubzoneInfo.cpp',
           'Scanner.cpp',
           'ClassicFEZoneWriter.cpp',
           'NodeMap.cpp',
           'SZLOrderedZoneWriter.cpp',
           'SZLFEPartitionWriter.cpp',
           'SZLFEZoneWriter.cpp',
           'ZoneVarMetadata.cpp',
           'MPINonBlockingCommunicationCollection.cpp',
           'SZLFEZoneHeaderWriter.cpp',
           'FileStreamWriter.cpp',
           'ClassicOrderedZoneWriter.cpp',
           'ZoneInfoCache.cpp',
           'readValueArray.cpp',
           'TecioSZL.cpp',
           'exportSubzonePlt.cpp',
           install : false,
	   dependencies: [mpi_dep], 
	   cpp_args: [default_warning_flags, tec_cxx_flags],
	   include_directories: 
	   teciompi_include)

teciompi_dep = declare_dependency(link_with: teciompi,
                                  include_directories: teciompi_include)
