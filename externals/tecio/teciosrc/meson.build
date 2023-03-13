tecio_cpp_flags = ['-fpermissive', '-DUSEENUM', '-DTHREED', '-DMAKEARCHIVE', '-DNO_ASSERTS', '-DNO_THIRD_PARTY_LIBS']

if (host_machine.system() == 'darwin')
  tecio_cpp_flags += ['-DDARWIN', '-DMAC64']
endif

if (host_machine.system() == 'linux')
  tecio_cpp_flags += ['-DLINUX', '-DLINUX64']
endif

if (host_machine.system() == 'windows')
  tecio_cpp_flags += ['-DMSWIN']
endif

if build_machine.system() != 'windows' and meson.get_compiler('cpp').get_id() != 'intel'
  desired_warnings = ['-Wno-misleading-indentation',
                      '-Wno-uninitialized',
                      '-Wno-placement-new',
                      '-Wno-pedantic',
                      '-Wno-error=deprecated-copy']
  foreach flag : desired_warnings
    if meson.get_compiler('cpp').has_argument(flag)
      tecio_cpp_flags += flag
    endif
  endforeach
endif

tecio_include = include_directories('../', './')
tecio = static_library('tecio',
           'ClassicZoneWriterAbstract.cpp',
           'szcombine.cpp',
           'FileSystem.cpp',
           'SZLOrderedZoneHeaderWriter.cpp',
           'ZoneWriterFactory.cpp',
           'importSzPltFile.cpp',
           'FaceNeighborGeneratorAbstract.cpp',
           'TranslatedString.cpp',
           'PartitionTecUtilDecorator.cpp',
           'tecio.cpp',
           'TecioData.cpp',
           'checkPercentDone.cpp',
           'Zone_s.cpp',
           'dataio4.cpp',
           'DataSetWriter.cpp',
           'ClassicZoneVariableWriter.cpp',
           'ClassicZoneFaceNeighborWriter.cpp',
           'UnicodeStringUtils.cpp',
           'ClassicZoneHeaderWriter.cpp',
           'SZLOrderedPartitionedZoneWriter.cpp',
           'SZLFEPartitionedZoneWriter.cpp',
           'writeValueArray.cpp',
           'TecioPLT.cpp',
           'FECellSubzoneCompressor.cpp',
           'fileStuff.cpp',
           'FileStreamReader.cpp',
           'strlist.cpp',
           'SZLOrderedPartitionWriter.cpp',
           'ZoneHeaderWriterAbstract.cpp',
           'TecioTecUtil.cpp',
           'ZoneWriterAbstract.cpp',
           'geom2.cpp',
           'NodeToElemMap_s.cpp',
           'FileIOStream.cpp',
           'SZLFEPartitionedZoneHeaderWriter.cpp',
           'ClassicFEZoneConnectivityWriter.cpp',
           'strutil.cpp',
           'q_msg.cpp',
           'auxdata.cpp',
           'MinMaxTree.cpp',
           'FieldData.cpp',
           'set.cpp',
           'SZLOrderedPartitionedZoneHeaderWriter.cpp',
           'AsciiOutputInfo.cpp',
           'FieldData_s.cpp',
           'dataio.cpp',
           'IntervalTree.cpp',
           'ClassicFEZoneFaceNeighborGenerator.cpp',
           'ClassicOrderedZoneFaceNeighborGenerator.cpp',
           'OrthogonalBisection.cpp',
           'ORBFESubzonePartitioner.cpp',
           'datautil.cpp',
           'IJKZoneInfo.cpp',
           'IJKSubzoneInfo.cpp',
           'Scanner.cpp',
           'ClassicFEZoneWriter.cpp',
           'NodeMap.cpp',
           'filestream.cpp',
           'SZLOrderedZoneWriter.cpp',
           'SZLFEPartitionWriter.cpp',
           'Geom_s.cpp',
           'SZLFEZoneWriter.cpp',
           'ZoneVarMetadata.cpp',
           'dataset.cpp',
           'SZLFEZoneHeaderWriter.cpp',
           'FileStreamWriter.cpp',
           'ClassicOrderedZoneWriter.cpp',
           'ZoneInfoCache.cpp',
           'arrlist.cpp',
           'dataset0.cpp',
           'readValueArray.cpp',
           'TecioSZL.cpp',
           'exportSubzonePlt.cpp',
           install : false,
	   cpp_args: [default_warning_flags, tecio_cpp_flags],
	   include_directories: tecio_include)

tecio_dep = declare_dependency(link_with: tecio,
                               include_directories: tecio_include)
