metis_c_args = ['-D_FILE_OFFSET_BITS=64', '-DNDEBUG', '-DNDEBUG2', '-DHAVE_GETLINE', '-DUSE_GKREGEX']

metis_default_warnings = []
if build_machine.system() != 'windows'
  metis_default_warnings += ['-Wno-implicit-function-declaration']
  if meson.get_compiler('cpp').get_id() != 'intel'
    desired_warnings = ['-Wno-unused-result',
                        '-Wno-unused-parameter',
                        '-Wno-unused-variable',
                        '-Wno-unused-but-set-variable',
                        '-Wno-unknown-pragmas',
                        '-Wno-sign-compare',
                        '-Wno-clobbered',
                        '-Wno-empty-body',
                        '-Wno-unused-label',
                        '-Wno-misleading-indentation',
                        '-Wno-maybe-uninitialized',
                        '-Wno-pedantic',
                        '-Wno-error=absolute-value',
                        '-Wno-absolute-value']
    foreach flag : desired_warnings
      if meson.get_compiler('cpp').has_argument(flag)
        metis_default_warnings += flag
      endif
    endforeach
  endif
endif

metis_src = ['GKlib/pqueue.c',
           'GKlib/blas.c',
           'GKlib/fkvkselect.c',
           'GKlib/csr.c',
           'GKlib/b64.c',
           'GKlib/memory.c',
           'GKlib/rw.c',
           'GKlib/util.c',
           'GKlib/timers.c',
           'GKlib/pdb.c',
           'GKlib/random.c',
           'GKlib/fs.c',
           'GKlib/gkregex.c',
           'GKlib/getopt.c',
           'GKlib/htable.c',
           'GKlib/mcore.c',
           'GKlib/sort.c',
           'GKlib/tokenizer.c',
           'GKlib/string.c',
           'GKlib/io.c',
           'GKlib/evaluate.c',
           'GKlib/seq.c',
           'GKlib/omp.c',
           'GKlib/itemsets.c',
           'GKlib/graph.c',
           'GKlib/error.c',
           'libmetis/frename.c',
           'libmetis/timing.c',
           'libmetis/meshpart.c',
           'libmetis/minconn.c',
           'libmetis/mesh.c',
           'libmetis/kwayfm.c',
           'libmetis/bucketsort.c',
           'libmetis/balance.c',
           'libmetis/contig.c',
           'libmetis/mmd.c',
           'libmetis/fm.c',
           'libmetis/util.c',
           'libmetis/refine.c',
           'libmetis/checkgraph.c',
           'libmetis/compress.c',
           'libmetis/parmetis.c',
           'libmetis/pmetis.c',
           'libmetis/stat.c',
           'libmetis/options.c',
           'libmetis/srefine.c',
           'libmetis/kwayrefine.c',
           'libmetis/sfm.c',
           'libmetis/separator.c',
           'libmetis/coarsen.c',
           'libmetis/kmetis.c',
           'libmetis/fortran.c',
           'libmetis/graph.c',
           'libmetis/mincover.c',
           'libmetis/auxapi.c',
           'libmetis/wspace.c',
           'libmetis/initpart.c',
           'libmetis/mcutil.c',
           'libmetis/debug.c',
           'libmetis/ometis.c',
           'libmetis/gklib.c']

metis_include = include_directories(['include','GKlib','libmetis'])
metis_lib     = static_library('metis', 
                               metis_src,
			       include_directories : metis_include,
			       c_args: metis_c_args + metis_default_warnings)

metis_dep = declare_dependency(link_with : metis_lib, 
                               include_directories : metis_include)

