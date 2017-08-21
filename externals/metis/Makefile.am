pkg_cppflags  = @SU2_METIS_CPPFLAGS@
pkg_cppflags += -I$(srcdir)/include  -I$(srcdir)/GKlib -I$(srcdir)/libmetis
pkg_sources = \
	include/metis.h \
	GKlib/b64.c \
	GKlib/blas.c \
	GKlib/csr.c \
	GKlib/error.c \
	GKlib/evaluate.c \
	GKlib/fkvkselect.c \
	GKlib/fs.c \
	GKlib/getopt.c \
	GKlib/gk_arch.h \
	GKlib/gk_defs.h \
	GKlib/gk_externs.h \
	GKlib/gk_getopt.h \
	GKlib/gk_macros.h \
	GKlib/gk_mkblas.h \
	GKlib/gk_mkmemory.h \
	GKlib/gk_mkpqueue.h \
	GKlib/gk_mkpqueue2.h \
	GKlib/gk_mkrandom.h \
	GKlib/gk_mksort.h \
	GKlib/gk_mkutils.h \
	GKlib/gk_proto.h \
	GKlib/gk_struct.h \
	GKlib/gk_types.h \
	GKlib/GKlib.h \
	GKlib/gkregex.c \
	GKlib/gkregex.h \
	GKlib/graph.c \
	GKlib/htable.c \
	GKlib/io.c \
	GKlib/itemsets.c \
	GKlib/mcore.c \
	GKlib/memory.c \
	GKlib/ms_inttypes.h \
	GKlib/ms_stat.h \
	GKlib/ms_stdint.h \
	GKlib/omp.c \
	GKlib/pdb.c \
	GKlib/pqueue.c \
	GKlib/random.c \
	GKlib/rw.c \
	GKlib/seq.c \
	GKlib/sort.c \
	GKlib/string.c \
	GKlib/timers.c \
	GKlib/tokenizer.c \
	GKlib/util.c \
	libmetis/auxapi.c \
	libmetis/balance.c \
	libmetis/bucketsort.c \
	libmetis/checkgraph.c \
	libmetis/coarsen.c \
	libmetis/compress.c \
	libmetis/contig.c \
	libmetis/debug.c \
	libmetis/defs.h \
	libmetis/fm.c \
	libmetis/fortran.c \
	libmetis/frename.c \
	libmetis/gklib.c \
	libmetis/gklib_defs.h \
	libmetis/gklib_rename.h \
	libmetis/graph.c \
	libmetis/initpart.c \
	libmetis/kmetis.c \
	libmetis/kwayfm.c \
	libmetis/kwayrefine.c \
	libmetis/macros.h \
	libmetis/mcutil.c \
	libmetis/mesh.c \
	libmetis/meshpart.c \
	libmetis/metislib.h \
	libmetis/minconn.c \
	libmetis/mincover.c \
	libmetis/mmd.c \
	libmetis/ometis.c \
	libmetis/options.c \
	libmetis/parmetis.c \
	libmetis/pmetis.c \
	libmetis/proto.h \
	libmetis/refine.c \
	libmetis/rename.h \
	libmetis/separator.c \
	libmetis/sfm.c \
	libmetis/srefine.c \
	libmetis/stat.c \
	libmetis/stdheaders.h \
	libmetis/struct.h \
	libmetis/timing.c \
	libmetis/util.c \
	libmetis/wspace.c

######################################################################
AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS  = $(pkg_cppflags) -fPIC # any defined package-specific include paths
AM_CFLAGS    =
AM_CXXFLAGS  =
AM_LDFLAGS   =

noinst_LIBRARIES = libmetis.a

libmetis_a_SOURCES  = $(pkg_sources)
libmetis_a_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
libmetis_a_CXXFLAGS = $(CXXFLAGS_DBG)
libmetis_a_CFLAGS   = $(CFLAGS_DBG)
