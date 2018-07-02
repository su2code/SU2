# include paths - note that parmetis needs metis!
pkg_cppflags  = @SU2_METIS_CPPFLAGS@ @SU2_PARMETIS_CPPFLAGS@
pkg_cppflags += -I$(top_srcdir)/externals/metis/GKlib \
                -I$(top_srcdir)/externals/metis/include \
                -I$(srcdir)/include \
                -I$(srcdir)/libparmetis
pkg_sources   = \
	include/parmetis.h \
	libparmetis/akwayfm.c \
	libparmetis/ametis.c \
	libparmetis/balancemylink.c \
	libparmetis/comm.c \
	libparmetis/csrmatch.c \
	libparmetis/ctrl.c \
	libparmetis/debug.c \
	libparmetis/defs.h \
	libparmetis/diffutil.c \
	libparmetis/frename.c \
	libparmetis/gkmetis.c \
	libparmetis/gkmpi.c \
	libparmetis/graph.c \
	libparmetis/initbalance.c \
	libparmetis/initmsection.c \
	libparmetis/initpart.c \
	libparmetis/kmetis.c \
	libparmetis/kwayrefine.c \
	libparmetis/macros.h \
	libparmetis/match.c \
	libparmetis/mdiffusion.c \
	libparmetis/mesh.c \
	libparmetis/mmetis.c \
	libparmetis/move.c \
	libparmetis/msetup.c \
	libparmetis/node_refine.c \
	libparmetis/ometis.c \
	libparmetis/parmetislib.h \
	libparmetis/proto.h \
	libparmetis/pspases.c \
	libparmetis/redomylink.c \
	libparmetis/remap.c \
	libparmetis/rename.h \
	libparmetis/renumber.c \
	libparmetis/rmetis.c \
	libparmetis/selectq.c \
	libparmetis/serial.c \
	libparmetis/stat.c \
	libparmetis/struct.h \
	libparmetis/timer.c \
	libparmetis/util.c \
	libparmetis/wave.c \
	libparmetis/weird.c \
	libparmetis/wspace.c \
	libparmetis/xyzpart.c

######################################################################
AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS  = $(pkg_cppflags) -fPIC # any defined package-specific include paths
AM_CFLAGS    =
AM_CXXFLAGS  =
AM_LDFLAGS   =

noinst_LIBRARIES = libparmetis.a

libparmetis_a_SOURCES  = $(pkg_sources)
libparmetis_a_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
libparmetis_a_CXXFLAGS = $(CXXFLAGS_DBG)
libparmetis_a_CFLAGS   = $(CFLAGS_DBG)
