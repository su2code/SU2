pkg_cppflags = @SU2_CGNS_CPPFLAGS@
pkg_cppflags += -I$(srcdir) -I$(srcdir)/adf
pkg_sources  = \
    cgnslib.h \
    cgns_header.h \
    cgns_io.h \
    cgnstypes.h \
    cgnsconfig.h \
    fortran_macros.h \
    adf/ADF.h \
    adf/ADF_fbind.h \
    adf/ADF_internals.h \
    cgns_error.c \
    cgns_internals.c \
    cgns_io.c \
    cgnslib.c \
    adf/ADF_interface.c \
    adf/ADF_internals.c


######################################################################
AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS  = $(pkg_cppflags) -fPIC # any defined package-specific include paths
AM_CFLAGS    =
AM_CXXFLAGS  =
AM_LDFLAGS   =

noinst_LIBRARIES = libcgns.a

libcgns_a_SOURCES  = $(pkg_sources)
libcgns_a_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
libcgns_a_CXXFLAGS = $(CXXFLAGS_DBG)
libcgns_a_CFLAGS   = $(CFLAGS_DBG)
