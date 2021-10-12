pkg_cppflags = @SU2_RAPIDJSON_CPPFLAGS@
pkg_cppflags += -I$(srcdir) -I$(srcdir)/error -I$(srcdir)/internal
pkg_sources  = \
    allocators.h\
    cursorstreamwrapper.h\
    document.h\
    encodedstream.h\
    encodings.h\
    filereadstream.h\
    filewritestream.h\
    fwd.h\
    istreamwrapper.h\
    memorybuffer.h\
    memorystream.h\
    ostreamwrapper.h\
    pointer.h\
    prettywriter.h\
    rapidjson.h\
    reader.h\
    schema.h\
    stream.h\
    stringbuffer.h\
    uri.h\
    writer.h\
    error/en.h\
    error/error.h\
    internal/biginteger.h\
    internal/clzll.h\
    internal/diyfp.h\
    internal/dtoa.h\
    internal/ieee754.h\    
    internal/itoa.h\
    internal/meta.h\
    internal/pow10.h\
    internal/regex.h\
    internal/stack.h\    
    internal/strfunc.h\
    internal/strtod.h\
    internal/swap.h


######################################################################
AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS  = $(pkg_cppflags) -fPIC # any defined package-specific include paths
AM_CFLAGS    =
AM_CXXFLAGS  =
AM_LDFLAGS   =

noinst_LIBRARIES = librapidjson.a

librapidjson_a_SOURCES  = $(pkg_sources)
librapidjson_a_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
librapidjson_a_CXXFLAGS = $(CXXFLAGS_DBG)
librapidjson_a_CFLAGS   = $(CFLAGS_DBG)
