#ifndef _TK_OGL
#define _TK_OGL

#include <tk.h>
#if defined(__WIN32__) || defined(_WIN32)
#   define WIN32_LEAN_AND_MEAN
#   include <windows.h>
#   undef WIN32_LEAN_AND_MEAN
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <assert.h>
#if defined(__WIN32__) || defined(_WIN32)
#   include "tkWinInt.h"
#   if defined(_MSC_VER)
#	    define EXPORT(a,b) __declspec(dllexport) a b
#	    define DllEntryPoint DllMain
#   else
#	    if defined(__BORLANDC__)
#               include <float.h>
#	        define EXPORT(a,b) a _export b
#	    else
#	        define EXPORT(a,b) a b
#	    endif
#   endif
#else
#   define EXPORT(a,b) a b
#   include <GL/glx.h>
#   include <X11/Xatom.h>		/* for XA_RGB_DEFAULT_MAP atom */
#   include <X11/Xmu/StdCmap.h>	/* for XmuLookupStandardColormap() */
#endif

typedef int (TkOGLExtProc) (Tcl_Interp* interp, int argc, char** argv);

EXPORT (int,RegisterTkOGLExtension) (Tcl_Interp* interp, char* extname, TkOGLExtProc* extproc);

#endif /* _TK_OGL */


