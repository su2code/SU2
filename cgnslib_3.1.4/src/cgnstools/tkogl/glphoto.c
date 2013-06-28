#include <tk.h>
#if defined(__WIN32__) || defined(_WIN32)
#   define WIN32_LEAN_AND_MEAN
#   include <windows.h>
#   undef WIN32_LEAN_AND_MEAN
#endif
#include <GL/gl.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "glphoto.h"
#include "tkOGL.h"

/*--------------------------------------------------------------------------
 *
 *  Main Procedure for processing glphoto commands
 *
 *--------------------------------------------------------------------------*/

int
glPhoto (Tcl_Interp *interp, int argc, char* argv [])
{

#define ERRMSG(msg) \
   { Tcl_AppendResult (interp, (msg), (char*) NULL);\
     result = TCL_ERROR;\
     goto done; }

#define ERRMSG2(msg1, msg2) \
   { Tcl_AppendResult (interp, (msg1), (msg2), (char*) NULL);\
     result = TCL_ERROR;\
     goto done; }

   int result = TCL_OK;
   Tk_PhotoHandle handle;
   Tk_PhotoImageBlock block;
   int len;

   if (argc < 1) ERRMSG ("No command specified for glphoto");

   len = strlen (argv[0]);
   if (strncmp (argv [0], "put", len) == 0) {
      if (argc < 2) ERRMSG ("No photo image name specified for 'put' option");
      handle = Tk_FindPhoto (interp,argv [1]);
      if (handle == NULL) ERRMSG2 ("Photo not defined: ", argv [1]);
      if (Tk_PhotoGetImage (handle, &block) != 1) 
	 ERRMSG2 ("Could not get image of photo ", argv [1]);
      printf ("width %d height %d pitch %d pixelsize %d offset %d %d %d\n",
	      block.width, block.height, block.pitch, block.pixelSize,
	      block.offset [0],block.offset [1],block.offset [2]);
	      
      if (block.pixelSize != 3 && block.pixelSize != 4) 
	 ERRMSG ("Image has invalid pixel size");
      glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
      glPixelZoom (1.0, -1.0);
      glRasterPos2i (0, 0);
      glDrawPixels (block.width, block.height, 
		    block.pixelSize == 3 ? GL_RGB : GL_RGBA,
		    GL_UNSIGNED_BYTE, block.pixelPtr);
   }
   else {
      ERRMSG2 ("Invalid glphoto command: ", argv [0]);
   }
done:
   return result;
}




	    
	    




