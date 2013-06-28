#include "tkogl.h"
#include "tkoglparse.h"
#include "feedback.h"
#include <string.h>
#include <stdlib.h>

#if defined(__WIN32__) || defined(_WIN32)
#   ifdef strncasecmp
#      undef strncasecmp
#   endif
#   if defined(__BORLANDC__) || defined(__CYGWIN__)
#      define strncasecmp(a,b,c) strnicmp(a,b,c)
#   else
#      define strncasecmp(a,b,c) _strnicmp(a,b,c)
#   endif
#endif

typedef int PrintVtx (Tcl_DString *, GLfloat* );

static int
ProcessBuffer (Tcl_Interp *interp, int size, GLfloat* buffer, PrintVtx * proc) 
{
   Tcl_DString list;
   int i, j, n;
   char buf [80];

   Tcl_DStringInit (&list);
   for (i = 0; i < size; i++) {
      int nvtx = 0;
      if (*buffer == GL_POINT_TOKEN) {
	 Tcl_DStringAppendElement (&list, "-point");
	 nvtx = 1;
      }
      else if (*buffer == GL_BITMAP_TOKEN) {
	 Tcl_DStringAppendElement (&list, "-bitmap");
	 nvtx = 1;
      }
      else if (*buffer == GL_DRAW_PIXEL_TOKEN) {
	 Tcl_DStringAppendElement (&list, "-drawpixel");
	 nvtx = 1;
      }
      else if (*buffer == GL_COPY_PIXEL_TOKEN) {
	 Tcl_DStringAppendElement (&list, "-copypixel");
	 nvtx = 1;
      }
      else if (*buffer == GL_LINE_TOKEN) {
	 Tcl_DStringAppendElement (&list, "-line");
	 nvtx = 2;
      }
      else if (*buffer == GL_LINE_RESET_TOKEN) {
	 Tcl_DStringAppendElement (&list, "-linereset");
	 nvtx = 2;
      }
      else if (*buffer == GL_POLYGON_TOKEN) {
	 Tcl_DStringAppendElement (&list, "-polygon");
	 buffer++;
	 i++;
	 nvtx = (int) *buffer;
      }
      else if (*buffer == GL_PASS_THROUGH_TOKEN) {
	 Tcl_DStringAppendElement (&list, "-passthrough");
	 buffer++;
	 i++;
	 sprintf (buf, "%d", (int) *buffer);
	 Tcl_DStringAppendElement (&list, buf);
      }
      if (nvtx > 0) {
	 ++buffer;
	 for (j = 0; j < nvtx; j++) {
	    n = (*proc) (&list, buffer);
	    buffer += n;
	    i += n;
	 }
      }
   }
   Tcl_DStringResult (interp, &list);
   Tcl_DStringFree (&list);

   return TCL_OK;
}
  
static int 
PrintVtx2D (Tcl_DString *list, GLfloat* buffer)
{
   char buf [80];
   int i;
   Tcl_DStringAppendElement (list, "-vertex");
   for (i = 0; i < 2; i++) {
      sprintf (buf, "%g", *buffer++);
      Tcl_DStringAppendElement (list, buf);
   }
   return 2;
}

static int 
PrintVtx3D (Tcl_DString *list, GLfloat* buffer)
{
   char buf [80];
   int i;
   Tcl_DStringAppendElement (list, "-vertex");
   for (i = 0; i < 3; i++) {
      sprintf (buf, "%g", *buffer++);
      Tcl_DStringAppendElement (list, buf);
   }
   return 3;
}

static int 
PrintVtx3DColor (Tcl_DString *list, GLfloat* buffer)
{
   char buf [80];
   int i;
   Tcl_DStringAppendElement (list, "-vertex");
   for (i = 0; i < 3; i++) {
      sprintf (buf, "%g", *buffer++);
      Tcl_DStringAppendElement (list, buf);
   }
   Tcl_DStringAppendElement (list, "-color");
   for (i = 0; i < 4; i++) {
      sprintf (buf, "%g", *buffer++);
      Tcl_DStringAppendElement (list, buf);
   }
   return 7;
}

static int 
PrintVtx3DColorTexture (Tcl_DString *list, GLfloat* buffer)
{
   char buf [80];
   int i;
   Tcl_DStringAppendElement (list, "-vertex");
   for (i = 0; i < 3; i++) {
      sprintf (buf, "%g", *buffer++);
      Tcl_DStringAppendElement (list, buf);
   }
   Tcl_DStringAppendElement (list, "-color");
   for (i = 0; i < 4; i++) {
      sprintf (buf, "%g", *buffer++);
      Tcl_DStringAppendElement (list, buf);
   }
   Tcl_DStringAppendElement (list, "-texture");
   for (i = 0; i < 4; i++) {
      sprintf (buf, "%g", *buffer++);
      Tcl_DStringAppendElement (list, buf);
   }
   return 11;
}

static int 
PrintVtx4DColorTexture (Tcl_DString *list, GLfloat* buffer)
{
   char buf [80];
   int i;
   Tcl_DStringAppendElement (list, "-vertex");
   for (i = 0; i < 4; i++) {
      sprintf (buf, "%g", *buffer++);
      Tcl_DStringAppendElement (list, buf);
   }
   Tcl_DStringAppendElement (list, "-color");
   for (i = 0; i < 4; i++) {
      sprintf (buf, "%g", *buffer++);
      Tcl_DStringAppendElement (list, buf);
   }
   Tcl_DStringAppendElement (list, "-texture");
   for (i = 0; i < 4; i++) {
      sprintf (buf, "%g", *buffer++);
      Tcl_DStringAppendElement (list, buf);
   }
   return 12;
}

int
Feedback (Tcl_Interp *interp, int argc, char* argv [])
{
   int result = TCL_OK;
   int narg;
   int len;
   int size;
   PrintVtx * proc;
   GLenum type;
   GLfloat *buffer;

   if (argc < 4) {
      Tcl_AppendResult (interp, "Not enough args", (char*) NULL);
      return TCL_ERROR;
   } 

   /* Parse 'size' */
   if (Tcl_GetInt (interp, argv [2], &size) != TCL_OK || size <= 0) {
      Tcl_AppendResult (interp, "invalid hit buffer size", argv [2], NULL);
      return TCL_ERROR;
   }

   /* Parse 'type' */
   len = (int)strlen (argv [3]);
   if (strncasecmp (argv [3], "2d", len) == 0) {
      type = GL_2D;
      proc = PrintVtx2D;
   } 
   else if (strncasecmp (argv [3], "3d", len) == 0) {	
      type = GL_3D;
      proc = PrintVtx3D;
   } 
   else if (strncasecmp (argv [3], "3dColor", len) == 0) {	
      type = GL_3D_COLOR;
      proc = PrintVtx3DColor;
   } 
   else if (strncasecmp (argv [3], "3dColorTexture", len) == 0) {	
      type = GL_3D_COLOR_TEXTURE;
      proc = PrintVtx3DColorTexture;
   } 
   else if (strncasecmp (argv [3], "4dColorTexture", len) == 0) {	
      type = GL_4D_COLOR_TEXTURE;
      proc = PrintVtx4DColorTexture;
   } 
   else {
      Tcl_AppendResult (interp, "Unknown type:", argv [3]);
      return TCL_ERROR;
   }
      
   /* Allocate feedback 'buffer' */
   buffer = (GLfloat*) malloc (sizeof (GLfloat)*size);
   assert (buffer != NULL);

   /* Render in feedback mode */
   glFeedbackBuffer (size, type, buffer);
   (void) glRenderMode (GL_FEEDBACK);
   argc -= 4; argv += 4;
   while (result == TCL_OK && argc > 0) {
      result = ParseGLFunc (interp, argc, argv, &narg);
      argc -= narg;
      argv += narg;
   }
   glFlush ();

   /* Generate feedback string */
   size = glRenderMode (GL_RENDER);
   if (result == TCL_OK) {
      result = ProcessBuffer (interp, size, buffer, proc);
   }
   free (buffer);
   return result;
}


