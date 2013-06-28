#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "tkogl.h"
#include "tess.h"

/*--------------------------------------------------------------------------
 *
 *  Main Procedure for generating tesselations
 *
 *--------------------------------------------------------------------------*/

static Tcl_Interp * globalinterp;
static int globalresult;

static void TessError (GLenum errno) {
   Tcl_AppendResult (globalinterp, "Tesselation error: ", 
		     gluErrorString(errno), " ", (char*) NULL);
   globalresult = TCL_ERROR;
}

#if defined(__WIN32__) || defined(_WIN32)
    typedef void (CALLBACK* TessCallback) ();
#else
    typedef void (* TessCallback) ();
#endif
    
int
Tesselate (Tcl_Interp *interp, int argc, char* argv [])
{
   int result = TCL_OK;
   int icoord = 0;
   int edgeflags = 1;
   int iarg;
   GLUtriangulatorObj* obj;
   int dlist = -1;
   GLdouble coord [3];
   GLfloat* vtx;
   GLfloat* vtxptr;

   globalinterp = interp;
   globalresult = TCL_OK;

   for (iarg = 2; iarg < argc; iarg++) {	
      int len = (int)strlen (argv [iarg]);
      if (strncmp (argv [iarg], "-displaylist", len) == 0) {
	 iarg++;
	 if (iarg == argc) {
            Tcl_AppendResult (interp, "not enough arguments", (char*)NULL);
            return TCL_ERROR;
         }
         if (strcmp (argv [iarg], "none") == 0) {
            dlist = 0;
         } 
         else if (Tcl_GetInt (interp, argv [iarg], &dlist) != TCL_OK) {
	    Tcl_AppendResult (interp,
	           "\nError parsing display list number", (char*) NULL);
	    return TCL_ERROR;
	 }
      } 
      else if (strncmp (argv [iarg], "-noedgeflags", len) == 0) {
	 edgeflags = 0;
      }
      else break;
   }

   if (argc - iarg < 9) {
      Tcl_AppendResult (interp,
	      "Not enough vertices", (char*) NULL);
      return TCL_ERROR;
   }

   obj = gluNewTess();
   vtx = (GLfloat*) malloc (sizeof (GLfloat) * (argc - iarg));
   vtxptr = vtx;

   assert (vtx != NULL);

   gluTessCallback(obj, GLU_BEGIN, glBegin);
   gluTessCallback(obj, GLU_VERTEX, glVertex3fv); 
   gluTessCallback(obj, GLU_END, glEnd);
   gluTessCallback(obj, GLU_ERROR, (TessCallback) TessError);
   if (edgeflags) gluTessCallback (obj, GLU_EDGE_FLAG, (TessCallback) glEdgeFlag);
   if (dlist == -1) dlist = glGenLists (1);   
   if (dlist != 0) glNewList (dlist, GL_COMPILE);
   gluBeginPolygon (obj);

   for (; iarg < argc; iarg++) {	
      int len = (int)strlen (argv [iarg]);
      if (strncmp (argv [iarg], "-contour", len) == 0) {
	 gluNextContour (obj, GLU_UNKNOWN);
      }
      else {
	 if (Tcl_GetDouble (interp, argv [iarg], &coord[icoord]) != TCL_OK) {
	    Tcl_AppendResult (interp,
	      "\nError parsing tesselation vertex coord", (char*) NULL);
	    result = TCL_ERROR;
	    break;
	 }
	 else {
	    icoord = (icoord+1)%3;
	    if (icoord == 0) {
	       *(vtxptr) = (GLfloat)coord [0];
	       *(vtxptr+1) = (GLfloat)coord [1];
	       *(vtxptr+2) = (GLfloat)coord [2];
	       gluTessVertex (obj, coord, vtxptr);
	       vtxptr += 3;
	    }
	 }
      }
   }
   
   gluEndPolygon (obj);
   gluDeleteTess (obj);
   free (vtx);
   if (dlist != 0) glEndList(); 

   if (result != TCL_OK || globalresult != TCL_OK) {
      if (dlist != 0) glDeleteLists (dlist, 1);
      return TCL_ERROR;
   }
   
   if (dlist != 0) {
     char tmp[128];
     sprintf (tmp, "%d", dlist);
     Tcl_SetResult(interp, tmp, TCL_VOLATILE);
   }
   return TCL_OK;
}




	    
	    




