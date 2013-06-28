#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "tkogl.h"
#include "quadric.h"
#include "tkoglparse.h"

/*--------------------------------------------------------------------------
 *
 *  Main Procedure for generating quadrics
 *
 *--------------------------------------------------------------------------*/

#define ERRMSG(msg) { Tcl_AppendResult(interp,msg,(char*)NULL); \
		      result = TCL_ERROR; goto done; }

#define ERRMSG2(msg1,msg2) { Tcl_AppendResult(interp,msg1,msg2,(char*)NULL); \
			     result = TCL_ERROR; goto done; }

int
Quadric (Tcl_Interp *interp, int argc, char* argv [])
{
   int result = TCL_OK;
   int icoord = 0;
   int iarg;
   GLUquadricObj* obj = gluNewQuadric();
   int dlist = -1;
   GLdouble val [6];

   for (iarg = 2; iarg < argc; iarg++) {
      int len = (int)strlen (argv [iarg]);
      if (strncmp (argv [iarg], "-normals", len) == 0) {
	 iarg++;
	 if (iarg == argc) ERRMSG ("No value for -normals");
	 if (strcmp (argv [iarg], "none") == 0) {
	    gluQuadricNormals (obj, GLU_NONE);
	 }
	 else if (strcmp (argv [iarg], "flat") == 0) {
	    gluQuadricNormals (obj, GLU_FLAT);
	 }
	 else if (strcmp (argv [iarg], "smooth") == 0) {
	    gluQuadricNormals (obj, GLU_SMOOTH);
	 }
	 else ERRMSG2 ("should be 'none', 'flat' or 'smooth':", argv [iarg]);
      }
      else if (strncmp (argv [iarg], "-drawstyle", len) == 0) {
	 iarg++;
	 if (iarg == argc) ERRMSG ("No value for -drawstyle");
	 if (strcmp (argv [iarg], "point") == 0) {
	    gluQuadricDrawStyle (obj, GLU_POINT);
	 }
	 else if (strcmp (argv [iarg], "line") == 0) {
	    gluQuadricDrawStyle (obj, GLU_LINE);
	 }
	 else if (strcmp (argv [iarg], "fill") == 0) {
	    gluQuadricDrawStyle (obj, GLU_FILL);
	 }
	 else if (strcmp (argv [iarg], "silhouette") == 0) {
	    gluQuadricDrawStyle (obj, GLU_SILHOUETTE);
	 }
	 else ERRMSG2 ("should be 'point', 'line', 'silhouette' or 'fill':",
		      argv [iarg]);
      }
      else if (strncmp (argv [iarg], "-orientation", len) == 0) {
	 iarg++;
	 if (iarg == argc) ERRMSG ("No value for -orientation");
	 if (strcmp (argv [iarg], "outside") == 0) {
	    gluQuadricOrientation (obj, GLU_OUTSIDE);
	 }
	 else if (strcmp (argv [iarg], "inside") == 0) {
	    gluQuadricOrientation (obj, GLU_INSIDE);
	 }
	 else ERRMSG2 ("should be 'outside' or 'inside':",
		       argv [iarg]);
      }
      else if (strncmp (argv [iarg], "-texture", len) == 0) {
	 int texture;
	 iarg++;
	 if (iarg == argc) ERRMSG ("No value for -texture");
	 result = Tcl_GetBoolean (interp, argv [iarg], &texture);
	 if (result != TCL_OK) goto done;
	 gluQuadricTexture (obj, (GLboolean) texture);
      }
      else if (strncmp (argv [iarg], "-displaylist", len) == 0) {
	 iarg++;
         if (strcmp (argv [iarg], "none") == 0) {
            dlist = 0;
         }
	 else {
	    result = Tcl_GetInt (interp, argv [iarg], &dlist);
	    if (result != TCL_OK) goto done;
	 }
      }
      else if (argv [iarg][0] == '-' && isalpha (argv [iarg][1])) {
	 ERRMSG2 ("Invalid option: ", argv [iarg])
      }
      else break;
   }
   for (; iarg < argc && icoord < 6; iarg++) {	
      if (Tcl_GetDouble (interp, argv [iarg], &val[icoord]) != TCL_OK) {
	 ERRMSG2 ("\nInvalid value in ", argv [1]);
      }
      icoord++;
   }
   if (dlist == -1) dlist = glGenLists (1);
   if (dlist != 0) glNewList (dlist, GL_COMPILE);
   switch (argv [1][0]) {
      case 'c': { /* Cylinder */
	 if (iarg != argc || icoord > 5) 
	    ERRMSG ("too many values for cylinder");
	 if (icoord < 5) 
	    ERRMSG ("too few values for cylinder");
	 gluCylinder (obj, val [0], val [1], val [2], (GLint) val [3],
		      (GLint) val [4]);
	 break;
      }
      case 'd': { /* Disk */
	 if (iarg != argc || icoord > 4) 
	    ERRMSG ("too many values for disk");
	 if (icoord < 4) 
	    ERRMSG ("too few values for disk");
	 gluDisk (obj, val [0], val [1], (GLint) val [2], (GLint) val [3]);
	 break;
      }
      case 'p': { /* PartialDisk */
	 if (iarg != argc || icoord > 6) 
	    ERRMSG ("too many values for partialdisk");
	 if (icoord < 6) 
	    ERRMSG ("too few values for partialdisk");
	 gluPartialDisk (obj, val [0], val [1], (GLint) val [2],
		      (GLint) val [3], val [4], val [5]);
	 break;
      }
      case 's': { /* Sphere */
	 if (iarg != argc || icoord > 3) 
	    ERRMSG ("too many values for sphere");
	 if (icoord < 3) 
	    ERRMSG ("too few values for sphere");
	 gluSphere (obj, val [0], (GLint) val [1], (GLint) val [2]);
	 break;
      }
   }
	    
done:
   gluDeleteQuadric (obj);
   if (dlist != 0) glEndList(); else return result;

   if (result == TCL_OK) {
     char tmp[128];
     sprintf (tmp, "%d", dlist);
     Tcl_SetResult(interp, tmp, TCL_VOLATILE);
   }
   else {
      glDeleteLists (dlist, 1);
   }
   return result;
}




	    
	    




