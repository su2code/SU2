#include <tk.h>
#if defined(__WIN32__) || defined(_WIN32)
#   define WIN32_LEAN_AND_MEAN
#   include <windows.h>
#   undef WIN32_LEAN_AND_MEAN
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "nurbs.h"
#include "tkogl.h"
#include "tkoglparse.h"

/*---------------------------------------------------------------------------
 *
 *  The routines below maintain a dynamic array of floats
 *
 *---------------------------------------------------------------------------*/

#define DEFAULTSIZE 8

typedef struct {
   int size;
   int count;
   GLfloat *value;
} FloatArrayStruct, * FloatArray;

static FloatArray 
NewFloatArray () 
{
   FloatArray ptr = malloc (sizeof(FloatArrayStruct));
   assert (ptr != NULL);
   ptr->size = DEFAULTSIZE;
   ptr->count = 0;
   ptr->value = malloc (sizeof (GLfloat)*ptr->size);
   assert (ptr->value != NULL);
   return ptr;
}

static void
DestroyFloatArray (FloatArray ptr)
{
   assert (ptr!=NULL);
   assert (ptr->value != NULL);
   free (ptr->value);
   free (ptr);
}

static void 
AddFloat (FloatArray ptr, GLfloat val)
{
   assert (ptr != NULL);
   assert (ptr->count <= ptr->size);
   if (ptr->count == ptr->size) {
      ptr->size *= 2;	
      ptr->value = realloc (ptr->value, sizeof (GLfloat)*ptr->size);
      assert (ptr->value != NULL);
   }
   ptr->value[ptr->count++] = val;
}

#if 0
static void 
SetFloat (FloatArray ptr, int index, GLfloat val) 
{
   assert (ptr != NULL);
   assert (ptr->count <= ptr->size);
   assert (index <= ptr->count);
   if (index == ptr->size) {
      ptr->size *= 2;	
      ptr->value = realloc (ptr->value, sizeof (GLfloat)*ptr->size);
      assert (ptr->value != NULL);
   }
   if (index == ptr->count) {
      ptr->count++;
   }
   ptr->value [index] = val;
}
#endif

/*--------------------------------------------------------------------------
 *
 *  Main Procedure for generating nurbs surfaces
 *
 *--------------------------------------------------------------------------*/

#define ERRMSG(msg) { Tcl_AppendResult(interp,msg,(char*)NULL); \
		      result = TCL_ERROR; goto done; }

#define ERRMSG2(msg1,msg2) { Tcl_AppendResult(interp,msg1,msg2,(char*)NULL); \
			     result = TCL_ERROR; goto done; }

int
NurbsSurface (Tcl_Interp *interp, int argc, char* argv [])
{
   int result = TCL_OK;
   GLint uOrder = 4;
   GLint vOrder = 4;
   GLenum type = GL_MAP2_VERTEX_3;
   int nCoords = 3;
   FloatArray uKnot = NewFloatArray ();
   FloatArray vKnot = NewFloatArray ();
   FloatArray cPoint = NewFloatArray ();
   GLfloat samplingTolerance = 50.0;
   GLfloat displayMode = GLU_FILL;
   GLfloat culling = GL_FALSE;
   int iarg;
   int dlist = 0;

   for (iarg = 2; iarg < argc; iarg++) {	
      int len = (int)strlen (argv [iarg]);
      if (strncmp (argv [iarg], "-uorder", len) == 0) {
	 int val;
	 iarg++;
	 if (iarg >= argc) ERRMSG ("No value given for -uorder");
	 if (Tcl_GetInt (interp, argv [iarg], &val) != TCL_OK ||
	     val < 2 || val > 8) 
	    ERRMSG2 ("\nInvalid value for -uorder:", argv [iarg]);
	 uOrder = val;
      } 
      else if (strncmp (argv [iarg], "-vorder", len) == 0) {
	 int val;
	 iarg++;
	 if (iarg >= argc) ERRMSG ("No value given for -vorder");
	 if (Tcl_GetInt (interp, argv [iarg], &val) != TCL_OK ||
	     val < 2 || val > 8) 
	    ERRMSG2 ("\nInvalid value for -vorder:", argv [iarg]);
	 vOrder = val;
      }
      else if (strncmp (argv [iarg], "-uknots", len) == 0) {
	 if (uKnot->count != 0) ERRMSG ("uknot values already given");
	 iarg++;
	 while (iarg < argc &&
		!(argv [iarg][0] == '-' && isalpha(argv [iarg][1]))) {
	    double val;
	    if (Tcl_GetDouble (interp, argv [iarg], &val) != TCL_OK) 
	       ERRMSG ("\nError parsing uknot value");
	    if (uKnot->count > 0 &&
		uKnot->value [uKnot->count-1] > val) 
	       ERRMSG ("uknot values not in non-descending order");
	    AddFloat (uKnot, (GLfloat)val);
	    iarg++;
	 }
	 iarg--;
      }	
      else if (strncmp (argv [iarg], "-vknots", len) == 0) {
	 if (vKnot->count != 0) ERRMSG ("vknot values already given");
	 iarg++;
	 while (iarg < argc &&
		!(argv [iarg][0] == '-' && isalpha(argv [iarg][1]))) {
	    double val;
	    if (Tcl_GetDouble (interp, argv [iarg], &val) != TCL_OK) 
	       ERRMSG ("\nError parsing uknot value");
	    if (vKnot->count > 0 &&
		vKnot->value [vKnot->count-1] > val) 
	       ERRMSG ("vknot values not in non-descending order");
	    AddFloat (vKnot, (GLfloat)val);
	    iarg++;
	 }
	 iarg--;
      }	
      else if (strncmp (argv [iarg], "-controlpoints", len) == 0) {
	 if (cPoint->count != 0) ERRMSG ("controlpoint values already given");
	 iarg++;
	 while (iarg < argc &&
		!(argv [iarg][0] == '-' && isalpha(argv [iarg][1]))) {
	    double val;
	    if (Tcl_GetDouble (interp, argv [iarg], &val) != TCL_OK) 
	       ERRMSG ("\nError parsing uknot value");
	    AddFloat (cPoint, (GLfloat)val);
	    iarg++;
	 }
	 iarg--;
      }	      
      else if (strncmp (argv [iarg], "-type", len) == 0) {
	 iarg++;
	 if (iarg >= argc) ERRMSG ("No -type value given");	 
	 if (strcmp (argv [iarg], "map2vertex3") ==0) {
	    type = GL_MAP2_VERTEX_3; nCoords = 3;
	 } else if (strcmp (argv [iarg], "map2vertex4") == 0) {
	    type = GL_MAP2_VERTEX_4; nCoords = 4;
	 } else if (strcmp (argv [iarg], "map2color4") == 0) {
	    type = GL_MAP2_COLOR_4; nCoords = 4;
	 } else if (strcmp (argv [iarg], "map2normal") == 0) {
	    type = GL_MAP2_NORMAL; nCoords = 3; 
	 } else if (strcmp (argv [iarg], "map2texturecoord1") == 0) {
	    type = GL_MAP2_TEXTURE_COORD_1; nCoords = 1; 
	 } else if (strcmp (argv [iarg], "map2texturecoord2") == 0) {
	    type = GL_MAP2_TEXTURE_COORD_2; nCoords = 2;
	 } else if (strcmp (argv [iarg], "map2texturecoord3") == 0) {
	    type = GL_MAP2_TEXTURE_COORD_3; nCoords = 3;
	 } else if (strcmp (argv [iarg], "map2texturecoord4") == 0) {
	    type = GL_MAP2_TEXTURE_COORD_4; nCoords = 4;
	 } else 
	    ERRMSG2 ("not a valid type:", argv [iarg]);
      }
      else if (strncmp (argv [iarg], "-samplingtolerance", len) == 0) {
	 double val;
	 iarg++;
	 if (iarg >= argc) ERRMSG ("No -samplingtolerance value given"); 
	 if (Tcl_GetDouble (interp, argv [iarg], &val) != TCL_OK) 
	    ERRMSG ("\nError parsing sampling tolerance");
	 samplingTolerance = (GLfloat)val;
      }
      else if (strncmp (argv [iarg], "-displaymode", len) == 0) {
	 iarg++;
	 if (iarg >= argc) ERRMSG ("No -displaymode value given");	 
	 if (strcmp (argv [iarg], "fill") == 0) {
	    displayMode = GLU_FILL;
	 } else if (strcmp (argv [iarg], "outlinepolygon") == 0) {
	    displayMode = GLU_OUTLINE_POLYGON;
	 } else if (strcmp (argv [iarg], "outlinepatch") == 0) {
	    displayMode = GLU_OUTLINE_PATCH;
	 } else {
	    ERRMSG2 ("not a valid display mode:", argv [iarg]);
	 }
      }
      else if (strncmp (argv [iarg], "-culling", len) == 0) {
	 int val;
	 iarg++;
	 if (iarg >= argc) ERRMSG ("No -culling value given");	 
	 if (Tcl_GetBoolean (interp, argv [iarg], &val) != TCL_OK) 
	    ERRMSG ("\nError parsing culling value");
	 culling = (GLfloat)val;
      }
      else {
	 ERRMSG2 ("invalid option:", argv [iarg]);
      }
   }

   if (vKnot->count == 0 || uKnot->count == 0 || cPoint->count == 0) 
      ERRMSG ("All of -uknot, -vknot and -cpoint options must be specified");

   /* Now try to guess the remaining arguments and call gluNurbsSurface */
   {
      GLint uKnotCount = uKnot->count;
      GLint vKnotCount = vKnot->count;
      GLint vStride = nCoords;
      GLint uStride = nCoords * (vKnotCount - vOrder);
      static GLUnurbsObj* obj = NULL;
      if (uStride * (uKnotCount - uOrder) != cPoint->count) {
	 char buf [80];
	 sprintf (buf, "%d", uStride * (uKnotCount - uOrder));
	 ERRMSG2 ("Incorrect number of controlpoint coordinates. Expected ",
		  buf);
      }
      /* Theoretically, a nurbs object could be allocated for each
	 invocation of NurbsSurface and then freed after the creation
	 of the display list. However, this produces a segmentation
	 violation on AIX OpenGL 1.0. Thus, only one nurbs object is
	 ever allocated and never freed.
      */
      if (obj == NULL) obj = gluNewNurbsRenderer();
      dlist = glGenLists (1);
      gluNurbsProperty (obj, GLU_SAMPLING_TOLERANCE, samplingTolerance);
      gluNurbsProperty (obj, GLU_DISPLAY_MODE, displayMode);
      gluNurbsProperty (obj, GLU_CULLING, culling);
      glNewList (dlist, GL_COMPILE); 
      gluBeginSurface (obj);
      gluNurbsSurface (obj, uKnotCount, uKnot->value, 
		       vKnotCount, vKnot->value,
		       uStride, vStride, cPoint->value,
		       uOrder, vOrder, type);
      gluEndSurface (obj);
      /* This is never used because of a bug in AIX OpenGL 1.0.
	 gluDeleteNurbsObj (obj);
      */
      glEndList(); 
      glFlush();
   }
      
done:

   DestroyFloatArray (uKnot);
   DestroyFloatArray (vKnot);
   DestroyFloatArray (cPoint);

   if (result == TCL_OK) {
     char tmp[128];
     sprintf (tmp, "%d", dlist);
     Tcl_SetResult(interp, tmp, TCL_VOLATILE);
   }

   return result;
}


