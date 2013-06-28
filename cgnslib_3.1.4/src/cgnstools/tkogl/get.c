#include <string.h>
#include <tk.h>
#if defined(__WIN32__) || defined(_WIN32)
#   define WIN32_LEAN_AND_MEAN
#   include <windows.h>
#   undef WIN32_LEAN_AND_MEAN
#endif
#include <GL/gl.h>

typedef enum {
   INTVAL,
   FLOATVAL,
   BOOLEANVAL,
   ISSETVAL
} valType;

typedef struct {
   char * name;
   GLenum code;
   valType type;
   int nVal;
} getStruct;

static getStruct getTable [] = {
   {"currentcolor", GL_CURRENT_COLOR, FLOATVAL, 4},
   {"currentnormal", GL_CURRENT_NORMAL, FLOATVAL, 4},
   {"modelviewmatrix", GL_MODELVIEW_MATRIX, FLOATVAL, 16},
   {"projectionmatrix", GL_PROJECTION_MATRIX, FLOATVAL, 16},
   {"viewport", GL_VIEWPORT, INTVAL, 4}
};


#define ERRMSG(msg) { Tcl_AppendResult(interp,msg,(char*)NULL); \
		      result = TCL_ERROR; goto done; }

#define ERRMSG2(msg1,msg2) { Tcl_AppendResult(interp,msg1,msg2,(char*)NULL); \
			     result = TCL_ERROR; goto done; }

int
GetGlVal (Tcl_Interp *interp, int argc, char* argv [])
{
   int i, j, len;
   char buf [80];
   GLfloat floatVal [16];
   GLint intVal [16];
 
    if (argc != 3) {
        Tcl_AppendResult (interp, "wrong # args", (char*) NULL);
        return TCL_ERROR;
    }
   
    len = (int)strlen (argv [2]);
    for (i = 0; i < sizeof(getTable)/sizeof(getStruct); i++) {
       if (strncmp (argv [2], getTable[i].name, len) == 0) goto found;
    }
    Tcl_AppendResult (interp, "not implemented", (char*) NULL);
    return TCL_ERROR;

found:
	 switch (getTable [i].code) {
	    case GL_CURRENT_COLOR: glGetFloatv (GL_CURRENT_COLOR, floatVal); break;
        case GL_CURRENT_NORMAL: glGetFloatv (GL_CURRENT_NORMAL, floatVal); break;
        case GL_MODELVIEW_MATRIX: glGetFloatv (GL_MODELVIEW_MATRIX, floatVal); break;
        case GL_PROJECTION_MATRIX: glGetFloatv (GL_PROJECTION_MATRIX, floatVal); break;
        case GL_VIEWPORT: glGetIntegerv (GL_VIEWPORT, intVal); break;
     }
     switch (getTable [i].type) {
        case FLOATVAL: {
		    for (j = 0; j < getTable [i].nVal; j++) {
                sprintf (buf, "%g", floatVal [j]);
		        Tcl_AppendElement (interp, buf);
	        }
	        break;
	    }
	    case INTVAL: {
	        for (j = 0; j < getTable [i].nVal; j++) {
		        sprintf (buf, "%d", intVal [j]);
		        Tcl_AppendElement (interp, buf);
	        }
	        break;
	    }
        default: break;
    }
   
    return TCL_OK;
}



