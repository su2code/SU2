/*
 * tkOGL.c --
 *
 *	This module implements "OGLwin" widgets. This is the TCL
 *	interface to OpenGL.
 *
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include "tkogl.h"
#include "tkoglparse.h"
#include "load3ds.h"
#include "tess.h"
#include "quadric.h"
#include "nurbs.h"
#include "get.h"
#include "gencyl.h"
#include "printstr.h"
#include "feedback.h"

#ifndef CONST
# define CONST
#endif

/*
 * A data structure of the following type is kept for each glxwin
 * widget managed by this file:
 */

typedef struct {
    Tk_Window tkwin;		/* Window that embodies the glxwin.  NULL
   		                     * means window has been deleted but
				 * widget record hasn't been cleaned
                                         * up yet. */
    Display *display;		/* X's token for the window's display. */
    Tcl_Interp *interp;		/* Interpreter associated with widget. */

    int width;			/* Width to request for window.  <= 0 means
				 * don't request any size. */
    int height;			/* Height to request for window.  <= 0 means
			           * don't request any size. */

    int absx, absy;                     /* Absolute x and y of window */

    char* context;	                    /* Name of gl window this window will share
                                         * display lists with */

    int flags;			/* Various flags;  see below for
				 * definitions. */

    int doubleBuffer,           /* Various options that control visual alloc. */
        depthSize, 
        stencilSize, 
        alphaSize,
        accumSize;

    int updatePending;		/* Set to 1 when redraw needed */

    int redrawList;		/* Number of redraw display list (or -1) */

    double aspectRatio; 	/* If not 0, represents a fraction width/height
				   which should be maintained when setting
				   the viewport during the event loop */

    Tk_3DBorder bgBorder;	/* Used for drawing background. */

#if defined(__WIN32__) || defined(_WIN32)

    /*
     * Information required to manage the OpenGl window in WIN32 environment
     */

    HPALETTE hPalette;
    HWND hwnd;
    HDC hdc;
    HGLRC hrc;
    TkWinDCState state;
#else
    GLXContext cx;              /* The GL X context */
#endif

} OGLwin;

/*
 * Information used for argv parsing.
 */

static Tk_ConfigSpec configSpecs[] = {
    {TK_CONFIG_PIXELS, "-height", "height", "Height",
       "300" , Tk_Offset(OGLwin, height), 0},

    {TK_CONFIG_PIXELS, "-width", "width", "Width",
       "300" , Tk_Offset(OGLwin, width), 0},

    {TK_CONFIG_STRING, "-context", "context", "Context",
	NULL, Tk_Offset (OGLwin, context), TK_CONFIG_NULL_OK},

    {TK_CONFIG_BOOLEAN, "-doublebuffer", "doublebuffer", "DoubleBuffer",
	"1", Tk_Offset (OGLwin, doubleBuffer), 0},

    {TK_CONFIG_INT, "-depthsize", "depthsize", "DepthSize",
	"16", Tk_Offset (OGLwin, depthSize), 0},

    {TK_CONFIG_INT, "-stencilsize", "stencilsize", "StencilSize",
	"0", Tk_Offset (OGLwin, stencilSize), 0},

    {TK_CONFIG_INT, "-alphasize", "alphasize", "AlphaSize",
	"0", Tk_Offset (OGLwin, alphaSize), 0},

    {TK_CONFIG_INT, "-accumsize", "accumsize", "AccumSize",
	"0", Tk_Offset (OGLwin, accumSize), 0},

    {TK_CONFIG_DOUBLE, "-aspectratio", "aspectratio", "AspectRatio",
	"0", Tk_Offset (OGLwin, aspectRatio), 0},


    {TK_CONFIG_BORDER, "-background", "background", "Background",
	"#d9d9d9", Tk_Offset(OGLwin, bgBorder), 0},


    {TK_CONFIG_END, (char *) NULL, (char *) NULL, (char *) NULL,
       (char *) NULL, 0, 0}
};

/*
 * Forward declarations for procedures defined later in this file:
 */

static int	OGLwinConfigure _ANSI_ARGS_((Tcl_Interp *interp,
			    OGLwin *glxwinPtr, int argc, char **argv,
			    int flags));
static void	OGLwinDestroy _ANSI_ARGS_((char* clientData));

static void	OGLwinEventProc _ANSI_ARGS_((ClientData clientData,
			    XEvent *eventPtr));

static int	OGLwinWidgetCmd _ANSI_ARGS_((ClientData clientData,
			    Tcl_Interp *, int argc, char **argv));


static void 	OGLwinRedraw _ANSI_ARGS_ ((ClientData clientData));

int             OGLwinCmd(ClientData, Tcl_Interp*, int, char**);

static int      UnusedDList ();

static int      FreeDisplayList (int) ;

static void     OGLwinViewport (OGLwin* oglwinPtr);

static void     MakeCurrent (OGLwin* oglwinPtr);

#if defined(__WIN32__) || defined(_WIN32)

    static LONG WINAPI WndProc (HWND, UINT, WPARAM, LPARAM);
    static int WinMakeWindowExist (OGLwin*);
    static void SetDCPixelFormat (OGLwin*);

#else
    static Colormap getColormap _ANSI_ARGS_((Display *dpy,
						 XVisualInfo *vi));
#endif

#define ERRMSG(msg) {\
    Tcl_AppendResult(interp,msg,(char*)NULL); \
    return TCL_ERROR; \
}
#define ERRMSG2(msg1,msg2) {\
    Tcl_AppendResult(interp,msg1,msg2,(char*)NULL); \
    return TCL_ERROR; \
}


#define ARRANGE_REDRAW(glxwinptr) \
       if (!(glxwinptr)->updatePending) {\
	        Tcl_DoWhenIdle ((Tcl_IdleProc *)OGLwinRedraw, \
			 (ClientData)(glxwinptr));\
	          (glxwinPtr)->updatePending = 1;\
       }



#ifdef __WIN32__

/*
 *----------------------------------------------------------------------
 *
 * DllEntryPoint --
 *
 *	This wrapper function is used by Windows to invoke the
 *	initialization code for the DLL.  If we are compiling
 *	with Visual C++, this routine will be renamed to DllMain.
 *	routine.
 *
 * Results:
 *	Returns TRUE;
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */


BOOL APIENTRY
DllEntryPoint(hInst, reason, reserved)
    HINSTANCE hInst;		/* Library instance handle. */
    DWORD reason;	    	/* Reason this function is being called. */
    LPVOID reserved;		/* Not used. */
{
    return TRUE;
}

static WNDCLASS childClass;	/* Tk Window class for child windows. */

static int createCount = 0;

static LONG WINAPI WndProc (HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
/*
 * This is the event proc for the opengl window
 *
 */
{
    LONG result;
    switch (msg) {
        case WM_WINDOWPOSCHANGED:
        case WM_MOVE:
	 break;
        default:
           return TkWinChildProc (hwnd, msg, wParam, lParam);
    }

    result = DefWindowProc(hwnd, msg, wParam, lParam);
    Tcl_ServiceAll();
    return result;
}

#endif

EXPORT(int,Tkogl_Init)(Tcl_Interp* interp)
{
    Tk_Window topLevel;
    topLevel = Tk_MainWindow(interp);

#ifdef __WIN32__

    /*** Fix floating point exception bug ***/
#   ifdef __BORLANDC__
    /* This is how to do it in Borland C */
    _control87(MCW_EM, MCW_EM);
#   else
    /* Don't know how to do it in Visual C++ - Need help here */
#   endif


#else
    /*** make sure OpenGL's GLX extension supported ***/
    if (!glXQueryExtension(Tk_Display(topLevel), NULL, NULL)) {
       ERRMSG ("X server has no OpenGL GLX extension");
    }
#endif

    /*** Initialize the GL function parse tables ***/
    InitHashTables ();

    /*** Create Tkogl main tcl commands ***/

    Tcl_CreateCommand(interp, "OGLwin", (Tcl_CmdProc *)OGLwinCmd,
        (ClientData)topLevel, (Tcl_CmdDeleteProc *)0);

    /*** Init standard extensions ***/
    if (RegisterTkOGLExtension (interp, "nurbssurface", NurbsSurface)!= TCL_OK)
       return TCL_ERROR;

    if (RegisterTkOGLExtension (interp, "load3ds", glLoad3DStudio) != TCL_OK)
       return TCL_ERROR;

    if (RegisterTkOGLExtension (interp, "tesselate", Tesselate) != TCL_OK)
       return TCL_ERROR;

    if (RegisterTkOGLExtension (interp, "cylinder", Quadric) != TCL_OK)
       return TCL_ERROR;

    if (RegisterTkOGLExtension (interp, "disk", Quadric) != TCL_OK)
       return TCL_ERROR;

    if (RegisterTkOGLExtension (interp, "partialdisk", Quadric) != TCL_OK)
       return TCL_ERROR;

    if (RegisterTkOGLExtension (interp, "sphere", Quadric) != TCL_OK)
       return TCL_ERROR;

    if (RegisterTkOGLExtension (interp, "get", GetGlVal) != TCL_OK)
       return TCL_ERROR;

    if (RegisterTkOGLExtension (interp, "gencyl", GenericCylinder) != TCL_OK)
       return TCL_ERROR;

    if (RegisterTkOGLExtension (interp, "printstring", PrintString) != TCL_OK)
       return TCL_ERROR;

    if (RegisterTkOGLExtension (interp, "loadbitmapfont", LoadBitmapFont) != TCL_OK)
       return TCL_ERROR;

    if (RegisterTkOGLExtension (interp, "feedback", Feedback) != TCL_OK)
       return TCL_ERROR;

    return Tcl_PkgProvide(interp, "Tkogl", "1.0");
}



#ifdef __WIN32__


static int WinMakeWindowExist (OGLwin* oglWinPtr)
/*         ==================
 *
 * Forces window to exist so that we can setup MS's version of OpenGL
 */
{
   static char* TkOGLClassName = "TkOGL Class";
   static int TkOGLClassInitted = 0;
   TkWindow *winPtr = (TkWindow *) oglWinPtr->tkwin;
   Display *dpy = Tk_Display (oglWinPtr->tkwin);
   Window parent;
   HWND hwnd, parentWin;
   HANDLE hInstance;
   WNDCLASS TkOGLClass;
   Tcl_HashEntry *hPtr;
   int new_flag;

   /* Destroy window if already exists */
   if (winPtr->window != None) {
      XDestroyWindow(dpy, winPtr->window);
   }

   /* Find parent of window */
   /* Necessary for creation */
   if ((winPtr->parentPtr == NULL) || (winPtr->flags & TK_TOP_LEVEL)) {
      parent = XRootWindow(winPtr->display, winPtr->screenNum);
   }
   else {
      if (winPtr->parentPtr->window == None) {
         Tk_MakeWindowExist((Tk_Window) winPtr->parentPtr);
      }
      parent = winPtr->parentPtr->window;
   }

   /* Create a window class for TkOGL windows if not done this yet */
   parentWin = Tk_GetHWND(parent);
   hInstance = Tk_GetHINSTANCE();
   if (!TkOGLClassInitted) {
       TkOGLClassInitted = 1;
       TkOGLClass.style = CS_HREDRAW | CS_VREDRAW;
       TkOGLClass.cbClsExtra = 0;
       TkOGLClass.cbWndExtra = 4;   /* to save struct OGLwin* */
       TkOGLClass.hInstance = hInstance;
       TkOGLClass.hbrBackground = NULL;
       TkOGLClass.lpszMenuName = NULL;
       TkOGLClass.lpszClassName = TkOGLClassName;
       TkOGLClass.lpfnWndProc = WndProc;
       TkOGLClass.hIcon = NULL;
       TkOGLClass.hCursor = NULL;
       if (!RegisterClass(&TkOGLClass)){
          Tcl_AppendResult(oglWinPtr->interp, "could not register TkOGL class",
    		 (char *) NULL);
          return TCL_ERROR;
       }
   }

   /* Create Window */
   hwnd = CreateWindow(TkOGLClassName, NULL, WS_CHILD | WS_CLIPCHILDREN
                       | WS_CLIPSIBLINGS, 0, 0, oglWinPtr->width, oglWinPtr->height,
                       parentWin, NULL, hInstance, NULL);
   SetWindowLong(hwnd, 0, (LONG) oglWinPtr);
   SetWindowPos(hwnd, HWND_TOP, 0, 0, 0, 0,
  	            SWP_NOACTIVATE | SWP_NOMOVE | SWP_NOSIZE);
   oglWinPtr->hwnd = hwnd;
   oglWinPtr->hdc = GetDC(hwnd);
   

   SetDCPixelFormat (oglWinPtr);
   winPtr->window = Tk_AttachHWND ((Tk_Window) winPtr, hwnd);   

   /* Don't know why the things below are necessary (copied
      from togl code) */

   hPtr = Tcl_CreateHashEntry(&winPtr->dispPtr->winTable,
                              (char *) winPtr->window, &new_flag);
   Tcl_SetHashValue(hPtr, winPtr);
   winPtr->dirtyAtts = 0;
   winPtr->dirtyChanges = 0;
   
   Tk_MapWindow (oglWinPtr->tkwin);
   /* XMapWindow (dpy, Tk_WindowId(oglWinPtr->tkwin)); */
   wglMakeCurrent (oglWinPtr->hdc, oglWinPtr->hrc);

   return TCL_OK;

} /* WinMakeWindowExist */


#endif 


int GetAbsXY (OGLwin *glxwinPtr)
/*  ========
 *
 * Updates the absolute x and y values for the upper left corner of
 * the window. Returns 1 if the recomputed values differ from the previous
 * one, or 0 otherwise
 */
{
   Tk_Window root = glxwinPtr->tkwin;
   int x = Tk_X (root);
   int y = Tk_Y (root);
   int modified;
   do {
      root = Tk_Parent (root);
      x += Tk_X (root);
      y += Tk_Y (root);
   } while (!Tk_IsTopLevel (root));
   modified = (x != glxwinPtr->absx || y != glxwinPtr->absy);
   glxwinPtr->absx = x;
   glxwinPtr->absy = y;
   return modified;
}

/*
 *--------------------------------------------------------------
 *
 * OGLwinCmd --
 *
 *	This procedure is invoked to process the "OGLwin" Tcl
 *	command.  It creates a new "OGLwin" widget.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	A new widget is created and configured.
 *
 *--------------------------------------------------------------
 */

int
OGLwinCmd(clientData, interp, argc, argv)
    ClientData clientData;	/* Main window associated with
				             * interpreter. */
    Tcl_Interp *interp;		/* Current interpreter. */
    int argc;			    /* Number of arguments. */
    char **argv;		    /* Argument strings. */
{
    Tk_Window mainwin = (Tk_Window) clientData;
    OGLwin *glxwinPtr;
    Tk_Window tkwin;

#ifndef __WIN32__
    Colormap cmap;
    XVisualInfo *vi;
    int configuration [50], *confPtr;
#endif

    if (argc < 2) {
    	Tcl_AppendResult(interp, "wrong # args:  should be \"",
    		argv[0], " pathName ?options?\"", (char *) NULL);
    	return TCL_ERROR;
    }

    tkwin = Tk_CreateWindowFromPath(interp, mainwin, argv[1], (char *) NULL);
    if (tkwin == NULL) {
       ERRMSG ("Could not create window");
    }


    /*
     * Allocate and initialize the widget record.
     */

    glxwinPtr = (OGLwin *) ckalloc(sizeof(OGLwin));
    glxwinPtr->tkwin = tkwin;
    glxwinPtr->display = Tk_Display(tkwin);
    glxwinPtr->interp = interp;
    glxwinPtr->width  = 300;
    glxwinPtr->height = 300;
    glxwinPtr->absx = glxwinPtr->absy = 0;
   
    glxwinPtr->context = (char*) NULL;
    glxwinPtr->updatePending = 0;
    glxwinPtr->redrawList = -1;
    glxwinPtr->aspectRatio = 0.0;
    glxwinPtr->doubleBuffer = 1;
    glxwinPtr->depthSize = 16;
    glxwinPtr->alphaSize = 0;
    glxwinPtr->accumSize = 0;
    glxwinPtr->stencilSize = 0;
    glxwinPtr->bgBorder = NULL;
#if defined(__WIN32__) || defined(_WIN32)
    glxwinPtr->hwnd = 0;
    glxwinPtr->hdc = 0;
    glxwinPtr->hrc = 0;
    glxwinPtr->hPalette = 0;
#endif
    
    /* create the widget itself */
    Tk_CreateEventHandler(glxwinPtr->tkwin, 
			  ExposureMask|StructureNotifyMask,
			  OGLwinEventProc,
			  (ClientData) glxwinPtr);

    Tcl_CreateCommand(interp, Tk_PathName(glxwinPtr->tkwin),
        (Tcl_CmdProc *)OGLwinWidgetCmd, (ClientData) glxwinPtr,
        (Tcl_CmdDeleteProc *)0);


    if (OGLwinConfigure(interp, glxwinPtr, argc-2, argv+2, 0) != TCL_OK) {
    	Tk_DestroyWindow(glxwinPtr->tkwin);
	    return TCL_ERROR;
    }

#ifdef __WIN32__
    if (WinMakeWindowExist (glxwinPtr) != TCL_OK) {
        Tcl_AppendResult(interp, "Could not make window exist", (char *) NULL);
    	return TCL_ERROR;
    }

#else
    /*
     * For OpenGL under X we may allocate the OpenGL X context and
     * visual here. First, we need to build a configuration
     * array corresponding to the the widget options.
     */
    confPtr = configuration;
    *confPtr++ = GLX_RGBA;
    if (glxwinPtr->doubleBuffer) {
       *confPtr++ = GLX_DOUBLEBUFFER; 
    }
    if (glxwinPtr->depthSize) {
       *confPtr++ = GLX_DEPTH_SIZE;
       *confPtr++ = glxwinPtr->depthSize;
    }
    if (glxwinPtr->stencilSize) {
       *confPtr++ = GLX_STENCIL_SIZE;
       *confPtr++ = glxwinPtr->stencilSize;
    }
    if (glxwinPtr->accumSize) {
       *confPtr++ = GLX_ACCUM_RED_SIZE;
       *confPtr++ = glxwinPtr->accumSize;
       *confPtr++ = GLX_ACCUM_GREEN_SIZE;
       *confPtr++ = glxwinPtr->accumSize;
       *confPtr++ = GLX_ACCUM_BLUE_SIZE;
       *confPtr++ = glxwinPtr->accumSize;
       if (glxwinPtr->alphaSize) {
    	  *confPtr++ = GLX_ACCUM_ALPHA_SIZE;
	      *confPtr++ = glxwinPtr->accumSize;
       }
    }
    if (glxwinPtr->alphaSize) {
       *confPtr++ = GLX_ALPHA_SIZE;
       *confPtr++ = glxwinPtr->alphaSize;
    }

    /*** find an appropriate visual and a colormap for it ***/
    /* first, try to find a 24 bit visual */
    confPtr [0] = GLX_RED_SIZE; 
    confPtr [1] = 8;
    confPtr [2] = GLX_GREEN_SIZE; 
    confPtr [3] = 8;
    confPtr [4] = GLX_BLUE_SIZE; 
    confPtr [5] = 8;
    confPtr [6] = None;

    vi = glXChooseVisual(glxwinPtr->display,
			 DefaultScreen(glxwinPtr->display), configuration);
    if (vi == NULL) {
       /* try an 8 bit visual */
       confPtr [0] = None;
       vi = glXChooseVisual(glxwinPtr->display,
	    		    DefaultScreen(glxwinPtr->display), configuration);

       if (vi == NULL) {
	        Tk_DestroyWindow(glxwinPtr->tkwin);
	        ERRMSG ("no appropriate RGB visual");
       }
    }
    cmap = getColormap(glxwinPtr->display, vi);

    if (!Tk_SetWindowVisual (tkwin, vi->visual, vi->depth, cmap)) {
       Tk_DestroyWindow(glxwinPtr->tkwin);
       ERRMSG ("Could not set window visual");
    }

    Tk_SetWindowColormap (tkwin, cmap);

    if ((Tk_Parent(tkwin) != NULL) &&
	(Tk_Colormap(tkwin) != Tk_Colormap (Tk_Parent(tkwin)))) {
       TkWmAddToColormapWindows(tkwin);
    } 

    /* See if this window will share display lists with another */
    if (glxwinPtr->context != NULL) {
        Tcl_CmdInfo info;
        if (!Tcl_GetCommandInfo (interp, glxwinPtr->context, &info) ||
	        info.proc != (Tcl_CmdProc *)OGLwinWidgetCmd) {
	        Tcl_AppendResult (interp, "Not a gl window: ", glxwinPtr->context,
			    (char*) NULL);	  
    	    Tk_DestroyWindow(glxwinPtr->tkwin);
	        return TCL_ERROR;
        }
        glxwinPtr->cx = glXCreateContext(Tk_Display(tkwin), vi,
		    			((OGLwin *)(info.clientData))->cx,
			    		/* direct rendering */ GL_TRUE);
     }
     else {
        glxwinPtr->cx = glXCreateContext(Tk_Display(tkwin), vi,
			    	     /* no sharing of display lists */ NULL,
				         /* direct rendering */ GL_TRUE);
     }

    if (glxwinPtr->cx == NULL) {
        Tcl_AppendResult (interp, "could not create rendering context", 
			  (char*) NULL);	  
        Tk_DestroyWindow(glxwinPtr->tkwin);
        return TCL_ERROR;
    }

    Tk_MapWindow (tkwin);

    XSetWMColormapWindows(glxwinPtr->display,
			  Tk_WindowId(tkwin), 
			  &(Tk_WindowId(tkwin)), 1 );

#endif
    MakeCurrent(glxwinPtr);
    OGLwinViewport (glxwinPtr);
    ARRANGE_REDRAW(glxwinPtr);
    GetAbsXY (glxwinPtr);
    
    Tcl_SetResult(interp, Tk_PathName(glxwinPtr->tkwin), TCL_VOLATILE);
    return TCL_OK;
}




/*---------------------------------------------------------------------------
 *
 * OGLwinViewport
 * 
 *    Called when window is resized. If a non-zero aspect ratio was defined
 *    for that window, a corresponding viewport is set by centering it
 *    inside the window.
 *
 *---------------------------------------------------------------------------
 */
static void
OGLwinViewport (OGLwin* oglWinPtr) 
{
   if (oglWinPtr->aspectRatio == 0.0) {
      glViewport (0, 0, oglWinPtr->width, oglWinPtr->height);
   }
   else {
      if (oglWinPtr->width/oglWinPtr->aspectRatio >  oglWinPtr->height) {
	 /* Should reduce the width of the viewport */
	 int adjWidth = (int) (oglWinPtr->height * oglWinPtr->aspectRatio);
	 glViewport ((oglWinPtr->width - adjWidth) / 2, 0, adjWidth, 
		     oglWinPtr->height);
      }
      else {
	 /* Should reduce the height of the viewport */
	 int adjHeight = (int) (oglWinPtr->width/oglWinPtr->aspectRatio);
	 glViewport (0, (oglWinPtr->height - adjHeight) / 2, 
		     oglWinPtr->width, adjHeight);
      }
   }
}

/*---------------------------------------------------------------------------
 *
 * MakeCurrent
 * 
 *    Called whenever it is necessary that the OpenGL context for a
 *    given window is current.
 *
 *---------------------------------------------------------------------------
 */
static void
MakeCurrent (OGLwin* oglwinPtr) 
{
   static OGLwin* previous = NULL;
   if (oglwinPtr == previous) return;
   if (oglwinPtr == NULL) return;
   previous = oglwinPtr;

#ifdef __WIN32__
    wglMakeCurrent (oglwinPtr->hdc, oglwinPtr->hrc);
#else
    glXMakeCurrent(oglwinPtr->display, Tk_WindowId(oglwinPtr->tkwin),
       		  oglwinPtr->cx);
#endif
}



/*---------------------------------------------------------------------------
 *
 * OGLwinRedraw
 * 
 *    Sample redraw routine. Redraws
 *    display list redrawList if defined
 *
 *---------------------------------------------------------------------------
 */
    
static void 
OGLwinRedraw (clientData)
     ClientData clientData;
{
   OGLwin *glxwinPtr = (OGLwin *) clientData;
   Tk_Window tkwin = glxwinPtr->tkwin;

   glxwinPtr->updatePending = 0;

   if (tkwin == NULL || !Tk_IsMapped(tkwin)) {
      return;
   }

#ifdef __WIN32__
    assert (glxwinPtr->hwnd != 0);
    assert (glxwinPtr->hdc != 0);

    MakeCurrent (glxwinPtr);

    if (glxwinPtr->redrawList != -1) {
       glCallList (glxwinPtr->redrawList);
    }

    if (glxwinPtr->doubleBuffer) {
        SwapBuffers (glxwinPtr->hdc);
    }
    else {
        glFlush ();
    }

#else
   MakeCurrent(glxwinPtr);
   
   if (glxwinPtr->redrawList != -1) {
      glCallList (glxwinPtr->redrawList);
   }

   if (glxwinPtr->doubleBuffer) {
      glXSwapBuffers(Tk_Display (glxwinPtr->tkwin),
		     Tk_WindowId(glxwinPtr->tkwin));	
      /* buffer swap does implicit glFlush */
   }
   else {
      glFlush();		
      /* explicit flush for single buffered case */
   }
#endif
}



/*
 *----------------------------------------------------------------------
 * 
 * Display list management.
 *
 *----------------------------------------------------------------------
 */

static int 
UnusedDList ()
{
   /* Returns an integer number corresponding to a free display list
    */
   return glGenLists (1);
}

static int
FreeDisplayList (displayList)
     int displayList;
{
   /* free the given display list. returns 0 if displayList is not in use or
    *  1 if it was correctly freed 
    */
   glDeleteLists (displayList, 1);
   return 1;
}



/*---------------------------------------------------------------------------
 *
 * Management of Extensions to the TKOGL widget
 *
 *---------------------------------------------------------------------------*/

static struct TkOGLExtStruct {
   char* name;
   TkOGLExtProc* proc;
   struct TkOGLExtStruct* next;
} *TkOGLExtTable = NULL;

EXPORT(int,RegisterTkOGLExtension) (interp, extname, extproc)
     Tcl_Interp* interp;
     char* extname;
     TkOGLExtProc* extproc;
{
   struct TkOGLExtStruct* ptr;
   ptr = (struct TkOGLExtStruct*) malloc (sizeof (struct TkOGLExtStruct));
   assert (ptr != NULL);
   ptr->name = (char*) malloc (strlen (extname)+1);
   assert (ptr->name != NULL);
   strcpy (ptr->name, extname);
   ptr->proc = extproc;
   ptr->next = TkOGLExtTable;
   TkOGLExtTable = ptr;
   return TCL_OK;
}


/*---------------------------------------------------------------------------
 *
 * Utility function to format the hit buffer used for selection as a
 * Tcl List
 *
 *---------------------------------------------------------------------------*/

static int
ProcessHits (interp, nhits, selectbuf)
     Tcl_Interp* interp;
     int nhits;
     GLuint * selectbuf;
{
   Tcl_DString hitlist;
   int i;
   unsigned int j;
   GLuint names, *ptr;

   Tcl_DStringInit (&hitlist);
   ptr = selectbuf;
   for (i = 0; i < nhits; i++) {
      char buf [80];
      names = *ptr++;
      Tcl_DStringStartSublist (&hitlist);
      sprintf (buf, "%u", names);
      Tcl_DStringAppendElement (&hitlist, buf);
      sprintf (buf, "%f", ((double) (*ptr++)) / (unsigned int) 0xffffffff );
      Tcl_DStringAppendElement (&hitlist, buf);
      sprintf (buf, "%f", ((double) (*ptr++)) / (unsigned int) 0xffffffff );
      Tcl_DStringAppendElement (&hitlist, buf);
      for (j = 0; j < names; j++) {
	 sprintf (buf, "%d", *ptr++);
	 Tcl_DStringAppendElement (&hitlist, buf);
      }
      Tcl_DStringEndSublist (&hitlist);
   }
   Tcl_DStringResult (interp, &hitlist);
   Tcl_DStringFree (&hitlist);

   return TCL_OK;
}



/*
 *--------------------------------------------------------------
 *
 * OGLwinWidgetCmd --
 *
 *	This procedure is invoked to process the Tcl command
 *	that corresponds to a widget managed by this module.
 *	See the user documentation for details on what it does.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	See the user documentation.
 *
 *--------------------------------------------------------------
 */

static int
OGLwinWidgetCmd(clientData, interp, argc, argv)
    ClientData clientData;		/* Information about glxwin widget. */
    Tcl_Interp *interp;			/* Current interpreter. */
    int argc;				/* Number of arguments. */
    char **argv;			/* Argument strings. */
{
   OGLwin *glxwinPtr = (OGLwin *) clientData;
   Tk_Window tkwin = glxwinPtr->tkwin;
   int result = TCL_OK;
   int length;
   char c;

#ifdef __WIN32__
#   define MAPWINDOW {\
        if (!Tk_IsMapped (tkwin)) Tk_MapWindow (tkwin); \
        assert (glxwinPtr->hdc != 0);\
        MakeCurrent (glxwinPtr);\
    }
#else
#   define MAPWINDOW { \
       if (!Tk_IsMapped (tkwin)) Tk_MapWindow (tkwin); \
       MakeCurrent(glxwinPtr);\
    }
#endif

   if (argc < 2) {
      Tcl_AppendResult(interp, "wrong # args: should be \"",
		       argv[0], " option ?arg arg ...?\"", (char *) NULL);
      return TCL_ERROR;
   }
   Tcl_Preserve((ClientData) glxwinPtr);
   c = argv[1][0];
   length = (int)strlen(argv[1]);

   if ((c == 'c') && (strncmp(argv[1], "configure", length) == 0)
       && (length >= 2)) {
        /* Standard "configure" tk command */
        if (argc == 2) {
           result = Tk_ConfigureInfo(interp, glxwinPtr->tkwin,
		         configSpecs, (char *) glxwinPtr, (char *) NULL, 0);
        }
        else if (argc == 3) {
	   result = Tk_ConfigureInfo(interp, glxwinPtr->tkwin,
			 configSpecs, (char *) glxwinPtr, argv[2], 0);
        }
        else {
	   result = OGLwinConfigure(interp, glxwinPtr, argc-2, argv+2,
				   TK_CONFIG_ARGV_ONLY);
      }
   }
   else if ((c == 'c') && (strncmp(argv[1], "cget", length) == 0)
	    && (length >= 2)) {
       /* Standard "cget" tk command */
       if (argc != 3) {
          Tcl_AppendResult(interp, "wrong # args: should be \"",
		    argv[0], " cget option\"",(char *) NULL);
	  goto error;
       }
       result = Tk_ConfigureValue(interp, glxwinPtr->tkwin, configSpecs,
		(char *) glxwinPtr, argv[2], 0);
   }
   else if ((c == 'm') && (length >= 4) &&
	(strncmp(argv[1], "mainlist", length) == 0)) {
        /* Establishes the main display list for the widget */
        int narg;
        MAPWINDOW;
        if (glxwinPtr->redrawList != -1) {
            glDeleteLists (glxwinPtr->redrawList, 1);
        } 
        else {
            glxwinPtr->redrawList = UnusedDList();
            if (glxwinPtr->redrawList == -1) {
                Tcl_AppendResult (interp, "Can't find unused display list",
          		       (char*) NULL);
                goto error;
            }
        }
        argc -= 2;
        argv += 2;
        glNewList (glxwinPtr->redrawList, GL_COMPILE);
        while (result == TCL_OK && argc > 0) {
	        result = ParseGLFunc (interp, argc, argv, &narg);
	        argc -= narg;
	        argv += narg;
        }
        glEndList();
        ARRANGE_REDRAW (glxwinPtr);
   }
   else if ((c == 'n') && (strncmp(argv[1], "newlist", length) == 0)) {
        /* Creates a new display list and returns its number */
        int newlist, narg;
        MAPWINDOW;
        if (argc > 2 && isdigit (argv [2][0])) {
	        if (Tcl_GetInt (interp, argv [2], &newlist) != TCL_OK) goto error;
	        argc -= 3;
	        argv += 3;
        } 
        else {	    
	        newlist = UnusedDList();
	        argc -= 2;
	        argv += 2;
        }
        glNewList (newlist, GL_COMPILE);
        while (result == TCL_OK && argc > 0) {
	        result = ParseGLFunc (interp, argc, argv, &narg);
	        argc -= narg;
	        argv += narg;
        }
        glEndList();
        if (result == TCL_OK) {
	  char tmp[128];
	  sprintf (tmp, "%d", newlist);
	  Tcl_SetResult(interp, tmp, TCL_VOLATILE);
	}
   }
   else if ((c == 'e') && (strncmp(argv[1], "eval", length) == 0)) {
      /* sends the gl commands directly */
      int narg;
      MAPWINDOW;
      argc -= 2;
      argv += 2;
      while (result == TCL_OK && argc > 0) {
    	  result = ParseGLFunc (interp, argc, argv, &narg);
	      argc -= narg;
	      argv += narg;
      }
      glFlush ();
   }
   else if ((c == 's') && (strncmp(argv[1], "select", length) == 0)) {
      /* sets rendermode to select, issue the commands, reset the rendermode
       * to render and returns the hits as a tcl list */
      int narg, nvals, nhits;
      GLuint *selectbuf;
      if (argc < 3) {
	        Tcl_AppendResult (interp, "no size for hit buffer specified", NULL);
	        goto error;
      }
      if (Tcl_GetInt (interp, argv [2], &nvals) != TCL_OK ||
	        nvals <= 0) {
	        Tcl_AppendResult (interp, "invalid hit buffer size", argv [2], NULL);
	        goto error;
      }
      selectbuf = (GLuint*) malloc (sizeof (GLuint)*nvals);
      assert (selectbuf != NULL);
      MAPWINDOW;      
      argc -= 3;
      argv += 3;
      glSelectBuffer (nvals, selectbuf);
      (void) glRenderMode (GL_SELECT);
      while (result == TCL_OK && argc > 0) {
	        result = ParseGLFunc (interp, argc, argv, &narg);
	        argc -= narg;
	        argv += narg;
      }
      glFlush ();
      nhits = glRenderMode (GL_RENDER);
      result = ProcessHits (interp, nhits, selectbuf);
      free (selectbuf);
   }
   else if ((c == 'd') && (strncmp(argv[1], "deletelist", length) == 0)) {
      /* Free a given display list */
      int listNumber;
      if (argc != 3) {
    	  Tcl_AppendResult (interp, "no display list number specified",
			    (char*) NULL);
	      goto error;
      }
      if (Tcl_GetInt (interp, argv [2], &listNumber) != TCL_OK) goto error;
      MAPWINDOW;
      if (!FreeDisplayList (listNumber)) {
    	  Tcl_AppendResult (interp, "not a display list: ", argv [2],
			    (char*) NULL);
	      goto error;
      }
   }
   else if ((c == 'p' && strncmp(argv[1], "project", length) == 0) ||
	    (c == 'u' && strncmp(argv[1], "unproject", length) == 0)) {
      /* Implementation of gluProject (gluUnProject) utility function */
      /* Given 3 world(window) coordinates, returns the corresponding */
      /* 3 window(world) coordinates */
      GLdouble x, y, z, modelmatrix[16], projectmatrix[16];
      GLint viewport [4], retval;
      if (argc != 5) {		
    	 Tcl_AppendResult (interp, argv [1], ": x y z coordinates expected",
			   (char*) NULL);
	     goto error;
      }
      if (Tcl_GetDouble (interp, argv [2], &x) != TCL_OK ||
	  Tcl_GetDouble (interp, argv [3], &y) != TCL_OK ||
	  Tcl_GetDouble (interp, argv [4], &z) != TCL_OK) goto error;
      MAPWINDOW;
      glGetDoublev (GL_MODELVIEW_MATRIX, modelmatrix);
      glGetDoublev (GL_PROJECTION_MATRIX, projectmatrix);
      glGetIntegerv (GL_VIEWPORT, viewport);
      if (c == 'p') {
    	 retval = gluProject (x, y, z, modelmatrix, projectmatrix, 
			      viewport, &x, &y, &z);
	 y = viewport [3] - y;
      } 
      else {
    	 y = viewport [3] - y;
	 retval = gluUnProject (x, y, z, modelmatrix, projectmatrix, 
				viewport, &x, &y, &z);	 
      }
      if (retval) {
	char tmp[128];
	sprintf (tmp, "%f %f %f", x, y, z);
	Tcl_SetResult(interp, tmp, TCL_VOLATILE);
      }
   } 
   else if ((c == 'r') && (strncmp(argv[1], "redraw", length) == 0)) {
      /* perform the redraw */
      ARRANGE_REDRAW(glxwinPtr);
   } 
   else {
      /* handle extensions */
      struct TkOGLExtStruct * extPtr = TkOGLExtTable;
      while (extPtr != NULL) {
	      if (strncmp(argv[1], extPtr->name, length) == 0) {
	            MAPWINDOW;
	            result = (*extPtr->proc) (interp, argc, argv);
	            break;
	      }
	      extPtr = extPtr->next;
      }
      if (extPtr == NULL) {
	        Tcl_AppendResult (interp, "bad command: ", argv[1], (char*) NULL);
	        goto error;
      }
   }

   Tcl_Release((ClientData) glxwinPtr);
   return result;

   error:
   Tcl_Release((ClientData) glxwinPtr);
   return TCL_ERROR;
}


/*
 *----------------------------------------------------------------------
 *
 * OGLwinConfigure --
 *
 *	This procedure is called to process an argv/argc list in
 *	conjunction with the Tk option database to configure (or
 *	reconfigure) a glxwin widget.
 *
 * Results:
 *	The return value is a standard Tcl result.  If TCL_ERROR is
 *	returned, then interp->result contains an error message.
 *
 * Side effects:
 *	Configuration information, such as colors, border width,
 *	etc. get set for glxwinPtr;  old resources get freed,
 *	if there were any.
 *
 *----------------------------------------------------------------------
 */

static int
OGLwinConfigure(interp, glxwinPtr, argc, argv, flags)
    Tcl_Interp *interp;			/* Used for error reporting. */
    OGLwin *glxwinPtr;			/* Information about widget. */
    int argc;				/* Number of valid entries in argv. */
    char **argv;			/* Arguments. */
    int flags;				/* Flags to pass to
					 * Tk_ConfigureWidget. */
{
    if (Tk_ConfigureWidget(interp, glxwinPtr->tkwin, configSpecs,
	    argc, (CONST char **)argv, (char *) glxwinPtr, flags) != TCL_OK) {
	return TCL_ERROR;
    }

    /*
     * Register the desired geometry for the window.  Then arrange for
     * the window to be redisplayed.
     */

/*
    Tk_SetWindowBackground(glxwinPtr->tkwin,
	    Tk_3DBorderColor(glxwinPtr->bgBorder)->pixel);
*/
    if ((glxwinPtr->width > 0) || (glxwinPtr->height > 0)) {

       Tk_GeometryRequest(glxwinPtr->tkwin, glxwinPtr->width,
			  glxwinPtr->height);
    }
    return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * OGLwinEventProc --
 *
 *	This procedure is invoked by the Tk dispatcher for various
 *	events on glxwins.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	When the window gets deleted, internal structures get
 *	cleaned up.  When it gets exposed, it is redisplayed.
 *
 *--------------------------------------------------------------
 */

static void
OGLwinEventProc(clientData, eventPtr)
    ClientData clientData;	/* Information about window. */
    XEvent *eventPtr;		/* Information about event. */
{
    OGLwin *glxwinPtr = (OGLwin *) clientData;
    
#ifdef __WIN32__

    assert (glxwinPtr->hwnd != 0);

    if (eventPtr->type == Expose) {
        if (GetAbsXY (glxwinPtr)) {
           /* This Kludge is necessary since Windows Tk does not pass move events
            * to child windows and hence our OpenGL window might be redrawn
            * in its old position
            */
           Tk_ResizeWindow (glxwinPtr->tkwin, glxwinPtr->width, glxwinPtr->height);
        } 
        ARRANGE_REDRAW(glxwinPtr);
    } 
    else if (eventPtr->type == ConfigureNotify) {
        GLsizei glnWidth = eventPtr->xconfigure.width; 
        GLsizei glnHeight = eventPtr->xconfigure.height; 
        glxwinPtr->width = glnWidth;
        glxwinPtr->height = glnHeight;
        OGLwinViewport (glxwinPtr);
    } 
#else
    if (eventPtr->type == Expose) {
       ARRANGE_REDRAW(glxwinPtr);
    } 
    else if (eventPtr->type == ConfigureNotify) {
       glxwinPtr->width = eventPtr->xconfigure.width;
       glxwinPtr->height = eventPtr->xconfigure.height;
       OGLwinViewport (glxwinPtr);
       ARRANGE_REDRAW(glxwinPtr);
    }
#endif /* WIN32 */
    else if (eventPtr->type == DestroyNotify) {
    	Tcl_DeleteCommand(glxwinPtr->interp, Tk_PathName(glxwinPtr->tkwin));
	    glxwinPtr->tkwin = NULL;
	    if (glxwinPtr->updatePending) {
    	   Tcl_CancelIdleCall(OGLwinRedraw, (ClientData) glxwinPtr);
        }
	    Tcl_EventuallyFree((ClientData) glxwinPtr, OGLwinDestroy);
    }

}


/*
 *----------------------------------------------------------------------
 *
 * OGLwinDestroy --
 *
 *	This procedure is invoked by Tcl_EventuallyFree or Tcl_Release
 *	to clean up the internal structure of a glxwin at a safe time
 *	(when no-one is using it anymore).
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Everything associated with the glxwin is freed up.
 *
 *----------------------------------------------------------------------
 */

static void
OGLwinDestroy(char* clientData)
{
    OGLwin *glxwinPtr = (OGLwin *) clientData;

#if defined(__WIN32__) || defined(_WIN32)
    if (glxwinPtr->hrc != 0) {
        wglMakeCurrent (NULL, NULL);
        wglDeleteContext (glxwinPtr->hrc);
        MakeCurrent (NULL);
        ReleaseDC (glxwinPtr->hwnd, glxwinPtr->hdc);
        if (glxwinPtr->hPalette != 0) {
            DeleteObject (glxwinPtr->hPalette);
        }
    }
#endif
    Tk_FreeOptions(configSpecs, (char *) glxwinPtr, glxwinPtr->display, 0);
    ckfree((char *) glxwinPtr);
}



#ifdef __WIN32__
/*
 *----------------------------------------------------------------------
 *
 *  SetDCPixelFormat sets the pixel format for a device context in
 *  preparation for creating a rendering context.
 *
 *  Input parameters:
 *    hdc = Device context handle
 *
 *  Returns:
 *    Nothing
 *
 *----------------------------------------------------------------------
 */
 
static void SetDCPixelFormat (OGLwin* glxwinPtr)
{
    HANDLE hHeap;
    int nColors, i;
    LPLOGPALETTE lpPalette;
    HDC hdc = glxwinPtr->hdc;

    BYTE byRedMask, byGreenMask, byBlueMask;

    PIXELFORMATDESCRIPTOR pfd = {
        sizeof (PIXELFORMATDESCRIPTOR),             // Size of this structure
        1,                                          // Version number
        PFD_DRAW_TO_WINDOW |                        // Flags
        PFD_SUPPORT_OPENGL,
        PFD_TYPE_RGBA,                              // RGBA pixel values
        24,                                         // 24-bit color
        0, 0, 0, 0, 0, 0,                           // Don't care about these
        0, 0,                                       // No alpha buffer
        0, 0, 0, 0, 0,                              // No accumulation buffer
        32,                                         // 32-bit depth buffer
        0,                                          // No stencil buffer
        0,                                          // No auxiliary buffers
        PFD_MAIN_PLANE,                             // Layer type
        0,                                          // Reserved (must be 0)
        0, 0, 0                                     // No layer masks
    };

    int nPixelFormat = 0;
    XVisualInfo VisInf;
    XVisualInfo *visinfo = &VisInf;
    TkWinColormap *cmap = (TkWinColormap *) ckalloc(sizeof(TkWinColormap));
    
    /* Just for portability, define the simplest visinfo */
    visinfo->visual = DefaultVisual(glxwinPtr->display, 
                                    DefaultScreen(glxwinPtr->display));
    visinfo->depth = visinfo->visual->bits_per_rgb;

    /* Set pfd fields according to the capabilities needed */
    if (glxwinPtr->doubleBuffer) {
       pfd.dwFlags |= PFD_DOUBLEBUFFER; 
    }
    else {
       pfd.dwFlags &= ~PFD_DOUBLEBUFFER;
    }

    pfd.cDepthBits = glxwinPtr->depthSize;
    pfd.cStencilBits = glxwinPtr->stencilSize;
    pfd.cAccumBits = glxwinPtr->accumSize;
    pfd.cAlphaBits = glxwinPtr->alphaSize;

    nPixelFormat = ChoosePixelFormat (hdc, &pfd);
    SetPixelFormat (hdc, nPixelFormat, &pfd);
    DescribePixelFormat (hdc, nPixelFormat, sizeof (PIXELFORMATDESCRIPTOR),
                         &pfd);

    glxwinPtr->depthSize = pfd.cDepthBits;
    glxwinPtr->stencilSize = pfd.cStencilBits;
    glxwinPtr->accumSize = pfd.cAccumBits;
    glxwinPtr->alphaSize = pfd.cAlphaBits;
    glxwinPtr->doubleBuffer = (pfd.dwFlags & PFD_DOUBLEBUFFER) != 0;

    if (pfd.dwFlags & PFD_NEED_PALETTE) {
        nColors = 1 << pfd.cColorBits;
        hHeap = GetProcessHeap ();

        (LPLOGPALETTE) lpPalette = HeapAlloc (hHeap, 0,
            sizeof (LOGPALETTE) + (nColors * sizeof (PALETTEENTRY)));
            
        lpPalette->palVersion = 0x300;
        lpPalette->palNumEntries = nColors;

        byRedMask = (1 << pfd.cRedBits) - 1;
        byGreenMask = (1 << pfd.cGreenBits) - 1;
        byBlueMask = (1 << pfd.cBlueBits) - 1;

        for (i=0; i<nColors; i++) {
            lpPalette->palPalEntry[i].peRed =
                (((i >> pfd.cRedShift) & byRedMask) * 255) / byRedMask;
            lpPalette->palPalEntry[i].peGreen =
                (((i >> pfd.cGreenShift) & byGreenMask) * 255) / byGreenMask;
            lpPalette->palPalEntry[i].peBlue =
                (((i >> pfd.cBlueShift) & byBlueMask) * 255) / byBlueMask;
            lpPalette->palPalEntry[i].peFlags = 0;
        }

        glxwinPtr->hPalette = CreatePalette (lpPalette);
        HeapFree (hHeap, 0, lpPalette);

        if (glxwinPtr->hPalette != NULL) {
            SelectPalette (hdc, glxwinPtr->hPalette, FALSE);
            RealizePalette (hdc);
        }
        
        cmap->palette = glxwinPtr->hPalette;
        cmap->size = nColors;
        cmap->stale = 0;
        
        /* Since this is a private colormap of a fix size, we do not need
           a valid hash table, but a dummy one */
        Tcl_InitHashTable(&cmap->refCounts, TCL_ONE_WORD_KEYS);
    }
    else {
        cmap = (TkWinColormap *) DefaultColormap(glxwinPtr->display,
                               DefaultScreen(glxwinPtr->display));
    }

    /* Now, allocate the OpenGL context */
    glxwinPtr->hrc = wglCreateContext (glxwinPtr->hdc);

    /* See if this window will share display lists with another */
    if (glxwinPtr->context != NULL) {
        Tcl_CmdInfo info;
        if (!Tcl_GetCommandInfo (glxwinPtr->interp, glxwinPtr->context,
            &info) || info.proc != (Tcl_CmdProc *)OGLwinWidgetCmd) {
            Tcl_ResetResult (glxwinPtr->interp);
            Tcl_AppendResult (glxwinPtr->interp,
                "Not a gl window: ", glxwinPtr->context, (char*) NULL);
	        Tcl_BackgroundError (glxwinPtr->interp);
        }
        else if (!wglShareLists(((OGLwin *)(info.clientData))->hrc,
                  glxwinPtr->hrc)) {
            Tcl_ResetResult (glxwinPtr->interp);
            Tcl_AppendResult (glxwinPtr->interp,
                "Cannot share lists with: ", glxwinPtr->context, (char*) NULL);
	        Tcl_BackgroundError (glxwinPtr->interp);
        }
    }
    /* Make sure Tk knows how to switch palettes */
    Tk_SetWindowVisual (glxwinPtr->tkwin, visinfo->visual,
                        visinfo->depth, (Colormap) cmap);
}
#else
/*---------------------------------------------------------------------------
 *
 * The following function allocates an X colormap appropriate for using
 * with OpenGL. This was taken from the example program 'glxdino' by
 * Mark Kilgard
 */

Colormap
getColormap (Display *dpy, XVisualInfo * vi)
{
    Status          status;
    XStandardColormap *standardCmaps;
    Colormap        cmap;
    int             i, numCmaps;

    /* if no standard colormap but TrueColor, just make an unshared one */
    status = XmuLookupStandardColormap(dpy, vi->screen, vi->visualid,
        vi->depth, XA_RGB_DEFAULT_MAP, /* replace */ False, /* retain */ True);
    if (status == 1 && (vi->class == DirectColor || vi->depth == 8)) {
	status = XGetRGBColormaps(dpy, RootWindow(dpy, vi->screen),
			     &standardCmaps, &numCmaps, XA_RGB_DEFAULT_MAP);
	if (status == 1)
	    for (i = 0; i < numCmaps; i++)
		if (standardCmaps[i].visualid == vi->visualid) {
		    cmap = standardCmaps[i].colormap;
		    XFree(standardCmaps);
		    return cmap;
		}
    }
    cmap = XCreateColormap(dpy, RootWindow(dpy, vi->screen),
        vi->visual, AllocNone);
    return cmap;
}
#endif




