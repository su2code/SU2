/*
 * winMain.c --
 *
 *	Main entry point for wish and other Tk-based applications.
 *
 * Copyright (c) 1995-1997 Sun Microsystems, Inc.
 * Copyright (c) 1998-1999 by Scriptics Corporation.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 */

#include <tk.h>
#include "tkWinInt.h"
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef WIN32_LEAN_AND_MEAN
#include <malloc.h>
#include <locale.h>

extern int CGIOtcl_Init _ANSI_ARGS_((Tcl_Interp *interp));
#ifdef USE_HTMLHELP
extern int WinHtml_Init _ANSI_ARGS_((Tcl_Interp *interp));
#endif

/*
 * Forward declarations for procedures defined later in this file:
 */

static void setargv _ANSI_ARGS_((int *argcPtr, char ***argvPtr));
static void WishPanic _ANSI_ARGS_(TCL_VARARGS(char *,format));
static int Proc_LoadIcon _ANSI_ARGS_((ClientData data,
    Tcl_Interp *interp, int argc, char **argv));

static HINSTANCE myInstance;

/*
 *----------------------------------------------------------------------
 *
 * WinMain --
 *
 *	Main entry point from Windows.
 *
 * Results:
 *	Returns false if initialization fails, otherwise it never
 *	returns.
 *
 * Side effects:
 *	Just about anything, since from here we call arbitrary Tcl code.
 *
 *----------------------------------------------------------------------
 */

int APIENTRY
WinMain(HINSTANCE hInstance,
        HINSTANCE hPrevInstance,
        LPSTR lpszCmdLine,
        int nCmdShow)
{
    char **argv;
    int argc;
    char buffer[MAX_PATH+1];
    char *p;

    Tcl_SetPanicProc((Tcl_PanicProc *)WishPanic);

    /*
     * Increase the application queue size from default value of 8.
     * At the default value, cross application SendMessage of WM_KILLFOCUS
     * will fail because the handler will not be able to do a PostMessage!
     * This is only needed for Windows 3.x, since NT dynamically expands
     * the queue.
     */

    SetMessageQueue(64);

    /*
     * Set up the default locale to be standard "C" locale so parsing
     * is performed correctly.
     */

    setlocale(LC_ALL, "C");
    setargv(&argc, &argv);

    /*
     * Replace argv[0] with full pathname of executable, and forward
     * slashes substituted for backslashes.
     */

    GetModuleFileName(NULL, buffer, sizeof(buffer));
    argv[0] = buffer;
    for (p = buffer; *p != '\0'; p++) {
	if (*p == '\\') {
	    *p = '/';
	}
    }
    myInstance = hInstance;

    Tk_Main(argc, argv, Tcl_AppInit);
    return 1;
}

/*
 *----------------------------------------------------------------------
 *
 * Tcl_AppInit --
 *
 *	This procedure performs application-specific initialization.
 *	Most applications, especially those that incorporate additional
 *	packages, will have their own version of this procedure.
 *
 * Results:
 *	Returns a standard Tcl completion code, and leaves an error
 *	message in the interp's result if an error occurs.
 *
 * Side effects:
 *	Depends on the startup script.
 *
 *----------------------------------------------------------------------
 */

int
Tcl_AppInit(Tcl_Interp *interp)
{
    if (Tcl_Init(interp) == TCL_ERROR) {
	goto error;
    }
    if (Tk_Init(interp) == TCL_ERROR) {
	goto error;
    }
    Tcl_StaticPackage(interp, "Tk", Tk_Init, Tk_SafeInit);

    if (CGIOtcl_Init(interp) == TCL_ERROR)
        goto error;
#ifdef USE_HTMLHELP
    if (WinHtml_Init(interp) == TCL_ERROR)
        goto error;
#endif
    Tcl_CreateCommand (interp, "LoadIcon", (Tcl_CmdProc *)Proc_LoadIcon,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);

    Tcl_SetVar(interp, "tcl_rcFileName", "~/wishrc.tcl", TCL_GLOBAL_ONLY);
    return TCL_OK;

error:
    MessageBeep(MB_ICONEXCLAMATION);
    MessageBox(NULL, Tcl_GetStringResult(interp), "Error in Wish",
	    MB_ICONSTOP | MB_OK | MB_TASKMODAL | MB_SETFOREGROUND);
    ExitProcess(1);
    /* we won't reach this, but we need the return */
    return TCL_ERROR;
}

/*
 *----------------------------------------------------------------------
 *
 * WishPanic --
 *
 *	Display a message and exit.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Exits the program.
 *
 *----------------------------------------------------------------------
 */

void
WishPanic TCL_VARARGS_DEF(char *,arg1)
{
    va_list argList;
    char buf[1024];
    char *format;

    format = TCL_VARARGS_START(char *,arg1,argList);
    vsprintf(buf, format, argList);

    MessageBeep(MB_ICONEXCLAMATION);
    MessageBox(NULL, buf, "Fatal Error in Wish",
	    MB_ICONSTOP | MB_OK | MB_TASKMODAL | MB_SETFOREGROUND);
#ifdef _MSC_VER
    DebugBreak();
#endif
    ExitProcess(1);
}

/*
 *-------------------------------------------------------------------------
 *
 * setargv --
 *
 *	Parse the Windows command line string into argc/argv.  Done here
 *	because we don't trust the builtin argument parser in crt0.
 *	Windows applications are responsible for breaking their command
 *	line into arguments.
 *
 *	2N backslashes + quote -> N backslashes + begin quoted string
 *	2N + 1 backslashes + quote -> literal
 *	N backslashes + non-quote -> literal
 *	quote + quote in a quoted string -> single quote
 *	quote + quote not in quoted string -> empty string
 *	quote -> begin quoted string
 *
 * Results:
 *	Fills argcPtr with the number of arguments and argvPtr with the
 *	array of arguments.
 *
 * Side effects:
 *	Memory allocated.
 *
 *--------------------------------------------------------------------------
 */

static void
setargv(int *argcPtr, char ***argvPtr)
{
    char *cmdLine, *p, *arg, *argSpace;
    char **argv;
    int argc, size, inquote, copy, slashes;

    cmdLine = GetCommandLine();	/* INTL: BUG */

    /*
     * Precompute an overly pessimistic guess at the number of arguments
     * in the command line by counting non-space spans.
     */

    size = 2;
    for (p = cmdLine; *p != '\0'; p++) {
	if ((*p == ' ') || (*p == '\t')) {	/* INTL: ISO space. */
	    size++;
	    while ((*p == ' ') || (*p == '\t')) { /* INTL: ISO space. */
		p++;
	    }
	    if (*p == '\0') {
		break;
	    }
	}
    }
    argSpace = (char *) Tcl_Alloc(
	    (unsigned) (size * sizeof(char *) + strlen(cmdLine) + 1));
    argv = (char **) argSpace;
    argSpace += size * sizeof(char *);
    size--;

    p = cmdLine;
    for (argc = 0; argc < size; argc++) {
	argv[argc] = arg = argSpace;
	while ((*p == ' ') || (*p == '\t')) {	/* INTL: ISO space. */
	    p++;
	}
	if (*p == '\0') {
	    break;
	}

	inquote = 0;
	slashes = 0;
	while (1) {
	    copy = 1;
	    while (*p == '\\') {
		slashes++;
		p++;
	    }
	    if (*p == '"') {
		if ((slashes & 1) == 0) {
		    copy = 0;
		    if ((inquote) && (p[1] == '"')) {
			p++;
			copy = 1;
		    } else {
			inquote = !inquote;
		    }
                }
                slashes >>= 1;
            }

            while (slashes) {
		*arg = '\\';
		arg++;
		slashes--;
	    }

	    if ((*p == '\0')
		    || (!inquote && ((*p == ' ') || (*p == '\t')))) { /* INTL: ISO space. */
		break;
	    }
	    if (copy != 0) {
		*arg = *p;
		arg++;
	    }
	    p++;
        }
	*arg = '\0';
	argSpace = arg + 1;
    }
    argv[argc] = NULL;

    *argcPtr = argc;
    *argvPtr = argv;
}

/*----------------------------------------------------------------------
 * set the toplevel window icon
 *----------------------------------------------------------------------*/

typedef struct TkWmInfo {
    TkWindow *winPtr;
    HWND hwnd;
} WmInfo;

static int Proc_LoadIcon(ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    TkWindow *winPtr;
    HANDLE hIcon;
    HWND hwnd = NULL;

    if (2 != argc) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
	    argv[0], " iconname\"", (char *) NULL);
        return TCL_ERROR;
    }

    winPtr = (TkWindow *) Tk_MainWindow (interp);
    if (NULL != winPtr && NULL != winPtr->wmInfoPtr)
        hwnd = winPtr->wmInfoPtr->hwnd;
    if (!hwnd) {
        Tcl_SetResult (interp,
            "toplevel window has not been created/mapped", TCL_STATIC);
        return TCL_ERROR;
    }
    if (argv[1][0] == '@')
        hIcon = LoadImage (NULL, &argv[1][1], IMAGE_ICON, 0, 0,
            LR_DEFAULTSIZE | LR_LOADFROMFILE);
    else
        hIcon = (HANDLE) LoadIcon (myInstance, argv[1]);
    if (!hIcon) {
        Tcl_SetResult (interp, "icon not found", TCL_STATIC);
        return TCL_ERROR;
    }
    SetClassLongPtr (hwnd, GCLP_HICON, (LONG_PTR)hIcon);
    Tcl_ResetResult (interp);
    return TCL_OK;
}

