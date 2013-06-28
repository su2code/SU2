#include <stdio.h>
#include <string.h>
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef WIN32_LEAN_AND_MEAN
#include <htmlhelp.h>
#include <io.h>
#include <tcl.h>

/*---------- WinHtml -----------------------------------------------
 * run windows HTML help
 *------------------------------------------------------------------*/

static int WinHtml (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int n;
    char *p;
    static char winname[65] = "";
    static char hlpfile[256] = "";
    static char usg_msg[] =
        "WinHtml file CHMfile ?window?\n"
        "WinHtml topic HTMLfile ?tag?\n"
        "WinHtml index\n"
        "WinHtml close";

    if (argc < 2) {
        Tcl_SetResult (interp, usg_msg, TCL_STATIC);
        return TCL_ERROR;
    }
    Tcl_ResetResult (interp);

    if (0 == strcmp ("close", argv[1])) {
        HtmlHelp (NULL, NULL, HH_CLOSE_ALL, 0);
        return TCL_OK;
    }

    if (0 == strcmp ("file", argv[1])) {
        if (argc < 3 || argc > 4) {
            Tcl_SetResult (interp,
                "usage: WinHtml file CHMfile [window]", TCL_STATIC);
            return TCL_ERROR;
        }
        p = strrchr (argv[2], '.');
        if (NULL == p || stricmp (p, ".chm") || access (argv[2], 0)) {
            *hlpfile = 0;
            Tcl_AppendResult (interp, "CHM file \"", argv[2],
                "\" does not exist", NULL);
            return TCL_ERROR;
        }
        if (strlen (argv[2]) >= sizeof(hlpfile)) {
            Tcl_SetResult (interp, "CHM file pathname too long", TCL_STATIC);
            return TCL_ERROR;
        }
        strcpy (hlpfile, argv[2]);
        if (argc == 4)
            strcpy (winname, argv[3]);
        else
            *winname = 0;
        return TCL_OK;
    }

    n = strlen (hlpfile);
    if (!n) {
        Tcl_SetResult (interp, "CHM help file not specified", TCL_STATIC);
        return TCL_ERROR;
    }

    if (0 == strcmp ("index", argv[1])) {
        if (*winname)
            sprintf (&hlpfile[n], ">%s", winname);
        HtmlHelp (NULL, hlpfile, HH_DISPLAY_TOPIC, 0);
        hlpfile[n] = 0;
        return TCL_OK;
    }

    if (strcmp ("topic", argv[1])) {
        Tcl_SetResult (interp, usg_msg, TCL_STATIC);
        return TCL_ERROR;
    }

    if (argc < 3 || argc > 4) {
        Tcl_SetResult (interp,
           "usage: WinHtml topic HTMLfile [tag]", TCL_STATIC);
        return TCL_ERROR;
    }
    sprintf (&hlpfile[n], "::/%s", argv[2]);
    p = hlpfile + strlen (hlpfile);
    if (argc == 4) {
        sprintf (p, "#%s", argv[3]);
        p = hlpfile + strlen (hlpfile);
    }
    if (*winname)
        sprintf (p, ">%s", winname);
    HtmlHelp (NULL, hlpfile, HH_DISPLAY_TOPIC, 0);
    hlpfile[n] = 0;
    return TCL_OK;
}

/*----------------------------------------------------------------------*/

#if defined(_WIN32) && defined(BUILD_DLL)
__declspec(dllexport)
#endif
int WinHtml_Init(Tcl_Interp *interp)
{
    Tcl_CreateCommand (interp, "WinHtml", (Tcl_CmdProc *)WinHtml,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    return TCL_OK;
}

