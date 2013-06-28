#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "tcl.h"
#include "calc.h"
#include "cgnslib.h"

#ifndef CONST
# define CONST
#endif

static Tcl_Interp *global_interp;
static char message[1024] = "";

/*---------- get_error -----------------------------------------------
 * gets error message on parsing error
 *--------------------------------------------------------------------*/

static void get_error (int errnum, char *errmsg, int pos, char *str)
{
    char *p = message;

    if (errnum < 0) {
        sprintf (message, "error_exit {%s}", errmsg);
        Tcl_Eval (global_interp, message);
    }
    if (NULL != str) {
        sprintf (p, "%s\n", str);
        p += strlen (p);
        while (pos-- > 0)
            *p++ = '-';
        *p++ = '^';
    }
    strcpy (p, errmsg);
}

/*---------- get_names ----------------------------------------------
 * collect symbol names into Tcl interp result
 *-------------------------------------------------------------------*/

static int get_names (VECSYM *sym, void *data)
{
    Tcl_Interp *interp = (Tcl_Interp *)data;

    if (VECSYM_EQUSTR == vecsym_type(sym)) {
        if (vecsym_nargs(sym) > 0)
            sprintf (message, "%s(%d)", vecsym_name(sym),
                vecsym_nargs(sym));
        else
            strcpy (message, vecsym_name(sym));
    }
    else if (VECSYM_FUNC == vecsym_type(sym)) {
        if (vecsym_nargs(sym) < 0)
            sprintf (message, "%s(...)", vecsym_name(sym));
        else
            sprintf (message, "%s(%d)", vecsym_name(sym), vecsym_nargs(sym));
    }
    else if (VECSYM_VECTOR == vecsym_type(sym))
        sprintf (message, "%s[%ld]", vecsym_name(sym),
            (long)vecsym_veclen(sym));
    else
        strcpy (message, vecsym_name(sym));
    Tcl_AppendElement (interp, message);
    return 0;
}

/*---------- CalcReset ---------------------------------------------
 * reset calculator (symbol table)
 *------------------------------------------------------------------*/

static int CalcReset (ClientData data, Tcl_Interp *interp,
                      int argc, char **argv)
{
    cgnsCalcReset ();
    return TCL_OK;
}

/*---------- CalcInit ----------------------------------------------
 * load CGNS file and initialize
 *------------------------------------------------------------------*/

static int CalcInit (ClientData data, Tcl_Interp *interp,
                     int argc, char **argv)
{
    int nb, idum, modify = 0;
    char name[33];

    Tcl_ResetResult (interp);
    if (argc != 3) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " mode CGNSfile\"", NULL);
        return TCL_ERROR;
    }
    if (argv[1][0] == 'w' || argv[1][0] == 'm') modify = 1;
    cgnsCalcInit (argv[2], modify, get_error);
    for (nb = 1; nb <= NumBases; nb++) {
        if (cg_base_read (cgnsFile, nb, name, &idum, &idum))
            cgnsCalcFatal ((char *)cg_get_error());
        Tcl_AppendElement (interp, name);
    }
    return TCL_OK;
}

/*---------- CalcDone ----------------------------------------------
 * close CGNS file
 *------------------------------------------------------------------*/

static int CalcDone (ClientData data, Tcl_Interp *interp,
                     int argc, char **argv)
{
    cgnsCalcDone ();
    return TCL_OK;
}

/*---------- CalcBase ----------------------------------------------
 * set base for calculations
 *------------------------------------------------------------------*/

static int CalcBase (ClientData data, Tcl_Interp *interp,
                     int argc, char **argv)
{
    int nb, nz;
    cgsize_t sizes[9];
    char name[33];

    Tcl_ResetResult (interp);
    if (2 != argc) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " base\"", NULL);
        return TCL_ERROR;
    }
    nb = atoi (argv[1]);
    if (nb < 1 || nb > NumBases) {
        sprintf (message, "base number %d out of range", nb);
        Tcl_SetResult (interp, message, TCL_STATIC);
        return TCL_ERROR;
    }
    if (nb != cgnsBase) cgnsCalcBase (nb);
    for (nz = 1; nz <= NumZones; nz++) {
        if (cg_zone_read (cgnsFile, cgnsBase, nz, name, sizes))
            cgnsCalcFatal ((char *)cg_get_error());
        Tcl_AppendElement (interp, name);
    }
    return TCL_OK;
}

/*---------- CalcBaseInfo ------------------------------------------
 * get current base information
 *------------------------------------------------------------------*/

static int CalcBaseInfo (ClientData data, Tcl_Interp *interp,
                         int argc, char **argv)
{
    char *p = message;

    Tcl_ResetResult (interp);
    if (!cgnsBase) return TCL_OK;
    sprintf (p, "Name      : %s\n", BaseName);
    p += strlen (p);
    sprintf (p, "Phy Dim   : %d\nCell Dim  : %d\n", PhyDim, CellDim);
    p += strlen (p);
    sprintf (p, "Zones     : %d\nRef Values: %d\n", NumZones, NumReference);
    p += strlen (p);
#if CGNS_VERSION >= 2500
    sprintf (p, "Data Class: %s\n", cg_DataClassName(BaseClass));
    p += strlen (p);
    sprintf (p, "Units     : %s %s %s %s %s",
        cg_MassUnitsName(BaseUnits[0]),
        cg_LengthUnitsName(BaseUnits[1]),
        cg_TimeUnitsName(BaseUnits[2]),
        cg_TemperatureUnitsName(BaseUnits[3]),
        cg_AngleUnitsName(BaseUnits[4]));
#else
    sprintf (p, "Data Class: %s\n", DataClassName[BaseClass]);
    p += strlen (p);
    sprintf (p, "Units     : %s %s %s %s %s",
        MassUnitsName[BaseUnits[0]],
        LengthUnitsName[BaseUnits[1]],
        TimeUnitsName[BaseUnits[2]],
        TemperatureUnitsName[BaseUnits[3]],
        AngleUnitsName[BaseUnits[4]]);
#endif
    Tcl_SetResult (interp, message, TCL_STATIC);
    return TCL_OK;
}

/*---------- CalcZone ----------------------------------------------
 * set zone for calculations
 *------------------------------------------------------------------*/

static int CalcZone (ClientData data, Tcl_Interp *interp,
                     int argc, char **argv)
{
    int nz, ns;
    char name[33];
    CGNS_ENUMT(GridLocation_t) location;

    Tcl_ResetResult (interp);
    if (2 != argc) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " zone\"", NULL);
        return TCL_ERROR;
    }
    nz = atoi (argv[1]);
    if (nz < 1 || nz > NumZones) {
        sprintf (message, "zone number %d out of range", nz);
        Tcl_SetResult (interp, message, TCL_STATIC);
        return TCL_ERROR;
    }
    if (nz != cgnsZone) cgnsCalcZone (nz);
    for (ns = 1; ns <= NumSolns; ns++) {
        if (cg_sol_info (cgnsFile, cgnsBase, cgnsZone, ns, name, &location))
            cgnsCalcFatal ((char *)cg_get_error());
        Tcl_AppendElement (interp, name);
    }
    return TCL_OK;
}

/*---------- CalcZoneInfo ------------------------------------------
 * get current zone information
 *------------------------------------------------------------------*/

static int CalcZoneInfo (ClientData data, Tcl_Interp *interp,
                         int argc, char **argv)
{
    char *p = message;
    int n, dim;

    Tcl_ResetResult (interp);
    if (!cgnsZone) return TCL_OK;
    sprintf (p, "Name      : %s\nType      : %s\n", ZoneName,
#if CGNS_VERSION >= 2500
        cg_ZoneTypeName(ZoneType));
#else
        ZoneTypeName[ZoneType]);
#endif
    p += strlen (p);
    dim = ZoneType == CGNS_ENUMV(Structured) ? CellDim : 1;
    sprintf (p, "Dimensions: %d", ZoneDims[0]);
    p += strlen(p);
    for (n = 1; n < dim; n++) {
        sprintf (p, " x %d", ZoneDims[n]);
        p += strlen(p);
    }
    sprintf (p, "\nCoords    : %d\nSolutions : %d\n",
        NumCoordinates, NumSolns);
    p += strlen (p);
#if CGNS_VERSION >= 2500
    sprintf (p, "Data Class: %s\n", cg_DataClassName(ZoneClass));
    p += strlen (p);
    sprintf (p, "Units     : %s %s %s %s %s",
        cg_MassUnitsName(ZoneUnits[0]),
        cg_LengthUnitsName(ZoneUnits[1]),
        cg_TimeUnitsName(ZoneUnits[2]),
        cg_TemperatureUnitsName(ZoneUnits[3]),
        cg_AngleUnitsName(ZoneUnits[4]));
#else
    sprintf (p, "Data Class: %s\n", DataClassName[ZoneClass]);
    p += strlen (p);
    sprintf (p, "Units     : %s %s %s %s %s",
        MassUnitsName[ZoneUnits[0]],
        LengthUnitsName[ZoneUnits[1]],
        TimeUnitsName[ZoneUnits[2]],
        TemperatureUnitsName[ZoneUnits[3]],
        AngleUnitsName[ZoneUnits[4]]);
#endif
    Tcl_SetResult (interp, message, TCL_STATIC);
    return TCL_OK;
}

/*---------- CalcSoln ----------------------------------------------
 * set solution for calculations
 *------------------------------------------------------------------*/

static int CalcSoln (ClientData data, Tcl_Interp *interp,
                     int argc, char **argv)
{
    int ns, nv;

    Tcl_ResetResult (interp);
    if (2 != argc) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " soln\"", NULL);
        return TCL_ERROR;
    }
    ns = atoi (argv[1]);
    if (ns < 1 || ns > NumSolns) {
        sprintf (message, "solution number %d out of range", ns);
        Tcl_SetResult (interp, message, TCL_STATIC);
        return TCL_ERROR;
    }
    if (ns != cgnsSoln) cgnsCalcSoln (ns);
    for (nv = 0; nv < NumVariables; nv++)
        Tcl_AppendElement (interp, variables[nv].name);
    return TCL_OK;
}

/*---------- CalcSolnInfo ------------------------------------------
 * get current solution information
 *------------------------------------------------------------------*/

static int CalcSolnInfo (ClientData data, Tcl_Interp *interp,
                         int argc, char **argv)
{
    char *p = message;
    int n;

    Tcl_ResetResult (interp);
    if (!cgnsSoln) return TCL_OK;
    sprintf (p, "Name      : %s\nLocation  : %s\n", SolnName,
#if CGNS_VERSION >= 2500
        cg_GridLocationName(SolnLocation));
#else
        GridLocationName[SolnLocation]);
#endif
    p += strlen (p);
    if (ZoneType == CGNS_ENUMV(Structured)) {
        sprintf (p, "Dimensions: %d", SolnDims[0]);
        p += strlen(p);
        for (n = 1; n < CellDim; n++) {
            sprintf (p, " x %d", SolnDims[n]);
            p += strlen(p);
        }
        sprintf (p, "\nRind Cells: %d %d %d %d %d %d\n",
            SolnRind[0], SolnRind[1], SolnRind[2],
            SolnRind[3], SolnRind[4], SolnRind[5]);
    }
    else
        sprintf (p, "Dimensions: %d\n", SolnDims[0]);
    p += strlen (p);
    sprintf (p, "Fields    : %d\n", NumVariables);
    p += strlen (p);
#if CGNS_VERSION >= 2500
    sprintf (p, "Data Class: %s\n", cg_DataClassName(SolnClass));
    p += strlen (p);
    sprintf (p, "Units     : %s %s %s %s %s",
        cg_MassUnitsName(SolnUnits[0]),
        cg_LengthUnitsName(SolnUnits[1]),
        cg_TimeUnitsName(SolnUnits[2]),
        cg_TemperatureUnitsName(SolnUnits[3]),
        cg_AngleUnitsName(SolnUnits[4]));
#else
    sprintf (p, "Data Class: %s\n", DataClassName[SolnClass]);
    p += strlen (p);
    sprintf (p, "Units     : %s %s %s %s %s",
        MassUnitsName[SolnUnits[0]],
        LengthUnitsName[SolnUnits[1]],
        TimeUnitsName[SolnUnits[2]],
        TemperatureUnitsName[SolnUnits[3]],
        AngleUnitsName[SolnUnits[4]]);
#endif
    Tcl_SetResult (interp, message, TCL_STATIC);
    return TCL_OK;
}

/*---------- CalcRegList -------------------------------------------
 * return list of regions
 *------------------------------------------------------------------*/

static int CalcRegList (ClientData data, Tcl_Interp *interp,
                        int argc, char **argv)
{
    int n, dim;
    char *p = message;

    Tcl_ResetResult (interp);
    if (!NumZones || !cgnsZone) return TCL_OK;
    dim = ZoneType == CGNS_ENUMV(Structured) ? CellDim : 1;
    sprintf (p, "<Zone>[%d", ZoneDims[0]);
    p += strlen(p);
    for (n = 1; n < dim; n++) {
        sprintf (p, ",%d", ZoneDims[n]);
        p += strlen(p);
    }
    strcpy (p, "]");
    Tcl_AppendElement (interp, message);
    return TCL_OK;
}

/*---------- CalcRegInfo -------------------------------------------
 * return region information
 *------------------------------------------------------------------*/

static int CalcRegInfo (ClientData data, Tcl_Interp *interp,
                        int argc, char **argv)
{
    Tcl_ResetResult (interp);
    return TCL_OK;
}

/*---------- CalcVarList -------------------------------------------
 * return list of variables
 *------------------------------------------------------------------*/

static int CalcVarList (ClientData data, Tcl_Interp *interp,
                        int argc, char **argv)
{
    int n, i;

    Tcl_ResetResult (interp);
#if 0
    for (n = 0; n < NumCoordinates; n++) {
        for (i = 0; i < NumVariables; i++) {
            if (0 == strcmp (coordinates[n].name, variables[i].name))
                break;
        }
        if (i < NumVariables)
            sprintf (message, "'%s[%d]", coordinates[n].name,
                coordinates[n].len);
        else
            sprintf (message, "%s[%d]", coordinates[n].name,
                coordinates[n].len);
        Tcl_AppendElement (interp, message);
    }
    for (n = 0; n < NumVariables; n++) {
        sprintf (message, "%s[%d]", variables[n].name,
            variables[n].len);
        Tcl_AppendElement (interp, message);
    }
    for (n = 0; n < NumReference; n++) {
        if (reference[n].len > 1)
            sprintf (message, "~%s[%d]", reference[n].name,
                reference[n].len);
        else
            sprintf (message, "~%s", reference[n].name);
        Tcl_AppendElement (interp, message);
    }
#else
    for (n = 0; n < NumCoordinates; n++) {
        for (i = 0; i < NumVariables; i++) {
            if (0 == strcmp (coordinates[n].name, variables[i].name))
                break;
        }
        if (i < NumVariables) {
            sprintf (message, "'%s", coordinates[n].name);
            Tcl_AppendElement (interp, message);
        }
        else
            Tcl_AppendElement (interp, coordinates[n].name);
    }
    for (n = 0; n < NumVariables; n++)
        Tcl_AppendElement (interp, variables[n].name);
    for (n = 0; n < NumReference; n++) {
        sprintf (message, "~%s", reference[n].name);
        Tcl_AppendElement (interp, message);
    }
#endif
    return TCL_OK;
}

/*---------- CalcVarInfo -------------------------------------------
 * return variable information
 *------------------------------------------------------------------*/

static int CalcVarInfo (ClientData data, Tcl_Interp *interp,
                        int argc, char **argv)
{
    char *p = message;
    Variable *var;

    Tcl_ResetResult (interp);
    if (2 != argc) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " varname\"", NULL);
        return TCL_ERROR;
    }
    var = cgnsCalcVarGet (argv[1]);
    if (var == NULL) {
        Tcl_AppendResult (interp, "variable \"", argv[1],
            "\" not found", NULL);
        return TCL_ERROR;
    }
    sprintf (p, "Name      : %s\n", var->name);
    p += strlen (p);
    if (var->type == 0) {
          sprintf (p, "Type      : Reference\n");
          p += strlen(p);
          if (var->len > 1)
              sprintf (p, "Size      : %d\n", var->len);
          else
              sprintf (p, "Value     : %g\n", var->vd->f.val);
    }
    else if (var->type == 1 || SolnLocation == CGNS_ENUMV(Vertex)) {
          sprintf (p, "Type      : %s\nLocation  : Vertex\n",
              var->type == 1 ? "Coordinate" : "Solution");
          p += strlen (p);
          if (ZoneType == CGNS_ENUMV(Structured))
              sprintf (p, "Size      : %d (%d x %d x %d)\n", var->len,
                  ZoneDims[0], ZoneDims[1], ZoneDims[2]);
          else
              sprintf (p, "Size      : %d\n", var->len);
    }
    else  {
          sprintf (p, "Type      : Solution\nLocation  : %s\n",
#if CGNS_VERSION >= 2500
              cg_GridLocationName(SolnLocation));
#else
              GridLocationName[SolnLocation]);
#endif
          p += strlen (p);
          if (ZoneType == CGNS_ENUMV(Structured)) {
              sprintf (p, "Size      : %d (%d x %d x %d)\n", var->len,
                  SolnDims[0], SolnDims[1], SolnDims[2]);
              p += strlen (p);
              sprintf (p, "Rind Cells: %d %d %d %d %d %d\n",
                  SolnRind[0], SolnRind[1], SolnRind[2],
                  SolnRind[3], SolnRind[4], SolnRind[5]);
          }
          else
              sprintf (p, "Size      : %d\n", var->len);
    }
    p += strlen (p);
    sprintf (p, "Data Class: %s\nData Type : %s\n",
#if CGNS_VERSION >= 2500
        cg_DataClassName(var->dataclass),
        cg_DataTypeName(var->datatype));
#else
        DataClassName[var->dataclass], DataTypeName[var->datatype]);
#endif
    p += strlen (p);
    if (var->hasunits)
        sprintf (p, "Units     : %s %s %s %s %s\n",
#if CGNS_VERSION >= 2500
            cg_MassUnitsName(var->units[0]),
            cg_LengthUnitsName(var->units[1]),
            cg_TimeUnitsName(var->units[2]),
            cg_TemperatureUnitsName(var->units[3]),
            cg_AngleUnitsName(var->units[4]));
#else
            MassUnitsName[var->units[0]],
            LengthUnitsName[var->units[1]],
            TimeUnitsName[var->units[2]],
            TemperatureUnitsName[var->units[3]],
            AngleUnitsName[var->units[4]]);
#endif
    else
        sprintf (p, "Units     : not specified\n");
    p += strlen (p);
    if (var->hasexp)
        sprintf (p, "Exponents : %g %g %g %g %g\n",
            var->exponent[0], var->exponent[1], var->exponent[2],
            var->exponent[3], var->exponent[4]);
    else
        sprintf (p, "Exponents : not specified\n");
    p += strlen (p);
    if (var->hasconv)
        sprintf (p, "Conversion: %g %g", var->dataconv[0], var->dataconv[1]);
    else
        sprintf (p, "Conversion: not specified");
    Tcl_SetResult (interp, message, TCL_STATIC);
    return TCL_OK;
}

/*---------- CalcSymList -------------------------------------------
 * return list of symbols
 *------------------------------------------------------------------*/

static int CalcSymList (ClientData data, Tcl_Interp *interp,
                        int argc, char **argv)
{
    Tcl_ResetResult (interp);
    sym_list (0, get_names, interp);
    return TCL_OK;
}

/*---------- CalcSymInfo -------------------------------------------
 * return symbol information
 *------------------------------------------------------------------*/

static int CalcSymInfo (ClientData data, Tcl_Interp *interp,
                        int argc, char **argv)
{
    char *p = message;
    VECSYM *sym;
    Units *units;

    Tcl_ResetResult (interp);
    if (2 != argc) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " symname\"", NULL);
        return TCL_ERROR;
    }
    sym = find_symbol (argv[1], 0);
    if (sym == NULL) {
        sym = find_symbol (argv[1], 1);
        if (sym == NULL) {
            Tcl_AppendResult (interp, "symbol \"", argv[1],
                "\" not found", NULL);
            return TCL_ERROR;
        }
    }
    sprintf (p, "Name      : %s\n", vecsym_name(sym));
    p += strlen (p);
    if (vecsym_type(sym) == VECSYM_EQUSTR) {
        sprintf (p, "Type      : equation\n");
        p += strlen(p);
        if (vecsym_nargs(sym) > 0)
            sprintf (p, "Arguments : %d", vecsym_nargs(sym));
        else
            sprintf (p, "Arguments : none");
    }
    else if (vecsym_type(sym) == VECSYM_FUNC) {
        sprintf (p, "Type      : function\n");
        p += strlen(p);
        if (vecsym_nargs(sym) < 0)
            sprintf (p, "Arguments : variable");
        else
            sprintf (p, "Arguments : %d", vecsym_nargs(sym));
    }
    else if (vecsym_type(sym) == VECSYM_VECTOR)
        sprintf (p, "Type      : vector\nSize      : %ld",
            (long)vecsym_veclen(sym));
    else
        sprintf (p, "Type      : value\nValue     : %g",
            vecsym_value(sym));
    if (vecsym_user(sym) != NULL) {
        units = (Units *)vecsym_user(sym);
        p += strlen(p);
        sprintf (p, "\nUnits     : %s %s %s %s %s\n",
#if CGNS_VERSION >= 2500
            cg_MassUnitsName(units->units[0]),
            cg_LengthUnitsName(units->units[1]),
            cg_TimeUnitsName(units->units[2]),
            cg_TemperatureUnitsName(units->units[3]),
            cg_AngleUnitsName(units->units[4]));
#else
            MassUnitsName[units->units[0]],
            LengthUnitsName[units->units[1]],
            TimeUnitsName[units->units[2]],
            TemperatureUnitsName[units->units[3]],
            AngleUnitsName[units->units[4]]);
#endif
        p += strlen(p);
        sprintf (p, "Exponents : %g %g %g %g %g",
            units->exps[0], units->exps[1], units->exps[2],
            units->exps[3], units->exps[4]);
    }
    Tcl_SetResult (interp, message, TCL_STATIC);
    return TCL_OK;
}

/*---------- CalcCommand -------------------------------------------
 * parse expression and return results
 *------------------------------------------------------------------*/

static int CalcCommand (ClientData data, Tcl_Interp *interp,
                        int argc, char **argv)
{
    size_t n;
    VECFLOAT *f, fmin, fmax;
    VECSYM *sym;

    Tcl_ResetResult (interp);
    if (2 != argc) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " command\"", NULL);
        return TCL_ERROR;
    }

    sym = cgnsCalcCommand (argv[1]);
    if (NULL == sym) {
        Tcl_SetResult (interp, message, TCL_STATIC);
        return TCL_ERROR;
    }
    if (VECSYM_EQUSTR == vecsym_type(sym)) {
        if (vecsym_nargs(sym) > 0)
            sprintf (message, "%s(%d)", vecsym_name(sym),
                vecsym_nargs(sym));
        else
            strcpy (message, vecsym_name(sym));
    }
    else if (VECSYM_FUNC == vecsym_type(sym)) {
        if (vecsym_nargs(sym) < 0)
            sprintf (message, "%s(...)", vecsym_name(sym));
        else
            sprintf (message, "%s(%d)", vecsym_name(sym), vecsym_nargs(sym));
    }
    else if (VECSYM_VALUE == vecsym_type(sym))
        sprintf (message, "%s = %g", vecsym_name(sym), vecsym_value(sym));
    else if (VECSYM_VECTOR == vecsym_type(sym)) {
        f = vecsym_vector(sym);
        fmin = fmax = f[0];
        for (n = 1; n < vecsym_veclen(sym); n++) {
            if (fmin > f[n]) fmin = f[n];
            if (fmax < f[n]) fmax = f[n];
        }
        sprintf (message, "%s -> len = %ld, min = %g, max = %g",
            vecsym_name(sym), (long)vecsym_veclen(sym), fmin, fmax);
    }
    else
        strcpy (message, vecsym_name(sym));

    Tcl_SetResult (interp, message, TCL_STATIC);
    return TCL_OK;
}

/*---------- CalcDelete --------------------------------------------
 * delete symbols
 *------------------------------------------------------------------*/

static int CalcDelete (ClientData data, Tcl_Interp *interp,
                       int argc, char **argv)
{
    int i, n, nargs;
    CONST char **args;

    for (n = 1; n < argc; n++) {
        if (TCL_OK == Tcl_SplitList (interp, argv[n], &nargs, &args)) {
            for (i = 0; i < nargs; i++)
                sym_delsym ((char *)args[i]);
            Tcl_Free ((char *)args);
        }
        else
            sym_delsym (argv[n]);
    }
    return TCL_OK;
}

/*---------- CalcSave ----------------------------------------------
 * save computed variables to results file
 *------------------------------------------------------------------*/

static int CalcSave (ClientData data, Tcl_Interp *interp,
                     int argc, char **argv)
{
    Tcl_SetResult (interp, "not implemented yet", TCL_STATIC);
    return TCL_ERROR;
}

/*---------- ADFtcl_Init ---------------------------------------
 * Initialize and create the commands
 *--------------------------------------------------------------*/

#if defined(_WIN32) && defined(BUILD_DLL)
__declspec(dllexport)
#endif
int Calctcl_Init(Tcl_Interp *interp)
{
    global_interp = interp;
    Tcl_CreateCommand (interp, "CalcReset", (Tcl_CmdProc *)CalcReset,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcInit", (Tcl_CmdProc *)CalcInit,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcDone", (Tcl_CmdProc *)CalcDone,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcBase", (Tcl_CmdProc *)CalcBase,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcBaseInfo", (Tcl_CmdProc *)CalcBaseInfo,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcZone", (Tcl_CmdProc *)CalcZone,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcZoneInfo", (Tcl_CmdProc *)CalcZoneInfo,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcSoln", (Tcl_CmdProc *)CalcSoln,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcSolnInfo", (Tcl_CmdProc *)CalcSolnInfo,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcRegList", (Tcl_CmdProc *)CalcRegList,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcRegInfo", (Tcl_CmdProc *)CalcRegInfo,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcVarList", (Tcl_CmdProc *)CalcVarList,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcVarInfo", (Tcl_CmdProc *)CalcVarInfo,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcSymList", (Tcl_CmdProc *)CalcSymList,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcSymInfo", (Tcl_CmdProc *)CalcSymInfo,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcCommand", (Tcl_CmdProc *)CalcCommand,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcDelete", (Tcl_CmdProc *)CalcDelete,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CalcSave", (Tcl_CmdProc *)CalcSave,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    return TCL_OK;
}

