#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "cgnslib.h"
#include "calc.h"
#include "vecerr.h"

#if CGNS_VERSION < 3000
# define Celsius Celcius
#endif

#ifndef CG_MODE_READ
# define CG_MODE_READ   MODE_READ
# define CG_MODE_MODIFY MODE_MODIFY
#endif

#ifdef READ_NODE
#include "cgns_header.h"
#include "cgns_io.h"
#endif

#ifndef CGNSTYPES_H
# define cgsize_t  int
# define cglong_t  long
# define cgulong_t unsigned long
#endif
#ifndef CGNS_ENUMT
# define CGNS_ENUMT(e) e
# define CGNS_ENUMV(e) e
#endif
#ifndef CG_MAX_INT32
# define CG_MAX_INT32 0x7FFFFFFF
#endif

int cgnsFile = 0;

/*--- base data ---*/

int NumBases = 0;
int cgnsBase;
char BaseName[33];
int CellDim, PhyDim;
int BaseClass, BaseUnits[5];

/*--- zone data ---*/

int NumZones = 0;
int cgnsZone;
char ZoneName[33];
int ZoneType;
int ZoneDims[6];
int ZoneClass, ZoneUnits[5];
int GridClass, GridUnits[5];

/*--- solution data ---*/

int NumSolns = 0;
int cgnsSoln;
char SolnName[33];
int SolnLocation;
int SolnDims[3], SolnRind[6];
int SolnClass, SolnUnits[5];

/*--- field data ---*/

int NumReference = 0;
Variable *reference;

int NumCoordinates = 0;
Variable *coordinates;

int NumVariables = 0;
Variable *variables;

/*--- local variables ---*/

static int cmdstrlen = 0;
static char *cmdstr;
static int VectorLen = 100;

/* unit specifications */

typedef struct {
    char *name;
    int type;
    int value;
} UnitSpec;

static UnitSpec unitspec[] = {
    {"cel", 3, CGNS_ENUMV(Celsius}),
    {"cen", 1, CGNS_ENUMV(Centimeter}),
    {"cm",  1, CGNS_ENUMV(Centimeter}),
    {"c",   3, CGNS_ENUMV(Celsius}),
    {"d",   4, CGNS_ENUMV(Degree}),
    {"fa",  3, CGNS_ENUMV(Fahrenheit}),
    {"fo",  1, CGNS_ENUMV(Foot}),
    {"ft",  1, CGNS_ENUMV(Foot}),
    {"f",   3, CGNS_ENUMV(Fahrenheit}),
    {"g",   0, CGNS_ENUMV(Gram}),
    {"in",  1, CGNS_ENUMV(Inch}),
    {"ke",  3, CGNS_ENUMV(Kelvin}),
    {"ki",  0, CGNS_ENUMV(Kilogram}),
    {"kg",  0, CGNS_ENUMV(Kilogram}),
    {"k",   3, CGNS_ENUMV(Kelvin}),
    {"lb",  0, CGNS_ENUMV(PoundMass}),
    {"me",  1, CGNS_ENUMV(Meter}),
    {"mi",  1, CGNS_ENUMV(Millimeter}),
    {"mm",  1, CGNS_ENUMV(Millimeter}),
    {"m",   1, CGNS_ENUMV(Meter}),
    {"p",   0, CGNS_ENUMV(PoundMass}),
    {"rad", 4, CGNS_ENUMV(Radian}),
    {"ran", 3, CGNS_ENUMV(Rankine}),
    {"r",   3, CGNS_ENUMV(Rankine}),
    {"se",  2, CGNS_ENUMV(Second}),
    {"sl",  0, CGNS_ENUMV(Slug}),
    {"s",   2, CGNS_ENUMV(Second})
};

#define NUM_UNITSPEC (sizeof(unitspec)/sizeof(UnitSpec))

/*---------- read_node ---------------------------------------------
 * read a node from the CGNS file
 *------------------------------------------------------------------*/

#ifdef READ_NODE

static VECDATA *read_node (char *nodename)
{
    cgns_file *cgfile;
    int n, bytes, dt, cgio;
    int ndim;
    cgsize_t np, dims[CGIO_MAX_DIMENSIONS];
    char *values;
    char type[CGIO_MAX_DATATYPE_LENGTH+1];
    char errmsg[CGIO_MAX_ERROR_LENGTH+1];
    double rootid, nodeid;
    VECDATA *vd;

    static struct dataTypes {
        char *name;
        int bytes;
    } data_types[6] = {
        {"I4", 4},
        {"I8", 8},
        {"U4", 4},
        {"U8", 8},
        {"R4", 4},
        {"R8", 8}
    };

    /* get node ID for node */

    cgfile = cgi_get_file (cgnsFile);
    cgio = cgfile->cgio;
    rootid = cgfile->rootid;

    if (cgio_get_node_id (cgio, rootid, nodename, &nodeid)) {
        cgio_error_message (errmsg);
        cgnsCalcFatal (errmsg);
    }

    /* get the type of data */

    if (cgio_get_data_type (cgio, nodeid, type)) {
        cgio_error_message (errmsg);
        cgnsCalcFatal (errmsg);
    }
    for (n = 0; n < CGIO_MAX_DATATYPE_LENGTH && type[n]; n++) {
        if (islower (type[n]))
            type[n] = toupper (type[n]);
    }
    for (bytes = 0, dt = 0; dt < 6; dt++) {
        if (0 == strncmp (type, data_types[dt].name, 2)) {
            bytes = data_types[dt].bytes;
            break;
        }
    }
    if (bytes == 0) {
        sprintf (errmsg, "can't handle data type %s", type);
        cgnsCalcFatal (errmsg);
    }

    /* get data dimensions */

    if (cgio_get_dimensions (cgio, nodeid, &ndim, dims)) {
        cgio_error_message (errmsg);
        cgnsCalcFatal (errmsg);
    }
    np = 0;
    if (ndim > 0) {
        for (np = 1, n = 0; n < ndim; n++)
            np *= dims[n];
    }
    if (np == 0)
        cgnsCalcFatal ("no data for node");
    if (np > CG_MAX_INT32)
        cgnsCalcFatal ("exceeded 32-bit integer");

    /* read the data */

    values = (char *) malloc ((size_t)(np * bytes));
    if (NULL == values)
        cgnsCalcFatal ("malloc failed for node data");

    if (cgio_read_all_data (cgio, nodeid, values)) {
        cgio_error_message (errmsg);
        cgnsCalcFatal (errmsg);
    }

    if (np == 1) {
        vd = vec_create (VEC_VALUE, 0, 1);
        if (dt == 0) {
            int *data = (int *)values;
            vd->f.val = (VECFLOAT)*data;
        }
        else if (dt == 1) {
            cglong_t *data = (cglong_t *)values;
            vd->f.val = (VECFLOAT)*data;
        }
        else if (dt == 2) {
            unsigned int *data = (unsigned int *)values;
            vd->f.val = (VECFLOAT)*data;
        }
        else if (dt == 3) {
            cgulong_t *data = (cgulong_t *)values;
            vd->f.val = (VECFLOAT)*data;
        }
        else if (dt == 4) {
            float *data = (float *)values;
            vd->f.val = (VECFLOAT)*data;
        }
        else {
            double *data = (double *)values;
            vd->f.val = (VECFLOAT)*data;
        }
    }
    else {
        vd = vec_create (VEC_VECTOR, (int)np, 1);
        if (dt == 0) {
            int *data = (int *)values;
            for (n = 0; n < np; n++)
                vd->f.vec[n] = (VECFLOAT)data[n];
        }
        else if (dt == 1) {
            cglong_t *data = (cglong_t *)values;
            for (n = 0; n < np; n++)
                vd->f.vec[n] = (VECFLOAT)data[n];
        }
        else if (dt == 2) {
            unsigned int *data = (unsigned int *)values;
            for (n = 0; n < np; n++)
                vd->f.vec[n] = (VECFLOAT)data[n];
        }
        else if (dt == 3) {
            cgulong_t *data = (cgulong_t *)values;
            for (n = 0; n < np; n++)
                vd->f.vec[n] = (VECFLOAT)data[n];
        }
        else if (dt == 4) {
            float *data = (float *)values;
            for (n = 0; n < np; n++)
                vd->f.vec[n] = (VECFLOAT)data[n];
        }
        else {
            double *data = (double *)values;
            for (n = 0; n < np; n++)
                vd->f.vec[n] = (VECFLOAT)data[n];
        }
    }

    free (values);
    return vd;
}

#endif

/*---------- read_units -----------------------------------------------
 * read unit specifications
 *---------------------------------------------------------------------*/

static int read_units (int units[5])
{
    int n;
    CGNS_ENUMT(MassUnits_t) mass;
    CGNS_ENUMT(LengthUnits_t) length;
    CGNS_ENUMT(TimeUnits_t) time;
    CGNS_ENUMT(TemperatureUnits_t) temp;
    CGNS_ENUMT(AngleUnits_t) angle;

    if (cg_units_read (&mass, &length, &time, &temp, &angle)) {
        for (n = 0; n < 5; n++)
            units[n] = 0;
        return 0;
    }
    units[0] = mass;
    units[1] = length;
    units[2] = time;
    units[3] = temp;
    units[4] = angle;
    return 1;
}

/*---------- read_class -----------------------------------------------
 * get data class, units and conversion factors
 *---------------------------------------------------------------------*/

static void read_class (Variable *var, int dataclass, int units[5])
{
    int i;
    CGNS_ENUMT(DataType_t) datatype;

    if (cg_dataclass_read ((CGNS_ENUMT(DataClass_t) *)&var->dataclass))
        var->dataclass = dataclass;
    var->hasunits = read_units (var->units);
    if (!var->hasunits) {
        for (i = 0; i < 5; i++)
            var->units[i] = units[i];
    }
    if (cg_conversion_info (&datatype) ||
       (datatype != CGNS_ENUMV(RealSingle) && datatype != CGNS_ENUMV(RealDouble))) {
        var->hasconv = 0;
        var->dataconv[0] = 1.0;
        var->dataconv[1] = 0.0;
    }
    else {
        var->hasconv = datatype;
        if (datatype == CGNS_ENUMV(RealSingle)) {
            float conv[2];
            if (cg_conversion_read (conv))
                cgnsCalcFatal ((char *)cg_get_error());
            for (i = 0; i < 2; i++)
                var->dataconv[i] = conv[i];
        }
        else {
            if (cg_conversion_read (var->dataconv))
                cgnsCalcFatal ((char *)cg_get_error());
        }
    }
    if (cg_exponents_info (&datatype) ||
       (datatype != CGNS_ENUMV(RealSingle) && datatype != CGNS_ENUMV(RealDouble))) {
        var->hasexp = 0;
        for (i = 0; i < 5; i++)
            var->exponent[i] = 0.0;
    }
    else {
        var->hasexp = datatype;
        if (datatype == CGNS_ENUMV(RealSingle)) {
            float exp[5];
            if (cg_exponents_read (exp))
                cgnsCalcFatal ((char *)cg_get_error());
            for (i = 0; i < 5; i++)
                var->exponent[i] = exp[i];
        }
        else {
            if (cg_exponents_read (var->exponent))
                cgnsCalcFatal ((char *)cg_get_error());
        }
    }
}

/*---------- read_reference --------------------------------------------
 * read reference conditions
 *----------------------------------------------------------------------*/

static void read_reference (void)
{
    int n, narrays, na, dim;
    cgsize_t vec[12];
    CGNS_ENUMT(DataType_t) datatype;
    char name[33];

    NumReference = 0;
    if (cg_goto (cgnsFile, cgnsBase, "ReferenceState_t", 1, "end") ||
        cg_narrays (&narrays) || narrays < 1)
        return;
    for (na = 1; na <= narrays; na++) {
        if (cg_array_info (na, name, &datatype, &dim, vec))
            cgnsCalcFatal ((char *)cg_get_error());
        if (datatype != CGNS_ENUMV(Character) && dim >= 1 && vec[0] >= 1)
            NumReference++;
    }
    if (!NumReference) return;
    reference = (Variable *) malloc (NumReference * sizeof(Variable));
    if (NULL == reference)
        cgnsCalcFatal ("malloc failed for reference variables");

    for (n = 0, na = 1; na <= narrays; na++) {
        if (cg_array_info (na, name, &datatype, &dim, vec))
            cgnsCalcFatal ((char *)cg_get_error());
        if (datatype != CGNS_ENUMV(Character) && dim >= 1 && vec[0] >= 1) {
            dim *= vec[0];
            strcpy (reference[n].name, name);
            reference[n].type = 0;
            reference[n].id = na;
            reference[n].len = dim;
            reference[n].valid = 1;
            reference[n].datatype = datatype;
            if (dim == 1) {
                reference[n].vd = vec_create (VEC_VALUE, 0, 0);
                if (cg_array_read_as (na, CGNS_ENUMV(RealDouble), &reference[n].vd->f.val))
                    cgnsCalcFatal ((char *)cg_get_error());
            }
            else {
                reference[n].vd = vec_create (VEC_VECTOR, dim, 0);
                if (cg_array_read_as (na, CGNS_ENUMV(RealDouble), reference[n].vd->f.vec))
                    cgnsCalcFatal ((char *)cg_get_error());
            }
            if (cg_goto (cgnsFile, cgnsBase, "ReferenceState_t", 1,
                    "DataArray_t", na, "end"))
                cgnsCalcFatal ((char *)cg_get_error());
            read_class (&reference[n], BaseClass, BaseUnits);
            cg_goto (cgnsFile, cgnsBase, "ReferenceState_t", 1, "end");
            n++;
        }
    }
}

/*---------- get_variable ----------------------------------------------
 * return variable data - read if neccessary
 *----------------------------------------------------------------------*/

static VECDATA *get_variable (Variable *var)
{
    if (var->vd == NULL) {
        int n;
        cgsize_t min[3], max[3];
        for (n = 0; n < 3; n++) {
            min[n] = 1;
            max[n] = SolnDims[n];
        }
        var->vd = vec_create (VEC_VECTOR, var->len, 0);
        if (cg_field_read (cgnsFile, cgnsBase, cgnsZone, cgnsSoln,
            var->name, CGNS_ENUMV(RealDouble), min, max, var->vd->f.vec))
            cgnsCalcFatal ((char *)cg_get_error());
    }
    return var->vd;
}

/*---------- get_coordinate --------------------------------------------
 * return coordinate data - read if neccessary
 *----------------------------------------------------------------------*/

static VECDATA *get_coordinate (Variable *var)
{
    if (var->vd == NULL) {
        int n;
        cgsize_t min[3], max[3];
        for (n = 0; n < 3; n++) {
            min[n] = 1;
            max[n] = ZoneDims[n];
        }
        var->vd = vec_create (VEC_VECTOR, var->len, 0);
        if (cg_coord_read (cgnsFile, cgnsBase, cgnsZone,
            var->name, CGNS_ENUMV(RealDouble), min, max, var->vd->f.vec))
            cgnsCalcFatal ((char *)cg_get_error());
    }
    return var->vd;
}

/*---------- print_error ---------------------------------------------
 * print error message on parsing error
 *--------------------------------------------------------------------*/

static void print_error (int errnum, char *errmsg, int pos, char *str)
{
    printf (errnum < 0 ? "FATAL:" : "ERROR:");
    if (NULL != str) {
        printf ("%s\n      ", str);
        while (pos-- > 0)
            putchar ('-');
        putchar ('^');
    }
    printf ("%s\n", errmsg);
}

/*---------- get_name -----------------------------------------------
 * get symbol name from string
 *-------------------------------------------------------------------*/

static char *get_name (char **str)
{
    int n;
    char *p = *str;
    static char name[SYMNAME_MAXLEN+1];

    if (*p == '"') {
        for (++p, n = 0; n < SYMNAME_MAXLEN && *p; n++) {
            if (*p == '"') break;
            name[n] = *p++;
        }
        if (*p++ != '"') return NULL;
    }
    else if (!isalpha (*p) && *p != '_')
        return (NULL);
    else {
        for (n = 0; n < SYMNAME_MAXLEN && *p; n++) {
            if (!isalnum(*p) && *p != '_') break;
            name[n] = *p++;
        }
    }
    name[n] = 0;
    *str = p;
    return name;
}

/*---------- print_symbols ------------------------------------------
 * print symbol names
 *-------------------------------------------------------------------*/

static int print_symbols (VECSYM *sym, void *data)
{
    FILE *fp = (FILE *)data;

    if (VECSYM_EQUSTR == vecsym_type(sym))
        fprintf (fp, "%s{%d}\n", vecsym_name(sym), vecsym_nargs(sym));
    else if (VECSYM_MACRO  == vecsym_type(sym))
        fprintf (fp, "%s<%d>\n", vecsym_name(sym), vecsym_nargs(sym));
    else if (VECSYM_FUNC == vecsym_type(sym)) {
        if (vecsym_nargs(sym) < 0)
            fprintf (fp, "%s(...)\n", vecsym_name(sym));
        else
            fprintf (fp, "%s(%d)\n", vecsym_name(sym), vecsym_nargs(sym));
    }
    else if (VECSYM_VECTOR == vecsym_type(sym))
        fprintf (fp, "%s[%ld]\n", vecsym_name(sym), (long)vecsym_veclen(sym));
    else
        fprintf (fp, "%s\n", vecsym_name(sym));
    return 0;
}

/*---------- delete_units -------------------------------------------
 * called when symbol deleted
 *-------------------------------------------------------------------*/

static void delete_units (VECSYM *sym)
{
    if (vecsym_user(sym) != NULL)
        free (vecsym_user(sym));
}

/*---------- callback -----------------------------------------------
 * callback function for vector parser
 *-------------------------------------------------------------------*/

static VECDATA *callback (int check, char **pp, char **err)
{
    int n, type = 0;
    char *p = *pp, *name;

    /* check for reading node data */

#ifdef READ_NODE
    if (*p == '{') {
        char nodename[257];
        for (n = 0; n < 256; n++) {
            if (!*++p || *p == '}') break;
            nodename[n] = *p;
        }
        if (n == 256) {
            *err = "internal node name length exceeded";
            return NULL;
        }
        if (!n || *p != '}') {
            *err = "incomplete node name specification";
            return NULL;
        }
        *pp = ++p;
        return read_node (nodename);
    }
#endif

    /* check for reference or coordinate data */

    if (*p == '~' || *p == '\'')
        type = *p++;

    /* get name */

    if ((name = get_name (&p)) != NULL) {

        /* check for variable */

        if (!type) {
            for (n = 0; n < NumVariables; n++) {
                if (0 == strcmp (variables[n].name, name)) {
                    *pp = p;
                    return get_variable (&variables[n]);
                }
            }
        }

        /* check for grid coordinates */

        if (type != '~') {
            for (n = 0; n < NumCoordinates; n++) {
                if (0 == strcmp (coordinates[n].name, name)) {
                    *pp = p;
                    return get_coordinate (&coordinates[n]);
                }
            }
        }

        /* check for reference quantity */

        if (type != '\'') {
            for (n = 0; n < NumReference; n++) {
                if (0 == strcmp (reference[n].name, name)) {
                    *pp = p;
                    return reference[n].vd;
                }
            }
        }
    }

    return (NULL);
}

/*---------- parse_units -------------------------------------------
 * get unit specification
 *------------------------------------------------------------------*/

static Units *parse_units (char **pp)
{
    int n, par, div;
    char *p = *pp, name[33];
    float exp;
    Units units, *u;
    UnitSpec *us;

    for (n = 0; n < 5; n++) {
        units.units[n] = 0;
        units.exps[n] = 0.0;
    }
    par = div = 0;
    while (1) {
        while (*p && isspace (*p))
            p++;
        if (*p == '*' || *p == '-') {
            p++;
            continue;
        }
        if (*p == '/') {
            div++;
            p++;
            continue;
        }
        if (*p == '(') {
            par++;
            p++;
            continue;
        }
        if (*p == ')') {
            if (!par--) return NULL;
            if (div) div--;
            p++;
            continue;
        }
        n = 0;
        while (*p && isalpha (*p)) {
            if (n < 32) name[n++] = tolower (*p);
            p++;
        }
        if (!n) break;
        name[n] = 0;
        for (n = 0; n < NUM_UNITSPEC; n++) {
            if (name[0] == unitspec[n].name[0]) break;
        }
        us = NULL;
        while (n < NUM_UNITSPEC) {
            if (name[0] != unitspec[n].name[0]) break;
            if (!strncmp (name, unitspec[n].name, strlen(unitspec[n].name))) {
                us = &unitspec[n];
                break;
            }
            n++;
        }
        if (us == NULL) return NULL;
        while (*p && isspace (*p))
            p++;
        if (*p == '^') {
            if (1 != sscanf (++p, "%f%n", &exp, &n)) return NULL;
            for (p += n; *p && isspace(*p); p++)
                ;
        }
        else
            exp = 1.0;
        units.units[us->type] = us->value;
        if (div) {
            units.exps[us->type] -= exp;
            if (*p != '-' && div > par) div--;
        }
        else
            units.exps[us->type] += exp;
    }
    u = (Units *) malloc (sizeof(Units));
    if (u == NULL)
        cgnsCalcFatal ("malloc failed for unit specification");
    for (n = 0; n < 5; n++) {
        u->units[n] = units.units[n];
        u->exps[n] = units.exps[n];
    }
    *pp = p;
    return u;
}

/*---------- free_all ----------------------------------------------
 * free all data
 *------------------------------------------------------------------*/

static void free_all (void)
{
    int n;

    if (NumReference) {
        for (n = 0; n < NumReference; n++)
            vec_destroy (reference[n].vd);
        free (reference);
        NumReference = 0;
    }
    if (NumCoordinates) {
        for (n = 0; n < NumCoordinates; n++)
            vec_destroy (coordinates[n].vd);
        free (coordinates);
        NumCoordinates = 0;
    }
    if (NumVariables) {
        for (n = 0; n < NumVariables; n++)
            vec_destroy (variables[n].vd);
        free (variables);
        NumVariables = 0;
    }
}

/*---------- cgnsCalcFatal -----------------------------------------
 * terminate with error message
 *------------------------------------------------------------------*/

void cgnsCalcFatal (char *errmsg)
{
    cgnsCalcDone ();
    if (NULL != errmsg && *errmsg) {
        if (NULL == vec_errhandler)
            print_error (-1, errmsg, 0, NULL);
        else
            (*vec_errhandler) (-1, errmsg, 0, NULL);
    }
    exit (-1);
}

/*---------- cgnsCalcError -----------------------------------------
 * print error message
 *------------------------------------------------------------------*/

void cgnsCalcError (char *errmsg)
{
    if (NULL != errmsg && *errmsg) {
        if (NULL == vec_errhandler)
            print_error (0, errmsg, 0, NULL);
        else
            (*vec_errhandler) (0, errmsg, 0, NULL);
    }
}

/*---------- cgnsCalcReset -----------------------------------------
 * reset calculator (symbol table)
 *------------------------------------------------------------------*/

void cgnsCalcReset (void)
{
    sym_free ();
#ifdef EXTERN_FUNCS
    add_funcs ();
#endif
    VectorLen = 100;
}

/*---------- cgnsCalcInit ------------------------------------------
 * load CGNS file and initialize
 *------------------------------------------------------------------*/

int cgnsCalcInit (char *cgnsfile, int modify,
                  void (*errhandler)(int,char *,int,char *))
{
    cgnsCalcDone ();

    /* set up error handler */

    if (NULL == errhandler)
        vec_errhandler = print_error;
    else
        vec_errhandler = errhandler;

    /* callback to delete unit data */

    sym_delfunc = delete_units;

    /* open CGNS file */

    if (modify) {
        if (cg_open (cgnsfile, CG_MODE_MODIFY, &cgnsFile))
            cgnsCalcFatal ("couldn't open file in modify mode");
    }
    else {
        if (cg_open (cgnsfile, CG_MODE_READ, &cgnsFile))
            cgnsCalcFatal ("couldn't open file in read mode");
    }
    if (cg_nbases (cgnsFile, &NumBases))
        cgnsCalcFatal ((char *)cg_get_error());
    if (NumBases < 1)
        cgnsCalcFatal ("no bases found in CGNS file");

    cgnsCalcBase (1);

    /* return number of bases */

    return NumBases;
}

/*---------- cgnsCalcDone ------------------------------------------
 * close CGNS file
 *------------------------------------------------------------------*/

void cgnsCalcDone (void)
{
    if (cgnsFile) {
        cg_close (cgnsFile);
        cgnsFile = 0;
    }
    free_all ();
    cgnsBase = NumBases = 0;
    cgnsZone = NumZones = 0;
    cgnsSoln = NumSolns = 0;
}

/*---------- cgnsCalcBase ------------------------------------------
 * set base for calculations
 *------------------------------------------------------------------*/

int cgnsCalcBase (int base)
{
    if (base < 1 || base > NumBases)
        cgnsCalcFatal ("invalid base specified");
    if (cg_base_read (cgnsFile, base, BaseName, &CellDim, &PhyDim))
        cgnsCalcFatal ((char *)cg_get_error());
    free_all ();
    cgnsBase = base;
    cgnsZone = NumZones = 0;
    cgnsSoln = NumSolns = 0;

    /* read base class and units */

    if (cg_goto (cgnsFile, cgnsBase, "end"))
        cgnsCalcFatal ((char *)cg_get_error());
    if (cg_dataclass_read ((CGNS_ENUMT(DataClass_t) *)&BaseClass))
        BaseClass = 0;
    read_units (BaseUnits);

    /* read reference conditions */

    read_reference ();

    /* get number of zones and initilize */

    if (cg_nzones (cgnsFile, cgnsBase, &NumZones))
        cgnsCalcFatal ((char *)cg_get_error());
    if (NumZones) cgnsCalcZone (1);

    return NumZones;
}

/*---------- cgnsCalcZone ------------------------------------------
 * set zone for calculations
 *------------------------------------------------------------------*/

int cgnsCalcZone (int zone)
{
    int n, nc, size = 1;
    cgsize_t dims[9];
    CGNS_ENUMT(DataType_t) datatype;
    char name[33];

    if (zone < 1 || zone > NumZones)
        cgnsCalcFatal ("invalid zone specified");
    if (cg_zone_read (cgnsFile, cgnsBase, zone, ZoneName, dims) ||
        cg_zone_type (cgnsFile, cgnsBase, zone, (CGNS_ENUMT(ZoneType_t) *)&ZoneType))
        cgnsCalcFatal ((char *)cg_get_error());
    cgnsZone = zone;
    cgnsSoln = NumSolns = 0;

    for (n = 0; n < 6; n++)
        ZoneDims[n] = 1;
    if (ZoneType == CGNS_ENUMV(Structured)) {
        for (n = 0; n < CellDim; n++) {
            ZoneDims[n] = dims[n];
            ZoneDims[n+3] = dims[n+CellDim];
            size *= dims[n];
        }
    }
    else if (ZoneType == CGNS_ENUMV(Unstructured)) {
        ZoneDims[0] = dims[0];
        ZoneDims[3] = dims[1];
        size = dims[0];
    }
    else
        cgnsCalcFatal ("invalid zone type");
    VectorLen = size;

    /* free-up previous data */

    if (NumCoordinates) {
        for (n = 0; n < NumCoordinates; n++)
            vec_destroy (coordinates[n].vd);
        free (coordinates);
        NumCoordinates = 0;
    }
    if (NumVariables) {
        for (n = 0; n < NumVariables; n++)
            vec_destroy (variables[n].vd);
        free (variables);
        NumVariables = 0;
    }

    /* read zone class and units */

    if (cg_goto (cgnsFile, cgnsBase, "Zone_t", cgnsZone, "end"))
        cgnsCalcFatal ((char *)cg_get_error());
    if (cg_dataclass_read ((CGNS_ENUMT(DataClass_t) *)&ZoneClass))
        ZoneClass = BaseClass;
    if (!read_units (ZoneUnits)) {
        for (n = 0; n < 5; n++)
            ZoneUnits[n] = BaseUnits[n];
    }
    GridClass = ZoneClass;
    for (n = 0; n < 5; n++)
        GridUnits[n] = ZoneUnits[n];

    /* get coordinate info */

    if (cg_ncoords (cgnsFile, cgnsBase, cgnsZone, &nc))
        cgnsCalcFatal ((char *)cg_get_error());
    if (nc > 0) {

        /* read grid class and units */

        if (cg_goto (cgnsFile, cgnsBase, "Zone_t", cgnsZone,
            "GridCoordinates_t", 1, "end"))
            cgnsCalcFatal ((char *)cg_get_error());
        if (cg_dataclass_read ((CGNS_ENUMT(DataClass_t) *)&GridClass))
            GridClass = ZoneClass;
        if (!read_units (GridUnits)) {
            for (n = 0; n < 5; n++)
                GridUnits[n] = ZoneUnits[n];
        }

        /* read coordinates */

        NumCoordinates = nc;
        coordinates = (Variable *) malloc (NumCoordinates * sizeof(Variable));
        if (NULL == coordinates)
            cgnsCalcFatal ("malloc failed for coordinate info");

        for (n = 0; n < NumCoordinates; n++) {
            if (cg_coord_info (cgnsFile, cgnsBase, cgnsZone, n+1,
                &datatype, name))
                cgnsCalcFatal ((char *)cg_get_error());
            strcpy (coordinates[n].name, name);
            coordinates[n].type = 1;
            coordinates[n].id = n + 1;
            coordinates[n].len = size;
            coordinates[n].valid = 1;
            coordinates[n].datatype = datatype;
            if (cg_goto (cgnsFile, cgnsBase, "Zone_t", cgnsZone,
                "GridCoordinates_t", 1, "DataArray_t", n+1, "end"))
                cgnsCalcFatal ((char *)cg_get_error());
            read_class (&coordinates[n], GridClass, GridUnits);
            coordinates[n].vd = NULL;
        }
    }

    /* get number of solutions and initilize */

    if (cg_nsols (cgnsFile, cgnsBase, cgnsZone, &NumSolns))
        cgnsCalcFatal ((char *)cg_get_error());
    if (NumSolns) cgnsCalcSoln (1);

    return NumSolns;
}

/*---------- cgnsCalcSoln ------------------------------------------
 * set solution for calculations
 *------------------------------------------------------------------*/

int cgnsCalcSoln (int soln)
{
    int i, n, size, nflds;
    CGNS_ENUMT(DataType_t) datatype;
    char name[33];

    if (soln < 1 || soln > NumSolns)
        cgnsCalcFatal ("invalid solution specified");
    if (cg_sol_info (cgnsFile, cgnsBase, cgnsZone, soln, SolnName,
        (CGNS_ENUMT(GridLocation_t) *)&SolnLocation))
        cgnsCalcFatal ((char *)cg_get_error());
    cgnsSoln = soln;

    for (n = 0; n < 3; n++)
        SolnDims[n] = 1;
    if (ZoneType == CGNS_ENUMV(Structured)) {
        size = 1;
        if (SolnLocation == CGNS_ENUMV(CellCenter)) {
            if (cg_goto (cgnsFile, cgnsBase, "Zone_t", cgnsZone,
                "FlowSolution_t", cgnsSoln, "end"))
                cgnsCalcFatal ((char *)cg_get_error());
            if (cg_rind_read (SolnRind)) {
                for (n = 0; n < 6; n++)
                    SolnRind[n] = 0;
            }
            for (i = 0, n = 0; n < CellDim; n++, i += 2) {
                SolnDims[n] = ZoneDims[n] - 1 + SolnRind[i] + SolnRind[i+1];
                size *= SolnDims[n];
            }
        }
        else {
            for (n = 0; n < 6; n++)
                SolnRind[n] = 0;
            for (n = 0; n < CellDim; n++) {
                SolnDims[n] = ZoneDims[n];
                size *= SolnDims[n];
            }
        }
    }
    else {
        for (n = 0; n < 6; n++)
            SolnRind[n] = 0;
        size = SolnLocation == CGNS_ENUMV(CellCenter) ? ZoneDims[3] : ZoneDims[0];
        SolnDims[0] = size;
    }
    VectorLen = size;

    /* free-up previous data */

    if (NumVariables) {
        for (n = 0; n < NumVariables; n++)
            vec_destroy (variables[n].vd);
        free (variables);
        NumVariables = 0;
    }

    /* read solution class and units */

    if (cg_goto (cgnsFile, cgnsBase, "Zone_t", cgnsZone,
        "FlowSolution_t", cgnsSoln, "end"))
        cgnsCalcFatal ((char *)cg_get_error());
    if (cg_dataclass_read ((CGNS_ENUMT(DataClass_t) *)&SolnClass))
        SolnClass = ZoneClass;
    if (!read_units (SolnUnits)) {
        for (n = 0; n < 5; n++)
            SolnUnits[n] = ZoneUnits[n];
    }

    /* get field info */

    if (cg_nfields (cgnsFile, cgnsBase, cgnsZone, cgnsSoln, &nflds))
        cgnsCalcFatal ((char *)cg_get_error());
    if (nflds > 0) {
        NumVariables = nflds;
        variables = (Variable *) malloc (NumVariables * sizeof(Variable));
        if (NULL == variables)
            cgnsCalcFatal ("malloc failed for field info");
        for (n = 0; n < NumVariables; n++) {
            if (cg_field_info (cgnsFile, cgnsBase, cgnsZone, cgnsSoln,
                    n+1, &datatype, name))
                cgnsCalcFatal ((char *)cg_get_error());
            strcpy (variables[n].name, name);
            variables[n].type = 2;
            variables[n].id = n + 1;
            variables[n].len = size;
            variables[n].valid = 1;
            variables[n].datatype = datatype;
            if (cg_goto (cgnsFile, cgnsBase, "Zone_t", cgnsZone,
                "FlowSolution_t", cgnsSoln, "DataArray_t", n+1, "end"))
                cgnsCalcFatal ((char *)cg_get_error());
            read_class (&variables[n], SolnClass, SolnUnits);
            variables[n].vd = NULL;
        }
    }

    return NumVariables;
}

/*---------- cgnsCalcCheck -----------------------------------------
 * check expression
 *------------------------------------------------------------------*/

int cgnsCalcCheck (char *expression)
{
    int length = (int)strlen (expression);
    char *p, *cmd, *name;
    Units *units;

    if (length > cmdstrlen) {
        if (cmdstrlen == 0) {
            cmdstrlen = length > 256 ? length : 256;
            cmdstr = (char *) malloc (cmdstrlen + 1);
        }
        else {
            cmdstrlen = length + 32;
            cmdstr = (char *) realloc (cmdstr, cmdstrlen + 1);
        }
        if (cmdstr == NULL)
            cgnsCalcFatal ("malloc failed for cmdstr");
    }
    cmd = strcpy (cmdstr, expression);

    for (p = cmd + strlen(cmd) - 1; p >= cmd && isspace(*p); p--)
        ;
    *++p = 0;
    while (*cmd && isspace(*cmd))
        cmd++;

    if (!*cmd) return (0);

    p = cmd;
    if ((name = get_name (&p)) != NULL) {
        int nargs = -1;
        char *equ;

        for (equ = p; *equ && isspace(*equ); equ++)
            ;

        /* check for units */

        if ('[' == *equ) {
            equ++;
            units = parse_units (&equ);
            if (units == NULL || *equ != ']') {
                cgnsCalcError ("bad units specification");
                if (units != NULL) free (units);
                return 0;
            }
            free (units);
            while (*++equ && isspace (*equ))
                ;
            if (!*equ) {
                if (find_symbol (name, 0) != NULL) return 1;
            }
        }

        /* check for equation */

        if ('(' == *equ) {
            char *arg = equ;
            while (*++arg && (isspace (*arg) || isdigit(*arg)))
                ;
            if (')' == *arg) {
                nargs = atoi (equ + 1);
                for (equ = arg+1; *equ && isspace (*equ); equ++)
                    ;
            }
        }

        if ('=' == *equ && '=' != *++equ) {
            for (cmd = equ; *cmd && isspace(*cmd); cmd++)
                ;
            if (nargs > 9) {
                cgnsCalcError ("invalid number of equation arguments");
                return 0;
            }
        }
    }

    return vec_check (cmd, VectorLen, callback);
}

/*---------- cgnsCalcCommand ---------------------------------------
 * parse expression and return results
 *------------------------------------------------------------------*/

VECSYM *cgnsCalcCommand (char *expression)
{
    int n, length = (int)strlen (expression);
    char *p, *cmd, *name, sym[SYMNAME_MAXLEN+1];
    VECDATA *vd;
    Units *units = NULL;

    if (length > cmdstrlen) {
        if (cmdstrlen == 0) {
            cmdstrlen = length > 256 ? length : 256;
            cmdstr = (char *) malloc (cmdstrlen + 1);
        }
        else {
            cmdstrlen = length + 32;
            cmdstr = (char *) realloc (cmdstr, cmdstrlen + 1);
        }
        if (cmdstr == NULL)
            cgnsCalcFatal ("malloc failed for cmdstr");
    }
    cmd = strcpy (cmdstr, expression);

    /* skip leading and trailing space */

    for (p = cmd + strlen(cmd) - 1; p >= cmd && isspace(*p); p--)
        ;
    *++p = 0;
    while (*cmd && isspace(*cmd))
        cmd++;

    /* empty string */

    if (!*cmd) return (NULL);

    /* check for defining a new symbol */

    p = cmd;
    strcpy (sym, "_temp_");

    if ((name = get_name (&p)) != NULL) {
        int nargs = -1;
        char *equ;

        for (equ = p; *equ && isspace(*equ); equ++)
            ;

        /* check for units */

        if ('[' == *equ) {
            equ++;
            units = parse_units (&equ);
            if (units == NULL || *equ != ']') {
                cgnsCalcError ("bad units specification");
                if (units != NULL) free (units);
                return (NULL);
            }
            while (*++equ && isspace (*equ))
                ;
            if (!*equ) {
                VECSYM *symbol = find_symbol (name, 0);
                if (symbol != NULL) {
                    if (vecsym_user(symbol) != NULL)
                        free (vecsym_user(symbol));
                    vecsym_user(symbol) = units;
                    return (symbol);
                }
            }
        }

        /* check for equation */

        if ('(' == *equ) {
            char *arg = equ;
            while (*++arg && (isspace (*arg) || isdigit(*arg)))
                ;
            if (')' == *arg) {
                nargs = atoi (equ + 1);
                for (equ = arg+1; *equ && isspace (*equ); equ++)
                    ;
            }
        }

        if ('=' == *equ && '=' != *++equ) {
            strcpy (sym, name);
            for (cmd = equ; *cmd && isspace(*cmd); cmd++)
                ;

            /* add equation as string */

            if (nargs >= 0) {
                if (nargs > 9) {
                    cgnsCalcError ("invalid number of equation arguments");
                    return (NULL);
                }
                n = sym_addequ (sym, nargs, cmd, units);
                if (n) {
                    cgnsCalcError (sym_errmsg (n));
                    return (NULL);
                }
                return (find_symbol (sym, 1));
            }
        }
    }

    vd = vec_parse (cmd, VectorLen, callback);

    if (NULL == vd) {
        if (NULL != units) free (units);
        return (NULL);
    }

    /* add to symbol table */

    if (VEC_VALUE == vd->type)
        n = sym_addval (sym, vd->f.val, units);
    else
        n = sym_addvec (sym, vd->len, vd->f.vec, units);
    vec_destroy (vd);
    if (!n)
        return (find_symbol (sym, 0));
    cgnsCalcError (sym_errmsg (n));
    return (NULL);
}

/*---------- cgnsCalcVarGet ----------------------------------------
 * return a variable
 *------------------------------------------------------------------*/

Variable *cgnsCalcVarGet (char *varname)
{
    int n, type = 0;
    char *name = varname;

    if (name == NULL || !*name) return NULL;
    if (*name == '~' || *name == '\'')
        type = *name++;

    /* check for variable */

    if (!type) {
        for (n = 0; n < NumVariables; n++) {
            if (0 == strcmp (variables[n].name, name))
                return &variables[n];
        }
    }

    /* check for grid coordinates */

    if (type != '~') {
        for (n = 0; n < NumCoordinates; n++) {
            if (0 == strcmp (coordinates[n].name, name))
                return &coordinates[n];
        }
    }

    /* check for reference quantity */

    if (type != '\'') {
        for (n = 0; n < NumReference; n++) {
            if (0 == strcmp (reference[n].name, name))
                return &reference[n];
        }
    }
    return NULL;
}

/*---------- cgnsCalcVarList ---------------------------------------
 * print variables
 *------------------------------------------------------------------*/

void cgnsCalcVarList (FILE *fp)
{
    int n;

    if (NULL == fp)
        fp = stdout;

    if (NumReference) {
        fprintf (fp, "=== Reference ===\n");
        for (n = 0; n < NumReference; n++) {
            if (reference[n].len > 1)
                fprintf (fp, "~%s[%d]\n", reference[n].name, reference[n].len);
            else
                fprintf (fp, "~%s\n", reference[n].name);
        }
    }

    if (NumCoordinates) {
        fprintf (fp, "=== Coordinates ===\n");
        for (n = 0; n < NumCoordinates; n++)
            fprintf (fp, "'%s[%d]\n", coordinates[n].name, coordinates[n].len);
    }

    if (NumVariables) {
        fprintf (fp, "=== Solution ===\n");
        for (n = 0; n < NumVariables; n++)
            fprintf (fp, "%s[%d]\n", variables[n].name, variables[n].len);
    }
}

/*---------- cgnsCalcSymList ---------------------------------------
 * print list of symbols
 *------------------------------------------------------------------*/

void cgnsCalcSymList (FILE *fp)
{
    if (NULL == fp)
        fp = stdout;
    fprintf (fp, "=== Symbols ===\n");
    sym_list (0, print_symbols, fp);
}

