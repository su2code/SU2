#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "tk.h"
#include "cgns_io.h"
#include "cgnslib.h" /* only needed for CGNS_VERSION */

#ifndef CGNSTYPES_H
# define cgsize_t  int
# define cglong_t  long
# define cgulong_t unsigned long
#endif

#ifndef CONST
# define CONST
#endif

/* these are the data types as used in CGIO */

#define I4 int
#define U4 unsigned int
#define I8 cglong_t
#define U8 cgulong_t
#define R4 float
#define R8 double
#define X4 float
#define X8 double
#define C1 char
#define B1 unsigned char

enum DataTags {
    MTdata = 0,
    I4data,
    I8data,
    U4data,
    U8data,
    R4data,
    R8data,
    X4data,
    X8data,
    C1data,
    B1data,
    LKdata
};

static struct DataType {
    char *name;
    int type;
    int bytes;
    int size;
} DataList[] = {
    {"MT", MTdata,  0, 0},
    {"I4", I4data,  4, sizeof(I4)},
    {"I8", I8data,  8, sizeof(I8)},
    {"U4", U4data,  4, sizeof(U4)},
    {"U8", U8data,  8, sizeof(U8)},
    {"R4", R4data,  4, sizeof(R4)},
    {"R8", R8data,  8, sizeof(R8)},
    {"X4", X4data,  8, sizeof(X4)*2},
    {"X8", X8data, 16, sizeof(X8)*2},
    {"C1", C1data,  1, 1},
    {"B1", B1data,  1, 1},
    {"LK", LKdata,  0, 0}
};

#define NumDataTypes (sizeof(DataList)/sizeof(struct DataType))

static int cgioNum = 0;
static double RootID = 0.0;

/*===================================================================
 *           local routines
 *===================================================================*/

/*---------- get_error ----------------------------------------------
 * save error message
 *-------------------------------------------------------------------*/

static int get_error (Tcl_Interp *interp, char *funcname)
{
    char errmsg[CGIO_MAX_ERROR_LENGTH+1];

    cgio_error_message (errmsg);
    Tcl_AppendResult (interp, /*funcname, ":",*/ errmsg, NULL);
    return TCL_ERROR;
}

/*---------- get_node ----------------------------------------------
 * get the node ID for a node
 *------------------------------------------------------------------*/

static int get_node (Tcl_Interp *interp, char *nodename, double *nodeid)
{
    if (!*nodename || 0 == strcmp (nodename, "/"))
        *nodeid = RootID;
    else {
        if (cgio_get_node_id (cgioNum, RootID, nodename, nodeid))
            return (get_error (interp, "cgio_get_node_id"));
    }
    return TCL_OK;
}

/*---------- get_parent ----------------------------------------------
 * get the parent node ID for a node
 *------------------------------------------------------------------*/

static char *get_parent (Tcl_Interp *interp, char *nodename, double *nodeid)
{
    char *node;

    if (!*nodename || 0 == strcmp (nodename, "/") ||
        NULL == (node = strrchr (nodename, '/'))) {
        Tcl_AppendResult (interp, "node does not have a parent", NULL);
        return NULL;
    }
    if (node == nodename)
        *nodeid = RootID;
    else {
        *node = 0;
        if (cgio_get_node_id (cgioNum, RootID, nodename, nodeid)) {
            get_error (interp, "cgio_get_node_id");
            return NULL;
        }
        *node = '/';
    }
    return (++node);
}

/*---------- data_size ---------------------------------------------
 * return number of data values for a node
 *------------------------------------------------------------------*/

static cgsize_t data_size (Tcl_Interp *interp, double nodeid)
{
    int n, ndim;
    cgsize_t np, dims[CGIO_MAX_DIMENSIONS];

    if (cgio_get_dimensions (cgioNum, nodeid, &ndim, dims)) {
        get_error (interp, "cgio_get_dimensions");
        return -1;
    }
    if (ndim == 0)
        return 0;
    for (np = 1, n = 0; n < ndim; n++)
        np *= dims[n];
    return np;
}

/*---------- data_type ---------------------------------------------
 * return number of data values for a node
 *------------------------------------------------------------------*/

static struct DataType * data_type (Tcl_Interp *interp, double nodeid)
{
    int n;
    char type[CGIO_MAX_DATATYPE_LENGTH+1];

    if (cgio_get_data_type (cgioNum, nodeid, type)) {
        get_error (interp, "cgio_get_data_type");
        return NULL;
    }
    for (n = 0; n < CGIO_MAX_DATATYPE_LENGTH && type[n]; n++) {
        if (islower (type[n]))
            type[n] = toupper (type[n]);
    }
    for (n = 0; n < NumDataTypes; n++) {
        if (0 == strncmp (type, DataList[n].name, 2))
            return (&DataList[n]);
    }
    Tcl_AppendResult (interp, "unrecognized data type:", type, NULL);
    return NULL;
}

/*---------- convert_C1 --------------------------------------------
 * convert data from C1
 *------------------------------------------------------------------*/

static char *convert_C1 (size_t np, char *olddata, struct DataType *dtnew)
{
    size_t n;
    C1 *d = (C1 *)olddata;
    char *newdata;

    if (np < 1 || dtnew->size == 0) return NULL;
    if (dtnew->type == C1data) return olddata;
    newdata = (char *) calloc (np, dtnew->size);
    if (NULL == newdata) return NULL;

    switch (dtnew->type) {
        case I4data: {
            I4 *i = (I4 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I4)*d++;
            break;
        }
        case U4data: {
            U4 *u = (U4 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U4)*d++;
            break;
        }
        case I8data: {
            I8 *i = (I8 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I8)*d++;
            break;
        }
        case U8data: {
            U8 *u = (U8 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U8)*d++;
            break;
        }
        case R4data: {
            R4 *r = (R4 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R4)*d++;
            break;
        }
        case R8data: {
            R8 *r = (R8 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R8)*d++;
            break;
        }
        case X4data: {
            X4 *x = (X4 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X4)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case X8data: {
            X8 *x = (X8 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X8)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case B1data: {
            B1 *b = (B1 *)newdata;
            for (n = 0; n < np; n++)
                *b++ = (B1)*d++;
            break;
        }
        default:
            free (newdata);
            newdata = NULL;
            break;
    }
    return newdata;
}

/*---------- convert_B1 --------------------------------------------
 * convert data from B1
 *------------------------------------------------------------------*/

static char *convert_B1 (size_t np, char *olddata, struct DataType *dtnew)
{
    size_t n;
    B1 *d = (B1 *)olddata;
    char *newdata;

    if (np < 1 || dtnew->size == 0) return NULL;
    if (dtnew->type == B1data) return olddata;
    newdata = (char *) calloc (np, dtnew->size);
    if (NULL == newdata) return NULL;

    switch (dtnew->type) {
        case I4data: {
            I4 *i = (I4 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I4)*d++;
            break;
        }
        case U4data: {
            U4 *u = (U4 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U4)*d++;
            break;
        }
        case I8data: {
            I8 *i = (I8 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I8)*d++;
            break;
        }
        case U8data: {
            U8 *u = (U8 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U8)*d++;
            break;
        }
        case R4data: {
            R4 *r = (R4 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R4)*d++;
            break;
        }
        case R8data: {
            R8 *r = (R8 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R8)*d++;
            break;
        }
        case X4data: {
            X4 *x = (X4 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X4)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case X8data: {
            X8 *x = (X8 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X8)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case C1data: {
            C1 *c = (C1 *)newdata;
            for (n = 0; n < np; n++)
                *c++ = (C1)*d++;
            break;
        }
        default:
            free (newdata);
            newdata = NULL;
            break;
    }
    return newdata;
}

/*---------- convert_I4 --------------------------------------------
 * convert data from I4
 *------------------------------------------------------------------*/

static char *convert_I4 (size_t np, char *olddata, struct DataType *dtnew)
{
    size_t n;
    I4 *d = (I4 *)olddata;
    char *newdata;

    if (np < 1 || dtnew->size == 0) return NULL;
    if (dtnew->type == I4data) return olddata;
    newdata = (char *) calloc (np, dtnew->size);
    if (NULL == newdata) return NULL;

    switch (dtnew->type) {
        case U4data: {
            U4 *u = (U4 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U4)*d++;
            break;
        }
        case I8data: {
            I8 *i = (I8 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I8)*d++;
            break;
        }
        case U8data: {
            U8 *u = (U8 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U8)*d++;
            break;
        }
        case R4data: {
            R4 *r = (R4 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R4)*d++;
            break;
        }
        case R8data: {
            R8 *r = (R8 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R8)*d++;
            break;
        }
        case X4data: {
            X4 *x = (X4 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X4)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case X8data: {
            X8 *x = (X8 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X8)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case C1data: {
            C1 *c = (C1 *)newdata;
            for (n = 0; n < np; n++)
                *c++ = (C1)*d++;
            break;
        }
        case B1data: {
            B1 *b = (B1 *)newdata;
            for (n = 0; n < np; n++)
                *b++ = (B1)*d++;
            break;
        }
        default:
            free (newdata);
            newdata = NULL;
            break;
    }
    return newdata;
}

/*---------- convert_U4 --------------------------------------------
 * convert data from U4
 *------------------------------------------------------------------*/

static char *convert_U4 (size_t np, char *olddata, struct DataType *dtnew)
{
    size_t n;
    U4 *d = (U4 *)olddata;
    char *newdata;

    if (np < 1 || dtnew->size == 0) return NULL;
    if (dtnew->type == U4data) return olddata;
    newdata = (char *) calloc (np, dtnew->size);
    if (NULL == newdata) return NULL;

    switch (dtnew->type) {
        case I4data: {
            I4 *i = (I4 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I4)*d++;
            break;
        }
        case I8data: {
            I8 *i = (I8 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I8)*d++;
            break;
        }
        case U8data: {
            U8 *u = (U8 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U8)*d++;
            break;
        }
        case R4data: {
            R4 *r = (R4 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R4)*d++;
            break;
        }
        case R8data: {
            R8 *r = (R8 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R8)*d++;
            break;
        }
        case X4data: {
            X4 *x = (X4 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X4)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case X8data: {
            X8 *x = (X8 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X8)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case C1data: {
            C1 *c = (C1 *)newdata;
            for (n = 0; n < np; n++)
                *c++ = (C1)*d++;
            break;
        }
        case B1data: {
            B1 *b = (B1 *)newdata;
            for (n = 0; n < np; n++)
                *b++ = (B1)*d++;
            break;
        }
        default:
            free (newdata);
            newdata = NULL;
            break;
    }
    return newdata;
}

/*---------- convert_I8 --------------------------------------------
 * convert data from I8
 *------------------------------------------------------------------*/

static char *convert_I8 (size_t np, char *olddata, struct DataType *dtnew)
{
    size_t n;
    I8 *d = (I8 *)olddata;
    char *newdata;

    if (np < 1 || dtnew->size == 0) return NULL;
    if (dtnew->type == I8data) return olddata;
    newdata = (char *) calloc (np, dtnew->size);
    if (NULL == newdata) return NULL;

    switch (dtnew->type) {
        case I4data: {
            I4 *i = (I4 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I4)*d++;
            break;
        }
        case U4data: {
            U4 *u = (U4 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U4)*d++;
            break;
        }
        case U8data: {
            U8 *u = (U8 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U8)*d++;
            break;
        }
        case R4data: {
            R4 *r = (R4 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R4)*d++;
            break;
        }
        case R8data: {
            R8 *r = (R8 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R8)*d++;
            break;
        }
        case X4data: {
            X4 *x = (X4 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X4)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case X8data: {
            X8 *x = (X8 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X8)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case C1data: {
            C1 *c = (C1 *)newdata;
            for (n = 0; n < np; n++)
                *c++ = (C1)*d++;
            break;
        }
        case B1data: {
            B1 *b = (B1 *)newdata;
            for (n = 0; n < np; n++)
                *b++ = (B1)*d++;
            break;
        }
        default:
            free (newdata);
            newdata = NULL;
            break;
    }
    return newdata;
}

/*---------- convert_U8 --------------------------------------------
 * convert data from U8
 *------------------------------------------------------------------*/

static char *convert_U8 (size_t np, char *olddata, struct DataType *dtnew)
{
    size_t n;
    U8 *d = (U8 *)olddata;
    char *newdata;

    if (np < 1 || dtnew->size == 0) return NULL;
    if (dtnew->type == U8data) return olddata;
    newdata = (char *) calloc (np, dtnew->size);
    if (NULL == newdata) return NULL;

    switch (dtnew->type) {
        case I4data: {
            I4 *i = (I4 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I4)*d++;
            break;
        }
        case U4data: {
            U4 *u = (U4 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U4)*d++;
            break;
        }
        case I8data: {
            I8 *i = (I8 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I8)*d++;
            break;
        }
        case R4data: {
            R4 *r = (R4 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R4)*d++;
            break;
        }
        case R8data: {
            R8 *r = (R8 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R8)*d++;
            break;
        }
        case X4data: {
            X4 *x = (X4 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X4)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case X8data: {
            X8 *x = (X8 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X8)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case C1data: {
            C1 *c = (C1 *)newdata;
            for (n = 0; n < np; n++)
                *c++ = (C1)*d++;
            break;
        }
        case B1data: {
            B1 *b = (B1 *)newdata;
            for (n = 0; n < np; n++)
                *b++ = (B1)*d++;
            break;
        }
        default:
            free (newdata);
            newdata = NULL;
            break;
    }
    return newdata;
}

/*---------- convert_R4 --------------------------------------------
 * convert data from R4
 *------------------------------------------------------------------*/

static char *convert_R4 (size_t np, char *olddata, struct DataType *dtnew)
{
    size_t n;
    R4 *d = (R4 *)olddata;
    char *newdata;

    if (np < 1 || dtnew->size == 0) return NULL;
    if (dtnew->type == R4data) return olddata;
    newdata = (char *) calloc (np, dtnew->size);
    if (NULL == newdata) return NULL;

    switch (dtnew->type) {
        case I4data: {
            I4 *i = (I4 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I4)*d++;
            break;
        }
        case U4data: {
            U4 *u = (U4 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U4)*d++;
            break;
        }
        case I8data: {
            I8 *i = (I8 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I8)*d++;
            break;
        }
        case U8data: {
            U8 *u = (U8 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U8)*d++;
            break;
        }
        case R8data: {
            R8 *r = (R8 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R8)*d++;
            break;
        }
        case X4data: {
            X4 *x = (X4 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X4)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case X8data: {
            X8 *x = (X8 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X8)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case C1data: {
            C1 *c = (C1 *)newdata;
            for (n = 0; n < np; n++)
                *c++ = (C1)*d++;
            break;
        }
        case B1data: {
            B1 *b = (B1 *)newdata;
            for (n = 0; n < np; n++)
                *b++ = (B1)*d++;
            break;
        }
        default:
            free (newdata);
            newdata = NULL;
            break;
    }
    return newdata;
}

/*---------- convert_R8 --------------------------------------------
 * convert data from R8
 *------------------------------------------------------------------*/

static char *convert_R8 (size_t np, char *olddata, struct DataType *dtnew)
{
    size_t n;
    R8 *d = (R8 *)olddata;
    char *newdata;

    if (np < 1 || dtnew->size == 0) return NULL;
    if (dtnew->type == R8data) return olddata;
    newdata = (char *) calloc (np, dtnew->size);
    if (NULL == newdata) return NULL;

    switch (dtnew->type) {
        case I4data: {
            I4 *i = (I4 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I4)*d++;
            break;
        }
        case U4data: {
            U4 *u = (U4 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U4)*d++;
            break;
        }
        case I8data: {
            I8 *i = (I8 *)newdata;
            for (n = 0; n < np; n++)
                *i++ = (I8)*d++;
            break;
        }
        case U8data: {
            U8 *u = (U8 *)newdata;
            for (n = 0; n < np; n++)
                *u++ = (U8)*d++;
            break;
        }
        case R4data: {
            R4 *r = (R4 *)newdata;
            for (n = 0; n < np; n++)
                *r++ = (R4)*d++;
            break;
        }
        case X4data: {
            X4 *x = (X4 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X4)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case X8data: {
            X8 *x = (X8 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X8)*d++;
                *x++ = 0.0;
            }
            break;
        }
        case C1data: {
            C1 *c = (C1 *)newdata;
            for (n = 0; n < np; n++)
                *c++ = (C1)*d++;
            break;
        }
        case B1data: {
            B1 *b = (B1 *)newdata;
            for (n = 0; n < np; n++)
                *b++ = (B1)*d++;
            break;
        }
        default:
            free (newdata);
            newdata = NULL;
            break;
    }
    return newdata;
}

/*---------- convert_X4 --------------------------------------------
 * convert data from X4
 *------------------------------------------------------------------*/

static char *convert_X4 (size_t np, char *olddata, struct DataType *dtnew)
{
    size_t n;
    X4 *d = (X4 *)olddata;
    char *newdata;

    if (np < 1 || dtnew->size == 0) return NULL;
    if (dtnew->type == X4data) return olddata;
    newdata = (char *) calloc (np, dtnew->size);
    if (NULL == newdata) return NULL;

    switch (dtnew->type) {
        case I4data: {
            I4 *i = (I4 *)newdata;
            for (n = 0; n < np; n++, d++)
                *i++ = (I4)*d++;
            break;
        }
        case U4data: {
            U4 *u = (U4 *)newdata;
            for (n = 0; n < np; n++, d++)
                *u++ = (U4)*d++;
            break;
        }
        case I8data: {
            I8 *i = (I8 *)newdata;
            for (n = 0; n < np; n++, d++)
                *i++ = (I8)*d++;
            break;
        }
        case U8data: {
            U8 *u = (U8 *)newdata;
            for (n = 0; n < np; n++, d++)
                *u++ = (U8)*d++;
            break;
        }
        case R4data: {
            R4 *r = (R4 *)newdata;
            for (n = 0; n < np; n++, d++)
                *r++ = (R4)*d++;
            break;
        }
        case R8data: {
            R8 *r = (R8 *)newdata;
            for (n = 0; n < np; n++, d++)
                *r++ = (R8)*d++;
            break;
        }
        case X8data: {
            X8 *x = (X8 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X8)*d++;
                *x++ = (X8)*d++;
            }
            break;
        }
        case C1data: {
            C1 *c = (C1 *)newdata;
            for (n = 0; n < np; n++, d++)
                *c++ = (C1)*d++;
            break;
        }
        case B1data: {
            B1 *b = (B1 *)newdata;
            for (n = 0; n < np; n++, d++)
                *b++ = (B1)*d++;
            break;
        }
        default:
            free (newdata);
            newdata = NULL;
            break;
    }
    return newdata;
}

/*---------- convert_X8 --------------------------------------------
 * convert data from C1
 *------------------------------------------------------------------*/

static char *convert_X8 (size_t np, char *olddata, struct DataType *dtnew)
{
    size_t n;
    X8 *d = (X8 *)olddata;
    char *newdata;

    if (np < 1 || dtnew->size == 0) return NULL;
    if (dtnew->type == X8data) return olddata;
    newdata = (char *) calloc (np, dtnew->size);
    if (NULL == newdata) return NULL;

    switch (dtnew->type) {
        case I4data: {
            I4 *i = (I4 *)newdata;
            for (n = 0; n < np; n++, d++)
                *i++ = (I4)*d++;
            break;
        }
        case U4data: {
            U4 *u = (U4 *)newdata;
            for (n = 0; n < np; n++, d++)
                *u++ = (U4)*d++;
            break;
        }
        case I8data: {
            I8 *i = (I8 *)newdata;
            for (n = 0; n < np; n++, d++)
                *i++ = (I8)*d++;
            break;
        }
        case U8data: {
            U8 *u = (U8 *)newdata;
            for (n = 0; n < np; n++, d++)
                *u++ = (U8)*d++;
            break;
        }
        case R4data: {
            R4 *r = (R4 *)newdata;
            for (n = 0; n < np; n++, d++)
                *r++ = (R4)*d++;
            break;
        }
        case R8data: {
            R8 *r = (R8 *)newdata;
            for (n = 0; n < np; n++, d++)
                *r++ = (R8)*d++;
            break;
        }
        case X4data: {
            X4 *x = (X4 *)newdata;
            for (n = 0; n < np; n++) {
                *x++ = (X4)*d++;
                *x++ = (X4)*d++;
            }
            break;
        }
        case C1data: {
            C1 *c = (C1 *)newdata;
            for (n = 0; n < np; n++, d++)
                *c++ = (C1)*d++;
            break;
        }
        case B1data: {
            B1 *b = (B1 *)newdata;
            for (n = 0; n < np; n++, d++)
                *b++ = (B1)*d++;
            break;
        }
        default:
            free (newdata);
            newdata = NULL;
            break;
    }
    return newdata;
}

/*===================================================================
 *           tcl routines
 *===================================================================*/

/*---------- CGNSversion -----------------------------------------------
 * get CGNS library version
 *----------------------------------------------------------------------*/

static int CGNSversion (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    char version[33];

    Tcl_ResetResult (interp);
    if (argc != 1) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], "\"", NULL);
        return TCL_ERROR;
    }
    sprintf (version, "%g", CGNS_DOTVERS);
    Tcl_AppendResult (interp, version, NULL);
    return TCL_OK;
}

/*---------- CGNSfile --------------------------------------------------
 * detect file type
 *----------------------------------------------------------------------*/

static int CGNSfile (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int file_type;

    Tcl_ResetResult (interp);
    if (argc != 2) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " filename\"", NULL);
        return TCL_ERROR;
    }
    if (cgio_check_file (argv[1], &file_type))
        return (get_error (interp, "cgio_check_file"));
    if (file_type == CGIO_FILE_ADF)
        Tcl_SetResult (interp, "adf", TCL_STATIC);
    else if (file_type == CGIO_FILE_HDF5)
        Tcl_SetResult (interp, "hdf5", TCL_STATIC);
    else {
        Tcl_SetResult (interp, "unknown file type", TCL_STATIC);
        return TCL_ERROR;
    }
    return TCL_OK;
}

/*---------- CGIOsupported ----------------------------------------------
 * get supported file types
 *----------------------------------------------------------------------*/

static int CGIOsupported (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int type = CGIO_FILE_NONE;

    Tcl_ResetResult (interp);
    if (argc < 1 || argc > 2) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " ?type?\"", NULL);
        return TCL_ERROR;
    }

    if (argc == 2) {
        switch (argv[1][0]) {
            case 'a':
            case 'A':
                type = CGIO_FILE_ADF;
                break;
            case 'h':
            case 'H':
                type = CGIO_FILE_HDF5;
                break;
            default:
                Tcl_SetResult (interp, "unknown file type", TCL_STATIC);
                return TCL_ERROR;
        }
    }
    if ((type == CGIO_FILE_NONE || type == CGIO_FILE_ADF) &&
            cgio_is_supported (CGIO_FILE_ADF) == 0)
        Tcl_AppendElement (interp, "adf");
    if ((type == CGIO_FILE_NONE || type == CGIO_FILE_HDF5) &&
            cgio_is_supported (CGIO_FILE_HDF5) == 0)
        Tcl_AppendElement (interp, "hdf5");
    return TCL_OK;
}

/*---------- CGIOversion ------------------------------------------------
 * get CGIO library version
 *----------------------------------------------------------------------*/

static int CGIOversion (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    char version[CGIO_MAX_VERSION_LENGTH+1];

    Tcl_ResetResult (interp);
    if (argc != 1) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], "\"", NULL);
        return TCL_ERROR;
    }

    if (cgioNum > 0) {
        if (cgio_library_version (cgioNum, version))
            return (get_error (interp, "cgio_library_version"));
        Tcl_AppendResult (interp, version, NULL);
    }
    return TCL_OK;
}

/*---------- CGIOopen ---------------------------------------------------
 * open CGIO database
 *----------------------------------------------------------------------*/

static int CGIOopen (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int mode, type, cgio_num;
    char rootname[CGIO_MAX_NAME_LENGTH+1];
    double root_id;

    Tcl_ResetResult (interp);
    if (argc < 2 || argc > 4) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " filename ?mode? ?type?\"", NULL);
        return TCL_ERROR;
    }

    if (cgioNum != 0) {
        cgio_close_file (cgioNum);
        cgioNum = 0;
    }

    mode = CGIO_MODE_READ;
    type = CGIO_FILE_NONE;
    if (argc > 2) {
        switch (argv[2][0]) {
            case 'r':
            case 'R':
                mode = CGIO_MODE_READ;
                break;
            case 'w':
            case 'W':
            case 'n':
            case 'N':
                mode = CGIO_MODE_WRITE;
                break;
            case 'm':
            case 'M':
            case 'a':
            case 'A':
            case 'o':
            case 'O':
                mode = CGIO_MODE_MODIFY;
                break;
            default:
                Tcl_SetResult (interp, "unknown file mode", TCL_STATIC);
                return TCL_ERROR;
        }
        if (argc > 3) {
            switch (argv[3][0]) {
                case 'a':
                case 'A':
                    type = CGIO_FILE_ADF;
                    break;
                case 'h':
                case 'H':
                    type = CGIO_FILE_HDF5;
                    break;
                default:
                    Tcl_SetResult (interp, "unknown file type", TCL_STATIC);
                    return TCL_ERROR;
            }
        }
    }

    if (cgio_open_file (argv[1], mode, type, &cgio_num))
        return (get_error (interp, "cgio_open_file"));

    if (cgio_get_root_id (cgio_num, &root_id))
        return (get_error (interp, "cgio_get_root_id"));
    if (cgio_get_name (cgio_num, root_id, rootname))
        return (get_error (interp, "cgio_get_name"));
    cgioNum = cgio_num;
    RootID = root_id;
    Tcl_AppendResult (interp, rootname, NULL);
    return TCL_OK;
}

/*---------- CGIOsave ---------------------------------------------------
 * save CGIO database
 *----------------------------------------------------------------------*/

static int CGIOsave (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int err, type, cgio_num;

    Tcl_ResetResult (interp);
    if (argc != 3) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " filename type\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    switch (argv[2][0]) {
        case 'a':
        case 'A':
            type = CGIO_FILE_ADF;
            break;
        case 'h':
        case 'H':
            type = CGIO_FILE_HDF5;
            break;
        default:
            Tcl_SetResult (interp, "unknown file type", TCL_STATIC);
            return TCL_ERROR;
    }

    if (cgio_open_file (argv[1], CGIO_MODE_WRITE, type, &cgio_num))
        return (get_error (interp, "cgio_open_file"));

    err = cgio_copy_file (cgioNum, cgio_num, 1);
    cgio_close_file (cgio_num);
    if (err)
        return (get_error (interp, "cgio_copy_file"));
    return TCL_OK;
}

/*---------- CGIOclose --------------------------------------------------
 * close CGIO database file
 *----------------------------------------------------------------------*/

static int CGIOclose (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int err;

    Tcl_ResetResult (interp);
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    err = cgio_close_file (cgioNum);
    cgioNum = 0;
    if (err)
        return (get_error (interp, "cgio_close_file"));
    return TCL_OK;
}

/*---------- CGIOcompress -----------------------------------------------
 * remove empty space from CGIO database file
 *----------------------------------------------------------------------*/

static int CGIOcompress (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int err;

    Tcl_ResetResult (interp);
    if (argc != 2) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " filename\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    err = cgio_compress_file (cgioNum, argv[1]);
    if (err)
        return (get_error (interp, "cgio_compress_file"));
    return TCL_OK;
}

/*---------- CGIOnode ---------------------------------------------------
 * get type of node
 *----------------------------------------------------------------------*/

static int CGIOnode (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int len;
    double node_id;

    Tcl_ResetResult (interp);
    if (argc != 2) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }
    if (0 == strcmp (argv[1], "/")) {
        Tcl_SetResult (interp, "root", TCL_STATIC);
        return TCL_OK;
    }
    if (cgio_get_node_id (cgioNum, RootID, argv[1], &node_id) == 0) {
        if (cgio_is_link (cgioNum, node_id, &len) == 0 && len > 0)
            Tcl_SetResult (interp, "link", TCL_STATIC);
        else
            Tcl_SetResult (interp, "node", TCL_STATIC);
    }
    return TCL_OK;
}

/*---------- CGIOname ---------------------------------------------------
 * get/set name for a node
 *----------------------------------------------------------------------*/

static int CGIOname (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    char name[CGIO_MAX_NAME_LENGTH+1];
    double parent_id, node_id;

    Tcl_ResetResult (interp);
    if (argc < 2 || argc > 3) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node ?newname?\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (get_node (interp, argv[1], &node_id))
        return TCL_ERROR;
    if (argc == 3) {
        if (NULL == get_parent (interp, argv[1], &parent_id))
            return TCL_ERROR;
        if (cgio_set_name (cgioNum, parent_id, node_id, argv[2]))
            return (get_error (interp, "cgio_set_name"));
    }
    if (cgio_get_name (cgioNum, node_id, name))
        return (get_error (interp, "cgio_get_name"));
    Tcl_AppendResult (interp, name, NULL);
    return TCL_OK;
}

/*---------- CGIOlabel --------------------------------------------------
 * get/set label for a node
 *----------------------------------------------------------------------*/

static int CGIOlabel (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    char label[CGIO_MAX_LABEL_LENGTH+1];
    double node_id;

    Tcl_ResetResult (interp);
    if (argc < 2 || argc > 3) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node ?newlabel?\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (get_node (interp, argv[1], &node_id))
        return TCL_ERROR;
    if (argc == 3) {
        if (cgio_set_label (cgioNum, node_id, argv[2]))
            return (get_error (interp, "cgio_set_label"));
    }
    if (cgio_get_label (cgioNum, node_id, label))
        return (get_error (interp, "cgio_get_label"));
    Tcl_AppendResult (interp, label, NULL);
    return TCL_OK;
}

/*---------- CGIOtype ---------------------------------------------------
 * get data type for a node
 *----------------------------------------------------------------------*/

static int CGIOtype (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int n, ndim, err;
    cgsize_t np, dims[CGIO_MAX_DIMENSIONS];
    char type[CGIO_MAX_DATATYPE_LENGTH+1];
    double node_id;
    struct DataType *dtnew, *dtold;
    char *olddata, *newdata;

    Tcl_ResetResult (interp);
    if (argc == 1) {
        for (ndim = 0; ndim < NumDataTypes; ndim++)
            Tcl_AppendElement (interp, DataList[ndim].name);
        return TCL_OK;
    }
    if (argc > 3) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node ?newtype?\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (get_node (interp, argv[1], &node_id))
        return TCL_ERROR;

    if (argc > 2) {

        /* get new data type */

        strncpy (type, argv[2], CGIO_MAX_DATATYPE_LENGTH);
        type[CGIO_MAX_DATATYPE_LENGTH] = 0;
        for (n = 0; n < CGIO_MAX_DATATYPE_LENGTH && type[n]; n++) {
            if (islower (type[n]))
                type[n] = toupper (type[n]);
        }
        dtnew = NULL;
        for (n = 0; n < NumDataTypes; n++) {
            if (0 == strncmp (type, DataList[n].name, 2)) {
                dtnew = &DataList[n];
                break;
            }
        }
        if (NULL == dtnew) {
            Tcl_AppendResult (interp, "unrecognized data type", NULL);
            return TCL_ERROR;
        }

        /* get old data type */

        dtold = data_type (interp, node_id);
        if (NULL == dtold)
            return TCL_ERROR;

        /* convert the data */

        if (dtnew->type != dtold->type) {

            /* get current data size */

            np = 0;
            if (cgio_get_dimensions (cgioNum, node_id, &ndim, dims))
                return (get_error (interp, "cgio_get_dimensions"));
            if (ndim > 0) {
                for (np = 1, n = 0; n < ndim; n++)
                    np *= dims[n];
            }
            newdata = NULL;

            if (np > 0 && dtold->size > 0) {
                olddata = (char *) malloc ((size_t)(np * dtold->size));
                if (NULL == olddata) {
                    Tcl_AppendResult (interp, "malloc failed for data", NULL);
                    return TCL_ERROR;
                }
                if (cgio_read_all_data (cgioNum, node_id, olddata)) {
                    free (olddata);
                    return (get_error (interp, "cgio_read_all_data"));
                }
                switch (dtold->type) {
                    case C1data:
                        newdata = convert_C1 ((size_t)np, olddata, dtnew);
                        break;
                    case B1data:
                        newdata = convert_B1 ((size_t)np, olddata, dtnew);
                        break;
                    case I4data:
                        newdata = convert_I4 ((size_t)np, olddata, dtnew);
                        break;
                    case U4data:
                        newdata = convert_U4 ((size_t)np, olddata, dtnew);
                        break;
                    case I8data:
                        newdata = convert_I8 ((size_t)np, olddata, dtnew);
                        break;
                    case U8data:
                        newdata = convert_U8 ((size_t)np, olddata, dtnew);
                        break;
                    case R4data:
                        newdata = convert_R4 ((size_t)np, olddata, dtnew);
                        break;
                    case R8data:
                        newdata = convert_R8 ((size_t)np, olddata, dtnew);
                        break;
                    case X4data:
                        newdata = convert_X4 ((size_t)np, olddata, dtnew);
                        break;
                    case X8data:
                        newdata = convert_X8 ((size_t)np, olddata, dtnew);
                        break;
                }
                if (newdata != olddata)
                    free (olddata);
            }

            /* rewrite the data */

            if (NULL == newdata) ndim = 0;
            if (cgio_set_dimensions (cgioNum, node_id,
                    dtnew->name, ndim, dims))
                return (get_error (interp, "cgio_set_dimensions"));
            if (NULL != newdata) {
                err = cgio_write_all_data (cgioNum, node_id, newdata);
                free (newdata);
                if (err)
                    return (get_error (interp, "cgio_write_all_data"));
            }
        }
    }
    if (cgio_get_data_type (cgioNum, node_id, type))
        return (get_error (interp, "cgio_get_data_type"));
    Tcl_AppendResult (interp, type, NULL);
    return TCL_OK;
}

/*---------- CGIOdimensions ---------------------------------------------
 * get/set data dimensions for a node
 *----------------------------------------------------------------------*/

static int CGIOdimensions (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int n;
    int ndim;
    cgsize_t dims[CGIO_MAX_DIMENSIONS];
    double node_id;
    CONST char **args;
    char str[65];

    Tcl_ResetResult (interp);
    if (argc < 2 || argc > 3) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node ?newdimensions?\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (get_node (interp, argv[1], &node_id))
        return TCL_ERROR;
    if (argc > 2) {
        if (cgio_get_data_type (cgioNum, node_id, str))
            return (get_error (interp, "cgio_get_data_type"));
        if (TCL_OK != Tcl_SplitList (interp, argv[2], &ndim, &args))
            return TCL_ERROR;
        if (ndim > CGIO_MAX_DIMENSIONS) {
            Tcl_Free ((char *)args);
            Tcl_AppendResult (interp, "invalid number of dimensions", NULL);
            return TCL_ERROR;
        }
        if (ndim) {
            for (n = 0; n < ndim; n++)
                dims[n] = atoi (args[n]);
            Tcl_Free ((char *)args);
        }
        if (cgio_set_dimensions (cgioNum, node_id, str, ndim, dims))
            return (get_error (interp, "cgio_set_dimensions"));
    }
    if (cgio_get_dimensions (cgioNum, node_id, &ndim, dims))
        return (get_error (interp, "cgio_get_dimensions"));
    if (ndim > 0) {
        for (n = 0; n < ndim; n++) {
            sprintf (str, "%ld", (long)dims[n]);
            Tcl_AppendElement (interp, str);
        }
    }
    return TCL_OK;
}

/*---------- CGIOsize ---------------------------------------------------
 * get number of bytes of data for a node
 *----------------------------------------------------------------------*/

static int CGIOsize (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    cgsize_t np;
    char str[65];
    double node_id;
    struct DataType *type;

    Tcl_ResetResult (interp);
    if (argc != 2) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (get_node (interp, argv[1], &node_id))
        return TCL_ERROR;

    np = data_size (interp, node_id);
    if (np == -1)
        return TCL_ERROR;
    type = data_type (interp, node_id);
    if (NULL == type)
        return TCL_ERROR;
    sprintf (str, "%u", (unsigned)np * (unsigned)type->bytes);
    Tcl_AppendResult (interp, str, NULL);
    return TCL_OK;
}

/*---------- CGIOread ---------------------------------------------------
 * read node data
 *----------------------------------------------------------------------*/

static int CGIOread (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int n;
    cgsize_t np;
    char *values, str[65];
    double node_id;
    struct DataType *dt;

    Tcl_ResetResult (interp);
    if (argc < 2 || argc > 14) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node ?range1 range2 ...?\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (get_node (interp, argv[1], &node_id))
        return TCL_ERROR;
    np = data_size (interp, node_id);
    if (np == -1)
        return TCL_ERROR;
    dt = data_type (interp, node_id);
    if (NULL == dt)
        return TCL_ERROR;
    if (np == 0 || dt->size == 0)
        return TCL_OK;

    values = (char *) malloc ((unsigned)np * (unsigned)dt->size + 1);
    if (NULL == values) {
        Tcl_AppendResult (interp, "malloc failed for data", NULL);
        return TCL_ERROR;
    }

    if (cgio_read_all_data (cgioNum, node_id, values)) {
        free (values);
        return (get_error (interp, "cgio_read_all_data"));
    }

    if (dt->type == C1data) {
        values[np] = 0;
        Tcl_AppendResult (interp, values, NULL);
    }
    else if (dt->type == B1data) {
        B1 *u = (B1 *)values;
        for (n = 0; n < np; n++, u++) {
            sprintf (str, "%d", (int)*u);
            Tcl_AppendElement (interp, str);
        }
    }
    else if (dt->type == I4data) {
        I4 *i = (I4 *)values;
        for (n = 0; n < np; n++, i++) {
            sprintf (str, "%d", *i);
            Tcl_AppendElement (interp, str);
        }
    }
    else if (dt->type == U4data) {
        U4 *u = (U4 *)values;
        for (n = 0; n < np; n++, u++) {
            sprintf (str, "%u", *u);
            Tcl_AppendElement (interp, str);
        }
    }
    else if (dt->type == I8data) {
        I8 *i = (I8 *)values;
        for (n = 0; n < np; n++, i++) {
            sprintf (str, "%ld", *i);
            Tcl_AppendElement (interp, str);
        }
    }
    else if (dt->type == U8data) {
        U8 *u = (U8 *)values;
        for (n = 0; n < np; n++, u++) {
            sprintf (str, "%lu", *u);
            Tcl_AppendElement (interp, str);
        }
    }
    else if (dt->type == R4data) {
        R4 *r = (R4 *)values;
        for (n = 0; n < np; n++, r++) {
            sprintf (str, "%g", *r);
            Tcl_AppendElement (interp, str);
        }
    }
    else if (dt->type == R8data) {
        R8 *r = (R8 *)values;
        for (n = 0; n < np; n++, r++) {
            sprintf (str, "%g", *r);
            Tcl_AppendElement (interp, str);
        }
    }
    else if (dt->type == X4data) {
        X4 *r = (X4 *)values;
        for (n = 0; n < np; n++, r++) {
            sprintf (str, "%g %g", *r, *(r+1));
            Tcl_AppendElement (interp, str);
            r++;
        }
    }
    else if (dt->type == X8data) {
        X8 *r = (X8 *)values;
        for (n = 0; n < np; n++, r++) {
            sprintf (str, "%g %g", *r, *(r+1));
            Tcl_AppendElement (interp, str);
            r++;
        }
    }
    else {
        free (values);
        Tcl_AppendResult (interp, "internal error - should not happen", NULL);
        return TCL_ERROR;
    }

    free (values);
    return TCL_OK;
}

/*---------- CGIOwrite --------------------------------------------------
 * write node data
 *----------------------------------------------------------------------*/

static int CGIOwrite (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int n, nv, ns;
    int ndim;
    cgsize_t np, dims[CGIO_MAX_DIMENSIONS];
    double node_id;
    CONST char **args;
    char *values, type[CGIO_MAX_DATATYPE_LENGTH+1];
    struct DataType *dt = NULL;

    Tcl_ResetResult (interp);
    if (argc < 3 || argc > 5) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node type ?dimensions? ?data?\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    /* get node ID for node */

    if (get_node (interp, argv[1], &node_id))
        return TCL_ERROR;

    /* get data type */

    strncpy (type, argv[2], CGIO_MAX_DATATYPE_LENGTH);
    type[CGIO_MAX_DATATYPE_LENGTH] = 0;
    for (n = 0; n < CGIO_MAX_DATATYPE_LENGTH && type[n]; n++) {
        if (islower (type[n]))
            type[n] = toupper (type[n]);
    }
    for (n = 0; n < NumDataTypes; n++) {
        if (0 == strncmp (type, DataList[n].name, 2)) {
            dt = &DataList[n];
            break;
        }
    }
    if (dt == NULL) {
        Tcl_AppendResult (interp, "data type not recognized", NULL);
        return TCL_ERROR;
    }

    /* get dimensions */

    ndim = 0;
    args = NULL;
    if (argc > 3 &&
        TCL_OK != Tcl_SplitList (interp, argv[3], &ndim, &args))
        return TCL_ERROR;
    if (ndim == 0) {
        if (cgio_set_dimensions (cgioNum, node_id, dt->name, ndim, dims))
            return (get_error (interp, "cgio_set_dimensions"));
        return TCL_OK;
    }

    if (ndim > CGIO_MAX_DIMENSIONS) {
        Tcl_Free ((char *)args);
        Tcl_AppendResult (interp, "invalid number of dimensions", NULL);
        return TCL_ERROR;
    }
    for (n = 0; n < ndim; n++)
        dims[n] = (cgsize_t) atol (args[n]);
    Tcl_Free ((char *)args);
    for (np = 1, n = 0; n < ndim; n++) {
        if (dims[n] < 1) {
            Tcl_AppendResult (interp, "invalid dimension", NULL);
            return TCL_ERROR;
        }
        np *= dims[n];
    }

    /* create data array */

    if (NULL == (values = (char *) calloc ((size_t)np, dt->size))) {
        Tcl_AppendResult (interp, "malloc failed for data", NULL);
        return TCL_ERROR;
    }

    /* get data */

    if (argc > 4) {
        if (dt->type == C1data && ndim == 1) {
            strncpy (values, argv[4], (size_t)np);
            for (ns = (int)strlen(argv[4]); ns < np; ns++)
                values[ns] = ' ';
        }
        else {
            if (TCL_OK != Tcl_SplitList (interp, argv[4], &nv, &args)) {
                free (values);
                return TCL_ERROR;
            }
            if (nv) {
                if (dt->type == C1data) np /= dims[0];
                if (np > nv) np = nv;

                if (dt->type == C1data) {
                    char *p = values;
                    for (n = 0; n < np; n++) {
                        strncpy (p, args[n], (size_t)dims[0]);
                        for (ns = (int)strlen(args[n]); ns < dims[0]; ns++)
                            p[ns] = ' ';
                        p += dims[0];
                    }
                }
                else if (dt->type == B1data) {
                    B1 *u = (B1 *)values;
                    for (n = 0; n < np; n++, u++)
                        *u = (B1) atoi (args[n]);
                }
                else if (dt->type == I4data) {
                    I4 *i = (I4 *)values;
                    for (n = 0; n < np; n++, i++)
                        *i = (I4) atoi (args[n]);
                }
                else if (dt->type == U4data) {
                    U4 *u = (U4 *)values;
                    for (n = 0; n < np; n++, u++)
                        *u = (U4) atoi (args[n]);
                }
                else if (dt->type == I8data) {
                    I8 *i = (I8 *)values;
                    for (n = 0; n < np; n++, i++)
                        *i = (I8) atol (args[n]);
                }
                else if (dt->type == U8data) {
                    U8 *u = (U8 *)values;
                    for (n = 0; n < np; n++, u++)
                        *u = (U8) atol (args[n]);
                }
                else if (dt->type == R4data) {
                    R4 *r = (R4 *)values;
                    for (n = 0; n < np; n++, r++)
                        *r = (R4) atof (args[n]);
                }
                else if (dt->type == R8data) {
                    R8 *r = (R8 *)values;
                    for (n = 0; n < np; n++, r++)
                        *r = (R8) atof (args[n]);
                }
                else if (dt->type == X4data) {
                    float fr, fi;
                    X4 *r = (X4 *)values;
                    for (n = 0; n < np; n++, r++) {
                        if (2 != sscanf (args[n], "%f %f", &fr, &fi)) {
                            fr = (float) atof (args[n]);
                            fi = 0.0;
                        }
                        *r = (R4) fr;
                        *++r = (R4) fi;
                    }
                }
                else if (dt->type == X8data) {
                    double dr, di;
                    X8 *r = (X8 *)values;
                    for (n = 0; n < np; n++, r++) {
                        if (2 != sscanf (args[n], "%lf %lf", &dr, &di)) {
                            dr = atof (args[n]);
                            di = 0.0;
                        }
                        *r = (X8) dr;
                        *++r = (X8) di;
                    }
                }
                else {
                    Tcl_Free ((char *)args);
                    free (values);
                    Tcl_AppendResult (interp, "internal error - should not happen", NULL);
                    return TCL_ERROR;
                }

                Tcl_Free ((char *)args);
            }
        }
    }

    if (cgio_set_dimensions (cgioNum, node_id, dt->name, ndim, dims)) {
        free (values);
        return (get_error (interp, "cgio_set_dimensions"));
    }
    n = cgio_write_all_data (cgioNum, node_id, values);
    free (values);
    if (n)
        return (get_error (interp, "cgio_write_all_data"));
    return TCL_OK;
}

/*---------- CGIOlink ---------------------------------------------------
 * create/retrieve link of a node
 *----------------------------------------------------------------------*/

static int CGIOlink (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int file_len, name_len;
    char *node, *name_in_file, *file_name;
    double parent_id, node_id;

    Tcl_ResetResult (interp);
    if (argc < 2 || argc > 4) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node ?linknode? ?linkfile?\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (argc > 2) {
        if (NULL == (node = get_parent (interp, argv[1], &parent_id)))
            return TCL_ERROR;
        name_in_file = argv[2];
        if (argc > 3)
            file_name = argv[3];
        else
            file_name = NULL;
        if (cgio_create_link (cgioNum, parent_id, node, file_name,
                name_in_file, &node_id))
            return (get_error (interp, "cgio_create_link"));
    }
    else {
        if (get_node (interp, argv[1], &node_id))
            return TCL_ERROR;
    }
    if (cgio_is_link (cgioNum, node_id, &name_len))
        return (get_error (interp, "cgio_is_link"));
    if (name_len > 0) {
        if (cgio_link_size (cgioNum, node_id, &file_len, &name_len))
            return (get_error (interp, "cgio_link_size"));
        file_name = (char *) malloc (file_len + name_len + 2);
        if (NULL == file_name) {
            Tcl_AppendResult (interp, "malloc failed for link", NULL);
            return TCL_ERROR;
        }
        name_in_file = file_name + file_len + 1;
        if (cgio_get_link (cgioNum, node_id, file_name, name_in_file)) {
            free (file_name);
            return (get_error (interp, "cgio_get_link"));
        }
        Tcl_AppendElement (interp, name_in_file);
        Tcl_AppendElement (interp, file_name);
        free (file_name);
    }
    return TCL_OK;
}

/*---------- CGIOchildren -----------------------------------------------
 * retrieve children of a node
 *----------------------------------------------------------------------*/

static int CGIOchildren (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int n, nchildren, len;
    char name[CGIO_MAX_NAME_LENGTH+1];
    double node_id, *ids;

    Tcl_ResetResult (interp);
    if (2 != argc) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (get_node (interp, argv[1], &node_id))
        return TCL_ERROR;
    if (cgio_number_children (cgioNum, node_id, &nchildren))
        return (get_error (interp, "cgio_number_children"));
    if (nchildren < 1)
        return TCL_OK;

    ids = (double *) malloc (nchildren * sizeof(double));
    if (NULL == ids) {
        Tcl_AppendResult (interp, "malloc failed for IDs", NULL);
        return TCL_ERROR;
    }
    if (cgio_children_ids (cgioNum, node_id, 1, nchildren, &len, ids)) {
        free (ids);
        Tcl_ResetResult (interp);
        return (get_error (interp, "cgio_children_ids"));
    }
    for (n = 0; n < nchildren; n++) {
        if (cgio_get_name (cgioNum, ids[n], name)) {
            free (ids);
            Tcl_ResetResult (interp);
            return (get_error (interp, "cgio_get_name"));
        }
        Tcl_AppendElement (interp, name);
    }
    free (ids);
    return TCL_OK;
}

/*---------- CGIOnumchild -----------------------------------------------
 * returns number of children for a node
 *----------------------------------------------------------------------*/

static int CGIOnumchild (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int nchildren;
    double node_id;
    char buf[33];

    Tcl_ResetResult (interp);
    if (2 != argc) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (get_node (interp, argv[1], &node_id))
        return TCL_ERROR;
    if (cgio_number_children (cgioNum, node_id, &nchildren))
        return (get_error (interp, "cgio_number_children"));
    sprintf (buf, "%d", nchildren);
    Tcl_AppendResult (interp, buf, NULL);
    return TCL_OK;
}

/*---------- CGIOchildname ----------------------------------------------
 * retrieve child name for child index
 *----------------------------------------------------------------------*/

static int CGIOchildname (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    int len;
    char name[CGIO_MAX_NAME_LENGTH+1];
    double node_id;

    Tcl_ResetResult (interp);
    if (3 != argc) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node childnum\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (get_node (interp, argv[1], &node_id))
        return TCL_ERROR;
    if (cgio_children_names (cgioNum, node_id, atoi(argv[2]), 1,
            CGIO_MAX_NAME_LENGTH+1, &len, name))
        return (get_error (interp, "cgio_children_names"));
    Tcl_AppendResult (interp, name, NULL);
    return TCL_OK;
}

/*---------- CGIOcreate -------------------------------------------------
 * create a new node
 *----------------------------------------------------------------------*/

static int CGIOcreate (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    char *node, label[CGIO_MAX_LABEL_LENGTH+1];
    double parent_id, node_id;

    Tcl_ResetResult (interp);
    if (argc < 2 || argc > 3) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node ?label?\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (NULL == (node = get_parent (interp, argv[1], &parent_id)))
        return TCL_ERROR;
    if (cgio_create_node (cgioNum, parent_id, node, &node_id))
        return (get_error (interp, "cgio_create_node"));
    if (argc > 2) {
        if (cgio_set_label (cgioNum, node_id, argv[2]))
            return (get_error (interp, "cgio_set_label"));
    }
    if (cgio_get_label (cgioNum, node_id, label))
        return (get_error (interp, "cgio_get_label"));
    Tcl_AppendResult (interp, label, NULL);
    return TCL_OK;
}

/*---------- CGIOdelete -------------------------------------------------
 * delete a node
 *----------------------------------------------------------------------*/

static int CGIOdelete (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    char *node;
    double parent_id, node_id;

    Tcl_ResetResult (interp);
    if (argc != 2) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (NULL == (node = get_parent (interp, argv[1], &parent_id)))
        return TCL_ERROR;
    if (cgio_get_node_id (cgioNum, parent_id, node, &node_id))
        return (get_error (interp, "cgio_get_node_id"));
    if (cgio_delete_node (cgioNum, parent_id, node_id))
        return (get_error (interp, "cgio_delete_node"));
    return TCL_OK;
}

/*---------- CGIOmove ---------------------------------------------------
 * move a node
 *----------------------------------------------------------------------*/

static int CGIOmove (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    char *node;
    double parent_id, node_id, new_parent_id;

    Tcl_ResetResult (interp);
    if (argc != 3) {
        Tcl_AppendResult (interp, "wrong # args: should be \"",
            argv[0], " node parent\"", NULL);
        return TCL_ERROR;
    }
    if (cgioNum == 0) {
        Tcl_AppendResult (interp, "no database is open", NULL);
        return TCL_ERROR;
    }

    if (NULL == (node = get_parent (interp, argv[1], &parent_id)))
        return TCL_ERROR;
    if (cgio_get_node_id (cgioNum, parent_id, node, &node_id))
        return (get_error (interp, "cgio_get_node_id"));
    if (get_node (interp, argv[2], &new_parent_id))
        return TCL_ERROR;
    if (cgio_move_node (cgioNum, parent_id, node_id, new_parent_id))
        return (get_error (interp, "cgio_move_node"));
    return TCL_OK;
}

#ifdef SINGLE_COMMAND

/*---------- CGIOcommand ----------------------------------------
 * process CGIO command
 *--------------------------------------------------------------*/

static int CGIOcommand (ClientData data, Tcl_Interp *interp,
    int argc, char **argv)
{
    static char usg_msg[] =
        "CGIO version\n"
        "CGIO supported ?type?\n"
        "CGIO open filename ?status? ?format?\n"
        "CGIO save filename format\n"
        "CGIO compress filename\n"
        "CGIO close\n"
        "CGIO node node\n"
        "CGIO name node ?newname?\n"
        "CGIO label node ?newlabel?\n"
        "CGIO type node ?newtype?\n"
        "CGIO dimensions node ?newdimensions?\n"
        "CGIO size node\n"
        "CGIO read node ?range1 range2 ...?\n"
        "CGIO write node type dimensions data\n"
        "CGIO link node ?linknode? ?linkfile?\n"
        "CGIO children node\n"
        "CGIO numchild node\n"
        "CGIO childname node childnum\n"
        "CGIO create node ?label?\n"
        "CGIO delete node\n"
        "CGIO move node newparent\n";

    if (argc < 2) {
        Tcl_SetResult (interp, usg_msg, TCL_STATIC);
        return TCL_ERROR;
    }
    if (0 == strcmp (argv[1], "version"))
        return CGIOversion (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "supported"))
        return CGIOsupported (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "open"))
        return CGIOopen (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "save"))
        return CGIOsave (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "close"))
        return CGIOclose (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "compress"))
        return CGIOcompress (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "node"))
        return CGIOnode (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "name"))
        return CGIOname (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "label"))
        return CGIOlabel (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "type"))
        return CGIOtype (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "dimensions"))
        return CGIOdimensions (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "size"))
        return CGIOsize (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "read"))
        return CGIOread (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "write"))
        return CGIOwrite (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "link"))
        return CGIOlink (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "children"))
        return CGIOchildren (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "numchild"))
        return CGIOnumchild (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "childname"))
        return CGIOchildname (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "create"))
        return CGIOcreate (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "delete"))
        return CGIOdelete (data, interp, argc-1, argv+1);
    if (0 == strcmp (argv[1], "move"))
        return CGIOmove (data, interp, argc-1, argv+1);

    Tcl_SetResult (interp, usg_msg, TCL_STATIC);
    return TCL_ERROR;
}

#endif

/*---------- CGIOtcl_Init ---------------------------------------
 * Initialize and create the commands
 *--------------------------------------------------------------*/

#if defined(_WIN32) && defined(BUILD_DLL)
__declspec(dllexport)
#endif
int CGIOtcl_Init(Tcl_Interp *interp)
{
#ifdef SINGLE_COMMAND
    Tcl_CreateCommand (interp, "CGIO", (Tcl_CmdProc *)CGIOcommand,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
#else
    Tcl_CreateCommand (interp, "CGIOversion", (Tcl_CmdProc *)CGIOversion,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOsupported", (Tcl_CmdProc *)CGIOsupported,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOopen", (Tcl_CmdProc *)CGIOopen,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOsave", (Tcl_CmdProc *)CGIOsave,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOclose", (Tcl_CmdProc *)CGIOclose,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOcompress", (Tcl_CmdProc *)CGIOcompress,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOnode", (Tcl_CmdProc *)CGIOnode,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOname", (Tcl_CmdProc *)CGIOname,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOlabel", (Tcl_CmdProc *)CGIOlabel,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOtype", (Tcl_CmdProc *)CGIOtype,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOdimensions", (Tcl_CmdProc *)CGIOdimensions,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOsize", (Tcl_CmdProc *)CGIOsize,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOread", (Tcl_CmdProc *)CGIOread,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOwrite", (Tcl_CmdProc *)CGIOwrite,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOlink", (Tcl_CmdProc *)CGIOlink,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOchildren", (Tcl_CmdProc *)CGIOchildren,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOnumchild", (Tcl_CmdProc *)CGIOnumchild,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOchildname", (Tcl_CmdProc *)CGIOchildname,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOcreate", (Tcl_CmdProc *)CGIOcreate,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOdelete", (Tcl_CmdProc *)CGIOdelete,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGIOmove", (Tcl_CmdProc *)CGIOmove,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
#endif
    Tcl_CreateCommand (interp, "CGNSversion", (Tcl_CmdProc *)CGNSversion,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGNSfile", (Tcl_CmdProc *)CGNSfile,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    return TCL_OK;
}

