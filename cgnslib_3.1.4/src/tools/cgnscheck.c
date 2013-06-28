/*
 * cgnscheck.c - check CGNS file
 */

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#ifdef _WIN32
# include <io.h>
# define access _access
# define unlink _unlink
#else
# include <unistd.h>
#endif
#include "getargs.h"
#include "hash.h"
#include "cgnslib.h"
#include "cgns_header.h"
#include "cgnames.h"

#if !defined(CGNS_VERSION) || CGNS_VERSION < 3100
# error You need at least CGNS Version 3.1
#endif

#define USE_MID_NODES

static int FileVersion;
static int LibraryVersion = CGNS_VERSION;

static int verbose = 0;
static int nwarn = 0, nerr = 0, totwarn = 0;
static int dowarn = 3, doerr = 1;
static int cgnsfn, cgnsbase, cgnszone;

static int CellDim, PhyDim;
static int BaseClass;
static int *pBaseUnits, BaseUnits[9];

static int BaseIter;
static int NumSteps;

typedef struct {
    float x, y, z;
} VERTEX;

typedef struct {
    cgsize_t e1, e2;
    int f1, f2;
    int nnodes;
    cgsize_t *nodes;
} FACE;

typedef struct {
    char name[33];
    CGNS_ENUMT(ElementType_t) type;
    cgsize_t is, ie;
    int ib;
    cgsize_t nv, ns, ne, nn;
    cgsize_t *elements;
    cgsize_t *parent;
    int rind[2];
    int invalid;
} ELEMSET;

typedef struct {
    char name[33];
    CGNS_ENUMT(ZoneType_t) type;
    int idim;
    cgsize_t dims[3][3];
    cgsize_t nnodes;
    cgsize_t maxnode;
    int dataclass;
    int *punits, units[9];
    int nsets;
    cgsize_t nv, ns, ne, nn;
    ELEMSET *sets;
    HASH *faces;
    cgsize_t nextnodes;
    cgsize_t *extnodes;
} ZONE;

static int MaxZones = 0;
static int NumZones = 0;
static ZONE *Zones;

typedef char CGNSNAME[33];

static int MaxFamily = 0;
static int NumFamily = 0;
static CGNSNAME *Family;

static int MaxFlowSolution = 0;
static int NumFlowSolution = 0;
static CGNSNAME *FlowSolution;

static int MaxGridCoordinate = 0;
static int NumGridCoordinate = 0;
static CGNSNAME *GridCoordinate;

static int MaxArbitraryGrid = 0;
static int NumArbitraryGrid = 0;
static CGNSNAME *ArbitraryGrid;

static int MaxRigidGrid = 0;
static int NumRigidGrid = 0;
static CGNSNAME *RigidGrid;

static int MaxZoneConn = 0;
static int NumZoneConn = 0;
static CGNSNAME *ZoneConn;

static int MaxZoneSubReg = 0;
static int NumZoneSubReg = 0;
static CGNSNAME *ZoneSubReg;

/* command line options */

static char options[] = "vVuUw:e";

static char *usgmsg[] = {
    "usage  : cgnscheck [options] CGNSfile [CGNSoutfile]",
    "options:",
    "   -v        : verbose output",
    "   -V        : more verbose - print descriptors",
    "   -u        : update CGNS file to CGNS Library Version and check",
    "   -U        : update CGNS file to CGNS Library Version only",
    "   -w<level> : warning level output (0 to 3)",
    "   -e        : don't print error",
    NULL
};

/*----------------------------------------------------------------------*/

static void warning (int level, char *format, ...)
{
    va_list arg;
    if (level <= dowarn) {
        va_start (arg, format);
        printf ("WARNING:");
        vprintf (format, arg);
        va_end(arg);
        putchar ('\n');
        nwarn++;
    }
    totwarn++;
}

static void error (char *format, ...)
{
    va_list arg;
    if (doerr) {
        va_start (arg, format);
        printf ("ERROR:");
        vprintf (format, arg);
        va_end(arg);
        putchar ('\n');
    }
    nerr++;
}

static void error_exit(char *func)
{
    printf("CGNSlib ERROR:");
    if (func != NULL && *func)
        printf("%s:", func);
    printf("%s\n", cg_get_error());
    exit(1);
}

static void fatal_error (char *format, ...)
{
    va_list arg;

    fflush(stdout);
    va_start (arg, format);
    fprintf (stderr, "INTERNAL ERROR:");
    vfprintf (stderr, format, arg);
    va_end(arg);
    exit(1);
}

/*----------------------------------------------------------------------*/

static void create_names (int cnt, int *maxcnt, CGNSNAME **namelist)
{
    CGNSNAME *names = *namelist;

    if (cnt > *maxcnt) {
        if (*maxcnt)
            names = (CGNSNAME *) realloc (names, cnt * sizeof(CGNSNAME));
        else
            names = (CGNSNAME *) malloc (cnt * sizeof(CGNSNAME));
        if (names == NULL)
            fatal_error ("malloc failed for cgns name list\n");
        *maxcnt = cnt;
        *namelist = names;
    }
}

/*=======================================================================*/

#define MAX_GOTO_DEPTH CG_MAX_GOTO_DEPTH

static char goLabel[MAX_GOTO_DEPTH][33];
static int goIndex[MAX_GOTO_DEPTH];
static int goDepth = 0;

/*-----------------------------------------------------------------------*/

static void goto_node ()
{
    int n, ier;
    char *labels[MAX_GOTO_DEPTH];

    for (n = 0; n < goDepth; n++)
        labels[n] = goLabel[n];
    ier = cg_golist (cgnsfn, cgnsbase, goDepth, labels, goIndex);
    if (ier) error_exit ("cg_golist");
}

/*-----------------------------------------------------------------------*/

static void go_absolute (char *dsname, ...)
{
    int num;
    char *name = dsname;
    va_list arg;

    va_start (arg, dsname);
    goDepth = 0;
    while (name != NULL) {
        num = va_arg (arg, int);
        if (goDepth == MAX_GOTO_DEPTH)
            fatal_error("maximum depth of goto exceeded\n");
        strncpy (goLabel[goDepth], name, 32);
        goLabel[goDepth][32] = 0;
        goIndex[goDepth] = num;
        goDepth++;
        name = va_arg (arg, char *);
    }
    va_end (arg);
    goto_node();
}

/*-----------------------------------------------------------------------*/

static void go_relative (char *dsname, ...)
{
    int num;
    char *name = dsname;
    va_list arg;

    va_start (arg, dsname);
    while (name != NULL) {
        num = va_arg (arg, int);
        if (num < 1) break;
        if (0 == strcmp (name, "..")) {
            goDepth -= num;
            if (goDepth < 0) goDepth = 0;
        }
        else if (strcmp (name, ".")) {
            if (goDepth == MAX_GOTO_DEPTH)
                fatal_error("maximum depth of goto exceeded\n");
            strncpy (goLabel[goDepth], name, 32);
            goLabel[goDepth][32] = 0;
            goIndex[goDepth++] = num;
        }
        name = va_arg (arg, char *);
    }
    va_end (arg);
    goto_node();
}

/*=======================================================================*/

static int check_node (char *label)
{
    int nchild;
    double pid, *ids;

    if (cgi_posit_id (&pid)) return CG_ERROR;
    if (cgi_get_nodes (pid, label, &nchild, &ids)) return CG_ERROR;
    if (nchild) {
        cg_free (ids);
        return CG_OK;
    }
    return CG_NODE_NOT_FOUND;
}

/*-----------------------------------------------------------------------*/

static int read_gridlocation (CGNS_ENUMT(GridLocation_t) *location)
{
    int ierr = check_node ("GridLocation_t");
    if (ierr == CG_OK)
        return cg_gridlocation_read (location);
    *location = CGNS_ENUMV(Vertex);
    return ierr;
}

/*-----------------------------------------------------------------------*/

static int read_ordinal (int *ordinal)
{
    int ierr = check_node ("Ordinal_t");
    if (ierr == CG_OK)
        return cg_ordinal_read (ordinal);
    return ierr;
}

/*-----------------------------------------------------------------------*/

static int read_rind (int *rind)
{
    int ierr = check_node ("Rind_t");
    if (ierr == CG_OK)
        return cg_rind_read (rind);
    return ierr;
}

/*-----------------------------------------------------------------------*/

static int check_interpolants (void)
{
    int n, na, ndim;
    cgsize_t dims[12];
    CGNS_ENUMT(DataType_t) dtype;
    char name[33];

    if (cg_narrays (&na)) error_exit ("cg_narrays");
    for (n = 1; n <= na; n++) {
        if (cg_array_info (n, name, &dtype, &ndim, dims))
            error_exit ("cg_array_info");
        if (0 == strcmp (name, "InterpolantsDonor")) return 1;
    }
    return 0;
}

/*=======================================================================*/

static char *temporary_file (char *basename)
{
    char *p, *temp;
    int n;

    if (basename == NULL || !*basename)
        basename = "cgnstmpfile";
    n = (int)strlen (basename);
    temp = (char *) malloc (n + 10);
    if (temp == NULL)
        fatal_error("malloc failed for temp filename\n");
    sprintf (temp, "%s.tmp", basename);
    p = temp + strlen(temp);
    for (n = 0; n < 1000; n++) {
        sprintf (p, "%3.3d~", n);
        if (access (temp, 0)) return temp;
    }
    fatal_error("failed to create temporary filename\n");
    return NULL;
}

/*-----------------------------------------------------------------------*/

static void copy_file (char *oldfile, char *newfile)
{
    int c;
    FILE *oldfp, *newfp;

    if (NULL == (oldfp = fopen (oldfile, "rb")))
        fatal_error("error opening input file for reading\n");
    if (NULL == (newfp = fopen (newfile, "w+b"))) {
        fclose (oldfp);
        fatal_error("error opening output file for writing\n");
    }
    while (EOF != (c = getc (oldfp)))
        putc (c, newfp);
    fclose (oldfp);
    fclose (newfp);
}

/*-----------------------------------------------------------------------*/

static char *update_version (char *cgnsfile, char *outfile)
{
    char *tempfile;
    float file_version;

    if (verbose) {
        puts ("checking file version");
        fflush (stdout);
    }
    if (cg_open (cgnsfile, CG_MODE_READ, &cgnsfn) ||
        cg_version (cgnsfn, &file_version) ||
        cg_close (cgnsfn))
        cg_error_exit ();
    if (LibraryVersion <= (int)(file_version * 1000.0 + 0.5)) {
        puts ("file version is current");
        return cgnsfile;
    }
    if (verbose) {
        printf ("creating a working copy of %s\n", cgnsfile);
        fflush (stdout);
    }
    tempfile = temporary_file (cgnsfile);
    copy_file (cgnsfile, tempfile);
    if (verbose) {
        printf ("updating version number for %s\n", tempfile);
        fflush (stdout);
    }
    if (cg_open (tempfile, CG_MODE_MODIFY, &cgnsfn) || cg_close (cgnsfn)) {
        unlink (tempfile);
        cg_error_exit ();
    }
    if (NULL == outfile || !*outfile) outfile = cgnsfile;
    if (verbose) {
        printf ("renaming %s -> %s\n", tempfile, outfile);
        fflush (stdout);
    }
    unlink (outfile);
    if (rename (tempfile, outfile))
        fatal_error("rename %s -> %s failed\n", tempfile, outfile);
    free (tempfile);
    return outfile;
}

/*===================================================================*/

static int sort_nodes (const void *nodes1, const void *nodes2)
{
    cgsize_t n1 = *((cgsize_t *)nodes1);
    cgsize_t n2 = *((cgsize_t *)nodes2);
    if (n1 == n2) return 0;
    return (n1 > n2 ? 1 : -1);
}

/*-------------------------------------------------------------------*/

static int compare_faces (void *v1, void *v2)
{
    int i;
    FACE *f1 = (FACE *)v1;
    FACE *f2 = (FACE *)v2;

    if (f1->nnodes != f2->nnodes)
        return (f1->nnodes - f2->nnodes);
    for (i = 0; i < f1->nnodes; i++) {
        if (f1->nodes[i] != f2->nodes[i])
            return (int)(f1->nodes[i] - f2->nodes[i]);
    }
    return (0);
}

/*-------------------------------------------------------------------*/

static size_t hash_face (void *v)
{
    FACE *f = (FACE *)v;
    int n;
    size_t hash = 0;

    for (n = 0; n < f->nnodes; n++)
        hash += (size_t)f->nodes[n];
    return (hash);
}

/*-----------------------------------------------------------------------*/

static cgsize_t find_extnode (ZONE *z, cgsize_t node)
{
    cgsize_t lo = 0, hi = z->nextnodes - 1, mid;

    if (node == z->extnodes[lo]) return node;
    if (node == z->extnodes[hi]) return node;
    while (lo <= hi) {
        mid = (lo + hi) >> 1;
        if (node == z->extnodes[mid]) return node;
        if (node < z->extnodes[mid])
            hi = mid - 1;
        else
            lo = mid + 1;
    }

    return 0;
}

/*-----------------------------------------------------------------------*/

static size_t get_maxnode (void *vface, void *vmaxnode)
{
    FACE *face = (FACE *)vface;
    int n;
    cgsize_t *maxnode = (cgsize_t *)vmaxnode;

    for (n = 0; n < face->nnodes; n++) {
        if (*maxnode < face->nodes[n]) *maxnode = face->nodes[n];
    }
    return 0;
}

/*-----------------------------------------------------------------------*/

static size_t get_extnodes (void *vface, void *vnodes)
{
    FACE *face = (FACE *)vface;

    if (face->e2 == 0) {
        int n;
        cgsize_t *nodes = (cgsize_t *)vnodes;
        for (n = 0; n < face->nnodes; n++)
            nodes[face->nodes[n]-1] = 1;
    }
    return 0;
}

/*-----------------------------------------------------------------------*/

static int valid_face (ZONE *z, cgsize_t elem)
{
    int ns, nn;
    cgsize_t n, ne, *pe;
    CGNS_ENUMT(ElementType_t) type;

    for (ns = 0; ns < z->nsets; ns++) {
        if (z->sets[ns].invalid || z->sets[ns].ns == 0) continue;
        if (elem >= z->sets[ns].is && elem <= z->sets[ns].ie) {
            type = z->sets[ns].type;
            pe = z->sets[ns].elements;
            ne = elem - z->sets[ns].is;
            if (type == CGNS_ENUMV(NGON_n)) {
                for (n = 0; n < ne; n++) {
                    nn = (int)*pe++;
                    pe += nn;
                }
                return (*pe < 3 ? 0 : 1);
            }
            if (type == CGNS_ENUMV(MIXED)) {
                for (n = 0; n < ne; n++) {
                    type = (CGNS_ENUMT(ElementType_t))*pe++;
                    if (type >= CGNS_ENUMV(NGON_n)) {
                        nn = (int)(type - CGNS_ENUMV(NGON_n));
                    }
                    else {
                        if (cg_npe (type, &nn) || nn <= 0) return 0;
                    }
                    pe += nn;
                }
                type = (CGNS_ENUMT(ElementType_t))*pe;
                if (type >= CGNS_ENUMV(NGON_n)+3) return 1;
            }
            if (type >= CGNS_ENUMV(TRI_3) &&
                type <= CGNS_ENUMV(QUAD_9)) return 1;
            return 0;
        }
    }
    return 0;            
}

/*-----------------------------------------------------------------------*/

static cgsize_t *find_element (ZONE *z, cgsize_t elem, int *dim, int *nnodes)
{
    int ns, nn;
    cgsize_t ne, *nodes;
    CGNS_ENUMT(ElementType_t) type;

    for (ns = 0; ns < z->nsets; ns++) {
        if (z->sets[ns].invalid) continue;
        if (elem >= z->sets[ns].is && elem <= z->sets[ns].ie) {
            ne = elem - z->sets[ns].is;
            nodes = z->sets[ns].elements;
            type = z->sets[ns].type;
            if (type == CGNS_ENUMV(NGON_n)) {
                while (ne-- > 0) {
                    nn = (int)*nodes++;
                    nodes += nn;
                }
                *dim = 2;
                *nnodes = (int)*nodes++;
                return nodes;
            }
            if (type == CGNS_ENUMV(NFACE_n)) {
                while (ne-- > 0) {
                    nn = (int)*nodes++;
                    nodes += nn;
                }
                *dim = 3;
                *nnodes = (int)*nodes++;
                return nodes;
            }
            cg_npe (type, &nn);
            if (nn) {
                nodes += (nn * ne);
            }
            else if (type == CGNS_ENUMV(MIXED)) {
                type = (CGNS_ENUMT(ElementType_t))*nodes++;
                while (ne-- > 0) {
                    if (type >= CGNS_ENUMV(NGON_n))
                        nn = (int)(type - CGNS_ENUMV(NGON_n));
                    else {
                        if (cg_npe (type, &nn) || nn <= 0)
                            return NULL;
                    }
                    nodes += nn;
                    type = (CGNS_ENUMT(ElementType_t))*nodes++;
                }
            }
            else {
                return NULL;
            }
            if (type == CGNS_ENUMV(NODE))
                *dim = 0;
            else if (type < CGNS_ENUMV(TRI_3))
                *dim = 1;
            else if (type < CGNS_ENUMV(TETRA_4) || type >= CGNS_ENUMV(NGON_n))
                *dim = 2;
            else {
                *dim = 3;
                if (type == CGNS_ENUMV(HEXA_27)) nn--;
            }
#ifndef USE_MID_NODES
            switch (type) {
                case CGNS_ENUMV(QUAD_9):
                case CGNS_ENUMV(PYRA_14):
                    nn--;
                    break;
                case CGNS_ENUMV(PENTA_18):
                    nn -= 3;
                    break;
                case CGNS_ENUMV(HEXA_27):
                    nn -= 6;
                    break;
            }
#endif            
            *nnodes = nn;
            return nodes;
        }
    }
    return NULL;
}

/*-----------------------------------------------------------------------*/

static FACE *new_face (int nnodes, cgsize_t *nodes)
{
    FACE *f = (FACE *) malloc (sizeof(FACE) + nnodes * sizeof(cgsize_t));
    if (f == NULL)
        fatal_error("malloc failed for a new face\n");
    f->e1 = f->e2 = 0;
    f->f1 = f->f2 = 0;
    f->nnodes = nnodes;
    f->nodes = (cgsize_t *)(f + 1);
    if (nodes != NULL) {
        int n;
        for (n = 0; n < nnodes; n++)
            f->nodes[n] = nodes[n];
        qsort (f->nodes, f->nnodes, sizeof(cgsize_t), sort_nodes);
    }
    return f;
}

/*-------------------------------------------------------------------*/

static int tetra_4[4][4] = {
    {3, 0, 2, 1},
    {3, 0, 1, 3},
    {3, 1, 2, 3},
    {3, 2, 0, 3}
};
static int tetra_10[4][7] = {
    {6, 0, 6, 2, 5, 1, 4},
    {6, 0, 4, 1, 8, 3, 7},
    {6, 1, 5, 2, 9, 3, 8},
    {6, 2, 6, 0, 7, 3, 9}
};
static int pyra_5[5][5] = {
    {4, 0, 3, 2, 1},
    {3, 0, 1, 4, 0},
    {3, 1, 2, 4, 0},
    {3, 2, 3, 4, 0},
    {3, 3, 0, 4, 0}
};
static int pyra_13[5][9] = {
    {8, 0, 8, 3, 7, 2, 6, 1, 5},
    {6, 0, 5, 1, 10, 4, 9, 0, 0},
    {6, 1, 6, 2, 11, 4, 10, 0, 0},
    {6, 2, 7, 3, 12, 4, 11, 0, 0},
    {6, 3, 8, 0, 9, 4, 12, 0, 0}
};
static int pyra_14[5][10] = {
    {9, 0, 8, 3, 7, 2, 6, 1, 5, 13},
    {6, 0, 5, 1, 10, 4, 9, 0, 0, 0},
    {6, 1, 6, 2, 11, 4, 10, 0, 0, 0},
    {6, 2, 7, 3, 12, 4, 11, 0, 0, 0},
    {6, 3, 8, 0, 9, 4, 12, 0, 0, 0}
};
static int penta_6[5][5] = {
    {4, 0, 1, 4, 3},
    {4, 1, 2, 5, 4},
    {4, 2, 0, 3, 5},
    {3, 0, 2, 1, 0},
    {3, 3, 4, 5, 0},
};
static int penta_15[5][9] = {
    {8, 0, 6, 1, 10, 4, 12, 3, 9},
    {8, 1, 7, 2, 11, 5, 13, 4, 10},
    {8, 2, 8, 0, 9, 3, 14, 5, 11},
    {6, 0, 8, 2, 7, 1, 6, 0, 0},
    {6, 3, 12, 4, 13, 5, 14, 0, 0},
};
static int penta_18[5][10] = {
    {9, 0, 6, 1, 10, 4, 12, 3, 9, 15},
    {9, 1, 7, 2, 11, 5, 13, 4, 10, 16},
    {9, 2, 8, 0, 9, 3, 14, 5, 11, 17},
    {6, 0, 8, 2, 7, 1, 6, 0, 0, 0},
    {6, 3, 12, 4, 13, 5, 14, 0, 0, 0},
};
static int hexa_8[6][5] = {
    {4, 0, 3, 2, 1},
    {4, 0, 1, 5, 4},
    {4, 1, 2, 6, 5},
    {4, 2, 3, 7, 6},
    {4, 0, 4, 7, 3},
    {4, 4, 5, 6, 7}
};
static int hexa_20[6][9] = {
    {8, 0, 11, 3, 10, 2, 9, 1, 8},
    {8, 0, 8, 1, 13, 5, 16, 4, 12},
    {8, 1, 9, 2, 14, 6, 17, 5, 13},
    {8, 2, 10, 3, 15, 7, 18, 6, 14},
    {8, 0, 12, 4, 19, 7, 15, 3, 11},
    {8, 4, 16, 5, 17, 6, 18, 7, 19}
};
static int hexa_27[6][10] = {
    {9, 0, 11, 3, 10, 2, 9, 1, 8, 20},
    {9, 0, 8, 1, 13, 5, 16, 4, 12, 21},
    {9, 1, 9, 2, 14, 6, 17, 5, 13, 22},
    {9, 2, 10, 3, 15, 7, 18, 6, 14, 23},
    {9, 0, 12, 4, 19, 7, 15, 3, 11, 24},
    {9, 4, 16, 5, 17, 6, 18, 7, 19, 25}
};

static FACE *element_face (ZONE *z, int fnum, CGNS_ENUMT(ElementType_t) type,
    cgsize_t *nodes)
{
    int n;
    int *nodemap;
    FACE *face;

    if (type == CGNS_ENUMV(NFACE_n)) {
        int dim;
        cgsize_t *face = find_element (z, abs(nodes[fnum]), &dim, &n);
        if (face == NULL || dim != 2)
            fatal_error("find_element returned invalid face\n");
        return new_face (n, face);
    }
    switch (type) {
        case CGNS_ENUMV(TETRA_4):
            nodemap = tetra_4[fnum];
            break;
        case CGNS_ENUMV(TETRA_10):
            nodemap = tetra_10[fnum];
            break;
        case CGNS_ENUMV(PYRA_5):
            nodemap = pyra_5[fnum];
            break;
        case CGNS_ENUMV(PYRA_13):
            nodemap = pyra_13[fnum];
            break;
        case CGNS_ENUMV(PYRA_14):
#ifdef USE_MID_NODES
            nodemap = pyra_14[fnum];
#else
            nodemap = pyra_13[fnum];
#endif            
            break;
        case CGNS_ENUMV(PENTA_6):
            nodemap = penta_6[fnum];
            break;
        case CGNS_ENUMV(PENTA_15):
            nodemap = penta_15[fnum];
            break;
        case CGNS_ENUMV(PENTA_18):
#ifdef USE_MID_NODES
            nodemap = penta_18[fnum];
#else
            nodemap = penta_15[fnum];
#endif            
            break;
        case CGNS_ENUMV(HEXA_8):
            nodemap = hexa_8[fnum];
            break;
        case CGNS_ENUMV(HEXA_20):
            nodemap = hexa_20[fnum];
            break;
        case CGNS_ENUMV(HEXA_27):
#ifdef USE_MID_NODES        
            nodemap = hexa_27[fnum];
#else
            nodemap = hexa_20[fnum];
#endif            
            break;
        default:
            fatal_error("invalid element type in element_face\n");
    }
    face = new_face (*nodemap++, NULL);
    for (n = 0; n < face->nnodes; n++)
        face->nodes[n] = nodes[nodemap[n]];
    qsort (face->nodes, face->nnodes, sizeof(cgsize_t), sort_nodes);
    return face;        
}

/*-----------------------------------------------------------------------*/

static int compare_elemsets (const void *e1, const void *e2)
{
    cgsize_t s1 = ((ELEMSET *)e1)->is;
    cgsize_t s2 = ((ELEMSET *)e2)->is;
    if (s1 == s2) return 0;
    return (s1 > s2 ? 1 : -1);
}

/*-----------------------------------------------------------------------*/

/* FINISH - allow for CellDim = 2. May also want to read and save
   coordinates to check connectivities */

static void read_zone (int nz)
{
    char name[33];
    int i, j, n;
    cgsize_t size[9];
    int ns, nsets, hasparent;
    cgsize_t ne, *pe;
    cgsize_t se, nelem, k;
    int nn, nf, ip, type, ierr;
    cgsize_t *nodes, maxnode;
    ELEMSET *es;
    CGNS_ENUMT(ZoneType_t) zonetype;
    ZONE *z = &Zones[nz++];
    FACE *face, *pf;

    if (cg_zone_read (cgnsfn, cgnsbase, nz, name, size))
        error_exit("cg_zone_read");
    if (cg_zone_type (cgnsfn, cgnsbase, nz, &zonetype))
        error_exit("cg_zone_type");

    printf ("reading zone \"%s\"\n", name);
    fflush (stdout);

    strcpy (z->name, name);
    z->type = zonetype;
    z->idim = 0;
    z->nnodes = 0;
    z->nsets = 0;
    z->nextnodes = 0;
    z->nv = z->ns = z->ne = z->nn = 0;
    z->faces = NULL;

    for (j = 0; j < 3; j++)
        for (i = 0; i < 3; i++)
            z->dims[j][i] = 0;

    if (zonetype == CGNS_ENUMV(Structured)) {
        z->idim = CellDim;
        for (n = 0, j = 0; j < 3; j++) {
            for (i = 0; i < CellDim; i++) {
                z->dims[j][i] = size[n++];
            }
        }
    }
    else if (zonetype == CGNS_ENUMV(Unstructured)) {
        z->idim = 1;
        for (n = 0; n < 3; n++)
            z->dims[n][0] = size[n];
    }
    else {
        return;
    }

    for (z->nnodes = 1, n = 0; n < z->idim; n++)
        z->nnodes *= z->dims[0][n];
    z->maxnode = z->nnodes;

    /* read element sets */

    if (cg_nsections (cgnsfn, cgnsbase, nz, &nsets))
        error_exit ("cg_nsections");
    if (z->type == CGNS_ENUMV(Structured)) {
        if (nsets)
            warning (1, "element sets are not used with Structured grid");
        return;
    }
    if (!nsets) {
        if (size[1]) error ("no element sets found");
        return;
    }
    z->nsets = nsets;
    z->sets = (ELEMSET *) malloc (nsets * sizeof(ELEMSET));
    if (NULL == z->sets)
        fatal_error("malloc failed for element sets\n");

    /* read element sets */

    for (es = z->sets, ns = 1; ns <= nsets; ns++, es++) {
        es->invalid = 0;
        if (cg_section_read (cgnsfn, cgnsbase, nz, ns, es->name,
                &es->type, &es->is, &es->ie, &es->ib, &hasparent))
            error_exit("cg_section_read");
        printf ("  reading element set \"%s\"\n", es->name);
        fflush (stdout);
        nelem = es->ie - es->is + 1;
        if (cg_ElementDataSize (cgnsfn, cgnsbase, nz, ns, &se))
            error_exit ("cg_ElementDataSize");
        if (se == 0) continue;
        es->elements = (cgsize_t *) malloc ((size_t)(se * sizeof(cgsize_t)));
        if (NULL == es->elements)
            fatal_error("malloc failed for elements\n");
        es->parent = NULL;
        if (hasparent) {
            es->parent = (cgsize_t *) malloc ((size_t)(4 * nelem * sizeof(cgsize_t)));
            if (NULL == es->parent)
                fatal_error("malloc failed for elemset parent data\n");
        }
        if (cg_elements_read (cgnsfn, cgnsbase, nz, ns, es->elements,
                es->parent)) error_exit ("cg_elements_read");

        go_absolute ("Zone_t", nz, "Elements_t", ns, NULL);
        ierr = read_rind (es->rind);
        if (ierr) {
            if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_rind_read");
            es->rind[0] = es->rind[1] = 0;
        }
        else {
            if (FileVersion < 2400)
                error ("rind not valid for elements");
        }

        es->nv = es->ns = es->ne = es->nn = 0;
        if (es->type < CGNS_ENUMV(NODE)) {
            es->invalid = -1;
            continue;
        }
        if (es->type == CGNS_ENUMV(MIXED)) {
            ip = ierr = 0;
            for (pe = es->elements, ne = 0; ne < nelem; ne++) {
                type = (int)*pe++;
                if (type < CGNS_ENUMV(NODE) || type == CGNS_ENUMV(MIXED)) {
                  if (type == CGNS_ENUMV(MIXED))
                        es->invalid = -2;
                    else
                        es->invalid = -1;
                    break;
                }
                if (type >= CGNS_ENUMV(NGON_n)) {
                    ip++;
                    nn = type - CGNS_ENUMV(NGON_n);
                    if (nn < 3) {
                        ierr++;
                    }
                    else {
                        (es->ns)++;
                        if (ne >= es->rind[0] && ne < nelem - es->rind[1])
                            (z->ns)++;
                    }
                    pe += nn;
                    continue;
                }
                if (type == CGNS_ENUMV(NODE)) {
                    (es->nn)++;
                    if (ne >= es->rind[0] && ne < nelem - es->rind[1])
                        (z->nn)++;
                }
                else if (type <= CGNS_ENUMV(BAR_3)) {
                    (es->ne)++;
                    if (ne >= es->rind[0] && ne < nelem - es->rind[1])
                        (z->ne)++;
                }
                else if (type <= CGNS_ENUMV(QUAD_9)) {
                    (es->ns)++;
                    if (ne >= es->rind[0] && ne < nelem - es->rind[1])
                        (z->ns)++;
                }
                else {
                    (es->nv)++;
                    if (ne >= es->rind[0] && ne < nelem - es->rind[1])
                        (z->nv)++;
                }
                if (cg_npe (type, &nn) || nn <= 0) {
                    es->invalid = -3;
                    break;
                }
                pe += nn;
            }
            if (ierr)
                warning (2, "%d NGON_n element with < 3 nodes in MIXED", ierr);
            if (ip) {
                if (FileVersion >= 3000)
                    warning (2,
                        "NGON_n in MIXED deprecated - use NGON_n element set");
            }
        }
        else if (es->type == CGNS_ENUMV(NGON_n)) {
            es->ns = nelem;
            z->ns += (nelem - es->rind[0] - es->rind[1]);
        }
        else if (es->type == CGNS_ENUMV(NFACE_n)) {
            es->nv = nelem;
            z->nv += (nelem - es->rind[0] - es->rind[1]);
        }
        else if (es->type == CGNS_ENUMV(NODE)) {
            es->nn = nelem;
            z->nn += (nelem - es->rind[0] - es->rind[1]);
        }
        else if (es->type <= CGNS_ENUMV(BAR_3)) {
            es->ne = nelem;
            z->ne += (nelem - es->rind[0] - es->rind[1]);
        }
        else if (es->type <= CGNS_ENUMV(QUAD_9)) {
            es->ns = nelem;
            z->ns += (nelem - es->rind[0] - es->rind[1]);
        }
        else {
            es->nv = nelem;
            z->nv += (nelem - es->rind[0] - es->rind[1]);
        }
    }

    /* sort element sets by element range */

    if (nsets > 1)
        qsort (z->sets, nsets, sizeof(ELEMSET), compare_elemsets);

    if (!z->nv) return;

    /* check element indices for volume face hash table */

    ierr = 0;
    for (es = z->sets, ns = 0; ns < nsets; ns++, es++) {
        if (es->invalid || es->nv == 0) continue;
        nelem = es->ie - es->is + 1 - es->rind[1];
        pe = es->elements;
        if (es->type == CGNS_ENUMV(MIXED)) {
            for (ne = 0; ne < nelem; ne++) {
                type = (int)*pe++;
                if (type >= CGNS_ENUMV(NGON_n))
                    nn = (int)(type - CGNS_ENUMV(NGON_n));
                else
                    cg_npe (type, &nn);
                if (nn <= 0) {
                    es->invalid = -3;
                    break;
                }
                if (ne >= es->rind[0]) {
                    for (i = 0; i < nn; i++) {
                        if (pe[i] < 1 || pe[i] > z->nnodes) {
                            ierr++;
                            (es->invalid)++;
                        }
                    }
                }
                pe += nn;
            }
        }
        else if (es->type == CGNS_ENUMV(NGON_n)) {
            for (ne = 0; ne < nelem; ne++) {
                nn = (int)*pe++;
                if (ne >= es->rind[0]) {
                    for (i = 0; i < nn; i++) {
                        if (pe[i] < 1 || pe[i] > z->maxnode) {
                            ierr++;
                            (es->invalid)++;
                        }
                    }
                }
                pe += nn;
            }
        }
        else if (es->type == CGNS_ENUMV(NFACE_n)) {
            for (ne = 0; ne < nelem; ne++) {
                nn = (int)*pe++;
                if (ne >= es->rind[0]) {
                    for (i = 0; i < nn; i++) {
                        if (!valid_face (z, abs(pe[i]))) {
                            ierr++;
                            (es->invalid)++;
                        }
                    }
                }
                pe += nn;
            }
        }
        else {
            if (cg_npe (es->type, &nn) || nn <= 0) {
                es->invalid = -3;
                continue;
            }
            if (es->rind[0] > 0) {
                nelem -= es->rind[0];
                pe += (nn * es->rind[0]);
            }
            for (k = 0; k < nn*nelem; k++) {
                if (pe[k] < 1 || pe[k] > z->nnodes) {
                    ierr++;
                    (es->invalid)++;
                }
            }
        }
    }

    /* build face hash table from volume elements */

    if (ierr) {
        printf("  skipping volume face checking since there were %d errors\n",
            ierr);
        return;
    }
    puts ("  building volume faces hash table...");
    fflush (stdout);
    z->faces = HashCreate ((size_t)(z->nv >> 2) + 1, compare_faces, hash_face);
    if (z->faces == NULL)
        fatal_error("malloc failed for face hash table\n");

    ierr = 0;
    for (es = z->sets, ns = 0; ns < nsets; ns++, es++) {
        if (es->invalid || es->nv == 0 || es->type < CGNS_ENUMV(TETRA_4))
            continue;
        if (es->type > CGNS_ENUMV(MIXED)) {
            if (es->type != CGNS_ENUMV(NFACE_n))
            continue;
        }
        nelem = es->ie - es->is + 1 - es->rind[1];
        type = es->type;
        pe = es->elements;
        cg_npe (es->type, &nn);

        for (ne = 0; ne < nelem; ne++) {
            if (es->type == CGNS_ENUMV(MIXED)) {
                type = (int)*pe++;
                if (type >= CGNS_ENUMV(NGON_n))
                    nn = (int)(type - CGNS_ENUMV(NGON_n));
                else
                    cg_npe (type, &nn);
            }
            switch (type) {
                case CGNS_ENUMV(TETRA_4):
                case CGNS_ENUMV(TETRA_10):
                    nf = 4;
                    break;
                case CGNS_ENUMV(PYRA_5):
                case CGNS_ENUMV(PYRA_13):
                case CGNS_ENUMV(PYRA_14):
                case CGNS_ENUMV(PENTA_6):
                case CGNS_ENUMV(PENTA_15):
                case CGNS_ENUMV(PENTA_18):
                    nf = 5;
                    break;
                case CGNS_ENUMV(HEXA_8):
                case CGNS_ENUMV(HEXA_20):
                case CGNS_ENUMV(HEXA_27):
                    nf = 6;
                    break;
                case CGNS_ENUMV(NFACE_n):
                    nf = (int)*pe++;
                    nn = nf;
                    break;
                default:
                    nf = 0;
                    break;
            }
            if (nf > 0 && ne >= es->rind[0]) {
                for (j = 0; j < nf; j++) {
                    face = element_face (z, j, type, pe);
                    pf = (FACE *) HashFind (z->faces, face);
                    if (pf == NULL) {
                        face->e1 = es->is + ne;
                        face->f1 = j + 1;
                        (void) HashAdd (z->faces, face);
                    }
                    else if (pf->e2) {
                        ierr++;
                    }
                    else {
                        pf->e2 = es->is + ne;
                        pf->f2 = j + 1;
                    }
                }
            }
            pe += nn;
        }
    }
    
    if (ierr)
        warning (1, "%d faces are shared by more than 2 volumes", ierr);

    /* get the list of exterior nodes */

    puts ("  finding exterior nodes...");
    fflush (stdout);
    maxnode = 0;
    HashList (z->faces, get_maxnode, &maxnode);
    nodes = (cgsize_t *) calloc ((size_t)maxnode, sizeof(cgsize_t));
    if (nodes == NULL)
        fatal_error("malloc failed for zone nodes\n");
    HashList (z->faces, get_extnodes, nodes);
    for (ne = 0, k = 0; k < z->nnodes; k++) {
        if (nodes[k]) ne++;
    }
    z->nextnodes = ne;
    z->extnodes = (cgsize_t *) malloc ((size_t)(ne * sizeof(cgsize_t)));
    if (z->extnodes == NULL)
        fatal_error("malloc failed for zone exterior nodes\n");
    for (ne = 0, k = 0; k < z->nnodes; k++) {
        if (nodes[k]) z->extnodes[ne++] = k + 1;
    }
    free (nodes);
}

/*=======================================================================*/

static cgsize_t get_data_size (ZONE *z, CGNS_ENUMT(GridLocation_t) location,
    int *rind)
{
    int n, i;
    cgsize_t datasize = 1;

    if (location == CGNS_ENUMV(Vertex)) {
        for (n = 0, i = 0; i < z->idim; i++) {
            datasize *= (z->dims[0][i] + rind[n] + rind[n+1]);
            n += 2;
        }
        return datasize;
    }

    if (location == CGNS_ENUMV(CellCenter)) {
        for (n = 0, i = 0; i < z->idim; i++) {
            datasize *= (z->dims[1][i] + rind[n] + rind[n+1]);
            n += 2;
        }
        return datasize;
    }

    if (z->type == CGNS_ENUMV(Unstructured)) {
        error ("grid location %s not valid for unstructured zone",
            cg_GridLocationName (location));
        return 0;
    }

    if (location == CGNS_ENUMV(FaceCenter)) {
        if (z->idim > 2) {
            error ("location is FaceCenter but index dimension > 2");
            return 0;
        }
        for (n = 0, i = 0; i < z->idim; i++) {
            datasize *= (z->dims[1][i] + rind[n] + rind[n+1]);
            n += 2;
        }
        return datasize;
    }

    if (location == CGNS_ENUMV(EdgeCenter)) {
        if (z->idim > 1) {
            error ("location is EdgeCenter but index dimension > 1");
            return 0;
        }
        datasize = z->dims[1][0] + rind[0] + rind[1];
        return datasize;
    }

    if (location == CGNS_ENUMV(IFaceCenter)) {
        for (n = 0, i = 1; i < z->idim; i++) {
            if (i == 0)
                datasize *= (z->dims[0][i] + rind[n] + rind[n+1]);
            else
                datasize *= (z->dims[1][i] + rind[n] + rind[n+1]);
            n += 2;
        }
        return datasize;
    }

    if (location == CGNS_ENUMV(JFaceCenter)) {
        if (z->idim < 2) {
            error ("location is JFaceCenter but index dimension < 2");
            return 0;
        }
        for (n = 0, i = 1; i < z->idim; i++) {
            if (i == 1)
                datasize *= (z->dims[0][i] + rind[n] + rind[n+1]);
            else
                datasize *= (z->dims[1][i] + rind[n] + rind[n+1]);
            n += 2;
        }
        return datasize;
    }

    if (location == CGNS_ENUMV(KFaceCenter)) {
        if (z->idim < 3) {
            error ("location is KFaceCenter but index dimension < 3");
            return 0;
        }
        for (n = 0, i = 1; i < z->idim; i++) {
            if (i == 2)
                datasize *= (z->dims[0][i] + rind[n] + rind[n+1]);
            else
                datasize *= (z->dims[1][i] + rind[n] + rind[n+1]);
            n += 2;
        }
        return datasize;
    }

    error ("grid location %s is invalid", cg_GridLocationName(location));
    return 0;
}

/*=======================================================================*/

static int read_dataclass (void)
{
    int ierr;
    CGNS_ENUMT(DataClass_t) dataclass;

    ierr = cg_dataclass_read (&dataclass);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_dataclass_read");
        return -1;
    }
    return (int)dataclass;
}

/*-----------------------------------------------------------------------*/

static int *read_units (int units[9])
{
    int n, ierr;
    CGNS_ENUMT(MassUnits_t) mass;
    CGNS_ENUMT(LengthUnits_t) length;
    CGNS_ENUMT(TimeUnits_t) time;
    CGNS_ENUMT(TemperatureUnits_t) temp;
    CGNS_ENUMT(AngleUnits_t) angle;
    CGNS_ENUMT(ElectricCurrentUnits_t) current;
    CGNS_ENUMT(SubstanceAmountUnits_t) amount;
    CGNS_ENUMT(LuminousIntensityUnits_t) intensity;

    for (n = 0; n < 9; n++)
        units[n] = 0;
    ierr = cg_unitsfull_read (&mass, &length, &time, &temp, &angle,
                              &current, &amount, &intensity);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_unitsfull_read");
        return NULL;
    }
    cg_nunits (&n);
    units[0] = mass;
    units[1] = length;
    units[2] = time;
    units[3] = temp;
    units[4] = angle;
    units[5] = current;
    units[6] = amount;
    units[7] = intensity;
    units[8] = n;
    return units;
}

/*-----------------------------------------------------------------------*/

static int read_exponents (float exps[9])
{
    int n, ierr;
    CGNS_ENUMT(DataType_t) type;

    for (n = 0; n < 9; n++)
        exps[n] = 0;
    ierr = cg_exponents_info (&type);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_exponents_info");
        return 0;
    }
    if (type == CGNS_ENUMV(RealSingle))
        ierr = cg_expfull_read (exps);
    else if (type == CGNS_ENUMV(RealDouble)) {
        double data[8];
        ierr = cg_expfull_read (data);
        if (ierr == CG_OK) {
            for (n = 0; n < 8; n++)
                exps[n] = (float)data[n];
        }
    }
    else {
        error("invalid data type for exponents");
        return 0;
    }
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND)
            error_exit("cg_expfull_read");
        return 0;
    }
    cg_nexponents (&n);
    exps[8] = (float)n;
    return 1;
}

/*-----------------------------------------------------------------------*/

static void print_indent (int indent)
{
    while (indent-- > 0) putchar (' ');
}

/*-----------------------------------------------------------------------*/

static void print_dataclass (int dataclass, int indent)
{
    print_indent (indent);
    printf ("Data Class=");
    if (dataclass < 0)
        puts ("<not specified>");
    else
        puts (cg_DataClassName(dataclass));
}

/*-----------------------------------------------------------------------*/

static void print_units (int *units, int indent)
{
    print_indent (indent);
    printf ("Units=[%s,%s,%s,%s,%s",
        cg_MassUnitsName(units[0]),
        cg_LengthUnitsName(units[1]),
        cg_TimeUnitsName(units[2]),
        cg_TemperatureUnitsName(units[4]),
        cg_AngleUnitsName(units[4]));
    if (units[8] > 5)
        printf (",%s,%s,%s",
            cg_ElectricCurrentUnitsName(units[5]),
            cg_SubstanceAmountUnitsName(units[6]),
            cg_LuminousIntensityUnitsName(units[7]));
    puts ("]");
}

/*-----------------------------------------------------------------------*/

static void print_exponents (float *exps, int indent)
{
    print_indent (indent);
    printf ("Exponents=[%g,%g,%g,%g,%g", exps[0], exps[1],
        exps[2], exps[3], exps[4]);
    if (exps[8] > 5.0)
        printf (",%g,%g,%g", exps[5], exps[6], exps[7]);
    puts ("]");
}

/*=======================================================================*/

static void check_quantity (int dnum, char *name,
    int parclass, int *parunits, int isref, int indent)
{
    int n, ne, dclass, *punits, hasexps, units[9];
    float defexps[8], exps[9];

    go_relative ("DataArray_t", dnum, NULL);
    dclass = read_dataclass ();
    punits = read_units (units);
    hasexps = read_exponents (exps);
    go_relative ("..", 1, NULL);

    if (verbose) {
        if (dclass >= 0) print_dataclass (dclass, indent);
        if (punits != NULL) print_units (punits, indent);
        if (hasexps) print_exponents (exps, indent);
    }
    if (dclass < 0) dclass = parclass;

    if (cg_get_identifier (name, &ne, defexps)) {
        if (isref > 0)
            warning (3, "not a CGNS data-name identifier");
        if (dclass < 0)
            warning (3, "dataclass is not given");
        else if ((CGNS_ENUMT(DataClass_t))dclass == CGNS_ENUMV(Dimensional)) {
            if (punits == NULL && parunits == NULL)
                warning (2, "units not given");
            if (!hasexps)
                warning (2, "exponents not given");
        }
        else {
            if (punits != NULL || hasexps)
                warning (2, "dataclass is %s, but units and/or exponents"
                    " are given", cg_DataClassName(dclass));
        }
        return;
    }

    if (isref < 0) return;

    if (dclass < 0)
        warning (2, "dataclass not given");
    else if ((CGNS_ENUMT(DataClass_t))dclass == CGNS_ENUMV(Dimensional) ||
             (CGNS_ENUMT(DataClass_t))dclass == CGNS_ENUMV(NormalizedByDimensional) ||
             (CGNS_ENUMT(DataClass_t))dclass == CGNS_ENUMV(NormalizedByUnknownDimensional)) {
        if (!ne)
            warning (2, "dataclass does not match CGNS specification");
        if ((CGNS_ENUMT(DataClass_t))dclass == CGNS_ENUMV(Dimensional) &&
            punits == NULL && parunits == NULL)
            warning (2, "units not given for dimensional quantity");
        if (ne > 0) {
            if (hasexps) {
                for (n = 0; n < ne; n++) {
                    if (exps[n] != defexps[n]) break;
                }
                if (n < ne)
                    warning (2, "exponents do not match CGNS specification");
            }
            else
                warning (3, "exponents not given");
        }
    }
    else if ((CGNS_ENUMT(DataClass_t))dclass == CGNS_ENUMV(NondimensionalParameter) ||
             (CGNS_ENUMT(DataClass_t))dclass == CGNS_ENUMV(DimensionlessConstant)) {
        if (ne)
            warning (2, "dataclass does not match CGNS specification");
        if (punits != NULL)
            warning (2, "units given for nondimensional quantity");
        if (hasexps)
            warning (2, "exponents given for nondimensional quantity");
    }
    else
        error ("invalid dataclass");
}

/*-----------------------------------------------------------------------*/

static void check_arrays (int parclass, int *parunits, int isref,
    cgsize_t length, int indent)
{
    int n, narrays, na, ndim;
    cgsize_t size, dims[12];
    int dataclass, *punits, units[9];
    CGNS_ENUMT(DataType_t) datatype;
    char name[33];

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, indent);
        if (punits) print_units (punits, indent);
    }
    if (dataclass < 0) dataclass = parclass;
    if (!punits) punits = parunits;

    if (cg_narrays (&narrays)) error_exit ("cg_narrays");

    for (na = 1; na <= narrays; na++) {
        if (cg_array_info (na, name, &datatype, &ndim, dims))
            error_exit ("cg_array_info");
        print_indent (indent);
        printf ("checking quantity \"%s\"\n", name);
        fflush (stdout);
        for (size = 1, n = 0; n < ndim; n++)
            size *= dims[n];
        if (ndim < 1 || size < 1)
            error ("invalid dimensions");
        if (length < 0 && size != 1 && size != -length)
            error ("array size not 1 or %d", -length);
        if (length > 0 && size != length)
            error ("array size not %d", length);
        check_quantity (na, name, dataclass, punits, isref, indent+2);
    }
}

/*-----------------------------------------------------------------------*/

static void check_user_data (int parclass, int *parunits, int indent)
{
    int n, nd, nu, nuser, na, ndim;
    cgsize_t dims[12];
    int dataclass, *punits, units[9];
    char name[33], *desc;
    CGNS_ENUMT(DataType_t) datatype;
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(PointSetType_t) ptype;
    int hasf, haso, hasl, hasp, ordinal;
    cgsize_t npnts;

    if (cg_nuser_data (&nuser)) error_exit ("cg_nuser_data");
    if (nuser <= 0) return;

    for (nu = 1; nu <= nuser; nu++) {
        if (cg_user_data_read (nu, name)) error_exit("cg_user_data_read");
        print_indent (indent);
        printf ("checking user data \"%s\"\n", name);
        fflush (stdout);
        go_relative ("UserDefinedData_t", nu, NULL);
        hasf = cg_famname_read (name);
        if (hasf && hasf != CG_NODE_NOT_FOUND)
            error_exit("cg_famname_read");
        haso = read_ordinal (&ordinal);
        if (haso && haso != CG_NODE_NOT_FOUND)
            error_exit("cg_ordinal_read");
        hasl = read_gridlocation (&location);
        if (hasl && hasl != CG_NODE_NOT_FOUND)
            error_exit("cg_gridlocation_read");
        /* ptset is only allowed below a zone */
        if (cgnszone == 0)
            hasp = CG_NODE_NOT_FOUND;
        else {
            hasp = cg_ptset_info (&ptype, &npnts);
            if (hasp && hasp != CG_NODE_NOT_FOUND)
                error_exit("cg_ptset_info");
        }
        if (verbose) {
            if (hasf == CG_OK) {
                print_indent (indent+2);
                printf ("Family Name=\"%s\"\n", name);
            }
            if (haso == CG_OK) {
                print_indent (indent+2);
                printf ("Ordinal=%d\n", ordinal);
            }
            if (hasl == CG_OK) {
                print_indent (indent+2);
                printf ("Grid Location=%s\n", cg_GridLocationName(location));
            }
            if (hasp == CG_OK) {
                print_indent (indent+2);
                printf ("Point Set Type=%s\n", cg_PointSetTypeName(ptype));
                print_indent (indent+2);
                printf ("Number Points=%ld\n", (long)npnts);
            }
        }
        if (hasf == CG_OK) {
            for (n = 0; n < NumFamily; n++) {
                if (0 == strcmp (name, Family[n])) break;
            }
            if (n == NumFamily)
                warning (2, "family name \"%s\" not found", name);
        }
        if (verbose > 1) {
            if (cg_ndescriptors (&nd)) error_exit ("cg_ndescriptors");
            for (n = 1; n <= nd; n++) {
                if (cg_descriptor_read (n, name, &desc))
                    error_exit("cg_descriptor_read");
                if (desc != NULL) {
                    print_indent (indent + 2);
                    printf ("Descriptor %s:\n%s\n", name, desc);
                    cg_free (desc);
                }
            }
        }

        dataclass = read_dataclass ();
        punits = read_units (units);
        if (verbose) {
            if (dataclass >= 0) print_dataclass (dataclass, indent + 2);
            if (punits) print_units (punits, indent + 2);
        }
        if (dataclass < 0) dataclass = parclass;
        if (!punits) punits = parunits;

        if (cg_narrays (&na)) error_exit ("cg_narrays");
        for (n = 1; n <= na; n++) {
            if (cg_array_info (n, name, &datatype, &ndim, dims))
                error_exit ("cg_array_info");
            print_indent (indent + 2);
            printf ("checking quantity \"%s\"\n", name);
            fflush (stdout);
            check_quantity (n, name, dataclass, punits, 0, indent+4);
        }

        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }
}

/*-----------------------------------------------------------------------*/

static void check_integral (int parclass, int *parunits, int indent)
{
    char *desc, name[33];
    int n, ni, nint, na, nd, ndim;
    cgsize_t dims[12];
    int dataclass, *punits, units[9];
    CGNS_ENUMT(DataType_t) datatype;

    if (cg_nintegrals (&nint)) error_exit ("cg_nintegrals");
    if (nint <= 0) return;
    print_indent (indent);
    puts ("checking integral data");
    fflush (stdout);

    for (ni = 1; ni <= nint; ni++) {
        if (cg_integral_read (ni, name)) error_exit("cg_integral_read");
        go_relative ("IntegralData_t", ni, NULL);
        if (verbose) {
            print_indent (indent + 2);
            printf ("integral data \"%s\"\n", name);
            if (verbose > 1) {
                if (cg_ndescriptors (&nd)) error_exit ("cg_ndescriptors");
                for (n = 1; n <= nd; n++) {
                    if (cg_descriptor_read (n, name, &desc))
                        error_exit("cg_descriptor_read");
                    if (desc != NULL) {
                        print_indent (indent + 4);
                        printf ("Descriptor %s:\n%s\n", name, desc);
                        cg_free (desc);
                    }
                }
            }
        }

        dataclass = read_dataclass ();
        punits = read_units (units);
        if (verbose) {
            if (dataclass >= 0) print_dataclass (dataclass, indent + 4);
            if (punits) print_units (punits, indent + 4);
        }
        if (dataclass < 0) dataclass = parclass;
        if (!punits) punits = parunits;

        if (cg_narrays (&na)) error_exit ("cg_narrays");
        for (n = 1; n <= na; n++) {
            if (cg_array_info (n, name, &datatype, &ndim, dims))
                error_exit ("cg_array_info");
            print_indent (indent + 4);
            printf ("checking quantity \"%s\"\n", name);
            fflush (stdout);
            if (ndim != 1 || dims[0] != 1)
                error ("dimension is not 1");
            check_quantity (n, name, dataclass, punits, 0, indent+6);
        }

        check_user_data (dataclass, punits, indent + 4);
        go_relative ("..", 1, NULL);
    }
}

/*-----------------------------------------------------------------------*/

static void check_rotating (float *point, float *vector,
    int parclass, int *parunits, int indent)
{
    char *desc, name[33];
    int n, na, nd, ndim;
    cgsize_t dims[12];
    int dataclass, *punits, units[9];
    CGNS_ENUMT(DataType_t) datatype;

    go_relative ("RotatingCoordinates_t", 1, NULL);

    if (verbose) {
        print_indent (indent);
        printf ("Center=[%g", point[0]);
        for (nd = 1; nd < PhyDim; nd++)
            printf (",%g", point[nd]);
        puts ("]");
        print_indent (indent);
        printf ("Rate Vector=[%g", vector[0]);
        for (nd = 1; nd < PhyDim; nd++)
            printf (",%g", vector[nd]);
        puts ("]");
    }

    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                print_indent (indent);
                printf ("Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, indent);
        if (punits) print_units (punits, indent);
    }
    if (dataclass < 0) dataclass = parclass;
    if (!punits) punits = parunits;

    if (cg_narrays (&na)) error_exit("cg_narrays");
    for (n = 1; n <= na; n++) {
        if (cg_array_info (n, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        print_indent (indent);
        printf ("checking rotating data \"%s\"\n", name);
        fflush (stdout);
        if (strcmp (name, "RotationCenter") &&
            strcmp (name, "RotationRateVector"))
            warning (1, "not valid as child of RotatingCoordinates");
        else
            check_quantity (n, name, dataclass, punits, 1, indent + 4);
    }

    check_user_data (dataclass, punits, indent + 2);

    go_relative ("..", 1, NULL);
}

/*-----------------------------------------------------------------------*/

static void check_convergence (int niter, char *NormDefs,
    int parclass, int *parunits, int indent)
{
    char name[33];
    int n, na, ndim;
    cgsize_t dims[12];
    int dataclass, *punits, units[9];
    CGNS_ENUMT(DataType_t) datatype;

    if (verbose) {
        print_indent (indent);
        printf ("Number Iterations=%d\n", niter);
        if (NormDefs != NULL) {
            print_indent (indent);
            printf ("Norm Definitions=%s\n", NormDefs);
        }
    }
    if (niter < 1)
        warning (2, "number of iterations is not > 0");

    if (verbose > 1) {
        int nd;
        char *desc;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL && strcmp (name, "NormDefinitions")) {
                print_indent (indent);
                printf ("Descriptor %s:\n%s\n", name, desc);
            }
            if (desc != NULL) cg_free (desc);
        }
    }

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, indent);
        if (punits) print_units (punits, indent);
    }
    if (dataclass < 0) dataclass = parclass;
    if (!punits) punits = parunits;

    if (cg_narrays (&na)) error_exit("cg_narrays");
    for (n = 1; n <= na; n++) {
        if (cg_array_info (n, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        print_indent (indent);
        printf ("checking convergence data \"%s\"\n", name);
        fflush (stdout);
        if (ndim != 1 || dims[0] != niter)
            error ("length of array is not the number of iterations");
        if (0 == strncmp (name, "RSD", 3) ||
            0 == strncmp (name, "CHG", 3))
            check_quantity (n, NULL, dataclass, punits, -1, indent + 2);
        else
            check_quantity (n, name, dataclass, punits, 1, indent + 2);
    }

    check_user_data (dataclass, punits, indent);
}

/*-----------------------------------------------------------------------*/

static void check_equation_set (int *flags, int parclass, int *parunits,
    int indent)
{
    char *desc, name[33];
    int n, nd, dataclass, ierr, *punits, units[9], ndiff, diff[6];
    CGNS_ENUMT(GoverningEquationsType_t) governing;
    CGNS_ENUMT(ModelType_t) model;
    int thermrelax, chemkin;
    int emelec, emmagn, emcond;

    ierr = cg_equationset_chemistry_read (&thermrelax, &chemkin);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND)
            error_exit("cg_equationset_chemistry_read");
        thermrelax = chemkin = 0;
    }
    ierr = cg_equationset_elecmagn_read (&emelec, &emmagn, &emcond);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND);
            error_exit("cg_equationset_elecmagn_read");
        emelec = emmagn = emcond = 0;
    }

    go_relative ("FlowEquationSet_t", 1, NULL);

    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                print_indent (indent);
                printf ("Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    if (LibraryVersion < 2000) {
        dataclass = -1;
        punits = NULL;
    }
    else {
        dataclass = read_dataclass ();
        punits = read_units (units);
    }
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, indent);
        if (punits) print_units (punits, indent);
        print_indent (indent);
        printf ("Equation Dimension=%d\n", flags[0]);
    }
    if (dataclass < 0) dataclass = parclass;
    if (!punits) punits = parunits;
    ndiff = (CellDim * (CellDim + 1)) / 2;

    if (flags[0] < 1)
        error ("equation dimension < 1");

    /* governing equations */

    if (flags[1]) {
        if (cg_governing_read (&governing))
            error_exit("cg_governing_read");
        go_relative ("GoverningEquations_t", 1, NULL);
        print_indent (indent);
        printf ("Governing Equation=%s\n",
            cg_GoverningEquationsTypeName(governing));
        if (LibraryVersion < 2200 && goDepth == 2) {
            warning (3, "can't get Diffusion Model for FlowEquations"
                " under CGNSBase_t\n        due to a bug in the CGNS"
                " library for versions prior to 2.2");
        }
        else {
            ierr = cg_diffusion_read (diff);
            if (ierr && ierr != CG_NODE_NOT_FOUND)
                error_exit("cg_diffusion_read");
            if (ierr == CG_OK) {
                if (verbose) {
                    print_indent (indent + 2);
                    printf ("Diffusion Model=[%d", diff[0]);
                    for (n = 1; n < ndiff; n++)
                        printf (",%d", diff[n]);
                    puts ("]");
                }
                for (n = 0; n < ndiff; n++) {
                    if (diff[n] < 0 || diff[n] > 1) break;
                }
                if (n < ndiff)
                    warning (1, "diffusion model terms not 0 or 1");
            }
        }
        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }

    /* gas model */

    if (flags[2]) {
        if (cg_model_read ("GasModel_t", &model))
            error_exit("cg_model_read");
        print_indent (indent);
        printf ("Gas Model=%s\n", cg_ModelTypeName(model));
        go_relative ("GasModel_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, indent + 2);
        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }

    /* viscosity model */

    if (flags[3]) {
        if (cg_model_read ("ViscosityModel_t", &model))
            error_exit("cg_model_read");
        print_indent (indent);
        printf ("Viscosity Model=%s\n", cg_ModelTypeName(model));
        go_relative ("ViscosityModel_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, indent + 2);
        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }

    /* thermal conductivity model */

    if (flags[4]) {
        if (cg_model_read ("ThermalConductivityModel_t", &model))
            error_exit("cg_model_read");
        print_indent (indent);
        printf ("Thermal Conductivity Model=%s\n",
            cg_ModelTypeName(model));
        go_relative ("ThermalConductivityModel_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, indent + 2);
        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }

    /* turbulence closure model */

    if (flags[5]) {
        if (cg_model_read ("TurbulenceClosure_t", &model))
            error_exit("cg_model_read");
        print_indent (indent);
        printf ("Turbulence Closure=%s\n", cg_ModelTypeName(model));
        go_relative ("TurbulenceClosure_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, indent + 2);
        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }

    /* turbulence model */

    if (flags[6]) {
        if (cg_model_read ("TurbulenceModel_t", &model))
            error_exit("cg_model_read");
        go_relative ("TurbulenceModel_t", 1, NULL);
        print_indent (indent);
        printf ("Turbulence Model=%s\n", cg_ModelTypeName(model));
        if (LibraryVersion < 2200 && goDepth == 2) {
            warning (3, "can't get Turbulent Diffusion Model for FlowEquations"
                " under CGNSBase_t\n        due to a bug in the CGNS library"
                " for versions prior to 2.2");
        }
        else {
            ierr = cg_diffusion_read (diff);
            if (ierr && ierr != CG_NODE_NOT_FOUND)
                error_exit("cg_diffusion_read");
            if (ierr == CG_OK) {
                if (verbose) {
                    print_indent (indent + 2);
                    printf ("Turbulent Diffusion Model=[%d", diff[0]);
                    for (n = 1; n < ndiff; n++)
                        printf (",%d", diff[n]);
                    puts ("]");
                }
                for (n = 0; n < ndiff; n++) {
                    if (diff[n] < 0 || diff[n] > 1) break;
                }
                if (n < ndiff)
                    warning (1, "turbulent diffusion model terms not 0 or 1");
            }
        }
        check_arrays (dataclass, punits, 1, 0, indent + 2);
        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }

    /* thermal relaxation model */

    if (thermrelax) {
        if (cg_model_read ("ThermalRelaxationModel_t", &model))
            error_exit ("cg_model_read");
        print_indent (indent);
        printf ("Thermal Relaxation Model=%s\n", cg_ModelTypeName(model));
        go_relative ("ThermalRelaxationModel_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, indent + 2);
        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }

    /* chemical kinetics model */


    if (chemkin) {
        if (cg_model_read ("ChemicalKineticsModel_t", &model))
            error_exit ("cg_model_read");
        print_indent (indent);
        printf ("Chemical Kinetics Model=%s\n", cg_ModelTypeName(model));
        go_relative ("ChemicalKineticsModel_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, indent + 2);
        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }

    /* EM electric field model */

    if (emelec) {
        if (cg_model_read ("EMElectricFieldModel_t", &model))
            error_exit ("cg_model_read");
        print_indent (indent);
        printf ("EM Electric Field Model=%s\n", cg_ModelTypeName(model));
        go_relative ("EMElectricFieldModel_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, indent + 2);
        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }

    /* EM magnetic field model */

    if (emmagn) {
        if (cg_model_read ("EMMagneticFieldModel_t", &model))
            error_exit ("cg_model_read");
        print_indent (indent);
        printf ("EM Magnetic Field Model=%s\n", cg_ModelTypeName(model));
        go_relative ("EMMagneticFieldModel_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, indent + 2);
        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }

    /* EM electric field model */

    if (emcond) {
        if (cg_model_read ("EMConductivityModel_t", &model))
            error_exit ("cg_model_read");
        print_indent (indent);
        printf ("EM Conductivity Model=%s\n", cg_ModelTypeName(model));
        go_relative ("EMConductivityModel_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, indent + 2);
        check_user_data (dataclass, punits, indent + 2);
        go_relative ("..", 1, NULL);
    }

    go_relative ("..", 1, NULL);
}

/*-----------------------------------------------------------------------*/

static void check_coordinates (int ng)
{
    char name[33];
    int ierr, n, rind[6];
    cgsize_t np, rmin[3], rmax[3];
    int nc, ncoords, mask, coordset[4];
    int *punits, units[9], dataclass;
    float *coord, cmin, cmax;
    CGNS_ENUMT(DataType_t) datatype;
    ZONE *z = &Zones[cgnszone-1];

    if (cg_grid_read (cgnsfn, cgnsbase, cgnszone, ng, name))
        error_exit("cg_grid_read");
    strcpy (GridCoordinate[ng-1], name);
    printf ("  checking coordinates \"%s\"\n", name);
    fflush (stdout);

    go_absolute ("Zone_t", cgnszone, "GridCoordinates_t", ng, NULL);

    if (verbose > 1) {
        int nd;
        char *desc;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 4);
        if (punits) print_units (punits, 4);
    }
    if (dataclass < 0) dataclass = z->dataclass;
    if (!punits) punits = z->punits;

    ierr = read_rind (rind);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_rind_read");
        for (n = 0; n < 6; n++)
            rind[n] = 0;
    }
    else {
        if (verbose) {
            printf ("    Rind=[%d", rind[0]);
            for (n = 1; n < 2 * z->idim; n++)
                printf (",%d", rind[n]);
            puts ("]");
        }
        if (z->type == CGNS_ENUMV(Unstructured) && FileVersion < 2400)
            error ("rind not valid for unstructured zones");
    }

    for (n = 0, np = 1; n < z->idim; n++) {
        rmin[n] = 1;
        rmax[n] = z->dims[0][n] + rind[2*n] + rind[2*n+1];
        np *= rmax[n];
    }
    if (NULL == (coord = (float *) malloc ((size_t)(np * sizeof(float)))))
        fatal_error("malloc failed for %ld coordinate values\n", (long)np);
    if (z->maxnode < np) z->maxnode = np;

    if (cg_ncoords (cgnsfn, cgnsbase, cgnszone, &ncoords))
        error_exit("cg_ncoords");
    if (ncoords < PhyDim)
        error ("number coordinates < physical dimensions");
    for (n = 0; n < 4; n++)
        coordset[n] = 0;

    for (nc = 1; nc <= ncoords; nc++) {
        if (cg_coord_info (cgnsfn, cgnsbase, cgnszone, nc, &datatype, name))
            error_exit("cg_coord_info");
        if (cg_coord_read (cgnsfn, cgnsbase, cgnszone, name, CGNS_ENUMV(RealSingle),
                rmin, rmax, coord))
            error_exit("cg_coord_read");
        printf ("    checking coordinate \"%s\"\n", name);
        fflush (stdout);
        cmin = cmax = coord[0];
        for (n = 1; n < np; n++) {
            if (cmin > coord[n]) cmin = coord[n];
            if (cmax < coord[n]) cmax = coord[n];
        }
        if (verbose)
            printf("      Coordinate Range=%g -> %g (%g)\n",
                cmin, cmax, cmax-cmin);
        if (cmin == cmax)
            warning(1, "coordinate range is 0");
        if (0 == strcmp (name, "CoordinateX"))
            coordset[0] |= 1;
        else if (0 == strcmp (name, "CoordinateY"))
            coordset[0] |= 2;
        else if (0 == strcmp (name, "CoordinateZ")) {
            coordset[0] |= 4;
            coordset[1] |= 4;
        }
        else if (0 == strcmp (name, "CoordinateR"))
            coordset[1] |= 1;
        else if (0 == strcmp (name, "CoordinateTheta")) {
            coordset[1] |= 2;
            coordset[2] |= 2;
        }
        else if (0 == strcmp (name, "CoordinatePhi"))
            coordset[2] |= 4;
        else if (0 == strcmp (name, "CoordinateXi"))
            coordset[3] |= 1;
        else if (0 == strcmp (name, "CoordinateEta"))
            coordset[3] |= 2;
        else if (0 == strcmp (name, "CoordinateZeta"))
            coordset[3] |= 4;
        check_quantity (nc, name, dataclass, punits, 1, 6);
    }
    free (coord);

    for (mask = 0, n = 0; n < PhyDim; n++)
        mask |= (1 << n);
    for (n = 0; n < 4; n++) {
        if ((coordset[n] & mask) == mask) break;
    }
    if (n == 4)
        error ("a complete coordinate system was not found");
}

/*-----------------------------------------------------------------------*/

static void check_elements (void)
{
    int nn, ns;
    int type, nf, np, nint, next;
    cgsize_t is, ne, nelem, *pe;
    ELEMSET *es;
    FACE *face, *pf;
    ZONE *z = &Zones[cgnszone-1];

    puts ("  checking elements");
    if (verbose) {
        printf ("    Number Element Sets=%d\n", z->nsets);
        printf ("    [0D,1D,2D,3D] Elements=[%ld,%ld,%ld,%ld]\n",
            (long)z->nn, (long)z->ne, (long)z->ns, (long)z->nv);
        if (z->faces)
            printf ("    Number Volume Faces=%ld\n", (long)HashSize(z->faces));
    }
    fflush (stdout);

    nelem = CellDim == 2 ? z->ns : z->nv;
    if (nelem < z->dims[1][0])
        error ("number %dD elements < specified in zone dimensions", CellDim);
    if (nelem > z->dims[1][0])
        warning (1, "number %dD elements > specified in zone dimensions",
            CellDim);
    if (CellDim == 3 && z->nv == 0)
        error ("no 3D elements exist and CellDim is 3");
    if (CellDim == 2) {
        if (z->nv)
            error ("3D elements exist and CellDim is 2");
        if (z->ns == 0)
            error ("no 2D elements exist and CellDim is 2");
    }

    is = z->sets->is;
    if (is != 1)
        warning (1, "element numbering does not start at 1");

    /* check element sets */

    for (es = z->sets, ns = 0; ns < z->nsets; ns++, es++) {
        printf ("  checking element set \"%s\"\n", es->name);
        if (verbose) {
            printf ("    Element Set Type=%s\n",
                cg_ElementTypeName(es->type));
            printf ("    Element Range=[%ld,%ld]\n", (long)es->is, (long)es->ie);
            if (es->rind[0] || es->rind[1])
                printf ("    Rind Elements=[%d,%d]\n",
                    es->rind[0], es->rind[1]);
            printf ("    [0D,1D,2D,3D] Elements=[%ld,%ld,%ld,%ld]\n",
                (long)es->nn, (long)es->ne, (long)es->ns, (long)es->nv);
        }
        if (ns && z->sets[ns].is != is)
            warning (1, "element numbers are not consecutative with \"%s\"",
                z->sets[ns-1].name);
        is = z->sets[ns].ie + 1;
        for (nn = 0; nn < z->nsets; nn++) {
            if (nn != ns) {
                if ((z->sets[nn].is >= z->sets[ns].is &&
                     z->sets[nn].is <= z->sets[ns].ie) ||
                    (z->sets[nn].ie >= z->sets[ns].is &&
                     z->sets[nn].ie <= z->sets[ns].ie))
                    error ("element numbers overlap those for \"%s\"",
                        z->sets[nn].name);
            }
        }

        if (es->invalid < 0) {
            if (es->invalid == -1)
                printf ("INTERNAL:can't handle %s element set\n",
                    cg_ElementTypeName(es->type));
            else if (es->invalid == -2)
                error ("MIXED elements not allowed in a MIXED element set");
            else if (es->invalid == -3)
                printf("INTERNAL: error with element size in MIXED\n");
            else
                printf("INTERNAL: error processing element set\n");
            continue;
        }
        if (es->invalid > 0) {
            if (es->type == CGNS_ENUMV(NFACE_n)) {
                error ("%d polyhedra faces are not valid face elements",
                    es->invalid);
                continue;
            }
            error ("%d element nodes are out of range", es->invalid);
            continue;
        }
        if (z->faces == NULL || es->ns == 0 || es->type < CGNS_ENUMV(TRI_3))
            continue;
        if (es->type > CGNS_ENUMV(QUAD_9) &&
            es->type != CGNS_ENUMV(MIXED) &&
            es->type != CGNS_ENUMV(NGON_n)) continue;

        nelem = es->ie - es->is + 1 - es->rind[1];
        type = es->type;
        pe = es->elements;
        nf = np = nint = next = 0;
        for (ne = 0; ne < nelem; ne++) {
            if (es->type == CGNS_ENUMV(MIXED)) {
                type = (int)*pe++;
                if (type >= CGNS_ENUMV(NGON_n)) {
                    nn = (int)(type - CGNS_ENUMV(NGON_n));
                    type = CGNS_ENUMV(NGON_n);
                }
                else {
                    cg_npe (type, &nn);
                }
            }
            else if (es->type == CGNS_ENUMV(NGON_n)) {
                nn = (int)*pe++;
            }
            else {
                cg_npe (type, &nn);
            }
            if (ne >= es->rind[0] && nn >= 3 && (type == CGNS_ENUMV(NGON_n) ||
                (type >= CGNS_ENUMV(TRI_3) && type <= CGNS_ENUMV(QUAD_9)))) {
#ifndef USE_MID_NODES
                if (type == CGNS_ENUMV(QUAD_9))
                    face = new_face (8, pe);
                else
#endif                               
                face = new_face (nn, pe);
                pf = (FACE *) HashFind (z->faces, face);
                free (face);
                if (pf == NULL)
                    nf++;
                else {
                    if (pf->e2)
                        nint++;
                    else
                        next++;
                    if (es->parent) {
                        if (es->parent[ne] == pf->e1) {
                            if (pf->f1 != es->parent[ne + 2*nelem] ||
                                pf->e2 != es->parent[ne + nelem] ||
                                pf->f2 != es->parent[ne + 3*nelem]) np++;
                        }
                        else if (es->parent[ne] == pf->e2) {
                            if (pf->f2 != es->parent[ne + 2*nelem] ||
                                pf->e1 != es->parent[ne + nelem] ||
                                pf->f1 != es->parent[ne + 3*nelem]) np++;
                        }
                        else
                            np++;
                    }
                }
            }
            pe += nn;
        }

        if (verbose)
            printf ("    [internal,external] Faces=[%d,%d]\n", nint, next);
        if (nf)
            warning (1, "%d faces are not faces of the volume elements", nf);
        if (np)
            error ("%d faces have invalid parent data", np);
        if (es->parent) {
            free (es->parent);
            es->parent = NULL;
        }
    }
}

/*-----------------------------------------------------------------------*/

static void check_struct_interface (ZONE *z, CGNS_ENUMT(PointSetType_t) ptype,
    cgsize_t npts, cgsize_t *pts, int bndry)
{
    int id, n1, n2, n3, nerr1, nerr2;
    cgsize_t n, np;

    nerr1 = nerr2 = 0;

    /* vertices */

    if (ptype == CGNS_ENUMV(PointList)) {
        for (n = 0, np = 0; np < npts; np++) {
            n1 = n2 = 0;
            for (id = 0; id < z->idim; id++) {
                if (pts[n] < 1 || pts[n] > z->dims[0][id])
                    n1++;
                else if (pts[n] == 1 || pts[n] == z->dims[0][id])
                    n2++;
                n++;
            }
            if (n1) nerr1++;
            if (!n2) nerr2++;
        }
        if (nerr1)
            error ("%d points are out of range", nerr1);
        if (nerr2 && bndry)
            warning (2, "%d points are not boundary points", nerr2);
        return;
    }

    /* faces */

    if (bndry) {
        for (n = 0, np = 0; np < npts; np++) {
            n1 = n2 = n3 = 0;
            for (id = 0; id < z->idim; id++) {
                if (pts[n] < 1 || pts[n] > z->dims[0][id])
                    n1++;
                else if (pts[n] == 1)
                    n2++;
                else if (pts[n] == z->dims[0][id])
                    n3++;
                n++;
            }
            if (n1 || n3 > 1) nerr1++;
            if (!n2 && !n3) nerr2++;
        }
    }

    /* cells */

    else {
        for (n = 0, np = 0; np < npts; np++) {
            n1 = 0;
            for (id = 0; id < z->idim; id++) {
                if (pts[n] < 1 || pts[n] >= z->dims[0][id])
                    n1++;
                n++;
            }
            if (n1) nerr1++;
        }
    }

    if (nerr1)
        error ("%d elements are out of range", nerr1);
    if (nerr2)
        warning (2, "%d elements are not boundary elements", nerr2);
}

/*-----------------------------------------------------------------------*/

static void check_unstruct_interface (ZONE *z, CGNS_ENUMT(PointSetType_t) ptype,
    cgsize_t npts, cgsize_t *pts, int bndry)
{
    int dim, nn, id, nerr1, nerr2, nerr3;
    cgsize_t n, *nodes;
    FACE *face, *pf;

    nerr1 = nerr2 = nerr3 = 0;

    /* vertices */

    if (ptype == CGNS_ENUMV(PointList)) {
        for (n = 0; n < npts; n++) {
            if (pts[n] < 1 || pts[n] > z->nnodes) nerr1++;
        }
        if (nerr1)
            error ("%d node numbers are out of range", nerr1);
        if (z->nextnodes && bndry) {
            for (n = 0; n < npts; n++) {
                if (!find_extnode (z, pts[n])) nerr2++;
            }
            if (nerr2)
                warning (2, "%d nodes are not boundary nodes", nerr2);
        }
        return;
    }

    /* elements */

    dim = bndry ? CellDim - 1 : CellDim;

    for (n = 0; n < npts; n++) {
        nodes = find_element (z, pts[n], &id, &nn);
        if (nodes == NULL)
            nerr1++;
        else {
            if (id != dim) nerr2++;
            if (z->faces != NULL && id == 2 && nn > 2) {
                face = new_face (nn, nodes);
                pf = (FACE *) HashFind (z->faces, face);
                if (pf == NULL || pf->e2) nerr3++;
                free (face);
            }
        }
    }
    if (nerr1)
        error ("%d element numbers are out of range", nerr1);
    if (nerr2)
        warning (1, "%d elements have invalid dimension", nerr2);
    if (nerr3 && bndry)
        warning (2, "%d elements are not boundary elements", nerr3);
}

/*-----------------------------------------------------------------------*/

static cgsize_t check_interface (ZONE *z, CGNS_ENUMT(PointSetType_t) ptype,
    CGNS_ENUMT(GridLocation_t) location, cgsize_t npts, cgsize_t *pts, int bndry)
{
    cgsize_t np, *p;

    if (ptype == CGNS_ENUMV(PointListDonor)) ptype = CGNS_ENUMV(PointList);
    if (ptype == CGNS_ENUMV(CellListDonor)) ptype = CGNS_ENUMV(ElementList);
    if (ptype != CGNS_ENUMV(PointRange)   && ptype != CGNS_ENUMV(PointList) &&
        ptype != CGNS_ENUMV(ElementRange) && ptype != CGNS_ENUMV(ElementList)) {
        error ("invalid point type");
        return 0;
    }
    if (location < CGNS_ENUMV(Vertex) || location > CGNS_ENUMV(EdgeCenter)) {
        error ("invalid grid location");
        return 0;
    }
    if (z->idim < 1 || z->idim > 3) {
        warning (1, "can't evaluate since zone \"%s\" in invalid", z->name);
        return 0;
    }

    if (ptype == CGNS_ENUMV(PointRange) && location != CGNS_ENUMV(Vertex))
        ptype = CGNS_ENUMV(ElementRange);
    if (ptype == CGNS_ENUMV(PointList)  && location != CGNS_ENUMV(Vertex))
        ptype = CGNS_ENUMV(ElementList);

    if (ptype == CGNS_ENUMV(PointRange) || ptype == CGNS_ENUMV(ElementRange)) {
        int n;
        cgsize_t i, j, k;
        cgsize_t pmin[3], pmax[3];
        for (np = 1, n = 0; n < z->idim; n++) {
            if (pts[n] < pts[n+z->idim]) {
                pmin[n] = pts[n];
                pmax[n] = pts[n+z->idim];
            }
            else {
                pmin[n] = pts[n+z->idim];
                pmax[n] = pts[n];
            }
            np *= (pmax[n] - pmin[n] + 1);
        }
        p = (cgsize_t *) malloc ((size_t)(np * z->idim * sizeof(cgsize_t)));
        if (p == NULL)
            fatal_error("malloc failed for point/element list\n");
        n = 0;
        if (z->idim == 1) {
            for (i = pmin[0]; i <= pmax[0]; i++)
                p[n++] = i;
        }
        else if (z->idim == 2) {
            for (i = pmin[0]; i <= pmax[0]; i++) {
                for (j = pmin[1]; j <= pmax[1]; j++) {
                    p[n++] = i;
                    p[n++] = j;
                }
            }
        }
        else {
            for (i = pmin[0]; i <= pmax[0]; i++) {
                for (j = pmin[1]; j <= pmax[1]; j++) {
                    for (k = pmin[2]; k <= pmax[2]; k++) {
                        p[n++] = i;
                        p[n++] = j;
                        p[n++] = k;
                    }
                }
            }
        }
        if (ptype == CGNS_ENUMV(PointRange))
            ptype = CGNS_ENUMV(PointList);
        else
            ptype = CGNS_ENUMV(ElementList);
    }
    else {
        np = npts;
        p = pts;
    }

    if (z->type == CGNS_ENUMV(Structured))
        check_struct_interface (z, ptype, np, p, bndry);
    else
        check_unstruct_interface (z, ptype, np, p, bndry);

    if (p != pts) free (p);
    return np;
}

/*-----------------------------------------------------------------------*/

static CGNS_ENUMT(GridLocation_t) check_location (ZONE *z,
    CGNS_ENUMT(PointSetType_t) ptype, CGNS_ENUMT(GridLocation_t) location)
{
    switch (location) {
        case CGNS_ENUMV(Vertex):
            break;
        case CGNS_ENUMV(EdgeCenter):
            if (CellDim == 1) {
                if (FileVersion >= 3140) {
                    warning (1, "use CellCenter for CellDim=1"
                                " instead of EdgeCenter");
                }
                return CGNS_ENUMV(CellCenter);
            }
            return location;
        case CGNS_ENUMV(FaceCenter):
            if (z->type == CGNS_ENUMV(Structured)) {
                warning (2, "use [IJK]FaceCenter with Structured grids");
                return location;
            }
            if (CellDim == 2) {
                if (FileVersion >= 3140) {
                    warning (1, "use CellCenter for CellDim=2"
                                " instead of FaceCenter");
                }
                return CGNS_ENUMV(CellCenter);
            }
            return location;
        case CGNS_ENUMV(IFaceCenter):
            if (z->type != CGNS_ENUMV(Structured)) {
                error ("IFaceCenter only valid for Structured grids");
                return CGNS_ENUMV(FaceCenter);
            }
            break;
        case CGNS_ENUMV(JFaceCenter):
            if (z->type != CGNS_ENUMV(Structured) || z->idim < 2) {
                error ("JFaceCenter only valid for Structured grids"
                    " with CellDim > 1");
                return CGNS_ENUMV(FaceCenter);
            }
            break;
        case CGNS_ENUMV(KFaceCenter):
            if (z->type != CGNS_ENUMV(Structured) || z->idim < 3) {
                error ("KFaceCenter only valid for Structured grids"
                    " with CellDim > 2");
                return CGNS_ENUMV(FaceCenter);
            }
            break;
        case CGNS_ENUMV(CellCenter):
            if (z->type == CGNS_ENUMV(Structured) && FileVersion >= 2300) {
                warning (2, "use [IJK]FaceCenter location rather"
                    " than CellCenter");
            }
            return location;
        default:
            error ("invalid grid location");
            break;
    }
    if (ptype == CGNS_ENUMV(ElementRange) ||
        ptype == CGNS_ENUMV(ElementList)) {
        return (CellDim == 1 ? CGNS_ENUMV(Vertex) :
               (CellDim == 2 ? CGNS_ENUMV(EdgeCenter) :
                               CGNS_ENUMV(FaceCenter)));
    }
    return location;
}

/*-----------------------------------------------------------------------*/

static void check_BCdata (CGNS_ENUMT(BCType_t) bctype, int dirichlet, int neumann,
    cgsize_t size, int parclass, int *parunits, int indent)
{
    char name[33], *desc;
    int ierr, n, nd;
    int *punits, units[9], dataclass;
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(PointSetType_t) ptype;
    int hasl, hasp;
    cgsize_t npnts;

    hasl = read_gridlocation (&location);
    if (hasl && hasl != CG_NODE_NOT_FOUND)
        error_exit("cg_gridlocation_read");
    /* ptset only valid below a zone */
    if (cgnszone == 0)
        hasp = CG_NODE_NOT_FOUND;
    else {
        hasp = cg_ptset_info (&ptype, &npnts);
        if (hasp && hasp != CG_NODE_NOT_FOUND)
            error_exit("cg_ptset_info");
        if (hasp == CG_OK && hasl != CG_OK) {
            location = CGNS_ENUMV(Vertex);
            hasl = CG_OK;
        }
    }

    if (verbose) {
        print_indent (indent);
        printf ("BC Type=%s\n", cg_BCTypeName(bctype));
        print_indent (indent);
        printf ("Dirichlet Data=%s\n", dirichlet ? "yes" : "no");
        print_indent (indent);
        printf ("Neumann Data=%s\n", neumann ? "yes" : "no");
        if (hasl == CG_OK) {
            print_indent (indent);
            printf ("Grid Location=%s\n", cg_GridLocationName(location));
        }
        if (hasp == CG_OK) {
            print_indent (indent);
            printf ("Point Set Type=%s\n", cg_PointSetTypeName(ptype));
            print_indent (indent);
            printf ("Number Points=%ld\n", (long)npnts);
        }
        if (verbose > 1) {
            if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
            for (n = 1; n <= nd; n++) {
                if (cg_descriptor_read (n, name, &desc))
                    error_exit("cg_descriptor_read");
                if (desc != NULL) {
                    print_indent (indent);
                    printf ("Descriptor %s:\n%s\n", name, desc);
                    cg_free (desc);
                }
            }
        }
    }

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, indent);
        if (punits) print_units (punits, indent);
    }
    if (dataclass < 0) dataclass = parclass;
    if (punits == NULL) punits = parunits;

    ierr = cg_state_read (&desc);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_state_read");
    if (ierr == CG_OK) {
        print_indent (indent);
        puts ("checking reference state");
        if (desc != NULL) {
            if (verbose > 1) {
                print_indent (indent+2);
                printf ("Descriptor:%s\n", desc);
            }
            cg_free (desc);
        }
        fflush (stdout);
        go_relative ("ReferenceState_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, indent+2);
        go_relative ("..", 1, NULL);
    }

    check_user_data (dataclass, punits, indent);

    if (hasp == CG_OK) {
        cgsize_t *pts;
        ZONE *z = &Zones[cgnszone-1];

        print_indent (indent);
        puts ("checking BCDataSet interface");
        fflush (stdout);
        location = check_location (z, ptype, location);
        if (npnts < 1) {
            error ("number of points for Point Set less than 1");
            size = 0;
        }
        else if (npnts != 2 &&
                 (ptype == CGNS_ENUMV(PointRange) ||
                  ptype == CGNS_ENUMV(ElementRange))) {
            error ("number of points not 2 for Point/Element Range");
            size = 0;
        }
        else {
            pts = (cgsize_t *) malloc ((size_t)(z->idim * npnts * sizeof(cgsize_t)));
            if (NULL == pts)
                fatal_error("malloc failed for BCDataSet points\n");
            if (cg_ptset_read (pts)) error_exit("cg_ptset_read");
            if (ptype == CGNS_ENUMV(PointRange) ||
                ptype == CGNS_ENUMV(ElementRange)) {
                for (n = 0; n < z->idim; n++) {
                    if (pts[n] > pts[n+z->idim]) {
                        warning (1, "start value > end value for range");
                        break;
                    }
                }
            }
            size = check_interface (z, ptype, location, npnts, pts, location != CellCenter);
            free (pts);
        }
    }
    else {
        if (hasp == CG_OK)
            warning (1, "grid location should be used only with point set");
    }

    if (size) size = -size;  /* allow array size of 1 */

    if (dirichlet) {
        print_indent (indent);
        puts ("checking Dirichlet data");
        go_relative ("BCData_t", CGNS_ENUMV(Dirichlet), NULL);
        check_arrays (dataclass, punits, 0, size, indent+2);
        go_relative ("..", 1, NULL);
    }

    if (neumann) {
        print_indent (indent);
        puts ("checking Neumann data");
        go_relative ("BCData_t", CGNS_ENUMV(Neumann), NULL);
        check_arrays (dataclass, punits, 0, size, indent+2);
        go_relative ("..", 1, NULL);
    }
}

/*-----------------------------------------------------------------------*/

static void check_BCdataset (int nb, int nd, cgsize_t size,
    int parclass, int *parunits)
{
    char name[33];
    int dirichlet, neumann;
    CGNS_ENUMT(BCType_t) bctype;

    if (cg_dataset_read (cgnsfn, cgnsbase, cgnszone, nb, nd,
        name, &bctype, &dirichlet, &neumann))
        error_exit("cg_dataset_read");
    printf ("    checking BC data set \"%s\"\n", name);
    fflush (stdout);
    go_absolute ("Zone_t", cgnszone, "ZoneBC_t", 1, "BC_t", nb,
        "BCDataSet_t", nd, NULL);
    check_BCdata (bctype, dirichlet, neumann, size, parclass, parunits, 6);
}

/*-----------------------------------------------------------------------*/

static void check_BC (int nb, int parclass, int *parunits)
{
    char name[33], *desc;
    int n, nd, ierr, ndataset, hasl;
    cgsize_t npts, *pts, nrmlflag;
    int nrmlindex[3];
    int *punits, units[9], dataclass;
    void *nrmllist;
    CGNS_ENUMT(BCType_t) bctype;
    CGNS_ENUMT(PointSetType_t) ptype;
    CGNS_ENUMT(DataType_t) datatype;
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(WallFunctionType_t) wtype;
    CGNS_ENUMT(AreaType_t) atype;
    float area;
    ZONE *z = &Zones[cgnszone-1];

    if (cg_boco_info (cgnsfn, cgnsbase, cgnszone, nb, name, &bctype,
            &ptype, &npts, nrmlindex, &nrmlflag, &datatype, &ndataset))
        error_exit("cg_boco_info");
    printf ("  checking BC \"%s\"\n", name);
    if (verbose) {
        printf ("    BC Type=%s\n", cg_BCTypeName(bctype));
        printf ("    Point Set Type=%s\n", cg_PointSetTypeName(ptype));
        printf ("    Number Points=%ld\n", (long)npts);
    }
    fflush (stdout);

    go_absolute ("Zone_t", cgnszone, "ZoneBC_t", 1, "BC_t", nb, NULL);

    if ((FileVersion >= 1270 && FileVersion <= 2200) || FileVersion >= 3100) {
        if (ptype != CGNS_ENUMV(PointRange) && ptype != CGNS_ENUMV(PointList)) {
            error ("point set type not PointRange or PointList");
            return;
        }
    }
    else {
        if (ptype != CGNS_ENUMV(PointRange)   && ptype != CGNS_ENUMV(PointList) &&
            ptype != CGNS_ENUMV(ElementRange) && ptype != CGNS_ENUMV(ElementList)) {
            error ("point set type not PointRange, PointList, ElementRange"
                " or ElementList");
            return;
        }
    }

#if CGNS_VERSION > 3100 && !defined(CG_ALLOW_ELEMENTLIST_RANGE)
    if (cg_boco_gridlocation_read(cgnsfn, cgnsbase, cgnszone, nb, &location))
        error_exit("cg_boco_gridlocation_read");
    hasl = CG_OK;
#else
    hasl = read_gridlocation (&location);
    if (hasl && hasl != CG_NODE_NOT_FOUND)
        error_exit("cg_gridlocation_read");
#endif
    if (verbose && hasl == CG_OK)
        printf ("    Grid Location=%s\n", cg_GridLocationName(location));
    location = check_location (z, ptype, location);

    if (npts < 1) {
        error ("number of points is less than 1");
        return;
    }
    if ((FileVersion >= 1270 && FileVersion <= 2200) || FileVersion >= 3100) {
        if (ptype == CGNS_ENUMV(PointRange) && npts != 2) {
            error ("number of points is not 2 for PointRange");
            return;
        }
    }
    else {
        if ((ptype == CGNS_ENUMV(PointRange) ||
             ptype == CGNS_ENUMV(ElementRange)) && npts != 2) {
            error ("number of points is not 2 for PointRange/ElementRange");
            return;
        }
    }
    if (z->idim == 0) return;

    ierr = check_node ("\"int[IndexDimension]\"");

    if (verbose) {
        if (ierr == CG_OK) {
            printf ("    Normal Index=[%d", nrmlindex[0]);
            for (n = 1; n < z->idim; n++)
                printf (",%d", nrmlindex[n]);
            puts ("]");
        }
        if (nrmlflag) {
            puts ("    Normals Defined=yes");
            printf ("    Number Normals=%ld\n", (long)(nrmlflag / PhyDim));
        }
        else
            puts ("    Normals Defined=no");
    }
    fflush (stdout);

    if (ierr == CG_OK && z->type != CGNS_ENUMV(Structured))
        error ("normal index is only valid for Structured grids");
    if (nrmlflag) {
        if (datatype != CGNS_ENUMV(RealSingle) &&
            datatype != CGNS_ENUMV(RealDouble)) {
            error ("normal data type is not RealSingle or RealDouble");
            if (LibraryVersion < 2200) return;
        }
    }

    pts = (cgsize_t *) malloc ((size_t)(z->idim * npts * sizeof(cgsize_t)));
    if (NULL == pts)
        fatal_error("malloc failed for BC points\n");
    nrmllist = NULL;
    if (nrmlflag && LibraryVersion < 2200) {
        int n = (datatype == CGNS_ENUMV(RealSingle) ? sizeof(float) : sizeof(double));
        nrmllist = (void *) malloc ((size_t)(nrmlflag * n));
        if (nrmllist == NULL)
            fatal_error("malloc failed for BC normals\n");
    }
    if (cg_boco_read (cgnsfn, cgnsbase, cgnszone, nb, pts, nrmllist))
        error_exit("cg_boco_read");
    if (nrmllist) free (nrmllist);

    if (verbose && (ptype == CGNS_ENUMV(PointRange) ||
                    ptype == CGNS_ENUMV(ElementRange))) {
        printf ("    Range=[%ld:%ld", (long)pts[0], (long)pts[z->idim]);
        for (n = 1; n < z->idim; n++)
            printf (",%ld:%ld", (long)pts[n], (long)pts[n+z->idim]);
        puts ("]");
    }

    ierr = cg_famname_read (name);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_famname_read");
    if (ierr == CG_NODE_NOT_FOUND) {
      if (bctype == CGNS_ENUMV(FamilySpecified))
            warning (1,
                "BC Type is FamilySpecified but no family name is given");
    }
    else {
        if (verbose) printf ("    Family=\"%s\"\n", name);
        for (n = 0; n < NumFamily; n++) {
            if (0 == strcmp (name, Family[n])) break;
        }
        if (n == NumFamily &&
            (FileVersion >= 1200 || strcmp(name, "ORPHAN")))
            warning (bctype == CGNS_ENUMV(FamilySpecified) ? 1 : 2,
                "family name \"%s\" not found", name);
    }

    ierr = read_ordinal (&n);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_ordinal_read");
    if (ierr == CG_OK && verbose)
        printf ("    Ordinal=%d\n", n);

    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 4);
        if (punits) print_units (punits, 4);
    }
    if (dataclass < 0) dataclass = parclass;
    if (punits == NULL) punits = parunits;

    ierr = cg_state_read (&desc);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_state_read");
    if (ierr == CG_OK) {
        puts ("    checking reference state");
        if (desc != NULL) {
            if (verbose > 1)
                printf ("      Descriptor:%s\n", desc);
            cg_free (desc);
        }
        fflush (stdout);
        go_relative ("ReferenceState_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, 6);
        go_relative ("..", 1, NULL);
    }

    check_user_data (dataclass, punits, 4);

    puts ("    checking BC interface");
    fflush (stdout);
    if (ptype == CGNS_ENUMV(PointRange) || ptype == CGNS_ENUMV(ElementRange)) {
        for (n = 0; n < z->idim; n++) {
            if (pts[n] > pts[n+z->idim]) {
                warning (1, "start value > end value for range");
                break;
            }
        }
    }
    npts = check_interface (z, ptype, location, npts, pts, location != CellCenter);
    free (pts);

    /* BCDataSet */

    for (nd = 1; nd <= ndataset; nd++)
        check_BCdataset (nb, nd, npts, dataclass, punits);

    /* BCProperty */

    if (check_node ("BCProperty_t") == CG_OK) {
        puts ("    checking BC property");
        fflush (stdout);
        go_absolute ("Zone_t", cgnszone, "ZoneBC_t", 1, "BC_t", nb,
            "BCProperty_t", 1, NULL);
        if (verbose > 1) {
            if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
            for (n = 1; n <= nd; n++) {
                if (cg_descriptor_read (n, name, &desc))
                    error_exit("cg_descriptor_read");
                if (desc != NULL) {
                    printf ("      Descriptor %s:\n%s\n", name, desc);
                    cg_free (desc);
                }
            }
        }
        check_user_data (dataclass, punits, 6);

        ierr = cg_bc_wallfunction_read (cgnsfn, cgnsbase, cgnszone, nb,
            &wtype);
        if (ierr && ierr != CG_NODE_NOT_FOUND)
            error_exit("cg_bc_wallfunction_read");
        if (ierr == CG_OK) {
            puts ("      checking Wall Function property");
            go_relative ("WallFunction_t", 1, NULL);
            if (verbose) {
                printf ("        Wall Function Type=%s\n",
                    cg_WallFunctionTypeName (wtype));
                if (verbose > 1) {
                    if (cg_ndescriptors (&nd))
                        error_exit("cg_ndescriptors");
                    for (n = 1; n <= nd; n++) {
                        if (cg_descriptor_read (n, name, &desc))
                            error_exit("cg_descriptor_read");
                        if (desc != NULL) {
                            printf ("        Descriptor %s:\n%s\n", name, desc);
                            cg_free (desc);
                        }
                    }
                }
            }
            check_user_data (dataclass, punits, 8);
            go_relative ("..", 1, NULL);
        }

        ierr = cg_bc_area_read (cgnsfn, cgnsbase, cgnszone, nb,
            &atype, &area, name);
        if (ierr && ierr != CG_NODE_NOT_FOUND)
            error_exit("cg_bc_area_read");
        if (ierr == CG_OK) {
            puts ("      checking Area property");
            go_relative ("Area_t", 1, NULL);
            if (verbose) {
                printf ("        Area Type=%s\n", cg_AreaTypeName (atype));
                printf ("        Area=%g\n", area);
                printf ("        Region=\"%s\"\n", name);
                if (verbose > 1) {
                    if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
                    for (n = 1; n <= nd; n++) {
                        if (cg_descriptor_read (n, name, &desc))
                            error_exit("cg_descriptor_read");
                        if (desc != NULL) {
                            printf ("        Descriptor %s:\n%s\n", name, desc);
                            cg_free (desc);
                        }
                    }
                }
            }
            check_quantity (1, "SurfaceArea", dataclass, punits, 1, 8);
            check_user_data (dataclass, punits, 8);
            go_relative ("..", 1, NULL);
        }
    }
}

/*-----------------------------------------------------------------------*/

static void check_zoneBC (void)
{
    char name[33], *desc;
    int n, nb, ierr;
    int *punits, units[9], dataclass;

    if (cg_nbocos (cgnsfn, cgnsbase, cgnszone, &nb))
        error_exit("cg_nbocos");
    if (nb < 1) return;
    puts ("  checking boundary conditions");
    fflush (stdout);
    go_absolute ("Zone_t", cgnszone, "ZoneBC_t", 1, NULL);

    if (verbose > 1) {
        int nd;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 4);
        if (punits) print_units (punits, 4);
    }
    if (dataclass < 0) dataclass = Zones[cgnszone-1].dataclass;
    if (punits == NULL) punits = Zones[cgnszone-1].punits;

    ierr = cg_state_read (&desc);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_state_read");
    if (ierr == CG_OK) {
        puts ("    checking reference state");
        if (desc != NULL) {
            if (verbose > 1)
                printf ("      Descriptor:%s\n", desc);
            cg_free (desc);
        }
        fflush (stdout);
        go_relative ("ReferenceState_t", 1, NULL);
        check_arrays (dataclass, punits, 1, 0, 6);
    }

    check_user_data (dataclass, punits, 4);

    for (n = 1; n <= nb; n++)
        check_BC (n, dataclass, punits);
}

/*-----------------------------------------------------------------------*/

static void check_1to1 (int nzc, int nc)
{
    char name[33], dname[33], *desc;
    int ierr, n, nd, trans[3];
    cgsize_t range[6], drange[6];
    ZONE *z = &Zones[cgnszone-1], *dz;
    int ndim;
    cgsize_t dims[12];
    float center[3], angle[3], translate[3];
    CGNS_ENUMT(AverageInterfaceType_t) average;
    CGNS_ENUMT(DataType_t) datatype;

    if (cg_1to1_read (cgnsfn, cgnsbase, cgnszone, nc, name,
        dname, range, drange, trans)) error_exit("cg_1to1_read");
    printf ("  checking 1to1 connectivity \"%s\"\n", name);
    if (verbose) {
        printf ("    Range=[%ld:%ld", (long)range[0], (long)range[z->idim]);
        for (n = 1; n < z->idim; n++)
            printf (",%ld:%ld", (long)range[n], (long)range[n+z->idim]);
        puts ("]");
    }

    for (n = 0; n < z->idim; n++) {
        if (range[n] > range[n+z->idim]) {
            warning (1, "start value > end value for range");
            break;
        }
    }
    puts ("    checking 1to1 interface");
    fflush (stdout);
    check_interface (z, CGNS_ENUMV(PointRange), CGNS_ENUMV(Vertex), 2, range, 1);

    for (dz = NULL, n = 0; n < NumZones; n++) {
        if (0 == strcmp (dname, Zones[n].name)) {
            dz = &Zones[n];
            break;
        }
    }
    if (verbose)
        printf ("    Donor Name=%s\n", dname);
    if (dz == NULL) {
        error ("donor zone \"%s\" not found", dname);
        if (verbose)
            puts ("    Donor Range=<can't evaluate>");
    }
    else {
        if (verbose) {
            printf ("    Donor Range=[%ld:%ld", (long)drange[0], (long)drange[dz->idim]);
            for (n = 1; n < dz->idim; n++)
                printf (",%ld:%ld", (long)drange[n], (long)drange[n+dz->idim]);
            puts ("]");
        }
        puts ("    checking donor interface");
        fflush (stdout);
        check_interface (dz, CGNS_ENUMV(PointRange), CGNS_ENUMV(Vertex), 2, drange, 1);
    }

    if (verbose) {
        printf ("    Transform=[%d", trans[0]);
        for (n = 1; n < z->idim; n++)
            printf (",%d", trans[n]);
        puts ("]");
    }

    if (dz != NULL) {
        int id;
        cgsize_t dp, np = 1, dnp = 1;
        for (n = 0; n < z->idim; n++) {
            dp = range[n+z->idim] - range[n];
            if (dp < 0) dp = -dp;
            np *= (dp + 1);
        }
        for (n = 0; n < z->idim; n++) {
            dp = drange[n+dz->idim] - drange[n];
            if (dp < 0) dp = -dp;
            dnp *= (dp + 1);
        }
        if (np != dnp)
            error ("number of points is not the same as for the donor zone");
        for (n = 0; n < z->idim; n++) {
            np = range[n+z->idim] - range[n];
            id = abs (trans[n]);
            if ((id == 0 && np) || id > dz->idim)
                dnp = np + 1;
            else {
                id--;
                dnp = drange[id+dz->idim] - drange[id];
                if (trans[n] < 0) dnp = -dnp;
            }
            if (np != dnp) {
                error ("transform specification is invalid");
                break;
            }
        }
    }

    go_absolute ("Zone_t", cgnszone, "ZoneGridConnectivity_t", nzc,
        "GridConnectivity1to1_t", nc, NULL);

    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    ierr = read_ordinal (&n);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_ordinal_read");
    if (ierr == CG_OK && verbose)
        printf ("    Ordinal=%d\n", n);

    check_user_data (z->dataclass, z->punits, 4);

    if (check_node ("GridConnectivityProperty_t")) return;
    puts ("    checking grid connectivity property");
    fflush (stdout);
    go_relative ("GridConnectivityProperty_t", 1, NULL);
    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("      Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }
    check_user_data (z->dataclass, z->punits, 6);

    ierr = cg_1to1_periodic_read (cgnsfn, cgnsbase, cgnszone, nc,
        center, angle, translate);
    if (ierr && ierr != CG_NODE_NOT_FOUND)
        error_exit("cg_1to1_periodic_read");
    if (ierr == CG_OK) {
        int *punits, units[9], dataclass;

        puts ("      checking periodic property");
        go_relative ("Periodic_t", 1, NULL);
        if (verbose) {
            printf ("        Center=[%g", center[0]);
            for (n = 1; n < PhyDim; n++)
                printf (",%g", center[n]);
            printf ("]\n        Angle=[%g", angle[0]);
            for (n = 1; n < PhyDim; n++)
                printf (",%g", angle[n]);
            printf ("]\n        Translation=[%g", translate[0]);
            for (n = 1; n < PhyDim; n++)
                printf (",%g", translate[n]);
            puts ("]");
            if (verbose > 1) {
                if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
                for (n = 1; n <= nd; n++) {
                    if (cg_descriptor_read (n, name, &desc))
                        error_exit("cg_descriptor_read");
                    if (desc != NULL) {
                        printf ("        Descriptor %s:\n%s\n", name, desc);
                        cg_free (desc);
                    }
                }
            }
        }

        dataclass = read_dataclass ();
        punits = read_units (units);
        if (verbose) {
            if (dataclass >= 0) print_dataclass (dataclass, 8);
            if (punits) print_units (punits, 8);
        }
        if (dataclass < 0) dataclass = z->dataclass;
        if (punits == NULL) punits = z->punits;

        if (cg_narrays (&nd)) error_exit("cg_narrays");
        for (n = 1; n <= nd; n++) {
            if (cg_array_info (n, name, &datatype, &ndim, dims))
                error_exit("cg_array_info");
            printf ("        checking periodic data %s\n", name);
            check_quantity (n, name, dataclass, punits, 1, 8);
        }
        go_relative ("..", 1, NULL);
    }

    ierr = cg_1to1_average_read (cgnsfn, cgnsbase, cgnszone, nc, &average);
    if (ierr && ierr != CG_NODE_NOT_FOUND)
        error_exit("cg_1to1_average_read");
    if (ierr == CG_OK) {
        puts ("      checking average interface property");
        fflush (stdout);
        go_relative ("AverageInterface_t", 1, NULL);
        if (verbose) {
            printf ("        Interface Type=%s\n",
                cg_AverageInterfaceTypeName (average));
            if (verbose > 1) {
                if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
                for (n = 1; n <= nd; n++) {
                    if (cg_descriptor_read (n, name, &desc))
                        error_exit("cg_descriptor_read");
                    if (desc != NULL) {
                        printf ("        Descriptor %s:\n%s\n", name, desc);
                        cg_free (desc);
                    }
                }
            }
        }
        check_user_data (z->dataclass, z->punits, 8);
    }
}

/*-----------------------------------------------------------------------*/

static void check_conn (int nzc, int nc)
{
    char name[33], dname[33], *desc;
    int ierr, n, nd, ndim;
    cgsize_t npts, dnpts, dims[12];
    int interp;
    cgsize_t *pts, *dpts;
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(GridConnectivityType_t) ctype;
    CGNS_ENUMT(PointSetType_t) ptype, dptype;
    CGNS_ENUMT(ZoneType_t) dztype;
    CGNS_ENUMT(DataType_t) datatype;
    float center[3], angle[3], trans[3];
    CGNS_ENUMT(AverageInterfaceType_t) average;
    ZONE *z = &Zones[cgnszone-1], *dz;
    CGNS_ENUMT(DataType_t) ddatatype;

    if (cg_conn_info (cgnsfn, cgnsbase, cgnszone, nc, name,
            &location, &ctype, &ptype, &npts, dname, &dztype,
            &dptype, &ddatatype, &dnpts)) error_exit("cg_conn_info");
    printf ("  checking connectivity \"%s\"\n", name);
    if (verbose) {
        printf ("    Connectivity Type=%s\n",
            cg_GridConnectivityTypeName(ctype));
        printf ("    Grid Location=%s\n",
            cg_GridLocationName(location));
        printf ("    Point Set Type=%s\n",
            cg_PointSetTypeName(ptype));
        printf ("    Number Points=%ld\n", (long)npts);
        printf ("    Donor Zone=\"%s\"\n", dname);
        printf ("    Donor Zone Type=%s\n", cg_ZoneTypeName(dztype));
        printf ("    Donor Point Set Type=%s\n",
            cg_PointSetTypeName(dptype));
        printf ("    Donor Number Points=%ld\n", (long)dnpts);
    }
    fflush (stdout);

    go_absolute ("Zone_t", cgnszone, "ZoneGridConnectivity_t", nzc,
        "GridConnectivity_t", nc, NULL);
    interp = check_interpolants ();

    ierr = 0;
    if (ctype == CGNS_ENUMV(Overset)) {
      if (location != CGNS_ENUMV(Vertex) && location != CGNS_ENUMV(CellCenter))
            warning (1, "grid location should be Vertex or CellCenter");
    }
    else if (ctype == CGNS_ENUMV(Abutting) || ctype == CGNS_ENUMV(Abutting1to1)) {
        if (FileVersion < 2200) {
          if (location != CGNS_ENUMV(Vertex) && location != CGNS_ENUMV(CellCenter))
                warning (1,
                    "grid location should be Vertex or CellCenter");
        }
        else if (FileVersion >= 2300) {
            switch (location) {
                case CGNS_ENUMV(Vertex):
                    break;
                case CGNS_ENUMV(FaceCenter):
                    if (z->type == CGNS_ENUMV(Structured))
                        warning (2,
                            "use [IJK]FaceCenter with Structured grids");
                    break;
                case CGNS_ENUMV(IFaceCenter):
                    if (z->type != CGNS_ENUMV(Structured))
                        error ("IFaceCenter only valid for Structured grids");
                    break;
                case CGNS_ENUMV(JFaceCenter):
                    if (z->type != CGNS_ENUMV(Structured) || z->idim < 2)
                        error ("JFaceCenter only valid for Structured grids"
                            " with CellDim > 1");
                    break;
                case CGNS_ENUMV(KFaceCenter):
                    if (z->type != CGNS_ENUMV(Structured) || z->idim < 3)
                        error ("KFaceCenter only valid for Structured grids"
                            " with CellDim > 2");
                    break;
                case CGNS_ENUMV(CellCenter):
                    if (z->type == CGNS_ENUMV(Structured))
                        warning (2, "use [IJK]FaceCenter location rather"
                            " than CellCenter");
                    else
                        warning (2, "use FaceCenter location rather than"
                            " CellCenter");
                    break;
                default:
                    error ("grid location not Vertex,CellCenter,FaceCenter"
                        " or [IJK]FaceCenter");
            }
        }
        else {
            if (location != CGNS_ENUMV(Vertex) && location != CGNS_ENUMV(CellCenter) &&
                location != CGNS_ENUMV(FaceCenter))
                warning (1, "grid location should be Vertex, FaceCenter"
                    " or CellCenter");
        }
    }
    else {
        error ("connectivity type not Overset,Abutting or Abutting1to1");
        ierr++;
    }
    if (location < CGNS_ENUMV(Vertex) || location >= CGNS_ENUMV(EdgeCenter)) ierr++;

    if (ptype == CGNS_ENUMV(PointRange)) {
        if (npts != 2) {
            error ("number of points is not 2 for PointRange");
            ierr++;
        }
    }
    else if (ptype == CGNS_ENUMV(PointList)) {
        if (npts < 1) {
            error ("number of points is less than 1");
            ierr++;
        }
    }
    else {
        error ("point set type is not PointList or PointRange");
        ierr++;
    }

    for (dz = NULL, n = 0; n < NumZones; n++) {
        if (0 == strcmp (dname, Zones[n].name)) {
            dz = &Zones[n];
            break;
        }
    }
    if (dz == NULL) {
        error ("donor zone \"%s\" not found", dname);
        ierr++;
    }
    else {
        if (dztype != dz->type) {
            error ("returned donor zone type does not match actual zone type");
            ierr++;
        }
    }

    if (dptype == CGNS_ENUMV(PointListDonor)) {
        if (ctype != CGNS_ENUMV(Abutting1to1) && FileVersion >= 2000)
            warning (1, "PointListDonor should only be used for Abutting1to1");
        if (interp)
            warning (1, "InterpolantsDonor given for PointListDonor");
    }
    else if (dptype == CGNS_ENUMV(CellListDonor)) {
        if (interp) {
            if (LibraryVersion < 2000)
                warning (1,
                    "InterpolantsDonor given but they are not readable\n"
                    "        due to a bug in the CGNS library prior to 2.0");
        }
        else
            warning (1, "InterpolantsDonor not given for CellListDonor");
    }
    else {
        error ("donor point set type is not PointListDonor or CellListDonor");
        ierr++;
    }

    if (dnpts < 1) {
        error ("donor number of points is less than 1");
        ierr++;
    }
    if (ierr) npts = 0;

    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    ierr = read_ordinal (&n);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_ordinal_read");
    if (ierr == CG_OK && verbose)
        printf ("    Ordinal=%d\n", n);

    check_user_data (z->dataclass, z->punits, 4);

    if (npts && dnpts) {
        pts = (cgsize_t *) malloc ((size_t)(npts * z->idim * sizeof(cgsize_t)));
        if (LibraryVersion < 2200) {
            /* a bug in version prior to 2.2 causes the base cell dimension */
            /* to be used here, instead of the donor zone index dimension */
            dpts = (cgsize_t *) malloc ((size_t)(dnpts * CellDim * sizeof(cgsize_t)));
        }
        else
            dpts = (cgsize_t *) malloc ((size_t)(dnpts * dz->idim * sizeof(cgsize_t)));
        if (NULL == pts || NULL == dpts)
            fatal_error("malloc failed for connectivity points\n");
        if (cg_conn_read (cgnsfn, cgnsbase, cgnszone, nc, pts,
                CGNS_ENUMV(Integer), dpts)) error_exit("cg_conn_read");

        puts ("    checking connectivity interface");
        fflush (stdout);
        check_interface (z, ptype, location, npts, pts, ctype != CGNS_ENUMV(Overset));
        free (pts);

        puts ("    checking donor zone interface");
        fflush (stdout);
        check_interface (dz, dptype, location, dnpts, dpts, ctype != CGNS_ENUMV(Overset));
        free (dpts);
    }

    if (check_node ("GridConnectivityProperty_t")) return;
    puts ("    checking grid connectivity property");
    fflush (stdout);
    go_relative ("GridConnectivityProperty_t", 1, NULL);
    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("      Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }
    check_user_data (z->dataclass, z->punits, 6);

    ierr = cg_conn_periodic_read (cgnsfn, cgnsbase, cgnszone, nc,
        center, angle, trans);
    if (ierr && ierr != CG_NODE_NOT_FOUND)
        error_exit("cg_conn_periodic_read");
    if (ierr == CG_OK) {
        int *punits, units[9], dataclass;

        puts ("      checking periodic property");
        go_relative ("Periodic_t", 1, NULL);
        if (verbose) {
            printf ("        Center=[%g", center[0]);
            for (n = 1; n < PhyDim; n++)
                printf (",%g", center[n]);
            printf ("]\n        Angle=[%g", angle[0]);
            for (n = 1; n < PhyDim; n++)
                printf (",%g", angle[n]);
            printf ("]\n        Translation=[%g", trans[0]);
            for (n = 1; n < PhyDim; n++)
                printf (",%g", trans[n]);
            puts ("]");
            if (verbose > 1) {
                if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
                for (n = 1; n <= nd; n++) {
                    if (cg_descriptor_read (n, name, &desc))
                        error_exit("cg_descriptor_read");
                    if (desc != NULL) {
                        printf ("        Descriptor %s:\n%s\n", name, desc);
                        cg_free (desc);
                    }
                }
            }
        }

        dataclass = read_dataclass ();
        punits = read_units (units);
        if (verbose) {
            if (dataclass >= 0) print_dataclass (dataclass, 8);
            if (punits) print_units (punits, 8);
        }
        if (dataclass < 0) dataclass = z->dataclass;
        if (punits == NULL) punits = z->punits;

        if (cg_narrays (&nd)) error_exit("cg_narrays");
        for (n = 1; n <= nd; n++) {
            if (cg_array_info (n, name, &datatype, &ndim, dims))
                error_exit("cg_array_info");
            printf ("        checking periodic data %s\n", name);
            check_quantity (n, name, dataclass, punits, 1, 8);
        }
        go_relative ("..", 1, NULL);
    }

    ierr = cg_conn_average_read (cgnsfn, cgnsbase, cgnszone, nc, &average);
    if (ierr && ierr != CG_NODE_NOT_FOUND)
        error_exit("cg_conn_average_read");
    if (ierr == CG_OK) {
        puts ("      checking average interface property");
        fflush (stdout);
        go_relative ("AverageInterface_t", 1, NULL);
        if (verbose) {
            printf ("        Interface Type=%s\n",
                cg_AverageInterfaceTypeName (average));
            if (verbose > 1) {
                if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
                for (n = 1; n <= nd; n++) {
                    if (cg_descriptor_read (n, name, &desc))
                        error_exit("cg_descriptor_read");
                    if (desc != NULL) {
                        printf ("        Descriptor %s:\n%s\n", name, desc);
                        cg_free (desc);
                    }
                }
            }
        }
        check_user_data (z->dataclass, z->punits, 8);
    }
}

/*-----------------------------------------------------------------------*/

static void check_hole (int nzc, int nh)
{
    char name[33];
    int ierr, n, nptsets;
    cgsize_t np, npts;
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(PointSetType_t) ptype;
    ZONE *z = &Zones[cgnszone-1];

    if (cg_hole_info (cgnsfn, cgnsbase, cgnszone, nh, name,
        &location, &ptype, &nptsets, &npts)) error_exit("cg_hole_info");
    printf ("  checking overset hole \"%s\"\n", name);
    if (verbose) {
        printf ("    Grid Location=%s\n",
            cg_GridLocationName(location));
        printf ("    Point Set Type=%s\n",
            cg_PointSetTypeName(ptype));
        printf ("    Number Point Sets=%d\n", nptsets);
    }
    fflush (stdout);

    ierr = 0;
    np = 0;
    if (location != CGNS_ENUMV(Vertex) && location != CGNS_ENUMV(CellCenter)) {
        error ("location not Vertex or CellCenter");
        ierr++;
    }
    if (ptype == CGNS_ENUMV(PointRange)) {
        if (nptsets < 1)
            error ("nptsets must be greater then 0 for PointRange");
        if (npts != 2 * nptsets)
            error ("npts not equal to 2 * nptsets for PointRange");
        np = 2 * nptsets;
    }
    else if (ptype == CGNS_ENUMV(PointList)) {
        if (nptsets != 1)
            error ("nptsets must be 1 for PointList");
        if (npts < 1)
            error ("npts is less than 1 for PointList");
        nptsets = 1;
        np = npts;
    }
    else {
        error ("point set type not PointList or PointRange");
        ierr++;
    }

    if (!ierr && np > 0) {
        cgsize_t *pnts = (cgsize_t *) malloc ((size_t)(np * z->idim * sizeof(cgsize_t)));
        if (pnts == NULL)
            fatal_error("malloc failed for hole data\n");
        if (cg_hole_read (cgnsfn, cgnsbase, cgnszone, nh, pnts))
            error_exit("cg_hole_read");
        if (ptype == CGNS_ENUMV(PointRange)) npts = 2;
        for (np = 0, n = 1; n <= nptsets; n++) {
            printf ("    checking point set %d interface\n", n);
            fflush (stdout);
            check_interface (z, ptype, location, npts, &pnts[np], CellDim);
            np += npts * z->idim;
        }
        free (pnts);
    }

    go_absolute ("Zone_t", cgnszone, "ZoneGridConnectivity_t", nzc,
        "OversetHoles_t", nh, NULL);

    if (verbose > 1) {
        int nd;
        char *desc;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    check_user_data (z->dataclass, z->punits, 4);
}

/*-----------------------------------------------------------------------*/

static void check_connectivity (int nzc)
{
    int ierr, n, nc;
    char zcname[33];

    if (cg_zconn_set(cgnsfn, cgnsbase, cgnszone, nzc))
        error_exit("cg_zconn_set");
    if (cg_zconn_read(cgnsfn, cgnsbase, cgnszone, nzc, zcname))
        error_exit("cg_zconn_set");
    strcpy (ZoneConn[nzc-1], zcname);
    printf ("  checking zone grid connectivty \"%s\"\n", zcname);
    fflush (stdout);

    ierr = cg_goto (cgnsfn, cgnsbase, "Zone_t", cgnszone,
        "ZoneGridConnectivity_t", nzc, "end");
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_goto");
        return;
    }
    go_absolute ("Zone_t", cgnszone, "ZoneGridConnectivity_t", nzc, NULL);

    if (verbose > 1) {
        int nd;
        char name[33], *desc;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    check_user_data (Zones[cgnszone-1].dataclass, Zones[cgnszone-1].punits, 4);

    if (cg_n1to1 (cgnsfn, cgnsbase, cgnszone, &nc)) error_exit("cg_n1to1");
    for (n = 1; n <= nc; n++)
         check_1to1 (nzc, n);

    if (cg_nconns (cgnsfn, cgnsbase, cgnszone, &nc)) error_exit("cg_nconns");
    for (n = 1; n <= nc; n++)
        check_conn (nzc, n);

    if (cg_nholes (cgnsfn, cgnsbase, cgnszone, &nc)) error_exit("cg_nholes");
    for (n = 1; n <= nc; n++)
        check_hole (nzc, n);
}

/*-----------------------------------------------------------------------*/

static void check_arbitrary_motion (int na)
{
    char name[33];
    int ierr, n, nd, id, rind[6];
    int ndim;
    cgsize_t datasize, size, dims[12];
    int *punits, units[9], dataclass;
    CGNS_ENUMT(DataType_t) datatype;
    CGNS_ENUMT(ArbitraryGridMotionType_t) type;
    CGNS_ENUMT(GridLocation_t) location;
    ZONE *z = &Zones[cgnszone-1];

    if (cg_arbitrary_motion_read (cgnsfn, cgnsbase, cgnszone,
            na, name, &type)) error_exit("cg_arbitrary_motion_read");
    strcpy (ArbitraryGrid[na-1], name);
    printf ("  checking arbitrary motion \"%s\"\n", name);
    if (verbose)
        printf ("    Arbitrary Motion Type=%s\n",
            cg_ArbitraryGridMotionTypeName (type));
    fflush (stdout);

    go_absolute ("Zone_t", cgnszone, "ArbitraryGridMotion_t", na, NULL);

    /* grid location */

    ierr = read_gridlocation (&location);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_gridlocation_read");
        location = CGNS_ENUMV(Vertex);
    }
    else {
        if (verbose)
            printf ("    Grid Location=%s\n", cg_GridLocationName(location));
    }

    /* rind */

    ierr = read_rind (rind);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_rind_read");
        for (n = 0; n < 6; n++)
            rind[n] = 0;
    }
    else {
        if (verbose) {
            printf ("    Rind=[%d", rind[0]);
            for (n = 1; n < 2 * z->idim; n++)
                printf (",%d", rind[n]);
            puts ("]");
        }
        if (z->type == CGNS_ENUMV(Unstructured) && FileVersion < 2400)
            error ("rind not valid for unstructured zones");
    }

    /* descriptors */

    if (verbose > 1) {
        char *desc;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    /* dataclass and dimensional units */

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 4);
        if (punits) print_units (punits, 4);
    }
    if (dataclass < 0) dataclass = z->dataclass;
    if (punits == NULL) punits = z->punits;

    /* get grid data */

    datasize = get_data_size (z, location, rind);

    if (cg_narrays (&nd)) error_exit("cg_narrays");
    if (nd == 0 && type != CGNS_ENUMV(DeformingGrid))
        warning (1, "grid velocity data is missing");

    for (n = 1; n <= nd; n++) {
        if (cg_array_info (n, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        printf ("    checking grid velocity \"%s\"\n", name);
        fflush (stdout);
        for (size = 1, id = 0; id < ndim; id++)
            size *= dims[id];
        if (ndim != z->idim || size < 1 ||
            (datasize && size != datasize))
            error ("bad dimension values");
        check_quantity (n, name, dataclass, punits, 1, 6);
    }

    check_user_data (dataclass, punits, 4);
}

/*-----------------------------------------------------------------------*/

static void check_rigid_motion (int nr)
{
    char name[33];
    int n, nd, i, ndim;
    cgsize_t size, dims[12];
    int dataclass, *punits, units[9];
    CGNS_ENUMT(RigidGridMotionType_t) type;
    CGNS_ENUMT(DataType_t) datatype;
    ZONE *z = &Zones[cgnszone-1];

    if (cg_rigid_motion_read (cgnsfn, cgnsbase, cgnszone,
            nr, name, &type)) error_exit("cg_rigid_motion_read");
    strcpy (RigidGrid[nr-1], name);
    printf ("  checking rigid motion \"%s\"\n", name);
    if (verbose)
        printf ("    Rigid Motion Type=%s\n",
            cg_RigidGridMotionTypeName (type));
    fflush (stdout);

    go_absolute ("Zone_t", cgnszone, "RigidGridMotion_t", nr, NULL);

    /* descriptors */

    if (verbose > 1) {
        char *desc;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    /* dataclass and dimensional units */

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 4);
        if (punits) print_units (punits, 4);
    }
    if (dataclass < 0) dataclass = z->dataclass;
    if (punits == NULL) punits = z->punits;

    if (cg_narrays (&nd)) error_exit("cg_narrays");
    for (n = 1; n <= nd; n++) {
        if (cg_array_info (n, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        printf ("    checking rigid grid data \"%s\"\n", name);
        fflush (stdout);
        if (0 == strcmp (name, "OriginLocation")) {
            if (ndim != 2 || dims[0] != PhyDim || dims[1] != 2)
                error ("bad dimension values");
            check_quantity (n, name, dataclass, punits, 1, 6);
        }
        else if (0 == strcmp (name, "RigidRotationAngle") ||
                 0 == strcmp (name, "RigidRotationRate") ||
                 0 == strcmp (name, "RigidVelocity")) {
            if (ndim != 1 || dims[0] != PhyDim)
                error ("invalid dimension");
            check_quantity (n, name, dataclass, punits, 1, 6);
        }
        else {
            for (size = 1, i = 0; i < ndim; i++)
                size *= dims[i];
            if (ndim < 1 || size < 1)
                error ("bad dimension values");
            check_quantity (n, name, dataclass, punits, -1, 6);
        }
    }

    check_user_data (dataclass, punits, 2);
}

/*-----------------------------------------------------------------------*/

static void check_discrete (int ndis)
{
    char name[33];
    int n, nd, id, ierr, rind[6];
    int ndim;
    cgsize_t datasize, size, dims[12];
    int *punits, units[9], dataclass;
    CGNS_ENUMT(DataType_t) datatype;
    CGNS_ENUMT(GridLocation_t) location;
    ZONE *z = &Zones[cgnszone-1];

    if (cg_discrete_read (cgnsfn, cgnsbase, cgnszone, ndis, name))
        error_exit("cg_discrete_read");
    printf ("  checking discrete data \"%s\"\n", name);
    fflush (stdout);

    go_absolute ("Zone_t", cgnszone, "DiscreteData_t", ndis, NULL);

    /* grid location */

    ierr = read_gridlocation (&location);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_gridlocation_read");
        location = CGNS_ENUMV(Vertex);
    }
    else {
        if (verbose)
            printf ("    Grid Location=%s\n", cg_GridLocationName(location));
    }

    /* rind */

    ierr = read_rind (rind);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_rind_read");
        for (n = 0; n < 6; n++)
            rind[n] = 0;
    }
    else {
        if (verbose) {
            printf ("    Rind=[%d", rind[0]);
            for (n = 1; n < 2 * z->idim; n++)
                printf (",%d", rind[n]);
            puts ("]");
        }
        if (z->type == CGNS_ENUMV(Unstructured) && FileVersion < 2400)
            error ("rind not valid for unstructured zones");
    }

    /* descriptors */

    if (verbose > 1) {
        char *desc;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    /* dataclass and dimensional units */

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 4);
        if (punits) print_units (punits, 4);
    }
    if (dataclass < 0) dataclass = z->dataclass;
    if (punits == NULL) punits = z->punits;

    /* get discrete data */

    datasize = get_data_size (z, location, rind);

    if (cg_narrays (&nd)) error_exit("cg_narrays");
    if (nd == 0)
        warning (2, "no discrete data arrays defined");

    for (n = 1; n <= nd; n++) {
        if (cg_array_info (n, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        printf ("    checking discrete data array \"%s\"\n", name);
        fflush (stdout);
        for (size = 1, id = 0; id < ndim; id++)
            size *= dims[id];
        if (ndim != z->idim || size < 1 ||
            (datasize && size != datasize))
            error ("bad dimension values");
        check_quantity (n, name, dataclass, punits, -1, 6);
    }

    check_user_data (dataclass, punits, 4);
}

/*-----------------------------------------------------------------------*/

static void check_solution (int ns)
{
    char name[33];
    int n, nf, id, ierr, rind[6];
    int ndim;
    cgsize_t datasize, size, dims[12];
    int *punits, units[9], dataclass;
    CGNS_ENUMT(DataType_t) datatype;
    CGNS_ENUMT(GridLocation_t) location;
    ZONE *z = &Zones[cgnszone-1];

    if (cg_sol_info (cgnsfn, cgnsbase, cgnszone, ns, name, &location))
        error_exit("cg_sol_info");
    strcpy (FlowSolution[ns-1], name);
    printf ("  checking solution \"%s\"\n", name);
    if (verbose)
        printf ("    Grid Location=%s\n", cg_GridLocationName (location));
    fflush (stdout);

#if 0
    if (location != CGNS_ENUMV(Vertex) && location != CGNS_ENUMV(CellCenter)) {
        if (z->type == CGNS_ENUMV(Unstructured))
            error ("grid location nust be Vertex or CellCenter for"
                " unstructured zones");
        else if (location != CGNS_ENUMV(IFaceCenter) &&
                 location != CGNS_ENUMV(JFaceCenter) &&
                 location != CGNS_ENUMV(KFaceCenter))
            error ("grid location not Vertex,CellCenter or [IJK]FaceCenter");
    }
#else
    location = check_location (z, CGNS_ENUMV(PointSetTypeNull), location);
#endif

    go_absolute ("Zone_t", cgnszone, "FlowSolution_t", ns, NULL);

    /* rind */

    ierr = read_rind (rind);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_rind_read");
        for (n = 0; n < 6; n++)
            rind[n] = 0;
    }
    else {
        if (verbose) {
            printf ("    Rind=[%d", rind[0]);
            for (n = 1; n < 2 * z->idim; n++)
                printf (",%d", rind[n]);
            puts ("]");
        }
        if (z->type == CGNS_ENUMV(Unstructured) && FileVersion < 2400)
            error ("rind not valid for unstructured zones");
    }

    /* descriptors */

    if (verbose > 1) {
        char *desc;
        int nd;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    /* dataclass and dimensional units */

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 4);
        if (punits) print_units (punits, 4);
    }
    if (dataclass < 0) dataclass = z->dataclass;
    if (punits == NULL) punits = z->punits;

    /* get solution data size */

    datasize = get_data_size (z, location, rind);

    /* read solution data as arrays to get size */

    if (cg_nfields (cgnsfn, cgnsbase, cgnszone, ns, &nf))
        error_exit("cg_nfields");
    if (nf == 0)
        warning (2, "no solution data arrays defined");

    for (n = 1; n <= nf; n++) {
        if (cg_array_info (n, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        printf ("    checking solution field \"%s\"\n", name);
        fflush (stdout);
        for (size = 1, id = 0; id < ndim; id++)
            size *= dims[id];
        if (ndim != z->idim || size < 1 ||
            (datasize && size != datasize))
            error ("bad dimension values");
        check_quantity (n, name, dataclass, punits, 1, 6);
    }

    /* user data */

    check_user_data (dataclass, punits, 4);
}

/*-----------------------------------------------------------------------*/

static cgsize_t subreg_size(int dim, int isBC, char *subname)
{
    char name[33], dname[33], *p;
    int n, nbc, ni, nc, nzc, nn, ns, nsets, trans[3];
    cgsize_t npnts, dnpnts, nsize, range[6], drange[6];
    CGNS_ENUMT(BCType_t) bctype;
    CGNS_ENUMT(PointSetType_t) ptype, dptype;
    CGNS_ENUMT(DataType_t) dtype;
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(GridConnectivityType_t) gctype;
    CGNS_ENUMT(ZoneType_t) dztype;

    if (isBC) {
        if (cg_nbocos(cgnsfn, cgnsbase, cgnszone, &nbc))
            error_exit("cg_nbocos");
        for (n = 1; n <= nbc; n++) {
            if (cg_boco_info(cgnsfn, cgnsbase, cgnszone, n, name,
                &bctype, &ptype, &npnts, trans, &nsize, &dtype, &nc))
                error_exit("cg_boco_info");
            if (0 == strcmp(name, subname)) {
                if (ptype == CGNS_ENUMV(PointRange) && npnts == 2) {
                    if (cg_boco_read(cgnsfn, cgnsbase, cgnszone,
                            n, range, NULL))
                        error_exit("cg_boco_read");
                    npnts = 1;
                    for (ni = 0; ni < dim; ni++)
                        npnts *= (range[ni+dim] - range[ni] + 1);
                }
                else if (ptype != CGNS_ENUMV(PointList)) {
                    error("BC for subregion not PointList/Range");
                    return 0;
                }
                if (npnts == 0) warning(1, "subregion size is 0");
                return npnts;
            }
        }
        error("BCRegionName not found");
        return 0;
    }

    /* check for ZoneConnectivity */

    if (cg_nzconns(cgnsfn, cgnsbase, cgnszone, &nzc))
        error_exit("cg_nzconns");
    if (nzc < 1) {
        error("GridConnectivetyRegionName not found");
        return 0;
    }

    if ((p = strchr(subname, '/')) == NULL) {
        p = subname;
        nzc = 1;
    }
    else {
        *p = 0;
        for (n = 1; n <= nzc; n++) {
            if (cg_zconn_read(cgnsfn, cgnsbase, cgnszone, n, name))
                error_exit("cg_zconn_read");
            if (0 == strcmp(subname, name)) break;
        }
        if (n > nzc) {
            error("ZoneConnectivity not found");
            return 0;
        }
        nzc = n;
        *p++ = '/';
    }

    if (cg_zconn_set(cgnsfn, cgnsbase, cgnszone, nzc))
        error_exit("cg_zconn_set");

    /* check 1to1 connections */

    if (cg_n1to1(cgnsfn, cgnsbase, cgnszone, &nc))
        error_exit("cg_n1to1");
    for (n = 1; n <= nc; n++) {
        if (cg_1to1_read(cgnsfn, cgnsbase, cgnszone, n, name,
                dname, range, drange, trans))
            error_exit("cg_1to1_read");
        if (0 == strcmp(p, name)) {
            npnts = 1;
            for (ni = 0; ni < dim; ni++)
                npnts *= (range[ni+dim] - range[ni] + 1);
            return npnts;
        }
    }

    /* check general connections */

    if (cg_nconns(cgnsfn, cgnsbase, cgnszone, &nc))
        error_exit("cg_nconns");
    for (n = 1; n <= nc; n++) {
        if (cg_conn_info(cgnsfn, cgnsbase, cgnszone, n, name,
                &location, &gctype, &ptype, &npnts, dname,
                &dztype, &dptype, &dtype, &dnpnts))
            error_exit("cg_conn_info");
        if (0 == strcmp(p, name)) {
            if (ptype == CGNS_ENUMV(PointRange) && npnts == 2) {
                if (cg_conn_read_short(cgnsfn, cgnsbase, cgnszone, n, range))
                    error_exit("cg_conn_read_short");
                npnts = 1;
                for (ni = 0; ni < dim; ni++)
                    npnts *= (range[ni+dim] - range[ni] + 1);
            }
            else if (ptype != CGNS_ENUMV(PointList)) {
                error("connection for subregion not PointList/Range");
                return 0;
            }
            if (npnts == 0) warning(1, "subregion size is 0");
            return npnts;
        }
    }

    /* check holes */

    if (cg_nholes(cgnsfn, cgnsbase, cgnszone, &nc))
        error_exit("cg_nholes");
    for (n = 1; n <= nc; n++) {
        if (cg_hole_info(cgnsfn, cgnsbase, cgnszone, n, name,
                &location, &ptype, &nsets, &npnts))
            error_exit("cg_hole_info");
        if (0 == strcmp(p, name)) {
            if (ptype == CGNS_ENUMV(PointRange) && npnts == 2) {
                cgsize_t *pnts;
                pnts = (cgsize_t *)malloc(2 * nsets * dim * sizeof(cgsize_t));
                if (pnts == NULL)
                    fatal_error("subreg_size:malloc failed for hole data\n");
                if (cg_hole_read(cgnsfn, cgnsbase, cgnszone, n, pnts))
                    error_exit("cg_hole_read");
                npnts = 1;
                nn = 0;
                for (ns = 0; ns < nsets; ns++) {
                    for (ni = 0; ni < dim; ni++) {
                        npnts *= (range[nn+dim] - range[nn] + 1);
                        nn++;
                    }
                    nn += dim;
                }
                free(pnts);
            }
            else if (ptype != CGNS_ENUMV(PointList)) {
                error("hole for subregion not PointList/Range");
                return 0;
            }
            if (npnts == 0) warning(1, "subregion size is 0");
            return npnts;
        }
    }

    error("GridConnectivetyRegionName not found");
    return 0;
}

/*-----------------------------------------------------------------------*/

static void check_subreg (int ns)
{
    char name[33], *sname;
    int dimension, bclen, gclen, ierr;
    int n, nd, rind[6];
    int *punits, units[9], dataclass;
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(PointSetType_t) ptype;
    cgsize_t npnts, datasize;
    ZONE *z = &Zones[cgnszone-1];

    if (cg_subreg_info (cgnsfn, cgnsbase, cgnszone, ns, name,
        &dimension, &location, &ptype, &npnts, &bclen, &gclen))
        error_exit("cg_subreg_info");
    strcpy (ZoneSubReg[ns-1], name);
    printf ("  checking subregion \"%s\"\n", name);

    if (verbose)
        printf ("    Dimension=%d\n", dimension);
    if (dimension < 1 || dimension > CellDim)
        error("subregion dimension invalid");

    if (verbose)
        printf ("    Grid Location=%s\n", cg_GridLocationName(location));
    if (dimension == 1) {
        if (CellDim == 1) {
            if (location != CGNS_ENUMV(Vertex) &&
                location != CGNS_ENUMV(CellCenter))
                error ("grid location nust be Vertex or CellCenter");
        }
        else {
            if (location != CGNS_ENUMV(Vertex) &&
                location != CGNS_ENUMV(EdgeCenter))
                error ("grid location nust be Vertex or EdgeCenter");
        }
    }
    else if (dimension == 2) {
        if (CellDim == 2) {
            if (location != CGNS_ENUMV(Vertex) &&
                location != CGNS_ENUMV(EdgeCenter) &&
                location != CGNS_ENUMV(CellCenter))
                error ("grid location nust be Vertex, EdgeCenter or CellCenter");
        }
        else {
            if (location != CGNS_ENUMV(Vertex) &&
                location != CGNS_ENUMV(EdgeCenter) &&
                location != CGNS_ENUMV(CellCenter) &&
                location != CGNS_ENUMV(FaceCenter) &&
                location != CGNS_ENUMV(IFaceCenter) &&
                location != CGNS_ENUMV(JFaceCenter) &&
                location != CGNS_ENUMV(KFaceCenter))
                error ("grid location nust be Vertex, EdgeCenter, CellCenter"
                       " or [IJK]FaceCenter");
        }
    }

    datasize = 0;
    if (ptype == CGNS_ENUMV(PointSetTypeNull)) {
        if (bclen) {
            if (gclen)
                error ("both BCRegionName and GridConnectivityName specified");
            else {
                sname = (char *)malloc(bclen + 1);
                if (sname == NULL)
                    fatal_error("malloc failed for BCRegionName\n");
                if (cg_subreg_bcname_read(cgnsfn, cgnsbase, cgnszone, ns, sname))
                    error_exit("cg_subreg_bcname_read");
                if (verbose)
                    printf ("    BC Region Name=%s\n", sname);
                datasize = subreg_size(z->idim, 1, sname);
                free(sname);
            }
        }
        else if (gclen) {
            sname = (char *)malloc(gclen + 1);
            if (sname == NULL)
                fatal_error("malloc failed for GridConnectivityName\n");
            if (cg_subreg_gcname_read(cgnsfn, cgnsbase, cgnszone, ns, sname))
                error_exit("cg_subreg_gcname_read");
            if (verbose)
                printf ("    Grid Connectivity Name=%s\n", sname);
            datasize = subreg_size(z->idim, 0, sname);
            free(sname);
        }
        else {
            error ("one of PointList, PointRange, BCRegionName or "
                   "GridConnectivityName not specified");
        }
    }
    else if (ptype == CGNS_ENUMV(PointRange) ||
             ptype == CGNS_ENUMV(PointList)) {
        if (verbose) {
            printf ("    Point Set Type=%s\n", cg_PointSetTypeName(ptype));
            printf ("    Number Points=%ld\n", (long)npnts);
        }
        if (bclen || gclen)
            error("both PointSetType and BCRegionName or GridConnectivityName specified");
        ierr = 0;
        if (ptype == CGNS_ENUMV(PointRange)) {
            if (npnts != 2) {
                ierr++;
                error ("npts not equal to 2 for PointRange");
            }
        }
        else {
            if (npnts < 1) {
                ierr++;
                error ("npts is less than 1 for PointList");
            }
        }
        if (!ierr) {
            cgsize_t *pnts = (cgsize_t *)malloc(npnts * z->idim * sizeof(cgsize_t));
            if (pnts == NULL)
                fatal_error("check_subreg:malloc failed for points\n");
            if (cg_subreg_ptset_read(cgnsfn, cgnsbase, cgnszone, ns, pnts))
                error_exit("cg_subreg_ptset_read");
            if (verbose && ptype == CGNS_ENUMV(PointRange)) {
                printf ("    Range=[%ld:%ld", (long)pnts[0], (long)pnts[z->idim]);
                for (n = 1; n < z->idim; n++)
                    printf (",%ld:%ld", (long)pnts[n], (long)pnts[n+z->idim]);
                puts ("]");
            }
            datasize = check_interface (z, ptype, location, npnts, pnts, 0);
            free(pnts);
        }
    }
    else {
        error ("point set type not PointList or PointRange");
    }

    go_absolute ("Zone_t", cgnszone, "ZoneSubRegion_t", ns, NULL);

    /* rind */

    ierr = read_rind (rind);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_rind_read");
        for (n = 0; n < 6; n++)
            rind[n] = 0;
    }
    else {
        if (verbose) {
            printf ("    Rind=[%d", rind[0]);
            for (n = 1; n < 2 * z->idim; n++)
                printf (",%d", rind[n]);
            puts ("]");
        }
        if (z->type == CGNS_ENUMV(Unstructured) && FileVersion < 2400)
            error ("rind not valid for unstructured zones");
    }

    /* family */

    ierr = cg_famname_read (name);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_famname_read");
    if (ierr == CG_OK) {
        if (verbose) printf ("    Family=\"%s\"\n", name);
        for (n = 0; n < NumFamily; n++) {
            if (0 == strcmp (name, Family[n])) break;
        }
        if (n == NumFamily &&
            (FileVersion >= 1200 || strcmp(name, "ORPHAN")))
            warning (2, "family name \"%s\" not found", name);
    }

    /* descriptors */

    if (verbose > 1) {
        char *desc;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    /* dataclass and dimensional units */

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 4);
        if (punits) print_units (punits, 4);
    }
    if (dataclass < 0) dataclass = z->dataclass;
    if (punits == NULL) punits = z->punits;

    check_user_data (dataclass, punits, 4);

    /* check data arrays */
    
    if (datasize) {
        if (verbose) printf("    Data Size=%ld\n", (long)datasize);
        if (cg_narrays (&nd)) error_exit("cg_narrays");
        if (nd) check_arrays (dataclass, punits, 0, datasize, 4);
    }
}

/*-----------------------------------------------------------------------*/

static void check_zone_iter (void)
{
    char *p, *desc, name[33], buff[33];
    int ierr, n, na, nd, nn, ndim;
    cgsize_t dims[12], size;
    int dataclass, *punits, units[9];
    CGNS_ENUMT(DataType_t) datatype;
    ZONE *z = &Zones[cgnszone-1];

    go_absolute ("Zone_t", cgnszone, "ZoneIterativeData_t", 1, NULL);

    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("    Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 4);
        if (punits) print_units (punits, 4);
    }
    if (dataclass < 0) dataclass = z->dataclass;
    if (!punits) punits = z->punits;

    if (cg_narrays (&na)) error_exit("cg_narrays");
    for (n = 1; n <= na; n++) {
        if (cg_array_info (n, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        printf ("    checking zone iterative data \"%s\"\n", name);
        fflush (stdout);
        for (size = 1, nd = 0; nd < ndim; nd++)
            size *= dims[nd];
        if (0 == strcmp (name, "ArbitraryGridMotionPointers") ||
            0 == strcmp (name, "FlowSolutionPointers") ||
            0 == strcmp (name, "GridCoordinatesPointers") ||
            0 == strcmp (name, "RigidGridMotionPointers") ||
            0 == strcmp (name, "ZoneGridConnectivityPointers") ||
            0 == strcmp (name, "ZoneSubRegionPointers")) {
            if (ndim != 2 || dims[0] != 32 || size < 1 ||
               (NumSteps && dims[1] != NumSteps))
                error ("invalid dimension values");
            else {
                desc = (char *) malloc ((size_t)size);
                if (desc == NULL)
                    fatal_error("malloc failed for zone iter data\n");
                if (cg_array_read (n, desc)) error_exit("cg_array_read");
                ierr = 0;
                if (0 == strcmp (name, "ArbitraryGridMotionPointers")) {
                    for (nd = 0; nd < dims[1]; nd++) {
                        strncpy (buff, &desc[nd<<5], 32);
                        buff[32] = 0;
                        p = buff + strlen(buff);
                        while (--p >= buff && isspace(*p))
                            ;
                        *++p = 0;
                        for (nn = 0; nn < NumArbitraryGrid; nn++) {
                            if (0 == strcmp (buff, ArbitraryGrid[nn])) break;
                        }
                        if (nn == NumArbitraryGrid) ierr++;
                    }
                }
                else if (0 == strcmp (name, "FlowSolutionPointers")) {
                    for (nd = 0; nd < dims[1]; nd++) {
                        strncpy (buff, &desc[nd<<5], 32);
                        buff[32] = 0;
                        p = buff + strlen(buff);
                        while (--p >= buff && isspace(*p))
                            ;
                        *++p = 0;
                        for (nn = 0; nn < NumFlowSolution; nn++) {
                            if (0 == strcmp (buff, FlowSolution[nn])) break;
                        }
                        if (nn == NumFlowSolution) ierr++;
                    }
                }
                else if (0 == strcmp (name, "GridCoordinatesPointers")) {
                    for (nd = 0; nd < dims[1]; nd++) {
                        strncpy (buff, &desc[nd<<5], 32);
                        buff[32] = 0;
                        p = buff + strlen(buff);
                        while (--p >= buff && isspace(*p))
                            ;
                        *++p = 0;
                        for (nn = 0; nn < NumGridCoordinate; nn++) {
                            if (0 == strcmp (buff, GridCoordinate[nn])) break;
                        }
                        if (nn == NumGridCoordinate) ierr++;
                    }
                }
                else if (0 == strcmp (name, "RigidGridMotionPointers")) {
                    for (nd = 0; nd < dims[1]; nd++) {
                        strncpy (buff, &desc[nd<<5], 32);
                        buff[32] = 0;
                        p = buff + strlen(buff);
                        while (--p >= buff && isspace(*p))
                            ;
                        *++p = 0;
                        for (nn = 0; nn < NumRigidGrid; nn++) {
                            if (0 == strcmp (buff, RigidGrid[nn])) break;
                        }
                        if (nn == NumRigidGrid) ierr++;
                    }
                }
                else if (0 == strcmp (name, "ZoneGridConnectivityPointers")) {
                    for (nd = 0; nd < dims[1]; nd++) {
                        strncpy (buff, &desc[nd<<5], 32);
                        buff[32] = 0;
                        p = buff + strlen(buff);
                        while (--p >= buff && isspace(*p))
                            ;
                        *++p = 0;
                        for (nn = 0; nn < NumZoneConn; nn++) {
                            if (0 == strcmp (buff, ZoneConn[nn])) break;
                        }
                        if (nn == NumZoneConn) ierr++;
                    }
                }
                else {
                    for (nd = 0; nd < dims[1]; nd++) {
                        strncpy (buff, &desc[nd<<5], 32);
                        buff[32] = 0;
                        p = buff + strlen(buff);
                        while (--p >= buff && isspace(*p))
                            ;
                        *++p = 0;
                        for (nn = 0; nn < NumZoneSubReg; nn++) {
                            if (0 == strcmp (buff, ZoneSubReg[nn])) break;
                        }
                        if (nn == NumZoneSubReg) ierr++;
                    }
                }
                free (desc);
                if (ierr)
                    error ("%d %s are invalid", ierr, name);
            }
        }
        else {
            if (ndim < 1 || size < 1)
                error ("invalid dimension values");
            check_quantity (n, name, dataclass, punits, 0, 6);
        }
    }

    check_user_data (dataclass, punits, 4);
}

/*-----------------------------------------------------------------------*/

static void check_zone (void)
{
    char name[33], *desc;
    int n, nd, niter, ierr, eqset[7];
    float point[3], vector[3];
    ZONE *z = &Zones[cgnszone-1];
    static char indexname[] = "IJK";

    printf ("\nchecking zone \"%s\"\n", z->name);
    fflush (stdout);

    if (z->type == CGNS_ENUMV(Structured)) {
        for (n = 0; n < CellDim; n++) {
            if (z->dims[0][n] < 2) {
                error ("number of points in %c-direction < 2", indexname[n]);
                z->idim = 0;
            }
            if (z->dims[1][n] != z->dims[0][n] - 1) {
                error ("number of cells in %c-direction is %d instead of %d",
                    indexname[n], z->dims[1][n], z->dims[0][n] - 1);
                z->dims[1][n] = z->dims[0][n] - 1;
            }
            if (z->dims[2][n] != 0)
                warning (1, "VertexSizeBoundary in %c-direction should be 0"
                    " for structured grid", indexname[n]);
        }
    }
    else if (z->type == CGNS_ENUMV(Unstructured)) {
        if (z->dims[0][0] < CellDim + 1) {
            error ("number of vertices < CellDim + 1");
            z->idim = 0;
        }
        if (z->dims[1][0] < 1) {
            warning (1, "number of cells < 1");
            z->idim = 0;
        }
        if (z->dims[2][0] > z->dims[0][0])
            error ("VertexBoundarySize > total number of vertices");
    }
    else
        error ("zone type is not Structured or Unstructured");
    if (z->idim == 0) return;

    /*----- descriptors -----*/

    go_absolute ("Zone_t", cgnszone, NULL);
    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("  Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    /*----- DataClass and DimensionalUnits -----*/

    z->dataclass = read_dataclass ();
    z->punits = read_units (z->units);
    if (verbose) {
        printf ("  Zone Type=%s\n", cg_ZoneTypeName(z->type));
        printf ("  Vertex Size=[%ld", (long)z->dims[0][0]);
        for (n = 1; n < z->idim; n++)
            printf (",%ld", (long)z->dims[0][n]);
        puts ("]");
        printf ("  Cell Size=[%ld", (long)z->dims[1][0]);
        for (n = 1; n < z->idim; n++)
            printf (",%ld", (long)z->dims[1][n]);
        puts ("]");
        if (z->dataclass >= 0) print_dataclass (z->dataclass, 2);
        if (z->punits) print_units (z->punits, 2);
    }
    if (z->dataclass < 0) z->dataclass = BaseClass;
    if (!z->punits) z->punits = pBaseUnits;

    /*----- FamilyName -----*/

    ierr = cg_famname_read (name);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_famname_read");
    if (ierr == CG_OK) {
        if (verbose) printf ("  Family=%s\n", name);
        for (n = 0; n < NumFamily; n++) {
            if (0 == strcmp (name, Family[n])) break;
        }
        if (n >= NumFamily &&
            (FileVersion >= 1200 || strcmp(name, "ORPHAN")))
            warning (2, "zone family name \"%s\" not found", name);
    }

    /*----- ReferenceState -----*/

    ierr = cg_state_read (&desc);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_state_read");
    if (ierr == CG_OK) {
        puts ("  checking reference state");
        if (desc != NULL) {
            if (verbose > 1)
                printf ("    Descriptor:%s\n", desc);
            cg_free (desc);
        }
        fflush (stdout);
        go_relative ("ReferenceState_t", 1, NULL);
        check_arrays (BaseClass, pBaseUnits, 1, 0, 4);
    }

    /*----- FlowEquationSet -----*/

    go_absolute ("Zone_t", cgnszone, NULL);
    ierr = cg_equationset_read (&eqset[0], &eqset[1], &eqset[2],
        &eqset[3], &eqset[4], &eqset[5], &eqset[6]);
    if (ierr && ierr != CG_NODE_NOT_FOUND)
        error_exit("cg_equationset_read");
    if (ierr == CG_OK) {
        puts ("  checking equation set");
        fflush (stdout);
        check_equation_set (eqset, z->dataclass, z->punits, 4);
    }

    /*----- Coordinates -----*/

    if (cg_ngrids (cgnsfn, cgnsbase, cgnszone,
        &NumGridCoordinate)) error_exit("cg_ngrids");
    if (!NumGridCoordinate)
        error ("no grid coordinates defined");
    else {
        create_names (NumGridCoordinate, &MaxGridCoordinate, &GridCoordinate);
        for (n = 1; n <= NumGridCoordinate; n++)
            check_coordinates (n);
    }

    /*----- Elements -----*/

    if (z->nsets) check_elements ();

    /*----- ZoneBC -----*/

    check_zoneBC ();

    /*----- ZoneGridConnectivity -----*/

    if (cg_nzconns(cgnsfn, cgnsbase, cgnszone, &NumZoneConn))
        error_exit("cg_nzconns");
    create_names (NumZoneConn, &MaxZoneConn, &ZoneConn);
    for (n = 1; n <= NumZoneConn; n++)
        check_connectivity (n);

    /*----- RotatingCoordinates -----*/

    go_absolute ("Zone_t", cgnszone, NULL);
    ierr = cg_rotating_read (vector, point);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_rotating_read");
    if (ierr == CG_OK) {
        puts ("  checking rotating coordinates");
        fflush (stdout);
        check_rotating (point, vector, z->dataclass, z->punits, 4);
    }

    /*----- ArbitraryGridMotion -----*/

    if (cg_n_arbitrary_motions (cgnsfn, cgnsbase, cgnszone,
        &NumArbitraryGrid)) error_exit("cg_n_arbitrary_motions");
    create_names (NumArbitraryGrid, &MaxArbitraryGrid, &ArbitraryGrid);
    for (n = 1; n <= NumArbitraryGrid; n++)
        check_arbitrary_motion (n);

    /*----- RigidGridMotion -----*/

    if (cg_n_rigid_motions (cgnsfn, cgnsbase, cgnszone, &NumRigidGrid))
        error_exit("cg_n_rigid_motions");
    create_names (NumRigidGrid, &MaxRigidGrid, &RigidGrid);
    for (n = 1; n <= NumRigidGrid; n++)
        check_rigid_motion (n);

    /*----- FlowSolution -----*/

    if (cg_nsols (cgnsfn, cgnsbase, cgnszone, &NumFlowSolution))
        error_exit("cg_nsols");
    create_names (NumFlowSolution, &MaxFlowSolution, &FlowSolution);
    for (n = 1; n <= NumFlowSolution; n++)
        check_solution (n);

    /*----- SubRegions -----*/

    if (cg_nsubregs (cgnsfn, cgnsbase, cgnszone, &NumZoneSubReg))
        error_exit("cg_nsubregs");
    create_names (NumZoneSubReg, &MaxZoneSubReg, &ZoneSubReg);
    for (n = 1; n <= NumZoneSubReg; n++)
        check_subreg (n);

    /*----- ConvergenceHistory -----*/

    go_absolute ("Zone_t", cgnszone, NULL);
    ierr = cg_convergence_read (&niter, &desc);
    if (ierr) {
        if (ierr != CG_NODE_NOT_FOUND) error_exit("cg_convergence_read");
        niter = 0;
    }
    else {
        puts ("  checking zone convergence history");
        fflush (stdout);
        go_relative ("ConvergenceHistory_t", 1, NULL);
        check_convergence (niter, desc, z->dataclass, z->punits, 4);
        if (desc != NULL) cg_free (desc);
    }

    /*----- ZoneIterativeData -----*/

    ierr = cg_ziter_read (cgnsfn, cgnsbase, cgnszone, name);
    /* prior to 2.3 returned ERROR instead of NODE_NOT_FOUND */
    if (ierr && ierr != CG_NODE_NOT_FOUND && LibraryVersion >= 2300)
        error_exit("cg_ziter_read");
    if (ierr == CG_OK) {
        printf ("  checking zone iterative data \"%s\"\n", name);
        fflush (stdout);
        if (BaseIter)
            error ("ZoneIterativeData requires BaseIterativeData");
        check_zone_iter ();
    }

    /*----- IntegralData -----*/

    go_absolute ("Zone_t", cgnszone, NULL);
    check_integral (z->dataclass, z->punits, 2);

    /*----- UserDefinedData -----*/

    go_absolute ("Zone_t", cgnszone, NULL);
    check_user_data (z->dataclass, z->punits, 2);

    /*----- DiscreteData -----*/

    if (cg_ndiscrete (cgnsfn, cgnsbase, cgnszone, &nd))
        error_exit("cg_ndiscrete");
    for (n = 1; n <= nd; n++)
        check_discrete (n);
}

/*-----------------------------------------------------------------------*/

static void check_axisymmetry (float *point, float *vector)
{
    char name[33];
    int n, na, ndim;
    cgsize_t dims[12];
    int dataclass, *punits, units[9];
    CGNS_ENUMT(DataType_t) datatype;

    if (verbose) {
        printf ("  Reference Point=[%g,%g]\n", point[0], point[1]);
        printf ("  Axis Vector=[%g,%g]\n", vector[0], vector[1]);
    }

    go_absolute ("Axisymmetry_t", 1, NULL);

    if (verbose > 1) {
        int nd;
        char *desc;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("  Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 2);
        if (punits) print_units (punits, 2);
    }
    if (dataclass < 0) dataclass = BaseClass;
    if (!punits) punits = pBaseUnits;

    if (cg_narrays (&na)) error_exit("cg_narrays");
    for (n = 1; n <= na; n++) {
        if (cg_array_info (n, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        printf ("  checking axisymmetry data \"%s\"\n", name);
        fflush (stdout);
        if (0 == strcmp (name, "AxisymmetryAxisVector") ||
            0 == strcmp (name, "AxisymmetryReferencePoint")) {
            if (ndim != 1 || dims[0] != 2)
                error ("bad dimension values");
            if (datatype != CGNS_ENUMV(RealSingle))
                error ("data type not real");
            check_quantity (n, name, dataclass, punits, 1, 4);
        }
        else if (0 == strcmp (name, "AxisymmetryAngle")) {
            if (ndim != 1 || dims[0] != 1)
                error ("bad dimension values");
            if (datatype != CGNS_ENUMV(RealSingle))
                error ("data type not real");
            check_quantity (n, name, dataclass, punits, 1, 4);
        }
        else if (0 == strcmp (name, "CoordinateNames")) {
            if (ndim != 2 || dims[0] != 32 || dims[1] != 2)
                error ("bad dimension values");
            if (datatype != CGNS_ENUMV(Character))
                error ("data type not character");
        }
        else
            warning (1, "not valid as child of Axisymmetry");
    }

    check_user_data (dataclass, punits, 2);
}

/*-----------------------------------------------------------------------*/

static void check_gravity (float *vector)
{
    char *desc, name[33];
    int n, na, nd, ndim;
    cgsize_t dims[12];
    int dataclass, *punits, units[9];
    CGNS_ENUMT(DataType_t) datatype;

    if (verbose) {
        printf ("  Vector=[%g]", vector[0]);
        for (nd = 1; nd < PhyDim; nd++)
            printf (",%g", vector[nd]);
        puts ("]");
    }

    go_absolute ("Gravity_t", 1, NULL);

    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("  Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 2);
        if (punits) print_units (punits, 2);
    }
    if (dataclass < 0) dataclass = BaseClass;
    if (!punits) punits = pBaseUnits;

    if (cg_narrays (&na)) error_exit("cg_narrays");
    for (n = 1; n <= na; n++) {
        if (cg_array_info (n, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        printf ("  checking gravity data \"%s\"\n", name);
        fflush (stdout);
        if (0 == strcmp (name, "GravityVector")) {
            if (ndim != 1 || dims[0] != PhyDim)
                error ("invalid dimension values");
            if (datatype != CGNS_ENUMV(RealSingle))
                error ("data type is not RealSingle");
            check_quantity (n, name, dataclass, punits, 1, 4);
        }
        else
            warning (1, "not valid as child of Gravity");
    }

    check_user_data (dataclass, punits, 2);
}

/*-----------------------------------------------------------------------*/

static void check_family (int fam)
{
    char famname[33], name[33], cad[33], *filename;
    int ierr, i, n, nbc, ngeo, nparts;
    CGNS_ENUMT(BCType_t) bctype;
    int nds, dirichlet, neumann;
    float point[3], vector[3];

    if (cg_family_read (cgnsfn, cgnsbase, fam, famname, &nbc, &ngeo))
        error_exit("cg_family_read");
    printf ("checking family \"%s\"\n", famname);

    go_absolute ("Family_t", fam, NULL);
    if (verbose > 1) {
        int nd;
        char *desc;
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("  Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    if (verbose)
        printf ("  Number BC=%d\n", nbc);
    if (nbc < 0 || nbc > 1)
        error ("number of BCs not 0 or 1 for %s", famname);
    for (n = 1; n <= nbc; n++) {
        if (cg_fambc_read (cgnsfn, cgnsbase, fam, n, name, &bctype))
            error_exit("cg_fambc_read");
        if (verbose) {
            printf ("    BC Name=\"%s\"\n", name);
            printf ("    BC Type=%s\n", cg_BCTypeName(bctype));
        }
        go_relative ("FamilyBC_t", n, NULL);
        ierr = cg_bcdataset_info (&nds);
        if (ierr && ierr != CG_NODE_NOT_FOUND)
            error_exit("cg_bcdataset_info");
        if (ierr == CG_OK && nds > 0) {
            for (i = 1; i <= nds; i++) {
                if (cg_bcdataset_read (i, name, &bctype, &dirichlet, &neumann))
                    error_exit("cg_bcdataset_read");
                printf ("  checking BC data set \"%s\"\n", name);
                fflush (stdout);
                if (FileVersion < 3100)
                    go_relative ("BCDataSet_t", i, NULL);
                else
                    go_relative ("FamilyBCDataSet_t", i, NULL);
                check_BCdata (bctype, dirichlet, neumann, 1,
                    BaseClass, pBaseUnits, 4);
                go_relative ("..", 1, NULL);
            }
        }
        go_relative ("..", 1, NULL);
    }

    if (verbose) printf ("  Number Geo=%d\n", ngeo);
    for (n = 1; n <= ngeo; n++) {
        if (cg_geo_read (cgnsfn, cgnsbase, fam, n, name, &filename,
                cad, &nparts)) error_exit("cg_geo_read");
        if (verbose) {
            printf ("    Geo Name=\"%s\"\n", name);
            printf ("    Geo File=\"%s\"\n",
                filename == NULL ? "<unknown>" : filename);
            printf ("    CAD=\"%s\"\n", cad);
            printf ("    Nparts=%d\n", nparts);
        }
        if (filename != NULL) cg_free (filename);
        for (i = 1; i <= nparts; i++) {
            if (cg_part_read (cgnsfn, cgnsbase, fam, n, i, name))
                error_exit("cg_part_read");
            if (verbose > 2)
                printf ("      Part %d Name=\"%s\"\n", i, name);
        }
    }

    ierr = read_ordinal (&i);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_ordinal_read");

    ierr = cg_rotating_read (vector, point);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_rotating_read");
    if (ierr == CG_OK) {
        puts ("  checking rotating coordinates");
        fflush (stdout);
        check_rotating (point, vector, BaseClass, pBaseUnits, 4);
    }

    check_user_data (BaseClass, pBaseUnits, 2);
}

/*-----------------------------------------------------------------------*/

static void check_base_iter (void)
{
    char *p, *desc, name[33];
    int ierr, n, na, ns, nd, nmax, ndim;
    cgsize_t dims[12];
    int dataclass, *punits, units[9];
    int *icnt, nnf = 0, nfp = 0, nnz = 0, nzp = 0;
    CGNS_ENUMT(DataType_t) datatype;

    if (verbose) printf ("  Number Steps=%d\n", NumSteps);
    if (NumSteps < 1)
        warning (2, "number of time steps is not > 0");

    go_absolute ("BaseIterativeData_t", 1, NULL);
    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("  Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    dataclass = read_dataclass ();
    punits = read_units (units);
    if (verbose) {
        if (dataclass >= 0) print_dataclass (dataclass, 2);
        if (punits) print_units (punits, 2);
    }
    if (dataclass < 0) dataclass = BaseClass;
    if (!punits) punits = pBaseUnits;

    if (cg_narrays (&na)) error_exit("cg_narrays");
    ns = 0;
    for (n = 1; n <= na; n++) {
        if (cg_array_info (n, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        if (0 == strcmp (name, "NumberOfFamilies")) {
            nnf = n;
            continue;
        }
        if (0 == strcmp (name, "FamilyPointers")) {
            nfp = n;
            continue;
        }
        if (0 == strcmp (name, "NumberOfZones")) {
            nnz = n;
            continue;
        }
        if (0 == strcmp (name, "ZonePointers")) {
            nzp = n;
            continue;
        }
        printf ("  checking \"%s\"\n", name);
        fflush (stdout);
        if (0 == strcmp (name, "IterationValues")) {
            if (ndim != 1 || dims[0] != NumSteps)
                error ("invalid dimension values");
            if (datatype != CGNS_ENUMV(Integer))
                error ("data type not integer");
            ns++;
        }
        else if (0 == strcmp (name, "TimeValues")) {
            if (ndim != 1 || dims[0] != NumSteps)
                error ("invalid dimension values");
            if (datatype != CGNS_ENUMV(RealSingle) &&
                datatype != CGNS_ENUMV(RealDouble))
                error ("data type not real or double");
            ns++;
        }
        else {
            check_quantity (n, name, dataclass, punits, 0, 4);
        }
    }
    if (!ns)
        error ("TimeValues and/or IterationValues is required");

    /* check family pointers */

    if (nnf && !nfp)
        error ("NumberOfFamilies given but not FamilyPointers");
    if (!nnf && nfp)
        error ("FamilyPointers given but not NumberOfFamilies");
    if (nnf && nfp) {
        ierr = nmax = 0;
        icnt = NULL;
        puts ("  checking \"NumberOfFamilies\"");
        fflush (stdout);
        if (cg_array_info (nnf, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        if (ndim != 1 || dims[0] != NumSteps) {
            error ("invalid dimension values");
            ierr = 1;
        }
        if (datatype != CGNS_ENUMV(Integer)) {
            error ("data type not integer");
            ierr = 1;
        }
        if (!ierr && NumSteps > 0) {
            icnt = (int *) malloc (NumSteps * sizeof(int));
            if (icnt == NULL)
                fatal_error("malloc failed for number of families\n");
            if (cg_array_read (nnf, icnt)) error_exit("cg_array_read");
            for (ns = 0; ns < NumSteps; ns++) {
                if (icnt[ns] < 0 || icnt[ns] > NumFamily) ierr++;
                if (nmax < icnt[ns]) nmax = icnt[ns];
            }
            if (ierr)
                error ("there are %d invalid entries", ierr);
            if (nmax == 0)
                error ("max number of families is 0");
        }

        puts ("  checking \"FamilyPointers\"");
        fflush (stdout);
        if (cg_array_info (nfp, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        if (ndim != 3 || dims[0] != 32 ||
            dims[1] != nmax || dims[2] != NumSteps) {
            error ("invalid dimension values");
            ierr = 1;
        }
        if (datatype != CGNS_ENUMV(Character)) {
            error ("data type not character");
            ierr = 1;
        }
        if (NumSteps > 0 && nmax > 0 && !ierr) {
            desc = (char *) malloc (32 * nmax * NumSteps * sizeof(char));
            if (NULL == desc)
                fatal_error("malloc failed for family pointers\n");
            if (cg_array_read (nfp, desc)) error_exit("cg_array_read");
            for (ierr = 0, n = 0, ns = 0; ns < NumSteps; ns++) {
                for (nd = 0; nd < icnt[ns]; nd++) {
                    strncpy (name, &desc[n + 32 * nd], 32);
                    name[32] = 0;
                    p = name + strlen(name);
                    while (--p >= name && isspace(*p))
                        ;
                    *++p = 0;
                    for (na = 0; na < NumFamily; na++) {
                        if (0 == strcmp (name, Family[na])) break;
                    }
                    if (na >= NumFamily) ierr++;
                }
                n += 32 * nmax;
            }
            free (desc);
            if (ierr)
                warning (1, "%d unknown families", ierr);
        }
        if (icnt != NULL) free (icnt);
    }

    /* check zone pointers */

    if (nnz && !nzp)
        error ("NumberOfZones given but not ZonePointers");
    if (!nnz && nzp)
        error ("ZonePointers given but not NumberOfZones");
    if (nnz && nzp) {
        ierr = nmax = 0;
        icnt = NULL;
        puts ("  checking \"NumberOfZones\"");
        fflush (stdout);
        if (cg_array_info (nnz, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        if (ndim != 1 || dims[0] != NumSteps) {
            error ("invalid dimension values");
            ierr = 1;
        }
        if (datatype != CGNS_ENUMV(Integer)) {
            error ("data type not integer");
            ierr = 1;
        }
        if (!ierr && NumSteps > 0) {
            icnt = (int *) malloc (NumSteps * sizeof(int));
            if (icnt == NULL)
                fatal_error("malloc failed for number of zones\n");
            if (cg_array_read (nnz, icnt)) error_exit("cg_array_read");
            for (ns = 0; ns < NumSteps; ns++) {
                if (icnt[ns] < 0 || icnt[ns] > NumZones) ierr++;
                if (nmax < icnt[ns]) nmax = icnt[ns];
            }
            if (ierr)
                error ("there are %d invalid entries", ierr);
            if (nmax == 0)
                error ("max number of zones is 0");
        }

        printf ("  checking \"ZonePointers\"\n");
        fflush (stdout);
        if (cg_array_info (nzp, name, &datatype, &ndim, dims))
            error_exit("cg_array_info");
        if (ndim != 3 || dims[0] != 32 ||
            dims[1] != nmax || dims[2] != NumSteps) {
            error ("invalid dimension values");
            ierr = 1;
        }
        if (datatype != CGNS_ENUMV(Character)) {
            error ("data type not character");
            ierr = 1;
        }
        if (NumSteps > 0 && nmax > 0 && !ierr) {
            desc = (char *) malloc (32 * nmax * NumSteps * sizeof(char));
            if (NULL == desc)
                fatal_error("malloc failed for zone pointers\n");
            if (cg_array_read (nzp, desc)) error_exit("cg_array_read");
            for (ierr = 0, n = 0, ns = 0; ns < NumSteps; ns++) {
                for (nd = 0; nd < icnt[ns]; nd++) {
                    strncpy (name, &desc[n + 32 * nd], 32);
                    name[32] = 0;
                    p = name + strlen(name);
                    while (--p >= name && isspace(*p))
                        ;
                    *++p = 0;
                    for (na = 0; na < NumZones; na++) {
                        if (0 == strcmp (name, Zones[na].name)) break;
                    }
                    if (na >= NumZones) ierr++;
                }
                n += 32 * nmax;
            }
            free (desc);
            if (ierr)
                warning (1, "%d unknown zones", ierr);
        }
        if (icnt != NULL) free (icnt);
    }

    check_user_data (dataclass, punits, 2);
}

/*=======================================================================*/

static void check_base (void)
{
    char basename[33], name[33], *desc;
    int n, nz, ierr, nd, nf, eqset[7];
    float point[3], vector[3];
    CGNS_ENUMT(SimulationType_t) simulation;

    /*----- base dimensions -----*/

    if (cg_base_read (cgnsfn, cgnsbase, basename, &CellDim, &PhyDim))
        error_exit("cg_base_read");
    printf ("\nreading base \"%s\"\n", basename);
    fflush (stdout);
    if (CellDim < 1) {
        error ("Cell Dimension < 1");
        return;
    }
    if (PhyDim < CellDim) {
        error ("Physical Dimension < Cell Dimension");
        return;
    }

    if (CellDim < 2 || CellDim > 3 || PhyDim < 2 || PhyDim > 3) {
        puts ("INTERNAL:can't handle CellDim and/or Phydim < 2 or > 3");
        return;
    }

    /*----- read zones -----*/

    for (nz = 0; nz < NumZones; nz++) {
        if (Zones[nz].nsets) {
            for (n = 0; n < Zones[nz].nsets; n++) {
                if (Zones[nz].sets[n].elements != NULL)
                    free (Zones[nz].sets[n].elements);
                if (Zones[nz].sets[n].parent != NULL)
                    free (Zones[nz].sets[n].parent);
            }
            free (Zones[nz].sets);
        }
        if (Zones[nz].faces)
            HashDestroy (Zones[nz].faces, free);
        if (Zones[nz].nextnodes)
            free (Zones[nz].extnodes);
    }

    if (cg_nzones (cgnsfn, cgnsbase, &NumZones)) error_exit("cg_nzones");
    if (NumZones > MaxZones) {
        if (MaxZones)
            Zones = (ZONE *) realloc (Zones, NumZones * sizeof(ZONE));
        else
            Zones = (ZONE *) malloc (NumZones * sizeof(ZONE));
        if (NULL == Zones)
            fatal_error("malloc failed for zones\n");
        MaxZones = NumZones;
    }

    for (nz = 0; nz < NumZones; nz++)
        read_zone (nz);
    cgnszone = 0;

    /*----- read families -----*/

    if (cg_nfamilies (cgnsfn, cgnsbase, &NumFamily))
        error_exit("cg_nfamilies");
    if (NumFamily) {
        puts ("reading families");
        fflush (stdout);
        create_names (NumFamily, &MaxFamily, &Family);
        for (nf = 0; nf < NumFamily; nf++) {
            if (cg_family_read (cgnsfn, cgnsbase, nf + 1, name, &n, &nd))
                error_exit("cg_family_read");
            strcpy (Family[nf], name);
        }
    }

    /*----- check base -----*/

    printf ("\nchecking base \"%s\"\n", basename);
    if (verbose)
        printf ("  Cell Dimension=%d\n  Physical Dimension=%d\n",
            CellDim, PhyDim);
    fflush (stdout);

    go_absolute (NULL);

    /*----- base descriptors -----*/

    if (verbose > 1) {
        if (cg_ndescriptors (&nd)) error_exit("cg_ndescriptors");
        for (n = 1; n <= nd; n++) {
            if (cg_descriptor_read (n, name, &desc))
                error_exit("cg_descriptor_read");
            if (desc != NULL) {
                printf ("  Descriptor %s:\n%s\n", name, desc);
                cg_free (desc);
            }
        }
    }

    /*----- SimulationType -----*/

    ierr = cg_simulation_type_read (cgnsfn, cgnsbase, &simulation);
    if (ierr && ierr != CG_NODE_NOT_FOUND)
        error_exit("cg_simulation_type_read");
    if (ierr == CG_OK && verbose)
        printf ("  Simulation Type=%s\n", cg_SimulationTypeName(simulation));

    /*----- data class, dimensional units -----*/

    BaseClass = read_dataclass ();
    pBaseUnits = read_units (BaseUnits);
    if (verbose) {
        if (BaseClass >= 0) print_dataclass (BaseClass, 2);
        if (pBaseUnits) print_units (pBaseUnits, 2);
    }

    /*----- ReferenceState -----*/

    ierr = cg_state_read (&desc);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_state_read");
    if (ierr == CG_OK) {
        puts ("checking reference state");
        if (desc != NULL) {
            if (verbose > 1)
                printf ("  Descriptor:%s\n", desc);
            cg_free (desc);
        }
        fflush (stdout);
        go_absolute ("ReferenceState_t", 1, NULL);
        check_arrays (BaseClass, pBaseUnits, 1, 0, 2);
    }

    /*----- FlowEquationSet -----*/

    go_absolute (NULL);
    ierr = cg_equationset_read (&eqset[0], &eqset[1], &eqset[2],
        &eqset[3], &eqset[4], &eqset[5], &eqset[6]);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_equationset_read");
    if (ierr == CG_OK) {
        puts ("checking equation set");
        fflush (stdout);
        check_equation_set (eqset, BaseClass, pBaseUnits, 2);
    }

    /*----- Families -----*/

    for (nf = 1; nf <= NumFamily; nf++)
        check_family (nf);

    /*----- Axisymmetry -----*/

    ierr = cg_axisym_read (cgnsfn, cgnsbase, point, vector);
    if (ierr && ierr != CG_NODE_NOT_FOUND) error_exit("cg_axisym_read");
    if (ierr == CG_OK) {
        if (PhyDim != 2)
            error ("axisymmetry is only valid for physical dimension of 2");
        else {
            puts ("checking axisymmetry");
            fflush (stdout);
            check_axisymmetry (point, vector);
        }
    }

    /*----- RotatingCoordinates -----*/

    go_absolute (NULL);
    ierr = cg_rotating_read (vector, point);
    if (ierr && ierr != CG_NODE_NOT_FOUND)
        error_exit("cg_rotating_read");
    if (ierr == CG_OK) {
        puts ("checking rotating coordinates");
        fflush (stdout);
        check_rotating (point, vector, BaseClass, pBaseUnits, 2);
    }

    /*----- Gravity -----*/

    ierr = cg_gravity_read (cgnsfn, cgnsbase, vector);
    if (ierr && ierr != CG_NODE_NOT_FOUND)
        error_exit("cg_gravity_read");
    if (ierr == CG_OK) {
        puts ("checking gravity");
        fflush (stdout);
        check_gravity (vector);
    }

    /*----- BaseIterativeData -----*/

    BaseIter = cg_biter_read (cgnsfn, cgnsbase, name, &NumSteps);
    if (BaseIter && BaseIter != CG_NODE_NOT_FOUND)
        error_exit("cg_biter_read");
    if (BaseIter == CG_OK) {
        printf ("checking base iterative data \"%s\"\n", name);
        fflush (stdout);
        check_base_iter ();
    }
    else
        NumSteps = 0;

    /*----- ConvergenceHistory -----*/

    go_absolute (NULL);
    ierr = cg_convergence_read (&nd, &desc);
    if (ierr && ierr != CG_NODE_NOT_FOUND)
        error_exit("cg_convergence_read");
    if (ierr == CG_OK && nd) {
        puts ("checking global convergence history");
        fflush (stdout);
        go_absolute ("ConvergenceHistory_t", 1, NULL);
        check_convergence (nd, desc, BaseClass, pBaseUnits, 2);
        if (desc != NULL) cg_free (desc);
    }

    /*=----- IntegralData -----*/

    go_absolute (NULL);
    check_integral (BaseClass, pBaseUnits, 0);

    /*----- UserDefinedData -----*/

    go_absolute (NULL);
    check_user_data (BaseClass, pBaseUnits, 0);

    /*----- Zones -----*/

    if (NumZones == 0) {
        warning (1, "no zones defined");
        return;
    }

    for (cgnszone = 1; cgnszone <= NumZones; cgnszone++)
        check_zone ();
}

/*=======================================================================*/

int main (int argc, char *argv[])
{
    char *cgnsfile;
    int n, nbases, update = 0;
    float file_version;

    if (argc < 2)
        print_usage (usgmsg, NULL);

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'v':
                verbose = 1;
                break;
            case 'V':
                verbose = 2;
                break;
            case 'u':
            case 'U':
                update = n;
                break;
            case 'w':
                dowarn = atoi (argarg);
                break;
            case 'e':
                doerr = 0;
                break;
        }
    }

    if (argind == argc)
        print_usage (usgmsg, "CGNSfile not given");
    cgnsfile = argv[argind++];

    /* update CGNS file by opening it in modify mode */

    if (update) {
        char *newfile = argind < argc ? argv[argind] : NULL;
        cgnsfile = update_version (cgnsfile, newfile);
        if (update == 'U')
            exit (0);
    }

    printf ("reading CGNS file %s\n", cgnsfile);
    fflush (stdout);
    if (cg_open (cgnsfile, CG_MODE_READ, &cgnsfn)) error_exit("cg_open");

    /* get version */

    if (cg_version (cgnsfn, &file_version)) error_exit("cg_version");
    FileVersion = (int)(1000.0 * file_version + 0.5);
    if (LibraryVersion < FileVersion)
        warning (1, "CGNS file version is more recent than library version");
    if (verbose) {
        int precision;
        printf ("File Version=%g\n", file_version);
        cg_precision (cgnsfn, &precision);
        if (precision) printf ("File Precision=%d\n", precision);
    }

    /* get number of bases */

    if (cg_nbases (cgnsfn, &nbases)) error_exit("cg_nbases");
    if (nbases < 1) warning (1, "no bases defined in CGNS file");
    for (cgnsbase = 1; cgnsbase <= nbases; cgnsbase++)
        check_base ();

    /* close CGNS file and exit */

    if (cg_close (cgnsfn)) error_exit("cg_close");
    puts ("\nchecking complete");
    if (totwarn) printf ("%d warnings (%d shown)\n", totwarn, nwarn);
    if (nerr) printf ("%d errors\n", nerr);
    return 0;
}

