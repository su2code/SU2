/*
 * cgns_to_vtk - read CGNS file and write VTK file
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
# include <io.h>
# include <direct.h>
# define chdir _chdir
#else
# include <unistd.h>
#endif

#include "getargs.h"
#include "cgnslib.h"

#ifndef CG_MODE_READ
# define CG_MODE_READ MODE_READ
#endif

#if defined(_WIN32) || defined(__CYGWIN__) || defined(__linux) || \
    defined(__alpha) || defined(__ultrix)
# ifndef BYTE_SWAPPED
#  define BYTE_SWAPPED
# endif
#endif

static char options[] = "b:z:s:vea";

static char *usgmsg[] = {
    "usage  : cgns_to_vtk [options] CGNSfile [outdir]",
    "options:",
    "   -b<base> = base number (default 1)",
    "   -z<zone> = zone number (default 0 - all)",
    "   -s<soln> = solution number (default 1)",
    "   -v       = verbose output",
    "   -e       = write element sets (unstructured zone)",
    "   -a       = write ascii instead of binary format",
    "<outdir> is the output directory for the VTK files.",
    "If not specified, it defaults to the current directory.",
    NULL
};

typedef float Node[3];
typedef struct {
    int cnt;
    char name[33];
} Variable;

static int nzones;
static int cgnsfn;
static int cgnsbase = 1;
static int cgnszone = 0;
static int cgnssol = 1;
static int CellDim, PhyDim;

static int nnodes;
static Node *nodes;

static CGNS_ENUMT(GridLocation_t) varloc;
static int nvars, ndata;
static Variable *vars;
static cgsize_t varrng[2][3];

static int verbose = 0;
static int ascii = 0;

/*---------- FATAL ----------------------------------------------------
 * exit with error message
 *---------------------------------------------------------------------*/

static void FATAL (char *errmsg)
{
    if (NULL == errmsg)
        fprintf (stderr, "CGNS error:%s\n", cg_get_error());
    else
        fprintf (stderr, "%s\n", errmsg);
    exit (1);
}

/*---------- create_filename -----------------------------------------
 * create valid filename
 *--------------------------------------------------------------------*/

static void create_filename (char *str, char *fname)
{
    int n = 0;
    char *p;

    for (p = str; *p; p++) {
#ifdef _WIN32
        if (strchr ("\\/:*?\"<>|", *p) == NULL)
#else
        if (isspace(*p)) continue;
        if (strchr ("\\/:*?\"<>|[]()", *p) == NULL)
#endif
            fname[n++] = *p;
        else
            fname[n++] = '_';
    }
    fname[n] = 0;
}

/*---------- fix_name ------------------------------------------------
 * remove invalid characters from variable name
 *--------------------------------------------------------------------*/

static void fix_name (char *str, char *name)
{
    int n = 0;
    char *p;

    for (p = str; *p; p++) {
        if (!isspace(*p))
            name[n++] = *p;
    }
    name[n] = 0;
}

/*---------- swap_bytes ----------------------------------------------
 * swap bytes
 *--------------------------------------------------------------------*/

#ifdef BYTE_SWAPPED

static void *swap_bytes (int bytes, void *data)
{
    static unsigned char buf[sizeof(double)];
    unsigned char *p = (unsigned char *)data;
    int n;

    p += bytes;
    for (n = 0; n < bytes; n++)
        buf[n] = *--p;
    return ((void *)buf);
}

#endif

/*---------- write_ints ----------------------------------------------
 * write integers to VTK file
 *--------------------------------------------------------------------*/

static void write_ints (FILE *fp, int cnt, int *data)
{
    if (ascii) {
        fprintf (fp, "%d", *data);
        while (--cnt > 0) {
            data++;
            fprintf (fp, " %d", *data);
        }
        putc ('\n', fp);
    }
    else {
#ifdef BYTE_SWAPPED
        while (cnt-- > 0) {
            fwrite (swap_bytes (sizeof(int), (void *)data),
                sizeof(int), 1, fp);
            data++;
        }
#else
        fwrite (data, sizeof(int), cnt, fp);
#endif
    }
}

/*---------- write_floats --------------------------------------------
 * write floats to VTK file
 *--------------------------------------------------------------------*/

static void write_floats (FILE *fp, int cnt, float *data)
{
    if (ascii) {
        fprintf (fp, "%g", *data);
        while (--cnt > 0) {
            data++;
            fprintf (fp, " %g", *data);
        }
        putc ('\n', fp);
    }
    else {
#ifdef BYTE_SWAPPED
        while (cnt-- > 0) {
            fwrite (swap_bytes (sizeof(float), (void *)data),
                sizeof(float), 1, fp);
            data++;
        }
#else
        fwrite (data, sizeof(float), cnt, fp);
#endif
    }
}

/*---------- get_nodes ------------------------------------------------
 * read zone nodes
 *---------------------------------------------------------------------*/

static int get_nodes (int nz, CGNS_ENUMT(ZoneType_t) zonetype, cgsize_t *sizes)
{
    int i, j, n, ncoords;
    int rind[6];
    cgsize_t nn, rng[2][3];
    CGNS_ENUMT(DataType_t) datatype;
    float *xyz;
    double rad, theta, phi;
    char name[33], coordtype[4];

    /* get number of coordinates */

    if (cg_ncoords (cgnsfn, cgnsbase, nz, &ncoords))
        FATAL (NULL);
    if (ncoords < PhyDim)
        FATAL ("less than PhyDim coordinates");

    /* check for rind */

    if (cg_goto (cgnsfn, cgnsbase, "Zone_t", nz,
        "GridCoordinates_t", 1, "end"))
        FATAL (NULL);
    if ((i = cg_rind_read (rind)) != CG_OK) {
        if (i != CG_NODE_NOT_FOUND)
            FATAL (NULL);
        for (n = 0; n < 6; n++)
            rind[n] = 0;
    }

    /* get grid coordinate range */

    if (zonetype == CGNS_ENUMV(Structured)) {
        for (n = 0; n < 3; n++) {
            rng[0][n] = 1;
            rng[1][n] = 1;
        }
        nn = 1;
        for (n = 0; n < CellDim; n++) {
            rng[0][n] = rind[2*n] + 1;
            rng[1][n] = rind[2*n] + sizes[n];
            nn *= sizes[n];
        }
    }
    else {
        nn = sizes[0] + rind[0] + rind[1];
        rng[0][0] = 1;
        rng[1][0] = nn;
    }
    nnodes = (int)nn;

    /* read the nodes */

    strcpy (coordtype, "   ");
    xyz = (float *) malloc (nnodes * sizeof(float));
    nodes = (Node *) malloc (nnodes * sizeof(Node));
    if (xyz == NULL || nodes == NULL)
        FATAL ("malloc failed for nodes");

    for (i = 1; i <= ncoords; i++) {
        if (cg_coord_info (cgnsfn, cgnsbase, nz, i, &datatype, name) ||
            cg_coord_read (cgnsfn, cgnsbase, nz, name,
                CGNS_ENUMV(RealSingle), rng[0], rng[1], xyz))
            FATAL (NULL);
        if (0 == strcmp (name, "CoordinateX") ||
            0 == strcmp (name, "CoordinateR"))
            j = 0;
        else if (0 == strcmp (name, "CoordinateY") ||
                 0 == strcmp (name, "CoordinateTheta"))
            j = 1;
        else if (0 == strcmp (name, "CoordinateZ") ||
                 0 == strcmp (name, "CoordinatePhi"))
            j = 2;
        else
            continue;
        if (coordtype[j] == ' ' || strchr ("XYZ", name[10]) != NULL)
            coordtype[j] = name[10];
        for (n = 0; n < nnodes; n++)
            nodes[n][j] = xyz[n];
    }
    free (xyz);

    /* change coordinate system to cartesian */

    if (0 == strncmp (coordtype, "RTZ", PhyDim)) {
        for (n = 0; n < nnodes; n++) {
            rad   = nodes[n][0];
            theta = nodes[n][1];
            nodes[n][0] = (float)(rad * cos (theta));
            nodes[n][1] = (float)(rad * sin (theta));
        }
    }
    else if (0 == strcmp (coordtype, "RTP")) {
        for (n = 0; n < nnodes; n++) {
            rad   = nodes[n][0];
            theta = nodes[n][1];
            phi   = nodes[n][2];
            nodes[n][0] = (float)(rad * sin (theta) * cos (phi));
            nodes[n][1] = (float)(rad * sin (theta) * sin (phi));
            nodes[n][2] = (float)(rad * cos (theta));
        }
    }
    else if (strncmp (coordtype, "XYZ", PhyDim))
        FATAL ("unknown coordinate system");

    return nnodes;
}

/*---------- sort_variables -----------------------------------------
 * sort variables by name
 *-------------------------------------------------------------------*/

static int sort_variables (const void *v1, const void *v2)
{
    Variable *var1 = (Variable *)v1;
    Variable *var2 = (Variable *)v2;

    return (strcmp (var1->name, var2->name));
}

/*---------- get_variables --------------------------------------------
 * get the solution vaiables
 *---------------------------------------------------------------------*/

static int get_variables (int nz, CGNS_ENUMT(ZoneType_t) zonetype, cgsize_t *sizes)
{
    char name[33];
    int n, len, nv, nsols;
    int rind[6];
    CGNS_ENUMT(DataType_t) datatype;

    nvars = 0;
    if (cg_nsols (cgnsfn, cgnsbase, nz, &nsols))
        FATAL (NULL);
    if (cgnssol < 1 || cgnssol > nsols) return 0;
    if (cg_sol_info (cgnsfn, cgnsbase, nz, cgnssol, name, &varloc) ||
        cg_nfields (cgnsfn, cgnsbase, nz, cgnssol, &nv))
        FATAL (NULL);
    if (nv < 1) return 0;
    if (varloc != CGNS_ENUMV(Vertex) && varloc != CGNS_ENUMV(CellCenter)) return 0;
    nvars = nv;

    /* check for rind */

    if (cg_goto (cgnsfn, cgnsbase, "Zone_t", nz,
        "FlowSolution_t", cgnssol, "end"))
        FATAL (NULL);
    if ((n = cg_rind_read (rind)) != CG_OK) {
        if (n != CG_NODE_NOT_FOUND)
            FATAL (NULL);
        for (n = 0; n < 6; n++)
            rind[n] = 0;
    }

    /* get solution data range */

    if (zonetype == CGNS_ENUMV(Structured)) {
        nv = varloc == CGNS_ENUMV(Vertex) ? 0 : CellDim;
        for (n = 0; n < 3; n++) {
            varrng[0][n] = 1;
            varrng[1][n] = 1;
        }
        ndata = 1;
        for (n = 0; n < CellDim; n++) {
            varrng[0][n] = rind[2*n] + 1;
            varrng[1][n] = rind[2*n] + sizes[n+nv];
            ndata *= (int)sizes[n+nv];
        }
    }
    else {
        nv = varloc == CGNS_ENUMV(Vertex) ? 0 : 1;
        ndata = (int)sizes[nv];
        varrng[0][0] = rind[0] + 1;
        varrng[1][0] = rind[0] + ndata;
    }

    /* get variable names */

    vars = (Variable *) malloc (nvars * sizeof(Variable));
    if (vars == NULL)
        FATAL ("malloc failed for variable names");

    for (nv = 0; nv < nvars; nv++) {
        if (cg_field_info (cgnsfn, cgnsbase, nz, cgnssol, nv+1,
                &datatype, name))
            FATAL (NULL);
        vars[nv].cnt = 1;
        strcpy (vars[nv].name, name);
    }
    qsort (vars, nvars, sizeof(Variable), sort_variables);

    /* get number of scalars and vectors */

    for (nv = 2; nv < nvars; nv++) {
        len = (int)strlen(vars[nv].name) - 1;
        if (vars[nv].name[len] == 'Z') {
            strcpy (name, vars[nv].name);
            name[len] = 'Y';
            if (0 == strcmp (name, vars[nv-1].name)) {
                name[len] = 'X';
                if (0 == strcmp (name, vars[nv-2].name)) {
                    vars[nv-2].cnt = 3;
                    vars[nv-1].cnt = 0;
                    vars[nv].cnt = 0;
                }
            }
        }
    }

    return nvars;
}

/*---------- write_solution -------------------------------------------
 * write solution to vtk file
 *---------------------------------------------------------------------*/

static void write_solution (FILE *fp, int nz, int *mask)
{
    char name[33];
    int n, nv, nd;
    int nscal, nvect;
    float *data, v[3];

    nscal = nvect = 0;
    for (nv = 0; nv < nvars; nv++) {
        if (vars[nv].cnt == 1)
            nscal++;
        else if (vars[nv].cnt == 3)
            nvect++;
    }

    if (nvect)
        data = (float *) malloc (3 * ndata * sizeof(float));
    else
        data = (float *) malloc (ndata * sizeof(float));
    if (data == NULL)
        FATAL ("malloc failed for solution data");

    if (mask == NULL)
        nd = ndata;
    else {
        for (nd = 0, n = 0; n < ndata; n++) {
            if (mask[n]) nd++;
        }
    }

    if (verbose) {
        printf ("  writing %d scalars and %d vectors as %s data\n",
            nscal, nvect, varloc == CGNS_ENUMV(Vertex) ? "point" : "cell");
        fflush (stdout);
    }
    fprintf (fp, "%s_DATA %d\n", varloc == CGNS_ENUMV(Vertex) ? "POINT" : "CELL", nd);

    if (nscal) {
        for (nv = 0; nv < nvars; nv++) {
            if (vars[nv].cnt != 1) continue;
            if (cg_field_read (cgnsfn, cgnsbase, nz, cgnssol,
                    vars[nv].name, CGNS_ENUMV(RealSingle),
                    varrng[0], varrng[1], data))
                FATAL (NULL);
            fix_name (vars[nv].name, name);
            fprintf (fp, "SCALARS %s float\n", name);
            fprintf (fp, "LOOKUP_TABLE default\n");
            for (n = 0; n < ndata; n++) {
                if (mask == NULL || mask[n])
                    write_floats (fp, 1, &data[n]);
            }
        }
    }

    if (nvect) {
        for (nv = 0; nv < nvars; nv++) {
            if (vars[nv].cnt != 3) continue;
            if (cg_field_read (cgnsfn, cgnsbase, nz, cgnssol,
                    vars[nv].name, CGNS_ENUMV(RealSingle),
                    varrng[0], varrng[1], data) ||
                cg_field_read (cgnsfn, cgnsbase, nz, cgnssol,
                    vars[nv+1].name, CGNS_ENUMV(RealSingle),
                    varrng[0], varrng[1], &data[ndata]) ||
                cg_field_read (cgnsfn, cgnsbase, nz, cgnssol,
                    vars[nv+2].name, CGNS_ENUMV(RealSingle),
                    varrng[0], varrng[1], &data[2*ndata]))
                FATAL (NULL);
            fix_name (vars[nv].name, name);
            name[strlen(name)-1] = 0;
            fprintf (fp, "VECTORS %s float\n", name);
            for (n = 0; n < ndata; n++) {
                if (mask == NULL || mask[n]) {
                    v[0] = data[n];
                    v[1] = data[n+ndata];
                    v[2] = data[n+2*ndata];
                    write_floats (fp, 3, v);
                }
            }
        }
    }

    free (data);
}

/*---------- write_volume_cells ---------------------------------------
 * write volume cell data to vtk file
 *---------------------------------------------------------------------*/

static void write_volume_cells (FILE *fp, int nz)
{
    int i, n, ns, nsect, nn, ip;
    int elemcnt, elemsize;
    int *types, cell[9];
    cgsize_t is, ie, nelems, maxsize, maxelems;
    cgsize_t size, *conn;
    CGNS_ENUMT(ElementType_t) elemtype, et;
    char name[33];

    if (cg_nsections (cgnsfn, cgnsbase, nz, &nsect))
        FATAL (NULL);
    if (nsect < 1) FATAL ("no sections defined");

    maxsize = maxelems = 0;
    for (ns = 1; ns <= nsect; ns++) {
        if (cg_section_read (cgnsfn, cgnsbase, nz, ns,
                name, &elemtype, &is, &ie, &nn, &ip) ||
            cg_ElementDataSize (cgnsfn, cgnsbase, nz, ns, &size))
            FATAL (NULL);
        nelems = ie - is + 1;
        if (maxelems < nelems) maxelems = nelems;
        if (maxsize < size) maxsize = size;
    }
    if (maxsize > CG_MAX_INT32) FATAL("too many elements for 32-bit integer");

    conn = (cgsize_t *) malloc ((size_t)maxsize * sizeof(cgsize_t));
    if (conn == NULL)
        FATAL ("malloc failed for element connectivity");

    /* count volume cells */

    elemcnt = elemsize = 0;
    for (ns = 1; ns <= nsect; ns++) {
        if (cg_section_read (cgnsfn, cgnsbase, nz, ns,
                name, &elemtype, &is, &ie, &nn, &ip))
            FATAL (NULL);
        if (elemtype < CGNS_ENUMV(TETRA_4) || elemtype > CGNS_ENUMV(MIXED)) continue;
        nelems = ie - is + 1;
        if (elemtype == CGNS_ENUMV(MIXED)) {
            if (cg_elements_read (cgnsfn, cgnsbase, nz, ns, conn, NULL))
                FATAL (NULL);
            for (i = 0, n = 0; n < nelems; n++) {
                et = (int)conn[i++];
                switch (et) {
                    case CGNS_ENUMV(TETRA_4):
                    case CGNS_ENUMV(TETRA_10):
                        elemcnt++;
                        elemsize += 5;
                        break;
                    case CGNS_ENUMV(PYRA_5):
                    case CGNS_ENUMV(PYRA_14):
                        elemcnt++;
                        elemsize += 6;
                        break;
                    case CGNS_ENUMV(PENTA_6):
                    case CGNS_ENUMV(PENTA_15):
                    case CGNS_ENUMV(PENTA_18):
                        elemcnt++;
                        elemsize += 7;
                        break;
                    case CGNS_ENUMV(HEXA_8):
                    case CGNS_ENUMV(HEXA_20):
                    case CGNS_ENUMV(HEXA_27):
                        elemcnt++;
                        elemsize += 9;
                        break;
                    default:
                        break;
                }
                if (cg_npe (et, &nn) || nn == 0)
                    FATAL ("invalid element type in MIXED");
                i += nn;
            }
        }
        else {
            switch (elemtype) {
                case CGNS_ENUMV(TETRA_4):
                case CGNS_ENUMV(TETRA_10):
                    nn = 5;
                    break;
                case CGNS_ENUMV(PYRA_5):
                case CGNS_ENUMV(PYRA_14):
                    nn = 6;
                    break;
                case CGNS_ENUMV(PENTA_6):
                case CGNS_ENUMV(PENTA_15):
                case CGNS_ENUMV(PENTA_18):
                    nn = 7;
                    break;
                case CGNS_ENUMV(HEXA_8):
                case CGNS_ENUMV(HEXA_20):
                case CGNS_ENUMV(HEXA_27):
                    nn = 9;
                    break;
                default:
                    nn = 0;
                    break;
            }
            if (nn) {
                elemcnt += (int)nelems;
                elemsize += (nn * (int)nelems);
            }
        }
    }

    if (elemcnt == 0) {
        free (conn);
        return;
    }

    types = (int *) malloc (elemcnt * sizeof(int));
    if (types == NULL)
        FATAL ("malloc failed for cell types");

    /* write the elements */

    if (verbose) {
        printf ("  writing %d cells\n", elemcnt);
        fflush (stdout);
    }
    fprintf (fp, "CELLS %d %d\n", elemcnt, elemsize);

    elemcnt = 0;
    for (ns = 1; ns <= nsect; ns++) {
        if (cg_section_read (cgnsfn, cgnsbase, nz, ns,
                name, &elemtype, &is, &ie, &nn, &ip))
            FATAL (NULL);
        if (elemtype < CGNS_ENUMV(TETRA_4) || elemtype > CGNS_ENUMV(MIXED)) continue;
        nelems = ie - is + 1;
        if (cg_elements_read (cgnsfn, cgnsbase, nz, ns, conn, NULL))
            FATAL (NULL);
        et = elemtype;
        for (i = 0, n = 0; n < nelems; n++) {
            if (elemtype == CGNS_ENUMV(MIXED)) et = (int)conn[i++];
            switch (et) {
                case CGNS_ENUMV(TETRA_4):
                case CGNS_ENUMV(TETRA_10):
                    nn = 4;
                    types[elemcnt++] = 10;
                    break;
                case CGNS_ENUMV(PYRA_5):
                case CGNS_ENUMV(PYRA_14):
                    nn = 5;
                    types[elemcnt++] = 14;
                    break;
                case CGNS_ENUMV(PENTA_6):
                case CGNS_ENUMV(PENTA_15):
                case CGNS_ENUMV(PENTA_18):
                    nn = 6;
                    types[elemcnt++] = 13;
                    break;
                case CGNS_ENUMV(HEXA_8):
                case CGNS_ENUMV(HEXA_20):
                case CGNS_ENUMV(HEXA_27):
                    nn = 8;
                    types[elemcnt++] = 12;
                    break;
                default:
                    nn = 0;
                    break;
            }
            if (nn) {
                cell[0] = nn;
                for (ip = 0; ip < nn; ip++)
                    cell[ip+1] = (int)conn[i+ip] - 1;
                write_ints (fp, nn+1, cell);
            }
            if (cg_npe (et, &nn) || nn == 0)
                FATAL ("invalid element type");
            i += nn;
        }
    }
    free (conn);

    /* write the element types */

    fprintf (fp, "CELL_TYPES %d\n", elemcnt);
    write_ints (fp, elemcnt, types);

    free (types);
}

/*---------- write_element_sets ---------------------------------------
 * write element sets as vtk files
 *---------------------------------------------------------------------*/

static void write_element_sets (int nz, cgsize_t *sizes)
{
    int i, n, ns, nsect, nn, ip;
    int elemcnt, elemsize, cell[9];
    int *nodemap, *types;
    cgsize_t is, ie, nelems;
    cgsize_t size, *conn;
    CGNS_ENUMT(ElementType_t) elemtype, et;
    char name[33], outfile[65], buff[33];
    FILE *fp;

    if (cg_nsections (cgnsfn, cgnsbase, nz, &nsect))
        FATAL (NULL);
    if (nsect < 1) FATAL ("no sections defined");

    nodemap = (int *) malloc (nnodes * sizeof(int));
    if (nodemap == NULL)
        FATAL ("malloc failed for node mapping data");

    for (ns = 1; ns <= nsect; ns++) {
        if (cg_section_read (cgnsfn, cgnsbase, nz, ns,
                name, &elemtype, &is, &ie, &nn, &ip) ||
            cg_ElementDataSize (cgnsfn, cgnsbase, nz, ns, &size))
            FATAL (NULL);
        if (size > CG_MAX_INT32)
            FATAL("element data too large for 32-bit integer");
        nelems = ie - is + 1;
        conn = (cgsize_t *) malloc ((size_t)size * sizeof(cgsize_t));
        if (conn == NULL)
            FATAL ("malloc failed for element connectivity");
        if (cg_elements_read (cgnsfn, cgnsbase, nz, ns, conn, NULL))
            FATAL (NULL);

        for (n = 0; n < nnodes; n++)
            nodemap[n] = 0;
        et = elemtype;
        elemcnt = elemsize = 0;
        for (is = 0, ie = 0; ie < nelems; ie++) {
            if (elemtype == CGNS_ENUMV(MIXED)) et = (int)conn[is++];
            switch (et) {
                case CGNS_ENUMV(NODE):
                    nn = 1;
                    break;
                case CGNS_ENUMV(BAR_2):
                case CGNS_ENUMV(BAR_3):
                    nn = 2;
                    break;
                case CGNS_ENUMV(TRI_3):
                case CGNS_ENUMV(TRI_6):
                    nn = 3;
                    break;
                case CGNS_ENUMV(QUAD_4):
                case CGNS_ENUMV(QUAD_8):
                case CGNS_ENUMV(QUAD_9):
                case CGNS_ENUMV(TETRA_4):
                case CGNS_ENUMV(TETRA_10):
                    nn = 4;
                    break;
                case CGNS_ENUMV(PYRA_5):
                case CGNS_ENUMV(PYRA_14):
                    nn = 5;
                    break;
                case CGNS_ENUMV(PENTA_6):
                case CGNS_ENUMV(PENTA_15):
                case CGNS_ENUMV(PENTA_18):
                    nn = 6;
                    break;
                case CGNS_ENUMV(HEXA_8):
                case CGNS_ENUMV(HEXA_20):
                case CGNS_ENUMV(HEXA_27):
                    nn = 8;
                    break;
                default:
                    nn = 0;
                    break;
            }
            if (nn) {
                elemcnt++;
                elemsize += (nn + 1);
                for (i = 0; i < nn; i++) {
                    n = (int)conn[is+i] - 1;
                    (nodemap[n])++;
                }
            }
            if (cg_npe (et, &nn) || nn == 0)
                FATAL ("invalid element type");
            is += nn;
        }
        if (elemcnt == 0) {
            free (conn);
            continue;
        }

        create_filename (name, outfile);
        if (nzones > 1) {
            sprintf (buff, ".%d", nz);
            strcat (outfile, buff);
        }
        strcat (outfile, ".vtk");
        printf ("writing element set %d to \"%s\"\n", ns, outfile);
        fflush (stdout);
        if ((fp = fopen (outfile, "w+b")) == NULL) {
            fprintf (stderr, "couldn't open <%s> for output\n", outfile);
            exit (1);
        }

        fprintf (fp, "# vtk DataFile Version 2.0\n");
        fprintf (fp, "zone %d, elemset %d - %s\n", nz, ns, name);
        fprintf (fp, "%s\n", ascii ? "ASCII" : "BINARY");
        fprintf (fp, "DATASET UNSTRUCTURED_GRID\n");

        /* write the points */

        for (nn = 0, n = 0; n < nnodes; n++) {
            if (nodemap[n]) {
                nodemap[n] = ++nn;
            }
        }

        if (verbose) {
            printf ("  writing %d points\n", nn);
            fflush (stdout);
        }
        fprintf (fp, "POINTS %d float\n", nn);
        for (n = 0; n < nnodes; n++) {
            if (nodemap[n]) {
                write_floats (fp, 3, nodes[n]);
            }
        }

        /* write the cells */

        types = (int *) malloc (elemcnt * sizeof(int));
        if (types == NULL)
            FATAL ("malloc failed for element types");

        if (verbose) {
            printf ("  writing %d cells\n", elemcnt);
            fflush (stdout);
        }
        fprintf (fp, "CELLS %d %d\n", elemcnt, elemsize);

        et = elemtype;
        elemcnt = 0;
        for (is = 0, ie = 0; ie < nelems; ie++) {
            if (elemtype == CGNS_ENUMV(MIXED)) et = (int)conn[is++];
            switch (et) {
                case CGNS_ENUMV(NODE):
                    nn = 1;
                    types[elemcnt++] = 1;
                    break;
                case CGNS_ENUMV(BAR_2):
                case CGNS_ENUMV(BAR_3):
                    nn = 2;
                    types[elemcnt++] = 3;
                    break;
                case CGNS_ENUMV(TRI_3):
                case CGNS_ENUMV(TRI_6):
                    nn = 3;
                    types[elemcnt++] = 5;
                    break;
                case CGNS_ENUMV(QUAD_4):
                case CGNS_ENUMV(QUAD_8):
                case CGNS_ENUMV(QUAD_9):
                    nn = 4;
                    types[elemcnt++] = 9;
                    break;
                case CGNS_ENUMV(TETRA_4):
                case CGNS_ENUMV(TETRA_10):
                    nn = 4;
                    types[elemcnt++] = 10;
                    break;
                case CGNS_ENUMV(PYRA_5):
                case CGNS_ENUMV(PYRA_14):
                    nn = 5;
                    types[elemcnt++] = 14;
                    break;
                case CGNS_ENUMV(PENTA_6):
                case CGNS_ENUMV(PENTA_15):
                case CGNS_ENUMV(PENTA_18):
                    nn = 6;
                    types[elemcnt++] = 13;
                    break;
                case CGNS_ENUMV(HEXA_8):
                case CGNS_ENUMV(HEXA_20):
                case CGNS_ENUMV(HEXA_27):
                    nn = 8;
                    types[elemcnt++] = 12;
                    break;
                default:
                    nn = 0;
                    break;
            }
            if (nn) {
                cell[0] = nn;
                for (i = 0; i < nn; i++)
                    cell[i+1] = nodemap[(int)conn[is+i]-1] - 1;
                write_ints (fp, nn+1, cell);
            }
            cg_npe (et, &nn);
            is += nn;
        }

        free (conn);

        /* write the cell types */

        fprintf (fp, "CELL_TYPES %d\n", elemcnt);
        write_ints (fp, elemcnt, types);

        free (types);

        /* write solution if Vertex */

        if (nvars) {
            if (varloc == CGNS_ENUMV(Vertex) && ndata == nnodes)
                write_solution (fp, nz, nodemap);
            else if (verbose) {
                printf ("  skipping solution - not Vertex\n");
                fflush (stdout);
            }
        }

        fclose (fp);
    }

    free (nodemap);
}

/*========== main =====================================================*/

int main (int argc, char *argv[])
{
    int n, nz;
    char name[33], outfile[37];
    int elemsets = 0;
    cgsize_t sizes[9];
    CGNS_ENUMT(ZoneType_t) zonetype;
    struct stat st;
    FILE *fp;

    if (argc < 2)
        print_usage (usgmsg, NULL);

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'b':
                cgnsbase = atoi (argarg);
                break;
            case 'z':
                cgnszone = atoi (argarg);
                break;
            case 's':
                cgnssol = atoi (argarg);
                break;
            case 'e':
                elemsets = 1;
                break;
            case 'v':
                verbose = 1;
                break;
            case 'a':
                ascii = 1;
                break;
        }
    }

    if (argind >= argc)
        print_usage (usgmsg, "filename not specified");
    if (stat (argv[argind], &st)) {
        fprintf (stderr, "can't stat <%s>\n", argv[argind]);
        exit (1);
    }
    if (S_IFREG != (st.st_mode & S_IFMT)) {
        fprintf (stderr, "<%s> is not a regular file\n", argv[argind]);
        exit (1);
    }

    /* open CGNS file */

    printf ("reading CGNS file from \"%s\"\n", argv[argind]);
    fflush (stdout);
    if (cg_open (argv[argind], CG_MODE_READ, &cgnsfn))
        FATAL (NULL);
    if (cg_base_read (cgnsfn, cgnsbase, name, &CellDim, &PhyDim))
        FATAL (NULL);
    printf ("  using base %d - %s\n", cgnsbase, name);
    fflush (stdout);
    if (PhyDim != 3 /*|| (CellDim != 1 && CellDim != 3)*/)
        FATAL ("cell and/or physical dimension invalid");

    if (cg_nzones (cgnsfn, cgnsbase, &nzones))
        FATAL (NULL);
    if (nzones == 0)
        FATAL ("no zones in the CGNS file");
    if (cgnszone && cgnszone > nzones)
        FATAL ("zone number invalid");

    /* file output directory */

    if (++argind < argc) {
        if (stat (argv[argind], &st) &&
#ifdef _WIN32
            _mkdir (argv[argind])) {
#else
            mkdir (argv[argind], S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH)) {
#endif
            cg_close (cgnsfn);
            fprintf (stderr, "couldn't create the directory <%s>\n",
                argv[argind]);
            exit (1);
        }
        if (chdir (argv[argind])) {
            cg_close (cgnsfn);
            fprintf (stderr, "couldn't chdir to <%s>\n", argv[argind]);
            exit (1);
        }
        printf ("writing %s VTK files to directory \"%s\"\n",
            ascii ? "ASCII" : "binary", argv[argind]);
    }
    else
        printf ("writing %s VTK files to current directory\n",
            ascii ? "ASCII" : "binary");

    for (nz = 1; nz <= nzones; nz++) {
        if (cgnszone && nz != cgnszone) continue;
        if (cg_zone_type (cgnsfn, cgnsbase, nz, &zonetype) ||
            cg_zone_read (cgnsfn, cgnsbase, nz, name, sizes))
            FATAL (NULL);
        if (zonetype == CGNS_ENUMV(Structured)) {
            if (sizes[0]*sizes[1]*sizes[2] > CG_MAX_INT32)
                FATAL("too many coordinates for 32-bit integer");
        }
        else if (zonetype == CGNS_ENUMV(Unstructured)) {
            if (sizes[0] > CG_MAX_INT32)
                FATAL("too many coordinates for 32-bit integer");
        }
        else
            FATAL ("invalid zone type");
        

        create_filename (name, outfile);
        strcat (outfile, ".vtk");
        printf ("writing zone %d as %s to \"%s\"\n", nz,
#if CGNS_VERSION >= 2500
            cg_ZoneTypeName(zonetype), outfile);
#else
            ZoneTypeName[zonetype], outfile);
#endif
        fflush (stdout);
        if ((fp = fopen (outfile, "w+b")) == NULL) {
            fprintf (stderr, "couldn't open <%s> for output\n", outfile);
            exit (1);
        }

        fprintf (fp, "# vtk DataFile Version 2.0\n");
        fprintf (fp, "zone %d - %s\n", nz, name);
        fprintf (fp, "%s\n", ascii ? "ASCII" : "BINARY");
        if (zonetype == CGNS_ENUMV(Structured)) {
            fprintf (fp, "DATASET STRUCTURED_GRID\n");
            fprintf (fp, "DIMENSIONS %d %d %d\n",
                (int)sizes[0], (int)sizes[1], (int)sizes[2]);
        }
        else
            fprintf (fp, "DATASET UNSTRUCTURED_GRID\n");

        get_nodes (nz, zonetype, sizes);
        if (verbose) {
            printf ("  writing %d points\n", nnodes);
            fflush (stdout);
        }
        fprintf (fp, "POINTS %d float\n", nnodes);
        for (n = 0; n < nnodes; n++)
            write_floats (fp, 3, nodes[n]);

        if (zonetype == CGNS_ENUMV(Unstructured))
            write_volume_cells (fp, nz);
        if (get_variables (nz, zonetype, sizes))
            write_solution (fp, nz, NULL);

        fclose (fp);

        if (elemsets && zonetype == CGNS_ENUMV(Unstructured))
            write_element_sets (nz, sizes);

        free (nodes);
        if (nvars) free (vars);
    }

    cg_close (cgnsfn);
    return 0;
}
