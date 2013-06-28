#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#ifdef _WIN32
# define WIN32_LEAN_AND_MEAN
# include <windows.h>
#endif
#include <GL/gl.h>

#include "tk.h"
#include "cgnslib.h"
#include "hash.h"

/* define this to exclude structured mesh boundaries */
/* imin, imax, jmin, jmax, kmin, kmax */

/*#define NO_MESH_BOUNDARIES*/

/* define this to disable a cutting plane */
/* the cutting plane requires more memory since all
   nodes and elements are saved */

/*#define NO_CUTTING_PLANE*/

/* this is the cosine of the angle between faces to define an interior edge */
/* it's set to 60 degrees - comment to remove interior edges */

#define EDGE_ANGLE 0.5

typedef float Node[3];

typedef struct {
    cgsize_t id;
    cgsize_t nodes[2];
} Edge;

typedef struct {
    cgsize_t id;
    int flags;
    int nnodes;
    cgsize_t *nodes;
    float normal[3];
} Face;

typedef struct {
    Face *face;
    cgsize_t num;
    int flags;
} PolyFace;

#ifndef NO_CUTTING_PLANE

typedef struct {
    int id;
    cgsize_t nodes[3];
    float ratio;
} CutNode;

typedef struct {
    cgsize_t nelems;
    cgsize_t *elems;
    cgsize_t nedges;
    Edge *edges;
} CutData;

typedef struct {
    float plane[4];
    cgsize_t nelems;
    cgsize_t nedges;
    cgsize_t nnodes;
    Node *nodes;
} CutPlane;

static float cutcolor[4] = {(float)0.8, (float)0.4, (float)0.8, (float)0.5};
static int usecutclr = 0;
static int ignorevis = 0;
static CutPlane cutplane;
static int CutDL = 0;
static int PlaneDL = 0;

#endif

typedef struct {
    char name[33];
    int type;
    int dim;
    cgsize_t data[10];
    char d_name[33];
    cgsize_t nedges;
    Edge *edges;
    cgsize_t nfaces;
    Face **faces;
#ifndef NO_CUTTING_PLANE
    CGNS_ENUMT(ElementType_t) elemtype;
    cgsize_t nelems;
    cgsize_t *elems;
    int npoly;
    Face **poly;
    CutData cut;
#endif
    float bbox[3][2];
    int dlist;
    int mode;
    float color[4];
    char errmsg[81];
} Regn;

typedef struct {
    char name[33];
    cgsize_t nnodes;
    Node *nodes;
    int nregs;
    Regn *regs;
} Zone;

static int cgnsfn = 0;
static int cgnsbase = 0;
static int cgnszone = 0;
static int nbases = 0;
static int nzones = 0;
static Zone *zones;
static int AxisDL = 0;

static char BaseName[33];
static int CellDim, PhyDim;

static Tcl_Interp *global_interp;

enum {
    REG_MESH,
    REG_ELEM,
    REG_1TO1,
    REG_CONN,
    REG_HOLE,
    REG_BOCO,
    REG_BNDS
};

/* mapping from elements to faces */

static int facenodes[22][5] = {
    /* tri */
    {3, 0, 1, 2, 0},
    /* quad */
    {4, 0, 1, 2, 3},
    /* tet */
    {3, 0, 2, 1, 0},
    {3, 0, 1, 3, 0},
    {3, 1, 2, 3, 0},
    {3, 2, 0, 3, 0},
    /* pyramid */
    {4, 0, 3, 2, 1},
    {3, 0, 1, 4, 0},
    {3, 1, 2, 4, 0},
    {3, 2, 3, 4, 0},
    {3, 3, 0, 4, 0},
    /* wedge */
    {4, 0, 1, 4, 3},
    {4, 1, 2, 5, 4},
    {4, 2, 0, 3, 5},
    {3, 0, 2, 1, 0},
    {3, 3, 4, 5, 0},
    /* hex */
    {4, 0, 3, 2, 1},
    {4, 0, 1, 5, 4},
    {4, 1, 2, 6, 5},
    {4, 2, 3, 7, 6},
    {4, 0, 4, 7, 3},
    {4, 4, 5, 6, 7}
};

/*===================================================================
 *           local routines
 *===================================================================*/

/*-------------------------------------------------------------------*/

static void FATAL (char *errmsg)
{
    char cmd[129];

    sprintf (cmd, "error_exit {%s}", errmsg);
    Tcl_Eval (global_interp, cmd);
    exit (1);
}

/*-------------------------------------------------------------------*/

static void *MALLOC (char *funcname, size_t bytes)
{
    void *data = calloc (bytes, 1);
    if (NULL == data) {
        char msg[128];
        if (funcname != NULL)
            sprintf (msg, "%s:malloc failed for %lu bytes", funcname,
                (unsigned long)bytes);
        else
            sprintf (msg, "malloc failed for %lu bytes",
                (unsigned long)bytes);
        FATAL (msg);
    }
    return data;
}

/*-------------------------------------------------------------------*/

static void *REALLOC (char *funcname, size_t bytes, void *old_data)
{
    void *new_data;

    if (NULL == old_data)
        return MALLOC (funcname, bytes);
    new_data = realloc (old_data, bytes);
    if (NULL == new_data) {
        char msg[128];
        if (funcname != NULL)
            sprintf (msg, "%s:realloc failed for %lu bytes", funcname,
                (unsigned long)bytes);
        else
            sprintf (msg, "realloc failed for %lu bytes",
                (unsigned long)bytes);
        FATAL (msg);
    }
    return new_data;
}

/*-------------------------------------------------------------------*/

static Face *new_face(char *funcname, int nnodes)
{
    Face *f = (Face *) calloc (nnodes * sizeof(cgsize_t) + sizeof(Face), 1);
    if (f == NULL) {
        char msg[128];
        if (funcname != NULL)
            sprintf (msg, "%s:malloc failed for face with %d nodes",
                funcname, nnodes);
        else
            sprintf (msg, "malloc failed for face with %d nodes",
                nnodes);
        FATAL (msg);
    }
    f->nnodes = nnodes;
    f->nodes  = (cgsize_t *)(f + 1);
    return f;
}

/*-------------------------------------------------------------------*/

static Face *copy_face(char *funcname, Face *face)
{
    int n;
    Face *f = new_face(funcname, face->nnodes);

    f->id = face->id;
    f->flags = face->flags;
    for (n = 0; n < face->nnodes; n++)
        f->nodes[n] = face->nodes[n];
    for (n = 0; n < 3; n++)
        f->normal[n] = face->normal[n];
    return f;
}

/*-------------------------------------------------------------------*/

static void zone_message (char *msg, char *name)
{
    char cmd[129];

    if (name != NULL)
        sprintf (cmd, "display_message {Zone %d : %s %s...}",
            cgnszone, msg, name);
    else
        sprintf (cmd, "display_message {Zone %d : %s...}", cgnszone, msg);
    Tcl_Eval (global_interp, cmd);
}

/*-------------------------------------------------------------------*/

static void free_all (void)
{
    int nz, nr, nf;

    if (!nzones) return;
    for (nz = 0; nz < nzones; nz++) {
        for (nr = 0; nr < zones[nz].nregs; nr++) {
            if (zones[nz].regs[nr].dlist)
                glDeleteLists (zones[nz].regs[nr].dlist, 1);
            if (zones[nz].regs[nr].nedges)
                free (zones[nz].regs[nr].edges);
            if (zones[nz].regs[nr].nfaces) {
                for (nf = 0; nf < zones[nz].regs[nr].nfaces; nf++) {
                    if (zones[nz].regs[nr].faces[nf])
                        free (zones[nz].regs[nr].faces[nf]);
                }
                free (zones[nz].regs[nr].faces);
            }
#ifndef NO_CUTTING_PLANE
            if (zones[nz].regs[nr].nelems)
                free (zones[nz].regs[nr].elems);
            if (zones[nz].regs[nr].npoly)
                free (zones[nz].regs[nr].poly);
            if (zones[nz].regs[nr].cut.nelems)
                free (zones[nz].regs[nr].cut.elems);
            if (zones[nz].regs[nr].cut.nedges)
                free (zones[nz].regs[nr].cut.edges);
#endif
        }
        if (zones[nz].nregs) free (zones[nz].regs);
        if (zones[nz].nnodes) free(zones[nz].nodes);
    }
    free (zones);
    nzones = 0;
#ifndef NO_CUTTING_PLANE
    if (cutplane.nnodes)
        free (cutplane.nodes);
    cutplane.nelems = 0;
    cutplane.nedges = 0;
    cutplane.nnodes = 0;
#endif
}

/*-------------------------------------------------------------------*/

static int compare_ints (const void *v1, const void *v2)
{
    return (int)(*((cgsize_t *)v1) - *((cgsize_t *)v2));
}

/*-------------------------------------------------------------------*/

static int find_int (cgsize_t value, cgsize_t nlist, cgsize_t *list)
{
    cgsize_t lo = 0, hi = nlist - 1, mid;

    if (!nlist || value < list[0]) return 0;
    if (value == list[0]) return 1;
    if (!hi || value > list[hi]) return 0;
    if (value == list[hi]) return 1;

    while (lo <= hi) {
        mid = (lo + hi) >> 1;
        if (value == list[mid]) return 1;
        if (value < list[mid])
            hi = mid - 1;
        else
            lo = mid + 1;
    }
    return 0;
}

/*====================================================================
 * structured grid regions
 *====================================================================*/

#define NODE_INDEX(I,J,K) ((I)+dim[0]*((J)+dim[1]*(K)))

/*-------------------------------------------------------------------*/

static int structured_range (Regn *reg, cgsize_t *dim, cgsize_t *ptrng,
                             CGNS_ENUMT(GridLocation_t) location)
{
    int n, i, j, nf;
    cgsize_t ii, jj, kk, nfaces, rng[3][2];
    Face *f;
    static char *funcname = "structured_range";

    if (location != CGNS_ENUMV(Vertex) &&
        location != CGNS_ENUMV(IFaceCenter) &&
        location != CGNS_ENUMV(JFaceCenter) &&
        location != CGNS_ENUMV(KFaceCenter)) {
        i = j = 0;
        for (n = 0; n < CellDim; n++) {
            if (ptrng[n] == ptrng[n+CellDim] &&
               (ptrng[n] == 1 || ptrng[n] == dim[n])) {
                if (ptrng[n] == 1)
                    i++;
                else if (j) {
                    j = 4;
                    break;
                }
                else
                    j = n + 1;
            }
        }
        if (!j && i == 1) {
            for (n = 0; n < CellDim; n++) {
                if (ptrng[n] == ptrng[n+CellDim] && ptrng[n] == 1) {
                    j = n + 1;
                    break;
                }
            }
        }
        if (j == 1)
            location = CGNS_ENUMV(IFaceCenter);
        else if (j == 2)
            location = CGNS_ENUMV(JFaceCenter);
        else if (j == 3)
            location = CGNS_ENUMV(KFaceCenter);
        else {
            strcpy (reg->errmsg,
                "unable to determine boundary - use [IJK]FaceCenter");
            return 0;
        }
    }

    nfaces = 1;
    if (location == CGNS_ENUMV(Vertex)) {
        for (n = 0, i = 0; i < CellDim; i++) {
            if (ptrng[i] < 1 || ptrng[i] > dim[i]) return 0;
            if (ptrng[i] == ptrng[i+CellDim]) {
                if (n || (ptrng[i] != 1 && ptrng[i] != dim[i]))
                    return 0;
                n = i + 1;
                rng[i][0] = rng[i][1] = ptrng[i] - 1;
            } else {
                if (ptrng[i] < ptrng[i+CellDim]) {
                    rng[i][0] = ptrng[i] - 1;
                    rng[i][1] = ptrng[i+CellDim] - 1;
                }
                else {
                    rng[i][0] = ptrng[i+CellDim] - 1;
                    rng[i][1] = ptrng[i] - 1;
                }
                nfaces *= (rng[i][1] - rng[i][0]);
            }
        }
    }
    else {
        if (location == CGNS_ENUMV(IFaceCenter))
            n = 0;
        else if (location == CGNS_ENUMV(JFaceCenter))
            n = 1;
        else
            n = 2;
        for (i = 0; i < CellDim; i++) {
            if (i == n) {
                if (ptrng[i] != ptrng[i+CellDim] ||
                   (ptrng[i] != 1 && ptrng[i] != dim[i])) return 0;
                rng[i][0] = rng[i][1] = ptrng[i] - 1;
            }
            else {
                if (ptrng[i] < ptrng[i+CellDim]) {
                    rng[i][0] = ptrng[i] - 1;
                    rng[i][1] = ptrng[i+CellDim];
                }
                else {
                    rng[i][0] = ptrng[i+CellDim] - 1;
                    rng[i][1] = ptrng[i];
                }
                if (rng[i][0] < 0 || rng[i][1] >= dim[i]) return 0;
                nfaces *= (rng[i][1] - rng[i][0]);
            }
        }
        n++;
    }
    if (!nfaces || n < 1 || n > CellDim) {
        strcpy (reg->errmsg, "couldn't find any exterior faces");
        return 0;
    }

    if (CellDim == 2) {
        reg->nedges = nfaces;
        reg->edges  = (Edge *) MALLOC (funcname,  (size_t)nfaces * sizeof(Edge));
        nf = 0;
        kk = 0;        

        if (n == 1) {
            ii = rng[0][0];
            for (jj = rng[1][0]; jj < rng[1][1]; jj++) {
                reg->edges[nf].nodes[0] = NODE_INDEX(ii, jj,   kk);
                reg->edges[nf].nodes[1] = NODE_INDEX(ii, jj+1, kk);
                nf++;
            }
        }
        else {
            jj = rng[1][0];
            for (ii = rng[0][0]; ii < rng[0][1]; ii++) {
                reg->edges[nf].nodes[0] = NODE_INDEX(ii,   jj, kk);
                reg->edges[nf].nodes[1] = NODE_INDEX(ii+1, jj, kk);
                nf++;
            }
        }
        return 1;
    }

    reg->nfaces = nfaces;
    reg->faces  = (Face **) MALLOC (funcname,  (size_t)nfaces * sizeof(Face *));
    for (nf = 0; nf < nfaces; nf++)
        reg->faces[nf] = new_face (funcname, 4);
    nf = 0;

    if (n == 1) {
        if ((ii = rng[0][0]) == 0) {
            for (kk = rng[2][0]; kk < rng[2][1]; kk++) {
                for (jj = rng[1][0]; jj < rng[1][1]; jj++) {
                    f = reg->faces[nf++];
                    f->nodes[0] = NODE_INDEX (ii, jj,   kk);
                    f->nodes[1] = NODE_INDEX (ii, jj,   kk+1);
                    f->nodes[2] = NODE_INDEX (ii, jj+1, kk+1);
                    f->nodes[3] = NODE_INDEX (ii, jj+1, kk);
                }
            }
        }
        else {
            for (kk = rng[2][0]; kk < rng[2][1]; kk++) {
                for (jj = rng[1][0]; jj < rng[1][1]; jj++) {
                    f = reg->faces[nf++];
                    f->nodes[0] = NODE_INDEX (ii, jj,   kk);
                    f->nodes[1] = NODE_INDEX (ii, jj+1, kk);
                    f->nodes[2] = NODE_INDEX (ii, jj+1, kk+1);
                    f->nodes[3] = NODE_INDEX (ii, jj,   kk+1);
                }
            }
        }
    }
    else if (n == 2) {
        if ((jj = rng[1][0]) == 0) {
            for (kk = rng[2][0]; kk < rng[2][1]; kk++) {
                for (ii = rng[0][0]; ii < rng[0][1]; ii++) {
                    f = reg->faces[nf++];
                    f->nodes[0] = NODE_INDEX (ii,   jj, kk);
                    f->nodes[1] = NODE_INDEX (ii+1, jj, kk);
                    f->nodes[2] = NODE_INDEX (ii+1, jj, kk+1);
                    f->nodes[3] = NODE_INDEX (ii,   jj, kk+1);
                }
            }
        }
        else {
            for (kk = rng[2][0]; kk < rng[2][1]; kk++) {
                for (ii = rng[0][0]; ii < rng[0][1]; ii++) {
                    f = reg->faces[nf++];
                    f->nodes[0] = NODE_INDEX (ii,   jj, kk);
                    f->nodes[1] = NODE_INDEX (ii,   jj, kk+1);
                    f->nodes[2] = NODE_INDEX (ii+1, jj, kk+1);
                    f->nodes[3] = NODE_INDEX (ii+1, jj, kk);
                }
            }
        }
    }
    else {
        if ((kk = rng[2][0]) == 0) {
            for (jj = rng[1][0]; jj < rng[1][1]; jj++) {
                for (ii = rng[0][0]; ii < rng[0][1]; ii++) {
                    f = reg->faces[nf++];
                    f->nodes[0] = NODE_INDEX (ii,   jj,   kk);
                    f->nodes[1] = NODE_INDEX (ii,   jj+1, kk);
                    f->nodes[2] = NODE_INDEX (ii+1, jj+1, kk);
                    f->nodes[3] = NODE_INDEX (ii+1, jj,   kk);
                }
            }
        }
        else {
            for (jj = rng[1][0]; jj < rng[1][1]; jj++) {
                for (ii = rng[0][0]; ii < rng[0][1]; ii++) {
                    f = reg->faces[nf++];
                    f->nodes[0] = NODE_INDEX (ii,   jj,   kk);
                    f->nodes[1] = NODE_INDEX (ii+1, jj,   kk);
                    f->nodes[2] = NODE_INDEX (ii+1, jj+1, kk);
                    f->nodes[3] = NODE_INDEX (ii,   jj+1, kk);
                }
            }
        }
    }
    return 2;
}

/*-------------------------------------------------------------------*/

static int structured_list (Regn *mesh, Regn *reg, cgsize_t *dim, cgsize_t npts,
                            cgsize_t *pts, CGNS_ENUMT(GridLocation_t) location)
{
    cgsize_t n, nn, nf, nfaces = 0, noff, nmax;
    Face **faces, *f;
    static char *funcname = "structured_list";

    if (location != CGNS_ENUMV(Vertex) &&
        location != CGNS_ENUMV(IFaceCenter) &&
        location != CGNS_ENUMV(JFaceCenter) &&
        location != CGNS_ENUMV(KFaceCenter)) {
        int i, j, k;
        cgsize_t rng[3][2];
        for (i = 0; i < CellDim; i++)
            rng[i][0] = rng[i][1] = pts[i];
        for (nf = 1; nf < npts; nf++) {
            nn = nf * CellDim;
            for (i = 0; i < CellDim; i++) {
                if (rng[i][0] > pts[nn+i]) rng[i][0] = pts[nn+i];
                if (rng[i][1] < pts[nn+i]) rng[i][1] = pts[nn+i];
            }
        }
        j = k = 0;
        for (i = 0; i < CellDim; i++) {
            if (rng[i][0] == rng[i][1] &&
                (rng[i][0] == 1 || rng[i][0] == dim[i])) {
                if (rng[i][0] == 1)
                    j++;
                else if (k) {
                    k = 4;
                    break;
                }
                else
                    k = i + 1;
            }
        }
        if (!k && j == 1) {
            for (i = 0; i < CellDim; i++) {
                if (rng[i][0] == rng[i][1] && rng[i][0] == 1) {
                    k = i + 1;
                    break;
                }
            }
        }
        if (k == 1)
            location = CGNS_ENUMV(IFaceCenter);
        else if (k == 2)
            location = CGNS_ENUMV(JFaceCenter);
        else if (k == 3)
            location = CGNS_ENUMV(KFaceCenter);
        else {
            strcpy (reg->errmsg,
                "unable to determine boundary - use [IJK]FaceCenter");
            return 0;
        }
    }
    nmax  = npts;

    if (CellDim == 2) {
        cgsize_t ii, jj, n0, n1, ne;
        Edge *edges = (Edge *) MALLOC (funcname, (size_t)nmax * sizeof(Edge));

        ne = 0;
        if (location == CGNS_ENUMV(Vertex)) {
            for (nn = 0, n = 0; n < npts; n++) {
                pts[n] = NODE_INDEX (pts[nn]-1, pts[nn+1]-1, 0);
                nn += 2;
            }
            for (n = 1; n < npts; n++) {
                if (pts[n] < pts[n-1]) {
                    qsort (pts, (size_t)npts, sizeof(cgsize_t), compare_ints);
                    break;
                }
            }

            ne = 0;
            jj = 0;
            for (ii = 1; ii < dim[0]; ii++) {
                n0 = NODE_INDEX(ii-1, jj, 0);
                n1 = NODE_INDEX(ii,   jj, 0);
                if (find_int(n0, npts, pts) &&
                    find_int(n1, npts, pts)) {
                    if (ne == nmax) {
                        nmax += 100;
                        edges = (Edge *) REALLOC (funcname,
                            (size_t)nmax * sizeof(Edge), edges);
                    }
                    edges[ne].nodes[0] = n0;
                    edges[ne].nodes[1] = n1;
                    ne++;
                }
            }

            jj = dim[1] - 1;
            for (ii = 1; ii < dim[0]; ii++) {
                n0 = NODE_INDEX(ii-1, jj, 0);
                n1 = NODE_INDEX(ii,   jj, 0);
                if (find_int(n0, npts, pts) &&
                    find_int(n1, npts, pts)) {
                    if (ne == nmax) {
                        nmax += 100;
                        edges = (Edge *) REALLOC (funcname,
                            (size_t)nmax * sizeof(Edge), edges);
                    }
                    edges[ne].nodes[0] = n0;
                    edges[ne].nodes[1] = n1;
                    ne++;
                }
            }

            ii = 0;
            for (jj = 1; jj < dim[1]; jj++) {
                n0 = NODE_INDEX(ii, jj-1, 0);
                n1 = NODE_INDEX(ii, jj,   0);
                if (find_int(n0, npts, pts) &&
                    find_int(n1, npts, pts)) {
                    if (ne == nmax) {
                        nmax += 100;
                        edges = (Edge *) REALLOC (funcname,
                            (size_t)nmax * sizeof(Edge), edges);
                    }
                    edges[ne].nodes[0] = n0;
                    edges[ne].nodes[1] = n1;
                    ne++;
                }
            }

            ii = dim[0] - 1;
            for (jj = 1; jj < dim[1]; jj++) {
                n0 = NODE_INDEX(ii, jj-1, 0);
                n1 = NODE_INDEX(ii, jj,   0);
                if (find_int(n0, npts, pts) &&
                    find_int(n1, npts, pts)) {
                    if (ne == nmax) {
                        nmax += 100;
                        edges = (Edge *) REALLOC (funcname,
                            (size_t)nmax * sizeof(Edge), edges);
                    }
                    edges[ne].nodes[0] = n0;
                    edges[ne].nodes[1] = n1;
                    ne++;
                }
            }
        }

        else if (location == CGNS_ENUMV(IFaceCenter)) {
            for (nn = 0, n = 0; n < npts; n++) {
                if ((pts[nn] == 1 || pts[nn] == dim[0]) &&
                    pts[nn+1] > 0 && pts[nn+1] < dim[1]) {
                    ii = pts[nn] - 1;
                    jj = pts[nn+1] - 1;
                    n0 = NODE_INDEX(ii, jj,   0);
                    n1 = NODE_INDEX(ii, jj+1, 0);
                    edges[ne].nodes[0] = n0;
                    edges[ne].nodes[1] = n1;
                    ne++;
                }
                nn += 2;
            }
        }

        else if (location == CGNS_ENUMV(JFaceCenter)) {
            for (nn = 0, n = 0; n < npts; n++) {
                if ((pts[nn+1] == 1 || pts[nn+1] == dim[1]) &&
                    pts[nn] > 0 && pts[nn] < dim[0]) {
                    ii = pts[nn] - 1;
                    jj = pts[nn+1] - 1;
                    n0 = NODE_INDEX(ii,   jj, 0);
                    n1 = NODE_INDEX(ii+1, jj, 0);
                    edges[ne].nodes[0] = n0;
                    edges[ne].nodes[1] = n1;
                    ne++;
                }
                nn += 2;
            }
        }
        
        if (ne == 0) {
            free(edges);
            strcpy (reg->errmsg, "couldn't find any exterior edges");
            return 0;
        }

        reg->nedges = ne;
        reg->edges  = edges;
        return 1;
    }

    faces = (Face **) MALLOC (funcname, (size_t)nmax * sizeof(Face *));

    if (location == CGNS_ENUMV(Vertex)) {
        for (nn = 0, n = 0; n < npts; n++) {
            pts[n] = NODE_INDEX (pts[nn]-1, pts[nn+1]-1, pts[nn+2]-1);
            nn += 3;
        }
        for (n = 1; n < npts; n++) {
            if (pts[n] < pts[n-1]) {
                qsort (pts, (size_t)npts, sizeof(cgsize_t), compare_ints);
                break;
            }
        }

        for (nf = 0; nf < mesh->nfaces; nf++) {
            f = mesh->faces[nf];
            for (nn = 0; nn < f->nnodes; nn++) {
                if (!find_int (f->nodes[nn], npts, pts))
                    break;
            }
            if (nn == f->nnodes) {
                if (nfaces == nmax) {
                    nmax += 100;
                    faces = (Face **) REALLOC (funcname,
                        (size_t)nmax * sizeof(Face *), faces);
                }
                faces[nfaces++] = copy_face (funcname, f);
            }
        }
    }

    else if (location == CGNS_ENUMV(IFaceCenter)) {
        for (n = 0; n < npts; n++) {
            nn = 3 * n;
            if ((pts[nn] == 1 || pts[nn] == dim[0]) &&
                pts[nn+1] > 0 && pts[nn+1] < dim[1] &&
                pts[nn+2] > 0 && pts[nn+2] < dim[2]) {
                nf = pts[nn+1]-1 + (pts[nn+2]-1) * (dim[1]-1);
                if (pts[nn] == dim[0])
                    nf += (dim[1]-1) * (dim[2]-1);
                faces[nfaces++] = copy_face (funcname, mesh->faces[nf]);
            }
        }
    }

    else if (location == CGNS_ENUMV(JFaceCenter)) {
        noff = 2 * (dim[1]-1) * (dim[2]-1);
        for (n = 0; n < npts; n++) {
            nn = 3 * n;
            if ((pts[nn+1] == 1 || pts[nn+1] == dim[1]) &&
                pts[nn] > 0 && pts[nn] < dim[0] &&
                pts[nn+2] > 0 && pts[nn+2] < dim[2]) {
                nf = noff + pts[nn]-1 + (pts[nn+2]-1) * (dim[0]-1);
                if (pts[nn+1] == dim[1])
                    nf += (dim[0]-1) * (dim[2]-1);
                faces[nfaces++] = copy_face (funcname, mesh->faces[nf]);
            }
        }
    }

    else  {
        noff = 2 * ((dim[1]-1) * (dim[2]-1) + (dim[0]-1) * (dim[2]-1));
        for (n = 0; n < npts; n++) {
            nn = 3 * n;
            if ((pts[nn+2] == 1 || pts[nn+2] == dim[2]) &&
                pts[nn] > 0 && pts[nn] < dim[0] &&
                pts[nn+1] > 0 && pts[nn+1] < dim[1]) {
                nf = noff + pts[nn]-1 + (pts[nn+1]-1) * (dim[0]-1);
                if (pts[nn+2] == dim[2])
                    nf += (dim[0]-1) * (dim[1]-1);
                faces[nfaces++] = copy_face (funcname, mesh->faces[nf]);
            }
        }
    }

    if (nfaces == 0) {
        free (faces);
        strcpy (reg->errmsg, "couldn't find any exterior faces");
        return 0;
    }

    reg->nfaces = nfaces;
    reg->faces = faces;
    return 2;
}

/*-------------------------------------------------------------------*/

static int structured_zone (Tcl_Interp *interp, cgsize_t *dim)
{
    char name[33], d_name[33];
    int nints, nconns, nbocos, nholes, ii, nn, nr, nsets;
    cgsize_t i, j, k, n, ni, nj, nk, nf, ne, fn;
    cgsize_t npts, *pts, ndpts;
    cgsize_t range[6], d_range[6];
    int transform[3];
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(GridConnectivityType_t) type;
    CGNS_ENUMT(PointSetType_t) ptype, d_ptype;
    CGNS_ENUMT(ZoneType_t) d_ztype;
    CGNS_ENUMT(DataType_t) datatype;
    CGNS_ENUMT(BCType_t) bctype;
    Face *f;
    Zone *z = &zones[cgnszone-1];
    static char *funcname = "structured_zone";

    zone_message ("finding exterior faces", NULL);
    if (cg_n1to1 (cgnsfn, cgnsbase, cgnszone, &nints) ||
        cg_nconns (cgnsfn, cgnsbase, cgnszone, &nconns) ||
        cg_nholes (cgnsfn, cgnsbase, cgnszone, &nholes) ||
        cg_nbocos (cgnsfn, cgnsbase, cgnszone, &nbocos)) {
        Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
        return 1;
    }
    z->nregs = nints + nconns + nholes + nbocos + 1;
#ifndef NO_MESH_BOUNDARIES
    z->nregs += (2 * CellDim);
#endif
    z->regs = (Regn *) MALLOC (funcname, z->nregs * sizeof(Regn));
    ni = dim[0] - 1;
    nj = dim[1] - 1;
    nk = dim[2] - 1;
    nr = 1;

    /* mesh boundaries */

    strcpy (z->regs[0].name, "<mesh>");
    z->regs[0].type = REG_MESH;
    for (n = 0; n < CellDim; n++) {
        z->regs[0].data[n] = dim[n];
        range[n] = 1;
        range[n+CellDim] = dim[n];
    }

    if (CellDim == 2) {
        z->regs[0].dim = 2;
        z->regs[0].nfaces = ni * nj;
        z->regs[0].faces = (Face **) MALLOC (funcname,
                           (size_t)z->regs[0].nfaces * sizeof(Face *));
        fn = 0;
        k  = 0;
        for (j = 0; j < nj; j++) {
            for (i = 0; i < ni; i++) {
                f = z->regs[0].faces[fn++] = new_face (funcname, 4);
                f->nodes[0] = NODE_INDEX (i, j, k);
                f->nodes[1] = NODE_INDEX (i, j+1, k);
                f->nodes[2] = NODE_INDEX (i+1, j+1, k);
                f->nodes[3] = NODE_INDEX (i+1, j, k);
            }
        }

#ifndef NO_CUTTING_PLANE
        nf = z->regs[0].nfaces;
        z->regs[0].elemtype = CGNS_ENUMV(QUAD_4);
        z->regs[0].nelems = nf;
        z->regs[0].elems = (cgsize_t *) MALLOC (funcname,
	                       (size_t)(4 * nf) * sizeof(cgsize_t));
        for (n = 0, j = 0; j < nf; j++) {
            f = z->regs[0].faces[j];
            for (i = 0; i < 4; i++)
                z->regs[0].elems[n++] = f->nodes[i];
        }
#endif

#ifndef NO_MESH_BOUNDARIES
        strcpy (z->regs[nr].name, "<imin>");
        z->regs[nr].dim = 1;
        z->regs[nr].type = REG_BNDS;
        for (nn = 0; nn < 4; nn++)
            z->regs[nr].data[nn] = range[nn];
        z->regs[nr].data[2] = 1;
        z->regs[nr].nedges = nj;
        z->regs[nr].edges = (Edge *) MALLOC (funcname,
                            (size_t)nj * sizeof(Edge));
        for (i = 0, j = 0; j < nj; j++) {
            z->regs[nr].edges[j].nodes[0] = NODE_INDEX(i, j,   k);
            z->regs[nr].edges[j].nodes[1] = NODE_INDEX(i, j+1, k);
        }
        nr++;

        strcpy (z->regs[nr].name, "<imax>");
        z->regs[nr].dim = 1;
        z->regs[nr].type = REG_BNDS;
        for (nn = 0; nn < 4; nn++)
            z->regs[nr].data[nn] = range[nn];
        z->regs[nr].data[0] = dim[0];
        z->regs[nr].nedges = nj;
        z->regs[nr].edges = (Edge *) MALLOC (funcname,
                            (size_t)nj * sizeof(Edge));
        for (i = ni, j = 0; j < nj; j++) {
            z->regs[nr].edges[j].nodes[0] = NODE_INDEX(i, j,   k);
            z->regs[nr].edges[j].nodes[1] = NODE_INDEX(i, j+1, k);
        }
        nr++;

        strcpy (z->regs[nr].name, "<jmin>");
        z->regs[nr].dim = 1;
        z->regs[nr].type = REG_BNDS;
        for (nn = 0; nn < 4; nn++)
            z->regs[nr].data[nn] = range[nn];
        z->regs[nr].data[3] = 1;
        z->regs[nr].nedges = ni;
        z->regs[nr].edges = (Edge *) MALLOC (funcname,
                            (size_t)ni * sizeof(Edge));
        for (j = 0, i = 0; i < ni; i++) {
            z->regs[nr].edges[i].nodes[0] = NODE_INDEX(i,   j, k);
            z->regs[nr].edges[i].nodes[1] = NODE_INDEX(i+1, j, k);
        }
        nr++;

        strcpy (z->regs[nr].name, "<jmax>");
        z->regs[nr].dim = 1;
        z->regs[nr].type = REG_BNDS;
        for (nn = 0; nn < 4; nn++)
            z->regs[nr].data[nn] = range[nn];
        z->regs[nr].data[1] = dim[1];
        z->regs[nr].nedges = ni;
        z->regs[nr].edges = (Edge *) MALLOC (funcname,
                            (size_t)ni * sizeof(Edge));
        for (j = nj, i = 0; i < ni; i++) {
            z->regs[nr].edges[i].nodes[0] = NODE_INDEX(i,   j, k);
            z->regs[nr].edges[i].nodes[1] = NODE_INDEX(i+1, j, k);
        }
        nr++;
#endif
    }

    else {
        z->regs[0].dim = 3;

#ifndef NO_CUTTING_PLANE
        z->regs[0].elemtype = CGNS_ENUMV(HEXA_8);
        z->regs[0].nelems = ni * nj * nk;
        z->regs[0].elems = pts = (cgsize_t *) MALLOC (funcname,
	        (size_t)(8 * z->regs[0].nelems) * sizeof(cgsize_t));

        for (n = 0, k = 0; k < nk; k++) {
            for (j = 0; j < nj; j++) {
                for (i = 0; i < ni; i++) {
                    pts[n++] = NODE_INDEX (i,   j,   k);
                    pts[n++] = NODE_INDEX (i+1, j,   k);
                    pts[n++] = NODE_INDEX (i+1, j+1, k);
                    pts[n++] = NODE_INDEX (i,   j+1, k);
                    pts[n++] = NODE_INDEX (i,   j,   k+1);
                    pts[n++] = NODE_INDEX (i+1, j,   k+1);
                    pts[n++] = NODE_INDEX (i+1, j+1, k+1);
                    pts[n++] = NODE_INDEX (i,   j+1, k+1);
                }
            }
        }
#endif

        z->regs[0].nfaces = 2 * (nj * nk + ni * nk + ni * nj);
        z->regs[0].faces = (Face **) MALLOC (funcname,
                           (size_t)z->regs[0].nfaces * sizeof(Face *));
        fn = 0;

        for (i = 0, k = 0; k < nk; k++) {
            for (j = 0; j < nj; j++) {
                f = z->regs[0].faces[fn++] = new_face (funcname, 4);
                f->nodes[0] = NODE_INDEX (i, j, k);
                f->nodes[1] = NODE_INDEX (i, j, k+1);
                f->nodes[2] = NODE_INDEX (i, j+1, k+1);
                f->nodes[3] = NODE_INDEX (i, j+1, k);
            }
        }
        for (i = ni, k = 0; k < nk; k++) {
            for (j = 0; j < nj; j++) {
                f = z->regs[0].faces[fn++] = new_face (funcname, 4);
                f->nodes[0] = NODE_INDEX (i, j, k);
                f->nodes[1] = NODE_INDEX (i, j+1, k);
                f->nodes[2] = NODE_INDEX (i, j+1, k+1);
                f->nodes[3] = NODE_INDEX (i, j, k+1);
            }
        }
        for (j = 0, k = 0; k < nk; k++) {
            for (i = 0; i < ni; i++) {
                f = z->regs[0].faces[fn++] = new_face (funcname, 4);
                f->nodes[0] = NODE_INDEX (i, j, k);
                f->nodes[1] = NODE_INDEX (i+1, j, k);
                f->nodes[2] = NODE_INDEX (i+1, j, k+1);
                f->nodes[3] = NODE_INDEX (i, j, k+1);
            }
        }
        for (j = nj, k = 0; k < nk; k++) {
            for (i = 0; i < ni; i++) {
                f = z->regs[0].faces[fn++] = new_face (funcname, 4);
                f->nodes[0] = NODE_INDEX (i, j, k);
                f->nodes[1] = NODE_INDEX (i, j, k+1);
                f->nodes[2] = NODE_INDEX (i+1, j, k+1);
                f->nodes[3] = NODE_INDEX (i+1, j, k);
            }
        }
        for (k = 0, j = 0; j < nj; j++) {
            for (i = 0; i < ni; i++) {
                f = z->regs[0].faces[fn++] = new_face (funcname, 4);
                f->nodes[0] = NODE_INDEX (i, j, k);
                f->nodes[1] = NODE_INDEX (i, j+1, k);
                f->nodes[2] = NODE_INDEX (i+1, j+1, k);
                f->nodes[3] = NODE_INDEX (i+1, j, k);
            }
        }
        for (k = nk, j = 0; j < nj; j++) {
            for (i = 0; i < ni; i++) {
                f = z->regs[0].faces[fn++] = new_face (funcname, 4);
                f->nodes[0] = NODE_INDEX (i, j, k);
                f->nodes[1] = NODE_INDEX (i+1, j, k);
                f->nodes[2] = NODE_INDEX (i+1, j+1, k);
                f->nodes[3] = NODE_INDEX (i, j+1, k);
            }
        }

#ifndef NO_MESH_BOUNDARIES
        fn = 0;
        nf = nj * nk;
        strcpy (z->regs[nr].name, "<imin>");
        z->regs[nr].dim = 2;
        z->regs[nr].type = REG_BNDS;
        for (nn = 0; nn < 6; nn++)
            z->regs[nr].data[nn] = range[nn];
        z->regs[nr].data[3] = 1;
        z->regs[nr].nfaces = nf;
        z->regs[nr].faces = (Face **) MALLOC (funcname,
                            (size_t)nf * sizeof(Face *));
        for (n = 0; n < nf; n++)
            z->regs[nr].faces[n] = copy_face (funcname, z->regs->faces[fn++]);
        nr++;

        strcpy (z->regs[nr].name, "<imax>");
        z->regs[nr].dim = 2;
        z->regs[nr].type = REG_BNDS;
        for (nn = 0; nn < 6; nn++)
            z->regs[nr].data[nn] = range[nn];
        z->regs[nr].data[0] = dim[0];
        z->regs[nr].nfaces = nf;
        z->regs[nr].faces = (Face **) MALLOC (funcname,
                            (size_t)nf * sizeof(Face *));
        for (n = 0; n < nf; n++)
            z->regs[nr].faces[n] = copy_face (funcname, z->regs->faces[fn++]);
        nr++;

        nf = ni * nk;
        strcpy (z->regs[nr].name, "<jmin>");
        z->regs[nr].dim = 2;
        z->regs[nr].type = REG_BNDS;
        for (nn = 0; nn < 6; nn++)
            z->regs[nr].data[nn] = range[nn];
        z->regs[nr].data[4] = 1;
        z->regs[nr].nfaces = nf;
        z->regs[nr].faces = (Face **) MALLOC (funcname,
                            (size_t)nf * sizeof(Face *));
        for (n = 0; n < nf; n++)
            z->regs[nr].faces[n] = copy_face (funcname, z->regs->faces[fn++]);
        nr++;

        strcpy (z->regs[nr].name, "<jmax>");
        z->regs[nr].dim = 2;
        z->regs[nr].type = REG_BNDS;
        for (nn = 0; nn < 6; nn++)
            z->regs[nr].data[nn] = range[nn];
        z->regs[nr].data[1] = dim[1];
        z->regs[nr].nfaces = nf;
        z->regs[nr].faces = (Face **) MALLOC (funcname,
                            (size_t)nf * sizeof(Face *));
        for (n = 0; n < nf; n++)
            z->regs[nr].faces[n] = copy_face (funcname, z->regs->faces[fn++]);
        nr++;

        nf = ni * nj;
        strcpy (z->regs[nr].name, "<kmin>");
        z->regs[nr].dim = 2;
        z->regs[nr].type = REG_BNDS;
        for (nn = 0; nn < 6; nn++)
            z->regs[nr].data[nn] = range[nn];
        z->regs[nr].data[5] = 1;
        z->regs[nr].nfaces = nf;
        z->regs[nr].faces = (Face **) MALLOC (funcname,
                            (size_t)nf * sizeof(Face *));
        for (n = 0; n < nf; n++)
            z->regs[nr].faces[n] = copy_face (funcname, z->regs->faces[fn++]);
        nr++;

        strcpy (z->regs[nr].name, "<kmax>");
        z->regs[nr].dim = 2;
        z->regs[nr].type = REG_BNDS;
        for (nn = 0; nn < 6; nn++)
            z->regs[nr].data[nn] = range[nn];
        z->regs[nr].data[2] = dim[2];
        z->regs[nr].nfaces = nf;
        z->regs[nr].faces = (Face **) MALLOC (funcname,
                            (size_t)nf * sizeof(Face *));
        for (n = 0; n < nf; n++)
            z->regs[nr].faces[n] = copy_face (funcname, z->regs->faces[fn++]);
        nr++;
#endif
    }

    /* 1 to 1 interfaces */

    for (nn = 1; nn <= nints; nn++) {
        if (cg_1to1_read (cgnsfn, cgnsbase, cgnszone, nn,
                name, d_name, range, d_range, transform)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        strcpy (z->regs[nr].name, name);
        z->regs[nr].type = REG_1TO1;
        z->regs[nr].data[0] = 2 * CellDim;
        for (ii = 0; ii < 2*CellDim; ii++)
            z->regs[nr].data[ii+1] = range[ii];
        strcpy (z->regs[nr].d_name, d_name);
        k = structured_range (&z->regs[nr], dim, range, CGNS_ENUMV(Vertex));
        z->regs[nr].dim = k;
        nr++;
    }

    /* general connectivities */

    for (nn = 1; nn <= nconns; nn++) {
        if (cg_conn_info (cgnsfn, cgnsbase, cgnszone, nn, name,
                &location, &type, &ptype, &npts, d_name, &d_ztype,
                &d_ptype, &datatype, &ndpts)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        pts = (cgsize_t *) MALLOC (funcname, (size_t)(3 * npts) * sizeof(cgsize_t));
        if (cg_conn_read_short (cgnsfn, cgnsbase, cgnszone, nn, pts)) {
            free (pts);
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        strcpy (z->regs[nr].name, name);
        z->regs[nr].type = REG_CONN;
        z->regs[nr].data[0] = type;
        z->regs[nr].data[1] = location;
        z->regs[nr].data[2] = ptype;
        z->regs[nr].data[3] = npts;
        if (ptype == CGNS_ENUMV(PointRange)) {
            z->regs[nr].data[3] = 2 * CellDim;
            for (ii = 0; ii < 2*CellDim; ii++)
                z->regs[nr].data[ii+4] = pts[ii];
        }
        strcpy (z->regs[nr].d_name, d_name);

        if (type == CGNS_ENUMV(Abutting1to1) || type == CGNS_ENUMV(Abutting)) {
            if (ptype == CGNS_ENUMV(PointRange))
                k = structured_range (&z->regs[nr], dim, pts, location);
            else if (ptype == CGNS_ENUMV(PointList))
                k = structured_list (z->regs, &z->regs[nr], dim, npts, pts, location);
            else {
                k = 0;
                strcpy (z->regs[nr].errmsg, "invalid point set type");
            }
        }
        else if (type == CGNS_ENUMV(Overset)) {
            k = 0;
            strcpy (z->regs[nr].errmsg, "Overset connectivity not implemented");
        }
        else {
            k = 0;
            strcpy (z->regs[nr].errmsg, "invalid connectivity type");
        }
        z->regs[nr].dim = k;
        free (pts);
        nr++;
    }

    /* holes */

    for (nn = 1; nn <= nholes; nn++) {
        if (cg_hole_info (cgnsfn, cgnsbase, cgnszone, nn, name,
                &location, &ptype, &nsets, &npts)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        pts = (cgsize_t *) MALLOC (funcname, (size_t)(3 * npts * nsets) * sizeof(cgsize_t));
        if (cg_hole_read (cgnsfn, cgnsbase, cgnszone, nn, pts)) {
            free (pts);
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        strcpy (z->regs[nr].name, name);
        z->regs[nr].type = REG_HOLE;
        z->regs[nr].data[0] = nsets;
        z->regs[nr].data[1] = location;
        z->regs[nr].data[2] = ptype;
        z->regs[nr].data[3] = npts;

        if (ptype == CGNS_ENUMV(PointRange)) {
            z->regs[nr].data[3] = 2 * CellDim;
            for (ii = 0; ii < 2*CellDim; ii++)
                z->regs[nr].data[ii+4] = pts[ii];
            z->regs[nr].dim = structured_range (&z->regs[nr], dim, pts, location);

            if (z->regs[nr].dim && nsets > 1) {
                Edge *edges = z->regs[nr].edges;
                Face **faces = z->regs[nr].faces;
                ne = z->regs[nr].nedges;
                nf = z->regs[nr].nfaces;
                for (ii = 1; ii < nsets; ii++) {
                    z->regs[nr].nedges = z->regs[nr].nfaces = 0;
                    k = structured_range (&z->regs[nr], dim,
                        &pts[ii*2*CellDim], location);
                    if (k && z->regs[nr].nedges) {
                        edges = (Edge *) REALLOC (funcname,
                            (ne + z->regs[nr].nedges) * sizeof(Edge), edges);
                        for (i = 0; i < z->regs[nr].nedges; i++) {
                            edges[ne].nodes[0] = z->regs[nr].edges[i].nodes[0];
                            edges[ne].nodes[1] = z->regs[nr].edges[i].nodes[1];
                            ne++;
                        }
                        free(z->regs[nr].edges);
                    }
                    if (k && z->regs[nr].nfaces) {
                        faces = (Face **) REALLOC (funcname,
                            (nf + z->regs[nr].nfaces) * sizeof(Face *), faces);
                        for (i = 0; i < z->regs[nr].nfaces; i++)
                            faces[nf++] = z->regs[nr].faces[i];
                        free(z->regs[nr].faces);
                    }
                }
                z->regs[nr].nedges = ne;
                z->regs[nr].edges = edges;
                z->regs[nr].nfaces = nf;
                z->regs[nr].faces = faces;
            }
        }
        else if (ptype == CGNS_ENUMV(PointList)) {
            z->regs[nr].dim = structured_list (z->regs, &z->regs[nr],
                              dim, npts, pts, location);
        }
        else {
            strcpy (z->regs[nr].errmsg, "invalid Point Set Type");
        }
        free(pts);
        nr++;
    }

    /* boundary conditions */

    for (nn = 1; nn <= nbocos; nn++) {
        if (cg_boco_info (cgnsfn, cgnsbase, cgnszone, nn, name,
                &bctype, &ptype, &npts, transform, &j, &datatype, &ii) ||
            cg_boco_gridlocation_read(cgnsfn, cgnsbase, cgnszone, nn,
                &location)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        pts = (cgsize_t *) MALLOC (funcname, (size_t)(3 * npts) * sizeof(cgsize_t));
        if (cg_boco_read (cgnsfn, cgnsbase, cgnszone, nn, pts, 0)) {
            free (pts);
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        strcpy (z->regs[nr].name, name);
        z->regs[nr].type = REG_BOCO;
        z->regs[nr].data[0] = bctype;
        z->regs[nr].data[1] = location;
        z->regs[nr].data[2] = ptype;
        z->regs[nr].data[3] = i;
        if (ptype == CGNS_ENUMV(PointRange) || ptype == CGNS_ENUMV(ElementRange)) {
            z->regs[nr].data[3] = 2 * CellDim;
            for (ii = 0; ii < 2*CellDim; ii++)
                z->regs[nr].data[ii+4] = pts[ii];
        }

        if (ptype == CGNS_ENUMV(PointRange) || ptype == CGNS_ENUMV(ElementRange))
            k = structured_range (&z->regs[nr], dim, pts, location);
        else if (ptype == CGNS_ENUMV(PointList) || ptype == CGNS_ENUMV(ElementList))
            k = structured_list (z->regs, &z->regs[nr], dim, npts, pts, location);
        else {
            k = 0;
            strcpy (z->regs[nr].errmsg, "invalid point set type");
        }
        z->regs[nr].dim = k;
        free (pts);
        nr++;
    }

    z->nregs = nr;

#ifndef NO_CUTTING_PLANE
    for (nr = 0; nr < z->nregs; nr++) {
        if (z->regs[nr].dim == 2 && z->regs[nr].nfaces) {
            nf = z->regs[nr].nfaces;
            z->regs[nr].elemtype = CGNS_ENUMV(QUAD_4);
            z->regs[nr].nelems = nf;
            z->regs[nr].elems = (cgsize_t *) MALLOC (funcname,
	        (size_t)(4 * nf) * sizeof(cgsize_t));
            for (n = 0, j = 0; j < nf; j++, f++) {
                f = z->regs[nr].faces[j];
                for (ii = 0; ii < 4; ii++)
                    z->regs[nr].elems[n++] = f->nodes[ii];
            }
        }
    }
#endif

    return 0;
}

/*====================================================================
 * unstructured grid regions
 *====================================================================*/

static int max_face_nodes = 0;
static cgsize_t *sort_face_nodes = NULL;

/*-------------------------------------------------------------------*/

static int compare_faces (void *v1, void *v2)
{
    Face *f1 = (Face *)v1;
    Face *f2 = (Face *)v2;
    int i, k;
    cgsize_t id, nn, *n1, *n2;

    if (f1->nnodes != f2->nnodes)
        return (f1->nnodes - f2->nnodes);

    if (f1->nnodes > max_face_nodes) {
        max_face_nodes += 10;
        sort_face_nodes = (cgsize_t *) REALLOC ("compare_faces",
            2 * max_face_nodes * sizeof(cgsize_t), sort_face_nodes);
    }
    n1 = sort_face_nodes;
    n2 = sort_face_nodes + max_face_nodes;

    for (i = 0; i < f1->nnodes; i++) {
        id = f1->nodes[i];
        for (k = 0; k < i; k++) {
            if (n1[k] > id) {
                nn = n1[k];
                n1[k] = id;
                id = nn;
            }
        }
        n1[i] = id;
    }
    for (i = 0; i < f2->nnodes; i++) {
        id = f2->nodes[i];
        for (k = 0; k < i; k++) {
            if (n2[k] > id) {
                nn = n2[k];
                n2[k] = id;
                id = nn;
            }
        }
        n2[i] = id;
    }

    for (i = 0; i < f1->nnodes; i++) {
        if (n1[i] != n2[i])
            return (int)(n1[i] - n2[i]);
    }
    return 0;
}

/*-------------------------------------------------------------------*/

static size_t hash_face (void *v)
{
    Face *f = (Face *)v;
    int n;
    size_t hash = 0;

    for (n = 0; n < f->nnodes; n++)
        hash += (size_t)f->nodes[n];
    return hash;
}

/*-------------------------------------------------------------------*/

static size_t get_faces (void *vf, void *vr)
{
    Face *f = (Face *)vf;
    Regn *r = (Regn *)vr;

    r->faces[r->nfaces] = f;
    (r->nfaces)++;
    return 1;
}

/*-------------------------------------------------------------------*/

static int compare_poly (void *v1, void *v2)
{
    PolyFace *p1 = (PolyFace *)v1;
    PolyFace *p2 = (PolyFace *)v2;

    return compare_faces(p1->face, p2->face);
}

/*-------------------------------------------------------------------*/

static size_t hash_poly (void *v)
{
    PolyFace *p = (PolyFace *)v;

    return hash_face(p->face);
}

/*-------------------------------------------------------------------*/

static cgsize_t nPolyFaces;

static size_t poly_faces (void *vp, void *vl)
{
    PolyFace *pf = (PolyFace *)vp;
    PolyFace **pl = (PolyFace **)vl;

    pl[nPolyFaces++] = pf;
    return 1;
}

/*-------------------------------------------------------------------*/

static int poly_sort (const void *v1, const void *v2)
{
    const PolyFace **p1 = (const PolyFace **)v1;
    const PolyFace **p2 = (const PolyFace **)v2;

    return (int)((*p1)->num - (*p2)->num);
}

/*-------------------------------------------------------------------*/

static Face *find_face (Zone *z, cgsize_t fnum)
{
    int nr;
    cgsize_t nf;

    for (nr = 0; nr < z->nregs; nr++) {
        if (z->regs[nr].type == REG_ELEM &&
            z->regs[nr].dim == 2 &&
            z->regs[nr].data[1] <= fnum &&
            z->regs[nr].data[2] >= fnum) {
            nf = fnum - z->regs[nr].data[1];
            return z->regs[nr].faces[nf];
        }
    }
    return NULL;
}

/*-------------------------------------------------------------------*/

static int element_dimension (CGNS_ENUMT(ElementType_t) elemtype)
{
    switch (elemtype) {
        case CGNS_ENUMV(NODE):
            return 0;
        case CGNS_ENUMV(BAR_2):
        case CGNS_ENUMV(BAR_3):
            return 1;
        case CGNS_ENUMV(TRI_3):
        case CGNS_ENUMV(TRI_6):
        case CGNS_ENUMV(QUAD_4):
        case CGNS_ENUMV(QUAD_8):
        case CGNS_ENUMV(QUAD_9):
        case CGNS_ENUMV(NGON_n):
            return 2;
        case CGNS_ENUMV(TETRA_4):
        case CGNS_ENUMV(TETRA_10):
        case CGNS_ENUMV(PYRA_5):
        case CGNS_ENUMV(PYRA_13):
        case CGNS_ENUMV(PYRA_14):
        case CGNS_ENUMV(PENTA_6):
        case CGNS_ENUMV(PENTA_15):
        case CGNS_ENUMV(PENTA_18):
        case CGNS_ENUMV(HEXA_8):
        case CGNS_ENUMV(HEXA_20):
        case CGNS_ENUMV(HEXA_27):
        case CGNS_ENUMV(NFACE_n):
            return 3;
        default:
            break;
    }
    return -1;
}

/*-------------------------------------------------------------------*/

static void edge_elements (Regn *r, cgsize_t *conn)
{
    int ip;
    cgsize_t istart, n, ne, nelems;

    istart = r->data[1];
    nelems = r->data[2] - istart + 1;
    cg_npe ((CGNS_ENUMT(ElementType_t))r->data[0], &ip);

    r->nedges = nelems;
    r->edges = (Edge *) MALLOC ("edge_elements", nelems * sizeof(Edge));
    for (n = 0, ne = 0; ne < nelems; ne++) {
        r->edges[ne].id = istart + ne;
        r->edges[ne].nodes[0] = conn[n] - 1;
        r->edges[ne].nodes[1] = conn[n+1] - 1;
        n += ip;
    }
}

/*-------------------------------------------------------------------*/

static void face_elements (Regn *r, cgsize_t *conn)
{
    int i, ip;
    cgsize_t ne, nn, istart, nelems;
    cgsize_t rind0, rind1;
    CGNS_ENUMT(ElementType_t) elemtype, type;
    static char *funcname = "face_elements";

    elemtype = type = (CGNS_ENUMT(ElementType_t))r->data[0];
    istart = r->data[1];
    nelems = r->data[2] - istart + 1;
    rind0  = r->data[3];
    rind1  = nelems - r->data[4];

    r->nfaces = nelems;
    r->faces  = (Face **) MALLOC (funcname, nelems * sizeof(Face *));

    for (nn = 0, ne = 0; ne < nelems; ne++) {
        if (elemtype == CGNS_ENUMV(MIXED))
            type = (CGNS_ENUMT(ElementType_t))conn[nn++];
        switch (type) {
            case CGNS_ENUMV(TRI_3):
            case CGNS_ENUMV(TRI_6):
                ip = 3;
                break;
            case CGNS_ENUMV(QUAD_4):
            case CGNS_ENUMV(QUAD_8):
            case CGNS_ENUMV(QUAD_9):
                ip = 4;
                break;
            case CGNS_ENUMV(NGON_n):
                ip = (int)conn[nn++];
                break;
            default:
                if (type < CGNS_ENUMV(NODE) || type > CGNS_ENUMV(HEXA_27))
                    FATAL ("face_elements:unknown element type");
                ip = 0;
                break;
        }
        if (ip) {
            r->faces[ne] = new_face (funcname, ip);
            r->faces[ne]->id = istart + ne;
            for (i = 0; i < ip; i++)
                r->faces[ne]->nodes[i] = conn[nn+i] - 1;
            if (ne < rind0 || ne >= rind1)
                r->faces[ne]->flags = 1;
        }
        else {
            r->faces[ne] = NULL;
        }
        if (type == CGNS_ENUMV(NGON_n))
            nn += ip;
        else {
            cg_npe (type, &i);
            nn += i;
        }
    }
}

/*-------------------------------------------------------------------*/

static void exterior_faces (Zone *z, Regn *r, cgsize_t *conn)
{
    int i, j, nf, ip, flag;
    cgsize_t ne, nn, istart, nelems;
    cgsize_t rind0, rind1;
    CGNS_ENUMT(ElementType_t) elemtype;
    HASH *facehash;
    Face *face, *pf;
    static char *funcname = "exterior_faces";

    elemtype = (CGNS_ENUMT(ElementType_t))r->data[0];
    istart = r->data[1];
    nelems = r->data[2] - istart + 1;
    rind0  = r->data[3];
    rind1  = nelems - r->data[4];

    facehash = HashCreate (nelems > 1024 ? (size_t)nelems / 3 : 127,
                           compare_faces, hash_face);
    if (NULL == facehash)
        FATAL ("exterior_faces:face hash table creation failed");

    if (elemtype == CGNS_ENUMV(NFACE_n)) {
        for (nn = 0, ne = 0; ne < nelems; ne++) {
            flag = (ne < rind0 || ne >= rind1) ? 1 : 0;
            nf = conn[nn++];
            for (j = 0; j < nf; j++) {
                face = find_face (z, abs(conn[nn + j]));
                if (face != NULL) {
                    pf = (Face *) HashFind (facehash, face);
                    if (NULL == pf) {
                        pf = copy_face (funcname, face);
                        pf->id = 0;
                        pf->flags = flag;
                        (void) HashAdd (facehash, pf);
                    }
                    else if (flag == pf->flags) {
                        HashDelete (facehash, pf);
                        free (pf);
                    }
                    else {
                        pf->flags = 0;
                    }
                }
            }
            nn += nf;
        }
    }
    else {
        CGNS_ENUMT(ElementType_t) type = elemtype;
        face = new_face(funcname, 4);
        for (nn = 0, ne = 0; ne < nelems; ne++) {
            flag = (ne < rind0 || ne >= rind1) ? 1 : 0;
            if (elemtype == CGNS_ENUMV(MIXED))
                type = (CGNS_ENUMT(ElementType_t))conn[nn++];
            switch (type) {
                case CGNS_ENUMV(TETRA_4):
                case CGNS_ENUMV(TETRA_10):
                    ip = 2;
                    nf = 4;
                    break;
                case CGNS_ENUMV(PYRA_5):
                case CGNS_ENUMV(PYRA_13):
                case CGNS_ENUMV(PYRA_14):
                    ip = 6;
                    nf = 5;
                    break;
                case CGNS_ENUMV(PENTA_6):
                case CGNS_ENUMV(PENTA_15):
                case CGNS_ENUMV(PENTA_18):
                    ip = 11;
                    nf = 5;
                    break;
                case CGNS_ENUMV(HEXA_8):
                case CGNS_ENUMV(HEXA_20):
                case CGNS_ENUMV(HEXA_27):
                    ip = 16;
                    nf = 6;
                    break;
                default:
                    if (type < CGNS_ENUMV(NODE) ||
                        type > CGNS_ENUMV(HEXA_27))
                        FATAL ("exterior_faces:unknown element type");
                    ip = 0;
                    nf = 0;
                    break;
            }
            for (j = 0; j < nf; j++) {
                face->nnodes = facenodes[ip+j][0];
                for (i = 0; i < face->nnodes; i++)
                    face->nodes[i] = conn[nn + facenodes[ip+j][i+1]] - 1;
                pf = (Face *) HashFind (facehash, face);
                if (NULL == pf) {
                    pf = copy_face (funcname, face);
                    pf->id = 0;
                    pf->flags = flag;
                    (void) HashAdd (facehash, pf);
                }
                else if (flag == pf->flags) {
                    HashDelete (facehash, pf);
                    free (pf);
                }
                else {
                    pf->flags = 0;
                }
            }
            cg_npe (type, &j);
            nn += j;
        }
        free(face);
    }

    r->nfaces = 0;
    ne = (cgsize_t) HashSize (facehash);
    if (ne) {
        r->faces = (Face **) MALLOC (funcname, ne * sizeof(Face *));
        HashList (facehash, get_faces, r);
    }
    else {
        strcpy (r->errmsg, "couldn't find any exterior faces");
    }
    HashDestroy (facehash, NULL);
}

/*-------------------------------------------------------------------*/

static void polyhedral_faces (Zone *z, Regn *r, cgsize_t *conn)
{
    int j, nf, flag;
    cgsize_t ne, nn, istart, nelems, nfaces, id;
    cgsize_t rind0, rind1;
    HASH *facehash;
    Face *face;
    PolyFace poly, *pf, **polylist;
    static char *funcname = "polyhedral_faces";

    istart = r->data[1];
    nelems = r->data[2] - istart + 1;
    rind0  = r->data[3];
    rind1  = nelems - r->data[4];

    facehash = HashCreate (nelems > 1024 ? (size_t)nelems / 3 : 127,
                           compare_poly, hash_poly);
    if (NULL == facehash)
        FATAL ("polyhedral_faces:face hash table creation failed");

    nfaces = 0;
    for (nn = 0, ne = 0; ne < nelems; ne++) {
        flag = (ne < rind0 || ne >= rind1) ? 1 : 0;
        nf = conn[nn++];
        for (j = 0; j < nf; j++) {
            id = conn[nn + j];
            face = find_face (z, abs(id));
            poly.face = face;
            pf = (PolyFace *) HashFind (facehash, &poly);
            if (NULL == pf) {
                pf = (PolyFace *)MALLOC(funcname, sizeof(PolyFace));
                pf->face = face;
                pf->num = ++nfaces;
                pf->flags = flag;
                (void) HashAdd (facehash, pf);
            }
            else {
                if ((pf->flags & 1) != flag)
                    pf->flags &= ~1;
                pf->flags |= 2;
            }
            conn[nn + j] = id < 0 ? -(pf->num) : pf->num;
        }
        nn += nf;
    }

    nfaces = (cgsize_t) HashSize (facehash);
    polylist = (PolyFace **) MALLOC (funcname, nfaces * sizeof(PolyFace *));
    nPolyFaces = 0;
    HashList (facehash, poly_faces, polylist);
    HashDestroy (facehash, NULL);

    qsort(polylist, nfaces, sizeof(PolyFace *), poly_sort);

    for (nn = 0, ne = 0; ne < nfaces; ne++) {
        if ((polylist[ne]->flags & 2) == 0) nn++;
    }
    r->nfaces = nn;
    r->faces = (Face **) MALLOC (funcname, nn * sizeof(Face *));
    for (nn = 0, ne = 0; ne < nfaces; ne++) {
        if ((polylist[ne]->flags & 2) == 0) {
            r->faces[nn] = copy_face(funcname, polylist[ne]->face);
            r->faces[nn]->id = 0;
            r->faces[nn]->flags = (polylist[ne]->flags & 1);
            nn++;
        }
    }

    r->npoly = nfaces;
    r->poly = (Face **) MALLOC (funcname, nfaces * sizeof(Face *));
    for (nn = 0; nn < nfaces; nn++) {
        face = polylist[nn]->face;
        face->flags |= (polylist[nn]->flags & 2);
        r->poly[nn] = face;
        free(polylist[nn]);
    }
    free(polylist);

    r->elemtype = CGNS_ENUMV(NFACE_n);
    r->nelems = nelems;
    r->elems = conn;
}

/*-------------------------------------------------------------------*/

static int sort_elemsets(const void *v1, const void *v2)
{
    return (((Regn *)v2)->dim - ((Regn *)v1)->dim);
}

/*-------------------------------------------------------------------*/

static cgsize_t unstructured_region (int nregs, Regn *regs, Regn *r,
                                     CGNS_ENUMT(PointSetType_t) ptype,
                                     cgsize_t nlist, cgsize_t *list)
{
    int nr, nn;
    cgsize_t nf, nfaces, maxfaces;
    Face **faces, *f;
    static char *funcname = "unstructured_region";

    if (ptype == CGNS_ENUMV(PointList) ||
        ptype == CGNS_ENUMV(ElementList)) {
        for (nf = 1; nf < nlist; nf++) {
            if (list[nf] < list[nf-1]) {
                qsort (list, (size_t)nlist, sizeof(cgsize_t), compare_ints);
                break;
            }
        }
        maxfaces = nlist;
    }
    else if (ptype == CGNS_ENUMV(PointRange) ||
             ptype == CGNS_ENUMV(ElementRange)) {
        if (list[0] > list[1]) {
            nf = list[0];
            list[0] = list[1];
            list[1] = nf;
        }
        maxfaces = list[1] - list[0] + 1;
    }
    else {
        strcpy (r->errmsg, "invalid point set type");
        maxfaces = 0;
    }

    if (maxfaces < 1) return 0;
    faces = (Face **)MALLOC(funcname, (size_t)maxfaces * sizeof(Face *));

    nfaces = 0;
    for (nr = 0; nr < nregs; nr++) {
        if (!regs[nr].nfaces) continue;
        for (nf = 0; nf < regs[nr].nfaces; nf++) {
            f = regs[nr].faces[nf];
            switch (ptype) {
                case CGNS_ENUMV(PointList):
                    if (f->id) continue;
                    for (nn = 0; nn < f->nnodes; nn++) {
                        if (!find_int (f->nodes[nn]+1, nlist, list))
                            break;
                    }
                    if (nn == f->nnodes) break;
                    continue;
                case CGNS_ENUMV(PointRange):
                    if (f->id) continue;
                    for (nn = 0; nn < f->nnodes; nn++) {
                        if (f->nodes[nn]+1 < list[0] ||
                            f->nodes[nn]+1 > list[1])
                            break;
                    }
                    if (nn == f->nnodes) break;
                    continue;
                case CGNS_ENUMV(ElementList):
                    if (f->id && find_int (f->id, nlist, list)) break;
                    continue;
                case CGNS_ENUMV(ElementRange):
                    if (f->id >= list[0] && f->id <= list[1]) break;
                    continue;
                default:
                    continue;
            }
            if (nfaces == maxfaces) {
                maxfaces += 100;
                faces = (Face **) REALLOC (funcname,
                    (size_t)maxfaces * sizeof(Face *), faces);
            }
            faces[nfaces++] = f;
        }
    }
    if (nfaces) {
        r->nfaces = nfaces;
        r->faces = (Face **)MALLOC(funcname, (size_t)nfaces * sizeof(Face *));
        for (nf = 0; nf < nfaces; nf++)
            r->faces[nf] = copy_face(funcname, faces[nf]);
    }
    else
        strcpy (r->errmsg, "couldn't find any exterior faces");
    free (faces);
    return nfaces;
}

/*-------------------------------------------------------------------*/

static int unstructured_zone (Tcl_Interp *interp)
{
    int i, ns, nb, ip, nr, haspoly, nsets;
    int nsect, nints, nconns, nholes, nbocos, nrmlindex[3];
    int transform[3], rind[2];
    cgsize_t is, ie, np, n, ne, nf;
    cgsize_t nelem, elemsize, *conn;
    cgsize_t range[6], d_range[6];
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(GridConnectivityType_t) type;
    CGNS_ENUMT(PointSetType_t) ptype, d_ptype;
    CGNS_ENUMT(ZoneType_t) d_ztype;
    CGNS_ENUMT(DataType_t) datatype;
    CGNS_ENUMT(BCType_t) bctype;
    CGNS_ENUMT(ElementType_t) elemtype;
    char name[33], d_name[33];
    Zone *z = &zones[cgnszone-1];
    static char *funcname = "unstructured_zone";
    static char *dspmsg = "finding exterior faces for";

    if (cg_nsections (cgnsfn, cgnsbase, cgnszone, &nsect) ||
        cg_n1to1 (cgnsfn, cgnsbase, cgnszone, &nints) ||
        cg_nconns (cgnsfn, cgnsbase, cgnszone, &nconns) ||
        cg_nholes (cgnsfn, cgnsbase, cgnszone, &nholes) ||
        cg_nbocos (cgnsfn, cgnsbase, cgnszone, &nbocos)) {
        Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
        return 1;
    }
    z->nregs = nsect + nints + nconns + nholes + nbocos;
    z->regs = (Regn *) MALLOC (funcname, z->nregs * sizeof(Regn));

    /* element sets */

    haspoly = 0;
    for (nr = 0, ns = 1; ns <= nsect; ns++, nr++) {
        if (cg_section_read (cgnsfn, cgnsbase, cgnszone, ns,
                name, &elemtype, &is, &ie, &nb, &ip) ||
            cg_ElementDataSize (cgnsfn, cgnsbase, cgnszone, ns, &elemsize)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        zone_message (dspmsg, name);
        strcpy (z->regs[nr].name, name);
        z->regs[nr].type = REG_ELEM;
        z->regs[nr].data[0] = elemtype;
        z->regs[nr].data[1] = is;
        z->regs[nr].data[2] = ie;

        if (cg_goto (cgnsfn, cgnsbase, "Zone_t", cgnszone,
                "Elements_t", ns, "end") ||
            cg_rind_read (rind)) {
            rind[0] = rind[1] = 0;
        }
        z->regs[nr].data[3] = rind[0];
        z->regs[nr].data[4] = rind[1];

        if (elemtype < CGNS_ENUMV(BAR_2) || elemtype > CGNS_ENUMV(NFACE_n)) {
            strcpy (z->regs[nr].errmsg, "invalid element type");
            continue;
        }

        /* do this after reading all the sections */

        if (elemtype == CGNS_ENUMV(NFACE_n)) {
            z->regs[nr].dim = 3;
            haspoly++;
            continue;
        }

        nelem = ie - is + 1;
        conn = (cgsize_t *) MALLOC (funcname, (size_t)elemsize * sizeof(cgsize_t));
        if (cg_elements_read (cgnsfn, cgnsbase, cgnszone, ns, conn, 0)) {
            free (conn);
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }

        /* check element indices */

        if (elemtype == CGNS_ENUMV(MIXED)) {
            int dim;
            CGNS_ENUMT(ElementType_t) type;
            z->regs[nr].dim = -1;
            for (n = 0, ne = 0; ne < nelem; ne++) {
                type = (CGNS_ENUMT(ElementType_t))conn[n++];
                if (cg_npe (type, &ip) || ip <= 0) {
                    strcpy(z->regs[nr].errmsg,
                        "unhandled element type found in MIXED");
                    break;
                }
                for (i = 0; i < ip; i++) {
                    if (conn[n] < 1 || conn[n] > z->nnodes) {
                        strcpy(z->regs[nr].errmsg, "invalid element index");
                        break;
                    }
                    n++;
                }
                if (i < ip) break;
                dim = element_dimension(type);
                if (z->regs[nr].dim < dim) z->regs[nr].dim = dim;
            }
        }
        else if (elemtype == CGNS_ENUMV(NGON_n)) {
            z->regs[nr].dim = 2;
            for (n = 0, ne = 0; ne < nelem; ne++) {
                ip = (int)conn[n++];
                for (i = 0; i < ip; i++) {
                    if (conn[n] < 1 || conn[n] > z->nnodes) {
                        strcpy(z->regs[nr].errmsg, "invalid element index");
                        break;
                    }
                    n++;
                }
                if (i < ip) break;
            }
        }
        else {
            z->regs[nr].dim = element_dimension(elemtype);
            cg_npe (elemtype, &ip);
            for (n = 0, ne = 0; ne < nelem; ne++) {
                for (i = 0; i < ip; i++) {
                    if (conn[n] < 1 || conn[n] > z->nnodes) {
                        strcpy(z->regs[nr].errmsg, "invalid element index");
                        break;
                    }
                    n++;
                }
                if (i < ip) break;
            }
        }
        if (ne == nelem && z->regs[nr].dim > 0) {

            if (z->regs[nr].dim == 1)
                edge_elements(&z->regs[nr], conn);
            else if (z->regs[nr].dim == 2)
                face_elements (&z->regs[nr], conn);
            else
                exterior_faces (z, &z->regs[nr], conn);

#ifndef NO_CUTTING_PLANE
            if (z->regs[nr].dim > 1) {
                z->regs[nr].elemtype = elemtype;
                z->regs[nr].nelems = nelem;
                z->regs[nr].elems = conn;
                /* fix element indexing */
                cg_npe (elemtype, &ip);
                for (n = 0, ne = 0; ne < nelem; ne++) {
                    if (elemtype == CGNS_ENUMT(MIXED)) {
                        nb = (int)conn[n++];
                        cg_npe ((CGNS_ENUMT(ElementType_t))nb, &ip);
                    }
                    else if (elemtype == CGNS_ENUMT(NGON_n)) {
                        ip = (int)conn[n++];
                    }
                    for (i = 0; i < ip; i++) {
                        (conn[n])--;
                        n++;
                    }
                }
            }
            else
#endif
                free (conn);
        }
    }

    /* process NFACE_n sections */

    if (haspoly) {
        for (ns = 0; ns < nsect; ns++) {
            if (z->regs[ns].data[0] != CGNS_ENUMV(NFACE_n)) continue;
            zone_message (dspmsg, z->regs[ns].name);
            nelem = z->regs[ns].data[2] - z->regs[ns].data[1] + 1;
            cg_ElementDataSize (cgnsfn, cgnsbase, cgnszone, ns+1, &elemsize);
            conn = (cgsize_t *) MALLOC (funcname, (size_t)elemsize * sizeof(cgsize_t));
            if (cg_elements_read (cgnsfn, cgnsbase, cgnszone, ns+1, conn, 0)) {
                free (conn);
                Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
                return 1;
            }

            /* check element indices */

            for (n = 0, ne = 0; ne < nelem; ne++) {
                ip = (int)conn[n++];
                for (i = 0; i < ip; i++) {
                    if (NULL == find_face(z, abs(conn[n++]))) {
                        strcpy(z->regs[ns].errmsg, "invalid face index");
                        break;
                    }
                }
                if (i < ip) break;
            }

            if (ne == nelem) {
#ifndef NO_CUTTING_PLANE
                polyhedral_faces (z, &z->regs[ns], conn);
#else
                exterior_faces (z, &z->regs[ns], conn);
                free (conn);
#endif
            }
        }
    }

    qsort(z->regs, nr, sizeof(Regn), sort_elemsets);

    /* 1to1 connectivities */

    for (ns = 1; ns <= nints; ns++) {
        if (cg_1to1_read (cgnsfn, cgnsbase, cgnszone, ns,
                name, d_name, range, d_range, transform)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        zone_message (dspmsg, name);
        strcpy (z->regs[nr].name, name);
        z->regs[nr].type = REG_1TO1;
        z->regs[nr].data[0] = 2;
        z->regs[nr].data[1] = range[0];
        z->regs[nr].data[2] = range[1];
        strcpy (z->regs[nr].d_name, d_name);
        if (unstructured_region (nsect, z->regs, &z->regs[nr],
                CGNS_ENUMV(PointRange), 2, range)) z->regs[nr].dim = 2;
        nr++;
    }

    /* general connectivities */

    for (ns = 1; ns <= nconns; ns++) {
        if (cg_conn_info (cgnsfn, cgnsbase, cgnszone, ns, name,
                &location, &type, &ptype, &np, d_name, &d_ztype,
                &d_ptype, &datatype, &ie)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        zone_message (dspmsg, name);
        conn = (cgsize_t *) MALLOC (funcname, (size_t)np * sizeof(cgsize_t));
        if (cg_conn_read_short (cgnsfn, cgnsbase, cgnszone, ns, conn)) {
            free (conn);
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        strcpy (z->regs[nr].name, name);
        z->regs[nr].type = REG_CONN;
        z->regs[nr].data[0] = type;
        z->regs[nr].data[1] = location;
        z->regs[nr].data[2] = ptype;
        z->regs[nr].data[3] = np;
        if (ptype == CGNS_ENUMT(PointRange)) {
            z->regs[nr].data[4] = conn[0];
            z->regs[nr].data[5] = conn[1];
        }
        strcpy (z->regs[nr].d_name, d_name);

        if (type != CGNS_ENUMV(Abutting) &&
            type != CGNS_ENUMV(Abutting1to1)) {
            strcpy(z->regs[nr].errmsg,
               "can only handle Abutting or Abutting1to1 currently");
        }
        else if (ptype != CGNS_ENUMV(PointList) &&
                 ptype != CGNS_ENUMV(PointRange)) {
            strcpy(z->regs[nr].errmsg,
               "point set type not PointList or PointRange");
        }
        else if (location != CGNS_ENUMV(Vertex) &&
                 location != CGNS_ENUMV(CellCenter) &&
                 location != CGNS_ENUMV(FaceCenter)) {
            strcpy(z->regs[nr].errmsg,
               "location not Vertex, CellCenter or FaceCenter");
        }
        else {
            if (location != CGNS_ENUMV(Vertex)) {
                ptype = (ptype == CGNS_ENUMV(PointRange) ?
                         CGNS_ENUMV(ElementRange) : CGNS_ENUMV(ElementList));
            }
            if (unstructured_region (nsect, z->regs, &z->regs[nr],
                    ptype, np, conn)) z->regs[nr].dim = 2;
        }
        free (conn);
        nr++;
    }

    /* holes */

    for (ns = 1; ns <= nholes; ns++) {
        if (cg_hole_info (cgnsfn, cgnsbase, cgnszone, ns, name,
                &location, &ptype, &nsets, &np)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        conn = (cgsize_t *) MALLOC (funcname, (size_t)(3 * np * nsets) * sizeof(cgsize_t));
        if (cg_hole_read (cgnsfn, cgnsbase, cgnszone, ns, conn)) {
            free (conn);
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        strcpy (z->regs[nr].name, name);
        z->regs[nr].type = REG_HOLE;
        z->regs[nr].data[0] = nsets;
        z->regs[nr].data[1] = location;
        z->regs[nr].data[2] = ptype;
        z->regs[nr].data[3] = np;

        if (ptype == CGNS_ENUMV(PointRange)) {
            z->regs[nr].data[4] = conn[0];
            z->regs[nr].data[5] = conn[1];
        }
        else if (ptype != CGNS_ENUMV(PointList)) {
            strcpy(z->regs[nr].errmsg,
               "point set type not PointList or PointRange");
        }
        else {
            if (location == CGNS_ENUMV(Vertex)) {
                d_ptype = ptype;
            }
            else {
                d_ptype = (ptype == CGNS_ENUMV(PointRange) ?
                         CGNS_ENUMV(ElementRange) : CGNS_ENUMV(ElementList));
            }
            if (unstructured_region (nsect, z->regs, &z->regs[nr],
                    d_ptype, np, conn)) z->regs[nr].dim = 2;
            if (z->regs[nr].dim && nsets > 1 &&
                ptype == CGNS_ENUMV(PointRange)) {
                Edge *edges = z->regs[nr].edges;
                Face **faces = z->regs[nr].faces;
                ne = z->regs[nr].nedges;
                nf = z->regs[nr].nfaces;
                for (ip = 1; ip < nsets; ip++) {
                    z->regs[nr].nedges = z->regs[nr].nfaces = 0;
                    is = unstructured_region (nsect, z->regs, &z->regs[nr],
                             d_ptype, np, &conn[ip*2]);
                    if (is && z->regs[nr].nedges) {
                        edges = (Edge *) REALLOC (funcname,
                            (ne + z->regs[nr].nedges) * sizeof(Edge), edges);
                        for (i = 0; i < z->regs[nr].nedges; i++) {
                            edges[ne].nodes[0] = z->regs[nr].edges[i].nodes[0];
                            edges[ne].nodes[1] = z->regs[nr].edges[i].nodes[1];
                            ne++;
                        }
                        free(z->regs[nr].edges);
                    }
                    if (is && z->regs[nr].nfaces) {
                        faces = (Face **) REALLOC (funcname,
                            (nf + z->regs[nr].nfaces) * sizeof(Face *), faces);
                        for (i = 0; i < z->regs[nr].nfaces; i++)
                            faces[nf++] = z->regs[nr].faces[i];
                        free(z->regs[nr].faces);
                    }
                }
                z->regs[nr].nedges = ne;
                z->regs[nr].edges = edges;
                z->regs[nr].nfaces = nf;
                z->regs[nr].faces = faces;
            }
        }
        free (conn);
        nr++;
    }

    /* boundary conditions */

    for (ns = 1; ns <= nbocos; ns++) {
        if (cg_boco_info (cgnsfn, cgnsbase, cgnszone, ns, name,
                &bctype, &ptype, &np, nrmlindex, &is, &datatype, &nb) ||
            cg_boco_gridlocation_read (cgnsfn, cgnsbase, cgnszone, ns,
                &location)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        zone_message (dspmsg, name);
        conn = (cgsize_t *) MALLOC (funcname, (size_t)np * sizeof(cgsize_t));
        if (cg_boco_read (cgnsfn, cgnsbase, cgnszone, ns, conn, 0)) {
            free (conn);
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return 1;
        }
        strcpy (z->regs[nr].name, name);
        z->regs[nr].type = REG_BOCO;
        z->regs[nr].data[0] = bctype;
        z->regs[nr].data[1] = location;
        z->regs[nr].data[2] = ptype;
        z->regs[nr].data[3] = np;
        if (ptype == CGNS_ENUMV(PointRange) ||
            ptype == CGNS_ENUMV(ElementRange)) {
            z->regs[nr].data[4] = conn[0];
            z->regs[nr].data[5] = conn[1];
        }

        if ((ptype == CGNS_ENUMV(PointRange) ||
             ptype == CGNS_ENUMV(PointList)) &&
            (location == CGNS_ENUMV(FaceCenter) ||
             location == CGNS_ENUMV(CellCenter))) {
            ptype = (ptype == CGNS_ENUMV(PointRange) ?
                     CGNS_ENUMV(ElementRange) : CGNS_ENUMV(ElementList));
        }
        if (unstructured_region (nsect, z->regs, &z->regs[nr],
            ptype, np, conn)) z->regs[nr].dim = 2;
        free (conn);
        nr++;
    }

    z->nregs = nr;
    return 0;
}

/*====================================================================
 * find region edges
 *====================================================================*/

/*-------------------------------------------------------------------*/

static int compare_edges (void *v1, void *v2)
{
    Edge *e1 = (Edge *)v1;
    Edge *e2 = (Edge *)v2;

    if (e1->nodes[0] < e1->nodes[1]) {
        if (e2->nodes[0] < e2->nodes[1])
            return (int)(e1->nodes[0] - e2->nodes[0]);
        return (int)(e1->nodes[0] - e2->nodes[1]);
    }
    if (e2->nodes[0] < e2->nodes[1])
        return (int)(e1->nodes[1] - e2->nodes[0]);
    return (int)(e1->nodes[1] - e2->nodes[1]);
}

/*-------------------------------------------------------------------*/

static size_t hash_edge (void *v)
{
    Edge *e = (Edge *)v;

    return ((size_t)(e->nodes[0] + e->nodes[1]));
}

/*-------------------------------------------------------------------*/

static size_t get_edges (void *ve, void *vr)
{
    Edge *e = (Edge *)ve;
    Regn *r = (Regn *)vr;

    r->edges[r->nedges].nodes[0] = e->nodes[0];
    r->edges[r->nedges].nodes[1] = e->nodes[1];
    (r->nedges)++;
    return 1;
}

/*-------------------------------------------------------------------*/

static void extract_edges (Regn *r)
{
    int i, k;
    cgsize_t j, n;
    size_t ne;
    Face *f;
    Edge edge, *ep;
    HASH *edgehash;
    float dot;
    static char *funcname = "extract_edges";

    if (!r->nfaces) return;
    edgehash = HashCreate ((size_t)r->nfaces, compare_edges, hash_edge);
    if (NULL == edgehash)
        FATAL ("edge hash table creation failed");
    for (j = 0; j < r->nfaces; j++) {
        f = r->faces[j];
        if (f->flags) continue;
        for (i = 0, k = f->nnodes-1; i < f->nnodes; k = i++) {
            if (f->nodes[i] == f->nodes[k]) continue;
            if (f->nodes[i] < f->nodes[k]) {
                edge.nodes[0] = f->nodes[i];
                edge.nodes[1] = f->nodes[k];
            }
            else {
                edge.nodes[0] = f->nodes[k];
                edge.nodes[1] = f->nodes[i];
            }
            ep = (Edge *) HashFind (edgehash, &edge);
            if (NULL == ep) {
                ep = (Edge *) MALLOC (funcname, sizeof(Edge));
                ep->nodes[0] = edge.nodes[0];
                ep->nodes[1] = edge.nodes[1];
                ep->id = j;
                (void) HashAdd (edgehash, ep);
            }
            else {
                n = ep->id;
                dot = r->faces[n]->normal[0] * f->normal[0] +
                      r->faces[n]->normal[1] * f->normal[1] +
                      r->faces[n]->normal[2] * f->normal[2];
                if (dot > EDGE_ANGLE) {
                    HashDelete (edgehash, ep);
                    free (ep);
                }
            }
        }
    }

    ne = HashSize (edgehash);
    if (ne) {
        r->nedges = 0;
        r->edges = (Edge *) MALLOC (funcname, ne * sizeof(Edge));
        HashList (edgehash, get_edges, r);
    }
    HashDestroy (edgehash, NULL);
}

/*===================================================================
 * region manipulation
 *===================================================================*/

/*-------------------------------------------------------------------*/

static float *compute_normal (Node n0, Node n1, Node n2, Node n3)
{
    int j;
    double xn, yn, zn, sn;
    double d1[3], d2[3];
    static float norm[3];

    /* triangle */

    if (NULL == n3) {
        for (j = 0; j < 3; j++) {
            d1[j] = n1[j] - n0[j];
            d2[j] = n2[j] - n0[j];
        }
        sn = 0.5;
    }

    /* quadrilateral */

    else {
        for (j = 0; j < 3; j++) {
            d1[j] = n2[j] - n0[j];
            d2[j] = n3[j] - n1[j];
        }
        sn = 1.0;
    }
    xn = sn * (d1[1] * d2[2] - d2[1] * d1[2]);
    yn = sn * (d1[2] * d2[0] - d2[2] * d1[0]);
    zn = sn * (d1[0] * d2[1] - d2[0] * d1[1]);
    sn = sqrt (xn * xn + yn * yn + zn * zn);
    if (sn == 0.0) sn = 1.0;
    norm[0] = (float)(xn / sn);
    norm[1] = (float)(yn / sn);
    norm[2] = (float)(zn / sn);
    return norm;
}

/*-------------------------------------------------------------------*/

static float *face_normal (Zone *z, int nnodes, cgsize_t *nodes)
{
    int i, n;
    float *n0, *n1, *n2, *n3;
    float *norm, sn;
    static float sum[3];

    if (nnodes < 3) {
        for (i = 0; i < 3; i++)
            sum[i] = 0.0;
        return sum;
    }
    if (nnodes <= 4) {
        n0 = z->nodes[nodes[0]];
        n1 = z->nodes[nodes[1]];
        n2 = z->nodes[nodes[2]];
        n3 = nnodes == 4 ? z->nodes[nodes[3]] : NULL;
        return compute_normal(n0, n1, n2, n3);
    }

    for (i = 0; i < 3; i++)
        sum[i] = 0.0;
    n0 = z->nodes[nodes[0]];
    n1 = z->nodes[nodes[1]];
    for (n = 2; n < nnodes; n++) {
        n2 = z->nodes[nodes[n]];
        norm = compute_normal(n0, n1, n2, NULL);
        for (i = 0; i < 3; i++)
            sum[i] += norm[i];
        n1 = n2;
    }
    sn = (float)sqrt(sum[0]*sum[0] + sum[1]*sum[1] + sum[2]*sum[2]);
    if (sn == 0.0) sn = 1.0;
    for (i = 0; i < 3; i++)
        sum[i] /= sn;
    return sum;
}

/*-------------------------------------------------------------------*/

static void region_normals (Zone *z, Regn *r)
{
    int i, n;
    Face *f;
    float *norm;

    for (n = 0; n < r->nfaces; n++) {
        f = r->faces[n];
        norm = face_normal(z, f->nnodes, f->nodes);
        for (i = 0; i < 3; i++)
            f->normal[i] = norm[i];
    }
}

/*-------------------------------------------------------------------*/

static void bounding_box (Zone *z, Regn *r)
{
    int i, j, n;

    if (r->nfaces) {
        Face *f = r->faces[0];
        for (j = 0; j < 3; j++)
            r->bbox[j][0] = r->bbox[j][1] = z->nodes[f->nodes[0]][j];
        for (n = 0; n < r->nfaces; n++) {
            f = r->faces[n];
            for (i = 0; i < f->nnodes; i++) {
                for (j = 0; j < 3; j++) {
                    if (r->bbox[j][0] > z->nodes[f->nodes[i]][j])
                        r->bbox[j][0] = z->nodes[f->nodes[i]][j];
                    if (r->bbox[j][1] < z->nodes[f->nodes[i]][j])
                        r->bbox[j][1] = z->nodes[f->nodes[i]][j];
                }
            }
        }
    }
    else if (r->nedges) {
        Edge *e = r->edges;
        for (j = 0; j < 3; j++)
            r->bbox[j][0] = r->bbox[j][1] = z->nodes[e->nodes[0]][j];
        for (n = 0; n < r->nedges; n++, e++) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < 3; j++) {
                    if (r->bbox[j][0] > z->nodes[e->nodes[i]][j])
                        r->bbox[j][0] = z->nodes[e->nodes[i]][j];
                    if (r->bbox[j][1] < z->nodes[e->nodes[i]][j])
                        r->bbox[j][1] = z->nodes[e->nodes[i]][j];
                }
            }
        }
    }
    else {
        for (j = 0; j < 3; j++)
            r->bbox[j][0] = r->bbox[j][1] = 0.0;
    }
}

/*-------------------------------------------------------------------*/

static void get_bounds (int all, float bbox[3][2])
{
    int nz, nr, n, first = 1;

    for (nz = 0; nz < nzones; nz++) {
        for (nr = 0; nr < zones[nz].nregs; nr++) {
            if (zones[nz].regs[nr].nfaces &&
               (all || zones[nz].regs[nr].mode)) {
                if (first) {
                    for (n = 0; n < 3; n++) {
                        bbox[n][0] = zones[nz].regs[nr].bbox[n][0];
                        bbox[n][1] = zones[nz].regs[nr].bbox[n][1];
                    }
                    first = 0;
                }
                else {
                    for (n = 0; n < 3; n++) {
                        if (bbox[n][0] > zones[nz].regs[nr].bbox[n][0])
                            bbox[n][0] = zones[nz].regs[nr].bbox[n][0];
                        if (bbox[n][1] < zones[nz].regs[nr].bbox[n][1])
                            bbox[n][1] = zones[nz].regs[nr].bbox[n][1];
                    }
                }
            }
        }
    }
    if (first) {
        for (n = 0; n < 3; n++) {
            bbox[n][0] = 0.0;
            bbox[n][1] = 1.0;
        }
    }
}

/*-------------------------------------------------------------------*/

static void draw_outlines (Zone *z, Regn *r)
{
    glDisable (GL_LIGHTING);
    glShadeModel (GL_FLAT);
    glBegin (GL_LINES);
    if (r->nedges) {
        int ne;
        for (ne = 0; ne < r->nedges; ne++) {
            glVertex3fv (z->nodes[r->edges[ne].nodes[0]]);
            glVertex3fv (z->nodes[r->edges[ne].nodes[1]]);
        }
    }
    else {
        glVertex3f (r->bbox[0][0], r->bbox[1][0], r->bbox[2][0]);
        glVertex3f (r->bbox[0][1], r->bbox[1][0], r->bbox[2][0]);
        glVertex3f (r->bbox[0][1], r->bbox[1][0], r->bbox[2][0]);
        glVertex3f (r->bbox[0][1], r->bbox[1][1], r->bbox[2][0]);
        glVertex3f (r->bbox[0][1], r->bbox[1][1], r->bbox[2][0]);
        glVertex3f (r->bbox[0][0], r->bbox[1][1], r->bbox[2][0]);
        glVertex3f (r->bbox[0][0], r->bbox[1][1], r->bbox[2][0]);
        glVertex3f (r->bbox[0][0], r->bbox[1][0], r->bbox[2][0]);
        glVertex3f (r->bbox[0][0], r->bbox[1][0], r->bbox[2][1]);
        glVertex3f (r->bbox[0][1], r->bbox[1][0], r->bbox[2][1]);
        glVertex3f (r->bbox[0][1], r->bbox[1][0], r->bbox[2][1]);
        glVertex3f (r->bbox[0][1], r->bbox[1][1], r->bbox[2][1]);
        glVertex3f (r->bbox[0][1], r->bbox[1][1], r->bbox[2][1]);
        glVertex3f (r->bbox[0][0], r->bbox[1][1], r->bbox[2][1]);
        glVertex3f (r->bbox[0][0], r->bbox[1][1], r->bbox[2][1]);
        glVertex3f (r->bbox[0][0], r->bbox[1][0], r->bbox[2][1]);
        glVertex3f (r->bbox[0][0], r->bbox[1][0], r->bbox[2][0]);
        glVertex3f (r->bbox[0][0], r->bbox[1][0], r->bbox[2][1]);
        glVertex3f (r->bbox[0][1], r->bbox[1][0], r->bbox[2][0]);
        glVertex3f (r->bbox[0][1], r->bbox[1][0], r->bbox[2][1]);
        glVertex3f (r->bbox[0][1], r->bbox[1][1], r->bbox[2][0]);
        glVertex3f (r->bbox[0][1], r->bbox[1][1], r->bbox[2][1]);
        glVertex3f (r->bbox[0][0], r->bbox[1][1], r->bbox[2][0]);
        glVertex3f (r->bbox[0][0], r->bbox[1][1], r->bbox[2][1]);
    }
    glEnd ();
}

/*-------------------------------------------------------------------*/

static void draw_mesh (Zone *z, Regn *r)
{
    int nf, nn;
    Face *f;

    glEnable (GL_LIGHTING);
    glShadeModel (GL_FLAT);
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    for (nf = 0; nf < r->nfaces; nf++) {
        f = r->faces[nf];
        if (f->flags || f->nnodes < 2) continue;
        if (f->nnodes == 2)
            glBegin (GL_LINES);
        else if (f->nnodes == 3)
            glBegin (GL_TRIANGLES);
        else if (f->nnodes == 4)
            glBegin (GL_QUADS);
        else
            glBegin (GL_POLYGON);
        glNormal3fv (f->normal);
        for (nn = 0; nn < f->nnodes; nn++)
            glVertex3fv (z->nodes[f->nodes[nn]]);
        glEnd ();
    }
}

/*-------------------------------------------------------------------*/

static void draw_shaded (Zone *z, Regn *r)
{
    int nf, nn;
    Face *f;

    glEnable (GL_LIGHTING);
    glShadeModel (GL_FLAT);
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    for (nf = 0; nf < r->nfaces; nf++) {
        f = r->faces[nf];
        if (f->flags || f->nnodes < 3) continue;
        if (f->nnodes == 3)
            glBegin (GL_TRIANGLES);
        else if (f->nnodes == 4)
            glBegin (GL_QUADS);
        else
            glBegin (GL_POLYGON);
        glNormal3fv (f->normal);
        for (nn = 0; nn < f->nnodes; nn++)
            glVertex3fv (z->nodes[f->nodes[nn]]);
        glEnd ();
    }
}

/*===================================================================
 *           tcl interface
 *===================================================================*/

/*---------- CGNSopen ----------------------------------------------
 * open a CGNS file - return bases
 *------------------------------------------------------------------*/

static int CGNSopen (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int fn, nb, idum;
    static char buff[33];

    if (argc != 2) {
        Tcl_SetResult (interp, "usage: CGNSopen filename", TCL_STATIC);
        return TCL_ERROR;
    }
    Tcl_ResetResult (interp);

    if (cgnsfn) {
        cg_close (cgnsfn);
        cgnsfn = 0;
    }
    if (cg_open (argv[1], CG_MODE_READ, &fn)) {
        Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
        return TCL_ERROR;
    }
    if (cg_nbases (fn, &nb)) {
        Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
        cg_close (fn);
        return TCL_ERROR;
    }
    if (nb < 1) {
        Tcl_SetResult (interp, "no bases defined", TCL_STATIC);
        cg_close (fn);
        return TCL_ERROR;
    }
    nbases = nb;
    for (nb = 1; nb <= nbases; nb++) {
        if (cg_base_read (fn, nb, buff, &idum, &idum)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            cg_close (fn);
            return TCL_ERROR;
        }
        Tcl_AppendElement (interp, buff);
    }
    cgnsfn = fn;
    return TCL_OK;
}

/*---------- CGNSclose ---------------------------------------------
 * close the open CGNS file
 *------------------------------------------------------------------*/

static int CGNSclose (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    if (cgnsfn && cg_close (cgnsfn)) {
        Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
        return TCL_ERROR;
    }
    cgnsfn = 0;
    free_all ();
    Tcl_ResetResult (interp);
    return TCL_OK;
}

/*---------- CGNSbase ----------------------------------------------
 * set the CGNS base - return zones
 *------------------------------------------------------------------*/

static int CGNSbase (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int base, nz;
    cgsize_t sizes[9];
    char buff[33];

    if (!cgnsfn) {
        Tcl_SetResult (interp, "CGNS file not open", TCL_STATIC);
        return TCL_ERROR;
    }
    if (argc != 2) {
        Tcl_SetResult (interp, "usage: CGNSbase basenum", TCL_STATIC);
        return TCL_ERROR;
    }
    base = atoi (argv[1]) + 1;
    if (base < 1 || base > nbases) {
        Tcl_SetResult (interp, "base number out of range", TCL_STATIC);
        return TCL_ERROR;
    }
    Tcl_ResetResult (interp);
    cgnsbase = base;
    if (cg_base_read (cgnsfn, cgnsbase, BaseName, &CellDim, &PhyDim) ||
        cg_nzones (cgnsfn, cgnsbase, &nz)) {
        Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
        return TCL_ERROR;
    }
    free_all ();
    if (CellDim < 2 || CellDim > 3 || PhyDim < 2 || PhyDim > 3) {
        Tcl_SetResult (interp, "CellDim and Phydim not 2 or 3", TCL_STATIC);
        return TCL_ERROR;
    }
    nzones = nz;
    zones = (Zone *) MALLOC ("CGNSbase", nzones * sizeof(Zone));
    for (nz = 1; nz <= nzones; nz++) {
        if (cg_zone_read (cgnsfn, cgnsbase, nz, buff, sizes)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return TCL_ERROR;
        }
        strcpy (zones[nz-1].name, buff);
        Tcl_AppendElement (interp, buff);
    }
    return TCL_OK;
}

/*---------- CGNSzone ----------------------------------------------
 * set the CGNS zone - return regions
 *------------------------------------------------------------------*/

static int CGNSzone (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int i, j, zone, nr, ncoords;
    cgsize_t sizes[9], rng[2][3], n, nnodes;
    int rind[6];
    CGNS_ENUMT(DataType_t) datatype;
    CGNS_ENUMT(ZoneType_t) zonetype;
    Node *nodes;
    float *xyz;
    double rad, theta, phi;
    char buff[65], coordtype[4];
    Zone *z;
#ifdef NO_CUTTING_PLANE
    int *tag, nf, nn;
#endif
    static char *funcname = "CGNSzone";

    if (!cgnsfn) {
        Tcl_SetResult (interp, "CGNS file not open", TCL_STATIC);
        return TCL_ERROR;
    }
    if (argc != 2) {
        Tcl_SetResult (interp, "usage: CGNSzone zonenum", TCL_STATIC);
        return TCL_ERROR;
    }
    zone = atoi (argv[1]) + 1;
    if (zone < 1 || zone > nzones) {
        Tcl_SetResult (interp, "zone number out of range", TCL_STATIC);
        return TCL_ERROR;
    }

    if (cg_zone_read (cgnsfn, cgnsbase, zone, buff, sizes) ||
        cg_zone_type (cgnsfn, cgnsbase, zone, &zonetype)) {
        Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
        return TCL_ERROR;
    }
    if (zonetype == CGNS_ENUMV(Structured)) {
        for (i = 0; i < CellDim; i++) {
            if (sizes[i] < 2) {
                Tcl_SetResult (interp, "zone dimension < 2", TCL_STATIC);
                return TCL_ERROR;
            }
        }
    }
    else if (zonetype != CGNS_ENUMV(Unstructured)) {
        Tcl_SetResult (interp, "invalid zone type", TCL_STATIC);
        return TCL_ERROR;
    }
    cgnszone = zone;
    z = &zones[zone-1];

    /* get number of coordinates */

    if (cg_ncoords (cgnsfn, cgnsbase, cgnszone, &ncoords)) {
        Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
        return TCL_ERROR;
    }
    if (ncoords < PhyDim) {
        Tcl_SetResult (interp, "less than PhyDim coordinates", TCL_STATIC);
        return TCL_ERROR;
    }

    /* check for rind */

    if (cg_goto (cgnsfn, cgnsbase, "Zone_t", zone,
        "GridCoordinates_t", 1, "end")) {
        Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
        return TCL_ERROR;
    }
    if ((i = cg_rind_read (rind)) != CG_OK) {
        if (i != CG_NODE_NOT_FOUND) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return TCL_ERROR;
        }
        for (i = 0; i < 6; i++)
            rind[i] = 0;
    }

    /* get grid coordinate range */

    if (zonetype == CGNS_ENUMV(Structured)) {
        for (i = 0; i < 3; i++) {
            rng[0][i] = 1;
            rng[1][i] = 1;
        }
        nnodes = 1;
        for (i = 0; i < CellDim; i++) {
            rng[0][i] = rind[2*i] + 1;
            rng[1][i] = rind[2*i] + sizes[i];
            nnodes *= sizes[i];
        }
    }
    else {
        nnodes = sizes[0] + rind[0] + rind[1];
        rng[0][0] = 1;
        rng[1][0] = nnodes;
    }

    /* read the nodes */

    strcpy (coordtype, "   ");
    zone_message ("reading coordinates", NULL);
    xyz = (float *) MALLOC (funcname, (size_t)nnodes * sizeof(float));
    nodes = (Node *) MALLOC (funcname, (size_t)nnodes * sizeof(Node));
    for (i = 1; i <= ncoords; i++) {
        if (cg_coord_info (cgnsfn, cgnsbase, cgnszone, i, &datatype, buff) ||
            cg_coord_read (cgnsfn, cgnsbase, cgnszone, buff,
                CGNS_ENUMV(RealSingle), rng[0], rng[1], xyz)) {
            free (xyz);
            free (nodes);
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return TCL_ERROR;
        }
        if (0 == strcmp (buff, "CoordinateX") ||
            0 == strcmp (buff, "CoordinateR"))
            j = 0;
        else if (0 == strcmp (buff, "CoordinateY") ||
                 0 == strcmp (buff, "CoordinateTheta"))
            j = 1;
        else if (0 == strcmp (buff, "CoordinateZ") ||
                 0 == strcmp (buff, "CoordinatePhi"))
            j = 2;
        else
            continue;
        if (coordtype[j] == ' ' || strchr ("XYZ", buff[10]) != NULL)
            coordtype[j] = buff[10];
        for (n = 0; n < nnodes; n++)
            nodes[n][j] = xyz[n];
    }
    free (xyz);
    if (0 == strncmp (coordtype, "RTZ", PhyDim)) {
        for (n = 0; n < nnodes; n++) {
            rad = nodes[n][0];
            theta = nodes[n][1];
            nodes[n][0] = (float)(rad * cos (theta));
            nodes[n][1] = (float)(rad * sin (theta));
        }
    }
    else if (0 == strcmp (coordtype, "RTP")) {
        for (n = 0; n < nnodes; n++) {
            rad = nodes[n][0];
            theta = nodes[n][1];
            phi = nodes[n][2];
            nodes[n][0] = (float)(rad * sin (theta) * cos (phi));
            nodes[n][1] = (float)(rad * sin (theta) * sin (phi));
            nodes[n][2] = (float)(rad * cos (theta));
        }
    }
    else if (strncmp (coordtype, "XYZ", PhyDim)) {
        free (nodes);
        Tcl_SetResult (interp, "unknown coordinate system", TCL_STATIC);
        return TCL_ERROR;
    }

    z->nnodes = nnodes;
    z->nodes = nodes;

    /* build regions */

    if (zonetype == CGNS_ENUMV(Structured)) {
        if (structured_zone (interp, sizes))
            return TCL_ERROR;
    }
    else {
        if (unstructured_zone (interp))
            return TCL_ERROR;
    }

#ifdef NO_CUTTING_PLANE

    tag = (int *) MALLOC (funcname, nnodes * sizeof(int));
    for (n = 0; n < nnodes; n++)
        tag[n] = -1;

    /* tag nodes which are actually used */

    for (nn = 0, nr = 0; nr < z->nregs; nr++) {
        for (nf = 0; nf < z->regs[nr].nfaces; nf++) {
            for (n = 0; n < z->regs[nr].faces[nf]->nnodes; n++) {
                i = z->regs[nr].faces[nf]->nodes[n];
                if (tag[i] < 0)
                    tag[i] = nn++;
            }
        }
    }

    nodes = (Node *) MALLOC (funcname, nn * sizeof(Node));
    for (n = 0; n < nnodes; n++) {
        if (tag[n] >= 0) {
            j = tag[n];
            for (i = 0; i < 3; i++)
                nodes[j][i] = z->nodes[n][i];
        }
    }

    free(z->nodes);
    z->nodes = nodes;
    z->nnodes = nn;

    /* re-index region faces */

    for (nr = 0; nr < z->nregs; nr++) {
        for (nf = 0; nf < z->regs[nr].nfaces; nf++) {
            for (n = 0; n < z->regs[nr].faces[nf]->nnodes; n++) {
                i = z->regs[nr].faces[nf]->nodes[n];
                z->regs[nr].faces[nf]->nodes[n] = tag[i];
            }
        }
    }

    free(tag);

#endif

    /* find region bounding box, edges and normals */

    zone_message ("finding normals and edges", NULL);
    for (nr = 0; nr < z->nregs; nr++) {
        if (z->regs[nr].nfaces) {
            bounding_box (z, &z->regs[nr]);
            region_normals (z, &z->regs[nr]);
            extract_edges (&z->regs[nr]);
        }
    }

    Tcl_ResetResult (interp);
    for (nr = 0; nr < z->nregs; nr++) {
        switch (z->regs[nr].type) {
            case REG_MESH:
                strcpy(buff, z->regs[nr].name);
                break;
            case REG_ELEM:
                sprintf(buff, "<Element Sections>/%s", z->regs[nr].name);
                break;
            case REG_1TO1:
                sprintf(buff, "<1to1 Connections>/%s", z->regs[nr].name);
                break;
            case REG_CONN:
                sprintf(buff, "<General Connections>/%s", z->regs[nr].name);
                break;
            case REG_HOLE:
                sprintf(buff, "<Overset Holes>/%s", z->regs[nr].name);
                break;
            case REG_BOCO:
                sprintf(buff, "<Boundary Conditions>/%s", z->regs[nr].name);
                break;
            case REG_BNDS:
                sprintf(buff, "<Mesh Boundaries>/%s", z->regs[nr].name);
                break;
        }
        Tcl_AppendElement (interp, buff);
    }
    return TCL_OK;
}

/*---------- CGNSsummary -------------------------------------------
 * return info summary string
 *------------------------------------------------------------------*/

static int CGNSsummary (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int n, nz;
    char *p, buff[128];
    Regn *r;

    if (!cgnsfn) {
        Tcl_SetResult (interp, "CGNS file not open", TCL_STATIC);
        return TCL_ERROR;
    }
    if (argc < 1 || argc > 3) {
        Tcl_SetResult (interp, "usage: CGNSsummary [zone [reg]]", TCL_STATIC);
        return TCL_ERROR;
    }
    Tcl_ResetResult (interp);
    if (argc == 1) {
        sprintf (buff, "Physical Dim = %d, Cell Dim = %d", PhyDim, CellDim);
        Tcl_AppendResult (interp, buff, NULL);
        return TCL_OK;
    }

    nz = atoi (argv[1]);
    if (nz < 0 || nz >= nzones) {
        Tcl_SetResult (interp, "zone number out of range", TCL_STATIC);
        return TCL_ERROR;
    }

    if (argc == 2) {
        cgsize_t sizes[9];
        CGNS_ENUMT(ZoneType_t) zonetype;
        if (cg_zone_read (cgnsfn, cgnsbase, nz+1, buff, sizes) ||
            cg_zone_type (cgnsfn, cgnsbase, nz+1, &zonetype)) {
            Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
            return TCL_ERROR;
        }
        Tcl_AppendResult (interp, cg_ZoneTypeName(zonetype),
            " Zone : ", NULL);
        if (zonetype == CGNS_ENUMV(Unstructured)) {
            sprintf (buff, "%ld vertices, %ld elements",
                (long)sizes[0], (long)sizes[1]);
            Tcl_AppendResult (interp, buff, NULL);
        }
        else {
            sprintf (buff, "%ld", (long)sizes[0]);
            for (n = 1; n < CellDim; n++) {
                p = buff + strlen(buff);
                sprintf (p, " x %ld", (long)sizes[n]);
            }
            Tcl_AppendResult (interp, buff, " vertices", NULL);
        }
        return TCL_OK;
    }

    n = atoi (argv[2]);
    if (n < 0 || n >= zones[nz].nregs) {
        Tcl_SetResult (interp, "region number out of range", TCL_STATIC);
        return TCL_ERROR;
    }
    r = &zones[nz].regs[n];

    switch (r->type) {
        case REG_MESH:
            if (CellDim == 2) {
                sprintf (buff, "%ld x %ld",
                    (long)r->data[0], (long)r->data[1]);
            }
            else {
                sprintf (buff, "%ld x %ld x %ld",
                    (long)r->data[0], (long)r->data[1], (long)r->data[2]);
            }
            Tcl_AppendResult (interp, "Structured Mesh : ",
                buff, " vertices", NULL);
            break;
        case REG_ELEM:
            sprintf (buff, "%ld", (long)(r->data[2] - r->data[1] + 1));
            Tcl_AppendResult (interp, cg_ElementTypeName(r->data[0]),
                " Element Set : ", buff, " elements", NULL);
            break;
        case REG_1TO1:
            if (r->data[0] == 2)
                sprintf (buff, "%ld", (long)(r->data[2] - r->data[1] + 1));
            else if (CellDim == 2) {
                sprintf (buff, "%ld x %ld",
                    (long)(r->data[3] - r->data[1] + 1),
                    (long)(r->data[4] - r->data[2] + 1));
            }
            else {
                sprintf (buff, "%ld x %ld x %ld",
                    (long)(r->data[4] - r->data[1] + 1),
                    (long)(r->data[5] - r->data[2] + 1),
                    (long)(r->data[6] - r->data[3] + 1));
            }
            Tcl_AppendResult (interp, "1to1 Connection : PointRange ",
                buff, " -> ", r->d_name, NULL);
            break;
        case REG_CONN:
            if (r->data[2] == CGNS_ENUMV(PointList) ||
                r->data[2] == CGNS_ENUMV(ElementList))
                sprintf (buff, " %ld", (long)r->data[3]);
            else if (r->data[3] == 2)
                sprintf (buff, " %ld", (long)(r->data[5] - r->data[4] + 1));
            else if (CellDim == 2) {
                sprintf (buff, " %ld x %ld",
                    (long)(r->data[6] - r->data[4] + 1),
                    (long)(r->data[7] - r->data[5] + 1));
            }
            else {
                sprintf (buff, " %ld x %ld x %ld",
                    (long)(r->data[7] - r->data[4] + 1),
                    (long)(r->data[8] - r->data[5] + 1),
                    (long)(r->data[9] - r->data[6] + 1));
            }
            Tcl_AppendResult (interp,
                cg_GridConnectivityTypeName(r->data[0]),
                " Connection : ", cg_PointSetTypeName(r->data[2]),
                buff, " -> ", r->d_name, NULL);
            break;
        case REG_HOLE:
            if (r->data[2] == CGNS_ENUMV(PointList) ||
                r->data[2] == CGNS_ENUMV(ElementList))
                sprintf (buff, " %ld", (long)r->data[3]);
            else if (r->data[3] == 2)
                sprintf (buff, " %ld", (long)(r->data[5] - r->data[4] + 1));
            else if (CellDim == 2) {
                sprintf (buff, " %ld x %ld",
                    (long)(r->data[6] - r->data[4] + 1),
                    (long)(r->data[7] - r->data[5] + 1));
            }
            else {
                sprintf (buff, " %ld x %ld x %ld",
                    (long)(r->data[7] - r->data[4] + 1),
                    (long)(r->data[8] - r->data[5] + 1),
                    (long)(r->data[9] - r->data[6] + 1));
            }
            Tcl_AppendResult (interp, "Overset Hole : ",
                cg_PointSetTypeName(r->data[2]), buff, NULL);
            break;
        case REG_BOCO:
            if (r->data[2] == CGNS_ENUMV(PointList) ||
                r->data[2] == CGNS_ENUMV(ElementList))
                sprintf (buff, " %ld", (long)r->data[3]);
            else if (r->data[3] == 2)
                sprintf (buff, " %ld", (long)(r->data[5] - r->data[4] + 1));
            else if (CellDim == 2) {
                sprintf (buff, " %ld x %ld",
                    (long)(r->data[6] - r->data[4] + 1),
                    (long)(r->data[7] - r->data[5] + 1));
            }
            else {
                sprintf (buff, " %ld x %ld x %ld",
                    (long)(r->data[7] - r->data[4] + 1),
                    (long)(r->data[8] - r->data[5] + 1),
                    (long)(r->data[9] - r->data[6] + 1));
            }
            Tcl_AppendResult (interp, cg_BCTypeName(r->data[0]),
                " Boundary Condition : ", cg_PointSetTypeName(r->data[2]),
                buff, NULL);
            break;
        case REG_BNDS:
            if (CellDim == 2) {
                sprintf (buff, "%ld x %ld",
                    (long)(r->data[2] - r->data[0] + 1),
                    (long)(r->data[3] - r->data[1] + 1));
            }
            else {
                sprintf (buff, "%ld x %ld x %ld",
                    (long)(r->data[3] - r->data[0] + 1),
                    (long)(r->data[4] - r->data[1] + 1),
                    (long)(r->data[5] - r->data[2] + 1));
            }
            Tcl_AppendResult (interp, "Mesh Boundary : ", buff,
                " vertices", NULL);
            break;
    }
    return TCL_OK;
}

/*---------- CGNSgetbase -------------------------------------------
 * get base properties
 *------------------------------------------------------------------*/

static int CGNSgetbase (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    char cd[16], pd[16], nz[16];

    if (!cgnsfn) {
        Tcl_SetResult (interp, "CGNS file not open", TCL_STATIC);
        return TCL_ERROR;
    }
    if (argc != 1) {
        Tcl_SetResult (interp, "usage: CGNSgetbase", TCL_STATIC);
        return TCL_ERROR;
    }
    Tcl_ResetResult (interp);

    sprintf (pd, "%d", PhyDim);
    sprintf (cd, "%d", CellDim);
    sprintf (nz, "%d", nzones);
    Tcl_AppendResult (interp,
          "Base Name   : ", BaseName,
        "\nPhysical Dim: ", pd,
        "\nCell Dim    : ", cd,
        "\nNumber Zones: ", nz, NULL);
    return TCL_OK;
}

/*---------- CGNSgetzone -------------------------------------------
 * get zone properties
 *------------------------------------------------------------------*/

static int CGNSgetzone (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int n, ndim, zone, cnts[4];
    cgsize_t sizes[9];
    CGNS_ENUMT(ZoneType_t) zonetype;
    char *p, buff[65];
    static char *cntname[] = {
        "Element Sections",
        "1to1 Connections",
        "General Connections",
        "Boundary Conditions"
    };

    if (!cgnsfn) {
        Tcl_SetResult (interp, "CGNS file not open", TCL_STATIC);
        return TCL_ERROR;
    }
    if (argc != 2) {
        Tcl_SetResult (interp, "usage: CGNSgetzone zonenum", TCL_STATIC);
        return TCL_ERROR;
    }
    zone = atoi (argv[1]) + 1;
    if (zone < 1 || zone > nzones) {
        Tcl_SetResult (interp, "zone number out of range", TCL_STATIC);
        return TCL_ERROR;
    }

    if (cg_zone_read (cgnsfn, cgnsbase, zone, buff, sizes) ||
        cg_zone_type (cgnsfn, cgnsbase, zone, &zonetype) ||
        cg_nsections (cgnsfn, cgnsbase, zone, &cnts[0]) ||
        cg_n1to1 (cgnsfn, cgnsbase, zone, &cnts[1]) ||
        cg_nconns (cgnsfn, cgnsbase, zone, &cnts[2]) ||
        cg_nbocos (cgnsfn, cgnsbase, zone, &cnts[3])) {
        Tcl_SetResult (interp, (char *)cg_get_error(), TCL_STATIC);
        return TCL_ERROR;
    }
    Tcl_ResetResult (interp);

    Tcl_AppendResult (interp, "Zone Name          : ", buff,
        "\nType of Zone       : ", cg_ZoneTypeName(zonetype),
        "\nVertex Dimensions  : ", NULL);

    ndim = zonetype == CGNS_ENUMV(Unstructured) ? 1 : CellDim;
    sprintf (buff, "%ld", (long)sizes[0]);
    for (n = 1; n < ndim; n++) {
        p = buff + strlen(buff);
        sprintf (p, " x %ld", (long)sizes[n]);
    }
    Tcl_AppendResult (interp, buff, "\nCell Dimensions    : ", NULL);

    sprintf (buff, "%ld", (long)sizes[ndim]);
    for (n = 1; n < ndim; n++) {
        p = buff + strlen(buff);
        sprintf (p, " x %ld", (long)sizes[n+CellDim]);
    }
    Tcl_AppendResult (interp, buff, NULL);

    for (n = 0; n < 4; n++) {
        sprintf (buff, "\n%-19s: %d", cntname[n], cnts[n]);
        Tcl_AppendResult (interp, buff, NULL);
    }
    return TCL_OK;
}

/*---------- CGNSgetregion -----------------------------------------
 * get region properties
 *------------------------------------------------------------------*/

static int CGNSgetregion (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int n;
    char buff[128];
    Zone *z;
    Regn *r;

    if (!cgnsfn) {
        Tcl_SetResult (interp, "CGNS file not open", TCL_STATIC);
        return TCL_ERROR;
    }
    if (argc != 3) {
        Tcl_SetResult (interp, "usage: CGNSgetregion zone reg", TCL_STATIC);
        return TCL_ERROR;
    }
    n = atoi (argv[1]);
    if (n < 0 || n >= nzones) {
        Tcl_SetResult (interp, "zone number out of range", TCL_STATIC);
        return TCL_ERROR;
    }
    z = &zones[n];
    n = atoi (argv[2]);
    if (n < 0 || n >= z->nregs) {
        Tcl_SetResult (interp, "region number out of range", TCL_STATIC);
        return TCL_ERROR;
    }
    r = &z->regs[n];

    Tcl_ResetResult (interp);
    switch (r->type) {
        case REG_MESH:
            if (CellDim == 2) {
                sprintf (buff, "%ld x %ld",
                    (long)r->data[0], (long)r->data[1]);
            }
            else {
                sprintf (buff, "%ld x %ld x %ld",
                    (long)r->data[0], (long)r->data[1], (long)r->data[2]);
            }
            Tcl_AppendResult (interp,
                  "Region Name    : ", r->name,
                "\nType of Region : Structured Mesh",
                "\nMesh Dimensions: ", buff, NULL);
            break;
        case REG_ELEM:
            sprintf (buff, "%ld -> %ld",
                (long)r->data[1], (long)r->data[2]);
            Tcl_AppendResult (interp,
                  "Region Name     : ", r->name,
                "\nType of Region  : Element Set",
                "\nElement Set Type: ", cg_ElementTypeName(r->data[0]),
                "\nElement Range   : ", buff, NULL);
            if (r->data[3] || r->data[4]) {
                sprintf (buff, "%ld %ld",
                    (long)r->data[3], (long)r->data[4]);
                Tcl_AppendResult (interp,
                    "\nRind Elements   : ", buff, NULL);
            }
            break;
        case REG_1TO1:
            Tcl_AppendResult (interp,
                  "Region Name   : ", r->name,
                "\nType of Region: 1to1 Connectivity",
                "\nPoint Set Type: PointRange", NULL);
            if (r->data[0] == 2) {
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[1], (long)r->data[2]);
                Tcl_AppendResult (interp,
                    "\nIndex Range   : ", buff, NULL);
            }
            else {
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[1], (long)r->data[1+CellDim]);
                Tcl_AppendResult (interp,
                    "\nI Index Range : ", buff, NULL);
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[2], (long)r->data[2+CellDim]);
                Tcl_AppendResult (interp,
                    "\nJ Index Range : ", buff, NULL);
                if (CellDim == 3) {
                    sprintf (buff, "%ld -> %ld",
                        (long)r->data[3], (long)r->data[6]);
                    Tcl_AppendResult (interp,
                        "\nK Index Range : ", buff, NULL);
                }
            }
            Tcl_AppendResult (interp, "\nDonor Zone    : ", r->d_name, NULL);
            break;
        case REG_CONN:
            Tcl_AppendResult (interp,
                  "Region Name      : ", r->name,
                "\nType of Region   : General Connectivity",
                "\nConnectivity Type: ",
                    cg_GridConnectivityTypeName(r->data[0]),
                "\nGrid Location    : ", cg_GridLocationName(r->data[1]),
                "\nPoint Set Type   : ", cg_PointSetTypeName(r->data[2]),
                NULL);
            if (r->data[2] == CGNS_ENUMV(PointList) || r->data[2] == CGNS_ENUMV(ElementList)) {
                sprintf (buff, "%ld", (long)r->data[3]);
                Tcl_AppendResult (interp, "\nNumber of Points : ", buff, NULL);
            }
            else if (r->data[3] == 2) {
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[4], (long)r->data[5]);
                Tcl_AppendResult (interp, "\nIndex Range      : ", buff, NULL);
            }
            else {
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[4], (long)r->data[4+CellDim]);
                Tcl_AppendResult (interp,
                    "\nI Index Range    : ", buff, NULL);
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[5], (long)r->data[5+CellDim]);
                Tcl_AppendResult (interp,
                    "\nJ Index Range    : ", buff, NULL);
                if (CellDim == 3) {
                    sprintf (buff, "%ld -> %ld",
                        (long)r->data[6], (long)r->data[9]);
                    Tcl_AppendResult (interp,
                        "\nK Index Range    : ", buff, NULL);
                }
            }
            Tcl_AppendResult (interp, "\nDonor Zone       : ", r->d_name, NULL);
            break;
        case REG_HOLE:
            sprintf(buff, "%ld", (long)r->data[0]);
            Tcl_AppendResult (interp,
                  "Region Name      : ", r->name,
                "\nType of Region   : Overset Hole",
                "\nNumber Sets      : ", buff,
                "\nGrid Location    : ", cg_GridLocationName(r->data[1]),
                "\nPoint Set Type   : ", cg_PointSetTypeName(r->data[2]),
                NULL);
            if (r->data[2] == CGNS_ENUMV(PointList) || r->data[2] == CGNS_ENUMV(ElementList)) {
                sprintf (buff, "%ld", (long)r->data[3]);
                Tcl_AppendResult (interp, "\nNumber of Points : ", buff, NULL);
            }
            else if (r->data[3] == 2) {
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[4], (long)r->data[5]);
                Tcl_AppendResult (interp, "\nIndex Range      : ", buff, NULL);
            }
            else {
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[4], (long)r->data[4+CellDim]);
                Tcl_AppendResult (interp,
                    "\nI Index Range    : ", buff, NULL);
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[5], (long)r->data[5+CellDim]);
                Tcl_AppendResult (interp,
                    "\nJ Index Range    : ", buff, NULL);
                if (CellDim == 3) {
                    sprintf (buff, "%ld -> %ld",
                        (long)r->data[6], (long)r->data[9]);
                    Tcl_AppendResult (interp,
                        "\nK Index Range    : ", buff, NULL);
                }
            }
            break;
        case REG_BOCO:
            Tcl_AppendResult (interp,
                  "Region Name     : ", r->name,
                "\nType of Region  : Boundary Condition",
                "\nType of BC      : ", cg_BCTypeName(r->data[0]),
                "\nGrid Location   : ", cg_GridLocationName(r->data[1]),
                "\nPoint Set Type  : ", cg_PointSetTypeName(r->data[2]),
                NULL);
            if (r->data[2] == CGNS_ENUMV(PointList) ||
                r->data[2] == CGNS_ENUMV(ElementList)) {
                sprintf (buff, "%ld", (long)r->data[3]);
                Tcl_AppendResult (interp, "\nNumber of Points: ", buff, NULL);
            }
            else if (r->data[3] == 2) {
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[4], (long)r->data[5]);
                Tcl_AppendResult (interp, "\nIndex Range     : ", buff, NULL);
            }
            else {
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[4], (long)r->data[4+CellDim]);
                Tcl_AppendResult (interp,
                    "\nI Index Range   : ", buff, NULL);
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[5], (long)r->data[5+CellDim]);
                Tcl_AppendResult (interp,
                    "\nJ Index Range   : ", buff, NULL);
                if (CellDim == 3) {
                    sprintf (buff, "%ld -> %ld",
                        (long)r->data[6], (long)r->data[9]);
                    Tcl_AppendResult (interp,
                        "\nK Index Range   : ", buff, NULL);
                }
            }
            break;
        case REG_BNDS:
            strcpy (buff, r->name);
            Tcl_AppendResult (interp,
                  "Region Name   : ", r->name,
                "\nType of Region: Mesh Boundary", NULL);
            sprintf (buff, "%ld -> %ld",
                (long)r->data[0], (long)r->data[CellDim]);
            Tcl_AppendResult (interp,
                "\nI Index Range : ", buff, NULL);
            sprintf (buff, "%ld -> %ld",
                (long)r->data[1], (long)r->data[1+CellDim]);
            Tcl_AppendResult (interp,
                "\nJ Index Range : ", buff, NULL);
            if (CellDim == 3) {
                sprintf (buff, "%ld -> %ld",
                    (long)r->data[2], (long)r->data[5]);
                Tcl_AppendResult (interp,
                    "\nK Index Range : ", buff, NULL);
            }
            break;
    }
    return TCL_OK;
}

/*---------- CGNSregiondim -----------------------------------------
 * return dimension of a region
 *------------------------------------------------------------------*/

static int CGNSregiondim (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int n;
    char buff[16];
    Zone *z;

    if (!cgnsfn) {
        Tcl_SetResult (interp, "CGNS file not open", TCL_STATIC);
        return TCL_ERROR;
    }
    if (argc != 3) {
        Tcl_SetResult (interp, "usage: CGNSregtype zone reg", TCL_STATIC);
        return TCL_ERROR;
    }
    n = atoi (argv[1]);
    if (n < 0 || n >= nzones) {
        Tcl_SetResult (interp, "zone number out of range", TCL_STATIC);
        return TCL_ERROR;
    }
    z = &zones[n];
    n = atoi (argv[2]);
    if (n < 0 || n >= z->nregs) {
        Tcl_SetResult (interp, "region number out of range", TCL_STATIC);
        return TCL_ERROR;
    }
    Tcl_ResetResult (interp);
    if (!z->regs[n].dim) {
        if (*(z->regs[n].errmsg))
            Tcl_SetResult (interp, z->regs[n].errmsg, TCL_STATIC);
        return TCL_ERROR;
    }
    sprintf (buff, "%d", z->regs[n].dim);
    Tcl_AppendResult (interp, buff, NULL);
    return TCL_OK;
}

/*---------- CGNSbounds --------------------------------------------
 * get bounding box
 *------------------------------------------------------------------*/

static void transform_bounds (float m[16], float bb[3][2])
{
    int i, j;
    float x, y, z, bbox[3][2];

    x = m[0] * bb[0][0] + m[4] * bb[1][0] +  m[8] * bb[2][0] + m[12];
    y = m[1] * bb[0][0] + m[5] * bb[1][0] +  m[9] * bb[2][0] + m[13];
    z = m[2] * bb[0][0] + m[6] * bb[1][0] + m[10] * bb[2][0] + m[14];
    bbox[0][0] = bbox[0][1] = x;
    bbox[1][0] = bbox[1][1] = y;
    bbox[2][0] = bbox[2][1] = z;

    x = m[0] * bb[0][1] + m[4] * bb[1][0] +  m[8] * bb[2][0] + m[12];
    y = m[1] * bb[0][1] + m[5] * bb[1][0] +  m[9] * bb[2][0] + m[13];
    z = m[2] * bb[0][1] + m[6] * bb[1][0] + m[10] * bb[2][0] + m[14];
    if (bbox[0][0] > x) bbox[0][0] = x;
    if (bbox[0][1] < x) bbox[0][1] = x;
    if (bbox[1][0] > y) bbox[1][0] = y;
    if (bbox[1][1] < y) bbox[1][1] = y;
    if (bbox[2][0] > z) bbox[2][0] = z;
    if (bbox[2][1] < z) bbox[2][1] = z;

    x = m[0] * bb[0][0] + m[4] * bb[1][1] +  m[8] * bb[2][0] + m[12];
    y = m[1] * bb[0][0] + m[5] * bb[1][1] +  m[9] * bb[2][0] + m[13];
    z = m[2] * bb[0][0] + m[6] * bb[1][1] + m[10] * bb[2][0] + m[14];
    if (bbox[0][0] > x) bbox[0][0] = x;
    if (bbox[0][1] < x) bbox[0][1] = x;
    if (bbox[1][0] > y) bbox[1][0] = y;
    if (bbox[1][1] < y) bbox[1][1] = y;
    if (bbox[2][0] > z) bbox[2][0] = z;
    if (bbox[2][1] < z) bbox[2][1] = z;

    x = m[0] * bb[0][1] + m[4] * bb[1][1] +  m[8] * bb[2][0] + m[12];
    y = m[1] * bb[0][1] + m[5] * bb[1][1] +  m[9] * bb[2][0] + m[13];
    z = m[2] * bb[0][1] + m[6] * bb[1][1] + m[10] * bb[2][0] + m[14];
    if (bbox[0][0] > x) bbox[0][0] = x;
    if (bbox[0][1] < x) bbox[0][1] = x;
    if (bbox[1][0] > y) bbox[1][0] = y;
    if (bbox[1][1] < y) bbox[1][1] = y;
    if (bbox[2][0] > z) bbox[2][0] = z;
    if (bbox[2][1] < z) bbox[2][1] = z;

    x = m[0] * bb[0][0] + m[4] * bb[1][0] +  m[8] * bb[2][1] + m[12];
    y = m[1] * bb[0][0] + m[5] * bb[1][0] +  m[9] * bb[2][1] + m[13];
    z = m[2] * bb[0][0] + m[6] * bb[1][0] + m[10] * bb[2][1] + m[14];
    if (bbox[0][0] > x) bbox[0][0] = x;
    if (bbox[0][1] < x) bbox[0][1] = x;
    if (bbox[1][0] > y) bbox[1][0] = y;
    if (bbox[1][1] < y) bbox[1][1] = y;
    if (bbox[2][0] > z) bbox[2][0] = z;
    if (bbox[2][1] < z) bbox[2][1] = z;

    x = m[0] * bb[0][1] + m[4] * bb[1][0] +  m[8] * bb[2][1] + m[12];
    y = m[1] * bb[0][1] + m[5] * bb[1][0] +  m[9] * bb[2][1] + m[13];
    z = m[2] * bb[0][1] + m[6] * bb[1][0] + m[10] * bb[2][1] + m[14];
    if (bbox[0][0] > x) bbox[0][0] = x;
    if (bbox[0][1] < x) bbox[0][1] = x;
    if (bbox[1][0] > y) bbox[1][0] = y;
    if (bbox[1][1] < y) bbox[1][1] = y;
    if (bbox[2][0] > z) bbox[2][0] = z;
    if (bbox[2][1] < z) bbox[2][1] = z;

    x = m[0] * bb[0][0] + m[4] * bb[1][1] +  m[8] * bb[2][1] + m[12];
    y = m[1] * bb[0][0] + m[5] * bb[1][1] +  m[9] * bb[2][1] + m[13];
    z = m[2] * bb[0][0] + m[6] * bb[1][1] + m[10] * bb[2][1] + m[14];
    if (bbox[0][0] > x) bbox[0][0] = x;
    if (bbox[0][1] < x) bbox[0][1] = x;
    if (bbox[1][0] > y) bbox[1][0] = y;
    if (bbox[1][1] < y) bbox[1][1] = y;
    if (bbox[2][0] > z) bbox[2][0] = z;
    if (bbox[2][1] < z) bbox[2][1] = z;

    x = m[0] * bb[0][1] + m[4] * bb[1][1] +  m[8] * bb[2][1] + m[12];
    y = m[1] * bb[0][1] + m[5] * bb[1][1] +  m[9] * bb[2][1] + m[13];
    z = m[2] * bb[0][1] + m[6] * bb[1][1] + m[10] * bb[2][1] + m[14];
    if (bbox[0][0] > x) bbox[0][0] = x;
    if (bbox[0][1] < x) bbox[0][1] = x;
    if (bbox[1][0] > y) bbox[1][0] = y;
    if (bbox[1][1] < y) bbox[1][1] = y;
    if (bbox[2][0] > z) bbox[2][0] = z;
    if (bbox[2][1] < z) bbox[2][1] = z;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 2; j++)
            bb[i][j] = bbox[i][j];
}

/*-------------------------------------------------------------------*/

static int CGNSbounds (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    float bbox[3][2], matrix[16];
    int n, all = 0;
    CONST char **args;
    char sbb[65];

    if (argc > 1) all = atoi(argv[1]);
    get_bounds (all, bbox);
    if (argc > 2) {
        if (TCL_OK != Tcl_SplitList (interp, argv[2], &n, &args))
            return TCL_ERROR;
        for (n = 0; n < 16; n++)
            matrix[n] = (float) atof (args[n]);
        Tcl_Free ((char *)args);
        transform_bounds (matrix, bbox);
    }
    Tcl_ResetResult (interp);
    for (n = 0; n < 3; n++) {
        sprintf (sbb, "%f %f", bbox[n][0], bbox[n][1]);
        Tcl_AppendElement (interp, sbb);
    }
    return TCL_OK;
}

/*---------- OGLregion ---------------------------------------------
 * create OGL display list for region
 *------------------------------------------------------------------*/

static int OGLregion (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int zone, regn, nc;
    CONST char **args;
    Zone *z;
    Regn *r;
    static char slist[17];

    if (argc != 5) {
        Tcl_SetResult (interp, "usage: OGLregion zone region mode color",
            TCL_STATIC);
        return TCL_ERROR;
    }
    zone = atoi (argv[1]);
    if (zone < 0 || zone >= nzones) {
        Tcl_SetResult (interp, "zone number out of range", TCL_STATIC);
        return TCL_ERROR;
    }
    z = &zones[zone];
    regn = atoi (argv[2]);
    if (regn < 0 || regn >= z->nregs) {
        Tcl_SetResult (interp, "region number out of range", TCL_STATIC);
        return TCL_ERROR;
    }
    r = &z->regs[regn];

    if (r->nfaces || r->nedges) {
        r->mode = atoi (argv[3]);
        if (TCL_OK != Tcl_SplitList (interp, argv[4], &nc, &args))
            return TCL_ERROR;
        if (nc != 3) {
            Tcl_Free ((char *)args);
            Tcl_SetResult (interp, "invalid color", TCL_STATIC);
            return TCL_ERROR;
        }
        for (nc = 0; nc < 3; nc++)
            r->color[nc] = (float)atof (args[nc]);
        r->color[3] = 1.0;
        Tcl_Free ((char *)args);

        if (!r->dlist) r->dlist = glGenLists (1);
        glNewList (r->dlist, GL_COMPILE);
        glColor3fv (r->color);
        glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, r->color);
        if (r->nfaces) {
            switch (r->mode) {
                case 1:
                    draw_outlines (z, r);
                    break;
                case 2:
                    draw_mesh (z, r);
                    break;
                case 3:
                    draw_shaded (z, r);
                    break;
                default:
                    r->mode = 0;
                    break;
            }
        }
        else if (r->mode < 1 || r->mode > 3) {
            r->mode = 0;
        }
        else {
            draw_outlines (z, r);
        }
        glEndList ();
    }

    sprintf (slist, "%d", r->dlist);
    Tcl_SetResult (interp, slist, TCL_STATIC);
    return TCL_OK;
}

/*---------- OGLaxis -----------------------------------------------
 * create OGL display list for axis
 *------------------------------------------------------------------*/

#define CHAR_W 8
#define CHAR_H 13

static GLubyte x_raster[] = {
    0x00, 0x00, 0xc3, 0x66, 0x66, 0x3c, 0x3c, 0x18,
    0x3c, 0x3c, 0x66, 0x66, 0xc3
};
static GLubyte y_raster[] = {
    0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18,
    0x3c, 0x3c, 0x66, 0x66, 0xc3
};
static GLubyte z_raster[] = {
    0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x7e,
    0x0c, 0x06, 0x03, 0x03, 0xff
};

/*-------------------------------------------------------------------*/

static int OGLaxis (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int vis;
    float bbox[3][2];
    static char slist[17];

    if (argc < 2 || argc > 3) {
        Tcl_SetResult (interp, "usage: OGLaxis visible [bounds]",
            TCL_STATIC);
        return TCL_ERROR;
    }
    vis = atoi (argv[1]);
    if (!AxisDL) AxisDL = glGenLists (1);

    glNewList (AxisDL, GL_COMPILE);
    if (vis) {
        if (argc == 3) {
            int nb, n = 0;
            CONST char **args;
            if (TCL_OK != Tcl_SplitList (interp, argv[2], &nb, &args))
                return TCL_ERROR;
            if (nb == 3) {
                for (n = 0; n < nb; n++) {
                    if (sscanf (args[n], "%f %f", &bbox[n][0], &bbox[n][1]) != 2)
                        break;
                }
            }
            Tcl_Free ((char *)args);
            if (n != 3) {
                Tcl_SetResult (interp, "invalid bounding box", TCL_STATIC);
                return TCL_ERROR;
            }
        }
        else
            get_bounds (0, bbox);
        glLineWidth (3.0);
        glDisable (GL_LIGHTING);
        glShadeModel (GL_FLAT);
        glBegin (GL_LINES);
        glColor3f (1.0, 0.0, 0.0);
        glVertex3f (bbox[0][0], bbox[1][0], bbox[2][0]);
        glVertex3f (bbox[0][1], bbox[1][0], bbox[2][0]);
        glColor3f (1.0, 1.0, 0.0);
        glVertex3f (bbox[0][0], bbox[1][0], bbox[2][0]);
        glVertex3f (bbox[0][0], bbox[1][1], bbox[2][0]);
        glColor3f (0.0, 1.0, 0.0);
        glVertex3f (bbox[0][0], bbox[1][0], bbox[2][0]);
        glVertex3f (bbox[0][0], bbox[1][0], bbox[2][1]);
        glEnd ();
        glLineWidth (1.0);

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glColor3f (1.0, 0.0, 0.0);
        glRasterPos3f (bbox[0][1], bbox[1][0], bbox[2][0]);
        glBitmap (CHAR_W, CHAR_H, -0.5 * (float)CHAR_W,
            0.5 * (float)CHAR_H, 0.0, 0.0, x_raster);
        glColor3f (1.0, 1.0, 0.0);
        glRasterPos3f (bbox[0][0], bbox[1][1], bbox[2][0]);
        glBitmap (CHAR_W, CHAR_H, -0.5 * (float)CHAR_W,
            0.5 * (float)CHAR_H, 0.0, 0.0, y_raster);
        glColor3f (0.0, 1.0, 0.0);
        glRasterPos3f (bbox[0][0], bbox[1][0], bbox[2][1]);
        glBitmap (CHAR_W, CHAR_H, -0.5 * (float)CHAR_W,
            0.5 * (float)CHAR_H, 0.0, 0.0, z_raster);
    }
    glEndList ();

    sprintf (slist, "%d", AxisDL);
    Tcl_SetResult (interp, slist, TCL_STATIC);
    return TCL_OK;
}

/*-------------------------------------------------------------------*/

static int OGLcolor (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int index, i, j;
    double r, g, b;
    double h, v, s;
    static char color[256];
    static double huemap[12] = {0,1,2,3,4,5,0.5,1.25,2.65,3.4,4.5,5.5};

    if (argc != 2) {
        Tcl_SetResult (interp, "usage: OGLcolor index",
            TCL_STATIC);
        return TCL_ERROR;
    }
    index = abs(atoi(argv[1])) % 132;
    h = huemap[index % 12];
    i = (int)h;
    h -= (double)i;
    j = index / 12;
    if ((j % 2) == 0) {
        v = 1.0;
        s = 1.0 - sqrt((double)j / 22.0);
    }
    else {
        v = 1.0 - sqrt((double)j / 44.0);
        s = 1.0;
    }
    r = g = b = 0.0;
    switch (i) {
        case 6:
            h = 0.0;
        case 0:
            r = v;
            g = v * (1.0 - (s * (1.0 - h)));
            b = v * (1.0 - s);
            break;
        case 1:
            r = v * (1.0 - (s * h));
            g = v;
            b = v * (1.0 - s);
            break;
        case 2:
            r = v * (1.0 - s);
            g = v;
            b = v * (1.0 - (s * (1.0 - h)));
            break;
        case 3:
            r = v * (1.0 - s);
            g = v * (1.0 - (s * h));
            b = v;
            break;
        case 4:
            r = v * (1.0 - (s * (1.0 - h)));
            g = v * (1.0 - s);
            b = v;
            break;
        case 5:
            r = v;
            g = v * (1.0 - s);
            b = v * (1.0 - (s * h));
            break;
    }

    if (r < 0.0) r = 0.0;
    if (r > 1.0) r = 1.0;
    if (g < 0.0) g = 0.0;
    if (g > 1.0) g = 1.0;
    if (b < 0.0) b = 0.0;
    if (b > 1.0) b = 1.0;

    sprintf(color, "%g %g %g", r, g, b);
    Tcl_SetResult (interp, color, TCL_STATIC);
    return TCL_OK;
}

/*==============================================================
 * cutting plane routines
 *==============================================================*/

#ifndef NO_CUTTING_PLANE

/*------------------------------------------------------------------*/

static void init_cutplane (float plane[4])
{
    int nz, nr;

    for (nz = 0; nz < nzones; nz++) {
        for (nr = 0; nr < zones[nz].nregs; nr++) {
            if (zones[nz].regs[nr].cut.nelems) {
                free (zones[nz].regs[nr].cut.elems);
                zones[nz].regs[nr].cut.nelems = 0;
            }
            if (zones[nz].regs[nr].cut.nedges) {
                free (zones[nz].regs[nr].cut.edges);
                zones[nz].regs[nr].cut.nedges = 0;
            }
        }
    }
    if (cutplane.nnodes) free (cutplane.nodes);
    cutplane.nelems = 0;
    cutplane.nedges = 0;
    cutplane.nnodes = 0;
    for (nr = 0; nr < 4; nr++)
        cutplane.plane[nr] = plane[nr];
}

/*------------------------------------------------------------------*/

static int classify_element (Zone *z, int nnodes, cgsize_t *nodeid)
{
    int n, index = 0;
    int mask = (1 << (nnodes - 1)) - 1;
    float *node;
    double s;

    for (n = 0; n < nnodes; n++) {
        node = z->nodes[nodeid[n]];
        s = node[0] * cutplane.plane[0] + node[1] * cutplane.plane[1] +
            node[2] * cutplane.plane[2] - cutplane.plane[3];
        if (s >= 0.0) index |= (1 << n);
    }
    if (index > mask) index ^= ((mask << 1) | 1);
    return index;
}

/*------------------------------------------------------------------*/

static int classify_polygon (Zone *z, int nnodes, cgsize_t *nodeid)
{
    int n, start, diff;
    float *node;
    double s;

    node = z->nodes[nodeid[0]];
    s = node[0] * cutplane.plane[0] + node[1] * cutplane.plane[1] +
        node[2] * cutplane.plane[2] - cutplane.plane[3];
    start = s >= 0.0 ? 1 : -1;

    for (n = 1; n < nnodes; n++) {
        node = z->nodes[nodeid[n]];
        s = node[0] * cutplane.plane[0] + node[1] * cutplane.plane[1] +
            node[2] * cutplane.plane[2] - cutplane.plane[3];
        diff = start * (s >= 0.0 ? 1 : -1);
        if (diff < 0) return 1;
    }
    return 0;
}

/*------------------------------------------------------------------*/

static cgsize_t find_elements ()
{
#define ELEM_INC 50
    int nz, nnodes, nn, nr, nf;
    cgsize_t n, ne, maxelems, nelems, *elems;
    CGNS_ENUMT(ElementType_t) type;
    Zone *z;
    Regn *r;
    Face *f;

    for (nz = 0; nz < nzones; nz++) {
        z = &zones[nz];
        for (nr = 0; nr < z->nregs; nr++) {
            r = &z->regs[nr];
            if (r->dim < 2 || (r->mode == 0 && !ignorevis)) continue;
            type = r->elemtype;
            cg_npe(type, &nnodes);
            maxelems = nelems = 0;
            elems = NULL;

            if (type == CGNS_ENUMV(NGON_n)) {
                for (n = 0, ne = 0; ne < r->nelems; ne++) {
                    nn = r->elems[n++];
                    if (nn > 2 &&
                        classify_polygon(z, nn, &r->elems[n])) {
                        if (nelems >= maxelems) {
                            maxelems += ELEM_INC;
                            elems = (cgsize_t *) REALLOC ("find_elements",
                                (size_t)maxelems * sizeof(cgsize_t), elems);
                        }
                        elems[nelems++] = n - 1;
                    }
                    n += nn;
                }
            }
            else if (type == CGNS_ENUMV(NFACE_n)) {
                for (n = 0, ne = 0; ne < r->nelems; ne++) {
                    nf = r->elems[n++];
                    for (nn = 0; nn < nf; nn++) {
                        f = r->poly[abs(r->elems[n+nn])-1];
                        if (f->nnodes > 2 &&
                            classify_polygon(z, f->nnodes, f->nodes)) {
                            if (nelems >= maxelems) {
                                maxelems += ELEM_INC;
                                elems = (cgsize_t *) REALLOC ("find_elements",
                                    (size_t)maxelems * sizeof(cgsize_t), elems);
                            }
                            elems[nelems++] = n - 1;
                            break;
                        }
                    }
                    n += nf;
                }
            }
            else {
                for (n = 0, ne = 0; ne < r->nelems; ne++) {
                    if (r->elemtype == CGNS_ENUMV(MIXED)) {
                        type = (CGNS_ENUMT(ElementType_t))r->elems[n++];
                        cg_npe(type, &nnodes);
                    }
                    switch (type) {
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
                        case CGNS_ENUMV(PYRA_13):
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
                    if (nn && classify_element(z, nn, &r->elems[n])) {
                        if (nelems >= maxelems) {
                            maxelems += ELEM_INC;
                            elems = (cgsize_t *) REALLOC ("find_elements",
                                (size_t)maxelems * sizeof(cgsize_t), elems);
                        }
                        if (r->elemtype == CGNS_ENUMV(MIXED))
                            elems[nelems] = n - 1;
                        else
                            elems[nelems] = n;
                        nelems++;
                    }
                    n += nnodes;
                }
            }
            r->cut.nelems = nelems;
            r->cut.elems = elems;
            cutplane.nelems += nelems;
        }
    }

    return cutplane.nelems;
}

/*-------------------------------------------------------------------*/

static int compare_cut_node (void *v1, void *v2)
{
    int i;
    CutNode *c1 = (CutNode *)v1;
    CutNode *c2 = (CutNode *)v2;

    for (i = 0; i < 3; i++) {
        if (c1->nodes[i] != c2->nodes[i])
            return (int)(c1->nodes[i] - c2->nodes[i]);
    }
    return 0;
}

/*-------------------------------------------------------------------*/

static size_t hash_cut_node (void *v)
{
    CutNode *c = (CutNode *)v;

    return ((size_t)(c->nodes[0] + c->nodes[1] + c->nodes[2]));
}

/*-------------------------------------------------------------------*/

static void get_cut_node (void *v)
{
    int i;
    CutNode *c = (CutNode *)v;
    Zone *z = &zones[c->nodes[0]];
    float *n, *n1, *n2;

    n  = cutplane.nodes[c->id];
    n1 = z->nodes[c->nodes[1]];
    n2 = z->nodes[c->nodes[2]];

    for (i = 0; i < 3; i++)
        n[i] = n1[i] + c->ratio * (n2[i] - n1[i]);
    free (c);
}

/*-------------------------------------------------------------------*/

static size_t get_cut_edge (void *ve, void *vc)
{
    Edge *e = (Edge *)ve;
    CutData *c = (CutData *)vc;

    c->edges[c->nedges].nodes[0] = e->nodes[0];
    c->edges[c->nedges].nodes[1] = e->nodes[1];
    (c->nedges)++;
    return 1;
}

/*----- tri elements -----*/

#define TRI_SIZE  3
#define TRI_EDGES 3

static int triCuts[TRI_SIZE+1][4] = {
    {0},
    {2,0,2, 0},
    {2,0,1, 0},
    {2,1,2, 0}
};
static int triEdges[TRI_EDGES][2] = {
    {0,1},
    {1,2},
    {2,0}
};

/*----- quad elements -----*/

#define QUAD_SIZE  7
#define QUAD_EDGES 4

static int quadCuts[QUAD_SIZE+1][4] = {
    {0},
    {2,0,3, 0},
    {2,0,1, 0},
    {2,1,3, 0},
    {2,1,2, 0},
    {2,0,3, 0},
    {2,0,2, 0},
    {2,2,3, 0}
};
static int quadEdges[QUAD_EDGES][2] = {
    {0,1},
    {1,2},
    {2,3},
    {3,0}
};

/*----- tet elements -----*/

#define TET_SIZE  7
#define TET_EDGES 6

static int tetCuts[TET_SIZE+1][6] = {
    {0},
    {3,0,3,2, 0},
    {3,0,1,4, 0},
    {4,1,4,3,2, 0},
    {3,1,2,5, 0},
    {4,0,3,5,1, 0},
    {4,0,2,5,4, 0},
    {3,3,5,4, 0}
};
static int tetEdges[TET_EDGES][2] = {
    {0,1},
    {1,2},
    {2,0},
    {0,3},
    {1,3},
    {2,3}
};

/*----- pyramid elements -----*/

#define PYR_SIZE  15
#define PYR_EDGES 8

static int pyrCuts[PYR_SIZE+1][9] = {
    {0},
    {3,0,4,3, 0},
    {3,0,1,5, 0},
    {4,1,5,4,3, 0},
    {3,1,2,6, 0},
    {3,0,4,3, 3,1,2,6, 0},
    {4,0,2,6,5, 0},
    {5,3,2,6,5,4, 0},
    {3,2,3,7, 0},
    {4,0,4,7,2, 0},
    {3,2,3,7, 3,0,1,5, 0},
    {5,1,5,4,7,2, 0},
    {4,1,3,7,6, 0},
    {5,0,4,7,6,1, 0},
    {5,0,3,7,6,5, 0},
    {4,4,7,6,5, 0},
};
static int pyrEdges[PYR_EDGES][2] = {
    {0,1},
    {1,2},
    {2,3},
    {3,0},
    {0,4},
    {1,4},
    {2,4},
    {3,4}
};

/*----- wedge elements -----*/

#define WDG_SIZE  31
#define WDG_EDGES 9

static int wdgCuts[WDG_SIZE+1][10] = {
    {0},
    {3,0,3,2, 0},
    {3,0,1,4, 0},
    {4,1,4,3,2, 0},
    {3,1,2,5, 0},
    {4,0,3,5,1, 0},
    {4,0,2,5,4, 0},
    {3,3,5,4, 0},
    {3,3,6,8, 0},
    {4,0,6,8,2, 0},
    {3,3,6,8, 3,0,1,4, 0},
    {5,1,4,6,8,2, 0},
    {3,3,6,8, 3,1,2,5, 0},
    {5,0,6,8,5,1, 0},
    {3,3,6,8, 4,0,2,5,4, 0},
    {4,4,6,8,5, 0},
    {3,4,7,6, 0},
    {3,4,7,6, 3,0,3,2, 0},
    {4,0,1,7,6, 0},
    {5,1,7,6,3,2, 0},
    {3,4,7,6, 3,1,2,5, 0},
    {3,4,7,6, 4,0,3,5,1, 0},
    {5,0,2,5,7,6, 0},
    {4,3,5,7,6, 0},
    {4,3,4,7,8, 0},
    {5,0,4,7,8,2, 0},
    {5,0,1,7,8,3, 0},
    {4,1,7,8,2, 0},
    {4,3,4,7,8, 3,1,2,5, 0},
    {3,0,4,1, 3,5,7,8, 0},
    {3,0,2,3, 3,5,7,8, 0},
    {3,5,7,8, 0}
};
static int wdgEdges[WDG_EDGES][2] = {
    {0,1},
    {1,2},
    {2,0},
    {0,3},
    {1,4},
    {2,5},
    {3,4},
    {4,5},
    {5,3}
};

/*----- hex elements -----*/

#define HEX_SIZE  127
#define HEX_EDGES 12

static int hexCuts[HEX_SIZE+1][17] = {
    {0},
    {3,3,0,4, 0},
    {3,0,1,5, 0},
    {4,3,1,5,4, 0},
    {3,1,2,6, 0},
    {3,3,0,4, 3,1,2,6, 0},
    {4,0,2,6,5, 0},
    {5,3,2,6,5,4, 0},
    {3,3,7,2, 0},
    {4,0,4,7,2, 0},
    {3,0,1,5, 3,3,7,2, 0},
    {5,4,7,2,1,5, 0},
    {4,3,7,6,1, 0},
    {5,0,4,7,6,1, 0},
    {5,0,3,7,6,5, 0},
    {4,7,6,5,4, 0},
    {3,4,11,8, 0},
    {4,0,8,11,3, 0},
    {3,0,1,5, 3,4,11,8, 0},
    {5,3,11,8,5,1, 0},
    {3,1,2,6, 3,4,11,8, 0},
    {4,0,8,11,3, 3,1,2,6, 0},
    {4,0,2,6,5, 3,4,11,8, 0},
    {6,3,2,6,5,8,11, 0},
    {3,3,7,2, 3,4,11,8, 0},
    {5,0,2,7,11,8, 0},
    {3,0,1,5, 3,3,7,2, 3,4,11,8, 0},
    {6,2,1,5,8,11,7, 0},
    {4,3,7,6,1, 3,4,11,8, 0},
    {6,0,8,11,7,6,1, 0},
    {5,0,3,7,6,5, 3,4,11,8, 0},
    {5,11,7,6,5,8, 0},
    {3,8,9,5, 0},
    {3,8,9,5, 3,3,0,4, 0},
    {4,0,8,9,1, 0},
    {5,4,3,1,9,8, 0},
    {3,8,9,5, 3,1,2,6, 0},
    {3,1,9,5, 3,3,0,4, 3,1,2,6, 0},
    {5,0,2,6,9,8, 0},
    {6,3,2,6,9,8,4, 0},
    {3,8,9,5, 3,3,7,2, 0},
    {4,0,4,7,2, 3,8,9,5, 0},
    {4,0,8,9,1, 3,3,7,2, 0},
    {6,4,7,2,1,9,8, 0},
    {4,3,7,6,1, 3,8,9,5, 0},
    {5,7,6,1,0,4, 3,8,9,5, 0},
    {6,0,3,7,6,9,8, 0},
    {5,4,7,6,9,8, 0},
    {4,4,11,9,5, 0},
    {5,0,3,11,9,5, 0},
    {5,0,4,11,9,1, 0},
    {4,3,1,9,11, 0},
    {4,4,11,9,5, 3,1,2,6, 0},
    {6,0,3,11,9,6,2, 0},
    {6,0,2,6,9,11,4, 0},
    {5,3,2,6,9,11, 0},
    {4,4,11,9,5, 3,3,7,2, 0},
    {6,0,2,7,11,9,5, 0},
    {5,11,9,1,0,4, 3,3,7,2, 0},
    {5,11,7,2,1,9, 0},
    {4,3,7,6,1, 4,4,11,9,5, 0},
    {4,11,7,6,9, 3,0,1,5, 0},
    {4,11,7,6,9, 3,3,0,4, 0},
    {4,11,7,6,9, 0},
    {3,9,10,6, 0},
    {3,9,10,6, 3,3,0,4, 0},
    {3,9,10,6, 3,0,1,5, 0},
    {4,4,3,1,5, 3,9,10,6, 0},
    {4,2,1,9,10, 0},
    {4,2,1,9,10, 3,3,0,4, 0},
    {5,0,2,10,9,5, 0},
    {6,3,2,10,9,5,4, 0},
    {3,3,7,2, 3,9,10,6, 0},
    {4,7,2,0,4, 3,9,10,6, 0},
    {3,0,1,5, 3,2,3,7, 3,9,10,6, 0},
    {4,4,7,10,9, 4,1,2,6,5, 0},
    {5,3,7,10,9,1, 0},
    {6,4,7,10,9,1,0, 0},
    {6,3,7,10,9,5,0, 0},
    {5,4,7,10,9,5, 0},
    {3,4,11,8, 3,9,10,6, 0},
    {4,0,8,11,3, 3,9,10,6, 0},
    {3,0,1,5, 3,4,11,8, 3,9,10,6, 0},
    {5,3,11,10,6,1, 3,1,9,5, 0},
    {4,1,2,10,9, 3,4,11,8, 0},
    {4,3,11,10,2, 4,0,8,9,1, 0},
    {5,4,11,10,2,0, 3,8,9,5, 0},
    {4,3,11,10,2, 3,8,9,5, 0},
    {3,3,7,2, 3,9,10,6, 3,4,11,8, 0},
    {5,0,2,6,9,8, 3,11,7,10, 0},
    {3,0,4,3, 3,11,7,10, 3,1,2,6, 3,8,9,5, 0},
    {3,1,9,5, 3,11,7,10, 3,1,2,6, 0},
    {5,3,1,9,8,4, 3,11,7,10, 0},
    {4,0,8,9,1, 3,11,7,10, 0},
    {3,11,7,10, 3,8,9,5, 3,3,0,4, 0},
    {3,11,7,10, 3,8,9,5, 0},
    {4,8,10,6,5, 0},
    {4,8,10,6,5, 3,3,0,4, 0},
    {5,0,8,10,6,1, 0},
    {6,3,4,8,10,6,1, 0},
    {5,10,2,1,5,8, 0},
    {5,3,2,10,8,4, 3,0,1,5, 0},
    {4,0,8,10,2, 0},
    {5,3,2,10,8,4, 0},
    {4,8,10,6,5, 3,3,7,2, 0},
    {4,4,7,10,8, 4,0,2,6,5, 0},
    {5,3,7,10,8,0, 3,1,2,6, 0},
    {4,4,7,10,8, 3,1,2,6, 0},
    {6,3,7,10,8,5,1, 0},
    {4,4,7,10,8, 3,0,1,5, 0},
    {5,3,7,10,8,0, 0},
    {4,4,7,10,8, 0},
    {5,4,11,10,6,5, 0},
    {6,0,3,11,10,6,5, 0},
    {6,4,11,10,6,1,0, 0},
    {5,3,11,10,6,1, 0},
    {6,4,11,10,2,1,5, 0},
    {4,3,11,10,2, 3,0,1,5, 0},
    {5,4,11,10,2,0, 0},
    {4,11,10,2,3, 0},
    {5,3,2,6,5,4, 3,11,7,10, 0},
    {4,0,2,6,5, 3,11,7,10, 0},
    {3,3,0,4, 3,1,2,6, 3,11,7,10, 0},
    {3,1,2,6, 3,11,7,10, 0},
    {4,4,3,1,5, 3,11,7,10, 0},
    {3,0,1,5, 3,11,7,10, 0},
    {3,11,7,10, 3,3,0,4, 0},
    {3,11,7,10, 0}
};
static int hexEdges[HEX_EDGES][2] = {
    {0,1},
    {1,2},
    {2,3},
    {3,0},
    {0,4},
    {1,5},
    {2,6},
    {3,7},
    {4,5},
    {5,6},
    {6,7},
    {7,4}
};

static int n_cut_nodes;
static HASH *cut_hash;

/*------------------------------------------------------------------*/

static void intersect_polygon (int zonenum, int nnodes, cgsize_t *nodeid,
                               HASH *edgehash)
{
    int i, n, nn, i1, i2;
    cgsize_t id = -1;
    float *node;
    double s1, s2;
    CutNode cnode, *cn;
    Edge edge, *ep;
    Zone *z = &zones[zonenum];
    static char *funcname = "intersect_polygon";

    if (nnodes < 3) return;
    node = z->nodes[nodeid[0]];
    for (s1 = 0.0, i = 0; i < 3; i++)
        s1 += (node[i] * cutplane.plane[i]);
    i1 = (s1 - cutplane.plane[3]) >= 0.0 ? 1 : -1;

    for (n = 1; n <= nnodes; n++) {
        nn = n % nnodes;
        node = z->nodes[nodeid[nn]];
        for (s2 = 0.0, i = 0; i < 3; i++)
            s2 += (node[i] * cutplane.plane[i]);
        i2 = (s2 - cutplane.plane[3]) >= 0.0 ? 1 : -1;
        if (i1 * i2 < 0) {
            cnode.nodes[0] = zonenum;
            if (nodeid[n-1] < nodeid[nn]) {
                cnode.nodes[1] = nodeid[n-1];
                cnode.nodes[2] = nodeid[nn];
                if (s1 == s2)
                    cnode.ratio = 0.0;
                else
                    cnode.ratio = (float)((cutplane.plane[3] - s1) / (s2 - s1));
            }
            else {
                cnode.nodes[1] = nodeid[nn];
                cnode.nodes[2] = nodeid[n-1];
                if (s1 == s2)
                    cnode.ratio = 0.0;
                else
                    cnode.ratio = (float)((cutplane.plane[3] - s2) / (s1 - s2));
            }
            cn = (CutNode *) HashFind (cut_hash, &cnode);
            if (cn == NULL) {
                cn = (CutNode *) MALLOC (funcname, sizeof(CutNode));
                for (i = 0; i < 3; i++)
                    cn->nodes[i] = cnode.nodes[i];
                cn->id = n_cut_nodes++;
                cn->ratio = cnode.ratio;
                (void) HashAdd (cut_hash, cn);
            }
            if (id >= 0) {
                edge.nodes[0] = id;
                edge.nodes[1] = cn->id;
                ep = (Edge *) HashFind (edgehash, &edge);
                if (NULL == ep) {
                    ep = (Edge *) MALLOC (funcname, sizeof(Edge));
                    ep->nodes[0] = edge.nodes[0];
                    ep->nodes[1] = edge.nodes[1];
                    (void) HashAdd (edgehash, ep);
                }
            }
            id = cn->id;
        }
        s1 = s2;
        i1 = i2;
    }
}

/*------------------------------------------------------------------*/

static void intersect_element (int zonenum, CGNS_ENUMT(ElementType_t) elemtype,
                               cgsize_t *nodeid, HASH *edgehash)
{
    int i, n, index, nn, nc, *cuts;
    int edgemap[HEX_EDGES][2], ids[7];
    cgsize_t swap;
    double s1, s2;
    float *node;
    CutNode cnode, *cn;
    Edge edge, *ep;
    Zone *z = &zones[zonenum];
    static char *funcname = "intersect_element";

    /* get intersection lookup table */

    switch (elemtype) {
        case CGNS_ENUMV(TRI_3):
        case CGNS_ENUMV(TRI_6):
            index = classify_element(z, 3, nodeid);
            if (index < 1 || index > TRI_SIZE) return;
            cuts = triCuts[index];
            for (n = 0; n < TRI_EDGES; n++) {
                edgemap[n][0] = triEdges[n][0];
                edgemap[n][1] = triEdges[n][1];
            }
            break;
        case CGNS_ENUMV(QUAD_4):
        case CGNS_ENUMV(QUAD_8):
        case CGNS_ENUMV(QUAD_9):
            index = classify_element(z, 4, nodeid);
            if (index < 1 || index > QUAD_SIZE) return;
            cuts = quadCuts[index];
            for (n = 0; n < QUAD_EDGES; n++) {
                edgemap[n][0] = quadEdges[n][0];
                edgemap[n][1] = quadEdges[n][1];
            }
            break;
        case CGNS_ENUMV(TETRA_4):
        case CGNS_ENUMV(TETRA_10):
            index = classify_element(z, 4, nodeid);
            if (index < 1 || index > TET_SIZE) return;
            cuts = tetCuts[index];
            for (n = 0; n < TET_EDGES; n++) {
                edgemap[n][0] = tetEdges[n][0];
                edgemap[n][1] = tetEdges[n][1];
            }
            break;
        case CGNS_ENUMV(PYRA_5):
        case CGNS_ENUMV(PYRA_13):
        case CGNS_ENUMV(PYRA_14):
            index = classify_element(z, 5, nodeid);
            if (index < 1 || index > PYR_SIZE) return;
            cuts = pyrCuts[index];
            for (n = 0; n < PYR_EDGES; n++) {
                edgemap[n][0] = pyrEdges[n][0];
                edgemap[n][1] = pyrEdges[n][1];
            }
            break;
        case CGNS_ENUMV(PENTA_6):
        case CGNS_ENUMV(PENTA_15):
        case CGNS_ENUMV(PENTA_18):
            index = classify_element(z, 6, nodeid);
            if (index < 1 || index > WDG_SIZE) return;
            cuts = wdgCuts[index];
            for (n = 0; n < WDG_EDGES; n++) {
                edgemap[n][0] = wdgEdges[n][0];
                edgemap[n][1] = wdgEdges[n][1];
            }
            break;
        case CGNS_ENUMV(HEXA_8):
        case CGNS_ENUMV(HEXA_20):
        case CGNS_ENUMV(HEXA_27):
            index = classify_element(z, 8, nodeid);
            if (index < 1 || index > HEX_SIZE) return;
            cuts = hexCuts[index];
            for (n = 0; n < HEX_EDGES; n++) {
                edgemap[n][0] = hexEdges[n][0];
                edgemap[n][1] = hexEdges[n][1];
            }
            break;
        default:
            return;
    }

    /* get the edge intersections */

    for (nc = 0; cuts[nc];) {
        nn = cuts[nc];
        for (n = 1; n <= nn; n++) {

            /* get edge nodes */

            cnode.nodes[0] = zonenum;
            cnode.nodes[1] = nodeid[edgemap[cuts[nc+n]][0]];
            cnode.nodes[2] = nodeid[edgemap[cuts[nc+n]][1]];
            if (cnode.nodes[2] < cnode.nodes[1]) {
                swap = cnode.nodes[1];
                cnode.nodes[1] = cnode.nodes[2];
                cnode.nodes[2] = swap;
            }
            cn = (CutNode *) HashFind (cut_hash, &cnode);

            /* add node to hash table if not there */

            if (NULL == cn) {
                cn = (CutNode *) MALLOC (funcname, sizeof(CutNode));
                cn->id = n_cut_nodes++;
                for (i = 0; i < 3; i++)
                    cn->nodes[i] = cnode.nodes[i];
                node = z->nodes[cn->nodes[1]];
                for (s1 = 0.0, i = 0; i < 3; i++)
                    s1 += (node[i] * cutplane.plane[i]);
                node = z->nodes[cn->nodes[2]];
                for (s2 = 0.0, i = 0; i < 3; i++)
                    s2 += (node[i] * cutplane.plane[i]);
                if (s1 == s2)
                    cn->ratio = 0.0;
                else
                    cn->ratio = (float)((cutplane.plane[3] - s1) / (s2 - s1));
                (void) HashAdd (cut_hash, cn);
            }
            ids[n-1] = cn->id;
        }
        ids[nn] = ids[0];

        /* add cutplane edge */

        for (n = 0; n < nn; n++) {
            edge.nodes[0] = ids[n];
            edge.nodes[1] = ids[n+1];
            ep = (Edge *) HashFind (edgehash, &edge);
            if (NULL == ep) {
                ep = (Edge *) MALLOC (funcname, sizeof(Edge));
                ep->nodes[0] = edge.nodes[0];
                ep->nodes[1] = edge.nodes[1];
                (void) HashAdd (edgehash, ep);
            }
        }

        /* next cut */

        nc += (nn + 1);
    }
}

/*------------------------------------------------------------------*/

static cgsize_t find_intersects ()
{
    int nz, nr, nf, nfaces, nnodes;
    cgsize_t n, ne;
    size_t nn;
    CGNS_ENUMT(ElementType_t) type;
    Face *f;
    Regn *r;
    HASH *edgehash;
    static char *funcname = "find_intersects";

    /* create hash table to store nodes at edge intersections */

    n_cut_nodes = 0;
    nn = (size_t)cutplane.nelems;
    cut_hash = HashCreate (nn > 1024 ? nn / 3 : 127,
        compare_cut_node, hash_cut_node);
    cutplane.nedges = 0;

    for (nz = 0; nz < nzones; nz++) {
        for (nr = 0; nr < zones[nz].nregs; nr++) {
            r = &zones[nz].regs[nr];
            if (r->cut.nelems == 0) continue;
            type = r->elemtype;

            nn = (size_t)r->cut.nelems;
            edgehash = HashCreate (nn > 1024 ? nn / 3: 127,
                compare_edges, hash_edge);

            for (n = 0, ne = 0; ne < r->cut.nelems; ne++) {
                n = r->cut.elems[ne];
                if (r->elemtype == CGNS_ENUMV(NGON_n)) {
                    nnodes = r->elems[n++];
                    intersect_polygon(nz, nnodes, &r->elems[n], edgehash);
                }
                else if (r->elemtype == CGNS_ENUMV(NFACE_n)) {
                    nfaces = r->elems[n++];
                    for (nf = 0; nf < nfaces; nf++) {
                        f = r->poly[abs(r->elems[n+nf])-1];
                        intersect_polygon(nz, f->nnodes, f->nodes, edgehash);
                    }
                }
                else {
                    cg_npe(type, &nnodes);
                    if (r->elemtype == CGNS_ENUMV(MIXED)) {
                        type = (CGNS_ENUMT(ElementType_t))r->elems[n++];
                        cg_npe(type, &nnodes);
                    }
                    intersect_element (nz, type, &r->elems[n], edgehash);
                }
            }

            r->cut.nedges = 0;
            nn = HashSize (edgehash);
            if (nn) {
                r->cut.edges = (Edge *) MALLOC (funcname, nn * sizeof(Edge));
                HashList (edgehash, get_cut_edge, &r->cut);
            }
            HashDestroy (edgehash, NULL);
            cutplane.nedges += r->cut.nedges;
        }
    }

    nn = HashSize (cut_hash);
    cutplane.nnodes = (cgsize_t)nn;
    cutplane.nodes = (Node *) MALLOC (funcname, nn * sizeof(Node));
    HashDestroy (cut_hash, get_cut_node);

    return cutplane.nelems;
}

/*------------------------------------------------------------------*/

static void draw_edges ()
{
    int nz, nr;
    cgsize_t ne, nn;
    Zone *z;
    Regn *r;

    glDisable(GL_LIGHTING);
    glShadeModel(GL_FLAT);
    glBegin(GL_LINES);
    glColor3fv(cutcolor);

    for (nz = 0; nz < nzones; nz++) {
        z = &zones[nz];
        for (nr = 0; nr < z->nregs; nr++) {
            r = &z->regs[nr];
            if (r->cut.nedges == 0) continue;
            if (!usecutclr) glColor3fv(r->color);
            for (ne = 0; ne < r->cut.nedges; ne++) {
                nn = r->cut.edges[ne].nodes[0];
                glVertex3fv(cutplane.nodes[nn]);
                nn = r->cut.edges[ne].nodes[1];
                glVertex3fv(cutplane.nodes[nn]);
            }
        }
    }

    glEnd();
}

/*------------------------------------------------------------------*/

static void draw_elements (int mode)
{
    int nz, nr, i, j;
    int ip, nf, nnodes;
    cgsize_t n, ne, nn;
    float *nodes[4], norm[3];
    CGNS_ENUMT(ElementType_t) type;
    Zone *z;
    Regn *r;
    Face *f;

    glEnable(GL_LIGHTING);
    glShadeModel(GL_FLAT);
    glPolygonMode(GL_FRONT_AND_BACK, mode);
    glColor3fv(cutcolor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, cutcolor);

    for (nz = 0; nz < nzones; nz++) {
        z = &zones[nz];
        for (nr = 0; nr < z->nregs; nr++) {
            r = &z->regs[nr];
            if (r->cut.nelems == 0) continue;
            type = r->elemtype;

            if (!usecutclr) {
                glColor3fv(r->color);
                glMaterialfv(GL_FRONT_AND_BACK,
                    GL_AMBIENT_AND_DIFFUSE, r->color);
            }

            if (type == CGNS_ENUMV(NGON_n)) {
                for (ne = 0; ne < r->cut.nelems; ne++) {
                    nn = r->cut.elems[ne];
                    nnodes = r->elems[nn++];
                    if (nnodes < 3) continue;
                    if (nnodes == 3)
                        glBegin(GL_TRIANGLES);
                    else if (nnodes == 4)
                        glBegin(GL_QUADS);
                    else
                        glBegin(GL_POLYGON);
                    glNormal3fv(face_normal(z, nnodes, &r->elems[nn]));
                    for (i = 0; i < nnodes; i++)
                        glVertex3fv(z->nodes[r->elems[nn++]]);
                    glEnd();
                }
            }
            else if (type == CGNS_ENUMV(NFACE_n)) {
                for (ne = 0; ne < r->cut.nelems; ne++) {
                    nn = r->cut.elems[ne];
                    nf = r->elems[nn++];
                    for (j = 0; j < nf; j++, nn++) {
                        f = r->poly[abs(r->elems[nn])-1];
                        if (f->nnodes < 3) continue;
                        if (f->nnodes == 3)
                            glBegin(GL_TRIANGLES);
                        else if (f->nnodes == 4)
                            glBegin(GL_QUADS);
                        else
                            glBegin(GL_POLYGON);
                        if (r->elems[nn] < 0) {
                            for (i = 0; i < 3; i++)
                                norm[i] = -f->normal[i];
                        }
                        else {
                            for (i = 0; i < 3; i++)
                                norm[i] = f->normal[i];
                        }
                        glNormal3fv(norm);
                        for (i = 0; i < f->nnodes; i++)
                            glVertex3fv(z->nodes[f->nodes[i]]);
                        glEnd();
                    }
                }
            }
            else {
                for (ne = 0; ne < r->cut.nelems; ne++) {
                    n = r->cut.elems[ne];
                    if (r->elemtype == CGNS_ENUMV(MIXED))
                        type = (CGNS_ENUMT(ElementType_t))r->elems[n++];
                    switch (type) {
                        case CGNS_ENUMV(TRI_3):
                        case CGNS_ENUMV(TRI_6):
                            ip = 0;
                            nf = 1;
                            break;
                        case CGNS_ENUMV(QUAD_4):
                        case CGNS_ENUMV(QUAD_8):
                        case CGNS_ENUMV(QUAD_9):
                            ip = 1;
                            nf = 1;
                            break;
                        case CGNS_ENUMV(TETRA_4):
                        case CGNS_ENUMV(TETRA_10):
                            ip = 2;
                            nf = 4;
                            break;
                        case CGNS_ENUMV(PYRA_5):
                        case CGNS_ENUMV(PYRA_14):
                            ip = 6;
                            nf = 5;
                            break;
                        case CGNS_ENUMV(PENTA_6):
                        case CGNS_ENUMV(PENTA_15):
                        case CGNS_ENUMV(PENTA_18):
                            ip = 11;
                            nf = 5;
                            break;
                        case CGNS_ENUMV(HEXA_8):
                        case CGNS_ENUMV(HEXA_20):
                        case CGNS_ENUMV(HEXA_27):
                            ip = 16;
                            nf = 6;
                            break;
                        default:
                            ip = 0;
                            nf = 0;
                            break;
                    }
                    for (j = 0; j < nf; j++) {
                        nnodes = facenodes[ip+j][0];
                        for (i = 0; i < nnodes; i++) {
                            nn = r->elems[n+facenodes[ip+j][i+1]];
                            nodes[i] = z->nodes[nn];
                        }
                        if (nnodes == 4) {
                            glBegin(GL_QUADS);
                            glNormal3fv(compute_normal(nodes[0], nodes[1],
                                nodes[2], nodes[3]));
                        }
                        else {
                            glBegin(GL_TRIANGLES);
                            glNormal3fv(compute_normal(nodes[0], nodes[1],
                                nodes[2], NULL));
                        }
                        for (i = 0; i < nnodes; i++)
                            glVertex3fv(nodes[i]);
                        glEnd();
                    }
                }
            }
        }
    }
}

/*---------- OGLcutplane -------------------------------------------
 * create OGL display list for a cutting plane
 *------------------------------------------------------------------*/

static int OGLcutplane (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int mode;
    float plane[4];
    static char slist[33];

    if (argc < 1 || argc > 3) {
        Tcl_SetResult (interp, "usage: OGLcutplane [mode] [plane]",
            TCL_STATIC);
        return TCL_ERROR;
    }

    /* create and return displaylist flag */

    if (argc == 1) {
        if (!CutDL) CutDL = glGenLists (1);
        sprintf (slist, "%d", CutDL);
        Tcl_SetResult (interp, slist, TCL_STATIC);
        return TCL_OK;
    }

    mode = atoi(argv[1]);

    if (argc == 3) {
        int np;
        CONST char **args;
        if (TCL_OK != Tcl_SplitList (interp, argv[2], &np, &args))
            return TCL_ERROR;
        if (np != 4) {
            Tcl_Free ((char *)args);
            Tcl_SetResult (interp, "invalid plane", TCL_STATIC);
            return TCL_ERROR;
        }
        for (np = 0; np < 4; np++)
            plane[np] = (float) atof (args[np]);
        Tcl_Free ((char *)args);
        init_cutplane(plane);
        find_elements();
    }

    if (!CutDL) CutDL = glGenLists (1);

    glNewList (CutDL, GL_COMPILE);
    if (mode && cutplane.nelems) {
        if (mode == 1) {
            if (cutplane.nedges == 0)
                find_intersects();
            draw_edges ();
        }
        else {
            draw_elements (mode > 2 ? GL_FILL : GL_LINE);
        }
    }
    glEndList ();

    sprintf (slist, "%ld %ld", (long)cutplane.nelems, (long)cutplane.nedges);
    Tcl_SetResult (interp, slist, TCL_STATIC);
    return TCL_OK;
}

/*---------- OGLdrawplane ------------------------------------------
 * draw the cutting plane
 *------------------------------------------------------------------*/

static int OGLdrawplane (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int n, np, i, j, k, index, n0, n1;
    CONST char **args;
    float plane[4], bbox[3][2], s[8], ds;
    float node[8][3], pnode[6][3];
    static char slist[17];

    if (argc < 1 || argc > 2) {
        Tcl_SetResult (interp, "usage: OGLdrawplane [plane]",
            TCL_STATIC);
        return TCL_ERROR;
    }

    if (!PlaneDL) PlaneDL = glGenLists (1);

    if (argc == 1) {
      glNewList (PlaneDL, GL_COMPILE);
      glEndList ();
      sprintf (slist, "%d", PlaneDL);
      Tcl_SetResult (interp, slist, TCL_STATIC);
      return TCL_OK;
    }

    if (TCL_OK != Tcl_SplitList (interp, argv[1], &np, &args))
        return TCL_ERROR;
    if (np != 4) {
        Tcl_Free ((char *)args);
        Tcl_SetResult (interp, "invalid plane", TCL_STATIC);
        return TCL_ERROR;
    }
    for (n = 0; n < np; n++)
        plane[n] = (float) atof (args[n]);
    Tcl_Free ((char *)args);

    get_bounds (ignorevis, bbox);
    index = n = 0;
    for (k = 0; k < 2; k++) {
        for (j = 0; j < 2; j++) {
            for (i = 0; i < 2; i++) {
                node[n][0] = bbox[0][(i+j)%2];
                node[n][1] = bbox[1][j];
                node[n][2] = bbox[2][k];
                s[n] = node[n][0] * plane[0] + node[n][1] * plane[1] +
                       node[n][2] * plane[2];
                if (s[n] >= plane[3]) index |= (1 << n);
                n++;
            }
        }
    }
    if (index > 0x7f) index ^= 0xff;
    if (index < 1 || index > HEX_SIZE) {
      Tcl_SetResult (interp, "plane doesn't intersect", TCL_STATIC);
      return TCL_ERROR;
    }

    np = hexCuts[index][0];
    for (n = 0; n < np; n++) {
        j = hexCuts[index][n+1];
        n0 = hexEdges[j][0];
        n1 = hexEdges[j][1];
        ds = s[n1] - s[n0];
        if (ds != 0.0)
            ds = (plane[3] - s[n0]) / ds;
        for (i = 0; i < 3; i++)
            pnode[n][i] = node[n0][i] + ds * (node[n1][i] - node[n0][i]);
    }

    glNewList (PlaneDL, GL_COMPILE);
    glEnable(GL_BLEND);
    glDepthMask(GL_FALSE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_LIGHTING);
    glShadeModel(GL_FLAT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor4fv(cutcolor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, cutcolor);

    glBegin(GL_TRIANGLE_FAN);
    glNormal3fv(plane);
    for (n = 0; n < np; n++)
        glVertex3fv(pnode[n]);
    glEnd();

    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
    glEndList ();

    sprintf (slist, "%d", PlaneDL);
    Tcl_SetResult (interp, slist, TCL_STATIC);
    return TCL_OK;
}

/*---------- OGLcutconfig ------------------------------------------
 * set the cutting plane color and operation
 *------------------------------------------------------------------*/

static int OGLcutconfig (ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
    int n, np;
    CONST char **args;

    if (argc < 2 || argc > 4) {
        Tcl_SetResult (interp, "usage: OGLcutconfig color [usecutclr] [ignorevis]",
            TCL_STATIC);
        return TCL_ERROR;
    }

    if (TCL_OK != Tcl_SplitList (interp, argv[1], &np, &args))
        return TCL_ERROR;
    if (np < 3 || np > 4) {
        Tcl_Free ((char *)args);
        Tcl_SetResult (interp, "invalid color", TCL_STATIC);
        return TCL_ERROR;
    }
    for (n = 0; n < np; n++)
        cutcolor[n] = (float) atof (args[n]);
    Tcl_Free ((char *)args);

    if (argc > 2) {
        usecutclr = atoi (argv[2]);
        if (argc > 3)
            ignorevis = atoi (argv[3]);
    }
    return TCL_OK;
}

#endif

/*---------- Cgnstcl_Init --------------------------------------
 * Initialize and create the commands
 *--------------------------------------------------------------*/

#if defined(_WIN32) && defined(BUILD_DLL)
__declspec(dllexport)
#endif
int Cgnstcl_Init(Tcl_Interp *interp)
{
    global_interp = interp;
    Tcl_CreateCommand (interp, "CGNSopen", (Tcl_CmdProc *)CGNSopen,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGNSclose", (Tcl_CmdProc *)CGNSclose,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGNSbase", (Tcl_CmdProc *)CGNSbase,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGNSzone", (Tcl_CmdProc *)CGNSzone,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGNSsummary", (Tcl_CmdProc *)CGNSsummary,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGNSgetbase", (Tcl_CmdProc *)CGNSgetbase,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGNSgetzone", (Tcl_CmdProc *)CGNSgetzone,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGNSgetregion", (Tcl_CmdProc *)CGNSgetregion,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGNSregiondim", (Tcl_CmdProc *)CGNSregiondim,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "CGNSbounds", (Tcl_CmdProc *)CGNSbounds,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "OGLregion", (Tcl_CmdProc *)OGLregion,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "OGLaxis", (Tcl_CmdProc *)OGLaxis,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "OGLcolor", (Tcl_CmdProc *)OGLcolor,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
#ifndef NO_CUTTING_PLANE
    Tcl_CreateCommand (interp, "OGLcutplane", (Tcl_CmdProc *)OGLcutplane,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "OGLdrawplane", (Tcl_CmdProc *)OGLdrawplane,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
    Tcl_CreateCommand (interp, "OGLcutconfig", (Tcl_CmdProc *)OGLcutconfig,
        (ClientData)0, (Tcl_CmdDeleteProc *)0);
#endif
    return TCL_OK;
}

