/*
 * interpolate_cgns.c - solution interpolation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#ifdef _WIN32
# define unlink _unlink
#else
# include <unistd.h>
#endif

#include "getargs.h"
#include "cgnslib.h"
#include "cgnsutil.h"

#ifndef CG_MODE_READ
# define CG_MODE_READ   MODE_READ
# define CG_MODE_MODIFY MODE_MODIFY
#endif

#define BBOX_PADDING 0.01
#define MAX_ELEMENTS 256
#define MAX_DEPTH    16

#define MAX_ITER     20
#define TOLERANCE    0.00001

/* command line options */

static char options[] = "c:b:B:S:wand:e:p:i:t:";

static char *usgmsg[] = {
    "usage  : interpolate_cgns [options] CGNSsol CGNSgrid [CGNSnew]",
    "         reads solution from CGNSsol and interpolates the",
    "         solution onto the grid from CGNSgrid.",
    "options:",
    "   -c<base>  = use CGNS base number <base> for CGNSsol (default 1)",
    "   -b<base>  = use CGNS base number <base> for CGNSgrid (default 1)",
    "   -B<name>  = write results to base <name> (default same as read)",
    "   -S<name>  = write results to solution <name>",
    "   -w        = use volume weighting",
    "   -a        = allow element extrapolation",
    "   -n        = use nearest point values",
    "   -d<depth> = max depth for octtree (default 16)",
    "   -e<nelem> = max number of elements in octtree branch (default 256)",
    "   -p<pad>   = bounding box padding fraction (default 0.01)",
    "   -i<iter>  = max newton iterations (default 20)",
    "   -t<tol>   = u,v,w tolerance (default 0.00001)",
    NULL
};

typedef struct _Element {
    int flag;
    int zone;
    int nnodes;
    cgsize_t nodes[8];
    double bbox[3][2];
} Element;

static cgsize_t num_elements;
static Element *elements;

typedef struct _OctTree {
    int depth;
    int subflag;
    struct _OctTree *parent;
    struct _OctTree *tree[8];
    double bbox[3][2];
    cgsize_t nelem;
    Element **elem;
} OctTree;

static OctTree root;
static int depths[3];
static cgsize_t counts[4];

static int nbasezones;
static ZONE *basezones;
static char *solname = NULL;
static int weighting = 0;
static int extrapolate = 0;
static int nearestpt = 0;
static int numout = 0;
static int numextrap = 0;
static int numconv = 0;

static int max_elements = MAX_ELEMENTS;
static int max_depth = MAX_DEPTH;
static float bbox_padding = (float)BBOX_PADDING;
static int max_iter = MAX_ITER;
static double tolerance = TOLERANCE;

static char buff[1024];

/*-------------------------------------------------------------------*/

static int sort_name (const void *v1, const void *v2)
{
    return strcmp (((FIELD *)v1)->name, ((FIELD *)v2)->name);
}

/*-------------------------------------------------------------------*/

static void check_solution (void)
{
    int nz, nf;
    ZONE *z;

    for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
        read_solution_field (nz, 1, 0);
        if (z->sols->nflds == 0) {
            sprintf (buff, "missing solution for zone %d", nz);
            FATAL ("check_solution", buff);
        }
        if (z->sols->location == CGNS_ENUMV(CellCenter))
            cell_vertex_solution (nz, 1, weighting);
        if (z->sols->nflds > 1)
            qsort (z->sols->flds, z->sols->nflds, sizeof(FIELD), sort_name);
        if (nz > 1) {
            if (z->sols->nflds != Zones->sols->nflds)
               FATAL ("check_solution",
                   "solution inconsistant between zones");
            for (nf = 0; nf < z->sols->nflds; nf++) {
                if (strcmp (z->sols->flds[nf].name, Zones->sols->flds[nf].name))
                   FATAL ("check_solution",
                       "solution inconsistant between zones");
            }
        }
    }
}

/*-------------------------------------------------------------------*/

static cgsize_t count_elements (int nz)
{
    int ns, et;
    cgsize_t n, nn, ne, nelem = 0;
    ZONE *z = &Zones[nz];

    for (ns = 0; ns < z->nesets; ns++) {
        ne = z->esets[ns].end - z->esets[ns].start + 1;
        et = z->esets[ns].type;
        if (et == CGNS_ENUMV(MIXED)) {
            for (n = 0, nn = 0; nn < ne; nn++) {
                et = (int)z->esets[ns].conn[n++];
                if (et < CGNS_ENUMV(NODE) || et > CGNS_ENUMV(HEXA_27))
                    FATAL ("count_elements", "unrecognized element type");
                if (et >= CGNS_ENUMV(TETRA_4) && et <= CGNS_ENUMV(HEXA_27)) nelem++;
                n += element_node_counts[et];
            }
        }
        else {
            if (et >= CGNS_ENUMV(TETRA_4) && et <= CGNS_ENUMV(HEXA_27)) nelem += ne;
        }
    }
    return nelem;
}

/*-------------------------------------------------------------------*/

static void add_elements (int nz)
{
    int i, nn, ns, et;
    cgsize_t n, j, ne, iv, nelem = 0;
    ZONE *z = &Zones[nz];
    VERTEX *v;
    Element *e = &elements[num_elements];

    for (ns = 0; ns < z->nesets; ns++) {
        ne = z->esets[ns].end - z->esets[ns].start + 1;
        et = z->esets[ns].type;
        if (et < CGNS_ENUMV(TETRA_4) || et > CGNS_ENUMV(MIXED)) continue;
        for (n = 0, j = 0; j < ne; j++) {
            if (z->esets[ns].type == CGNS_ENUMV(MIXED))
                et = (int)z->esets[ns].conn[n++];
            switch (et) {
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
                e[nelem].zone = nz;
                e[nelem].nnodes = nn;
                for (i = 0; i < nn; i++) {
                    iv = z->esets[ns].conn[n+i] - 1;
                    e[nelem].nodes[i] = iv;
                    v = &z->verts[iv];
                    if (i) {
                        if (e[nelem].bbox[0][0] > v->x)
                            e[nelem].bbox[0][0] = v->x;
                        if (e[nelem].bbox[0][1] < v->x)
                            e[nelem].bbox[0][1] = v->x;
                        if (e[nelem].bbox[1][0] > v->y)
                            e[nelem].bbox[1][0] = v->y;
                        if (e[nelem].bbox[1][1] < v->y)
                            e[nelem].bbox[1][1] = v->y;
                        if (e[nelem].bbox[2][0] > v->z)
                            e[nelem].bbox[2][0] = v->z;
                        if (e[nelem].bbox[2][1] < v->z)
                            e[nelem].bbox[2][1] = v->z;
                    }
                    else {
                        e[nelem].bbox[0][0] = e[nelem].bbox[0][1] = v->x;
                        e[nelem].bbox[1][0] = e[nelem].bbox[1][1] = v->y;
                        e[nelem].bbox[2][0] = e[nelem].bbox[2][1] = v->z;
                    }
                }
                nelem++;
            }
            n += element_node_counts[et];
        }
    }
    num_elements += nelem;
}

/*-------------------------------------------------------------------*/

static int contains_point (OctTree *tree, VERTEX *vert)
{
    if (vert->x < tree->bbox[0][0] || vert->x > tree->bbox[0][1] ||
        vert->y < tree->bbox[1][0] || vert->y > tree->bbox[1][1] ||
        vert->z < tree->bbox[2][0] || vert->z > tree->bbox[2][1])
        return 0;
    return 1;
}

#if 0
/*-------------------------------------------------------------------*/

static int contains_edge (OctTree *tree, VERTEX *v1, VERTEX *v2)
{
    int n;
    double s, x, y, z;

    if ((v1->x < tree->bbox[0][0] && v2->x < tree->bbox[0][0]) ||
        (v1->x > tree->bbox[0][1] && v2->x > tree->bbox[0][1]) ||
        (v1->y < tree->bbox[1][0] && v2->y < tree->bbox[1][0]) ||
        (v1->y > tree->bbox[1][1] && v2->y > tree->bbox[1][1]) ||
        (v1->z < tree->bbox[2][0] && v2->z < tree->bbox[2][0]) ||
        (v1->z > tree->bbox[2][1] && v2->z > tree->bbox[2][1]))
        return 0;

    if (v1->x != v2->x) {
        for (n = 0; n < 2; n++) {
            s = (tree->bbox[0][n] - v1->x) / (v2->x - v1->x);
            if (s >= 0.0 && s <= 1.0) {
                y = v1->y + s * (v2->y - v1->y);
                z = v1->z + s * (v2->z - v1->z);
                if (y >= tree->bbox[1][0] && y <= tree->bbox[1][1] &&
                    z >= tree->bbox[2][0] && z <= tree->bbox[2][1])
                    return 1;
            }
        }
    }

    if (v1->y != v2->y) {
        for (n = 0; n < 2; n++) {
            s = (tree->bbox[1][n] - v1->y) / (v2->y - v1->y);
            if (s >= 0.0 && s <= 1.0) {
                x = v1->x + s * (v2->x - v1->x);
                z = v1->z + s * (v2->z - v1->z);
                if (x >= tree->bbox[0][0] && x <= tree->bbox[0][1] &&
                    z >= tree->bbox[2][0] && z <= tree->bbox[2][1])
                    return 1;
            }
        }
    }

    if (v1->z != v2->z) {
        for (n = 0; n < 2; n++) {
            s = (tree->bbox[2][n] - v1->z) / (v2->z - v1->z);
            if (s >= 0.0 && s <= 1.0) {
                x = v1->x + s * (v2->x - v1->x);
                y = v1->y + s * (v2->y - v1->y);
                if (x >= tree->bbox[0][0] && x <= tree->bbox[0][1] &&
                    y >= tree->bbox[1][0] && y <= tree->bbox[1][1])
                    return 1;
            }
        }
    }

    return 0;
}
#endif

/*-------------------------------------------------------------------*/

static int contains_element (OctTree *tree, Element *elem)
{
    int n;
#if 0
    ZONE *z = &Zones[elem->zone];
    VERTEX *v[8];
#endif

    for (n = 0; n < 3; n++) {
        if (elem->bbox[n][0] > tree->bbox[n][1] ||
            elem->bbox[n][1] < tree->bbox[n][0])
            return 0;
    }
#if 0
    for (n = 0; n < elem->nnodes; n++) {
        v[n] = &z->verts[elem->nodes[n]];
        if (contains_point (tree, v[n]))
            return 1;
    }

    if (elem->nnodes == 4) {
        return (contains_edge (tree, v[0], v[1]) ||
                contains_edge (tree, v[0], v[2]) ||
                contains_edge (tree, v[0], v[3]) ||
                contains_edge (tree, v[1], v[2]) ||
                contains_edge (tree, v[1], v[3]) ||
                contains_edge (tree, v[2], v[3]));
    }
    if (elem->nnodes == 5) {
        return (contains_edge (tree, v[0], v[1]) ||
                contains_edge (tree, v[0], v[3]) ||
                contains_edge (tree, v[0], v[4]) ||
                contains_edge (tree, v[1], v[2]) ||
                contains_edge (tree, v[1], v[4]) ||
                contains_edge (tree, v[2], v[3]) ||
                contains_edge (tree, v[2], v[4]) ||
                contains_edge (tree, v[3], v[4]));
    }
    if (elem->nnodes == 6) {
        return (contains_edge (tree, v[0], v[1]) ||
                contains_edge (tree, v[0], v[2]) ||
                contains_edge (tree, v[0], v[3]) ||
                contains_edge (tree, v[1], v[2]) ||
                contains_edge (tree, v[1], v[4]) ||
                contains_edge (tree, v[2], v[5]) ||
                contains_edge (tree, v[3], v[4]) ||
                contains_edge (tree, v[3], v[5]) ||
                contains_edge (tree, v[4], v[5]));
    }
    return (contains_edge (tree, v[0], v[1]) ||
            contains_edge (tree, v[0], v[3]) ||
            contains_edge (tree, v[0], v[4]) ||
            contains_edge (tree, v[1], v[2]) ||
            contains_edge (tree, v[1], v[5]) ||
            contains_edge (tree, v[2], v[3]) ||
            contains_edge (tree, v[2], v[6]) ||
            contains_edge (tree, v[3], v[7]) ||
            contains_edge (tree, v[4], v[5]) ||
            contains_edge (tree, v[4], v[7]) ||
            contains_edge (tree, v[5], v[6]) ||
            contains_edge (tree, v[6], v[7]));
#endif
    /*
     * add if bounding boxes overlap - this does not imply that
     * the element is actually contained in the tree bounding box
     * and may result in additional elements that must be searched.
     * This does, however, guarantee that all possible elements
     * will be included in the tree branch.
     */
    return 1;
}

/*-------------------------------------------------------------------*/

static void subdivide (OctTree *parent)
{
    int i, j, k, n, ne;
    double bbox[3][3];
    OctTree *tree;

    if (parent->nelem <= max_elements || parent->depth >= max_depth) {
        if (depths[0] > parent->depth) depths[0] = parent->depth;
        if (depths[1] < parent->depth) depths[1] = parent->depth;
        depths[2] += parent->depth;
        if (counts[0] > parent->nelem) counts[0] = parent->nelem;
        if (counts[1] < parent->nelem) counts[1] = parent->nelem;
        counts[2] += parent->nelem;
        (counts[3])++;
        return;
    }

    parent->subflag = 1;
    for (n = 0; n < 3; n++) {
        bbox[n][0] = parent->bbox[n][0];
        bbox[n][1] = 0.5 * (parent->bbox[n][0] + parent->bbox[n][1]);
        bbox[n][2] = parent->bbox[n][1];
    }

    for (n = 0; n < 8; n++) {
        parent->tree[n] = tree = (OctTree *) malloc (sizeof(OctTree));
        if (NULL == tree)
            FATAL ("subdivide", "malloc failed for octtree branch");
        i = (n & 1);
        j = (n & 2) >> 1;
        k = (n & 4) >> 2;
        tree->bbox[0][0] = bbox[0][i];
        tree->bbox[0][1] = bbox[0][i+1];
        tree->bbox[1][0] = bbox[1][j];
        tree->bbox[1][1] = bbox[1][j+1];
        tree->bbox[2][0] = bbox[2][k];
        tree->bbox[2][1] = bbox[2][k+1];

        tree->depth = parent->depth + 1;
        tree->subflag = 0;
        tree->parent = parent;
        for (ne = 0, i = 0; i < parent->nelem; i++) {
            if (contains_element (tree, parent->elem[i])) {
                parent->elem[i]->flag = 1;
                ne++;
            }
            else
                parent->elem[i]->flag = 0;
        }

        tree->nelem = ne;
        if (ne) {
            tree->elem = (Element **) malloc (ne * sizeof(Element *));
            if (NULL == tree->elem)
                FATAL ("subdivide", "malloc failed for octtree elements");
            for (ne = 0, i = 0; i < parent->nelem; i++) {
                if (parent->elem[i]->flag)
                    tree->elem[ne++] = parent->elem[i];
            }
        }
    }
    free (parent->elem);
    parent->nelem = 0;

    for (n = 0; n < 8; n++)
        subdivide (parent->tree[n]);
}

/*-------------------------------------------------------------------*/

static void build_octree (void)
{
    int i, nz;
    cgsize_t n, ne;
    double diff;
    ZONE *z = Zones;

    root.depth = 0;
    root.subflag = 0;
    root.parent = NULL;

    root.bbox[0][0] = root.bbox[0][1] = z->verts->x;
    root.bbox[1][0] = root.bbox[1][1] = z->verts->y;
    root.bbox[2][0] = root.bbox[2][1] = z->verts->z;
    for (ne = 0, nz = 0; nz < nZones; nz++, z++) {
        for (n = 0; n < z->nverts; n++) {
            if (root.bbox[0][0] > z->verts[n].x)
                root.bbox[0][0] = z->verts[n].x;
            if (root.bbox[0][1] < z->verts[n].x)
                root.bbox[0][1] = z->verts[n].x;
            if (root.bbox[1][0] > z->verts[n].y)
                root.bbox[1][0] = z->verts[n].y;
            if (root.bbox[1][1] < z->verts[n].y)
                root.bbox[1][1] = z->verts[n].y;
            if (root.bbox[2][0] > z->verts[n].z)
                root.bbox[2][0] = z->verts[n].z;
            if (root.bbox[2][1] < z->verts[n].z)
                root.bbox[2][1] = z->verts[n].z;
        }
        ne += count_elements (nz);
    }
    if (!ne) FATAL ("build_octree", "no volume elements found");

    /* add buffer around bounding box */

    for (i = 0; i < 3; i++) {
        diff = bbox_padding * (root.bbox[i][1] - root.bbox[i][0]);
        root.bbox[i][0] -= diff;
        root.bbox[i][1] += diff;
    }

    /* build element list */

    elements = (Element *) malloc ((size_t)ne * sizeof(Element));
    if (NULL == elements) FATAL ("build_octree", "malloc failed for elements");
    num_elements = 0;
    for (nz = 0; nz < nZones; nz++)
        add_elements (nz);
    if (num_elements != ne) FATAL ("build_octree", "mismatch in element count");

    root.nelem = num_elements;
    root.elem = (Element **) malloc ((size_t)num_elements * sizeof(Element *));
    if (NULL == root.elem)
        FATAL ("build_octree", "malloc failed for element pointers");
    for (ne = 0; ne < num_elements; ne++)
        root.elem[ne] = &elements[ne];

    depths[0] = max_depth;
    counts[0] = num_elements;
    for (i = 1; i < 3; i++) {
        depths[i] = 0;
        counts[i] = 0;
    }
    counts[3] = 0;

    subdivide (&root);
}

/*-------------------------------------------------------------------*/

#define SWAP(A,B) {temp=(A);(A)=(B);(B)=temp;}

static int invert3x3 (double a[3][3], double b[3])
{
    int i, j, k, irow = 0, icol = 0;
    int indxc[3], indxr[3], ipiv[3];
    double big, temp, pivinv;

    for (j = 0; j < 3; j++)
        ipiv[j] = 0;
    for (i = 0; i < 3; i++) {
        big = 0.0;
        for (j = 0; j < 3; j++) {
            if (ipiv[j] != 1) {
                for (k = 0; k < 3; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs (a[j][k]) >= big) {
                            big = fabs (a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else {
                        if (ipiv[k] > 1) return 1;
                    }
                }
            }
        }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (j = 0; j < 3; j++)
                SWAP (a[irow][j], a[icol][j]);
            SWAP (b[irow], b[icol]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0) return 2;
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for (j = 0; j < 3; j++) a[icol][j] *= pivinv;
        b[icol] *= pivinv;
        for (k = 0; k < 3; k++) {
            if (k != icol) {
                temp = a[k][icol];
                a[k][icol] = 0.0;
                for (j = 0; j < 3; j++) a[k][j] -= (a[icol][j] * temp);
                b[k] -= (b[icol] * temp);
            }
        }
    }
    for (j = 2; j >= 0; j--) {
        if (indxr[j] != indxc[j]) {
            for (k = 0; k < 3; k++)
                SWAP (a[k][indxr[j]], a[k][indxc[j]]);
        }
    }
    return 0;
}

/*-------------------------------------------------------------------*/

static void compute_shapef (int nnodes, double uvw[3], double shapef[8],
                            double deriv[8][3])
{
    if (nnodes == 4) {
        shapef[0] = 1.0 - uvw[0] - uvw[1] - uvw[2];
        shapef[1] = uvw[0];
        shapef[2] = uvw[1];
        shapef[3] = uvw[2];
    }
    else if (nnodes == 5) {
        shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
        shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
        shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
        shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
        shapef[4] = uvw[2];
    }
    else if (nnodes == 6) {
        shapef[0] = (1.0 - uvw[0] - uvw[1]) * (1.0 - uvw[2]);
        shapef[1] = uvw[0] * (1.0 - uvw[2]);
        shapef[2] = uvw[1] * (1.0 - uvw[2]);
        shapef[3] = (1.0 - uvw[0] - uvw[1]) * uvw[2];
        shapef[4] = uvw[0] * uvw[2];
        shapef[5] = uvw[1] * uvw[2];
    }
    else if (nnodes == 8) {
        shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
        shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
        shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
        shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
        shapef[4] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * uvw[2];
        shapef[5] = uvw[0] * (1.0 - uvw[1]) * uvw[2];
        shapef[6] = uvw[0] * uvw[1] * uvw[2];
        shapef[7] = (1.0 - uvw[0]) * uvw[1] * uvw[2];
    }
    else
        FATAL ("compute_shapef", "invalid number of nodes for element");

    if (deriv != NULL) {
        if (nnodes == 4) {
            deriv[0][0] = -1.0;
            deriv[0][1] = -1.0;
            deriv[0][2] = -1.0;
            deriv[1][0] =  1.0;
            deriv[1][1] =  0.0;
            deriv[1][2] =  0.0;
            deriv[2][0] =  0.0;
            deriv[2][1] =  1.0;
            deriv[2][2] =  0.0;
            deriv[3][0] =  0.0;
            deriv[3][1] =  0.0;
            deriv[3][2] =  1.0;
        }
        else if (nnodes == 5) {
            deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
            deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
            deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
            deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
            deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
            deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
            deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
            deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
            deriv[2][2] = -uvw[0] * uvw[1];
            deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
            deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
            deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
            deriv[4][0] =  0.0;
            deriv[4][1] =  0.0;
            deriv[4][2] =  1.0;
        }
        else if (nnodes == 6) {
            deriv[0][0] = -(1.0 - uvw[2]);
            deriv[0][1] = -(1.0 - uvw[2]);
            deriv[0][2] = -(1.0 - uvw[0] - uvw[1]);
            deriv[1][0] =  (1.0 - uvw[2]);
            deriv[1][1] =  0.0;
            deriv[1][2] = -uvw[0];
            deriv[2][0] =  0.0;
            deriv[2][1] =  (1.0 - uvw[2]);
            deriv[2][2] = -uvw[1];
            deriv[3][0] = -uvw[2];
            deriv[3][1] = -uvw[2];
            deriv[3][2] =  (1.0 - uvw[0] - uvw[1]);
            deriv[4][0] =  uvw[2];
            deriv[4][1] =  0.0;
            deriv[4][2] =  uvw[0];
            deriv[5][0] =  0.0;
            deriv[5][1] =  uvw[2];
            deriv[5][2] =  uvw[1];
        }
        else {
            deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
            deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
            deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
            deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
            deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
            deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
            deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
            deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
            deriv[2][2] = -uvw[0] * uvw[1];
            deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
            deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
            deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
            deriv[4][0] = -(1.0 - uvw[1]) * uvw[2];
            deriv[4][1] = -(1.0 - uvw[0]) * uvw[2];
            deriv[4][2] =  (1.0 - uvw[0]) * (1.0 - uvw[1]);
            deriv[5][0] =  (1.0 - uvw[1]) * uvw[2];
            deriv[5][1] = -uvw[0] * uvw[2];
            deriv[5][2] =  uvw[0] * (1.0 - uvw[1]);
            deriv[6][0] =  uvw[1] * uvw[2];
            deriv[6][1] =  uvw[0] * uvw[2];
            deriv[6][2] =  uvw[0] * uvw[1];
            deriv[7][0] = -uvw[1] * uvw[2];
            deriv[7][1] =  (1.0 - uvw[0]) * uvw[2];
            deriv[7][2] =  (1.0 - uvw[0]) * uvw[1];
        }
    }
}

/*-------------------------------------------------------------------*/

static int compute_uvw (Element *elem, VERTEX *pt, double uvw[3])
{
    int i, j, n;
    double dist;
    double a[3][3], b[3], shapef[8], dw[8][3];
    VERTEX *v[8];
    ZONE *z = &basezones[elem->zone];

    /* check if point is within element bounding box */

    if (!extrapolate) {
        if (pt->x < elem->bbox[0][0] || pt->x > elem->bbox[0][1] ||
            pt->y < elem->bbox[1][0] || pt->y > elem->bbox[1][1] ||
            pt->z < elem->bbox[2][0] || pt->z > elem->bbox[2][1])
            return 0;
    }

    for (n = 0; n < elem->nnodes; n++)
        v[n] = &z->verts[elem->nodes[n]];

    /* for tetrahedron, direct solution */

    if (elem->nnodes == 4) {
        for (i = 0; i < 3; i++) {
            a[0][i] = v[i+1]->x - v[0]->x;
            a[1][i] = v[i+1]->y - v[0]->y;
            a[2][i] = v[i+1]->z - v[0]->z;
        }
        b[0] = pt->x - v[0]->x;
        b[1] = pt->y - v[0]->y;
        b[2] = pt->z - v[0]->z;
        if (invert3x3 (a, b)) return 0;
        for (n = 0; n < 3; n++)
            uvw[n] = b[n];
        if (!extrapolate) {
            for (n = 0; n < 3; n++) {
                if (uvw[n] < 0.0) {
                    if (uvw[n] < -tolerance) return 0;
                    uvw[n] = 0.0;
                }
                else if (uvw[n] > 1.0) {
                    if (uvw[n] > 1.0+tolerance) return 0;
                    uvw[n] = 1.0;
                }
            }
        }
        return 1;
    }

    /* all other elements, use Newton iteration */

    for (n = 0; n < 3; n++)
        uvw[n] = 0.5;
    for (n = 0; n < max_iter; n++) {
        compute_shapef (elem->nnodes, uvw, shapef, dw);
        b[0] = -pt->x;
        b[1] = -pt->y;
        b[2] = -pt->z;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                a[i][j] = 0.0;
            }
        }
        for (i = 0; i < elem->nnodes; i++) {
            b[0] += (shapef[i] * v[i]->x);
            b[1] += (shapef[i] * v[i]->y);
            b[2] += (shapef[i] * v[i]->z);
            for (j = 0; j < 3; j++) {
                a[0][j] -= (dw[i][j] * v[i]->x);
                a[1][j] -= (dw[i][j] * v[i]->y);
                a[2][j] -= (dw[i][j] * v[i]->z);
            }
        }
        if (invert3x3 (a, b)) return 0;
        for (dist = 0.0, i = 0; i < 3; i++) {
            dist += b[i] * b[i];
            uvw[i] += b[i];
        }
        if (dist <= tolerance * tolerance) {
            if (!extrapolate) {
                for (i = 0; i < 3; i++) {
                    if (uvw[i] < 0.0) {
                        if (uvw[i] < -tolerance) return 0;
                        uvw[i] = 0.0;
                    }
                    else if (uvw[i] > 1.0) {
                        if (uvw[i] > 1.0+tolerance) return 0;
                        uvw[i] = 1.0;
                    }
                }
            }
            return 1;
        }
    }
#if 0
    printf ("%g (%g,%g,%g) (%g,%g,%g)\n", dist, b[0], b[1], b[2],
        uvw[0], uvw[1], uvw[2]);
#endif
    numconv++;
    return 0;
}

/*-------------------------------------------------------------------*/

static Element *closest_point (OctTree *tree, VERTEX *pt, double shapef[8])
{
    int ne, nn, index = 0, nt, ntree;
    double dist, min_dist = 1.0e32;
    OctTree **tp;
    Element *e, *elem = NULL;
    VERTEX *v;
    ZONE *z;

    if (tree->nelem == 0) {
        tp = tree->parent->tree;
        ntree = 8;
    }
    else {
        tp = &tree;
        ntree = 1;
    }
    for (nt = 0; nt < ntree; nt++) {
        for (ne = 0; ne < tp[nt]->nelem; ne++) {
            e = tp[nt]->elem[ne];
            z = &basezones[e->zone];
            for (nn = 0; nn < e->nnodes; nn++) {
                v = &z->verts[e->nodes[nn]];
                v->id = 0;
            }
        }
    }
    for (nt = 0; nt < ntree; nt++) {
        for (ne = 0; ne < tp[nt]->nelem; ne++) {
            e = tp[nt]->elem[ne];
            z = &basezones[e->zone];
            for (nn = 0; nn < e->nnodes; nn++) {
                v = &z->verts[e->nodes[nn]];
                if (!v->id) {
                    v->id = 1;
                    dist = (pt->x - v->x) * (pt->x - v->x) +
                           (pt->y - v->y) * (pt->y - v->y) +
                           (pt->z - v->z) * (pt->z - v->z);
                    if (elem == NULL || dist < min_dist) {
                        elem = e;
                        min_dist = dist;
                        index = nn;
                    }
                }
            }
        }
    }
    for (nn = 0; nn < elem->nnodes; nn++)
        shapef[nn] = 0.0;
    shapef[index] = 1.0;
    numout++;
    return elem;
}

/*-------------------------------------------------------------------*/

static Element *find_element (VERTEX *pt, double shapef[8])
{
    int n, ne, nt, ntree;
    double xm, ym, zm, dist, min_dist;
    double uvw[3], min_uvw[3];
    OctTree **tp, *tree = &root;
    Element *elem = NULL;

    /* find branch containing point */

    if (!contains_point (tree, pt)) {
        sprintf (buff, "vertex %g,%g,%g outside bounds",
            pt->x, pt->y, pt->z);
        FATAL ("find_element", buff);
    }
    while (tree->subflag) {
        xm = 0.5 * (tree->bbox[0][0] + tree->bbox[0][1]);
        ym = 0.5 * (tree->bbox[1][0] + tree->bbox[1][1]);
        zm = 0.5 * (tree->bbox[2][0] + tree->bbox[2][1]);
        n = (pt->x < xm ? 0 : 1) |
            (pt->y < ym ? 0 : 2) |
            (pt->z < zm ? 0 : 4);
        tree = tree->tree[n];
    }

    if (nearestpt)
        return closest_point (tree, pt, shapef);

    /* find element containing point */

    if (tree->nelem == 0 && extrapolate) {
        tp = tree->parent->tree;
        ntree = 8;
    }
    else {
        tp = &tree;
        ntree = 1;
    }
    for (nt = 0; nt < ntree; nt++) {
        for (ne = 0; ne < tp[nt]->nelem; ne++) {
            if (compute_uvw (tp[nt]->elem[ne], pt, uvw)) {
                for (dist = 0.0, n = 0; n < 3; n++)
                    dist += (uvw[n] - 0.5) * (uvw[n] - 0.5);
                if (elem == NULL || dist < min_dist) {
                    elem = tp[nt]->elem[ne];
                    min_dist = dist;
                    for (n = 0; n < 3; n++)
                        min_uvw[n] = uvw[n];
                }
            }
        }
    }

    /* if no element was found use nearest point */

    if (elem == NULL)
        elem = closest_point (tree, pt, shapef);
    else
        compute_shapef (elem->nnodes, min_uvw, shapef, NULL);
    return elem;
}

/*-------------------------------------------------------------------*/

static void build_solution (int nz)
{
    int n, ns, nv, nf;
    SOLUTION *sol = NULL;
    ZONE *bz, *z = &Zones[nz];
    Element *elem;
    double wsum, fsum, shapef[8];

    if (z->nsols) {
        if (solname == NULL)
            sol = z->sols;
        else {
            for (ns = 0; ns < z->nsols; ns++) {
                if (!strcmp (solname, z->sols[ns].name)) {
                    sol = &z->sols[ns];
                    break;
                }
            }
            if (sol == NULL) {
                z->sols = (SOLUTION *) realloc (z->sols,
                    (z->nsols + 1) * sizeof(SOLUTION));
                if (NULL == z->sols)
                    FATAL ("build_solution",
                        "realloc failed for a new solution");
                sol = &z->sols[(z->nsols)++];
            }
        }
    }
    else {
        z->nsols = 1;
        z->sols = sol = new_solution (1);
    }

    if (solname != NULL) {
        strncpy (sol->name, solname, 32);
        sol->name[32] = 0;
    }
    sol->location = CGNS_ENUMV(Vertex);
    sol->size = z->nverts;
    for (n = 0; n < 5; n++)
        sol->units[n] = basezones->sols->units[n];
    sol->dataclass = basezones->sols->dataclass;
    sol->nflds = basezones->sols->nflds;
    sol->flds = new_field (sol->nflds, z->nverts);

    for (nf = 0; nf < sol->nflds; nf++) {
        strcpy (sol->flds[nf].name, basezones->sols->flds[nf].name);
        sol->flds[nf].datatype = basezones->sols->flds[nf].datatype;
        for (n = 0; n < 5; n++)
            sol->flds[nf].units[n] = basezones->sols->flds[nf].units[n];
        sol->flds[nf].dataclass = basezones->sols->flds[nf].dataclass;
        sol->flds[nf].convtype = basezones->sols->flds[nf].convtype;
        for (n = 0; n < 2; n++)
        sol->flds[nf].dataconv[n] = basezones->sols->flds[nf].dataconv[n];
        sol->flds[nf].exptype = basezones->sols->flds[nf].exptype;
        for (n = 0; n < 5; n++)
            sol->flds[nf].exponent[n] = basezones->sols->flds[nf].exponent[n];
    }

    for (nv = 0; nv < z->nverts; nv++) {
        elem = find_element (&z->verts[nv], shapef);
        for (n = 0; n < elem->nnodes; n++) {
            if (shapef[n] < 0.0 || shapef[n] > 1.0) {
                numextrap++;
                break;
            }
        }
        bz = &basezones[elem->zone];
        for (wsum = 0.0, n = 0; n < elem->nnodes; n++) {
            shapef[n] *= bz->verts[elem->nodes[n]].w;
            wsum += shapef[n];
        }
        if (wsum == 0.0) wsum = 1.0;
        for (nf = 0; nf < sol->nflds; nf++) {
            for (fsum = 0.0, n = 0; n < elem->nnodes; n++)
                fsum += shapef[n] * bz->sols->flds[nf].data[elem->nodes[n]];
            sol->flds[nf].data[nv] = fsum / wsum;
        }
    }
}

/*-------------------------------------------------------------------*/

#if 0

/* this is for debugging - solution and grid files must have same nodes */

static void compare_solution ()
{
    int nz, nf, n;
    double f, diff, fmax, favg, fsum;

    for (nz = 0; nz < nZones; nz++) {
        for (nf = 0; nf < Zones[nz].sols->nflds; nf++) {
            f = basezones[nz].sols->flds[nf].data[0];
            favg = fmax = fabs (Zones[nz].sols->flds[nf].data[0] - f);
            fsum = fabs (f);
            for (n = 1; n < Zones[nz].sols->size; n++) {
                f = basezones[nz].sols->flds[nf].data[n];
                diff = fabs (Zones[nz].sols->flds[nf].data[n] - f);
                if (fmax < diff) fmax = diff;
                favg += diff;
                fsum += fabs (f);
            }
            favg /= (double)n;
            fsum /= (double)n;
            printf ("%d %s %g %g %g\n", nz+1,
                Zones[nz].sols->flds[nf].name, fsum, favg, fmax);
        }
    }
}

#endif

/*-------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
    int n, nz, dim;
    int base1 = 1, base2 = 1;
    char *tmpfile, *newbase = NULL, basename[33];

    if (argc < 2)
        print_usage (usgmsg, NULL);

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'c':
                base1 = atoi (argarg);
                break;
            case 'b':
                base2 = atoi (argarg);
                break;
            case 'B':
                newbase = argarg;
                break;
            case 'S':
                solname = argarg;
                break;
            case 'a':
                extrapolate = 1;
                break;
            case 'w':
                weighting = 1;
                break;
            case 'n':
                nearestpt = 1;
                break;
            case 'd':
                max_depth = atoi (argarg);
                break;
            case 'e':
                max_elements = atoi (argarg);
                break;
            case 'p':
                bbox_padding = (float)atof (argarg);
                break;
            case 'i':
                max_iter = atoi (argarg);
                break;
            case 't':
                tolerance = atof (argarg);
                break;
        }
    }

    if (argind >= argc + 1)
        print_usage (usgmsg, "CGNSbase and/or CGNSfile not given");
    if (!file_exists (argv[argind]) || !file_exists (argv[argind+1]))
        FATAL (NULL, "CGNSbase and/or CGNSfile do not exist or is not a file");

    /* read solution CGNS file */

    printf ("reading CGNS solution file from %s\n", argv[argind]);
    if (cg_open (argv[argind], CG_MODE_READ, &cgnsfn) ||
        cg_base_read (cgnsfn, base1, basename, &dim, &dim))
        FATAL (NULL, NULL);
    cgnsbase = base1;
    printf ("reading zone information for base %d - %s\n",
        cgnsbase, basename);

    read_zones ();
    for (nz = 1; nz <= nZones; nz++) {
        read_zone_grid (nz);
        if (Zones[nz-1].type == CGNS_ENUMV(Structured))
            structured_elements (nz);
        else
            read_zone_element (nz);
        if (weighting)
            vertex_volumes (nz);
        read_zone_solution (nz);
        if (!Zones[nz-1].nsols) {
            sprintf (buff, "zone %d does not contain a solution", nz);
            FATAL (NULL, buff);
        }
    }

    puts ("checking zone solutions...");
    fflush (stdout);
    check_solution ();

    puts ("building octtree...");
    fflush (stdout);
    build_octree ();
    printf ("               min     max     avg\n");
    printf (" depth:   %8d%8d%8d\n",
        depths[0], depths[1], depths[2] / (int)counts[3]);
    printf (" elements:%8ld%8ld%8ld\n",
        (long)counts[0], (long)counts[1], (long)(counts[2] / counts[3]));

    /* save zone information */
    nbasezones = nZones;
    basezones = Zones;
    cg_close (cgnsfn);

    /* create a working copy */

    printf ("\ncreating a working copy of %s\n", argv[++argind]);
    tmpfile = temporary_file (argv[argind]);
    copy_file (argv[argind], tmpfile);

    printf ("reading CGNS file from %s\n", tmpfile);
    if (cg_open (tmpfile, CG_MODE_MODIFY, &cgnsfn) ||
        cg_base_read (cgnsfn, base2, basename, &dim, &dim))
        FATAL (NULL, NULL);
    cgnsbase = base2;

    /* conversion may leave temp file (older CGNS versions) */

    sprintf (buff, "%s.temp", tmpfile);
    unlink (buff);

    printf ("reading zone information for base %d - %s\n",
        cgnsbase, basename);
    read_cgns ();
    printf ("interpolating solution using %s averaging...\n",
        weighting ? "volume" : "simple");
    for (nz = 0; nz < nZones; nz++) {
        printf ("  zone %d, %ld vertices...\n", nz+1, (long)Zones[nz].nverts);
        fflush (stdout);
        build_solution (nz);
    }
#if 0
    if (numconv)
        printf ("Newton iteration failed for %d points\n", numconv);
#endif
    if (numout)
        printf ("%d points were set to the nearest existing point\n", numout);
    if (numextrap)
        printf ("%d points were extrapolated from element values\n", numextrap);

    if (newbase != NULL) {
        strncpy (basename, newbase, 32);
        basename[32] = 0;
        if (cg_base_write (cgnsfn, basename, 3, 3, &cgnsbase))
            FATAL (NULL, NULL);
    }
    printf ("writing data to base %d - %s...\n", cgnsbase, basename);
    fflush (stdout);
    write_cgns ();
    cg_close (cgnsfn);

    if (argind + 1 < argc) argind++;
    printf ("renaming %s to %s\n", tmpfile, argv[argind]);
    unlink (argv[argind]);
    if (rename (tmpfile, argv[argind])) {
        char msg[512];
        sprintf (msg, "rename %s -> %s failed", tmpfile, argv[argind]);
        FATAL (NULL, msg);
        exit (1);
    }
    free (tmpfile);
    return 0;
}
