/*
 * import routines for CGNS
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "cgnsImport.h"
#include "hash.h"

#ifndef CG_MODE_READ
# define CG_MODE_READ   MODE_READ
# define CG_MODE_WRITE  MODE_WRITE
# define CG_MODE_MODIFY MODE_MODIFY
#endif

#ifndef DOUBLE      /* data type for coordinates */
#define DOUBLE      double
#endif

#ifndef TOLERANCE   /* tolerance for duplicate checking */
#define TOLERANCE   1.e-03
#endif

#ifndef REGION_BASE /* base name for unnamed regions */
#define REGION_BASE "Region"
#endif

/*--- variables ---*/

static int num_vars = 0;
static char **var_names = 0;

/*--- node bit flags ---*/

#define USED_BIT    1
#define REGN_BIT    2
#define REFS_BIT    4

/*--- node structure ---*/

typedef struct {
    cgsize_t id;    /* node ID */
    int flags;      /* references to node */
    DOUBLE x, y, z; /* coordinates */
    DOUBLE dist;    /* distance from origin */
    DOUBLE *vars;   /* variables */
} cgnsNODE;

/*--- node mapping data ---*/

#define NODEMAP_INC 50  /* realloc this many at a time */

typedef struct {
    cgsize_t nodeid;    /* node id */
    cgsize_t mapped;    /* set when mapped */
    cgnsNODE *node;     /* pointer to node data */
} NODEMAP;

static cgsize_t num_nodes = 0; /* number of nodes */
static cgsize_t max_nodes = 0; /* number of nodes malloced */
static NODEMAP *nodemap;       /* list of nodes */

/*--- duplicate node checking ---*/

#define ENTRY_INC   50

static int no_check = 0;
static double def_tol = TOLERANCE;
static double tolerance = TOLERANCE;
DOUBLE xmin, xmax, ymin, ymax, zmin, zmax;

static cgsize_t num_entries = 0;
static cgsize_t max_entries = 0;
static cgnsNODE **node_list;

/*--- element data ---*/

#define ELEMENT_INC 50  /* realloc this many at a time */

typedef struct {
    cgsize_t elemid;  /* element ID */
    int elemtype;     /* element type (number of nodes) */
    cgsize_t *nodeid; /* node ID's for element */
    char facemap[6];  /* remapping of faces */
} cgnsELEM;

static cgsize_t num_elements = 0;/* number of elements */
static cgsize_t max_elements = 0;/* number of elements malloced */
static cgnsELEM *elemlist;       /* list of elements */

/*--- region data ---*/

#define REGION_INC  50          /* step increment for realloc */

static char region_name[33];      /* region name */
static int region_type;           /* type of region */
static cgsize_t region_id = 0;    /* region ID */
static cgsize_t region_max = 0;   /* malloced size of region_list */
static cgsize_t region_nobjs = 0; /* number objects in region */
static cgsize_t *region_list;     /* list of nodes in region */

typedef struct {
    char name[33];   /* region name */
    int type;        /* region type */
    cgsize_t nobjs;  /* number of objects */
    cgsize_t *objid; /* object ID's */
} cgnsREGN;

static int num_regions = 0; /* number of regions */
static cgnsREGN *reglist;     /* list of regions */

/*--- external faces */

typedef struct  {
    cgsize_t faceid;
    int nnodes;
    cgsize_t nodeid[4];
    int flags;
} cgnsFACE;

static cgsize_t num_faces = 0;
static cgnsFACE **facelist;

/*--- CGNS data ---*/

int cgnsFile = 0;
int cgnsBase = 0;
int cgnsZone = 0;
char cgnsZoneName[33] = "";

/*--- error handling callback ---*/

static void (*errcallback)( /* callback pointer to user routine */
    char *errmsg            /* error message */
) = NULL;

/*======================================================================
 * Node routines
 *======================================================================*/

/*---------- NewNode ---------------------------------------------------
 * create a new node
 *----------------------------------------------------------------------*/

static cgnsNODE *NewNode (cgsize_t id, double x, double y, double z)
{
    cgnsNODE *node = (cgnsNODE *) malloc (sizeof(cgnsNODE));

    if (NULL != node) {
        node->id = id;
        node->flags = 0;
        node->x = x;
        node->y = y;
        node->z = z;
        node->dist = 0.0;
        node->vars = 0;
        if (num_vars)
            node->vars = (DOUBLE *) calloc (num_vars, sizeof(DOUBLE));
    }
    return (node);
}

/*---------- GetNode ---------------------------------------------------
 * return the node for a given node id
 *----------------------------------------------------------------------*/

static cgnsNODE *GetNode (cgsize_t nodeid, cgsize_t *pos)
{
    cgsize_t lo = 0, hi = num_nodes - 1;

    *pos = 0;
    if (!num_nodes || nodeid < nodemap[0].nodeid)
        return (NULL);
    if (nodeid == nodemap[0].nodeid)
        return (nodemap[0].node);
    if (!hi || nodeid > nodemap[hi].nodeid) {
        *pos = num_nodes;
        return (NULL);
    }
    if (nodeid == nodemap[hi].nodeid) {
        *pos = hi;
        return (nodemap[hi].node);
    }

    while (1) {
        *pos = (lo + hi) >> 1;
        if (nodeid == nodemap[*pos].nodeid)
            return (nodemap[*pos].node);
        if (hi - lo <= 1)
            break;
        if (nodeid < nodemap[*pos].nodeid)
            hi = *pos;
        else
            lo = *pos;
    }
    *pos = hi;
    return (NULL);
}

/*---------- CompareNodes ----------------------------------------------
 * compare two nodes, returns 0 if nodes are the same within
 * the specifed tolerance, else 1
 *----------------------------------------------------------------------*/

static int CompareNodes (cgnsNODE *node1, cgnsNODE *node2)
{
    double dist = (node2->x - node1->x) * (node2->x - node1->x) +
                  (node2->y - node1->y) * (node2->y - node1->y) +
                  (node2->z - node1->z) * (node2->z - node1->z);
    return (dist < (tolerance * tolerance) ? 0 : 1);
}

/*======================================================================
  duplicate node checking routines

The nodes are stored in a sorted list based on radius from the origin.
Once the position in the list is determined, then a scan backwards and
forwards in the list is done for a matching node.
========================================================================*/

/*---------- FindPosition ----------------------------------------------
 * bisection search to locate position for node
 *----------------------------------------------------------------------*/

static cgsize_t FindPosition (cgnsNODE *node)
{
    cgsize_t mid, lo = 0, hi = num_entries - 1;

    if (!num_entries || node->dist <= node_list[0]->dist)
        return (0);

    if (!hi || node->dist > node_list[hi]->dist)
        return (num_entries);
    if (node->dist == node_list[hi]->dist)
        return (hi);

    while ((hi - lo) > 1) {
        mid = (lo + hi) >> 1;
        if (node->dist == node_list[mid]->dist)
            return (mid);
        if (node->dist < node_list[mid]->dist)
            hi = mid;
        else
            lo = mid;
    }
    return (hi);
}

/*---------- FindNode --------------------------------------------------
 * search for matching node in Node List
 *----------------------------------------------------------------------*/

static cgnsNODE *FindNode (cgnsNODE *node, cgsize_t *pos)
{
    cgsize_t n;

    *pos = FindPosition (node);

    for (n = *pos - 1; n >= 0; n--) {
        if (fabs (node->dist - node_list[n]->dist) >= tolerance)
            break;
        if (!CompareNodes (node, node_list[n]))
            return (node_list[n]);
    }
    for (n = *pos; n < num_entries; n++) {
        if (fabs (node->dist - node_list[n]->dist) >= tolerance)
            break;
        if (!CompareNodes (node, node_list[n]))
            return (node_list[n]);
    }
    return (NULL);
}

/*---------- AddNode ---------------------------------------------------
 * add a new node to the duplicate node checking list
 *----------------------------------------------------------------------*/

static void AddNode (cgnsNODE *node, cgsize_t pos)
{
    cgsize_t n;

    if (num_entries == max_entries) {
        if (!max_entries)
            node_list = (cgnsNODE **) malloc (ENTRY_INC * sizeof(cgnsNODE *));
        else
            node_list = (cgnsNODE **) realloc (node_list,
                (size_t)(max_entries + ENTRY_INC) * sizeof(cgnsNODE *));
        if (NULL == node_list)
            cgnsImportFatal (
            "AddNode:malloc failed for new node entry in duplicate node list");
        max_entries += ENTRY_INC;
    }
    for (n = num_entries; n > pos; n--)
        node_list[n] = node_list[n-1];
    node_list[pos] = node;
    num_entries++;
}

/*======================================================================
 * Element routines
 *======================================================================*/

/*---------- GetElement ------------------------------------------------
 * return the element for a given element id
 *----------------------------------------------------------------------*/

static cgnsELEM *GetElement (cgsize_t elemid, cgsize_t *pos)
{
    cgsize_t lo = 0, hi = num_elements - 1;

    *pos = 0;
    if (!num_elements || elemid < elemlist[0].elemid)
        return (NULL);
    if (elemid == elemlist[0].elemid)
        return (&elemlist[0]);
    if (!hi || elemid > elemlist[hi].elemid) {
        *pos = num_elements;
        return (NULL);
    }
    if (elemid == elemlist[hi].elemid) {
        *pos = hi;
        return (&elemlist[hi]);
    }

    while (1) {
        *pos = (lo + hi) >> 1;
        if (elemid == elemlist[*pos].elemid)
            return (&elemlist[*pos]);
        if (hi - lo <= 1)
            break;
        if (elemid < elemlist[*pos].elemid)
            hi = *pos;
        else
            lo = *pos;
    }
    *pos = hi;
    return (NULL);
}

/*---------- NewElement ------------------------------------------------
 * add new element to list of elements
 *----------------------------------------------------------------------*/

static cgnsELEM *NewElement (cgsize_t pos)
{
    int i;
    cgsize_t n;

    /* malloc/realloc if needed */

    if (num_elements == max_elements) {
        if (!max_elements)
            elemlist = (cgnsELEM *) malloc (ELEMENT_INC * sizeof(cgnsELEM));
        else
            elemlist = (cgnsELEM *) realloc (elemlist,
                (size_t)(max_elements + ELEMENT_INC) * sizeof(cgnsELEM));
        if (NULL == elemlist)
            cgnsImportFatal ("AddElement:malloc failed for element list");
        max_elements += ELEMENT_INC;
    }

    /* insert new element */

    for (n = num_elements; n > pos; n--) {
        elemlist[n].elemid   = elemlist[n-1].elemid;
        elemlist[n].elemtype = elemlist[n-1].elemtype;
        elemlist[n].nodeid   = elemlist[n-1].nodeid;
        for (i = 0; i < 6; i++)
            elemlist[n].facemap[i] = elemlist[n-1].facemap[i];
    }
    num_elements++;
    return (&elemlist[pos]);
}

/*======================================================================
 * Region routines
 *======================================================================*/

/*---------- GetRegion -------------------------------------------------
 * return a region for a given region name
 *----------------------------------------------------------------------*/

static cgnsREGN *GetRegion (char *name, int *pos)
{
    int cmp, lo = 0, hi = num_regions - 1;

    *pos = 0;
    if (!num_regions || (cmp = strcmp (name, reglist[0].name)) < 0)
        return (NULL);
    if (0 == cmp)
        return (&reglist[0]);
    if (!hi || (cmp = strcmp (name, reglist[hi].name)) > 0) {
        *pos = num_regions;
        return (NULL);
    }
    if (0 == cmp) {
        *pos = hi;
        return (&reglist[hi]);
    }

    while (1) {
        *pos = (lo + hi) >> 1;
        if (0 == (cmp = strcmp (name, reglist[*pos].name)))
            return (&reglist[*pos]);
        if (hi - lo <= 1)
            break;
        if (cmp < 0)
            hi = *pos;
        else
            lo = *pos;
    }
    *pos = hi;
    return (NULL);
}

/*---------- NewRegion -------------------------------------------------
 * add a new region to region list
 *----------------------------------------------------------------------*/

static cgnsREGN *NewRegion (char *name, int pos)
{
    int n;
    static char *errmsg = "NewRegion:malloc failed for region list";

    if (!num_regions)
        reglist = (cgnsREGN *) malloc (sizeof(cgnsREGN));
    else
        reglist = (cgnsREGN *) realloc (reglist,
            (num_regions + 1) * sizeof(cgnsREGN));
    if (NULL == reglist)
        cgnsImportFatal (errmsg);
    for (n = num_regions; n > pos; n--) {
        strcpy (reglist[n].name, reglist[n-1].name);
        reglist[n].type  = reglist[n-1].type;
        reglist[n].nobjs = reglist[n-1].nobjs;
        reglist[n].objid = reglist[n-1].objid;
    }
    strncpy (reglist[pos].name, name, 32);
    reglist[pos].name[32] = 0;
    reglist[pos].type  = 0;
    reglist[pos].nobjs = 0;
    reglist[pos].objid = NULL;
    num_regions++;
    return (&reglist[pos]);
}

/*======================================================================
 * external face regions
 *======================================================================*/

/*---------- get_face_nodes -----------------------------------------
 * get nodes for an element face
 *-------------------------------------------------------------------*/

static int get_face_nodes (cgsize_t faceid, cgsize_t nodeid[4])
{
    cgnsELEM *elem;
    cgsize_t elemid = faceid >> 3;
    int facenum = (int)(faceid & 7);
    int n, nfaces = 0, noff = 0, nnodes;
    static int facenodes[20][5] = {
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

    if (elemid < 0 || elemid >= num_elements)
        cgnsImportFatal ("get_face_nodes:invalid element number");
    elem = &elemlist[elemid];
    switch (elem->elemtype) {
        case cgnsELEM_TET:
            noff = 0;
            nfaces = 4;
            break;
        case cgnsELEM_PYR:
            noff = 4;
            nfaces = 5;
            break;
        case cgnsELEM_WDG:
            noff = 9;
            nfaces = 5;
            break;
        case cgnsELEM_HEX:
            noff = 14;
            nfaces = 6;
            break;
        default:
            cgnsImportFatal ("get_face_nodes:invalid element type");
    }
    if (facenum < 0 || facenum >= nfaces)
        return (0);
    noff += (int)elem->facemap[facenum];
    nnodes = facenodes[noff][0];
    for (n = 0; n < nnodes; n++)
        nodeid[n] = elem->nodeid[facenodes[noff][n+1]];
    return (nnodes);
}

/*---------- compare_faces -------------------------------------------
 * face comparison routine
 *--------------------------------------------------------------------*/

static int compare_faces (void *v1, void *v2)
{
    cgnsFACE *f1 = (cgnsFACE *)v1;
    cgnsFACE *f2 = (cgnsFACE *)v2;
    int n;

    if (f1->nnodes != f2->nnodes)
        return (f1->nnodes - f2->nnodes);

    /* the following assumes nodes have been sorted */

    for (n = 0; n < f1->nnodes; n++) {
        if (f1->nodeid[n] != f2->nodeid[n])
            return (int)(f1->nodeid[n] - f2->nodeid[n]);
    }
    return (0);
}

/*---------- get_faces ----------------------------------------------
 * get the exterior faces
 *-------------------------------------------------------------------*/

static void get_faces (void *v)
{
    facelist[num_faces++] = (cgnsFACE *)v;
}

/*---------- hash_face -----------------------------------------------
 * face hash function
 *--------------------------------------------------------------------*/

static size_t hash_face (void *v)
{
    cgnsFACE *face = (cgnsFACE *)v;
    int n;
    size_t hash = 0;

    for (n = 0; n < face->nnodes; n++)
        hash += (size_t)face->nodeid[n];
    return (hash);
}

/*---------- sortfaces -------------------------------------------------
 * called by qsort to sort the list of faces
 *----------------------------------------------------------------------*/

static int sortfaces (const void *f1, const void *f2)
{
    return (int)((*((cgnsFACE **)f1))->faceid - (*((cgnsFACE **)f2))->faceid);
}

/*---------- exterior_faces -----------------------------------------
 * find exterior faces
 *-------------------------------------------------------------------*/

static void exterior_faces (void)
{
    int i, j, k, nfaces;
    cgsize_t nodeid[4];
    cgsize_t n, nn, id, faceid;
    HASH FaceList;
    cgnsFACE *pf, face;

    FaceList = HashCreate (2047, compare_faces, hash_face);
    if (NULL == FaceList)
        cgnsImportFatal ("exterior_faces:malloc failed for face hash table");

    for (n = 0; n < num_elements; n++) {
        switch (elemlist[n].elemtype) {
            case cgnsELEM_WDG:
                nfaces = 5;
                break;
            case cgnsELEM_HEX:
                nfaces = 6;
                break;
            default:
                nfaces = elemlist[n].elemtype;
                break;
        }

        /* loop over element faces */

        for (j = 0; j < nfaces; j++) {

            /* get face nodes and sort */

            faceid = (n << 3) | j;
            face.nnodes = get_face_nodes (faceid, nodeid);
            for (i = 0; i < face.nnodes; i++) {
                id = nodeid[i];
                for (k = 0; k < i; k++) {
                    if (face.nodeid[k] > id) {
                        nn = face.nodeid[k];
                        face.nodeid[k] = id;
                        id = nn;
                    }
                }
                face.nodeid[i] = id;
            }

            if (NULL == (pf = (cgnsFACE *) HashFind (FaceList, &face))) {

                /* create new face and add to list */

                if (NULL == (pf = (cgnsFACE *) malloc (sizeof(cgnsFACE))))
                    cgnsImportFatal ("exterior_faces:malloc failed for new face");
                pf->faceid = faceid;
                pf->flags = 0;
                pf->nnodes = face.nnodes;
                for (i = 0; i < face.nnodes; i++)
                    pf->nodeid[i] = face.nodeid[i];
                (void) HashAdd (FaceList, pf);
            }

            /* else already exists */

            else {
                HashDelete (FaceList, pf);
                free (pf);
            }
        }
    }

    facelist = (cgnsFACE **) malloc (HashSize (FaceList) * sizeof(cgnsFACE *));
    if (NULL == facelist)
        cgnsImportFatal ("exterior_faces:malloc failed for exterior face list");
    num_faces = 0;
    HashDestroy (FaceList, get_faces);

    /* check if faces need sorting */

    for (n = 1; n < num_faces; n++) {
        if (facelist[n]->faceid < facelist[n-1]->faceid)
            break;
    }
    if (n < num_faces)
        qsort (facelist, (size_t)num_faces, sizeof(cgnsFACE *), sortfaces);

    /* get face nodes in the correct order */

    for (n = 0; n < num_faces; n++) {
        get_face_nodes (facelist[n]->faceid, nodeid);
        for (i = 0; i < 4; i++)
            facelist[n]->nodeid[i] = nodeid[i];
    }
}

/*===================================================================
 * write regions to cgns file
 *===================================================================*/

/*---------- sortnodes -------------------------------------------------
 * called by qsort to sort list of node ID mappings
 *----------------------------------------------------------------------*/

static int sortnodes (const void *n1, const void *n2)
{
    return (int)(*((cgsize_t *)n1) - *((cgsize_t *)n2));
}

/*---------- write_node_region --------------------------------------
 * write region from node list
 *-------------------------------------------------------------------*/

static cgsize_t write_node_region (cgnsREGN *reg, cgsize_t offset)
{
    int nn, isect;
    cgsize_t i, j, mid, lo, hi, pos;
    cgsize_t nfaces, nc, *conns;
    CGNS_ENUMT(ElementType_t) elemtype = CGNS_ENUMV(ElementTypeNull);
    cgnsNODE *node;

    /* get exterior faces */

    if (num_faces == 0) exterior_faces ();
    for (j = 0; j < num_faces; j++)
        facelist[j]->flags = 0;

    /* sort region nodes */

    for (i = 1; i < reg->nobjs; i++) {
        if (reg->objid[i] < reg->objid[i-1])
            break;
    }
    if (i < reg->nobjs)
        qsort (reg->objid, (size_t)reg->nobjs, sizeof(cgsize_t), sortnodes);

    /* scan list of exterior faces */

    nfaces = nc = 0;
    for (j = 0; j < num_faces; j++) {
        if (facelist[j]->flags) continue;
        for (nn = 0; nn < facelist[j]->nnodes; nn++) {
            lo = 0;
            hi = reg->nobjs - 1;
            while (lo <= hi) {
                mid = (lo + hi) >> 1;
                if (facelist[j]->nodeid[nn] == reg->objid[mid])
                    break;
                if (facelist[j]->nodeid[nn] < reg->objid[mid])
                    hi = mid - 1;
                else
                    lo = mid + 1;
            }
            if (lo > hi)
                break;
        }
        if (nn == facelist[j]->nnodes) {
            nfaces++;
            facelist[j]->flags = 1;
            if (nc != nn) {
                elemtype = nc ? CGNS_ENUMV(MIXED) :
                           (nn == 3 ? CGNS_ENUMV(TRI_3) : CGNS_ENUMV(QUAD_4));
                nc = nn;
            }
        }
    }
    if (!nfaces) return 0;

    conns = (cgsize_t *) malloc ((size_t)(5 * nfaces) * sizeof(cgsize_t));
    if (NULL == conns)
        cgnsImportFatal ("write_node_region:malloc failed for connectivity");

    /* write face connectivities */

    for (nc = 0, j = 0; j < num_faces; j++) {
        if (facelist[j]->flags) {
            if (elemtype == CGNS_ENUMV(MIXED))
                conns[nc++] = facelist[j]->nnodes == 3 ?
                              CGNS_ENUMV(TRI_3) : CGNS_ENUMV(QUAD_4);
            for (nn = 0; nn < facelist[j]->nnodes; nn++) {
                if (NULL == (node = GetNode (facelist[j]->nodeid[nn], &pos)))
                    cgnsImportFatal ("write_node_region:missing element node");
                conns[nc++] = pos + 1;
            }
        }
    }

    if (cg_section_write (cgnsFile, cgnsBase, cgnsZone, reg->name,
            elemtype, offset, offset + nfaces - 1, 0, conns, &isect))
        cgnsImportFatal ((char *)cg_get_error());

    /* create parent cell mapping */

    for (nc = 0, j = 0; j < num_faces; j++) {
        if (facelist[j]->flags)
            conns[nc++] = (facelist[j]->faceid >> 3) + 1;
    }
    for (j = 0; j < nfaces; j++)
        conns[nc++] = 0;
    for (j = 0; j < num_faces; j++) {
        if (facelist[j]->flags)
            conns[nc++] = (facelist[j]->faceid & 7) + 1;
    }
    for (j = 0; j < nfaces; j++)
        conns[nc++] = 0;
    if (cg_parent_data_write (cgnsFile, cgnsBase, cgnsZone, isect, conns))
        cgnsImportFatal ((char *)cg_get_error());

    free (conns);
    return nfaces;
}

/*---------- write_face_region --------------------------------------
 * write region from face list
 *-------------------------------------------------------------------*/

static cgsize_t write_face_region (cgnsREGN *reg, cgsize_t offset)
{
    int nn, facenum, i, isect;
    cgsize_t elemid, nodeid[4];
    cgsize_t n, nc, pos, *conns;
    CGNS_ENUMT(ElementType_t) elemtype = CGNS_ENUMV(ElementTypeNull);
    cgnsELEM *elem;
    cgnsNODE *node;

    if (!reg->nobjs) return 0;
    conns = (cgsize_t *) malloc ((size_t)(5 * reg->nobjs) * sizeof(cgsize_t));
    if (NULL == conns)
        cgnsImportFatal ("write_face_region:malloc failed for connectivity");

    for (i = 0, n = 0; n < reg->nobjs; n++) {
        elemid = reg->objid[n] >> 3;
        facenum = (int)(reg->objid[n] & 7) - 1;
        if (NULL == (elem = GetElement (elemid, &pos)))
            cgnsImportFatal ("write_face_region:region element not found");
        nn = get_face_nodes ((pos << 3) | facenum, nodeid);
        if (i != nn) {
            if (i) {
                elemtype = CGNS_ENUMV(MIXED);
                break;
            }
            i = nn;
            elemtype = nn == 3 ? CGNS_ENUMV(TRI_3) : CGNS_ENUMV(QUAD_4);
        }
    }

    for (nc = 0, n = 0; n < reg->nobjs; n++) {
        elemid = reg->objid[n] >> 3;
        facenum = (int)(reg->objid[n] & 7) - 1;
        elem = GetElement (elemid, &pos);
        nn = get_face_nodes ((pos << 3) | facenum, nodeid);
        if (elemtype == CGNS_ENUMV(MIXED))
            conns[nc++] = nn == 3 ? CGNS_ENUMV(TRI_3) : CGNS_ENUMV(QUAD_4);
        for (i = 0; i < nn; i++) {
            if (NULL == (node = GetNode (nodeid[i], &pos)))
                cgnsImportFatal ("write_face_region:missing element node");
            conns[nc++] = pos + i;
        }
    }

    if (cg_section_write (cgnsFile, cgnsBase, cgnsZone, reg->name,
            elemtype, offset, offset + reg->nobjs - 1, 0, conns, &isect))
        cgnsImportFatal ((char *)cg_get_error());

    free (conns);
    return reg->nobjs;
}

/*---------- write_elem_region --------------------------------------
 * write elements as region
 *-------------------------------------------------------------------*/

static cgsize_t write_elem_region (cgnsREGN *reg, cgsize_t offset)
{
    return 0;
}

/*======================================================================
 * API routines
 *======================================================================*/

/*---------- cgnsImportError -------------------------------------------
 * setup error handler call back
 *----------------------------------------------------------------------*/

void cgnsImportError (void (*callback)(char *msg))
{
    errcallback = callback;
}

/*---------- cgnsImportFatal -------------------------------------------
 * write error message and exit
 *----------------------------------------------------------------------*/

void cgnsImportFatal (char *errmsg)
{
    if (NULL != errcallback)
        (*errcallback) (errmsg);
    else if (NULL != errmsg && *errmsg)
        fprintf (stderr, "%s\n", errmsg);
    exit (-1);
}

/*---------- cgnsImportSetTol ------------------------------------------
 * setup tolerance for duplicate node checking
 *----------------------------------------------------------------------*/

double cgnsImportSetTol (double tol)
{
    tolerance = tol >= 0.0 ? tol : TOLERANCE;
    tol = def_tol;
    def_tol = tolerance;
    return (tol);
}

/*---------- cgnsImportGetTol ------------------------------------------
 * return tolerance for duplicate node checking
 *----------------------------------------------------------------------*/

double cgnsImportGetTol (int rel)
{
    double tol = def_tol;

    if (rel && num_nodes) {
        double avgvol = (xmax-xmin) * (ymax-ymin) * (zmax-zmin) /
            (DOUBLE)num_nodes;
        tol *= pow (avgvol, 0.33333);
    }
    return (tol);
}

/*---------- cgnsImportSetCheck ----------------------------------------
 * set duplicate node checking on/off
 *----------------------------------------------------------------------*/

void cgnsImportSetCheck (int set)
{
    no_check = set;
}

/*---------- cgnsImportRange --------------------------------------------
 * gets bounding box of node coordinates
 *-----------------------------------------------------------------------*/

cgsize_t cgnsImportRange (double *x1, double *y1, double *z1,
                     double *x2, double *y2, double *z2)
{
    *x1 = xmin; *y1 = ymin; *z1 = zmin;
    *x2 = xmax; *y2 = ymax; *z2 = zmax;
    return num_nodes;
}

/*---------- cgnsImportCheck --------------------------------------------
 * check for and remove duplicate nodes
 *-----------------------------------------------------------------------*/

cgsize_t cgnsImportCheck (int rel)
{
    cgsize_t n, pos, dup_cnt = 0;
    cgnsNODE *node;

    if (num_nodes < 2)
        return (0);

    /* set tolerance */

    tolerance = cgnsImportGetTol (rel);

    /* scan list of nodes, and remove duplicates */

    for (n = 0; n < num_nodes; n++) {
        if (!nodemap[n].mapped) {
            nodemap[n].node->dist = sqrt ((double)(
                (nodemap[n].node->x - xmin) * (nodemap[n].node->x - xmin) +
                (nodemap[n].node->y - ymin) * (nodemap[n].node->y - ymin) +
                (nodemap[n].node->z - zmin) * (nodemap[n].node->z - zmin)));
            node = FindNode (nodemap[n].node, &pos);
            if (NULL == node)
                AddNode (nodemap[n].node, pos);
            else if (node != nodemap[n].node) {
                if (REFS_BIT == (nodemap[n].node->flags & REFS_BIT)) {
                    for (pos = 0; pos < num_nodes; pos++) {
                        if (nodemap[pos].mapped &&
                            nodemap[pos].node == nodemap[n].node)
                            nodemap[pos].node = node;
                    }
                }
                node->flags |=
                    ((nodemap[n].node->flags & USED_BIT) | REFS_BIT);
                free (nodemap[n].node);
                nodemap[n].node = node;
                dup_cnt++;
            }
            nodemap[n].mapped = 1;
        }
    }

    /* free duplicate node checking list */

    free (node_list);
    num_entries = max_entries = 0;
    return (dup_cnt);
}

/*---------- cgnsImportMap ----------------------------------------------
 * map a node explictly to another
 *-----------------------------------------------------------------------*/

cgsize_t cgnsImportMap (cgsize_t nodeid, cgsize_t mapid)
{
    cgsize_t p1, p2, ret;
    cgnsNODE *n1 = GetNode (nodeid, &p1);
    cgnsNODE *n2 = GetNode (mapid, &p2);

    if (NULL == n1 || NULL == n2)
        return (0);
    if (n1 == n2)
        return (n1->id);
    ret = CompareNodes (n1, n2) ? -1 : n1->id;
    if (REFS_BIT == (n2->flags & REFS_BIT)) {
        cgsize_t n;
        for (n = 0; n < num_nodes; n++) {
            if (nodemap[n].node == n2) {
                nodemap[n].node = n1;
                nodemap[n].mapped = 1;
            }
        }
    }
    else {
        nodemap[p2].node = n1;
        nodemap[p2].mapped = 1;
    }
    n1->flags |= ((n2->flags & USED_BIT) | REFS_BIT);
    free (n2);
    return (ret);
}

/*---------- cgnsImportNode ---------------------------------------------
 * import a node
 *-----------------------------------------------------------------------*/

cgsize_t cgnsImportNode (cgsize_t nodeid, double x, double y, double z)
{
    cgsize_t n, pos;

    if (nodeid <= 0)
        return (0);

    /* find min/max bounds */

    if (!num_nodes) {
        xmin = xmax = x;
        ymin = ymax = y;
        zmin = zmax = z;
    }
    else {
        if (xmin > x) xmin = x;
        if (xmax < x) xmax = x;
        if (ymin > y) ymin = y;
        if (ymax < y) ymax = y;
        if (zmin > z) zmin = z;
        if (zmax < z) zmax = z;
    }

    /* find position to place node id */

    if (NULL != GetNode (nodeid, &pos)) {
        nodemap[pos].node->x = (DOUBLE)x;
        nodemap[pos].node->y = (DOUBLE)y;
        nodemap[pos].node->z = (DOUBLE)z;
        return (-1);
    }

    /* malloc/realloc if needed */

    if (num_nodes == max_nodes) {
        if (!max_nodes)
            nodemap = (NODEMAP *) malloc (NODEMAP_INC * sizeof(NODEMAP));
        else
            nodemap = (NODEMAP *) realloc (nodemap,
                (size_t)(max_nodes + NODEMAP_INC) * sizeof(NODEMAP));
        if (NULL == nodemap)
            cgnsImportFatal (
                "cgnsImportNode:malloc failed for node mapping data");
        max_nodes += NODEMAP_INC;
    }

    /* insert new node */

    for (n = num_nodes; n > pos; n--) {
        nodemap[n].nodeid = nodemap[n-1].nodeid;
        nodemap[n].mapped = nodemap[n-1].mapped;
        nodemap[n].node = nodemap[n-1].node;
    }
    nodemap[pos].nodeid = nodeid;
    nodemap[pos].mapped = no_check;
    nodemap[pos].node = NewNode (nodeid, x, y, z);
    if (NULL == nodemap[pos].node)
        cgnsImportFatal ("cgnsImportNode:malloc failed for a new node");
    num_nodes++;
    return (nodeid);
}

/*---------- cgnsImportSetNode ------------------------------------------
 * set node coordinates for node ID
 *-----------------------------------------------------------------------*/

cgsize_t cgnsImportSetNode (cgsize_t nodeid, double x, double y, double z)
{
    cgsize_t n;
    cgnsNODE *node = GetNode (nodeid, &n);

    if (NULL != node) {
        node->x = (DOUBLE)x;
        node->y = (DOUBLE)y;
        node->z = (DOUBLE)z;
        return (node->id);
    }
    return (0);
}

/*---------- cgnsImportGetNode ------------------------------------------
 * return node coordinates for node ID
 *-----------------------------------------------------------------------*/

cgsize_t cgnsImportGetNode (cgsize_t nodeid, double *x, double *y, double *z)
{
    cgsize_t n;
    cgnsNODE *node = GetNode (nodeid, &n);

    if (NULL != node) {
        *x = node->x;
        *y = node->y;
        *z = node->z;
        return (node->id);
    }
    return (0);
}

/*---------- cgnsImportNodeList -----------------------------------------
 * return list of all node ID's
 *-----------------------------------------------------------------------*/

cgsize_t *cgnsImportNodeList (void)
{
    cgsize_t n, *nodeids;
    
    nodeids = (cgsize_t *) malloc ((size_t)(num_nodes + 1) * sizeof(cgsize_t));
    if (NULL == nodeids)
        cgnsImportFatal ("cgnsImportNodeList:malloc failed for node ID list");
    nodeids[0] = num_nodes;
    for (n = 0; n < num_nodes; n++)
        nodeids[n+1] = nodemap[n].nodeid;
    return (nodeids);
}

/*---------- cgnsImportAddVariable -------------------------------------
 * create a node variable
 *----------------------------------------------------------------------*/

int cgnsImportAddVariable (char *varname)
{
    int n;
    cgnsNODE *node;

    if (varname == NULL || !*varname) return -1;
    for (n = 0; n < num_vars; n++) {
        if (0 == strcmp(varname, var_names[n])) return n;
    }
    if (num_vars)
        var_names = (char **) realloc (var_names, (num_vars + 1) * sizeof(char *));
    else
        var_names = (char **) malloc (sizeof(char *));
    if (var_names == NULL)
        cgnsImportFatal ("AddVariable:malloc/realloc failed for variable name list");
    var_names[num_vars] = (char *) malloc (strlen(varname) + 1);
    if (var_names[num_vars] == NULL)
        cgnsImportFatal ("AddVariable:malloc failed for variable name");
    strcpy(var_names[num_vars++], varname);

    for (n = 0; n < num_entries; n++) {
        node = node_list[n];
        if (num_vars == 1)
            node->vars = (DOUBLE *) malloc (sizeof(DOUBLE));
        else
            node->vars = (DOUBLE *) realloc (node->vars, num_vars * sizeof(DOUBLE));
        if (node->vars == NULL)
            cgnsImportFatal ("AddVariable:malloc failed for node variables");
        node->vars[num_vars-1] = 0.0;
    }
    return num_vars-1;
}

/*---------- cgnsImportGetVariable -------------------------------------
 * get the variable number for a node variable
 *----------------------------------------------------------------------*/

int cgnsImportGetVariable (char *varname)
{
    int n;

    if (varname != NULL && *varname) {
        for (n = 0; n < num_vars; n++) {
            if (0 == strcmp(varname, var_names[n])) return n;
        }
    }
    return -1;
}

/*---------- cgnsImportVariable ----------------------------------------
 * set the value of a variable at a node
 *----------------------------------------------------------------------*/

cgsize_t cgnsImportVariable (cgsize_t nodeid, int varnum, double val)
{
    cgsize_t n;
    cgnsNODE *node = GetNode (nodeid, &n);

    if (NULL == node || varnum < 0 || varnum >= num_vars) return 0;
    node->vars[varnum] = (DOUBLE)val;
    return node->id;
}

/*---------- cgnsImportElement ------------------------------------------
 * import an element
 *-----------------------------------------------------------------------*/

int cgnsImportElement (cgsize_t elemid, int elemtype, cgsize_t *nodelist)
{
    int n, ret;
    cgsize_t pos;
    cgnsNODE *node;
    cgnsELEM *elem;

    if (elemid <= 0 ||
       (elemtype != cgnsELEM_TET && elemtype != cgnsELEM_PYR &&
        elemtype != cgnsELEM_WDG && elemtype != cgnsELEM_HEX))
        return (0);

    /* element not found */

    if (NULL == (elem = GetElement (elemid, &pos))) {
        ret = elemtype;
        elem = NewElement (pos);
    }

    /* element already exists */

    else {
        ret = -1;
        free (elem->nodeid);
    }

    /* set element values */

    elem->elemid   = elemid;
    elem->elemtype = elemtype;
    elem->nodeid   = (cgsize_t *) malloc (elemtype * sizeof(cgsize_t));
    if (NULL == elem->nodeid)
        cgnsImportFatal (
            "cgnsImportElement:malloc failed for a new element");

    for (n = 0; n < elemtype; n++) {
        if (NULL == (node = GetNode (nodelist[n], &pos))) {
            char errmsg[50];
            sprintf (errmsg, "cgnsImportElement:element node %ld not found",
                (long)nodelist[n]);
            cgnsImportFatal (errmsg);
        }
        elem->nodeid[n] = node->id;
        node->flags |= USED_BIT;
    }
    for (n = 0; n < 6; n++)
        elem->facemap[n] = n;

    return (ret);
}

/*---------- cgnsImportGetElement ---------------------------------------
 * return element for element ID
 *-----------------------------------------------------------------------*/

int cgnsImportGetElement (cgsize_t elemid, cgsize_t nodeid[])
{
    int n;
    cgsize_t pos;
    cgnsELEM *elem = GetElement (elemid, &pos);

    if (NULL != elem) {
        for (n = 0; n < elem->elemtype; n++)
            nodeid[n] = elem->nodeid[n];
        return (elem->elemtype);
    }
    return (0);
}

/*---------- cgnsImportElementList --------------------------------------
 * return list of all element ID's
 *-----------------------------------------------------------------------*/

cgsize_t *cgnsImportElementList (void)
{
    cgsize_t n, *elemids;
    
    elemids = (cgsize_t *) malloc ((size_t)(num_elements + 1) * sizeof(int));
    if (NULL == elemids)
        cgnsImportFatal (
            "cgnsImportElementList:malloc failed for element ID list");
    elemids[0] = num_elements;
    for (n = 0; n < num_elements; n++)
        elemids[n+1] = elemlist[n].elemid;
    return (elemids);
}

/*---------- cgnsImportGetFace ------------------------------------------
 * return element face node ID's
 *-----------------------------------------------------------------------*/

int cgnsImportGetFace (cgsize_t elemid, int facenum, cgsize_t nodeid[])
{
    int nfaces;
    cgsize_t pos;
    cgnsELEM *elem = GetElement (elemid, &pos);

    if (NULL == elem)
        return (-1);
    switch (elem->elemtype) {
        case cgnsELEM_WDG:
            nfaces = 5;
            break;
        case cgnsELEM_HEX:
            nfaces = 6;
            break;
        default:
            nfaces = elem->elemtype;
            break;
    }
    if (--facenum < 0 || facenum >= nfaces)
        return (0);
    return get_face_nodes ((pos << 3) | facenum, nodeid);
}

/*---------- cgnsImportFindFace -----------------------------------------
 * return element face number given face node ID's
 *-----------------------------------------------------------------------*/

int cgnsImportFindFace (cgsize_t elemid, int nnodes, cgsize_t nodeid[])
{
    int i, j, nfaces = 0, noff = 0, mask = 0;
    cgsize_t pos;
    cgnsELEM *elem = GetElement (elemid, &pos);
    static int facemask[4][6] = {
        /* tet */
        { 7,  11,  14,  13,   0,   0},
        /* pyramid */
        {15,  19,  22,  28,  25,   0},
        /* wedge */
        {27,  54,  45,   7,  56,   0},
        /* hex */
        {15,  51, 102, 204, 153, 240}
    };

    if (NULL == elem || NULL == nodeid)
        return (-1);

    switch (elem->elemtype) {
        case cgnsELEM_TET:
            if (nnodes != 3)
                return (-1);
            noff = 0;
            nfaces = 4;
            break;
        case cgnsELEM_PYR:
            if (nnodes < 3 || nnodes > 4)
                return (-1);
            noff = 1;
            nfaces = 5;
            break;
        case cgnsELEM_WDG:
            if (nnodes < 3 || nnodes > 4)
                return (-1);
            noff = 2;
            nfaces = 5;
            break;
        case cgnsELEM_HEX:
            if (nnodes != 4)
                return (-1);
            noff = 3;
            nfaces = 6;
            break;
        default:
            cgnsImportFatal ("cgnsImportFindFace:invalid element type");
    }

    for (j = 0; j < nnodes; j++) {
        for (i = 0; i < elem->elemtype; i++) {
            if (nodeid[j] == elem->nodeid[i])
                break;
        }
        if (i == elem->elemtype)
            return (0);
        mask |= (1 << i);
    }
    for (i = 0; i < nfaces; i++) {
        if (mask == facemask[noff][i]) {
            for (j = 0; j < 6; j++) {
                if (i == (int)elem->facemap[j])
                    return (j + 1);
            }
        }
    }
    return (0);
}

/*---------- cgnsImportBegReg --------------------------------------------
 * begin a region specification
 *-----------------------------------------------------------------------*/

cgsize_t cgnsImportBegReg (char *regname, int regtype)
{
    int n;
    cgnsREGN *reg;

    if (region_id)
        cgnsImportEndReg ();

    /* initialize region node list */

    if (0 == region_max) {
        region_list = (cgsize_t *) malloc (REGION_INC * sizeof(cgsize_t));
        if (NULL == region_list)
            cgnsImportFatal (
                "cgnsImportBegReg:malloc failed for region node list");
        region_max = REGION_INC;
    }

    /* initialize region data */

    region_id = num_regions + 1;
    region_type = regtype;
    if (NULL == regname || !*regname)
        sprintf (region_name, "%s%ld", REGION_BASE, (long)region_id);
    else {
        strncpy (region_name, regname, sizeof(region_name));
        region_name[sizeof(region_name)-1] = 0;
    }
    region_nobjs = 0;

    if (NULL == (reg = GetRegion (region_name, &n)))
        return (0);
    if (reg->type != regtype)
        cgnsImportFatal ("cgnsImportBegReg:only 1 type allowed for a region");
    return (reg->nobjs);
}

/*---------- cgnsImportAddReg -------------------------------------------
 * add nodes to the region
 *-----------------------------------------------------------------------*/

cgsize_t cgnsImportAddReg (cgsize_t numobjs, cgsize_t *objlist)
{
    cgsize_t n, pos;
    char errmsg[50];

    if (!region_id)
        cgnsImportFatal ("cgnsImportAddReg:region not defined");

    /* realloc region list array if needed */

    if (region_nobjs + numobjs > region_max) {
        n = region_nobjs + numobjs - region_max;
        if (n < REGION_INC) n = REGION_INC;
        region_list = (cgsize_t *) realloc (region_list,
            (size_t)(region_max + n) * sizeof(cgsize_t));
        if (NULL == region_list)
            cgnsImportFatal (
                "cgnsImportAddReg:malloc failed for region node list");
        region_max += n;
    }

    /* node region */

    if (region_type == cgnsREG_NODES) {
        cgnsNODE *node;
        for (n = 0; n < numobjs; n++) {
            if (NULL == (node = GetNode (objlist[n], &pos))) {
                sprintf (errmsg, "cgnsImportAddReg:region node %ld not found",
                    (long)objlist[n]);
                cgnsImportFatal (errmsg);
            }
            region_list[region_nobjs++] = node->id;
            node->flags |= USED_BIT;
        }
    }

    /* face region */

    else if (region_type == cgnsREG_FACES) {
        int facenum, nfaces;
        cgsize_t elemid;
        cgnsELEM *elem;
        for (n = 0; n < numobjs; n++) {
            elemid = objlist[n] >> 3;
            facenum = (int)(objlist[n] & 7);
            if (NULL == (elem = GetElement (elemid, &pos))) {
                sprintf (errmsg, "cgnsImportAddReg:region element %ld not found",
                    (long)elemid);
                cgnsImportFatal (errmsg);
            }
            if (elem->elemtype == cgnsELEM_WDG)
                nfaces = 5;
            else if (elem->elemtype == cgnsELEM_HEX)
                nfaces = 6;
            else
                nfaces = elem->elemtype;
            if (facenum < 1 || facenum > nfaces)
                cgnsImportFatal ("cgnsImportAddReg:region face number out of range");
            region_list[region_nobjs++] = objlist[n];
        }
    }

    /* element region */

    else if (region_type == cgnsREG_ELEMS) {
        for (n = 0; n < numobjs; n++) {
            if (NULL == GetElement (objlist[n], &pos)) {
                sprintf (errmsg, "cgnsImportAddReg:region element %ld not found",
                    (long)objlist[n]);
                cgnsImportFatal (errmsg);
            }
            region_list[region_nobjs++] = objlist[n];
        }
    }

    else
        cgnsImportFatal ("cgnsImportAddReg:undefined region type");

    return (region_nobjs);
}

/*---------- cgnsImportEndReg -------------------------------------------
 * end region definition and import region
 *-----------------------------------------------------------------------*/

cgsize_t cgnsImportEndReg (void)
{
    int pos;
    cgsize_t n;
    cgnsREGN *reg;

    if (!region_id || !region_nobjs)
        return (region_id = 0);
    region_id = 0;

    /* create a new region */

    if (NULL == (reg = GetRegion (region_name, &pos))) {
        reg = NewRegion (region_name, pos);
        reg->type = region_type;
    }
    if (0 == reg->nobjs)
        reg->objid = (cgsize_t *) malloc ((size_t)region_nobjs * sizeof(cgsize_t));
    else
        reg->objid = (cgsize_t *) realloc (reg->objid,
            (size_t)(reg->nobjs + region_nobjs) * sizeof(cgsize_t));
    if (NULL == reg->objid)
        cgnsImportFatal (
            "cgnsImportRegion:malloc failed for the region object list");

    for (n = 0; n < region_nobjs; n++)
        reg->objid[n+reg->nobjs] = region_list[n];
    reg->nobjs += region_nobjs;
    return (reg->nobjs);
}

/*---------- cgnsImportRegion -------------------------------------------
 * import a named region
 *-----------------------------------------------------------------------*/

cgsize_t cgnsImportRegion (char *regname, int regtype,
                           cgsize_t numobjs, cgsize_t *objlist)
{
    cgnsImportBegReg (regname, regtype);
    cgnsImportAddReg (numobjs, objlist);
    return (cgnsImportEndReg ());
}

/*---------- cgnsImportRegionList ---------------------------------------
 * return a list of all region names
 *-----------------------------------------------------------------------*/

char **cgnsImportRegionList (void)
{
    int n, len = 0;
    char **namelist, *names;

    for (n = 0; n < num_regions; n++)
        len += ((int)strlen (reglist[n].name) + 1);
    n = num_regions + 1;
    namelist = (char **) malloc (len + n * sizeof(char *));
    if (NULL == namelist)
        cgnsImportFatal (
            "cgnsImportRegionList:malloc failed for region name list");
    names = (char *) (namelist + n);
    for (n = 0; n < num_regions; n++) {
        namelist[n] = names;
        strcpy (names, reglist[n].name);
        names += (strlen (reglist[n].name) + 1);
    }
    namelist[num_regions] = NULL;
    return (namelist);
}

/*---------- cgnsImportGetRegion ----------------------------------------
 * get node ID's for a region
 *-----------------------------------------------------------------------*/

cgsize_t *cgnsImportGetRegion (char *regname)
{
    int pos;
    cgsize_t n, *objlist;
    cgnsREGN *reg;

    if (NULL == regname || !*regname ||
        NULL == (reg = GetRegion (regname, &pos)))
        return (NULL);
    objlist = (cgsize_t *) malloc ((size_t)(reg->nobjs + 2) * sizeof(cgsize_t));
    if (NULL == objlist)
        cgnsImportFatal (
            "cgnsImportGetRegion:malloc failed for region object ID list");
    objlist[0] = reg->type;
    objlist[1] = reg->nobjs;
    for (n = 0; n < reg->nobjs; n++)
        objlist[n+2] = reg->objid[n];
    return (objlist);
}

/*---------- cgnsImportOpen ---------------------------------------------
 * open CGNS file
 *-----------------------------------------------------------------------*/

int cgnsImportOpen (char *filename)
{
    cgnsImportClose ();
    if (cg_open (filename, CG_MODE_MODIFY, &cgnsFile) &&
        cg_open (filename, CG_MODE_WRITE, &cgnsFile))
        cgnsImportFatal ((char *)cg_get_error());
    return cgnsFile;
}

/*---------- cgnsImportBase ---------------------------------------------
 * set CGNS base
 *-----------------------------------------------------------------------*/

int cgnsImportBase (char *basename)
{
    if (cg_base_write (cgnsFile, basename, 3, 3, &cgnsBase))
        cgnsImportFatal ((char *)cg_get_error());
    return cgnsBase;
}

/*---------- cgnsImportZone ---------------------------------------------
 * set CGNS zone
 *-----------------------------------------------------------------------*/

void cgnsImportZone (char *zonename)
{
    int n;

    for (n = 0; n < num_nodes; n++) {
        if (nodemap[n].node != NULL)
            free (nodemap[n].node);
    }
    for (n = 0; n < num_elements; n++) {
        if (elemlist[n].nodeid != NULL)
            free (elemlist[n].nodeid);
    }
    if (num_regions) {
        for (n = 0; n < num_regions; n++) {
            if (reglist[n].objid != NULL)
                free (reglist[n].objid);
        }
        free (reglist);
    }
    if (num_faces) {
        for (n = 0; n < num_faces; n++) {
            if (facelist[n] != NULL) free (facelist[n]);
        }
        free (facelist);
    }
    num_nodes = num_elements = num_faces = 0;
    num_regions = 0;

    strncpy (cgnsZoneName, zonename, 32);
    cgnsZoneName[32] = 0;
}

/*---------- cgnsImportWrite --------------------------------------------
 * write data to the CGNS file
 *-----------------------------------------------------------------------*/

int cgnsImportWrite (void)
{
    int icoord, isect;
    cgsize_t n, nn, nnodes, sizes[3];
    cgsize_t nc, nconn, *conns, pos;
    CGNS_ENUMT(ElementType_t) elemtype = CGNS_ENUMV(ElementTypeNull);
#ifdef DOUBLE_PRECISION
    CGNS_ENUMT(DataType_t) datatype = CGNS_ENUMV(RealDouble);
    double *xyz;
#else
    CGNS_ENUMT(DataType_t) datatype = CGNS_ENUMV(RealSingle);
    float *xyz;
#endif
    cgnsNODE *node;
    cgnsELEM *elem;
    cgnsREGN *regn;

    if (!cgnsFile)
        cgnsImportFatal ("cgnsImportWrite:CGNS file not open");
    if (region_id)
        cgnsImportEndReg ();

    /* count the nodes */

    nnodes = 0;
    for (n = 0; n < num_nodes; n++) {
        if (nodemap[n].nodeid == nodemap[n].node->id &&
            USED_BIT == (nodemap[n].node->flags & USED_BIT))
            nnodes++;
    }
    if (!nnodes) return 0;

    if (!cgnsBase) cgnsImportBase ("Base");
    if (!*cgnsZoneName) strcpy (cgnsZoneName, "Zone");

    sizes[0] = nnodes;
    sizes[1] = num_elements;
    sizes[2] = 0;

    if (cg_zone_write (cgnsFile, cgnsBase, cgnsZoneName,
            sizes, CGNS_ENUMV(Unstructured), &cgnsZone))
        cgnsImportFatal ((char *)cg_get_error());

    /* write the node list */

#ifdef DOUBLE_PRECISION
    xyz = (double *) malloc ((size_t)nnodes * sizeof(double));
#else
    xyz = (float *) malloc ((size_t)nnodes * sizeof(float));
#endif
    if (NULL == xyz)
        cgnsImportFatal ("cgnsImportWrite:malloc failed for nodes");

    for (nn = 0, n = 0; n < num_nodes; n++) {
        if (nodemap[n].nodeid == nodemap[n].node->id &&
            USED_BIT == (nodemap[n].node->flags & USED_BIT))
#ifdef DOUBLE_PRECISION
            xyz[nn++] = (double)nodemap[n].node->x;
#else
            xyz[nn++] = (float)nodemap[n].node->x;
#endif
    }
    if (cg_coord_write (cgnsFile, cgnsBase, cgnsZone, datatype,
            "CoordinateX", (void *)xyz, &icoord))
        cgnsImportFatal ((char *)cg_get_error());

    for (nn = 0, n = 0; n < num_nodes; n++) {
        if (nodemap[n].nodeid == nodemap[n].node->id &&
            USED_BIT == (nodemap[n].node->flags & USED_BIT))
#ifdef DOUBLE_PRECISION
            xyz[nn++] = (double)nodemap[n].node->y;
#else
            xyz[nn++] = (float)nodemap[n].node->y;
#endif
    }
    if (cg_coord_write (cgnsFile, cgnsBase, cgnsZone, datatype,
            "CoordinateY", (void *)xyz, &icoord))
        cgnsImportFatal ((char *)cg_get_error());

    for (nn = 0, n = 0; n < num_nodes; n++) {
        if (nodemap[n].nodeid == nodemap[n].node->id &&
            USED_BIT == (nodemap[n].node->flags & USED_BIT))
#ifdef DOUBLE_PRECISION
            xyz[nn++] = (double)nodemap[n].node->z;
#else
            xyz[nn++] = (float)nodemap[n].node->z;
#endif
    }
    if (cg_coord_write (cgnsFile, cgnsBase, cgnsZone, datatype,
            "CoordinateZ", (void *)xyz, &icoord))
        cgnsImportFatal ((char *)cg_get_error());

    /* write variables */

    if (num_vars) {
        int isol, ifld, nv;
        if (cg_sol_write(cgnsFile, cgnsBase, cgnsZone,
                "NodeVariables", CGNS_ENUMV(Vertex), &isol))
            cgnsImportFatal ((char *)cg_get_error());
        for (nv = 0; nv < num_vars; nv++) {
            for (nn = 0, n = 0; n < num_nodes; n++) {
                if (nodemap[n].nodeid == nodemap[n].node->id &&
                    USED_BIT == (nodemap[n].node->flags & USED_BIT))
#ifdef DOUBLE_PRECISION
                    xyz[nn++] = (double)nodemap[n].node->vars[nv];
#else
                    xyz[nn++] = (float)nodemap[n].node->vars[nv];
#endif
            }
            if (strlen(var_names[nv]) > 32) var_names[nv][32] = 0;
            if (cg_field_write(cgnsFile, cgnsBase, cgnsZone, isol,
                    datatype, var_names[nv], xyz, &ifld))
                cgnsImportFatal ((char *)cg_get_error());
        }
    }

    free (xyz);

    /* write the element list */

    switch (elemlist->elemtype) {
        case cgnsELEM_TET:
            elemtype = CGNS_ENUMV(TETRA_4);
            break;
        case cgnsELEM_PYR:
            elemtype = CGNS_ENUMV(PYRA_5);
            break;
        case cgnsELEM_WDG:
            elemtype = CGNS_ENUMV(PENTA_6);
            break;
        case cgnsELEM_HEX:
            elemtype = CGNS_ENUMV(HEXA_8);
            break;
    }
    for (n = 0, elem = elemlist; n < num_elements; n++, elem++) {
        if (elem->elemtype != elemlist->elemtype) {
            elemtype = CGNS_ENUMV(MIXED);
            break;
        }
    }

    if (elemtype == CGNS_ENUMV(MIXED)) {
        nconn = 0;
        for (n = 0, elem = elemlist; n < num_elements; n++, elem++)
            nconn += (1 + elem->elemtype);
    }
    else
        nconn = num_elements * elemlist->elemtype;
    conns = (cgsize_t *) malloc ((size_t)nconn * sizeof(cgsize_t));
    if (NULL == conns)
        cgnsImportFatal ("cgnsImportWrite:malloc failed for element data");

    nc = 0;
    for (n = 0, elem = elemlist; n < num_elements; n++, elem++) {
        if (elemtype == CGNS_ENUMV(MIXED)) {
            switch (elem->elemtype) {
                case cgnsELEM_TET :
                    conns[nc] = CGNS_ENUMV(TETRA_4);
                    break;
                case cgnsELEM_PYR:
                    conns[nc] = CGNS_ENUMV(PYRA_5);
                    break;
                case cgnsELEM_WDG:
                    conns[nc] = CGNS_ENUMV(PENTA_6);
                    break;
                case cgnsELEM_HEX:
                    conns[nc] = CGNS_ENUMV(HEXA_8);
                    break;
            }
            nc++;
        }
        for (nn = 0; nn < elem->elemtype; nn++) {
            if (NULL == (node = GetNode (elem->nodeid[nn], &pos)))
                cgnsImportFatal ("cgnsImportWrite:missing element node");
            conns[nc++] = pos + 1;
        }
    }

    if (cg_section_write (cgnsFile, cgnsBase, cgnsZone, "GridElements",
            elemtype, 1, num_elements, 0, conns, &isect))
        cgnsImportFatal ((char *)cg_get_error());

    free (conns);

    /* write the regions */

    nn = num_elements + 1;
    for (n = 0, regn = reglist; n < num_regions; n++, regn++) {
        if (regn->type == cgnsREG_NODES)
            nn += write_node_region (regn, nn);
        else if (regn->type == cgnsREG_FACES)
            nn += write_face_region (regn, nn);
        else if (regn->type == cgnsREG_ELEMS)
            nn += write_elem_region (regn, nn);
    }

    return cgnsZone;
}

/*---------- cgnsImportClose --------------------------------------------
 * close the CGNS file
 *-----------------------------------------------------------------------*/

void cgnsImportClose (void)
{
    if (cgnsFile) {
        cg_close (cgnsFile);
        cgnsFile = 0;
    }
}

