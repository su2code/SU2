/* include file for CGNS imports */

#ifndef _cgnsIMPORT_H_
#define _cgnsIMPORT_H_

#include "cgnslib.h"

#ifndef CGNSTYPES_H
# define cgsize_t int
#endif
#ifndef CGNS_ENUMT
# define CGNS_ENUMT(e) e
# define CGNS_ENUMV(e) e
#endif

/*--- allowable element types ---*/

#define cgnsELEM_TET     4   /* tet element (4 nodes ) */
#define cgnsELEM_PYR     5   /* pyramid element (5 nodes ) */
#define cgnsELEM_WDG     6   /* wedge element (6 nodes) */
#define cgnsELEM_HEX     8   /* hex element (8 nodes) */

/*--- allowable region types ---*/

#define cgnsREG_NODES    1
#define cgnsREG_FACES    2
#define cgnsREG_ELEMS    3

/*--- convert element ID/face number to face ID */

#define cgnsFACEID(elemid,facenum) ((elemid << 3) | (facenum & 7))

/*--- function prototypes ---*/

#ifdef __cplusplus
extern "C" {
#endif

void cgnsImportError (   /* define callback for errors */
    void (*callback)(    /* user-supplied call back routine */
        char *errmsg     /* error message */
    )
);

void cgnsImportFatal (   /* terminate with error message */
    char *errmsg         /* error message */
);

double cgnsImportSetTol (/* setup node checking tolerance */
    double tol           /* tolerance for comparisons */
);

double cgnsImportGetTol (/* get duplicate node tolerance */
    int rel_tol          /* 0 - absolute tol, else relative tol */
);

void cgnsImportSetCheck (/* set dup checking for new nodes */
    int set              /* 0 - allow dup checking, else disallow */
);

cgsize_t cgnsImportRange (/* returns bounding box of nodes */
    double *xmin,        /* lower x limit */
    double *ymin,        /* lower y limit */
    double *zmin,        /* lower z limit */
    double *xmax,        /* upper x limit */
    double *ymax,        /* upper y limit */
    double *zmax         /* upper z limit */
);

cgsize_t cgnsImportCheck (/* check for duplicate nodes */
    int rel_tol          /* 0 - absolute tol, else relative tol */
);

cgsize_t cgnsImportMap ( /* explicitly map 2 nodes */
    cgsize_t nodeid,     /* reference node id */
    cgsize_t mapid       /* node id to map */
);

cgsize_t cgnsImportNode (/* import a node */
    cgsize_t nodeid,     /* node number */
    double x,            /* x coordinate */
    double y,            /* y coordinate */
    double z             /* z coordinate */
);

cgsize_t cgnsImportSetNode (/* set node coordinates */
    cgsize_t nodeid,     /* node ID */
    double x,            /* coordinates */
    double y,
    double z
);

cgsize_t cgnsImportGetNode (/* get node coordinates */
    cgsize_t nodeid,     /* node ID */
    double *x,           /* returned coordinates */
    double *y,
    double *z
);

cgsize_t *cgnsImportNodeList (/* return list of all node ID's */
    void
);

int cgnsImportAddVariable (/*add a variable to the nodes */
    char *varname        /* name of the variable */
);

int cgnsImportGetVariable (/* return variable number */
    char *varname        /* name of the variable */
);

cgsize_t cgnsImportVariable (/* set variable value at node */
    cgsize_t nodeid,     /* node number */
    int varnum,          /* variable number */
    double val           /* variable value */
);

int cgnsImportElement (  /* import an element */
    cgsize_t elemid,     /* element number */
    int elemtype,        /* element type - tet,pyr,wdg or hex */
    cgsize_t *nodelist   /* node numbers defining element */
);

int cgnsImportGetElement (/* get element nodes */
    cgsize_t elemid,     /* element ID */
    cgsize_t nodeid[]    /* returned node IDs */
);

cgsize_t *cgnsImportElementList (/* return list of all element ID's */
    void
);

int cgnsImportGetFace (  /* get element face nodes */
    cgsize_t elemid,     /* element ID */
    int facenum,         /* element face number */
    cgsize_t nodeid[]    /* face nodes */
);

int cgnsImportFindFace ( /* get element face number */
    cgsize_t elemid,     /* element ID */
    int nnodes,          /* number of nodes */
    cgsize_t nodeid[]    /* face nodes */
);

cgsize_t cgnsImportBegReg (/* begin a region specification */
    char *regname,       /* region name */
    int regtype          /* type of region (nodes,faces or elements) */
);

cgsize_t cgnsImportAddReg (/* add nodes to a region specification */
    cgsize_t numobjs,    /* number of objects to add */
    cgsize_t *objlist    /* object list for region */
);

cgsize_t cgnsImportEndReg (/* end region specification */
    void
);

cgsize_t cgnsImportRegion (/* import region of nodes */
    char *regname,       /* region name */
    int regtype,         /* region type */
    cgsize_t numobjs,    /* number of objects in region */
    cgsize_t *objlist    /* object IDs in region */
);

char **cgnsImportRegionList (/* get list of region names */
    void
);

cgsize_t *cgnsImportGetRegion (/* get region object ID's */
    char *regname        /* region name */
);

int cgnsImportOpen (     /* open CGNS file */
    char *filename       /* name of the file */
);

int cgnsImportBase (     /* set CGNS base */
    char *basename       /* name for base */
);

void cgnsImportZone (    /* set CGNS zone */
    char *zonename       /* name for zone */
);

int cgnsImportWrite (    /* write data to CGNS file */
    void
);

void cgnsImportClose (   /* close the CGNS file */
    void
);

#ifdef __cplusplus
}
#endif

#endif  /* _cgnsIMPORT_H_ */
