/*
 * cgnsutil.c - CGNS utility routines
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
# include <io.h>
# include <direct.h>
# define access _access
# define unlink _unlink
# define getcwd _getcwd
#else
# include <unistd.h>
#endif

#include "cgnslib.h"
#include "cgnsutil.h"

#define MAXDIRLEN 256

#ifdef _WIN32
# define PATH_DELIM ';'
# include <direct.h>
# include <io.h>
#else
# define PATH_DELIM ':'
# include <unistd.h>
#endif

#ifndef CG_MODE_READ
# define CG_MODE_READ   MODE_READ
# define CG_MODE_WRITE  MODE_WRITE
# define CG_MODE_MODIFY MODE_MODIFY
#endif

int nZones = 0;
ZONE *Zones;

int baseclass = 0;
int baseunits[5] = {0, 0, 0, 0, 0};
int cgnsbase = 1;
int cgnsfn = 0;

int element_node_counts[] = {
    0, 0,      /* ElementTypeNull, ElementTypeUserDefined */
    1, 2, 3,   /* NODE, BAR_2, BAR_3 */
    3, 6,      /* TRI_3, TRI_6 */
    4, 8, 9,   /* QUAD_4, QUAD_8, QUAD_9 */
    4, 10,     /* TETRA_4, TETRA_10 */
    5, 14,     /* PYRA_5, PYRA_14 */
    6, 15, 18, /* PENTA_6, PENTA_15, PENTA_18 */
    8, 20, 27, /* HEXA_8, HEXA_20, HEXA_27 */
    0, 0       /* MIXED, NGON_n */
};

static char cgnstemp[16] = "";

/*---------- FATAL ----------------------------------------------------
 * exit with error message
 *---------------------------------------------------------------------*/

void FATAL (char *procname, char *errmsg)
{
    char *msg = errmsg;

    if (NULL == msg) {
        msg = (char *)cg_get_error ();
        if (NULL == msg || !*msg)
            msg = "unknown error";
    }
    fflush (stdout);
    if (NULL == procname || !*procname)
        fprintf (stderr, "%s\n", msg);
    else
        fprintf (stderr, "%s:%s\n", procname, msg);
    if (cgnsfn) cg_close (cgnsfn);
    if (cgnstemp[0]) unlink (cgnstemp);
    exit (1);
}

/*---------- new_zone -------------------------------------------------
 * create new zone(s)
 *---------------------------------------------------------------------*/

ZONE *new_zone (int count)
{
    int n;
    ZONE *z;

    z = (ZONE *) calloc (count, sizeof(ZONE));
    if (NULL == z)
        FATAL ("new_zone", "calloc failed for new zones");
    for (n = 0; n < count; n++) {
        z[n].id = n + 1;
        sprintf (z[n].name, "Zone%d", n + 1);
        z[n].type = CGNS_ENUMV(Structured);
        z[n].vertflags = 0;
        z[n].datatype = CGNS_ENUMV(RealDouble);
    }
    return z;
}

/*---------- new_vertex -----------------------------------------------
 * create coordinate array for a zone
 *---------------------------------------------------------------------*/

VERTEX *new_vertex (cgsize_t nverts)
{
    cgsize_t n;
    VERTEX *verts;

    verts = (VERTEX *) calloc ((size_t)nverts, sizeof(VERTEX));
    if (NULL == verts)
        FATAL ("new_vertex", "calloc failed for new vertex array");
    for (n = 0; n < nverts; n++) {
        verts[n].id = n + 1;
        verts[n].w = 1.0;
    }
    return verts;
}

/*---------- new_elemset ----------------------------------------------
 * create element set array
 *---------------------------------------------------------------------*/

ELEMSET *new_elemset (int nesets)
{
    int n;
    ELEMSET *esets;

    esets = (ELEMSET *) calloc (nesets, sizeof(ELEMSET));
    if (NULL == esets)
        FATAL ("new_elemset", "calloc failed for new element set array");
    for (n = 0; n < nesets; n++) {
        esets[n].id = n + 1;
        sprintf (esets[n].name, "ElemSet%d", n + 1);
    }
    return esets;
}

/*---------- new_interface --------------------------------------------
 * create grid 1to1 interface array for a zone
 *---------------------------------------------------------------------*/

INTERFACE *new_interface (int nints)
{
    int n;
    INTERFACE *ints;

    ints = (INTERFACE *) calloc (nints, sizeof(INTERFACE));
    if (NULL == ints)
        FATAL ("new_interface", "calloc failed for new interface array");
    for (n = 0; n < nints; n++) {
        ints[n].id = n + 1;
        sprintf (ints[n].name, "Interface%d", n + 1);
    }
    return ints;
}

/*---------- new_connect ----------------------------------------------
 * create grid connectivity array for a zone
 *---------------------------------------------------------------------*/

CONNECT *new_connect (int nconns)
{
    int n;
    CONNECT *conns;

    conns = (CONNECT *) calloc (nconns, sizeof(CONNECT));
    if (NULL == conns)
        FATAL ("new_connect", "calloc failed for new connectivity array");
    for (n = 0; n < nconns; n++) {
        conns[n].id = n + 1;
        sprintf (conns[n].name, "Connectivity%d", n + 1);
    }
    return conns;
}

/*---------- new_boco -------------------------------------------------
 * create boundary condition array
 *---------------------------------------------------------------------*/

BOCO *new_boco (int nbocos)
{
    int n;
    BOCO *bocos;

    bocos = (BOCO *) calloc (nbocos, sizeof(BOCO));
    if (NULL == bocos)
        FATAL ("new_boco", "calloc failed for new boco array");
    for (n = 0; n < nbocos; n++) {
        bocos[n].id = n + 1;
        sprintf (bocos[n].name, "Boco%d", n + 1);
    }
    return bocos;
}

/*---------- new_solution ---------------------------------------------
 * create solution array for a zone
 *---------------------------------------------------------------------*/

SOLUTION *new_solution (int nsols)
{
    int n;
    SOLUTION *sols;

    sols = (SOLUTION *) calloc (nsols, sizeof(SOLUTION));
    if (NULL == sols)
        FATAL ("new_solution", "calloc failed for new solution array");
    for (n = 0; n < nsols; n++) {
        sols[n].id = n + 1;
        sprintf (sols[n].name, "FlowSolution%d", n + 1);
    }
    return sols;
}

/*---------- new_field ------------------------------------------------
 * create solution variable array for a zone
 *---------------------------------------------------------------------*/

FIELD *new_field (int nflds, cgsize_t size)
{
    int n;
    FIELD *flds;

    flds = (FIELD *) calloc (nflds, sizeof(FIELD));
    if (NULL == flds)
        FATAL ("new_field", "calloc failed for new field array");
    for (n = 0; n < nflds; n++) {
        flds[n].id = n + 1;
        sprintf (flds[n].name, "Field%d", n+1);
        flds[n].datatype = CGNS_ENUMV(RealDouble);
        if (size > 0) {
            flds[n].data = (double *) calloc ((size_t)size, sizeof(double));
            if (NULL == flds[n].data)
                FATAL ("new_field", "calloc failed for field data array");
        }
    }
    return flds;
}

/*---------- new_desc -------------------------------------------------
 * create descriptor array
 *---------------------------------------------------------------------*/

DESC *new_desc (int ndesc)
{
    int n;
    DESC *desc;

    desc = (DESC *) calloc (ndesc, sizeof(DESC));
    if (NULL == desc)
        FATAL ("new_desc", "calloc failed for new descriptor array");
    for (n = 0; n < ndesc; n++) {
        desc[n].id = n + 1;
        sprintf (desc[n].name, "Descriptor%d", n + 1);
    }
    return desc;
}

/*---------- vertex_index ---------------------------------------------
 * get index in vertex array for structured grid
 *---------------------------------------------------------------------*/

cgsize_t vertex_index (ZONE *z, cgsize_t i, cgsize_t j, cgsize_t k)
{
    return i - 1 + z->dim[0] * ((j - 1) + z->dim[1] * (k - 1));
}

/*---------- cell_index -----------------------------------------------
 * get index in cell array for structured grid
 *---------------------------------------------------------------------*/

cgsize_t cell_index (ZONE *z, cgsize_t i, cgsize_t j, cgsize_t k)
{
    return i - 1 + (z->dim[0] - 1) * ((j - 1) + (z->dim[1] - 1) * (k - 1));
}

/*---------- solution_index -------------------------------------------
 * get index in solution array for structured grid
 *---------------------------------------------------------------------*/

cgsize_t solution_index (ZONE *z, SOLUTION *s, cgsize_t i, cgsize_t j, cgsize_t k)
{
    cgsize_t ni, nj;

    ni = z->dim[0] - 1 + s->rind[0][0] + s->rind[0][1];
    nj = z->dim[1] - 1 + s->rind[1][0] + s->rind[1][1];
    i += s->rind[0][0];
    j += s->rind[1][0];
    k += s->rind[2][0];
    return i - 1 + ni * ((j - 1) + nj * (k - 1));
}

/*---------- file_exists ----------------------------------------------
 * check if a file exists
 *---------------------------------------------------------------------*/

int file_exists (char *file)
{
    struct stat st;

    if (access (file, 0) || stat (file, &st) ||
        S_IFREG != (st.st_mode & S_IFMT)) return 0;
    return 1;
}

/*---------- is_executable -----------------------------------------------
 * checks if pathname exists and is executable
 *------------------------------------------------------------------------*/

int is_executable (char *file)
{
    struct stat st;

    /* needs to be executable and not a directory */

    if (access (file, 1) || stat (file, &st) ||
        S_IFDIR == (st.st_mode & S_IFMT))
        return (0);
    return (1);
}

#ifdef _WIN32

/*---------- check_extensions -------------------------------------------
 * check for DOS/Windows executable
 *-----------------------------------------------------------------------*/

static int check_extensions (char *pathname)
{
    int n;
    char *p;
    static char *exts[] = {".com", ".exe", ".bat"};

    /* fix path */

    for (p = pathname; *p; p++) {
        if (*p == '/')
            *p = '\\';
    }
    if (is_executable (pathname))
        return (1);
    for (n = 0; n < 3; n++) {
        strcpy (p, exts[n]);
        if (is_executable (pathname))
            return (1);
    }
    *p = 0;
    return (0);
}

#endif

/*---------- find_executable --------------------------------------------
 * locate and build pathname to executable
 *-----------------------------------------------------------------------*/

char *find_executable (char *exename)
{
    int n;
    char *p, *s;
    static char exepath[MAXDIRLEN+1];

    if (exename == NULL || !*exename)
        return (NULL);

#ifdef _WIN32

    /* full path */

    if (*exename == '/' || *exename == '\\' || *(exename+1) == ':') {
        strcpy (exepath, exename);
        return (check_extensions (exepath) ? exepath : NULL);
    }

    /* get current directory */

    if (NULL == getcwd (exepath, MAXDIRLEN))
        FATAL ("find_executable", "couldn't get working directory");
    p = exepath + strlen(exepath);
    if (*(p-1) == '\\')
        *--p = 0;

    /* relative path */

    if (0 == strncmp (exename, ".\\", 2) ||
        0 == strncmp (exename, "..\\", 3) ||
        0 == strncmp (exename, "./", 2) ||
        0 == strncmp (exename, "../", 3)) {
        if (exename[1] != '.')
            strcpy (p, &exename[1]);
        else {
            if (NULL == (p = strrchr (exepath, '\\')))
                p = exepath;
            strcpy (p, &exename[2]);
        }
        return (check_extensions (exepath) ? exepath : NULL);
    }

    /* current directory */

    *p++ = '\\';
    strcpy (p, exename);
    if (check_extensions (exepath))
        return (exepath);

#else

    /* full path */

    if (*exename == '/') {
        if (is_executable (exename))
            return (strcpy (exepath, exename));
        return (NULL);
    }

    /* relative path */

    if (0 == strncmp (exename, "./", 2) ||
        0 == strncmp (exename, "../", 3)) {
        if (NULL == getcwd (exepath, MAXDIRLEN))
            FATAL ("find_executable", "couldn't get working directory");
        p = exepath + strlen(exepath);
        if (*(p-1) == '/')
            *--p = 0;
        if (exename[1] != '.')
            strcpy (p, &exename[1]);
        else {
            if (NULL == (p = strrchr (exepath, '/')))
                p = exepath;
            strcpy (p, &exename[2]);
        }
        return (is_executable (exepath) ? exepath : NULL);
    }

#endif

    /* scan $PATH environment variable */

    if (NULL == (p = getenv ("PATH")))
        return (NULL);
    while (*p) {
        if (NULL == (s = strchr (p, PATH_DELIM))) {
            strcpy (exepath, p);
            n = (int)strlen (exepath);
        }
        else {
            n = (int)(s++ - p);
            strncpy (exepath, p, n);
            exepath[n] = 0;
        }
        if (n) {
            p = exepath + n;
#ifdef _WIN32
            if (*(p-1) != '\\')
                *p++ = '\\';
            strcpy (p, exename);
            if (check_extensions (exepath))
                return (exepath);
#else
            if (*(p-1) != '/')
                *p++ = '/';
            strcpy (p, exename);
            if (is_executable (exepath))
                return (exepath);
#endif
        }
        if (NULL == (p = s)) break;
    }
    return (NULL);
}

/*---------- find_file ------------------------------------------------
 * get pathname to a file
 *---------------------------------------------------------------------*/

char *find_file (char *filename, char *exename)
{
    char *p;
    static char pathname[MAXDIRLEN+1];

    if (file_exists (filename))
        return strcpy (pathname, filename);
    if ((p = find_executable (exename)) == NULL)
        return NULL;
    strcpy (pathname, p);
    while ((p = strrchr (pathname, '/')) != NULL ||
           (p = strrchr (pathname, '\\')) != NULL) {
        strcpy (p+1, filename);
        if (file_exists (pathname))
            return pathname;
        *p = 0;
    }
    return NULL;
}

/*---------- same_file ------------------------------------------------
 * check if 2 files are the same
 *---------------------------------------------------------------------*/

int same_file (char *file1, char *file2)
{
    int n = file_exists (file1) | (file_exists (file2) << 1);
#ifdef _WIN32
    char path1[257], path2[257];

    if (n == 1 || n == 2) return 0;
    if (_fullpath (path1, file1, 256) != NULL &&
        _fullpath (path2, file2, 256) != NULL)
        return (_stricmp (path1, path2) == 0);
    return (_stricmp (file1, file2) == 0);
#else
    if (n == 3) {
        struct stat st1, st2;
        stat (file1, &st1);
        stat (file2, &st2);
        return (st1.st_ino == st2.st_ino);
    }
    if (n == 0)
        return (strcmp (file1, file2) == 0);
    return 0;
#endif
}

/*---------- temporary_file -------------------------------------------
 * create a temporary file
 *---------------------------------------------------------------------*/

char *temporary_file (char *basename)
{
    char *p, *temp;
    int n;

    if (basename == NULL || !*basename)
        basename = "cgnstmpfile";
    n = (int)strlen (basename);
    temp = (char *) malloc (n + 10);
    if (temp == NULL)
        FATAL ("temporary_file", "malloc failed for temp filename");
    sprintf (temp, "%s.tmp", basename);
    p = temp + strlen(temp);
    for (n = 0; n < 1000; n++) {
        sprintf (p, "%3.3d~", n);
        if (access (temp, 0)) return temp;
    }
    FATAL ("temporary_file", "failed to create temporary filename");
    return 0;
}

/*---------- copy_file ------------------------------------------------
 * make a copy of a file
 *---------------------------------------------------------------------*/

void copy_file (char *oldfile, char *newfile)
{
    int c;
    FILE *oldfp, *newfp;

    if (NULL == (oldfp = fopen (oldfile, "rb")))
        FATAL ("copy_file", "error opening input file for reading");
    if (NULL == (newfp = fopen (newfile, "w+b"))) {
        fclose (oldfp);
        FATAL ("copy_file", "error opening output file for writing");
    }
    while (EOF != (c = getc (oldfp)))
        putc (c, newfp);
    fclose (oldfp);
    fclose (newfp);
}

/*---------- open_cgns ------------------------------------------------
 * open a CGNS file
 *---------------------------------------------------------------------*/

int open_cgns (char *cgnsfile, int read_only)
{
    int nbases;

    if (read_only) {
        if (cg_open (cgnsfile, CG_MODE_READ, &cgnsfn) ||
            cg_nbases (cgnsfn, &nbases))
            FATAL ("open_cgns", NULL);
        return nbases;
    }
    if (cg_open (cgnsfile, CG_MODE_MODIFY, &cgnsfn)) {
        if (cg_open (cgnsfile, CG_MODE_WRITE, &cgnsfn))
            FATAL ("open_cgns", NULL);
        return 0;
    }
    if (cg_nbases (cgnsfn, &nbases))
        FATAL ("open_cgns", NULL);
    return nbases;
}

/*---------- find_base ------------------------------------------------
 * find base id from base name
 *---------------------------------------------------------------------*/

int find_base (char *basename)
{
    int nbases, nb, idum;
    char buff[33];

    if (cg_nbases (cgnsfn, &nbases))
        FATAL ("find_base", NULL);
    for (nb = 1; nb <= nbases; nb++) {
        if (cg_base_read (cgnsfn, nb, buff, &idum, &idum))
            FATAL ("find_base", NULL);
        if (!strcmp (buff, basename)) return nb;
    }
    return 0;
}

/*---------- read_cgns ------------------------------------------------
 * read the CGNS file
 *---------------------------------------------------------------------*/

void read_cgns (void)
{
    int nz, ns;

    read_zones ();
    for (nz = 1; nz <= nZones; nz++) {
        read_zone_grid (nz);
        read_zone_element (nz);
        read_zone_interface (nz);
        read_zone_connect (nz);
        read_zone_boco (nz);
        read_zone_solution (nz);
        for (ns = 1; ns <= Zones[nz-1].nsols; ns++)
            read_solution_field (nz, ns, 0);
    }
}

/*---------- read_zones -----------------------------------------------
 * read zone information from CGNS file
 *---------------------------------------------------------------------*/

int read_zones (void)
{
    int n, nz, nd, celldim;
    cgsize_t sizes[9];
    CGNS_ENUMT(ZoneType_t) zonetype;
    char buff[33];

    if (cg_goto (cgnsfn, cgnsbase, "end"))
        FATAL ("read_zones", NULL);
    read_units (baseunits);
    if (cg_dataclass_read ((CGNS_ENUMT(DataClass_t) *)&baseclass))
        baseclass = 0;

    if (cg_base_read (cgnsfn, cgnsbase, buff, &celldim, &nd) ||
        cg_nzones (cgnsfn, cgnsbase, &nZones))
        FATAL ("read_zones", NULL);
    Zones = new_zone (nZones);

    /* read the zone information */

    for (nz = 0; nz < nZones; nz++) {
        if (cg_zone_read (cgnsfn, cgnsbase, nz+1, buff, sizes) ||
            cg_zone_type (cgnsfn, cgnsbase, nz+1, &zonetype))
            FATAL ("read_zones", NULL);
        if (zonetype != CGNS_ENUMV(Structured) && zonetype != CGNS_ENUMV(Unstructured))
            FATAL ("read_zones", "invalid zone type");
        Zones[nz].id = nz + 1;
        strcpy (Zones[nz].name, buff);
        Zones[nz].type = zonetype;
        Zones[nz].idim = zonetype == CGNS_ENUMV(Structured) ? celldim : 1;
        for (n = 0; n < 3; n++)
            Zones[nz].dim[n] = sizes[n];

        /* get units */

        if (cg_goto (cgnsfn, cgnsbase, "Zone_t", nz+1, "end"))
            FATAL ("read_zones", NULL);
        if (!read_units (Zones[nz].units)) {
            for (n = 0; n < 5; n++)
                Zones[nz].units[n] = baseunits[n];
        }
        if (cg_dataclass_read ((CGNS_ENUMT(DataClass_t) *)&Zones[nz].dataclass))
            Zones[nz].dataclass = baseclass;

        /* get descriptors */

        if (cg_ndescriptors (&nd))
            FATAL ("cg_ndescriptors", NULL);
        if (nd) {
            Zones[nz].ndesc = nd;
            Zones[nz].desc = new_desc (nd);
            for (n = 0; n < nd; n++) {
                if (cg_descriptor_read (n+1, Zones[nz].desc[n].name,
                    &Zones[nz].desc[n].desc))
                    FATAL ("cg_descriptor_read", NULL);
            }
        }

        /* get zone counts */

        if (cg_nsections (cgnsfn, cgnsbase, nz+1, &Zones[nz].nesets) ||
            cg_n1to1 (cgnsfn, cgnsbase, nz+1, &Zones[nz].nints) ||
            cg_nconns (cgnsfn, cgnsbase, nz+1, &Zones[nz].nconns) ||
            cg_nbocos (cgnsfn, cgnsbase, nz+1, &Zones[nz].nbocos) ||
            cg_nsols (cgnsfn, cgnsbase, nz+1, &Zones[nz].nsols))
            FATAL ("read_zones", NULL);
    }
    return nZones;
}

/*---------- read_zone_data -------------------------------------------
 * read all zone data
 *---------------------------------------------------------------------*/

void read_zone_data (int nz)
{
    int ns, nsols;

    read_zone_grid (nz);
    read_zone_element (nz);
    read_zone_interface (nz);
    read_zone_connect (nz);
    read_zone_boco (nz);
    nsols = read_zone_solution (nz);
    for (ns = 1; ns <= nsols; ns++)
        read_solution_field (nz, ns, 0);
}

/*---------- read_zone_grid -------------------------------------------
 * read zone grid coordinates
 *---------------------------------------------------------------------*/

cgsize_t read_zone_grid (int nz)
{
    int n, nc, ncoords;
    int rind[6];
    cgsize_t nn, nverts, rng[2][3];
    CGNS_ENUMT(DataType_t) datatype;
    char buff[33];
    double *xyz;
    ZONE *z = &Zones[nz-1];

    if (cg_goto (cgnsfn, cgnsbase, "Zone_t", nz,
                 "GridCoordinates_t", 1, NULL))
        FATAL ("read_zone_grid", NULL);
    if (cg_rind_read (rind) == CG_OK) {
        for (n = 0; n < 2 * z->idim; n++) {
            if (rind[n])
                FATAL ("read_zone_grid",
                    "currently can't handle coordinate rind");
        }
    }

    if (z->type == CGNS_ENUMV(Structured)) {
        nverts = z->dim[0] * z->dim[1] * z->dim[2];
        for (n = 0; n < 3; n++) {
            rng[0][n] = 1;
            rng[1][n] = z->dim[n];
        }
    }
    else {
        nverts = z->dim[0];
        for (n = 0; n < 3; n++) {
            rng[0][n] = 1;
            rng[1][n] = nverts;
        }
    }

    xyz = (double *) malloc ((size_t)nverts * sizeof(double));
    if (NULL == xyz)
        FATAL ("read_zone_grid", "malloc failed for coordinate working array");
    z->vertflags = 0;
    z->datatype = 0;
    z->nverts = nverts;
    z->verts = new_vertex (nverts);


    /* read the nodes */

    if (cg_ncoords (cgnsfn, cgnsbase, nz, &ncoords))
        FATAL ("read_zone_grid", NULL);
    for (nc = 1; nc <= ncoords; nc++) {
        if (cg_coord_info (cgnsfn, cgnsbase, nz, nc, &datatype, buff) ||
            cg_coord_read (cgnsfn, cgnsbase, nz, buff, CGNS_ENUMV(RealDouble),
                rng[0], rng[1], xyz))
            FATAL ("read_zone_grid", NULL);
        if (z->datatype < datatype) z->datatype = datatype;
        if (0 == strcmp (buff, "CoordinateX")) {
            z->vertflags |= 1;
            for (nn = 0; nn < nverts; nn++)
                z->verts[nn].x = xyz[nn];
        }
        else if (0 == strcmp (buff, "CoordinateY")) {
            z->vertflags |= 2;
            for (nn = 0; nn < nverts; nn++)
                z->verts[nn].y = xyz[nn];
        }
        else if (0 == strcmp (buff, "CoordinateZ")) {
            z->vertflags |= 4;
            for (nn = 0; nn < nverts; nn++)
                z->verts[nn].z = xyz[nn];
        }
        else if (0 == strcmp (buff, "CoordinateR")) {
            z->vertflags |= 8;
            for (nn = 0; nn < nverts; nn++)
                z->verts[nn].z = xyz[nn];
        }
        else if (0 == strcmp (buff, "CoordinateTheta")) {
            z->vertflags |= 16;
            for (nn = 0; nn < nverts; nn++)
                z->verts[nn].z = xyz[nn];
        }
        else if (0 == strcmp (buff, "CoordinatePhi")) {
            z->vertflags |= 32;
            for (nn = 0; nn < nverts; nn++)
                z->verts[nn].z = xyz[nn];
        }
        else {
            FATAL ("read_zone_grid", "unknown coordinate type");
        }
    }
    free (xyz);

    /* cylindrical coordinates */

    if (z->vertflags == 24 || z->vertflags == 28) {
        double r, t;
        for (nn = 0; nn < nverts; nn++) {
            r = z->verts[nn].x;
            t = z->verts[nn].y;
            z->verts[nn].x = r * cos (t);
            z->verts[nn].y = r * sin (t);
        }
        z->vertflags |= 3;
    }
    else if (z->vertflags == 56) {      /* spherical coordinates */
        double r, t, p;
        for (nn = 0; nn < nverts; nn++) {
            r = z->verts[nn].x;
            t = z->verts[nn].y;
            p = z->verts[nn].z;
            z->verts[nn].x = r * sin (t) * cos (p);
            z->verts[nn].y = r * sin (t) * sin (p);
            z->verts[nn].z = r * cos (t);
        }
        z->vertflags |= 7;
    }
    else {
        if (z->vertflags < 1 || z->vertflags > 7)
            FATAL ("read_zone_grid", "unknown coordinate system");
    }
    return nverts;
}

/*---------- read_zone_element ----------------------------------------
 * read zone element sets and elements
 *---------------------------------------------------------------------*/

int read_zone_element (int nz)
{
    int ns, iparent;
    cgsize_t size;
    int rind[6];
    ZONE *z = &Zones[nz-1];
    ELEMSET *eset;

    if (cg_nsections (cgnsfn, cgnsbase, nz, &z->nesets))
        FATAL ("read_zone_element", NULL);
    if (z->nesets) {
        z->esets = eset = new_elemset (z->nesets);
        for (ns = 1; ns <= z->nesets; ns++, eset++) {
            if (cg_goto (cgnsfn, cgnsbase, "Zone_t", nz,
                         "Elements_t", ns, NULL))
                FATAL ("read_zone_element", NULL);
            if (cg_rind_read (rind) == CG_OK && (rind[0] || rind[1]))
                FATAL ("read_zone_element",
                    "currently can't handle element set rind");
            if (cg_section_read (cgnsfn, cgnsbase, nz, ns,
                    eset->name, (CGNS_ENUMT(ElementType_t) *)&eset->type,
                    &eset->start, &eset->end, &eset->nbndry, &iparent) ||
                cg_ElementDataSize (cgnsfn, cgnsbase, nz, ns, &size))
                FATAL ("read_zone_element", NULL);
            eset->conn = (cgsize_t *) malloc ((size_t)size * sizeof(cgsize_t));
            if (NULL == eset->conn)
                FATAL ("read_zone_element",
                    "malloc failed for element connectivity");
            if (iparent) {
                size = 4 * (eset->end - eset->start + 1);
                eset->parent = (cgsize_t *) malloc ((size_t)size * sizeof(cgsize_t));
                if (NULL == eset->conn)
                    FATAL ("read_zone_element","malloc failed for parent data");
            }
            if (cg_elements_read (cgnsfn, cgnsbase, nz, ns,
                    eset->conn, eset->parent))
                FATAL ("read_zone_element", NULL);
        }
    }
    return z->nesets;
}

/*---------- structured_elements --------------------------------------
 * build elements from structured zone
 *---------------------------------------------------------------------*/

int structured_elements (int nz)
{
    cgsize_t i, j, k, n, nelems;
    ZONE *z = &Zones[nz-1];
    ELEMSET *eset;

    if (z->type != CGNS_ENUMV(Structured)) return 0;
    nelems = (z->dim[0] - 1) * (z->dim[1] - 1) * (z->dim[2] - 1);
    if (nelems == 0) return 0;
    z->nesets = 1;
    z->esets = eset = new_elemset (1);
    strcpy (eset->name, "StructuredGridElements");
    eset->type = CGNS_ENUMV(HEXA_8);
    eset->start = 1;
    eset->end = nelems;
    eset->nbndry = 0;
    eset->conn = (cgsize_t *) malloc ((size_t)(8 * nelems) * sizeof(cgsize_t));
    if (NULL == eset->conn)
        FATAL ("structured_elements", "malloc failed for element connectivity");
    eset->parent = NULL;
    for (n = 0, k = 1; k < z->dim[2]; k++) {
        for (j = 1; j < z->dim[1]; j++) {
            for (i = 1; i < z->dim[0]; i++) {
                eset->conn[n++] = vertex_index (z, i,   j,   k) + 1;
                eset->conn[n++] = vertex_index (z, i+1, j,   k) + 1;
                eset->conn[n++] = vertex_index (z, i+1, j+1, k) + 1;
                eset->conn[n++] = vertex_index (z, i,   j+1, k) + 1;
                eset->conn[n++] = vertex_index (z, i,   j,   k+1) + 1;
                eset->conn[n++] = vertex_index (z, i+1, j,   k+1) + 1;
                eset->conn[n++] = vertex_index (z, i+1, j+1, k+1) + 1;
                eset->conn[n++] = vertex_index (z, i,   j+1, k+1) + 1;
            }
        }
    }
    return 1;
}

/*---------- read_zone_interface --------------------------------------
 * read zone 1 to 1 interfaces
 *---------------------------------------------------------------------*/

int read_zone_interface (int nz)
{
    int i, j, n, ni;
    cgsize_t range[2][3], d_range[2][3];
    ZONE *z = &Zones[nz-1];
    INTERFACE *ints;

    if (cg_n1to1 (cgnsfn, cgnsbase, nz, &z->nints))
        FATAL ("read_zone_interface", NULL);
    if (z->nints) {
        z->ints = ints = new_interface (z->nints);
        for (ni = 1; ni <= z->nints; ni++, ints++) {
            if (cg_1to1_read (cgnsfn, cgnsbase, nz, ni,
                    ints->name, ints->d_name, (cgsize_t *)range,
                    (cgsize_t *)d_range, (int *)ints->transform))
                FATAL ("read_zone_interface", NULL);
            for (j = 0; j < 2; j++) {
                for (i = 0; i < 3; i++) {
                    ints->range[i][j] = range[j][i];
                    ints->d_range[i][j] = d_range[j][i];
                }
            }
            for (n = 0; n < nZones; n++) {
                if (!strcmp (Zones[n].name, ints->d_name)) {
                    ints->d_zone = n + 1;
                    break;
                }
            }
        }
    }
    return z->nints;
}

/*---------- read_zone_connect ----------------------------------------
 * read zone connectivities
 *---------------------------------------------------------------------*/

int read_zone_connect (int nz)
{
    int n, nc;
    cgsize_t npnts;
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(GridConnectivityType_t) type;
    CGNS_ENUMT(PointSetType_t) ptype, d_ptype;
    CGNS_ENUMT(ZoneType_t) d_ztype;
    CGNS_ENUMT(DataType_t) d_datatype;
    ZONE *z = &Zones[nz-1];
    CONNECT *conns;

    if (cg_nconns (cgnsfn, cgnsbase, nz, &z->nconns))
        FATAL ("read_zone_connect", NULL);
    if (z->nconns) {
        z->conns = conns = new_connect (z->nconns);
        for (nc = 1; nc <= z->nconns; nc++, conns++) {
            if (cg_conn_info (cgnsfn, cgnsbase, nz, nc,
                    conns->name, &location, &type, &ptype,
                    &conns->npnts, conns->d_name, &d_ztype,
                    &d_ptype, &d_datatype, &conns->d_npnts))
                FATAL ("read_zone_connect", NULL);
            conns->location = location;
            conns->type = type;
            conns->ptype = ptype;
            conns->d_ztype = d_ztype;
            conns->d_ptype = d_ptype;
            npnts = conns->npnts * z->idim;
            conns->pnts = (cgsize_t *) calloc ((size_t)npnts, sizeof(cgsize_t));
            npnts = conns->d_npnts * z->idim;
            conns->d_pnts = (cgsize_t *) calloc ((size_t)npnts, sizeof(cgsize_t));
            if (NULL == conns->pnts || NULL == conns->d_pnts)
                FATAL ("read_zone_connect",
                    "malloc failed for connectivity point arrays");
            if (cg_conn_read (cgnsfn, cgnsbase, nz, nc,
                    conns->pnts, CGNS_ENUMV(Integer), conns->d_pnts))
                FATAL ("read_zone_connect", NULL);
            for (n = 0; n < nZones; n++) {
                if (!strcmp (Zones[n].name, conns->d_name)) {
                    conns->d_zone = n + 1;
                    break;
                }
            }
        }
    }
    return z->nconns;
}

/*---------- read_zone_boco  ------------------------------------------
 * read zone boundary conditions
 *---------------------------------------------------------------------*/

int read_zone_boco (int nz)
{
    int nb, ndatasets;
    cgsize_t npnts;
    CGNS_ENUMT(BCType_t) bctype;
    CGNS_ENUMT(PointSetType_t) ptype;
    CGNS_ENUMT(DataType_t) datatype;
    ZONE *z = &Zones[nz-1];
    BOCO *bocos;

    if (cg_nbocos (cgnsfn, cgnsbase, nz, &z->nbocos))
        FATAL ("read_zone_boco", NULL);
    if (z->nbocos) {
        z->bocos = bocos = new_boco (z->nbocos);
        for (nb = 1; nb <= z->nbocos; nb++, bocos++) {
            if (cg_boco_info (cgnsfn, cgnsbase, nz, nb, bocos->name,
                    &bctype, &ptype, &bocos->npnts, bocos->n_index,
                    &bocos->n_cnt, &datatype, &ndatasets))
                FATAL ("read_zone_boco", NULL);
            bocos->type = bctype;
            bocos->ptype = ptype;
            bocos->n_type = datatype;
            npnts = bocos->npnts * z->idim;
            bocos->pnts = (cgsize_t *) calloc ((size_t)npnts, sizeof(cgsize_t));
            if (NULL == bocos->pnts)
                FATAL ("read_zone_boco",
                    "calloc failed for boco point arrays");
            if (bocos->n_cnt) {
                bocos->n_list = (double *) calloc ((size_t)bocos->n_cnt, sizeof(double));
                if (NULL == bocos->n_list)
                    FATAL ("read_zone_boco",
                        "calloc failed for boco normal list");
            }
            if (cg_boco_read (cgnsfn, cgnsbase, nz, nb,
                    bocos->pnts, bocos->n_list))
                FATAL ("read_zone_boco", NULL);
        }
    }
    return z->nbocos;
}

/*---------- read_zone_solution ---------------------------------------
 * read zone solution
 *---------------------------------------------------------------------*/

int read_zone_solution (int nz)
{
    int i, j, ns, nd;
    CGNS_ENUMT(DataType_t) datatype;
    CGNS_ENUMT(GridLocation_t) location;
    ZONE *z = &Zones[nz-1];
    SOLUTION *sols;

    if (cg_nsols (cgnsfn, cgnsbase, nz, &z->nsols))
        FATAL ("read_zone_solution", NULL);
    if (z->nsols) {
        z->sols = sols = new_solution (z->nsols);
        for (ns = 1; ns <= z->nsols; ns++, sols++) {
            if (cg_sol_info (cgnsfn, cgnsbase, nz, ns,
                    sols->name, &location))
                FATAL ("read_zone_solution", NULL);
            sols->location = location;
            if (z->type == CGNS_ENUMV(Structured)) {
                if (sols->location == CGNS_ENUMV(Vertex)) {
                    for (i = 0; i < 3; i++)
                        for (j = 0; j < 2; j++)
                            sols->rind[i][j] = 0;
                    sols->size = z->dim[0] * z->dim[1] * z->dim[2];
                }
                else if (sols->location == CGNS_ENUMV(CellCenter)) {
                    if (cg_goto (cgnsfn, cgnsbase, "Zone_t", nz,
                          "FlowSolution_t", ns, "end"))
                        FATAL ("read_zone_solution", NULL);
                    if (cg_rind_read ((int *)sols->rind)) {
                        for (i = 0; i < 3; i++)
                            for (j = 0; j < 2; j++)
                                sols->rind[i][j] = 0;
                    }
                    sols->size = 1;
                    for (i = 0; i < 3; i++) {
                        sols->size *= (z->dim[i] - 1 +
                            sols->rind[i][0] + sols->rind[i][1]);
                    }
                }
                else
                    FATAL ("read_zone_solution",
                        "solution location not Vertex or CellCenter");
            }
            else {
                sols->size = sols->location == CGNS_ENUMV(Vertex) ? z->dim[0] : z->dim[1];
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 2; j++)
                        sols->rind[i][j] = 0;
            }
            if (cg_nfields (cgnsfn, cgnsbase, nz, ns, &sols->nflds))
                FATAL ("read_zone_solution", NULL);
            if (sols->nflds) {
                sols->flds = new_field (sols->nflds, 0);
                for (i = 0; i < sols->nflds; i++) {
                    if (cg_field_info (cgnsfn, cgnsbase, nz, ns, i+1,
                            &datatype, sols->flds[i].name))
                        FATAL ("read_zone_solution", NULL);
                }
            }

            /* get units */

            if (cg_goto (cgnsfn, cgnsbase, "Zone_t", nz,
                    "FlowSolution_t", ns, "end"))
                FATAL ("read_zone_solution", NULL);
            if (!read_units (sols->units)) {
                for (i = 0; i < 5; i++)
                    sols->units[i] = z->units[i];
            }
            if (cg_dataclass_read ((CGNS_ENUMT(DataClass_t) *)&sols->dataclass))
                sols->dataclass = z->dataclass;

            /* get descriptors */

            if (cg_ndescriptors (&nd))
                FATAL ("cg_ndescriptors", NULL);
            if (nd) {
                sols->ndesc = nd;
                sols->desc = new_desc (nd);
                for (i = 0; i < nd; i++) {
                    if (cg_descriptor_read (i+1, sols->desc[i].name,
                            &sols->desc[i].desc))
                        FATAL ("cg_descriptor_read", NULL);
                }
            }
        }
    }
    return z->nsols;
}

/*---------- read_solution_field --------------------------------------
 * read solution field data for a zone
 *---------------------------------------------------------------------*/

cgsize_t read_solution_field (int nz, int ns, int nf)
{
    int n, is, ie;
    cgsize_t min[3], max[3];
    CGNS_ENUMT(DataType_t) datatype;
    ZONE *z = &Zones[nz-1];
    SOLUTION *s = &z->sols[ns-1];
    FIELD *f;

    if (z->type == CGNS_ENUMV(Structured)) {
        for (n = 0; n < 3; n++) {
            min[n] = 1;
            max[n] = z->dim[n];
        }
        if (s->location == CGNS_ENUMV(CellCenter)) {
            for (n = 0; n < 3; n++)
                max[n] += s->rind[n][0] + s->rind[n][1] - 1;
        }
    }
    else {
        for (n = 0; n < 3; n++) {
            min[n] = 1;
            max[n] = s->size;
        }
    }
    if (nf) {
        is = ie = nf;
    }
    else {
        is = 1;
        ie = s->nflds;
    }
    f = &s->flds[is-1];
    for (nf = is; nf <= ie; nf++, f++) {
        if (cg_field_info (cgnsfn, cgnsbase, nz, ns, nf,
                &datatype, f->name))
            FATAL ("read_solution_field", NULL);
        f->id = nf;
        f->datatype = datatype;
        f->data = (double *) malloc ((size_t)s->size * sizeof(double));
        if (NULL == f->data)
            FATAL ("read_solution_field",
                "malloc failed for solution field data");
        if (cg_field_read (cgnsfn, cgnsbase, nz, ns, f->name,
                CGNS_ENUMV(RealDouble), min, max, f->data))
            FATAL ("read_solution_field", NULL);

        if (cg_goto (cgnsfn, cgnsbase, "Zone_t", nz,
                "FlowSolution_t", ns, "DataArray_t", nf, "end"))
            FATAL ("read_solution_field", NULL);
        if (!read_units (f->units)) {
            for (n = 0; n < 5; n++)
                f->units[n] = s->units[n];
        }

        /* read data class, conversion and exponents */

        if (cg_dataclass_read ((CGNS_ENUMT(DataClass_t) *)&f->dataclass))
            f->dataclass = s->dataclass;

        if (cg_conversion_info (&datatype))
            f->dataconv[0] = 1.0;
        else {
            f->convtype = datatype;
            if (datatype == CGNS_ENUMV(RealSingle)) {
                float conv[2];
                if (cg_conversion_read (conv))
                    FATAL ("read_solution_field", NULL);
                for (n = 0; n < 2; n++)
                    f->dataconv[n] = conv[n];
            }
            else if (datatype == CGNS_ENUMV(RealDouble)) {
                if (cg_conversion_read (f->dataconv))
                    FATAL ("read_solution_field", NULL);
            }
            else
                FATAL ("cg_conversion_info", "invalid data type");
        }

        if (!cg_exponents_info (&datatype)) {
            f->exptype = datatype;
            if (datatype == CGNS_ENUMV(RealSingle)) {
                float exp[5];
                if (cg_exponents_read (exp))
                    FATAL ("read_solution_field", NULL);
                for (n = 0; n < 5; n++)
                    f->exponent[n] = exp[n];
            }
            else if (datatype == CGNS_ENUMV(RealDouble)) {
                if (cg_exponents_read (f->exponent))
                    FATAL ("read_solution_field", NULL);
            }
            else
                FATAL ("cg_exponents_info", "invalid data type");
        }
    }
    return s->size;
}

/*---------- read_units -----------------------------------------------
 * read unit specifications
 *---------------------------------------------------------------------*/

int read_units (int units[5])
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

/*---------- write_cgns -----------------------------------------------
 * write the CGNS file
 *---------------------------------------------------------------------*/

void write_cgns (void)
{
    int nz, ns;

    write_zones ();
    for (nz = 1; nz <= nZones; nz++) {
        write_zone_grid (nz);
        write_zone_element (nz);
        write_zone_interface (nz);
        write_zone_connect (nz);
        write_zone_boco (nz);
        for (ns = 1; ns <= Zones[nz-1].nsols; ns++) {
            write_zone_solution (nz, ns);
            write_solution_field (nz, ns, 0);
        }
    }
}

/*---------- write_zones ----------------------------------------------
 * write zone information to CGNS file
 *---------------------------------------------------------------------*/

void write_zones (void)
{
    int n, nz;
    cgsize_t sizes[3][3];
    ZONE *z = Zones;

    for (n = 0; n < 5; n++) {
        if (baseunits[n]) {
            if (cg_goto (cgnsfn, cgnsbase, "end") ||
                cg_units_write ((CGNS_ENUMT(MassUnits_t))baseunits[0],
                                (CGNS_ENUMT(LengthUnits_t))baseunits[1],
                                (CGNS_ENUMT(TimeUnits_t))baseunits[2],
                                (CGNS_ENUMT(TemperatureUnits_t))baseunits[3], 
                                (CGNS_ENUMT(AngleUnits_t))baseunits[4]))
                FATAL ("write_zones", NULL);
            break;
        }
    }
    if (baseclass) {
        if (cg_goto (cgnsfn, cgnsbase, "end") ||
            cg_dataclass_write ((CGNS_ENUMT(DataClass_t))baseclass))
            FATAL ("write_zones", NULL);
    }

    /* write the zone information */

    for (nz = 0; nz < nZones; nz++, z++) {
        if (!z->id) continue;
        for (n = 0; n < 3; n++) {
            sizes[0][n] = z->dim[n];
            sizes[1][n] = z->dim[n] - 1;
            sizes[2][n] = 0;
        }
        if (cg_zone_write (cgnsfn, cgnsbase, z->name,
                (cgsize_t *)sizes, (CGNS_ENUMT(ZoneType_t))z->type, &z->id))
            FATAL ("write_zones", NULL);

        for (n = 0; n < 5; n++) {
            if (z->units[n] && z->units[n] != baseunits[n]) {
                if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id, "end") ||
                    cg_units_write ((CGNS_ENUMT(MassUnits_t))z->units[0],
                                    (CGNS_ENUMT(LengthUnits_t))z->units[1],
                                    (CGNS_ENUMT(TimeUnits_t))z->units[2],
                                    (CGNS_ENUMT(TemperatureUnits_t))z->units[3],
                                    (CGNS_ENUMT(AngleUnits_t))z->units[4]))
                    FATAL ("write_zones", NULL);
                break;
            }
        }

        if (z->dataclass && z->dataclass != baseclass) {
            if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id, "end") ||
                cg_dataclass_write ((CGNS_ENUMT(DataClass_t))z->dataclass))
            FATAL ("write_zones", NULL);
        }

        if (z->ndesc) {
            if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id, "end"))
                FATAL ("write_zones", NULL);
            for (n = 0; n < z->ndesc; n++) {
                if (cg_descriptor_write (z->desc[n].name, z->desc[n].desc))
                    FATAL ("cg_descriptor_write", NULL);
            }
        }
    }
}

/*---------- write_zone_data ------------------------------------------
 * write all zone data
 *---------------------------------------------------------------------*/

void write_zone_data (int nz)
{
    int n, ns;
    cgsize_t sizes[3][3];
    ZONE *z = &Zones[nz-1];

    /* write the zone information */

    for (n = 0; n < 3; n++) {
        sizes[0][n] = z->dim[n];
        sizes[1][n] = z->dim[n] - 1;
        sizes[2][n] = 0;
    }
    if (cg_zone_write (cgnsfn, cgnsbase, z->name,
            (cgsize_t *)sizes, (CGNS_ENUMT(ZoneType_t))z->type, &z->id))
        FATAL ("write_zone_data", NULL);

    for (n = 0; n < 5; n++) {
        if (z->units[n] && z->units[n] != baseunits[n]) {
            if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id, "end") ||
                cg_units_write ((CGNS_ENUMT(MassUnits_t))z->units[0],
                                (CGNS_ENUMT(LengthUnits_t))z->units[1],
                                (CGNS_ENUMT(TimeUnits_t))z->units[2],
                                (CGNS_ENUMT(TemperatureUnits_t))z->units[3],
                                (CGNS_ENUMT(AngleUnits_t))z->units[4]))
                FATAL ("write_zone_data", NULL);
            break;
        }
    }

    write_zone_grid (nz);
    write_zone_element (nz);
    write_zone_interface (nz);
    write_zone_connect (nz);
    write_zone_boco (nz);
    for (ns = 1; ns <= Zones[nz-1].nsols; ns++) {
        write_zone_solution (nz, ns);
        write_solution_field (nz, ns, 0);
    }
}

/*---------- write_zone_grid ------------------------------------------
 * write zone grid coordinates
 *---------------------------------------------------------------------*/

void write_zone_grid (int nz)
{
    int nc;
    cgsize_t n;
    ZONE *z = &Zones[nz-1];

    if (z->verts == NULL || (z->vertflags & 7) == 0) return;

    if (z->datatype == CGNS_ENUMV(RealSingle)) {
        float *xyz = (float *) malloc ((size_t)z->nverts * sizeof(float));
        if (NULL == xyz)
            FATAL ("write_zone_grid",
                "malloc failed for coordinate working array");
        if ((z->vertflags & 1) == 1) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = (float)z->verts[n].x;
            if (cg_coord_write (cgnsfn, cgnsbase, z->id, CGNS_ENUMV(RealSingle),
                "CoordinateX", xyz, &nc)) FATAL ("write_zone_grid", NULL);
        }
        if ((z->vertflags & 2) == 2) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = (float)z->verts[n].y;
            if (cg_coord_write (cgnsfn, cgnsbase, z->id, CGNS_ENUMV(RealSingle),
                "CoordinateY", xyz, &nc)) FATAL ("write_zone_grid", NULL);
        }
        if ((z->vertflags & 4) == 4) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = (float)z->verts[n].z;
            if (cg_coord_write (cgnsfn, cgnsbase, z->id, CGNS_ENUMV(RealSingle),
                "CoordinateZ", xyz, &nc)) FATAL ("write_zone_grid", NULL);
        }
        free (xyz);
    }
    else {
        double *xyz = (double *) malloc ((size_t)z->nverts * sizeof(double));
        if (NULL == xyz)
            FATAL ("write_zone_grid",
                "malloc failed for coordinate working array");
        if ((z->vertflags & 1) == 1) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = z->verts[n].x;
            if (cg_coord_write (cgnsfn, cgnsbase, z->id, CGNS_ENUMV(RealDouble),
                "CoordinateX", xyz, &nc)) FATAL ("write_zone_grid", NULL);
        }
        if ((z->vertflags & 2) == 2) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = z->verts[n].y;
            if (cg_coord_write (cgnsfn, cgnsbase, z->id, CGNS_ENUMV(RealDouble),
                "CoordinateY", xyz, &nc)) FATAL ("write_zone_grid", NULL);
        }
        if ((z->vertflags & 4) == 4) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = z->verts[n].z;
            if (cg_coord_write (cgnsfn, cgnsbase, z->id, CGNS_ENUMV(RealDouble),
                "CoordinateZ", xyz, &nc)) FATAL ("write_zone_grid", NULL);
        }
        free (xyz);
    }
}

/*---------- write_zone_element ---------------------------------------
 * write zone element sets and elements
 *---------------------------------------------------------------------*/

void write_zone_element (int nz)
{
    int ns;
    ZONE *z = &Zones[nz-1];
    ELEMSET *eset = z->esets;

    if (eset == NULL || z->type == CGNS_ENUMV(Structured)) return;

    for (ns = 1; ns <= z->nesets; ns++, eset++) {
        if (eset->id) {
            if (cg_section_write (cgnsfn, cgnsbase, z->id,
                    eset->name, (CGNS_ENUMT(ElementType_t))eset->type,
                    eset->start, eset->end, eset->nbndry,
                    eset->conn, &eset->id))
                FATAL ("write_zone_element", NULL);
            if (eset->parent != NULL &&
                cg_parent_data_write (cgnsfn, cgnsbase, z->id,
                    eset->id, eset->parent))
                FATAL ("write_zone_element", NULL);
        }
    }
}

/*---------- write_zone_interface -------------------------------------
 * write zone 1 to 1 interfaces
 *---------------------------------------------------------------------*/

void write_zone_interface (int nz)
{
    int i, j, ni;
    cgsize_t range[2][3], d_range[2][3];
    ZONE *z = &Zones[nz-1];
    INTERFACE *ints = z->ints;

    if (ints == NULL) return;

    for (ni = 1; ni <= z->nints; ni++, ints++) {
        if (ints->id) {
            for (j = 0; j < 2; j++) {
                for (i = 0; i < 3; i++) {
                    range[j][i] = ints->range[i][j];
                    d_range[j][i] = ints->d_range[i][j];
                }
            }
            if (cg_1to1_write (cgnsfn, cgnsbase, z->id,
                    ints->name, ints->d_name, (cgsize_t *)range,
                    (cgsize_t *)d_range, ints->transform, &ints->id))
                FATAL ("write_zone_interface", NULL);
        }
    }
}

/*---------- write_zone_connect ---------------------------------------
 * write zone connectivities
 *---------------------------------------------------------------------*/

void write_zone_connect (int nz)
{
    int nc;
    ZONE *z = &Zones[nz-1];
    CONNECT *conns = z->conns;

    if (conns == NULL) return;

    for (nc = 1; nc <= z->nconns; nc++, conns++) {
        if (conns->id &&
            cg_conn_write (cgnsfn, cgnsbase, z->id,
                conns->name, (CGNS_ENUMT(GridLocation_t))conns->location,
                (CGNS_ENUMT(GridConnectivityType_t))conns->type,
                (CGNS_ENUMT(PointSetType_t))conns->ptype, conns->npnts,
                conns->pnts, conns->d_name,
                (CGNS_ENUMT(ZoneType_t))conns->d_ztype,
                (CGNS_ENUMT(PointSetType_t))conns->d_ptype,
                CGNS_ENUMV(Integer), conns->d_npnts, conns->d_pnts,
                &conns->id))
            FATAL ("write_zone_connect", NULL);
    }
}

/*---------- write_zone_bocos -----------------------------------------
 * write zone boundary conditions
 *---------------------------------------------------------------------*/

void write_zone_boco (int nz)
{
    int nb;
    ZONE *z = &Zones[nz-1];
    BOCO *bocos = z->bocos;

    if (bocos == NULL) return;

    for (nb = 1; nb <= z->nbocos; nb++, bocos++) {
        if (bocos->id &&
           (cg_boco_write (cgnsfn, cgnsbase, z->id,
                bocos->name, (CGNS_ENUMT(BCType_t))bocos->type,
                (CGNS_ENUMT(PointSetType_t))bocos->ptype, bocos->npnts,
                bocos->pnts, &bocos->id) ||
            cg_boco_normal_write (cgnsfn, cgnsbase, z->id, bocos->id,
                bocos->n_index, (int)bocos->n_cnt,
                (CGNS_ENUMT(DataType_t))bocos->n_type, bocos->n_list)))
            FATAL ("write_zone_boco", NULL);
    }
}

/*---------- write_zone_solution --------------------------------------
 * write zone solution
 *---------------------------------------------------------------------*/

void write_zone_solution (int nz, int ns)
{
    int n;
    ZONE *z = &Zones[nz-1];
    SOLUTION *s;

    if (z->sols == NULL || ns < 1 || ns > z->nsols) return;
    s = &z->sols[ns-1];

    if (cg_sol_write (cgnsfn, cgnsbase, z->id, s->name,
        (CGNS_ENUMT(GridLocation_t))s->location, &s->id))
        FATAL ("write_zone_solution", NULL);
    if (z->type == CGNS_ENUMV(Structured) && s->location == CGNS_ENUMV(CellCenter)) {
        if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id,
            "FlowSolution_t", s->id, "end") ||
            cg_rind_write ((int *)s->rind))
            FATAL ("write_zone_solution", NULL);
    }

    for (n = 0; n < 5; n++) {
        if (s->units[n] && s->units[n] != z->units[n]) {
            if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "end") ||
                cg_units_write ((CGNS_ENUMT(MassUnits_t))s->units[0],
                                (CGNS_ENUMT(LengthUnits_t))s->units[1],
                                (CGNS_ENUMT(TimeUnits_t))s->units[2],
                                (CGNS_ENUMT(TemperatureUnits_t))s->units[3],
                                (CGNS_ENUMT(AngleUnits_t))s->units[4]))
                FATAL ("write_zone_solution", NULL);
            break;
        }
    }

    if (s->dataclass && s->dataclass != z->dataclass) {
        if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id,
            "FlowSolution_t", s->id, "end") ||
            cg_dataclass_write ((CGNS_ENUMT(DataClass_t))s->dataclass))
        FATAL ("write_zone_solution", NULL);
    }

    if (s->ndesc) {
        if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id,
            "FlowSolution_t", s->id, "end"))
            FATAL ("write_zone_solution", NULL);
        for (n = 0; n < s->ndesc; n++) {
            if (cg_descriptor_write (s->desc[n].name, s->desc[n].desc))
                FATAL ("cg_descriptor_write", NULL);
        }
    }
}

/*---------- write_solution_field -------------------------------------
 * write solution field data for a zone
 *---------------------------------------------------------------------*/

void write_solution_field (int nz, int ns, int nf)
{
    int n, is, ie;
    float *data = NULL;
    ZONE *z = &Zones[nz-1];
    SOLUTION *s = &z->sols[ns-1];
    FIELD *f;

    if (nf) {
        is = ie = nf;
    }
    else {
        is = 1;
        ie = s->nflds;
    }
    f = &s->flds[is-1];

    for (nf = is; nf <= ie; nf++, f++) {
        if (f->data == NULL) continue;
        if (f->datatype == CGNS_ENUMV(RealSingle)) {
            if (data == NULL) {
                data = (float *) malloc ((size_t)s->size * sizeof(float));
                if (NULL == data)
                    FATAL ("write_solution_field",
                        "malloc failed for working array");
            }
            for (n = 0; n < s->size; n++)
                data[n] = (float)f->data[n];
            if (cg_field_write (cgnsfn, cgnsbase, z->id, s->id,
                CGNS_ENUMV(RealSingle), f->name, data, &f->id))
                FATAL ("write_solution_field", NULL);
        }
        else {
            if (cg_field_write (cgnsfn, cgnsbase, z->id, s->id,
                CGNS_ENUMV(RealDouble), f->name, f->data, &f->id))
                FATAL ("write_solution_field", NULL);
        }
        for (n = 0; n < 5; n++) {
            if (f->units[n] && f->units[n] != s->units[n]) {
                if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                    cg_units_write ((CGNS_ENUMT(MassUnits_t))f->units[0],
                                    (CGNS_ENUMT(LengthUnits_t))f->units[1],
                                    (CGNS_ENUMT(TimeUnits_t))f->units[2],
                                    (CGNS_ENUMT(TemperatureUnits_t))f->units[3],
                                    (CGNS_ENUMT(AngleUnits_t))f->units[4]))
                    FATAL ("write_solution_field", NULL);
                break;
            }
        }
        if (f->dataclass && f->dataclass != s->dataclass) {
            if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                cg_dataclass_write ((CGNS_ENUMT(DataClass_t))f->dataclass))
            FATAL ("write_solution_field", NULL);
        }
        if (f->convtype == CGNS_ENUMV(RealSingle)) {
            float conv[2];
            for (n = 0; n < 2; n++)
                conv[n] = (float)f->dataconv[n];
            if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                cg_conversion_write (CGNS_ENUMV(RealSingle), conv))
                FATAL ("write_solution_field", NULL);
        }
        else if (f->convtype == CGNS_ENUMV(RealDouble)) {
            if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                cg_conversion_write (CGNS_ENUMV(RealDouble), f->dataconv))
                FATAL ("write_solution_field", NULL);
        }
        else {}
        if (f->exptype == CGNS_ENUMV(RealSingle)) {
            float exp[5];
            for (n = 0; n < 5; n++)
                exp[n] = (float)f->dataconv[n];
            if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                cg_exponents_write (CGNS_ENUMV(RealSingle), exp))
                FATAL ("write_solution_field", NULL);
        }
        else if (f->exptype == CGNS_ENUMV(RealDouble)) {
            if (cg_goto (cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                cg_exponents_write (CGNS_ENUMV(RealDouble), f->exponent))
                FATAL ("write_solution_field", NULL);
        }
        else {}
    }
    if (data != NULL) free (data);
}

/*---------- volume_tet -----------------------------------------------
 * compute volume of a tetrahedron
 *---------------------------------------------------------------------*/

double volume_tet (VERTEX *v1, VERTEX *v2, VERTEX *v3, VERTEX *v4)
{
    double vol =
        ((v4->x - v1->x) * ((v2->y - v1->y) * (v3->z - v1->z) -
                            (v2->z - v1->z) * (v3->y - v1->y)) +
         (v4->y - v1->y) * ((v2->z - v1->z) * (v3->x - v1->x) -
                            (v2->x - v1->x) * (v3->z - v1->z)) +
         (v4->z - v1->z) * ((v2->x - v1->x) * (v3->y - v1->y) -
                            (v2->y - v1->y) * (v3->x - v1->x))) / 6.0;
    return vol;
}

/*---------- volume_pyr -----------------------------------------------
 * compute volume of a pyramid
 *---------------------------------------------------------------------*/

double volume_pyr (VERTEX *v1, VERTEX *v2, VERTEX *v3, VERTEX *v4, VERTEX *v5)
{
    VERTEX p;
    double vol;

    p.x = 0.25 * (v1->x + v2->x + v3->x + v4->x);
    p.y = 0.25 * (v1->y + v2->y + v3->y + v4->y);
    p.z = 0.25 * (v1->z + v2->z + v3->z + v4->z);

    vol = (v5->x - p.x) * ((v3->y - v1->y) * (v4->z - v2->z) -
                           (v3->z - v1->z) * (v4->y - v2->y)) +
          (v5->y - p.y) * ((v3->z - v1->z) * (v4->x - v2->x) -
                           (v3->x - v1->x) * (v4->z - v2->z)) +
          (v5->z - p.z) * ((v3->x - v1->x) * (v4->y - v2->y) -
                           (v3->y - v1->y) * (v4->x - v2->x));
    return vol;
}

/*---------- volume_wdg -----------------------------------------------
 * compute volume of a wedge (prism)
 *---------------------------------------------------------------------*/

double volume_wdg (VERTEX *v1, VERTEX *v2, VERTEX *v3,
                   VERTEX *v4, VERTEX *v5, VERTEX *v6)
{
    VERTEX p;
    double vol;

    p.x = (v1->x + v2->x + v3->x + v4->x + v5->x + v6->x) / 6.0;
    p.y = (v1->y + v2->y + v3->y + v4->y + v5->y + v6->y) / 6.0;
    p.z = (v1->z + v2->z + v3->z + v4->z + v5->z + v6->z) / 6.0;

    vol = volume_tet (v1, v2, v3, &p) +
          volume_tet (v4, v6, v5, &p) +
          volume_pyr (v1, v4, v5, v2, &p) +
          volume_pyr (v2, v5, v6, v3, &p) +
          volume_pyr (v1, v3, v6, v4, &p);
    return vol;
}

/*---------- volume_hex -----------------------------------------------
 * compute volume of a hexahedron
 *---------------------------------------------------------------------*/

double volume_hex (VERTEX *v1, VERTEX *v2, VERTEX *v3, VERTEX *v4,
                   VERTEX *v5, VERTEX *v6, VERTEX *v7, VERTEX *v8)
{
    VERTEX p;
    double vol;

    p.x = 0.125 * (v1->x + v2->x + v3->x + v4->x +
                   v5->x + v6->x + v7->x + v8->x);
    p.y = 0.125 * (v1->y + v2->y + v3->y + v4->y +
                   v5->y + v6->y + v7->y + v8->y);
    p.z = 0.125 * (v1->z + v2->z + v3->z + v4->z +
                   v5->z + v6->z + v7->z + v8->z);

    vol = volume_pyr (v1, v2, v3, v4, &p) +
          volume_pyr (v1, v5, v6, v2, &p) +
          volume_pyr (v2, v6, v7, v3, &p) +
          volume_pyr (v3, v7, v8, v4, &p) +
          volume_pyr (v4, v8, v5, v1, &p) +
          volume_pyr (v5, v8, v7, v6, &p);
    return vol;
}

/*---------- volume_element -------------------------------------------
 * compute volume of an element
 *---------------------------------------------------------------------*/

double volume_element (int nn, VERTEX *v[])
{
    switch (nn) {
        case 4: return volume_tet (v[0], v[1], v[2], v[3]);
        case 5: return volume_pyr (v[0], v[1], v[2], v[3], v[4]);
        case 6: return volume_wdg (v[0], v[1], v[2], v[3], v[4], v[5]);
        case 8: return volume_hex (v[0], v[1], v[2], v[3],
                                   v[4], v[5], v[6], v[7]);
    }
    return 0.0;
}

/*---------- cell_volume ----------------------------------------------
 * compute cell volume in a zone
 *---------------------------------------------------------------------*/

double cell_volume (ZONE *z, cgsize_t i, cgsize_t j, cgsize_t k)
{
    cgsize_t ni, nj;
    VERTEX *v;
    double vol;

    ni = z->dim[0];
    nj = ni * z->dim[1];
    v = &z->verts[vertex_index (z, i, j, k)];

    vol = volume_hex (v, v + 1, v + ni + 1, v + ni,
        v + nj, v + nj + 1, v + nj + ni + 1, v + nj + ni);
    return fabs (vol);
}

/*---------- vertex_volumes -------------------------------------------
 * compute cell volumes at vertices
 *---------------------------------------------------------------------*/

void vertex_volumes (int nz)
{
    int ns, nn, et, nv, ii;
    cgsize_t i, j, n, ne, ni;
    ZONE *z = &Zones[nz-1];
    VERTEX *v[8];
    double vol;

    if (z->verts == NULL) read_zone_grid (nz);
    if (z->esets == NULL) {
        if (z->type == CGNS_ENUMV(Structured))
            structured_elements (nz);
        else
            read_zone_element (nz);
    }
    for (i = 0; i < z->nverts; i++)
        z->verts[i].w = 0.0;

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
                    nv = 4;
                    break;
                case CGNS_ENUMV(PYRA_5):
                case CGNS_ENUMV(PYRA_14):
                    nv = 5;
                    break;
                case CGNS_ENUMV(PENTA_6):
                case CGNS_ENUMV(PENTA_15):
                case CGNS_ENUMV(PENTA_18):
                    nv = 6;
                    break;
                case CGNS_ENUMV(HEXA_8):
                case CGNS_ENUMV(HEXA_20):
                case CGNS_ENUMV(HEXA_27):
                    nv = 8;
                    break;
                default:
                    nv = 0;
                    break;
            }
            nn = element_node_counts[et];
            if (nv) {
                for (ii = 0; ii < nv; ii++)
                    v[ii] = &z->verts[z->esets[ns].conn[n+ii]-1];
                vol = fabs (volume_element (nv, v)) / (double)nn;
                for (ii = 0; ii < nn; ii++) {
                    ni = z->esets[ns].conn[n+ii] - 1;
                    z->verts[ni].w += vol;
                }
            }
            n += nn;
        }
    }
}

/*---------- to_cell_vertex -------------------------------------------
 * convert cell-centered field to cell-vertex
 *---------------------------------------------------------------------*/

static void to_cell_vertex (ZONE *z, SOLUTION *s, FIELD *f, double *w)
{
    int nes, et, nv, ii;
    cgsize_t n, i, j, k, ni, nj, nk;
    cgsize_t np, nn, ns, nw;
    double *data;
    double *wsum;

    if (z->type == CGNS_ENUMV(Structured))
        np = z->dim[0] * z->dim[1] * z->dim[2];
    else
        np = z->dim[0];
    data = (double *) malloc ((size_t)np * sizeof(double));
    wsum = (double *) malloc ((size_t)np * sizeof(double));
    if (NULL == data || NULL == wsum)
        FATAL ("to_cell_center", "malloc failed for field data");
    for (n = 0; n < np; n++) {
        data[n] = 0.0;
        wsum[n] = 0.0;
    }

    if (z->type == CGNS_ENUMV(Structured)) {
        ni = z->dim[0];
        nj = z->dim[1];
        nk = z->dim[2];
        for (k = 1; k < nk; k++) {
            for (j = 1; j < nj; j++) {
                for (i = 1; i < ni; i++) {
                    ns = solution_index (z, s, i, j, k);
                    nw = cell_index (z, i, j, k);

                    nn = vertex_index (z, i, j, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, j, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i, j+1, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, j+1, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i, j, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, j, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i, j+1, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, j+1, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];
                }
            }
        }

        if (s->rind[0][0]) {
            for (k = 1; k < nk; k++) {
                for (j = 1; j < nj; j++) {
                    ns = solution_index (z, s, 0, j, k);
                    nw = cell_index (z, 1, j, k);

                    nn = vertex_index (z, 1, j, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, 1, j+1, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, 1, j, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, 1, j+1, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];
                }
            }
        }

        if (s->rind[0][1]) {
            for (k = 1; k < nk; k++) {
                for (j = 1; j < nj; j++) {
                    ns = solution_index (z, s, ni, j, k);
                    nw = cell_index (z, ni-1, j, k);

                    nn = vertex_index (z, ni, j, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, ni, j+1, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, ni, j, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, ni, j+1, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];
                }
            }
        }

        if (s->rind[1][0]) {
            for (k = 1; k < nk; k++) {
                for (i = 1; i < ni; i++) {
                    ns = solution_index (z, s, i, 0, k);
                    nw = cell_index (z, i, 1, k);

                    nn = vertex_index (z, i, 1, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, 1, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i, 1, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, 1, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];
                }
            }
        }

        if (s->rind[1][1]) {
            for (k = 1; k < nk; k++) {
                for (i = 1; i < ni; i++) {
                    ns = solution_index (z, s, i, nj, k);
                    nw = cell_index (z, i, nj-1, k);

                    nn = vertex_index (z, i, nj, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, nj, k);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i, nj, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, nj, k+1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];
                }
            }
        }

        if (s->rind[2][0]) {
            for (j = 1; j < nj; j++) {
                for (i = 1; i < ni; i++) {
                    ns = solution_index (z, s, i, j, 0);
                    nw = cell_index (z, i, j, 1);

                    nn = vertex_index (z, i, j, 1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, j, 1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i, j+1, 1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, j+1, 1);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];
                }
            }
        }

        if (s->rind[2][1]) {
            for (j = 1; j < nj; j++) {
                for (i = 1; i < ni; i++) {
                    ns = solution_index (z, s, i, j, nk);
                    nw = cell_index (z, i, j, nk-1);

                    nn = vertex_index (z, i, j, nk);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, j, nk);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i, j+1, nk);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];

                    nn = vertex_index (z, i+1, j+1, nk);
                    data[nn] += (w[nw] * f->data[ns]);
                    wsum[nn] += w[nw];
                }
            }
        }
    }

    else {
        for (nes = 0; nes < z->nesets; nes++) {
            et = z->esets[nes].type;
            if (et < CGNS_ENUMV(TETRA_4) || et > CGNS_ENUMV(MIXED)) continue;
            for (n = 0, j = z->esets[nes].start-1; j < z->esets[nes].end; j++) {
                if (z->esets[nes].type == CGNS_ENUMV(MIXED))
                    et = (int)z->esets[nes].conn[n++];
                nv = element_node_counts[et];
                if (j < s->size && et >= CGNS_ENUMV(TETRA_4) && et <= CGNS_ENUMV(HEXA_27)) {
                    for (ii = 0; ii < nv; ii++) {
                        nn = z->esets[nes].conn[n+ii] - 1;
                        data[nn] += (w[j] * f->data[j]);
                        wsum[nn] += w[j];
                    }
                }
                n += nv;
            }
        }
    }

    for (n = 0; n < np; n++) {
        if (wsum[n] != 0.0)
            data[n] /= wsum[n];
    }

    free (wsum);
    free (f->data);
    f->data = data;
}

/*---------- cell_vertex_solution -------------------------------------
 * convert cell-center solution to cell-vertex
 *---------------------------------------------------------------------*/

void cell_vertex_solution (int nz, int ns, int weighting)
{
    int nf, nes, et, nv, ii, jj;
    cgsize_t n, i, j, k, ni, nj, nk, nc;
    double *w;
    ZONE *z = &Zones[nz-1];
    SOLUTION *s = &z->sols[ns-1];

    if (s->location == CGNS_ENUMV(Vertex) || !s->size || !s->nflds) return;
    if (s->location != CGNS_ENUMV(CellCenter))
        FATAL ("cell_vertex_solution", "solution not cell-centered");
    if (z->type == CGNS_ENUMV(Structured)) {
        if (z->dim[0] < 2 || z->dim[1] < 2 || z->dim[2] < 2)
            FATAL ("cell_vertex_solution",
                "can only handle 3-dimensional structured zones");
        ni = z->dim[0] - 1;
        nj = z->dim[1] - 1;
        nk = z->dim[2] - 1;
        nc = ni * nj * nk;
    }
    else {
        ni = nc = s->size;
        nj = nk = 1;
        if (z->esets == NULL && !read_zone_element (nz))
            FATAL ("cell_vertex_solution", "no element sets defined");
    }
    for (nf = 0; nf < s->nflds; nf++) {
        if (s->flds[nf].data != NULL) break;
    }
    if (nf == s->nflds) return;

    w = (double *) malloc ((size_t)nc * sizeof(double));
    if (NULL == w)
        FATAL ("cell_vertex_solution",
            "malloc failed for weighting function");

    if (weighting) {
        if (z->verts == NULL) read_zone_grid (nz);
        if (z->type == CGNS_ENUMV(Structured)) {
            n = 0;
            for (k = 1; k <= nk; k++) {
                for (j = 1; j <= nj; j++) {
                    for (i = 1; i <= ni; i++) {
                        w[n++] = fabs (cell_volume (z, i, j, k));
                    }
                }
            }
        }
        else {
            VERTEX *v[8];
            for (n = 0; n < nc; n++)
                w[n] = 0.0;
            for (nes = 0; nes < z->nesets; nes++) {
                et = z->esets[nes].type;
                if (et < CGNS_ENUMV(TETRA_4) || et > CGNS_ENUMV(MIXED)) continue;
                for (n = 0, j = z->esets[nes].start-1; j < z->esets[nes].end; j++) {
                    if (z->esets[nes].type == CGNS_ENUMV(MIXED))
                        et = (int)z->esets[nes].conn[n++];
                    switch (et) {
                        case CGNS_ENUMV(TETRA_4):
                        case CGNS_ENUMV(TETRA_10):
                            nv = 4;
                            break;
                        case CGNS_ENUMV(PYRA_5):
                        case CGNS_ENUMV(PYRA_14):
                            nv = 5;
                            break;
                        case CGNS_ENUMV(PENTA_6):
                        case CGNS_ENUMV(PENTA_15):
                        case CGNS_ENUMV(PENTA_18):
                            nv = 6;
                            break;
                        case CGNS_ENUMV(HEXA_8):
                        case CGNS_ENUMV(HEXA_20):
                        case CGNS_ENUMV(HEXA_27):
                            nv = 8;
                            break;
                        default:
                            nv = 0;
                            break;
                    }
                    if (nv && j < nc) {
                        for (ii = 0; ii < nv; ii++)
                            v[ii] = &z->verts[z->esets[nes].conn[n+ii]-1];
                        w[j] = fabs (volume_element (nv, v));
                    }
                    n += element_node_counts[et];
                }
            }
        }
    }
    else {
        for (n = 0; n < nc; n++)
            w[n] = 1.0;
    }

    for (nf = 0; nf < s->nflds; nf++) {
        if (s->flds[nf].data != NULL)
            to_cell_vertex (z, s, &s->flds[nf], w);
    }
    free (w);

    s->location = CGNS_ENUMV(Vertex);
    for (jj = 0; jj < 3; jj++)
        for (ii = 0; ii < 2; ii++)
            s->rind[jj][ii] = 0;
    s->size = z->nverts;
}

/*---------- to_cell_center -------------------------------------------
 * convert cell-vertex field to cell-centered
 *---------------------------------------------------------------------*/

static void to_cell_center (ZONE *z, SOLUTION *s, FIELD *f, double *w)
{
    int nes, et, nv, ii;
    cgsize_t n, i, j, k, ni, nj, nk, np, nn, ns;
    double *data;
    double *wsum;

    if (z->type == CGNS_ENUMV(Structured)) {
        ni = z->dim[0];
        nj = z->dim[1];
        nk = z->dim[2];
        np = (ni - 1 + s->rind[0][0] + s->rind[0][1]) *
             (nj - 1 + s->rind[1][0] + s->rind[1][1]) *
             (nk - 1 + s->rind[2][0] + s->rind[2][1]);
    }
    else
        np = z->dim[1];

    data = (double *) malloc ((size_t)np * sizeof(double));
    wsum = (double *) malloc ((size_t)np * sizeof(double));
    if (NULL == data || NULL == wsum)
        FATAL ("to_cell_center", "malloc failed for field data");
    for (n = 0; n < np; n++) {
        data[n] = 0.0;
        wsum[n] = 0.0;
    }

    if (z->type == CGNS_ENUMV(Structured)) {
        for (k = 1; k < nk; k++) {
            for (j = 1; j < nj; j++) {
                for (i = 1; i < ni; i++) {
                    ns = solution_index (z, s, i, j, k);

                    nn = vertex_index (z, i, j, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, j, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i, j+1, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, j+1, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i, j, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, j, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i, j+1, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, j+1, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];
                }
            }
        }

        if (s->rind[0][0]) {
            for (k = 1; k < nk; k++) {
                for (j = 1; j < nj; j++) {
                    ns = solution_index (z, s, 0, j, k);

                    nn = vertex_index (z, 1, j, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, 1, j+1, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, 1, j, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, 1, j+1, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];
                }
            }
        }

        if (s->rind[0][1]) {
            for (k = 1; k < nk; k++) {
                for (j = 1; j < nj; j++) {
                    ns = solution_index (z, s, ni, j, k);

                    nn = vertex_index (z, ni, j, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, ni, j+1, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, ni, j, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, ni, j+1, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];
                }
            }
        }

        if (s->rind[1][0]) {
            for (k = 1; k < nk; k++) {
                for (i = 1; i < ni; i++) {
                    ns = solution_index (z, s, i, 0, k);

                    nn = vertex_index (z, i, 1, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, 1, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i, 1, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, 1, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];
                }
            }
        }

        if (s->rind[1][1]) {
            for (k = 1; k < nk; k++) {
                for (i = 1; i < ni; i++) {
                    ns = solution_index (z, s, i, nj, k);

                    nn = vertex_index (z, i, nj, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, nj, k);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i, nj, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, nj, k+1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];
                }
            }
        }

        if (s->rind[2][0]) {
            for (j = 1; j < nj; j++) {
                for (i = 1; i < ni; i++) {
                    ns = solution_index (z, s, i, j, 0);

                    nn = vertex_index (z, i, j, 1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, j, 1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i, j+1, 1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, j+1, 1);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];
                }
            }
        }

        if (s->rind[2][1]) {
            for (j = 1; j < nj; j++) {
                for (i = 1; i < ni; i++) {
                    ns = solution_index (z, s, i, j, nk);

                    nn = vertex_index (z, i, j, nk);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, j, nk);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i, j+1, nk);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];

                    nn = vertex_index (z, i+1, j+1, nk);
                    data[ns] += (w[nn] * f->data[nn]);
                    wsum[ns] += w[nn];
                }
            }
        }

        if (s->rind[0][0] && s->rind[1][0]) {
            for (k = 1; k < nk; k++) {
                ns = solution_index (z, s, 0, 0, k);

                nn = solution_index (z, s, 0, 1, k);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, 1, 0, k);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[0][0] && s->rind[1][1]) {
            for (k = 1; k < nk; k++) {
                ns = solution_index (z, s, 0, nj, k);

                nn = solution_index (z, s, 0, nj-1, k);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, 1, nj, k);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[0][1] && s->rind[1][0]) {
            for (k = 1; k < nk; k++) {
                ns = solution_index (z, s, ni, 0, k);

                nn = solution_index (z, s, ni, 1, k);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, ni-1, 0, k);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[0][1] && s->rind[1][1]) {
            for (k = 1; k < nk; k++) {
                ns = solution_index (z, s, ni, nj, k);

                nn = solution_index (z, s, ni, nj-1, k);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, ni-1, nj, k);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[0][0] && s->rind[2][0]) {
            for (j = 1; j < nj; j++) {
                ns = solution_index (z, s, 0, j, 0);

                nn = solution_index (z, s, 0, j, 1);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, 1, j, 0);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[0][0] && s->rind[2][1]) {
            for (j = 1; j < nj; j++) {
                ns = solution_index (z, s, 0, j, nk);

                nn = solution_index (z, s, 0, j, nk-1);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, 1, j, nk);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[0][1] && s->rind[2][0]) {
            for (j = 1; j < nj; j++) {
                ns = solution_index (z, s, ni, j, 0);

                nn = solution_index (z, s, ni, j, 1);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, ni-1, j, 0);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[0][1] && s->rind[2][1]) {
            for (j = 1; j < nj; j++) {
                ns = solution_index (z, s, ni, j, nk);

                nn = solution_index (z, s, ni, j, nk-1);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, ni-1, j, nk);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[1][0] && s->rind[2][0]) {
            for (i = 1; i < ni; i++) {
                ns = solution_index (z, s, i, 0, 0);

                nn = solution_index (z, s, i, 0, 1);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, i, 1, 0);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[1][0] && s->rind[2][1]) {
            for (i = 1; i < ni; i++) {
                ns = solution_index (z, s, i, 0, nk);

                nn = solution_index (z, s, i, 0, nk-1);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, i, 1, nk);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[1][1] && s->rind[2][0]) {
            for (i = 1; i < ni; i++) {
                ns = solution_index (z, s, i, nj, 0);

                nn = solution_index (z, s, i, nj, 1);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, i, nj-1, 0);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[1][1] && s->rind[2][1]) {
            for (i = 1; i < ni; i++) {
                ns = solution_index (z, s, i, nj, nk);

                nn = solution_index (z, s, i, nj, nk-1);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];

                nn = solution_index (z, s, i, nj-1, nk);
                data[ns] += data[nn];
                wsum[ns] += wsum[nn];
            }
        }

        if (s->rind[0][0] && s->rind[1][0] && s->rind[2][0]) {
            ns = solution_index (z, s, 0, 0, 0);

            nn = solution_index (z, s, 1, 0, 0);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, 0, 1, 0);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, 0, 0, 1);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];
        }

        if (s->rind[0][1] && s->rind[1][0] && s->rind[2][0]) {
            ns = solution_index (z, s, ni, 0, 0);

            nn = solution_index (z, s, ni-1, 0, 0);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, ni, 1, 0);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, ni, 0, 1);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];
        }

        if (s->rind[0][0] && s->rind[1][1] && s->rind[2][0]) {
            ns = solution_index (z, s, 0, nj, 0);

            nn = solution_index (z, s, 1, nj, 0);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, 0, nj-1, 0);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, 0, nj, 1);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];
        }

        if (s->rind[0][1] && s->rind[1][1] && s->rind[2][0]) {
            ns = solution_index (z, s, ni, nj, 0);

            nn = solution_index (z, s, ni-1, nj, 0);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, ni, nj-1, 0);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, ni, nj, 1);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];
        }

        if (s->rind[0][0] && s->rind[1][0] && s->rind[2][1]) {
            ns = solution_index (z, s, 0, 0, nk);

            nn = solution_index (z, s, 1, 0, nk);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, 0, 1, nk);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, 0, 0, nk-1);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];
        }

        if (s->rind[0][1] && s->rind[1][0] && s->rind[2][1]) {
            ns = solution_index (z, s, ni, 0, nk);

            nn = solution_index (z, s, ni-1, 0, nk);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, ni, 1, nk);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, ni, 0, nk-1);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];
        }

        if (s->rind[0][0] && s->rind[1][1] && s->rind[2][1]) {
            ns = solution_index (z, s, 0, nj, nk);

            nn = solution_index (z, s, 1, nj, nk);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, 0, nj-1, nk);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, 0, nj, nk-1);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];
        }

        if (s->rind[0][1] && s->rind[1][1] && s->rind[2][1]) {
            ns = solution_index (z, s, ni, nj, nk);

            nn = solution_index (z, s, ni-1, nj, nk);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, ni, nj-1, nk);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];

            nn = solution_index (z, s, ni, nj, nk-1);
            data[ns] += data[nn];
            wsum[ns] += wsum[nn];
        }
    }

    else {
        for (nes = 0; nes < z->nesets; nes++) {
            et = z->esets[nes].type;
            if (et < CGNS_ENUMV(TETRA_4) || et > CGNS_ENUMV(MIXED)) continue;
            for (n = 0, j = z->esets[nes].start-1; j < z->esets[nes].end; j++) {
                if (z->esets[nes].type == CGNS_ENUMV(MIXED))
                    et = (int)z->esets[nes].conn[n++];
                nv = element_node_counts[et];
                if (j < np && et >= CGNS_ENUMV(TETRA_4) && et <= CGNS_ENUMV(HEXA_27)) {
                    for (ii = 0; ii < nv; ii++) {
                        nn = z->esets[nes].conn[n+ii] - 1;
                        data[j] += (w[nn] * f->data[nn]);
                        wsum[j] += w[nn];
                    }
                }
                n += nv;
            }
        }
    }

    for (n = 0; n < np; n++) {
        if (wsum[n] != 0.0)
            data[n] /= wsum[n];
    }

    free (wsum);
    free (f->data);
    f->data = data;
}

/*---------- cell_center_solution -------------------------------------
 * convert cell-vertex solution to cell-center
 *---------------------------------------------------------------------*/

void cell_center_solution (int nz, int ns, int weighting)
{
    int ii, jj, nf;
    cgsize_t n, i, j, k, ni, nj, nk, np, size;
    double *w, vol;
    ZONE *z = &Zones[nz-1];
    SOLUTION *s = &z->sols[ns-1];

    if (s->location == CGNS_ENUMV(CellCenter) || !s->size || !s->nflds) return;
    if (s->location != CGNS_ENUMV(Vertex))
        FATAL ("cell_center_solution", "solution not at cell-vertex");
    if (z->type == CGNS_ENUMV(Structured)) {
        if (z->dim[0] < 2 || z->dim[1] < 2 || z->dim[2] < 2)
            FATAL ("cell_vertex_solution",
                "can only handle 3-dimensional structured zones");
        ni = z->dim[0];
        nj = z->dim[1];
        nk = z->dim[2];
        np = ni * nj * nk;
        for (ii = 0; ii < 3; ii++) {
            for (jj = 0; jj < 2; jj++) {
                if (s->rind[ii][jj] < 0 || s->rind[ii][jj] > 1)
                    FATAL ("cell_center_solution", "rind must be 0 or 1");
            }
        }
        size = (ni - 1 + s->rind[0][0] + s->rind[0][1]) *
               (nj - 1 + s->rind[1][0] + s->rind[1][1]) *
               (nk - 1 + s->rind[2][0] + s->rind[2][1]);
    }
    else {
        ni = np = z->dim[0];
        nj = nk = 1;
        if (z->esets == NULL && !read_zone_element (nz))
            FATAL ("cell_vertex_solution", "no element sets defined");
        size = z->dim[1];
    }
    for (nf = 0; nf < s->nflds; nf++) {
        if (s->flds[nf].data != NULL) break;
    }
    if (nf == s->nflds) return;

    w = (double *) malloc ((size_t)np * sizeof(double));
    if (NULL == w)
        FATAL ("cell_center_solution", "malloc failed for weighting function");

    if (weighting) {
        if (z->verts == NULL) read_zone_grid (nz);
        if (z->type == CGNS_ENUMV(Structured)) {
            int *cnt = (int *) malloc ((size_t)np * sizeof(int));
            if (NULL == cnt)
                FATAL ("cell_center_solution", "malloc failed for cnt array");
            for (n = 0; n < np; n++) {
                w[n] = 0.0;
                cnt[n] = 0;
            }

            for (k = 1; k < nk; k++) {
                for (j = 1; j < nj; j++) {
                    for (i = 1; i < ni; i++) {
                        vol = cell_volume (z, i, j, k);
                        n = vertex_index (z, i, j, k);
                        w[n] += vol;
                        cnt[n]++;
                        n = vertex_index (z, i+1, j, k);
                        w[n] += vol;
                        cnt[n]++;
                        n = vertex_index (z, i, j+1, k);
                        w[n] += vol;
                        cnt[n]++;
                        n = vertex_index (z, i+1, j+1, k);
                        w[n] += vol;
                        cnt[n]++;
                        n = vertex_index (z, i, j, k+1);
                        w[n] += vol;
                        cnt[n]++;
                        n = vertex_index (z, i+1, j, k+1);
                        w[n] += vol;
                        cnt[n]++;
                        n = vertex_index (z, i, j+1, k+1);
                        w[n] += vol;
                        cnt[n]++;
                        n = vertex_index (z, i+1, j+1, k+1);
                        w[n] += vol;
                        cnt[n]++;
                    }
                }
            }
            for (n = 0; n < np; n++)
                w[n] /= (double)cnt[n];
            free (cnt);
        }
        else {
            for (n = 0; n < np; n++)
                w[n] = z->verts[n].w;
            vertex_volumes (nz);
            for (n = 0; n < np; n++) {
                vol = z->verts[n].w;
                z->verts[n].w = w[n];
                w[n] = vol;
            }
        }
    }
    else {
        for (n = 0; n < np; n++)
            w[n] = 1.0;
    }

    for (nf = 0; nf < s->nflds; nf++) {
        if (s->flds[nf].data != NULL)
            to_cell_center (z, s, &s->flds[nf], w);
    }
    free (w);

    s->location = CGNS_ENUMV(CellCenter);
    s->size = size;
}
