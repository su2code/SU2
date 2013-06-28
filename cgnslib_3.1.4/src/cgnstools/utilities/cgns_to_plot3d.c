/*
 * cgns_to_plot3d.c - read CGNS file and write Plot3d file
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
# define unlink _unlink
#else
# include <unistd.h>
#endif

#include "getargs.h"
#include "cgnslib.h"
#include "cgnsutil.h"

#ifdef USE_FORTRAN
#include "pd3fint.h"
#else
extern void OPENF   (int *,char *,int);
extern void CLOSEF  (void);
extern void WRITEIF (int *,int *,int *);
extern void WRITEFF (int *,float *,int *);
extern void WRITEDF (int *,double *,int *);
extern void WRITEGFF(int *,float *,int *,int *);
extern void WRITEGDF(int *,double *,int *,int *);
#endif

/* command line options */

static char options[] = "spnfudb:B:g:wS:";

static char *usgmsg[] = {
    "usage  : cgns_to_plot3d [options] CGNSfile XYZfile [Qfile]",
    "options:",
    "   -s       = write Plot3d as single block",
    "   -p       = write Plot3d as planar format",
    "   -n       = don't write iblank array",
    "   -f       = write formatted (ASCII) Plot3d file",
    "   -u       = write Fortran unformatted Plot3d file",
    "   -d       = use double-precision (64-bit)",
    "   -b<base> = use CGNS base number <base> (default 1)",
    "   -B<name> = use CGNS base named <name>",
    "   -S<sol>  = solution to use if multiple (default 1)",
    "   -w       = use volume weighting",
    "   -g<gamma>= gamma for data conversions (default 1.4)",
    "Default is multi-block binary format with iblank data.",
    NULL
};

static int format = 'b';
static int mblock = 1;
static int whole = 1;
static int use_iblank = 1;
static int weighting = 0;
static int usesol = 1;
static int use_double = 0;

static double reference[4];
static double gamma = 1.4;

static int nblocks = 0;
static int *iblank;
static double *q[5];

/*---------- compute_iblank -------------------------------------------
 * fill in iblank values for a zone
 *---------------------------------------------------------------------*/

static void compute_iblank (int nz)
{
    cgsize_t n, ni, i, j, k, nk, nj, nn;
    cgsize_t ns[3], ne[3];
    ZONE *zone = &Zones[nz];
    INTERFACE *ints = zone->ints;

    for (n = 0; n < zone->nverts; n++)
        iblank[n] = 1;
    nj = zone->dim[0];
    nk = zone->dim[1] * nj;
    for (ni = 0; ni < zone->nints; ni++, ints++) {
        if (ints->d_zone) {
            for (n = 0; n < 3; n++) {
                if (ints->range[n][1] < ints->range[n][0]) {
                    ns[n] = ints->range[n][1] - 1;
                    ne[n] = ints->range[n][0] - 1;
                }
                else {
                    ns[n] = ints->range[n][0] - 1;
                    ne[n] = ints->range[n][1] - 1;
                }
            }
            for (k = ns[2]; k <= ne[2]; k++) {
                for (j = ns[1]; j <= ne[1]; j++) {
                    nn = k * nk + j * nj + ns[0];
                    for (i = ns[0]; i <= ne[0]; i++)
                        iblank[nn++] = -ints->d_zone;
                }
            }
        }
    }
}

/*---------- write_xyz_binary -----------------------------------------
 * write binary Plot3d file
 *---------------------------------------------------------------------*/

static void write_xyz_binary (char *xyzfile)
{
    int nz, i, dims[3];
    int n, k, nk, np;
    VERTEX *verts;
    void *xyz;
    float *xyzf;
    double *xyzd;
    FILE *fp;

    printf ("\nwriting binary XYZ file to %s\n", xyzfile);
    printf ("  in %s-precision", use_double ? "double" : "single");
    if (use_iblank) printf (" with iblank array");
    putchar ('\n');

    for (np = 0, nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            nk = whole ? (int)Zones[nz].nverts :
                (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
            if (np < nk) np = nk;
        }
    }
    if (use_double) {
        xyz = (void *) malloc (np * sizeof(double));
        xyzd = (double *)xyz;
    }
    else {
        xyz = (void *) malloc (np * sizeof(float));
        xyzf = (float *)xyz;
    }
    if (NULL == xyz)
        FATAL ("write_xyx_binary", "malloc failed for working arrays");

    if (NULL == (fp = fopen (xyzfile, "w+b"))) {
        fprintf (stderr, "couldn't open <%s> for writing\n", xyzfile);
        exit (1);
    }
    if (mblock || nblocks > 1)
        fwrite (&nblocks, sizeof(int), 1, fp);
    for (nz = 0; nz < nZones; nz++) {
	if (Zones[nz].type == CGNS_ENUMV(Structured)) {
	    for (i = 0; i < 3; i++)
		dims[i] = (int)Zones[nz].dim[i];
            fwrite (dims, sizeof(int), 3, fp);
	}
    }
    for (nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            printf ("  zone %d - %d x %d x %d ... ", nz+1,
                (int)Zones[nz].dim[0], (int)Zones[nz].dim[1],
		(int)Zones[nz].dim[2]);
            fflush (stdout);
            if (use_iblank) compute_iblank (nz);
            if (whole) {
                np = (int)Zones[nz].nverts;
                nk = 1;
            }
            else {
                np = (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
                nk = (int)Zones[nz].dim[2];
            }
            for (k = 0; k < nk; k++) {
                verts = &Zones[nz].verts[k*np];
                if (use_double) {
                    for (n = 0; n < np; n++)
                        xyzd[n] = verts[n].x;
                    fwrite (xyzd, sizeof(double), np, fp);
                    for (n = 0; n < np; n++)
                        xyzd[n] = verts[n].y;
                    fwrite (xyzd, sizeof(double), np, fp);
                    for (n = 0; n < np; n++)
                        xyzd[n] = verts[n].z;
                    fwrite (xyzd, sizeof(double), np, fp);
                }
                else {
                    for (n = 0; n < np; n++)
                        xyzf[n] = (float)verts[n].x;
                    fwrite (xyzf, sizeof(float), np, fp);
                    for (n = 0; n < np; n++)
                        xyzf[n] = (float)verts[n].y;
                    fwrite (xyzf, sizeof(float), np, fp);
                    for (n = 0; n < np; n++)
                        xyzf[n] = (float)verts[n].z;
                    fwrite (xyzf, sizeof(float), np, fp);
                }
                if (use_iblank)
                    fwrite (&iblank[k*np], sizeof(int), np, fp);
            }
            puts ("done");
        }
    }
    fclose (fp);
    free (xyz);
}

/*---------- write_xyz_unformatted ------------------------------------
 * write unformatted Plot3d file
 *---------------------------------------------------------------------*/

static void write_xyz_unformatted (char *xyzfile)
{
    int n, i, ierr, *indices;
    int nz, k, nk, np;
    char buff[129];
    VERTEX *verts;
    void *xyz;
    float *xyzf;
    double *xyzd;

    printf ("\nwriting unformatted XYZ file to %s\n", xyzfile);
    printf ("  in %s-precision", use_double ? "double" : "single");
    if (use_iblank) printf (" with iblank array");
    putchar ('\n');

    for (np = 0, nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            nk = whole ? (int)Zones[nz].nverts :
                (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
            if (np < nk) np = nk;
        }
    }
    indices = (int *) malloc (3 * nblocks * sizeof(int));
    if (use_double) {
        xyz = (void *) malloc (3 * np * sizeof(double));
        xyzd = (double *)xyz;
    }
    else {
        xyz = (void *) malloc (3 * np * sizeof(float));
        xyzf = (float *)xyz;
    }
    if (NULL == indices || NULL == xyz)
        FATAL ("write_xyx_unformatted", "malloc failed for working arrays");

    unlink (xyzfile);
    strcpy (buff, xyzfile);
    for (n = (int)strlen(buff); n < 128; n++)
        buff[n] = ' ';
    buff[128] = 0;
    n = 0;
    OPENF (&n, buff, 128);

    if (mblock || nblocks > 1) {
        n = 1;
        WRITEIF (&n, &nblocks, &ierr);
    }
    for (n = 0, nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            for (i = 0; i < 3; i++)
                indices[n++] = (int)Zones[nz].dim[i];
        }
    }
    WRITEIF (&n, indices, &ierr);

    for (nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            printf ("  zone %d - %d x %d x %d ... ", nz+1,
                (int)Zones[nz].dim[0], (int)Zones[nz].dim[1],
		(int)Zones[nz].dim[2]);
            fflush (stdout);
            if (use_iblank) compute_iblank (nz);
            if (whole) {
                np = (int)Zones[nz].nverts;
                nk = 1;
            }
            else {
                np = (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
                nk = (int)Zones[nz].dim[2];
            }
            for (k = 0; k < nk; k++) {
                verts = &Zones[nz].verts[k*np];
                if (use_double) {
                    for (i = 0, n = 0; n < np; n++)
                        xyzd[i++] = verts[n].x;
                    for (n = 0; n < np; n++)
                        xyzd[i++] = verts[n].y;
                    for (n = 0; n < np; n++)
                        xyzd[i++] = verts[n].z;
                    if (use_iblank)
                        WRITEGDF (&np, xyzd, &iblank[k*np], &ierr);
                    else
                        WRITEDF (&i, xyzd, &ierr);
                }
                else {
                    for (i = 0, n = 0; n < np; n++)
                        xyzf[i++] = (float)verts[n].x;
                    for (n = 0; n < np; n++)
                        xyzf[i++] = (float)verts[n].y;
                    for (n = 0; n < np; n++)
                        xyzf[i++] = (float)verts[n].z;
                    if (use_iblank)
                        WRITEGFF (&np, xyzf, &iblank[k*np], &ierr);
                    else
                        WRITEFF (&i, xyzf, &ierr);
                }
            }
            puts ("done");
        }
    }
    CLOSEF ();
    free (indices);
    free (xyz);
}

/*---------- write_xyz_formatted --------------------------------------
 * write formatted Plot3d file
 *---------------------------------------------------------------------*/

static void write_xyz_formatted (char *xyzfile)
{
    int n, k, nk, nz, np, *ib;
    VERTEX *verts;
    FILE *fp;

    printf ("\nwriting formatted XYZ file to %s\n", xyzfile);
    if (use_iblank) printf ("  with iblank array\n");

    if (NULL == (fp = fopen (xyzfile, "w+"))) {
        fprintf (stderr, "couldn't open <%s> for writing\n", xyzfile);
        exit (1);
    }
    if (mblock || nblocks > 1)
        fprintf (fp, "%d\n", nblocks);
    for (nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured))
            fprintf (fp, "%d %d %d\n", (int)Zones[nz].dim[0],
                (int)Zones[nz].dim[1], (int)Zones[nz].dim[2]);
    }
    for (nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            printf ("  zone %d - %d x %d x %d ... ", nz+1,
                (int)Zones[nz].dim[0], (int)Zones[nz].dim[1],
		(int)Zones[nz].dim[2]);
            fflush (stdout);
            if (use_iblank) compute_iblank (nz);
            if (whole) {
                np = (int)Zones[nz].nverts;
                nk = 1;
            }
            else {
                np = (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
                nk = (int)Zones[nz].dim[2];
            }
            for (k = 0; k < nk; k++) {
                verts = &Zones[nz].verts[k*np];
                for (n = 0; n < np; n++) {
                    if (n) putc ((n % 5) == 0 ? '\n' : ' ', fp);
                    fprintf (fp, "%#g", verts[n].x);
                }
                putc ('\n', fp);
                for (n = 0; n < np; n++) {
                    if (n) putc ((n % 5) == 0 ? '\n' : ' ', fp);
                    fprintf (fp, "%#g", verts[n].y);
                }
                putc ('\n', fp);
                for (n = 0; n < np; n++) {
                    if (n) putc ((n % 5) == 0 ? '\n' : ' ', fp);
                    fprintf (fp, "%#g", verts[n].z);
                }
                putc ('\n', fp);
                if (use_iblank) {
                    ib = &iblank[k*np];
                    for (n = 0; n < np; n++, ib++) {
                        if (n && (n % 10) == 0)
                            putc ('\n', fp);
                        fprintf (fp, "%5d", *ib);
                    }
                    putc ('\n', fp);
                }
            }
            puts ("done");
        }
    }
    fclose (fp);
}

/*---------- check_solution -------------------------------------------
 * check zone for a complete solution
 *---------------------------------------------------------------------*/

static int check_solution (int nz)
{
    int nf, flags = 0;
    SOLUTION *sol;

    if (!read_zone_solution (nz+1) ||
        usesol > Zones[nz].nsols) return 0;
    sol = &Zones[nz].sols[usesol-1];
    if (sol->nflds < 5 || (sol->location != CGNS_ENUMV(Vertex) &&
        sol->location != CGNS_ENUMV(CellCenter))) return 0;
    for (nf = 0; nf < sol->nflds; nf++) {
        if (!strcmp (sol->flds[nf].name, "Density")) {
            flags |= 0x01;
            continue;
        }
        if (!strcmp (sol->flds[nf].name, "VelocityX") ||
            !strcmp (sol->flds[nf].name, "MomentumX")) {
            flags |= 0x02;
            continue;
        }
        if (!strcmp (sol->flds[nf].name, "VelocityY") ||
            !strcmp (sol->flds[nf].name, "MomentumY")) {
            flags |= 0x04;
            continue;
        }
        if (!strcmp (sol->flds[nf].name, "VelocityZ") ||
            !strcmp (sol->flds[nf].name, "MomentumZ")) {
            flags |= 0x08;
            continue;
        }
        if (!strcmp (sol->flds[nf].name, "Pressure") ||
            !strcmp (sol->flds[nf].name, "EnergyStagnationDensity")) {
            flags |= 0x10;
            continue;
        }
    }
    return (flags & 0x1f) == 0x1f ? 1 : 0;
}

/*---------- get_reference --------------------------------------------
 * get the reference conditions
 *---------------------------------------------------------------------*/

static void get_reference (void)
{
    int n, narrays, na, dim;
    cgsize_t vec[12];
    CGNS_ENUMT(DataType_t) datatype;
    CGNS_ENUMT(AngleUnits_t) angle;
    int aloc = 0, units[5];
    char name[33];
    static char *refnames[4] = {
        "Mach",
        "AngleofAttack",
        "Reynolds",
        "TimeLatest"
    };

    reference[0] = 0.5;   /* Mach Number */
    reference[1] = 0.0;   /* angle of attack */
    reference[2] = 1.0e6; /* Reynolds Number */
    reference[3] = 0.0;   /* time */

    if (cg_goto (cgnsfn, cgnsbase, "ReferenceState_t", 1, "end") ||
        cg_narrays (&narrays) || narrays < 1)
        return;
    for (na = 1; na <= narrays; na++) {
        if (cg_array_info (na, name, &datatype, &dim, vec))
            FATAL ("get_reference", NULL);
        if (dim != 1 || vec[0] != 1) continue;
        for (n = 0; n < 4; n++) {
            if (!strcmp (refnames[n], name)) {
                if (cg_array_read_as (na, CGNS_ENUMV(RealDouble), &reference[n]))
                    FATAL ("get_reference", NULL);
                if (n == 1) aloc = na;
                break;
            }
        }
    }

    /* angle of attack units */

    if (aloc) {
        if (cg_goto (cgnsfn, cgnsbase, "ReferenceState_t", 1,
            "DataArray_t", aloc, "end"))
            FATAL ("get_reference", NULL);
        if (read_units (units))
            angle = (CGNS_ENUMT(AngleUnits_t))units[4];
        else
            angle = (CGNS_ENUMT(AngleUnits_t))baseunits[4];
        if (angle == CGNS_ENUMV(Radian))
            reference[1] *= 57.29578;
    }
}

/*---------- compute_solution -----------------------------------------
 * compute solution for a zone
 *---------------------------------------------------------------------*/

static void compute_solution (int nz)
{
    int n, nf, loc[5], con[5];
    double vel2;
    SOLUTION *sol = &Zones[nz].sols[usesol-1];

    for (n = 0; n < 5; n++)
        loc[n] = -1;
    for (nf = 0; nf < sol->nflds; nf++) {
        if (!strcmp (sol->flds[nf].name, "Density")) {
            loc[0] = nf;
            con[0] = 0;
        }
        else if (!strcmp (sol->flds[nf].name, "VelocityX")) {
            if (loc[1] >= 0) continue;
            loc[1] = nf;
            con[1] = 1;
        }
        else if (!strcmp (sol->flds[nf].name, "MomentumX")) {
            loc[1] = nf;
            con[1] = 0;
        }
        else if (!strcmp (sol->flds[nf].name, "VelocityY")) {
            if (loc[2] >= 0) continue;
            loc[2] = nf;
            con[2] = 1;
        }
        else if (!strcmp (sol->flds[nf].name, "MomentumY")) {
            loc[2] = nf;
            con[2] = 0;
        }
        else if (!strcmp (sol->flds[nf].name, "VelocityZ")) {
            if (loc[3] >= 0) continue;
            loc[3] = nf;
            con[3] = 1;
        }
        else if (!strcmp (sol->flds[nf].name, "MomentumZ")) {
            loc[3] = nf;
            con[3] = 0;
        }
        else if (!strcmp (sol->flds[nf].name, "Pressure")) {
            if (loc[4] >= 0) continue;
            loc[4] = nf;
            con[4] = 1;
        }
        else if (!strcmp (sol->flds[nf].name, "EnergyStagnationDensity")) {
            loc[4] = nf;
            con[4] = 0;
        }
        else
            continue;
        read_solution_field (nz+1, usesol, nf+1);
    }
    if (sol->location != CGNS_ENUMV(Vertex))
        cell_vertex_solution (nz+1, usesol, weighting);

    for (nf = 0; nf < 5; nf++) {
        for (n = 0; n < sol->size; n++)
            q[nf][n] = sol->flds[loc[nf]].data[n];
    }
    for (nf = 1; nf <= 3; nf++) {
        if (con[nf]) {
            for (n = 0; n < sol->size; n++)
                q[nf][n] *= q[0][n];
        }
    }
    if (con[4]) {
        for (n = 0; n < sol->size; n++) {
            vel2 = 0.0;
            for (nf = 1; nf <= 3; nf++)
                vel2 += (q[nf][n] * q[nf][n]);
            q[4][n] = q[4][n] / (gamma - 1.0) + 0.5 * vel2 / q[0][n];
        }
    }
}

/*---------- write_q_binary -------------------------------------------
 * write binary Plot3d file
 *---------------------------------------------------------------------*/

static void write_q_binary (char *qfile)
{
    int i, j, n, k, nk, nz, np, dim[3];
    float qf[4];
    FILE *fp;

    printf ("\nwriting binary Q file to %s\n", qfile);
    printf ("  in %s-precision\n", use_double ? "double" : "single");

    if (NULL == (fp = fopen (qfile, "w+b"))) {
        fprintf (stderr, "couldn't open <%s> for writing\n", qfile);
        exit (1);
    }
    if (mblock || nblocks > 1)
        fwrite (&nblocks, sizeof(int), 1, fp);
    for (nz = 0; nz < nZones; nz++) {
	if (Zones[nz].type == CGNS_ENUMV(Structured)) {
	    for (i = 0; i < 3; i++)
		dim[i] = (int)Zones[nz].dim[i];
            fwrite (dim, sizeof(int), 3, fp);
	}
    }
    for (nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            printf ("  zone %d ... ", nz+1);
            fflush (stdout);
            compute_solution (nz);
            if (whole) {
                np = (int)Zones[nz].nverts;
                nk = 1;
            }
            else {
                np = (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
                nk = (int)Zones[nz].dim[2];
            }
            if (use_double) {
                fwrite (reference, sizeof(double), 4, fp);
                for (k = 0; k < nk; k++) {
                    for (i = 0; i < 5; i++)
                        fwrite (&q[i][k*np], sizeof(double), np, fp);
                }
            }
            else {
                 for (n = 0; n < 4; n++)
                     qf[n] = (float)reference[n];
                fwrite (qf, sizeof(float), 4, fp);
                for (k = 0; k < nk; k++) {
                    for (i = 0; i < 5; i++) {
                        j = k * np;
                        for (n = 0; n < np; n++, j++) {
                            qf[0] = (float)q[i][j];
                            fwrite (qf, sizeof(float), 1, fp);
                        }
                    }
                }
            }
            puts ("done");
        }
    }
    fclose (fp);
}

/*---------- write_q_unformatted --------------------------------------
 * write unformatted Plot3d file
 *---------------------------------------------------------------------*/

static void write_q_unformatted (char *qfile)
{
    int np, nk, nz, nq, ierr;
    int i, j, k, n, *indices;
    char buff[129];
    void *qdata;
    float *qf = 0;
    double *qd = 0;

    printf ("\nwriting unformatted Q file to %s\n", qfile);
    printf ("  in %s-precision\n", use_double ? "double" : "single");

    for (np = 0, nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            nk = whole ? (int)Zones[nz].nverts :
                (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
            if (np < nk) np = nk;
        }
    }
    indices = (int *) malloc (3 * nblocks * sizeof(int));
    if (use_double) {
        qdata = (void *) malloc (5 * np * sizeof(double));
        qd = (double *)qdata;
    }
    else {
        qdata = (void *) malloc (5 * np * sizeof(float));
        qf = (float *)qdata;
    }
    if (NULL == indices || NULL == qdata)
        FATAL ("write_q_unformatted", "malloc failed for working arrays");

    unlink (qfile);
    strcpy (buff, qfile);
    for (n = (int)strlen(buff); n < 128; n++)
        buff[n] = ' ';
    buff[128] = 0;
    n = 0;
    OPENF (&n, buff, 128);

    if (mblock || nblocks > 1) {
        n = 1;
        WRITEIF (&n, &nblocks, &ierr);
    }
    for (np = 0, nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            for (n = 0; n < 3; n++)
                indices[np++] = (int)Zones[nz].dim[n];
        }
    }
    WRITEIF (&np, indices, &ierr);

    for (nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            printf ("  zone %d ... ", nz+1);
            fflush (stdout);
            compute_solution (nz);
            if (whole) {
                np = (int)Zones[nz].nverts;
                nk = 1;
            }
            else {
                np = (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
                nk = (int)Zones[nz].dim[2];
            }
            if (use_double) {
                n = 4;
                WRITEDF (&n, reference, &ierr);
                for (k = 0; k < nk; k++) {
                    for (nq = 0, j = 0; j < 5; j++) {
                        i = k * np;
                        for (n = 0; n < np; n++, i++)
                            qd[nq++] = q[j][i];
                    }
                    WRITEDF (&nq, qd, &ierr);
                }
            }
            else {
                for (n = 0; n < 4; n++)
                    qf[n] = (float)reference[n];
                WRITEFF (&n, qf, &ierr);
                for (k = 0; k < nk; k++) {
                    for (nq = 0, j = 0; j < 5; j++) {
                        i = k * np;
                        for (n = 0; n < np; n++, i++)
                            qf[nq++] = (float)q[j][i];
                    }
                    WRITEFF (&nq, qf, &ierr);
                }
            }
            puts ("done");
        }
    }
    CLOSEF ();
    free (indices);
    free (qdata);
}

/*---------- write_q_formatted ----------------------------------------
 * write formatted Plot3d file
 *---------------------------------------------------------------------*/

static void write_q_formatted (char *qfile)
{
    int nz, i, j, n, k, nk, np;
    FILE *fp;

    printf ("\nwriting formatted Q file to %s\n", qfile);

    if (NULL == (fp = fopen (qfile, "w+"))) {
        fprintf (stderr, "couldn't open <%s> for writing\n", qfile);
        exit (1);
    }
    if (mblock || nblocks > 1)
        fprintf (fp, "%d\n", nblocks);
    for (nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured))
            fprintf (fp, "%d %d %d\n", (int)Zones[nz].dim[0],
                (int)Zones[nz].dim[1], (int)Zones[nz].dim[2]);
    }
    for (nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            printf ("  zone %d ... ", nz+1);
            fflush (stdout);
            fprintf (fp, "%#g %#g %#g %#g\n", reference[0], reference[1],
                reference[2], reference[3]);
            compute_solution (nz);
            if (whole) {
                np = (int)Zones[nz].nverts;
                nk = 1;
            }
            else {
                np = (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
                nk = (int)Zones[nz].dim[2];
            }
            for (k = 0; k < nk; k++) {
                i = k * np;
                for (j = 0; j < 5; j++) {
                    for (n = 0; n < np; n++) {
                        if (n) putc ((n % 5) == 0 ? '\n' : ' ', fp);
                        fprintf (fp, "%#g", q[j][i+n]);
                    }
                    putc ('\n', fp);
                }
            }
            puts ("done");
        }
    }
    fclose (fp);
}

/*========== main =====================================================*/

int main (int argc, char *argv[])
{
    int n, ib, nb, is, nz, celldim, phydim;
    cgsize_t imax;
    char basename[33];

    if (argc < 2)
        print_usage (usgmsg, NULL);

    ib = 0;
    basename[0] = 0;
    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 's':
                mblock = 0;
                break;
            case 'p':
                whole = 0;
                break;
            case 'n':
                use_iblank = 0;
                break;
            case 'f':
            case 'u':
                format = n;
                break;
            case 'd':
                use_double = 1;
                break;
            case 'b':
                ib = atoi (argarg);
                break;
            case 'B':
                strncpy (basename, argarg, 32);
                basename[32] = 0;
                break;
            case 'g':
                gamma = atof (argarg);
                if (gamma <= 1.0)
                    FATAL (NULL, "invalid value for gamma");
                break;
            case 'w':
                weighting = 1;
                break;
            case 'S':
                usesol = atoi (argarg);
                break;
        }
    }

    if (argind > argc - 2)
        print_usage (usgmsg, "CGNSfile and/or XYZfile not given");
    if (!file_exists (argv[argind]))
        FATAL (NULL, "CGNSfile does not exist or is not a file");

    /* open CGNS file */

    printf ("reading CGNS file from %s\n", argv[argind]);
    nb = open_cgns (argv[argind], 1);
    if (!nb)
        FATAL (NULL, "no bases found in CGNS file");
    if (*basename && 0 == (ib = find_base (basename)))
        FATAL (NULL, "specified base not found");
    if (ib > nb) FATAL (NULL, "base index out of range");
    cgnsbase = ib ? ib : 1;
    if (cg_base_read (cgnsfn, cgnsbase, basename, &celldim, &phydim))
        FATAL (NULL, NULL);
    if (celldim != 3 || phydim != 3)
        FATAL (NULL, "cell and/or physical dimension must be 3");
    printf ("  using base %d - %s\n", cgnsbase, basename);

    read_zones ();
    for (nz = 0; nz < nZones; nz++) {
	if (Zones[nz].type == CGNS_ENUMV(Structured)) {
	    /* verify we can write out using ints */
	    for (n = 0; n < 3; n++) {
		if (Zones[nz].dim[n] > CG_MAX_INT32)
		    FATAL(NULL, "zone dimensions too large for integer");
	    }
	    if (whole) {
		if (Zones[nz].nverts > CG_MAX_INT32)
		    FATAL(NULL, "zone too large to write as whole using an integer");
	    }
	    else {
	        if (Zones[nz].dim[0]*Zones[nz].dim[1] > CG_MAX_INT32)
		    FATAL(NULL, "zone too large to write using an integer");
	    }
	    nblocks++;
	}
    }
    if (!nblocks) FATAL (NULL, "no structured zones found");

    /* read the nodes */

    printf ("reading %d zones\n", nblocks);
    ib = is = 0;
    imax = 0;
    for (nz = 0; nz < nZones; nz++) {
        if (Zones[nz].type == CGNS_ENUMV(Structured)) {
            printf ("  zone %d - %s ... ", nz+1, Zones[nz].name);
            fflush (stdout);
            read_zone_grid (nz+1);
            ib += read_zone_interface (nz+1);
            is += check_solution (nz);
            if (imax < Zones[nz].nverts) imax = Zones[nz].nverts;
            puts ("done");
        }
    }

    if (!ib) use_iblank = 0;
    if (use_iblank) {
        iblank = (int *) malloc ((size_t)imax * sizeof(int));
        if (NULL == iblank)
            FATAL (NULL, "malloc failed for iblank array");
    }

    /* write Plot3d XYZ file */

    if (format == 'f')
        write_xyz_formatted (argv[++argind]);
    else if (format == 'u')
        write_xyz_unformatted (argv[++argind]);
    else
        write_xyz_binary (argv[++argind]);

    if (use_iblank) free (iblank);

    /* write solution file */

    if (++argind < argc) {
        if (is != nblocks) {
            fprintf (stderr, "solution file is not being written since not\n");
            fprintf (stderr, "all the blocks contain a complete solution\n");
            cg_close (cgnsfn);
            exit (1);
        }
        for (n = 0; n < 5; n++) {
            q[n] = (double *) malloc ((size_t)imax * sizeof(double));
            if (NULL == q[n])
                FATAL (NULL, "malloc failed for solution working array");
        }
        get_reference ();
        if (format == 'f')
            write_q_formatted (argv[argind]);
        else if (format == 'u')
            write_q_unformatted (argv[argind]);
        else
            write_q_binary (argv[argind]);
    }

    cg_close (cgnsfn);
    return 0;
}
