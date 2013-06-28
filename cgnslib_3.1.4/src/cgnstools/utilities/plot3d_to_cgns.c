/*
 * plot3d_to_cgns.c - convert a PLOT3D file to CGNS
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "getargs.h"
#include "binaryio.h"
#include "cgnslib.h"
#include "cgnsutil.h"

/* command line options */

static char options[] = "spinfuM:db:B:g:c";

static char *usgmsg[] = {
    "usage: plot3d_to_cgns [options] XYZfile [Qfile] CGNSfile",
    "options are :",
    "   -s       = single block file (default multi-block)",
    "   -p       = planar grid format (default whole format)",
    "   -i       = has iblank array",
    "   -n       = read iblank array, but ignore it",
    "   -f       = formatted (ASCII) Plot3d file format",
    "   -u       = Fortran unformatted Plot3d file format",
    "   -M<mach> = machine type for binary or unformatted - one of:",
    "                ieee  bsieee  iris  alpha   hp  ibm",
    "                sun   dec     cray  convex  nt  linux",
    "              case is not significant and only the first",
    "              2 characters are needed.",
    "   -d       = use double-precision (64-bit)",
    "   -b<base> = use CGNS base index <base>",
    "   -B<name> = set CGNS base name to <name>",
    "   -c       = convert solution to primitive variables",
    "   -g<gamma>= gamma for data conversions (default 1.4)",
    "Default is multi-block binary format with no iblank data.",
    NULL
};

static int mblock = 1;
static int whole = 1;
static int has_iblank = 0;
static int use_iblank = 0;
static int is_double = 0;

static double reference[4];
static double gamma = 1.4;
static int convert = 0;

/*---------- get_machine --------------------------------------------
 * get machine type from command line option
 *-------------------------------------------------------------------*/

static int get_machine (char *name)
{
    switch (*name++) {
        case 'i':
        case 'I':
            if ('r' == *name || 'R' == *name)
                return (MACH_IRIS);
            if ('b' == *name || 'B' == *name)
                return (MACH_IBM);
            return (MACH_IEEE);
        case 'b':
        case 'B':
            return (MACH_BSIEEE);
        case 's':
        case 'S':
            return (MACH_SUN);
        case 'h':
        case 'H':
            return (MACH_HP);
        case 'a':
        case 'A':
            return (MACH_ALPHA);
        case 'd':
        case 'D':
            return (MACH_DEC);
        case 'c':
        case 'C':
            if ('r' == *name || 'R' == *name)
                return (MACH_CRAY);
            return (MACH_CONVEX);
        case 'l':
        case 'L':
            return (MACH_LINUX);
        case 'n':
            return (MACH_WIN32);
        default:
            break;
    }
    print_usage (usgmsg, "unknown machine name - option (-M)");
    return (0); /* quite compiler */
}

/*---------- read_xyz -----------------------------------------------
 * read PLOT3D XYZ file
 *-------------------------------------------------------------------*/

static void read_xyz (BINARYIO *bf)
{
    int i, k, n, nk, nz, np, nmax;
    int *indices, *iblank;
    void *xyz;
    VERTEX *verts;

    /* get number of grid blocks */

    if (mblock) {
        bf_getints (bf, 1, &nZones);
        if (nZones < 1 || nZones > 100000) {
            fprintf (stderr, "found %d blocks\n", nZones);
            fprintf (stderr, "file type and/or format is probably incorrect\n");
            bf_close (bf);
            exit (1);
        }
    }
    else
        nZones = 1;
    printf ("reading %d grid blocks\n", nZones);

    /* read indices for grids */

    indices = (int *) malloc (3 * nZones * sizeof(int));
    if (NULL == indices)
        FATAL ("read_xyz", "malloc failed for grid indices");
    bf_getints (bf, 3 * nZones, indices);

    /* create zone structures */

    Zones = new_zone (nZones);
    for (nmax = 0, nz = 0; nz < nZones; nz++) {
        Zones[nz].type = CGNS_ENUMV(Structured);
        for (np = 1, n = 0; n < 3; n++) {
            nk = indices[3 * nz + n];
            Zones[nz].dim[n] = nk;
            np *= nk;
        }
        Zones[nz].vertflags = 7;
        Zones[nz].datatype = is_double ? CGNS_ENUMV(RealDouble) : CGNS_ENUMV(RealSingle);
        Zones[nz].nverts = np;
        Zones[nz].verts = new_vertex (np);
        nk = whole ? (int)Zones[nz].nverts : (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
        if (nmax < nk) nmax = nk;
    }

    free (indices);

    if (is_double)
        xyz = (void *) malloc (3 * nmax * sizeof(double));
    else
        xyz = (void *) malloc (3 * nmax * sizeof(float));
    if (NULL == xyz)
        FATAL ("read_xyz", "malloc failed for coordinate working array");
    if (has_iblank) {
        iblank = (int *) malloc (nmax * sizeof(int));
        if (NULL == iblank)
            FATAL ("read_xyz", "malloc failed for iblank array");
    }
    else
        use_iblank = 0;

    /* read the grid blocks */

    for (nz = 0; nz < nZones; nz++) {
        printf ("reading block %d grid %dx%dx%d ...", nz+1,
            (int)Zones[nz].dim[0], (int)Zones[nz].dim[1],
            (int)Zones[nz].dim[2]);
        fflush (stdout);
        if (whole) {
            np = (int)Zones[nz].nverts;
            nk = 1;
        }
        else {
            np = (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
            nk = (int)Zones[nz].dim[2];
        }
        verts = Zones[nz].verts;
        for (k = 0; k < nk; k++) {
            if (is_double)
                bf_getdoubles (bf, 3 * np, xyz);
            else
                bf_getfloats (bf, 3 * np, xyz);
            if (has_iblank)
                bf_getints (bf, np, iblank);
            if (is_double) {
                for (i = 0, n = 0; n < np; n++, i++)
                    verts[n].x = ((double *)xyz)[i];
                for (n = 0; n < np; n++, i++)
                    verts[n].y = ((double *)xyz)[i];
                for (n = 0; n < np; n++, i++)
                    verts[n].z = ((double *)xyz)[i];
            }
            else {
                for (i = 0, n = 0; n < np; n++, i++)
                    verts[n].x = ((float *)xyz)[i];
                for (n = 0; n < np; n++, i++)
                    verts[n].y = ((float *)xyz)[i];
                for (n = 0; n < np; n++, i++)
                    verts[n].z = ((float *)xyz)[i];
            }
            for (n = 0; n < np; n++, verts++)
                verts->id = use_iblank ? iblank[n] : 1;
        }
        puts (" done");
    }

    free (xyz);
    if (has_iblank) free (iblank);
}

/*---------- read_q -------------------------------------------------
 * read PLOT3D solution file
 *-------------------------------------------------------------------*/

static void read_q (BINARYIO *bf)
{
    int i, k, n, nk, nz, np, nv, nmax;
    int *indices;
    void *data;
    SOLUTION *sol;
    FIELD *flds;
    double qq, rho;
    static char *fldnames[] = {
        "Density",
        "MomentumX",
        "MomentumY",
        "MomentumZ",
        "EnergyStagnationDensity",
        "VelocityX",
        "VelocityY",
        "VelocityZ",
        "Pressure"
    };

    /* get number of grid blocks */

    if (mblock) {
        bf_getints (bf, 1, &nz);
        if (nz != nZones)
            FATAL ("read_q", "number of blocks not the same as the XYZ file");
    }

    /* read indices for grids */

    indices = (int *) malloc (3 * nZones * sizeof(int));
    if (NULL == indices)
        FATAL ("read_q", "malloc failed for grid indices");
    bf_getints (bf, 3 * nZones, indices);
    for (nz = 0; nz < nZones; nz++) {
        for (n = 0; n < 3; n++) {
            if (indices[3*nz+n] != Zones[nz].dim[n])
                FATAL ("read_q", "mismatch in block sizes");
        }
    }
    free (indices);

    /* create solution data arrays */

    for (nmax = 0, nz = 0; nz < nZones; nz++) {
        np = whole ? (int)Zones[nz].nverts : (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
        if (nmax < np) nmax = np;
        sol = new_solution (1);
        strcpy (sol->name, "FlowSolution");
        sol->location = CGNS_ENUMV(Vertex);
        sol->size = Zones[nz].nverts;
        sol->nflds = 5;
        sol->flds = new_field (5, sol->size);
        for (nv = 0; nv < 5; nv++) {
            strcpy (sol->flds[nv].name, fldnames[nv]);
            sol->flds[nv].datatype = is_double ? CGNS_ENUMV(RealDouble) : CGNS_ENUMV(RealSingle);
        }
        Zones[nz].nsols = 1;
        Zones[nz].sols = sol;
    }

    if (nmax < 4) nmax = 4;
    if (is_double)
        data = (void *) malloc (nmax * sizeof(double));
    else
        data = (void *) malloc (nmax * sizeof(float));
    if (NULL == data)
        FATAL ("read_q", "malloc failed for solution working array");

    /* read the solution data */

    for (nz = 0; nz < nZones; nz++) {
        printf ("reading block %d solution ...", nz+1);
        fflush (stdout);
        if (is_double) {
            bf_getdoubles (bf, 4, data);
            if (0 == nz) {
                for (n = 0; n < 4; n++)
                    reference[n] = ((double *)data)[n];
            }
        }
        else {
            bf_getfloats (bf, 4, data);
            if (0 == nz) {
                for (n = 0; n < 4; n++)
                    reference[n] = ((float *)data)[n];
            }
        }

        if (whole) {
            np = (int)Zones[nz].nverts;
            nk = 1;
        }
        else {
            np = (int)(Zones[nz].dim[0] * Zones[nz].dim[1]);
            nk = (int)Zones[nz].dim[2];
        }
        flds = Zones[nz].sols->flds;
        for (k = 0; k < nk; k++) {
            i = k * np;
            for (nv = 0; nv < 5; nv++) {
                if (is_double) {
                    bf_getdoubles (bf, np, data);
                    for (n = 0; n < np; n++)
                        flds[nv].data[n+i] = ((double *)data)[n];
                }
                else {
                    bf_getfloats (bf, np, data);
                    for (n = 0; n < np; n++)
                        flds[nv].data[n+i] = ((float *)data)[n];
                }
            }
        }

        if (convert) {
            for (nv = 1; nv < 5; nv++)
                strcpy (flds[nv].name, fldnames[4+nv]);
            for (n = 0; n < Zones[nz].nverts; n++) {
                rho = flds[0].data[n];
                for (qq = 0.0, nv = 1; nv < 4; nv++) {
                    flds[nv].data[n] /= rho;
                    qq += flds[nv].data[n] * flds[nv].data[n];
                }
                flds[4].data[n] = (gamma - 1.0) *
                    (flds[4].data[n] - 0.5 * rho * qq);
            }
        }
        puts (" done");
    }

    free (data);
}

/*---------- build_interfaces -----------------------------------------
 * create interfaces from iblank data
 *---------------------------------------------------------------------*/

static void build_interfaces (void)
{
}

/*---------- write_reference ------------------------------------------
 * write reference conditions to CGNS file
 *---------------------------------------------------------------------*/

static void write_reference (void)
{
    int n;
    cgsize_t cnt = 1;
    CGNS_ENUMT(DataType_t) datasize;
    float ref[4];
    void *mach, *alpha, *rey, *time;

    printf ("writing reference state...");
    fflush (stdout);
    if (cg_goto (cgnsfn, cgnsbase, "end") ||
        cg_state_write ("PLOT3D reference state"))
        FATAL ("write_reference", NULL);

    if (is_double) {
        datasize = CGNS_ENUMV(RealDouble);
        mach  = (void *)&reference[0];
        alpha = (void *)&reference[1];
        rey   = (void *)&reference[2];
        time  = (void *)&reference[3];
    }
    else {
        for (n = 0; n < 4; n++)
            ref[n] = (float)reference[n];
        datasize = CGNS_ENUMV(RealSingle);
        mach  = (void *)&ref[0];
        alpha = (void *)&ref[1];
        rey   = (void *)&ref[2];
        time  = (void *)&ref[3];
    }

    if (cg_goto (cgnsfn, cgnsbase, "ReferenceState_t", 1, "end") ||
        cg_array_write ("Mach", datasize, 1, &cnt, mach) ||
        cg_goto (cgnsfn, cgnsbase, "ReferenceState_t", 1,
            "DataArray_t", 1, "end") ||
        cg_dataclass_write (CGNS_ENUMV(NondimensionalParameter)))
        FATAL ("write_reference", NULL);

    if (cg_goto (cgnsfn, cgnsbase, "ReferenceState_t", 1, "end") ||
        cg_array_write ("AngleofAttack", datasize, 1, &cnt, alpha) ||
        cg_goto (cgnsfn, cgnsbase, "ReferenceState_t", 1,
            "DataArray_t", 2, "end") ||
        cg_dataclass_write (CGNS_ENUMV(Dimensional)) ||
        cg_units_write (CGNS_ENUMV(MassUnitsNull), CGNS_ENUMV(LengthUnitsNull), CGNS_ENUMV(TimeUnitsNull),
            CGNS_ENUMV(TemperatureUnitsNull), CGNS_ENUMV(Degree)))
        FATAL ("write_reference", NULL);

    if (cg_goto (cgnsfn, cgnsbase, "ReferenceState_t", 1, "end") ||
        cg_array_write ("Reynolds", datasize, 1, &cnt, rey) ||
        cg_goto (cgnsfn, cgnsbase, "ReferenceState_t", 1,
            "DataArray_t", 3, "end") ||
        cg_dataclass_write (CGNS_ENUMV(NondimensionalParameter)))
        FATAL ("write_reference", NULL);

    if (cg_goto (cgnsfn, cgnsbase, "ReferenceState_t", 1, "end") ||
        cg_array_write ("TimeLatest", datasize, 1, &cnt, time) ||
        cg_goto (cgnsfn, cgnsbase, "ReferenceState_t", 1,
            "DataArray_t", 4, "end") ||
        cg_dataclass_write (CGNS_ENUMV(Dimensional)) ||
        cg_units_write (CGNS_ENUMV(MassUnitsNull), CGNS_ENUMV(LengthUnitsNull),
            CGNS_ENUMV(Second), CGNS_ENUMV(TemperatureUnitsNull), CGNS_ENUMV(AngleUnitsNull)))
        FATAL ("write_reference", NULL);

    puts (" done");
}

/*========== main ===================================================*/

int main (int argc, char *argv[])
{
    int n, ib = 0, nb, flags = 0, has_q = 0;
    BINARYIO *bf;
    static char basename[33] = "Base";

    if (argc < 3)
        print_usage (usgmsg, NULL);

    /* get options */

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'f':
                flags &= ~OPEN_FORTRAN;
                flags |= OPEN_ASCII;
                break;
            case 'u':
                flags &= ~OPEN_ASCII;
                flags |= OPEN_FORTRAN;
                break;
            case 's':
                mblock = 0;
                break;
            case 'p':
                whole = 0;
                break;
            case 'i':
                use_iblank = 1;
                /* fall through */
            case 'n':
                has_iblank = 1;
                break;
            case 'd':
                is_double = 1;
                break;
            case 'M':
                flags &= ~MACH_UNKNOWN;
                flags |= get_machine (argarg);
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
            case 'c':
                convert = 1;
                break;
        }
    }

    if (argind > argc - 2)
        print_usage (usgmsg, "XYZfile and/or CGNSfile not given");

    /* read Plot3d file */

    printf ("reading PLOT3D grid file %s\n", argv[argind]);
    printf ("  as %s-block %s", mblock ? "multi" : "single",
        flags == OPEN_ASCII ? "ASCII" :
        (flags == OPEN_FORTRAN ? "FORTRAN unformatted" : "binary"));
    if (has_iblank) printf (" with iblank array");
    putchar ('\n');
    if (!file_exists (argv[argind]))
        FATAL (NULL, "XYZ file does not exist or is not a file");
    if (NULL == (bf = bf_open (argv[argind], flags | OPEN_READ))) {
        fprintf (stderr, "can't open <%s> for reading", argv[argind]);
        exit (1);
    }
    read_xyz (bf);
    bf_close (bf);

    if (use_iblank) build_interfaces ();

    /* read solution file if given */

    if (++argind < argc-1) {
        printf ("\nreading PLOT3D solution file %s\n", argv[argind]);
        if (!file_exists (argv[argind]))
            FATAL (NULL, "Solution file does not exist or is not a file");

        if (NULL == (bf = bf_open (argv[argind], flags | OPEN_READ))) {
            fprintf (stderr, "can't open <%s> for reading", argv[argind]);
            exit (1);
        }
        read_q (bf);
        bf_close (bf);
        argind++;
        has_q = 1;
    }

    /* open CGNS file */

    printf ("\nwriting CGNS file to %s\n", argv[argind]);
    nb = open_cgns (argv[argind], 0);
    if (ib) {
        if (ib > nb)
            FATAL (NULL, "specified base index out of range");
        if (cg_base_read (cgnsfn, ib, basename, &n, &n))
            FATAL (NULL, NULL);
    }
    if (cg_base_write (cgnsfn, basename, 3, 3, &cgnsbase) ||
        cg_goto (cgnsfn, cgnsbase, "end") ||
        cg_dataclass_write (CGNS_ENUMV(NormalizedByUnknownDimensional)))
        FATAL (NULL, NULL);
    printf ("  output to base %d - %s\n", cgnsbase, basename);

    write_zones ();
    for (n = 1; n <= nZones; n++) {
        printf ("writing zone %d ... grid", n);
        fflush (stdout);
        write_zone_grid (n);
        write_zone_interface (n);
        if (has_q) {
            printf (", solution");
            fflush (stdout);
            write_zone_solution (n, 1);
            write_solution_field (n, 1, 0);
        }
        puts (" done");
    }
    if (has_q) write_reference ();

    cg_close (cgnsfn);

    return 0;
}
