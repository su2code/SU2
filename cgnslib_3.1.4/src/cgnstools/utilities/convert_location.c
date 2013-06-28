/*
 * convert_location.c - convert between vertex and cell-center
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

#ifndef CG_MODE_MODIFY
# define CG_MODE_MODIFY MODE_MODIFY
#endif

/* command line options */

#if defined(CELL_TO_VERTEX)

static char options[] = "wb:z:s:S:";

static char *usgmsg[] = {
    "usage  : cell_to_vertex [options] CGNSfile [newCGNSfile]",
    "options:",
    "   -w       = use volume weighting",
    "   -b<base> = use CGNS base number <base> (default 1)",
    "   -z<zone> = read zone number <zone> (default all)",
    "   -s<sol>  = read solution number <sol> (default all)",
    "   -S<name> = write to solution <name> (default same as read)",
    NULL
};

#elif defined(VERTEX_TO_CELL)

static char options[] = "wb:z:s:S:ijkIJK";

static char *usgmsg[] = {
    "usage  : vertex_to_cell [options] CGNSfile [newCGNSfile]",
    "options:",
    "   -w       = use volume weighting",
    "   -b<base> = use CGNS base number <base> (default 1)",
    "   -z<zone> = read zone number <zone> (default all)",
    "   -s<sol>  = read solution number <sol> (default all)",
    "   -S<name> = write to solution <name> (default same as read)",
    "   -i       = add rind cell at imin",
    "   -I       = add rind cell at imax",
    "   -j       = add rind cell at jmin",
    "   -J       = add rind cell at jmax",
    "   -k       = add rind cell at kmin",
    "   -K       = add rind cell at kmax",
    NULL
};

#else

static char options[] = "cvwb:z:s:S:ijkIJK";

static char *usgmsg[] = {
    "usage  : convert_location [options] CGNSfile [newCGNSfile]",
    "options:",
    "   -c       = convert to cell-center",
    "   -v       = convert to vertex",
    "   -w       = use volume weighting",
    "   -b<base> = use CGNS base number <base> (default 1)",
    "   -z<zone> = read zone number <zone> (default all)",
    "   -s<sol>  = read solution number <sol> (default all)",
    "   -S<name> = write to solution <name> (default same as read)",
    "   for conversions to cell-center:",
    "   -i       = add rind cell at imin",
    "   -I       = add rind cell at imax",
    "   -j       = add rind cell at jmin",
    "   -J       = add rind cell at jmax",
    "   -k       = add rind cell at kmin",
    "   -K       = add rind cell at kmax",
    NULL
};

#endif

static int weighting = 0;
static int izone = 0;
static int isol = 0;

/*-------------------------------------------------------------------*/

static void check_zones (int location)
{
    int iz, nz, ns, icnt = 0;
    ZONE *z;
    SOLUTION *s;

    for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
        iz = 0;
        if (!izone || nz == izone) {
            read_zone_solution (nz);
            for (s = z->sols, ns = 1; ns <= z->nsols; ns++, s++) {
                if ((!isol || ns == isol) && s->location == location) {
                    icnt += s->nflds;
                    iz = z->id;
                }
                else
                    s->id = 0;
            }
        }
        z->id = iz;
    }

    if (!icnt) {
        printf ("nothing to do\n");
        exit (0);
    }
}

/*-------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
    int i, j, n, nz, ns, celldim, phydim;
    int rind[3][2], location = 0;
    char basename[33], *solname = NULL, *tmpfile;
    ZONE *z;
    SOLUTION *s;

    if (argc < 2)
        print_usage (usgmsg, NULL);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 2; j++)
            rind[i][j] = 0;

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'c':
                location += CGNS_ENUMV(CellCenter);
                break;
            case 'v':
                location += CGNS_ENUMV(Vertex);
                break;
            case 'w':
                weighting = 1;
                break;
            case 'b':
                cgnsbase = atoi (argarg);
                break;
            case 'z':
                izone = atoi (argarg);
                break;
            case 's':
                isol = atoi (argarg);
                break;
            case 'S':
                solname = argarg;
                break;
            case 'i':
            case 'j':
            case 'k':
                rind[n-'i'][0] = 1;
                break;
            case 'I':
            case 'J':
            case 'K':
                rind[n-'I'][1] = 1;
                break;
        }
    }

    if (argind == argc)
        print_usage (usgmsg, "CGNSfile not given");
    if (!file_exists (argv[argind]))
        FATAL (NULL, "CGNSfile does not exist or is not a file");

#if defined(CELL_TO_VERTEX)
    location = CGNS_ENUMV(Vertex);
#elif defined(VERTEX_TO_CELL)
    location = CGNS_ENUMV(CellCenter);
#else
    if (!location)
        print_usage (usgmsg, "please select either the -c or -v option");
    if (location != CGNS_ENUMV(CellCenter) && location != CGNS_ENUMV(Vertex))
        print_usage (usgmsg, "please select only one of -c or -v options");
#endif

    /* create a working copy */

    printf ("creating a working copy of %s\n", argv[argind]);
    tmpfile = temporary_file (argv[argind]);
    copy_file (argv[argind], tmpfile);

    /* read CGNS file */

    printf ("reading CGNS file from %s\n", tmpfile);
    if (cg_open (tmpfile, CG_MODE_MODIFY, &cgnsfn) ||
        cg_base_read (cgnsfn, cgnsbase, basename, &celldim, &phydim))
        FATAL (NULL, NULL);
    if (celldim != 3 || phydim != 3)
        FATAL (NULL, "cell and/or physical dimension must be 3");
    printf ("reading zone information for base %d - %s\n",
        cgnsbase, basename);
    printf ("converting solution location to %s\n",
        location == CGNS_ENUMV(CellCenter) ? "CellCenter" : "Vertex");

    read_zones ();
    check_zones (location == CGNS_ENUMV(CellCenter) ? CGNS_ENUMV(Vertex) : CGNS_ENUMV(CellCenter));

    /* convert solution location */

    for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
        if (z->id) {
            for (s = z->sols, ns = 1; ns <= z->nsols; ns++, s++) {
                if (s->id) {
                    printf ("converting zone %d, solution %d ... ", nz, ns);
                    fflush (stdout);
                    read_solution_field (nz, ns, 0);
                    if (location == CGNS_ENUMV(CellCenter)) {
                        for (i = 0; i < 3; i++)
                            for (j = 0; j < 2; j++)
                                s->rind[i][j] = rind[i][j];
                        cell_center_solution (nz, ns, weighting);
                    }
                    else
                        cell_vertex_solution (nz, ns, weighting);
                    puts ("done");
                }
            }
        }
    }

    /* write CGNS file */

    for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
        if (z->id) {
            for (s = z->sols, ns = 1; ns <= z->nsols; ns++, s++) {
                if (s->id) {
                    printf ("writing zone %d, solution %d ... ", nz, ns);
                    fflush (stdout);
                    if (solname != NULL) {
                        if (z->nsols == 1)
                            strncpy (s->name, solname, 32);
                        else if (strlen (solname) > 30)
                            sprintf (s->name, "%30.30s%d", solname, s->id);
                        else
                            sprintf (s->name, "%s%d", solname, s->id);
                        s->name[32] = 0;
                    }
                    write_zone_solution (nz, ns);
                    write_solution_field (nz, ns, 0);
                    puts ("done");
                }
            }
        }
    }

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
