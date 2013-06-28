/*
 * cgns_info.c - print CGNS information
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "getargs.h"
#include "cgnslib.h"
#include "cgnsutil.h"

/* command line options */

static char options[] = "vb:";

static char *usgmsg[] = {
    "usage  : cgns_info [options] CGNSfile",
    "options:",
    "   -v       = verbose output",
    "   -b<base> = use CGNS base number <base> (default 1)",
    NULL
};

static void print_interface (ZONE *zone)
{
    int ni;
    INTERFACE *ints;

    read_zone_interface (zone->id);
    ints = zone->ints;

    for (ni = 0; ni < zone->nints; ni++) {
        printf ("  interface %d - %s\n", ints[ni].id, ints[ni].name);
        printf ("    range       = %d %d %d %d %d %d\n",
            (int)ints[ni].range[0][0], (int)ints[ni].range[0][1],
            (int)ints[ni].range[1][0], (int)ints[ni].range[1][1],
            (int)ints[ni].range[2][0], (int)ints[ni].range[2][1]);
        printf ("    donor name  = %s\n", ints[ni].d_name);
        printf ("    donor range = %d %d %d %d %d %d\n",
            (int)ints[ni].d_range[0][0], (int)ints[ni].d_range[0][1],
            (int)ints[ni].d_range[1][0], (int)ints[ni].d_range[1][1],
            (int)ints[ni].d_range[2][0], (int)ints[ni].d_range[2][1]);
        printf ("    transform   = %d %d %d\n", ints[ni].transform[0],
            ints[ni].transform[1], ints[ni].transform[2]);
        printf ("    donor zone  = %d\n", ints[ni].d_zone);
    }
}

static void print_connect (ZONE *zone)
{
    int nc;
    CONNECT *conns;

    read_zone_connect (zone->id);
    conns = zone->conns;

    for (nc = 0; nc < zone->nconns; nc++) {
        printf ("  connectivity %d - %s\n", conns[nc].id, conns[nc].name);
        printf ("    type          = %d\n", conns[nc].type);
        printf ("    location      = %d\n", conns[nc].location);
        printf ("    pt type       = %d\n", conns[nc].ptype);
        printf ("    points        = %d\n", (int)conns[nc].npnts);
        printf ("    donor name    = %s\n", conns[nc].d_name);
        printf ("    donor pt type = %d\n", conns[nc].d_ptype);
        printf ("    donor points  = %d\n", (int)conns[nc].d_npnts);
        printf ("    donor zone    = %d\n", conns[nc].d_zone);
    }
}

static void print_solution (ZONE *zone)
{
    int ns, nf;
    SOLUTION *sols;

    read_zone_solution (zone->id);
    sols = zone->sols;

    for (ns = 0; ns < zone->nsols; ns++) {
        printf ("  solution %d - %s\n", sols[ns].id, sols[ns].name);
        printf ("    location = %d\n", sols[ns].location);
        printf ("    rind     = %d %d %d %d %d %d\n",
            sols[ns].rind[0][0], sols[ns].rind[0][1],
            sols[ns].rind[1][0], sols[ns].rind[1][1],
            sols[ns].rind[2][0], sols[ns].rind[2][1]);
        printf ("    size     = %d\n", (int)sols[ns].size);
        printf ("    fields   = %d\n", sols[ns].nflds);
        for (nf = 0; nf < sols[ns].nflds; nf++)
            printf ("      %s\n", sols[ns].flds[nf].name);
    }
}

int main (int argc, char *argv[])
{
    int n, nz, nbases, celldim, phydim;
    int verbose = 0;
    char basename[33];
    ZONE *z;

    if (argc < 2)
        print_usage (usgmsg, NULL);

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'v':
                verbose = 1;
                break;
            case 'b':
                cgnsbase = atoi (argarg);
                break;
        }
    }

    if (argind == argc)
        print_usage (usgmsg, "CGNSfile not given");
    if (!file_exists (argv[argind]))
        FATAL (NULL, "CGNSfile does not exist or is not a file");

    /* read CGNS file */

    printf ("reading CGNS file from %s\n", argv[argind]);
    fflush (stdout);
    nbases = open_cgns (argv[argind], 1);
    if (!nbases) FATAL (NULL, "no bases in CGNS file");
    if (cgnsbase < 1 || cgnsbase > nbases)
        FATAL (NULL, "invalid base index");

    if (cg_base_read (cgnsfn, cgnsbase, basename, &celldim, &phydim))
        FATAL (NULL, NULL);
    printf ("using base %d - %s\n", cgnsbase, basename);
    printf ("cell dimension     = %d\n", celldim);
    printf ("physical dimension = %d\n", phydim);

    read_zones ();

    for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
        printf ("\nzone %d - %s\n", z->id, z->name);
        printf ("type      = %d\n", z->type);
        printf ("dimension = %d x %d x %d\n", (int)z->dim[0],
            (int)z->dim[1], (int)z->dim[2]);
        printf ("1to1      = %d\n", z->nints);
        if (verbose) print_interface (z);
        printf ("connects  = %d\n", z->nconns);
        if (verbose) print_connect (z);
        printf ("solutions = %d\n", z->nsols);
        if (verbose) print_solution (z);
    }

    if (cg_close (cgnsfn))
        FATAL (NULL, NULL);

    return 0;
}
