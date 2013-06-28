#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
# include <io.h>
# define unlink _unlink
#else
# include <unistd.h>
#endif
#include "utils.h"

int cgfile, cgbase, cgzone, cgcoord;
cgsize_t size[9];

#define NUM_SIDE 5

float coord[NUM_SIDE*NUM_SIDE*NUM_SIDE];

int main (int argc, char **argv)
{
    int n, nz, nzones = 50;
    double start, finish;
    char name[33], linkpath[33];
    char fname[33], linkfile[33];
    int cgftop,cgfchild;

    for (n = 0; n < 3; n++) {
        size[n]   = NUM_SIDE;
        size[n+3] = NUM_SIDE - 1;
        size[n+6] = 0;
    }
    if (argc > 1)
        nzones = atoi (argv[1]);
    printf ("number of zones  = %d\n", nzones);

    for (nz = 1; nz <= nzones; nz++) {
        sprintf (fname, "zone%d.cgns", nz);
        unlink (fname);
    }

    printf ("creating zones ...");
    fflush (stdout);
    start = elapsed_time ();
    for (nz = 1; nz <= nzones; nz++) {
        sprintf (fname, "zone%d.cgns", nz);
        if (cg_open (fname, CG_MODE_WRITE, &cgfile) ||
            cg_base_write (cgfile, "Base", 3, 3, &cgbase) ||
            cg_zone_write (cgfile, cgbase, "Zone", size, CGNS_ENUMV( Structured ),
                &cgzone) ||
            cg_coord_write(cgfile, cgbase, cgzone, CGNS_ENUMV( RealSingle ),
                "CoordinateX", coord, &cgcoord) ||
            cg_coord_write(cgfile, cgbase, cgzone, CGNS_ENUMV( RealSingle ),
                "CoordinateY", coord, &cgcoord) ||
            cg_coord_write(cgfile, cgbase, cgzone, CGNS_ENUMV( RealSingle ),
                "CoordinateZ", coord, &cgcoord))
            cg_error_exit();
        if (cg_close(cgfile)) cg_error_exit();
    }
    finish = elapsed_time ();
    printf (" %g secs\n", finish - start);

    strcpy (fname, "links.cgns");
    unlink (fname);
    strcpy (linkpath, "/Base/Zone");
    printf ("creating link file ...");
    fflush (stdout);
    start = elapsed_time ();
    if (cg_open (fname, CG_MODE_WRITE, &cgfile) ||
        cg_base_write (cgfile, "Base", 3, 3, &cgbase))
        cg_error_exit();
    for (nz = 1; nz <= nzones; nz++) {
        sprintf (name, "Zone%d", nz);
        sprintf (linkfile, "zone%d.cgns", nz);
        if (cg_goto (cgfile, cgbase, "end") ||
            cg_link_write (name, linkfile, linkpath))
            cg_error_exit();
    }
    cg_close (cgfile);
    finish = elapsed_time ();
    printf (" %g secs\n", finish - start);
    printf ("file size        = %g Mb\n", file_size(fname));

    printf ("opening top link file  ...");
    fflush (stdout);
    if (cg_open (fname, CG_MODE_READ, &cgftop)) cg_error_exit();
    cg_nzones(cgftop,1,&n);
    cg_close (cgftop);
    printf ("close top link file\n");

    printf ("opening child link file ...");
    fflush (stdout);
    if (cg_open ("zone1.cgns", CG_MODE_READ, &cgfchild)) cg_error_exit();
    cg_nzones(cgfchild,1,&n);
    printf (" [%d] zones\n",n);

    printf ("opening top link file (child still open) ...");
    fflush (stdout);
    if (cg_open (fname, CG_MODE_READ, &cgftop)) cg_error_exit();
    cg_nzones(cgftop,1,&n);
    printf (" [%d] zones...",n);
    cg_close (cgftop);
    printf ("close top link file\n");

    printf ("read child link file ...");
    cg_nzones(cgfchild,1,&n);
    printf (" [%d] zones...",n);
    cg_close (cgfchild);
    printf ("close child link file\n");

    return 0;
}

