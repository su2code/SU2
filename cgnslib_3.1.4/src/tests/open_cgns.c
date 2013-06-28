#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

int main (int argc, char **argv)
{
    double start, finish;
    int cgfile, mode = CG_MODE_READ;

    if (argc < 2 || argc > 3) {
        fprintf (stderr, "open_cgns [-m] CGNSfile\n");
        exit (1);
    }
    if (argc > 2) {
        mode = CG_MODE_MODIFY;
        cg_configure (CG_CONFIG_COMPRESS, (void *)1);
    }

    printf ("opening cgns file <%s> ...", argv[argc-1]);
    fflush (stdout);
    start = elapsed_time ();
    if (cg_open (argv[argc-1], mode, &cgfile)) cg_error_exit();
    finish = elapsed_time ();
    printf (" %g secs\n", finish - start);

    printf ("closing cgns file ...");
    fflush (stdout);
    start = elapsed_time ();
    cg_close (cgfile);
    finish = elapsed_time ();
    printf (" %g secs\n", finish - start);

    return 0;
}

