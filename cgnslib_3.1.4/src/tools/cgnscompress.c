#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "cgns_io.h"

int main (int argc, char **argv)
{
    char *inpfile, *outfile;
    int inpcg;
    size_t inpsize, outsize;
    struct stat st;

    if (argc < 2 || argc > 3) {
        fprintf(stderr, "usage: cgnscompress InputFile [OutputFile]\n");
        exit(1);
    }
    inpfile = argv[1];
    outfile = argv[argc-1];

    if (stat(inpfile, &st)) {
        fprintf (stderr, "can't stat %s\n", inpfile);
        exit (1);
    }
    inpsize = st.st_size;

    if (cgio_open_file (inpfile, 'r', CGIO_FILE_NONE, &inpcg))
        cgio_error_exit("cgio_open_file");
    if (cgio_compress_file (inpcg, outfile))
        cgio_error_exit("cgio_compress_file");

    if (stat(outfile, &st)) {
        fprintf (stderr, "can't stat %s\n", outfile);
        exit (1);
    }
    outsize = st.st_size;

    if (outsize < inpsize)
        printf ("file size reduced from %.3g Mb to %.3g Mb (%d percent)\n",
            (double)inpsize / 1048576.0, (double)outsize / 1048576.0,
            (int)(100.0 * (double)(inpsize - outsize) / (double)inpsize));
    else
        printf ("no reduction in file size\n");
    return 0;
}

