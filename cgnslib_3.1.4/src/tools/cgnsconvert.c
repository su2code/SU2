#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
# include <io.h>
# define unlink _unlink
#else
# include <unistd.h>
#endif

#include "cgns_io.h"
#include "getargs.h"

static char options[] = "ahfl";
static char *usgmsg[] = {
    "usage  : cgnsconvert [options] InputFile [OutputFile]",
    "options:",
    "   -a : write ADF file",
    "   -h : write HDF5 file",
    "   -f : force output if input format is same as output",
    "   -l : expand links in ouput file",
    NULL
};

int main (int argc, char **argv)
{
    char *inpfile, *outfile;
    char tempfile[1024];
    int n, inptype, outtype = CGIO_FILE_NONE;
    int inpcg, outcg, links = 0;
    int force = 0;
    struct stat inpst, outst;
    time_t ts, te;
    static char *FileType[] = {"NONE", "ADF", "HDF5"};

    if (argc < 2)
        print_usage (usgmsg, NULL);
    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'a':
                outtype = CGIO_FILE_ADF;
                break;
            case 'h':
                outtype = CGIO_FILE_HDF5;
                break;
            case 'f':
                force = 1;
                break;
            case 'l':
                links = 1;
                break;
        }
    }

    if (argind == argc)
        print_usage (usgmsg, "InputFile not given");
    if (outtype != CGIO_FILE_NONE && cgio_is_supported(outtype)) {
        fprintf(stderr, "output type %s not supported\n",
            FileType[outtype]);
        exit(1);
    }

    inpfile = argv[argind++];
    if (argind < argc)
        outfile = argv[argind];
    else
        outfile = inpfile;
    sprintf(tempfile, "%s.temp", outfile);
    unlink(tempfile);

    if (stat (inpfile, &inpst)) {
        fprintf (stderr, "can't stat %s\n", inpfile);
        exit (1);
    }

    if (cgio_open_file (inpfile, 'r', CGIO_FILE_NONE, &inpcg))
        cgio_error_exit("cgio_open_file");
    if (cgio_get_file_type (inpcg, &inptype))
        cgio_error_exit("cgio_get_file_type");

    if (outtype == CGIO_FILE_NONE) {
        if (inptype == CGIO_FILE_ADF)
            outtype = CGIO_FILE_HDF5;
        else
            outtype = CGIO_FILE_ADF;
    }
    if (!force && outtype == inptype) {
        cgio_close_file(inpcg);
        fputs("input and output formats the same: use -f to force write\n",
            stderr);
        return 1;
    }

    printf("converting %s file %s to %s file %s\n",
        FileType[inptype], inpfile, FileType[outtype], outfile);
    if (links)
        printf ("links will be included in output file\n");
    fflush(stdout);

    ts = time (NULL);
    if (cgio_open_file (tempfile, 'w', outtype, &outcg))
        cgio_error_exit("cgio_open_file");
    if (cgio_copy_file (inpcg, outcg, links))
        cgio_error_exit("cgio_copy_file");
    if (cgio_close_file (inpcg) || cgio_close_file (outcg))
        cgio_error_exit("cgio_close_file");
    te = time (NULL);

    unlink (outfile);
    if (rename (tempfile, outfile)) {
        fprintf (stderr, "rename %s -> %s failed", tempfile, outfile);
        exit (1);
    }

    if (stat (outfile, &outst)) {
        fprintf (stderr, "can't stat %s\n", outfile);
        exit (1);
    }

    printf ("%-4s input  file size  = %ld bytes\n",
        FileType[inptype], (long)inpst.st_size);
    printf ("%-4s output file size  = %ld bytes\n",
        FileType[outtype], (long)outst.st_size);
    printf ("conversion time = %d secs\n", (int)(te - ts));
    return 0;
}

