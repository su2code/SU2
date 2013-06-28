#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

void int_to_a (int num, char *s, int len)
{
    int i, j, k;

    for (i = 1; i * 10 <= num; i *= 10)
        ;
    k = 0;
    while (i > 0 && k < len - 1) {
        j = num / i;
        num %= i;
        i /= 10;
        s[k] = j + 48;
        k++;
    }
    s[k] = 0;
}

int main (int argc, char **argv)
{
    char *dbname = "dbtest.cgns";
    char buf[30];
    int i, numzones, numvalues;
    cgsize_t isize[3][1];
    int index_file, index_base, index_zone, index_coord;
    float *values;
    double start_time, end_time;

    i = 1;
    if (argc > 1 && strchr ("-ah", argv[1][0]) != NULL) {
        int type = argv[1][0];
        if (type == '-') type = argv[1][1];
        if (type == 'a')
            type = cg_set_file_type(CG_FILE_ADF);
        else if (type == 'h')
            type = cg_set_file_type(CG_FILE_HDF5);
        else {
            fputs ("unknown file type\n", stderr);
            exit (1);
        }
        if (type) {
            fprintf (stderr, "cg_set_file_type:%s\n", cg_get_error());
            exit (1);
        }
        i++;
    }
    if (i > argc - 2) {
        fputs ("usage: dbtest [a|h] numzones numvalues [CGNSfile]\n", stderr);
        exit (1);
    }
    numzones = atoi(argv[i++]);
    numvalues = atoi(argv[i++]);
    if (i < argc)
        dbname = argv[i];

    values = (float *) malloc (numvalues * sizeof(float));
    if (values == NULL) {
        perror ("malloc");
        exit (-1);
    }
    for (i = 0; i < numvalues; i++)
        values[i] = (float)i;

    start_time = elapsed_time();

    cg_open (dbname, CG_MODE_WRITE, &index_file);
    cg_base_write (index_file, "Base", 1, 1, &index_base);

    isize[0][0] = numvalues;
    isize[1][0] = isize[0][0] - 1;
    isize[2][0] = 0;

    for (i = 0; i < numzones; i++) {
        int_to_a (i, buf, sizeof(buf));
        if (cg_zone_write (index_file, index_base, buf, *isize,
                CGNS_ENUMV(Structured), &index_zone) ||
            cg_coord_write (index_file, index_base, index_zone,
                CGNS_ENUMV(RealSingle), "CoordinateX", values, &index_coord))
            cg_error_exit();
    }

    cg_close (index_file);
    end_time = elapsed_time();

#if 0
    printf ("numzones:%d numvalues:%d time:%f [s] size:%f [Mb]\n",
        numzones, numvalues, end_time - start_time,
        file_size (dbname));
#else
    printf ("%d %d %g %g\n", numzones, numvalues,
        end_time - start_time, file_size (dbname));
#endif

    return 0;
}

