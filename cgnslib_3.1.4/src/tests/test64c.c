#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef _WIN32
# include <io.h>
# define unlink _unlink
#else
# include <unistd.h>
#endif
#include "utils.h"

#if CGNS_VERSION < 3100
# error need CGNS version 3.1 or greater
#endif
#if CG_SIZEOF_SIZE != 64
# error need 64-bit build
#endif

#define MAXMEM 16000000000

void print_usage()
{
    printf("usage: test64c [options]\n");
    printf("options:\n");
    printf("  --                  - this message\n");
    printf("  [-]he[lp]           - this message\n");
    printf("  [-]a[df]            - use ADF\n");
    printf("  [-]h[df5]           - use HDF5\n");
    printf("  [-]t[est]           - read and test data\n");
    printf("  [-]m bytes[k|m|g]   - set max memory (default=16G)\n");
    printf("  [-]c ncoords[k|m|b] - set num coords (default=maxmem/4) \n");
    printf("  [-]e nelems[k|m|b]] - set num elements (default=maxmem/32)\n");
}

char *format_size(cgsize_t size)
{
    static char buff[64];

    if (size >= 1000000000)
        sprintf(buff, "%gb", (double)size / 1000000000);
    else if (size >= 1000000)
        sprintf(buff, "%gm", (double)size / 1000000);
    else if (size >= 1000)
        sprintf(buff, "%gk", (double)size / 1000);
    else
        sprintf(buff, "%ld", (long)size);
    return buff;
}

char *format_time(int secs)
{
    int h, m, s;
    static char buff[64];

    h = secs / 3600;
    s = secs - (h * 3600);
    m = s / 60;
    s -= (m * 60);
    if (h)
        sprintf(buff, "%dh%dm%ds", h, m, s);
    else if (m)
        sprintf(buff, "%dm%ds", m, s);
    else
        sprintf(buff, "%ds", s);
    return buff;
}

int main (int argc, char **argv)
{
    char *dbname = "test64.cgns";
    int i, j, k, celldim = 3, phydim = 3;
    cgsize_t n, sizes[3];
    cgsize_t ns, ne, nn, count;
    cgsize_t maxmem = 0;
    cgsize_t nnodes = -1;
    cgsize_t nelems = -1;
    int index_file, index_base, index_zone;
    int index_coord, index_sect, test_data = 0;
    void *data;
    float *nodes;
    cgsize_t *elems;
    double start_time, end_time, dsize;
    char bname[33], zname[33], cname[33];
    cgsize_t nerr, errs = 0;
    CGNS_ENUMT(ZoneType_t) ztype;
    CGNS_ENUMT(ElementType_t) etype;
    time_t tstart, tend;

    for (i = 1; i < argc; i++) {
        j = 0;
        if (argv[i][j] == '-') j++;
        switch (argv[i][j]) {
            case 'a':
            case 'A':
                if (cg_set_file_type(CG_FILE_ADF)) {
                    fprintf(stderr, "cg_set_file_type:%s\n",
                        cg_get_error());
                    exit(1);
                }
                break;
            case 'h':
            case 'H':
                if (argv[i][j+1] == 'e' || argv[i][j+1] == 'E') {
                    print_usage();
                    return 0;
                }
                if (cg_set_file_type(CG_FILE_HDF5)) {
                    fprintf(stderr, "cg_set_file_type:%s\n",
                        cg_get_error());
                    exit(1);
                }
                break;
            case 't':
                test_data = 1;
                break;
            case 'm':
            case 'c':
            case 'e':
                k = argv[i][j];
                if (argv[i][++j] == 0) {
                    if (++i >= argc) {
                        fprintf(stderr, "missing argument to %c\n", k);
                        exit(1);
                    }
                    j = 0;
                }
                dsize = atof(&argv[i][j]);
                if (dsize < 0.0) {
                    fprintf(stderr, "bad size given to %c\n", k);
                    exit(1);
                }
                if (strchr(&argv[i][j], 'k') ||
                    strchr(&argv[i][j], 'K')) {
                    count = (cgsize_t)(1000.0 * dsize);
                }
                else if (strchr(&argv[i][j], 'm') ||
                         strchr(&argv[i][j], 'M')) {
                    count = (cgsize_t)(1000000.0 * dsize);
                }
                else if (strchr(&argv[i][j], 'b') ||
                         strchr(&argv[i][j], 'B') ||
                         strchr(&argv[i][j], 'g') ||
                         strchr(&argv[i][j], 'G')) {
                    count = (cgsize_t)(1000000000.0 * dsize);
                }
                else {
                    count = (cgsize_t)dsize;
                }
                if (k == 'm')
                    maxmem = count;
                else if (k == 'c')
                    nnodes = count;
                else
                    nelems = count;
                break;
            default:
                fprintf(stderr, "unknown option %s\n", argv[i]);
            case '?':
            case '-':
                print_usage();
                exit(argv[i][j] != '?' && argv[i][j] != '-');
        }
    }

    if (maxmem == 0) maxmem = MAXMEM;
    if (nnodes < 0) nnodes = maxmem / sizeof(float);
    if (nelems < 0) nelems = maxmem / (4 * sizeof(cgsize_t));

    printf("sizeof int      = %d\n", (int)sizeof(int));
    printf("sizeof long     = %d\n", (int)sizeof(long));
    printf("sizeof size_t   = %d\n", (int)sizeof(size_t));
    printf("sizeof cgsize_t = %d\n", (int)sizeof(cgsize_t));
    printf("sizeof cglong_t = %d\n", (int)sizeof(cglong_t));
    printf("sizeof void *   = %d\n", (int)sizeof(void *));
    printf("memory size     = %s\n", format_size(maxmem));
    printf("num coords      = %s\n", format_size(nnodes));
    printf("num elements    = %s\n", format_size(nelems));

    data = malloc ((size_t)maxmem);
    if (data == NULL) {
        printf("failed to allocate %s bytes of memory\n",
            format_size(maxmem));
        exit(1);
    }

    sizes[0] = nnodes ? nnodes : 1;
    sizes[1] = nelems ? nelems : 1;
    sizes[2] = 0;

    tstart = time(NULL);
    start_time = elapsed_time();

    unlink(dbname);
    printf("creating %s...", dbname);
    fflush(stdout);
    if (cg_open (dbname, CG_MODE_WRITE, &index_file) ||
        cg_base_write (index_file, "Base", celldim, phydim, &index_base) ||
        cg_zone_write (index_file, index_base, "Zone", sizes,
            CGNS_ENUMV(Unstructured), &index_zone))
        cg_error_exit();
    puts("finished");
    if (cg_get_file_type(index_file, &i))
        cg_error_exit();
    printf("using %s file type\n", i == CG_FILE_ADF ? "ADF" :
        (i == CG_FILE_HDF5 ? "HDF5" : "unknown"));
    fflush(stdout);

    if (nnodes > 0) {
        count = maxmem / sizeof(float);
        if (count > nnodes) count = nnodes;
        nodes = (float *)data;

        if (count == nnodes) {
            printf("writing %s coordinates...", format_size(nnodes));
            fflush(stdout);
            for (n = 0; n < nnodes; n++)
                nodes[n] = (float)(n + 1);
            if (cg_coord_write (index_file, index_base, index_zone,
                    CGNS_ENUMV(RealSingle), "Coordinates", nodes,
                    &index_coord))
                cg_error_exit();
            puts("finished");
        }
        else {
            ns = 1;
            while (ns < nnodes) {
                ne = ns + count - 1;
                if (ne > nnodes) ne = nnodes;
                printf("writing coordinates %s", format_size(ns));
                printf(" -> %s...", format_size(ne));
                fflush(stdout);
                for (n = 0, nn = ns; nn <= ne; nn++, n++)
                    nodes[n] = (float)nn;
                if (cg_coord_partial_write (index_file, index_base,
                        index_zone, CGNS_ENUMV(RealSingle), "Coordinates",
                        &ns, &ne, nodes, &index_coord))
                    cg_error_exit();
                ns = ne + 1;
                puts("finished");
            }
        }
    }

    if (nelems > 0) {
        count = maxmem / (4 * sizeof(cgsize_t));
        if (count > nelems) count = nelems;
        elems = (cgsize_t *)data;

        if (count == nelems) {
            printf("writing %s tetrahedra...", format_size(nelems));
            fflush(stdout);
            for (n = 0, nn = 1; nn <= nelems; nn++) {
                for (i = 0; i < 4; i++) {
                    elems[n++] = nn;
                }
            }
            if (cg_section_write(index_file, index_base, index_zone,
                    "Elements",  CGNS_ENUMV(TETRA_4), 1, nelems,
                    0, elems, &index_sect))
                cg_error_exit();
            puts("finished");
        }
        else {
            if (cg_section_partial_write(index_file, index_base,
                    index_zone, "Elements",  CGNS_ENUMV(TETRA_4),
                    1, nelems, 0, &index_sect))
                cg_error_exit();
            ns = 1;
            while (ns < nelems) {
                ne = ns + count - 1;
                if (ne > nelems) ne = nelems;
                printf("writing tetrahedra %s", format_size(ns));
                printf(" -> %s...", format_size(ne));
                fflush(stdout);
                for (n = 0, nn = ns; nn <= ne; nn++) {
                    for (i = 0; i < 4; i++) {
                        elems[n++] = nn;
                    }
                }
                if (cg_elements_partial_write(index_file, index_base,
                        index_zone, index_sect, ns, ne, elems))
                    cg_error_exit();
                ns = ne + 1;
                puts("finished");
            }
        }
    }

    printf("closing file...\n");
    cg_close (index_file);
    tend = time(NULL);
    end_time = elapsed_time();
    count = 1000000 * (cgsize_t)(0.5 + file_size(dbname));
    printf("elapsed time=%s, ", format_time((int)(tend - tstart)));
    printf("cpu time=%s, file size=%s\n",
        format_time((int)(end_time - start_time)),
        format_size(count));

    if (!test_data) return 0;

    tstart = time(NULL);
    start_time = elapsed_time();

    printf("opening in read only mode...\n");
    index_base = index_zone = 1;
    if (cg_open (dbname, CG_MODE_READ, &index_file) ||
        cg_base_read (index_file, index_base, bname, &celldim, &phydim) ||
        cg_zone_type (index_file, index_base, index_zone, &ztype) ||
        cg_zone_read (index_file, index_base, index_zone, zname, sizes))
        cg_error_exit();

    printf("reading and checking data...\n");
    if (strcmp (bname, "Base")) {
        errs++;
        printf("bad base name=%s\n", bname);
    }
    if (celldim != 3 || phydim != 3) {
        errs++;
        printf("bad cell dim=%d or phy dim=%d\n", celldim, phydim);
    }
    if (strcmp (zname, "Zone")) {
        errs++;
        printf("bad zone name=%s\n", zname);
    }
    if (ztype != CGNS_ENUMV(Unstructured)) {
        errs++;
        printf("bad zone type=%d\n", ztype);
    }
    if (nnodes && sizes[0] != nnodes) {
        errs++;
        printf("bad num points=%ld\n", (long)sizes[0]);
    }
    if (nelems && sizes[1] != nelems) {
        errs++;
        printf("bad num elements=%ld\n", (long)sizes[1]);
    }
    if (sizes[2] != 0) {
        errs++;
        printf("bad num bndry=%ld\n", sizes[2]);
    }

    if (nnodes > 0) {
        count = maxmem / sizeof(float);
        if (count > nnodes) count = nnodes;
        nodes = (float *)data;

        ns = 1;
        while (ns < nnodes) {
            nerr = 0;
            ne = ns + count - 1;
            if (ne > nnodes) ne = nnodes;
            printf("reading coordinates %s ->", format_size(ns));
            printf(" %s...", format_size(ne));
            fflush(stdout);
            if (cg_coord_read(index_file, index_base, index_zone,
                    "Coordinates", CGNS_ENUMV(RealSingle),
                    &ns, &ne, nodes))
                cg_error_exit();
            for (n = 0, nn = ns; nn <= ne; nn++, n++) {
                if (nodes[n] != (float)nn) nerr++;
            }
            printf("%s errors\n", nerr ? format_size(nerr) : "no");
            fflush(stdout);
            errs += nerr;
            ns = ne + 1;
        }
    }

    if (nelems > 0) {
        nerr = 0;
        if (cg_section_read(index_file, index_base, index_zone,
                index_sect, cname, &etype, &ns, &ne, &i, &j))
            cg_error_exit();
        if (strcmp (cname, "Elements")) {
            errs++;
            printf("bad section name=%s\n", cname);
        }
        if (ns != 1 || ne != nelems) {
            nerr++;
            errs++;
            printf("bad start=%ld or end=%ld element range\n",
                (long)ns, (long)ne);
        }
        if (etype != CGNS_ENUMV(TETRA_4)) {
            nerr++;
            errs++;
            printf("bad element type=%d\n", etype);
        }

        if (nerr == 0) {
            count = maxmem / (4 * sizeof(cgsize_t));
            if (count > nelems) count = nelems;
            elems = (cgsize_t *)data;

            if (count == nelems) {
                nerr = 0;
                printf("reading tetrahedra 1 -> %s...",
                    format_size(nelems));
                fflush(stdout);
                if (cg_elements_read(index_file, index_base,
                        index_zone, index_sect, elems, 0))
                    cg_error_exit();
                for (n = 0, nn = 1; nn <= nelems; nn++) {
                    for (i = 0; i < 4; i++) {
                        if (elems[n++] != nn) nerr++;
                    }
                }
                printf("%s errors\n", nerr ? format_size(nerr) : "no");
                fflush(stdout);
                errs += nerr;
            }
            else {
                ns = 1;
                while (ns < nelems) {
                    nerr = 0;
                    ne = ns + count - 1;
                    if (ne > nelems) ne = nelems;
                    printf("reading tetrahedra %s", format_size(ns));
                    printf(" -> %s...", format_size(ne));
                    fflush(stdout);
                    if (cg_elements_partial_read(index_file, index_base,
                            index_zone, index_sect, ns, ne, elems, 0))
                        cg_error_exit();
                    for (n = 0, nn = ns; nn <= ne; nn++) {
                        for (i = 0; i < 4; i++) {
                            if (elems[n++] != nn) nerr++;
                        }
                    }
                    printf("%s errors\n", nerr ? format_size(nerr) : "no");
                    fflush(stdout);
                    errs += nerr;
                    ns = ne + 1;
                }
            }
        }
    }

    cg_close (index_file);
    tend = time(NULL);
    end_time = elapsed_time();
    free(data);
    printf ("%s errors, elapsed time=%s, ",
        errs ? format_size(errs) : "no",
        format_time((int)(tend - tstart)));
    printf("cpu time=%s\n", format_time((int)(end_time - start_time)));

    return 0;
}

