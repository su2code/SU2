#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgns_io.h"
#include "getargs.h"

#ifndef CGNSTYPES_H
# define cgsize_t int
#endif

static int nocase = 0;
static int nospace = 0;
static int quiet = 0;
static int follow_links = 0;
static int recurse = 0;
static int node_data = 0;
static double tol = 0.0;

static int cgio1, cgio2;

static char options[] = "ciqrfdt:";
static char *usgmsg[] = {
    "usage  : cgnsdiff [options] CGNSfile1 [dataset1] CGNSfile2 [dataset2]",
    "options:",
    "   -c      : case insensitive names",
    "   -i      : ignore white space in names",
    "   -q      : print only if they differ",
    "   -f      : follow links",
    "   -r      : recurse (used only when dataset given)",
    "   -d      : compare node data also",
    "   -t<tol> : tolerance for comparing floats/doubles (default 0)",
    NULL
};

static void err_exit (char *msg, char *name)
{
    char errmsg[128];

    fflush (stdout);
    if (cgio_error_message (errmsg)) {
        if (msg != NULL)
            fprintf (stderr, "%s:", msg);
        if (name != NULL)
            fprintf (stderr, "%s:", name);
        fprintf (stderr, "%s\n", errmsg);
    }
    else {
        if (msg != NULL)
            fprintf (stderr, "%s", msg);
        if (name != NULL) {
            if (msg != NULL) putc (':', stderr);
            fprintf (stderr, "%s", name);
        }
        putc ('\n', stderr);
    }
    cgio_cleanup ();
    exit (1);
}

static size_t data_size (char *type, int ndim, cgsize_t *dims, int *size)
{
    int n;
    size_t bytes;

    *size = 0;
    if (ndim < 1) return 0;
    if (0 == strcmp (type, "C1") ||
        0 == strcmp (type, "B1"))
        bytes = sizeof(char);
    else if (0 == strcmp (type, "I4") ||
             0 == strcmp (type, "U4"))
        bytes = sizeof(int);
    else if (0 == strcmp (type, "I8") ||
             0 == strcmp (type, "U8"))
        bytes = sizeof(cglong_t);
    else if (0 == strcmp (type, "R4")) {
        *size = 4;
        bytes = sizeof(float);
    }
    else if (0 == strcmp (type, "R8")) {
        *size = 8;
        bytes = sizeof(double);
    }
    else if (0 == strcmp (type, "X4")) {
        *size = 4;
        bytes = 2 * sizeof(float);
    }
    else if (0 == strcmp (type, "X8")) {
        *size = 8;
        bytes = 2 * sizeof(double);
    }
    else
        return 0;

    for (n = 0; n < ndim; n++)
        bytes *= (size_t)dims[n];
    return bytes;
}

static int compare_bytes (size_t cnt, unsigned char *d1, unsigned char *d2)
{
    size_t n;

    for (n = 0; n < cnt; n++) {
        if (d1[n] != d2[n]) return 1;
    }
    return 0;
}

static int compare_floats (size_t cnt, float *d1, float *d2)
{
    size_t n;

    for (n = 0; n < cnt; n++) {
        if (fabs(d1[n] - d2[n]) > tol) return 1;
    }
    return 0;
}

static int compare_doubles (size_t cnt, double *d1, double *d2)
{
    size_t n;

    for (n = 0; n < cnt; n++) {
        if (fabs(d1[n] - d2[n]) > tol) return 1;
    }
    return 0;
}

static void compare_data (char *name1, double id1, char *name2, double id2)
{
    int n, err;
    size_t bytes;
    char label1[CGIO_MAX_NAME_LENGTH+1];
    char type1[CGIO_MAX_NAME_LENGTH+1];
    int ndim1, ndim2;
    cgsize_t dims1[CGIO_MAX_DIMENSIONS];
    cgsize_t dims2[CGIO_MAX_DIMENSIONS];
    char label2[CGIO_MAX_NAME_LENGTH+1];
    char type2[CGIO_MAX_NAME_LENGTH+1];
    void *data1, *data2;

    /* compare labels */

    if (cgio_get_label (cgio1, id1, label1))
        err_exit (name1, "cgio_get_label");
    if (cgio_get_label (cgio2, id2, label2))
        err_exit (name2, "cgio_get_label");
    if (strcmp (label1, label2)) {
        printf ("%s <> %s : labels differ\n", name1, name2);
        return;
    }

    /* compare data types */

    if (cgio_get_data_type (cgio1, id1, type1))
        err_exit (name1, "cgio_get_data_type");
    if (cgio_get_data_type (cgio2, id2, type2))
        err_exit (name2, "cgio_get_data_type");
    if (strcmp (type1, type2)) {
        printf ("%s <> %s : data types differ\n", name1, name2);
        return;
    }

    /* compare number of dimensions */

    if (cgio_get_dimensions (cgio1, id1, &ndim1, dims1))
        err_exit (name1, "cgio_get_dimensions");
    if (cgio_get_dimensions (cgio2, id2, &ndim2, dims2))
        err_exit (name2, "cgio_get_dimensions");
    if (ndim1 != ndim2) {
        printf ("%s <> %s : number of dimensions differ\n", name1, name2);
        return;
    }

    /* compare dimensions */

    if (ndim1 > 0) {
        for (n = 0; n < ndim1; n++) {
            if (dims1[n] != dims2[n]) {
                printf ("%s <> %s : dimensions differ\n", name1, name2);
                return;
            }
        }
    }

    if (!node_data || ndim1 <= 0) return;

    /* compare data */

    bytes = data_size (type1, ndim1, dims1, &n);
    if (bytes > 0) {
        data1 = malloc (bytes);
        if (data1 == NULL) {
            fprintf (stderr, "%s:malloc failed for node data\n", name1);
            exit (1);
        }
        data2 = malloc (bytes);
        if (data2 == NULL) {
            fprintf (stderr, "%s:malloc failed for node data\n", name2);
            exit (1);
        }
        if (cgio_read_all_data (cgio1, id1, data1))
            err_exit (name1, "cgio_read_all_data");
        if (cgio_read_all_data (cgio2, id2, data2))
            err_exit (name2, "cgio_read_all_data");
        if (tol > 0.0 && n) {
            if (n == 4)
                err = compare_floats (bytes >> 2, data1, data2);
            else
                err = compare_doubles (bytes >> 3, data1, data2);
        }
        else {
            err = compare_bytes (bytes, (unsigned char *)data1,
                  (unsigned char *)data2);
        }
        free (data1);
        free (data2);
        if (err)
            printf ("%s <> %s : data values differ\n", name1, name2);
    }
}

static void copy_name (char *name, char *newname)
{
    int n1, n2;

    if (nospace) {
        if (nocase) {
            for (n1 = 0, n2 = 0; name[n1]; n1++) {
                if (!isspace (name[n1]))
                    newname[n2++] = tolower (name[n1]);
            }
        }
        else {
            for (n1 = 0, n2 = 0; name[n1]; n1++) {
                if (!isspace (name[n1]))
                    newname[n2++] = name[n1];
            }
        }
        newname[n2] = 0;
    }
    else if (nocase) {
        for (n1 = 0; name[n1]; n1++)
            newname[n1] = tolower (name[n1]);
        newname[n1] = 0;
    }
    else
        strcpy (newname, name);
}

static int sort_children (const void *v1, const void *v2)
{
    char p1[33], p2[33];

    copy_name ((char *)v1, p1);
    copy_name ((char *)v2, p2);
    return strcmp (p1, p2);
}

static int find_name (char *name, int nlist, char *namelist)
{
    int cmp, mid, lo = 0, hi = nlist - 1;
    char p1[33], p2[33];

    copy_name (name, p1);
    copy_name (namelist, p2);
    if (0 == strcmp (p1, p2)) return 0;
    copy_name (&namelist[33*hi], p2);
    if (0 == strcmp (p1, p2)) return hi;

    while (lo <= hi) {
        mid = (lo + hi) >> 1;
        copy_name (&namelist[33*mid], p2);
        cmp = strcmp (p1, p2);
        if (0 == cmp) return mid;
        if (cmp > 0)
            lo = mid + 1;
        else
            hi = mid - 1;
    }
    return -1;
}

static void compare_nodes (char *name1, double id1, char *name2, double id2)
{
    int n1, n2, nc1, nc2, nret;
    char *p, *children1 = NULL, *children2 = NULL;
    char path1[1024], path2[1024];
    double cid1, cid2;

    compare_data (name1, id1, name2, id2);
    if (!recurse) return;
    if (!follow_links) {
        if (cgio_is_link (cgio1, id1, &nret))
            err_exit (name1, "cgio_is_link");
        if (nret > 0) return;
        if (cgio_is_link (cgio2, id2, &nret))
            err_exit (name2, "cgio_is_link");
        if (nret > 0) return;
    }

    if (cgio_number_children (cgio1, id1, &nc1))
        err_exit (name1, "cgio_number_children");
    if (nc1) {
        children1 = (char *) malloc (33 * nc1);
        if (children1 == NULL) {
            fprintf (stderr, "%s:malloc failed for children names\n", name1);
            exit (1);
        }
        if (cgio_children_names (cgio1, id1, 1, nc1, 33,
                &nret, children1))
            err_exit (name1, "cgio_children_names");
        if (nc1 > 1)
            qsort (children1, nc1, 33, sort_children);
    }

    if (cgio_number_children (cgio2, id2, &nc2))
        err_exit (name2, "cgio_number_children");
    if (nc2) {
        children2 = (char *) malloc (33 * nc2);
        if (children2 == NULL) {
            fprintf (stderr, "%s:malloc failed for children names\n", name2);
            exit (1);
        }
        if (cgio_children_names (cgio2, id2, 1, nc2, 33,
                &nret, children2))
            err_exit (name2, "cgio_children_names");
        if (nc2 > 1)
            qsort (children2, nc2, 33, sort_children);
    }

    if (0 == strcmp (name1, "/")) name1 = "";
    if (0 == strcmp (name2, "/")) name2 = "";
    if (nc1 == 0) {
        if (nc2 == 0) return;
        for (n2 = 0; n2 < nc2; n2++)
            printf ("> %s/%s\n", name2, &children2[33*n2]);
        free (children2);
        return;
    }
    if (nc2 == 0) {
        for (n1 = 0; n1 < nc1; n1++)
            printf ("< %s/%s\n", name1, &children1[33*n1]);
        free (children1);
        return;
    }

    for (n1 = 0, n2 = 0; n1 < nc1; n1++) {
        p = &children1[33*n1];
        nret = find_name (p, nc2, children2);
        if (nret < 0) {
            printf ("< %s/%s\n", name1, p);
            continue;
        }
        while (n2 < nret) {
            printf ("> %s/%s\n", name2, &children2[33*n2]);
            n2++;
        }
        if (cgio_get_node_id (cgio1, id1, p, &cid1))
            err_exit (name1, "cgio_get_node_id");
        sprintf (path1, "%s/%s", name1, p);
        p = &children2[33*n2];
        if (cgio_get_node_id (cgio2, id2, p, &cid2))
            err_exit (name2, "cgio_get_node_id");
        sprintf (path2, "%s/%s", name2, p);
        compare_nodes (path1, cid1, path2, cid2);
        n2++;
    }
    while (n2 < nc2-1) {
        printf ("> %s/%s\n", name2, &children2[33*n2]);
        n2++;
    }
    free (children1);
    free (children2);
}

int main (int argc, char *argv[])
{
    double root1, root2;
    double node1, node2;
    int n;
    char *file1, *file2;
    char *ds1 = NULL;
    char *ds2 = NULL;

    if (argc < 3)
        print_usage (usgmsg, NULL);
    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'c':
                nocase = 1;
                break;
            case 'i':
                nospace = 1;
                break;
            case 'q':
                quiet = 1;
                break;
            case 'f':
                follow_links = 1;
                break;
            case 'r':
                recurse = 1;
                break;
            case 'd':
                node_data = 1;
                break;
            case 't':
                tol = atof (argarg);
                break;
        }
    }

    if (argind > argc - 2)
        print_usage (usgmsg, "CGNSfile1 and/or CGNSfile2 not given");

    file1 = argv[argind++];
    if (cgio_open_file (file1, 'r', CGIO_FILE_NONE, &cgio1))
        err_exit (file1, "cgio_open_file");
    if (cgio_get_root_id (cgio1, &root1))
        err_exit (file1, "cgio_get_root_id");
    if (argind < argc - 1)
        ds1 = argv[argind++];

    file2 = argv[argind++];
    if (cgio_open_file (file2, 'r', CGIO_FILE_NONE, &cgio2))
        err_exit (file2, "cgio_open_file");
    if (cgio_get_root_id (cgio2, &root2))
        err_exit (file2, "cgio_get_root_id");
    if (argind < argc)
        ds2 = argv[argind++];

    if (ds1 == NULL) {
        recurse = 1;
        compare_nodes ("/", root1, "/", root2);
    }
    else {
        if (cgio_get_node_id (cgio1, root1, ds1, &node1))
            err_exit (ds1, "cgio_get_node_id");
        if (ds2 == NULL) ds2 = ds1;
        if (cgio_get_node_id (cgio2, root2, ds2, &node2))
            err_exit (ds2, "cgio_get_node_id");
        compare_nodes (ds1, node1, ds2, node2);
    }

    cgio_close_file (cgio1);
    cgio_close_file (cgio2);
    return 0;
}

