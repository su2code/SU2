#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "cgns_io.h"
#include "getargs.h"

#ifndef CGNSTYPES_H
# define cgsize_t int
#endif

#define MAX_LEADER 1024
#define INDENT     2

static char leader[MAX_LEADER+1];
static int leader_len;
static int indent = INDENT;
static int follow_links = 0;
static int out_flags = 0;

static char options[] = "bi:faltds";
static char *usgmsg[] = {
    "usage  : cgnslist [options] CGNSfile [node]",
    "options:",
    "   -b      = brief - file summary only",
    "   -i<cnt> = set indent level (default 2)",
    "   -f      = follow links",
    "   -l      = print node label",
    "   -t      = print node data type",
    "   -d      = print node dimensions",
    "   -s      = print node size in bytes",
    "   -a      = print all -ltds",
    NULL
};

static void print_node (int cgio, double node_id)
{
    int n, ndim;
    char label[CGIO_MAX_NAME_LENGTH+1];
    char type[CGIO_MAX_NAME_LENGTH+1];
    cgsize_t bytes, dims[CGIO_MAX_DIMENSIONS];

    if ((out_flags & 1) != 0) {
        if (cgio_get_label (cgio, node_id, label))
            cgio_error_exit ("cgio_get_label");
        printf (" %s", label);
    }
    if ((out_flags & 10) != 0) {
        if (cgio_get_data_type (cgio, node_id, type))
            cgio_error_exit ("cgio_get_data_type");
        if ((out_flags & 2) != 0)
            printf (" %s", type);
    }
    if ((out_flags & 12) != 0) {
        if (cgio_get_dimensions (cgio, node_id, &ndim, dims))
            cgio_error_exit ("cgio_get_data_type");
        if ((out_flags & 4) != 0) {
            printf (" (");
            if (ndim > 0) {
                printf ("%ld", (long)dims[0]);
                for (n = 1; n < ndim; n++)
                    printf (",%ld", (long)dims[n]);
            }
            putchar (')');
        }
        if ((out_flags & 8) != 0) {
            if (ndim < 1 || NULL != strchr ("LlMm", type[0]))
                bytes = 0;
            else if (NULL != strchr ("CcBb", type[0]))
                bytes = 1;
            else if (type[0] == 'X' || type[0] == 'x')
                bytes = type[1] == '8' ? 16 : 8;
            else
                bytes = type[1] == '8' ? 8 : 4;
            for (n = 0; n < ndim; n++)
                bytes *= dims[n];
            printf (" %ld", (long)bytes);
        }
    }
}

static void print_children (int cgio, double parent_id)
{
    int nc, nchildren, len_ret;
    char *p = leader + leader_len;
    char name[CGIO_MAX_NAME_LENGTH+1];
    char name_in_file[CGIO_MAX_LINK_LENGTH+1];
    char file_name[CGIO_MAX_FILE_LENGTH+1];
    double child_id;

    if (cgio_number_children (cgio, parent_id, &nchildren))
        cgio_error_exit ("cgio_number_children");
    if (!nchildren) return;

    if (leader_len + indent > MAX_LEADER) {
        fprintf (stderr, "nesting is too deep\n");
        exit (1);
    }
    leader_len += indent;
    for (nc = 0; nc < indent; nc++)
        p[nc] = ' ';
    p[indent] = 0;

    for (nc = 1; nc <= nchildren; nc++) {
        if (cgio_children_ids (cgio, parent_id, nc, 1, &len_ret,
               &child_id))
            cgio_error_exit ("cgio_children_ids");
        if (cgio_get_name (cgio, child_id, name))
            cgio_error_exit ("cgio_get_name");
        if (cgio_is_link (cgio, child_id, &len_ret))
            cgio_error_exit ("cgio_is_link");

        *p = 0;
        if (len_ret > 0) {
            if (cgio_get_link (cgio, child_id, file_name, name_in_file))
                cgio_error_exit ("cgio_get_link");
            if (*file_name)
                printf ("%s+-%s  -> %s @ %s\n", leader, name,
                    name_in_file, file_name);
            else
                printf ("%s+-%s  -> %s\n", leader, name, name_in_file);
        }
        else if (out_flags) {
            printf ("%s+-%s  --", leader, name);
            print_node (cgio, child_id);
            putchar ('\n');
        }
        else
            printf ("%s+-%s\n", leader, name);

        if (follow_links || len_ret <= 0) {
            *p = (char)(nc < nchildren ? '|' : ' ');
            print_children (cgio, child_id);
        }
    }
    *p = 0;
    leader_len -= indent;
}

int main (int argc, char *argv[])
{
    double root_id, node_id;
    float cgns_version;
    int n = 1, cgio, file_type, brief = 0;
    char *name, rootname[CGIO_MAX_NAME_LENGTH+1];
    struct stat st;
    char version[CGIO_MAX_NAME_LENGTH+1];
    char created[CGIO_MAX_NAME_LENGTH+1];
    char modified[CGIO_MAX_NAME_LENGTH+1];
    static char *FileType[] = {"NONE", "ADF", "HDF5"};

    if (argc < 2)
        print_usage (usgmsg, NULL);
    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'b':
                brief = 1;
                break;
            case 'i':
                indent = atoi (argarg);
                if (indent < 1) {
                    fprintf (stderr, "indent must be > 0\n");
                    exit (1);
                }
                break;
            case 'f':
                follow_links = 1;
                break;
            case 'l':
                out_flags |= 1;
                break;
            case 't':
                out_flags |= 2;
                break;
            case 'd':
                out_flags |= 4;
                break;
            case 's':
                out_flags |= 8;
                break;
            case 'a':
                out_flags |= 15;
                break;
        }
    }

    if (argind == argc)
        print_usage (usgmsg, "CGNSfile not given");

    if (stat (argv[argind], &st)) {
        fprintf (stderr, "can't stat %s\n", argv[argind]);
        exit (1);
    }

    if (cgio_open_file (argv[argind], 'r', CGIO_FILE_NONE, &cgio))
        cgio_error_exit ("cgio_open_file");
    if (cgio_get_root_id (cgio, &root_id))
        cgio_error_exit ("cgio_get_root_id");

    if (brief) {
        if (cgio_get_file_type (cgio, &file_type))
            cgio_error_exit ("cgio_get_file_type");
        if (cgio_file_version (cgio, version, created, modified))
            cgio_error_exit ("cgio_file_version");
        if (0 == cgio_get_node_id (cgio, root_id,
                     "CGNSLibraryVersion",&node_id) &&
            0 == cgio_read_all_data (cgio, node_id, &cgns_version))
            printf ("CGNS version  : %4.2f\n", cgns_version);
        else
            printf ("CGNS version  : not defined\n");
        printf ("file type     : %s\n", FileType[file_type]);
        printf ("file version  : %s\n", version);
        printf ("file size     : %ld bytes\n", (long)st.st_size);
        printf ("creation date : %s\n", created);
        printf ("modified date : %s\n", modified);
        if (cgio_close_file (cgio))
            cgio_error_exit ("cgio_close_file");
        return 0;
    }

    if (++argind < argc) {
        name = argv[argind];
        if (cgio_get_node_id (cgio, root_id, name, &node_id))
            cgio_error_exit ("cgio_get_root_id");
    }
    else {
        if (cgio_get_name (cgio, root_id, rootname))
            cgio_error_exit ("cgio_get_name");
        node_id = root_id;
        name = rootname;
    }

    for (n = 0; n < indent; n++)
        leader[n] = ' ';
    leader[indent] = 0;
    leader_len = indent;

    if (out_flags) {
        printf ("%s  --", name);
        print_node (cgio, node_id);
        putchar ('\n');
    }
    else
        printf ("%s\n", name);
    print_children (cgio, node_id);

    if (cgio_close_file (cgio))
        cgio_error_exit ("cgio_close_file");
    return 0;
}

