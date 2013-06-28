/*
 * patran_to_cgns - Patran Neutral File Import
 * reads packets 1 (nodes), 2 (elements) and 21 (named groups)
 * and optionally packet 6 (loads), and writes CGNS
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "cgnsImport.h"
#include "getargs.h"

#define getline()   \
    {if (NULL == fgets (buffer, sizeof(buffer), fp)) \
        cgnsImportFatal ("premature EOF"); \
    lineno++;}

static char options[] = "ldDt:B:";

static char *usgmsg[] = {
    "usage  : patran_to_cgns [options] Patranfile CGNSfile",
    "options:",
    "   -l       = process packet 6 - distributed loads",
    "   -d       = duplicate node checking with rel tol",
    "   -D       = duplicate node checking with abs tol",
    "   -t<tol>  = duplicate node checking tolerance",
    "   -B<name> = set CGNS base name to <name>",
    NULL
};

/*---------- print_error -------------------------------------------
 * print error message and line number
 *------------------------------------------------------------------*/

static int lineno = 0;

static void print_error (char *errmsg)
{
    fprintf (stderr, "%s on line %d\n", errmsg, lineno);
}

/*---------- add_face -----------------------------------------------
 * add an element face to the region list
 *-------------------------------------------------------------------*/

static void add_face (int elemid, char *data)
{
    int n, nodes[8], nnodes;
    int elemtype, faceid;
    cgsize_t elemnodes[8], nodeid[8];
    char errmsg[81];
    static int facemap[5][7] = {
        {0, 1, 2, 3, 4, 0, 0},
        {0, 1, 2, 3, 4, 5, 0},
        {0, 4, 5, 1, 3, 2, 0},
        {0, 0, 0, 0, 0, 0, 0},
        {0, 2, 4, 1, 3, 6, 5}
    };

    /* check for node flags set */

    for (nnodes = 0, n = 0; n < 8; n++) {
        if ('1' == data[n])
            nodes[nnodes++] = n;
    }

    elemtype = cgnsImportGetElement (elemid, elemnodes);
    if (!elemtype) {
        sprintf (errmsg, "element %d not found for packet 6\n", elemid);
        cgnsImportFatal (errmsg);
    }

    /* if node flags set, use the node values */

    if (nnodes) {
        for (n = 0; n < nnodes; n++) {
            if (nodes[n] >= elemtype) {
                sprintf (errmsg,
                    "invalid node flags for element %d\n", elemid);
                cgnsImportFatal (errmsg);
            }
            nodeid[n] = elemnodes[nodes[n]];
        }
    }

    /* else get nodes from face number */

    else {
        faceid = atoi (&data[8]);
        if (faceid < 1 || faceid > 6) {
            sprintf (errmsg, "invalid faceid for element %d\n", elemid);
            cgnsImportFatal (errmsg);
        }
        faceid = facemap[elemtype-4][faceid];
        nnodes = cgnsImportGetFace (elemid, faceid, nodeid);
        if (nnodes < 0) {
            sprintf (errmsg,
                "element %d not found for packet 6\n", elemid);
            cgnsImportFatal (errmsg);
        }
        if (0 == nnodes) {
            sprintf (errmsg,
                "invalid face number for element %d\n", elemid);
            cgnsImportFatal (errmsg);
        }
    }

    cgnsImportAddReg (nnodes, nodeid);
}

/*========== main ===================================================*/

int main (int argc, char *argv[])
{
    int n, packet, nlines, nodeid;
    int nnodes, elemid;
    cgsize_t nodes[8];
    int lastid = -1, loadid;
    int do_loads = 0, do_chk = 0;
    double xyz[3];
    char *p, buffer[256], *basename = NULL;
    FILE *fp;

    if (argc < 2)
        print_usage (usgmsg, NULL);

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'l':
                do_loads = 1;
                break;
            case 'D':
            case 'd':
                do_chk = n;
                break;
            case 't':
                cgnsImportSetTol (atof (argarg));
                break;
            case 'B':
                basename = argarg;
                break;
        }
    }

    if (argind + 1 >= argc)
        print_usage (usgmsg, "Patran and/or CGNS filename not specified\n");

    if (NULL == (fp = fopen (argv[argind], "r"))) {
        sprintf (buffer, "can't open <%s> for reading", argv[argind]);
        cgnsImportFatal (buffer);
    }
    printf ("reading Patran Neutral file from %s\n", argv[argind]);
    cgnsImportError (print_error);

    getline ();

    /* header - packet 25 */

    if ((packet = atoi (buffer)) == 25) {
        getline ();
        fputs (buffer, stdout);
        getline ();
        packet = atoi (buffer);
    }

    /* summary - packet 26 */

    if (packet == 26) {
        getline ();
        getline ();
        packet = atoi (buffer);
    }

    /* get remaining packets */

    while (packet != 99) {

        /* node */

        if (packet == 1) {
            nodeid = atoi (&buffer[2]);
            getline ();
            p = buffer + 48;
            for (n = 2; n >= 0; n--) {
                *p = 0;
                p -= 16;
                xyz[n] = atof (p);
            }
            getline ();
            cgnsImportNode (nodeid, xyz[0], xyz[1], xyz[2]);
        }

        /* element */

        else if (packet == 2) {
            elemid = atoi (&buffer[2]);
            n = atoi (&buffer[10]);
            nlines = atoi (&buffer[18]);
            if (n == 5 || n == 7 || n == 8) {
                getline ();
                nnodes = n == 8 ? n : n-1;
                lineno++;
                for (n = 0; n < nnodes; n++) {
                    if (1 != fscanf (fp, "%d", &nodeid) || nodeid < 1)
                        cgnsImportFatal ("missing or invalid node ID");
                    nodes[n] = nodeid;
                }
                while (getc (fp) != '\n')
                    ;
                nlines -= 2;
                cgnsImportElement (elemid, nnodes, nodes);
            }
            while (nlines-- > 0)
                getline ();
        }

        /* distributed loads */

        else if (packet == 6 && do_loads) {
            elemid = atoi (&buffer[2]);
            loadid = atoi (&buffer[10]);
            nlines = atoi (&buffer[18]);

            if (loadid != lastid) {
                sprintf (buffer, "PatranLoad%d", loadid);
                cgnsImportBegReg (buffer, cgnsREG_NODES);
                lastid = loadid;
            }
            getline ();

            /* add if element load flag is set */

            if ('1' == buffer[0])
                add_face (elemid, &buffer[9]);
            while (--nlines > 0)
                getline ();
        }

        /* named component */

        else if (packet == 21) {
            int cnt, type, id;
            elemid = atoi (&buffer[2]);
            nnodes = atoi (&buffer[10]) / 2;
            getline ();

            /* strip leading and trailing spaces */

            buffer[sizeof(buffer)-1] = 0;
            p = buffer + strlen (buffer);
            while (--p >= buffer && isspace(*p))
                ;
            *++p = 0;
            for (p = buffer; *p && isspace(*p); p++)
                ;
            cgnsImportBegReg (p, cgnsREG_NODES);

            /* currently only handle type 5 (nodes) in groups */

            for (n = 0, cnt = 0; n < nnodes; n++) {
                if (0 == (n % 5))
                    lineno++;
                fscanf (fp, "%d%d", &type, &id);
                if (5 == type) {
                    nodes[cnt++] = id;
                    if (8 == cnt) {
                        cgnsImportAddReg (8, nodes);
                        cnt = 0;
                    }
                }
            }
            while (getc (fp) != '\n')
                ;
            if (cnt)
                cgnsImportAddReg (cnt, nodes);
            cgnsImportEndReg ();
        }

        /* all others */

        else {
            nlines = atoi (&buffer[18]);
            while (nlines--)
                getline ();
        }

        if (NULL == fgets (buffer, sizeof(buffer), fp))
            break;
        lineno++;
        packet = atoi (buffer);
    }

    fclose (fp);
    cgnsImportError (NULL);

    /* check for duplicate nodes */

    if (do_chk) {
        printf ("checking for duplicate nodes...\n");
        cgnsImportCheck (do_chk == 'd');
    }

    /* output to CGNS file */

    printf ("writing CGNS data to %s\n", argv[++argind]);
    cgnsImportOpen (argv[argind]);
    if (basename != NULL)
        cgnsImportBase (basename);
    cgnsImportWrite ();
    cgnsImportClose ();
    return 0;
}
