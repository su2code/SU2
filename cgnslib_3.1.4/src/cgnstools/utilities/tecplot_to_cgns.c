/*
 * tecplot_to_cgns - convert Tecplot to CGNS
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#ifdef _WIN32
# define strcasecmp(A,B)    _stricmp(A,B)
# define strncasecmp(A,B,C) _strnicmp(A,B,C)
#endif

#include "cgnsImport.h"
#include "getargs.h"

static char options[] = "fdDt:B:";

static char *usgmsg[] = {
    "usage  : tecplot_to_cgns [options] Tecplotfile CGNSfile",
    "options:",
    "   -f       = fix degenerate brick elements",
    "   -d       = duplicate node checking with rel tol",
    "   -D       = duplicate node checking with abs tol",
    "   -t<tol>  = duplicate node checking tolerance",
    "   -B<name> = set CGNS base name to <name>",
    NULL
};

static int nvar = 3;
static int xloc = 0;
static int yloc = 1;
static int zloc = 2;

static char buffer[257], name[33];

/*---------- check_ascii --------------------------------------------
 * check first 256 bytes for non-ascii characters
 *-------------------------------------------------------------------*/

static void check_ascii (char *fname)
{
    int n, np;
    FILE *fp = fopen (fname, "rb");

    if (fp == NULL) {
        fprintf (stderr, "can't open <%s> for reading\n", fname);
        exit (1);
    }
    np = (int)fread (buffer, 1, sizeof(buffer), fp);
    fclose (fp);
    for (n = 0; n < np; n++) {
        if (!buffer[n] || !isascii (buffer[n]) ||
           (buffer[n] < ' ' && !isspace (buffer[n]))) {
            fprintf (stderr, "can only read ASCII Tecplot files\n");
            exit (1);
        }
    }
}

/*---------- get_line -----------------------------------------------
 * get next non-blank line
 *-------------------------------------------------------------------*/

static char *get_line (FILE *fp)
{
    char *p;

    while (fgets (buffer, sizeof(buffer), fp) != NULL) {
        buffer[sizeof(buffer)-1] = 0;
        for (p = buffer+strlen(buffer)-1; p >= buffer && isspace(*p); p--)
            ;
        *++p = 0;
        for (p = buffer; *p && isspace (*p); p++)
            ;
        if (*p) return p;
    }
    return NULL;
}

/*---------- getvar -------------------------------------------------
 * get the next variable name
 *-------------------------------------------------------------------*/

static char *getvar (char **pp)
{
    int n = 0;
    char *p;

    for (p = *pp; *p && (isspace (*p) || *p == ','); p++)
        ;
    if (*p != '"') return NULL;
    for (++p; *p && *p != '"'; p++) {
        if (n < sizeof(name))
            name[n++] = *p;
    }
    name[n] = 0;
    if (*p == '"') p++;
    *pp = p;
    return name;
}

/*---------- getzone ------------------------------------------------
 * get zone description data
 *-------------------------------------------------------------------*/

static char *getzone (char **pp, int *tag)
{
    int n, match;
    char *p, what[5];

    for (p = *pp; *p && (isspace (*p) || *p == ','); p++)
        ;
    if (!*p) return NULL;
    for (n = 0; *p && !isspace(*p) && *p != '='; p++) {
        if (n < sizeof(what))
            what[n++] = *p;
    }
    what[n] = 0;
    while (*p && isspace (*p))
        p++;
    if (*p != '=')
        cgnsImportFatal ("zone specification error");
    while (*++p && isspace (*p))
        ;
    if (*p == '"' || *p == '(') {
        match = *p == '"' ? '"' : ')';
        for (n = 0, ++p; *p && *p != match; p++) {
            if (n < sizeof(name))
                name[n++] = *p;
        }
        if (*p == match) p++;
    }
    else if (*p) {
        for (n = 0; *p && !isspace (*p) && *p != ','; p++) {
            if (n < sizeof(name))
                name[n++] = *p;
        }
    }
    else
        cgnsImportFatal ("zone specification error");
    name[n] = 0;

    switch (tolower (what[0])) {
        case 't': *tag = 0; break;
        case 'f': *tag = 1; break;
        case 'i': *tag = 2; break;
        case 'j': *tag = 3; break;
        case 'k': *tag = 4; break;
        case 'n': *tag = 5; break;
        case 'e': *tag = tolower(what[1]) == 't' ? 7 : 6; break;
        default: *tag = -1; break;
    }
    *pp = p;
    return name;
}

/*---------- block_nodes --------------------------------------------
 * read nodes in block format
 *-------------------------------------------------------------------*/

static void block_nodes (FILE *fp, int nnodes)
{
    int nv, nn, vnum;
    double *x, *y, *z, v;

    for (nn = 0; nn < nnodes; nn++)
        cgnsImportNode (nn+1, 0.0, 0.0, 0.0);

    x = (double *) calloc (3 * nnodes, sizeof(double));
    if (x == NULL)
        cgnsImportFatal ("malloc failed for variables");
    y = x + nnodes;
    z = y + nnodes;

    for (vnum = 0, nv = 0; nv < nvar; nv++) {
        if (nv == xloc) {
            for (nn = 0; nn < nnodes; nn++) {
                if (1 != fscanf (fp, "%lf", &x[nn]))
                    cgnsImportFatal ("error reading variables");
            }
        }
        else if (nv == yloc) {
            for (nn = 0; nn < nnodes; nn++) {
                if (1 != fscanf (fp, "%lf", &y[nn]))
                    cgnsImportFatal ("error reading variables");
            }
        }
        else if (nv == zloc) {
            for (nn = 0; nn < nnodes; nn++) {
                if (1 != fscanf (fp, "%lf", &z[nn]))
                    cgnsImportFatal ("error reading variables");
            }
        }
        else {
            for (nn = 0; nn < nnodes; nn++) {
                if (1 != fscanf (fp, "%lf", &v))
                    cgnsImportFatal ("error reading variables");
                cgnsImportVariable(nn+1, vnum++, v);
            }
        }
    }
    for (nn = 0; nn < nnodes; nn++)
        cgnsImportSetNode (nn+1, x[nn], y[nn], z[nn]);
    free (x);
}

/*---------- point_nodes --------------------------------------------
 * read nodes in point format
 *-------------------------------------------------------------------*/

static void point_nodes (FILE *fp, int nnodes)
{
    int nv, nn, vnum;
    double x = 0.0, y = 0.0, z = 0.0;
    double *v;

    v = (double *) malloc (nvar * sizeof(double));
    if (v == NULL)
        cgnsImportFatal ("malloc failed for variables");

    for (nn = 0; nn < nnodes; nn++) {
        for (nv = 0; nv < nvar; nv++) {
            if (1 != fscanf (fp, "%lf", &v[nv]))
                cgnsImportFatal ("error reading variables");
        }
        if (xloc >= 0 && xloc < nvar) x = v[xloc];
        if (yloc >= 0 && yloc < nvar) y = v[yloc];
        if (zloc >= 0 && zloc < nvar) z = v[zloc];
        cgnsImportNode (nn+1, x, y, z);
        for (vnum = 0, nv = 0; nv < nvar; nv++) {
            if (nv != xloc && nv != yloc && nv != zloc)
                cgnsImportVariable(nn+1, vnum++, v[nv]);
        }
    }
    free (v);
}

/*========== main ===================================================*/

int main (int argc, char *argv[])
{
    int n, i, j, k, ni, nj, nk;
    int nn, ne, et, nz, block;
    cgsize_t nodes[8];
    int do_chk = 0, fix_bricks = 0;
    char elemname[33];
    char zonename[33], *p, *s, *basename = NULL;
    FILE *fp;

    if (argc < 2)
        print_usage (usgmsg, NULL);

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'f':
                fix_bricks = 1;
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
        print_usage (usgmsg, "Tecplot and/or CGNS filename not specified\n");

    check_ascii (argv[argind]);
    if (NULL == (fp = fopen (argv[argind], "r"))) {
        fprintf (stderr, "can't open <%s> for reading\n", argv[argind]);
        exit (1);
    }
    printf ("reading Tecplot file from %s\n", argv[argind++]);
    printf ("writing CGNS data to %s\n", argv[argind]);
    cgnsImportOpen (argv[argind]);
    if (basename != NULL)
        cgnsImportBase (basename);

    nz = 0;
    p = get_line (fp);

    while (p != NULL) {

        /* VARIABLES */

        if (0 == strncasecmp (p, "variables", 9)) {
            for (p += 9; *p && isspace(*p); p++)
                ;
            if (*p++ != '=') {
                p = get_line (fp);
                continue;
            }
            nvar = 0;
            xloc = yloc = zloc = -1;
            while (1) {
                if ((s = getvar (&p)) == NULL) {
                    p = get_line (fp);
                    if (p == NULL || (s = getvar (&p)) == NULL) break;
                }
                if (0 == strcasecmp ("x", s))
                    xloc = nvar++;
                else if (0 == strcasecmp ("y", s))
                    yloc = nvar++;
                else if (0 == strcasecmp ("z", s))
                    zloc = nvar++;
                else {
                    cgnsImportAddVariable(s);
                    nvar++;
                }
            }
            if (xloc == -1 && yloc == -1 && zloc == -1)
                cgnsImportFatal("X, Y and Z variables missing");
            continue;
        }

        /* ZONE */

        if (0 == strncasecmp (p, "zone", 4)) {
            p += 4;
            ni = nj = nk = nn = ne = et = -1;
            block = 0;
            sprintf (zonename, "Zone%d", ++nz);
            while (1) {
                if ((s = getzone (&p, &n)) == NULL) {
                    while ((n = fgetc (fp)) != EOF &&
                        (isspace(n) || n == ','))
                        ;
                    if (n == EOF)
                        cgnsImportFatal ("end of file while reading zone");
                    ungetc (n, fp);
                    if (!isalpha (n)) break;
                    p = get_line (fp);
                    continue;
                }
                switch (n) {
                    case 0:     /* T */
                        strncpy (zonename, s, 32);
                        zonename[32] = 0;
                        break;
                    case 1:     /* F */
                        if (0 == strncasecmp (s, "block", 5) ||
                            0 == strncasecmp (s, "feblock", 7))
                            block = 1;
                        break;
                    case 2:     /* I */
                        ni = atoi (s);
                        break;
                    case 3:     /* J */
                        nj = atoi (s);
                        break;
                    case 4:     /* K */
                        nk = atoi (s);
                        break;
                    case 5:     /* N */
                        nn = atoi (s);
                        break;
                    case 6:     /* E */
                        ne = atoi (s);
                        break;
                    case 7:     /* ET */
                        strncpy (elemname, s, 32);
                        elemname[32] = 0;
                        if (0 == strncasecmp (elemname, "tri", 3))
                            et = 0;
                        else if (0 == strncasecmp (elemname, "qua", 3))
                            et = 1;
                        else if (0 == strncasecmp (elemname, "tet", 3))
                            et = 2;
                        else if (0 == strncasecmp (elemname, "bri", 3))
                            et = 3;
                        else
                            printf("unhandled element type %s", elemname);
                        break;
                    default:
                        break;
                }
            }

            printf("zone %s:", zonename);

            if (nn == -1) {
                if (ni < 2 || nj < 2 || nk < 2) {
                    printf("missing I, J and K - skipping zone\n");
                    p = get_line(fp);
                    continue;
                }
                nn = ni * nj * nk;
                ne = (ni - 1) * (nj - 1) * (nk - 1);
                et = 5;
            }
            else {
                if (nn < 3 || ne < 1 || et < 0) {
                    printf("%d nodes, %d %s elements - skipping zone\n",
                        nn, ne, elemname);
                    p = get_line(fp);
                    continue;
                }
                if (et < 2) {
                    printf ("%s elements - skipping zone\n", elemname);
                    p = get_line (fp);
                    continue;
                }
            }

            printf ("%d nodes...", nn);
            fflush (stdout);

            cgnsImportZone (zonename);
            if (block)
                block_nodes (fp, nn);
            else
                point_nodes (fp, nn);

            printf (" %d %s elements...", ne, elemname);
            fflush (stdout);

            if (et == 5) {
                for (n = 1, k = 1; k < nk; k++) {
                    for (j = 1; j < nj; j++) {
                        for (i = 1; i < ni; i++) {
                            n = i + ni * ((j - 1) + nj * (k - 1));
                            nodes[0] = n;
                            nodes[1] = n + 1;
                            nodes[2] = n + 1 + ni;
                            nodes[3] = n + ni;
                            n = i + ni * ((j - 1) + nj * k);
                            nodes[4] = n;
                            nodes[5] = n + 1;
                            nodes[6] = n + 1 + ni;
                            nodes[7] = n + ni;
                            cgnsImportElement (n++, 8, nodes);
                        }
                    }
                }
            }
            else {
                j = 1 << et;
                for (n = 1; n <= ne; n++) {
                    for (i = 0; i < j; i++) {
                        if (1 != fscanf (fp, "%d", &k))
                            cgnsImportFatal ("error reading elements");
                        nodes[i] = k;
                    }
                    k = j;
                    if (fix_bricks && k == 8 && nodes[6] == nodes[7]) {
                        if (nodes[4] == nodes[6] && nodes[5] == nodes[6]) {
                            if (nodes[2] == nodes[3]) {
                                k = 4;
                                nodes[3] = nodes[4];
                            }
                            else
                                k = 5;
                        }
                        else if (nodes[2] == nodes[3]) {
                            k = 6;
                            nodes[3] = nodes[4];
                            nodes[4] = nodes[5];
                            nodes[5] = nodes[6];
                        }
                    }
                    cgnsImportElement (n, k, nodes);
                }
            }
            puts (" done");

            if (do_chk) {
                printf ("checking for duplicate nodes...\n");
                cgnsImportCheck (do_chk == 'd');
            }

            cgnsImportWrite ();
        }

        p = get_line (fp);
    }

    fclose (fp);
    cgnsImportClose ();
    return 0;
}
