#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "getargs.h"
#include "cgnslib.h"
#include "cgnsutil.h"

static char options[] = "ab:B:S:w";

static char *usgmsg[] = {
    "usage: cgns_to_tecplot [options] CGNSfile Tecplotfile",
    "  options are:",
    "   -a       = write ascii instead of binary format",
    "   -b<base> = use CGNS base number <base> (default 1)",
    "   -B<name> = use CGNS base named <name>",
    "   -S<sol>  = solution to use, 0 for no solution (default 1)",
    "   -w       = use volume weighting",
    NULL
};

static int weighting = 0;
static int usesol = -1;
static int ascii = 0;

/*--------------------------------------------------------------------*/

static int check_solution (void)
{
    int nz, nf, nvar = -1;

    if (!usesol) return 0;

    for (nz = 0; nz < nZones; nz++) {
        read_zone_solution (nz+1);
        if (usesol < 0) {
            if (Zones[nz].nsols == 0) {
                usesol = 0;
                return 0;
            }
            nf = Zones[nz].sols[0].nflds;
        }
        else {
            if (usesol > Zones[nz].nsols)
                FATAL (NULL, "solution index out of range");
            nf = Zones[nz].sols[usesol-1].nflds;
        }
        if (nvar < 0)
            nvar = nf;
        else {
            if (nvar != nf) {
                if (usesol < 0) {
                    usesol = 0;
                    return 0;
                }
                FATAL (NULL, "number solution fields not the same for all zones");
            }
        }
    }
    if (usesol < 0) usesol = 1;
    return nvar;
}

/*--------------------------------------------------------------------*/

static void write_string (FILE *fp, char *data)
{
    int c;

    if (ascii)
        fprintf (fp, "\"%s\"\n", data);
    else {
        do {
            c = (int)*data;
            fwrite (&c, 4, 1, fp);
        }
        while (*data++);
    }
}

/*--------------------------------------------------------------------*/

static void write_ints (FILE *fp, int cnt, int *data)
{
    if (ascii) {
        int nout = 0;
        while (cnt-- > 0) {
            if (nout == 10) {
                putc ('\n', fp);
                nout = 0;
            }
            if (nout++)
                putc (' ', fp);
            fprintf (fp, "%d", *data);
            data++;
        }
        putc ('\n', fp);
    }
    else {
        fwrite (data, sizeof(int), cnt, fp);
    }
}

/*--------------------------------------------------------------------*/

static void write_floats (FILE *fp, int cnt, float *data)
{
    if (ascii) {
        int nout = 0;
        while (cnt-- > 0) {
            if (nout == 5) {
                putc ('\n', fp);
                nout = 0;
            }
            if (nout++)
                putc (' ', fp);
            fprintf (fp, "%#g", *data);
            data++;
        }
        putc ('\n', fp);
    }
    else
        fwrite (data, sizeof(float), cnt, fp);
}

/*--------------------------------------------------------------------*/

static void count_elements (int nz, cgsize_t *nelems, int *elemtype)
{
    int ns, et;
    cgsize_t n, nn, ne, nb = 0, nt = 0;
    ZONE *z = &Zones[nz];

    if (z->esets == NULL)
        read_zone_element (nz+1);
    for (ns = 0; ns < z->nesets; ns++) {
        ne = z->esets[ns].end - z->esets[ns].start + 1;
        et = z->esets[ns].type;
        if (et == CGNS_ENUMV(MIXED)) {
            for (n = 0, nn = 0; nn < ne; nn++) {
                et = (int)z->esets[ns].conn[n++];
                if (et < CGNS_ENUMV(NODE) || et > CGNS_ENUMV(HEXA_27))
                    FATAL (NULL, "unrecognized element type");
                if (et >= CGNS_ENUMV(TETRA_4)) {
                    if (et > CGNS_ENUMV(TETRA_10))
                        nb++;
                    else
                        nt++;
                }
                n += element_node_counts[et];
            }
            continue;
        }
        if (et >= CGNS_ENUMV(TETRA_4) && et <= CGNS_ENUMV(HEXA_27)) {
            if (et > CGNS_ENUMV(TETRA_10))
                nb += ne;
            else
                nt += ne;
        }
    }

    if (nt + nb == 0)
        FATAL (NULL, "no volume elements in the zone");
    *nelems = nt + nb;
    *elemtype = nb ? 3 : 2;
}

/*--------------------------------------------------------------------*/

static int *volume_elements (int nz, cgsize_t *nelems, int *elemtype)
{
    int i, np, ns, nt, et;
    int *elems;
    cgsize_t n, nn, ne;
    ZONE *z = &Zones[nz];

    count_elements (nz, &ne, &nt);
    *nelems = ne;
    *elemtype = nt;

    elems = (int *) malloc ((size_t)ne * (1 << nt) * sizeof(int));
    if (NULL == elems)
        FATAL (NULL, "malloc failed for elements");

    for (np = 0, ns = 0; ns < z->nesets; ns++) {
        ne = z->esets[ns].end - z->esets[ns].start + 1;
        et = z->esets[ns].type;
        if (et < CGNS_ENUMV(TETRA_4) || et > CGNS_ENUMV(MIXED)) continue;
        for (n = 0, nn = 0; nn < ne; nn++) {
            if (z->esets[ns].type == CGNS_ENUMV(MIXED))
                et = (int)z->esets[ns].conn[n++];
            switch (et) {
                case CGNS_ENUMV(TETRA_4):
                case CGNS_ENUMV(TETRA_10):
                    if (nt == 2) {
                        for (i = 0; i < 4; i++)
                            elems[np++] = (int)z->esets[ns].conn[n+i];
                    }
                    else {
                        for (i = 0; i < 3; i++)
                            elems[np++] = (int)z->esets[ns].conn[n+i];
                        elems[np++] = (int)z->esets[ns].conn[n+2];
                        for (i = 0; i < 4; i++)
                            elems[np++] = (int)z->esets[ns].conn[n+3];
                    }
                    break;
                case CGNS_ENUMV(PYRA_5):
                case CGNS_ENUMV(PYRA_14):
                    for (i = 0; i < 4; i++)
                        elems[np++] = (int)z->esets[ns].conn[n+i];
                    for (i = 0; i < 4; i++)
                        elems[np++] = (int)z->esets[ns].conn[n+4];
                    break;
                case CGNS_ENUMV(PENTA_6):
                case CGNS_ENUMV(PENTA_15):
                case CGNS_ENUMV(PENTA_18):
                    for (i = 0; i < 3; i++)
                        elems[np++] = (int)z->esets[ns].conn[n+i];
                    elems[np++] = (int)z->esets[ns].conn[n+2];
                    for (i = 3; i < 6; i++)
                        elems[np++] = (int)z->esets[ns].conn[n+i];
                    elems[np++] = (int)z->esets[ns].conn[n+5];
                    break;
                case CGNS_ENUMV(HEXA_8):
                case CGNS_ENUMV(HEXA_20):
                case CGNS_ENUMV(HEXA_27):
                    for (i = 0; i < 8; i++)
                        elems[np++] = (int)z->esets[ns].conn[n+i];
                    break;
            }
            n += element_node_counts[et];
        }
    }

    return elems;
}

/*--------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
    int i, n, ib, nb, nz, nv, celldim, phydim;
    int nn, type, *elems = 0, idata[5];
    cgsize_t ne;
    char *p, basename[33], title[65];
    float value, *var;
    SOLUTION *sol;
    FILE *fp;

    if (argc < 2)
        print_usage (usgmsg, NULL);

    ib = 0;
    basename[0] = 0;
    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'a':
                ascii = 1;
                break;
            case 'b':
                ib = atoi (argarg);
                break;
            case 'B':
                strncpy (basename, argarg, 32);
                basename[32] = 0;
                break;
            case 'w':
                weighting = 1;
                break;
            case 'S':
                usesol = atoi (argarg);
                break;
        }
    }

    if (argind > argc - 2)
        print_usage (usgmsg, "CGNSfile and/or Tecplotfile not given");
    if (!file_exists (argv[argind]))
        FATAL (NULL, "CGNSfile does not exist or is not a file");

    /* open CGNS file */

    printf ("reading CGNS file from %s\n", argv[argind]);
    nb = open_cgns (argv[argind], 1);
    if (!nb)
        FATAL (NULL, "no bases found in CGNS file");
    if (*basename && 0 == (ib = find_base (basename)))
        FATAL (NULL, "specified base not found");
    if (ib > nb) FATAL (NULL, "base index out of range");
    cgnsbase = ib ? ib : 1;
    if (cg_base_read (cgnsfn, cgnsbase, basename, &celldim, &phydim))
        FATAL (NULL, NULL);
    if (celldim != 3 || phydim != 3)
        FATAL (NULL, "cell and physical dimension must be 3");
    printf ("  using base %d - %s\n", cgnsbase, basename);

    if (NULL == (p = strrchr (argv[argind], '/')) &&
        NULL == (p = strrchr (argv[argind], '\\')))
        strncpy (title, argv[argind], sizeof(title));
    else
        strncpy (title, ++p, sizeof(title));
    title[sizeof(title)-1] = 0;
    if ((p = strrchr (title, '.')) != NULL)
        *p = 0;

    read_zones ();
    if (!nZones)
        FATAL (NULL, "no zones in the CGNS file");
    
    /* verify dimensions fit in an integer */

    for (nz = 0; nz < nZones; nz++) {
        if (Zones[nz].nverts > CG_MAX_INT32)
	    FATAL(NULL, "zone size too large to write with integers");
	if (Zones[nz].type == CGNS_ENUMV(Unstructured)) {
            count_elements (nz, &ne, &type);
            if (ne > CG_MAX_INT32)
	        FATAL(NULL, "too many elements to write with integers");
        }
     }

    nv = 3 + check_solution ();

    /* open Tecplot file */

    printf ("writing %s Tecplot data to <%s>\n",
        ascii ? "ASCII" : "binary", argv[++argind]);
    if (NULL == (fp = fopen (argv[argind], ascii ? "w+" : "w+b")))
        FATAL (NULL, "couldn't open Tecplot output file");

    /* write file header */

    if (ascii)
        fprintf (fp, "TITLE = \"%s\"\n", title);
    else {
        fwrite ("#!TDV75 ", 1, 8, fp);
        i = 1;
        write_ints (fp, 1, &i);
        write_string (fp, title);
    }

    /* write variables */

    if (ascii) {
        fprintf (fp, "VARIABLES = \"X\", \"Y\", \"Z\"");
        if (usesol) {
            sol = Zones->sols;
            for (n = 0; n < sol->nflds; n++)
                fprintf (fp, ",\n\"%s\"", sol->flds[n].name);
        }
    }
    else {
        write_ints (fp, 1, &nv);
        write_string (fp, "X");
        write_string (fp, "Y");
        write_string (fp, "Z");
        if (usesol) {
            sol = Zones->sols;
            for (n = 0; n < sol->nflds; n++)
                write_string (fp, sol->flds[n].name);
        }
    }

    /* write zones */

    if (!ascii) {
        for (nz = 0; nz < nZones; nz++) {
            if (Zones[nz].type == CGNS_ENUMV(Structured)) {
                idata[0] = 0;          /* BLOCK */
                idata[1] = -1;         /* color not specified */
                idata[2] = (int)Zones[nz].dim[0];
                idata[3] = (int)Zones[nz].dim[1];
                idata[4] = (int)Zones[nz].dim[2];
            }
            else {
                count_elements (nz, &ne, &type);
                idata[0] = 2;          /* FEBLOCK */
                idata[1] = -1;         /* color not specified */
                idata[2] = (int)Zones[nz].dim[0];
                idata[3] = (int)ne;
                idata[4] = type;
            }
            value = 299.0;
            write_floats (fp, 1, &value);
            write_string (fp, Zones[nz].name);
            write_ints (fp, 5, idata);
        }
        value = 357.0;
        write_floats (fp, 1, &value);
    }

    for (nz = 0; nz < nZones; nz++) {
        printf ("  zone %d...", nz+1);
        fflush (stdout);
        read_zone_grid (nz+1);
        nn = (int)Zones[nz].nverts;
        var = (float *) malloc (nn * sizeof(float));
        if (NULL == var)
            FATAL (NULL, "malloc failed for temp float array");
        if (Zones[nz].type == CGNS_ENUMV(Unstructured))
            elems = volume_elements (nz, &ne, &type);

        if (ascii) {
            if (Zones[nz].type == CGNS_ENUMV(Structured))
                fprintf (fp, "\nZONE T=\"%s\", I=%d, J=%d, K=%d, F=BLOCK\n",
                    Zones[nz].name, (int)Zones[nz].dim[0],
                    (int)Zones[nz].dim[1], (int)Zones[nz].dim[2]);
            else
                fprintf (fp, "\nZONE T=\"%s\", N=%d, E=%d, F=FEBLOCK, ET=%s\n",
                    Zones[nz].name, nn, (int)ne, type == 2 ? "TETRAHEDRON" : "BRICK");
        }
        else {
            value = 299.0;
            write_floats (fp, 1, &value);
            i = 0;
            write_ints (fp, 1, &i);
            i = 1;
            for (n = 0; n < nv; n++)
                write_ints (fp, 1, &i);
        }

        for (n = 0; n < nn; n++)
            var[n] = (float)Zones[nz].verts[n].x;
        write_floats (fp, nn, var);
        for (n = 0; n < nn; n++)
            var[n] = (float)Zones[nz].verts[n].y;
        write_floats (fp, nn, var);
        for (n = 0; n < nn; n++)
            var[n] = (float)Zones[nz].verts[n].z;
        write_floats (fp, nn, var);

        if (usesol) {
            read_solution_field (nz+1, usesol, 0);
            sol = &Zones[nz].sols[usesol-1];
            if (sol->location != CGNS_ENUMV(Vertex))
                cell_vertex_solution (nz+1, usesol, weighting);
            for (nv = 0; nv < sol->nflds; nv++) {
                for (n = 0; n < nn; n++)
                    var[n] = (float)sol->flds[nv].data[n];
                write_floats (fp, nn, var);
            }
        }

        free (var);

        if (Zones[nz].type == CGNS_ENUMV(Unstructured)) {
            if (!ascii) {
                i = 0;
                write_ints (fp, 1, &i);
            }
            nn = 1 << type;
            for (i = 0, n = 0; n < ne; n++, i += nn)
                write_ints (fp, nn, &elems[i]);
            free (elems);
        }
        puts ("done");
    }

    fclose (fp);
    cg_close (cgnsfn);
    return 0;
}
