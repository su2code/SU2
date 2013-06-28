#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
# include <io.h>
# define unlink _unlink
#else
# include <unistd.h>
#endif
#include "cgnslib.h"

#ifndef CGNSTYPES_H
# define cgsize_t int
#endif
#ifndef CGNS_ENUMT
# define CGNS_ENUMT(e) e
# define CGNS_ENUMV(e) e
#endif

#define NUM_SIDE 5
#define NUM_RANDOM 10

float *xcoord;
float *ycoord;
float *zcoord;
float *fbuf;

cgsize_t *elements;
cgsize_t *faces;
cgsize_t *parent;
cgsize_t *ptmp;
cgsize_t *ibuf, *pbuf;

#define NODE_INDEX(I,J,K) ((I)+NUM_SIDE*(((J)-1)+NUM_SIDE*((K)-1)))
#define CELL_INDEX(I,J,K) ((I)+(NUM_SIDE-1)*(((J)-1)+(NUM_SIDE-1)*((K)-1)))

int irandom(int imin, int imax) {
    float r = (float)rand() / (float)RAND_MAX;
    int i = imin + (int)(r * (float)(imax - imin));
    if (i < imin) return imin;
    if (i > imax) return imax;
    return i;
}

void get_parent(int rmin, int rmax, int nelems, int nfaces)
{
    int i, j, k, n, nn, np;

    np = rmax - rmin + 1;
    nn = rmin - nelems - 1;
    for (n = 0, j = 0; j < 4; j++) {
        k = j * nfaces + nn;
        for (i = 0; i < np; i++)
            ptmp[n++] = parent[k++];
    }
}

int mixed_offset(int num, int nelems) {
    int offset;
    int nmixed = nelems << 1;

    if (--num < nmixed) {
        int i = num >> 1;
        offset = 14 * i;
        if (num != (i << 1)) offset += 5;
    }
    else
        offset = 14 * nelems + 5 * (num - nmixed);
    return offset;
}

int main (int argc, char **argv)
{
    int n, i, j, k, nn, nf, np;
    int nnodes, nelems, nfaces;
    int cgfile, cgbase, cgzone, cgcoord, cgsol, cgfld;
    int cgsect, cgelems, cgfaces;
    cgsize_t size[3], rmin, rmax, nr;
    cgsize_t is, ie;
    char name[33];
    CGNS_ENUMT( ElementType_t )  type;
    static char *fname = "partial.cgns";

    if (argc > 1) {
        n = 0;
        if (argv[1][n] == '-') n++;
        if (argv[1][n] == 'a' || argv[1][n] == 'A')
            nn = CG_FILE_ADF;
        else if (argv[1][n] == 'h' || argv[1][n] == 'H')
            nn = CG_FILE_HDF5;
        else {
            fprintf(stderr, "unknown option\n");
            exit (1);
        }
        if (cg_set_file_type(nn))
            cg_error_exit();
    }

    nnodes = NUM_SIDE * NUM_SIDE * NUM_SIDE;
    xcoord = (float *) malloc (4 * nnodes * sizeof(float));
    ycoord = xcoord + nnodes;
    zcoord = ycoord + nnodes;
    fbuf   = zcoord + nnodes;

    for (n = 0, k = 1; k <= NUM_SIDE; k++) {
        for (j = 1; j <= NUM_SIDE; j++) {
            for (i = 1; i <= NUM_SIDE; i++) {
                xcoord[n] = (float)i;
                ycoord[n] = (float)j;
                zcoord[n] = (float)k;
                n++;
            }
        }
    }

    nelems = (NUM_SIDE - 1) * (NUM_SIDE - 1) * (NUM_SIDE - 1);
    elements = (cgsize_t *) malloc (16 * nelems * sizeof(cgsize_t));
    ibuf = elements + 8 * nelems;

    for (n = 0, k = 1; k < NUM_SIDE; k++) {
        for (j = 1; j < NUM_SIDE; j++) {
            for (i = 1; i < NUM_SIDE; i++) {
                nn = NODE_INDEX(i, j, k);
                elements[n++] = nn;
                elements[n++] = nn + 1;
                elements[n++] = nn + 1 + NUM_SIDE;
                elements[n++] = nn + NUM_SIDE;
                nn += NUM_SIDE * NUM_SIDE;
                elements[n++] = nn;
                elements[n++] = nn + 1;
                elements[n++] = nn + 1 + NUM_SIDE;
                elements[n++] = nn + NUM_SIDE;
            }
        }
    }

    nfaces = 6 * (NUM_SIDE - 1) * (NUM_SIDE - 1);
    faces = (cgsize_t *) malloc (18 * nfaces * sizeof(cgsize_t));
    parent = faces + 4 * nfaces;
    ptmp = parent + 4 * nfaces;
    pbuf = ptmp + 4 * nfaces;

    for (n = 0; n < 4 * nfaces; n++)
        parent[n] = 0;
    np = nf = 0;
    n = 2 * nfaces;

    i = 1;
    for (k = 1; k < NUM_SIDE; k++) {
        for (j = 1; j < NUM_SIDE; j++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + NUM_SIDE * NUM_SIDE;
            faces[nf++]  = nn + NUM_SIDE * (NUM_SIDE + 1);
            faces[nf++]  = nn + NUM_SIDE;
            parent[np]   = CELL_INDEX(i, j, k);
            parent[np+n] = 5;
            np++;
        }
    }
    i = NUM_SIDE;
    for (k = 1; k < NUM_SIDE; k++) {
        for (j = 1; j < NUM_SIDE; j++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + NUM_SIDE;
            faces[nf++]  = nn + NUM_SIDE * (NUM_SIDE + 1);
            faces[nf++]  = nn + NUM_SIDE * NUM_SIDE;
            parent[np]   = CELL_INDEX(i-1, j, k);
            parent[np+n] = 3;
            np++;
        }
    }
    j = 1;
    for (k = 1; k < NUM_SIDE; k++) {
        for (i = 1; i < NUM_SIDE; i++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + 1;
            faces[nf++]  = nn + 1 + NUM_SIDE * NUM_SIDE;
            faces[nf++]  = nn + NUM_SIDE * NUM_SIDE;
            parent[np]   = CELL_INDEX(i, j, k);
            parent[np+n] = 2;
            np++;
        }
    }
    j = NUM_SIDE;
    for (k = 1; k < NUM_SIDE; k++) {
        for (i = 1; i < NUM_SIDE; i++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + NUM_SIDE * NUM_SIDE;
            faces[nf++]  = nn + 1 + NUM_SIDE * NUM_SIDE;
            faces[nf++]  = nn + 1;
            parent[np]   = CELL_INDEX(i, j-1, k);
            parent[np+n] = 4;
            np++;
        }
    }
    k = 1;
    for (j = 1; j < NUM_SIDE; j++) {
        for (i = 1; i < NUM_SIDE; i++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + NUM_SIDE;
            faces[nf++]  = nn + NUM_SIDE + 1;
            faces[nf++]  = nn + 1;
            parent[np]   = CELL_INDEX(i, j, k);
            parent[np+n] = 1;
            np++;
        }
    }
    k = NUM_SIDE;
    for (j = 1; j < NUM_SIDE; j++) {
        for (i = 1; i < NUM_SIDE; i++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + 1;
            faces[nf++]  = nn + NUM_SIDE + 1;
            faces[nf++]  = nn + NUM_SIDE;
            parent[np]   = CELL_INDEX(i, j, k-1);
            parent[np+n] = 6;
            np++;
        }
    }

    unlink (fname);
    printf ("creating CGNS file %s\n", fname);
    if (cg_open (fname, CG_MODE_WRITE, &cgfile) ||
        cg_base_write (cgfile, "Base", 3, 3, &cgbase))
        cg_error_exit ();

    /* create zone and sections */

    size[0] = nnodes;
    size[1] = nelems;
    size[2] = 0;

    if (cg_zone_write (cgfile, cgbase, "Zone", size,
            CGNS_ENUMV(Unstructured), &cgzone) ||
        cg_section_partial_write(cgfile, cgbase, cgzone, "Elements",
	    CGNS_ENUMV(HEXA_8), (cgsize_t)1, (cgsize_t)nelems, 0, &cgelems) ||
        cg_section_partial_write(cgfile, cgbase, cgzone, "Faces",
	    CGNS_ENUMV(QUAD_4), (cgsize_t)(nelems + 1),
	    (cgsize_t)(nelems + nfaces), 0, &cgfaces))
        cg_error_exit();

    /* write zone with partial write */

    puts("\nwriting zone with partial write");

    /* write every other coordinate plane */

    np = NUM_SIDE * NUM_SIDE;
    for (k = 0; k < NUM_SIDE; k += 2) {
        n    = k * np;
        rmin = n + 1;
        rmax = n + np;
        printf("coordinates %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_coord_partial_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
                "CoordinateX", &rmin, &rmax, &xcoord[n], &cgcoord) ||
            cg_coord_partial_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
                "CoordinateY", &rmin, &rmax, &ycoord[n], &cgcoord) ||
            cg_coord_partial_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
                "CoordinateZ", &rmin, &rmax, &zcoord[n], &cgcoord))
            cg_error_exit();
    }

    /* write every other cell plane */

    np = (NUM_SIDE - 1) * (NUM_SIDE - 1);
    for (k = 1; k < NUM_SIDE; k += 2) {
        nn   = (k - 1) * np;
        n    = nn << 3;
        rmin = nn + 1;
        rmax = nn + np;
        printf("elements %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_elements_partial_write(cgfile, cgbase, cgzone,
                cgelems, rmin, rmax, &elements[n]))
            cg_error_exit();
    }

    /* write every other face */

    for (k = 0; k < 6; k += 2) {
        nn   = k * np;
        n    = nn << 2;
        rmin = nn + 1 + nelems;
        rmax = nn + np + nelems;
        get_parent((int)rmin, (int)rmax, nelems, nfaces);
        printf("faces %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_elements_partial_write(cgfile, cgbase, cgzone,
                cgfaces, rmin, rmax, &faces[n]) ||
            cg_parent_data_partial_write(cgfile, cgbase, cgzone,
                cgfaces, rmin, rmax, ptmp))
            cg_error_exit();
    }

    /* write every other solution value */

    if (cg_sol_write(cgfile, cgbase, cgzone, "Solution", CGNS_ENUMV(Vertex), &cgsol))
        cg_error_exit();

    puts("field -> 1,3,5,7 ...");
    for (n = 0; n < nnodes; n += 2) {
        rmin = n + 1;
        rmax = n + 1;
        if (cg_field_partial_write(cgfile, cgbase, cgzone, cgsol,
	        CGNS_ENUMV(RealSingle), "Field", &rmin, &rmax, &xcoord[n], &cgfld))
            cg_error_exit();
    }

    puts ("closing and reopening in modify mode");
    cg_close (cgfile);

    if (cg_open (fname, CG_MODE_MODIFY, &cgfile))
        cg_error_exit ();

    /* fill in missing coordinate planes */

    np = NUM_SIDE * NUM_SIDE;
    for (k = 1; k < NUM_SIDE; k += 2) {
        n    = k * np;
        rmin = n + 1;
        rmax = n + np;
        printf("coordinates %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_coord_partial_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
                "CoordinateX", &rmin, &rmax, &xcoord[n], &cgcoord) ||
            cg_coord_partial_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
                "CoordinateY", &rmin, &rmax, &ycoord[n], &cgcoord) ||
            cg_coord_partial_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
                "CoordinateZ", &rmin, &rmax, &zcoord[n], &cgcoord))
            cg_error_exit();
    }

    /* fill in missing cell planes */

    np = (NUM_SIDE - 1) * (NUM_SIDE - 1);
    for (k = 2; k < NUM_SIDE; k += 2) {
        nn   = (k - 1) * np;
        n    = nn << 3;
        rmin = nn + 1;
        rmax = nn + np;
        printf("elements %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_elements_partial_write(cgfile, cgbase, cgzone,
                cgelems, rmin, rmax, &elements[n]))
            cg_error_exit();
    }

    /* fill in missing faces */

    for (k = 1; k < 6; k += 2) {
        nn   = k * np;
        n    = nn << 2;
        rmin = nn + 1 + nelems;
        rmax = nn + np + nelems;
        get_parent((int)rmin, (int)rmax, nelems, nfaces);
        printf("faces %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_elements_partial_write(cgfile, cgbase, cgzone,
                cgfaces, rmin, rmax, &faces[n]) ||
            cg_parent_data_partial_write(cgfile, cgbase, cgzone,
                cgfaces, rmin, rmax, ptmp))
            cg_error_exit();
    }

    /* fill in missing solution values */

    puts("field -> 2,4,6,8 ...");
    for (n = 1; n < nnodes; n += 2) {
        rmin = n + 1;
        rmax = n + 1;
        if (cg_field_partial_write(cgfile, cgbase, cgzone, cgsol,
		CGNS_ENUMV(RealSingle), "Field", &rmin, &rmax, &xcoord[n], &cgfld))
            cg_error_exit();
    }

#if NUM_RANDOM > 0

    /* do random insertion */

    printf("doing random modification\n");

    for (k = 0; k < NUM_RANDOM; k++) {
        rmin = irandom(1, nnodes);
        rmax = irandom(1, nnodes);
        if (rmin > rmax) {
            nr = rmin;
            rmin = rmax;
            rmax = nr;
        }
        n = (int)(rmin - 1);
        printf("coordinates %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_coord_partial_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
                "CoordinateX", &rmin, &rmax, &xcoord[n], &cgcoord) ||
            cg_coord_partial_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
                "CoordinateY", &rmin, &rmax, &ycoord[n], &cgcoord) ||
            cg_coord_partial_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
                "CoordinateZ", &rmin, &rmax, &zcoord[n], &cgcoord))
            cg_error_exit();
    }

    for (k = 0; k < NUM_RANDOM; k++) {
        rmin = irandom(1, nelems);
        rmax = irandom(1, nelems);
        if (rmin > rmax) {
            nr = rmin;
            rmin = rmax;
            rmax = nr;
        }
        n = (int)(rmin - 1) << 3;
        printf("elements %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_elements_partial_write(cgfile, cgbase, cgzone,
                cgelems, rmin, rmax, &elements[n]))
            cg_error_exit();
    }

    for (k = 0; k < NUM_RANDOM; k++) {
        rmin = irandom(1, nfaces);
        rmax = irandom(1, nfaces);
        if (rmin > rmax) {
            nr = rmin;
            rmin = rmax;
            rmax = nr;
        }
        n = (int)(rmin - 1) << 2;
        rmin += nelems;
        rmax += nelems;
        get_parent((int)rmin, (int)rmax, nelems, nfaces);
        printf("faces %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_elements_partial_write(cgfile, cgbase, cgzone,
                cgfaces, rmin, rmax, &faces[n]) ||
            cg_parent_data_partial_write(cgfile, cgbase, cgzone,
                cgfaces, rmin, rmax, ptmp))
            cg_error_exit();
    }

    for (k = 0; k < NUM_RANDOM; k++) {
        rmin = irandom(1, nnodes);
        rmax = irandom(1, nnodes);
        if (rmin > rmax) {
            nr = rmin;
            rmin = rmax;
            rmax = nr;
        }
        n = (int)rmin - 1;
        printf("field %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_field_partial_write(cgfile, cgbase, cgzone, cgsol,
	        CGNS_ENUMV(RealSingle), "Field", &rmin, &rmax, &xcoord[n], &cgfld))
            cg_error_exit();
    }

#endif

    /* check the data */

    puts("checking the data");

    nn = 0;

    /* check coordinates */

    rmin = 1;
    rmax = nnodes;
    if (cg_coord_read(cgfile, cgbase, cgzone, "CoordinateX",
            CGNS_ENUMV(RealSingle), &rmin, &rmax, fbuf))
        cg_error_exit();
    for (np = 0, n = 0; n < nnodes; n++) {
        if (fbuf[n] != xcoord[n]) np++;
    }
    nn += np;
    if (np) printf("%d differences in CoordinateX\n", np);

    if (cg_coord_read(cgfile, cgbase, cgzone, "CoordinateY",
	    CGNS_ENUMV(RealSingle), &rmin, &rmax, fbuf))
        cg_error_exit();
    for (np = 0, n = 0; n < nnodes; n++) {
        if (fbuf[n] != ycoord[n]) np++;
    }
    nn += np;
    if (np) printf("%d differences in CoordinateY\n", np);

    if (cg_coord_read(cgfile, cgbase, cgzone, "CoordinateZ",
	    CGNS_ENUMV(RealSingle), &rmin, &rmax, fbuf))
        cg_error_exit();
    for (np = 0, n = 0; n < nnodes; n++) {
        if (fbuf[n] != zcoord[n]) np++;
    }
    nn += np;
    if (np) printf("%d differences in CoordinateZ\n", np);

    /* check elements */

    if (cg_section_read(cgfile, cgbase, cgzone, 1, name,
            &type, &is, &ie, &k, &n) ||
        cg_elements_read(cgfile, cgbase, cgzone, 1, ibuf, NULL))
        cg_error_exit();
    if (strcmp(name, "Elements") || type != CGNS_ENUMV(HEXA_8) || is != 1 ||
        ie != nelems || k != 0 || n != 0) {
        nn++;
        puts("differences in Elements");
    }
    for (np = 0, n = 0; n < 8*nelems; n++) {
        if (elements[n] != ibuf[n]) np++;
    }
    nn += np;
    if (np) printf("%d differences in Elements connectivity\n", np);

    /* check faces */

    if (cg_section_read(cgfile, cgbase, cgzone, 2, name,
            &type, &is, &ie, &k, &n) ||
        cg_elements_read(cgfile, cgbase, cgzone, 2, ibuf, pbuf))
        cg_error_exit();
    if (strcmp(name, "Faces") || type != CGNS_ENUMV(QUAD_4) || is != (nelems+1) ||
        ie != (nelems+nfaces) || k != 0 || n != 1) {
        nn++;
        puts("differences in Faces");
    }
    for (np = 0, n = 0; n < 4*nfaces; n++) {
        if (faces[n] != ibuf[n]) np++;
    }
    nn += np;
    if (np) printf("%d differences in Faces connectivity\n", np);

    for (np = 0, n = 0; n < 4*nfaces; n++) {
        if (parent[n] != pbuf[n]) np++;
    }
    nn += np;
    if (np) printf("%d differences in Faces parent data\n", np);

    /* check field */

    if (cg_field_read(cgfile, cgbase, cgzone, 1, "Field",
	    CGNS_ENUMV(RealSingle), &rmin, &rmax, fbuf))
        cg_error_exit();
    for (np = 0, n = 0; n < nnodes; n++) {
        if (fbuf[n] != xcoord[n]) np++;
    }
    nn += np;
    if (np) printf("%d differences in Field\n", np);

    if (nn == 0) puts("no diferences");

#if NUM_RANDOM > 0

    /* do random checks */

    puts("doing random reads/checks");

    for (k = 0; k < NUM_RANDOM; k++) {
        rmin = irandom(1, nnodes);
        rmax = irandom(1, nnodes);
        if (rmin > rmax) {
            nr = rmin;
            rmin = rmax;
            rmax = nr;
        }
        printf("coordinates %d -> %d\n", (int)rmin, (int)rmax);
        n = (int)rmin - 1;
        np = 0;
        nn = (int)(rmax - rmin) + 1;
        if (cg_coord_read(cgfile, cgbase, cgzone, "CoordinateX",
	        CGNS_ENUMV(RealSingle), &rmin, &rmax, fbuf))
            cg_error_exit();
        for (i = 0; i < nn; i++) {
            if (fbuf[i] != xcoord[n+i]) np++;
        }

        if (cg_coord_read(cgfile, cgbase, cgzone, "CoordinateY",
	        CGNS_ENUMV(RealSingle), &rmin, &rmax, fbuf))
            cg_error_exit();
        for (i = 0; i < nn; i++) {
            if (fbuf[i] != ycoord[n+i]) np++;
        }

        if (cg_coord_read(cgfile, cgbase, cgzone, "CoordinateZ",
	        CGNS_ENUMV(RealSingle), &rmin, &rmax, fbuf))
            cg_error_exit();
        for (i = 0; i < nn; i++) {
            if (fbuf[i] != zcoord[n+i]) np++;
        }
        if (np) printf("%d differences in Coordinates\n", np);
    }

    for (k = 0; k < NUM_RANDOM; k++) {
        rmin = irandom(1, nelems);
        rmax = irandom(1, nelems);
        if (rmin > rmax) {
            nr = rmin;
            rmin = rmax;
            rmax = nr;
        }
        n = (int)(rmin - 1) << 3;
        nn = (int)(rmax - rmin + 1) << 3;
        printf("elements %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_ElementPartialSize(cgfile, cgbase, cgzone, 1,
                rmin, rmax, &nr) ||
            cg_elements_partial_read(cgfile, cgbase, cgzone, 1,
                rmin, rmax, ibuf, NULL))
            cg_error_exit();
        if (nr != nn) puts("diference in element data size");
        for (np = 0, i = 0; i < nn; i++) {
            if (ibuf[i] != elements[n+i]) np++;
        }
        if (np) printf("%d differences in element connectivity\n", np);
    }

    for (k = 0; k < NUM_RANDOM; k++) {
        rmin = irandom(1, nfaces);
        rmax = irandom(1, nfaces);
        if (rmin > rmax) {
            nr = rmin;
            rmin = rmax;
            rmax = nr;
        }
        n = (int)(rmin - 1) << 2;
        nn = (int)(rmax - rmin + 1) << 2;
        rmin += nelems;
        rmax += nelems;
        get_parent((int)rmin, (int)rmax, nelems, nfaces);
        printf("faces %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_ElementPartialSize(cgfile, cgbase, cgzone, 2,
                rmin, rmax, &nr) ||
            cg_elements_partial_read(cgfile, cgbase, cgzone, 2,
                rmin, rmax, ibuf, pbuf))
            cg_error_exit();
        if (nr != nn) puts("diference in face data size");
        for (np = 0, i = 0; i < nn; i++) {
            if (ibuf[i] != faces[n+i]) np++;
        }
        if (np) printf("%d differences in face connectivity\n", np);
        for (np = 0, i = 0; i < nn; i++) {
            if (pbuf[i] != ptmp[i]) np++;
        }
        if (np) printf("%d differences in face parent data\n", np);
    }

    for (k = 0; k < NUM_RANDOM; k++) {
        rmin = irandom(1, nnodes);
        rmax = irandom(1, nnodes);
        if (rmin > rmax) {
            nr = rmin;
            rmin = rmax;
            rmax = nr;
        }
        n = (int)rmin - 1;
        nn = (int)(rmax - rmin) + 1;
        printf("field %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_field_read(cgfile, cgbase, cgzone, 1, "Field",
	        CGNS_ENUMV(RealSingle), &rmin, &rmax, fbuf))
            cg_error_exit();
        for (np = 0, i = 0; i < nn; i++) {
            if (fbuf[i] != xcoord[n+i]) np++;
        }
        if (np) printf("%d differences in field data\n", np);
    }

#endif

    puts("deleting Elements and Faces and creating Mixed");

    /* delete sections and create as mixed */

    if (cg_goto(cgfile, cgbase, "Zone_t", 1, NULL) ||
        cg_delete_node("Elements") || cg_delete_node("Faces"))
        cg_error_exit();

    if (cg_section_partial_write(cgfile, cgbase, cgzone, "Mixed",
	    CGNS_ENUMV(MIXED), (cgsize_t)1, (cgsize_t)(nelems + nfaces),
	    0, &cgsect))
        cg_error_exit();

    /* create mixed element connectivity */

    nn = (nelems << 3) + nelems + (nfaces << 2) + nfaces;
    ptmp = (cgsize_t *) malloc (2 * nn * sizeof(cgsize_t));

    i = j = n = 0;
    for (nf = 0; nf < nelems; nf++) {
        ptmp[n++] = CGNS_ENUMV(QUAD_4);
        for (k = 0; k < 4; k++)
            ptmp[n++] = faces[j++];
        ptmp[n++] = CGNS_ENUMV(HEXA_8);
        for (k = 0; k < 8; k++)
            ptmp[n++] = elements[i++];
    }
    while (nf++ < nfaces) {
        ptmp[n++] = CGNS_ENUMV(QUAD_4);
        for (k = 0; k < 4; k++)
            ptmp[n++] = faces[j++];
    }

    free (elements);
    elements = ptmp;
    ibuf = elements + nn;

    /* create parent data */

    np = nelems + nfaces;
    nn = np << 2;
    ptmp = (cgsize_t *) malloc (3 * nn * sizeof(cgsize_t));

    for (n = 0; n < nfaces; n++)
        parent[n] <<= 1;

    for (n = 0; n < nn; n++)
        ptmp[n] = 0;
    for (n = 0, j = 0; j < 4; j++) {
        k = j * np;
        for (i = 0; i < nelems; i++) {
            ptmp[k] = parent[n++];
            k += 2;
        }
        while (i++ < nfaces)
            ptmp[k++] = parent[n++];
    }

    free(faces);
    parent = ptmp;
    ptmp = parent + nn;
    pbuf = ptmp + nn;

    rmin = 2 * nelems + 1;
    rmax = np;
    n = mixed_offset((int)rmin, nelems);
    get_parent((int)rmin, (int)rmax, 0, np);

    printf("mixed %d -> %d\n", (int)rmin, (int)rmax);
    if (cg_elements_partial_write(cgfile, cgbase, cgzone,
            cgsect, rmin, rmax, &elements[n]) ||
        cg_parent_data_partial_write(cgfile, cgbase, cgzone,
            cgsect, rmin, rmax, ptmp))
        cg_error_exit();

    printf("mixed %d -> %d (2 at a time)\n", 1, nelems << 1);
    for (i = 0; i < nelems; i++) {
        rmin = (i << 1) + 1;
        rmax = rmin + 1;
        n = mixed_offset((int)rmin, nelems);
        get_parent((int)rmin, (int)rmax, 0, np);
        if (cg_elements_partial_write(cgfile, cgbase, cgzone,
                cgsect, rmin, rmax, &elements[n]) ||
            cg_parent_data_partial_write(cgfile, cgbase, cgzone,
                cgsect, rmin, rmax, ptmp))
            cg_error_exit();
    }

#if NUM_RANDOM > 0

    /* do random writes */

    puts("doing random writes");

    for (k = 0; k < NUM_RANDOM; k++) {
        rmin = irandom(1, np);
        rmax = irandom(1, np);
        if (rmin > rmax) {
            nr = rmin;
            rmin = rmax;
            rmax = nr;
        }
        n = mixed_offset((int)rmin, nelems);
        get_parent((int)rmin, (int)rmax, 0, np);
        printf("mixed %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_elements_partial_write(cgfile, cgbase, cgzone,
                cgsect, rmin, rmax, &elements[n]) ||
            cg_parent_data_partial_write(cgfile, cgbase, cgzone,
                cgsect, rmin, rmax, ptmp))
            cg_error_exit();
    }

#endif

    puts("checking the data");

    if (cg_section_read(cgfile, cgbase, cgzone, cgsect, name,
            &type, &is, &ie, &k, &n) ||
        cg_elements_read(cgfile, cgbase, cgzone, cgsect, ibuf, pbuf))
        cg_error_exit();
    if (strcmp(name, "Mixed") || type != CGNS_ENUMV(MIXED) || is != 1 ||
        ie != np || k != 0 || n != 1) {
        puts("differences in Mixed");
    }
    nn = mixed_offset(np, nelems);
    for (i = 0, n = 0; n < nn; n++) {
        if (elements[n] != ibuf[n]) i++;
    }
    if (i) printf("%d differences in Mixed connectivity\n", i);

    nn = (nelems + nfaces) << 2;
    for (i = 0, n = 0; n < nn; n++) {
        if (parent[n] != pbuf[n]) i++;
    }
    if (i) printf("%d differences in Mixed parent data\n", i);

#if NUM_RANDOM > 0

    /* do random checks */

    puts("doing random reads/checks");

    for (k = 0; k < NUM_RANDOM; k++) {
        rmin = irandom(1, np);
        rmax = irandom(1, np);
        if (rmin > rmax) {
            nr = rmin;
            rmin = rmax;
            rmax = nr;
        }
        n = mixed_offset((int)rmin, nelems);
        nn = mixed_offset((int)(rmax + 1), nelems) - n;
        get_parent((int)rmin, (int)rmax, 0, np);
        printf("mixed %d -> %d\n", (int)rmin, (int)rmax);
        if (cg_ElementPartialSize(cgfile, cgbase, cgzone, cgsect,
                rmin, rmax, &nr) ||
            cg_elements_partial_read(cgfile, cgbase, cgzone, cgsect,
                rmin, rmax, ibuf, pbuf))
            cg_error_exit();
        if (nr != nn) puts("diference in mixed data size");
        for (nf = 0, i = 0; i < nn; i++) {
            if (ibuf[i] != elements[n+i]) nf++;
        }
        if (nf) printf("%d differences in mixed connectivity\n", nf);
        nn = (int)(rmax - rmin + 1) << 2;
        for (nf = 0, i = 0; i < nn; i++) {
            if (pbuf[i] != ptmp[i]) nf++;
        }
        if (nf) printf("%d differences in mixed parent data\n", nf);
    }

#endif

    puts ("closing file");
    cg_close (cgfile);

    return 0;
}

