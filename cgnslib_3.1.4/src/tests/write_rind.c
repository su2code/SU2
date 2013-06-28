#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

int CellDim = 3, PhyDim = 3;
int cgfile, cgbase, cgzone, cggrid, cgsect, cgsol, cgfld;

cgsize_t size[9] = {5, 5, 5,  4, 4, 4,  0, 0, 0};
int rind[6] = {0, 0,  0, 0,  1, 1};

#define NUM_I (size[0] + rind[0] + rind[1])
#define NUM_J (size[1] + rind[2] + rind[3])
#define NUM_K (size[2] + rind[4] + rind[5])

#define INDEX(I,J,K) (int)((I) + NUM_I * ((J) + NUM_J * (K)))

cgsize_t num_coord;
float *xcoord, *ycoord, *zcoord;
cgsize_t *nmap;

cgsize_t num_element, *elements;

float *solution;

static void compute_coord (n, i, j, k)
{
    xcoord[n] = (float)(i - rind[0]);
    ycoord[n] = (float)(j - rind[2]);
    zcoord[n] = (float)(k - rind[4]);
}

static void compute_element (ne, i, j, k)
{
    int n = ne << 3;

    elements[n++] = nmap[INDEX(i,   j,   k)];
    elements[n++] = nmap[INDEX(i+1, j,   k)];
    elements[n++] = nmap[INDEX(i+1, j+1, k)];
    elements[n++] = nmap[INDEX(i,   j+1, k)];
    elements[n++] = nmap[INDEX(i,   j,   k+1)];
    elements[n++] = nmap[INDEX(i+1, j,   k+1)];
    elements[n++] = nmap[INDEX(i+1, j+1, k+1)];
    elements[n]   = nmap[INDEX(i,   j+1, k+1)];
}

int main (int argc, char *argv[])
{
    int n, nn, i, j, k, ni, nj, nk;
    cgsize_t dims[3];
    int grind[2], erind[2];

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

    num_coord = NUM_I * NUM_J * NUM_K;
    num_element = (NUM_I - 1) * (NUM_J - 1) * (NUM_K - 1);
    xcoord = (float *) malloc((size_t)(4 * num_coord * sizeof(float)));
    nmap = (cgsize_t *) malloc((size_t)((num_coord + 8 * num_element) * sizeof(cgsize_t)));
    if (NULL == xcoord || NULL == nmap) {
        fprintf(stderr, "malloc failed for data\n");
        exit(1);
    }
    ycoord = xcoord + num_coord;
    zcoord = ycoord + num_coord;
    solution = zcoord + num_coord;
    elements = nmap + num_coord;

    for (n = 0; n < num_coord; n++)
        solution[n] = (float)n;

    unlink("rind.cgns");
    if (cg_open("rind.cgns", CG_MODE_WRITE, &cgfile))
        cg_error_exit();

    /*--- structured grid with rind ---*/

    printf ("writing structured base with rind\n");
    fflush (stdout);

    for (n = 0, k = 0; k < NUM_K; k++) {
        for (j = 0; j < NUM_J; j++) {
            for (i = 0; i < NUM_I; i++) {
                compute_coord(n++, i, j, k);
            }
        }
    }

    if (cg_base_write(cgfile, "Structured", CellDim, PhyDim, &cgbase) ||
        cg_goto(cgfile, cgbase, "end") ||
        cg_dataclass_write(CGNS_ENUMV( NormalizedByUnknownDimensional )) ||
        cg_zone_write(cgfile, cgbase, "Zone", size, CGNS_ENUMV( Structured ), &cgzone))
        cg_error_exit();

    /* can't use cg_coord_write to write rind coordinates
       need to use cg_grid_write to create the node, cg_goto to set
       position at the node, then write rind and coordinates as an array */

    dims[0] = NUM_I;
    dims[1] = NUM_J;
    dims[2] = NUM_K;

    if (cg_grid_write(cgfile, cgbase, cgzone, "GridCoordinates", &cggrid) ||
        cg_goto(cgfile, cgbase, "Zone_t", cgzone,
            "GridCoordinates_t", cggrid, "end") ||
        cg_rind_write(rind) ||
        cg_array_write("CoordinateX", CGNS_ENUMV( RealSingle ), 3, dims, xcoord) ||
        cg_array_write("CoordinateY", CGNS_ENUMV( RealSingle ), 3, dims, ycoord) ||
        cg_array_write("CoordinateZ", CGNS_ENUMV( RealSingle ), 3, dims, zcoord))
        cg_error_exit();

    /* a similiar technique is used for the solution with rind,
       but we use cg_field_write instead of cg_array_write
       and the solution dimensions come from the zone sizes */

    if (cg_sol_write(cgfile, cgbase, cgzone, "VertexSolution",
		     CGNS_ENUMV( Vertex ), &cgsol) ||
        cg_goto(cgfile, cgbase, "Zone_t", cgzone,
            "FlowSolution_t", cgsol, "end") ||
        cg_rind_write(rind) ||
        cg_field_write(cgfile, cgbase, cgzone, cgsol, CGNS_ENUMV( RealSingle ),
            "Density", solution, &cgfld))
        cg_error_exit();

    /*--- unstructured with rind ---*/

    printf ("writing unstructured base with rind\n");
    fflush (stdout);

    /* rind here has dimension rind[2], so need to put all the
       rind coordinates at the beginning and/or end of the array.
       Just for grins, I'll put some at both ends, although
       it's probably best to put the rind coordinates at the end.
       I'll use the nmap array for building elements */

    ni = (int)(size[0] + rind[0]);
    nj = (int)(size[1] + rind[2]);
    nk = (int)(size[2] + rind[4]);

    for (n = 0, i = 0; i < rind[0]; i++) {
        for (k = 0; k < NUM_K; k++) {
            for (j = 0; j < NUM_J; j++) {
                compute_coord (n, i, j, k);
                nn = INDEX(i, j, k);
                nmap[nn] = ++n;
            }
        }
    }
    for (j = 0; j < rind[2]; j++) {
        for (k = 0; k < NUM_K; k++) {
            for (i = rind[0]; i < ni; i++) {
                compute_coord (n, i, j, k);
                nn = INDEX(i, j, k);
                nmap[nn] = ++n;
            }
        }
    }
    for (k = 0; k < rind[4]; k++) {
        for (j = rind[2]; j < nj; j++) {
            for (i = rind[0]; i < ni; i++) {
                compute_coord (n, i, j, k);
                nn = INDEX(i, j, k);
                nmap[nn] = ++n;
            }
        }
    }
    grind[0] = n;

    for (k = rind[4]; k < nk; k++) {
        for (j = rind[2]; j < nj; j++) {
            for (i = rind[0]; i < ni; i++) {
                compute_coord (n, i, j, k);
                nn = INDEX(i, j, k);
                nmap[nn] = ++n;
            }
        }
    }
    grind[1] = (int)(num_coord - n);

    for (i = ni; i < NUM_I; i++) {
        for (k = 0; k < NUM_K; k++) {
            for (j = 0; j < NUM_J; j++) {
                compute_coord (n, i, j, k);
                nn = INDEX(i, j, k);
                nmap[nn] = ++n;
            }
        }
    }
    for (j = nj; j < NUM_J; j++) {
        for (k = 0; k < NUM_K; k++) {
            for (i = rind[0]; i < ni; i++) {
                compute_coord (n, i, j, k);
                nn = INDEX(i, j, k);
                nmap[nn] = ++n;
            }
        }
    }
    for (k = nk; k < NUM_K; k++) {
        for (j = rind[2]; j < nj; j++) {
            for (i = rind[0]; i < ni; i++) {
                compute_coord (n, i, j, k);
                nn = INDEX(i, j, k);
                nmap[nn] = ++n;
            }
        }
    }

    /* rind elements are like the coordinates, they need to go
       at the beginning and/or end of the element array, although
       at the end is probably best */

    for (n = 0, i = 0; i < rind[0]; i++) {
        for (k = 0; k < NUM_K - 1; k++) {
            for (j = 0; j < NUM_J - 1; j++) {
                compute_element(n++, i, j, k);
            }
        }
    }
    for (j = 0; j < rind[2]; j++) {
        for (k = 0; k < NUM_K - 1; k++) {
            for (i = rind[0]; i < ni - 1; i++) {
                compute_element(n++, i, j, k);
            }
        }
    }
    for (k = 0; k < rind[4]; k++) {
        for (j = rind[2]; j < nj - 1; j++) {
            for (i = rind[0]; i < ni - 1; i++) {
                compute_element(n++, i, j, k);
            }
        }
    }
    erind[0] = n;

    for (k = rind[4]; k < nk - 1; k++) {
        for (j = rind[2]; j < nj - 1; j++) {
            for (i = rind[0]; i < ni - 1; i++) {
                compute_element(n++, i, j, k);
            }
        }
    }
    erind[1] = (int)(num_element - n);

    for (i = ni - 1; i < NUM_I - 1; i++) {
        for (k = 0; k < NUM_K - 1; k++) {
            for (j = 0; j < NUM_J - 1; j++) {
                compute_element(n++, i, j, k);
            }
        }
    }
    for (j = nj - 1; j < NUM_J - 1; j++) {
        for (k = 0; k < NUM_K - 1; k++) {
            for (i = rind[0]; i < ni - 1; i++) {
                compute_element(n++, i, j, k);
            }
        }
    }
    for (k = nk - 1; k < NUM_K - 1; k++) {
        for (j = rind[2]; j < nj - 1; j++) {
            for (i = rind[0]; i < ni - 1; i++) {
                compute_element(n++, i, j, k);
            }
        }
    }

    /* create base, zone and write coordinates.
       As for the structured case, the rind coordinates
       and elements are not included in the zone totals */

    dims[0] = num_coord - grind[0] - grind[1];
    dims[1] = num_element - erind[0] - erind[1];
    dims[2] = 0;

    if (cg_base_write(cgfile, "Unstructured", CellDim, PhyDim, &cgbase) ||
        cg_goto(cgfile, cgbase, "end") ||
        cg_dataclass_write(CGNS_ENUMV( NormalizedByUnknownDimensional )) ||
        cg_zone_write(cgfile, cgbase, "Zone", dims, CGNS_ENUMV( Unstructured ), &cgzone) ||
        cg_grid_write(cgfile, cgbase, cgzone, "GridCoordinates", &cggrid) ||
        cg_goto(cgfile, cgbase, "Zone_t", cgzone,
            "GridCoordinates_t", cggrid, "end") ||
        cg_rind_write(grind) ||
        cg_array_write("CoordinateX", CGNS_ENUMV( RealSingle ), 1, &num_coord, xcoord) ||
        cg_array_write("CoordinateY", CGNS_ENUMV( RealSingle ), 1, &num_coord, ycoord) ||
        cg_array_write("CoordinateZ", CGNS_ENUMV( RealSingle ), 1, &num_coord, zcoord))
        cg_error_exit();

    /* to write the elements with rind elements,
       write all the elements, then use goto to add rind */

    if (cg_section_write(cgfile, cgbase, cgzone, "Elements", CGNS_ENUMV( HEXA_8 ),
            1, num_element, 0, elements, &cgsect) ||
        cg_goto(cgfile, cgbase, "Zone_t", cgzone,
            "Elements_t", cgsect, "end") ||
        cg_rind_write(erind))
        cg_error_exit();

    /* write solution a vertices with rind */

    if (cg_sol_write(cgfile, cgbase, cgzone, "VertexSolution",
		     CGNS_ENUMV( Vertex ), &cgsol) ||
        cg_goto(cgfile, cgbase, "Zone_t", cgzone,
            "FlowSolution_t", cgsol, "end") ||
        cg_rind_write(grind) ||
        cg_field_write(cgfile, cgbase, cgzone, cgsol, CGNS_ENUMV( RealSingle ),
            "Density", solution, &cgfld))
        cg_error_exit();

    cg_close(cgfile);
    return 0;
}

