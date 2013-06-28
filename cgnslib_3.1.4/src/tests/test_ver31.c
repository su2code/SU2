#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#if defined(_WIN32) && !defined(__NUTC__)
# include <io.h>
# define unlink _unlink
#else
# include <unistd.h>
#endif
#include "cgnslib.h"

/* set UNSTRUCTURED_1TO1 to write the unstructured case
   zone connectivities as 1to1 instead of abutting1to1
#define UNSTRUCTURED_1TO1
*/

/* set ABUTTING1TO1_FACES to write the unstructured case
   zone connectivites as abutting1to1 with faces(elements)
   instead of points.
#define ABUTTING1TO1_FACES
*/

int CellDim = 3, PhyDim = 3;

int cgfile, cgbase, cgzone;
int CellDim, PhyDim;
cgsize_t size[9];

#define NUM_ZCONN  5
#define ZCONN_NAME "ZoneConnectivity"

#define NUM_SIDE 5
#define NODE_INDEX(I,J,K) ((I)+NUM_SIDE*(((J)-1)+NUM_SIDE*((K)-1)))
#define CELL_INDEX(I,J,K) ((I)+(NUM_SIDE-1)*(((J)-1)+(NUM_SIDE-1)*((K)-1)))

int num_coord;
float *xcoord, *ycoord, *zcoord;
int num_element, num_face;
cgsize_t *elements, *faces, *parent;

int npts;
cgsize_t *pts, *d_pts;

char errmsg[256];

static float exponents[5] = {0, 1, 0, 0, 0};

void init_data();
void write_structured(), write_unstructured();
void test_zconn(), test_subreg();

void error_exit (char *where)
{
    fprintf (stderr, "ERROR:%s:%s\n", where, cg_get_error());
    exit (1);
}

int main (int argc, char *argv[])
{
    char outfile[32];
    int nbases, nzones;

    strcpy (outfile, "ver31.cgns");
    if (argc > 1) {
        int type = 0;
        int n = 0;
        if (argv[1][n] == '-') n++;
        if (argv[1][n] == 'a' || argv[1][n] == 'A') {
            type = CG_FILE_ADF;
            strcpy (outfile, "ver31.cga");
        }
        else if (argv[1][n] == 'h' || argv[1][n] == 'H') {
            type = CG_FILE_HDF5;
            strcpy (outfile, "ver31.cgh");
        }
        else {
            fprintf(stderr, "unknown option\n");
            exit (1);
        }
        if (cg_set_file_type(type))
            error_exit("cg_set_file_type");
    }
    init_data();
    unlink(outfile);
    if (cg_open(outfile, CG_MODE_WRITE, &cgfile)) error_exit("cg_open");
    write_structured();
    write_unstructured();
    if (cg_close(cgfile)) error_exit("cg_close");

    if (cg_open(outfile, CG_MODE_READ, &cgfile)) error_exit("cg_open");
    if (cg_nbases(cgfile, &nbases)) error_exit("cg_nbases");
    for (cgbase = 1; cgbase <= nbases; cgbase++) {
        if (cg_nzones(cgfile, cgbase, &nzones)) error_exit("cg_nzones");
        for (cgzone = 1; cgzone <= nzones; cgzone++) {
            test_zconn();
            test_subreg();
        }
    }
    if (cg_close(cgfile)) error_exit("cg_close");
        
    return 0;
}

void init_data()
{
    int n, i, j, k, nn, nf, np;

    /* compute coordinates - make it twice as big for use with cylindrical */

    num_coord = NUM_SIDE * NUM_SIDE * NUM_SIDE;
    xcoord = (float *) malloc (6 * num_coord * sizeof(float));
    if (NULL == xcoord) {
        fprintf(stderr, "malloc failed for coordinates\n");
        exit(1);
    }
    ycoord = xcoord + 2 * num_coord;
    zcoord = ycoord + 2 * num_coord;
    for (n = 0, k = 0; k < NUM_SIDE; k++) {
        for (j = 0; j < NUM_SIDE; j++) {
            for (i = 0; i < NUM_SIDE; i++, n++) {
                xcoord[n] = (float)i;
                ycoord[n] = (float)j;
                zcoord[n] = (float)k;
            }
        }
    }

    /* compute elements */

    num_element = (NUM_SIDE - 1) * (NUM_SIDE - 1) * (NUM_SIDE - 1);
    elements = (cgsize_t *) malloc (8 * num_element * sizeof(cgsize_t));
    if (NULL == elements) {
        fprintf(stderr, "malloc failed for elements");
        exit(1);
    }
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

    /* compute outside face elements */

    num_face = 6 * (NUM_SIDE - 1) * (NUM_SIDE - 1);
    faces = (cgsize_t *) malloc (4 * num_face * sizeof(cgsize_t));
    parent = (cgsize_t *) malloc (4 * num_face * sizeof(cgsize_t));
    if (NULL == faces || NULL == parent) {
        fprintf(stderr, "malloc failed for elements");
        exit(1);
    }
    for (n = 0; n < 4*num_face; n++)
        parent[n] = 0;
    nf = np = 0;
    n = 2 * num_face;
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

    /* connectivity points - make it big enough to hold 4 surfaces */

    npts = NUM_SIDE * NUM_SIDE;
    pts = (cgsize_t *) malloc (12 * npts * sizeof(cgsize_t));
    if (NULL == pts) {
        fprintf(stderr, "malloc failed for connectivity points");
        exit(1);
    }
    d_pts = pts + 6 * npts;
}

void write_coords(int nz)
{
    int k, nn, n, nij, koff, cgcoord;

    koff = nz == 1 ? 1 - NUM_SIDE : 0;
    nij = NUM_SIDE * NUM_SIDE;
    for (n = 0, k = 0; k < NUM_SIDE; k++) {
        for (nn = 0; nn < nij; nn++)
            zcoord[n++] = (float)(k + koff);
    }

    if (cg_coord_write(cgfile, cgbase, nz, CGNS_ENUMV(RealSingle),
            "CoordinateX", xcoord, &cgcoord) ||
        cg_coord_write(cgfile, cgbase, nz, CGNS_ENUMV(RealSingle),
            "CoordinateY", ycoord, &cgcoord) ||
        cg_coord_write(cgfile, cgbase, nz, CGNS_ENUMV(RealSingle),
            "CoordinateZ", zcoord, &cgcoord)) {
        sprintf (errmsg, "zone %d coordinates", nz);
        error_exit(errmsg);
    }
    for (n = 1; n <= 3; n++) {
        if (cg_goto(cgfile, cgbase, "Zone_t", nz, "GridCoordinates_t", 1,
                "DataArray_t", n, NULL) ||
            cg_exponents_write(CGNS_ENUMV(RealSingle), exponents))
            error_exit("coordinate exponents");
    }
}

void write_elements(int nz)
{
    int cgsect;

    if (cg_section_write(cgfile, cgbase, nz, "Elements", CGNS_ENUMV(HEXA_8),
            1, num_element, 0, elements, &cgsect) ||
        cg_section_write(cgfile, cgbase, nz, "Faces", CGNS_ENUMV(QUAD_4),
            num_element+1, num_element+num_face, 0, faces, &cgsect) ||
        cg_parent_data_write(cgfile, cgbase, nz, cgsect, parent)) {
        sprintf (errmsg, "zone %d elements", nz);
        error_exit(errmsg);
    }
}

/*------------------------------------------------------------------------
 * structured grid
 *------------------------------------------------------------------------*/

void write_structured()
{
    int i, j, n, nc, cgfam, cgz, cgbc, cgconn, cghole, cgsr;
    cgsize_t range[6], d_range[6], dims[2];
    int transform[3], rind[6], iter[NUM_ZCONN];
    char name[33], cname[33], hname[33], sname[33], zcname[33];
    char cpointers[32*NUM_ZCONN+1], spointers[32*NUM_ZCONN+1];

    if (cg_base_write(cgfile, "Structured", CellDim, PhyDim, &cgbase) ||
        cg_goto(cgfile, cgbase, "end") ||
        cg_descriptor_write("Descriptor", "Multi-block Structured Grid") ||
        cg_dataclass_write(CGNS_ENUMV(Dimensional)) ||
        cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter),
            CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)))
        error_exit("structured base");
    if (cg_simulation_type_write(cgfile, cgbase, CGNS_ENUMV(TimeAccurate)))
        error_exit("simulation type");

    /* write a FamilyBCDataset node and children */

    if (cg_family_write(cgfile, cgbase, "Family", &cgfam) ||
        cg_fambc_write(cgfile, cgbase, cgfam, "Wall",
            CGNS_ENUMV(BCWall), &cgbc) ||
        cg_goto(cgfile, cgbase, "Family_t", 1, "FamilyBC_t", 1, NULL) ||
        cg_bcdataset_write("Dataset", CGNS_ENUMV(BCWall),
            CGNS_ENUMV(Dirichlet)) ||
        cg_gorel(cgfile, "Dataset", 0, NULL) ||
/*        cg_goto(cgfile, cgbase, "Family_t", 1, "FamilyBC_t", 1,
            "FamilyBCDataSet_t", 1, NULL) ||*/
        cg_descriptor_write("descriptor", "this is a descriptor") ||
        cg_dataclass_write(CGNS_ENUMV(NormalizedByDimensional)) ||
        cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter),
            CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)) ||
        cg_state_write("reference state") ||
        cg_user_data_write("Userdata") ||
        cg_gorel(cgfile, "UserDefinedData_t", 1, NULL) ||
        cg_famname_write("Family"))
        error_exit("family bc");

    /* write zones */

    for (n = 0; n < 3; n++) {
        size[n]   = NUM_SIDE;
        size[n+3] = NUM_SIDE - 1;
        size[n+6] = 0;
    }
    for (n = 1; n <= 2; n++) {
        sprintf(name, "Zone%d", n);
        if (cg_zone_write(cgfile, cgbase, name, size,
                CGNS_ENUMV(Structured), &cgzone)) {
            sprintf (errmsg, "structured zone %d", n);
            error_exit(errmsg);
        }
        write_coords(n);
    }

    /* write inlet BC (zone 1) as point range */

    for (n = 0; n < 3; n++) {
        range[n] = 1;
        range[n+3] = NUM_SIDE;
    }
    range[5] = 1;
    if (cg_boco_write(cgfile, cgbase, 1, "Inlet", CGNS_ENUMV(BCInflow),
            CGNS_ENUMV(PointRange), 2, range, &cgbc))
        error_exit("inlet boco");

    /* write outlet BC (zone 2) as point list */

    for (n = 0, j = 1; j <= NUM_SIDE; j++) {
        for (i = 1; i <= NUM_SIDE; i++) {
            pts[n++] = i;
            pts[n++] = j;
            pts[n++] = NUM_SIDE;
        }
    }
    if (cg_boco_write(cgfile, cgbase, 2, "Outlet", CGNS_ENUMV(BCOutflow),
            CGNS_ENUMV(PointList), npts, pts, &cgbc))
        error_exit("outlet boco");

    /* write connectivities in multiple ZoneGridConnectivity_t nodes */
    
    for (nc = 1; nc <= NUM_ZCONN; nc++) {
        /* create ZoneGridConnectivity_t node */
        sprintf(zcname, "%s%d", ZCONN_NAME, nc);
        if (cg_zconn_write(cgfile, cgbase, 1, zcname, &cgz) ||
            cg_zconn_write(cgfile, cgbase, 2, zcname, &cgz))
            error_exit(name);
        sprintf(cname, "conn%d", nc);
        sprintf(hname, "hole%d", nc);

        /* write zone 1 to zone 2 connectivity as 1to1 */

        for (n = 0; n < 3; n++) {
            range[n] = d_range[n] = 1;
            range[n+3] = d_range[n+3] = NUM_SIDE;
            transform[n] = n + 1;
        }
        range[2] = NUM_SIDE;
        d_range[5] = 1;
        if (cg_1to1_write(cgfile, cgbase, 1, cname, "Zone2",
                range, d_range, transform, &cgconn)) {
            sprintf(errmsg, "Zone1/%s/%s", zcname, cname);
            error_exit(errmsg);
        }

        /* write zone 2 to zone 1 connectivity as Abutting1to1 */

        for (n = 0; n < 3; n++) {
            range[n] = 1;
            range[n+3] = NUM_SIDE;
        }
        range[5] = 1;
        for (n = 0, j = 1; j <= NUM_SIDE; j++) {
            for (i = 1; i <= NUM_SIDE; i++) {
                d_pts[n++] = i;
                d_pts[n++] = j;
                d_pts[n++] = 1;
            }
        }
        if (cg_conn_write(cgfile, cgbase, 2, cname,
                CGNS_ENUMV(Vertex), CGNS_ENUMV(Abutting1to1),
                CGNS_ENUMV(PointRange), 2, range, "Zone1",
                CGNS_ENUMV(Structured), CGNS_ENUMV(PointListDonor),
                CGNS_ENUMV(Integer), npts, d_pts, &cgconn)) {
            sprintf(errmsg, "Zone2/%s/%s", zcname, cname);
            error_exit(errmsg);
        }

        /* write hole in zone 1 as PointList */
        
        if (cg_hole_write(cgfile, cgbase, 1, hname, CGNS_ENUMV(Vertex),
                CGNS_ENUMV(PointList), 1, npts, d_pts, &cghole)) {
            sprintf(errmsg, "Zone1/%s/%s", zcname, hname);
            error_exit(errmsg);
        }
        
        /* write hole in zone 2 as PointRange */
        
        if (cg_hole_write(cgfile, cgbase, 2, hname, CGNS_ENUMV(Vertex),
                CGNS_ENUMV(PointRange), 1, 2, range, &cghole)) {
            sprintf(errmsg, "Zone2/%s/%s", zcname, hname);
            error_exit(errmsg);
        }
    }

    /* write SubRegion 1 as BCRegionName */

    if (cg_subreg_bcname_write(cgfile, cgbase, 1, "SubRegion1",
            2, "Inlet", &cgsr))
        error_exit("cg_subreg_bcname_write(Inlet)");
    if (cg_goto(cgfile, cgbase, "Zone_t", 1, "ZoneSubRegion_t", 1, NULL))
        error_exit("cg_goto ZoneSubRegion_t 1");
    if (cg_dataclass_write(CGNS_ENUMV(Dimensional)))
        error_exit("cg_dataclass_write Dimensional");
    if (cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter),
            CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Degree)))
        error_exit("cg_units_write");
    if (cg_descriptor_write("descriptor", "this is zone 1, subregion 1"))
        error_exit("cg_descriptor_write");

    if (cg_subreg_bcname_write(cgfile, cgbase, 2, "SubRegion1",
            2, "Outlet", &cgsr))
        error_exit("cg_subreg_bcname_write(Outlet)");
    if (cg_goto(cgfile, cgbase, "Zone_t", 2, "SubRegion1", 0, NULL))
        error_exit("cg_goto SubRegion1");
    if (cg_famname_write("Family"))
        error_exit("cg_famname_write SubRegion1");
    if (cg_user_data_write("Userdata"))
        error_exit("cg_user_data_write");
    if (cg_gorel(cgfile, "UserDefinedData_t", 1, NULL))
        error_exit("cg_gorel UserDefinedData");
    if (cg_famname_write("Family"))
        error_exit("cg_famname_write SubRegion1/UserDefinedData");

    /* write SubRegion 2 as ConnectivityRegionName */

    sprintf(cname, "%s2/conn2", ZCONN_NAME);
    strcpy(sname, "SubRegion2");
    if (cg_subreg_gcname_write(cgfile, cgbase, 1, sname, 2, cname, &cgsr) ||
        cg_subreg_gcname_write(cgfile, cgbase, 2, sname, 2, cname, &cgsr)) {
        sprintf(errmsg, "cg_subreg_gcname_write(%s)", sname);
        error_exit(errmsg);
    }

    /* write subregion as PointRange in Zone1 and PointList in Zone2 */

    for (n = 0; n < 3; n++) {
        range[n] = 1;
        range[n+3] = NUM_SIDE;
    }
    range[5] = 1;
    for (n = 0, j = 1; j <= NUM_SIDE; j++) {
        for (i = 1; i <= NUM_SIDE; i++) {
            pts[n++] = i;
            pts[n++] = j;
            pts[n++] = NUM_SIDE;
        }
    }
    for (n = 0; n < 6; n++)
        rind[n] = 1;
    dims[0] = NUM_SIDE * NUM_SIDE;

    for (nc = 3; nc <= NUM_ZCONN; nc++) {
        sprintf(sname, "SubRegion%d", nc);
        if (cg_subreg_ptset_write(cgfile, cgbase, 1, sname, 2,
                CGNS_ENUMV(Vertex), CGNS_ENUMV(PointRange),
                2, range, &cgsr)) {
            sprintf(errmsg, "cg_subreg_ptset_write(%s PointRange)", sname);
            error_exit(errmsg);
        }
        if (cg_goto(cgfile, cgbase, "Zone_t", 1, sname, 0, NULL)) {
            sprintf(errmsg, "cg_goto %s", sname);
            error_exit(errmsg);
        }
        if (cg_rind_write(rind)) {
            sprintf(errmsg, "cg_rind_write %s", sname);
            error_exit(errmsg);
        }
        if (cg_array_write("DataArray", CGNS_ENUMV(RealSingle), 1,
                dims, xcoord)) {
            sprintf(errmsg, "cg_array_write %s", sname);
            error_exit(errmsg);
        }
        if (cg_gorel(cgfile, "DataArray", 0, NULL) ||
            cg_exponents_write(CGNS_ENUMV(RealSingle), exponents)) {
            sprintf(errmsg, "cg_exponents_write %s", sname);
            error_exit(errmsg);
        }
        if (cg_subreg_ptset_write(cgfile, cgbase, 2, sname, 2,
                CGNS_ENUMV(EdgeCenter), CGNS_ENUMV(PointList),
                npts, pts, &cgsr)) {
            sprintf(errmsg, "cg_subreg_ptset_write(%s PointList)", sname);
            error_exit(errmsg);
        }
        if (cg_goto(cgfile, cgbase, "Zone_t", 2, sname, 0, NULL)) {
            sprintf(errmsg, "cg_goto %s", sname);
            error_exit(errmsg);
        }
        if (cg_array_write("DataArray", CGNS_ENUMV(RealSingle), 1,
                dims, xcoord)) {
            sprintf(errmsg, "cg_array_write %s", sname);
            error_exit(errmsg);
        }
        if (cg_gorel(cgfile, "DataArray", 0, NULL) ||
            cg_exponents_write(CGNS_ENUMV(RealSingle), exponents)) {
            sprintf(errmsg, "cg_exponents_write %s", sname);
            error_exit(errmsg);
        }
    }

    /* create BaseIterativeData_t node */

    for (n = 0; n < NUM_ZCONN; n++)
        iter[n] = n + 1;
    dims[0] = NUM_ZCONN;

    if (cg_biter_write(cgfile, cgbase, "BaseIterativeData", NUM_ZCONN))
        error_exit("cg_biter_write");
    if (cg_goto(cgfile, cgbase, "BaseIterativeData", 0, NULL))
        error_exit("cg_goto BaseIterativeData");
    if (cg_array_write("IterationValues", CGNS_ENUMV(Integer),
        1, dims, iter)) error_exit("cg_array_write IterationValues");

    /* create the ZoneIterativeData_t node */

    dims[0] = 32;
    dims[1] = NUM_ZCONN;
    n = 32 * NUM_ZCONN;
    for (i = 0; i < n; i++) {
        cpointers[i] = ' ';
        spointers[i] = ' ';
    }
    cpointers[n] = 0;
    spointers[n] = 0;

    for (n = 0, nc = 1; nc <= NUM_ZCONN; nc++, n+=32) {
        sprintf(zcname, "%s%d", ZCONN_NAME, nc);
        i = strlen(zcname);
        strcpy(&cpointers[n], zcname);
        cpointers[n+i] = ' ';
        sprintf(sname, "SubRegion%d", nc);
        i = strlen(sname);
        strcpy(&spointers[n], sname);
        spointers[n+i] = ' ';
    }

    for (n = 1; n <= 2; n++) {
        sprintf(name, "Zone%dInterativeData", n);
        if (cg_ziter_write(cgfile, cgbase, n, name)) {
            sprintf(errmsg, "%s:cg_ziter_write", name);
            error_exit(errmsg);
        }
        if (cg_goto(cgfile, cgbase, "Zone_t", n, name, 0, NULL)) {
            sprintf(errmsg, "%s:cg_goto", name);
            error_exit(errmsg);
        }
        if (cg_array_write("ZoneGridConnectivityPointers",
                CGNS_ENUMV(Character), 2, dims, cpointers)) {
            sprintf(errmsg, "%s:cg_array_write ZoneGridConnectivityPointers",
                name);
            error_exit(errmsg);
        }
        if (cg_array_write("ZoneSubRegionPointers",
                CGNS_ENUMV(Character), 2, dims, spointers)) {
            sprintf(errmsg, "%s:cg_array_write ZoneSubRegionPointers",
                name);
            error_exit(errmsg);
        }
    }
}

/*------------------------------------------------------------------------
 * unstructured grid
 *------------------------------------------------------------------------*/

void write_unstructured()
{
    int n, nelem, nc, cgconn, cgz, cghole;
    int iter[NUM_ZCONN];
    cgsize_t range[2], dims[2];
    char name[33], zcname[33], pointers[32*NUM_ZCONN+1];
#ifdef UNSTRUCTURED_1TO1
    int d_range[2], transform;
#else
# ifdef ABUTTING1TO1_FACES
    GridLocation_t location;
# else
    int i, j;
# endif
#endif

    if (cg_base_write(cgfile, "Unstructured", CellDim, PhyDim, &cgbase) ||
        cg_goto(cgfile, cgbase, "end") ||
        cg_descriptor_write("Descriptor", "Multi-block Unstructured Grid") ||
        cg_dataclass_write(CGNS_ENUMV(NormalizedByDimensional)) ||
        cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter),
            CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)))
        error_exit("unstructured base");

    /* write zones */

    for (n = 0; n < 9; n++)
        size[n] = 0;
    size[0] = num_coord;
    size[1] = num_element;
    for (n = 1; n <= 2; n++) {
        sprintf(name, "Zone%d", n);
        if (cg_zone_write(cgfile, cgbase, name, size,
                CGNS_ENUMV(Unstructured), &cgzone)) {
            sprintf (errmsg, "unstructured zone %d", n);
            error_exit(errmsg);
        }
        write_coords(n);
        write_elements(n);
    }
    nelem = (NUM_SIDE - 1) * (NUM_SIDE - 1);

    /* write connectivities in multiple ZoneGridConnectivity_t nodes */

    for (nc = 1; nc <= NUM_ZCONN; nc++) {
        /* create ZoneGridConnectivity_t node */
        sprintf(zcname, "%s%d", ZCONN_NAME, nc);
        if (cg_zconn_write(cgfile, cgbase, 1, zcname, &cgz) ||
            cg_zconn_write(cgfile, cgbase, 2, zcname, &cgz))
            error_exit(name);

        sprintf(name, "conn%d", nc);

#ifdef UNSTRUCTURED_1TO1

        /* write connectivities as 1to1 */

        range[0] = NODE_INDEX(1, 1, NUM_SIDE);
        range[1] = NODE_INDEX(NUM_SIDE, NUM_SIDE, NUM_SIDE);
        d_range[0] = NODE_INDEX(1, 1, 1);
        d_range[1] = NODE_INDEX(NUM_SIDE, NUM_SIDE, 1);
        transform = 1;
        if (cg_1to1_write(cgfile, cgbase, 1, name, "Zone2",
                range, d_range, &transform, &cgconn)) {
            sprintf(errmsg, "Zone1/%s/%s", name);
            error_exit(errmsg);
        }

        if (cg_1to1_write(cgfile, cgbase, 2, name, "Zone1",
                d_range, range, &transform, &cgconn)) {
            sprintf(errmsg, "Zone2/%s/%s", zcname, name);
            error_exit(errmsg);
        }

#else
# ifdef ABUTTING1TO1_FACES

        /* zone 1 to zone 2 connectivity as Abutting1to1 with element range */

        range[0] = num_element + num_face - nelem + 1;
        range[1] = num_element + num_face;
        for (n = 0; n < nelem; n++)
            d_pts[n] = range[0] + n - nelem;
        location = FaceCenter;

        if (cg_conn_write(cgfile, cgbase, 1, name,
                location, CGNS_ENUMV(Abutting1to1),
                CGNS_ENUMV(PointRange), 2, range, "Zone2",
                CGNS_ENUMV(Unstructured), CGNS_ENUMV(PointListDonor),
                CGNS_ENUMV(Integer), nelem, d_pts, &cgconn)) {
            sprintf(errmsg, "Zone1/%s/%s", zcname, name);
            error_exit(errmsg);
        }

        /* zone 2 to zone 1 connectivity as Abutting1to1 with element list */

        for (n = 0; n < nelem; n++) {
            pts[n]   = num_element + num_face - 2 * nelem + 1 + n;
            d_pts[n] = pts[n] + nelem;
        }
        if (cg_conn_write(cgfile, cgbase, 2, name,
                location, CGNS_ENUMV(Abutting1to1),
                CGNS_ENUMV(PointList), nelem, pts, "Zone1",
                CGNS_ENUMV(Unstructured), CGNS_ENUMV(PointListDonor),
                CGNS_ENUMV(Integer), nelem, d_pts, &cgconn)) {
            sprintf(errmsg, "Zone2/%s/%s", zcname, name);
            error_exit(errmsg);
        }

# else

        /* zone 1 to zone 2 connectivity as Abutting1to1 with point range */

        range[0] = NODE_INDEX(1, 1, NUM_SIDE);
        range[1] = NODE_INDEX(NUM_SIDE, NUM_SIDE, NUM_SIDE);
        for (n = 0, j = 1; j <= NUM_SIDE; j++) {
            for (i = 1; i <= NUM_SIDE; i++)
                d_pts[n++] = NODE_INDEX(i, j, 1);
        }
        if (cg_conn_write(cgfile, cgbase, 1, name,
                CGNS_ENUMV(Vertex), CGNS_ENUMV(Abutting1to1),
                CGNS_ENUMV(PointRange), 2, range, "Zone2",
                CGNS_ENUMV(Unstructured), CGNS_ENUMV(PointListDonor),
                CGNS_ENUMV(Integer), npts, d_pts, &cgconn)) {
            sprintf(errmsg, "Zone1/%s/%s", zcname, name);
            error_exit(errmsg);
        }

        /* zone 2 to zone 1 connectivity as Abutting1to1 with point list */

        for (n = 0, j = 1; j <= NUM_SIDE; j++) {
            for (i = 1; i <= NUM_SIDE; i++) {
                pts[n] = d_pts[n];
                d_pts[n++] = NODE_INDEX(i, j, NUM_SIDE);
            }
        }
        if (cg_conn_write(cgfile, cgbase, 2, name,
                CGNS_ENUMV(Vertex), CGNS_ENUMV(Abutting1to1),
                CGNS_ENUMV(PointList), npts, pts, "Zone1",
                CGNS_ENUMV(Unstructured), CGNS_ENUMV(PointListDonor),
                CGNS_ENUMV(Integer), npts, d_pts, &cgconn)) {
            sprintf(errmsg, "Zone2/%s/%s", zcname, name);
            error_exit(errmsg);
        }
        
        sprintf(name, "hole%d", nc);

        /* write hole in zone 1 as PointRange */

        if (cg_hole_write(cgfile, cgbase, 1, name, CGNS_ENUMV(Vertex),
                CGNS_ENUMV(PointRange), 1, 2, range, &cghole)) {
            sprintf(errmsg, "Zone1/%s/%s", zcname, name);
            error_exit(errmsg);
        }
        
        /* write hole in zone 2 as PointList */
        
        if (cg_hole_write(cgfile, cgbase, 2, name, CGNS_ENUMV(Vertex),
                CGNS_ENUMV(PointList), 1, npts, pts, &cghole)) {
            sprintf(errmsg, "Zone2/%s/%s", zcname, name);
            error_exit(errmsg);
        }

# endif
#endif
    }

    /* create BaseIterativeData_t node */
    

    for (n = 0; n < NUM_ZCONN; n++)
        iter[n] = n + 1;
    dims[0] = NUM_ZCONN;

    if (cg_biter_write(cgfile, cgbase, "BaseIterativeData", NUM_ZCONN))
        error_exit("cg_biter_write");
    if (cg_goto(cgfile, cgbase, "BaseIterativeData", 0, NULL))
        error_exit("cg_goto BaseIterativeData");
    if (cg_array_write("IterationValues", CGNS_ENUMV(Integer),
        1, dims, iter)) error_exit("cg_array_write IterationValues");
    
    /* create the ZoneIterativeData_t node */
    
    dims[0] = 32;
    dims[1] = NUM_ZCONN;
    n = 32 * NUM_ZCONN;
    for (i = 0; i < n; i++)
        pointers[i] = ' ';
    pointers[n] = 0;

    for (n = 0, nc = 1; nc <= NUM_ZCONN; nc++, n+=32) {
        sprintf(zcname, "%s%d", ZCONN_NAME, nc);
        i = strlen(zcname);
        strcpy(&pointers[n], zcname);
        pointers[n+i] = ' ';
    }

    for (n = 1; n <= 2; n++) {
        sprintf(name, "Zone%dInterativeData", n);
        if (cg_ziter_write(cgfile, cgbase, n, name)) {
            sprintf(errmsg, "%s:cg_ziter_write", name);
            error_exit(errmsg);
        }
        if (cg_goto(cgfile, cgbase, "Zone_t", n, name, 0, NULL)) {
            sprintf(errmsg, "%s:cg_goto", name);
            error_exit(errmsg);
        }
        if (cg_array_write("ZoneGridConnectivityPointers",
                CGNS_ENUMV(Character), 2, dims, pointers)) {
            sprintf(errmsg, "%s:cg_array_write", name);
            error_exit(errmsg);
        }
    }
}

/*------------------------------------------------------------------------
 * test Zone connectivity data
 *------------------------------------------------------------------------*/

void set_zconn(const char *zcname)
{
    int n;
    char name[33];

    for (n = 1; n <= NUM_ZCONN; n++) {
        if (cg_zconn_read(cgfile, cgbase, cgzone, n, name)) {
            sprintf(errmsg, "%s:cg_zconn_read", zcname);
            error_exit(errmsg);
        }
        if (0 == strcmp(name, zcname)) {
            if (cg_zconn_set(cgfile, cgbase, cgzone, n)) {
                sprintf(errmsg, "%s:cg_zconn_set", zcname);
                error_exit(errmsg);
            }
            return;
        }
    }
    sprintf(errmsg, "zconn %s not found", zcname);
    error_exit(errmsg);
}

void test_zconn()
{
    int n, nz, nc;
    int nconn, n1to1, nholes;
    char name[33], zcname[33], expected[33], pointers[32*NUM_ZCONN+1];
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(PointSetType_t) ptype;
    cgsize_t np;

    if (cg_nzconns(cgfile, cgbase, cgzone, &nz) || nz != NUM_ZCONN)
        error_exit("cg_nzconns");
    if (cg_ziter_read(cgfile, cgbase, cgzone, name))
        error_exit("cg_ziter_read");
    if (cg_goto(cgfile, cgbase, "Zone_t", cgzone, name, 0, NULL))
        error_exit("cg_goto ziter");
    if (cg_array_read(1, pointers)) error_exit("cg_array_read");
    
    for (nz = 0; nz < NUM_ZCONN; nz++) {
        strncpy(zcname, &pointers[32*nz], 32);
        zcname[32] = 0;
        for (n = strlen(zcname)-1; n >= 0 && zcname[n] == ' '; n--)
            ;
        zcname[++n] = 0;
        set_zconn(zcname);
        n = strlen(ZCONN_NAME);
        nc = atoi(&zcname[n]);
        
        if (cg_nconns(cgfile, cgbase, cgzone, &nconn)) {
            sprintf(errmsg, "%s:cg_nconns", zcname);
            error_exit(errmsg);
        }
        if (cg_n1to1(cgfile, cgbase, cgzone, &n1to1)) {
            sprintf(errmsg, "%s:cg_n1to1", zcname);
            error_exit(errmsg);
        }
        if ((nconn + n1to1) != 1) {
            sprintf(errmsg, "%s:nconns=%d", zcname, (nconn + n1to1));
            error_exit(errmsg);
        }
        if (n1to1) {
            cgsize_t range[6], d_range[6];
            int transform[3];
            if (cg_1to1_read(cgfile, cgbase, cgzone, 1, name, expected,
                    range, d_range, transform)) {
                sprintf(errmsg, "%s:cg_1to1_read", zcname);
                error_exit(errmsg);
            }
        }
        else {
            CGNS_ENUMT(GridConnectivityType_t) type;
            CGNS_ENUMT(PointSetType_t) d_ptype;
            CGNS_ENUMT(ZoneType_t) d_ztype;
            CGNS_ENUMT(DataType_t) d_dtype;
            cgsize_t d_npnts;
           if (cg_conn_info(cgfile, cgbase, cgzone, 1, name, &location,
                    &type, &ptype, &np, expected, &d_ztype, &d_ptype,
                    &d_dtype, &d_npnts)) {
                sprintf(errmsg, "%s:cg_conn_info", zcname);
                error_exit(errmsg);
            }
        }
        sprintf(expected, "conn%d", nc);
        if (strcmp(expected, name)) {
            sprintf(errmsg, "%s:conn name %s != %s", zcname, name, expected);
            error_exit(errmsg);
        }

        if (cg_nholes(cgfile, cgbase, cgzone, &nholes)) {
            sprintf(errmsg, "%s:cg_zconn_set", zcname);
            error_exit(errmsg);
        }
        if (nholes != 1) {
            sprintf(errmsg, "%s:nholes=%d", zcname, nholes);
            error_exit(errmsg);
        }
        if (cg_hole_info(cgfile, cgbase, cgzone, 1, name, &location,
                &ptype, &n, &np)) {
            sprintf(errmsg, "%s:cg_hole_info", zcname);
            error_exit(errmsg);
        }
        sprintf(expected, "hole%d", nc);
        if (strcmp(expected, name)) {
            sprintf(errmsg, "%s:hole name %s != %s", zcname, name, expected);
            error_exit(errmsg);
        }
    }
}

/*------------------------------------------------------------------------
 * test sub region data
 *------------------------------------------------------------------------*/

void test_subreg()
{
    int ns, dim, rind[6];
    int bclen, gclen, ndescr;
    char srname[33], name[33], text[65], *descr;
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(PointSetType_t) ptype;
    CGNS_ENUMT(DataClass_t) dclass;
    CGNS_ENUMT(MassUnits_t) mass;
    CGNS_ENUMT(LengthUnits_t) length;
    CGNS_ENUMT(TimeUnits_t) time;
    CGNS_ENUMT(TemperatureUnits_t) temp;
    CGNS_ENUMT(AngleUnits_t) ang;
    cgsize_t np;
    
    if (cg_nsubregs(cgfile, cgbase, cgzone, &ns))
        error_exit("cg_nsubregs");
    if (ns == 0 && cgbase == 2) return;
    if (ns != NUM_ZCONN) error_exit("invalid nsubregs");

    for (ns = 1; ns < NUM_ZCONN; ns++) {
        if (cg_subreg_info(cgfile, cgbase, cgzone, ns, srname, &dim,
                &location, &ptype, &np, &bclen, &gclen))
            error_exit("cg_subreg_info");
        if (dim != 2) error_exit("dim not 2");
        
        if (bclen) {
            if (gclen || ptype != CGNS_ENUMV(PointSetTypeNull) || np)
                error_exit("bclen");
            if (strcmp(srname, "SubRegion1")) error_exit("SubRegion1");
            /* these should be errors */
            if (!cg_subreg_gcname_read(cgfile, cgbase, cgzone, ns, text) ||
                !cg_subreg_ptset_read(cgfile, cgbase, cgzone, ns, pts))
                error_exit("cg_subreg_bcname_read:multiple");
            /* valid */
            if (cg_subreg_bcname_read(cgfile, cgbase, cgzone, ns, text))
                error_exit("gc_subreg_bcname_read");

            if (cgzone == 1) {
                if (strcmp(text, "Inlet")) error_exit("Inlet");
                if (cg_goto(cgfile, cgbase, "Zone_t", 1, "ZoneSubRegion_t", ns, NULL))
                    error_exit("cg_goto ZoneSubRegion_t");
                if (cg_dataclass_read(&dclass) ||
                    cg_units_read(&mass, &length, &time, &temp, &ang) ||
                    cg_ndescriptors(&ndescr))
                    error_exit("dataclass,units or ndescriptors");
                if (ndescr != 1) error_exit("ndescriptors not 1");
                if (cg_descriptor_read(1, name, &descr))
                    error_exit("cg_descriptor_read");
                if (strcmp(name, "descriptor") ||
                    strcmp(descr, "this is zone 1, subregion 1"))
                    error_exit("descriptor");
                cg_free(descr);
            }
            else {
                if (strcmp(text, "Outlet")) error_exit("Outlet");
                if (cg_goto(cgfile, cgbase, "Zone_t", 2, srname, 0, NULL))
                    error_exit("cg_goto SubRegion1");
                if (cg_famname_read(name) || cg_user_data_read(1, text))
                    error_exit("famname or user data");
                if (strcmp(name, "Family")) error_exit("Family");
                if (strcmp(text, "Userdata")) error_exit("Userdata");
                if (cg_gorel(cgfile, "UserDefinedData_t", 1, NULL))
                    error_exit("cg_gorel SubRegion1/UserDefinedData");
                if (cg_famname_read(name) || strcmp(name, "Family"))
                    error_exit("Userdata/Family");
            }
        }
        
        else if (gclen) {
            if (bclen || ptype != CGNS_ENUMV(PointSetTypeNull) || np)
                error_exit("gclen");
            if (strcmp(srname, "SubRegion2")) error_exit("SubRegion2");
            /* these should be errors */
            if (!cg_subreg_bcname_read(cgfile, cgbase, cgzone, ns, text) ||
                !cg_subreg_ptset_read(cgfile, cgbase, cgzone, ns, pts))
                error_exit("cg_subreg_bcname_read:multiple");
            /* valid */
            if (cg_subreg_gcname_read(cgfile, cgbase, cgzone, ns, text))
                error_exit("cg_subreg_gcname_read");
            if (strcmp(text, "ZoneConnectivity2/conn2"))
                error_exit("ZoneConnectivity2/conn2");
        }

        else {
            if (bclen || gclen) error_exit("ptset");
            if (strcmp(srname, "SubRegion3") &&
                strcmp(srname, "SubRegion4") &&
                strcmp(srname, "SubRegion5")) error_exit("SubRegion[345]");
            /* these should be errors */
            if (!cg_subreg_bcname_read(cgfile, cgbase, cgzone, ns, text) ||
                !cg_subreg_gcname_read(cgfile, cgbase, cgzone, ns, text))
                error_exit("cg_subreg_ptset_read:multiple");

            if (cgzone == 1) {
                if (location != CGNS_ENUMV(Vertex) ||
                    ptype != CGNS_ENUMV(PointRange) || np != 2)
                    error_exit("not Vertex,PointRange or 2");
                if (cg_goto(cgfile, cgbase, "Zone_t", 1, srname, 0, NULL))
                    error_exit("cg_goto SubRegion[345]");
                if (cg_rind_read(rind)) error_exit("cg_rind_read");
            }
            else {
                if (location != CGNS_ENUMV(EdgeCenter) ||
                    ptype != CGNS_ENUMV(PointList) || np != npts)
                    error_exit("not EdgeCenter,PointList or npts");
            }
            if (cg_subreg_ptset_read(cgfile, cgbase, cgzone, ns, pts))
                error_exit("cg_subreg_ptset_read");
        }
    }
}

