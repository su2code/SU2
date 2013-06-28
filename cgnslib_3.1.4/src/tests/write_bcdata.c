#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
# include <io.h>
# define unlink _unlink
#else
# include <unistd.h>
#endif
#include "utils.h"

char *fname = "bcdata.cgns";

int cgfile, cgbase, cgzone, cgcoord, cgsect, cgbc, cgfam, cgdset;
char name[33];

cgsize_t size[] = {4, 1, 0};
float coord[3][4] = {
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}
};
cgsize_t tris[12] = {
    1, 3, 2,
    1, 2, 4,
    1, 4, 3,
    2, 3, 4
};
cgsize_t tets[4] = {1, 2, 3, 4};

cgsize_t dim, npnts = 4;
cgsize_t pnts[4] = {1, 2, 3, 4};
float data[4] = {1, 2, 3, 4};

float start, finish;

int main (int argc, char **argv)
{
    int i, j, nb = 5, nv = 5;

    if (argc > 1) {
        nb = atoi (argv[1]);
        if (argc > 2)
            nv = atoi (argv[2]);
    }

    unlink (fname);
    printf ("creating file ...");
    fflush (stdout);
    start = (float)elapsed_time();

    if (cg_open (fname, CG_MODE_WRITE, &cgfile) ||
        cg_base_write (cgfile, "Base", 3, 3, &cgbase) ||
        cg_zone_write (cgfile, cgbase, "Zone", size,
            CGNS_ENUMV(Unstructured), &cgzone) ||
        cg_coord_write (cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
            "CoordinateX", coord[0], &cgcoord) ||
        cg_coord_write (cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
            "CoordinateY", coord[1], &cgcoord) ||
        cg_coord_write (cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
            "CoordinateZ", coord[2], &cgcoord) ||
        cg_section_write (cgfile, cgbase, cgzone, "Tris",
            CGNS_ENUMV(TRI_3), 1, 4, 0, tris, &cgsect) ||
        cg_section_write (cgfile, cgbase, cgzone, "Tets",
            CGNS_ENUMV(TETRA_4), 5, 5, 0, tets, &cgsect))
        cg_error_exit();

    dim = npnts;
    for (j = 1; j <= nb; j++) {
        sprintf (name, "BC%d", j);
        if (cg_boco_write (cgfile, cgbase, cgzone, name, CGNS_ENUMV(BCWall),
                CGNS_ENUMV(ElementList), npnts, pnts, &cgbc))
            cg_error_exit();

        for (i = 1; i <= nv; i++) {
            sprintf (name, "BCData%d", i);
            if (cg_dataset_write (cgfile, cgbase, cgzone, cgbc, name,
				  CGNS_ENUMV(BCWall), &cgdset) ||
                cg_bcdata_write (cgfile, cgbase, cgzone, cgbc, cgdset,
				 CGNS_ENUMV(Dirichlet)) ||
                cg_goto (cgfile, cgbase, "Zone_t", 1, "ZoneBC_t", 1,
                    "BC_t", cgbc, "BCDataSet_t", cgdset, "end") ||
                cg_descriptor_write ("label", name) ||
                cg_descriptor_write ("basis", "intensive") ||
                cg_goto (cgfile, cgbase, "Zone_t", 1, "ZoneBC_t", 1,
                    "BC_t", cgbc, "BCDataSet_t", cgdset,
		    "BCData_t", CGNS_ENUMV(Dirichlet), "end") ||
                cg_array_write ("Data", CGNS_ENUMV(RealSingle), 1, &dim, data))
                cg_error_exit();
        }
    }

    dim = 1;
    for (j = 1; j <= nb; j++) {
        sprintf (name, "Family%d", j);
        if (cg_family_write (cgfile, cgbase, name, &cgfam) ||
            cg_fambc_write (cgfile, cgbase, cgfam, "FamilyBC",
                CGNS_ENUMV(BCWall), &cgbc))
            cg_error_exit();

        for (i = 1; i <= nv; i++) {
            sprintf (name, "FamilyBCData%d", i);
            if (cg_goto (cgfile, cgbase, "Family_t", j, "FamilyBC_t", 1, NULL) ||
                cg_bcdataset_write (name, CGNS_ENUMV(BCWall),
                    CGNS_ENUMV(Dirichlet)) ||
                cg_gorel (cgfile, "FamilyBCDataSet_t", i, NULL) ||
                cg_descriptor_write ("label", name) ||
                cg_descriptor_write ("basis", "intensive") ||
                cg_gorel (cgfile, "DirichletData", 0, NULL) ||
                cg_array_write ("Data", CGNS_ENUMV(RealSingle), 1, &dim, data))
                cg_error_exit();
        }
    }

    finish = (float)elapsed_time();
    printf (" %.2f secs\n", finish - start);

    printf ("closing file  ...");
    fflush (stdout);
    start = (float)elapsed_time();
    if (cg_close(cgfile)) cg_error_exit();
    finish = (float)elapsed_time();
    printf (" %.2f secs\n", finish - start);

    printf ("opening file  ...");
    fflush (stdout);
    start = (float)elapsed_time();
    if (cg_open (fname, CG_MODE_MODIFY, &cgfile)) cg_error_exit();
    finish = (float)elapsed_time();
    printf (" %.2f secs\n", finish - start);
    cg_close (cgfile);

    return 0;
}

