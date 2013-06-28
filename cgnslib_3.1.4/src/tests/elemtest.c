#include <stdio.h>
#include <stdlib.h>
#if defined(_WIN32) && !defined(__NUTC__)
# include <io.h>
# define unlink _unlink
#else
# include <unistd.h>
#endif
#include "cgnslib.h"

static double xc[42] = {
    0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 1.0, 1.0, 1.0,
    1.0, 2.0, 1.0, 0.0, 0.0, 2.0, 2.0, 0.0, 1.0, 2.0, 1.0, 0.0,
    1.0, 1.0, 2.0, 1.0, 0.0, 1.0, 1.0,
    1.5, 0.5, 1.0, 1.5, 0.5, 1.5, 0.5,
    0.5, 1.5, 1.5, 0.5,
    1.0
};
static double yc[42] = {
    0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 4.0, 4.0, 2.0,
    0.0, 1.0, 2.0, 1.0, 0.0, 0.0, 2.0, 2.0, 0.0, 1.0, 2.0, 1.0,
    1.0, 0.0, 1.0, 2.0, 1.0, 1.0, 1.0,
    3.0, 3.0, 4.0, 3.0, 3.0, 3.0, 3.0,
    1.0, 1.0, 2.0, 2.0,
    3.0
};
static double zc[42] = {
    0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 0.0, 2.0, 4.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0,
    0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0,
    0.0, 0.0, 1.0, 2.0, 2.0, 1.0, 1.0,
    3.0, 3.0, 3.0, 3.0,
    3.0
};

static cgsize_t node[1]   = {1};
static cgsize_t bar[3]    = {1, 2, 12};
static cgsize_t tri[6]    = {4, 3, 9, 14, 31, 32};
static cgsize_t quad[9]   = {1, 2, 3, 4, 12, 13, 14, 15, 24};
static cgsize_t tetra[10] = {8, 7, 10, 11, 22, 35, 36, 41, 40, 42};
static cgsize_t pyra[14]  = {5, 6, 7, 8, 11, 20, 21, 22, 23, 38, 39,
                             40, 41, 29};
static cgsize_t penta[18] = {4, 3, 9, 8, 7, 10, 14, 31, 32, 19, 18, 35,
                             22, 34, 35, 27, 36, 37};
static cgsize_t hexa[27]  = {1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15,
                             16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                             26, 27, 28, 29, 30};

static int npoly = 16;
static cgsize_t poly[] = {4, 1, 4, 3, 2,
                          4, 2, 3, 7, 6,
                          4, 3, 4, 8, 7,
                          4, 1, 5, 8, 4,
                          4, 1, 5, 6, 2,
                          4, 5, 6, 7, 8,
                          3, 3, 4, 9,
                          3, 8, 7, 10,
                          4, 3, 9, 10, 7,
                          4, 4, 8, 10, 9,
                          3, 5, 6, 11,
                          3, 6, 7, 11,
                          3, 7, 8, 11,
                          3, 5, 11, 8,
                          3, 7, 10, 11,
                          3, 8, 11, 10};
static int nface = 4;
static cgsize_t face[] = {6, 1, 2, 3, 4, 5, 6,
                          5, -3, 7, 8, 9, 10,
                          5, -6, 11, 12, 13, 14,
                          4, -8, -13, 15, 16};

static cgsize_t elems[256];

int main (int argc, char *argv[])
{
    int j, n;
    cgsize_t ne;
    int fnum, bnum, znum, snum, cnum;
    cgsize_t size[3];
    char *outfile = "elemtest.cgns";

    unlink (outfile);

    if (cg_open (outfile, CG_MODE_WRITE, &fnum) ||
        cg_base_write (fnum, "Base", 3, 3, &bnum))
        cg_error_exit ();

    /* zone with linear elements */

    size[0] = 11;
    size[1] = 12;
    size[2] = 0;
    if (cg_zone_write (fnum, bnum, "Linear", size,
            CGNS_ENUMV(Unstructured), &znum) ||
        cg_coord_write (fnum, bnum, znum, CGNS_ENUMV(RealDouble),
            "CoordinateX", xc, &cnum) ||
        cg_coord_write (fnum, bnum, znum, CGNS_ENUMV(RealDouble),
            "CoordinateY", yc, &cnum) ||
        cg_coord_write (fnum, bnum, znum, CGNS_ENUMV(RealDouble),
            "CoordinateZ", zc, &cnum))
        cg_error_exit ();
    ne = j = 0;

    /* NGON_n first so polyhedra face references are correct */
 
    if (cg_section_write (fnum, bnum, znum, "NGON_n", CGNS_ENUMV(NGON_n),
            ne+1, ne+npoly, 0, poly, &snum))
        cg_error_exit();
    ne += npoly;

    /* NODE */

    if (cg_section_write (fnum, bnum, znum, "NODE", CGNS_ENUMV(NODE),
            ne+1, ne+1, 0, node, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(NODE);
    elems[j++] = node[0];

    /* BAR_2 */

    if (cg_section_write (fnum, bnum, znum, "BAR_2", CGNS_ENUMV(BAR_2),
            ne+1, ne+1, 0, bar, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(BAR_2);
    for (n = 0; n < 2; n++)
        elems[j++] = bar[n];
    
    /* TRI_3 */

    if (cg_section_write (fnum, bnum, znum, "TRI_3", CGNS_ENUMV(TRI_3),
            ne+1, ne+1, 0, tri, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(TRI_3);
    for (n = 0; n < 3; n++)
        elems[j++] = tri[n];
    
    /* QUAD_4 */

    if (cg_section_write (fnum, bnum, znum, "QUAD_4", CGNS_ENUMV(QUAD_4),
            ne+1, ne+1, 0, quad, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(QUAD_4);
    for (n = 0; n < 4; n++)
        elems[j++] = quad[n];
    
    /* TETRA_4 */
    
    if (cg_section_write (fnum, bnum, znum, "TETRA_4", CGNS_ENUMV(TETRA_4),
            ne+1, ne+1, 0, tetra, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(TETRA_4);
    for (n = 0; n < 4; n++)
        elems[j++] = tetra[n];
    
    /* PYRA_5 */

    if (cg_section_write (fnum, bnum, znum, "PYRA_5", CGNS_ENUMV(PYRA_5),
            ne+1, ne+1, 0, pyra, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(PYRA_5);
    for (n = 0; n < 5; n++)
        elems[j++] = pyra[n];
 
    /* PENTA_6 */

    if (cg_section_write (fnum, bnum, znum, "PENTA_6", CGNS_ENUMV(PENTA_6),
            ne+1, ne+1, 0, penta, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(PENTA_6);
    for (n = 0; n < 6; n++)
        elems[j++] = penta[n];
   
    /* HEXA_8 */

    if (cg_section_write (fnum, bnum, znum, "HEXA_8", CGNS_ENUMV(HEXA_8),
            ne+1, ne+1, 0, hexa, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(HEXA_8);
    for (n = 0; n < 8; n++)
        elems[j++] = hexa[n];

    /* MIXED */
 
    if (cg_section_write (fnum, bnum, znum, "MIXED", CGNS_ENUMV(MIXED),
            ne+1, ne+8, 0, elems, &snum))
        cg_error_exit ();
    ne += 8;

    /* NFACE_n */

    if (cg_section_write (fnum, bnum, znum, "NFACE_n", CGNS_ENUMV(NFACE_n),
            ne+1, ne+nface, 0, face, &snum))
        cg_error_exit ();

    /* zone with quadratic elements */

    size[0] = 42;
    size[1] = 8;
    size[2] = 0;
    if (cg_zone_write (fnum, bnum, "Quadratic", size,
            CGNS_ENUMV(Unstructured), &znum) ||
        cg_coord_write (fnum, bnum, znum, CGNS_ENUMV(RealDouble),
            "CoordinateX", xc, &cnum) ||
        cg_coord_write (fnum, bnum, znum, CGNS_ENUMV(RealDouble),
            "CoordinateY", yc, &cnum) ||
        cg_coord_write (fnum, bnum, znum, CGNS_ENUMV(RealDouble),
            "CoordinateZ", zc, &cnum))
        cg_error_exit ();
    ne = j = 0;

    /* BAR_3 */

    if (cg_section_write (fnum, bnum, znum, "BAR_3", CGNS_ENUMV(BAR_3),
            ne+1, ne+1, 0, bar, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(BAR_3);
    for (n = 0; n < 3; n++)
        elems[j++] = bar[n];
    
    /* TRI_6 */

    if (cg_section_write (fnum, bnum, znum, "TRI_6", CGNS_ENUMV(TRI_6),
            ne+1, ne+1, 0, tri, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(TRI_6);
    for (n = 0; n < 6; n++)
        elems[j++] = tri[n];
    
    /* QUAD_8 */

    if (cg_section_write (fnum, bnum, znum, "QUAD_8", CGNS_ENUMV(QUAD_8),
            ne+1, ne+1, 0, quad, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(QUAD_8);
    for (n = 0; n < 8; n++)
        elems[j++] = quad[n];
    
    /* TETRA_10 */
    
    if (cg_section_write (fnum, bnum, znum, "TETRA_10", CGNS_ENUMV(TETRA_10),
            ne+1, ne+1, 0, tetra, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(TETRA_10);
    for (n = 0; n < 10; n++)
        elems[j++] = tetra[n];
    
    /* PYRA_13 */

    if (cg_section_write (fnum, bnum, znum, "PYRA_13", CGNS_ENUMV(PYRA_13),
            ne+1, ne+1, 0, pyra, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(PYRA_13);
    for (n = 0; n < 13; n++)
        elems[j++] = pyra[n];
 
    /* PENTA_15 */

    if (cg_section_write (fnum, bnum, znum, "PENTA_15", CGNS_ENUMV(PENTA_15),
            ne+1, ne+1, 0, penta, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(PENTA_15);
    for (n = 0; n < 15; n++)
        elems[j++] = penta[n];
   
    /* HEXA_20 */

    if (cg_section_write (fnum, bnum, znum, "HEXA_20", CGNS_ENUMV(HEXA_20),
            ne+1, ne+1, 0, hexa, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(HEXA_20);
    for (n = 0; n < 20; n++)
        elems[j++] = hexa[n];

    /* MIXED */
 
    if (cg_section_write (fnum, bnum, znum, "MIXED", CGNS_ENUMV(MIXED),
            ne+1, ne+7, 0, elems, &snum))
        cg_error_exit ();
    ne += 7;

    /* zone with quadratic elements with mid-nodes */

    size[0] = 42;
    size[1] = 8;
    size[2] = 0;
    if (cg_zone_write (fnum, bnum, "Quadratic with mid-nodes", size,
            CGNS_ENUMV(Unstructured), &znum) ||
        cg_coord_write (fnum, bnum, znum, CGNS_ENUMV(RealDouble),
            "CoordinateX", xc, &cnum) ||
        cg_coord_write (fnum, bnum, znum, CGNS_ENUMV(RealDouble),
            "CoordinateY", yc, &cnum) ||
        cg_coord_write (fnum, bnum, znum, CGNS_ENUMV(RealDouble),
            "CoordinateZ", zc, &cnum))
        cg_error_exit ();
    ne = j = 0;

    /* QUAD_9 */

    if (cg_section_write (fnum, bnum, znum, "QUAD_9", CGNS_ENUMV(QUAD_9),
            ne+1, ne+1, 0, quad, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(QUAD_9);
    for (n = 0; n < 9; n++)
        elems[j++] = quad[n];
    
    /* TETRA_10 */
    
    if (cg_section_write (fnum, bnum, znum, "TETRA_10", CGNS_ENUMV(TETRA_10),
            ne+1, ne+1, 0, tetra, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(TETRA_10);
    for (n = 0; n < 10; n++)
        elems[j++] = tetra[n];
    
    /* PYRA_14 */

    if (cg_section_write (fnum, bnum, znum, "PYRA_14", CGNS_ENUMV(PYRA_14),
            ne+1, ne+1, 0, pyra, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(PYRA_14);
    for (n = 0; n < 14; n++)
        elems[j++] = pyra[n];
 
    /* PENTA_18 */

    if (cg_section_write (fnum, bnum, znum, "PENTA_18", CGNS_ENUMV(PENTA_18),
            ne+1, ne+1, 0, penta, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(PENTA_18);
    for (n = 0; n < 18; n++)
        elems[j++] = penta[n];
   
    /* HEXA_27 */

    if (cg_section_write (fnum, bnum, znum, "HEXA_27", CGNS_ENUMV(HEXA_27),
            ne+1, ne+1, 0, hexa, &snum))
        cg_error_exit ();
    ne++;

    elems[j++] = (int)CGNS_ENUMV(HEXA_27);
    for (n = 0; n < 27; n++)
        elems[j++] = hexa[n];

    /* MIXED */
 
    if (cg_section_write (fnum, bnum, znum, "MIXED", CGNS_ENUMV(MIXED),
            ne+1, ne+5, 0, elems, &snum))
        cg_error_exit ();
    ne += 5;
 
    if (cg_close (fnum)) cg_error_exit ();
    return 0;
}
