/* This example illustrates using separate grid
 * and solution files.
 */

#include "TECIO.h"
#include "MASTER.h" /* for defintion of NULL */

int main()
{
    /* DOCSTART:gridsolution_grid_tecini.txt*/
    INTEGER4 I; /* use to check return values */

    INTEGER4 Debug      = 1;
    INTEGER4 VIsDouble  = 0;
    INTEGER4 FileType   = 1; /* 1 = grid file. */
    INTEGER4 FileFormat = 0; // 0 == PLT, 1 == SZPLT; only PLT is currently
                             // supported for grid/solution files

    I = TECINI142((char*)"Example: Separate grid and solution files",
                  (char*)"X Y Z",  /* Defines the variables for the data file.
                                    * Each zone must contain each of the vars
                                    * listed here. The order of the variables
                                    * in the list is used to define the
                                    * variable number (e.g. X is Variable 1).
                                    * When referring to variables in other
                                    * TecIO functions, you will refer to the
                                    * variable by its number.
                                    */
                  (char*)"gridsolution-grid.plt",
                  (char*)".",       /* scratch directory */
                  &FileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);
    /* DOCEND */

    /* DOCSTART:gridsolution_grid_teczne.txt*/
    /* TECZNE Parameters */
    INTEGER4 ZoneType           = 7;   /* FE Polyhedron */
    INTEGER4 NumPts             = 20;  /* the number of unique
                                        * nodes in the zone.
                                        */
    INTEGER4 NumElems           = 1;
    INTEGER4 NumFaces           = 12;  /* the number of unique
                                        * faces in the zone.
                                        */
    INTEGER4 ICellMax           = 0;   /* not used */
    INTEGER4 JCellMax           = 0;   /* not used */
    INTEGER4 KCellMax           = 0;   /* not used */
    double   SolutionTime       = 0.0;
    INTEGER4 StrandID           = 1;   /* time strand for
                                        * unsteady solution.
                                        */
    INTEGER4 ParentZone         = 0;
    INTEGER4 IsBlock            = 1;
    INTEGER4 NumFaceConnections = 0;
    INTEGER4 FaceNeighborMode   = 1;
    INTEGER4 SharConn           = 0;

    /* For this zone, the total number of face nodes is
     * five times number of faces, because each face
     * is a pentagon.
     */
    INTEGER4 TotalNumFaceNodes  = 5 * NumFaces;

    /* This zone has no connected boundary faces.
     */
    INTEGER4 TotalNumBndryFaces = 0;
    INTEGER4 TotalNumBndryConns = 0;

    I = TECZNE142((char*)"Dodecahedron", /* Name of the zone. */
                  &ZoneType,
                  &NumPts,
                  &NumElems,
                  &NumFaces,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &SolutionTime,
                  &StrandID,
                  &ParentZone,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  &TotalNumFaceNodes,
                  &TotalNumBndryFaces,
                  &TotalNumBndryConns,
                  NULL,
                  NULL, /* All nodal variables */
                  NULL,
                  &SharConn);
    /* DOCEND */

    /* DOCSTART:gridsolution_grid_tecdat.txt*/
    /* TECDAT Parameters */
    double  Phi = 0.5 * (1.0 + sqrt(5.0));
    double  Pi  = 3.141592653578;
    double *X  = new double[NumPts];
    double *Y  = new double[NumPts];
    double *Z  = new double[NumPts];
    int     Count = 0;

    for(int J = 0; J <= 4; J++)
    {
        X[Count] = 2.0 * cos(2.0 / 5.0 * Pi * J);
        Y[Count] = 2.0 * sin(2.0 / 5.0 * Pi * J);
        Z[Count] = Phi + 1.0;
        Count++;

        X[Count] = -X[Count - 1];
        Y[Count] = -Y[Count - 1];
        Z[Count] = -Z[Count - 1];
        Count++;

        X[Count] = 2.0 * Phi * cos(2.0 / 5.0 * Pi * J);
        Y[Count] = 2.0 * Phi * sin(2.0 / 5.0 * Pi * J);
        Z[Count] = Phi - 1.0;
        Count++;

        X[Count] = -X[Count - 1];
        Y[Count] = -Y[Count - 1];
        Z[Count] = -Z[Count - 1];
        Count++;
    }

    INTEGER4  IsDouble = 1;
    
    I = TECDAT142(&NumPts, X, &IsDouble);
    I = TECDAT142(&NumPts, Y, &IsDouble);
    I = TECDAT142(&NumPts, Z, &IsDouble);
    
    delete X;
    delete Y;
    delete Z;
    /* DOCEND */

    /* DOCSTART:gridsolution_grid_facenodes.txt*/
    /* TecPolyFace Parameters */

    /* Create a FaceNodes array, dimensioned by the total number
     * of face nodes in the zone.
     */
    INTEGER4 *FaceNodes = new INTEGER4[TotalNumFaceNodes];
    int n = 0;

    /* Face Nodes for face 1 of the dodecahedron */
    FaceNodes[n++] = 2;
    FaceNodes[n++] = 6;
    FaceNodes[n++] = 10;
    FaceNodes[n++] = 14;
    FaceNodes[n++] = 18;

    /* Face Nodes for face 2 */
    FaceNodes[n++] = 6;
    FaceNodes[n++] = 8;
    FaceNodes[n++] = 19;
    FaceNodes[n++] = 12;
    FaceNodes[n++] = 10;

    /* Face Nodes for face 3 */
    FaceNodes[n++] = 3;
    FaceNodes[n++] = 12;
    FaceNodes[n++] = 10;
    FaceNodes[n++] = 14;
    FaceNodes[n++] = 16;

    /* Face Nodes for face 4 */
    FaceNodes[n++] = 7;
    FaceNodes[n++] = 16;
    FaceNodes[n++] = 14;
    FaceNodes[n++] = 18;
    FaceNodes[n++] = 20;

    /* Face Nodes for face 5 */
    FaceNodes[n++] = 2;
    FaceNodes[n++] = 4;
    FaceNodes[n++] = 11;
    FaceNodes[n++] = 20;
    FaceNodes[n++] = 18;

    /* Face Nodes for face 6 */
    FaceNodes[n++] = 2;
    FaceNodes[n++] = 4;
    FaceNodes[n++] = 15;
    FaceNodes[n++] = 8;
    FaceNodes[n++] = 6;

    /* Face Nodes for face 7 */
    FaceNodes[n++] = 1;
    FaceNodes[n++] = 3;
    FaceNodes[n++] = 12;
    FaceNodes[n++] = 19;
    FaceNodes[n++] = 17;

    /* Face Nodes for face 8 */
    FaceNodes[n++] = 1;
    FaceNodes[n++] = 3;
    FaceNodes[n++] = 16;
    FaceNodes[n++] = 7;
    FaceNodes[n++] = 5;

    /* Face Nodes for face 9 */
    FaceNodes[n++] = 5;
    FaceNodes[n++] = 7;
    FaceNodes[n++] = 20;
    FaceNodes[n++] = 11;
    FaceNodes[n++] = 9;

    /* Face Nodes for face 10 */
    FaceNodes[n++] = 4;
    FaceNodes[n++] = 11;
    FaceNodes[n++] = 9;
    FaceNodes[n++] = 13;
    FaceNodes[n++] = 15;

    /* Face Nodes for face 11 */
    FaceNodes[n++] = 8;
    FaceNodes[n++] = 15;
    FaceNodes[n++] = 13;
    FaceNodes[n++] = 17;
    FaceNodes[n++] = 19;

    /* Face Nodes for face 12 */
    FaceNodes[n++] = 1;
    FaceNodes[n++] = 5;
    FaceNodes[n++] = 9;
    FaceNodes[n++] = 13;
    FaceNodes[n++] = 17;
    /* DOCEND */

    /* Specify the number of nodes for each face, and the right and
     * left neighboring elements. The neighboring elements can be
     * determined using the right-hand rule. For each face, curl
     * the fingers of your right hand in the direction of
     * incrementing node numbers (i.e. from Node 1 to Node 2 and
     * so on). Your thumb will point toward the right element.
     * A value of zero indicates that there is no
     * neighboring element on that side.  A negative value
     * indicates that the neighboring element is in another zone.
     * In that case, the number is a pointer into the
     * FaceBndryConnectionElems and FaceBndryConnectionZones arrays.
     */

    /* DOCSTART:gridsolution_grid_tecpoly.txt*/
    INTEGER4 *FaceNodeCounts = new INTEGER4[NumFaces];
    INTEGER4 *FaceLeftElems  = new INTEGER4[NumFaces];
    INTEGER4 *FaceRightElems = new INTEGER4[NumFaces];

    /* For this particular zone, each face has the 5 nodes. */
    for(int J = 0; J < NumFaces; J++)
        FaceNodeCounts[J] = 5;

    /* Set the right and left elements for each face. */
    FaceRightElems[0]  = 1;
    FaceRightElems[1]  = 1;
    FaceRightElems[2]  = 0;
    FaceRightElems[3]  = 0;
    FaceRightElems[4]  = 0;
    FaceRightElems[5]  = 1;
    FaceRightElems[6]  = 1;
    FaceRightElems[7]  = 0;
    FaceRightElems[8]  = 0;
    FaceRightElems[9]  = 1;
    FaceRightElems[10] = 1;
    FaceRightElems[11] = 0;

    FaceLeftElems[0]  = 0;
    FaceLeftElems[1]  = 0;
    FaceLeftElems[2]  = 1;
    FaceLeftElems[3]  = 1;
    FaceLeftElems[4]  = 1;
    FaceLeftElems[5]  = 0;
    FaceLeftElems[6]  = 0;
    FaceLeftElems[7]  = 1;
    FaceLeftElems[8]  = 1;
    FaceLeftElems[9]  = 0;
    FaceLeftElems[10] = 0;
    FaceLeftElems[11] = 1;

    I = TECPOLYFACE142(&NumFaces,
                       FaceNodeCounts,
                       FaceNodes,
                       FaceLeftElems,
                       FaceRightElems);

    delete FaceNodes;
    delete FaceLeftElems;
    delete FaceRightElems;
    /* DOCEND */

    /* DOCSTART:gridsolution_grid_tecend.txt*/
    I = TECEND142();
    /* DOCEND */

    for(int J = 0; J < 5; J++)
    {
        char SolutionFileName[128];
        sprintf(SolutionFileName, "gridsolution-solution%d.plt", J);

        /* DOCSTART:gridsolution_solution_tecini.txt*/
        FileType  = 2; /* 2 = solution file. */

        I = TECINI142((char*)"Example: Separate grid and solution files",
                      (char*)"P T",  /* Defines the variables for the solution file.
                                      * Note that these are different variables from
                                      * the grid file.
                                      */
                      SolutionFileName,
                      (char*)".",       /* scratch directory */
                      &FileFormat,
                      &FileType,
                      &Debug,
                      &VIsDouble);
        /* DOCEND */

        /* DOCSTART:gridsolution_solution_teczne.txt*/
        /* TECZNE Parameters are mostly unchanged from creation of the grid file. */
        TotalNumFaceNodes = 0;
        SolutionTime = J;

        char ZoneName[128];
        sprintf(ZoneName, "Dodecahedron Time=%g", SolutionTime);
        I = TECZNE142(ZoneName,
                      &ZoneType,
                      &NumPts,
                      &NumElems,
                      &NumFaces,
                      &ICellMax,
                      &JCellMax,
                      &KCellMax,
                      &SolutionTime,
                      &StrandID,
                      &ParentZone,
                      &IsBlock,
                      &NumFaceConnections,
                      &FaceNeighborMode,
                      &TotalNumFaceNodes,
                      &TotalNumBndryFaces,
                      &TotalNumBndryConns,
                      NULL,
                      NULL, /* All nodal variables */
                      NULL,
                      &SharConn);
        /* DOCEND */

        /* DOCSTART:gridsolution_solution_tecdat.txt*/
        /* TECDAT Parameters */
        double *P = new double[NumPts];
        double *T = new double[NumPts];

        for(int K = 0; K < NumPts; K++)
        {
            P[K] = (double)(K + J);
            T[K] = 1.0 + K + K;
        }

        I = TECDAT142(&NumPts, P, &IsDouble);
        I = TECDAT142(&NumPts, T, &IsDouble);

        delete P;
        delete T;
        /* DOCEND */

        /* DOCSTART:gridsolution_solution_tecend.txt*/
        I = TECEND142();
        /* DOCEND */
    }
        
    return 0;
}
