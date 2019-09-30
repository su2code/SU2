/* This example illustrates using TecPolyFace and TecPolyBConn to create
 * two polygonal zones. The first zone contains six hexagons, and the
 * second zone contains a hexagon and an octagon. Refer to the Data
 * Format Guide for a picture of the configuration, including node
 * and face numbers.
 */

#include "TECIO.h"
#include "MASTER.h" /* for defintion of NULL */

int main()
{
    /* DOCSTART:hexagons_tecini.txt*/
    INTEGER4 I; /* use to check return values */

    INTEGER4 Debug      = 1;
    INTEGER4 VIsDouble  = 0;
    INTEGER4 FileType   = 0;
    INTEGER4 FileFormat = 0; // 0 == PLT, 1 == SZPLT

    I = TECINI142((char*)"Example: Multiple polygonal zones",
                  (char*)"X Y P",  /* Defines the variables for the data file.
                                    * Each zone must contain each of the vars
                                    * listed here. The order of the variables
                                    * in the list is used to define the
                                    * variable number (e.g. X is Variable 1).
                                    * When referring to variables in other
                                    * TecIO functions, you will refer to the
                                    * variable by its number.
                                    */
                  (char*)"multiplepolygons-HexAndOct.plt",
                  (char*)".",       /* scratch directory */
                  &FileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);
    /* DOCEND */

    /* DOCSTART:hexagons_zone1_teczne.txt*/
    /* TECZNE Parameters */
    INTEGER4 ZoneType           = 6;   /* FE Polygon */
    INTEGER4 NumPts_Z1          = 13;  /* the number of unique
                                        * nodes in the zone.
                                        */
    INTEGER4 NumElems_Z1        = 3;
    INTEGER4 NumFaces_Z1        = 15;  /* the number of unique
                                        * faces in the zone.
                                        */
    INTEGER4 ICellMax           = 0;   /* not used */
    INTEGER4 JCellMax           = 0;   /* not used */
    INTEGER4 KCellMax           = 0;   /* not used */
    double   SolutionTime       = 0.0;
    INTEGER4 StrandID           = 0;
    INTEGER4 ParentZone         = 0;
    INTEGER4 IsBlock            = 1;
    INTEGER4 NumFaceConnections = 0;
    INTEGER4 FaceNeighborMode   = 1;
    INTEGER4 SharConn           = 0;

    INTEGER4 ValueLocation[3]   = { 1, 1, 0 };

    /* For a polygonal zone, the total number of face nodes is
     * twice the total number of faces. This is because each face
     * is composed of exactly two nodes.
     */
    INTEGER4 TotalNumFaceNodes_Z1  = 2 * NumFaces_Z1;

    /* A boundary face is a face that is neighbored by an element
     * or elements in another zone or zone(s). In Zone 1, Face 9,
     * Face 10 and Face 12 have a neighbor in Zone 2. Therefore,
     * the total number of boundary faces is “3”.
     */
    INTEGER4 TotalNumBndryFaces_Z1 = 3;

    /* Each boundary face has one or more boundary connections. A
     * boundary connection is defined as another element in another
     * zone. Face 9 has a boundary connection with Element 1 in
     * Zone 2. In this example, each boundary face is connected to
     * one other element, so the total number of boundary
     * connections is equivalent to the total number of boundary
     * faces (3).
     */
    INTEGER4 TotalNumBndryConns_Z1 = 3;

    I = TECZNE142((char*)"Zone 1: 3 Hexagons", /* Specifies the name of
                                                * the entire dataset. When
                                                * the file is loaded into
                                                * Tecplot, the value is
                                                * available via the Data
                                                * Set Info dialog.
                                                */
                  &ZoneType,
                  &NumPts_Z1,
                  &NumElems_Z1,
                  &NumFaces_Z1,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &SolutionTime,
                  &StrandID,
                  &ParentZone,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  &TotalNumFaceNodes_Z1,
                  &TotalNumBndryFaces_Z1,
                  &TotalNumBndryConns_Z1,
                  NULL,
                  ValueLocation,
                  NULL,
                  &SharConn);
    /* DOCEND */

    /* DOCSTART:hexagons_zone1_tecdat.txt*/
    /* TECDAT Parameters */
    double *X_Z1 = new double[NumPts_Z1];
    double *Y_Z1 = new double[NumPts_Z1];

    X_Z1[0]  = 1;
    Y_Z1[0]  = 6;

    X_Z1[1]  = 2;
    Y_Z1[1]  = 6;

    X_Z1[2]  = 3;
    Y_Z1[2]  = 5;

    X_Z1[3]  = 2;
    Y_Z1[3]  = 4;

    X_Z1[4]  = 1;
    Y_Z1[4]  = 4;

    X_Z1[5]  = 0;
    Y_Z1[5]  = 5;

    X_Z1[6]  = 4;
    Y_Z1[6]  = 5;

    X_Z1[7]  = 5;
    Y_Z1[7]  = 4;

    X_Z1[8]  = 4;
    Y_Z1[8]  = 3;

    X_Z1[9]  = 3;
    Y_Z1[9]  = 3;

    X_Z1[10] = 2;
    Y_Z1[10] = 2;

    X_Z1[11] = 1;
    Y_Z1[11] = 2;

    X_Z1[12] = 0;
    Y_Z1[12] = 3;

    double *P_Z1 = new double[NumElems_Z1];
    P_Z1[0] = 2;
    P_Z1[1] = 4;
    P_Z1[2] = 5;

    INTEGER4  IsDouble = 1;
    I = TECDAT142(&NumPts_Z1,   X_Z1, &IsDouble);
    I = TECDAT142(&NumPts_Z1,   Y_Z1, &IsDouble);
    I = TECDAT142(&NumElems_Z1, P_Z1, &IsDouble);
    delete X_Z1;
    delete Y_Z1;
    delete P_Z1;
    /* DOCEND */

    /* DOCSTART:hexagons_zone1_facenodes.txt*/
    /* TecPolyFace Parameters */

    /* Create a FaceNodes array, dimensioned by the total number
     * of face nodes in the zone.
     */
    INTEGER4 *FaceNodes_Z1 = new INTEGER4[TotalNumFaceNodes_Z1];

    /* Face Nodes for Element 1 */
    FaceNodes_Z1[0]  = 1;
    FaceNodes_Z1[1]  = 2;

    FaceNodes_Z1[2]  = 2;
    FaceNodes_Z1[3]  = 3;

    FaceNodes_Z1[4]  = 3;
    FaceNodes_Z1[5]  = 4;

    FaceNodes_Z1[6]  = 4;
    FaceNodes_Z1[7]  = 5;

    FaceNodes_Z1[8]  = 5;
    FaceNodes_Z1[9]  = 6;

    FaceNodes_Z1[10] = 6;
    FaceNodes_Z1[11] = 1;

    /* Face Nodes for Element 2 */
    FaceNodes_Z1[12] = 3;
    FaceNodes_Z1[13] = 7;

    FaceNodes_Z1[14] = 7;
    FaceNodes_Z1[15] = 8;

    FaceNodes_Z1[16] = 8;
    FaceNodes_Z1[17] = 9;

    FaceNodes_Z1[18] = 9;
    FaceNodes_Z1[19] = 10;

    FaceNodes_Z1[20] = 10;
    FaceNodes_Z1[21] = 4;

    /* Face Nodes for Element 3 */
    FaceNodes_Z1[22] = 10;
    FaceNodes_Z1[23] = 11;

    FaceNodes_Z1[24] = 11;
    FaceNodes_Z1[25] = 12;

    FaceNodes_Z1[26] = 12;
    FaceNodes_Z1[27] = 13;

    FaceNodes_Z1[28] = 13;
    FaceNodes_Z1[29] = 5;
    /* DOCEND */

    /* Specify the right and left neighboring elements.
     * The neighboring elements can be determined using the
     * right-hand rule. For each face, place your right-hand along
     * the face with your fingers pointing the direction of
     * incrementing node numbers (i.e. from Node 1 to Node 2). The
     * right side of your hand will indicate the right element,
     * and the left side of your hand will indicate the left
     * element. A value of zero indicates that there is no
     * neighboring element on that side.  A negative value
     * indicates that the neighboring element is in another zone.
     * The number is a pointer into the FaceBndryConnectionElems
     * and FaceBndryConnectionZones arrays.
     */

    /* DOCSTART:hexagons_zone1_neighbors.txt*/
    INTEGER4 *FaceLeftElems_Z1  = new INTEGER4[NumFaces_Z1];
    INTEGER4 *FaceRightElems_Z1 = new INTEGER4[NumFaces_Z1];

    /* Left Face Elems for Element 1 */
    FaceLeftElems_Z1[0]  = 0;
    FaceLeftElems_Z1[1]  = 0;
    FaceLeftElems_Z1[2]  = 2;
    FaceLeftElems_Z1[3]  = 3;
    FaceLeftElems_Z1[4]  = 0;

    /* Left Face Elems for Element 2 */
    FaceLeftElems_Z1[5]  =  0;
    FaceLeftElems_Z1[6]  =  0;
    FaceLeftElems_Z1[7]  =  0;
    FaceLeftElems_Z1[8]  = -1;
    FaceLeftElems_Z1[9]  = -2;
    FaceLeftElems_Z1[10] =  3;

    /* Left Face Elems for Element 3 */
    FaceLeftElems_Z1[11]  = -3;
    FaceLeftElems_Z1[12]  =  0;
    FaceLeftElems_Z1[13]  =  0;
    FaceLeftElems_Z1[14]  =  0;

    /* Set Right Face Elems.  Because of the way we numbered the
     * nodes and faces, the right element for every face is the
     * element itself.
     */
    for (INTEGER4 ii = 0; ii < 6; ii++)
        FaceRightElems_Z1[ii] = 1;

    for (INTEGER4 ii = 6; ii < 11; ii++)
        FaceRightElems_Z1[ii] = 2;

    for (INTEGER4 ii = 11; ii <= 14; ii++)
        FaceRightElems_Z1[ii] = 3;

    I = TECPOLYFACE142(&NumFaces_Z1,
                       NULL,         /* Not used for polygon zones */
                       FaceNodes_Z1,
                       FaceLeftElems_Z1,
                       FaceRightElems_Z1);

    delete FaceNodes_Z1;
    delete FaceLeftElems_Z1;
    delete FaceRightElems_Z1;
    /* DOCEND */

    /* DOCSTART:hexagons_zone1_tecpoly.txt */
    /* TecPolyBConn Parameters */

    /* The FaceBndryConnectionCounts array is used to define the
     * number of boundary connections for each face that has a
     * boundary connection. For example, if a zone has three
     * boundary connections in total (NumConnectedBoundaryFaces),
     * two of those boundary connections are in one face, and the
     * remaining boundary connection is in a second face, the
     * FaceBndryConnectionCounts array would be: [2 1].
     *
     * In this example, the total number of connected boundary
     * faces (specified via TECZNE) is equal to three. Each
     * boundary face is connected to only one other element,
     * so the FaceBoundaryConnectionCounts array is (1, 1, 1).
     */
    INTEGER4 FaceBndryConnectionCounts_Z1[3]  = {1, 1, 1};

    /* The value(s) in the FaceBndryConnectionElems and
     * FaceBndryConnectionZones arrays specifies the element number
     * and zone number, respectively, that a given boundary
     * connection is connected to. In this case, the first
     * boundary connection face is connected to Element 1 in Zone 2
     * and the remaining connection is to Element 2 in Zone 2.
     */
    INTEGER4 FaceBndryConnectionElems_Z1[3]   = {1, 2, 2};
    INTEGER4 FaceBndryConnectionZones_Z1[3]   = {2, 2, 2};

    I = TECPOLYBCONN142(&TotalNumBndryFaces_Z1,
                        FaceBndryConnectionCounts_Z1,
                        FaceBndryConnectionElems_Z1,
                        FaceBndryConnectionZones_Z1);
    /* DOCEND */

    /* Define Zone 2.  Zone 2 contains an octagon and a hexagon. */
    /* TECZNE Parameters */
    /* DOCSTART:hexagons_zone2_teczne.txt*/
    INTEGER4 NumPts_Z2               = 12; /* number of unique
                                            * nodes in the zone
                                            */
    INTEGER4 NumElems_Z2             = 2;
    INTEGER4 NumFaces_Z2             = 13; /* number of unique
                                            * faces in the zone
                                            */
    INTEGER4 NumFaceConnections_Z2   = 0;
    /* In polygonal zones, each face has exactly two nodes */
    INTEGER4 TotalNumFaceNodes_Z2    = NumFaces_Z2 * 2;

    /* A boundary face is a face that is neighbored by an element or
     * elements from another zone or zone(s). In Zone 2, Face 6,
     * Face 7 and Face 13 have a neighbor in Zone 1. Therefore, the
     * total number of boundary faces is “3”.
     */
    INTEGER4 TotalNumBndryFaces_Z2   = 3;

    /* Each boundary face has one or more boundary connections. In
     * this example, each boundary face is connected to one other
     * element (i.e. the number of boundary faces and the number of
     * boundary connections is one-to-one).
     */
    INTEGER4 TotalNumBndryConns_Z2   = 3;

    I = TECZNE142((char*)"Zone 2: 1 Hexagon and 1 Octagon",
                  &ZoneType,
                  &NumPts_Z2,
                  &NumElems_Z2,
                  &NumFaces_Z2,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &SolutionTime,
                  &StrandID,
                  &ParentZone,
                  &IsBlock,
                  &NumFaceConnections_Z2,
                  &FaceNeighborMode,
                  &TotalNumFaceNodes_Z2,
                  &TotalNumBndryFaces_Z2,
                  &TotalNumBndryConns_Z2,
                  NULL,
                  ValueLocation,
                  NULL,
                  &SharConn);
    /* DOCEND */

    /* TECDAT Parameters */
    /* DOCSTART:hexagons_zone2_tecdat.txt*/
    double *X_Z2 = new double[NumPts_Z2];
    double *Y_Z2 = new double[NumPts_Z2];

    X_Z2[0]      = 5;
    Y_Z2[0]      = 4;

    X_Z2[1]      = 6;
    Y_Z2[1]      = 4;

    X_Z2[2]      = 7;
    Y_Z2[2]      = 3;

    X_Z2[3]      = 6;
    Y_Z2[3]      = 2;

    X_Z2[4]      = 5;
    Y_Z2[4]      = 2;

    X_Z2[5]      = 4;
    Y_Z2[5]      = 3;

    X_Z2[6]      = 3;
    Y_Z2[6]      = 3;

    X_Z2[7]      = 5;
    Y_Z2[7]      = 1;

    X_Z2[8]      = 4;
    Y_Z2[8]      = 0;

    X_Z2[9]      = 3;
    Y_Z2[9]      = 0;

    X_Z2[10]     = 2;
    Y_Z2[10]     = 1;

    X_Z2[11]     = 2;
    Y_Z2[11]     = 2;

    /* In the call to TecZne, P was set to a cell centered variable.
     * As such, only two values need to be defined.
     */
    double  *P_Z2 = new double[NumPts_Z2];

    P_Z2[0] = 8;
    P_Z2[1] = 6;

    I = TECDAT142(&NumPts_Z2,   X_Z2, &IsDouble);
    I = TECDAT142(&NumPts_Z2,   Y_Z2, &IsDouble);
    I = TECDAT142(&NumElems_Z2, P_Z2, &IsDouble);

    delete X_Z2;
    delete Y_Z2;
    delete P_Z2;
    /* DOCEND */

    /* TecPolyFace Parameters */
    /* DOCSTART:hexagons_zone2_facemap.txt*/
    INTEGER4 *FaceNodes_Z2;
    FaceNodes_Z2 = new INTEGER4[TotalNumFaceNodes_Z2];

    /* Face Nodes for Element 1 */
    FaceNodes_Z2[0]  = 1;
    FaceNodes_Z2[1]  = 2;

    FaceNodes_Z2[2]  = 2;
    FaceNodes_Z2[3]  = 3;

    FaceNodes_Z2[4]  = 3;
    FaceNodes_Z2[5]  = 4;

    FaceNodes_Z2[6]  = 4;
    FaceNodes_Z2[7]  = 5;

    FaceNodes_Z2[8]  = 5;
    FaceNodes_Z2[9]  = 6;

    FaceNodes_Z2[10] = 6;
    FaceNodes_Z2[11] = 1;

    /* Face Nodes for Element 2 */
    FaceNodes_Z2[12] = 7;
    FaceNodes_Z2[13] = 6;

    FaceNodes_Z2[14] = 5;
    FaceNodes_Z2[15] = 8;

    FaceNodes_Z2[16] = 8;
    FaceNodes_Z2[17] = 9;

    FaceNodes_Z2[18] = 9;
    FaceNodes_Z2[19] = 10;

    FaceNodes_Z2[20] = 10;
    FaceNodes_Z2[21] = 11;

    FaceNodes_Z2[22] = 11;
    FaceNodes_Z2[23] = 12;

    FaceNodes_Z2[24] = 12;
    FaceNodes_Z2[25] = 7;
    /* DOCEND */

    /* DOCSTART:hexagons_zone2_tecpolyface.txt*/
    /* Specify the right and left neighboring elements.
     * The neighboring elements can be determined using the
     * right-hand rule. For each face, place your right-hand along
     * the face with your fingers pointing the direction of
     * incrementing node numbers (i.e. from Node 1 to Node 2). The
     * right side of your hand will indicate the right element,
     * and the left side of your hand will indicate the left
     * element. A value of zero indicates that there is no
     * neighboring element on that side.  A negative value
     * indicates that the neighboring element is in another zone.
     * The number is a pointer into the FaceBndryConnectionElems
     * and FaceBndryConnectionZones arrays.
     */

    INTEGER4 *FaceLeftElems_Z2  = new INTEGER4[NumFaces_Z2];
    INTEGER4 *FaceRightElems_Z2 = new INTEGER4[NumFaces_Z2];

    /* Left Face Elems for Element 1 */
    FaceLeftElems_Z2[0]  =  0;
    FaceLeftElems_Z2[1]  =  0;
    FaceLeftElems_Z2[2]  =  0;
    FaceLeftElems_Z2[3]  =  0;
    FaceLeftElems_Z2[4]  =  2;
    FaceLeftElems_Z2[5]  = -1;

    /* Left Face Elems for Element 2 */
    FaceLeftElems_Z2[6]   = -2;
    FaceLeftElems_Z2[7]   =  0;
    FaceLeftElems_Z2[8]   =  0;
    FaceLeftElems_Z2[9]   =  0;
    FaceLeftElems_Z2[10]  =  0;
    FaceLeftElems_Z2[11]  =  0;
    FaceLeftElems_Z2[12]  = -3;

    /* Set Right Face Elems.  Because of the way we numbered the
     * nodes and faces, the right element for every face is the
     * element itself. */
    for (INTEGER4 ii = 0; ii < 6; ii++)
        FaceRightElems_Z2[ii]  = 1;

    for (INTEGER4 ii = 6; ii < 13; ii++)
        FaceRightElems_Z2[ii]  = 2;

    I = TECPOLYFACE142(&NumFaces_Z2,
                       NULL,
                       FaceNodes_Z2,
                       FaceLeftElems_Z2,
                       FaceRightElems_Z2);

    delete FaceNodes_Z2;
    delete FaceLeftElems_Z2;
    delete FaceRightElems_Z2;
    /* DOCEND */

    /* DOCSTART:hexagons_zone2_tecpoly.txt*/
    /* The FaceBndryConnectionCounts array is used to define the
     * number of boundary connections for each face that has a
     * boundary connection. In this example, the total number of
     * connected boundary faces (specified via TECZNE) is equal to
     * three. Each boundary face is connected to only one other
     * element, so the FaceBoundaryConnectionCounts array is
     * (1, 1, 1).
     */
    INTEGER4 FaceBndryConnectionCounts_Z2[3]  = {1, 1, 1};

    /* The value(s) in the FaceBndryConnectionElems and
     * FaceBndryConnectionZones arrays specifies that element
     * number and zone number, respectively, that a given boundary
     * connection is connected to. In this case, the first boundary
     * connection face is connected to Element 2 in Zone 1 and the
     * remaining connections are Element 3 in Zone 1.
     */
    INTEGER4 FaceBndryConnectionElems_Z2[3]   = {2, 3, 3};
    INTEGER4 FaceBndryConnectionZones_Z2[3]   = {1, 1, 1};

    I = TECPOLYBCONN142(&TotalNumBndryFaces_Z2,
                        FaceBndryConnectionCounts_Z2,
                        FaceBndryConnectionElems_Z2,
                        FaceBndryConnectionZones_Z2);
    /* DOCEND */

    /* DOCSTART:hexagons_tecend.txt*/
    I = TECEND142();
    /* DOCEND */

    return 0;
}
