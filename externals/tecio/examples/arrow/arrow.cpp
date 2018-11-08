/* This example creates two simple polyhedral zones in the shape
 * of a three-dimensional arrow.  Obscured boundary faces are used.
 */

#include <stdio.h>
#include "TECIO.h"

int main()
{
    /* DOCSTART:arrow_tecini.txt*/
    INTEGER4 Debug     = 1;
    INTEGER4 VIsDouble = 1;
    INTEGER4 FileFormat = 0; // 0 == PLT, 1 == SZPLT; Only PLT is currently
                             // supported for ployhedral zones
    INTEGER4 FileType  = 0;
    INTEGER4 I;

    /* Open the file and write the Tecplot datafile
     * header information
     */
    I = TECINI142((char*)"Multiple polyhedral zones", /* Name of the entire
                                                       * dataset.
                                                       */
                  (char*)"X Y Z P",  /* Defines the variables for the data
                                      * file. Each zone must contain each of
                                      * the variables listed here. The order
                                      * of the variables in the list is used
                                      * to define the variable number (e.g.
                                      * X is Var 1).
                                      */
                  (char*)"Arrow.plt",
                  (char*)".",         /* Scratch Directory */
                  &FileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);
    /* DOCEND */

    /* After TECINI is called, call TECZNE to create one or more
     * zones for your data file. In this example, Zone 1 contains a
     * single rectangular solid created as a face-based finite-element
     * (i.e. polyhedral zone). The zone has eight points (or nodes),
     * six faces and one element.
     */
    /* DOCSTART:arrow_teczne_rect.txt*/
    /* TECZNE Parameters */
    INTEGER4 ZoneType                    = 7; /* sets the zone type
                                             * to polyhedral */
    INTEGER4 NumPts_Rect                 = 8;
    INTEGER4 NumElems_Rect               = 1;
    INTEGER4 NumFaces_Rect               = 6;
    INTEGER4 ICellMax                    = 0; /* not used */
    INTEGER4 JCellMax                    = 0; /* not used */
    INTEGER4 KCellMax                    = 0; /* not used */
    double   SolutionTime                = 0.0;
    INTEGER4 StrandID                    = 0;
    INTEGER4 ParentZone                  = 0;
    INTEGER4 IsBlock                     = 1;
    INTEGER4 NumFaceConnections          = 0; /* not used */
    INTEGER4 FaceNeighborMode            = 1; /* not used */
    INTEGER4 SharConn                    = 0;

    /* In a rectangular solid, each face is composed of four nodes.
     * As such, the total number of face nodes is twenty-four (four
     * nodes for each of the six faces).
     */
    INTEGER4 TotalNumFaceNodes_Rect      = 24;

    /* There is one connected boundary face in this zone (the face on
     * the rectangle adjacent to the arrowhead). Refer to the Data
     * Format Guide for additional information. */
    INTEGER4 NumConnBndryFaces_Rect      = 1;

    /* The connected boundary face has one connection, the face on
     * the bottom of the arrowhead. A connection is an element-zone
     * tuple that indicates a neighboring element (and its zone) when
     * the neighboring element is in a different zone. Generally,
     * there will be one boundary connection for each boundary face.
     */
    INTEGER4 TotalNumBndryConns_Rect     = 1;

    /* For illustrative purposes, the grid variables (X, Y, and Z)
     * are nodal variables (i.e. ValueLocation = 1), and the pressure
     * variable (P) is a cell-centered variable (i.e.
     * ValueLocation = 0).
     */
    INTEGER4 ValueLocation[4] = { 1, 1, 1, 0 };

    I = TECZNE142((char*)"Zone 1: Rectangular Solid",
                  &ZoneType,
                  &NumPts_Rect,
                  &NumElems_Rect,
                  &NumFaces_Rect,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &SolutionTime,
                  &StrandID,
                  &ParentZone,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  &TotalNumFaceNodes_Rect,
                  &NumConnBndryFaces_Rect,
                  &TotalNumBndryConns_Rect,
                  NULL,
                  ValueLocation,
                  NULL,
                  &SharConn);
    /* DOCEND */

    /* DOCSTART:arrow_tecdat_rect.txt*/
    //set variable values (X_Rect, Y_Rect, Z_Rect & P_Rect)
    double *X_Rect = new double[NumPts_Rect];
    double *Y_Rect = new double[NumPts_Rect];
    double *Z_Rect = new double[NumPts_Rect];
    double *P_Rect = new double[NumElems_Rect];

    for (INTEGER4 ii = 0; ii <= NumPts_Rect / 2; ii += 4)
    {
        X_Rect[ii]   = 0;
        X_Rect[ii+1] = 3;
        X_Rect[ii+2] = 3;
        X_Rect[ii+3] = 0;

        Y_Rect[ii]   = 3;
        Y_Rect[ii+1] = 3;
        Y_Rect[ii+2] = 1;
        Y_Rect[ii+3] = 1;
    }

    for (INTEGER4 ii = 0; ii < 4; ii++)
        Z_Rect[ii] = 0;

    for (INTEGER4 ii = 4; ii < NumPts_Rect; ii++)
        Z_Rect[ii] = -2;

    P_Rect[0] = 10;

    INTEGER4 IsDouble = 1;
    I = TECDAT142(&NumPts_Rect,   X_Rect, &IsDouble);
    I = TECDAT142(&NumPts_Rect,   Y_Rect, &IsDouble);
    I = TECDAT142(&NumPts_Rect,   Z_Rect, &IsDouble);
    I = TECDAT142(&NumElems_Rect, P_Rect, &IsDouble);
    /* DOCEND */

    /* DOCSTART:arrow_facenodes_rect.txt*/
    /* The FaceNodeCounts array is used to describe the number of
     * nodes in each face of the zone. The first value in the array
     * is the number of nodes in Face 1, the second value is the
     * number of nodes in Face 2 and so forth. In this example, each
     * face of the zone has four nodes.
     */

    INTEGER4 *FaceNodeCounts_Rect = new INTEGER4[NumFaces_Rect];
    //For this particular zone, each face has the 4 nodes
    for (INTEGER4 ii = 0; ii < NumFaces_Rect; ii++)
        FaceNodeCounts_Rect[ii] = 4;

    /* The FaceNodes array is used to specify the nodes that compose
     * each face.  For each face (n of N), the number of nodes used
     * to define the face is specified by the nth value in the
     * FaceNodeCounts array. For example, if the first value in the
     * FaceNodeCounts array is 4 (indicating Face 1 is composed of
     * four nodes), the first four values in the FaceNodes array are
     * the node numbers of the nodes in Face 1.
     *
     * ------------
     * WARNING
     * When providing the node numbers for each face, you must
     * provide the node numbers in a consistent order (either
     * clockwise or counter-clockwise. Providing the node numbers
     * out of order results in contorted faces.
     * ------------
     */

    INTEGER4 *FaceNodes_Rect  = new INTEGER4[TotalNumFaceNodes_Rect];

    //Nodes for Face 1
    FaceNodes_Rect[0]   = 1;
    FaceNodes_Rect[1]   = 2;
    FaceNodes_Rect[2]   = 3;
    FaceNodes_Rect[3]   = 4;

    //Nodes for Face 2
    FaceNodes_Rect[4]   = 1;
    FaceNodes_Rect[5]   = 4;
    FaceNodes_Rect[6]   = 8;
    FaceNodes_Rect[7]   = 5;

    //Nodes for Face 3
    FaceNodes_Rect[8]   = 5;
    FaceNodes_Rect[9]   = 8;
    FaceNodes_Rect[10]  = 7;
    FaceNodes_Rect[11]  = 6;

    //Nodes for Face 4
    FaceNodes_Rect[12]  = 2;
    FaceNodes_Rect[13]  = 6;
    FaceNodes_Rect[14]  = 7;
    FaceNodes_Rect[15]  = 3;

    //Nodes for Face 5
    FaceNodes_Rect[16]  = 6;
    FaceNodes_Rect[17]  = 2;
    FaceNodes_Rect[18]  = 1;
    FaceNodes_Rect[19]  = 5;

    //Nodes for Face 6
    FaceNodes_Rect[20]  = 3;
    FaceNodes_Rect[21]  = 7;
    FaceNodes_Rect[22]  = 8;
    FaceNodes_Rect[23]  = 4;
    /* DOCEND */

    /* DOCSTART:arrow_neighbors_rect.txt*/
    INTEGER4 *FaceLeftElems_Rect  = new INTEGER4[NumFaces_Rect];
    INTEGER4 *FaceRightElems_Rect = new INTEGER4[NumFaces_Rect];

    /* Since this zone has just one element, all leftelems are
     * NoNeighboring Element and all right elems are itself
     */
    for (INTEGER4 ii = 0; ii < NumFaces_Rect; ii++)
    {
        FaceRightElems_Rect[ii] = 1;
        FaceLeftElems_Rect[ii]  = 0;
    }

    /* The negative value in the FaceLeftElems array indicates that
     * the face is connected to an element in another zone. In this
     * case, Face 4 is connected to a face in Zone 2 (to be defined
     * later in the example). The FaceBoundaryConnectionElems array
     * lists all of the element numbers in other zones that the
     * current zone shares boundary connections with. Similarly, the
     * FaceBoundaryConnectionZones array lists all of the zone numbers
     * with which the current zone shares boundaries.  A negative
     * value in the FaceLeftElems or FaceRightElems array indicates
     * the position within these arrays that defines the neighboring
     * element and zone for a face.
     *
     * For example, if the FaceBoundaryConnectionElems array is:
     * [1 8 2] and the FaceBoundaryConnectionZones array is: [2 5 3],
     * a FaceLeftElems or FaceRightElems value of -2 indicates that
     * the face in question has a boundary connection with Element 8
     * in Zone 5.
     */
    FaceLeftElems_Rect[3]    = -1;

    I = TECPOLYFACE142(&NumFaces_Rect,
                       FaceNodeCounts_Rect,
                       FaceNodes_Rect,
                       FaceLeftElems_Rect,
                       FaceRightElems_Rect);
    /* DOCEND */

    /* DOCSTART:arrow_tecpoly_rect.txt*/
    /* The FaceBndryConnectionCounts array is used to define the
     * number of boundary connections for each face that has a
     * boundary connection. For example, if a zone has three boundary
     * connections in total (NumConnectedBoundaryFaces), two of those
     * boundary connections are in one face, and the remaining
     * boundary connection is in a second face, the
     * FaceBndryConnectionCounts array would be: [2 1].
     * In this example, the total number of connected boundary faces
     * (specified via TECZNE) is equal to one, so the
     * FaceBoundaryConnectionCounts array contains a single value (1).
     */
    INTEGER4 *FaceBndryConnCounts_Rect = new INTEGER4[NumConnBndryFaces_Rect];
    FaceBndryConnCounts_Rect[0] = 1;

    /* The value(s) in the FaceBndryConnectionElems and
     * FaceBndryConnectionZones arrays specify the element number and
     * zone number, respectively, that a given boundary connection is
     * connected to. In this case, the boundary connection face is
     * connected to Element 1 in Zone 2.
     */
    INTEGER4 *FaceBndryConnElems_Rect = new INTEGER4[TotalNumBndryConns_Rect];
    INTEGER4 *FaceBndryConnZones_Rect = new INTEGER4[TotalNumBndryConns_Rect];

    FaceBndryConnElems_Rect[0]  = 1;
    FaceBndryConnZones_Rect[0]  = 2;

    I = TECPOLYBCONN142(&NumConnBndryFaces_Rect,
                        FaceBndryConnCounts_Rect,
                        FaceBndryConnElems_Rect,
                        FaceBndryConnZones_Rect);

    /* cleanup */
    delete X_Rect;
    delete Y_Rect;
    delete Z_Rect;
    delete P_Rect;
    delete FaceNodeCounts_Rect;
    delete FaceNodes_Rect;
    delete FaceLeftElems_Rect;
    delete FaceRightElems_Rect;
    delete FaceBndryConnCounts_Rect;
    delete FaceBndryConnElems_Rect;
    delete FaceBndryConnZones_Rect;
    /* DOCEND */

    /* The data for Zone 1 has been written to the data file, so we
     * are ready to create Zone 2. For simplicity, we will reuse many
     * of the variables created for the rectangular zone that are not
     * relevant to this tutorial. */

    /* Zone 2 (the arrowhead or prism) has a single element composed
     * of six nodes and five faces.
     */

    /* DOCSTART:arrow_teczne_prism.txt*/
    //TECZNE Parameters
    INTEGER4 NumPts_Prism             = 6;
    INTEGER4 NumElems_Prism           = 1;
    INTEGER4 NumFaces_Prism           = 5;

    /* The prism is composed of two triangular faces and three
     * rectangular faces.  The total number of face nodes is the sum
     * of the nodes in each triangular face (2 times 3) and the nodes
     * in each rectangular face (3 times 4).
     */
    INTEGER4 TotalNumFaceNodes_Prism  = 18;

    /* As with Zone 1, Zone 2 has one connected boundary face, the
     * face that is connected to Zone 1.
     */
    INTEGER4 NumConnBndryFaces_Prism  = 1;

    /* In this case, we have set the total number of boundary
     * connections for the connected face to two. The first boundary
     * connection is the connection to Zone 1.  The second boundary
     * connection is used to indicate that the face is only partially
     * obscured by the face from Zone 1. If we omitted the second
     * boundary connection, the connected face of the prism would
     * disappear if the rectangular zone was deactivated.
     */
    INTEGER4 TotalNumBndryConns_Prism = 2;

    I = TECZNE142((char*)"Zone 2: Prism",
                  &ZoneType,
                  &NumPts_Prism,
                  &NumElems_Prism,
                  &NumFaces_Prism,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &SolutionTime,
                  &StrandID,
                  &ParentZone,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  &TotalNumFaceNodes_Prism,
                  &NumConnBndryFaces_Prism,
                  &TotalNumBndryConns_Prism,
                  NULL,
                  ValueLocation,
                  NULL,
                  &SharConn);
    /* DOCEND */

    /* DOCSTART:arrow_tecdat_prism.txt*/
    double *X_Prism = new double[NumPts_Prism];
    double *Y_Prism = new double[NumPts_Prism];
    double *Z_Prism = new double[NumPts_Prism];

    /* Set the X and Y variable values, one z-plane at a time */
    double ZVal = 0;
    for (INTEGER4 ii = 0; ii < 2; ii++)
    {
        // triangle in Z=ZVal plane
        X_Prism[3*ii]   = 3;
        Y_Prism[3*ii]   = 4;
        Z_Prism[3*ii]   = ZVal;

        X_Prism[3*ii+1] = 7;
        Y_Prism[3*ii+1] = 2;
        Z_Prism[3*ii+1] = ZVal;

        X_Prism[3*ii+2] = 3;
        Y_Prism[3*ii+2] = 0;
        Z_Prism[3*ii+2] = ZVal;

        ZVal = ZVal - 2;
    }

    /* When we called TecZne, we specified that the variable 4
     * (pressure) is cell-centered. As such, only NumElements number
     * of values needs to be written to the data file for the pressure
     * variable.
     */
    double *P_Prism = new double[NumElems_Prism];
    P_Prism[0] = 20;

    I = TECDAT142(&NumPts_Prism,  X_Prism, &IsDouble);
    I = TECDAT142(&NumPts_Prism,  Y_Prism, &IsDouble);
    I = TECDAT142(&NumPts_Prism,  Z_Prism, &IsDouble);
    I = TECDAT142(&NumElems_Prism, P_Prism, &IsDouble);
    /* DOCEND */

    /* DOCSTART:arrow_facemap_prism.txt*/
    INTEGER4 *FaceNodeCounts_Prism = new INTEGER4[NumFaces_Prism];
    INTEGER4 *FaceNodes_Prism    = new INTEGER4[TotalNumFaceNodes_Prism];

    /* Because of the way we chose to number our faces, the first
     * three faces are rectangular and the last two are triangular.
     * The numbering of the faces is arbitrary, but the faces must
     * be referred to consistently.
     */
    for (INTEGER4 ii = 0; ii < 3; ii++)
        FaceNodeCounts_Prism[ii] = 4;

    for (INTEGER4 ii = 3; ii < NumFaces_Prism; ii++)
        FaceNodeCounts_Prism[ii] = 3;

    //Nodes for Face 1
    FaceNodes_Prism[0]   = 1;
    FaceNodes_Prism[1]   = 3;
    FaceNodes_Prism[2]   = 6;
    FaceNodes_Prism[3]   = 4;

    //Nodes for Face 2
    FaceNodes_Prism[4]   = 1;
    FaceNodes_Prism[5]   = 4;
    FaceNodes_Prism[6]   = 5;
    FaceNodes_Prism[7]   = 2;

    //Nodes for Face 3
    FaceNodes_Prism[8]   = 3;
    FaceNodes_Prism[9]   = 2;
    FaceNodes_Prism[10]  = 5;
    FaceNodes_Prism[11]  = 6;

    //Nodes for Face 4
    FaceNodes_Prism[12]  = 5;
    FaceNodes_Prism[13]  = 4;
    FaceNodes_Prism[14]  = 6;

    //Nodes for Face 5
    FaceNodes_Prism[15]  = 1;
    FaceNodes_Prism[16]  = 2;
    FaceNodes_Prism[17]  = 3;
    /* DOCEND */

    /* DOCSTART:arrow_neighbors_prism.txt*/
    /* Since this zone has just one element, all leftelems are
     * NoNeighboring Element and all right elems are itself.
     */
    INTEGER4 *FaceLeftElems_Prism  =  new INTEGER4[NumFaces_Prism];
    INTEGER4 *FaceRightElems_Prism =  new INTEGER4[NumFaces_Prism];

    for (INTEGER4 ii = 0; ii < NumFaces_Prism; ii++)
    {
        FaceRightElems_Prism[ii] = 1;
        FaceLeftElems_Prism[ii]  = 0;
    }

    /* The negative value in the FaceLeftElems array indicates that
     * the face is connected to an element in another zone. In this
     * case, Face 1 is connected to a face in Zone 1 (as indicated in
     * Line 6). The FaceBoundaryConnectionElems array lists all of
     * the element numbers in other zones that the current zone shares
     * boundary connections with. Similarly, the
     * FaceBoundaryConnectionZones array lists all of the zone numbers
     * with which the current zone shares boundaries.  A negative
     * value in the FaceLeftElems or FaceRightElems array indicates
     * the position within these arrays that defines the neighboring
     * element and zone for a face.
     */
    FaceLeftElems_Prism[0]  = -1;

    I = TECPOLYFACE142(&NumFaces_Prism,
                       FaceNodeCounts_Prism,
                       FaceNodes_Prism,
                       FaceLeftElems_Prism,
                       FaceRightElems_Prism);
    /* DOCEND */

    /* DOCSTART:arrow_tecpoly_prism.txt*/
    INTEGER4 *FaceBndryConnCounts_Prism = new INTEGER4[NumConnBndryFaces_Prism];
    FaceBndryConnCounts_Prism[0] = 2;

    INTEGER4 *FaceBndryConnElems_Prism = new INTEGER4[TotalNumBndryConns_Prism];
    INTEGER4 *FaceBndryConnZones_Prism = new INTEGER4[TotalNumBndryConns_Prism];

    /* As previously mentioned, a connected boundary face is a face
     * that has either multiple neighboring faces or neighbor(s) that
     * belong to another zone.  Those cases are sufficient when the
     * combination of all of the face’s neighbors completely cover the
     * face. However, there are some cases (such as the bottom of the
     * arrowhead) where the face is not completely covered by its
     * neighbors. In those cases the face is referred to as “partially
     * obscured”. A partially obscured face is indicated by
     * incrementing the value in TotalNumConnectedBoundaryFaces and
     * entering a value of 0 in both the FaceBndryConnectionElems and
     * FaceBoundaryConnectionZones arrays for the boundary connection
     * for the partially obscured face.
     */
    FaceBndryConnElems_Prism[0]  = 0;
    FaceBndryConnZones_Prism[0]  = 0;

    /* Indicates that Face 1 is connected to Element 1 in Zone 1. */
    FaceBndryConnElems_Prism[1]  = 1;
    FaceBndryConnZones_Prism[1]  = 1;

    I = TECPOLYBCONN142(&NumConnBndryFaces_Prism,
                        FaceBndryConnCounts_Prism,
                        FaceBndryConnElems_Prism,
                        FaceBndryConnZones_Prism);

    /* cleanup */
    delete X_Prism;
    delete Y_Prism;
    delete Z_Prism;
    delete P_Prism;
    delete FaceNodeCounts_Prism;
    delete FaceNodes_Prism;
    delete FaceLeftElems_Prism;
    delete FaceRightElems_Prism;
    delete FaceBndryConnCounts_Prism;
    delete FaceBndryConnElems_Prism;
    delete FaceBndryConnZones_Prism;
    /* DOCEND */

    /* DOCSTART:arrow_tecend.txt*/
    I = TECEND142();
    /* DOCEND */

    return 0;
}
