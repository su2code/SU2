/* This example illustrates how to create two simple
 * FE-quadilateral zones and create a face neighbor
 * connection between the two zones.  In order to keep the
 * example as simple as possible, error checking is not included.
 */

#include "TECIO.h"
#include "MASTER.h"

int main()
{
    /* Initialize the Data File using TECINI.  TECINI is required
     * for all data files.  It is used to: open the data file and
     * initialize the file header information (name the data file,
     * the variables for the data file, and the file type).
     */

    /* DOCSTART:faceneighbors_tecini.txt*/

    INTEGER4 Debug      = 1;
    INTEGER4 VIsDouble  = 0;
    INTEGER4 FileType   = 0;
    INTEGER4 FileFormat = 0; // 0 == PLT, 1 == SZPLT
    INTEGER4 I          = 0;  /* Used to track return codes */

    I = TECINI142((char*)"Face Neighbors Example", /* Specifies the name
                                                    * of the entire
                                                    * dataset
                                                    */
                  (char*)"X Y P",                  /* Defines the
                                                    * variables for the
                                                    * data file.  Each
                                                    * zone must contain
                                                    * each of the vars
                                                    * listed.  The order
                                                    * of the variables in
                                                    * the list is used to
                                                    * define the variable
                                                    * number (e.g. X is
                                                    * Var 1.)
                                                    */

                  (char*)"FaceNeighbors.plt",      /* Specifies the
                                                    * file name.
                                                    */
                  (char*)".",
                  &FileFormat,
                  &FileType,                       /* The FileType is set to
                                                    * zero, indicating it is
                                                    * a full file (containing
                                                    * both grid and solution
                                                    * data).
                                                    */
                  &Debug,
                  &VIsDouble);
    /* DOCEND */

    /* After TECINI is called, call TECZNE to create one or
     * more zones for your data file.
     */
    /* DOCSTART:faceneighbors_teczne1.txt*/
    INTEGER4 ZoneType   = 3;     /* set the zone type to
                                * FEQuadrilateral
                                */
    INTEGER4 NumPts     = 6;
    INTEGER4 NumElems   = 2;
    INTEGER4 NumFaces   = 8;
    INTEGER4 ICellMax   = 0;     /* not used */
    INTEGER4 JCellMax   = 0;     /* not used */
    INTEGER4 KCellMax   = 0;     /* not used */
    double   SolTime    = 360.0;
    INTEGER4 StrandID   = 0;     /* StaticZone */
    INTEGER4 ParentZn   = 0;
    INTEGER4 IsBlock    = 1;     /* Block */
    INTEGER4 NFConns    = 1;     /* Specify the number of Face
                                * Neighbor Connections in the
                                * Zone.  When this value is
                                * greater than zero, TECFACE must
                                * be called prior to creating the
                                * next zone or ending the file.
                                */

    /* Specify the Face Neighbor Mode.
     * A value of 2 indicated that the face neighbor mode is global
     * one-to-one.  The scope of the face neighbors (local or
     * global) is with respect to the zones. A value of global
     * indicates that the face neighbor(s) is/are shared aross zones;
     * a value of local indicates that the face neighbor(s) are
     * shared within the current zone.  The terms one-to-one and
     * one-to-many are used to indicate whether the face in question
     * is shared with one cell or several cells.
     * For example, if your data is arranged as follows:

                   -----------------------
                  |       |       |       |
                  |   1   |   2   |   3   |
                  |       |       |       |
                   -----------------------
                  |       |               |
                  |   4   |      5        |
                  |       |               |
                   -----------------------
     * The face between 1 & 4 is local-one-to-one. The face between
     * 5 and (2 & 3) is local one-to-many.
     */

    INTEGER4 FNMode                       = 2;

    INTEGER4 TotalNumFaceNodes            = 1;  /* Not used for
                                               * FEQuad zones*/
    INTEGER4 NumConnectedBoundaryFaces    = 1;  /* Not used for
                                               * FEQuad zones*/
    INTEGER4 TotalNumBoundaryConnections  = 1;  /* Not used for
                                               * FEQuad zones*/
    INTEGER4 ShrConn                      = 0;

    INTEGER4 ValueLocation[3] = {1, 1, 1};  /* Specify the variable
                                         * values at the nodes.
                                         * NOTE:  Because all of
                                         * the variables are
                                         * defined at the nodes,
                                         * we can just pass
                                         * NULL for this array.
                                         * We are providing the
                                         * array for illustration
                                         * purposes.
                                         */

    I = TECZNE142((char*)"Zone 1",
                  &ZoneType,
                  &NumPts,
                  &NumElems,
                  &NumFaces,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &SolTime,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NFConns,
                  &FNMode,
                  &TotalNumFaceNodes,
                  &NumConnectedBoundaryFaces,
                  &TotalNumBoundaryConnections,
                  NULL,
                  ValueLocation,
                  NULL,
                  &ShrConn);
    /* DOCEND */

    /* Set up the variable values.  The variable values will be
     * written to the file using TECDAT.  Because we are specifying
     * nodal variables (as specified via the ValueLocation
     * parameter in TECZNE, each variable is dimensioned by the
     * number of points (NumPts) in the Zone.  You have the option
     * to specify some variables with nodal values and some with
     * cell-centered values.  Refer to the Binary Function
     * Reference for details.
     */

    /* DOCSTART:faceneighbors_tecdat1.txt*/
    float *X = new float[NumPts];
    float *Y = new float[NumPts];
    float *P = new float[NumPts];

    /* For this example, we will create 2 rectangular cells in Zone
     * 1.  Before defining your variables, you must establish a
     * consistent node numbering scheme for your data.  Once the
     * node numbers are defined, supply the variable values in the
     * node numbering order.  In this example, node 1 is defined at
     * X = 0 and Y = 0.  As such, the first value supplied for X
     * (i.e. X[0]) is 0.  Similarly, the first value supplied for Y
     * is 0.
     *
     * It is important that you refer to node numbers consistently.
     * The node numbers will be used later to define the
     * connectivity for each element.
     */

    X[0] = 0;
    X[1] = 0;
    X[2] = 1;
    X[3] = 1;
    X[4] = 2;
    X[5] = 2;

    Y[0] = 0;
    Y[1] = 1;
    Y[2] = 0;
    Y[3] = 1;
    Y[4] = 0;
    Y[5] = 1;

    for (INTEGER4 ii = 0; ii < NumPts; ii++)
        P[ii] = (float)(NumPts - ii);

    INTEGER4 DIsDouble =  0;  /* Set DIsDouble to zero to use
                             * variables in float format.
                             */

    /* Call TECDAT once for each variable */
    I   = TECDAT142(&NumPts, &X[0], &DIsDouble);
    I   = TECDAT142(&NumPts, &Y[0], &DIsDouble);
    I   = TECDAT142(&NumPts, &P[0], &DIsDouble);
    /* DOCEND */

    /*  Define the face neighbors connections.
     *  The Connectivity List is used to specify the nodes that
     *  compose each element.  When working with nodal variables, the
     *  numbering of the nodes is implicitly defined when the
     *  variables are declared.  The first value of each variable is
     *  for node one, the second value for node two, and so on.
     *
     *  Because this zone contains two quadilateral elements, we must
     *  supply 8 values in the connectivity list.  The first four
     *  values define the nodes that form element 1. Similarly, the
     *  second four values define the nodes that form element 2.
     */

    /* DOCSTART:faceneighbors_tecnod1.txt*/
    INTEGER4 ConnList[8] = {1, 3, 4, 2,
                            3, 5, 6, 4
                           };
    I   = TECNOD142(ConnList);
    /* DOCEND */

    /*  TIP! It is important to provide the node list in either a
     *  clockwise or counter-clockwise order.  Otherwise, your
     *  elements will be misformed.  For example, if the first two
     *  numbers in the above connectivity list, the zone would
     *  appear as follows:
     */

    /*  Now that TECNOD has been called, the creation of Zone 1
     *  is complete.  However, in this example, we will define a
     *  face neighbor between Zone 1 and Zone 2 (to be created
     *  later in the example).  Face Neighbor connections are used
     *  to define connections that are not created via the
     *  connectivity list.  For example, local face neighbors may
     *  need to be defined when a zone wraps itself and global face
     *  neighbors may need to be defined to smooth edges across
     *  zones.  Face Neighbors are used when deriving variables and
     *  drawing contours.
     *
     *  In this example, we are creating a face neighbor connection
     *  between cell 2 in Zone 1 and cell 1 in Zone 2.  The
     *  information required when specifying face  neighbors
     *  depends upon the type of connection.
     *
     *  In this case, we must supply (in this order):
     *    - the cell number in the current zone that contains the
     *    - the number of the face in that cell that contains the
     *      face neighbor
     *    - the number of the other zone to which the face is
     *      connected
     *    - the number of the cell in the other zone to which the
     *      face is connected
     *  The face numbering for cell-based finite elements is
     *  defined using the picture displayed in the Data Format
     *  Guide.  In this example, face 2 in cell 2 in the current
     *  zone is connected to cell 1 in zone 2.
     */

    /* DOCSTART:faceneighbors_tecface1.txt*/
    INTEGER4 FaceConn[4] = {2, 2, 2, 1};
    I   = TECFACE142(FaceConn);
    /* DOCEND */

    /* The creation of Zone 1 is complete.  We are ready to create
     * Zone 2.  For simplicity, Zone 2 is a copy of Zone 1 shifted
     * along the X-axis.  As such, many of the variables used to
     * create Zone 1 are re-used here.
     */
    /* DOCSTART:faceneighbors_teczne2.txt*/
    /* Call TECZNE to create Zone 2 */
    I = TECZNE142((char*)"Zone 2",
                  &ZoneType,
                  &NumPts,
                  &NumElems,
                  &NumFaces,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &SolTime,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NFConns,
                  &FNMode,
                  &TotalNumFaceNodes,
                  &NumConnectedBoundaryFaces,
                  &TotalNumBoundaryConnections,
                  NULL,
                  ValueLocation,
                  NULL,
                  &ShrConn);
    /* DOCEND */

    /* Define the variables for Zone 2.  Because Zone 2 is a copy
     * of Zone 1, shifted along the X-axis, we can share the Y
     * variable definition used to Zone.  We will also create a
     * second pressure variable for Zone 2 (P2).
     */

    /* DOCSTART:faceneighbors_tecdat2.txt*/
    float *X2 = new float[NumPts];
    float *P2 = new float[NumPts];

    for (INTEGER4 ii = 0; ii < NumPts; ii++)
    {
        X2[ii] = X[ii] + 2;
        P2[ii] = 2 * (float)ii;
    }

    I   = TECDAT142(&NumPts, &X2[0], &DIsDouble);
    I   = TECDAT142(&NumPts, &Y[0], &DIsDouble);
    I   = TECDAT142(&NumPts, &P2[0], &DIsDouble);

    delete X;
    delete Y;
    delete P;
    delete X2;
    delete P2;
    /* DOCEND */

    /* As with Zone 1, we must define the connectivity list for
     * Zone 2. Because, the node numbering restarts at one for each
     * new zone and the nodal arrangement is identical between the
     * two zones, we may reuse the connectivity list from Zone 1.
     */

    /* DOCSTART:faceneighbors_tecnod2.txt*/
    I   = TECNOD142(ConnList);
    /* DOCEND */

    /* We will now specify the face neighbor connection with
     * respect to our new current zone of Zone 2.
     */

    /* DOCSTART:faceneighbors_tecface2.txt*/
    INTEGER4 FaceConn2[4] = {1, 4, 1, 2};  /* cell 1, face 4 in
                                       * current zone is a
                                       * neighbor to cell 2 in
                                       * zone 1.
                                       */
    I   = TECFACE142(FaceConn2);
    /* DOCEND */

    /* Call TECEND to close the file */
    /* DOCSTART:faceneighbors_tecend.txt*/
    I = TECEND142();
    /* DOCEND */

    return 0;
}
