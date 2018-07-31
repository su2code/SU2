/* This example illustrates how to create a zone with a single
 * polygonal cell.  Please refer to the Data Format Guide for
 * additional information, including figures that display node
 * numbering.
 */

#include "TECIO.h"
#include "MASTER.h" /* for defintion of NULL */

int main()
{
    /* DOCSTART:octagon_tecini.txt*/
    INTEGER4 Debug      = 1;
    INTEGER4 VIsDouble  = 0;
    INTEGER4 FileType   = 0;
    INTEGER4 FileFormat = 0; // 0 == PLT, 1 == SZPLT; Only PLT is currently
                             // supported for polygonal zones.
    INTEGER4 I;             /* used to check return codes */

    /*
     * Open the file and write the Tecplot datafile
     * header information
     */

    I = TECINI142((char*)"Octagon",
                  (char*)"X Y P",   /* Defines the variables for the data
                                     * file. Each zone must contain each
                                     * of the vars listed here. The order
                                     * of the variables in the list is
                                     * used to define the variable number
                                     * (e.g. X is Variable 1). When
                                     * referring to variables in other
                                     * TecIO functions, you will refer to
                                     * thevariable by its number.
                                     */
                  (char*)"octagon.plt",
                  (char*)".",       /* scratch directory */
                  &FileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);
    /* DOCEND */

    /* Declare TECZNE variables */
    /* DOCSTART:octagon_teczne.txt*/
    /* In this example, we will create a single octagonal cell in
     * Tecplot 360's polyhedral file format.
     */
    INTEGER4 ZoneType  = 6;     /* FEPolygon */
    INTEGER4 NumNodes  = 8;     /* Number of nodes in the octagon.*/
    INTEGER4 NumElems  = 1;     /* Number of octagonal elements. */
    INTEGER4 NumFaces  = 8;     /* Number of faces in the octagon.*/
    INTEGER4 ICellMax  = 0;     /* Not Used */
    INTEGER4 JCellMax  = 0;     /* Not Used */
    INTEGER4 KCellMax  = 0;     /* Not Used */
    double   SolTime   = 360.0;
    INTEGER4 StrandID  = 0;     /* Static Zone */
    INTEGER4 ParentZn  = 0;     /* No Parent Zone */
    INTEGER4 IsBlock   = 1;     /* Block */
    INTEGER4 NFConns   = 0;
    INTEGER4 FNMode    = 0;


    /* For polygonal zones, the total number of face nodes is equal
     * to twice the number of nodes.  This is because, each face
     * has exactly two nodes.
     */
    INTEGER4 NumFaceNodes    = 2 * NumNodes;
    /* Boundary Faces and Boundary Connections are not used in this
     * example.
     */
    INTEGER4 NumBFaces       = 0;
    INTEGER4 NumBConnections = 0;

    INTEGER4 ShrConn         = 0;

    I = TECZNE142((char*)"Octagonal Zone",
                  &ZoneType,
                  &NumNodes,
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
                  &NumFaceNodes,
                  &NumBFaces,
                  &NumBConnections,
                  NULL,
                  NULL,  /* When Value Location is not specified,
                          * Tecplot will treat all variables as
                          * nodal variables.
                          */
                  NULL,
                  &ShrConn);
    /* DOCEND */

    /* Establish numbering.
     * For this example, we will a single octagonal cell.
     * Before defining your variables, you must establish a
     * consistent node numbering scheme for your  data.  Once the
     * node numbers are defined, supply the variable values in the
     * node numbering order. In this example, node 1 is defined at
     * X = 0 and Y = 0. As such, the first value supplied for X
     * (i.e. X[0]) is 0. Similarly, the first value supplied for
     * Y is 0.
     *
     * It is important that you refer to node numbers consistently.
     * The node numbers will be used later to define the
     * connectivity for each element.
     */

    /* TECDAT Variables */
    /* Set up the variable values.  The variable values will be
     * written to the file using TECDAT.  Because we are specifying
     * nodal variables (as specified via the ValueLocation
     * parameter in TECZNE, each variable is dimensioned by the
     * number of points (NumPts) in the Zone.  You have the option
     * to specify some variables with nodal values and some with
     * cell-centered values. Refer to the Binary Function Reference
     * for details.
     */
    /* DOCSTART:octagon_tecdat.txt*/
    float *X = new float[NumNodes];
    float *Y = new float[NumNodes];
    float *P = new float[NumNodes];

    //Define the grid values.
    X[0] = 0.25;
    Y[0] = 0.0;

    X[1] = 0.75;
    Y[1] = 0.0;

    X[2] = 1.0;
    Y[2] = 0.25;

    X[3] = 1.0;
    Y[3] = 0.75;

    X[4] = 0.75;
    Y[4] = 1.0;

    X[5] = 0.25;
    Y[5] = 1.0;

    X[6] = 0.0;
    Y[6] = 0.75;

    X[7] = 0.0;
    Y[7] = 0.25;

    for (INTEGER4 ii = 0; ii < 8; ii++)
        P[ii] = .5;

    /* Write out the field data using TECDAT */
    INTEGER4 DIsDouble = 0;  /* set IsDouble to 0 to use float
                            * variables. */

    I   = TECDAT142(&NumNodes,  X,  &DIsDouble);
    I   = TECDAT142(&NumNodes,  Y,  &DIsDouble);
    I   = TECDAT142(&NumNodes,  P,  &DIsDouble);

    delete X;
    delete Y;
    delete P;
    /* DOCEND */

    /* Define the Face Nodes.

     * The FaceNodes array is used to indicate which nodes define
     * which face. As mentioned earlier, the number of the nodes is
     * implicitly defined by the order in which the nodal data is
     * provided.  The first value of each nodal variable describes
     * node 1, the second value describes node 2, and so on.
     *
     * The face numbering is implicitly defined.  Because there are
     * two nodes in each face, the first two nodes provided define
     * face 1, the next two define face 2 and so on.  If there was
     * a variable number of nodes used to define the faces, the
     * array would be more complicated.
     */
    /* DOCSTART:octagon_facenodes.txt*/
    INTEGER4 *FaceNodes = new INTEGER4[NumFaceNodes];

    /*
     * Loop over number of sides, and set each side to two
     * consecutive nodes.
     */
    for (INTEGER4 ii = 0; ii < 8; ii++)
    {
        FaceNodes[2*ii] = ii + 1;
        FaceNodes[2*ii+1] = ii + 2;
    }
    FaceNodes[15] = 1;
    /* DOCEND */
    /* Define the right and left elements of each face.

     * The last step for writing out the polyhedral data is to
     * define the right and left neighboring elements for each
     * face.  The neighboring elements can be determined using the
     * right-hand rule.  For each face, place your right-hand along
     * the face which your fingers pointing the direction of
     * incrementing node numbers (i.e. from node 1 to node 2).
     * Your right thumb will point towards the right element; the
     * element on the other side of your hand is the left element.
     *
     * The number zero is used to indicate that there isn't an
     * element on that side of the face.
     *
     * Because of the way we numbered the nodes and faces, the
     * right element for every face is the element itself
     * (element 1) and the left element is "no-neighboring element"
     * (element 0).
     */
    /* DOCSTART:octagon_rightleft.txt*/
    INTEGER4 *FaceLeftElems  = new INTEGER4[NumFaces];
    INTEGER4 *FaceRightElems = new INTEGER4[NumFaces];

    for (INTEGER4 ii = 0; ii < NumFaces; ii++)
    {
        FaceLeftElems[ii]   = 0;
        FaceRightElems[ii]  = 1;
    }
    /* DOCEND */
    /* Write the polyhedral data to the file.  */
    /* DOCSTART:octagon_tecpoly.txt*/
    I = TECPOLYFACE142(&NumFaces,
                       NULL,
                       FaceNodes,
                       FaceLeftElems,
                       FaceRightElems);

    delete FaceNodes;
    delete FaceLeftElems;
    delete FaceRightElems;
    /* DOCEND */
    /* DOCSTART:octagon_tecend.txt*/
    I = TECEND142();
    /* DOCEND */

    return 0;
}
