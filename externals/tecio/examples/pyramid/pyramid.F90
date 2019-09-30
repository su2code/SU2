! This example creates a zone with a single polyhedral cell.

program pyramid
    implicit none
    include 'tecio.f90'

    POINTER(NullPtr, NULL)
    INTEGER(4), dimension(*) :: NULL
    INTEGER(4) :: FileType   = 0   ! 0 for full file
    INTEGER(4) :: FileFormat = 0   ! 0 == PLT, 1 == SZPLT Only PLT is currently supported for polyhedral zones
    INTEGER(4) :: Debug      = 0
    INTEGER(4) :: VIsdouble  = 1
    INTEGER(4) :: I          = 0      ! use to check return codes
    INTEGER(4) :: ZoneType   = 7      ! 7 for FEPolyhedron
    INTEGER(4) :: NumNodes   = 5      ! number of unique nodes
    INTEGER(4) :: NumElems   = 1      ! number of elements
    INTEGER(4) :: NumFaces   = 5      ! number of unique faces
    INTEGER(4) :: ICellMax   = 0      ! Not Used, set to zero
    INTEGER(4) :: JCellMax   = 0      ! Not Used, set to zero
    INTEGER(4) :: KCellMax   = 0      ! Not Used, set to zero
    double precision :: SolTime    = 12.65d0 ! solution time  
    INTEGER(4) :: StrandID   = 0      ! static zone    
    INTEGER(4) :: ParentZone = 0      ! no parent zone 
    INTEGER(4) :: IsBlock    = 1      ! block format 
    INTEGER(4) :: NFConns    = 0      ! not used for FEPolyhedron zones                             
    INTEGER(4) :: FNMode     = 0      ! not used for FEPolyhedron zones
    INTEGER(4) :: ShrConn    = 0
    INTEGER(4) ::   NumFaceNodes = 16
    INTEGER(4) ::   NumBConns    = 0  ! No Boundary Connections
    INTEGER(4) ::   NumBItems    = 0  ! No Boundary Items
    double precision, dimension(5) :: X
    double precision, dimension(5) :: Y
    double precision, dimension(5) :: Z
    INTEGER(4), dimension(5)  :: FaceNodeCounts
    INTEGER(4), dimension(16) :: FaceNodes
    INTEGER(4), dimension(5)  :: FaceRightElems
    INTEGER(4), dimension(5)  :: FaceLeftElems
    
    NullPtr = 0

    I = TECINI142('Pyramid' // char(0), &        ! Data Set Title
                  'X Y Z' // char(0), &          ! Variable List
                  'pyramidf90.plt' // char(0), & ! File Name
                  '.' // char(0), &              ! Scratch Directory
                  FileFormat, &
                  FileType, &
                  Debug, &
                  VIsdouble)

    ! The number of face nodes in the zone.  This example creates
    ! a zone with a single pyramidal cell.  This cell has four
    ! triangular faces and one rectangular face, yielding a total
    ! of 16 face nodes.

    I = TECZNE142('Polyhedral Zone (Octahedron)'//char(0), &
                  ZoneType, &
                  NumNodes, &
                  NumElems, &
                  NumFaces, &
                  ICellMax, &
                  JCellMax, &
                  KCellMax, &
                  SolTime, &
                  StrandID, &
                  ParentZone, &
                  IsBlock, &
                  NFConns, &
                  FNMode, &
                  NumFaceNodes, &
                  NumBConns, &
                  NumBItems, &
                  NULL, & ! PassiveVarArray
                  NULL, & ! ValueLocArray
                  NULL, & ! VarShareArray
                  ShrConn)

    ! Initialize arrays of nodal data
    X(1) = 0
    Y(1) = 0
    Z(1) = 0

    X(2) = 1
    Y(2) = 1
    Z(2) = 2

    X(3) = 2
    Y(3) = 0
    Z(3) = 0

    X(4) = 2
    Y(4) = 2
    Z(4) = 0

    X(5) = 0
    Y(5) = 2
    Z(5) = 0

    ! Write the data.
    I = TECDATD142(NumNodes, X)
    I = TECDATD142(NumNodes, Y)
    I = TECDATD142(NumNodes, Z)

    ! Define the Face Nodes.
    !
    ! The FaceNodes array is used to indicate which nodes define
    ! which face. As mentioned earlier, the number of the nodes is
    ! implicitly defined by the order in which the nodal data is
    ! provided.  The first value of each nodal variable describes
    ! node 1, the second value describes node 2, and so on.
    !
    ! The face numbering is implicitly defined.  Because there are
    ! two nodes in each face, the first two nodes provided define
    ! face 1, the next two define face 2 and so on.  If there was
    ! a variable number of nodes used to define the faces, the
    ! array would be more complicated.

    ! The first four faces are triangular, i.e. have three nodes.
    ! The fifth face is rectangular, i.e. has four nodes.
    FaceNodeCounts(1) = 3
    FaceNodeCounts(2) = 3
    FaceNodeCounts(3) = 3
    FaceNodeCounts(4) = 3
    FaceNodeCounts(5) = 4

    ! Face Nodes for Face 1
    FaceNodes(1)  = 1
    FaceNodes(2)  = 2
    FaceNodes(3)  = 3

    ! Face Nodes for Face 2
    FaceNodes(4)  = 3
    FaceNodes(5)  = 2
    FaceNodes(6)  = 4

    ! Face Nodes for Face 3
    FaceNodes(7)  = 5
    FaceNodes(8)  = 2
    FaceNodes(9)  = 4

    ! Face Nodes for Face 4
    FaceNodes(10)  = 1
    FaceNodes(11) = 2
    FaceNodes(12) = 5

    ! Face Nodes for Face 5
    FaceNodes(13) = 1
    FaceNodes(14) = 5
    FaceNodes(15) = 4
    FaceNodes(16) = 3

    ! Define the right and left elements of each face.
    !
    ! The last step for writing out the polyhedral data is to
    ! define the right and left neighboring elements for each
    ! face.  The neighboring elements can be determined using the
    ! right-hand rule.  For each face, place your right-hand along
    ! the face which your fingers pointing the direction of
    ! incrementing node numbers (i.e. from node 1 to node 2).
    ! Your right thumb will point towards the right element the
    ! element on the other side of your hand is the left element.
    !
    ! The number zero is used to indicate that there isn't an
    ! element on that side of the face.
    !
    ! Because of the way we numbered the nodes and faces, the
    ! right element for every face is the element itself
    ! (element 1) and the left element is "no-neighboring element"
    ! (element 0).

    FaceLeftElems(1) = 1
    FaceLeftElems(2) = 1
    FaceLeftElems(3) = 0
    FaceLeftElems(4) = 0
    FaceLeftElems(5) = 0

    FaceRightElems(1) = 0
    FaceRightElems(2) = 0
    FaceRightElems(3) = 1
    FaceRightElems(4) = 1
    FaceRightElems(5) = 1
    
    ! Write the face map (created above) using TECPOLYFACE142.
    I = TECPOLYFACE142(NumFaces, &
                       FaceNodeCounts, &  ! The face node counts array
                       FaceNodes, &       ! The face nodes array
                       FaceLeftElems, &   ! The left elements array 
                       FaceRightElems)    ! The right elements array

    I = TECEND142()
end program pyramid
! DOCEND
