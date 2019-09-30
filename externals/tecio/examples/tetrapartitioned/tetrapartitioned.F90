
program tetrapartitioned
    implicit none
#if defined TECIOMPI
    include 'mpif.h'
#endif
    !
    ! Example C++ program to write a partitioned
    ! binary datafile for tecplot. This example
    ! does the following:
    !
    !   1.  Open a datafile called "tetrapartitioned.szplt"
    !   2.  Assign values for x, y, z and p.
    !   3.  Write out a tetrahedral zone in 3 partitions.
    !   4.  Close the datafile.
    !
    ! If TECIOMPI is #defined, this program may be executed with mpiexec with
    ! up to 3 MPI ranks (processes). In this case, it must be linked
    ! with the MPI version of TECIO.
    !

    include "tecio.f90"

    ! The zone will appear the same as an ordered zone with these dimensions:
    integer(4), parameter :: XDIM = 10
    integer(4), parameter :: YDIM = 9
    integer(4), parameter :: ZDIM = 8

    integer(4), parameter :: FULL = 0
    integer(4), parameter :: GRID = 1
    integer(4), parameter :: SOLUTION = 2

    !
    ! Open the file and write the tecplot datafile
    ! header information
    !
    integer(4) :: fileFormat = 1 ! 1 == SZPLT (cannot output partitioned zones to PLT)
    integer(4) :: fileType   = FULL
    integer(4) :: debug = 1
    integer(4) :: vIsDouble = 0
    integer(4) :: returnValue = 0;
    integer(4) :: zoneType      = 4                   ! Tetrahedral
    integer(4) :: nNodes        = XDIM * YDIM * ZDIM  ! Overall zone dimensions
    integer(4) :: nCells        = (XDIM - 1) * (YDIM - 1) * (ZDIM - 1) * 5
    integer(4) :: nFaces        = 4                   ! Not used
    integer(4) :: iCellMax      = 0
    integer(4) :: jCellMax      = 0
    integer(4) :: kCellMax      = 0
    double precision :: solTime = 360.0
    integer(4) :: strandID      = 0                   ! StaticZone
    integer(4) :: parentZn      = 0                   ! No Parent
    integer(4) :: isBlock       = 1                   ! Block
    integer(4) :: nFConns       = 0
    integer(4) :: fNMode        = 0
    integer(4), dimension(4) :: valueLocations = (/1, 1, 1, 0/)    ! 1 = Nodal, 0 = Cell-Centered
    integer(4) :: shrConn   = 0
    integer(4) :: zone = 1
    integer(4) :: dIsDouble = 0
    integer(4), dimension(3) :: iMin, iMax, jMin, jMax, kMin, kMax
    integer(4), dimension(3) :: iDim, jDim, kDim
    integer(4), dimension(3) :: pNNodes, pNCells ! Partition node and cell counts, including ghost items
    integer(4) :: ptn, connectivityCount
    integer(4), dimension(:, :), allocatable :: connectivity
    real(4), dimension(:, :), allocatable ::  x, y, z, p
    integer :: allocateStatus
    integer :: maxCells, maxNodes
    integer :: i, j, k, index
    integer(4), dimension(3) :: nGNodes, nGCells ! Partition ghost node and ghost cell counts
    integer(4), dimension(:, :), allocatable, target :: ghostNodes, gNPartitions, gNPNodes, ghostCells
    ! Upper bound on counts of ghost nodes and cells
    integer :: maxGCells = XDIM * ZDIM + YDIM * ZDIM * 5
    integer :: maxGNodes = 2 * (XDIM * ZDIM + YDIM * ZDIM)
    POINTER   (NullPtr,Null)
    integer(4) Null(*)
    integer(4), dimension(8) :: brickNodes
    integer, dimension(4, 5, 2) :: tetraCorners
    integer :: evenOdd, whichTet, whichCorner
    
    interface
        subroutine GatherGhostNodesAndCells( &
            nGNodes, ghostNodes, gNPartitions, gNPNodes, &
            nGCells, ghostCells, iDim, jDim, kDim, jMin, jMax, maxGNodes, maxGCells)
            implicit none
            integer(4), dimension(maxGNodes,3), target :: ghostNodes, gNPartitions, gNPNodes
            integer(4), dimension(maxGCells, 3), target :: ghostCells
            integer(4) :: nGNodes(3), nGCells(3), iDim(3), jDim(3), kDim(3), jMin(3), jMax(3)
            integer(4) :: maxGNodes, maxGCells
        end subroutine GatherGhostNodesAndCells
    end interface

#if defined TECIOMPI
    integer(4) :: commSize
    integer(4) :: commRank
    integer(4) :: mainRank = 0
    integer(4) :: numPartitions = 3
    integer(4), dimension(3) :: partitionOwners
    integer(4) :: ierr
    INTEGER(4) :: mpiComm = MPI_COMM_WORLD
    call MPI_Init(ierr)
    call MPI_Comm_size(mpiComm, commSize, ierr)
    call MPI_Comm_rank(mpiComm, commRank, ierr)
#endif
    returnValue = TECINI142("SIMPLE DATASET"//char(0), &
                  "x y z p"//char(0), &
                  "tetrapartitionedf90.plt"//char(0), &
                  "."//char(0), &
                  fileFormat, &
                  fileType, &
                  debug, &
                  vIsDouble)
#if defined TECIOMPI
    if (returnValue .eq. 0) then
        returnValue = TECMPIINIT142(mpiComm, mainRank)
    endif
#endif

    /*
     * Zone
     */
    NullPtr = 0
    if (returnValue .eq. 0) &
        returnValue = TECZNE142("Partitioned Zone"//char(0), &
                      zoneType, &
                      nNodes, &
                      nCells, &
                      nFaces, &
                      iCellMax, &
                      jCellMax, &
                      kCellMax, &
                      solTime, &
                      strandID, &
                      parentZn, &
                      isBlock, &
                      nFConns, &
                      fNMode, &
                      0, &              ! TotalNumFaceNodes
                      0, &              ! NumConnectedBoundaryFaces
                      0, &              ! TotalNumBoundaryConnections
                      NULL, &           ! PassiveVarList
                      valueLocations, &
                      NULL, &           ! SharVarFromZone
                      shrConn)

    ! Divide the zone into 3 partitions, identified by the index ranges
    ! of an equivalent unpartitioned ordered zone.

    ! Partition 1 node range, which will include one layer of ghost cells on the IMAX boundary:
    iMin(1) = 0 ! We use zero-based indices here because it simplifies index arithmetic in C.
    iMax(1) = XDIM / 2 + 2
    jMin(1) = 0
    jMax(1) = YDIM
    kMin(1) = 0
    kMax(1) = ZDIM

    ! Partition 2; ghost cells on IMIN and JMAX boundaries:
    iMin(2) = iMax(1) - 3
    iMax(2) = XDIM
    jMin(2) = jMin(1);
    jMax(2) = YDIM / 2 + 2
    kMin(2) = kMin(1)
    kMax(2) = kMax(1)

    ! Partition 3; ghost cells on IMIN and JMIN boundaries:
    iMin(3) = iMin(2)
    iMax(3) = iMax(2)
    jMin(3) = jMax(2) - 3
    jMax(3) = YDIM
    kMin(3) = kMin(2)
    kMax(3) = kMax(2)

    ! Local partition dimensions (of equivalent ordered zones)
    do i = 1, 3
        iDim(i) = iMax(i) - iMin(i)
        jDim(i) = jMax(i) - jMin(i)
        kDim(i) = kMax(i) - kMin(i)
    enddo

    ! Allocate memory for connectivity and variable values
    do ptn = 1, 3
        pNNodes(ptn) = iDim(ptn) * jDim(ptn) * kDim(ptn)
        pNCells(ptn) = (iDim(ptn) - 1) * (jDim(ptn) - 1) * (kDim(ptn) - 1) * 5
    enddo
    
    maxCells = max(pNCells(1), pNCells(2), pNCells(3))
    maxNodes = max(pNNodes(1), pNNodes(2), pNNodes(3))
    allocate(connectivity(4 * maxCells, 3), stat = allocateStatus)
    if (allocateStatus .eq. 0) allocate(x(maxNodes, 3), stat = allocateStatus)
    if (allocateStatus .eq. 0) allocate(y(maxNodes, 3), stat = allocateStatus)
    if (allocateStatus .eq. 0) allocate(z(maxNodes, 3), stat = allocateStatus)
    if (allocateStatus .eq. 0) allocate(p(maxCells, 3), stat = allocateStatus)
    if (allocateStatus .ne. 0) stop 'Unable to allocate memory'
    
    ! Calculate variable and connectivity values for partitions.
    do ptn = 1, 3
        ! Nodes
        do k = 0, kDim(ptn) - 1
            do j = 0, jDim(ptn) - 1
                do i = 0, iDim(ptn) - 1
                    index = (k * jDim(ptn) + j) * iDim(ptn) + i + 1;
                    x(index, ptn) = real(i + iMin(ptn) + 1)
                    y(index, ptn) = real(j + jMin(ptn) + 1)
                    z(index, ptn) = real(k + kMin(ptn) + 1)
                enddo
            enddo
        enddo
        ! p (cell-centered) and connectivity
        do k = 0, kDim(ptn) - 2
            do j = 0, jDim(ptn) - 2
                do i = 0, iDim(ptn) - 2
                    !
                    ! Start with a brick and divide it into 5 tets.
                    ! Need to mirror the subdivision for neighboring
                    ! bricks to avoid internal faces.
                    !
                    tetraCorners = reshape((/ &
                        ! "Even" bricks
                        1,2,3,6, &
                        1,3,4,8, &
                        5,8,6,1, &
                        6,8,7,3, &
                        1,3,8,6, &
                        ! "Odd" bricks
                        1,2,4,5, &
                        2,3,4,7, &
                        5,7,6,2, &
                        5,8,7,4, &
                        2,4,5,7/), shape(tetraCorners))
                    evenOdd = MOD(i + iMin(ptn) + j + jMin(ptn) + k + kMin(ptn), 2) + 1

                    brickNodes(1) = (k * jDim(ptn) + j) * iDim(ptn) + i + 1
                    brickNodes(2) = brickNodes(1) + 1
                    brickNodes(3) = brickNodes(1) + iDim(ptn) + 1
                    brickNodes(4) = brickNodes(1) + iDim(ptn)
                    brickNodes(5) = brickNodes(1) + iDim(ptn) * jDim(ptn)
                    brickNodes(6) = brickNodes(2) + iDim(ptn) * jDim(ptn)
                    brickNodes(7) = brickNodes(3) + iDim(ptn) * jDim(ptn)
                    brickNodes(8) = brickNodes(4) + iDim(ptn) * jDim(ptn)

                    do whichTet = 1, 5
                        index = ((k * (jDim(ptn) - 1) + j) * (iDim(ptn) - 1) + i) * 5 + whichTet
                        p(index, ptn) = real((i + iMin(ptn) + 1) * (j + jMin(ptn) + 1) * (k + kMin(ptn) + 1))
                        do whichCorner = 1, 4
                            connectivity((index - 1) * 4 + whichCorner, ptn) = &
                                brickNodes(tetraCorners(whichCorner, whichTet, evenOdd))
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    allocate(ghostNodes(maxGNodes, 3), stat = allocateStatus)
    if (allocateStatus .eq. 0) allocate(gNPartitions(maxGNodes, 3), stat = allocateStatus)
    if (allocateStatus .eq. 0) allocate(gNPNodes(maxGNodes, 3), stat = allocateStatus)
    if (allocateStatus .eq. 0) allocate(ghostCells(maxGCells, 3), stat = allocateStatus)
    if (allocateStatus .ne. 0) stop 'Unable to allocate memory for ghost items'
    call GatherGhostNodesAndCells(nGNodes, ghostNodes, gNPartitions, gNPNodes, nGCells, ghostCells, &
        iDim, jDim, kDim, jMin, jMax, maxGNodes, maxGCells)

    ! Output partitions
#if defined TECIOMPI
    do ptn = 0, 2
        partitionOwners(ptn + 1) = mod(ptn, commSize)
    enddo
    if (returnValue .eq. 0) &
        returnValue = TECZNEMAP142(numPartitions, partitionOwners)
#endif

    do ptn = 1, 3
#if defined TECIOMPI
        if (partitionOwners(ptn) .eq. commRank) then
#endif
        if (returnValue .eq. 0) &
            returnValue = TECFEPTN142(ptn,                  &
                                      pNNodes(ptn),         &
                                      pNCells(ptn),         &
                                      nGNodes(ptn),         &
                                      ghostNodes(1, ptn),   &
                                      gNPartitions(1, ptn), &
                                      gNPNodes(1, ptn),     &
                                      nGCells(ptn),         &
                                      ghostCells(1, ptn))
        /*
         * Write out the field data.
         */
        if (returnValue == 0) returnValue = TECDAT142(pNNodes(ptn), x(1, ptn), dIsDouble)
        if (returnValue == 0) returnValue = TECDAT142(pNNodes(ptn), y(1, ptn), dIsDouble)
        if (returnValue == 0) returnValue = TECDAT142(pNNodes(ptn), z(1, ptn), dIsDouble)
        if (returnValue == 0) returnValue = TECDAT142(pNCells(ptn), p(1, ptn), dIsDouble)

        connectivityCount = 4 * pNCells(ptn)
        if (returnValue == 0) returnValue = TECNODE142(connectivityCount, connectivity(1, ptn))

#if defined TECIOMPI
        endif
#endif
    enddo

    if (returnValue == 0) returnValue = TECEND142()

    deallocate(ghostNodes, gNPartitions, gNPNodes, ghostCells, x, y, z, p)

#if defined TECIOMPI
    call MPI_Finalize(ierr)
#endif

end program tetrapartitioned

subroutine appendGhostNodes( &
    nGhosts, ghosts, gPartitions, gPGhosts, ptn, gIDim, gJDim, &
    gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd, oIDim, oJDim, &
    oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd)
    
    integer(4) :: nGhosts, ptn
    integer(4), dimension(:) :: ghosts, gPartitions, gPGhosts
    integer(4) :: gIDim, gJDim, gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd
    integer(4) :: oIDim, oJDim, oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd

    integer :: i, j, k
    integer(4) oI, oJ, oK
    do i = gIStart, gIEnd - 1
        do j = gJStart, gJEnd - 1
            do k = gKStart, gKEnd - 1
                nGhosts = nGhosts + 1
                ghosts(nGhosts) = (k * gJDim + j) * gIDim + i + 1
                gPartitions(nGhosts) = ptn
                oI = i - gIStart + oIStart
                oJ = j - gJStart + oJStart
                oK = k - gKStart + oKStart
                gPGhosts(nGhosts) = (oK * oJDim + oJ) * oIDim + oI + 1;
            enddo
        enddo
    enddo
end

subroutine appendGhostCells( &
    nGhosts, ghosts, gIDim, gJDim, &
    gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd)
    integer(4) nGhosts
    integer(4), dimension(:) :: ghosts
    integer(4) :: gIDim, gJDim, gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd

    integer :: i, j, k, whichTet

    do i = gIStart, gIEnd - 1
        do j = gJStart, gJEnd - 1
            do k = gKStart, gKEnd - 1
                do whichTet = 1, 5
                    nGhosts = nGhosts + 1
                    ghosts(nGhosts) = ((k * gJDim + j) * gIDim + i) * 5 + whichTet
                enddo
            enddo
        enddo
    enddo
end


subroutine GatherGhostNodesAndCells( &
    nGNodes, ghostNodes, gNPartitions, gNPNodes, &
    nGCells, ghostCells, iDim, jDim, kDim, jMin, jMax, maxGNodes, maxGCells)
    implicit none
    integer(4), dimension(maxGNodes,3), target :: ghostNodes, gNPartitions, gNPNodes
    integer(4), dimension(maxGCells, 3), target :: ghostCells
    integer(4) :: nGNodes(3), nGCells(3), iDim(3), jDim(3), kDim(3), jMin(3), jMax(3)
    integer(4) :: maxGNodes, maxGCells
    !
    ! Assemble lists of ghost nodes and cells--nodes and cells near partition boundaries that
    ! coincide with those "owned" by neighboring partitions. For each set of coincident nodes
    ! or cells, exactly one partition must own the node or cell, and all other involved partitions
    ! must report it as a ghost node or cell.
    !
    ! Arbitrarily, we say that the first partition owns any nodes that do not overlay the interior
    ! of neighboring partitions. That is, it owns any nodes that its "real" (non-ghost) cells use.
    ! So only a single layer of nodes and cells--on its IMax boundary--are ghosts, and are owned
    ! by the second and third partitions. We use the same logic to assign ownership for nodes
    ! shared by partitions 2 and 3.
    !
    integer(4), pointer, dimension(:) :: gnPtr, gnpPtr, gnpnPtr, gcPtr
    POINTER   (NullPtr,Null)
    Integer(4) Null(1)
    interface
        subroutine appendGhostNodes( &
        nGhosts, ghosts, gPartitions, gPGhosts, ptn, gIDim, gJDim, &
        gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd, oIDim, oJDim, &
        oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd)
    
        integer(4) :: nGhosts, ptn
        integer(4), dimension(:) :: ghosts, gPartitions, gPGhosts
        integer(4) :: gIDim, gJDim, gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd
        integer(4) :: oIDim, oJDim, oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd
        end subroutine appendGhostNodes
        
        subroutine appendGhostCells( &
        nGhosts, ghosts, gIDim, gJDim, &
        gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd)
        
        integer(4) nGhosts
        integer(4), dimension(:) :: ghosts
        integer(4) :: gIDim, gJDim, gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd
        end subroutine appendGhostCells
    end interface
    nGNodes(1) = 0
    nGCells(1) = 0
    gnPtr => ghostNodes(1:maxGNodes, 1)
    gnpPtr => gnPartitions(1:maxGNodes, 1)
    gnpnPtr => gnpNodes(1:maxGNodes, 1)
    gcPtr => ghostCells(1:maxGCells, 1)
    ! Nodes owned by the second partition:
    call appendGhostNodes( &
        nGNodes(1), gnPtr, gnpPtr, gnpnPtr, 2, &
        iDim(1), jDim(1), &                                 ! I- and J-dimensions
        iDim(1) - 1, iDim(1), 0, jDim(2) - 1, 0, kDim(1), & ! local index ranges
        iDim(2), jDim(2), &                                 ! I- and J-dimensions
        2, 3, 0, jDim(2) - 1, 0, kDim(2))                   ! local index ranges
    ! Nodes owned by the third partition:
    call appendGhostNodes( &
        nGNodes(1), gnPtr, gnpPtr, gnpnPtr, 3, &
        iDim(1), jDim(1), &                                       ! I- and J-dimensions
        iDim(1) - 1, iDim(1), jMin(3) + 2, jMax(3), 0, kDim(1), & ! local index ranges
        iDim(3), jDim(3), &                                       ! I- and J-dimensions
        2, 3, 0, jDim(3), 0, kDim(3))                             ! local index ranges
    ! Cells owned by the second partition:
    call appendGhostCells( &
        nGCells(1), gcptr, &
        iDim(1) - 1, jDim(1) - 1, &                               ! I- and J-dimensions
        iDim(1) - 2, iDim(1) - 1, 0, jDim(2) - 2, 0, kDim(1) - 1) ! local index ranges
    ! Cells owned by the third partition:
    call appendGhostCells( &
        nGCells(1), gcptr, &
        iDim(1) - 1, jDim(1) - 1, &                                         ! I- and J-dimensions
        iDim(1) - 2, iDim(1) - 1, jDim(2) - 2, jDim(1) - 1, 0, kDim(1) - 1) ! local index ranges

    nGNodes(2) = 0
    nGCells(2) = 0
    gnPtr => ghostNodes(1:maxGNodes, 2)
    gnpPtr => gnPartitions(1:maxGNodes, 2)
    gnpnPtr => gnpNodes(1:maxGNodes, 2)
    gcPtr => ghostCells(1:maxGCells, 2)
    ! Nodes owned by the first partition.
    call appendGhostNodes( &
        nGNodes(2), gnPtr, gnpPtr, gnpnPtr, 1, &
        iDim(2), jDim(2), &                               ! I- and J-dimensions
        0, 2, 0, jDim(2), 0, kDim(2), &                   ! local index ranges
        iDim(1), jDim(1), &                               ! I- and J-dimensions
        iDim(1) - 3, iDim(1) - 1, 0, jDim(2), 0, kDim(1)) ! local index ranges
    ! Nodes owned by the third partition.
    call appendGhostNodes( &
        nGNodes(2), gnPtr, gnpPtr, gnpnPtr, 3, &
        iDim(2), jDim(2), &                             ! I- and J-dimensions
        2, iDim(2), jDim(2) - 1, jDim(2), 0, kDim(2), & ! local index ranges
        iDim(3), jDim(3), &                             ! I- and J-dimensions
        2, iDim(3), 2, 3, 0, kDim(3))                   ! local index ranges
    ! Cells owned by the first partition.
    call appendGhostCells( &
        nGCells(2), gcptr, &
        iDim(2) - 1, jDim(2) - 1, &           ! I- and J-dimensions
        0, 1, 0, jDim(2) - 1, 0, kDim(2) - 1) ! local index ranges
    ! Cells owned by the third partition.
    call appendGhostCells( &
        nGCells(2), gcptr, &
        iDim(2) - 1, jDim(2) - 1, &                               ! I- and J-dimensions
        1, iDim(2) - 1, jDim(2) - 2, jDim(2) - 1, 0, kDim(2) - 1) ! local index ranges

    nGNodes(3) = 0
    nGCells(3) = 0
    gnPtr => ghostNodes(1:maxGNodes, 3)
    gnpPtr => gnPartitions(1:maxGNodes, 3)
    gnpnPtr => gnpNodes(1:maxGNodes, 3)
    gcPtr => ghostCells(1:maxGCells, 3)
    ! Nodes owned by the first partition
    call appendGhostNodes( &
        nGNodes(3), gnPtr, gnpPtr, gnpnPtr, 1, &
        iDim(3), jDim(3), &                                     ! I- and J-dimensions
        0, 2, 0, jDim(3), 0, kDim(3), &                         ! local index ranges
        iDim(1), jDim(1), &                                     ! I- and J-dimensions
        iDim(1) - 3, iDim(1) - 1, jMin(3), jMax(3), 0, kDim(1)) ! local index ranges
    ! Nodes owned by the second partition.
    call appendGhostNodes( &
        nGNodes(3), gnPtr, gnpPtr, gnpnPtr, 2, &
        iDim(3), jDim(3), &                               ! I- and J-dimensions
        2, iDim(3), 0, 2, 0, kDim(3), &                   ! local index ranges
        iDim(2), jDim(2), &                               ! I- and J-dimensions
        2, iDim(2), jDim(2) - 3, jDim(2) - 1, 0, kDim(2)) ! local index ranges
    ! Cells owned by the first partition
    call appendGhostCells( &
        nGCells(3), gcPtr, &
        iDim(3) - 1, jDim(3) - 1, &           ! I- and J-dimensions
        0, 1, 0, jDim(3) - 1, 0, kDim(3) - 1) ! local index ranges
    ! Nodes owned by the second partition.
    call appendGhostCells( &
        nGCells(3), gcPtr, &
        iDim(3) - 1, jDim(3) - 1, &           ! I- and J-dimensions
        1, iDim(3) - 1, 0, 1, 0, kDim(3) - 1) ! local index ranges
end
