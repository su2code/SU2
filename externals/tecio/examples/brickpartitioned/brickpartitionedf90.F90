program brickpartitioned
    implicit none
#if defined TECIOMPI
    include "mpif.h"
#endif
    ! internal testing
    ! RUNFLAGS:none
    !
    ! Example FORTRAN program to write a partitioned binary grid and solution
    ! SZPLT data file for Tecplot. This example does the following:
    !   1.  Open a datafile called "brickpartitioned_grid.szplt"
    !       and "brickpartitioned_solution.szplt"
    !   2.  Assign values for x, y, z to the grid and p to the solution.
    !   3.  Write out a hexahedral (brick) zone in 3 partitions.
    !   4.  Close the datafiles
    !
    ! If TECIOMPI is #defined, this program may be executed with mpiexec with
    ! up to 3 MPI ranks (processes). In this case, it must be linked
    ! with the MPI version of TECIO.
    !

    include "tecio.f90"

    interface
        integer(4) function initializeFile( &
#           if defined TECIOMPI
            mpiComm, &
#           endif
            fOffset, numSolutionFiles, numPartitions, solTime)
            implicit none
#           if defined TECIOMPI
            integer(4) :: mpiComm
#           endif
            integer(4) :: fOffset
            integer(4) :: numSolutionFiles
            integer(4) :: numPartitions
            double precision :: solTime
        end function initializeFile

        integer(4) function createZone( &
            fOffset, numPartitions, partitionOwners, solTime)
            implicit none
            integer(4) :: fOffset
            integer(4) :: numPartitions
            integer(4), dimension(3) :: partitionOwners
            double precision :: solTime
        end function createZone

        integer(4) function createData( &
#           if defined TECIOMPI
            commRank, &
#           endif
            fOffset, numSolutionFiles, numPartitions, partitionOwners, solTime)
            implicit none
#           if defined TECIOMPI
            integer(4) :: commRank
#           endif
            integer(4) :: fOffset
            integer(4) :: numSolutionFiles
            integer(4) :: numPartitions
            integer(4), dimension(3) :: partitionOwners
            double precision :: solTime
        end function createData

        integer(4) function finalizeFile()
        end function finalizeFile
    end interface

    integer(4) :: returnValue = 0

    integer(4), dimension(3) :: partitionOwners
    integer(4) :: numPartitions = 3
    double precision :: solTime

    ! Set numSolutionFiles to 0 for a FULL file or 1 to 9 for grid/solution file(s).
    ! The solution file names will be decorated with the time step. A value of 0
    ! creates a single file, brickpartitioned.szplt, while a value of 1 or greater creates
    ! two or more files, brickpartitioned_grid.szplt and brickpartitioned_solution1.szplt,
    ! brickpartitioned_solution2.szplt, etc.
    integer(4) :: numSolutionFiles = 3

    integer(4) :: numOutputFiles
    integer(4) :: fOffset

#   if defined TECIOMPI
    integer(4) :: ptn
    integer(4) :: commSize
    integer(4) :: commRank
    integer(4) :: mainRank = 0
    integer(4) :: ierr
    INTEGER(4) :: mpiComm = MPI_COMM_WORLD
    call MPI_Init(ierr)
    call MPI_Comm_size(mpiComm, commSize, ierr)
    call MPI_Comm_rank(mpiComm, commRank, ierr)
    do ptn = 0, 2
        partitionOwners(ptn + 1) = mod(ptn, commSize)
    enddo
#   endif

    numOutputFiles = 1 + numSolutionFiles
    
    do fOffset = 0, numOutputFiles-1
        if (fOffset .eq. 0) then
            solTime = 360
        else
            solTime = 360 + fOffset-1
        endif

        if (returnValue .eq. 0) then
            returnValue = initializeFile( &
#               if defined TECIOMPI
                mpiComm, &
#               endif
                fOffset, numSolutionFiles, numPartitions, solTime)
        endif

        ! Create zone
        if (returnValue .eq. 0) then
            returnValue = createZone(fOffset, numPartitions, partitionOwners, solTime)
        endif

        ! Create the connectivity and variable data
        if (returnValue .eq. 0) then
            returnValue = createData( &
#               if defined TECIOMPI
                commRank, &
#               endif
                fOffset, numSolutionFiles, numPartitions, partitionOwners, solTime)
        endif

        if (returnValue .eq. 0) then
            returnValue = finalizeFile()
        endif
    enddo
    
#   if defined TECIOMPI
    call MPI_Finalize(ierr)
#   endif

end program brickpartitioned

integer(4) function initializeFile( &
#   if defined TECIOMPI
    mpiComm, &
#   endif
    fOffset, numSolutionFiles, numPartitions, solTime)
    implicit none
#   if defined TECIOMPI
    integer(4) :: mpiComm
#   endif
    integer(4) :: fOffset
    integer(4) :: numSolutionFiles
    integer(4) :: numPartitions
    double precision :: solTime

    integer(4) :: fileFormat = 1
    integer(4) :: debug = 1
    integer(4) :: vIsDouble = 0
    character(128) :: fileName
    character(24) :: varNames
    integer(4) :: dataFileType
#   if defined TECIOMPI
    integer(4) :: mainRank = 0
#   endif

    integer(4), parameter :: FULL = 0
    integer(4), parameter :: GRID = 1
    integer(4), parameter :: SOLUTION = 2

    include "tecio.f90"

    initializeFile = 0

    if (numSolutionFiles .eq. 0) then
        dataFileType = FULL
        if (numPartitions .gt. 0) then
            fileName = "brickpartitionedf90.szplt"//char(0)
        else
            fileName = "brickf90.szplt"//char(0)
        endif
        varNames = "x y z p"//char(0)
    else if (fOffset .eq. 0) then
        dataFileType = GRID
        if (numPartitions .gt. 0) then
            fileName = "brickpartitionedf90_grid.szplt"//char(0)
        else
            fileName = "brickf90_grid.szplt"//char(0)
        endif
        varNames = "x y z"//char(0)
    else
        dataFileType = SOLUTION
        if (numPartitions .gt. 0) then
            fileName = ' '
            write(fileName, fmt='(a,i1,a)') 'brickpartitionedf90_solution', &
                foffset, '.szplt'
            fileName = trim(fileName) // char(0)
        else
            fileName = "brickf90_solution.szplt"//char(0)
        endif
        varNames = "p"//char(0)
    endif

    initializeFile = TECINI142("SIMPLE DATASET"//char(0), &
                  varNames, &
                  fileName, &
                  "."//char(0), &
                  fileFormat, &
                  dataFileType, &
                  debug, &
                  vIsDouble)

#   if defined TECIOMPI
        if (initializeFile .eq. 0) then
            initializeFile = TECMPIINIT142(mpiComm, mainRank)
        endif
#   endif

    return
end

integer(4) function createZone( &
    fOffset, numPartitions, partitionOwners, solTime)
    implicit none
    integer(4) :: fOffset
    integer(4) :: numPartitions
    integer(4), dimension(3) :: partitionOwners
    double precision :: solTime

    ! These should really be either declared globally or passed to the function
    integer(4), parameter :: XDIM = 10
    integer(4), parameter :: YDIM = 9
    integer(4), parameter :: ZDIM = 8

    integer(4) :: zoneType  = 5      ! Brick
    integer(4) :: nNodes    = XDIM * YDIM * ZDIM ! Overall zone dimensions
    integer(4) :: nCells    = (XDIM - 1) * (YDIM - 1) * (ZDIM - 1)
    integer(4) :: nFaces    = 6     ! Not used
    integer(4) :: iCellMax  = 0
    integer(4) :: jCellMax  = 0
    integer(4) :: kCellMax  = 0
    integer(4) :: strandID  = 1
    integer(4) :: parentZn  = 0     ! No Parent
    integer(4) :: isBlock   = 1     ! Block
    integer(4) :: nFConns   = 0
    integer(4) :: fNMode    = 0
    integer(4), dimension(4) :: gridAndSolValLocs = (/1, 1, 1, 0/) ! 1 = Nodal, 0 = Cell-Centered
    integer(4) :: shrConn   = 0
    character(1024) :: zoneTitle
    integer(4) baseVarOffset
    POINTER   (NullPtr,Null)
    Integer(4) Null(*)

    include "tecio.f90"

    createZone = 0

    NullPtr = 0
    if (numPartitions .gt. 0) then
        zoneTitle = "partitioned Zone"//char(0)
    else
        zoneTitle = "Zone"//char(0)
    endif

    if (fOffset .eq. 0) then
        createZone = TECZNE142( &
            zoneTitle, zoneType, nNodes, nCells, nFaces, iCellMax, jCellMax, kCellMax, &
            solTime, strandID, parentZn, isBlock, nFConns, fNMode, &
            0, &              ! TotalNumFaceNodes
            0, &              ! NumConnectedBoundaryFaces
            0, &              ! TotalNumBoundaryConnections
            NULL, &           ! PassiveVarList
            gridAndSolValLocs(1), &
            NULL, &           ! SharVarFromZone
            shrConn)
    else
        createZone = TECZNE142( &
            zoneTitle, zoneType, nNodes, nCells, nFaces, iCellMax, jCellMax, kCellMax, &
            solTime, strandID, parentZn, isBlock, nFConns, fNMode, &
            0, &              ! TotalNumFaceNodes
            0, &              ! NumConnectedBoundaryFaces
            0, &              ! TotalNumBoundaryConnections
            NULL, &           ! PassiveVarList
            gridAndSolValLocs(4), &
            NULL, &           ! SharVarFromZone
            shrConn)
    endif

#   if defined TECIOMPI
        ! Output partitions
        if (createZone .eq. 0 .and. numPartitions .gt. 0) then
            createZone = TECZNEMAP142(numPartitions, partitionOwners)
        endif
#   endif

    return
end

integer(4) function createData( &
#   if defined TECIOMPI
    commRank, &
#   endif
    fOffset, numSolutionFiles, numPartitions, partitionOwners, solTime)
    implicit none
#   if defined TECIOMPI
    integer(4) :: commRank
#   endif
    integer(4) :: fOffset
    integer(4) :: numSolutionFiles
    integer(4) :: numPartitions
    integer(4), dimension(3) :: partitionOwners
    double precision :: solTime

    ! These should really be either declared globally or passed to the function
    integer(4), parameter :: XDIM = 10
    integer(4), parameter :: YDIM = 9
    integer(4), parameter :: ZDIM = 8

    integer(4), dimension(3) :: iMin, iMax, jMin, jMax, kMin, kMax
    integer(4), dimension(3) :: iDim, jDim, kDim
    integer(4), dimension(3) :: pNNodes, pNCells ! Partition node and cell counts, including ghost items
    integer(4) :: ptn, connectivityCount
    integer(4) :: f
    integer(4), dimension(:, :), allocatable :: connectivity
    real(4), dimension(:, :), allocatable ::  x, y, z, p
    integer :: allocateStatus
    integer :: maxCells, maxNodes
    integer :: i, j, k, index
    integer(4), dimension(3) :: nGNodes, nGCells ! Partition ghost node and ghost cell counts
    integer(4), dimension(:, :), allocatable, target :: ghostNodes, gNPartitions, gNPNodes, ghostCells
    ! Upper bound on counts of ghost nodes and cells
    integer :: maxGCells = XDIM * ZDIM + YDIM * ZDIM
    integer :: maxGNodes = 2 * (XDIM * ZDIM + YDIM * ZDIM)
    integer(4) :: effectiveNumPartitions
    character(40) :: auxDataValue

    include "tecio.f90"

    interface
        subroutine GatherGhostNodesAndCells( &
            nGNodes, ghostNodes, gNPartitions, gNPNodes, &
            nGCells, ghostCells, iDim, jDim, kDim, jMin, jMax, maxGNodes, maxGCells)
            implicit none
            integer(4), dimension(maxGNodes,3), target :: ghostNodes, gNPartitions, gNPNodes
            integer(4), dimension(maxGCells,3), target :: ghostCells
            integer(4), dimension(3) :: nGNodes, nGCells, iDim, jDim, kDim, jMin, jMax
            integer(4) :: maxGNodes, maxGCells
        end subroutine GatherGhostNodesAndCells
    end interface

    createData = 0

    ! Add aux data to solution files; for MPI, only the main output rank may do this.
    if (numSolutionFiles == 0 .or. fOffset > 0) then
        write(auxDataValue, fmt='(i10)') fOffset
#       if defined TECIOMPI
            if (commRank == 0) then
#       endif
        createData = TECZAUXSTR142("TimeStep"//char(0), trim(adjustl(auxDataValue))//char(0))
#       if defined TECIOMPI
            endif
#       endif
    endif

    if (numPartitions .gt. 0) then
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
        jMin(2) = jMin(1)
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

        effectiveNumPartitions = numPartitions
    else
        iMin(1) = 0
        iMax(1) = XDIM
        jMin(1) = 0
        jMax(1) = YDIM
        kMin(1) = 0
        kMax(1) = ZDIM

        iMin(2) = 0
        iMax(2) = 0
        jMin(2) = 0
        jMax(2) = 0
        kMin(2) = 0
        kMax(2) = 0

        iMin(3) = 0
        iMax(3) = 0
        jMin(3) = 0
        jMax(3) = 0
        kMin(3) = 0
        kMax(3) = 0

        effectiveNumPartitions = 1
    endif

    ! Local partition dimensions (of equivalent ordered zones)
    do ptn = 1, 3
        if (ptn .le. effectiveNumPartitions) then
            iDim(ptn) = iMax(ptn) - iMin(ptn)
            jDim(ptn) = jMax(ptn) - jMin(ptn)
            kDim(ptn) = kMax(ptn) - kMin(ptn)
        else
            iDim(ptn) = 0
            jDim(ptn) = 0
            kDim(ptn) = 0
        endif
    enddo

    ! Allocate memory for connectivity and variable values
    do ptn = 1, effectiveNumPartitions
        if (ptn .le. effectiveNumPartitions) then
            pNNodes(ptn) = iDim(ptn) * jDim(ptn) * kDim(ptn)
            pNCells(ptn) = (iDim(ptn) - 1) * (jDim(ptn) - 1) * (kDim(ptn) - 1)
        else
            pNNodes(ptn) = 0
            pNCells(ptn) = 0
        endif
    enddo
    
    maxCells = max(pNCells(1), pNCells(2), pNCells(3))
    maxNodes = max(pNNodes(1), pNNodes(2), pNNodes(3))
    allocate(connectivity(8 * maxCells, 3), stat = allocateStatus)
    if (allocateStatus .eq. 0) allocate(x(maxNodes, 3), stat = allocateStatus)
    if (allocateStatus .eq. 0) allocate(y(maxNodes, 3), stat = allocateStatus)
    if (allocateStatus .eq. 0) allocate(z(maxNodes, 3), stat = allocateStatus)
    if (allocateStatus .eq. 0) allocate(p(maxCells, 3), stat = allocateStatus)
    if (allocateStatus .ne. 0) stop 'Unable to allocate memory'
    
    ! Calculate variable and connectivity values for partitions.
    do ptn = 1, effectiveNumPartitions
        if (fOffset .eq. 0) then
            ! Create grid variables
            do k = 0, kDim(ptn) - 1
                do j = 0, jDim(ptn) - 1
                    do i = 0, iDim(ptn) - 1
                        index = (k * jDim(ptn) + j) * iDim(ptn) + i + 1
                        x(index, ptn) = real(i + iMin(ptn) + 1)
                        y(index, ptn) = real(j + jMin(ptn) + 1)
                        z(index, ptn) = real(k + kMin(ptn) + 1)
                    enddo
                enddo
            enddo
        else
            x(1, ptn) = 0.0
            y(1, ptn) = 0.0
            z(1, ptn) = 0.0
        endif

        ! p (cell-centered) and connectivity
        do k = 0, kDim(ptn) - 2
            do j = 0, jDim(ptn) - 2
                do i = 0, iDim(ptn) - 2
                    index = (k * (jDim(ptn) - 1) + j) * (iDim(ptn) - 1) + i
                    p(index + 1, ptn) = real((i + iMin(ptn) + 1) * (j + jMin(ptn) + 1) * (k + kMin(ptn) + 1) + solTime)

                    if (fOffset .eq. 0) then
                        connectivity(8 * index + 1, ptn) = (k * jDim(ptn) + j) * iDim(ptn) + i + 1
                        connectivity(8 * index + 2, ptn) = connectivity(8 * index + 1, ptn) + 1
                        connectivity(8 * index + 3, ptn) = connectivity(8 * index + 1, ptn) + iDim(ptn) + 1
                        connectivity(8 * index + 4, ptn) = connectivity(8 * index + 1, ptn) + iDim(ptn)
                        connectivity(8 * index + 5, ptn) = connectivity(8 * index + 1, ptn) + iDim(ptn) * jDim(ptn)
                        connectivity(8 * index + 6, ptn) = connectivity(8 * index + 2, ptn) + iDim(ptn) * jDim(ptn)
                        connectivity(8 * index + 7, ptn) = connectivity(8 * index + 3, ptn) + iDim(ptn) * jDim(ptn)
                        connectivity(8 * index + 8, ptn) = connectivity(8 * index + 4, ptn) + iDim(ptn) * jDim(ptn)
                    endif
                enddo
            enddo
        enddo
    enddo

    if (numPartitions .gt. 0) then
        allocate(ghostNodes(maxGNodes, 3), stat = allocateStatus)
        if (allocateStatus .eq. 0) allocate(gNPartitions(maxGNodes, 3), stat = allocateStatus)
        if (allocateStatus .eq. 0) allocate(gNPNodes(maxGNodes, 3), stat = allocateStatus)
        if (allocateStatus .eq. 0) allocate(ghostCells(maxGCells, 3), stat = allocateStatus)
        if (allocateStatus .ne. 0) stop 'Unable to allocate memory for ghost items'
        call GatherGhostNodesAndCells(nGNodes, ghostNodes, gNPartitions, gNPNodes, nGCells, ghostCells, &
            iDim, jDim, kDim, jMin, jMax, maxGNodes, maxGCells)
    endif

    do ptn = 1, effectiveNumPartitions
#       if defined TECIOMPI
        if (numPartitions .eq. 0 .or. partitionOwners(ptn) .eq. commRank) then
#       endif
        if (numPartitions .gt. 0) then
            if (createData .eq. 0) &
                createData = TECFEPTN142(ptn, pNNodes(ptn), pNCells(ptn), &
                    nGNodes(ptn), ghostNodes(1, ptn), gNPartitions(1, ptn), gNPNodes(1, ptn), &
                    nGCells(ptn), ghostCells(1, ptn))
        endif

        if (fOffset .eq. 0) then
            if (createData .eq. 0) then
                createData = TECDATF142(pNNodes(ptn), x(1, ptn))
            endif
            if (createData .eq. 0) then
                createData = TECDATF142(pNNodes(ptn), y(1, ptn))
            endif
            if (createData .eq. 0) then
                createData = TECDATF142(pNNodes(ptn), z(1, ptn))
            endif
        endif

        if (fOffset .ne. 0 .or. numSolutionFiles .eq. 0) then
            ! Write out the solution variable.
            if (createData .eq. 0) then
                createData = TECDATF142(pNCells(ptn), p(1, ptn))
            else
            endif
        endif

        if (fOffset .eq. 0) then
            ! Write out the connectivityCount
            connectivityCount = 8 * pNCells(ptn)
            if (createData .eq. 0) then
                createData = TECNODE142(connectivityCount, connectivity(1, ptn))
            endif
        endif

#       if defined TECIOMPI
        endif
#       endif
    enddo

    if (numPartitions .gt. 0) then
        deallocate(ghostNodes, gNPartitions, gNPNodes, ghostCells)
    endif

    deallocate(x, y, z, p)

    return
end

integer(4) function finalizeFile()
    include "tecio.f90"
    finalizeFile = TECEND142()
    return
end

subroutine appendGhostItems( &
    nGhosts, ghosts, gPartitions, gPGhosts, ptn, gIDim, gJDim, &
    gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd, oIDim, oJDim, &
    oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd, doPartition)
    
    integer(4) :: nGhosts, ptn
    integer(4), dimension(:) :: ghosts, gPartitions, gPGhosts
    integer(4) :: gIDim, gJDim, gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd
    integer(4) :: oIDim, oJDim, oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd
    logical :: doPartition

    integer :: i, j, k
    integer(4) oI, oJ, oK
    do i = gIStart, gIEnd - 1
        do j = gJStart, gJEnd - 1
            do k = gKStart, gKEnd - 1
                nGhosts = nGhosts + 1
                ghosts(nGhosts) = (k * gJDim + j) * gIDim + i + 1
                if (doPartition) then
                    gPartitions(nGhosts) = ptn
                    oI = i - gIStart + oIStart
                    oJ = j - gJStart + oJStart
                    oK = k - gKStart + oKStart
                    gPGhosts(nGhosts) = (oK * oJDim + oJ) * oIDim + oI + 1
                endif
            enddo
        enddo
    enddo
end


subroutine GatherGhostNodesAndCells( &
    nGNodes, ghostNodes, gNPartitions, gNPNodes, &
    nGCells, ghostCells, iDim, jDim, kDim, jMin, jMax, maxGNodes, maxGCells)
    implicit none
    integer(4), target :: ghostNodes(maxGNodes,3), gNPartitions(maxGNodes, 3), gNPNodes(maxGNodes, 3), ghostCells(maxGCells, 3)
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
        subroutine appendGhostItems( &
        nGhosts, ghosts, gPartitions, gPGhosts, ptn, gIDim, gJDim, &
        gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd, oIDim, oJDim, &
        oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd, doPartition)
    
        integer(4) :: nGhosts, ptn
        integer(4), dimension(:) :: ghosts, gPartitions, gPGhosts
        integer(4) :: gIDim, gJDim, gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd
        integer(4) :: oIDim, oJDim, oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd
        logical :: doPartition
        end subroutine appendGhostItems
    end interface

    NullPtr = 0

    nGNodes(1) = 0
    nGCells(1) = 0
    gnPtr => ghostNodes(1:maxGNodes, 1)
    gnpPtr => gnPartitions(1:maxGNodes, 1)
    gnpnPtr => gnpNodes(1:maxGNodes, 1)
    gcPtr => ghostCells(1:maxGCells, 1)
    ! Nodes owned by the second partition:
    call appendGhostItems( &
        nGNodes(1), gnPtr, gnpPtr, gnpnPtr, 2, &
        iDim(1), jDim(1), &                                 ! I- and J-dimensions
        iDim(1) - 1, iDim(1), 0, jDim(2) - 1, 0, kDim(1), & ! local index ranges
        iDim(2), jDim(2), &                                 ! I- and J-dimensions
        2, 3, 0, jDim(2) - 1, 0, kDim(2), .true.)           ! local index ranges
    ! Nodes owned by the third partition:
    call appendGhostItems( &
        nGNodes(1), gnPtr, gnpPtr, gnpnPtr, 3, &
        iDim(1), jDim(1), &                                       ! I- and J-dimensions
        iDim(1) - 1, iDim(1), jMin(3) + 2, jMax(3), 0, kDim(1), & ! local index ranges
        iDim(3), jDim(3), &                                       ! I- and J-dimensions
        2, 3, 0, jDim(3), 0, kDim(3), .true.)                     ! local index ranges
    ! Cells owned by the second partition:
    call appendGhostItems( &
        nGCells(1), gcptr, NULL, NULL, 2, &
        iDim(1) - 1, jDim(1) - 1, &                                 ! I- and J-dimensions
        iDim(1) - 2, iDim(1) - 1, 0, jDim(2) - 2, 0, kDim(1) - 1, & ! local index ranges
        iDim(2) - 1, jDim(2) - 1, &                                 ! I- and J-dimensions
        1, 2, 0, jDim(2) - 2, 0, kDim(2) - 1, .false.)              ! local index ranges
    ! Cells owned by the third partition:
    call appendGhostItems( &
        nGCells(1), gcptr, NULL, NULL, 3, &
        iDim(1) - 1, jDim(1) - 1, &                                           ! I- and J-dimensions
        iDim(1) - 2, iDim(1) - 1, jDim(2) - 2, jDim(1) - 1, 0, kDim(1) - 1, & ! local index ranges
        iDim(3) - 1, jDim(3) - 1, &                                           ! I- and J-dimensions
        1, 2, 0, jDim(3) - 1, 0, kDim(3), .false.)                            ! local index ranges

    nGNodes(2) = 0
    nGCells(2) = 0
    gnPtr => ghostNodes(1:maxGNodes, 2)
    gnpPtr => gnPartitions(1:maxGNodes, 2)
    gnpnPtr => gnpNodes(1:maxGNodes, 2)
    gcPtr => ghostCells(1:maxGCells, 2)
    ! Nodes owned by the first partition.
    call appendGhostItems( &
        nGNodes(2), gnPtr, gnpPtr, gnpnPtr, 1, &
        iDim(2), jDim(2), &                                       ! I- and J-dimensions
        0, 2, 0, jDim(2), 0, kDim(2), &                           ! local index ranges
        iDim(1), jDim(1), &                                       ! I- and J-dimensions
        iDim(1) - 3, iDim(1) - 1, 0, jDim(2), 0, kDim(1), .true.) ! local index ranges
    ! Nodes owned by the third partition.
    call appendGhostItems( &
        nGNodes(2), gnPtr, gnpPtr, gnpnPtr, 3, &
        iDim(2), jDim(2), &                             ! I- and J-dimensions
        2, iDim(2), jDim(2) - 1, jDim(2), 0, kDim(2), & ! local index ranges
        iDim(3), jDim(3), &                             ! I- and J-dimensions
        2, iDim(3), 2, 3, 0, kDim(3), .true.)           ! local index ranges
    ! Cells owned by the first partition.
    call appendGhostItems( &
        nGCells(2), gcptr, NULL, NULL, 1, &
        iDim(2) - 1, jDim(2) - 1, &                                        ! I- and J-dimensions
        0, 1, 0, jDim(2) - 1, 0, kDim(2) - 1, &                            ! local index ranges
        iDim(1) - 1, jDim(1) - 1, &                                        ! I- and J-dimensions
        iDim(1) - 3, iDim(1) - 2, 0, jDim(2) - 1, 0, kDim(1) - 1, .false.) ! local index ranges
    ! Cells owned by the third partition.
    call appendGhostItems( &
        nGCells(2), gcptr, NULL, NULL, 3, &
        iDim(2) - 1, jDim(2) - 1, &                                 ! I- and J-dimensions
        1, iDim(2) - 1, jDim(2) - 2, jDim(2) - 1, 0, kDim(2) - 1, & ! local index ranges
        iDim(3) - 1, jDim(3) - 1, &                                 ! I- and J-dimensions
        1, iDim(3) - 1, 1, 2, 0, kDim(3) - 1, .false.)              ! local index ranges

    nGNodes(3) = 0
    nGCells(3) = 0
    gnPtr => ghostNodes(1:maxGNodes, 3)
    gnpPtr => gnPartitions(1:maxGNodes, 3)
    gnpnPtr => gnpNodes(1:maxGNodes, 3)
    gcPtr => ghostCells(1:maxGCells, 3)
    ! Nodes owned by the first partition
    call appendGhostItems( &
        nGNodes(3), gnPtr, gnpPtr, gnpnPtr, 1, &
        iDim(3), jDim(3), &                                             ! I- and J-dimensions
        0, 2, 0, jDim(3), 0, kDim(3), &                                 ! local index ranges
        iDim(1), jDim(1), &                                             ! I- and J-dimensions
        iDim(1) - 3, iDim(1) - 1, jMin(3), jMax(3), 0, kDim(1), .true.) ! local index ranges
    ! Nodes owned by the second partition.
    call appendGhostItems( &
        nGNodes(3), gnPtr, gnpPtr, gnpnPtr, 2, &
        iDim(3), jDim(3), &                                       ! I- and J-dimensions
        2, iDim(3), 0, 2, 0, kDim(3), &                           ! local index ranges
        iDim(2), jDim(2), &                                       ! I- and J-dimensions
        2, iDim(2), jDim(2) - 3, jDim(2) - 1, 0, kDim(2), .true.) ! local index ranges
    ! Cells owned by the first partition
    call appendGhostItems( &
        nGCells(3), gcPtr, NULL, NULL, 1, &
        iDim(3) - 1, jDim(3) - 1, &                                              ! I- and J-dimensions
        0, 1, 0, jDim(3) - 1, 0, kDim(3) - 1, &                                  ! local index ranges
        iDim(1) - 1, jDim(1) - 1, &                                              ! I- and J-dimensions
        iDim(1) - 2, iDim(1) - 1, jMin(3), jMax(3) - 1, 0, kDim(1) - 1, .false.) ! local index ranges
    ! Nodes owned by the second partition.
    call appendGhostItems( &
        nGCells(3), gcPtr, NULL, NULL, 2, &
        iDim(3) - 1, jDim(3) - 1, &                                        ! I- and J-dimensions
        1, iDim(3) - 1, 0, 1, 0, kDim(3) - 1, &                            ! local index ranges
        iDim(2) - 1, jDim(2) - 1, &                                        ! I- and J-dimensions
        1, iDim(2) - 1, jDim(2) - 2, jDim(2) - 1, 0, kDim(2) - 1, .false.) ! local index ranges
end
