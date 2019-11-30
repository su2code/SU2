/*
 * Example C++ program to write a partitioned
 * binary datafile for tecplot. This example
 * does the following:
 *
 *   1.  Open a datafile called "tetrapartitioned.szplt"
 *   2.  Assign values for x, y, z and p.
 *   3.  Write out a tetrahedral zone in 3 partitions.
 *   4.  Close the datafile.
 *
 * If TECIOMPI is #defined, this program may be executed with mpiexec with
 * up to 3 MPI ranks (processes). In this case, it must be linked
 * with the MPI version of TECIO.
 */

#if defined TECIOMPI
#include <mpi.h>
#endif

#include "TECIO.h"
#include <cstdlib>
#include <cstdio>

#ifndef NULL
#define NULL 0
#endif

// The zone will appear the same as an ordered zone with these dimensions:
#define XDIM 10
#define YDIM 9
#define ZDIM 8

enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };

void GatherGhostNodesAndCells(
    INTEGER4* nGNodes, INTEGER4** ghostNodes, INTEGER4** gNPartitions, INTEGER4** gNPNodes,
    INTEGER4* nGCells, INTEGER4** ghostCells, int* iDim, int* jDim, int* kDim, int* jMin, int* jMax);

int main(int argc, char** argv)
{
    /*
     * Open the file and write the tecplot datafile
     * header information
     */
    INTEGER4 fileFormat = 1; // 1 == SZPLT (cannot output partitioned zones to PLT)
    INTEGER4 fileType   = FULL;
    INTEGER4 debug = 1;
    INTEGER4 vIsDouble = 0;
    int returnValue = 0;
#if defined TECIOMPI
    int commSize;
    int commRank;
    int mainRank = 0;
    INTEGER4 numPartitions = 3;
    INTEGER4 partitionOwners[3];
    MPI_Init(&argc, &argv);
    MPI_Comm mpiComm = MPI_COMM_WORLD;
    MPI_Comm_size(mpiComm, &commSize);
    MPI_Comm_rank(mpiComm, &commRank);
#endif
    returnValue = TECINI142((char*)"SIMPLE DATASET",
                  (char*)"x y z p",
                  (char*)"tetrapartitioned.szplt",
                  (char*)".",
                  &fileFormat,
                  &fileType,
                  &debug,
                  &vIsDouble);
#if defined TECIOMPI
    if (returnValue == 0)
        returnValue = TECMPIINIT142(&mpiComm, &mainRank);
#endif

    /*
     * Zone
     */
    INTEGER4 zoneType  = 4;      // Brick
    INTEGER4 nNodes    = XDIM * YDIM * ZDIM; // Overall zone dimensions
    INTEGER4 nCells    = (XDIM - 1) * (YDIM - 1) * (ZDIM - 1) * 5;
    INTEGER4 nFaces    = 6;         // Not used
    INTEGER4 iCellMax  = 0;
    INTEGER4 jCellMax  = 0;
    INTEGER4 kCellMax  = 0;
    double solTime     = 360.0;
    INTEGER4 strandID  = 0;      // StaticZone
    INTEGER4 parentZn  = 0;      // No Parent
    INTEGER4 isBlock   = 1;      // Block
    INTEGER4 nFConns   = 0;
    INTEGER4 fNMode    = 0;
    INTEGER4 valueLocations[4] = {1, 1, 1, 0}; // 1 = Nodal, 0 = Cell-Centered
    INTEGER4 shrConn   = 0;

    if (returnValue == 0)
        returnValue = TECZNE142((char*)"Partitioned Zone",
                      &zoneType,
                      &nNodes,
                      &nCells,
                      &nFaces,
                      &iCellMax,
                      &jCellMax,
                      &kCellMax,
                      &solTime,
                      &strandID,
                      &parentZn,
                      &isBlock,
                      &nFConns,
                      &fNMode,
                      0,              // TotalNumFaceNodes
                      0,              // NumConnectedBoundaryFaces
                      0,              // TotalNumBoundaryConnections
                      NULL,           // PassiveVarList
                      valueLocations,
                      NULL,           // SharVarFromZone
                      &shrConn);

    // Divide the zone into 3 partitions, identified by the index ranges
    // of an equivalent unpartitioned ordered zone.
    int iMin[3], iMax[3], jMin[3], jMax[3], kMin[3], kMax[3];

    // Partition 1 node range, which will include one layer of ghost cells on the IMAX boundary:
    iMin[0] = 0; // We use zero-based indices here because it simplifies index arithmetic in C.
    iMax[0] = XDIM / 2 + 2;
    jMin[0] = 0;
    jMax[0] = YDIM;
    kMin[0] = 0;
    kMax[0] = ZDIM;

    // Partition 2; ghost cells on IMIN and JMAX boundaries:
    iMin[1] = iMax[0] - 3;
    iMax[1] = XDIM;
    jMin[1] = jMin[0];
    jMax[1] = YDIM / 2 + 2;
    kMin[1] = kMin[0];
    kMax[1] = kMax[0];

    // Partition 3; ghost cells on IMIN and JMIN boundaries:
    iMin[2] = iMin[1];
    iMax[2] = iMax[1];
    jMin[2] = jMax[1] - 3;
    jMax[2] = YDIM;
    kMin[2] = kMin[1];
    kMax[2] = kMax[1];

    // Local partition dimensions (of equivalent ordered zones)
    int iDim[3], jDim[3], kDim[3];
    for(int ii = 0; ii < 3; ++ii)
    {
        iDim[ii] = iMax[ii] - iMin[ii];
        jDim[ii] = jMax[ii] - jMin[ii];
        kDim[ii] = kMax[ii] - kMin[ii];
    }

    // Calculate variable and connectivity values for partitions.
    int *connectivity[3];
    float *x[3], *y[3], *z[3], *p[3];
    INTEGER4 pNNodes[3], pNCells[3]; // Partition node and cell counts, including ghost items
    for (INTEGER4 ptn = 0; ptn < 3; ++ptn)
    {
        pNNodes[ptn] = iDim[ptn] * jDim[ptn] * kDim[ptn];
        x[ptn] = (float*)malloc(pNNodes[ptn] * sizeof(float));
        y[ptn] = (float*)malloc(pNNodes[ptn] * sizeof(float));
        z[ptn] = (float*)malloc(pNNodes[ptn] * sizeof(float));
        // Nodes
        for (int k = 0; k < kDim[ptn]; ++k)
            for (int j = 0; j < jDim[ptn]; ++j)
                for (int i = 0; i < iDim[ptn]; ++i)
                {
                    int index = (k * jDim[ptn] + j) * iDim[ptn] + i;
                    x[ptn][index] = (float)(i + iMin[ptn] + 1);
                    y[ptn][index] = (float)(j + jMin[ptn] + 1);
                    z[ptn][index] = (float)(k + kMin[ptn] + 1);
                }
        // p (cell-centered) and connectivity
        pNCells[ptn] = (iDim[ptn] - 1) * (jDim[ptn] - 1) * (kDim[ptn] - 1) * 5;
        p[ptn] = (float*)malloc(pNCells[ptn] * sizeof(float));
        int connectivityCount = 4 * pNCells[ptn];
        connectivity[ptn] = (INTEGER4*)malloc(connectivityCount * sizeof(INTEGER4));
        for (int k = 0; k < kDim[ptn] - 1; ++k)
            for (int j = 0; j < jDim[ptn] - 1; ++j)
                for (int i = 0; i < iDim[ptn] - 1; ++i)
                {
                    /*
                     * Start with a brick and divide it into 5 tets.
                     * Need to mirror the subdivision for neighboring
                     * bricks to avoid internal faces.
                     */
                    INTEGER4 brickNodes[8];
                    int tetraCorners[2][5][4] = {
                        { /* "Even" bricks */
                            {0,1,2,5},
                            {0,2,3,7},
                            {4,7,5,0},
                            {5,7,6,2},
                            {0,2,7,5}
                        },
                        { /* "Odd" bricks */
                            {0,1,3,4},
                            {1,2,3,6},
                            {4,6,5,1},
                            {4,7,6,3},
                            {1,3,4,6}
                        }
                    };
                    int evenOdd = (i + iMin[ptn] + j + jMin[ptn] + k + kMin[ptn]) % 2;
                    int whichTet;
                    int whichCorner;

                    brickNodes[0] = (k * jDim[ptn] + j) * iDim[ptn] + i + 1;
                    brickNodes[1] = brickNodes[0] + 1;
                    brickNodes[2] = brickNodes[0] + iDim[ptn] + 1;
                    brickNodes[3] = brickNodes[0] + iDim[ptn];
                    brickNodes[4] = brickNodes[0] + iDim[ptn] * jDim[ptn];
                    brickNodes[5] = brickNodes[1] + iDim[ptn] * jDim[ptn];
                    brickNodes[6] = brickNodes[2] + iDim[ptn] * jDim[ptn];
                    brickNodes[7] = brickNodes[3] + iDim[ptn] * jDim[ptn];

                    for(whichTet = 0; whichTet < 5; ++whichTet)
                    {
                        int index = ((k * (jDim[ptn] - 1) + j) * (iDim[ptn] - 1) + i) * 5 + whichTet;
                        p[ptn][index] = (float)((i + iMin[ptn] + 1) * (j + jMin[ptn] + 1) * (k + kMin[ptn] + 1));
                        for(whichCorner = 0; whichCorner < 4; ++whichCorner)
                        {
                            connectivity[ptn][index * 4 + whichCorner] = brickNodes[tetraCorners[evenOdd][whichTet][whichCorner]];
                        }
                    }
                }
    }

    INTEGER4 nGNodes[3], nGCells[3]; // Partition ghost node and ghost cell counts
    INTEGER4 *ghostNodes[3], *gNPartitions[3], *gNPNodes[3], *ghostCells[3];
    GatherGhostNodesAndCells(nGNodes, ghostNodes, gNPartitions, gNPNodes, nGCells, ghostCells, iDim, jDim, kDim, jMin, jMax);

    // Output partitions
#if defined TECIOMPI
    for(int ptn = 0; ptn < 3; ++ptn)
        partitionOwners[ptn] = (ptn % commSize);
    if (returnValue == 0)
        returnValue = TECZNEMAP142(&numPartitions, &partitionOwners[0]);
#endif

    INTEGER4 dIsDouble = 0;
    for(INTEGER4 ptn = 1; ptn <= 3; ++ptn)
    {
#if defined TECIOMPI
        if (partitionOwners[ptn - 1] == commRank)
        {
#endif
        if (returnValue == 0)
            returnValue = TECFEPTN142(&ptn, &pNNodes[ptn - 1], &pNCells[ptn - 1],
                &nGNodes[ptn - 1], ghostNodes[ptn - 1], gNPartitions[ptn - 1], gNPNodes[ptn - 1],
                &nGCells[ptn - 1], ghostCells[ptn - 1]);
        if (ghostNodes[ptn - 1])
            free(ghostNodes[ptn - 1]);
        if (gNPartitions[ptn - 1])
            free(gNPartitions[ptn - 1]);
        if (gNPNodes[ptn - 1])
            free(gNPNodes[ptn - 1]);
        if (ghostCells[ptn - 1])
            free(ghostCells[ptn - 1]);
        /*
         * Write out the field data.
         */
        if (returnValue == 0)
            returnValue = TECDAT142(&pNNodes[ptn - 1], x[ptn - 1], &dIsDouble);
        if (returnValue == 0)
            returnValue = TECDAT142(&pNNodes[ptn - 1], y[ptn - 1], &dIsDouble);
        if (returnValue == 0)
            returnValue = TECDAT142(&pNNodes[ptn - 1], z[ptn - 1], &dIsDouble);
        if (returnValue == 0)
            returnValue = TECDAT142(&pNCells[ptn - 1], p[ptn - 1], &dIsDouble);

        INTEGER4 connectivityCount = 4 * pNCells[ptn - 1];
        if (returnValue == 0)
            returnValue = TECNODE142(&connectivityCount, connectivity[ptn - 1]);

        free(x[ptn - 1]);
        free(y[ptn - 1]);
        free(z[ptn - 1]);
        free(p[ptn - 1]);
        free(connectivity[ptn - 1]);
#if defined TECIOMPI
        }
#endif
    }

    if (returnValue == 0)
        returnValue = TECEND142();

#if defined TECIOMPI
    MPI_Finalize();
#endif

    return returnValue;
}

void appendGhostNodes(
    INTEGER4* nGhosts, INTEGER4** ghosts, INTEGER4** gPartitions, INTEGER4** gPGhosts, INTEGER4 ptn,
    INTEGER4 gIDim, INTEGER4 gJDim,
    INTEGER4 gIStart, INTEGER4 gIEnd, INTEGER4 gJStart, INTEGER4 gJEnd, INTEGER4 gKStart, INTEGER4 gKEnd,
    INTEGER4 oIDim, INTEGER4 oJDim,
    INTEGER4 oIStart, INTEGER4 oIEnd, INTEGER4 oJStart, INTEGER4 oJEnd, INTEGER4 oKStart, INTEGER4 oKEnd)
{
    INTEGER4 localNumGhosts = (gIEnd - gIStart) * (gJEnd - gJStart) * (gKEnd - gKStart);
    INTEGER4 gOffset = *nGhosts;
    *nGhosts += localNumGhosts;
    *ghosts = (INTEGER4*)realloc(*ghosts, *nGhosts * sizeof(INTEGER4));
    *gPartitions = (INTEGER4*)realloc(*gPartitions, *nGhosts * sizeof(INTEGER4));
    *gPGhosts = (INTEGER4*)realloc(*gPGhosts, *nGhosts * sizeof(INTEGER4));

    for(INTEGER4 i = gIStart; i < gIEnd; ++i)
        for(INTEGER4 j = gJStart; j < gJEnd; ++j)
            for(INTEGER4 k = gKStart; k < gKEnd; ++k)
            {
                (*ghosts)[gOffset] = (k * gJDim + j) * gIDim + i + 1;
                (*gPartitions)[gOffset] = ptn;
                INTEGER4 oI = i - gIStart + oIStart;
                INTEGER4 oJ = j - gJStart + oJStart;
                INTEGER4 oK = k - gKStart + oKStart;
                (*gPGhosts)[gOffset] = (oK * oJDim + oJ) * oIDim + oI + 1;
                ++gOffset;
            }
}

void appendGhostCells(
    INTEGER4* nGhosts, INTEGER4** ghosts, INTEGER4 gIDim, INTEGER4 gJDim,
    INTEGER4 gIStart, INTEGER4 gIEnd, INTEGER4 gJStart, INTEGER4 gJEnd, INTEGER4 gKStart, INTEGER4 gKEnd)
{
    INTEGER4 localNumGhosts = (gIEnd - gIStart) * (gJEnd - gJStart) * (gKEnd - gKStart) * 5;
    INTEGER4 gOffset = *nGhosts;
    *nGhosts += localNumGhosts;
    *ghosts = (INTEGER4*)realloc(*ghosts, *nGhosts * sizeof(INTEGER4));

    for(INTEGER4 i = gIStart; i < gIEnd; ++i)
        for(INTEGER4 j = gJStart; j < gJEnd; ++j)
            for(INTEGER4 k = gKStart; k < gKEnd; ++k)
                for(INTEGER4 whichTet = 0; whichTet < 5; ++whichTet)
                {
                    (*ghosts)[gOffset] = ((k * gJDim + j) * gIDim + i) * 5 + whichTet + 1;
                    ++gOffset;
                }
}


void GatherGhostNodesAndCells(
    INTEGER4* nGNodes, INTEGER4** ghostNodes, INTEGER4** gNPartitions, INTEGER4** gNPNodes,
    INTEGER4* nGCells, INTEGER4** ghostCells, int* iDim, int* jDim, int* kDim, int* jMin, int* jMax)
{
    /*
     * Assemble lists of ghost nodes and cells--nodes and cells near partition boundaries that
     * coincide with those "owned" by neighboring partitions. For each set of coincident nodes
     * or cells, exactly one partition must own the node or cell, and all other involved partitions
     * must report it as a ghost node or cell.
     *
     * Arbitrarily, we say that the first partition owns any nodes that do not overlay the interior
     * of neighboring partitions. That is, it owns any nodes that its "real" (non-ghost) cells use.
     * So only a single layer of nodes and cells--on its IMax boundary--are ghosts, and are owned
     * by the second and third partitions. We use the same logic to assign ownership for nodes
     * shared by partitions 2 and 3.
     */
    nGNodes[0] = 0;
    ghostNodes[0] = NULL;
    gNPartitions[0] = NULL;
    gNPNodes[0] = NULL;
    nGCells[0] = 0;
    ghostCells[0] = NULL;
    // Nodes owned by the second partition:
    appendGhostNodes(
        &nGNodes[0], &ghostNodes[0], &gNPartitions[0], &gNPNodes[0], 2,
        iDim[0], jDim[0],                                 // I- and J-dimensions
        iDim[0] - 1, iDim[0], 0, jDim[1] - 1, 0, kDim[0], // local index ranges
        iDim[1], jDim[1],                                 // I- and J-dimensions
        2, 3, 0, jDim[1] - 1, 0, kDim[1]);                // local index ranges
    // Nodes owned by the third partition:
    appendGhostNodes(
        &nGNodes[0], &ghostNodes[0], &gNPartitions[0], &gNPNodes[0], 3,
        iDim[0], jDim[0],                                       // I- and J-dimensions
        iDim[0] - 1, iDim[0], jMin[2] + 2, jMax[2], 0, kDim[0], // local index ranges
        iDim[2], jDim[2],                                       // I- and J-dimensions
        2, 3, 0, jDim[2], 0, kDim[2]);                          // local index ranges
    // Cells owned by the second partition:
    appendGhostCells(
        &nGCells[0], &ghostCells[0],
        iDim[0] - 1, jDim[0] - 1,                                  // I- and J-dimensions
        iDim[0] - 2, iDim[0] - 1, 0, jDim[1] - 2, 0, kDim[0] - 1); // local index ranges
    // Cells owned by the third partition:
    appendGhostCells(
        &nGCells[0], &ghostCells[0],
        iDim[0] - 1, jDim[0] - 1,                                            // I- and J-dimensions
        iDim[0] - 2, iDim[0] - 1, jDim[1] - 2, jDim[0] - 1, 0, kDim[0] - 1); // local index ranges

    nGNodes[1] = 0;
    ghostNodes[1] = NULL;
    gNPartitions[1] = NULL;
    gNPNodes[1] = NULL;
    nGCells[1] = 0;
    ghostCells[1] = NULL;
    // Nodes owned by the first partition.
    appendGhostNodes(
        &nGNodes[1], &ghostNodes[1], &gNPartitions[1], &gNPNodes[1], 1,
        iDim[1], jDim[1],                                  // I- and J-dimensions
        0, 2, 0, jDim[1], 0, kDim[1],                      // local index ranges
        iDim[0], jDim[0],                                  // I- and J-dimensions
        iDim[0] - 3, iDim[0] - 1, 0, jDim[1], 0, kDim[0]); // local index ranges
    // Nodes owned by the third partition.
    appendGhostNodes(
        &nGNodes[1], &ghostNodes[1], &gNPartitions[1], &gNPNodes[1], 3,
        iDim[1], jDim[1],                             // I- and J-dimensions
        2, iDim[1], jDim[1] - 1, jDim[1], 0, kDim[1], // local index ranges
        iDim[2], jDim[2],                             // I- and J-dimensions
        2, iDim[2], 2, 3, 0, kDim[2]);                // local index ranges
    // Cells owned by the first partition.
    appendGhostCells(
        &nGCells[1], &ghostCells[1],
        iDim[1] - 1, jDim[1] - 1,              // I- and J-dimensions
        0, 1, 0, jDim[1] - 1, 0, kDim[1] - 1); // local index ranges
    // Cells owned by the third partition.
    appendGhostCells(
        &nGCells[1], &ghostCells[1],
        iDim[1] - 1, jDim[1] - 1,                                  // I- and J-dimensions
        1, iDim[1] - 1, jDim[1] - 2, jDim[1] - 1, 0, kDim[1] - 1); // local index ranges

    nGNodes[2] = 0;
    ghostNodes[2] = NULL;
    gNPartitions[2] = NULL;
    gNPNodes[2] = NULL;
    nGCells[2] = 0;
    ghostCells[2] = NULL;
    // Nodes owned by the first partition
    appendGhostNodes(
        &nGNodes[2], &ghostNodes[2], &gNPartitions[2], &gNPNodes[2], 1,
        iDim[2], jDim[2],                                        // I- and J-dimensions
        0, 2, 0, jDim[2], 0, kDim[2],                            // local index ranges
        iDim[0], jDim[0],                                        // I- and J-dimensions
        iDim[0] - 3, iDim[0] - 1, jMin[2], jMax[2], 0, kDim[0]); // local index ranges
    // Nodes owned by the second partition.
    appendGhostNodes(
        &nGNodes[2], &ghostNodes[2], &gNPartitions[2], &gNPNodes[2], 2,
        iDim[2], jDim[2],                                  // I- and J-dimensions
        2, iDim[2], 0, 2, 0, kDim[2],                      // local index ranges
        iDim[1], jDim[1],                                  // I- and J-dimensions
        2, iDim[1], jDim[1] - 3, jDim[1] - 1, 0, kDim[1]); // local index ranges
    // Cells owned by the first partition
    appendGhostCells(
        &nGCells[2], &ghostCells[2],
        iDim[2] - 1, jDim[2] - 1,              // I- and J-dimensions
        0, 1, 0, jDim[2] - 1, 0, kDim[2] - 1); // local index ranges
    // Nodes owned by the second partition.
    appendGhostCells(
        &nGCells[2], &ghostCells[2],
        iDim[2] - 1, jDim[2] - 1,              // I- and J-dimensions
        1, iDim[2] - 1, 0, 1, 0, kDim[2] - 1); // local index ranges
}
