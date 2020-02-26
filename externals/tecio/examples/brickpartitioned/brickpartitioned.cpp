/*
 * Example C++ program to write a partitioned binary grid and solution
 * datafiles for tecplot. This example does the following:
 *
 *   1.  Open a datafile called "brickpartitioned.szplt"
 *       and "brickpartitioned_solution.szplt"
 *   2.  Assign values for x, y, z to the grid and p to the solution.
 *   3.  Write out a hexahedral (brick) zone in 3 partitions.
 *   4.  Close the datafiles
 *
 * If TECIOMPI is #defined, this program may be executed with mpiexec with
 * up to 3 MPI ranks (processes). In this case, it must be linked
 * with the MPI version of TECIO.
 */

#if defined TECIOMPI
#include <mpi.h>
#endif

#include "TECIO.h"
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#if defined _WIN32 && !defined snprintf
#define snprintf _snprintf
#endif

// internal testing
// RUNFLAGS:none
// RUNFLAGS:--num-solution-files 5

#ifndef NULL
#define NULL 0
#endif

/* The zone will appear the same as an ordered zone with these dimensions: */
#define XDIM 10
#define YDIM 9
#define ZDIM 8

enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };

/**
 * Open the file and write the tecplot datafile header information.
 */
static INTEGER4 initializeFile(
    #if defined TECIOMPI
    MPI_Comm mpiComm,
    #endif
    int      fOffset,
    int      numSolutionFiles,
    int      numPartitions,
    double   solTime)
{
    INTEGER4 returnValue = 0;

    INTEGER4 fileFormat = 1; /* 1 == SZPLT (cannot output partitioned zones to PLT) */
    INTEGER4 debug = 1;
    INTEGER4 vIsDouble = 0;

    if (returnValue == 0)
    {
        INTEGER4 dataFileType;
        size_t const fileNameSize = 1024;
        char fileName[fileNameSize];
        char const* baseFileName = numPartitions > 0 ? "brickpartitioned" : "brick";
        char const* varNames;
        if (numSolutionFiles == 0)
        {
            dataFileType = FULL;
            snprintf(fileName, fileNameSize, "%s.szplt", baseFileName);
            varNames = "x y z p";
        }
        else if (fOffset == 0)
        {
            dataFileType = GRID;
            snprintf(fileName, fileNameSize, "%s_grid.szplt", baseFileName);
            varNames = "x y z";
        }
        else
        {
            dataFileType = SOLUTION;
            if (numSolutionFiles == 1)
                snprintf(fileName, fileNameSize, "%s_solution.szplt", baseFileName);
            else
                snprintf(fileName, fileNameSize, "%s_solution_%04d.szplt", baseFileName, static_cast<int>(std::ceil(solTime)));
            varNames = "p";
        }

        returnValue = TECINI142((char*)"SIMPLE DATASET",
                                (char*)varNames,
                                (char*)fileName,
                                (char*)".",
                                &fileFormat,
                                &dataFileType,
                                &debug,
                                &vIsDouble);
    }

    #if defined TECIOMPI
    {
        if (returnValue == 0)
        {
            INTEGER4 mainRank = 0;
            returnValue = TECMPIINIT142(&mpiComm, &mainRank);
        }
    }
    #endif

    return returnValue;
}

/**
 */
static INTEGER4 createZones(
    int             fOffset,
    INTEGER4        numPartitions,
    INTEGER4 const* partitionOwners,
    double          solTime)
{
    INTEGER4 returnValue = 0;

    INTEGER4 zoneType  = 5;      /* Brick */
    INTEGER4 nNodes    = XDIM * YDIM * ZDIM; /* Overall zone dimensions */
    INTEGER4 nCells    = (XDIM - 1) * (YDIM - 1) * (ZDIM - 1);
    INTEGER4 nFaces    = 6;         /* Not used */
    INTEGER4 iCellMax  = 0;
    INTEGER4 jCellMax  = 0;
    INTEGER4 kCellMax  = 0;
    INTEGER4 strandID  = 1;
    INTEGER4 parentZn  = 0;      /* No Parent */
    INTEGER4 isBlock   = 1;      /* Block */
    INTEGER4 nFConns   = 0;
    INTEGER4 fNMode    = 0;
    INTEGER4 valueLocations[4] = {1, 1, 1, 0}; /* 1 = Nodal, 0 = Cell-Centered */
    INTEGER4 shrConn   = 0;

    size_t const zoneTitleSize = 1024;
    char zoneTitle[zoneTitleSize];
    char const* baseTitle = numPartitions > 0 ? "Partitioned Zone" : "Zone";
    snprintf(zoneTitle, zoneTitleSize, "%s Time=%d", baseTitle, static_cast<int>(std::ceil(solTime)));
    int const baseVarOffset = fOffset == 0 ? 0 : 3;

    if (returnValue == 0)
        returnValue = TECZNE142(zoneTitle,
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
                                0,              /* TotalNumFaceNodes */
                                0,              /* NumConnectedBoundaryFaces */
                                0,              /* TotalNumBoundaryConnections */
                                NULL,           /* PassiveVarList */
                                &valueLocations[baseVarOffset],
                                NULL,           /* SharVarFromZone */
                                &shrConn);

    if (returnValue == 0)
    {
        #if defined TECIOMPI
        {
            /* Output partitions */
            if (numPartitions > 0 && returnValue == 0)
                returnValue = TECZNEMAP142(&numPartitions, partitionOwners);
        }
        #endif
    }

    return returnValue;
}

/**
 */
static void appendGhostItems(
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
    if (gPartitions)
        *gPartitions = (INTEGER4*)realloc(*gPartitions, *nGhosts * sizeof(INTEGER4));
    if (gPGhosts)
        *gPGhosts = (INTEGER4*)realloc(*gPGhosts, *nGhosts * sizeof(INTEGER4));

    for(INTEGER4 i = gIStart; i < gIEnd; ++i)
    for(INTEGER4 j = gJStart; j < gJEnd; ++j)
    for(INTEGER4 k = gKStart; k < gKEnd; ++k)
    {
        (*ghosts)[gOffset] = (k * gJDim + j) * gIDim + i + 1;
        if (gPartitions)
            (*gPartitions)[gOffset] = ptn;
        if (gPGhosts)
        {
            INTEGER4 oI = i - gIStart + oIStart;
            INTEGER4 oJ = j - gJStart + oJStart;
            INTEGER4 oK = k - gKStart + oKStart;
            (*gPGhosts)[gOffset] = (oK * oJDim + oJ) * oIDim + oI + 1;
        }
        ++gOffset;
    }
}


/**
 */
static void GatherGhostNodesAndCells(
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
    /* Nodes owned by the second partition: */
    appendGhostItems(
        &nGNodes[0], &ghostNodes[0], &gNPartitions[0], &gNPNodes[0], 2,
        iDim[0], jDim[0],                                 /* I- and J-dimensions */
        iDim[0] - 1, iDim[0], 0, jDim[1] - 1, 0, kDim[0], /* local index ranges */
        iDim[1], jDim[1],                                 /* I- and J-dimensions */
        2, 3, 0, jDim[1] - 1, 0, kDim[1]);                /* local index ranges */
    /* Nodes owned by the third partition: */
    appendGhostItems(
        &nGNodes[0], &ghostNodes[0], &gNPartitions[0], &gNPNodes[0], 3,
        iDim[0], jDim[0],                                       /* I- and J-dimensions */
        iDim[0] - 1, iDim[0], jMin[2] + 2, jMax[2], 0, kDim[0], /* local index ranges */
        iDim[2], jDim[2],                                       /* I- and J-dimensions */
        2, 3, 0, jDim[2], 0, kDim[2]);                          /* local index ranges */
    /* Cells owned by the second partition: */
    appendGhostItems(
        &nGCells[0], &ghostCells[0], NULL, NULL, 2,
        iDim[0] - 1, jDim[0] - 1,                                 /* I- and J-dimensions */
        iDim[0] - 2, iDim[0] - 1, 0, jDim[1] - 2, 0, kDim[0] - 1, /* local index ranges */
        iDim[1] - 1, jDim[1] - 1,                                 /* I- and J-dimensions */
        1, 2, 0, jDim[1] - 2, 0, kDim[1] - 1);                    /* local index ranges */
    /* Cells owned by the third partition: */
    appendGhostItems(
        &nGCells[0], &ghostCells[0], NULL, NULL, 3,
        iDim[0] - 1, jDim[0] - 1,                                           /* I- and J-dimensions */
        iDim[0] - 2, iDim[0] - 1, jDim[1] - 2, jDim[0] - 1, 0, kDim[0] - 1, /* local index ranges */
        iDim[2] - 1, jDim[2] - 1,                                           /* I- and J-dimensions */
        1, 2, 0, jDim[2] - 1, 0, kDim[2]);                                  /* local index ranges */

    nGNodes[1] = 0;
    ghostNodes[1] = NULL;
    gNPartitions[1] = NULL;
    gNPNodes[1] = NULL;
    nGCells[1] = 0;
    ghostCells[1] = NULL;
    /* Nodes owned by the first partition. */
    appendGhostItems(
        &nGNodes[1], &ghostNodes[1], &gNPartitions[1], &gNPNodes[1], 1,
        iDim[1], jDim[1],                                  /* I- and J-dimensions */
        0, 2, 0, jDim[1], 0, kDim[1],                      /* local index ranges */
        iDim[0], jDim[0],                                  /* I- and J-dimensions */
        iDim[0] - 3, iDim[0] - 1, 0, jDim[1], 0, kDim[0]); /* local index ranges */
    /* Nodes owned by the third partition. */
    appendGhostItems(
        &nGNodes[1], &ghostNodes[1], &gNPartitions[1], &gNPNodes[1], 3,
        iDim[1], jDim[1],                             /* I- and J-dimensions */
        2, iDim[1], jDim[1] - 1, jDim[1], 0, kDim[1], /* local index ranges */
        iDim[2], jDim[2],                             /* I- and J-dimensions */
        2, iDim[2], 2, 3, 0, kDim[2]);                /* local index ranges */
    /* Cells owned by the first partition. */
    appendGhostItems(
        &nGCells[1], &ghostCells[1], NULL, NULL, 1,
        iDim[1] - 1, jDim[1] - 1,                                  /* I- and J-dimensions */
        0, 1, 0, jDim[1] - 1, 0, kDim[1] - 1,                      /* local index ranges */
        iDim[0] - 1, jDim[0] - 1,                                  /* I- and J-dimensions */
        iDim[0] - 3, iDim[0] - 2, 0, jDim[1] - 1, 0, kDim[0] - 1); /* local index ranges */
    /* Cells owned by the third partition. */
    appendGhostItems(
        &nGCells[1], &ghostCells[1], NULL, NULL, 3,
        iDim[1] - 1, jDim[1] - 1,                                 /* I- and J-dimensions */
        1, iDim[1] - 1, jDim[1] - 2, jDim[1] - 1, 0, kDim[1] - 1, /* local index ranges */
        iDim[2] - 1, jDim[2] - 1,                                 /* I- and J-dimensions */
        1, iDim[2] - 1, 1, 2, 0, kDim[2] - 1);                    /* local index ranges */

    nGNodes[2] = 0;
    ghostNodes[2] = NULL;
    gNPartitions[2] = NULL;
    gNPNodes[2] = NULL;
    nGCells[2] = 0;
    ghostCells[2] = NULL;
    /* Nodes owned by the first partition */
    appendGhostItems(
        &nGNodes[2], &ghostNodes[2], &gNPartitions[2], &gNPNodes[2], 1,
        iDim[2], jDim[2],                                        /* I- and J-dimensions */
        0, 2, 0, jDim[2], 0, kDim[2],                            /* local index ranges */
        iDim[0], jDim[0],                                        /* I- and J-dimensions */
        iDim[0] - 3, iDim[0] - 1, jMin[2], jMax[2], 0, kDim[0]); /* local index ranges */
    /* Nodes owned by the second partition. */
    appendGhostItems(
        &nGNodes[2], &ghostNodes[2], &gNPartitions[2], &gNPNodes[2], 2,
        iDim[2], jDim[2],                                  /* I- and J-dimensions */
        2, iDim[2], 0, 2, 0, kDim[2],                      /* local index ranges */
        iDim[1], jDim[1],                                  /* I- and J-dimensions */
        2, iDim[1], jDim[1] - 3, jDim[1] - 1, 0, kDim[1]); /* local index ranges */
    /* Cells owned by the first partition */
    appendGhostItems(
        &nGCells[2], &ghostCells[2], NULL, NULL, 1,
        iDim[2] - 1, jDim[2] - 1,                                        /* I- and J-dimensions */
        0, 1, 0, jDim[2] - 1, 0, kDim[2] - 1,                            /* local index ranges */
        iDim[0] - 1, jDim[0] - 1,                                        /* I- and J-dimensions */
        iDim[0] - 2, iDim[0] - 1, jMin[2], jMax[2] - 1, 0, kDim[0] - 1); /* local index ranges */
    /* Nodes owned by the second partition. */
    appendGhostItems(
        &nGCells[2], &ghostCells[2], NULL, NULL, 2,
        iDim[2] - 1, jDim[2] - 1,                                  /* I- and J-dimensions */
        1, iDim[2] - 1, 0, 1, 0, kDim[2] - 1,                      /* local index ranges */
        iDim[1] - 1, jDim[1] - 1,                                  /* I- and J-dimensions */
        1, iDim[1] - 1, jDim[1] - 2, jDim[1] - 1, 0, kDim[1] - 1); /* local index ranges */
}

/**
 */
static INTEGER4 createData(
    #if defined TECIOMPI
    int             commRank,
    #endif
    int             fOffset,
    int             numSolutionFiles,
    INTEGER4        numPartitions,
    INTEGER4 const* partitionOwners,
    double          solTime)
{
    INTEGER4 returnValue = 0;

    // Add aux data to solution files; for MPI, only the main output rank may do this.
    if (numSolutionFiles == 0 || fOffset > 0)
    {
        char auxDataValue[40];
        snprintf(auxDataValue, sizeof(auxDataValue), "%d", fOffset);
        #if defined TECIOMPI
            if (commRank == 0)
        #endif
        returnValue = TECZAUXSTR142("TimeStep", auxDataValue);
    }

    /*
     * If the file isn't partitioned, we still allocate out array resources as if there was a single
     * partition as it makes the code between partitioned and non-partitioned the same.
     */
    INTEGER4 const effectiveNumPartitions = numPartitions > 0 ? numPartitions : 1;

    /*
     * Divide the zone into number of partitions, identified by the index ranges
     * of an equivalent unpartitioned ordered zone.
     */
    int* iMin = (int*)malloc(effectiveNumPartitions * sizeof(int));
    int* iMax = (int*)malloc(effectiveNumPartitions * sizeof(int));
    int* jMin = (int*)malloc(effectiveNumPartitions * sizeof(int));
    int* jMax = (int*)malloc(effectiveNumPartitions * sizeof(int));
    int* kMin = (int*)malloc(effectiveNumPartitions * sizeof(int));
    int* kMax = (int*)malloc(effectiveNumPartitions * sizeof(int));

    assert(numPartitions == 0 || numPartitions == 3); /* ...because this example only handles on or the other */
    if (numPartitions > 0)
    {
        /* Partition 1 node range, which will include one layer of ghost cells on the IMAX boundary: */
        iMin[0] = 0; /* We use zero-based indices here because it simplifies index arithmetic in C. */
        iMax[0] = XDIM / 2 + 2;
        jMin[0] = 0;
        jMax[0] = YDIM;
        kMin[0] = 0;
        kMax[0] = ZDIM;

        /* Partition 2; ghost cells on IMIN and JMAX boundaries: */
        iMin[1] = iMax[0] - 3;
        iMax[1] = XDIM;
        jMin[1] = jMin[0];
        jMax[1] = YDIM / 2 + 2;
        kMin[1] = kMin[0];
        kMax[1] = kMax[0];

        /* Partition 3; ghost cells on IMIN and JMIN boundaries: */
        iMin[2] = iMin[1];
        iMax[2] = iMax[1];
        jMin[2] = jMax[1] - 3;
        jMax[2] = YDIM;
        kMin[2] = kMin[1];
        kMax[2] = kMax[1];
    }
    else
    {
        iMin[0] = 0;
        iMax[0] = XDIM;
        jMin[0] = 0;
        jMax[0] = YDIM;
        kMin[0] = 0;
        kMax[0] = ZDIM;
    }

    /* Local partition dimensions (of equivalent ordered zones) */
    int* iDim = (int*)malloc(effectiveNumPartitions * sizeof(int));
    int* jDim = (int*)malloc(effectiveNumPartitions * sizeof(int));
    int* kDim = (int*)malloc(effectiveNumPartitions * sizeof(int));
    for(int ii = 0; ii < effectiveNumPartitions; ++ii)
    {
        iDim[ii] = iMax[ii] - iMin[ii];
        jDim[ii] = jMax[ii] - jMin[ii];
        kDim[ii] = kMax[ii] - kMin[ii];
    }

    /* Calculate variable and connectivity values for partitions. */
    int** connectivity = (int**)malloc(effectiveNumPartitions * sizeof(int*));
    float** x = (float**)malloc(effectiveNumPartitions * sizeof(float*));
    float** y = (float**)malloc(effectiveNumPartitions * sizeof(float*));
    float** z = (float**)malloc(effectiveNumPartitions * sizeof(float*));
    float** p = (float**)malloc(effectiveNumPartitions * sizeof(float*));

    /* Partition node and cell counts, including ghost items */
    INTEGER4* pNNodes = (INTEGER4*)malloc(effectiveNumPartitions * sizeof(INTEGER4));
    INTEGER4* pNCells = (INTEGER4*)malloc(effectiveNumPartitions * sizeof(INTEGER4));

    for (INTEGER4 ptn = 0; ptn < effectiveNumPartitions; ++ptn)
    {
        pNNodes[ptn] = iDim[ptn] * jDim[ptn] * kDim[ptn];
        if (fOffset == 0)
        {
            /* create grid variables */
            x[ptn] = (float*)malloc(pNNodes[ptn] * sizeof(float));
            y[ptn] = (float*)malloc(pNNodes[ptn] * sizeof(float));
            z[ptn] = (float*)malloc(pNNodes[ptn] * sizeof(float));
            for (int k = 0; k < kDim[ptn]; ++k)
            for (int j = 0; j < jDim[ptn]; ++j)
            for (int i = 0; i < iDim[ptn]; ++i)
            {
                int index = (k * jDim[ptn] + j) * iDim[ptn] + i;
                x[ptn][index] = (float)(i + iMin[ptn] + 1);
                y[ptn][index] = (float)(j + jMin[ptn] + 1);
                z[ptn][index] = (float)(k + kMin[ptn] + 1);
            }
        }
        else
        {
            x[ptn] = 0;
            y[ptn] = 0;
            z[ptn] = 0;
        }

        /* p (cell-centered) and connectivity */
        pNCells[ptn] = (iDim[ptn] - 1) * (jDim[ptn] - 1) * (kDim[ptn] - 1);
        p[ptn] = (float*)malloc(pNCells[ptn] * sizeof(float));
        if (fOffset == 0)
        {
            int connectivityCount = 8 * pNCells[ptn];
            connectivity[ptn] = (INTEGER4*)malloc(connectivityCount * sizeof(INTEGER4));
        }
        else
        {
            connectivity[ptn] = 0;
        }

        for (int k = 0; k < kDim[ptn] - 1; ++k)
        for (int j = 0; j < jDim[ptn] - 1; ++j)
        for (int i = 0; i < iDim[ptn] - 1; ++i)
        {
            int index = (k * (jDim[ptn] - 1) + j) * (iDim[ptn] - 1) + i;
            p[ptn][index] = (float)((i + iMin[ptn] + 1) * (j + jMin[ptn] + 1) * (k + kMin[ptn] + 1) + solTime);

            if (fOffset == 0)
            {
                connectivity[ptn][8 * index] = (k * jDim[ptn] + j) * iDim[ptn] + i + 1; /* One-based to feed TECIO */
                connectivity[ptn][8 * index + 1] = connectivity[ptn][8 * index] + 1;
                connectivity[ptn][8 * index + 2] = connectivity[ptn][8 * index] + iDim[ptn] + 1;
                connectivity[ptn][8 * index + 3] = connectivity[ptn][8 * index] + iDim[ptn];
                connectivity[ptn][8 * index + 4] = connectivity[ptn][8 * index] + iDim[ptn] * jDim[ptn];
                connectivity[ptn][8 * index + 5] = connectivity[ptn][8 * index + 1] + iDim[ptn] * jDim[ptn];
                connectivity[ptn][8 * index + 6] = connectivity[ptn][8 * index + 2] + iDim[ptn] * jDim[ptn];
                connectivity[ptn][8 * index + 7] = connectivity[ptn][8 * index + 3] + iDim[ptn] * jDim[ptn];
            }
        }
    }

    INTEGER4*  nGNodes      = NULL;
    INTEGER4*  nGCells      = NULL;
    INTEGER4** ghostNodes   = NULL;
    INTEGER4** gNPartitions = NULL;
    INTEGER4** gNPNodes     = NULL;
    INTEGER4** ghostCells   = NULL;

    if (numPartitions > 0)
    {
        /* Partition ghost node and ghost cell counts */
        nGNodes      = (INTEGER4*)malloc(numPartitions * sizeof(INTEGER4));
        nGCells      = (INTEGER4*)malloc(numPartitions * sizeof(INTEGER4));

        ghostNodes   = (INTEGER4**)malloc(numPartitions * sizeof(INTEGER4*));
        gNPartitions = (INTEGER4**)malloc(numPartitions * sizeof(INTEGER4*));
        gNPNodes     = (INTEGER4**)malloc(numPartitions * sizeof(INTEGER4*));
        ghostCells   = (INTEGER4**)malloc(numPartitions * sizeof(INTEGER4*));

        GatherGhostNodesAndCells(nGNodes, ghostNodes, gNPartitions, gNPNodes, nGCells, ghostCells, iDim, jDim, kDim, jMin, jMax);
    }

    INTEGER4 dIsDouble = 0;
    for(INTEGER4 ptn = 1; ptn <= effectiveNumPartitions; ++ptn)
    {
        #if defined TECIOMPI
        if (numPartitions == 0 || partitionOwners[ptn - 1] == commRank)
        #endif
        {
            if (numPartitions > 0)
            {
                if (returnValue == 0)
                    returnValue = TECFEPTN142(&ptn, &pNNodes[ptn - 1], &pNCells[ptn - 1],
                        &nGNodes[ptn - 1], ghostNodes[ptn - 1], gNPartitions[ptn - 1], gNPNodes[ptn - 1],
                        &nGCells[ptn - 1], ghostCells[ptn - 1]);

                free(ghostNodes[ptn - 1]);
                free(gNPartitions[ptn - 1]);
                free(gNPNodes[ptn - 1]);
                free(ghostCells[ptn - 1]);
            }

            if (fOffset == 0)
            {
                /* write out the grid field data */
                if (returnValue == 0)
                    returnValue = TECDAT142(&pNNodes[ptn - 1], x[ptn - 1], &dIsDouble);
                if (returnValue == 0)
                    returnValue = TECDAT142(&pNNodes[ptn - 1], y[ptn - 1], &dIsDouble);
                if (returnValue == 0)
                    returnValue = TECDAT142(&pNNodes[ptn - 1], z[ptn - 1], &dIsDouble);
            }

            if (fOffset != 0 || numSolutionFiles == 0)
            {
                /* write out the solution variable */
                if (returnValue == 0)
                    returnValue = TECDAT142(&pNCells[ptn - 1], p[ptn - 1], &dIsDouble);
            }

            if (fOffset == 0)
            {
                /* write out the connectivity */
                INTEGER4 connectivityCount = 8 * pNCells[ptn - 1];
                if (returnValue == 0)
                    returnValue = TECNODE142(&connectivityCount, connectivity[ptn - 1]);
            }

            if (fOffset == 0)
            {
                free(x[ptn - 1]);
                free(y[ptn - 1]);
                free(z[ptn - 1]);
                free(connectivity[ptn - 1]);
            }
            free(p[ptn - 1]);
        }
    }

    free(iMin);
    free(iMax);
    free(jMin);
    free(jMax);
    free(kMin);
    free(kMax);

    free(iDim);
    free(jDim);
    free(kDim);

    free(connectivity);
    free(x);
    free(y);
    free(z);
    free(p);

    free(pNNodes);
    free(pNCells);

    free(nGNodes);
    free(nGCells);
    free(ghostNodes);
    free(gNPartitions);
    free(gNPNodes);
    free(ghostCells);

    return returnValue;
}

/**
 */
static INTEGER4 finalizeFile()
{
    return TECEND142();
}

/**
 */
int main(int argc, char** argv)
{
    int returnValue = 0;

    #if defined TECIOMPI
    {
        MPI_Init(&argc, &argv);
    }
    #endif

    char const* NUM_SOLUTION_FILES_FLAG = "--num-solution-files";
    char* endptr = 0;
    int numSolutionFiles = 0;
    if (argc == 3 && (strcmp(argv[1], NUM_SOLUTION_FILES_FLAG) == 0 &&
                      (strtol(argv[2], &endptr, 10) >= 0 && *endptr == 0)))
    {
        numSolutionFiles = static_cast<int>(strtol(argv[2], &endptr, 10));
    }
    else if (argc != 1)
    {
        returnValue = 1;
        fprintf(stderr, "%s: Expected \"%s <count>, where <count> is zero or greater.\n",
                argv[0], NUM_SOLUTION_FILES_FLAG);
    }

    /*
     * Open the file(s) and write the tecplot datafile header information.
     */
    INTEGER4 const numPartitions = 3; /* 0: non-partitioned file; 3: partitioned */
    assert(numPartitions == 0 || numPartitions == 3); /* ...because this example only handles on or the other */
    INTEGER4* partitionOwners = NULL;
    if (numPartitions > 0)
        partitionOwners = (INTEGER4*)malloc(numPartitions * sizeof(INTEGER4));

    #if defined TECIOMPI
    MPI_Comm mpiComm = MPI_COMM_WORLD;
    int commRank = 0;
    if (returnValue == 0)
    {
        int commSize = 0;
        MPI_Comm_size(mpiComm, &commSize);
        MPI_Comm_rank(mpiComm, &commRank);
        if (commSize == 0)
        {
            returnValue = 1; /* error */
        }
        else
        {
            /* create partition owners */
            for(int ptn = 0; ptn < numPartitions; ++ptn)
                partitionOwners[ptn] = (ptn % commSize);
        }
    }
    #endif

    int const numOutputFiles = 1 + numSolutionFiles;
    for (int fOffset = 0; returnValue == 0 && fOffset < numOutputFiles; ++fOffset)
    {
        double const solTime = static_cast<double>(360 + (fOffset == 0 ? 0 : fOffset-1));

        if (returnValue == 0)
            returnValue = initializeFile(
                #if defined TECIOMPI
                mpiComm,
                #endif
                fOffset, numSolutionFiles, numPartitions, solTime);

        /*
         * Create zones.
         */
        if (returnValue == 0)
            returnValue = createZones(fOffset, numPartitions, partitionOwners, solTime);

        /*
         * Create the connectivity and variable data.
         */
        if (returnValue == 0)
            returnValue = createData(
                #if defined TECIOMPI
                commRank,
                #endif
                fOffset, numSolutionFiles, numPartitions, partitionOwners, solTime);

        if (returnValue == 0)
            returnValue = finalizeFile();
    }

    #if defined TECIOMPI
    {
        MPI_Finalize();
    }
    #endif

    free(partitionOwners);

    return returnValue;
}
