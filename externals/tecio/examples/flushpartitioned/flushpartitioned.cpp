/*
 * Simple example c program to write a
 * binary datafile for tecplot.  This example
 * does the following:
 *
 *   1.  Open a data file called "flush.szplt"
 *   2.  Assign values for x, y, and z.
 *   3.  Loop over time steps
 *       a. Assign values for p
 *       b. Flush the data to file for the time step.
 *   4.  Close the datafile.
 *
 * If TECIOMPI is #defined, this program may be executed with mpiexec with
 * up to 3 MPI ranks (processes). In this case, it must be linked
 * with the MPI version of TECIO.
 *
 * If OLD_API is defined, the old TECXXX routines will be used. Otherwise,
 * the new tec... routines will be used.
 *
 * Note that for the new API, there are 64-bit versions of several functions. To enable
 * us to demonstrate both old and new APIs here, we'll use the 32-bit versions.
 */

// Internal testing flags
// RUNFLAGS:none

#if defined TECIOMPI
    #include <mpi.h>
#endif

#include "TECIO.h"
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

// Uncomment the following line to use the old TECXXX routines:
//#define OLD_API

#if defined _WIN32 && !defined snprintf
    #define snprintf _snprintf
#endif

#ifndef NULL
    #define NULL 0
#endif

/* The zone will appear the same as an ordered zone with these dimensions: */
#define XDIM 10
#define YDIM 9
#define ZDIM 8

int32_t initializeFile(
    #if defined TECIOMPI
        MPI_Comm mpiComm,
    #endif
    int numPartitions,
    void** fileHandle);

int32_t createZone(
    int             timeStep,
    int32_t         numPartitions,
    int32_t const*  partitionOwners,
    double          solTime,
    void*           fileHandle,
    int32_t*        zone);

int32_t createData(
    #if defined TECIOMPI
        int             commRank,
    #endif
    int             timeStep,
    int32_t         numPartitions,
    int32_t const*  partitionOwners,
    double          solTime,
    void*           fileHandle,
    int32_t         zone);

int32_t finalizeFile(void** fileHandle);

/**
 */
int main(int argc, char** argv)
{
    int returnValue = 0;

    #if defined TECIOMPI
        MPI_Init(&argc, &argv);
    #endif

    /*
     * Open the file(s) and write the tecplot datafile header information.
     */
    int32_t const numPartitions = 3; /* 0: non-partitioned file; 3: partitioned */
    assert(numPartitions == 0 || numPartitions == 3); /* ...because this example handles only on or the other */
    int32_t* partitionOwners = NULL;
    if (numPartitions > 0)
        partitionOwners = (int32_t*)malloc(numPartitions * sizeof(int32_t));

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

    void* fileHandle = NULL;
    if (returnValue == 0)
        returnValue = initializeFile(
            #if defined TECIOMPI
                mpiComm,
            #endif
            numPartitions, &fileHandle);

    if (returnValue == 0)
    {
        // Add custom labels to the file
        #if defined TECIOMPI
            if (commRank == 0)
        #endif
            if (returnValue == 0)
            #if defined OLD_API
                returnValue = TECLAB142("Sunday,Monday,Tuesday,Wednesday,Thursday,Friday,Saturday");
            #else
                returnValue = tecCustomLabelsAddSet(fileHandle, "Sunday,Monday,Tuesday,Wednesday,Thursday,Friday,Saturday");
            #endif
    }

    int const numTimeSteps = 5;
    int32_t   zone = 0;
    for (int timeStep = 0; returnValue == 0 && timeStep < numTimeSteps; ++timeStep)
    {
        double const solTime = static_cast<double>(360 + timeStep);

        /*
         * Create zone.
         */
        if (returnValue == 0)
            returnValue = createZone(timeStep, numPartitions, partitionOwners, solTime, fileHandle, &zone);

        /*
         * Create the connectivity and variable data.
         */
        if (returnValue == 0)
            returnValue = createData(
                #if defined TECIOMPI
                    commRank,
                #endif
                timeStep, numPartitions, partitionOwners, solTime, fileHandle, zone);

    }
    if (returnValue == 0)
        returnValue = finalizeFile(&fileHandle);

    #if defined TECIOMPI
        MPI_Finalize();
    #endif

    free(partitionOwners);

    return returnValue;
}

/**
 * Open the file and write the tecplot datafile header information.
 */
int32_t initializeFile(
    #if defined TECIOMPI
        MPI_Comm mpiComm,
    #endif
    int numPartitions,
    void** fileHandle)
{
    size_t const fileNameSize = 1024;
    char fileName[fileNameSize];
    char const* baseFileName = numPartitions > 0 ? "flushpartitioned" : "flush";
    char const* varNames;
    snprintf(fileName, fileNameSize, "%s.szplt", baseFileName);
    varNames = "x y z p";

    #if defined OLD_API
        int32_t fileFormat = FILEFORMAT_SZL;
        int32_t fileType = FILETYPE_FULL;
        int32_t debug = 1;
        int32_t vIsDouble = 1;
        int32_t returnValue = TECINI142("Flush Example Data Set", (char*)varNames, fileName, ".", &fileFormat, &fileType, &debug, &vIsDouble);
    #else
        int32_t returnValue = tecFileWriterOpen(fileName, "Flush Example Data Set", varNames, FILEFORMAT_SZL, FILETYPE_FULL, FieldDataType_Double, NULL, fileHandle);
        if (returnValue == 0)
            returnValue = tecFileSetDiagnosticsLevel(*fileHandle, 1);
    #endif

    #if defined TECIOMPI
        if (returnValue == 0)
        {
            int32_t mainRank = 0;
            #if defined OLD_API
                returnValue = TECMPIINIT142(&mpiComm, &mainRank);
            #else
                returnValue = tecMPIInitialize(*fileHandle, mpiComm, mainRank);
            #endif
        }
    #endif

    return returnValue;
}

/**
 */
int32_t createZone(
    int             timeStep,
    int32_t         numPartitions,
    int32_t const*  partitionOwners,
    double          solTime,
    void*           fileHandle,
    int32_t*        zone)
{
    #if defined OLD_API
        int32_t nNodes    = XDIM * YDIM * ZDIM; /* Overall zone dimensions */
        int32_t nCells    = (XDIM - 1) * (YDIM - 1) * (ZDIM - 1);
    #else
        int64_t nNodes    = XDIM * YDIM * ZDIM; /* Overall zone dimensions */
        int64_t nCells    = (XDIM - 1) * (YDIM - 1) * (ZDIM - 1);
    #endif
    int32_t strandID  = 1;
    int32_t valueLocations[4] = {1, 1, 1, 0}; /* 1 = Nodal, 0 = Cell-Centered */
    int32_t shareVarFromZone[4] = {1, 1, 1, 0};
    int32_t shrConn   = (timeStep == 0 ? 0 : 1);

    size_t const zoneTitleSize = 1024;
    char zoneTitle[zoneTitleSize];
    char const* baseTitle = numPartitions > 0 ? "Partitioned Zone" : "Zone";
    snprintf(zoneTitle, zoneTitleSize, "%s Time=%d", baseTitle, static_cast<int>(std::ceil(solTime)));
    int32_t* shareVarPtr = (timeStep == 0 ? NULL : shareVarFromZone);

    #if defined OLD_API
        int32_t zoneType  = 5;      /* Brick */
        int32_t nFaces    = 6;         /* Not used */
        int32_t iCellMax  = 0;
        int32_t jCellMax  = 0;
        int32_t kCellMax  = 0;
        int32_t parentZn  = 0;      /* No Parent */
        int32_t isBlock   = 1;      /* Block */
        int32_t nFConns   = 0;
        int32_t fNMode    = 0;
        int32_t zero = 0;
        int32_t returnValue = TECZNE142(zoneTitle, &zoneType, &nNodes, &nCells, &nFaces, &iCellMax, &jCellMax, &kCellMax, &solTime, &strandID, &parentZn, &isBlock,
            &nFConns, &fNMode, &zero, &zero, &zero, NULL, valueLocations, shareVarPtr, &shrConn);
        ++*zone; 
    #else
        // Note: There are 64-bit versions of this and related routines for zones that overflow 32-bit limits
        int32_t returnValue = tecZoneCreateFE(fileHandle, zoneTitle, ZONETYPE_FEBRICK, nNodes, nCells, NULL, shareVarPtr, valueLocations, NULL, shrConn, 0, 0, zone);
        if (returnValue == 0)
            returnValue = tecZoneSetUnsteadyOptions(fileHandle, *zone, solTime, strandID);
    #endif

    if (returnValue == 0)
    {
        #if defined TECIOMPI
            /* Output partitions */
            if (numPartitions > 0 && returnValue == 0)
                #if defined OLD_API
                    returnValue = TECZNEMAP142(&numPartitions, partitionOwners);
                #else
                    returnValue = tecZoneMapPartitionsToMPIRanks(fileHandle, *zone, numPartitions, partitionOwners);
                #endif
        #endif
    }

    return returnValue;
}

/**
 */
static void appendGhostItems(
    int32_t* nGhosts, int32_t** ghosts, int32_t** gPartitions, int32_t** gPGhosts, int32_t ptn,
    int32_t gIDim, int32_t gJDim,
    int32_t gIStart, int32_t gIEnd, int32_t gJStart, int32_t gJEnd, int32_t gKStart, int32_t gKEnd,
    int32_t oIDim, int32_t oJDim,
    int32_t oIStart, int32_t oIEnd, int32_t oJStart, int32_t oJEnd, int32_t oKStart, int32_t oKEnd)
{
    int32_t localNumGhosts = (gIEnd - gIStart) * (gJEnd - gJStart) * (gKEnd - gKStart);
    int32_t gOffset = *nGhosts;
    *nGhosts += localNumGhosts;
    *ghosts = (int32_t*)realloc(*ghosts, *nGhosts * sizeof(int32_t));
    if (gPartitions)
        *gPartitions = (int32_t*)realloc(*gPartitions, *nGhosts * sizeof(int32_t));
    if (gPGhosts)
        *gPGhosts = (int32_t*)realloc(*gPGhosts, *nGhosts * sizeof(int32_t));

    for(int32_t i = gIStart; i < gIEnd; ++i)
    for(int32_t j = gJStart; j < gJEnd; ++j)
    for(int32_t k = gKStart; k < gKEnd; ++k)
    {
        (*ghosts)[gOffset] = (k * gJDim + j) * gIDim + i + 1;
        if (gPartitions)
            (*gPartitions)[gOffset] = ptn;
        if (gPGhosts)
        {
            int32_t oI = i - gIStart + oIStart;
            int32_t oJ = j - gJStart + oJStart;
            int32_t oK = k - gKStart + oKStart;
            (*gPGhosts)[gOffset] = (oK * oJDim + oJ) * oIDim + oI + 1;
        }
        ++gOffset;
    }
}


/**
 */
static void GatherGhostNodesAndCells(
    int32_t* nGNodes, int32_t** ghostNodes, int32_t** gNPartitions, int32_t** gNPNodes,
    int32_t* nGCells, int32_t** ghostCells, int32_t* iDim, int32_t* jDim, int32_t* kDim, int32_t* jMin, int32_t* jMax)
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
int32_t createData(
    #if defined TECIOMPI
    int         commRank,
    #endif
    int             timeStep,
    int32_t         numPartitions,
    int32_t const*  partitionOwners,
    double          solTime,
    void*           fileHandle,
    int32_t         zone)
{
    int32_t returnValue = 0;

    /*
     * Add zone aux data.
     */
    char timeStepString[40];
    snprintf(timeStepString, sizeof(timeStepString), "%d", timeStep);
    #if defined TECIOMPI
        if (commRank == 0)
    #endif
    #if defined OLD_API
        returnValue = TECZAUXSTR142("TimeStep", timeStepString);
    #else
        returnValue = tecZoneAddAuxData(fileHandle, zone, "TimeStep", timeStepString);
    #endif

    /*
     * Add 3D text attached to the zone
     */
    char timeStepLabel[40];
    snprintf(timeStepLabel, sizeof(timeStepLabel), "Time Step = %d", timeStep);
    #if defined TECIOMPI
        if (commRank == 0)
    #endif
    {
    #if defined OLD_API
        double      xPos             = 5.0;
        double      yPos             = 0.0;
        double      zPos             = 8.0;
        INTEGER4    posCoordMode     = 6; // Grid3D
        INTEGER4    attachToZone     = 1;
        INTEGER4    font             = 0; // Helvetica
        INTEGER4    fontHeightUnits  = 2; // Points
        double      fontHeight       = 18.0;
        INTEGER4    boxType          = 0; // No box
        double      boxMargin        = 5.0;
        double      boxLineThickness = .01;
        INTEGER4    boxColor         = 0;
        INTEGER4    boxFillColor     = 7;
        double      angle            = 0.0;
        INTEGER4    anchor           = 0;
        double      lineSpacing      = 1.0;
        INTEGER4    textColor        = 0;
        INTEGER4    scope            = 0;
        INTEGER4    clipping         = 0;
        char const* MFC              = "";

        returnValue = TECTXT142(&xPos, &yPos, &zPos, &posCoordMode, &attachToZone, &zone,
            &font, &fontHeightUnits, &fontHeight,
            &boxType, &boxMargin, &boxLineThickness, &boxColor, &boxFillColor,
            &angle, &anchor, &lineSpacing,
            &textColor, &scope, &clipping, timeStepLabel, MFC);
    #else
        if (returnValue == 0)
            returnValue = tecText3DBegin(fileHandle, timeStepLabel, 5.0, 0.0, 8.0, 18.0, 2);
        if (returnValue == 0)
            returnValue = tecTextAttachToZone(fileHandle, zone);
        if (returnValue == 0)
            returnValue = tecTextEnd(fileHandle);
    #endif
    }

    /*
     * Add a circle attached to the zone
     */
    #if defined TECIOMPI
        if (commRank == 0)
    #endif
    {
    #if defined OLD_API
        double      xPos                = 3.0;
        double      yPos                = 2.0 + 0.1 * timeStep;
        double      zPos                = 0.0; // Unused for circles
        INTEGER4    posCoordMode        = 0; // Grid
        INTEGER4    attachToZone        = 1;
        INTEGER4    color               = 0;
        INTEGER4    fillColor           = 7;
        INTEGER4    isFilled            = 0;
        INTEGER4    geomType            = 3; // Circle
        INTEGER4    linePattern         = 0; // Solid
        double      patternLength       = 1.0;
        double      lineThicknessness   = 0.1;
        INTEGER4    numEllipsePts       = 72;
        INTEGER4    arrowheadStyle      = 0;
        INTEGER4    arrowheadAttachment = 0;
        double      arrowheadSize       = 5.0;
        double      arrowheadAngle      = 12.0;
        INTEGER4    scope               = 1;
        INTEGER4    clipping            = 0;
        INTEGER4    numSegments         = 0;
        INTEGER4*   numSegPts           = NULL;
        float       xOrThetaGeomData[]  = { 1.0 };
        float*      yOrRGeomData        = NULL;
        float*      zGeomData           = NULL;
        char const* MFC                 = "";
        returnValue = TECGEO142(&xPos, &yPos, &zPos, &posCoordMode, &attachToZone, &zone,
            &color, &fillColor, &isFilled, &geomType, &linePattern, &patternLength, &lineThicknessness,
            &numEllipsePts, &arrowheadStyle, &arrowheadAttachment, &arrowheadSize, &arrowheadAngle,
            &scope, &clipping, &numSegments, numSegPts, xOrThetaGeomData, yOrRGeomData, zGeomData, MFC);
    #else
        if (returnValue == 0)
            returnValue = tecGeomCircleBegin(fileHandle, 3.0, 2.0 + 0.1 * timeStep, 1.0, 0 /* grid-mode */);
        if (returnValue == 0)
            returnValue = tecGeomAttachToZone(fileHandle, zone);
        if (returnValue == 0)
            returnValue = tecGeomEnd(fileHandle);
    #endif
    }

    /*
     * If the file isn't partitioned, we still allocate out array resources as if there was a single
     * partition as it makes the code between partitioned and non-partitioned the same.
     */
    int32_t const effectiveNumPartitions = numPartitions > 0 ? numPartitions : 1;

    /*
     * Divide the zone into number of partitions, identified by the index ranges
     * of an equivalent unpartitioned ordered zone.
     */
    int32_t* iMin = (int32_t*)malloc(effectiveNumPartitions * sizeof(int32_t));
    int32_t* iMax = (int32_t*)malloc(effectiveNumPartitions * sizeof(int32_t));
    int32_t* jMin = (int32_t*)malloc(effectiveNumPartitions * sizeof(int32_t));
    int32_t* jMax = (int32_t*)malloc(effectiveNumPartitions * sizeof(int32_t));
    int32_t* kMin = (int32_t*)malloc(effectiveNumPartitions * sizeof(int32_t));
    int32_t* kMax = (int32_t*)malloc(effectiveNumPartitions * sizeof(int32_t));

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
    int32_t* iDim = (int32_t*)malloc(effectiveNumPartitions * sizeof(int32_t));
    int32_t* jDim = (int32_t*)malloc(effectiveNumPartitions * sizeof(int32_t));
    int32_t* kDim = (int32_t*)malloc(effectiveNumPartitions * sizeof(int32_t));
    for(int ii = 0; ii < effectiveNumPartitions; ++ii)
    {
        iDim[ii] = iMax[ii] - iMin[ii];
        jDim[ii] = jMax[ii] - jMin[ii];
        kDim[ii] = kMax[ii] - kMin[ii];
    }

    /* Calculate variable and connectivity values for partitions. */
    int32_t** connectivity = (int32_t**)malloc(effectiveNumPartitions * sizeof(int32_t*));
    float** x = (float**)malloc(effectiveNumPartitions * sizeof(float*));
    float** y = (float**)malloc(effectiveNumPartitions * sizeof(float*));
    float** z = (float**)malloc(effectiveNumPartitions * sizeof(float*));
    float** p = (float**)malloc(effectiveNumPartitions * sizeof(float*));

    /* Partition node and cell counts, including ghost items */
    int32_t* pNNodes = (int32_t*)malloc(effectiveNumPartitions * sizeof(int32_t));
    int32_t* pNCells = (int32_t*)malloc(effectiveNumPartitions * sizeof(int32_t));

    for (int32_t ptn = 0; ptn < effectiveNumPartitions; ++ptn)
    {
        pNNodes[ptn] = iDim[ptn] * jDim[ptn] * kDim[ptn];
        if (timeStep == 0)
        {
            /* create grid variables */
            x[ptn] = (float*)malloc(pNNodes[ptn] * sizeof(float));
            y[ptn] = (float*)malloc(pNNodes[ptn] * sizeof(float));
            z[ptn] = (float*)malloc(pNNodes[ptn] * sizeof(float));
            for (int32_t k = 0; k < kDim[ptn]; ++k)
            for (int32_t j = 0; j < jDim[ptn]; ++j)
            for (int32_t i = 0; i < iDim[ptn]; ++i)
            {
                int32_t index = (k * jDim[ptn] + j) * iDim[ptn] + i;
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
        if (timeStep == 0)
        {
            int32_t connectivityCount = 8 * pNCells[ptn];
            connectivity[ptn] = (int32_t*)malloc(connectivityCount * sizeof(int32_t));
        }
        else
        {
            connectivity[ptn] = 0;
        }

        for (int32_t k = 0; k < kDim[ptn] - 1; ++k)
        for (int32_t j = 0; j < jDim[ptn] - 1; ++j)
        for (int32_t i = 0; i < iDim[ptn] - 1; ++i)
        {
            int32_t index = (k * (jDim[ptn] - 1) + j) * (iDim[ptn] - 1) + i;
            p[ptn][index] = (float)(solTime + (i + iMin[ptn] + 1) * (j + jMin[ptn] + 1) * (k + kMin[ptn] + 1));

            if (timeStep == 0)
            {
                connectivity[ptn][8 * index] = (k * jDim[ptn] + j) * iDim[ptn] + i + 1; /* One-based to feed TECIO old API */
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

    int32_t*  nGNodes      = NULL;
    int32_t*  nGCells      = NULL;
    int32_t** ghostNodes   = NULL;
    int32_t** gNPartitions = NULL;
    int32_t** gNPNodes     = NULL;
    int32_t** ghostCells   = NULL;

    if (numPartitions > 0)
    {
        /* Partition ghost node and ghost cell counts */
        nGNodes      = (int32_t*)malloc(numPartitions * sizeof(int32_t));
        nGCells      = (int32_t*)malloc(numPartitions * sizeof(int32_t));

        ghostNodes   = (int32_t**)malloc(numPartitions * sizeof(int32_t*));
        gNPartitions = (int32_t**)malloc(numPartitions * sizeof(int32_t*));
        gNPNodes     = (int32_t**)malloc(numPartitions * sizeof(int32_t*));
        ghostCells   = (int32_t**)malloc(numPartitions * sizeof(int32_t*));

        GatherGhostNodesAndCells(nGNodes, ghostNodes, gNPartitions, gNPNodes, nGCells, ghostCells, iDim, jDim, kDim, jMin, jMax);
    }

    int32_t dIsDouble = 0;
    for(int32_t ptn = 1; returnValue == 0 && ptn <= effectiveNumPartitions; ++ptn)
    {
        #if defined TECIOMPI
            if (numPartitions == 0 || partitionOwners[ptn - 1] == commRank)
        #endif
        {
            if (numPartitions > 0)
            {
                if (returnValue == 0)
                    #if defined OLD_API
                        returnValue = TECFEPTN142(
                            &ptn,
                            &pNNodes[ptn - 1],
                            &pNCells[ptn - 1],
                            &nGNodes[ptn - 1],
                            ghostNodes[ptn - 1],
                            gNPartitions[ptn - 1],
                            gNPNodes[ptn - 1],
                            &nGCells[ptn - 1],
                            ghostCells[ptn - 1]);
                    #else
                        returnValue = tecFEPartitionCreate32(
                            fileHandle,
                            zone,
                            ptn,
                            pNNodes[ptn - 1],
                            pNCells[ptn - 1],
                            nGNodes[ptn - 1],
                            ghostNodes[ptn - 1],
                            gNPartitions[ptn - 1],
                            gNPNodes[ptn - 1],
                            nGCells[ptn - 1],
                            ghostCells[ptn - 1]);
                    #endif

                free(ghostNodes[ptn - 1]);
                free(gNPartitions[ptn - 1]);
                free(gNPNodes[ptn - 1]);
                free(ghostCells[ptn - 1]);
            }

            if (timeStep == 0)
            {
                /* write out the grid field data */
                #if defined OLD_API
                    if (returnValue == 0)
                        returnValue = TECDAT142(&pNNodes[ptn - 1], x[ptn - 1], &dIsDouble);
                    if (returnValue == 0)
                        returnValue = TECDAT142(&pNNodes[ptn - 1], y[ptn - 1], &dIsDouble);
                    if (returnValue == 0)
                        returnValue = TECDAT142(&pNNodes[ptn - 1], z[ptn - 1], &dIsDouble);
                #else
                    if (returnValue == 0)
                        returnValue = tecZoneVarWriteFloatValues(fileHandle, zone, 1, ptn, pNNodes[ptn - 1], x[ptn - 1]);
                    if (returnValue == 0)
                        returnValue = tecZoneVarWriteFloatValues(fileHandle, zone, 2, ptn, pNNodes[ptn - 1], y[ptn - 1]);
                    if (returnValue == 0)
                        returnValue = tecZoneVarWriteFloatValues(fileHandle, zone, 3, ptn, pNNodes[ptn - 1], z[ptn - 1]);
                #endif
            }

            /* write out the solution variable */
            if (returnValue == 0)
                #if defined OLD_API
                    returnValue = TECDAT142(&pNCells[ptn - 1], p[ptn - 1], &dIsDouble);
                #else
                    returnValue = tecZoneVarWriteFloatValues(fileHandle, zone, 4, ptn, pNCells[ptn - 1], p[ptn - 1]);
                #endif

            if (timeStep == 0)
            {
                /* write out the connectivity */
                #if defined OLD_API
                    int32_t connectivityCount = 8 * pNCells[ptn - 1];
                    if (returnValue == 0)
                        returnValue = TECNODE142(&connectivityCount, connectivity[ptn - 1]);
                #else
                    int32_t connectivityCount = 8 * pNCells[ptn - 1];
                    if (returnValue == 0)
                        returnValue = tecZoneNodeMapWrite32(fileHandle, zone, ptn, 1, connectivityCount, connectivity[ptn - 1]);
                #endif
            }

            if (timeStep == 0)
            {
                free(x[ptn - 1]);
                free(y[ptn - 1]);
                free(z[ptn - 1]);
                free(connectivity[ptn - 1]);
            }
            free(p[ptn - 1]);
        }
    }
    int32_t numZonesToRetain = 1;
    int32_t zonesToRetain[] = {1};
    if (returnValue == 0)
        #if defined OLD_API
            returnValue = TECFLUSH142(&numZonesToRetain, zonesToRetain);
        #else
            returnValue = tecFileWriterFlush(fileHandle, numZonesToRetain, zonesToRetain);
        #endif

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
int32_t finalizeFile(void** fileHandle)
{
    #if defined OLD_API
        return TECEND142();
    #else
        return tecFileWriterClose(fileHandle);
    #endif
}
