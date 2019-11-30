// This example creates an IJK-ordered zone in 4 partitions.

#if defined TECIOMPI
#include "mpi.h"
#endif

#include <iostream>
#include <vector>

#include "TECIO.h"

#ifndef NULL
#define NULL 0
#endif

#define XDIM 10
#define YDIM 9
#define ZDIM 8

//#define IS_DOUBLE

#if IS_DOUBLE
#define VAR_TYPE double
#else
#define VAR_TYPE float
#endif

INTEGER4 outputVarData(
    VAR_TYPE* var,
    INTEGER4 iDim, INTEGER4 jDim, INTEGER4 kDim,
    INTEGER4 iMin, INTEGER4 jMin, INTEGER4 kMin,
    INTEGER4 iMax, INTEGER4 jMax, INTEGER4 kMax);

int main(int argc, char** argv)
{
    INTEGER4 Debug      = 1;
    INTEGER4 VIsDouble  = 0;
    INTEGER4 FileType   = 0;
    INTEGER4 FileFormat = 1; // SZPLT; .PLT not supported for partitioned zones
    INTEGER4 I          = 0; // Used to track return codes
    #if defined TECIOMPI
        MPI_Init(&argc, &argv);
        MPI_Comm mpiComm = MPI_COMM_WORLD;
        int commSize;
        MPI_Comm_size(mpiComm, &commSize);
        int commRank;
        MPI_Comm_rank(mpiComm, &commRank);
        #if defined _DEBUG
            if (commRank == 0)
            {
                std::cout << "Press return to continue" << std::endl;
                std::cin.get();
            }
        #endif
    #endif

    /*
     * Open the file and write the tecplot datafile
     * header information
     */
    I = TECINI142((char*)"IJK Ordered Zone",
                  (char*)"X Y Z P",
                  (char*)"ijkpartitioned.szplt",
                  (char*)".",
                  &FileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);

    #if defined TECIOMPI
        INTEGER4 mainRank = 0;
        I = TECMPIINIT142(&mpiComm, &mainRank);
    #endif

    float* x = new float[XDIM * YDIM * ZDIM];
    float* y = new float[XDIM * YDIM * ZDIM];
    float* z = new float[XDIM * YDIM * ZDIM];
    float* p = new float[(XDIM - 1) * (YDIM - 1) * (ZDIM - 1)];

    INTEGER4 IMax                     = XDIM;
    INTEGER4 JMax                     = YDIM;
    INTEGER4 KMax                     = ZDIM;
    INTEGER4 ICellMax                 = 0;
    INTEGER4 JCellMax                 = 0;
    INTEGER4 KCellMax                 = 0;
    INTEGER4 DIsDouble                = 0;
    double   SolTime                  = 360.0;
    INTEGER4 StrandID                 = 0;      /* StaticZone */
    INTEGER4 ParentZn                 = 0;
    INTEGER4 IsBlock                  = 1;      /* Block */
    INTEGER4 NFConns                  = 0;
    INTEGER4 FNMode                   = 0;
    INTEGER4 TotalNumFaceNodes        = 1;
    INTEGER4 TotalNumBndryFaces       = 1;
    INTEGER4 TotalNumBndryConnections = 1;
    INTEGER4 valueLocations[]         = {1, 1, 1, 0};
    INTEGER4 ShrConn                  = 0;

    // Ordered Zone Parameters
    for(int i = 0; i < XDIM; ++i)
        for(int j = 0; j < YDIM; ++j)
            for(int k = 0; k < ZDIM; ++k)
            {
                int index = (k * YDIM + j) * XDIM + i;
                x[index] = (float)i;
                y[index] = (float)j;
                z[index] = (float)k;
                if (i < XDIM - 1 && j < YDIM - 1 && k < ZDIM - 1)
                {
                    int cindex = (k * (YDIM - 1) + j) * (XDIM - 1) + i;
                    p[cindex] = (x[index] + 0.5f) * (y[index] + 0.5f) * (z[index] + 0.5f);
                }
            }

    //  Create the ordered Zone
    INTEGER4 ZoneType = 0;
    I = TECZNE142((char*)"Ordered Zone",
                  &ZoneType,
                  &IMax,
                  &JMax,
                  &KMax,
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
                  &TotalNumBndryFaces,
                  &TotalNumBndryConnections,
                  NULL, // PassiveVarList
                  valueLocations, // ValueLocation
                  NULL, // ShareVarFromZone
                  &ShrConn);

    /*
     * The nodal index ranges of each partition.
     * Interior boundary index ranges must coincide. For example, if one partition's I index
     * range goes from 1 to 3, its neighboring partition must start with an I index of 3.
     */
    INTEGER4 partitionIndices[4][6] = {
        { 1, 1, 1, XDIM / 3 + 1, YDIM, ZDIM },
        { XDIM / 3 + 1, 1, 1, XDIM, YDIM / 3 + 1, ZDIM },
        { XDIM / 3 + 1, YDIM / 3 + 1, 1, XDIM, YDIM, ZDIM / 2 + 1 },
        { XDIM / 3 + 1, YDIM / 3 + 1, ZDIM / 2 + 1, XDIM, YDIM, ZDIM }
    };

    // Output partitions
    #if defined TECIOMPI
        INTEGER4 numPartitions = 4;
        std::vector<INTEGER4> partitionOwners;
        for (INTEGER4 ptn = 0; ptn < numPartitions; ++ptn)
            partitionOwners.push_back(ptn % commSize);
        TECZNEMAP142(&numPartitions, &partitionOwners[0]);
    #endif

    for (INTEGER4 partition = 1; partition <= 4; ++partition)
    {
        #if defined TECIOMPI
            if (partitionOwners[partition - 1] == commRank)
            {
        #endif
        INTEGER4 partitionIMin = partitionIndices[partition - 1][0];
        INTEGER4 partitionJMin = partitionIndices[partition - 1][1];
        INTEGER4 partitionKMin = partitionIndices[partition - 1][2];
        INTEGER4 partitionIMax = partitionIndices[partition - 1][3];
        INTEGER4 partitionJMax = partitionIndices[partition - 1][4];
        INTEGER4 partitionKMax = partitionIndices[partition - 1][5];
        I = TECIJKPTN142(&partition, &partitionIMin, &partitionJMin, &partitionKMin, &partitionIMax, &partitionJMax, &partitionKMax);
        I = outputVarData(x, XDIM, YDIM, ZDIM, partitionIMin, partitionJMin, partitionKMin, partitionIMax, partitionJMax, partitionKMax);
        I = outputVarData(y, XDIM, YDIM, ZDIM, partitionIMin, partitionJMin, partitionKMin, partitionIMax, partitionJMax, partitionKMax);
        I = outputVarData(z, XDIM, YDIM, ZDIM, partitionIMin, partitionJMin, partitionKMin, partitionIMax, partitionJMax, partitionKMax);
        I = outputVarData(p, XDIM - 1, YDIM - 1, ZDIM - 1, partitionIMin, partitionJMin, partitionKMin, partitionIMax - 1, partitionJMax - 1, partitionKMax - 1);
        #if defined TECIOMPI
            }
        #endif
    }

    // The I=IMax boundary is output as three unpartitioned surface zones.
    // For MPI, each of these is output by the MPI rank that also output
    // the partition that it boundaries (partitions 2 through 4 touch the
    // I=IMax boundary).
    for (INTEGER4 partition = 2; partition <= 4; ++partition)
    {
        #if defined TECIOMPI
            // Only the main rank and the rank that owns the zone need call TECZNE.
            // The others may, but it becomes a no-op for them.
            if (commRank == mainRank || partitionOwners[partition - 1] == commRank)
            {
        #endif
        INTEGER4 surfaceZoneJMin = partitionIndices[partition - 1][1];
        INTEGER4 surfaceZoneKMin = partitionIndices[partition - 1][2];
        INTEGER4 surfaceZoneJMax = partitionIndices[partition - 1][4];
        INTEGER4 surfaceZoneKMax = partitionIndices[partition - 1][5];
        IMax = 1;
        JMax = surfaceZoneJMax - surfaceZoneJMin + 1;
        KMax = surfaceZoneKMax - surfaceZoneKMin + 1;
        I = TECZNE142((char*)"Ordered Zone",
            &ZoneType,
            &IMax,
            &JMax,
            &KMax,
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
            &TotalNumBndryFaces,
            &TotalNumBndryConnections,
            NULL, // PassiveVarList
            valueLocations, // ValueLocation
            NULL, // ShareVarFromZone
            &ShrConn);
        #if defined TECIOMPI
            INTEGER4 one = 1;
            I = TECZNEMAP142(&one, &partitionOwners[partition - 1]);
        #endif
        I = outputVarData(x, XDIM, YDIM, ZDIM,
            XDIM, surfaceZoneJMin, surfaceZoneKMin, XDIM, surfaceZoneJMax, surfaceZoneKMax);
        I = outputVarData(y, XDIM, YDIM, ZDIM,
            XDIM, surfaceZoneJMin, surfaceZoneKMin, XDIM, surfaceZoneJMax, surfaceZoneKMax);
        I = outputVarData(z, XDIM, YDIM, ZDIM,
            XDIM, surfaceZoneJMin, surfaceZoneKMin, XDIM, surfaceZoneJMax, surfaceZoneKMax);
        // Output the cell-centered values for the iMax - .5 plane of cells
        I = outputVarData(p, XDIM - 1, YDIM - 1, ZDIM - 1,
            XDIM - 1, surfaceZoneJMin, surfaceZoneKMin, XDIM - 1, surfaceZoneJMax - 1, surfaceZoneKMax - 1);
        #if defined TECIOMPI
            }
        #endif
    }

    I = TECEND142();
    #if defined TECIOMPI
        MPI_Finalize();
    #endif

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] p;

    return 0;
}

INTEGER4 outputVarData(
    VAR_TYPE* var,
    INTEGER4 iDim, INTEGER4 jDim, INTEGER4 kDim,
    INTEGER4 iMin, INTEGER4 jMin, INTEGER4 kMin,
    INTEGER4 iMax, INTEGER4 jMax, INTEGER4 kMax)
{
    INTEGER4 count = iMax - iMin + 1;
#if defined IS_DOUBLE
    INTEGER4 isDouble = 1;
#else
    INTEGER4 isDouble = 0;
#endif
    INTEGER4 result = 0;
    for (INTEGER4 k = kMin; result == 0 && k <= kMax; ++k)
    {
        for (INTEGER4 j = jMin; result == 0 && j <= jMax; ++j)
        {
            INTEGER4 index = ((k - 1) * jDim + j - 1) * iDim + iMin - 1;
            result = TECDAT142(&count, &var[index], &isDouble);
        }
    }
    return result;
}

