/*
 * Simple example c program to write a
 * binary datafile for tecplot.  This example
 * does the following:
 *
 *   1.  Open a datafile called "t.plt"
 *   2.  Assign values for x,y, and p.
 *   3.  Write out a tetrahedral zone.
 *   4.  Close the datafile.
 */

// Internal testing flags
// RUNFLAGS:none
// RUNFLAGS:--szl

#include "TECIO.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef NULL
#define NULL 0
#endif

/* The zone will have the same nodes as an ordered zone with these dimensions: */
#define XDIM 5
#define YDIM 4
#define ZDIM 3

enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };

int main(int argc, const char *argv[])
{
    float *x, *y, *z, *p;
    int *connectivity;
    double solTime;
    INTEGER4 debug, i, j, k, dIsDouble, vIsDouble, zoneType, strandID, parentZn, isBlock;
    INTEGER4 iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType;


    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    if (argc == 2 && strncmp(argv[1],"--szl",5) == 0)
        fileFormat = 1; 
    else
        fileFormat = 0; 

    INTEGER4 nNodes, nCells, nFaces, connectivityCount, index;

    debug     = 1;
    vIsDouble = 0;
    dIsDouble = 0;
    nNodes = XDIM * YDIM * ZDIM;
    nCells = (XDIM - 1) * (YDIM - 1) * (ZDIM - 1) * 5;
    nFaces = 6; /* Not used */
    zoneType  = 4;      /* Tetrahedral */
    solTime   = 360.0;
    strandID  = 0;     /* StaticZone */
    parentZn  = 0;      /* No Parent */
    isBlock   = 1;      /* Block */
    iCellMax  = 0;
    jCellMax  = 0;
    kCellMax  = 0;
    nFConns   = 0;
    fNMode    = 0;
    shrConn   = 0;
    fileType  = FULL;

    /*
     * Open the file and write the tecplot datafile
     * header information
     */
    i = TECINI142((char*)"SIMPLE DATASET",
                  (char*)"x y z p",
                  (char*)"tetra.plt",
                  (char*)".",
                  &fileFormat,
                  &fileType,
                  &debug,
                  &vIsDouble);

    x = (float*)malloc(nNodes * sizeof(float));
    y = (float*)malloc(nNodes * sizeof(float));
    z = (float*)malloc(nNodes * sizeof(float));
    p = (float*)malloc(nNodes * sizeof(float));
    for (k = 0; k < ZDIM; k++)
        for (j = 0; j < YDIM; j++)
            for (i = 0; i < XDIM; i++)
            {
                index = (k * YDIM + j) * XDIM + i;
                x[index] = (float)(i + 1);
                y[index] = (float)(j + 1);
                z[index] = (float)(k + 1);
                p[index] = (float)((i + 1) * (j + 1) * (k + 1));
            }

    connectivityCount = 4 * nCells;
    connectivity = (INTEGER4*)malloc(connectivityCount * sizeof(INTEGER4));
    for (k = 0; k < ZDIM - 1; k++)
        for (j = 0; j < YDIM - 1; j++)
            for (i = 0; i < XDIM - 1; i++)
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
                int evenOdd = (i + j + k) % 2;
                int whichTet;
                int whichCorner;

                brickNodes[0] = (k * YDIM + j) * XDIM + i + 1;
                brickNodes[1] = brickNodes[0] + 1;
                brickNodes[2] = brickNodes[0] + XDIM + 1;
                brickNodes[3] = brickNodes[0] + XDIM;
                brickNodes[4] = brickNodes[0] + XDIM * YDIM;
                brickNodes[5] = brickNodes[1] + XDIM * YDIM;
                brickNodes[6] = brickNodes[2] + XDIM * YDIM;
                brickNodes[7] = brickNodes[3] + XDIM * YDIM;

                index = ((k * (YDIM - 1) + j) * (XDIM - 1) + i) * 5 * 4;
                for(whichTet = 0; whichTet < 5; ++whichTet)
                {
                    for(whichCorner = 0; whichCorner < 4; ++whichCorner)
                    {
                        connectivity[index++] = brickNodes[tetraCorners[evenOdd][whichTet][whichCorner]];
                    }
                }
                evenOdd = 1 - evenOdd;
            }
    /*
     * Write the zone header information.
     */
    i = TECZNE142((char*)"Simple Zone",
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
                  NULL,           /* ValueLocation = Nodal */
                  NULL,           /* SharVarFromZone */
                  &shrConn);
    /*
     * Write out the field data.
     */
    i = TECDAT142(&nNodes, x, &dIsDouble);
    i = TECDAT142(&nNodes, y, &dIsDouble);
    i = TECDAT142(&nNodes, z, &dIsDouble);
    i = TECDAT142(&nNodes, p, &dIsDouble);
    free(x);
    free(y);
    free(z);
    free(p);

    i = TECNODE142(&connectivityCount, connectivity);
    free(connectivity);

    i = TECEND142();

    return 0;
}
