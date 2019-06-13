/*
 * Simple example c program to write a
 * binary datafile for tecplot.  This example
 * does the following:
 *
 *   1.  Open a datafile called "t.plt"
 *   2.  Assign values for x,y, and p.
 *   3.  Write out a hexahedral (brick) zone.
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

/* The zone will appear the same as an ordered zone with these dimensions: */
#define XDIM 5
#define YDIM 4
#define ZDIM 3

enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };

int main(int argc, const char *argv[])
{
    float *x, *y, *z, *p, *cc;
    int *connectivity;
    double solTime;
    INTEGER4 debug, i, j, k, dIsDouble, vIsDouble, zoneType, strandID, parentZn, isBlock;
    INTEGER4 iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType;
    INTEGER4 nNodes, nCells, nFaces, connectivityCount, index;
    int valueLocation[] = {1,1,1,1,0}; 

    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    if (argc == 2 && strncmp(argv[1],"--szl",5) == 0)
        fileFormat = 1; 
    else
        fileFormat = 0; 

    debug     = 1;
    vIsDouble = 0;
    dIsDouble = 0;
    nNodes = XDIM * YDIM * ZDIM;
    nCells = (XDIM - 1) * (YDIM - 1) * (ZDIM - 1);
    nFaces = 6; /* Not used */
    zoneType  = 5;      /* Brick */
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
                  (char*)"x y z p cc",  // NOTE: Make sure and change valueLocation above if this changes.
                  (char*)"brick.plt",
                  (char*)".",
                  &fileFormat,
                  &fileType,
                  &debug,
                  &vIsDouble);

    x  = (float*)malloc(nNodes * sizeof(float));
    y  = (float*)malloc(nNodes * sizeof(float));
    z  = (float*)malloc(nNodes * sizeof(float));
    p  = (float*)malloc(nNodes * sizeof(float));
    cc = (float*)malloc(nCells * sizeof(float));
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

    for (i = 0; i < nCells; ++i)
        cc[i] = (float)i;

    connectivityCount = 8 * nCells;
    connectivity = (INTEGER4*)malloc(connectivityCount * sizeof(INTEGER4));
    for (k = 0; k < ZDIM - 1; k++)
        for (j = 0; j < YDIM - 1; j++)
            for (i = 0; i < XDIM - 1; i++)
            {
                index = ((k * (YDIM - 1) + j) * (XDIM - 1) + i) * 8;
                connectivity[index] = (k * YDIM + j) * XDIM + i + 1;
                connectivity[index + 1] = connectivity[index] + 1;
                connectivity[index + 2] = connectivity[index] + XDIM + 1;
                connectivity[index + 3] = connectivity[index] + XDIM;
                connectivity[index + 4] = connectivity[index] + XDIM * YDIM;
                connectivity[index + 5] = connectivity[index + 1] + XDIM * YDIM;
                connectivity[index + 6] = connectivity[index + 2] + XDIM * YDIM;
                connectivity[index + 7] = connectivity[index + 3] + XDIM * YDIM;
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
                  valueLocation,  /* ValueLocation = Nodal */
                  NULL,           /* SharVarFromZone */
                  &shrConn);
    /*
     * Write out the field data.
     */
    i = TECDAT142(&nNodes, x, &dIsDouble);
    i = TECDAT142(&nNodes, y, &dIsDouble);
    i = TECDAT142(&nNodes, z, &dIsDouble);
    i = TECDAT142(&nNodes, p, &dIsDouble);
    i = TECDAT142(&nCells, cc, &dIsDouble);
    free(x);
    free(y);
    free(z);
    free(p);
    free(cc);

    i = TECNODE142(&connectivityCount, connectivity);
    free(connectivity);

    i = TECEND142();

    return 0;
}
