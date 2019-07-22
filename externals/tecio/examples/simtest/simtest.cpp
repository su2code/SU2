/*
 * Simple example c program to write a
 * binary datafile for tecplot.  This example
 * does the following:
 *
 *   1.  Open a datafile called "t.plt"
 *   2.  Assign values for X,Y, and P
 *   3.  Write out a zone dimensioned 4x5
 *   4.  Close the datafile.
 */

// Internal testing flags
// RUNFLAGS:none
// RUNFLAGS:--szl

#include "TECIO.h"
#include <string.h>

#ifndef NULL
#define NULL 0
#endif

enum FileType { FULL = 0, GRID = 1, SOLUTION = 2 };

int main(int argc, const char *argv[])
{
    float X[5][4], Y[5][4], P[5][4];
    double SolTime;
    INTEGER4 Debug, I, J, III, DIsDouble, VIsDouble, IMax, JMax, KMax, ZoneType, StrandID, ParentZn, IsBlock;
    INTEGER4 ICellMax, JCellMax, KCellMax, NFConns, FNMode, ShrConn, FileType;

    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    if (argc == 2 && strncmp(argv[1],"--szl",5) == 0)
        fileFormat = 1; 
    else
        fileFormat = 0; 


    Debug     = 1;
    VIsDouble = 0;
    DIsDouble = 0;
    IMax      = 4;
    JMax      = 5;
    KMax      = 1;
    ZoneType  = 0;      /* Ordered */
    SolTime   = 360.0;
    StrandID  = 0;     /* StaticZone */
    ParentZn  = 0;      /* No Parent */
    IsBlock   = 1;      /* Block */
    ICellMax  = 0;
    JCellMax  = 0;
    KCellMax  = 0;
    NFConns   = 0;
    FNMode    = 0;
    ShrConn   = 0;
    FileType  = FULL;

    /*
     * Open the file and write the tecplot datafile
     * header information
     */
    I = TECINI142((char*)"SIMPLE DATASET",
                  (char*)"X Y P",
                  (char*)"simtestcpp-t.plt",
                  (char*)".",
                  &fileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);

    for (J = 0; J < 5; J++)
        for (I = 0; I < 4; I++)
        {
            X[J][I] = (float)(I + 1);
            Y[J][I] = (float)(J + 1);
            P[J][I] = (float)((I + 1) * (J + 1));
        }
    /*
     * Write the zone header information.
     */
    I = TECZNE142((char*)"Simple Zone",
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
                  0,              /* TotalNumFaceNodes */
                  0,              /* NumConnectedBoundaryFaces */
                  0,              /* TotalNumBoundaryConnections */
                  NULL,           /* PassiveVarList */
                  NULL,           /* ValueLocation = Nodal */
                  NULL,           /* SharVarFromZone */
                  &ShrConn);
    /*
     * Write out the field data.
     */
    III = IMax * JMax;
    I   = TECDAT142(&III, &X[0][0], &DIsDouble);
    I   = TECDAT142(&III, &Y[0][0], &DIsDouble);
    I   = TECDAT142(&III, &P[0][0], &DIsDouble);

    I = TECEND142();

    return 0;
}
