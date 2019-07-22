/* This example creates a simple set of IJ-ordered zones */
/* DOCSTART:ij_ordered.txt*/


// Internal testing flags
// RUNFLAGS:none
// RUNFLAGS:--szl


#include "TECIO.h"
#include "MASTER.h" /* for defintion of NULL */
#include <string.h>

int main(int argc, const char *argv[])
{
    INTEGER4 Debug      = 1;
    INTEGER4 VIsDouble  = 0;
    INTEGER4 FileType   = 0;
    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    if (argc == 2 && strncmp(argv[1],"--szl",5) == 0)
        fileFormat = 1; 
    else
        fileFormat = 0; 

    INTEGER4 I          = 0; /* Used to track return codes */

    /*
     * Open the file and write the tecplot datafile
     * header information
     */
    I = TECINI142((char*)"IJ Ordered Zones", /* Name of the entire
                                              * dataset.
                                              */
                  (char*)"X Y P",  /* Defines the variables for the data
                                    * file. Each zone must contain each of
                                    * the variables listed here. The order
                                    * of the variables in the list is used
                                    * to define the variable number (e.g.
                                    * X is Var 1).
                                    */
                  (char*)"ij_ordered.plt",
                  (char*)".",      /* Scratch Directory */
                  &fileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);

    float X1[4];
    float Y1[4];
    float P1[4];
    float X2[4];
    float Y2[4];
    float P2[4];

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
    INTEGER4 ShrConn                  = 0;

    /*Ordered Zone Parameters*/
    INTEGER4 IMax = 2;
    INTEGER4 JMax = 2;
    INTEGER4 KMax = 1;

    X1[0] = .125;
    Y1[0] = .5;
    P1[0] = 5;

    X1[1] = .625;
    Y1[1] = .5;
    P1[1] = 7.5;

    X1[2] = .125;
    Y1[2] = .875;
    P1[2] = 10;

    X1[3] = .625;
    Y1[3] = .875;
    P1[3] = 7.5;

    X2[0] = .375;
    Y2[0] = .125;
    P2[0] = 5;

    X2[1] = .875;
    Y2[1] = .125;
    P2[1] = 7.5;

    X2[2] = .375;
    Y2[2] = .5;
    P2[2] = 10;

    X2[3] = .875;
    Y2[3] = .5;
    P2[3] = 7.5;

    /*  Ordered Zone */
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
                  NULL,
                  NULL,
                  NULL,
                  &ShrConn);
    INTEGER4 III = IMax * JMax * KMax;
    I   = TECDAT142(&III, X1, &DIsDouble);
    I   = TECDAT142(&III, Y1, &DIsDouble);
    I   = TECDAT142(&III, P1, &DIsDouble);

    I = TECZNE142((char*)"Ordered Zone2",
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
                  NULL,
                  NULL,
                  NULL,
                  &ShrConn);

    I   = TECDAT142(&III, X2, &DIsDouble);
    I   = TECDAT142(&III, Y2, &DIsDouble);
    I   = TECDAT142(&III, P2, &DIsDouble);

    I = TECEND142();
    return 0;
}
/* DOCEND */
