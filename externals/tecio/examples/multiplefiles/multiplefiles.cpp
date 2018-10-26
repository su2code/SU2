/* This example illustrates working with TecFil to create multiple
 * plt files simultaneously.
 */
#if defined _MSC_VER
#pragma warning (disable: 4996) /* Windows strcpy warning off */
#endif

// Internal testing flags
// RUNFLAGS:none
// RUNFLAGS:--szl

/* DOCSTART:mulitplefiles.txt */
#include "TECIO.h"
#include "MASTER.h" /* for defintion of NULL */
#include <string.h>

int main(int argc, const char *argv[])
{
    /*
     * Open the file and write the tecplot datafile
     * header information
     */
    INTEGER4 Debug      = 1;
    INTEGER4 VIsDouble  = 0;
    INTEGER4 FileType   = 0;

    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    if (argc == 2 && strncmp(argv[1],"--szl",5) == 0)
        fileFormat = 1; 
    else
        fileFormat = 0; 

    INTEGER4 I          = 0;   /* Used to check the return value */

    I = TECINI142((char*)"SIMPLE DATASET", /* Name of the entire dataset.*/

                  (char*)"X1 Y1 P1",  /* Defines the variables for the data
                                       * file. Each zone must contain each of
                                       * the variables listed here. The order
                                       * of the variables in the list is used
                                       * to define the variable number (e.g.
                                       * X1 is Var 1).
                                       */
                  (char*)"multiplefiles-file1.plt",
                  (char*)".",      /* Scratch Directory */
                  &fileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);

    /* Set the parameters for TecZne */
    INTEGER4 ZoneType           = 0; /* sets the zone type to
                                    * ordered
                                    */
    INTEGER4 IMax               = 2; /* Create an IJ-ordered zone,
                                    * by using IMax and JMax
                                    * values that are greater
                                    * than one, and setting KMax
                                    * to one.
                                    */
    INTEGER4 JMax               = 2;
    INTEGER4 KMax               = 1;

    double   SolTime            = 0;
    INTEGER4 StrandID           = 0; /* StaticZone */
    INTEGER4 ParentZn           = 0; /* used for surface streams */

    INTEGER4 ICellMax           = 0; /* not used */
    INTEGER4 JCellMax           = 0; /* not used */
    INTEGER4 KCellMax           = 0; /* not used */

    INTEGER4 IsBlock            = 1; /* Block */

    INTEGER4 NFConns            = 0; /* this example does not use
                                    * face neighbors */
    INTEGER4 FNMode             = 0;
    INTEGER4 TotalNumFaceNodes  = 1;
    INTEGER4 TotalNumBndryFaces = 1;
    INTEGER4 TotalNumBndryConn  = 1;
    INTEGER4 ShrConn            = 0;


    /* Create an Ordered Zone */
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
                  &TotalNumBndryConn,
                  NULL,
                  NULL,
                  NULL,
                  &ShrConn);

    /* Set the variable values for the ordered zone. */
    float X1[4];
    float Y1[4];
    float P1[4];

    X1[0] = 0.125;
    Y1[0] = 0.5;
    P1[0] = 7.5;

    X1[1] = 0.625;
    Y1[1] = 0.5;
    P1[1] = 10.0;

    X1[2] = 0.125;
    Y1[2] = 0.875;
    P1[2] = 5.0;

    X1[3] = 0.625;
    Y1[3] = 0.875;
    P1[3] = 7.5;

    INTEGER4 DIsDouble   = 0;  /* set DIsDouble to 0, for float
                              * values.
                              */

    INTEGER4 III = IMax * JMax * KMax;
    I   = TECDAT142(&III, X1, &DIsDouble);
    I   = TECDAT142(&III, Y1, &DIsDouble);
    I   = TECDAT142(&III, P1, &DIsDouble);

    /* Open a new data file.  note: the first file is still open
     * because TecEnd was not called.
     */
    I = TECINI142((char*)"Auxiliary Data",
                  (char*)"X1 Y1 P1",
                  (char*)"multiplefiles-file2.plt",
                  (char*)".",
                  &fileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);

    /* Switch the active file to the newly created data file
     * (file2.plt) which is the second file opened with TECINI142
     * so we use 2.
     */
    INTEGER4 WhichFile = 2;
    I = TECFIL142(&WhichFile);

    /* Create a second zone, using many of the values from the first
     * zone, and write it to the second data file.
     */

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
                  &TotalNumBndryConn,
                  NULL,
                  NULL,
                  NULL,
                  &ShrConn);
    /* set the variable values for the second zone */
    float X2[4];
    float Y2[4];
    float P2[4];

    X2[0] = 0.375;
    Y2[0] = 0.125;
    P2[0] = 5;

    X2[1] = 0.875;
    Y2[1] = 0.125;
    P2[1] = 7.5;

    X2[2] = 0.375;
    Y2[2] = 0.5;
    P2[2] = 10;

    Y2[3] = 0.5;
    X2[3] = 0.875;
    P2[3] = 7.5;

    III = IMax * JMax * KMax;
    I   = TECDAT142(&III, X2, &DIsDouble);
    I   = TECDAT142(&III, Y2, &DIsDouble);
    I   = TECDAT142(&III, P2, &DIsDouble);

    /* Switch to the first file. */
    WhichFile = 1;
    I = TECFIL142(&WhichFile);

    /* Create an auxiliary data value and write it to the file */
    char DeformationValue[128];
    strcpy(DeformationValue, "0.98");

    I = TECAUXSTR142((char*)"DeformationValue",
                     DeformationValue);
    /* Close the first file */
    I = TECEND142();

    /* The remaining file will become the active file.  As such,
     * TecFil does not need to be called again to close the second
     * file.
     */
    I = TECEND142();

    return 0;
}

/* DOCEND */
