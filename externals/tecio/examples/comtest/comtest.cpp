/*
 * Complex example C program to write a
 * binary data file for Tecplot. This example
 * does the following:
 *
 *   1.  Open a data file called "field.plt."
 *   2.  Open a data file called "line.plt."
 *   3.  Assign values for X, Y and P. These will be used
 *       in both the ordered and finite-element data files.
 *   4.  Write out an ordered zone dimensioned 4 x 5 to "field.plt."
 *   5.  Assign values for XL and YL arrays.
 *   6.  Write out data for line plot to "line.plt." Make the data
 *       use double precision.
 *   7.  Write out a finite-element zone to "field.plt."
 *   8.  Write out a text record to "field.plt."
 *   9.  Write out a geometry (circle) record to "field.plt."
 *  10.  Close file 1.
 *  11.  Close file 2.
 */

// Internal testing flags
// RUNFLAGS:none
// RUNFLAGS:--szl

#include <stdio.h>
#include <math.h>
#include "TECIO.h"
#include <string.h>

int main(int argc, const char *argv[])
{
    float    X[5][4], Y[5][4], P[5][4];
    double   XL[50], YL[50];
    double   SolTime;
    INTEGER4 Debug, I, J, K, L, III, NPts, NElm, DIsDouble, VIsDouble, IMax, JMax, KMax;
    INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, Clipping;
    INTEGER4 StrandID, ParentZn, FieldFileType, LineFileType;

    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    if (argc == 2 && strncmp(argv[1],"--szl",5) == 0)
        fileFormat = 1; 
    else
        fileFormat = 0; 

    INTEGER4 SharingZone[3] = {0, 0, 0};
    INTEGER4 IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
    INTEGER4 NM[12][4];
    double   XP, YP, ZP, FH, LineSpacing, PatternLength;
    double   BoxMargin, BoxLineThickness, TextAngle;
    INTEGER4 AttachToZone, Zone, Scope, PositionCoordSys, FontType, HeightUnits;
    INTEGER4 IsFilled, GeomType, LinePattern, NumEllipsePts;
    INTEGER4 Anchor, BoxType, BoxColor, BoxFillColor, TextColor, Color, FillColor;
    INTEGER4 ArrowheadStyle, ArrowheadAttachment, NumSegments, NumSegPts[1];
    double   LineThickness, ArrowheadSize, ArrowheadAngle;
    float    XGeomData[1], YGeomData[1], ZGeomData[1];
    enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };

    Debug         = 2;
    VIsDouble     = 0;
    DIsDouble     = 0;
    FieldFileType = FULL;
    LineFileType  = FULL;
    /*
     * Open order.plt and write the header information.
     */
    I = TECINI142((char*)"DATASET WITH ONE ORDERED ZONE AND ONE FE-QUAD ZONE OVER 2 TIME STEPS",
                  (char*)"X Y P",
                  (char*)"comtestcpp-field.plt",
                  (char*)".",
                  &fileFormat,
                  &FieldFileType,
                  &Debug,
                  &VIsDouble);
    /*
     * Open line.plt and write the header information.
     */
    VIsDouble = 1;
    I = TECINI142((char*)"DATASET WITH ONE I-ORDERED ZONE",
                  (char*)"X Y",
                  (char*)"comtestcpp-line.plt",
                  (char*)".",
                  &fileFormat,
                  &LineFileType,
                  &Debug,
                  &VIsDouble);

    /*
     * Calculate values for the field variables.
     */
    for (J = 0; J < 5; J++)
        for (I = 0; I < 4; I++)
        {
            X[J][I] = (float)(I + 1);
            Y[J][I] = (float)(J + 1);
            P[J][I] = (float)((I + 1) * (J + 1));
        }

    /*
     * Make sure writing to file #1.
     */
    III = 1;
    I = TECFIL142(&III);

    /*
     * Write the zone header information for the ordered zone.
     */
    IMax      = 4;
    JMax      = 5;
    KMax      = 1;
    ICellMax  = 0;
    JCellMax  = 0;
    KCellMax  = 0;
    ZoneType  = 0;
    SolTime   = 10.0;
    StrandID  = 1;
    ParentZn  = 0;
    IsBlock   = 1;
    NumFaceConnections = 0;
    FaceNeighborMode   = 0;
    ShareConnectivityFromZone = 0;
    I = TECZNE142((char*)"Ordered Zone 1",
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
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  NULL,      /* ShareVarFromZone */
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  &ShareConnectivityFromZone);
    /*
     * Write out the field data for the ordered zone.
     */
    III = IMax * JMax;
    I   = TECDAT142(&III, &X[0][0], &DIsDouble);
    I   = TECDAT142(&III, &Y[0][0], &DIsDouble);
    I   = TECDAT142(&III, &P[0][0], &DIsDouble);

    /*
     * Calculate values for the I-ordered zone.
     */

    for (I = 0; I < 50; I++)
    {
        XL[I] = I + 1;
        YL[I] = sin((double)(I + 1) / 20.0);
    }
    /*
     * Switch to the "line.plt" file (file number 2)
     * and write out the line plot data.
     */

    III = 2;
    I = TECFIL142(&III);

    /*
     * Write the zone header information for the XY-data.
     */
    IMax = 50;
    JMax = 1;
    KMax = 1;
    SolTime = 0.0;
    StrandID = 0; /* StaticZone */
    I = TECZNE142((char*)"XY Line plot",
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
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  NULL,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
    /*
     * Write out the line plot.
     */
    DIsDouble = 1;
    III = IMax;
    I   = TECDAT142(&III, (float *) & XL[0], &DIsDouble);
    I   = TECDAT142(&III, (float *) & YL[0], &DIsDouble);

    /*
     * Switch back to the field plot file and write out
     * the finite-element zone.
     */
    III = 1;
    I = TECFIL142(&III);

    /*
     * Move the coordinates so this zone's not on top of the other
     */
    for (J = 0; J < 5; J++)
        for (I = 0; I < 4; I++)
        {
            X[J][I] = (float)(I + 6);
            Y[J][I] = (float)(J + 1);
            P[J][I] = (float)((I + 1) * (J + 1));
        }

    /*
     * Write the zone header information for the finite-element zone.
     */
    ZoneType  = 3;  /* FEQuad */
    NPts      = 20; /* Number of points */
    NElm      = 12; /* Number of elements */
    KMax      = 0;  /* Unused */
    SolTime   = 10.0;
    StrandID  = 2;
    I = TECZNE142((char*)"Finite Zone 1",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &SolTime,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  NULL,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
    /*
     * Write out the field data for the finite-element zone.
     */
    IMax      = 4;
    JMax      = 5;
    III       = IMax * JMax;
    DIsDouble = 0;
    I    = TECDAT142(&III, &X[0][0], &DIsDouble);
    I    = TECDAT142(&III, &Y[0][0], &DIsDouble);
    I    = TECDAT142(&III, &P[0][0], &DIsDouble);

    /*
     * Calculate and then write out the connectivity list.
     * Note: The NM array references cells starting with
     *       offset of 1.
     */

    for (I = 1; I < IMax; I++)
        for (J = 1; J < JMax; J++)
        {
            K = I + (J - 1) * (IMax - 1);
            L = I + (J - 1) * IMax;
            NM[K-1][0] = L;
            NM[K-1][1] = L + 1;
            NM[K-1][2] = L + IMax + 1;
            NM[K-1][3] = L + IMax;
        }

    I = TECNOD142((INTEGER4 *)NM);

    /*
     * Calculate values for the new solution variable.
     */
    for (J = 0; J < 5; J++)
        for (I = 0; I < 4; I++)
        {
            P[J][I] = (float)(2.0 * (I + 1) * (J + 1));
        }

    /*
     * Write the zone header information for time step 2
     */
    ZoneType = 0;
    IMax      = 4;
    JMax      = 5;
    KMax      = 1;
    SolTime = 20.0;
    StrandID = 1;
    SharingZone[0] = 1;
    SharingZone[1] = 1;
    SharingZone[2] = 0; /* solution variable is not shared */
    ShareConnectivityFromZone = 0;

    I = TECZNE142((char*)"Ordered Zone 2",
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
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,
                  NULL,
                  SharingZone,
                  &ShareConnectivityFromZone);

    /*
     * Write out the solution variable the grid variables are shared.
     */
    IMax      = 4;
    JMax      = 5;
    III       = IMax * JMax;
    DIsDouble = 0;
    I    = TECDAT142(&III, &P[0][0], &DIsDouble);

    /*
     * Calculate values for the new solution variable.
     */
    for (J = 0; J < 5; J++)
        for (I = 0; I < 4; I++)
        {
            P[J][I] = (float)(3.0 * (I + 1) * (J + 1));
        }

    /*
     * Write another time step for the FEZone and share from the first
     */
    ZoneType = 3;
    SolTime = 20.0;
    StrandID = 2;
    KMax = 0;
    SharingZone[0] = 2;
    SharingZone[1] = 2;
    SharingZone[2] = 0; /* solution variable is not shared */
    ShareConnectivityFromZone = 2;
    I = TECZNE142((char*)"Finite Zone 2",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &SolTime,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  SharingZone,
                  &ShareConnectivityFromZone);

    /*
     * Write out the solution variable the grid variables are shared.
     */
    IMax      = 4;
    JMax      = 5;
    III       = IMax * JMax;
    DIsDouble = 0;
    I    = TECDAT142(&III, &P[0][0], &DIsDouble);

    /*
     * Prepare to write out text record. Text is positioned
     * at 0.5, 0.5 in frame units and has a height
     * of 0.05 frame units.
     */
    XP               = 50.0;
    YP               = 50.0;
    ZP               = 0.0;
    FH               = 5.0;
    Scope            = 1; /* Local */
    Clipping         = 1; /* Clip to frame */
    PositionCoordSys = 1; /* Frame */
    FontType         = 1; /* Helv Bold */
    HeightUnits      = 1; /* Frame */
    AttachToZone     = 0;
    Zone             = 0;
    BoxType          = 0; /* None */
    BoxMargin        = 5.0;
    BoxLineThickness = 0.5;
    BoxColor         = 3;
    BoxFillColor     = 7;
    TextAngle        = 0.0;
    Anchor           = 0; /* Left */
    LineSpacing      = 1.0;
    TextColor        = 0; /* Black */

    III =  TECTXT142(&XP,
                     &YP,
                     &ZP,
                     &PositionCoordSys,
                     &AttachToZone,
                     &Zone,
                     &FontType,
                     &HeightUnits,
                     &FH,
                     &BoxType,
                     &BoxMargin,
                     &BoxLineThickness,
                     &BoxColor,
                     &BoxFillColor,
                     &TextAngle,
                     &Anchor,
                     &LineSpacing,
                     &TextColor,
                     &Scope,
                     &Clipping,
                     (char*)"Hi Mom",
                     (char*)"");

    /*
     * Prepare to write out geometry record (circle). Circle is
     * positioned at 25, 25 (in frame units) and has a radius of
     * 20 percent. Circle is drawn using a dashed line.
     */


    XP                  = 25.0;
    YP                  = 25.0;
    ZP                  = 0.0;
    IsFilled            = 0;
    Color               = 0;
    FillColor           = 7;
    GeomType            = 3; /* Circle */
    LinePattern         = 1; /* Dashed */
    LineThickness       = 0.3;
    PatternLength       = 1.5;
    NumEllipsePts       = 72;
    ArrowheadStyle      = 0;
    ArrowheadAttachment = 0;
    ArrowheadSize       = 0.1;
    ArrowheadAngle      = 15.0;
    NumSegments         = 1;
    NumSegPts[0]        = 1;

    XGeomData[0] = 20.0;
    YGeomData[0] = 0.0;
    ZGeomData[0] = 0.0;


    III =  TECGEO142(&XP,
                     &YP,
                     &ZP,
                     &PositionCoordSys,
                     &AttachToZone,
                     &Zone,
                     &Color,
                     &FillColor,
                     &IsFilled,
                     &GeomType,
                     &LinePattern,
                     &PatternLength,
                     &LineThickness,
                     &NumEllipsePts,
                     &ArrowheadStyle,
                     &ArrowheadAttachment,
                     &ArrowheadSize,
                     &ArrowheadAngle,
                     &Scope,
                     &Clipping,
                     &NumSegments,
                     NumSegPts,
                     &XGeomData[0],
                     &YGeomData[0],
                     &ZGeomData[0],
                     (char*)"");

    /*
     * Close out file 1.
     */
    I = TECEND142();

    /*
     * Close out file 2.
     */
    III = 2;
    I = TECFIL142(&III);
    I = TECEND142();

    return 0;
}
