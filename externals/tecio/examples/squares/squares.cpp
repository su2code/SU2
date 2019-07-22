/* This example creates a group of square geometries, each with a
 * different fill color */
#if defined _MSC_VER
#pragma warning (disable: 4996) /* Windows strcpy warning off */
#endif

// Internal testing flags
// RUNFLAGS:none
// RUNFLAGS:--szl

/* DOCSTART:tecgeo.txt*/
#include "TECIO.h"
#include <string.h>

int main(int argc, const char *argv[])
{
    INTEGER4 Debug      = 1;
    INTEGER4 VIsDouble  = 0;

    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    if (argc == 2 && strncmp(argv[1],"--szl",5) == 0)
        fileFormat = 1; 
    else
        fileFormat = 0; 

    INTEGER4 FileType   = 0;
    INTEGER4 I          = 0; /* use to check return values */


    /* Open the file and write the tecplot datafile
     * header information
     */
    I = TECINI142((char*)"Square Geometries",
                  (char*)"X Y P",
                  (char*)"squares.plt",
                  (char*)".",
                  &fileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);

    double ZPos = 0.0;  /* N/A for squares */
    double XPos;
    double YPos;

    INTEGER4 PosCoordMode  =  0; /* use grid coordinates */

    /* opt not to attach the text to a given zone.  When text is
     * attached to a given zone, it is displayed only when the zone
     * is displayed.
     */
    INTEGER4 AttachToZone  =  0;
    INTEGER4 Zone          =  1;

    /* Set the Geometry Style Values */
    INTEGER4 Color         =  0;   /* set the outline color to
                                  * black
                                  */
    INTEGER4 IsFilled      =  1;
    INTEGER4 GeomType      =  2;   /* set the geometry type to
                                  * square
                                  */
    INTEGER4 LinePattern   =  5;   /* set the line pattern to
                                  * DashDotDot
                                  */
    double   PatternLength =  .1;
    double   LineThick     =  .2;

    /* N/A for square geometries */
    INTEGER4 NumPts        = 100;
    INTEGER4 ArrowStyle    = 1;
    INTEGER4 ArrowAttach   = 0;
    double   ArrowSize     = 1;
    double   ArrowAngle    = 30;
    INTEGER4 NumSegments   = 15;
    INTEGER4 NumSegPts     = 25;


    INTEGER4 Scope         =  1;  /* set the text to "local", i.e.
                                 * available in the current frame
                                 * only.
                                 */
    INTEGER4 Clipping      =  1;

    /* Specify the length of a side of the square.  The units used
     * are those defined with PosCoordMode.
     */
    float XGeomData     = 2.5;

    float YGeomData     = 0;  /* N/A for square geometries */
    float ZGeomData     = 0;  /* N/A for square geometries */

    char * MFC = new char[128];
    strcpy(MFC, "SQUARE");

    for (INTEGER4 ii = 0; ii <= 7; ii++)
    {
        INTEGER4 FillColor     =  ii;
        XPos                   = (double) ii;
        YPos                   = (double) ii;

        I = TECGEO142(&XPos,
                      &YPos,
                      &ZPos,
                      &PosCoordMode,
                      &AttachToZone,
                      &Zone,
                      &Color,
                      &FillColor,
                      &IsFilled,
                      &GeomType,
                      &LinePattern,
                      &PatternLength,
                      &LineThick,
                      &NumPts,
                      &ArrowStyle,
                      &ArrowAttach,
                      &ArrowSize,
                      &ArrowAngle,
                      &Scope,
                      &Clipping,
                      &NumSegments,
                      &NumSegPts,
                      &XGeomData,
                      &YGeomData,
                      &ZGeomData,
                      MFC);
    }

    I = TECEND142();

    delete MFC;
    return 0;
}

/* DOCEND */
