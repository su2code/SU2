/* This example demonstrates adding a text object to a Tecplot
 * data file.
 */
#if defined _MSC_VER
#pragma warning (disable: 4996) /* Windows strcpy warning off */
#endif

// Internal testing flags
// RUNFLAGS:none
// RUNFLAGS:--szl

/* DOCSTART:tectxt.txt*/
#include "TECIO.h"
#include <string.h>

int main(int argc, const char *argv[])
{
    /* Open the file & write the datafile header information */
    INTEGER4 Debug     = 1;
    INTEGER4 VIsDouble = 0;

    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    if (argc == 2 && strncmp(argv[1],"--szl",5) == 0)
        fileFormat = 1; 
    else
        fileFormat = 0; 

    INTEGER4 FileType  = 0;
    INTEGER4 I         = 0;  /* used to check the return value */

    I = TECINI142((char*)"Text",
                  (char*)"X Y P",
                  (char*)"text.plt",
                  (char*)".",
                  &fileFormat,
                  &FileType,
                  &Debug,
                  &VIsDouble);

    /* Specify the X, Y and Z position of the anchor point */
    double   XPos              = 0.0;
    double   YPos              = 1.0;
    double   ZPos              = 0.0; /* N/A for 2D text */

    INTEGER4 PosCoordMode      = 0;   /* use grid coordinates */

    /* opt not to attach the text to a given zone.  When text is
     * attached to a given zone, it is displayed only when the zone
     * is displayed.
     */
    INTEGER4 AttachToZone      = 0;
    INTEGER4 Zone              = 2;


    /* Specify the font values */
    INTEGER4 Font              = 1;  /* Helvetica Bold */
    INTEGER4 FontHeightUnits   = 2;  /* in grid coordinates */
    double   FontHeight        = 18;

    /* Set the box style parameters */
    INTEGER4 BoxType           = 1;     /* filled box */
    double   BoxMargin         = .5;    /* margin between the text
                                       * and the text box
                                       */
    double   BoxLineThickness  = .1;
    INTEGER4 BoxColor          = 0;     /* set the box line color
                                       * to black.
                                       */
    INTEGER4 BoxFillColor      = 1;     /* set the box fill color
                                       * to red.
                                       */

    /* set the font properties */
    double   Angle             = 30;    /* angle of the text */
    INTEGER4 Anchor            = 1;     /* set the anchor point to
                                       * the center of the text
                                       * box.
                                       */
    double   LineSpacing       = 1.5;
    INTEGER4 TextColor         = 7;     /* set the font color to
                                       * white
                                       */

    INTEGER4 Scope             = 1;     /* set the text to "local",
                                       * i.e. available in the
                                       * current frame only.
                                       */
    INTEGER4 Clipping          = 1;


    char     Text[60];
    char     MFC[24];
    strcpy(Text, "Sample Text");
    strcpy(MFC, "My Macro");

    I = TECTXT142(&XPos,
                  &YPos,
                  &ZPos,
                  &PosCoordMode,
                  &AttachToZone,
                  &Zone,
                  &Font,
                  &FontHeightUnits,
                  &FontHeight,
                  &BoxType,
                  &BoxMargin,
                  &BoxLineThickness,
                  &BoxColor,
                  &BoxFillColor,
                  &Angle,
                  &Anchor,
                  &LineSpacing,
                  &TextColor,
                  &Scope,
                  &Clipping,
                  Text,
                  MFC);

    I = TECEND142();

    return 0;
}

/* DOCEND */
