#include "tkogl.h"

#if defined(__WIN32__) || defined(_WIN32)
#   if defined(__BORLANDC__) || defined(__CYGWIN__)
#      define strcasecmp(a,b) stricmp(a,b)
#   else
#      define strcasecmp(a,b) _stricmp(a,b)
#   endif
#endif

typedef struct {
   unsigned char x, y;
} CPoint;

typedef struct {
   int nPoint;
   CPoint * point;
} Stroke;

typedef struct {
   CPoint box;
   CPoint org;
   CPoint inc;
   int nStroke;
   Stroke * stroke;
} CDesc;

typedef CDesc * StrokeFont [256];

typedef StrokeFont* StrokeFontPtr;

/*----------------------------------------------------------------------
 *
 *  NewFont
 *  =======
 *
 *  Returns a pointer to a newly allocated empty font 
 */

static StrokeFontPtr NewFont () 
{
   StrokeFontPtr fontPtr = (StrokeFontPtr) malloc (sizeof (StrokeFont));
   int glyph;
   for (glyph = 0; glyph < 256; glyph++) {
      (*fontPtr) [glyph] = NULL;
   }
   return fontPtr;
}

/*----------------------------------------------------------------------
 *
 *  FreeFont
 *  =======
 *
 *  Given a pointer to a font, deallocates it.
 */

static void FreeFont (StrokeFontPtr fontPtr)
{
   int glyph;

   for (glyph = 0; glyph < 256; glyph++) {
      if ((*fontPtr) [glyph] != NULL) {
	 Stroke * stroke = ((*fontPtr) [glyph]) -> stroke;
	 int istroke = ((*fontPtr) [glyph]) -> nStroke;
	 while (istroke--) {
	    free (stroke->point);
	    stroke++;
	 }
	 free (((*fontPtr) [glyph]) -> stroke);
      }
   }
   free (fontPtr);
}

/*----------------------------------------------------------------------
 *
 *  LoadFont
 *  ========
 *
 *  Given the name of a text file containing a description of a stroke font,
 *  loads it into the supplied (empty) font. Returns a standard Tcl
 *  result indicating the success of the operation.
 */

static int LoadFont (Tcl_Interp * interp, char * filename, StrokeFont font) 
{
   int glyph, x, y;
   CDesc *desc;
   int nstroke;
   int npoint;
   Stroke strokebuf [50];
   CPoint pointbuf [1000];
   char buf [80];
   FILE *f;
   int ch;

#  define ERROR(s1,s2) { \
      Tcl_AppendResult (interp, s1, s2, (char*) NULL); \
      return TCL_ERROR; \
   }

#  define STORESTROKES {\
      desc->nStroke = nstroke; \
      desc->stroke = (Stroke*) malloc (sizeof (Stroke) * nstroke); \
      memcpy (desc->stroke, strokebuf, sizeof (Stroke) * nstroke); \
   }

   if ((f = fopen (filename, "r")) == 0) ERROR("Can't open ", filename);

   for (glyph = 0; glyph < 256; glyph++) {
      font [glyph] = NULL;
   }

   nstroke = 0;
   glyph = -1;

   while (!feof (f)) {
      if (fscanf (f, "%79s", buf) != 1) 
	 ERROR("Error reading ", filename);
      if (strcasecmp (buf, "glyph") == 0) {
	 if (glyph != -1) font [glyph] = desc;
	 if (nstroke > 0) STORESTROKES;
	 if (fscanf (f, "%d", &glyph)!= 1) 
	    ERROR ("Error in glyph ", filename);
	 desc = (CDesc*) malloc (sizeof (CDesc));
	 desc->box.x = desc->box.y = 0;
	 desc->org.x = desc->org.y = 0;
	 desc->inc.x = desc->inc.y = 0;
	 desc->nStroke = nstroke = 0;
	 desc->stroke = NULL;
      }
      else if (strcasecmp (buf, "blackbox") == 0) {
	 if (fscanf (f, "%d %d", &x, &y) != 2) 
	    ERROR ("Error in blackbox", filename);
	 desc->box.x = (unsigned char) x;
	 desc->box.y = (unsigned char) y;
      }
      else if (strcasecmp (buf, "origin") == 0) {
	 if (fscanf (f, "%d %d", &x, &y) != 2) 
	    ERROR ("Error in origin", filename);
	 desc->org.x = (unsigned char) x;
	 desc->org.y = (unsigned char) y;
      }
      else if (strcasecmp (buf, "cellinc") == 0) {
	 if (fscanf (f, "%d %d", &x, &y) != 2) 
	    ERROR ("Error in cellinc", filename);
	 desc->inc.x = (unsigned char) x;
	 desc->inc.y = (unsigned char) y;
      }
      else if (strcasecmp (buf, "stroke") == 0) {
	 npoint = 0;
	 printf ("stroke\n");
	 for (;;) {
	    do { 
	       ch = getc (f);
	    } while (!feof (f) && isspace (ch));
	    ungetc (ch, f);
	    if (isalpha (ch)) break;
	    if (fscanf (f, "%d %d", &x, &y) != 2) 
	       ERROR ("Error reading stroke in ", filename);
	    pointbuf [npoint].x = (unsigned char) x;
	    pointbuf [npoint].y = (unsigned char) y;
	    npoint++;
	    printf ("%d %d\n", x, y);
	 }
	 strokebuf [nstroke].nPoint = npoint;
	 strokebuf [nstroke].point = (CPoint*) malloc (sizeof (CPoint)*npoint);
	 memcpy (strokebuf [nstroke].point, pointbuf, sizeof (CPoint)*npoint);
	 nstroke++;
      }
   }
   
   if (glyph != -1) font [glyph] = desc;
   if (nstroke > 0) STORESTROKES;

   fclose (f);

   return TCL_OK;
}


/*----------------------------------------------------------------------
 *
 *  StrokeFontExt
 *  =============
 *
 *  Main Tkogl extension command for handling stroked fonts.
 */

int StrokeFontExt (Tcl_Interp *interp, int argc, char* argv [])
{   
   int iarg;
   static StrokeFontPtr fontPtr = NULL;
   
   for (iarg = 2; iarg < argc; iarg++) {	
      int len = strlen (argv [iarg]);
      if (strncmp (argv [iarg], "-load", len) == 0) {
	 iarg++;
	 if (iarg == argc) {
            Tcl_AppendResult (interp, "not enough arguments", (char*)NULL);
            return TCL_ERROR;
         }
	 if (fontPtr != NULL) FreeFont (fontPtr);
	 fontPtr = NewFont ();
	 return LoadFont (interp, argv [iarg], *fontPtr);
      }
      else {
	 Tcl_AppendResult (interp, "unknown font option: ", argv [iarg]);
	 return TCL_ERROR;
      }
   }
}


