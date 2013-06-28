#include <tk.h>
#if defined(__WIN32__) || defined(_WIN32)
#   define WIN32_LEAN_AND_MEAN
#   include <windows.h>
#   undef WIN32_LEAN_AND_MEAN
#endif
#include <GL/gl.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "load3ds.h"
#include "tkogl.h"

/*---------------------------------------------------------------------------*
 *
 *  The following are data structures corresponding to the 
 *  chunk/chunktype overall format of 3DStudio files
 *
 *---------------------------------------------------------------------------*/

typedef unsigned char  byte;
typedef unsigned short word;
typedef unsigned long  dword;

typedef struct {
    word    id;
    dword   len;
} TChunkHeader, *PChunkHeader;

enum {
    CHUNK_RGBF      = 0x0010,
    CHUNK_RGBB      = 0x0011,
    CHUNK_PERCI	    = 0x0030,
    CHUNK_PERCF     = 0x0031,
    CHUNK_MAIN      = 0x4D4D,
        CHUNK_OBJMESH   = 0x3D3D,
            CHUNK_BKGCOLOR  = 0x1200,
            CHUNK_AMBCOLOR  = 0x2100,
            CHUNK_OBJBLOCK  = 0x4000,
                CHUNK_TRIMESH   = 0x4100,
                    CHUNK_VERTLIST  = 0x4110,
                    CHUNK_FACELIST  = 0x4120,
                    CHUNK_FACEMAT   = 0x4130,
                    CHUNK_MAPLIST   = 0x4140,
                    CHUNK_SMOOLIST  = 0x4150,
                    CHUNK_TRMATRIX  = 0x4160,
                CHUNK_LIGHT     = 0x4600,
                    CHUNK_SPOTLIGHT = 0x4610,
                CHUNK_CAMERA    = 0x4700,
        CHUNK_MATERIAL  = 0xAFFF,
            CHUNK_MATNAME   = 0xA000,
            CHUNK_AMBIENT   = 0xA010,
            CHUNK_DIFFUSE   = 0xA020,
            CHUNK_SPECULAR  = 0xA030,
            CHUNK_SHININESS = 0xA040,
            CHUNK_SHININESS1= 0xA041,
            CHUNK_SHININESS2= 0xA042,
            CHUNK_TEXTURE   = 0xA200,
            CHUNK_BUMPMAP   = 0xA230,
            CHUNK_MAPFILE   = 0xA300,
        CHUNK_KEYFRAMER = 0xB000,
            CHUNK_FRAMES    = 0xB008,

};

/*---------------------------------------------------------------------------
 *
 * The following are forward declarations for
 * procedures for reading each chunk type
 *
 *---------------------------------------------------------------------------*/

typedef void ChunkReadProc (Tcl_DString* desc, FILE *f, long p);

static ChunkReadProc  
   ChunkReader,
   TriMeshReader,
   RGBFReader,
   RGBBReader,   
   PercentIReader,
   PercentFReader,
   ObjBlockReader,
   VertListReader,
   FaceListReader,
   FaceMatReader, 
   MapListReader, 
   SmooListReader,
   TrMatrixReader,
   LightReader,   
   SpotLightReader,
   CameraReader,  
   MatNameReader, 
   MapFileReader, 
   FramesReader;

/*---------------------------------------------------------------------------
 *
 * The following is a table that associates each chunk type with the proper
 * chunk reader function and a string to be used in the property list
 * to be returned in a Tcl result
 *
 *---------------------------------------------------------------------------*/

struct {
    word id;
    char *name;
    ChunkReadProc* func;
} ChunkNames[] = {
    {CHUNK_RGBF,        "rgb",              RGBFReader},
    {CHUNK_RGBB,        "rgb",              RGBBReader},
    {CHUNK_PERCI, 	"percent",          PercentIReader},
    {CHUNK_PERCF,       "percent",          PercentFReader},
    {CHUNK_MAIN,        "main",             NULL},
    {CHUNK_OBJMESH,     "object-mesh",      NULL},
    {CHUNK_BKGCOLOR,    "background",       NULL},
    {CHUNK_AMBCOLOR,    "ambient",          NULL},
    {CHUNK_OBJBLOCK,    "object-block",     ObjBlockReader},
    {CHUNK_TRIMESH,     "tri-mesh",         TriMeshReader},
    {CHUNK_VERTLIST,    "vertex-list",      VertListReader},
    {CHUNK_FACELIST,    "face-list",        FaceListReader},
    {CHUNK_FACEMAT,     "face-material",    FaceMatReader},
    {CHUNK_MAPLIST,     "mappings-list",    MapListReader},
    {CHUNK_SMOOLIST,    "smoothings",       SmooListReader},
    {CHUNK_TRMATRIX,    "matrix",           TrMatrixReader},
    {CHUNK_LIGHT,       "light",            LightReader},
    {CHUNK_SPOTLIGHT,   "spotlight",        SpotLightReader},
    {CHUNK_CAMERA,      "camera",           CameraReader},

    {CHUNK_MATERIAL,    "material",         NULL},
    {CHUNK_MATNAME,     "name",   	    MatNameReader},
    {CHUNK_AMBIENT,     "ambient",          NULL},
    {CHUNK_DIFFUSE,     "diffuse",          NULL},
    {CHUNK_SPECULAR,    "specular",         NULL},
    {CHUNK_SHININESS,   "shininess",        NULL},
    {CHUNK_SHININESS1,  "shininess1",       NULL},
    {CHUNK_SHININESS2,  "shininess2",       NULL},
    
    {CHUNK_TEXTURE,     "texture-map",      NULL},
    {CHUNK_BUMPMAP,     "bump-map",         NULL},
    {CHUNK_MAPFILE,     "name",     	    MapFileReader},

    {CHUNK_KEYFRAMER,   "keyframer-data",   NULL},
    {CHUNK_FRAMES,      "frames",           FramesReader},  

};

/*---------------------------------------------------------------------------
 *
 * The following are data structures for storing the geometry of a
 * Tri-Mesh object
 *
 *---------------------------------------------------------------------------*/

typedef GLfloat Vector [3];

typedef struct {
   int ivtx [3];
   Vector normal;
   dword group;
   word flags;
} FaceInfo;
   
typedef struct SmoothInfoNode {
   Vector normal;
   dword group;
   struct SmoothInfoNode *next;
} *SmoothList;

typedef struct {
   Vector coord;
   char texFlag;
   GLfloat texCoord [2];
   SmoothList smoothList;
} VertexInfo;

typedef struct FaceListNode {
   int idx;
   struct FaceListNode* next;
} *FaceList;

typedef struct {
   char name [80];
   int displayList;
   FaceList faceList;
   int n;
} MaterialFaceInfo;

static int nmat = 0;
static int nvtx = 0;
static int nface = 0;

static VertexInfo * vtx = NULL;
static FaceInfo * face = NULL;
static MaterialFaceInfo * matface = NULL;

/*---------------------------------------------------------------------------
 *
 *  The following are procedures to search and modify some of the lists
 *  used in keeping track of faces, vertices, materials, and smoothing groups
 *
 *---------------------------------------------------------------------------*/

static void 
AllocMatFace ()
{
   /* Allocate and initialize vector matface */
   int i;
   
   assert (nmat >= 0);
   assert (matface == NULL);
   
   matface = (MaterialFaceInfo*) malloc (sizeof (MaterialFaceInfo) * nmat);
   assert (matface != NULL);

   for (i = 0; i < nmat; ++i) {
      matface [i].name [0] = '\0';
      matface [i].displayList = -1;
      matface [i].n = 0;
      matface [i].faceList = (FaceList) NULL;
   }
}

static void 
FreeMatFace ()
{
   /* Frees the matface table */
   int i, j;
   FaceList flist, tmp;

   if (matface == NULL) return;

   for (i = 0; i < nmat; ++i) {
      flist = matface [i].faceList;
      for (j = 0; j < matface [i].n; ++j) {
	 assert (flist != NULL);
	 tmp = flist->next;
	 free (flist);
	 flist = tmp;
      }
      assert (flist == NULL);
   }
   free (matface);
   matface = NULL;
}

static void 
AllocFace ()
{
   /* Allocates and initializes the face table */
   assert (face == NULL);
   assert (nface > 0);

   face = (FaceInfo*) malloc (sizeof (FaceInfo) * nface);
   assert (face != NULL);
}

static void 
FreeFace ()
{
   /* Frees memory associated with the face table */
   free (face);
   face = NULL;
   nface = 0;
}

static void
AllocVtx ()
{
   /* Allocates and initializes the vtx table */
   int i;
   assert (vtx == NULL);
   assert (nvtx > 0);
   
   vtx = (VertexInfo*) malloc (sizeof (VertexInfo) * nvtx);
   assert (vtx != NULL);
   
   for (i = 0; i < nvtx; i++) {
      vtx [i].smoothList = (SmoothList) NULL;
   }
}

static void
FreeVtx ()
{
   /* Deallocates heap space associated with the vtx array */
   int i;
   for (i = 0; i < nvtx; i++) {
      SmoothList ptr = vtx [i].smoothList, aux;
      while (ptr != NULL) {
	 aux = ptr->next;
	 free (ptr);
	 ptr = aux;
      }
   }

   free (vtx);
   vtx = NULL;
   nvtx = 0;
}
	 
static int 
FindMatFace ()
{
   /* Finds an unused matface entry in 'matface' and returns its index */
   int i; 

   if (matface == NULL) {
      AllocMatFace ();
   }

   for (i = 0; i < nmat; i++) {
      if (matface [i].faceList == NULL) return i;
   }
   return -1;
}

   
static void
FindSmooth (int ivtx, dword group)
{
   /* Finds the in the list of normals/smoothing groups of vertex ivtx
    * the entry corresponding to 'group'. If it does not exist yet, one
    * is created. The entry is moved to the head of the list in order
    * to speed up new searches
    */
   SmoothList *ptr;

   assert (vtx != NULL);
   assert (ivtx < nvtx);
   
   ptr = &(vtx [ivtx].smoothList);
   while (*ptr != NULL) {
      if ((*ptr)->group == group) break;
      ptr = &((*ptr)->next);
   }
   if (*ptr == NULL) {
      *ptr = (SmoothList) malloc (sizeof (struct SmoothInfoNode));
      (*ptr)->group = group;
      (*ptr)->normal [0] = (*ptr)->normal [1] = (*ptr)->normal [2] = 0.0;
      (*ptr)->next = (SmoothList) NULL;
   }
   if (ptr != &(vtx [ivtx].smoothList)) {
      /* Move to front */
      SmoothList aux = *ptr;
      *ptr = aux->next;
      aux->next = vtx [ivtx].smoothList;
      vtx [ivtx].smoothList = aux;
   }
}

static void
SmoothVertices ()
{
   /* Average all vertex normals */
   int i;
   SmoothList ptr;

   assert (vtx != NULL);
   for (i = 0; i < nvtx; i++) {
      for (ptr = vtx [i].smoothList; ptr != NULL; ptr = ptr->next) {
	 float hyp = (float)sqrt (ptr->normal [0]*ptr->normal [0]+
			   ptr->normal [1]*ptr->normal [1]+
			   ptr->normal [2]*ptr->normal [2]);
	 if (hyp == 0.0) continue;
	 ptr->normal [0] /= hyp;
	 ptr->normal [1] /= hyp;
	 ptr->normal [2] /= hyp;
      }
   }
}

static int 
RenderFaceList (FaceList ptr)
{
   /* Renders all triangular faces in list 'ptr'. Returns
    * the number of the display list used.
    */
   int disp = glGenLists (1);
   assert (disp >= 0);
	
   glNewList (disp, GL_COMPILE);
   glBegin (GL_TRIANGLES);

   for (;ptr != NULL; ptr = ptr->next) {
      int iface = ptr->idx;
      dword group = face [iface].group;
      int i;
      for (i = 0; i < 3; i++) {
	 int ivtx = face [iface].ivtx [i];
	 FindSmooth (ivtx, group);	 
	 glNormal3fv (vtx [ivtx].smoothList->normal);
	 if (vtx [ivtx].texFlag) {
	    glTexCoord2fv (vtx [ivtx].texCoord);
	 }
	 glVertex3fv (vtx [ivtx].coord);
      }
   }
   glEnd ();
   glEndList();
   return disp;
}

static void 
RenderAllFaces (Tcl_DString* desc)
{
   int imatface, iface;
   char buf [80];

   SmoothVertices ();

   if (matface == NULL) {
      /* No material/face list was given: render all faces with
       * flat normals and material = {} */
      int disp = glGenLists(1);
      assert (disp >= 0);

      Tcl_DStringStartSublist (desc);
      Tcl_DStringAppendElement (desc, "material-face");

      Tcl_DStringStartSublist (desc);
      Tcl_DStringAppendElement (desc, "material");
      Tcl_DStringAppendElement (desc, "");
      Tcl_DStringEndSublist (desc);

      Tcl_DStringStartSublist (desc);
      Tcl_DStringAppendElement (desc, "nfaces");
      sprintf (buf, "%d", nface);
      Tcl_DStringAppendElement (desc, buf);
      Tcl_DStringEndSublist (desc);

      Tcl_DStringStartSublist (desc);
      Tcl_DStringAppendElement (desc, "displaylist");
      sprintf (buf, "%d", disp);
      Tcl_DStringAppendElement (desc, buf);
      Tcl_DStringEndSublist (desc);

      Tcl_DStringEndSublist (desc);

      glNewList (disp, GL_COMPILE);
      glBegin (GL_TRIANGLES);
      for (iface = 0; iface < nface; iface++) {
	 int i;
	 glNormal3fv (face [iface].normal);
	 for (i = 0; i < 3; i++) {
	    int ivtx = face [iface].ivtx [i];
	    glVertex3fv (vtx [ivtx].coord);
	    if (vtx [ivtx].texFlag) {
	       glTexCoord2fv (vtx [ivtx].texCoord);
	    }
	 }
      }
      glEnd ();
      glEndList();
   }
   else {
      /* Render faces in separate display lists according to matface */
      for (imatface = 0; imatface < nmat; imatface++) {
	 FaceList flist = matface [imatface].faceList;
	 if (flist == NULL) continue;
	 matface [imatface].displayList = RenderFaceList (flist);

	 Tcl_DStringStartSublist (desc);
	 Tcl_DStringAppendElement (desc, "material-face");

	 Tcl_DStringStartSublist (desc);
	 Tcl_DStringAppendElement (desc, "material");
	 Tcl_DStringAppendElement (desc, matface [imatface].name);
	 Tcl_DStringEndSublist (desc);

	 Tcl_DStringStartSublist (desc);
	 Tcl_DStringAppendElement (desc, "nfaces");
	 sprintf (buf, "%d", matface [imatface].n);
	 Tcl_DStringAppendElement (desc, buf);
	 Tcl_DStringEndSublist (desc);

	 Tcl_DStringStartSublist (desc);
	 Tcl_DStringAppendElement (desc, "displaylist");
	 sprintf (buf, "%d", matface [imatface].displayList);
	 Tcl_DStringAppendElement (desc, buf);
	 Tcl_DStringEndSublist (desc);

	 Tcl_DStringEndSublist (desc);
      }
   }
}

/*---------------------------------------------------------------------------
 *
 *  The following are procedures for reading simple types from
 *  a file stream. They require a fair amount of byte swapping in order
 *  to cope with the different way Intel microprocessors store ints, floats,
 *  etc. All procs return 1 on success and 0 if fail.
 *
 *---------------------------------------------------------------------------*/

static int 
ReadWord (FILE *f, word *wptr) 
{
#if defined(__WIN32__) || defined(_WIN32) || defined(LINUX)
   unsigned char b [2];
   if (fread (b, 2, 1, f) != 1) return 0;
#else
   unsigned char b [2], tmp;
   if (fread (b, 2, 1, f) != 1) return 0;
   tmp = b [0];
   b [0] = b [1];
   b [1] = tmp;
#endif

   memcpy ((char*) wptr, b, 2);

   return 1;
}

static int 
ReadWords (FILE *f, word *wptr, int nwords) 
{
   while (nwords--) {
      if (ReadWord (f, wptr++) != 1) return 0;
   }
   return 1;
}

static int 
ReadDWord (FILE *f, dword *wptr) 
{
#if defined(__WIN32__) || defined(_WIN32) || defined(LINUX)
   word b [2];
   if (ReadWords (f, b, 2) != 1) return 0;
#else
   word b [2], tmp;
   if (ReadWords (f, b, 2) != 1) return 0;
   tmp = b [0];
   b [0] = b [1];
   b [1] = tmp; 
#endif

   memcpy ((char*) wptr, (char*) b, 4);
   return 1;
}

static int 
ReadFloat (FILE *f, float* fptr) 
{
#if defined(__WIN32__) || defined(_WIN32) || defined(LINUX)
   word b [2];
   if (ReadWords (f, b, 2) != 1) return 0;
#else
   word b [2], tmp;
   if (ReadWords (f, b, 2) != 1) return 0;
   tmp = b [0];
   b [0] = b [1];
   b [1] = tmp; 
#endif

   memcpy ((char*) fptr, b, 4);
   return 1;
}

static int 
ReadFloats (FILE *f, float *fptr, int nfloats)
{
   while (nfloats--) {
      if (ReadFloat (f, fptr++) != 1) return 0;
   }
   return 1;
}


/*---------------------------------------------------------------------------
 *
 *  These are the various procedures for reading each chnk type. The
 *  default behaviour is to read the chunk and 
 *  store the appropriate information
 *  about the chunk in the Tcl Desc result string 
 *
 *---------------------------------------------------------------------------*/

static int 
FindChunk(word id) 
{
    int i;
    for (i = 0; i < sizeof(ChunkNames)/sizeof(ChunkNames[0]); i++)
        if (id == ChunkNames[i].id)
            return i;
    return -1;
}

static void 
ChunkReader(Tcl_DString* desc, FILE *f, long p) 
{
    TChunkHeader h;
    int n;
    long pc;
/*    char buf [80];*/

    while (ftell(f) < p) {
       pc = ftell(f);
       if (ReadWord (f, &h.id) != 1 ||
	   ReadDWord (f, &h.len) != 1) return;
/*
       printf ("Chunk %x len = %d\n", h.id, h.len); fflush (stdout);
*/
       n = FindChunk(h.id);
       if (n < 0) {
/*
           printf ("Unknown id\n"); fflush (stdout);
	   sprintf (buf, "unknown {id 0x%04X} {offset 0x%lX} {size %d}",
		    h.id, pc, h.len);
	   Tcl_DStringAppendElement (desc, buf);
*/
	  if (h.len == 0) return;	     
	  fseek(f, pc + h.len, SEEK_SET);
       } else {
	  Tcl_DStringStartSublist (desc);
	  Tcl_DStringAppendElement (desc, ChunkNames[n].name);
	  pc = pc + h.len;
	  if (ChunkNames[n].func != NULL)
	     ChunkNames[n].func(desc, f, pc);
	  else
	     ChunkReader(desc, f, pc);
	  Tcl_DStringEndSublist (desc);
	  fseek(f, pc, SEEK_SET);
       }
       if (ferror(f))
	  break;
    }
}

static void 
TriMeshReader(Tcl_DString* desc, FILE *f, long p) 
{
   TChunkHeader h;
   int n;
   long pc;

   while (ftell(f) < p) {
      pc = ftell(f);
      if (ReadWord (f, &h.id) != 1 ||
	  ReadDWord (f, &h.len) != 1) return;
      n = FindChunk(h.id);
      if (n < 0) {
	 fseek(f, pc + h.len, SEEK_SET);
      } else {
	 Tcl_DStringStartSublist (desc);
	 Tcl_DStringAppendElement (desc, ChunkNames[n].name);
	 pc = pc + h.len;
	 if (ChunkNames[n].func != NULL)
	    ChunkNames[n].func(desc, f, pc);
	 else
	    ChunkReader(desc, f, pc);
	 
	 Tcl_DStringEndSublist (desc);
	 fseek(f, pc, SEEK_SET);
      }
      if (ferror(f))
	 break;
   }
    
   RenderAllFaces (desc);
   FreeVtx ();
   FreeMatFace ();
   FreeFace ();
}


static void 
RGBFReader (Tcl_DString* desc, FILE *f, long p) 
{
   float c[3];
   char buf [80];
   int i;
   if (ReadFloats (f, c, 3) != 1) return;
   for (i = 0; i < 3; ++i) {
      sprintf (buf, "%f", c[i]);
      Tcl_DStringAppendElement (desc, buf);
   }
}

static void 
RGBBReader (Tcl_DString* desc, FILE *f, long p) 
{
   byte c[3];
   char buf [80];
   int i;
   if (fread(&c, sizeof(c), 1, f) != 1) return;
   for (i = 0; i < 3; ++i) {
      sprintf (buf, "%f", c[i]/255.0);
      Tcl_DStringAppendElement (desc, buf);
   }
}

static void 
PercentIReader (Tcl_DString* desc, FILE *f, long p) 
{
   word perc;
   char buf [80];
   if (ReadWord (f, &perc) != 1) return;
   sprintf (buf, "%d", perc);
   Tcl_DStringAppendElement (desc, buf);
}

static void 
PercentFReader (Tcl_DString* desc, FILE *f, long p) 
{
   float perc;
   char buf [80];
   if (ReadFloat (f, &perc) != 1) return;
   sprintf (buf, "%f", perc);
   Tcl_DStringAppendElement (desc, buf);
}

static void 
ObjBlockReader (Tcl_DString* desc, FILE *f, long p) 
{
   int i, c;
   char buf [80];
    
   /* Read ASCIIZ object name */
   for (i = 0; i < 80; i++) {
       c = fgetc(f);
       if (c == EOF || c == '\0') break;
       buf [i] = c;
   }
   buf [i] = '\0';
   
   Tcl_DStringStartSublist (desc);
   Tcl_DStringAppendElement (desc, "name");
   Tcl_DStringAppendElement (desc, buf);
   Tcl_DStringEndSublist (desc);
    
   /* Read rest of chunks inside this one. */
   ChunkReader(desc, f, p);
}

static void
VertListReader (Tcl_DString* desc, FILE *f, long p)
{
   word ivtx, nv, i;
   float c[3];
   char buf [80];
   Vector max;
   Vector min;

   if (ReadWord (f, &nv) != 1) return;

   sprintf (buf, "%d", nv);

   Tcl_DStringStartSublist (desc);
   Tcl_DStringAppendElement (desc, "vertices");
   Tcl_DStringAppendElement (desc, buf);
   Tcl_DStringEndSublist (desc);

   nvtx = nv;
   AllocVtx ();

   for (ivtx = 0; ivtx < nv; ivtx++) {
      if (ReadFloats (f, c, 3) != 1) assert (0);
      vtx [ivtx].texFlag = 0;
      if (ivtx == 0) {
         for (i = 0; i < 3; i++) {
            vtx [ivtx].coord [i] = c [i];
            min [i] = c [i];
            max [i] = c [i];
         }
      }
      else {
         for (i = 0; i < 3; i++) {
	    vtx [ivtx].coord [i] = c [i];
	    if (c [i] > max [i]) max [i] = c [i];
	    if (c [i] < min [i]) min [i] = c [i];
         }
      }
   }

   Tcl_DStringStartSublist (desc);
   Tcl_DStringAppendElement (desc, "min");
   for (i = 0; i < 3; i++) {
      sprintf (buf, "%f", min [i]);
      Tcl_DStringAppendElement (desc, buf);
   }
   Tcl_DStringEndSublist (desc);

   Tcl_DStringStartSublist (desc);
   Tcl_DStringAppendElement (desc, "max");
   for (i = 0; i < 3; i++) {
      sprintf (buf, "%f", max [i]);
      Tcl_DStringAppendElement (desc, buf);
   }
   Tcl_DStringEndSublist (desc);

}

static void
FaceListReader (Tcl_DString* desc, FILE *f, long p) 
{
   word iface, nv;
   word c[4];
   char buf [80];

   if (ReadWord (f, &nv) != 1) return;

   sprintf (buf, "%d", nv);
   Tcl_DStringStartSublist (desc);
   Tcl_DStringAppendElement (desc, "faces");
   Tcl_DStringAppendElement (desc, buf);
   Tcl_DStringEndSublist (desc);

   nface = nv;
   AllocFace ();
   for (iface = 0; iface < nface; iface++) {
      Vector e1, e2;
      float hypot;
      if (ReadWords (f, c, 4) != 1) assert(0);
      face [iface].ivtx [0] = c[0];
      face [iface].ivtx [1] = c[1];
      face [iface].ivtx [2] = c[2];
      e1 [0] = vtx[c[1]].coord[0]-vtx[c[0]].coord[0];
      e1 [1] = vtx[c[1]].coord[1]-vtx[c[0]].coord[1];
      e1 [2] = vtx[c[1]].coord[2]-vtx[c[0]].coord[2];
      e2 [0] = vtx[c[2]].coord[0]-vtx[c[1]].coord[0];
      e2 [1] = vtx[c[2]].coord[1]-vtx[c[1]].coord[1];
      e2 [2] = vtx[c[2]].coord[2]-vtx[c[1]].coord[2];
      face [iface].normal [0] = e1[1]*e2[2]-e1[2]*e2[1];
      face [iface].normal [1] = e2[0]*e1[2]-e2[2]*e1[0];
      face [iface].normal [2] = e1[0]*e2[1]-e1[1]*e2[0];
      hypot = (float)sqrt (face [iface].normal [0] * face [iface].normal [0] +
		    face [iface].normal [1] * face [iface].normal [1] +
		    face [iface].normal [2] * face [iface].normal [2]);
      if (hypot != 0.0) {
	 face [iface].normal [0] /= hypot;
	 face [iface].normal [1] /= hypot;
	 face [iface].normal [2] /= hypot;
      }
      face [iface].flags = c[3];
      face [iface].group = 0;
   }

   /* Read rest of chunks inside this one. */
   ChunkReader(desc, f, p);
}


static void 
FaceMatReader (Tcl_DString* desc, FILE *f, long p) 
{
   int c, i, imat;
   word n, nf;
   FaceList *ptr;
   char buf[80];

   /* Read ASCIIZ object name */
   for (i = 0; i < 80; i++) {
      c = fgetc(f);
      if (c == EOF || c == '\0') break;
      buf [i] = c;
   }
   assert (i < 80);
   buf [i] = '\0';

   imat = FindMatFace (buf);
   assert (imat >= 0 && imat < nmat);

   strcpy (matface [imat].name, buf);   
   if (ReadWord (f, &n) != 1) assert (0);
   matface [imat].n = n;

   ptr = &matface [imat].faceList;
   while (n-- > 0) {
      if (ReadWord (f, &nf) != 1) assert (0);
      *ptr = (FaceList) malloc (sizeof (struct FaceListNode));
      assert (*ptr != NULL);
      (*ptr)->next = NULL;
      (*ptr)->idx = nf;
      ptr = &(*ptr)->next;
   }
}

static void 
MapListReader (Tcl_DString* desc, FILE *f, long p) 
{
   word ivtx, nv;
   float c[2];
   char buf[80];

   if (ReadWord (f, &nv) != 1) return;


   Tcl_DStringStartSublist (desc);
   Tcl_DStringAppendElement (desc, "vertexmap");
   sprintf (buf, "%d", nv);
   Tcl_DStringAppendElement (desc, buf);
   Tcl_DStringEndSublist (desc);

   for (ivtx = 0; ivtx < nv; ivtx++) {
      if (ReadFloats (f, c, 2) != 1) return;
      vtx [ivtx].texFlag = 1;
      vtx [ivtx].texCoord [0] = c [0];
      vtx [ivtx].texCoord [1] = c [1];
   }
}

static void 
SmooListReader (Tcl_DString* desc, FILE *f, long p) 
{
   dword s;
   int i, j, n, ivtx;

   n = 0;
   while (ftell(f) < p) {
      if (ReadDWord (f, &s) != 1) assert(0);
      assert (n <= nface);
      face [n].group = s;
      for (i = 0; i < 3; i++) {
	 ivtx = face [n].ivtx[i];
	 FindSmooth (ivtx, s);
	 for (j = 0; j < 3; j++) {
	    vtx [ivtx].smoothList->normal[j] += face [n].normal[j];
	 }
      }
      n++;
   }
   assert (n == nface);
}

static void 
TrMatrixReader(Tcl_DString* desc, FILE *f, long p) 
{
   float rot[12];
   int i;
   char buf [80];

   if (ReadFloats (f, rot, 12) != 1) return;

   for (i = 0; i < 12; ++i) {
      sprintf (buf, "%f", rot [i]);
      Tcl_DStringAppendElement (desc, buf);
   }
}

static void 
LightReader(Tcl_DString* desc, FILE *f, long p) 
{
   float c[3];
   int i;
   char buf [80];

   if (ReadFloats (f, c, 3) != 1) return;

   Tcl_DStringStartSublist (desc);
   Tcl_DStringAppendElement (desc, "position");
   for (i = 0; i < 3; ++i) {
      sprintf (buf, "%f", c[i]);
      Tcl_DStringAppendElement (desc, buf);
   }
   Tcl_DStringEndSublist (desc);

   /* Read rest of chunks inside this one. */
   ChunkReader(desc, f, p);
}

static void 
SpotLightReader(Tcl_DString* desc, FILE *f, long p)
{
   float c[5];
   char buf [80];
   int i;

   if (ReadFloats (f, c, 5) != 1) return;

   Tcl_DStringStartSublist (desc);
   Tcl_DStringAppendElement (desc, "target");
   for (i = 0; i < 3; i++) {
      sprintf (buf, "%f", c[i]);
      Tcl_DStringAppendElement (desc, buf);
   }
   Tcl_DStringAppend (desc, "} {hotspot", -1);
   sprintf (buf, "%f", c[3]);
   Tcl_DStringAppendElement (desc, buf);
   Tcl_DStringAppend (desc, "} {falloff", -1);
   sprintf (buf, "%f", c[4]);
   Tcl_DStringAppendElement (desc, buf);
   Tcl_DStringEndSublist (desc);

}
 
static void 
CameraReader(Tcl_DString* desc, FILE *f, long p)
{
   float c[8];
   char buf [80];
   int i;

   if (ReadFloats (f, c, 8) != 1) return;

   Tcl_DStringStartSublist (desc);
   Tcl_DStringAppendElement (desc, "position");
   for (i = 0; i < 3; i++) {
      sprintf (buf, "%f", c[i]);
      Tcl_DStringAppendElement (desc, buf);
   }
   Tcl_DStringAppend (desc, "} {target", -1);
   for (i = 3; i < 6; i++) {
      sprintf (buf, "%f", c[i]);
      Tcl_DStringAppendElement (desc, buf);
   }

   Tcl_DStringAppend (desc, "} {bank", -1);
   sprintf (buf, "%f", c[6]);
   Tcl_DStringAppendElement (desc, buf);

   Tcl_DStringAppend (desc, "} {lens", -1);
   sprintf (buf, "%f", c[7]);
   Tcl_DStringAppendElement (desc, buf);
   Tcl_DStringEndSublist (desc);

}

static void 
MatNameReader (Tcl_DString* desc, FILE *f, long p)
{
   int i, c;
   char buf [80];

   /* Read ASCIIZ object name */
   for (i = 0; i < 80; i++) {
      c = fgetc(f);
      if (c == EOF || c == '\0') break;
      buf [i] = c;
   }
   buf [i] = '\0';

   nmat++;

   Tcl_DStringAppendElement (desc, buf);
}

static void 
MapFileReader (Tcl_DString* desc, FILE *f, long p)
{
   int i, c;
   char buf [80];

   /* Read ASCIIZ object name */
   for (i = 0; i < 80; i++) {
      c = fgetc(f);
      if (c == EOF || c == '\0') break;
      buf [i] = c;
   }
   buf [i] = '\0';

   Tcl_DStringAppendElement (desc, buf);
}

static void 
FramesReader (Tcl_DString* desc, FILE *f, long p)
{
   dword c[2];
   char buf[80];

   if (ReadDWord (f, &c[0]) != 1 || 
       ReadDWord (f, &c[1]) != 1) return;

   Tcl_DStringStartSublist (desc);
   Tcl_DStringAppend (desc, "start", -1);
   sprintf (buf, "%d", (int)c[0]);
   Tcl_DStringAppendElement (desc, buf);

   Tcl_DStringAppend (desc, "} {end", -1);
   sprintf (buf, "%d", (int)c[1]);
   Tcl_DStringAppendElement (desc, buf);
   Tcl_DStringEndSublist (desc);
}



/*
 *--------------------------------------------------------------------------- 
 *
 *  This procedure loads a 3ds file created by the 3D-Studio software
 *  package and pipes any triangles it may find into the GL pipeline
 *
 *  The syntax is:
 *
 *  	<pathName> "load3ds" <filename> 
 *
 *  where 
 *	<filename> is the name of an ASC file
 *
 *  Result: If OK, a list where each element is a property list of 
 * 	objects described in the file. Otherwise, an error message.
 *
 *---------------------------------------------------------------------------
 */

int
glLoad3DStudio (interp, argc, argv)
     Tcl_Interp *interp;		/* Current interpreter. */
     int argc;			/* Arg count */
     char *argv [];		/* Argument strings. */
{
   int result = TCL_OK;
   long p;
   FILE * file3D;
   Tcl_DString desclist;
   char* filename;

   if (argc != 3) {
      Tcl_AppendResult (interp, "wrong # args", (char*) NULL);
      return TCL_ERROR;
   }

   filename = argv [2];

#if defined(__WIN32__) || defined(_WIN32)
   file3D = fopen (filename, "rb");
#else
   file3D = fopen (filename, "r");
#endif
    
   if (file3D == NULL) { 
      Tcl_AppendResult (interp, "Could not read ", filename, (char*) NULL);
      return TCL_ERROR;
   }

   /* Find file size. */
   fseek(file3D, 0, SEEK_END);
   p = ftell(file3D);
   fseek(file3D, 0, SEEK_SET);

   /* Allocate a dynamic string to hold the result */
   Tcl_DStringInit (&desclist);

   /* Go! */
   ChunkReader(&desclist, file3D, p);

   fclose (file3D);

   Tcl_DStringResult (interp, &desclist);
   Tcl_DStringFree (&desclist);

   return result;
}


