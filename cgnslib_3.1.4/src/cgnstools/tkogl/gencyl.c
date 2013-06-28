#include <tk.h>
#if defined(__WIN32__) || defined(_WIN32)
#   define WIN32_LEAN_AND_MEAN
#   include <windows.h>
#   undef WIN32_LEAN_AND_MEAN
#endif
#include <GL/gl.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "gencyl.h"
#include "tkogl.h"

#ifndef PI
/* The ubiquitous PI constant */
#define PI 3.14159265358979323846264338327950288
#endif

/* Constants that define the various flags for rendering a generic
 * cylinder */

#define STITCH_ENDS 0x01 /* Stitch together the ends of the cyl */
#define STITCH_LOOPS 0x02 /* Close the cross-sections of the cyl */
#define SHADE_SMOOTH_ROWS 0x04 /* Average normals across rows */
#define SHADE_SMOOTH_COLS 0x08 /* Average normals across columns */
#define TEXTURE_GEN 0x10 /* Generate texture coordinates */
#define ADAPTIVE 0x20 /* Average only normals which form a small angle */
#define CLOSE_FIRST 0x40 /* Close first cross section */
#define CLOSE_LAST 0x80 /* Close last cross section */

/*---------------------------------------------------------------------------
 *
 *  General operations with vectors and transformation matrices
 *
 *---------------------------------------------------------------------------*/

typedef GLfloat Vector [3];
typedef GLfloat Matrix [4][4];

static void 
MatrixMult (Matrix m1, Matrix m2, Matrix result) 
{
   /* Multiplies m1 by m2 storing the result in result */
   int i, j, k;
   Matrix tmp;

   for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) { 
         tmp [i][j] = 0.0;
         for (k = 0; k < 4; k++) {
            tmp [i][j] += m1 [i][k] * m2 [k][j];
         }
      }
   }
   for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
         result [i][j] = tmp [i][j];
      }
   }
}

static void 
TransformVector (const Vector v, Matrix m, Vector result)
{
   /* Applies affine linear transformation m to v storing the result in 
    * result */
   int i, j;
   Vector tmp;

   for (i = 0; i < 3; i++) {
      tmp [i] = m [i][3];
      for (j = 0; j < 3; j++) {
         tmp [i] += m [i][j] * v [j];
      }
   }
   for (i = 0; i < 3; i++) {
      result [i] = tmp [i];
   }
}

static void 
NormalizeVector (Vector v)
{
   /* Makes v a unit vector */
   GLfloat hypot = (GLfloat) sqrt (v [0] * v [0] +
                                   v [1] * v [1] +
                                   v [2] * v [2]);
   if (hypot == 0.0) {
      return;
   }
   v [0] /= hypot;
   v [1] /= hypot;
   v [2] /= hypot;
}

static void 
AddVector (const Vector v1, const Vector v2, Vector result)
{
   /* Adds v1 to v2 and stores in result */
   result [0] = v1 [0] + v2 [0];
   result [1] = v1 [1] + v2 [1];
   result [2] = v1 [2] + v2 [2];
}

#define InitVector(v,a,b,c) {v[0]=a; v[1]=b; v[2]=c;}
#define CopyVector(dst,src) {dst[0]=src[0];dst[1]=src[1];dst[2]=src[2];}
#define DotProduct(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

#if 0
static void
ComputeTriangleNormal2 (const Vector v1, const Vector v2, const Vector v3,
                       Vector normal)
{
   /* Computes the normal of the triangle given by the three vertices v1, v2, v3
    * and stores it in the vector given by normal */
   Vector e1, e2; 
   InitVector (e1, v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]);
   InitVector (e2, v3[0] - v2[0], v3[1] - v2[1], v3[2] - v2[2]);
   normal [0] = e1[1]*e2[2]-e1[2]*e2[1];
   normal [1] = e2[0]*e1[2]-e2[2]*e1[0];
   normal [2] = e1[0]*e2[1]-e1[1]*e2[0];
}

static void
ComputeQuadNormal (const Vector v1, const Vector v2, const Vector v3, 
                   const Vector v4, Vector normal)
{
   /* Computes the normal at vertex v1 of the quadrilateral given
    * by v1-v2-v3-v4 and stores it in the vector given by normal */
   Vector normal1, normal2;
   GLfloat size1, size2;
   ComputeTriangleNormal2 (v1, v2, v3, normal1);
   ComputeTriangleNormal2 (v3, v4, v1, normal2);
   size1 = normal1[0]*normal1[0]+normal1[1]*normal1[1]+normal1[2]*normal1[2];
   size2 = normal2[0]*normal2[0]+normal2[1]*normal2[1]+normal2[2]*normal2[2];
   if (size1 > 100*size2) {
      normal [0] = normal1 [0]/size1;
      normal [1] = normal1 [1]/size1;
      normal [2] = normal1 [1]/size1;
   } else if (size2 > 100*size1) {
      normal [0] = normal2 [0]/size2;
      normal [1] = normal2 [1]/size2;
      normal [2] = normal2 [1]/size2;
   } else {
      AddVector (normal1, normal2, normal);
      NormalizeVector (normal);
   }
}
#endif

static void                                                                               
ComputeQuadNormal2 (const Vector v1, const Vector v2, const Vector v3,                    
                    const Vector v4, Vector normal)  
{
   /* Computes the squared normal at vertex v1 of the quadrilateral given
    * by v1-v2-v3-v4 and stores it in the vector given by normal */

   normal [0] = (v4 [1] - v1 [1]) * (v4 [2] + v1[2]);
   normal [1] = (v4 [2] - v1 [2]) * (v4 [0] + v1[0]);
   normal [2] = (v4 [0] - v1 [0]) * (v4 [1] + v1[1]);

   normal [0] += (v1 [1] - v2 [1]) * (v1 [2] + v2[2]);
   normal [1] += (v1 [2] - v2 [2]) * (v1 [0] + v2[0]);
   normal [2] += (v1 [0] - v2 [0]) * (v1 [1] + v2[1]);

   normal [0] += (v2 [1] - v3 [1]) * (v2 [2] + v3[2]);
   normal [1] += (v2 [2] - v3 [2]) * (v2 [0] + v3[0]);
   normal [2] += (v2 [0] - v3 [0]) * (v2 [1] + v3[1]);

   normal [0] += (v3 [1] - v4 [1]) * (v3 [2] + v4[2]);
   normal [1] += (v3 [2] - v4 [2]) * (v3 [0] + v4[0]);
   normal [2] += (v3 [0] - v4 [0]) * (v3 [1] + v4[1]);

}


/*---------------------------------------------------------------------------
 *
 *  Management of cross sections
 *
 *---------------------------------------------------------------------------*/

typedef struct {
   Vector * vtx; /* Vertex Coordinates */
   Vector * normalUL,  /* Normal to the Up Left Side of the vertex */
          * normalUR,  /* Normal to the Up Right Side of the vertex */
          * normalDL,  /* Normal to the Down Left Side of the vertex */  
          * normalDR ; /* Normal to the Down Right Side of the vertex */ 
   int nvtx,     /* Number of vertices in cross section */
       vtxsize;  /* Number of vertices allocated in vtx */
} CrossSection;

static CrossSection* 
NewCrossSection ()
{
   /* Allocates a new cross section structure */
   CrossSection * result = (CrossSection*) malloc (sizeof (CrossSection));
   assert (result != NULL);
   result->nvtx = 0;
   result->vtxsize = 8;
   result->vtx = (Vector*) malloc (sizeof (Vector)*result->vtxsize);
   result->normalUL = result->normalUR = result->normalDL = result->normalDR = NULL;
   assert (result->vtx != NULL);
   return result;
}

static void 
FreeCrossSection (CrossSection * s)
{
   /* Deallocates the memory associated with cross-section s */
   assert (s != NULL);
   assert (s->vtx != NULL);
   if (s->normalUL != NULL) free(s->normalUL);
   if (s->normalUR != NULL) free(s->normalUR);
   if (s->normalDL != NULL) free(s->normalDL);
   if (s->normalDR != NULL) free(s->normalDR);
   free (s->vtx);
   free (s);
}

static void
AddCrossSectionVtx (CrossSection * s, GLfloat x, GLfloat y, GLfloat z) 
{
   /* Stores another vertex with coords (x,y,0) as the next vertex in
    * cross section 's' */
   assert (s != NULL);
   if (s->vtxsize == s->nvtx) {
      s->vtxsize *= 2;
      s->vtx = (Vector*) realloc (s->vtx, sizeof(Vector)*s->vtxsize);
      assert (s->vtx != NULL);
   }
   s->vtx [s->nvtx][0] = x;
   s->vtx [s->nvtx][1] = y;
   s->vtx [s->nvtx][2] = z;
   s->nvtx++;
}

static CrossSection* 
PolygonCrossSection (GLfloat radius, int nsides)
{
   /* Returns a cross section which is a regular polygon with 'nsides' sides
    * and radius 'radius' */
   GLfloat x, y, ang, incAng;
   int i;
   CrossSection *cross = NewCrossSection ();
   incAng = (GLfloat)(PI * 2 / nsides);
   ang = 0.0;
   for (i = 0; i < nsides; i++) {
      x = (GLfloat)(radius * sin (ang));
      y = (GLfloat)(radius * cos (ang));
      AddCrossSectionVtx (cross, x, y, 0.0);
      ang += incAng;
   }
   return cross;
}

static void 
TransformCrossSection (CrossSection* s, Matrix m)
{
   /* Apply transformation m to cross section s */
   int i;
   for (i = 0; i < s->nvtx; i++) {
      TransformVector (s->vtx [i], m, s->vtx [i]);
   }
}

static void 
DupCrossSection (const CrossSection* src, CrossSection* dst)
{
   /* Make dst an exact copy of src */
   if (dst->vtxsize < src->vtxsize) {
      dst->vtxsize = src->vtxsize;
      dst->vtx = (Vector*) realloc (dst->vtx, sizeof(Vector)*dst->vtxsize);
      assert (dst->vtx != NULL);
   }
   dst->nvtx = src->nvtx;
   memcpy (dst->vtx, src->vtx, sizeof (Vector)*dst->nvtx);
}

static void 
ExpandCrossSection (const CrossSection* src, CrossSection* dst, int ndst)
{
   /* Make dst a copy of src expanded to have n vertices */
   int i, nsrc;
   nsrc = src->nvtx;
   assert (nsrc < ndst);
   if (dst->vtxsize < ndst) {
      dst->vtxsize = ndst;
      dst->vtx = (Vector*) realloc (dst->vtx, sizeof (Vector)*ndst);
   }
   dst->nvtx = ndst;
   for (i = 0; i < ndst; i++) {
      memcpy (&(dst->vtx [i]), &(src->vtx [i*nsrc/ndst]), sizeof(Vector));
   }
}

static void
CrossSectionNormal2 (const CrossSection* cross, Vector normal)
{
   /* Computes a non-normalized vector that is normal to the cross section */
   int ivtx, prev;
   int nvtx = cross->nvtx;
   Vector * vtx = cross->vtx;
   prev = nvtx;
   normal [0] = normal [1] = normal [2] = 0.0;
   for (ivtx = 0; ivtx < cross->nvtx; ivtx++) {
      normal [0] += (vtx[prev] [1] - vtx[ivtx] [1]) * (vtx[prev] [2] + vtx[ivtx][2]);
      normal [1] += (vtx[prev] [2] - vtx[ivtx] [2]) * (vtx[prev] [0] + vtx[ivtx][0]);
      normal [2] += (vtx[prev] [0] - vtx[ivtx] [0]) * (vtx[prev] [1] + vtx[ivtx][1]);
      prev = ivtx;
   }
}

/*---------------------------------------------------------------------------
 *
 * Model management
 *
 *---------------------------------------------------------------------------*/

typedef struct {
   double smin, smax, tmin, tmax;
   int flags;
   int ncross, sizecross;
   double adaptivethreshold;
   CrossSection** cross;
} Model;

static Model * 
NewModel () 
{
   /* Allocates a new model and returns a pointer to it */
   Model* result;
   
   result = (Model*) malloc (sizeof (Model));
   result->adaptivethreshold = 0;
   result->ncross = 0;
   result->sizecross = 8;
   result->cross = (CrossSection**) malloc (sizeof (CrossSection*) * 
                                            result->sizecross);
   assert (result->cross != NULL);
   return result;
}

static void 
FreeModel (Model* model)
{      
   /* Deallocates all memory associated with model model */
   int i;
   for (i = 0; i < model->ncross; i++) {
      FreeCrossSection (model->cross [i]);
   }
   free (model->cross);
}

static void
AddModelCrossSection (Model* model, CrossSection* cross)
{
   /* Adds another cross section to the model */
   if (model->sizecross == model->ncross) {
      model->sizecross *= 2;
      model->cross = (CrossSection**) realloc (model->cross,
                                sizeof (CrossSection*) * model->sizecross);
   }
   model->cross [model->ncross++] = cross;
}

static void
UniformCrossSectionLengths (Model * model) 
{
   /* Force all CrossSections to be of uniform length */
   int icross, maxlength;
   maxlength = 0;
   for (icross = 0; icross < model->ncross; icross++) {   
      if (model->cross [icross]->nvtx > maxlength) {
         maxlength = model->cross [icross]->nvtx;
      }
   }
   for (icross = 0; icross < model->ncross; icross++) {   
      if (model->cross [icross]->nvtx < maxlength) {
         CrossSection* tmp = NewCrossSection ();
         ExpandCrossSection (model->cross [icross], tmp, maxlength);
         FreeCrossSection (model->cross [icross]);
         model->cross [icross] = tmp;
      }
   }
}

static void
ComputeModelNormals (Model* model)
{
   /* Computes normals for each vertex of the model */
   int icross, ivtx, prevvtx;
   int flags = model->flags;
   CrossSection *thisCross, *nextCross, *prevCross = NULL;
   Vector *a, *b, *c, *d;
   int nvtx = model->cross [0]->nvtx; /* Assume every cross section has the 
                                         same number of vertices */
   assert (model->ncross > 1);
   assert (nvtx > 1);

   /* First compute Up Right normals (face normals) */
   for (icross = 0; icross < model->ncross; icross++) {
      thisCross = model->cross [icross];
      assert (thisCross->nvtx == nvtx);
      thisCross->normalUR = (Vector*) malloc (sizeof (Vector) * nvtx);
      if (icross+1 == model->ncross) {
         if (flags&STITCH_ENDS) {
            /* Assume last cross section wraps with first cross section */
            nextCross = model->cross [0];
         } 
         else {
            /* Last Cross section repeats normals at right from previous cross 
               sections */
            assert (prevCross != NULL);
            memcpy (thisCross->normalUR, prevCross->normalUR,
                    sizeof (Vector) * nvtx);
            break;
         }
      } 
      else {
         nextCross = model->cross [icross+1];
      }
      for (ivtx = 0; ivtx < nvtx; ivtx++) {
         if (ivtx+1 == nvtx) {
            if (flags&STITCH_LOOPS) {
               /* Assume last vertex wraps with first */
               b = &(thisCross->vtx [0]);
               c = &(nextCross->vtx [0]);
            } 
            else {
               /* Last Vertex repeats normal above from previous vertex */
               CopyVector (thisCross->normalUR[ivtx],
                           thisCross->normalUR[ivtx-1]);
               break;
            }
         }
         else {
            b = &(thisCross->vtx [ivtx+1]);
            c = &(nextCross->vtx [ivtx+1]);
         }
         a = &(thisCross->vtx [ivtx]);
         d = &(nextCross->vtx [ivtx]);
         ComputeQuadNormal2 (*a, *d, *c, *b, thisCross->normalUR [ivtx]);
         NormalizeVector (thisCross->normalUR [ivtx]);
      }
      prevCross = thisCross;
   }

   /* If only face normals are needed, return here */
   if ((flags&(SHADE_SMOOTH_ROWS | SHADE_SMOOTH_COLS)) == 0) return;
   
   /* Copy normals to the remaining 3 directions */
   for (icross = 0; icross < model->ncross; icross++) {
      thisCross = model->cross [icross];
      thisCross->normalUL = (Vector*) malloc (sizeof (Vector) * nvtx);
      thisCross->normalDR = (Vector*) malloc (sizeof (Vector) * nvtx);
      thisCross->normalDL = (Vector*) malloc (sizeof (Vector) * nvtx);
      if (icross == 0) {
         if (flags&STITCH_ENDS) {
            /* Assume first cross section wraps with last cross section */
            prevCross = model->cross [model->ncross-1];
         } 
         else {
            /* First Cross section repeats normals at left from the right*/ 
            prevCross = thisCross;
         }
      } 
      else {
         prevCross = model->cross [icross-1];
      }
      for (ivtx = 0; ivtx < nvtx; ivtx++) {
         if (ivtx == 0) {
            if (flags&STITCH_LOOPS) {
               /* Assume last vertex wraps with first */
               prevvtx = nvtx-1;
            } 
            else {
               /* First Vertex repeats normal below from above */
               prevvtx = 0;
            }
         }
         else {
            prevvtx = ivtx-1;
         }
         CopyVector (thisCross->normalUL [ivtx], prevCross->normalUR [ivtx]);
         CopyVector (thisCross->normalDR [ivtx], thisCross->normalUR [prevvtx]);
         CopyVector (thisCross->normalDL [ivtx], prevCross->normalUR [prevvtx]);
      }
   }

   /* Smooth Normals */
   for (icross = 0; icross < model->ncross; icross++) {
      thisCross = model->cross [icross];
      for (ivtx = 0; ivtx < nvtx; ivtx++) {
         if ((flags & SHADE_SMOOTH_ROWS)) {
            Vector tmp;
            if (!(flags&ADAPTIVE)||
                DotProduct(thisCross->normalUL [ivtx], 
                           thisCross->normalUR [ivtx])
                > model->adaptivethreshold) {
               AddVector (thisCross->normalUL [ivtx], 
                          thisCross->normalUR [ivtx], tmp);
               NormalizeVector (tmp);
               CopyVector (thisCross->normalUL [ivtx], tmp);
               CopyVector (thisCross->normalUR [ivtx], tmp);
            }
            if (!(flags&ADAPTIVE)||
                DotProduct(thisCross->normalDL [ivtx], 
                           thisCross->normalDR [ivtx])
                > model->adaptivethreshold) {
               AddVector (thisCross->normalDL [ivtx], 
                          thisCross->normalDR [ivtx], tmp);
               NormalizeVector (tmp);
               CopyVector (thisCross->normalDL [ivtx], tmp);
               CopyVector (thisCross->normalDR [ivtx], tmp);
            }
         }
         if (flags & SHADE_SMOOTH_COLS) {
            Vector tmp;
            if (!(flags&ADAPTIVE)||
                DotProduct(thisCross->normalUL [ivtx], 
                           thisCross->normalDL [ivtx])
                > model->adaptivethreshold) {
               AddVector (thisCross->normalUL [ivtx], 
                          thisCross->normalDL [ivtx], tmp);
               NormalizeVector (tmp);
               CopyVector (thisCross->normalUL [ivtx], tmp);
               CopyVector (thisCross->normalDL [ivtx], tmp);
            }
            if (!(flags&ADAPTIVE)||
                DotProduct(thisCross->normalUR [ivtx], 
                           thisCross->normalDR [ivtx])
                > model->adaptivethreshold) {
               AddVector (thisCross->normalUR [ivtx], 
                          thisCross->normalDR [ivtx], tmp);
               NormalizeVector (tmp);
               CopyVector (thisCross->normalUR [ivtx], tmp);
               CopyVector (thisCross->normalDR [ivtx], tmp);
            }
         }
      }
   }
}

static void
ComputeModelBBox (Tcl_Interp* interp, Model* model)
{
   /* Computes the model's bounding box and stores it as two elements
      in the Interp's result: one for the minimum corner and another for
      the maximum corner */
   Vector min;
   Vector max;
   char buf [200];
   int icross, ivtx, icoord, ifirst = 1;
   for (icross = 0; icross < model->ncross; icross++) {
      CrossSection *cross = model->cross [icross];
      for (ivtx = 0; ivtx < cross->nvtx; ivtx++) {
         if (ifirst) {
            for (icoord = 0; icoord < 3; icoord++) {
                  min [icoord] = max [icoord] = cross->vtx [ivtx][icoord];
            }
            ifirst = 0;
         }
         else {
            for (icoord = 0; icoord < 3; icoord++) {
               if (cross->vtx [ivtx][icoord] < min[icoord]) {
                  min [icoord] = cross->vtx [ivtx][icoord];
               }
               if (cross->vtx [ivtx][icoord] > max[icoord]) {
                  max [icoord] = cross->vtx [ivtx][icoord];
               }
            }
         }
      }
   }
   sprintf (buf, "min %f %f %f", min[0], min[1], min[2]);
   Tcl_AppendElement (interp, buf);
   sprintf (buf, "max %f %f %f", max[0], max[1], max[2]);
   Tcl_AppendElement (interp, buf);
}


static void
RenderModel (Model* model)
{
   /* Renders the surface between each pair of cross sections */

   int icross, previcross, ivtx, previvtx;
   int ncross = model->ncross;
   int flags = model->flags;
   int nvtx = model->cross [0]->nvtx;
   int gentex = (flags & TEXTURE_GEN) != 0;
   GLfloat texcoord [2];
   GLfloat texincr [2];
   texcoord [0] = (GLfloat)model->tmin;
   if (flags&STITCH_ENDS) {
      icross = 0;
      previcross = ncross-1;
      texincr [0] = (GLfloat)((model->tmax-model->tmin) / ncross);
   }
   else {
      icross = 1;
      previcross = 0;
      texincr [0] = (GLfloat)((model->tmax-model->tmin) / (ncross-1));
   }
   if (flags&STITCH_LOOPS) {
      texincr [1] = (GLfloat)((model->smax-model->smin) / nvtx);
   } 
   else {
      texincr [1] = (GLfloat)((model->smax-model->smin) / (nvtx-1));
   }
   while (icross < ncross) {
      CrossSection *a = model->cross [previcross];
      CrossSection *b = model->cross [icross];     
      assert (a->nvtx == b->nvtx);
      if (flags & STITCH_LOOPS) {
         ivtx = 0;
         previvtx = nvtx-1;
      } 
      else {
         ivtx = 1;
         previvtx = 0;
      }
      texcoord [1] = (GLfloat)model->smin;

      if (flags & (SHADE_SMOOTH_COLS|SHADE_SMOOTH_ROWS)) {
         /* One normal per vertex */
         for (; ivtx < nvtx; ivtx++) {
            glBegin (GL_TRIANGLE_STRIP);

            if (gentex) glTexCoord2f (texcoord [0], texcoord [1]);
            glNormal3fv (a->normalUR [previvtx]);
            glVertex3fv (a->vtx [previvtx]);

            if (gentex) glTexCoord2f (texcoord [0]+texincr[0], texcoord [1]);
            glNormal3fv (b->normalUL [previvtx]);
            glVertex3fv (b->vtx [previvtx]);

            if (gentex) glTexCoord2f (texcoord [0], texcoord [1]+texincr[1]);      
            glNormal3fv (a->normalDR [ivtx]);
            glVertex3fv (a->vtx [ivtx]);

            if (gentex) glTexCoord2f (texcoord [0]+texincr[0], texcoord [1]+texincr[1]);      
            glNormal3fv (b->normalDL [ivtx]);
            glVertex3fv (b->vtx [ivtx]);

            previvtx = ivtx;
            if (gentex) texcoord [1] += texincr [1];
            glEnd ();
         }
      }
      else {
         for (; ivtx < nvtx; ivtx++) {
            glBegin (GL_TRIANGLE_STRIP);

            if (gentex) glTexCoord2f (texcoord [0], texcoord [1]);
            glNormal3fv (a->normalUR [previvtx]);
            glVertex3fv (a->vtx [previvtx]);

            if (gentex) glTexCoord2f (texcoord [0]+texincr[0], texcoord [1]);
            glVertex3fv (b->vtx [previvtx]);

            if (gentex) glTexCoord2f (texcoord [0], texcoord [1]+texincr[1]);      
            glVertex3fv (a->vtx [ivtx]);

            if (gentex) glTexCoord2f (texcoord [0]+texincr[0], texcoord [1]+texincr[1]);      
            glVertex3fv (b->vtx [ivtx]);

            previvtx = ivtx;
            if (gentex) texcoord [1] += texincr [1];
            glEnd();
         }
      }

      previcross = icross;
      if (gentex) texcoord [0] += texincr [0];
      icross++;
   }
   if (flags&(CLOSE_FIRST|CLOSE_LAST)) {
      GLUtriangulatorObj* obj;      
      Vector normal;
      GLdouble v [3];
      obj = gluNewTess();
      gluTessCallback(obj, GLU_BEGIN, glBegin);
      gluTessCallback(obj, GLU_VERTEX, glVertex3fv); 
      gluTessCallback(obj, GLU_END, glEnd);
      if (flags&CLOSE_FIRST) {
         CrossSection *a = model->cross [0];
         CrossSectionNormal2 (a, normal);
         NormalizeVector (normal);
         glNormal3fv (normal);
         gluBeginPolygon (obj);
         for (ivtx = 0; ivtx < nvtx; ivtx++) {
            v [0] = a->vtx [ivtx][0];
            v [1] = a->vtx [ivtx][1];
            v [2] = a->vtx [ivtx][2];
            gluTessVertex (obj, v, &(a->vtx [ivtx][0]));
         }
         gluEndPolygon (obj);
      }
      if (flags&CLOSE_LAST) {
         CrossSection *a = model->cross [model->ncross-1];
         CrossSectionNormal2 (a, normal);
         NormalizeVector (normal);
         normal [0] = -normal [0];
         normal [1] = -normal [1];
         normal [2] = -normal [2];
         glNormal3fv (normal);
         gluBeginPolygon (obj);
         for (ivtx = nvtx-1; ivtx >= 0; ivtx--) {
            v [0] = a->vtx [ivtx][0];
            v [1] = a->vtx [ivtx][1];
            v [2] = a->vtx [ivtx][2];
            gluTessVertex (obj, v, &(a->vtx [ivtx][0]));
         }
         gluEndPolygon (obj);
      }
      gluDeleteTess (obj);
   }
}

/*--------------------------------------------------------------------------
 *
 *  Main Procedure for generating generic cylinders
 *
 *--------------------------------------------------------------------------*/

int
GenericCylinder (Tcl_Interp *interp, int argc, char* argv [])
{

#define ERRMSG(msg) \
   { Tcl_AppendResult (interp, (msg), (char*) NULL);\
     result = TCL_ERROR;\
     goto done; }

#define ERRMSG2(msg1, msg2) \
   { Tcl_AppendResult (interp, (msg1), (msg2), (char*) NULL);\
     result = TCL_ERROR;\
     goto done; }

#define I { \
      {1, 0, 0, 0}, \
      {0, 1, 0, 0}, \
      {0, 0, 1, 0}, \
      {0, 0, 0, 1}  \
   }

   int result = TCL_OK;

   Matrix transf = I;
   CrossSection *currentCross = NewCrossSection ();
   Model* model = NewModel ();
   int iarg;
   int dlist = -1;

   model->flags = STITCH_LOOPS | SHADE_SMOOTH_ROWS | SHADE_SMOOTH_COLS;
   model->adaptivethreshold = 0;

   AddCrossSectionVtx (currentCross, -1.0, -1.0, 0.0);
   AddCrossSectionVtx (currentCross, -1.0, 1.0, 0.0);
   AddCrossSectionVtx (currentCross, 1.0, 1.0, 0.0);
   AddCrossSectionVtx (currentCross, 1.0, -1.0, 0.0);

   for (iarg = 2; iarg < argc; iarg++) {        
      int len = (int)strlen (argv [iarg]);
      if (strncmp (argv [iarg], "-plot", len) == 0) {
         CrossSection* cross = NewCrossSection ();
         DupCrossSection (currentCross, cross);
         TransformCrossSection (cross, transf);
         AddModelCrossSection (model, cross);
      }
      else if (strncmp (argv [iarg], "-displaylist", len) == 0) {
         iarg++;
         if (strcmp (argv [iarg], "none") == 0) {
            dlist = 0;
         }
         else {
            result = Tcl_GetInt (interp, argv [iarg], &dlist);
            if (result != TCL_OK) goto done;
         }
      }
      else if (strncmp (argv [iarg], "-cross", len) == 0) {
         FreeCrossSection (currentCross);
         currentCross = NewCrossSection ();
         while (iarg+3 < argc && !isalpha (argv [iarg+1][1])) {
            double x, y, z;
            if (Tcl_GetDouble (interp, argv [iarg+1], &x) != TCL_OK ||
                Tcl_GetDouble (interp, argv [iarg+2], &y) != TCL_OK ||
                Tcl_GetDouble (interp, argv [iarg+3], &z) != TCL_OK) {
               ERRMSG ("\n error in -cross");
            }
            AddCrossSectionVtx (currentCross, (GLfloat)x, (GLfloat)y, (GLfloat)z);
            iarg += 3;
         }
      }
      else if (strncmp (argv [iarg], "-polygon", len) == 0) {
         double radius;
         int nsides;
         if (iarg+2 >= argc ||
             Tcl_GetDouble (interp, argv [iarg+1], &radius) != TCL_OK ||
             Tcl_GetInt (interp, argv [iarg+2], &nsides) != TCL_OK) {
            ERRMSG ("\nError in -polygon");
         }
         iarg += 2;
         FreeCrossSection (currentCross);
         currentCross = PolygonCrossSection ((GLfloat)radius, nsides);
      }
      else if (strncmp (argv [iarg], "-stitch", len) == 0) {
         if (iarg+1 >= argc) ERRMSG ("No value for -stitch");
         iarg++;
         len = (int)strlen (argv [iarg]);
         if (strncmp (argv [iarg], "both", len) == 0) {
            model->flags |= STITCH_LOOPS | STITCH_ENDS;
         } else if (strncmp (argv [iarg], "none", len) == 0) {
            model->flags &= ~(STITCH_LOOPS | STITCH_ENDS);
         } else if (strncmp (argv [iarg], "ends", len) == 0) {
            model->flags |= STITCH_ENDS;
            model->flags &= ~STITCH_LOOPS;
         } else if (strncmp (argv [iarg], "loops", len) == 0) {
            model->flags &= ~STITCH_ENDS;
            model->flags |= STITCH_LOOPS;
         } else {
            ERRMSG2 ("Should be 'both', 'none', 'loops' or 'ends':",
                     argv [iarg]);
         }
      }
      else if (strncmp (argv [iarg], "-close", len) == 0) {
         if (iarg+1 >= argc) ERRMSG ("No value for -close");
         iarg++;
         len = (int)strlen (argv [iarg]);
         if (strncmp (argv [iarg], "both", len) == 0) {
            model->flags |= CLOSE_FIRST | CLOSE_LAST;
         } else if (strncmp (argv [iarg], "none", len) == 0) {
            model->flags &= ~(CLOSE_FIRST | CLOSE_LAST);
         } else if (strncmp (argv [iarg], "first", len) == 0) {
            model->flags |= CLOSE_FIRST;
            model->flags &= ~CLOSE_LAST;
         } else if (strncmp (argv [iarg], "last", len) == 0) {
            model->flags &= ~CLOSE_LAST;
            model->flags |= CLOSE_FIRST;
         } else {
            ERRMSG2 ("Should be 'both', 'none', 'first' or 'last':",
                     argv [iarg]);
         }
      }
      else if (strncmp (argv [iarg], "-shade", len) == 0) {
         if (iarg+1 >= argc) ERRMSG ("No value for -shade");
         iarg++;
         len = (int)strlen (argv [iarg]);
         if (strncmp (argv [iarg], "smooth", len) == 0) {
            model->flags |= SHADE_SMOOTH_COLS | SHADE_SMOOTH_ROWS;
         } else if (strncmp (argv [iarg], "flat", len) == 0) {
            model->flags &= ~(SHADE_SMOOTH_COLS | SHADE_SMOOTH_ROWS);
         } else if (strncmp (argv [iarg], "smoothrows", len) == 0) {
            model->flags &= ~SHADE_SMOOTH_COLS;
            model->flags |= SHADE_SMOOTH_ROWS;
         } else if (strncmp (argv [iarg], "smoothcols", len) == 0) {
            model->flags &= ~SHADE_SMOOTH_ROWS;
            model->flags |= SHADE_SMOOTH_COLS;
         } else {
            ERRMSG2 ("Should be 'flat', 'smooth', 'smoothrows' or 'smoothcols':",
                     argv [iarg]);
         }
      }
      else if (strncmp (argv [iarg], "-adaptive", len) == 0) {
         double ang;
         model->flags |= ADAPTIVE;
         if (iarg+1 >= argc ||
             Tcl_GetDouble (interp, argv [iarg+1], &ang) != TCL_OK) 
            ERRMSG ("\nError in -adaptive");
         model->adaptivethreshold = cos(ang*PI/180);
         iarg++;
      }
      else if (strncmp (argv [iarg], "-texgen", len) == 0) {
         double smin, smax, tmin, tmax;
         model->flags |= TEXTURE_GEN;
         if (iarg+4 >= argc ||
             Tcl_GetDouble (interp, argv [iarg+1], &smin) != TCL_OK ||
             Tcl_GetDouble (interp, argv [iarg+2], &smax) != TCL_OK ||
             Tcl_GetDouble (interp, argv [iarg+3], &tmin) != TCL_OK ||
             Tcl_GetDouble (interp, argv [iarg+4], &tmax) != TCL_OK)
            ERRMSG ("\nError in -rotate");
         iarg += 4;
         model->smin = smin;
         model->smax = smax;
         model->tmin = tmin;
         model->tmax = tmax;
      }
      else if (strncmp (argv [iarg], "-identity", len) == 0) {
         Matrix ident = I;
         memcpy (transf, ident, sizeof (Matrix));
      }
      else if (strncmp (argv [iarg], "-rotate", len) == 0) {
         double angle, x, y, z, norm, sint, cost;
         Matrix m = I;
         if (iarg+4 >= argc ||
             Tcl_GetDouble (interp, argv [iarg+1], &angle) != TCL_OK ||
             Tcl_GetDouble (interp, argv [iarg+2], &x) != TCL_OK ||
             Tcl_GetDouble (interp, argv [iarg+3], &y) != TCL_OK ||
             Tcl_GetDouble (interp, argv [iarg+4], &z) != TCL_OK)
            ERRMSG ("\nError in -rotate");
         iarg += 4;
         norm = sqrt (x*x+y*y+z*z);
         if (norm == 0.0) ERRMSG ("Null Vector");
         x /= norm; y /= norm; z /= norm;
         angle *= PI/180;
         sint = sin(angle);
         cost = cos(angle);
         m [0][0] = (GLfloat)(x*x + cost * (1 - x*x) + sint * 0); 
         m [0][1] = (GLfloat)(x*y + cost * (0 - x*y) + sint * (-z)); 
         m [0][2] = (GLfloat)(x*z + cost * (0 - x*z) + sint * (y));
         m [1][0] = (GLfloat)(y*x + cost * (0 - y*x) + sint * (z));
         m [1][1] = (GLfloat)(y*y + cost * (1 - y*y) + sint * 0);
         m [1][2] = (GLfloat)(y*z + cost * (0 - y*z) + sint * (-x));
         m [2][0] = (GLfloat)(z*x + cost * (0 - z*x) + sint * (-y));
         m [2][1] = (GLfloat)(z*y + cost * (0 - z*y) + sint * (x));
         m [2][2] = (GLfloat)(z*z + cost * (1 - z*z) + sint * 0); 
         MatrixMult (m, transf, transf);
      }
      else if (strncmp (argv [iarg], "-translate", len) == 0) {
         double x, y, z;
         Matrix m = I;
         if (iarg+3 >= argc ||
             Tcl_GetDouble (interp, argv [iarg+1], &x) != TCL_OK ||
             Tcl_GetDouble (interp, argv [iarg+2], &y) != TCL_OK ||
             Tcl_GetDouble (interp, argv [iarg+3], &z) != TCL_OK)
            ERRMSG ("\nError in -translate");
         iarg += 3;
         m [0][3] = (GLfloat)x;  m [1][3] = (GLfloat)y; m [2][3] = (GLfloat)z;
         MatrixMult (m, transf, transf);
      }
      else if (strncmp (argv [iarg], "-scale", len) == 0) {
         double x, y, z;
         Matrix m = I;
         if (iarg+3 >= argc ||
             Tcl_GetDouble (interp, argv [iarg+1], &x) != TCL_OK ||
             Tcl_GetDouble (interp, argv [iarg+2], &y) != TCL_OK ||
             Tcl_GetDouble (interp, argv [iarg+3], &z) != TCL_OK)
            ERRMSG ("\nError in -scale");
         iarg += 3;
         m [0][0] = (GLfloat)x;  m [1][1] = (GLfloat)y; m [2][2] = (GLfloat)z;
         MatrixMult (m, transf, transf);
      }
      else {
         ERRMSG2 ("No such option: ", argv [iarg]);
      }
   }

done:
   if (result == TCL_OK) {
      char buf [100];
      if (dlist == -1) dlist = glGenLists (1);
      if (dlist != 0) glNewList (dlist, GL_COMPILE);
      UniformCrossSectionLengths (model);
      ComputeModelNormals (model);
      RenderModel (model);
      if (dlist != 0) {
         glEndList ();
         sprintf (buf, "displaylist %d", dlist);
      }
      Tcl_AppendElement (interp, buf);
      ComputeModelBBox (interp, model);
   }
   FreeModel (model);
   return result;
}




            
            




