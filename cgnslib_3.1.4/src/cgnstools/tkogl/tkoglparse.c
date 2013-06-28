#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "tkogl.h"
#include "tkoglparse.h"

#define ERRMSG(msg) { Tcl_AppendResult(interp,msg,(char*)NULL); \
		      result = TCL_ERROR; goto done; }

#define ERRMSG2(msg1,msg2) { Tcl_AppendResult(interp,msg1,msg2,(char*)NULL); \
			     result = TCL_ERROR; goto done; }

#define MAXARGS 200

/*---------------------------------------------------------------------------
 *
 * The following is a table that implements the binding between
 * GL enum constants and their tcl names
 *
 *---------------------------------------------------------------------------*/

typedef struct {
   char * name;
   GLenum code;
} parseItem;

static parseItem enumTable [] = {
   {"accum", GL_ACCUM},
   {"accumbuffer", GL_ACCUM_BUFFER_BIT},
   {"accumbufferbit", GL_ACCUM_BUFFER_BIT},
   {"add", GL_ADD},
   {"alphatest", GL_ALPHA_TEST},
   {"always", GL_ALWAYS},
   {"allattrib", GL_ALL_ATTRIB_BITS},
   {"allattribbits", GL_ALL_ATTRIB_BITS},
   {"ambient", GL_AMBIENT},
   {"ambientanddiffuse", GL_AMBIENT_AND_DIFFUSE},
   {"autonormal", GL_AUTO_NORMAL},
   {"aux0", GL_AUX0},
   {"aux1", GL_AUX1},
   {"aux2", GL_AUX2},
   {"aux3", GL_AUX3},
   {"back", GL_BACK},
   {"backleft", GL_BACK_LEFT},
   {"backright", GL_BACK_RIGHT},
   {"blend", GL_BLEND},
   {"bluebias", GL_BLUE_BIAS},
   {"bluescale", GL_BLUE_SCALE},
   {"ccw", GL_CCW},
   {"clamp", GL_CLAMP},
   {"clipplane0", GL_CLIP_PLANE0},
   {"clipplane1", GL_CLIP_PLANE1},
   {"clipplane2", GL_CLIP_PLANE2},
   {"clipplane3", GL_CLIP_PLANE3},
   {"clipplane4", GL_CLIP_PLANE4},
   {"clipplane5", GL_CLIP_PLANE5},
   {"colorbuffer", GL_COLOR_BUFFER_BIT},
   {"colorbufferbit", GL_COLOR_BUFFER_BIT},
   {"colorindex", GL_COLOR_INDEX},
   {"colormaterial", GL_COLOR_MATERIAL},
   {"compile", GL_COMPILE},
   {"compileandaxecute", GL_COMPILE_AND_EXECUTE},
   {"constantattenuation", GL_CONSTANT_ATTENUATION},
   {"cullface", GL_CULL_FACE},
   {"current", GL_CURRENT_BIT},
   {"currentbit", GL_CURRENT_BIT},
   {"cw", GL_CW},
   {"decal", GL_DECAL},
   {"decr", GL_DECR},
   {"depthbuffer", GL_DEPTH_BUFFER_BIT},
   {"depthbufferbit", GL_DEPTH_BUFFER_BIT},
   {"depthtest", GL_DEPTH_TEST},
   {"diffuse", GL_DIFFUSE},
   {"dither", GL_DITHER},
   {"dstalpha", GL_DST_ALPHA},
   {"dstcolor", GL_DST_COLOR},
   {"enable", GL_ENABLE_BIT},
   {"enablebit", GL_ENABLE_BIT},
   {"emission", GL_EMISSION},
   {"equal", GL_EQUAL},
   {"eval", GL_EVAL_BIT},
   {"evalbit", GL_EVAL_BIT},
   {"exp", GL_EXP},
   {"exp2", GL_EXP2},
   {"eyelinear", GL_EYE_LINEAR},
   {"eyeplane", GL_EYE_PLANE},
   {"fill", GL_FILL},
   {"flat", GL_FLAT},
   {"fog", GL_FOG},
   {"fogbit", GL_FOG_BIT},
   {"fogcolor", GL_FOG_COLOR},
   {"fogdensity", GL_FOG_DENSITY},
   {"fogend", GL_FOG_END},
   {"fogmode", GL_FOG_MODE},
   {"fogstart", GL_FOG_START},
   {"front", GL_FRONT},
   {"frontandback", GL_FRONT_AND_BACK},
   {"frontleft", GL_FRONT_LEFT},
   {"frontright", GL_FRONT_RIGHT},
   {"gequal", GL_GEQUAL},
   {"greater", GL_GREATER},
   {"greenbias", GL_GREEN_BIAS},
   {"greenscale", GL_GREEN_SCALE},
   {"hint", GL_HINT_BIT},
   {"hintbit", GL_HINT_BIT},
   {"incr", GL_INCR},
   {"invert", GL_INVERT},
   {"keep", GL_KEEP},
   {"left", GL_LEFT},
   {"lequal", GL_LEQUAL},
   {"less", GL_LESS},
   {"light0", GL_LIGHT0},
   {"light1", GL_LIGHT1},
   {"light2", GL_LIGHT2},
   {"light3", GL_LIGHT3},
   {"light4", GL_LIGHT4},
   {"light5", GL_LIGHT5},
   {"light6", GL_LIGHT6},
   {"light7", GL_LIGHT7},
   {"lighting", GL_LIGHTING},
   {"lightingbit", GL_LIGHTING_BIT},
   {"lightmodelambient", GL_LIGHT_MODEL_AMBIENT},
   {"lightmodellocalviewer", GL_LIGHT_MODEL_LOCAL_VIEWER},
   {"lightmodeltwoside", GL_LIGHT_MODEL_TWO_SIDE},
   {"line", GL_LINE},
   {"linebit", GL_LINE_BIT},
   {"linear", GL_LINEAR},
   {"linearattenuation", GL_LINEAR_ATTENUATION},
   {"lineloop", GL_LINE_LOOP},
   {"lines", GL_LINES},
   {"linesmooth", GL_LINE_SMOOTH},
   {"linestipple", GL_LINE_STIPPLE},
   {"linestrip", GL_LINE_STRIP},
   {"list", GL_LIST_BIT},
   {"listbit", GL_LIST_BIT},
   {"load", GL_LOAD},
   {"map1color4", GL_MAP1_COLOR_4},
   {"map1normal", GL_MAP1_NORMAL},
   {"map1texturecoord1", GL_MAP1_TEXTURE_COORD_1},
   {"map1texturecoord2", GL_MAP1_TEXTURE_COORD_2},
   {"map1texturecoord3", GL_MAP1_TEXTURE_COORD_3},
   {"map1texturecoord4", GL_MAP1_TEXTURE_COORD_4},
   {"map1vertex3", GL_MAP1_VERTEX_3},
   {"map1vertex4", GL_MAP1_VERTEX_4},
   {"map2color4", GL_MAP2_COLOR_4},
   {"map2normal", GL_MAP2_NORMAL},
   {"map2texturecoord1", GL_MAP2_TEXTURE_COORD_1},
   {"map2texturecoord2", GL_MAP2_TEXTURE_COORD_2},
   {"map2texturecoord3", GL_MAP2_TEXTURE_COORD_3},
   {"map2texturecoord4", GL_MAP2_TEXTURE_COORD_4},
   {"map2vertex3", GL_MAP2_VERTEX_3},
   {"map2vertex4", GL_MAP2_VERTEX_4},
   {"modelview", GL_MODELVIEW},
   {"modulate", GL_MODULATE},
   {"mult", GL_MULT},
   {"nearest", GL_NEAREST},
   {"never", GL_NEVER},
   {"none", GL_NONE},
   {"normalize", GL_NORMALIZE},
   {"notequal", GL_NOTEQUAL},
   {"objectlinear", GL_OBJECT_LINEAR},
   {"objectplane", GL_OBJECT_PLANE},
   {"one", GL_ONE},
   {"oneminusdstalpha", GL_ONE_MINUS_DST_ALPHA},
   {"oneminusdstcolor", GL_ONE_MINUS_DST_COLOR},
   {"oneminussrcalpha", GL_ONE_MINUS_SRC_ALPHA},
   {"oneminussrccolor", GL_ONE_MINUS_SRC_COLOR},
   {"packalignment", GL_PACK_ALIGNMENT},
   {"packlsbfirst", GL_PACK_LSB_FIRST},
   {"packrowlength", GL_PACK_ROW_LENGTH},
   {"packskippixels", GL_PACK_SKIP_PIXELS},
   {"packskiprows", GL_PACK_SKIP_ROWS},
   {"packswapbytes", GL_PACK_SWAP_BYTES},
   {"pixelmode", GL_PIXEL_MODE_BIT},
   {"pixelmodebit", GL_PIXEL_MODE_BIT},
   {"point", GL_POINT},
   {"pointbit", GL_POINT_BIT},
   {"points", GL_POINTS},
   {"polygon", GL_POLYGON},
   {"polygonbit", GL_POLYGON_BIT},
   {"polygonstipple", GL_POLYGON_STIPPLE_BIT},
   {"polygonstipplebit", GL_POLYGON_STIPPLE_BIT},
   {"position", GL_POSITION},
   {"projection", GL_PROJECTION},
   {"q", GL_Q},
   {"quadraticattenuation", GL_QUADRATIC_ATTENUATION},
   {"quads", GL_QUADS},
   {"quadstrip", GL_QUAD_STRIP},
   {"r", GL_R},
   {"redbias", GL_RED_BIAS},
   {"redscale", GL_RED_SCALE},
   {"repeat", GL_REPEAT},
   {"replace", GL_REPLACE},
   {"return", GL_RETURN},
   {"right", GL_RIGHT},
   {"s", GL_S},
   {"scissor", GL_SCISSOR_TEST},
   {"scissorbit", GL_SCISSOR_BIT},
   {"shininess", GL_SHININESS},
   {"smooth", GL_SMOOTH},
   {"specular", GL_SPECULAR},
   {"spheremap", GL_SPHERE_MAP},
   {"spotcutoff", GL_SPOT_CUTOFF},
   {"spotdirecion", GL_SPOT_DIRECTION},
   {"spotexponent", GL_SPOT_EXPONENT},
   {"srcalpha", GL_SRC_ALPHA},
   {"srcalphasaturate", GL_SRC_ALPHA_SATURATE},
   {"srccolor", GL_SRC_COLOR},
   {"stenciltest", GL_STENCIL_TEST},
   {"stencilbuffer", GL_STENCIL_BUFFER_BIT},
   {"stencilbufferbit", GL_STENCIL_BUFFER_BIT},
   {"t", GL_T},
   {"texture", GL_TEXTURE},
   {"texture1d", GL_TEXTURE_1D},
   {"texture2d", GL_TEXTURE_2D},
   {"texturebit", GL_TEXTURE_BIT},
   {"texturebordercolor", GL_TEXTURE_BORDER_COLOR},
   {"textureenv", GL_TEXTURE_ENV},
   {"textureenvcolor", GL_TEXTURE_ENV_COLOR},
   {"textureenvmode", GL_TEXTURE_ENV_MODE},
   {"texturegenmode", GL_TEXTURE_GEN_MODE},
   {"texturegens", GL_TEXTURE_GEN_S},
   {"texturegent", GL_TEXTURE_GEN_T},
   {"texturemagfilter", GL_TEXTURE_MAG_FILTER},
   {"textureminfilter", GL_TEXTURE_MIN_FILTER},
   {"texturewraps", GL_TEXTURE_WRAP_S},
   {"texturewrapt", GL_TEXTURE_WRAP_T},
   {"transform", GL_TRANSFORM_BIT},
   {"transformbit", GL_TRANSFORM_BIT},
   {"triangles", GL_TRIANGLES},
   {"trianglefan", GL_TRIANGLE_FAN},
   {"trianglestrip", GL_TRIANGLE_STRIP},
   {"unpackalignment", GL_UNPACK_ALIGNMENT},
   {"unpacklsbfirst", GL_UNPACK_LSB_FIRST},
   {"unpackrowlength", GL_UNPACK_ROW_LENGTH},
   {"unpackskippixels", GL_UNPACK_SKIP_PIXELS},
   {"unpackskiprows", GL_UNPACK_SKIP_ROWS},
   {"unpackswapbytes", GL_UNPACK_SWAP_BYTES},
   {"viewport", GL_VIEWPORT_BIT},
   {"viewportbit", GL_VIEWPORT_BIT},
   {"zero", GL_ZERO}
};

/*---------------------------------------------------------------------------
 *
 * The following are tables of enum constants that are required as
 * arguments for various OpenGL functions 
 *
 *---------------------------------------------------------------------------*/

static GLenum primTable [] = {
   GL_POINTS, GL_LINES, GL_LINE_STRIP, GL_LINE_LOOP, GL_POLYGON,
   GL_QUADS, GL_QUAD_STRIP, GL_TRIANGLES, GL_TRIANGLE_STRIP,
   GL_TRIANGLE_FAN
};

/* Names for capabilities set/reset by glEnable/glDisable */
static GLenum capabTable [] = { 
   GL_ALPHA_TEST, GL_AUTO_NORMAL, GL_BLEND, GL_COLOR_MATERIAL,
   GL_CLIP_PLANE0, GL_CLIP_PLANE1, GL_CLIP_PLANE2, GL_CLIP_PLANE3,
   GL_CLIP_PLANE4, GL_CLIP_PLANE5, GL_CULL_FACE, GL_DEPTH_TEST,
   GL_DITHER, GL_FOG, GL_LINE_SMOOTH, GL_LIGHT0, GL_LIGHT1, GL_LIGHT2,
   GL_LIGHT3, GL_LIGHT4, GL_LIGHT5, GL_LIGHT6, GL_LIGHT7, GL_LIGHTING,
   GL_LINE_STIPPLE, GL_MAP1_VERTEX_3, GL_MAP1_VERTEX_4, GL_MAP1_COLOR_4,
   GL_MAP1_NORMAL, GL_MAP1_TEXTURE_COORD_1, GL_MAP1_TEXTURE_COORD_2,
   GL_MAP1_TEXTURE_COORD_3, GL_MAP1_TEXTURE_COORD_4, GL_MAP2_VERTEX_3,
   GL_MAP2_VERTEX_4, GL_MAP2_COLOR_4, GL_MAP2_NORMAL,
   GL_MAP2_TEXTURE_COORD_1, GL_MAP2_TEXTURE_COORD_2,
   GL_MAP2_TEXTURE_COORD_3, GL_MAP2_TEXTURE_COORD_4, GL_NORMALIZE,
   GL_SCISSOR_TEST, GL_STENCIL_TEST,
   GL_TEXTURE_1D, GL_TEXTURE_2D, GL_TEXTURE_GEN_S, GL_TEXTURE_GEN_T 
};

/* Names for glMatrixMode */
static GLenum matrixModeTable [] = {
   GL_MODELVIEW, GL_PROJECTION, GL_TEXTURE,
};
  

/* Names describing face sides (for glPolygonMode and Material) */
static GLenum faceTable [] = {
   GL_BACK, GL_FRONT, GL_FRONT_AND_BACK,
};

/* Names describing face sides (for glCullFace) */
static GLenum cullFaceTable [] = {
   GL_BACK, GL_FRONT
};

/* Names for glPolygonMode */
static GLenum pModeTable [] = {
   GL_FILL, GL_LINE, GL_POINT
};

/* Names for glLightModel */
static GLenum lightModelTable [] = {
   GL_LIGHT_MODEL_AMBIENT, GL_LIGHT_MODEL_LOCAL_VIEWER,
   GL_LIGHT_MODEL_TWO_SIDE
};

/* Names for Material attribute (glMaterial) */
static GLenum matParmTable[] = {
   GL_AMBIENT, GL_SPECULAR, GL_DIFFUSE, GL_EMISSION, GL_SHININESS,
   GL_AMBIENT_AND_DIFFUSE, GL_COLOR_INDEX
};

/* Names for mode of glColorMaterial */
static GLenum modeColorMatTable[] = {
   GL_AMBIENT, GL_SPECULAR, GL_DIFFUSE, GL_EMISSION, 
   GL_AMBIENT_AND_DIFFUSE
};

/* Names for Light numbers (glLight) */
static GLenum lightNumberTable[] = {
   GL_LIGHT0, GL_LIGHT1, GL_LIGHT2, GL_LIGHT3, GL_LIGHT4, GL_LIGHT5,
   GL_LIGHT6, GL_LIGHT7
};


/* Names for Light attributes (glLight) */
static GLenum lightParmTable[] = {
   GL_AMBIENT, GL_CONSTANT_ATTENUATION, GL_DIFFUSE,
   GL_LINEAR_ATTENUATION, GL_POSITION, GL_QUADRATIC_ATTENUATION,
   GL_SPECULAR, GL_SPOT_CUTOFF, GL_SPOT_DIRECTION, GL_SPOT_EXPONENT
};

/* Names for glShadeModel */
static GLenum shadeModelTable [] = {
   GL_FLAT, GL_SMOOTH
};

/* Names for glClear */
static GLenum bufferTable [] = {
   GL_ACCUM_BUFFER_BIT, GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT,
   GL_STENCIL_BUFFER_BIT
};

/* Names for glPixelStore */
#if 0
static GLenum pixelStoreTable [] = {
   GL_PACK_ALIGNMENT, GL_PACK_LSB_FIRST, GL_PACK_ROW_LENGTH,
   GL_PACK_SKIP_PIXELS, GL_PACK_SKIP_ROWS, GL_PACK_SWAP_BYTES,
   GL_UNPACK_ALIGNMENT, GL_UNPACK_LSB_FIRST,
   GL_UNPACK_ROW_LENGTH, GL_UNPACK_SKIP_PIXELS,
   GL_UNPACK_SKIP_ROWS, GL_UNPACK_SWAP_BYTES
};
#endif

/* Param names for glPixelTransfer */
static GLenum pixelTransferTable [] = {
   GL_RED_SCALE, GL_GREEN_SCALE, GL_BLUE_SCALE, GL_RED_BIAS,
   GL_GREEN_BIAS, GL_BLUE_BIAS
};

/* Target names for glTexParameter */
static GLenum targetTexParamTable [] = {
   GL_TEXTURE_1D, GL_TEXTURE_2D
};

/* pnames for glTexParameter */
static GLenum nameTexParamTable [] = {
   GL_TEXTURE_BORDER_COLOR, GL_TEXTURE_MAG_FILTER,
   GL_TEXTURE_MIN_FILTER, GL_TEXTURE_WRAP_S, GL_TEXTURE_WRAP_T
};


/* values for glTexParameter */
static GLenum valueTexParamTable [] = {
   GL_NEAREST, GL_LINEAR, GL_CLAMP, GL_REPEAT
};

/* target names for glTexEnv */
static GLenum targetTexEnvTable [] = {
   GL_TEXTURE_ENV
};

/* pnames for glTexEnv */
static GLenum nameTexEnvTable [] = {
   GL_TEXTURE_ENV_MODE, GL_TEXTURE_ENV_COLOR
};

/* value names for glTexEnv */
static GLenum valueTexEnvTable [] = {
   GL_MODULATE, GL_DECAL, GL_BLEND
};


/* coord names for glTexGen */
static GLenum coordTexGenTable [] = {
   GL_S, GL_T, GL_R, GL_Q
};

/* pnames for glTexGen */ 
static GLenum nameTexGenTable [] = {
   GL_TEXTURE_GEN_MODE, GL_OBJECT_PLANE, GL_EYE_PLANE
};

/* values for glTexGen */
static GLenum valueTexGenTable [] = {
   GL_OBJECT_LINEAR, GL_EYE_LINEAR, GL_SPHERE_MAP
};

/* target names for glMap1 */
static GLenum targetMap1Table [] = {
   GL_MAP1_VERTEX_3, GL_MAP1_VERTEX_4, GL_MAP1_COLOR_4,
   GL_MAP1_NORMAL, GL_MAP1_TEXTURE_COORD_1, GL_MAP1_TEXTURE_COORD_2,
   GL_MAP1_TEXTURE_COORD_3, GL_MAP1_TEXTURE_COORD_4

};

/* target names for glMap2 */
static GLenum targetMap2Table [] = {
   GL_MAP2_VERTEX_3, GL_MAP2_VERTEX_4, GL_MAP2_COLOR_4,
   GL_MAP2_NORMAL, GL_MAP2_TEXTURE_COORD_1, GL_MAP2_TEXTURE_COORD_2,
   GL_MAP2_TEXTURE_COORD_3, GL_MAP2_TEXTURE_COORD_4
};

/* Names for glEvalMesh1 */
static GLenum modeMesh1Table [] = {
   GL_LINE, GL_POINT
};

/* Names for glEvalMesh2 */
static GLenum modeMesh2Table [] = {
   GL_FILL, GL_LINE, GL_POINT
};

/* Names for glFog */
static GLenum nameFogTable [] = {
   GL_FOG_MODE, GL_FOG_DENSITY, GL_FOG_START, GL_FOG_END, GL_FOG_COLOR
};

/* Mode names for glFog (GL_FOG_MODE ...) */
static GLenum fogFogModeTable [] = {
   GL_EXP, GL_EXP2, GL_LINEAR
};

/* Names for glFrontFace */
static GLenum nameFrontFaceTable [] = {
   GL_CW, GL_CCW
};

/* mode names for glDrawBuffer */
static GLenum modeDrawBufferTable [] = {
   GL_FRONT, GL_BACK, GL_RIGHT, GL_LEFT, GL_FRONT_RIGHT, GL_FRONT_LEFT,
   GL_BACK_RIGHT, GL_BACK_LEFT, GL_FRONT_AND_BACK, 
   GL_AUX0, GL_AUX1, GL_AUX2, GL_AUX3, GL_NONE
};

/* mode names for glReadBuffer */
static GLenum modeReadBufferTable [] = {
   GL_FRONT, GL_BACK, GL_RIGHT, GL_LEFT, GL_FRONT_RIGHT, GL_FRONT_LEFT,
   GL_BACK_RIGHT, GL_BACK_LEFT, 
   GL_AUX0, GL_AUX1, GL_AUX2, GL_AUX3
};

/* func names for glAlphaFunc, glDepthFunc and glStencilFunc */
static GLenum funcAlphaStencilTable [] = {
   GL_NEVER, GL_ALWAYS, GL_LESS, GL_LEQUAL, GL_EQUAL, GL_GEQUAL, GL_GREATER,
   GL_NOTEQUAL
};

/* op names for glStencilOp */
static GLenum opStencilTable [] = {
   GL_KEEP, GL_ZERO, GL_REPLACE, GL_INCR, GL_DECR, GL_INVERT
};

/* op names for glAccum */
static GLenum opAccumTable [] = {
   GL_ACCUM, GL_LOAD, GL_ADD, GL_MULT, GL_RETURN
};


/* sfactor names for glBlendFunc */
static GLenum sfactorBlendTable [] = {
   GL_ZERO, GL_ONE, GL_DST_COLOR, 
   GL_ONE_MINUS_DST_COLOR,  GL_SRC_ALPHA,
   GL_ONE_MINUS_SRC_ALPHA, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA,
   GL_SRC_ALPHA_SATURATE
};

/* dfactor names for glBlendFunc */
static GLenum dfactorBlendTable [] = {
   GL_ZERO, GL_ONE, GL_SRC_COLOR,
   GL_ONE_MINUS_SRC_COLOR, GL_SRC_ALPHA,
   GL_ONE_MINUS_SRC_ALPHA, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA,
};

/* clipplane names for glClipPlane */
static GLenum clipPlaneTable [] = {
   GL_CLIP_PLANE0, GL_CLIP_PLANE1,  GL_CLIP_PLANE2, 
   GL_CLIP_PLANE3, GL_CLIP_PLANE4,  GL_CLIP_PLANE5
};
 
/* bit names for glPushAttrib */
static GLenum attribNameTable [] = {
   GL_ACCUM_BUFFER_BIT, GL_ALL_ATTRIB_BITS, GL_COLOR_BUFFER_BIT,
   GL_CURRENT_BIT, GL_DEPTH_BUFFER_BIT, GL_ENABLE_BIT, GL_EVAL_BIT,
   GL_FOG_BIT, GL_HINT_BIT, GL_LIGHTING_BIT, GL_LINE_BIT, GL_LIST_BIT,
   GL_PIXEL_MODE_BIT, GL_POINT_BIT, GL_POLYGON_BIT,
   GL_POLYGON_STIPPLE_BIT, GL_SCISSOR_BIT, GL_STENCIL_BUFFER_BIT,
   GL_TEXTURE_BIT, GL_TRANSFORM_BIT, GL_VIEWPORT_BIT
};

/*  names for glNewList */
static GLenum newListTable [] = {
   GL_COMPILE, GL_COMPILE_AND_EXECUTE
};

/*  target names for glHint */
static GLenum hintTargetTable [] = {
   GL_FOG_HINT, GL_LINE_SMOOTH_HINT, GL_PERSPECTIVE_CORRECTION_HINT, 
   GL_POINT_SMOOTH_HINT, GL_POLYGON_SMOOTH_HINT
};

/*  mode names for glHint */
static GLenum hintModeTable [] = {
   GL_FASTEST, GL_NICEST, GL_DONT_CARE
};


/*---------------------------------------------------------------------------
 *
 *  The following data structures describe what arguments each GL function
 *  expects.
 *
 *---------------------------------------------------------------------------*/

typedef enum {
   ARG_ENUM,	  /* Name of parse table follows */
   ARG_BIT_FIELD, /* Name of parse table follows */
   ARG_FLOAT, 	  /* Number of floats follow */
   ARG_INT, 	  /* Number of ints follow */
   ARG_BOOLEAN,	  /* Number of booleans follow */
   ARG_VAR_FLOAT,  /* Min # args, max # args, defaults follow */
   ARG_STRING,      /* A char string */
   ARG_DONT_CARE    /* All remaining args are unprocessed initially */
} ArgType;

typedef struct ArgDetailStruct {
   ArgType type;
   union {
      struct {
	 GLenum * enumTable;
	 int enumSize;
      } enumDetail;
      int nValues;
      struct {
	 int minArgs, maxArgs;
	 GLfloat def [3];
      } varFloatDetail;
   } detail;
   struct ArgDetailStruct * next;
} * ArgDetailList, ArgDetail;



/*----------------------------------------------------------------------------
 * 
 * Each of the following procs generates an appropriate Arg Detail description
 *
 *---------------------------------------------------------------------------*/

#define ArgEnum(enumTable,next) \
   ArgEnumTable (enumTable, sizeof(enumTable)/sizeof(GLenum), ARG_ENUM, next)

#define ArgBitField(enumTable,next) \
   ArgEnumTable (enumTable, sizeof(enumTable)/sizeof(GLenum), ARG_BIT_FIELD,\
		 next)

static ArgDetailList 
ArgEnumTable (enumTable, enumTableSize, type, next) 
   GLenum * enumTable;
   int enumTableSize;
   ArgDetailList next;
   ArgType type;
{
   ArgDetailList arg = (ArgDetailList) malloc (sizeof (ArgDetail));
   assert (arg != NULL);
   arg->type = type;
   arg->detail.enumDetail.enumTable = enumTable;
   arg->detail.enumDetail.enumSize = enumTableSize;
   arg->next = next;
   return arg;
}

static ArgDetailList
ArgFloat (nValues, next)
     int nValues;
     ArgDetailList next;
{
   ArgDetailList arg = (ArgDetailList) malloc (sizeof (ArgDetail));
   assert (arg != NULL);
   arg->type = ARG_FLOAT;
   arg->detail.nValues = nValues;
   arg->next = next;
   return arg;
}

static ArgDetailList
ArgInt (nValues, next)
     int nValues;
     ArgDetailList next;
{
   ArgDetailList arg = (ArgDetailList) malloc (sizeof (ArgDetail));
   assert (arg != NULL);
   arg->type = ARG_INT;
   arg->detail.nValues = nValues;
   arg->next = next;
   return arg;
}

static ArgDetailList
ArgBoolean (nValues, next)
     int nValues;
     ArgDetailList next;
{
   ArgDetailList arg = (ArgDetailList) malloc (sizeof (ArgDetail));
   assert (arg != NULL);
   arg->type = ARG_BOOLEAN;
   arg->detail.nValues = nValues;
   arg->next = next;
   return arg;
}

static ArgDetailList
ArgString (next)
     ArgDetailList next;
{
   ArgDetailList arg = (ArgDetailList) malloc (sizeof (ArgDetail));
   assert (arg != NULL);
   arg->type = ARG_STRING;
   arg->next = next;
   return arg;
}

static ArgDetailList
ArgDontCare (next)
     ArgDetailList next;
{
   ArgDetailList arg = (ArgDetailList) malloc (sizeof (ArgDetail));
   assert (arg != NULL);
   arg->type = ARG_DONT_CARE;
   arg->next = next;
   return arg;
}

static ArgDetailList
ArgVarFloat (minArgs, maxArgs, def0, def1, def2, next)
     int minArgs, maxArgs;
     GLfloat def0, def1, def2;
     ArgDetailList next;
{
   ArgDetailList arg = (ArgDetailList) malloc (sizeof (ArgDetail));
   assert (arg != NULL);
   assert (maxArgs - minArgs <= 3);
   arg->type = ARG_VAR_FLOAT;
   arg->detail.varFloatDetail.minArgs = minArgs;
   arg->detail.varFloatDetail.maxArgs = maxArgs;
   arg->detail.varFloatDetail.def[0] = (GLfloat)def0;
   arg->detail.varFloatDetail.def[1] = (GLfloat)def1;
   arg->detail.varFloatDetail.def[2] = (GLfloat)def2;
   arg->next = next;
   return arg;
}

/*---------------------------------------------------------------------------
 *
 * Data structures that define how to parse command lines that generate
 * calls to GL functions. Each GL function has an entry in funcDescTable
 * which contains the corresponding tcl binding (that has to be preceded
 * by a dash in the style of tcl options),  an encoded description of its 
 * argument list and a function which is to be called if the parse is 
 * successful
 *
 *---------------------------------------------------------------------------*/

typedef int TkOGLFunc _ANSI_ARGS_((Tcl_Interp* interp, void** args, int nargs));

static TkOGLFunc 
    TkAccum, 	  
    TkAlphaFunc,   
    TkBegin,        
    TkBlendFunc,
    TkCall,        
    TkClear,
    TkClearAccum,   
    TkClearColor,  
    TkClearDepth,  
    TkClearStencil, 
    TkClipPlane,   
    TkColor,
    TkColorMask,   
    TkColorMaterial,
    TkCopyPixels,
    TkCullFace,
    TkDepthFunc,
    TkDepthMask,
    TkDisable,     
    TkDrawBuffer,     
    TkDrawPixels,
    TkEdgeFlag,
    TkEnable,      
    TkEnd,          
    TkEndList,
    TkEvalCoord1,
    TkEvalCoord2,
    TkEvalMesh1,   
    TkEvalMesh2,
    TkFlush,
    TkFinish,
    TkFog, 	 
    TkFrontFace, 	  
    TkFrustum,      
    TkHint, 
    TkInitNames,
    TkLight,       
    TkLightModel,   
    TkLineStipple, 
    TkLineWidth,
    TkLoadIdentity, 
    TkLoadMatrix,  
    TkLoadName, 
    TkLookAt,
    TkMap1,      
    TkMap2,      
    TkMapGrid1,  
    TkMapGrid2,
    TkMaterial,    
    TkMatrixMode,
    TkMultMatrix,  
    TkNewList, 	  
    TkNormal,      
    TkOrtho,        
    TkPerspective,
    TkPickMatrix,  
    TkPixelTransfer, 
    TkPixelZoom,     
    TkPointSize, 	  
    TkPolygonMode, 
    TkPopAttrib, 	 
    TkPopMatrix,   
    TkPopName,  
    TkPushAttrib,  
    TkPushMatrix,
    TkPushName,   
    TkRasterPos,
    TkReadBuffer,   
    TkReadPixels, 
    TkRect,  
    TkRotate,      
    TkScale,
    TkScissor,
    TkShadeModel,  
    TkStencilFunc, 
    TkStencilMask,	 
    TkStencilOp,    
    TkTexCoord,      
    TkTexEnv,        
    TkTexGen,        
    TkTexImage1D,    
    TkTexImage2D,    
    TkTexParameter,  
    TkTranslate,   
    TkVertex;

    
typedef struct funcDescStruct {
   TkOGLFunc *func;
   ArgDetailList argList;
} FuncDesc;

/*---------------------------------------------------------------------------
 *
 *  Initialization of the Hash Tables that govern the parsing 
 *  subsystem
 *
 *---------------------------------------------------------------------------*/

Tcl_HashTable 
   funcDescHashTable,
   enumHashTable;

static void
InsDesc (tclName, func, argList)
     char * tclName;
     TkOGLFunc *func;
     ArgDetailList argList;
{
   Tcl_HashEntry *descEntry;
   int newEntry = 0;
   FuncDesc * where;
   descEntry = Tcl_CreateHashEntry (&funcDescHashTable, tclName, &newEntry);
   assert (newEntry);
   where = (FuncDesc *) malloc (sizeof (FuncDesc));
   assert (where!=NULL);
   (where)->func = func;
   (where)->argList = argList;
   Tcl_SetHashValue (descEntry, (ClientData) where);
}

void
InitHashTables () 
{
   int i, newEntry;
   Tcl_HashEntry *enumEntry;
   Tcl_InitHashTable (&funcDescHashTable, TCL_STRING_KEYS);
   Tcl_InitHashTable (&enumHashTable, TCL_STRING_KEYS);

   for (i = 0; i < sizeof (enumTable)/sizeof(parseItem); i++) {
      enumEntry = Tcl_CreateHashEntry (&enumHashTable, 
				       enumTable [i].name, &newEntry);
      Tcl_SetHashValue (enumEntry, (ClientData)((size_t)(enumTable [i].code)));
      assert (newEntry);
   }	       

   InsDesc ("matrixmode" , TkMatrixMode,  ArgEnum(matrixModeTable, NULL));
   InsDesc ("pushmatrix" , TkPushMatrix,  NULL);
   InsDesc ("popmatrix"  , TkPopMatrix,   NULL);
   InsDesc ("loadidentity",TkLoadIdentity,NULL);
   InsDesc ("multmatrix" , TkMultMatrix,  ArgFloat(16, NULL));
   InsDesc ("loadmatrix" , TkLoadMatrix,  ArgFloat(16, NULL));
   InsDesc ("translate"  , TkTranslate,   ArgFloat(3, NULL));
   InsDesc ("rotate"     , TkRotate,      ArgFloat(4, NULL));
   InsDesc ("scale"      , TkScale,       ArgFloat(3, NULL));
   InsDesc ("call"	 , TkCall,        ArgInt(1, NULL));
   InsDesc ("calllist"	 , TkCall,        ArgInt(1, NULL));
   InsDesc ("vertex"     , TkVertex,      ArgVarFloat(2, 4, 0., 1., 0., NULL));
   InsDesc ("normal"     , TkNormal,      ArgFloat(3, NULL));
   InsDesc ("color"      , TkColor,       ArgVarFloat(3, 4, 1., 0., 0., NULL));
   InsDesc ("clear"      , TkClear,       ArgBitField(bufferTable, NULL));
   InsDesc ("clearcolor" , TkClearColor,  ArgVarFloat(3, 4, 1., 0., 0., NULL));
   InsDesc ("cleardepth" , TkClearDepth,  ArgFloat(1, NULL));
   InsDesc ("clearstencil",TkClearStencil,ArgInt(1, NULL));
   InsDesc ("enable"     , TkEnable,      ArgEnum(capabTable, NULL));
   InsDesc ("disable"    , TkDisable,     ArgEnum(capabTable, NULL));
   InsDesc ("polygonmode", TkPolygonMode, ArgEnum(faceTable, 
					  ArgEnum (pModeTable, NULL)));
   InsDesc ("material"   , TkMaterial,    ArgEnum(faceTable, 
					  ArgEnum (matParmTable, 
					  ArgVarFloat(1,4, 0., 0., 1.,NULL))));
   InsDesc ("light"      , TkLight,       ArgEnum(lightNumberTable, 
					  ArgEnum (lightParmTable, 
					  ArgVarFloat(1,4, 0., 0., 1.,NULL))));
   InsDesc ("lightmodel" , TkLightModel,  ArgEnum(lightModelTable, 
					  ArgVarFloat(1,4, .2, .2, 1.,NULL)));
   InsDesc ("shademodel" , TkShadeModel,  ArgEnum(shadeModelTable, NULL));
   InsDesc ("begin"      , TkBegin,       ArgEnum(primTable, NULL));
   InsDesc ("end"	 , TkEnd,         NULL);
   InsDesc ("perspective", TkPerspective, ArgFloat(4,NULL));
   InsDesc ("lookat"	 , TkLookAt,	  ArgFloat(9,NULL));
   InsDesc ("ortho"	 , TkOrtho,	  ArgFloat(6,NULL));
   InsDesc ("frustum"	 , TkFrustum,	  ArgFloat(6,NULL));
   InsDesc ("readpixels" , TkReadPixels,  ArgInt(2,
					  ArgString(NULL)));
   InsDesc ("drawpixels" , TkDrawPixels,  ArgString(NULL));
   InsDesc ("copypixels" , TkCopyPixels,  ArgInt (4, NULL));
   InsDesc ("pixelzoom"  , TkPixelZoom,   ArgFloat (2, NULL));
   InsDesc ("rasterpos"  , TkRasterPos,   ArgVarFloat (2, 4, 0., 1., 0.,NULL));
   InsDesc ("pixeltransfer", TkPixelTransfer, ArgEnum(pixelTransferTable, 
					      ArgFloat(1, NULL)));
   InsDesc ("texcoord"   , TkTexCoord,    ArgVarFloat (1, 4, 0., 0., 1.,NULL));
   InsDesc ("teximage2d" , TkTexImage2D,  ArgInt(2, 
					  ArgString(NULL)));
   InsDesc ("teximage1d" , TkTexImage1D,  ArgInt(2, 
					  ArgString(NULL)));
   InsDesc ("texparameter",TkTexParameter,ArgEnum(targetTexParamTable, 
					  ArgEnum(nameTexParamTable,
					  ArgDontCare(NULL))));
   InsDesc ("texenv",      TkTexEnv,      ArgEnum(targetTexEnvTable, 
					  ArgEnum(nameTexEnvTable,
					  ArgDontCare(NULL))));
   InsDesc ("texgen",      TkTexGen,      ArgEnum(coordTexGenTable, 
					  ArgEnum(nameTexGenTable,
					  ArgDontCare(NULL))));
   InsDesc ("initnames",   TkInitNames,   NULL);
   InsDesc ("loadname",    TkLoadName,    ArgInt(1,NULL));
   InsDesc ("pushname",    TkPushName,    ArgInt(1,NULL));
   InsDesc ("popname",     TkPopName,     NULL);
   InsDesc ("pickmatrix",  TkPickMatrix,  ArgFloat(4,NULL));
   InsDesc ("map1",	   TkMap1,        ArgEnum(targetMap1Table,
					  ArgFloat(2,
					  ArgInt(2,
					  ArgDontCare(NULL)))));
   InsDesc ("map2",	   TkMap2,        ArgEnum(targetMap2Table,
					  ArgFloat(2,
					  ArgInt(2,
					  ArgFloat(2,
					  ArgInt(2,
					  ArgDontCare(NULL)))))));
   InsDesc ("evalcoord1",  TkEvalCoord1,  ArgFloat (1, NULL));
   InsDesc ("evalcoord2",  TkEvalCoord2,  ArgFloat (2, NULL));
   InsDesc ("mapgrid1",    TkMapGrid1,    ArgInt (1,
					  ArgFloat (2, NULL)));
   InsDesc ("mapgrid2",    TkMapGrid2,    ArgInt (1,
					  ArgFloat (2, 
					  ArgInt (1,
					  ArgFloat (2, NULL)))));
   InsDesc ("evalmesh1",   TkEvalMesh1,   ArgEnum (modeMesh1Table,
					  ArgInt (2, NULL)));
   InsDesc ("evalmesh2",   TkEvalMesh2,   ArgEnum (modeMesh2Table,
					  ArgInt (4, NULL)));
   InsDesc ("pointsize",   TkPointSize,   ArgFloat (1, NULL));
   InsDesc ("linewidth",   TkLineWidth,   ArgFloat (1, NULL));
   InsDesc ("linestipple", TkLineStipple, ArgInt (2, NULL));
   InsDesc ("fog", 	   TkFog,	  ArgEnum (nameFogTable,
					  ArgDontCare(NULL)));
   InsDesc ("frontface",   TkFrontFace,   ArgEnum (nameFrontFaceTable, NULL));
   InsDesc ("colormaterial", TkColorMaterial, ArgEnum (faceTable,
					  ArgEnum (modeColorMatTable, NULL)));
   InsDesc ("clearaccum",  TkClearAccum,  ArgFloat (4,NULL));
   InsDesc ("drawbuffer",  TkDrawBuffer,  ArgEnum (modeDrawBufferTable, NULL));
   InsDesc ("colormask",   TkColorMask,   ArgBoolean (4, NULL));
   InsDesc ("depthmask",   TkDepthMask,   ArgBoolean (1, NULL));
   InsDesc ("stencilmask", TkStencilMask, ArgInt (1, NULL));
   InsDesc ("scissor", 	   TkScissor, 	  ArgInt (4, NULL));
   InsDesc ("alphafunc",   TkAlphaFunc,   ArgEnum (funcAlphaStencilTable,
					  ArgFloat (1, NULL)));
   InsDesc ("depthfunc",   TkDepthFunc,   ArgEnum (funcAlphaStencilTable,
					  NULL));
   InsDesc ("stencilfunc", TkStencilFunc, ArgEnum (funcAlphaStencilTable,
					  ArgInt(2, NULL)));
   InsDesc ("stencilop",   TkStencilOp,   ArgEnum (opStencilTable,
					  ArgEnum (opStencilTable,
					  ArgEnum (opStencilTable, NULL))));
   InsDesc ("accum",       TkAccum,       ArgEnum (opAccumTable,
					  ArgFloat (1, NULL)));
   InsDesc ("readbuffer",  TkReadBuffer,  ArgEnum (modeReadBufferTable, NULL));
   InsDesc ("blendfunc",   TkBlendFunc,   ArgEnum (sfactorBlendTable, 
					  ArgEnum (dfactorBlendTable, NULL)));
   InsDesc ("clipplane",   TkClipPlane,   ArgEnum (clipPlaneTable,
					  ArgFloat (4, NULL)));
   InsDesc ("pushattrib",  TkPushAttrib,  ArgBitField(attribNameTable, NULL));
   InsDesc ("popattrib",   TkPopAttrib,   NULL);
   InsDesc ("newlist",     TkNewList,     ArgInt (1, 
					  ArgEnum (newListTable, NULL)));
   InsDesc ("endlist",     TkEndList,     NULL);
   InsDesc ("edgeflag",    TkEdgeFlag,    ArgBoolean (1, NULL));
   InsDesc ("flush",       TkFlush,       NULL);
   InsDesc ("finish",      TkFinish,      NULL);
   InsDesc ("hint",        TkHint,        ArgEnum (hintTargetTable, 
                                          ArgEnum (hintModeTable, NULL)));
   InsDesc ("rect",        TkRect,        ArgFloat (4, NULL));
   InsDesc ("cullface",    TkCullFace, 	  ArgEnum (cullFaceTable, NULL));
}

     
/*---------------------------------------------------------------------------
 *
 * Main routines that parse a tcl binding for a GL function
 *
 *---------------------------------------------------------------------------*/

int 
SearchEnumVal (interp, name, val)
     Tcl_Interp* interp;
     char * name;
     GLenum* val;
{
   Tcl_HashEntry* entry;
   char buf [80];
   int i;
   for (i = 0; i < 80; i++) {
      buf [i] = tolower(name [i]);
      if (name [i] == '\0') break;
   }
   entry = Tcl_FindHashEntry (&enumHashTable, buf);
   if (entry == NULL) {
      Tcl_AppendResult (interp, "Not a valid enum:", name, (char*) NULL);
      return TCL_ERROR;
   }
   *val = (GLenum)((size_t)(Tcl_GetHashValue (entry)));
   return TCL_OK;
}

int 
SearchEnumName (interp, val, name)
     Tcl_Interp* interp;
     GLenum val;
     char ** name;
{
   Tcl_HashEntry* entry;
   Tcl_HashSearch search;
   char buf [20];
   entry = Tcl_FirstHashEntry (&enumHashTable, &search);
   while (entry != NULL) {
      if ((GLenum)((size_t)(Tcl_GetHashValue (entry))) == val) {
	 *name = Tcl_GetHashKey (&enumHashTable, entry);
	 return TCL_OK;
      }
      entry = Tcl_NextHashEntry (&search);
   }
   sprintf (buf, "%d", val);
   Tcl_AppendResult (interp, "not a known enum value: ", buf, (char*) NULL);
   return TCL_ERROR;
}

static int
SearchEnum (interp, table, tablesize, name, val) 
     Tcl_Interp* interp;
     GLenum* table;
     int tablesize;
     char * name;
     GLenum* val;
{
   int i;
   if (SearchEnumVal (interp, name, val) != TCL_OK) return TCL_ERROR;

   for (i = 0; i < tablesize; i++) {
      if (table [i] == *val) return TCL_OK;
   }

   Tcl_AppendResult (interp, "Invalid enum for this function: ", 
		     name, "\n One of the following was expected: " ,
		     (char*) NULL);

   for (i = 0; i < tablesize; i++) {
      char* enumName;
      int result = SearchEnumName (interp, table [i], &enumName);
      assert (result == TCL_OK);
      Tcl_AppendResult (interp, enumName, ", ", (char*)NULL);
   }

   return TCL_ERROR;
}


FuncDesc*
SearchFuncDesc (tclName)
     char * tclName;
{
   Tcl_HashEntry* entry;
   char buf [80];
   int i;
   for (i = 0; i < 80; i++) {
      buf [i] = tolower(tclName [i]);
      if (tclName [i] == '\0') break;
   }
   entry = Tcl_FindHashEntry (&funcDescHashTable, buf);
   if (entry == NULL) return NULL;
   return (FuncDesc*) (Tcl_GetHashValue (entry));
}

int 
ParseGLFunc (interp, argc, argv, nArg)
     Tcl_Interp *interp;
     int argc;
     char *argv [];
     int *nArg;
{
   static GLfloat floatArgs [MAXARGS];
   static GLenum enumArgs [MAXARGS];
   static GLint intArgs [MAXARGS];
   static void *argVal [MAXARGS]; 
   ArgDetailList argList;
   FuncDesc* funcDesc = NULL;
   char *glcmdstring = argv [0];
   int iarg = 0;
   int ival = 0;
   int result = TCL_OK;
   assert (argc > 0);
   if (argv [0][0] != '-') {
      ERRMSG2 ("GL command should start with '-': ", argv [0]);
   }
   funcDesc = SearchFuncDesc (&argv [0][1]);
   if (funcDesc == NULL) {
      ERRMSG2 ("Invalid GL command: ", argv [0]);
   }
   argc--;
   argv++;
   argList = funcDesc->argList;
   for (iarg = 0; 
	argList != NULL && 
	iarg < argc && 
	!(argv [iarg][0] == '-' && isalpha(argv [iarg][1])); 
	iarg++, argList = argList->next) {      
      assert (ival < MAXARGS);
      switch (argList->type) {
	 case ARG_ENUM: {
	    result = SearchEnum (interp, argList->detail.enumDetail.enumTable, 
				 argList->detail.enumDetail.enumSize,
				 argv [iarg], &enumArgs [ival]);
	    if (result != TCL_OK) goto done;
	    argVal [ival] = &enumArgs [ival];
	    ival++;
	    break;
	 }
	 case ARG_BIT_FIELD: {
	    int i;
	    int mask = 0;
	    GLenum val;
	    for (i = 0; 
		 iarg < argc &&
		 !(argv [iarg][0] == '-' && isalpha(argv [iarg][1])); 
		 iarg++, i++) {
	       result = SearchEnum (interp, 
				    argList->detail.enumDetail.enumTable, 
				    argList->detail.enumDetail.enumSize,
				    argv [iarg], &val);
	       if (result != TCL_OK) goto done;
	       mask |= val;
	    }
	    intArgs [ival] = mask;
	    argVal [ival] = &intArgs [ival];
	    ival++;
	    iarg--;
	    break;
	 }
	 case ARG_FLOAT: {
	    int i;
	    for (i = 0; 
		 iarg < argc &&
		 i < argList->detail.nValues &&
		 !(argv [iarg][0] == '-' && isalpha(argv [iarg][1])); 
		 iarg++, i++) {
	       double d;
	       assert (iarg < MAXARGS);
	       if (Tcl_GetDouble (interp, argv [iarg], &d) != TCL_OK) {
		  result = TCL_ERROR;
		  goto done;
	       }
	       floatArgs [ival] = (GLfloat)d;
	       argVal [ival] = &floatArgs [ival];
	       ival++;
	    }
	    iarg--;
	    if (i != argList->detail.nValues) {
	       ERRMSG ("Not enough arguments were specified");
	    }
	    break;
	 }
	 case ARG_DONT_CARE: {
	    int i;
	    for (i = 0; 
		 iarg < argc &&
		 !(argv [iarg][0] == '-' && isalpha(argv [iarg][1])); 
		 iarg++, i++) {
	       assert (iarg < MAXARGS);
	       argVal [ival] = argv [iarg];
	       ival++;
	    }
	    if (i != 0) iarg--;
	    break;
	 }
	 case ARG_STRING: {
	    if (iarg >= argc) {
	       ERRMSG ("Not enough arguments were specified");
	    }
	    argVal [ival] = argv [iarg];
	    ival++;
	    break;
	 }
	 case ARG_INT: {
	    int i;
	    for (i = 0; 
		 iarg < argc &&
		 i < argList->detail.nValues &&
		 !(argv [iarg][0] == '-' && isalpha(argv [iarg][1])); 
		 iarg++, i++) {
	       int j;
	       assert (iarg < MAXARGS);
	       if (Tcl_GetInt (interp, argv [iarg], &j) != TCL_OK) {
		  result = TCL_ERROR;
		  goto done;
	       }
	       intArgs [ival] = j;
	       argVal [ival] = &intArgs [ival];
	       ival++;
	    }
	    iarg--;
	    if (i != argList->detail.nValues) {
	       ERRMSG ("Not enough arguments were specified");
	    }
	    break;
	 }
	 case ARG_BOOLEAN: {
	    int i;
	    for (i = 0;
		 iarg < argc &&
		 i < argList->detail.nValues &&
		 !(argv [iarg][0] == '-' && isalpha(argv [iarg][1]));
		 iarg++, i++) {
	       int j;
	       assert (iarg < MAXARGS);
	       if (Tcl_GetBoolean (interp, argv [iarg], &j) != TCL_OK) {
		  result = TCL_ERROR;
		  goto done;
	       }
	       intArgs [ival] = j;
	       argVal [ival] = &intArgs [ival];
	       ival++;
	    }
	    iarg--;
	    if (i != argList->detail.nValues) {
	       ERRMSG ("Not enough arguments were specified");
	    }
	    break;
	 }
	 case ARG_VAR_FLOAT: {
	    int i, j;
	    for (i = 0;
		 iarg < argc &&
		 i < argList->detail.varFloatDetail.maxArgs &&
		 !(argv [iarg][0] == '-' && isalpha(argv [iarg][1])); 
		 iarg++, i++) {
	       double d;
	       assert (iarg < MAXARGS);
	       if (Tcl_GetDouble (interp, argv [iarg], &d) != TCL_OK) {
		  result = TCL_ERROR;
		  goto done;
	       }
	       floatArgs [ival] = (GLfloat)d;
	       argVal [ival] = &floatArgs [ival];
	       ival++;
	    }
	    if (i < argList->detail.varFloatDetail.minArgs) {
	       ERRMSG ("Not enough arguments were specified");
	    }
	    for (j = i-argList->detail.varFloatDetail.minArgs; 
		 i < argList->detail.varFloatDetail.maxArgs; 
		 i++,j++) {
	       floatArgs [ival] = argList->detail.varFloatDetail.def [j];
	       argVal [ival] = &floatArgs [ival];
	       ival++;
	    }
	    iarg--;		 
	    break;
	 }
      }  
   }
   if (argList != NULL) {
      ERRMSG2 ("Not enough arguments for command: ", argv [-1]);
   }
done:
   if (result == TCL_OK) {
      result = (*(funcDesc->func)) (interp, argVal, ival);
      *nArg = iarg+1;
   } 
   if (result == TCL_ERROR) {
      Tcl_AppendResult (interp, "\nError processing gl command ", 
			glcmdstring, (char*)NULL);
   }
   return result;
}

/*---------------------------------------------------------------------------
 *
 *  Functions that implement the tcl bindings of each GL function
 *
 *---------------------------------------------------------------------------*/

static int
TkMatrixMode (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs == 1);
   glMatrixMode (*((int*)args[0]));
   return TCL_OK;
}

static int 
TkPushMatrix (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs == 0);
   glPushMatrix();
   return TCL_OK;
}

static int 
TkPopMatrix (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs == 0);
   glPopMatrix();
   return TCL_OK;
}

static int 
TkLoadIdentity (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs == 0);
   glLoadIdentity();
   return TCL_OK;
}

static int 
TkMultMatrix (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs == 16);
   glMultMatrixf ((GLfloat*) args [0]);
   return TCL_OK;
}

static int 
TkLoadMatrix (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs == 16);
   glLoadMatrixf ((GLfloat*) args [0]);
   return TCL_OK;
}

static int 
TkTranslate (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{  
   GLfloat* coord;
   assert (nargs == 3);
   coord = (GLfloat*) args [0];
   glTranslatef (coord [0], coord [1], coord [2]);
   return TCL_OK;
}

static int 
TkRotate (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* coord;
   assert (nargs == 4);
   coord = (GLfloat*) args [0];
   glRotatef (coord [0], coord [1], coord [2], coord [3]);
   return TCL_OK;
}

static int 
TkScale (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* coord = (GLfloat*) args [0];
   glScalef (coord [0], coord [1], coord [2]);
   return TCL_OK;
}

static int 
TkCall (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs == 1);
   glCallList (*((int*) args [0]));
   return TCL_OK;
}

static int 
TkVertex (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* coord;
   assert (nargs >= 2 && nargs <= 4);
   coord = (GLfloat*) args [0];
   switch (nargs) {
      case 2: glVertex2fv (coord); break;
      case 3: glVertex3fv (coord); break;
      case 4: glVertex4fv (coord); break;
   }
   return TCL_OK;
}

static int 
TkNormal (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* coord;
   assert (nargs == 3);
   coord = (GLfloat*) args [0];
   glNormal3fv (coord);
   return TCL_OK;
}

static int 
TkColor (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* coord;
   assert (nargs >= 3 && nargs <= 4);
   coord = (GLfloat*) args [0];
   switch (nargs) {
      case 3: glColor3fv (coord); break;
      case 4: glColor4fv (coord); break;
   }
   return TCL_OK;
}

static int 
TkClear (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glClear (*((int*) args [0]));
   return TCL_OK;
}

static int 
TkClearColor (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* coord;
   assert (nargs == 3 || nargs == 4);
   coord = (GLfloat*) args [0];
   glClearColor (coord [0], coord [1], coord [2], coord [3]);
   return TCL_OK;
}

static int 
TkClearDepth (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glClearDepth ((double) *((GLfloat*) args [0]));
   return TCL_OK;
}

static int 
TkClearStencil (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glClearStencil (*((GLint*) args[0]));
   return TCL_OK;
}

static int 
TkEnable (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glEnable (*((GLenum*) args[0]));
   return TCL_OK;
}

static int 
TkDisable (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs == 1);
   glDisable (*((GLenum*) args[0]));
   return TCL_OK;
}

static int 
TkPolygonMode (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs == 2);
   glPolygonMode (*((GLenum*) args[0]),*((GLenum*) args[1]));
   return TCL_OK;
}

static int 
TkMaterial (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum face, parm; 
   GLfloat* val;
   assert (nargs>2);
   face = *((GLenum*) args [0]);
   parm = *((GLenum*) args [1]);
   val = ((GLfloat*) args [2]);
   glMaterialfv (face, parm, val);
   return TCL_OK;
}

static int 
TkLight (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum light, parm; 
   GLfloat* val;
   assert (nargs>2);
   light = *((GLenum*) args [0]);
   parm = *((GLenum*) args [1]);
   val = ((GLfloat*) args [2]);
   glLightfv (light, parm, val);
   return TCL_OK;
}

static int 
TkLightModel (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum parm; 
   GLfloat* val;
   assert (nargs>=2);
   parm = *((GLenum*) args [0]);
   val = ((GLfloat*) args [1]);
   glLightModelfv (parm, val);
   return TCL_OK;
}

static int 
TkShadeModel (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs==1);
   glShadeModel (*((GLenum*) args [0]));
   return TCL_OK;
}

static int 
TkPerspective (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* val = (GLfloat*) args [0];
   assert (nargs == 4);
   gluPerspective ((double) val [0], (double) val [1],
		   (double) val [2], (double) val [3]);
   return TCL_OK;
}

static int 
TkLookAt (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* val = (GLfloat*) args [0];
   assert (nargs == 9);
   gluLookAt ((double) val [0], (double) val [1], (double) val [2],
	      (double) val [3], (double) val [4], (double) val [5],
	      (double) val [6], (double) val [7], (double) val [8]);
   return TCL_OK;
}

static int 
TkOrtho (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* val = (GLfloat*) args [0];
   assert (nargs == 6);
   glOrtho ((double) val [0], (double) val [1], (double) val [2],
	    (double) val [3], (double) val [4], (double) val [5]);
   return TCL_OK;
}

static int 
TkFrustum (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* val = (GLfloat*) args [0];
   assert (nargs == 6);
   glFrustum ((double) val [0], (double) val [1], (double) val [2],
	      (double) val [3], (double) val [4], (double) val [5]);
   return TCL_OK;
}

static int 
TkBegin (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs == 1);
   glBegin (*((GLenum*) args [0]));
   return TCL_OK;
}

static int 
TkEnd (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   assert (nargs == 0);
   glEnd ();
   return TCL_OK;
}

static int 
TkReadPixels (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   int result = TCL_OK;
   char* imagename;
   int x, y;
   Tk_PhotoHandle handle;
   Tk_PhotoImageBlock block;
   assert (nargs == 3);
   x = *((int *) args [0]);
   y = *((int *) args [1]);
   imagename = (char*) args [2];
   handle = Tk_FindPhoto (interp, imagename);
   if (handle == NULL) ERRMSG2 ("Photo not defined: ", imagename);
   if (Tk_PhotoGetImage (handle, &block) != 1) 
      ERRMSG2 ("Could not get image of photo ", imagename);
   if (block.pixelSize != 3 && block.pixelSize != 4) 
      ERRMSG ("Image has invalid pixel size");
   switch (block.pitch - block.width * block.pixelSize) {
      case 0: {
	 glPixelStorei (GL_PACK_ALIGNMENT, 1);
	 break;
      }
      case 1: {
	 glPixelStorei (GL_PACK_ALIGNMENT, 2);
	 break;
      }
      case 2:
      case 3: {
	 glPixelStorei (GL_PACK_ALIGNMENT, 4);
	 break;
      }	 
      default: 
         printf ("unknown alignment\n");
   }
   glReadPixels (x, y, block.width, block.height,
		 block.pixelSize == 3 ? GL_RGB : GL_RGBA,
		 GL_UNSIGNED_BYTE, block.pixelPtr);
   /* Swap rows so that image will not end upside down */
   {
      int iy;
      char *tmp;
      tmp = (char*) malloc (block.pitch);
      for (iy = 0; iy < block.height/2; iy++) {
	 char *from = (char*)block.pixelPtr+iy*block.pitch;
	 char *to = (char*)block.pixelPtr+(block.height-iy-1)*block.pitch;
	 memcpy (tmp, from, block.pitch);
	 memcpy (from, to, block.pitch);
	 memcpy (to, tmp, block.pitch);
      }
      free (tmp);
   }
done:
   return result;
}   



static int 
TkDrawPixels (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   int result = TCL_OK;
   char* imagename;
   Tk_PhotoHandle handle;
   Tk_PhotoImageBlock block;
   assert (nargs == 1);
   imagename = (char*) args [0];
   handle = Tk_FindPhoto (interp, imagename);
   if (handle == NULL) ERRMSG2 ("Photo not defined: ", imagename);
   if (Tk_PhotoGetImage (handle, &block) != 1) 
      ERRMSG2 ("Could not get image of photo ", imagename);
   if (block.pixelSize != 3 && block.pixelSize != 4) 
      ERRMSG ("Image has invalid pixel size");
   glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
   glDrawPixels (block.width, block.height,
		 block.pixelSize == 3 ? GL_RGB : GL_RGBA,
		 GL_UNSIGNED_BYTE, block.pixelPtr);
done:
   return result;
}   


static int 
TkCopyPixels (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLint *val = (GLint*) args [0];
   assert (nargs == 4);
   glCopyPixels (val [0], val [1], val [2], val [3], GL_COLOR);
   return TCL_OK;
}

static int 
TkRasterPos (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat *val = (GLfloat*) args [0];
   assert (nargs >= 2);
   glRasterPos4fv (val);
   return TCL_OK;
}

static int 
TkPixelZoom (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat *val = (GLfloat*) args [0];
   assert (nargs == 2);
   glPixelZoom (val [0], val [1]);
   return TCL_OK;
}

static int 
TkPixelTransfer (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum parm;
   GLfloat val;
   assert (nargs == 2);
   parm = *((GLenum*) args [0]);
   val = *((GLfloat*) args [1]);
   glPixelTransferf (parm, val);
   return TCL_OK;
}

static int 
TkTexCoord (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* coord;
   assert (nargs >= 2 && nargs <= 4);
   coord = (GLfloat*) args [0];
   switch (nargs) {
      case 2: glTexCoord2fv (coord); break;
      case 3: glTexCoord3fv (coord); break;
      case 4: glTexCoord4fv (coord); break;
   }
   return TCL_OK;
}

#ifndef NDEBUG

/*	Create checkerboard texture	*/
#define	checkImageWidth 64
#define	checkImageHeight 64
GLubyte checkImage[checkImageWidth][checkImageHeight][3];

void makeCheckImage(void)
{
    int i, j, c;
    
    for (i = 0; i < checkImageWidth; i++) {
	for (j = 0; j < checkImageHeight; j++) {
	    c = (((i&0x8)==0)^((j&0x8)==0))*255;
	    checkImage[i][j][0] = (GLubyte) c;
	    checkImage[i][j][1] = (GLubyte) c;
	    checkImage[i][j][2] = (GLubyte) c;
	}
    }
}
#endif

static int 
TkTexImage2D (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   int result = TCL_OK;
   int level, border, i, n;
   char* imagename;
   Tk_PhotoHandle handle;
   Tk_PhotoImageBlock block;
   assert (nargs == 3);
   level = *((int*) args [0]);
   border = *((int*) args [1]);
   imagename = (char*) args [2];
   handle = Tk_FindPhoto (interp, imagename);
   if (handle == NULL) ERRMSG2 ("Photo not defined: ", imagename);
   if (Tk_PhotoGetImage (handle, &block) != 1) 
      ERRMSG2 ("Could not get image of photo ", imagename);
   if (block.pixelSize != 3 && block.pixelSize != 4) 
      ERRMSG ("Image has invalid pixel size");
   n = block.width - border;
   for (i = 0; i < 16; i++) {
      if (n == (1 << i)) break;
   }
   if (i == 16) {
      char buf [10];
      sprintf (buf, "%d", block.width);
      ERRMSG2 ("image width must be a power of 2", buf);
   }
   n = block.height - border;
   for (i = 0; i < 16; i++) {
      if (n == (1 << i)) break;
   }
   if (i == 16) {
      char buf [20];
      sprintf (buf, "%d", block.height);
      ERRMSG2 ("image height must be a power of 2", buf);
   }
   glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
   glTexImage2D (GL_TEXTURE_2D, level, 
		 block.pixelSize, 
		 block.width, block.height, border,
		 block.pixelSize == 3 ? GL_RGB : GL_RGBA,
		 GL_UNSIGNED_BYTE, block.pixelPtr);
done:
   return result;
}

static int 
TkTexImage1D (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   int result = TCL_OK;
   int level, border, i, n;
   char* imagename;
   Tk_PhotoHandle handle;
   Tk_PhotoImageBlock block;
   assert (nargs == 3);
   level = *((int*) args [0]);
   border = *((int*) args [1]);
   imagename = (char*) args [2];
   handle = Tk_FindPhoto (interp, imagename);
   if (handle == NULL) ERRMSG2 ("Photo not defined: ", imagename);
   if (Tk_PhotoGetImage (handle, &block) != 1) 
      ERRMSG2 ("Could not get image of photo ", imagename);
   if (block.pixelSize != 3 && block.pixelSize != 4) 
      ERRMSG ("Image has invalid pixel size");
   n = block.width - border;
   for (i = 0; i < 16; i++) {
      if (n == (1 << i)) break;
   }
   if (i == 16) {
      char buf [20];
      sprintf (buf, "%d", block.width);
      ERRMSG2 ("image width must be a power of 2", buf);
   }
   glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
   glTexImage1D (GL_TEXTURE_1D, level, 
		 block.pixelSize, 
		 block.width, border,
		 block.pixelSize == 3 ? GL_RGB : GL_RGBA,
		 GL_UNSIGNED_BYTE, block.pixelPtr);
done:
   return result;
}

static int 
TkTexParameter (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   int result = TCL_OK;
   GLenum target, pname, pval;
   GLfloat val [4] = {0.0, 0.0, 0.0, 0.0};
   assert (nargs >= 3 && nargs <= 6);
   target = *((GLenum*) args [0]);
   pname = *((GLenum*) args [1]);
   if (pname == GL_TEXTURE_BORDER_COLOR) {
      int i;
      for (i = 0; i+2 < nargs; i++) {
	 double d;
	 if (Tcl_GetDouble (interp, (char*) args [i+2], &d) != TCL_OK) {
	    ERRMSG ("Invalid color coordinate");
	 }
	 val [i] = (GLfloat)d;
      }
      glTexParameterfv (target, pname, val);
   }
   else {
      result = SearchEnum (interp, valueTexParamTable, 
			   sizeof (valueTexParamTable)/sizeof(GLenum),
			   (char*) args [2], &pval);
      if (result != TCL_OK) goto done;
      glTexParameteri (target, pname, pval);
   }
done:
   return result;
}

static int 
TkTexEnv (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   int result = TCL_OK;
   GLenum target, pname, pval;
   GLfloat val [4] = {0.0, 0.0, 0.0, 0.0};
   assert (nargs >= 3 && nargs <= 6);
   target = *((GLenum*) args [0]);
   pname = *((GLenum*) args [1]);
   if (pname == GL_TEXTURE_ENV_COLOR) {
      int i;
      for (i = 0; i+2 < nargs; i++) {
	 double d;
	 if (Tcl_GetDouble (interp, (char*) args [i+2], &d) != TCL_OK) {
	    ERRMSG ("Invalid color coordinate");
	 }
	 val [i] = (GLfloat)d;
      }
      glTexEnvfv (target, pname, val);
   }
   else {
      result = SearchEnum (interp, valueTexEnvTable, 
			   sizeof(valueTexEnvTable)/sizeof(GLenum),
			   (char*) args [2], &pval);
      if (result != TCL_OK) goto done;
      glTexEnvi (target, pname, pval);
   }
done:
   return result;
}

static int 
TkTexGen (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   int result = TCL_OK;
   GLenum coord, pname, pval;
   GLfloat val [4] = {0.0, 0.0, 0.0, 0.0};
   assert (nargs >= 3 && nargs <= 6);
   coord = *((GLenum*) args [0]);
   pname = *((GLenum*) args [1]);
   if (pname != GL_TEXTURE_GEN_MODE) {
      int i;
      for (i = 0; i+2 < nargs; i++) {
	 double d;
	 if (Tcl_GetDouble (interp, (char*) args [i+2], &d) != TCL_OK) {
	    ERRMSG ("Invalid texgen coordinate");
	 }
	 val [i] = (GLfloat)d;
      }
      glTexGenfv (coord, pname, val);
   }
   else {
      result = SearchEnum (interp, valueTexGenTable,
			   sizeof(valueTexGenTable)/sizeof(GLenum),
			   (char*) args [2], &pval);
      if (result != TCL_OK) goto done;
      glTexGeni (coord, pname, pval);
   }
done:
   return result;
}

static int 
TkInitNames (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glInitNames();
   return TCL_OK;
}

static int 
TkLoadName (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glLoadName(*((GLint*) args [0]));
   return TCL_OK;
}


static int 
TkPushName (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glPushName(*((GLint*) args [0]));
   return TCL_OK;
}


static int 
TkPopName (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glPopName ();
   return TCL_OK;
}

static int 
TkPickMatrix (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat *val = ((GLfloat*) args [0]);
   GLint viewport [4];

   glGetIntegerv (GL_VIEWPORT, viewport);
   gluPickMatrix ((double) val [0], (double) (viewport [3] - val [1]),
		  val [2], val [3], viewport);
   
   return TCL_OK;
}

static int
TkMap1 (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum target = *((GLenum*) args [0]);
   GLfloat *u = ((GLfloat*) args [1]);
   GLint stride = *((GLint*) args [3]);
   GLint order = *((GLint*) args [4]);
   int result = TCL_OK;
   int i;
   GLfloat *pt;

   if (nargs - 5 != order * stride) {
      char buf [20];
      sprintf (buf, "%d", order * stride);
      ERRMSG2 (buf, " coordinate values were expected");
   }

   pt = (GLfloat*) malloc (sizeof (GLfloat) * stride * order);
   assert (pt != NULL);

   for (i = 0; i+5 < nargs; i++) {
      double d;
      if (Tcl_GetDouble (interp, (char*) args [i+5], &d) != TCL_OK) {
	 free (pt);
	 ERRMSG ("Invalid map1 coordinate");
      }
      pt [i] = (GLfloat)d;
   }
   glMap1f (target, u [0], u [1], stride, order, pt);
   free (pt);

done:
   return result;
}

static int
TkMap2 (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum target = *((GLenum*) args [0]);
   GLfloat *u = ((GLfloat*) args [1]);
   GLint ustride = *((GLint*) args [3]);
   GLint uorder = *((GLint*) args [4]);
   GLfloat *v = ((GLfloat*) args [5]);
   GLint vstride = *((GLint*) args [7]);
   GLint vorder = *((GLint*) args [8]);
   int result = TCL_OK;
   int i;
   GLfloat *pt;

   if (nargs - 9 != uorder * vorder * ustride) {
      char buf [20];
      sprintf (buf, "%d", uorder * vorder * ustride);
      ERRMSG2 (buf, " coordinate values were expected");
   }

   if (vstride != ustride * uorder) {
      char buf [20];
      sprintf (buf, "%d", uorder * ustride);
      ERRMSG2 (" vstride parameter should be ", buf);
   }


   pt = (GLfloat*) malloc (sizeof (GLfloat) * ustride * vorder * uorder);
   assert (pt != NULL);

   for (i = 0; i+9 < nargs; i++) {
      double d;
      if (Tcl_GetDouble (interp, (char*) args [i+9], &d) != TCL_OK) {
	 free (pt);
	 ERRMSG ("Invalid map1 coordinate");
      }
      pt [i] = (GLfloat)d;
   }
   glMap2f (target, u [0], u [1], ustride, uorder, 
	    v [0], v [1], vstride, vorder, pt);
   free (pt);

done:
   return result;
}


static int
TkEvalCoord1 (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat *u = ((GLfloat*) args [0]);
   glEvalCoord1f (u[0]);
   return TCL_OK;
}

static int
TkEvalCoord2 (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat *u = ((GLfloat*) args [0]);
   glEvalCoord2f (u[0], u[1]);
   return TCL_OK;
}

static int
TkMapGrid1 (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLint n = *((GLint*) args [0]);
   GLfloat *u = ((GLfloat*) args [1]);
   glMapGrid1f (n, u[0], u[1]);
   return TCL_OK;
}

static int
TkMapGrid2 (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLint nu = *((GLint*) args [0]);
   GLfloat *u = ((GLfloat*) args [1]);
   GLint nv = *((GLint*) args [3]);
   GLfloat *v = ((GLfloat*) args [4]);
   glMapGrid2f (nu, u[0], u[1], nv, v[0], v[1]);
   return TCL_OK;
}

static int
TkEvalMesh1 (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum mode = *((GLenum*) args [0]);
   GLint *p = ((GLint*) args [1]);
   glEvalMesh1 (mode, p[0], p[1]);
   return TCL_OK;
}

static int
TkEvalMesh2 (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum mode = *((GLenum*) args [0]);
   GLint *p = ((GLint*) args [1]);
   glEvalMesh2 (mode, p[0], p[1], p[2], p[3]);
   return TCL_OK;
}


static int
TkPointSize (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat size = *((GLfloat*) args [0]);
   glPointSize (size);
   return TCL_OK;
}

static int
TkLineWidth (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat width = *((GLfloat*) args [0]);
   glLineWidth (width);
   return TCL_OK;
}

static int
TkLineStipple (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLint* arg = ((GLint*) args [0]);
   glLineStipple (arg [0], (GLushort) arg [1]);
   return TCL_OK;
}

static int
TkFog (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   int result = TCL_OK;
   GLenum mode = *((GLenum*) args [0]);
   if (mode == GL_FOG_MODE) {
      GLenum fogmode;
      result = SearchEnum (interp, fogFogModeTable,
			   sizeof(fogFogModeTable)/sizeof(GLenum),
			   (char*) args [1], &fogmode);
      if (result != TCL_OK) goto done;
      glFogi (mode, fogmode);
   }
   else if (mode == GL_FOG_COLOR) {
      int i;
      GLfloat fogcolor[4] = {0.0, 0.0, 0.0, 1.0};
      for (i = 0; i+1 < 4 && i+1 < nargs; i++) {
	 double d;
	 if (Tcl_GetDouble (interp, (char*)args[i+1], &d) !=TCL_OK) {
	    ERRMSG ("\nInvalid fog color coordinate");
	 }
	 fogcolor [i] = (GLfloat)d;
      }
      glFogfv (mode, fogcolor);
   }
   else {
      double d;
      if (Tcl_GetDouble (interp, (char*) args [1], &d) != TCL_OK) {
	 ERRMSG ("\nInvalid fog argument");
      }
      glFogf (mode, (GLfloat)d);
   }
done:
   return result;
}

static int
TkFrontFace (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum mode = *((GLenum*) args [0]);
   glFrontFace (mode);
   return TCL_OK;
}

static int
TkColorMaterial (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum face = *((GLenum*) args [0]);
   GLenum mode = *((GLenum*) args [1]);
   glColorMaterial (face, mode);
   return TCL_OK;
}

static int
TkClearAccum (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* s = ((GLfloat*) args [0]);
   glClearAccum (s[0], s[1], s[2], s[3]);
   return TCL_OK;
}

static int
TkDrawBuffer (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum mode = *((GLenum*) args [0]);
   glDrawBuffer (mode);
   return TCL_OK;
}

static int
TkColorMask (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLboolean* s = ((GLboolean*) args [0]);
   glColorMask (s[0], s[1], s[2], s[3]);
   return TCL_OK;
}

static int
TkDepthMask (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLboolean* s = ((GLboolean*) args [0]);
   glDepthMask (s[0]);
   return TCL_OK;
}

static int
TkStencilMask (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLint* s = ((GLint*) args [0]);
   glStencilMask (s[0]);
   return TCL_OK;
}

static int
TkScissor (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLint* s = ((GLint*) args [0]);
   glScissor (s[0],s[1],s[2],s[3]);
   return TCL_OK;
}

static int
TkAlphaFunc (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum func = *((GLenum*) args [0]);
   GLfloat ref = *((GLfloat*) args [1]);
   glAlphaFunc (func, ref);
   return TCL_OK;
}

static int
TkDepthFunc (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum func = *((GLenum*) args [0]);
   glDepthFunc (func);
   return TCL_OK;
}

static int
TkStencilFunc (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum func = *((GLenum*) args [0]);
   GLint ref = *((GLint*) args [1]);
   GLint mask = *((GLint*) args [2]);
   glStencilFunc (func, ref, mask);
   return TCL_OK;
}

static int
TkStencilOp (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum* op = ((GLenum*) args [0]);
   glStencilOp (op [0], op [1], op [2]);
   return TCL_OK;
}

static int
TkAccum (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum op = *((GLenum*) args [0]);
   GLfloat val = *((GLfloat*) args [1]);
   glAccum (op, val);
   return TCL_OK;
}

static int
TkReadBuffer(interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum mode = *((GLenum*) args [0]);
   glReadBuffer (mode);
   return TCL_OK;
}

static int
TkBlendFunc (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum sfactor = *((GLenum*) args [0]);
   GLenum dfactor = *((GLenum*) args [1]);
   glBlendFunc (sfactor, dfactor);
   return TCL_OK;
}

static int
TkClipPlane (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLenum clipplane = *((GLenum*) args [0]);
   GLdouble eq [4];
   eq [0]= *((GLfloat*) args [1]);
   eq [1]= *((GLfloat*) args [2]);
   eq [2]= *((GLfloat*) args [3]);
   eq [3]= *((GLfloat*) args [4]);
   glClipPlane (clipplane, eq);
   return TCL_OK;
}

static int
TkPushAttrib (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glPushAttrib(*((GLint*) args [0]));
   return TCL_OK;
}

static int
TkPopAttrib (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glPopAttrib();
   return TCL_OK;
}


static int
TkNewList (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glNewList(*((GLint*) args [0]), *((GLenum*) args [1]));
   return TCL_OK;
}

static int
TkEndList (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glEndList ();
   return TCL_OK;
}

static int
TkEdgeFlag (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLboolean* s = ((GLboolean*) args [0]);
   glEdgeFlag (s[0]);
   return TCL_OK;
}

static int
TkFlush (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glFlush ();
   return TCL_OK;
}

static int
TkFinish (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glFinish ();
   return TCL_OK;
}

static int
TkHint (interp, args, nargs)
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
      glHint (*((GLenum*) args [0]), *((GLenum*) args [1]));
   return TCL_OK;
}

static int 
TkRect (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   GLfloat* coord;
   coord = (GLfloat*) args [0];
   glRectf (coord [0], coord [1], coord [2], coord [3]);
   return TCL_OK;
}

static int 
TkCullFace (interp, args, nargs) 
     Tcl_Interp* interp;
     void** args;
     int nargs;
{
   glCullFace (*((GLenum*) args [0]));
   return TCL_OK;
}


