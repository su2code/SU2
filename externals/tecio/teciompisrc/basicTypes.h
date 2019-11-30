 #pragma once
#include "MASTER.h"
#include "GLOBAL.h"
#include "TASSERT.h"
#include "ItemAddress.h"
namespace tecplot { namespace ___3933 { typedef int32_t ___4636; ___4636 const MAX_NUM_ZONES = ___2391; ___4636 const BAD_ZONE_INDEX = ___4636(0x7FFFFFFF); typedef int32_t ___4352; ___4352 const MAX_NUM_VARS = ___2391; ___4352 const BAD_VAR_INDEX = ___4352(0x7FFFFFFF); typedef int64_t ___81; ___81 const MAX_ANY_INDEX = ___81(0x7FFFFFFFFFFFFFFEll); ___81 const BAD_ANY_INDEX = ___81(0x7FFFFFFFFFFFFFFFll); typedef int64_t ___465; ___465 const MAX_NUM_CELLS = ___465(0x7FFFFFFFFFFFFFFEll); ___465 const BAD_CELL_INDEX = ___465(0x7FFFFFFFFFFFFFFFll); typedef int64_t ___2718; ___2718 const MAX_NUM_NODES = ___2718(0x7FFFFFFFFFFFFFFEll); ___2718 const BAD_NODE_INDEX = ___2718(0x7FFFFFFFFFFFFFFFll); typedef uint8_t ___682; ___682 const MAX_NUM_CELL_CORNERS = 8; ___682 const BAD_CORNER_INDEX = ___682(-1); ___682 const NUM_IJK_CELL_CORNERS = 8; ___682 const NUM_TETRA_CELL_CORNERS = 4; ___682 const NUM_BRICK_CELL_CORNERS = 8; typedef uint8_t FaceIndex_t; FaceIndex_t const MAX_NUM_CELL_FACES = 6; FaceIndex_t const BAD_FACE_INDEX = FaceIndex_t(-1); FaceIndex_t const NUM_IJK_CELL_FACES = 6; FaceIndex_t const NUM_TETRA_CELL_FACES = 4; FaceIndex_t const NUM_BRICK_CELL_FACES = 6; typedef uint64_t ___1393; ___1393 const ___330 = ___1393(-1);
 #define VALID_FILE_LOC(fileLoc) ( (fileLoc) != ___330 )
typedef uint16_t RefSubzoneOffset_t; RefSubzoneOffset_t const BAD_REFSZ_INDEX = RefSubzoneOffset_t(-1); }}
