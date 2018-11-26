#include "stdafx.h"
#include "MASTER.h"
 #define ___1404
#include "GLOBAL.h"
#include "TASSERT.h"
#include "ALLOC.h"
#include "SYSTEM.h"
#include "FILESTREAM.h"
___1405 *___1402(FILE      *File, ___372  ___2007) { REQUIRE(VALID_REF(File) || File == NULL); ___1405 *___3359 = ALLOC_ITEM(___1405, "FileStream"); if (___3359 != NULL) { ___3359->File              = File; ___3359->___2007 = ___2007; } ENSURE(VALID_REF(___3359) || ___3359 == NULL); return ___3359; } void ___1403(___1405 **___1401) { REQUIRE(VALID_REF(___1401)); REQUIRE(VALID_REF(*___1401) || *___1401 == NULL); if (*___1401 != NULL) { ___1531(*___1401, "FileStream"); *___1401 = NULL; } ENSURE(*___1401 == NULL); }
