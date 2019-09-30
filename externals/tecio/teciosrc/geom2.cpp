#include "stdafx.h"
#include "MASTER.h"
 #define ___1559
#include "GLOBAL.h"
#include "TASSERT.h"
#include "ALLOC.h"
#include "GEOM.h"
#include "TEXT.h"
#include "STRUTIL.h"
#include "GEOM2.h"
#include "DATASET0.h"
FieldDataType_e ___1747(___1632 const* ___1556) { FieldDataType_e ___3359; REQUIRE(VALID_REF(___1556)); REQUIRE(VALID_REF(___1556->___1573.___1548.___4293)); ___3359 = ___1556->___907; ENSURE(VALID_GEOM_FIELD_DATA_TYPE(___3359)); ENSURE(IMPLICATION(VALID_REF(___1556->___1573.___1548.___4293), ___3359 == ___1726(___1556->___1573.___1548.___4293))); ENSURE(IMPLICATION(VALID_REF(___1556->___1573.___1548.___4295), ___3359 == ___1726(___1556->___1573.___1548.___4295))); ENSURE(IMPLICATION(VALID_REF(___1556->___1573.___1548.___4297), ___3359 == ___1726(___1556->___1573.___1548.___4297))); return ___3359; }
