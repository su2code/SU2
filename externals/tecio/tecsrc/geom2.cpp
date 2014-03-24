#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE

/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2010 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#define GEOM2MODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "ALLOC.h"

#include "GEOM.h"
#include "TEXT.h"
#include "STRUTIL.h"
#include "GEOM2.h"

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#include "DATASET0.h"


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# if 0 /* 3D geometry arrowheads are not drawn at this time. */
#endif
# if 0 /* 3D geometry arrowheads are not drawn at this time. */
# endif
#if 0
#endif
#             ifndef NO_ASSERTS
#             endif
#             ifndef NO_ASSERTS
#             endif
#endif /* TECPLOTKERNEL */


FieldDataType_e GetGeomFieldDataType(Geom_s const* Geom)
{
    FieldDataType_e Result;

    REQUIRE(VALID_REF(Geom));
    REQUIRE(VALID_REF(Geom->GeomData.Generic.V1Base));

    Result = Geom->DataType;

    ENSURE(VALID_GEOM_FIELD_DATA_TYPE(Result));
    /*
     * Check that the geom's field data arrays (if they exist)
     * have the same type as the geometry.
     */
    ENSURE(IMPLICATION(VALID_REF(Geom->GeomData.Generic.V1Base), Result == GetFieldDataType(Geom->GeomData.Generic.V1Base)));
    ENSURE(IMPLICATION(VALID_REF(Geom->GeomData.Generic.V2Base), Result == GetFieldDataType(Geom->GeomData.Generic.V2Base)));
    ENSURE(IMPLICATION(VALID_REF(Geom->GeomData.Generic.V3Base), Result == GetFieldDataType(Geom->GeomData.Generic.V3Base)));

    return Result;
}
