/****************************************************************************
 *
 * CAEP Plugin example
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#define SHOW_PWP_MESSAGES
#define SHOW_PWGM_MESSAGES
#include "apiGridModel.h"

#define MODEL      PWGM_HGRIDMODEL_GMIMPL(model)
#define BLKMODEL   PWGM_HBLOCK_GMIMPL(block)
#define DOMMODEL   PWGM_HDOMAIN_GMIMPL(domain)
#define VERTMODEL  PWGM_HVERTEX_GMIMPL(vertex)
#define ELEMMODEL  PWGM_HELEMENT_GMIMPL(element)
#define BNDRYMODEL PWGM_HBNDRY_GMIMPL(boundary)
#define CNXNMODEL  PWGM_HCNXN_GMIMPL(connection)

const PWGM_HBLOCK badBlock = PWGM_HBLOCK_INIT;
const PWGM_HDOMAIN badDomain = PWGM_HDOMAIN_INIT;
const PWGM_HVERTEX badVertex = PWGM_HVERTEX_INIT;
const PWGM_HELEMENT badElement = PWGM_HELEMENT_INIT;
const PWGM_HBNDRY badBoundary = PWGM_HBNDRY_INIT;
const PWGM_HCNXN badConnection = PWGM_HCNXN_INIT;


/**************************************/
PWP_UINT32
PwModBlockCount(PWGM_HGRIDMODEL model)
{
    return (MODEL && MODEL->PwModBlockCountCB) ?
        MODEL->PwModBlockCountCB(model) : 0;
}

/**************************************/
PWGM_HBLOCK
PwModEnumBlocks(PWGM_HGRIDMODEL model, PWP_UINT32 ndx)
{
    return (MODEL && MODEL->PwModEnumBlocksCB) ?
        MODEL->PwModEnumBlocksCB(model, ndx) : badBlock;
}



#if !defined(PWGM_HIDE_UNSTRUCTURED_API)

/**************************************/
PWP_UINT32
PwModDomainCount(PWGM_HGRIDMODEL model)
{
    return (MODEL && MODEL->PwModDomainCountCB) ?
        MODEL->PwModDomainCountCB(model) : 0;
}

/**************************************/
PWGM_HDOMAIN
PwModEnumDomains(PWGM_HGRIDMODEL model, PWP_UINT32 ndx)
{
    return (MODEL && MODEL->PwModEnumDomainsCB) ?
        MODEL->PwModEnumDomainsCB(model, ndx) : badDomain;
}

/**************************************/
PWGM_HVERTEX
PwModEnumVertices(PWGM_HGRIDMODEL model, PWP_UINT32 ndx)
{
    return (MODEL && MODEL->PwModEnumVerticesCB) ?
        MODEL->PwModEnumVerticesCB(model, ndx) : badVertex;
}

/**************************************/
PWP_UINT32
PwModVertexCount(PWGM_HGRIDMODEL model)
{
    return (MODEL && MODEL->PwModVertexCountCB) ?
        MODEL->PwModVertexCountCB(model) : 0;
}

/**************************************/
PWP_UINT32
PwBlkElementCount(PWGM_HBLOCK block, PWGM_ELEMCOUNTS *pCounts)
{
    return (PWGM_HBLOCK_ISVALID(block) && BLKMODEL->PwBlkElementCountCB) ?
        BLKMODEL->PwBlkElementCountCB(block, pCounts) : 0;
}

/**************************************/
PWGM_HELEMENT
PwBlkEnumElements(PWGM_HBLOCK block, PWP_UINT32 ndx)
{
    return (PWGM_HBLOCK_ISVALID(block) && BLKMODEL->PwBlkEnumElementsCB) ?
        BLKMODEL->PwBlkEnumElementsCB(block, ndx) : badElement;
}

/**************************************/
PWP_BOOL
PwBlkCondition(PWGM_HBLOCK block, PWGM_CONDDATA *pCondData)
{
    return (pCondData && PWGM_HBLOCK_ISVALID(block) &&
            BLKMODEL->PwBlkConditionCB) ?
        BLKMODEL->PwBlkConditionCB(block, pCondData) : PWP_FALSE;
}

/**************************************/
PWP_UINT32
PwDomElementCount(PWGM_HDOMAIN domain, PWGM_ELEMCOUNTS *pCounts)
{
    return (PWGM_HDOMAIN_ISVALID(domain) && DOMMODEL->PwDomElementCountCB) ?
        DOMMODEL->PwDomElementCountCB(domain, pCounts) : 0;
}

/**************************************/
PWGM_HELEMENT
PwDomEnumElements(PWGM_HDOMAIN domain, PWP_UINT32 ndx)
{
    return (PWGM_HDOMAIN_ISVALID(domain) && DOMMODEL->PwDomEnumElementsCB) ?
        DOMMODEL->PwDomEnumElementsCB(domain, ndx) : badElement;
}

/**************************************/
PWP_BOOL
PwDomCondition(PWGM_HDOMAIN domain, PWGM_CONDDATA *pCondData)
{
    return (pCondData && PWGM_HDOMAIN_ISVALID(domain) &&
            DOMMODEL->PwDomConditionCB) ?
        DOMMODEL->PwDomConditionCB(domain, pCondData) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwVertDataMod(PWGM_HVERTEX vertex, PWGM_VERTDATA *pVertData)
{
    return (pVertData && PWGM_HVERTEX_ISVALID(vertex) &&
            VERTMODEL->PwVertDataModCB) ?
        VERTMODEL->PwVertDataModCB(vertex, pVertData) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwVertIndexMod(PWGM_HVERTEX vertex, PWP_UINT32 *pIndex)
{
    return (pIndex && PWGM_HVERTEX_ISVALID(vertex) &&
            VERTMODEL->PwVertIndexModCB) ?
        VERTMODEL->PwVertIndexModCB(vertex, pIndex) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwVertXyzVal(PWGM_HVERTEX vertex, PWGM_ENUM_XYZ which, PWGM_XYZVAL *pVal)
{
    return (pVal && PWGM_HVERTEX_ISVALID(vertex) &&
            VERTMODEL->PwVertXyzValCB) ?
        VERTMODEL->PwVertXyzValCB(vertex, which, pVal) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwElemDataMod(PWGM_HELEMENT element, PWGM_ELEMDATA *pElemData)
{
    return (pElemData && PWGM_HELEMENT_ISVALID(element) &&
            ELEMMODEL->PwElemDataModCB) ?
        ELEMMODEL->PwElemDataModCB(element, pElemData) : PWP_FALSE;
}

#endif // !defined(PWGM_HIDE_UNSTRUCTURED_API)



//***************************************************************************
//***************************************************************************
//***************************************************************************
//***************************************************************************

#if !defined(PWGM_HIDE_STRUCTURED_API)

/**************************************/
PWP_UINT32
PwModBoundaryCount(PWGM_HGRIDMODEL model)
{
    return (MODEL && MODEL->PwModBoundaryCountCB) ?
        MODEL->PwModBoundaryCountCB(model) : 0;
}

/**************************************/
PWGM_HBNDRY
PwModEnumBoundaries(PWGM_HGRIDMODEL model, PWP_UINT32 ndx)
{
    return (MODEL && MODEL->PwModEnumBoundariesCB) ?
        MODEL->PwModEnumBoundariesCB(model, ndx) : badBoundary;
}

/**************************************/
PWP_UINT32
PwModConnectionCount(PWGM_HGRIDMODEL model)
{
    return (MODEL && MODEL->PwModConnectionCountCB) ?
        MODEL->PwModConnectionCountCB(model) : 0;
}

/**************************************/
PWGM_HCNXN
PwModEnumConnections(PWGM_HGRIDMODEL model, PWP_UINT32 ndx)
{
    return (MODEL && MODEL->PwModEnumConnectionsCB) ?
        MODEL->PwModEnumConnectionsCB(model, ndx) : badConnection;
}

/**************************************/
PWP_BOOL
PwModNdxBoundary(PWGM_HGRIDMODEL model, PWP_UINT32 ndx,
    PWGM_BNDRYDATA *pBndryData)
{
    return (MODEL && MODEL->PwModNdxBoundaryCB) ?
            MODEL->PwModNdxBoundaryCB(model, ndx, pBndryData) :
            PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwModNdxBoundaryAndCondition(PWGM_HGRIDMODEL model, PWP_UINT32 ndx,
    PWGM_BNDRYDATA *pBndryData, PWGM_CONDDATA *pCondData)
{
    return (MODEL && MODEL->PwModNdxBoundaryAndConditionCB &&
            (pBndryData || pCondData)) ?
        MODEL->PwModNdxBoundaryAndConditionCB(model, ndx, pBndryData,
            pCondData) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwModNdxConnection(PWGM_HGRIDMODEL model, PWP_UINT32 ndx,
    PWGM_CNXNDATA *pCnxnData)
{
    return (MODEL && MODEL->PwModNdxConnectionCB && pCnxnData) ?
        MODEL->PwModNdxConnectionCB(model, ndx, pCnxnData) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwBlkSize(PWGM_HBLOCK block, PWGM_STR_SIZE *pSize)
{
    return (BLKMODEL && BLKMODEL->PwBlkSizeCB) ?
        BLKMODEL->PwBlkSizeCB(block, pSize) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwBlkNdxVertData(PWGM_HBLOCK block, PWGM_INDEX3 ndx3, PWGM_VERTDATA *pVertData)
{
    return (BLKMODEL && BLKMODEL->PwBlkNdxVertDataCB) ?
        BLKMODEL->PwBlkNdxVertDataCB(block, ndx3, pVertData) : PWP_FALSE;
}

/**************************************/
PWP_UINT32
PwBlkBoundaryCount(PWGM_HBLOCK block)
{
    return (BLKMODEL && BLKMODEL->PwBlkBoundaryCountCB) ?
        BLKMODEL->PwBlkBoundaryCountCB(block) : 0;
}

/**************************************/
PWGM_HBNDRY
PwBlkEnumBoundaries(PWGM_HBLOCK block, PWP_UINT32 ndx)
{
    return (BLKMODEL && BLKMODEL->PwBlkEnumBoundariesCB) ?
        BLKMODEL->PwBlkEnumBoundariesCB(block, ndx) : badBoundary;
}

/**************************************/
PWP_UINT32
PwBlkConnectionCount(PWGM_HBLOCK block)
{
    return (BLKMODEL && BLKMODEL->PwBlkConnectionCountCB) ?
        BLKMODEL->PwBlkConnectionCountCB(block) : 0;
}

/**************************************/
PWGM_HCNXN
PwBlkEnumConnections(PWGM_HBLOCK block, PWP_UINT32 ndx)
{
    return (BLKMODEL && BLKMODEL->PwBlkEnumConnectionsCB) ?
        BLKMODEL->PwBlkEnumConnectionsCB(block, ndx) : badConnection;
}

/**************************************/
PWP_BOOL
PwBlkNdxBoundary(PWGM_HBLOCK block, PWP_UINT32 ndx, PWGM_BNDRYDATA *pBndryData)
{
    return (BLKMODEL && BLKMODEL->PwBlkNdxBoundaryCB && pBndryData) ?
        BLKMODEL->PwBlkNdxBoundaryCB(block, ndx, pBndryData) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwBlkNdxBoundaryAndCondition(PWGM_HBLOCK block, PWP_UINT32 ndx,
    PWGM_BNDRYDATA *pBndryData, PWGM_CONDDATA *pCondData)
{
    return (BLKMODEL && BLKMODEL->PwBlkNdxBoundaryAndConditionCB &&
            (pBndryData || pCondData)) ?
        BLKMODEL->PwBlkNdxBoundaryAndConditionCB(block, ndx, pBndryData,
            pCondData) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwBlkNdxConnection(PWGM_HBLOCK block, PWP_UINT32 ndx, PWGM_CNXNDATA *pCnxnData)
{
    return (BLKMODEL && BLKMODEL->PwBlkNdxConnectionCB && pCnxnData) ?
        BLKMODEL->PwBlkNdxConnectionCB(block, ndx, pCnxnData) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwBoundary(PWGM_HBNDRY boundary, PWGM_BNDRYDATA *pBndryData)
{
    return (BNDRYMODEL && BNDRYMODEL->PwBoundaryCB && pBndryData) ?
        BNDRYMODEL->PwBoundaryCB(boundary, pBndryData) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwBndryCondition(PWGM_HBNDRY boundary, PWGM_CONDDATA *pCondData)
{
    return (BNDRYMODEL && BNDRYMODEL->PwBndryConditionCB) ?
        BNDRYMODEL->PwBndryConditionCB(boundary, pCondData) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwConnection(PWGM_HCNXN connection, PWGM_CNXNDATA *pCnxnData)
{
    return (CNXNMODEL && CNXNMODEL->PwConnectionCB && pCnxnData) ?
        CNXNMODEL->PwConnectionCB(connection, pCnxnData) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwBlock(PWGM_HBLOCK block, PWGM_BLOCKDATA *pBlockData)
{
    return (BLKMODEL && BLKMODEL->PwBlockCB && pBlockData) ?
        BLKMODEL->PwBlockCB(block, pBlockData) : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwXform2to3(const PWGM_INDEX_XFORM2 *pX2, PWGM_INDEX_XFORM *pX3)
{
    // XFORM2       XFORM
    // ------       -------
    // A B C        A B 0 C
    // D E F   ==>  D E 0 F
    //              0 0 1 0
    PWP_BOOL ret = (pX2 && pX3) ? PWP_TRUE : PWP_FALSE;
    if (ret) {
        // row 0
        pX3->m[0][0] = pX2->m[0][0];
        pX3->m[0][1] = pX2->m[0][1];
        pX3->m[0][2] = 0;
        pX3->m[0][0] = pX2->m[0][2];
        // row 1
        pX3->m[1][0] = pX2->m[1][0];
        pX3->m[1][1] = pX2->m[1][1];
        pX3->m[1][2] = 0;
        pX3->m[1][0] = pX2->m[1][2];
        // row 2
        pX3->m[2][0] = 0;
        pX3->m[2][1] = 0;
        pX3->m[2][2] = 1;
        pX3->m[2][0] = 0;
    }
    return ret;
}

/**************************************/
PWP_BOOL
PwXform3to2(const PWGM_INDEX_XFORM *pX3, PWGM_INDEX_XFORM2 *pX2)
{
    // XFORM         XFORM2
    // -------       ------
    // A B C D       A B D
    // E F G H  ==>  E F H
    // I J K L       
    PWP_BOOL ret = (pX2 && pX3) ? PWP_TRUE : PWP_FALSE;
    if (ret) {
        // row 0
        pX2->m[0][0] = pX3->m[0][0];
        pX2->m[0][1] = pX3->m[0][1];
        pX2->m[0][2] = pX3->m[0][3];
        // row 1
        pX2->m[1][0] = pX3->m[1][0];
        pX2->m[1][1] = pX3->m[1][1];
        pX2->m[1][2] = pX3->m[1][3];
    }
    return ret;
}

/**************************************/
PWGM_INDEX3
PwXformApply(const PWGM_INDEX_XFORM *pX3, PWGM_INDEX3 ijk)
{
    PWGM_INDEX3 ret;
    ret.i = pX3->m[0][0] * ijk.i + pX3->m[0][1] * ijk.j + pX3->m[0][2] * ijk.k
                + pX3->m[0][3];
    ret.j = pX3->m[1][0] * ijk.i + pX3->m[1][1] * ijk.j + pX3->m[1][2] * ijk.k
                + pX3->m[1][3];
    ret.k = pX3->m[2][0] * ijk.i + pX3->m[2][1] * ijk.j + pX3->m[2][2] * ijk.k
                + pX3->m[2][3];
    return ret;
}

/**************************************/
PWGM_ENUM_IJK
PwXformFollows(const PWGM_INDEX_XFORM *pX3, PWGM_ENUM_IJK localAxis,
    PWP_BOOL *pFlipped)
{
    PWGM_ENUM_IJK transformedAxis = PWGM_IJK_SIZE; /* invalid */
    PWGM_INDEX3 ijk = { 0, 0, 0 };
    PWGM_INDEX3 ijk2 = { 0, 0, 0 };
    if (pFlipped) {
        *pFlipped = PWP_FALSE;
    }
    /* construct axis vec for localAxis that will be rotated below */
    switch (localAxis) {
    case PWGM_IJK_I: ijk.i = 1; break;
    case PWGM_IJK_J: ijk.j = 1; break;
    case PWGM_IJK_K: ijk.k = 1; break;
    }
    /* apply rotation only */
    ijk2.i = pX3->m[0][0] * ijk.i + pX3->m[0][1] * ijk.j + pX3->m[0][2] * ijk.k;
    ijk2.j = pX3->m[1][0] * ijk.i + pX3->m[1][1] * ijk.j + pX3->m[1][2] * ijk.k;
    ijk2.k = pX3->m[2][0] * ijk.i + pX3->m[2][1] * ijk.j + pX3->m[2][2] * ijk.k;
    /* determine which axis it is in the transformed coord system */
    if ((0 != ijk2.i) && (0 == ijk2.j) && (0 == ijk2.k)) {
        transformedAxis = PWGM_IJK_I;
        if (pFlipped) {
            *pFlipped = ((1 == ijk2.i) ? PWP_FALSE : PWP_TRUE);
        }
    }
    else if ((0 != ijk2.j) && (0 == ijk2.i) && (0 == ijk2.k)) {
        transformedAxis = PWGM_IJK_J;
        if (pFlipped) {
            *pFlipped = ((1 == ijk2.j) ? PWP_FALSE : PWP_TRUE);
        }
    }
    else if ((0 != ijk2.k) && (0 == ijk2.i) && (0 == ijk2.j)) {
        transformedAxis = PWGM_IJK_K;
        if (pFlipped) {
            *pFlipped = ((1 == ijk2.k) ? PWP_FALSE : PWP_TRUE);
        }
    }
    return transformedAxis;
}

/**************************************/
PWGM_INDEX3
PwXform2Apply(const PWGM_INDEX_XFORM2 *pX2, PWGM_INDEX3 ijk)
{
    PWGM_INDEX3 ret;
    ret.i = pX2->m[0][0] * ijk.i + pX2->m[0][1] * ijk.j + pX2->m[0][2];
    ret.j = pX2->m[1][0] * ijk.i + pX2->m[1][1] * ijk.j + pX2->m[1][2];
    ret.k = ijk.k;
    return ret;
}

/**************************************/
PWGM_ENUM_IJK
PwXform2Follows(const PWGM_INDEX_XFORM2 *pX2, PWGM_ENUM_IJK localAxis,
    PWP_BOOL *pFlipped)
{
    PWGM_ENUM_IJK ret = PWGM_IJK_SIZE;
    PWGM_INDEX_XFORM x3;
    if (PwXform2to3(pX2, &x3)) {
        ret = PwXformFollows(&x3, localAxis, pFlipped);
    }
    return ret;
}

/**************************************/
PWP_BOOL
PwInRange(PWGM_INDEX3 ijk, const PWGM_STR_RANGE *pRange)
{
    return (ijk.i >= pRange->beg.i) && (ijk.i <= pRange->end.i) &&
           (ijk.j >= pRange->beg.j) && (ijk.j <= pRange->end.j) &&
           (ijk.k >= pRange->beg.k) && (ijk.k <= pRange->end.k);
}

#endif // !defined(PWGM_HIDE_STRUCTURED_API)
