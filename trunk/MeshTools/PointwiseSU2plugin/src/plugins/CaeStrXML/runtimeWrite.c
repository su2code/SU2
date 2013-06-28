/****************************************************************************
 *
 * CAEP Plugin example - PwCaeGridWrite implementation
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#include "apiCAEP.h"
#include "apiCAEPUtils.h"
#include "apiGridModel.h"
#include "apiPWP.h"
#include "runtimeWrite.h"
#include "pwpPlatform.h"

#ifdef __cplusplus
} /* extern "C" */
#endif


const char *faceIdStr[] = {
    "KMIN",  /* PWGM_FACE_KMIN */
    "KMAX",  /* PWGM_FACE_KMAX */
    "IMIN",  /* PWGM_FACE_IMIN */
    "IMAX",  /* PWGM_FACE_IMAX */
    "JMIN",  /* PWGM_FACE_JMIN */
    "JMAX",  /* PWGM_FACE_JMAX */
    "UNSTR"  /* PWGM_FACE_UNSTR */
};

//------------------------------------------------
static void
writeComment(CAEP_RTITEM *pRti, const char *pComment)
{
    if (pRti && pRti->fp) {
        fprintf(pRti->fp, "<!-- %s -->\n", (pComment? pComment : ""));
    }
}

//------------------------------------------------
static void
writeIndex3(CAEP_RTITEM *pRti, const PWGM_INDEX3 *pNdx)
{
    if (CAEPU_RT_DIM_3D(pRti)) {
        fprintf(pRti->fp, "<index i='%lu' j='%lu' k='%lu' />\n",
            (unsigned long)pNdx->i, (unsigned long)pNdx->j,
            (unsigned long)pNdx->k);
    }
    else {
        fprintf(pRti->fp, "<index i='%lu' j='%lu' />\n",
            (unsigned long)pNdx->i, (unsigned long)pNdx->j);
    }
}

//------------------------------------------------
static void
writeRange(CAEP_RTITEM *pRti, const PWGM_STR_RANGE *pRange)
{
    if (CAEPU_RT_DIM_3D(pRti)) {
        fprintf(pRti->fp,
            "<range imin='%lu' jmin='%lu' kmin='%lu' "
                   "imax='%lu' jmax='%lu' kmax='%lu' />\n",
            (unsigned long)pRange->beg.i, (unsigned long)pRange->beg.j,
            (unsigned long)pRange->beg.k,
            (unsigned long)pRange->end.i, (unsigned long)pRange->end.j,
            (unsigned long)pRange->end.k);
    }
    else {
        fprintf(pRti->fp,
            "<range imin='%lu' jmin='%lu' imax='%lu' jmax='%lu' />\n",
            (unsigned long)pRange->beg.i, (unsigned long)pRange->beg.j,
            (unsigned long)pRange->end.i, (unsigned long)pRange->end.j);
    }
}

//------------------------------------------------
static void
writeTransform(CAEP_RTITEM *pRti, const PWGM_INDEX_XFORM *pXform)
{
    /* PWP_INT m[3][4];    m[row][col] */
    int row;
    int col;
    fputs("<xform", pRti->fp);
    if (CAEPU_RT_DIM_3D(pRti)) {
        for (row=0; row < 3; ++row) {
            if (row) {
                fputs("\n      ", pRti->fp);
            }
            for (col=0; col < 4; ++col) {
                fprintf(pRti->fp, " m%u%u='%d'", row, col,
                    pXform->m[row][col]);
            }
        }
    }
    else {
        PWGM_INDEX_XFORM2 x2;
        if (PwXform3to2(pXform, &x2)) {
            for (row=0; row < 2; ++row) {
                if (row) {
                    fputs("\n      ", pRti->fp);
                }
                for (col=0; col < 3; ++col) {
                    fprintf(pRti->fp, " m%u%u='%d'", row, col, x2.m[row][col]);
                }
            }
        }
    }
    fputs(" />\n", pRti->fp);
}

//------------------------------------------------
static PWP_VOID
writeLink(CAEP_RTITEM *pRti, PWP_BOOL fwd, PWGM_HBLOCK block,
    PWGM_FACE_ID face, const PWGM_STR_RANGE *pRange,
    const PWGM_INDEX_XFORM *pXform)
{
    fprintf(pRti->fp, "<link dir='%s' blk='%lu' face='%s'>\n",
        (fwd ? "fwd" : "rev"), (unsigned long)PWGM_HBLOCK_ID(block),
        faceIdStr[face]);
    writeRange(pRti, pRange);
    writeTransform(pRti, pXform);
    fputs("</link>\n", pRti->fp);
}

//------------------------------------------------
static PWP_VOID
writeConnectionData(CAEP_RTITEM *pRti, const PWGM_CNXNDATA *pCnxnData)
{
    if (pCnxnData) {
        const char *name = (pCnxnData->name ? pCnxnData->name : "");
        fprintf(pRti->fp, "<cnxn name='%s'>\n", name);
        writeLink(pRti, PWP_TRUE, pCnxnData->block1, pCnxnData->face1,
                  &pCnxnData->range1, &pCnxnData->from1to2);
        writeLink(pRti, PWP_FALSE, pCnxnData->block2, pCnxnData->face2,
                  &pCnxnData->range2, &pCnxnData->from2to1);
        fputs("</cnxn>\n", pRti->fp);
    }
}

//------------------------------------------------
static PWP_BOOL
writeConnectionNdx(CAEP_RTITEM *pRti, PWP_UINT32 ndx)
{
    PWGM_CNXNDATA cnxnData;
    PWP_BOOL ret = PwModNdxConnection(pRti->model, ndx, &cnxnData);
    if (ret) {
        writeConnectionData(pRti, &cnxnData);
    }
    return ret;
}

//------------------------------------------------
static PWP_VOID
writeConditionData(CAEP_RTITEM *pRti, const PWGM_CONDDATA *pCondData)
{
    if (pCondData) {
        fprintf(pRti->fp,
            "<condition name='%s' id='%lu' type='%s' typeid='%lu' />\n",
            pCondData->name, (unsigned long)pCondData->id, pCondData->type,
            (unsigned long)pCondData->tid);
    }
}

//------------------------------------------------
static PWP_VOID
writeBoundaryAndConditionData(CAEP_RTITEM *pRti,
    const PWGM_BNDRYDATA *pBndryData, const PWGM_CONDDATA *pCondData)
{
    if (pBndryData) {
        const char *name = (pBndryData->name ? pBndryData->name : "");
        fprintf(pRti->fp, "<boundary name='%s' blk='%lu' face='%s'>\n", name,
            (unsigned long)PWGM_HBLOCK_ID(pBndryData->block),
            faceIdStr[pBndryData->face]);
        writeRange(pRti, &pBndryData->range);
        writeConditionData(pRti, pCondData);
        fputs("</boundary>\n", pRti->fp);
    }
}

//------------------------------------------------
static PWP_BOOL
writeBoundaryNdx(CAEP_RTITEM *pRti, PWP_UINT32 ndx)
{
    PWGM_BNDRYDATA bndryData;
    PWGM_CONDDATA condData;
    PWP_BOOL ret = PwModNdxBoundaryAndCondition(pRti->model, ndx, &bndryData, &condData);
    if (ret) {
        writeBoundaryAndConditionData(pRti, &bndryData, &condData);
    }
    return ret;
}

//------------------------------------------------
static void
writeVertexData(CAEP_RTITEM *pRti, const PWGM_INDEX3 *pIjk,
    const PWGM_VERTDATA *pVertData)
{
    if (CAEPU_RT_DIM_3D(pRti)) {
        fprintf(pRti->fp,
            "<vertex x='%.5g' y='%.5g' z='%.5g' i='%lu' j='%lu' k='%lu' />\n",
            pVertData->x, pVertData->y, pVertData->z,
            (unsigned long)pIjk->i, (unsigned long)pIjk->j,
            (unsigned long)pIjk->k);
    }
    else {
        fprintf(pRti->fp, "<vertex x='%.5g' y='%.5g' i='%lu' j='%lu' />\n",
            pVertData->x, pVertData->y,
            (unsigned long)pIjk->i,
            (unsigned long)pIjk->j);
    }
}

//------------------------------------------------
static PWP_BOOL
writeBlkVertices(CAEP_RTITEM *pRti, PWGM_HBLOCK block)
{
    PWGM_STR_SIZE size;
    PWP_BOOL ret = PwBlkSize(block, &size);
    if (ret) {
        PWGM_VERTDATA vertData;
        PWGM_INDEX3 ijk;
        fprintf(pRti->fp, "<vertices cnt='%lu'>\n",
            (unsigned long)(size.i * size.j * size.k));
        for (ijk.i=0; ret && (ijk.i < size.i); ++ijk.i) {
            for (ijk.j=0; ret && (ijk.j < size.j); ++ijk.j) {
                for (ijk.k=0; ret && (ijk.k < size.k); ++ijk.k) {
                    ret = PwBlkNdxVertData(block, ijk, &vertData) &&
                          !pRti->opAborted;
                    if (ret) {
                        writeVertexData(pRti, &ijk, &vertData);
                    }
                }
            }
        }
        fputs("</vertices>\n", pRti->fp);
    }
    return ret;
}

//------------------------------------------------
static PWP_BOOL
writeBlkConnections(CAEP_RTITEM *pRti, PWGM_HBLOCK block)
{
    PWP_BOOL ret = PWP_TRUE;

    PWP_UINT32 ndx = 0;
    PWP_UINT32 cnt = PwBlkConnectionCount(block);
    PWGM_CNXNDATA cnxnData;
    fprintf(pRti->fp, "<connections cnt='%lu'>\n", (unsigned long)cnt);
    while (PwBlkNdxConnection(block, ndx++, &cnxnData)) {
        writeConnectionData(pRti, &cnxnData);
    }
    fputs("</connections>\n", pRti->fp);
    return ret;
}

//------------------------------------------------
static PWP_BOOL
writeBlkBoundaries(CAEP_RTITEM *pRti, PWGM_HBLOCK block)
{
    PWP_BOOL ret = PWP_TRUE;
    PWP_UINT32 ndx = 0;
    PWGM_BNDRYDATA bndryData;
    PWGM_CONDDATA condData;
    fprintf(pRti->fp, "<boundaries cnt='%lu'>\n",
        (unsigned long)PwBlkBoundaryCount(block));
    while (PwBlkNdxBoundaryAndCondition(block, ndx++, &bndryData, &condData)) {
        writeBoundaryAndConditionData(pRti, &bndryData, &condData);
    }
    fputs("</boundaries>\n", pRti->fp);
    return ret;
}

//------------------------------------------------
static PWP_BOOL
writeBlock(CAEP_RTITEM *pRti, PWGM_HBLOCK block)
{
    PWP_BOOL ret = PWP_FALSE;
    PWGM_STR_SIZE size;
    PWGM_BLOCKDATA blkData;
    if (pRti && PWGM_HBLOCK_ISVALID(block) && PwBlkSize(block, &size) &&
            PwBlock(block, &blkData)) {
        if (0 == blkData.name) {
            blkData.name = "";
        }
        if (CAEPU_RT_DIM_3D(pRti)) {
            fprintf(pRti->fp,
                "<block name='%s' id='%lu' i='%lu' j='%lu' k='%lu'>\n",
                blkData.name,
                (unsigned long)PWGM_HBLOCK_ID(block), (unsigned long)size.i,
                (unsigned long)size.j, (unsigned long)size.k);
        }
        else {
            fprintf(pRti->fp,
                "<block name='%s' id='%lu' i='%lu' j='%lu'>\n",
                blkData.name, (unsigned long)PWGM_HBLOCK_ID(block),
                (unsigned long)size.i, (unsigned long)size.j);
        }
        ret = !pRti->opAborted && writeBlkVertices(pRti, block);
        ret = !pRti->opAborted && writeBlkConnections(pRti, block);
        ret = !pRti->opAborted && writeBlkBoundaries(pRti, block);
        fputs("</block>\n", pRti->fp);
    }
    return ret;
}

/**************************************/
static void
stepWriteModelBlocks(CAEP_RTITEM *pRti)
{
    PWP_UINT32 cnt = PwModBlockCount(pRti->model);
    if (caeuProgressBeginStep(pRti, cnt)) {
        PWP_UINT32 ndx = 0;
        PWGM_HBLOCK block = PwModEnumBlocks(pRti->model, ndx);
        fprintf(pRti->fp, "<blocks cnt='%lu'>\n", (unsigned long)cnt);
        while (!pRti->opAborted && writeBlock(pRti, block)) {
            caeuProgressIncr(pRti);
            block = PwModEnumBlocks(pRti->model, ++ndx);
        }
        fputs("</blocks>\n", pRti->fp);
        caeuProgressEndStep(pRti);
        if (ndx != cnt) {
            pRti->opAborted = PWP_TRUE;
        }
    }
}

/**************************************/
static void
stepWriteModelConnections(CAEP_RTITEM *pRti)
{
    /* write model-centric connection data
    */
    PWP_UINT32 cnt = PwModConnectionCount(pRti->model);
    if (caeuProgressBeginStep(pRti, cnt)) {
        PWP_UINT32 ndx = 0;
        fprintf(pRti->fp, "<connections cnt='%lu'>\n", (unsigned long)cnt);
        while (!pRti->opAborted && writeConnectionNdx(pRti, ndx)) {
            caeuProgressIncr(pRti);
            ++ndx;
        }
        fputs("</connections>\n", pRti->fp);
        caeuProgressEndStep(pRti);
        if (ndx != cnt) {
            pRti->opAborted = PWP_TRUE;
        }
    }
}

/**************************************/
static void
stepWriteModelBoundaries(CAEP_RTITEM *pRti)
{
    /* write model-centric boundary data
    */
    PWP_UINT32 cnt = PwModBoundaryCount(pRti->model);
    if (caeuProgressBeginStep(pRti, cnt)) {
        PWP_UINT32 ndx = 0;
        fprintf(pRti->fp, "<boundaries cnt='%lu'>\n", (unsigned long)cnt);
        while (!pRti->opAborted && writeBoundaryNdx(pRti, ndx)) {
            caeuProgressIncr(pRti);
            ++ndx;
        }
        fputs("</boundaries>\n", pRti->fp);
        caeuProgressEndStep(pRti);
        if (ndx != cnt) {
            pRti->opAborted = PWP_TRUE;
        }
    }
}

//------------------------------------------------
static void
writeInfo(CAEP_RTITEM *pRti)
{
    char buf[2048];
    char strTime[256];
    time_t szClock;
    PWP_PLUGININFO pluginInfo;
    time(&szClock);
    strftime(strTime, sizeof(strTime), "%Y-%m-%d %H:%M:%S", localtime(&szClock));
    PwpGetPluginInfo(&pluginInfo);
    sprintf(buf,
        "\n**************************************************************\n"
        "\tcreated by: \"%s %s\" CAE exporter plugin v%lu.%lu\n"
        "\tcreated on: %s\n"
        "\tspec: %s\n"
        "\tplugin id: %lu / %lu:%lu / 0x%08lX\n"
        "\tauthor: %s %s\n"
        "\tsupport: %s\n"
        "\tfile: %s\n"
        "\tbc only: %s\n"
        "\tencoding: %s\n"
        "\tprecision: %s\n"
        "\tdimension: %s\n"
        "**************************************************************\n",
        pRti->FormatInfo.group, pRti->FormatInfo.name,
            (unsigned long)pluginInfo.libVer.major,
            (unsigned long)pluginInfo.libVer.minor,
        strTime, //asctime(localtime(&szClock)),
        pRti->pApiData->apiInfo.name,
        (unsigned long)pRti->FormatInfo.id,
            (unsigned long)((pRti->FormatInfo.id >> 16) & 0x0000FFFFUL),
            (unsigned long)(pRti->FormatInfo.id & 0x0000FFFFUL),
            (unsigned long)pRti->FormatInfo.id,
        pluginInfo.author, pluginInfo.copyright,
        pluginInfo.support,
        pRti->pWriteInfo->fileDest,
        pRti->pWriteInfo->conditionsOnly ? "true" : "false",
        CAEPU_RT_ENCODING_TEXT(pRti),
        CAEPU_RT_PREC_TEXT(pRti),
        CAEPU_RT_DIM_TEXT(pRti));
    writeComment(pRti, buf);
}

/**************************************/
PWP_BOOL
runtimeWrite(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model,
             const CAEP_WRITEINFO *pWriteInfo)
{
    PWP_BOOL ret = PWP_FALSE;
    if (pRti && pWriteInfo && PWGM_HGRIDMODEL_ISVALID(model)) {
        PWP_UINT32 cnt = 3;  /* the # of MAJOR progress steps */
        if (caeuProgressInit(pRti, cnt)) {
            fputs("<?xml version='1.0' encoding='UTF-8' ?>\n", pRti->fp);
            writeInfo(pRti);
            fprintf(pRti->fp, "<pointwise dim='%s' type='structured'>\n",
                        CAEPU_RT_DIM_TEXT(pRti));
            stepWriteModelBlocks(pRti);
            stepWriteModelConnections(pRti);
            stepWriteModelBoundaries(pRti);
            fputs("</pointwise>\n", pRti->fp);
        }
        ret = !pRti->opAborted;
    }
    return ret;
}
