/****************************************************************************
 *
 * CAEP Plugin example - PwCaeGridWrite implementation
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

//#include <stdio.h>
//#include <time.h>
//
#include "apiCAEP.h"
#include "apiCAEPUtils.h"
#include "apiPWP.h"
#include "runtimeWrite.h"


// There are 5 predefined entity references in XML:
//
// &lt; 	< 	less than
// &gt; 	> 	greater than
// &amp; 	& 	ampersand 
// &apos; 	' 	apostrophe
// &quot; 	" 	quotation mark
//
// Note: Only the characters "<" and "&" are strictly illegal in XML. The
// greater than character is legal, but it is a good habit to replace it. 

//------------------------------------------------
static void
writeComment(CAEP_RTITEM *pRti, const char *pComment)
{
    if (pRti) {
        fprintf(pRti->fp, "<!-- %s -->\n", (pComment? pComment : ""));
    }
}

//------------------------------------------------
static void
writeCondData(CAEP_RTITEM *pRti, PWGM_CONDDATA *pCondData)
{
    if (pRti && pCondData) {
        fputs("\t\t\t<condition>\n", pRti->fp);
        fprintf(pRti->fp, "\t\t\t\t<physical id='%ld' name='%s' />\n",
            (long)pCondData->tid, pCondData->type);
        fprintf(pRti->fp, "\t\t\t\t<user id='%ld' name='%s' />\n",
            (long)pCondData->id, pCondData->name);
        fputs("\t\t\t</condition>\n", pRti->fp);
    }
}

//------------------------------------------------
static const char *
elemType2Str(PWGM_ENUM_ELEMTYPE type)
{
    const char * ret = 0;
    switch (type) {
        case PWGM_ELEMTYPE_BAR:     ret = "bar"; break;
        case PWGM_ELEMTYPE_HEX:     ret = "hex"; break;
        case PWGM_ELEMTYPE_QUAD:    ret = "quad"; break;
        case PWGM_ELEMTYPE_TRI:     ret = "tri"; break;
        case PWGM_ELEMTYPE_TET:     ret = "tet"; break;
        case PWGM_ELEMTYPE_WEDGE:   ret = "wedge"; break;
        case PWGM_ELEMTYPE_PYRAMID: ret = "pyramid"; break;
        default: ret = "ERROR"; break;
    }
    return ret;
}

//------------------------------------------------
static void
writeElemData(CAEP_RTITEM *pRti, PWGM_ELEMDATA *pElemData)
{
    if (pRti && pElemData) {
        PWP_UINT32 ii;
        fprintf(pRti->fp, "\t\t\t<element type='%s' count='%lu'>\n",
            elemType2Str(pElemData->type), (unsigned long)pElemData->vertCnt);
        for (ii=0; ii < pElemData->vertCnt; ++ii) {
            fprintf(pRti->fp, "\t\t\t\t<index type='vertex' val='%lu' />\n",
                (unsigned long)pElemData->index[ii]);
        }
        fputs("\t\t\t</element>\n", pRti->fp);
    }
}

//------------------------------------------------
static void
writeElemCounts(CAEP_RTITEM *pRti, PWGM_ELEMCOUNTS *pCnts)
{
    if (pRti && pCnts) {
        fputs("\t\t\t<namedcounts type='elementtype'>\n", pRti->fp);
        fprintf(pRti->fp, "\t\t\t\t<count name='bar' val='%lu' />\n",
            (unsigned long)PWGM_ECNT_Bar(*pCnts));
        fprintf(pRti->fp, "\t\t\t\t<count name='hex' val='%lu' />\n",
            (unsigned long)PWGM_ECNT_Hex(*pCnts));
        fprintf(pRti->fp, "\t\t\t\t<count name='quad' val='%lu' />\n",
            (unsigned long)PWGM_ECNT_Quad(*pCnts));
        fprintf(pRti->fp, "\t\t\t\t<count name='tri' val='%lu' />\n",
            (unsigned long)PWGM_ECNT_Tri(*pCnts));
        fprintf(pRti->fp, "\t\t\t\t<count name='tet' val='%lu' />\n",
            (unsigned long)PWGM_ECNT_Tet(*pCnts));
        fprintf(pRti->fp, "\t\t\t\t<count name='wedge' val='%lu' />\n",
            (unsigned long)PWGM_ECNT_Wedge(*pCnts));
        fprintf(pRti->fp, "\t\t\t\t<count name='pyramid' val='%lu' />\n",
            (unsigned long)PWGM_ECNT_Pyramid(*pCnts));
        fputs("\t\t\t</namedcounts>\n", pRti->fp);
    }
}

//------------------------------------------------
static int
writeBlock(CAEP_RTITEM *pRti, PWGM_HBLOCK hBlk, PWP_BOOL condOnly)
{
    int ret = 0;
    if (pRti && PWGM_HBLOCK_ISVALID(hBlk)) {
        PWGM_ELEMDATA eData;
        PWGM_CONDDATA CondData;
        PWGM_ELEMCOUNTS elemCnts;
        PWP_UINT32 eCnt = PwBlkElementCount(hBlk, &elemCnts);

        if (caeuProgressBeginStep(pRti, condOnly ? 1 : eCnt)) {
            fprintf(pRti->fp, "\t\t<block id='%lu' count='%lu'>\n",
                (unsigned long)PWGM_HDOMAIN_ID(hBlk), (unsigned long)eCnt);
            if (PwBlkCondition(hBlk, &CondData)) {
                writeCondData(pRti, &CondData);
            }
            if (condOnly) {
                caeuProgressIncr(pRti);
            }
            else {
                writeElemCounts(pRti, &elemCnts);
                eCnt = 0;
                while (PwElemDataMod(PwBlkEnumElements(hBlk, eCnt++), &eData)) {
                    writeElemData(pRti, &eData);
                    if (!caeuProgressIncr(pRti)) {
                        break;
                    }
                }
            }
            fputs("\t\t</block>\n", pRti->fp);
            caeuProgressEndStep(pRti);
            ret = !pRti->opAborted;
        }
    }
    return ret;
}

//------------------------------------------------
static void
writeBlocks(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model, PWP_BOOL condOnly)
{
    if (!pRti->opAborted) {
        PWP_UINT32 cnt = PwModBlockCount(model);
        writeComment(pRti,
            "\n**************************************************************\n"
            "\tGRID BLOCKS\n"
            "**************************************************************\n"
        );
        fprintf(pRti->fp, "\t<blocks count='%lu'>\n", (unsigned long)cnt);
        cnt = 0;
        while (writeBlock(pRti, PwModEnumBlocks(model, cnt++), condOnly)) {
        }
        fputs("\t</blocks>\n\n", pRti->fp);
    }
}

//------------------------------------------------
static int
writeDomain(CAEP_RTITEM *pRti, PWGM_HDOMAIN hDom, PWP_BOOL condOnly)
{
    int ret = 0;
    if (pRti && PWGM_HDOMAIN_ISVALID(hDom)) {
        PWGM_ELEMDATA eData;
        PWGM_CONDDATA CondData;
        PWGM_ELEMCOUNTS elemCnts;
        PWP_UINT32 eCnt = PwDomElementCount(hDom, &elemCnts);

        if (caeuProgressBeginStep(pRti, condOnly ? 1 : eCnt)) {
            fprintf(pRti->fp, "\t\t<domain id='%lu' count='%lu'>\n",
                (unsigned long)PWGM_HDOMAIN_ID(hDom), (unsigned long)eCnt);
            if (PwDomCondition(hDom, &CondData)) {
                writeCondData(pRti, &CondData);
            }
            if (condOnly) {
                caeuProgressIncr(pRti);
            }
            else {
                eCnt = 0;
                while (PwElemDataMod(PwDomEnumElements(hDom, eCnt++), &eData)) {
                    writeElemData(pRti, &eData);
                    if (!caeuProgressIncr(pRti)) {
                        break;
                    }
                }
            }
            fputs("\t\t</domain>\n", pRti->fp);
            caeuProgressEndStep(pRti);
            ret = !pRti->opAborted;
        }
    }
    return ret;
}

//------------------------------------------------
static void
writeDomains(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model, PWP_BOOL condOnly)
{
    if (!pRti->opAborted) {
        PWP_UINT32 cnt = PwModDomainCount(model);
        writeComment(pRti,
            "\n**************************************************************\n"
            "\tGRID DOMAINS\n"
            "**************************************************************\n"
        );
        fprintf(pRti->fp, "\t<domains count='%lu'>\n", (unsigned long)cnt);
        cnt = 0;
        while (writeDomain(pRti, PwModEnumDomains(model, cnt++), condOnly)) {
        }
        fputs("\t</domains>\n", pRti->fp);
    }
}

//------------------------------------------------
static int
writeVertex(CAEP_RTITEM *pRti, PWGM_HVERTEX vertex)
{
    int ret = 0;
    if (pRti && PWGM_HVERTEX_ISVALID(vertex)) {
        PWGM_VERTDATA VertData;
        int prec = CAEPU_RT_PREC_SINGLE(pRti) ? 8 : 16;
        PwVertDataMod(vertex, &VertData);
        if (CAEPU_RT_DIM_3D(pRti)) {
            fprintf(pRti->fp,
                "\t\t<vertex i='%lu' x='%.*g' y='%.*g' z='%.*g' />\n",
                (unsigned long)VertData.i, prec, VertData.x, prec, VertData.y,
                prec, VertData.z);
        }
        else {
            fprintf(pRti->fp, "\t\t<vertex i='%lu' x='%.*g' y='%.*g' />\n",
                (unsigned long)VertData.i, prec, VertData.x, prec, VertData.y);
        }
        ret = !pRti->opAborted;
    }
    return ret;
}

//------------------------------------------------
static void
writeVertices(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model, PWP_BOOL condOnly)
{
    if (!condOnly && !pRti->opAborted) {
        PWP_UINT32 cnt = PwModVertexCount(model);
        if (caeuProgressBeginStep(pRti, cnt)) {
            writeComment(pRti,
                "\n**************************************************************\n"
                "\tGRID VERTICES\n"
                "**************************************************************\n"
            );
            fprintf(pRti->fp, "\t<vertices count='%lu'>\n",
                (unsigned long)cnt);
            cnt=0;
            while (writeVertex(pRti, PwModEnumVertices(model, cnt++))) {
                caeuProgressIncr(pRti);
            }
            fputs("\t</vertices>\n", pRti->fp);
            caeuProgressEndStep(pRti);
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
        pRti->FormatInfo.group,
            pRti->FormatInfo.name,
            (unsigned long)pluginInfo.libVer.major,
            (unsigned long)pluginInfo.libVer.minor,
        strTime, //asctime(localtime(&szClock)),
        pRti->pApiData->apiInfo.name,
        (unsigned long)pRti->FormatInfo.id,
            (unsigned long)(pRti->FormatInfo.id >> 16) & 0x0000FFFFUL,
            (unsigned long)(pRti->FormatInfo.id & 0x0000FFFFUL),
            (unsigned long)pRti->FormatInfo.id,
        pluginInfo.author, pluginInfo.copyright,
        pluginInfo.support,
        pRti->pWriteInfo->fileDest,
        pRti->pWriteInfo->conditionsOnly ? "true" : "false",
        CAEPU_RT_ENCODING_TEXT(pRti) ,
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

    if (pRti && model) {
        // Treat all verts as 1 step. If condsOnly, skip verts
        PWP_UINT32 cnt = (pWriteInfo->conditionsOnly ? 0 : 1)
                       + PwModBlockCount(model)
                       + PwModDomainCount(model);

        if (caeuProgressInit(pRti, cnt)) {
            fputs("<?xml version='1.0' encoding='UTF-8' ?>\n", pRti->fp);
            writeInfo(pRti);
            fprintf(pRti->fp, "<pointwise dim='%s' type='unstructured'>\n",
                    CAEPU_RT_DIM_TEXT(pRti));
            writeVertices(pRti, model, pWriteInfo->conditionsOnly);
            writeDomains(pRti, model, pWriteInfo->conditionsOnly);
            writeBlocks(pRti, model, pWriteInfo->conditionsOnly);
            fputs("</pointwise>\n", pRti->fp);
            caeuFileClose(pRti, pWriteInfo);
            caeuProgressEnd(pRti, ret);
            ret = !pRti->opAborted;
        }
    }
    return ret;
}
