/****************************************************************************
 *
 * CAEP Plugin example
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#include <stdio.h>
#include <stddef.h>
#include <time.h>

#define SHOW_PWP_MESSAGES
#include "apiCAEP.h"

#include "apiCAEPUtils.h"
#include "apiPWPUtils.h"
#include "apiUtils.h"
#include "runtimeWrite.h"


//************************************************
// impl-defined CAE format support data
//************************************************
#include "rtCaepSupportData.h"

/*------------------------------------*/
CAEP_RTITEM caepRtItem[] = {
    //************************************************
    // impl-defined CAE format data
    //************************************************
#   include "rtCaepInitItems.h"
};
PWP_UINT32 caepFormatCnt = ARRAYSIZE(caepRtItem);

/**************************************/
/* Safely casts the incoming CAEP_EXPORTER handle to a CAEP_RTITEM*.
   NOTE: caeuH2Rti() needs to be here to gain compile-time access
   to caepRtItem array.
*/
static CAEP_RTITEM *
caeuH2Rti(CAEP_EXPORTER handle)
{
    CAEP_RTITEM *ret = 0;
    if (handle) {
        ptrdiff_t diff = ((char*)handle) - (char*)caepRtItem;
        ptrdiff_t mod = diff % sizeof(caepRtItem[0]);
        if ((0 == mod) && (0 <= diff) &&
            (diff < (ptrdiff_t)sizeof(caepRtItem))) {
            ret = ((CAEP_RTITEM*)handle);
        }
    }
    return ret;
}


/**************************************/
CAEP_EXPORTER
PwCreateCaeById(PWP_UINT32 id)
{
    CAEP_RTITEM* ret = caeuFindFormatById(id);
    return (CAEP_EXPORTER)ret;
}

/**************************************/
CAEP_EXPORTER
PwCreateCaeByName(const char name[])
{
    CAEP_RTITEM* ret = caeuFindFormatByName(name);
    return (CAEP_EXPORTER)ret;
}

/**************************************/
PWP_VOID
PwDestroyCae(CAEP_EXPORTER *handle)
{
    if (handle && caeuH2Rti(*handle)) {
        *handle = 0;
    }
}

/**************************************/
const char*
PwEnumCaeFormat(PWP_UINT32 ndx, CAEP_FORMATINFO *pFormatInfo)
{
    const char* ret = (ndx < caepFormatCnt) ? caepRtItem[ndx].FormatInfo.name : 0;
    if (ret && pFormatInfo) {
        *pFormatInfo = caepRtItem[ndx].FormatInfo;
    }
    return ret;
}

/**************************************/
PWP_UINT32
PwGetCaeFormatCount()
{
    return caepFormatCnt;
}

/**************************************/
const char*
PwCaeFormat(CAEP_EXPORTER handle, CAEP_FORMATINFO *pFormatInfo)
{
    CAEP_RTITEM *pRti = caeuH2Rti(handle);
    const char* ret = pRti ? pRti->FormatInfo.name : 0;
    if (ret && pFormatInfo) {
        *pFormatInfo = pRti->FormatInfo;
    }
    return ret;
}

/**************************************/
PWP_BOOL
PwCaeElementType(CAEP_EXPORTER handle, PWGM_ENUM_ELEMTYPE which)
{
    PWP_BOOL ret = PWP_FALSE;
    CAEP_RTITEM *pRti = caeuH2Rti(handle);
    if (pRti) {
        switch (which) {
            case PWGM_ELEMTYPE_BAR:
            case PWGM_ELEMTYPE_HEX:
            case PWGM_ELEMTYPE_QUAD:
            case PWGM_ELEMTYPE_TRI:
            case PWGM_ELEMTYPE_TET:
            case PWGM_ELEMTYPE_WEDGE:
            case PWGM_ELEMTYPE_PYRAMID: ret = pRti->elemType[which]; break;
        }
    }
    return ret;
}

/**************************************/
const char*
PwCaeEnumBCs(CAEP_EXPORTER handle, PWP_UINT32 ndx, CAEP_BCINFO *pBCInfo)
{
    const char* ret = 0;
    CAEP_RTITEM *pRti = caeuH2Rti(handle);
    if (pRti) {
        ret = (ndx < pRti->BCCnt) ? pRti->pBCInfo[ndx].phystype : 0;
        if (ret && pBCInfo) {
            *pBCInfo = pRti->pBCInfo[ndx];
        }
    }
    return ret;
}

/**************************************/
const char*
PwCaeEnumFileExt(CAEP_EXPORTER handle, PWP_UINT32 ndx)
{
    const char* ret = 0;
    CAEP_RTITEM *pRti = caeuH2Rti(handle);
    ret = (pRti && (ndx < pRti->ExtCnt)) ? pRti->pFileExt[ndx] : 0;
    return ret;
}

/**************************************/
const char*
PwCaeEnumVCs(CAEP_EXPORTER handle, PWP_UINT32 ndx, CAEP_VCINFO *pVCInfo)
{
    const char* ret = 0;
    CAEP_RTITEM *pRti = caeuH2Rti(handle);
    if (pRti && pRti->FormatInfo.allowedVolumeConditions) {
        ret = (ndx < pRti->VCCnt) ? pRti->pVCInfo[ndx].phystype : 0;
        if (ret && pVCInfo) {
            *pVCInfo = pRti->pVCInfo[ndx];
        }
    }
    return ret;
}

/**************************************/
PWP_UINT32
PwCaeGetBCCount(CAEP_EXPORTER handle)
{
    PWP_UINT32 ret = 0;
    CAEP_RTITEM *pRti = caeuH2Rti(handle);
    if (pRti) {
        ret = pRti->BCCnt;
    }
    return ret;
}

/**************************************/
PWP_UINT32
PwCaeGetFileExtCount(CAEP_EXPORTER handle)
{
    PWP_UINT32 ret = 0;
    CAEP_RTITEM *pRti = caeuH2Rti(handle);
    if (pRti) {
        ret = pRti->ExtCnt;
    }
    return ret;
}

/**************************************/
PWP_UINT32
PwCaeGetVCCount(CAEP_EXPORTER handle)
{
    PWP_UINT32 ret = 0;
    CAEP_RTITEM *pRti = caeuH2Rti(handle);
    if (pRti && pRti->FormatInfo.allowedVolumeConditions) {
        ret = pRti->VCCnt;
    }
    return ret;
}

/**************************************/
PWP_BOOL
PwCaeGridWrite(CAEP_EXPORTER handle, PWGM_HGRIDMODEL model,
               const CAEP_WRITEINFO *pWriteInfo)
{
    PWP_BOOL ret = PWP_FALSE;
    CAEP_RTITEM *pRti = caeuH2Rti(handle);
    if (model && caeuFileOpen(pRti, pWriteInfo)) {
        pRti->model = model;
        pRti->pWriteInfo = pWriteInfo;
        // give impl control!
        ret = runtimeWrite(pRti, model, pWriteInfo);
        // we opened it! we close it!
        caeuFileClose(pRti, pWriteInfo);
        // force cwd back to original location
        while (0 == pwpCwdPop());
    }
    return ret;
}
