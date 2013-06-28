/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#include <string.h>

#include "apiPWP.h"
#include "apiCAEPUtils.h"
#include "pwpPlatform.h"


/**************************************/
CAEP_RTITEM*
caeuFindFormatById (PWP_UINT32 id)
{
    CAEP_RTITEM *ret = 0;
    PWP_UINT32 ii;
    for (ii=0; ii < caepFormatCnt; ++ii) {
        if (id == caepRtItem[ii].FormatInfo.id) {
            ret = &(caepRtItem[ii]);
            break;
        }
    }
    return ret;
}

/**************************************/
CAEP_RTITEM*
caeuFindFormatByName (const char name[])
{
    CAEP_RTITEM *ret = 0;
    PWP_UINT32 ii;
    for (ii=0; ii < caepFormatCnt; ++ii) {
        if (0 == strcmp(caepRtItem[ii].FormatInfo.name, name)) {
            ret = &(caepRtItem[ii]);
            break;
        }
    }
    return ret;
}

/**************************************/
PWP_BOOL
caeuProgressInit(CAEP_RTITEM *pRti, PWP_UINT32 cnt)
{
    PWP_BOOL ret = PWP_FALSE;
    if (pRti) {
        pRti->progTotal = 0;
        pRti->progComplete  = 0;
        pRti->clocks[CAEPU_CLKS_PROGUPDATE] =
            pRti->clocks[CAEPU_CLKS_PROGINIT] = clock();
        pRti->opAborted = PWP_FALSE;
        ret = PwuProgressBegin(pRti->pApiData->apiInfo.name, cnt);
        pRti->opAborted |= !ret;
    }
    return ret;
}

/**************************************/
PWP_BOOL
caeuProgressBeginStep(CAEP_RTITEM *pRti, PWP_UINT32 total)
{
    PWP_BOOL ret = PWP_FALSE;
    if (pRti && (total > 0)) {
        pRti->progComplete = 0;
        pRti->progTotal = total;
        ret = !pRti->opAborted;
        pRti->clocks[CAEPU_CLKS_BEGSTEP] = clock();
    }
    return ret;
}

/**************************************/
PWP_BOOL
caeuProgressIncr(CAEP_RTITEM *pRti)
{
    PWP_BOOL ret = PWP_FALSE;
    if (pRti && !pRti->opAborted) {
        // send update a max of 2 times per second
        const clock_t DELAY = (CLOCKS_PER_SEC / 2);
        pRti->clocks[CAEPU_CLKS_PROGINCR] = clock();
        ++pRti->progComplete;
        if (DELAY <= CAEPU_RT_CLKS_DIFF(pRti, CAEPU_CLKS_PROGUPDATE,
                        CAEPU_CLKS_PROGINCR)) {
            ret = PwuProgressStatus (pRti->pApiData->apiInfo.name,
                pRti->progComplete, pRti->progTotal);
            pRti->clocks[CAEPU_CLKS_PROGUPDATE] =
                pRti->clocks[CAEPU_CLKS_PROGINCR];
            pRti->opAborted |= !ret;
        }
        else {
            ret = PWP_TRUE;
        }
    }
    return ret;
}

/**************************************/
PWP_BOOL
caeuProgressEndStep(CAEP_RTITEM *pRti)
{
    PWP_BOOL ret = PWP_FALSE;
    if (pRti) {
        pRti->clocks[CAEPU_CLKS_ENDSTEP] = clock();
        pRti->progComplete = 0;
        pRti->progTotal = 0;
        ret = PwuProgressNextStep(pRti->pApiData->apiInfo.name);
        pRti->opAborted |= !ret;
    }
    return ret;
}

/**************************************/
void
caeuProgressEnd(CAEP_RTITEM *pRti, PWP_BOOL ok)
{
    if (pRti) {
        pRti->clocks[CAEPU_CLKS_PROGEND] = clock();
        pRti->progTotal = 0;
        pRti->progComplete  = 0;
        pRti->clocks[CAEPU_CLKS_PROGUPDATE] = 0;
        PwuProgressEnd(pRti->pApiData->apiInfo.name, ok);
    }
}

/**************************************/
int
caeuFileClose(CAEP_RTITEM* pRti, const CAEP_WRITEINFO *pWriteInfo)
{
    int ret = 0;
    if (pRti) {
        switch (pRti->FormatInfo.fileDest) {
            case PWP_FILEDEST_FILENAME:
                if (pRti->fp) {
                    if (0 == pwpFileClose(pRti->fp)) {
                        if (pRti->opAborted) {
                            pwpFileDelete(pWriteInfo->fileDest);
                        }
                        else {
                            ret = 1;
                        }
                    }
                    pRti->fp = 0;
                    if (0 != pwpCwdPop()) {
                        ret = 0;
                    }
                }
                break;

            case PWP_FILEDEST_BASENAME:
            case PWP_FILEDEST_FOLDER:
                /* restore cwd to the original folder.
                */
                if (0 == pwpCwdPop()) {
                    ret = 1;
                }
                break;
        }
    }
    return ret;
}

/*------------------------------------*/
static int
openFileName(CAEP_RTITEM* pRti, const CAEP_WRITEINFO *pWriteInfo)
{
    int ret = 0;
    if (pRti && pWriteInfo && pWriteInfo->fileDest && pWriteInfo->fileDest[0]) {
        caeuFileClose(pRti, pWriteInfo); // sanity check
        switch (pWriteInfo->encoding) {
            case PWP_ENCODING_ASCII:
                pRti->fp = pwpFileOpen(pWriteInfo->fileDest, (sysFILEMODE)(pwpWrite | pwpAscii));
                break;
            case PWP_ENCODING_BINARY:
                pRti->fp = pwpFileOpen(pWriteInfo->fileDest, (sysFILEMODE)(pwpWrite | pwpBinary));
                break;
            case PWP_ENCODING_UNFORMATTED:
                pRti->fp = pwpFileOpen(pWriteInfo->fileDest, (sysFILEMODE)(pwpWrite | pwpUnformatted));
                break;
        }
        ret = (0 != pRti->fp);
    }
    return ret;
}

/*------------------------------------*/
static int
fileDestCwdPush(const char *fileDest)
{
    int ret = 0;
    if (fileDest && fileDest[0]) {
        char dir[FILENAME_MAX];
        char *p;
        strcpy(dir, fileDest);
        p = strrchr(dir, '/');
        if (p) {
            // strip file/base name from end of path
            *p = '\0';
            if (0 == pwpCwdPush(dir)) {
                ret = 1;
            }
        }
    }
    return ret;
}

/**************************************/
int
caeuFileOpen(CAEP_RTITEM* pRti, const CAEP_WRITEINFO *pWriteInfo)
{
    int ret = 0;
    if (pRti) {
        switch (pRti->FormatInfo.fileDest) {
            /* Plugin wants a full "/path/to/filename.ext"
            */
            case PWP_FILEDEST_FILENAME:
                if (fileDestCwdPush(pWriteInfo->fileDest) &&
                    !openFileName(pRti, pWriteInfo)) {
                    pwpCwdPop();
                }
                else {
                    ret = 1;
                }
                break;

            /* Plugin wants a base "/path/to/basefilename" (no ext)
            */
            case PWP_FILEDEST_BASENAME:
                /* Set the cwd to the given folder.
                   It is the plugin's job to open the appropriate files.
                */
                ret = fileDestCwdPush(pWriteInfo->fileDest);
                break;

            /* Plugin wants a folder "/path/to/folder/" (no file or ext)
            */
            case PWP_FILEDEST_FOLDER:
                /* Set the cwd to the given folder.
                   It is the plugin's job to open the appropriate files.
                */
                if (0 == pwpCwdPush(pWriteInfo->fileDest)) {
                    ret = 1;
                }
                break;
        }
    }
    return ret;
}

static const char *invalid = "!invalid";

/**************************************/
const char *
caeuEncodeToText(CAEP_ENUM_ENCODING enc)
{
    switch (enc) {
    case PWP_ENCODING_ASCII: return "ascii";
    case PWP_ENCODING_BINARY: return "binary";
    case PWP_ENCODING_UNFORMATTED: return "unformatted";
    }
    return invalid;
}

/**************************************/
const char *
caeuPrecisionToText(CAEP_ENUM_PRECISION prec)
{
    switch (prec) {
    case PWP_PRECISION_SINGLE: return "single";
    case PWP_PRECISION_DOUBLE: return "double";
    }
    return invalid;
}

/**************************************/
const char *
caeuDimensionToText(CAEP_ENUM_DIMENSION dim)
{
    switch (dim) {
    case PWP_DIMENSION_2D: return "2D";
    case PWP_DIMENSION_3D: return "3D";
    }
    return invalid;
}
