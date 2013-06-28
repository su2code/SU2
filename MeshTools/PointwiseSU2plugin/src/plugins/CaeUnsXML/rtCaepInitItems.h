/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#ifndef _RTCAEPINITITEMS_H_
#define _RTCAEPINITITEMS_H_

/*.................................................
    initialize caepRtItem[0]
*/
#define ID_CaeUnsXML  10
{
    /*== CAEP_FORMATINFO FormatInfo */
    {   PWP_SITE_GROUPNAME,   /* const char *group */
        "XML-Unstructured",      /* const char *name */
        MAKEGUID(ID_CaeUnsXML),  /* PWP_UINT32 id */

        PWP_FILEDEST_FILENAME,  /* CAEP_ENUM_FILEDEST fileDest */

        PWP_TRUE,               /* PWP_BOOL allowedExportConditionsOnly */
        PWP_TRUE,               /* PWP_BOOL allowedVolumeConditions */

        PWP_TRUE,               /* PWP_BOOL allowedFileFormatASCII */
        PWP_FALSE,              /* PWP_BOOL allowedFileFormatBinary */
        PWP_FALSE,              /* PWP_BOOL allowedFileFormatUnformatted */

        PWP_TRUE,               /* PWP_BOOL allowedDataPrecisionSingle */
        PWP_TRUE,               /* PWP_BOOL allowedDataPrecisionDouble */

        PWP_TRUE,               /* PWP_BOOL allowedDimension2D */
        PWP_TRUE                /* PWP_BOOL allowedDimension3D */
    },

    &pwpRtItem[1],  /* PWU_RTITEM* */

    /*== CAEP_BCINFO*    pBCInfo;    -- array of format BCs or NULL */
    /*== PWP_UINT32      BCCnt;      -- # format BCs */
    CaeUnsXMLBCInfo,            /* CAEP_BCINFO* */
    ARRAYSIZE(CaeUnsXMLBCInfo), /* PWP_UINT32 BCCnt */

    /*== CAEP_VCINFO*    pVCInfo;    -- array of format VCs or NULL */
    /*== PWP_UINT32      VCCnt;      -- # format VCs */
    CaeUnsXMLVCInfo,            /* CAEP_VCINFO* pVCInfo */
    ARRAYSIZE(CaeUnsXMLVCInfo), /* PWP_UINT32 VCCnt */

    /*== const char**    pFileExt;   -- array of valid file extensions */
    /*== PWP_UINT32      ExtCnt;     -- # valid file extensions */
    CaeUnsXMLFileExt,            /* const char **pFileExt */
    ARRAYSIZE(CaeUnsXMLFileExt), /* PWP_UINT32 ExtCnt */

    /*== PWP_BOOL  elemType[PWGM_ELEMTYPE_SIZE]; -- un/supported elem */
    {   PWP_TRUE,              /* elemType[PWGM_ELEMTYPE_BAR] */
        PWP_TRUE,              /* elemType[PWGM_ELEMTYPE_HEX] */
        PWP_TRUE,              /* elemType[PWGM_ELEMTYPE_QUAD] */
        PWP_TRUE,              /* elemType[PWGM_ELEMTYPE_TRI] */
        PWP_TRUE,              /* elemType[PWGM_ELEMTYPE_TET] */
        PWP_TRUE,              /* elemType[PWGM_ELEMTYPE_WEDGE] */
        PWP_TRUE },            /* elemType[PWGM_ELEMTYPE_PYRAMID] */

    0,  /* FILE *fp */

    /* PWU_UNFDATA UnfData */
    {   0,          /* PWP_UINT32 status */
        0,          /* FILE *fp */
        0,          /* sysFILEPOS fPos */
        PWP_FALSE,  /* PWP_BOOL hadError */
        PWP_FALSE,  /* PWP_BOOL inRec */
        0,          /* PWP_UINT32 recBytes */
        0,          /* PWP_UINT32 totRecBytes */
        0    },     /* PWP_UINT32 recCnt */

    0,  /* PWGM_HGRIDMODEL model */

    0,  /* const CAEP_WRITEINFO *pWriteInfo */

    0,   /* PWP_UINT32 progTotal */
    0,   /* PWP_UINT32 progComplete */
    {0}, /* clock_t clocks[] */
    0,   /* PWP_BOOL opAborted */

    /* if you added any custom data in rtCaepInstanceData.h,
       you need to initialize it here. The init below matches the 
       example MY_CAEP_DATA struct given in rtCaepInstanceData.h */
    /*
    {   0,
        0,
        0.0,
        "string" },
    */
},

#endif /* _RTCAEPINITITEMS_H_ */
