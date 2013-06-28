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

/*! \cond */
/*.................................................
    initialize caepRtItem[0]
*/
#define ID_CaeUnsSU2  30
{
    /*== CAEP_FORMATINFO FormatInfo */
    {   PWP_SITE_GROUPNAME,     /* const char *group */
        "SU2",             /* const char *name */
        MAKEGUID(ID_CaeUnsSU2),  /* PWP_UINT32 id */

        PWP_FILEDEST_FILENAME,  /* CAEP_ENUM_FILEDEST fileDest */

        PWP_FALSE,              /* PWP_BOOL allowedExportConditionsOnly */
        PWP_FALSE,              /* PWP_BOOL allowedVolumeConditions */

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
    0, //CaeUnsSU2BCInfo,            /* CAEP_BCINFO* */
    0, //ARRAYSIZE(CaeUnsSU2BCInfo), /* PWP_UINT32 BCCnt */

    /*== CAEP_VCINFO*    pVCInfo;    -- array of format VCs or NULL */
    /*== PWP_UINT32      VCCnt;      -- # format VCs */
    0, //CaeUnsSU2VCInfo,            /* CAEP_VCINFO* pVCInfo */
    0, //ARRAYSIZE(CaeUnsSU2VCInfo), /* PWP_UINT32 VCCnt */

    /*== const char**    pFileExt;   -- array of valid file extensions */
    /*== PWP_UINT32      ExtCnt;     -- # valid file extensions */
    CaeUnsSU2FileExt,            /* const char **pFileExt */
    ARRAYSIZE(CaeUnsSU2FileExt), /* PWP_UINT32 ExtCnt */

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

    0,  /* PWP_UINT32 progTotal */
    0,  /* PWP_UINT32 progComplete */
    {0},  /* clock_t clocks[CAEPU_CLKS_SIZE]; */
    0,  /* PWP_BOOL opAborted */

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
/*! \endcond */

/************************************************************************/
/*! \file
\brief Static Initialization Data for the CAEP_RTITEM Array

The file \sf{%rtCaepInitItems.h} defines the static, compile-time initialization
of the global CAEP_RTITEM \ref caepRtItem[] array. The CAE Plugin SDK uses this
array to implement the functions and behaviors required by the
\ref DOXGRP_APICAEP. If you want to see the SDK implementation details,
look in the \sf{/shared/CAEP/apiCAEP.c} file.

The SDK file \sf{/shared/CAEP/apiCAEP.c} includes \sf{%rtCaepInitItems.h} as shown
below.
\par
\dontinclude apiCAEP.c
\skip caepRtItem[] =
\until };

The format of \sf{%rtCaepInitItems.h} must be valid for the static
initialization of an array of C-struct's. It is important to note that some of
\ref CAEP_RTITEM's data members are also structs. This will require
curly-braces \p {} around these nested data members. If you are not familiar
with static initialization, see the \ref example_cstruct_init page.

When copied from the \sf{src/plugins/templates/CAEP/} folder to your plugins project
folder, \sf{%rtCaepInitItems.h} will contain example initilization data for 3
CAEP_RTITEM array items. This example data must be culled and edited to
define the settings appropriate for your plugin's implementation.

\note
The global \ref caepRtItem[] is an array so that a plugin can implement
multiple CAE exporters in a single binary. However, due to a limitation in
the "Export-CAE/1.0" implementation, only one CAE exporter is allowed at
this time. As a result, the \ref caepRtItem[] array must be of size 1 for
"Export-CAE/1.0" plugins. Pointwise hopes to remove this limitation in
future releases.

\note
If you add custom data members to CAEP_RTITEM using rtCaepInstanceData.h, be
sure to add the additional static initializers when editing \sf{%rtCaepInitItems.h}
to prevent compiler warnings or errors!
*/

#endif /* _RTCAEPINITITEMS_H_ */
