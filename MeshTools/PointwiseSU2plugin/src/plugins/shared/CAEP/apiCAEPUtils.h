/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#ifndef _APICAEPUTILS_H_
#define _APICAEPUTILS_H_

#include <stdio.h>
#include <time.h>

#include "apiCAEP.h"
#include "apiPWP.h"
#include "apiPWPUtils.h"
#include "apiGridModel.h"

#ifdef __cplusplus
extern "C" {
#endif

/************************************************************************/
/*! \file
    \brief CAEP utilities.

    A collection of helpful data and functions useful to CAEP
    compliant plugins.

    \par "Progress functions"

    Functions used to maintain progress status during an export.
    
    These functions access data stored in a CAEP_RTITEM
    instance to track the progress and forward calls the appropriate
    PwuProgressBegin(), PwuProgressEnd(), PwuProgressStatus(),
    PwuProgressNextStep(), PwuProgressQuit() functions.

    \par "File functions"

    Functions used to open and close a CAE export file using
    CAEP_WRITEINFO runtime settings.

    \par "Search functions"

    Functions used to find an entry in the CAEP_RTITEM array.
*/

// include impl-defined cae instance runtime data types.
// Must also define the macro CAEP_RUNTIME_INSTDATADECL
#include "rtCaepInstanceData.h"

/*! \def CAEP_RUNTIME_INSTDATADECL
* \brief Implementation defined CAE runtime instance data macro.
*
* By default, this macro resolves to nothing (no data).
*
* If a plugin needs additional runtime instance data,
* the developer may add data members to the CAEP_RTITEM structure by
* editing the rtCaepInstanceData.h header.
*/
#if !defined(CAEP_RUNTIME_INSTDATADECL)
#   define CAEP_RUNTIME_INSTDATADECL
#endif

/*---------------------------------------------------------*/
/*! \brief Supported CAEPU clock id values. 
*
* A collection of timer values that can be used for measuring export timings.
* Typically, these values are used for debugging and performance purposes.
*
* All times are initialized to 0.
*
* \sa CAEPU_IS_CLKS_ID, CAEPU_RT_CLKS_ID, CAEPU_RT_CLKS_DIFF,
*   CAEPU_RT_CLKS_POLL, CAEPU_RT_CLKS_POLL_STEP, CAEPU_RT_CLKS_DIFF_STEP,
*   CAEPU_RT_CLKS_POLL_TOTAL, CAEPU_RT_CLKS_DIFF_TOTAL, CAEPU_CLKS_TO_MSECS,
*   CAEPU_CLKS_TO_SECS, CAEPU_CLKS_TO_MINS, CAEPU_CLKS_TO_HOURS,
*   CAEPU_CLKS_TO_FSECS, CAEP_RTITEM::clocks
*/
typedef enum CAEPU_ENUM_CLOCKS_e {
    CAEPU_CLKS_PROGUPDATE, /*!< last time a progress update msg was sent */
    CAEPU_CLKS_PROGINIT,   /*!< time caeuProgressInit() was called */
    CAEPU_CLKS_BEGSTEP,    /*!< time caeuProgressBeginStep() was called */
    CAEPU_CLKS_PROGINCR,   /*!< time caeuProgressIncr() was called */
    CAEPU_CLKS_ENDSTEP,    /*!< time caeuProgressEndStep() was called */
    CAEPU_CLKS_PROGEND,    /*!< time caeuProgressEnd() was called */
    /* add new group enum values before this line */

    /*! \cond sdkINTERNALS */
    /*! for Bookkeeping only - NOT valid "id" values
    @{ */
    CAEPU_CLKS_SIZE,
    CAEPU_CLKS_LAST = CAEPU_CLKS_SIZE-1
    /* @} */
    /*! \endcond sdkINTERNALS */
}
CAEPU_ENUM_CLOCKS;


/************************************************************************/
/*! \brief The data representing a CAE exporter instance.
*/
typedef struct CAEP_RTITEM_t {
    /*! \brief The CAE Plugin format data.
    *
    * Accessed by the low-level API call PwEnumCaeFormat().
    */
    CAEP_FORMATINFO FormatInfo;

    /*! \brief Pointer to the associated PWU_RTITEM structure.
    */
    PWU_RTITEM *pApiData;

    /*! \brief Pointer to an array of supported BC definitions.
    *
    * Accessed by the CAEP-API call PwCaeEnumBCs(). If BCs are not
    * supported, set this to null.
    */
    CAEP_BCINFO *pBCInfo;

    /*! \brief The number of BC definitions.
    *
    * If BCs are not supported, set this to 0.
    */
    PWP_UINT32 BCCnt;

    /*! \brief Pointer to an array of supported VC definitions.
    *
    * Accessed by the CAEP-API call PwCaeEnumVCs(). If VCs are not supported,
    * set this to null.
    */
    CAEP_VCINFO *pVCInfo;

    /*! \brief The number of VC definitions.
    *
    * Accessed by the low-level API call PwCaeGetFileExtCount(). If VCs are
    * not supported, set this to 0.
    */
    PWP_UINT32 VCCnt;

    /*! \brief Pointer to an array of valid file extensions.
    *
    * Accessed by the CAEP-API call PwCaeEnumFileExt(). If not supported, set
    * this to null.
    */
    const char **pFileExt;

    /*! \brief The number of valid file extensions.
    *
    * If not supported, set this to 0.
    */
    PWP_UINT32 ExtCnt;

    /*! \brief Array of supported element-type flags.
    *
    * The exporter supports element type T if elemType[T] is PWP_TRUE. Where T
    * is one of the PWGM_ENUM_ELEMTYPE values.
    */
    PWP_BOOL elemType[PWGM_ELEMTYPE_SIZE];

    /*! \brief Runtime FILE pointer.
    *
    * This data member may be used by the plugin.
    *
    * When a plugin defines FormatInfo.fileDest as PWP_FILEDEST_FILENAME, fp is
    * automatically initialized by the SDK before runtimeWrite() is called and
    * closed by the SDK after runtimeWrite() returns.
    *
    * For all other FormatInfo.fileDest settings, fp will be null. The CAE
    * exporter is allowed to initialize fp using the pwpFileOpen() function as
    * needed.
    *
    * \sa pwpFileOpen(), pwpFileClose(), pwpFileEof(), pwpFileFlush(),
    *     pwpFileGetpos(), pwpFileSetpos(), pwpFileRead(), pwpFileWrite(),
    *     pwpFileWriteStr(), pwpFileRewind(), pwpFileDelete(), pwpCwdPush(),
    *     pwpCwdPop()
    */
    FILE *fp;

    /*! \brief Unformatted file I/O data.
    *
    * Used by the PwuUnfXxxx() helper functions. This data member is for
    * internal use only.
    *
    * \sa PwuUnfFileSetEndianness(), PwuUnfFileGetEndianness(),
    *     PwuUnfFileBegin(), PwuUnfRecBegin(), PwuUnfRecWriteArr(),
    *     PwuUnfRecWriteBuf(), PwuUnfRecWriteEndianBuf()
    */
    PWU_UNFDATA unfData;

    /*! \brief Runtime grid model handle to export.
    *
    * This is the same handle as passed to runtimeWrite().
    */
    PWGM_HGRIDMODEL model;

    /*! \brief Runtime export CAEP_WRITEINFO data.
    *
    * This is the same pointer passed to runtimeWrite().
    */
    const CAEP_WRITEINFO *pWriteInfo;

/*! \cond sdkINTERNALS */

    /*! \brief Internal progress handling data.
    *
    * Used by the caeuProgressXxxx() helper functions. This data member is
    * for internal use only.
    */
    PWP_UINT32 progTotal;

    /*! \brief Internal progress handling data.
    *
    * Used by the caeuProgressXxxx() helper functions. This data member is
    * for internal use only.
    */
    PWP_UINT32 progComplete;

    /*! \brief Internal progress handling data.
    *
    * Used by the caeuProgressXxxx() helper functions. This data member is
    * for internal use only.
    */
    clock_t clocks[CAEPU_CLKS_SIZE];

    /*! \brief Internal progress handling data.
    *
    * Used by the caeuProgressXxxx() helper functions. This data member is
    * for internal use only.
    */
    PWP_BOOL opAborted;

/*! \endcond sdkINTERNALS */

    /*! \brief Implementation defined CAE runtime instance data macro.
    *
    * By default, this macro resolves to nothing (no data).
    *
    * If a plugin needs additional runtime instance data,
    * the developer may add data members to the CAEP_RTITEM structure by
    * editing the rtCaepInstanceData.h header.
    */
    CAEP_RUNTIME_INSTDATADECL
}
CAEP_RTITEM;


/************************************************************************/
/*! \brief The runtime array of CAEP_RTITEM items.
*
* This array is initialized in apiCAEP.c by including the plugin-defined data
* from the rtCaepSupportData.h and rtCaepInitItems.h header files.
*
* There will be one entry in this array for each CAE Exporter supported by
* the plugin.
*
* The items in this array are used to generate CAEP_FORMATINFO values for
* PwEnumCaeFormat() and PwCaeFormat().
*
* \sa PwEnumCaeFormat(), PwCaeFormat(), CAEP_FORMATINFO
*/
extern CAEP_RTITEM caepRtItem[];


/************************************************************************/
/*! \brief The number of entries in caepRtItem[] array.
*
* Accessed by the low-level API call PwGetCaeFormatCount(). This value MUST be
* 1 or greater.
*/
extern PWP_UINT32 caepFormatCnt;


/*************************************************
                 helper functions
**************************************************/

/************************************************************************/
/*! \brief Find an item in caepRtItem[] by it's id.
*
* Find caepRtItem[] item with caepRtItem[n].FormatInfo.id == id.
*
* \param id
*    CAE format GUID
*
* \return Pointer to matching item or null if not found.
*/
CAEP_RTITEM* caeuFindFormatById (PWP_UINT32 id);


/************************************************************************/
/*! \brief Find an item in caepRtItem[] by it's name.
*
* Find caepRtItem[] item with caepRtItem[n].FormatInfo.name == name.
*
* \param name
*    CAE format name
*
* \return Pointer to matching item or null if not found
*/
CAEP_RTITEM* caeuFindFormatByName (const char name[]);


/************************************************************************/
/*! \brief Initializes a progress tracking session.
*
* Called once by a plugin before an export begins.
*
* \param pRti
*    The CAEP_RTITEM pointer passed into runtimeWrite().
*
* \param cnt
*    The number of major steps to be used for this export session.
*
* \par "Sample usage:"
* \code
*   // Set a total of 3 MAJOR progress steps.
*   if (caeuProgressInit(pRti, 3)) {
*       PWP_BOOL ok;
*
*       // begin MAJOR step 1 composed of 44 sub-steps
*       if (ok = caeuProgressBeginStep(pRti, 44)) {
*           int ii;
*           for (ii=0; ii < 44 && ok; ++ii) {
*               // repeat for the 44 sub-steps
*               ok = caeuProgressIncr(pRti);
*           }
*
*           // This call is optional. It is valid to call
*           // caeuProgressBeginStep() or caeuProgressEnd() instead.
*           if (ok) {
*               ok = caeuProgressEndStep(pRti);
*           }
*       }
*
*       // begin MAJOR step 2 composed of 22 sub-steps
*       if (ok && ok = caeuProgressBeginStep(pRti, 22)) {
*           ...snip...
*       }
*
*       // begin MAJOR step 3 composed of 10 sub-steps
*       if (ok && ok = caeuProgressBeginStep(pRti, 10)) {
*           ...snip...
*       }
*       caeuProgressEnd(pRti, ok);
*   }
* \endcode
*/
PWP_BOOL caeuProgressInit(CAEP_RTITEM *pRti, PWP_UINT32 cnt);


/************************************************************************/
/*! \brief Begins a progress tracking step.
*
* Called by a plugin before begining a major export step.
*
* \param pRti
*    The CAEP_RTITEM pointer passed into runtimeWrite().
*
* \param total
*    The number of sub-steps in this major step.
*
* \note See example usage in caeuProgressInit().
*/
PWP_BOOL caeuProgressBeginStep(CAEP_RTITEM *pRti, PWP_UINT32 total);


/************************************************************************/
/*! \brief Completes a progress tracking sub-step.
*
* Called by a plugin after each export sub-step is completed.
*
* \param pRti
*    The CAEP_RTITEM pointer passed into runtimeWrite().
*
* \note See example usage in caeuProgressInit().
*/
PWP_BOOL caeuProgressIncr(CAEP_RTITEM *pRti);


/************************************************************************/
/*! \brief Completes a progress tracking major step.
*
* Called by a plugin after each export major step is completed.
*
* \param pRti
*    The CAEP_RTITEM pointer passed into runtimeWrite().
*
* \note See example usage in caeuProgressInit().
*/
PWP_BOOL caeuProgressEndStep(CAEP_RTITEM *pRti);


/************************************************************************/
/*! \brief Ends all progress tracking.
*
* Called by a plugin after all major steps are completed or when exiting early
* for an export error.
*
* \param pRti
*    The CAEP_RTITEM pointer passed into runtimeWrite().
*
* \param ok
*    Set to PWP_TRUE/PWP_FALSE to indicate sucess/failure.
*
* \note See example usage in caeuProgressInit().
*/
void caeuProgressEnd(CAEP_RTITEM *pRti, PWP_BOOL ok);


/************************************************************************/
/*! \brief Closes pRti for file I/O as specified by pWriteInfo.
*
* Prior to this call, pRti must have been opened with a prior call to
* caeuFileOpen(). The current working directory is restored to the location
* in effect when caeuFileOpen() was called.
*
* \param pRti
*    The CAEP_RTITEM pointer passed into runtimeWrite().
*
* \param pWriteInfo
*    The CAEP_WRITEINFO pointer passed into runtimeWrite().
*
* \return non-zero value on success.
*/
int caeuFileClose(CAEP_RTITEM* pRti, const CAEP_WRITEINFO *pWriteInfo);

/************************************************************************/
/*! \brief Prepare pRti for file I/O as specified by pWriteInfo.
*
* \param pRti
*    The CAEP_RTITEM pointer passed into runtimeWrite().
*
* \param pWriteInfo
*    The CAEP_WRITEINFO pointer passed into runtimeWrite().
*
* \return non-zero value on success.
* \note 
*   The behavior of this call is determined by the value of
*   pRti->FormatInfo.fileDest:
*
*  \li PWP_FILEDEST_FILENAME - pRti->fp is opened to file
*      specified by pWriteInfo->fileDest.
*  \li PWP_FILEDEST_BASENAME - pRti->fp is null.
*  \li PWP_FILEDEST_FOLDER - pRti->fp is null.
*
*  In all cases, the current working directory is set to the path specified in
*  pWriteInfo->fileDest.
*/
int caeuFileOpen(CAEP_RTITEM* pRti, const CAEP_WRITEINFO *pWriteInfo);

/************************************************************************/
/*! \brief Converts a CAEP_ENUM_ENCODING value to a text string
*   representation.
*
* \param enc
*    The CAEP_ENUM_ENCODING value.
*
* \return The text string or "!invalid".
*/
const char * caeuEncodeToText(CAEP_ENUM_ENCODING enc);

/************************************************************************/
/*! \brief Converts a CAEP_ENUM_PRECISION value to a text string
*   representation.
*
* \param prec
*    The CAEP_ENUM_PRECISION value.
*
* \return The text string or "!invalid".
*/
const char * caeuPrecisionToText(CAEP_ENUM_PRECISION prec);

/************************************************************************/
/*! \brief Converts a CAEP_ENUM_DIMENSION value to a text string
*   representation.
*
* \param dim
*    The CAEP_ENUM_DIMENSION value.
*
* \return The text string or "!invalid".
*/
const char * caeuDimensionToText(CAEP_ENUM_DIMENSION dim);



/************************************************************************/
/*! \brief Get the export dimension from CAEP_RTITEM data pointed to by rti.
*/
#define CAEPU_RT_DIM(rti)       (rti)->pWriteInfo->dimension

/************************************************************************/
/*! \brief Get the export dimension text from CAEP_RTITEM data pointed to by
*   rti.
*/
#define CAEPU_RT_DIM_TEXT(rti) \
                caeuDimensionToText((rti)->pWriteInfo->dimension)

/************************************************************************/
/*! \brief Returns PWP_TRUE if export dimension is d.
*/
#define CAEPU_RT_DIM_IS(rti, d) PWP_CAST_BOOL((d) == CAEPU_RT_DIM(rti))

/************************************************************************/
/*! \brief Returns PWP_TRUE if export dimension is PWP_DIMENSION_3D.
*/
#define CAEPU_RT_DIM_3D(rti)     CAEPU_RT_DIM_IS(rti, PWP_DIMENSION_3D)

/************************************************************************/
/*! \brief Returns PWP_TRUE if export dimension is PWP_DIMENSION_2D.
*/
#define CAEPU_RT_DIM_2D(rti)     CAEPU_RT_DIM_IS(rti, PWP_DIMENSION_2D)


/************************************************************************/
/*! \brief Get the export precision from CAEP_RTITEM data pointed to by rti.
*/
#define CAEPU_RT_PREC(rti)       (rti)->pWriteInfo->precision

/************************************************************************/
/*! \brief Get the export precision text from CAEP_RTITEM data pointed to by
*   rti.
*/
#define CAEPU_RT_PREC_TEXT(rti) \
                caeuPrecisionToText((rti)->pWriteInfo->precision)

/************************************************************************/
/*! \brief Returns PWP_TRUE if export precision is prec.
*/
#define CAEPU_RT_PREC_IS(rti, prec) PWP_CAST_BOOL((prec) == CAEPU_RT_PREC(rti))

/************************************************************************/
/*! \brief Returns PWP_TRUE if export precision is PWP_PRECISION_SINGLE.
*/
#define CAEPU_RT_PREC_SINGLE(rti)     CAEPU_RT_PREC_IS(rti, PWP_PRECISION_SINGLE)

/************************************************************************/
/*! \brief Returns PWP_TRUE if export dimension is PWP_PRECISION_DOUBLE.
*/
#define CAEPU_RT_PREC_DOUBLE(rti)     CAEPU_RT_PREC_IS(rti, PWP_PRECISION_DOUBLE)


/************************************************************************/
/*! \brief Get the export encoding from CAEP_RTITEM data pointed to by rti.
*/
#define CAEPU_RT_ENCODING(rti)       (rti)->pWriteInfo->encoding

/************************************************************************/
/*! \brief Get the export precision text from CAEP_RTITEM data pointed to by
*   rti.
*/
#define CAEPU_RT_ENCODING_TEXT(rti) \
                caeuEncodeToText((rti)->pWriteInfo->encoding)

/************************************************************************/
/*! \brief Returns PWP_TRUE if export encoding is enc.
*/
#define CAEPU_RT_ENC_IS(rti, e) PWP_CAST_BOOL((e) == CAEPU_RT_ENCODING(rti))

/************************************************************************/
/*! \brief Returns PWP_TRUE if export encoding is PWP_ENCODING_ASCII.
*/
#define CAEPU_RT_ENC_ASCII(rti)     CAEPU_RT_ENC_IS(rti, PWP_ENCODING_ASCII)

/************************************************************************/
/*! \brief Returns PWP_TRUE if export encoding is PWP_ENCODING_BINARY.
*/
#define CAEPU_RT_ENC_BINARY(rti)    CAEPU_RT_ENC_IS(rti, PWP_ENCODING_BINARY)

/************************************************************************/
/*! \brief Returns PWP_TRUE if export encoding is PWP_ENCODING_UNFORMATTED.
*/
#define CAEPU_RT_ENC_UNFORMATTED(rti)  \
                CAEPU_RT_ENC_IS(rti, PWP_ENCODING_UNFORMATTED)


/************************************************************************/
/*! \brief Returns PWP_TRUE if the export has been aborted by the user or
        the plugin itself.
*/
#define CAEPU_RT_IS_ABORTED(rti)    PWP_CAST_BOOL((rti)->opAborted)

/************************************************************************/
/*! \brief Used by plugin to mark the export operation as aborted.
*/
#define CAEPU_RT_ABORT(rti)    ((rti)->opAborted = PWP_TRUE)

/************************************************************************/
/*! \brief Attempt theOp and abort the export if needed. Returns PWP_FALSE if
*       the export operation is marked as aborted.
*
*   The export operation is marked as aborted if and only if theOp returns 0.
*
*   Though theOp can be any "boolean" expression, it is typically a function
*   call.
*
* \note If theOp is a function call, it will only be called if the export has
*       not already been marked as aborted.
*
* \par "Sample usage:"
* \code
*   // assume:
*   //      int someFunc(int n);
*   // returns 0 on failure
*
*   caeuProgressInit(pRti, 3);
*
*   caeuProgressBeginStep(pRti, 44); // 44 substeps
*
*   if (CAEPU_RT_TRY(pRti, someFunc(1))) {
*       // someFunc(1) succeeded (returned a non-zero value)
*       ...snip...
*   }
*
*   caeuProgressBeginStep(pRti, 55); // 55 substeps
*
*   // The call to someFunc(2) will NOT be made if someFunc(1) failed or if
*   // CAEPU_RT_IS_ABORTED(pRti) is already true for some other reason.
*
*   if (CAEPU_RT_TRY(pRti, someFunc(2))) {
*       // someFunc(2) succeeded (returned a non-zero value)
*       ...snip...
*   }
*
*   caeuProgressEnd(pRti, !CAEPU_RT_IS_ABORTED(pRti));
* \endcode
*/
#define CAEPU_RT_TRY(rti, theOp) \
            ((rti)->opAborted ? PWP_FALSE : \
                (((rti)->opAborted = !(theOp)), !CAEPU_RT_IS_ABORTED(rti)))


/************************************************************************/
/*! \brief Returns PWP_TRUE if id is a valid CAEPU_ENUM_CLOCKS id.
*/
#define CAEPU_IS_CLKS_ID(id) \
            PWP_CAST_BOOL(((id) >= 0) && ((id) <= CAEPU_CLKS_LAST))

/************************************************************************/
/*! \brief Returns the clock time value for the given id. The id is validated.
*/
#define CAEPU_RT_CLKS_ID(rti, id) (CAEPU_IS_CLKS_ID(id)? (rti)->clocks[id]: 0)

/************************************************************************/
/*! \brief Returns the clock time difference between startId and endId as
*       clocks[endId] - clocks[startId]. The ids are validated.
*/
#define CAEPU_RT_CLKS_DIFF(rti, startId, endId) \
    ((CAEPU_IS_CLKS_ID(startId) && CAEPU_IS_CLKS_ID(endId)) ? \
        ((rti)->clocks[endId] - (rti)->clocks[startId]) : 0)

/************************************************************************/
/*! \brief Returns the current clock time difference between 
*       id and CAEPU_CLKS_PROGINCR.
*/
#define CAEPU_RT_CLKS_POLL(rti, id) \
            CAEPU_RT_CLKS_DIFF(rti, id, CAEPU_CLKS_PROGINCR)

/************************************************************************/
/*! \brief Returns the current clock time difference between 
*       CAEPU_CLKS_BEGSTEP and CAEPU_CLKS_PROGINCR.
*/
#define CAEPU_RT_CLKS_POLL_STEP(rti) \
            CAEPU_RT_CLKS_POLL(rti, CAEPU_CLKS_BEGSTEP)

/************************************************************************/
/*! \brief Returns the current clock time difference between 
*       CAEPU_CLKS_BEGSTEP and CAEPU_CLKS_ENDSTEP.
*/
#define CAEPU_RT_CLKS_DIFF_STEP(rti) \
            CAEPU_RT_CLKS_DIFF(rti, CAEPU_CLKS_BEGSTEP, CAEPU_CLKS_ENDSTEP)

/************************************************************************/
/*! \brief Returns the current clock time difference between 
*       CAEPU_CLKS_PROGINIT and CAEPU_CLKS_PROGINCR.
*/
#define CAEPU_RT_CLKS_POLL_TOTAL(rti) \
            CAEPU_RT_CLKS_POLL(rti, CAEPU_CLKS_PROGINIT)

/************************************************************************/
/*! \brief Returns the current clock time difference between 
*       CAEPU_CLKS_PROGINIT and CAEPU_CLKS_PROGEND.
*/
#define CAEPU_RT_CLKS_DIFF_TOTAL(rti) \
            CAEPU_RT_CLKS_POLL(rti, CAEPU_CLKS_PROGINIT, CAEPU_CLKS_PROGEND)

/************************************************************************/
/*! \brief Returns the clock value c as milli seconds (PWP_INT32). Only whole
*       ms values are possible.
*/
#define CAEPU_CLKS_TO_MSECS(c) ((PWP_INT32)((c) * 1000) / CLOCKS_PER_SEC)

/************************************************************************/
/*! \brief Returns the clock value c as seconds (PWP_INT32). Only whole second
*       values are possible.
*/
#define CAEPU_CLKS_TO_SECS(c) ((PWP_INT32)(c) / CLOCKS_PER_SEC)

/************************************************************************/
/*! \brief Returns the clock value c as minutes (PWP_INT32). Only whole minute
*       values are possible.
*/
#define CAEPU_CLKS_TO_MINS(c) (CAEPU_CLKS_TO_SECS(c) / 60)

/************************************************************************/
/*! \brief Returns the clock value c as hours (PWP_INT32). Only whole hour
*       values are possible.
*/
#define CAEPU_CLKS_TO_HOURS(c) (CAEPU_CLKS_TO_MINS(c) / 60)

/************************************************************************/
/*! \brief Returns the clock value c as seconds (PWP_FLOAT).
*/
#define CAEPU_CLKS_TO_FSECS(c) ((PWP_FLOAT)(c) / (PWP_FLOAT)CLOCKS_PER_SEC)


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  /* _APICAEPUTILS_H_ */
