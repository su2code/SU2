/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#ifndef _APIPWPUTILS_H_
#define _APIPWPUTILS_H_

#include "apiPWP.h"
#include "pwpPlatform.h"


#ifdef __cplusplus
extern "C" {
#endif

/*! \file
    \brief Data and functions useful to PWP-API compliant plugins.
*/

/***********************************************************/
/***********************************************************/
/*! \addtogroup DOXGRP_APIPWP
*  @{
*/

/*! \defgroup DOXGRP_APIPWP_UTILS PWP-API Utilities
*
*  A collection of SDK helper/utility data and functions that implement the
*  PWP-API.
*
*  @{
*/

//***************************************************************************
//***************************************************************************
/*! \defgroup DOXGRP_APIPWP_UTILS_DATA PWP-API/SDK Data Types and Instances
*
*  The SDK data types and instances that implement the PWP-API.
*
*  @{
*/
/*! \brief The runtime data representing a PWP-API suported by a plugin.
*/
typedef struct PWP_RTITEM_t {
    /*! \brief The PWP-API instance information.
    *
    * This data is returned by PwpEnumAPIs(). Edit the
    * .../YourPlugin/rtPwpInitItems.h header file to initialize this data.
    */
    PWP_APIINFO apiInfo;

    /*! \brief The API's assigned message callback.
    *
    * The callback is PWP_NULL by default. The callback is set by the
    * framework using PwpSetMessageCallback().
    */
    PWP_MESSAGECB msgCB;
}
PWU_RTITEM;


/*! \brief The runtime array of PWU_RTITEM items.
*
* This array is initialized in apiPWP.c by including the plugin-defined data
* from the .../YourPlugin/rtPwpVersions.h and .../YourPlugin/rtPwpInitItems.h
* header files.
*
* There will be one entry in this array for each published API supported by
* the plugin and one entry for each "special" unpublished APIs used internally
* by the SDK. The unpublished APIs are always last in the array.
*
* The published items in this array are used to generate return values for
* PwpEnumAPIs().
*
* \note A typical CAE plugin will define 2 items:
*  - The base PWP-API, "Plugin-PWP/1.0"
*  - The CAEP-API, "Export-CAE/1.0"
*
* \sa PwpEnumAPIs(), PWP_MESSAGECB_DEFAULT, PWP_MESSAGECB_SPY
*/
extern PWU_RTITEM pwpRtItem[];


/*! \brief The total # of published and unpublished entries in pwpRtItem[]
*
* \sa PwpEnumAPIs(), PWP_MESSAGECB_DEFAULT, PWP_MESSAGECB_SPY
*/
extern PWP_UINT32 totalApiCnt;


/*! \brief The total # of published entries in pwpRtItem[]
*
* \sa PwpEnumAPIs()
*/
extern PWP_UINT32 publishedApiCnt;

/*! @} */ /* DOXGRP_APIPWP_UTILS_DATA */


//***************************************************************************
//***************************************************************************
/*! \defgroup DOXGRP_APIPWP_UTILS_SDKPWP PWP/SDK Data Access Functions
*
*  These calls are used by the SDK to access the SDK data structures. They
*  are avaliable for use by a plugin, but probably will not be needed.
*
*  @{
*/
//-------------------------------------------------------------------------
/*! \brief Find any api in pwpRtItem[].
* \param api The API to search for.
* \return Pointer to the found PWU_RTITEM or PWP_NULL.
* \sa PWP_MESSAGECB_DEFAULT, PWP_MESSAGECB_SPY
* \note The PWP SDK maintains unpublished API's to support the "special"
*       message callbacks.
*/
PWU_RTITEM* PwuFindTotalAPI (const char api[]);

//-------------------------------------------------------------------------
/*! \brief Find a published api in pwpRtItem[].
*
* \param api The API to search for.
* \return Pointer to the found PWU_RTITEM or PWP_NULL.
*/
PWU_RTITEM* PwuFindPublishedAPI (const char api[]);

//-------------------------------------------------------------------------
/*! \brief Search pwpRtItem[] for an API's messageCB.
*
* \param api The API to search for.
* \return The message callback function pointer or PWP_NULL if not found or
*         callback is not set.
*/
PWP_MESSAGECB PwuFindApiMsgCB (const char api[]);

/*! @} */ /* DOXGRP_APIPWP_UTILS_SDKPWP */


//***************************************************************************
//***************************************************************************
/*! \defgroup DOXGRP_APIPWP_UTILS_MSGS Message Handling
*
*  Bundle and send a text message back to framework.
*
*  The SDK implements the PwpSetMessageCallback() and PwpGetMessageCallback()
*  API functions. If any framework message callbacks are registered,
*  these calls will bundle and route the text messages as required by the
*  PWP-API.
*
*  \sa PWP_MESSAGECB, PWP_MSG_TEXT
*
*  @{
*/
//-------------------------------------------------------------------------
/*! \brief Send a message from an api.
*
*  Generic form of the function. Other, more specific forms are provided for
*  common messages.
* \param api The API sending the message. One of the values returned from
*            PwpEnumAPIs().
* \param id The message type id.
* \param pMsg The message type specific data.
*/
PWP_UINT32 PwuSendMsg (const char api[], PWP_ENUM_MSGID id, void *pMsg);

//-------------------------------------------------------------------------
/*! \brief Send a debug text message (PWP_MSGID_DEBUG) to the framework.
*
* \param api The API sending the message. One of the values returned from
*            PwpEnumAPIs().
* \param txt The message text.
* \param code The API-defined message code.
*/
void PwuSendDebugMsg (const char api[], const char txt[], PWP_UINT32 code);

//-------------------------------------------------------------------------
/*! \brief Send an info text message (PWP_MSGID_INFO) to the framework.
*
* \param api The API sending the message. One of the values returned from
*            PwpEnumAPIs().
* \param txt The message text.
* \param code The API-defined message code.
*/
void PwuSendInfoMsg (const char api[], const char txt[], PWP_UINT32 code);

//-------------------------------------------------------------------------
/*! \brief Send a warning text message (PWP_MSGID_WARNING) to the framework.
*
* \param api The API sending the message. One of the values returned from
*            PwpEnumAPIs().
* \param txt The message text.
* \param code The API-defined message code.
*/
void PwuSendWarningMsg (const char api[], const char txt[], PWP_UINT32 code);

//-------------------------------------------------------------------------
/*! \brief Send an error text message (PWP_MSGID_ERROR) to the framework.
*
* \param api The API sending the message. One of the values returned from
*            PwpEnumAPIs().
* \param txt The message text.
* \param code The API-defined message code.
*/
void PwuSendErrorMsg (const char api[], const char txt[], PWP_UINT32 code);

/*! @} */ /* DOXGRP_APIPWP_UTILS_MSGS */


//***************************************************************************
//***************************************************************************
/*! \defgroup DOXGRP_APIPWP_UTILS_PROGRESS Progress Reporting
*
*  Bundle and send a progress message back to framework.
*
*  The SDK implements the PwpSetMessageCallback() and PwpGetMessageCallback()
*  API functions. If any framework message callbacks are registered,
*  these calls will bundle and route the progress messages as required by the
*  PWP-API.
*
*  \sa PWP_MESSAGECB, PWP_MSG_PROGRESS
*
*  @{
*/
//-------------------------------------------------------------------------
/*! \brief Send a progress begin message (PWP_MSGID_PROGBEGIN) to the
*   framework.
*
* \param api The API sending the message. One of the values returned from
*            PwpEnumAPIs().
* \param totalSteps The total number of major steps.
* \return PWP_TRUE if framework wants to continue. PWP_FALSE if framework has
*         aborted op.
*/
PWP_BOOL PwuProgressBegin (const char api[], PWP_UINT32 totalSteps);

//-------------------------------------------------------------------------
/*! \brief Send a progress end message (PWP_MSGID_PROGEND) to the framework.
*
* \param api The API sending the message. One of the values returned from
*            PwpEnumAPIs().
* \param ok The progress status; Set to PWP_TRUE to continue. Set to PWP_FALSE
*           to indicate operation failure and stop.
* \return PWP_TRUE if framework wants to continue. PWP_FALSE if framework has
*         aborted op.
*/
void PwuProgressEnd (const char api[], PWP_BOOL ok);

//-------------------------------------------------------------------------
/*! \brief Send a progress status message (PWP_MSGID_PROGSTATUS, value >= 0)
*   to the framework.
*
* \param api The API sending the message. One of the values returned from
*            PwpEnumAPIs().
* \param complete The number of completed sub-steps.
* \param total The total number of sub-steps in the current step.
* \return PWP_TRUE if framework wants to continue. PWP_FALSE if framework has
*         aborted op.
*/
PWP_BOOL PwuProgressStatus (const char api[], PWP_UINT32 complete,
                              PWP_UINT32 total);

//-------------------------------------------------------------------------
/*! \brief Send a progress "next step" message (PWP_MSGID_PROGSTATUS,
*   value = -1) to the framework.
*
* \param api The API sending the message. One of the values returned from
*            PwpEnumAPIs().
* \return PWP_TRUE if framework wants to continue. PWP_FALSE if framework has
*         aborted op.
*/
PWP_BOOL PwuProgressNextStep (const char api[]);

//-------------------------------------------------------------------------
/*! \brief Send a progress query-quit message (PWP_MSGID_PROGQUIT) to the
*   framework.
*
* \param api The API sending the message. One of the values returned from
*            PwpEnumAPIs().
* \return PWP_TRUE if framework has aborted op. PWP_FALSE if framework wants
*         to continue.
*/
PWP_BOOL PwuProgressQuit (const char api[]);

/*! @} */ /* DOXGRP_APIPWP_UTILS_PROGRESS */


//***************************************************************************
//***************************************************************************
/*! \defgroup DOXGRP_APIPWP_UTILS_ENDIAN Platform Endian Handling
*
*  Platform Endianness detection.
*
*  @{
*/

//--------------------------------------------
/*! \brief Flags used to indicate endianness or control endian behaviors in
*   functions.
*
* \sa PwuUnfFileSetEndianness(), PwuUnfFileGetEndianness()
*/
typedef enum PWU_ENDIANNESS_t {
    PWU_ENDIAN_ERROR,   /*!< error indicator */
    PWU_ENDIAN_LITTLE,  /*!< force little-endian */
    PWU_ENDIAN_BIG,     /*!< force big-endian */
    PWU_ENDIAN_NATIVE,  /*!< use native platform endian */
    PWU_ENDIAN_FOREIGN  /*!< force opposite of native platform endian */
}
PWU_ENDIANNESS;

//--------------------------------------------
/*! \brief Query the OS's native endianness.
*
* \return One of PWU_ENDIAN_LITTLE or PWU_ENDIAN_BIG
*/
PWU_ENDIANNESS PwuGetOsEndianness(void);

/*! @} */ /* DOXGRP_APIPWP_UTILS_ENDIAN */


//***************************************************************************
//***************************************************************************
/*! \defgroup DOXGRP_APIPWP_UTILS_UNFORMATTED Unformatted FORTRAN File Handling
*
*  Unformatted FORTRAN file I/O utility functions and data-types.
*
*  @{
*/

//--------------------------------------------
/*! \brief Unformatted file data block.
*
*   Data block used by the \c PwuUnfFileXxx() functions to track and control
*   unformated file I/O operations.
*
* \note This data is intended for internal use only and should
*       not be accessed by a plugin directly.
*/
typedef struct UNFDATA_t {
    PWP_UINT32     status;      /*!< current file status */
    FILE           *fp;         /*!< file pointer */
    sysFILEPOS     fPos;        /*!< file position value */
    PWP_BOOL       hadError;    /*!< error flag */
    PWP_BOOL       inRec;       /*!< "in record" flag */
    PWP_UINT32     recBytes;    /*!< # bytes written to current record */
    PWP_UINT32     totRecBytes; /*!< total # bytes written to all records */
    PWP_UINT32     recCnt;      /*!< # of records written */
    PWU_ENDIANNESS endianness;  /*!< write data using this endianness */
}
PWU_UNFDATA;

//-------------------------------------------------------------------------
/*! \brief Set the output endianness.
*
*   The endianness will apply to all (future) writes using this PWU_UNFDATA
*   block.
*
* \param pUData Pointer to an I/O control block initialized by
*               PwuUnfFileBegin().
* \param endianness The endianness flag. See notes.
* \return Previous endianness setting.
* \note If endianness may be one of:
*   \li \c PWU_ENDIAN_LITTLE - force output to little endian
*   \li \c PWU_ENDIAN_BIG - force output to big endian
*   \li \c PWU_ENDIAN_NATIVE - force output to OS's native endianness
*   \li \c PWU_ENDIAN_FOREIGN - force output to OS's non-native endianness
*/
PWU_ENDIANNESS PwuUnfFileSetEndianness (PWU_UNFDATA *pUData,
                                        PWU_ENDIANNESS endianness);

//-------------------------------------------------------------------------
/*! \brief Get the output endianness setting for this PWU_UNFDATA block.
*
* \param pUData Pointer to an I/O control block initialized by a call to
*               PwuUnfFileBegin().
* \return Current endianness setting.
*/
PWU_ENDIANNESS PwuUnfFileGetEndianness (PWU_UNFDATA *pUData);

//-------------------------------------------------------------------------
/*! \brief Prepares a PWU_UNFDATA block for a new unformatted file I/O session.
*
* \param fp The file pointer used for writing. The file must be open and
*           ready for write before calling this function. The file position is
*           not changed by this call.
* \param pUData Pointer to an uninitialized I/O control block.
* \return PWP_TRUE on success. PWP_FALSE on failure.
* \sa pwpFileOpen(), PwuUnfFileSetEndianness()
* \note PWU_ENDIAN_NATIVE is the default endianness.
*/
PWP_BOOL PwuUnfFileBegin (FILE *fp, PWU_UNFDATA *pUData);

//-------------------------------------------------------------------------
/*! \brief Prepares a PWU_UNFDATA block for writing a new unformatted data
*   record.
*
* \param pUData Pointer to an I/O control block initialized by a call to
*               PwuUnfFileBegin().
* \return PWP_TRUE on success. PWP_FALSE on failure.
* \sa PwuUnfRecEnd()
* \note If a record is currently active, it will be ended automatically before
*       the new record is prepared.
*/
PWP_BOOL PwuUnfRecBegin (PWU_UNFDATA *pUData);

//-------------------------------------------------------------------------
/*! \brief Write an array of data to the current record.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param arr Pointer to the array data.
* \param itemSize Size of each array item (bytes).
* \param itemCnt The number of array items.
* \return PWP_TRUE on success. PWP_FALSE on failure.
* \note Data is written as-is. Endianness is not enforced.
*/
PWP_BOOL PwuUnfRecWriteArr (PWU_UNFDATA *pUData, const void *arr,
                            size_t itemSize, size_t itemCnt);

//-------------------------------------------------------------------------
/*! \brief Write a data buffer to the current record.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param buf Pointer to the data buffer.
* \param size Size of the data buffer (bytes).
* \return PWP_TRUE on success. PWP_FALSE on failure.
* \note Data is written as-is. Endianness is not enforced.
*/
PWP_BOOL PwuUnfRecWriteBuf (PWU_UNFDATA *pUData, const void *buf, size_t size);

//-------------------------------------------------------------------------
/*! \brief Write a data buffer to the current record with endian order applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param buf Pointer to the data buffer.
* \param size Size of the data buffer (bytes). Size must be an even value
*             between 2 and 32.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteEndianBuf(PWU_UNFDATA *pUData, const void *buf,
                                 size_t size);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_INT value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteINT (PWU_UNFDATA *pUData, PWP_INT val);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_UINT value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteUINT (PWU_UNFDATA *pUData, PWP_UINT val);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_INT8 value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteINT8 (PWU_UNFDATA *pUData, PWP_INT8 val);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_UINT8 value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteUINT8 (PWU_UNFDATA *pUData, PWP_UINT8 val);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_INT16 value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteINT16 (PWU_UNFDATA *pUData, PWP_INT16 val);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_UINT16 value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteUINT16 (PWU_UNFDATA *pUData, PWP_UINT16 val);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_INT32 value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteINT32 (PWU_UNFDATA *pUData, PWP_INT32 val);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_UINT32 value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteUINT32 (PWU_UNFDATA *pUData, PWP_UINT32 val);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_INT64 value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteINT64 (PWU_UNFDATA *pUData, PWP_INT64 val);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_UINT64 value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteUINT64 (PWU_UNFDATA *pUData, PWP_UINT64 val);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_FLOAT value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteFLOAT (PWU_UNFDATA *pUData, PWP_FLOAT val);

//-------------------------------------------------------------------------
/*! \brief Write a PWP_REAL value to the current record with endian order
*   applied.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \param val The value to write.
* \return PWP_TRUE on success. PWP_FALSE on failure.
*/
PWP_BOOL PwuUnfRecWriteREAL (PWU_UNFDATA *pUData, PWP_REAL val);

//-------------------------------------------------------------------------
/*! \brief Finalize the current record write.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfRecBegin().
* \return PWP_TRUE on success. PWP_FALSE on failure.
* \sa PwuUnfRecBegin()
*/
PWP_BOOL PwuUnfRecEnd (PWU_UNFDATA *pUData);

//-------------------------------------------------------------------------
/*! \brief Finalize an unformatted file I/O session.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfFileBegin().
* \return PWP_TRUE on success. PWP_FALSE on failure.
* \sa pwpFileClose()
* \note If a record is currently active, it will be ended automatically before
*       the I/O session is finalized.
* \note It is the callers reponsibility to close the file pointer originally
*       passed to PwuUnfFileBegin().
*/
PWP_BOOL PwuUnfFileEnd (PWU_UNFDATA *pUData);

//-------------------------------------------------------------------------
/*! \brief Check if an unformatted file I/O session has detected any errors.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfFileBegin().
* \return PWP_TRUE if an error was detected. PWP_FALSE if an error was not
*         detected.
* \note It is valid to call this function after PwuUnfFileEnd(). The internals
*       are only cleared on a call to PwuUnfFileBegin().
*/
PWP_BOOL PwuUnfHadError (PWU_UNFDATA *pUData);

//-------------------------------------------------------------------------
/*! \brief Get the running total number of bytes written to the current record
*   during an unformatted file I/O session.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfFileBegin().
* \return The number of bytes written to the current record so far.
* \note It is valid to call this function after PwuUnfRecEnd(). The record
*       internals are only cleared on a call to PwuUnfRecBegin().
*/
PWP_UINT32 PwuUnfRecBytes (PWU_UNFDATA *pUData);

//-------------------------------------------------------------------------
/*! \brief Get the running total number of bytes written to all records during
*   an unformatted file I/O session.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfFileBegin().
* \return The total number of bytes written to all records so far.
* \note It is valid to call this function after PwuUnfFileEnd(). The internals
*       are only cleared on a call to PwuUnfFileBegin().
*/
PWP_UINT32 PwuUnfTotRecBytes (PWU_UNFDATA *pUData);

//-------------------------------------------------------------------------
/*! \brief Get the running total number of finalized records written during an
*   unformatted file I/O session.
*
* \param pUData Pointer to an I/O control block already prepared by a call to
*               PwuUnfFileBegin().
* \return The total number of bytes written to all records so far.
* \note It is valid to call this function after PwuUnfFileEnd(). The internals
*       are only cleared on a call to PwuUnfFileBegin().
*/
PWP_UINT32 PwuUnfRecCount (PWU_UNFDATA *pUData);

/*! @} */ /* DOXGRP_APIPWP_UTILS_UNFORMATTED */

//***************************************************************************
//***************************************************************************
/*! \defgroup DOXGRP_APIPWP_UTILS_HEAP Private Heap Manager Tools
*
*  Implements a private heap manager subsystem.
*
*  Very useful for plugins that may need to allocate a HUGE number of small
*  memory allocations. On some platforms, there is a significant overhead per
*  allocation.
*
*  @{
*/

/***************************************************************************/
/*! \brief Heap statistics data structure passed to PWU_HEAP_CALCSTATS.
*/
typedef struct PwuHeapStats_t {
	size_t totalUnused;   /*!< Total unused items in all segments */
	size_t totalSize;     /*!< Total available items in all segments */
    float totalUnusedPct; /*!< Total unused items as a 0 to 100 percentage */
	size_t segCnt;        /*!< Total number of segments */
	size_t fullSegCnt;    /*!< Number of full segments */
} PwuHeapStats;


/*! \cond sdkINTERNALS */

// fwd declare
typedef struct PwuHeap_t PwuHeap;
typedef struct PwuHeapSeg_t PwuHeapSeg;

/***************************************************************************/
struct PwuHeapSeg_t {
	PwuHeap     *pHeap;    /*!< ptr to the segment's heap */
	PwuHeapSeg  *pNextSeg; /*!< ptr to next segment in heap chain */
	size_t      size;      /*!< # items in segment */
	size_t      offset;    /*!< offset into heap of last used item */
};

/***************************************************************************/
struct PwuHeap_t {
	const char*    name;        /*!< the heap datatype's name */
	size_t         hint;        /*!< suggested # items per heap segment */
	size_t         itemSize;    /*!< bytes per item */
	size_t         segHdrSize;  /*!< bytes per item */
	PwuHeapSeg     *pSeg;       /*!< start of active segment chain */
	PwuHeapSeg     *pFullSeg;   /*!< start of full segment chain */
};

void* PwuHeapAllocitems(PwuHeap *pHeap, size_t num);
void PwuHeapDtor(PwuHeap *pHeap);
void PwuHeapCalcStats(PwuHeap *pHeap, PwuHeapStats *pStats);

/*! \endcond */ /* sdkINTERNALS */


/***************************************************************************/
/*                        PUBLIC INTERFACE TO HEAP                         */
/***************************************************************************/

/***************************************************************************/
/*! \brief The heap-type suitable for function parameters, etc.
*
* \param datatype  The heap item type-name.
*
* \par
* \code
*  // C-function, pass heap by ptr
*  void myCFunc(PWU_HEAP_TYPE(MyStruct) *pHeap) {
*       MyStruct *pChunk = PWU_HEAP_ALLOCITEMS(*pHeap, numItems);
*  }
*  // C++-function, pass heap by ref
*  void myCppFunc(PWU_HEAP_TYPE(MyStruct) &heap) {
*       MyStruct *pChunk = PWU_HEAP_ALLOCITEMS(heap, numItems);
*  }
* \endcode
*/
#define PWU_HEAP_TYPE(datatype) \
	PwuHeap_##datatype


/***************************************************************************/
/*! \brief Declares the strictly-typed, public interface for a given heap-type.
*
* \param datatype  The heap item type-name.
*
* \note This macro must be used in conjunction with a corresponding
*       PWU_HEAP_BODY statement.
*
* \note The framework already declares heap datatype's for PWP_INT8, PWP_UINT8,
*       PWP_INT16, PWP_UINT16, PWP_INT32, PWP_UINT32, PWP_INT64, PWP_UINT64,
*       PWP_FLOAT, and PWP_REAL. You can use PWU_HEAP_CREATE for these item
*       types.
*
* \par
* \code
*  PWU_HEAP_CREATE(int8Heap, PWP_INT8, 100000);
*  PWU_HEAP_CREATE(uint8Heap, PWP_UINT8, 100000);
*  PWU_HEAP_CREATE(int16Heap, PWP_INT16, 100000);
*  PWU_HEAP_CREATE(uint16Heap, PWP_UINT16, 100000);
*  ..snip ..
*  PWU_HEAP_CREATE(floatHeap, PWP_FLOAT, 100000);
*  PWU_HEAP_CREATE(realHeap, PWP_REAL, 100000);
* \endcode
*
* \sa PWU_HEAP_BODY, PWU_HEAP_CREATE, PWU_HEAP_ALLOCITEMS
*/
#define PWU_HEAP_HEADER(datatype) \
	/*---------------------------------------------------------------------*/ \
	typedef struct PwuHeap_##datatype##_t PwuHeap_##datatype; \
	typedef	datatype* allocItemsFunc_##datatype(PwuHeap_##datatype *pHeap, size_t num); \
	struct PwuHeap_##datatype##_t { \
		PwuHeap heap; /* base heap struct */ \
		allocItemsFunc_##datatype *allocItems; \
	}; \
	extern datatype* PwuHeap_##datatype##_allocitems(PwuHeap_##datatype *pHeap, size_t num); \
	/*---------------------------------------------------------------------*/ \
	typedef struct PwuHeap_##datatype##_seg_t PwuHeap_##datatype##_seg; \
	struct PwuHeap_##datatype##_seg_t { \
		PwuHeapSeg seg; /* base seg struct */ \
		datatype   items[1]; /* this segment's item array */ \
	}


/***************************************************************************/
/*! \brief External declaration for a distinct instance of a typed heap.
*
* \param heapName  The heap instance name.
* \param datatype  The heap item type-name.
*
* \note This macro must be used after PWU_HEAP_HEADER and in conjunction with
*       a corresponding PWU_HEAP_CREATE statement.
*
* \note Typically placed in a .h file when a heap is shared among multiple
*       compilation units.
*
* \par
* \code
*  //--------------------------------------------
*  // in myHeap.h
*  //--------------------------------------------
*  // my heap datatype
*  typedef struct MyStruct_t {
*     ..snip..
*  } MyStruct;
*
*  // declare public heap interface
*  PWU_HEAP_HEADER(MyStruct);
*
*  // extern the heap instance
*  PWU_HEAP_EXTERN(myHeap, MyStruct);
*
*  //--------------------------------------------
*  // in myHeap.c
*  //--------------------------------------------
*  // These could also be moved to either file1.c or file2.c if you do not
*  // want a seperate .c file for the heap instance.
*  PWU_HEAP_BODY(MyStruct);
*  PWU_HEAP_CREATE(myHeap, MyStruct, 100000);
*
*  //--------------------------------------------
*  // in file1.c
*  //--------------------------------------------
*  #include "myHeap.h"
*  void doHeapStuff1() {
*       MyStruct *pChunk = PWU_HEAP_ALLOCITEMS(myHeap, numItems);
*  }
*
*  //--------------------------------------------
*  // in file2.c
*  //--------------------------------------------
*  #include "myHeap.h"
*  void doHeapStuff2() {
*       MyStruct *pChunk = PWU_HEAP_ALLOCITEMS(myHeap, numItems);
*  }
* \endcode
*
* \sa PWU_HEAP_BODY, PWU_HEAP_CREATE, PWU_HEAP_ALLOCITEMS
*/
#define PWU_HEAP_EXTERN(heapName, datatype) \
	extern PwuHeap_##datatype heapName;


/***************************************************************************/
/*! \brief Implements the heap functions for a given heap-type.
*
* \param datatype  The heap item type-name.
*
* \note This macro must be used in conjunction with a corresponding
*       PWU_HEAP_HEADER statement.
*
* \par
* \code
*  // usually placed in a .h file
*  PWU_HEAP_HEADER(MyStruct);
*
*  // placed in a .c file
*  PWU_HEAP_BODY(MyStruct);
*
*  // declares 2 distinct heaps of this item-type
*  PWU_HEAP_CREATE(myHeap1, MyStruct, 100000);
*  PWU_HEAP_CREATE(myHeap2, MyStruct, 200000);
*
* \endcode
*
* \sa PWU_HEAP_HEADER, PWU_HEAP_CREATE, PWU_HEAP_ALLOCITEMS
*/
#define PWU_HEAP_BODY(datatype) \
	/*---------------------------------------------------------------------*/ \
	datatype* PwuHeap_##datatype##_allocitems(PwuHeap_##datatype *pHeap, size_t num) \
	{ \
		return (datatype*)PwuHeapAllocitems((PwuHeap*)pHeap, num); \
	}


/***************************************************************************/
/*! \brief Creates a distinct instance of a typed heap.
*
* \param heapName  The heap instance name. Must be unique for a given scope.
* \param datatype  The heap item type-name.
* \param hint The suggested number of heap items per page. This value can be
*             tuned to trade off between unused page-items and the number of
*             pages in a heap. As the number of pages gets large, wasted space
*             will be minimal, but allocations can slow down. Also, the page
*             size should be significantly larger than the average number of
*             items requested per heap allocation.
*
* \note This macro must be preceeded with corresponding PWU_HEAP_HEADER and
*       PWU_HEAP_BODY statements.
*
* \note When using a strict C-compiler, PWU_HEAP_CREATE() can only be used
*       where variable declarations are allowed.
*
* \note "hint" example; I expect to allocate about 10 million PWP_UINT32
*       items in contiguous chunks of 4 to 8 items. I will use a heap segment
*       hint size of 500,000 items (== 10M / 20). This will give a worst-case
*       waste of about 5%. See code snippet.
*
* \par
* \code
*  // declare a PWP_UINT32 heap
*  PWU_HEAP_CREATE(heap, PWP_UINT32, 500000);
*
*  PWP_UINT32 *pChunk;
*  for(int ii=0; ii < 10000000; ++ii) {
*      numItems = calcNumItems(4, 8);
*      pChunk = PWU_HEAP_ALLOCITEMS(heap, numItems);
*      // do something with pChunk
*  }
*
*  // All done! Release heap memory
*  PWU_HEAP_DESTROY(heap);
* \endcode
*
* \sa PWU_HEAP_HEADER, PWU_HEAP_BODY, PWU_HEAP_ALLOCITEMS
*/
#define PWU_HEAP_CREATE(heapName, datatype, hint) \
	PwuHeap_##datatype heapName = { \
		{ \
			#heapName,                        /* PwuHeap.name */ \
			hint,                             /* PwuHeap.hint */ \
			sizeof(datatype),                 /* PwuHeap.itemSize */ \
			sizeof(PwuHeap_##datatype##_seg), /* PwuHeap.segHdrSize */ \
			0,                                /* PwuHeap.pSeg */ \
            0,                                /* PwuHeap.pFullSeg */ \
		}, \
		PwuHeap_##datatype##_allocitems /* allocItemsFunc_##datatype.allocItems */ \
	}


/***************************************************************************/
/*! \brief Allocates a contiguous block of items from a heap.
*
* \param heapName  The heap instance name.
* \param num  The number of contiguous items to allocate in heap.
*
* \note This macro must be preceeded with corresponding PWU_HEAP_HEADER,
*       PWU_HEAP_BODY, and PWU_HEAP_CREATE statements.
*
* \note The returned pointer is guaranteed to stay valid for the heap's
*       lifetime.
*
* \sa PWU_HEAP_HEADER, PWU_HEAP_BODY, PWU_HEAP_CREATE
*/
#define PWU_HEAP_ALLOCITEMS(heapName, num) \
	(*heapName.allocItems)(&heapName, num)


/***************************************************************************/
/*! \brief Changes a heap's suggested items per page.
*
* \param heapName  The heap instance name.
* \param h The suggested number of heap items per page. Only pages allocated
*          after this call will use this value. Existing pages will NOT be
*          resized.
*
* \note This macro can be used any time after a heap has been created using
*       PWU_HEAP_CREATE.
*
* \sa PWU_HEAP_HEADER, PWU_HEAP_BODY, PWU_HEAP_CREATE
*/
#define PWU_HEAP_SET_HINT(heapName, h) \
	heapName.heap.hint = h


/***************************************************************************/
/*! \brief Release all memory associated with a heap.
*
* \param heapName  The heap instance name.
*
* \note Any pointers obtained with PWU_HEAP_ALLOCITEMS will be invalid.
*
* \sa PWU_HEAP_HEADER, PWU_HEAP_BODY, PWU_HEAP_CREATE
*/
#define PWU_HEAP_DESTROY(heapName) \
	PwuHeapDtor((PwuHeap*)&heapName)


/***************************************************************************/
/*! \brief Obtain heap usage statistics.
*
* \param heapName  The heap instance name.
* \param stats The heap stats struct of type PwuHeapStats.
*
* \note This is intended for debugging use only.
*
* \par
* \code
*  PwuHeapStats stats;
*  PWU_HEAP_CREATE(heap, PWP_UINT32, 500000);
*
*  for(looping) {
*      ..snip.. // use heap here
*  }
*
*  // get heap statistics before destroying it...
*  PWU_HEAP_CALCSTATS(heap, stats);
*
*  PWU_HEAP_DESTROY(heap);
* \endcode
*
* \sa PWU_HEAP_HEADER, PWU_HEAP_BODY, PWU_HEAP_CREATE
*/
#define PWU_HEAP_CALCSTATS(heapName, stats) \
    PwuHeapCalcStats((PwuHeap*)&heapName, &stats)


/***************************************************************************/
/*                     DECLARE STANDARD HEAP TYPES                         */
/***************************************************************************/
/*! \cond sdkINTERNALS */
PWU_HEAP_HEADER(PWP_INT8);
PWU_HEAP_HEADER(PWP_UINT8);
PWU_HEAP_HEADER(PWP_INT16);
PWU_HEAP_HEADER(PWP_UINT16);
PWU_HEAP_HEADER(PWP_INT32);
PWU_HEAP_HEADER(PWP_UINT32);
PWU_HEAP_HEADER(PWP_INT64);
PWU_HEAP_HEADER(PWP_UINT64);
PWU_HEAP_HEADER(PWP_FLOAT);
PWU_HEAP_HEADER(PWP_REAL);
/*! \endcond */ /* sdkINTERNALS */

/*! @} */ /* DOXGRP_APIPWP_UTILS_HEAP */

/*! @} */ /* DOXGRP_APIPWP_UTILS */

/*! @} */ /* DOXGRP_APIPWP */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  /* _APIPWPUTILS_H_ */
