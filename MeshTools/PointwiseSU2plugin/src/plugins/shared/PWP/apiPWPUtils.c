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
#include "apiPWPUtils.h"


/*------------------------------------*/
static PWU_RTITEM*
PwuFindAPI(const char api[], PWP_UINT32 cnt)
{
    PWU_RTITEM *ret = 0;
    if (api && api[0]) {
        PWP_UINT32 ii;
        for (ii=0; ii < cnt; ++ii) {
            if (0 == strcmp(pwpRtItem[ii].apiInfo.name, api)) {
                ret = &(pwpRtItem[ii]);
                break;
            }
        }
    }
    return ret;
}

/**************************************/
PWU_RTITEM*
PwuFindTotalAPI(const char api[])
{
    return PwuFindAPI(api, totalApiCnt);
}

/**************************************/
PWU_RTITEM*
PwuFindPublishedAPI(const char api[])
{
    return PwuFindAPI(api, publishedApiCnt);
}

/**************************************/
PWP_MESSAGECB
PwuFindApiMsgCB(const char api[])
{
    PWU_RTITEM *pApiInfo = PwuFindTotalAPI(api);
    return (pApiInfo? pApiInfo->msgCB : 0);
}


//***********************************************************************
//***********************************************************************
//        send a generic message back to framework
//***********************************************************************
//***********************************************************************

/**************************************/
PWP_UINT32
PwuSendMsg(const char api[], PWP_ENUM_MSGID id, void *pMsg)
{
    PWP_UINT32 ret = 0;
    PWP_MESSAGECB msgCB = PwuFindApiMsgCB(api);
    if (!msgCB) {
        // api does NOT have a CB, use default CB
        msgCB = PwuFindApiMsgCB(PWP_MESSAGECB_DEFAULT);
    }

    // invoke CB if defined
    if (msgCB) {
        ret = msgCB(api, id, pMsg);
    }

    // always invoke spy CB if defined. ignore ret.
    msgCB = PwuFindApiMsgCB(PWP_MESSAGECB_SPY);
    if (msgCB) {
        msgCB(api, id, pMsg);
    }
    return ret;
}

//***********************************************************************
//***********************************************************************
//        bundle and send a text message back to framework
//***********************************************************************
//***********************************************************************

/*------------------------------------*/
static void
PwuSendTextMsg(const char api[], PWP_ENUM_MSGID id, const char txt[], PWP_UINT32 code)
{
    PWP_MSG_TEXT msg;
    msg.code = code;
    msg.text = (txt ? txt : "");
    PwuSendMsg(api, id, (void*)&msg);
}

/**************************************/
void
PwuSendDebugMsg(const char api[], const char txt[], PWP_UINT32 code)
{
    PwuSendTextMsg(api, PWP_MSGID_DEBUG, txt, code);
}

/**************************************/
void
PwuSendInfoMsg(const char api[], const char txt[], PWP_UINT32 code)
{
    PwuSendTextMsg(api, PWP_MSGID_INFO, txt, code);
}

/**************************************/
void
PwuSendWarningMsg(const char api[], const char txt[], PWP_UINT32 code)
{
    PwuSendTextMsg(api, PWP_MSGID_WARNING, txt, code);
}

/**************************************/
void
PwuSendErrorMsg(const char api[], const char txt[], PWP_UINT32 code)
{
    PwuSendTextMsg(api, PWP_MSGID_ERROR, txt, code);
}

//***********************************************************************
//***********************************************************************
//        bundle and send a progress message back to framework
//***********************************************************************
//***********************************************************************

/*------------------------------------*/
static PWP_UINT32
PwuSendProgressMsg (const char api[], PWP_ENUM_MSGID id, PWP_UINT32 value)
{
    PWP_MSG_PROGRESS msg;
    msg.value = value;
    return PwuSendMsg(api, id, (void*)&msg);
}

/**************************************/
PWP_BOOL
PwuProgressBegin (const char api[], PWP_UINT32 totalSteps)
{
    //For Message id         msg.value =        return
    //--------------------   ----------------   ---------------------
    //PWP_MSGID_PROGBEGIN    total steps        >0 if op can continue
    return PwuSendProgressMsg (api, PWP_MSGID_PROGBEGIN, totalSteps)
             ? PWP_TRUE : PWP_FALSE;
}

/**************************************/
void
PwuProgressEnd (const char api[], PWP_BOOL ok)
{
    //For Message id         msg.value =        return
    //--------------------   ----------------   ---------------------
    //PWP_MSGID_PROGEND      0=fail/>0=ok       nothing
    PwuSendProgressMsg (api, PWP_MSGID_PROGEND, (ok ? 1 : 0));
}

/**************************************/
PWP_BOOL
PwuProgressStatus (const char api[], PWP_UINT32 complete, PWP_UINT32 total)
{
    //For Message id         msg.value =        return
    //--------------------   ----------------   ---------------------
    //PWP_MSGID_PROGSTATUS   0-100              >0 if op can continue
    PWP_UINT32 percent = 100;
    if ((0 != total) && (complete < total)) {
        percent = (PWP_UINT32)((PWP_FLOAT)100 * (PWP_FLOAT)complete
                                              / (PWP_FLOAT)total);
    }
    return PwuSendProgressMsg (api, PWP_MSGID_PROGSTATUS, percent)
             ? PWP_TRUE : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwuProgressNextStep (const char api[])
{
    //For Message id         msg.value =        return
    //--------------------   ----------------   ---------------------
    //PWP_MSGID_PROGSTATUS   -1                 >0 if op can continue
    return PwuSendProgressMsg (api, PWP_MSGID_PROGSTATUS, (PWP_UINT32)-1)
             ? PWP_TRUE : PWP_FALSE;
}

/**************************************/
PWP_BOOL
PwuProgressQuit (const char api[])
{
    //For Message id         msg.value =        return
    //--------------------   ----------------   ---------------------
    //PWP_MSGID_PROGQUIT     not used           >0 if op canceled
    return PwuSendProgressMsg (api, PWP_MSGID_PROGQUIT, 0)
             ? PWP_TRUE : PWP_FALSE;
}



//***********************************************************************
//***********************************************************************
//          unformatted file I/O helper functions
//***********************************************************************
//***********************************************************************

#define STAT_MASK    ((PWP_UINT32)(0xFFFFFFF0))
#define STAT_BASE    ((PWP_UINT32)(0x9ABCDEF0))
#define STAT_OPEN    ((PWP_UINT32)(STAT_BASE | 0x1))
#define STAT_CLOSED  ((PWP_UINT32)(STAT_BASE | 0x0))

#define UDATA_ISINIT(pUD)   (STAT_BASE == (STAT_MASK & (pUD)->status))
#define UDATA_ISOPEN(pUD)   ((pUD) && (STAT_OPEN == (pUD)->status) && (pUD)->fp && !(pUD)->hadError)
#define UDATA_INREC(pUD)    (UDATA_ISOPEN(pUD) && (pUD)->inRec)
#define UDATA_ISCLOSED(pUD) ((pUD) && (STAT_CLOSED == (pUD)->status))


/*------------------------------------*/
static const void*
pwuApplyEndian(PWU_ENDIANNESS endianness, const void *buf, size_t size)
{
#   define ENDIAN_MAXSIZE 32
    const void *ret = 0;
    // buf size must be even
    if (buf && !(size & 0x1) && (ENDIAN_MAXSIZE >= size)) {
        PWU_ENDIANNESS osEndianness = PwuGetOsEndianness();
        switch (endianness) {
            case PWU_ENDIAN_LITTLE:
            case PWU_ENDIAN_BIG:
                if (osEndianness == endianness) {
                    // no need to do swap
                    ret = buf;
                    break;
                }
                // drop through

            case PWU_ENDIAN_FOREIGN: {
                // must swap data order
                static unsigned char tmpBuf[ENDIAN_MAXSIZE];
                const char *cbuf = (const char *)buf;
                size_t hNdx = 0;
                size_t tNdx = size-1;
                while (hNdx < tNdx) {
                    tmpBuf[hNdx]   = cbuf[tNdx];
                    tmpBuf[tNdx--] = cbuf[hNdx++];
                }
                // ret only valid until next call to pwuApplyEndian()
                ret = tmpBuf;
                break; }

            case PWU_ENDIAN_NATIVE:
                // no need to do swap
                ret = buf;
                break;

            case PWU_ENDIAN_ERROR:
            default:
                // huh?
                break;
        }
    }
    return ret;
#   undef ENDIAN_MAXSIZE
}

/*------------------------------------*/
static PWP_BOOL
unfUDataInit(PWU_UNFDATA *pUData)
{
    if (pUData) {
        memset(pUData, 0, sizeof(PWU_UNFDATA));
        pUData->status = STAT_OPEN;
        pUData->endianness = PWU_ENDIAN_NATIVE;
    }
    return pUData ? PWP_TRUE : PWP_FALSE;
}

/*------------------------------------*/
static PWP_BOOL
unfHdrLenWrite(PWU_UNFDATA *pUData)
{
    PWP_BOOL ret = PWP_FALSE;
    if (UDATA_INREC(pUData)) {
        const void *p = pwuApplyEndian(pUData->endianness, &(pUData->recBytes),
                                       sizeof(pUData->recBytes));
        size_t cnt = pwpFileWrite(p, sizeof(pUData->recBytes), 1, pUData->fp);
        if (1 == cnt) {
            ret = PWP_TRUE;
        }
        if (!ret) {
            pUData->hadError = PWP_TRUE;
        }
    }
    return ret;
}

/**************************************/
PWU_ENDIANNESS
PwuGetOsEndianness(void)
{
    union endian_test_t {
        PWP_UINT32    uint;
        unsigned char ch[sizeof(PWP_UINT32)];
    }
    const etest = { 0xAABBCCDD };
    return (0xAA == etest.ch[0]) ? PWU_ENDIAN_BIG : PWU_ENDIAN_LITTLE;
}


/***************************************************************************/
/***************************************************************************/
/*                        UNFORMATTED FORTRAN FILE UTILS                   */
/***************************************************************************/
/***************************************************************************/

/**************************************/
PWU_ENDIANNESS
PwuUnfFileSetEndianness (PWU_UNFDATA *pUData, PWU_ENDIANNESS endianness)
{
    PWU_ENDIANNESS ret = PWU_ENDIAN_ERROR;
    if (UDATA_ISOPEN(pUData) && (PWU_ENDIAN_ERROR != endianness)) {
        ret = pUData->endianness;
        pUData->endianness = endianness;
    }
    return ret;
}

/**************************************/
PWU_ENDIANNESS
PwuUnfFileGetEndianness (PWU_UNFDATA *pUData)
{
    return pUData ? pUData->endianness : PWU_ENDIAN_ERROR;
}

/**************************************/
PWP_BOOL
PwuUnfFileBegin(FILE *fp, PWU_UNFDATA *pUData)
{
    PWP_BOOL ret = PWP_FALSE;
    if (unfUDataInit(pUData) && fp) {
        pUData->fp = fp;
        ret = PWP_TRUE;
    }
    return ret;
}

/**************************************/
PWP_BOOL
PwuUnfRecEnd(PWU_UNFDATA *pUData)
{
    PWP_BOOL ret = PWP_FALSE;
    if (UDATA_INREC(pUData)) {
        // write final rec len value to footer
        if (unfHdrLenWrite(pUData)) {
            // save current pos. will restore it later
            sysFILEPOS fPosSave;
            if (pwpFileGetpos(pUData->fp, &fPosSave)) {
                // fail
            }
            // set pos to rec header location saved in PwuUnfRecBegin()
            else if (pwpFileSetpos(pUData->fp, &(pUData->fPos))) {
                // fail
            }
            // write final rec len value to header
            else if (!unfHdrLenWrite(pUData)) {
                // fail
            }
            // set pos to saved location
            else if (pwpFileSetpos(pUData->fp, &fPosSave)) {
                // fail
            }
            else {
                ret = PWP_TRUE;
            }
        }
        // update running totals
        pUData->totRecBytes += pUData->recBytes;
        ++pUData->recCnt;

        if (!ret) {
            pUData->hadError = PWP_TRUE;
        }
    }

    if (pUData) {
        // clear for next rec
        pUData->inRec = PWP_FALSE;
        // do NOT reset recBytes so that PwuUnfRecBytes()
        // can be used to get last record's byte count.
        //pUData->recBytes = 0;
    }
    return ret;
}

/**************************************/
PWP_BOOL
PwuUnfRecBegin(PWU_UNFDATA *pUData)
{
    PWP_BOOL ret = PWP_FALSE;
    if (UDATA_ISOPEN(pUData)) {
        // auto-close rec (if any)
        PwuUnfRecEnd(pUData);
        // PwuUnfRecEnd() does NOT reset recBytes
        pUData->recBytes = 0;
        // save current file pos. will write the final rec length later
        if (0 == pwpFileGetpos (pUData->fp, &(pUData->fPos))) {
            // must be true for call to unfHdrLenWrite()
            pUData->inRec = PWP_TRUE;
            // write dummy value... will set final value in PwuUnfRecEnd()
            ret = pUData->inRec = unfHdrLenWrite(pUData);
        }
        if (!ret) {
            pUData->hadError = PWP_TRUE;
        }
    }
    return ret;
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteArr(PWU_UNFDATA *pUData, const void *arr, size_t itemSize,
                  size_t itemCnt)
{
    PWP_BOOL ret = PWP_FALSE;
    if (UDATA_INREC(pUData)) {
        if (itemCnt == pwpFileWrite(arr, itemSize, itemCnt, pUData->fp)) {
            pUData->recBytes += (PWP_UINT32)(itemCnt * itemSize);
            ret = PWP_TRUE;
        }
    }
    if (!ret) {
        pUData->hadError = PWP_TRUE;
    }
    return ret;
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteBuf(PWU_UNFDATA *pUData, const void *buf, size_t size)
{
    return PwuUnfRecWriteArr(pUData, buf, size, 1);
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteEndianBuf(PWU_UNFDATA *pUData, const void *buf, size_t size)
{
    const void *p = pwuApplyEndian(pUData->endianness, buf, size);
    return PwuUnfRecWriteArr(pUData, p, size, 1);
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteINT(PWU_UNFDATA *pUData, PWP_INT val)
{
    return PwuUnfRecWriteEndianBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteUINT(PWU_UNFDATA *pUData, PWP_UINT val)
{
    return PwuUnfRecWriteEndianBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteINT8(PWU_UNFDATA *pUData, PWP_INT8 val)
{
    return PwuUnfRecWriteBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteUINT8(PWU_UNFDATA *pUData, PWP_UINT8 val)
{
    return PwuUnfRecWriteBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteINT16(PWU_UNFDATA *pUData, PWP_INT16 val)
{
    return PwuUnfRecWriteEndianBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteUINT16(PWU_UNFDATA *pUData, PWP_UINT16 val)
{
    return PwuUnfRecWriteEndianBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteINT32(PWU_UNFDATA *pUData, PWP_INT32 val)
{
    return PwuUnfRecWriteEndianBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteUINT32(PWU_UNFDATA *pUData, PWP_UINT32 val)
{
    return PwuUnfRecWriteEndianBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteINT64(PWU_UNFDATA *pUData, PWP_INT64 val)
{
    return PwuUnfRecWriteEndianBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteUINT64(PWU_UNFDATA *pUData, PWP_UINT64 val)
{
    return PwuUnfRecWriteEndianBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteFLOAT(PWU_UNFDATA *pUData, PWP_FLOAT val)
{
    return PwuUnfRecWriteEndianBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfRecWriteREAL(PWU_UNFDATA *pUData, PWP_REAL val)
{
    return PwuUnfRecWriteEndianBuf(pUData, &val, sizeof(val));
}

/**************************************/
PWP_BOOL
PwuUnfFileEnd(PWU_UNFDATA *pUData)
{
    PWP_BOOL ret = PWP_FALSE;
    if (pUData) {
        // make sure rec is closed
        PwuUnfRecEnd(pUData);
        pUData->fp = 0;
        pUData->status = STAT_CLOSED;
        ret = !pUData->hadError;
    }
    return ret;
}

/**************************************/
PWP_BOOL
PwuUnfHadError(PWU_UNFDATA *pUData)
{
    return UDATA_ISINIT(pUData) && pUData->hadError;
}

/**************************************/
PWP_UINT32
PwuUnfRecBytes(PWU_UNFDATA *pUData)
{
    return UDATA_ISINIT(pUData) ? pUData->recBytes : 0;
}

/**************************************/
PWP_UINT32
PwuUnfTotRecBytes(PWU_UNFDATA *pUData)
{
    return UDATA_ISINIT(pUData) ? pUData->totRecBytes : 0;
}

/**************************************/
PWP_UINT32
PwuUnfRecCount(PWU_UNFDATA *pUData)
{
    return UDATA_ISINIT(pUData) ? pUData->recCnt : 0;
}


/***************************************************************************/
/***************************************************************************/
/*                           HEAP MANAGER SUBSYSTEM                        */
/***************************************************************************/
/***************************************************************************/

/***************************************************************************/
/*                     IMPLEMENT STANDARD HEAP TYPES                       */
/***************************************************************************/
PWU_HEAP_BODY(PWP_INT8);
PWU_HEAP_BODY(PWP_UINT8);
PWU_HEAP_BODY(PWP_INT16);
PWU_HEAP_BODY(PWP_UINT16);
PWU_HEAP_BODY(PWP_INT32);
PWU_HEAP_BODY(PWP_UINT32);
PWU_HEAP_BODY(PWP_INT64);
PWU_HEAP_BODY(PWP_UINT64);
PWU_HEAP_BODY(PWP_FLOAT);
PWU_HEAP_BODY(PWP_REAL);


/***************************************************************************/
/*                        PRIVATE INTERFACE TO HEAP                        */
/***************************************************************************/

/*---------------------------------------------------------------------*/
static int
PwuHeapNewseg(PwuHeap *pHeap, size_t size)
{
    PwuHeapSeg *pNewSeg = 0;
    if (pHeap) {
        int ii;
        // allow up to 4 attempts if out of memory
        for (ii = 0; (0 < size) && (ii < 4); ++ii) {
            /* subtracting pHeap->itemSize accounts for the single items[1]
               present in the typed segment header. See the PWU_HEAP_HEADER
               macro. */
            size_t bytes =
                pHeap->segHdrSize + (size * pHeap->itemSize) - pHeap->itemSize;
            pNewSeg = (PwuHeapSeg*)malloc(bytes);
            if (pNewSeg) {
                memset(pNewSeg, 0, bytes);
                pNewSeg->pHeap = pHeap;
                pNewSeg->pNextSeg = pHeap->pSeg;
                pNewSeg->size = size;
                pNewSeg->offset = size;
                pHeap->pSeg = pNewSeg;
                break; // all done
            }
            // if we get here, alloc failed. Try again with a smaller size.
            size /= 2;
        }
    }
    return 0 != pNewSeg;
}

/*---------------------------------------------------------------------*/
static void
PwuCheckForFullSegs(PwuHeapSeg **pOwnersSegPtr)
{
    // This logic may be a little obscure... so read carefully:
    //
    //   We have an "owner" that contains a PwuHeapSeg ptr data member. The
    //   pOwnersSegPtr arg points to this storage. The segment pointed to by
    //   *pOwnersSegPtr (aka pThisSeg) is checked. If pThisSeg is full, we
    //   want to "unown" pThisSeg by making the owner point to
    //   pThisSeg->pNextSeg. This removes pThisSeg from the owner's seg chain.
    //   Then we take pThisSeg and make it the new head of the pHeap->pFullSeg
    //   chain.
    //
    while (pOwnersSegPtr && *pOwnersSegPtr) {
        PwuHeapSeg *pThisSeg = (*pOwnersSegPtr);
        PwuHeap *pHeap = pThisSeg->pHeap;
        // is pThisSeg full?
        if (0 == pThisSeg->offset) {
            // Make pThisSeg's owner point to next seg in owner's chain.
            // pThisSeg is now unowned (dangling).
            *pOwnersSegPtr = pThisSeg->pNextSeg;
            // Chain pThisSeg to the current head of the fullseg chain.
            pThisSeg->pNextSeg = pHeap->pFullSeg;
            // Make pThisSeg the new head of the fullseg chain.
            // pThisSeg is no longer dangling.
            pHeap->pFullSeg = pThisSeg;
            // Do NOT change pOwnersSegPtr here. On next pass, we will check
            // the newly owned seg (pThisSeg->pNextSeg) on the next pass.
        }
        else {
            // pThisSeg is NOT full. Make pThisSeg the current "owner" and make
            // another pass.
            pOwnersSegPtr = &pThisSeg->pNextSeg;
        }
    }
}

/*---------------------------------------------------------------------*/
static void*
PwuHeapFindsegitems(PwuHeapSeg *pSeg, size_t num)
{
    void *ret = 0;
    while (pSeg) {
        if (pSeg->offset >= num) {
            PwuHeap *pHeap = pSeg->pHeap;
            size_t segHdrSize = pHeap->segHdrSize;
            size_t itemSize = pHeap->itemSize;
            char *items = (char*)(pSeg) + segHdrSize - itemSize;
            pSeg->offset -= num;
            ret = items + (pSeg->offset * itemSize);
            if (0 == pSeg->offset) {
                PwuCheckForFullSegs(&pHeap->pSeg);
            }
            break;
        }
        // try next seg *without* recursion
        pSeg = pSeg->pNextSeg;
    }
    return ret;
}

/*---------------------------------------------------------------------*/
static void
PwuHeapFreesegChain(PwuHeapSeg *pSeg)
{
    while (pSeg) {
        // save for next pass
        PwuHeapSeg *pNextSeg = pSeg->pNextSeg;
        pSeg->pNextSeg = 0;
        free(pSeg);
        pSeg = pNextSeg;
    }
}

/*---------------------------------------------------------------------*/
static size_t
PwuHeapMergeStats(PwuHeapSeg *pSeg, PwuHeapStats *pStats, size_t *pChainSegCnt)
{
    size_t ret = 0; /* local seg count */
    if (pSeg && pStats) {
        ret = pStats->segCnt; /* save */
        while (pSeg) {
            ++pStats->segCnt;
            pStats->totalUnused += pSeg->offset;
            pStats->totalSize += pSeg->size;
            pSeg = pSeg->pNextSeg;
        }
        ret = pStats->segCnt - ret; /* calc delta */
    }
    if (pChainSegCnt) {
        *pChainSegCnt = ret;
    }
    return ret;
}


/***************************************************************************/
/*                        PUBLIC INTERFACE TO HEAP                         */
/***************************************************************************/

/***************************************************************************/
void*
PwuHeapAllocitems(PwuHeap *pHeap, size_t num)
{
    void *ret = 0;
    ret = PwuHeapFindsegitems(pHeap->pSeg, num);
    if (!ret) {
        size_t segsize = (num > pHeap->hint) ? num : pHeap->hint;
        if (PwuHeapNewseg(pHeap, segsize)) {
            ret = PwuHeapFindsegitems(pHeap->pSeg, num);
        }
    }
    return ret;
}

/***************************************************************************/
void
PwuHeapDtor(PwuHeap *pHeap)
{
    if (pHeap) {
        PwuHeapFreesegChain(pHeap->pSeg);
        pHeap->pSeg = 0;
        PwuHeapFreesegChain(pHeap->pFullSeg);
        pHeap->pFullSeg = 0;
        // do NOT reset other values so pHeap can be reused.
    }
}

/***************************************************************************/
void
PwuHeapCalcStats(PwuHeap *pHeap, PwuHeapStats *pStats)
{
    if (pHeap && pStats) {
        pStats->totalUnused = 0;
        pStats->totalSize = 0;
        pStats->totalUnusedPct = 0;
        pStats->segCnt = 0;
        pStats->fullSegCnt = 0;

        PwuHeapMergeStats(pHeap->pSeg, pStats, 0);
        PwuHeapMergeStats(pHeap->pFullSeg, pStats, &pStats->fullSegCnt);

        if (0 < pStats->totalSize) {
            pStats->totalUnusedPct =
                (float)(100.0 * pStats->totalUnused) / pStats->totalSize;
        }
    }
}
