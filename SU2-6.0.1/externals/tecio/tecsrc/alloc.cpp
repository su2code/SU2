#include "stdafx.h"
#include "MASTER.h"

#define TECPLOTENGINEMODULE

/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/

#define ALLOCMODULE
#include "GLOBAL.h"
#include "ALLOC.h"
#include "TASSERT.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TRACK_MEMORY_USAGE
static size_t memInUse = 0;
static size_t memTotalHighMark = 0;
static size_t memCurrentHighMark = 0;
static size_t memSavedHighMark = 0;
Mutex_pa memMutex;

void initMemoryUsageTracking(void)
{
    REQUIRE(!Thread_ThreadingIsInitialized());
    Thread_InitMutex(&memMutex);
}

void cleanUpMemoryUsageTracking(void)
{
    REQUIRE(!Thread_ThreadingIsInitialized());
    Thread_FreeMutex(&memMutex);
}

void trackMemoryClearHighMark(void)
{
    memCurrentHighMark = memInUse;
}

void trackMemorySaveHighMark(void)
{
    memSavedHighMark = memCurrentHighMark;
}

void trackMemoryAlloc(size_t size)
{
    REQUIRE(memInUse >= 0);

    if (Thread_ThreadingIsInitialized())
        Thread_LockMutex(memMutex);

    memInUse += size;
    if (memInUse > memTotalHighMark)
        memTotalHighMark = memInUse;
    if (memInUse > memCurrentHighMark)
        memCurrentHighMark = memInUse;

    if (Thread_ThreadingIsInitialized())
        Thread_UnlockMutex(memMutex);
}

void trackMemoryFree(size_t size)
{
    if (Thread_ThreadingIsInitialized())
        Thread_LockMutex(memMutex);

    memInUse -= size;

    if (Thread_ThreadingIsInitialized())
        Thread_UnlockMutex(memMutex);

    ENSURE(memInUse >= 0);
}

void getMemoryUsage(size_t* memoryInUse,
                    size_t* memoryCurrentHighMark,
                    size_t* memorySavedHighMark,
                    size_t* memoryTotalHighMark)
{
    REQUIRE(VALID_REF_OR_NULL(memoryInUse));
    REQUIRE(VALID_REF_OR_NULL(memoryCurrentHighMark));
    REQUIRE(VALID_REF_OR_NULL(memorySavedHighMark));
    REQUIRE(VALID_REF_OR_NULL(memoryTotalHighMark));
    if (memoryInUse != NULL)
        *memoryInUse = memInUse;
    if (memoryCurrentHighMark != NULL)
        *memoryCurrentHighMark = memCurrentHighMark;
    if (memorySavedHighMark != NULL)
        *memorySavedHighMark = memSavedHighMark;
    if (memoryTotalHighMark != NULL)
        *memoryTotalHighMark = memTotalHighMark;
}
#endif

#if defined MSWIN && defined ALLOC_HEAP
#define HEAPMIN 512
#endif

#if defined MSWIN && defined ALLOC_HEAP
/**
 */
void *MSWinAlloc(DWORD nSize)
{
    long *pMem = NULL;
    if (nSize < HEAPMIN)
        pMem = (long *)malloc(sizeof(long) + nSize);
    else
        pMem = (long *)HeapAlloc(GetProcessHeap(), NULL, sizeof(long) + nSize);
    if (pMem)
        pMem[0] = nSize;
    return (void *)&(pMem[1]);
}
#endif

#if defined MSWIN && defined ALLOC_HEAP
/**
 */
void MSWinFree(void *pMem)
{
    REQUIRE(VALID_REF(pMem));
    if (pMem)
    {
        long *pMemLong = &(((long *)pMem)[-1]);
        if (pMemLong[0] < HEAPMIN)
            free((void *)pMemLong);
        else
            HeapFree(GetProcessHeap(), NULL, (void *)pMemLong);
    }
}
#endif
