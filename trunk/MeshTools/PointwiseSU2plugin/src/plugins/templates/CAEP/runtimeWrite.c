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


#if 0
/**************************************/
static void
stepN(CAEP_RTITEM *pRti)
{
    PWP_UINT32 cnt = 1; /* the # of MINOR progress sub-steps */
    if (caeuProgressBeginStep(pRti, cnt)) {
        while (cnt--) {
            /*
            // PUT YOUR SUB-STEP OUTPUT LOGIC HERE
            */
            /* incr to next sub-step */
            caeuProgressIncr(pRti);
        }
        caeuProgressEndStep(pRti);
    }
}
#endif


/**************************************/
PWP_BOOL
runtimeWrite(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model,
             const CAEP_WRITEINFO *pWriteInfo)
{
    PWP_BOOL ret = PWP_FALSE;

    if (pRti && model && pWriteInfo) {
        PWP_UINT32 cnt = 1; /* the # of MAJOR progress steps */

        if (caeuProgressInit(pRti, cnt)) {
            /*
            // PUT YOUR MAJOR-STEP OUTPUT LOGIC HERE
            step1(pRti);
            step2(pRti);
            ...
            stepN(pRti);
            */
            caeuProgressEnd(pRti, ret);
            ret = !pRti->opAborted;
        }
    }
    return ret;
}
