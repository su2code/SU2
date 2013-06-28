/****************************************************************************
 *
 * PWP Plugin example
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#include <string.h>

#define SHOW_PWP_MESSAGES
#include "apiPWP.h"
#include "apiCAEP.h"

#include "apiPWPUtils.h"
#include "apiUtils.h"

//************************************************
// impl-defined plugin version values
//************************************************
#include "rtPwpVersions.h"


/*------------------------------------*/
PWU_RTITEM pwpRtItem[] = {

// include the impl-defined PWP runtime item array data
#include "rtPwpInitItems.h"

    // special items NOT returned bt PwpEnumAPIs() - always LAST!
    /*............................*/
    {
        {PWP_MESSAGECB_DEFAULT, {0,0}},
        0,
    },
    /*............................*/
    {
        {PWP_MESSAGECB_SPY, {0,0}},
        0,
    },
};

PWP_UINT32 totalApiCnt     = ARRAYSIZE(pwpRtItem);
PWP_UINT32 publishedApiCnt = ARRAYSIZE(pwpRtItem) - 2;


/**************************************/
PWP_BOOL PwpInitialize()
{
    PWP_BOOL ret = PWP_TRUE;
    return ret;
}

/**************************************/
PWP_BOOL PwpActivateAPI(const char api[])
{
    return PwuFindPublishedAPI(api) ? PWP_TRUE : PWP_FALSE;
}

/**************************************/
PWP_VOID PwpDestroy()
{
}

/**************************************/
const char* PwpEnumAPIs(PWP_UINT32 ndx, PWP_APIINFO *pInfo)
{
    const char* ret = (ndx < publishedApiCnt) ? pwpRtItem[ndx].apiInfo.name : 0;
    if (ret && pInfo) {
        *pInfo = pwpRtItem[ndx].apiInfo;
    }
    return ret;
}

/**************************************/
PWP_UINT32 PwpGetAPICount()
{
    return publishedApiCnt;
}

/**************************************/
PWP_MESSAGECB PwpGetMessageCallback(const char api[])
{
    return PwuFindApiMsgCB(api);
}

/**************************************/
PWP_VERSIONVAL PwpGetPluginInfo(PWP_PLUGININFO *pInfo)
{
    PWP_VERSIONVAL ret = VERSION_PWP_MAJOR;
    PWP_PLUGININFO info = {
#       include "rtPwpPluginInfo.h"
    };

    info.apiCount = publishedApiCnt;
    info.defCB    = PwuFindApiMsgCB(PWP_MESSAGECB_DEFAULT);
    info.spyCB    = PwuFindApiMsgCB(PWP_MESSAGECB_SPY);

    if (pInfo) {
        *pInfo = info;
    }
    return ret;
}

/**************************************/
PWP_BOOL PwpIsLicensed(const char api[], const PWP_LICENSEDATA *pLicenseData)
{
    return (pLicenseData && PwuFindPublishedAPI(api)) ? PWP_TRUE : PWP_FALSE;
}

/**************************************/
PWP_MESSAGECB PwpSetMessageCallback(const char api[], PWP_MESSAGECB msgCallback)
{
    PWP_MESSAGECB ret = 0;
    PWU_RTITEM *pApiItem = PwuFindTotalAPI(api);
    if (pApiItem) {
        ret = pApiItem->msgCB;
        PwuSendDebugMsg(api, "removing callback", 0);
        pApiItem->msgCB = msgCallback;
        PwuSendDebugMsg(api, "adding callback", 1);
    }
    return ret;
}

/**************************************/
PWP_VOID PwpSetLanguage(const char language[])
{
    if (language && language[0]) {
    }
}
