/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#ifndef _APIUTILS_H_
#define _APIUTILS_H_

#include "site.h"


#ifdef __cplusplus
extern "C" {
#endif


/*! \file
    \brief Base plugin utilities.

    Defines a set of helper/utility macros useful to all plugins.
*/


/*************************************************
               GUID macros
**************************************************/

/*-----------------------------------------------------*/
/*! \brief Builds a 32-bit id value.
*
* This macro is intended to build GUIDs. But, it can be used to combine any 2,
* 16-bit values into a single 32-bit value.
*
* \param hi16
*    Any unsigned 16-bit value. These bits are shifted to the upper word of
*    the resulting 32-bit value.
*
* \param low16
*    Any unsigned 16-bit value. These bits are placed in the lower word of
*    the resulting 32-bit value.
*
* \sa MAKEGUID()
*/
#define MAKEGUID2(hi16,low16) \
            (PWP_UINT32)( ( ((hi16) << 16) & 0xFFFF0000UL) | \
            ((low16) & 0x0000FFFFUL) )


/*-----------------------------------------------------*/
/*! \brief The first group id value available for local, internal use.
*
* Group ids issued by Pointwise will always be \b less than this value.
* The Maximum allowed group id value is 65535.
*
* Unless specified in site.h, \ref PWP_SITE_GROUPID will default to this value.
*
* \sa PWP_SITE_GROUPID
*/
#define PWP_GRPID_USER   (PWP_UINT32)(63000UL)

/*-----------------------------------------------------*/
/*! \brief The default site group name.
*
* Unless specified in site.h, \ref PWP_SITE_GROUPNAME will default to this
* value.
*
* \sa PWP_SITE_GROUPNAME
*/
#define PWP_GROUPNAME_DEFAULT "GroupUndefined"


/*! \cond sdkINTERNALS */
/*-----------------------------------------------------*/
/*! \def PWP_SITE_GROUPID
* \brief The default site group id.
*
* If you plan on releasing a plugin outside your company, you will need to
* obtain a PWP_SITE_GROUPID value from Pointwise support.
*
* You can define your PWP_SITE_GROUPID value in the
* PluginSDK\src\plugins\site.h header file. If a value is not defined, it will
* default to PWP_GRPID_USER.
*
* \sa site.h
*/
#if !defined(PWP_SITE_GROUPID)
    #define PWP_SITE_GROUPID PWP_GRPID_USER
#endif
/*! \endcond sdkINTERNALS */


/*-----------------------------------------------------*/
/*! \def MAKEGUID(usr_id)
*   \brief Builds a 32-bit Globally Unique Id (GUID).
*
* This macro uses the PWP_SITE_GROUPID for the upper, 16-bits of the
* resulting GUID. See PWP_SITE_GROUPID for important information.
*
* \param usr_id
*    Any locally unique, unsigned 16-bit value. The developer is responsible
*    for the local uniqueness of the usr_id values.
*
* \sa PWP_SITE_GROUPID, MAKEGUID2()
*/
#define MAKEGUID(usr_id) MAKEGUID2(PWP_SITE_GROUPID,usr_id)


/*! \cond sdkINTERNALS */
/*-----------------------------------------------------*/
/*! \def PWP_SITE_GROUPNAME
*   \brief The default company or group name.
*
* If you define a PWP_SITE_GROUPNAME value in the PluginSDK\src\plugins\site.h
* header file, it will be used for building the plugin. If a value is not
* defined, it will default to PWP_GROUPNAME_DEFAULT.
*
* This value is returned by PwpGetPluginInfo() in the PWP_PLUGININFO::author
* data member.
*
* \sa PWP_PLUGININFO, PwpGetPluginInfo()
*/
#if !defined(PWP_SITE_GROUPNAME)
    #define PWP_SITE_GROUPNAME PWP_GROUPNAME_DEFAULT
#endif
/*! \endcond sdkINTERNALS */


/*-----------------------------------------------------*/
/*! \brief Calculates the size of a statically declared array
*
* \note The array MUST be statically declared in the scope where this macro is used.
*
* \par "Sample GOOD usage:"
* \code
*   static int globalArray[22];
*
*   void func()
*   {
*       // globalArray is still in scope here
*       // hence, sz == 22
*       int sz = ARRAYSIZE(globalArray);
*   }
* \endcode
*
* \par "Sample BAD usage:"
* \code
*   void func()
*   {
*       int localArray[22];
*       func2(localArray);
*   }
*   
*   void func2(int *array)
*   {
*       // localArray is NOT in scope here
*       // hence, sz != 22
*       int sz = ARRAYSIZE(array);
*   }
* \endcode
*/
#define ARRAYSIZE(arrname)  (PWP_UINT32)(sizeof(arrname)/sizeof(arrname[0]))



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  /* _APIUTILS_H_ */
