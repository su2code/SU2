#ifndef _APIPWP_H_
#define _APIPWP_H_
/****************************************************************************
 *
 * Pointwise Plugin API (PWP-API) v1.0
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

/*! \cond sdkINTERNALS */

//*********************************************
//   platform-dependent export defines
//*********************************************
#if defined(WINDOWS)
#  define winmsgSTRING2(ss)  #ss
#  define winmsgSTRING(ss)  winmsgSTRING2(ss)
#  define PW_DLLEXPORT __declspec(dllexport)
#  define PW_DLLIMPORT __declspec(dllimport)
#else
#  define PW_DLLEXPORT
#  define PW_DLLIMPORT
#endif /* WINDOWS */

#undef PWP_PROTOTYPE_DECL

#ifdef BUILD_PWPLUGIN_DYNLIB

#include "apiUtils.h"

// Building a plugin dyn lib module.
// Don't need this stuff when building GgPlugin
// for Pointwise distro.

#if defined(WINDOWS)
#  if defined(SHOW_PWP_MESSAGES)
#    pragma message ("PWP_SITE_GROUPNAME='" PWP_SITE_GROUPNAME "'")
#    pragma message ("PWP_SITE_GROUPID='" winmsgSTRING(PWP_SITE_GROUPID) "'")
#  endif /* SHOW_PWP_MESSAGES */
#endif /* WINDOWS */

//*********************************************
/*! \brief Cross-platform function export declaration.
*/
# define PWP_PROTOTYPE_DECL PW_DLLEXPORT

#else // building pointwise distro

/*! \brief Turn API prototypes into function typedefs so we can
    use as function pointers in Pointwise.
*/
# define PWP_PROTOTYPE_DECL typedef

#endif

/*! \endcond */ /* sdkINTERNALS */


#if defined(WINDOWS) && defined(SHOW_PWP_MESSAGES)
#  pragma message ("PWP_PROTOTYPE_DECL='" winmsgSTRING(PWP_PROTOTYPE_DECL) "'")
#endif /* WINDOWS */


#ifdef __cplusplus
extern "C" {
#endif

/*! \file
 \brief Pointwise Plugin API (PWP-API)

 This defines the minimum API required for all Pointwise plugins.

 \sa \link DOXGRP_APIPWP PWP-API specification \endlink
 \sa \link DOXGRP_APICAEP CAEP-API specification \endlink
*/


/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWP Pointwise Plugin API Specification (PWP-API)
*
*  For a plugin to be compatible with Pointwise, a dynamic library MUST
*  implement the PWP-API. There will be additional API requirements for
*  more specialized plugins.
*
* \sa \link DOXGRP_APICAEP CAEP-API specification \endlink
*
*  @{
*/

/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWP_BASENAMES Supported Plugin API Basenames
*  @{
*   The supported Plugin API specification base names. A full API spec
*   name is in the format: basename/majorVer.minorVer
*   (i.e. "Plugin-PWP/1.0")
*
* \sa PwpActivateAPI(), PwpIsLicensed(), PwpGetMessageCallback(),
*     PwpSetMessageCallback(), PWP_MESSAGECB
*/

/*! \brief The PWP API specification base name.
*/
#define PWP_API_PLUGIN      "Plugin-PWP"

/*! @} */ /* DOXGRP_APIPWP_BASENAMES */


/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWP_XPLATFORMTYPES PWP-API Cross-Platform Base Types
*  @{
*   Cross-platform "C" datatypes. Cross-platform data issues can be minimized
*   by using the following datatypes internally for plugins. All API calls use
*   these types.
*/
/*! \brief integer same size as void* */
typedef long                PWP_INT;
/*! \brief unsigned integer same size as void* */
typedef unsigned long       PWP_UINT;
/*! \brief maximum valid PWP_UINT value */
#define PWP_UINT_MAX        ((PWP_UINT)(~0))
/*! \brief "undefined" PWP_UINT value */
#define PWP_UINT_UNDEF      PWP_UINT_MAX

/* Bit length noted ints */
/*! \brief 8-bit integer */
typedef signed char         PWP_INT8;  
/*! \brief 8-bit unsigned integer */
typedef unsigned char       PWP_UINT8; 
/*! \brief maximum valid PWP_UINT8 value */
#define PWP_UINT8_MAX       ((PWP_UINT8)(~0))
/*! \brief "undefined" PWP_UINT8 value */
#define PWP_UINT8_UNDEF     PWP_UINT8_MAX

/*! \brief 16-bit integer */
typedef short               PWP_INT16; 
/*! \brief 16-bit unsigned integer */
typedef unsigned short      PWP_UINT16;
/*! \brief maximum valid PWP_UINT16 value */
#define PWP_UINT16_MAX      ((PWP_UINT16)(~0))
/*! \brief "undefined" PWP_UINT16 value */
#define PWP_UINT16_UNDEF    PWP_UINT16_MAX

/*! \brief 32-bit integer */
typedef int                 PWP_INT32; 
/*! \brief 32-bit unsigned integer */
typedef unsigned int        PWP_UINT32;
/*! \brief maximum valid PWP_UINT32 value */
#define PWP_UINT32_MAX      ((PWP_UINT32)(~0))
/*! \brief "undefined" PWP_UINT32 value */
#define PWP_UINT32_UNDEF    PWP_UINT32_MAX

/*! \brief 64-bit integer */
typedef long long           PWP_INT64; 
/*! \brief 64-bit unsigned integer */
typedef unsigned long long  PWP_UINT64;
/*! \brief maximum valid PWP_UINT64 value */
#define PWP_UINT64_MAX      ((PWP_UINT64)(~0))
/*! \brief "undefined" PWP_UINT64 value */
#define PWP_UINT64_UNDEF    PWP_UINT64_MAX

/*! \brief 32-bit real */
typedef float               PWP_FLOAT; 
/*! \brief 64-bit real */
typedef double              PWP_REAL;  

/*! \brief logical value */
typedef int                 PWP_BOOL;  
/*! \brief PWP_BOOL logical "false" value */
#define PWP_FALSE (0)         
/*! \brief PWP_BOOL logical "true" value */
#define PWP_TRUE  (!PWP_FALSE)

/*! \brief Cast a value to a PWP_BOOL value (PWP_TRUE or PWP_FALSE) */
#define PWP_CAST_BOOL(v)    ((v) ? PWP_TRUE : PWP_FALSE)

/*! \brief no value */
typedef void                PWP_VOID; 

/*! @} */ /* DOXGRP_APIPWP_XPLATFORMTYPES */


/*! \cond sdkINTERNALS */

/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWP_DATAHANDLE_HELPERS PWP-API Data Handle Helper Macros
*  @{
*   Base data handle helper macros. The are used to implement opaque,
*   API-specific, data handles. See PWGM_HGRIDMODEL for an example.
*/

/*---------------------------------------------------------*/
/*! \brief Declares a root-level, strongly-typed data handle type.
*
* \param name  The handle type-name.
*/
#define PWP_DECLARE_HANDLE(name) \
                typedef struct name##_t { \
                    /*! \cond */ \
                    int unused; \
                    /*! \endcond */ \
                } * name
/*! \brief Test the validity of a PWP_DECLARE_HANDLE() handle. */
#define PWP_HANDLE_ISVALID(h)    (0 != h)
/*! \brief Static init value for a PWP_DECLARE_HANDLE() handle. */
#define PWP_HANDLE_INIT          0
/*! \brief Runtime set of a PWP_DECLARE_HANDLE() handle. */
#define PWP_HANDLE_SET(h, v)     h=v

/*! \brief Bad id value. */
#define PWP_BADID     (~((PWP_UINT32)0))
/*! \brief Bad type value. */
#define PWP_BADTYPE   ((unsigned char)~0)


/*---------------------------------------------------------*/
/*! \brief Declares a parented, strongly-typed, element group data handle type.
*
* \param pname Parent type-name declared using PWP_DECLARE_HANDLE(pname)
* \param name  The element group handle type-name.
*/
#define PWP_DECLARE_HELEMGROUP(pname,name) \
                typedef struct name##_t { \
                    /*! \cond */ \
                    pname hP; \
                    PWP_UINT32 id; \
                    /*! \endcond */ \
                } name
/*! \brief Test the validity of a PWP_DECLARE_HELEMGROUP() handle. */
#define PWP_HEGRP_ISVALID(h)    (PWP_HANDLE_ISVALID((h).hP) && \
                                    (PWP_BADID != (h).id))
/*! \brief Static init value for a PWP_DECLARE_HELEMGROUP() handle. */
#define PWP_HEGRP_INIT          {0,PWP_BADID}
/*! \brief Runtime set of a PWP_DECLARE_HELEMGROUP() handle. */
#define PWP_HEGRP_SET(h, p, v)  { (h).hP=(p); (h).id=(v); }
/*! \brief Extract the parent handle from a PWP_DECLARE_HELEMGROUP() handle. */
#define PWP_HEGRP_H(h)          ((h).hP)
/*! \brief Extract the id from a PWP_DECLARE_HELEMGROUP() handle. */
#define PWP_HEGRP_ID(h)         ((h).id)


/*---------------------------------------------------------*/
/*! \brief Declares a sub-element group data handle type.
*
* \param sname Parent type-name declared using
*              PWP_DECLARE_HELEMGROUP(pname,sname)
* \param name  The sub-element group handle type-name.
*/
#define PWP_DECLARE_HEGRPITEM(sname,name) \
                typedef struct name##_t { \
                    /*! \cond */ \
                    sname parent; \
                    unsigned char ptype; \
                    PWP_UINT32 id; \
                    /*! \endcond */ \
                } name
/*! \brief Test the validity of a PWP_DECLARE_HEGRPITEM() handle. */
#define PWP_HEGI_ISVALID(h)       (PWP_HANDLE_ISVALID((h).parent.hP) && \
                                        (PWP_BADID != (h).parent.id) && \
                                        (PWP_BADTYPE != (h).ptype) && \
                                        (PWP_BADID != (h).id))
/*! \brief Static init value for a PWP_DECLARE_HEGRPITEM() handle. */
#define PWP_HEGI_INIT             {{0, PWP_BADID}, PWP_BADTYPE, PWP_BADID}
/*! \brief Runtime set of a PWP_DECLARE_HEGRPITEM() handle. */
#define PWP_HEGI_SET(h, p, pt, pid, v)  { (h).parent.hP=(p); \
                                                (h).parent.id=(pid); \
                                                (h).ptype=(pt); \
                                                (h).id=(v); }
/*! \brief Extract the parent PWP_DECLARE_HELEMGROUP() handle. */
#define PWP_HEGI_H(h)      ((h).parent.hP)
/*! \brief Extract the parent PWP_DECLARE_HELEMGROUP() id. */
#define PWP_HEGI_PID(h)    ((h).parent.id)
/*! \brief Extract the parent-type id from a PWP_DECLARE_HEGRPITEM() handle. */
#define PWP_HEGI_PTYPE(h)  ((h).ptype)
/*! \brief Extract the item id from a PWP_DECLARE_HEGRPITEM() handle. */
#define PWP_HEGI_ID(h)     ((h).id)

/*! @} */ /* DOXGRP_APIPWP_DATAHANDLE_HELPERS */

/*! \endcond */ /* sdkINTERNALS */


/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWP_TYPES Data Types
*  @{
*/

/*---------------------------------------------------------*/
/*! \brief Version data component value.
*
* \sa PWP_VERSION
*/
typedef PWP_UINT32  PWP_VERSIONVAL;

/*---------------------------------------------------------*/
/*! \brief Version data.
*
* \sa PWP_APIINFO, PWP_PLUGININFO
*/
typedef struct PWP_VERSION_t {
    /*! \brief the major version value
    */
    PWP_VERSIONVAL major;

    /*! \brief the minor version value
    */
    PWP_VERSIONVAL minor;
}
PWP_VERSION;

/*---------------------------------------------------------*/
/*! \brief Installation's license data. NOT IMPLEMENTED YET.
*
* \sa PwpIsLicensed
*/
typedef struct PWP_LICENSEDATA_t {
    /*! \cond sdkINTERNALS */
    /*! \brief TBD
    */
    int dummy;
    /*! \endcond sdkINTERNALS */
}
PWP_LICENSEDATA;

/*---------------------------------------------------------*/
/*! \brief Supported PWP-API message ids.
*
* \sa PWP_MESSAGECB, PwpSetMessageCallback, PwpGetMessageCallback
*/
typedef enum PWP_ENUM_MSGID_e {
    /*! \brief Debug text message id (see PWP_MSG_TEXT)
    */
    PWP_MSGID_DEBUG,

    /*! \brief Information text message id (see PWP_MSG_TEXT)
    */
    PWP_MSGID_INFO,

    /*! \brief Non-fatal error text message id (see PWP_MSG_TEXT)
    */
    PWP_MSGID_WARNING,

    /*! \brief Fatal error text message id (see PWP_MSG_TEXT)
    */
    PWP_MSGID_ERROR,

    // PWP_MSG_PROGRESS ids
    /*! \brief Begin progress message id (see PWP_MSG_PROGRESS)
    */
    PWP_MSGID_PROGBEGIN,

    /*! \brief End progress message id (see PWP_MSG_PROGRESS)
    */
    PWP_MSGID_PROGEND,

    /*! \brief Status progress message id (see PWP_MSG_PROGRESS)
    */
    PWP_MSGID_PROGSTATUS,

    /*! \brief Query progress quit message id (see PWP_MSG_PROGRESS)
    */
    PWP_MSGID_PROGQUIT,

    // add new group enum values before this line

    /*! \cond sdkINTERNALS */
    //! for bookkeeping only - NOT a valid "id" value
    PWP_MSGID_SIZE,
    //! for bookkeeping only - NOT a valid "id" value
    PWP_MSGID_LAST = PWP_MSGID_SIZE-1
    /*! \endcond sdkINTERNALS */
}
PWP_ENUM_MSGID;


/*---------------------------------------------------------*/
/*! \brief Message handler callback function signature.
*
* \sa PWP_ENUM_MSGID, PwpSetMessageCallback, PwpGetMessageCallback
*/
typedef PWP_UINT32 (*PWP_MESSAGECB)(const char api[], PWP_ENUM_MSGID id,
                                    void *pMsg);

/*! \brief Special API name used to register the default message handler.
*
*   This special API name can be passed to PwpSetMessageCallback(api).
*
*   The default handler is only invoked if an API specific handler
*   was not set.
*
* \sa PWP_ENUM_MSGID, PWP_MESSAGECB, PwpSetMessageCallback,
*     PwpGetMessageCallback
*/
#define PWP_MESSAGECB_DEFAULT  "@@default"

/*! \brief Special API name used to register the spy message handler.
*
*   This special API name can be passed to PwpSetMessageCallback(api).
*
*   The spy handler is always invoked after the API specific or
*   default handler is invoked. The spy handler cannot modify the
*   return value.
*
* \sa PWP_ENUM_MSGID, PWP_MESSAGECB, PwpSetMessageCallback,
*     PwpGetMessageCallback
*/
#define PWP_MESSAGECB_SPY      "@@spy"


/*---------------------------------------------------------*/
/*! \brief The data sent by plugins for progress messages.
*
* \note The framework supports a 2-level progress interface. The progress
*       message sequence is shown below:
*
* \par
* \code
*   PWP_MSGID_PROGBEGIN(value=3)      // 3 level1 steps
*     PWP_MSGID_PROGSTATUS(value=10)  // step 1 at 10%
*     PWP_MSGID_PROGSTATUS(value=20)  // step 1 at 20%
*        ...
*     PWP_MSGID_PROGSTATUS(value=100) // step 1 at 100%
*
*     PWP_MSGID_PROGSTATUS(value=-1)  // start next step 2of3
*     PWP_MSGID_PROGSTATUS(value=10)  // step 2 at 10%
*     PWP_MSGID_PROGSTATUS(value=20)  // step 2 at 20%
*        ...
*     PWP_MSGID_PROGSTATUS(value=100) // step 2 at 100%
*
*     PWP_MSGID_PROGSTATUS(value=-1)  // start next step 3of3
*     PWP_MSGID_PROGSTATUS(value=10)  // step 3 at 10%
*     PWP_MSGID_PROGSTATUS(value=20)  // step 3 at 20%
*        ...
*     PWP_MSGID_PROGSTATUS(value=100) // step 3 at 100%
*   PWP_MSGID_PROGEND(value=0 or 1)
* \endcode
*/
typedef struct PWP_MSG_PROGRESS_t {
    /*! \brief The progress value
    *
    * \note The "meaning" of value depends on the which message is being used.
    *       The usage is shown in the table below:
    *
    * \par
    * \code
    *                          plugin will set      framework will
    *   For PWP_MSGID_x        msg.value =          return
    *   --------------------   ------------------   ---------------------
    *   PWP_MSGID_PROGBEGIN    total L1 steps       >0 if op can continue
    *   PWP_MSGID_PROGSTATUS   0-100 level2 prog    >0 if op can continue
    *   PWP_MSGID_PROGSTATUS   -1 next level1       >0 if op can continue
    *   PWP_MSGID_PROGEND      0 = fail, !0 = ok    nothing
    *   PWP_MSGID_PROGQUIT     not used             >0 if op canceled
    * \endcode
    */
    PWP_UINT32 value;
}
PWP_MSG_PROGRESS;


/*---------------------------------------------------------*/
/*! \brief The data sent by plugins for text messages.
*
*   The framework will decide how to display the messages.
*/
typedef struct PWP_MSG_TEXT_t {
    /*! \brief API defined message code
    */
    PWP_UINT32 code; 
    /*! \brief API defined message string
    */
    const char *text;
}
PWP_MSG_TEXT;


/*---------------------------------------------------------*/
/*! \brief The API information returned by plugins for each supported API.
*
*   This information is used to integrate the plugin into the framework.
*/
typedef struct PWP_APIINFO_t {
    /*! \brief full API spec name
    */
    const char *name;   
    /*! \brief API spec version
    */
    PWP_VERSION ver;    
}
PWP_APIINFO;


/*---------------------------------------------------------*/
/*! \brief Provides general information about a plugin.
*/
typedef struct PWP_PLUGININFO_t {
    /*! \brief plugin conforms to this PWP-API version
    */
    PWP_VERSION pwpVer;   
    /*! \brief software library release version
    */
    PWP_VERSION libVer;   
    /*! \brief company/author description
    */
    const char* author;   
    /*! \brief support description (phone, web-link).
    */
    const char* support;  
    /*! \brief copyright description
    */
    const char* copyright;
    /*! \brief number of APIs implemented by this plugin
    */
    PWP_UINT32 apiCount;  
    /*! \brief assigned default message callback
    */
    PWP_MESSAGECB defCB;  
    /*! \brief assigned spy message callback
    */
    PWP_MESSAGECB spyCB;  
}
PWP_PLUGININFO;

/*! @} */ /* DOXGRP_APIPWP_TYPES */


/***********************************************************/
/***********************************************************/
/*! \defgroup DOXGRP_APIPWP_FUNCTIONS Functions
*  @{
*/

/*---------------------------------------------------------*/
/*! \brief Initializes the plugin.
*
* Called once and first by the framework immediately after the plugin library
* is opened.
*/
PWP_PROTOTYPE_DECL PWP_BOOL
PwpInitialize();


/*---------------------------------------------------------*/
/*! \brief Activates the plugin for a given API spec.
*
* An API will be activated by the framework before any API specific calls are
* made to the plugin.
*
* \param api 
*    The target api. Must be one of the values returned from
*    PwpEnumAPIs().
*/
PWP_PROTOTYPE_DECL PWP_BOOL
PwpActivateAPI(const char api[]);


/*---------------------------------------------------------*/
/*! \brief Called by framework just before plugin library is closed.
*/
PWP_PROTOTYPE_DECL PWP_VOID
PwpDestroy();


/*---------------------------------------------------------*/
/*! \brief Enumerates the APIs supported by this plugin.
*
* \param ndx
*    The api index starting with 0.
*
* \param pInfo
*    Pointer to the destination PWP_APIINFO buffer.
*
* \return The API spec name (same as pInfo->name). NULL if ndx is not valid.
*
* \sa PwpGetAPICount()
*
* \par "Sample usage:"
* \code
*    PWP_APIINFO apiInfo;
*    PWP_UINT32 ndx = 0;
*    while (caep->PwpEnumAPIs(ndx++, &apiInfo)) {
*        printf ("api: '%s' v%lu.%lu\n",
*            apiInfo.name, apiInfo.ver.major, apiInfo.ver.minor);
*    }
*    // output:
*    // api: 'Exporter-CAE/1.0' v1.0
* \endcode
*/
PWP_PROTOTYPE_DECL const char*
PwpEnumAPIs(PWP_UINT32 ndx, PWP_APIINFO *pInfo);


/*---------------------------------------------------------*/
/*! \brief Get the number of APIs supported by this plugin.
*
* \return The number of published APIs.
*/
PWP_PROTOTYPE_DECL PWP_UINT32
PwpGetAPICount();


/*---------------------------------------------------------*/
/*! \brief Gets the current message callback.
*
* \param api 
*    The target api. Must be one of the values returned from
*    PwpEnumAPIs().
*
* \sa PWP_MESSAGECB, PwpSetMessageCallback
*/
PWP_PROTOTYPE_DECL PWP_MESSAGECB
PwpGetMessageCallback(const char api[]);


/*---------------------------------------------------------*/
/*! \brief Get information about this plugin.
*
* \param pInfo
*    Pointer to the PWP_PLUGININFO buffer.
*
* \return The major PWP-API spec version number.
*           Same as "pInfo->pwpVer.major".
* \sa PWP_PLUGININFO
*/
PWP_PROTOTYPE_DECL PWP_VERSIONVAL
PwpGetPluginInfo(PWP_PLUGININFO *pInfo);


/*---------------------------------------------------------*/
/*! \brief Determines if plugin api is licensed for use on this machine.
*
* \param api 
*    The target api. Must be one of the values returned from
*    PwpEnumAPIs().
*
* \param pLicenseData
*    Pointer to the PWP_LICENSEDATA buffer.
*
* \note This call is for future use. API should always return PWP_TRUE.
* \sa PWP_LICENSEDATA
*/
PWP_PROTOTYPE_DECL PWP_BOOL
PwpIsLicensed(const char api[], const PWP_LICENSEDATA *pLicenseData);


/*---------------------------------------------------------*/
/*! \brief Sets the message callback function for the given api.
*
* The plugin uses this callback to send api specific messages to the framework.
*
* \param api 
*    The target api. Must be one of the values returned from
*    PwpEnumAPIs() or one of the special api targets listed in notes.
*
* \param msgCB
*    Pointer to a function with the signature defined by PWP_MESSAGECB.
*
* \return The previous message callback function.
*
* \sa PwpEnumAPIs, PWP_MESSAGECB_DEFAULT, PWP_MESSAGECB_SPY
*
* \note The available special api targets are:
* \par 
*      PWP_MESSAGECB_DEFAULT - Any messages not routed to an api-specific callback will be
*      sent to this handler.
* \par 
*      PWP_MESSAGECB_SPY - All messages are copied and routed to this handler after any
*      api-specific callbacks are done. Any values returned from the
*      spy are ignored.
*/
PWP_PROTOTYPE_DECL PWP_MESSAGECB
PwpSetMessageCallback(const char api[], PWP_MESSAGECB msgCB);


/*---------------------------------------------------------*/
/*! \brief Sets the active language.
*
* If the specified language is not supported, the plugin should default to
* "us-english".
*
* \param language
*    Pointer to a function with the signature defined by PWP_MESSAGECB.
*
* \note The supported language values are:
* \par
*    "us-english" - United States English
*/
PWP_PROTOTYPE_DECL PWP_VOID
PwpSetLanguage(const char language[]);

/*! @} */ /* DOXGRP_APIPWP_FUNCTIONS */

/*! @} */ /* DOXGRP_APIPWP */


#ifdef __cplusplus
} // extern "C"
#endif

#endif /* !_APIPWP_H_ */
