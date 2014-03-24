/*****************************************************************
 *****************************************************************
 *******                                                  ********
 ****** Copyright (C) 1988-2010 Tecplot, Inc.              *******
 *******                                                  ********
 *****************************************************************
 *****************************************************************/
/* BEGINREMOVEFROMADDON */
/* NOTE: All code contained between comments that look like
 *             BEGINREMOVEFROMADDON
 *             ENDREMOVEFROMADDON
 * are pulled out to create the MASTER.h file used in addons.
 */
/* ENDREMOVEFROMADDON */

#ifndef _MASTER_H_
#define _MASTER_H_

/*
 * Annotations that specify the life cycle of objects returned from functions
 * and input and output parameters sent as function parameters. The following
 * table specifies the meaning in their context. The annotations provide code
 * generation tools with information for building language bindings to various
 * Tecplot 360 and Tecplot SDK related libraries.
 *
 * For purposes of this table the client is one making the call and the service
 * is the recipient.
 *
 * +==================+=========================+=================================================================+
 * | Function Context | Annotation              | Meaning                                                         |
 * |   Result or      |                         |                                                                 |
 * |   Parameter      |                         |                                                                 |
 * |==================+=========================+=================================================================|
 * | Result           | TP_OUT                  | Default for a function return value that does not transfer      |
 * |                  |                         | ownership. Because this is the most common scenario this        |
 * |                  |                         | annotation is implied and never explicitly used in this         |
 * |                  |                         | context.                                                        |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Scalar Result    | TP_GIVES                | Annotates a function scalar return value as one who's ownership |
 * |                  |                         | is transfered to the client. The client is responsible for      |
 * |                  |                         | properly disposing the value.                                   |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Array Result     | TP_ARRAY_GIVES          | Annotates a function array return value as one who's ownership  |
 * |                  |                         | is transfered to the client. The client is responsible for      |
 * |                  |                         | properly disposing the value.                                   |
 * |==================+=========================+=================================================================|
 * | Parameter        | TP_IN                   | Default for a function input parameter value sent to the        |
 * |                  |                         | service. Because this is the most common scenario this          |
 * |                  |                         | annotation is implied and never explicitly used.                |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Parameter        | TP_ACQUIRES             | Annotates a function parameter as one that sends a value to     |
 * |                  |                         | the service through the parameter and acquires shared           |
 * |                  |                         | ownership of the input value with the client. The service is    |
 * |                  |                         | not responsible for disposing the value however it is           |
 * |                  |                         | expected that a symmetric API exists that "releases" the        |
 * |                  |                         | library of this shared ownership. For example:                  |
 * |                  |                         |   void addListener(TP_ACQUIRES Listener& listener);             |
 * |                  |                         |   void removeListener(TP_RELEASES Listener& listener);          |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Parameter        | TP_RELEASES             | Annotates a function parameter as one that sends a value to     |
 * |                  |                         | the service through the parameter and releases previously       |
 * |                  |                         | shared ownership of the                                         |
 * |                  |                         | input value with the client. The service is not responsible     |
 * |                  |                         | for disposing the value however it is expected that a           |
 * |                  |                         | symmetric API exists that "releases" the library of this        |
 * |                  |                         | shared ownership. For example:                                  |
 * |                  |                         |   void addListener(TP_ACQUIRES Listener& listener);             |
 * |                  |                         |   void removeListener(TP_RELEASES Listener& listener);          |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Scalar Parameter | TP_OUT                  | Annotates a function scalar parameter as one that returns a     |
 * |                  |                         | value to the client through the parameter but does not          |
 * |                  |                         | transfer ownership of the output value to the client.           |
 * |                  |                         | The client is not responsible for disposing the value.          |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Scalar Parameter | TP_IN_OUT               | Annotates a function scalar parameter as one that both sends    |
 * |                  |                         | a value to the service and returns a value to the client        |
 * |                  |                         | through the parameter. Ownership of the input value is not      |
 * |                  |                         | transfered to the service nor is ownership of the output value  |
 * |                  |                         | transfered to the client. The service is not responsible for    |
 * |                  |                         | disposing  the input value and the client is not responsible    |
 * |                  |                         | for disposing the output value.                                 |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Array Parameter  | TP_ARRAY_OUT            | Annotates a function array parameter as one that returns a      |
 * |                  |                         | value to the client through the parameter but does not          |
 * |                  |                         | transfer ownership of the output value to the client.           |
 * |                  |                         | The client is not responsible for disposing the value.          |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Array Parameter  | TP_ARRAY_IN_OUT         | Annotates a function array parameter as one that both sends     |
 * |                  |                         | a value to the service and returns a value to the client        |
 * |                  |                         | through the parameter. Ownership of the input value is not      |
 * |                  |                         | transfered to the service nor is ownership of the output value  |
 * |                  |                         | transfered to the client. The service is not responsible for    |
 * |                  |                         | disposing  the input value and the client is not responsible    |
 * |                  |                         | for disposing the output value.                                 |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Scalar Parameter | TP_GIVES                | Annotates a function scalar parameter as one that returns a     |
 * |                  |                         | value to the client through the parameter and transfers         |
 * |                  |                         | ownership of the output value to the client. The client is      |
 * |                  |                         | responsible for properly disposing the value.                   |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Scalar Parameter | TP_RECEIVES             | Annotates a function scalar parameter as one that sends a value |
 * |                  |                         | to the service through the parameter and transfers ownership    |
 * |                  |                         | of the input value to the service. The service is responsible   |
 * |                  |                         | for properly disposing the value.                               |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Scalar Parameter | TP_RECEIVES_GIVES       | Annotates a function scalar parameter as one that both sends    |
 * |                  |                         | a value to the service and returns a value to the client        |
 * |                  |                         | through the  parameter. Ownership of the input value is         |
 * |                  |                         | transfered to the service and ownership of the output value is  |
 * |                  |                         | transfered to the client. The service is responsible for        |
 * |                  |                         | properly disposing the input value and the client is            |
 * |                  |                         | responsible for properly disposing the output value.            |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Array Parameter  | TP_ARRAY_GIVES          | Annotates a function array parameter as one that returns a      |
 * |                  |                         | value to the client through the parameter and transfers         |
 * |                  |                         | ownership of the output value to the client. The client is      |
 * |                  |                         | responsible for properly disposing the value.                   |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Array Parameter  | TP_ARRAY_RECEIVES       | Annotates a function array parameter as one that sends a value  |
 * |                  |                         | to the service through the parameter and transfers ownership    |
 * |                  |                         | of the input value to the service. The service is responsible   |
 * |                  |                         | for properly disposing the value.                               |
 * |------------------+-------------------------+-----------------------------------------------------------------|
 * | Array Parameter  | TP_ARRAY_RECEIVES_GIVES | Annotates a function array parameter as one that both sends     |
 * |                  |                         | a value to the service and returns a value to the client        |
 * |                  |                         | through the  parameter. Ownership of the input value is         |
 * |                  |                         | transfered to the service and ownership of the output value is  |
 * |                  |                         | transfered to the client. The service is responsible for        |
 * |                  |                         | properly disposing the input value and the client is            |
 * |                  |                         | responsible for properly disposing the output value.            |
 * |==================+===================+=======================================================================|
 */

/*
 * First check to make sure that our life-cycle keywords are not in conflict with any system defines.
 */
#if defined TP_ACQUIRES             || \
    defined TP_RELEASES             || \
    defined TP_OUT                  || \
    defined TP_IN_OUT               || \
    defined TP_ARRAY_OUT            || \
    defined TP_ARRAY_IN_OUT         || \
    defined TP_GIVES                || \
    defined TP_RECEIVES             || \
    defined TP_RECEIVES_GIVES       || \
    defined TP_ARRAY_GIVES          || \
    defined TP_ARRAY_RECEIVES       || \
    defined TP_ARRAY_RECEIVES_GIVES
        #error "Tecplot's parameter life-cycle keywords are in direct conflict with other meanings."
#endif

#if defined INCLUDE_OBJECT_LIFECYCLE_ANNOTATIONS
    #define TP_ACQUIRES             __attribute((gccxml("acquires","in")))
    #define TP_RELEASES             __attribute((gccxml("releases","in")))
    #define TP_OUT                  __attribute((gccxml("out")))
    #define TP_IN_OUT               __attribute((gccxml("in","out")))
    #define TP_ARRAY_OUT            __attribute((gccxml("array","out")))
    #define TP_ARRAY_IN_OUT         __attribute((gccxml("array","in","out")))
    #define TP_GIVES                __attribute((gccxml("gives","out")))
    #define TP_RECEIVES             __attribute((gccxml("receives","in")))
    #define TP_RECEIVES_GIVES       __attribute((gccxml("receives","in","gives","out")))
    #define TP_ARRAY_GIVES          __attribute((gccxml("array","gives","out")))
    #define TP_ARRAY_RECEIVES       __attribute((gccxml("array","receives","in")))
    #define TP_ARRAY_RECEIVES_GIVES __attribute((gccxml("array","receives","in","gives","out")))
#else
    #define TP_ACQUIRES
    #define TP_RELEASES
    #define TP_OUT
    #define TP_IN_OUT
    #define TP_ARRAY_OUT
    #define TP_ARRAY_IN_OUT
    #define TP_GIVES
    #define TP_RECEIVES
    #define TP_RECEIVES_GIVES
    #define TP_ARRAY_GIVES         
    #define TP_ARRAY_RECEIVES      
    #define TP_ARRAY_RECEIVES_GIVES
#endif

/* BEGINREMOVEFROMADDON */
#ifdef NO_ASSERTS /* obfuscate names */
#define ShutDownLicensing FreeAllExtraMapData
#define ProcessYMapInXDirection
#endif /* NO_ASSERTS */


/**************************************
 * LICENSING
 **************************************/
#if defined TECPLOTKERNEL && !defined ENGINE
/* CORE SOURCE CODE REMOVED */
#if defined FLEXLM && defined RLM
#endif
#if !defined FLEXLM && !defined RLM
#endif
#endif

#include "stdafx.h"

#if defined MSWIN && defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#include <string>
#include <map>
#include <vector>
#include <queue>

#include "TranslatedString.h"

/*
 * The following is a temporary fix for figuring out which product is
 * running.  In the future when Focus and 360 use the same code base,
 * we will have to do this dynamically (either with flags on the compiler
 * or variables within Tecplot).
 */
/* ENDREMOVEFROMADDON */

#if defined _WIN32

#if !defined TECPLOTKERNEL

#if !defined MSWIN
#define MSWIN
#endif /* !MSWIN */

/* For the sake of some older add-ons,
   defined _WINDOWS, WINDOWS, and WIN32
   New code should always use MSWIN */

#if !defined WINDOWS
#define WINDOWS
#endif /* WINDOWS */

#if !defined _WINDOWS
#define _WINDOWS
#endif /* !_WINDOWS */

#if !defined WIN32
#define WIN32
#endif /* !WIN32 */

#if defined _DEBUG
#if !defined DEBUG
#define DEBUG
#endif
#elif defined CHECKED_BUILD
#if defined NO_ASSERTS
#undef NO_ASSERTS
#endif
#if !defined NDEBUG
#define NDEBUG
#endif
#else /* RELEASE */
#if !defined NDEBUG
#define NDEBUG
#endif
#if !defined NO_ASSERTS
#define NO_ASSERTS
#endif
#endif /* _DEBUG */
#endif /* TECPLOTKERNEL */

#if _MSC_VER >= 1400
#define VS_2005 /* Using VS2005 Compiler */
#endif

#if !defined TECPLOTKERNEL && defined VS_2005
/* Suppress the warnings about the
     deprecated c runtime functions. */

#if !defined _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif
#endif /* !TECPLOTKERNEL && VS_2005 */

#endif /* MSWIN */

#ifdef NDEBUG
# ifdef _DEBUG
#   error "Both NDEBUG and _DEBUG defined"
# endif
#elif defined TECPLOTKERNEL
# ifndef _DEBUG
#   define _DEBUG
# endif
#endif

/* Now a requirement */
#define USE_3D_HARDWARE

#ifndef THREED
#  define THREED
#endif

#include <stdio.h>
#include <ctype.h>
#include <math.h>

#if defined QUICKDEMO
#define DEMO
#endif

#if defined MicrosoftC
#define DOS
#endif

#if defined CRAYX
#define CRAY
#endif

#if defined IRISX
#define IRIS
#endif

#if defined HPX
#define HPUX
#define HP
#endif

#if defined IBMRS6000X
#define IBMRS6000
#endif

#if defined COMPAQALPHAX
#define COMPAQALPHA
#define COMPAQX
#define COMPAQ
#endif

#if defined DECALPHAX
#define DECALPHA
#define DECX
#endif

#if defined DECX
#define DEC
#endif

#if defined SUNSOLARISX || defined SUNSOLARIS86X
#define SUNX
#endif

#if defined SUNX
#define SUN
#endif

#if defined IRISX || defined CRAYX || defined HPX || defined SUNX || defined CONVEXX
#define UNIXX
#define SYSV
#endif

#if defined DECX || defined LINUX || defined IBMRS6000X || defined COMPAQX || defined DARWIN
#define UNIXX
#endif

/* BEGINREMOVEFROMADDON */
#include <stdarg.h>


/* A bit of OEM stuff */
#define OEM_INVALID_CHECKSUM (LgIndex_t) -1

/* Hide the name of the checksum function */
#if defined NDEBUG
# define DECRYPTTIMEDCODE          FixupPlot
# define CHECKHASHEDCODE           ExpandPlot
# define UPDATECLASSICOEMEHCKSUM   ToggleQuadrants
# define UPDATEOEMCHECKSUM         ComputeAngleFromQuatrant
# define InitOemSettings           InitAngleQuatrantSettings
#endif

#if defined MSWIN
#define USE_TRUETYPEFONTS
#endif
/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */

#ifdef __cplusplus // STL

#ifdef MSWIN

#pragma warning(push, 1) /* warning disabling bellow doesn't actually have any effect on compiler warning.
* It appears that Microsft STL enables all the warning right back on.
* Therefore, the only way to hide them is to push existing warning level,
* lower the level for the time while STL headers are included and then restore
    * previous warning level with a "pragma warning(pop)"
    */

#pragma warning(disable: 4018)  // signed/unsigned mismatch
#pragma warning(disable: 4100)  // unreferenced formal parameter
#pragma warning(disable: 4146)  // unary minus operator applied to unsigned type, 
    // result still unsigned
#pragma warning(disable: 4244)  // 'conversion' conversion from 'type1' to 'type2', 
    // possible loss of data
#pragma warning(disable: 4245)  // conversion from 'type1' to 'type2', signed/unsigned 
    // mismatch
#pragma warning(disable: 4511)  // 'class' : copy constructor could not be generated
#pragma warning(disable: 4512)  // 'class' : assignment operator could not be generated
#pragma warning(disable: 4663)  // C++ language change: to explicitly specialize class 
    // template 'vector'
#pragma warning(disable: 4710)  // 'function' : function not inlined
#pragma warning(disable: 4786)  // identifier was truncated to 'number' characters 
    // in the debug information
#endif

#ifdef MSWIN
#pragma warning(pop) //Restore old warning state.
#endif //MSWIN

#endif //__cplusplus

    /* ENDREMOVEFROMADDON */

#ifdef MSWIN
    /* BEGINREMOVEFROMADDON */
#ifdef TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#ifdef _DEBUG
#endif
#endif /* TECPLOTKERNEL */
    /* ENDREMOVEFROMADDON */

#ifndef TECPLOTKERNEL
#if defined VS_2005
#define Widget LONG_PTR /* correct for 32 & 64 bit builds */
#else
#define Widget long
#endif
#endif



#endif /* MSWIN */


#if defined UNIXX && defined ENGINE
    typedef void *Widget;
#endif


#include <string.h>

#if !defined SYSV && !defined MSWIN
#include <strings.h>
#endif

#if defined (MicrosoftC)
#include <stdlib.h>
#define EXECOS
#ifndef FAR
#define FAR
#endif
#define VOID       void
#endif

#include <sys/types.h>
#include <stdlib.h>

#if defined UNIXX
#if !defined ENGINE
#define X11
#define MOTIF
#endif
#define FAR
#define NEAR
#include <unistd.h>
#endif

/* BEGINREMOVEFROMADDON */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined THREADS_BY_PTHREADS && !defined THREADS_BY_WINAPI
#endif
#if defined THREADS_BY_PTHREADS
#endif
#endif
/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
/* OPENGL currently a must have */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if 0 /* including GLEW header file is currently defining GLAPI to "extern" which causes problems with a kludged Mesa header GLwDrawA.h external glwMDrawingAreaWidgetClass */
    #if defined USE_VBOs
    #endif
#endif
#  if !defined ENGINE
#    if defined UNIXX
#    endif
#  endif
#endif
/* ENDREMOVEFROMADDON */
/*
 * If not building the tecplot kernel then at least
 * include the X Instrinsics.  This will make most
 * development for addons etc work.
 */

/* NOTE: MOTIF not defined if ENGINE is defined */
#if defined MOTIF
#  if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#    if XmVERSION == 1 && XmREVISION == 0
#    endif
#  else
#    include <X11/Intrinsic.h>
#  endif
#endif

#if defined MOTIF
#define CREATE_DIALOG_PARAMS Widget W
typedef Widget ComboBoxWidget_t;
typedef Widget DropDownListWidget_t;
typedef Widget FileDialogWidget_t;
typedef Widget LabelWidget_t;
typedef Widget ListWidget_t;
typedef Widget OptionMenuWidget_t;
typedef Widget PullDownMenuWidget_t;
typedef Widget ScaleWidget_t;
typedef Widget TextFieldWidget_t;
typedef Widget ToggleWidget_t;
typedef Widget ButtonWidget_t;
typedef Widget GridWidget_t;
#endif
#if defined MSWIN
#include <windows.h>
#define CREATE_DIALOG_PARAMS     CWnd *, LaunchDialogMode_e
typedef Widget ComboBoxWidget_t;
typedef Widget DropDownListWidget_t;
typedef Widget FileDialogWidget_t;
typedef Widget LabelWidget_t;
typedef Widget ListWidget_t;
typedef Widget OptionMenuWidget_t;
typedef Widget PullDownMenuWidget_t;
typedef Widget ScaleWidget_t;
typedef Widget TextFieldWidget_t;
typedef Widget ToggleWidget_t;
typedef Widget ButtonWidget_t;
typedef Widget GridWidget_t;
#endif

/* BEGINREMOVEFROMADDON */
#if defined MSWIN && defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined TRACE
#endif
#if defined TRACE0
#endif
#if defined TRACE1
#endif
#if defined TRACE2
#endif
#if defined TRACE3
#endif
#if defined NDEBUG
#else
#endif
#endif /* MSWIN */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */
/* ENDREMOVEFROMADDON */

/* Assume that if TRACE is not defined, then none of the TRACE macros are */
#if !defined (TRACE)
/* TRACE is not used by non-debug builds */
#if defined NDEBUG
#if defined MSWIN
#define TRACE              __noop
#define TRACE0(s)          __noop
#define TRACE1(S,a1)       __noop
#define TRACE2(s,a1,a2)    __noop
#define TRACE3(s,a1,a2,a3) __noop
#else
#define TRACE(str)           ((void)0)
#define TRACE0(str)          ((void)0)
#define TRACE1(str,a1)       ((void)0)
#define TRACE2(str,a1,a2)    ((void)0)
#define TRACE3(str,a1,a2,a3) ((void)0)
#endif /* MSWIN */
#else /* DEBUG */
#if defined MSWIN
/* If the add-on is running in debug mode but does not
 * use MFC, then no TRACE macro is available. Thus, to make tracing available,
 * map TRACE to the win32 OutpuDebugString() function.
 */
# define TRACE(str)           do { OutputDebugStringA(str); } while (0)
# define TRACE1(str,a1)       do { char s[5000]; sprintf(s,str,a1);       OutputDebugStringA(s); } while (0)
# define TRACE2(str,a1,a2)    do { char s[5000]; sprintf(s,str,a1,a2);    OutputDebugStringA(s); } while (0)
# define TRACE3(str,a1,a2,a3) do { char s[5000]; sprintf(s,str,a1,a2,a3); OutputDebugStringA(s); } while (0)
# define TRACE0(str) TRACE(str)
#else
#define TRACE  printf
#define TRACE0 printf
#define TRACE1 printf
#define TRACE2 printf
#define TRACE3 printf
#endif /* MSWIN */
#endif /* NDEBUG */
#endif /* !defined (TRACE) */


/*
  Platform independent way for add-ons to know how much space
  to allocate for a filename.
*/
#if !defined MAX_SIZEOFUTF8CHAR
#define MAX_SIZEOFUTF8CHAR 1
#endif

#if !defined (MaxCharsFilePath)
# if defined (MSWIN)
#   define MaxCharsFilePath (_MAX_PATH*MAX_SIZEOFUTF8CHAR+1) /* Includes traling '\0' */
# else
#   define MaxCharsFilePath 2047 /* ...not really a hard limit for Linux/Unix */
# endif /* MSWIN */
#endif /* !MaxCharsFilePath */

/* BEGINREMOVEFROMADDON */

/*
 * Under Windows, if we are doing a release build (NDEBUG) that is not a CHECKED_BUILD
 * then NO_ASSERTS should be defined
 */
#if defined MSWIN && defined NDEBUG && !defined NO_ASSERTS && !defined CHECKED_BUILD
/* intentionally break the compile */
#  error "define NO_ASSERTS for release builds"
#endif

/*
 * Under Windows, if we are doing a CHECKED_BUILD then it should
 * also be a release build (NDEBUG)
 */
#if defined MSWIN && defined CHECKED_BUILD && !defined NDEBUG
#  error "CHECKED_BUILDS must also be release builds! NDEBUG should be defined but isn't."
#endif


#if defined NO_ASSERTS
#  if !defined USE_MACROS_FOR_FUNCTIONS
#    define USE_MACROS_FOR_FUNCTIONS
#  endif
#endif
/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
/*
 * Under Linux the definition of NULL has a cast that conflicts with our own
 * casting causing warnings that make it tough to find real problems.
 */
#if defined LINUX && defined NULL
# undef NULL
# define NULL 0
#endif

/*
 */
#if !defined MSWIN && !defined ENGINE && !defined ISMESA
#define DISALLOW_OFFSCREEN_EXPORT_IN_BATCH
#endif

/* indentify the platforms capable of using FFMPEG for encoding video formats */
#if defined MSWIN || defined LINUX || defined DARWIN
    #define HAVE_FFMPEG
#endif
/* ENDREMOVEFROMADDON */

/* In windows min and max are being redefined in windef.h.
 * As we want to use the ones provided by the STL we undefined them
 */
#if defined MSWIN && defined max
# undef max
#endif

#if defined MSWIN && defined min
# undef min
#endif

#endif /* _MASTER_H_ */
