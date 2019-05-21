 #ifndef _MASTER_H_
 #define _MASTER_H_
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
 #if defined ___1918
 #define TP_ACQUIRES             __attribute((___1546("acquires","in")))
 #define TP_RELEASES             __attribute((___1546("releases","in")))
 #define TP_OUT                  __attribute((___1546("out")))
 #define TP_IN_OUT               __attribute((___1546("in","out")))
 #define TP_ARRAY_OUT            __attribute((___1546("array","out")))
 #define TP_ARRAY_IN_OUT         __attribute((___1546("array","in","out")))
 #define TP_GIVES                __attribute((___1546("gives","out")))
 #define TP_RECEIVES             __attribute((___1546("receives","in")))
 #define TP_RECEIVES_GIVES       __attribute((___1546("receives","in","gives","out")))
 #define TP_ARRAY_GIVES          __attribute((___1546("array","gives","out")))
 #define TP_ARRAY_RECEIVES       __attribute((___1546("array","receives","in")))
 #define TP_ARRAY_RECEIVES_GIVES __attribute((___1546("array","receives","in","gives","out")))
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
 #if defined TP_QUERY
 #error "Tecplot's parameter annotation keywords are in direct conflict with other meanings."
 #endif
 #define TP_QUERY
 #ifdef NO_ASSERTS 
 #define ___3587 ___1529
 #define ___3233
 #endif 
#include "stdafx.h"
#include <string>
#include <map>
#include <vector>
#include <queue>
 #if defined _WIN32
 #if !defined TECPLOTKERNEL
 #if !defined MSWIN
 #define MSWIN
 #endif 
 #if !defined WINDOWS
 #define WINDOWS
 #endif 
 #if !defined _WINDOWS
 #define _WINDOWS
 #endif 
 #if !defined WIN32
 #define WIN32
 #endif 
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
 #else 
 #if !defined NDEBUG
 #define NDEBUG
 #endif
 #if !defined NO_ASSERTS
 #define NO_ASSERTS
 #endif
 #endif 
 #endif 
 #if _MSC_VER >= 1400
 #define ___4444 
 #endif
 #if !defined TECPLOTKERNEL && defined ___4444
 #if !defined _CRT_SECURE_NO_DEPRECATE
 #define _CRT_SECURE_NO_DEPRECATE
 #endif
 #endif 
 #endif 
 #ifdef NDEBUG
 # ifdef _DEBUG
 #   error "Both NDEBUG and _DEBUG defined"
 # endif
 #elif defined TECPLOTKERNEL
 # ifndef _DEBUG
 #   define _DEBUG
 # endif
 #endif
#include "TranslatedString.h"
 #define ___4281
 #ifndef THREED
 #  define THREED
 #endif
#include <stdio.h>
#include <ctype.h>
#include <math.h>
 #if defined ___3260
 #define ___961
 #endif
 #if defined ___2467
 #define ___1100
 #endif
 #if defined CRAYX
 #define CRAY
 #endif
 #if defined ___1995
 #define ___1994
 #endif
 #if defined HPX
 #define HPUX
 #define ___1831
 #endif
 #if defined IBMRS6000X
 #define ___1833
 #endif
 #if defined COMPAQALPHAX
 #define ___534
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
 #if defined ___3892 || defined ___3891
 #define ___3893
 #endif
 #if defined ___3893
 #define ___3886
 #endif
 #if defined ___1995 || defined CRAYX || defined HPX || defined ___3893 || defined ___657
 #define UNIXX
 #define ___3922
 #endif
 #if defined DECX || defined LINUX || defined IBMRS6000X || defined COMPAQX || defined DARWIN
 #define UNIXX
 #endif
#include <stdarg.h>
 #define OEM_INVALID_CHECKSUM (___2227) -1
 #if defined MSWIN
 #define USE_TRUETYPEFONTS
 #endif
 #ifdef MSWIN
 #if defined ___4444
 #define Widget ___2322 
 #else
 #define Widget long
 #endif
 #endif 
 #if defined UNIXX
typedef void *Widget;
 #endif
#include <string.h>
 #if !defined ___3922 && !defined MSWIN
#include <strings.h>
 #endif
 #if defined (___2467)
#include <stdlib.h>
 #define ___1199
 #ifndef ___1306
 #define ___1306
 #endif
 #define VOID       void
 #endif
#include <sys/types.h>
#include <stdlib.h>
 #if defined UNIXX
 #define ___1306
 #define ___2690
#include <unistd.h>
 #endif
 #if defined MSWIN
#include <windows.h>
 #endif
 #if !defined (TRACE)
 #if defined NDEBUG
 #if defined MSWIN
 #define TRACE                       __noop
 #define TRACE0(s)                   __noop
 #define TRACE1(S,a1)                __noop
 #define TRACE2(s,a1,a2)             __noop
 #define TRACE3(s,a1,a2,a3)          __noop
 #define TRACE4(s,a1,a2,a3,a4)       __noop
 #define TRACE5(s,a1,a2,a3,a4,a5)    __noop
 #define TRACE6(s,a1,a2,a3,a4,a5,a6) __noop
 #else
 #define TRACE(str)                      ((void)0)
 #define TRACE0(str)                     ((void)0)
 #define TRACE1(str,a1)                  ((void)0)
 #define TRACE2(str,a1,a2)               ((void)0)
 #define TRACE3(str,a1,a2,a3)            ((void)0)
 #define TRACE4(str,a1,a2,a3,a4)         ((void)0)
 #define TRACE5(str,a1,a2,a3,a4,a5)      ((void)0)
 #define TRACE6(str,a1,a2,a3,a4,a5,a6)   ((void)0)
 #endif 
 #else 
 #if defined MSWIN
 # define TRACE(str)                    do {                                                 OutputDebugStringA(str); }    while (0)
 # define TRACE1(str,a1)                do { char s[5000]; sprintf(s,str,a1);                OutputDebugStringA(s);   }    while (0)
 # define TRACE2(str,a1,a2)             do { char s[5000]; sprintf(s,str,a1,a2);             OutputDebugStringA(s);   }    while (0)
 # define TRACE3(str,a1,a2,a3)          do { char s[5000]; sprintf(s,str,a1,a2,a3);          OutputDebugStringA(s);   }    while (0)
 # define TRACE4(str,a1,a2,a3,a4)       do { char s[5000]; sprintf(s,str,a1,a2,a3,a4);       OutputDebugStringA(s);   }    while (0)
 # define TRACE5(str,a1,a2,a3,a4,a5)    do { char s[5000]; sprintf(s,str,a1,a2,a3,a4,a5);    OutputDebugStringA(s);   }    while (0)
 # define TRACE6(str,a1,a2,a3,a4,a5,a6) do { char s[5000]; sprintf(s,str,a1,a2,a3,a4,a5,a6); OutputDebugStringA(s);   }    while (0)
 # define TRACE0(str) TRACE(str)
 #else
 #define TRACE  printf
 #define TRACE0 printf
 #define TRACE1 printf
 #define TRACE2 printf
 #define TRACE3 printf
 #define TRACE4 printf
 #define TRACE5 printf
 #define TRACE6 printf
 #endif 
 #endif 
 #endif 
 #if !defined MAX_SIZEOFUTF8CHAR
 #define MAX_SIZEOFUTF8CHAR 1
 #endif
 #if !defined (MaxCharsFilePath)
 # if defined (MSWIN)
 #   define MaxCharsFilePath (_MAX_PATH*MAX_SIZEOFUTF8CHAR+1) 
 # else
 #   define MaxCharsFilePath 2047 
 # endif 
 #endif 
 #if defined MSWIN && defined NDEBUG && !defined NO_ASSERTS && !defined CHECKED_BUILD
 #  error "define NO_ASSERTS for release builds"
 #endif
 #if defined MSWIN && defined CHECKED_BUILD && !defined NDEBUG
 #  error "CHECKED_BUILDS must also be release builds! NDEBUG should be defined but isn't."
 #endif
 #if defined NO_ASSERTS
 #  if !defined USE_MACROS_FOR_FUNCTIONS
 #    define USE_MACROS_FOR_FUNCTIONS
 #  endif
 #endif
 #if defined LINUX && defined NULL
 # undef NULL
 # define NULL 0
 #endif
 #if defined MSWIN || defined LINUX || defined DARWIN
 #define ___1823
 #endif
 #if defined __GNUC__ && !defined ___1545
 #define ___1545 (__GNUC__ * 10000 + \
 __GNUC_MINOR__ * 100 + \
 __GNUC_PATCHLEVEL__)
 #endif
 #if defined MSWIN && defined max
 # undef max
 #endif
 #if defined MSWIN && defined min
 # undef min
 #endif
 #endif 
