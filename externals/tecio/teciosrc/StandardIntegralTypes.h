 #ifndef TECPLOT_STANDARDINTEGRALTYPES_H
 #define TECPLOT_STANDARDINTEGRALTYPES_H
 #if !defined __STDC_CONSTANT_MACROS
 #define __STDC_CONSTANT_MACROS
 #endif
 #if !defined __STDC_FORMAT_MACROS
 #define __STDC_FORMAT_MACROS
 #endif
 #if defined _MSC_VER && _MSC_VER < 1800 
 #ifndef _MSC_VER 
 #error "Use this header only with Microsoft Visual C++ compilers!"
 #endif 
 #ifndef _MSC_STDINT_H_ 
 #define _MSC_STDINT_H_
 #if _MSC_VER > 1000
 #pragma once
 #endif
#include <limits.h>
 #ifdef __cplusplus
extern "C" {
 #endif
#  include <wchar.h>
 #ifdef __cplusplus
}
 #endif
 #ifndef _W64
 #  if !defined(__midl) && (defined(_X86_) || defined(_M_IX86)) && _MSC_VER >= 1300
 #     define _W64 __w64
 #  else
 #     define _W64
 #  endif
 #endif
 #if (_MSC_VER < 1300)
typedef signed char       int8_t; typedef signed short      int16_t; typedef signed int        int32_t; typedef unsigned char     uint8_t; typedef unsigned short    uint16_t; typedef unsigned int      uint32_t;
 #else
typedef signed __int8     int8_t; typedef signed __int16    int16_t; typedef signed __int32    int32_t; typedef unsigned __int8   uint8_t; typedef unsigned __int16  uint16_t; typedef unsigned __int32  uint32_t;
 #endif
typedef signed __int64       int64_t; typedef unsigned __int64     uint64_t; typedef int8_t    int_least8_t; typedef int16_t   int_least16_t; typedef int32_t   int_least32_t; typedef int64_t   int_least64_t; typedef uint8_t   uint_least8_t; typedef uint16_t  uint_least16_t; typedef uint32_t  uint_least32_t; typedef uint64_t  uint_least64_t; typedef int8_t    int_fast8_t; typedef int16_t   int_fast16_t; typedef int32_t   int_fast32_t; typedef int64_t   int_fast64_t; typedef uint8_t   uint_fast8_t; typedef uint16_t  uint_fast16_t; typedef uint32_t  uint_fast32_t; typedef uint64_t  uint_fast64_t;
 #ifdef _WIN64 
typedef signed __int64    intptr_t; typedef unsigned __int64  uintptr_t;
 #else 
typedef _W64 signed int   intptr_t; typedef _W64 unsigned int uintptr_t;
 #endif 
typedef int64_t   intmax_t; typedef uint64_t  uintmax_t;
 #if !defined(__cplusplus) || defined(__STDC_LIMIT_MACROS) 
 #define INT8_MIN     ((int8_t)_I8_MIN)
 #define INT8_MAX     _I8_MAX
 #define INT16_MIN    ((int16_t)_I16_MIN)
 #define INT16_MAX    _I16_MAX
 #define INT32_MIN    ((int32_t)_I32_MIN)
 #define INT32_MAX    _I32_MAX
 #define INT64_MIN    ((int64_t)_I64_MIN)
 #define INT64_MAX    _I64_MAX
 #define UINT8_MAX    _UI8_MAX
 #define UINT16_MAX   _UI16_MAX
 #define UINT32_MAX   _UI32_MAX
 #define UINT64_MAX   _UI64_MAX
 #define INT_LEAST8_MIN    INT8_MIN
 #define INT_LEAST8_MAX    INT8_MAX
 #define INT_LEAST16_MIN   INT16_MIN
 #define INT_LEAST16_MAX   INT16_MAX
 #define INT_LEAST32_MIN   INT32_MIN
 #define INT_LEAST32_MAX   INT32_MAX
 #define INT_LEAST64_MIN   INT64_MIN
 #define INT_LEAST64_MAX   INT64_MAX
 #define UINT_LEAST8_MAX   UINT8_MAX
 #define UINT_LEAST16_MAX  UINT16_MAX
 #define UINT_LEAST32_MAX  UINT32_MAX
 #define UINT_LEAST64_MAX  UINT64_MAX
 #define INT_FAST8_MIN    INT8_MIN
 #define INT_FAST8_MAX    INT8_MAX
 #define INT_FAST16_MIN   INT16_MIN
 #define INT_FAST16_MAX   INT16_MAX
 #define INT_FAST32_MIN   INT32_MIN
 #define INT_FAST32_MAX   INT32_MAX
 #define INT_FAST64_MIN   INT64_MIN
 #define INT_FAST64_MAX   INT64_MAX
 #define UINT_FAST8_MAX   UINT8_MAX
 #define UINT_FAST16_MAX  UINT16_MAX
 #define UINT_FAST32_MAX  UINT32_MAX
 #define UINT_FAST64_MAX  UINT64_MAX
 #ifdef _WIN64 
 #  define INTPTR_MIN   INT64_MIN
 #  define INTPTR_MAX   INT64_MAX
 #  define UINTPTR_MAX  UINT64_MAX
 #else 
 #  define INTPTR_MIN   INT32_MIN
 #  define INTPTR_MAX   INT32_MAX
 #  define UINTPTR_MAX  UINT32_MAX
 #endif 
 #define INTMAX_MIN   INT64_MIN
 #define INTMAX_MAX   INT64_MAX
 #define UINTMAX_MAX  UINT64_MAX
 #ifdef _WIN64 
 #  define PTRDIFF_MIN  _I64_MIN
 #  define PTRDIFF_MAX  _I64_MAX
 #else  
 #  define PTRDIFF_MIN  _I32_MIN
 #  define PTRDIFF_MAX  _I32_MAX
 #endif  
 #define SIG_ATOMIC_MIN  INT_MIN
 #define SIG_ATOMIC_MAX  INT_MAX
 #ifndef SIZE_MAX 
 #  ifdef _WIN64 
 #     define SIZE_MAX  _UI64_MAX
 #  else 
 #     define SIZE_MAX  _UI32_MAX
 #  endif 
 #endif 
 #ifndef WCHAR_MIN 
 #  define WCHAR_MIN  0
 #endif  
 #ifndef WCHAR_MAX 
 #  define WCHAR_MAX  _UI16_MAX
 #endif  
 #define WINT_MIN  0
 #define WINT_MAX  _UI16_MAX
 #endif 
 #if !defined(__cplusplus) || defined(__STDC_CONSTANT_MACROS) 
 #define INT8_C(___4298)  ___4298##i8
 #define INT16_C(___4298) ___4298##i16
 #define INT32_C(___4298) ___4298##i32
 #define INT64_C(___4298) ___4298##i64
 #define UINT8_C(___4298)  ___4298##ui8
 #define UINT16_C(___4298) ___4298##ui16
 #define UINT32_C(___4298) ___4298##ui32
 #define UINT64_C(___4298) ___4298##ui64
 #if defined INTMAX_C
 #undef INTMAX_C
 #endif
 #if defined UINTMAX_C
 #undef UINTMAX_C
 #endif
 #define INTMAX_C   INT64_C
 #define UINTMAX_C  UINT64_C
 #endif 
 #endif 
 #define PRId32 "I32d"
 #define PRIu32 "I32u"
 #define PRIx32 "I32x"
 #define PRIu64 "I64u"
 #define PRIx64 "I64x"
 #else 
#include <inttypes.h>
 #endif
typedef int32_t ___1172;
 #ifdef INDEX_16_BIT
typedef  int16_t         ___2227;
 #else
typedef  int64_t         ___2227;
 #endif
static int const InvalidEnumValue = 255;
 #if defined MSWIN
typedef long long          ldfmt_t; typedef unsigned long long lufmt_t;
 #else
typedef long int          ldfmt_t; typedef unsigned long int lufmt_t;
 #endif
 #endif
