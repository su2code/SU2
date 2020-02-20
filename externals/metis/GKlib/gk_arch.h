/*!
\file gk_arch.h
\brief This file contains various architecture-specific declerations

\date   Started 3/27/2007
\author George
\version\verbatim $Id: gk_arch.h 10711 2011-08-31 22:23:04Z karypis $ \endverbatim
*/

#ifndef _GK_ARCH_H_
#define _GK_ARCH_H_

/*************************************************************************
* Architecture-specific differences in header files
**************************************************************************/
#ifdef __linux__
#if !defined(__USE_XOPEN)
#define __USE_XOPEN
#endif
#if !defined(_XOPEN_SOURCE)
#define _XOPEN_SOURCE 600
#endif
#if !defined(__USE_XOPEN2K)
#define __USE_XOPEN2K
#endif
#endif


#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif


#ifdef __MSC__
  #include "ms_stdint.h"
  #include "ms_inttypes.h"
  #include "ms_stat.h"
#else
#ifndef SUNOS
  #include <stdint.h>
#endif
  #include <inttypes.h>
  #include <sys/types.h>
#ifndef __MINGW32__
  #include <sys/resource.h>
#endif
  #include <sys/time.h>
#endif


/*************************************************************************
* Architecture-specific modifications
**************************************************************************/
#if defined WIN32 || defined(_MSC_VER)
typedef ptrdiff_t ssize_t;
#endif


#ifdef SUNOS
#define PTRDIFF_MAX  INT64_MAX
#endif

#if (__STDC_VERSION__ >= 199901L)
#define C99
#endif

#if defined __MSC__ && defined C99
/* MSC does not have rint() function */
#define rint(x) ((int)((x)+0.5))

/* MSC does not have INFINITY defined */
#ifndef INFINITY
#define INFINITY FLT_MAX
#endif
#endif

#endif
