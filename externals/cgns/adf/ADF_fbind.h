/*****************************************************************************/
/*                     McDonnell Douglas Corporation                         */
/*****************************************************************************/
/*                                                                           */
/*****************************************************************************/
/*      'typedefs' and macro definitions for interfacing FORTRAN and C       */
/*****************************************************************************/

#if defined(SYSTYPE)
#undef SYSTYPE
#endif

#ifndef USE_ADF_MACROS
/* use the CGNS Fortran macros */
#include "fortran_macros.h"
#define SYSTYPE
typedef	VINTEGER		Fint;
typedef	VFLOAT			Freal;
typedef VREAL			Fdouble;
typedef VCHARACTER		*Fchar;
#define F2CP(FCHAR) STR_PTR(FCHAR)
#define FNAME(A,B) FMNAME(A,B)

#else
#if defined(__convexc__)
/* Convex */
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A ## _

#elif defined(cray)
/* Cray */
#define SYSTYPE
#include <fortran.h>
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef _fcd			Fchar;
#define F2CP(FCHAR) (_fcdtocp(FCHAR))
#define FNAME(A,B) B

#elif defined(__hpux)
/* Hewlett Packard HP-UX */
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A

#elif defined(_AIX)
/* IBM RS/6000 */
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A

#elif defined(__PARAGON__)
/* Intel Paragon */
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A ## _

#elif defined(sgi)
/* Silicon Graphics */
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A ## _

#elif defined(VMS)
/* DEC VAX/VMS */
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef struct dsc$descriptor_s	*Fchar;
#define F2CP(FCHAR) ((FCHAR)->dsc$a_pointer)
#define FNAME(A,B) A
#include <descrip.h>

#elif defined(__alpha)
/* DEC ALPHA OSF/1 */
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A ## _

#elif defined(PPRO)
/* Pentium Pro (P6) using the Intel Reference Compiler */
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A

#elif defined(sun)
/* Sun */
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A ## _

#elif defined(_WIN32)
/* WIN32  */
#define SYSTYPE
typedef	int				Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) B

#elif defined(_CX_UX)
/* Harris Nighthawk */
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A ## _

#elif defined(m88k)
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A/**/_

#elif defined(__linux__)
/* LINUX on Intel */
#define SYSTYPE
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A ## _

#else
typedef	int			Fint;
typedef	float			Freal;
typedef double			Fdouble;
typedef char			*Fchar;
#define F2CP(FCHAR) (FCHAR)
#define FNAME(A,B) A
#endif
#endif
