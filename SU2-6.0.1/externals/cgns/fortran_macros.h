/*-------------------------------------------------------------------------
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from
the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not
   be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution.
-------------------------------------------------------------------------*/

#ifndef	FORTRAN_MACROS_H
#define	FORTRAN_MACROS_H

/**************************************************************/
/* Include file used to interface with FORTRAN */
/**************************************************************/

/* FMNAME is used when you prototype a FORTRAN routine for calling from C   */
/*        (which you HAVE to do), or when you create the subroutine         */
/*	  INTEGER FMNAME(abc,ABC) (INTEGER ival, REAL rval);                */
/*	  void FMNAME(abc,ABC) (INTEGER ival, REAL rval)                    */
/*	  { *ival = 1; *rval = 2.0; return; }                               */
/* FMCALL is used when you call a FORTRAN subroutine from C                 */
/*	  VINTEGER ival;                                                    */
/*	  VREAL rval;                                                       */
/*	  FMCALL(abc,ABC) (&ival, &rval);                                   */

/* STR_PSTR is used in receiving arguments from FORTRAN                     */
/* STR_PLEN is used in end of arguments from FORTRAN (no comma before it)   */
/* STR_PTR  is used to get the address of a string from FORTRAN             */
/* STR_LEN  is used to get the length of a string from FORTRAN              */
/*	    INTEGER FMNAME(abc,ABC) (STR_PSTR(str), INTEGER ival            */
/*				     STR_PLEN(str))                         */
/*	    { char *pointer = STR_PTR(str); int length = STR_LEN(str); }    */

#ifdef NO_CONCATENATION
# define IDENTITY(x)      x
# define CONCATENATE(a,b) IDENTITY(a)b
#else
# define CONCATENATE(a,b) a##b
#endif

/* upper case Fortran name */

#if defined(UPPERCASE)
#	define FMNAME(lname,uname) uname
#	define FMCALL(lname,uname) uname

/* upper case Fortran name with trailing underscore */

#elif defined(UPPERCASE_)
#	define FMNAME(lname,uname) CONCATENATE(uname,_)
#	define FMCALL(lname,uname) CONCATENATE(uname,_)

/* upper case Fortran name with 2 trailing underscores */

#elif defined(UPPERCASE__)
#	define FMNAME(lname,uname) CONCATENATE(uname,__)
#	define FMCALL(lname,uname) CONCATENATE(uname,__)

/* lower case Fortran name */

#elif defined(LOWERCASE)
#	define FMNAME(lname,uname) lname
#	define FMCALL(lname,uname) lname

/* lower case Fortran name with trailing underscore */

#elif defined(LOWERCASE_)
#	define FMNAME(lname,uname) CONCATENATE(lname,_)
#	define FMCALL(lname,uname) CONCATENATE(lname,_)

/* lower case Fortran name with 2 trailing underscores */

#elif defined(LOWERCASE__)
#	define FMNAME(lname,uname) CONCATENATE(lname,__)
#	define FMCALL(lname,uname) CONCATENATE(lname,__)

/* Cray Super Computer */

#elif defined(_CRAY) || defined(cray)
#       include <fortran.h>
#	define FMNAME(lname,uname) uname
#	define FMCALL(lname,uname) uname
#	define STR_PSTR(str)       _fcd str
#	define STR_PLEN(str)
#	define STR_PTR(str)	   _fcdtocp (str)
#	define STR_LEN(str)	   _fcdlen (str)

/* Vax VMS */

#elif defined(VMS)
#       include <descrip.h>
#       define FMNAME(lname,uname) uname
#       define FMCALL(lname,uname) uname
#       define STR_PSTR(str)       struct dsc$descriptor_s *str
#       define STR_PLEN(str)
#       define STR_PTR(str)        str->dsc$a_pointer
#       define STR_LEN(str)        str->dsc$w_length

/* MS Windows */

#elif defined(_WIN32)
# if defined(__NUTC__) || defined(__WIN32_BINARY__)
#	define FMNAME(lname,uname) __cdecl lname
#       define FMCALL(lname,uname) lname
#       define STR_PSTR(str)       char *str
#       define STR_PLEN(str)       , int CONCATENATE(Len,str)
# else
#       ifndef WIN32_FORTRAN
#         define WIN32_FORTRAN
#       endif
#       define FMNAME(lname,uname) __stdcall uname
#       define FMCALL(lname,uname) uname
#       define STR_PSTR(str)       char *str, int CONCATENATE(Len,str)
#       define STR_PLEN(str)
# endif
#       define STR_PTR(str)        str
#       define STR_LEN(str)        CONCATENATE(Len,str)

/* assume lower case Fortran name with trailing underscore */

#else
#	define FMNAME(lname,uname) CONCATENATE(lname,_)
#	define FMCALL(lname,uname) CONCATENATE(lname,_)

#endif

#ifndef STR_PSTR
#	define STR_PSTR(str)       char *str
#	define STR_PLEN(str)       , int CONCATENATE(Len,str)
#	define STR_PTR(str)        str
#	define STR_LEN(str)        CONCATENATE(Len,str)
#endif

/*************/
/* Datatypes */
/*************/

typedef	char		VCHARACTER;
typedef int		VINTEGER;
typedef double		VREAL;
typedef float		VFLOAT;

typedef VCHARACTER	* CHARACTER;
typedef	VINTEGER	* INTEGER;
typedef VREAL		* REAL;
#if !defined(_WIN32) || !defined(_WINDEF_)
typedef VFLOAT		* PFLOAT;
#endif

/**************/
/* Prototypes */
/**************/

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************/
/*  tocstr - convert a fortran character string into a fortran integer array */
/*	     that has a trailing null character to look like a c-string	     */
/*									     */
/*  character	str			(in) Fortran chacter string	     */
/*  integer	icstr(*)		(out) Fortran integer array	     */
/*									     */
/*  notes:								     */
/*    1) Trailing blanks are removed, and leading/trailing '"' are also.     */
/*    2) To keep the trailing blanks, quote them '"'			     */
/*    3) It is a little faster to call this routine without any trailing     */
/*       blanks; ex. call tocstr (abc[1:notblk], icstr)                      */
/*****************************************************************************/

void FMNAME(tocstr,TOCSTR) (
STR_PSTR(str),				/* (in)  Fortran character string */
CHARACTER icstr				/* (out) Fortran integer array */
STR_PLEN(str)				/* (in)  Compiler passed len of str */
			   );

/*****************************************************************************/
/*  frcstr - convert a fortran integer array into a fortran character string */
/*									     */
/*  integer	icstr(*)		(in) Fortran integer array	     */
/*  character	str			(out) Fortran chacter string	     */
/*****************************************************************************/

void FMNAME(frcstr,FRCSTR) (
CHARACTER icstr,			/* (in)  Fortran integer array */
STR_PSTR(str)				/* (out) Fortran character string */
STR_PLEN(str)				/* (in)  Compiler passed len of str */
			   );

/*************************************************************************/
/*  fstr_to_cstr - convert a fortran character string into a c character */
/*		   string						 */
/*************************************************************************/

void fstr_to_cstr (
char *str,				/* (in) Pointer to character string */
int   ilen,				/* (in) Max length of str */
char *icstr				/* (out) C character string */
		  );

/*************************************************************************/
/*  cstr_to_fstr - convert a c character string into a fortran character */
/*		   string						 */
/*************************************************************************/

void cstr_to_fstr (
char *icstr,				/* (in) C character string */
int   ilen,				/* (in) Max length of str */
char *str				/* (out) Pointer to character string */
		  );

#ifdef __cplusplus
}
#endif

#endif	/* FORTRAN_MACROS_H */
