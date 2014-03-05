#ifndef Q_MSG_H
#define Q_MSG_H
/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
#if defined EXTERN
#undef EXTERN
#endif
#if defined Q_MSGMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

#define MAX_STATUS_LINE_MSG_LEN 255

#include "TranslatedString.h"

EXTERN Boolean_t WrapString(const char  *OldString,
                            char       **NewString);
EXTERN void Warning(tecplot::strutil::TranslatedString format,
                    ...); /* zero or more arguments */
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
EXTERN void ErrMsg(tecplot::strutil::TranslatedString format,
                   ...); /* zero or more arguments */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined ENGINE
#endif
#if !defined ENGINE
#if defined MOTIF
#endif
#endif
#if !defined ENGINE
#endif
#if defined Q_MSGMODULE
#else
#endif
#endif // TECPLOTKERNEL

#endif // Q_MSG_H
