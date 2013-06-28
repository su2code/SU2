#ifndef TECPLOT_SYSTEM_H
#define TECPLOT_SYSTEM_H

/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2010 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/
#if defined EXTERN
#undef EXTERN
#endif
#if defined SYSTEMMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

/*
 * Windows has trouble with large chunks across SAMBA mounted file systems (not sure
 * where the limit is).
 */
#define MAX_BYTES_PER_CHUNK 131072L // ...== 2^17

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
EXTERN int       OpenFileListGetCount(void);
EXTERN char     *GetLongFileName(const char *FileName);
EXTERN Boolean_t VerifyToOverwriteFile(const char *FName);
EXTERN Boolean_t IsValidDirectory(const char *FName);
EXTERN Boolean_t FileExists(const char *F,
                            Boolean_t   ShowErr);
EXTERN void ErrFName(const char *FName);
EXTERN Boolean_t IsValidFileName(const char *FileName,
                                 Boolean_t   IsReading,
                                 Boolean_t   ShowError);
EXTERN Boolean_t ResizeFile(FILE   *File,
                            Int64_t Length);
EXTERN Boolean_t Close_File(FILE     **F,
                            Boolean_t  ShowErr);
EXTERN Boolean_t Open_File(FILE       **F,
                           const char *FName,
                           Boolean_t  IsReading,
                           Boolean_t  IsAppending,
                           Boolean_t  ForceOpen,
                           Boolean_t  ShowErr,
                           Boolean_t  IsAscii);

#endif
