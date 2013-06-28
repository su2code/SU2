/****************************************************************************
 *
 * Pointwise Plugin cross platform functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#ifndef _PWPPLATFORM_H_
#define _PWPPLATFORM_H_

#include <stdio.h>

// defines memory management functions malloc/calloc/free
#if defined(WINDOWS)
#   include <malloc.h>
#else
#   include <stdlib.h>
#endif /* WINDOWS */


#ifdef __cplusplus
extern "C" {
#endif

/*! \file
 \brief Cross Platform Functions

 Certain operations and data-types required by a plugin are platform dependent.
 To simplify access to these operations, the following set of cross-platform
 functions have been implemented in the SDK.

 If a plugin uses these functions and data-types, it should compile and run on
 all Pointwise-supported platforms.
*/


/*-----------------------------------------------------*/
/*! \def sysEOF
* \brief End-of-file value.
*
* \typedef sysFILEPOS
* \brief File position data type.
*/
#if defined(WINDOWS)
#   define sysEOF   EOF
    typedef fpos_t  sysFILEPOS;

#else /* *nix or mac */
#   include <unistd.h>
#   include <stddef.h>
#   define sysEOF   EOF
    typedef fpos_t  sysFILEPOS;

#endif /* WINDOWS */


/*! \brief Bit flags used for pwpFileOpen(mode)
*
* To compose a full mode value, combine one mode flag with one format flag.
*/
typedef enum sysFILEMODE_e {
    // base flags.
    // do not use these directly.
    /*! \cond sdkINTERNALS */
    pwpRead_       = 0x01,
    pwpWrite_      = 0x02,
    pwpAppend_     = 0x04,
    pwpPlus_       = 0x08,

    pwpMaskBase_   = 0x0F, /*!< bits used for base flags */
    pwpMaskFormat_ = 0xF0, /*! bits used for format flags */
    /*! \endcond sdkINTERNALS */

    // combine one of these modes:
    /*! mode: Read access */
    pwpRead        = pwpRead_,
    /*! mode: Write access */
    pwpWrite       = pwpWrite_,
    /*! mode: Append access */
    pwpAppend      = pwpAppend_,
    /*! mode: Read/Write access. File must exist. */
    pwpRWExists    = (pwpRead_ | pwpPlus_),
    /*! mode: Read/Write access. File will be truncated. */
    pwpRWTruncate  = (pwpWrite_ | pwpPlus_),
    /*! mode: Read with Write append access. */
    pwpReadAppend  = (pwpAppend_ | pwpPlus_),

    // with one of these formats:
    /*! format: Generic binary. */
    pwpBinary       = 0x10,
    /*! format: Formatted FORTRAN (ASCII). */
    pwpFormatted    = 0x20,
    /*! format: Unformatted FORTRAN (binary records). */
    pwpUnformatted  = 0x40,
    /*! format: Generic ASCII. */
    pwpAscii        = 0x80
}
sysFILEMODE;


/*---------------------------------------------------------*/
/*! \brief Opens a file for I/O.
*
* \param filename
*    The name of file to open.
*
* \param mode
*    Controls how the file is opened. This value is built from the bitwise
*    OR of 2 sysFILEMODE flags.
*
* \return Pointer to the file or 0 if file could not be opened.
*
* \sa sysFILEMODE, pwpFileClose()
*/
FILE *
pwpFileOpen (const char *filename, int mode);


/*---------------------------------------------------------*/
/*! \brief Closes a file opened with pwpFileOpen().
*
* \param fp
*    A file pointer obtained from pwpFileOpen().
*
* \return 0 on success, sysEOF on error.
*
* \sa pwpFileOpen()
*/
int   
pwpFileClose (FILE *fp);


/*---------------------------------------------------------*/
/*! \brief Queries end-of-file status.
*
* \param fp
*    A file pointer obtained from pwpFileOpen().
*
* \return !0 if fp is at end-of-file.
*/
int   
pwpFileEof (FILE *fp);


/*---------------------------------------------------------*/
/*! \brief Flush a file to disk.
*
* \param fp
*    A file pointer obtained from pwpFileOpen().
*
* \return 0 on success, sysEOF on error.
*/
int   
pwpFileFlush (FILE *fp);


/*---------------------------------------------------------*/
/*! \brief Query the current file position.
*
* \param fp
*    A file pointer obtained from pwpFileOpen().
*
* \param pos
*    Pointer to a sysFILEPOS object.
*
* \return 0 on success, !0 on error.
*
* \sa pwpFileSetpos()
*/
int
pwpFileGetpos (FILE *fp, sysFILEPOS *pos);


/*---------------------------------------------------------*/
/*! \brief Set the current file position.
*
* \param fp
*    A file pointer obtained from pwpFileOpen().
*
* \param pos
*    Pointer to an sysFILEPOS object obtained from a call to pwpFileGetpos().
*
* \return 0 on success, !0 on error.
*
* \sa pwpFileGetpos()
*/
int   
pwpFileSetpos (FILE *fp, const sysFILEPOS *pos);


/*---------------------------------------------------------*/
/*! \brief Read an collection of data items from a file.
*
* \param buf
*    Pointer to the read buffer.
*
* \param size
*    The size of each item.
*
* \param count
*    The number of items to read.
*
* \param fp
*    A file pointer obtained from pwpFileOpen().
*
* \return The number of full items read. Will be < count on an error.
*/
size_t
pwpFileRead (void *buf, size_t size, size_t count, FILE *fp);


/*---------------------------------------------------------*/
/*! \brief Write an collection of data items to a file.
*
* \param buf
*    Pointer to the buffer to write.
*
* \param size
*    The size of each item.
*
* \param count
*    The number of items to write.
*
* \param fp
*    A file pointer obtained from pwpFileOpen().
*
* \return The number of full items written. Will be < count on an error.
*/
size_t
pwpFileWrite (const void *buf, size_t size, size_t count, FILE *fp);


/*---------------------------------------------------------*/
/*! \brief Write a null-terminated string to a file.
*
* \param str
*    The string to write. The null-byte is not written.
*
* \param fp
*    A file pointer obtained from pwpFileOpen().
*
* \return 1 on success, 0 on error.
*/
size_t
pwpFileWriteStr (const char *str, FILE *fp);


/*---------------------------------------------------------*/
/*! \brief Reset position to the beginning of the file.
*
* \param fp
*    A file pointer obtained from pwpFileOpen().
*
* \sa pwpFileSetpos(), pwpFileGetpos()
*/
void  
pwpFileRewind (FILE *fp);


/*---------------------------------------------------------*/
/*! \brief Delete a file.
*
* \param filename
*    The name of an existing file to delete.
*/
int   
pwpFileDelete (const char *filename);


/*---------------------------------------------------------*/
/*! \brief Change the current directory.
*
* Before changing the cwd, the current directory is pushed onto a stack and
* can be restored by calling pwpCwdPop().
*
* \param dir
*    The directory to change to.
*
* \return 0 on success, -1 on error.
*
* \sa pwpCwdPop()
*/
int   
pwpCwdPush (const char *dir);


/*---------------------------------------------------------*/
/*! \brief Restore the current directory.
*
* This call is only valid after calling pwpCwdPush().
*
* \return 0 on success, -1 on error.
*
* \sa pwpCwdPush()
*/
int   
pwpCwdPop (void);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  /* _PWPPLATFORM_H_ */
