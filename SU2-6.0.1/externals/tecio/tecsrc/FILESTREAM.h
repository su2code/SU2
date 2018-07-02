/*
 *****************************************************************
 *****************************************************************
 *******                                                  ********
 ******     Copyright (C) 1988-2010 Tecplot, Inc.          *******
 *******                                                  ********
 *****************************************************************
 *****************************************************************
*/
#if !defined FILESTREAM_h
#define FILESTREAM_h

#if defined EXTERN
#  undef EXTERN
#endif
#if defined FILESTREAMMODULE
#  define EXTERN
#else
#  define EXTERN extern
#endif

typedef struct
{
    FILE      *File;
    Boolean_t  IsByteOrderNative;
} FileStream_s;

/**
 * Creates a structure for associating an open file stream with its byte
 * order. The byte order can changed at any time.
 *
 * @param File
 *     Handle to a file which can be NULL.
 * @param IsByteOrderNative
 *     TRUE if the file's byte order is native, FALSE if foreign.
 *
 * @return
 *     An allocated structure associating an open file to its byte order.
 */
EXTERN FileStream_s *FileStreamAlloc(FILE      *File,
                                     Boolean_t  IsByteOrderNative);

/**
 * Deallocates the structure associating the file stream with the byte order.
 * This function does NOT close the file.
 *
 * @param FileStream
 *     Pointer to an open file stream or a pointer to NULL.
 */
EXTERN void FileStreamDealloc(FileStream_s **FileStream);

#endif
