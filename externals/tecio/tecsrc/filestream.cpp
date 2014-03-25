#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE

/*
 *****************************************************************
 *****************************************************************
 *******                                                  ********
 ******     Copyright (C) 1988-2010 Tecplot, Inc.          *******
 *******                                                  ********
 *****************************************************************
 *****************************************************************
*/

#define FILESTREAMMODULE

#include "GLOBAL.h"
#include "TASSERT.h"
#include "ALLOC.h"
#include "SYSTEM.h"
#include "FILESTREAM.h"

/**
 */
FileStream_s *FileStreamAlloc(FILE      *File,
                              Boolean_t  IsByteOrderNative)
{
    REQUIRE(VALID_REF(File) || File == NULL);

    FileStream_s *Result = ALLOC_ITEM(FileStream_s, "FileStream");
    if (Result != NULL)
    {
        Result->File              = File;
        Result->IsByteOrderNative = IsByteOrderNative;
    }

    ENSURE(VALID_REF(Result) || Result == NULL);
    return Result;
}

/**
 */
void FileStreamDealloc(FileStream_s **FileStream)
{
    REQUIRE(VALID_REF(FileStream));
    REQUIRE(VALID_REF(*FileStream) || *FileStream == NULL);

    if (*FileStream != NULL)
    {
        FREE_ITEM(*FileStream, "FileStream");
        *FileStream = NULL;
    }

    ENSURE(*FileStream == NULL);
}
