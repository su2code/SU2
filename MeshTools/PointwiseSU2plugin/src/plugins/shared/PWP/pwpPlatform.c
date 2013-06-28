/****************************************************************************
 *
 * runtime Plugin platform dependent function calls
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#include <string.h>

#if defined(WINDOWS)
//#   include <io.h>
//#   include <windows.h>
#   include <direct.h>
#   include <malloc.h>
#else
//#   include <unistd.h>
//#   include <limits.h>
#   include <stdlib.h>
#endif /* WINDOWS */

#include "pwpPlatform.h"

#if defined(WINDOWS)
#   define sysTextMode    't'
#   define sysBinaryMode  'b'
#else /* *nix or mac */
#   define sysTextMode    '\0'
#   define sysBinaryMode  'b'
#endif /* WINDOWS */


#define strOK(ss)       ((ss) && (ss[0]))
#define chrAPPEND(p,c)  ((p && c) ? (*p++ = c) : 0)

#define isMODE(allBits,modeBits) (((allBits) & (modeBits)) == modeBits)

#if !defined(FILENAME_MAX)
#   define FILENAME_MAX   512
#endif


/*************************************************************/
FILE *
pwpFileOpen(const char *filename, int mode)
{
    FILE *ret = 0;
    while (strOK(filename)) {
        char modeStr[4] = {0,0,0,0};
        char *p = modeStr;
        // get base mode
        if (isMODE(mode,pwpRead)) {
            chrAPPEND(p, 'r');
        }
        else if (isMODE(mode,pwpWrite)) {
            chrAPPEND(p, 'w');
        }
        else if (isMODE(mode,pwpAppend)) {
            chrAPPEND(p, 'a');
        }
        else {
            break; // error
        }

        // get extended mode
        if (isMODE(mode,pwpPlus_)) {
            chrAPPEND(p, '+');
        }

        // get format
        if (isMODE(mode,pwpBinary) || isMODE(mode,pwpUnformatted)) {
            chrAPPEND(p, sysBinaryMode);
        }
        else if (isMODE(mode,pwpFormatted) || isMODE(mode,pwpAscii)) {
            chrAPPEND(p, sysTextMode);
        }
        else {
            break; // error
        }

        ret = fopen(filename, modeStr);
        break; // force exit from while()
    }
    return ret;
}

/*************************************************************/
int
pwpFileClose(FILE *fp)
{
    int ret = sysEOF;
    if (fp) {
        ret = fclose(fp);
    }
    return ret;
}

/*************************************************************/
int
pwpFileEof(FILE *fp)
{
    return fp ? feof(fp) : 0;
}

/*************************************************************/
int
pwpFileFlush(FILE *fp)
{
    return fp ? fflush(fp) : -1;
}

/*************************************************************/
int
pwpFileGetpos(FILE *fp, sysFILEPOS *pos)
{
    int ret = -1;
    if (fp && pos) {
        ret = fgetpos(fp, pos);
    }
    return ret;
}

/*************************************************************/
int
pwpFileSetpos(FILE *fp, const sysFILEPOS *pos)
{
    int ret = -1;
    if (fp && pos) {
        ret = fsetpos(fp, pos);
    }
    return ret;
}

/*************************************************************/
size_t
pwpFileRead(void *buf, size_t size, size_t count, FILE *fp)
{
    size_t ret = 0;
    if (buf && size && count && fp) {
        ret = fread(buf, size, count, fp);
    }
    return ret;
}

/*************************************************************/
size_t
pwpFileWrite(const void *buf, size_t size, size_t count, FILE *fp)
{
    size_t ret = 0;
    if (buf && size && count && fp) {
        ret = fwrite(buf, size, count, fp);
    }
    return ret;
}

/*************************************************************/
size_t
pwpFileWriteStr(const char *str, FILE *fp)
{
    return str ? pwpFileWrite(str, strlen(str), 1, fp) : 0;
}

/*************************************************************/
void
pwpFileRewind(FILE *fp)
{
    if (fp) {
        rewind(fp);
    }
}

/*************************************************************/
int
pwpFileDelete(const char *filename)
{
    int ret = 0;
    if (strOK(filename)) {
        ret = unlink(filename);
    }
    return ret;
}


#define MAXSTACK 50
static char * dirStack[MAXSTACK+1] = { 0 };
static int stackCnt = 0;

/*************************************************************/
int
pwpCwdPush(const char *dir)
{
    int ret = -1;
    if (dir && dir[0] && (stackCnt < MAXSTACK)) {
        char * pCwd = getcwd(0, FILENAME_MAX);
        if (pCwd) {
            if (0 == chdir(dir)) {
                dirStack[stackCnt++] = pCwd;
                ret = 0;
            }
            else {
                free(pCwd);
            }
        }
    }
    return ret;
}

/*************************************************************/
int
pwpCwdPop(void)
{
    int ret = -1;
    if ((0 < stackCnt) && (0 != dirStack[--stackCnt])) {
        ret = chdir(dirStack[stackCnt]);
        free(dirStack[stackCnt]);
        dirStack[stackCnt] = 0;
    }
    return ret;
}
