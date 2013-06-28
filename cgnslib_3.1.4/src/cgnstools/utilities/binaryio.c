/*
 * binaryio.c - reads C or FORTRAN binary files
 */

#include <stdio.h>
#include <stdlib.h>     /* for malloc */
#include <string.h>     /* for memcpy, etc */
#include <ctype.h>
#include <math.h>       /* for atof() */
#include "binaryio.h"
#ifdef _WIN32
#include <fcntl.h>      /* for O_BINARY */
#include <io.h>         /* for setmode() */
# define setmode _setmode
# define fileno  _fileno
#endif

/*define STRING_SKIP_SPACE*/

/*----- forward references -----*/

static void swapbytes(unsigned char *, int);

static unsigned char *copy2(unsigned char *);
static unsigned char *copy4(unsigned char *);
static unsigned char *copy8(unsigned char *);

static unsigned char *swap2(unsigned char *);
static unsigned char *swap4(unsigned char *);
static unsigned char *swap8(unsigned char *);

static unsigned char *swap2to4(unsigned char *);
static unsigned char *swap4to2(unsigned char *);

static unsigned char *swap4to8(unsigned char *);
static unsigned char *swap8to4(unsigned char *);

static unsigned char *cray_short(unsigned char *);
static unsigned char *cray_long(unsigned char *);
static unsigned char *cray_float(unsigned char *);
static unsigned char *cray_double(unsigned char *);
static unsigned char *short_cray(unsigned char *);
static unsigned char *long_cray(unsigned char *);
static unsigned char *float_cray(unsigned char *);
static unsigned char *double_cray(unsigned char *);

static unsigned char *convex_float(unsigned char *);
static unsigned char *convex_double(unsigned char *);
static unsigned char *float_convex(unsigned char *);
static unsigned char *double_convex(unsigned char *);

static int rdlongs(BINARYIO *,int,long *);

/*----- error codes and messages -----*/

enum ERRCODES {
    BIOERR_WRITE = 1,
    BIOERR_READ,
    BIOERR_NOEOR,
    BIOERR_BADEOR
};

static char *errmsg[] = {
    "tried to read from file opened for write",
    "tried to write to file opened for read",
    "Fortran EOR mark not found",
    "misplaced Fortran EOR mark"
};

void (*binaryio_error)(char *msg) = NULL;

/*----- macros for local data conversions -----*/

#if ARCH_LOCAL == ARCH_IEEE
# define ieee2short(S)  *((short *)(S))
# define ieee2int(I)    *((int *)(I))
# define ieee2long(L)   *((long *)(L))
# define ieee2float(F)  *((float *)(F))
# define ieee2double(D) *((double *)(D))
# define short2ieee(S)  (unsigned char *)(S)
# define int2ieee(I)    (unsigned char *)(I)
# define long2ieee(L)   (unsigned char *)(L)
# define float2ieee(F)  (unsigned char *)(F)
# define double2ieee(D) (unsigned char *)(D)
#endif

#if ARCH_LOCAL == ARCH_BSIEEE
# define ieee2short(S)  *((short *)swap2((unsigned char *)(S)))
# define ieee2int(I)    *((int *)swap4((unsigned char *)(I)))
# define ieee2long(L)   *((long *)swap4((unsigned char *)(L)))
# define ieee2float(F)  *((float *)swap4((unsigned char *)(F)))
# define ieee2double(D) *((double *)swap8((unsigned char *)(D)))
# define short2ieee(S)  swap2((unsigned char *)(S))
# define int2ieee(I)    swap4((unsigned char *)(I))
# define long2ieee(L)   swap4((unsigned char *)(L))
# define float2ieee(F)  swap4((unsigned char *)(F))
# define double2ieee(D) swap8((unsigned char *)(D))
#endif

#if ARCH_LOCAL == ARCH_CRAY
# define ieee2short(S)  *((short *)short_cray((unsigned char *)(S)))
# define ieee2int(I)    *((int *)long_cray((unsigned char *)(I)))
# define ieee2long(L)   *((long *)long_cray((unsigned char *)(L)))
# define ieee2float(F)  *((float *)float_cray((unsigned char *)(F)))
# define ieee2double(D) *((double *)double_cray((unsigned char *)(D)))
# define short2ieee(S)  cray_short((unsigned char *)(S))
# define int2ieee(I)    cray_long((unsigned char *)(I))
# define long2ieee(L)   cray_long((unsigned char *)(L))
# define float2ieee(F)  cray_float((unsigned char *)(F))
# define double2ieee(D) cray_double((unsigned char *)(D))
#endif

#if ARCH_LOCAL == ARCH_CONVEX
# define ieee2short(S)  *((short *)(S))
# define ieee2int(I)    *((int *)(I))
# define ieee2long(L)   *((long *)(L))
# define ieee2float(F)  *((float *)float_convex((unsigned char *)(F)))
# define ieee2double(D) *((double *)double_convex((unsigned char *)(D)))
# define short2ieee(S)  (unsigned char *)(S)
# define int2ieee(I)    (unsigned char *)(I)
# define long2ieee(L)   (unsigned char *)(L)
# define float2ieee(F)  convex_float((unsigned char *)(F))
# define double2ieee(D) convex_double((unsigned char *)(D))
#endif

#if ARCH_LOCAL == ARCH_ALPHA
# define ieee2short(S)  *((short *)swap2((unsigned char *)(S)))
# define ieee2int(I)    *((int *)swap4((unsigned char *)(I)))
# define ieee2long(L)   *((long *)swap4to8((unsigned char *)(L)))
# define ieee2float(F)  *((float *)swap4((unsigned char *)(F)))
# define ieee2double(D) *((double *)swap8((unsigned char *)(D)))
# define short2ieee(S)  swap2((unsigned char *)(S))
# define int2ieee(I)    swap4((unsigned char *)(I))
# define long2ieee(L)   swap8to4((unsigned char *)(L))
# define float2ieee(F)  swap4((unsigned char *)(F))
# define double2ieee(D) swap8((unsigned char *)(D))
#endif

#if ARCH_LOCAL == ARCH_DOS
# define ieee2short(S)  *((short *)swap2((unsigned char *)(S)))
# define ieee2int(I)    *((int *)swap4to2((unsigned char *)(I)))
# define ieee2long(L)   *((long *)swap4((unsigned char *)(L)))
# define ieee2float(F)  *((float *)swap4((unsigned char *)(F)))
# define ieee2double(D) *((double *)swap8((unsigned char *)(D)))
# define short2ieee(S)  swap2((unsigned char *)(S))
# define int2ieee(I)    swap2to4((unsigned char *)(I))
# define long2ieee(L)   swap4((unsigned char *)(L))
# define float2ieee(F)  swap4((unsigned char *)(F))
# define double2ieee(D) swap8((unsigned char *)(D))
#endif

#ifdef NO_MEMCPY

/*---------- memcpy ------------------------------------------------
 * copys memory buffers
 *------------------------------------------------------------------*/

char *memcpy (dest, s, cnt)
char *dest;
register char *s;
register int cnt;
{
    register char *d = dest;

    if (d != NULL) {
        while (cnt--)
            *d++ = *s++;
    }
    return (dest);
}

#endif  /* NO_MEMCPY */

#ifdef NO_MEMSET

/*---------- memset ------------------------------------------------
 * sets memory to a given character
 *------------------------------------------------------------------*/

char *memset (dest, c, cnt)
char *dest;
register int c;
register int cnt;
{
    register char *d = dest;

    if (d != NULL) {
        while (cnt--)
            *d++ = c;
    }
    return (dest);
}

#endif  /* NO_MEMSET */

/*===================================================================
 * utility routines
 *===================================================================*/

/*---------- set_flags ----------------------------------------------
 * set up I/O flags and conversion routines
 *-------------------------------------------------------------------*/

static int set_flags (BINARYIO *bf, int flags)
{
    /* can't write fortran */

    if ((OPEN_WRITE|OPEN_FORTRAN) == (flags & (OPEN_WRITE|OPEN_FORTRAN)))
        return (0);

    /* initialize */

    if (MACH_DEFAULT == (flags & MACH_UNKNOWN)) {
        switch (flags & ARCH_IEEE) {
            case ARCH_IEEE:
                flags |= MACH_IEEE;
                break;
            case ARCH_BSIEEE:
                flags |= MACH_BSIEEE;
                break;
            case ARCH_CRAY:
                flags |= MACH_CRAY;
                break;
            case ARCH_CONVEX:
                flags |= MACH_CONVEX;
                break;
            case ARCH_ALPHA:
                flags |= MACH_ALPHA;
                break;
            case ARCH_DOS:
                if (OPEN_FORTRAN == (flags & OPEN_FORTRAN))
                    flags |= MACH_DOS32;
                else
                    flags |= MACH_DOS16;
                break;
            default:
                flags |= MACH_LOCAL;
                break;
        }
    }

    bf->short_size  = 2;
    bf->int_size    = 4;
    bf->long_size   = 4;
    bf->float_size  = 4;
    bf->double_size = 8;

    bf->rec_num  = 0L;
    bf->rec_size = 0L;
    bf->rec_read = 0L;

    /* set values based on machine type */

    switch (flags & MACH_UNKNOWN) {
        case MACH_SUN:
        case MACH_IRIS:
        case MACH_HP:
        case MACH_IBM:
        case MACH_IEEE:
        case MACH_UNKNOWN:
            bf->arch       = ARCH_IEEE;
            bf->fromshort  = copy2;
            bf->fromint    = copy4;
            bf->fromlong   = copy4;
            bf->fromfloat  = copy4;
            bf->fromdouble = copy8;
            bf->toshort    = copy2;
            bf->toint      = copy4;
            bf->tolong     = copy4;
            bf->tofloat    = copy4;
            bf->todouble   = copy8;
            break;
        case MACH_DOS32:
        case MACH_DEC:
        case MACH_BSIEEE:
        case MACH_WIN32:
        case MACH_LINUX:
            bf->arch       = ARCH_BSIEEE;
            bf->fromshort  = swap2;
            bf->fromint    = swap4;
            bf->fromlong   = swap4;
            bf->fromfloat  = swap4;
            bf->fromdouble = swap8;
            bf->toshort    = swap2;
            bf->toint      = swap4;
            bf->tolong     = swap4;
            bf->tofloat    = swap4;
            bf->todouble   = swap8;
            break;
        case MACH_CRAY:
            bf->arch        = ARCH_CRAY;
            bf->short_size  = 8;
            bf->int_size    = 8;
            bf->long_size   = 8;
            bf->float_size  = 8;
            bf->double_size = 8;
            bf->fromshort   = cray_short;
            bf->fromint     = cray_long;
            bf->fromlong    = cray_long;
            bf->fromfloat   = cray_float;
            bf->fromdouble  = cray_double;
            bf->toshort     = short_cray;
            bf->toint       = long_cray;
            bf->tolong      = long_cray;
            bf->tofloat     = float_cray;
            bf->todouble    = double_cray;
            break;
        case MACH_CONVEX:
            bf->arch       = ARCH_CONVEX;
            bf->fromshort  = copy2;
            bf->fromint    = copy4;
            bf->fromlong   = copy4;
            bf->fromfloat  = convex_float;
            bf->fromdouble = convex_double;
            bf->toshort    = copy2;
            bf->toint      = copy4;
            bf->tolong     = copy4;
            bf->tofloat    = float_convex;
            bf->todouble   = double_convex;
            break;
        case MACH_ALPHA:
            bf->arch       = ARCH_ALPHA;
            bf->long_size  = 8;
            bf->fromshort  = swap2;
            bf->fromint    = swap4;
            bf->fromlong   = swap8to4;
            bf->fromfloat  = swap4;
            bf->fromdouble = swap8;
            bf->toshort    = swap2;
            bf->toint      = swap4;
            bf->tolong     = swap4to8;
            bf->tofloat    = swap4;
            bf->todouble   = swap8;
            break;
        case MACH_DOS16:
            bf->arch       = ARCH_DOS;
            bf->int_size   = 2;
            bf->fromshort  = swap2;
            bf->fromint    = swap2to4;
            bf->fromlong   = swap4;
            bf->fromfloat  = swap4;
            bf->fromdouble = swap8;
            bf->toshort    = swap2;
            bf->toint      = swap4to2;
            bf->tolong     = swap4;
            bf->tofloat    = swap4;
            bf->todouble   = swap8;
            break;
        default:
            return (0);
    }

    bf->flags = flags & (OPEN_WRITE|OPEN_FORTRAN|OPEN_ASCII|MACH_UNKNOWN);
    return (bf->flags);
}

/*---------- fatal_error --------------------------------------------
 * default fatal error handler
 *-------------------------------------------------------------------*/

static int fatal_error (char *funcname, int errcode)
{
    char msg[81];

    sprintf (msg, "%s: %s", funcname, errmsg[errcode-1]);
    if (NULL == binaryio_error) {
        fprintf (stderr, "%s\n", msg);
        exit (errcode);
    }
    (*binaryio_error) (msg);
    return (-1);
}

/*---------- beg_record ----------------------------------------------
 * reads record size in bytes for FORTRAN data
 *--------------------------------------------------------------------*/

static int beg_record (BINARYIO *bf)
{
    long rec_size;

    switch (bf->flags & MACH_UNKNOWN) {

        /*
         * most machines - integer byte count with
         * begin and end records marks
         */

        case MACH_SUN:
        case MACH_IRIS:
        case MACH_HP:
        case MACH_IBM:
        case MACH_DEC:
        case MACH_ALPHA:
        case MACH_CONVEX:
        case MACH_WIN32:
        case MACH_LINUX:
        case MACH_IEEE:
        case MACH_BSIEEE:
        case MACH_UNKNOWN:
            if (1 != rdlongs (bf, 1, &rec_size))
                return (0);
            bf->rec_size = rec_size;
            break;

        /*
         * Cray FORTRAN - record count is given in words (64-bit)
         * and needs to be masked by 0x1FFF. No end record mark
         */

        case MACH_CRAY:
            if (1 != rdlongs (bf, 1, &rec_size))
                return (0);
            rec_size &= 0x01FF;
            bf->rec_size = rec_size << 3;
            break;

        /*
         * Microsoft FORTRAN - first byte of file is 0x4B
         * each record is 128 bytes or less, with a byte for
         * the count at beginning and end of record
         * It also appears that a record mark of 129 indicates
         * a partial record, and 130 indicates EOF.
        */

        case MACH_DOS16:
        case MACH_DOS32:
            if (0L == bf->rec_num)
                (void) getc (bf->fp);
            rec_size = (long) getc (bf->fp);
            if ((long)EOF == rec_size || 130L == rec_size)
                return (0);
            bf->rec_size = rec_size > 128L ? 128L : rec_size;
            break;

        default:
            return (0);
    }

    bf->rec_num++;
    bf->rec_read = 0L;
    return (1);
}

/*---------- end_record ---------------------------------------------
 * read end of record mark
 *-------------------------------------------------------------------*/

static int end_record (BINARYIO *bf)
{
    long rec_size;

    switch (bf->flags & MACH_UNKNOWN) {

        case MACH_SUN:
        case MACH_IRIS:
        case MACH_HP:
        case MACH_IBM:
        case MACH_DEC:
        case MACH_ALPHA:
        case MACH_CONVEX:
        case MACH_WIN32:
        case MACH_LINUX:
        case MACH_IEEE:
        case MACH_BSIEEE:
        case MACH_UNKNOWN:
            if (1 != rdlongs (bf, 1, &rec_size) ||
                bf->rec_size != rec_size)
                return (0);
            break;

        case MACH_CRAY:
            break;

        case MACH_DOS16:
        case MACH_DOS32:
            rec_size = (long) getc (bf->fp);
            if (rec_size > 128L) rec_size = 128L;
            if (bf->rec_size != rec_size)
                return (0);
            break;

        default:
            return (0);
    }
    bf->rec_size = 0L;
    bf->rec_read = 0L;
    return (1);
}

/*---------- skip_space ---------------------------------------------
 * skip over white space and commas
 *-------------------------------------------------------------------*/

static int skip_space (BINARYIO *bf)
{
    int c;

    while (EOF != (c = getc (bf->fp))) {
        if (!isspace (c) && ',' != c)
            return (EOF == ungetc (c, bf->fp));
        if ('\n' == c) {
            bf->rec_num++;
            bf->rec_read = 0L;
        }
        else
            bf->rec_read++;
    }
    return (EOF);
}

/*---------- read_double --------------------------------------------
 * reads a floating point number from ASCII stream
 *-------------------------------------------------------------------*/

static int read_double (BINARYIO *bf, double *dval)
{
    int c, n = 0, point = 0;
    char str[128];

    *dval = 0.0;
    if (skip_space (bf))
        return (-1);

    /* get to first digit */

    if (EOF == (c = getc (bf->fp)))
        return (EOF);
    if ('-' == c || '+' == c) {
        str[n++] = c;
        bf->rec_read++;
        if (EOF == (c = getc (bf->fp)))
            return (EOF);
    }
    if ('.' == c) {
        str[n++] = c;
        bf->rec_read++;
        point = 1;
        if (EOF == (c = getc (bf->fp)))
            return (EOF);
    }

    /* next character needs to be a digit */

    if (!isdigit (c))
        return (EOF == ungetc (c, bf->fp));

    /* get the number */

    do {
        str[n++] = c;
        bf->rec_read++;
        if (EOF == (c = getc (bf->fp)))
            return (EOF);
        if (point && '.' == c) {
            str[n] = 0;
            *dval = atof (str);
            return (EOF == ungetc (c, bf->fp));
        }
    } while (isdigit (c) || '.' == c);

    /* check for an exponent */

    if (NULL != strchr ("DEde", c)) {
        str[n++] = 'e';
        bf->rec_read++;
        if (EOF == (c = getc (bf->fp)))
            return (EOF);
        if ('-' == c || '+' == c) {
            str[n++] = c;
            bf->rec_read++;
            if (EOF == (c = getc (bf->fp)))
                return (EOF);
        }
        if (isdigit (c)) {
            do {
                str[n++] = c;
                bf->rec_read++;
                if (EOF == (c = getc (bf->fp)))
                    return (EOF);
            } while (isdigit (c));
        }
    }

    /* convert number */

    str[n] = 0;
    *dval = atof (str);
    return (EOF == ungetc (c, bf->fp));
}

/*===================================================================
 * C binary reads
 *===================================================================*/

/*---------- rdshorts -----------------------------------------------
 * read 16-bit integers from a file with data conversions
 *-------------------------------------------------------------------*/

static int rdshorts (BINARYIO *bf, int count, short *data)
{
    int n;

    if (ARCH_LOCAL == bf->arch && sizeof(short) == bf->short_size)
        n = (int)fread (data, sizeof(short), count, bf->fp);
    else {
        unsigned char buf[8];
        for (n = 0; n < count; n++, data++) {
            if (1 != fread (buf, bf->short_size, 1, bf->fp))
                break;
            *data = ieee2short ((*bf->fromshort) (buf));
        }
    }
    return (n);
}

/*---------- rdints -------------------------------------------------
 * read integers from a file with data conversions
 *-------------------------------------------------------------------*/

static int rdints (BINARYIO *bf, int count, int *data)
{
    int n;

    if (ARCH_LOCAL == bf->arch && sizeof(int) == bf->int_size)
        n = (int)fread (data, sizeof(int), count, bf->fp);
    else {
        unsigned char buf[8];
        for (n = 0; n < count; n++, data++) {
            if (1 != fread (buf, bf->int_size, 1, bf->fp))
                break;
            *data = ieee2int ((*bf->fromint) (buf));
        }
    }
    return (n);
}

/*---------- rdINTS -------------------------------------------------
 * read Fortran INTEGER*4 integers from a file with data conversions
 *-------------------------------------------------------------------*/

static int rdINTS (BINARYIO *bf, int count, int *data)
{
    int n;
#if ARCH_LOCAL == ARCH_DOS
    unsigned char buf[4];

    for (n = 0; n < count; n++, data++) {
        if (1 != fread (buf, 4, 1, bf->fp))
            break;
        *data = ieee2int ((*bf->fromint) (buf));
    }
#else
    n = rdints (bf, count, data);
#endif
    return (n);
}

/*---------- rdlongs ------------------------------------------------
 * read 32-bit integers from a file with data conversions
 *-------------------------------------------------------------------*/

static int rdlongs (BINARYIO *bf, int count, long *data)
{
    int n;

    if (ARCH_LOCAL == bf->arch && sizeof(long) == bf->long_size)
        n = (int)fread (data, sizeof(long), count, bf->fp);
    else {
        unsigned char buf[8];
        for (n = 0; n < count; n++, data++) {
            if (1 != fread (buf, bf->long_size, 1, bf->fp))
                break;
            *data = ieee2long ((*bf->fromlong) (buf));
        }
    }
    return (n);
}

/*---------- rdLONGS ------------------------------------------------
 * read Fortran INTEGER*8 integers from a file with data conversions
 *-------------------------------------------------------------------*/

static int rdLONGS (BINARYIO *bf, int count, long *data)
{
    int n;
    unsigned char buf[8];

    for (n = 0; n < count; n++, data++) {
        if (1 != fread (buf, 8, 1, bf->fp))
            break;
        if (ARCH_BSIEEE == bf->arch)
            *data = ieee2long ((*bf->fromlong) (buf));
        else
            *data = ieee2long ((*bf->fromlong) (&buf[4]));
    }
    return (n);
}

/*---------- rdfloats -----------------------------------------------
 * read 32-bit floats from a file with data conversions
 *-------------------------------------------------------------------*/

static int rdfloats (BINARYIO *bf, int count, float *data)
{
    int n;

    if (ARCH_LOCAL == bf->arch && sizeof(float) == bf->float_size)
        n = (int)fread (data, sizeof(float), count, bf->fp);
    else {
        unsigned char buf[8];
        for (n = 0; n < count; n++, data++) {
            if (1 != fread (buf, bf->float_size, 1, bf->fp))
                break;
            *data = ieee2float ((*bf->fromfloat) (buf));
        }
    }
    return (n);
}

/*---------- rddoubles ----------------------------------------------
 * read 64-bit floats from a file with data conversions
 *-------------------------------------------------------------------*/

static int rddoubles (BINARYIO *bf, int count, double *data)
{
    int n;

    if (ARCH_LOCAL == bf->arch && sizeof(double) == bf->double_size)
        n = (int)fread (data, sizeof(double), count, bf->fp);
    else {
        unsigned char buf[8];
        for (n = 0; n < count; n++, data++) {
            if (1 != fread (buf, bf->double_size, 1, bf->fp))
                break;
            *data = ieee2double ((*bf->fromdouble) (buf));
        }
    }
    return (n);
}

/*===================================================================
 * C binary writes
 *===================================================================*/

/*---------- wrtshorts ----------------------------------------------
 * write 16-bit integers to a file with data conversions
 *-------------------------------------------------------------------*/

static int wrtshorts (BINARYIO *bf, int count, short *data)
{
    int n;

    if (ARCH_LOCAL == bf->arch && sizeof(short) == bf->short_size)
        n = (int)fwrite (data, sizeof(short), count, bf->fp);
    else {
        unsigned char *buf;
        for (n = 0; n < count; n++, data++) {
            buf = (*bf->toshort) (short2ieee (data));
            if (1 != fwrite (buf, bf->short_size, 1, bf->fp))
                break;
        }
    }
    return (n);
}

/*---------- wrtints ------------------------------------------------
 * write integers to a file with data conversions
 *-------------------------------------------------------------------*/

static int wrtints (BINARYIO *bf, int count, int *data)
{
    int n;

    if (ARCH_LOCAL == bf->arch && sizeof(int) == bf->int_size)
        n = (int)fwrite (data, sizeof(int), count, bf->fp);
    else {
        unsigned char *buf;
        for (n = 0; n < count; n++, data++) {
            buf = (*bf->toint) (int2ieee (data));
            if (1 != fwrite (buf, bf->int_size, 1, bf->fp))
                break;
        }
    }
    return (n);
}

/*---------- wrtlongs -----------------------------------------------
 * write 32-bit integers to a file with data conversions
 *-------------------------------------------------------------------*/

static int wrtlongs (BINARYIO *bf, int count, long *data)
{
    int n;

    if (ARCH_LOCAL == bf->arch && sizeof(long) == bf->long_size)
        n = (int)fwrite (data, sizeof(long), count, bf->fp);
    else {
        unsigned char *buf;
        for (n = 0; n < count; n++, data++) {
            buf = (*bf->tolong) (long2ieee (data));
            if (1 != fwrite (buf, bf->long_size, 1, bf->fp))
                break;
        }
    }
    return (n);
}

/*---------- wrtfloats ----------------------------------------------
 * write 32-bit floats to a file with data conversions
 *-------------------------------------------------------------------*/

static int wrtfloats (BINARYIO *bf, int count, float *data)
{
    int n;

    if (ARCH_LOCAL == bf->arch && sizeof(float) == bf->float_size)
        n = (int)fwrite (data, sizeof(float), count, bf->fp);
    else {
        unsigned char *buf;
        for (n = 0; n < count; n++, data++) {
            buf = (*bf->tofloat) (float2ieee (data));
            if (1 != fwrite (buf, bf->float_size, 1, bf->fp))
                break;
        }
    }
    return (n);
}

/*---------- wrtdoubles ---------------------------------------------
 * write 64-bit floats to a file with data conversions
 *-------------------------------------------------------------------*/

static int wrtdoubles (BINARYIO *bf, int count, double *data)
{
    int n;

    if (ARCH_LOCAL == bf->arch && sizeof(double) == bf->double_size)
        n = (int)fwrite (data, sizeof(double), count, bf->fp);
    else {
        unsigned char *buf;
        for (n = 0; n < count; n++, data++) {
            buf = (*bf->todouble) (double2ieee (data));
            if (1 != fwrite (buf, bf->double_size, 1, bf->fp))
                break;
        }
    }
    return (n);
}

/*===================================================================
 * file manipulation
 *===================================================================*/

/*---------- bf_new -------------------------------------------------
 * create the BINARYIO structure and set flags
 *-------------------------------------------------------------------*/

BINARYIO *bf_new (FILE *fp, int flags)
{
    BINARYIO *bf;

    if (NULL == (bf = (BINARYIO *) malloc (sizeof(BINARYIO))))
        return (NULL);
    if (!set_flags (bf, flags)) {
        free (bf);
        return (NULL);
    }
    bf->fp = fp;
    bf->did_open = 0;
#if defined(MSDOS) || defined(__MSDOS__) || defined(_WIN32)
    if (OPEN_ASCII != (bf->flags & OPEN_ASCII))
        setmode (fileno (bf->fp), O_BINARY);
#endif
    return (bf);
}

/*---------- bf_open ------------------------------------------------
 * open binary file
 *-------------------------------------------------------------------*/

BINARYIO *bf_open (char *fname, int flags)
{
    BINARYIO *bf;

    if (NULL == fname || !*fname ||
        NULL == (bf = (BINARYIO *) malloc (sizeof(BINARYIO))))
        return (NULL);
    if (!set_flags (bf, flags)) {
        free (bf);
        return (NULL);
    }
#if defined(MSDOS) || defined(__MSDOS__) || defined(_WIN32)
    if (OPEN_ASCII == (bf->flags & OPEN_ASCII))
        bf->fp = fopen (fname, bf->flags & OPEN_WRITE ? "w+" : "r");
    else
        bf->fp = fopen (fname, bf->flags & OPEN_WRITE ? "w+b" : "rb");
#else
    bf->fp = fopen (fname, bf->flags & OPEN_WRITE ? "w+" : "r");
#endif
    if (NULL == bf->fp) {
        free (bf);
        return (NULL);
    }
    bf->did_open = 1;
    return (bf);
}

/*---------- bf_close -----------------------------------------------
 * close the file
 *-------------------------------------------------------------------*/

void bf_close (BINARYIO *bf)
{
    if (bf->did_open)
        fclose (bf->fp);
    free (bf);
}

/*---------- bf_rewind ----------------------------------------------
 * rewind the file
 *-------------------------------------------------------------------*/

void bf_rewind (BINARYIO *bf)
{
    rewind (bf->fp);
    bf->rec_num  = 0L;
    bf->rec_size = 0L;
    bf->rec_read = 0L;
}

/*---------- bf_tell ------------------------------------------------
 * get current file offset
 *-------------------------------------------------------------------*/

long bf_tell (BINARYIO *bf)
{
    long offset = ftell (bf->fp);

    if (0L < offset &&
        OPEN_FORTRAN == (bf->flags & (OPEN_FORTRAN|OPEN_ASCII))) {
        int mach = bf->flags & MACH_UNKNOWN;

        /* Microsoft FORTRAN */

        if (MACH_DOS16 == mach || MACH_DOS32 == mach) {
            offset -= (bf->rec_num << 1) + 1L;
            if (bf->rec_size)
                offset--;
        }

        /* Cray FORTRAN */

        else if (MACH_CRAY == mach)
            offset -= (bf->rec_num << 3);

        /* most everybody else */

        else {
            offset -= (bf->rec_num << 2);
            if (bf->rec_size)
                offset -= 4L;
        }
    }
    return (offset);
}

/*---------- bf_seek ------------------------------------------------
 * set position in file
 *-------------------------------------------------------------------*/

int bf_seek (BINARYIO *bf, long offset)
{
    long cur_off = 0L;

    if (0L > offset)
        return (-1);
    if (0L == offset) {
        bf_rewind (bf);
        return (0);
    }

    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        int c;
        bf_rewind (bf);
        while (EOF != (c = getc (bf->fp))) {
            if ('\n' == c) {
                bf->rec_num++;
                bf->rec_read = 0L;
            }
            else
                bf->rec_read++;
            if (++cur_off == offset)
                return (0);
        }
        return (-1);
    }

    if (OPEN_FORTRAN == (bf->flags & OPEN_FORTRAN)) {
        bf_rewind (bf);
        while (1) {
            if (!beg_record (bf))
                return (-1);
            if (offset < cur_off + bf->rec_size)
                break;
            if (fseek (bf->fp, bf->rec_size, 1))
                return (-1);
            cur_off += bf->rec_size;
            if (!end_record (bf))
                return (-1);
        }
        offset -= cur_off;
        bf->rec_read = offset;
        return (fseek (bf->fp, offset, 1));
    }

    return (fseek (bf->fp, offset, 0));
}

/*---------- bf_unget ------------------------------------------------
 * unget last byte
 *--------------------------------------------------------------------*/

int bf_unget (BINARYIO *bf, int c)
{
    if (EOF == ungetc (c, bf->fp))
        return (0);
    if (bf->rec_read)
        bf->rec_read--;
    return (1);
}

/*---------- bf_nextrec ----------------------------------------------
 * seek to start of next record
 *--------------------------------------------------------------------*/

int bf_nextrec (BINARYIO *bf)
{
    int cnt = 0;

    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        int c;
        while (EOF != (c = getc (bf->fp))) {
            cnt++;
            if ('\n' == c) {
                bf->rec_num++;
                bf->rec_read = 0L;
                return (cnt);
            }
        }
        return (0);
    }
    if (OPEN_FORTRAN == (bf->flags & OPEN_FORTRAN)) {
        if (0L == bf->rec_size)
            beg_record (bf);
        cnt = bf->rec_size - bf->rec_read;
        if (bf->rec_read < bf->rec_size)
            fseek (bf->fp, bf->rec_size - bf->rec_read, 1);
        end_record (bf);
        beg_record (bf);
    }
    return (cnt);
}

/*---------- bf_reclen -----------------------------------------------
 * return the current record len
 *--------------------------------------------------------------------*/

int bf_reclen (BINARYIO *bf)
{
    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        int c;
        long len = bf->rec_read;
        while (EOF != (c = getc (bf->fp))) {
            len++;
            if ('\n' == c)
                break;
        }
        fseek (bf->fp, bf->rec_read - len, 1);
        return ('\n' == c ? (int)(len - 1L) : (int)len);
    }
    if (OPEN_FORTRAN == (bf->flags & OPEN_FORTRAN)) {
        if (0L == bf->rec_size)
            beg_record (bf);
        return ((int)bf->rec_size);
    }
    return (0);
}

/*---------- bf_skipspace --------------------------------------------
 * skips white space and commas for ASCII file
 *--------------------------------------------------------------------*/

int bf_skipspace (BINARYIO *bf)
{
    if (OPEN_ASCII == (bf->flags & OPEN_ASCII))
        return (skip_space (bf));
    return (0);
}

/*=====================================================================
 * information
 *=====================================================================*/

/*---------- bf_machname ----------------------------------------------
 * return machine name for given type
 *---------------------------------------------------------------------*/

static char *machnames[] = {
    "unknown",
    "Sun",
    "Iris",
    "HP",
    "IBM",
    "DEC",
    "Alpha",
    "Cray",
    "Convex",
    "DOS16",
    "DOS32",
    "Win32",
    "Linux",
    "generic BSIEEE",
    "generic IEEE"
};

#define NUM_MACH    (sizeof(machnames)/sizeof(char *))

char *bf_machname (int mach)
{
    if (MACH_DEFAULT == mach)
        mach = MACH_LOCAL;
    if (mach < 0 || mach > NUM_MACH)
        mach = 0;
    return (machnames[mach]);
}

/*---------- bf_archname ---------------------------------------------
 * return the architecture name
 *--------------------------------------------------------------------*/

static char *arch_name[] = {
    "32-bit IEEE",
    "32-bit byteswapped IEEE",
    "64-bit Cray",
    "32-bit native Convex",
    "64-bit DEC Alpha",
    "16-bit DOS",
};

char *bf_archname (int mach)
{
    int arch;

    if (MACH_DEFAULT == mach)
        mach = MACH_LOCAL;
    switch (mach) {
        case MACH_SUN:
        case MACH_IRIS:
        case MACH_HP:
        case MACH_IBM:
        case MACH_IEEE:
        case MACH_UNKNOWN:
        default:
            arch = 0;
            break;
        case MACH_DOS32:
        case MACH_DEC:
        case MACH_WIN32:
        case MACH_LINUX:
        case MACH_BSIEEE:
            arch = 1;
            break;
        case MACH_CRAY:
            arch = 2;
            break;
        case MACH_CONVEX:
            arch = 3;
            break;
        case MACH_ALPHA:
            arch = 4;
            break;
        case MACH_DOS16:
            arch = 5;
            break;
    }
    return (arch_name[arch]);
}

/*===================================================================
 * reads
 *===================================================================*/

/*---------- bf_getbytes --------------------------------------------
 * read bytes from a file
 *-------------------------------------------------------------------*/

int bf_getbytes (BINARYIO *bf, int count, unsigned char *data)
{
    int n = 0;
    static char *fname = "bf_getbytes";

    if (OPEN_WRITE == (bf->flags & OPEN_WRITE))
        return (fatal_error (fname, BIOERR_WRITE));

    /* ASCII read */

    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        n = (int)fread (data, 1, count, bf->fp);
        for (count = 0; count < n; count++, data++) {
            if ('\n' == *data) {
                bf->rec_num++;
                bf->rec_read = 0L;
            }
            else
                bf->rec_read++;
        }
    }

    /* Fortran binary read */

    else if (OPEN_FORTRAN == (bf->flags & OPEN_FORTRAN)) {
        int nread;
        while (n < count) {
            if (0L == bf->rec_size) {
                if (!beg_record (bf))
                    return (n);
                continue;
            }
            if (bf->rec_read == bf->rec_size) {
                if (!end_record (bf))
                    return (fatal_error (fname, BIOERR_NOEOR));
                continue;
            }

            nread = (int)(bf->rec_size - bf->rec_read);
            if (nread > count - n)
                nread = count - n;

            if (nread != fread (data, 1, nread, bf->fp))
                return (n + nread);
            bf->rec_read += (long)nread;
            data += nread;
            n += nread;
        }
    }

    /* C binary read */

    else
        n = (int)fread (data, 1, count, bf->fp);

    return (n);
}

/*---------- bf_getstring -------------------------------------------
 * read bytes from a file up to len or end of record
 *-------------------------------------------------------------------*/

int bf_getstring (BINARYIO *bf, int count, char *data)
{
    int n = 0;
    static char *fname = "bf_getstring";

    if (OPEN_WRITE == (bf->flags & OPEN_WRITE))
        return (fatal_error (fname, BIOERR_WRITE));

    /* ASCII read */

    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        int c;
#ifdef STRING_SKIP_SPACE
        if (skip_space (bf))
            return (n);
#endif
        for (n = 0; n < count; n++) {
            if (EOF == (c = getc (bf->fp)))
                return (EOF);
            if ('\n' == c) {
                bf->rec_num++;
                bf->rec_read = 0;
                return (n);
            }
            bf->rec_read++;
            *data++ = c;
        }
    }

    /* Fortran binary read */

    else if (OPEN_FORTRAN == (bf->flags & OPEN_FORTRAN)) {
        if (0L == bf->rec_size) {
            if (!beg_record (bf))
                return (EOF);
        }
        if (bf->rec_read == bf->rec_size) {
            if (!end_record (bf))
                return (fatal_error (fname, BIOERR_NOEOR));
            if (!beg_record (bf))
                return (EOF);
        }
        n = (int)(bf->rec_size - bf->rec_read);
        if (n > count) n = count;
        n = (int)fread (data, 1, n, bf->fp);
        bf->rec_read += (long)n;
    }

    /* C binary read */

    else {
        int c;
        for (n = 0; n < count; n++) {
            if (EOF == (c = getc (bf->fp)))
                return (EOF);
            if ('\0' == c)
                return (n);
            *data++ = c;
        }
    }

    return (n);
}

/*---------- bf_getshorts -------------------------------------------
 * read 16-bit integers from a file
 *-------------------------------------------------------------------*/

int bf_getshorts (BINARYIO *bf, int count, short *data)
{
    int n = 0, cnt;
    static char *fname = "bf_getshorts";

    if (OPEN_WRITE == (bf->flags & OPEN_WRITE))
        return (fatal_error (fname, BIOERR_WRITE));

    /* ASCII read */

    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        for (n = 0; n < count; n++, data++) {
            if (skip_space(bf) || 1 != fscanf (bf->fp, "%hd%n", data, &cnt))
                break;
            bf->rec_read += (long)cnt;
        }
    }

    /* Fortran binary read of INTEGER*2 */

    else if (OPEN_FORTRAN == (bf->flags & OPEN_FORTRAN)) {
        int nread, shift = (ARCH_CRAY == bf->arch ? 3 : 1);
        long nleft;
        while (n < count) {
            if (0L == bf->rec_size) {
                if (!beg_record (bf))
                    return (n);
                continue;
            }
            if (bf->rec_read == bf->rec_size) {
                if (!end_record (bf))
                    return (fatal_error (fname, BIOERR_NOEOR));
                continue;
            }

            nleft = bf->rec_size - bf->rec_read;
            nread = (int)(nleft >> shift);
            if (nread >= count - n) {
                nread = count - n;
                nleft = (long)nread << shift;
            }
            else if (nleft != ((long)nread << shift))
                return (fatal_error (fname, BIOERR_BADEOR));

            if (nread != rdshorts (bf, nread, data))
                return (n + nread);
            bf->rec_read += nleft;
            data += nread;
            n += nread;
        }
    }

    /* C binary read */

    else
        n = rdshorts (bf, count, data);

    return (n);
}

/*---------- bf_getints ---------------------------------------------
 * read integers from a file
 *-------------------------------------------------------------------*/

int bf_getints (BINARYIO *bf, int count, int *data)
{
    int n = 0, cnt;
    static char *fname = "bf_getints";

    if (OPEN_WRITE == (bf->flags & OPEN_WRITE))
        return (fatal_error (fname, BIOERR_WRITE));

    /* ASCII read */

    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        for (n = 0; n < count; n++, data++) {
            if (skip_space(bf) || 1 != fscanf (bf->fp, "%d%n", data, &cnt))
                break;
            bf->rec_read += (long)cnt;
        }
    }

    /* Fortran binary read of INTEGER*4 */

    else if (OPEN_FORTRAN == (bf->flags & OPEN_FORTRAN)) {
        int nread, shift = (ARCH_CRAY == bf->arch ? 3 : 2);
        long nleft;
        while (n < count) {
            if (0L == bf->rec_size) {
                if (!beg_record (bf))
                    return (n);
                continue;
            }
            if (bf->rec_read == bf->rec_size) {
                if (!end_record (bf))
                    return (fatal_error (fname, BIOERR_NOEOR));
                continue;
            }

            nleft = bf->rec_size - bf->rec_read;
            nread = (int)(nleft >> shift);
            if (nread >= count - n) {
                nread = count - n;
                nleft = (long)nread << shift;
            }
            else if (nleft != ((long)nread << shift))
                return (fatal_error (fname, BIOERR_BADEOR));

            if (nread != rdINTS (bf, nread, data))
                return (n + nread);
            bf->rec_read += nleft;
            data += nread;
            n += nread;
        }
    }

    /* C binary read */

    else
        n = rdints (bf, count, data);

    return (n);
}

/*---------- bf_getlongs --------------------------------------------
 * read 32/64-bit integers from a file
 *-------------------------------------------------------------------*/

int bf_getlongs (BINARYIO *bf, int count, long *data)
{
    int n = 0, cnt;
    static char *fname = "bf_getlongs";

    if (OPEN_WRITE == (bf->flags & OPEN_WRITE))
        return (fatal_error (fname, BIOERR_WRITE));

    /* ASCII read */

    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        for (n = 0; n < count; n++, data++) {
            if (skip_space(bf) || 1 != fscanf (bf->fp, "%ld%n", data, &cnt))
                break;
            bf->rec_read += (long)cnt;
        }
    }

    /* Fortran binary read of INTEGER*8 */

    else if (OPEN_FORTRAN == (bf->flags & OPEN_FORTRAN)) {
        int nread;
        long nleft;
        while (n < count) {
            if (0L == bf->rec_size) {
                if (!beg_record (bf))
                    return (n);
                continue;
            }
            if (bf->rec_read == bf->rec_size) {
                if (!end_record (bf))
                    return (fatal_error (fname, BIOERR_NOEOR));
                continue;
            }

            nleft = bf->rec_size - bf->rec_read;
            nread = (int)(nleft >> 3);
            if (nread >= count - n) {
                nread = count - n;
                nleft = (long)nread << 3;
            }
            else if (nleft != ((long)nread << 3))
                return (fatal_error (fname, BIOERR_BADEOR));

            if (nread != rdLONGS (bf, nread, data))
                return (n + nread);
            bf->rec_read += nleft;
            data += nread;
            n += nread;
        }
    }

    /* C binary read */

    else
        n = rdlongs (bf, count, data);

    return (n);
}

/*---------- bf_getfloats -------------------------------------------
 * read 32-bit floats from a file
 *-------------------------------------------------------------------*/

int bf_getfloats (BINARYIO *bf, int count, float *data)
{
    int n = 0;
    double dval;
    static char *fname = "bf_getfloats";

    if (OPEN_WRITE == (bf->flags & OPEN_WRITE))
        return (fatal_error (fname, BIOERR_WRITE));

    /* ASCII read */

    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        for (n = 0; n < count; n++) {
            if (read_double (bf, &dval))
                break;
            *data++ = (float)dval;
        }
    }

    /* Fortran binary read of REAL*4 */

    else if (OPEN_FORTRAN == (bf->flags & OPEN_FORTRAN)) {
        int nread, shift = (ARCH_CRAY == bf->arch ? 3 : 2);
        long nleft;
        while (n < count) {
            if (0L == bf->rec_size) {
                if (!beg_record (bf))
                    return (n);
                continue;
            }
            if (bf->rec_read == bf->rec_size) {
                if (!end_record (bf))
                    return (fatal_error (fname, BIOERR_NOEOR));
                continue;
            }

            nleft = bf->rec_size - bf->rec_read;
            nread = (int)(nleft >> shift);
            if (nread >= count - n) {
                nread = count - n;
                nleft = (long)nread << shift;
            }
            else if (nleft != ((long)nread << shift))
                return (fatal_error (fname, BIOERR_BADEOR));

            if (nread != rdfloats (bf, nread, data))
                return (n + nread);
            bf->rec_read += nleft;
            data += nread;
            n += nread;
        }
    }

    /* C binary read */

    else
        n = rdfloats (bf, count, data);

    return (n);
}

/*---------- bf_getdoubles ------------------------------------------
 * read 64-bit floats from a file
 *-------------------------------------------------------------------*/

int bf_getdoubles (BINARYIO *bf, int count, double *data)
{
    int n = 0;
    static char *fname = "bf_getdoubles";

    if (OPEN_WRITE == (bf->flags & OPEN_WRITE))
        return (fatal_error (fname, BIOERR_WRITE));

    /* ASCII read */

    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        for (n = 0; n < count; n++, data++) {
            if (read_double (bf, data))
                break;
        }
    }

    /* Fortran binary read of REAL*8 */

    else if (OPEN_FORTRAN == (bf->flags & OPEN_FORTRAN)) {
        int nread;
        long nleft;
        while (n < count) {
            if (0L == bf->rec_size) {
                if (!beg_record (bf))
                    return (n);
                continue;
            }
            if (bf->rec_read == bf->rec_size) {
                if (!end_record (bf))
                    return (fatal_error (fname, BIOERR_NOEOR));
                continue;
            }

            nleft = bf->rec_size - bf->rec_read;
            nread = (int)(nleft >> 3);
            if (nread >= count - n) {
                nread = count - n;
                nleft = (long)nread << 3;
            }
            else if (nleft != ((long)nread << 3))
                return (fatal_error (fname, BIOERR_BADEOR));

            if (nread != rddoubles (bf, nread, data))
                return (n + nread);
            bf->rec_read += nleft;
            data += nread;
            n += nread;
        }
    }

    /* C binary read */

    else
        n = rddoubles (bf, count, data);

    return (n);
}

/*===================================================================
 * writes
 *===================================================================*/

/*---------- bf_putbytes --------------------------------------------
 * write bytes to a file
 *-------------------------------------------------------------------*/

int bf_putbytes (BINARYIO *bf, int count, unsigned char *data)
{
    if (OPEN_WRITE != (bf->flags & OPEN_WRITE))
        return (fatal_error ("bf_putbytes", BIOERR_READ));
    return ((int)fwrite (data, 1, count, bf->fp));
}

/*---------- bf_putshorts -------------------------------------------
 * write 16-bit integers to a file
 *-------------------------------------------------------------------*/

int bf_putshorts (BINARYIO *bf, int count, short *data)
{
    if (OPEN_WRITE != (bf->flags & OPEN_WRITE))
        return (fatal_error ("bf_putshorts", BIOERR_READ));
    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        int n;
        for (n = 0; n < count; n++, data++) {
            if (-1 == fprintf (bf->fp, " %hd", *data))
                break;
        }
        putc ('\n', bf->fp);
        return (n);
    }
    return (wrtshorts (bf, count, data));
}

/*---------- bf_putints ---------------------------------------------
 * write integers to a file
 *-------------------------------------------------------------------*/

int bf_putints (BINARYIO *bf, int count, int *data)
{
    if (OPEN_WRITE != (bf->flags & OPEN_WRITE))
        return (fatal_error ("bf_putints", BIOERR_READ));
    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        int n;
        for (n = 0; n < count; n++, data++) {
            if (-1 == fprintf (bf->fp, " %d", *data))
                break;
        }
        putc ('\n', bf->fp);
        return (n);
    }
    return (wrtints (bf, count, data));
}

/*---------- bf_putlongs --------------------------------------------
 * write 32-bit integers to a file
 *-------------------------------------------------------------------*/

int bf_putlongs (BINARYIO *bf, int count, long *data)
{
    if (OPEN_WRITE != (bf->flags & OPEN_WRITE))
        return (fatal_error ("bf_putlongs", BIOERR_READ));
    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        int n;
        for (n = 0; n < count; n++, data++) {
            if (-1 == fprintf (bf->fp, " %ld", *data))
                break;
        }
        putc ('\n', bf->fp);
        return (n);
    }
    return (wrtlongs (bf, count, data));
}

/*---------- bf_putfloats -------------------------------------------
 * write 32-bit floats to a file
 *-------------------------------------------------------------------*/

int bf_putfloats (BINARYIO *bf, int count, float *data)
{
    if (OPEN_WRITE != (bf->flags & OPEN_WRITE))
        return (fatal_error ("bf_putfloats", BIOERR_READ));
    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        int n;
        for (n = 0; n < count; n++, data++) {
            if (-1 == fprintf (bf->fp, " %#g", *data))
                break;
        }
        putc ('\n', bf->fp);
        return (n);
    }
    return (wrtfloats (bf, count, data));
}

/*---------- bf_putdoubles ------------------------------------------
 * write 64-bit floats to a file
 *-------------------------------------------------------------------*/

int bf_putdoubles (BINARYIO *bf, int count, double *data)
{
    if (OPEN_WRITE != (bf->flags & OPEN_WRITE))
        return (fatal_error ("bf_putdoubles", BIOERR_READ));
    if (OPEN_ASCII == (bf->flags & OPEN_ASCII)) {
        int n;
        for (n = 0; n < count; n++, data++) {
            if (-1 == fprintf (bf->fp, " %#lg", *data))
                break;
        }
        putc ('\n', bf->fp);
        return (n);
    }
    return (wrtdoubles (bf, count, data));
}

/*===================================================================
 * data conversion to and from 32-bit IEEE
 *
 * byte ordering used is MSB first
 * byte structure of floating point data types
 *
 *    short  - 2 bytes
 *    int    - 4 bytes
 *    long   - 4 bytes
 *
 *    float  - 4 bytes  <sign><exponent><mantissa>  exp bias = 127
 *                        1       8         23
 *    double - 8 bytes  <sign><exponent><mantissa>  exp bias = 1023
 *                        1      11         52
 *
 * both floats and doubles are normalized to a value between 1 and 2, thus
 * the high order bit of the mantissa is always 1, and is not stored.
 *===================================================================*/

/*-------------------------------------------------------------------
 * IEEE to IEEE, just do copy
 *-------------------------------------------------------------------*/

static unsigned char *copy2 (unsigned char *bytes)
{
    static unsigned char val[2];

    if (bytes != val)
        memcpy (val, bytes, 2);
    return (val);
}

static unsigned char *copy4 (unsigned char *bytes)
{
    static unsigned char val[4];

    if (bytes != val)
        memcpy (val, bytes, 4);
    return (val);
}

static unsigned char *copy8 (unsigned char *bytes)
{
    static unsigned char val[8];

    if (bytes != val)
        memcpy (val, bytes, 8);
    return (val);
}

/*-------------------------------------------------------------------
 * byte-swapped IEEE
 *-------------------------------------------------------------------*/

static void swapbytes (unsigned char *bytes, int nbytes)
{
    int i = 0, j = nbytes - 1;
    unsigned char tmp;

    while (i < j) {
        tmp = bytes[i];
        bytes[i++] = bytes[j];
        bytes[j--] = tmp;
    }
}

static unsigned char *swap2 (unsigned char *bytes)
{
    static unsigned char val[2];

    if (bytes == val)
        swapbytes (val, 2);
    else {
        val[1] = *bytes++;
        val[0] = *bytes;
    }
    return (val);
}

static unsigned char *swap4 (unsigned char *bytes)
{
    static unsigned char val[4];

    if (bytes == val)
        swapbytes (val, 4);
    else {
        int n;
        for (n = 3; n >= 0; n--)
            val[n] = *bytes++;
    }
    return (val);
}

static unsigned char *swap8 (unsigned char *bytes)
{
    static unsigned char val[8];

    if (bytes == val)
        swapbytes (val, 8);
    else {
        int n;
        for (n = 7; n >= 0; n--)
            val[n] = *bytes++;
    }
    return (val);
}

/*-----------------------------------------------------------------------
 * conversions for 16-bit DOS - ints are 2 bytes
 *-----------------------------------------------------------------------*/

static unsigned char *swap2to4 (unsigned char *bytes)
{
    static unsigned char val[4];

    if (bytes == val)
        swapbytes (val, 4);
    else {
        val[3] = *bytes++;
        val[2] = *bytes;
        if (0 == (*bytes & 0x80))
            val[1] = val[0] = 0;
        else
            val[1] = val[0] = 0xFF;
    }
    return (val);
}

static unsigned char *swap4to2 (unsigned char *bytes)
{
    static unsigned char val[2];

    if (bytes == val)
        swapbytes (val, 2);
    else {
        bytes += 2;
        val[1] = *bytes++;
        val[0] = *bytes;
    }
    return (val);
}

/*-----------------------------------------------------------------------
 * conversions for DEC Alpha - longs are 8 bytes
 *-----------------------------------------------------------------------*/

static unsigned char *swap4to8 (unsigned char *bytes)
{
    static unsigned char val[8];

    if (bytes == val)
        swapbytes (val, 8);
    else {
        int n;
        for (n = 3; n >= 0; n--)
            val[n] = *bytes++;
        if (0 == (val[3] & 0x80)) {
            for (n = 4; n < 8 ; n++)
                val[n] = 0;
        }
        else {
            for (n = 4; n < 8 ; n++)
                val[n] = 0xFF;
        }
    }
    return (val);
}

static unsigned char *swap8to4 (unsigned char *bytes)
{
    static unsigned char val[4];

    if (bytes == val)
        swapbytes (val, 4);
    else {
        int n;
        for (n = 3; n >= 0; n--)
            val[n] = *bytes++;
    }
    return (val);
}

/*-----------------------------------------------------------------------
 * data conversions for Cray
 *
 * all data types on the Cray, except for char, are 64 bit quantities
 *
 * floating point values are stored as:
 *
 *    <sign><exponent><implicit bit><mantissa>
 *      1      15           1          47
 *-----------------------------------------------------------------------*/

/*----- IEEE short -> Cray short -----*/

static unsigned char *short_cray (unsigned char *bytes)
{
    static unsigned char val[8];

    memset (val, *bytes & 0x80 ? 0xFF : 0, 8);
    memcpy (&val[6], bytes, 2);
    return (val);
}

/*----- IEEE long -> Cray long -----*/

static unsigned char *long_cray (unsigned char *bytes)
{
    static unsigned char val[8];

    memset (val, *bytes & 0x80 ? 0xFF : 0, 8);
    memcpy (&val[4], bytes, 4);
    return (val);
}

/*----- IEEE float -> Cray float -----*/

static unsigned char *float_cray (unsigned char *bytes)
{
    unsigned exp;
    static unsigned char val[8];

    memset (val, 0, 8);
    exp = ((unsigned)(bytes[0] & 0x7F) << 1) |
          ((unsigned)(bytes[1] & 0x80) >> 7);
    if (exp) {
        exp += 0x3F82;
        val[0] = (bytes[0] & 0x80) | ((exp >> 8) & 0x7F);
        val[1] = exp & 0xFF;
        val[2] = (bytes[1] & 0x7F) | 0x80;
        memcpy (&val[3], &bytes[2], 2);
    }
    return (val);
}

/*----- IEEE double -> Cray double -----*/

static unsigned char *double_cray (unsigned char *bytes)
{
    unsigned exp;
    static unsigned char val[8];

    memset (val, 0, 8);
    exp = ((unsigned)(bytes[0] & 0x7F) << 4) |
          ((unsigned)(bytes[1] & 0xF0) >> 4);
    if (exp) {
        int n;
        exp += 0x3C02;
        val[0] = (bytes[0] & 0x80) | ((exp >> 8) & 0x7F);
        val[1] = exp & 0xFF;
        val[2] = ((bytes[1] & 0x0F) << 3) | ((bytes[2] & 0xE0) >> 5) | 0x80;
        for (n = 3; n < 8; n++)
            val[n] = ((bytes[n-1] & 0x1F) << 3) | ((bytes[n] & 0xE0) >> 5);
    }
    return (val);
}

/*----- Cray short -> IEEE short -----*/

static unsigned char *cray_short (unsigned char *bytes)
{
    static unsigned char val[2];

    memcpy (val, bytes + 6, 2);
    return (val);
}

/*----- Cray long -> IEEE long -----*/

static unsigned char *cray_long (unsigned char *bytes)
{
    static unsigned char val[4];

    memcpy (val, bytes + 4, 4);
    return (val);
}

/*----- Cray float -> IEEE float -----*/

static unsigned char *cray_float (unsigned char *bytes)
{
    unsigned exp;
    static unsigned char val[4];

    exp = ((unsigned)(bytes[0] & 0x7F) << 8) | (unsigned)bytes[1];
    if (exp < 0x3F82 || exp > 0x4081)
        memset (val, exp < 0x3F82 ? 0 : 0xFF, 4);
    else {
        exp -= 0x3F82;
        val[0] = (bytes[0] & 0x80) | ((exp & 0xFF) >> 1);
        val[1] = (bytes[2] & 0x7F) | ((exp & 0x01) << 7);
        memcpy (&val[2], &bytes[3], 2);
    }
    return (val);
}

/*----- Cray double -> IEEE double -----*/

static unsigned char *cray_double (unsigned char *bytes)
{
    unsigned exp;
    static unsigned char val[8];

    exp = ((unsigned)(bytes[0] & 0x7F) << 8) | (unsigned)bytes[1];
    if (exp < 0x3C02 || exp > 0x4401)
        memset (val, exp < 0x3C02 ? 0 : 0xFF, 8);
    else {
        int n;
        exp -= 0x3C02;
        val[0] = (bytes[0] & 0x80) | ((exp >> 4) & 0x7F);
        val[1] = ((bytes[2] & 0x78) >> 3) | ((exp & 0x0F) << 4);
        for (n = 2; n < 7; n++)
            val[n] = ((bytes[n] & 0x07) << 5) | ((bytes[n+1] & 0xF8) >> 3);
        val[7] = (bytes[7] & 0x07) << 5;
    }
    return (val);
}

/*------------------------------------------------------------------------
 * floating point data conversions for Convex
 *
 * native Convex floating point representation is basically
 * IEEE with an exponent bias of 129 instead of 127
 *------------------------------------------------------------------------*/

/*----- IEEE float -> Convex float -----*/

static unsigned char *float_convex (unsigned char *bytes)
{
    unsigned exp;
    static unsigned char val[4];

    exp = ((unsigned)(bytes[0] & 0x7F) << 8) | (bytes[1] & 0x80);
    if (exp == 0)
        memset (val, 0, 4);
    else if (exp < 0x7E80) {
        exp += 0x0100;
        val[0] = (bytes[0] & 0x80) | (exp >> 8);
        val[1] = (bytes[1] & 0x7F) | (exp & 0x80);
        memcpy (&val[2], &bytes[2], 2);
    }
    else {
        val[0] = (bytes[0] & 0x80) | 0x7F;
        val[1] = 0x80;
        memset (&val[2], 0, 2);
    }
    return (val);
}

/*----- IEEE double -> Convex double -----*/

static unsigned char *double_convex (unsigned char *bytes)
{
    unsigned exp;
    static unsigned char val[8];

    exp = ((unsigned)(bytes[0] & 0x7F) << 8) | (bytes[1] & 0xF0);
    if (exp == 0)
        memset (val, 0, 8);
    else if (exp < 0x7FD0) {
        exp += 0x0020;
        val[0] = (bytes[0] & 0x80) | (exp >> 8);
        val[1] = (bytes[1] & 0x0F) | (exp & 0xF0);
        memcpy (&val[2], &bytes[2], 6);
    }
    else {
        val[0] = (bytes[0] & 0x80) | 0x7F;
        val[1] = 0xF0;
        memset (&val[2], 0, 6);
    }
    return (val);
}

/*----- Convex float -> IEEE float -----*/

static unsigned char *convex_float (unsigned char *bytes)
{
    unsigned exp;
    static unsigned char val[4];

    exp = ((unsigned)(bytes[0] & 0x7F) << 8) | (bytes[1] & 0x80);
    if (exp < 0x0100)
        memset (val, 0, 4);
    else if (exp < 0x7F80) {
        exp -= 0x0100;
        val[0] = (bytes[0] & 0x80) | (exp >> 8);
        val[1] = (bytes[1] & 0x7F) | (exp & 0x80);
        memcpy (&val[2], &bytes[2], 2);
    }
    else {
        val[0] = (bytes[0] & 0x80) | 0x7F;
        val[1] = 0x80;
        memset (&val[2], 0, 2);
    }
    return (val);
}

/*----- Convex double -> IEEE double -----*/

static unsigned char *convex_double (unsigned char *bytes)
{
    unsigned exp;
    static unsigned char val[8];

    exp = ((unsigned)(bytes[0] & 0x7F) << 8) | (bytes[1] & 0xF0);
    if (exp < 0x0200)
        memset (val, 0, 8);
    else if (exp < 0x7FF0) {
        exp -= 0x0020;
        val[0] = (bytes[0] & 0x80) | (exp >> 8);
        val[1] = (bytes[1] & 0x0F) | (exp & 0xF0);
        memcpy (&val[2], &bytes[2], 6);
    }
    else {
        val[0] = (bytes[0] & 0x80) | 0x7F;
        val[1] = 0xF0;
        memset (&val[2], 0, 6);
    }
    return (val);
}
