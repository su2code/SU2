#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "binaryio.h"

static FILE *p3dout;

void OPENF (int *read, char *fname, int length)
{
    char *p, *buff;

    if (*read) {
        fprintf (stderr, "reading not supported\n");
        exit (1);
    }
#if MACH_LOCAL == MACH_CRAY  || \
    MACH_LOCAL == MACH_DOS16 || \
    MACH_LOCAL == MACH_DOS32 || \
    MACH_LOCAL == MACH_UNKNOWN
    fprintf (stderr,
        "Fortran unformatted output not supported for %s machine\n",
        bf_machname (MACH_LOCAL);
    exit (1);
#else
    buff = (char *) malloc (length + 1);
    if (NULL == buff) {
        fprintf (stderr, "malloc failed for filename working buffer\n");
        exit (1);
    }
    strncpy (buff, fname, length);
    buff[length] = 0;
    for (p = buff+strlen(buff)-1; p >= buff && isspace(*p); p--)
        ;
    *++p = 0;
    for (p = buff; *p && isspace(*p); p++)
        ;
    if (NULL == (p3dout = fopen (p, "w+b"))) {
        fprintf (stderr, "couldn't open <%s> for writing\n", fname);
        exit (1);
    }
    free (buff);
#endif
}

void CLOSEF (void)
{
    fclose (p3dout);
}

void WRITEIF (int *icnt, int *idata, int *ierr)
{
    unsigned int reclen = (unsigned int)*icnt * sizeof(int);

    if (fwrite (&reclen, sizeof(int), 1, p3dout) != 1 ||
        fwrite (idata, sizeof(int), *icnt, p3dout) != *icnt ||
        fwrite (&reclen, sizeof(int), 1, p3dout) != 1)
        *ierr = 1;
    else
        *ierr = 0;
    return;
}

void WRITEFF (int *icnt, float *rdata, int *ierr)
{
    unsigned int reclen = (unsigned int)*icnt * sizeof(float);

    if (fwrite (&reclen, sizeof(int), 1, p3dout) != 1 ||
        fwrite (rdata, sizeof(float), *icnt, p3dout) != *icnt ||
        fwrite (&reclen, sizeof(int), 1, p3dout) != 1)
        *ierr = 1;
    else
        *ierr = 0;
    return;
}

void WRITEDF (int *icnt, double *rdata, int *ierr)
{
    unsigned int reclen = (unsigned int)*icnt * sizeof(double);

    if (fwrite (&reclen, sizeof(int), 1, p3dout) != 1 ||
        fwrite (rdata, sizeof(double), *icnt, p3dout) != *icnt ||
        fwrite (&reclen, sizeof(int), 1, p3dout) != 1)
        *ierr = 1;
    else
        *ierr = 0;
    return;
}

void WRITEGFF(int *icnt, float *rdata, int *idata, int *ierr)
{
    unsigned int rcnt = 3 * (unsigned int)*icnt;
    unsigned int reclen = rcnt * sizeof(float) +
                          (unsigned int)*icnt * sizeof(int);

    if (fwrite (&reclen, sizeof(int), 1, p3dout) != 1 ||
        fwrite (rdata, sizeof(float), rcnt, p3dout) != rcnt ||
        fwrite (idata, sizeof(int), *icnt, p3dout) != *icnt ||
        fwrite (&reclen, sizeof(int), 1, p3dout) != 1)
        *ierr = 1;
    else
        *ierr = 0;
    return;
}

void WRITEGDF(int *icnt, double *rdata, int *idata, int *ierr)
{
    unsigned int rcnt = 3 * (unsigned int)*icnt;
    unsigned int reclen = rcnt * sizeof(double) +
                          (unsigned int)*icnt * sizeof(int);

    if (fwrite (&reclen, sizeof(int), 1, p3dout) != 1 ||
        fwrite (rdata, sizeof(double), rcnt, p3dout) != rcnt ||
        fwrite (idata, sizeof(int), *icnt, p3dout) != *icnt ||
        fwrite (&reclen, sizeof(int), 1, p3dout) != 1)
        *ierr = 1;
    else
        *ierr = 0;
    return;
}

