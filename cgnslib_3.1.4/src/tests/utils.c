#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#ifndef CLK_TCK
# define CLK_TCK CLOCKS_PER_SEC
#endif

/* keep the old code in here just in case */
#if 0
#if !defined(_WIN32) || defined(__NUTC__)
#include <sys/times.h>
#endif

double elapsed_time (void)
{
#if defined(_WIN32) && !defined(__NUTC__)
    return (double)clock() / (double)CLK_TCK;
#else
    struct tms clock;
    double usrtime;
    double systime;

    times(&clock);
    usrtime = (double) clock.tms_utime / (double)CLK_TCK;
    systime = (double) clock.tms_stime / (double)CLK_TCK;
    return (usrtime + systime);
#endif
}
#else
double elapsed_time (void)
{
    return (double)clock() / (double)CLK_TCK;
}
#endif

double file_size (char *fname)
{
    struct stat st;

    if (stat (fname, &st)) return 0.0;
    return (double)st.st_size / 1048576.0;
}

