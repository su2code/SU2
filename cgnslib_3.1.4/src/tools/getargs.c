#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "getargs.h"

/*---------- usage --------------------------------------------------
 * display usage message and exit
 *-------------------------------------------------------------------*/

void print_usage (char **usgmsg, char *errmsg)
{
    int n;

    if (NULL != errmsg)
        fprintf (stderr, "ERROR: %s\n", errmsg);
    for (n = 0; NULL != usgmsg[n]; n++)
        fprintf (stderr, "%s\n", usgmsg[n]);
    exit (NULL != errmsg);
}

/*---------- getargs ---------------------------------------------------
 * get option letter from argument vector or terminates on error
 * this is similar to getopt()
 *----------------------------------------------------------------------*/

int argind = 0;  /* index into argv array */
int argerr = 1;  /* error output flag */
char *argarg;    /* pointer to argument string */

int getargs (int argc, char **argv, char *ostr)
{
    int argopt;
    char *oli;
    static char *place;
    static int nextarg;

    /* initialization */

    if (!argind)
        nextarg = 1;

    if (nextarg) {  /* update scanning pointer */
        if (argind >= argc || ++argind == argc) {
            argarg = NULL;
            return (-1);
        }
        if ('-' != argv[argind][0]) {
            argarg = argv[argind];
            return (0);
        }
        place = argarg = &argv[argind][1];
        if (!*place) {
            if (++argind == argc) {
                argarg = NULL;
                return (-1);
            }
            argarg = argv[argind];
            return (0);
        }
        nextarg = 0;
    }

    /* check for valid option */

    if ((argopt = *place++) == ':' || argopt == ';' ||
        (oli = strchr (ostr, argopt)) == NULL) {
        if (argerr) {
            fprintf (stderr, "invalid option - `%c'\n", argopt);
            exit (-1);
        }
        return (argopt);
    }

    /* don't need argument */

    if (*++oli != ':') {
        if (*place && *oli == ';') {    /* optional argument */
            argarg = place;
            nextarg = 1;
        }
        else {
            argarg = NULL;
            if (!*place)
                nextarg = 1;
        }
        return (argopt);
    }

    /* get argument */

    if (!*place) {
        if (++argind >= argc) {
            if (!argerr) return (':');
            fprintf (stderr, "missing argument for option `%c'\n", argopt);
            exit (1);
        }
        place = argv[argind];
    }
    argarg = place;
    nextarg = 1;
    return (argopt);
}
