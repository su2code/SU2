#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "calc.h"

static int verbose = 0;
static int recurs = 0;
static char buff[1025];
static char errmsg[513];

/*=====================================================================
 * main
 *=====================================================================*/

#define PROMPT "\n$ "

static char *commands[] = {
    "cmd        = execute command 'cmd'",
    "> cmd      = execute command 'cmd' with output",
    "< cmdfile  = read command file 'cmdfile'",
    "& cgnsfile = load CGNS file 'cgnsfile'",
    "%>         = turn on output",
    "%<         = turn off output",
    "%%         = reset calculator",
    "%b base    = set base number 'base'",
    "%z zone    = set zone number 'zone'",
    "%s soln    = set solution number 'soln'",
    "?n         = list number of bases,zones and solutions",
    "?[rcvsi]   = list Reference,Coordinates,Variables,Symbols,Intrinsics",
    "?%         = show base,zone,solution and output state",
    "**         = exit",
    NULL
};

/*---------- next_line ----------------------------------------------
 * return next input line, handling blanks, comments and continuations
 *-------------------------------------------------------------------*/

static char *next_line (FILE *fp)
{
    int n = 0, len;
    char *p, line[257];

    while (fgets (line, sizeof(line), fp) != NULL) {
        line[sizeof(line)-1] = 0;
        p = line + strlen(line);
        while (--p >= line && isspace(*p))
            ;
        *++p = 0;
        for (p = line; *p && isspace(*p); p++)
            ;
        if (!*p) continue;
        strcpy (buff, p);
        n = strlen (buff);
        while (buff[n-1] == '\\') {
            for (n -= 2; n >= 0 && isspace(buff[n]); n--)
                ;
            buff[++n] = 0;
            if (fgets (line, sizeof(line), fp) == NULL) break;
            line[sizeof(line)-1] = 0;
            p = line + strlen(line);
            while (--p >= line && isspace(*p))
                ;
            *++p = 0;
            for (p = line; *p && isspace(*p); p++)
                ;
            if (!*p) break;
            len = strlen (p);
            if (n + len >= sizeof(buff))
                cgnsCalcFatal ("internal command buffer length exceeded");
            strcpy (&buff[n], p);
            n += len;
        }
        if ((p = strchr (buff, '#')) != NULL)
            *p = 0;
        for (p = buff+strlen(buff)-1; p >= buff && isspace(*p); p--)
            ;
        *++p = 0;
        for (p = buff; *p && isspace(*p); p++)
            ;
        if (*p) return (p);
    }
    return (NULL);
}

/*---------- process_command ---------------------------------------
 * extract prefix/suffix and process command
 *------------------------------------------------------------------*/

static void process_command (char *expression)
{
    int echo = verbose;
    char *p, *cmd = expression;
    VECSYM *sym;

    /* skip leading and trailing space */

    for (p = cmd + strlen(cmd) - 1; p >= cmd && isspace(*p); p--)
        ;
    *++p = 0;
    while (*cmd && isspace(*cmd))
        cmd++;
    if (!*cmd) return;

    /* check for echo */

    if (*cmd == '>') {
        echo = 1;
        while (*++cmd && isspace (*cmd))
            ;
    }

    /* empty string */

    if (!*cmd) {
        if (echo) putchar ('\n');
        return;
    }

    /* evaluate the expression */

    if (NULL == (sym = cgnsCalcCommand (cmd)) || !echo)
        return;

    /* value */

    if (vecsym_type(sym) == VECSYM_VALUE) {
        printf ("%s = %g\n", vecsym_name(sym), vecsym_value(sym));
        return;
    }

    /* vector */

    if (vecsym_type(sym) == VECSYM_VECTOR) {
        int n, len = vecsym_veclen(sym);
        VECFLOAT vmin, vmax, *vec = vecsym_vector(sym);
        vmin = vmax = *vec;
        for (n = 1; n < len; n++) {
            if (vmin > vec[n]) vmin = vec[n];
            if (vmax < vec[n]) vmax = vec[n];
        }
        printf ("%s = %d:%g->%g\n", vecsym_name(sym), len, vmin, vmax);
        return;
    }

    /* equation or macro */

    if (vecsym_type(sym) == VECSYM_EQUSTR ||
        vecsym_type(sym) == VECSYM_MACRO) {
        printf ("%s(%d) %s %s\n", vecsym_name(sym), vecsym_nargs(sym),
            vecsym_type(sym) == VECSYM_MACRO ? ":=" : "=", vecsym_equstr(sym));
    }
}

/*---------- parse_commands ---------------------------------------------
 * process command lines
 *-----------------------------------------------------------------------*/

static void parse_commands (FILE *fp, char *pmt)
{
    int n, len;
    char *p;

    if (++recurs == 10)
        cgnsCalcFatal ("too many recursions");

    while (1) {
        if (NULL != pmt && *pmt)
            printf (pmt);
        if ((p = next_line (fp)) == NULL) break;

        /* quit */

        if (*p == '*') {
            if (*++p == '*') {
                cgnsCalcDone ();
                exit (0);
            }
            continue;
        }

        /* set output */

        else if (*p == '%') {
            switch (*++p) {
                case '>':
                    verbose = 1;
                    break;
                case '<':
                    verbose = 0;
                    break;
                case '%':
                    cgnsCalcReset ();
                    break;
                case 'b':
                    n = atoi (++p);
                    if (n < 1 || n > NumBases) {
                        sprintf (errmsg, "base %d invalid\n", n);
                        cgnsCalcError (errmsg);
                    }
                    else
                        cgnsCalcBase (n);
                    break;
                case 'z':
                    n = atoi (++p);
                    if (n < 1 || n > NumZones) {
                        sprintf (errmsg, "zone %d invalid\n", n);
                        cgnsCalcError (errmsg);
                    }
                    else
                        cgnsCalcZone (n);
                    break;
                case 's':
                    n = atoi (++p);
                    if (n < 1 || n > NumSolns) {
                        sprintf (errmsg, "solution %d invalid\n", n);
                        cgnsCalcError (errmsg);
                    }
                    else
                        cgnsCalcSoln (n);
                    break;
            }
            continue;
        }

        /* help */

        if (*p == '?') {
            if (!*++p) {
                for (n = 0; commands[n] != NULL; n++)
                    printf ("%s\n", commands[n]);
                continue;
            }
            while (*p) {
                switch (*p) {
                    case 'n':
                        printf ("number of bases = %d\n", NumBases);
                        printf ("          zones = %d\n", NumZones);
                        printf ("          solns = %d\n", NumSolns);
                        break;
                    case 'r':
                        printf ("=== Reference ===\n");
                        for (n = 0; n < NumReference; n++) {
                            if (reference[n].len > 1)
                                printf ("%s[%d]\n", reference[n].name,
                                    reference[n].len);
                            else
                                printf ("%s\n", reference[n].name);
                        }
                        break;
                    case 'c':
                        printf ("=== Coordinates ===\n");
                        for (n = 0; n < NumCoordinates; n++)
                            printf ("%s[%d]\n", coordinates[n].name,
                                coordinates[n].len);
                        break;
                    case 'v':
                        printf ("=== Solution ===\n");
                        for (n = 0; n < NumVariables; n++)
                            printf ("%s[%d]\n", variables[n].name,
                                variables[n].len);
                        break;
                    case 'i': {
                        char **in = vec_list ();
                        printf ("=== Intrinsics ===\n");
                        for (n = 0; in[n] != NULL; n++)
                            printf ("%s\n", in[n]);
                        break;
                    }
                    case 's':
                        cgnsCalcSymList (stdout);
                        break;
                    case '%':
                        printf ("current base = %d\n", cgnsBase);
                        printf ("        zone = %d\n", cgnsZone);
                        printf ("        soln = %d\n", cgnsSoln);
                        printf ("output is %s\n", verbose ? "on" : "off");
                        break;
                    default:
                        break;
                }
                p++;
            }
            continue;
        }

        /* load CGNS file */

        if (*p == '&') {
            while (*++p && isspace(*p))
                ;
            if (!*p) continue;
            if (access (p, 0)) {
                sprintf (errmsg, "results file <%s> does not exist", p);
                cgnsCalcError (errmsg);
            }
            else
                cgnsCalcInit (p, 0, NULL);
            continue;
        }

        /* load command file */

        if (*p == '<') {
            FILE *file;
            while (*++p && isspace (*p))
                ;
            if (!*p) continue;
            if (access (p, 0)) {
                sprintf (errmsg, "command file <%s> does not exist", p);
                cgnsCalcError (errmsg);
            }
            else if ((file = fopen (p, "r")) == NULL) {
                sprintf (errmsg, "couldn't open command file <%s>", p);
                cgnsCalcError (errmsg);
            }
            else {
                parse_commands (file, NULL);
                fclose (file);
            }
            continue;
        }

        /* process command */

        process_command (p);
    }
    recurs--;
}

/*================================================================*/

int main (int argc, char *argv[])
{
    int n = 1;

    if (n < argc) {
        cgnsCalcInit (argv[n++], 0, NULL);
        if (n < argc) {
            FILE *fp = fopen (argv[n], "r");
            if (NULL == fp) {
                sprintf (errmsg, "couldn't open command file <%s>", argv[n]);
                cgnsCalcFatal (errmsg);
            }
            parse_commands (fp, NULL);
        }
    }
    if (isatty(fileno(stdin)) && isatty(fileno(stdout))) {
        putchar ('\n');
        for (n = 0; commands[n] != NULL; n++)
            printf ("%s\n", commands[n]);
        parse_commands (stdin, PROMPT);
    }
    else
        parse_commands (stdin, NULL);

    cgnsCalcDone ();
    return 0;
}

