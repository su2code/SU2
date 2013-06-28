/*
 * convert_variables.c - convert between primitive and conserved variables
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef _WIN32
# define unlink _unlink
#else
# include <unistd.h>
#endif

#include "getargs.h"
#include "cgnslib.h"
#include "cgnsutil.h"
#include "vec.h"
#include "vecerr.h"
#include "vecsym.h"

#ifndef CG_MODE_MODIFY
# define CG_MODE_MODIFY MODE_MODIFY
#endif

/* command line options */

static char options[] = "pcf:b:z:s:S:g:";

static char *usgmsg[] = {
    "usage  : convert_variables [options] CGNSfile [newCGNSfile]",
    "options:",
    "   -p       = convert to primitive variables",
    "   -c       = convert to conserved variables",
    "   -f<file> = read conversion expressions from <file>",
    "   -b<base> = use CGNS base number <base> (default 1)",
    "   -z<zone> = read zone number <zone> (default all)",
    "   -s<sol>  = read solution number <sol> (default all)",
    "   -S<name> = write to solution <name> (default same as read)",
    "   -g<gamma>= value of gamma for conversions (default 1.4)",
    NULL
};

static double gamma = 1.4;

typedef struct _VAR {
    char name[33];
    int valid;
    VECDATA vd;
} VAR;

static int varlen = 0;
static int maxvars = 0;
static int numvars = 0;
static VAR *vars;

static int numrefs = 0;
static VAR *refs;

static char buff[1024];

/*-------------------------------------------------------------------*/

static void get_error (int errnum, char *errmsg, int pos, char *str)
{
    if (pos >= 0 && NULL != str) {
        fprintf (stderr, "\n%s\n", str);
        while (pos-- > 0)
            putc ('-', stderr);
        putc ('^', stderr);
    }
    FATAL (NULL, errmsg);
}

/*-------------------------------------------------------------------*/

static void get_reference (void)
{
    int n, narrays, na, dim;
    cgsize_t len, vec[12];
    CGNS_ENUMT(DataType_t) datatype;
    char name[33];

    numrefs = 0;
    if (cg_goto (cgnsfn, cgnsbase, "ReferenceState_t", 1, "end") ||
        cg_narrays (&narrays) || narrays < 1)
        return;
    for (na = 1; na <= narrays; na++) {
        if (cg_array_info (na, name, &datatype, &dim, vec))
            FATAL ("get_reference", NULL);
        if (dim >= 1 && vec[0] >= 1) numrefs++;
    }
    if (!numrefs) return;
    refs = (VAR *) malloc (numrefs * sizeof(VAR));
    if (NULL == refs)
        FATAL ("get_reference", "malloc failed for reference variables");

    for (n = 0, na = 1; na <= narrays; na++) {
        if (cg_array_info (na, name, &datatype, &dim, vec))
            FATAL ("get_reference", NULL);
        if (dim < 1 || vec[0] < 1) continue;
        strcpy (refs[n].name, name);
        refs[n].valid = 1;
        len = dim * vec[0];
        if (len == 1) {
            refs[n].vd.type = VEC_VALUE;
            refs[n].vd.len = 0;
            if (cg_array_read_as (na, CGNS_ENUMV(RealDouble), &refs[n].vd.f.val))
                FATAL ("get_reference", NULL);
        }
        else {
            refs[n].vd.type = VEC_VECTOR;
            refs[n].vd.len = len;
            refs[n].vd.f.vec = (VECFLOAT *) malloc (len * sizeof(VECFLOAT));
            if (NULL == refs[n].vd.f.vec)
                FATAL ("get_reference", "malloc failed for reference data");
            if (cg_array_read_as (na, CGNS_ENUMV(RealDouble), refs[n].vd.f.vec))
                FATAL ("get_reference", NULL);
        }
        n++;
    }
}

/*-------------------------------------------------------------------*/

static void init_conversions (ZONE *z, SOLUTION *s)
{
    int n;

    sym_free ();
    sym_addval ("gamma", gamma, NULL);
    read_solution_field (z->id, s->id, 0);
    if (s->nflds > maxvars) {
        n = s->nflds + 10;
        if (maxvars)
            vars = (VAR *) realloc (vars, n * sizeof(VAR));
        else
            vars = (VAR *) malloc (n * sizeof(VAR));
        if (vars == NULL)
            FATAL ("init_conversions", "malloc failed for variables");
        maxvars = n;
    }
    numvars = s->nflds;
    for (n = 0; n < s->nflds; n++) {
        strcpy (vars[n].name, s->flds[n].name);
        vars[n].valid = 1;
        vars[n].vd.type = VEC_VECTOR;
        vars[n].vd.len = s->size;
        vars[n].vd.f.vec = s->flds[n].data;
    }
    varlen = s->size;
}

/*-------------------------------------------------------------------*/

static void update_solution (ZONE *z, SOLUTION *s)
{
    int n, i, nflds = 0;

    for (n = 0; n < numvars; n++) {
        if (vars[n].valid) nflds++;
    }
    s->nflds = nflds;
    if (nflds)
        s->flds = new_field (nflds, 0);
    for (i = 0, n = 0; n < numvars; n++) {
        if (vars[n].valid) {
            strcpy (s->flds[i].name, vars[n].name);
            s->flds[i].data = vars[n].vd.f.vec;
            i++;
        }
        else {
            if (vars[n].vd.f.vec != NULL)
                free (vars[n].vd.f.vec);
        }
    }
}

/*-------------------------------------------------------------------*/

static char *next_token (char *str, char token[33])
{
    int n;
    char *p = str;

    while (*p && (isspace(*p) || *p == ',')) p++;
    if (!*p) return NULL;
    if (*p == '"') {
        p++;
        for (n = 0; n < 32; n++) {
            if (!*p || *p == '"') break;
            token[n] = *p++;
        }
        if (*p++ != '"') FATAL ("missing quote", str);
    }
    else {
        for (n = 0; n < 32; n++) {
            if (!*p || (!isalnum(*p) && *p != '_')) break;
            token[n] = *p++;
        }
    }
    if (!n) return NULL;
    token[n] = 0;
    while (*p && isspace(*p)) p++;
    return p;
}

/*-------------------------------------------------------------------*/

static VECDATA *callback (int check, char **pp, char **err)
{
    int n, nr = 0;
    char name[33], *p = *pp;

    if (*p == '~') {
        nr = 1;
        p++;
    }
    if ((p = next_token (p, name)) == NULL) return NULL;
    if (nr) {
        for (n = 0; n < numrefs; n++) {
            if (0 == strcmp (refs[n].name, name)) {
                *pp = p;
                return &refs[n].vd;
            }
        }
        *err = "reference variable not found";
    }
    else {
        for (n = 0; n < numvars; n++) {
            if (0 == strcmp (vars[n].name, name)) {
                *pp = p;
                return &vars[n].vd;
            }
        }
    }
    return NULL;
}

/*-------------------------------------------------------------------*/

static int delete_var (char *name)
{
    int n;

    for (n = 0; n < numvars; n++) {
        if (0 == strcmp (vars[n].name, name)) {
            vars[n].valid = 0;
            return 1;
        }
    }
    return 0;
}

/*-------------------------------------------------------------------*/

static int add_var (char *name)
{
    int n;
    VECDATA *vd;

    for (n = 0; n < numvars; n++) {
        if (0 == strcmp (vars[n].name, name))
            return 0;
    }
    vd = vec_parse (name, varlen, callback);
    if (vd == NULL || vd->type != VEC_VECTOR || vd->len != varlen)
        FATAL ("error computing", name);
    if (numvars == maxvars) {
        maxvars += 10;
        vars = (VAR *) realloc (vars, maxvars * sizeof(VAR));
        if (vars == NULL)
            FATAL (NULL, "realloc failed for variables");
    }
    strcpy (vars[numvars].name, name);
    vars[numvars].valid = 1;
    vars[numvars].vd.len = vd->len;
    vars[numvars].vd.len = vd->len;
    vars[numvars].vd.f.vec = vd->f.vec;
    numvars++;
    return 1;
}

/*-------------------------------------------------------------------*/

static int parse_command (char *cmd)
{
    int n, nargs = 0;
    char *p, name[33];

    for (p = cmd; *p && isspace(*p); p++)
        ;
    if (!*p) return 0;
    if (*p == '+') {
        strcpy (name, "add");
        p++;
    }
    else if (*p == '-') {
        strcpy (name, "rem");
        p++;
    }
    else {
        p = next_token (p, name);
    }
    if (p == NULL || !*p) return 0;

    if (*p == '(') {
        if (*++p && *p != ')') {
            nargs = atoi (p);
            if (nargs < 0 || nargs > 9)
                get_error (0, "invalid number of arguments",
                    (int)(p - cmd), cmd);
            while (*p && *p != ')') p++;
        }
        if (*p != ')')
            get_error (0, "missing ')'", (int)(p - cmd), cmd);
        while (*++p && isspace(*p));
    }

    if (*p == '=') {
        while (*++p && isspace(*p));
        if (*p)
            sym_addequ (name, nargs, p, NULL);
        return 0;
    }

    if (0 == strncmp (name, "rem", 3) ||
        0 == strncmp (name, "del", 3)) {
        n = 0;
        while ((p = next_token (p, name)) != NULL)
            n += delete_var (name);
        return n;
    }

    if (0 == strncmp (name, "add", 3) ||
        0 == strncmp (name, "sav", 3)) {
        n = 0;
        while ((p = next_token (p, name)) != NULL)
            n += add_var (name);
        return n;
    }

    FATAL ("unknown cmd", name);
    return 0; /* quite compiler */
}

/*-------------------------------------------------------------------*/

char *next_line (FILE *fp)
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
        n = (int)strlen (buff);
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
            len = (int)strlen (p);
            if (n + len >= sizeof(buff))
                FATAL ("next_line", "internal command buffer length exceeded");
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

/*-------------------------------------------------------------------*/

static int parse_file (char *cmdfile)
{
    int changes = 0;
    char *cmd;
    FILE *fp = fopen (cmdfile, "r");

    if (fp == NULL) {
        sprintf (buff, "couldn't open <%s> for reading", cmdfile);
        FATAL ("parse_file", buff);
    }
    while ((cmd = next_line (fp)) != NULL)
        changes += parse_command (cmd);
    fclose (fp);
    return changes;
}

/*-------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
    int n, iz, is, nz, ns, dim;
    int izone = 0, isol = 0;
    char *p, *tmpfile, basename[33];
    char *solname = NULL, *convfile = NULL;
    ZONE *z;
    SOLUTION *s;

    if (argc < 2)
        print_usage (usgmsg, NULL);

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'p':
                convfile = "primitive.cnv";
                break;
            case 'c':
                convfile = "conserved.cnv";
                break;
            case 'f':
                convfile = argarg;
                break;
            case 'b':
                cgnsbase = atoi (argarg);
                break;
            case 'z':
                izone = atoi (argarg);
                break;
            case 's':
                isol = atoi (argarg);
                break;
            case 'S':
                solname = argarg;
                break;
            case 'g':
                gamma = atof (argarg);
                if (gamma <= 1.0)
                    FATAL (NULL, "invalid value for gamma");
                break;
        }
    }

    if (argind == argc)
        print_usage (usgmsg, "CGNSfile not given");
    if (!file_exists (argv[argind]))
        FATAL (NULL, "CGNSfile does not exist or is not a file");

    /* get the conversion expression file */

    if (convfile == NULL)
        print_usage (usgmsg, "need to select one of -p, -c or -f options");
    if ((p = find_file (convfile, argv[0])) == NULL)
        FATAL ("can't locate", convfile);
    convfile = p;
    printf ("converting variables using expressions from %s\n", convfile);

    /* create a working copy */

    printf ("creating a working copy of %s\n", argv[argind]);
    tmpfile = temporary_file (argv[argind]);
    copy_file (argv[argind], tmpfile);

    /* read CGNS file */

    printf ("reading CGNS file from %s\n", tmpfile);
    if (cg_open (tmpfile, CG_MODE_MODIFY, &cgnsfn) ||
        cg_base_read (cgnsfn, cgnsbase, basename, &dim, &dim))
        FATAL (NULL, NULL);
    printf ("reading zone information for base %d - %s\n",
        cgnsbase, basename);

    vec_errhandler = get_error;
    get_reference ();
    read_zones ();

    for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
        iz = 0;
        if (!izone || nz == izone) {
            read_zone_solution (nz);
            for (s = z->sols, ns = 1; ns <= z->nsols; ns++, s++) {
                is = 0;
                if ((!isol || ns == isol) && s->size && s->nflds >= 5) {
                    printf ("checking zone %d, solution %d ... ", nz, ns);
                    fflush (stdout);
                    init_conversions (z, s);
                    if (parse_file (convfile)) {
                        update_solution (z, s);
                        iz = z->id;
                        is = s->id;
                        puts ("updated");
                    }
                    else
                        puts ("ok");
                }
                s->id = is;
            }
        }
        z->id = iz;
    }

    /* write CGNS file */

    for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
        if (z->id) {
            for (s = z->sols, ns = 1; ns <= z->nsols; ns++, s++) {
                if (s->id) {
                    printf ("writing zone %d, solution %d ... ", nz, ns);
                    fflush (stdout);
                    if (solname != NULL) {
                        if (z->nsols == 1)
                            strncpy (s->name, solname, 32);
                        else if (strlen (solname) > 30)
                            sprintf (s->name, "%30.30s%d", solname, s->id);
                        else
                            sprintf (s->name, "%s%d", solname, s->id);
                        s->name[32] = 0;
                    }
                    write_zone_solution (nz, ns);
                    write_solution_field (nz, ns, 0);
                    puts ("done");
                }
            }
        }
    }

    cg_close (cgnsfn);

    if (argind + 1 < argc) argind++;
    printf ("renaming %s to %s\n", tmpfile, argv[argind]);
    unlink (argv[argind]);
    if (rename (tmpfile, argv[argind])) {
        char msg[512];
        sprintf (msg, "rename %s -> %s failed", tmpfile, argv[argind]);
        FATAL (NULL, msg);
        exit (1);
    }
    free (tmpfile);
    return 0;
}
