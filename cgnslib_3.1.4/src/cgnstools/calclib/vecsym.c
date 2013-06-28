/*
 * vecsym.c - symbol table for vector processor
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "vecsym.h"

#define HASH_SIZE   127     /* size of hash table */

static VECSYM *hash_table[HASH_SIZE];/* user symbol hash table */

static int num_symbols = 0;         /* number symbols in hash table */
static int max_namelen = 0;         /* max length of user symbol name */

/*----- error mesages -----*/

static char *errmsg[] = {
    "no error",
    "NULL or blank symbol name",
    "NULL or blank equation string",
    "failed to add new symbol to symbol table",
    "malloc failed - out of space",
    "NULL funtion pointer",
    "user function has too many arguments",
    "0 length or NULL vector",
    "symbol name too long",
    "invalid error code"
};

/*----- callback for deleting symbol -----*/

void (*sym_delfunc) (struct symbol_ *) = NULL;

/*----- local functions -----*/

static void init_hash(void);
static int  hash_name(char *name);
static VECSYM *add_symbol(char *name);

/*==================================================================
 * these routines are needed by vec.c
 *==================================================================*/

/*---------- find_symbol -------------------------------------------
 * find a symbol in the hash table
 *------------------------------------------------------------------*/

VECSYM *find_symbol (char *name, int is_func)
{
    VECSYM *sym;

    if (!num_symbols)
        return (NULL);
    sym = hash_table[hash_name (name)];
    while (sym != NULL) {
        if (strcmp (name, vecsym_name(sym)) == 0) {
            if (vecsym_type(sym) == VECSYM_FUNC)
                return (is_func ? sym : NULL);
            if (vecsym_type(sym) == VECSYM_EQUSTR ||
                vecsym_type(sym) == VECSYM_MACRO) {
                if (vecsym_nargs(sym) > 0)
                    return (is_func ? sym : NULL);
                return (sym);
            }
            return (is_func ? NULL : sym);
        }
        sym = sym->next;
    }
    return (NULL);
}

/*==================================================================
 * hashing routines for user symbols
 *==================================================================*/

/*---------- init_hash ---------------------------------------------
 * initialize hash table
 *------------------------------------------------------------------*/

static void init_hash (void)
{
    int n;

    for (n = 0; n < HASH_SIZE; n++)
        hash_table[n] = NULL;
}

/*---------- hash_name ---------------------------------------------
 * hashing routine
 *------------------------------------------------------------------*/

static int hash_name (char *name)
{
    unsigned long hash = 0;

    while (*name)
        hash += *name++;
    return ((int)(hash % HASH_SIZE));
}

/*---------- add_symbol --------------------------------------------
 * add a symbol to the hash table
 *------------------------------------------------------------------*/

static VECSYM *add_symbol (char *name)
{
    int hash = hash_name (name);
    VECSYM *sym = hash_table[hash];

    while (sym != NULL) {
        if (strcmp (name, vecsym_name(sym)) == 0) {
            if (vecsym_type(sym) == VECSYM_EQUSTR)
                free (vecsym_equstr(sym));
            else if (vecsym_type(sym) == VECSYM_MACRO)
                free (vecsym_macro(sym));
            else if (vecsym_type(sym) == VECSYM_VECTOR)
                free (vecsym_vector(sym));
            return (sym);
        }
        sym = sym->next;
    }

    sym = (VECSYM *) malloc (sizeof(VECSYM));
    if (sym != NULL) {
        int len = (int)strlen (name);
        strcpy (vecsym_name(sym), name);
        sym->next = hash_table[hash];
        hash_table[hash] = sym;
        num_symbols++;
        if (max_namelen < len)
            max_namelen = len;
    }
    return (sym);
}

/*---------- sym_free ----------------------------------------------
 * free symbol table resources
 *------------------------------------------------------------------*/

void sym_free (void)
{
    int n;
    VECSYM *sym, *next;

    if (!num_symbols)
        return;
    for (n = 0; n < HASH_SIZE; n++) {
        sym = hash_table[n];
        while (sym != NULL) {
            next = sym->next;
            if (sym_delfunc != NULL) (*sym_delfunc) (sym);
            if (vecsym_type(sym) == VECSYM_EQUSTR)
                free (vecsym_equstr(sym));
            else if (vecsym_type(sym) == VECSYM_MACRO)
                free (vecsym_macro(sym));
            else if (vecsym_type(sym) == VECSYM_VECTOR)
                free (vecsym_vector(sym));
            free (sym);
            sym = next;
        }
        hash_table[n] = NULL;
    }
    num_symbols = 0;
    max_namelen = 0;
}

/*---------- sym_errmsg ---------------------------------------------
 * returns error message
 *-------------------------------------------------------------------*/

char *sym_errmsg (int errnum)
{
    if (errnum >= 0 && errnum < SYMERR_INVALID)
        return (errmsg[errnum]);
    return (errmsg[SYMERR_INVALID]);
}

/*---------- sym_addval --------------------------------------------
 * adds a value to user symbol table
 *------------------------------------------------------------------*/

int sym_addval (char *name, VECFLOAT val, void *user)
{
    VECSYM *sym;

    if (name == NULL || !*name)
        return (SYMERR_NONAME);
    if (strlen (name) > SYMNAME_MAXLEN)
        return (SYMERR_TOOLONG);

    if (!num_symbols)
        init_hash ();
    if ((sym = add_symbol (name)) == NULL)
        return (SYMERR_SYMTABLE);
    vecsym_type(sym)   = VECSYM_VALUE;
    vecsym_nargs(sym)  = 0;
    vecsym_veclen(sym) = 0;
    vecsym_value(sym)  = val;
    vecsym_user(sym)   = user;
    return (0);
}

/*---------- sym_addvec --------------------------------------------
 * adds a vector to user symbol table
 *------------------------------------------------------------------*/

int sym_addvec (char *name, size_t len, VECFLOAT *vec, void *user)
{
    VECSYM *sym;

    if (name == NULL || !*name)
        return (SYMERR_NONAME);
    if (strlen (name) > SYMNAME_MAXLEN)
        return (SYMERR_TOOLONG);
    if (len < 1 || vec == NULL)
        return (SYMERR_NOVEC);

    if (!num_symbols)
        init_hash ();
    if ((sym = add_symbol (name)) == NULL)
        return (SYMERR_SYMTABLE);
    vecsym_vector(sym) = (VECFLOAT *) malloc (len * sizeof(VECFLOAT));
    if (vecsym_vector(sym) == NULL) {
        vecsym_type(sym) = VECSYM_VALUE;
        sym_delsym (name);
        return (SYMERR_MALLOC);
    }
    vecsym_type(sym)   = VECSYM_VECTOR;
    vecsym_nargs(sym)  = 0;
    vecsym_veclen(sym) = len;
    vecsym_user(sym)   = user;
    memcpy (vecsym_vector(sym), vec, len * sizeof(VECFLOAT));
    return (0);
}

/*---------- sym_addequ --------------------------------------------
 * adds an equation to user symbol table
 *------------------------------------------------------------------*/

int sym_addequ (char *name, int nargs, char *equ, void *user)
{
    VECSYM *sym;

    if (name == NULL || !*name)
        return (SYMERR_NONAME);
    if (strlen (name) > SYMNAME_MAXLEN)
        return (SYMERR_TOOLONG);
    if (equ == NULL)
        return (SYMERR_NOEXPR);
    while (*equ && isspace (*equ))
        equ++;
    if (!*equ)
        return (SYMERR_NOEXPR);

    if (!num_symbols)
        init_hash ();
    if ((sym = add_symbol (name)) == NULL)
        return (SYMERR_SYMTABLE);
    vecsym_equstr(sym) = (char *) malloc (strlen (equ) + 1);
    if (vecsym_equstr(sym) == NULL) {
        vecsym_type(sym) = VECSYM_VALUE;
        sym_delsym (name);
        return (SYMERR_MALLOC);
    }
    vecsym_type(sym)  = VECSYM_EQUSTR;
    vecsym_nargs(sym) = nargs;
    vecsym_user(sym)  = user;
    strcpy (vecsym_equstr(sym), equ);
    return (0);
}

/*---------- sym_addfunc -------------------------------------------
 * adds a function to user symbol table
 *------------------------------------------------------------------*/

int sym_addfunc (char *name, int nargs, VECFUNC func, void *user)
{
    VECSYM *sym;

    if (name == NULL || !*name)
        return (SYMERR_NONAME);
    if (strlen (name) > SYMNAME_MAXLEN)
        return (SYMERR_TOOLONG);
    if (func == NULL)
        return (SYMERR_NOFUNC);
    if (nargs > FUNC_MAXARGS)
        return (SYMERR_MAXARGS);

    if (!num_symbols)
        init_hash ();
    if ((sym = add_symbol (name)) == NULL)
        return (SYMERR_SYMTABLE);
    vecsym_type(sym)  = VECSYM_FUNC;
    vecsym_nargs(sym) = nargs;
    vecsym_func(sym)  = func;
    vecsym_user(sym)  = user;
    return (0);
}

/*---------- sym_addmacro ------------------------------------------
 * adds a macro to user symbol table
 *------------------------------------------------------------------*/

int sym_addmacro (char *name, int nargs, char *macro, void *user)
{
    VECSYM *sym;

    if (name == NULL || !*name)
        return (SYMERR_NONAME);
    if (strlen (name) > SYMNAME_MAXLEN)
        return (SYMERR_TOOLONG);
    if (macro == NULL)
        return (SYMERR_NOEXPR);
    while (*macro && isspace (*macro))
        macro++;
    if (!*macro)
        return (SYMERR_NOEXPR);

    if (!num_symbols)
        init_hash ();
    if ((sym = add_symbol (name)) == NULL)
        return (SYMERR_SYMTABLE);
    vecsym_macro(sym) = (char *) malloc (strlen (macro) + 1);
    if (vecsym_macro(sym) == NULL) {
        vecsym_type(sym) = VECSYM_VALUE;
        sym_delsym (name);
        return (SYMERR_MALLOC);
    }
    vecsym_type(sym)  = VECSYM_MACRO;
    vecsym_nargs(sym) = nargs;
    vecsym_user(sym)  = user;
    strcpy (vecsym_macro(sym), macro);
    return (0);
}

/*---------- sym_adddata -------------------------------------------
 * adds user data to user symbol table
 *------------------------------------------------------------------*/

int sym_adddata (char *name, void *data, void *user)
{
    VECSYM *sym;

    if (name == NULL || !*name)
        return (SYMERR_NONAME);
    if (strlen (name) > SYMNAME_MAXLEN)
        return (SYMERR_TOOLONG);

    if (!num_symbols)
        init_hash ();
    if ((sym = add_symbol (name)) == NULL)
        return (SYMERR_SYMTABLE);
    vecsym_type(sym)  = VECSYM_DATA;
    vecsym_nargs(sym) = 0;
    vecsym_data(sym)  = data;
    vecsym_user(sym)  = user;
    return (0);
}

/*---------- sym_delsym ----------------------------------------------
 * delete a symbol from the symbol table
 *--------------------------------------------------------------------*/

void sym_delsym (char *name)
{
    VECSYM *sym, *prev;
    int hash;

    if (!num_symbols || name == NULL || !*name)
        return;
    hash = hash_name (name);
    if ((prev = hash_table[hash]) == NULL)
        return;
    if (strcmp (name, vecsym_name(prev)) == 0) {
        hash_table[hash] = prev->next;
        if (sym_delfunc != NULL) (*sym_delfunc) (prev);
        if (vecsym_type(prev) == VECSYM_EQUSTR)
            free (vecsym_equstr(prev));
        else if (vecsym_type(prev) == VECSYM_MACRO)
            free (vecsym_macro(prev));
        else if (vecsym_type(prev) == VECSYM_VECTOR)
            free (vecsym_vector(prev));
        free (prev);
        return;
    }
    sym = prev->next;
    while (sym != NULL) {
        if (strcmp (name, vecsym_name(sym)) == 0) {
            prev->next = sym->next;
            if (sym_delfunc != NULL) (*sym_delfunc) (sym);
            if (vecsym_type(sym) == VECSYM_EQUSTR)
                free (vecsym_equstr(sym));
            else if (vecsym_type(sym) == VECSYM_MACRO)
                free (vecsym_macro(sym));
            else if (vecsym_type(sym) == VECSYM_VECTOR)
                free (vecsym_vector(sym));
            free (sym);
            return;
        }
        prev = sym;
        sym = prev->next;
    }
}

/*---------- sym_count ---------------------------------------------
 * return number of symbols of given type
 *------------------------------------------------------------------*/

int sym_count (int type)
{
    int n, cnt = 0;
    VECSYM *sym;

    if (num_symbols) {
        for (n = 0; n < HASH_SIZE; n++) {
            sym = hash_table[n];
            while (sym != NULL) {
                if (!type || vecsym_type(sym) == type)
                    cnt++;
                sym = sym->next;
            }
        }
    }
    return (cnt);
}

/*---------- sym_names ---------------------------------------------
 * return list of symbol names
 *------------------------------------------------------------------*/

char **sym_names (int type)
{
    int n, cnt = sym_count (type);
    VECSYM *sym;
    char **names = NULL;

    if (cnt) {
        names = (char **) malloc ((cnt+1) * sizeof(char *));
        if (NULL == names)
            return (NULL);
        for (n = 0, cnt = 0; n < HASH_SIZE; n++) {
            sym = hash_table[n];
            while (sym != NULL) {
                if (!type || vecsym_type(sym) == type)
                    names[cnt++] = vecsym_name(sym);
                sym = sym->next;
            }
        }
        names[cnt] = NULL;
    }
    return (names);
}

/*---------- sym_list ----------------------------------------------
 * lists symbols
 *------------------------------------------------------------------*/

int sym_list (int type, int (*func)(struct symbol_ *, void *), void *userdata)
{
    int n, cnt = 0;
    VECSYM *sym;

    if (num_symbols && NULL != func) {
        for (n = 0; n < HASH_SIZE; n++) {
            sym = hash_table[n];
            while (sym != NULL) {
                if (!type || vecsym_type(sym) == type)
                    cnt += (*func) (sym, userdata);
                sym = sym->next;
            }
        }
    }
    return (cnt);
}

/*---------- sym_print ---------------------------------------------
 * print symbol list
 *------------------------------------------------------------------*/

void sym_print (FILE *fp)
{
    int i, fld;
    VECSYM *sym;

    if (!num_symbols) {
        fprintf (fp, "\nNo Symbols defined\n");
        return;
    }
    fld = max_namelen > 5 ? max_namelen : 5;
    fprintf (fp, "Defined Symbols:\n");
    fprintf (fp, "%-*s  Type        Nargs/Expression/Value/Length\n",
        fld, "Name");
    i = fld + 43;
    while (i--)
        putc ('-', fp);
    putc ('\n', fp);
    for (i = 0; i < HASH_SIZE; i++) {
        sym = hash_table[i];
        while (sym != NULL) {
            fprintf (fp, "%-*.*s  ", fld, fld, vecsym_name(sym));
            if (vecsym_type(sym) == VECSYM_VALUE)
                fprintf (fp, "value       %g\n", vecsym_value(sym));
            else if (vecsym_type(sym) == VECSYM_VECTOR)
                fprintf (fp, "vector      %ld\n", (long)vecsym_veclen(sym));
            else if (vecsym_type(sym) == VECSYM_EQUSTR)
                fprintf (fp, "equation    %s\n", vecsym_equstr(sym));
            else if (vecsym_type(sym) == VECSYM_MACRO)
                fprintf (fp, "macro       %s\n", vecsym_macro(sym));
            else if (vecsym_type(sym) == VECSYM_FUNC) {
                fprintf (fp, "function    ");
                if (vecsym_nargs(sym) < 0)
                    fprintf (fp, "variable\n");
                else
                    fprintf (fp, "%d\n", vecsym_nargs(sym));
            }
            else if (vecsym_type(sym) == VECSYM_DATA)
                fprintf (fp, "user data\n");
            else
                fprintf (fp, "unknown type\n");
            sym = sym->next;
        }
    }
}
