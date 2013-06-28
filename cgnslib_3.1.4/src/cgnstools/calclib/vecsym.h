/*
 * vecsym.h - symbol table interface for vec.c
 */

#ifndef _VECSYM_H_
#define _VECSYM_H_

/*===================================================================
 * none of this is used directly by vec.c
 * The symbol structure is used indirectly, however, through
 * the macros defined below
 *===================================================================*/

#include <stdio.h>
#ifndef _VEC_H_
#include "vec.h"    /* for definition of VECDATA */
#endif

/*----- error return codes -----*/

enum SymErrors {
    SYMERR_NOERR=0, /* no error */
    SYMERR_NONAME,  /* NULL or blank symbol name */
    SYMERR_NOEXPR,  /* NULL or blank equation string */
    SYMERR_SYMTABLE,/* failed to add new symbol to symbol table */
    SYMERR_MALLOC,  /* malloc failed - out of space */
    SYMERR_NOFUNC,  /* NULL function pointer */
    SYMERR_MAXARGS, /* user function has too many arguments */
    SYMERR_NOVEC,   /* 0 length or NULL vector */
    SYMERR_TOOLONG, /* symbol name too long */
    SYMERR_INVALID  /* invalid error code */
};

/*----- symbol table structure -----*/

#define SYMNAME_MAXLEN 32

struct symbol_ {
    char name[SYMNAME_MAXLEN+1];
    int type;
    int nargs;
    size_t len;
    union {
        VECFLOAT val;
        VECFLOAT *vec;
        char *equ;
        VECFUNC func;
        void *data;
    } op;
    void *user;
    struct symbol_ *next;
};

/*----- callback for deleting a symbol -----*/

extern void (*sym_delfunc) (struct symbol_ *);

/*----- function prototypes -----*/

void sym_free (void); /* free symbol table */

char *sym_errmsg (  /* get error message */
    int errnum      /* error number */
);

int sym_addval (    /* add value symbol */
    char *name,     /* symbol name */
    VECFLOAT val,   /* symbol value */
    void *user      /* user data */
);

int sym_addvec (    /* add vector symbol */
    char *name,     /* symbol name */
    size_t len,     /* vector length */
    VECFLOAT *vec,  /* vector */
    void *user      /* user data */
);

int sym_addequ (    /* add equation symbol */
    char *name,     /* symbol name */
    int nargs,      /* number of arguments */
    char *equ,      /* equation string */
    void *user      /* user data */
);

int sym_addfunc (   /* add function symbol */
    char *name,     /* symbol name */
    int nargs,      /* number of arguments */
    VECFUNC func,   /* function pointer */
    void *user      /* user data */
);

int sym_addmacro (  /* add macro symbol */
    char *name,     /* symbol name */
    int nargs,      /* number of arguments */
    char *equ,      /* macro string */
    void *user      /* user data */
);

int sym_adddata (   /* add user data symbol */
    char *name,     /* symbol name */
    void *data,     /* symbol data */
    void *user      /* user data */
);

void sym_delsym (   /* remove symbol */
    char *name      /* symbol name */
);

int sym_count (     /* count number of symbols */
    int type        /* symbol type (0 for all) */
);

char **sym_names (  /* return NULL terminated list of symbol names */
    int type        /* symbol type (0 for all) */
);

int sym_list (      /* list symbols */
    int type,       /* symbol type (0 for all) */
    int (*func) (   /* callback function for list */
        struct symbol_ *,   /* symbol */
        void *      /* user data */
    ),
    void *userdata  /* user data */
);

void sym_print (    /* print list of symbols */
    FILE *fp        /* where output goes */
);

/*===================================================================
 * the following need to be defined for use by vec.c
 *===================================================================*/

/*----- symbol definition -----*/

typedef struct symbol_ VECSYM;

/*----- maximum arguments to a function -----*/

#define FUNC_MAXARGS    5

/*----- symbol types -----*/

#define VECSYM_VALUE    1   /* constant value */
#define VECSYM_VECTOR   2   /* vector */
#define VECSYM_EQUSTR   3   /* equation string */
#define VECSYM_FUNC     4   /* user function */
#define VECSYM_MACRO    5   /* macro string */
#define VECSYM_DATA     6   /* user data */

/*----- macros for accessing symbols -----*/

#define vecsym_name(SYM)    ((SYM)->name)
#define vecsym_type(SYM)    ((SYM)->type)
#define vecsym_user(SYM)    ((SYM)->user)

#define vecsym_value(SYM)   ((SYM)->op.val)

#define vecsym_veclen(SYM)  ((SYM)->len)
#define vecsym_vector(SYM)  ((SYM)->op.vec)

#define vecsym_nargs(SYM)   ((SYM)->nargs)
#define vecsym_equstr(SYM)  ((SYM)->op.equ)
#define vecsym_macro(SYM)   ((SYM)->op.equ)
#define vecsym_func(SYM)    ((SYM)->op.func)

#define vecsym_data(SYM)    ((SYM)->op.data)

/*----- function prototypes -----*/

VECSYM *find_symbol (   /* find a symbol */
    char *name,         /* symbol name */
    int is_func         /* non-zero if symbol is a function */
);

#endif  /* _VECSYM_H_ */
