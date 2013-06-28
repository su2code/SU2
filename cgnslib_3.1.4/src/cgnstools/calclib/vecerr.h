#ifndef _VECERR_H_
#define _VECERR_H_

extern int vec_maxerr;      /* max error code */
extern int vec_errnum;      /* error number */

/*----- user error handler -----*/

extern void (*vec_errhandler) ( /* user error handler */
    int errnum,             /* error number */
    char *errmsg,           /* error message */
    int offset,             /* offset in parse string */
    char *exp               /* expression being parsed */
);

/*----- error codes -----*/

enum VecErrors {
    VECERR_NOERR=0, /* sucessfull evaluation */
    VECERR_SYMBOL,  /* unrecognized or invalid symbol */
    VECERR_BADEXPR, /* bad expression */
    VECERR_NOEXPR,  /* NULL or blank expression string argument */
    VECERR_FUNC,    /* function name not followed by '(' */
    VECERR_DIVIDE,  /* math error: divide by 0 */
    VECERR_PAREN,   /* '(' without ')' */
    VECERR_IFELSE,  /* if (?) without else (:) */
    VECERR_EQUALS,  /* '=' instead of '==' */
    VECERR_NOTEQS,  /* '!' instead of '!=' */
    VECERR_MATH,    /* general math error */
    VECERR_MISSING, /* missing operator */
    VECERR_POINTER, /* user function did not update pointer */
    VECERR_CALLS,   /* too many recursive calls */
    VECERR_BADNUM,  /* single '.' not preceeded/followed by digit */
    VECERR_RECURSE, /* recursion in expression detected */
    VECERR_ARGS,    /* function argument list error */
    VECERR_MAXARGS, /* user function has too many arguments */
    VECERR_MALLOC,  /* malloc failed */
    VECERR_OPERR,   /* internal error - op code */
    VECERR_PUSH,    /* out of space for variable stack */
    VECERR_POP,     /* tried to pop empty variable stack */
    VECERR_OP_STK,  /* couldn't malloc/realloc operation stack */
    VECERR_VAR_STK, /* couldn't malloc variable stack */
    VECERR_BRACKET, /* '[' without ']' */
    VECERR_NOVECT,  /* couldn't malloc vector data */
    VECERR_VECLEN,  /* vector lengths not the same */
    VECERR_INDEX,   /* index not associated with a variable */
    VECERR_BADINDX, /* index out of range of vector data */
    VECERR_VECTYPE, /* bad vector data type */
    VECERR_ZEROLEN, /* vector len < 1 */
    VECERR_EQUARG,  /* equation argument out of range */
    VECERR_ACOSH,   /* argument to acosh < 1 */
    VECERR_ATANH,   /* argument to atanh <= -1 or >= 1 */
    VECERR_BADVEC,  /* bad vector specification */
    VECERR_USER,    /* user call-back function error */
    VECERR_INVALID  /* invalid error message number */
};

#ifdef INCLUDE_ERRDATA

/*----- error messages -----*/

static char *err_msg[] = {
    "sucessfull evaluation",
    "unrecognized or invalid symbol",
    "bad expression",
    "NULL or blank expression string argument",
    "function name not followed by '('",
    "math error: divide by 0",
    "'(' without ')'",
    "if (?) without else (:)",
    "'=' instead of '=='",
    "'!' instead of '!='",
    "                                             ",
    "missing operator",
    "user function did not update pointer",
    "too many recursive calls",
    "single '.' not preceeded or followed by digit",
    "recursion in expression detected",
    "function argument list error",
    "user function has too many arguments",
    "malloc failed",
    "internal error - op code",
    "out of space for variable stack",
    "tried to pop empty variable stack",
    "couldn't malloc/realloc operation stack",
    "couldn't malloc variable stack",
    "'[' without ']'",
    "couldn't malloc vector data",
    "vector lengths not the same",
    "index not associated with a variable",
    "index out of range of vector data",
    "bad vector data type",
    "vector len < 1",
    "equation argument out of range",
    "argument to acosh() < 1",
    "argument to atanh() <= -1 or >= 1",
    "bad vector specification",
    NULL,
    "invalid error message number"
};

#endif  /* INCLUDE_ERRDATA */

#endif  /* _VECERR_H_ */
