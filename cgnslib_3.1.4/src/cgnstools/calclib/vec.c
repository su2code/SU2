/*
 * vec.c - array processor for an equation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <setjmp.h>
#include <signal.h>
#include <time.h>

#include "vec.h"
#include "vecsym.h"

#define INCLUDE_ERRDATA
#include "vecerr.h"

/*----- machine dependent defines -----*/

#ifndef SIGNAL
#define SIGNAL      void
#endif

/*
#ifndef sgi
#define USE_MATHERR
#endif
*/

#ifndef SIG_ERR
# ifdef BADSIG
#  define SIG_ERR   BADSIG
# else
#  define SIG_ERR   ((SIGNAL (*)())-1)
# endif
#endif

/*----- defined data values -----*/

#define STACK_INC   50      /* amount to increase size of op stack */
#define MAX_ARGS    9       /* maximum function arguments */

#define E           2.71828182845904523536
#ifndef PI
# define PI         3.14159265358979323846
#endif
#define TOLERANCE   0.0000000001

#define OPERATORS   "^*/%+-<>=!&|?:),]"

enum Tokens {

    /* operation tokens */

    OP_NEGATE = 1,
    OP_POWER,
    OP_MULTIPLY,
    OP_DIVIDE,
    OP_MODULUS,
    OP_PLUS,
    OP_MINUS,
    OP_LESS,
    OP_LESS_EQ,
    OP_GREATER,
    OP_GREATER_EQ,
    OP_EQUALS,
    OP_NOT_EQUAL,
    OP_AND,
    OP_OR,
    OP_NOT,
    OP_IF,
    OP_ELSE,
    OP_ENDIF,
    OP_BEGIN_FUNC,
    OP_NOP,

    OP_VECTOR,
    OP_NUMBER,
    OP_INTEGER,
    OP_INDEX,
    OP_COUNTER,

    /* function tokens */

    FUNC_ABS,
    FUNC_ACOS,
    FUNC_ACOSH,
    FUNC_ASIN,
    FUNC_ASINH,
    FUNC_ATAN,
    FUNC_ATAN2,
    FUNC_ATANH,
    FUNC_CEIL,
    FUNC_COS,
    FUNC_COSH,
    FUNC_DEG,
    FUNC_EXP,
    FUNC_FACT,
    FUNC_FLOOR,
    FUNC_LOG,
    FUNC_LOG10,
    FUNC_MOD,
    FUNC_POW,
    FUNC_POW10,
    FUNC_RAD,
    FUNC_RAND,
    FUNC_ROUND,
    FUNC_SIN,
    FUNC_SINH,
    FUNC_SQR,
    FUNC_SQRT,
    FUNC_TAN,
    FUNC_TANH,

    FUNC_AVG,
    FUNC_MAX,
    FUNC_MIN,
    FUNC_SUM,
    FUNC_SIZE,
    FUNC_LEN
};

/*----- operation stack -----*/

struct op_stack_ {
    int op;
    size_t len;
    union {
        int ival;
        double num;
        VECFLOAT *vec;
    } op_val;
};

static struct op_stack_ *op_start;  /* operation stack */
static int op_size = 0;             /* current size of operation stack */
static int op_pos;                  /* current position in stack */
static int op_index;                /* counter for processing */

/*----- variable stack -----*/

static double *var_stack;           /* processing stack */
static int var_size = 0;            /* size of processing stack */
static int var_pos;                 /* current position in stack */

/*----- constant values -----*/

static double constants[] = {
#define CONSTANT_E      0
    E,
#define CONSTANT_PI     1
    PI,
#define CONSTANT_TOL    2
    TOLERANCE
};
#define CONSTANT_DIM    3

/*----- intrinsic symbol table -----*/

#define VEC_ACROSS      0
#define VEC_ALONG       1
#define VEC_VARIABLE    2

static struct intrinsic_ {
    char *name;
    int type;
    int nargs;
    int op;
}
intrinsics[] = {
    {"abs",     VEC_ACROSS,     1,  FUNC_ABS    },
    {"acos",    VEC_ACROSS,     1,  FUNC_ACOS   },
    {"acosh",   VEC_ACROSS,     1,  FUNC_ACOSH  },
    {"asin",    VEC_ACROSS,     1,  FUNC_ASIN   },
    {"asinh",   VEC_ACROSS,     1,  FUNC_ASINH  },
    {"atan",    VEC_ACROSS,     1,  FUNC_ATAN   },
    {"atan2",   VEC_ACROSS,     2,  FUNC_ATAN2  },
    {"atanh",   VEC_ACROSS,     1,  FUNC_ATANH  },
    {"avg",     VEC_ALONG,      1,  FUNC_AVG    },
    {"ceil",    VEC_ACROSS,     1,  FUNC_CEIL   },
    {"cos",     VEC_ACROSS,     1,  FUNC_COS    },
    {"cosh",    VEC_ACROSS,     1,  FUNC_COSH   },
    {"deg",     VEC_ACROSS,     1,  FUNC_DEG    },
    {"dim",     CONSTANT_DIM,   0,  OP_NUMBER   },
    {"e",       CONSTANT_E,     0,  OP_NUMBER   },
    {"exp",     VEC_ACROSS,     1,  FUNC_EXP    },
    {"fact",    VEC_ACROSS,     1,  FUNC_FACT   },
    {"floor",   VEC_ACROSS,     1,  FUNC_FLOOR  },
    {"index",   VEC_VARIABLE,   0,  OP_COUNTER  },
    {"len",     VEC_ALONG,      1,  FUNC_LEN    },
    {"log",     VEC_ACROSS,     1,  FUNC_LOG    },
    {"log10",   VEC_ACROSS,     1,  FUNC_LOG10  },
    {"mag",     VEC_ALONG,      1,  FUNC_SIZE   },
    {"max",     VEC_ALONG,      1,  FUNC_MAX    },
    {"min",     VEC_ALONG,      1,  FUNC_MIN    },
    {"mod",     VEC_ACROSS,     2,  FUNC_MOD    },
    {"pi",      CONSTANT_PI,    0,  OP_NUMBER   },
    {"pow",     VEC_ACROSS,     2,  FUNC_POW    },
    {"pow10",   VEC_ACROSS,     1,  FUNC_POW10  },
    {"rad",     VEC_ACROSS,     1,  FUNC_RAD    },
    {"rand",    VEC_ACROSS,     0,  FUNC_RAND   },
    {"round",   VEC_ACROSS,     1,  FUNC_ROUND  },
    {"sin",     VEC_ACROSS,     1,  FUNC_SIN    },
    {"sinh",    VEC_ACROSS,     1,  FUNC_SINH   },
    {"sqr",     VEC_ACROSS,     1,  FUNC_SQR    },
    {"sqrt",    VEC_ACROSS,     1,  FUNC_SQRT   },
    {"sum",     VEC_ALONG,      1,  FUNC_SUM    },
    {"tan",     VEC_ACROSS,     1,  FUNC_TAN    },
    {"tanh",    VEC_ACROSS,     1,  FUNC_TANH   },
    {"tol",     CONSTANT_TOL,   0,  OP_NUMBER   }
};

#define NUM_INTRINSICS  sizeof(intrinsics)/sizeof(struct intrinsic_)

/*----- list of intrinsics -----*/

static char *in_list[] = {
    "Intrinsic operators in order of precedence:",
    "Operator Type   Operators",
    "-------------------------",
    "Expression      ()      const   func()",
    "Unary           -       !",
    "Exponential     ^       **",
    "Multiplicative  *       /       %",
    "Additive        +       -",
    "Relational      <       <=      >       >=      ==      !=",
    "Logical         &[&]    |[|]    ?:",
    "",
    "Defined constants:",
    "Name            Value",
    "---------------------",
    "pi              3.14159265358979323846",
    "e               2.71828182845904523536",
    "tol             0.0000000001",
    "index           current array index",
    "dim             minimum array len",
    "<number>        any constant",
    "",
    "Defined functions:",
    "Name            Operation",
    "-------------------------",
    "abs()           absolute value",
    "acos()          arc cosine (radians)",
    "acosh()         arc hyperbolic cosine",
    "asin()          arc sine (radians)",
    "asinh()         arc hyperbolic sine",
    "atan()          arc tangent (radians)",
    "atanh()         arc hyperbolic tangent",
    "avg()           average value in a vector",
    "ceil()          next largest integer",
    "cos()           cosine (radians)",
    "cosh()          hyperbolic cosine",
    "deg()           convert from radians to degrees",
    "exp()           power of e",
    "fact()          factorial of a number",
    "floor()         next smallest integer",
    "len()           length of a vector",
    "log()           natural log",
    "log10()         base 10 log",
    "size()          magnitude of a vector",
    "max()           maximum value in a vector",
    "min()           minimum value in a vector",
    "pow10()         power of 10",
    "rad()           convert from degrees to radians",
    "rand()          random number between 0 and 1",
    "round()         round a number to integer",
    "sin()           sine (radians)",
    "sinh()          hyperbolic sine",
    "sqr()           square a number",
    "sqrt()          take square root of a number",
    "sum()           sum of values in a vector",
    "tan()           tangent (radians)",
    "tanh()          hyperbolic tangent",
    NULL
};

/*----- error handling -----*/

int vec_maxerr = VECERR_INVALID;
int vec_errnum = 0;

void (*vec_errhandler)(int,char *,int,char *) = NULL;
typedef SIGNAL (*SIGHANDLER)(int);

/*----- temporary vectors -----*/

typedef struct veclist_ {
    VECDATA *v;
    struct veclist_ *prev;
} VECLIST;

static VECLIST *veclist = NULL;

/*----- defined macros -----*/

#define advance(amt)    exp_ptr+=(amt)
#define skip_space(p)   while(*(p) && isspace(*(p)))(p)++
#define skip_int(p)     while(*(p) && isdigit(*(p)))(p)++

/*----- parsing data -----*/

static int num_calls = 0;       /* recursion counter */
static int max_calls = VEC_MAXCALLS;/* maximum recursive calls */
static size_t max_len = 0;      /* maximum vector length */
static size_t counter_len = 0;  /* min length of counter */
static VECCALL usr_func;        /* pointer to user supplied function */
static char *exp_start = NULL;  /* start of string being parsed */
static char *exp_ptr;           /* location in string being parsed */
static SIGHANDLER old_fpe_sig;  /* old SIGFPE handler */
static jmp_buf env_ptr;         /* enviroment for return from longjmp() */

static int nest_level;          /* nesting count for if-then-else */
static int vec_op;              /* doing operations on array lengths */

static int num_args = 0;          /* number arguments for equation */
static VECDATA *arglist[MAX_ARGS];/* equation argument list */

static int checking = 0;        /* set if checking equation only */

/*----- function prototypes -----*/

static VECDATA *parse();
static void *handle_error();
static void clean_up();

static int  do_intrinsic();
static int  do_symbol();
static int  get_args();
static void chk_index();
static void add_op();

static void logical();
static void relational();
static void additive();
static void multiplicative();
static void exponential();
static void unary();
static void expression();

static void push_var();
static void pop_var();
static double process_op();

static SIGNAL fpe_err (int signum);

/*==================================================================
 * user callable routines
 *==================================================================*/

/*---------- vec_create --------------------------------------------
 * creates a vector data structure
 *------------------------------------------------------------------*/

VECDATA *vec_create (int type, size_t len, int temp)
{
    VECDATA *vdata;

    if (type == VEC_VECTOR) {
        if (len < 1 && temp >= 0)
            return (handle_error (VECERR_ZEROLEN));
    }
    else if (type == VEC_VALUE)
        len = 0;
    else
        return (handle_error (VECERR_VECTYPE));

    vdata = (VECDATA *) calloc (1, sizeof(VECDATA));
    if (vdata == NULL)
        return (handle_error (VECERR_MALLOC));

    if (len && temp >= 0) {
        vdata->f.vec = (VECFLOAT *) calloc (len, sizeof(VECFLOAT));
        if (vdata->f.vec == NULL) {
            free (vdata);
            return (handle_error (VECERR_MALLOC));
        }
    }
    else
        vdata->f.vec = NULL;

    vdata->type = type;
    vdata->len  = len;

    /* add to temporary list */

    if (temp) {
        VECLIST *vl = (VECLIST *) malloc (sizeof(VECLIST));
        if (vl == NULL) {
            if (len)
                free (vdata->f.vec);
            free (vdata);
            return (handle_error (VECERR_MALLOC));
        }
        vl->v = vdata;
        vl->prev = veclist;
        veclist = vl;
    }

    return (vdata);
}

/*---------- vec_destroy -------------------------------------------
 * destroys a vector data structure
 *------------------------------------------------------------------*/

void vec_destroy (VECDATA *vdata)
{
    if (vdata != NULL) {
        vec_remove (vdata);
        if (vdata->type == VEC_VECTOR && vdata->f.vec != NULL)
            free (vdata->f.vec);
        free (vdata);
    }
}

/*---------- vec_add -----------------------------------------------
 * add vector to temporary vector list
 *------------------------------------------------------------------*/

int vec_add (VECDATA *vd)
{
    VECLIST *vl = veclist;

    if (vd == NULL)
        return (0);

    /* check if already in list */

    while (vl != NULL) {
        if (vl->v == vd)
            return (1);
        vl = vl->prev;
    }

    /* it's not, so add it */

    if ((vl = (VECLIST *) malloc (sizeof(VECLIST))) == NULL)
        return (0);
    vl->v = vd;
    vl->prev = veclist;
    veclist = vl;
    return (1);
}

/*---------- vec_remove --------------------------------------------
 * remove temporary vector from list
 *------------------------------------------------------------------*/

void vec_remove (VECDATA *vd)
{
    VECLIST *prev, *vlist = veclist;

    if (vd == NULL || vlist == NULL)
        return;
    if (vlist->v == vd) {
        veclist = vlist->prev;
        free (vlist);
        return;
    }
    prev = vlist->prev;
    while (prev != NULL) {
        if (prev->v == vd) {
            vlist->prev = prev->prev;
            free (prev);
            return;
        }
        vlist = prev;
        prev = vlist->prev;
    }
}

/*---------- vec_free ----------------------------------------------
 * releases internal resources
 *------------------------------------------------------------------*/

void vec_free (void)
{
    /* free operation stack */

    if (op_size) {
        op_size = 0;
        free (op_start);
    }

    /* free variable stack */

    if (var_size) {
        var_size = 0;
        free (var_stack);
    }
}

/*---------- vec_maxcalls ------------------------------------------
 * sets maximum recursive calls
 *------------------------------------------------------------------*/

void vec_maxcalls (int val)
{
    max_calls = val > 0 ? val : VEC_MAXCALLS;
}

/*---------- vec_errmsg --------------------------------------------
 * returns error message
 *------------------------------------------------------------------*/

char *vec_errmsg (void)
{
    if (vec_errnum < 0 || vec_errnum >= VECERR_INVALID)
        vec_errnum = VECERR_INVALID;
    return (err_msg[vec_errnum]);
}

/*---------- vec_list ----------------------------------------------
 * return list of intrinsics
 *------------------------------------------------------------------*/

char **vec_list (void)
{
    return (in_list);
}

/*---------- vec_parse ---------------------------------------------
 * parse and process equation string
 *------------------------------------------------------------------*/

VECDATA *vec_parse (char *str, int min_len, VECCALL func)
{
    VECDATA *vd;
    size_t old_len = counter_len;
    VECCALL old_func = usr_func;
    char *old_ptr = exp_ptr;
    char *old_start = exp_start;

    /* is there a string ? */

    exp_start = NULL;
    if (str != NULL)
        skip_space (str);
    if (str == NULL || !*str)
        return (handle_error (VECERR_NOEXPR));

    /* create operation stack if not done */

    if (op_size == 0) {
        op_start = (struct op_stack_ *)
            malloc (STACK_INC * sizeof(struct op_stack_));
        if (op_start == NULL)
            return (handle_error (VECERR_OP_STK));
        op_size = STACK_INC;
    }

    /* create processing stack if not already done */

    if (var_size == 0) {
        var_stack = (double *) malloc (STACK_INC * sizeof(double));
        if (var_stack == NULL)
            return (handle_error (VECERR_VAR_STK));
        var_size = STACK_INC;
    }

    /* install the SIGFPE handler */

    if (num_calls == 0) {
        checking = op_pos = var_pos = op_index = 0;
        max_len = 0;
        old_fpe_sig = signal (SIGFPE, fpe_err);
        vec_randinit ((int) time (NULL));

        /* jump here on parsing error */

        vec_errnum = setjmp (env_ptr);
        if (vec_errnum != 0) {
            clean_up ();
            return (handle_error (vec_errnum));
        }
    }

    usr_func = func;
    exp_ptr = str;
    exp_start = str;
    counter_len = min_len > 0 ? min_len : 1;

    vd = parse ();
    if (*exp_ptr)
        longjmp (env_ptr, VECERR_MISSING);
    vec_remove (vd);

    if (num_calls == 0)
        clean_up ();
    else {
        counter_len = old_len;
        usr_func  = old_func;
        exp_ptr   = old_ptr;
        exp_start = old_start;
    }
    return (vd);
}

/*---------- vec_check ---------------------------------------------
 * parse equation string and check for errors
 *------------------------------------------------------------------*/

size_t vec_check (char *str, int min_len, VECCALL func)
{
    VECDATA *vd;
    size_t len;
    size_t old_len = counter_len;
    VECCALL old_func = usr_func;
    char *old_ptr = exp_ptr;
    char *old_start = exp_start;

    /* is there a string ? */

    exp_start = NULL;
    if (str != NULL)
        skip_space (str);
    if (str == NULL || !*str) {
        handle_error (VECERR_NOEXPR);
        return (0);
    }

    /* create operation stack if not done */

    if (op_size == 0) {
        op_start = (struct op_stack_ *)
            malloc (STACK_INC * sizeof(struct op_stack_));
        if (op_start == NULL) {
            handle_error (VECERR_OP_STK);
            return (0);
        }
        op_size = STACK_INC;
    }

    /* create processing stack if not already done */

    if (var_size == 0) {
        var_stack = (double *) malloc (STACK_INC * sizeof(double));
        if (var_stack == NULL) {
            handle_error (VECERR_VAR_STK);
            return (0);
        }
        var_size = STACK_INC;
    }

    /* install the SIGFPE handler */

    if (num_calls == 0) {
        checking = 1;
        max_len = op_pos = var_pos = op_index = 0;
        old_fpe_sig = signal (SIGFPE, fpe_err);
        vec_randinit ((int) time (NULL));

        /* jump here on parsing error */

        vec_errnum = setjmp (env_ptr);
        if (vec_errnum != 0) {
            clean_up ();
            handle_error (vec_errnum);
            return (0);
        }
    }

    usr_func = func;
    exp_ptr = str;
    exp_start = str;
    counter_len = min_len > 0 ? min_len : 1;

    vd = parse ();
    if (*exp_ptr)
        longjmp (env_ptr, VECERR_MISSING);
    len = vd->len ? vd->len : 1;

    if (num_calls == 0)
        clean_up ();
    else {
        counter_len = old_len;
        usr_func  = old_func;
        exp_ptr   = old_ptr;
        exp_start = old_start;
    }
    return (len);
}

/*===================================================================
 * support routines callable by user
 *===================================================================*/

/*---------- vec_number --------------------------------------------
 * gets a double and advances string pointer
 * This is similar to 'strtod, but I know what it's doing
 * Expects a number of the form :
 * [whitespace][{+|-}][digits][.digits][{d|D|e|E}[sign]digits]
 *------------------------------------------------------------------*/

int vec_number (double *dval, char **sp)
{
    register char *p, c;
    char *start;
    int point = 0;

    p = start = *sp;
    *dval = 0.0;

    /* skip whitespace */

    while (*p && isspace (*p))
        p++;

    /* get to first required digit */

    if (*p == '-' || *p == '+')
        p++;
    if (*p == '.') {
        point++;
        p++;
    }

    /* next character needs to be digit */

    if (!*p || !isdigit (*p))
        return (0);

    /* get past the number */

    while (*p && (isdigit (*p) || *p == '.')) {
        if (*p == '.') {
            if (point) {
                *p = 0;
                *dval = (double) atof (*sp);
                *p = '.';
                *sp = p;
                return (1);
            }
            point++;
        }
        p++;
    }
    *sp = p;

    /* check for an exponent */

    if (*p == 'd' || *p == 'D' || *p == 'e' || *p == 'E') {
        if (*++p == '+' || *p == '-')
            p++;
        if (isdigit (*p)) {
            while (*++p && isdigit (*p))
                ;
            *sp = p;
        }
    }
    c = **sp;
    **sp = 0;
    *dval = atof (start);
    **sp = c;
    return (1);
}

/*-------------------------------------------------------------------
 * system independent random number generator from Kluth
 *-------------------------------------------------------------------*/

#define MBIG    1000000000
#define MSEED   161803398

static int Rand_ext = 0;
static int Rand_extp = 31;
static long Rand_ma[56] = {
            0,473269125,708359770,742015427, 65602533,275831218, 24831226,
    897336983,199633574,555982024,514994949,279315436,156002375,334914356,
    270182445,993586635,796323046,222779050,530019863,240216188,246247465,
    251350083, 27923555, 17362715,286349234,561741882, 61883183, 25293241,
    182316584,384320687, 97284929,343171996,939345275,385350967,340911449,
    606343026,885561620,105803147,288011295,407490891,632823362,921005485,
    393546951,638589878,430524660,  1651896,884594510,251198018,883223679,
    254238950,266438063,664142955,409571047,306976444,378529592,649134132
};

/*---------- vec_randinit -------------------------------------------
 * initialize the random number generator
 *-------------------------------------------------------------------*/

void vec_randinit (int seed)
{
    int i, k;
    long mj, mk;

    mj = MSEED - (seed < 0 ? -seed : seed);
    mj %= MBIG;
    Rand_ma[55] = mj;
    mk = 1;
    for (i = 1; i <= 54; i++) {
        k = (21 * i) % 55;
        Rand_ma[k] = mk;
        mk = mj - mk;
        if (mk < 0)
            mk += MBIG;
        mj = Rand_ma[k];
    }
    for (k = 1; k <= 4; k++) {
        for (i = 1; i <= 55; i++) {
            Rand_ma[i] -= Rand_ma[1+(i+30) % 55];
            if (Rand_ma[i] < 0)
                Rand_ma[i] += MBIG;
        }
    }
    Rand_ext = 0;
    Rand_extp = 31;
}

/*---------- vec_rand -----------------------------------------------
 * return random number between 0 and 1
 *-------------------------------------------------------------------*/

double vec_rand (void)
{
    long mj;

    if (++Rand_ext == 56)
        Rand_ext = 1;
    if (++Rand_extp == 56)
        Rand_extp = 1;
    mj = Rand_ma[Rand_ext] - Rand_ma[Rand_extp];
    if (mj < 0)
        mj += MBIG;
    Rand_ma[Rand_ext] = mj;
    return ((double) mj / (double) MBIG);
}

/*==================================================================
 * parsing routines
 *==================================================================*/

/*---------- parse --------------------------------------------------
 * does the parsing and evaluation of the string
 *-------------------------------------------------------------------*/

static VECDATA *parse ()
{
    size_t n;
    size_t cur_len = max_len;
    int cur_op    = op_pos;
    int cur_index = op_index;
    int cur_var   = var_pos;
    VECDATA *vd;

    if (++num_calls == max_calls)
        longjmp (env_ptr, VECERR_CALLS);

    max_len = 0;
    logical ();
    skip_space (exp_ptr);

    if (checking) {
        if (max_len == 1)
            vd = vec_create (VEC_VALUE, 0, -1);
        else
            vd = vec_create (VEC_VECTOR, max_len, -1);
    }

    else {
        vd = vec_create (VEC_VECTOR, max_len, 1);

        /* process the variables */

        for (n = 0; n < max_len; n++) {
            var_pos  = cur_var;
            op_index = cur_op;
            nest_level = vec_op = 0;
            vd->f.vec[n] = (VECFLOAT)process_op (n);
        }

        if (max_len == 1) {
            VECFLOAT f = vd->f.vec[0];
            free (vd->f.vec);
            vd->type  = VEC_VALUE;
            vd->len   = 0;
            vd->f.val = f;
        }
    }

    max_len  = cur_len;
    op_pos   = cur_op;
    op_index = cur_index;
    var_pos  = cur_var;
    num_calls--;
    return (vd);
}

/*---------- handle_error -------------------------------------------
 * handle error message
 *-------------------------------------------------------------------*/

static void *handle_error (err)
int err;
{
    int pos = -1;
    char *str = exp_start;

    if (num_calls)
        longjmp (env_ptr, vec_errnum);
    if (err < 0 || err >= VECERR_INVALID)
        vec_errnum = VECERR_INVALID;
    else
        vec_errnum = err;
    exp_start = NULL;
    if (str != NULL) {
        pos = (int) (exp_ptr - str);
        if (pos < 0 || pos > (int)strlen (str)) {
            pos = -1;
            str = NULL;
        }
    }
    if (vec_errhandler != NULL)
        (*vec_errhandler) (vec_errnum, err_msg[vec_errnum], pos, str);
    else {
        if (str != NULL) {
            fprintf (stderr, "%s\n", str);
            while (pos-- > 0)
                putc ('-', stderr);
            putc ('^', stderr);
        }
        fprintf (stderr, "%s\n", err_msg[vec_errnum]);
    }
    return (NULL);
}

/*---------- reset_recurs ------------------------------------------
 * reset recursion flags for equation symbols - called from sym_list
 *------------------------------------------------------------------*/

static int reset_recurs (VECSYM *sym, void *userdata)
{
    if (vecsym_nargs(sym) < 0)
        vecsym_nargs(sym) = 0;
    return (0);
}

/*---------- clean_up ----------------------------------------------
 * clean up on error or parsing completion
 *------------------------------------------------------------------*/

static void clean_up ()
{
    VECLIST *prev;

    if (old_fpe_sig != SIG_ERR)
        signal (SIGFPE, old_fpe_sig);
    (void) sym_list (VECSYM_EQUSTR, reset_recurs, NULL);
    while (veclist != NULL) {
        prev = veclist->prev;
        if (veclist->v->type == VEC_VECTOR &&
            veclist->v->f.vec != NULL)
            free (veclist->v->f.vec);
        free (veclist->v);
        free (veclist);
        veclist = prev;
    }
    num_calls = num_args = 0;
}

/*---------- do_intrinsic ------------------------------------------
 * process an intrinsic variable
 *------------------------------------------------------------------*/

static int do_intrinsic ()
{
    char c, *p = exp_ptr;
    int lo, hi, n, m;
    struct intrinsic_ *ins;

    if (!isalpha (*p))
        return (0);
    while (*++p && (isalnum (*p) || *p == '_'))
        ;
    c = *p;
    *p = 0;

    lo  = 0;
    hi  = NUM_INTRINSICS - 1;
    n   = hi / 2;
    ins = &intrinsics[n];
    while (1) {
        if ((m = strcmp (exp_ptr, ins->name)) == 0)
            break;
        if (lo >= hi)
            break;
        if (m < 0)
            hi = n - 1;
        else
            lo = n + 1;
        n = (lo + hi) / 2;
        ins = &intrinsics[n];
    }
    *p = c;

    if (m)
        return (0);
    exp_ptr = p;

    if (ins->op == OP_NUMBER) {
        if (ins->type == CONSTANT_DIM) {
            double dval = counter_len;
            add_op (OP_NUMBER, &dval);
        }
        else
            add_op (OP_NUMBER, &constants[ins->type]);
    }
    else if (ins->type == VEC_VARIABLE)
        add_op (ins->op, NULL);
    else if (ins->type == VEC_ALONG) {
        size_t msave = max_len;
        add_op (OP_BEGIN_FUNC, NULL);
        (void) get_args (ins->nargs);
        add_op (ins->op, NULL);
        max_len = msave ? msave : 1;
    }
    else {
        (void) get_args (ins->nargs);
        add_op (ins->op, NULL);
    }
    return (1);
}

/*---------- get_args ----------------------------------------------
 * extracts arguments for intrinsic function
 *------------------------------------------------------------------*/

static int get_args (nargs)
int nargs;
{
    int n = 0;

    skip_space (exp_ptr);
    if (*exp_ptr != '(')
        longjmp (env_ptr, VECERR_FUNC);

    advance (1);
    skip_space (exp_ptr);
    if (*exp_ptr == ')') {
        if (nargs > 0)
            longjmp (env_ptr, VECERR_ARGS);
        advance (1);
        return (0);
    }

    while (1) {
        skip_space (exp_ptr);
        logical ();
        n++;
        skip_space (exp_ptr);
        if (*exp_ptr == ')')
            break;
        if (*exp_ptr == ',') {
            if (n == nargs)
                longjmp (env_ptr, VECERR_ARGS);
            if (n == FUNC_MAXARGS)
                longjmp (env_ptr, VECERR_MAXARGS);
            advance (1);
        }
        else
            longjmp (env_ptr, VECERR_ARGS);
    }
    if (nargs > 0 && n != nargs)
        longjmp (env_ptr, VECERR_ARGS);
    advance (1);
    return (n);
}

/*---------- do_symbol ---------------------------------------------
 * processes a user symbol
 *------------------------------------------------------------------*/

static int do_symbol ()
{
    int c, is_func;
    char symname[SYMNAME_MAXLEN+2], *p = exp_ptr;
    VECSYM *sym;

    /* check for replaceable argument */

    if (*p == '$') {
        if (*++p && isdigit (*p)) {
            c = *p - '1';
            if (c < 0 || c >= num_args)
                longjmp (env_ptr, VECERR_EQUARG);
            exp_ptr = p;
            advance (1);
            if (arglist[c]->type == VEC_VALUE) {
                double dval = arglist[c]->f.val;
                add_op (OP_NUMBER, &dval);
            }
            else if (arglist[c]->type == VEC_VECTOR)
                add_op (OP_VECTOR, arglist[c]);
            else
                longjmp (env_ptr, VECERR_VECTYPE);
            return (1);
        }
        return (0);
    }

    /* get symbol name */

    if (*p == '\\') p++;
    if (!*p || (!isalpha (*p) && *p != '_' && *p != '"'))
        return (0);

    symname[0] = '.';
    if (*p == '"') {
        for (++p, c = 1; c <= SYMNAME_MAXLEN; c++) {
            if (!*p || *p == '"') break;
            symname[c] = *p++;
        }
        if (*p++ != '"') return 0;
    }
    else {
        for (c = 1; c <= SYMNAME_MAXLEN; c++) {
            symname[c] = *p++;
            if (!*p || (!isalnum (*p) && *p != '_')) break;
        }
    }
    if (c > SYMNAME_MAXLEN) return (0);
    symname[++c] = 0;

    /* check if function */

    while (*p && isspace (*p))
        p++;
    is_func = *p == '(';

    /* find symbol */

    if (*exp_ptr == '\\' || (sym = find_symbol (symname, is_func)) == NULL)
        sym = find_symbol (&symname[1], is_func);

    /*--- symbol not found ---*/

    if (sym == NULL)
        return (0);

    /*--- value ---*/

    if (vecsym_type(sym) == VECSYM_VALUE) {
        double val = vecsym_value(sym);
        add_op (OP_NUMBER, &val);
        exp_ptr = p;
        return (1);
    }

    /*--- vector ---*/

    if (vecsym_type(sym) == VECSYM_VECTOR) {
        VECDATA vdat;
        exp_ptr = p;
        vdat.len   = vecsym_veclen(sym);
        vdat.f.vec = vecsym_vector(sym);
        chk_index (&vdat);
        return (1);
    }

    /*--- equation ---*/

    if (vecsym_type(sym) == VECSYM_EQUSTR) {
        char *s;
        int nargs = 0;
        VECDATA *vd[9], *vtmp;
        if (vecsym_nargs(sym) < 0)
            longjmp (env_ptr, VECERR_RECURSE);
        if (vecsym_nargs(sym) == 0) {
            vecsym_nargs(sym) = -1;
            if (is_func) {
                p++;
                skip_space (p);
                if (*p != ')')
                    longjmp (env_ptr, VECERR_PAREN);
                p++;
            }
        }
        else {
            exp_ptr = p;
            if (*exp_ptr++ != '(')
                longjmp (env_ptr, VECERR_FUNC);
            skip_space (exp_ptr);
            while (*exp_ptr && *exp_ptr != ')') {
                if (nargs == vecsym_nargs(sym))
                    longjmp (env_ptr, VECERR_MAXARGS);
                vd[nargs++] = parse ();
                if (*exp_ptr == ',')
                    advance (1);
                skip_space (exp_ptr);
            }
            if (*exp_ptr != ')' || nargs != vecsym_nargs(sym))
                longjmp (env_ptr, VECERR_ARGS);
            advance (1);
            p = exp_ptr;
        }
        for (c = 0; c < 9; c++) {
            vtmp = arglist[c];
            arglist[c] = vd[c];
            vd[c] = vtmp;
        }
        c = nargs;
        nargs = num_args;
        num_args = c;
        s = exp_start;
        exp_start = exp_ptr = vecsym_equstr(sym);
        logical ();
        if (*exp_ptr)
            longjmp (env_ptr, VECERR_MISSING);
        if (vecsym_nargs(sym) < 0)
            vecsym_nargs(sym) = 0;
        exp_start = s;
        exp_ptr = p;
        num_args = nargs;
        for (c = 0; c < num_args; c++)
            arglist[c] = vd[c];
        return (1);
    }

    /*--- function ---*/

    if (vecsym_type(sym) == VECSYM_FUNC) {
        int nargs = 0;
        VECDATA *vd[FUNC_MAXARGS+1];
        exp_ptr = p;
        if (*exp_ptr++ != '(')
            longjmp (env_ptr, VECERR_FUNC);
        skip_space (exp_ptr);
        if (*exp_ptr != ')') {
            while (1) {
                vd[nargs] = parse ();
                if (++nargs > FUNC_MAXARGS)
                    longjmp (env_ptr, VECERR_MAXARGS);
                if (*exp_ptr == ',')
                    advance (1);
                else if (*exp_ptr != ')')
                    longjmp (env_ptr, VECERR_ARGS);
                else
                    break;
            }
        }
        advance (1);
        if (vecsym_nargs(sym) >= 0 && vecsym_nargs(sym) != nargs)
            longjmp (env_ptr, VECERR_ARGS);

        p = NULL;
        vd[nargs] = (*vecsym_func(sym)) (checking, nargs, vd, &p);
        if (vd[nargs] == NULL) {
            if (p == NULL || !*p)
                longjmp (env_ptr, VECERR_SYMBOL);
            else {
                err_msg[VECERR_USER] = p;
                longjmp (env_ptr, VECERR_USER);
            }
        }
        if (vd[nargs]->type == VEC_VALUE) {
            double dval = vd[nargs]->f.val;
            add_op (OP_NUMBER, &dval);
        }
        else if (vd[nargs]->type == VEC_VECTOR)
            add_op (OP_VECTOR, vd[nargs]);
        else
            longjmp (env_ptr, VECERR_VECTYPE);
        for (c = 0; c < nargs; c++)
            vec_destroy (vd[c]);
        return (1);
    }

    /*--- unknown ---*/

    return (0);
}

/*---------- chk_index ---------------------------------------------
 * check for a vector index
 *------------------------------------------------------------------*/

static void chk_index (vdata)
VECDATA *vdata;
{
    if (*exp_ptr == '[') {
        advance (1);
        logical ();
        skip_space (exp_ptr);
        if (*exp_ptr != ']')
            longjmp (env_ptr, VECERR_BRACKET);
        advance (1);
        add_op (OP_INDEX, NULL);
        if (max_len <= 1) {
            add_op (OP_VECTOR, vdata);
            max_len = 1;
            return;
        }
    }
    add_op (OP_VECTOR, vdata);
}

/*---------- do_userfunc -------------------------------------------
 * call user routine
 *------------------------------------------------------------------*/

static int do_userfunc ()
{
    char *p, *s;
    double dval;
    VECDATA *vd;

    if (usr_func == NULL)
        return (0);
    p = exp_ptr;
    s = NULL;
    if ((vd = (*usr_func) (checking, &p, &s)) == NULL) {
        if (s == NULL || !*s)
            return (0);
        err_msg[VECERR_USER] = s;
        longjmp (env_ptr, VECERR_USER);
    }
    if (p <= exp_ptr)
        longjmp (env_ptr, VECERR_POINTER);
    exp_ptr = p;
    if (vd->type == VEC_VALUE) {
        dval = vd->f.val;
        add_op (OP_NUMBER, (void *)&dval);
    }
    else if (vd->type == VEC_VECTOR) {
        skip_space (exp_ptr);
        chk_index (vd);
    }
    else
        longjmp (env_ptr, VECERR_VECTYPE);
    return (1);
}

/*---------- add_op ------------------------------------------------
 * adds operation to the stack
 *------------------------------------------------------------------*/

static void add_op (op, op_val)
int op;
void *op_val;
{
    if (op_pos >= op_size) {
        op_size += STACK_INC;
        op_start = (struct op_stack_ *)
            realloc (op_start, op_size * sizeof(struct op_stack_));
        if (op_start == NULL)
            longjmp (env_ptr, VECERR_OP_STK);
    }

    switch (op) {
        case OP_INTEGER:
            op_start[op_pos].len = 1;
            op_start[op_pos].op_val.ival = *((int *)op_val);
            break;
        case OP_NUMBER:
            op_start[op_pos].len = 1;
            op_start[op_pos].op_val.num = *((double *)op_val);
            break;
        case OP_VECTOR:
        {
            VECDATA *vdat = (VECDATA *)op_val;
            op_start[op_pos].len = vdat->len;
            op_start[op_pos].op_val.vec = vdat->f.vec;
/*            if (counter_len < vdat->len)*/
                counter_len = vdat->len;
            break;
        }
        case OP_INDEX:
        case FUNC_RAND:
            op_start[op_pos].len = 1;
            break;
        case OP_COUNTER:
            op_start[op_pos].len = counter_len;
            break;
        default:
            op_start[op_pos].len = 0;
            break;
    }

    op_start[op_pos].op = op;
    if (max_len < op_start[op_pos].len)
        max_len = op_start[op_pos].len;
    op_pos++;
}

/*---------- logical -----------------------------------------------
 * evaluates logical operations
 *------------------------------------------------------------------*/

static void logical ()
{
    relational ();

    while (1) {
        skip_space (exp_ptr);

        if (*exp_ptr == '?') {
            advance (1);
            add_op (OP_IF, NULL);
            logical ();
            skip_space (exp_ptr);
            if (*exp_ptr != ':')
                longjmp (env_ptr, VECERR_IFELSE);
            advance (1);
            add_op (OP_ELSE, NULL);
            logical ();
            add_op (OP_ENDIF, NULL);
        }
        else if (*exp_ptr == '&') {
            advance (1);
            if (*exp_ptr == '&')
                advance (1);
            relational ();
            add_op (OP_AND, NULL);
        }
        else if (*exp_ptr == '|') {
            advance (1);
            if (*exp_ptr == '|')
                advance (1);
            relational ();
            add_op (OP_OR, NULL);
        }
        else
            break;
    }
}

/*---------- relational --------------------------------------------
 * evaluates relational operators
 *------------------------------------------------------------------*/

static void relational ()
{
    additive ();

    while (1) {
        skip_space (exp_ptr);

        if (*exp_ptr == '<') {
            advance (1);
            if (*exp_ptr == '=') {
                advance (1);
                additive ();
                add_op (OP_LESS_EQ, NULL);
            }
            else {
                additive ();
                add_op (OP_LESS, NULL);
            }
        }
        else if (*exp_ptr == '>') {
            advance (1);
            if (*exp_ptr == '=') {
                advance (1);
                additive ();
                add_op (OP_GREATER_EQ, NULL);
            }
            else {
                additive ();
                add_op (OP_GREATER, NULL);
            }
        }
        else if (*exp_ptr == '=') {
            advance (1);
            if (*exp_ptr != '=')
                longjmp (env_ptr, VECERR_EQUALS);
            advance (1);
            additive ();
            add_op (OP_EQUALS, NULL);
        }
        else if (*exp_ptr == '!') {
            advance (1);
            if (*exp_ptr != '=')
                longjmp (env_ptr, VECERR_NOTEQS);
            advance (1);
            additive ();
            add_op (OP_NOT_EQUAL, NULL);
        }
        else
            break;
    }
}

/*---------- additive ----------------------------------------------
 * evaluates additive operators
 *------------------------------------------------------------------*/

static void additive ()
{
    multiplicative ();

    while (1) {
        skip_space (exp_ptr);

        if (*exp_ptr == '+') {
            advance (1);
            multiplicative ();
            add_op (OP_PLUS, NULL);
        }
        else if (*exp_ptr == '-') {
            advance (1);
            multiplicative ();
            add_op (OP_MINUS, NULL);
        }
        else
            break;
    }
}

/*---------- multiplicative ----------------------------------------
 * evaluates multiplicative operators
 *------------------------------------------------------------------*/

static void multiplicative ()
{
    exponential ();

    while (1) {
        skip_space (exp_ptr);

        if (*exp_ptr == '*') {
            advance (1);
            exponential ();
            add_op (OP_MULTIPLY, NULL);
        }
        else if (*exp_ptr == '/') {
            advance (1);
            exponential ();
            add_op (OP_DIVIDE, NULL);
        }
        else if (*exp_ptr == '%') {
            advance (1);
            exponential ();
            add_op (OP_MODULUS, NULL);
        }
        else if (strchr (OPERATORS, *exp_ptr) == NULL) {
            exponential ();
            add_op (OP_MULTIPLY, NULL);
        }
        else
            break;
    }
}

/*---------- exponential -------------------------------------------
 * evaluates exponential operators
 *------------------------------------------------------------------*/

static void exponential ()
{
    unary ();

    while (1) {
        skip_space (exp_ptr);

        if (*exp_ptr == '^') {
            advance (1);
            unary ();
            add_op (OP_POWER, NULL);
        }
        else if (*exp_ptr == '*' && *(exp_ptr+1) == '*') {
            advance (2);
            unary ();
            add_op (OP_POWER, NULL);
        }
        else
            break;
    }
}

/*---------- unary -------------------------------------------------
 * evaluates unary operations
 *------------------------------------------------------------------*/

static void unary ()
{
    skip_space (exp_ptr);

    if (*exp_ptr == '-') {
        advance (1);
        expression ();
        add_op (OP_NEGATE, NULL);
    }
    else if (*exp_ptr == '!') {
        advance (1);
        expression ();
        add_op (OP_NOT, NULL);
    }
    else
        expression ();
}

/*---------- expression --------------------------------------------
 * evaluates an expression
 *------------------------------------------------------------------*/

static void expression ()
{
    char *p;
    int i;
    size_t ii;
    double dval;
    VECDATA *vd;

    skip_space (exp_ptr);

    if (*exp_ptr == '(') {
        advance (1);
        logical ();
        skip_space (exp_ptr);
        if (*exp_ptr != ')')
            longjmp (env_ptr, VECERR_PAREN);
        advance (1);
    }

    else if (*exp_ptr == '[') {
        advance (1);
        p = exp_ptr;
        i = 0;
        while (*exp_ptr && *exp_ptr != ']' &&
               vec_number (&dval, &exp_ptr)) {
            i++;
            skip_space (exp_ptr);
            if (*exp_ptr == ',')
                advance(1);
        }
        if (!i || *exp_ptr != ']')
            longjmp (env_ptr, VECERR_BADVEC);
        advance (1);
        vd = vec_create (VEC_VECTOR, i, 1);
        for (ii = 0; ii < vd->len; ii++) {
            vec_number (&dval, &p);
            vd->f.vec[ii] = (VECFLOAT)dval;
            skip_space (p);
            if (*p == ',')
                p++;
        }
        add_op (OP_VECTOR, (void *)vd);
    }

    else if (isdigit (*exp_ptr) || *exp_ptr == '.') {
        if (!vec_number (&dval, &exp_ptr))
            longjmp (env_ptr, VECERR_BADNUM);
        add_op (OP_NUMBER, (void *)&dval);
    }

    else if (*exp_ptr == '\\') {
        advance (1);
        if (!do_userfunc() && !do_symbol() && !do_intrinsic())
            longjmp (env_ptr, VECERR_SYMBOL);
    }

    else if (!do_intrinsic() && !do_symbol() && !do_userfunc())
        longjmp (env_ptr, VECERR_SYMBOL);
}

/*==================================================================
 * processing routines
 *==================================================================*/

/*---------- push_var ----------------------------------------------
 * push a variable on the stack
 *------------------------------------------------------------------*/

static void push_var (v)
double v;
{
    if (var_pos >= var_size) {
        var_size += STACK_INC;
        var_stack = (double *)
             realloc (var_stack, var_size * sizeof(double));
        if (var_stack == NULL)
            longjmp (env_ptr, VECERR_PUSH);
    }
    var_stack[var_pos++] = v;
}

/*---------- pop_var ------------------------------------------------
 * pop a variable from the stack
 *------------------------------------------------------------------*/

static void pop_var (v)
double *v;
{
    if (--var_pos < 0)
        longjmp (env_ptr, VECERR_POP);
    *v = var_stack[var_pos];
}

/*---------- process_op --------------------------------------------
 * processes the operation stack for given variables
 *------------------------------------------------------------------*/

static double process_op (n)
int n;
{
    int m, done = 0, level, index = 0;
    long j;
    struct op_stack_ *op;
    double f1, f2;

    while (op_index < op_pos && !done) {
        op = &op_start[op_index++];
        switch (op->op) {
            case OP_NOP:
            case OP_BEGIN_FUNC:
                break;
            case OP_INTEGER:
                push_var ((double)op->op_val.ival);
                break;
            case OP_COUNTER:
                push_var ((double)(n + 1));
                break;
            case OP_INDEX:
                pop_var (&f1);
                index = (int)(f1 + TOLERANCE);
                if (index < 1)
                    longjmp (env_ptr, VECERR_BADINDX);
                if (op_index == op_pos)
                    longjmp (env_ptr, VECERR_POP);
                op = &op_start[op_index++];
                if (op->op != OP_VECTOR)
                     longjmp (env_ptr, VECERR_INDEX);
            case OP_VECTOR:
            {
                VECFLOAT *vec = op->op_val.vec;
                if (!index)
                    index = n + 1;
                if (--index < 0 || index >= (int)op->len)
                    longjmp (env_ptr, VECERR_BADINDX);
                push_var ((double)(vec[index]));
                index = 0;
                break;
            }
            case OP_NUMBER:
                push_var (op->op_val.num);
                break;

            /* unary operations */

            case OP_NEGATE:
                pop_var (&f1);
                push_var (-f1);
                break;
            case OP_NOT:
                pop_var (&f1);
                push_var ((double)(f1 == 0.0 ? 1 : 0));
                break;

            /* arithmetic operations */

            case OP_POWER:
                pop_var (&f2);
                pop_var (&f1);
                push_var (pow (f1, f2));
                break;
            case OP_MULTIPLY:
                pop_var (&f2);
                pop_var (&f1);
                push_var (f1 * f2);
                break;
            case OP_DIVIDE:
                pop_var (&f2);
                pop_var (&f1);
                if (f2 == 0.0)
                    longjmp (env_ptr, VECERR_DIVIDE);
                push_var (f1 / f2);
                break;
            case OP_MODULUS:
                pop_var (&f2);
                pop_var (&f1);
                push_var (fmod (f1, f2));
                break;
            case OP_PLUS:
                pop_var (&f2);
                pop_var (&f1);
                push_var (f1 + f2);
                break;
            case OP_MINUS:
                pop_var (&f2);
                pop_var (&f1);
                push_var (f1 - f2);
                break;

            /* comparisons */

            case OP_LESS:
                pop_var (&f2);
                pop_var (&f1);
                push_var ((double)(f1 < f2 ? 1 : 0));
                break;
            case OP_LESS_EQ:
                pop_var (&f2);
                pop_var (&f1);
                push_var ((double)(f1 <= f2 ? 1 : 0));
                break;
            case OP_GREATER:
                pop_var (&f2);
                pop_var (&f1);
                push_var ((double)(f1 > f2 ? 1 : 0));
                break;
            case OP_GREATER_EQ:
                pop_var (&f2);
                pop_var (&f1);
                push_var ((double)(f1 >= f2 ? 1 : 0));
                break;
            case OP_EQUALS:
                pop_var (&f2);
                pop_var (&f1);
                push_var ((double)(f1 == f2 ? 1 : 0));
                break;
            case OP_NOT_EQUAL:
                pop_var (&f2);
                pop_var (&f1);
                push_var ((double)(f1 != f2 ? 1 : 0));
                break;
            case OP_AND:
                pop_var (&f2);
                pop_var (&f1);
                push_var ((double)(f1 && f2 ? 1 : 0));
                break;
            case OP_OR:
                pop_var (&f2);
                pop_var (&f1);
                push_var ((double)(f1 || f2 ? 1 : 0));
                break;

            /* if-then operation */

            case OP_IF:
                level = ++nest_level;
                pop_var (&f1);
                if (f1 == 0.0) {    /* skip if(TRUE) operations */
                    while (1) {
                        op = &op_start[op_index++];
                        if (op->op == OP_ELSE) {
                            if (level == nest_level)
                                break;
                            level--;
                        }
                        else if (op->op == OP_IF)
                            level++;
                    }
                }
                push_var (process_op (n));
                nest_level--;
                break;
            case OP_ELSE:
                level = nest_level;
                while (1) {
                    op = &op_start[op_index++];
                    if (op->op == OP_ENDIF) {
                        if (level == nest_level)
                            break;
                        level--;
                    }
                    else if (op->op == OP_IF)
                        level++;
                }
            case OP_ENDIF:
                if (nest_level)
                    done = 1;
                break;

            /* functions which operate on the vectors */

            case FUNC_MAX:
            case FUNC_MIN:
            case FUNC_SUM:
            case FUNC_AVG:
            case FUNC_SIZE:
            case FUNC_LEN:

                /* if currently processing vector, return */

                if (vec_op) {
                    done = 1;
                    break;
                }

                /* else, first time here - find start of operations */
                /* also check that all vectors are the same length */

                index = op_index - 1;
                level = m = 0;
                while (op_start[--op_index].op != OP_BEGIN_FUNC) {
                    if (op_start[op_index].op == OP_VECTOR) {
                        if (!level)
                            level = op_start[op_index].len;
                        else if (level != op_start[op_index].len)
                            longjmp (env_ptr, VECERR_VECLEN);
                    }
                    else if (op_start[op_index].op == OP_COUNTER)
                        m = op_start[op_index].len;
                }
                if (!level && m)
                    level = m;

                /* set BEGIN_FUNC op to NOP and set start of processing */

                op_start[op_index].op = OP_NOP;
                vec_op = ++op_index;

                /* evaluate intrinsic function */

                pop_var (&f1);
                if (op->op == FUNC_LEN)
                    f1 = level ? level : 1;
                else {
                    for (m = 1; m < level; m++) {
                        op_index = vec_op;
                        f2 = process_op (m);
                        switch (op->op) {
                            case FUNC_MAX:
                                if (f1 < f2)
                                    f1 = f2;
                                break;
                            case FUNC_MIN:
                                if (f1 > f2)
                                    f1 = f2;
                                break;
                            case FUNC_AVG:
                            case FUNC_SUM:
                                f1 += f2;
                                break;
                            case FUNC_SIZE:
                                f1 += (f2 * f2);
                                break;
                        }
                    }
                    if (op->op == FUNC_AVG)
                        f1 /= (double)(level ? level : 1);
                    else if (op->op == FUNC_SIZE)
                        f1 = sqrt (f1);
                }

                /* set stack ops to NOP */

                for (m = vec_op; m < index; m++)
                    op_start[m].op = OP_NOP;
                op = &op_start[index];
                op_index = index + 1;
                vec_op = index = 0;

                /* put final result on stack */

                op->op = OP_NUMBER;
                op->op_val.num = f1;
                push_var (f1);
                break;

            /* function operations */

            case FUNC_ABS:
                pop_var (&f1);
                push_var (fabs (f1));
                break;
            case FUNC_ACOS:
                pop_var (&f1);
                push_var (acos (f1));
                break;
            case FUNC_ACOSH:
                pop_var (&f1);
                if (f1 < 1.0)
                    longjmp (env_ptr, VECERR_ACOSH);
                push_var (log (f1 + sqrt (f1 * f1 - 1.0)));
                break;
            case FUNC_ASIN:
                pop_var (&f1);
                push_var (asin (f1));
                break;
            case FUNC_ASINH:
                pop_var (&f1);
                push_var (log (f1 + sqrt (f1 * f1 + 1.0)));
                break;
            case FUNC_ATAN:
                pop_var (&f1);
                push_var (atan (f1));
                break;
            case FUNC_ATAN2:
                pop_var (&f2);
                pop_var (&f1);
                push_var (atan2 (f1, f2));
                break;
            case FUNC_ATANH:
                pop_var (&f1);
                if (f1 <= -1.0 || f1 >= 1.0)
                    longjmp (env_ptr, VECERR_ATANH);
                push_var (0.5 * log ((1.0 + f1) / (1.0 - f1)));
                break;
            case FUNC_CEIL:
                pop_var (&f1);
                push_var (ceil (f1));
                break;
            case FUNC_COSH:
                pop_var (&f1);
                push_var (cosh (f1));
                break;
            case FUNC_COS:
                pop_var (&f1);
                push_var (cos (f1));
                break;
            case FUNC_DEG:
                pop_var (&f1);
                push_var (f1 * 180.0 / PI);
                break;
            case FUNC_EXP:
                pop_var (&f1);
                push_var (exp (f1));
                break;
            case FUNC_FACT:
                pop_var (&f1);
                j = (int)f1;
                f2 = j;
                if (j == 0)
                    f2 = 1;
                else if (j > 0) {
                    while (--j)
                        f2 *= (double)j;
                }
                else {
                    while (++j)
                        f2 *= (double)j;
                }
                push_var (f2);
                break;
            case FUNC_FLOOR:
                pop_var (&f1);
                push_var (floor (f1));
                break;
            case FUNC_LOG10:
                pop_var (&f1);
                push_var (log10 (f1));
                break;
            case FUNC_LOG:
                pop_var (&f1);
                push_var (log (f1));
                break;
            case FUNC_MOD:
                pop_var (&f2);
                pop_var (&f1);
                push_var (fmod (f1, f2));
                break;
            case FUNC_POW:
                pop_var (&f2);
                pop_var (&f1);
                push_var (pow (f1, f2));
                break;
            case FUNC_POW10:
                pop_var (&f1);
                push_var (pow (10.0, f1));
                break;
            case FUNC_RAD:
                pop_var (&f1);
                push_var (f1 * PI / 180.0);
                break;
            case FUNC_RAND:
                push_var (vec_rand ());
                break;
            case FUNC_ROUND:
                pop_var (&f1);
                if (fmod (f1, 1.0) >= 0.5)
                    push_var (ceil (f1));
                else
                    push_var (floor (f1));
                break;
            case FUNC_SINH:
                pop_var (&f1);
                push_var (sinh (f1));
                break;
            case FUNC_SIN:
                pop_var (&f1);
                push_var (sin (f1));
                break;
            case FUNC_SQRT:
                pop_var (&f1);
                push_var (sqrt (f1));
                break;
            case FUNC_SQR:
                pop_var (&f1);
                push_var (f1 * f1);
                break;
            case FUNC_TANH:
                pop_var (&f1);
                push_var (tanh (f1));
                break;
            case FUNC_TAN:
                pop_var (&f1);
                push_var (tan (f1));
                break;

            default:
                longjmp (env_ptr, VECERR_OPERR);
        }
    }
    pop_var (&f1);
    return (f1);
}

/*===================================================================
 * math exception handling
 *===================================================================*/

#ifdef _WIN32
#include <float.h>
#endif

/*---------- fpe_err -----------------------------------------------
 * signal handler for arithmetic exceptions
 *------------------------------------------------------------------*/

static SIGNAL fpe_err (int signum)
{
    strcpy (err_msg[VECERR_MATH], "arithmetic exception (SIGFPE)");
#if defined(__MSDOS__) || defined(MSDOS) || defined(_WIN32)
    _fpreset ();
#endif
    longjmp (env_ptr, VECERR_MATH);
}

#ifdef USE_MATHERR

#ifdef _WIN32
#define matherr   _matherr
#define exception _exception
#endif

/*---------- matherr -----------------------------------------------
 * use this math error handler, instead of library
 * This handles errors generated by math function calls
 *------------------------------------------------------------------*/

int matherr (e)
struct exception *e;
{
    switch (e->type) {
        case DOMAIN:
            sprintf (err_msg[VECERR_MATH],
                "%s : argument domain error", e->name);
            break;
        case SING:
            sprintf (err_msg[VECERR_MATH],
                "%s : argument singularity", e->name);
            break;
        case OVERFLOW:
            sprintf (err_msg[VECERR_MATH],
                "%s : overflow range error", e->name);
            break;
        case UNDERFLOW:
            sprintf (err_msg[VECERR_MATH],
                "%s : underflow range error", e->name);
            break;
        case PLOSS:
            sprintf (err_msg[VECERR_MATH],
                "%s : partial loss of significance", e->name);
            break;
        case TLOSS:
            sprintf (err_msg[VECERR_MATH],
                "%s : total loss of significance", e->name);
            break;
        default:
            sprintf (err_msg[VECERR_MATH],
                "%s : unknown error", e->name);
            break;
    }
    longjmp (env_ptr, VECERR_MATH);
    return (0); /* quite compiler */
}

#endif  /* USE_MATHERR */

/*========== print_op ==============================================
 * debug print out of operation stack
 *==================================================================*/

#ifdef VEC_DEBUG

static char HEADER1[] = "op code                 pop   push     stack\n";
static char HEADER2[] = "-------                 ---   ----     -----\n";
static char FORMAT[]  = "%-20s%7d%7d%10d\n";
static char NUMFORM[] = "%-20.3g%7d%7d%10d\n";

void print_op (fp)
FILE *fp;
{
    int m, push, pop, stack_size = 0;
    char *opname;
    struct op_stack_ *op;

    if (op_size == 0)
        fprintf (fp, "OP STACK not created\n");
    else if (op_pos == 0)
        fprintf (fp, "OP STACK is empty\n");
    else {
        fprintf (fp, HEADER1);
        fprintf (fp, HEADER2);
        op_index = 0;

        while (op_index < op_pos) {
            op = &op_start[op_index++];
            push = pop = 1;
            switch (op->op) {
                case OP_NOP:
                    push = pop = 0;
                    opname = "nop";
                    break;
                case OP_BEGIN_FUNC:
                    push = pop = 0;
                    opname = "beginfunc";
                    break;
                case OP_INTEGER:
                    pop = 0;
                    opname = "integer";
                    break;
                case OP_COUNTER:
                    pop = 0;
                    opname = "counter";
                    break;
                case OP_INDEX:
                    push = 0;
                    opname = "index";
                    break;
                case OP_VECTOR:
                    pop = 0;
                    opname = "vector";
                    break;
                case OP_NUMBER:
                    pop = 0;
                    opname = "number";
                    break;

                /* unary operations */

                case OP_NEGATE:
                    opname = "negate";
                    break;
                case OP_NOT:
                    opname = "not";
                    break;

                /* arithmetic operations */

                case OP_POWER:
                    pop = 2;
                    opname = "power";
                    break;
                case OP_MULTIPLY:
                    pop = 2;
                    opname = "multiply";
                    break;
                case OP_DIVIDE:
                    pop = 2;
                    opname = "divide";
                    break;
                case OP_MODULUS:
                    pop = 2;
                    opname = "modulus";
                    break;
                case OP_PLUS:
                    pop = 2;
                    opname = "plus";
                    break;
                case OP_MINUS:
                    pop = 2;
                    opname = "minus";
                    break;

                /* comparisons */

                case OP_LESS:
                    pop = 2;
                    opname = "<";
                    break;
                case OP_LESS_EQ:
                    pop = 2;
                    opname = "<=";
                    break;
                case OP_GREATER:
                    pop = 2;
                    opname = ">";
                    break;
                case OP_GREATER_EQ:
                    pop = 2;
                    opname = ">=";
                    break;
                case OP_EQUALS:
                    pop = 2;
                    opname = "==";
                    break;
                case OP_NOT_EQUAL:
                    pop = 2;
                    opname = "!=";
                    break;
                case OP_AND:
                    pop = 2;
                    opname = "and";
                    break;
                case OP_OR:
                    pop = 2;
                    opname = "or";
                    break;

                /* if-then operation */

                case OP_IF:
                    opname = "if";
                    break;
                case OP_ELSE:
                    push = pop = 0;
                    opname = "else";
                    break;
                case OP_ENDIF:
                    push = pop = 0;
                    opname = "endif";
                    break;

                /* functions which operate on the vectors */

                case FUNC_MAX:
                    opname = "max()";
                    break;
                case FUNC_MIN:
                    opname = "min()";
                    break;
                case FUNC_SUM:
                    opname = "sum()";
                    break;
                case FUNC_AVG:
                    opname = "avg()";
                    break;
                case FUNC_SIZE:
                    opname = "size()";
                    break;
                case FUNC_LEN:
                    opname = "len()";
                    break;

                /* function operations */
            
                case FUNC_ABS:
                    opname = "abs()";
                    break;
                case FUNC_ACOS:
                    opname = "acos()";
                    break;
                case FUNC_ACOSH:
                    opname = "acosh()";
                    break;
                case FUNC_ASIN:
                    opname = "asin()";
                    break;
                case FUNC_ASINH:
                    opname = "asinh()";
                    break;
                case FUNC_ATAN:
                    opname = "atan()";
                    break;
                case FUNC_ATANH:
                    opname = "atanh()";
                    break;
                case FUNC_CEIL:
                    opname = "ceil()";
                    break;
                case FUNC_COSH:
                    opname = "cosh()";
                    break;
                case FUNC_COS:
                    opname = "cos()";
                    break;
                case FUNC_DEG:
                    opname = "deg()";
                    break;
                case FUNC_EXP:
                    opname = "exp()";
                    break;
                case FUNC_FACT:
                    opname = "fact()";
                    break;
                case FUNC_FLOOR:
                    opname = "floor()";
                    break;
                case FUNC_LOG10:
                    opname = "log10()";
                    break;
                case FUNC_LOG:
                    opname = "log()";
                    break;
                case FUNC_POW10:
                    opname = "pow10()";
                    break;
                case FUNC_RAD:
                    opname = "rad()";
                    break;
                case FUNC_RAND:
                    pop = 0;
                    opname = "rand()";
                    break;
                case FUNC_ROUND:
                    opname = "round()";
                    break;
                case FUNC_SINH:
                    opname = "sinh()";
                    break;
                case FUNC_SIN:
                    opname = "sin()";
                    break;
                case FUNC_SQRT:
                    opname = "sqrt()";
                    break;
                case FUNC_SQR:
                    opname = "sqr()";
                    break;
                case FUNC_TANH:
                    opname = "tanh()";
                    break;
                case FUNC_TAN:
                    opname = "tan()";
                    break;

                default:
                    fprintf (fp, "UNREGOGNIZED op code\n");
                    return;
            }
            stack_size += (push - pop);
            fprintf (fp, FORMAT, opname, pop, push, stack_size);
        }

        fprintf (fp, FORMAT, "RESULT", 1, 0, --stack_size);
    }

    if (veclist == NULL)
        fprintf (fp, "Vector list is clear\n");
    else {
        VECLIST *vl = veclist;
        m = 0;
        while (vl != NULL) {
            m++;
            vl = vl->prev;
        }
        fprintf (fp, "%d Vectors in vector list\n");
    }
}

#endif  /* VEC_DEBUG */
