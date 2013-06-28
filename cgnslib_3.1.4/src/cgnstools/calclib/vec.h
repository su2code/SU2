/*
 * vec.h
 */

#ifndef _VEC_H_
#define _VEC_H_

#define VEC_MAXCALLS    10  /* max recursive calls */

/*----- vector data structure -----*/

#define VEC_VALUE       0
#define VEC_VECTOR      1
#define VEC_COMPLEX     2

typedef double VECFLOAT;

typedef struct {
    int type;          /* type of data */
    size_t len;        /* length of vector */
    union {
        VECFLOAT val;  /* single data value */
        VECFLOAT *vec; /* vector data */
    } f;
} VECDATA;

/*--- function prototypes ---*/

#ifdef __cplusplus
extern "C" {
#endif

typedef VECDATA *(*VECFUNC)(
    int checking,   /* set if checking string instead of parsing */
    int nv,         /* number of arguments */
    VECDATA **v,    /* arguments */
    char **errmsg   /* returned error message */
);
typedef VECDATA *(*VECCALL)(
    int checking,   /* set if checking string instead of parsing */
    char **pos,     /* pointer to position in equation string */
    char **errmsg   /* returned error message */
);

VECDATA *vec_create (   /* create a vector data structure */
    int type,           /* type of data */
    size_t len,         /* length of vector */
    int temp            /* set if temporary */
);

void vec_destroy (      /* destroy a vector data structure */
    VECDATA *vdata      /* structure to destroy */
);

int vec_add (           /* add vector to temporary list */
    VECDATA *vdata      /* data to add */
);

void vec_remove (       /* remove vector from temporary list */
    VECDATA *vdata      /* data to remove */
);

void vec_free (         /* free internal memory */
    void
);

void vec_maxcalls (     /* sets maximum number recursive calls */
    int maxcalls        /* max recursive calls */
);

char *vec_errmsg (      /* returns error messages */
    void
);

char **vec_list (       /* return list of intrinsics */
    void
);

VECDATA *vec_parse (    /* parse and process equation */
    char *str,          /* equation string */
    int min_len,        /* minimum length of vector for counter */
    VECCALL func        /* user call-back for unknown symbols */
);

size_t vec_check (      /* parse equation and check for errors */
    char *str,          /* equation string */
    int min_len,        /* minimum length of vector for counter */
    VECCALL func        /* user call-back for unknown symbols */
);

int vec_number (        /* get floating point number from string */
    double *dval,       /* returned number */
    char **sp           /* updated string pointer */
);

void vec_randinit (     /* initialize random number generator */
    int seed            /* initialization seed */
);

double vec_rand (       /* return random floating point number */
    void
);

#ifdef __cplusplus
}
#endif

#endif      /* _VEC_H_ */
