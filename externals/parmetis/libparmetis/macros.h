/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * macros.h
 *
 * This file contains macros used in multilevel
 *
 * Started 9/25/94
 * George
 *
 * $Id: macros.h 10578 2011-07-14 18:10:15Z karypis $
 *
 */


/* The following macro returns a random number in the specified range */
#define AND(a, b) ((a) < 0 ? ((-(a))&(b)) : ((a)&(b)))
#define OR(a, b) ((a) < 0 ? -((-(a))|(b)) : ((a)|(b)))
#define XOR(a, b) ((a) < 0 ? -((-(a))^(b)) : ((a)^(b)))

#define HASHFCT(key, size) ((key)%(size))

/* set/reset the current workspace core */
#define WCOREPUSH    do {PASSERT(ctrl,ctrl->mcore!=NULL); gk_mcorePush(ctrl->mcore);}while(0)
#define WCOREPOP     do {PASSERT(ctrl,ctrl->mcore!=NULL); gk_mcorePop(ctrl->mcore);}while(0)


/* Timer macros */
#define cleartimer(tmr) (tmr = 0.0)
#define starttimer(tmr) (tmr -= MPI_Wtime())
#define stoptimer(tmr) (tmr += MPI_Wtime())
#define gettimer(tmr) (tmr)
#define STARTTIMER(ctrl, tmr) \
   do { \
     IFSET((ctrl)->dbglvl, DBG_TIME, gkMPI_Barrier((ctrl)->gcomm));\
     IFSET((ctrl)->dbglvl, DBG_TIME, starttimer((tmr))); \
   } while (0)
#define STOPTIMER(ctrl, tmr) \
   do { \
     IFSET((ctrl)->dbglvl, DBG_TIME, gkMPI_Barrier((ctrl)->gcomm));\
     IFSET((ctrl)->dbglvl, DBG_TIME, stoptimer((tmr))); \
   } while (0)



/* Debugging macros */
#ifndef NDEBUG
#   define PASSERT(ctrl, expr)                                          \
    if (!(expr)) {                                               \
       myprintf(ctrl, "***ASSERTION failed on line %d of file %s: " #expr "\n", \
            __LINE__, __FILE__);                               \
       assert(expr);                                           \
    }

#   define PASSERTP(ctrl, expr, msg)                                          \
    if (!(expr)) {                                               \
        myprintf(ctrl, "***ASSERTION failed on line %d of file %s:" #expr "\n", \
              __LINE__, __FILE__);                               \
        myprintf msg ; \
        assert(expr); \
    }

#else
#   define PASSERT(ctrl, expr) ;
#   define PASSERTP(ctrl, expr,msg) ;
#endif 


/*************************************************************************
 * * These macros insert and remove nodes from the boundary list
 * **************************************************************************/
#define BNDInsert(nbnd, bndind, bndptr, vtx) \
   do { \
	bndind[nbnd] = vtx; \
	bndptr[vtx] = nbnd++;\
      } while(0)

#define BNDDelete(nbnd, bndind, bndptr, vtx) \
   do { \
        bndind[bndptr[vtx]] = bndind[--nbnd]; \
	bndptr[bndind[nbnd]] = bndptr[vtx]; \
	bndptr[vtx] = -1; \
      } while(0)


