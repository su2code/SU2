#ifndef __CG_MALLOC_H__
#define __CG_MALLOC_H__

#define malloc(S)    cgmalloc(S)
#define realloc(P,S) cgrealloc(P,S)
#define calloc(N,S)  cgcalloc(N,S)
#define free(P)      cgfree(P)

extern void  *cgmalloc(size_t);
extern void  *cgrealloc(void *,size_t);
extern void  *cgcalloc(size_t,size_t);
extern void   cgfree(void *);

extern size_t cgmemnow();
extern size_t cgmemmax();
extern size_t cgalloccalls();
extern size_t cgfreecalls();

#endif
