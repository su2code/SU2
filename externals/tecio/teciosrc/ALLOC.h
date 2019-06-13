 #ifndef ALLOC_H
 #define ALLOC_H
#include "TASSERT.h"
 #if defined __cplusplus
#include <new>
 #endif
 #if !defined __cplusplus
 #define ___23(N,___4236,str) (___4236 *)malloc((N)*sizeof(___4236))
 #define ALLOC_ITEM(___4236,str)    (___4236 *)malloc(sizeof(___4236))
 #ifdef _DEBUG
 #define ___1530(X,str)  do { free((void *)(X)); *((void **)&(X)) = (void *)0xFFFF; } while (0)
 #define ___1531(X,str)   do { free((void *)(X)); *((void **)&(X)) = (void *)0xFFFF; } while (0)
 #else
 #define ___1530(X,str)  free((void *)(X))
 #define ___1531(X,str)   free((void *)(X))
 #endif
 #else
 #ifdef TRACK_MEMORY_USAGE
extern void ___1936(void); extern void ___491(void); extern void ___4207(size_t size); extern void ___4209(size_t size); extern void ___4208(void); extern void ___4210(void); extern void ___1758(size_t* ___2407, size_t* ___2406, size_t* ___2408, size_t* ___2409);
 #endif
 #if defined MSWIN && defined _DEBUG && defined TRACK_MEMORY_USAGE
template <typename T> inline T *nonExceptionNew(size_t      ___2812, const char* ___1394, int         lineNumber) { REQUIRE(___2812 > 0); REQUIRE(VALID_REF(___1394)); REQUIRE(lineNumber > 0); T* ___3358 = NULL; try {
 #ifdef DEBUG_NEW
 #ifdef new
 #undef new
 #define USING_DEBUG_NEW
 #endif
___3358 = new(___1394, lineNumber) T[___2812];
 #ifdef USING_DEBUG_NEW
 #undef USING_DEBUG_NEW
 #endif
 #else
___3358 = new T[___2812];
 #endif
} catch (std::bad_alloc&) { ___3358 = NULL; }
 #ifdef TRACK_MEMORY_USAGE
if (___3358 != NULL) {
 #ifdef MSWIN
___4207(_msize(___3358));
 #else
___4207(malloc_usable_size(___3358));
 #endif
}
 #endif
ENSURE(VALID_REF_OR_NULL(___3358)); return ___3358; }
 #define ___23(N,___4236,str) nonExceptionNew<___4236>((N),__FILE__,__LINE__)
 #else
template <typename T> inline T *nonExceptionNew(size_t ___2812) { REQUIRE(___2812 > 0); T *___3358 = NULL; try { ___3358 = new T[___2812]; } catch (std::bad_alloc&) { ___3358 = NULL; }
 #ifdef TRACK_MEMORY_USAGE
if (___3358 != NULL) {
 #ifdef MSWIN
___4207(_msize(___3358));
 #else
___4207(malloc_usable_size(___3358));
 #endif
}
 #endif
return ___3358; }
 #define ___23(N,___4236,str) nonExceptionNew<___4236>((N))
 #endif
 #define ALLOC_ITEM(___4236,str)    ___23(1,___4236,str)
template <typename T> inline void nonExceptionDelete(T* &___3251) { REQUIRE(VALID_REF(___3251));
 #if defined TRACK_MEMORY_USAGE
{ if (___3251 != NULL) {
 #ifdef MSWIN
___4209(_msize(___3251));
 #else
___4209(malloc_usable_size(___3251));
 #endif
} }
 #endif
delete [] ___3251;
 #if !defined NO_ASSERTS
___3251 = (T*)(void*)0xFFFF;
 #endif
} template <typename T> inline void nonExceptionDelete(T* const& ___3251) { nonExceptionDelete(const_cast<T*&>(___3251)); }
 #define ___1530(___3251,str)  nonExceptionDelete((___3251))
 #define ___1531(___3251,str)   ___1530(___3251,str)
 #endif
struct ___956 { template<typename T> void operator()(T*& object) { delete object; object = NULL; } };
 #endif 
