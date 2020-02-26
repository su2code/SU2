 #ifndef TASSERT_H
 #define TASSERT_H
 #if !defined TECPLOTKERNEL && !defined STD_ASSERTS  && !defined CHECKED_BUILD
 #define STD_ASSERTS
 #endif
#  include <assert.h>
 #  if !defined ASSERT
 #    define ASSERT assert
 #  endif
 #if defined MSWIN
 #if defined ENSURE
 #undef ENSURE
 #endif 
 #endif 
 #define INVALID_REF       ((void *)0x0000FFFF)
 #define UNINITIALIZED_REF ((void *)0x11111111)
 #define INVALID_FN_REF    ((void *)NULL)
 #ifdef UNIXX
 #if !defined TECPLOTKERNEL
 #  if !defined VALID_REF
 #    define VALID_REF(p)      ( (p)  != NULL )
 #  endif
 #  if !defined VALID_FN_REF
 #    define VALID_FN_REF(___1481)  ( (___1481) != NULL )
 #  endif
 #endif 
 # define VALID_WIDGET(widget)       VALID_REF((widget))
 # define VALID_MENU_WIDGET(widget)   VALID_REF((widget))
 #endif 
 #ifdef MSWIN
 #if defined NO_ASSERTS
 #  if !defined VALID_REF
 #    define VALID_REF(p)      ((p)  != NULL)
 #  endif
 #  if !defined VALID_FN_REF
 #    define VALID_FN_REF(___3002)  ((___3002) != NULL)
 #  endif
 #else
 #  if !defined VALID_REF
 #    define VALID_REF(p)      ((p)  != NULL && !IsBadReadPtr((const void *)(p), 1))
 #  endif
 #  if !defined VALID_FN_REF
 #    define VALID_FN_REF(___3002)  ((___3002) != NULL && !IsBadReadPtr((const void *)(___3002),(UINT_PTR)sizeof(const void*)))
 #  endif
 #endif
 #   define VALID_WIDGET(widget)       ((widget) != NULL)
 # define VALID_MENU_WIDGET(widget)  ((widget)!=NULL)
 #endif 
 #define VALID_HANDLE(handle)       ((handle)!=0)
 #if !defined VALID_BOOLEAN
 #  define VALID_BOOLEAN(b)           ((b) == ___4226 || (b) == ___1305)
 #endif
 #if !defined VALID_ENUM
 #  define VALID_ENUM(___4314, type)    (0 <= (int)(___4314) && \
 (int)(___4314) < END_##type)
 #endif
 #if !defined VALID_REF_OR_NULL
 #  define VALID_REF_OR_NULL(___3251) IMPLICATION((___3251) != NULL, VALID_REF(___3251))
 #endif
 #if !defined VALID_FN_REF_OR_NULL
 #  define VALID_FN_REF_OR_NULL(___3251) IMPLICATION((___3251) != NULL, VALID_FN_REF(___3251))
 #endif
 #define VALID_TRANSLATED_STRING(___4228) (!(___4228).___2035())
struct ___802; namespace tecplot { class ___2090; } bool VALID_FE_CLASSIC_CELL_INDEX( ___802 const* ___800, ___2227       ___462); bool VALID_FE_CELL_INDEX( ___802 const* ___800, ___2227       ___462); bool VALID_FE_CELL_INDEX( ___802 const*             ___800, tecplot::___2090 const& ___451);
 #define VALID_IPLANE_CELL_INDEX(___801,___463) \
 (  \
 (___463) >= 0 && \
 ___1842((___801),___463) <= MAX((___801)->___2811,1) && \
 ___2112((___801),___463) <  MAX((___801)->___2816,1) && \
 ___2157((___801),___463) <  MAX((___801)->___2819,1))
 #define VALID_JPLANE_CELL_INDEX(___801,___463) \
 (  \
 (___463) >= 0 && \
 ___1842((___801),___463) <  MAX((___801)->___2811,1) && \
 ___2112((___801),___463) <= MAX((___801)->___2816,1) && \
 ___2157((___801),___463) <  MAX((___801)->___2819,1))
 #define VALID_KPLANE_CELL_INDEX(___801,___463) \
 (  \
 (___463) >= 0 && \
 ___1842((___801),___463) <  MAX((___801)->___2811,1) && \
 ___2112((___801),___463) <  MAX((___801)->___2816,1) && \
 ___2157((___801),___463) <= MAX((___801)->___2819,1))
 #define VALID_ORDERED_CELL_INDEX(___801, ___463, ___3095) \
 (  \
 ((IJKPlanes_e)(___3095) == ___1867 || \
 (IJKPlanes_e)(___3095) == ___1872 || \
 (IJKPlanes_e)(___3095) == ___1874 || \
 (IJKPlanes_e)(___3095) == ___1876) && \
 \
   \
 (IMPLICATION(((IJKPlanes_e)(___3095) == ___1867 || \
 (IJKPlanes_e)(___3095) == ___1876), \
 VALID_IPLANE_CELL_INDEX((___801),___463)) && \
 IMPLICATION(((IJKPlanes_e)(___3095) == ___1872 || \
 (IJKPlanes_e)(___3095) == ___1876), \
 VALID_JPLANE_CELL_INDEX((___801),___463)) && \
 IMPLICATION(((IJKPlanes_e)(___3095) == ___1874 || \
 (IJKPlanes_e)(___3095) == ___1876), \
 VALID_KPLANE_CELL_INDEX((___801),___463))))
bool VALID_CELL_INDEX( ___802 const* ___800, ___2227       ___462, IJKPlanes_e     ___1865); bool VALID_CELL_INDEX( ___802 const*             ___800, tecplot::___2090 const& ___451, IJKPlanes_e                 ___1865);
 #define VALID_DATASET(___882,___484) (((___882) != NULL) && \
 IMPLICATION((___484),(___882)->___2847 >= 1))
 #ifdef MSWIN
 # define VALID_WND(___3257) (::___2083((___3257)->___1771()))
 #else 
 # define VALID_WND(___3257) 
 #endif 
 #if !defined VALID_NON_ZERO_LEN_STR
 #  if defined MSWIN
 #    if defined NO_ASSERTS
 #      define VALID_NON_ZERO_LEN_STR(str) (VALID_REF(str) && !___2017(str))
 #    else
 #      define VALID_NON_ZERO_LEN_STR(str) \
 (VALID_REF(str)                                                            && \
 !IsBadReadPtr((const void*)(str),(UINT_PTR)(1+strlen((const char*)(str)))) && \
 !___2017(str))
 #    endif
 #  else
 #    define VALID_NON_ZERO_LEN_STR(str) (VALID_REF(str) && !___2017(str))
 #  endif
 #endif
 #if !defined VALID_SET_INDEX
 #  define VALID_SET_INDEX(___3492) (((___3493)___3492)>=(___3493)1)
 #endif
 #if !defined VALID_FILE_HANDLE
 #  define VALID_FILE_HANDLE(___3792) ((___3792) != NULL)
 #endif
 #define VALID_BASIC_COLOR(___351) \
 (___1420<=(___351) && (___351)<=___2195)
 #define VALID_CONTOUR_COLOR(Color) \
 (___614<=(Color) && \
 (Color)<___614+___1547.___2241.___2379+1)
 #define VALID_PLOTTING_COLOR(Color) \
 (VALID_BASIC_COLOR(Color) || VALID_CONTOUR_COLOR(Color))
 #define VALID_INTERFACE_SPECIFIC_COLOR(___351) \
 (___1423<=(___351) && (___351)<=___2200)
 #define VALID_INTERFACE_COLOR(Color) \
 (VALID_PLOTTING_COLOR(Color) || VALID_INTERFACE_SPECIFIC_COLOR(Color))
 #define VALID_MULTICOLOR_COLOR(Color) \
 (((Color) == ___2662) || ((Color) == ___2655) || \
 ((Color) == ___2656) || ((Color) == ___2657) || \
 ((Color) == ___2658) || ((Color) == ___2659) || \
 ((Color) == ___2660) || ((Color) == ___2661))
 #define VALID_RGB_COLOR(Color) \
 ((Color) == ___3375)
 #define VALID_ASSIGNABLE_COLOR(C) \
 (VALID_BASIC_COLOR(C)      || \
 VALID_MULTICOLOR_COLOR(C) || \
 VALID_RGB_COLOR(C))
 #define VALID_PEN_OFFSET(___3000) \
 (___364<=(___3000) && (___3000)<=___2826)
 #define VALID_PEN_OFFSET_FOR_OBJECT(___3000) \
 (___1424<=(___3000) && (___3000)<=___2202)
 #define VALID_NAME(___2686, ___2376) \
 (VALID_REF(___2686) && \
 (___2017(___2686) || \
 (!tecplot::isspace((___2686)[0]) && !tecplot::isspace((___2686)[strlen(___2686)-1]))) && \
 strlen(___2686) <= (___2376))
 #define VALID_ZONE_NAME(___2686) VALID_NAME((___2686), ___2358)
 #define VALID_VAR_NAME(___2686)  VALID_NAME((___2686), ___2356)
 #define VALID_LIGHTINGEFFECT(___2165) \
 (((___2165) == ___2239) || ((___2165) == ___2236))
typedef void (*TAssertFailureNotifyFunc)( const char *___1246, const char *___1396, int        line);
 #if !defined STD_ASSERTS
extern void TAssert( const char *___1246, const char *___1396, int        line); extern TAssertFailureNotifyFunc ___1957( TAssertFailureNotifyFunc ___2698);
 #endif 
 #if defined NO_ASSERTS
 #   define TASSERT(___1245)
 #   if !defined INVARIANT
 #     define INVARIANT(___1245)
 #   endif
 #   if !defined REQUIRE
 #     define REQUIRE(___1245)
 #   endif
 #   if !defined ENSURE
 #     define ENSURE(___1245)
 #   endif
 #   if !defined ___478
 #     define ___478(___1245)
 #   endif
 #   ifdef VERIFY
 #     undef VERIFY
 #   endif
 #   define VERIFY(___1245)    ((void)(___1245))
 #   if !defined ASSERT_ONLY
 #     define ASSERT_ONLY(___2972)
 #   endif
 #   if !defined NOT_IMPLEMENTED
 #     if defined ___1840
 #       define NOT_IMPLEMENTED() ___478(___1305)
 #     else
 #       if defined MSWIN
 #         define NOT_IMPLEMENTED(x)  TAssert("Not Implemented", __FILE__, __LINE__)
 #       endif
 #       if defined UNIXX
 #         define NOT_IMPLEMENTED()  not ___1907 
 #       endif
 #     endif
 #   endif
 #elif defined STD_ASSERTS
 #   define TASSERT(___1245)         assert(___1245)
 #   if !defined INVARIANT
 #     define INVARIANT(___1245)       assert(___1245)
 #   endif
 #   if !defined REQUIRE
 #     define REQUIRE(___1245)         assert(___1245)
 #   endif
 #   if !defined ENSURE
 #     define ENSURE(___1245)          assert(___1245)
 #   endif
 #   if !defined ___478
 #     define ___478(___1245)           assert(___1245)
 #   endif
 #   ifdef VERIFY
 #     undef VERIFY
 #   endif
 #   ifndef VERIFY
 #     if defined NDEBUG
 #       define VERIFY(___1245) ((void)(___1245))
 #     else
 #       define VERIFY(___1245) assert(___1245)
 #     endif
 #   endif 
 #   if !defined NOT_IMPLEMENTED
 #     define NOT_IMPLEMENTED()     assert(!("Not Implemented"))
 #   endif
 #   if !defined ASSERT_ONLY
 #     define ASSERT_ONLY(___2972) ___2972
 #   endif
 #else
 #if defined (MSWIN)
 #if defined CHECKED_BUILD
#include <string>
#include <vector>
#include <algorithm>
class ___212 { public: static void ___1934(const std::string &___1394); static bool ___2031(); static void ___5(const std::string &___2432); private: static void ___4540(const std::string &___2432); private: static bool                     ___2317; static std::string              ___2318; static std::vector<std::string> ___213; };
 #define TASSERT(___1245)\
 do { if (!(___1245)) {   } } while (0)
 #else
 #define TASSERT(___1245) ASSERT(___1245)
 #endif 
 #else
 #define TASSERT(___1245) (void)((___1245) || (TAssert(#___1245, __FILE__, __LINE__), 0))
 #endif
 #   if !defined INVARIANT
 #   if defined NO_INVARIANTS
 #     define INVARIANT(___1245)
 #   else
 #     define INVARIANT(___1245) TASSERT(___1245)
 #   endif
 #   endif
 #   if !defined REQUIRE
 #   if defined ___2753
 #     define REQUIRE(___1245)
 #   else
 #     define REQUIRE(___1245) TASSERT(___1245)
 #   endif
 #   endif
 #   if !defined ENSURE
 #   if defined ___2752
 #     define ENSURE(___1245)
 #   else
 #     define ENSURE(___1245) TASSERT(___1245)
 #   endif
 #   endif
 #   if !defined ___478
 #   if defined NO_CHECKS
 #     define ___478(___1245)
 #   else
 #     define ___478(___1245)  TASSERT(___1245)
 #   endif
 #   endif
 #   if !defined VERIFY
 #   if defined NO_CHECKS
 #     define VERIFY(___1245)  ((void)(___1245))
 #   else
 #     if defined NDEBUG
 #       define VERIFY(___1245) ((void)(___1245))
 #     else
 #       define VERIFY(___1245) TASSERT(___1245)
 #     endif
 #   endif
 #   endif
 #   if defined NICE_NOT_IMPLEMENTED
 #     define NOT_IMPLEMENTED() ___2706()
 #   else
 #     define NOT_IMPLEMENTED() TASSERT(!("Not Implemented"))
 #   endif
 #   if !defined ASSERT_ONLY
 #     define ASSERT_ONLY(___2972) ___2972
 #   endif
 #endif
 #if !defined ASSERT_ONLY_PARAM && defined ASSERT_ONLY
 #   define ASSERT_ONLY_PARAM(___2972) ASSERT_ONLY(___2972)
 #endif
 #if !defined STD_ASSERTS
extern void ___4027(void);
 #endif 
 #if defined NICE_NOT_IMPLEMENTED
extern void ___2706(void);
 #endif
 #if !defined IMPLICATION
 #  define IMPLICATION(___2894,___3258) (!(___2894) || (___3258))
 #endif
 #if !defined EQUIVALENCE
 #  define EQUIVALENCE(___2894,___3258) ((___2894) == (___3258))
 #endif
 #endif 
