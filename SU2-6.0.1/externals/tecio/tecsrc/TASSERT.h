/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
/*
 * Provide four levels of assertion control. Assertions provide a mechanism
 * to enforce a contract between a client and service provider. The assertions
 * are listed in order of highest to lowest priority. Assertions can be turned
 * off individually by defining the appropriate name (see preprossessor
 * definitions below), however, lower priority assertions should be turned
 * off prior to higher ones. As confidence in the code increases all assertions
 * can be turned off by defining NO_ASSERTS.
 *
 * The assertions defined below have the following meanings:
 *
 *     INVARIANT - Asserts that a property's state is invariant throughout the
 *                 life of the property's scope. Stating invariant properties
 *                 of an application provides a deeper understanding of the
 *                 application's state.  These statements are usually
 *                 positioned just ahead of the preconditions and just after
 *                 the postconditions.
 *
 *     REQUIRE   - Asserts that a method's preconditions are within their
 *                 valid domains. Preconditions are conditions placed upon
 *                 any state information relied upon for the call. These
 *                 statements should be as close to the top of the method
 *                 as possible (except for assertions on invariant properties).
 *
 *     ENSURE    - Asserts that a method's postconditions are within their
 *                 valid ranges. Postconditions are conditions placed upon
 *                 any state information modified by the call. These
 *                 statements should be as close to the bottom of the method
 *                 (presumably there is only one exit point) as possible
 *                 (except for assertions on invariant properties).
 *
 *     CHECK     - Any other assertion not covered by the above assertions.
 *                 These are often added within a method body to specify
 *                 something that may not be immediately obvious to the reader
 *                 or to validate your assumptions about a call to a 3rd party
 *                 method that does not use runtime assertions for its
 *                 preconditions or postconditions. Obviously if the 3rd party
 *                 method uses assertions then there is no need for the CHECK.
 *
 * Additionally a convenience macro is available to place in code that is
 * pending implementation.
 *
 *     NOT_IMPLEMENTED - Assertion that always fails during runtime for debug
 *                       builds and always fails at compile time for release
 *                       builds.
 */
#if !defined TASSERT_H
#define TASSERT_H

#if defined (MSWIN)
# include <assert.h>
#endif /* MSWIN */

#if !defined TECPLOTKERNEL && !defined STD_ASSERTS
#define STD_ASSERTS
#endif

#if !defined (MSWIN)
#  include <assert.h>
#  if !defined ASSERT
#    define ASSERT assert
#  endif
#endif

#if defined MSWIN
/* MFC .NET defines ENSURE, so we undefine it here */
#if defined ENSURE
#undef ENSURE
#endif /* ENSURE */
#endif /* MSWIN */

/* BEGINREMOVEFROMADDON */
#define INVALID_REF       ((void *)0x0000FFFF)
/*
 * Chances are low the address 0x11111111 will be used, so we'll risk asserting
 * against it (see unitialized assignment in newmalloc).
 */
#define UNINITIALIZED_REF ((void *)0x11111111)
#define INVALID_FN_REF    ((void *)NULL)
/* ENDREMOVEFROMADDON */

#ifdef UNIXX
/* BEGINREMOVEFROMADDON */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#  if defined NO_ASSERTS
#  else
#  endif
#endif /* TECPLOTKERNAL */
/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
#if !defined TECPLOTKERNEL
/* For add-ons, there is a problem with VALID_REF, so just test for non-NULL */
/* ENDREMOVEFROMADDON */
#  if !defined VALID_REF
#    define VALID_REF(p)      ( (p)  != NULL )
#  endif
#  if !defined VALID_FN_REF
#    define VALID_FN_REF(fp)  ( (fp) != NULL )
#  endif
/* BEGINREMOVEFROMADDON */
#endif /* !defined TECPLOTKERNAL */
/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
/* Widgets are pointers under Motif */
# define VALID_WIDGET(widget)       VALID_REF((widget))
/* Menu widgets are pointers too */
# define VALID_MENU_WIDGET(widget)   VALID_REF((widget))
/* ENDREMOVEFROMADDON */
#endif /* UNIXX */

#ifdef MSWIN
/* BEGINREMOVEFROMADDON */
/* Don't use AfxIsValidAddress()! See Bug <7245>.
   1/4/08, dto. */
/* ENDREMOVEFROMADDON */

#if defined NO_ASSERTS
/* release build in TecUtil layer uses these for TUASSERT */
#  if !defined VALID_REF
#    define VALID_REF(p)      ((p)  != NULL)
#  endif
#  if !defined VALID_FN_REF
#    define VALID_FN_REF(pf)  ((pf) != NULL)
#  endif
#else
#  if !defined VALID_REF
#    define VALID_REF(p)      ((p)  != NULL && !IsBadReadPtr((const void *)(p), 1))
#  endif
#  if !defined VALID_FN_REF
#    define VALID_FN_REF(pf)  ((pf) != NULL && !IsBadReadPtr((const void *)(pf),(UINT_PTR)sizeof(const void*)))
#  endif
#endif

/* BEGINREMOVEFROMADDON */
/* Widgets are numbers under Windows, so we decode it with GetWindowFromWidget */
# if defined ENGINE
#   define VALID_WIDGET(widget)       ((widget) != NULL)
# else
#   define VALID_WIDGET(widget)       ((widget) != NULL && GetWindowFromWidget((widget))!=NULL)
# endif // ENGINE

/* Menu widgets are numbers too, so we just check against zero */
# define VALID_MENU_WIDGET(widget)  ((widget)!=NULL)
/* ENDREMOVEFROMADDON */
#endif /* MSWIN */
/* BEGINREMOVEFROMADDON */
/* handles are not pointers to memory, so the only test we can */
/* perform is to check for 0 */
#define VALID_HANDLE(handle)       ((handle)!=0)

#if defined FLEXLM
#define VALID_FLEX_JOB_HANDLE(handle) ((handle) != NULL)
#define VALID_FLEX_ERROR_CODE(ErrorCode)(ErrorCode <= 0)
#endif /* FLEXLM */

/* ENDREMOVEFROMADDON */
/* other useful validity checks */
#if !defined VALID_BOOLEAN
#  define VALID_BOOLEAN(b)           ((b) == TRUE || (b) == FALSE)
#endif
#if !defined VALID_ENUM
#  define VALID_ENUM(value, type)    (0 <= (int)(value) && \
                                           (int)(value) < END_##type)
#endif

/* Test a parameter than can be NULL or a valid pointer */
#if !defined VALID_REF_OR_NULL
#  define VALID_REF_OR_NULL(ptr) IMPLICATION((ptr) != NULL, VALID_REF(ptr))
#endif
#if !defined VALID_FN_REF_OR_NULL
#  define VALID_FN_REF_OR_NULL(ptr) IMPLICATION((ptr) != NULL, VALID_FN_REF(ptr))
#endif

/* BEGINREMOVEFROMADDON */
#define VALID_TRANSLATED_STRING(ts) (!(ts).isNull())

/**
 * These macros are a little complicated but it allows one to
 * write a simple assertion regardless of the zone type or
 * selected plane:
 *
 *   REQUIRE(VALID_CELL_INDEX(CZData, CellIndex, Plane)));
 *
 * Prior to using the macros a call to SetupXxx,
 * or at a minimum SetupCZData, must be called to setup
 * the globals defining the dataset structure.
 */
#define VALID_FE_CELL_INDEX(CZData, CellIndex) \
            (/* CellIndex range test */ \
             0 <= (CellIndex) && \
                  (CellIndex) < (CZData)->NumElements)

#define VALID_IPLANE_CELL_INDEX(CZData,CellIndex) \
            (/* CellIndex range test */ \
             (CellIndex) >= 0 && \
             IINDEX((CZData),CellIndex) <= MAX((CZData)->NumIPtsM1,1) && \
             JINDEX((CZData),CellIndex) <  MAX((CZData)->NumJPtsM1,1) && \
             KINDEX((CZData),CellIndex) <  MAX((CZData)->NumKPtsM1,1))

#define VALID_JPLANE_CELL_INDEX(CZData,CellIndex) \
            (/* CellIndex range test */ \
             (CellIndex) >= 0 && \
             IINDEX((CZData),CellIndex) <  MAX((CZData)->NumIPtsM1,1) && \
             JINDEX((CZData),CellIndex) <= MAX((CZData)->NumJPtsM1,1) && \
             KINDEX((CZData),CellIndex) <  MAX((CZData)->NumKPtsM1,1))

#define VALID_KPLANE_CELL_INDEX(CZData,CellIndex) \
            (/* CellIndex range test */ \
             (CellIndex) >= 0 && \
             IINDEX((CZData),CellIndex) <  MAX((CZData)->NumIPtsM1,1) && \
             JINDEX((CZData),CellIndex) <  MAX((CZData)->NumJPtsM1,1) && \
             KINDEX((CZData),CellIndex) <= MAX((CZData)->NumKPtsM1,1))

#define VALID_ORDERED_CELL_INDEX(CZData, CellIndex, Plane) \
            (/* macro preconditions */ \
             ((IJKPlanes_e)(Plane) == IJKPlanes_I || \
              (IJKPlanes_e)(Plane) == IJKPlanes_J || \
              (IJKPlanes_e)(Plane) == IJKPlanes_K || \
              (IJKPlanes_e)(Plane) == IJKPlanes_Volume) && \
\
             /* CellIndex range test */ \
             (IMPLICATION(((IJKPlanes_e)(Plane) == IJKPlanes_I || \
                           (IJKPlanes_e)(Plane) == IJKPlanes_Volume), \
                          VALID_IPLANE_CELL_INDEX((CZData),CellIndex)) && \
              IMPLICATION(((IJKPlanes_e)(Plane) == IJKPlanes_J || \
                           (IJKPlanes_e)(Plane) == IJKPlanes_Volume), \
                          VALID_JPLANE_CELL_INDEX((CZData),CellIndex)) && \
              IMPLICATION(((IJKPlanes_e)(Plane) == IJKPlanes_K || \
                           (IJKPlanes_e)(Plane) == IJKPlanes_Volume), \
                          VALID_KPLANE_CELL_INDEX((CZData),CellIndex))))

#define VALID_CELL_INDEX(CZData, CellIndex, Plane) \
           (((CZData)->NM != NULL || (CZData)->FM != NULL) ? \
              VALID_FE_CELL_INDEX((CZData), (CellIndex)) : \
              VALID_ORDERED_CELL_INDEX((CZData), (CellIndex), (Plane)))

#define VALID_DATASET(dataSet,checkNumZones) (((dataSet) != NULL) && \
                                              IMPLICATION((checkNumZones),(dataSet)->NumZones >= 1))



#ifdef MSWIN
/* Here is a more specific check in Windows for a valid
   pointer to an MFC Window object.
   Note that GetSafeHwnd() works even if pWnd is NULL, because
   it checks the 'this' pointer first */
# define VALID_WND(pWnd) (::IsWindow((pWnd)->GetSafeHwnd()))

#else /* !MSWIN */
# define VALID_WND(pWnd) /* Should not be used in Motif */
#endif /* MSWIN */
/* ENDREMOVEFROMADDON */

/* Check for a non-zero length string */
#if !defined VALID_NON_ZERO_LEN_STR
#  if defined MSWIN
#    if defined NO_ASSERTS
#      define VALID_NON_ZERO_LEN_STR(str) (VALID_REF(str) && !ISEMPTYSTRING(str))
#    else
#      define VALID_NON_ZERO_LEN_STR(str) \
           (VALID_REF(str)                                                            && \
           !IsBadReadPtr((const void*)(str),(UINT_PTR)(1+strlen((const char*)(str)))) && \
           !ISEMPTYSTRING(str))
#    endif
#  else
#    define VALID_NON_ZERO_LEN_STR(str) (VALID_REF(str) && !ISEMPTYSTRING(str))
#  endif
#endif

#if !defined VALID_SET_INDEX
#  define VALID_SET_INDEX(setIndex) (((SetIndex_t)setIndex)>=(SetIndex_t)1)
#endif

/* Check for valid stdio file handle */
#if !defined VALID_FILE_HANDLE
#  define VALID_FILE_HANDLE(stream) ((stream) != NULL)
#endif

/* To check colors and pen numbers */
/* BEGINREMOVEFROMADDON */
#define VALID_BASIC_COLOR(BColor) \
          (FirstBasicColor<=(BColor) && (BColor)<=LastBasicColor)
#define VALID_CONTOUR_COLOR(Color) \
          (ContourColorOffset<=(Color) && \
           (Color)<ContourColorOffset+GeneralBase.Limits.MaxNumContourLevels+1)
#define VALID_PLOTTING_COLOR(Color) \
          (VALID_BASIC_COLOR(Color) || VALID_CONTOUR_COLOR(Color))
#define VALID_INTERFACE_SPECIFIC_COLOR(BColor) \
          (FirstInterfaceColor<=(BColor) && (BColor)<=LastInterfaceColor)
#define VALID_INTERFACE_COLOR(Color) \
          (VALID_PLOTTING_COLOR(Color) || VALID_INTERFACE_SPECIFIC_COLOR(Color))
#define VALID_MULTICOLOR_COLOR(Color) \
          (((Color) == MultiColor_C) || ((Color) == MultiColor2_C) || \
           ((Color) == MultiColor3_C) || ((Color) == MultiColor4_C) || \
           ((Color) == MultiColor5_C) || ((Color) == MultiColor6_C) || \
           ((Color) == MultiColor7_C) || ((Color) == MultiColor8_C))
#define VALID_RGB_COLOR(Color) \
          ((Color) == RGBColor_C)
#define VALID_ASSIGNABLE_COLOR(C) \
        (VALID_BASIC_COLOR(C)      || \
         VALID_MULTICOLOR_COLOR(C) || \
         VALID_RGB_COLOR(C))
#define VALID_PEN_OFFSET(PenOffset) \
          (Black_C<=(PenOffset) && (PenOffset)<=NumPlotterPens)
#define VALID_PEN_OFFSET_FOR_OBJECT(PenOffset) \
          (FirstObjectPen<=(PenOffset) && (PenOffset)<=LastObjectPen)


/* to check FE cells */
#define VALID_ELEMENT_TYPE(element_type) \
          ((element_type) == ZoneType_FETriangle || \
           (element_type) == ZoneType_FEQuad     || \
           (element_type) == ZoneType_FETetra    || \
           (element_type) == ZoneType_FEBrick    || \
           (element_type) == ZoneType_FELineSeg)



/*
 * Test validity of zone and variable names. A valid name is one that has a
 * valid reference, is not padded with spaces and is within the maximum
 * specified length.
 */
#define VALID_NAME(Name, MaxLength) \
          (VALID_REF(Name) && \
           (ISEMPTYSTRING(Name) || \
            (!tecplot::isspace((Name)[0]) && !tecplot::isspace((Name)[strlen(Name)-1]))) && \
           strlen(Name) <= (MaxLength))
#define VALID_ZONE_NAME(Name) VALID_NAME((Name), MaxChrsZnTitle)
#define VALID_VAR_NAME(Name)  VALID_NAME((Name), MaxChrsVarName)


/* Special test for lighting effect (don't allow "none" in some cases) */
#define VALID_LIGHTINGEFFECT(L) \
          (((L) == LightingEffect_Paneled) || ((L) == LightingEffect_Gouraud))


/* type definition for assert failure notification function */
typedef void (*TAssertFailureNotifyFunc)(
    const char *expression, /* text representation of the assertion */
    const char *file_name,  /* name of the file containing the assertion */
    int        line);       /* line number in the file of the assertion */

#if !defined STD_ASSERTS
/* external function prototypes */
extern void TAssert(
    const char *expression, /* text representation of the assertion */
    const char *file_name,  /* name of the file containing the assertion */
    int        line);       /* line number in the file of the assertion */

extern TAssertFailureNotifyFunc InstallTAssertFailureNotify(
    TAssertFailureNotifyFunc new_function); /* new notification function */
#endif /* !STD_ASSERTS */
/* ENDREMOVEFROMADDON */

#if defined NO_ASSERTS
/* BEGINREMOVEFROMADDON */
#   define TASSERT(EXPR)
/* ENDREMOVEFROMADDON */
#   if !defined INVARIANT
#     define INVARIANT(EXPR)
#   endif
#   if !defined REQUIRE
#     define REQUIRE(EXPR)
#   endif
#   if !defined ENSURE
#     define ENSURE(EXPR)
#   endif
#   if !defined CHECK
#     define CHECK(EXPR)
#   endif
#   ifdef VERIFY
#     undef VERIFY
#   endif
#   define VERIFY(EXPR)    ((void)(EXPR))

/*
 * If a project is compiled with "warnings treated as errors" we need a way to
 * to remove parameters that are only used for assertions if assertions are
 * turned off.
 *
 * This macro is used in the implementation as follows:
 *
 *     void someFunction(int const ASSERT_ONLY_PARAM(paramOnlyUsedInAssertions))
 *     {
 *        ...
 *     }
 */
#   if !defined ASSERT_ONLY_PARAM
#     define ASSERT_ONLY_PARAM(PARAMETER)
#   endif

/*
 * Only define IGNORENOTIMPLEMENTED if building a "test" release build
 * that you are fully aware may contain unimplemented features.
 */
#   if !defined NOT_IMPLEMENTED
#     if defined IGNORENOTIMPLEMENTED
#       define NOT_IMPLEMENTED() CHECK(FALSE)
#     else
#       if defined MSWIN
/*
 * NOT_IMPLEMENTED is defined using a parameter, but should be called with none,
 * this will then throw a warning and not break the compile. Unix doesn't pick
 * up this warning, so break the compile under Unix
 */
#         define NOT_IMPLEMENTED(x)  TAssert("Not Implemented", __FILE__, __LINE__)
#       endif
#       if defined UNIXX
#         define NOT_IMPLEMENTED()  not implemented /* intentionally break the compile */
#       endif
#     endif
#   endif
#elif defined STD_ASSERTS
/* BEGINREMOVEFROMADDON */
#   define TASSERT(EXPR)         assert(EXPR)
/* ENDREMOVEFROMADDON */
#   if !defined INVARIANT
#     define INVARIANT(EXPR)       assert(EXPR)
#   endif
#   if !defined REQUIRE
#     define REQUIRE(EXPR)         assert(EXPR)
#   endif
#   if !defined ENSURE
#     define ENSURE(EXPR)          assert(EXPR)
#   endif
#   if !defined CHECK
#     define CHECK(EXPR)           assert(EXPR)
#   endif
#   ifdef VERIFY
#     undef VERIFY
#   endif
#   ifndef VERIFY
#     if defined NDEBUG
#       define VERIFY(EXPR) ((void)(EXPR))
#     else
#       define VERIFY(EXPR) assert(EXPR)
#     endif
#   endif /* VERIFY */
#   if !defined NOT_IMPLEMENTED
#     define NOT_IMPLEMENTED()     assert(!("Not Implemented"))
#   endif
    /*
     * See note above for this macro.
     */
#   if !defined ASSERT_ONLY_PARAM
#     define ASSERT_ONLY_PARAM(PARAMETER) PARAMETER
#   endif
#else
/* BEGINREMOVEFROMADDON */
#if defined (MSWIN)
#if defined CHECKED_BUILD
#include <string>
#include <vector>
#include <algorithm>

class AssertionLog
{
public:
    static void initializeAssertLog(const std::string &fileName);
    static bool isLoggingAssertions();
    static void addAssertion(const std::string &message);
private:
    static void writeOutAssertion(const std::string &message);
private:
    static bool                     logAssertions;
    static std::string              logFileName;
    static std::vector<std::string> assertList;
};

extern void TWinCheckedFailedLine(const char *Expr,
                                  const char *FileName,
                                  int LineNum);

#define TASSERT(EXPR)\
        do { if (!(EXPR)) { TWinCheckedFailedLine(#EXPR,__FILE__,__LINE__); } } while (0)
#else
#define TASSERT(EXPR) ASSERT(EXPR) /* MFC assert. 
Works in both release & debug builds */
#endif /* CHECKED_BUILD */
#else
#define TASSERT(EXPR) (void)((EXPR) || (TAssert(#EXPR, __FILE__, __LINE__), 0))
#endif

#   if defined NO_INVARIANTS
#     define INVARIANT(EXPR)
#   else
#     define INVARIANT(EXPR) TASSERT(EXPR)
#   endif

#   if defined NO_PRECONDITIONS
#     define REQUIRE(EXPR)
#   else
#     define REQUIRE(EXPR) TASSERT(EXPR)
#   endif

#   if defined NO_POSTCONDITIONS
#     define ENSURE(EXPR)
#   else
#     define ENSURE(EXPR) TASSERT(EXPR)
#   endif

#   if defined VERIFY
#     undef VERIFY
#   endif

#   if defined NO_CHECKS
#     define CHECK(EXPR)
#     define VERIFY(EXPR)  ((void)(EXPR))
#   else
#     define CHECK(EXPR)  TASSERT(EXPR)
#     if defined NDEBUG
#       define VERIFY(EXPR) ((void)(EXPR))
#     else
#       define VERIFY(EXPR) TASSERT(EXPR)
#     endif
#   endif

#   if defined NICE_NOT_IMPLEMENTED
#     define NOT_IMPLEMENTED() NiceNotImplemented()
#   else
#     define NOT_IMPLEMENTED() TASSERT(!("Not Implemented"))
#   endif
    /*
     * See note above for this macro.
     */
#   define ASSERT_ONLY_PARAM(PARAMETER) PARAMETER
/* ENDREMOVEFROMADDON */
#endif
/* BEGINREMOVEFROMADDON */
#if !defined STD_ASSERTS
extern void TecplotMopupOnAssert(void);
#endif /* !STD_ASSERTS */

#if defined NICE_NOT_IMPLEMENTED
extern void NiceNotImplemented(void);
#endif
/* ENDREMOVEFROMADDON */

/* convenience macros for implication, P -> Q, and equivalence, P <-> Q. */
#if !defined IMPLICATION
#  define IMPLICATION(P,Q) (!(P) || (Q))
#endif
#if !defined EQUIVALENCE
#  define EQUIVALENCE(P,Q) ((P) == (Q))
#endif

/* BEGINREMOVEFROMADDON */
#if defined RLM
#define VALID_RLM_HANDLE(h) ((h) != NULL)
#endif /* RLM */
/* ENDREMOVEFROMADDON */


#endif /* TASSERT_H */
