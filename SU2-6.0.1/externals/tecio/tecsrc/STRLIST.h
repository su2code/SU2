/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2010 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/
#if defined EXTERN
#  undef EXTERN
#endif
#if defined STRLISTMODULE
#  define EXTERN
#else
#  define EXTERN extern
#endif

#if !defined ARRLIST_h
#  error "Include ARRLIST.h before including STRLIST.h"
#endif

/*
 *
 * For building pltview.exe under Windows, we use
 * tecio.dll (which is linked to pltview).
 * Since pltview.exe uses a few of the
 * functions here, they need to be exported into
 * the tecio.dll, thus "TECXXX.h" is included for the
 * LIBFUNCTION & LIBCALL keywords. They are not
 * documented with the other TECXXX() functions,
 * however.
 *
 * If pltview requires other string functions
 * in the future, they can be added to the dll
 * by adding LIBFUNCTION & LIBCALL as in
 * StringListDealloc(), etc. below.
 *
 * When building the tecplot kernel, LIBFUNCTION
 * and LIBCALL are nop's.
 *
 */
#include "TECXXX.h"

EXTERN Boolean_t     StringListValid(StringList_pa StringList);
EXTERN void          StringListClear(StringList_pa StringList);
EXTERN void          StringListRemoveStrings(StringList_pa StringList,
                                             LgIndex_t     StringOffset,
                                             LgIndex_t     Count);
EXTERN void          StringListRemoveString(StringList_pa StringList,
                                            LgIndex_t     StringOffset);
LIBFUNCTION void LIBCALL StringListDealloc(StringList_pa* StringList);
EXTERN StringList_pa StringListAlloc(void);
EXTERN Boolean_t     StringListAppendString(StringList_pa StringList,
                                            char   const* String);
LIBFUNCTION LgIndex_t LIBCALL StringListCount(StringList_pa StringList);
LIBFUNCTION char* LIBCALL StringListGetString(StringList_pa StringList,
                                              LgIndex_t     StringOffset);

#if defined USE_MACROS_FOR_FUNCTIONS
#  define StringListGetStringRef StringListGetStringRef_MACRO
#else
#  define StringListGetStringRef StringListGetStringRef_FUNC
#endif

#if !defined USE_MACROS_FOR_FUNCTIONS
EXTERN char const* StringListGetStringRef_FUNC(StringList_pa StringList,
                                               LgIndex_t     StringOffset);
#endif
/**
 * To maintain the string list's integrity the result is cast to a
 * (char const*) to minimize the risk of users passing the result
 * to FREE_ARRAY.
 */
#define StringListGetStringRef_MACRO(StringList, StringOffset) \
   static_cast<char const*>(ArrayListGetCharPtr(reinterpret_cast<ArrayList_pa>(StringList), StringOffset))

EXTERN Boolean_t     StringListSetString(StringList_pa StringList,
                                         LgIndex_t     StringOffset,
                                         char const*   String);
EXTERN Boolean_t     StringListInsertString(StringList_pa StringList,
                                            LgIndex_t     StringOffset,
                                            char const*   String);
EXTERN StringList_pa StringListCopy(StringList_pa StringList);
EXTERN Boolean_t     StringListAppend(StringList_pa Target,
                                      StringList_pa Source);

EXTERN char*         StringListToNLString(StringList_pa StringList);
EXTERN StringList_pa StringListFromNLString(char const* String);
EXTERN char**        StringListToArray(StringList_pa StringList);
EXTERN StringList_pa StringListFromArray(char const** StringArray,
                                         LgIndex_t    Count);
EXTERN StringList_pa StringListFromCompound(char const* String);
EXTERN char*         StringListToCompound(StringList_pa StringList,
                                           char          GroupJoinCharacter,
                                           char const*   CharsToEscape);
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
