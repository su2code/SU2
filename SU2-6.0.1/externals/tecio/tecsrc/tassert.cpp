#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE

/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2010 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#define TASSERTMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined (MSWIN)
#endif
#endif

#include "STRUTIL.h"

using namespace tecplot::strutil;
using namespace std;

#define MAX_ERRMSG_LENGTH 2096

/* the mopup from assert and the writing out of crash.lay are */
/* used by TUASSERT and thus are needed even if NO_ASSERTS */
/* is set */
#if !defined NO_TU_ASSERTS || !defined NO_ASSERTS

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined MSWIN /* ...Unix/Linux calls this via signal handlers */
#endif
#endif /* TECPLOTKERNEL */

#endif /* Mopup function needed ... */



#if !defined STD_ASSERTS
/*
 * Define the final assertion notification function.
 */
#if defined UNIXX && !defined NO_ASSERTS
/*******************************************************************
 *                                                                 *
 *                          UNIX                                   *
 *                                                                 *
 *******************************************************************/


#  if defined NDEBUG
/*
 * if NDEBUG is defined __assert is NOT defined so we must supply
 * our own assertion notification function.....
 */
# define ASSERT assert
static void UnixAssert(const char *expression,
                       const char *file_name,
                       int        line)
{
    fprintf(stderr, "Assertion: %s\n"
            "Tecplot version: %s\n"
            "File Name: %s\n"
            "Line Number: %d\n",
            expression, TecVersionId, file_name, line);
    exit(ExitCode_AssertionFailure);
}
static TAssertFailureNotifyFunc assert_failure_notify = UnixAssert;
#  else
/*
 * NDEBUG is not defined so __assert is available....
 */
#    if defined LINUX
#      define LOWLEVELASSERTFUNCTION __linuxassertproxy
/*
 * In linux, __assert does not exist but rather
 * __assert_fail which has a differnt API.  Thus
 * a proxy is provided
 */
static void __linuxassertproxy(const char *__assertion,
                               const char *__file,
                               int         __line)
{
    __assert_fail(__assertion, __file, __line, __ASSERT_FUNCTION);
}
#    elif defined DARWIN
#      define LOWLEVELASSERTFUNCTION __darwinassertproxy
/*
 * In Darwin (Mac OS X), __assert is #defined to a call to __eprintf,
 * which also has a different API. Another proxy...
 */
static void __darwinassertproxy(const char *__assertion,
                                const char *__file,
                                int         __line)
{
    __eprintf("Assertion: %s\n"
              "Tecplot version: %s\n"
              "File Name: %s\n"
              "Line Number: %d\n",
              __assertion, TecVersionId, __file, (unsigned)__line);
}
#   else
#     define LOWLEVELASSERTFUNCTION __assert
#   endif

static TAssertFailureNotifyFunc assert_failure_notify = (TAssertFailureNotifyFunc) LOWLEVELASSERTFUNCTION;

#  endif
#endif /* UNIXX */

#if defined UNIXX && !defined NO_ASSERTS
/*
 * Replace the current assert failure notification function and
 * return the current one.
 *
 * Assumptions:
 *     new function points to valid function (not null) that
 *     conforms to the specified prototype
 *
 * Guarantees:
 *     result is a pointer to the previously installed
 *     function (not null)
 */
TAssertFailureNotifyFunc InstallTAssertFailureNotify(
    TAssertFailureNotifyFunc new_function) /* new notification function */
{
    TAssertFailureNotifyFunc result = 0; /* old function address */

    ASSERT(new_function != 0);

    result = assert_failure_notify;
    assert_failure_notify = new_function;

    ASSERT(result != 0);
    return result;
}






/*
 * Perform the installed assert failure notification action.
 */
void TAssert(const char *expression, /* text representation of the assertion */
             const char *file_name,  /* name of the file containing the assertion */
             int        line)        /* line number in the file of the assertion */
{
    static Boolean_t InTAssert = FALSE;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    char Message[MAX_ERRMSG_LENGTH+1];
    ASSERT(expression != 0 && strlen(expression) != 0);
    ASSERT(file_name != 0 && strlen(file_name) != 0);
    ASSERT(line >= 1);

    /* check for recursion */
    if (InTAssert)
    {
        fprintf(stderr, "Already in assert!\n");
        fprintf(stderr, "Assertion: %s\n"
                "Tecplot version: %s\n"
                "File Name: %s\n"
                "Line Number: %d\n",
                expression, TecVersionId, file_name, line);
        PrintCurBacktrace(stderr, 100);
        ASSERT(FALSE); /*... really exit */
    }

    InTAssert = TRUE;

    sprintf(Message, "Assertion: %s\n"
            "Tecplot version: %s\n"
            "File Name: %s\n"
            "Line Number: %d\n",
            expression, TecVersionId, file_name, line);

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# else
    fprintf(stderr, "%s", Message);
# endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

        (*assert_failure_notify)(expression, file_name, line);
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

    InTAssert = FALSE; /* just in case assert_failure_notify has an ignore */
}
#endif /* defined UNIXX && !defined NO_ASSERTS */
#endif /* STD_ASSERTS */


#if defined MSWIN && defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined CHECKED_BUILD
#endif
#if !defined ENGINE
#  if defined CHECKED_BUILD
#  endif // CHECKED_BUILD
#endif //!ENGINE
#endif /* MSWIN  */


#if defined NICE_NOT_IMPLEMENTED
static Boolean_t NotImplementedCalled = FALSE;
void NiceNotImplemented(void)
{
    if (!NotImplementedCalled)
    {
        Warning("Not Implemented!");
        NotImplementedCalled = TRUE;
    }
}
#endif
