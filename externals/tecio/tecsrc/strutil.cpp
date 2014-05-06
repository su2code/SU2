#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE

/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/

#define STRUTILMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#include "ARRLIST.h"
#include "STRLIST.h"
#include "CHARTYPE.h"
#include "STRUTIL.h"
#include "ALLOC.h"

#include "Q_MSG.h"

#include <algorithm>
#include <cctype> // ...needed to find std::tolower and std::toupper
#include <limits.h>
#include "TranslatedString.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

using std::string;
using tecplot::strutil::translate;
using tecplot::strutil::dontTranslate;
using tecplot::strutil::TranslatedString;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#ifdef MSWIN
# pragma warning (disable : 4786) /* STL warning about trucated identifiers */
#endif

/* END HEADER */

/**
 */
#define           INITIAL_FORMAT_BUFFER_SIZE 16384*3
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
static char      *FormatStringBuffer = NULL;
static int        FormatStringBufferSize = INITIAL_FORMAT_BUFFER_SIZE;

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
    #if defined MSWIN
    #else
    #endif /* !MSWIN */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/**
 * This should be one of the last functions called by Tecplot while mopping up.
 */
void FormatStringBufferCleanup(void)
{
    /*
     * NOTE: We use free instead of FREE_ARRAY for the scratch buffer because in
     *       debug mode FREE_ARRAY uses ErrMsg which uses vFormatString causing
     *       infinite recursion.
     */
    if (FormatStringBuffer != NULL)
        free(FormatStringBuffer);
    FormatStringBuffer = NULL;
}

/**
 */
char *vFormatString(const char *Format,
                    va_list     Arguments)
{
    char *Result = NULL;

    REQUIRE(VALID_REF(Format));

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

    /*
     * NOTE: We use malloc instead of ALLOC_ARRAY for the scratch buffer because
     *       in debug mode ALLOC_ARRAY uses ErrMsg which uses vFormatString
     *       causing infinite recursion.
     */
    if (FormatStringBuffer == NULL)
        FormatStringBuffer = (char *)malloc(FormatStringBufferSize);

    if (FormatStringBuffer != NULL)
    {
        Boolean_t TryAgain = FALSE;
        do
        {
            /*
             * Assign a value other than '\0' to the end of the buffer so that we
             * can determine if the buffer needs to be expanded. If after we call
             * vsnprintf the end of the buffer has a '\0' we need to expand it.
             */
            FormatStringBuffer[FormatStringBufferSize - 1] = (char)!'\0';

#         if defined MSWIN
            memset(FormatStringBuffer, 0, FormatStringBufferSize - 1);

            TryAgain =
                _vsnprintf(FormatStringBuffer,
                           FormatStringBufferSize,
                           Format,
                           Arguments) == -1;
#         elif defined IRIX62
            vsprintf(FormatStringBuffer,
                     Format,
                     Arguments);
            CHECK(strlen(FormatStringBuffer) < FormatStringBufferSize);
#         else
            vsnprintf(FormatStringBuffer,
                      FormatStringBufferSize,
                      Format,
                      Arguments);
#         endif

#ifndef MSWIN
            TryAgain = (FormatStringBuffer[FormatStringBufferSize - 1] == '\0');
#endif
            if (TryAgain)
            {
                /*
                 * Reallocate the buffer and try again.
                 *
                 * NOTE: We use malloc/free instead of ALLOC/FREE_ARRAY for the
                 *       scratch buffer because in debug mode ALLOC/FREE_ARRAY
                 *       uses ErrMsg which uses vFormatString causing infinite
                 *       recursion.
                 */
                free(FormatStringBuffer);
                FormatStringBufferSize += MAX(1, FormatStringBufferSize / 2);
                FormatStringBuffer = (char *)malloc(FormatStringBufferSize);
                TryAgain = (FormatStringBuffer != NULL);
                if (!TryAgain)
                    FormatStringBufferSize = INITIAL_FORMAT_BUFFER_SIZE;
            }
        }
        while (TryAgain);

        if (FormatStringBuffer != NULL)
            Result = DupString(dontTranslate(FormatStringBuffer));
    }

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

    ENSURE(VALID_REF(Result) || Result == NULL);
    return Result;
}

/**
 */
char *FormatString(TranslatedString Format,
                   ...) /* 0 or more variable arguments */
{
    REQUIRE(!Format.isNull());

    va_list Arguments;
    va_start(Arguments, Format);
    char *Result = vFormatString(Format.c_str(), Arguments);
    va_end(Arguments);

    ENSURE(VALID_REF(Result) || Result == NULL);
    return Result;
}

/**
 */
int FormatString(string&           Buffer,
                 TranslatedString  Format
                 ...) /* 0 or more variable arguments */
{
    REQUIRE(!Format.isNull());

    va_list Arguments;
    va_start(Arguments, Format);
    char *FormattedString = vFormatString(Format.c_str(), Arguments);
    va_end(Arguments);

    int Result;
    if (FormattedString != NULL)
    {
        Buffer.assign(FormattedString);
        Result = (int)Buffer.size();
        FREE_ARRAY(FormattedString, "FormattedString");
    }
    else
        Result = -1;


    ENSURE(Result == -1 || Result >= 0);
    return Result;
}

/**
 * Returns a duplicate of the string or NULL if sufficient memory is not
 * available.
 *
 * NOTE: This function was created because ResetString(...) does not
 *       duplicate zero length strings but returns NULL instead.
 */
char *DupString(TranslatedString String)
{
    REQUIRE(VALID_TRANSLATED_STRING(String));

    char *Result = ALLOC_ARRAY(strlen(String.c_str()) + 1, char, "duplicate string");
    if (Result != NULL)
        strcpy(Result, String.c_str());

    ENSURE(Result == NULL || (VALID_REF(Result) && strcmp(Result, String.c_str()) == 0));
    return Result;
}


/*
 * Copy up to 'Count' characters from the 'Source' string beginning at
 * position 'Index' to the 'Target' string. The actual number of characters
 * copied may be less than 'Count' if a '\0' was encountered in the
 * 'Source' string before 'Count' characters were copied.
 *
 * NOTE: The 'Target' and 'Source' strings may overlap.
 */
void CopySubString(char       *Target,
                   const char *Source,
                   int         Index,
                   int         Count)
{
    LgIndex_t Length = 0;

    REQUIRE(VALID_REF(Target));
    REQUIRE("Target string is sized to accommodate a string who's length "
            "is at least MIN(strlen(&Source[Index]), Count) characters.");
    REQUIRE(VALID_REF(Source));
    REQUIRE(0 <= Index && Index <= (LgIndex_t)strlen(Source));
    REQUIRE(Count >= 0);

    Length = MIN((LgIndex_t)strlen(&Source[Index]), Count);
    memmove(Target, &Source[Index], Length);
    Target[Length] = '\0';

    ENSURE(VALID_REF(Target) && (LgIndex_t)strlen(Target) == Length);
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

/*
 * Remove any leading white space from the string and return
 * a reference to it. NOTE: The input string is modified.
 */
char *StringFlushLeft(char *String)
{
    char *Result = String;
    char *Start = String;

    REQUIRE(VALID_REF(String));

    /* move the substring beginning at the first non-whitespace */
    /* character to the head of the string                      */
    while (tecplot::isspace(*Start))
        Start++;
    if (Start != String)
        memmove(String, Start, strlen(Start) + 1);

    ENSURE(VALID_REF(Result) && Result == String);
    return Result;
}


/*
 * Remove any trailing white space from the string and return
 * a reference to it. NOTE: The input string is modified.
 */
static char *StringFlushRight(char *String)
{
    char *Result = String;
    char *End = NULL;

    REQUIRE(VALID_REF(String));

    for (End = EndOfString(String); End != String && tecplot::isspace(End[-1]); End--)
        End[-1] = '\0';

    ENSURE(VALID_REF(Result) && Result == String);
    return Result;
}


/*
 * Remove any leading and trailing white space from the string
 * and return a reference to it.  The return value is not
 * absolutely necessary since the input string is modified
 * but it is convenient sometimes.
 * NOTE: The input string is modified but no memory is
 * allocated nor deallocated.
 */
char *TrimLeadAndTrailSpaces(char *String)
{
    REQUIRE((String == NULL) || VALID_REF(String));
    if (String)
        return (StringFlushLeft(StringFlushRight(String)));
    else
        return String;
}


/*
 * If the specified string is longer than the maximum specified length
 * truncate it and return a reference to it.
 *
 * String
 *   String to truncate if necessary.
 * MaxLength
 *   Length at which to truncate the specified string if exceeded.
 *
 * Return
 *   Reference to the input string.
 */

// Okay for UTF-8
char *StringTruncate(char      *String,
                     LgIndex_t MaxLength)
{
    REQUIRE(VALID_REF(String));
    REQUIRE(MaxLength >= 0);

    if ((LgIndex_t)strlen(String) > MaxLength)
        String[MaxLength] = '\0';/* UTF8_SetAt(String,'\0',MaxLength); */

    ENSURE(VALID_REF(String));
    ENSURE((LgIndex_t)strlen(String) <= MaxLength);
    return String;
}


/*
 * Trim and truncate the specified string such that its trimmed length
 * does not exceed the specified length and return a reference to it.
 *
 * String
 *   String to trim and truncate if necessary.
 * MaxLength
 *   Length at which to truncate the trimmed string if exceeded.
 *
 * Return
 *   Reference to the input string.
 */
char *StringTrimAndTruncate(char      *String,
                            LgIndex_t MaxLength)
{
    REQUIRE(VALID_REF(String));
    REQUIRE(MaxLength >= 0);

    /*
     * Note that we are careful to truncate the string after trimming
     * whitespace from the left side but then trim whitespace from the end
     * after truncating to make sure we don't return a string that has
     * whitespace simply because it truncated on a word break.
     */
    StringFlushLeft(String);
    StringTruncate(String,MaxLength);
    StringFlushRight(String);

    ENSURE(VALID_REF(String));
    ENSURE((LgIndex_t)strlen(String) <= MaxLength);
    ENSURE(IMPLICATION(strlen(String) != 0,
                       (!tecplot::isspace(String[0]) &&
                        !tecplot::isspace(String[strlen(String)-1]))));
    return String;
}

/**
 */

#ifndef MSWIN
StringList_pa LineBreakString(const char *String,
                              UInt32_t    WrapMargin)
{
    REQUIRE(VALID_REF(String));

    StringList_pa Result = StringListAlloc();
    if (Result != NULL)
    {
        Boolean_t IsOk = TRUE;
        if (strlen(String) > WrapMargin)
        {
            char *StringCopy = DupString(dontTranslate(String));
            IsOk = (StringCopy != NULL);
            if (IsOk)
            {
                char *CPtr = StringCopy;
                char *SubString = StringCopy;
                UInt32_t SubStringLen = 0;
                while (*CPtr != '\0' && IsOk)
                {
                    while (*CPtr != '\0' && SubStringLen < WrapMargin)
                    {
                        /* look for a hard break */
                        if (*CPtr == '\n')
                        {
                            *CPtr = '\0'; /* replace the newline */
                            CPtr++;
                            break;
                        }

                        CPtr++;
                        SubStringLen++;
                    }

                    /*
                     * If we didn't find a hard break or the end of the string
                     * then we need to back up and find the closest space.
                     */
                    if (*CPtr != '\0' && SubStringLen == WrapMargin)
                    {
                        /* find the closes space from the right */
                        if (*CPtr != ' ')
                        {
                            while (CPtr != SubString && *CPtr != ' ')
                                CPtr--;
                            if (*CPtr != ' ')
                            {
                                /*
                                 * Bummer, this line will exceed the wrap margin.
                                 * Search forward for the next space or newline.
                                 */
                                while (*CPtr != '\0' && *CPtr != ' ' && *CPtr != '\n')
                                    CPtr++;
                                while (*CPtr != '\0' && *CPtr == ' ')
                                    CPtr++; /* skip over the white space */
                            }
                        }

                        if (*CPtr != '\0')
                        {
                            *CPtr = '\0';
                            CPtr++;
                        }
                        StringFlushRight(SubString);
                    }

                    IsOk = StringListAppendString(Result, SubString);
                    SubString = CPtr;
                    SubStringLen = 0;
                }

                FREE_ARRAY(StringCopy, "StringCopy");
            }
        }
        else
            IsOk = StringListAppendString(Result, String);

        if (!IsOk)
            StringListDealloc(&Result);
    }

    ENSURE(Result == NULL || VALID_REF(Result));
    return Result;
}
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */


/**
 * Lexicographically compares, at most, the first 'Len' characters of
 * s1 and s2.
 *
 * param s1
 *     First string or NULL.
 * param s2
 *     Second string or NULL.
 * param Len
 *     Maximum number of characters to compare.
 * return
 *     Integer value greater than, equal to, or less than zero according
 *     to whether the first 'Len' characters of 's1' are greater than,
 *     equal to, or less than 's2'.
 */

// Okay for UTF-8
int ustrncmp(const  char *s1,
             const  char *s2,
             size_t Len)
{
    REQUIRE((s1 == NULL) || VALID_REF(s1));
    REQUIRE((s2 == NULL) || VALID_REF(s2));

    char *t1;
    char *t2;
    char ct1;
    char ct2;
    size_t I = 0;
    if ((s1 == NULL) && (s2 == NULL))
        return 0;
    if (s1 == NULL)
        return -1;
    else if (s2 == NULL)
        return 1;

    t1 = (char*)s1;
    t2 = (char*)s2;

    while (*t1 && *t2 && (I < Len))
    {
        ct1 = CAPITAL(*t1);
        ct2 = CAPITAL(*t2);
        if (ct1 != ct2)
            return (ct1 - ct2);
        t1++;
        t2++;
        I++;
    }
    if ((I == Len) ||
        ((*t1 == '\0') && (*t2 == '\0')))
        return 0;
    else
        return CAPITAL(*t1) - CAPITAL(*t2);


}


/**
 * Lexicographically compares the characters of s1 and s2.
 *
 * param s1
 *     First string or NULL.
 * param s2
 *     Second string or NULL.
 * return
 *     Integer value greater than, equal to, or less than zero according to
 *     whether the characters of 's1' are greater than, equal to, or less
 *     than 's2'.
 */
int ustrcmp(const char *s1,
            const char *s2)
{
    REQUIRE((s1 == NULL) || VALID_REF(s1));
    REQUIRE((s2 == NULL) || VALID_REF(s2));

    return (ustrncmp(s1, s2, INT_MAX));
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined NO_ASSERTS && defined DEBUG_ALLOC
#endif
#endif /* TECPLOTKERNEL */


/*
 * The problem with passing file names for release builds is that
 * the full path name is used (i.e., c:\user\craig\v7.5\tecplot\alloc.c)
 */

// Okay for UTF-8

#if !defined NO_ASSERTS
Boolean_t InternalResetString(char       **SBase,
                              const char  *NewString,
                              Boolean_t    IssueErrMsg,
                              const char  *FileName,
                              int          LineNumber)
#else
Boolean_t InternalResetString(char       **SBase,
                              const char  *NewString,
                              Boolean_t    IssueErrMsg)
#endif
{
    REQUIRE(VALID_REF(SBase));
    REQUIRE(*SBase == NULL || VALID_REF(*SBase));
    REQUIRE(NewString == NULL || VALID_REF(NewString));
    REQUIRE(IMPLICATION(VALID_REF(*SBase), *SBase != NewString)); /* Prevent calling with same string. */
    REQUIRE(VALID_BOOLEAN(IssueErrMsg));
    REQUIRE(VALID_NON_ZERO_LEN_STR(FileName));
    REQUIRE(LineNumber >= 1);

    if (*SBase)
    {
#if !defined NO_ASSERTS && defined DEBUG_ALLOC
        char S[80+1];
        MakeDebugRecord(FileName, LineNumber, "releasing", *SBase, S, 80);
        FREE_ARRAY(*SBase, S);
#else
        FREE_ARRAY(*SBase, "");
#endif
    }
    if (NewString == NULL)
    {
        *SBase = NULL;
        return (TRUE);
    }
    else
    {
#if !defined NO_ASSERTS && defined DEBUG_ALLOC
        char S[80+1];
        MakeDebugRecord(FileName, LineNumber, "duplicating", NewString, S, 80);
        *SBase = ALLOC_ARRAY(strlen(NewString) + 1, char, S);
#else
#  if defined MSWIN && defined _DEBUG && !defined(MAKEARCHIVE) && !defined(NO_ASSERTS)
        /* Allow the MFC memory leak detection to report the leak at the
         * calling programs location, and not here (which is fairly useless).
         * But first, we have to turn off the preprocessor definition of new
         * and get to the original. */
#     undef new
        *SBase = new(FileName, LineNumber) char[strlen(NewString)+1];
#     define new DEBUG_NEW
#  else
        *SBase = ALLOC_ARRAY(strlen(NewString) + 1, char, "");
#  endif
#endif
        if (*SBase)
        {
            strcpy(*SBase, NewString);
            return (TRUE);
        }
        else
        {
            if (IssueErrMsg)
                ErrMsg(translate("Out of memory"));
            return (FALSE);
        }
    }
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif



/*
 * Unfortunately, this routine uses the interface of DeleteStringToAdd
 * forcing StringToAdd to be non-const.  Another copy of this routine
 * for const char *'s is below.  Eventually we should get rid of
 * the two routines and DeleteStringToAdd, always using the const version
 * and deleting the string if necessary in the code that calls TackOnString.
 */
Boolean_t TackOnString(char       **SBase,
                       const char  *StringToAdd,
                       Boolean_t    DeleteStringToAdd,
                       Boolean_t    ConvertNewlinesToAscii)
{
    size_t      CurLen;
    size_t      NewLen;
    int         NumNewlines = 0;
    char       *NewString;
    const char *CPtr = StringToAdd;
    Boolean_t   IsOk = TRUE;

    REQUIRE(VALID_REF(SBase));
    REQUIRE((StringToAdd == NULL) || VALID_REF(StringToAdd));
    REQUIRE(VALID_BOOLEAN(DeleteStringToAdd));
    REQUIRE(VALID_BOOLEAN(ConvertNewlinesToAscii));

    if ((StringToAdd  == NULL) ||
        (*StringToAdd == '\0'))
    {
        if (StringToAdd            &&
            (*StringToAdd == '\0') &&
            DeleteStringToAdd)
        {
            char *TMP = (char *)StringToAdd;
            FREE_ARRAY(TMP, "empty string to add");
        }
    }
    else
    {

        if (*SBase == NULL)
            CurLen = 0;
        else
            CurLen = strlen(*SBase);

        while (*CPtr)
            if (*CPtr++ == '\n')
                NumNewlines++;

        NewLen = CurLen + strlen(StringToAdd) + 1 + NumNewlines;

        NewString = ALLOC_ARRAY(NewLen, char, StringToAdd);

        if (NewString == NULL)
        {
            if (DeleteStringToAdd)
            {
                char *TMP = (char *)StringToAdd;
                FREE_ARRAY(TMP, StringToAdd);
            }
            IsOk = FALSE;
        }
        else
        {
            if (*SBase)
            {
                strcpy(NewString, *SBase);
                FREE_ARRAY(*SBase, (CurLen > 0 ? *SBase : "previous text"));
            }
            else
                *NewString = '\0';

            {
                char *NPtr = EndOfString(NewString);
                const char *APtr = StringToAdd;
                while (*APtr)
                {
                    if ((*APtr == '\n') && ConvertNewlinesToAscii)
                    {
                        *NPtr++ = '\\';
                        *NPtr++ = 'n';
                    }
                    else
                        *NPtr++ = *APtr; //UTF8_AssignAndIncrement(&NPtr,(char**)&APtr,TRUE,FALSE); //*NPtr++ = *APtr;
                    APtr++;
                }
                *NPtr = '\0';
            }

            if (DeleteStringToAdd)
            {
                char *TMP = (char *)StringToAdd;
                FREE_ARRAY(TMP, StringToAdd);
            }

            *SBase = NewString;
        }
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return (IsOk);
}


/*
 * See TackOnString for discussion.
 */

// Okay for UTF-8
Boolean_t TackOnConstString(char      **SBase,
                            const char *StringToAdd,
                            Boolean_t   ConvertNewlinesToAscii)
{
    size_t      CurLen;
    size_t      NewLen;
    int         NumNewlines = 0;
    char       *NewString;
    const char *CPtr = StringToAdd;
    Boolean_t   IsOk = TRUE;

    REQUIRE(VALID_REF(SBase));
    REQUIRE((StringToAdd == NULL) || VALID_REF(StringToAdd));
    REQUIRE(VALID_BOOLEAN(ConvertNewlinesToAscii));

    if ((StringToAdd  != NULL) &&
        (*StringToAdd != '\0'))
    {
        if (*SBase == NULL)
            CurLen = 0;
        else
            CurLen = strlen(*SBase);

        while (*CPtr)
            if (*CPtr++ == '\n')
                NumNewlines++;

        NewLen = CurLen + strlen(StringToAdd) + 1 + NumNewlines;

        NewString = ALLOC_ARRAY(NewLen, char, StringToAdd);

        if (NewString == NULL)
        {
            IsOk = FALSE;
        }
        else
        {
            if (*SBase)
            {
                strcpy(NewString, *SBase);
                FREE_ARRAY(*SBase, (CurLen > 0 ? *SBase : "previous text"));
            }
            else
                *NewString = '\0';

            {
                char *NPtr = EndOfString(NewString);
                const char *APtr = StringToAdd;
                while (*APtr)
                {
                    if ((*APtr == '\n') && ConvertNewlinesToAscii)
                    {
                        *NPtr++ = '\\';
                        *NPtr++ = 'n';
                    }
                    else
                        *NPtr++ = *APtr; // UTF8_AssignAndIncrement(&NPtr,(char**)&APtr,TRUE,FALSE);
                    APtr++;
                }
                *NPtr = '\0';
            }
            *SBase = NewString;
        }
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return (IsOk);
}


// Okay for UTF-8
Boolean_t TackOnChar(char  **SBase,
                     char    CharToAdd)
{
    REQUIRE(VALID_REF(SBase));

    char S[2];
    S[0] = CharToAdd;
    S[1] = '\0';

    return (TackOnString(SBase, S, FALSE, FALSE));
}



/**
 * Converts all one character new line characters in the allocated string
 * to a two character "\n" sequence.
 *
 * param String
 *     String to scan and convert if necessary. Note that the string will
 *     be reallocated if any new line characters are discovered.
 *
 * return
 *     TRUE if the request was successfull, FALSE otherwise.
 */

// Okay for UTF-8
Boolean_t ReplaceNewlineWithBackslashN(char **String)
{
    size_t    I;
    LgIndex_t NewlineCount;
    size_t    Length;
    char     *Replacement;

    REQUIRE(VALID_REF(String));
    REQUIRE(VALID_REF(*String));

    /* count how many new line character are present */
    NewlineCount = 0;
    Length = strlen(*String);
    for (I = 0; I < Length; I++)
        if ((*String)[I] == '\n')
            NewlineCount++;

    if (NewlineCount != 0)
    {
        /* allocate a new string and convert */
        Replacement = ALLOC_ARRAY(Length + NewlineCount + 1, char,
                                  "replacement string");
        if (Replacement != NULL)
        {
            size_t J;
            for (I = J = 0; I < Length + 1; I++, J++)
            {
                if ((*String)[I] == '\n')
                {
                    Replacement[J] = '\\';
                    J++;
                    Replacement[J] = 'n';
                }
                else
                {
                    Replacement[J] = (*String)[I];
                }
            }

            /* sanity check */
            CHECK(I == Length + 1);
            CHECK(J == Length + NewlineCount + 1);
        }

        /* release the old string and record the new one */
        FREE_ARRAY(*String, "original string");
        *String = Replacement;
    }

    ENSURE(*String == NULL || VALID_REF(*String));
    return (*String != NULL);
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined TECPLOTKERNEL
#if !defined NO_ASSERTS
#endif /* !NO_ASSERTS */
#endif /*TECPLOTKERNEL*/
#endif
