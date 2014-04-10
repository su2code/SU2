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

#define STRLISTMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "STRUTIL.h"
#include "ALLOC.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#include "ARRLIST.h"
#include "STRLIST.h"

/* END HEADER */

using tecplot::strutil::dontTranslate;

/*
 * This set of functions provide a wrapper around the array list utilities
 * thereby making it aware of item allocation and deallocation. All strings
 * given to the string list and returned to the client are copies. Therefore
 * it is the client's responsibility to deallocate string results when no
 * longer needed.
 */


/*
 * Destructor for cleaning up string allocations.
 *
 * param ItemRef
 *     Reference to the string item to destroy.
 * param ClientData
 *     Any client data needed for destroying the string.
 *
 * return
 *     TRUE is a requirement
 */
static Boolean_t StringListItemDestructor(void*      ItemRef,
                                          ArbParam_t ClientData)
{
    REQUIRE(VALID_REF(ItemRef));
    REQUIRE(VALID_REF(*static_cast<char**>(ItemRef)) || *static_cast<char**>(ItemRef) == NULL);
    UNUSED(ClientData);

    char** StringRef = static_cast<char**>(ItemRef);
    if (*StringRef != NULL)
    {
        FREE_ARRAY(*StringRef, "string");
        *StringRef = NULL;
    }

    ENSURE(*StringRef == NULL);
    return TRUE;
}

/*
 * String item duplicator.
 *
 * param TargetItemRef
 *     Reference to the string list item to receive the duplicate.
 * param SourceItemRef
 *     Reference to the string list item to duplicate.
 * param ClientData
 *     Any client data required for duplication.
 *
 * return
 *     TRUE if the duplication was a success
 *     FALSE otherwise. If the duplication failed it
 *     is the client's responsibility to cleanup any
 *     partial duplication
 */
static Boolean_t StringListItemDuplicator(void*      TargetItemRef,
                                          void*      SourceItemRef,
                                          ArbParam_t ClientData)
{
    REQUIRE(VALID_REF(TargetItemRef));
    REQUIRE(VALID_REF(SourceItemRef));
    REQUIRE(VALID_REF(*static_cast<char**>(SourceItemRef)) || *static_cast<char**>(SourceItemRef) == NULL);
    UNUSED(ClientData);

    Boolean_t IsOk = TRUE;
    char** TargetStringRef = static_cast<char**>(TargetItemRef);
    char** SourceStringRef = static_cast<char**>(SourceItemRef);

    if (*SourceStringRef != NULL)
        IsOk = ((*TargetStringRef = DupString(dontTranslate(*SourceStringRef))) != NULL);
    else
        *TargetStringRef = NULL;

    ENSURE(VALID_REF(*TargetStringRef) || *TargetStringRef == NULL);
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/*
 * Determine if the string list handle and its members are sane.
 */
Boolean_t StringListValid(StringList_pa stringList)
{
    Boolean_t isValid = ArrayListIsValid(reinterpret_cast<ArrayList_pa>(stringList));
    if (isValid)
    {
        LgIndex_t stringCount = ArrayListGetCount(reinterpret_cast<ArrayList_pa>(stringList));

        #if defined PERFORM_EXPENSIVE_STRLIST_TESTS
        {
            for (LgIndex_t index = 0; index < stringCount; index++)
            {
                char* string = ArrayListGetCharPtr(reinterpret_cast<ArrayList_pa>(stringList), index);
                if (string != NULL && !VALID_REF(string))
                {
                    isValid = FALSE;
                    break;
                }
            }
        }
        #else
        {
            /* Check first and last only */
            if (stringCount > 0)
            {
                char* string = ArrayListGetCharPtr(reinterpret_cast<ArrayList_pa>(stringList), 0);
                if (string != NULL && !VALID_REF(string))
                {
                    isValid = FALSE;
                }
            }
            if (isValid && stringCount > 1)
            {
                char* string = ArrayListGetCharPtr(reinterpret_cast<ArrayList_pa>(stringList), stringCount - 1);
                if (string != NULL && !VALID_REF(string))
                {
                    isValid = FALSE;
                }
            }
        }
        #endif /* PERFORM_SKIP_EXPENSIVE_STRLIST_TESTS */
    }

    ENSURE(VALID_BOOLEAN(isValid));
    return isValid;
}


/*
 * Remove all members of the string list.
 */
void StringListClear(StringList_pa StringList)
{
    REQUIRE(StringListValid(StringList));

    ArrayListDeleteAllItems(reinterpret_cast<ArrayList_pa>(StringList), StringListItemDestructor, 0);

    ENSURE(StringListValid(StringList) && StringListCount(StringList) == 0);
}


/*
 * Remove 'Count' strings from the list beginning at the specified offset.
 * The members following the items removed are shifted down accordingly to
 * fill the vacated space.
 */
void StringListRemoveStrings(StringList_pa StringList,
                             LgIndex_t     StringOffset,
                             LgIndex_t     Count)
{
    REQUIRE(StringListValid(StringList));
    REQUIRE(0 <= StringOffset && StringOffset <= StringListCount(StringList) - 1);
    REQUIRE(1 <= Count && StringOffset + Count <= StringListCount(StringList));

    ArrayListDeleteItems(reinterpret_cast<ArrayList_pa>(StringList), StringOffset, Count,
                         StringListItemDestructor, 0);

    ENSURE(StringListValid(StringList));
}


/*
 * Remove the string from the list at the specified offset. The members
 * following the item removed are shifted down accordingly to fill the
 * vacated space.
 */
void StringListRemoveString(StringList_pa StringList,
                            LgIndex_t     StringOffset)
{
    REQUIRE(StringListValid(StringList));
    REQUIRE(0 <= StringOffset && StringOffset <= StringListCount(StringList) - 1);

    ArrayListDeleteItems(reinterpret_cast<ArrayList_pa>(StringList), StringOffset, 1,
                         StringListItemDestructor, 0);

    ENSURE(StringListValid(StringList));
}


/*
 * Deallocate the string list members and handle and set the handle to NULL.
 */
void LIBCALL StringListDealloc(StringList_pa* StringList)
{
    REQUIRE(VALID_REF(StringList));
    REQUIRE(*StringList == NULL || StringListValid(*StringList));

    if (*StringList != NULL)
        ArrayListDealloc(reinterpret_cast<ArrayList_pa*>(StringList), StringListItemDestructor, 0);

    ENSURE(*StringList == NULL);
}


/*
 * Allocate a string list handle. A handle of NULL is
 * returned if sufficient memory is not available.
 */
StringList_pa StringListAlloc(void)
{
    StringList_pa Result = reinterpret_cast<StringList_pa>(ArrayListAlloc(0, ArrayListType_CharPtr, NULL, 0));

    ENSURE(Result == NULL || StringListValid(Result));
    return Result;
}


/*
 * Append a copy of the string to the string list. The string list will be
 * expanded to accommodate the additional item. A return value of TRUE
 * indicates the operation was successful otherwise FALSE is returned
 * indicating that sufficient memory was not available for the additional
 * item.
 */
Boolean_t StringListAppendString(StringList_pa StringList,
                                 char const*   String)
{
    REQUIRE(StringListValid(StringList));
    REQUIRE(String == NULL || VALID_REF(String));

    Boolean_t IsOk = StringListSetString(StringList, StringListCount(StringList), String);

    ENSURE(StringListValid(StringList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/*
 * Return the number of strings currently in the string list.
 */
LgIndex_t LIBCALL StringListCount(StringList_pa StringList)
{
    REQUIRE(StringListValid(StringList));

    LgIndex_t Result = ArrayListGetCount(reinterpret_cast<ArrayList_pa>(StringList));

    ENSURE(Result >= 0);
    return Result;
}


/*
 * Return a copy of the string at the specified offset in the string list.
 */
char* LIBCALL StringListGetString(StringList_pa StringList,
                                  LgIndex_t     StringOffset)
{
    REQUIRE(StringListValid(StringList));
    REQUIRE(0 <= StringOffset && StringOffset <= StringListCount(StringList) - 1);

    char* Result;
    char const* StringRef = StringListGetStringRef(StringList, StringOffset);
    if (StringRef == NULL)
        Result = NULL;
    else
        Result = DupString(dontTranslate(StringRef));

    ENSURE(Result == NULL || VALID_REF(Result));
    return Result;
}


#if !defined USE_MACROS_FOR_FUNCTIONS
/*
 * Returns actual string at the specified offset in the string list.  Do not
 * attempt to free this string.  Changing this string should be done with
 * utmost caution.
 */
char const* StringListGetStringRef_FUNC(StringList_pa StringList,
                                        LgIndex_t     StringOffset)
{
    REQUIRE(StringListValid(StringList));
    REQUIRE(0 <= StringOffset && StringOffset <= StringListCount(StringList) - 1);

    char const* Result = StringListGetStringRef_MACRO(StringList, StringOffset);

    ENSURE(Result == NULL || VALID_REF(Result));
    return Result;
}
#endif


/*
 * Place a copy of the specified string at the specified offset. If the offset
 * is beyond the end of the string list it is sized accordingly and the
 * intervening string references between the last item of the original
 * state and the last item of the new state are assigned NULL. If a string
 * already exists at the specified location its resources are released.
 * A return value of TRUE indicates the operation was successful otherwise
 * FALSE is returned indicating that sufficient memory was not available
 * for the additional item at the specified offset.
 */
Boolean_t StringListSetString(StringList_pa StringList,
                              LgIndex_t     StringOffset,
                              char const*   String)
{
    REQUIRE(StringListValid(StringList));
    REQUIRE(StringOffset >= 0);
    REQUIRE(String == NULL || VALID_REF(String));

    Boolean_t       IsOk;
    ArrayListItem_u ItemCopy;

    if (String != NULL)
    {
        ItemCopy.CharPtr = DupString(dontTranslate(String));
        IsOk = (ItemCopy.CharPtr != NULL);
    }
    else
    {
        ItemCopy.CharPtr = NULL;
        IsOk = TRUE;
    }

    if (IsOk)
        IsOk = ArrayListSetItem(reinterpret_cast<ArrayList_pa>(StringList), StringOffset, ItemCopy,
                                StringListItemDestructor, 0);

    ENSURE(StringListValid(StringList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/*
 * Insert a copy of the string into the string list at the specified offset.
 * The string list will be expanded to accommodate the additional item.
 * A return value of TRUE indicates the operation was successful otherwise
 * FALSE is returned indicating that sufficient memory was not available
 * for the additional item.
 */
Boolean_t StringListInsertString(StringList_pa StringList,
                                 LgIndex_t     StringOffset,
                                 char const*   String)
{
    REQUIRE(StringListValid(StringList));
    REQUIRE(StringOffset >= 0);
    REQUIRE(String == NULL || VALID_REF(String));

    Boolean_t       IsOk;
    ArrayListItem_u ItemCopy;

    if (String != NULL)
    {
        ItemCopy.CharPtr = DupString(dontTranslate(String));
        IsOk = (ItemCopy.CharPtr != NULL);
    }
    else
    {
        ItemCopy.CharPtr = NULL;
        IsOk = TRUE;
    }

    if (IsOk)
        IsOk = ArrayListInsertItem(reinterpret_cast<ArrayList_pa>(StringList),
                                   StringOffset, ItemCopy);

    ENSURE(StringListValid(StringList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/*
 * Return a handle to a duplicate of the specified string list and its contents.
 * A handle of NULL is returned if sufficient memory is not available.
 */
StringList_pa StringListCopy(StringList_pa StringList)
{
    REQUIRE(StringListValid(StringList));

    StringList_pa Result =
        reinterpret_cast<StringList_pa>(ArrayListCopy(reinterpret_cast<ArrayList_pa>(StringList),
                                                      StringListItemDuplicator, 0));

    ENSURE(Result == NULL ||
           (StringListValid(Result) &&
            StringListCount(Result) == StringListCount(StringList)));
    return Result;
}



/*
 * Append a copy of the contents of the source list to the target list.
 * A return value of TRUE indicates the operation was successful otherwise
 * FALSE is returned indicating that sufficient memory was not available
 * for the request.
 */
Boolean_t StringListAppend(StringList_pa Target,
                           StringList_pa Source)
{
    REQUIRE(StringListValid(Target));
    REQUIRE(StringListValid(Source));

    StringList_pa SourceCopy = StringListCopy(Source);
    Boolean_t IsOk = (SourceCopy != NULL);
    if (IsOk)
    {
        ArrayListAppend(reinterpret_cast<ArrayList_pa>(Target), reinterpret_cast<ArrayList_pa>(SourceCopy));
        /* deallocate the list but not the string items since Target now owns them */
        ArrayListDealloc(static_cast<ArrayList_pa*>(static_cast<void*>(&SourceCopy)), NULL, 0);
    }

    ENSURE(StringListValid(Target));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}



/*
 * Return a new line, '\n', separated string representation of the string list.
 * Caller is responsible for de-allocating the result.
 */
char* StringListToNLString(StringList_pa StringList)
{
    REQUIRE(StringListValid(StringList));

    size_t Length = 0;

    /* determine the resulting new line, '\n', separated string length */
    LgIndex_t Count = StringListCount(StringList);
    if (Count >= 1)
    {
        LgIndex_t Index;
        for (Index = 0, Length = strlen("\n") * (Count - 1);
             Index < Count;
             Index++)
        {
            char* String = ArrayListGetCharPtr(reinterpret_cast<ArrayList_pa>(StringList), Index);
            if (String != NULL)
                Length += strlen(String);
        }
    }

    /* create a new line, '\n', separated string */
    char* Result = ALLOC_ARRAY(Length + 1, char, "new line separated string");
    if (Result != NULL)
    {
        LgIndex_t Index;
        for (Index = 0, strcpy(Result, "");
             Index < Count;
             Index++)
        {
            char* String = ArrayListGetCharPtr(reinterpret_cast<ArrayList_pa>(StringList), Index);

            if (Index != 0)
                strcat(Result, "\n");

            if (String != NULL)
                strcat(Result, String);
        }
    }

    ENSURE(Result == NULL || VALID_REF(Result));
    return Result;
}


/*
 * Create a string list from the new line, '\n', separated string. The string
 * is copied and therefore owned and managed by the caller.
 */
StringList_pa StringListFromNLString(char const* String)
{
    REQUIRE(VALID_REF(String));

    /* create the string list and scan the entire string */
    StringList_pa Result = StringListAlloc();
    LgIndex_t     StartIndex;
    LgIndex_t     EndIndex;
    for (StartIndex = EndIndex = 0; Result != NULL; EndIndex++)
    {
        /* end of sub-string ? */
        if (String[EndIndex] == '\n' || String[EndIndex] == '\0')
        {
            /* extract the sub-string and append it to the string list */
            LgIndex_t Length = EndIndex - StartIndex;
            char*     SubString = ALLOC_ARRAY(Length + 1, char, "sub string");
            if (SubString != NULL)
            {
                CopySubString(SubString, String, StartIndex, Length);
                StringListAppendString(Result, SubString);

                FREE_ARRAY(SubString, "sub string");

                if (String[EndIndex] != '\0')
                    StartIndex = EndIndex + 1;
                else
                    break; /* nothing left to scan */
            }
            else
            {
                /* memory allocation failure: bail out */
                StringListDealloc(&Result);
                Result = NULL;
                break;
            }
        }
    }

    ENSURE(Result == NULL || StringListValid(Result));
    return Result;
}


/*
 * Return a 'C' string array representation of the string list.
 * Caller is responsible for de-allocating the result.
 */
char** StringListToArray(StringList_pa StringList)
{
    REQUIRE(StringListValid(StringList));

    char** Result = static_cast<char**>(ArrayListToArray(reinterpret_cast<ArrayList_pa>(StringList),
                                                         StringListItemDuplicator, 0));

    ENSURE(Result == NULL || VALID_REF(Result));
    return Result;
}



/*
 * Create a string list from the 'C' string array. The string array
 * is copied and therefore owned and managed by the caller.
 */
StringList_pa StringListFromArray(char const** StringArray,
                                  LgIndex_t    Count)
{
    REQUIRE((Count == 0 && StringArray == NULL) ||
            (Count >= 1 && VALID_REF(StringArray)));

    StringList_pa Result = reinterpret_cast<StringList_pa>(ArrayListFromArray(static_cast<void*>(StringArray),
                                                                              Count, ArrayListType_CharPtr,
                                                                              StringListItemDuplicator, 0));

    ENSURE(Result == NULL || StringListValid(Result));
    return Result;
}



#define ISJOINCHAR(c) ((c == ';') || (c == '+'))

static void SkipWhiteSpaceOrComma(char const** CPtr)
{
    REQUIRE(VALID_REF(CPtr) && VALID_REF(*CPtr));
    while (ISWHITESPACE(**CPtr) || (**CPtr == ','))
        (*CPtr)++;
}

/*
 * Obtain the next sub-string.  This can be of the form:
 *
 *  [del]any-character-sequence[del]
 *
 *            or
 *
 *   limited-character-sequence
 *
 *  where a limited-character-sequence cannot contain
 * any of the following characters:  +;,<space>
 *
 */
static Boolean_t GetNextSubString(char const** OriginalCPtr,
                                  char**       NextSubString)
{
    REQUIRE(VALID_REF(OriginalCPtr) && (VALID_REF(*OriginalCPtr)));
    REQUIRE(VALID_REF(NextSubString));

    Boolean_t IsOk = TRUE;

    *NextSubString = NULL;

    char const* CPtr = *OriginalCPtr;
    SkipWhiteSpaceOrComma(&CPtr);

    char InsideDelimiter = '\0';
    if (*CPtr == '"' || *CPtr == '\'')
    {
        InsideDelimiter = *CPtr;
        CPtr++;
    }

    char const* CStart = CPtr;

    while (*CPtr &&
           ((InsideDelimiter && (*CPtr != InsideDelimiter)) ||
            (!InsideDelimiter && (*CPtr != ',')       &&
             !ISJOINCHAR(*CPtr)  &&
             !ISWHITESPACE(*CPtr))))
    {
        if (InsideDelimiter  &&
            (*CPtr == '\\')  &&
            (*(CPtr + 1) == InsideDelimiter))
            CPtr += 2;
        else
            CPtr++;
    }

    if (InsideDelimiter && (*CPtr != InsideDelimiter))
        IsOk = FALSE;


    if (IsOk && CStart < CPtr)
    {
        size_t StrLen = static_cast<size_t>(CPtr - CStart);
        *NextSubString = ALLOC_ARRAY(StrLen + 1, char, "GetNextSubString: NextSubString");
        if (*NextSubString)
        {
            char* NPtr = *NextSubString;
            /*
             * Don't just copy the string because escaped delimiters need to have
             * the escape removed...
             */
            while (CStart < CPtr)
            {
                if ((*CStart == '\\') && (*(CStart + 1) == InsideDelimiter))
                    CStart++;
                *NPtr++ = *CStart++;
            }
            *NPtr = '\0';
        }
        else
            IsOk = FALSE;
    }

    if (IsOk)
    {
        if (InsideDelimiter)
            CPtr++;
        SkipWhiteSpaceOrComma(&CPtr);
        *OriginalCPtr = CPtr;
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}




/*
 * Return a string list representation of a compound string.
 *
 * The compound String parameter has the following form:
 *
 * [del]<character-sequence>[del] [GroupJoinCharacter] [del]<character-sequence>[del] [GroupJoinCharacter] .....
 *                       or
 * <nospace-character-sequence> <nospace-character-sequence> ...
 *
 * where:
 *   [del] is an optional single quote or a double quote.  [del] must be used
 *   if <character-sequence> contains spaces, commas, or the plus symbol.
 *
 *   GroupJoinCharacter can be either a "+" or a ";"
 *
 * The GroupJoinCharacter symbol is used to separate character sequences that
 * are to be grouped together. If the GroupJoinCharacter symbol is omitted then
 * a new group is started.
 *
 * Internally, the original string is converted to a list of strings where
 * each string uses newlines to separate one sub-string from the next.
 */
StringList_pa StringListFromCompound(char const* String)
{
    REQUIRE(VALID_REF(String));
    SkipWhiteSpaceOrComma(&String);
    REQUIRE(!ISJOINCHAR(*String));

    Boolean_t IsOk = TRUE;

    /* extract character sequences */
    StringList_pa Result = StringListAlloc();
    char const*   CPtr   = String;

    char* CurString = NULL;
    while (IsOk && *CPtr != '\0')
    {
        char*     NextSubString = NULL;
        Boolean_t WantsToJoin   = FALSE;

        if (ISJOINCHAR(*CPtr))
        {
            WantsToJoin = TRUE;
            CPtr++;
            SkipWhiteSpaceOrComma(&CPtr);
        }

        IsOk = GetNextSubString(&CPtr,
                                &NextSubString);

        if (IsOk)
        {
            /*
             * Tack on the sub-string to the running string.
             */
            if (WantsToJoin)
                TackOnChar(&CurString, '\n');
            if (NextSubString != NULL && strlen(NextSubString) != 0)
                IsOk = TackOnString(&CurString, NextSubString, FALSE, FALSE);
            else if (CurString == NULL)
                CurString = DupString(dontTranslate(""));
        }

        if (NextSubString != NULL)
            FREE_ARRAY(NextSubString, "StringListFromCompound: NextSubString");

        /*
         * If this is the end of processing or if the next character is
         * not a join character then add the current string to the stringlist.
         */

        if (IsOk && !ISJOINCHAR(*CPtr))
        {
            StringListAppendString(Result, CurString);
            if (CurString != NULL)
                FREE_ARRAY(CurString, "current string");
            CurString = NULL;
        }
    }

    if (CurString != NULL)
        FREE_ARRAY(CurString, "current string");

    if (!IsOk)
        StringListDealloc(&Result);

    ENSURE(Result == NULL || StringListValid(Result));
    return Result;
}


/*
 * Return a compound string representation of a string list.
 *
 * One common usage in Tecplot:
 *   The $!OpenLayout command in tecplot has the sub-option
 *   ALTDATALOADINSTRUCTIONS that has the form:
 *     '"instr-string1" [GroupJoinCharacter] "instr-string2" [+] ...'
 */
char *StringListToCompound(StringList_pa StringList,
                           char          GroupJoinCharacter,
                           char const*   CharsToEscape)
{
    REQUIRE(StringListValid(StringList));
    REQUIRE(StringListCount(StringList) >= 1);
    REQUIRE(ISJOINCHAR(GroupJoinCharacter));
    REQUIRE(VALID_REF(CharsToEscape));

    char* Result = NULL;

    Boolean_t IsOk = TRUE;
    LgIndex_t Index;
    LgIndex_t Count;
    for (Index = 0, Count = StringListCount(StringList), IsOk = TRUE;
         Index < Count && IsOk;
         Index++)
    {
        char* String = StringListGetString(StringList, Index);

        if (String != NULL && strlen(String) != 0)
        {
            char*       CStart = NULL;
            char*       CEnd = NULL;
            char*       EscapedString = NULL;
            char const* EscChar = NULL;
            char*       StrChar = NULL;

            /* First scan the string and escape any specified characters.  */
            /* Note that the Escape sequence is a double backslash because */
            /* it the first escape escapes the escape for variable usage.  */
            for (StrChar = String; *StrChar != '\0'; StrChar++)
            {
                for (EscChar = CharsToEscape; *EscChar != '\0'; EscChar++)
                    if (*StrChar == *EscChar)
                    {
                        IsOk = TackOnChar(&EscapedString, '\\');
                        IsOk = TackOnChar(&EscapedString, '\\');
                        break;
                    }
                IsOk = TackOnChar(&EscapedString, *StrChar);
            }

            CEnd = EscapedString;
            while (IsOk && *CEnd != '\0')
            {
                int Len = 0;
                CStart = CEnd;
                while (*CEnd != '\0' && *CEnd != '\n')
                {
                    Len++;
                    if (*CEnd == '"')
                        Len++;
                    CEnd++;
                }

                char* TString = ALLOC_ARRAY(Len + 4, char, "temp compound sub-string");
                if (TString != NULL)
                {
                    /* prepend the new string with either   */
                    /* a space character or the plus symbol */
                    if (CStart == EscapedString)
                    {
                        if (Index != 0)
                            IsOk = TackOnChar(&Result, ' ');
                    }
                    else
                    {
                        IsOk = TackOnChar(&Result, GroupJoinCharacter);
                    }

                    /* stuff TString and append the new string */
                    char* TStr = TString;
                    *TStr++ = '"';
                    while (CStart != CEnd)
                    {
                        if (*CStart == '"')
                            *TStr++ = '\\';
                        *TStr++ = *CStart++;
                    }
                    *TStr++ = '"';
                    *TStr = '\0';

                    TackOnString(&Result, TString, FALSE, FALSE);
                    FREE_ARRAY(TString, "StringListToCompound");
                    TString = NULL;
                    if (*CEnd)
                        CEnd++;
                }
                else
                {
                    IsOk = FALSE;
                }
            }

            if (EscapedString != NULL)
                FREE_ARRAY(EscapedString, "escaped string");
        }
        else
        {
            /* a null pointer or length of zero indicates an empty sub-string */
            if (Index == 0)
                TackOnString(&Result, "\"\"", FALSE, FALSE);
            else
                TackOnString(&Result, " \"\"", FALSE, FALSE);
        }

        if (String != NULL)
            FREE_ARRAY(String, "string list item");
    }

    if (!IsOk)
    {
        if (Result != NULL)
        {
            FREE_ARRAY(Result, "StringListToCompound");
            Result = NULL;
        }
    }

    ENSURE(Result == NULL || VALID_REF(Result));
    return Result;
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
