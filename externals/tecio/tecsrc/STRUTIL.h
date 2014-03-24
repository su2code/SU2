#ifndef TECPLOT_STRUTIL
#define TECPLOT_STRUTIL
/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/

#if defined EXTERN
#undef EXTERN
#endif
#if defined STRUTILMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

#include <string>

namespace tecplot
{
namespace strutil
{
class Scanner;
}
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif



EXTERN void FormatStringBufferCleanup(void);

/*
 * This is a helper function for FormatString or any other functions that want
 * to format a string based on a format string followed by a set of arguments.
 * See FormatString or ErrMsg functions for example usage.
 *
 * @param Format
 *   C format string.
 * @param Arguments
 *   Variable argument list already fetched using va_start().
 *
 * @return
 *   Allocated string with the formatted string or NULL if it failed.
 */
EXTERN char *vFormatString(const char *Format,
                           va_list     Arguments);

/**
 * Formats a string using the specified C format string.
 *
 * @param Format
 *   C format string.
 * @param ...
 *   Any arguments needed by the C format string.
 *
 * @return
 *   Allocated string with the formatted string or NULL if it failed. The
 *   client is responsible for deallocating the resource.
 */
EXTERN char *FormatString(tecplot::strutil::TranslatedString Format,
                          ...); /* 0 or more variable arguments */

/**
 * Formats a string using the specified C format string and places the result
 * in the string buffer.
 *
 * @param Buffer
 *   String buffer to receive the formatted string.
 * @param Format
 *   C format string.
 * @param ...
 *   Any arguments needed by the C format string.
 *
 * @return
 *   Upon successful return, these functions return the number of characters
 *   printed, not including the trailing '\0' used to end output to strings. If
 *   unsuccessful -1 is returned.
 */
EXTERN int FormatString(std::string&                       Buffer,
                        tecplot::strutil::TranslatedString Format
                        ...); /* 0 or more variable arguments */
EXTERN char *DupString(tecplot::strutil::TranslatedString String);
EXTERN void CopySubString(char       *Target,
                          const char *Source,
                          int         Index,
                          int         Count);

#if !defined MSWIN

EXTERN void ReplaceCharInString(char *S,
                                short OldChar,
                                short NewChar);
#endif

EXTERN void MakeStringLowerCase(char *str);
EXTERN void MakeStringUpperCase(char *str);
EXTERN char *TrimLeadAndTrailSpaces(char *String);
EXTERN char *StringFlushLeft(char *String);
EXTERN char *StringTruncate(char      *String,
                            LgIndex_t MaxLength);
EXTERN char *StringTrimAndTruncate(char      *String,
                                   LgIndex_t MaxLength);

#ifndef MSWIN
EXTERN StringList_pa LineBreakString(const char *String,
                                     UInt32_t    WrapMargin);
#endif

EXTERN void RemoveSeparator(const char **CPtr);
EXTERN void SkipWhiteSpace(const char **CPtr);
EXTERN void SkipNonWhiteSpace(char **CPtr);
EXTERN const char *ustrstr(const char *s1,
                           const char *s2);
EXTERN int  ustrncmp(const char *s1,
                     const char *s2,
                     size_t      Len);
EXTERN int  ustrcmp(const char *s1,
                    const char *s2);

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/* public access */
/* InternalResetString should not be used directly (use ResetString macro) */
#if !defined NO_ASSERTS
EXTERN Boolean_t InternalResetString(char       **SBase,
                                     const char  *NewString,
                                     Boolean_t    IssueErrMsg,
                                     const char  *FileName,
                                     int          LineNumber);
#  define ResetString(SBase, NewString, IssueErrMsg) InternalResetString( \
                                                       SBase, \
                                                       NewString, \
                                                       IssueErrMsg, \
                                                       __FILE__, __LINE__)
#else
EXTERN Boolean_t InternalResetString(char       **SBase,
                                     const char  *NewString,
                                     Boolean_t    IssueErrMsg);
#  define ResetString(SBase, NewString, IssueErrMsg) InternalResetString( \
                                                       SBase, \
                                                       NewString, \
                                                       IssueErrMsg)
#endif

EXTERN Boolean_t ScanForString(tecplot::strutil::Scanner &scanner,
                               std::string               &DestString,
                               Boolean_t                  GrabEntireStringIncludingDelimiters);
EXTERN Boolean_t TackOnString(char      **SBase,
                              const char *StringToAdd,
                              Boolean_t   DeleteStringToAdd,
                              Boolean_t   ConvertNewlineToAscii);
EXTERN Boolean_t TackOnConstString(char      **SBase,
                                   const char *StringToAdd,
                                   Boolean_t   ConvertNewlineToAscii);
EXTERN Boolean_t TackOnChar(char     **SBase,
                            char       CharToAdd);
EXTERN Boolean_t ReplaceNewlineWithBackslashN(char **String);
EXTERN Boolean_t ReplaceBackslashNWithNewline(char **S);

EXTERN Boolean_t EscapeOutDelimitersInString(char **S,
                                             char Delimiter);
EXTERN Boolean_t ScanForSymbol(tecplot::strutil::Scanner  &scanner,
                               char                        Symbol,
                               Boolean_t                   OnlySkipWhiteSpace);


/* Newline Delimited Strings */
EXTERN char *ConvertStringToNewlineDelimitedString(const char *OriginalString);
EXTERN char *ConvertNewlineDelimitedStringToQuotedString(const char *NewlineDelimitedString,
                                                         Boolean_t   SeparateInstructionsWithPlus);



EXTERN char *InsertNameAtPlaceHolder(char *BaseString,
                                     char *NameToInsert);
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined NO_ASSERTS
#endif /* !NO_ASSERTS */
#endif //TECPLOTKERNEL

inline char* EndOfString(char* str)
{
    return str + strlen(str);
}
inline char const* EndOfString(char const* str)
{
    return str + strlen(str);
}

#endif
