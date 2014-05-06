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

#define Q_MSGMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "ALLOC.h"
#include "ARRLIST.h"

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined ENGINE
#if defined MOTIF
#endif
#if defined MSWIN
#endif
#endif
#endif

#include "STRUTIL.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

using std::string;
using namespace tecplot::strutil;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#define MAXCHARACTERSPERLINE 60
/*
 * Wrap a string so it contains at most CharactersPerLine
 * characters. Embedded newlines are left alone. Spaces
 * following newlines are also left alone.
 */
Boolean_t WrapString(const char  *OldString,
                     char       **NewString)
{
    size_t L;
    if (OldString == NULL)
        return (FALSE);

    /*
     *  Assume Old string has ample spaces.  Will only be
     *  replacing some spaces with newlines and removing
     *  other spaces.  New string can be allocated to be
     *  same length as old string.
     */

    L = strlen(OldString);
    *NewString = ALLOC_ARRAY(L + 1, char, "new error message string");
    if (*NewString == NULL)
        return (FALSE);

    strcpy(*NewString, OldString);

    if (L > MAXCHARACTERSPERLINE)
    {
        char *LineStart = *NewString;
        char *LastWord  = *NewString;
        char *WPtr      = *NewString;
        while (WPtr && (*WPtr != '\0'))
        {
            size_t CurLen;
            /*
             * Find next hard newline.  If there is one befor the
             * line should be chopped then reset the Last Word to
             * be at the first word after the newline.
             */
            WPtr = strchr(LineStart, '\n');
            if (WPtr && ((WPtr - LineStart) < MAXCHARACTERSPERLINE))
            {
                WPtr++;
                while (*WPtr == '\n')
                    WPtr++;
                LineStart = WPtr;
                /*
                 * Skip over trailing spaces.  Special handling to
                 * allow indent after hard newline.
                 */
                while (*WPtr == ' ')
                    WPtr++;
                LastWord  = WPtr;
                continue;
            }
            /*
             *  Find next "word"
             */
            WPtr = strchr(LastWord, ' ');
            if (WPtr != NULL)
            {
                while (*WPtr == ' ')
                    WPtr++;
            }

            if (WPtr == NULL)
            {
                CurLen = strlen(LineStart);
            }
            else
            {
                CurLen = WPtr - LineStart;
            }

            if (CurLen > MAXCHARACTERSPERLINE)
            {
                /*
                 * Line is too long.  Back up to previous
                 * word and replace preceeding space with
                 * a newline.
                 */
                if (LastWord == LineStart)
                {
                    /*
                     * Bad news, line has very long word.
                     */
                    if (WPtr && (*WPtr != '\0'))
                    {
                        *(WPtr - 1) = '\n';
                        LastWord = WPtr;
                    }
                }
                else
                {
                    *(LastWord - 1) = '\n';
                }
                LineStart = LastWord;
            }
            else
                LastWord = WPtr;
        }
    }
    return (TRUE);
}


static void SendWarningToFile(FILE *F,
                              const char *S)
{
    char *S2;
    REQUIRE(VALID_REF(F));
    REQUIRE(VALID_REF(S));
    if (WrapString(S, &S2))
    {
        fprintf(F, "Warning: %s\n", S2);
        FREE_ARRAY(S2, "temp warning string");
    }
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 * Show the warning message. Note that the format string can be the only
 * argument, in which case it is essentially the warning message itself.
 *
 * param format
 *     C format string or a simple message.
 * param ...
 *     Zero or more variable arguments. The number of arguments must correspond
 *     to the placeholders in the format string.
 */
void Warning(TranslatedString format,
             ...) /* zero or more arguments */
{
    REQUIRE(!format.isNull());

    static Boolean_t InWarning = FALSE; /* ...used to prevent recursive deadlock */
    if (!InWarning)
    {
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

        InWarning = TRUE;
        {
            /*
             * Attempt to format the string. Failing that simply use the original format
             * string argument which, if we ran out of memory while formatting, is
             * probably just an warning message saying that we ran out of memory in some
             * previous operation anyway.
             */
            Boolean_t cleanUp = TRUE;

            va_list  Arguments;
            va_start(Arguments, format);
            char* message = vFormatString(format.c_str(), Arguments);
            va_end(Arguments);

            if (message == NULL)
            {
                cleanUp = FALSE; // ...this boolean allows us to "carefully" cast away the const'ness
                message = (char*)format.c_str();
            }

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#ifdef MSWIN
#endif
#if defined UNIXX
#endif
#if defined MSWIN
#endif
#else /* !TECPLOTKERNEL */
            {
                SendWarningToFile(stderr, message);
            }
#endif

            if (cleanUp)
                FREE_ARRAY(message, "message");
        }
        InWarning = FALSE;

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    }
}


static void SendErrToFile(FILE       *File,
                          const char *Msg)
{
    char *formattedMsg;
    REQUIRE(VALID_REF(File));
    REQUIRE(VALID_REF(Msg));
    if (WrapString(Msg, &formattedMsg))
    {
        fprintf(File, "Err: %s\n", formattedMsg);
        FREE_ARRAY(formattedMsg, "temp error string");
    }
    else
        fprintf(File, "Err: %s\n", Msg);
}


/* Fall-back ErrMsg procedure when nothing else works */
static void DefaultErrMsg(const char *Msg)
{
    REQUIRE(VALID_REF(Msg));

#ifdef MSWIN
#ifdef TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
    MessageBox(NULL, Msg, "Error", MB_OK | MB_ICONERROR);
#endif
#else
    SendErrToFile(stderr, Msg);
#endif
}

static void PostErrorMessage(TranslatedString format,
                             va_list          Arguments)
{
    REQUIRE(!format.isNull());

    /*
     * Attempt to format the string. Failing that simply use the original format
     * string argument which, if we ran out of memory while formatting, is
     * probably just an error message saying that we ran out of memory in some
     * previous operation anyway.
     */
    Boolean_t cleanUp = TRUE;
    char* messageString = vFormatString(format.c_str(), Arguments);
    if (messageString == NULL)
    {
        cleanUp = FALSE; // ...this boolean allows us to "carefully" cast away the const'ness
        messageString = (char*)format.c_str();
    }

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#ifdef MSWIN
#endif
#if defined UNIXX
#if !defined ENGINE
#endif
#endif
#else /* !TECPLOTKERNEL */
    {
        DefaultErrMsg(messageString);
    }
#endif

    /* cleanup if we allocated the string */
    if (cleanUp)
        FREE_ARRAY(messageString, "messageString");
}

/*
 * NOTES:
 *   This function is thread safe in that it may be safely called by multiple
 *   threads however when running interactively only the first error message is
 *   queued for display on idle. In batch mode all messages are sent to the
 *   batch log file.
 */
void vErrMsg(TranslatedString format,
             va_list          Arguments)
{
    REQUIRE(!format.isNull());

    static Boolean_t InErrMsg = FALSE; /* ...used to prevent recursive deadlock */
    if (!InErrMsg)
    {
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

        InErrMsg = TRUE;
        {
            PostErrorMessage(format, Arguments);
        }
        InErrMsg = FALSE;

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    }
}

/**
 * Show the error message. Note that the format string can be the only
 * argument, in which case it is essentially the error message itself.
 *
 * @param format
 *   C format string or a simple message.
 * @param ...
 *   Zero or more variable arguments. The number of arguments must correspond
 *   to the placeholders in the format string.
 */
void ErrMsg(TranslatedString format,
            ...) /* zero or more arguments */
{
    REQUIRE(!format.isNull());

    va_list Arguments;
    va_start(Arguments, format);
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
    PostErrorMessage(format, Arguments);
#endif
    va_end(Arguments);
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined ENGINE
#endif
#if !defined ENGINE
#if defined MOTIF
#endif
#if defined MSWIN
#endif
#endif
#if !defined ENGINE
#if defined MOTIF
#endif
#if defined MSWIN
#endif
#endif
#if !defined ENGINE
#if defined MOTIF
#endif /* MOTIF */
#if defined MSWIN
#endif
#if defined MOTIF
#endif /* MOTIF */
#endif
#if !defined ENGINE
#endif
#endif /* TECPLOTKERNEL */

