/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2010 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#ifndef TECPLOT_STRUTIL_TRANSLATEDSTRING
#define TECPLOT_STRUTIL_TRANSLATEDSTRING

#if defined MSWIN
#pragma warning(disable : 4181)
#endif

namespace tecplot
{
namespace strutil
{

/**
 * Class responsible for translating strings for internationalization. This
 * class is used both to perform the translation and to identify which strings
 * are/aren't in need of translation.
 *
 * With the exception of the empty constructor all translated strings are
 * created via static methods or inline helper functions named translate() and
 * dontTranslate(). Translated strings created with a call to translate() are
 * flagged as needing human translation and return the translated value at
 * runtime. Translated strings created with a call to dontTranslate() are
 * flagged as not needing human translation and return the non-translated value
 * at runtime. Examples:
 *
 *   ErrMsg(translate("Can I have %d cookies?", numCookies);
 *   ErrMsg(dontTranslate("%s"), errMsgString);
 *
 * Conversion methods exists for std::string so that they operate well
 * together. Therefore you can pass translated strings to methods or functions
 * that receive std::string or std::string references. Additionally you can
 * perform std::string operations by casting the translated string to a
 * std::string. For example:
 *
 *   if (string(dialogTitle).size() != 0)
 *     ...do something useful
 *
 * We have intentionally not provided conversion methods for "const char*" (an
 * internal representation of the object's string) because of the risk of the
 * client using the pointer after the translated string object goes out of
 * scope.
 *
 * @author David Ossorio
 */
class TranslatedString
{
public:
    /**
     * Enumeration describing available translation modes.
     */
    typedef enum { DontTranslate, DoTranslate } Mode;

    /**
     * Constructs an empty translated string. It is equivalent to calling
     * TranslatedString::null().
     */
    TranslatedString();

    /**
     * Convenience function for creating a NULL translated string.
     *
     * @return
     *     NULL translated string.
     */
    static TranslatedString null();

    /**
     * Creates a translated string object and marks it as a string needing
     * translation.
     *
     * @param str
     *     Character string to translate.
     * @param translatorNotes
     *     Optional notes for the human translator describing the meaning
     *     of the string.
     */
    static TranslatedString translate(const char* str,
                                      const char* translatorNotes = NULL);

    /**
     * Creates a translated string object and marks it as a string not needing
     * translation.
     *
     * @param str
     *     Character string to translate. The str parameter can be a NULL
     *     pointer.
     */
    static TranslatedString dontTranslate(const char* str);

    /**
     * Destructor.
     */
    virtual ~TranslatedString();

    /**
     * Indicates if the object's string is NULL.
     *
     * @return
     *     true if the object's string is NULL, false otherwise.
     */
    virtual bool isNull() const;

    /**
     * Indicates if the object's string is NULL or zero length.
     *
     * @return
     *     true if the object's string is NULL or zero length, false otherwise.
     */
    virtual bool isNullOrZeroLength() const;

    /**
     * Returns the internal representation of the object's string. Use this
     * result carefully. When this object goes out of scope so does this
     * references.
     *
     * @return
     *     Pointer to the internal representation of the object's string.
     */
    virtual const char* c_str();

#if defined MSWIN
    /**
     * Returns the internal representation of the wide character version of the
     * object's string. Use this result carefully. When this object goes out of
     * scope so does this references.
     *
     * @return
     *     Pointer to the internal representation of the object's string.
     */
    virtual const wchar_t* c_wstr();
#endif

    /**
     * Returns a copy of the object's string by this object. The result is
     * translated or not depending on the platform and how it was created.
     *
     * @return
     *     Copy of the object's string.
     */
    virtual operator std::string();

#if defined MSWIN
    /**
     * Returns a copy of the wide character version of the object's string.
     * The result is translated or not depending on the platform and how it was
     * created.
     *
     * @return
     *     Copy of the wide character version of the object's string.
     */
    virtual operator std::wstring();
#endif

    /**
     * Assignment operator.
     */
    virtual TranslatedString& operator =(const TranslatedString& other);

    /**
     * Copy constructor.
     */
    TranslatedString(const TranslatedString& other);

#if !defined NO_ASSERTS
    /**
     * Used only for assertions.
     */
    virtual bool isValid() const;
#endif /* !NO_ASSERTS */

private:
    /**
     * Constructs a translated string. This is declared private to make sure
     * clients use translate() or dontTranslate() so that Tecplot's translation
     * parser can easily extract strings from the source code.
     *
     * @param mode
     *     Indicates if this string is to be translated or not at runtime.
     * @param str
     *     Character string to translate.
     * @param translatorNotes
     *     Optional notes for the human translator describing the meaning
     *     of the string.
     */
    TranslatedString(TranslatedString::Mode mode,
                     const char*            str,
                     const char*            translatorNotes);

    /**
     * Convenience method for initialize a translated string.
     *
     * @param mode
     *     Indicates if this string is to be translated or not at runtime.
     * @param str
     *     Character string to translate.
     * @param translatorNotes
     *     Optional notes for the human translator describing the meaning
     *     of the string.
     */
    void init(TranslatedString::Mode mode,
              const char*            str,
              const char*            translatorNotes);

    /**
     * Private instance data.
     */
    TranslatedString::Mode  m_mode;
    bool                    m_isNull;
    std::string             m_string;
    std::string            *m_utf8String;
#if defined MSWIN
    std::wstring           *m_wideString;
#endif
};

/**
 * Convenience function for creating a translated string object and marking it
 * as a string needing translation.
 *
 * @param str
 *     Character string to translate.
 * @param translatorNotes
 *     Optional notes for the human translator describing the meaning
 *     of the string.
 */
inline TranslatedString translate(const char* str,
                                  const char* translatorNotes = NULL)
{
    return TranslatedString::translate(str, translatorNotes);
}

/**
 * Convenience function for creating a translated string object and marking it
 * as a string not needing translation.
 *
 * @param str
 *     Character string to translate. The str parameter can be a NULL
 *     pointer.
 */
inline TranslatedString dontTranslate(const char* str)
{
    return TranslatedString::dontTranslate(str);
}

/**
 * Convenience function for creating a translated string object and marks it as
 * a string not needing translation.
 *
 * @param str
 *     String to translate.
 */
inline TranslatedString dontTranslate(const std::string& str)
{
    return TranslatedString::dontTranslate(str.c_str());
}

}
}

#endif
