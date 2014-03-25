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

#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"


using namespace std;

namespace tecplot
{
namespace strutil
{

#if defined MSWIN && !defined TECPLOTKERNEL
/**
 * Stub function for non-TECPLOTKERNEL
 */
string LookUpTranslation(string& str)
{
    return string(str);
}
#endif

/**
 * Convenience function for creating Utf8 string translations.
 *
 * @param str
 *     String to translate.
 *
 * @return
 *     A new Utf8 translated string.
 */
static inline string* createUtf8StringTranslation(string& str)
{
#if defined MSWIN
    string *result = new string(LookUpTranslation(str));
#else
    string *result = new string(str);
#endif
    ENSURE(VALID_REF(result));
    return result;
}

#if defined MSWIN
/**
 * Convenience function for creating wide string translations.
 *
 * @param str
 *     String to translate.
 *
 * @return
 *     A new wide translated string.
 */
static inline wstring* createWideStringTranslation(string& str)
{
    wstring *result = new wstring;
    *result = StringToWString(LookUpTranslation(str));

    ENSURE(VALID_REF(result));
    return result;
}
#endif

#if defined MSWIN
/**
 * Convenience function for creating wide string with the given mode.
 *
 * @param mode
 *     Indicates if this string is to be translated or not.
 * @param str
 *     String to translate.
 *
 * @return
 *     A new wide translated string.
 */
static inline wstring* createWideString(TranslatedString::Mode mode,
                                        string&                str)
{
    REQUIRE(mode == TranslatedString::DoTranslate || mode == TranslatedString::DontTranslate);

    wstring* result;
    if (mode == TranslatedString::DoTranslate)
        result = createWideStringTranslation(str);
    else
        result = new wstring(StringToWString(str));

    return result;
}
#endif

/**
 */
void TranslatedString::init(TranslatedString::Mode mode,
                            const char*            str,
                            const char*            translatorNotes)
{
    REQUIRE(mode == DoTranslate || mode == DontTranslate);
    REQUIRE(VALID_REF_OR_NULL(str));
    REQUIRE(VALID_REF_OR_NULL(translatorNotes));

    m_mode   = mode;
    m_isNull = (str == NULL);
    if (!m_isNull)
        m_string = str;
    m_utf8String = NULL; // ...on demand resource
#if defined MSWIN
    m_wideString = NULL; // ...on demand resource
#endif
}

/**
 */
TranslatedString::TranslatedString()
{
    init(DontTranslate, (const char*)NULL, (const char*)NULL);
    ENSURE(this->isValid());
}

/**
 */
TranslatedString TranslatedString::null()
{
    return dontTranslate(NULL);
}

/**
 */
TranslatedString::TranslatedString(TranslatedString::Mode mode,
                                   const char*            str,
                                   const char*            translatorNotes)
{

    REQUIRE(mode == DoTranslate || mode == DontTranslate);
    REQUIRE(VALID_REF_OR_NULL(str));
    REQUIRE(VALID_REF_OR_NULL(translatorNotes));

    init(mode, str, translatorNotes);
    ENSURE(this->isValid());
}

/**
 */
TranslatedString::~TranslatedString()
{
    delete m_utf8String;
#if defined MSWIN
    delete m_wideString;
#endif
}

#if !defined NO_ASSERTS
/**
 */
bool TranslatedString::isValid() const
{
    CHECK(IMPLICATION(m_isNull, m_string.length() == 0));
#if 0
    /* TODO(DTO/RMS/CAM): 11/2/2007
    *   Code currently exists in Tecplot where we call translate() on a
    *   variable. This seems wrong and at times (PleaseWait() in
                                                 *   particular) the variable passed is a NULL pointer which causes
    *   this assertion to fail. There is not enough time before v11.2
    *   release to remove all translate() calls to non-literal strings so
    *   we'll have to do this as a post release cleanup. For now just
    *   deactivate this assertion.
    */
    CHECK(IMPLICATION(m_isNull, m_mode == DontTranslate));
#endif

    return true;
}
#endif

/**
 */
bool TranslatedString::isNull() const
{
    INVARIANT(this->isValid());
    return m_isNull;
}

/**
 */
bool TranslatedString::isNullOrZeroLength() const
{
    INVARIANT(this->isValid());
    return m_isNull || m_string.length() == 0;
}

/**
 */
const char* TranslatedString::c_str()
{
    INVARIANT(this->isValid());

    const char* result = NULL;
    if (!isNull())
    {
        if (m_mode == DoTranslate)
        {
            if (m_utf8String == NULL)
                m_utf8String = createUtf8StringTranslation(m_string);
            result = m_utf8String->c_str();
        }
        else // ...if we aren't translating don't bother creating another Utf8 copy just use m_string
            result = m_string.c_str();
    }

    ENSURE(result == NULL || VALID_REF(result));
    return result;
}

#if defined MSWIN
/**
 */
const wchar_t *TranslatedString::c_wstr()
{
    INVARIANT(this->isValid());

    const wchar_t *result = NULL;
    if (!isNull())
    {
        if (m_wideString == NULL)
            m_wideString = createWideString(m_mode, m_string);
        result = m_wideString->c_str();
    }

    ENSURE(result == NULL || VALID_REF(result));
    return result;
}
#endif

/**
 */
TranslatedString::operator string()
{
    INVARIANT(this->isValid());
    REQUIRE(!isNull());

    string* result;
    if (m_mode == DoTranslate)
    {
        if (m_utf8String == NULL)
            m_utf8String = createUtf8StringTranslation(m_string);
        result = m_utf8String;
    }
    else // ...if we aren't translating don't bother creating another Utf8 copy just use m_string
        result = &m_string;

    return *result;
}

#if defined MSWIN
/**
 */
TranslatedString::operator wstring()
{
    INVARIANT(this->isValid());
    REQUIRE(!isNull());

    if (m_wideString == NULL)
        m_wideString = createWideString(m_mode, m_string);

    return *m_wideString;
}
#endif

/**
 */
TranslatedString& TranslatedString::operator =(const TranslatedString& other)
{
    REQUIRE(other.isValid());

    if (this != &other) // ...only perform if not self assignment
    {
        m_mode       = other.m_mode;
        m_isNull     = other.m_isNull;
        m_string     = other.m_string;
        m_utf8String = (other.m_utf8String != NULL ? new string(*other.m_utf8String) : NULL);
#if defined MSWIN
        m_wideString = (other.m_wideString != NULL ? new wstring(*other.m_wideString) : NULL);
#endif
    }

    ENSURE(this->isValid());
    return *this;
}

/**
 */
TranslatedString::TranslatedString(const TranslatedString& other)
{
    REQUIRE(other.isValid());

    m_mode       = other.m_mode;
    m_isNull     = other.m_isNull;
    m_string     = other.m_string;
    m_utf8String = (other.m_utf8String != NULL ? new string(*other.m_utf8String) : NULL);
#if defined MSWIN
    m_wideString = (other.m_wideString != NULL ? new wstring(*other.m_wideString) : NULL);
#endif

    ENSURE(this->isValid());
}

/**
 */
TranslatedString TranslatedString::translate(const char* str,
                                             const char* translatorNotes)
{
    REQUIRE(VALID_REF_OR_NULL(str));
    REQUIRE(VALID_REF_OR_NULL(translatorNotes));

    return TranslatedString(DoTranslate, str, translatorNotes);
}

/**
 */
TranslatedString TranslatedString::dontTranslate(const char* str)
{
    REQUIRE(VALID_REF_OR_NULL(str));

    return TranslatedString(DontTranslate, str, NULL);
}

}
}
