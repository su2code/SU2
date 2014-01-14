/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/


#if !defined Q_UNICODE_H_
# define Q_UNICODE_H_

#if defined EXTERN
#undef EXTERN
#endif
#if defined Q_UNICODEMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

namespace tecplot
{
namespace strutil
{

// functions
EXTERN Boolean_t IsValidUtf8LeadByte(Byte_t ch);
EXTERN Boolean_t IsValidUtf8ContinuingByte(Byte_t ch);
EXTERN Boolean_t IsValidUtf8Byte(Byte_t ch);

EXTERN Boolean_t IsPrintable8BitAsciiChar(wchar_t wChar);

EXTERN Boolean_t IsValidUtf8String(const char *str);
EXTERN Boolean_t ShouldConvertWideStringToUtf8String(const wchar_t *str);
EXTERN void InitTranslatedStrings();
EXTERN void CleanUpTranslatedStrings();

EXTERN Boolean_t IsNullOrZeroLengthString(const char *S);
EXTERN Boolean_t IsNullOrZeroLengthString(tecplot::strutil::TranslatedString TS);

EXTERN Boolean_t IsEmptyString(const char *S);
EXTERN Boolean_t IsEmptyString(tecplot::strutil::TranslatedString S);
EXTERN Boolean_t IsEmptyString(const wchar_t* S);

EXTERN std::string AsciiToUtf8String(unsigned char asciiChar);

#if defined MSWIN

EXTERN std::string  LookUpTranslation(std::string& strEnglish);
EXTERN void MsWinInitTranslatedStrings();

EXTERN std::string    WStringToString(std::wstring str);
EXTERN std::wstring   StringToWString(std::string str);

EXTERN std::wstring   MultiByteToWideChar(const char*  Utf8Str,
                                          unsigned int CodePage);

EXTERN std::string    WideCharToMultiByte(const wchar_t* WideStr,
                                          unsigned int   CodePage);

// Conversion
EXTERN std::string    WideCharToUtf8(const wchar_t* str);
EXTERN std::wstring   Utf8ToWideChar(const char *str);
EXTERN char *getenv(const char *str);

#endif

}
}

#endif
