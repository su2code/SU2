#include "UnicodeStringUtils.h"
#include "CodeContract.h"
#include "ThirdPartyHeadersBegin.h"
 #if !defined NO_THIRD_PARTY_LIBS
#  include <boost/format.hpp>
#  include <cstring>
#  include <QByteArray>
#  include <QDebug>
#  include <QString>
#  include <QTextCodec>
 #endif
 #if defined MSWIN
#  include <windows.h>
 #endif
#include "ThirdPartyHeadersEnd.h"
namespace tecplot {
 #if defined MSWIN
std::wstring utf8ToWideString(std::string const& utf8String) {
 #if !defined NO_THIRD_PARTY_LIBS
REQUIRE_VALID_UTF8_STR(utf8String.c_str());
 #endif
const unsigned int codePage = CP_UTF8; std::string  UTF8str(utf8String); std::wstring ___3861; if (UTF8str.length() > 0) { size_t wLen = static_cast<size_t>(::MultiByteToWideChar(codePage, 0, UTF8str.c_str(), -1, NULL, 0)); if (wLen > 0) { wchar_t* pBuffer = new wchar_t[wLen]; VERIFY(::MultiByteToWideChar(codePage, 0, UTF8str.c_str(), static_cast<int>(UTF8str.length() + 1), pBuffer, static_cast<int>(wLen)) != 0); ___3861 = pBuffer; delete[] pBuffer; } else { ___478(false); } } ENSURE("strResult is a valid wstring"); return ___3861; } std::string  wideStringToUtf8(const wchar_t* wideString) { REQUIRE(VALID_REF(wideString)); const unsigned int codePage = CP_UTF8; std::wstring wString(wideString); std::string ___3861; if (wString.length() > 0) { size_t nLen = static_cast<size_t>(::WideCharToMultiByte(codePage, 0, wString.c_str(), -1, NULL, 0, NULL, NULL)); if (nLen > 0) { char* tmpBuffer = new char[nLen]; VERIFY(::WideCharToMultiByte(codePage, 0, wideString, static_cast<int>(wString.length() + 1), tmpBuffer, static_cast<int>(nLen), NULL, NULL) != 0); ___3861 = tmpBuffer; delete[] tmpBuffer; } else { ___478(___1305); } }
 #if !defined NO_THIRD_PARTY_LIBS
ENSURE_VALID_UTF8_STR(___3861.c_str());
 #endif
return ___3861; }
 #endif 
 #if !defined NO_THIRD_PARTY_LIBS
bool isUtf8(char const* stringToTest) { REQUIRE(VALID_REF(stringToTest)); QByteArray const byteArray(stringToTest); QTextCodec::ConverterState state; QTextCodec* const codec = QTextCodec::codecForName("UTF-8"); ___478(VALID_REF(codec)); codec->toUnicode(byteArray.constData(), byteArray.size(), &state); return state.invalidChars == 0; } bool isUtf8Checked(char const* stringToTest, std::string const& sourceFile  , int line  ) { REQUIRE(VALID_REF(stringToTest)); REQUIRE(!sourceFile.empty()); REQUIRE(line > 0); if (!isUtf8(stringToTest)) { std::string msg = boost::str(boost::format( "Internal error: Invalid UTF-8 string:\n\n" "___3813: \"%1%\"\n" "File: %2%\n" "Line: %3%\n") % stringToTest % sourceFile % line); qDebug() << msg.c_str(); ___478(false); } return true; } char* strToUtf8(std::string const& inputString) { REQUIRE(!inputString.empty()); QByteArray utf8String(QString(inputString.c_str()).toUtf8()); ___478(utf8String.constData()[utf8String.length()] == '\0'); char* ___3358 = new char[utf8String.length() + 1]; std::strncpy(___3358, utf8String.constData(),utf8String.length()); ___3358[utf8String.length()] = '\0'; ENSURE_VALID_UTF8_STR(___3358); return ___3358; }
 #endif 
}
