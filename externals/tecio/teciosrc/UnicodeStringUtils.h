 #pragma once
#include "ThirdPartyHeadersBegin.h"
#  include <string>
#include "ThirdPartyHeadersEnd.h"
 #if defined MAKEARCHIVE
 #define common_strutil_API 
 #else
#include "common_strutil_Exports.h" 
 #endif
namespace tecplot {
 #if defined MSWIN
common_strutil_API std::wstring utf8ToWideString(std::string const& utf8String); common_strutil_API std::string wideStringToUtf8(const wchar_t* wideString);
 #endif
 #if !defined NO_THIRD_PARTY_LIBS
common_strutil_API bool isUtf8(char const* stringToTest); common_strutil_API bool isUtf8Checked(char const* stringToTest, std::string const& sourceFile, int line); common_strutil_API char* strToUtf8(std::string const& inputString);
 #endif 
}
 #if !defined NO_THIRD_PARTY_LIBS
 #define REQUIRE_VALID_UTF8_STR(str) REQUIRE(tecplot::isUtf8Checked((str), __FILE__, __LINE__))
 #define CHECK_VALID_UTF8_STR(str)   ___478(tecplot::isUtf8Checked((str),   __FILE__, __LINE__))
 #define ENSURE_VALID_UTF8_STR(str)  ENSURE(tecplot::isUtf8Checked((str),  __FILE__, __LINE__))
 #define REQUIRE_VALID_UTF8_STR_OR_NULL(str) REQUIRE((str) == nullptr || tecplot::isUtf8Checked((str), __FILE__, __LINE__))
 #define CHECK_VALID_UTF8_STR_OR_NULL(str)   ___478((str)   == nullptr || tecplot::isUtf8Checked((str), __FILE__, __LINE__))
 #define ENSURE_VALID_UTF8_STR_OR_NULL(str)  ENSURE((str)  == nullptr || tecplot::isUtf8Checked((str), __FILE__, __LINE__))
 #endif 
