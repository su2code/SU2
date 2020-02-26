 #ifndef TECPLOT_CHARTYPE
 #define TECPLOT_CHARTYPE
#include <cctype>
#include <locale>
 #if defined EXTERN
 #undef EXTERN
 #endif
 #if defined ___477
 #define EXTERN
 #else
 #define EXTERN extern
 #endif
namespace tecplot { template <typename CHAR_TYPE> inline bool isspace(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return ::isspace(static_cast<int>(ch));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return std::use_facet<std::ctype<wchar_t> >(___2313).is(std::ctype_base::space, static_cast<wchar_t>(ch));
 #endif
} template <typename CHAR_TYPE> inline bool ___2050(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return ::___2050(static_cast<int>(ch));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return std::use_facet<std::ctype<wchar_t> >(___2313).is(std::ctype_base::print, static_cast<wchar_t>(ch));
 #endif
} template <typename CHAR_TYPE> inline bool ___2009(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return ::___2009(static_cast<int>(ch));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return std::use_facet<std::ctype<wchar_t> >(___2313).is(std::ctype_base::cntrl, static_cast<wchar_t>(ch));
 #endif
} template <typename CHAR_TYPE> inline bool isupper(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return ::isupper(static_cast<int>(ch));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return std::use_facet<std::ctype<wchar_t> >(___2313).is(std::ctype_base::upper, static_cast<wchar_t>(ch));
 #endif
} template <typename CHAR_TYPE> inline bool ___2032(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return ::___2032(static_cast<int>(ch));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return std::use_facet<std::ctype<wchar_t> >(___2313).is(std::ctype_base::lower, static_cast<wchar_t>(ch));
 #endif
} template <typename CHAR_TYPE> inline bool ___1998(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return ::___1998(static_cast<int>(ch));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return std::use_facet<std::ctype<wchar_t> >(___2313).is(std::ctype_base::alpha, static_cast<wchar_t>(ch));
 #endif
} template <typename CHAR_TYPE> inline bool ___2012(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return ::___2012(static_cast<int>(ch));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return std::use_facet<std::ctype<wchar_t> >(___2313).is(std::ctype_base::digit, static_cast<wchar_t>(ch));
 #endif
} template <typename CHAR_TYPE> inline bool ___2052(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return ::___2052(static_cast<int>(ch));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return std::use_facet<std::ctype<wchar_t> >(___2313).is(std::ctype_base::punct, static_cast<wchar_t>(ch));
 #endif
} template <typename CHAR_TYPE> inline bool ___2084(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return ::___2084(static_cast<int>(ch));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return std::use_facet<std::ctype<wchar_t> >(___2313).is(std::ctype_base::xdigit, static_cast<wchar_t>(ch));
 #endif
} template <typename CHAR_TYPE> inline bool ___1997(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return ::___1997(static_cast<int>(ch));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return std::use_facet<std::ctype<wchar_t> >(___2313).is(std::ctype_base::alnum, static_cast<wchar_t>(ch));
 #endif
} template <typename CHAR_TYPE> inline bool ___2024(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return ::___2024(static_cast<int>(ch));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return std::use_facet<std::ctype<wchar_t> >(___2313).is(std::ctype_base::graph, static_cast<wchar_t>(ch));
 #endif
} template <typename CHAR_TYPE> inline CHAR_TYPE toupper(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return static_cast<CHAR_TYPE>(::toupper(static_cast<int>(ch)));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return static_cast<CHAR_TYPE>( std::use_facet<std::ctype<wchar_t> >(___2313).toupper(static_cast<wchar_t>(ch)));
 #endif
} template <typename CHAR_TYPE> inline CHAR_TYPE tolower(CHAR_TYPE          ch, std::locale const& ___2313 = std::locale::classic()) {
 #if defined ___3893
REQUIRE(___2313 == std::locale::classic()); return static_cast<CHAR_TYPE>(::tolower(static_cast<int>(ch)));
 #else
REQUIRE(std::has_facet<std::ctype<wchar_t> >(___2313)); return static_cast<CHAR_TYPE>( std::use_facet<std::ctype<wchar_t> >(___2313).tolower(static_cast<wchar_t>(ch)));
 #endif
} }
 #endif
