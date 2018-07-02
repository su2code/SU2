#ifndef TECPLOT_CHARTYPE
#define TECPLOT_CHARTYPE
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
#if defined CHARTYPEMODULE
#define EXTERN
#else
#define EXTERN extern
#endif


/**
 */
namespace tecplot
{

/**
 * Same as std::isspace() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline bool isspace(CHAR_TYPE          ch,
                    std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return ::isspace(static_cast<int>(ch));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return std::use_facet<std::ctype<wchar_t> >(loc).is(std::ctype_base::space, static_cast<wchar_t>(ch));
    #endif
}

/**
 * Same as std::isprint() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline bool isprint(CHAR_TYPE          ch,
                    std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return ::isprint(static_cast<int>(ch));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return std::use_facet<std::ctype<wchar_t> >(loc).is(std::ctype_base::print, static_cast<wchar_t>(ch));
    #endif
}

/**
 * Same as std::iscntrl() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline bool iscntrl(CHAR_TYPE          ch,
                    std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return ::iscntrl(static_cast<int>(ch));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return std::use_facet<std::ctype<wchar_t> >(loc).is(std::ctype_base::cntrl, static_cast<wchar_t>(ch));
    #endif
}

/**
 * Same as std::isupper() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline bool isupper(CHAR_TYPE          ch,
                    std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return ::isupper(static_cast<int>(ch));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return std::use_facet<std::ctype<wchar_t> >(loc).is(std::ctype_base::upper, static_cast<wchar_t>(ch));
    #endif
}

/**
 * Same as std::islower() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline bool islower(CHAR_TYPE          ch,
                    std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return ::islower(static_cast<int>(ch));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return std::use_facet<std::ctype<wchar_t> >(loc).is(std::ctype_base::lower, static_cast<wchar_t>(ch));
    #endif
}

/**
 * Same as std::isalpha() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline bool isalpha(CHAR_TYPE          ch,
                    std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return ::isalpha(static_cast<int>(ch));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return std::use_facet<std::ctype<wchar_t> >(loc).is(std::ctype_base::alpha, static_cast<wchar_t>(ch));
    #endif
}

/**
 * Same as std::isdigit() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline bool isdigit(CHAR_TYPE          ch,
                    std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return ::isdigit(static_cast<int>(ch));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return std::use_facet<std::ctype<wchar_t> >(loc).is(std::ctype_base::digit, static_cast<wchar_t>(ch));
    #endif
}

/**
 * Same as std::ispunct() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline bool ispunct(CHAR_TYPE          ch,
                    std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return ::ispunct(static_cast<int>(ch));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return std::use_facet<std::ctype<wchar_t> >(loc).is(std::ctype_base::punct, static_cast<wchar_t>(ch));
    #endif
}

/**
 * Same as std::isxdigit() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline bool isxdigit(CHAR_TYPE          ch,
                     std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return ::isxdigit(static_cast<int>(ch));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return std::use_facet<std::ctype<wchar_t> >(loc).is(std::ctype_base::xdigit, static_cast<wchar_t>(ch));
    #endif
}

/**
 * Same as std::isalnum() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline bool isalnum(CHAR_TYPE          ch,
                    std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return ::isalnum(static_cast<int>(ch));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return std::use_facet<std::ctype<wchar_t> >(loc).is(std::ctype_base::alnum, static_cast<wchar_t>(ch));
    #endif
}

/**
 * Same as std::isgraph() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline bool isgraph(CHAR_TYPE          ch,
                    std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return ::isgraph(static_cast<int>(ch));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return std::use_facet<std::ctype<wchar_t> >(loc).is(std::ctype_base::graph, static_cast<wchar_t>(ch));
    #endif
}

/**
 * Same as std::toupper() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline CHAR_TYPE toupper(CHAR_TYPE          ch,
                         std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return static_cast<CHAR_TYPE>(::toupper(static_cast<int>(ch)));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return static_cast<CHAR_TYPE>(
        std::use_facet<std::ctype<wchar_t> >(loc).toupper(static_cast<wchar_t>(ch)));
    #endif
}

/**
 * Same as std::tolower() with the classic default C locale for portable
 * comparisons.
 */
template <typename CHAR_TYPE>
inline CHAR_TYPE tolower(CHAR_TYPE          ch,
                         std::locale const& loc = std::locale::classic())
{
    #if defined SUNX
    REQUIRE(loc == std::locale::classic());
    return static_cast<CHAR_TYPE>(::tolower(static_cast<int>(ch)));
    #else
    REQUIRE(std::has_facet<std::ctype<wchar_t> >(loc));
    return static_cast<CHAR_TYPE>(
        std::use_facet<std::ctype<wchar_t> >(loc).tolower(static_cast<wchar_t>(ch)));
    #endif
}

}

#endif

