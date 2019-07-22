 #if defined(_MSC_VER)
 #pragma warning (push, 0)
 #pragma warning (disable:4191)
 #pragma warning (disable:4242)
 #pragma warning (disable:4263)
 #pragma warning (disable:4264)
 #pragma warning (disable:4265)
 #pragma warning (disable:4266)
 #pragma warning (disable:4302)
 #pragma warning (disable:4702)
 #pragma warning (disable:4715)
 #pragma warning (disable:4905)
 #pragma warning (disable:4906)
 #pragma warning (disable:4928)
 #if (_MSC_VER > 1600)
#include <codeanalysis\warnings.h>
 #pragma warning (push)
 #pragma warning (disable : ALL_CODE_ANALYSIS_WARNINGS)
 #endif
 #elif defined(__clang__)
 #pragma clang diagnostic push
 #pragma clang diagnostic ignored "-Wall"
 #pragma clang diagnostic ignored "-Wextra"
 #elif defined(__GNUC__)
 #if ((__GNUC__  >= 4) && (__GNUC_MINOR__ >= 6)) || (__GNUC__ >= 5) 
 #pragma GCC diagnostic push
 #pragma GCC diagnostic ignored "-Wall"
 #pragma GCC diagnostic ignored "-Wextra"
 #pragma GCC diagnostic ignored "-Waddress"
 #pragma GCC diagnostic ignored "-Wc++11-compat"
 #pragma GCC diagnostic ignored "-Wchar-subscripts"
 #pragma GCC diagnostic ignored "-Wcomment"
 #pragma GCC diagnostic ignored "-Wenum-compare"
 #pragma GCC diagnostic ignored "-Wformat"
 #pragma GCC diagnostic ignored "-Winit-self"
 #pragma GCC diagnostic ignored "-Wmain"
 #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
 #pragma GCC diagnostic ignored "-Wmissing-braces"
 #pragma GCC diagnostic ignored "-Wnarrowing"
 #pragma GCC diagnostic ignored "-Wnonnull"
 #pragma GCC diagnostic ignored "-Wparentheses"
 #pragma GCC diagnostic ignored "-Wreorder"
 #pragma GCC diagnostic ignored "-Wreturn-type"
 #pragma GCC diagnostic ignored "-Wsequence-point"
 #pragma GCC diagnostic ignored "-Wsign-compare"
 #pragma GCC diagnostic ignored "-Wsizeof-pointer-memaccess"
 #pragma GCC diagnostic ignored "-Wstrict-aliasing"
 #pragma GCC diagnostic ignored "-Wswitch"
 #pragma GCC diagnostic ignored "-Wtrigraphs"
 #pragma GCC diagnostic ignored "-Wuninitialized"
 #pragma GCC diagnostic ignored "-Wunknown-pragmas"
 #pragma GCC diagnostic ignored "-Wunused-function"
 #pragma GCC diagnostic ignored "-Wunused-label"
 #pragma GCC diagnostic ignored "-Wunused-value"
 #pragma GCC diagnostic ignored "-Wunused-variable"
 #pragma GCC diagnostic ignored "-Wvolatile-register-var"
 #pragma GCC diagnostic ignored "-Wclobbered"
 #pragma GCC diagnostic ignored "-Wempty-body"
 #pragma GCC diagnostic ignored "-Wignored-qualifiers"
 #pragma GCC diagnostic ignored "-Wmissing-field-initializers"
 #pragma GCC diagnostic ignored "-Wsign-compare"
 #pragma GCC diagnostic ignored "-Wtype-limits"
 #pragma GCC diagnostic ignored "-Wuninitialized"
 #pragma GCC diagnostic ignored "-Wunused-parameter"
 #pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
 #pragma GCC diagnostic ignored "-Wliteral-suffix"
 #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
 #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
 #endif 
 #endif
