 #pragma once
 #if defined EXTERN
 #undef EXTERN
 #endif
 #if defined ___3862
 #define EXTERN
 #else
 #define EXTERN extern
 #endif
#include "ThirdPartyHeadersBegin.h"
#include <string>
#include <set>
#include <vector>
#include "ThirdPartyHeadersEnd.h"
namespace tecplot { class ___3441; } EXTERN void ___1476(void); EXTERN char *___4410(const char *Format, va_list     ___93); EXTERN char *___1475(tecplot::___4218 Format, ...); EXTERN int ___1475(std::string&                       ___417, tecplot::___4218 Format ...); EXTERN char *___1135(tecplot::___4218 ___3813); EXTERN void ___678(char       *___3946, const char *___3642, int         ___1926, int         ___684);
 #if !defined MSWIN
EXTERN void ___3352(char *S, short ___2871, short ___2696);
 #endif
EXTERN void ___2335(char *str); EXTERN void ___2336(char *str); EXTERN char *___4225(char *___3813); EXTERN char *___3818(char *___3813); EXTERN char *___3859(char *___3813, int   ___2376); char* ___3858( char* ___3813, int   ___2376);
 #ifndef MSWIN
EXTERN ___3839 ___2245(const char *___3813, uint32_t    ___4477);
 #endif
EXTERN void SkipSeparator(char const** ___685); EXTERN void ___3609(char const** ___685); EXTERN void ___3608(char const** ___685); EXTERN const char *ustrstr(const char *___3430, const char *___3431); EXTERN int  ustrncmp(const char *___3430, const char *___3431, size_t      ___2223); EXTERN int  ustrcmp(const char *___3430, const char *___3431);
 #if !defined NO_ASSERTS
EXTERN ___372 ___1982(char       **___3433, const char  *___2700, ___372    ___2063, const char  *___1395, int          ___2262);
 #  define ___3355(___3433, ___2700, ___2063) ___1982( \
 ___3433, \
 ___2700, \
 ___2063, \
 __FILE__, __LINE__)
 #else
EXTERN ___372 ___1982(char       **___3433, const char  *___2700, ___372    ___2063);
 #  define ___3355(___3433, ___2700, ___2063) ___1982( \
 ___3433, \
 ___2700, \
 ___2063)
 #endif
EXTERN ___372 ___3941( char**      ___3433, char const* ___3856, ___372   ___958, ___372   convertNewlineToEscSeq); EXTERN ___372 ___3940( char**      ___3433, char const* ___3856, ___372   convertNewlineToEscSeq); EXTERN ___372 ___3939(char     **___3433, char       ___476); EXTERN ___372 ___3353(char **___3813); EXTERN ___372 ___3351(char **S); EXTERN ___372 ___1187(char **S, char ___959); EXTERN ___372 ___3439(tecplot::___3441  &___3440, char                        ___3915, ___372                   ___2875); EXTERN char *___656(const char *___2884); EXTERN char *___654(const char *___2699, ___372   ___3474); EXTERN char *___1954(char *___339, char *___2687); inline char* EndOfString(char* str) { return str + strlen(str); }; inline char const* EndOfString(char const* str) { return str + strlen(str); };
