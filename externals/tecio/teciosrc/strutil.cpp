#include "stdafx.h"
#include "MASTER.h"
 #define ___3862
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "ARRLIST.h"
#include "STRLIST.h"
#include "CHARTYPE.h"
#include "STRUTIL.h"
#include "ALLOC.h"
#include "Q_MSG.h"
#include "ThirdPartyHeadersBegin.h"
#include "ThirdPartyHeadersEnd.h"
#include <algorithm>
#include <cctype> 
#include <limits.h>
#include "TranslatedString.h"
#include "stringformat.h"
using std::string; using tecplot::___4217; using tecplot::___1097; using tecplot::___4218;
 #ifdef MSWIN
 # pragma warning (disable : 4786) 
 #endif
 #define           FORMAT_BUFFER_SIZE 16384
static char      *FormatStringBuffer = NULL; static int        FormatStringBufferSize = FORMAT_BUFFER_SIZE; void ___1476(void) { if (FormatStringBuffer != NULL) free(FormatStringBuffer); FormatStringBuffer = NULL; } char *___4410(const char *Format, va_list     ___93) { char *___3359 = NULL; REQUIRE(VALID_REF(Format)); if (FormatStringBuffer == NULL) FormatStringBuffer = (char *)malloc(FormatStringBufferSize); if (FormatStringBuffer != NULL) {
 # if defined MSWIN
ASSERT_ONLY(int nChars =) _vsnprintf(FormatStringBuffer, FormatStringBufferSize, Format, ___93);
 # else
ASSERT_ONLY(int nChars =) vsnprintf(FormatStringBuffer, FormatStringBufferSize, Format, ___93);
 # endif
FormatStringBuffer[FormatStringBufferSize-1] = '\0'; ___478(nChars > 0); ___3359 = ___1135(___1097(FormatStringBuffer)); } ENSURE(VALID_REF(___3359) || ___3359 == NULL); return ___3359; } char *___1475(___4218 Format, ...) { REQUIRE(!Format.___2035()); va_list ___93; va_start(___93, Format); char *___3359 = ___4410(Format.c_str(), ___93); va_end(___93); ENSURE(VALID_REF(___3359) || ___3359 == NULL); return ___3359; } int ___1475(string&           ___417, ___4218  Format ...) { REQUIRE(!Format.___2035()); va_list ___93; va_start(___93, Format); char *FormattedString = ___4410(Format.c_str(), ___93); va_end(___93); int ___3359; if (FormattedString != NULL) { ___417.assign(FormattedString); ___3359 = (int)___417.size(); ___1530(FormattedString, "FormattedString"); } else ___3359 = -1; ENSURE(___3359 == -1 || ___3359 >= 0); return ___3359; } char *___1135(___4218 ___3813) { REQUIRE(VALID_TRANSLATED_STRING(___3813)); char *___3359 = ___23(strlen(___3813.c_str()) + 1, char, "duplicate string"); if (___3359 != NULL) strcpy(___3359, ___3813.c_str()); ENSURE(___3359 == NULL || (VALID_REF(___3359) && strcmp(___3359, ___3813.c_str()) == 0)); return ___3359; } void ___678(char       *___3946, const char *___3642, int         ___1926, int         ___684) { int       ___2224 = 0; REQUIRE(VALID_REF(___3946)); REQUIRE("Target string is sized to accommodate a string who's length " "is at least MIN(strlen(&Source[Index]), Count) characters."); REQUIRE(VALID_REF(___3642)); REQUIRE(0 <= ___1926 && ___1926 <= (int)strlen(___3642)); REQUIRE(___684 >= 0); ___2224 = MIN((int)strlen(&___3642[___1926]), ___684); memmove(___3946, &___3642[___1926], ___2224); ___3946[___2224] = '\0'; ENSURE(VALID_REF(___3946) && (int)strlen(___3946) == ___2224); } char *___3818(char *___3813) { char *___3359 = ___3813; char *Start = ___3813; REQUIRE(VALID_REF(___3813)); while (tecplot::isspace(*Start)) Start++; if (Start != ___3813) memmove(___3813, Start, strlen(Start) + 1); ENSURE(VALID_REF(___3359) && ___3359 == ___3813); return ___3359; } static char* StringFlushRight(char *___3813) { REQUIRE(VALID_REF(___3813)); char* ___3359 = ___3813; char* End = ___3813 + strlen(___3813) - 1; while (End > ___3813 && tecplot::isspace(*End)) End--; *(End + 1) = '\0'; ENSURE(VALID_REF(___3359) && ___3359 == ___3813); return ___3359; } char *___4225(char *___3813) { REQUIRE((___3813 == NULL) || VALID_REF(___3813)); if (___3813) return (___3818(StringFlushRight(___3813))); else return ___3813; } char *___3859(char  *___3813, int    ___2376) { REQUIRE(VALID_REF(___3813)); REQUIRE(___2376 >= 0); if ((int)strlen(___3813) > ___2376) ___3813[___2376] = '\0'; ENSURE(VALID_REF(___3813)); ENSURE((int)strlen(___3813) <= ___2376); return ___3813; } char* ___3858( char* ___3813, int   ___2376) { REQUIRE(VALID_REF(___3813)); REQUIRE(___2376 >= 0); ___3818(___3813); ___3859(___3813,___2376); StringFlushRight(___3813); ENSURE(VALID_REF(___3813)); ENSURE((int)strlen(___3813) <= ___2376); ENSURE(IMPLICATION(strlen(___3813) != 0, (!tecplot::isspace(___3813[0]) && !tecplot::isspace(___3813[strlen(___3813)-1])))); return ___3813; }
 #ifndef MSWIN
___3839 ___2245(const char *___3813, uint32_t    ___4477) { REQUIRE(VALID_REF(___3813)); ___3839 ___3359 = ___3821(); if (___3359 != NULL) { ___372 ___2040 = ___4226; if (strlen(___3813) > ___4477) { char *StringCopy = ___1135(___1097(___3813)); ___2040 = (StringCopy != NULL); if (___2040) { char *___685 = StringCopy; char *SubString = StringCopy; uint32_t SubStringLen = 0; while (*___685 != '\0' && ___2040) { while (*___685 != '\0' && SubStringLen < ___4477) { if (*___685 == '\n') { *___685 = '\0'; ___685++; break; } ___685++; SubStringLen++; } if (*___685 != '\0' && SubStringLen == ___4477) { if (*___685 != ' ') { while (___685 != SubString && *___685 != ' ') ___685--; if (*___685 != ' ') { while (*___685 != '\0' && *___685 != ' ' && *___685 != '\n') ___685++; while (*___685 != '\0' && *___685 == ' ') ___685++; } } if (*___685 != '\0') { *___685 = '\0'; ___685++; } StringFlushRight(SubString); } ___2040 = ___3823(___3359, SubString); SubString = ___685; SubStringLen = 0; } ___1530(StringCopy, "StringCopy"); } } else ___2040 = ___3823(___3359, ___3813); if (!___2040) ___3828(&___3359); } ENSURE(___3359 == NULL || VALID_REF(___3359)); return ___3359; }
 #endif
int ustrncmp(const  char *___3430, const  char *___3431, size_t ___2223) { REQUIRE((___3430 == NULL) || VALID_REF(___3430)); REQUIRE((___3431 == NULL) || VALID_REF(___3431)); char *t1; char *t2; char ct1; char ct2; size_t ___1832 = 0; if ((___3430 == NULL) && (___3431 == NULL)) return 0; if (___3430 == NULL) return -1; else if (___3431 == NULL) return 1; t1 = (char*)___3430; t2 = (char*)___3431; while (*t1 && *t2 && (___1832 < ___2223)) { ct1 = ___444(*t1); ct2 = ___444(*t2); if (ct1 != ct2) return (ct1 - ct2); t1++; t2++; ___1832++; } if ((___1832 == ___2223) || ((*t1 == '\0') && (*t2 == '\0'))) return 0; else return ___444(*t1) - ___444(*t2); } int ustrcmp(const char *___3430, const char *___3431) { REQUIRE((___3430 == NULL) || VALID_REF(___3430)); REQUIRE((___3431 == NULL) || VALID_REF(___3431)); return (ustrncmp(___3430, ___3431, INT_MAX)); }
 #if !defined NO_ASSERTS
___372 ___1982(char       **___3433, const char  *___2700, ___372    ___2063, const char  *___1395, int          ___2262)
 #else
___372 ___1982(char       **___3433, const char  *___2700, ___372    ___2063)
 #endif
{ REQUIRE(VALID_REF(___3433)); REQUIRE(*___3433 == NULL || VALID_REF(*___3433)); REQUIRE(___2700 == NULL || VALID_REF(___2700)); REQUIRE(IMPLICATION(VALID_REF(*___3433), *___3433 != ___2700)); REQUIRE(VALID_BOOLEAN(___2063)); REQUIRE(VALID_NON_ZERO_LEN_STR(___1395)); REQUIRE(___2262 >= 1); if (*___3433) {
 #if !defined NO_ASSERTS && defined DEBUG_ALLOC
char S[80+1]; MakeDebugRecord(___1395, ___2262, "releasing", *___3433, S, 80); ___1530(*___3433, S);
 #else
___1530(*___3433, "");
 #endif
} if (___2700 == NULL) { *___3433 = NULL; return (___4226); } else {
 #if !defined NO_ASSERTS && defined DEBUG_ALLOC
char S[80+1]; MakeDebugRecord(___1395, ___2262, "duplicating", ___2700, S, 80); *___3433 = ___23(strlen(___2700) + 1, char, S);
 #else
 #  if defined MSWIN && defined _DEBUG && !defined(MAKEARCHIVE) && !defined(NO_ASSERTS)
 #     undef new
*___3433 = new(___1395, ___2262) char[strlen(___2700)+1];
 #     define new DEBUG_NEW
 #  else
*___3433 = ___23(strlen(___2700) + 1, char, "");
 #  endif
 #endif
if (*___3433) { strcpy(*___3433, ___2700); return (___4226); } else { if (___2063) ___1177(___4217("Out of memory")); return (___1305); } } } ___372 ___3941( char**      ___3433, char const* ___3856, ___372   ___958, ___372   convertNewlineToEscSeq) { size_t      CurLen; size_t      NewLen; int         NumNewlines = 0; char       *___2700; const char *___685 = ___3856; ___372   ___2040 = ___4226; REQUIRE(VALID_REF(___3433)); REQUIRE((___3856 == NULL) || VALID_REF(___3856)); REQUIRE(VALID_BOOLEAN(___958)); REQUIRE(VALID_BOOLEAN(convertNewlineToEscSeq)); if ((___3856  == NULL) || (*___3856 == '\0')) { if (___3856            && (*___3856 == '\0') && ___958) { char *TMP = (char *)___3856; ___1530(TMP, "empty string to add"); } } else { if (*___3433 == NULL) CurLen = 0; else CurLen = strlen(*___3433); while (*___685) if (*___685++ == '\n') NumNewlines++; NewLen = CurLen + strlen(___3856) + 1 + NumNewlines; ___2700 = ___23(NewLen, char, ___3856); if (___2700 == NULL) { if (___958) { char *TMP = (char *)___3856; ___1530(TMP, ___3856); } ___2040 = ___1305; } else { if (*___3433) { strcpy(___2700, *___3433); ___1530(*___3433, (CurLen > 0 ? *___3433 : "previous text")); } else *___2700 = '\0'; { char *NPtr = EndOfString(___2700); const char *APtr = ___3856; while (*APtr) { if ((*APtr == '\n') && convertNewlineToEscSeq) { *NPtr++ = '\\'; *NPtr++ = 'n'; } else *NPtr++ = *APtr; APtr++; } *NPtr = '\0'; } if (___958) { char *TMP = (char *)___3856; ___1530(TMP, ___3856); } *___3433 = ___2700; } } ENSURE(VALID_BOOLEAN(___2040)); return (___2040); } ___372 ___3940( char**      ___3433, char const* ___3856, ___372   convertNewlineToEscSeq) { size_t      CurLen; size_t      NewLen; int         NumNewlines = 0; char       *___2700; const char *___685 = ___3856; ___372   ___2040 = ___4226; REQUIRE(VALID_REF(___3433)); REQUIRE((___3856 == NULL) || VALID_REF(___3856)); REQUIRE(VALID_BOOLEAN(convertNewlineToEscSeq)); if ((___3856  != NULL) && (*___3856 != '\0')) { if (*___3433 == NULL) CurLen = 0; else CurLen = strlen(*___3433); while (*___685) if (*___685++ == '\n') NumNewlines++; NewLen = CurLen + strlen(___3856) + 1 + NumNewlines; ___2700 = ___23(NewLen, char, ___3856); if (___2700 == NULL) { ___2040 = ___1305; } else { if (*___3433) { strcpy(___2700, *___3433); ___1530(*___3433, (CurLen > 0 ? *___3433 : "previous text")); } else *___2700 = '\0'; { char *NPtr = EndOfString(___2700); const char *APtr = ___3856; while (*APtr) { if ((*APtr == '\n') && convertNewlineToEscSeq) { *NPtr++ = '\\'; *NPtr++ = 'n'; } else *NPtr++ = *APtr; APtr++; } *NPtr = '\0'; } *___3433 = ___2700; } } ENSURE(VALID_BOOLEAN(___2040)); return (___2040); } ___372 ___3939(char  **___3433, char    ___476) { REQUIRE(VALID_REF(___3433)); char S[2]; S[0] = ___476; S[1] = '\0'; return (___3941(___3433, S, ___1305, ___1305)); } ___372 ___3353(char **___3813) { size_t    ___1832; int       NewlineCount; size_t    ___2224; char     *Replacement; REQUIRE(VALID_REF(___3813)); REQUIRE(VALID_REF(*___3813)); NewlineCount = 0; ___2224 = strlen(*___3813); for (___1832 = 0; ___1832 < ___2224; ___1832++) if ((*___3813)[___1832] == '\n') NewlineCount++; if (NewlineCount != 0) { Replacement = ___23(___2224 + NewlineCount + 1, char, "replacement string"); if (Replacement != NULL) { size_t ___2106; for (___1832 = ___2106 = 0; ___1832 < ___2224 + 1; ___1832++, ___2106++) { if ((*___3813)[___1832] == '\n') { Replacement[___2106] = '\\'; ___2106++; Replacement[___2106] = 'n'; } else { Replacement[___2106] = (*___3813)[___1832]; } } ___478(___1832 == ___2224 + 1); ___478(___2106 == ___2224 + NewlineCount + 1); } ___1530(*___3813, "original string"); *___3813 = Replacement; } ENSURE(*___3813 == NULL || VALID_REF(*___3813)); return (*___3813 != NULL); }
