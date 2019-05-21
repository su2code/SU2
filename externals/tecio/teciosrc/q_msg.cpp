#include "stdafx.h"
#include "MASTER.h"
 #define ___3259
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
 #if defined MSWIN
#include "UnicodeStringUtils.h"
 #endif
#include "ALLOC.h"
#include "STRUTIL.h"
#include "stringformat.h"
using std::string; using namespace tecplot;
 #define MAXCHARACTERSPERLINE 60
___372 ___4478(const char  *___2873, char       **___2700) { size_t ___2165; if (___2873 == NULL) return (___1305); ___2165 = strlen(___2873); *___2700 = ___23(___2165 + 1, char, "new error message string"); if (*___2700 == NULL) return (___1305); strcpy(*___2700, ___2873); if (___2165 > MAXCHARACTERSPERLINE) { char *LineStart = *___2700; char *LastWord  = *___2700; char *WPtr      = *___2700; while (WPtr && (*WPtr != '\0')) { size_t CurLen; WPtr = strchr(LineStart, '\n'); if (WPtr && ((WPtr - LineStart) < MAXCHARACTERSPERLINE)) { WPtr++; while (*WPtr == '\n') WPtr++; LineStart = WPtr; while (*WPtr == ' ') WPtr++; LastWord  = WPtr; continue; } WPtr = strchr(LastWord, ' '); if (WPtr != NULL) { while (*WPtr == ' ') WPtr++; } if (WPtr == NULL) { CurLen = strlen(LineStart); } else { CurLen = WPtr - LineStart; } if (CurLen > MAXCHARACTERSPERLINE) { if (LastWord == LineStart) { if (WPtr && (*WPtr != '\0')) { *(WPtr - 1) = '\n'; LastWord = WPtr; } } else { *(LastWord - 1) = '\n'; } LineStart = LastWord; } else LastWord = WPtr; } } return (___4226); } namespace { void SendMessageToFile( FILE*            F, char const*      S, MessageBoxType_e messageBoxType) { REQUIRE(VALID_REF(F)); REQUIRE(VALID_REF(S)); REQUIRE(VALID_ENUM(messageBoxType, MessageBoxType_e)); REQUIRE(messageBoxType == ___2443 || messageBoxType == ___2447 || messageBoxType == ___2444); char *S2; if (___4478(S, &S2)) { switch (messageBoxType) { case ___2443:       fprintf(F, "Error: %s\n",   S2); break; case ___2447:     fprintf(F, "Warning: %s\n", S2); break; case ___2444: fprintf(F, "Info: %s\n",    S2); break; default: ___478(___1305); break; } } ___1530(S2, "temp info string"); } } namespace { void defaultErrMsg(char const* Msg) { REQUIRE(VALID_REF(Msg));
 #if defined MSWIN
::MessageBoxA(NULL, Msg, "Error", MB_OK | MB_ICONERROR);
 #else
SendMessageToFile(stderr, Msg, ___2443);
 #endif
} } namespace { void postMessage( std::string const& ___2432, MessageBoxType_e   messageBoxType) { { if (messageBoxType == ___2443) defaultErrMsg(___2432.c_str()); else SendMessageToFile(stderr, ___2432.c_str(), messageBoxType); } } } namespace { void Message(const char* ___2432, MessageBoxType_e messageBoxType) { REQUIRE(VALID_NON_ZERO_LEN_STR(___2432)); REQUIRE(VALID_ENUM(messageBoxType, MessageBoxType_e)); REQUIRE(messageBoxType == ___2443 || messageBoxType == ___2447 || messageBoxType == ___2444); static ___372 InMessage = ___1305; if (!InMessage) { { SendMessageToFile(stderr, ___2432, messageBoxType); } } } } void Information(___4218 format, ...) { REQUIRE(!format.___2035()); ___372 cleanUp = ___4226; va_list  ___93; va_start(___93, format); char* ___2432 = ___4410(format.c_str(), ___93); va_end(___93); if (___2432 == NULL) { cleanUp = ___1305; ___2432 = (char*)format.c_str(); } Message(___2432, ___2444); if (cleanUp) ___1530(___2432, "message"); } void ___4447(___4218 format, ...) { REQUIRE(!format.___2035()); ___372 cleanUp = ___4226; va_list  ___93; va_start(___93, format); char* ___2432 = ___4410(format.c_str(), ___93); va_end(___93); if (___2432 == NULL) { cleanUp = ___1305; ___2432 = (char*)format.c_str(); } Message(___2432, ___2447); if (cleanUp) ___1530(___2432, "message"); } namespace { void PostErrorMessage(___4218 format, va_list          ___93) { REQUIRE(!format.___2035()); ___372 cleanUp = ___4226; char* ___2455 = ___4410(format.c_str(), ___93); if (___2455 == NULL) { cleanUp = ___1305; ___2455 = (char*)format.c_str(); } postMessage(___2455, ___2443); if (cleanUp) ___1530(___2455, "messageString"); } } void ___4407(___4218 format, va_list          ___93) { REQUIRE(!format.___2035()); static ___372 InErrMsg = ___1305; if (!InErrMsg) { PostErrorMessage(format, ___93); } } void ___1177(___4218 format, ...) { REQUIRE(!format.___2035()); va_list ___93; va_start(___93, format); PostErrorMessage(format, ___93); va_end(___93); }
