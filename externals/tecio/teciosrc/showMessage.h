 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <cstdarg>
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
extern bool _showMessage(MessageBoxType_e messageBoxType, char const* ___2432); namespace tecplot { namespace ___3933 { inline bool messageBox(MessageBoxType_e messageBoxType, char const* format, va_list args) { REQUIRE(VALID_ENUM(messageBoxType, MessageBoxType_e)); REQUIRE(VALID_REF(format)); char ___2455[5000]; vsprintf(___2455, format, args); return _showMessage(messageBoxType, ___2455); } inline bool ___1186(char const* format, ...) { va_list args; va_start(args, format); (void)messageBox(___2443, format, args); va_end(args); return false; } inline bool ___3247(char const* format, ...) { va_list args; va_start(args, format); bool ___2039 = (messageBox(___2449, format, args) == ___4226); va_end(args); return ___2039; } inline bool ___1931(FILE* file, char const* format, ...) { va_list args; va_start(args, format); bool ___2039; if (file == NULL) ___2039 = messageBox(___2444, format, args); else ___2039 = (vfprintf(file, format, args) >= 0); va_end(args); return ___2039; } inline bool ___1931(char const* format, ...) { va_list args; va_start(args, format); bool ___2039 = messageBox(___2444, format, args); va_end(args); return ___2039; } }}
