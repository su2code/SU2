#include "stdafx.h"
#include "MASTER.h"
#include "GLOBAL.h"
#include "TASSERT.h"
#include "UnicodeStringUtils.h"
#include "ThirdPartyHeadersBegin.h"
#include <string>
#include "ThirdPartyHeadersEnd.h"
using std::string; using std::wstring; namespace tecplot { static inline string* createUtf8StringTranslation(string const& str) { string *___3358 = new string(str); ENSURE(VALID_REF(___3358)); return ___3358; }
 #if defined MSWIN
static inline wstring* createWideStringTranslation(string const& str) { wstring *___3358 = new wstring; *___3358 = utf8ToWideString(str.c_str()); ENSURE(VALID_REF(___3358)); return ___3358; }
 #endif
 #if defined MSWIN
static inline wstring* createWideString(___4218::___2505 ___2504, string const&          str) { REQUIRE(___2504 == ___4218::___1104 || ___2504 == ___4218::___1098); wstring* ___3358; if (___2504 == ___4218::___1104) ___3358 = createWideStringTranslation(str); else ___3358 = new wstring(utf8ToWideString(str.c_str())); return ___3358; }
 #endif
void ___4218::___1932(___4218::___2505 ___2504, char const*            str, char const*            ASSERT_ONLY(___4219)) { REQUIRE(___2504 == ___1104 || ___2504 == ___1098); REQUIRE(VALID_REF_OR_NULL(str)); REQUIRE(VALID_REF_OR_NULL(___4219)); ___2493   = ___2504; ___2487 = (str == NULL); if (!___2487) ___2623 = str; ___2666 = NULL;
 #if defined MSWIN
___2675 = NULL;
 #endif
} ___4218::___4218() { ___1932(___1098, NULL, NULL); ENSURE(this->___2067()); } ___4218 ___4218::___2765() { return ___1097(NULL); } ___4218::___4218(___4218::___2505 ___2504, char const*            str, char const*            ___4219) { REQUIRE(___2504 == ___1104 || ___2504 == ___1098); REQUIRE(VALID_REF_OR_NULL(str)); REQUIRE(VALID_REF_OR_NULL(___4219)); ___1932(___2504, str, ___4219); ENSURE(this->___2067()); } ___4218::~___4218() { delete ___2666;
 #if defined MSWIN
delete ___2675;
 #endif
}
 #if !defined NO_ASSERTS
bool ___4218::___2067() const { ___478(IMPLICATION(___2487, ___2623.length() == 0));
 #if 0
___478(IMPLICATION(___2487, ___2493 == ___1098));
 #endif
return true; }
 #endif
bool ___4218::___2035() const { INVARIANT(this->___2067()); return ___2487; } bool ___4218::___2036() const { INVARIANT(this->___2067()); return ___2487 || ___2623.length() == 0; } char const* ___4218::c_str() const { INVARIANT(this->___2067()); char const* ___3358 = NULL; if (!___2035()) { if (___2493 == ___1104) { if (___2666 == NULL) ___2666 = createUtf8StringTranslation(___2623); ___3358 = ___2666->c_str(); } else ___3358 = ___2623.c_str(); } ENSURE(___3358 == NULL || VALID_REF(___3358)); return ___3358; }
 #if defined MSWIN && !defined MAKEARCHIVE
wchar_t const* ___4218::___798() const { INVARIANT(this->___2067()); wchar_t const* ___3358 = NULL; if (!___2035()) { if (___2675 == NULL) ___2675 = createWideString(___2493, ___2623); ___3358 = ___2675->c_str(); } ENSURE(___3358 == NULL || VALID_REF(___3358)); return ___3358; }
 #endif
___4218::operator string() { INVARIANT(this->___2067()); REQUIRE(!___2035()); string* ___3358; if (___2493 == ___1104) { if (___2666 == NULL) ___2666 = createUtf8StringTranslation(___2623); ___3358 = ___2666; } else ___3358 = &___2623; return *___3358; }
 #if defined MSWIN && !defined MAKEARCHIVE
___4218::operator wstring() { INVARIANT(this->___2067()); REQUIRE(!___2035()); if (___2675 == NULL) ___2675 = createWideString(___2493, ___2623); return *___2675; }
 #endif
___4218& ___4218::operator=(___4218 const& ___2888) { REQUIRE(___2888.___2067()); if (this != &___2888) { ___2493       = ___2888.___2493; ___2487     = ___2888.___2487; ___2623     = ___2888.___2623; delete ___2666; ___2666 = (___2888.___2666 != NULL ? new string(*___2888.___2666) : NULL);
 #if defined MSWIN
delete ___2675; ___2675 = (___2888.___2675 != NULL ? new wstring(*___2888.___2675) : NULL);
 #endif
} ENSURE(this->___2067()); return *this; } ___4218::___4218(___4218 const& ___2888) { REQUIRE(___2888.___2067()); ___2493       = ___2888.___2493; ___2487     = ___2888.___2487; ___2623     = ___2888.___2623; ___2666 = (___2888.___2666 != NULL ? new string(*___2888.___2666) : NULL);
 #if defined MSWIN
___2675 = (___2888.___2675 != NULL ? new wstring(*___2888.___2675) : NULL);
 #endif
ENSURE(this->___2067()); } ___4218 ___4218::___4217(char const* str, char const* ___4219) { REQUIRE(VALID_REF_OR_NULL(str)); REQUIRE(VALID_REF_OR_NULL(___4219)); return ___4218(___1104, str, ___4219); } ___4218 ___4218::___1097(char const* str) { REQUIRE(VALID_REF_OR_NULL(str)); return ___4218(___1098, str, NULL); } }
