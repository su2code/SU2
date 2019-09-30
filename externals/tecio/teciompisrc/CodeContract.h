 #pragma once
#include <cstdlib>
#include <iostream>
 #ifndef VALID_REF
 #define VALID_REF(p)      ((p)  != 0)
 #define VALID_FN_REF(___3002)  ((___3002) != 0)
 #endif
 #if defined NO_ASSERTS
 #ifndef ASSERT
 #define ASSERT(___1245)
 #endif
 #ifndef ASSERT_ONLY
 #define ASSERT_ONLY(___1245)
 #endif
 #else   
 #ifndef ASSERT_ONLY
 #define ASSERT_ONLY(___1245) ___1245
 #endif
 #if defined NDEBUG    
 #ifndef ASSERT
template<class T> inline bool checkedAssert(T const& expr) { return (expr ? true : false); } inline bool checkedAssert(char const* expr) { return expr != 0; }
 #define ASSERT(___1245) \
 do \
 { \
 if (!checkedAssert(___1245)) \
 { \
 std::cerr << __FILE__ << ':' << __LINE__ << ':' << "Assertion '" << #___1245 << "' failed."; \
 abort(); \
 } \
 } while (0);
 #endif  
 #else   
 #ifndef ASSERT
#include <assert.h>
 #define ASSERT(___1245) assert(___1245)
 #endif
 #endif  
 #endif
 #if !defined ASSERT_ONLY_PARAM
 #define ASSERT_ONLY_PARAM(___2972) ASSERT_ONLY(___2972)
 #endif
 #if !defined INVARIANT
 #define INVARIANT(___1245) ASSERT(___1245)
 #endif
 #if !defined REQUIRE
 #define REQUIRE(___1245)   ASSERT(___1245)
 #endif
 #if !defined ENSURE
 #define ENSURE(___1245)    ASSERT(___1245)
 #endif
 #if !defined ___478
 #define ___478(___1245)     ASSERT(___1245)
 #endif
 #ifdef VERIFY
 #undef VERIFY
 #endif
 #ifndef VERIFY
 #if defined NO_ASSERTS
 #define VERIFY(___1245) ((void)(___1245))
 #elif defined NDEBUG
 #define VERIFY(___1245) \
 do \
 { \
 if ((___1245) == 0) \
 { \
 std::cerr << __FILE__ << ':' << __LINE__ << ':' << "Assertion '" << #___1245 << "' failed."; \
 abort(); \
 } \
 } while(0);
 #else
 #define VERIFY(___1245) assert(___1245)
 #endif
 #endif 
 #if !defined IMPLICATION
 #define IMPLICATION(___2894,___3258) (!(___2894) || (___3258))
 #endif
 #if !defined EQUIVALENCE
 #define EQUIVALENCE(___2894,___3258) ((___2894) == (___3258))
 #endif
 #define VALID_MAP_KEY(key,map) (map.find(key) != map.end())
 #if !defined VALID_REF_OR_NULL
 #define VALID_REF_OR_NULL(p)       (VALID_REF(p) || p == 0)
 #endif
 #if !defined VALID_BOOLEAN
 #define VALID_BOOLEAN(b)           ((b) == 1 || (b) == 0)
 #endif
 #if !defined VALID_ENUM
 #define VALID_ENUM(___4314, type)    (0 <= (___4314) && (___4314) < END_##type)
 #endif
 #if !defined VALID_CLASS_ENUM
 #define VALID_CLASS_ENUM(e) (static_cast<std::decay<decltype((e))>::type>(0) <= (e) && (e) < std::decay<decltype((e))>::type::END_ENUM)
 #endif
