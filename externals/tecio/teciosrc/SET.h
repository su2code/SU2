 #if defined EXTERN
 #undef EXTERN
 #endif
 #if defined ___3500
 #define EXTERN
 #else
 #define EXTERN extern
 #endif
 #ifndef _SET_H_INCLUDED
 #define _SET_H_INCLUDED
#include "ThirdPartyHeadersBegin.h"
#include <algorithm>
#include <set>
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "RawArray.h"
#include "CodeContract.h"
 #define ___2896(X,Y)           ((int)(((X)-1)/(Y)+1)*(Y))
 #define ___3479            (8*sizeof(___3483))
 #define ___3498            (((___3480)1)<<(___3479-1))
 #if defined _DEBUG
 #  define USE_FUNCTIONS_FOR_SETS
 #endif
struct ___3502 { ___3493 size; ___3481 data; };
 #define ___2056(___3476) ((___3476)==NULL)
inline size_t ___3482(___3501 ___3476) { REQUIRE(VALID_REF(___3476)); return ___3476->size / ___3479 * sizeof(___3483); } EXTERN ___3501 ___29(___372 ___3574); EXTERN void ___939(___3501 *___3476); EXTERN ___372 ___3496(void       *___2098, ___90  ___494); EXTERN ___372 ___1201(___3501     ___3476, ___3493 max_val, ___372  ___3574); EXTERN ___372 ___676(___3501    ___1121, ___3501    ___3656, ___372 ___3574); EXTERN ___372 ___83(___3501 ___1121, ___3501 ___3656, ___372 ___3574); EXTERN void ___493(___3501 ___3476);
 #if defined USE_FUNCTIONS_FOR_SETS
EXTERN ___372 ___17(___3501     ___3476, ___3493 ___2401, ___372  ___3574);
 #else
 #if defined __cplusplus
inline ___372 ___17(___3501     ___3476, ___3493 ___2401, ___372  ___3574) { if (___3476 && (___2401 + 1 <= ___3476->size || ___1201(___3476, ___2401 + 1, ___3574))) { ___3493 word = ___2401 / ___3479; ___3480 bit = (___3480)1 << (___2401 % ___3479); ___3476->data[word] |= bit; return ___4226; } else return ___1305; }
 #elif defined TECPLOTKERNEL
 #define ___17(___3476,___2401,___3574) \
 (((___3476) && \
 ((___2401)+1 <= (___3476)->size || \
 ___1201((___3476), (___2401)+1, (___3574)))) \
 ? (((___3476)->data[(___2401) / ___3479].___1346((___3480)1 << ((___2401) % ___3479))), ___4226) \
 : ___1305)
 #else
 #define ___17(___3476,___2401,___3574) \
 (((___3476) && \
 ((___2401)+1 <= (___3476)->size || \
 ___1201((___3476), (___2401)+1, (___3574)))) \
 ? (((___3476)->data[(___2401) / ___3479] |= (___3480)1 << ((___2401) % ___3479)), ___4226) \
 : ___1305)
 #endif
 #endif
EXTERN void ___3334(___3501     ___3476, ___3493 ___2401); EXTERN void ___957(___3501     ___3476, ___3493 ___2402); EXTERN ___372 ___1955(___3501     ___3476, ___3493 ___2402, ___372  ___3571);
 #if defined USE_FUNCTIONS_FOR_SETS
EXTERN ___372 ___1956(___3501     ___3476, ___3493 ___2401);
 #else
 #if defined __cplusplus
inline ___372 ___1956(___3501     ___3476, ___3493 ___2401) { if (___3476 && (0 <= ___2401 && ___2401 < ___3476->size)) { ___3493 word = ___2401 / ___3479; ___3480 bit = (___3480)1 << (___2401 % ___3479); return (___3476->data[word]&bit) != 0; } else return ___1305; }
 #elif defined TECPLOTKERNEL
 #define ___1956(___3476,___2401)  ((___3476 && (0<=(___2401) && (___2401)<(___3476)->size)) \
 ? ((___3476)->data[(___2401)/___3479].load()&((___3480)1<<((___2401)%___3479)))!=0 \
 : ___1305)
 #else
 #define ___1956(___3476,___2401)  ((___3476 && (0<=(___2401) && (___2401)<(___3476)->size)) \
 ? ((___3476)->data[(___2401)/___3479]&((___3480)1<<((___2401)%___3479)))!=0 \
 : ___1305)
 #endif
 #endif
EXTERN ___372 ___2015(___3501 ___3476); EXTERN ___372 ___1822(___3501 ___3476); EXTERN ___3493 ___2403(___3501 ___3476); EXTERN ___372 ___2033(___3501 ___3476); EXTERN ___3493 ___1761(___3501     ___3476, ___3493 ___3682); EXTERN ___3493 ___1769(___3501     ___3476, ___3493 ___3682); EXTERN ___372 ___1175(___3501  ___3477, ___3501  ___3478); ___3501 intersection( ___3501 ___3477, ___3501 ___3478); EXTERN ___372 ___2062(___3501 ___486, ___3501 ___2973); EXTERN ___3493 ___2404(___3501 ___3476, ___3493    ___2402); EXTERN ___3493 ___2867(___3501 ___3476, ___3493    ___2866); EXTERN ___372 ___677(___3501     ___1126, ___3493 ___1125, ___3501     ___3663, ___3493 ___3662); EXTERN void ___3560(___3501     ___3476, ___3493 ___3558, ___3493 ___3559, ___3493 ___3556);
 #define ___1746(___3476) (___1761((___3476), ___333))
 #define ___1751(___3476)  (___1769((___3476), ___333))
 #define ___1472(___2402, ___3476) \
 for (___2402 = ___1746((___3476)); \
 ___2402 != ___333; \
 ___2402 = ___1761((___3476), (___2402)))
 #define ForAllMembersInEntIndexSet(___2402, ___3476) \
 for (___2402 = static_cast<___1172>(___1746((___3476))); \
 ___2402 != static_cast<___1172>(___333); \
 ___2402 = static_cast<___1172>(___1761((___3476), (___2402))))
 #define ___1471(___2402, ___3476) \
 for (___2402 = ___1751((___3476)); \
 ___2402 != ___333; \
 ___2402 = ___1769((___3476), (___2402)))
namespace tecplot { template <typename T> std::vector<T> ___4194(___3501 ___2100) { REQUIRE(VALID_REF(___2100) || ___2100 == 0); std::vector<T> ___3358; size_t const count = ___2403(___2100); if (count != 0) { ___3358.reserve(count); ___3493 ___2085; ___1472(___2085,___2100) ___3358.push_back(static_cast<T>(___2085)); } return ___3358; } template <typename T> inline std::set<T> ___4186(___3501 const set) { REQUIRE(VALID_REF_OR_NULL(set)); ___1172 ___4314; std::set<T> ___3358; if (set != NULL) { ___1472(___4314, set) { ___3358.insert(static_cast<T>(___4314)); } } return ___3358; } template <typename CONTAINER> ___3501 ___4186( CONTAINER const& ___2099, bool             isSorted = true) { REQUIRE(IMPLICATION(isSorted && !___2099.empty(), ___2099[___2099.size()-1] == ___333 || ___2099[___2099.size()-1] == *std::max_element(&___2099[0], &___2099[0]+___2099.size()))); ___3501 ___3358 = ___29(___1305); if (___3358 == NULL) throw std::bad_alloc(); if (!___2099.empty()) { typename CONTAINER::value_type largestMember = static_cast<typename CONTAINER::value_type>(___333); if (isSorted) { for (typename CONTAINER::value_type const* iter = &___2099[___2099.size()-1]; iter >= &___2099[0]; --iter) if ((largestMember = *iter) != static_cast<typename CONTAINER::value_type>(___333)) break; } else { largestMember = *std::max_element(&___2099[0], &___2099[0]+___2099.size()); } if (largestMember != static_cast<typename CONTAINER::value_type>(___333)) { if (!___1201(___3358, static_cast<___3493>(largestMember + 1), ___1305)) { ___939(&___3358); throw std::bad_alloc(); } typename CONTAINER::value_type const* itemsArray = &___2099[0]; size_t const ___2812 = ___2099.size(); for (size_t ___1992 = 0; ___1992 < ___2812; ++___1992) if (itemsArray[___1992] != static_cast<typename CONTAINER::value_type>(___333)) (void)___17(___3358,static_cast<___3493>(itemsArray[___1992]),___1305); } } ENSURE(VALID_REF(___3358)); return ___3358; } template <typename T> void ___4185( ___3501       ___2100, ___3269<T>& ___3358) { REQUIRE(VALID_REF(___2100) || ___2100 == 0); size_t const count = ___2403(___2100); if (count != 0) { ___3358.reserve(count); ___3358.___3503(count); T* ___3360 = &___3358[0]; size_t ___2865 = 0; ___3493 ___2085; ___1472(___2085,___2100) ___3360[___2865++] = static_cast<T>(___2085); } else { ___3358.___3503(0); } } }
 #endif 
