 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <algorithm>
#include "ThirdPartyHeadersEnd.h"
#include "CodeContract.h"
class ___2479 : public std::pair<double, double> {
 #define INVALID_MINMAX_MIN_VALUE (10.*___2179) 
 #define INVALID_MINMAX_MAX_VALUE (-10.*___2179)
public: ___2479() { invalidate(); } ___2479(double newMin, double newMax) { first  = newMin; second = newMax; } explicit ___2479(double ___4298) { first  = ___4298; second = first; } inline void swap(___2479& ___2888) { using std::swap; swap(first, ___2888.first); swap(second, ___2888.second); } inline double minValue(void) const { return first; } inline double maxValue(void) const { return second; } inline void ___3499(double newMin, double newMax) { REQUIRE(newMin <= newMax); first  = newMin; second = newMax; } inline void ___3499(___2479 const& ___2888) { first  = ___2888.first; second = ___2888.second; } inline void include(double ___4298) { first  = std::min(first, ___4298); second = std::max(second, ___4298); } inline void include(___2479 const& minMax) { REQUIRE(minMax.___2067()); first  = std::min(first, minMax.first); second = std::max(second, minMax.second); } inline bool containsValue(double ___4298) const { return ( first <= ___4298 && ___4298 <= second ); } inline void invalidate(void) { first  = INVALID_MINMAX_MIN_VALUE; second = INVALID_MINMAX_MAX_VALUE; } inline bool ___2067(void) const { REQUIRE(*(uint64_t *)&first != 0xcdcdcdcdcdcd || *(uint64_t *)&second != 0xcdcdcdcdcdcd); bool const valid = (first <= second); ___478(IMPLICATION(!valid, first==INVALID_MINMAX_MIN_VALUE && second==INVALID_MINMAX_MAX_VALUE)); return valid; } };
