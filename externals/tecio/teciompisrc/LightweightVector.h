 #pragma once
#include "ThirdPartyHeadersBegin.h"
#  include <algorithm>
#  include <new>
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
#include "ClassMacros.h"
#include "CodeContract.h"
#include "MinMax.h"
#include "ItemAddress.h"
inline uint64_t numBytesForNumBits(uint64_t numBits) { return (numBits + 7) / 8; } namespace { template <typename T> uint64_t const* getRealMemPtr(T* const& ___3251) { REQUIRE(VALID_REF(___3251)); uint64_t const* uint64Array = (uint64_t const*)___3251; return uint64Array-1; } } namespace { template <typename T> uint64_t* getRealMemPtr(T*& ___3251) { REQUIRE(VALID_REF(___3251)); uint64_t* uint64Array = (uint64_t*)___3251; return uint64Array-1; } } namespace tecplot { namespace ___3933 { template<typename T> class ___2240 { public:
 #if (defined _MSC_VER && __cplusplus >= 199711L) || __cplusplus >= 201103L
___2240(___2240 const&) = delete; ___2240& operator=(___2240 const&) = delete; ___2240(___2240&& ___2888) : m_ptr(std::move(___2888.m_ptr))
 #if !defined NO_ASSERTS
, m_size(std::move(___2888.m_size))
 #endif
{ } ___2240& operator=(___2240&& ___3392) { if (this != &___3392) { ___937(); m_ptr = std::move(___3392.m_ptr); ___3392.m_ptr = nullptr;
 #if !defined NO_ASSERTS
m_size = std::move(___3392.m_size); ___3392.m_size = 0;
 #endif
} return *this; }
 #else
private: ___2240(___2240 const&); ___2240& operator=(___2240 const&); public:
 #endif
void ___2319( char const* container, size_t numElements) {
 #ifdef LARGE_ARRAY_MEMORY_LOGGING
size_t const MEMTRACK_CUTOFF = size_t(1000)*size_t(1000); if (numElements * sizeof(T) >= MEMTRACK_CUTOFF) { FILE *file = fopen("memtrack.txt", "at"); if (file) { fprintf(file, "%s\t%" "I64u" "\t%" "I64u" "\t%s\n", container, numElements, sizeof(T), typeid(T).___2685()); fclose(file); } else throw std::bad_alloc(); }
 #else
___4278(container); ___4278(numElements);
 #endif
} private: T* m_ptr;
 #ifndef NO_ASSERTS
uint64_t m_size;
 #endif
inline bool allocRawData(T*& ___3251, uint64_t requestedSize) { REQUIRE(___3251 == NULL); REQUIRE(requestedSize>0); uint64_t const totalBytesRequired = sizeof(uint64_t) + requestedSize*uint64_t(sizeof(T)); uint64_t* mem = NULL; if ( sizeof(size_t) == 4 && totalBytesRequired > uint32_t(-1) ) mem = NULL; else { ___2319("LightweightVector", requestedSize); mem = (uint64_t *)malloc( static_cast<size_t>(totalBytesRequired) ); } bool ___2039 = (mem != NULL); if ( ___2039 ) { mem[0] = requestedSize; ___3251 = (T *)&mem[1]; ___2039 = true; } return ___2039; } inline void freeRawData(T*& ___3251) { REQUIRE(VALID_REF(___3251)); free(getRealMemPtr(___3251)); ___3251 = NULL; } public: ___2240() : m_ptr(NULL)
 #ifndef NO_ASSERTS
, m_size(0)
 #endif
{ } ~___2240() { ___937(); } inline void swap(___2240<T>& ___2888) { using std::swap; swap(m_ptr, ___2888.m_ptr);
 #ifndef NO_ASSERTS
swap(m_size, ___2888.m_size);
 #endif
} inline uint64_t size() const { uint64_t ___3358; if ( empty() ) ___3358 = 0; else ___3358 = getRealMemPtr(m_ptr)[0]; ENSURE(___3358 == m_size); return ___3358; } inline uint64_t numBytesAllocated(uint64_t knownSize) const { REQUIRE(IMPLICATION(empty(),knownSize==0)); REQUIRE(IMPLICATION(!empty(),knownSize==size())); if ( empty() ) return 0; else return sizeof(uint64_t) + knownSize*uint64_t(sizeof(T)); } inline bool empty() const { INVARIANT(EQUIVALENCE(m_ptr==NULL,m_size==0)); return ( m_ptr == NULL ); } inline bool alloc(uint64_t requestedSize) { REQUIRE(empty()); REQUIRE(requestedSize>0); bool ___2039; if ( !empty() ) ___2039 = false; else { ___2039 = allocRawData(m_ptr, requestedSize); if ( ___2039 ) {
 #ifndef NO_ASSERTS
m_size = requestedSize;
 #endif
uint64_t pos = 0; try { for ( pos = 0; pos < requestedSize; pos++ ) ::new(&m_ptr[pos]) T; } catch (...) { for ( uint64_t pos2 = 0; pos2 < pos; pos++ ) m_ptr[pos].~T(); freeRawData(m_ptr);
 #ifndef NO_ASSERTS
m_size = 0;
 #endif
___2039 = false; } } } return ___2039; } inline bool alloc(uint64_t requestedSize, T const& padVal) { REQUIRE(empty()); REQUIRE(requestedSize>0); bool ___2039; if ( !empty() ) ___2039 = false; else { ___2039 = allocRawData(m_ptr, requestedSize); if ( ___2039 ) {
 #ifndef NO_ASSERTS
m_size = requestedSize;
 #endif
uint64_t pos = 0; try { for ( pos = 0; pos < requestedSize; pos++ ) ::new(&m_ptr[pos]) T(padVal); } catch (...) { for ( uint64_t pos2 = 0; pos2 < pos; pos++ ) m_ptr[pos].~T(); freeRawData(m_ptr);
 #ifndef NO_ASSERTS
m_size = 0;
 #endif
___2039 = false; } } } return ___2039; } inline bool reallocate(uint64_t requestedSize) { bool ___2039 = true; if (empty()) { ___2039 = alloc(requestedSize); } else if (size() != requestedSize) { uint64_t const origSize = size(); T* newPtr = 0; ___2039 = allocRawData(newPtr, requestedSize); if ( ___2039 ) { uint64_t pos = 0; try { for (pos = 0; pos < requestedSize; ++pos) ::new(&newPtr[pos]) T; uint64_t const numItemsToSwap = std::min(origSize,requestedSize); for (pos = 0; pos < numItemsToSwap; ++pos) newPtr[pos].swap(m_ptr[pos]); for (pos = 0; pos < origSize; ++pos) m_ptr[pos].~T(); freeRawData(m_ptr); m_ptr = newPtr;
 #ifndef NO_ASSERTS
m_size = requestedSize;
 #endif
} catch (...) { for (uint64_t pos2 = 0; pos2 < pos; ++pos) newPtr[pos].~T(); freeRawData(newPtr); ___2039 = false; } } } return ___2039; } inline void ___937() { if ( m_ptr ) { uint64_t const ___2812 = size(); for ( uint64_t pos = 0; pos < ___2812; pos++ ) m_ptr[pos].~T(); freeRawData(m_ptr);
 #ifndef NO_ASSERTS
m_size = 0;
 #endif
} ENSURE(empty()); } inline T& operator[](uint64_t pos) { REQUIRE(pos<size()); return m_ptr[pos]; } inline T const& operator[](uint64_t pos) const { REQUIRE(pos<size()); return m_ptr[pos]; } inline T* data() { return m_ptr; } inline T const* data() const { return m_ptr; } typedef T* iterator; inline iterator begin() { return m_ptr; } inline iterator end(uint64_t knownSize) { REQUIRE(size() == knownSize); return m_ptr+knownSize; } typedef T const* const_iterator; inline const_iterator begin() const { return m_ptr; } inline const_iterator end(uint64_t knownSize) const { REQUIRE(size() == knownSize); return m_ptr+knownSize; } };
 #define LWV_SPECIALIZE_PLAIN_DATA_VECTORS
 #ifdef LWV_SPECIALIZE_PLAIN_DATA_VECTORS
template<typename T> class ___3094 { public:
 #if (defined _MSC_VER && __cplusplus >= 199711L) || __cplusplus >= 201103L
___3094(___3094 const&) = delete; ___3094& operator=(___3094 const&) = delete;
 #if defined _MSC_VER && _MSC_VER <= 1800 
 #if !defined NO_ASSERTS
___3094(___3094&& ___2888) : m_ptr(std::move(___2888.m_ptr)) , m_size(std::move(___2888.m_size)) {} ___3094& operator=(___3094&& ___3392) { if (this != &___3392) { m_ptr = std::move(___3392.m_ptr); m_size = std::move(___3392.m_size); } return *this; }
 #else
___3094(___3094&& ___2888) : m_ptr(std::move(___2888.m_ptr)) {} ___3094& operator=(___3094&& ___3392) { if (this != &___3392) { m_ptr = std::move(___3392.m_ptr); } return *this; }
 #endif
 #else
___3094(___3094&&) = default; ___3094& operator=(___3094&&) = default;
 #endif
 #else
private: ___3094(___3094 const&); ___3094& operator=(___3094 const&); public:
 #endif
void ___2319( char const* container, size_t numElements) {
 #ifdef LARGE_ARRAY_MEMORY_LOGGING
size_t const MEMTRACK_CUTOFF = size_t(1000)*size_t(1000); if (numElements * sizeof(T) >= MEMTRACK_CUTOFF) { FILE *file = fopen("memtrack.txt", "at"); if (file) { fprintf(file, "%s\t%" "I64u" "\t%" "I64u" "\t%s\n", container, numElements, sizeof(T), typeid(T).___2685()); fclose(file); } else throw std::bad_alloc(); }
 #else
___4278(container); ___4278(numElements);
 #endif
} private: T* m_ptr;
 #ifndef NO_ASSERTS
uint64_t m_size;
 #endif
public: ___3094() : m_ptr(NULL)
 #ifndef NO_ASSERTS
, m_size(0)
 #endif
{ } ~___3094() { ___937(); } inline void swap(___3094<T>& ___2888) { using std::swap; swap(m_ptr, ___2888.m_ptr);
 #ifndef NO_ASSERTS
swap(m_size, ___2888.m_size);
 #endif
} inline uint64_t numBytesAllocated(uint64_t knownSize) const { REQUIRE(knownSize==m_size); REQUIRE(EQUIVALENCE(empty(),knownSize==0)); if ( empty() ) return 0; else return knownSize*sizeof(T); } inline bool empty() const { REQUIRE(EQUIVALENCE(m_ptr==NULL,m_size==0)); return ( m_ptr == NULL ); }
 #ifndef NO_ASSERTS
inline uint64_t size() const { return m_size; }
 #endif
inline bool alloc(uint64_t requestedSize) { REQUIRE(empty()); REQUIRE(requestedSize>0); bool ___2039; if ( !empty() || requestedSize == 0 ) ___2039 = false; else { uint64_t const totalBytesRequired = requestedSize*uint64_t(sizeof(T)); if ( sizeof(size_t) == 4 && totalBytesRequired > uint32_t(-1) ) m_ptr = NULL; else { ___2319("PlainDataVector", requestedSize); m_ptr = (T *)malloc( static_cast<size_t>(totalBytesRequired) ); } ___2039 = (m_ptr != NULL);
 #ifndef NO_ASSERTS
if ( ___2039 ) m_size = requestedSize;
 #endif
} ENSURE(EQUIVALENCE(___2039, VALID_REF(m_ptr))); return ___2039; } inline bool alloc(uint64_t requestedSize, T padVal) { REQUIRE(empty()); bool ___2039 = alloc(requestedSize); if ( ___2039 ) { try { for ( uint64_t pos = 0; pos < requestedSize; pos++ ) new(&m_ptr[pos]) T(padVal); } catch (...) { ___937(); ___2039 = false; } } return ___2039; } inline bool reallocate( uint64_t origSize, uint64_t requestedSize) { REQUIRE(origSize == size()); bool ___2039 = true; if (empty()) { ___2039 = alloc(requestedSize); } else if (origSize != requestedSize) { T* newPtr = 0; uint64_t const totalBytesRequired = requestedSize*uint64_t(sizeof(T)); if ( sizeof(size_t) == 4 && totalBytesRequired > uint32_t(-1) ) newPtr = NULL; else { ___2319("PlainDataVector", requestedSize); newPtr = (T *)malloc( static_cast<size_t>(totalBytesRequired) ); } ___2039 = (newPtr != NULL); if ( ___2039 ) { uint64_t const bytesToCopy = std::min(origSize,requestedSize)*uint64_t(sizeof(T)); uint8_t const* const sourceBytePtr = reinterpret_cast<uint8_t const*>(m_ptr); uint8_t* const targetBytePtr = reinterpret_cast<uint8_t*>(newPtr); std::copy(sourceBytePtr, sourceBytePtr+bytesToCopy, targetBytePtr); free(m_ptr); m_ptr = newPtr;
 #ifndef NO_ASSERTS
m_size = requestedSize;
 #endif
} } return ___2039; } inline void ___937() { if ( m_ptr ) { free(m_ptr); m_ptr = NULL;
 #ifndef NO_ASSERTS
m_size = 0;
 #endif
} ENSURE(empty()); } inline T& operator[](uint64_t pos) { REQUIRE(!empty() && pos < size()); return m_ptr[pos]; } inline T const& operator[](uint64_t pos) const { REQUIRE(!empty() && pos < size()); return m_ptr[pos]; } inline T* data() { return m_ptr; } inline T const* data() const { return m_ptr; } typedef T* iterator; inline iterator begin() { return m_ptr; } inline iterator end(uint64_t knownSize) { REQUIRE(knownSize == size()); return m_ptr+knownSize; } typedef T const* const_iterator; inline const_iterator begin() const { return m_ptr; } inline const_iterator end(uint64_t knownSize) const { REQUIRE(knownSize == size()); return m_ptr+knownSize; } }; template<> class ___2240<___2090> : public ___3094<___2090> {}; template<> class ___2240<___2090::SubzoneAddress> : public ___3094<___2090::SubzoneAddress> {}; template<> class ___2240<double> : public ___3094<double> {}; template<> class ___2240<float> : public ___3094<float> {}; template<> class ___2240<uint64_t> : public ___3094<uint64_t> {}; template<> class ___2240<int64_t> : public ___3094<int64_t> {}; template<> class ___2240<uint32_t> : public ___3094<uint32_t> {}; template<> class ___2240<int32_t> : public ___3094<int32_t> {}; template<> class ___2240<uint16_t> : public ___3094<uint16_t> {}; template<> class ___2240<int16_t> : public ___3094<int16_t> {}; template<> class ___2240<uint8_t> : public ___3094<uint8_t> {}; template<> class ___2240<int8_t> : public ___3094<int8_t> {}; template<typename T> class ___2240<T*> : public ___3094<T*> {}; template<> class ___2240<___2479> : public ___3094<___2479> {};
 #endif
template<typename T> inline bool ___3356(___2240<___2240<T> >& twoDLwVector, uint64_t newDim1, uint64_t newDim2) { REQUIRE(newDim1>0 && newDim2>0); bool ___2039 = false; try { twoDLwVector.alloc(newDim1); for ( uint64_t ___1841 = 0; ___1841 < newDim1; ___1841++ ) twoDLwVector[___1841].alloc(newDim2); ___2039 = true; } catch (...) { if ( !twoDLwVector.empty() ) { ___478(twoDLwVector.size() == newDim1); for ( uint64_t ___1841 = 0; ___1841 < newDim1; ___1841++ ) twoDLwVector[___1841].___937(); twoDLwVector.___937(); } ___2039 = false; } ENSURE(IMPLICATION(___2039, twoDLwVector.size()==newDim1)); ENSURE(IMPLICATION(___2039, twoDLwVector[0].size()==newDim2)); ENSURE(IMPLICATION(___2039, twoDLwVector[newDim1/2].size()==newDim2)); ENSURE(IMPLICATION(___2039, twoDLwVector[newDim1-1].size()==newDim2)); return ___2039; } template<typename T> inline bool ___3356(___2240<___2240<T> >& twoDLwVector, uint64_t newDim1, uint64_t newDim2, T padValue) { REQUIRE(newDim1>0 && newDim2>0); bool ___2039 = false; try { twoDLwVector.alloc(newDim1); for ( uint64_t ___1841 = 0; ___1841 < newDim1; ___1841++ ) twoDLwVector[___1841].alloc(newDim2, padValue); ___2039 = true; } catch (...) { if ( !twoDLwVector.empty() ) { ___478(twoDLwVector.size() == newDim1); for ( uint64_t ___1841 = 0; ___1841 < newDim1; ___1841++ ) twoDLwVector[___1841].___937(); twoDLwVector.___937(); } ___2039 = false; } ENSURE(IMPLICATION(___2039, twoDLwVector.size()==newDim1)); ENSURE(IMPLICATION(___2039, twoDLwVector[0].size()==newDim2)); ENSURE(IMPLICATION(___2039, twoDLwVector[newDim1/2].size()==newDim2)); ENSURE(IMPLICATION(___2039, twoDLwVector[newDim1-1].size()==newDim2)); return ___2039; } }}
