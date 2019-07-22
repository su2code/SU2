 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <stdexcept>
#include <typeinfo>
#include "ThirdPartyHeadersEnd.h"
#include "LightweightVector.h"
namespace tecplot { namespace ___3933 { class AnyTypeLightweightVector { private: struct VectorHolderBase { virtual ~VectorHolderBase(){}
 #if !defined NO_ASSERTS
virtual uint64_t size() const = 0;
 #endif
virtual void ___3504(uint64_t i, double ___4314) = 0; virtual double operator[] (uint64_t i) const = 0; }; template <typename T, bool isBitArray = false> struct VectorHolder : public VectorHolderBase { ___2240<T> lightweightVector; explicit VectorHolder(uint64_t size) { REQUIRE(size > 0); if (!lightweightVector.alloc(size)) throw std::bad_alloc(); } VectorHolder(uint64_t size, T const& initialValue) { REQUIRE(size > 0); if (!lightweightVector.alloc(size, initialValue)) throw std::bad_alloc(); } virtual ~VectorHolder() {}
 #if !defined NO_ASSERTS
virtual uint64_t size() const { return lightweightVector.size(); }
 #endif
virtual void ___3504(uint64_t i, double ___4314) { lightweightVector[i] = static_cast<T>(___4314); } virtual double operator[](uint64_t i) const { return static_cast<double>(lightweightVector[i]); } }; VectorHolderBase* m_vectorHolder; UNCOPYABLE_CLASS(AnyTypeLightweightVector); public: AnyTypeLightweightVector() : m_vectorHolder(NULL) {} ~AnyTypeLightweightVector() { if (m_vectorHolder) delete m_vectorHolder; }
 #if (defined _MSC_VER && __cplusplus >= 199711L) || __cplusplus >= 201103L
AnyTypeLightweightVector(AnyTypeLightweightVector&& that) { m_vectorHolder = that.m_vectorHolder; that.m_vectorHolder = nullptr; } AnyTypeLightweightVector& operator=(AnyTypeLightweightVector&& that) { if (this != &that) { ___937(); m_vectorHolder = that.m_vectorHolder; that.m_vectorHolder = nullptr; } return *this; }
 #else
void moveAssignment(AnyTypeLightweightVector& that) { if (this != &that) { ___937(); m_vectorHolder = that.m_vectorHolder; that.m_vectorHolder = NULL; } }
 #endif
template <typename T, bool isBitArray> ___372 alloc(uint64_t size) { REQUIRE(empty()); REQUIRE(size > 0); try { m_vectorHolder = new VectorHolder<T, isBitArray>(size); } catch (...) { return ___1305; } return ___4226; } ___372 alloc(FieldDataType_e ___1363, uint64_t size) { REQUIRE(VALID_ENUM(___1363, FieldDataType_e)); REQUIRE(size > 0); REQUIRE(empty()); switch (___1363) { case FieldDataType_Float: return alloc<float, false>(size); case FieldDataType_Double: return alloc<double, false>(size); case FieldDataType_Int32: return alloc<int32_t, false>(size); case FieldDataType_Int16: return alloc<int16_t, false>(size); case FieldDataType_Byte: return alloc<uint8_t, false>(size); case ___1365: return alloc<uint8_t, true>(size); default: ___478(___1305); return ___1305; } } ___372 alloc(FieldDataType_e ___1363, uint64_t size, double initialValue) { REQUIRE(VALID_ENUM(___1363, FieldDataType_e)); REQUIRE(size > 0); REQUIRE(empty()); switch (___1363) { case FieldDataType_Float: return alloc<float, false>(size, static_cast<float>(initialValue)); case FieldDataType_Double: return alloc<double, false>(size, initialValue); case FieldDataType_Int32: return alloc<int32_t, false>(size, static_cast<int32_t>(initialValue)); case FieldDataType_Int16: return alloc<int16_t, false>(size, static_cast<int16_t>(initialValue)); case FieldDataType_Byte: return alloc<uint8_t, false>(size, static_cast<uint8_t>(initialValue)); case ___1365: return alloc<uint8_t, true>(size, (initialValue >= 1.0 ? 1 : 0)); default: ___478(___1305); return ___1305; } } template <typename T, bool isBitArray> ___372 alloc(uint64_t size, T const& initialValue) { REQUIRE(size > 0); REQUIRE(empty()); try { m_vectorHolder = new VectorHolder<T, isBitArray>(size, initialValue); } catch (...) { return ___1305; } return ___4226; } void ___937() { delete m_vectorHolder; m_vectorHolder = 0; }
 #if !defined NO_ASSERTS
uint64_t size() const { return m_vectorHolder ? m_vectorHolder->size() : 0;  }
 #endif
bool empty() const { return !m_vectorHolder;  } template <typename T, bool isBitArray> ___2240<T>& get() const { REQUIRE(typeid(*m_vectorHolder) == typeid(VectorHolder<T, isBitArray>)); return (static_cast<VectorHolder<T>*>(m_vectorHolder))->lightweightVector; } void ___3504(uint64_t i, double ___4314) { REQUIRE(i < size()); m_vectorHolder->___3504(i, ___4314); } double operator[](uint64_t i) const { REQUIRE(i < size());  return m_vectorHolder->operator[](i); } }; template <> struct AnyTypeLightweightVector::VectorHolder<uint8_t, true> : public AnyTypeLightweightVector::VectorHolderBase { ___2240<uint8_t> lightweightVector; explicit VectorHolder(uint64_t size) { REQUIRE(size > 0); if (!lightweightVector.alloc(numBytesForNumBits(size))) throw std::bad_alloc(); } VectorHolder(uint64_t size, uint8_t const& initialValue) { REQUIRE(size > 0); if (!lightweightVector.alloc(numBytesForNumBits(size), initialValue)) throw std::bad_alloc(); } virtual ~VectorHolder() {}
 #if !defined NO_ASSERTS
virtual uint64_t size() const { return 8 * lightweightVector.size(); }
 #endif
virtual void ___3504(uint64_t i, double ___4314) { REQUIRE(i < size()); if (___4314 >= 1.0) lightweightVector[i / 8] |= (1 << (i % 8)); else lightweightVector[i / 8] &= ~(1 << (i % 8)); } virtual double operator[](uint64_t i) const { REQUIRE(i < size()); return static_cast<double>((lightweightVector[i / 8] >> (i % 8)) & (uint8_t)0x1); } }; }}
