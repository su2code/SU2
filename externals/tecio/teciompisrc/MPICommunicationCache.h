 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <boost/scoped_ptr.hpp>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include "ThirdPartyHeadersEnd.h"
#include "SimpleVector.h"
#include "StandardIntegralTypes.h"
namespace tecplot { namespace teciompi { class MPICommunicationCache { public: class Error : public std::logic_error { public: Error(std::string const& ___1186) : std::logic_error(___1186) {} }; MPICommunicationCache(); explicit MPICommunicationCache(SimpleVector<uint8_t> const& serializedData); MPICommunicationCache(uint8_t const* serializedData, int serializedDataSize); ~MPICommunicationCache(); template <typename T> void addScalar(T ___4298, int32_t tag); template <typename T> void retrieveScalar(T& ___4298, int32_t tag) const; template <typename T> void addVector(SimpleVector<T> const& vec, int32_t tag); template <typename T> void retrieveVector(SimpleVector<T>& vec, int32_t tag) const; SimpleVector<uint8_t> data() const; private: struct Impl; boost::scoped_ptr<Impl> const m_impl; }; }}
