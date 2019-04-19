#include "MPINonBlockingCommunicationCollection.h"
#include "ThirdPartyHeadersBegin.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/ref.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
#include "CodeContract.h"
#include "MinMax.h"
#include "mpiDatatype.h"
#include "MPIError.h"
namespace tecplot { namespace teciompi { namespace { void throwMPIError(char const* method, char const* mpiRoutine, int errorCode) { std::ostringstream ___2432; ___2432 << "Error in " << method << ": " << mpiRoutine << " returned error code " << errorCode; throw(MPIError(___2432.str())); } } template <typename T> MPINonBlockingCommunicationCollection::MPINonBlockingSendScalar<T>::MPINonBlockingSendScalar( T const& ___4298, int dest, int tag, MPI_Comm comm)
 #if defined _DEBUG
: MPINonBlockingCommunication(dest, tag, "send")
 #else
: MPINonBlockingCommunication(dest, tag)
 #endif
, m_complete(false) { int ___3358 = MPI_Isend(const_cast<T*>(&___4298), 1, mpiDatatype<T>(), dest, tag, comm, &m_request); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingSendScalar", "MPI_Isend", ___3358); } template <> MPINonBlockingCommunicationCollection::MPINonBlockingSendScalar<___2479>::MPINonBlockingSendScalar( ___2479 const& ___4298, int dest, int tag, MPI_Comm comm)
 #if defined _DEBUG
: MPINonBlockingCommunication(dest, tag, "send")
 #else
: MPINonBlockingCommunication(dest, tag)
 #endif
, m_complete(false) { int ___3358 = MPI_Isend(const_cast<___2479*>(&___4298), 2, MPI_DOUBLE, dest, tag, comm, &m_request); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingSendScalar", "MPI_Isend", ___3358); } template <typename T> MPINonBlockingCommunicationCollection::MPINonBlockingSendScalar<T>::~MPINonBlockingSendScalar() { if (!isComplete() && m_request != MPI_REQUEST_NULL) MPI_Request_free(&m_request); } template <typename T> bool tecplot::teciompi::MPINonBlockingCommunicationCollection::MPINonBlockingSendScalar<T>::isComplete() { if (!m_complete) { int completeFlag; int ___3358 = MPI_Test(&m_request, &completeFlag, MPI_STATUS_IGNORE); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingSendScalar::isComplete", "MPI_Test", ___3358); m_complete = (completeFlag == 1); } return m_complete; } template <typename T> void tecplot::teciompi::MPINonBlockingCommunicationCollection::MPINonBlockingSendScalar<T>::___4446() { int ___3358 = MPI_Wait(&m_request, MPI_STATUS_IGNORE); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingSendScalar::waitUntilComplete", "MPI_Wait", ___3358); m_complete = true; } template <typename T> MPINonBlockingCommunicationCollection::MPINonBlockingSendScalarCopy<T>::MPINonBlockingSendScalarCopy( T ___4298, int dest, int tag, MPI_Comm comm)
 #if defined _DEBUG
: MPINonBlockingCommunication(dest, tag, "send")
 #else
: MPINonBlockingCommunication(dest, tag)
 #endif
, m_val(___4298) , m_sendScalar(m_val, dest, tag, comm) {} template <typename T> MPINonBlockingCommunicationCollection::MPINonBlockingSendScalarCopy<T>::~MPINonBlockingSendScalarCopy() { if (!isComplete()) {
 #if defined _DEBUG
std::cerr << "Warning: Uncompleted MPINonBlockingSendScalarCopy going out of scope." << std::endl;
 #endif
} } template <typename T> bool tecplot::teciompi::MPINonBlockingCommunicationCollection::MPINonBlockingSendScalarCopy<T>::isComplete() { return m_sendScalar.isComplete(); } template <typename T> void tecplot::teciompi::MPINonBlockingCommunicationCollection::MPINonBlockingSendScalarCopy<T>::___4446() { m_sendScalar.___4446(); } template <typename T> tecplot::teciompi::MPINonBlockingCommunicationCollection::MPINonBlockingReceiveScalar<T>:: MPINonBlockingReceiveScalar(T& ___4298, int ___3656, int tag, MPI_Comm comm)
 #if defined _DEBUG
: MPINonBlockingCommunication(___3656, tag, "recv")
 #else
: MPINonBlockingCommunication(___3656, tag)
 #endif
, m_val(___4298) , m_complete(false) , m_request(MPI_REQUEST_NULL) { int ___3358 = MPI_Irecv(&m_val, 1, mpiDatatype<T>(), ___3656, tag, comm, &m_request); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingReceiveScalar", "MPI_Irecv", ___3358); } template <> tecplot::teciompi::MPINonBlockingCommunicationCollection::MPINonBlockingReceiveScalar<___2479>:: MPINonBlockingReceiveScalar(___2479& ___4298, int ___3656, int tag, MPI_Comm comm)
 #if defined _DEBUG
: MPINonBlockingCommunication(___3656, tag, "recv")
 #else
: MPINonBlockingCommunication(___3656, tag)
 #endif
, m_val(___4298) , m_complete(false) , m_request(MPI_REQUEST_NULL) { int ___3358 = MPI_Irecv(&m_val, 2, MPI_DOUBLE, ___3656, tag, comm, &m_request); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingReceiveScalar", "MPI_Irecv", ___3358); } template <typename T> MPINonBlockingCommunicationCollection::MPINonBlockingReceiveScalar<T>::~MPINonBlockingReceiveScalar() { if (!isComplete()) { if (m_request != MPI_REQUEST_NULL) MPI_Request_free(&m_request);
 #if defined _DEBUG
std::cerr << "Warning: Uncompleted MPINonBlockingReceiveScalar going out of scope." << std::endl;
 #endif
} } template <typename T> bool tecplot::teciompi::MPINonBlockingCommunicationCollection::MPINonBlockingReceiveScalar<T>::isComplete() { if (!m_complete) { int completeFlag; int ___3358 = MPI_Test(&m_request, &completeFlag, MPI_STATUS_IGNORE); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingReceiveScalar::isComplete", "MPI_Test", ___3358); m_complete = (completeFlag == 1); } return m_complete; } template <typename T> void tecplot::teciompi::MPINonBlockingCommunicationCollection::MPINonBlockingReceiveScalar<T>::___4446() { int ___3358 = MPI_Wait(&m_request, MPI_STATUS_IGNORE); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingReceiveScalar::waitUntilComplete", "MPI_Wait", ___3358); m_complete = true; } template <typename T> MPINonBlockingCommunicationCollection::MPINonBlockingSendVector<T>::MPINonBlockingSendVector( SimpleVector<T> const& vec, int dest, int sizeTag, int vecTag, MPI_Comm comm)
 #if defined _DEBUG
: MPINonBlockingCommunication(dest, sizeTag, "send")
 #else
: MPINonBlockingCommunication(dest, sizeTag)
 #endif
, m_vec(vec) , m_complete(false) , m_sizeRequest(MPI_REQUEST_NULL) , m_vecRequest(MPI_REQUEST_NULL) { int ___3358 = MPI_Isend(const_cast<int*>(&m_vec.size()), 1, MPI_INT, dest, sizeTag, comm, &m_sizeRequest); if (___3358 == MPI_SUCCESS) ___3358 = MPI_Isend(const_cast<T*>(m_vec.begin()), m_vec.size(), mpiDatatype<T>(), dest, vecTag, comm, &m_vecRequest); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingSendVectorCopy", "MPI_Isend", ___3358); } template <> MPINonBlockingCommunicationCollection::MPINonBlockingSendVector<___2479>::MPINonBlockingSendVector( SimpleVector<___2479> const& vec, int dest, int sizeTag, int vecTag, MPI_Comm comm)
 #if defined _DEBUG
: MPINonBlockingCommunication(dest, sizeTag, "send")
 #else
: MPINonBlockingCommunication(dest, sizeTag)
 #endif
, m_vec(vec) , m_complete(false) , m_sizeRequest(MPI_REQUEST_NULL) , m_vecRequest(MPI_REQUEST_NULL) { int ___3358 = MPI_Isend(const_cast<int*>(&m_vec.size()), 1, MPI_INT, dest, sizeTag, comm, &m_sizeRequest); if (___3358 == MPI_SUCCESS) ___3358 = MPI_Isend(const_cast<___2479*>(m_vec.begin()), 2 * m_vec.size(), MPI_DOUBLE, dest, vecTag, comm, &m_vecRequest); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingSendVectorCopy", "MPI_Isend", ___3358); } template <typename T> tecplot::teciompi::MPINonBlockingCommunicationCollection::MPINonBlockingSendVector<T>::~MPINonBlockingSendVector() { if (!isComplete()) { if (m_sizeRequest != MPI_REQUEST_NULL) MPI_Request_free(&m_sizeRequest); if (m_vecRequest != MPI_REQUEST_NULL) MPI_Request_free(&m_vecRequest); } } template <typename T> bool MPINonBlockingCommunicationCollection::MPINonBlockingSendVector<T>::isComplete() { if (!m_complete) { int sizeFlag; int vecFlag = 0; int ___3358 = MPI_Test(&m_sizeRequest, &sizeFlag, MPI_STATUS_IGNORE); if (___3358 == MPI_SUCCESS) ___3358 = MPI_Test(&m_vecRequest, &vecFlag, MPI_STATUS_IGNORE); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingSendVectorCopy::isComplete", "MPI_Test", ___3358); m_complete = (sizeFlag == 1 && vecFlag == 1); } return m_complete; } template <typename T> void MPINonBlockingCommunicationCollection::MPINonBlockingSendVector<T>::___4446() { int ___3358 = MPI_Wait(&m_sizeRequest, MPI_STATUS_IGNORE); if (___3358 == MPI_SUCCESS) ___3358 = MPI_Wait(&m_vecRequest, MPI_STATUS_IGNORE); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingSendVectorCopy::waitUntilComplete", "MPI_Wait", ___3358); m_complete = true; } template <typename T> MPINonBlockingCommunicationCollection::MPINonBlockingSendVectorCopy<T>::MPINonBlockingSendVectorCopy( SimpleVector<T> const& vec, int dest, int sizeTag, int vecTag, MPI_Comm comm)
 #if defined _DEBUG
: MPINonBlockingCommunication(dest, sizeTag, "send")
 #else
: MPINonBlockingCommunication(dest, sizeTag)
 #endif
, m_vec(vec) , m_sendVector(m_vec, dest, sizeTag, vecTag, comm) {} template <typename T> MPINonBlockingCommunicationCollection::MPINonBlockingSendVectorCopy<T>::MPINonBlockingSendVectorCopy( std::vector<T> const& vec, int dest, int sizeTag, int vecTag, MPI_Comm comm)
 #if defined _DEBUG
: MPINonBlockingCommunication(dest, sizeTag, "send")
 #else
: MPINonBlockingCommunication(dest, sizeTag)
 #endif
, m_vec(vec.begin(), vec.end()) , m_sendVector(m_vec, dest, sizeTag, vecTag, comm) {} template <typename T> tecplot::teciompi::MPINonBlockingCommunicationCollection::MPINonBlockingSendVectorCopy<T>::~MPINonBlockingSendVectorCopy() { if (!isComplete()) {
 #if defined _DEBUG
std::cerr << "Warning: Uncompleted MPINonBlockingSendVectorCopy going out of scope." << std::endl;
 #endif
} } template <typename T> bool MPINonBlockingCommunicationCollection::MPINonBlockingSendVectorCopy<T>::isComplete() { return m_sendVector.isComplete(); } template <typename T> void MPINonBlockingCommunicationCollection::MPINonBlockingSendVectorCopy<T>::___4446() { m_sendVector.___4446(); } template <typename T> MPINonBlockingCommunicationCollection::MPINonBlockingReceiveVector<T>::MPINonBlockingReceiveVector( SimpleVector<T>& vec, int ___3656, int sizeTag, int vecTag, MPI_Comm comm)
 #if defined _DEBUG
: MPINonBlockingCommunication(___3656, sizeTag, "recv")
 #else
: MPINonBlockingCommunication(___3656, sizeTag)
 #endif
, m_size(0) , m_vec(vec) , m_src(___3656) , m_vecTag(vecTag) , m_comm(comm) , m_complete(false) , m_sizeRequest(MPI_REQUEST_NULL) , m_vecRequest(MPI_REQUEST_NULL) { int ___3358 = MPI_Irecv(&m_size, 1, MPI_INT, ___3656, sizeTag, comm, &m_sizeRequest); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingReceiveVector", "MPI_Irecv", ___3358); } template <typename T> tecplot::teciompi::MPINonBlockingCommunicationCollection::MPINonBlockingReceiveVector<T>::~MPINonBlockingReceiveVector() { if (!isComplete()) { if (m_sizeRequest != MPI_REQUEST_NULL) MPI_Request_free(&m_sizeRequest); if (m_vecRequest != MPI_REQUEST_NULL) MPI_Request_free(&m_vecRequest);
 #if defined _DEBUG
std::cerr << "Warning: Uncompleted MPINonBlockingReceiveVector going out of scope." << std::endl;
 #endif
} } template <typename T> bool MPINonBlockingCommunicationCollection::MPINonBlockingReceiveVector<T>::tryToComplete(CompleteFunction completeFunction) { if (!m_complete) { if (m_sizeRequest != MPI_REQUEST_NULL) { int sizeFlag = 1; completeFunction(&m_sizeRequest, &sizeFlag, MPI_STATUS_IGNORE); if (sizeFlag == 1) { ___478(m_sizeRequest == MPI_REQUEST_NULL); m_vec.allocate(m_size); int ___3358 = MPI_Irecv(m_vec.begin(), m_vec.size(), mpiDatatype<T>(), m_src, m_vecTag, m_comm, &m_vecRequest); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingReceiveVector::tryToComplete", "MPI_Irecv", ___3358); } } if (m_sizeRequest == MPI_REQUEST_NULL && m_vecRequest != MPI_REQUEST_NULL) { int vecFlag; completeFunction(&m_vecRequest, &vecFlag, MPI_STATUS_IGNORE); } if (m_sizeRequest == MPI_REQUEST_NULL && m_vecRequest == MPI_REQUEST_NULL) m_complete = true; } return m_complete; } template <> bool MPINonBlockingCommunicationCollection::MPINonBlockingReceiveVector<___2479>::tryToComplete(CompleteFunction completeFunction) { if (!m_complete) { if (m_sizeRequest != MPI_REQUEST_NULL) { int sizeFlag = 1; completeFunction(&m_sizeRequest, &sizeFlag, MPI_STATUS_IGNORE); if (sizeFlag == 1) { ___478(m_sizeRequest == MPI_REQUEST_NULL); m_vec.allocate(m_size); int ___3358 = MPI_Irecv(m_vec.begin(), 2 * m_vec.size(), MPI_DOUBLE, m_src, m_vecTag, m_comm, &m_vecRequest); if (___3358 != MPI_SUCCESS) throwMPIError("MPINonBlockingReceiveVector::tryToComplete", "MPI_Irecv", ___3358); } } if (m_sizeRequest == MPI_REQUEST_NULL && m_vecRequest != MPI_REQUEST_NULL) { int vecFlag; completeFunction(&m_vecRequest, &vecFlag, MPI_STATUS_IGNORE); } if (m_sizeRequest == MPI_REQUEST_NULL && m_vecRequest == MPI_REQUEST_NULL) m_complete = true; } return m_complete; } template <typename T> bool MPINonBlockingCommunicationCollection::MPINonBlockingReceiveVector<T>::isComplete() { return tryToComplete(boost::bind(MPI_Test, _1, _2, _3)); } template <typename T> void MPINonBlockingCommunicationCollection::MPINonBlockingReceiveVector<T>::___4446() { tryToComplete(boost::bind(MPI_Wait, _1, _3)); } MPINonBlockingCommunicationCollection::MPINonBlockingCommunicationCollection(MPI_Comm communicator, int numRequests  ) : m_communicator(communicator) { REQUIRE(numRequests >= 0); if (numRequests > 0) m_requests.reserve(numRequests); } template<typename T> void tecplot::teciompi::MPINonBlockingCommunicationCollection::sendScalar(T const& ___4298, int dest, int tag) { m_requests.push_back(boost::make_shared<MPINonBlockingSendScalar<T> > (___4298, dest, tag, m_communicator)); } template<typename T> void tecplot::teciompi::MPINonBlockingCommunicationCollection::sendScalarCopy(T ___4298, int dest, int tag) { m_requests.push_back(boost::make_shared<MPINonBlockingSendScalarCopy<T> > (___4298, dest, tag, m_communicator)); } template<typename T> void tecplot::teciompi::MPINonBlockingCommunicationCollection::receiveScalar(T& ___4298, int ___3656, int tag) { m_requests.push_back(boost::make_shared<MPINonBlockingReceiveScalar<T> > (boost::ref(___4298), ___3656, tag, m_communicator)); } template<typename T> void MPINonBlockingCommunicationCollection::sendVector(SimpleVector<T> const& vec, int dest, int sizeTag, int vecTag) { REQUIRE(sizeTag != vecTag); m_requests.push_back(boost::make_shared<MPINonBlockingSendVector<T> > (boost::ref(vec), dest, sizeTag, vecTag, m_communicator)); } template<typename T> void MPINonBlockingCommunicationCollection::sendVectorCopy(SimpleVector<T> const& vec, int dest, int sizeTag, int vecTag) { REQUIRE(sizeTag != vecTag); m_requests.push_back(boost::make_shared<MPINonBlockingSendVectorCopy<T> > (boost::ref(vec), dest, sizeTag, vecTag, m_communicator)); } template<typename T> void tecplot::teciompi::MPINonBlockingCommunicationCollection::sendVectorCopy(std::vector<T> const& vec, int dest, int sizeTag, int vecTag) { REQUIRE(sizeTag != vecTag); m_requests.push_back(boost::make_shared<MPINonBlockingSendVectorCopy<T> > (boost::ref(vec), dest, sizeTag, vecTag, m_communicator)); } void MPINonBlockingCommunicationCollection::sendStringCopy(std::string const& str, int dest, int sizeTag, int vecTag) { REQUIRE(sizeTag != vecTag); SimpleVector<char> vec(str); m_requests.push_back(boost::make_shared<MPINonBlockingSendVectorCopy<char> > (boost::ref(vec), dest, sizeTag, vecTag, m_communicator)); } template<typename T> void MPINonBlockingCommunicationCollection::receiveVector(SimpleVector<T>& vec, int ___3656, int sizeTag, int vecTag) { REQUIRE(sizeTag != vecTag); m_requests.push_back(boost::make_shared<MPINonBlockingReceiveVector<T> > (boost::ref(vec), ___3656, sizeTag, vecTag, m_communicator)); } bool MPINonBlockingCommunicationCollection::isComplete() {
 #if defined MPI_INSTRUMENTATION
std::ofstream outputFile("iscomplete.txt"); bool isComplete = true; BOOST_FOREACH(boost::shared_ptr<MPINonBlockingCommunicationCollection::MPINonBlockingCommunication>& request, m_requests) { outputFile << request->m_sendOrReceive << ", " << request->m_other << ", "; outputFile.setf(std::ios::hex, std::ios::basefield); outputFile.setf(std::ios::showbase); outputFile << request->m_tag; outputFile.setf(std::ios::dec, std::ios::basefield); outputFile.unsetf(std::ios::showbase); if (request->isComplete()) { outputFile << ", true"; } else { outputFile << ", false"; isComplete = false; } outputFile << std::endl; } return isComplete;
 #else
BOOST_FOREACH(boost::shared_ptr<MPINonBlockingCommunicationCollection::MPINonBlockingCommunication>& request, m_requests) if (!request->isComplete()) return false; return true;
 #endif
} void MPINonBlockingCommunicationCollection::___4446() { size_t numRequests = m_requests.size(); for (size_t i = 0; i < numRequests; ++i) { bool complete = true; for (size_t i = 0; i < numRequests; ++i) complete = m_requests[i]->isComplete() && complete; if (complete) break; } BOOST_FOREACH(boost::shared_ptr<MPINonBlockingCommunicationCollection::MPINonBlockingCommunication>& request, m_requests) request->___4446(); m_requests.clear(); }
 #define INSTANTIATE_FOR_TYPE(T) \
 template void MPINonBlockingCommunicationCollection::sendScalar<T>(T const& ___4298, int dest, int tag); \
 template void MPINonBlockingCommunicationCollection::sendScalarCopy<T>(T ___4298, int dest, int tag); \
 template void MPINonBlockingCommunicationCollection::receiveScalar<T>(T& ___4298, int dest, int tag); \
 template void MPINonBlockingCommunicationCollection::sendVector<T>(SimpleVector<T> const& vec, int dest, int sizeTag, int vecTag); \
 template void MPINonBlockingCommunicationCollection::sendVectorCopy<T>(SimpleVector<T> const& vec, int dest, int sizeTag, int vecTag); \
 template void MPINonBlockingCommunicationCollection::sendVectorCopy<T>(std::vector<T> const& vec, int dest, int sizeTag, int vecTag); \
 template void MPINonBlockingCommunicationCollection::receiveVector<T>(SimpleVector<T>& vec, int ___3656, int sizeTag, int vecTag);
INSTANTIATE_FOR_TYPE(char) INSTANTIATE_FOR_TYPE(uint8_t) INSTANTIATE_FOR_TYPE(uint16_t) INSTANTIATE_FOR_TYPE(int32_t) INSTANTIATE_FOR_TYPE(uint32_t) INSTANTIATE_FOR_TYPE(int64_t) INSTANTIATE_FOR_TYPE(uint64_t) INSTANTIATE_FOR_TYPE(float) INSTANTIATE_FOR_TYPE(double) INSTANTIATE_FOR_TYPE(___2479)
 #undef INSTANTIATE_FOR_TYPE
}}
