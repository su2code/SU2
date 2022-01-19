 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <boost/scoped_ptr.hpp>
#include <mpi.h>
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
#include "basicTypes.h"
#include "ClassMacros.h"
#include "FileIOStreamInterface.h"
namespace tecplot { namespace ___3933 { class FileIOStatistics; } } namespace tecplot { namespace teciompi { class MPIFileIOStream : public ___3933::FileIOStreamInterface { public: MPIFileIOStream(std::string const& ___1394, MPI_Comm comm); virtual ~MPIFileIOStream(); virtual ___372 ___2041() const; virtual ___372 close(bool ___3361); virtual ___3933::___1393 fileLoc(); virtual ___372 ___3460(); virtual ___372 ___3459(___3933::___1393 fileLoc); virtual ___372 seekToFileEnd(); virtual std::string const& ___1394() const; virtual void ___3494(___372 ___2002); virtual ___372 ___2002() const; virtual void setDataFileType(DataFileType_e ___844); virtual DataFileType_e ___844() const; virtual ___3933::FileIOStatistics& statistics(); ___372 open(int ___2504); MPI_File fileHandle() const; MPI_Comm comm() const; private: struct Impl; boost::scoped_ptr<Impl> m_impl; UNCOPYABLE_CLASS(MPIFileIOStream) }; }}
