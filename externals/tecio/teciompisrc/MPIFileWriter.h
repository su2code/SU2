 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <boost/scoped_ptr.hpp>
#include <mpi.h>
#include <string>
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
#include "basicTypes.h"
#include "ClassMacros.h"
#include "FileWriterInterface.h"
namespace tecplot { namespace teciompi { class MPIFileWriter : public ___3933::FileWriterInterface { public: MPIFileWriter( std::string const& ___1394, MPI_Comm           comm, size_t             bufferSizeInMB = 8); virtual ~MPIFileWriter(); virtual ___372 ___2041() const; virtual ___372 close(bool ___3361); virtual ___3933::___1393 fileLoc(); virtual ___372 ___3460(); virtual ___372 ___3459(___3933::___1393 fileLoc); virtual ___372 seekToFileEnd(); virtual std::string const& ___1394() const; virtual void ___3494(___372 ___2002); virtual ___372 ___2002() const; virtual void setDataFileType(DataFileType_e ___844); virtual DataFileType_e ___844() const; virtual ___3933::FileIOStatistics& statistics(); virtual ___372 open(bool update); virtual size_t fwrite(void const* ___416, size_t size, size_t count); virtual int fprintf(char const* format, ...); class ScopedCaching { public: ScopedCaching(MPIFileWriter& fileWriter, size_t cacheSize); virtual ~ScopedCaching(); private: MPIFileWriter& m_fileWriter; }; private: struct Impl; boost::scoped_ptr<Impl> const m_impl; UNCOPYABLE_CLASS(MPIFileWriter) }; }}
