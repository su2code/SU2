 #pragma once
#include "ClassMacros.h"
#include "FileIOStatistics.h"
#include "FileIOStreamInterface.h"
namespace tecplot { namespace ___3933 { class FileIOStream : public FileIOStreamInterface { public: explicit FileIOStream(std::string const& ___1394); virtual ~FileIOStream(); virtual ___372 ___2041() const; virtual ___372 close(bool ___3361); virtual ___1393 fileLoc(); virtual ___372 ___3460(); virtual ___372 ___3459(___1393 fileLoc); virtual ___372 seekToFileEnd(); virtual std::string const& ___1394() const; virtual void ___3494(___372 ___2002); virtual ___372 ___2002() const; virtual void setDataFileType(DataFileType_e ___844); virtual DataFileType_e ___844() const; virtual class FileIOStatistics& statistics(); ___372 open(std::string const& ___2504); FILE* handle() const; private: FileIOStatistics  m_statistics; FILE*             m_fileHandle; std::string const ___2461; bool              m_isAscii; DataFileType_e    m_dataFileType; UNCOPYABLE_CLASS(FileIOStream) }; }}
