 #pragma once
#include "ThirdPartyHeadersBegin.h"
#   include <string>
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
#include "basicTypes.h"
#include "FileIOStream.h"
#include "FileWriterInterface.h"
namespace tecplot { namespace ___3933 { class FileStreamWriter : public FileWriterInterface { public: explicit FileStreamWriter(std::string const& ___1394); virtual ~FileStreamWriter(); virtual ___372 ___2041() const; virtual ___372 close(bool ___3361); virtual ___1393 fileLoc(); virtual ___372 ___3460(); virtual ___372 ___3459(___1393 fileLoc); virtual ___372 seekToFileEnd(); virtual std::string const& ___1394() const; virtual void ___3494(___372 ___2002); virtual ___372 ___2002() const; virtual void setDataFileType(DataFileType_e ___844); virtual DataFileType_e ___844() const; virtual class FileIOStatistics& statistics(); virtual ___372 open(bool update); virtual size_t fwrite(void const* ___416, size_t size, size_t count); virtual int fprintf(char const* format, ...); private: FileIOStream m_fileIOStream; UNCOPYABLE_CLASS(FileStreamWriter) }; }}
