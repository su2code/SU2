 #pragma once
#include "ThirdPartyHeadersBegin.h"
#   include <cstdio>
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
#include "basicTypes.h"
#include "FileIOStream.h"
#include "FileReaderInterface.h"
namespace tecplot { namespace ___3933 { class FileStreamReader : public ___1399 { public: explicit FileStreamReader(std::string const& ___1394); virtual ~FileStreamReader(); virtual ___372 ___2041() const; virtual ___372 close(bool ___3361); virtual ___1393 fileLoc(); virtual ___372 ___3460(); virtual ___372 ___3459(___1393 fileLoc); virtual ___372 seekToFileEnd(); virtual std::string const& ___1394() const; virtual void ___3494(___372 ___2002); virtual ___372 ___2002() const; virtual void setDataFileType(DataFileType_e ___844); virtual DataFileType_e ___844() const; virtual class FileIOStatistics& statistics(); virtual ___372 open(); virtual size_t fread(void* ___416, size_t size, size_t count); virtual char* fgets(char* s, int size); virtual int feof(); virtual int getc(); virtual int ungetc(int c); virtual int fscanf(char const* format, void* ___3251); virtual int fscanf(char const* format, void* ptr1, void* ptr2); virtual int fscanf(char const* format, void* ptr1, void* ptr2, void* ptr3); private: FileIOStream m_fileIOStream; UNCOPYABLE_CLASS(FileStreamReader); }; }}
