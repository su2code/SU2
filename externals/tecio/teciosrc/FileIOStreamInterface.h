 #pragma once
#include "ThirdPartyHeadersBegin.h"
#   include <string>
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
#include "CodeContract.h"
#include "basicTypes.h"
namespace tecplot { namespace ___3933 { class FileIOStreamInterface { public: virtual ___372 ___2041() const = 0; virtual ___372 close(bool ___3361) = 0; virtual ___1393 fileLoc() = 0; virtual ___372 ___3460() = 0; virtual ___372 ___3459(___1393 fileLoc) = 0; virtual ___372 seekToFileEnd() = 0; virtual std::string const& ___1394() const = 0; virtual void ___3494(___372 ___2002) = 0; virtual ___372 ___2002() const = 0; virtual void setDataFileType(DataFileType_e ___844) = 0; virtual DataFileType_e ___844() const = 0; virtual class FileIOStatistics& statistics() = 0; virtual ~FileIOStreamInterface() {} }; }}
