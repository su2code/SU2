 #pragma once
#include "ThirdPartyHeadersBegin.h"
#   include <cstddef>
#include "ThirdPartyHeadersEnd.h"
 #if defined getc
 #   undef getc 
 #endif
#include "FileIOStreamInterface.h"
namespace tecplot { namespace ___3933 { class FileWriterInterface : public FileIOStreamInterface { public: virtual ___372 open(bool update) = 0; virtual size_t fwrite(void const* ___416, size_t size, size_t count) = 0; virtual int fprintf(char const* format, ...) = 0; virtual ~FileWriterInterface() {} }; }}
