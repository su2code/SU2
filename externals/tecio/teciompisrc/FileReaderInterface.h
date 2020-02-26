 #pragma once
#include "ThirdPartyHeadersBegin.h"
#   include <cstddef>
#include "ThirdPartyHeadersEnd.h"
 #if defined getc
 #   undef getc 
 #endif
#include "FileIOStreamInterface.h"
namespace tecplot { namespace ___3933 { class ___1399 : public FileIOStreamInterface { public: virtual ___372 open() = 0; virtual size_t fread(void* ___416, size_t size, size_t count) = 0; virtual char* fgets(char* s, int size) = 0; virtual int feof() = 0; virtual int getc() = 0; virtual int ungetc(int c) = 0; virtual int fscanf(char const* format, void* ___3251) = 0; virtual int fscanf(char const* format, void* ptr1, void* ptr2) = 0; virtual int fscanf(char const* format, void* ptr1, void* ptr2, void* ptr3) = 0; virtual ~___1399() {} }; }}
