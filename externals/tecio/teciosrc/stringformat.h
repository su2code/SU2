 #pragma once
#include "ThirdPartyHeadersBegin.h"
#  include <iomanip>
#  include <sstream>
#  include <stdio.h>
#  include <string>
#include "ThirdPartyHeadersEnd.h"
 #ifdef _MSC_VER
 #define snprintf std_snprintf
namespace { inline int std_vsnprintf(char* str, size_t size, const char* format, va_list ap) { int count = -1; if (size != 0) count = _vsnprintf_s(str, size, _TRUNCATE, format, ap); if (count == -1) count = _vscprintf(format, ap); return count; } inline int std_snprintf(char* str, size_t size, const char* format, ...) { int count; va_list ap; va_start(ap, format); count = std_vsnprintf(str, size, format, ap); va_end(ap); return count; } }
 #endif 
namespace tecplot { template <typename T> std::string ___4187( T   ___4314, int precision = 6) { std::stringstream ss; ss << std::setprecision(precision) << ___4314; return ss.str(); } } namespace tecplot { template <typename T> std::string toString(T ___4314) { return tecplot::___4187(___4314); } } template <typename T> std::string ___4187(T ___4314) { return tecplot::___4187(___4314); }
