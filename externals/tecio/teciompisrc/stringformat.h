 #pragma once
#include "ThirdPartyHeadersBegin.h"
#  include <iomanip>
#  include <sstream>
#  include <stdio.h>
#  include <string>
#include "ThirdPartyHeadersEnd.h"
namespace tecplot { template <typename T> std::string ___4187( T   ___4314, int precision = 6) { std::stringstream ss; ss << std::setprecision(precision) << ___4314; return ss.str(); } } namespace tecplot { template <typename T> std::string toString(T ___4314) { return tecplot::___4187(___4314); } } template <typename T> std::string ___4187(T ___4314) { return tecplot::___4187(___4314); }
