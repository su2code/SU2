 #pragma once
#include "ThirdPartyHeadersBegin.h"
#  include <stdio.h>
#  include <string>
#include "ThirdPartyHeadersEnd.h"
#include "StandardIntegralTypes.h"
namespace boost { namespace filesystem { class path; }} namespace tecplot { namespace filesystem { FILE* fileOpen(std::string const& ___1394, std::string const& ___2504); FILE* fileReopen(std::string const& ___1394, std::string const& ___2504, FILE* file); int   fileRename(std::string const& ___1394, std::string const& newFileName); int   fileRemove(std::string const& ___1394); bool  fileSize(char const* ___1394, uint64_t& sizeResult); bool  fileSize(std::string const& ___1394, uint64_t& sizeResult);
 #if !defined NO_THIRD_PARTY_LIBS
bool fileExists(boost::filesystem::path const& filePath); bool dirExists(boost::filesystem::path const& dirPath);
 #endif
}}
