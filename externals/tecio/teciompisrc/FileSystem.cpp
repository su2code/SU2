#include "FileSystem.h"
#include "CodeContract.h"
 #if defined MSWIN
#   include "UnicodeStringUtils.h"
 #endif
 #if !defined NO_THIRD_PARTY_LIBS
#  include "ThirdPartyHeadersBegin.h"
#    include <boost/filesystem.hpp>
#    include <boost/system/error_code.hpp>
#  include "ThirdPartyHeadersEnd.h"
 #endif
 #if !defined NO_THIRD_PARTY_LIBS
namespace bs  = boost::system; namespace bfs = boost::filesystem;
 #endif
namespace tecplot { namespace filesystem { FILE* fileOpen(std::string const& ___1394, std::string const& ___2504) { REQUIRE(!___1394.empty()); REQUIRE(!___2504.empty());
 #if defined MSWIN
return _wfopen(tecplot::utf8ToWideString(___1394).c_str(), tecplot::utf8ToWideString(___2504).c_str());
 #else
return fopen(___1394.c_str(), ___2504.c_str());
 #endif
} FILE* fileReopen(std::string const& ___1394, std::string const& ___2504, FILE* file) { REQUIRE(!___1394.empty()); REQUIRE(!___2504.empty()); REQUIRE(VALID_REF(file));
 #if defined MSWIN
return _wfreopen(tecplot::utf8ToWideString(___1394).c_str(), tecplot::utf8ToWideString(___2504).c_str(), file);
 #else
return freopen(___1394.c_str(), ___2504.c_str(), file);
 #endif
} int fileRename(std::string const& ___1394, std::string const& newFileName) { REQUIRE(!___1394.empty()); REQUIRE(!newFileName.empty());
 #if defined MSWIN
return _wrename(tecplot::utf8ToWideString(___1394).c_str(), tecplot::utf8ToWideString(newFileName).c_str());
 #else
return rename(___1394.c_str(), newFileName.c_str());
 #endif
} int fileRemove(std::string const& ___1394) { REQUIRE(!___1394.empty());
 #if defined MSWIN
return _wremove(tecplot::utf8ToWideString(___1394).c_str());
 #else
return remove(___1394.c_str());
 #endif
}
 #if !defined NO_THIRD_PARTY_LIBS
bool fileExists(boost::filesystem::path const& filePath) { REQUIRE("filepath can be valid or empty"); bool ___3358 = filePath.has_filename(); if (___3358) { bs::error_code errorCode; bfs::file_status fileStatus; fileStatus = bfs::status(filePath,errorCode); ___3358 = bfs::exists(filePath, errorCode) && !bfs::is_directory(filePath, errorCode); } return ___3358; } bool dirExists(boost::filesystem::path const& dirPath) { REQUIRE("dirpath can be anything - including empty"); bool ___3358 = dirPath.has_filename(); if (___3358) { bs::error_code errorCode; ___3358 = bfs::exists(dirPath, errorCode) && bfs::is_directory(dirPath, errorCode); } return ___3358; }
 #endif
}}
