#include "MASTER.h"
#include "GLOBAL.h"
#include "CodeContract.h"
#include "FileIOStream.h"
#include "FileSystem.h"
namespace tecplot { namespace ___3933 { FileIOStream::FileIOStream(std::string const& ___1394) : m_fileHandle(NULL) , ___2461(___1394) , m_isAscii(false) , m_dataFileType(___845) { REQUIRE(!___2461.empty()); } FileIOStream::~FileIOStream() { close(true); } ___372 FileIOStream::___2041() const { return m_fileHandle != NULL; } ___372 FileIOStream::close(bool ___3361) { ___372 ___3358 = ___4226; if (___2041()) { ___4195(m_fileHandle); m_fileHandle = NULL; if (!___3361) { ___478(!___1394().empty()); if (remove(___1394().c_str()) != 0) ___3358 = ___1305; }
 #ifdef PROFILE_FILE_ACCESS
m_statistics.numFSeeksPerformed = 0; m_statistics.numReadWritesPerformed = 0; m_statistics.___2780 = 0;
 #endif
} else { ___3358 = ___1305; } ENSURE(VALID_BOOLEAN(___3358)); return ___3358; } ___1393 FileIOStream::fileLoc() { REQUIRE(___2041()); ___1393 fileLoc = ___1393(___4201(m_fileHandle)); ENSURE(fileLoc != ___330); return fileLoc; } ___372 FileIOStream::___3460() { REQUIRE(___2041());
 #ifdef PROFILE_FILE_ACCESS
m_statistics.numFSeeksPerformed++;
 #endif
return ___4200(m_fileHandle, ___1393(0), SEEK_SET) == 0; } ___372 FileIOStream::___3459(___1393 fileLoc) { REQUIRE(___2041()); REQUIRE(fileLoc != ___330);
 #ifdef PROFILE_FILE_ACCESS
m_statistics.numFSeeksPerformed++;
 #endif
return ___4200(m_fileHandle, fileLoc, SEEK_SET) == 0; } ___372 FileIOStream::seekToFileEnd() { REQUIRE(___2041());
 #ifdef PROFILE_FILE_ACCESS
m_statistics.numFSeeksPerformed++;
 #endif
return ___4200(m_fileHandle, ___1393(0), SEEK_END) == 0; } std::string const& FileIOStream::___1394() const { return ___2461; } void FileIOStream::___3494(___372 ___2002) { REQUIRE(VALID_BOOLEAN(___2002)); m_isAscii = (___2002 == ___4226); } ___372 FileIOStream::___2002() const { return m_isAscii; } void FileIOStream::setDataFileType(DataFileType_e ___844) { REQUIRE(VALID_ENUM(___844, DataFileType_e)); m_dataFileType = ___844; } DataFileType_e FileIOStream::___844() const { return m_dataFileType; } class FileIOStatistics& FileIOStream::statistics() { return m_statistics; } ___372 FileIOStream::open(std::string const& ___2504) { REQUIRE(!___1394().empty()); REQUIRE(!___2041()); m_fileHandle = tecplot::filesystem::fileOpen(___1394(), ___2504); ___372 ___2039 = m_fileHandle != NULL;
 #ifdef PROFILE_FILE_ACCESS
___478(m_statistics.numFSeeksPerformed == 0 && m_statistics.numReadWritesPerformed == 0 && m_statistics.___2780 == 0);
 #endif
return ___2039; } FILE* FileIOStream::handle() const { return m_fileHandle; } }}
