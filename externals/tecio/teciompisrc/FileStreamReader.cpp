#include "MASTER.h"
#include "GLOBAL.h"
#include "CodeContract.h"
#include "showMessage.h"
#include "FileStreamReader.h"
namespace tecplot { namespace ___3933 { FileStreamReader::FileStreamReader(std::string const& ___1394) : m_fileIOStream(___1394) { } FileStreamReader::~FileStreamReader() { m_fileIOStream.close(true); } ___372 FileStreamReader::___2041() const { return m_fileIOStream.___2041(); } ___372 FileStreamReader::open() { REQUIRE(!___2041()); ___372 ___2039 = m_fileIOStream.open("rb"); if (!___2039) ___1186("Cannot read file %s", ___1394().c_str());
 #ifdef PROFILE_FILE_ACCESS
___478(statistics().numFSeeksPerformed == 0 && statistics().numReadWritesPerformed == 0 && statistics().___2780 == 0);
 #endif
return ___2039; } ___372 FileStreamReader::close(bool ___3361) { return m_fileIOStream.close(___3361); } ___1393 FileStreamReader::fileLoc() { REQUIRE(___2041()); return m_fileIOStream.fileLoc(); } ___372 FileStreamReader::___3460() { REQUIRE(___2041()); return m_fileIOStream.___3460(); } ___372 FileStreamReader::___3459(___1393 fileLoc) { REQUIRE(___2041()); REQUIRE(fileLoc != ___330); return m_fileIOStream.___3459(fileLoc); } ___372 FileStreamReader::seekToFileEnd() { REQUIRE(___2041()); return m_fileIOStream.seekToFileEnd(); } std::string const& FileStreamReader::___1394() const { return m_fileIOStream.___1394(); } void FileStreamReader::___3494(___372 ___2002) { REQUIRE(VALID_BOOLEAN(___2002)); m_fileIOStream.___3494(___2002); } ___372 FileStreamReader::___2002() const { return m_fileIOStream.___2002(); } void FileStreamReader::setDataFileType(DataFileType_e ___844) { REQUIRE(VALID_ENUM(___844, DataFileType_e)); m_fileIOStream.setDataFileType(___844); } DataFileType_e FileStreamReader::___844() const { return m_fileIOStream.___844(); } class FileIOStatistics& FileStreamReader::statistics() { return m_fileIOStream.statistics(); } size_t FileStreamReader::fread(void* ___416, size_t size, size_t count) { REQUIRE(VALID_REF(___416)); REQUIRE(___2041()); size_t ___3358 = ::fread(___416, size, count, m_fileIOStream.handle());
 #ifdef PROFILE_FILE_ACCESS
if ( ___3358 > 0 ) { statistics().numReadWritesPerformed++; statistics().___2780 += ___3358*size; }
 #endif
return ___3358; } char* FileStreamReader::fgets(char* s, int size) { REQUIRE(VALID_REF(s)); REQUIRE(size >= 0); REQUIRE(___2041()); char* ___3358 = ::fgets(s, size, m_fileIOStream.handle());
 #ifdef PROFILE_FILE_ACCESS
if ( ___3358 != NULL ) { statistics().numReadWritesPerformed++; statistics().___2780 += ___3358 ? strlen(s) : 0; }
 #endif
return ___3358; } int FileStreamReader::feof() { REQUIRE(___2041()); return ::feof(m_fileIOStream.handle()); } int FileStreamReader::getc() { REQUIRE(___2041()); int ___3358 = ::getc(m_fileIOStream.handle());
 #ifdef PROFILE_FILE_ACCESS
statistics().numReadWritesPerformed++; statistics().___2780++;
 #endif
return ___3358; } int FileStreamReader::ungetc(int c) { REQUIRE(___2041()); int ___3358 = ::ungetc(c, m_fileIOStream.handle());
 #ifdef PROFILE_FILE_ACCESS
___478(statistics().___2780>0); statistics().___2780--;
 #endif
return ___3358; } int FileStreamReader::fscanf(char const* format, void* ___3251) { REQUIRE(___2041()); REQUIRE(VALID_NON_ZERO_LEN_STR(format)); REQUIRE(VALID_REF(___3251));
 #ifdef PROFILE_FILE_ACCESS
___1393 startLoc = fileLoc();
 #endif
int ___3358 = ::fscanf(m_fileIOStream.handle(), format, ___3251);
 #ifdef PROFILE_FILE_ACCESS
___1393 endLoc = fileLoc(); statistics().numReadWritesPerformed++; statistics().___2780 += (endLoc-startLoc);
 #endif
return ___3358; } int FileStreamReader::fscanf(char const* format, void* ptr1, void* ptr2) { REQUIRE(___2041()); REQUIRE(VALID_NON_ZERO_LEN_STR(format)); REQUIRE(VALID_REF(ptr1)); REQUIRE(VALID_REF(ptr2));
 #ifdef PROFILE_FILE_ACCESS
___1393 startLoc = fileLoc();
 #endif
int ___3358 = ::fscanf(m_fileIOStream.handle(), format, ptr1, ptr2);
 #ifdef PROFILE_FILE_ACCESS
___1393 endLoc = fileLoc(); statistics().numReadWritesPerformed++; statistics().___2780 += (endLoc-startLoc);
 #endif
return ___3358; } int FileStreamReader::fscanf(char const* format, void* ptr1, void* ptr2, void* ptr3) { REQUIRE(___2041()); REQUIRE(VALID_NON_ZERO_LEN_STR(format)); REQUIRE(VALID_REF(ptr1)); REQUIRE(VALID_REF(ptr2)); REQUIRE(VALID_REF(ptr3));
 #ifdef PROFILE_FILE_ACCESS
___1393 startLoc = fileLoc();
 #endif
int ___3358 = ::fscanf(m_fileIOStream.handle(), format, ptr1, ptr2, ptr3);
 #ifdef PROFILE_FILE_ACCESS
___1393 endLoc = fileLoc(); statistics().numReadWritesPerformed++; statistics().___2780 += (endLoc-startLoc);
 #endif
return ___3358; } }}
