#include "MASTER.h"
#include "GLOBAL.h"
#include "CodeContract.h"
#include "showMessage.h"
#include "FileStreamWriter.h"
namespace tecplot { namespace ___3933 { FileStreamWriter::FileStreamWriter(std::string const& ___1394) : m_fileIOStream(___1394) { } FileStreamWriter::~FileStreamWriter() { m_fileIOStream.close(true); } ___372 FileStreamWriter::___2041() const { return m_fileIOStream.___2041(); } ___372 FileStreamWriter::open(bool update) { REQUIRE(!___2041()); ___372 ___2039; if (update) { ___2039 = m_fileIOStream.open("rb+"); if (!___2039) ___2039 = m_fileIOStream.open("wb+"); } else { ___2039 = m_fileIOStream.open("wb+"); } if (___2039) {
 #if !defined NDEBUG
{ if (m_fileIOStream.___2002()) if (setvbuf(m_fileIOStream.handle(), NULL, _IONBF, 0) != 0) ___2039 = ___1305; }
 #endif
} else { ___1186( "Cannot write to file %s", m_fileIOStream.___1394().c_str()); }
 #ifdef PROFILE_FILE_ACCESS
___478(statistics().numFSeeksPerformed == 0 && statistics().numReadWritesPerformed == 0 && statistics().___2780 == 0);
 #endif
return ___2039; } ___372 FileStreamWriter::close(bool ___3361) { return m_fileIOStream.close(___3361); } ___1393 FileStreamWriter::fileLoc() { REQUIRE(___2041()); return m_fileIOStream.fileLoc(); } ___372 FileStreamWriter::___3460() { REQUIRE(___2041()); return m_fileIOStream.___3460(); } ___372 FileStreamWriter::___3459(___1393 fileLoc) { REQUIRE(___2041()); REQUIRE(fileLoc != ___330); return m_fileIOStream.___3459(fileLoc); } ___372 FileStreamWriter::seekToFileEnd() { REQUIRE(___2041()); return m_fileIOStream.seekToFileEnd(); } std::string const& FileStreamWriter::___1394() const { return m_fileIOStream.___1394(); } void FileStreamWriter::___3494(___372 ___2002) { REQUIRE(VALID_BOOLEAN(___2002)); m_fileIOStream.___3494(___2002); } ___372 FileStreamWriter::___2002() const { return m_fileIOStream.___2002(); } void FileStreamWriter::setDataFileType(DataFileType_e ___844) { REQUIRE(VALID_ENUM(___844, DataFileType_e)); m_fileIOStream.setDataFileType(___844); } DataFileType_e FileStreamWriter::___844() const { return m_fileIOStream.___844(); } class FileIOStatistics& FileStreamWriter::statistics() { return m_fileIOStream.statistics(); } size_t FileStreamWriter::fwrite(void const* ___416, size_t size, size_t count) { REQUIRE(___2041()); REQUIRE(VALID_REF(___416)); size_t ___3358 = ::fwrite(___416, size, count, m_fileIOStream.handle());
 #ifdef PROFILE_FILE_ACCESS
if ( ___3358 > 0 ) { statistics().numReadWritesPerformed++; statistics().___2780 += ___3358*size; }
 #endif
return ___3358; } int FileStreamWriter::fprintf(char const* format, ...) { REQUIRE(___2041()); REQUIRE(VALID_NON_ZERO_LEN_STR(format)); va_list args; va_start(args, format); int ___3358 = ::vfprintf(m_fileIOStream.handle(), format, args); va_end (args);
 #ifdef PROFILE_FILE_ACCESS
if (___3358 > 0) { statistics().numReadWritesPerformed++; statistics().___2780 += ___3358; }
 #endif
return ___3358; } }}
