#include "ZoneWriterAbstract.h"
#include "AltTecUtil.h"
#include "ThirdPartyHeadersBegin.h"
#include <stdexcept>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
 #if defined OUTPUT_TIMES
 #if defined _WINDOWS
#include <sys\types.h>
#include <sys\timeb.h>
 #else
#include <sys/time.h>
 #endif
 #endif
#include "ThirdPartyHeadersEnd.h"
#include "ItemSetIterator.h"
#include "writeValueArray.h"
namespace tecplot { namespace ___3933 { ___4709::___4709( ItemSetIterator&              varIter, ___4636                   zone, ___4636                   ___341, std::vector<___372> const& ___4564, ___372                     ___4499, ___37&                   ___36) : m_varIter(varIter) , ___2677(zone) , m_baseZone(___341) , m_writeVariables(___4564) , m_writeConnectivity(___4499) , ___2337(___36) , m_zoneHeaderFilePosition(___330) , m_zoneFileSize(0) { ___478(0 <= m_baseZone && m_baseZone <= ___2677); if (!___2337.___896()) throw std::logic_error("Attempted to write zone for non-existent data set."); else if (___2677 < 0 || !___2337.___4638(___2677 + 1)) throw std::logic_error("Attempted to write non-existent zone."); } ___4709::~___4709() {} uint64_t ___4709::zoneFileSize(bool ___2002) { if (!m_zoneFileSize) m_zoneFileSize = zoneConnectivityFileSize(___2002) + zoneDataFileSize(___2002) + zoneHeaderFileSize(___2002); return m_zoneFileSize; } ___372 ___4709::writeZone( FileWriterInterface& szpltFile, ___1393            fileLoc) { ASSERT_ONLY(uint64_t fileSize1 = fileLoc); ASSERT_ONLY(uint64_t predictedFileSize2 = fileSize1 + zoneConnectivityFileSize(szpltFile.___2002() == ___4226)); szpltFile.___3459(fileLoc); if (writeZoneConnectivity(szpltFile)) { ASSERT_ONLY(uint64_t fileSize2 = szpltFile.fileLoc()); ___478(fileSize2 == predictedFileSize2); ASSERT_ONLY(uint64_t predictedFileSize3 = fileSize2 + zoneDataFileSize(szpltFile.___2002() == ___4226)); if (writeZoneData(szpltFile)) { ASSERT_ONLY(uint64_t fileSize3 = szpltFile.fileLoc()); ___478(fileSize3 == predictedFileSize3); ASSERT_ONLY(uint64_t predictedFileSize4 = fileSize3 + zoneHeaderFileSize(szpltFile.___2002() == ___4226)); m_zoneHeaderFilePosition = szpltFile.fileLoc(); ___372 ___3358 = writeZoneHeader(szpltFile); ASSERT_ONLY(uint64_t fileSize4 = szpltFile.fileLoc()); ___478(fileSize4 == predictedFileSize4); return ___3358; } } return ___1305; } uint64_t ___4709::getZoneHeaderFilePosition() const { if (m_zoneHeaderFilePosition == ___330) throw std::runtime_error("Internal error: getZoneHeaderFilePosition() called before writeZone()."); return m_zoneHeaderFilePosition; } ___2479 ___4709::varMinMax(___4352 datasetVar) { REQUIRE(m_varIter.baseItem() <= datasetVar && datasetVar - m_varIter.baseItem() < m_varIter.___2812()); double minValue; double maxValue; ___2337.___913(___2677 + 1, datasetVar + 1, &minValue, &maxValue); return ___2479(minValue, maxValue); }
 #if defined OUTPUT_TIMES
uint64_t ___4709::___717(void) {
 #if defined UNIXX
 #  if defined LINUXALPHA
{ struct timeb now; ftime(&now); return now.time * 1000L + (long)now.millitm; }
 #  elif defined LINUX
{ struct timeval now; gettimeofday(&now, NULL); return now.tv_sec * 1000L + now.tv_usec / 1000L; }
 #  else
{ struct timeval now; struct timezone ___4179; gettimeofday(&now, &___4179); return now.tv_sec * 1000L + now.tv_usec / 1000L; }
 #  endif
 #endif 
 #if defined _WIN32
{ struct _timeb curTimeB; _ftime(&curTimeB); return curTimeB.time * 1000LL + curTimeB.millitm; }
 #endif 
}
 #endif
}}
