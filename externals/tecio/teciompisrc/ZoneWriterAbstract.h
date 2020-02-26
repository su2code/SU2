 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <string>
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
#include "MinMax.h"
namespace tecplot { namespace ___3933 { class ___37; class FileWriterInterface; class ItemSetIterator; class ___4709 { public: ___4709( ItemSetIterator&              varIter, ___4636                   zone, ___4636                   ___341, std::vector<___372> const& ___4564, ___372                     ___4499, ___37&                   ___36); virtual ~___4709(); uint64_t zoneFileSize(bool ___2002); ___372 writeZone( FileWriterInterface& szpltFile, ___1393            fileLoc); virtual uint64_t getZoneHeaderFilePosition() const; virtual ___2479 varMinMax(___4352 ___4336);
 #if defined OUTPUT_TIMES
static uint64_t ___717();
 #endif
protected: ItemSetIterator&       m_varIter; ___4636 const      ___2677; ___4636 const      m_baseZone; std::vector<___372> m_writeVariables; ___372 const        m_writeConnectivity; ___37&            ___2337; uint64_t               m_zoneHeaderFilePosition; uint64_t               m_zoneFileSize; private: virtual uint64_t zoneConnectivityFileSize(bool ___2002) = 0; virtual uint64_t zoneDataFileSize(bool ___2002) = 0; virtual uint64_t zoneHeaderFileSize(bool ___2002) = 0; virtual ___372 writeZoneConnectivity(FileWriterInterface& szpltFile) = 0; virtual ___372 writeZoneData(FileWriterInterface& szpltFile) = 0; virtual ___372 writeZoneHeader(FileWriterInterface& szpltFile) = 0; }; }}
