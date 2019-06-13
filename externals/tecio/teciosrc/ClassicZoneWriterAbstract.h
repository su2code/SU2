 #pragma once
#include "ClassicZoneHeaderWriter.h"
#include "ClassicZoneVariableWriter.h"
#include "ZoneWriterAbstract.h"
#include "ClassicZoneFileLocations.h"
namespace tecplot { namespace ___3933 { class ClassicZoneWriterAbstract : public ___4709 { public: ClassicZoneWriterAbstract( ItemSetIterator&              varIter, ___4636                   zone, ___4636                   ___341, std::vector<___372> const& ___4564, ___372                     ___4499, ___37&                   ___36); virtual ~ClassicZoneWriterAbstract(); protected: ClassicZoneVariableWriter m_variableWriter; ClassicZoneHeaderWriter m_headerWriter; ClassicZoneFileLocations m_zoneFileLocations; private: virtual uint64_t zoneConnectivityFileSize(bool ___2002) = 0; virtual uint64_t zoneDataFileSize(bool ___2002); virtual uint64_t zoneHeaderFileSize(bool ___2002); virtual ___372 writeZoneData(FileWriterInterface& szpltFile); virtual ___372 writeZoneConnectivity(FileWriterInterface& szpltFile) = 0; virtual ___372 writeZoneHeader(FileWriterInterface& szpltFile); }; }}
