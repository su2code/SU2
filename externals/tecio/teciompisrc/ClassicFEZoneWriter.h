 #pragma once
#include "ClassicFEZoneConnectivityWriter.h"
#include "ClassicFEZoneFaceNeighborGenerator.h"
#include "ClassicZoneFaceNeighborWriter.h"
#include "ClassicZoneWriterAbstract.h"
namespace tecplot { namespace ___3933 { class ClassicFEZoneWriter : public ClassicZoneWriterAbstract { UNCOPYABLE_CLASS(ClassicFEZoneWriter) public: ClassicFEZoneWriter( ItemSetIterator&              varIter, ___4636                   zone, ___4636                   ___341, std::vector<___372> const& ___4564, ___372                     ___4499, ___37&                   ___36); virtual ~ClassicFEZoneWriter(); private: virtual uint64_t zoneConnectivityFileSize(bool ___2002); virtual ___372 writeZoneConnectivity(FileWriterInterface& szpltFile); ClassicFEZoneConnectivityWriter m_connectivityWriter; ClassicFEZoneFaceNeighborGenerator m_faceNeighborGenerator; ClassicZoneFaceNeighborWriter m_faceNeighborWriter; std::string m_zoneNumberLabel; }; }}
