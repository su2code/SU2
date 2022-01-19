 #pragma once
#include "ClassicZoneFaceNeighborWriter.h"
#include "ClassicOrderedZoneFaceNeighborGenerator.h"
#include "ClassMacros.h"
#include "ClassicZoneWriterAbstract.h"
namespace tecplot { namespace ___3933 { class ItemSetIterator; class ClassicOrderedZoneWriter : public ClassicZoneWriterAbstract { UNCOPYABLE_CLASS(ClassicOrderedZoneWriter) public: ClassicOrderedZoneWriter( ItemSetIterator&              varIter, ___4636                   zone, ___4636                   ___341, std::vector<___372> const& ___4564, ___372                     ___4499, ___37&                   ___36); virtual ~ClassicOrderedZoneWriter(); private: virtual uint64_t zoneConnectivityFileSize(bool ___2002); virtual ___372 writeZoneConnectivity(FileWriterInterface& szpltFile); ClassicOrderedZoneFaceNeighborGenerator m_faceNeighborGenerator; ClassicZoneFaceNeighborWriter m_faceNeighborWriter; ___372 ___4514( FileWriterInterface& file, ValueLocation_e      ___4326, ___4352           ___4336); }; }}
