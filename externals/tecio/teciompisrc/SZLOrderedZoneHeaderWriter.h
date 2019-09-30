 #pragma once
#include "SzlFileLoader.h"
#include "ZoneHeaderWriterAbstract.h"
namespace tecplot { namespace ___3933 { class ___37; class ___1881; class ItemSetIterator; class SZLOrderedZoneHeaderWriter : public ZoneHeaderWriterAbstract { public: SZLOrderedZoneHeaderWriter( ItemSetIterator&    varIter, ___4636         zone, ___4636         ___341, ___37&         ___36, ___1881 const&  ijkZoneInfo, ___1392 const& varFileLocs); virtual ~SZLOrderedZoneHeaderWriter(); virtual uint64_t sizeInFile(bool ___2002) const; virtual ___372 write(FileWriterInterface& fileWriter) const; private: ___1881 const& m_ijkZoneInfo; ___1392 const& ___2673; }; }}
