 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <boost/scoped_ptr.hpp>
#include <mpi.h>
#include <map>
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "SZLOrderedPartitionedZoneWriter.h"
namespace tecplot { namespace ___3933 { class ___1881; } } namespace tecplot { namespace teciompi { class SZLOrderedPartitionedZoneWriterMPI : public ___3933::SZLOrderedPartitionedZoneWriter { public: SZLOrderedPartitionedZoneWriterMPI( ___3933::ItemSetIterator&        varIter, ___3933::___4636             zone, ___3933::___4636             ___341, std::vector<___372> const&  ___4564, ___372                      ___4499, ___3933::___37&             ___36, ___3933::ZoneInfoCache&          zoneInfoCache, MPI_Comm                       communicator, int                            mainProcess); virtual ~SZLOrderedPartitionedZoneWriterMPI(); virtual ___2479 varMinMax(___3933::___4352 datasetVar); private: virtual uint64_t zoneDataFileSize(bool ___2002); virtual uint64_t zoneHeaderFileSize(bool ___2002); virtual ___372 writeZoneData(___3933::FileWriterInterface& szpltFile); virtual ___372 writeZoneHeader(___3933::FileWriterInterface& szpltFile); struct Impl; boost::scoped_ptr<Impl> m_impl; }; }}
