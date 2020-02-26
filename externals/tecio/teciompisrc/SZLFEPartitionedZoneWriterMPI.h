 #pragma once
#include "SZLFEPartitionedZoneWriter.h"
#include "ThirdPartyHeadersBegin.h"
#include <boost/scoped_ptr.hpp>
#include <mpi.h>
#include <vector>
#include "ThirdPartyHeadersEnd.h"
namespace tecplot { namespace ___3933 { class ___37; class ZoneInfoCache; class ItemSetIterator; }} namespace tecplot { namespace teciompi { class SZLFEPartitionedZoneWriterMPI: public ___3933::SZLFEPartitionedZoneWriter { public: SZLFEPartitionedZoneWriterMPI( ___3933::ItemSetIterator&       varIter, ___3933::___4636            zone, ___3933::___4636            ___341, std::vector<___372> const& ___4564, ___372                     ___4499, ___3933::___37&            ___36, ___3933::ZoneInfoCache&         zoneInfoCache, MPI_Comm                      communicator, int                           mainProcess); virtual ~SZLFEPartitionedZoneWriterMPI(); virtual ___2479 varMinMax(___3933::___4352 datasetVar); private: virtual uint64_t zoneDataFileSize(bool ___2002); virtual uint64_t zoneHeaderFileSize(bool ___2002); virtual ___372 writeZoneData(___3933::FileWriterInterface& szpltFile); virtual ___372 writeZoneHeader(___3933::FileWriterInterface& szpltFile); void createPartitionWriters(); void collateAndReturnResponses( class MPINonBlockingCommunicationCollection& communicationCollection); void exchangeGhostInfo(std::map<int, boost::shared_ptr<___3933::___1350> >& partitionInfoMap); struct Impl; boost::scoped_ptr<Impl> const m_impl; }; }}
