 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <mpi.h>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
#include "ZoneWriterFactory.h"
namespace tecplot { namespace ___3933 { class ___37; class ZoneInfoCache; class ___4709; class ItemSetIterator; }} namespace tecplot { namespace teciompi { class ZoneWriterFactoryMPI : public ___3933::___4711 { public: ZoneWriterFactoryMPI( ___3933::ZoneInfoCache& zoneInfoCache, ___3933::___37&    ___36, MPI_Comm              communicator, int                   mainProcess); boost::shared_ptr<___3933::___4709> ___4708( ___3933::ItemSetIterator&       varIter, ___3933::___4636            zone, ___3933::___4636            ___341, std::vector<___372> const& ___4564, ___372                     ___4499); private: MPI_Comm m_communicator; int m_mainProcess; }; }}
