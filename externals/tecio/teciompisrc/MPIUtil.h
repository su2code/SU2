 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <mpi.h>
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "AltTecUtil.h"
#include "basicTypes.h"
#include "StandardIntegralTypes.h"
namespace tecplot { namespace teciompi { bool everyRankAppearsOnce(MPI_Comm comm, std::vector<int32_t> const& ranks); bool everyRankOwnsOnePartition(MPI_Comm comm, ___3933::___37& ___36, ___3933::___4636 ___4658); void gatherScatterPartitionFileLocs(___3933::___1393& finalFileLoc, ___3933::___1393& localPartitionFileLoc, ___3933::___37 &___36, ___3933::___4636 const zone, int32_t localProcess, uint64_t localPartitionFileSize, int32_t mainProcess, MPI_Comm comm); }}
