 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
namespace tecplot { namespace ___3933 { class FaceNeighborGeneratorAbstract; class ClassicZoneFaceNeighborWriter { public: ClassicZoneFaceNeighborWriter( FaceNeighborGeneratorAbstract& faceNeighborGenerator, ___4636 zone, ___4636 ___341); ___372 write(class FileWriterInterface& file); uint64_t sizeInFile(bool ___2002); private: FaceNeighborGeneratorAbstract const& m_faceNeighborGenerator; ___4636 const ___2677; ___4636 const m_baseZone; std::string const m_zoneNumberLabel; std::vector<int32_t> m_userFaceNeighbors; }; }}
