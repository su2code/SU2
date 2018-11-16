 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
#include "FaceNeighborGeneratorAbstract.h"
namespace tecplot { namespace ___3933 { class ___37; class ClassicOrderedZoneFaceNeighborGenerator : public FaceNeighborGeneratorAbstract { public: explicit ClassicOrderedZoneFaceNeighborGenerator(___37& ___36); ___372 gatherUserFaceNeighbors(std::vector<int32_t>& userFaceNeighbors, ___4636 zone) const; }; }}
