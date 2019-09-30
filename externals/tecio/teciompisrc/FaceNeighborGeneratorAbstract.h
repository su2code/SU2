 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
namespace tecplot { namespace ___3933 { class FaceNeighborGeneratorAbstract { public: explicit FaceNeighborGeneratorAbstract(class ___37& ___36); virtual ~FaceNeighborGeneratorAbstract() {} virtual ___372 gatherUserFaceNeighbors(std::vector<int32_t>& userFaceNeighbors, ___4636 zone) const = 0; protected: void appendUserFaceNeighborsForCellFace( std::vector<int32_t>& userFaceNeighbors, ___1292 ___1274, ___3501 ___1153, ___4636 zone, int32_t ___449, int32_t face) const; class ___37& ___2337; }; }}
