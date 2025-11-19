#pragma once
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"
#include <vector>
#include <cmath>
#include "../../../Common/include/CConfig.hpp"

inline su2double ComputeLinSysResRMS(const CSolver* solver, const CGeometry* geometry) {
    unsigned short nVar = solver->GetnVar();
    unsigned long nPointDomain = geometry->GetnPointDomain();
    std::vector<su2double> sumRes(nVar, 0.0);
    su2double localSum = 0.0;
    for (unsigned long i = 0; i < nPointDomain; ++i) {
        const su2double* res = solver->LinSysRes.GetBlock(i);
        for (unsigned short v = 0; v < nVar; ++v) {
            sumRes[v] += res[v] * res[v];
        }
    }
    for (unsigned short v = 0; v < nVar; ++v) localSum += sumRes[v];
    su2double globalSum = 0.0;
    SU2_MPI::Allreduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    unsigned long globalNPointDomain = 0;
    SU2_MPI::Allreduce(&nPointDomain, &globalNPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
    if (globalNPointDomain == 0) return 0.0;
    return std::sqrt(globalSum / (globalNPointDomain * nVar));
}
