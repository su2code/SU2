#pragma once
#include "COutputModule.hpp"

class CMeshDeformModule final : public CSolverOutputModule {

public:
  explicit CMeshDeformModule(CConfig* config, int nDim) : CSolverOutputModule(nDim, config->GetDeform_Mesh()){}

  void LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                       const IterationInfo& iterationInfo) override;

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields) override;

};
