#pragma once
#include "COutputModule.hpp"

class CMeshDeformModule final : public CSolverOutputModule {

public:
  explicit CMeshDeformModule(CConfig* config, int nDim) : CSolverOutputModule(nDim, config->GetDeform_Mesh()){}

  void LoadHistoryData(CHistoryOutFieldManager& historyFields);

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields);

};
