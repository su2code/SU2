#pragma once
#include "COutputModule.hpp"

class CMeshDeformModule final : public CSolverOutputModule {

public:
  explicit CMeshDeformModule(CConfig* config) : CSolverOutputModule(config->GetDeform_Mesh()){}

  void LoadHistoryData(CHistoryOutFieldManager& historyFields);

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields);

};
