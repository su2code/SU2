#pragma once
#include "COutputModule.hpp"

class CMeshDeformModule final : public CSolverOutputModule {

public:
  explicit CMeshDeformModule(CConfig* config) : CSolverOutputModule(config->GetDeform_Mesh()){}

  void LoadHistoryData(COutFieldCollection& fieldCollection);

  void DefineHistoryFields(COutFieldCollection& fieldCollection);

};
