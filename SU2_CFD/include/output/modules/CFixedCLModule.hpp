#include "COutputModule.hpp"

class CFixedCLModule final : public CSolverOutputModule {

public:
  explicit CFixedCLModule(CConfig* config, int nDim) : CSolverOutputModule(nDim, config->GetFixed_CL_Mode()){}

  void LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                       const IterationInfo& iterationInfo) override;

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields) override;

};
