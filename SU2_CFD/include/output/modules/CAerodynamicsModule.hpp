#pragma once
#include "COutputModule.hpp"

class CAerodynamicsModule : public CSolverOutputModule{

protected:
  su2double factor, Alpha, Beta, RefLength, Gas_Constant, Gamma, Prandtl_Lam;
  bool axisymmetric;

  su2double Force[3], Moment[3], SkinFriction[3], Y_Plus;

public:

  explicit CAerodynamicsModule(CConfig *config, int nDim);

  void LoadHistoryDataPerSurface(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                                 const IterationInfo& iterationInfo) override;

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields) override;

  void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) override;

  void LoadSurfaceData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                       const IterationInfo& iterationInfo, const PointInfo& pointInfo) override;

};

