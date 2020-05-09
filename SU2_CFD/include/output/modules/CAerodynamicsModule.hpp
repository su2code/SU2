#pragma once
#include "COutputModule.hpp"

class CAerodynamicsModule : public CSolverOutputModule{

protected:
  su2double factor, Alpha, Beta, RefArea, RefLength, Gas_Constant, Gamma, Prandtl_Lam;
  bool axisymmetric;

  su2double Force[3], Moment[3], SkinFriction[3], Heatflux, Y_Plus;

public:

  explicit CAerodynamicsModule(CConfig *config);

  void LoadHistoryDataPerSurface(COutFieldCollection& fieldCollection) override;

  void DefineHistoryFields(COutFieldCollection& fieldCollection) override;

  void DefineVolumeFields(COutFieldCollection& fieldCollection) override;

  void LoadSurfaceData(COutFieldCollection& fieldCollection) override;

};

