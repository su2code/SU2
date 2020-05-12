#pragma once
#include "CSolverModule.hpp"
#include "../../../include/solvers/CEulerSolver.hpp"

class CFixedCLModule : public CSolverModule{

  su2double Total_CD_Prev = 0.0;
  su2double Total_CL_Prev = 0.0;
  su2double Total_CMx_Prev = 0.0;
  su2double Total_CMy_Prev = 0.0;
  su2double Total_CMz_Prev = 0.0;
  su2double AoA_Prev;
  su2double dCL_dAlpha;
  su2double AoA_inc;

  bool Start_AoA_FD, End_AoA_FD;
  unsigned long Iter_Update_AoA;


public:

  explicit CFixedCLModule(CConfig *config) : CSolverModule(config->GetFixed_CL_Mode()) {};

  void PreIterationHook(COutput *output, CConfig *config, CGeometry *geometry, CSolver **solver) override;

  void PostIterationHook(COutput *output, CConfig *config, CGeometry *geometry, CSolver **solver) override;

  void SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                       CConfig *config, bool Output);

  void SetCoefficient_Gradients(CConfig *config, CSolver **solver);

  bool FixedCL_Convergence(CConfig* config, CSolver **solver, bool convergence);

  void Monitor_CL_Convergence(COutput *output, CConfig *config, CGeometry *geometry, CSolver **solver);
};
