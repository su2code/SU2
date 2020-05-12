#include "../../../include/iteration/modules/CFixedCLModule.hpp"
#include "../../../include/output/COutput.hpp"


void CFixedCLModule::PreIterationHook(COutput *output, CConfig *config, CGeometry *geometry, CSolver **solver){
  SetFarfield_AoA(geometry, solver, config, true);
}
void CFixedCLModule::PostIterationHook(COutput *output, CConfig *config, CGeometry *geometry, CSolver **solver){
  Monitor_CL_Convergence(output, config, geometry, solver);
}

void CFixedCLModule::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, bool Output) {

  su2double AoA = 0.0, Vel_Infty[3], Vel_Infty_Mag;
  const unsigned long InnerIter = config->GetInnerIter();
  const su2double Beta = config->GetAoS();

  const int iMesh = geometry->GetiMesh();

  const int nDim = geometry->GetnDim();
  /* --- Initialize values at first iteration --- */

  if (InnerIter == 0) {
    Total_CD_Prev = 0.0;
    Total_CL_Prev = 0.0;
    Total_CMx_Prev = 0.0;
    Total_CMy_Prev = 0.0;
    Total_CMz_Prev = 0.0;
    AoA_Prev = config->GetAoA();
    dCL_dAlpha = config->GetdCL_dAlpha();
    AoA_inc = 0.0;
    Iter_Update_AoA = 0;
    Start_AoA_FD = false;
    End_AoA_FD = false;
  }

  /*--- Retrieve the AoA (degrees) ---*/

  AoA = config->GetAoA();

  /* --- Set new AoA if needed --- */

  if (fabs(AoA_inc) > 0.0 && Output) {
    cout << "SETTING FARFIELD AOA" << endl;

    unsigned short iDim;

    /* --- Update *_Prev values with current coefficients --- */

    SetCoefficient_Gradients(config, solver_container);

    Total_CD_Prev = solver_container[FLOW_SOL]->GetTotal_CD();
    Total_CL_Prev = solver_container[FLOW_SOL]->GetTotal_CL();
    Total_CMx_Prev = solver_container[FLOW_SOL]->GetTotal_CMx();
    Total_CMy_Prev = solver_container[FLOW_SOL]->GetTotal_CMy();
    Total_CMz_Prev = solver_container[FLOW_SOL]->GetTotal_CMz();
    AoA_Prev = AoA;

    /*--- Compute a new value for AoA on the fine mesh only (degrees)---*/

    if (iMesh == MESH_0) AoA = AoA + AoA_inc;
    else { AoA = config->GetAoA(); }

    /*--- Only the fine mesh stores the updated values for AoA in config ---*/

    if (iMesh == MESH_0) {
      config->SetAoA(AoA);
    }

    /*--- Update the freestream velocity vector at the farfield ---*/

    for (iDim = 0; iDim < nDim; iDim++)
      Vel_Infty[iDim] = solver_container[FLOW_SOL]->GetVelocity_Inf(iDim);

    /*--- Compute the magnitude of the free stream velocity ---*/

    Vel_Infty_Mag = 0;
    for (iDim = 0; iDim < nDim; iDim++)
      Vel_Infty_Mag += Vel_Infty[iDim]*Vel_Infty[iDim];
    Vel_Infty_Mag = sqrt(Vel_Infty_Mag);

    /*--- Compute the new freestream velocity with the updated AoA ---*/

    if (nDim == 2) {
      Vel_Infty[0] = cos(AoA*PI_NUMBER/180.0)*Vel_Infty_Mag;
      Vel_Infty[1] = sin(AoA*PI_NUMBER/180.0)*Vel_Infty_Mag;
    }
    if (nDim == 3) {
      Vel_Infty[0] = cos(AoA*PI_NUMBER/180.0)*cos(Beta*PI_NUMBER/180.0)*Vel_Infty_Mag;
      Vel_Infty[1] = sin(Beta)*Vel_Infty_Mag;
      Vel_Infty[2] = sin(AoA*PI_NUMBER/180.0)*cos(Beta*PI_NUMBER/180.0)*Vel_Infty_Mag;
    }

    /*--- Store the new freestream velocity vector for the next iteration ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      solver_container[FLOW_SOL]->GetVelocity_Inf()[iDim] = Vel_Infty[iDim];
    }

    /*--- Only the fine mesh stores the updated values for velocity in config ---*/

    if (iMesh == MESH_0) {
      for (iDim = 0; iDim < nDim; iDim++)
        config->SetVelocity_FreeStreamND(Vel_Infty[iDim], iDim);
    }

  }
}

void CFixedCLModule::SetCoefficient_Gradients(CConfig *config, CSolver **solver){
  su2double dCL_dAlpha_, dCD_dCL_, dCMx_dCL_, dCMy_dCL_, dCMz_dCL_;
  su2double AoA = config->GetAoA();

  if (AoA != AoA_Prev) {
  /* --- Calculate gradients of coefficients w.r.t. CL --- */

  dCL_dAlpha_ = (solver[FLOW_SOL]->GetTotal_CL()-Total_CL_Prev)/(AoA - AoA_Prev);
  dCD_dCL_    = (solver[FLOW_SOL]->GetTotal_CD()-Total_CD_Prev)/(solver[FLOW_SOL]->GetTotal_CL()-Total_CL_Prev);
  dCMx_dCL_   = (solver[FLOW_SOL]->GetTotal_CMx()-Total_CMx_Prev)/(solver[FLOW_SOL]->GetTotal_CL()-Total_CL_Prev);
  dCMy_dCL_   = (solver[FLOW_SOL]->GetTotal_CMy()-Total_CMy_Prev)/(solver[FLOW_SOL]->GetTotal_CL()-Total_CL_Prev);
  dCMz_dCL_   = (solver[FLOW_SOL]->GetTotal_CMz()-Total_CMz_Prev)/(solver[FLOW_SOL]->GetTotal_CL()-Total_CL_Prev);

  /*--- Set the value of the  dOF/dCL in the config file ---*/

  config->SetdCD_dCL(dCD_dCL_);
  config->SetdCMx_dCL(dCMx_dCL_);
  config->SetdCMy_dCL(dCMy_dCL_);
  config->SetdCMz_dCL(dCMz_dCL_);
  config->SetdCL_dAlpha(dCL_dAlpha_);
  }
}

bool CFixedCLModule::FixedCL_Convergence(CConfig *config, CSolver **solver, bool convergence){
  su2double Target_CL = config->GetTarget_CL();
  unsigned long curr_iter = config->GetInnerIter();
  unsigned long Iter_dCL_dAlpha = config->GetIter_dCL_dAlpha();
  bool fixed_cl_conv = false;
  AoA_inc = 0.0;


  /*--- if in Fixed CL mode, before finite differencing --- */

  if (!Start_AoA_FD){
    if (convergence){

      /* --- C_L and solution are converged, start finite differencing --- */

      if (fabs(solver[FLOW_SOL]->GetTotal_CL()-Target_CL) < (config->GetCauchy_Eps()/2)) {

        /* --- If no finite differencing required --- */

        if (Iter_dCL_dAlpha == 0){
          fixed_cl_conv = true;
          return fixed_cl_conv;
        }

        /* --- Else, set up finite differencing routine ---*/

        Iter_Update_AoA = curr_iter;
        Start_AoA_FD = true;
        fixed_cl_conv = false;
        AoA_inc = 0.001;
      }

      /* --- C_L is not converged to target value and some iterations
          have passed since last update, so update AoA --- */

      else if ((curr_iter - Iter_Update_AoA) > config->GetStartConv_Iter()){
        Iter_Update_AoA = curr_iter;
        fixed_cl_conv = false;
        if (fabs(solver[FLOW_SOL]->GetTotal_CL() - Target_CL) > (config->GetCauchy_Eps()/2)) {
          AoA_inc = (1.0/dCL_dAlpha)*(Target_CL - solver[FLOW_SOL]->GetTotal_CL());
        }
      }
    }

    /* --- If the iteration limit between AoA updates is met, so update AoA --- */

    else if ((curr_iter - Iter_Update_AoA) == config->GetUpdate_AoA_Iter_Limit()) {
      Iter_Update_AoA = curr_iter;
      fixed_cl_conv = false;
      if (fabs(solver[FLOW_SOL]->GetTotal_CL()-Target_CL) > (config->GetCauchy_Eps()/2)) {
        AoA_inc = (1.0/dCL_dAlpha)*(Target_CL - solver[FLOW_SOL]->GetTotal_CL());
      }
    }

    /* --- If the total iteration limit is reached, start finite differencing --- */

    if (curr_iter == config->GetnInner_Iter() - Iter_dCL_dAlpha){
      if (Iter_dCL_dAlpha == 0){
        End_AoA_FD = true;
      }
      Iter_Update_AoA = curr_iter;
      Start_AoA_FD = true;
      fixed_cl_conv = false;
      AoA_inc = 0.001;
    }
  }

  /* --- If Finite Difference Mode has ended, end simulation --- */

  if (End_AoA_FD){
    //fixed_cl_conv = true;
    return true;
  }

  /* --- If starting Finite Difference Mode --- */

  if (Start_AoA_FD){

    /* --- Disable history writing --- */

    config->SetHistory_Wrt_Freq(2, 0);

    /* --- End Finite Difference Mode if iteration limit is reached, so simualtion is converged --- */

    End_AoA_FD = ((curr_iter - Iter_Update_AoA - 2) == Iter_dCL_dAlpha ||
      curr_iter == config->GetnInner_Iter()- 2 );

    if (convergence && (curr_iter - Iter_Update_AoA) > config->GetStartConv_Iter())
      End_AoA_FD = true;


    /* --- If Finite Difference mode is ending, reset AoA and calculate Coefficient Gradients --- */

    if (End_AoA_FD){
      SetCoefficient_Gradients(config, solver);
      config->SetAoA(AoA_Prev);
    }
  }

  return fixed_cl_conv;

}

void CFixedCLModule::Monitor_CL_Convergence(COutput *output, CConfig *config, CGeometry *geometry, CSolver **solver){

  bool fixed_cl_convergence = FixedCL_Convergence(config, solver, output->GetConvergence());

  /* --- If Fixed CL mode has ended and Finite Differencing has started: --- */

  if (Start_AoA_FD && Iter_Update_AoA == config->GetInnerIter()){

    /* --- Print convergence history and volume files since fixed CL mode has converged--- */
    if (SU2_MPI::GetRank() == MASTER_NODE) output->PrintConvergenceSummary();

    output->SetResult_Files(geometry, config, solver,
                            config->GetInnerIter(), true);

    /* --- Set finite difference mode in config (disables output) --- */
    config->SetFinite_Difference_Mode(true);
  }

  output->SetConvergence(fixed_cl_convergence);

}
