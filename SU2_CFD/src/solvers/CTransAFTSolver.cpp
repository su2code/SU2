/*!
 * \file CTransAFTSolver.cpp
 * \brief Main subroutines for Langtry-Menter Transition model solver.
 * \author S. Kang.
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/solvers/CTransAFTSolver.hpp"
#include "../../include/variables/CTransAFTVariable.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../include/variables/CTurbSAVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

/*---  This is the implementation of the Amplification Factor Transport(SA-AFT2019b) transition model.
       The main reference for this model is:Coder, In AIAA scitech 2019 forum.
       DOI: https://doi.org/10.2514/6.2019-0039 ---*/

// Note: TransAFT seems to use rho*AF, rho*ln(Intermittency) as Solution variables, thus Conservative=true

CTransAFTSolver::CTransAFTSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTurbSolver(geometry, config, true) {
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;

  bool multizone = config->GetMultizone_Problem();

  /*--- Dimension of the problem --> 2 Transport equations (AF, ln(intermittency)) ---*/
  nVar = 2;
  nPrimVar = 2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Define variables needed for transition from config file */
  options = config->GetAFTParsedOptions();
  TransCorrelations.SetOptions(options);

  /*--- Single grid simulation ---*/

  if (iMesh == MESH_0) {

    /*--- Define some auxiliary vector related with the residual ---*/

    Residual_RMS.resize(nVar,0.0);
    Residual_Max.resize(nVar,0.0);
    Point_Max.resize(nVar,0);
    Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (AFT transition model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    System.SetxIsZero(true);

    if (ReducerStrategy)
      EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS.resize(nVar,0.0);
      Residual_Max_BGS.resize(nVar,0.0);
      Point_Max_BGS.resize(nVar,0);
      Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
    }

  }

  /*--- Initialize lower and upper limits---*/
  lowerlimit[0] = 0.0;
  upperlimit[0] = 100.0;

  lowerlimit[1] = -20.0;
  upperlimit[1] = 1.0e-12;

  /*--- Far-field flow state quantities and initialization. ---*/
  const su2double AF_Inf = 0.0;
  const su2double lnIntermittency_Inf = 0.0;

  Solution_Inf[0] = AF_Inf;
  Solution_Inf[1] = lnIntermittency_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/
  nodes = new CTransAFTVariable(AF_Inf, lnIntermittency_Inf, nPoint, nDim, nVar, config);  
  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
      SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
    }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
    due to arbitrary number of turbulence variables ---*/

  Inlet_TurbVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_TurbVars[iMarker].resize(nVertex[iMarker],nVar);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
      Inlet_TurbVars[iMarker](iVertex,0) = AF_Inf;
      Inlet_TurbVars[iMarker](iVertex,1) = lnIntermittency_Inf;
    }
  }

   const su2double CFL = config->GetCFL(MGLevel)*config->GetCFLRedCoeff_Turb();
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
    nodes->SetAuxVar(iPoint, 0, geometry->nodes->GetWall_Distance(iPoint));
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name. ---*/
  SolverName = "AFT model";

}

void CTransAFTSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
         unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)

  /*--- Upwind second order reconstruction and gradients ---*/
  CommonPreprocessing(geometry, config, Output);
}

void CTransAFTSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

  /*--- Compute LM model gradients. ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetSolution_Gradient_GG(geometry, config, -1);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetSolution_Gradient_LS(geometry, config, -1);
  }

  AD::StartNoSharedReading();
  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());
  auto* turbNodes = su2staticcast_p<CTurbVariable*>(solver_container[TURB_SOL]->GetNodes());

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    // Here the nodes already have the new solution, thus I have to compute everything from scratch
    const su2double AF = nodes->GetSolution(iPoint,0);
    const su2double lnIntermittency = nodes->GetSolution(iPoint,1);
    const su2double Critical_N_Factor = config->GetN_Critical();
    /* Ampification Factor term */
    const su2double vel_u = flowNodes->GetVelocity(iPoint, 0);
    const su2double vel_v = flowNodes->GetVelocity(iPoint, 1);
    const su2double vel_w = (nDim == 3) ? flowNodes->GetVelocity(iPoint, 2) : 0.0;
    const su2double Velocity_Mag = sqrt(vel_u * vel_u + vel_v * vel_v + vel_w * vel_w);
    const su2double dist_i = geometry->nodes->GetWall_Distance(iPoint);
    const su2double StrainMag_i = flowNodes->GetStrainMag(iPoint);
    const su2double Eddy_Viscosity_i = turbNodes->GetmuT(iPoint);
    const su2double Laminar_Viscosity_i = flowNodes->GetLaminarViscosity(iPoint);
    const su2double Density_i = flowNodes->GetDensity(iPoint);
    const su2double c_1 = 100.0;
    const su2double c_2 = 0.06;
    const su2double c_3 = 50.0;
    su2double Temp1 = nodes->GetAuxVarGradient(iPoint, 0, 0);
    su2double Temp2 = nodes->GetAuxVarGradient(iPoint, 0, 1);
    su2double Temp3 = flowNodes->GetDensity(iPoint) * flowNodes->GetVelocity(iPoint, 0) * nodes->GetAuxVarGradient(iPoint, 0, 0);
    su2double Temp4 = flowNodes->GetDensity(iPoint) * flowNodes->GetVelocity(iPoint, 1) * nodes->GetAuxVarGradient(iPoint, 0, 1);

    if(nDim == 2) {
      nodes->SetAuxVar(iPoint, 1, Temp3 + Temp4);
    }
    else {
      nodes->SetAuxVar(iPoint, 1, Temp3 + Temp4 + flowNodes->GetDensity(iPoint) * flowNodes->GetVelocity(iPoint, 2) * nodes->GetAuxVarGradient(iPoint, 0, 2));
    }

    su2double VorticityMag = 0.0;
    su2double HLGradTerm = 0.0;
    HLGradTerm = nodes->GetAuxVarGradient(iPoint, 1, 0) * nodes->GetAuxVarGradient(iPoint, 0, 0) + nodes->GetAuxVarGradient(iPoint, 1, 1) * nodes->GetAuxVarGradient(iPoint, 0, 1);

    if(nDim == 2) {
      VorticityMag = sqrt(flowNodes->GetVorticity(iPoint)[0] * flowNodes->GetVorticity(iPoint)[0] + flowNodes->GetVorticity(iPoint)[1] * flowNodes->GetVorticity(iPoint)[1] );
      
    }
    else {
      VorticityMag = sqrt(flowNodes->GetVorticity(iPoint)[0] * flowNodes->GetVorticity(iPoint)[0] + flowNodes->GetVorticity(iPoint)[1] * flowNodes->GetVorticity(iPoint)[1]
                    + flowNodes->GetVorticity(iPoint)[2] * flowNodes->GetVorticity(iPoint)[2] );
      HLGradTerm += nodes->GetAuxVarGradient(iPoint, 1, 2) * nodes->GetAuxVarGradient(iPoint, 0, 2);
    }

    /*--- Cal H12, Hk, dNdRet, Ret0 ---*/
    const su2double HL = dist_i * dist_i / Laminar_Viscosity_i * HLGradTerm;
    const su2double H12 = TransCorrelations.H12_Correlations(HL);
    const su2double dNdRet = TransCorrelations.dNdRet_Correlations(H12);
    const su2double Ret0 = TransCorrelations.Ret0_Correlations(H12);
    const su2double D_H12 = TransCorrelations.D_Correlations(H12);
    const su2double l_H12 = TransCorrelations.l_Correlations(H12);
    const su2double m_H12 = TransCorrelations.m_Correlations(H12, l_H12);
    const su2double kv = TransCorrelations.kv_Correlations(H12);
    const su2double Rev = Density_i * dist_i * dist_i * StrainMag_i / (Laminar_Viscosity_i + Eddy_Viscosity_i);
    const su2double Rev0 = kv * Ret0;

    /*--- production term of the amplification factor ---*/
    const su2double F_growth = max(D_H12 * (1.0 + m_H12) / 2.0 * l_H12, 0.0);
    const su2double F_crit = (Rev < Rev0) ? 0.0 : 1.0;
    const su2double PAF = Density_i * VorticityMag * F_crit * F_growth * dNdRet;

    /*--- production term of the intermittency(Gamma) ---*/
    const su2double RT = Eddy_Viscosity_i / Laminar_Viscosity_i;
    const su2double F_onset1 = AF / Critical_N_Factor;
    const su2double F_onset2 = min(F_onset1, 2.0);
    const su2double F_onset3 = max(1.0 - pow(RT / 3.5, 3), 0.0);
    const su2double F_onset = max(F_onset2 - F_onset3, 0.0);
    const su2double F_turb = exp(-pow(RT / 2.0,4));
    const su2double Pg = c_1 * Density_i * StrainMag_i * F_onset * (1.0 - exp(lnIntermittency));

    /*-- destruction term of Intermeittency(Gamma) --*/
    const su2double Dg = c_2 * Density_i * VorticityMag * F_turb * (c_3 * exp(lnIntermittency) - 1.0);
    nodes -> SetAFT_Wonder_Func(iPoint, HL, H12, dNdRet, Ret0, D_H12, l_H12, m_H12, kv, Rev, Rev0, F_growth, F_crit, PAF, Pg, Dg);
  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

  switch (config->GetKind_Gradient_Method()) {
    case GREEN_GAUSS:
      SetAuxVar_Gradient_GG(geometry, config);
      break;
    case WEIGHTED_LEAST_SQUARES:
      SetAuxVar_Gradient_LS(geometry, config);
    default:
      break;
    }
}


void CTransAFTSolver::Viscous_Residual(const unsigned long iEdge, const CGeometry* geometry, CSolver** solver_container,
                                     CNumerics* numerics, const CConfig* config) {

  /*--- Define an object to set solver specific numerics contribution. ---*/

  auto SolverSpecificNumerics = [&](unsigned long iPoint, unsigned long jPoint) {};

  /*--- Now instantiate the generic implementation with the functor above. ---*/

  Viscous_Residual_impl(SolverSpecificNumerics, iEdge, geometry, solver_container, numerics, config);
}


void CTransAFTSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());
  //auto* turbNodes = su2staticcast_p<CFlowVariable*>(solver_container[TURB_SOL]->GetNodes());
  CVariable* turbNodes = solver_container[TURB_SOL]->GetNodes();

  /*--- Pick one numerics object per thread. ---*/
  auto* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Loop over all points. ---*/

  AD::StartNoSharedReading();

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {



    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint), nullptr);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
    /*--- ScalarVar & ScalarVarGradient : Turbulence model solution(k&w) ---*/

    numerics->SetScalarVar(turbNodes->GetSolution(iPoint), nullptr);
    numerics->SetScalarVarGradient(turbNodes->GetGradient(iPoint), nullptr);

    /*--- Transition variables w/o reconstruction, and its gradient ---*/

    numerics->SetTransVar(nodes->GetSolution(iPoint), nullptr);
    numerics->SetTransVarGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Set distance to the surface ---*/

    numerics->SetDistance(geometry->nodes->GetWall_Distance(iPoint), 0.0);

    /*--- Set vorticity and strain rate magnitude ---*/

    numerics->SetVorticity(flowNodes->GetVorticity(iPoint), nullptr);

    numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint), 0.0);

    /*--- Set coordinate (for debugging) ---*/
    numerics->SetCoord(geometry->nodes->GetCoord(iPoint), nullptr);
    numerics->SetAuxVarGrad(nodes->GetAuxVarGradient(iPoint), nullptr);

    su2double Temp1 = flowNodes->GetDensity(iPoint) * flowNodes->GetVelocity(iPoint, 0) * nodes->GetAuxVarGradient(iPoint, 0, 0);
    su2double Temp2 = flowNodes->GetDensity(iPoint) * flowNodes->GetVelocity(iPoint, 1) * nodes->GetAuxVarGradient(iPoint, 0, 1);
    su2double Temp3 = nodes->GetAuxVarGradient(iPoint, 0, 0);
    su2double Temp4 = nodes->GetAuxVarGradient(iPoint, 0, 1);
    su2double Temp5 = nodes->GetAuxVarGradient(iPoint, 1, 0);
    su2double Temp6 = nodes->GetAuxVarGradient(iPoint, 1, 1);
    

    /*--- Compute the source term ---*/
    auto residual = numerics->ComputeResidual(config);

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, residual);
    if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

}

void CTransAFTSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                       CConfig *config, unsigned short iMesh) {
}

void CTransAFTSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Allocate the value at the infinity ---*/

      auto V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at the wall, and at infinity ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Solution_Inf);

      /*--- Set Normal (it is necessary to change the sign) ---*/
      /*--- It's mean wall normal zero flux. */

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Grid Movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute residuals and Jacobians ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Add residuals and Jacobians ---*/

      LinSysRes.AddBlock(iPoint, residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
    }
  }
  END_SU2_OMP_FOR

}

void CTransAFTSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}

void CTransAFTSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                unsigned short val_marker) {
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/

      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Non-dimensionalize Inlet_TurbVars if Inlet-Files are used. ---*/
      su2double Inlet_Vars[MAXNVAR];
      Inlet_Vars[0] = Inlet_TurbVars[val_marker][iVertex][0];
      Inlet_Vars[1] = Inlet_TurbVars[val_marker][iVertex][1];
      if (config->GetInlet_Profile_From_File()) {
        Inlet_Vars[0] /= pow(config->GetVelocity_Ref(), 2);
        Inlet_Vars[1] *= config->GetViscosity_Ref() / (config->GetDensity_Ref() * pow(config->GetVelocity_Ref(), 2));
      }

      /*--- Set the LM variable states. ---*/
      /*--- Load the inlet transition LM model variables (uniform by default). ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Inlet_Vars);

      /*--- Set various other quantities in the solver class ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
    }
  }
  END_SU2_OMP_FOR

}

void CTransAFTSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                               CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  BC_Far_Field(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CTransAFTSolver::LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter,
                                  bool val_update_geo) {

  const string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

    if (config->GetRead_Binary_Restart()) {
      Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
    } else {
      Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
    }

    /*--- Skip flow variables ---*/

    unsigned short skipVars = nDim + solver[MESH_0][FLOW_SOL]->GetnVar() + solver[MESH_0][TURB_SOL] ->GetnVar();

    /*--- Adjust the number of solution variables in the incompressible
     restart. We always carry a space in nVar for the energy equation in the
     mean flow solver, but we only write it to the restart if it is active.
     Therefore, we must reduce skipVars here if energy is inactive so that
     the turbulent variables are read correctly. ---*/

    const bool incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
    const bool energy = config->GetEnergy_Equation();
    const bool weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

    if (incompressible && ((!energy) && (!weakly_coupled_heat))) skipVars--;

    /*--- Load data from the restart into correct containers. ---*/

    unsigned long counter = 0;
    for (auto iPoint_Global = 0ul; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      const auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {
        /*--- We need to store this point's data, so jump to the correct
         offset in the buffer of data from the restart file and load it. ---*/

        const auto index = counter * Restart_Vars[1] + skipVars;
        for (auto iVar = 0u; iVar < nVar; iVar++) nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index + iVar]);
        nodes ->SetAFT_Wonder_Func(iPoint_Local, Restart_Data[index + 2], Restart_Data[index + 3]
              , Restart_Data[index + 4], Restart_Data[index + 5], Restart_Data[index + 6], Restart_Data[index + 7]
              , Restart_Data[index + 8], Restart_Data[index + 9], Restart_Data[index + 10], Restart_Data[index + 11]
              , Restart_Data[index + 12], Restart_Data[index + 13], Restart_Data[index + 14], Restart_Data[index + 15]
              , Restart_Data[index + 16]);

        /*--- Increment the overall counter for how many points have been loaded. ---*/
        counter++;
      }
    }

    /*--- Detect a wrong solution file ---*/

    if (counter != nPointDomain) {
      SU2_MPI::Error(string("The solution file ") + restart_filename + string(" does not match with the mesh file!\n") +
                         string("This can be caused by empty lines at the end of the file."),
                     CURRENT_FUNCTION);
    }

  }  // end SU2_OMP_MASTER, pre and postprocessing are thread-safe.
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- MPI solution and compute the eddy viscosity ---*/

  solver[MESH_0][TRANS_SOL]->InitiateComms(geometry[MESH_0], config, MPI_QUANTITIES::SOLUTION);
  solver[MESH_0][TRANS_SOL]->CompleteComms(geometry[MESH_0], config, MPI_QUANTITIES::SOLUTION);

  /*--- For turbulent+species simulations the solver Pre-/Postprocessing is done by the species solver. ---*/
  if (config->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
    solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER,
                                            RUNTIME_FLOW_SYS, false);
    solver[MESH_0][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
    solver[MESH_0][TRANS_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
  }

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (auto iMesh = 1u; iMesh <= config->GetnMGLevels(); iMesh++) {

    MultigridRestriction(*geometry[iMesh - 1], solver[iMesh - 1][TRANS_SOL]->GetNodes()->GetSolution(),
                         *geometry[iMesh], solver[iMesh][TRANS_SOL]->GetNodes()->GetSolution());
    solver[iMesh][TRANS_SOL]->InitiateComms(geometry[iMesh], config, MPI_QUANTITIES::SOLUTION);
    solver[iMesh][TRANS_SOL]->CompleteComms(geometry[iMesh], config, MPI_QUANTITIES::SOLUTION);

    if (config->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
      solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS,
                                            false);
      solver[iMesh][TRANS_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
    }
  }

  /*--- Go back to single threaded execution. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Delete the class memory that is used to load the restart. ---*/

    Restart_Vars.clear();
    Restart_Data.clear();
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

}