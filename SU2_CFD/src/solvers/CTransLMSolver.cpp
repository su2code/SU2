/*!
 * \file CTransLMSolver.cpp
 * \brief Main subroutines for Langtry-Menter Transition model solver.
 * \author A. Aranake, A. Rausa, M. Cerabona
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/solvers/CTransLMSolver.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

template class CScalarSolver<CTransLMVariable>;

/*---  This is the implementation of the Langtry-Menter transition model.
       The main reference for this model is: https://turbmodels.larc.nasa.gov/langtrymenter_4eqn.html ---*/


CTransLMSolver::CTransLMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTransSolver(geometry, config, true) {

  unsigned short nLineLets;
  unsigned long iPoint;

  bool multizone = config->GetMultizone_Problem();

  /*--- Dimension of the problem --> 2 Transport equations (intermittency, Reth) ---*/
  nVar = 2;
  nPrimVar = 2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Single grid simulation ---*/

  if (iMesh == MESH_0) {

    /*--- Define some auxiliary vector related with the residual ---*/

    Residual_RMS.resize(nVar,0.0);
    Residual_Max.resize(nVar,0.0);
    Point_Max.resize(nVar,0);
    Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (LM transition model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

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

  /*--- Initialize value for model constants ---*/
  constants[0] = 2.0;    //c_a1
  constants[1] = 0.06;   //c_a2
  constants[2] = 1.0;    //c_e1
  constants[3] = 50.0;   //c_e2
  constants[4] = 0.03;   //c_thetat
  constants[5] = 2.0;   // s_1
  constants[6] = 1.0;    //sigma_f
  constants[7] = 2.0;    //sigma_thetat
  constants[8] = 0.6;    //c_CF


  /*--- Initialize lower and upper limits---*/
  lowerlimit[0] = 1.0e-4;
  upperlimit[0] = 5.0;

  lowerlimit[1] = 1.0e-4;
  upperlimit[1] = 1.0e15;

  /*--- Far-field flow state quantities and initialization. ---*/
  su2double TU_inf;

  TU_inf = config->GetTurbulenceIntensity_FreeStream()*100;
  TU_inf = max(TU_inf, 0.027); // Limit the turbulence intensity to 0.027

  su2double intermittency_Inf = 1.0;
  su2double re_thetat_Inf = 100.0;


  /*--- Momentum thickness Reynolds number, initialized from freestream turbulent intensity ---*/
  if (TU_inf <= 1.3) {
      re_thetat_Inf = (1173.51-589.428*TU_inf+0.2196/(TU_inf*TU_inf));
  }
  else {
    re_thetat_Inf = 331.5*pow(TU_inf-0.5658,-0.671);
  }


  Solution_Inf[0] = intermittency_Inf;
  Solution_Inf[1] = re_thetat_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/
  nodes = new CTransLMVariable(intermittency_Inf, re_thetat_Inf, constants, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
          SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
          SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
      }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTransSolver),
  due to arbitrary number of transition variables ---*/

  Inlet_TransVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
      Inlet_TransVars[iMarker].resize(nVertex[iMarker],nVar);
      for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
          Inlet_TransVars[iMarker](iVertex,0) = intermittency_Inf;
          Inlet_TransVars[iMarker](iVertex,1) = re_thetat_Inf;
      }
  }

  /*--- Store the initial CFL number for all grid points. ---*/

  const su2double CFL = config->GetCFL(MGLevel)*config->GetCFLRedCoeff_Trans();
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
      nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "LM model";   

}

void CTransLMSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                   unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {


  config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);

  /*--- Upwind second order reconstruction and gradients ---*/
  CommonPreprocessing(geometry, config, Output);
}

void CTransLMSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container,
                                    CConfig *config, unsigned short iMesh) {

    /*--- Compute transition gradients. ---*/

    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
        SetSolution_Gradient_GG(geometry, config);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
        SetSolution_Gradient_LS(geometry, config);
    }

    AD::StartNoSharedReading();


    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

        nodes->SetGammaEff(iPoint);

    }
    END_SU2_OMP_FOR

    AD::EndNoSharedReading();
}

void CTransLMSolver::Viscous_Residual(unsigned long iEdge, CGeometry* geometry, CSolver** solver_container,
                                      CNumerics* numerics, CConfig* config) {


    /*--- Define an object to set solver specific numerics contribution. ---*/
    auto SolverSpecificNumerics = [&](unsigned long iPoint, unsigned long jPoint) { };

    /*--- Now instantiate the generic implementation with the functor above. ---*/

    Viscous_Residual_impl(SolverSpecificNumerics, iEdge, geometry, solver_container, numerics, config);
    
}

void CTransLMSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {


    bool axisymmetric = config->GetAxisymmetric();

    const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

    auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());
    auto* turbNodes = su2staticcast_p<CTurbVariable*>(solver_container[TURB_SOL]->GetNodes());

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

        /*--- Transition variables w/o reconstruction, and its gradient ---*/

        numerics->SetScalarVar(nodes->GetSolution(iPoint), nullptr);

        numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nullptr);

        /*--- Set volume ---*/

        numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

        /*--- Set distance to the surface ---*/

        numerics->SetDistance(geometry->nodes->GetWall_Distance(iPoint), 0.0);

        if(geometry->nodes->GetWall_Distance(iPoint) > 1e-10) {

          su2double VorticityMag = GeometryToolbox::Norm(3, flowNodes->GetVorticity(iPoint));
          VorticityMag = max(VorticityMag, 1e-12);  

          // Aggiunto da me, sono da settare
          su2double Velocity[nDim];
          Velocity[0] = flowNodes->GetVelocity(iPoint, 0);
          Velocity[1] = flowNodes->GetVelocity(iPoint, 1);
          if (nDim == 3) Velocity[2] = flowNodes->GetVelocity(iPoint, 2);

          su2double VelocityMag = GeometryToolbox::Norm(3, Velocity);
          VelocityMag = max(VelocityMag, 1e-12); 


          if(TurbModelFamily(config->GetKind_Turb_Model()) == TURB_FAMILY::SA && config->GetConvertSA2SST())
            nodes->Set_kAndw(iPoint, VorticityMag, VelocityMag, turbNodes->GetmuT(iPoint),  flowNodes->GetStrainMag(iPoint),  geometry->nodes->GetWall_Distance(iPoint),
                             flowNodes->GetLaminarViscosity(iPoint),  config);


          nodes->SetF_onset(iPoint, flowNodes->GetDensity(iPoint), flowNodes->GetStrainMag(iPoint),
                            geometry->nodes->GetWall_Distance(iPoint), flowNodes->GetLaminarViscosity(iPoint),
                            turbNodes->GetSolution(iPoint), turbNodes->GetmuT(iPoint), VorticityMag, VelocityMag, config);

          nodes->SetF_length(Velocity, iPoint, flowNodes->GetDensity(iPoint), turbNodes->GetSolution(iPoint),
                             geometry->nodes->GetWall_Distance(iPoint), flowNodes->GetLaminarViscosity(iPoint), config);
          nodes->SetF_turb(iPoint);
          nodes->Setrethetat_eq(iPoint, Velocity, VelocityMag, flowNodes->GetVelocityGradient(iPoint), turbNodes->GetSolution(iPoint),
                                flowNodes->GetLaminarViscosity(iPoint), flowNodes->GetDensity(iPoint), config);
          nodes->SetT(iPoint, VelocityMag, flowNodes->GetLaminarViscosity(iPoint), flowNodes->GetDensity(iPoint), geometry->nodes->GetMaxLength(iPoint), turbNodes->GetmuT(iPoint));
          nodes->SetF_thetat(iPoint, flowNodes->GetDensity(iPoint), turbNodes->GetSolution(iPoint),
                             geometry->nodes->GetWall_Distance(iPoint), flowNodes->GetLaminarViscosity(iPoint),
                             VorticityMag, VelocityMag, constants, config);

          if (config->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM2015) {
            nodes->SetF_thetat_2(iPoint, geometry->nodes->GetWall_Distance(iPoint));
            nodes->SetReThetat_SCF(iPoint, geometry->nodes->GetWall_Distance(iPoint), flowNodes->GetDensity(iPoint),
                                   Velocity, VelocityMag, config->GethRoughness(), flowNodes->GetVorticity(iPoint),
                                   flowNodes->GetLaminarViscosity(iPoint), turbNodes->GetmuT(iPoint));
          }
          nodes->Setgamma_sep(iPoint, flowNodes->GetDensity(iPoint), turbNodes->GetSolution(iPoint),
                              flowNodes->GetLaminarViscosity(iPoint), geometry->nodes->GetWall_Distance(iPoint),
                              flowNodes->GetStrainMag(iPoint), VorticityMag, constants);

        }



        /*--- Quantities for the LM model contributions ---*/


        /*--- Set quantities useful for the source residual and jacobian computations ---*/
        numerics->SetF_length(nodes->GetF_length(iPoint), 0.0);
        numerics->SetF_onset(nodes->GetF_onset(iPoint), 0.0);
        numerics->SetF_thetat(nodes->GetF_thetat(iPoint), 0.0);
        numerics->SetF_turb(nodes->GetF_turb(iPoint), 0.0);
        numerics->SetF_wake(nodes->GetF_wake(iPoint), 0.0);
        numerics->Setrethetat_eq(nodes->Getrethetat_eq(iPoint), 0.0);
        numerics->SetT_param(nodes->GetT(iPoint), 0.0);
        numerics->Setdelta_param(nodes->Getdelta_param(iPoint), 0.0);

        if (config->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM2015) {
          numerics->SetReThetat_SCF(nodes->GetReThetat_SCF(iPoint), 0.0);
          numerics->SetF_thetat_2(nodes->GetF_thetat_2(iPoint), 0.0);
        }


        /*--- Set vorticity and strain rate magnitude ---*/

        numerics->SetVorticity(flowNodes->GetVorticity(iPoint), nullptr);

        numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint), 0.0);


        if (axisymmetric){
            /*--- Set y coordinate ---*/
            numerics->SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(iPoint));
        }

        /*--- Compute the source term ---*/
        auto residual = numerics->ComputeResidual(config);

        /*--- Subtract residual and the Jacobian ---*/

        LinSysRes.SubtractBlock(iPoint, residual);


        if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

    }
    END_SU2_OMP_FOR


    AD::EndNoSharedReading();

}

void CTransLMSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {
}

// Modifica del 16-04-22 alle 19:47
void CTransLMSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      conv_numerics->SetPrimitive(V_domain, V_domain);

      /*--- Set turbulent variable at the wall, and at infinity ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));

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

// Modifica del 16-04-22 alle 19:47 : Sostituito stessa function scritta sopra con questa.
// Sembra non cambiare niente
//void CTransLMSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
//                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
//
//  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
//
//  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
//  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//
//    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
//
//    if (geometry->nodes->GetDomain(iPoint)) {
//
//      /*--- Retrieve solution at the farfield boundary node ---*/
//
//      /*--- distance to closest neighbor ---*/
//      const auto jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
//
//      /*--- Set the solution values and zero the residual ---*/
//      nodes->SetSolution_Old(iPoint,nodes->GetSolution(jPoint));
//      nodes->SetSolution(iPoint,nodes->GetSolution(jPoint));
//      LinSysRes.SetBlock_Zero(iPoint);
//
//      if (implicit) {
//        /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
//        Jacobian.DeleteValsRowi(iPoint*nVar);
//        Jacobian.DeleteValsRowi(iPoint*nVar+1);
//      }
//    }
//  }
//  END_SU2_OMP_FOR
//
//}

void CTransLMSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}



void CTransLMSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

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

void CTransLMSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
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

      /*--- Non-dimensionalize Inlet_TransVars if Inlet-Files are used. ---*/
      su2double Inlet_Vars[MAXNVAR];
      Inlet_Vars[0] = Inlet_TransVars[val_marker][iVertex][0];
      Inlet_Vars[1] = Inlet_TransVars[val_marker][iVertex][1];
      // Non-dimensionalization should not be used since they are no-dimensional quantities
//      if (config->GetInlet_Profile_From_File()) {
//        Inlet_Vars[0] /= pow(config->GetVelocity_Ref(), 2);
//        Inlet_Vars[1] *= config->GetViscosity_Ref() / (config->GetDensity_Ref() * pow(config->GetVelocity_Ref(), 2));
//      }

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

void CTransLMSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                               CConfig *config, unsigned short val_marker) {

  BC_Far_Field(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}




