/*!
 * \file CPassiveScalarSolver.cpp
 * \brief Main subroutines for solving transported scalar class
 * \author T. Economon, N. Beishuizen
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/solvers/CPassiveScalarSolver.hpp"
#include "../../include/variables/CPassiveScalarVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

CPassiveScalarSolver::CPassiveScalarSolver(void) : CScalarLegacySolver() {}

CPassiveScalarSolver::CPassiveScalarSolver(CGeometry *geometry,
                                           CConfig *config,
                                           unsigned short iMesh)
: CScalarLegacySolver(geometry, config) {
  unsigned short nLineLets;

  const bool turbulent = ((config->GetKind_Solver() == RANS) ||
                         (config->GetKind_Solver() == DISC_ADJ_RANS));
  const bool turb_SST  = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  const bool turb_SA   = ((turbulent) && (config->GetKind_Turb_Model() == SA));
  bool multizone = config->GetMultizone_Problem();

  /*--- Dimension of the problem --> passive scalar will only ever
   have a single equation. Other child classes of CScalarLegacySolver
   can have variable numbers of equations. ---*/
 
  nVar     = 1;
  nPrimVar = 1;

  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();
  
  /*--- Fluid model pointer initialization ---*/
  
  FluidModel = nullptr;
  

  /*--- Single grid simulation ---*/

  /*--- Define some auxiliary vector related with the solution ---*/
  Solution = new su2double[nVar];
  Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];

  /*--- do not see a reason to use only single grid ---*/
  //if (iMesh == MESH_0) {

    /*--- Define some auxiliary vector related with the residual ---*/

    Residual_RMS.resize(nVar,0.0);
    Residual_Max.resize(nVar,0.0);
    Point_Max.resize(nVar,0);
    Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (passive scalar model)." << endl;
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

  //} //iMESH_0

  /*--- Initialize lower and upper limits---*/

  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  if (config->GetScalar_Clipping()){
    for (auto iVar=0u;iVar<nVar;iVar++){
      lowerlimit[iVar] = config->GetScalar_Clipping_Min(iVar);
      upperlimit[iVar] = config->GetScalar_Clipping_Max(iVar);
    }
  }
  else {
    for (auto iVar=0u;iVar<nVar;iVar++){
      lowerlimit[iVar] = -1.0e15;
      upperlimit[iVar] =  1.0e15;
    }
  }
  /*--- Far-field flow state quantities and initialization. ---*/
  //su2double Density_Inf, Viscosity_Inf;
  //Density_Inf   = config->GetDensity_FreeStreamND();
  //Viscosity_Inf = config->GetViscosity_FreeStreamND();

  /*--- Set up fluid model for the diffusivity ---*/

  su2double Diffusivity_Ref = 1.0;
  su2double DiffusivityND = config->GetDiffusivity_Constant()/Diffusivity_Ref;
  config->SetDiffusivity_ConstantND(DiffusivityND);

  FluidModel = new CFluidModel();
  FluidModel->SetMassDiffusivityModel(config);

  /*--- Scalar variable state at the far-field. ---*/

  Scalar_Inf = new su2double[nVar];
  for (auto iVar = 0u; iVar < nVar; iVar++){
    Scalar_Inf[iVar] = config->GetScalar_Init(iVar);
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CPassiveScalarVariable(Scalar_Inf, nPoint, nDim, nVar, config);


  /*--- initialize the mass diffusivity ---*/
  for (auto iVar = 0u; iVar < nVar; iVar++){
    auto M = FluidModel->GetMassDiffusivity(); // returns a su2double, note that for more species this should be a vector
    // loop over all points and set diffusivity
    // why construct the entire diffusivity matrix?
   for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
     nodes->SetDiffusivity(iPoint, M, iVar);
  }

  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Initialize quantities for SlidingMesh Interface ---*/

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
      SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
    }
  }


  /*-- Allocation of inlets has to happen in derived classes
   (not CScalarLegacySolver), due to arbitrary number of scalar variables.
   First, we also set the column index for any inlet profiles. ---*/
  
  Inlet_Position = nDim*2+2;
  if (turbulent) {
    if (turb_SA) Inlet_Position += 1;
    else if (turb_SST) Inlet_Position += 2;
  }

 /*-- Allocation of inlets has to happen in derived classes
   (not CScalarLegacySolver), due to arbitrary number of scalar variables.
   First, we also set the column index for any inlet profiles. ---*/
  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
   * due to arbitrary number of turbulence variables ---*/

  Inlet_ScalarVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++){
    Inlet_ScalarVars[iMarker].resize(nVertex[iMarker],nVar);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++){
        Inlet_ScalarVars[iMarker](iVertex,iVar) = Scalar_Inf[iVar];
      }
    }
  }

  /*--- The turbulence models are always solved implicitly, so set the
  implicit flag in case we have periodic BCs. ---*/

  SetImplicitPeriodic(true);

  /*--- Store the initial CFL number for all grid points. ---*/
 
  const su2double CFL = config->GetCFL(MGLevel)*config->GetCFLRedCoeff_Scalar();
  for (auto iPoint = 0u; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "SCALAR";

}

CPassiveScalarSolver::~CPassiveScalarSolver(void) {
  if (FluidModel != nullptr) delete FluidModel;
}


void CPassiveScalarSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
         unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool muscl = config->GetMUSCL_Scalar();
  const bool limiter = (config->GetKind_SlopeLimit_Scalar() != NO_LIMITER) &&
                       (config->GetInnerIter() <= config->GetLimiterIter());

  /*--- Clear residual and system matrix, not needed for
   * reducer strategy as we write over the entire matrix. ---*/
  if (!ReducerStrategy && !Output) {
    LinSysRes.SetValZero();
    if (implicit) Jacobian.SetValZero();
    else {SU2_OMP_BARRIER}
  }

  /*--- Upwind second order reconstruction and gradients ---*/

  if (config->GetReconstructionGradientRequired()) {
    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS)
      SetSolution_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES)
      SetSolution_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES)
      SetSolution_Gradient_LS(geometry, config, true);
  }

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
    SetSolution_Gradient_GG(geometry, config);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
    SetSolution_Gradient_LS(geometry, config);

  if (limiter && muscl) SetSolution_Limiter(geometry, config);

}

void CPassiveScalarSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container,
                                    CConfig *config, unsigned short iMesh) {

/*--- your postprocessing goes here ---*/

}

void CPassiveScalarSolver::SetInitialCondition(CGeometry **geometry,
                                               CSolver ***solver_container,
                                               CConfig *config,
                                               unsigned long ExtIter) {
  bool Restart   = (config->GetRestart() || config->GetRestart_Flow());
  
  
  if ((!Restart) && ExtIter == 0) {
    if (rank == MASTER_NODE){
      cout << "Initializing passive scalar (initial condition)." << endl;
      cout << "initialization = " << nVar << " " << config->GetScalar_Init(0)<<endl;
    }  

    su2double* scalar_init = new su2double[nVar];

    for (unsigned long i_mesh = 0; i_mesh <= config->GetnMGLevels(); i_mesh++) {

      for (unsigned long i_point = 0; i_point < geometry[i_mesh]->GetnPoint(); i_point++) {
        
        for (unsigned long i_var = 0; i_var < nVar; i_var++)
          Solution[i_var] = 0.0;
        
        for(int i_scalar = 0;i_scalar < nVar;i_scalar++){
          scalar_init[i_scalar] = config->GetScalar_Init(i_scalar);
        }
        
        solver_container[i_mesh][SCALAR_SOL]->GetNodes()->SetSolution(i_point, scalar_init);

      }

      solver_container[i_mesh][SCALAR_SOL]->InitiateComms(geometry[i_mesh], config, SOLUTION);
      solver_container[i_mesh][SCALAR_SOL]->CompleteComms(geometry[i_mesh], config, SOLUTION);

      solver_container[i_mesh][FLOW_SOL]->InitiateComms(geometry[i_mesh], config, SOLUTION);
      solver_container[i_mesh][FLOW_SOL]->CompleteComms(geometry[i_mesh], config, SOLUTION);

      solver_container[i_mesh][FLOW_SOL]->Preprocessing( geometry[i_mesh], solver_container[i_mesh], config, i_mesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      
    }


  delete[] scalar_init;
  }


}


void CPassiveScalarSolver::SetPreconditioner(CGeometry *geometry, CSolver **solver_container,  CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  
  su2double  BetaInc2, Density, dRhodT, dRhodC, Temperature, Delta;
  
  bool variable_density = (config->GetKind_DensityModel() == INC_DENSITYMODEL::VARIABLE);
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Access the primitive variables at this node. ---*/
    
    Density     = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    BetaInc2    = solver_container[FLOW_SOL]->GetNodes()->GetBetaInc2(iPoint);
    Temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
    
    unsigned short nVar_Flow = solver_container[FLOW_SOL]->GetnVar();
    
    su2double SolP = solver_container[FLOW_SOL]->LinSysSol[iPoint*nVar_Flow+0];
    su2double SolT = solver_container[FLOW_SOL]->LinSysSol[iPoint*nVar_Flow+nDim+1];
    
    /*--- We need the derivative of the equation of state to build the
     preconditioning matrix. For now, the only option is the ideal gas
     law, but in the future, dRhodT should be in the fluid model. ---*/
    
    if (variable_density) {
      dRhodT = -Density/Temperature;
    } else {
      dRhodT = 0.0;
    }
    
    /*--- Passive scalars have no impact on the density. ---*/
    
    dRhodC = 0.0;
    
    /*--- Modify matrix diagonal with term including volume and time step. ---*/
    
    su2double Vol = geometry->nodes->GetVolume(iPoint);
    Delta = Vol / (config->GetCFLRedCoeff_Scalar()*
                   solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint));
    
    /*--- Calculating the inverse of the preconditioning matrix
     that multiplies the time derivative during time integration. ---*/
    
    if (implicit) {
      // nijso: do we need to wipe the entire jacobian for preconditioning?
      //for (int i_var = 0; i_var < nVar; i_var++) {
      //  for (int j_var = 0; j_var < nVar; j_var++) {
      //    Jacobian_i[i_var][j_var] = 0.0;
      //  }
      //}
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        total_index = iPoint*nVar+iVar;
        
        su2double c = nodes->GetSolution(iPoint,iVar);
        
        /*--- Compute the lag terms for the decoupled linear system from
         the mean flow equations and add to the residual for the scalar.
         In short, we are effectively making these terms explicit. ---*/
        
        su2double artcompc1 = SolP * c/(Density*BetaInc2);
        su2double artcompc2 = SolT * dRhodT * c/(Density);
        
        LinSysRes[total_index] += artcompc1 + artcompc2;
        
        /*--- Add the extra Jacobian term to the scalar system. ---*/
        
        su2double Jaccomp = c * dRhodC + Density; //This is Gamma
        su2double JacTerm = Jaccomp*Delta;
        
        Jacobian.AddVal2Diag(iPoint, JacTerm);
        
      }
      
    }
    
  }

}

void CPassiveScalarSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                           CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  bool axisymmetric = config->GetAxisymmetric();

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  vector<su2double> zero_sources(nVar,0.);

  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Set scalar sources to zero ---*/
  numerics->SetScalarSources(&zero_sources[0]);

  /*--- Loop over all points. ---*/

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint), nullptr);
    
    /*--- Scalar variables w/o reconstruction, and its gradient ---*/
    
    numerics->SetScalarVarLegacy(nodes->GetSolution(iPoint), nullptr);
    numerics->SetScalarVarLegacyGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Mass diffusivity coefficients. ---*/
    su2double *M;
    M = nodes->GetDiffusivity(iPoint);
    /*--- for multiple species, we need an array---*/
    numerics->SetDiffusionCoeff(M, M);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Axisymmetry source term for the scalar equation. ---*/

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
}

void CPassiveScalarSolver::BC_Inlet(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CNumerics *conv_numerics,
                                    CNumerics *visc_numerics,
                                    CConfig *config,
                                    unsigned short val_marker) {

  su2double *inlet_scalar = new su2double[nVar]{};

  //bool grid_movement  = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  //su2double   temp_inlet    = config->GetInlet_Ttotal       (Marker_Tag);
  inlet_scalar  = config->GetInlet_ScalarVal    (Marker_Tag);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      nodes->SetSolution_Old(iPoint, inlet_scalar);

      LinSysRes.SetBlock_Zero(iPoint);

      for (auto iVar = 0u; iVar < nVar; iVar++) {
        nodes->SetVal_ResTruncError_Zero(iPoint, iVar);
      }

      /*--- Includes 1 in the diagonal ---*/
      for (auto iVar = 0u; iVar < nVar; iVar++) {
        auto total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    }
  }
  END_SU2_OMP_FOR

}

void CPassiveScalarSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                               CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /* strong zero flux Neumann boundary condition at the outlet */
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {
      
        /*--- Allocate the value at the outlet ---*/
        auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor(); 
          
        for (auto iVar = 0u; iVar < nVar; iVar++)
          Solution[iVar] = nodes->GetSolution(Point_Normal, iVar);
        nodes->SetSolution_Old(iPoint, Solution);
    
        LinSysRes.SetBlock_Zero(iPoint);

        for (auto iVar = 0u; iVar < nVar; iVar++){
          nodes->SetVal_ResTruncError_Zero(iPoint, iVar);
        }

        /*--- Includes 1 in the diagonal ---*/
        for (auto iVar = 0u; iVar < nVar; iVar++) {
          auto total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
    }
  }
  END_SU2_OMP_FOR

}

void CPassiveScalarSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  

}

void CPassiveScalarSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}


