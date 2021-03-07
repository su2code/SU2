/*!
 * \file CPassiveScalarSolver.cpp
 * \brief Main subroutines for solving transported scalar class
 * \author T. Economon, N. Beishuizen
 * \version 7.1.0 "Blackbird"
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

CPassiveScalarSolver::CPassiveScalarSolver(void) : CScalarSolver() {}

CPassiveScalarSolver::CPassiveScalarSolver(CGeometry *geometry,
                                           CConfig *config,
                                           unsigned short iMesh)
: CScalarSolver(geometry, config) {
  
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  su2double Density_Inf, Viscosity_Inf;
  
  bool turbulent = ((config->GetKind_Solver() == RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool turb_SST  = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool turb_SA   = ((turbulent) && (config->GetKind_Turb_Model() == SA));

  /*--- Dimension of the problem --> passive scalar will only ever
   have a single equation. Other child classes of CScalarSolver
   can have variable numbers of equations. ---*/
 
  nVar     = 1;
  nPrimVar = 1;

  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim = geometry->GetnDim();
  
  /*--- Fluid model pointer initialization ---*/
  
  FluidModel = NULL;
  
  /*--- Define some auxiliar vector related with the residual ---*/
  
  Residual     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual    [iVar] = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar] = 0.0;
  Residual_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i  [iVar] = 0.0;
  Residual_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j  [iVar] = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar] = 0.0;
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliar vector related with the solution ---*/
  
  Solution = new su2double[nVar];
  Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
  
  /*--- Define some auxiliar vector related with the geometry ---*/
  
  Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
  
  /*--- Define some auxiliar vector related with the flow solution ---*/
  
  FlowPrimVar_i = new su2double [nDim+9]; FlowPrimVar_j = new su2double [nDim+9];
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  /*--- Initialization of the structure of the whole Jacobian ---*/
  
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Passive Scalar). MG level: " << iMesh <<"." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  
  if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
    if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
  }
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Computation of gradients by least squares ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*transpose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    
    /*--- c vector := transpose(WA)*(Wb) ---*/
    
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Initialize lower and upper limits---*/
  
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  
  lowerlimit[0] = -1.0e15;
  upperlimit[0] =  1.0e15;
  
  /*--- Read farfield conditions from config ---*/
  
  Density_Inf   = config->GetDensity_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();

  /*--- Set up fluid model for the diffusivity ---*/
  
  su2double Diffusivity_Ref = 1.0;
  
  su2double DiffusivityND = config->GetDiffusivity_Constant()/Diffusivity_Ref;
  
  config->SetDiffusivity_ConstantND(DiffusivityND);
  
  FluidModel = new CFluidModel();
  FluidModel->SetMassDiffusivityModel(config);
  
  /*--- Scalar variable state at the far-field. ---*/
  
  Scalar_Inf = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++){
    Scalar_Inf[iVar] = config->GetScalar_Init(iVar);
  }
  /*--- Initialize the solution to the far-field state everywhere. ---*/
  
  nodes = new CPassiveScalarVariable(Scalar_Inf, nPoint, nDim, nVar, config);
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
   (not CScalarSolver), due to arbitrary number of scalar variables.
   First, we also set the column index for any inlet profiles. ---*/
  
  Inlet_Position = nDim*2+2;
  if (turbulent) {
    if (turb_SA) Inlet_Position += 1;
    else if (turb_SST) Inlet_Position += 2;
  }

 /*-- Allocation of inlets has to happen in derived classes
   (not CScalarSolver), due to arbitrary number of scalar variables.
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
  
}

CPassiveScalarSolver::~CPassiveScalarSolver(void) {
  
  unsigned long iMarker, iVertex;
  unsigned short iVar;
  
  if (FluidModel != NULL) delete FluidModel;

}

void CPassiveScalarSolver::Preprocessing(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CConfig *config,
                                         unsigned short iMesh,
                                         unsigned short iRKStep,
                                         unsigned short RunTime_EqSystem,
                                         bool Output) {
  
  unsigned long ErrorCounter = 0;
  unsigned long InnerIter = config->GetInnerIter();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool limiter_flow     = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) &&
                           (InnerIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_scalar   = ((config->GetKind_SlopeLimit_Scalar() != NO_LIMITER) &&
                           (InnerIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (!disc_adjoint) Jacobian.SetValZero();
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  
  /*--- Upwind second order reconstruction ---*/
  
  if ((limiter_scalar) && (iMesh == MESH_0)) SetSolution_Limiter(geometry, config);
  
  if ((limiter_flow) && (iMesh == MESH_0)) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);
  
}

void CPassiveScalarSolver::Postprocessing(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CConfig *config,
                                          unsigned short iMesh) { }

unsigned long CPassiveScalarSolver::SetPrimitive_Variables(CSolver **solver_container,
                                                           CConfig *config,
                                                           bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double Density, Temperature, lam_visc = 0.0, eddy_visc = 0.0, Cp, diff;
  unsigned short iVar, turb_model = config->GetKind_Turb_Model();
  CFluidModel *FluidModel= solver_container[FLOW_SOL]->GetFluidModel();
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Retrieve the density, temperature, Cp, and laminar viscosity. ---*/

    Density     = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    Cp          = solver_container[FLOW_SOL]->GetNodes()->GetSpecificHeatCp(iPoint);
    Temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
    lam_visc    = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    
    /*--- Retrieve the value of the kinetic energy (if needed) ---*/
    
    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->GetNodes()->GetmuT(iPoint);
    }

    /*--- Compute and store the mass diffusivity. ---*/
    // nijso: necessary? does not seem to matter!
    FluidModel->SetDiffusivityState(Temperature, Density, lam_visc, eddy_visc, Cp);


    for (int i_scalar = 0; i_scalar < nVar; i_scalar++)
      nodes->SetDiffusivity(iPoint, FluidModel->GetMassDiffusivity(), i_scalar);
    
    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  return ErrorCounter;
  
}

void CPassiveScalarSolver::SetInitialCondition(CGeometry **geometry,
                                               CSolver ***solver_container,
                                               CConfig *config,
                                               unsigned long ExtIter) {
  su2double *coords;
  bool Restart   = (config->GetRestart() || config->GetRestart_Flow());
  
  
  if ((!Restart) && ExtIter == 0) {
    if (rank == MASTER_NODE){
      cout << "Initializing passive scalar (initial condition)." << endl;
      cout << "initialization = " << nVar << " " << config->GetScalar_Init(0)<<endl;
    }  

    su2double *scalar_init    = new su2double[nVar];
    CFluidModel *fluid_model_local;

    for (unsigned long i_mesh = 0; i_mesh <= config->GetnMGLevels(); i_mesh++) {

      fluid_model_local = solver_container[i_mesh][FLOW_SOL]->GetFluidModel();

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
  
  su2double  BetaInc2, Density, dRhodT, dRhodC, Temperature, Cp, Delta;
  
  bool variable_density = (config->GetKind_DensityModel() == VARIABLE);
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Access the primitive variables at this node. ---*/
    
    Density     = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    BetaInc2    = solver_container[FLOW_SOL]->GetNodes()->GetBetaInc2(iPoint);
    Cp          = solver_container[FLOW_SOL]->GetNodes()->GetSpecificHeatCp(iPoint);
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
        
        su2double c = nodes->GetSolution(iPoint,0);
        
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

void CPassiveScalarSolver::Source_Residual(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CNumerics *numerics,
                                      CNumerics *second_numerics,
                                      CConfig *config,
                                      unsigned short iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool axisymmetric   = config->GetAxisymmetric();
  bool viscous        = config->GetViscous();
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint), NULL);
    
    /*--- Gradient of the primitive and conservative variables ---*/
    
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint), NULL);
    
    /*--- Scalar variables w/o reconstruction, and its gradient ---*/
    
    numerics->SetScalarVar(nodes->GetSolution(iPoint), NULL);
    numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), NULL);
    
    /*--- Set volume ---*/
    
    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));
    
    /*--- Compute the source term ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);
    
    /*--- Subtract residual and the Jacobian ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    
  }
  
  /*--- Axisymmetry source term for the scalar equation. ---*/
  
  if (axisymmetric) {
    
    /*--- Zero out Jacobian structure ---*/
    
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }
    
    /*--- loop over points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Primitive variables w/o reconstruction ---*/
      
      second_numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint), NULL);
      
      /*--- Scalar variables w/o reconstruction ---*/
      
      second_numerics->SetScalarVar(nodes->GetSolution(iPoint), NULL);
      
      /*--- Mass diffusivity coefficients. ---*/
      
      second_numerics->SetDiffusionCoeff(nodes->GetDiffusivity(iPoint),
                                         NULL);
      
      /*--- Set control volume ---*/
      
      second_numerics->SetVolume(geometry->nodes->GetVolume(iPoint));
      
      /*--- Set y coordinate ---*/
      
      second_numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                                NULL);
      
      /*--- If viscous, we need gradients for extra terms. ---*/
      
      if (viscous) {
        
        /*--- Gradient of the scalar variables ---*/
        
        second_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), NULL);
        
      }
      
      /*--- Compute Source term Residual ---*/
      
      second_numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }

}





void CPassiveScalarSolver::BC_Inlet(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CNumerics *conv_numerics,
                                    CNumerics *visc_numerics,
                                    CConfig *config,
                                    unsigned short val_marker) {


 
  unsigned short iDim, iVar;
  unsigned long iVertex, iPoint,total_index;
  su2double *V_inlet, *V_domain, *Normal, *inlet_scalar;
  
  Normal = new su2double[nDim];
  inlet_scalar = new su2double[nVar];
  for (int i=0;i<nVar;i++)
    inlet_scalar[i]=0.0;

  bool grid_movement  = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  su2double   temp_inlet    = config->GetInlet_Ttotal       (Marker_Tag);
  inlet_scalar  = config->GetInlet_ScalarVal    (Marker_Tag);

  
  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->nodes->GetDomain(iPoint)) {

      nodes->SetSolution_Old(iPoint, inlet_scalar);

      LinSysRes.SetBlock_Zero(iPoint);

      for (iVar = 0; iVar < nVar; iVar++) {
        nodes->SetVal_ResTruncError_Zero(iPoint, iVar);
      }

      /*--- Includes 1 in the diagonal ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }      
    }
  }

   /*--- Free locally allocated memory ---*/
  delete[] Normal;
  //delete[] inlet_scalar;
  
}

void CPassiveScalarSolver::BC_Outlet(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config,
                                     unsigned short val_marker) {
 
 unsigned long iPoint, iVertex, Point_Normal, total_index;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;
  
  bool grid_movement  = config->GetGrid_Movement();
  CFluidModel *fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();  
  Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /* strong zero flux Neumann boundary condition at the outlet */
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->nodes->GetDomain(iPoint)) {
      
        /*--- Allocate the value at the outlet ---*/
        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor(); 
          
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = nodes->GetSolution(Point_Normal, iVar);
        nodes->SetSolution_Old(iPoint, Solution);
    
        LinSysRes.SetBlock_Zero(iPoint);

        for (iVar = 0; iVar < nVar; iVar++){
          nodes->SetVal_ResTruncError_Zero(iPoint, iVar);
        }

        /*--- Includes 1 in the diagonal ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;

}

void CPassiveScalarSolver::BC_Isothermal_Wall(CGeometry *geometry,
                                              CSolver **solver_container,
                                              CNumerics *conv_numerics,
                                              CNumerics *visc_numerics,
                                              CConfig *config,
                                              unsigned short val_marker) {
  
  /*--- Convective fluxes across viscous walls are equal to zero. ---*/
  
}
