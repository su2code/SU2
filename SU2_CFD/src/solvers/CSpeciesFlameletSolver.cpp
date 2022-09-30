/*!
 * \file CSpeciesFlameletSolver.cpp
 * \brief Main subroutines of CSpeciesFlameletSolver class
 * \author D. Mayer, T. Economon, N. Beishuizen
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

#include "../../include/solvers/CSpeciesSolver.hpp"

#include "../../include/solvers/CSpeciesFlameletSolver.hpp"
#include "../../include/variables/CSpeciesFlameletVariable.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/fluid/CFluidFlamelet.hpp"

CSpeciesFlameletSolver::CSpeciesFlameletSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh)
    : CSpeciesSolver(geometry, config, true) {

  /*--- Store if an implicit scheme is used, for use during periodic boundary conditions. ---*/
  SetImplicitPeriodic(config->GetKind_TimeIntScheme_Species() == EULER_IMPLICIT);

  /*--- Dimension of the problem. ---*/

  // nijso: what about getnspecies, use that?

  // we set nscalars in fluidflamelet
  nVar = config->GetNScalars();
  // flamelet: only thing that changes
  //nVar = 2;
  nPrimVar = nVar;

  if (nVar > MAXNVAR)
    SU2_MPI::Error("Increase static array size MAXNVAR for CSpeciesVariable and proceed.", CURRENT_FUNCTION);

  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();


  /* ****************************************************************************** */  
  /* ***** additional stuff compared to CSpeciesSolver                        ***** */

  /*--- Single grid simulation ---*/

  /*--- Define some auxiliary vector related with the solution ---*/
  Solution   = new su2double[nVar];
  Solution_i = new su2double[nVar];
  Solution_j = new su2double[nVar];
  
  /*--- Allocates a 3D array with variable "middle" sizes and init to 0. ---*/

  auto Alloc3D = [](unsigned long M, const vector<unsigned long>& N, unsigned long P, vector<su2activematrix>& X) {
    X.resize(M);
    for (unsigned long i = 0; i < M; ++i) X[i].resize(N[i], P) = su2double(0.0);
  };

  /*--- Store the values of the temperature and the heat flux density at the boundaries,
   used for coupling with a solid donor cell ---*/
  constexpr auto n_conjugate_var = 4u;

  Alloc3D(nMarker, nVertex, n_conjugate_var, conjugate_var);
  for (auto& x : conjugate_var) x = config->GetTemperature_FreeStreamND();

  /* ****************************************************************************** */  
  
    /*--- Single grid simulation ---*/

  if (iMesh == MESH_0 || config->GetMGCycle() == FULLMG_CYCLE) {

    /*--- Define some auxiliary vector related with the residual ---*/

    Residual_RMS.resize(nVar, 0.0);
    Residual_Max.resize(nVar, 0.0);
    Point_Max.resize(nVar, 0);
    Point_Max_Coord.resize(nVar, nDim) = su2double(0.0);

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (config->GetMultizone_Problem()) {
      Residual_BGS.resize(nVar, 0.0);
      Residual_Max_BGS.resize(nVar, 0.0);
      Point_Max_BGS.resize(nVar, 0);
      Point_Max_Coord_BGS.resize(nVar, nDim) = su2double(0.0);
    }

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (flamelet model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      const auto nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE)
        cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    System.SetxIsZero(true);

    if (ReducerStrategy) EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
  }

  /*--- Initialize lower and upper limits---*/

  if (config->GetSpecies_Clipping()) {
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      lowerlimit[iVar] = config->GetSpecies_Clipping_Min(iVar);
      upperlimit[iVar] = config->GetSpecies_Clipping_Max(iVar);
    }
  } else {
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      /*--- we fix the lower limit to 0 ---*/  
      lowerlimit[iVar] = -1.0e15;
      upperlimit[iVar] = 1.0e15;
    }
  }

  //nijso: we probably do not need this for cflamelet
  //FluidModel = new CFluidModel();
  //FluidModel->SetMassDiffusivityModel(config);

  /*--- Scalar variable state at the far-field. ---*/

  for (auto iVar = 0u; iVar < nVar; iVar++) {
    Solution_Inf[iVar] = config->GetSpecies_Init()[iVar];
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CSpeciesFlameletVariable(Solution_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Initialize the mass diffusivity. Nondimensionalization done in the flow solver. ---*/
  //SU2_OMP_FOR_STAT(omp_chunk_size)
  //for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
  //  for (auto iVar = 0u; iVar < nVar; iVar++) {
  //    const auto MassDiffusivity = config->GetDiffusivity_ConstantND();
  //    nodes->SetDiffusivity(iPoint, MassDiffusivity, iVar);
  //  }
  //}
  //END_SU2_OMP_FOR

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Set the column number for species in inlet-files.
   * e.g. Coords(nDim), Temp(1), VelMag(1), Normal(nDim), Turb(1 or 2), Species(arbitrary) ---*/
  Inlet_Position = nDim + 2 + nDim + config->GetnTurbVar();

  /*-- Allocation of inlet-values. Will be filled either by an inlet files,
   * or uniformly by a uniform boundary condition. ---*/

  Inlet_SpeciesVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_SpeciesVars[iMarker].resize(nVertex[iMarker], nVar);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        Inlet_SpeciesVars[iMarker](iVertex, iVar) = Solution_Inf[iVar];
      }
    }
  }

  /*--- Store the initial CFL number for all grid points. ---*/

  const su2double CFL = config->GetCFL(MGLevel) * config->GetCFLRedCoeff_Species();
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0u; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  END_SU2_OMP_FOR
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "FLAMELET";
}

void CSpeciesFlameletSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
         unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long n_not_in_domain     = 0;
  unsigned long global_table_misses = 0;

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool muscl    =  config->GetMUSCL_Species();
  //const bool limiter  = (config->GetKind_SlopeLimit_Species() != SLOPE_LIMITER::NO_LIMITER) &&
  //                      (config->GetInnerIter() <= config->GetLimiterIter());

  for (auto i_point = 0u; i_point < nPoint; i_point++) {

    CFluidModel * fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();
    su2double * scalars = nodes->GetSolution(i_point);

    fluid_model_local->SetTDState_T(0.0,scalars);

    // nijso TODO n_not_in_domain has to come from somewhere else
    //n_not_in_domain += fluid_model_local->SetTDState_T(0.0,scalars); /*--- first argument (temperature) is not used ---*/
    fluid_model_local->SetTDState_T(0.0,scalars); /*--- first argument (temperature) is not used ---*/

    /*--- set the diffusivity in the fluid model to the diffusivity obtained from the lookup table ---*/  
    for(auto i_scalar = 0u; i_scalar < config->GetNScalars(); ++i_scalar){
      nodes->SetDiffusivity(i_point, fluid_model_local->GetMassDiffusivity(), i_scalar);
    }

    if (!Output) LinSysRes.SetBlock_Zero(i_point);

  }

  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Reduce(&n_not_in_domain, &global_table_misses, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
    if (rank == MASTER_NODE) {
      SetNTableMisses(global_table_misses);
    } 
  }

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

  //if (limiter && muscl) SetSolution_Limiter(geometry, config);
}

void CSpeciesFlameletSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container,
                                    CConfig *config, unsigned short iMesh) {

/*--- your postprocessing goes here ---*/

}

void CSpeciesFlameletSolver::SetInitialCondition(CGeometry **geometry,
                                          CSolver ***solver_container,
                                          CConfig *config,
                                          unsigned long ExtIter) {
  bool Restart   = (config->GetRestart() || config->GetRestart_Flow());
  
  
  if ((!Restart) && ExtIter == 0) {
    if (rank == MASTER_NODE){
      cout << "Initializing progress variable and total enthalpy (using temperature)" << endl;
    }

    su2double *scalar_init  = new su2double[nVar];
    su2double *flame_offset = config->GetFlameOffset();
    su2double *flame_normal = config->GetFlameNormal();

    su2double prog_burnt;
    su2double prog_unburnt    = 0.0;
    su2double flame_thickness = config->GetFlameThickness();
    su2double burnt_thickness = config->GetFlameBurntThickness();
    su2double flamenorm       = sqrt( flame_normal[0]*flame_normal[0]
                                     +flame_normal[1]*flame_normal[1]
                                     +flame_normal[2]*flame_normal[2]);

   
    su2double temp_inlet = config->GetInc_Temperature_Init(); // should do reverse lookup of enthalpy
    su2double prog_inlet = config->GetSpecies_Init()[I_PROGVAR]; 
    su2double enth_inlet = config->GetSpecies_Init()[I_ENTH]; 
    if (rank == MASTER_NODE){
      cout << "initial condition: T = " << temp_inlet << endl; 
      cout << "initial condition: c = " << prog_inlet << endl; 
      cout << "initial condition: h = " << enth_inlet << endl; 
    }

    //su2double enth_inlet;
    su2double point_loc;
    unsigned long n_not_iterated  = 0;
    unsigned long n_not_in_domain = 0;

    CFluidModel *fluid_model_local;

    vector<string>     look_up_tags;
    vector<su2double*> look_up_data;
    string name_enth = config->GetLUTScalarName(I_ENTH);
    string name_prog = config->GetLUTScalarName(I_PROGVAR);

    for (unsigned long i_mesh = 0; i_mesh <= config->GetnMGLevels(); i_mesh++) {

      fluid_model_local = solver_container[i_mesh][FLOW_SOL]->GetFluidModel();
    
      prog_burnt = fluid_model_local->GetLookUpTable()->GetTableLimitsProg().second;
      for (unsigned long i_point = 0; i_point < nPointDomain; i_point++) {
        
        for (unsigned long i_var = 0; i_var < nVar; i_var++)
          Solution[i_var] = 0.0;
        
        auto coords = geometry[i_mesh]->nodes->GetCoord(i_point);

        /* determine if our location is above or below the plane, assuming the normal 
           is pointing towards the burned region*/ 
        point_loc = 0.0;
        for (unsigned short i_dim =0; i_dim < geometry[i_mesh]->GetnDim(); i_dim++) {
          point_loc +=  flame_normal[i_dim]*(coords[i_dim]-flame_offset[i_dim]); 
        }

        /* compute the exact distance from point to plane */          
        point_loc = point_loc/flamenorm;

        /* --- unburnt region upstream of the flame --- */
        if (point_loc <= 0){
          scalar_init[I_PROGVAR] = prog_unburnt;

         /* --- flame zone --- */
        } else if ( (point_loc > 0) && (point_loc <= flame_thickness) ){
          scalar_init[I_PROGVAR] = prog_unburnt + point_loc * (prog_burnt - prog_unburnt)/flame_thickness;

        /* --- burnt region --- */
        } else if ( (point_loc > flame_thickness) && (point_loc <= flame_thickness + burnt_thickness) ){
          scalar_init[I_PROGVAR] = prog_burnt;

        /* --- unburnt region downstream of the flame --- */
        } else {
          scalar_init[I_PROGVAR] = prog_unburnt;
        }

        n_not_in_domain += fluid_model_local->GetLookUpTable()->LookUp_ProgEnth(look_up_tags, look_up_data, scalar_init[I_PROGVAR], scalar_init[I_ENTH],name_prog,name_enth);
        n_not_iterated += fluid_model_local->GetEnthFromTemp(&enth_inlet, prog_inlet, temp_inlet);
        scalar_init[I_ENTH] = enth_inlet;

        // initialize other transported scalars  (not pv and enth)
        // skip progress variable and enthalpy
        // we can make an init based on the lookup table. 
        for(int i_scalar = 2; i_scalar < config->GetNScalars(); ++i_scalar){
          //if ( (i_scalar != I_ENTH) && (i_scalar != I_PROGVAR) ){
            scalar_init[i_scalar] = config->GetSpecies_Init()[i_scalar];
          //}
        }

        solver_container[i_mesh][SPECIES_SOL]->GetNodes()->SetSolution(i_point, scalar_init);

      }

      solver_container[i_mesh][SPECIES_SOL]->InitiateComms(geometry[i_mesh], config, SOLUTION);
      solver_container[i_mesh][SPECIES_SOL]->CompleteComms(geometry[i_mesh], config, SOLUTION);

      solver_container[i_mesh][FLOW_SOL]->InitiateComms(geometry[i_mesh], config, SOLUTION);
      solver_container[i_mesh][FLOW_SOL]->CompleteComms(geometry[i_mesh], config, SOLUTION);

      solver_container[i_mesh][FLOW_SOL]->Preprocessing( geometry[i_mesh], solver_container[i_mesh], config, i_mesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      
    }

    delete[] scalar_init;

    if (rank == MASTER_NODE && (n_not_in_domain > 0 || n_not_iterated > 0))
      cout << endl;
    
    if (rank == MASTER_NODE && n_not_in_domain > 0)
      cout << " !!! Initial condition: Number of points outside of table domain: " << n_not_in_domain << " !!!" << endl;

    if (rank == MASTER_NODE && n_not_iterated > 0)
      cout << " !!! Initial condition: Number of points in which enthalpy could not be iterated: " << n_not_iterated << " !!!" << endl;

    if (rank == MASTER_NODE && (n_not_in_domain > 0 || n_not_iterated > 0))
      cout << endl;
  }
}

void CSpeciesFlameletSolver::SetPreconditioner(CGeometry *geometry, CSolver **solver_container,  CConfig *config) {
  
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
    Delta = Vol / (config->GetCFLRedCoeff_Species()*
                   solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint));
    
    /*--- Calculating the inverse of the preconditioning matrix
     that multiplies the time derivative during time integration. ---*/
    
    if (implicit) {
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        total_index = iPoint*nVar+iVar;
        
        su2double scalar = nodes->GetSolution(iPoint,iVar);
        
        /*--- Compute the lag terms for the decoupled linear system from
         the mean flow equations and add to the residual for the scalar.
         In short, we are effectively making these terms explicit. ---*/
        
        su2double artcompc1 = SolP * scalar/(Density*BetaInc2);
        su2double artcompc2 = SolT * dRhodT * scalar/(Density);
        
        LinSysRes[total_index] += artcompc1 + artcompc2;
        
        /*--- Add the extra Jacobian term to the scalar system. ---*/
        
        su2double Jaccomp = scalar * dRhodC + Density; //This is Gamma
        su2double JacTerm = Jaccomp*Delta;
        
        Jacobian.AddVal2Diag(iPoint, JacTerm);
        
      }
      
    }
    
  }

}

void CSpeciesFlameletSolver::Source_Residual(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CNumerics **numerics_container,
                                      CConfig *config,
                                      unsigned short iMesh) {

  bool implicit            = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool axisymmetric        = config->GetAxisymmetric();

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

  auto* first_numerics  = numerics_container[SOURCE_FIRST_TERM  + omp_get_thread_num()*MAX_TERMS];

  CFluidModel *fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();
  
  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (auto i_point = 0u; i_point < nPointDomain; i_point++) {

    /*--- Set primitive variables w/o reconstruction ---*/

    first_numerics->SetPrimitive(flowNodes->GetPrimitive(i_point), nullptr);

    /*--- Set scalar variables w/o reconstruction ---*/

    first_numerics->SetScalarVar(nodes->GetSolution(i_point), nullptr);

    first_numerics->SetDiffusionCoeff(nodes->GetDiffusivity(i_point), nodes->GetDiffusivity(i_point));

    /*--- Set volume of the dual cell. ---*/

    first_numerics->SetVolume(geometry->nodes->GetVolume(i_point));

    /*--- Update scalar sources in the fluidmodel ---*/

    /*--- Axisymmetry source term for the scalar equation. ---*/
    if (axisymmetric){
      /*--- Set y coordinate ---*/
      first_numerics->SetCoord(geometry->nodes->GetCoord(i_point), geometry->nodes->GetCoord(i_point));
      /*-- gradients necessary for axisymmetric flows only? ---*/
      first_numerics->SetScalarVarGradient(nodes->GetGradient(i_point), nullptr);
   }

    fluid_model_local->SetScalarSources(nodes->GetSolution(i_point));

    /*--- Retrieve scalar sources from fluidmodel and update numerics class data. ---*/
    first_numerics->SetScalarSources(fluid_model_local->GetScalarSources());

    auto residual = first_numerics->ComputeResidual(config);

    /*--- Add Residual ---*/
    
    LinSysRes.SubtractBlock(i_point, residual);
    
    /*--- Implicit part ---*/
    
    if (implicit) Jacobian.SubtractBlock2Diag(i_point, residual.jacobian_i);
    
  }
  END_SU2_OMP_FOR
}

void CSpeciesFlameletSolver::BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                              CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  // nijso: remember that we now have progress variable and enthalpy as inlet values, but we want temperature as inlet
  // (could be an additional keyword)
  
  su2double enth_inlet;
  su2double temp_inlet    = config->GetInlet_Ttotal(Marker_Tag);
  const su2double  *inlet_scalar_original  = config->GetInlet_SpeciesVal(Marker_Tag);
  su2double  *inlet_scalar  = const_cast<su2double*>(inlet_scalar_original);

  CFluidModel  *fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();
  
  fluid_model_local->GetEnthFromTemp(&enth_inlet, inlet_scalar[I_PROGVAR], temp_inlet);
  inlet_scalar[I_ENTH] = enth_inlet;

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    if (config->GetSpecies_StrongBC()) {
      nodes->SetSolution_Old(iPoint, inlet_scalar);

      LinSysRes.SetBlock_Zero(iPoint);

      for (auto iVar = 0; iVar < nVar; iVar++) {
        nodes->SetVal_ResTruncError_Zero(iPoint, iVar);
      }

      /*--- Includes 1 in the diagonal ---*/
      for (auto iVar = 0u; iVar < nVar; iVar++) {
        auto total_index = iPoint * nVar + iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    } else {  // weak BC
      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++) Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/

      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the species variable state at the inlet. ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Inlet_SpeciesVars[val_marker][iVertex]);
      //conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), inlet_scalar);

      /*--- Set various other quantities in the solver class ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/
      const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR
}


void CSpeciesFlameletSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
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

void CSpeciesFlameletSolver::BC_HeatFlux_Wall(CGeometry *geometry, 
                                       CSolver **solver_container, 
                                       CNumerics *conv_numerics,
                                       CNumerics *visc_numerics, 
                                       CConfig *config, 
                                       unsigned short val_marker) {
}

void CSpeciesFlameletSolver::BC_Isothermal_Wall(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CNumerics *conv_numerics,
                                         CNumerics *visc_numerics,
                                         CConfig *config,
                                         unsigned short val_marker) {

  unsigned long iVertex, iPoint, total_index;

  bool implicit                   = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  string Marker_Tag               = config->GetMarker_All_TagBound(val_marker);
  su2double temp_wall             = config->GetIsothermal_Temperature(Marker_Tag);
  CFluidModel *fluid_model_local  = solver_container[FLOW_SOL]->GetFluidModel();    
  su2double enth_wall, prog_wall;
  unsigned long n_not_iterated    = 0;

  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {
      /*--- Set enthalpy on the wall ---*/

      prog_wall = solver_container[SPECIES_SOL]->GetNodes()->GetSolution(iPoint)[I_PROGVAR];
      n_not_iterated += fluid_model_local->GetEnthFromTemp(&enth_wall, prog_wall, temp_wall);

      /*--- Impose the value of the enthalpy as a strong boundary
      condition (Dirichlet) and remove any
      contribution to the residual at this node. ---*/

      nodes->SetSolution(iPoint, I_ENTH, enth_wall);
      nodes->SetSolution_Old(iPoint, I_ENTH, enth_wall);

      LinSysRes(iPoint, I_ENTH) = 0.0;

      nodes->SetVal_ResTruncError_Zero(iPoint, I_ENTH);

      if (implicit) {
        total_index = iPoint * nVar + I_ENTH;

        Jacobian.DeleteValsRowi(total_index);
      }
    }
  }
  if (rank == MASTER_NODE && n_not_iterated > 0){
    cout << " !!! Isothermal wall bc ("  << Marker_Tag << "): Number of points in which enthalpy could not be iterated: " << n_not_iterated << " !!!" << endl;
  }

}

void CSpeciesFlameletSolver::BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                              CConfig *config, unsigned short val_marker) {

  unsigned long iVertex, iPoint, total_index;

  bool implicit                   = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  string Marker_Tag               = config->GetMarker_All_TagBound(val_marker);
  su2double temp_wall             = config->GetIsothermal_Temperature(Marker_Tag);
  CFluidModel *fluid_model_local  = solver_container[FLOW_SOL]->GetFluidModel();    
  su2double enth_wall, prog_wall;
  unsigned long n_not_iterated    = 0;

  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();


    temp_wall = GetConjugateHeatVariable(val_marker, iVertex, 0);

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {
      /*--- Set enthalpy on the wall ---*/

      prog_wall = solver_container[SPECIES_SOL]->GetNodes()->GetSolution(iPoint)[I_PROGVAR];
      n_not_iterated += fluid_model_local->GetEnthFromTemp(&enth_wall, prog_wall, temp_wall);

      /*--- Impose the value of the enthalpy as a strong boundary
      condition (Dirichlet) and remove any
      contribution to the residual at this node. ---*/

      nodes->SetSolution(iPoint, I_ENTH, enth_wall);
      nodes->SetSolution_Old(iPoint, I_ENTH, enth_wall);

      LinSysRes(iPoint, I_ENTH) = 0.0;

      nodes->SetVal_ResTruncError_Zero(iPoint, I_ENTH);

      if (implicit) {
        total_index = iPoint * nVar + I_ENTH;

        Jacobian.DeleteValsRowi(total_index);
      }
    }
  }
  if (rank == MASTER_NODE && n_not_iterated > 0){
    cout << " !!! CHT interface ("  << Marker_Tag << "): Number of points in which enthalpy could not be iterated: " << n_not_iterated << " !!!" << endl;
  }
}

