/*!
 * \file CSpeciesFlameletSolver.cpp
 * \brief Main subroutines of CSpeciesFlameletSolver class
 * \author T. Kattmann
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
//#include "../../include/variables/CSpeciesFlameletVariable.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/fluid/CFluidFlamelet.hpp"

CSpeciesFlameletSolver::CSpeciesFlameletSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh)
    : CSpeciesSolver(geometry, config, true) {

  /*--- Store if an implicit scheme is used, for use during periodic boundary conditions. ---*/
  SetImplicitPeriodic(config->GetKind_TimeIntScheme_Species() == EULER_IMPLICIT);

  /*--- Dimension of the problem. ---*/

  //nVar = config->GetnSpecies();
  // flamelet: only thing that changes
  nVar = 2;
  nPrimVar = nVar;

  if (nVar > MAXNVAR)
    SU2_MPI::Error("Increase static array size MAXNVAR for CSpeciesVariable and proceed.", CURRENT_FUNCTION);

  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

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
      lowerlimit[iVar] = -1.0e15;
      upperlimit[iVar] = 1.0e15;
    }
  }

  /*--- Scalar variable state at the far-field. ---*/

  for (auto iVar = 0u; iVar < nVar; iVar++) {
    Solution_Inf[iVar] = config->GetSpecies_Init()[iVar];
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CSpeciesVariable(Solution_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Initialize the mass diffusivity. Nondimensionalization done in the flow solver. ---*/
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      const auto MassDiffusivity = config->GetDiffusivity_ConstantND();
      nodes->SetDiffusivity(iPoint, MassDiffusivity, iVar);
    }
  }
  END_SU2_OMP_FOR

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


void CSpeciesFlameletSolver::SetInitialCondition(CGeometry **geometry,
                                          CSolver ***solver_container,
                                          CConfig *config,
                                          unsigned long ExtIter) {
  bool Restart   = (config->GetRestart() || config->GetRestart_Flow());
  
  
  if ((!Restart) && ExtIter == 0) {
    if (rank == MASTER_NODE){
      cout << "Initializing progress variable and temperature (initial condition)." << endl;
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
    if (rank == MASTER_NODE){
      cout << "initial condition: T = " << temp_inlet << endl; 
      cout << "initial condition: pv = " << prog_inlet << endl; 
    }

    su2double enth_inlet;
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

      n_not_iterated         += fluid_model_local->GetEnthFromTemp(&enth_inlet, prog_inlet,temp_inlet);
      scalar_init[I_ENTH] = enth_inlet;

      /*--- the burnt value of the progress variable is set to a value slightly below the maximum value ---*/

      prog_burnt = 0.95*fluid_model_local->GetLookUpTable()->GetTableLimitsProg().second;
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

        n_not_in_domain        += fluid_model_local->GetLookUpTable()->LookUp_ProgEnth(look_up_tags, look_up_data, scalar_init[I_PROGVAR], scalar_init[I_ENTH],name_prog,name_enth);

        // skip progress variable and enthalpy
        // we can make an init based on the lookup table. 
        for(int i_scalar = 0; i_scalar < config->GetnSpecies(); ++i_scalar){
          if ( (i_scalar != I_ENTH) && (i_scalar != I_PROGVAR) )
            scalar_init[i_scalar] = config->GetSpecies_Init()[i_scalar];
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