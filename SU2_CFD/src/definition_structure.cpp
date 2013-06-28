/*!
 * \file definition_structure.cpp
 * \brief Main subroutines used by SU2_CFD.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/definition_structure.hpp"

void Geometrical_Definition(CGeometry **geometry, CConfig *config) {
	unsigned short iMGlevel;
	int rank = MASTER_NODE;
	
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
	
	/*--- Single grid geometry structure generation ---*/
	
	if (config->GetKind_Solver() != NO_SOLVER) {
				
		if ((rank == MASTER_NODE) || (rank == 1)) 
			cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
		
		/*--- Compute elements surrounding points, points surrounding points, and elements surrounding elements ---*/
		if ((rank == MASTER_NODE) || (rank == 1)) cout << "Setting local point and element connectivity." <<endl; 
		geometry[MESH_0]->SetEsuP(); 
		geometry[MESH_0]->SetPsuP(); 
		geometry[MESH_0]->SetEsuE();
		
		/*--- Check the orientation before computing geometrical quantities ---*/
		if ((rank == MASTER_NODE) || (rank == 1)) cout << "Checking the numerical grid orientation." <<endl; 
		geometry[MESH_0]->SetBoundVolume();
		geometry[MESH_0]->Check_Orientation(config);
		
		/*--- Create the edge structure ---*/
		if ((rank == MASTER_NODE) || (rank == 1)) cout << "Identifying edges and vertices." <<endl; 
		geometry[MESH_0]->SetEdges();
		geometry[MESH_0]->SetVertex(config);
		
		/*--- Compute center of gravity ---*/
		if ((rank == MASTER_NODE) || (rank == 1)) cout << "Computing centers of gravity." << endl; 
		geometry[MESH_0]->SetCG();
		
		/*--- Create the control volume structures ---*/
		if ((rank == MASTER_NODE) || (rank == 1)) cout << "Setting the control volume structure." << endl; 
		geometry[MESH_0]->SetControlVolume(config, ALLOCATE);
		geometry[MESH_0]->SetBoundControlVolume(config, ALLOCATE);		
		
		/*--- For a rotating frame, set the velocity due to rotation at each mesh point ---*/
		if (config->GetRotating_Frame()) geometry[MESH_0]->SetRotationalVelocity(config);
		if ((config->GetMGLevels() != 0) && ((rank == MASTER_NODE) || (rank == 1))) cout << "Setting the multigrid structure." <<endl; 

	}
	
#ifndef NO_MPI
	/*--- Synchronization point after the multigrid stuff ---*/
	MPI::COMM_WORLD.Barrier();
#endif	
	
	/*--- Multigrid geometry generation ---*/
	if (config->GetKind_Solver() != NO_SOLVER) {

		/*--- Loop over all the new grid ---*/
		for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
			
			/*--- Create main aglomeration structure (ingluding MPI stuff) ---*/
			geometry[iMGlevel] = new CMultiGridGeometry(geometry[iMGlevel-1], config, iMGlevel);
			
			/*--- Compute element surronding points, point surronding point, and element surronding elements ---*/
			geometry[iMGlevel]->SetPsuP(geometry[iMGlevel-1]);
			
			/*--- Create the edge structure ---*/
			geometry[iMGlevel]->SetEdges();
			geometry[iMGlevel]->SetVertex(geometry[iMGlevel-1], config);
			
			/*--- Create the control volume structures ---*/
			geometry[iMGlevel]->SetControlVolume(config,geometry[iMGlevel-1], ALLOCATE);
			geometry[iMGlevel]->SetBoundControlVolume(config,geometry[iMGlevel-1], ALLOCATE);
			geometry[iMGlevel]->SetCoord(geometry[iMGlevel-1]);
			
			/*--- For a rotating frame, set the velocity due to rotation at each mesh point ---*/
			if (config->GetRotating_Frame()) geometry[iMGlevel]->SetRotationalVelocity(config);
			
		}
	}
	
}

void Solver_Definition(CNumerics ****solver_container, CSolution ***solution_container, CIntegration **integration_container, 
											 CGeometry **geometry, CConfig *config) {
	unsigned short iMGlevel, iSol, nDim, nVar_Template = 0, nVar_Flow = 0, nVar_Adj_Flow = 0, nVar_Plasma = 0, nVar_LevelSet = 0, nVar_Combustion = 0, nVar_Turb = 0, nVar_Adj_Turb = 0, nVar_Eikonal = 0, 
	nVar_Elec = 0, nVar_Adj_Elec = 0, nVar_Lin_Flow = 0, nVar_Lin_Elec = 0, nSpecies = 0, nFluids = 0, nDiatomics = 0, nMonatomics = 0;
	bool no_solver, potential, euler, navierstokes, combustion, plasma, plasma_monatomic, plasma_diatomic, levelset, adj_pot, lin_pot, adj_euler, lin_euler, adj_ns, lin_ns, turbulent, 
	adj_turb, lin_turb, electric, adj_elec, lin_elec, spalart_allmaras, menter_sst, template_solver;
	
	bool incompressible = config->GetIncompressible();
	
	/*--- Initialize some useful booleans ---*/
	no_solver = false;
	potential = false;	euler = false;		navierstokes = false;	combustion = false; turbulent = false;	electric = false;	plasma_monatomic = false;	
	plasma_diatomic = false; levelset = false; plasma = false;
	adj_pot = false;	adj_euler = false;	adj_ns = false;			adj_turb = false;	adj_elec = false; 	spalart_allmaras = false;
	lin_pot = false;	lin_euler = false;	lin_ns = false;			lin_turb = false;	lin_elec = false; 	menter_sst = false;
	template_solver = false;
	
	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
		case NO_SOLVER: no_solver = true; break;
		case TEMPLATE_SOLVER: template_solver = true; break;
		case POTENTIAL_FLOW: potential = true; break;
		case EULER : euler = true; break;
		case NAVIER_STOKES: navierstokes = true; break;
		case TWO_PHASE_FLOW: euler = true; navierstokes = false; levelset = true; break;
		case RANS : navierstokes = true; turbulent = true; break;
		case MULTI_SPECIES_NAVIER_STOKES : plasma = true; electric = true; break;
		case COMBUSTION: euler = true; combustion = true; break;
		case ELECTRIC_POTENTIAL: electric = true; break;
		case ADJ_EULER : case ONE_SHOT_ADJ_EULER : euler = true; adj_euler = true; break;
		case ADJ_NAVIER_STOKES : case ONE_SHOT_ADJ_NAVIER_STOKES : navierstokes = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
		case ADJ_RANS : navierstokes = true; turbulent = true; adj_ns = true; adj_turb = true; break;
		case ADJ_ELECTRIC_POTENTIAL: electric = true; adj_elec = true; break;
		case LIN_EULER: euler = true; lin_euler = true; break;
		case LIN_ELECTRIC_POTENTIAL: electric = true; lin_elec = true; break;
	}
	
	/*--- Assign turbulence model booleans --- */
	if (turbulent)
		switch (config->GetKind_Turb_Model()){
			case SA: case SA_COMP: spalart_allmaras = true; break;
			case SST: menter_sst = true; break;
			default: cout << "Specified turbulence model unavailable or none selected" << endl; cin.get(); break;
		}
	
	if (plasma)
		switch (config->GetKind_GasModel()){
			case AIR7: plasma_diatomic = true; break;
			case ARGON: plasma_monatomic = true; break;
			default: cout << "Specified plasma model unavailable or none selected" << endl; cin.get(); break;
		}
	
	if (no_solver) {
		integration_container[TEMPLATE_SOL] = new CIntegration(config);
		integration_container[FLOW_SOL] = new CIntegration(config);
		integration_container[TURB_SOL] = new CIntegration(config);
		integration_container[ELEC_SOL] = new CIntegration(config);
		integration_container[ADJFLOW_SOL] = new CIntegration(config);
		integration_container[ADJTURB_SOL] = new CIntegration(config);
		integration_container[ADJELEC_SOL] = new CIntegration(config);
		integration_container[LINFLOW_SOL] = new CIntegration(config);
		integration_container[LINELEC_SOL] = new CIntegration(config);
		integration_container[EIKONAL_SOL] = new CIntegration(config);
		integration_container[PLASMA_SOL] = new CIntegration(config);
		integration_container[LEVELSET_SOL] = new CIntegration(config);
		integration_container[COMBUSTION_SOL] = new CIntegration(config);
		return;
	}
	
	/*--- Definition of the Class for the solution: solution_container[MESH_LEVEL][EQUATION]. Note that euler, navierstokes 
	 and potential are incompatible, they use the same position in sol container ---*/
	for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
		
		/*--- Allocate solution for a template problem ---*/
		if (template_solver) {
			solution_container[iMGlevel][TEMPLATE_SOL] = new CTemplateSolution(geometry[iMGlevel], config);
			integration_container[TEMPLATE_SOL] = new CSingleGridIntegration(config);
		}
		
		/*--- Allocate solution for direct problem ---*/
		if (potential) {
			solution_container[iMGlevel][FLOW_SOL] = new CPotentialSolution(geometry[iMGlevel], config);
			integration_container[FLOW_SOL] = new CPotentialIntegration(config);
		}
		if (euler) {
			solution_container[iMGlevel][FLOW_SOL] = new CEulerSolution(geometry[iMGlevel], config);
			integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
		}
		if (combustion) {
			solution_container[iMGlevel][COMBUSTION_SOL] = new CCombustionSolution(geometry[iMGlevel], config);
			integration_container[COMBUSTION_SOL] = new CSingleGridIntegration(config);
		}
		if (navierstokes) {
			solution_container[iMGlevel][FLOW_SOL] = new CNSSolution(geometry[iMGlevel], config);
			integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
		}
		if (turbulent) {
			if (spalart_allmaras) solution_container[iMGlevel][TURB_SOL] = new CTurbSASolution(geometry[iMGlevel], config);
			else if (menter_sst) solution_container[iMGlevel][TURB_SOL] = new CTurbSSTSolution(geometry[iMGlevel], config);
			integration_container[TURB_SOL] = new CSingleGridIntegration(config);
			solution_container[iMGlevel][EIKONAL_SOL] = new CEikonalSolution(geometry[iMGlevel], config);
			integration_container[EIKONAL_SOL] = new CEikonalIntegration(config);
		}
		if (electric) {
			solution_container[iMGlevel][ELEC_SOL] = new CElectricSolution(geometry[iMGlevel], config);
			integration_container[ELEC_SOL] = new CPotentialIntegration(config);
		}
		if (plasma) {
			if (plasma_monatomic) solution_container[iMGlevel][PLASMA_SOL] = new CPlasmaMonatomicSolution(geometry[iMGlevel], config);
			else if (plasma_diatomic) solution_container[iMGlevel][PLASMA_SOL] = new CPlasmaDiatomicSolution(geometry[iMGlevel], config);
			integration_container[PLASMA_SOL] = new CMultiGridIntegration(config);
		}
		if (levelset) {
			solution_container[iMGlevel][LEVELSET_SOL] = new CLevelSetSolution(geometry[iMGlevel], config);
			integration_container[LEVELSET_SOL] = new CSingleGridIntegration(config);
		}
		
		/*--- Allocate solution for adjoint problem ---*/
		if (adj_pot) { 
			cout <<"Equation not implemented." << endl; cin.get(); break; 
		}
		if (adj_euler) {
			solution_container[iMGlevel][ADJFLOW_SOL] = new CAdjEulerSolution(geometry[iMGlevel], config);
			integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
		}
		if (adj_ns) {
			solution_container[iMGlevel][ADJFLOW_SOL] = new CAdjNSSolution(geometry[iMGlevel], config);
			integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
		}
		if (adj_turb) {
			solution_container[iMGlevel][ADJTURB_SOL] = new CAdjTurbSolution(geometry[iMGlevel], config);
			integration_container[ADJTURB_SOL] = new CSingleGridIntegration(config);
		}
		if (adj_elec) {
			solution_container[iMGlevel][ADJELEC_SOL] = new CAdjElectricSolution(geometry[iMGlevel], config);
			integration_container[ADJELEC_SOL] = new CPotentialIntegration(config);
		}
		
		/*--- Allocate solution for linear problem (at the moment we use the same scheme as the adjoint problem) ---*/
		if (lin_pot) { 
			cout <<"Equation not implemented." << endl; cin.get(); break; 
		}
		if (lin_euler) {
			solution_container[iMGlevel][LINFLOW_SOL] = new CLinEulerSolution(geometry[iMGlevel], config);
			integration_container[LINFLOW_SOL] = new CMultiGridIntegration(config);
		}
		if (lin_ns) { 
			cout <<"Equation not implemented." << endl; cin.get(); break; 
		}		
		if (lin_elec) {
			solution_container[iMGlevel][LINELEC_SOL] = new CLinElectricSolution(geometry[iMGlevel], config);
			integration_container[LINELEC_SOL] = new CPotentialIntegration(config);
		}
		
	}
	
	/*--- Number of variables for the template ---*/
	if (template_solver) nVar_Flow = solution_container[MESH_0][FLOW_SOL]->GetnVar();

	/*--- Number of variables for direct problem ---*/
	if (potential)		nVar_Flow = solution_container[MESH_0][FLOW_SOL]->GetnVar();
	if (euler)				nVar_Flow = solution_container[MESH_0][FLOW_SOL]->GetnVar();
	if (navierstokes)	nVar_Flow = solution_container[MESH_0][FLOW_SOL]->GetnVar();
	if (turbulent) {	
		nVar_Turb     = solution_container[MESH_0][TURB_SOL]->GetnVar();
		nVar_Eikonal  = solution_container[MESH_0][EIKONAL_SOL]->GetnVar(); 
	}
	if (electric)		nVar_Elec = solution_container[MESH_0][ELEC_SOL]->GetnVar();
	if (plasma)	{ 
		nVar_Plasma = solution_container[MESH_0][PLASMA_SOL]->GetnVar();
		nSpecies    = solution_container[MESH_0][PLASMA_SOL]->GetnSpecies();
		nFluids     = solution_container[MESH_0][PLASMA_SOL]->GetnFluids();
		nDiatomics  = solution_container[MESH_0][PLASMA_SOL]->GetnDiatomics();
		nMonatomics = solution_container[MESH_0][PLASMA_SOL]->GetnMonatomics();
	}
	if (levelset)		nVar_LevelSet = solution_container[MESH_0][LEVELSET_SOL]->GetnVar();
	if (combustion)		nVar_Combustion = solution_container[MESH_0][COMBUSTION_SOL]->GetnVar();
	
	/*--- Number of variables for adjoint problem ---*/
	if (adj_pot)		nVar_Adj_Flow = solution_container[MESH_0][ADJFLOW_SOL]->GetnVar();
	if (adj_euler)	nVar_Adj_Flow = solution_container[MESH_0][ADJFLOW_SOL]->GetnVar();
	if (adj_ns)			nVar_Adj_Flow = solution_container[MESH_0][ADJFLOW_SOL]->GetnVar();
	if (adj_turb)		nVar_Adj_Turb = solution_container[MESH_0][ADJTURB_SOL]->GetnVar();
	if (adj_elec)		nVar_Adj_Elec = solution_container[MESH_0][ADJELEC_SOL]->GetnVar();
	
	/*--- Number of variables for the linear problem ---*/
	if (lin_euler)	nVar_Lin_Flow = solution_container[MESH_0][LINFLOW_SOL]->GetnVar();
	if (lin_elec)		nVar_Lin_Elec = solution_container[MESH_0][LINELEC_SOL]->GetnVar();
	
	/*--- Number of dimensions ---*/
	nDim = geometry[MESH_0]->GetnDim();
	
	/*--- Definition of the Class for the numerical method: solver_container[MESH_LEVEL][EQUATION][EQ_TERM] ---*/
	for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
		solver_container[iMGlevel] = new CNumerics** [MAX_SOLS];
		for (iSol = 0; iSol < MAX_SOLS; iSol++)
			solver_container[iMGlevel][iSol] = new CNumerics* [MAX_TERMS];
	}
	
	/*--- Solver definition for the template problem ---*/
	if (template_solver) {
		
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Template()) {
			case SPACE_CENTRED : case SPACE_UPWIND :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][TEMPLATE_SOL][CONV_TERM] = new CConvective_Template(nDim, nVar_Template, config);
				break;
			default : cout << "Convective scheme not implemented." << endl; cin.get(); break;
		}
		
		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Template()) {
			case AVG_GRAD : case AVG_GRAD_CORRECTED : case GALERKIN :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][TEMPLATE_SOL][VISC_TERM] = new CViscous_Template(nDim, nVar_Template, config);
				break;
			default : cout << "Viscous scheme not implemented." << endl; cin.get(); break; 
		}
		
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Template()) {
			case PIECEWISE_CONSTANT :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][TEMPLATE_SOL][SOUR_TERM] = new CSource_Template(nDim, nVar_Template, config);
				break;
			default : cout << "Source term not implemented." << endl; cin.get(); break;
		}
		
		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
			solver_container[iMGlevel][TEMPLATE_SOL][BOUND_TERM] = new CConvective_Template(nDim, nVar_Template, config);	
		}
		
	}
	
	/*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
	if ((potential) || (euler) || (navierstokes)) {
		
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Flow()) {
			case NO_CONVECTIVE :
				cout << "No convective scheme." << endl; cin.get();
				break;
				
			case SPACE_CENTRED :
				if (incompressible) {
					/*--- Incompressible flow, use artificial compressibility method ---*/
					switch (config->GetKind_Centred_Flow()) {
						case NO_CENTRED : cout << "No centered scheme." << endl; break;	
						case LAX : solver_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim, nVar_Flow, config); break;
						case JST : solver_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJSTArtComp_Flow(nDim, nVar_Flow, config); break;
						default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
					}
					for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
						solver_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim,nVar_Flow, config);
					
					/*--- Definition of the boundary condition method ---*/
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
						solver_container[iMGlevel][FLOW_SOL][BOUND_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config);
					
				}
				else {
					/*--- Compressible flow ---*/
					switch (config->GetKind_Centred_Flow()) {
						case NO_CENTRED : cout << "No centered scheme." << endl; break;	
						case LAX : solver_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLaxComp_Flow(nDim,nVar_Flow, config); break;
						case JST : solver_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim,nVar_Flow, config); break;
						default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
					}
					for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
						solver_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLaxComp_Flow(nDim,nVar_Flow, config);

					/*--- Definition of the boundary condition method ---*/
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
						solver_container[iMGlevel][FLOW_SOL][BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);					
				
				}
				break;
			case SPACE_UPWIND :
				if (incompressible) {
					/*--- Incompressible flow, use artificial compressibility method ---*/
					switch (config->GetKind_Upwind_Flow()) {
						case NO_UPWIND : cout << "No upwind scheme." << endl; break;	
						case ROE_1ST : case ROE_2ND :
							solver_container[MESH_0][FLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config); break;
						default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
					}
					for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
						solver_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config);

					/*--- Definition of the boundary condition method ---*/
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
						solver_container[iMGlevel][FLOW_SOL][BOUND_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config);					
				
				}
				else {
					/*--- Compressible flow ---*/
					switch (config->GetKind_Upwind_Flow()) {
						case NO_UPWIND : cout << "No upwind scheme." << endl; break;	
						case ROE_1ST : case ROE_2ND :
							for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
								solver_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
								solver_container[iMGlevel][FLOW_SOL][BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
							}
							break;
							
						case AUSM_1ST : case AUSM_2ND :
							for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
								solver_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
								solver_container[iMGlevel][FLOW_SOL][BOUND_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
							}
							break;
						
						case HLLC_1ST : case HLLC_2ND :
							for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
								solver_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
								solver_container[iMGlevel][FLOW_SOL][BOUND_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
							}
							break;
							
						default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
					}

				}
				break;
				
			default :
				cout << "Convective scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Flow()) {
			case NONE :
				break;
			case AVG_GRAD :
				if (incompressible) {
					/*--- Incompressible flow, use artificial compressibility method ---*/
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
						solver_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
				}
				else {
					/*--- Compressible flow ---*/
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
						solver_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
				}
				break;
			case AVG_GRAD_CORRECTED :
				if (incompressible) {
					/*--- Incompressible flow, use artificial compressibility method ---*/
					solver_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_Flow(nDim, nVar_Flow, config);
					for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
						solver_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
				}
				else {
					/*--- Compressible flow ---*/
					solver_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrected_Flow(nDim, nVar_Flow, config);
					for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
						solver_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
				}
				break;
			case GALERKIN :
				if (potential) solver_container[MESH_0][FLOW_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Flow, config);
				if (euler || navierstokes) cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
			default :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Flow()) {
			case NONE :
				break;
			case PIECEWISE_CONSTANT :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					
					if (config->GetRotating_Frame() == YES)
						solver_container[iMGlevel][FLOW_SOL][SOUR_TERM] = new CSourceRotationalFrame_Flow(nDim, nVar_Flow, config);
					else
						solver_container[iMGlevel][FLOW_SOL][SOUR_TERM] = new CSourceNothing(nDim, nVar_Flow, config);

					solver_container[iMGlevel][FLOW_SOL][CONS_SOUR_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
				}
				
				break;
			default :
				cout << "Source term not implemented." << endl; cin.get();
				break;
		}
		
	}
	
	/*--- Solver definition for the turbulent model problem ---*/
	if (turbulent) {
		
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Turb()) {
			case NONE :
				break;
			case SPACE_UPWIND :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
					if (spalart_allmaras) solver_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
					else if (menter_sst) solver_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
				}
				break;
			default :
				cout << "Convective scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Turb()) {	
			case NONE :
				break;
			case AVG_GRAD :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
					if (spalart_allmaras) solver_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, config);
					else if (menter_sst) solver_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSST(nDim, nVar_Turb, config);
				}
				break;
			case AVG_GRAD_CORRECTED :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
					if (spalart_allmaras) solver_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSA(nDim, nVar_Turb, config);
					else if (menter_sst) solver_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSST(nDim, nVar_Turb, config);
				}
				break;
			case GALERKIN :
				cout << "Viscous scheme not implemented." << endl;
				cin.get(); break;
			default :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Turb()) {
			case NONE :
				break;
			case PIECEWISE_CONSTANT :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					if (spalart_allmaras) solver_container[iMGlevel][TURB_SOL][SOUR_TERM] = new CSourcePieceWise_TurbSA(nDim, nVar_Turb, config);
					else if (menter_sst) solver_container[iMGlevel][TURB_SOL][SOUR_TERM] = new CSourcePieceWise_TurbSST(nDim, nVar_Turb, config);
					solver_container[iMGlevel][TURB_SOL][CONS_SOUR_TERM] = new CSourceNothing(nDim, nVar_Turb, config);
				}
				break;
			default :
				cout << "Source term not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
			if (spalart_allmaras) solver_container[iMGlevel][TURB_SOL][BOUND_TERM] = new CUpwLin_TurbSA(nDim, nVar_Turb, config);
			else if (menter_sst) solver_container[iMGlevel][TURB_SOL][BOUND_TERM] = new CUpwLin_TurbSST(nDim, nVar_Turb, config);
		}
	}
	
	/*--- Solver definition for the multi species plasma model problem ---*/
	if (plasma) {
		
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Plasma()) {
			case NONE :
				break;
			case SPACE_UPWIND :
				if (plasma_monatomic)
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
						solver_container[iMGlevel][PLASMA_SOL][CONV_TERM] = new CUpwRoe_Plasma(nDim, nVar_Plasma, nSpecies, nFluids, config);
						solver_container[iMGlevel][PLASMA_SOL][BOUND_TERM] = new CUpwRoe_Plasma(nDim, nVar_Plasma, nSpecies, nFluids,config);	
					}
				if (plasma_diatomic)
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
						solver_container[iMGlevel][PLASMA_SOL][CONV_TERM] = new CUpwRoe_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
						solver_container[iMGlevel][PLASMA_SOL][BOUND_TERM] = new CUpwRoe_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
					}
				break;
			default :
				cout << "Convective scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Plasma()) {
			case NONE :	
				break;
			case AVG_GRAD :
				if (plasma_monatomic) {
//					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
//						solver_container[iMGlevel][PLASMA_SOL][VISC_TERM] = new CAvgGrad_PlasmaMonatomic(nDim, nVar_Plasma, config);
				}
				if (plasma_diatomic) {
//					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
//						solver_container[iMGlevel][PLASMA_SOL][VISC_TERM] = new CAvgGrad_PlasmaDiatomic(nDim, nVar_Plasma, nDiatomics, nMonatomics, config);
				}
				break;
			default : 
				cout << "Viscous scheme not implemented." << endl; cin.get(); 
				break; 
		}

		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Plasma()) {
			case NONE :
				break;
			case PIECEWISE_CONSTANT :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					if (plasma_monatomic)
						solver_container[iMGlevel][PLASMA_SOL][SOUR_TERM] = new CSourcePieceWise_Plasma(nDim, nVar_Plasma, config);
					if (plasma_diatomic)
						solver_container[iMGlevel][PLASMA_SOL][SOUR_TERM] = new CSourcePieceWise_PlasmaDiatomic(nDim, nVar_Plasma, config);
					solver_container[iMGlevel][PLASMA_SOL][CONS_SOUR_TERM] = new CSourceNothing(nDim, nVar_Plasma, config);
				}
				
				break;
			default :
				cout << "Source term not implemented." << endl; cin.get();
				break;
		}		
	}
	
	/*--- Solver definition for the electric potential problem ---*/
	if (electric) {
		
		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Elec()) {	
			case GALERKIN :
				solver_container[MESH_0][ELEC_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Elec, config);
				break;
			default : cout << "Viscous scheme not implemented." << endl; cin.get(); break;
		}
		
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Elec()) {
			case NONE :
				break;
			case PIECEWISE_CONSTANT :
				solver_container[MESH_0][ELEC_SOL][SOUR_TERM] = new CSourcePieceWise_Elec(nDim, nVar_Elec, config);
				solver_container[MESH_0][ELEC_SOL][CONS_SOUR_TERM] = new CSourceNothing(nDim, nVar_Elec, config);
				cout << " done with electric source residuals" << endl;
				break;
			default :
				cout << "Source term not implemented." << endl; cin.get();
				break;
		}
	}
	
	/*--- Solver definition for the level set model problem ---*/
	if (levelset) {

		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_LevelSet()) {
			case NO_CONVECTIVE : cout << "No convective scheme." << endl; cin.get(); break;
			case SPACE_CENTRED :
				switch (config->GetKind_Centred_LevelSet()) {
					case NO_UPWIND : cout << "No centered scheme." << endl; cin.get(); break;
					default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
				}
				break;
			case SPACE_UPWIND :
				switch (config->GetKind_Upwind_LevelSet()) {
					case NO_UPWIND : cout << "No upwind scheme." << endl; cin.get(); break;
					case ROE_1ST : case ROE_2ND :
						for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
							solver_container[iMGlevel][LEVELSET_SOL][CONV_TERM] = new CUpwLin_LevelSet(nDim, nVar_LevelSet, config);
						}
						break; 
					default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
				}
				break;
			default : cout << "Convective scheme not implemented." << endl; cin.get(); break;
		}
		
		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
			solver_container[iMGlevel][LEVELSET_SOL][BOUND_TERM] = new CUpwLin_LevelSet(nDim, nVar_LevelSet, config);	
		}	
	}
	
	/*--- Solver definition for the level set model problem ---*/
	if (combustion) {
		
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Combustion()) {
			case NONE :
				break;
			case SPACE_UPWIND :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][COMBUSTION_SOL][CONV_TERM] = new CUpwLin_Combustion(nDim, nVar_Combustion, config);
				break;
			default :
				cout << "Convective scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Combustion()) {
			case NONE :
				break;
			case PIECEWISE_CONSTANT :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][COMBUSTION_SOL][SOUR_TERM] = new CSourcePieceWise_Combustion(nDim, nVar_Combustion, config);
				break;
			default :
				cout << "Source term not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
			solver_container[iMGlevel][COMBUSTION_SOL][BOUND_TERM] = new CUpwLin_Combustion(nDim, nVar_Combustion, config);	
		}
		
	}
	
	
	/*--- Solver definition for the flow adjoint problem ---*/
	if ((adj_pot)||(adj_euler)||(adj_ns)) {
		
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_AdjFlow()) {
			case NONE :
				break;
			case SPACE_CENTRED :
				switch (config->GetKind_Centred_AdjFlow()) {
					case LAX : solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config); break;
					case JST : solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentJST_AdjFlow(nDim, nVar_Adj_Flow, config); break;
					default : cout << "Centred scheme not implemented." << endl; cin.get(); break;
				}
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config);
				break;
			case SPACE_UPWIND : 
				switch (config->GetKind_Upwind_AdjFlow()) {	
					case ROE_1ST : solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config); break;
					case ROE_2ND : solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config); break;
					default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
				}
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
				break;
			default :
				cout << "Convective scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_AdjFlow()) {
			case NONE :
				break;
			case AVG_GRAD :
				solver_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
				break;
			case AVG_GRAD_CORRECTED :
				solver_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrected_AdjFlow(nDim, nVar_Adj_Flow, config);
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
				break;
			default :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_AdjFlow()) {
			case NONE :
				break;
			case PIECEWISE_CONSTANT :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					solver_container[iMGlevel][ADJFLOW_SOL][SOUR_TERM] = new CSourcePieceWise_AdjFlow(nDim, nVar_Adj_Flow, config);
					solver_container[iMGlevel][ADJFLOW_SOL][CONS_SOUR_TERM] = new CSourceConservative_AdjFlow(nDim, nVar_Adj_Flow, config);
				}
				break;
			default :
				cout << "Source term not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
			solver_container[iMGlevel][ADJFLOW_SOL][BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
		
	}
	
	/*--- Solver definition for the linearized flow problem ---*/
	if (lin_euler) {
		
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_LinFlow()) {
			case NONE :
				break;
			case SPACE_CENTRED :
				switch (config->GetKind_Centred_LinFlow()) {
					case LAX : solver_container[MESH_0][LINFLOW_SOL][CONV_TERM] = new CCentLax_LinFlow(nDim, nVar_Lin_Flow, config); break;
					case JST : solver_container[MESH_0][LINFLOW_SOL][CONV_TERM] = new CCentJST_LinFlow(nDim, nVar_Lin_Flow, config); break;
					default : cout << "Centred scheme not implemented." << endl; cin.get(); break;
				}
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][LINFLOW_SOL][CONV_TERM] = new CCentLax_LinFlow(nDim, nVar_Lin_Flow, config);
				break;
			default : 
				cout << "Convective scheme not implemented." << endl; cin.get(); 
				break; 
		}
		
		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
			solver_container[iMGlevel][LINFLOW_SOL][BOUND_TERM] = new CCentLax_LinFlow(nDim, nVar_Lin_Flow, config);
	}
	
	/*--- Solver definition for the turbulent adjoint problem ---*/
	if (adj_turb) {
		
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_AdjTurb()) {
			case NONE :
				break;
			case SPACE_CENTRED :
				cout << "Convective scheme not implemented." << endl; cin.get();
				break;
			case SPACE_UPWIND :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJTURB_SOL][CONV_TERM] = new CUpwSca_AdjTurb(nDim, nVar_Adj_Turb, config);
				break;
			default :
				cout << "Convective scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_AdjTurb()) {
			case NONE :
				break;
			case AVG_GRAD :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
			case AVG_GRAD_CORRECTED :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJTURB_SOL][VISC_TERM] = new CAvgGradCorrected_AdjTurb(nDim, nVar_Adj_Turb, config);
				break;
			default :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Turb()) {	
			case NONE :
				break;
			case PIECEWISE_CONSTANT :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					solver_container[iMGlevel][ADJTURB_SOL][SOUR_TERM] = new CSourcePieceWise_AdjTurb(nDim, nVar_Adj_Turb, config);
					solver_container[iMGlevel][ADJTURB_SOL][CONS_SOUR_TERM] = new CSourceConservative_AdjTurb(nDim, nVar_Adj_Turb, config);
				}
				break;
			default :
				cout << "Source term not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
			solver_container[iMGlevel][ADJTURB_SOL][BOUND_TERM] = new CUpwLin_AdjTurb(nDim, nVar_Adj_Turb, config);
	}
	
	/*--- Solver definition for the adjoint electric potential problem ---*/
	if (adj_elec) {
		
		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Elec()) {
			case NONE :
				break;
			case AVG_GRAD :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
			case AVG_GRAD_CORRECTED :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
			case GALERKIN :
				solver_container[MESH_0][ADJELEC_SOL][VISC_TERM] = new CGalerkin_Flow(nDim,nVar_Elec, config);
				break;
			default :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Elec()) {	
			case NONE :
				break;
			case PIECEWISE_CONSTANT :
				cout << "Source term not implemented." << endl; cin.get();
				break;
			case CHARGE_DIST :
				solver_container[MESH_0][ADJELEC_SOL][SOUR_TERM] = new CSourcePieceWise_AdjElec(nDim,nVar_Elec, config);
				break;
			default :
				cout << "Source term not implemented." << endl; cin.get();
				break;
		}
	}
	
	/*--- Solver definition for the linearized electric potential problem ---*/
	if (lin_elec) {
		
		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Elec()) {
			case NONE :
				break;
			case AVG_GRAD :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
			case AVG_GRAD_CORRECTED :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
			case GALERKIN :
				solver_container[MESH_0][LINELEC_SOL][VISC_TERM] = new CGalerkin_Flow(nDim,nVar_Lin_Elec, config);
				break;
			default :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
		}
		
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Elec()) {	
			case NONE :
				break;
			case PIECEWISE_CONSTANT :
				cout << "Source term not implemented." << endl; cin.get();
				break;
			case CHARGE_DIST :
				solver_container[MESH_0][LINELEC_SOL][SOUR_TERM] = new CSourcePieceWise_LinElec(nDim,nVar_Lin_Elec, config);
				break;
			default : cout << "Source term not implemented." << endl; cin.get(); break;
				break;
		}
	}
}

void Geometrical_Deallocation(CGeometry **geometry, CConfig *config) {

	/*--- Deallocation of geometry pointer ---*/
	/*if (config->GetKind_Solver() != NO_SOLVER)
	 for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete geometry[iMGlevel];
	 delete geometry[MESH_0];
	 delete [] geometry;*/
	
}

void Solver_Deallocation(CNumerics ****solver_container, CSolution ***solution_container, CIntegration **integration_container,
											COutput *output, CGeometry **geometry, CConfig *config){
	
	unsigned short iMGlevel;
	bool no_solver, potential, euler, navierstokes, plasma, adj_pot, lin_pot, adj_euler, lin_euler, adj_ns, lin_ns, turbulent,
	adj_turb, lin_turb, electric, adj_elec, lin_elec, spalart_allmaras, menter_sst;
	
	/*--- Initialize some useful booleans ---*/
	no_solver = false;
	potential = false;	euler = false;		navierstokes = false;	turbulent = false;	electric = false;	plasma = false;
	adj_pot = false;	adj_euler = false;	adj_ns = false;			adj_turb = false;	adj_elec = false; 	spalart_allmaras = false;
	lin_pot = false;	lin_euler = false;	lin_ns = false;			lin_turb = false;	lin_elec = false; 	menter_sst = false;
	
	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
		case NO_SOLVER: no_solver = true; break;
		case POTENTIAL_FLOW: potential = true; break;
		case EULER : euler = true; break;
		case NAVIER_STOKES: navierstokes = true; break;
		case RANS : navierstokes = true; turbulent = true; break;
		case MULTI_SPECIES_NAVIER_STOKES : plasma = true; break;//electric = true; break;
		case ELECTRIC_POTENTIAL: electric = true; break;
		case ADJ_EULER : case ONE_SHOT_ADJ_EULER : euler = true; adj_euler = true; break;
		case ADJ_NAVIER_STOKES : case ONE_SHOT_ADJ_NAVIER_STOKES : navierstokes = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
		case ADJ_RANS : navierstokes = true; turbulent = true; adj_ns = true; adj_turb = true; break;
		case ADJ_ELECTRIC_POTENTIAL: electric = true; adj_elec = true; break;
		case LIN_EULER: euler = true; lin_euler = true; break;
		case LIN_ELECTRIC_POTENTIAL: electric = true; lin_elec = true; break;
	}
	
	/*--- Assign turbulence model booleans --- */
	if (turbulent)
		switch (config->GetKind_Turb_Model()){
			case SA: case SA_COMP: spalart_allmaras = true; break;
			case SST: menter_sst = true; break;
		}
	
	/*--- Deallocation of integration_container---*/
	if (no_solver) {
		delete integration_container[0];
		delete integration_container[1];
		delete integration_container[2];
		delete integration_container[3];
		delete integration_container[4];
		return;
	}
	if (potential)		delete integration_container[FLOW_SOL];
	if (euler)          delete integration_container[FLOW_SOL];
	if (navierstokes)   delete integration_container[FLOW_SOL];
	if (turbulent) {    delete integration_container[TURB_SOL];    delete integration_container[EIKONAL_SOL];	}
	if (electric) 	    delete integration_container[ELEC_SOL];
	if (plasma) 	    delete integration_container[PLASMA_SOL];
	if (adj_euler)      delete integration_container[ADJFLOW_SOL];
	if (adj_ns) 	    delete integration_container[ADJFLOW_SOL];
	if (adj_turb) 	    delete integration_container[ADJTURB_SOL];
	if (adj_elec) 	    delete integration_container[ADJELEC_SOL];
	if (lin_euler) 	    delete integration_container[LINFLOW_SOL];
	if (lin_elec) 	    delete integration_container[LINELEC_SOL];
	
	delete [] integration_container;
	
	/*--- Deallocation of solution_container---*/
	for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
		//if (potential) 		delete solution_container[iMGlevel][FLOW_SOL];
		if (euler) 		    delete solution_container[iMGlevel][FLOW_SOL];
		//if (navierstokes)   delete solution_container[iMGlevel][FLOW_SOL];
		//if (turbulent) {
		//	if (spalart_allmaras) delete solution_container[iMGlevel][TURB_SOL];
		//	else if (menter_sst)  delete solution_container[iMGlevel][TURB_SOL];
		//	delete solution_container[iMGlevel][EIKONAL_SOL];
		//}
		//if (electric) 		delete solution_container[iMGlevel][ELEC_SOL];
		//if (plasma) 		delete solution_container[iMGlevel][PLASMA_SOL];
		//if (adj_euler) 		delete solution_container[iMGlevel][ADJFLOW_SOL];
		//if (adj_ns) 		delete solution_container[iMGlevel][ADJFLOW_SOL];
		//if (adj_turb)		delete solution_container[iMGlevel][ADJTURB_SOL];
		//if (adj_elec) 		delete solution_container[iMGlevel][ADJELEC_SOL];
		//if (lin_euler) 		delete solution_container[iMGlevel][LINFLOW_SOL];
		//if (lin_elec) 		delete solution_container[iMGlevel][LINELEC_SOL];
		delete [] solution_container[iMGlevel];
	}
	delete [] solution_container;
	
	
	/*--- Deallocation of solver_container---*/
	/*--- Deallocation for the Potential, Euler, Navier-Stokes problems ---*/
	if ((potential) || (euler) || (navierstokes)) {
		
		//--- Deallocation of the convective scheme for each equation and mesh level ---
		switch (config->GetKind_ConvNumScheme_Flow()) {
			case NONE :
				break;
			case SPACE_CENTRED :
				switch (config->GetKind_Centred_Flow()) {
					case LAX : delete solver_container[MESH_0][FLOW_SOL][CONV_TERM]; break;
					case JST : delete solver_container[MESH_0][FLOW_SOL][CONV_TERM]; break;
				}
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					delete solver_container[iMGlevel][FLOW_SOL][CONV_TERM];
				break;
			case SPACE_UPWIND :
				switch (config->GetKind_Upwind_Flow()) {
					case ROE_1ST : case ROE_2ND : case AUSM_1ST : case AUSM_2ND : case HLLC_1ST : case HLLC_2ND :
						delete solver_container[MESH_0][FLOW_SOL][CONV_TERM]; break;
				}
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					delete solver_container[iMGlevel][FLOW_SOL][CONV_TERM];
				break;
		} /*
			 
			 //--- Deallocation of the viscous scheme for each equation and mesh level ---
			 switch (config->GetKind_ViscNumScheme_Flow()) {
			 case NONE :
			 break;
			 case DIVERGENCE_THEOREM :
			 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
			 delete solver_container[iMGlevel][FLOW_SOL][VISC_TERM];
			 break;
			 case DIVERGENCE_THEOREM_WEISS :
			 delete solver_container[MESH_0][FLOW_SOL][VISC_TERM];
			 for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
			 delete solver_container[iMGlevel][FLOW_SOL][VISC_TERM];
			 break;
			 case GALERKIN :
			 if (potential) delete solver_container[MESH_0][FLOW_SOL][VISC_TERM];
			 break;
			 }
			 
			 //--- Deallocation of the source term integration scheme for each equation and mesh level ---
			 switch (config->GetKind_SourNumScheme_Flow()) {
			 case NONE :
			 break;
			 case PIECEWISE_CONSTANT :
			 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
			 delete solver_container[iMGlevel][FLOW_SOL][SOUR_TERM];
			 delete solver_container[iMGlevel][FLOW_SOL][CONS_SOUR_TERM];
			 }
			 break;
			 }
			 */
		//--- Deallocation of the boundary condition method ---
		//for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
		//	delete solver_container[iMGlevel][FLOW_SOL][BOUND_TERM];
	}
	/*
	 //--- Solver definition for the turbulent model problem ---
	 if (turbulent) {
	 
	 //--- Definition of the convective scheme for each equation and mesh level ---
	 switch (config->GetKind_ConvNumScheme_Turb()) {
	 case NONE :
	 break;
	 case SPACE_UPWIND :
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
	 if (spalart_allmaras) delete solver_container[iMGlevel][TURB_SOL][CONV_TERM];
	 else if (menter_sst)  delete solver_container[iMGlevel][TURB_SOL][CONV_TERM];
	 }
	 break;
	 }
	 
	 //--- Definition of the viscous scheme for each equation and mesh level ---
	 switch (config->GetKind_ViscNumScheme_Turb()) {
	 case NONE :
	 break;
	 case DIVERGENCE_THEOREM :
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
	 if (spalart_allmaras) delete solver_container[iMGlevel][TURB_SOL][VISC_TERM];
	 else if (menter_sst)  delete solver_container[iMGlevel][TURB_SOL][VISC_TERM];
	 }
	 break;
	 case DIVERGENCE_THEOREM_WEISS :
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
	 if (spalart_allmaras) delete solver_container[iMGlevel][TURB_SOL][VISC_TERM];
	 else if (menter_sst)  delete solver_container[iMGlevel][TURB_SOL][VISC_TERM];
	 }
	 break;
	 }
	 
	 //--- Definition of the source term integration scheme for each equation and mesh level ---
	 switch (config->GetKind_SourNumScheme_Turb()) {
	 case NONE :
	 break;
	 case PIECEWISE_CONSTANT :
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
	 if (spalart_allmaras) delete solver_container[iMGlevel][TURB_SOL][SOUR_TERM];
	 else if (menter_sst)  delete solver_container[iMGlevel][TURB_SOL][SOUR_TERM];
	 delete solver_container[iMGlevel][TURB_SOL][CONS_SOUR_TERM];
	 }
	 break;
	 }
	 
	 //--- Definition of the boundary condition method ---
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
	 if (spalart_allmaras) delete solver_container[iMGlevel][TURB_SOL][BOUND_TERM];
	 else if (menter_sst)  delete solver_container[iMGlevel][TURB_SOL][BOUND_TERM];
	 }
	 }
	 
	 //--- Solver definition for the multi species plasma model problem ---
	 if (plasma) {
	 
	 //--- Definition of the convective scheme for each equation and mesh level ---
	 switch (config->GetKind_ConvNumScheme_Plasma()) {
	 case NONE :
	 break;
	 case SPACE_UPWIND :
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
	 delete solver_container[iMGlevel][PLASMA_SOL][CONV_TERM];
	 }
	 break;
	 default :
	 cout << "Convective scheme not implemented." << endl; cin.get();
	 break;
	 }
	 
	 //--- Definition of the viscous scheme for each equation and mesh level ---
	 //This was commented.
	 //	switch (config->GetKind_ViscNumScheme_Plasma()) {
	 //		case NONE :
	 //			break;
	 //		case DIVERGENCE_THEOREM :
	 //			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 //				solver_container[iMGlevel][PLASMA_SOL][VISC_TERM] = new CDivergence_Plasma(nDim, nVar_Plasma, config);
	 //			break;
	 //		default :
	 //			cout << "Viscous scheme not implemented." << endl; cin.get();
	 //			break;
	 //	}
	 
	 //--- Definition of the source term integration scheme for each equation and mesh level ---
	 
	 //	switch (config->GetKind_SourNumScheme_Plasma()) {
	 //		case NONE :
	 //			break;
	 //		case PIECEWISE_CONSTANT :
	 //			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
	 //				solver_container[iMGlevel][PLASMA_SOL][SOUR_TERM] = new CSourcePieceWise_Plasma(nDim, nVar_Plasma, config);
	 //				solver_container[iMGlevel][PLASMA_SOL][CONS_SOUR_TERM] = new CSourceNothing(nDim, nVar_Plasma);
	 //			}
	 //			break;
	 //		default :
	 //			cout << "Source term not implemented." << endl; cin.get();
	 //			break;
	 //	}
	 
	 //Until here.
	 //--- Definition of the boundary condition method ---
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
	 delete solver_container[iMGlevel][PLASMA_SOL][BOUND_TERM];
	 }
	 }
	 
	 //--- Solver definition for the electric potential problem ---
	 if (electric) {
	 
	 //--- Definition of the viscous scheme for each equation and mesh level ---
	 switch (config->GetKind_ViscNumScheme_Elec()) {
	 case GALERKIN :
	 delete solver_container[MESH_0][ELEC_SOL][VISC_TERM];
	 break;
	 default : cout << "Viscous scheme not implemented." << endl; cin.get(); break;
	 }
	 
	 //--- Definition of the source term integration scheme for each equation and mesh level ---
	 switch (config->GetKind_SourNumScheme_Elec()) {
	 case NONE :
	 break;
	 case PIECEWISE_CONSTANT :
	 delete solver_container[MESH_0][ELEC_SOL][SOUR_TERM];
	 delete solver_container[MESH_0][ELEC_SOL][CONS_SOUR_TERM];
	 break;
	 default :
	 cout << "Source term not implemented." << endl; cin.get();
	 break;
	 }
	 }
	 
	 //--- Solver definition for the flow adjoint problem ---
	 if ((adj_pot)||(adj_euler)||(adj_ns)) {
	 
	 //--- Definition of the convective scheme for each equation and mesh level ---
	 switch (config->GetKind_ConvNumScheme_AdjFlow()) {
	 case NONE :
	 break;
	 case SPACE_CENTRED :
	 switch (config->GetKind_Centred_AdjFlow()) {
	 case LAX : delete solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM]; break;
	 case JST : delete solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM]; break;
	 default : cout << "Centred scheme not implemented." << endl; cin.get(); break;
	 }
	 for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete solver_container[iMGlevel][ADJFLOW_SOL][CONV_TERM];
	 break;
	 case SPACE_UPWIND :
	 switch (config->GetKind_Upwind_AdjFlow()) {
	 case ROE_1ST : delete solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM]; break;
	 case ROE_2ND : delete solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM]; break;
	 default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
	 }
	 for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete solver_container[iMGlevel][ADJFLOW_SOL][CONV_TERM];
	 break;
	 default :
	 cout << "Convective scheme not implemented." << endl; cin.get();
	 break;
	 }
	 
	 //--- Definition of the viscous scheme for each equation and mesh level ---
	 switch (config->GetKind_ViscNumScheme_AdjFlow()) {
	 case NONE :
	 break;
	 case DIVERGENCE_THEOREM :
	 delete solver_container[MESH_0][ADJFLOW_SOL][VISC_TERM];
	 for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete solver_container[iMGlevel][ADJFLOW_SOL][VISC_TERM];
	 break;
	 case DIVERGENCE_THEOREM_WEISS :
	 delete solver_container[MESH_0][ADJFLOW_SOL][VISC_TERM];
	 for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete solver_container[iMGlevel][ADJFLOW_SOL][VISC_TERM];
	 break;
	 default :
	 cout << "Viscous scheme not implemented." << endl; cin.get();
	 break;
	 }
	 
	 //--- Definition of the source term integration scheme for each equation and mesh level ---
	 switch (config->GetKind_SourNumScheme_AdjFlow()) {
	 case NONE :
	 break;
	 case PIECEWISE_CONSTANT :
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
	 delete solver_container[iMGlevel][ADJFLOW_SOL][SOUR_TERM];
	 delete solver_container[iMGlevel][ADJFLOW_SOL][CONS_SOUR_TERM];
	 }
	 break;
	 default :
	 cout << "Source term not implemented." << endl; cin.get();
	 break;
	 }
	 
	 //--- Definition of the boundary condition method ---
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete solver_container[iMGlevel][ADJFLOW_SOL][BOUND_TERM];
	 
	 }
	 
	 //--- Solver definition for the linearized flow problem ---
	 if (lin_euler) {
	 
	 //--- Definition of the convective scheme for each equation and mesh level ---
	 switch (config->GetKind_ConvNumScheme_LinFlow()) {
	 case NONE :
	 break;
	 case SPACE_CENTRED :
	 switch (config->GetKind_Centred_LinFlow()) {
	 case LAX : delete solver_container[MESH_0][LINFLOW_SOL][CONV_TERM]; break;
	 case JST : delete solver_container[MESH_0][LINFLOW_SOL][CONV_TERM]; break;
	 default : cout << "Centred scheme not implemented." << endl; cin.get(); break;
	 }
	 for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete solver_container[iMGlevel][LINFLOW_SOL][CONV_TERM];
	 break;
	 default :
	 cout << "Convective scheme not implemented." << endl; cin.get();
	 break;
	 }
	 
	 //--- Definition of the boundary condition method ---
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete solver_container[iMGlevel][LINFLOW_SOL][BOUND_TERM];
	 }
	 
	 //--- Solver definition for the turbulent adjoint problem ---
	 if (adj_turb) {
	 
	 //--- Definition of the convective scheme for each equation and mesh level ---
	 switch (config->GetKind_ConvNumScheme_AdjTurb()) {
	 case NONE :
	 break;
	 case SPACE_CENTRED :
	 cout << "Convective scheme not implemented." << endl; cin.get();
	 break;
	 case SPACE_UPWIND :
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete solver_container[iMGlevel][ADJTURB_SOL][CONV_TERM];
	 break;
	 default :
	 cout << "Convective scheme not implemented." << endl; cin.get();
	 break;
	 }
	 
	 //--- Definition of the viscous scheme for each equation and mesh level ---
	 switch (config->GetKind_ViscNumScheme_AdjTurb()) {
	 case NONE :
	 break;
	 case DIVERGENCE_THEOREM :
	 cout << "Viscous scheme not implemented." << endl; cin.get();
	 break;
	 case DIVERGENCE_THEOREM_WEISS :
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete solver_container[iMGlevel][ADJTURB_SOL][VISC_TERM];
	 break;
	 default :
	 cout << "Viscous scheme not implemented." << endl; cin.get();
	 break;
	 }
	 
	 //--- Definition of the source term integration scheme for each equation and mesh level ---
	 switch (config->GetKind_SourNumScheme_Turb()) {
	 case NONE :
	 break;
	 case PIECEWISE_CONSTANT :
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
	 delete solver_container[iMGlevel][ADJTURB_SOL][SOUR_TERM];
	 delete solver_container[iMGlevel][ADJTURB_SOL][CONS_SOUR_TERM];
	 }
	 break;
	 default :
	 cout << "Source term not implemented." << endl; cin.get();
	 break;
	 }
	 
	 //--- Definition of the boundary condition method ---
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete solver_container[iMGlevel][ADJTURB_SOL][BOUND_TERM];
	 }
	 
	 //--- Solver definition for the adjoint electric potential problem ---
	 if (adj_elec) {
	 
	 //--- Definition of the viscous scheme for each equation and mesh level ---
	 switch (config->GetKind_ViscNumScheme_Elec()) {
	 case NONE :
	 break;
	 case DIVERGENCE_THEOREM :
	 cout << "Viscous scheme not implemented." << endl; cin.get();
	 break;
	 case DIVERGENCE_THEOREM_WEISS :
	 cout << "Viscous scheme not implemented." << endl; cin.get();
	 break;
	 case GALERKIN :
	 delete solver_container[MESH_0][ADJELEC_SOL][VISC_TERM];
	 break;
	 default :
	 cout << "Viscous scheme not implemented." << endl; cin.get();
	 break;
	 }
	 
	 //--- Definition of the source term integration scheme for each equation and mesh level ---
	 switch (config->GetKind_SourNumScheme_Elec()) {
	 case NONE :
	 break;
	 case PIECEWISE_CONSTANT :
	 cout << "Source term not implemented." << endl; cin.get();
	 break;
	 case CHARGE_DIST :
	 delete solver_container[MESH_0][ADJELEC_SOL][SOUR_TERM];
	 break;
	 default :
	 cout << "Source term not implemented." << endl; cin.get();
	 break;
	 }
	 }
	 
	 //--- Solver definition for the linearized electric potential problem ---
	 if (lin_elec) {
	 
	 //--- Definition of the viscous scheme for each equation and mesh level ---
	 switch (config->GetKind_ViscNumScheme_Elec()) {
	 case NONE :
	 break;
	 case DIVERGENCE_THEOREM :
	 cout << "Viscous scheme not implemented." << endl; cin.get();
	 break;
	 case DIVERGENCE_THEOREM_WEISS :
	 cout << "Viscous scheme not implemented." << endl; cin.get();
	 break;
	 case GALERKIN :
	 delete solver_container[MESH_0][LINELEC_SOL][VISC_TERM];
	 break;
	 default :
	 cout << "Viscous scheme not implemented." << endl; cin.get();
	 break;
	 }
	 
	 //--- Definition of the source term integration scheme for each equation and mesh level ---
	 switch (config->GetKind_SourNumScheme_Elec()) {
	 case NONE :
	 break;
	 case PIECEWISE_CONSTANT :
	 cout << "Source term not implemented." << endl; cin.get();
	 break;
	 case CHARGE_DIST :
	 delete solver_container[MESH_0][LINELEC_SOL][SOUR_TERM];
	 break;
	 default : cout << "Source term not implemented." << endl; cin.get(); break;
	 break;
	 }
	 }
	 
	 
	 //--- Definition of the Class for the numerical method: solver_container[MESH_LEVEL][EQUATION][EQ_TERM] ---
	 for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
	 solver_container[iMGlevel] = new CNumerics** [MAX_SOLS];
	 for (iSol = 0; iSol < MAX_SOLS; iSol++)
	 solver_container[iMGlevel][iSol] = new CNumerics* [MAX_TERMS];
	 }*/
	
	/*--- Deallocation of output pointer ---*/
	delete output;
	
	/*--- Deallocation of config pointer ---*/
	delete config;
	
	cout << "Deallocation completed." << endl;
}

