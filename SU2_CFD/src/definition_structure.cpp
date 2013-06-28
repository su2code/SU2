/*!
 * \file definition_structure.cpp
 * \brief Main subroutines used by SU2_CFD.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.1.
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


unsigned short GetnDomain(CConfig *config, string val_mesh_filename, unsigned short val_format) {

	string text_line, Marker_Tag;
	ifstream mesh_file;
	short nDomain = 1;
	bool isFound = false;
	char cstr[200];
	string::size_type position;
	int rank = MASTER_NODE;

#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();	
#endif

	if (rank == MASTER_NODE)
		cout << endl <<"---------------------- Read grid file information -----------------------" << endl;

	switch (val_format) {
	case SU2:

		/*--- Open grid file ---*/
		strcpy (cstr, val_mesh_filename.c_str());
		mesh_file.open(cstr, ios::in);
		if (mesh_file.fail()) {
			cout << "There is no geometry file!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
#ifdef NO_MPI
			exit(1);
#else
			MPI::COMM_WORLD.Abort(1);
			MPI::Finalize();
#endif
		}
		/*--- Read SU2 mesh file ---*/
		while (getline (mesh_file,text_line)) {
			/*--- Search for the "NDOM" keyword to see if there are multiple domains ---*/
			position = text_line.find ("NDOM=",0);
			if (position != string::npos) {
				text_line.erase (0,5); nDomain = atoi(text_line.c_str()); isFound = true;
				if (rank == MASTER_NODE) {
					if (nDomain == 1) cout << "SU2 mesh file format with a single domain." << endl;
					else if (nDomain >  1) cout << "SU2 mesh file format with " << nDomain << " domains." << endl;
					else if (nDomain <= 0) {
						cout << "Error: Number of mesh domains is less than 1 !!!" << endl;
						cout << "Press any key to exit..." << endl;
						cin.get();
#ifdef NO_MPI
						exit(1);
#else
						MPI::COMM_WORLD.Abort(1);
						MPI::Finalize();
#endif
					}
				}
			}
		}
		/*--- If the "NDOM" keyword was not found, assume this is an ordinary
       simulation on a single domain ---*/
		if (!isFound) {
			nDomain = 1;
			if (rank == MASTER_NODE) cout << "SU2 mesh file format with a single domain." << endl;
		}
		break;

	case CGNS:

		nDomain = 1;
		if (rank == MASTER_NODE) cout << "CGNS mesh file format with a single domain." << endl;
		break;

	case NETCDF_ASCII:

		nDomain = 1;
		if (rank == MASTER_NODE) cout << "NETCDF mesh file format with a single domain." << endl;
		break;

	}

	return (unsigned short) nDomain;
}



void Geometrical_Definition(CGeometry ***geometry, CConfig *config, unsigned short val_nDomain) {
	unsigned short iMGlevel, iDomain;
	int rank = MASTER_NODE;

#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	if (rank == MASTER_NODE)
		cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;

	for (iDomain = 0; iDomain < val_nDomain; iDomain++) {

		if ((rank == MASTER_NODE) && val_nDomain > 1)
			cout << "  Domain " << iDomain+1 << ":" << endl;

		/*--- Single grid geometry structure generation ---*/

		/*--- Compute elements surrounding points, points surrounding points, and elements surrounding elements ---*/
		if (rank == MASTER_NODE) cout << "Setting local point and element connectivity." <<endl; 
		geometry[iDomain][MESH_0]->SetEsuP(); 
		geometry[iDomain][MESH_0]->SetPsuP(); 
		geometry[iDomain][MESH_0]->SetEsuE();

		/*--- Check the orientation before computing geometrical quantities ---*/
		if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." <<endl; 
		geometry[iDomain][MESH_0]->SetBoundVolume();
		geometry[iDomain][MESH_0]->Check_Orientation(config);

		/*--- Create the edge structure ---*/
		if (rank == MASTER_NODE) cout << "Identifying edges and vertices." <<endl; 
		geometry[iDomain][MESH_0]->SetEdges();
		geometry[iDomain][MESH_0]->SetVertex(config);

		/*--- Compute center of gravity ---*/
		if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl; 
		geometry[iDomain][MESH_0]->SetCG();

		/*--- Create the control volume structures ---*/
		if (rank == MASTER_NODE) cout << "Setting the control volume structure." << endl; 
		geometry[iDomain][MESH_0]->SetControlVolume(config, ALLOCATE);
		geometry[iDomain][MESH_0]->SetBoundControlVolume(config, ALLOCATE);		

		/*--- Identify closest normal neighbor ---*/
		if (rank == MASTER_NODE) cout << "Searching for closest normal neighbor on the surface." << endl; 
		geometry[iDomain][MESH_0]->FindClosestNeighbor(config);

		/*--- Find any sharp edges ---*/
		if (rank == MASTER_NODE) cout << "Searching for sharp corners on the geometry." << endl; 
		geometry[iDomain][MESH_0]->FindSharpEdges(config);

		/*--- For a rotating frame, set the velocity due to rotation at each mesh point ---*/
		if (config->GetRotating_Frame()) geometry[iDomain][MESH_0]->SetRotationalVelocity(config);
		if ((config->GetMGLevels() != 0) && (rank == MASTER_NODE)) cout << "Setting the multigrid structure." <<endl; 

#ifndef NO_MPI
		/*--- Synchronization point after the multigrid stuff ---*/
		MPI::COMM_WORLD.Barrier();
#endif	

		/*--- Loop over all the new grid ---*/
		for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++) {

			/*--- Create main aglomeration structure (ingluding MPI stuff) ---*/
			geometry[iDomain][iMGlevel] = new CMultiGridGeometry(geometry[iDomain][iMGlevel-1], config, iMGlevel);

			/*--- Compute element surronding points, point surronding point, and element surronding elements ---*/
			geometry[iDomain][iMGlevel]->SetPsuP(geometry[iDomain][iMGlevel-1]);

			/*--- Create the edge structure ---*/
			geometry[iDomain][iMGlevel]->SetEdges();
			geometry[iDomain][iMGlevel]->SetVertex(geometry[iDomain][iMGlevel-1], config);

			/*--- Create the control volume structures ---*/
			geometry[iDomain][iMGlevel]->SetControlVolume(config,geometry[iDomain][iMGlevel-1], ALLOCATE);
			geometry[iDomain][iMGlevel]->SetBoundControlVolume(config,geometry[iDomain][iMGlevel-1], ALLOCATE);
			geometry[iDomain][iMGlevel]->SetCoord(geometry[iDomain][iMGlevel-1]);

			/*--- Find closest neighbor to a surface point ---*/
			geometry[iDomain][iMGlevel]->FindClosestNeighbor(config);

			/*--- For a rotating frame, set the velocity due to rotation at each mesh point ---*/
			if (config->GetRotating_Frame()) geometry[iDomain][iMGlevel]->SetRotationalVelocity(config);

		}

	}
}

void Solution_Definition(CSolution ***solution_container, CGeometry **geometry, CConfig *config) {

	unsigned short iMGlevel;
	bool potential, euler, navierstokes, combustion, plasma,
	plasma_monatomic, plasma_diatomic, levelset, adj_pot, lin_pot, adj_euler,
	lin_euler, adj_ns, lin_ns, turbulent, adj_turb, lin_turb, electric, wave, adj_levelset, adj_plasma,
	spalart_allmaras, menter_sst, template_solver;

	/*--- Initialize some useful booleans ---*/
	potential = false;	euler = false;		navierstokes = false;	combustion = false; turbulent = false;	electric = false;	plasma_monatomic = false;
	plasma_diatomic = false; levelset = false; plasma = false;
	adj_pot = false;	adj_euler = false;	adj_ns = false;			adj_turb = false;	wave = false; 	adj_levelset = false;	 spalart_allmaras = false;
	lin_pot = false;	lin_euler = false;	lin_ns = false;			lin_turb = false;	menter_sst = false; adj_plasma = false;
	template_solver = false;

	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
	case TEMPLATE_SOLVER: template_solver = true; break;
	case EULER : euler = true; break;
	case NAVIER_STOKES: navierstokes = true; break;
	case RANS : navierstokes = true; turbulent = true; break;
	case FREE_SURF_EULER: euler = true; levelset = true; break;
	case FREE_SURF_NAVIER_STOKES: navierstokes = true; levelset = true; break;
	case FREE_SURF_RANS: navierstokes = true; turbulent = true; levelset = true; break;
	case NS_PLASMA : plasma = true; break;
	case ELECTRIC_POTENTIAL: electric = true; break;
	case WAVE: wave = true; break;
	case ADJ_EULER : euler = true; adj_euler = true; break;
	case ADJ_NAVIER_STOKES : navierstokes = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
	case ADJ_RANS : navierstokes = true; turbulent = true; adj_ns = true; adj_turb = true; break;
	case ADJ_FREE_SURF_EULER: euler = true; adj_euler = true; levelset = true; adj_levelset = true; break;
	case ADJ_FREE_SURF_NAVIER_STOKES: navierstokes = true; adj_ns = true; levelset = true; adj_levelset = true; break;
	case ADJ_FREE_SURF_RANS: navierstokes = true; adj_ns = true; turbulent = true; adj_turb = true; levelset = true; adj_levelset = true; break;
	case ADJ_NS_PLASMA : plasma = true; adj_plasma = true; break;
	case LIN_EULER: euler = true; lin_euler = true; break;
	}

	/*--- Assign turbulence model booleans --- */
	if (turbulent)
		switch (config->GetKind_Turb_Model()){
		case SA: case SA_COMP: spalart_allmaras = true; break;
		case SST: menter_sst = true; break;
		default: cout << "Specified turbulence model unavailable or none selected" << endl; cin.get(); break;
		}

	if (plasma) {
		switch (config->GetKind_GasModel()){
		case AIR7: plasma_diatomic = true; break;
		case ARGON: plasma_monatomic = true; break;
		case AIR21: plasma_diatomic = true; break;
		default: cout << "Specified plasma model unavailable or none selected" << endl; cin.get(); break;
		}
		if (config->GetElectricSolver()) electric  = true;
	}


	/*--- Definition of the Class for the solution: solution_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, navierstokes
   and potential are incompatible, they use the same position in sol container ---*/
	for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {

		/*--- Allocate solution for a template problem ---*/
		if (template_solver) {
			solution_container[iMGlevel][TEMPLATE_SOL] = new CTemplateSolution(geometry[iMGlevel], config);
		}

		/*--- Allocate solution for direct problem ---*/
		if (potential) {
			solution_container[iMGlevel][FLOW_SOL] = new CPotentialSolution(geometry[iMGlevel], config);
		}
		if (euler) {
			solution_container[iMGlevel][FLOW_SOL] = new CEulerSolution(geometry[iMGlevel], config);
		}
		if (navierstokes) {
			solution_container[iMGlevel][FLOW_SOL] = new CNSSolution(geometry[iMGlevel], config);
		}
		if (turbulent) {
			if (spalart_allmaras) solution_container[iMGlevel][TURB_SOL] = new CTurbSASolution(geometry[iMGlevel], config);
			else if (menter_sst) solution_container[iMGlevel][TURB_SOL] = new CTurbSSTSolution(geometry[iMGlevel], config);
		}
		if (electric) {
			solution_container[iMGlevel][ELEC_SOL] = new CElectricSolution(geometry[iMGlevel], config);
		}
		if (plasma) {
			solution_container[iMGlevel][PLASMA_SOL] = new CPlasmaSolution(geometry[iMGlevel], config);
		}
		if (levelset) {
			solution_container[iMGlevel][LEVELSET_SOL] = new CLevelSetSolution(geometry[iMGlevel], config);
		}		
		if (wave) {
			solution_container[iMGlevel][WAVE_SOL] = new CWaveSolution(geometry[iMGlevel], config);
		}

		/*--- Allocate solution for adjoint problem ---*/
		if (adj_pot) {
			cout <<"Equation not implemented." << endl; cin.get(); break;
		}
		if (adj_euler) {
			solution_container[iMGlevel][ADJFLOW_SOL] = new CAdjEulerSolution(geometry[iMGlevel], config);
		}
		if (adj_ns) {
			solution_container[iMGlevel][ADJFLOW_SOL] = new CAdjNSSolution(geometry[iMGlevel], config);
		}
		if (adj_turb) {
			solution_container[iMGlevel][ADJTURB_SOL] = new CAdjTurbSolution(geometry[iMGlevel], config);
		}
		if (adj_levelset) {
			solution_container[iMGlevel][ADJLEVELSET_SOL] = new CAdjLevelSetSolution(geometry[iMGlevel], config);
		}
		if (adj_plasma) {
			solution_container[iMGlevel][ADJPLASMA_SOL] = new CAdjPlasmaSolution(geometry[iMGlevel], config);
		}

		/*--- Allocate solution for linear problem (at the moment we use the same scheme as the adjoint problem) ---*/
		if (lin_pot) {
			cout <<"Equation not implemented." << endl; cin.get(); break;
		}
		if (lin_euler) {
			solution_container[iMGlevel][LINFLOW_SOL] = new CLinEulerSolution(geometry[iMGlevel], config);
		}
		if (lin_ns) {
			cout <<"Equation not implemented." << endl; cin.get(); break;
		}

	}

}


void Integration_Definition(CIntegration **integration_container, CGeometry **geometry, CConfig *config) {

	bool potential, euler, navierstokes, combustion, plasma, plasma_monatomic, plasma_diatomic, levelset, adj_plasma, adj_pot, lin_pot, adj_euler, lin_euler, adj_ns, lin_ns, turbulent,
	adj_turb, lin_turb, electric, wave, spalart_allmaras, menter_sst, template_solver, adj_levelset;

	/*--- Initialize some useful booleans ---*/
	potential = false;	euler = false;		navierstokes = false;	combustion = false; turbulent = false;	electric = false;	plasma_monatomic = false;
	plasma_diatomic = false; levelset = false; plasma = false;
	adj_pot = false;	adj_euler = false;	adj_ns = false;			adj_turb = false;	wave = false; adj_levelset = false;	spalart_allmaras = false;
	lin_pot = false;	lin_euler = false;	lin_ns = false;			lin_turb = false;	menter_sst = false; adj_plasma = false;
	template_solver = false;

	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
	case TEMPLATE_SOLVER: template_solver = true; break;
	case EULER : euler = true; break;
	case NAVIER_STOKES: navierstokes = true; break;
	case FREE_SURF_EULER: euler = true; levelset = true; break;
	case FREE_SURF_NAVIER_STOKES: navierstokes = true; levelset = true; break;
	case RANS : navierstokes = true; turbulent = true; break;
	case NS_PLASMA : plasma = true; break;
	case ELECTRIC_POTENTIAL: electric = true; break;
	case WAVE: wave = true; break;
	case ADJ_EULER : euler = true; adj_euler = true; break;
	case ADJ_NAVIER_STOKES : navierstokes = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
	case ADJ_RANS : navierstokes = true; turbulent = true; adj_ns = true; adj_turb = true; break;
	case ADJ_NS_PLASMA : plasma = true; adj_plasma = true; break;
	case ADJ_FREE_SURF_EULER: euler = true; levelset = true; adj_euler = true; adj_levelset = true; break;
	case ADJ_FREE_SURF_NAVIER_STOKES: navierstokes = true; levelset = true; adj_ns = true; adj_levelset = true; break;
	case LIN_EULER: euler = true; lin_euler = true; break;
	}

	/*--- Assign turbulence model booleans --- */
	if (turbulent)
		switch (config->GetKind_Turb_Model()){
		case SA: case SA_COMP: spalart_allmaras = true; break;
		case SST: menter_sst = true; break;
		default: cout << "Specified turbulence model unavailable or none selected" << endl; cin.get(); break;
		}

	if (plasma) {
		switch (config->GetKind_GasModel()){
		case AIR7: plasma_diatomic = true; break;
		case ARGON: plasma_monatomic = true; break;
		case AIR21: plasma_diatomic = true; break;
		default: cout << "Specified plasma model unavailable or none selected" << endl; cin.get(); break;
		}
		if (config->GetElectricSolver()) electric  = true;
	}

	/*--- Allocate solution for a template problem ---*/
	if (template_solver) integration_container[TEMPLATE_SOL] = new CSingleGridIntegration(config);

	/*--- Allocate solution for direct problem ---*/
	if (potential) integration_container[FLOW_SOL] = new CPotentialIntegration(config);
	if (euler) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
	if (navierstokes) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
	if (turbulent) integration_container[TURB_SOL] = new CSingleGridIntegration(config);
	if (electric) integration_container[ELEC_SOL] = new CPotentialIntegration(config);
	if (plasma) integration_container[PLASMA_SOL] = new CMultiGridIntegration(config);
	if (levelset) integration_container[LEVELSET_SOL] = new CSingleGridIntegration(config);
	if (wave) integration_container[WAVE_SOL] = new CSingleGridIntegration(config);

	/*--- Allocate solution for adjoint problem ---*/
	if (adj_pot) { cout <<"Equation not implemented." << endl; cin.get(); }
	if (adj_euler) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
	if (adj_ns) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
	if (adj_turb) integration_container[ADJTURB_SOL] = new CSingleGridIntegration(config);
	if (adj_plasma) integration_container[ADJPLASMA_SOL] = new CMultiGridIntegration(config);
	if (adj_levelset) integration_container[ADJLEVELSET_SOL] = new CSingleGridIntegration(config);

	/*--- Allocate solution for linear problem (at the moment we use the same scheme as the adjoint problem) ---*/
	if (lin_pot) { cout <<"Equation not implemented." << endl; cin.get(); }
	if (lin_euler) integration_container[LINFLOW_SOL] = new CMultiGridIntegration(config);
	if (lin_ns) { cout <<"Equation not implemented." << endl; cin.get(); }

}


void Solver_Definition(CNumerics ****solver_container, CSolution ***solution_container, CGeometry **geometry, CConfig *config) {

	unsigned short iMGlevel, iSol, nDim, nVar_Template = 0, nVar_Flow = 0, nVar_Adj_Flow = 0, nVar_Plasma = 0, nVar_LevelSet = 0, nVar_Turb = 0, nVar_Adj_Turb = 0, 
			nVar_Elec = 0, nVar_Wave = 0, nVar_Lin_Flow = 0, nVar_Adj_LevelSet = 0, nVar_Adj_Plasma = 0, nSpecies = 0, nFluids = 0, nDiatomics = 0, nMonatomics = 0;
	bool potential, euler, navierstokes, combustion, plasma, plasma_monatomic, plasma_diatomic, levelset, adj_pot, adj_plasma, adj_levelset, lin_pot, adj_euler, lin_euler, adj_ns, lin_ns, turbulent, 
	adj_turb, lin_turb, electric, wave, spalart_allmaras, menter_sst, template_solver;

	bool incompressible = config->GetIncompressible();

	/*--- Initialize some useful booleans ---*/
	potential = false;	euler = false;		navierstokes = false;	combustion = false; turbulent = false;	electric = false;	plasma_monatomic = false;	
	plasma_diatomic = false; levelset = false; plasma = false;
	adj_pot = false;	adj_euler = false;	adj_ns = false;			adj_turb = false;	wave = false; adj_levelset = false;	spalart_allmaras = false;
	lin_pot = false;	lin_euler = false;	lin_ns = false;			lin_turb = false;	menter_sst = false; adj_plasma = false;
	template_solver = false;

	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
	case TEMPLATE_SOLVER: template_solver = true; break;
	case EULER : euler = true; break;
	case NAVIER_STOKES: navierstokes = true; break;
	case RANS : navierstokes = true; turbulent = true; break;
	case FREE_SURF_EULER: euler = true; levelset = true; break;
	case FREE_SURF_NAVIER_STOKES: navierstokes = true; levelset = true; break;
	case FREE_SURF_RANS: navierstokes = true; turbulent = true; levelset = true; break;
	case NS_PLASMA : plasma = true; break;
	case ELECTRIC_POTENTIAL: electric = true; break;
	case WAVE: wave = true; break;
	case ADJ_EULER : euler = true; adj_euler = true; break;
	case ADJ_NAVIER_STOKES : navierstokes = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
	case ADJ_RANS : navierstokes = true; turbulent = true; adj_ns = true; adj_turb = true; break;
	case ADJ_FREE_SURF_EULER: euler = true; adj_euler = true; levelset = true; adj_levelset = true; break;
	case ADJ_FREE_SURF_NAVIER_STOKES: navierstokes = true; adj_ns = true; levelset = true; adj_levelset = true; break;
	case ADJ_FREE_SURF_RANS: navierstokes = true; adj_ns = true; turbulent = true; adj_turb = true; levelset = true; adj_levelset = true; break;
	case ADJ_NS_PLASMA : plasma = true; adj_plasma = true; break;
	case LIN_EULER: euler = true; lin_euler = true; break;
	}

	/*--- Assign turbulence model booleans --- */
	if (turbulent)
		switch (config->GetKind_Turb_Model()){
		case SA: case SA_COMP: spalart_allmaras = true; break;
		case SST: menter_sst = true; break;
		default: cout << "Specified turbulence model unavailable or none selected" << endl; cin.get(); break;
		}

	if (plasma) {
		switch (config->GetKind_GasModel()){
		case AIR7: plasma_diatomic = true; break;
		case ARGON: plasma_monatomic = true; break;
		case AIR21: plasma_diatomic = true; break;
		default: cout << "Specified plasma model unavailable or none selected" << endl; cin.get(); break;
		}
		if (config->GetElectricSolver()) electric  = true;
	}

	/*--- Number of variables for the template ---*/
	if (template_solver) nVar_Flow = solution_container[MESH_0][FLOW_SOL]->GetnVar();

	/*--- Number of variables for direct problem ---*/
	if (potential)		nVar_Flow = solution_container[MESH_0][FLOW_SOL]->GetnVar();
	if (euler)				nVar_Flow = solution_container[MESH_0][FLOW_SOL]->GetnVar();
	if (navierstokes)	nVar_Flow = solution_container[MESH_0][FLOW_SOL]->GetnVar();
	if (turbulent)		nVar_Turb     = solution_container[MESH_0][TURB_SOL]->GetnVar();
	if (electric)			nVar_Elec = solution_container[MESH_0][ELEC_SOL]->GetnVar();
	if (plasma)	{ 
		nVar_Plasma = solution_container[MESH_0][PLASMA_SOL]->GetnVar();
		nSpecies    = solution_container[MESH_0][PLASMA_SOL]->GetnSpecies();
		nFluids     = solution_container[MESH_0][PLASMA_SOL]->GetnFluids();
		nDiatomics  = solution_container[MESH_0][PLASMA_SOL]->GetnDiatomics();
		nMonatomics = solution_container[MESH_0][PLASMA_SOL]->GetnMonatomics();
	}
	if (levelset)		nVar_LevelSet = solution_container[MESH_0][LEVELSET_SOL]->GetnVar();
	if (wave)				nVar_Wave = solution_container[MESH_0][WAVE_SOL]->GetnVar();

	/*--- Number of variables for adjoint problem ---*/
	if (adj_pot)		  nVar_Adj_Flow = solution_container[MESH_0][ADJFLOW_SOL]->GetnVar();
	if (adj_euler)  	nVar_Adj_Flow = solution_container[MESH_0][ADJFLOW_SOL]->GetnVar();
	if (adj_ns)			  nVar_Adj_Flow = solution_container[MESH_0][ADJFLOW_SOL]->GetnVar();
	if (adj_turb)		  nVar_Adj_Turb = solution_container[MESH_0][ADJTURB_SOL]->GetnVar();
	if (adj_levelset)	nVar_Adj_LevelSet = solution_container[MESH_0][ADJLEVELSET_SOL]->GetnVar();
	if (adj_plasma)		nVar_Adj_Plasma = solution_container[MESH_0][ADJPLASMA_SOL]->GetnVar();

	/*--- Number of variables for the linear problem ---*/
	if (lin_euler)	nVar_Lin_Flow = solution_container[MESH_0][LINFLOW_SOL]->GetnVar();

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
				case LAX : solver_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim,nVar_Flow, config); break;
				case JST : solver_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim,nVar_Flow, config); break;
				default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
				}
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim,nVar_Flow, config);

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
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
						solver_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config);
						solver_container[iMGlevel][FLOW_SOL][BOUND_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config);
					}
					break;
				default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
				}					

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
				else if (config->GetAxisymmetric() == YES)
					solver_container[iMGlevel][FLOW_SOL][SOUR_TERM] = new CSourceAxisymmetric_Flow(nDim,nVar_Flow, config);
				else if (config->GetGravityForce() == YES)
					solver_container[iMGlevel][FLOW_SOL][SOUR_TERM] = new CSourcePieceWise_Gravity(nDim, nVar_Flow, config);
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
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
				if (config->GetKind_GasModel() == ARGON) {
					solver_container[iMGlevel][PLASMA_SOL][CONV_TERM] = new CUpwRoe_Plasma(nDim, nVar_Plasma, nSpecies, nFluids, config);
					solver_container[iMGlevel][PLASMA_SOL][BOUND_TERM] = new CUpwRoe_Plasma(nDim, nVar_Plasma, nSpecies, nFluids,config);
				}
				if (config->GetKind_GasModel() == AIR7) {
					solver_container[iMGlevel][PLASMA_SOL][CONV_TERM] = new CUpwRoe_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
					solver_container[iMGlevel][PLASMA_SOL][BOUND_TERM] = new CUpwRoe_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				}
				if (config->GetKind_GasModel() == AIR21) {
					solver_container[iMGlevel][PLASMA_SOL][CONV_TERM] = new CUpwRoe_PlasmaDiatomic_AM(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
					solver_container[iMGlevel][PLASMA_SOL][BOUND_TERM] = new CUpwRoe_PlasmaDiatomic_AM(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);

				}
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
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][PLASMA_SOL][VISC_TERM] = new CAvgGrad_Plasma(nDim, nVar_Plasma, nSpecies, nFluids, config);
			}
			if (plasma_diatomic) {
				if (config->GetKind_GasModel() == AIR21) {

					//			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					//			solver_container[iMGlevel][PLASMA_SOL][VISC_TERM] = new CAvgGrad_Plasma(nDim, nVar_Plasma, nDiatomics, nMonatomics, config);
				}
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
				if (config->GetKind_GasModel() == ARGON)
					solver_container[iMGlevel][PLASMA_SOL][SOUR_TERM] = new CSourcePieceWise_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				if (config->GetKind_GasModel() == AIR7)
					solver_container[iMGlevel][PLASMA_SOL][SOUR_TERM] = new CSourcePieceWise_PlasmaDiatomic(nDim, nVar_Plasma, config);
				if (config->GetKind_GasModel() == AIR21)
					solver_container[iMGlevel][PLASMA_SOL][SOUR_TERM] = new CSourcePieceWise_Plasma_Air(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);

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
				case SCALAR_UPWIND_1ST : case SCALAR_UPWIND_2ND :
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

		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_LevelSet()) {
		case NONE :
			break;
		case PIECEWISE_CONSTANT :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
				solver_container[iMGlevel][LEVELSET_SOL][SOUR_TERM] = new CSourcePieceWise_LevelSet(nDim, nVar_LevelSet, config);
			break;
		default :
			cout << "Source term not implemented." << endl; cin.get();
			break;
		}
	}

	/*--- Solver definition for the flow adjoint problem ---*/
	if ((adj_pot)||(adj_euler)||(adj_ns)) {

		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_AdjFlow()) {
		case NO_CONVECTIVE :
			cout << "No convective scheme." << endl; cin.get();
			break;
		case SPACE_CENTRED :
			if (incompressible) {
				/*--- Incompressible flow, use artificial compressibility method ---*/
				switch (config->GetKind_Centred_AdjFlow()) {
				case NO_CENTRED : cout << "No centered scheme." << endl; break;
				case LAX : solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLaxArtComp_AdjFlow(nDim, nVar_Adj_Flow, config); break;
				case JST : solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentJSTArtComp_AdjFlow(nDim, nVar_Adj_Flow, config); break;
				default : cout << "Centred scheme not implemented." << endl; cin.get(); break;
				}
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLaxArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);

				/*--- Definition of the boundary condition method ---*/
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][BOUND_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);

			}
			else {
				/*--- Compressible flow ---*/
				switch (config->GetKind_Centred_AdjFlow()) {
				case NO_CENTRED : cout << "No centered scheme." << endl; break;
				case LAX : solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config); break;
				case JST : solver_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentJST_AdjFlow(nDim, nVar_Adj_Flow, config); break;
				default : cout << "Centred scheme not implemented." << endl; cin.get(); break;
				}
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config);

				/*--- Definition of the boundary condition method ---*/
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
			}
			break;
		case SPACE_UPWIND :
			if (incompressible) {
				/*--- Incompressible flow, use artificial compressibility method ---*/
				switch (config->GetKind_Upwind_AdjFlow()) {
				case NO_UPWIND : cout << "No upwind scheme." << endl; break;
				case ROE_1ST : case ROE_2ND :
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
						solver_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
						solver_container[iMGlevel][ADJFLOW_SOL][BOUND_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
					}
					break;
				default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
				}
			}
			else {
				/*--- Compressible flow ---*/
				switch (config->GetKind_Upwind_AdjFlow()) {
				case NO_UPWIND : cout << "No upwind scheme." << endl; break;
				case ROE_1ST : case ROE_2ND :
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
						solver_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
						solver_container[iMGlevel][ADJFLOW_SOL][BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
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
		switch (config->GetKind_ViscNumScheme_AdjFlow()) {
		case NONE :
			break;
		case AVG_GRAD :
			if (incompressible) {
				/*--- Incompressible flow, use artificial compressibility method ---*/
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
			}
			else {
				/*--- Compressible flow ---*/
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
			}
			break;
		case AVG_GRAD_CORRECTED :
			if (incompressible) {
				/*--- Incompressible flow, use artificial compressibility method ---*/
				solver_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
			}
			else {
				/*--- Compressible flow ---*/
				solver_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrected_AdjFlow(nDim, nVar_Adj_Flow, config);
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					solver_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
			}
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
				if (config->GetRotating_Frame() == YES)
					solver_container[iMGlevel][ADJFLOW_SOL][SOUR_TERM] = new CSourceRotationalFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
				else
					solver_container[iMGlevel][ADJFLOW_SOL][SOUR_TERM] = new CSourcePieceWise_AdjFlow(nDim, nVar_Adj_Flow, config);

				solver_container[iMGlevel][ADJFLOW_SOL][CONS_SOUR_TERM] = new CSourceConservative_AdjFlow(nDim, nVar_Adj_Flow, config);
			}
			break;
		default :
			cout << "Source term not implemented." << endl; cin.get();
			break;
		}
	}

	/*--- Solver definition for the multi species plasma model problem ---*/
	if (adj_plasma) {

		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Plasma()) {
		case NONE :
			break;
		case SPACE_UPWIND :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
				if (config->GetKind_GasModel() == AIR7) {
					solver_container[iMGlevel][PLASMA_SOL][CONV_TERM] = new CUpwRoe_AdjPlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
					solver_container[iMGlevel][PLASMA_SOL][BOUND_TERM] = new CUpwRoe_AdjPlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				}
			}
			break;
		default :
			cout << "Convective scheme not implemented." << endl; cin.get();
			break;
		}


		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		/*		switch (config->GetKind_SourNumScheme_Plasma()) {
			case NONE :
				break;
			case PIECEWISE_CONSTANT :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					if (config->GetKind_GasModel() == AIR7)
						solver_container[iMGlevel][PLASMA_SOL][SOUR_TERM] = new CSourcePieceWise_ADjPlasmaDiatomic(nDim, nVar_Plasma, config);

					solver_container[iMGlevel][PLASMA_SOL][CONS_SOUR_TERM] = new CSourceNothing(nDim, nVar_Plasma, config);
				}
				break;
			default :
				cout << "Source term not implemented." << endl; cin.get();
				break;
		} */
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
		switch (config->GetKind_SourNumScheme_AdjTurb()) {
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

	/*--- Solver definition for the level set model problem ---*/
	if (adj_levelset) {

		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_AdjLevelSet()) {
		case NO_CONVECTIVE : cout << "No convective scheme." << endl; cin.get(); break;
		case SPACE_CENTRED :
			switch (config->GetKind_Centred_AdjLevelSet()) {
			case NO_UPWIND : cout << "No centered scheme." << endl; cin.get(); break;
			default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
			}
			break;
			case SPACE_UPWIND :
				switch (config->GetKind_Upwind_AdjLevelSet()) {
				case NO_UPWIND : cout << "No upwind scheme." << endl; cin.get(); break;
				case SCALAR_UPWIND_1ST : case SCALAR_UPWIND_2ND :
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
						solver_container[iMGlevel][ADJLEVELSET_SOL][CONV_TERM] = new CUpwLin_AdjLevelSet(nDim, nVar_Adj_LevelSet, config);
						solver_container[iMGlevel][ADJLEVELSET_SOL][BOUND_TERM] = new CUpwLin_AdjLevelSet(nDim, nVar_Adj_LevelSet, config);
					}
					break;
				default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
				}
				break;
				default : cout << "Convective scheme not implemented." << endl; cin.get(); break;
		}

		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_AdjLevelSet()) {
		case NONE :
			break;
		case PIECEWISE_CONSTANT :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
				solver_container[iMGlevel][ADJLEVELSET_SOL][SOUR_TERM] = new CSourcePieceWise_AdjLevelSet(nDim, nVar_Turb, config);
			break;
		default :
			cout << "Source term not implemented." << endl; cin.get();
			break;
		}
	}

	/*--- Solver definition for the adjoint electric potential problem ---*/
	if (wave) {

		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Wave()) {
		case NONE :
			break;
		case AVG_GRAD :
			cout << "Viscous scheme not implemented." << endl; cin.get();
			break;
		case AVG_GRAD_CORRECTED :
			cout << "Viscous scheme not implemented." << endl; cin.get();
			break;
		case GALERKIN :
			solver_container[MESH_0][WAVE_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Wave, config);
			break;
		default :
			cout << "Viscous scheme not implemented." << endl; cin.get();
			break;
		}

		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Wave()) {
		case NONE : break;
		case PIECEWISE_CONSTANT : break;
		default : break;
		}
	}

}

void Geometrical_Deallocation(CGeometry **geometry, CConfig *config) {

	/*--- Deallocation of geometry pointer ---*/
	/*if (config->GetKind_Solver() != NO_SOLVER)
	 for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
	 delete geometry[iMGlevel];
	 delete geometry[iDomain][MESH_0];
	 delete [] geometry;*/

}

void Solver_Deallocation(CNumerics ****solver_container, CSolution ***solution_container, CIntegration **integration_container,
		COutput *output, CGeometry **geometry, CConfig *config){

	unsigned short iMGlevel;
	bool potential, euler, navierstokes, plasma, adj_pot, lin_pot, adj_euler, lin_euler, adj_ns, lin_ns, turbulent,
	adj_turb, lin_turb, electric, wave, spalart_allmaras, menter_sst;

	/*--- Initialize some useful booleans ---*/
	potential = false;	euler = false;		navierstokes = false;	turbulent = false;	electric = false;	plasma = false;
	adj_pot = false;	adj_euler = false;	adj_ns = false;			adj_turb = false;	wave = false; 	spalart_allmaras = false;
	lin_pot = false;	lin_euler = false;	lin_ns = false;			lin_turb = false;	menter_sst = false;

	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
	case EULER : euler = true; break;
	case NAVIER_STOKES: navierstokes = true; break;
	case RANS : navierstokes = true; turbulent = true; break;
	case NS_PLASMA : plasma = true; break;
	case ELECTRIC_POTENTIAL: electric = true; break;
	case WAVE: wave = true; break;
	case ADJ_EULER : euler = true; adj_euler = true; break;
	case ADJ_NAVIER_STOKES : navierstokes = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
	case ADJ_RANS : navierstokes = true; turbulent = true; adj_ns = true; adj_turb = true; break;
	case LIN_EULER: euler = true; lin_euler = true; break;
	}

	/*--- Assign turbulence model booleans --- */
	if (turbulent)
		switch (config->GetKind_Turb_Model()){
		case SA: case SA_COMP: spalart_allmaras = true; break;
		case SST: menter_sst = true; break;
		}

	/*--- Deallocation of integration_container---*/
	if (potential)		delete integration_container[FLOW_SOL];
	if (euler)          delete integration_container[FLOW_SOL];
	if (navierstokes)   delete integration_container[FLOW_SOL];
	if (turbulent)     delete integration_container[TURB_SOL]; 
	if (electric) 	    delete integration_container[ELEC_SOL];
	if (plasma) 	    delete integration_container[PLASMA_SOL];
	if (wave) 	    delete integration_container[WAVE_SOL];
	if (adj_euler)      delete integration_container[ADJFLOW_SOL];
	if (adj_ns) 	    delete integration_container[ADJFLOW_SOL];
	if (adj_turb) 	    delete integration_container[ADJTURB_SOL];
	if (lin_euler) 	    delete integration_container[LINFLOW_SOL];

	delete [] integration_container;

	/*--- Deallocation of solution_container---*/
	for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
		//if (potential) 		delete solution_container[iMGlevel][FLOW_SOL];
		if (euler) 		    delete solution_container[iMGlevel][FLOW_SOL];
		//if (navierstokes)   delete solution_container[iMGlevel][FLOW_SOL];
		//if (turbulent) {
		//	if (spalart_allmaras) delete solution_container[iMGlevel][TURB_SOL];
		//	else if (menter_sst)  delete solution_container[iMGlevel][TURB_SOL];
		//}
		//if (wave) 		delete solution_container[iMGlevel][WAVE_SOL];
		//if (electric) 		delete solution_container[iMGlevel][ELEC_SOL];
		//if (plasma) 		delete solution_container[iMGlevel][PLASMA_SOL];
		//if (adj_euler) 		delete solution_container[iMGlevel][ADJFLOW_SOL];
		//if (adj_ns) 		delete solution_container[iMGlevel][ADJFLOW_SOL];
		//if (adj_turb)		delete solution_container[iMGlevel][ADJTURB_SOL];
		//if (lin_euler) 		delete solution_container[iMGlevel][LINFLOW_SOL];
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
	 if (wave) {

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
	 delete solver_container[MESH_0][WAVE_SOL][VISC_TERM];
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
	 default :
	 cout << "Source term not implemented." << endl; cin.get();
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

