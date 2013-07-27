/*!
 * \file definition_structure.cpp
 * \brief Main subroutines used by SU2_CFD.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
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


unsigned short GetnZone(string val_mesh_filename, unsigned short val_format, CConfig *config) {
	string text_line, Marker_Tag;
	ifstream mesh_file;
	short nZone = 1;
	bool isFound = false;
	char cstr[200];
	string::size_type position;
	int rank = MASTER_NODE;

#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
	if (MPI::COMM_WORLD.Get_size() != 1) {
		val_mesh_filename.erase (val_mesh_filename.end()-4, val_mesh_filename.end());
		val_mesh_filename = val_mesh_filename + "_1.su2";
	}
#endif

	/*--- Search the mesh file for the 'NZONE' keyword. ---*/
	switch (val_format) {
	case SU2:

		/*--- Open grid file ---*/
		strcpy (cstr, val_mesh_filename.c_str());
		mesh_file.open(cstr, ios::in);
		if (mesh_file.fail()) {
			cout << "cstr=" << cstr << endl;
			cout << "There is no geometry file (GetnZone))!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
#ifdef NO_MPI
			exit(1);
#else
			MPI::COMM_WORLD.Abort(1);
			MPI::Finalize();
#endif
		}

		/*--- Open the SU2 mesh file ---*/
		while (getline (mesh_file,text_line)) {

			/*--- Search for the "NZONE" keyword to see if there are multiple Zones ---*/
			position = text_line.find ("NZONE=",0);
			if (position != string::npos) {
				text_line.erase (0,6); nZone = atoi(text_line.c_str()); isFound = true;
				if (rank == MASTER_NODE) {
					//					if (nZone == 1) cout << "SU2 mesh file format with a single zone." << endl;
					//					else if (nZone >  1) cout << "SU2 mesh file format with " << nZone << " zones." << endl;
					//					else
					if (nZone <= 0) {
						cout << "Error: Number of mesh zones is less than 1 !!!" << endl;
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
		/*--- If the "NZONE" keyword was not found, assume this is an ordinary
       simulation on a single Zone ---*/
		if (!isFound) {
			nZone = 1;
			//			if (rank == MASTER_NODE) cout << "SU2 mesh file format with a single zone." << endl;
		}
		break;

	case CGNS:

		nZone = 1;
		//		if (rank == MASTER_NODE) cout << "CGNS mesh file format with a single zone." << endl;
		break;

	case NETCDF_ASCII:

		nZone = 1;
		//		if (rank == MASTER_NODE) cout << "NETCDF mesh file format with a single zone." << endl;
		break;

	}

	/*--- For time spectral integration, nZones = nTimeInstances. ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		nZone = config->GetnTimeInstances();
	}

	return (unsigned short) nZone;
}

unsigned short GetnDim(string val_mesh_filename, unsigned short val_format) {

	string text_line, Marker_Tag;
	ifstream mesh_file;
	short nDim = 3;
	bool isFound = false;
	char cstr[200];
	string::size_type position;

#ifndef NO_MPI
	if (MPI::COMM_WORLD.Get_size() != 1) {
		val_mesh_filename.erase (val_mesh_filename.end()-4, val_mesh_filename.end());
		val_mesh_filename = val_mesh_filename + "_1.su2";
	}
#endif

	switch (val_format) {
	case SU2:

		/*--- Open grid file ---*/
		strcpy (cstr, val_mesh_filename.c_str());
		mesh_file.open(cstr, ios::in);

		/*--- Read SU2 mesh file ---*/
		while (getline (mesh_file,text_line)) {
			/*--- Search for the "NDIM" keyword to see if there are multiple Zones ---*/
			position = text_line.find ("NDIME=",0);
			if (position != string::npos) {
				text_line.erase (0,6); nDim = atoi(text_line.c_str()); isFound = true;
			}
		}
		break;

	case CGNS:
		nDim = 3;
		break;

	case NETCDF_ASCII:
		nDim = 3;
		break;
	}
	return (unsigned short) nDim;
}



void Geometrical_Preprocessing(CGeometry ***geometry, CConfig **config, unsigned short val_nZone) {


	unsigned short iMGlevel, iZone;
	unsigned long iPoint;
	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	for (iZone = 0; iZone < val_nZone; iZone++) {

		/*--- Compute elements surrounding points, points surrounding points,
     and elements surrounding elements ---*/
		if (rank == MASTER_NODE) cout << "Setting local point and element connectivity." <<endl;
		geometry[iZone][MESH_0]->SetEsuP();
		geometry[iZone][MESH_0]->SetPsuP();
		geometry[iZone][MESH_0]->SetEsuE();
        
		/*--- Check the orientation before computing geometrical quantities ---*/
		if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." <<endl;
		geometry[iZone][MESH_0]->SetBoundVolume();
		geometry[iZone][MESH_0]->Check_Orientation(config[iZone]);

		/*--- Create the edge structure ---*/
		if (rank == MASTER_NODE) cout << "Identifying edges and vertices." <<endl;
		geometry[iZone][MESH_0]->SetEdges();
		geometry[iZone][MESH_0]->SetVertex(config[iZone]);

		/*--- Compute center of gravity ---*/
		if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
		geometry[iZone][MESH_0]->SetCG();

		/*--- Create the control volume structures ---*/
		if (rank == MASTER_NODE) cout << "Setting the control volume structure." << endl;
		geometry[iZone][MESH_0]->SetControlVolume(config[iZone], ALLOCATE);
		geometry[iZone][MESH_0]->SetBoundControlVolume(config[iZone], ALLOCATE);

		/*--- Identify closest normal neighbor ---*/
		if (rank == MASTER_NODE) cout << "Searching for closest normal neighbor on the surface." << endl;
		geometry[iZone][MESH_0]->FindNormal_Neighbor(config[iZone]);

		/*--- For a rotating frame, set the velocity due to rotation at each mesh point ---*/
		if (config[iZone]->GetRotating_Frame())
			geometry[iZone][MESH_0]->SetRotationalVelocity(config[iZone]);

		if ((config[iZone]->GetMGLevels() != 0) && (rank == MASTER_NODE))
			cout << "Setting the multigrid structure." <<endl;

	}

#ifndef NO_MPI
	/*--- Synchronization point after the multigrid stuff ---*/
	MPI::COMM_WORLD.Barrier();
#endif

	/*--- Loop over all the new grid ---*/
	for (iMGlevel = 1; iMGlevel <= config[ZONE_0]->GetMGLevels(); iMGlevel++) {

		/*--- Loop over all zones at each grid level. ---*/
		for (iZone = 0; iZone < val_nZone; iZone++) {

			/*--- Create main aglomeration structure (ingluding MPI stuff) ---*/
			geometry[iZone][iMGlevel] = new CMultiGridGeometry(geometry, config, iMGlevel, iZone);

			/*--- Compute points surrounding points. ---*/
			geometry[iZone][iMGlevel]->SetPsuP(geometry[iZone][iMGlevel-1]);

			/*--- Create the edge structure ---*/
			geometry[iZone][iMGlevel]->SetEdges();
			geometry[iZone][iMGlevel]->SetVertex(geometry[iZone][iMGlevel-1], config[iZone]);

			/*--- Create the control volume structures ---*/
			geometry[iZone][iMGlevel]->SetControlVolume(config[iZone],geometry[iZone][iMGlevel-1], ALLOCATE);
			geometry[iZone][iMGlevel]->SetBoundControlVolume(config[iZone],geometry[iZone][iMGlevel-1], ALLOCATE);
			geometry[iZone][iMGlevel]->SetCoord(geometry[iZone][iMGlevel-1]);

			/*--- Find closest neighbor to a surface point ---*/
			geometry[iZone][iMGlevel]->FindNormal_Neighbor(config[iZone]);

			/*--- For a rotating frame, set the velocity due to rotation at each mesh point ---*/
			if (config[iZone]->GetRotating_Frame())
				geometry[iZone][iMGlevel]->SetRotationalVelocity(config[iZone]);

		}
	}

	/*--- For unsteady simulations, initialize the grid volumes
   and coordinates for previous solutions. Loop over all zones/grids. ---*/
	for (iZone = 0; iZone < val_nZone; iZone++) {
		if (config[iZone]->GetGrid_Movement()) {
			for (iMGlevel = 0; iMGlevel <= config[iZone]->GetMGLevels(); iMGlevel++) {
				for (iPoint = 0; iPoint < geometry[iZone][iMGlevel]->GetnPoint(); iPoint++) {
					geometry[iZone][iMGlevel]->node[iPoint]->SetVolume_n();
					geometry[iZone][iMGlevel]->node[iPoint]->SetVolume_nM1();

					geometry[iZone][iMGlevel]->node[iPoint]->SetCoord_n();
					geometry[iZone][iMGlevel]->node[iPoint]->SetCoord_n1();
				}
			}
		}
	}

}

void Solver_Preprocessing(CSolver ***solver_container, CGeometry **geometry, CConfig *config, unsigned short iZone) {

	unsigned short iMGlevel;
	bool euler, ns, plasma_euler, plasma_ns,
	plasma_monatomic, plasma_diatomic, levelset, adj_euler,
	lin_euler, adj_ns, lin_ns, turbulent, adj_turb, lin_turb, electric, wave, fea, adj_levelset, 
	adj_plasma_euler, adj_plasma_ns, spalart_allmaras, menter_sst, template_solver, transition;


	/*--- Initialize some useful booleans ---*/
	euler = false;		ns = false; turbulent = false;	electric = false;	plasma_monatomic = false;
	plasma_diatomic = false; levelset = false; plasma_euler = false; plasma_ns = false; transition = false;
	adj_euler = false;	adj_ns = false;			adj_turb = false;	wave = false; fea = false; 	adj_levelset = false;	 spalart_allmaras = false;
	lin_euler = false;	lin_ns = false;			lin_turb = false;	menter_sst = false; adj_plasma_euler = false;	adj_plasma_ns = false;
	template_solver = false;

	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
	case TEMPLATE_SOLVER: template_solver = true; break;
	case EULER : euler = true; break;
	case NAVIER_STOKES: ns = true; break;
	case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
	case FREE_SURFACE_EULER: euler = true; levelset = true; break;
	case FREE_SURFACE_NAVIER_STOKES: ns = true; levelset = true; break;
	case FREE_SURFACE_RANS: ns = true; turbulent = true; levelset = true; break;
	case FLUID_STRUCTURE_EULER: euler = true; fea = true; break;
	case FLUID_STRUCTURE_NAVIER_STOKES: ns = true; fea = true; break;
	case FLUID_STRUCTURE_RANS: ns = true; turbulent = true; fea = true; break;
	case AEROACOUSTIC_NAVIER_STOKES: ns = true; wave = true; break;
	case AEROACOUSTIC_RANS: ns = true; turbulent = true; wave = true; break;
	case ELECTRIC_POTENTIAL: electric = true; break;
	case WAVE_EQUATION: wave = true; break;
	case LINEAR_ELASTICITY: fea = true; break;
	case ADJ_EULER : euler = true; adj_euler = true; break;
	case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
	case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = ( ((config->GetKind_Adjoint() == CONTINUOUS) && (!config->GetFrozen_Visc())) || (config->GetKind_Adjoint() == HYBRID) ); break;
	case ADJ_FREE_SURFACE_EULER: euler = true; adj_euler = true; levelset = true; adj_levelset = true; break;
	case ADJ_FREE_SURFACE_NAVIER_STOKES: ns = true; adj_ns = true; levelset = true; adj_levelset = true; break;
	case ADJ_FREE_SURFACE_RANS: ns = true; adj_ns = true; turbulent = true; adj_turb = (((config->GetKind_Adjoint() != HYBRID) && (!config->GetFrozen_Visc())) || (config->GetKind_Adjoint() == HYBRID)); levelset = true; adj_levelset = true; break;
	case ADJ_PLASMA_EULER : plasma_euler = true; adj_plasma_euler = true; break;
	case ADJ_PLASMA_NAVIER_STOKES : plasma_ns = true; adj_plasma_ns = true; break;
	case LIN_EULER: euler = true; lin_euler = true; break;

	/*--- Specify by zone for the aeroacoustic problem ---*/
	case AEROACOUSTIC_EULER:
		if (iZone == ZONE_0) {
			euler = true;
		} else if (iZone == ZONE_1) {
			wave = true;
		}
		break;
	case ADJ_AEROACOUSTIC_EULER:
		if (iZone == ZONE_0) {
			euler = true; adj_euler = true;
		} else if (iZone == ZONE_1) {
			wave = true;
		}
		break;
	case PLASMA_EULER:
		if (iZone == ZONE_0) {
			plasma_euler = true;
		} else if (iZone == ZONE_1) {
			electric = true;
		}
		break;
	case PLASMA_NAVIER_STOKES:
		if (iZone == ZONE_0) {
			plasma_ns = true;
		} else if (iZone == ZONE_1) {
			electric = true;
		}
		break;
	}
	/*--- Assign turbulence model booleans --- */
	if (turbulent)
		switch (config->GetKind_Turb_Model()){
		case SA: spalart_allmaras = true; break;
		case SST: menter_sst = true; break;
		default: cout << "Specified turbulence model unavailable or none selected" << endl; cin.get(); break;
		}

	if (plasma_euler || plasma_ns) {
		switch (config->GetKind_GasModel()){
		case AIR7: plasma_diatomic = true; break;
		case O2: plasma_diatomic = true; break;
		case N2: plasma_diatomic = true; break;
		case AIR5: plasma_diatomic = true; break;
		case ARGON: plasma_monatomic = true; break;
		case ARGON_SID: plasma_diatomic = true; break;
		default: cout << "Specified plasma model unavailable or none selected" << endl; cin.get(); break;
		}
		//if (config->GetElectricSolver()) electric  = true;
	}

	/*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/
	for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {

		/*--- Allocate solution for a template problem ---*/
		if (template_solver) {
			solver_container[iMGlevel][TEMPLATE_SOL] = new CTemplateSolver(geometry[iMGlevel], config);
		}

		/*--- Allocate solution for direct problem, and run the preprocessing and postprocessing ---*/
		if (euler) {
			solver_container[iMGlevel][FLOW_SOL] = new CEulerSolver(geometry[iMGlevel], config, iMGlevel);
		}
		if (ns) {
			solver_container[iMGlevel][FLOW_SOL] = new CNSSolver(geometry[iMGlevel], config, iMGlevel);
		}
		if (turbulent) {
			if (spalart_allmaras) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbSASolver(geometry[iMGlevel], config, iMGlevel);
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
			else if (menter_sst) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbSSTSolver(geometry[iMGlevel], config, iMGlevel);
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
			if (transition) {
        solver_container[iMGlevel][TRANS_SOL] = new CTransLMSolver(geometry[iMGlevel], config, iMGlevel);
      }
		}
		if (electric) {
			solver_container[iMGlevel][ELEC_SOL] = new CElectricSolver(geometry[iMGlevel], config);
		}
		if (plasma_euler || plasma_ns) {
			solver_container[iMGlevel][PLASMA_SOL] = new CPlasmaSolver(geometry[iMGlevel], config);
		}
		if (levelset) {
			solver_container[iMGlevel][LEVELSET_SOL] = new CLevelSetSolver(geometry[iMGlevel], config, iMGlevel);
		}		
		if (wave) {
			solver_container[iMGlevel][WAVE_SOL] = new CWaveSolver(geometry[iMGlevel], config);
		}
		if (fea) {
			solver_container[iMGlevel][FEA_SOL] = new CFEASolver(geometry[iMGlevel], config);
		}

		/*--- Allocate solution for adjoint problem ---*/
		if (adj_euler) {
			solver_container[iMGlevel][ADJFLOW_SOL] = new CAdjEulerSolver(geometry[iMGlevel], config, iMGlevel);
		}
		if (adj_ns) {
			solver_container[iMGlevel][ADJFLOW_SOL] = new CAdjNSSolver(geometry[iMGlevel], config, iMGlevel);
		}
		if (adj_turb) {
			solver_container[iMGlevel][ADJTURB_SOL] = new CAdjTurbSolver(geometry[iMGlevel], config);
		}
		if (adj_levelset) {
			solver_container[iMGlevel][ADJLEVELSET_SOL] = new CAdjLevelSetSolver(geometry[iMGlevel], config, iMGlevel);
		}
		if (adj_plasma_euler || adj_plasma_ns) {
			solver_container[iMGlevel][ADJPLASMA_SOL] = new CAdjPlasmaSolver(geometry[iMGlevel], config);
		}

		/*--- Allocate solution for linear problem (at the moment we use the same scheme as the adjoint problem) ---*/
		if (lin_euler) {
			solver_container[iMGlevel][LINFLOW_SOL] = new CLinEulerSolver(geometry[iMGlevel], config);
		}
		if (lin_ns) {
			cout <<"Equation not implemented." << endl; cin.get(); break;
		}

	}

}


void Integration_Preprocessing(CIntegration **integration_container, CGeometry **geometry, CConfig *config, unsigned short iZone) {

	bool euler, ns, plasma_euler, plasma_ns, plasma_monatomic, plasma_diatomic, levelset, adj_plasma_euler, adj_plasma_ns,
	adj_euler, lin_euler, adj_ns, lin_ns, turbulent, adj_turb, lin_turb, electric, wave, fea, spalart_allmaras, menter_sst, template_solver, adj_levelset, transition;

	/*--- Initialize some useful booleans ---*/
	euler = false;		ns = false; turbulent = false;	electric = false;	plasma_monatomic = false;
	plasma_diatomic = false; levelset = false; plasma_euler = false; plasma_ns = false;
	adj_euler = false;	adj_ns = false;			adj_turb = false;	wave = false; fea = false; adj_levelset = false;	spalart_allmaras = false;
	lin_euler = false;	lin_ns = false;			lin_turb = false;	menter_sst = false; adj_plasma_euler = false; adj_plasma_ns = false; transition = false;
	template_solver = false;

	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
        case TEMPLATE_SOLVER: template_solver = true; break;
        case EULER : euler = true; break;
        case NAVIER_STOKES: ns = true; break;
        case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
        case FREE_SURFACE_EULER: euler = true; levelset = true; break;
        case FREE_SURFACE_NAVIER_STOKES: ns = true; levelset = true; break;
        case FREE_SURFACE_RANS: ns = true; turbulent = true; levelset = true; break;
        case FLUID_STRUCTURE_EULER: euler = true; fea = true; break;
        case FLUID_STRUCTURE_NAVIER_STOKES: ns = true; fea = true; break;
        case FLUID_STRUCTURE_RANS: ns = true; turbulent = true; fea = true; break;
        case AEROACOUSTIC_NAVIER_STOKES: ns = true; wave = true; break;
        case AEROACOUSTIC_RANS: ns = true; turbulent = true; wave = true; break;
        case ELECTRIC_POTENTIAL: electric = true; break;
        case WAVE_EQUATION: wave = true; break;
        case LINEAR_ELASTICITY: fea = true; break;
        case ADJ_EULER : euler = true; adj_euler = true; break;
        case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
        case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = ( ((config->GetKind_Adjoint() == CONTINUOUS) && (!config->GetFrozen_Visc())) || (config->GetKind_Adjoint() == HYBRID) ); break;
        case ADJ_PLASMA_EULER : plasma_euler = true; adj_plasma_euler = true; break;
        case ADJ_PLASMA_NAVIER_STOKES : plasma_ns = true; adj_plasma_ns = true; break;
        case ADJ_FREE_SURFACE_EULER: euler = true; levelset = true; adj_euler = true; adj_levelset = true; break;
        case ADJ_FREE_SURFACE_NAVIER_STOKES: ns = true; levelset = true; adj_ns = true; adj_levelset = true; break;
        case ADJ_FREE_SURFACE_RANS : ns = true; turbulent = true; levelset = true; adj_ns = true; adj_turb = (((config->GetKind_Adjoint() != HYBRID) && (!config->GetFrozen_Visc())) || (config->GetKind_Adjoint() == HYBRID)); adj_levelset = true; break;
        case LIN_EULER: euler = true; lin_euler = true; break;
            
            /*--- Specify by zone for the aeroacoustic problem ---*/
        case AEROACOUSTIC_EULER:
            if (iZone == ZONE_0) {
                euler = true;
            } else if (iZone == ZONE_1) {
                wave = true;
            }
            break;
        case ADJ_AEROACOUSTIC_EULER:
            if (iZone == ZONE_0) {
                euler = true; adj_euler = true;
            } else if (iZone == ZONE_1) {
                wave = true;
            }
            break;
        case PLASMA_EULER:
            if (iZone == ZONE_0) {
                plasma_euler = true;
            } else if (iZone == ZONE_1) {
                electric = true;
            }
            break;
        case PLASMA_NAVIER_STOKES:
            if (iZone == ZONE_0) {
                plasma_ns = true;
            } else if (iZone == ZONE_1) {
                electric = true;
            }
            break;
	}
    
	/*--- Assign turbulence model booleans --- */
	if (turbulent) {
		switch (config->GetKind_Turb_Model()) {
      case SA: spalart_allmaras = true; break;
      case SST: menter_sst = true; break;
      default: cout << "Specified turbulence model unavailable or none selected" << endl; cin.get(); break;
		}
  }
  
	if (plasma_euler || plasma_ns) {
		switch (config->GetKind_GasModel()) {
      case AIR7: plasma_diatomic = true; break;
      case O2: plasma_diatomic = true; break;
      case N2: plasma_diatomic = true; break;
      case AIR5: plasma_diatomic = true; break;
      case ARGON: plasma_monatomic = true; break;
      case AIR21: plasma_diatomic = true; break;
      case ARGON_SID: plasma_diatomic = true; break;
      default: cout << "Specified plasma model unavailable or none selected" << endl; cin.get(); break;
		}
	}

	/*--- Allocate solution for a template problem ---*/
	if (template_solver) integration_container[TEMPLATE_SOL] = new CSingleGridIntegration(config);

	/*--- Allocate solution for direct problem ---*/
	if (euler) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
	if (ns) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
	if (turbulent) integration_container[TURB_SOL] = new CSingleGridIntegration(config);
	if (transition) integration_container[TRANS_SOL] = new CSingleGridIntegration(config);
	if (electric) integration_container[ELEC_SOL] = new CPotentialIntegration(config);
	if (plasma_euler) integration_container[PLASMA_SOL] = new CMultiGridIntegration(config);
	if (plasma_ns) integration_container[PLASMA_SOL] = new CMultiGridIntegration(config);
	if (levelset) integration_container[LEVELSET_SOL] = new CSingleGridIntegration(config);
	if (wave) integration_container[WAVE_SOL] = new CSingleGridIntegration(config);
	if (fea) integration_container[FEA_SOL] = new CSingleGridIntegration(config);

	/*--- Allocate solution for adjoint problem ---*/
	if (adj_euler) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
	if (adj_ns) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
	if (adj_turb) integration_container[ADJTURB_SOL] = new CSingleGridIntegration(config);
	if (adj_plasma_euler) integration_container[ADJPLASMA_SOL] = new CMultiGridIntegration(config);
	if (adj_plasma_ns) integration_container[ADJPLASMA_SOL] = new CMultiGridIntegration(config);
	if (adj_levelset) integration_container[ADJLEVELSET_SOL] = new CSingleGridIntegration(config);

	/*--- Allocate solution for linear problem (at the moment we use the same scheme as the adjoint problem) ---*/
	if (lin_euler) integration_container[LINFLOW_SOL] = new CMultiGridIntegration(config);
	if (lin_ns) { cout <<"Equation not implemented." << endl; cin.get(); }

}


void Numerics_Preprocessing(CNumerics ****numerics_container, CSolver ***solver_container, CGeometry **geometry, 
		CConfig *config, unsigned short iZone) {

	unsigned short iMGlevel, iSol, nDim, nVar_Template = 0, nVar_Flow = 0, nVar_Trans = 0, nVar_Adj_Flow = 0, nVar_Plasma = 0,
			nVar_LevelSet = 0, nVar_Turb = 0, nVar_Adj_Turb = 0, nVar_Elec = 0, nVar_FEA = 0, nVar_Wave = 0, nVar_Lin_Flow = 0,
			nVar_Adj_LevelSet = 0, nVar_Adj_Plasma = 0, nSpecies = 0, nDiatomics = 0, nMonatomics = 0;

	double *constants = NULL;

	bool euler, ns, plasma_euler, plasma_ns, plasma_monatomic, plasma_diatomic,
	levelset, adj_plasma_euler, adj_plasma_ns, adj_levelset, adj_euler, lin_euler, adj_ns, lin_ns, turbulent,
	adj_turb, lin_turb, electric, wave, fea, spalart_allmaras, menter_sst, template_solver, transition;

	bool incompressible = config->GetIncompressible();

	/*--- Initialize some useful booleans ---*/
	euler          = false;   ns               = false;   turbulent        = false;
	electric       = false;   plasma_monatomic = false;   plasma_diatomic  = false;  levelset         = false;   plasma_euler     = false;
	plasma_ns      = false;   adj_euler        = false;	  adj_ns           = false;	 adj_turb         = false;
	wave           = false;   fea              = false;   adj_levelset     = false;	 spalart_allmaras = false; 
	lin_euler      = false;   lin_ns           = false;   lin_turb         = false;	 menter_sst       = false;   adj_plasma_euler = false;
	adj_plasma_ns  = false;   transition       = false;   template_solver  = false;

	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
	case TEMPLATE_SOLVER: template_solver = true; break;
	case EULER : euler = true; break;
	case NAVIER_STOKES: ns = true; break;
	case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
	case FREE_SURFACE_EULER: euler = true; levelset = true; break;
	case FREE_SURFACE_NAVIER_STOKES: ns = true; levelset = true; break;
	case FREE_SURFACE_RANS: ns = true; turbulent = true; levelset = true; break;
	case FLUID_STRUCTURE_EULER: euler = true; fea = true; break;
	case FLUID_STRUCTURE_NAVIER_STOKES: ns = true; fea = true; break;
	case FLUID_STRUCTURE_RANS: ns = true; turbulent = true; fea = true; break;
	case AEROACOUSTIC_NAVIER_STOKES: ns = true; wave = true; break;
	case AEROACOUSTIC_RANS: ns = true; turbulent = true; wave = true; break;
	case ELECTRIC_POTENTIAL: electric = true; break;
	case WAVE_EQUATION: wave = true; break;
	case LINEAR_ELASTICITY: fea = true; break;
	case ADJ_EULER : euler = true; adj_euler = true; break;
	case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
	case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = ( ((config->GetKind_Adjoint() == CONTINUOUS) && (!config->GetFrozen_Visc())) || (config->GetKind_Adjoint() == HYBRID) ); break;
	case ADJ_FREE_SURFACE_EULER: euler = true; adj_euler = true; levelset = true; adj_levelset = true; break;
	case ADJ_FREE_SURFACE_NAVIER_STOKES: ns = true; adj_ns = true; levelset = true; adj_levelset = true; break;
	case ADJ_FREE_SURFACE_RANS: ns = true; adj_ns = true; turbulent = true; adj_turb = (((config->GetKind_Adjoint() != HYBRID) && (!config->GetFrozen_Visc())) || (config->GetKind_Adjoint() == HYBRID)); levelset = true; adj_levelset = true; break;
	case ADJ_PLASMA_EULER : plasma_euler = true; adj_plasma_euler = true; break;
	case ADJ_PLASMA_NAVIER_STOKES : plasma_ns = true; adj_plasma_ns = true; break;
	case LIN_EULER: euler = true; lin_euler = true; break;

	/*--- Specify by zone for the aeroacoustic problem ---*/
	case AEROACOUSTIC_EULER:
		if (iZone == ZONE_0) {
			euler = true;
		} else if (iZone == ZONE_1) {
			wave = true;
		}
		break;
	case ADJ_AEROACOUSTIC_EULER:
		if (iZone == ZONE_0) {
			euler = true; adj_euler = true;
		} else if (iZone == ZONE_1) {
			wave = true;
		}
		break;
	case PLASMA_EULER:
		if (iZone == ZONE_0) {
			plasma_euler = true;
		} else if (iZone == ZONE_1) {
			electric = true;
		}
		break;
	case PLASMA_NAVIER_STOKES:
		if (iZone == ZONE_0) {
			plasma_ns = true;
		} else if (iZone == ZONE_1) {
			electric = true;
		}
		break;
	}

	/*--- Assign turbulence model booleans --- */
	if (turbulent)
		switch (config->GetKind_Turb_Model()){
		case SA: spalart_allmaras = true; break;
		case SST: menter_sst = true; constants = solver_container[MESH_0][TURB_SOL]->GetConstants(); break;
		default: cout << "Specified turbulence model unavailable or none selected" << endl; cin.get(); break;
		}

	if (plasma_euler || plasma_ns) {
		switch (config->GetKind_GasModel()){
		case AIR7: plasma_diatomic = true; break;
		case O2: plasma_diatomic = true; break;
		case N2: plasma_diatomic = true; break;
		case AIR5: plasma_diatomic = true; break;
		case ARGON: plasma_monatomic = true; break;
		case AIR21: plasma_diatomic = true; break;
		case ARGON_SID: plasma_diatomic = true; break;
		default: cout << "Specified plasma model unavailable or none selected" << endl; cin.get(); break;
		}
		//if (config->GetElectricSolver()) electric  = true;
	}

	/*--- Number of variables for the template ---*/
	if (template_solver) nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();

	/*--- Number of variables for direct problem ---*/
	if (euler)				nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
	if (ns)	nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
	if (turbulent)		nVar_Turb = solver_container[MESH_0][TURB_SOL]->GetnVar();
	if (transition)		nVar_Trans = solver_container[MESH_0][TRANS_SOL]->GetnVar();
	if (electric)			nVar_Elec = solver_container[MESH_0][ELEC_SOL]->GetnVar();
	if (plasma_euler || plasma_ns)	{ 
		nVar_Plasma = solver_container[MESH_0][PLASMA_SOL]->GetnVar();
		nSpecies    = solver_container[MESH_0][PLASMA_SOL]->GetnSpecies();
		nDiatomics  = solver_container[MESH_0][PLASMA_SOL]->GetnDiatomics();
		nMonatomics = solver_container[MESH_0][PLASMA_SOL]->GetnMonatomics();
	}
	if (levelset)		nVar_LevelSet = solver_container[MESH_0][LEVELSET_SOL]->GetnVar();
	if (wave)				nVar_Wave = solver_container[MESH_0][WAVE_SOL]->GetnVar();
	if (fea)				nVar_FEA = solver_container[MESH_0][FEA_SOL]->GetnVar();

	/*--- Number of variables for adjoint problem ---*/
	if (adj_euler)  	nVar_Adj_Flow = solver_container[MESH_0][ADJFLOW_SOL]->GetnVar();
	if (adj_ns)			  nVar_Adj_Flow = solver_container[MESH_0][ADJFLOW_SOL]->GetnVar();
	if (adj_turb)		  nVar_Adj_Turb = solver_container[MESH_0][ADJTURB_SOL]->GetnVar();
	if (adj_levelset)	nVar_Adj_LevelSet = solver_container[MESH_0][ADJLEVELSET_SOL]->GetnVar();
	if (adj_plasma_euler || adj_plasma_ns)		nVar_Adj_Plasma = solver_container[MESH_0][ADJPLASMA_SOL]->GetnVar();

	/*--- Number of variables for the linear problem ---*/
	if (lin_euler)	nVar_Lin_Flow = solver_container[MESH_0][LINFLOW_SOL]->GetnVar();

	/*--- Number of dimensions ---*/
	nDim = geometry[MESH_0]->GetnDim();

	/*--- Definition of the Class for the numerical method: numerics_container[MESH_LEVEL][EQUATION][EQ_TERM] ---*/
	for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
		numerics_container[iMGlevel] = new CNumerics** [MAX_SOLS];
		for (iSol = 0; iSol < MAX_SOLS; iSol++)
			numerics_container[iMGlevel][iSol] = new CNumerics* [MAX_TERMS];
	}

	/*--- Solver definition for the template problem ---*/
	if (template_solver) {

		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Template()) {
		case SPACE_CENTERED : case SPACE_UPWIND :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
				numerics_container[iMGlevel][TEMPLATE_SOL][CONV_TERM] = new CConvective_Template(nDim, nVar_Template, config);
			break;
		default : cout << "Convective scheme not implemented (template_solver)." << endl; cin.get(); break;
		}

		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Template()) {
		case AVG_GRAD : case AVG_GRAD_CORRECTED : case GALERKIN :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
				numerics_container[iMGlevel][TEMPLATE_SOL][VISC_TERM] = new CViscous_Template(nDim, nVar_Template, config);
			break;
		default : cout << "Viscous scheme not implemented." << endl; cin.get(); break;
		}

		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Template()) {
		case PIECEWISE_CONSTANT :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
				numerics_container[iMGlevel][TEMPLATE_SOL][SOURCE_FIRST_TERM] = new CSource_Template(nDim, nVar_Template, config);
			break;
		default : cout << "Source term not implemented." << endl; cin.get(); break;
		}

		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
			numerics_container[iMGlevel][TEMPLATE_SOL][CONV_BOUND_TERM] = new CConvective_Template(nDim, nVar_Template, config);	
		}

	}

	/*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
	if ((euler) || (ns)) {

		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Flow()) {
		case NO_CONVECTIVE :
			cout << "No convective scheme." << endl; cin.get();
			break;

		case SPACE_CENTERED :
			if (incompressible) {
				/*--- Incompressible flow, use artificial compressibility method ---*/
				switch (config->GetKind_Centered_Flow()) {
				case NO_CENTERED : cout << "No centered scheme." << endl; break;
				case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim, nVar_Flow, config); break;
				case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJSTArtComp_Flow(nDim, nVar_Flow, config); break;
				default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
				}
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim,nVar_Flow, config);

				/*--- Definition of the boundary condition method ---*/
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config);

			}
			else {
				/*--- Compressible flow ---*/
				switch (config->GetKind_Centered_Flow()) {
				case NO_CENTERED : cout << "No centered scheme." << endl; break;
				case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim,nVar_Flow, config); break;
				case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim,nVar_Flow, config); break;
				default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
				}

        if (!config->GetLowFidelitySim()) {
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
        }
        else {
          numerics_container[MESH_1][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim, nVar_Flow, config);
          for (iMGlevel = 2; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
        }

				/*--- Definition of the boundary condition method ---*/
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);

			}
			break;
		case SPACE_UPWIND :
			if (incompressible) {
				/*--- Incompressible flow, use artificial compressibility method ---*/
				switch (config->GetKind_Upwind_Flow()) {
				case NO_UPWIND : cout << "No upwind scheme." << endl; break;
				case ROE_1ST : case ROE_2ND :
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
						numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config);
						numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config);
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
						numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
						numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
					}
					break;

				case AUSM_1ST : case AUSM_2ND :
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
						numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
						numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
					}
					break;

				case ROE_TURKEL_1ST : case ROE_TURKEL_2ND :
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
						numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoe_Turkel_Flow(nDim, nVar_Flow, config);
						numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Turkel_Flow(nDim, nVar_Flow, config);
					}
					break;

				case HLLC_1ST : case HLLC_2ND :
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
						numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
						numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
					}
					break;

				default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
				}

			}
			break;

		default :
			cout << "Convective scheme not implemented (euler and ns)." << endl; cin.get();
			break;
		}

		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Flow()) {
		case NONE :
			break;
		case AVG_GRAD :
			if (incompressible) {
				/*--- Incompressible flow, use artificial compressibility method ---*/
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
          numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
        }
			}
			else {
				/*--- Compressible flow ---*/
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
          numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
        }
			}
			break;
		case AVG_GRAD_CORRECTED :
			if (incompressible) {
				/*--- Incompressible flow, use artificial compressibility method ---*/
				numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_Flow(nDim, nVar_Flow, config);
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
        
        /*--- Definition of the boundary condition method ---*/
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
			}
			else {
				/*--- Compressible flow ---*/
				numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrected_Flow(nDim, nVar_Flow, config);
				for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
        
        /*--- Definition of the boundary condition method ---*/
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
					numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
			}
			break;
		case GALERKIN :
			cout << "Galerkin viscous scheme not implemented." << endl; cin.get(); exit(1);
			break;
		default :
			cout << "Numerical viscous scheme not recognized." << endl; cin.get(); exit(1);
			break;
		}

		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Flow()) {
		case NONE :
			break;
		case PIECEWISE_CONSTANT :

			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {

				if (config->GetRotating_Frame() == YES)
					numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceRotatingFrame_Flow(nDim, nVar_Flow, config);
				else if (config->GetAxisymmetric() == YES)
					numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceAxisymmetric_Flow(nDim,nVar_Flow, config);
				else if (config->GetGravityForce() == YES)
					numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_Gravity(nDim, nVar_Flow, config);
				else if (config->GetMagnetic_Force() == YES)
					numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSource_Magnet(nDim, nVar_Flow, config);
				else if (config->GetJouleHeating() == YES)
					numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSource_JouleHeating(nDim, nVar_Flow, config);
				else
					numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Flow, config);

				numerics_container[iMGlevel][FLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
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
				if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
				else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
			}
			break;
		default :
			cout << "Convective scheme not implemented (turbulent)." << endl; cin.get();
			break;
		}

		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Turb()) {	
		case NONE :
			break;
		case AVG_GRAD :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
				if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, config);
				else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSST(nDim, nVar_Turb, config);
			}
			break;
		case AVG_GRAD_CORRECTED :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
				if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSA(nDim, nVar_Turb, config);
				else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSST(nDim, nVar_Turb, constants, config);
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
				if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA(nDim, nVar_Turb, config);
				else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSST(nDim, nVar_Turb, constants, config);
				numerics_container[iMGlevel][TURB_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Turb, config);
			}
			break;
		default :
			cout << "Source term not implemented." << endl; cin.get();
			break;
		}

		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
			if (spalart_allmaras) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, config);
      }
			else if (menter_sst) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSST(nDim, nVar_Turb, config);
      }
		}
	}

	/*--- Solver definition for the transition model problem ---*/
	if (transition) {

		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Turb()) {
		case NONE :
			break;
		case SPACE_UPWIND :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
				numerics_container[iMGlevel][TRANS_SOL][CONV_TERM] = new CUpwSca_TransLM(nDim, nVar_Trans, config);
			}
			break;
		default :
			cout << "Convective scheme not implemented (transition)." << endl; cin.get();
			break;
		}

		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Turb()) {	
		case NONE :
			break;
		case AVG_GRAD :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
				numerics_container[iMGlevel][TRANS_SOL][VISC_TERM] = new CAvgGrad_TransLM(nDim, nVar_Trans, config);
			}
			break;
		case AVG_GRAD_CORRECTED :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
				numerics_container[iMGlevel][TRANS_SOL][VISC_TERM] = new CAvgGradCorrected_TransLM(nDim, nVar_Trans, config);
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
				numerics_container[iMGlevel][TRANS_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TransLM(nDim, nVar_Trans, config);
				numerics_container[iMGlevel][TRANS_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Trans, config);
			}
			break;
		default :
			cout << "Source term not implemented." << endl; cin.get();
			break;
		}

		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
			numerics_container[iMGlevel][TRANS_SOL][CONV_BOUND_TERM] = new CUpwLin_TransLM(nDim, nVar_Trans, config);
		}
	}

	/*--- Solver definition for the multi species plasma model problem ---*/
	if (plasma_euler || plasma_ns) {

		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Plasma()) {
		case NONE :
			break;
		case SPACE_UPWIND :

			switch (config->GetKind_Upwind_Plasma()) {
			case NO_UPWIND : cout << "No upwind scheme." << endl; break;
			case ROE_1ST : case ROE_2ND :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					switch (config->GetKind_GasModel()) {
					case ARGON:
						numerics_container[iMGlevel][PLASMA_SOL][CONV_TERM]  = new CUpwRoe_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
						numerics_container[iMGlevel][PLASMA_SOL][CONV_BOUND_TERM] = new CUpwRoe_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
						break;
					case O2: case N2:	case AIR5:	case AIR7: case ARGON_SID:
						numerics_container[iMGlevel][PLASMA_SOL][CONV_TERM]  = new CUpwRoe_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
						numerics_container[iMGlevel][PLASMA_SOL][CONV_BOUND_TERM] = new CUpwRoe_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
						break;
					}
				}
				break;
			case HLLC_1ST :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					switch (config->GetKind_GasModel()) {
					case O2 : case N2 : case AIR5 : case AIR7 : case ARGON_SID:
						numerics_container[iMGlevel][PLASMA_SOL][CONV_TERM]  = new CUpwHLLC_PlasmaDiatomic(nDim, nVar_Flow, config);
						numerics_container[iMGlevel][PLASMA_SOL][CONV_BOUND_TERM] = new CUpwHLLC_PlasmaDiatomic(nDim, nVar_Flow, config);
						break;
					default:
						cout << "HLLC Upwind scheme not implemented for the selected gas chemistry model..." << endl; cin.get(); break;
					}
				}
				break;
			case ROE_TURKEL_1ST : case ROE_TURKEL_2ND :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					numerics_container[iMGlevel][PLASMA_SOL][CONV_TERM]  = new CUpwRoe_Turkel_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
					numerics_container[iMGlevel][PLASMA_SOL][CONV_BOUND_TERM] = new CUpwRoe_Turkel_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				}
				break;
			case SW_1ST : case SW_2ND :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					switch (config->GetKind_GasModel()) {
					case O2 : case N2 : case AIR5 : case AIR7 : case ARGON_SID:
						numerics_container[iMGlevel][PLASMA_SOL][CONV_TERM]  = new CUpwSW_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
						numerics_container[iMGlevel][PLASMA_SOL][CONV_BOUND_TERM] = new CUpwSW_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
						break;
					default:
						cout << "Steger-Warming Upwind scheme not implemented for the selected gas chemistry model..." << endl; cin.get(); break;
					}
				}
				break;
      case MSW_1ST : case MSW_2ND :
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            switch (config->GetKind_GasModel()) {
              case O2 : case N2 : case AIR5 : case AIR7 : case ARGON_SID:
                numerics_container[iMGlevel][PLASMA_SOL][CONV_TERM]  = new CUpwMSW_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
                numerics_container[iMGlevel][PLASMA_SOL][CONV_BOUND_TERM] = new CUpwMSW_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
                break;
              default:
                cout << "Modified Steger-Warming Upwind scheme not implemented for the selected gas chemistry model..." << endl; cin.get(); break;
            }
          }
          break;
			default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
			}
			break;

			case SPACE_CENTERED :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					switch (config->GetKind_Centered_Plasma()) {
					case JST:
						switch (config->GetKind_GasModel()) {
						case ARGON:
							numerics_container[iMGlevel][PLASMA_SOL][CONV_TERM] = new CCentJST_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
							numerics_container[iMGlevel][PLASMA_SOL][CONV_BOUND_TERM] = new CUpwRoe_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
							break;
						case O2: case N2: case AIR5: case AIR7: case ARGON_SID:
							numerics_container[iMGlevel][PLASMA_SOL][CONV_TERM] = new CCentJST_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
							numerics_container[iMGlevel][PLASMA_SOL][CONV_BOUND_TERM] = new CUpwRoe_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
							break;
						default:
							break;
						}
						break;

						case LAX:
							switch (config->GetKind_GasModel()) {
							case ARGON:
								cout << "Not implemented..." << endl; cin.get();
								break;
							case O2: case N2: case AIR5: case AIR7: case ARGON_SID:
								numerics_container[iMGlevel][PLASMA_SOL][CONV_TERM] = new CCentLax_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
								numerics_container[iMGlevel][PLASMA_SOL][CONV_BOUND_TERM] = new CUpwRoe_PlasmaDiatomic(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
								break;
							default:
								break;
							}
							break;
					}
				}
				break;
		}

		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		if (plasma_ns) {
			switch (config->GetKind_ViscNumScheme_Plasma()) {
			case NONE :
				break;
			case AVG_GRAD :
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][PLASMA_SOL][VISC_TERM] = new CAvgGrad_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics,config);
            numerics_container[iMGlevel][PLASMA_SOL][VISC_BOUND_TERM] = new CAvgGrad_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics,config);
          }
				break;
			case AVG_GRAD_CORRECTED :
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][PLASMA_SOL][VISC_TERM] = new CAvgGradCorrected_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics,config);
            numerics_container[iMGlevel][PLASMA_SOL][VISC_BOUND_TERM] = new CAvgGrad_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics,config);
          }
				break;
			default :
				cout << "Viscous scheme not implemented." << endl; cin.get();
				break;
			}
		}

		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Plasma()) {
		case NONE :
			break;
		case PIECEWISE_CONSTANT :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
				if (config->GetKind_GasModel() == ARGON)
					numerics_container[iMGlevel][PLASMA_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				if (config->GetKind_GasModel() == AIR7) {
					numerics_container[iMGlevel][PLASMA_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				}
				if (config->GetKind_GasModel() == AIR21)
					cout << "ERROR: 21 Species air gas chemistry model not implemented!!!" << endl;
				if (config->GetKind_GasModel() == O2)
					numerics_container[iMGlevel][PLASMA_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				if (config->GetKind_GasModel() == N2)
					numerics_container[iMGlevel][PLASMA_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				if (config->GetKind_GasModel() == AIR5) 
					numerics_container[iMGlevel][PLASMA_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				if (config->GetKind_GasModel() == ARGON_SID)
					numerics_container[iMGlevel][PLASMA_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_Plasma(nDim, nVar_Plasma, nSpecies, nDiatomics, nMonatomics, config);



				numerics_container[iMGlevel][PLASMA_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Plasma, config);
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
			numerics_container[MESH_0][ELEC_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Elec, config);
			break;
		default : cout << "Viscous scheme not implemented." << endl; cin.get(); break;
		}

		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Elec()) {
		case NONE :
			break;
		case PIECEWISE_CONSTANT :
			numerics_container[MESH_0][ELEC_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_Elec(nDim, nVar_Elec, config);
			numerics_container[MESH_0][ELEC_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Elec, config);
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
      case NONE :
        break;
			case SPACE_UPWIND :
				switch (config->GetKind_Upwind_LevelSet()) {
          case SCALAR_UPWIND_1ST : case SCALAR_UPWIND_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][LEVELSET_SOL][CONV_TERM] = new CUpwLin_LevelSet(nDim, nVar_LevelSet, config);
            }
            break;
          default :
            cout << "Upwind scheme not implemented." << endl; cin.get();
            break;
				}
				break;
      default : cout << "Convective scheme not implemented (levelset)." << endl; cin.get(); break;
		}
    
		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
			numerics_container[iMGlevel][LEVELSET_SOL][CONV_BOUND_TERM] = new CUpwLin_LevelSet(nDim, nVar_LevelSet, config);
		}
    
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_LevelSet()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][LEVELSET_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_LevelSet(nDim, nVar_LevelSet, config);
        break;
      default :
        cout << "Source term not implemented." << endl; cin.get();
        break;
		}
	}

	/*--- Solver definition for the flow adjoint problem ---*/
	if (adj_euler || adj_ns) {
    
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_AdjFlow()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; cin.get();
        break;
      case SPACE_CENTERED :
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Centered_AdjFlow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLaxArtComp_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            case JST : cout << "Centered scheme not implemented." << endl; cin.get(); break;
            default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLaxArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
          
        }
        else {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Centered_AdjFlow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            case JST : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentJST_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
        }
        break;
      case SPACE_UPWIND :
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_AdjFlow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE_1ST : case ROE_2ND :
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
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
                  numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
                  numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
          }
        }
        break;
        
      default :
        cout << "Convective scheme not implemented (adj_euler and adj_ns)." << endl; cin.get();
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
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
        }
        else {
          /*--- Compressible flow ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
          }
        }
        break;
      case AVG_GRAD_CORRECTED :
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          numerics_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
        }
        else {
          /*--- Compressible flow ---*/
          numerics_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrected_AdjFlow(nDim, nVar_Adj_Flow, config);
          numerics_container[MESH_0][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
          }
          
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
          
          /*--- Note that RANS is incompatible with Axisymmetric or Rotational (Fix it!) ---*/ 
          if ((adj_ns) && (!incompressible)) {
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceViscous_AdjFlow(nDim, nVar_Adj_Flow, config);
            if (config->GetRotating_Frame() == YES)
              numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
            else
              numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceConservative_AdjFlow(nDim, nVar_Adj_Flow, config);

          }
          else {
            if (config->GetRotating_Frame() == YES)
              numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
            else
              numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);

            if (config->GetAxisymmetric() == YES)
              numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceAxisymmetric_AdjFlow(nDim, nVar_Adj_Flow, config);
            else
              numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
            
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
          }
        
        }
        break;
      default :
        cout << "Source term not implemented." << endl; cin.get();
        break;
		}
	}

	/*--- Solver definition for the multi species plasma model problem ---*/
	if (adj_plasma_euler || adj_plasma_ns) {

		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_AdjPlasma()) {
		case NONE :
			break;
		case SPACE_UPWIND :

			switch (config->GetKind_Upwind_AdjPlasma()) {
			case NO_UPWIND : cout << "No upwind scheme." << endl; break;
			case ROE_1ST :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					numerics_container[iMGlevel][ADJPLASMA_SOL][CONV_TERM] = new CUpwRoe_AdjPlasmaDiatomic(nDim, nVar_Adj_Plasma, nSpecies, nDiatomics, nMonatomics, config);
					numerics_container[iMGlevel][ADJPLASMA_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjPlasmaDiatomic(nDim, nVar_Adj_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				}
				break;
			case SW_1ST:
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					numerics_container[iMGlevel][ADJPLASMA_SOL][CONV_TERM] = new CUpwSW_AdjPlasmaDiatomic(nDim, nVar_Adj_Plasma, nSpecies, nDiatomics, nMonatomics, config);
					numerics_container[iMGlevel][ADJPLASMA_SOL][CONV_BOUND_TERM] = new CUpwSW_AdjPlasmaDiatomic(nDim, nVar_Adj_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				}
				break;
			}
			break;
			case SPACE_CENTERED :
				for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
					numerics_container[iMGlevel][ADJPLASMA_SOL][CONV_TERM] = new CCentLax_AdjPlasmaDiatomic(nDim, nVar_Adj_Plasma, nSpecies, nDiatomics, nMonatomics, config);
					numerics_container[iMGlevel][ADJPLASMA_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjPlasmaDiatomic(nDim, nVar_Adj_Plasma, nSpecies, nDiatomics, nMonatomics, config);
				}
				break;
			default :
				cout << "Convective scheme not implemented (adj_plasma_euler and adj_plasma_ns)." << endl; cin.get();
				break;
		}


		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_AdjPlasma()) {
		case NONE :
			break;
		case PIECEWISE_CONSTANT :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
				numerics_container[iMGlevel][ADJPLASMA_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_Plasma(nDim, nVar_Adj_Plasma, nSpecies, nDiatomics, nMonatomics, config);
			}
			break;
		default :
			cout << "Source term not implemented." << endl; cin.get();
			break;
		}
	}

	/*--- Solver definition for the linearized flow problem ---*/
	if (lin_euler) {

		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_LinFlow()) {
		case NONE :
			break;
		case SPACE_CENTERED :
			switch (config->GetKind_Centered_LinFlow()) {
			case LAX : numerics_container[MESH_0][LINFLOW_SOL][CONV_TERM] = new CCentLax_LinFlow(nDim, nVar_Lin_Flow, config); break;
			case JST : numerics_container[MESH_0][LINFLOW_SOL][CONV_TERM] = new CCentJST_LinFlow(nDim, nVar_Lin_Flow, config); break;
			default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
			}
			for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
				numerics_container[iMGlevel][LINFLOW_SOL][CONV_TERM] = new CCentLax_LinFlow(nDim, nVar_Lin_Flow, config);
			break;
			default :
				cout << "Convective scheme not implemented (lin_euler)." << endl; cin.get();
				break;
		}

		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
			numerics_container[iMGlevel][LINFLOW_SOL][CONV_BOUND_TERM] = new CCentLax_LinFlow(nDim, nVar_Lin_Flow, config);
	}

	/*--- Solver definition for the turbulent adjoint problem ---*/
	if (adj_turb) {
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_AdjTurb()) {
		case NONE :
			break;
		case SPACE_UPWIND :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
				if (spalart_allmaras) {
						numerics_container[iMGlevel][ADJTURB_SOL][CONV_TERM] = new CUpwSca_AdjTurb(nDim, nVar_Adj_Turb, config);
				}
				else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; cin.get();}
			break;
		default :
			cout << "Convective scheme not implemented (adj_turb)." << endl; cin.get();
			break;
		}

		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_AdjTurb()) {
		case NONE :
			break;
		case AVG_GRAD :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
			if (spalart_allmaras){
				cout << "Viscous scheme not implemented." << endl; cin.get();
			}
			else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; cin.get();}
			break;
		case AVG_GRAD_CORRECTED :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
				if (spalart_allmaras){
					numerics_container[iMGlevel][ADJTURB_SOL][VISC_TERM] = new CAvgGradCorrected_AdjTurb(nDim, nVar_Adj_Turb, config);
				}
				else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; cin.get();}
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
				if (spalart_allmaras) {
						numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_AdjTurb(nDim, nVar_Adj_Turb, config);
						numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_SECOND_TERM] = new CSourceConservative_AdjTurb(nDim, nVar_Adj_Turb, config);
				}
				else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; cin.get();}
			}
			break;
		default :
			cout << "Source term not implemented." << endl; cin.get();
			break;
		}

		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
			if (spalart_allmaras) numerics_container[iMGlevel][ADJTURB_SOL][CONV_BOUND_TERM] = new CUpwLin_AdjTurb(nDim, nVar_Adj_Turb, config);
			else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; cin.get();}
		}

	}

	/*--- Solver definition for the level set model problem ---*/
	if (adj_levelset) {

		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_AdjLevelSet()) {
		case NO_CONVECTIVE : cout << "No convective scheme." << endl; cin.get(); break;
		case SPACE_CENTERED :
			switch (config->GetKind_Centered_AdjLevelSet()) {
			case NO_UPWIND : cout << "No centered scheme." << endl; cin.get(); break;
			default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
			}
			break;
			case SPACE_UPWIND :
				switch (config->GetKind_Upwind_AdjLevelSet()) {
				case NO_UPWIND : cout << "No upwind scheme." << endl; cin.get(); break;
				case SCALAR_UPWIND_1ST : case SCALAR_UPWIND_2ND :
					for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
						numerics_container[iMGlevel][ADJLEVELSET_SOL][CONV_TERM] = new CUpwLin_AdjLevelSet(nDim, nVar_Adj_LevelSet, config);
						numerics_container[iMGlevel][ADJLEVELSET_SOL][CONV_BOUND_TERM] = new CUpwLin_AdjLevelSet(nDim, nVar_Adj_LevelSet, config);
					}
					break;
				default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
				}
				break;
				default : cout << "Convective scheme not implemented (adj_levelset)." << endl; cin.get(); break;
		}

		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_AdjLevelSet()) {
		case NONE :
			break;
		case PIECEWISE_CONSTANT :
			for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
				numerics_container[iMGlevel][ADJLEVELSET_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_AdjLevelSet(nDim, nVar_Turb, config);
			break;
		default :
			cout << "Source term not implemented." << endl; cin.get();
			break;
		}
	}

	/*--- Solver definition for the wave problem ---*/
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
			numerics_container[MESH_0][WAVE_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Wave, config);
			break;
		default :
			cout << "Viscous scheme not implemented." << endl; cin.get();
			break;
		}

		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Wave()) {
		case NONE : break;
		case PIECEWISE_CONSTANT :
			break;
		default : break;
		}
	}

	/*--- Solver definition for the FEA problem ---*/
	if (fea) {

		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_FEA()) {
		case NONE :
			break;
		case AVG_GRAD :
			cout << "Viscous scheme not implemented." << endl; cin.get();
			break;
		case AVG_GRAD_CORRECTED :
			cout << "Viscous scheme not implemented." << endl; cin.get();
			break;
		case GALERKIN :
			numerics_container[MESH_0][FEA_SOL][VISC_TERM] = new CGalerkin_FEA(nDim, nVar_Wave, config);
			break;
		default :
			cout << "Viscous scheme not implemented." << endl; cin.get();
			break;
		}

		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Wave()) {
		case NONE : break;
		case PIECEWISE_CONSTANT :
			break;
		default : break;
		}
	}

}
