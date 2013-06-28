/*!
 * \file SU2_MAC.cpp
 * \brief Main file of Mesh Adaptation Code (SU2_MAC).
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

#include "../include/SU2_MAC.hpp"
using namespace std;

int main(int argc, char *argv[]) {
	
	/*--- Variable definitions ---*/
	char file_name[200];
  unsigned short nZone = 1;

	/*-- Definition of the Class for the definition of the problem ---*/
	CConfig *config;
	if (argc == 2) config = new CConfig(argv[1], SU2_MAC);
	else {
		strcpy (file_name, "default.cfg");
		config = new CConfig(file_name, SU2_MAC);
	}
	
	/*-- Definition of the Class for the geometry ---*/
	CGeometry *geometry; 
	geometry = new CGeometry;
	geometry = new CPhysicalGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat(), DOMAIN_0, nZone);
	
  /*--- Perform the non-dimensionalization, in case any values are needed ---*/
  config->SetNondimensionalization(geometry->GetnDim(), MASTER_NODE, nZone);
  
	cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
	
	/*--- Compute elements surrounding points, points surrounding points, and elements surronding elements ---*/
	cout << "Setting local point and element connectivity." <<endl; 
	geometry->SetEsuP(); geometry->SetPsuP(); geometry->SetEsuE();
	
	/*--- Check the orientation before computing geometrical quantities ---*/
	cout << "Check numerical grid orientation." <<endl; 
	geometry->SetBoundVolume(); geometry->Check_Orientation(config);
	
	/*--- Create the edge structure ---*/
	cout << "Identify edges and vertices." <<endl; 
	geometry->SetEdges(); geometry->SetVertex(config); geometry->SetCG();
	
	/*--- Create the control volume structures ---*/
	cout << "Set control volume structure." << endl; 
	geometry->SetControlVolume(config, ALLOCATE); geometry->SetBoundControlVolume(config, ALLOCATE);
	
	if (config->GetKind_Adaptation() != NONE) {
		
		cout << endl <<"--------------------- Start numerical grid adaptation -------------------" << endl;	
		
		/*-- Definition of the Class for grid adaptation ---*/
		CGridAdaptation *grid_adaptation; 
		grid_adaptation = new CGridAdaptation(geometry, config);
		
		/*--- Read the flow solution and/or the adjoint solution 
		 and choose the elements to adapt ---*/
		if ((config->GetKind_Adaptation() != NONE) && (config->GetKind_Adaptation() != FULL) 
				&& (config->GetKind_Adaptation() != WAKE) && (config->GetKind_Adaptation() != TWOPHASE) 
				&& (config->GetKind_Adaptation() != SMOOTHING) && (config->GetKind_Adaptation() != HORIZONTAL_PLANE) 
				&& (config->GetKind_Adaptation() != CURVE_SURFACE) && (config->GetKind_Adaptation() != DOUBLE_SURFACE)
				&& (config->GetKind_Adaptation() != SUPERSONIC_SHOCK)) 
			grid_adaptation->GetFlowSolution(geometry, config);
		
		switch (config->GetKind_Adaptation()) {
			case NONE:
				break;
			case SMOOTHING:
				config->SetSmoothNumGrid(true);
				grid_adaptation->SetNo_Refinement(geometry, 1);
				break;
			case FULL:
				grid_adaptation->SetComplete_Refinement(geometry, 1);
				break;
			case WAKE:
				grid_adaptation->SetWake_Refinement(geometry, 1);
				break;
			case TWOPHASE:
				grid_adaptation->SetTwoPhase_Refinement(geometry, 1);
				break;
			case SUPERSONIC_SHOCK:
				grid_adaptation->SetSupShock_Refinement(geometry, config);
				break;	
			case HORIZONTAL_PLANE:
				grid_adaptation->SetNearField_Refinement(geometry, config);
				break;
			case CURVE_SURFACE: case DOUBLE_SURFACE:
				break;
			case FULL_FLOW: 
				grid_adaptation->SetComplete_Refinement(geometry, 1);
				break;
			case FULL_ADJOINT:
				grid_adaptation->GetAdjSolution(geometry, config);
				grid_adaptation->SetComplete_Refinement(geometry, 1);
				break;
			case FULL_LINEAR:
				grid_adaptation->GetLinSolution(geometry, config);
				grid_adaptation->SetComplete_Refinement(geometry, 1);
				break;			
			case GRAD_FLOW:
				grid_adaptation->SetIndicator_Flow(geometry, config, 1);
				break;
			case GRAD_ADJOINT:
				grid_adaptation->GetAdjSolution(geometry, config);
				grid_adaptation->SetIndicator_Adj(geometry, config, 1);
				break;
			case GRAD_FLOW_ADJ:
				grid_adaptation->GetAdjSolution(geometry, config);
				grid_adaptation->SetIndicator_FlowAdj(geometry, config);
				break;
			case COMPUTABLE:
				grid_adaptation->GetAdjSolution(geometry, config);
				grid_adaptation->GetFlowResidual(geometry, config);
				grid_adaptation->SetIndicator_Computable(geometry, config);			
				break;
			case REMAINING:
				cout << "Adaptation method not implemented."<< endl;
				cout << "Press any key to exit..." << endl;
				cin.get();
				exit(1);	
				break;
			case ROBUST:
				grid_adaptation->GetFlowResidual(geometry, config);
				grid_adaptation->GetAdjResidual(geometry, config);
				grid_adaptation->SetIndicator_Robust(geometry, config);
				break;
			case COMPUTABLE_ROBUST:
				grid_adaptation->GetAdjSolution(geometry, config);
				grid_adaptation->GetLinResidual(geometry, config);
				grid_adaptation->SetIndicator_Computable_Robust(geometry, config);
				break;
				
			default :
				cout << "The adaptation is not defined" << endl;
		}
		
		/*--- Perform an homothetic adaptation of the grid ---*/
		CPhysicalGeometry *geo_adapt; geo_adapt = new CPhysicalGeometry;
		if ((config->GetKind_Adaptation() != CURVE_SURFACE) && 
				(config->GetKind_Adaptation() != DOUBLE_SURFACE)) {
			grid_adaptation->SetHomothetic_Adaptation(geometry, geo_adapt, config);
		}
		
		/*--- Identify the nearfield surface using the main lines of the grid ---*/
		if (config->GetKind_Adaptation() == CURVE_SURFACE) {
			grid_adaptation->SetDomain_Interface(geometry, geo_adapt, config);
		}
		
		/*--- Smooth the numerical grid coordinates ---*/
		if (config->GetSmoothNumGrid()) {
			cout << "Preprocessing for doing the implicit smoothing." << endl;
			geo_adapt->SetEsuP(); geo_adapt->SetPsuP(); geo_adapt->SetEsuE();
			geo_adapt->SetEdges(); geo_adapt->SetVertex(config);
			cout << "Implicit smoothing of the numerical grid coordinates." << endl;
			geo_adapt->SetCoord_Smoothing(5, 1.5, config);
		}
		
		/*--- Original and adapted grid ---*/
		if (config->GetOutput_FileFormat() == PARAVIEW)  {
			strcpy (file_name, "original_grid.vtk");
			geometry->SetParaView(file_name);
		}
		if (config->GetOutput_FileFormat() == TECPLOT) {
			strcpy (file_name, "original_grid.plt");
			geometry->SetTecPlot(file_name); 
		}
		
		/*--- Write the adaptation sensor ---*/
		grid_adaptation->WriteAdaptSensor(geometry, file_name);
		
		if (config->GetOutput_FileFormat() == PARAVIEW) {
			strcpy (file_name, "adapted_grid.vtk");
			geo_adapt->SetParaView(file_name);
			strcpy (file_name, "adapted_surface.vtk");
			geo_adapt->SetBoundParaView(config,file_name);
		}
		if (config->GetOutput_FileFormat() == TECPLOT) {
			strcpy (file_name, "adapted_grid.plt");		
			geo_adapt->SetTecPlot(file_name);
			strcpy (file_name, "adapted_surface.plt");
			geo_adapt->SetBoundTecPlot(config,file_name);
		}
		
		/*--- Write adapted grid ---*/
		cout << "Write the mesh file." << endl;
		if ((config->GetKind_Adaptation() != HORIZONTAL_PLANE) 
				&& (config->GetKind_Adaptation() != CURVE_SURFACE) 
				&& (config->GetKind_Adaptation() != DOUBLE_SURFACE)) {
			/*--- Write the new adapted grid, including the modified boundaries surfaces ---*/
			geo_adapt->SetMeshFile(config, config->GetMesh_Out_FileName());
		}
		else {
			/*--- Write the new grid, including two nearfield boundaries ---*/
			if ((config->GetKind_Adaptation() == HORIZONTAL_PLANE) 
					&& (config->GetKind_Adaptation() == CURVE_SURFACE))
				geo_adapt->SetMeshFile_IntSurface(config, config->GetMesh_Out_FileName());
			if (config->GetKind_Adaptation() == DOUBLE_SURFACE)
				geometry->SetMeshFile_IntSurface(config, config->GetMesh_Out_FileName());
		}

		/*--- Write the restart file ---*/
		if ((config->GetKind_Adaptation() != SMOOTHING) && (config->GetKind_Adaptation() != FULL) && 
				(config->GetKind_Adaptation() != WAKE) && (config->GetKind_Adaptation() != TWOPHASE) && 
				(config->GetKind_Adaptation() != HORIZONTAL_PLANE) && (config->GetKind_Adaptation() != CURVE_SURFACE) && 
				(config->GetKind_Adaptation() != SUPERSONIC_SHOCK) && (config->GetKind_Adaptation() != DOUBLE_SURFACE))
			grid_adaptation->SetReStart_FlowSolution(geometry, config, config->GetReStart_FlowFileName());
		
		if ((config->GetKind_Adaptation() == GRAD_FLOW_ADJ) || (config->GetKind_Adaptation() == GRAD_ADJOINT) 
				|| (config->GetKind_Adaptation() == FULL_ADJOINT) || (config->GetKind_Adaptation() == ROBUST) 
				|| (config->GetKind_Adaptation() == COMPUTABLE) || (config->GetKind_Adaptation() == COMPUTABLE_ROBUST) || 
				(config->GetKind_Adaptation() == REMAINING))
			grid_adaptation->SetReStart_AdjSolution(geometry, config, config->GetReStart_AdjFileName());		
		
		if ((config->GetKind_Adaptation() == FULL_LINEAR) || (config->GetKind_Adaptation() == COMPUTABLE_ROBUST)) {
			grid_adaptation->SetReStart_LinSolution(geometry, config, config->GetReStart_LinFileName());
		}
	}
	else {
		if (config->GetOutput_FileFormat() == PARAVIEW) {
			strcpy (file_name, "original_grid.vtk");
			geometry->SetParaView(file_name);
		}
		if (config->GetOutput_FileFormat() == TECPLOT) {
			strcpy (file_name, "original_grid.plt");
			geometry->SetTecPlot(file_name);
		}
		geometry->SetMeshFile (config, config->GetMesh_Out_FileName());
	}

	
	return 1;
}

