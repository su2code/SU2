/*!
 * \file SU2_PBC.cpp
 * \brief Main file of Periodic Boundary Code (SU2_PBC).
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

#include "../include/SU2_PBC.hpp"
using namespace std;

int main(int argc, char *argv[]) {
	
	char buffer_vtk[100], buffer_plt[100]; 
	string MeshFile;
	unsigned short nZone = 1;
    
#ifndef NO_MPI
    MPI::Init(argc, argv);
#endif
  
	/*--- Definition of the class for the definition of the problem ---*/
	CConfig *config;
	if (argc == 2) config = new CConfig(argv[1], SU2_PBC, ZONE_0, nZone, VERB_HIGH);
	else {
		char grid_file[200];
		strcpy (grid_file, "default.cfg");
		config = new CConfig(grid_file, SU2_PBC, ZONE_0, nZone, VERB_HIGH);
	}

	/*--- Definition of the class for the geometry ---*/
	CGeometry *geometry; geometry = new CGeometry;
	geometry = new CPhysicalGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat(), ZONE_0, nZone);
  
  /*--- Perform the non-dimensionalization, in case any values are needed ---*/
  config->SetNondimensionalization(geometry->GetnDim(), ZONE_0);
  
	cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
	
	/*--- Compute elements surrounding points, points surrounding points, and elements surrounding elements ---*/
	cout << "Setting local point and element connectivity." <<endl; 
	geometry->SetEsuP(); geometry->SetPsuP(); geometry->SetEsuE();
	
	/*--- Check the orientation before computing geometrical quantities ---*/
	cout << "Checking the numerical grid orientation." <<endl;
	geometry->SetBoundVolume(); geometry->Check_Orientation(config);
	
	/*--- Create the edge structure ---*/
	cout << "Identifying edges and vertices." <<endl; 
	geometry->SetEdges(); geometry->SetVertex(config); 
	
	/*--- Compute center of gravity ---*/
	cout << "Computing centers of gravity." << endl;
	geometry->SetCG();
  
	/*--- Create the control volume structures ---*/
	cout << "Setting the control volume structure." << endl; 
	geometry->SetControlVolume(config, ALLOCATE);
	geometry->SetBoundControlVolume(config, ALLOCATE);
	
	cout << endl <<"-------------------- Setting the periodic boundaries --------------------" << endl;

	/*--- Set periodic boundary conditions ---*/
	geometry->SetPeriodicBoundary(config);
	
  /*--- Original grid for debugging purposes ---*/
  strcpy (buffer_plt, "periodic_original.plt"); geometry->SetTecPlot(buffer_plt);
 
	/*--- Create a new grid with the right periodic boundary ---*/
	CGeometry *periodic; periodic = new CPeriodicGeometry(geometry, config);
	periodic->SetPeriodicBoundary(geometry, config);
	periodic->SetMeshFile(geometry, config, config->GetMesh_Out_FileName());
	
	/*--- Output of the grid for debuging purposes ---*/
  strcpy (buffer_plt, "periodic_halo.plt"); periodic->SetTecPlot(buffer_plt);
	
#ifndef NO_MPI
	MPI::Finalize();
#endif
    
	/*--- End solver ---*/
	cout << endl <<"------------------------- Exit Success (SU2_PBC) ------------------------" << endl << endl;
	
	return EXIT_SUCCESS;
}
