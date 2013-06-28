/*!
 * \file SU2_PBC.cpp
 * \brief Main file of Periodic Boundary Code (SU2_PBC).
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

#include "../include/SU2_PBC.hpp"

using namespace std;

int main(int argc, char *argv[]) {
	
	char buffer_vtk[100]; 
	string MeshFile;
	
	/*--- Definition of the class for the definition of the problem ---*/
	CConfig *config;
	if (argc == 2) config = new CConfig(argv[1], SU2_PBC);
	else {
		char grid_file[200];
		strcpy (grid_file, "default.cfg");
		config = new CConfig(grid_file, SU2_PBC);
	}
	
	/*--- Definition of the class for the geometry ---*/
	CGeometry *geometry; geometry = new CGeometry;
	geometry = new CPhysicalGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat());
  
	cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
	
	/*--- Compute elements surrounding points, points surrounding points, and elements surronding elements ---*/
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
	
	/*--- Set periodic boundary conditions ---*/
	geometry->SetPeriodicBoundary(config);
	
  /*--- Original grid in .vtk format for debugigng purposes ---*/
  strcpy (buffer_vtk, "periodic_orig.vtk"); geometry->SetParaView(buffer_vtk);	
  
	/*--- Create a new grid with the right periodic boundary ---*/
	CGeometry *periodic; periodic = new CPeriodicGeometry(geometry, config);
	periodic->SetPeriodicBoundary(geometry, config);
	periodic->SetMeshFile(config, config->GetMesh_Out_FileName());
	
	/*--- Output of the grid in .vtk format for debuging purposes ---*/
	strcpy (buffer_vtk, "periodic.vtk"); periodic->SetParaView(buffer_vtk);	
	
	return 1;
}
