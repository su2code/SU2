/*!
 * \file SU2_GDC.cpp
 * \brief Main file of Geometry Design Code (SU2_GDC).
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

#include "../include/SU2_GDC.hpp"

using namespace std;

int main(int argc, char *argv[]) {
	
	/*--- Variable definitions ---*/
	char file_name[200];
	
	/*-- Definition of the Class for the definition of the problem ---*/
	CConfig *config;
	if (argc == 2) config = new CConfig(argv[1], SU2_GDC);
	else {
		strcpy (file_name, "default.cfg");
		config = new CConfig(file_name, SU2_GDC);
	}
	
	/*-- Definition of the Class for the geometry ---*/
	CGeometry *geometry; 
	geometry = new CGeometry;
	geometry = new CPhysicalGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat());
	
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
	geometry->SetControlVolume(config, ALLOCATE);
	geometry->SetBoundControlVolume(config, ALLOCATE);
	
	return 1;
}
