/*!
 * \file SU2_DDC.cpp
 * \brief Main file of Domain Decomposition Code (SU2_DDC).
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

#include "../include/SU2_DDC.hpp"

using namespace std;

int main(int argc, char *argv[]) {
	
	unsigned short iDomain;
	char buffer_su2[8], buffer_vtk[8]; 
	string MeshFile;
	
	/*--- Definition of the Class for the definition of the problem ---*/
	CConfig *config;
	if (argc == 2) config = new CConfig(argv[1], SU2_DDC);
	else {
		char grid_file[200];
		strcpy (grid_file, "default.cfg");
		config = new CConfig(grid_file, SU2_DDC);
	}
	
	/*--- Definition of the Class for the geometry ---*/
	CGeometry *geometry; geometry = new CGeometry;
	
	geometry = new CPhysicalGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat());

	cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
	
	/*--- Compute element surrounding points, point surrounding points, and element surrounding elements ---*/
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
	
	/*--- Set domains for parallel computation (if any) ---*/
	if (config->GetnDomain() > 1) {
		
		/*--- Color the initial grid and set the send-receive domains ---*/
		geometry->SetColorGrid(config, config->GetnDomain());
		geometry->SetSendReceive(config, config->GetnDomain());
		
		/*--- Write .su2 and .vtk format ---*/
		CGeometry **domain; domain = new CGeometry *[config->GetnDomain()];
		MeshFile = config->GetMesh_FileName();
		MeshFile.erase (MeshFile.end()-4, MeshFile.end());
		
		for (iDomain = 0; iDomain < config->GetnDomain(); iDomain++) {
			domain[iDomain] = new CDomainGeometry(geometry, config, iDomain);
			domain[iDomain]->SetSendReceive(geometry, config, iDomain);
			
			if (config->GetVisualize_Partition()) {
				sprintf (buffer_vtk, "_%d.vtk", int(iDomain+1));
				string MeshFile_vtk = MeshFile + buffer_vtk;
				char *cstr_vtk = strdup(MeshFile_vtk.c_str());
				domain[iDomain]->SetParaView(cstr_vtk);
			}
			
			sprintf (buffer_su2, "_%d.su2", int(iDomain+1));
			string MeshFile_su2 = MeshFile + buffer_su2;
			char *cstr_su2 = strdup(MeshFile_su2.c_str());
			domain[iDomain]->SetMeshFile(config, cstr_su2);
			
		}		
		
	}
	
	return 1;
}
