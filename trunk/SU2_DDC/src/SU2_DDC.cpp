/*!
 * \file SU2_DDC.cpp
 * \brief Main file of Domain Decomposition Code (SU2_DDC).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.1
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
	
	unsigned short iDomain, nZone = 1;
	char buffer_su2[8], buffer_vtk[8], buffer_plt[8], file_name[200];
	string MeshFile;
	
#ifndef NO_MPI
	MPI::Init(argc, argv);
#endif
	
	/*--- Definition of some important class ---*/
	CConfig *config = NULL;
	CSurfaceMovement *surface_mov = NULL;
	CFreeFormChunk** chunk = NULL;
	unsigned short nChunk = MAX_NUMBER_CHUNK;

	/*--- Definition of the config problem ---*/
	if (argc == 2) config = new CConfig(argv[1], SU2_DDC, ZONE_0, nZone, VERB_HIGH);
	else { strcpy (file_name, "default.cfg"); config = new CConfig(file_name, SU2_DDC, ZONE_0, nZone, VERB_HIGH); }
	
	/*--- Definition of the Class for the geometry ---*/
	CGeometry *geometry; geometry = new CGeometry;
  geometry = new CPhysicalGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat(), ZONE_0, nZone);
  		
	cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
		
	/*--- Create the edge structure ---*/
	cout << "Identifying vertices." <<endl; 
	geometry->SetVertex(config);
		
	/*--- Read the FFD information from the grid file (3D problems)---*/
	if (geometry->GetnDim() == 3) {
		cout << endl <<"--------------------------- Read FFD information ------------------------" << endl;
		geometry->SetVertex(config);
		chunk = new CFreeFormChunk*[nChunk];
		surface_mov = new CSurfaceMovement();
		surface_mov->ReadFFDInfo(config, geometry, chunk, config->GetMesh_FileName());
	}
	
	/*--- Set domains for parallel computation (if any) ---*/
	if (config->GetnDomain() > 1) {
		
		/*--- Color the initial grid and set the send-receive domains ---*/
		geometry->SetColorGrid(config);
		geometry->SetSendReceive(config);
				
		/*--- Write the new subgrid ---*/
		CGeometry **domain; domain = new CGeometry *[config->GetnDomain()];
		MeshFile = config->GetMesh_FileName();
		MeshFile.erase (MeshFile.end()-4, MeshFile.end());
		
		cout << endl <<"----------------------------- Write mesh files --------------------------" << endl;

		for (iDomain = 0; iDomain < config->GetnDomain(); iDomain++) {
			
			domain[iDomain] = new CDomainGeometry(geometry, config, iDomain);
			domain[iDomain]->SetSendReceive(geometry, config, iDomain);

			/*--- Write tecplot and paraview files ---*/
			if (config->GetVisualize_Partition()) {
				if(config->GetOutput_FileFormat() == PARAVIEW) {
					sprintf (buffer_vtk, "_%d.vtk", int(iDomain+1));
					string MeshFile_vtk = MeshFile + buffer_vtk;
					char *cstr_vtk = strdup(MeshFile_vtk.c_str());
					domain[iDomain]->SetParaView(cstr_vtk);
				}
				if(config->GetOutput_FileFormat() == TECPLOT) {
					sprintf (buffer_plt, "_%d.plt", int(iDomain+1));
					string MeshFile_plt = MeshFile + buffer_plt;
					char *cstr_plt = strdup(MeshFile_plt.c_str());
					domain[iDomain]->SetTecPlot(cstr_plt);
				}
			}
			
			/*--- Write .su2 file ---*/
			sprintf (buffer_su2, "_%d.su2", int(iDomain+1));
			string MeshFile_su2 = MeshFile + buffer_su2;
			char *cstr_su2 = strdup(MeshFile_su2.c_str());
			domain[iDomain]->SetMeshFile(geometry, config, cstr_su2, iDomain);
			
			/*--- Write the FFD information (3D problems)---*/
			if (geometry->GetnDim() == 3)
				surface_mov->WriteFFDInfo(geometry, domain[iDomain], config, chunk, cstr_su2);
			
		}		
	}
	
#ifndef NO_MPI
	MPI::Finalize();
#endif
	
	/*--- End solver ---*/
	cout << endl <<"------------------------- Exit Success (SU2_DDC) ------------------------" << endl << endl;
	
	return EXIT_SUCCESS;
	
}
