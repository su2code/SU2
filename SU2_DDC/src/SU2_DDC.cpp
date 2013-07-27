/*!
 * \file SU2_DDC.cpp
 * \brief Main file of Domain Decomposition Code (SU2_DDC).
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

#include "../include/SU2_DDC.hpp"
using namespace std;

int main(int argc, char *argv[]) {
	
	unsigned short nZone = 1;
	char buffer_su2[8], buffer_vtk[8], buffer_plt[8], file_name[200];
	string MeshFile;
  
  int rank = MASTER_NODE;
  int size = 1;
  
#ifndef NO_MPI
	/*--- MPI initialization, and buffer setting ---*/
  static char buffer[MAX_MPI_BUFFER]; // buffer size in bytes
  
  void *ptr;
	MPI::Init(argc, argv);
	MPI::Attach_buffer(buffer, MAX_MPI_BUFFER);
  
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();
#endif
	
	/*--- Definition of some important class ---*/
	CConfig *config = NULL;
  CGeometry *geometry = NULL;
	CSurfaceMovement *surface_mov = NULL;
	CFreeFormDefBox** FFDBox = NULL;
  
	/*--- Definition of the config problem ---*/
	if (argc == 2) { config = new CConfig(argv[1], SU2_DDC, ZONE_0, nZone, VERB_HIGH); }
	else { strcpy (file_name, "default.cfg"); config = new CConfig(file_name, SU2_DDC, ZONE_0, nZone, VERB_HIGH); }
  
  if (rank == MASTER_NODE) {
    
    /*--- Definition of the Class for the geometry ---*/
    geometry = new CPhysicalGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat(), ZONE_0, nZone);
    
  }
  
#ifndef NO_MPI
  MPI::COMM_WORLD.Barrier();
#endif
  
	/*--- Set domains for parallel computation (if any) ---*/
	if (size > 1) {
		
    /*--- Write the new subgrid ---*/
    MeshFile = config->GetMesh_FileName();
    MeshFile.erase (MeshFile.end()-4, MeshFile.end());
    
    if (rank == MASTER_NODE) {
      
      cout << endl <<"------------------------ Divide the numerical grid ----------------------" << endl;
      
      /*--- Color the initial grid and set the send-receive domains ---*/
      geometry->SetColorGrid(config);
      
    }
    
#ifndef NO_MPI
    MPI::COMM_WORLD.Barrier();
#endif
    
    /*--- Allocate the memory of the current domain, and
     divide the grid between the nodes ---*/
    CDomainGeometry *domain = new CDomainGeometry(geometry, config);
    
    /*--- Add the Send/Receive boundaries ---*/
    domain->SetSendReceive(config);
    
#ifndef NO_MPI
    MPI::COMM_WORLD.Barrier();
#endif
    
    if (rank == MASTER_NODE)
      cout << endl <<"----------------------------- Write mesh files --------------------------" << endl;
    
#ifndef NO_MPI
    MPI::COMM_WORLD.Barrier();
#endif
    
    /*--- Write tecplot files ---*/
    if (config->GetVisualize_Partition()) {
      sprintf (buffer_plt, "_%d.dat", int(rank+1));
      string MeshFile_plt = MeshFile + buffer_plt;
      char *cstr_plt = strdup(MeshFile_plt.c_str());
      domain->SetTecPlot(cstr_plt);
    }
    
    /*--- Write .su2 file ---*/
    sprintf (buffer_su2, "_%d.su2", int(rank+1));
    string MeshFile_su2 = MeshFile + buffer_su2;
    char *cstr_su2 = strdup(MeshFile_su2.c_str());
    domain->SetMeshFile(config, cstr_su2);
    
#ifndef NO_MPI
    MPI::COMM_WORLD.Barrier();
#endif
    
    cout << "Domain " << rank <<": Mesh writing done (" << MeshFile_su2 <<")." << endl;
    
#ifndef NO_MPI
    MPI::COMM_WORLD.Barrier();
#endif
    
    /*--- Write the FFD information (3D problems)---*/
    if (domain->GetnDim() == 3) {
      
#ifndef NO_MPI
      MPI::COMM_WORLD.Barrier();
#endif
      
      if (rank == MASTER_NODE)
        cout << endl <<"---------------------- Read and write FFD information -------------------" << endl;
      
#ifndef NO_MPI
      MPI::COMM_WORLD.Barrier();
#endif
      
      FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];
      surface_mov = new CSurfaceMovement();
      surface_mov->ReadFFDInfo(domain, config, FFDBox, config->GetMesh_FileName(), false);
      surface_mov->WriteFFDInfo(domain, config, FFDBox, cstr_su2);
      
    }
    
#ifndef NO_MPI
    /*--- Finalize MPI parallelization ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::Detach_buffer(ptr);
    MPI::Finalize();
#endif
    
  }
  
	/*--- End solver ---*/
	if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_DDC) ------------------------" << endl << endl;
	
	return EXIT_SUCCESS;
	
}
