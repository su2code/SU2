/*!
 * \file SU2_DDC.cpp
 * \brief Main file of Domain Decomposition Code (SU2_DDC).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.10
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/SU2_DDC.hpp"
using namespace std;

int main(int argc, char *argv[]) {
	
	unsigned short nZone = 1;
	char buffer_su2[8], buffer_plt[8], file_name[200];
	string MeshFile;
  
  int rank = MASTER_NODE;
  int size = 1;
  
#ifndef NO_MPI
	/*--- MPI initialization, and buffer setting ---*/
  static char buffer[MAX_MPI_BUFFER]; // buffer size in bytes
  
  void *ptr;

#ifdef WINDOWS
	MPI_Init(&argc,&argv);
	MPI_Buffer_attach(buffer,MAX_MPI_BUFFER);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
#else
	MPI::Init(argc, argv);
	MPI::Attach_buffer(buffer, MAX_MPI_BUFFER);
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();
#endif	
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
    geometry = new CPhysicalGeometry(config, ZONE_0, nZone);
    
  }
  
#ifndef NO_MPI
#ifdef WINDOWS
  MPI_Barrier(MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Barrier();
#endif
#endif
  
	/*--- Set domains for parallel computation (if any) ---*/
	if (size > 1) {
		
    /*--- Write the new subgrid, and remove the extension ---*/
    MeshFile = config->GetMesh_FileName();
    unsigned short lastindex = MeshFile.find_last_of(".");
    MeshFile = MeshFile.substr(0, lastindex);

    if (rank == MASTER_NODE) {
      
      cout << endl <<"------------------------ Divide the numerical grid ----------------------" << endl;
      
      /*--- Color the initial grid and set the send-receive domains ---*/
      geometry->SetColorGrid(config);
      
    }
    
#ifndef NO_MPI
#ifdef WINDOWS
  MPI_Barrier(MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Barrier();
#endif
#endif
    
    /*--- Allocate the memory of the current domain, and
     divide the grid between the nodes ---*/
    CDomainGeometry *domain = new CDomainGeometry(geometry, config);
    
    /*--- Add the Send/Receive boundaries ---*/
    domain->SetSendReceive(config);
    
#ifndef NO_MPI
#ifdef WINDOWS
  MPI_Barrier(MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Barrier();
#endif
#endif
    
    if (rank == MASTER_NODE)
      cout << endl <<"----------------------------- Write mesh files --------------------------" << endl;
    
#ifndef NO_MPI
#ifdef WINDOWS
  MPI_Barrier(MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Barrier();
#endif
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
#ifdef WINDOWS
  MPI_Barrier(MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Barrier();
#endif
#endif
    
    if (rank == MASTER_NODE)
      cout << "Mesh writing done (" << MeshFile <<")." << endl;
    
    /*--- Write the FFD information (3D problems)---*/
    if (domain->GetnDim() == 3) {
      
#ifndef NO_MPI
#ifdef WINDOWS
  MPI_Barrier(MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Barrier();
#endif
#endif
      
      if (rank == MASTER_NODE)
        cout << endl <<"---------------------- Read and write FFD information -------------------" << endl;
      
#ifndef NO_MPI
#ifdef WINDOWS
  MPI_Barrier(MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Barrier();
#endif
#endif
      
      FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];
      surface_mov = new CSurfaceMovement();
      surface_mov->ReadFFDInfo(domain, config, FFDBox, config->GetMesh_FileName(), false);
      surface_mov->WriteFFDInfo(domain, config, FFDBox, cstr_su2);
      
    }
    
#ifndef NO_MPI
    /*--- Finalize MPI parallelization ---*/
#ifdef WINDOWS
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Buffer_detach(buffer,NULL);
	MPI_Finalize();
#else
    MPI::COMM_WORLD.Barrier();
    MPI::Detach_buffer(ptr);
    MPI::Finalize();
#endif
#endif
    
  }
  
	/*--- End solver ---*/
	if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_DDC) ------------------------" << endl << endl;
	
	return EXIT_SUCCESS;
	
}
