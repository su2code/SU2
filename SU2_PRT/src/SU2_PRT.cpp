/*!
 * \file SU2_PRT.cpp
 * \brief Main file of Domain Decomposition Code (SU2_PRT).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.2 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

#include "../include/SU2_PRT.hpp"
using namespace std;

int main(int argc, char *argv[]) {
	
	unsigned short nZone = 1;
  double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
	char buffer_su2[8], buffer_plt[8], file_name[MAX_STRING_SIZE];
	string MeshFile;
  
  int rank = MASTER_NODE;
  int size = 1;
  
#ifdef HAVE_MPI
	/*--- MPI initialization ---*/
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
	
	/*--- Definition of some important class ---*/
	CConfig *config = NULL;
	CGeometry *geometry = NULL;
	CSurfaceMovement *surface_mov = NULL;
	CFreeFormDefBox** FFDBox = NULL;
  
	/*--- Definition of the config problem ---*/
	if (argc == 2) { config = new CConfig(argv[1], SU2_PRT, ZONE_0, nZone, 0, VERB_HIGH); }
	else { strcpy (file_name, "default.cfg"); config = new CConfig(file_name, SU2_PRT, ZONE_0, nZone, 0, VERB_HIGH); }
  
  if (rank == MASTER_NODE) {
    
    /*--- Definition of the Class for the geometry ---*/
    geometry = new CPhysicalGeometry(config, ZONE_0, nZone);
    
  }
  
#ifndef HAVE_MPI
  StartTime = double(clock())/double(CLOCKS_PER_SEC);
#else
  MPI_Barrier(MPI_COMM_WORLD);
  StartTime = MPI_Wtime();
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
    
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    /*--- Allocate the memory of the current domain, and
     divide the grid between the nodes ---*/
    CPhysicalGeometry *domain = new CPhysicalGeometry(geometry, config);

    /*--- Add the Send/Receive boundaries ---*/
    domain->SetSendReceive(config);
    
    /*--- Setting the right order for the MPI boundaries ---*/
    domain->SetBoundaries(config);

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if (rank == MASTER_NODE)
      cout << endl <<"----------------------------- Write mesh files --------------------------" << endl;
    
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

    /*--- Write tecplot files ---*/
    if (config->GetVisualize_Partition()) {
      sprintf (buffer_plt, "_%d.dat", int(rank+1));
      string MeshFile_plt = MeshFile + buffer_plt;
      char *cstr_plt = strdup(MeshFile_plt.c_str());
      domain->SetTecPlot(cstr_plt, true);
    }

    /*--- Write .su2 file ---*/
    sprintf (buffer_su2, "_%d.su2", int(rank+1));
    string MeshFile_su2 = MeshFile + buffer_su2;
    char *cstr_su2 = strdup(MeshFile_su2.c_str());
    domain->SetMeshFile(config, cstr_su2);
    
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if (rank == MASTER_NODE)
      cout << "Mesh writing done (" << MeshFile <<")." << endl;
    
    /*--- Write the FFD information ---*/
    
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if (rank == MASTER_NODE)
      cout << endl <<"---------------------- Read and write FFD information -------------------" << endl;
    
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];
    surface_mov = new CSurfaceMovement();
    surface_mov->ReadFFDInfo(domain, config, FFDBox, config->GetMesh_FileName(), false);
    surface_mov->WriteFFDInfo(domain, config, FFDBox, cstr_su2);
    
  }
  
  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/
  
#ifndef HAVE_MPI
  StopTime = double(clock())/double(CLOCKS_PER_SEC);
#else
  MPI_Barrier(MPI_COMM_WORLD);
  StopTime = MPI_Wtime();
#endif
  
  /*--- Compute/print the total time for performance benchmarking. ---*/
  
  UsedTime = StopTime-StartTime;
  if (rank == MASTER_NODE) {
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on "<< size;
    if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
  }
  
  
	/*--- End solver ---*/
  
	if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_PRT) ------------------------" << endl << endl;
	
  
#ifdef HAVE_MPI
  /*--- Finalize MPI parallelization ---*/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
  
	return EXIT_SUCCESS;
	
}
