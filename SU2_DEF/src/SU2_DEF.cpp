/*!
 * \file SU2_DEF.cpp
 * \brief Main file of Mesh Deformation Code (SU2_DEF).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.0 "eagle"
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

#include "../include/SU2_DEF.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  
  double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  unsigned short nZone = 1;
  char buffer_char[50], out_file[MAX_STRING_SIZE], in_file[MAX_STRING_SIZE], mesh_file[MAX_STRING_SIZE];
  int rank = MASTER_NODE, size = SINGLE_NODE;
  string str;
  
#ifdef HAVE_MPI
  /*--- MPI initialization ---*/
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
#endif
  
  /*--- Pointer to different structures that will be used throughout the entire code ---*/
  
  CConfig **config = NULL;
  CPhysicalGeometry **geometry = NULL;
  CSurfaceMovement *surface_movement = NULL;
  CVolumetricMovement *grid_movement = NULL;
  
  /*--- Definition of the containers by zone (currently only one zone is
   allowed, but this can be extended if necessary). ---*/
  
  config   = new CConfig*[1];
  geometry = new CPhysicalGeometry*[1];
  
  /*--- Definition of the configuration class, and open the config file ---*/
  
  if (argc == 2) config[ZONE_0] = new CConfig(argv[1], SU2_DEF, ZONE_0, nZone, 0, VERB_HIGH);
  else {
    strcpy (mesh_file, "default.cfg");
    config[ZONE_0] = new CConfig(mesh_file, SU2_DEF, ZONE_0, nZone, 0, VERB_HIGH);
  }
  
#ifdef HAVE_MPI
  
  /*--- Change the name of the input-output files for the parallel computation ---*/
  
  config[ZONE_0]->SetFileNameDomain(rank+1);
  
#endif
  
  /*--- Definition of the geometry class ---*/
  
  geometry[ZONE_0] = new CPhysicalGeometry(config[ZONE_0], ZONE_0, nZone);
  
  /*--- Set up a timer for performance benchmarking (preprocessing time is not included) ---*/
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  StartTime = MPI_Wtime();
#else
  StartTime = double(clock())/double(CLOCKS_PER_SEC);
#endif
  
  /*--- Computational grid preprocesing ---*/
  
  if (rank == MASTER_NODE) cout << endl << "----------------------- Preprocessing computations ----------------------" << endl;
  
  /*--- Compute elements surrounding points, points surrounding points ---*/
  
  if (rank == MASTER_NODE) cout << "Setting local point connectivity." <<endl;
  geometry[ZONE_0]->SetPoint_Connectivity();
  
  /*--- Check the orientation before computing geometrical quantities ---*/
  
  if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." <<endl;
  geometry[ZONE_0]->SetBoundVolume();
  geometry[ZONE_0]->Check_Orientation(config[ZONE_0]);
  
  /*--- Create the edge structure ---*/
  
  if (rank == MASTER_NODE) cout << "Identify edges and vertices." <<endl;
  geometry[ZONE_0]->SetEdges(); geometry[ZONE_0]->SetVertex(config[ZONE_0]);
  
  /*--- Compute center of gravity ---*/
  
  if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
  geometry[ZONE_0]->SetCG();
  
  /*--- Create the dual control volume structures ---*/
  
  if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
  geometry[ZONE_0]->SetBoundControlVolume(config[ZONE_0], ALLOCATE);
  
  /*--- Output original grid for visualization, if requested (surface and volumetric) ---*/
  
  if (config[ZONE_0]->GetVisualize_Deformation()) {
    
    if (rank == MASTER_NODE) cout << "Writing an Tecplot file of the volumetric mesh." << endl;
    if (size > 1) sprintf (buffer_char, "_%d.plt", rank+1); else sprintf (buffer_char, ".plt");
    strcpy (out_file, "Volumetric_Grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetTecPlot(out_file, true);
    
    if (rank == MASTER_NODE) cout << "Writing an Tecplot file of the surface mesh." << endl;
    if (size > 1) sprintf (buffer_char, "_%d.plt", rank+1); else sprintf (buffer_char, ".plt");
    strcpy (out_file, "Surface_Grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundTecPlot(out_file, true, config[ZONE_0]);
    
    if (rank == MASTER_NODE) cout << "Writing an STL file of the surface mesh." << endl;
    if (size > 1) sprintf (buffer_char, "_%d.stl", rank+1); else sprintf (buffer_char, ".stl");
    strcpy (out_file, "Surface_Grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundSTL(out_file, true, config[ZONE_0]);
    
  }
  
  /*--- Surface grid deformation using design variables ---*/
  
  if (rank == MASTER_NODE) cout << endl << "------------------------- Surface grid deformation ----------------------" << endl;
  
  /*--- Definition and initialization of the surface deformation class ---*/
  
  surface_movement = new CSurfaceMovement();
  surface_movement->CopyBoundary(geometry[ZONE_0], config[ZONE_0]);
  
  /*--- Surface grid deformation ---*/
  
  if (rank == MASTER_NODE) cout << "Performing the deformation of the surface grid." << endl;
  surface_movement->SetSurface_Deformation(geometry[ZONE_0], config[ZONE_0]);
  
#ifdef HAVE_MPI
  /*--- MPI syncronization point ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- Volumetric grid deformation ---*/
  
  if (config[ZONE_0]->GetDesign_Variable(0) != FFD_SETTING) {
    
    if (rank == MASTER_NODE) cout << endl << "----------------------- Volumetric grid deformation ---------------------" << endl;
    
    /*--- Definition of the Class for grid movement ---*/
    
    grid_movement = new CVolumetricMovement(geometry[ZONE_0]);
    
    if (rank == MASTER_NODE) cout << "Performing the deformation of the volumetric grid." << endl;
    
    grid_movement->SetVolume_Deformation(geometry[ZONE_0], config[ZONE_0], false);
    
  }
  
  /*--- Computational grid preprocesing ---*/
  
  if (rank == MASTER_NODE) cout << endl << "----------------------- Write deformed grid files -----------------------" << endl;
  
  /*--- Output deformed grid for visualization, if requested (surface and volumetric) ---*/
  
  if (config[ZONE_0]->GetVisualize_Deformation()) {
    
    if (rank == MASTER_NODE) cout << "Writing a Tecplot file of the volumetric mesh." << endl;
    if (size > 1) sprintf (buffer_char, "_%d.plt", rank+1); else sprintf (buffer_char, ".plt");
    strcpy (out_file, "Volumetric_Grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetTecPlot(out_file, false);
    
    if (rank == MASTER_NODE) cout << "Writing a Tecplot file of the surface mesh." << endl;
    if (size > 1) sprintf (buffer_char, "_%d.plt", rank+1); else sprintf (buffer_char, ".plt");
    strcpy (out_file, "Surface_Grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundTecPlot(out_file, false, config[ZONE_0]);

    if (rank == MASTER_NODE) cout << "Writing a STL file of the surface mesh." << endl;
    if (size > 1) sprintf (buffer_char, "_%d.stl", rank+1); else sprintf (buffer_char, ".stl");
    strcpy (out_file, "Surface_Grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundSTL(out_file, false, config[ZONE_0] );
    
  }
  
  /*--- Write the new SU2 native mesh after deformation. ---*/
  
  if (rank == MASTER_NODE) cout << "Writing a SU2 file of the volumetric mesh." << endl;
  
  if (size > 1) sprintf (buffer_char, "_%d.su2", rank+1); else sprintf (buffer_char, ".su2");

  str = config[ZONE_0]->GetMesh_Out_FileName(); str.erase (str.end()-4, str.end());
  strcpy (out_file, str.c_str()); strcat(out_file, buffer_char);
  
  str = config[ZONE_0]->GetMesh_FileName();
  strcpy (in_file, str.c_str());
  
  geometry[ZONE_0]->SetMeshFile(config[ZONE_0], out_file, in_file);
  
  /*--- Write the the free-form deformation boxes after deformation. ---*/
  
  if (rank == MASTER_NODE) cout << "Adding FFD information to the SU2 file." << endl;

  surface_movement->WriteFFDInfo(geometry[ZONE_0], config[ZONE_0], out_file);
  
  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  StopTime = MPI_Wtime();
#else
  StopTime = double(clock())/double(CLOCKS_PER_SEC);
#endif
  
  /*--- Compute/print the total time for performance benchmarking. ---*/
  
  UsedTime = StopTime-StartTime;
  if (rank == MASTER_NODE) {
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on "<< size;
    if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
  }
  
  /*--- Exit the solver cleanly ---*/
  
  if (rank == MASTER_NODE)
  cout << endl << "------------------------- Exit Success (SU2_DEF) ------------------------" << endl << endl;
  
#ifdef HAVE_MPI
  /*--- Finalize MPI parallelization ---*/
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  
  return EXIT_SUCCESS;
  
}
