/*!
 * \file SU2_MDC.cpp
 * \brief Main file of Mesh Deformation Code (SU2_MDC).
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

#include "../include/SU2_MDC.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  unsigned short nZone = 1, iZone;
  char buffer_char[50], out_file[200], in_file[200], mesh_file[200];
  int rank = MASTER_NODE, size = SINGLE_NODE;
  
#ifndef NO_MPI
  /*--- MPI initialization ---*/
#ifdef WINDOWS
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
#else  
  MPI::Init(argc,argv);
  rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
#endif
#endif
  
  /*--- Pointer to different structures that will be used throughout the entire code ---*/
  
  CConfig **config = NULL;
  CPhysicalGeometry **geometry = NULL;
  CSurfaceMovement *surface_movement = NULL;
  CVolumetricMovement *grid_movement = NULL;
  
  /*--- Definition of the containers by zone (currently only one zone is
   allowed, but this can be extended if necessary). ---*/
  
  config   = new CConfig*[nZone];
  geometry = new CPhysicalGeometry*[nZone];
  
  /*--- Instantiate the config and geometry objects. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Definition of the configuration class, and open the config file ---*/
    
    if (argc == 2) config[iZone] = new CConfig(argv[1], SU2_MDC, iZone, nZone, VERB_HIGH);
    else {
      strcpy (mesh_file, "default.cfg");
      config[iZone] = new CConfig(mesh_file, SU2_MDC, iZone, nZone, VERB_HIGH);
    }
    
#ifndef NO_MPI
    /*--- Change the name of the input-output files for the parallel computation ---*/
    
    config[iZone]->SetFileNameDomain(rank+1);
#endif
    
    /*--- Definition of the geometry class ---*/
    
    geometry[iZone] = new CPhysicalGeometry(config[iZone], iZone+1, nZone);
    
  }
  
  /*--- Computational grid preprocesing ---*/
  
  if (rank == MASTER_NODE) cout << endl << "----------------------- Preprocessing computations ----------------------" << endl;
  
  /*--- Compute elements surrounding points, points surrounding points ---*/
  
  if (rank == MASTER_NODE) cout << "Setting local point and element connectivity." <<endl;
  geometry[ZONE_0]->SetEsuP();
  geometry[ZONE_0]->SetPsuP();
  
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
    
    if (rank == MASTER_NODE) cout << "Writing original grid in Tecplot format." << endl;
    if (size > 1) sprintf (buffer_char, "_%d.plt", rank+1);
    else sprintf (buffer_char, ".plt");
    strcpy (out_file, "original_volumetric_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetTecPlot(out_file);
    strcpy (out_file, "original_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundTecPlot(config[ZONE_0], out_file);
    
    if(config[ZONE_0]->GetOutput_FileFormat() == STL) {
      if (rank == MASTER_NODE) cout << "Writing an STL file of the original surface mesh." << endl;
      if (size > 1) sprintf (buffer_char, "_%d.stl", rank+1);
      else sprintf (buffer_char, ".stl");
      strcpy (out_file, "original_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundSTL(config[ZONE_0],out_file);
    }
    
  }
  
  /*--- Surface grid deformation using design variables ---*/
  
  if (rank == MASTER_NODE) cout << endl << "------------------------- Surface grid deformation ----------------------" << endl;
  
  /*--- Definition and initialization of the surface deformation class ---*/
  
  surface_movement = new CSurfaceMovement();
  surface_movement->CopyBoundary(geometry[ZONE_0], config[ZONE_0]);
  
  /*--- Surface grid deformation ---*/
  
  if (rank == MASTER_NODE) cout << "Performing the deformation of the surface grid." << endl;
  surface_movement->SetSurface_Deformation(geometry[ZONE_0], config[ZONE_0]);
  
#ifndef NO_MPI
  /*--- MPI syncronization point ---*/
#ifdef WINDOWS
  MPI_Barrier(MPI_COMM_WORLD);
#else 
  MPI::COMM_WORLD.Barrier();
#endif
#endif
  
  /*--- Volumetric grid deformation ---*/
  
  if (rank == MASTER_NODE) cout << endl << "----------------------- Volumetric grid deformation ---------------------" << endl;

  /*--- Definition of the Class for grid movement ---*/
  
  grid_movement = new CVolumetricMovement(geometry[ZONE_0]);
  
  if (config[ZONE_0]->GetDesign_Variable(0) != NO_DEFORMATION) {
    if (rank == MASTER_NODE) cout << "Performing the deformation of the volumetric grid." << endl;
    grid_movement->SetVolume_Deformation(geometry[ZONE_0], config[ZONE_0], false);
  }
  
  /*--- Computational grid preprocesing ---*/
  
  if (rank == MASTER_NODE) cout << endl << "----------------------- Write deformed grid files -----------------------" << endl;
  
  /*--- Output deformed grid for visualization, if requested (surface and volumetric) ---*/
  
  if (config[ZONE_0]->GetVisualize_Deformation()) {
    
    if (rank == MASTER_NODE) cout << "Writing deformed grid in Tecplot format." << endl;
    if (size > 1) sprintf (buffer_char, "_%d.plt", rank+1);
    else sprintf (buffer_char, ".plt");
    strcpy (out_file, "deformed_volumetric_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetTecPlot(out_file);
    strcpy (out_file, "deformed_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundTecPlot(config[ZONE_0], out_file);
    
    if(config[ZONE_0]->GetOutput_FileFormat() == STL) {
      if (rank == MASTER_NODE) cout << "Writing an STL file of the deformed surface mesh." << endl;
      if (size > 1) sprintf (buffer_char, "_%d.stl", rank+1);
      else sprintf (buffer_char, ".stl");
      strcpy (out_file, "deformed_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundSTL(config[ZONE_0], out_file);
    }
  }
  
  /*--- Write the new SU2 native mesh after deformation. ---*/
  
  if (rank == MASTER_NODE) cout << "Writing the new .su2 mesh after deformation." << endl;
  
  if (size > 1) sprintf (buffer_char, "_%d.su2", rank+1);
  else sprintf (buffer_char, ".su2");
  
  string str = config[ZONE_0]->GetMesh_Out_FileName();
  str.erase (str.end()-4, str.end()); strcpy (out_file, str.c_str()); strcat(out_file, buffer_char);
  
  str = config[ZONE_0]->GetMesh_FileName();
  strcpy (in_file, str.c_str());
  
  geometry[ZONE_0]->SetMeshFile(config[ZONE_0], out_file, in_file);
  
  if (geometry[ZONE_0]->GetnDim() == 3)
  surface_movement->WriteFFDInfo(geometry[ZONE_0], config[ZONE_0], out_file);
  
  
#ifndef NO_MPI
  /*--- Finalize MPI parallelization ---*/
#ifdef WINDOWS
  MPI_Finalize();
#else  
  MPI::Finalize();
#endif
#endif
  
  /*--- End solver ---*/
  
  if (rank == MASTER_NODE)
  cout << endl << "------------------------- Exit Success (SU2_MDC) ------------------------" << endl << endl;
  
  return EXIT_SUCCESS;
  
}
