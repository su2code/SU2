/*!
 * \file SU2_MDC.cpp
 * \brief Main file of Mesh Deformation Code (SU2_MDC).
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

#include "../include/SU2_MDC.hpp"
using namespace std;

int main(int argc, char *argv[]) {
	unsigned short nZone = 1, iZone;
	char buffer_char[50], out_file[200], mesh_file[200];
	int rank = MASTER_NODE, size = 1;
	
#ifndef NO_MPI
	/*--- MPI initialization ---*/
	MPI::Init(argc,argv);
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();
#endif
	
  /*--- Check the number of zones in the mesh file ---*/
	CConfig *config_container = NULL;
	if (argc == 2) config_container = new CConfig(argv[1]);
	else { strcpy (mesh_file, "default.cfg"); config_container = new CConfig(mesh_file); }
	nZone = GetnZone(config_container->GetMesh_FileName(), config_container->GetMesh_FileFormat(), config_container);
  
	/*--- Pointer to different structures that will be used throughout the entire code ---*/
	CConfig **config = NULL;
	CPhysicalGeometry **geometry = NULL;
	CSurfaceMovement *surface_movement = NULL;
  CVolumetricMovement *grid_movement = NULL;

  /*--- Definition of the containers per zones ---*/
  config = new CConfig*[nZone];
	geometry = new CPhysicalGeometry*[nZone];
  
  /*--- Instantiate the config and geometry objects for all zones. We will
   assume that the mesh deformation only occurs in the first zone for now
   but this can be extended in the future. ---*/
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
    geometry[iZone] = new CPhysicalGeometry(config[iZone], config[iZone]->GetMesh_FileName(),
                                            config[iZone]->GetMesh_FileFormat(), iZone+1, nZone);
    
  }
  
  /*--- COmputational grid preprocesing ---*/
	if (rank == MASTER_NODE)
		cout << endl << "----------------------- Preprocessing computations ----------------------" << endl;
  
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
	
	/*--- Create the control volume structures ---*/
	if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
	geometry[ZONE_0]->SetBoundControlVolume(config[ZONE_0], ALLOCATE);
  
  /*--- Output original grid (surface and volumetric) ---*/
	if (config[ZONE_0]->GetVisualize_Deformation()) {
    
    if (rank == MASTER_NODE) cout << "Writting original grid in Tecplot format." << endl;

    if (size > 1) sprintf (buffer_char, "_%d.plt", rank+1);
    else sprintf (buffer_char, ".plt");
    strcpy (out_file, "original_volumetric_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetTecPlot(out_file);
    strcpy (out_file, "original_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundTecPlot(config[ZONE_0], out_file);

		if(config[ZONE_0]->GetOutput_FileFormat() == STL) {
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
	MPI::COMM_WORLD.Barrier();
#endif
	
  
  /*--- Volumetric grid deformation ---*/
	if (rank == MASTER_NODE) cout << endl << "----------------------- Volumetric grid deformation ---------------------" << endl;
  
	/*--- Definition of the Class for grid movement ---*/
	grid_movement = new CVolumetricMovement(geometry[ZONE_0]);
  
	if (config[ZONE_0]->GetDesign_Variable(0) != NO_DEFORMATION) {
		if (rank == MASTER_NODE) cout << "Performing the deformation of the volumetric grid." << endl;
		grid_movement->SetVolume_Deformation(geometry[ZONE_0], config[ZONE_0], false);
  }
	
  /*--- Output deformed grid (surface and volumetric) ---*/
	if (config[ZONE_0]->GetVisualize_Deformation()) {
    
    if (rank == MASTER_NODE) cout << "Writting deformed grid in Tecplot format." << endl;

    if (size > 1) sprintf (buffer_char, "_%d.plt", rank+1);
    else sprintf (buffer_char, ".plt");
    strcpy (out_file, "deformed_volumetric_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetTecPlot(out_file);
    strcpy (out_file, "deformed_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundTecPlot(config[ZONE_0], out_file);
    
		if(config[ZONE_0]->GetOutput_FileFormat() == STL) {
			if (size > 1) sprintf (buffer_char, "_%d.stl", rank+1);
			else sprintf (buffer_char, ".stl");
			strcpy (out_file, "deformed_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundSTL(config[ZONE_0], out_file);
		}
	}
	
  /*--- Check for multiple zones and call appropriate mesh file export routine. ---*/
  if (nZone == 1) SetSingleZone_MeshFile(geometry[ZONE_0], config[ZONE_0], surface_movement);
	else SetMultiZone_MeshFile(geometry, config, nZone);
	
#ifndef NO_MPI
	/*--- Finalize MPI parallelization ---*/
	MPI::Finalize();
#endif
	
	/*--- End solver ---*/
	if (rank == MASTER_NODE)
	  cout << endl << "------------------------- Exit Success (SU2_MDC) ------------------------" << endl << endl;
	
	return EXIT_SUCCESS;
	
}

void SetSingleZone_MeshFile(CPhysicalGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement) {
	char buffer_char[50], out_file[200], in_file[200];
	int rank = MASTER_NODE, size = 1;
	
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();
#endif
  
  if (rank == MASTER_NODE) cout << "Writting .su2 mesh." << endl;
  
  if (size > 1) sprintf (buffer_char, "_%d.su2", rank+1);
  else sprintf (buffer_char, ".su2");
  
  string str = config->GetMesh_Out_FileName();
  str.erase (str.end()-4, str.end()); strcpy (out_file, str.c_str()); strcat(out_file, buffer_char);
  
  str = config->GetMesh_FileName();
  strcpy (in_file, str.c_str());
  
  geometry->SetMeshFile(config, out_file, in_file);
  
  if (geometry->GetnDim() == 3)
		surface_movement->WriteFFDInfo(geometry, config, out_file);
  
}

unsigned short GetnZone(string val_mesh_filename, unsigned short val_format, CConfig *config) {
	string text_line, Marker_Tag;
	ifstream mesh_file;
	short nZone = 1;
	bool isFound = false;
	char cstr[200];
	string::size_type position;
	int rank = MASTER_NODE;
	
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
	if (MPI::COMM_WORLD.Get_size() != 1) {
		val_mesh_filename.erase (val_mesh_filename.end()-4, val_mesh_filename.end());
		val_mesh_filename = val_mesh_filename + "_1.su2";
	}
#endif
	
	/*--- Search the mesh file for the 'NZONE' keyword. ---*/
	switch (val_format) {
		case SU2:
			
			/*--- Open grid file ---*/
			strcpy (cstr, val_mesh_filename.c_str());
			mesh_file.open(cstr, ios::in);
			if (mesh_file.fail()) {
				cout << "There is no geometry file (GetnZone))!" << endl;
				cout << "Press any key to exit..." << endl;
				cin.get();
#ifdef NO_MPI
				exit(1);
#else
				MPI::COMM_WORLD.Abort(1);
				MPI::Finalize();
#endif
			}
			
			/*--- Open the SU2 mesh file ---*/
			while (getline (mesh_file,text_line)) {
				
				/*--- Search for the "NZONE" keyword to see if there are multiple Zones ---*/
				position = text_line.find ("NZONE=",0);
				if (position != string::npos) {
					text_line.erase (0,6); nZone = atoi(text_line.c_str()); isFound = true;
					if (rank == MASTER_NODE) {
						if (nZone <= 0) {
							cout << "Error: Number of mesh zones is less than 1 !!!" << endl;
							cout << "Press any key to exit..." << endl;
							cin.get();
#ifdef NO_MPI
							exit(1);
#else
							MPI::COMM_WORLD.Abort(1);
							MPI::Finalize();
#endif
						}
					}
				}
			}
			/*--- If the "NZONE" keyword was not found, assume this is an ordinary
       simulation on a single Zone ---*/
			if (!isFound) {
				nZone = 1;
			}
			break;
			
		case CGNS:
			nZone = 1;
			break;
			
		case NETCDF_ASCII:
			nZone = 1;
			break;
			
	}
  
  /*--- For time spectral integration, nZones = nTimeInstances. ---*/
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    nZone = 1; //config->GetnTimeInstances();
  }
	
	return (unsigned short) nZone;
}

void SetMultiZone_MeshFile(CPhysicalGeometry **geometry, CConfig **config, unsigned short nZone) {
  
  /*--- Local variables ---*/
	unsigned long iElem, iPoint, iElem_Bound;
	unsigned short iMarker, iNodes, iDim, iZone;
	unsigned short iPeriodic, nPeriodic = 0;
  unsigned short nDim = geometry[ZONE_0]->GetnDim();
	ofstream output_file;
	string Grid_Marker;
	char *cstr;
	double *center, *angles, *transl;
  string val_mesh_out_filename = config[ZONE_0]->GetMesh_Out_FileName();
	cstr = new char [val_mesh_out_filename.size()+1];
	strcpy (cstr, val_mesh_out_filename.c_str());
  
	/*--- Open .su2 grid file ---*/
	output_file.precision(15);
	output_file.open(cstr, ios::out);
  
  /*--- Write the total number of zones first, then loop through each. ---*/
  output_file << "NZONE= " << nZone << endl;
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Write dimension, number of elements and number of points ---*/
    output_file << "IZONE= " << iZone+1 << endl;
    output_file << "NDIME= " << nDim << endl;
    output_file << "NELEM= " << geometry[iZone]->GetnElem() << endl;
    for (iElem = 0; iElem < geometry[iZone]->GetnElem(); iElem++) {
      output_file << geometry[iZone]->elem[iElem]->GetVTK_Type();
      for (iNodes = 0; iNodes < geometry[iZone]->elem[iElem]->GetnNodes(); iNodes++)
        output_file << "\t" << geometry[iZone]->elem[iElem]->GetNode(iNodes);
      output_file << "\t"<<iElem<<endl;
    }
    
    /*--- Write the node coordinates ---*/
    output_file << "NPOIN= " << geometry[iZone]->GetnPoint() << "\t" << geometry[iZone]->GetnPointDomain() << endl;
    output_file.precision(15);
    for (iPoint = 0; iPoint < geometry[iZone]->GetnPoint(); iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++)
        output_file << scientific << "\t" << geometry[iZone]->node[iPoint]->GetCoord(iDim) ;
#ifdef NO_MPI
      output_file << "\t" << iPoint << endl;
#else
      output_file << "\t" << iPoint << "\t" << geometry[iZone]->node[iPoint]->GetGlobalIndex() << endl;
#endif
    }
    
    /*--- Loop through and write the boundary info ---*/
    output_file << "NMARK= " << config[iZone]->GetnMarker_All() << endl;
    for (iMarker = 0; iMarker < config[iZone]->GetnMarker_All(); iMarker++) {
      
      /*--- Ignore SEND_RECEIVE for the moment ---*/
      if (geometry[iZone]->bound[iMarker][0]->GetVTK_Type() != VERTEX) {
        
        Grid_Marker = config[iZone]->GetMarker_All_Tag(iMarker);
        output_file << "MARKER_TAG= " << Grid_Marker <<endl;
        output_file << "MARKER_ELEMS= " << geometry[iZone]->nElem_Bound[iMarker]<< endl;
        
        if (nDim == 2) {
          for (iElem_Bound = 0; iElem_Bound < geometry[iZone]->nElem_Bound[iMarker]; iElem_Bound++) {
            output_file << geometry[iZone]->bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
            for (iNodes = 0; iNodes < geometry[iZone]->bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
              output_file << geometry[iZone]->bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
            output_file	<< iElem_Bound << endl;
          }
        }
        
        if (nDim == 3) {
          for (iElem_Bound = 0; iElem_Bound < geometry[iZone]->nElem_Bound[iMarker]; iElem_Bound++) {
            output_file << geometry[iZone]->bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
            for (iNodes = 0; iNodes < geometry[iZone]->bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
              output_file << geometry[iZone]->bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
            output_file	<< iElem_Bound << endl;
          }
        }
        
      } else if (geometry[iZone]->bound[iMarker][0]->GetVTK_Type() == VERTEX) {
        output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
        output_file << "MARKER_ELEMS= " << geometry[iZone]->nElem_Bound[iMarker]<< endl;
        if (config[iZone]->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config[iZone]->GetMarker_All_SendRecv(iMarker) << endl;
        if (config[iZone]->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config[iZone]->GetMarker_All_SendRecv(iMarker) << endl;
        
        for (iElem_Bound = 0; iElem_Bound < geometry[iZone]->nElem_Bound[iMarker]; iElem_Bound++) {
          output_file << geometry[iZone]->bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" <<
          geometry[iZone]->bound[iMarker][iElem_Bound]->GetNode(0) << "\t" <<
          geometry[iZone]->bound[iMarker][iElem_Bound]->GetRotation_Type() << "\t" <<
          geometry[iZone]->bound[iMarker][iElem_Bound]->GetMatching_Zone()<< endl;
        }
        
      }
    }
    
    /*--- Get the total number of periodic transformations ---*/
    nPeriodic = config[iZone]->GetnPeriodicIndex();
    output_file << "NPERIODIC= " << nPeriodic << endl;
    
    /*--- From iPeriodic obtain the iMarker ---*/
    for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
      
      /*--- Retrieve the supplied periodic information. ---*/
      center = config[iZone]->GetPeriodicCenter(iPeriodic);
      angles = config[iZone]->GetPeriodicRotation(iPeriodic);
      transl = config[iZone]->GetPeriodicTranslate(iPeriodic);
      
      output_file << "PERIODIC_INDEX= " << iPeriodic << endl;
      output_file << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
      output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << endl;
      output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << endl;
      
    }
  }
	output_file.close();
  
}
