/*!
 * \file SU2_MDC.cpp
 * \brief Main file of Mesh Deformation Code (SU2_MDC).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.3
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
	unsigned short iChunk, nZone = 1, iDV, iLevel, iChild, iParent, jChunk, iZone;
	char buffer_char[50], out_file[200], grid_file[200];
	int rank = MASTER_NODE, size = 1, iExtIter = 0;
	string ChunkTag;
	
#ifndef NO_MPI
	/*--- MPI initialization, and buffer setting ---*/
	void *buffer, *old_buffer;
	int bufsize;
	bufsize = MAX_MPI_BUFFER;
	buffer = new char[bufsize];
	MPI::Init(argc,argv);
	MPI::Attach_buffer(buffer, bufsize);
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();
#endif
	
  /*--- Check the number of zones in the mesh file ---*/
	CConfig *config_container = NULL;
	if (argc == 2) config_container = new CConfig(argv[1]);
	else { strcpy (grid_file, "default.cfg"); config_container = new CConfig(grid_file); }
	nZone = GetnZone(config_container->GetMesh_FileName(), config_container->GetMesh_FileFormat(), config_container);
  
	/*--- Pointer to different structures that will be used throughout the entire code ---*/
	CConfig **config = NULL;
	CPhysicalGeometry **geometry = NULL;
	CFreeFormChunk** chunk = NULL;
	CSurfaceMovement *surface_mov = NULL;
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
      strcpy (grid_file, "default.cfg");
      config[iZone] = new CConfig(grid_file, SU2_MDC, iZone, nZone, VERB_HIGH);
    }
    
#ifndef NO_MPI
    /*--- Change the name of the input-output files for the
     parallel computation ---*/
    config[iZone]->SetFileNameDomain(rank+1);
#endif
    
    /*--- Definition of the geometry class ---*/
    geometry[iZone] = new CPhysicalGeometry(config[iZone], config[iZone]->GetMesh_FileName(),
                                            config[iZone]->GetMesh_FileFormat(), iZone+1, nZone);
  }
  
  /*--- Perform all actions in SU2_MDC on ZONE_0 only. ---*/
	if (rank == MASTER_NODE)
		cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
  
	/*--- Compute elements surrounding points, points surrounding points, and elements surronding elements ---*/
	if (rank == MASTER_NODE) cout << "Setting local point and element connectivity." <<endl;
	geometry[ZONE_0]->SetEsuP(); geometry[ZONE_0]->SetPsuP();
	
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
  
	if (rank == MASTER_NODE)
	  cout << endl <<"------------------------- Surface grid deformation ----------------------" << endl;
  
	/*--- Definition and initialization of the surface deformation class ---*/
	surface_mov = new CSurfaceMovement();
	surface_mov->CopyBoundary(geometry[ZONE_0], config[ZONE_0]);
  
	/*--- Definition of the FFD deformation class ---*/
	unsigned short nChunk = MAX_NUMBER_CHUNK;
	chunk = new CFreeFormChunk*[nChunk];
  
	/*--- Output original computational grids ---*/
	if (config[ZONE_0]->GetVisualize_Deformation()) {
		if(config[ZONE_0]->GetOutput_FileFormat() == PARAVIEW) {
			if (size > 1) sprintf (buffer_char, "_%d.vtk", rank+1);
			else sprintf (buffer_char, ".vtk");
			strcpy (out_file, "original_volumetric_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetParaView(out_file);
			strcpy (out_file, "original_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundParaView(config[ZONE_0], out_file);
		}
		if(config[ZONE_0]->GetOutput_FileFormat() == TECPLOT) {
			if (size > 1) sprintf (buffer_char, "_%d.plt", rank+1);
			else sprintf (buffer_char, ".plt");
			strcpy (out_file, "original_volumetric_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetTecPlot(out_file);
			strcpy (out_file, "original_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundTecPlot(config[ZONE_0], out_file);
		}
		if(config[ZONE_0]->GetOutput_FileFormat() == STL) {
			if (size > 1) sprintf (buffer_char, "_%d.stl", rank+1);
			else sprintf (buffer_char, ".stl");
			strcpy (out_file, "original_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundSTL(config[ZONE_0],out_file);
		}
	}
  
	/*--- Bump deformation for 2D problems ---*/
	if (geometry[ZONE_0]->GetnDim() == 2) {
		
		if (rank == MASTER_NODE) {
			if (config[ZONE_0]->GetDesign_Variable(0) != NO_DEFORMATION)
				cout << "Perform 2D deformation of the surface." << endl;
			else
				cout << "No deformation of the surface grid." << endl;
		}
    
		/*--- Apply the design variables to the control point position ---*/
		for (unsigned short iDV = 0; iDV < config[ZONE_0]->GetnDV(); iDV++) {
			switch ( config[ZONE_0]->GetDesign_Variable(iDV) ) {
				case HICKS_HENNE : surface_mov->SetHicksHenne(geometry[ZONE_0], config[ZONE_0], iDV, false); break;
				case DISPLACEMENT : surface_mov->SetDisplacement(geometry[ZONE_0], config[ZONE_0], iDV, false); break;
				case ROTATION : surface_mov->SetRotation(geometry[ZONE_0], config[ZONE_0], iDV, false); break;
				case NACA_4DIGITS : surface_mov->SetNACA_4Digits(geometry[ZONE_0], config[ZONE_0]); break;
				case PARABOLIC : surface_mov->SetParabolic(geometry[ZONE_0], config[ZONE_0]); break;
				case OBSTACLE : surface_mov->SetObstacle(geometry[ZONE_0], config[ZONE_0]); break;
				case STRETCH : surface_mov->SetStretch(geometry[ZONE_0], config[ZONE_0]); break;
        case SURFACE_FILE : surface_mov->SetExternal_Deformation(geometry[ZONE_0], config[ZONE_0], ZONE_0, iExtIter); break;
			}
		}
	}
  
	/*--- Free Form deformation or surface file input for 3D problems ---*/
	if (geometry[ZONE_0]->GetnDim() == 3) {
		
    if (config[ZONE_0]->GetDesign_Variable(0) == SURFACE_FILE) {
      
      /*--- Check whether a surface file exists for input ---*/
      ofstream Surface_File;
      string filename = config[ZONE_0]->GetMotion_FileName();
      Surface_File.open(filename.c_str(), ios::in);
      
      /*--- A surface file does not exist, so write a new one for the
       markers that are specified as part of the motion. ---*/
      if (Surface_File.fail()) {
        
        if (rank == MASTER_NODE)
          cout << "No surface file found. Writing a new file: " << filename << "." << endl;
        Surface_File.open(filename.c_str(), ios::out);
        Surface_File.precision(15);
        unsigned long iMarker, jPoint, GlobalIndex, iVertex; double *Coords;
        for (iMarker = 0; iMarker < config[ZONE_0]->GetnMarker_All(); iMarker++) {
          if (config[ZONE_0]->GetMarker_All_Moving(iMarker) == YES) {
            for(iVertex = 0; iVertex < geometry[ZONE_0]->nVertex[iMarker]; iVertex++) {
              jPoint = geometry[ZONE_0]->vertex[iMarker][iVertex]->GetNode();
              GlobalIndex = geometry[ZONE_0]->node[jPoint]->GetGlobalIndex();
              Coords = geometry[ZONE_0]->node[jPoint]->GetCoord();
              Surface_File << GlobalIndex << "\t" << Coords[0] << "\t" << Coords[1];
              if (geometry[ZONE_0]->GetnDim() == 2) Surface_File << endl;
              else Surface_File << "\t" << Coords[2] << endl;
              //else Surface_File << "\t" << Coords[2] + 2.5/30.0*Coords[1] << endl;
              //cout << Coords[2] << "  " <<  2.5/30.0*Coords[1]<< endl;
            }
          }
        }
        Surface_File.close();
        
        /*--- A surface file exists, so read in the coordinates ---*/
      } else {
        Surface_File.close();
        if (rank == MASTER_NODE)
          cout << "Updating the surface coordinates from the input file." << endl;
        surface_mov->SetExternal_Deformation(geometry[ZONE_0], config[ZONE_0], ZONE_0, iExtIter);
      }
      
		} else {
      
      if (rank == MASTER_NODE) {
        if (config[ZONE_0]->GetDesign_Variable(0) != NO_DEFORMATION)
          cout << "Performing the deformation of the surface grid." << endl;
        else
          cout << "No deformation of the surface grid." << endl;
      }
      
      /*--- Read the FFD information fron the grid file ---*/
      surface_mov->ReadFFDInfo(config[ZONE_0], geometry[ZONE_0], chunk, config[ZONE_0]->GetMesh_FileName());
      
      /*--- If the chunk was not defined in the input file ---*/
      if (!surface_mov->GetChunkDefinition()) {
        
        if (rank == MASTER_NODE)
          cout << endl <<"----------------- FFD technique (cartesian -> parametric) ---------------" << endl;
        
        /*--- Create a unitary chunk as baseline for other chunks shapes ---*/
        CFreeFormChunk chunk_unitary(1,1,1);
        chunk_unitary.SetUnitCornerPoints();
        
        /*--- Compute the control points of the unitary box, in this case the degree is 1 and the order is 2 ---*/
        chunk_unitary.SetControlPoints_Parallelepiped();
        
        for (iChunk = 0; iChunk < surface_mov->GetnChunk(); iChunk++) {
          /*--- Compute the support control points for the final FFD using the unitary box ---*/
          chunk_unitary.SetSupportCP(chunk[iChunk]);
          
          /*--- Compute control points in the support box ---*/
          chunk_unitary.SetSupportCPChange(chunk[iChunk]);
          
          /*--- Compute the parametric coordinates, it also find the points in
           the chunk using the parametrics coordinates ---*/
          surface_mov->SetParametricCoord(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], iChunk);
          
        }
        
      }
      
      /*--- Output original FFD chunk ---*/
      for (iChunk = 0; iChunk < surface_mov->GetnChunk(); iChunk++) {
        if(config[ZONE_0]->GetOutput_FileFormat() == PARAVIEW) {
          sprintf (buffer_char, "original_chunk.vtk");
          if (iChunk == 0) chunk[iChunk]->SetParaView(buffer_char, true);
          else chunk[iChunk]->SetParaView(buffer_char, false);
        }
        if(config[ZONE_0]->GetOutput_FileFormat() == TECPLOT) {
          sprintf (buffer_char, "original_chunk.plt");
          if (iChunk == 0) chunk[iChunk]->SetTecplot(buffer_char, true);
          else chunk[iChunk]->SetTecplot(buffer_char, false);
        }
      }
      
      if (rank == MASTER_NODE)
        cout << endl <<"----------------- FFD technique (parametric -> cartesian) ---------------" << endl;
      
      /*--- Loop over all the FFD boxes levels ---*/
      for (iLevel = 0; iLevel < surface_mov->GetnLevel(); iLevel++) {
        
        /*--- Loop over all FFD chunks ---*/
        for (iChunk = 0; iChunk < surface_mov->GetnChunk(); iChunk++) {
          
          /*--- Check the level of the FFD box ---*/
          if(chunk[iChunk]->GetLevel() == iLevel) {
            
            /*--- Compute the parametric coordinates of the child box
             control points (using the parent chunk)  ---*/
            for (iChild = 0; iChild < chunk[iChunk]->GetnChildChunk(); iChild++) {
              ChunkTag = chunk[iChunk]->GetChildChunkTag(iChild);
              for (jChunk = 0; jChunk < surface_mov->GetnChunk(); jChunk++)
                if (ChunkTag == chunk[jChunk]->GetTag()) break;
              surface_mov->SetParametricCoordCP(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], chunk[jChunk]);
            }
            
            /*--- Update the parametric coordinates if it is a child chunk ---*/
            if (iLevel > 0) surface_mov->UpdateParametricCoord(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], iChunk);
            
            /*--- Apply the design variables to the control point position ---*/
            for (iDV = 0; iDV < config[ZONE_0]->GetnDV(); iDV++) {
              switch ( config[ZONE_0]->GetDesign_Variable(iDV) ) {
                case FFD_CONTROL_POINT : surface_mov->SetFFDCPChange(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], iChunk, iDV, false); break;
                case FFD_DIHEDRAL_ANGLE : surface_mov->SetFFDDihedralAngle(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], iChunk, iDV, false); break;
                case FFD_TWIST_ANGLE : surface_mov->SetFFDTwistAngle(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], iChunk, iDV, false); break;
                case FFD_ROTATION : surface_mov->SetFFDRotation(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], iChunk, iDV, false); break;
                case FFD_CAMBER : surface_mov->SetFFDCamber(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], iChunk, iDV, false); break;
                case FFD_THICKNESS : surface_mov->SetFFDThickness(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], iChunk, iDV, false); break;
                case FFD_VOLUME : surface_mov->SetFFDVolume(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], iChunk, iDV, false); break;
              }
            }
            
            /*--- Recompute cartesian coordinates using the new control point location ---*/
            surface_mov->SetCartesianCoord(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], iChunk);
            
            /*--- Reparametrization of the parent FFD box ---*/
            for (iParent = 0; iParent < chunk[iChunk]->GetnParentChunk(); iParent++) {
              ChunkTag = chunk[iChunk]->GetParentChunkTag(iParent);
              for (jChunk = 0; jChunk < surface_mov->GetnChunk(); jChunk++)
                if (ChunkTag == chunk[jChunk]->GetTag()) break;
              surface_mov->UpdateParametricCoord(geometry[ZONE_0], config[ZONE_0], chunk[jChunk], jChunk);
            }
            
            /*--- Compute the new location of the control points of the child boxes
             (using the parent chunk) ---*/
            for (iChild = 0; iChild < chunk[iChunk]->GetnChildChunk(); iChild++) {
              ChunkTag = chunk[iChunk]->GetChildChunkTag(iChild);
              for (jChunk = 0; jChunk < surface_mov->GetnChunk(); jChunk++)
                if (ChunkTag == chunk[jChunk]->GetTag()) break;
              surface_mov->GetCartesianCoordCP(geometry[ZONE_0], config[ZONE_0], chunk[iChunk], chunk[jChunk]);
            }
          }
          
        }
      }
      
      /*--- Output the deformed chunks ---*/
      for (iChunk = 0; iChunk < surface_mov->GetnChunk(); iChunk++) {
        if (config[ZONE_0]->GetDesign_Variable(0) != NO_DEFORMATION) {
          if(config[ZONE_0]->GetOutput_FileFormat() == PARAVIEW) {
            sprintf (buffer_char, "deformed_chunk.vtk");
            if (iChunk == 0) chunk[iChunk]->SetParaView(buffer_char, true);
            else chunk[iChunk]->SetParaView(buffer_char, false);
          }
          if(config[ZONE_0]->GetOutput_FileFormat() == TECPLOT) {
            sprintf (buffer_char, "deformed_chunk.plt");
            if (iChunk == 0) chunk[iChunk]->SetTecplot(buffer_char, true);
            else chunk[iChunk]->SetTecplot(buffer_char, false);
          }
        }
      }
    }
	}
	
  /*--- MPI sincronization point ---*/
#ifndef NO_MPI
	MPI::COMM_WORLD.Barrier();
#endif
	
	if (rank == MASTER_NODE)
	  cout << endl <<"----------------------- Volumetric grid deformation ---------------------" << endl;
  
	/*--- Definition of the Class for grid movement ---*/
	grid_movement = new CVolumetricMovement(geometry[ZONE_0]);
  
	/*--- Grid deformation ---*/
	if (config[ZONE_0]->GetDesign_Variable(0) != NO_DEFORMATION) {
		if (rank == MASTER_NODE)
			cout << "Performing the deformation of the volumetric grid." << endl;
		if (config[ZONE_0]->GetKind_GridDef_Method() == ALGEBRAIC) grid_movement->AlgebraicMethod(geometry[ZONE_0], config[ZONE_0], false);
		if (config[ZONE_0]->GetKind_GridDef_Method() == SPRING) grid_movement->SpringMethod(geometry[ZONE_0], config[ZONE_0], false);
		if (config[ZONE_0]->GetKind_GridDef_Method() == TORSIONAL_SPRING) grid_movement->TorsionalSpringMethod(geometry[ZONE_0], config[ZONE_0], false);
	}
	else {
		if (rank == MASTER_NODE)
			cout << "No deformation deformation of the volumetric grid." << endl;
	}
	
	/*--- Output grid ---*/
	if (rank == MASTER_NODE)
		cout << "End and write output files." << endl;
	
	if (config[ZONE_0]->GetVisualize_Deformation() &&
			(config[ZONE_0]->GetDesign_Variable(0) != NO_DEFORMATION)) {
		if(config[ZONE_0]->GetOutput_FileFormat() == PARAVIEW) {
			if (size > 1) sprintf (buffer_char, "_%d.vtk", rank+1);
			else sprintf (buffer_char, ".vtk");
			strcpy (out_file, "deformed_volumetric_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetParaView(out_file);
			strcpy (out_file, "deformed_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundParaView(config[ZONE_0], out_file);
		}
		if(config[ZONE_0]->GetOutput_FileFormat() == TECPLOT) {
			if (size > 1) sprintf (buffer_char, "_%d.plt", rank+1);
			else sprintf (buffer_char, ".plt");
			strcpy (out_file, "deformed_volumetric_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetTecPlot(out_file);
			strcpy (out_file, "deformed_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundTecPlot(config[ZONE_0], out_file);
		}
		if(config[ZONE_0]->GetOutput_FileFormat() == STL) {
			if (size > 1) sprintf (buffer_char, "_%d.stl", rank+1);
			else sprintf (buffer_char, ".stl");
			strcpy (out_file, "deformed_surface_grid"); strcat(out_file, buffer_char); geometry[ZONE_0]->SetBoundSTL(config[ZONE_0], out_file);
		}
	}
	
  /*--- Check for multiple zones and call appropriate mesh file export routine. ---*/
  if (nZone == 1) {
    if (size > 1) sprintf (buffer_char, "_%d.su2", rank+1);
    else sprintf (buffer_char, ".su2");
    string str = config[ZONE_0]->GetMesh_Out_FileName();
    str.erase (str.end()-4, str.end()); strcpy (out_file, str.c_str()); strcat(out_file, buffer_char);
    geometry[ZONE_0]->SetMeshFile(config[ZONE_0], out_file);
	}
	else {
		/*--- Call special write routine for more than one zone. ---*/
    SetMultiZone_MeshFile(geometry, config, nZone);
  }
	
	if (geometry[ZONE_0]->GetnDim() == 3)
		surface_mov->WriteFFDInfo(geometry[ZONE_0], config[ZONE_0], chunk, out_file);
	
#ifndef NO_MPI
	/*--- Finalize MPI parallelization ---*/
	old_buffer = buffer;
	MPI::Detach_buffer(old_buffer);
	//	delete [] buffer;
	MPI::Finalize();
#endif
	
	/*--- End solver ---*/
	if (rank == MASTER_NODE)
	  cout << endl <<"------------------------- Exit Success (SU2_MDC) ------------------------" << endl << endl;
	
	return EXIT_SUCCESS;
	
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
						//					if (nZone == 1) cout << "SU2 mesh file format with a single zone." << endl;
						//					else if (nZone >  1) cout << "SU2 mesh file format with " << nZone << " zones." << endl;
						//					else
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
				//			if (rank == MASTER_NODE) cout << "SU2 mesh file format with a single zone." << endl;
			}
			break;
			
		case CGNS:
			
			nZone = 1;
			//		if (rank == MASTER_NODE) cout << "CGNS mesh file format with a single zone." << endl;
			break;
			
		case NETCDF_ASCII:
			
			nZone = 1;
			//		if (rank == MASTER_NODE) cout << "NETCDF mesh file format with a single zone." << endl;
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
