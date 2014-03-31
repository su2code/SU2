/*!
 * \file SU2_UMC.cpp
 * \brief Main file for UMarc Coupling (SU2_UMC).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.9
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

#include "../include/SU2_UMC.hpp"

int main(int argc, char *argv[]) {
		  
  /*--- Local variables ---*/
  unsigned short mode, nZone = 1;
  char grid_filename[200];
  bool complete;

  int rank = MASTER_NODE;
#ifndef NO_MPI
  /*--- MPI initialization, and buffer setting ---*/
  MPI::Init(argc,argv);
  rank = MPI::COMM_WORLD.Get_rank();
#endif

  /*--- Declare pointers to class objects ---*/
  CConfig *config = NULL;
  CBoundaryGeometry *boundary = NULL;
 
  /*--- Instatiate an object of the config class
        based on the name of the config file given ---*/
  strcpy(grid_filename, argv[1]);
  config = new CConfig(grid_filename, SU2_GDC, ZONE_0, nZone, 
                         0, VERB_HIGH);

  /*--- Instantiate an object of the boundary-geometry class ---*/
  boundary = new CBoundaryGeometry(config, 
                                   config->GetMesh_FileName(), 
                                   config->GetMesh_FileFormat());
  
  /*--- determine the mode of operation ---*/
  if (argc < 3) {
    if (argc == 2) {
      cout << "\t ------------------------------------------ \t" << endl;
      cout << "\t *** ERROR! PLEASE ENTER A RUNNING MODE *** \t" << endl;
      cout << "\t ------------------------------------------ \t" << endl;
      return EXIT_SUCCESS;
    }
    if (argc == 1) {
      cout << "\t -------------------------------------------------- \t" << endl;
      cout << "\t *** ERROR! PLEASE ENTER A CONFIG FILE AND MODE *** \t" << endl;
      cout << "\t -------------------------------------------------- \t" << endl;
      return EXIT_SUCCESS;
    }
  }
  mode = atoi(argv[2]);
  cout << endl;
  
  switch (mode) {
  
    case 1:
      cout << "\t ------------------------------------ \t" << endl;
      cout << "\t *** AIRFOIL-SECTION-CUTTING MODE *** \t" << endl;
      cout << "\t ------------------------------------ \t" << endl;
      
      /*--- if they exist, delete old versions of files "span," 
            "tracked_points," "Airfoil_Sections.plt," and 
            "coordinates.csv" ---*/
      cout << endl;
      clean_up("span");
      clean_up("tracked_points");
      clean_up("Airfoil_Sections.plt");
      clean_up("coordinates.csv");
      
      /*--- make the sectional cuts along the surfaces of interest and
            write the files "span," "tracked_points," and 
            "Airfoil_Sections.plt" to the current directory ---*/
      Make_Sectional_Cuts(config, boundary);
  
      /*--- find the coordinates of the sections and 
            print to the "coordinate.csv" file ---*/
      Tracked_Points_to_Coords(boundary);
      
      break;
      
    case 2:
      cout << "\t ------------------------------------------ \t" << endl;
      cout << "\t *** SECTIONAL-FORCES MODE (SU2->UMARC) *** \t" << endl;
      cout << "\t ------------------------------------------ \t" << endl;
      
      /*--- recompute the coordinates and rewrite "coordinates.csv" 
            file, just in case the geometry has been deformed since 
            running mode 1. N.B. The original cuts can still be found
            in "Airfoil_Sections.plt" ---*/
      cout << endl;
      clean_up("coordinates.csv");
      Tracked_Points_to_Coords(boundary);
      
      /*--- read the new coordinates, find the corresponding values 
            of pressure from "surface_flow.csv," and write the 
            sectional-forces files cl, cd, and cm ---*/
      complete = Compute_Sectional_Forces();
     

 
      
      
      if (complete == false) {
        cout << endl << "Sectional forces have not been computed properly." << endl;
        return EXIT_SUCCESS;
      }
      
      break;
      
    case 3:
      cout << "\t --------------------------------------------------- \t" << endl;
      cout << "\t *** DEFORMATION-COMMUNICATION MODE (UMARC->SU2) *** \t" << endl;
      cout << "\t --------------------------------------------------- \t" << endl;

      break;
      
    default:
      cout << "\t -------------------------------------------- \t" << endl;
      cout << "\t *** ERROR! PLEASE RUN IN MODE 1, 2, OR 3 *** \t" << endl;
      cout << "\t -------------------------------------------- \t" << endl;
      return EXIT_SUCCESS;
  }
                                   
  /*--- Delete dynamically allocated memory ---*/
	if (config != NULL) delete config;
	if (boundary != NULL) delete boundary;
	
	/*--- End routine ---*/
	cout << endl <<"------------------------- Exit Success (SU2_UMC) ------------------------" << endl << endl;
  
	return EXIT_SUCCESS;
	
}

void Make_Sectional_Cuts(CConfig *config, CBoundaryGeometry *boundary) {
  
  /*--- Local variables ---*/
  unsigned short nDim, nMarker, nSection, iDim, iMarker, iSection, 
                 Boundary;
  unsigned long iNode, nNode;
  vector<unsigned long> **point1_Airfoil, **point2_Airfoil;
  double root, tip, minPlane, maxPlane, radial_station, interval,
         starting_angle, starting_rad, rotation, magnitude;
  double starting_vec [3]; 
  double *Plane_P0, *Plane_Normal;
  vector<double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, 
                 **weight1_Airfoil;
  bool Monitoring, CCW_orientation;
  string Marker_Tag, filename;
  ofstream span, tracked_points;
  
  /*--- intialize the various pointers to null ---*/
  Xcoord_Airfoil = NULL;
  Ycoord_Airfoil = NULL;
  Zcoord_Airfoil = NULL;
  point1_Airfoil = NULL;
  point2_Airfoil = NULL;
  weight1_Airfoil = NULL;
  Plane_P0 = NULL;
  Plane_Normal = NULL;
  
  /*--- find the number of dimensions in the problem ---*/
  nDim = boundary->GetnDim();
  
  /*--- find the number of markers in the mesh ---*/
  nMarker = boundary->GetnMarker();
  
  if (nDim == 3) {
    
    /*--- USER-DEFINED: set the total number 
          of slices for all markers ---*/
    nSection = 15;
    
    /*--- USER-DEFINED: provide the locations, in the spanwise 
          direction, of the blade/wing root and tip ---*/
    root = 0.0;
    tip = 14.0;
     
    /*--- USER-DEFINED: set the starting and ending
          points of the slices (N.B. These need not 
          be the same as the root and tip locations) ---*/
    minPlane = root;
    maxPlane = tip;
		
    /*--- initialize coordinate vectors --*/
    Xcoord_Airfoil = new vector<double> [nSection];
    Ycoord_Airfoil = new vector<double> [nSection];
    Zcoord_Airfoil = new vector<double> [nSection];
    
    /*--- iniitialize the 2D arrays of vectors ---*/
    point1_Airfoil = new vector<unsigned long> *[nMarker];
    point2_Airfoil = new vector<unsigned long> *[nMarker];
    weight1_Airfoil = new vector<double> *[nMarker];
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      point1_Airfoil[iMarker] = new vector<unsigned long> [nSection];
      point2_Airfoil[iMarker] = new vector<unsigned long> [nSection];
      weight1_Airfoil[iMarker] = new vector<double> [nSection];
    }
    
    /*--- For now, the slicing planes are hard-coded for each 
          marker, but this will be added to the config file in the 
          future. ---*/
    Plane_P0 = new double [3];
    Plane_Normal = new double [3];
   
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      Boundary   = config->GetMarker_All_Boundary(iMarker);
      Monitoring = config->GetMarker_All_Monitoring(iMarker);
      if ((Boundary == EULER_WALL) ||
          (Boundary == HEAT_FLUX) ||
          (Boundary == ISOTHERMAL) || 
          (Boundary == NEARFIELD_BOUNDARY)) {
        if (Monitoring) {
          
          /*--- retreive the name of the current marker ---*/
          Marker_Tag = config->GetMarker_All_Tag(iMarker);

          /*--- USER-DEFINED: angle (in degrees) 
                the blade/wing is originally at 
                N.B. This only works if you're in the X-Y plane---*/
          
          /*--- ONERA M6 WING ---*/
	        if (Marker_Tag == "WING") {
            starting_angle = 90.0;
            root = 0.0;
            tip = 14.0;
            minPlane = root;
            maxPlane = tip;
            CCW_orientation = 1;
          }
          
          /*--- UH-60 FIRST blade ---*/
	        if (Marker_Tag == "Blade1") {
	          starting_angle = 0.0;
	          root = 1.0;
	          tip = 8.2;
	          minPlane = 1.5;
	          maxPlane = 7.5;
	          CCW_orientation = 1;
	        }
          
          /*--- UH-60 SECOND blade ---*/
	        if (Marker_Tag == "Blade2") {
	          starting_angle = 90.0;
	          root = 1.0;
	          tip = 8.2;
	          minPlane = 1.5;
	          maxPlane = 7.5;
	          CCW_orientation = 1;
	        }
          
          /*--- UH-60 THIRD blade ---*/
	        if (Marker_Tag == "Blade3") {
	          starting_angle = 180.0;
	          root = 1.0;
	          tip = 8.2;
	          minPlane = 1.5;
	          maxPlane = 7.5;
	          CCW_orientation = 1;
	        }
          
          /*--- UH-60 FOURTH blade ---*/
	        if (Marker_Tag == "Blade4") {
	          starting_angle = 270.0;
	          root = 1.0;
	          tip = 8.2;
	          minPlane = 1.5;
	          maxPlane = 7.5;
	          CCW_orientation = 1;
	        }
	  
          /*--- Caradonna-Tung FIRST blade ---*/
          if (Marker_Tag == "blade_1") {
            starting_angle = 90.0;
            root = 0.15;
            tip = 1.14;
            minPlane = root;
            maxPlane = tip;
            CCW_orientation = 0;
          }
          
          /*--- Caradonna-Tung SECOND blade ---*/
	        if (Marker_Tag == "blade_2") {
            starting_angle = 270.0;
            root = 0.15;
            tip = 1.14;
            minPlane = root;
            maxPlane = tip;
            CCW_orientation = 0;
          }
					
	        /*--- USER INPUTS ABOUT GEOMETRY END HERE ---*/
	        
          /*--- comupte the distance between sections ---*/
          interval = (maxPlane-minPlane)/double(nSection-1);
          
          /*--- write the "span" file, listing the radial 
                stations along the geometry being sliced ---*/
          span.precision(15);
          span.open("span", ios::out);
          for (iSection = nSection; iSection > 0; iSection--) {
            radial_station = minPlane + (iSection-1)*interval;
            span << (radial_station-root)/(tip-root) << endl;
          }
          span.close();
          
          /*--- using the user-defined inputs, compute the
                points and planes that define the sections ---*/
          
          /*--- convert the starting angle to radians ---*/
          starting_rad = starting_angle*(PI_NUMBER/180.0);
	        
          /*--- section cutting ---*/
          
          /*--- compute the vector pointing 
                along the blade/wing ---*/
	        starting_vec[0] = cos(starting_rad);
	        starting_vec[1] = sin(starting_rad);
	        starting_vec[2] = 0.0;
          
          /*--- USER DEFINED: if, for whatever reason, the mesh is 
                rotated from the starting angles supplied above, 
                enter that rotation in degrees here. ---*/
          rotation = 0.0;
          
          /*--- find the vector pointing along this blade  ---*/
	        Plane_Normal[0] = cos(rotation)*starting_vec[0] - 
	                          sin(rotation)*starting_vec[1] + 
	                          0.0*starting_vec[2];
	        Plane_Normal[1] = sin(rotation)*starting_vec[0] + 
	                          cos(rotation)*starting_vec[1] + 
	                          0.0*starting_vec[2];
	        Plane_Normal[2] = 0.0*starting_vec[0] + 
	                          0.0*starting_vec[1] + 
	                          1.0*starting_vec[2];
	        
          /*--- turn it into a unit vector ---*/
	        magnitude = sqrt(pow(Plane_Normal[0],2.0) + 
	                         pow(Plane_Normal[1],2.0) + 
	                         pow(Plane_Normal[2],2.0));
	        for (iDim = 0; iDim < nDim; iDim++) {
	          Plane_Normal[iDim] /= magnitude;
	          if (fabs(Plane_Normal[iDim]) < 1e-5) {
	            Plane_Normal[iDim] = 0.0;
	          }
	        }
	        
	        /*--- print marker information to the screen ---*/
          cout << endl << "Marker: " << Marker_Tag << endl;
          cout << "  Points along: {";
          cout << Plane_Normal[0] << ", ";
          cout << Plane_Normal[1] << ", ";
          cout << Plane_Normal[2] << "}" << endl;
          cout << "  Sectional cuts pass through the points:" << endl;
          
          
          
          /*--- walk along each marker, making sectional cuts ---*/
          for (iSection = 0; iSection < nSection; iSection++) {
            
            /*--- find the point that defines the cutting plane ---*/
	          for (iDim = 0; iDim < nDim; iDim++) {
	            Plane_P0[iDim] = Plane_Normal[iDim]*(minPlane +
	                             (double)iSection*interval);
	          }
	          
	          cout << "    Section #" << iSection << ": {";
	          cout << Plane_P0[0] << ", "; 
	          cout << Plane_P0[1] << ", ";
	          cout << Plane_P0[2] << "}" << endl;
            
//            /*--- find the points defining the airfoil section --*/
//	          boundary->ComputeAirfoil_Section(Plane_P0, Plane_Normal,
//	                                          iSection, config,
//                                            Xcoord_Airfoil[iSection],
//                                            Ycoord_Airfoil[iSection],
//                                            Zcoord_Airfoil[iSection],
//                                   point1_Airfoil[iMarker][iSection],
//                                   point2_Airfoil[iMarker][iSection],
//                                  weight1_Airfoil[iMarker][iSection],
//                                              true, CCW_orientation);
            
            /*--- writing of the tracked-points file.this file 
                  records the node number of the edges being 
                  interescted as well as the correspond "weight" of 
                  the first point, based on distance. ---*/
            
            /*--- write the header ---*/    
	          if (iMarker == 0 && iSection == 0) {  
              tracked_points.precision(15);
	            filename = "tracked_points";
	            tracked_points.open(filename.c_str(), ios::out);
	            tracked_points <<  "Marker_Name, Marker_Number, ";
	            tracked_points << "Section, ";
	            tracked_points << "point1, point2, weight1" << endl;
	            tracked_points.close();
            }
            
	          /*--- find the number of nodes for 
                  the current airfoil section ---*/ 
	          nNode = point1_Airfoil[iMarker][iSection].size();
            
	          /*--- print information to filefile for a given node 
                      at a given section on a given marker ---*/
            tracked_points.open(filename.c_str(), ios::out | ios::app);
	          for (iNode = 0; iNode < nNode; iNode++) {
	            tracked_points << Marker_Tag << ", " << iMarker << ", ";
	            tracked_points << iSection << ", ";
	            tracked_points << point1_Airfoil[iMarker][iSection].at(iNode) << ", ";
              tracked_points << point2_Airfoil[iMarker][iSection].at(iNode) << ", ";
	            tracked_points << weight1_Airfoil[iMarker][iSection].at(iNode) << endl;
            }
            
            /*--- close the file tracking points---*/
            tracked_points.close();
            
          } /*--- end loop over sections ---*/
          
	      }  
      }    	    
    }

    /*--- delete the dynamically allocated memory  ---*/
    if (Xcoord_Airfoil != NULL) {delete [] Xcoord_Airfoil;}
    if (Ycoord_Airfoil != NULL) {delete [] Ycoord_Airfoil;}
    if (Zcoord_Airfoil != NULL) {delete [] Zcoord_Airfoil;}
    if (Plane_P0 != NULL) {delete [] Plane_P0;}
    if (Plane_Normal != NULL) {delete [] Plane_Normal;}
    if (point1_Airfoil != NULL) {
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        for (iSection = 0; iSection < nSection; iSection++) {
          point1_Airfoil[iMarker][iSection].clear();
        }
        delete [] point1_Airfoil[iMarker];
      }
      delete [] point1_Airfoil;
    }
    if (point2_Airfoil != NULL) {
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        for (iSection = 0; iSection < nSection; iSection++) {
          point2_Airfoil[iMarker][iSection].clear();
        }
        delete [] point2_Airfoil[iMarker];
      }
      delete [] point2_Airfoil;
    }
    if (weight1_Airfoil != NULL) {
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        for (iSection = 0; iSection < nSection; iSection++) {
          weight1_Airfoil[iMarker][iSection].clear();
        }
			  delete [] weight1_Airfoil[iMarker];
		  }
		  delete [] weight1_Airfoil;
    }
    
  }
	else {
		cout << "WARNING: SU2_UMC cannot be run for one- or two-dimensional problems" << endl;
	}
	
}

void Tracked_Points_to_Coords(CBoundaryGeometry *boundary) {
  
  /*--- local variables ---*/
  fstream tp_file;
  ofstream coords_file;
  string line;
  unsigned short counter, iDim, nDim;
  unsigned long point1, point2;
  double weight1, weight2, x_coord, y_coord, z_coord;
  double *coord1, *coord2;
  bool dummy;
  
  /*--- preliminaries ---*/
  nDim = boundary->GetnDim();
  coord1 = NULL;
  coord2 = NULL;
  
  /*--- dynamically allocate memory ---*/
  coord1 = new double [nDim];
  coord2 = new double [nDim];
  
  /*--- read the tracked points file ---*/
  tp_file.open("tracked_points");
  coords_file.open("coordinates.csv", ios::out | ios::app);
  coords_file.precision(15);
  
  if (tp_file.is_open() && coords_file.is_open()) {
    
    /*--- read the header ---*/
    getline(tp_file,line);
    
    /*--- read through each line ---*/
    while (getline(tp_file,line)) {
    
      
      /*--- turn it into a char array ---*/
      char tmp [500];
      strcpy(tmp,line.c_str());
      
      /*--- read the line, storing point1, point2, and weight1 ---*/
      char *pointer_to = NULL;
      pointer_to = strtok(tmp," (,)\r\n");
      
      counter = 1;
      while (pointer_to != NULL) {
        
        switch (counter) {
          case 4:
            point1 = strtoul(pointer_to,NULL,0);
            break;
          case 5:
            point2 = strtoul(pointer_to,NULL,0);
            break;
          case 6:
            weight1 = strtod(pointer_to,NULL);
            break;
          default:
            dummy = 0;
        }
        pointer_to = strtok(NULL," (,)\r\n");
        counter++;
      }
        
      /*--- convert the points to coordinates ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        coord1[iDim] = boundary->node[point1]->GetCoord(iDim);
        coord2[iDim] = boundary->node[point2]->GetCoord(iDim);
      }
      
      /*--- retreive & compute the weight of the second node ---*/
      weight2 = 1.0 - weight1;
      
      /*--- find the new coordinates of the section---*/
      x_coord = weight1*coord1[0] + weight2*coord2[0];
      y_coord = weight1*coord1[1] + weight2*coord2[1];
      z_coord = weight1*coord1[2] + weight2*coord2[2];
        
      /*--- write to the coordinates file ---*/
      coords_file << x_coord << ", ";
      coords_file << y_coord << ", ";
      coords_file << z_coord << endl;
        
    }
    
    /*--- close the files ---*/
    tp_file.close();
    coords_file.close();
    
  }
  else {
    if (tp_file.is_open())
      cout << "Unable to open \"coordinates.csv\" file!" << endl;
    if (coords_file.is_open())
      cout << "Unable to open \"tracked points\" file!" << endl;
  }
  
  /*--- delete dynamically allocated memory ---*/
  if (coord1 != NULL) {delete [] coord1;}
  if (coord2 != NULL) {delete [] coord2;}
  
}

void clean_up(const char *filename) {

  /*--- if the exists, delete it ---*/
  if (fexists(filename) == true) {
    if(remove(filename) != 0) {
      cout << "Error: Old copy of \"" << filename; 
      cout << "\" not deleted!" << endl;
    }
    else {
      cout << "N.B. Old copy of \"" << filename;
      cout << "\" deleted." << endl;
    }
  }

}

bool fexists(const char *filename) {

  /*--- TRUE if the file exists ---*/
  ifstream ifile(filename);
  return ifile;

}

bool Compute_Sectional_Forces() {
  
  /*--- local variables ---*/
  fstream span_file, tp_file, coords_file;
  unsigned short nSection, nMarker, old_size, new_size,
                 i, iMarker, iSection, counter, iDim, nDim=3,
                 iSegment, nSegment;
  unsigned long the_section, point1, point2, nNodes, iNode;
  vector<unsigned long> *point1_Airfoil, *point2_Airfoil;
  double weight1, coord, maxDistance, distance, chord,
         x_node, y_node, z_node, first_line_x, first_line_y, first_line_z;
  double *normal_force, *chord_force, *pitching_moment, 
         *cl, *cd, *cm;
  double trailing_node [3], leading_node [3], chord_line [3], 
         q_chord [3], first_line [3], section_normal [3], 
         normal_line [3];
  vector<double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, 
                 *weight1_Airfoil;
  string line, previous, current_marker, row, the_marker;
  vector <string> markers;
  char tmp[200];
  char *pointer_to;
  int line_number;
  bool dummy;
  
  
  /*--- intialize the various pointers to null ---*/
  pointer_to = NULL;
  Xcoord_Airfoil = NULL;
  Ycoord_Airfoil = NULL;
  Zcoord_Airfoil = NULL;
  point1_Airfoil = NULL;
  point2_Airfoil = NULL;
  weight1_Airfoil = NULL;
  normal_force = NULL;
  chord_force = NULL;
  pitching_moment = NULL;
  cl = NULL;
  cd = NULL;
  cm = NULL;
    
  /*--- check for the "span" file ---*/
  if (fexists("span") == false) {
    cout << endl;
    cout << "*** ERROR! \"span\" does not exist! ***" << endl;
    cout << "*** Please run mode 1 before running mode 2 ***" << endl;
    return 0;
  }
  
  /*--- check for the "tracked_points" file ---*/
  if (fexists("tracked_points") == false) {
    cout << endl;
    cout << "*** ERROR! \"tracked_points\" does not exist! ***" << endl;
    cout << "*** Please run mode 1 before running mode 2 ***" << endl;
    return 0;
  }
  
  /*--- check for the "coordinates.csv" file ---*/
  if (fexists("coordinates.csv") == false) {
    cout << endl;
    cout << "*** ERROR! \"coordinates.csv\" does not exist! ***" << endl;
    cout << "*** Please check subroutine Tracked_Points_to_Coords ***" << endl;
    return 0;
  }
  
  /*--- find the number of sections by reading the "span" file ---*/
  nSection = 0;
  span_file.open("span");
  if (span_file.is_open()) {
    while (getline(span_file,line)) {
      nSection++;
    }
  }
  else
    cout << endl << "Unable to open \"span\" file!" << endl;
  span_file.close();
  
  /*--- print number of cuts to the screen ---*/
  if (nSection == 1)
    cout << endl << "There is one sectional cut ";
  else
    cout << endl << "There are " << nSection << " sectional cuts ";
  cout << "along each marker." << endl;
  
  /*--- initialize coordinate vectors --*/
  Xcoord_Airfoil = new vector<double> [nSection];
  Ycoord_Airfoil = new vector<double> [nSection];
  Zcoord_Airfoil = new vector<double> [nSection];
  
  /*--- initialize the interpolating vectors ---*/
  point1_Airfoil = new vector<unsigned long> [nSection];
  point2_Airfoil = new vector<unsigned long> [nSection];
  weight1_Airfoil = new vector<double> [nSection];
  
  /*--- dynamically-defined arrays for airfoil normal
     force, chord force, and pitching moment  -
     and their nondimensional coefficients ---*/
  normal_force = new double [nSection];
  chord_force = new double [nSection];
  pitching_moment = new double [nSection];
  cl = new double [nSection];
  cd = new double [nSection];
  cm = new double [nSection];
  
  /*--- find the number of markers in the "tracked_points" file ---*/
  tp_file.open("tracked_points");
  if (tp_file.is_open()) {

    /*--- read through the file, reading the marker names ---*/
    cout << endl << "Markers found in \"tracked_points.\"" << endl;
    line_number = 1;
    old_size = 0;
    while (getline(tp_file,line)) {
      
      /*--- if the marker name changes, make a record of it ---*/
      if (line_number != 1) {
        
        strcpy(tmp,line.c_str());
        pointer_to = strtok(tmp," ,\r\n");
        
        if (line_number == 2)
          markers.push_back(pointer_to);
        else {
          previous = markers.at(markers.size()-1);
          if (pointer_to != previous)
            markers.push_back(pointer_to);
        }
        
        /*--- increment the line number ---*/
        line_number++;
      }
      else {
        /*--- skip the header ---*/
        line_number++;
      }
      
      /*--- print marker-name information to the screen --*/
      new_size = markers.size();
      if (old_size != new_size) {
        cout << "  -at line #" << line_number << ": ";
        cout << "( ";
        for (i = 0; i < markers.size(); i++ ) {
          cout << markers.at(i);
          if (i < markers.size()-1)
            cout << ", ";
          else
            cout << " )" << endl;
        }
        old_size = new_size;
      }
    }
  }
  else
    cout << endl << "Unable to open \"tracked_points\" file!" << endl;
  tp_file.close();
  nMarker = markers.size();
  
  
  /*--- loop over the markers ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
        
    /*--- identify the name of the marker we're interested in ---*/
    current_marker = markers.at(iMarker);
    
    /*--- since we've started working with a new marker, 
          clear the arrays of coordinate vectors ---*/
    for (iSection = 0; iSection < nSection; iSection++) {
      Xcoord_Airfoil[iSection].clear();
      Ycoord_Airfoil[iSection].clear();
      Zcoord_Airfoil[iSection].clear();
    }
    
    /*--- read from the "tracked_points" file and from the 
          "coordinates.csv" to populate the coordinate vectors ---*/
    tp_file.open("tracked_points");
    coords_file.open("coordinates.csv");
    if (tp_file.is_open() && coords_file.is_open()) {
      
      /*--- read the header of "tracked_points" ---*/
      getline(tp_file,line);
      
      /*--- read through the rows in "coordinates.csv" ---*/
      while (getline(coords_file,row)) {
      
        /*--- read the corresponding line in "tracked_points" ---*/
        getline(tp_file,line);
        
        /*--- turn it into a char array ---*/
        strcpy(tmp,line.c_str());
        
        /*--- read the line, storing the marker name, 
              section number, point1, point2, and weight1 ---*/
        pointer_to = strtok(tmp," ,\r\n");
        
        counter = 1;
        while (pointer_to != NULL) {
          
          switch (counter) {
            case 1:
              the_marker = pointer_to;
              break;
            case 3:
              the_section = strtoul(pointer_to,NULL,0);
              break;
            case 4:
              point1 = strtoul(pointer_to,NULL,0);
              point1_Airfoil[the_section].push_back(point1);
              break;
            case 5:
              point2 = strtoul(pointer_to,NULL,0);
              point2_Airfoil[the_section].push_back(point2);
              break;
            case 6:
              weight1 = strtod(pointer_to,NULL);
              weight1_Airfoil[the_section].push_back(weight1);
              break;
            default:
              dummy = 0;
          }
          
          /*--- move to the next entry in the line ---*/
          pointer_to = strtok(NULL," ,\r\n");
          
          /*--- increment the counter ---*/
          counter++;
          
        } /*--- end reading line in "tracked_points" ---*/
        
        /*--- for all the intersection points along this marker, 
              start filling up the arrays of coordinate vectors ---*/
        if (the_marker == current_marker) {
          
          /*--- turn the row into a char array ---*/
          strcpy(tmp,row.c_str());
          
          /*--- read the line, storing the marker name, 
                section number, point1, point2, and weight1 ---*/
          pointer_to = strtok(tmp," ,\r\n");
          
          counter = 1;
          while (pointer_to != NULL) {
            
            /*--- convert to a double ---*/
            coord = strtod(pointer_to,NULL);
            
            switch (counter) {
              case 1:
                Xcoord_Airfoil[the_section].push_back(coord);
                break;
              case 2:
                Ycoord_Airfoil[the_section].push_back(coord);
                break;
              case 3:
                Zcoord_Airfoil[the_section].push_back(coord);
                break;
              default:
                dummy = 0;
            }
            
            /*--- move to the next entry in the row ---*/
            pointer_to = strtok(NULL," ,\r\n");
            
            /*--- increment the counter ---*/
            counter++;
            
          } /*--- end reading line in "coordinates.csv" ---*/
        
        } /*--- end if (the_marker == current_marker) ---*/
        
      } /*--- end reading rows in "coordinates.csv" ---*/
    
    } /*--- end if (the files are open) ---*/
    else {
      if (tp_file.is_open())
        cout << "Unable to open \"coordinates.csv\" file!" << endl;
      if (coords_file.is_open())
        cout << "Unable to open \"tracked points\" file!" << endl;
    }
    tp_file.close();
    coords_file.close();
    
    /*--- loop over the sections and find the integrated forces ---*/
    for (iSection = 0; iSection < nSection; iSection++) {
    
      /*--- find the number of intersection points 
            defining this arifoil section ---*/
      nNodes = Xcoord_Airfoil[iSection].size();
      
      /*--- find the node that corresponds to the trailing edge ---*/
      trailing_node[0] = Xcoord_Airfoil[iSection].at(0);
      trailing_node[1] = Ycoord_Airfoil[iSection].at(0);
      trailing_node[2] = Zcoord_Airfoil[iSection].at(0);
      
      /*--- find the node that corresponds to the leading edge ---*/
      maxDistance = 0;
      for (iNode = 1; iNode < nNodes; iNode++) {
      
        /*--- compute the distance from the trailing edge ---*/
        x_node = Xcoord_Airfoil[iSection].at(iNode);
        y_node = Ycoord_Airfoil[iSection].at(iNode);
        z_node = Zcoord_Airfoil[iSection].at(iNode);
        distance = sqrt(pow(trailing_node[0]-x_node,2.0) 
                   + pow(trailing_node[1]-y_node,2.0) 
                   + pow(trailing_node[2]-z_node,2.0));
              
        /*--- assume the leading edge is the point
              farthest away from the trailing edge ---*/
        if (distance > maxDistance) {
          maxDistance = distance;
          leading_node[0] = x_node;
          leading_node[1] = y_node;
          leading_node[2] = z_node;
        }
      }
      
      /*--- compute the length of the chord ---*/
      chord = sqrt(pow(trailing_node[0]-leading_node[0],2.0) 
                   + pow(trailing_node[1]-leading_node[1],2.0) 
                   + pow(trailing_node[2]-leading_node[2],2.0));
      
      /*--- define vector pointing from the trailing 
            edge to the leading edge ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        chord_line[iDim] = leading_node[iDim]-trailing_node[iDim];
      }
      
      /*--- find the coordinates of the quarter-chord point
            (move three-quarters of the chord
            length along the chord line) ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        q_chord[iDim] = trailing_node[iDim]+0.75*chord_line[iDim];
      }
      
      /*--- find the vector that points from the trailing node to
             the first point along the top of the airfoil surface
             (assume that this is always the
             second point of the coord vectors.) ---*/
      first_line_x = Xcoord_Airfoil[iSection].at(1);
      first_line_y = Ycoord_Airfoil[iSection].at(1);
      first_line_z = Zcoord_Airfoil[iSection].at(1);
      first_line[0] = first_line_x - trailing_node[0];
      first_line[1] = first_line_y - trailing_node[1];
      first_line[2] = first_line_z - trailing_node[2];
            
      /*--- find the vector that lies normal to the plane of the 
            airfoil (we are assuming, that after a deformation, the 
            plane's normal will no longer be oriented in the 
            conventional +Y direction.) ---*/
      section_normal[0] = chord_line[1]*first_line[2] - chord_line[2]*first_line[1];
      section_normal[1] = chord_line[2]*first_line[0] - chord_line[0]*first_line[2];
      section_normal[2] = chord_line[0]*first_line[1] - chord_line[1]*first_line[0];
            
      /*--- find the direction of the normal force (take a cross 
            product between the section normal and the chordline) ---*/
      normal_line[0] = section_normal[1]*chord_line[2] - section_normal[2]*chord_line[1];
      normal_line[1] = section_normal[2]*chord_line[0] - section_normal[0]*chord_line[2];
      normal_line[2] = section_normal[0]*chord_line[1] - section_normal[1]*chord_line[0];
      
      /*--- define the number of segments along the section ---*/
      nSegment = Xcoord_Airfoil[iSection].size()-1;
      
      /*--- Information about each segment along the airfoil section 
            will be stored in a structure of type segment_t ---*/
      struct segment_t {
        double endpoint1 [3];          // coordinates of the first endpoint
        double endpoint2 [3];          // coordinates of the second endpoint
        double length;                 // segment length
        double midpoint [3];           // coordinates of the segment midpoint
        double normal [3];             // components of the outward-pointing unit normal
        double normal_force;           // force per unit length projected 90 degs above chordline
        double chord_force;            // force per unit length projected along the chordline
        double pitching_moment;        // force/unit length times lever arm (midpoint to quarter chord)
        double pressure_endpoint1;     // pressure at point 1
        double pressure_endpoint2;     // pressure at point 2
      };
      
      segment_t *segments;                  // pointer of type segment_t
      segments = new segment_t [nSegment];  // array of structures
      
      /*--- append the first entry (corresponding to the trailing
            point) to the ends of the coord vectors and the 
            interpolating vectors so that we have easily identifiable 
            segments all around the airfoil section. ---*/
      vector<double> Xcoord_around, Ycoord_around, Zcoord_around, weight1_around;
      vector<unsigned long> point1_around, point2_around;
      
      Xcoord_around = Xcoord_Airfoil[iSection];
      Xcoord_around.push_back(Xcoord_Airfoil[iSection].at(0));
      
      Ycoord_around = Ycoord_Airfoil[iSection];
      Ycoord_around.push_back(Ycoord_Airfoil[iSection].at(0));
      
      Zcoord_around = Zcoord_Airfoil[iSection];
      Zcoord_around.push_back(Zcoord_Airfoil[iSection].at(0));
      
      point1_around = point1_Airfoil[iSection];
      point1_around.push_back(point1_Airfoil[iSection].at(0));
      
      point2_around = point2_Airfoil[iSection];
      point2_around.push_back(point2_Airfoil[iSection].at(0));
      
      weight1_around = weight1_Airfoil[iSection];
      weight1_around.push_back(weight1_Airfoil[iSection].at(0));
      
      /*--- for each segment around the airfoil, populate a data structure ---*/
      for (iSegment = 0; iSegment < nSegment; iSegment++) {
        
        /*--- coordinates of first endpoint ---*/
        segments[iSegment].endpoint1[0] = Xcoord_around.at(iSegment);
        segments[iSegment].endpoint1[1] = Ycoord_around.at(iSegment);
        segments[iSegment].endpoint1[2] = Zcoord_around.at(iSegment);
        
        /*--- coordinates of second endpoint ---*/
        segments[iSegment].endpoint2[0] = Xcoord_around.at(iSegment+1);
        segments[iSegment].endpoint2[1] = Ycoord_around.at(iSegment+1);
        segments[iSegment].endpoint2[2] = Zcoord_around.at(iSegment+1);
        
        /*--- segment length ---*/
        double deltaX = segments[iSegment].endpoint2[0] - segments[iSegment].endpoint1[0];
        double deltaY = segments[iSegment].endpoint2[1] - segments[iSegment].endpoint1[1];
        double deltaZ = segments[iSegment].endpoint2[2] - segments[iSegment].endpoint1[2];
        segments[iSegment].length = sqrt(pow(deltaX,2.0) + pow(deltaY,2.0) + pow(deltaZ,2.0));
        
        /*--- coordinates of the segment midpoint ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          segments[iSegment].midpoint[iDim] = 0.5*(segments[iSegment].endpoint1[iDim] 
                                                 + segments[iSegment].endpoint2[iDim]);
        }
        
        /*--- compute the unit "tangent" vector, i.e.
               the vector pointing along the segment ---*/
        double tangent [3];
        tangent[0] = (1.0/segments[iSegment].length)*(deltaX);
        tangent[1] = (1.0/segments[iSegment].length)*(deltaY);
        tangent[2] = (1.0/segments[iSegment].length)*(deltaZ);
        
        /*--- compute the outward-pointing normal ---*/
        segments[iSegment].normal[0] = section_normal[1]*tangent[2] - section_normal[2]*tangent[1];
        segments[iSegment].normal[1] = section_normal[2]*tangent[0] - section_normal[0]*tangent[2];
        segments[iSegment].normal[2] = section_normal[0]*tangent[1] - section_normal[1]*tangent[0];
        
        /*--- compute the outward-pointing unit normal ---*/
        double norm_dot_norm = pow(segments[iSegment].normal[0],2.0) 
                               + pow(segments[iSegment].normal[1],2.0) 
                               + pow(segments[iSegment].normal[2],2.0);
        double normal_length = sqrt(norm_dot_norm);
        for (iDim = 0; iDim < nDim; iDim++) {
          segments[iSegment].normal[iDim] = segments[iSegment].normal[iDim]/normal_length;
        }
        
        /*--- find the pressure at endpoint1 of this segment by interpolating
               the pressures at the two nodes of the edge it sits on ---*/
        double pressure_endpoint1_p1, pressure_endpoint1_p2;
        string input_file = "surface_flow.csv";
        unsigned short column = 6;
        char convert_me[32];
        string entry;
        
        // for endpoint1, find the pressure at poi.nt1
        //sprintf(convert_me,"%d",point1[iSection].at(iSegment));
        //entry = string(convert_me);
        //pressure_endpoint1_p1 = strtod(file_find(input_file, column, entry),NULL);
        
        // for endpoint1, find the pressure at point2
        //pressure_endpoint1_p2 = strtod(file_find(input_file, column, point1[iSection].at(iSegment+1)),NULL);
        
        
               
               
               // for endpoint1, compute weight2
               
               // using the two weights, interpolate to find the pressure at endpoint1,
               // i.e. segments[iSegment].pressure_endpoint1
               
               // rinse and repeat for endpoint2
               
               //segments[iSegment].pressure_endpoint1 = (1.0/distance_total)*(distance_p2*pressure_endpoint1_p1 + distance_p1*pressure_endpoint1_p2);
              
              /*--- find the pressure at endpoint2 of this segment by interpolating
               the pressures at the two nodes of the edge it sits on ---*/
              // segments[iSegment].pressure_endpoint2 = (1.0/distance_total)*(distance_p2*pressure_endpoint2_p1 + distance_p1*pressure_endpoint2_p2);
              
              /*--- find the value of pressure at the segment's midpoint
               by averaging the values found at the endpoints ---*/
              //double midPressure = 0.5*(segments[iSegment].pressure_endpoint1 + segments[iSegment].pressure_endpoint2);
              
              /*--- find the force per unit length carried by the segment ---*/
              //double midForce = midPressure*segments[iSegment].length;
              
              /*--- create a vector associated with the force, pointing
               in the direction of the inward-facing unit-normal ---*/
              //double force_vector [3];
              //for (iDim = 0; iDim < nDim; iDim++) {
              //  force_vector[iDim] = midForce*(-1.0*segments[iSegment].normal[iDim]);
              //}
      
     
            } // end loop over the segments
    
    
    
    
          
          
          
    } /*--- end loop over the sections ---*/
        
  } /*--- end loop over the markers ---*/
    
    
  
  /*--- delete the dynamically allocated memory  ---*/
  if (Xcoord_Airfoil != NULL) delete [] Xcoord_Airfoil;
  if (Ycoord_Airfoil != NULL) delete [] Ycoord_Airfoil;
  if (Zcoord_Airfoil != NULL) delete [] Zcoord_Airfoil;
  if (point1_Airfoil != NULL) delete [] point1_Airfoil;
  if (point2_Airfoil != NULL) delete [] point2_Airfoil;
  if (weight1_Airfoil != NULL) delete [] weight1_Airfoil;
  if (normal_force != NULL) delete [] normal_force;
  if (chord_force != NULL) delete [] chord_force;
  if (pitching_moment != NULL) delete [] pitching_moment;
  if (cl != NULL) delete [] cl;
  if (cd != NULL) delete [] cd;
  if (cm != NULL) delete [] cm;
  
  return 1;
}

string file_find(string input_file, unsigned short column, string entry) {
  
  /*--- local variables ---*/
  fstream file;
  string found_value, line;
  char tmp[400];
  char *pointer_to;
  unsigned short counter;
  
  /*--- initializations ---*/
  found_value = "";
  
  /*--- check to make sure the file exists ---*/
  if (fexists(input_file.c_str()) == false) {
    cout << endl;
    cout << "*** ERROR! \"" << input_file << "\" does not exist! ***" << endl;
    cout << "*** Please run mode 1 before running mode 2 ***" << endl;
    found_value = "file not found!";
    return found_value;
  }
  
  /*--- open the file ---*/
  file.open(input_file.c_str());
  
  /*--- check if it opened ---*/
  if (file.is_open()) {
  
    /*--- start reading lines ---*/
    while (getline(file,line) && (found_value != "")) {
      
      /*--- turn the line into a char array ---*/
      strcpy(tmp,line.c_str());
      
      /*--- start reading the line ---*/
      pointer_to = strtok(tmp," ,\r\n");
      
      /*--- does it match what we want? ---*/
      if (pointer_to == entry) {
        
        /*--- read the rest of the line ---*/
        counter = 1;
        while (pointer_to != NULL) {
          
          /*--- if it's the desired column ---*/
          if (counter == column) {
          
            /*--- record it ---*/
            found_value = pointer_to;
          }
        }
        
        /*--- if no vlaue assigned, print error ---*/
        if (found_value == "") {
          cout << "SOMETHING WENT WRONG! ";
          cout << "ENTRY FOUND BUT VALUE NOT ASSIGNED!" << endl;
        }
      }
    }
    
    /*--- entry not found, print error ---*/
    if (found_value == "") {
      cout << "SOMETHING WENT WRONG! ";
      cout << "DESIRED ENTRY NOT FOUND!" << endl;
    }
  }
  else
    cout << "Oh no! Unable to open the \"" << input_file << "\" file!" << endl;
  file.close();
  
  return found_value;
}
