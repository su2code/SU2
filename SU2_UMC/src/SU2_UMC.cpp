/*!
 * \file SU2_UMC.cpp
 * \brief Main file for UMarc Coupling (SU2_UMC).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.8
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
                         VERB_HIGH);

  /*--- Instantiate an object of the boundary-geometry class ---*/
  boundary = new CBoundaryGeometry(config, 
                                   config->GetMesh_FileName(), 
                                   config->GetMesh_FileFormat());
  
  /*--- determine the mode of operation ---*/
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
      
      /*--- Make the sectional cuts along the surfaces of interest and
            write the files "span," "tracked points," and 
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
	if (config != NULL) {delete config;}
	if (boundary != NULL) {delete boundary;}
	
	/*--- End routine ---*/
	cout << endl <<"------------------------- Exit Success (SU2_UMC) ------------------------" << endl << endl;
  
	return EXIT_SUCCESS;
	
}

void Make_Sectional_Cuts(CConfig *config, CBoundaryGeometry *boundary) {
  
  /*--- Local variables ---*/
  unsigned short nDim, nMarker, nSection, iDim, iMarker, iSection;
  unsigned long iNode, nNode;
  vector<unsigned long> **point1_Airfoil, **point2_Airfoil;
  double root, tip, minPlane, maxPlane, radial_station, interval,
         starting_angle, starting_rad, rotation, magnitude;
  double starting_vec [3]; 
  double *Plane_P0, *Plane_Normal;
  vector<double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, 
                 **weight1_Airfoil;
  bool Boundary, Monitoring;
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
	        if (Marker_Tag == "WING") {starting_angle = 90.0;}
          
	        /*--- UH-60 FIRST blade ---*/
	        if (Marker_Tag == "Blade1") {starting_angle = 0.0;}
          
	        /*--- UH-60 SECOND blade ---*/
	        if (Marker_Tag == "Blade2") {starting_angle = 90.0;}
          
	        /*--- UH-60 THIRD blade ---*/
	        if (Marker_Tag == "Blade3") {starting_angle = 180.0;}
          
	        /*--- UH-60 FOURTH blade ---*/
	        if (Marker_Tag == "Blade4") {starting_angle = 270.0;}
	        
	        /*--- Caradonna-Tung FIRST blade ---*/
          if (Marker_Tag == "blade_1") {
            starting_angle = 90.0;
            root = 0.15;
            tip = 1.14;
            minPlane = root;
            maxPlane = tip;
          }
          
	        /*--- Caradonna-Tung SECOND blade ---*/
	        if (Marker_Tag == "blade_2") {
            starting_angle = 270.0;
            root = 0.15;
            tip = 1.14;
            minPlane = root;
            maxPlane = tip;
          }
	        
	        /*--- USER INPUTS ABOUT GEOMETRY END HERE ---*/
          /*--- comupte the distance between sections ---*/
          interval = (maxPlane-minPlane)/double(nSection-1);
          
          /*--- write the span file, listing the radial 
                stations along the geometry being sliced ---*/
          span.precision(15);
          span.open("span", ios::out);
          for (iSection = nSection; iSection > 0; iSection--) {
            radial_station = minPlane + (iSection-1)*interval;
            span << radial_station/(tip-root) << endl;
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
            
            /*--- find the points defining the airfoil section --*/
	          boundary->ComputeAirfoil_Section(Plane_P0, Plane_Normal,
	                                          iSection, config,
                                            Xcoord_Airfoil[iSection],
                                            Ycoord_Airfoil[iSection],
                                            Zcoord_Airfoil[iSection],
                                   point1_Airfoil[iMarker][iSection],
                                   point2_Airfoil[iMarker][iSection],
                                  weight1_Airfoil[iMarker][iSection],
                                                               true);
            
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
            
	          /*--- print information to file for a given node 
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
      cout << "Unable to open coordinates.csv file!" << endl;
    if (coords_file.is_open())
      cout << "Unable to open tracked points file!" << endl;
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

