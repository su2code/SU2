/*!
 * \file output_su2.cpp
 * \brief Main subroutines for output solver information.
 * \author F. Palacios, T. Economon, M. Colonno
 * \version 3.2.7.2 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (fpalacios@stanford.edu).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/output_structure.hpp"

void COutput::SetSU2_MeshASCII(CConfig *config, CGeometry *geometry) {
  
  ofstream SU2_File;
  char cstr[MAX_STRING_SIZE], out_file[MAX_STRING_SIZE], in_file[MAX_STRING_SIZE];
  string str;
  unsigned long iElem, iPoint, iElem_Bound, nElem_, nElem_Bound_, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], vnodes_tetra[4], vnodes_hexa[8], vnodes_wedge[6], vnodes_pyramid[5];
  unsigned short iMarker, iDim, nDim = geometry->GetnDim(), iChar, iPeriodic, nPeriodic = 0, VTK_Type, nDim_, nMarker_;
  double *center, *angles, *transl;
  ofstream output_file;
  ifstream input_file;
  string Grid_Marker, text_line, Marker_Tag;
  string::size_type position;

  /*--- Read the name of the output and input file ---*/
  
  str = config->GetMesh_Out_FileName();
  strcpy (out_file, str.c_str());
  strcpy (cstr, out_file);
  output_file.precision(15);
  output_file.open(cstr, ios::out);

  str = config->GetMesh_FileName();
  strcpy (in_file, str.c_str());
  strcpy (cstr, in_file);
  input_file.open(cstr, ios::out);
  
  /*--- Read grid file with format SU2 ---*/
  
  while (getline (input_file, text_line)) {
    
    /*--- Read the dimension of the problem ---*/
    
    position = text_line.find ("NDIME=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nDim_ = atoi(text_line.c_str());
      output_file << "NDIME= " << nDim_ << endl;
    }
    
    /*--- Read the information about inner elements ---*/
    
    position = text_line.find ("NELEM=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nElem_ = atoi(text_line.c_str());
      output_file << "NELEM= " << nElem_ << endl;
      
      /*--- Loop over all the volumetric elements ---*/
      
      for (iElem = 0; iElem < nElem_;  iElem++) {
        getline(input_file, text_line);
        istringstream elem_line(text_line);
        
        elem_line >> VTK_Type;
        output_file << VTK_Type;
        
        switch(VTK_Type) {
          case TRIANGLE:
            elem_line >> vnodes_triangle[0]; elem_line >> vnodes_triangle[1]; elem_line >> vnodes_triangle[2];
            output_file << "\t" << vnodes_triangle[0] << "\t" << vnodes_triangle[1] << "\t" << vnodes_triangle[2] << endl;
            break;
          case RECTANGLE:
            elem_line >> vnodes_quad[0]; elem_line >> vnodes_quad[1]; elem_line >> vnodes_quad[2]; elem_line >> vnodes_quad[3];
            output_file << "\t" << vnodes_quad[0] << "\t" << vnodes_quad[1] << "\t" << vnodes_quad[2] << "\t" << vnodes_quad[3] << endl;
            break;
          case TETRAHEDRON:
            elem_line >> vnodes_tetra[0]; elem_line >> vnodes_tetra[1]; elem_line >> vnodes_tetra[2]; elem_line >> vnodes_tetra[3];
            output_file << "\t" << vnodes_tetra[0] << "\t" << vnodes_tetra[1] << "\t" << vnodes_tetra[2] << "\t" << vnodes_tetra[3] << endl;
            break;
          case HEXAHEDRON:
            elem_line >> vnodes_hexa[0]; elem_line >> vnodes_hexa[1]; elem_line >> vnodes_hexa[2];
            elem_line >> vnodes_hexa[3]; elem_line >> vnodes_hexa[4]; elem_line >> vnodes_hexa[5];
            elem_line >> vnodes_hexa[6]; elem_line >> vnodes_hexa[7];
            output_file << "\t" << vnodes_hexa[0] << "\t" << vnodes_hexa[1] << "\t" << vnodes_hexa[2] << "\t" << vnodes_hexa[3] << "\t" << vnodes_hexa[4] << "\t" << vnodes_hexa[5] << "\t" << vnodes_hexa[6] << "\t" << vnodes_hexa[7] << endl;
            break;
          case WEDGE:
            elem_line >> vnodes_wedge[0]; elem_line >> vnodes_wedge[1]; elem_line >> vnodes_wedge[2];
            elem_line >> vnodes_wedge[3]; elem_line >> vnodes_wedge[4]; elem_line >> vnodes_wedge[5];
            output_file << "\t" << vnodes_wedge[0] << "\t" << vnodes_wedge[1] << "\t" << vnodes_wedge[2] << "\t" << vnodes_wedge[3] << "\t" << vnodes_wedge[4] << "\t" << vnodes_wedge[5] << endl;
            break;
          case PYRAMID:
            elem_line >> vnodes_pyramid[0]; elem_line >> vnodes_pyramid[1]; elem_line >> vnodes_pyramid[2];
            elem_line >> vnodes_pyramid[3]; elem_line >> vnodes_pyramid[4];
            output_file << "\t" << vnodes_pyramid[0] << "\t" << vnodes_pyramid[1] << "\t" << vnodes_pyramid[2] << "\t" << vnodes_pyramid[3] << "\t" << vnodes_pyramid[4] << endl;
            break;
        }
      }
    }
    
    /*--- Coordinates ---*/
    
    position = text_line.find ("NPOIN=",0);
    if (position != string::npos) {

      /*--- Skip the lines about the points ---*/
      
      for (iPoint = 0; iPoint < nGlobal_Poin;  iPoint++) {
        getline(input_file, text_line);
      }
      
      /*--- Add the new coordinates ---*/
      
      output_file << "NPOIN= " << nGlobal_Poin << endl;
      
      /*--- Write surface and volumetric solution data. ---*/
      
      for (iPoint = 0; iPoint < nGlobal_Poin; iPoint++) {
        
        /*--- Write the node coordinates ---*/
        
        for(iDim = 0; iDim < nDim; iDim++)
          output_file << scientific << Coords[iDim][iPoint] << "\t";
        
        output_file << iPoint << endl;
        
      }
      
    }
    
    /*--- Write the physical boundaries ---*/
    
    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      
      text_line.erase (0,6); nMarker_ = atoi(text_line.c_str());
      output_file << "NMARK= " << nMarker_ << endl;
      
      for (iMarker = 0 ; iMarker < nMarker_; iMarker++) {
        
        getline (input_file,text_line);
        text_line.erase (0,11);
        string::size_type position;
        for (iChar = 0; iChar < 20; iChar++) {
          position = text_line.find( " ", 0 );
          if(position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\r", 0 );
          if(position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\n", 0 );
          if(position != string::npos) text_line.erase (position,1);
        }
        Marker_Tag = text_line.c_str();
        
        /*--- Standart physical boundary ---*/
        
          getline (input_file, text_line);
          
          text_line.erase (0,13); nElem_Bound_ = atoi(text_line.c_str());
          output_file << "MARKER_TAG= " << Marker_Tag << endl;
          output_file << "MARKER_ELEMS= " << nElem_Bound_<< endl;
          
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
            
            getline(input_file, text_line);
            istringstream bound_line(text_line);
            
            bound_line >> VTK_Type;
            output_file << VTK_Type;
            
            switch(VTK_Type) {
              case LINE:
                bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                output_file << "\t" << vnodes_edge[0] << "\t" << vnodes_edge[1] << endl;
                break;
              case TRIANGLE:
                bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                output_file << "\t" << vnodes_triangle[0] << "\t" << vnodes_triangle[1] << "\t" << vnodes_triangle[2] << endl;
                break;
              case RECTANGLE:
                bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                output_file << "\t" << vnodes_quad[0] << "\t" << vnodes_quad[1] << "\t" << vnodes_quad[2] << "\t" << vnodes_quad[3] << endl;
                break;
            }
          }
      }
    }
  }
  
  
  /*--- Get the total number of periodic transformations ---*/
  
  nPeriodic = config->GetnPeriodicIndex();
  output_file << "NPERIODIC= " << nPeriodic << endl;
  
  /*--- From iPeriodic obtain the iMarker ---*/
  
  for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
    
    /*--- Retrieve the supplied periodic information. ---*/
    
    center = config->GetPeriodicCenter(iPeriodic);
    angles = config->GetPeriodicRotation(iPeriodic);
    transl = config->GetPeriodicTranslate(iPeriodic);
    
    output_file << "PERIODIC_INDEX= " << iPeriodic << endl;
    output_file << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
    output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << endl;
    output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << endl;
    
  }
  
  input_file.close();
  output_file.close();
  
}

void COutput::SetSU2_MeshBinary(CConfig *config, CGeometry *geometry) { }
