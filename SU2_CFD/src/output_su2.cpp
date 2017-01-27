/*!
 * \file output_su2.cpp
 * \brief Main subroutines for output solver information.
 * \author F. Palacios, T. Economon, M. Colonno
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
  
  char cstr[MAX_STRING_SIZE], out_file[MAX_STRING_SIZE];
  unsigned long iElem, iPoint, iElem_Bound, nElem_Bound_, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], iNode, nElem;
  unsigned short iMarker, iDim, nDim = geometry->GetnDim(), iChar, iPeriodic, nPeriodic = 0, VTK_Type, nMarker_;
  su2double *center, *angles, *transl;
  ofstream output_file;
  ifstream input_file;
  string Grid_Marker, text_line, Marker_Tag, str;
  string::size_type position;

  /*--- Read the name of the output and input file ---*/
  
  str = config->GetMesh_Out_FileName();
  strcpy (out_file, str.c_str());
  strcpy (cstr, out_file);
  output_file.precision(15);
  output_file.open(cstr, ios::out);
  
  /*--- Write dimensions data. ---*/

  output_file << "NDIME= " << nDim << endl;
  
  /*--- Write the angle of attack offset. ---*/
  
  output_file << "AOA_OFFSET= " << config->GetAoA_Offset() << endl;
  
  /*--- Write the angle of attack offset. ---*/
  
  output_file << "AOS_OFFSET= " << config->GetAoS_Offset() << endl;

  /*--- Write connectivity data. ---*/
  
  nElem = nGlobal_Tria+nGlobal_Quad+nGlobal_Tetr+nGlobal_Hexa+nGlobal_Pris+nGlobal_Pyra;
  
  output_file << "NELEM= " << nElem<< endl;
  
  nElem = 0;
  
  for (iElem = 0; iElem < nGlobal_Tria; iElem++) {
    iNode = iElem*N_POINTS_TRIANGLE;
    output_file << "5\t";
    output_file << Conn_Tria[iNode+0]-1 << "\t"; output_file << Conn_Tria[iNode+1]-1 << "\t";
    output_file << Conn_Tria[iNode+2]-1 << "\t";
    output_file << nElem << "\n"; nElem++;
  }
  
  for (iElem = 0; iElem < nGlobal_Quad; iElem++) {
    iNode = iElem*N_POINTS_QUADRILATERAL;
    output_file << "9\t";
    output_file << Conn_Quad[iNode+0]-1 << "\t"; output_file << Conn_Quad[iNode+1]-1 << "\t";
    output_file << Conn_Quad[iNode+2]-1 << "\t"; output_file << Conn_Quad[iNode+3]-1 << "\t";
    output_file << nElem << "\n"; nElem++;
  }
  
  for (iElem = 0; iElem < nGlobal_Tetr; iElem++) {
    iNode = iElem*N_POINTS_TETRAHEDRON;
    output_file << "10\t";
    output_file << Conn_Tetr[iNode+0]-1 << "\t" << Conn_Tetr[iNode+1]-1 << "\t";
    output_file << Conn_Tetr[iNode+2]-1 << "\t" << Conn_Tetr[iNode+3]-1 << "\t";
    output_file << nElem << "\n"; nElem++;
  }
  
  for (iElem = 0; iElem < nGlobal_Hexa; iElem++) {
    iNode = iElem*N_POINTS_HEXAHEDRON;
    output_file << "12\t";
    output_file << Conn_Hexa[iNode+0]-1 << "\t" << Conn_Hexa[iNode+1]-1 << "\t";
    output_file << Conn_Hexa[iNode+2]-1 << "\t" << Conn_Hexa[iNode+3]-1 << "\t";
    output_file << Conn_Hexa[iNode+4]-1 << "\t" << Conn_Hexa[iNode+5]-1 << "\t";
    output_file << Conn_Hexa[iNode+6]-1 << "\t" << Conn_Hexa[iNode+7]-1 << "\t";
    output_file << nElem << "\n"; nElem++;
  }
  
  for (iElem = 0; iElem < nGlobal_Pris; iElem++) {
    iNode = iElem*N_POINTS_PRISM;
    output_file << "13\t";
    output_file << Conn_Pris[iNode+0]-1 << "\t" << Conn_Pris[iNode+1]-1 << "\t";
    output_file << Conn_Pris[iNode+2]-1 << "\t" << Conn_Pris[iNode+3]-1 << "\t";
    output_file << Conn_Pris[iNode+4]-1 << "\t" << Conn_Pris[iNode+5]-1 << "\t";
    output_file << nElem << "\n"; nElem++;
  }
  
  for (iElem = 0; iElem < nGlobal_Pyra; iElem++) {
    iNode = iElem*N_POINTS_PYRAMID;
    output_file << "14\t";
    output_file << Conn_Pyra[iNode+0]-1 << "\t" << Conn_Pyra[iNode+1]-1 << "\t";
    output_file << Conn_Pyra[iNode+2]-1 << "\t" << Conn_Pyra[iNode+3]-1 << "\t";
    output_file << Conn_Pyra[iNode+4]-1 << "\t";
    output_file << nElem << "\n"; nElem++;
  }
  
  /*--- Write the node coordinates ---*/
  
  output_file << "NPOIN= " << nGlobal_Doma;
  if (geometry->GetGlobal_nPointDomain() != nGlobal_Doma)
    output_file << "\t" << geometry->GetGlobal_nPointDomain();
  output_file << endl;

  for (iPoint = 0; iPoint < nGlobal_Doma; iPoint++) {
     for (iDim = 0; iDim < nDim; iDim++)
      output_file << scientific << Coords[iDim][iPoint] << "\t";
    output_file << iPoint << endl;
  }
  
  /*--- Read the boundary information ---*/
  
  input_file.open("boundary.su2", ios::out);
  
  /*--- Read grid file with format SU2 ---*/
  
  while (getline (input_file, text_line)) {

    /*--- Write the physical boundaries ---*/
    
    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      
      text_line.erase (0,6); nMarker_ = atoi(text_line.c_str());
      output_file << "NMARK= " << nMarker_ << endl;
      
      for (iMarker = 0 ; iMarker < nMarker_; iMarker++) {
        
        getline (input_file, text_line);
        text_line.erase (0,11);
        string::size_type position;
        for (iChar = 0; iChar < 20; iChar++) {
          position = text_line.find( " ", 0 );
          if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\r", 0 );
          if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\n", 0 );
          if (position != string::npos) text_line.erase (position,1);
        }
        Marker_Tag = text_line.c_str();
        
        /*--- Standart physical boundary ---*/
        
          getline (input_file, text_line);
          
          text_line.erase (0,13); nElem_Bound_ = atoi(text_line.c_str());
          output_file << "MARKER_TAG= " << Marker_Tag << endl;
          output_file << "MARKER_ELEMS= " << nElem_Bound_<< endl;
          
          if (Marker_Tag == "SEND_RECEIVE") {
            if (config->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
            if (config->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
          }
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
              case QUADRILATERAL:
                bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                output_file << "\t" << vnodes_quad[0] << "\t" << vnodes_quad[1] << "\t" << vnodes_quad[2] << "\t" << vnodes_quad[3] << endl;
                break;
              case VERTEX:
                bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                output_file << "\t" << vnodes_edge[0] <<  "\t" << vnodes_edge[1] <<endl;
                break;
            }
          }
      }
    }
    
  }
  
  input_file.close();

  remove("boundary.su2");

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
  
  output_file.close();
  
}

void COutput::SetSU2_MeshBinary(CConfig *config, CGeometry *geometry) { }
