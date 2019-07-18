/*!
 * \file output_su2.cpp
 * \brief Main subroutines for output solver information.
 * \author F. Palacios, T. Economon, M. Colonno
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

void COutput::SetSU2_MeshASCII(CConfig *config, CGeometry *geometry, unsigned short val_iZone, ofstream& output_file) {
  
  unsigned long iElem, iPoint, iElem_Bound, nElem_Bound_, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], iNode, nElem;
  unsigned short iMarker, iDim, nDim = geometry->GetnDim(), iChar, VTK_Type, nMarker_;
  short SendTo;
  ifstream input_file;
  string Grid_Marker, text_line, Marker_Tag, str;
  string::size_type position;

  if (config->GetnZone() > 1){
    output_file << "IZONE= " << val_iZone+1 << endl;
  }

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

  str = "boundary.dat";

  str = config->GetMultizone_FileName(str, val_iZone);

  input_file.open(str.c_str(), ios::out);
  
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
          getline (input_file, text_line);

          text_line.erase (0,8); SendTo = atoi(text_line.c_str());

          if (Marker_Tag == "SEND_RECEIVE"){
            output_file << "SEND_TO= " << SendTo << endl;
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

  remove(str.c_str());

}

void COutput::SetSU2_MeshBinary(CConfig *config, CGeometry *geometry) { }

void COutput::WriteCoordinates_Binary(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;
  char cstr[200];
  
  /*--- Prepare the file name. ---*/
  
  strcpy(cstr, "coordinates");
  if (config->GetnZone() > 1){
    char appstr[200];
    SPRINTF(appstr, "_%u", val_iZone);
    strcat(cstr, appstr);
  }
  strcat(cstr,".dat");
  
  /*--- Prepare the first ints containing the counts. The first is a
   the total number of points written. The second is the dimension.
   We know the rest of the file will contain the coords (nPoints*nDim). ---*/
  
  int var_buf_size = 2;
  int var_buf[2]   = {(int)nGlobal_Poin, nDim};
  
  /*--- Prepare the 1D data buffer on this rank. ---*/
  
  passivedouble *buf = new passivedouble[nGlobal_Poin*nDim];
  
  /*--- For now, create a temp 1D buffer to load up the data for writing.
   This will be replaced with a derived data type most likely. ---*/
  
  for (iPoint = 0; iPoint < nGlobal_Poin; iPoint++)
    for (iDim = 0; iDim < nDim; iDim++)
      buf[iPoint*nDim+iDim] = SU2_TYPE::GetValue(Coords[iDim][iPoint]);
  
  /*--- Write the binary file. Everything is done in serial, as only the
   master node has the data and enters this routine. ---*/
  
  FILE* fhw;
  fhw = fopen(cstr, "wb");
  
  /*--- Error check for opening the file. ---*/
  
  if (!fhw) {
    SU2_MPI::Error(string("Unable to open binary coordinates file ") +
                   string(cstr), CURRENT_FUNCTION);
  }
  
  /*--- First, write the number of variables and points. ---*/
  
  fwrite(var_buf, var_buf_size, sizeof(int), fhw);
  
  /*--- Call to write the entire restart file data in binary in one shot. ---*/
  
  fwrite(buf, nGlobal_Poin*nDim, sizeof(passivedouble), fhw);
  
  /*--- Close the file. ---*/
  
  fclose(fhw);
  
  /*--- Release buffer memory. ---*/
  
  delete [] buf;
  
}

void COutput::WriteProjectedSensitivity(CConfig *config,
                                        CGeometry *geometry,
                                        unsigned short val_iZone,
                                        unsigned short val_nZone) {
  
  unsigned short iVar;
  unsigned long iPoint, iElem, iNode;
  unsigned long *LocalIndex = NULL;
  bool *SurfacePoint = NULL;
  string filename, fieldname;
  
  filename = config->GetDV_Sens_Filename();
  
  ofstream Sens_File;
  Sens_File.open(filename.c_str(), ios::out);
  Sens_File.precision(15);
  
  /*--- This is surface output, so print only the points
   that are in the element list. Change the numbering. ---*/
  
  LocalIndex   = new unsigned long [nGlobal_Poin+1];
  SurfacePoint = new bool [nGlobal_Poin+1];
  
  for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++)
    SurfacePoint[iPoint] = false;
  
  for (iElem = 0; iElem < nGlobal_Line; iElem++) {
    iNode = iElem*N_POINTS_LINE;
    SurfacePoint[Conn_Line[iNode+0]] = true;
    SurfacePoint[Conn_Line[iNode+1]] = true;
  }
  for (iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
    iNode = iElem*N_POINTS_TRIANGLE;
    SurfacePoint[Conn_BoundTria[iNode+0]] = true;
    SurfacePoint[Conn_BoundTria[iNode+1]] = true;
    SurfacePoint[Conn_BoundTria[iNode+2]] = true;
  }
  for (iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
    iNode = iElem*N_POINTS_QUADRILATERAL;
    SurfacePoint[Conn_BoundQuad[iNode+0]] = true;
    SurfacePoint[Conn_BoundQuad[iNode+1]] = true;
    SurfacePoint[Conn_BoundQuad[iNode+2]] = true;
    SurfacePoint[Conn_BoundQuad[iNode+3]] = true;
  }
  
  nSurf_Poin = 0;
  for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) {
    LocalIndex[iPoint] = 0;
    if (SurfacePoint[iPoint]) {nSurf_Poin++; LocalIndex[iPoint] = nSurf_Poin;}
  }
  
  /*--- Write surface x,y,z and surface dJ/dx, dJ/dy, dJ/dz data. ---*/
  
  for (iPoint = 0; iPoint < nGlobal_Poin; iPoint++) {
    
    if (LocalIndex[iPoint+1] != 0) {
      
      /*--- Write the node coordinates and the sensitivities. Note that
       we subtract 2 from the fields list to ignore the initial ID and
       final sens.normal value in the Data array. ---*/
      
      for (iVar = 0; iVar < config->fields.size()-2; iVar++)
      Sens_File << scientific << Data[iVar][iPoint] << "\t";
      Sens_File << scientific << "\n";
      
    }
  }
  
}

