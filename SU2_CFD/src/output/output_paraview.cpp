/*!
 * \file output_paraview.cpp
 * \brief Main subroutines for the output of ParaView visualization files.
 * \author F. Palacios, T. Economon, E. van der Weide
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

#include "../../include/output/COutput.hpp"

/*--- Subroutine to swap bytes, in case we need to convert to
 big endian, which is expected for ParaView binary legacy format. ---*/

void SwapBytes(char *buffer,
               size_t nBytes,
               unsigned long nVar) {
  
  /*--- Store half the number of bytes in kk. ---*/
  
  const int kk = (int)nBytes/2;
  
  /*--- Loop over the number of variables in the buffer. ---*/
  
  for (int j = 0; j < (int)nVar; j++) {
    
    /*--- Initialize ii and jj, which are used to store the
     indices of the bytes to be swapped. ---*/
    
    int ii = j*(int)nBytes;
    int jj = ii + (int)nBytes - 1;
    
    /*--- Swap the bytes. ---*/
    
    for (int i = 0; i < kk; i++) {
      char tmp   = buffer[jj];
      buffer[jj] = buffer[ii];
      buffer[ii] = tmp;
      
      ii++;
      jj--;
      
    }
  }
}

void COutput::WriteParaViewASCII_Parallel(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol) {

  unsigned short iDim, nDim = geometry->GetnDim();

  unsigned long iPoint, iElem, iNode;

  unsigned long nSurf_Elem_Storage;
  unsigned long nGlobal_Elem_Storage;

  string filename, fieldname;
  ofstream Paraview_File;

  int iProcessor;
  
  if (surf_sol) filename = config->GetFilename(SurfaceFilename, ".vtk");
  else filename          = config->GetFilename(VolumeFilename, ".vtk");


  /*--- Open Paraview ASCII file and write the header. ---*/

    if (rank == MASTER_NODE) {
  Paraview_File.open(filename.c_str(), ios::out);
  Paraview_File.precision(6);
  Paraview_File << "# vtk DataFile Version 3.0\n";
  Paraview_File << "vtk output\n";
  Paraview_File << "ASCII\n";
  Paraview_File << "DATASET UNSTRUCTURED_GRID\n";

  /*--- Write the header ---*/
  if (surf_sol) Paraview_File << "POINTS "<< nGlobal_Surf_Poin <<" double\n";
  else Paraview_File << "POINTS "<< nGlobal_Poin_Par <<" double\n";

    }

  Paraview_File.close();

#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  /*--- Each processor opens the file. ---*/

  Paraview_File.open(filename.c_str(), ios::out | ios::app);

  /*--- Write surface and volumetric point coordinates. ---*/

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {

      /*--- Write the node data from this proc ---*/

      if (surf_sol) {
        for (iPoint = 0; iPoint < nSurf_Poin_Par; iPoint++) {
          for (iDim = 0; iDim < nDim; iDim++)
            Paraview_File << scientific << Parallel_Surf_Data[iDim][iPoint] << "\t";
          if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
        }
      } else {

        for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
          for (iDim = 0; iDim < nDim; iDim++)
            Paraview_File << scientific << Parallel_Data[iDim][iPoint] << "\t";
          if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
        }
      }
    }
    Paraview_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }

  /*--- Reduce the total number of each element. ---*/

  unsigned long nTot_Line, nTot_BoundTria, nTot_BoundQuad, nTot_Tria, nTot_Quad, nTot_Tetr, nTot_Hexa, nTot_Pris, nTot_Pyra;
#ifdef HAVE_MPI
  SU2_MPI::Reduce(&nParallel_Line, &nTot_Line, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_BoundTria, &nTot_BoundTria, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_BoundQuad, &nTot_BoundQuad, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);

  SU2_MPI::Reduce(&nParallel_Tria, &nTot_Tria, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_Quad, &nTot_Quad, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_Tetr, &nTot_Tetr, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_Hexa, &nTot_Hexa, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_Pris, &nTot_Pris, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_Pyra, &nTot_Pyra, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
  nTot_Line      = nParallel_Line;
  nTot_BoundTria = nParallel_BoundTria;
  nTot_BoundQuad = nParallel_BoundQuad;

  nTot_Tria = nParallel_Tria;
  nTot_Quad = nParallel_Quad;
  nTot_Tetr = nParallel_Tetr;
  nTot_Hexa = nParallel_Hexa;
  nTot_Pris = nParallel_Pris;
  nTot_Pyra = nParallel_Pyra;
#endif

  if (rank == MASTER_NODE) {

  /*--- Write the header ---*/
  nSurf_Elem_Storage = nTot_Line*3 +nTot_BoundTria*4 + nTot_BoundQuad*5;
  nGlobal_Elem_Storage = nTot_Tria*4 + nTot_Quad*5 + nTot_Tetr*5 + nTot_Hexa*9 + nTot_Pris*7 + nTot_Pyra*6;

  if (surf_sol) Paraview_File << "\nCELLS " << nSurf_Elem_Par << "\t" << nSurf_Elem_Storage << "\n";
  else Paraview_File << "\nCELLS " << nGlobal_Elem_Par << "\t" << nGlobal_Elem_Storage << "\n";

  }

  Paraview_File.flush();
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  /*--- Write connectivity data. ---*/

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {

  if (surf_sol) {

    for (iElem = 0; iElem < nParallel_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      Paraview_File << N_POINTS_LINE << "\t";
      Paraview_File << Conn_BoundLine_Par[iNode+0]-1 << "\t";
      Paraview_File << Conn_BoundLine_Par[iNode+1]-1 << "\t";
    }

    for (iElem = 0; iElem < nParallel_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      Paraview_File << N_POINTS_TRIANGLE << "\t";
      Paraview_File << Conn_BoundTria_Par[iNode+0]-1 << "\t";
      Paraview_File << Conn_BoundTria_Par[iNode+1]-1 << "\t";
      Paraview_File << Conn_BoundTria_Par[iNode+2]-1 << "\t";
    }

    for (iElem = 0; iElem < nParallel_BoundQuad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      Paraview_File << N_POINTS_QUADRILATERAL << "\t";
      Paraview_File << Conn_BoundQuad_Par[iNode+0]-1 << "\t";
      Paraview_File << Conn_BoundQuad_Par[iNode+1]-1 << "\t";
      Paraview_File << Conn_BoundQuad_Par[iNode+2]-1 << "\t";
      Paraview_File << Conn_BoundQuad_Par[iNode+3]-1 << "\t";
    }

  }
  else {

    for (iElem = 0; iElem < nParallel_Tria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      Paraview_File << N_POINTS_TRIANGLE << "\t";
      Paraview_File << Conn_Tria_Par[iNode+0]-1 << "\t";
      Paraview_File << Conn_Tria_Par[iNode+1]-1 << "\t";
      Paraview_File << Conn_Tria_Par[iNode+2]-1 << "\t";
    }

    for (iElem = 0; iElem < nParallel_Quad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      Paraview_File << N_POINTS_QUADRILATERAL << "\t";
      Paraview_File << Conn_Quad_Par[iNode+0]-1 << "\t";
      Paraview_File << Conn_Quad_Par[iNode+1]-1 << "\t";
      Paraview_File << Conn_Quad_Par[iNode+2]-1 << "\t";
      Paraview_File << Conn_Quad_Par[iNode+3]-1 << "\t";
    }

    for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
      iNode = iElem*N_POINTS_TETRAHEDRON;
      Paraview_File << N_POINTS_TETRAHEDRON << "\t";
      Paraview_File << Conn_Tetr_Par[iNode+0]-1 << "\t" << Conn_Tetr_Par[iNode+1]-1 << "\t";
      Paraview_File << Conn_Tetr_Par[iNode+2]-1 << "\t" << Conn_Tetr_Par[iNode+3]-1 << "\t";
    }

    for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
      iNode = iElem*N_POINTS_HEXAHEDRON;
      Paraview_File << N_POINTS_HEXAHEDRON << "\t";
      Paraview_File << Conn_Hexa_Par[iNode+0]-1 << "\t" << Conn_Hexa_Par[iNode+1]-1 << "\t";
      Paraview_File << Conn_Hexa_Par[iNode+2]-1 << "\t" << Conn_Hexa_Par[iNode+3]-1 << "\t";
      Paraview_File << Conn_Hexa_Par[iNode+4]-1 << "\t" << Conn_Hexa_Par[iNode+5]-1 << "\t";
      Paraview_File << Conn_Hexa_Par[iNode+6]-1 << "\t" << Conn_Hexa_Par[iNode+7]-1 << "\t";
    }

    for (iElem = 0; iElem < nParallel_Pris; iElem++) {
      iNode = iElem*N_POINTS_PRISM;
      Paraview_File << N_POINTS_PRISM << "\t";
      Paraview_File << Conn_Pris_Par[iNode+0]-1 << "\t" << Conn_Pris_Par[iNode+1]-1 << "\t";
      Paraview_File << Conn_Pris_Par[iNode+2]-1 << "\t" << Conn_Pris_Par[iNode+3]-1 << "\t";
      Paraview_File << Conn_Pris_Par[iNode+4]-1 << "\t" << Conn_Pris_Par[iNode+5]-1 << "\t";
    }

    for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
      iNode = iElem*N_POINTS_PYRAMID;
      Paraview_File << N_POINTS_PYRAMID << "\t";
      Paraview_File << Conn_Pyra_Par[iNode+0]-1 << "\t" << Conn_Pyra_Par[iNode+1]-1 << "\t";
      Paraview_File << Conn_Pyra_Par[iNode+2]-1 << "\t" << Conn_Pyra_Par[iNode+3]-1 << "\t";
      Paraview_File << Conn_Pyra_Par[iNode+4]-1 << "\t";
    }
  }
    }    Paraview_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }

    if (rank == MASTER_NODE) {

  /*--- Write the header ---*/
  if (surf_sol) Paraview_File << "\nCELL_TYPES " << nSurf_Elem_Par << "\n";
  else Paraview_File << "\nCELL_TYPES " << nGlobal_Elem_Par << "\n";
    }

  Paraview_File.flush();
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      if (surf_sol) {
        for (iElem = 0; iElem < nParallel_Line; iElem++) Paraview_File << "3\t";
        for (iElem = 0; iElem < nParallel_BoundTria; iElem++) Paraview_File << "5\t";
        for (iElem = 0; iElem < nParallel_BoundQuad; iElem++) Paraview_File << "9\t";
      }
      else {
        for (iElem = 0; iElem < nParallel_Tria; iElem++) Paraview_File << "5\t";
        for (iElem = 0; iElem < nParallel_Quad; iElem++) Paraview_File << "9\t";
        for (iElem = 0; iElem < nParallel_Tetr; iElem++) Paraview_File << "10\t";
        for (iElem = 0; iElem < nParallel_Hexa; iElem++) Paraview_File << "12\t";
        for (iElem = 0; iElem < nParallel_Pris; iElem++) Paraview_File << "13\t";
        for (iElem = 0; iElem < nParallel_Pyra; iElem++) Paraview_File << "14\t";
      }
    }    Paraview_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }
  
    if (rank == MASTER_NODE) {
  /*--- Write the header ---*/
  if (surf_sol) Paraview_File << "\nPOINT_DATA "<< nGlobal_Surf_Poin <<"\n";
  else Paraview_File << "\nPOINT_DATA "<< nGlobal_Poin_Par <<"\n";

    }

  Paraview_File.flush();
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  unsigned short varStart = 2;
  if (nDim == 3) varStart++;

  /*--- Need to adjust container location to avoid PointID tag and coords. ---*/
  unsigned short VarCounter = varStart;

  for (unsigned short iField = varStart; iField < Variable_Names.size(); iField++) {

    fieldname = Variable_Names[iField];

    fieldname.erase(remove(fieldname.begin(), fieldname.end(), '"'), fieldname.end());

    bool output_variable = true, isVector = false;
    size_t found = Variable_Names[iField].find("_x");
    if (found!=string::npos) {
      output_variable = true;
      isVector = true;
    }
    found = Variable_Names[iField].find("_y");
    if (found!=string::npos) {
      output_variable = false;
      //skip
      Paraview_File.flush();
#ifdef HAVE_MPI
      SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
      VarCounter++;
    }
found = Variable_Names[iField].find("_z");
    if (found!=string::npos) {
      output_variable = false;
      //skip
      Paraview_File.flush();
#ifdef HAVE_MPI
      SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
      VarCounter++;
    }

    if (output_variable && isVector) {

      fieldname.erase(fieldname.end()-2,fieldname.end());

      if (rank == MASTER_NODE) {
        Paraview_File << "\nVECTORS " << fieldname << " double\n";
      }

      Paraview_File.flush();
#ifdef HAVE_MPI
      SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

      /*--- Write surface and volumetric point coordinates. ---*/

      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        if (rank == iProcessor) {

          /*--- Write the node data from this proc ---*/

          if (surf_sol) {
            for (iPoint = 0; iPoint < nSurf_Poin_Par; iPoint++) {
              /*--- Loop over the vars/residuals and write the values to file ---*/
              Paraview_File << scientific << Parallel_Surf_Data[VarCounter+0][iPoint] << "\t" << Parallel_Surf_Data[VarCounter+1][iPoint] << "\t";
              if (nDim == 3) Paraview_File << scientific << Parallel_Surf_Data[VarCounter+2][iPoint] << "\t";
              if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
            }
          } else {
            for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
              Paraview_File << scientific << Parallel_Data[VarCounter+0][iPoint] << "\t" << Parallel_Data[VarCounter+1][iPoint] << "\t";
              if (nDim == 3) Paraview_File << scientific << Parallel_Data[VarCounter+2][iPoint] << "\t";
              if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
            }
          }
        }
        Paraview_File.flush();
#ifdef HAVE_MPI
        SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
      }

      VarCounter++;

    } else if (output_variable) {

      if (rank == MASTER_NODE) {

        Paraview_File << "\nSCALARS " << fieldname << " double 1\n";
        Paraview_File << "LOOKUP_TABLE default\n";
      }

      Paraview_File.flush();
#ifdef HAVE_MPI
      SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

      /*--- Write surface and volumetric point coordinates. ---*/

      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        if (rank == iProcessor) {

          /*--- Write the node data from this proc ---*/

          if (surf_sol) {
            for (iPoint = 0; iPoint < nSurf_Poin_Par; iPoint++) {
              Paraview_File << scientific << Parallel_Surf_Data[VarCounter][iPoint] << "\t";
            }
          } else {
            for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
              Paraview_File << scientific << Parallel_Data[VarCounter][iPoint] << "\t";
            }
          }
        }
        Paraview_File.flush();
#ifdef HAVE_MPI
        SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
      }
      
      VarCounter++;
    }

  }

  Paraview_File.close();
  
}

void COutput::WriteParaViewBinary_Parallel(CConfig *config,
                                           CGeometry *geometry,
                                           unsigned short val_iZone,
                                           unsigned short val_nZone,
                                           bool surf_sol) {
  
  unsigned short iDim, nDim = geometry->GetnDim();
  
  unsigned long iPoint, iElem, iNode;
  
  string filename, fieldname;
  ofstream Paraview_File;
    
  int MAX_STRING_LENGTH = 255;
  char str_buf[MAX_STRING_LENGTH], fname[100];
  
  const int NCOORDS = 3;
  
  if (surf_sol) filename = config->GetFilename(SurfaceFilename, ".vtk");
  else filename          = config->GetFilename(VolumeFilename, ".vtk");

  strcpy(fname, filename.c_str());
  
  /*--- Check for big endian. We have to swap bytes otherwise. ---*/
  
  bool BigEndian;
  union {int i; char c[4];} val;
  val.i = 0x76543210;
  if (val.c[0] == 0x10) BigEndian = false;
  else BigEndian = true;
  
  /*--- Serial implementation in case we have not compiled with MPI. ---*/
  
#ifndef HAVE_MPI
  
  FILE* fhw;
  fhw = fopen(fname, "wb");
  
  unsigned long iNode2;
  unsigned long nSurf_Elem_Storage;
  unsigned long nGlobal_Elem_Storage;
  
  /*--- Error check for opening the file. ---*/
  
  if (!fhw) {
    SU2_MPI::Error(string("Unable to open VTK binary legacy file ") +
                   filename, CURRENT_FUNCTION);
  }
  
  /*--- File header written in ASCII. ---*/
  
  strcpy(str_buf, "# vtk DataFile Version 3.0\n");
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  
  strcpy(str_buf, "vtk output\n");
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  
  strcpy(str_buf, "BINARY\n");
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  
  strcpy(str_buf, "DATASET UNSTRUCTURED_GRID\n");
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  
  /*--- Write the point coordinates. ---*/
  
  unsigned long GlobalPoint;
  su2double **Data;
  if (surf_sol) {
    GlobalPoint  = nSurf_Poin_Par;
    Data         = Parallel_Surf_Data;
  } else {
    GlobalPoint  = nGlobal_Poin_Par;
    Data         = Parallel_Data;
  }
  
  SPRINTF(str_buf, "POINTS %i float\n", (int)GlobalPoint);
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  
  /*--- Load/write the 1D buffer of point coordinates. ---*/
  
  float *coord_buf = new float[GlobalPoint*NCOORDS];
  for (iPoint = 0; iPoint < GlobalPoint; iPoint++) {
    for (iDim = 0; iDim < NCOORDS; iDim++) {
      if (nDim == 2 && iDim == 2) {
        coord_buf[iPoint*NCOORDS + iDim] = 0.0;
      } else {
        float val = (float)SU2_TYPE::GetValue(Data[iDim][iPoint]);
        coord_buf[iPoint*NCOORDS + iDim] = val;
      }
    }
  }
  if (!BigEndian) SwapBytes((char *)coord_buf, sizeof(float), 3*GlobalPoint);
  
  fwrite(coord_buf, sizeof(float), 3*GlobalPoint, fhw);
  delete [] coord_buf;
  
  /*--- Write the connectivity data. ---*/
  
  unsigned long nTot_Line, nTot_BoundTria, nTot_BoundQuad;
  nTot_Line      = nParallel_Line;
  nTot_BoundTria = nParallel_BoundTria;
  nTot_BoundQuad = nParallel_BoundQuad;
  nSurf_Elem_Storage = nTot_Line*3 +nTot_BoundTria*4 + nTot_BoundQuad*5;
  
  unsigned long nTot_Tria, nTot_Quad;
  unsigned long nTot_Tetr, nTot_Hexa, nTot_Pris, nTot_Pyra;
  nTot_Tria = nParallel_Tria;
  nTot_Quad = nParallel_Quad;
  nTot_Tetr = nParallel_Tetr;
  nTot_Hexa = nParallel_Hexa;
  nTot_Pris = nParallel_Pris;
  nTot_Pyra = nParallel_Pyra;
  nGlobal_Elem_Storage = (nTot_Tria*4 + nTot_Quad*5 + nTot_Tetr*5 +
                          nTot_Hexa*9 + nTot_Pris*7 + nTot_Pyra*6);
  
  int *conn_buf = NULL;
  
  if (surf_sol) {
    SPRINTF (str_buf, "\nCELLS %i %i\n", (int)nSurf_Elem_Par,
             (int)nSurf_Elem_Storage);
    fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
    conn_buf = new int[nSurf_Elem_Par*(N_POINTS_QUADRILATERAL+1)];
  } else {
    SPRINTF (str_buf, "\nCELLS %i %i\n", (int)nGlobal_Elem_Par,
             (int)nGlobal_Elem_Storage);
    fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
    conn_buf = new int[nGlobal_Elem_Par*(N_POINTS_HEXAHEDRON+1)];
  }
  
  /*--- Load/write 1D buffers for the connectivity of each element type. ---*/
  
  if (surf_sol) {
    
    for (iElem = 0; iElem < nParallel_Line; iElem++) {
      iNode  = iElem*N_POINTS_LINE;
      iNode2 = iElem*(N_POINTS_LINE+1);
      conn_buf[iNode2+0] = N_POINTS_LINE;
      conn_buf[iNode2+1] = Conn_BoundLine_Par[iNode+0]-1;
      conn_buf[iNode2+2] = Conn_BoundLine_Par[iNode+1]-1;
    }
    if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                nParallel_Line*(N_POINTS_LINE+1));
    fwrite(conn_buf, sizeof(int),
           nParallel_Line*(N_POINTS_LINE+1), fhw);
    
    for (iElem = 0; iElem < nParallel_BoundTria; iElem++) {
      iNode  = iElem*N_POINTS_TRIANGLE;
      iNode2 = iElem*(N_POINTS_TRIANGLE+1);
      conn_buf[iNode2+0] = N_POINTS_TRIANGLE;
      conn_buf[iNode2+1] = Conn_BoundTria_Par[iNode+0]-1;
      conn_buf[iNode2+2] = Conn_BoundTria_Par[iNode+1]-1;
      conn_buf[iNode2+3] = Conn_BoundTria_Par[iNode+2]-1;
    }
    if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                nParallel_BoundTria*(N_POINTS_TRIANGLE+1));
    fwrite(conn_buf, sizeof(int),
           nParallel_BoundTria*(N_POINTS_TRIANGLE+1), fhw);
    
    for (iElem = 0; iElem < nParallel_BoundQuad; iElem++) {
      iNode  = iElem*N_POINTS_QUADRILATERAL;
      iNode2 = iElem*(N_POINTS_QUADRILATERAL+1);
      conn_buf[iNode2+0] = N_POINTS_QUADRILATERAL;
      conn_buf[iNode2+1] = Conn_BoundQuad_Par[iNode+0]-1;
      conn_buf[iNode2+2] = Conn_BoundQuad_Par[iNode+1]-1;
      conn_buf[iNode2+3] = Conn_BoundQuad_Par[iNode+2]-1;
      conn_buf[iNode2+4] = Conn_BoundQuad_Par[iNode+3]-1;
    }
    if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                nParallel_BoundQuad*(N_POINTS_QUADRILATERAL+1));
    fwrite(conn_buf, sizeof(int),
           nParallel_BoundQuad*(N_POINTS_QUADRILATERAL+1), fhw);
    
  } else {
    
    for (iElem = 0; iElem < nParallel_Tria; iElem++) {
      iNode  = iElem*N_POINTS_TRIANGLE;
      iNode2 = iElem*(N_POINTS_TRIANGLE+1);
      conn_buf[iNode2+0] = N_POINTS_TRIANGLE;
      conn_buf[iNode2+1] = Conn_Tria_Par[iNode+0]-1;
      conn_buf[iNode2+2] = Conn_Tria_Par[iNode+1]-1;
      conn_buf[iNode2+3] = Conn_Tria_Par[iNode+2]-1;
    }
    if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                              nParallel_Tria*(N_POINTS_TRIANGLE+1));
    fwrite(conn_buf, sizeof(int),
           nParallel_Tria*(N_POINTS_TRIANGLE+1), fhw);
    
    for (iElem = 0; iElem < nParallel_Quad; iElem++) {
      iNode  = iElem*N_POINTS_QUADRILATERAL;
      iNode2 = iElem*(N_POINTS_QUADRILATERAL+1);
      conn_buf[iNode2+0] = N_POINTS_QUADRILATERAL;
      conn_buf[iNode2+1] = Conn_Quad_Par[iNode+0]-1;
      conn_buf[iNode2+2] = Conn_Quad_Par[iNode+1]-1;
      conn_buf[iNode2+3] = Conn_Quad_Par[iNode+2]-1;
      conn_buf[iNode2+4] = Conn_Quad_Par[iNode+3]-1;
    }
    if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                              nParallel_Quad*(N_POINTS_QUADRILATERAL+1));
    fwrite(conn_buf, sizeof(int),
           nParallel_Quad*(N_POINTS_QUADRILATERAL+1), fhw);
    
    for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
      iNode  = iElem*N_POINTS_TETRAHEDRON;
      iNode2 = iElem*(N_POINTS_TETRAHEDRON+1);
      conn_buf[iNode2+0] = N_POINTS_TETRAHEDRON;
      conn_buf[iNode2+1] = Conn_Tetr_Par[iNode+0]-1;
      conn_buf[iNode2+2] = Conn_Tetr_Par[iNode+1]-1;
      conn_buf[iNode2+3] = Conn_Tetr_Par[iNode+2]-1;
      conn_buf[iNode2+4] = Conn_Tetr_Par[iNode+3]-1;
    }
    if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                              nParallel_Tetr*(N_POINTS_TETRAHEDRON+1));
    fwrite(conn_buf, sizeof(int),
           nParallel_Tetr*(N_POINTS_TETRAHEDRON+1), fhw);
    
    for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
      iNode  = iElem*N_POINTS_HEXAHEDRON;
      iNode2 = iElem*(N_POINTS_HEXAHEDRON+1);
      conn_buf[iNode2+0] = N_POINTS_HEXAHEDRON;
      conn_buf[iNode2+1] = Conn_Hexa_Par[iNode+0]-1;
      conn_buf[iNode2+2] = Conn_Hexa_Par[iNode+1]-1;
      conn_buf[iNode2+3] = Conn_Hexa_Par[iNode+2]-1;
      conn_buf[iNode2+4] = Conn_Hexa_Par[iNode+3]-1;
      conn_buf[iNode2+5] = Conn_Hexa_Par[iNode+4]-1;
      conn_buf[iNode2+6] = Conn_Hexa_Par[iNode+5]-1;
      conn_buf[iNode2+7] = Conn_Hexa_Par[iNode+6]-1;
      conn_buf[iNode2+8] = Conn_Hexa_Par[iNode+7]-1;
    }
    if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                              nParallel_Hexa*(N_POINTS_HEXAHEDRON+1));
    fwrite(conn_buf, sizeof(int),
           nParallel_Hexa*(N_POINTS_HEXAHEDRON+1), fhw);
    
    for (iElem = 0; iElem < nParallel_Pris; iElem++) {
      iNode  = iElem*N_POINTS_PRISM;
      iNode2 = iElem*(N_POINTS_PRISM+1);
      conn_buf[iNode2+0] = N_POINTS_PRISM;
      conn_buf[iNode2+1] = Conn_Pris_Par[iNode+0]-1;
      conn_buf[iNode2+2] = Conn_Pris_Par[iNode+1]-1;
      conn_buf[iNode2+3] = Conn_Pris_Par[iNode+2]-1;
      conn_buf[iNode2+4] = Conn_Pris_Par[iNode+3]-1;
      conn_buf[iNode2+5] = Conn_Pris_Par[iNode+4]-1;
      conn_buf[iNode2+6] = Conn_Pris_Par[iNode+5]-1;
    }
    if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                              nParallel_Pris*(N_POINTS_PRISM+1));
    fwrite(conn_buf, sizeof(int),
           nParallel_Pris*(N_POINTS_PRISM+1), fhw);
    
    for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
      iNode  = iElem*N_POINTS_PYRAMID;
      iNode2 = iElem*(N_POINTS_PYRAMID+1);
      conn_buf[iNode2+0] = N_POINTS_PYRAMID;
      conn_buf[iNode2+1] = Conn_Pyra_Par[iNode+0]-1;
      conn_buf[iNode2+2] = Conn_Pyra_Par[iNode+1]-1;
      conn_buf[iNode2+3] = Conn_Pyra_Par[iNode+2]-1;
      conn_buf[iNode2+4] = Conn_Pyra_Par[iNode+3]-1;
      conn_buf[iNode2+5] = Conn_Pyra_Par[iNode+4]-1;
    }
    if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                              nParallel_Pyra*(N_POINTS_PYRAMID+1));
    fwrite(conn_buf, sizeof(int),
           nParallel_Pyra*(N_POINTS_PYRAMID+1), fhw);
    
  }
  
  if (conn_buf != NULL) delete [] conn_buf;
  
  /*--- Load/write the cell type for all elements in the file. ---*/
  
  if (surf_sol) {
    SPRINTF (str_buf, "\nCELL_TYPES %i\n", SU2_TYPE::Int(nSurf_Elem_Par));
    fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  } else {
    SPRINTF (str_buf, "\nCELL_TYPES %i\n", SU2_TYPE::Int(nGlobal_Elem_Par));
    fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  }
  
  int *type_buf = NULL;
  if (surf_sol) {
    
    type_buf = new int[nSurf_Elem_Par];
    
    for (iElem = 0; iElem < nParallel_Line; iElem++) {
      type_buf[iElem] = LINE;
    }
    if (!BigEndian)
      SwapBytes((char *)type_buf, sizeof(int), nParallel_Line);
    fwrite(type_buf, sizeof(int), nParallel_Line, fhw);
    
    for (iElem = 0; iElem < nParallel_BoundTria; iElem++) {
      type_buf[iElem] = TRIANGLE;
    }
    if (!BigEndian)
      SwapBytes((char *)type_buf, sizeof(int), nParallel_BoundTria);
    fwrite(type_buf, sizeof(int), nParallel_BoundTria, fhw);
    
    for (iElem = 0; iElem < nParallel_BoundQuad; iElem++) {
      type_buf[iElem] = QUADRILATERAL;
    }
    if (!BigEndian)
      SwapBytes((char *)type_buf, sizeof(int), nParallel_BoundQuad);
    fwrite(type_buf, sizeof(int), nParallel_BoundQuad, fhw);
    
  } else {
    
    type_buf = new int[nGlobal_Elem_Par];
    
    for (iElem = 0; iElem < nParallel_Tria; iElem++) {
      type_buf[iElem] = TRIANGLE;
    }
    if (!BigEndian)
      SwapBytes((char *)type_buf, sizeof(int), nParallel_Tria);
    fwrite(type_buf, sizeof(int), nParallel_Tria, fhw);
    
    for (iElem = 0; iElem < nParallel_Quad; iElem++) {
      type_buf[iElem] = QUADRILATERAL;
    }
    if (!BigEndian)
      SwapBytes((char *)type_buf, sizeof(int), nParallel_Quad);
    fwrite(type_buf, sizeof(int), nParallel_Quad, fhw);
    
    for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
      type_buf[iElem] = TETRAHEDRON;
    }
    if (!BigEndian)
      SwapBytes((char *)type_buf, sizeof(int), nParallel_Tetr);
    fwrite(type_buf, sizeof(int), nParallel_Tetr, fhw);
    
    for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
      type_buf[iElem] = HEXAHEDRON;
    }
    if (!BigEndian)
      SwapBytes((char *)type_buf, sizeof(int), nParallel_Hexa);
    fwrite(type_buf, sizeof(int), nParallel_Hexa, fhw);
    
    for (iElem = 0; iElem < nParallel_Pris; iElem++) {
      type_buf[iElem] = PRISM;
    }
    if (!BigEndian)
      SwapBytes((char *)type_buf, sizeof(int), nParallel_Pris);
    fwrite(type_buf, sizeof(int), nParallel_Pris, fhw);
    
    for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
      type_buf[iElem] = PYRAMID;
    }
    if (!BigEndian)
      SwapBytes((char *)type_buf, sizeof(int), nParallel_Pyra);
    fwrite(type_buf, sizeof(int), nParallel_Pyra, fhw);
    
  }
  
  if (type_buf != NULL) delete [] type_buf;
  
  /*--- Now write the scalar and vector data (reuse the counts above). ---*/
  
  SPRINTF (str_buf, "\nPOINT_DATA %i\n", (int)GlobalPoint);
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  
  unsigned short varStart = 2;
  if (nDim == 3) varStart++;
  
  /*--- Need to adjust container location to avoid PointID tag and coords. ---*/
  
  unsigned short iField, VarCounter = varStart;
  for (iField = varStart; iField < Variable_Names.size(); iField++) {
    
    fieldname = Variable_Names[iField];
    fieldname.erase(remove(fieldname.begin(), fieldname.end(), '"'),
                    fieldname.end());
    
    bool output_variable = true, isVector = false;
    size_t found = Variable_Names[iField].find("_x");
    if (found!=string::npos) {
      output_variable = true;
      isVector = true;
    }
    found = Variable_Names[iField].find("_y");
    if (found!=string::npos) {
      //skip
      output_variable = false;
      VarCounter++;
    }
    found = Variable_Names[iField].find("_z");
    if (found!=string::npos) {
      //skip
      output_variable = false;
      VarCounter++;
    }
    
    if (output_variable && isVector) {
      
      fieldname.erase(fieldname.end()-2,fieldname.end());
      SPRINTF (str_buf, "\nVECTORS %s float\n", fieldname.c_str());
      fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
      
      /*--- Prepare the 1D data buffer on this rank. ---*/
      
      float *vec_buf = new float[GlobalPoint*NCOORDS];
      
      /*--- For now, create a temp 1D buffer to load up the data for writing.
       This will be replaced with a derived data type most likely. ---*/
      
      float val = 0.0;
      for (iPoint = 0; iPoint < GlobalPoint; iPoint++)
        for (iDim = 0; iDim < NCOORDS; iDim++) {
          if (nDim == 2 && iDim == 2) {
            vec_buf[iPoint*NCOORDS + iDim] = 0.0;
          } else {
            val = (float)SU2_TYPE::GetValue(Data[VarCounter+iDim][iPoint]);
            vec_buf[iPoint*NCOORDS + iDim] = val;
          }
        }
      if (!BigEndian)
        SwapBytes((char *)vec_buf, sizeof(float), NCOORDS*GlobalPoint);
      fwrite(vec_buf, sizeof(float), NCOORDS*GlobalPoint, fhw);
      
      delete [] vec_buf;
      
      VarCounter++;
      
    } else if (output_variable) {
      
      SPRINTF (str_buf, "\nSCALARS %s float 1\n", fieldname.c_str());
      fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
      
      SPRINTF (str_buf, "LOOKUP_TABLE default\n");
      fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
      
      /*--- Prepare the 1D data buffer on this rank. ---*/
      
      float *scalar_buf = new float[GlobalPoint];
      
      /*--- For now, create a temp 1D buffer to load up the data for writing.
       This will be replaced with a derived data type most likely. ---*/
      
      for (iPoint = 0; iPoint < GlobalPoint; iPoint++) {
        float val = (float)SU2_TYPE::GetValue(Data[VarCounter][iPoint]);
        scalar_buf[iPoint] = val;
      }
      if (!BigEndian)
        SwapBytes((char *)scalar_buf, sizeof(float), GlobalPoint);
      fwrite(scalar_buf, sizeof(float), GlobalPoint, fhw);
      
      delete [] scalar_buf;
      
      VarCounter++;
    }
    
  }
  
  /*--- Close the file. ---*/
  
  fclose(fhw);
  
#else
  
  /*--- Parallel binary output using MPI I/O. ---*/
  
  MPI_File fhw;
  SU2_MPI::Status status;
  MPI_Datatype etype, filetype;
  MPI_Offset disp, disp2;
  int ierr;
  
  /*--- All ranks open the file using MPI. Here, we try to open the file with
   exclusive so that an error is generated if the file exists. We always want
   to write a fresh output file, so we delete any existing files and create
   a new one. ---*/
  
  ierr = MPI_File_open(MPI_COMM_WORLD, fname,
                       MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                       MPI_INFO_NULL, &fhw);
  if (ierr != MPI_SUCCESS)  {
    MPI_File_close(&fhw);
    if (rank == 0)
      MPI_File_delete(fname, MPI_INFO_NULL);
    ierr = MPI_File_open(MPI_COMM_WORLD, fname,
                         MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                         MPI_INFO_NULL, &fhw);
  }
  
  /*--- Error check opening the file. ---*/
  
  if (ierr) {
    SU2_MPI::Error(string("Unable to open VTK binary legacy file ") +
                   string(fname), CURRENT_FUNCTION);
  }
  
  /*--- Set pointer to our output data for simplicity. ---*/
  
  su2double **Data;
  if (surf_sol) {
    Data = Parallel_Surf_Data;
  } else {
    Data = Parallel_Data;
  }
  
  /*--- Write the initial strings to the file. Only the master will
   write the header lines, but all ranks will store the offsets. ---*/
  
  disp = 0;
  strcpy(str_buf, "# vtk DataFile Version 3.0\n");
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  
  strcpy(str_buf, "vtk output\n");
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  
  strcpy(str_buf, "BINARY\n");
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  
  strcpy(str_buf, "DATASET UNSTRUCTURED_GRID\n");
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  
  /*--- Communicate the number of total points that will be
   written by each rank. After this communication, each proc knows how
   many poinnts will be written before its location in the file and the
   offsets can be correctly set. ---*/
  
  unsigned long myPoint, GlobalPoint;
  if (surf_sol) {
    GlobalPoint = nGlobal_Surf_Poin;
    myPoint     = nSurf_Poin_Par;
  } else {
    GlobalPoint = nGlobal_Poin_Par;
    myPoint     = nParallel_Poin;
  }
  
  int *nPoint_Snd = new int[size+1];
  int *nPoint_Cum = new int[size+1];
  
  nPoint_Snd[0] = 0; nPoint_Cum[0] = 0;
  for (int ii=1; ii < size; ii++) {
    nPoint_Snd[ii] = myPoint; nPoint_Cum[ii] = 0;
  }
  nPoint_Snd[size] = myPoint; nPoint_Cum[size] = 0;
  
  /*--- Communicate the local counts to all ranks for building offsets. ---*/
  
  SU2_MPI::Alltoall(&(nPoint_Snd[1]), 1, MPI_INT,
                    &(nPoint_Cum[1]), 1, MPI_INT, MPI_COMM_WORLD);
  
  /*--- Put the counters into cumulative storage format. ---*/
  
  for (int ii = 0; ii < size; ii++) {
    nPoint_Cum[ii+1] += nPoint_Cum[ii];
  }
  
  SPRINTF(str_buf, "POINTS %i float\n", SU2_TYPE::Int(GlobalPoint));
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  
  /*--- Load/write the 1D buffer of point coordinates. Note that we
   always have 3 coordinate dimensions, even for 2D problems. ---*/
  
  float *coord_buf = new float[myPoint*NCOORDS];
  for (iPoint = 0; iPoint < myPoint; iPoint++) {
    for (iDim = 0; iDim < NCOORDS; iDim++) {
      if (nDim == 2 && iDim == 2) {
        coord_buf[iPoint*NCOORDS + iDim] = 0.0;
      } else {
        float val = (float)SU2_TYPE::GetValue(Data[iDim][iPoint]);
        coord_buf[iPoint*NCOORDS + iDim] = val;
      }
    }
  }
  if (!BigEndian) SwapBytes((char *)coord_buf, sizeof(float), myPoint*NCOORDS);

  /*--- We will write the point coordinates as floats. ---*/
  
  etype = MPI_FLOAT;
  
  /*--- Define a derived datatype for this ranks contiguous
   chunk of data that will be placed in the file. ---*/
  
  MPI_Type_contiguous(myPoint*NCOORDS, MPI_FLOAT, &filetype);
  MPI_Type_commit(&filetype);
  
  /*--- Compute the offset for this rank's linear partition of the
   data in bytes. ---*/
  
  disp2 = disp + NCOORDS*nPoint_Cum[rank]*sizeof(float);
  
  /*--- Set the view for the MPI file write, i.e., describe the
   location in the file that this rank "sees" for writing its
   piece of the file. ---*/
  
  MPI_File_set_view(fhw, disp2, etype, filetype,
                    (char*)"native", MPI_INFO_NULL);
  
  /*--- Collective call for all ranks to write simultaneously. ---*/
  
  MPI_File_write_all(fhw, coord_buf, myPoint*NCOORDS, MPI_FLOAT, &status);
  
  /*--- Update the displacement position for MPI IO. ---*/
  
  disp += NCOORDS*nPoint_Cum[size]*sizeof(float);
  
  /*--- Free the derived datatype and coordinate array. ---*/
  
  MPI_Type_free(&filetype);
  delete [] coord_buf;
  
  /*--- Compute our local number of elements, the required storage,
   and reduce the total number of elements and storage globally. ---*/
  
  unsigned long nTot_Line, nTot_BoundTria, nTot_BoundQuad;
  unsigned long nTot_Tria, nTot_Quad;
  unsigned long nTot_Tetr, nTot_Hexa, nTot_Pris, nTot_Pyra;
  unsigned long myElem, myElemStorage, GlobalElem, GlobalElemStorage;
  
  if (surf_sol) {
    
    SU2_MPI::Allreduce(&nParallel_Line,      &nTot_Line,      1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nParallel_BoundTria, &nTot_BoundTria, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nParallel_BoundQuad, &nTot_BoundQuad, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    
    myElem        = (nParallel_Line   + nParallel_BoundTria   +
                     nParallel_BoundQuad);
    myElemStorage = (nParallel_Line*3 + nParallel_BoundTria*4 +
                     nParallel_BoundQuad*5);
    
    GlobalElem        = nTot_Line   + nTot_BoundTria   + nTot_BoundQuad;
    GlobalElemStorage = nTot_Line*3 + nTot_BoundTria*4 + nTot_BoundQuad*5;
    
  } else {
    
    SU2_MPI::Allreduce(&nParallel_Tria, &nTot_Tria, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nParallel_Quad, &nTot_Quad, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nParallel_Tetr, &nTot_Tetr, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nParallel_Hexa, &nTot_Hexa, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nParallel_Pris, &nTot_Pris, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nParallel_Pyra, &nTot_Pyra, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    
    myElem        = (nParallel_Tria + nParallel_Quad + nParallel_Tetr +
                     nParallel_Hexa + nParallel_Pris + nParallel_Pyra);
    myElemStorage = (nParallel_Tria*4 + nParallel_Quad*5 + nParallel_Tetr*5 +
                     nParallel_Hexa*9 + nParallel_Pris*7 + nParallel_Pyra*6);
    
    GlobalElem        = (nTot_Tria   + nTot_Quad   + nTot_Tetr   +
                         nTot_Hexa   + nTot_Pris   + nTot_Pyra);
    GlobalElemStorage = (nTot_Tria*4 + nTot_Quad*5 + nTot_Tetr*5 +
                         nTot_Hexa*9 + nTot_Pris*7 + nTot_Pyra*6);
    
  }
  
  /*--- Communicate the number of total cells/storage that will be
   written by each rank. After this communication, each proc knows how
   many cells will be written before its location in the file and the
   offsets can be correctly set. ---*/
  
  int *nElem_Snd = new int[size+1]; int *nElemStorage_Snd = new int[size+1];
  int *nElem_Cum = new int[size+1]; int *nElemStorage_Cum = new int[size+1];
  
  nElem_Snd[0] = 0; nElemStorage_Snd[0] = 0;
  nElem_Cum[0] = 0; nElemStorage_Cum[0] = 0;
  for (int ii=1; ii < size; ii++) {
    nElem_Snd[ii] = myElem; nElemStorage_Snd[ii] = myElemStorage;
    nElem_Cum[ii] = 0;      nElemStorage_Cum[ii] = 0;
  }
  nElem_Snd[size] = myElem; nElemStorage_Snd[size] = myElemStorage;
  nElem_Cum[size] = 0;      nElemStorage_Cum[size] = 0;
  
  /*--- Communicate the local counts to all ranks for building offsets. ---*/
  
  SU2_MPI::Alltoall(&(nElem_Snd[1]), 1, MPI_INT,
                    &(nElem_Cum[1]), 1, MPI_INT, MPI_COMM_WORLD);
  
  SU2_MPI::Alltoall(&(nElemStorage_Snd[1]), 1, MPI_INT,
                    &(nElemStorage_Cum[1]), 1, MPI_INT, MPI_COMM_WORLD);
  
  /*--- Put the counters into cumulative storage format. ---*/
  
  for (int ii = 0; ii < size; ii++) {
    nElem_Cum[ii+1]        += nElem_Cum[ii];
    nElemStorage_Cum[ii+1] += nElemStorage_Cum[ii];
  }
  
  /*--- Reset the file view before writing the next ASCII line for cells. ---*/
  
  MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE,
                    (char*)"native", MPI_INFO_NULL);
  SPRINTF(str_buf, "\nCELLS %i %i\n", SU2_TYPE::Int(GlobalElem),
          SU2_TYPE::Int(GlobalElemStorage));
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  
  /*--- Load/write 1D buffers for the connectivity of each element type. ---*/
  
  int *conn_buf = new int[myElemStorage];
  unsigned long iStorage = 0;
  
  if (surf_sol) {
    
    for (iElem = 0; iElem < nParallel_Line; iElem++) {
      iNode  = iElem*N_POINTS_LINE;
      conn_buf[iStorage+0] = N_POINTS_LINE;
      conn_buf[iStorage+1] = Conn_BoundLine_Par[iNode+0]-1;
      conn_buf[iStorage+2] = Conn_BoundLine_Par[iNode+1]-1;
      iStorage += (N_POINTS_LINE+1);
    }
    
    for (iElem = 0; iElem < nParallel_BoundTria; iElem++) {
      iNode  = iElem*N_POINTS_TRIANGLE;
      conn_buf[iStorage+0] = N_POINTS_TRIANGLE;
      conn_buf[iStorage+1] = Conn_BoundTria_Par[iNode+0]-1;
      conn_buf[iStorage+2] = Conn_BoundTria_Par[iNode+1]-1;
      conn_buf[iStorage+3] = Conn_BoundTria_Par[iNode+2]-1;
      iStorage += (N_POINTS_TRIANGLE+1);
    }
    
    for (iElem = 0; iElem < nParallel_BoundQuad; iElem++) {
      iNode  = iElem*N_POINTS_QUADRILATERAL;
      conn_buf[iStorage+0] = N_POINTS_QUADRILATERAL;
      conn_buf[iStorage+1] = Conn_BoundQuad_Par[iNode+0]-1;
      conn_buf[iStorage+2] = Conn_BoundQuad_Par[iNode+1]-1;
      conn_buf[iStorage+3] = Conn_BoundQuad_Par[iNode+2]-1;
      conn_buf[iStorage+4] = Conn_BoundQuad_Par[iNode+3]-1;
      iStorage += (N_POINTS_QUADRILATERAL+1);
    }
    
  } else {
    
    for (iElem = 0; iElem < nParallel_Tria; iElem++) {
      iNode  = iElem*N_POINTS_TRIANGLE;
      conn_buf[iStorage+0] = N_POINTS_TRIANGLE;
      conn_buf[iStorage+1] = Conn_Tria_Par[iNode+0]-1;
      conn_buf[iStorage+2] = Conn_Tria_Par[iNode+1]-1;
      conn_buf[iStorage+3] = Conn_Tria_Par[iNode+2]-1;
      iStorage += (N_POINTS_TRIANGLE+1);
    }
    
    for (iElem = 0; iElem < nParallel_Quad; iElem++) {
      iNode  = iElem*N_POINTS_QUADRILATERAL;
      conn_buf[iStorage+0] = N_POINTS_QUADRILATERAL;
      conn_buf[iStorage+1] = Conn_Quad_Par[iNode+0]-1;
      conn_buf[iStorage+2] = Conn_Quad_Par[iNode+1]-1;
      conn_buf[iStorage+3] = Conn_Quad_Par[iNode+2]-1;
      conn_buf[iStorage+4] = Conn_Quad_Par[iNode+3]-1;
      iStorage += (N_POINTS_QUADRILATERAL+1);
      
    }
    
    for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
      iNode  = iElem*N_POINTS_TETRAHEDRON;
      conn_buf[iStorage+0] = N_POINTS_TETRAHEDRON;
      conn_buf[iStorage+1] = Conn_Tetr_Par[iNode+0]-1;
      conn_buf[iStorage+2] = Conn_Tetr_Par[iNode+1]-1;
      conn_buf[iStorage+3] = Conn_Tetr_Par[iNode+2]-1;
      conn_buf[iStorage+4] = Conn_Tetr_Par[iNode+3]-1;
      iStorage += (N_POINTS_TETRAHEDRON+1);
      
    }
    
    for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
      iNode  = iElem*N_POINTS_HEXAHEDRON;
      conn_buf[iStorage+0] = N_POINTS_HEXAHEDRON;
      conn_buf[iStorage+1] = Conn_Hexa_Par[iNode+0]-1;
      conn_buf[iStorage+2] = Conn_Hexa_Par[iNode+1]-1;
      conn_buf[iStorage+3] = Conn_Hexa_Par[iNode+2]-1;
      conn_buf[iStorage+4] = Conn_Hexa_Par[iNode+3]-1;
      conn_buf[iStorage+5] = Conn_Hexa_Par[iNode+4]-1;
      conn_buf[iStorage+6] = Conn_Hexa_Par[iNode+5]-1;
      conn_buf[iStorage+7] = Conn_Hexa_Par[iNode+6]-1;
      conn_buf[iStorage+8] = Conn_Hexa_Par[iNode+7]-1;
      iStorage += (N_POINTS_HEXAHEDRON+1);
    }
    
    for (iElem = 0; iElem < nParallel_Pris; iElem++) {
      iNode  = iElem*N_POINTS_PRISM;
      conn_buf[iStorage+0] = N_POINTS_PRISM;
      conn_buf[iStorage+1] = Conn_Pris_Par[iNode+0]-1;
      conn_buf[iStorage+2] = Conn_Pris_Par[iNode+1]-1;
      conn_buf[iStorage+3] = Conn_Pris_Par[iNode+2]-1;
      conn_buf[iStorage+4] = Conn_Pris_Par[iNode+3]-1;
      conn_buf[iStorage+5] = Conn_Pris_Par[iNode+4]-1;
      conn_buf[iStorage+6] = Conn_Pris_Par[iNode+5]-1;
      iStorage += (N_POINTS_PRISM+1);
    }
    
    for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
      iNode  = iElem*N_POINTS_PYRAMID;
      conn_buf[iStorage+0] = N_POINTS_PYRAMID;
      conn_buf[iStorage+1] = Conn_Pyra_Par[iNode+0]-1;
      conn_buf[iStorage+2] = Conn_Pyra_Par[iNode+1]-1;
      conn_buf[iStorage+3] = Conn_Pyra_Par[iNode+2]-1;
      conn_buf[iStorage+4] = Conn_Pyra_Par[iNode+3]-1;
      conn_buf[iStorage+5] = Conn_Pyra_Par[iNode+4]-1;
      iStorage += (N_POINTS_PYRAMID+1);
    }
    
  }
  if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int), myElemStorage);
  
  /*--- We write the connectivity with MPI_INTs. ---*/
  
  etype = MPI_INT;
  
  /*--- Define a derived datatype for this ranks contiguous
   chunk of data that will be placed in the file. ---*/
  
  MPI_Type_contiguous(myElemStorage, MPI_INT, &filetype);
  MPI_Type_commit(&filetype);
  
  /*--- Compute the offset for this rank's linear partition of the
   data in bytes. ---*/
  
  disp2 = (disp + nElemStorage_Cum[rank]*sizeof(int));
  
  /*--- Set the view for the MPI file write, i.e., describe the
   location in the file that this rank "sees" for writing its
   piece of the file. ---*/
  
  MPI_File_set_view(fhw, disp2, etype, filetype,
                    (char*)"native", MPI_INFO_NULL);
  
  /*--- Collective call for all ranks to write simultaneously. ---*/
  
  MPI_File_write_all(fhw, conn_buf, myElemStorage, MPI_INT, &status);
  
  /*--- Update the displacement position for MPI IO. ---*/
  
  disp += nElemStorage_Cum[size]*sizeof(int);
  
  /*--- Free the derived datatype. ---*/
  
  MPI_Type_free(&filetype);
  delete [] conn_buf;
  
  /*--- Load/write the cell type for all elements in the file. ---*/
  
  MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE,
                    (char*)"native", MPI_INFO_NULL);
  SPRINTF (str_buf, "\nCELL_TYPES %i\n", SU2_TYPE::Int(GlobalElem));
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  
  int *type_buf = new int[myElem];
  unsigned long jElem = 0;
  
  if (surf_sol) {
    for (iElem = 0; iElem < nParallel_Line; iElem++) {
      type_buf[jElem] = LINE; jElem++;
    }
    for (iElem = 0; iElem < nParallel_BoundTria; iElem++) {
      type_buf[jElem] = TRIANGLE; jElem++;
    }
    for (iElem = 0; iElem < nParallel_BoundQuad; iElem++) {
      type_buf[jElem] = QUADRILATERAL; jElem++;
    }
  } else {
    for (iElem = 0; iElem < nParallel_Tria; iElem++) {
      type_buf[jElem] = TRIANGLE; jElem++;
    }
    for (iElem = 0; iElem < nParallel_Quad; iElem++) {
      type_buf[jElem] = QUADRILATERAL; jElem++;
    }
    for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
      type_buf[jElem] = TETRAHEDRON; jElem++;
    }
    for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
      type_buf[jElem] = HEXAHEDRON; jElem++;
    }
    for (iElem = 0; iElem < nParallel_Pris; iElem++) {
      type_buf[jElem] = PRISM; jElem++;
    }
    for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
      type_buf[jElem] = PYRAMID; jElem++;
    }
  }
  if (!BigEndian) SwapBytes((char *)type_buf, sizeof(int), myElem);

  /*--- We write the cell types with MPI_INTs. ---*/
  
  etype = MPI_INT;
  
  /*--- Define a derived datatype for this ranks contiguous
   chunk of data that will be placed in the file. ---*/
  
  MPI_Type_contiguous(myElem, MPI_INT, &filetype);
  MPI_Type_commit(&filetype);
  
  /*--- Compute the offset for this rank's linear partition of the
   data in bytes. ---*/
  
  disp2 = (disp + nElem_Cum[rank]*sizeof(int));
  
  /*--- Set the view for the MPI file write, i.e., describe the
   location in the file that this rank "sees" for writing its
   piece of the file. ---*/
  
  MPI_File_set_view(fhw, disp2, etype, filetype,
                    (char*)"native", MPI_INFO_NULL);
  
  /*--- Collective call for all ranks to write simultaneously. ---*/
  
  MPI_File_write_all(fhw, type_buf, myElem, MPI_INT, &status);
  
  /*--- Update the displacement position for MPI IO. ---*/
  
  disp += nElem_Cum[size]*sizeof(int);
  
  /*--- Free the derived datatype. ---*/
  
  MPI_Type_free(&filetype);
  if (type_buf != NULL) delete [] type_buf;
  
  /*--- Now write the scalar and vector point data. ---*/
  
  MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE,
                    (char*)"native", MPI_INFO_NULL);
  SPRINTF (str_buf, "\nPOINT_DATA %i\n", SU2_TYPE::Int(GlobalPoint));
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  
  /*--- Adjust container start location to avoid point coords. ---*/
  
  unsigned short varStart = 2;
  if (nDim == 3) varStart++;
  
  /*--- Loop over all variables that have been registered in the output. ---*/
  
  unsigned short iField, VarCounter = varStart;
  for (iField = varStart; iField < Variable_Names.size(); iField++) {
    
    fieldname = Variable_Names[iField];
    fieldname.erase(remove(fieldname.begin(), fieldname.end(), '"'),
                    fieldname.end());
    
    /*--- Check whether this field is a vector or scalar. ---*/
    
    bool output_variable = true, isVector = false;
    size_t found = Variable_Names[iField].find("_x");
    if (found!=string::npos) {
      output_variable = true;
      isVector        = true;
    }
    found = Variable_Names[iField].find("_y");
    if (found!=string::npos) {
      /*--- We have found a vector, so skip the Y component. ---*/
      output_variable = false;
      VarCounter++;
    }
    found = Variable_Names[iField].find("_z");
    if (found!=string::npos) {
      /*--- We have found a vector, so skip the Z component. ---*/
      output_variable = false;
      VarCounter++;
    }
    
    /*--- Write the point data as an <X,Y,Z> vector or a scalar. ---*/
    
    if (output_variable && isVector) {
      
      /*--- Adjust the string name to remove the leading "X-" ---*/
      
      fieldname.erase(fieldname.end()-2,fieldname.end());
      MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE,
                        (char*)"native", MPI_INFO_NULL);
      SPRINTF (str_buf, "\nVECTORS %s float\n", fieldname.c_str());
      if (rank == MASTER_NODE)
        MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                          MPI_CHAR, MPI_STATUS_IGNORE);
      disp += strlen(str_buf)*sizeof(char);
      
      /*--- Prepare the 1D data buffer on this rank. ---*/
      
      float *vec_buf = new float[myPoint*NCOORDS];
      
      /*--- Load up the buffer for writing this rank's vector data. ---*/
      
      float val = 0.0;
      for (iPoint = 0; iPoint < myPoint; iPoint++) {
        for (iDim = 0; iDim < NCOORDS; iDim++) {
          if (nDim == 2 && iDim == 2) {
            vec_buf[iPoint*NCOORDS + iDim] = 0.0;
          } else {
            val = (float)SU2_TYPE::GetValue(Data[VarCounter+iDim][iPoint]);
            vec_buf[iPoint*NCOORDS + iDim] = val;
          }
        }
      }
      if (!BigEndian)
        SwapBytes((char *)vec_buf, sizeof(float), myPoint*NCOORDS);

      /*--- We will write the point data as floats. ---*/
      
      etype = MPI_FLOAT;
      
      /*--- Define a derived datatype for this ranks contiguous
       chunk of data that will be placed in the file. ---*/
      
      MPI_Type_contiguous(myPoint*NCOORDS, MPI_FLOAT, &filetype);
      MPI_Type_commit(&filetype);
      
      /*--- Compute the offset for this rank's linear partition of the
       data in bytes. ---*/
      
      disp2 = disp + NCOORDS*nPoint_Cum[rank]*sizeof(float);
      
      /*--- Set the view for the MPI file write, i.e., describe the
       location in the file that this rank "sees" for writing its
       piece of the file. ---*/
      
      MPI_File_set_view(fhw, disp2, etype, filetype,
                        (char*)"native", MPI_INFO_NULL);
      
      /*--- Collective call for all ranks to write simultaneously. ---*/
      
      MPI_File_write_all(fhw, vec_buf, myPoint*NCOORDS, MPI_FLOAT, &status);
      
      /*--- Update the displacement position for MPI IO. ---*/
      
      disp += NCOORDS*nPoint_Cum[size]*sizeof(float);
      
      /*--- Free the derived datatype and coordinate array. ---*/
      
      MPI_Type_free(&filetype);
      delete [] vec_buf; vec_buf = NULL;
      
      VarCounter++;
      
    } else if (output_variable) {
      
      MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE,
                        (char*)"native", MPI_INFO_NULL);
      SPRINTF (str_buf, "\nSCALARS %s float 1\n", fieldname.c_str());
      if (rank == MASTER_NODE)
        MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                          MPI_CHAR, MPI_STATUS_IGNORE);
      disp += strlen(str_buf)*sizeof(char);
      
      MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE,
                        (char*)"native", MPI_INFO_NULL);
      SPRINTF (str_buf, "LOOKUP_TABLE default\n");
      if (rank == MASTER_NODE)
        MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                          MPI_CHAR, MPI_STATUS_IGNORE);
      disp += strlen(str_buf)*sizeof(char);
      
      /*--- Prepare the 1D data buffer on this rank. ---*/
      
      float *scalar_buf = new float[myPoint];
      
      /*--- For now, create a temp 1D buffer to load up the data for writing.
       This will be replaced with a derived data type most likely. ---*/
      
      for (iPoint = 0; iPoint < myPoint; iPoint++) {
        float val = (float)SU2_TYPE::GetValue(Data[VarCounter][iPoint]);
        scalar_buf[iPoint] = val;
      }
      if (!BigEndian) SwapBytes((char *)scalar_buf, sizeof(float), myPoint);
      
      /*--- We will write the point data as floats. ---*/
      
      etype = MPI_FLOAT;
      
      /*--- Define a derived datatype for this ranks contiguous
       chunk of data that will be placed in the file. ---*/
      
      MPI_Type_contiguous(myPoint, MPI_FLOAT, &filetype);
      MPI_Type_commit(&filetype);
      
      /*--- Compute the offset for this rank's linear partition of the
       data in bytes. ---*/
      
      disp2 = disp + nPoint_Cum[rank]*sizeof(float);
      
      /*--- Set the view for the MPI file write, i.e., describe the
       location in the file that this rank "sees" for writing its
       piece of the file. ---*/
      
      MPI_File_set_view(fhw, disp2, etype, filetype,
                        (char*)"native", MPI_INFO_NULL);
      
      /*--- Collective call for all ranks to write simultaneously. ---*/
      
      MPI_File_write_all(fhw, scalar_buf, myPoint, MPI_FLOAT, &status);
      
      /*--- Update the displacement position for MPI IO. ---*/
      
      disp += nPoint_Cum[size]*sizeof(float);
      
      /*--- Free the derived datatype and coordinate array. ---*/
      
      MPI_Type_free(&filetype);
      delete [] scalar_buf; scalar_buf = NULL;
      
      VarCounter++;
    }
    
  }
  
  /*--- All ranks close the file after writing. ---*/
  
  MPI_File_close(&fhw);
  
  /*--- Delete the offset counters that we needed for MPI IO. ---*/
  
  delete [] nElem_Snd;        delete [] nElem_Cum;
  delete [] nElemStorage_Snd; delete [] nElemStorage_Cum;
  delete [] nPoint_Snd;       delete [] nPoint_Cum;
  
#endif
  
}
