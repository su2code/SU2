#include "../../../include/output/filewriter/CParaviewBinaryFileWriter.hpp"

const string CParaviewBinaryFileWriter::fileExt = ".vtk";

CParaviewBinaryFileWriter::CParaviewBinaryFileWriter(vector<string> fields, unsigned short nDim, string fileName, 
                                                     CParallelDataSorter *dataSorter) : 
  CFileWriter(std::move(fields), std::move(fileName), dataSorter, fileExt, nDim){}


CParaviewBinaryFileWriter::~CParaviewBinaryFileWriter(){
  
}

void CParaviewBinaryFileWriter::Write_Data(){
    
  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }
  
  unsigned short iDim;
  
  unsigned long iPoint, iElem;
  
  ofstream Paraview_File;
    
  const int MAX_STRING_LENGTH = 255;
  char str_buf[MAX_STRING_LENGTH], fname[100];
  
  const int NCOORDS = 3;

  strcpy(fname, fileName.c_str());
  
  /* Check for big endian. We have to swap bytes otherwise. 
   * Since size of character is 1 byte when the character pointer
   *  is de-referenced it will contain only first byte of integer. ---*/
  
  bool BigEndian = false;  
  unsigned int i = 1;  
  char *c = (char*)&i;
  if (*c) BigEndian = false;
  else BigEndian = true;
  
  file_size = 0.0;
  
  /*--- Set a timer for the file writing. ---*/
  
#ifndef HAVE_MPI
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StartTime = MPI_Wtime();
#endif
  
  /*--- Serial implementation in case we have not compiled with MPI. ---*/
  
#ifndef HAVE_MPI
  
  FILE* fhw;
  fhw = fopen(fname, "wb");
  
  unsigned long iNode2;
  unsigned long nGlobal_Elem_Storage;
  
  /*--- Error check for opening the file. ---*/
  
  if (!fhw) {
    SU2_MPI::Error(string("Unable to open VTK binary legacy file ") +
                   fileName, CURRENT_FUNCTION);
  }
  
  /*--- File header written in ASCII. ---*/
  
  strcpy(str_buf, "# vtk DataFile Version 3.0\n");
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  file_size += sizeof(char)*strlen(str_buf);
  
  strcpy(str_buf, "vtk output\n");
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  file_size += sizeof(char)*strlen(str_buf);
  
  strcpy(str_buf, "BINARY\n");
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  file_size += sizeof(char)*strlen(str_buf);
  
  strcpy(str_buf, "DATASET UNSTRUCTURED_GRID\n");
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  file_size += sizeof(char)*strlen(str_buf);
  
  /*--- Write the point coordinates. ---*/
  
  unsigned long GlobalPoint = dataSorter->GetnPointsGlobal();
  
  SPRINTF(str_buf, "POINTS %i float\n", (int)GlobalPoint);
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  file_size += sizeof(char)*strlen(str_buf);
  
  /*--- Load/write the 1D buffer of point coordinates. ---*/
  
  float *coord_buf = new float[GlobalPoint*NCOORDS];
  for (iPoint = 0; iPoint < GlobalPoint; iPoint++) {
    for (iDim = 0; iDim < NCOORDS; iDim++) {
      if (nDim == 2 && iDim == 2) {
        coord_buf[iPoint*NCOORDS + iDim] = 0.0;
      } else {
        float val = (float)dataSorter->GetData(iDim,iPoint);
        coord_buf[iPoint*NCOORDS + iDim] = val;
      }
    }
  }
  if (!BigEndian) SwapBytes((char *)coord_buf, sizeof(float), 3*GlobalPoint);
  
  fwrite(coord_buf, sizeof(float), 3*GlobalPoint, fhw);
  file_size += sizeof(char)*3*GlobalPoint;
  
  delete [] coord_buf;
  
  /*--- Write the connectivity data. ---*/
  
  unsigned long nTot_Line;
  unsigned long nTot_Tria, nTot_Quad;
  unsigned long nTot_Tetr, nTot_Hexa, nTot_Pris, nTot_Pyra;
  nTot_Line = dataSorter->GetnElem(LINE);  
  nTot_Tria = dataSorter->GetnElem(TRIANGLE);
  nTot_Quad = dataSorter->GetnElem(QUADRILATERAL);
  nTot_Tetr = dataSorter->GetnElem(TETRAHEDRON);
  nTot_Hexa = dataSorter->GetnElem(HEXAHEDRON);
  nTot_Pris = dataSorter->GetnElem(PRISM);
  nTot_Pyra = dataSorter->GetnElem(PYRAMID);
  nGlobal_Elem_Storage = (nTot_Line*3 + nTot_Tria*4 + nTot_Quad*5 + nTot_Tetr*5 +
                          nTot_Hexa*9 + nTot_Pris*7 + nTot_Pyra*6);
  
  int *conn_buf = NULL;
  
  SPRINTF (str_buf, "\nCELLS %i %i\n", (int)dataSorter->GetnElem(),
           (int)nGlobal_Elem_Storage);
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  file_size += sizeof(char)*strlen(str_buf);
  
  conn_buf = new int[dataSorter->GetnElem()*(N_POINTS_HEXAHEDRON+1)];
  
  
  /*--- Load/write 1D buffers for the connectivity of each element type. ---*/
  
  
  for (iElem = 0; iElem < nTot_Line; iElem++) {
    iNode2 = iElem*(N_POINTS_LINE+1);
    conn_buf[iNode2+0] = N_POINTS_LINE;
    conn_buf[iNode2+1] = dataSorter->GetElem_Connectivity(LINE, iElem, 0)-1;
    conn_buf[iNode2+2] = dataSorter->GetElem_Connectivity(LINE, iElem, 1)-1;
  }
  if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                            nTot_Line*(N_POINTS_LINE+1));
  fwrite(conn_buf, sizeof(int),
         nTot_Line*(N_POINTS_LINE+1), fhw);
  
  file_size += sizeof(int)*nTot_Line*(N_POINTS_LINE+1);
  
  for (iElem = 0; iElem < nTot_Tria; iElem++) {
    iNode2 = iElem*(N_POINTS_TRIANGLE+1);
    conn_buf[iNode2+0] = N_POINTS_TRIANGLE;
    conn_buf[iNode2+1] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 0)-1;
    conn_buf[iNode2+2] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 1)-1;
    conn_buf[iNode2+3] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 2)-1;
  }
  if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                            nTot_Tria*(N_POINTS_TRIANGLE+1));
  fwrite(conn_buf, sizeof(int),
         nTot_Tria*(N_POINTS_TRIANGLE+1), fhw);
  file_size += sizeof(int)*nTot_Tria*(N_POINTS_TRIANGLE+1);
  
  for (iElem = 0; iElem < nTot_Quad; iElem++) {
    iNode2 = iElem*(N_POINTS_QUADRILATERAL+1);
    conn_buf[iNode2+0] = N_POINTS_QUADRILATERAL;
    conn_buf[iNode2+1] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 0)-1;
    conn_buf[iNode2+2] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 1)-1;
    conn_buf[iNode2+3] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 2)-1;
    conn_buf[iNode2+4] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 3)-1;
  }
  if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                            nTot_Quad*(N_POINTS_QUADRILATERAL+1));
  fwrite(conn_buf, sizeof(int),
         nTot_Quad*(N_POINTS_QUADRILATERAL+1), fhw);
  file_size += sizeof(int)*nTot_Quad*(N_POINTS_QUADRILATERAL+1);
  
  for (iElem = 0; iElem < nTot_Tetr; iElem++) {
    iNode2 = iElem*(N_POINTS_TETRAHEDRON+1);
    conn_buf[iNode2+0] = N_POINTS_TETRAHEDRON;
    conn_buf[iNode2+1] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 0)-1;
    conn_buf[iNode2+2] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 1)-1;
    conn_buf[iNode2+3] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2)-1;
    conn_buf[iNode2+4] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3)-1;
  }
  if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                            nTot_Tetr*(N_POINTS_TETRAHEDRON+1));
  fwrite(conn_buf, sizeof(int),
         nTot_Tetr*(N_POINTS_TETRAHEDRON+1), fhw);
  file_size += sizeof(int)*nTot_Tetr*(N_POINTS_TETRAHEDRON+1);
   
  for (iElem = 0; iElem < nTot_Hexa; iElem++) {
    iNode2 = iElem*(N_POINTS_HEXAHEDRON+1);
    conn_buf[iNode2+0] = N_POINTS_HEXAHEDRON;
    conn_buf[iNode2+1] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 0)-1;
    conn_buf[iNode2+2] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 1)-1;
    conn_buf[iNode2+3] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 2)-1;
    conn_buf[iNode2+4] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 3)-1;
    conn_buf[iNode2+5] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 4)-1;
    conn_buf[iNode2+6] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 5)-1;
    conn_buf[iNode2+7] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 6)-1;
    conn_buf[iNode2+8] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 7)-1;
  }
  if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                            nTot_Hexa*(N_POINTS_HEXAHEDRON+1));
  fwrite(conn_buf, sizeof(int),
         nTot_Hexa*(N_POINTS_HEXAHEDRON+1), fhw);
  file_size += sizeof(int)*nTot_Hexa*(N_POINTS_HEXAHEDRON+1);
  
  for (iElem = 0; iElem < nTot_Pris; iElem++) {
    iNode2 = iElem*(N_POINTS_PRISM+1);
    conn_buf[iNode2+0] = N_POINTS_PRISM;
    conn_buf[iNode2+1] = dataSorter->GetElem_Connectivity(PRISM, iElem, 0)-1;
    conn_buf[iNode2+2] = dataSorter->GetElem_Connectivity(PRISM, iElem, 1)-1;
    conn_buf[iNode2+3] = dataSorter->GetElem_Connectivity(PRISM, iElem, 2)-1;
    conn_buf[iNode2+4] = dataSorter->GetElem_Connectivity(PRISM, iElem, 3)-1;
    conn_buf[iNode2+5] = dataSorter->GetElem_Connectivity(PRISM, iElem, 4)-1;
    conn_buf[iNode2+6] = dataSorter->GetElem_Connectivity(PRISM, iElem, 5)-1;
  }
  if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                            nTot_Pris*(N_POINTS_PRISM+1));
  fwrite(conn_buf, sizeof(int),
         nTot_Pris*(N_POINTS_PRISM+1), fhw);
  file_size += sizeof(int)*nTot_Pris*(N_POINTS_PRISM+1);
  
  for (iElem = 0; iElem < nTot_Pyra; iElem++) {
    iNode2 = iElem*(N_POINTS_PYRAMID+1);
    conn_buf[iNode2+0] = N_POINTS_PYRAMID;
    conn_buf[iNode2+1] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 0)-1;
    conn_buf[iNode2+2] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 1)-1;
    conn_buf[iNode2+3] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 2)-1;
    conn_buf[iNode2+4] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 3)-1;
    conn_buf[iNode2+5] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 4)-1;
  }
  if (!BigEndian) SwapBytes((char *)conn_buf, sizeof(int),
                            nTot_Pyra*(N_POINTS_PYRAMID+1));
  fwrite(conn_buf, sizeof(int),
         nTot_Pyra*(N_POINTS_PYRAMID+1), fhw);
  file_size += sizeof(int)*nTot_Pyra*(N_POINTS_PYRAMID+1);
  
  
  if (conn_buf != NULL) delete [] conn_buf;
  
  /*--- Load/write the cell type for all elements in the file. ---*/
  

  SPRINTF (str_buf, "\nCELL_TYPES %i\n", SU2_TYPE::Int(dataSorter->GetnElem()));
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  file_size += sizeof(char)*strlen(str_buf);
  
  int *type_buf = NULL;
  
  type_buf = new int[dataSorter->GetnElem()];
  
  for (iElem = 0; iElem < nTot_Line; iElem++) {
    type_buf[iElem] = LINE;
  }
  if (!BigEndian)
    SwapBytes((char *)type_buf, sizeof(int), nTot_Line);
  fwrite(type_buf, sizeof(int), nTot_Line, fhw);
  file_size += sizeof(int)*nTot_Line;
    
  for (iElem = 0; iElem < nTot_Tria; iElem++) {
    type_buf[iElem] = TRIANGLE;
  }
  if (!BigEndian)
    SwapBytes((char *)type_buf, sizeof(int), nTot_Tria);
  fwrite(type_buf, sizeof(int), nTot_Tria, fhw);
  file_size += sizeof(int)*nTot_Tria;
  
  for (iElem = 0; iElem < nTot_Quad; iElem++) {
    type_buf[iElem] = QUADRILATERAL;
  }
  if (!BigEndian)
    SwapBytes((char *)type_buf, sizeof(int), nTot_Quad);
  fwrite(type_buf, sizeof(int), nTot_Quad, fhw);
  file_size += sizeof(int)*nTot_Quad;
  
  for (iElem = 0; iElem < nTot_Tetr; iElem++) {
    type_buf[iElem] = TETRAHEDRON;
  }
  if (!BigEndian)
    SwapBytes((char *)type_buf, sizeof(int), nTot_Tetr);
  fwrite(type_buf, sizeof(int), nTot_Tetr, fhw);
  file_size += sizeof(int)*nTot_Tetr;
  
  for (iElem = 0; iElem < nTot_Hexa; iElem++) {
    type_buf[iElem] = HEXAHEDRON;
  }
  if (!BigEndian)
    SwapBytes((char *)type_buf, sizeof(int), nTot_Hexa);
  fwrite(type_buf, sizeof(int), nTot_Hexa, fhw);
  file_size += sizeof(int)*nTot_Hexa;
  
  for (iElem = 0; iElem < nTot_Pris; iElem++) {
    type_buf[iElem] = PRISM;
  }
  if (!BigEndian)
    SwapBytes((char *)type_buf, sizeof(int), nTot_Pris);
  fwrite(type_buf, sizeof(int), nTot_Pris, fhw);
  file_size += sizeof(int)*nTot_Pris;
  
  for (iElem = 0; iElem < nTot_Pyra; iElem++) {
    type_buf[iElem] = PYRAMID;
  }
  if (!BigEndian)
    SwapBytes((char *)type_buf, sizeof(int), nTot_Pyra);
  fwrite(type_buf, sizeof(int), nTot_Pyra, fhw);
  file_size += sizeof(int)*nTot_Pyra;
  
  
  if (type_buf != NULL) delete [] type_buf;
  
  /*--- Now write the scalar and vector data (reuse the counts above). ---*/
  
  SPRINTF (str_buf, "\nPOINT_DATA %i\n", (int)GlobalPoint);
  fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
  file_size += sizeof(char)* strlen(str_buf);
  
  unsigned short varStart = 2;
  if (nDim == 3) varStart++;
  
  /*--- Need to adjust container location to avoid PointID tag and coords. ---*/
  
  unsigned short iField, VarCounter = varStart;
  for (iField = varStart; iField < fieldnames.size(); iField++) {
    
    string fieldname = fieldnames[iField];
    fieldname.erase(remove(fieldname.begin(), fieldname.end(), '"'),
                    fieldname.end());
    
    bool output_variable = true, isVector = false;
    size_t found = fieldnames[iField].find("_x");
    if (found!=string::npos) {
      output_variable = true;
      isVector = true;
    }
    found = fieldnames[iField].find("_y");
    if (found!=string::npos) {
      //skip
      output_variable = false;
      VarCounter++;
    }
    found = fieldnames[iField].find("_z");
    if (found!=string::npos) {
      //skip
      output_variable = false;
      VarCounter++;
    }
    
    if (output_variable && isVector) {
      
      fieldname.erase(fieldname.end()-2,fieldname.end());
      SPRINTF (str_buf, "\nVECTORS %s float\n", fieldname.c_str());
      fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
      file_size += sizeof(char)* strlen(str_buf);
      
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
            val = (float)dataSorter->GetData(VarCounter+iDim,iPoint);
            vec_buf[iPoint*NCOORDS + iDim] = val;
          }
        }
      if (!BigEndian)
        SwapBytes((char *)vec_buf, sizeof(float), NCOORDS*GlobalPoint);
      fwrite(vec_buf, sizeof(float), NCOORDS*GlobalPoint, fhw);
      file_size += sizeof(float)*NCOORDS*GlobalPoint;
      
      delete [] vec_buf;
      
      VarCounter++;
      
    } else if (output_variable) {
      
      SPRINTF (str_buf, "\nSCALARS %s float 1\n", fieldname.c_str());
      fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
      file_size += sizeof(char)* strlen(str_buf);
      
      SPRINTF (str_buf, "LOOKUP_TABLE default\n");
      fwrite(str_buf, sizeof(char), strlen(str_buf), fhw);
      file_size += sizeof(char)* strlen(str_buf);
      
      /*--- Prepare the 1D data buffer on this rank. ---*/
      
      float *scalar_buf = new float[GlobalPoint];
      
      /*--- For now, create a temp 1D buffer to load up the data for writing.
       This will be replaced with a derived data type most likely. ---*/
      
      for (iPoint = 0; iPoint < GlobalPoint; iPoint++) {
        float val = (float)dataSorter->GetData(VarCounter,iPoint);
        scalar_buf[iPoint] = val;
      }
      if (!BigEndian)
        SwapBytes((char *)scalar_buf, sizeof(float), GlobalPoint);
      fwrite(scalar_buf, sizeof(float), GlobalPoint, fhw);
      file_size += sizeof(float)*GlobalPoint;
      
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
    
  /*--- Write the initial strings to the file. Only the master will
   write the header lines, but all ranks will store the offsets. ---*/
  
  disp = 0;
  strcpy(str_buf, "# vtk DataFile Version 3.0\n");
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  file_size += sizeof(char)*strlen(str_buf);
  
  
  strcpy(str_buf, "vtk output\n");
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  file_size += sizeof(char)*strlen(str_buf);
  
  strcpy(str_buf, "BINARY\n");
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  file_size += sizeof(char)*strlen(str_buf);
  
  strcpy(str_buf, "DATASET UNSTRUCTURED_GRID\n");
  if (rank == MASTER_NODE)
    MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  disp += strlen(str_buf)*sizeof(char);
  file_size += sizeof(char)*strlen(str_buf);
  
  /*--- Communicate the number of total points that will be
   written by each rank. After this communication, each proc knows how
   many poinnts will be written before its location in the file and the
   offsets can be correctly set. ---*/
  
  unsigned long myPoint, GlobalPoint;
  
  GlobalPoint = dataSorter->GetnPointsGlobal();
  myPoint     = dataSorter->GetnPoints();
  
  
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
  file_size += sizeof(char)*strlen(str_buf);
  
  /*--- Load/write the 1D buffer of point coordinates. Note that we
   always have 3 coordinate dimensions, even for 2D problems. ---*/
  
  float *coord_buf = new float[myPoint*NCOORDS];
  for (iPoint = 0; iPoint < myPoint; iPoint++) {
    for (iDim = 0; iDim < NCOORDS; iDim++) {
      if (nDim == 2 && iDim == 2) {
        coord_buf[iPoint*NCOORDS + iDim] = 0.0;
      } else {
        float val = (float)dataSorter->GetData(iDim, iPoint);
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
  file_size += sizeof(float)*myPoint*NCOORDS;
  
  /*--- Free the derived datatype and coordinate array. ---*/
  
  MPI_Type_free(&filetype);
  delete [] coord_buf;
  
  /*--- Compute our local number of elements, the required storage,
   and reduce the total number of elements and storage globally. ---*/
  
  unsigned long nTot_Line;
  unsigned long nTot_Tria, nTot_Quad;
  unsigned long nTot_Tetr, nTot_Hexa, nTot_Pris, nTot_Pyra;
  unsigned long myElem, myElemStorage, GlobalElem, GlobalElemStorage;
  
  unsigned long nParallel_Line = dataSorter->GetnElem(LINE),
                nParallel_Tria = dataSorter->GetnElem(TRIANGLE),
                nParallel_Quad = dataSorter->GetnElem(QUADRILATERAL),
                nParallel_Tetr = dataSorter->GetnElem(TETRAHEDRON),
                nParallel_Hexa = dataSorter->GetnElem(HEXAHEDRON),
                nParallel_Pris = dataSorter->GetnElem(PRISM),
                nParallel_Pyra = dataSorter->GetnElem(PYRAMID);
  
  SU2_MPI::Allreduce(&nParallel_Line, &nTot_Line, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
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
  
  myElem        = (nParallel_Line + nParallel_Tria + nParallel_Quad + nParallel_Tetr +
                   nParallel_Hexa + nParallel_Pris + nParallel_Pyra);
  myElemStorage = (nParallel_Line*3 + nParallel_Tria*4 + nParallel_Quad*5 + nParallel_Tetr*5 +
                   nParallel_Hexa*9 + nParallel_Pris*7 + nParallel_Pyra*6);
  
  GlobalElem        = (nTot_Line + nTot_Tria   + nTot_Quad   + nTot_Tetr   +
                       nTot_Hexa   + nTot_Pris   + nTot_Pyra);
  GlobalElemStorage = (nTot_Line*3 + nTot_Tria*4 + nTot_Quad*5 + nTot_Tetr*5 +
                       nTot_Hexa*9 + nTot_Pris*7 + nTot_Pyra*6);
  
  
  
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
  file_size += sizeof(char)*strlen(str_buf);
  
  /*--- Load/write 1D buffers for the connectivity of each element type. ---*/
  
  int *conn_buf = new int[myElemStorage];
  unsigned long iStorage = 0;
  
  for (iElem = 0; iElem < nParallel_Line; iElem++) {
    conn_buf[iStorage+0] = N_POINTS_LINE;
    conn_buf[iStorage+1] = dataSorter->GetElem_Connectivity(LINE, iElem, 0)-1;
    conn_buf[iStorage+2] = dataSorter->GetElem_Connectivity(LINE, iElem, 1)-1;
    iStorage += (N_POINTS_LINE+1);
  }
  
  for (iElem = 0; iElem < nParallel_Tria; iElem++) {
    conn_buf[iStorage+0] = N_POINTS_TRIANGLE;
    conn_buf[iStorage+1] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 0)-1;
    conn_buf[iStorage+2] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 1)-1;
    conn_buf[iStorage+3] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 2)-1 ;
    iStorage += (N_POINTS_TRIANGLE+1);
  }
  
  for (iElem = 0; iElem < nParallel_Quad; iElem++) {
    conn_buf[iStorage+0] = N_POINTS_QUADRILATERAL;
    conn_buf[iStorage+1] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 0)-1;
    conn_buf[iStorage+2] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 1)-1;
    conn_buf[iStorage+3] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 2)-1;
    conn_buf[iStorage+4] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 3)-1;
    iStorage += (N_POINTS_QUADRILATERAL+1);
    
  }
  
  for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
    conn_buf[iStorage+0] = N_POINTS_TETRAHEDRON;
    conn_buf[iStorage+1] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 0)-1;
    conn_buf[iStorage+2] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 1)-1;
    conn_buf[iStorage+3] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2)-1;
    conn_buf[iStorage+4] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3)-1;
    iStorage += (N_POINTS_TETRAHEDRON+1);
    
  }
  
  for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
    conn_buf[iStorage+0] = N_POINTS_HEXAHEDRON;
    conn_buf[iStorage+1] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 0)-1;
    conn_buf[iStorage+2] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 1)-1;
    conn_buf[iStorage+3] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 2)-1;
    conn_buf[iStorage+4] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 3)-1;
    conn_buf[iStorage+5] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 4)-1;
    conn_buf[iStorage+6] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 5)-1;
    conn_buf[iStorage+7] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 6)-1;
    conn_buf[iStorage+8] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 7)-1;
    iStorage += (N_POINTS_HEXAHEDRON+1);
  }
  
  for (iElem = 0; iElem < nParallel_Pris; iElem++) {
    conn_buf[iStorage+0] = N_POINTS_PRISM;
    conn_buf[iStorage+1] = dataSorter->GetElem_Connectivity(PRISM, iElem, 0)-1;
    conn_buf[iStorage+2] = dataSorter->GetElem_Connectivity(PRISM, iElem, 1)-1;
    conn_buf[iStorage+3] = dataSorter->GetElem_Connectivity(PRISM, iElem, 2)-1;
    conn_buf[iStorage+4] = dataSorter->GetElem_Connectivity(PRISM, iElem, 3)-1;
    conn_buf[iStorage+5] = dataSorter->GetElem_Connectivity(PRISM, iElem, 4)-1;
    conn_buf[iStorage+6] = dataSorter->GetElem_Connectivity(PRISM, iElem, 5)-1;
    iStorage += (N_POINTS_PRISM+1);
  }
  
  for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
    conn_buf[iStorage+0] = N_POINTS_PYRAMID;
    conn_buf[iStorage+1] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 0)-1;
    conn_buf[iStorage+2] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 1)-1;
    conn_buf[iStorage+3] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 2)-1;
    conn_buf[iStorage+4] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 3)-1;
    conn_buf[iStorage+5] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 4)-1;
    iStorage += (N_POINTS_PYRAMID+1);
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
  
  file_size += sizeof(int)*myElemStorage;  
  
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
  file_size += sizeof(char)*strlen(str_buf);  
  
  int *type_buf = new int[myElem];
  unsigned long jElem = 0;
  
  for (iElem = 0; iElem < nParallel_Line; iElem++) {
    type_buf[jElem] = LINE; jElem++;
  }
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
  
  file_size += sizeof(int)*myElem;    
  
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
  file_size += sizeof(char)*strlen(str_buf);  
  
  /*--- Adjust container start location to avoid point coords. ---*/
  
  unsigned short varStart = 2;
  if (nDim == 3) varStart++;
  
  /*--- Loop over all variables that have been registered in the output. ---*/
  
  unsigned short iField, VarCounter = varStart;
  for (iField = varStart; iField < fieldnames.size(); iField++) {
    
    string fieldname = fieldnames[iField];
    fieldname.erase(remove(fieldname.begin(), fieldname.end(), '"'),
                    fieldname.end());
    
    /*--- Check whether this field is a vector or scalar. ---*/
    
    bool output_variable = true, isVector = false;
    size_t found = fieldnames[iField].find("_x");
    if (found!=string::npos) {
      output_variable = true;
      isVector        = true;
    }
    found = fieldnames[iField].find("_y");
    if (found!=string::npos) {
      /*--- We have found a vector, so skip the Y component. ---*/
      output_variable = false;
      VarCounter++;
    }
    found = fieldnames[iField].find("_z");
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
      file_size += sizeof(char)*strlen(str_buf);  
      
      /*--- Prepare the 1D data buffer on this rank. ---*/
      
      float *vec_buf = new float[myPoint*NCOORDS];
      
      /*--- Load up the buffer for writing this rank's vector data. ---*/
      
      float val = 0.0;
      for (iPoint = 0; iPoint < myPoint; iPoint++) {
        for (iDim = 0; iDim < NCOORDS; iDim++) {
          if (nDim == 2 && iDim == 2) {
            vec_buf[iPoint*NCOORDS + iDim] = 0.0;
          } else {
            val = (float)dataSorter->GetData(VarCounter+iDim,iPoint);
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
      
      file_size += sizeof(float)*myPoint*NCOORDS;  
      
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
      file_size += sizeof(char)*strlen(str_buf);  
      
      MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE,
                        (char*)"native", MPI_INFO_NULL);
      SPRINTF (str_buf, "LOOKUP_TABLE default\n");
      if (rank == MASTER_NODE)
        MPI_File_write_at(fhw, disp, str_buf, strlen(str_buf),
                          MPI_CHAR, MPI_STATUS_IGNORE);
      disp += strlen(str_buf)*sizeof(char);
      file_size += sizeof(char)*strlen(str_buf);  
      
      /*--- Prepare the 1D data buffer on this rank. ---*/
      
      float *scalar_buf = new float[myPoint];
      
      /*--- For now, create a temp 1D buffer to load up the data for writing.
       This will be replaced with a derived data type most likely. ---*/
      
      for (iPoint = 0; iPoint < myPoint; iPoint++) {
        float val = (float)dataSorter->GetData(VarCounter,iPoint);
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
      
      file_size += sizeof(float)*myPoint;  
      
      /*--- Free the derived datatype and coordinate array. ---*/
      
      MPI_Type_free(&filetype);
      delete [] scalar_buf; scalar_buf = NULL;
      
      VarCounter++;
    }
    
  }
  
  /*--- All ranks close the file after writing. ---*/
  
  MPI_File_close(&fhw);
  
  
  /*--- Compute and store the write time. ---*/
  
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  UsedTime = StopTime-StartTime;
  
  /*--- Communicate the total file size for the restart ---*/
  
#ifdef HAVE_MPI
  su2double my_file_size = file_size;
  SU2_MPI::Allreduce(&my_file_size, &file_size, 1,
                     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Compute and store the bandwidth ---*/
  
  Bandwidth = file_size/(1.0e6)/UsedTime;
  
  /*--- Delete the offset counters that we needed for MPI IO. ---*/
  
  delete [] nElem_Snd;        delete [] nElem_Cum;
  delete [] nElemStorage_Snd; delete [] nElemStorage_Cum;
  delete [] nPoint_Snd;       delete [] nPoint_Cum;
  
#endif
}


/*--- Subroutine to swap bytes, in case we need to convert to
 big endian, which is expected for ParaView binary legacy format. ---*/

void CParaviewBinaryFileWriter::SwapBytes(char *buffer, size_t nBytes, unsigned long nVar) {
  
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
