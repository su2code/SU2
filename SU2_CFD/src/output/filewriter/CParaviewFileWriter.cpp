#include "../../../include/output/filewriter/CParaviewFileWriter.hpp"

const string CParaviewFileWriter::fileExt = ".vtk";

CParaviewFileWriter::CParaviewFileWriter(vector<string> fields, unsigned short nDim, 
                                         string fileName, CParallelDataSorter *dataSorter) : 
  CFileWriter(std::move(fields), std::move(fileName), dataSorter, fileExt, nDim){}


CParaviewFileWriter::~CParaviewFileWriter(){}

void CParaviewFileWriter::Write_Data(){
    
  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }
  
  unsigned short iDim;

  unsigned long iPoint, iElem;

  unsigned long nGlobal_Elem_Storage;

  ofstream Paraview_File;

  int iProcessor;
  
  /*--- Set a timer for the file writing. ---*/
  
#ifndef HAVE_MPI
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StartTime = MPI_Wtime();
#endif
  
  /*--- Open Paraview ASCII file and write the header. ---*/
  
  if (rank == MASTER_NODE) {
    Paraview_File.open(fileName.c_str(), ios::out);
    Paraview_File.precision(6);
    Paraview_File << "# vtk DataFile Version 3.0\n";
    Paraview_File << "vtk output\n";
    Paraview_File << "ASCII\n";
    Paraview_File << "DATASET UNSTRUCTURED_GRID\n";
    
    /*--- Write the header ---*/
    Paraview_File << "POINTS "<< dataSorter->GetnPointsGlobal() <<" double\n";
    
  }

  Paraview_File.close();

#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  /*--- Each processor opens the file. ---*/

  Paraview_File.open(fileName.c_str(), ios::out | ios::app);

  /*--- Write surface and volumetric point coordinates. ---*/

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      
      /*--- Write the node data from this proc ---*/
      
      
      for (iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++)
          Paraview_File << scientific << dataSorter->GetData(iDim, iPoint) << "\t";
        if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
      }
    }
    
    Paraview_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }

  /*--- Reduce the total number of each element. ---*/

  unsigned long nTot_Line, nTot_Tria, nTot_Quad, nTot_Tetr, nTot_Hexa, nTot_Pris, nTot_Pyra;
  unsigned long nParallel_Line = dataSorter->GetnElem(LINE),
                nParallel_Tria = dataSorter->GetnElem(TRIANGLE),
                nParallel_Quad = dataSorter->GetnElem(QUADRILATERAL),
                nParallel_Tetr = dataSorter->GetnElem(TETRAHEDRON),
                nParallel_Hexa = dataSorter->GetnElem(HEXAHEDRON),
                nParallel_Pris = dataSorter->GetnElem(PRISM),
                nParallel_Pyra = dataSorter->GetnElem(PYRAMID);
#ifdef HAVE_MPI
  SU2_MPI::Reduce(&nParallel_Line, &nTot_Line, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_Tria, &nTot_Tria, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_Quad, &nTot_Quad, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_Tetr, &nTot_Tetr, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_Hexa, &nTot_Hexa, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_Pris, &nTot_Pris, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&nParallel_Pyra, &nTot_Pyra, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
  nTot_Line      = nParallel_Line;

  nTot_Tria = nParallel_Tria;
  nTot_Quad = nParallel_Quad;
  nTot_Tetr = nParallel_Tetr;
  nTot_Hexa = nParallel_Hexa;
  nTot_Pris = nParallel_Pris;
  nTot_Pyra = nParallel_Pyra;
#endif
  
  if (rank == MASTER_NODE) {
    
    /*--- Write the header ---*/
    nGlobal_Elem_Storage = nTot_Line*3 + nTot_Tria*4 + nTot_Quad*5 + nTot_Tetr*5 + nTot_Hexa*9 + nTot_Pris*7 + nTot_Pyra*6;
    
    Paraview_File << "\nCELLS " << dataSorter->GetnElem() << "\t" << nGlobal_Elem_Storage << "\n";
    
  }

  Paraview_File.flush();
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- Write connectivity data. ---*/
  
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      
      
      for (iElem = 0; iElem < nParallel_Line; iElem++) {
        Paraview_File << N_POINTS_LINE << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(LINE, iElem, 0)-1 << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(LINE, iElem, 1)-1 << "\t";
      }
      
      for (iElem = 0; iElem < nParallel_Tria; iElem++) {
        Paraview_File << N_POINTS_TRIANGLE << "\t";
        Paraview_File <<  dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 0)-1 << "\t";
        Paraview_File <<  dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 1)-1 << "\t";
        Paraview_File <<  dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 2)-1 << "\t";
      }
      
      for (iElem = 0; iElem < nParallel_Quad; iElem++) {
        Paraview_File << N_POINTS_QUADRILATERAL << "\t";
        Paraview_File <<  dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 0)-1 << "\t";
        Paraview_File <<  dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 1)-1 << "\t";
        Paraview_File <<  dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 2)-1 << "\t";
        Paraview_File <<  dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 3)-1 << "\t";
      }
      
      
      for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
        Paraview_File << N_POINTS_TETRAHEDRON << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 0)-1 << "\t" 
                      << dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 1)-1 << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2)-1 << "\t" 
                      << dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3)-1 << "\t";
      }
      
      for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
        Paraview_File << N_POINTS_HEXAHEDRON << "\t"; 
        Paraview_File << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 0)-1 << "\t" 
                      << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 1)-1 << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 2)-1 << "\t" 
                      << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 3)-1 << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 4)-1 << "\t" 
                      << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 5)-1 << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 6)-1 << "\t" 
                      << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 7)-1 << "\t";
      }
      
      for (iElem = 0; iElem < nParallel_Pris; iElem++) {
        Paraview_File << N_POINTS_PRISM << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(PRISM, iElem, 0)-1 << "\t" 
                      << dataSorter->GetElem_Connectivity(PRISM, iElem, 1)-1 << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(PRISM, iElem, 2)-1 << "\t" 
                      << dataSorter->GetElem_Connectivity(PRISM, iElem, 3)-1 << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(PRISM, iElem, 4)-1 << "\t" 
                      << dataSorter->GetElem_Connectivity(PRISM, iElem, 5)-1 << "\t";
      }
      
      for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
        Paraview_File << N_POINTS_PYRAMID << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(PYRAMID, iElem, 0)-1 << "\t" 
                      << dataSorter->GetElem_Connectivity(PYRAMID, iElem, 1)-1 << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(PYRAMID, iElem, 2)-1 << "\t" 
                      << dataSorter->GetElem_Connectivity(PYRAMID, iElem, 3)-1 << "\t";
        Paraview_File << dataSorter->GetElem_Connectivity(PYRAMID, iElem, 4)-1 << "\t";
      }
      
    }    Paraview_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }
  
  if (rank == MASTER_NODE) {
    
    /*--- Write the header ---*/
    Paraview_File << "\nCELL_TYPES " << dataSorter->GetnElem() << "\n";
    
  }

  Paraview_File.flush();
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iElem = 0; iElem < nParallel_Line; iElem++) Paraview_File << "3\t";
      for (iElem = 0; iElem < nParallel_Tria; iElem++) Paraview_File << "5\t";
      for (iElem = 0; iElem < nParallel_Quad; iElem++) Paraview_File << "9\t";
      for (iElem = 0; iElem < nParallel_Tetr; iElem++) Paraview_File << "10\t";
      for (iElem = 0; iElem < nParallel_Hexa; iElem++) Paraview_File << "12\t";
      for (iElem = 0; iElem < nParallel_Pris; iElem++) Paraview_File << "13\t";
      for (iElem = 0; iElem < nParallel_Pyra; iElem++) Paraview_File << "14\t";  
    }   
    Paraview_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }
  
  if (rank == MASTER_NODE) {
    /*--- Write the header ---*/
    Paraview_File << "\nPOINT_DATA "<< dataSorter->GetnPointsGlobal() <<"\n";
    
  }

  Paraview_File.flush();
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  unsigned short varStart = 2;
  if (nDim == 3) varStart++;

  /*--- Need to adjust container location to avoid PointID tag and coords. ---*/
  unsigned short VarCounter = varStart;

  for (unsigned short iField = varStart; iField < fieldnames.size(); iField++) {

    string fieldname = fieldnames[iField];

    fieldname.erase(remove(fieldname.begin(), fieldname.end(), '"'), fieldname.end());

    bool output_variable = true, isVector = false;
    size_t found = fieldnames[iField].find("_x");
    if (found!=string::npos) {
      output_variable = true;
      isVector = true;
    }
    found = fieldnames[iField].find("_y");
    if (found!=string::npos) {
      output_variable = false;
      //skip
      Paraview_File.flush();
#ifdef HAVE_MPI
      SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
      VarCounter++;
    }
found = fieldnames[iField].find("_z");
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
          
          for (iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
            Paraview_File << scientific << dataSorter->GetData(VarCounter+0, iPoint) << "\t" << dataSorter->GetData(VarCounter+1, iPoint) << "\t";
            if (nDim == 3) Paraview_File << scientific << dataSorter->GetData(VarCounter+2, iPoint) << "\t";
            if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
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
          
          for (iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
            Paraview_File << scientific << dataSorter->GetData(VarCounter, iPoint) << "\t";
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
  
  
  /*--- Compute and store the write time. ---*/
  
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  UsedTime = StopTime-StartTime;
  
  file_size = Determine_Filesize(fileName);
  
  /*--- Compute and store the bandwidth ---*/
  
  Bandwidth = file_size/(1.0e6)/UsedTime;
}

