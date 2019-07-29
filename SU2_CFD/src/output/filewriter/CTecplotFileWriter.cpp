#include "../../../include/output/filewriter/CTecplotFileWriter.hpp"

CTecplotFileWriter::CTecplotFileWriter(vector<string> fields, unsigned short nDim, unsigned long time_iter, su2double timestep) : 
  CFileWriter(fields, nDim){

  file_ext = ".dat";
    
  this->time_iter = time_iter;
  
  this->timestep = timestep;
  
}


CTecplotFileWriter::~CTecplotFileWriter(){
  
}

void CTecplotFileWriter::Write_Data(string filename, CParallelDataSorter *data_sorter){
  
  filename += file_ext;
  
  if (!data_sorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }
  
  unsigned short iVar;
  
  unsigned long iPoint, iElem;
 
  int iProcessor;

  ofstream Tecplot_File;
  
  file_size = 0.0;
  
  /*--- Set a timer for the file writing. ---*/
  
#ifndef HAVE_MPI
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StartTime = MPI_Wtime();
#endif
  
  /*--- Reduce the total number of each element. ---*/

  unsigned long nTot_Line, nTot_Tria, nTot_Quad, nTot_Tetr, nTot_Hexa, nTot_Pris, nTot_Pyra;
  unsigned long nParallel_Line = data_sorter->GetnElem(LINE),
                nParallel_Tria = data_sorter->GetnElem(TRIANGLE),
                nParallel_Quad = data_sorter->GetnElem(QUADRILATERAL),
                nParallel_Tetr = data_sorter->GetnElem(TETRAHEDRON),
                nParallel_Hexa = data_sorter->GetnElem(HEXAHEDRON),
                nParallel_Pris = data_sorter->GetnElem(PRISM),
                nParallel_Pyra = data_sorter->GetnElem(PYRAMID);
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nParallel_Line, &nTot_Line, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Tria, &nTot_Tria, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Quad, &nTot_Quad, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Tetr, &nTot_Tetr, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Hexa, &nTot_Hexa, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Pris, &nTot_Pris, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Pyra, &nTot_Pyra, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTot_Line      = nParallel_Line;

  nTot_Tria = nParallel_Tria;
  nTot_Quad = nParallel_Quad;
  nTot_Tetr = nParallel_Tetr;
  nTot_Hexa = nParallel_Hexa;
  nTot_Pris = nParallel_Pris;
  nTot_Pyra = nParallel_Pyra;
#endif
    
  /*--- Open Tecplot ASCII file and write the header. ---*/
  
  if (rank == MASTER_NODE) {
    Tecplot_File.open(filename.c_str(), ios::out);
    Tecplot_File.precision(6);
    Tecplot_File << "TITLE = \"Visualization of the solution\"" << endl;
    
    Tecplot_File << "VARIABLES = ";
    for (iVar = 0; iVar < fieldnames.size()-1; iVar++) {
      Tecplot_File << "\"" << fieldnames[iVar] << "\",";
    }
    Tecplot_File << "\"" << fieldnames[fieldnames.size()-1] << "\"" << endl;
    
    /*--- Write the header ---*/
    
    Tecplot_File << "ZONE ";

    if (timestep > 0.0){
      Tecplot_File << "STRANDID="<<SU2_TYPE::Int(time_iter+1)<<", SOLUTIONTIME="<< time_iter*timestep <<", ";      
    }
    
//    if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
//      Tecplot_File << "STRANDID="<<SU2_TYPE::Int(iExtIter+1)<<", SOLUTIONTIME="<<config->GetDelta_UnstTime()*iExtIter<<", ";
//    } else if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
//      /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
//      su2double period = config->GetHarmonicBalance_Period();
//      su2double deltaT = period/(su2double)(config->GetnTimeInstances());
//      Tecplot_File << "STRANDID="<<SU2_TYPE::Int(val_iZone+1)<<", SOLUTIONTIME="<<deltaT*val_iZone<<", ";
//    }
//    if (nDim == 2) {
//      if (surf_sol) Tecplot_File << "NODES= "<< nGlobal_Surf_Poin <<", ELEMENTS= "<< nSurf_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
//      else Tecplot_File << "NODES= "<< nGlobal_Poin_Par <<", ELEMENTS= "<< nGlobal_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
//    } else {
//      if (surf_sol) Tecplot_File << "NODES= "<< nGlobal_Surf_Poin <<", ELEMENTS= "<< nSurf_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
//      else Tecplot_File << "NODES= "<< nGlobal_Poin_Par <<", ELEMENTS= "<< nGlobal_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
//    }
    
    Tecplot_File << "NODES= "<< data_sorter->GetnPointsGlobal() <<", ELEMENTS= "<< data_sorter->GetnElem();
    
    if (nDim == 3){
      if ((nTot_Quad > 0 || nTot_Tria > 0) && (nTot_Hexa + nTot_Pris + nTot_Pyra + nTot_Tetr == 0)){
        Tecplot_File << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;
      }
      else {
        Tecplot_File <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
      }
    }
    else {
      if (nTot_Line > 0 && (nTot_Tria + nTot_Quad == 0)){
        Tecplot_File << ", DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;       
      }
      else{
        Tecplot_File << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;   
      }
    }
    Tecplot_File.close();    
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- Each processor opens the file. ---*/
  
  Tecplot_File.open(filename.c_str(), ios::out | ios::app);
  
  /*--- Write surface and volumetric solution data. ---*/
  
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      
      /*--- Write the node data from this proc ---*/
      
      
      for (iPoint = 0; iPoint < data_sorter->GetnPoints(); iPoint++) {
        for (iVar = 0; iVar < fieldnames.size(); iVar++)
          Tecplot_File << scientific << data_sorter->GetData(iVar, iPoint) << "\t";
        Tecplot_File << endl;
      }
    }
    
    Tecplot_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }

  
  /*--- Write connectivity data. ---*/
  
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      
      for (iElem = 0; iElem < nParallel_Line; iElem++) {
        Tecplot_File << data_sorter->GetElem_Connectivity(LINE, iElem, 0) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(LINE, iElem, 1)<< "\n";
      }
      
      
      for (iElem = 0; iElem < nParallel_Tria; iElem++) {
        Tecplot_File << data_sorter->GetElem_Connectivity(TRIANGLE, iElem, 0) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(TRIANGLE, iElem, 1) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(TRIANGLE, iElem, 2) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(TRIANGLE, iElem, 2) << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Quad; iElem++) {
        Tecplot_File << data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 0) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 1) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 2) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 3) << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
        Tecplot_File << data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 0) << "\t" << data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 1) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2) << "\t" << data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3) << "\t" << data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3) << "\t" << data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3) << "\n";        
      }
      
      for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
        Tecplot_File << data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 0) << "\t" << data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 1) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 2) << "\t" << data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 3) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 4) << "\t" << data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 5) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 6) << "\t" << data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 7) << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Pris; iElem++) {
        Tecplot_File << data_sorter->GetElem_Connectivity(PRISM, iElem, 0) << "\t" << data_sorter->GetElem_Connectivity(PRISM, iElem, 1) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(PRISM, iElem, 1) << "\t" << data_sorter->GetElem_Connectivity(PRISM, iElem, 2) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(PRISM, iElem, 3) << "\t" << data_sorter->GetElem_Connectivity(PRISM, iElem, 4) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(PRISM, iElem, 4) << "\t" << data_sorter->GetElem_Connectivity(PRISM, iElem, 5) << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
        Tecplot_File << data_sorter->GetElem_Connectivity(PYRAMID, iElem, 0) << "\t" << data_sorter->GetElem_Connectivity(PYRAMID, iElem, 1) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(PYRAMID, iElem, 2) << "\t" << data_sorter->GetElem_Connectivity(PYRAMID, iElem, 3) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(PYRAMID, iElem, 4) << "\t" << data_sorter->GetElem_Connectivity(PYRAMID, iElem, 4) << "\t";
        Tecplot_File << data_sorter->GetElem_Connectivity(PYRAMID, iElem, 4) << "\t" << data_sorter->GetElem_Connectivity(PYRAMID, iElem, 4) << "\n";
      }
      
      
    }
    Tecplot_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }
  
  Tecplot_File.close();
  
  /*--- Compute and store the write time. ---*/
  
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  UsedTime = StopTime-StartTime;
  
  file_size = Determine_Filesize(filename);
  
  /*--- Compute and store the bandwidth ---*/
  
  Bandwidth = file_size/(1.0e6)/UsedTime;
}


