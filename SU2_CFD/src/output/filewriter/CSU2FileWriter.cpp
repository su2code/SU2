#include "../../../include/output/filewriter/CSU2FileWriter.hpp"




CSU2FileWriter::CSU2FileWriter(vector<string> fields, unsigned short nDim) : 
  CFileWriter(fields, nDim){

  file_ext = ".dat";
    
}


CSU2FileWriter::~CSU2FileWriter(){
  
}

void CSU2FileWriter::Write_Data(string filename, CParallelDataSorter *data_sorter){
  
  filename += file_ext;
  
  /*--- Local variables ---*/
  
  unsigned short iVar;
  unsigned long iPoint;

  ofstream restart_file;
  
  int iProcessor;
  
  /*--- Set a timer for the file writing. ---*/
  
#ifndef HAVE_MPI
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StartTime = MPI_Wtime();
#endif
  
  /*--- Only the master node writes the header. ---*/
  
  if (rank == MASTER_NODE) {
    restart_file.open(filename.c_str(), ios::out);
    restart_file.precision(15);
    restart_file << "\"PointID\"";
    for (iVar = 0; iVar < fieldnames.size()-1; iVar++)
      restart_file << "\t\"" << fieldnames[iVar] << "\"";
    restart_file << "\t\"" << fieldnames[fieldnames.size()-1] << "\"" << endl;
    restart_file.close();
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- All processors open the file. ---*/
  
  restart_file.open(filename.c_str(), ios::out | ios::app);
  restart_file.precision(15);
  
  /*--- Write the restart file in parallel, processor by processor. ---*/
  
  unsigned long myPoint = 0, offset = 0, Global_Index;
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint < data_sorter->GetnPoints(); iPoint++) {
        
        /*--- Global Index of the current point. (note outer loop over procs) ---*/
        
        Global_Index = data_sorter->GetGlobalIndex(iPoint);
        
        /*--- Write global index. (note outer loop over procs) ---*/
        
        restart_file << Global_Index << "\t";
        myPoint++;
        
        /*--- Loop over the variables and write the values to file ---*/
        
        for (iVar = 0; iVar < fieldnames.size(); iVar++) {
          restart_file << scientific << data_sorter->GetData(iVar, iPoint) << "\t";
        }
        restart_file << "\n";
      }
      
    }
    /*--- Flush the file and wait for all processors to arrive. ---*/
    restart_file.flush();
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
    
  }
  
  /*--- Compute and store the write time. ---*/
  
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  UsedTime = StopTime-StartTime;
  
  /*--- Determine the file size ---*/
  
  file_size = Determine_Filesize(filename);
  
  /*--- Compute and store the bandwidth ---*/
  
  Bandwidth = file_size/(1.0e6)/UsedTime;
  
  /*--- Write the metadata (master rank alone) ----*/

//  if (rank == MASTER_NODE) {
//    if (dual_time)
//      restart_file <<"EXT_ITER= " << config->GetExtIter() + 1 << endl;
//    else
//      restart_file <<"EXT_ITER= " << config->GetExtIter() + config->GetExtIter_OffSet() + 1 << endl;
//    restart_file <<"AOA= " << config->GetAoA() - config->GetAoA_Offset() << endl;
//    restart_file <<"SIDESLIP_ANGLE= " << config->GetAoS() - config->GetAoS_Offset() << endl;
//    restart_file <<"INITIAL_BCTHRUST= " << config->GetInitial_BCThrust() << endl;
//    restart_file <<"DCD_DCL_VALUE= " << config->GetdCD_dCL() << endl;
//    restart_file <<"DCMX_DCL_VALUE= " << config->GetdCMx_dCL() << endl;
//    restart_file <<"DCMY_DCL_VALUE= " << config->GetdCMy_dCL() << endl;
//    restart_file <<"DCMZ_DCL_VALUE= " << config->GetdCMz_dCL() << endl;

//    if (( config->GetKind_Solver() == DISC_ADJ_EULER ||
//          config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES ||
//          config->GetKind_Solver() == DISC_ADJ_RANS ) && adjoint) {
//      restart_file << "SENS_AOA=" << solver[ADJFLOW_SOL]->GetTotal_Sens_AoA() * PI_NUMBER / 180.0 << endl;
//    }
//  }

  /*--- All processors close the file. ---*/

  restart_file.close();
}
