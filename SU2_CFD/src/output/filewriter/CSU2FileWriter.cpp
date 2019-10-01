#include "../../../include/output/filewriter/CSU2FileWriter.hpp"

const string CSU2FileWriter::fileExt = ".csv";

CSU2FileWriter::CSU2FileWriter(vector<string> fields, unsigned short nDim, 
                               string fileName, CParallelDataSorter *dataSorter) : 
  CFileWriter(std::move(fields), std::move(fileName), dataSorter, fileExt, nDim){}


CSU2FileWriter::~CSU2FileWriter(){
  
}

void CSU2FileWriter::Write_Data(){
    
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
    restart_file.open(fileName.c_str(), ios::out);
    restart_file.precision(15);
    restart_file << "\"PointID\"";
    for (iVar = 0; iVar < fieldnames.size()-1; iVar++)
      restart_file << ",\"" << fieldnames[iVar] << "\"";
    restart_file << ",\"" << fieldnames[fieldnames.size()-1] << "\"" << endl;
    restart_file.close();
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- All processors open the file. ---*/
  
  restart_file.open(fileName.c_str(), ios::out | ios::app);
  restart_file.precision(15);
  
  /*--- Write the restart file in parallel, processor by processor. ---*/
  
  unsigned long myPoint = 0, Global_Index;
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
        
        /*--- Global Index of the current point. (note outer loop over procs) ---*/
        
        Global_Index = dataSorter->GetGlobalIndex(iPoint);
        
        /*--- Write global index. (note outer loop over procs) ---*/
        
        restart_file << Global_Index << ", ";
        myPoint++;
        
        /*--- Loop over the variables and write the values to file ---*/
        
        for (iVar = 0; iVar < fieldnames.size()-1; iVar++) {
          restart_file << scientific << dataSorter->GetData(iVar, iPoint) << ", ";
        }
        restart_file << scientific << dataSorter->GetData(fieldnames.size()-1, iPoint) << "\n";        
      }
      
    }
    /*--- Flush the file and wait for all processors to arrive. ---*/
    restart_file.flush();
#ifdef HAVE_MPI
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
  
  file_size = Determine_Filesize(fileName);
  
  /*--- Compute and store the bandwidth ---*/
  
  Bandwidth = file_size/(1.0e6)/UsedTime;
  
  /*--- All processors close the file. ---*/

  restart_file.close();
}
