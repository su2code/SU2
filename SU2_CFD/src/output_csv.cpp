#include "../include/output_structure.hpp"

void COutput::WriteSurface_CSV(CConfig *config, CGeometry *geometry){
  
  unsigned short iVar;  
  unsigned long iPoint, index, offset, myPoint;
  ofstream Surf_file;
  
  int iProcessor;
  
  string filename;
  
  filename = config->GetFilename(SurfaceFilename, ".csv");
  
  Surf_file.precision(15);
  
  
  if (rank == MASTER_NODE) {
    Surf_file.open(filename.c_str(), ios::out);
    Surf_file << "\"Point\",";
    
    for (iVar = 0; iVar < Variable_Names.size()-1; iVar++) {
      Surf_file << "\"" << Variable_Names[iVar] << "\",";
    }
    Surf_file << "\"" << Variable_Names[Variable_Names.size()-1] << "\"" << endl;
    Surf_file.close();
  }
  
  /*--- Open file on every processor ---*/
  
  Surf_file.open(filename.c_str(), ios::out| ios::app);
  
  /*--- Write surface and volumetric solution data. ---*/
  
  offset = 0; myPoint = 0;
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint < nSurf_Poin_Par; iPoint++) {
        
        index = iPoint + offset;
        
        Surf_file << Renumber2Global[index+1] << ", ";
        
        for (iVar = 0; iVar < GlobalField_Counter; iVar++){
          Surf_file << scientific << Parallel_Surf_Data[iVar][iPoint];
          if (iVar != GlobalField_Counter -1) Surf_file << ", ";
        }
        Surf_file << endl;
        
        myPoint++;
      }
    }
    Surf_file.flush();
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);    
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }
}