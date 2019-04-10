#include "../include/output_structure.hpp"

void COutput::WriteSurface_CSV(CConfig *config, CGeometry *geometry){
  
  unsigned short iVar;
  
  unsigned long iPoint;
  unsigned long iExtIter = config->GetExtIter();
  char cstr[200], buffer[50];
  ofstream Surf_file;
  
  int iProcessor;
  
  string filename;
  
  switch (config->GetKind_Solver()) {
  case EULER : case NAVIER_STOKES : case RANS :
    filename = config->GetSurfFlowCoeff_FileName();
    break;
  case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS :
  case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    filename = config->GetSurfAdjCoeff_FileName();
    break;
  default: break;
  }
  
  strcpy (cstr, filename.c_str());
  
  if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(config->GetiInst()));
    
  } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.csv",  SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.csv",   SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.csv",    SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));
  }
  else
    SPRINTF (buffer, ".csv");
  
  strcat(cstr, buffer);
  Surf_file.precision(15);
  
  
  if (rank == MASTER_NODE) {
    Surf_file.open(cstr, ios::out);
    Surf_file << "\"Point\",";
    
    for (iVar = 0; iVar < Variable_Names.size()-1; iVar++) {
      Surf_file << "\"" << Variable_Names[iVar] << "\",";
    }
    Surf_file << "\"" << Variable_Names[Variable_Names.size()-1] << "\"" << endl;
    
  }
  
  /*--- Write surface and volumetric solution data. ---*/
  
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint < nSurf_Poin_Par; iPoint++) {
        
        Surf_file << Renumber2Global[iPoint+1] << ", ";
        
        for (iVar = 0; iVar < GlobalField_Counter; iVar++){
          Surf_file << scientific << Parallel_Surf_Data[iVar][iPoint];
          if (iVar != GlobalField_Counter -1) Surf_file << ", ";
        }
        Surf_file << endl;
      }
    }
    Surf_file.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }
}