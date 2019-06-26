#include "../../../include/output/filewriter/CParallelDataSorter.hpp"
#include "../../../Common/include/fem_geometry_structure.hpp"

CParallelDataSorter::CParallelDataSorter(CConfig *config, unsigned short nFields){
  
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  GlobalField_Counter = nFields;
  
  nParallel_Hexa = 0;
  nParallel_Line = 0;
  nParallel_Quad = 0;
  nParallel_Tetr = 0;
  nParallel_Pris = 0;
  nParallel_Pyra = 0;
  nParallel_Tria = 0;
  
  Conn_Line_Par = NULL;
  Conn_Hexa_Par = NULL;
  Conn_Pris_Par = NULL;
  Conn_Quad_Par = NULL;
  Conn_Tetr_Par = NULL;
  Conn_Tria_Par = NULL;
  Conn_Pyra_Par = NULL;
  
  Parallel_Data = NULL;
}

CParallelDataSorter::~CParallelDataSorter(){
  
  DeallocateConnectivity();
  
  DeallocateData();
  
}


unsigned long CParallelDataSorter::GetnElem(GEO_TYPE type){
  
  switch (type) {
    case LINE:
      return nParallel_Line;
      break;
    case TRIANGLE:
      return nParallel_Tria;
      break;
    case QUADRILATERAL:
      return nParallel_Quad;
      break;
    case TETRAHEDRON:
      return nParallel_Tetr;
      break;
    case HEXAHEDRON:
      return nParallel_Hexa;
      break;
    case PRISM:
      return nParallel_Pris;
      break;
    case PYRAMID:
      return nParallel_Pyra;
      break;
    default:
      break;
  }
  
  SU2_MPI::Error("GEO_TYPE not found", CURRENT_FUNCTION);
  
  return 0;
}


unsigned long CParallelDataSorter::GetElem_Connectivity(GEO_TYPE type, unsigned long iElem, unsigned long iNode) {
  
  switch (type) {
    case LINE:
      return Conn_Line_Par[iElem*N_POINTS_LINE + iNode];
      break;
    case TRIANGLE:
      return Conn_Tria_Par[iElem*N_POINTS_TRIANGLE + iNode];
      break;
    case QUADRILATERAL:
      return Conn_Quad_Par[iElem*N_POINTS_QUADRILATERAL + iNode];
      break;
    case TETRAHEDRON:
      return Conn_Tetr_Par[iElem*N_POINTS_TETRAHEDRON + iNode];
      break;
    case HEXAHEDRON:
      return Conn_Hexa_Par[iElem*N_POINTS_HEXAHEDRON + iNode];
      break;
    case PRISM:
      return Conn_Pris_Par[iElem*N_POINTS_PRISM + iNode];
      break;
    case PYRAMID:
      return Conn_Pyra_Par[iElem*N_POINTS_PYRAMID + iNode];
      break;
    default:
      break;
  }
  
  SU2_MPI::Error("GEO_TYPE not found", CURRENT_FUNCTION);
  
  return 0;
}

void CParallelDataSorter::DeallocateConnectivity() {
  
  /*--- Deallocate memory for connectivity data on each processor. ---*/
  
  if (nParallel_Line > 0 && Conn_Line_Par != NULL) delete [] Conn_Line_Par;
  if (nParallel_Tria > 0 && Conn_Tria_Par != NULL) delete [] Conn_Tria_Par;
  if (nParallel_Quad > 0 && Conn_Quad_Par != NULL) delete [] Conn_Quad_Par;
  if (nParallel_Tetr > 0 && Conn_Tetr_Par != NULL) delete [] Conn_Tetr_Par;
  if (nParallel_Hexa > 0 && Conn_Hexa_Par != NULL) delete [] Conn_Hexa_Par;
  if (nParallel_Pris > 0 && Conn_Pris_Par != NULL) delete [] Conn_Pris_Par;
  if (nParallel_Pyra > 0 && Conn_Pyra_Par != NULL) delete [] Conn_Pyra_Par;
  
  
}

void CParallelDataSorter::DeallocateData() {
  
  /*--- Deallocate memory for solution data ---*/
  
  for (unsigned short iVar = 0; iVar < GlobalField_Counter; iVar++) {
    if (Parallel_Data[iVar] != NULL) delete [] Parallel_Data[iVar];
  }
  if (Parallel_Data != NULL) delete [] Parallel_Data;

}
