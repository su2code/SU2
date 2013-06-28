/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.2
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2013 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/output_structure.hpp"

void COutput::CleanUp(CConfig *config, CGeometry *geometry) {
  
	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Local variables and initilization ---*/
  
	unsigned short iVar, iDim;
	unsigned short nDim = geometry->GetnDim();
  
	/*--- The master node alone owns all data found in this routine. ---*/
	if (rank == MASTER_NODE) {
    
    /*--- Release data associated with geometry merging.  We'll take these
     boolean guards off soon. ---*/
    if (config->GetWrt_Sol_Tec_ASCII()  ||
        config->GetWrt_Sol_Tec_Binary() ||
        config->GetWrt_Sol_CGNS()) {
      /*--- Deallocate memory for coordinate data ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        delete [] Coords[iDim];
      }
      delete [] Coords;
      
      /*--- Deallocate memory for connectivity data ---*/
      if (nGlobal_Tria > 0) delete [] Conn_Tria;
      if (nGlobal_Quad > 0) delete [] Conn_Quad;
      if (nGlobal_Tetr > 0) delete [] Conn_Tetr;
      if (nGlobal_Hexa > 0) delete [] Conn_Hexa;
      if (nGlobal_Wedg > 0) delete [] Conn_Wedg;
      if (nGlobal_Pyra > 0) delete [] Conn_Pyra;
      if (nGlobal_Line > 0) delete [] Conn_Line;
      
    }
    
    /*--- Same with solution data. We'll take the boolean guards off here soon
     too. ---*/
    if (config->GetWrt_Sol_Tec_ASCII()  ||
        config->GetWrt_Sol_Tec_Binary() ||
        config->GetWrt_Sol_CGNS()) {
      /*--- Deallocate memory for solution data ---*/
      for (iVar = 0; iVar < nVar_Total; iVar++) {
        delete [] Data[iVar];
      }
      delete [] Data;

    }
    
	}
  
}