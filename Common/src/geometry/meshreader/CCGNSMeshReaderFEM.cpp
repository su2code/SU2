/*!
 * \file CCGNSMeshReaderFEM.cpp
 * \brief Class that reads a single zone of a CGNS mesh file from disk into
 *        linear partitions across all ranks.
 * \author T. Economon, E. van der Weide
 * \version 7.0.5 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../../include/geometry/meshreader/CCGNSMeshReaderFEM.hpp"

CCGNSMeshReaderFEM::CCGNSMeshReaderFEM(CConfig        *val_config,
                                       unsigned short val_iZone,
                                       unsigned short val_nZone)
: CCGNSMeshReaderBase(val_config, val_iZone, val_nZone) {

#ifdef HAVE_CGNS
  
  /*--- The CGNS reader currently assumes a single database. ---*/
  cgnsBase = 1;
  OpenCGNSFile(config->GetMesh_FileName());
  
  /*--- Read the basic information about the database and zone(s). ---*/
  ReadCGNSDatabaseMetadata();
  ReadCGNSZoneMetadata();

  /*--- Read the basic information about the sections. ---*/
  ReadCGNSSectionMetadata();

  SU2_MPI::Error(string("Not fully implemented yet"), CURRENT_FUNCTION);
  
  /*--- We have extracted all CGNS data. Close the CGNS file. ---*/
  if (cg_close(cgnsFileID)) cg_error_exit();

#endif
}

CCGNSMeshReaderFEM::~CCGNSMeshReaderFEM(void) { }

