/*!
 * \file CSU2MeshBinaryFileWriter.cpp
 * \brief Filewriter class SU2 binary mesh format.
 * \author E. van der Weide
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../../include/output/filewriter/CSU2MeshBinaryFileWriter.hpp"
#include "../../../../Common/include/toolboxes/printing_toolbox.hpp"

const string CSU2MeshBinaryFileWriter::fileExt = ".su2b";

CSU2MeshBinaryFileWriter::CSU2MeshBinaryFileWriter(CParallelDataSorter *valDataSorter,
                                                   unsigned short valiZone, unsigned short valnZone) :
   CFileWriter(valDataSorter, fileExt), iZone(valiZone), nZone(valnZone) {}

void CSU2MeshBinaryFileWriter::WriteData(string val_filename) {
  SU2_MPI::Error("Function not implemented yet.", CURRENT_FUNCTION);
}
