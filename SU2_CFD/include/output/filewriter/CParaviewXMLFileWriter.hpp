/*!
 * \file CParaviewXMLFileWriter.hpp
 * \brief Headers fo paraview binary file writer class.
 * \author T. Albring
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#pragma once

#include "CFileWriter.hpp"

class CParaviewXMLFileWriter final: public CFileWriter{

private:

  /*!
   * \brief Enum that defines the VTK datatypes
   */
  enum class VTKDatatype {
    FLOAT32,
    INT32,
    UINT8
  };

  /*!
   * \brief Boolean storing whether we are on a big or little endian machine
   */
  bool bigEndian;

  /*!
   * \brief The current data offset that is used to find data in the binary blob at the end of the file
   */
  unsigned long dataOffset;

public:

  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] val_filename - The name of the file
   * \param[in] valDataSorter - The parallel sorted data to write
   */
  CParaviewXMLFileWriter(string val_filename, CParallelDataSorter* valDataSorter);

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] valDataSorter - The parallel sorted data to write
   */
  CParaviewXMLFileWriter(CParallelDataSorter* valDataSorter);

  /*!
   * \brief Destructor
   */
  ~CParaviewXMLFileWriter() override;

  /*!
   * \brief Write sorted data to file in paraview binary file format
   */
  void WriteData(string val_filename) override ;

private:

  /*!
   * \brief Add a new data array definition to the vtu file.
   * \param[in] type - The vtk datatype
   * \param[in] name - The name of the array
   * \param[in] nComponents - The number of components
   * \param[in] size        - The total size of the array
   * \param[in] globalSize  - The global size of the array over all processors
   */
  void AddDataArray(VTKDatatype type, string name, unsigned short nComponents, unsigned long size, unsigned long globalSize);

  /*!
   * \brief Write an array that has previously been defined with ::AddDataArray to the vtu file in binary format
   * \param[in] data - Pointer to the data
   * \param[in] type - The vtk datatype
   * \param[in] size - The total size of the array
   * \param[in] globalSize - The global size of the array over all processors
   * \param[in] offset - The displacement in the file view for the current processor
   */
  void WriteDataArray(void *data, VTKDatatype type, unsigned long size, unsigned long globalSize, unsigned long offset);

  /*!
   * \brief Get the type string and size of a VTK datatype
   * \param[in]  type - The VTK datatype
   * \param[out] typeStr - The string name of the type
   * \param[out] typeSize - The size in bytes of the type
   */
  inline void GetTypeInfo(const VTKDatatype type, string &typeStr, unsigned long &typeSize) const {
    switch (type) {
      case VTKDatatype::FLOAT32:
        typeStr = "\"Float32\"";
        typeSize = sizeof(float);
        break;
      case VTKDatatype::INT32:
        typeStr = "\"Int32\"";
        typeSize = sizeof(int);
        break;
      case VTKDatatype::UINT8:
        typeStr = "\"UInt8\"";
        typeSize = sizeof(char);
        break;
      default:
        SU2_MPI::Error("Unknown Type", CURRENT_FUNCTION);
        break;
    }
  }
};

