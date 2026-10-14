/*!
 * \file CParaviewPVTUFileWriter.cpp
 * \brief Filewriter class for the Paraview parallel XML format (.pvtu and one .vtu piece per rank).
 * \author A. Rausa
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

#include "../../../include/output/filewriter/CParaviewPVTUFileWriter.hpp"

#include <algorithm>
#include <limits>
#include <numeric>

#if defined(_WIN32)
#include <direct.h>
#else
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace {

constexpr unsigned short NCOORDS = 3;

/*--- The element types in the order in which they are written. ---*/
const std::array<GEO_TYPE, 7> ElemTypes = {LINE,       TRIANGLE, QUADRILATERAL, TETRAHEDRON,
                                           HEXAHEDRON, PRISM,    PYRAMID};

/*--- Write a data array in the "appended raw" format: the size in bytes of the array followed by its
 values, converted to the type the file declares. Returns the number of bytes written. ---*/
template <class T, class U>
size_t WriteDataArray(FILE* file, const vector<U>& values) {
  const auto nBytes = static_cast<uint64_t>(values.size() * sizeof(T));
  fwrite(&nBytes, sizeof(nBytes), 1, file);

  vector<T> converted(values.begin(), values.end());
  fwrite(converted.data(), sizeof(T), converted.size(), file);

  return sizeof(nBytes) + nBytes;
}

/*--- Open a new file for writing. On POSIX systems it is created with the permissions of the other output files
 (owner read/write, group and others read), fopen would create it with 0666 before the umask. ---*/
FILE* OpenFile(const string& name, const char* mode) {
#if defined(_WIN32)
  return fopen(name.c_str(), mode);
#else
  const int fd = open(name.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fd < 0) return nullptr;
  FILE* file = fdopen(fd, mode);
  if (!file) close(fd);
  return file;
#endif
}

/*--- Create a folder, if it does not exist yet. ---*/
void MakeFolder(const string& name) {
#if defined(_WIN32)
#ifdef __MINGW32__
  mkdir(name.c_str());
#else
  _mkdir(name.c_str());
#endif
#else
  mkdir(name.c_str(), 0777);
#endif
}

/*--- Size in bytes of a data array in the file, including the size written in front of it. ---*/
size_t DataArraySize(size_t nValues, size_t typeSize) { return sizeof(uint64_t) + nValues * typeSize; }

/*--- Name of a field without the component suffix of a vector, e.g. "Momentum_x" -> "Momentum". ---*/
string FieldBaseName(string name) {
  name.erase(remove(name.begin(), name.end(), '"'), name.end());
  return name.substr(0, name.size() - 2);
}

}  // namespace

const string CParaviewPVTUFileWriter::fileExt = ".pvtu";

CParaviewPVTUFileWriter::CParaviewPVTUFileWriter(CParallelDataSorter* valDataSorter, bool valDoublePrecision)
    : CFileWriter(valDataSorter, fileExt), doublePrecision(valDoublePrecision) {
  /* Check for big endian. We have to swap bytes otherwise. */

  bigEndian = false;
  unsigned int i = 1;
  char* c = (char*)&i;
  bigEndian = *c == 0;
}

CParaviewPVTUFileWriter::~CParaviewPVTUFileWriter() = default;

void CParaviewPVTUFileWriter::WriteData(string val_filename) {
  if (!dataSorter->GetConnectivitySorted()) {
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  if (bigEndian) {
    SU2_MPI::Error("The Paraview parallel writer is only implemented for little endian machines.", CURRENT_FUNCTION);
  }

  startTime = SU2_MPI::Wtime();

  /*--- The pieces are written into a folder next to the .pvtu file, one file per rank. ---*/

  const string folderName = val_filename;

  if (rank == MASTER_NODE) MakeFolder(folderName);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

  const auto lastSlash = val_filename.find_last_of("/\\");
  const string baseName = (lastSlash == string::npos) ? val_filename : val_filename.substr(lastSlash + 1);

  vector<string> pieceNames(size);
  for (int i = 0; i < size; i++) pieceNames[i] = baseName + "/" + baseName + "_" + to_string(i) + ".vtu";

  /*--- Make the data of this rank self-contained and write its piece. ---*/

  const auto nPointsPiece = CollectExternalPoints();

  /*--- All pieces must declare the same types, so the connectivity type is that of the largest piece. ---*/

  unsigned long localConn = dataSorter->GetnConn(), maxConn = 0;
  SU2_MPI::Allreduce(&localConn, &maxConn, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());
  const bool connInt64 = maxConn > static_cast<unsigned long>(std::numeric_limits<int32_t>::max());

  WritePiece(folderName + "/" + baseName + "_" + to_string(rank) + ".vtu", nPointsPiece, connInt64);

  /*--- The master node writes the file that lists the pieces. ---*/

  WriteIndexFile(val_filename + fileExt, pieceNames, connInt64);

  stopTime = SU2_MPI::Wtime();
  usedTime = stopTime - startTime;

  su2double myFileSize = fileSize;
  SU2_MPI::Allreduce(&myFileSize, &fileSize, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  bandwidth = fileSize / (1.0e6) / usedTime;
}

unsigned long CParaviewPVTUFileWriter::CollectExternalPoints() {
  const auto nFields = dataSorter->GetFieldNames().size();
  const auto firstPoint = dataSorter->GetnPointCumulative(rank);

  /*--- The piece holds exactly the points used by the elements of this rank, wherever they live. ---*/

  usedPoints.clear();
  usedPoints.reserve(dataSorter->GetnConn());

  for (auto type : ElemTypes) {
    const auto nNodes = nPointsOfElementType(type);
    for (unsigned long iElem = 0; iElem < dataSorter->GetnElem(type); iElem++)
      for (unsigned short iNode = 0; iNode < nNodes; iNode++)
        usedPoints.push_back(dataSorter->GetElemConnectivity(type, iElem, iNode) - 1);
  }

  sort(usedPoints.begin(), usedPoints.end());
  usedPoints.erase(unique(usedPoints.begin(), usedPoints.end()), usedPoints.end());
  usedPoints.shrink_to_fit();

  if (usedPoints.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
    SU2_MPI::Error("Too many points in one piece of the Paraview parallel output.", CURRENT_FUNCTION);

  /*--- Ask the rank that owns each point for its data. The points are sorted by their global index, and
   the ranks own ranges of increasing global indices, so they are already grouped by rank. ---*/

  vector<int> nRequest(size, 0), nAnswer(size, 0), requestOffset(size + 1, 0), answerOffset(size + 1, 0);

  for (auto iPoint : usedPoints) nRequest[dataSorter->FindProcessor(iPoint)]++;

  SU2_MPI::Alltoall(nRequest.data(), 1, MPI_INT, nAnswer.data(), 1, MPI_INT, SU2_MPI::GetComm());

  for (int i = 0; i < size; i++) {
    requestOffset[i + 1] = requestOffset[i] + nRequest[i];
    answerOffset[i + 1] = answerOffset[i] + nAnswer[i];
  }

  vector<unsigned long> answeredPoints(answerOffset[size]);

  SU2_MPI::Alltoallv(usedPoints.data(), nRequest.data(), requestOffset.data(), MPI_UNSIGNED_LONG,
                     answeredPoints.data(), nAnswer.data(), answerOffset.data(), MPI_UNSIGNED_LONG,
                     SU2_MPI::GetComm());

  /*--- Answer with the data of the points this rank owns, one field at a time so that the MPI counts are numbers
   of points. The data is passive, so the passive MPI wrapper is used also in AD builds. ---*/

  const auto nPointsPiece = usedPoints.size();
  pieceData.resize(nPointsPiece * nFields);

  vector<passivedouble> answerData(answeredPoints.size());

  for (size_t iField = 0; iField < nFields; iField++) {
    for (size_t i = 0; i < answeredPoints.size(); i++)
      answerData[i] = dataSorter->GetData(iField, answeredPoints[i] - firstPoint);

    SelectMPIWrapper<passivedouble>::W::Alltoallv(answerData.data(), nAnswer.data(), answerOffset.data(), MPI_DOUBLE,
                                                  &pieceData[iField * nPointsPiece], nRequest.data(),
                                                  requestOffset.data(), MPI_DOUBLE, SU2_MPI::GetComm());
  }

  return usedPoints.size();
}

void CParaviewPVTUFileWriter::WritePiece(const string& val_filename, unsigned long nPointsPiece, bool connInt64) {
  const auto& fieldNames = dataSorter->GetFieldNames();
  const auto nFields = fieldNames.size();
  const auto nDim = dataSorter->GetnDim();

  /*--- Local index of a point of the piece. ---*/

  auto localIndex = [&](unsigned long iPoint) {
    return static_cast<unsigned long>(lower_bound(usedPoints.begin(), usedPoints.end(), iPoint) - usedPoints.begin());
  };

  /*--- Connectivity and cell types of the piece. ---*/

  const auto nElem = dataSorter->GetnElem();
  const auto nConn = dataSorter->GetnConn();

  vector<unsigned long> connectivity, offsets;
  vector<uint8_t> types;
  connectivity.reserve(nConn);
  offsets.reserve(nElem);
  types.reserve(nElem);

  for (auto type : ElemTypes) {
    const auto nNodes = nPointsOfElementType(type);
    for (unsigned long iElem = 0; iElem < dataSorter->GetnElem(type); iElem++) {
      for (unsigned short iNode = 0; iNode < nNodes; iNode++)
        connectivity.push_back(localIndex(dataSorter->GetElemConnectivity(type, iElem, iNode) - 1));
      offsets.push_back(connectivity.size());
      types.push_back(static_cast<uint8_t>(type));
    }
  }

  /*--- Sizes of the arrays in the file, needed for the offsets of the appended data. ---*/

  const size_t realSize = doublePrecision ? sizeof(double) : sizeof(float);
  const string realType = doublePrecision ? "Float64" : "Float32";

  const size_t connSize = connInt64 ? sizeof(int64_t) : sizeof(int32_t);
  const string connType = connInt64 ? "Int64" : "Int32";

  /*--- Header of the piece, the offset of each array in the appended data is written with it. ---*/

  ostringstream header;
  size_t offset = 0;

  header << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" "
            "header_type=\"UInt64\">\n";
  header << "<UnstructuredGrid>\n";
  header << "<Piece NumberOfPoints=\"" << nPointsPiece << "\" NumberOfCells=\"" << nElem << "\">\n";

  header << "<Points>\n<DataArray type=\"" << realType << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\""
         << offset << "\"/>\n</Points>\n";
  offset += DataArraySize(nPointsPiece * NCOORDS, realSize);

  header << "<Cells>\n";
  header << "<DataArray type=\"" << connType << "\" Name=\"connectivity\" format=\"appended\" offset=\"" << offset
         << "\"/>\n";
  offset += DataArraySize(connectivity.size(), connSize);
  header << "<DataArray type=\"" << connType << "\" Name=\"offsets\" format=\"appended\" offset=\"" << offset
         << "\"/>\n";
  offset += DataArraySize(offsets.size(), connSize);
  header << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << offset << "\"/>\n";
  offset += DataArraySize(types.size(), sizeof(uint8_t));
  header << "</Cells>\n";

  /*--- The fields start after the coordinates, a field whose name ends in _x is written as a vector
   together with the fields that follow it. ---*/

  const size_t varStart = (nDim == 3) ? 3 : 2;

  header << "<PointData>\n";
  for (size_t iField = varStart; iField < nFields; iField++) {
    if (fieldNames[iField].find("_y") != string::npos || fieldNames[iField].find("_z") != string::npos) continue;

    const bool isVector = fieldNames[iField].find("_x") != string::npos;
    string name = fieldNames[iField];
    if (isVector) {
      name = FieldBaseName(name);
    } else {
      name.erase(remove(name.begin(), name.end(), '"'), name.end());
    }

    header << "<DataArray type=\"" << realType << "\" Name=\"" << name << "\" NumberOfComponents=\""
           << (isVector ? NCOORDS : 1) << "\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += DataArraySize(nPointsPiece * (isVector ? NCOORDS : 1), realSize);
  }
  header << "</PointData>\n</Piece>\n</UnstructuredGrid>\n<AppendedData encoding=\"raw\">\n_";

  /*--- Each rank writes its own file, no rank waits for another one. The arrays are written one at a time. ---*/

  FILE* file = OpenFile(val_filename, "wb");
  if (!file) SU2_MPI::Error(string("Unable to open file ") + val_filename, CURRENT_FUNCTION);

  const auto headerString = header.str();
  fwrite(headerString.data(), sizeof(char), headerString.size(), file);
  fileSize = headerString.size();

  vector<passivedouble> buffer;

  auto writeField = [&](size_t iField, unsigned short nComponents) {
    buffer.assign(nPointsPiece * nComponents, 0.0);
    for (unsigned short iComp = 0; iComp < nComponents; iComp++) {
      if ((nDim == 2) && (nComponents == NCOORDS) && (iComp == 2)) continue;
      const auto* field = &pieceData[(iField + iComp) * nPointsPiece];
      for (unsigned long iPoint = 0; iPoint < nPointsPiece; iPoint++)
        buffer[iPoint * nComponents + iComp] = field[iPoint];
    }
    fileSize += doublePrecision ? WriteDataArray<double>(file, buffer) : WriteDataArray<float>(file, buffer);
  };

  auto writeConn = [&](const vector<unsigned long>& values) {
    fileSize += connInt64 ? WriteDataArray<int64_t>(file, values) : WriteDataArray<int32_t>(file, values);
  };

  writeField(0, NCOORDS);

  writeConn(connectivity);
  writeConn(offsets);
  fileSize += WriteDataArray<uint8_t>(file, types);

  for (size_t iField = varStart; iField < nFields; iField++) {
    if (fieldNames[iField].find("_y") != string::npos || fieldNames[iField].find("_z") != string::npos) continue;

    const bool isVector = fieldNames[iField].find("_x") != string::npos;
    writeField(iField, isVector ? NCOORDS : 1);
  }

  const string tail = "\n</AppendedData>\n</VTKFile>\n";
  fwrite(tail.data(), sizeof(char), tail.size(), file);
  fclose(file);

  fileSize += tail.size();
}

void CParaviewPVTUFileWriter::WriteIndexFile(const string& val_filename, const vector<string>& pieceNames,
                                             bool connInt64) {
  if (rank != MASTER_NODE) return;

  const auto& fieldNames = dataSorter->GetFieldNames();
  const auto nDim = dataSorter->GetnDim();

  const string realType = doublePrecision ? "Float64" : "Float32";
  const string connType = connInt64 ? "Int64" : "Int32";

  ostringstream file;

  file << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
  file << "<PUnstructuredGrid GhostLevel=\"0\">\n";
  file << "<PPoints>\n<PDataArray type=\"" << realType << "\" NumberOfComponents=\"3\"/>\n</PPoints>\n";
  file << "<PCells>\n";
  file << "<PDataArray type=\"" << connType << "\" Name=\"connectivity\"/>\n";
  file << "<PDataArray type=\"" << connType << "\" Name=\"offsets\"/>\n";
  file << "<PDataArray type=\"UInt8\" Name=\"types\"/>\n";
  file << "</PCells>\n";

  const size_t varStart = (nDim == 3) ? 3 : 2;

  file << "<PPointData>\n";
  for (size_t iField = varStart; iField < fieldNames.size(); iField++) {
    if (fieldNames[iField].find("_y") != string::npos || fieldNames[iField].find("_z") != string::npos) continue;

    const bool isVector = fieldNames[iField].find("_x") != string::npos;
    string name = fieldNames[iField];
    if (isVector) {
      name = FieldBaseName(name);
    } else {
      name.erase(remove(name.begin(), name.end(), '"'), name.end());
    }

    file << "<PDataArray type=\"" << realType << "\" Name=\"" << name << "\" NumberOfComponents=\""
         << (isVector ? NCOORDS : 1) << "\"/>\n";
  }
  file << "</PPointData>\n";

  for (const auto& piece : pieceNames) file << "<Piece Source=\"" << piece << "\"/>\n";

  file << "</PUnstructuredGrid>\n</VTKFile>\n";

  const auto fileString = file.str();

  FILE* indexFile = OpenFile(val_filename, "w");
  if (!indexFile) SU2_MPI::Error(string("Unable to open file ") + val_filename, CURRENT_FUNCTION);
  fwrite(fileString.data(), sizeof(char), fileString.size(), indexFile);
  fclose(indexFile);

  fileSize += fileString.size();
}
