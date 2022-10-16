/*!
 * \file CFileReaderLUT.hpp
 * \brief reading lookup table for tabulated fluid properties
 * \author D. Mayer, T. Economon
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../Common/include/containers/CFileReaderLUT.hpp"
#include "../../Common/include/option_structure.hpp"
#include "../../Common/include/parallelization/mpi_structure.hpp"
#include "../../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../../Common/include/toolboxes/CSquareMatrixCM.hpp"

using namespace std;

void CFileReaderLUT::ReadRawDRG(const string& file_name) {
  version_reader = "1.0.1";

  /*--- Store MPI rank. ---*/
  rank = SU2_MPI::GetRank();

  string line;
  string word;

  istringstream stream_names_var;

  ifstream file_stream;

  bool eoHeader = false;
  bool eoData = false;
  bool eoConnectivity = false;
  bool eoHull = false;

  file_stream.open(file_name.c_str(), ifstream::in);

  if (!file_stream.is_open()) {
    SU2_MPI::Error(string("There is no look-up-table file called ") + file_name, CURRENT_FUNCTION);
  }

  /* Read header */
  line = SkipToFlag(&file_stream, "<Header>");

  while (getline(file_stream, line) && !eoHeader) {

    /*--- first strip any possible carriage returns ---*/
    //if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r' )) {
    //  line.erase(line.length()-1);
    //}
    line = GetStrippedLine(line);

    /* number of points in LUT */
    if (line.compare("[Version]") == 0) {
      getline(file_stream, line);

      //if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r' )) {
      //  line.erase(line.length()-1);
      //}
      line = GetStrippedLine(line);
      version_lut = line;
    }

    /* number of points in LUT */
    if (line.compare("[Number of points]") == 0) {
      getline(file_stream, line);
      line = GetStrippedLine(line);
      n_points = stoi(line);
    }

    /* number of triangles in LUT */
    if (line.compare("[Number of triangles]") == 0) {
      getline(file_stream, line);
      line = GetStrippedLine(line);
      n_triangles = stoi(line);
    }

    /* number of points on the hull */
    if (line.compare("[Number of hull points]") == 0) {
      getline(file_stream, line);
      line = GetStrippedLine(line);
      n_hull_points = stoi(line);
    }

    /* number of variables in LUT */
    if (line.compare("[Number of variables]") == 0) {
      getline(file_stream, line);
      line = GetStrippedLine(line);
      n_variables = stoi(line);
    }

    /* variable names */
    if (line.compare("[Variable names]") == 0) {
      
      for (unsigned long i = 0; i < n_variables; i++){

        /*--- grab a single line ---*/
        getline(file_stream, line);
        line = GetStrippedLine(line);
        names_var.push_back(line.substr(line.find(":")+1)); 
      }
    }

    // check if end of header is reached
    if (line.compare("</Header>") == 0) eoHeader = true;
  }

  // check version_lut
  if (version_lut.compare(version_reader) != 0)
    SU2_MPI::Error("Version conflict between Dragon reader and Dragon library file.", CURRENT_FUNCTION);

  // check header quantities
  if (n_points == 0 || n_triangles == 0 || n_variables == 0 || n_hull_points == 0)
    SU2_MPI::Error(
        "Number of points, triangles, hull points, or variables in Dragon "
        "library header is zero.",
        CURRENT_FUNCTION);

  // check if number of variables is consistent
  if (n_variables != names_var.size())
    SU2_MPI::Error(
        "Number of read variables does not match number of "
        "variables specified in Dragon "
        "library header.",
        CURRENT_FUNCTION);

  /* now that n_variables, n_points, n_hull_points and n_variables is available,
   * allocate memory */
  if (rank == MASTER_NODE) cout << "allocating memory for the data" << endl;
  table_data.resize(GetNVariables(), GetNPoints());


  if (rank == MASTER_NODE) cout << "allocating memory for the triangles" << endl;
  triangles.resize(GetNTriangles(), 3);

  if (rank == MASTER_NODE) cout << "allocating memory for the hull points" << endl;
  hull.resize(GetNHullPoints());

  /* flush any cout */
  if (rank == MASTER_NODE) cout << endl;

  // read data block
  if (rank == MASTER_NODE) cout << "loading data block" << endl;

  line = SkipToFlag(&file_stream, "<Data>");

  unsigned long pointCounter = 0;
    while (getline(file_stream, line) && !eoData) {
    if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r' )) {
      line.erase(line.length()-1);
    }
    // check if end of data is reached
    if (line.compare("</Data>") == 0) eoData = true;

    if (!eoData) {
      // one line contains values for one point for all variables
      istringstream streamDataLine(line);

      // add next line to table array
      for (unsigned long iVar = 0; iVar < n_variables; iVar++) {
        streamDataLine >> word;
        passivedouble tmp = stod(word);
        table_data[iVar][pointCounter] = (su2double)tmp;
      }
    }
    pointCounter++;
  }

  if (n_points != pointCounter - 1)
    SU2_MPI::Error(
        "Number of read points does not match number of points "
        "specified in Dragon library header.",
        CURRENT_FUNCTION);

  // read connectivity
  if (rank == MASTER_NODE) cout << "loading connectivity block" << endl;

  line = SkipToFlag(&file_stream, "<Connectivity>");

  unsigned long triCounter = 0;
  while (getline(file_stream, line) && !eoConnectivity) {
    if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r' )) {
      line.erase(line.length()-1);
    }
    // check if end of data is reached
    if (line.compare("</Connectivity>") == 0) eoConnectivity = true;

    if (!eoConnectivity) {
      // one line contains values for one triangle (3 points)
      istringstream streamTriLine(line);

      // add next line to triangles
      for (int iPoint = 0; iPoint < 3; iPoint++) {
        streamTriLine >> word;
        // Dragon table index starts with 1, convert to c++ indexing starting
        // with 0:
        triangles[triCounter][iPoint] = stol(word) - 1;
      }
    }
    triCounter++;
  }

  if (n_triangles != triCounter - 1)
    SU2_MPI::Error(
        "Number of read triangles does not match number of points "
        "specified in Dragon library header.",
        CURRENT_FUNCTION);

  // read hull points
  if (rank == MASTER_NODE) cout << "loading hull block" << endl;

  line = SkipToFlag(&file_stream, "<Hull>");

  unsigned long hullCounter = 0;
  while (getline(file_stream, line) && !eoHull) {
    if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r' )) {
      line.erase(line.length()-1);
    }
    // check if end of data is reached
    if (line.compare("</Hull>") == 0) eoHull = true;

    if (!eoHull) {
      // one line contains one point ID for one point on the hull
      istringstream streamHullLine(line);

      streamHullLine >> word;

      // Dragon table indices start with 1, convert to c++ indexing starting
      // with 0:
      hull[hullCounter] = stol(word) - 1;
    }
    hullCounter++;
  }

  if (n_hull_points != hullCounter - 1)
    SU2_MPI::Error(
        "Number of read hull points does not match number of points "
        "specified in Dragon library header.",
        CURRENT_FUNCTION);

  file_stream.close();

}

string CFileReaderLUT::SkipToFlag(ifstream* file_stream, const string& flag) {
  string line;
  getline(*file_stream, line);

  while (line.find(flag) == std::string::npos && !(*file_stream).eof()) {
    getline(*file_stream, line);
  }

  if ((*file_stream).eof()) SU2_MPI::Error("Flag " + flag + " not found in file", CURRENT_FUNCTION);

  return line;
}
