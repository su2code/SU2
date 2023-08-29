/*!
 * \file CFileReaderLUT.cpp
 * \brief reading lookup table for tabulated fluid properties
 * \author D. Mayer, T. Economon
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

#include "../../Common/include/containers/CFileReaderLUT.hpp"

#include "../../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../../Common/include/toolboxes/CSquareMatrixCM.hpp"
#include "../../Common/include/option_structure.hpp"
#include "../../Common/include/parallelization/mpi_structure.hpp"

using namespace std;

void CFileReaderLUT::ReadRawLUT(const string& file_name) {
  version_reader = "1.1.0";

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

  /*--- Read header ---*/
  SkipToFlag(file_stream, line, "<Header>");
  table_dim = 2;
  bool found_level_count = false;
  while (GetNextNonEmptyLine(file_stream, line) && !eoHeader) {
    /*--- check version_lut ---*/
    if (line.compare("[Version]") == 0) {
      GetNextNonEmptyLine(file_stream, line);
      version_lut = line;

      if (version_lut.compare("1.0.1") == 0) {
        table_dim = 2;
        n_levels = 1;
        found_level_count = true;
        n_points.resize(n_levels);
        n_triangles.resize(n_levels);
        n_hull_points.resize(n_levels);
      } else if (version_lut.compare("1.1.0") == 0) {
        table_dim = 3;
      } else {
        SU2_MPI::Error("Version conflict between LUT reader and LUT library file.", CURRENT_FUNCTION);
      }
    }

    if (line.compare("[Number of table levels]") == 0) {
      found_level_count = true;
      GetNextNonEmptyLine(file_stream, line);
      n_levels = stoul(line);

      n_points.resize(n_levels);
      n_triangles.resize(n_levels);
      n_hull_points.resize(n_levels);
      table_levels.resize(n_levels);
    }

    /*--- number of points in LUT ---*/
    if (line.compare("[Number of points]") == 0) {
      if (!found_level_count)
        SU2_MPI::Error("Number of points provided before specifying level count.", CURRENT_FUNCTION);
      for (unsigned long i_level = 0; i_level < n_levels; i_level++) {
        GetNextNonEmptyLine(file_stream, line);
        n_points[i_level] = stoi(line);
      }
    }

    /*--- number of triangles in LUT ---*/
    if (line.compare("[Number of triangles]") == 0) {
      if (!found_level_count)
        SU2_MPI::Error("Number of triangles provided before specifying level count.", CURRENT_FUNCTION);
      for (unsigned long i_level = 0; i_level < n_levels; i_level++) {
        GetNextNonEmptyLine(file_stream, line);
        n_triangles[i_level] = stoi(line);
      }
    }

    /*--- number of points on the hull ---*/
    if (line.compare("[Number of hull points]") == 0) {
      if (!found_level_count)
        SU2_MPI::Error("Number of hull points provided before specifying level count.", CURRENT_FUNCTION);
      for (unsigned long i_level = 0; i_level < n_levels; i_level++) {
        GetNextNonEmptyLine(file_stream, line);
        n_hull_points[i_level] = stoi(line);
      }
    }

    if (line.compare("[Table levels]") == 0) {
      for (unsigned long i_level = 0; i_level < n_levels; i_level++) {
        GetNextNonEmptyLine(file_stream, line);
        table_levels[i_level] = stod(line);
      }

      /*--- Checking whether table levels are ordered and appearance of duplicates ---*/
      for (auto i_level = 0u; i_level < n_levels - 1; i_level++) {
        if (table_levels[i_level] > table_levels[i_level + 1])
          SU2_MPI::Error("Table levels are not provided in ascending order.", CURRENT_FUNCTION);
      }
      auto duplicate = adjacent_find(table_levels.begin(), table_levels.end()) != table_levels.end();
      if (duplicate) SU2_MPI::Error("Duplicate table levels are present in LUT file", CURRENT_FUNCTION);
    }

    /*--- number of variables in LUT ---*/
    if (line.compare("[Number of variables]") == 0) {
      GetNextNonEmptyLine(file_stream, line);
      n_variables = stoi(line);
    }

    /*--- variable names ---*/
    if (line.compare("[Variable names]") == 0) {
      names_var.resize(n_variables);
      for (unsigned long i = 0; i < n_variables; i++) {
        /*--- grab a single line ---*/
        GetNextNonEmptyLine(file_stream, line);
        names_var[i] = line.substr(line.find(':') + 1);
      }
    }

    /*--- check if end of header is reached ---*/
    if (line.compare("</Header>") == 0) eoHeader = true;
  }

  /*--- check header quantities ---*/
  for (unsigned long i_level = 0; i_level < n_levels; i_level++) {
    if (n_points[i_level] == 0 || n_triangles[i_level] == 0 || n_variables == 0 || n_hull_points[i_level] == 0)
      SU2_MPI::Error(
          "Number of points, triangles, hull points, or variables in lookup table "
          "library header is zero.",
          CURRENT_FUNCTION);
  }

  /*--- check if number of variables is consistent ---*/
  if (n_variables != names_var.size())
    SU2_MPI::Error(
        "Number of read variables does not match number of "
        "variables specified in lookup table "
        "library header.",
        CURRENT_FUNCTION);

  /*--- now that n_variables, n_points, n_hull_points and n_variables are available, allocate memory ---*/
  if (rank == MASTER_NODE)
    cout << "allocating memory for the data, size = ( " << GetNVariables() << " , " << GetNPoints() << " )" << endl;
  table_data.resize(n_levels);
  for (unsigned long i_level = 0; i_level < n_levels; i_level++)
    table_data[i_level].resize(GetNVariables(), GetNPoints(i_level));

  if (rank == MASTER_NODE) cout << "allocating memory for the triangles, size = " << GetNTriangles() << endl;
  triangles.resize(n_levels);
  for (unsigned long i_level = 0; i_level < n_levels; i_level++) triangles[i_level].resize(GetNTriangles(i_level), 3);

  if (rank == MASTER_NODE) cout << "allocating memory for the hull points, size = " << GetNHullPoints() << endl;
  hull.resize(n_levels);
  for (unsigned long i_level = 0; i_level < n_levels; i_level++) hull[i_level].resize(GetNHullPoints(i_level));

  /*--- flush any cout ---*/
  if (rank == MASTER_NODE) cout << endl;

  /*--- read data block ---*/
  if (rank == MASTER_NODE) cout << "loading data block" << endl;

  SkipToFlag(file_stream, line, "<Data>");
  bool eoLevel = true;
  unsigned long pointCounter = 0;
  unsigned long levelCounter = 0;
  while (GetNextNonEmptyLine(file_stream, line) && !eoData) {
    /*--- check if end of data is reached ---*/
    if (line.compare("</Data>") == 0) eoData = true;

    if (!eoData) {
      if (line.compare("<Level>") == 0) {
        eoLevel = false;
        pointCounter = 0;
        GetNextNonEmptyLine(file_stream, line);
      }
      if (line.compare("</Level>") == 0) {
        eoLevel = true;
        if (n_points[levelCounter] != pointCounter)
          SU2_MPI::Error("Number of read points on level " + std::to_string(levelCounter) +
                             " does not match number of points "
                             "specified in lookup table library header.",
                         CURRENT_FUNCTION);
        levelCounter++;
      }
      if (!eoLevel || (table_dim == 2)) {
        /*--- one line contains values for one point for all variables ---*/
        istringstream streamDataLine(line);

        /*--- add next line to table array ---*/
        for (unsigned long iVar = 0; iVar < n_variables; iVar++) {
          streamDataLine >> word;
          passivedouble tmp = stod(word);
          table_data[levelCounter][iVar][pointCounter] = (su2double)tmp;
        }
        pointCounter++;
      }
    }
  }

  /*--- read connectivity ---*/
  if (rank == MASTER_NODE) cout << "loading connectivity block" << endl;

  SkipToFlag(file_stream, line, "<Connectivity>");
  unsigned long triCounter = 0;
  eoLevel = true;
  levelCounter = 0;
  while (GetNextNonEmptyLine(file_stream, line) && !eoConnectivity) {
    if (!line.empty() && (line[line.length() - 1] == '\n' || line[line.length() - 1] == '\r')) {
      line.erase(line.length() - 1);
    }
    /*--- check if end of data is reached ---*/
    if (line.compare("</Connectivity>") == 0) eoConnectivity = true;

    if (!eoConnectivity) {
      if (line.compare("<Level>") == 0) {
        eoLevel = false;
        triCounter = 0;
        GetNextNonEmptyLine(file_stream, line);
      }
      if (line.compare("</Level>") == 0) {
        eoLevel = true;
        if (n_triangles[levelCounter] != triCounter)
          SU2_MPI::Error("Number of read triangles on level " + std::to_string(levelCounter) +
                             " does not match number of triangles "
                             "specified in lookup table library header.",
                         CURRENT_FUNCTION);
        levelCounter++;
      }
      if (!eoLevel || (table_dim == 2)) {
        /*--- one line contains values for one triangle (3 points) ---*/
        istringstream streamTriLine(line);

        /*--- add next line to triangles ---*/
        for (int iPoint = 0; iPoint < 3; iPoint++) {
          streamTriLine >> word;

          /*--- lookup table index starts with 1, convert to c++ indexing starting with 0: ---*/
          triangles[levelCounter][triCounter][iPoint] = stol(word) - 1;
        }
        triCounter++;
      }
    }
  }

  /*--- read hull points ---*/
  if (rank == MASTER_NODE) cout << "loading hull block" << endl;

  SkipToFlag(file_stream, line, "<Hull>");
  unsigned long hullCounter = 0;
  eoLevel = true;
  levelCounter = 0;
  while (GetNextNonEmptyLine(file_stream, line) && !eoHull) {
    if (!line.empty() && (line[line.length() - 1] == '\n' || line[line.length() - 1] == '\r')) {
      line.erase(line.length() - 1);
    }
    /*--- check if end of data is reached ---*/
    if (line.compare("</Hull>") == 0) eoHull = true;

    if (!eoHull) {
      if (line.compare("<Level>") == 0) {
        eoLevel = false;
        hullCounter = 0;
        GetNextNonEmptyLine(file_stream, line);
      }
      if (line.compare("</Level>") == 0) {
        eoLevel = true;
        if (n_hull_points[levelCounter] != hullCounter)
          SU2_MPI::Error("Number of read hull nodes on level " + std::to_string(levelCounter) +
                             " does not match number of hull nodes "
                             "specified in lookup table library header.",
                         CURRENT_FUNCTION);
        levelCounter++;
      }
      if (!eoLevel || (table_dim == 2)) {
        /*--- one line contains one point ID for one point on the hull  ---*/
        istringstream streamHullLine(line);

        streamHullLine >> word;

        /*--- Lookup table indices start with 1, convert to c++ indexing starting with 0: ---*/
        hull[levelCounter][hullCounter] = stol(word) - 1;
        hullCounter++;
      }
    }
  }

  file_stream.close();
}

void CFileReaderLUT::SkipToFlag(ifstream& file_stream, const string& current_line, const string& flag) const {
  string next_line;

  /*--- compare current line to flag and return if equal ---*/
  if (current_line.compare(flag) == 0) return;

  /*--- else, search for flag ---*/
  while (next_line.find(flag) == string::npos && !(file_stream).eof()) {
    GetStrippedLine(file_stream, next_line);
  }

  /*--- throw error if end of file reached ---*/
  if ((file_stream).eof()) SU2_MPI::Error("Flag " + flag + " not found in file", CURRENT_FUNCTION);
}

bool CFileReaderLUT::GetNextNonEmptyLine(ifstream& file_stream, string& line) const {
  /*--- get next line and save return value ---*/
  bool return_value = GetStrippedLine(file_stream, line);

  /*--- skip empty lines ---*/
  while (line.empty() && !(file_stream).eof()) {
    return_value = GetStrippedLine(file_stream, line);
  }

  /*--- return true if line is not empty, else return false ---*/
  return return_value;
}

bool CFileReaderLUT::GetStrippedLine(ifstream& file_stream, string& line) const {
  /*--- get next line and save return value ---*/
  getline(file_stream, line);

  /*--- find last non-control-character character ---*/
  size_t end = line.find_last_not_of(" \n\r\t\f\v");

  /*--- clean up line ---*/
  line = (end == string::npos) ? "" : line.substr(0, end + 1);

  /*--- return true if line is not empty, else return false ---*/
  return !line.empty();
}
