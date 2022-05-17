/*!
 * \file CFileReaderLUT.hpp
 * \brief reading lookup table for tabulated fluid properties
 * \author D. Mayer, T. Economon
 * \version 7.3.1 "Blackbird"
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
#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "../../Common/include/parallelization/mpi_structure.hpp"

using namespace std;

class CFileReaderLUT {
 protected:
  int rank;

  string type_lut;
  string version_lut;
  string version_reader;
  unsigned long n_points;
  unsigned long n_triangles;
  unsigned long n_hull_points;
  unsigned long n_variables;

  /* !brief
   * Holds the variable names stored in the table file. Order is in sync with
   * tableFlamelet.
   */
  vector<string> names_var;

  /* !brief
   * Holds all data stored in the table. First index addresses the variable
   * while second index addresses the point.
   */
  vector<vector<su2double> > table_data;

  vector<vector<unsigned long> > triangles;

  vector<unsigned long> hull;

  string SkipToFlag(ifstream* file_stream, string flag);

  inline void SetTypeLUT(string value) { type_lut = value; }
  inline void SetVersionLUT(string value) { version_lut = value; }
  inline void SetNPoints(unsigned long value) { n_points = value; }
  inline void SetNTriangles(unsigned long value) { n_triangles = value; }
  inline void SetNHullPoints(unsigned long value) { n_hull_points = value; }
  inline void SetNVariables(unsigned long value) { n_variables = value; }

  inline void PushNameVar(string value) { names_var.push_back(value); }
  inline void PopNameVar() { names_var.pop_back(); }

  inline void AllocMemData() { table_data.resize(GetNVariables(), vector<su2double>(GetNPoints())); }

  inline void AllocMemTriangles() { triangles.resize(GetNTriangles(), vector<unsigned long>(3)); }

  inline void AllocMemHull() { hull.resize(GetNHullPoints()); }

 public:
  CFileReaderLUT();

  inline string GetTypeLUT() { return type_lut; }
  inline string GetVersionLUT() { return version_lut; }
  inline string GetVersionReader() { return version_reader; }
  inline unsigned long GetNPoints() { return n_points; }
  inline unsigned long GetNTriangles() { return n_triangles; }
  inline unsigned long GetNHullPoints() { return n_hull_points; }
  inline unsigned long GetNVariables() { return n_variables; }

  inline const vector<string>& GetNamesVar() const { return names_var; }

  inline const vector<vector<su2double> >& GetTableData() const { return table_data; }

  inline const vector<vector<unsigned long> >& GetTriangles() const { return triangles; };

  inline const vector<unsigned long>& GetHull() const { return hull; };

  void ReadRawDRG(string file_name);
};
