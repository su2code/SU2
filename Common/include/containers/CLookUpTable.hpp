/*!
 * \file CLookupTable.hpp
 * \brief tabulation of fluid properties
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
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../../Common/include/option_structure.hpp"
#include "../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../Common/include/toolboxes/CSquareMatrixCM.hpp"
#include "CFileReaderLUT.hpp"
#include "CTrapezoidalMap.hpp"

using namespace std;

class CLookUpTable {
 protected:
  int rank; /*!< \brief MPI Rank. */

  string file_name_lut;
  string type_lut;
  string version_lut;
  string version_reader;
  unsigned long n_points;
  unsigned long n_triangles;
  unsigned long n_variables;
  unsigned long n_hull_points;

  su2double limits_table_enth[2];
  su2double limits_table_prog[2];

  /*! \brief Holds the variable names stored in the table file.
   * Order is in sync with data
   */
  vector<string> names_var;

  /*! \brief
   * Holds all data stored in the table. First index addresses the variable
   * while second index addresses the point.
   */
  su2activematrix table_data;

  su2matrix<unsigned long> triangles;

  /* we do not know this size in advance until we go through the entire lookup table */
  vector<vector<unsigned long> > edges;
  vector<vector<unsigned long> > edge_to_triangle;

  vector<unsigned long> hull;

  CTrapezoidalMap trap_map_prog_enth;

  vector<su2activematrix> interp_mat_inv_prog_enth;

  inline int GetIndexOfVar(string nameVar) {
    int index = find(names_var.begin(), names_var.end(), nameVar) - names_var.begin();
    if (index == int(names_var.size())) {
      index = -1;
      string error_msg = "Variable '";
      error_msg.append(nameVar);
      error_msg.append("' is not in the lookup table.");
      SU2_MPI::Error(error_msg, CURRENT_FUNCTION);
    }
    return index;
  }

/* get the pointer to the column data of the table (density, temperature, source terms, ...) */
inline const su2double* GetDataP(string name_var) {
    int ix_var = GetIndexOfVar(name_var);

    su2double* tableDataRow(table_data[ix_var]); 

    return tableDataRow;
  }

  /*!
   * \brief find the table limits, i.e. the minimum and maximum values of the 2 independent
   * controlling variables (progress variable and enthalpy). We put the values in the variables
   * limits_table_prog[2] and limit_table_enth[2]. 
   * \param[in] name_prog - the string name for the first controlling variable
   * \param[in] name_enth - the string name of the second controlling variable 
   */
  void FindTableLimits(string name_prog, string name_enth);

  void IdentifyUniqueEdges();

  void LoadTableRaw(string file_name_lut);

  void ComputeInterpCoeffs(string name_prog, string name_enth);

  void GetInterpMatInv(const su2double* vec_x, const su2double* vec_y, std::array<unsigned long,3>& point_ids,
                       su2activematrix& interp_mat_inv);

  void GetInterpCoeffs(su2double val_x, su2double val_y, su2activematrix& interp_mat_inv,
                       std::array<su2double,3>& interp_coeffs);

  su2double Interpolate(const su2double* val_samples, std::array<unsigned long,3>& val_triangle,
                        std::array<su2double,3>& val_interp_coeffs);

  unsigned long FindNearestNeighborOnHull(su2double val_enth, su2double val_prog, string name_prog, string name_enth);

  bool IsInTriangle(su2double val_x, su2double val_y, unsigned long val_id_triangle, string name_prog,
                    string name_enth);

  inline su2double TriArea(su2double x1, su2double y1, su2double x2, su2double y2, su2double x3, su2double y3) {
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) * 0.5);
  }

 public:
  CLookUpTable(string file_name_lut, string name_prog, string name_enth);

  void PrintTableInfo();

  unsigned long LookUp_ProgEnth(string val_name_var, su2double* val_var, su2double val_prog, su2double val_enth,
                                string name_prog, string name_enth);

  unsigned long LookUp_ProgEnth(vector<string>& val_names_var, vector<su2double*>& val_vars, su2double val_prog,
                                su2double val_enth, string name_prog, string name_enth);

  unsigned long LookUp_ProgEnth(vector<string>& val_names_var, vector<su2double>& val_vars, su2double val_prog,
                                su2double val_enth, string name_prog, string name_enth);

  inline pair<su2double, su2double> GetTableLimitsEnth() {
    return make_pair(limits_table_enth[0], limits_table_enth[1]);
  }

  inline pair<su2double, su2double> GetTableLimitsProg() {
    return make_pair(limits_table_prog[0], limits_table_prog[1]);
  }
};
