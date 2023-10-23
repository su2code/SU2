/*!
 * \file CLookupTable.cpp
 * \brief tabulation of fluid properties
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

#include <utility>

#include "../../../Common/include/containers/CLookUpTable.hpp"

#include "../../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../../Common/include/toolboxes/CSquareMatrixCM.hpp"

using namespace std;

CLookUpTable::CLookUpTable(const string& var_file_name_lut, string name_CV1_in, string name_CV2_in)
    : file_name_lut{var_file_name_lut}, name_CV1{std::move(name_CV1_in)}, name_CV2{std::move(name_CV2_in)} {
  rank = SU2_MPI::GetRank();

  LoadTableRaw(var_file_name_lut);

  FindTableLimits(name_CV1, name_CV2);

  if (rank == MASTER_NODE)
    cout << "Detecting all unique edges and setting edge to triangle connectivity "
            "..."
         << endl;

  IdentifyUniqueEdges();

  if (rank == MASTER_NODE) cout << " done." << endl;

  PrintTableInfo();

  if (rank == MASTER_NODE) switch (table_dim) {
      case 2:
        cout << "Building a trapezoidal map for the (" + name_CV1 + ", " + name_CV2 +
                    ") "
                    "space ..."
             << endl;
        break;
      case 3:
        cout << "Building trapezoidal map stack for the (" + name_CV1 + ", " + name_CV2 +
                    ") "
                    "space ..."
             << endl;
        break;
      default:
        break;
    }

  trap_map_x_y.resize(n_table_levels);
  su2double startTime = SU2_MPI::Wtime();
  unsigned short barwidth = 65;
  bool display_map_info = (n_table_levels < 2);
  double tmap_memory_footprint = 0;
  for (auto i_level = 0ul; i_level < n_table_levels; i_level++) {
    trap_map_x_y[i_level] =
        CTrapezoidalMap(GetDataP(name_CV1, i_level), GetDataP(name_CV2, i_level), table_data[i_level].cols(),
                        edges[i_level], edge_to_triangle[i_level], display_map_info);
    tmap_memory_footprint += trap_map_x_y[i_level].GetMemoryFootprint();
    /* Display a progress bar to monitor table generation process */
    if (rank == MASTER_NODE) {
      su2double progress = su2double(i_level) / n_table_levels;
      auto completed = floor(progress * barwidth);
      auto to_do = barwidth - completed;
      cout << "[" << setfill('=') << setw(completed);
      cout << '>';
      cout << setfill(' ') << setw(to_do) << std::right << "] " << 100 * progress << "%\r";
      cout.flush();
    }
  }
  su2double stopTime = SU2_MPI::Wtime();

  if (rank == MASTER_NODE) {
    switch (table_dim) {
      case 2:
        cout << "\nConstruction of trapezoidal map took " << stopTime - startTime << " seconds\n" << endl;
        break;
      case 3:
        cout << "\nConstruction of trapezoidal map stack took " << stopTime - startTime << " seconds\n" << endl;
        break;
      default:
        break;
    }
    cout << "Trapezoidal map memory footprint: " << tmap_memory_footprint << " MB\n";
    cout << "Table data memory footprint: " << memory_footprint_data << " MB\n" << endl;
  }

  ComputeInterpCoeffs();

  if (rank == MASTER_NODE) cout << "LUT fluid model ready for use" << endl;
}

void CLookUpTable::LoadTableRaw(const string& var_file_name_lut) {
  CFileReaderLUT file_reader;

  if (rank == MASTER_NODE) cout << "Loading lookup table, filename = " << var_file_name_lut << " ..." << endl;

  file_reader.ReadRawLUT(var_file_name_lut);
  table_dim = file_reader.GetTableDim();
  n_table_levels = file_reader.GetNLevels();

  n_points.resize(n_table_levels);
  n_triangles.resize(n_table_levels);
  n_hull_points.resize(n_table_levels);
  table_data.resize(n_table_levels);
  hull.resize(n_table_levels);
  triangles.resize(n_table_levels);
  interp_mat_inv_x_y.resize(n_table_levels);
  edges.resize(n_table_levels);
  edge_to_triangle.resize(n_table_levels);

  for (unsigned long i_level = 0; i_level < n_table_levels; i_level++) {
    n_points[i_level] = file_reader.GetNPoints(i_level);
    n_triangles[i_level] = file_reader.GetNTriangles(i_level);
    n_hull_points[i_level] = file_reader.GetNHullPoints(i_level);
    table_data[i_level] = file_reader.GetTableData(i_level);
    triangles[i_level] = file_reader.GetTriangles(i_level);
    hull[i_level] = file_reader.GetHull(i_level);
    memory_footprint_data += n_points[i_level] * sizeof(su2double);
  }
  memory_footprint_data /= 1e6;

  n_variables = file_reader.GetNVariables();
  version_lut = file_reader.GetVersionLUT();
  version_reader = file_reader.GetVersionReader();
  names_var = file_reader.GetNamesVar();

  if (table_dim == 3) {
    z_values_levels.resize(n_table_levels);
    for (unsigned long i_level = 0; i_level < n_table_levels; i_level++) {
      z_values_levels[i_level] = file_reader.GetTableLevel(i_level);
    }
  }
  if (rank == MASTER_NODE) cout << " done." << endl;
}

void CLookUpTable::FindTableLimits(const string& name_cv1, const string& name_cv2) {
  int idx_CV1 = GetIndexOfVar(name_cv1);
  int idx_CV2 = GetIndexOfVar(name_cv2);
  limits_table_x.resize(n_table_levels);
  limits_table_y.resize(n_table_levels);

  /* we find the lowest and highest value of y and x in the table */
  for (auto i_level = 0u; i_level < n_table_levels; i_level++) {
    limits_table_y[i_level] =
        minmax_element(table_data[i_level][idx_CV2], table_data[i_level][idx_CV2] + table_data[i_level].cols());
    limits_table_x[i_level] =
        minmax_element(table_data[i_level][idx_CV1], table_data[i_level][idx_CV1] + table_data[i_level].cols());
  }

  if (table_dim == 3) {
    limits_table_z = minmax_element(z_values_levels.data(), z_values_levels.data() + z_values_levels.size());
  }
}

void CLookUpTable::PrintTableInfo() {
  if (rank == MASTER_NODE) {
    cout << setfill(' ');
    cout << endl;
    cout << "+------------------------------------------------------------------+\n";
    cout << "|                     Look-Up-Table (LUT) info                     |\n";
    cout << "+------------------------------------------------------------------+" << endl;
    cout << "| File name:" << setw(54) << right << file_name_lut << " |" << endl;
    cout << "| Table version:" << setw(50) << right << version_lut << " |" << endl;
    cout << "| Table reader version:" << setw(43) << right << version_reader << " |" << endl;
    cout << "| Table dimension:" << setw(48) << right << table_dim << " |" << endl;
    cout << "| Number of variables:" << setw(44) << right << n_variables << " |" << endl;

    su2double n_points_av = 0, n_tria_av = 0, n_edges_av = 0, min_x = 1 / EPS, max_x = -1 / EPS, min_y = 1 / EPS,
              max_y = -1 / EPS;
    for (auto i_level = 0ul; i_level < n_table_levels; i_level++) {
      n_points_av += n_points[i_level] / n_table_levels;
      n_tria_av += n_triangles[i_level] / n_table_levels;
      n_edges_av += edges[i_level].size() / n_table_levels;
      min_x = min(min_x, *limits_table_x[i_level].first);
      min_y = min(min_y, *limits_table_y[i_level].first);
      max_x = max(max_x, *limits_table_x[i_level].second);
      max_y = max(max_y, *limits_table_y[i_level].second);
    }
    switch (table_dim) {
      case 2:
        cout << "| Number of points:" << setw(47) << right << n_points_av << " |" << endl;
        cout << "| Number of triangles:" << setw(44) << right << n_tria_av << " |" << endl;
        cout << "| Number of edges:" << setw(48) << right << n_edges_av << " |" << endl;
        cout << "+------------------------------------------------------------------+" << endl;
        cout << "| Minimum Y:" << setw(54) << right << min_y << " |" << endl;
        cout << "| Maximum Y:" << setw(54) << right << max_y << " |" << endl;
        cout << "| Minimum X:" << setw(54) << right << min_x << " |" << endl;
        cout << "| Maximum X:" << setw(54) << right << max_x << " |" << endl;
        cout << "+------------------------------------------------------------------+" << endl;
        cout << "| Variable names:                                                  |" << endl;
        cout << "|                                                                  |" << endl;

        for (unsigned long i_var = 0; i_var < names_var.size(); i_var++)
          cout << "| " << right << setw(3) << i_var << ": " << left << setw(59) << names_var[i_var] << " |" << endl;

        cout << "+------------------------------------------------------------------+" << endl;

        cout << endl;
        break;
      case 3:

        cout << "| Number of table levels:" << setw(41) << right << n_table_levels << " |" << endl;
        cout << "| Average number of points:" << setw(39) << right << floor(n_points_av) << " |" << endl;
        cout << "| Average number of triangles:" << setw(36) << right << floor(n_tria_av) << " |" << endl;
        cout << "| Average number of edges:" << setw(40) << right << floor(n_edges_av) << " |" << endl;
        cout << "+------------------------------------------------------------------+" << endl;
        cout << "| Minimum Z:" << setw(54) << right << *limits_table_z.first << " |" << endl;
        cout << "| Maximum Z:" << setw(54) << right << *limits_table_z.second << " |" << endl;
        cout << "| Minimum Y:" << setw(54) << right << min_y << " |" << endl;
        cout << "| Maximum Y:" << setw(54) << right << max_y << " |" << endl;
        cout << "| Minimum X:" << setw(54) << right << min_x << " |" << endl;
        cout << "| Maximum X:" << setw(54) << right << max_x << " |" << endl;
        cout << "+------------------------------------------------------------------+" << endl;
        cout << "| Variable names:                                                  |" << endl;
        cout << "|                                                                  |" << endl;

        for (unsigned long i_var = 0; i_var < names_var.size(); i_var++)
          cout << "| " << right << setw(3) << i_var << ": " << left << setw(59) << names_var[i_var] << " |" << endl;

        cout << "+------------------------------------------------------------------+" << endl;
        break;
      default:
        break;
    }
  }
}

void CLookUpTable::IdentifyUniqueEdges() {
  for (auto i_level = 0u; i_level < n_table_levels; i_level++) {
    /* loop through elements and store the vector of element IDs (neighbors)
       for each of the points in the table */

    vector<vector<unsigned long> > neighborElemsOfPoint;

    neighborElemsOfPoint.resize(n_points[i_level]);
    for (unsigned long iElem = 0; iElem < n_triangles[i_level]; iElem++) {
      /* loop over 3 points per triangle */
      for (unsigned long iPoint = 0; iPoint < N_POINTS_TRIANGLE; iPoint++) {
        /* get the global ID of the current point */
        const unsigned long GlobalIndex = triangles[i_level][iElem][iPoint];

        /* add the current element ID to the neighbor list for this point */
        neighborElemsOfPoint[GlobalIndex].push_back(iElem);
      }
    }

    /* remove duplicates from the neighboring element lists*/
    vector<unsigned long>::iterator vecIt;
    for (unsigned long iPoint = 0; iPoint < n_points[i_level]; iPoint++) {
      /* sort neighboring elements for each point */
      sort(neighborElemsOfPoint[iPoint].begin(), neighborElemsOfPoint[iPoint].end());

      /* uniquify list of neighboring elements */
      vecIt = unique(neighborElemsOfPoint[iPoint].begin(), neighborElemsOfPoint[iPoint].end());

      /* adjust size of vector */
      neighborElemsOfPoint[iPoint].resize(vecIt - neighborElemsOfPoint[iPoint].begin());
    }

    /* loop through all neighboring elements of each point and store
       the point IDs that are neighboring points */
    vector<vector<unsigned long> > neighborPointsOfPoint;
    neighborPointsOfPoint.resize(n_points[i_level]);
    for (unsigned long iPoint = 0; iPoint < n_points[i_level]; iPoint++) {
      for (unsigned long iElem = 0; iElem < neighborElemsOfPoint[iPoint].size(); iElem++) {
        /* loop over element points */
        for (unsigned long jPoint = 0; jPoint < N_POINTS_TRIANGLE; jPoint++) {
          /* get the global ID of the current point */
          const unsigned long GlobalIndex = triangles[i_level][neighborElemsOfPoint[iPoint][iElem]][jPoint];

          /* add the current element ID to the neighbor list for this point */
          if (GlobalIndex != iPoint) neighborPointsOfPoint[iPoint].push_back(GlobalIndex);
        }
      }
    }

    /* remove duplicates from the neighboring points lists */
    for (unsigned long iPoint = 0; iPoint < n_points[i_level]; iPoint++) {
      /* sort neighboring points for each point */
      sort(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());

      /* uniquify list of neighboring elements */
      vecIt = unique(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());

      /* adjust size of vector */
      neighborPointsOfPoint[iPoint].resize(vecIt - neighborPointsOfPoint[iPoint].begin());
    }

    /* loop through our point neighbors and fill the vector of the unique
       point pairs making up each edge in the grid. We impose a requirement
       that the smaller global index is in the first position and the larger
       is in the second to make for a unique set, so there's no need to
       remove duplicates. */
    for (unsigned long iPoint = 0; iPoint < n_points[i_level]; iPoint++) {
      for (unsigned long jPoint = 0; jPoint < neighborPointsOfPoint[iPoint].size(); jPoint++) {
        /* Store the neighbor index more clearly. */
        const unsigned long GlobalIndex = neighborPointsOfPoint[iPoint][jPoint];

        /* Store the edge so that the lower index of the pair is always first. */
        if (iPoint < GlobalIndex) {
          std::array<unsigned long, 2> edge;
          edge[0] = iPoint;
          edge[1] = GlobalIndex;
          edges[i_level].push_back(edge);
        }
      }
    }

    /* Loop over our edges data structure. For the first point in each
     pair, loop through the neighboring elements and store the two
     elements that contain the second point in the edge ('left' and 'right' of the edge). */
    edge_to_triangle[i_level].resize(edges[i_level].size());
    for (unsigned long iEdge = 0; iEdge < edges[i_level].size(); iEdge++) {
      /* Store the two points of the edge more clearly. */
      const unsigned long iPoint = edges[i_level][iEdge][0];
      const unsigned long jPoint = edges[i_level][iEdge][1];

      /* Loop over all neighboring elements to iPoint. */
      for (unsigned long iElem = 0; iElem < neighborElemsOfPoint[iPoint].size(); iElem++) {
        /* loop over 3 points per triangle */
        for (unsigned long kPoint = 0; kPoint < N_POINTS_TRIANGLE; kPoint++) {
          /* Get the global ID of the current point. */
          const unsigned long GlobalIndex = triangles[i_level][neighborElemsOfPoint[iPoint][iElem]][kPoint];

          /* Add the current element ID to the neighbor list for this point. */
          if (GlobalIndex == jPoint) edge_to_triangle[i_level][iEdge].push_back(neighborElemsOfPoint[iPoint][iElem]);
        }
      }
    }
  }
}

void CLookUpTable::ComputeInterpCoeffs() {
  for (auto i_level = 0ul; i_level < n_table_levels; i_level++) {
    /* build KD tree for y, x space */

    std::array<unsigned long, 3> next_triangle;

    const su2double* val_CV1 = GetDataP(name_CV1, i_level);
    const su2double* val_CV2 = GetDataP(name_CV2, i_level);

    /* calculate weights for each triangle (basically a distance function) and
     * build inverse interpolation matrices */
    interp_mat_inv_x_y[i_level].resize(n_triangles[i_level]);
    for (unsigned long i_triangle = 0; i_triangle < n_triangles[i_level]; i_triangle++) {
      for (int p = 0; p < 3; p++) {
        next_triangle[p] = triangles[i_level][i_triangle][p];
      }

      su2activematrix x_interp_mat_inv(3, 3);
      GetInterpMatInv(val_CV1, val_CV2, next_triangle, x_interp_mat_inv);
      interp_mat_inv_x_y[i_level][i_triangle] = x_interp_mat_inv;
    }
  }
}

void CLookUpTable::GetInterpMatInv(const su2double* vec_x, const su2double* vec_y,
                                   std::array<unsigned long, 3>& point_ids, su2activematrix& interp_mat_inv) {
  const unsigned int M = 3;
  CSquareMatrixCM global_M(3);

  /* setup LHM matrix for the interpolation */
  for (unsigned int i_point = 0; i_point < M; i_point++) {
    su2double x = vec_x[point_ids[i_point]];
    su2double y = vec_y[point_ids[i_point]];

    global_M(i_point, 0) = SU2_TYPE::GetValue(1.0);
    global_M(i_point, 1) = SU2_TYPE::GetValue(x);
    global_M(i_point, 2) = SU2_TYPE::GetValue(y);
  }

  global_M.Invert();
  global_M.Transpose();

  for (unsigned int i = 0; i < M; i++) {
    for (unsigned int j = 0; j < M; j++) {
      interp_mat_inv[i][j] = global_M(i, j);
    }
  }
}

std::pair<unsigned long, unsigned long> CLookUpTable::FindInclusionLevels(const su2double val_CV3) {
  /*--- Find the table levels with constant z-values directly below and above the query value val_CV3 ---*/

  /* Check if val_CV3 lies outside table bounds */
  if (val_CV3 >= *limits_table_z.second) {
    return std::make_pair(n_table_levels - 1, n_table_levels - 1);
  } else if (val_CV3 <= *limits_table_z.first) {
    return std::make_pair(0, 0);
  } else {
    std::pair<unsigned long, unsigned long> inclusion_levels;
    /* Loop over table levels to find the levels directly above and below the query value */
    for (auto i_level = 0ul; i_level < n_table_levels - 1; i_level++) {
      if ((val_CV3 >= z_values_levels[i_level]) && (val_CV3 < z_values_levels[i_level + 1])) {
        inclusion_levels = std::make_pair(i_level, i_level + 1);
      }
    }
    return inclusion_levels;
  }
}

unsigned long CLookUpTable::LookUp_XYZ(const std::string& val_name_var, su2double* val_var, su2double val_CV1,
                                       su2double val_CV2, su2double val_CV3) {
  /*--- Perform quasi-3D interpolation for a single variable named val_name_var on a query point
        with coordinates val_CV1, val_CV2, and val_CV3 ---*/

  /* 1: Find table levels directly above and below the query point (the levels that sandwhich val_CV3) */
  std::pair<unsigned long, unsigned long> inclusion_levels = FindInclusionLevels(val_CV3);
  bool within_z_limits = (inclusion_levels.first != inclusion_levels.second);

  if (within_z_limits) {
    /* 2: Determine val_CV1 and val_CV2 for the inclusion table levels. */
    std::array<std::array<su2double, 2>, 2> lower_upper_CV1_2 =
        ComputeNormalizedXY(inclusion_levels, val_CV1, val_CV2, val_CV3);
    su2double val_CV1_lower = lower_upper_CV1_2[0][0], val_CV2_lower = lower_upper_CV1_2[0][1],
              val_CV1_upper = lower_upper_CV1_2[1][0], val_CV2_upper = lower_upper_CV1_2[1][1];

    /* 3: Perform 2D interpolations on upper and lower inclusion levels */
    unsigned long lower_level = inclusion_levels.first;
    unsigned long upper_level = inclusion_levels.second;

    su2double val_var_lower, val_var_upper;
    unsigned long exit_code_lower = LookUp_XY(val_name_var, &val_var_lower, val_CV1_lower, val_CV2_lower, lower_level);
    unsigned long exit_code_upper = LookUp_XY(val_name_var, &val_var_upper, val_CV1_upper, val_CV2_upper, upper_level);

    /* 4: Perform linear interpolation along the z-direction using the x-y interpolation results
             from upper and lower trapezoidal maps */
    Linear_Interpolation(val_CV3, lower_level, upper_level, val_var_lower, val_var_upper, val_var);

    return max(exit_code_lower, exit_code_upper);
  } else {
    /* Perform single, 2D interpolation when val_CV3 lies outside table bounds */
    unsigned long bound_level = inclusion_levels.first;
    LookUp_XY(val_name_var, val_var, val_CV1, val_CV2, bound_level);
    return 1;
  }
}
unsigned long CLookUpTable::LookUp_XYZ(const std::vector<std::string>& val_names_var, std::vector<su2double>& val_vars,
                                       su2double val_CV1, su2double val_CV2, su2double val_CV3) {
  /*--- Perform quasi-3D interpolation for a vector of variables with names val_names_var
        on a query point with coordinates val_CV1, val_CV2, and val_CV3 ---*/

  /* 1: Find table levels directly above and below the query point (the levels that sandwhich val_CV3) */
  std::pair<unsigned long, unsigned long> inclusion_levels = FindInclusionLevels(val_CV3);
  bool within_z_limits = (inclusion_levels.first != inclusion_levels.second);

  if (within_z_limits) {
    /* 2: Determine val_CV1 and val_CV2 for the inclusion table levels. */
    std::array<std::array<su2double, 2>, 2> lower_upper_CV1_2 =
        ComputeNormalizedXY(inclusion_levels, val_CV1, val_CV2, val_CV3);
    su2double val_CV1_lower = lower_upper_CV1_2[0][0], val_CV2_lower = lower_upper_CV1_2[0][1],
              val_CV1_upper = lower_upper_CV1_2[1][0], val_CV2_upper = lower_upper_CV1_2[1][1];

    /* 3: Perform 2D interpolations on upper and lower inclusion levels */
    unsigned long lower_level = inclusion_levels.first;
    unsigned long upper_level = inclusion_levels.second;

    std::vector<su2double> val_vars_lower, val_vars_upper;
    val_vars_lower.resize(val_vars.size());
    val_vars_upper.resize(val_vars.size());
    unsigned long exit_code_lower = LookUp_XY(val_names_var, val_vars_lower, val_CV1_lower, val_CV2_lower, lower_level);
    unsigned long exit_code_upper = LookUp_XY(val_names_var, val_vars_upper, val_CV1_upper, val_CV2_upper, upper_level);

    /* 4: Perform linear interpolation along the z-direction using the x-y interpolation results
             from upper and lower trapezoidal maps */
    Linear_Interpolation(val_CV3, lower_level, upper_level, val_vars_lower, val_vars_upper, val_vars);

    return max(exit_code_lower, exit_code_upper);
  } else {
    /* Perform single, 2D interpolation when val_CV3 lies outside table bounds */
    unsigned long bound_level = inclusion_levels.first;
    LookUp_XY(val_names_var, val_vars, val_CV1, val_CV2, bound_level);
    return 1;
  }
}

void CLookUpTable::Linear_Interpolation(const su2double val_CV3, const unsigned long lower_level,
                                        const unsigned long upper_level, su2double& lower_value, su2double& upper_value,
                                        su2double*& var_val) const {
  /* Perform linear interpolation along the z-direction of the table for a single variable */

  /* Retrieve constant z-values of inclusion levels */
  su2double val_z_lower = z_values_levels[lower_level];
  su2double val_z_upper = z_values_levels[upper_level];

  /* Compute linear interpolation coefficients */
  su2double factor_upper = (val_CV3 - val_z_lower) / (val_z_upper - val_z_lower);
  su2double factor_lower = (val_z_upper - val_CV3) / (val_z_upper - val_z_lower);

  /* Perform linear interpolation */
  *var_val = lower_value * factor_lower + upper_value * factor_upper;
}

void CLookUpTable::Linear_Interpolation(const su2double val_CV3, const unsigned long lower_level,
                                        const unsigned long upper_level, std::vector<su2double>& lower_values,
                                        std::vector<su2double>& upper_values, std::vector<su2double>& var_vals) const {
  /* Perform linear interpolation along the z-direction of the table for multiple variables */

  /* Retrieve constant z-values of inclusion levels */
  su2double val_z_lower = z_values_levels[lower_level];
  su2double val_z_upper = z_values_levels[upper_level];

  /* Compute linear interpolation coefficients */
  su2double factor_upper = (val_CV3 - val_z_lower) / (val_z_upper - val_z_lower);
  su2double factor_lower = (val_z_upper - val_CV3) / (val_z_upper - val_z_lower);

  /* Perform linear interpolation */
  for (size_t iVar = 0; iVar < var_vals.size(); iVar++) {
    var_vals[iVar] = lower_values[iVar] * factor_lower + upper_values[iVar] * factor_upper;
  }
}

std::array<std::array<su2double, 2>, 2> CLookUpTable::ComputeNormalizedXY(
    std::pair<unsigned long, unsigned long>& inclusion_levels, const su2double val_CV1, const su2double val_CV2,
    const su2double val_CV3) {
  /* Determine the values of the controlling variables of the upper and lower table levels corresponding to the
     normalized query point coordinates. */

  unsigned long lower_level = inclusion_levels.first;
  unsigned long upper_level = inclusion_levels.second;

  /* Get the upper and lower limits of the first and second controlling variables on the inclusion levels. */
  su2double x_min_lower = *limits_table_x[lower_level].first, x_max_lower = *limits_table_x[lower_level].second;
  su2double x_min_upper = *limits_table_x[upper_level].first, x_max_upper = *limits_table_x[upper_level].second;
  su2double y_min_lower = *limits_table_y[lower_level].first, y_max_lower = *limits_table_y[lower_level].second;
  su2double y_min_upper = *limits_table_y[upper_level].first, y_max_upper = *limits_table_y[upper_level].second;

  /* Interpolate the table limits to the query point coordinates. */
  su2double x_min_q = x_min_lower + (x_min_upper - x_min_lower) * (val_CV3 - z_values_levels[lower_level]) /
                                        (z_values_levels[upper_level] - z_values_levels[lower_level] + EPS);
  su2double x_max_q = x_max_lower + (x_max_upper - x_max_lower) * (val_CV3 - z_values_levels[lower_level]) /
                                        (z_values_levels[upper_level] - z_values_levels[lower_level] + EPS);
  su2double y_min_q = y_min_lower + (y_min_upper - y_min_lower) * (val_CV3 - z_values_levels[lower_level]) /
                                        (z_values_levels[upper_level] - z_values_levels[lower_level] + EPS);
  su2double y_max_q = y_max_lower + (y_max_upper - y_max_lower) * (val_CV3 - z_values_levels[lower_level]) /
                                        (z_values_levels[upper_level] - z_values_levels[lower_level] + EPS);

  /* Compute the normalized first and second controlling variable values. */
  su2double x_norm_q = (val_CV1 - x_min_q) / (x_max_q - x_min_q);
  su2double y_norm_q = (val_CV2 - y_min_q) / (y_max_q - y_min_q);

  /* De-normalize the normalized coordinates for the upper and lower table levels. */
  su2double val_CV1_upper = x_min_upper + (x_max_upper - x_min_upper) * x_norm_q;
  su2double val_CV2_upper = y_min_upper + (y_max_upper - y_min_upper) * y_norm_q;
  su2double val_CV1_lower = x_min_lower + (x_max_lower - x_min_lower) * x_norm_q;
  su2double val_CV2_lower = y_min_lower + (y_max_lower - y_min_lower) * y_norm_q;

  /* Store coordinates in output */
  std::array<std::array<su2double, 2>, 2> lower_upper_CVs;
  lower_upper_CVs[0][0] = val_CV1_lower;
  lower_upper_CVs[0][1] = val_CV2_lower;
  lower_upper_CVs[1][0] = val_CV1_upper;
  lower_upper_CVs[1][1] = val_CV2_upper;

  return lower_upper_CVs;
}

unsigned long CLookUpTable::LookUp_XY(const string& val_name_var, su2double* val_var, su2double val_CV1,
                                      su2double val_CV2, unsigned long i_level) {
  unsigned long exit_code = 1;

  if (noSource(val_name_var)) {
    *val_var = 0.0;
    exit_code = 0;
    return exit_code;
  }

  /* check if x and y value is in table range */
  if ((val_CV1 >= *limits_table_x[i_level].first && val_CV1 <= *limits_table_x[i_level].second) &&
      (val_CV2 >= *limits_table_y[i_level].first && val_CV2 <= *limits_table_y[i_level].second)) {
    /* find the triangle that holds the (x, y) point */
    unsigned long id_triangle = trap_map_x_y[i_level].GetTriangle(val_CV1, val_CV2);

    if (IsInTriangle(val_CV1, val_CV2, id_triangle, i_level)) {
      /* get interpolation coefficients for point on triangle */
      std::array<su2double, 3> interp_coeffs{0};
      GetInterpCoeffs(val_CV1, val_CV2, interp_mat_inv_x_y[i_level][id_triangle], interp_coeffs);

      /* first, copy the single triangle from the large triangle list*/
      std::array<unsigned long, 3> triangle{0};
      for (int p = 0; p < 3; p++) triangle[p] = triangles[i_level][id_triangle][p];

      *val_var = Interpolate(GetDataP(val_name_var, i_level), triangle, interp_coeffs);
      exit_code = 0;
    }
  }
  if (exit_code == 1) InterpolateToNearestNeighbors(val_CV1, val_CV2, val_name_var, val_var, i_level);

  return exit_code;
}

unsigned long CLookUpTable::LookUp_XY(const vector<string>& val_names_var, vector<su2double>& val_vars,
                                      su2double val_CV1, su2double val_CV2, unsigned long i_level) {
  vector<su2double*> look_up_data(val_names_var.size());

  for (long unsigned int i_var = 0; i_var < val_vars.size(); ++i_var) {
    look_up_data[i_var] = &val_vars[i_var];
  }

  unsigned long exit_code = LookUp_XY(val_names_var, look_up_data, val_CV1, val_CV2, i_level);

  return exit_code;
}

unsigned long CLookUpTable::LookUp_XY(const vector<string>& val_names_var, vector<su2double*>& val_vars,
                                      su2double val_CV1, su2double val_CV2, unsigned long i_level) {
  unsigned long exit_code = 1;
  unsigned long id_triangle = 0;
  std::array<su2double, 3> interp_coeffs{0};
  std::array<unsigned long, 3> triangle{0};

  /* check if x value is in table x-dimension range
   * and if y is in table y-dimension table range */
  if ((val_CV1 >= *limits_table_x[i_level].first && val_CV1 <= *limits_table_x[i_level].second) &&
      (val_CV2 >= *limits_table_y[i_level].first && val_CV2 <= *limits_table_y[i_level].second)) {
    /* if so, try to find the triangle that holds the (prog, enth) point */
    id_triangle = trap_map_x_y[i_level].GetTriangle(val_CV1, val_CV2);

    /* check if point is inside a triangle (if table domain is non-rectangular,
     * the previous range check might be true but the point could still be outside of the domain) */
    if (IsInTriangle(val_CV1, val_CV2, id_triangle, i_level)) {
      /* if so, get interpolation coefficients for point in  the triangle */
      GetInterpCoeffs(val_CV1, val_CV2, interp_mat_inv_x_y[i_level][id_triangle], interp_coeffs);

      /* exit_code 0 means point was in triangle */
      exit_code = 0;
    }
  }

  /* loop over variable names and interpolate / get values */
  for (long unsigned int i_var = 0; i_var < val_names_var.size(); ++i_var) {
    if (noSource(val_names_var[i_var])) {
      *val_vars[i_var] = 0.0;
      exit_code = 0;
    } else {
      if (exit_code == 0) {
        /* first, copy the single triangle from the large triangle list*/
        for (int p = 0; p < 3; p++) triangle[p] = triangles[i_level][id_triangle][p];
        *val_vars[i_var] = Interpolate(GetDataP(val_names_var[i_var], i_level), triangle, interp_coeffs);
      } else {
        InterpolateToNearestNeighbors(val_CV1, val_CV2, val_names_var[i_var], val_vars[i_var], i_level);
      }
    }
  }

  return exit_code;
}

void CLookUpTable::GetInterpCoeffs(su2double val_CV1, su2double val_CV2, su2activematrix& interp_mat_inv,
                                   std::array<su2double, 3>& interp_coeffs) {
  std::array<su2double, 3> query_vector = {1, val_CV1, val_CV2};

  su2double d;
  for (int i = 0; i < 3; i++) {
    d = 0;
    for (int j = 0; j < 3; j++) {
      d = d + interp_mat_inv[i][j] * query_vector[j];
    }
    interp_coeffs[i] = d;
  }
}

su2double CLookUpTable::Interpolate(const su2double* val_samples, std::array<unsigned long, 3>& val_triangle,
                                    std::array<su2double, 3>& val_interp_coeffs) {
  su2double result = 0;
  su2double z = 0;

  for (int i_point = 0; i_point < N_POINTS_TRIANGLE; i_point++) {
    z = val_samples[val_triangle[i_point]];
    result += val_interp_coeffs[i_point] * z;
  }

  return result;
}

unsigned long CLookUpTable::FindNearestNeighborOnHull(su2double val_CV1, su2double val_CV2, unsigned long i_level) {
  su2double min_distance = 1.e99;
  su2double next_distance = 1.e99;
  su2double next_x_norm;
  su2double next_y_norm;
  unsigned long neighbor_id = 0;

  const su2double* x_table = GetDataP(name_CV1, i_level);
  const su2double* y_table = GetDataP(name_CV2, i_level);

  su2double norm_coeff_x = 1. / (limits_table_x[i_level].second - limits_table_x[i_level].first);
  su2double norm_coeff_y = 1. / (limits_table_y[i_level].second - limits_table_y[i_level].first);
  su2double val_CV1_norm = val_CV1 / (limits_table_x[i_level].second - limits_table_x[i_level].first);
  su2double val_CV2_norm = val_CV2 / (limits_table_y[i_level].second - limits_table_y[i_level].first);

  for (unsigned long i_point = 0; i_point < n_hull_points[i_level]; ++i_point) {
    next_x_norm = x_table[hull[i_level][i_point]] * norm_coeff_x;
    next_y_norm = y_table[hull[i_level][i_point]] * norm_coeff_y;

    next_distance = sqrt(pow(val_CV1_norm - next_x_norm, 2) + pow(val_CV2_norm - next_y_norm, 2));

    if (next_distance < min_distance) {
      min_distance = next_distance;
      neighbor_id = hull[i_level][i_point];
    }
  }
  return neighbor_id;
}

void CLookUpTable::InterpolateToNearestNeighbors(const su2double val_CV1, const su2double val_CV2,
                                                 const std::vector<std::string>& names_var,
                                                 std::vector<su2double*>& var_vals, const unsigned long i_level) {
  /* Interpolate data using distance-weighted averaging on the two nearest table nodes. */

  su2double min_distance = 1e99, second_distance = 1e99;
  su2double norm_coeff_x = 1. / (*limits_table_x[i_level].second - *limits_table_x[i_level].first);
  su2double norm_coeff_y = 1. / (*limits_table_y[i_level].second - *limits_table_y[i_level].first);
  su2double val_CV1_norm = val_CV1 / (*limits_table_x[i_level].second - *limits_table_x[i_level].first);
  su2double val_CV2_norm = val_CV2 / (*limits_table_y[i_level].second - *limits_table_y[i_level].first);

  const su2double* x_table = GetDataP(name_CV1, i_level);
  const su2double* y_table = GetDataP(name_CV2, i_level);
  unsigned long i_nearest = 0, i_second_nearest = 0;

  for (unsigned long i_point = 0; i_point < n_hull_points[i_level]; ++i_point) {
    su2double next_x_norm = x_table[hull[i_level][i_point]] * norm_coeff_x;
    su2double next_y_norm = y_table[hull[i_level][i_point]] * norm_coeff_y;

    su2double next_distance = pow(val_CV1_norm - next_x_norm, 2) + pow(val_CV2_norm - next_y_norm, 2);

    /* Check if distance to node is lower than current nearest or second nearest node. */
    if (next_distance < min_distance) {
      second_distance = min_distance;
      min_distance = next_distance;

      i_second_nearest = i_nearest;
      i_nearest = hull[i_level][i_point];
    } else if ((next_distance > min_distance) && (next_distance < second_distance)) {
      i_second_nearest = hull[i_level][i_point];

      second_distance = next_distance;
    }
  }

  min_distance = sqrt(min_distance);
  second_distance = sqrt(second_distance);

  /* Interpolate data using distance-weighted averaging */
  su2double delimiter = (1.0 / min_distance) + (1.0 / second_distance);
  for (auto iVar = 0u; iVar < var_vals.size(); iVar++) {
    su2double data_nearest = GetDataP(names_var[iVar], i_level)[i_nearest],
              data_second_nearest = GetDataP(names_var[iVar], i_level)[i_second_nearest];
    *var_vals[iVar] = (data_nearest * (1.0 / min_distance) + data_second_nearest * (1.0 / second_distance)) / delimiter;
  }
}

void CLookUpTable::InterpolateToNearestNeighbors(const su2double val_CV1, const su2double val_CV2,
                                                 const std::string& name_var, su2double* var_vals,
                                                 const unsigned long i_level) {
  std::vector<std::string> names_var = {name_var};
  std::vector<su2double*> val_names_var = {var_vals};
  InterpolateToNearestNeighbors(val_CV1, val_CV2, names_var, val_names_var, i_level);
}

bool CLookUpTable::IsInTriangle(su2double val_CV1, su2double val_CV2, unsigned long val_id_triangle,
                                unsigned long i_level) {
  su2double tri_x_0 = GetDataP(name_CV1, i_level)[triangles[i_level][val_id_triangle][0]];
  su2double tri_y_0 = GetDataP(name_CV2, i_level)[triangles[i_level][val_id_triangle][0]];

  su2double tri_x_1 = GetDataP(name_CV1, i_level)[triangles[i_level][val_id_triangle][1]];
  su2double tri_y_1 = GetDataP(name_CV2, i_level)[triangles[i_level][val_id_triangle][1]];

  su2double tri_x_2 = GetDataP(name_CV1, i_level)[triangles[i_level][val_id_triangle][2]];
  su2double tri_y_2 = GetDataP(name_CV2, i_level)[triangles[i_level][val_id_triangle][2]];

  su2double area_tri = TriArea(tri_x_0, tri_y_0, tri_x_1, tri_y_1, tri_x_2, tri_y_2);

  su2double area_0 = TriArea(val_CV1, val_CV2, tri_x_1, tri_y_1, tri_x_2, tri_y_2);
  su2double area_1 = TriArea(tri_x_0, tri_y_0, val_CV1, val_CV2, tri_x_2, tri_y_2);
  su2double area_2 = TriArea(tri_x_0, tri_y_0, tri_x_1, tri_y_1, val_CV1, val_CV2);

  return (abs(area_tri - (area_0 + area_1 + area_2)) < area_tri * 1e-10);
}
