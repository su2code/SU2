/*!
 * \file CLookupTable.cpp
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

#include "../../../Common/include/containers/CLookUpTable.hpp"
#include "../../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../../Common/include/toolboxes/CSquareMatrixCM.hpp"

using namespace std;

CLookUpTable::CLookUpTable(const string& var_file_name_lut, const string& name_prog, const string& name_enth) {
  file_name_lut = var_file_name_lut;

  rank = SU2_MPI::GetRank();

  LoadTableRaw(var_file_name_lut);

  FindTableLimits(name_prog, name_enth);

  if (rank == MASTER_NODE)
    cout << "Detecting all unique edges and setting edge to triangle connectivity "
            "..."
         << endl;

  IdentifyUniqueEdges();

  if (rank == MASTER_NODE) cout << " done." << endl;

  PrintTableInfo();

  if (rank == MASTER_NODE)
    cout << "Building a trapezoidal map for the (progress variable, enthalpy) "
            "space ..."
         << endl;

  trap_map_prog_enth = CTrapezoidalMap(GetDataP(name_prog), GetDataP(name_enth), table_data.cols(), edges, edge_to_triangle);

  if (rank == MASTER_NODE) cout << " done." << endl;

  if (rank == MASTER_NODE) {
    cout << "Precomputing interpolation coefficients..." << endl;
  }

  ComputeInterpCoeffs(name_prog, name_enth);

  if (rank == MASTER_NODE) cout << "LUT fluid model ready for use" << endl;
}

void CLookUpTable::LoadTableRaw(const string& var_file_name_lut) {
  CFileReaderLUT file_reader;

  if (rank == MASTER_NODE) cout << "Loading look-up-table-file " << var_file_name_lut << " ..." << endl;

  file_reader.ReadRawDRG(var_file_name_lut);

  n_points = file_reader.GetNPoints();
  n_triangles = file_reader.GetNTriangles();
  n_variables = file_reader.GetNVariables();
  n_hull_points = file_reader.GetNHullPoints();
  type_lut = file_reader.GetTypeLUT();
  version_lut = file_reader.GetVersionLUT();
  version_reader = file_reader.GetVersionReader();
  names_var = file_reader.GetNamesVar();
  table_data = file_reader.GetTableData();
  triangles = file_reader.GetTriangles();
  hull = file_reader.GetHull();

  if (rank == MASTER_NODE) cout << " done." << endl;
}

void CLookUpTable::FindTableLimits(const string& name_prog, const string& name_enth) {
  int ixEnth = GetIndexOfVar(name_enth);
  int ixProg = GetIndexOfVar(name_prog);

  /* we find the lowest and highest value of enthalpy and progress variable in the table */
  limits_table_enth[0] = *min_element(&table_data[ixEnth][0], &table_data[ixEnth][0]+table_data.cols());
  limits_table_enth[1] = *max_element(&table_data[ixEnth][0], &table_data[ixEnth][0]+table_data.cols());
  limits_table_prog[0] = *min_element(&table_data[ixProg][0], &table_data[ixProg][0]+table_data.cols());
  limits_table_prog[1] = *max_element(&table_data[ixProg][0], &table_data[ixProg][0]+table_data.cols());

}

void CLookUpTable::PrintTableInfo() {
  if (rank == MASTER_NODE) {
    cout << setfill(' ');
    cout << endl;
    cout << "+------------------------------------------------------------------+\n";
    cout << "|                     Look-Up-Table (LUT) info                     |\n";
    cout << "+------------------------------------------------------------------+" << endl;
    cout << "| File name:" << setw(54) << right << file_name_lut << " |" << endl;
    cout << "| Table type:" << setw(53) << right << type_lut << " |" << endl;
    cout << "| Table version:" << setw(50) << right << version_lut << " |" << endl;
    cout << "| Table reader version:" << setw(43) << right << version_reader << " |" << endl;
    cout << "| Number of variables:" << setw(44) << right << n_variables << " |" << endl;
    cout << "| Number of points:" << setw(47) << right << n_points << " |" << endl;
    cout << "| Number of triangles:" << setw(44) << right << n_triangles << " |" << endl;
    cout << "| Number of edges:" << setw(48) << right << edges.size() << " |" << endl;
    cout << "+------------------------------------------------------------------+" << endl;
    cout << "| Minimum enthalpy:" << setw(47) << right << limits_table_enth[0] << " |" << endl;
    cout << "| Maximum enthalpy:" << setw(47) << right << limits_table_enth[1] << " |" << endl;
    cout << "| Minimum progress variable:" << setw(38) << right << limits_table_prog[0] << " |" << endl;
    cout << "| Maximum progress variable:" << setw(38) << right << limits_table_prog[1] << " |" << endl;
    cout << "+------------------------------------------------------------------+" << endl;
    cout << "| Variable names:                                                  |" << endl;
    cout << "|                                                                  |" << endl;

    for (unsigned long i_var = 0; i_var < names_var.size(); i_var++)
      cout << "| " << right << setw(3) << i_var << ": " << left << setw(59) << names_var[i_var] << " |" << endl;

    cout << "+------------------------------------------------------------------+" << endl;

    cout << endl;
  }
}

void CLookUpTable::IdentifyUniqueEdges() {

  /* loop through elements and store the element ID as
     a neighbor for each point in the element */

  vector<vector<unsigned long> > neighborElemsOfPoint;
  
  neighborElemsOfPoint.resize(n_points);
  for (unsigned long iElem = 0; iElem < n_triangles; iElem++) {
    /* loop over 3 points per triangle */
    for (unsigned long iPoint = 0; iPoint < N_POINTS_TRIANGLE; iPoint++) {
      /* get the global ID of the current point */
      const unsigned long GlobalIndex = triangles[iElem][iPoint];

      /* add the current element ID to the neighbor list for this point */
      neighborElemsOfPoint[GlobalIndex].push_back(iElem);
    }
  }

  /* remove duplicates from the neighboring element lists*/
  vector<unsigned long>::iterator vecIt;
  for (unsigned long iPoint = 0; iPoint < n_points; iPoint++) {
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
  neighborPointsOfPoint.resize(n_points);
  for (unsigned long iPoint = 0; iPoint < n_points; iPoint++) {
    for (unsigned long iElem = 0; iElem < neighborElemsOfPoint[iPoint].size(); iElem++) {
      /* loop over element points */
      for (unsigned long jPoint = 0; jPoint < N_POINTS_TRIANGLE; jPoint++) {
        /* get the global ID of the current point */
        const unsigned long GlobalIndex = triangles[neighborElemsOfPoint[iPoint][iElem]][jPoint];

        /* add the current element ID to the neighbor list for this point */
        if (GlobalIndex != iPoint) neighborPointsOfPoint[iPoint].push_back(GlobalIndex);
      }
    }
  }

  /* remove duplicates from the neighboring points lists */
  for (unsigned long iPoint = 0; iPoint < n_points; iPoint++) {
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
  for (unsigned long iPoint = 0; iPoint < n_points; iPoint++) {
    for (unsigned long jPoint = 0; jPoint < neighborPointsOfPoint[iPoint].size(); jPoint++) {
      /* Store the neighbor index more clearly. */
      const unsigned long GlobalIndex = neighborPointsOfPoint[iPoint][jPoint];

      /* Store the edge so that the lower index of the pair is always first. */
      if (iPoint < GlobalIndex) {
        vector<unsigned long> edge(2);
        edge[0] = iPoint;
        edge[1] = GlobalIndex;
        edges.push_back(edge);
      }
    }
  }

  /* Loop over our edges data structure. For the first point in each
   pair, loop through the neighboring elements and store the two
   elements that contain the second point in the edge. */
  edge_to_triangle.resize(edges.size());
  for (unsigned long iEdge = 0; iEdge < edges.size(); iEdge++) {
    /* Store the two points of the edge more clearly. */
    const unsigned long iPoint = edges[iEdge][0];
    const unsigned long jPoint = edges[iEdge][1];

    /* Loop over all neighobring elements to iPoint. */
    for (unsigned long iElem = 0; iElem < neighborElemsOfPoint[iPoint].size(); iElem++) {
      /* loop over 3 points per triangle */
      for (unsigned long kPoint = 0; kPoint < N_POINTS_TRIANGLE; kPoint++) {
        /* Get the global ID of the current point. */
        const unsigned long GlobalIndex = triangles[neighborElemsOfPoint[iPoint][iElem]][kPoint];

        /* Add the current element ID to the neighbor list for this point. */
        if (GlobalIndex == jPoint) edge_to_triangle[iEdge].push_back(neighborElemsOfPoint[iPoint][iElem]);
      }
    }
  }
}

void CLookUpTable::ComputeInterpCoeffs(const string& name_prog, const string& name_enth) {
  /* build KD tree for enthalpy, progress variable space */

  std::array<unsigned long, 3> next_triangle;

  const su2double* prog = GetDataP(name_prog);
  const su2double* enth = GetDataP(name_enth);

  /* calculate weights for each triangle (basically a distance function) and
   * build inverse interpolation matrices */
  for (unsigned long i_triangle = 0; i_triangle < n_triangles; i_triangle++) {

    for (int p = 0; p < 3; p++) {
      next_triangle[p] = triangles[i_triangle][p];
    }

    su2activematrix prog_interp_mat_inv(3, 3);
    GetInterpMatInv(prog, enth, next_triangle, prog_interp_mat_inv);
    interp_mat_inv_prog_enth.push_back(prog_interp_mat_inv);
  }
}

void CLookUpTable::GetInterpMatInv(const su2double* vec_x, const su2double* vec_y,
                                   std::array<unsigned long,3>& point_ids, su2activematrix& interp_mat_inv) {
  unsigned int M = 3;
  CSquareMatrixCM global_M(M);

  /* setup LHM matrix for the interpolation */
  for (int i_point = 0; i_point < 3; i_point++) {
    su2double x = vec_x[point_ids[i_point]];
    su2double y = vec_y[point_ids[i_point]];

    global_M(i_point,0) = SU2_TYPE::GetValue(1.0);
    global_M(i_point,1) = SU2_TYPE::GetValue(x);
    global_M(i_point,2) = SU2_TYPE::GetValue(y);
  }

  global_M.Invert();
  global_M.Transpose();
 
  /* convert back into vector<vector>> for now */
  for (unsigned int i=0; i<M; i++){
    for (unsigned int j=0; j<M; j++){
      interp_mat_inv[i][j] = global_M(i,j);
    }
  }

}

unsigned long CLookUpTable::LookUp_ProgEnth(const string& val_name_var, su2double *val_var, su2double val_prog,
                                            su2double val_enth, const string& name_prog, const string& name_enth) {
  unsigned long exit_code = 0;

  if (val_name_var.compare("NULL") == 0) {
    *val_var = 0.0;
    exit_code = 0;
    return exit_code;
  }

  /* check if progress variable and enthalpy value is in table range */
  if (val_prog >= limits_table_prog[0] && val_prog <= limits_table_prog[1]) {
    //&&
    // val_enth >= limits_table_enth[0] && val_enth <= limits_table_enth[1] ){

    /* find the triangle that holds the (prog, enth) point */
    unsigned long id_triangle = trap_map_prog_enth.GetTriangle(val_prog, val_enth);

    if (IsInTriangle(val_prog, val_enth, id_triangle, name_prog, name_enth)) {
      /* get interpolation coefficients for point on triangle */
      std::array<su2double,3> interp_coeffs{0};
      GetInterpCoeffs(val_prog, val_enth, interp_mat_inv_prog_enth[id_triangle], interp_coeffs);

      /* first, copy the single triangle from the large triangle list*/
      std::array<unsigned long,3> triangle{0};
      for (int p = 0; p < 3; p++) 
        triangle[p] = triangles[id_triangle][p]; 

      //*val_var = Interpolate(GetData(val_name_var), triangles[id_triangle], interp_coeffs);
      *val_var = Interpolate(GetDataP(val_name_var), triangle, interp_coeffs);
      exit_code = 0;
    } else {
      /* in bounding box but outside of table */
      unsigned long nearest_neighbor = FindNearestNeighborOnHull(val_prog, val_enth, name_prog, name_enth);
      *val_var = GetDataP(val_name_var)[nearest_neighbor];
      exit_code = 1;
    }
  } else {
    /* lookup is outside of table bounding box */
    unsigned long nearest_neighbor = FindNearestNeighborOnHull(val_prog, val_enth, name_prog, name_enth);
    *val_var = GetDataP(val_name_var)[nearest_neighbor];
    exit_code = 1;
  }
  return exit_code;
}

unsigned long CLookUpTable::LookUp_ProgEnth(vector<string>& val_names_var, vector<su2double>& val_vars,
                                            su2double val_prog, su2double val_enth, const string& name_prog,
                                            const string& name_enth) {
  vector<su2double*> look_up_data;

  for (long unsigned int i_var = 0; i_var < val_vars.size(); ++i_var) {
    look_up_data.push_back(&val_vars[i_var]);
  }

  unsigned long exit_code = LookUp_ProgEnth(val_names_var, look_up_data, val_prog, val_enth, name_prog, name_enth);

  return exit_code;
}

unsigned long CLookUpTable::LookUp_ProgEnth(vector<string>& val_names_var, vector<su2double*>& val_vars,
                                            su2double val_prog, su2double val_enth, const string& name_prog,
                                            const string& name_enth) {
  unsigned long exit_code = 0;
  unsigned long nearest_neighbor = 0;
  unsigned long id_triangle;
  std::array<su2double,3> interp_coeffs{0};
  std::array<unsigned long,3> triangle{0};

  /* check if progress variable value is in progress variable table range
   * and if enthalpy is in enthalpy table range */
  if (val_prog >= limits_table_prog[0] && val_prog <= limits_table_prog[1] && val_enth >= limits_table_enth[0] &&
      val_enth <= limits_table_enth[1]) {
    /* if so, try to find the triangle that holds the (prog, enth) point */
    id_triangle = trap_map_prog_enth.GetTriangle(val_prog, val_enth);

    /* check if point is inside a triangle (if table domain is non-rectangular,
     * the previous range check might be true but the point is still outside of the domain) */
    if (IsInTriangle(val_prog, val_enth, id_triangle, name_prog, name_enth)) {
      /* if so, get interpolation coefficients for point in  the triangle */
      GetInterpCoeffs(val_prog, val_enth, interp_mat_inv_prog_enth[id_triangle], interp_coeffs);

      /* exit_code 0 means point was in triangle */
      exit_code = 0;

    } else {
      // cout << "lookup_progenth: outside table range, c,h = "<< val_prog<< " "<<val_enth<<endl;
      /* if point is not inside a triangle (outside table domain) search nearest neighbor */
      nearest_neighbor = FindNearestNeighborOnHull(val_prog, val_enth, name_prog, name_enth);
      exit_code = 1;
    }

  } else {
    // if (rank == MASTER_NODE) cout << "WARNING: LookUp_ProgEnth: lookup is outside of table bounding box, c,h = "<<

    /* if point is outside of table range, search nearest neighbor */
    nearest_neighbor = FindNearestNeighborOnHull(val_prog, val_enth, name_prog, name_enth);
    exit_code = 1;
  }

  /* loop over variable names and interpolate / get values */
  for (long unsigned int i_var = 0; i_var < val_names_var.size(); ++i_var) {
    if (val_names_var[i_var].compare("NULL") == 0) {
      *val_vars[i_var] = 0.0;
    } else {
      if (exit_code == 0){

        /* first, copy the single triangle from the large triangle list*/
        for (int p = 0; p < 3; p++) 
          triangle[p] = triangles[id_triangle][p]; 
        *val_vars[i_var] = Interpolate(GetDataP(val_names_var[i_var]), triangle, interp_coeffs);
      }
      else
        *val_vars[i_var] = GetDataP(val_names_var[i_var])[nearest_neighbor];
    }
  }
  return exit_code;
}

void CLookUpTable::GetInterpCoeffs(su2double val_x, su2double val_y, su2activematrix& interp_mat_inv,
                                   std::array<su2double,3>& interp_coeffs) {

  std::array<su2double,3> query_vector = {1,val_x,val_y};

  su2double d;
  for (int i = 0; i < 3; i++) {
    d = 0;
    for (int j = 0; j < 3; j++) {
      d = d + interp_mat_inv[i][j] * query_vector[j];
    }
    interp_coeffs[i] = d;
  }
}

su2double CLookUpTable::Interpolate(const su2double* val_samples, std::array<unsigned long,3>& val_triangle,
                                    std::array<su2double,3>& val_interp_coeffs) {
  su2double result = 0;
  su2double z = 0;

  for (int i_point = 0; i_point < N_POINTS_TRIANGLE; i_point++) {
    z = val_samples[val_triangle[i_point]];
    result += val_interp_coeffs[i_point] * z;
  }

  return result;
}

unsigned long CLookUpTable::FindNearestNeighborOnHull(su2double val_prog, su2double val_enth, const string& name_prog,
                                                      const string& name_enth) {
  su2double min_distance = 1.e99;
  su2double next_distance = 1.e99;
  su2double next_prog_norm;
  su2double next_enth_norm;
  unsigned long neighbor_id = 0;

  const su2double* prog_table = GetDataP(name_prog);
  const su2double* enth_table = GetDataP(name_enth);

  su2double norm_coeff_prog = 1. / (limits_table_prog[1] - limits_table_prog[0]);
  su2double norm_coeff_enth = 1. / (limits_table_enth[1] - limits_table_enth[0]);
  su2double val_prog_norm = val_prog / (limits_table_prog[1] - limits_table_prog[0]);
  su2double val_enth_norm = val_enth / (limits_table_enth[1] - limits_table_enth[0]);

  for (unsigned long i_point = 0; i_point < n_hull_points; ++i_point) {
    next_prog_norm = prog_table[hull[i_point]] * norm_coeff_prog;
    next_enth_norm = enth_table[hull[i_point]] * norm_coeff_enth;

    next_distance = sqrt(pow(val_prog_norm - next_prog_norm, 2) + pow(val_enth_norm - next_enth_norm, 2));

    if (next_distance < min_distance) {
      min_distance = next_distance;
      neighbor_id = hull[i_point];
    }
  }
  return neighbor_id;
}

bool CLookUpTable::IsInTriangle(su2double val_prog, su2double val_enth, unsigned long val_id_triangle, const string& name_prog,
                                const string& name_enth) {
  su2double tri_prog_0 = GetDataP(name_prog)[triangles[val_id_triangle][0]];
  su2double tri_enth_0 = GetDataP(name_enth)[triangles[val_id_triangle][0]];
  su2double tri_prog_1 = GetDataP(name_prog)[triangles[val_id_triangle][1]];
  su2double tri_enth_1 = GetDataP(name_enth)[triangles[val_id_triangle][1]];
  su2double tri_prog_2 = GetDataP(name_prog)[triangles[val_id_triangle][2]];
  su2double tri_enth_2 = GetDataP(name_enth)[triangles[val_id_triangle][2]];

  su2double area_tri = TriArea(tri_prog_0, tri_enth_0, tri_prog_1, tri_enth_1, tri_prog_2, tri_enth_2);
  su2double area_0 = TriArea(val_prog, val_enth, tri_prog_1, tri_enth_1, tri_prog_2, tri_enth_2);
  su2double area_1 = TriArea(tri_prog_0, tri_enth_0, val_prog, val_enth, tri_prog_2, tri_enth_2);
  su2double area_2 = TriArea(tri_prog_0, tri_enth_0, tri_prog_1, tri_enth_1, val_prog, val_enth);

  return (abs(area_tri - (area_0 + area_1 + area_2)) < area_tri * 1e-10);
}
