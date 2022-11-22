/*!
 * \file CLookupTable.cpp
 * \brief tabulation of fluid properties
 * \author D. Mayer, T. Economon, E.C.Bunschoten
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

#include "../../../Common/include/containers/CLookUpTable.hpp"
#include "../../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../../Common/include/toolboxes/CSquareMatrixCM.hpp"

using namespace std;

CLookUpTable::CLookUpTable(const string& var_file_name_lut, const string& name_x, const string& name_y) {
  file_name_lut = var_file_name_lut;

  rank = SU2_MPI::GetRank();

  name_CV1 = name_x;
  name_CV2 = name_y;

  LoadTableRaw(var_file_name_lut);

  FindTableLimits(name_x, name_y);

  if (rank == MASTER_NODE)
    cout << "Detecting all unique edges and setting edge to triangle connectivity "
            "..."
         << endl;

  IdentifyUniqueEdges();

  if (rank == MASTER_NODE) 
    cout << " done." << endl;

  PrintTableInfo();

  if (rank == MASTER_NODE)
    cout << "Building a trapezoidal map for the (" + name_x + ", " + name_y + ") "
            "space ..."
         << endl;
  trap_map_x_y = new CTrapezoidalMap[n_table_levels];
  for(auto i_level = 0ul; i_level<n_table_levels; i_level++)
    trap_map_x_y[i_level] = CTrapezoidalMap(GetDataP(name_x, i_level), GetDataP(name_y, i_level), table_data[i_level].cols(), edges[i_level], edge_to_triangle[i_level]);

  if (rank == MASTER_NODE) 
    cout << " done." << endl;

  if (rank == MASTER_NODE) 
    cout << "Precomputing interpolation coefficients..." << endl;
  

  ComputeInterpCoeffs();

  if (rank == MASTER_NODE) 
    cout << "LUT fluid model ready for use" << endl;
}

CLookUpTable::~CLookUpTable()
{
  delete [] trap_map_x_y;
  delete [] n_points;
  delete [] n_triangles;
  delete [] n_hull_points;
  delete [] table_data;
  delete [] hull;
  delete [] triangles;
  delete [] interp_mat_inv_x_y;
  delete [] edges;
  delete [] edge_to_triangle;
}

void CLookUpTable::LoadTableRaw(const string& var_file_name_lut) {
  CFileReaderLUT file_reader;

  if (rank == MASTER_NODE) 
    cout << "Loading lookup table, filename = " << var_file_name_lut << " ..." << endl;

  file_reader.ReadRawLUT(var_file_name_lut);
  unsigned short table_dim = file_reader.GetTableDim();
  n_table_levels = file_reader.GetNLevels();

  n_points            = new unsigned long[n_table_levels];
  n_triangles         = new unsigned long[n_table_levels];
  n_hull_points       = new unsigned long[n_table_levels];
  table_data          = new su2activematrix[n_table_levels];
  hull                = new std::vector<unsigned long>[n_table_levels];
  triangles           = new su2matrix<unsigned long>[n_table_levels];
  interp_mat_inv_x_y  = new std::vector<su2activematrix>[n_table_levels];
  edges               = new std::vector<std::vector<unsigned long> >[n_table_levels];
  edge_to_triangle    = new std::vector<std::vector<unsigned long> >[n_table_levels];

  for(unsigned long i_level=0; i_level<n_table_levels; i_level++){
    n_points[i_level] = file_reader.GetNPoints(i_level);
    n_triangles[i_level] = file_reader.GetNTriangles(i_level);
    n_hull_points[i_level] = file_reader.GetNHullPoints(i_level);
    table_data[i_level] = file_reader.GetTableData(i_level);
    triangles[i_level] = file_reader.GetTriangles(i_level);
    hull[i_level] = file_reader.GetHull(i_level);
  }
  
  n_variables = file_reader.GetNVariables();
  version_lut = file_reader.GetVersionLUT();
  version_reader = file_reader.GetVersionReader();
  names_var = file_reader.GetNamesVar();

  if(table_dim == 3){
    z_values_levels = new su2double[n_table_levels];
    for(unsigned long i_level=0; i_level<n_table_levels; i_level++){
      z_values_levels[i_level] = file_reader.GetTableLevel(i_level);
    }
  }
  if (rank == MASTER_NODE) 
    cout << " done." << endl;
}

void CLookUpTable::FindTableLimits(const string& name_x, const string& name_y) {
  int idx_X = GetIndexOfVar(name_x);
  int idx_Y = GetIndexOfVar(name_y);
  limits_table_x = new pair<su2double, su2double>[n_table_levels];
  limits_table_y = new pair<su2double, su2double>[n_table_levels];

  /* we find the lowest and highest value of y and x in the table */
  for(auto i_level=0u; i_level<n_table_levels; i_level++)
  {
    su2double min_x, max_x, min_y, max_y;
    min_y = *min_element(&table_data[i_level][idx_Y][0], &table_data[i_level][idx_Y][0]+table_data[i_level].cols());
    min_x = *min_element(&table_data[i_level][idx_X][0], &table_data[i_level][idx_X][0]+table_data[i_level].cols());
    max_y = *max_element(&table_data[i_level][idx_Y][0], &table_data[i_level][idx_Y][0]+table_data[i_level].cols());
    max_x = *max_element(&table_data[i_level][idx_X][0], &table_data[i_level][idx_X][0]+table_data[i_level].cols());
    limits_table_x[i_level] = make_pair(min_x, max_x);
    limits_table_y[i_level] = make_pair(min_y, max_y);
  }
}

void CLookUpTable::PrintTableInfo() {
  //TODO: make compatible with 3D table
  if (rank == MASTER_NODE) {
    cout << setfill(' ');
    cout << endl;
    cout << "+------------------------------------------------------------------+\n";
    cout << "|                     Look-Up-Table (LUT) info                     |\n";
    cout << "+------------------------------------------------------------------+" << endl;
    cout << "| File name:" << setw(54) << right << file_name_lut << " |" << endl;
    cout << "| Table version:" << setw(50) << right << version_lut << " |" << endl;
    cout << "| Table reader version:" << setw(43) << right << version_reader << " |" << endl;
    cout << "| Number of variables:" << setw(44) << right << n_variables << " |" << endl;
    cout << "| Number of points:" << setw(47) << right << n_points[0] << " |" << endl;
    cout << "| Number of triangles:" << setw(44) << right << n_triangles[0] << " |" << endl;
    cout << "| Number of edges:" << setw(48) << right << edges[0].size() << " |" << endl;
    cout << "+------------------------------------------------------------------+" << endl;
    cout << "| Minimum Y:" << setw(47) << right << limits_table_y[0].first << " |" << endl;
    cout << "| Maximum Y:" << setw(47) << right << limits_table_y[0].second << " |" << endl;
    cout << "| Minimum X:" << setw(38) << right << limits_table_x[0].first << " |" << endl;
    cout << "| Maximum X:" << setw(38) << right << limits_table_x[0].second << " |" << endl;
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

  for(auto i_level=0u; i_level<n_table_levels; i_level++){

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
        vector<unsigned long> edge(2);
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
  for(auto i_level=0ul; i_level<n_table_levels; i_level++){
    /* build KD tree for y, x space */

    std::array<unsigned long, 3> next_triangle;

    const su2double* val_x = GetDataP(name_CV1, i_level);
    const su2double* val_y = GetDataP(name_CV2, i_level);

    /* calculate weights for each triangle (basically a distance function) and
    * build inverse interpolation matrices */
    for (unsigned long i_triangle = 0; i_triangle < n_triangles[i_level]; i_triangle++) {

      for (int p = 0; p < 3; p++) {
        next_triangle[p] = triangles[i_level][i_triangle][p];
      }

      su2activematrix x_interp_mat_inv(3, 3);
      GetInterpMatInv(val_x, val_y, next_triangle, x_interp_mat_inv);
      interp_mat_inv_x_y[i_level].push_back(x_interp_mat_inv);
    }
  }
}

void CLookUpTable::GetInterpMatInv(const su2double* vec_x, const su2double* vec_y,
                                   std::array<unsigned long,3>& point_ids, su2activematrix& interp_mat_inv) {
  const unsigned int M = 3;
  CSquareMatrixCM global_M(3);

  /* setup LHM matrix for the interpolation */
  for (unsigned int i_point = 0; i_point < M; i_point++) {
    su2double x = vec_x[point_ids[i_point]];
    su2double y = vec_y[point_ids[i_point]];

    global_M(i_point,0) = SU2_TYPE::GetValue(1.0);
    global_M(i_point,1) = SU2_TYPE::GetValue(x);
    global_M(i_point,2) = SU2_TYPE::GetValue(y);
  }

  global_M.Invert();
  global_M.Transpose();
 
  for (unsigned int i=0; i<M; i++){
    for (unsigned int j=0; j<M; j++){
      interp_mat_inv[i][j] = global_M(i,j);
    }
  }

}

std::pair<unsigned long, unsigned long> CLookUpTable::FindInclusionLevels(su2double val_z, bool*within_table_limits)
{
  std::pair<unsigned long, unsigned long> inclusion_levels;
  if(val_z >= limits_table_z.second){
    *within_table_limits = false;
    inclusion_levels = std::make_pair(n_table_levels-1, n_table_levels-1);
  }else if(val_z <= limits_table_z.first){
    *within_table_limits = false;
    inclusion_levels = std::make_pair(0, 0);
  }else{
    for(auto i_level = 0ul; i_level < n_table_levels-1; i_level++){
      if((val_z >= z_values_levels[i_level]) && (val_z < z_values_levels[i_level+1])){
        *within_table_limits = true;
        inclusion_levels = std::make_pair(i_level, i_level+1);
      }
    }
  }
  return inclusion_levels;
}

unsigned long CLookUpTable::LookUp_XYZ(const std::string& val_name_var, su2double* val_var, su2double val_x,
                                su2double val_y, su2double val_z)
{
  bool within_z_limits, exit_code;
  unsigned long lower_level, upper_level;
  std::pair<unsigned long, unsigned long> inclusion_levels;
  if(val_z == 0){
    inclusion_levels = std::make_pair(0, 0);
    within_z_limits = false;
  }else{
    inclusion_levels = FindInclusionLevels(val_z, &within_z_limits);
  }
  
  if(within_z_limits){
    lower_level = inclusion_levels.first;
    upper_level = inclusion_levels.second;
    su2double val_var_lower, val_var_upper;
    unsigned long exit_code_lower = LookUp_XY(val_name_var, &val_var_lower, val_x, val_y, lower_level);
    unsigned long exit_code_upper = LookUp_XY(val_name_var, &val_var_upper, val_x, val_y, upper_level);

    Linear_Interpolation(val_z, inclusion_levels, val_var_lower, val_var_upper, val_var);

    return max(exit_code_lower, exit_code_upper);
  }else{
    lower_level = inclusion_levels.first;
    exit_code = LookUp_XY(val_name_var, val_var, val_x, val_y, lower_level);
    if( val_z == 0 ){
      return exit_code;
    }else{
      return 1;
    }
  }
}
unsigned long CLookUpTable::LookUp_XYZ(const std::vector<std::string>& val_names_var, std::vector<su2double*>& val_vars, su2double val_x,
                                su2double val_y, su2double val_z)
{
  bool within_z_limits, exit_code;
  unsigned long lower_level, upper_level;
  std::pair<unsigned long, unsigned long> inclusion_levels;
  if(val_z == 0){
    inclusion_levels = std::make_pair(0, 0);
    within_z_limits = false;
  }else{
    inclusion_levels = FindInclusionLevels(val_z, &within_z_limits);
  }
  
  if(within_z_limits){
    lower_level = inclusion_levels.first;
    upper_level = inclusion_levels.second;
    std::vector<su2double> val_vars_lower, val_vars_upper;
    val_vars_lower.resize(val_vars.size());
    val_vars_upper.resize(val_vars.size());
    unsigned long exit_code_lower = LookUp_XY(val_names_var, val_vars_lower, val_x, val_y, lower_level);
    unsigned long exit_code_upper = LookUp_XY(val_names_var, val_vars_upper, val_x, val_y, upper_level);

    Linear_Interpolation(val_z, inclusion_levels, val_vars_lower, val_vars_upper, val_vars);

    return max(exit_code_lower, exit_code_upper);
  }else{
    lower_level = inclusion_levels.first;
    exit_code = LookUp_XY(val_names_var, val_vars, val_x, val_y, lower_level);
    if( val_z == 0 ){
      return exit_code;
    }else{
      return 1;
    }
  }
}

void CLookUpTable::Linear_Interpolation(su2double val_z, std::pair<unsigned long, unsigned long> &inclusion_levels, su2double&lower_value,su2double&upper_value, su2double*&var_val)
{
  su2double val_z_lower = z_values_levels[inclusion_levels.first];
  su2double val_z_upper = z_values_levels[inclusion_levels.second];
  su2double factor_upper = (val_z - val_z_lower) / (val_z_upper - val_z_lower);
  su2double factor_lower = (val_z_upper - val_z) / (val_z_upper - val_z_lower);

  *var_val = lower_value * factor_lower + upper_value * factor_upper;
  
}

void CLookUpTable::Linear_Interpolation(su2double val_z, std::pair<unsigned long, unsigned long> &inclusion_levels, std::vector<su2double>&lower_values,std::vector<su2double>&upper_values, std::vector<su2double*>&var_vals)
{
  su2double val_z_lower = z_values_levels[inclusion_levels.first];
  su2double val_z_upper = z_values_levels[inclusion_levels.second];
  su2double factor_upper = (val_z - val_z_lower) / (val_z_upper - val_z_lower);
  su2double factor_lower = (val_z_upper - val_z) / (val_z_upper - val_z_lower);

  for(unsigned short iVar=0; iVar<var_vals.size(); iVar++){
    *var_vals[iVar] = lower_values[iVar] * factor_lower + upper_values[iVar] * factor_upper;
  }
}

unsigned long CLookUpTable::LookUp_XY(const string& val_name_var, su2double *val_var, su2double val_x,
                                            su2double val_y, unsigned long i_level) {
  unsigned long exit_code = 0;

  if (val_name_var.compare("NULL") == 0) {
    *val_var = 0.0;
    exit_code = 0;
    return exit_code;
  }

  /* check if x and y value is in table range */
  if ((val_x >= limits_table_x[i_level].first && val_x <= limits_table_x[i_level].second) &&
      (val_y >= limits_table_y[i_level].first && val_y <= limits_table_y[i_level].second)){

    /* find the triangle that holds the (x, y) point */
    unsigned long id_triangle = trap_map_x_y[i_level].GetTriangle(val_x, val_y);

    if (IsInTriangle(val_x, val_y, id_triangle, i_level)) {
      /* get interpolation coefficients for point on triangle */
      std::array<su2double,3> interp_coeffs{0};
      GetInterpCoeffs(val_x, val_y, interp_mat_inv_x_y[i_level][id_triangle], interp_coeffs);

      /* first, copy the single triangle from the large triangle list*/
      std::array<unsigned long,3> triangle{0};
      for (int p = 0; p < 3; p++) 
        triangle[p] = triangles[i_level][id_triangle][p]; 

      *val_var = Interpolate(GetDataP(val_name_var), triangle, interp_coeffs);
      exit_code = 0;
    } else {
      /* in bounding box but outside of table */
      unsigned long nearest_neighbor = FindNearestNeighborOnHull(val_x, val_y, i_level);
      *val_var = GetDataP(val_name_var, i_level)[nearest_neighbor];
      exit_code = 1;
    }
  } else {
    /* lookup is outside of table bounding box */
    unsigned long nearest_neighbor = FindNearestNeighborOnHull(val_x, val_y, i_level);
    *val_var = GetDataP(val_name_var, i_level)[nearest_neighbor];
    exit_code = 1;
  }
  return exit_code;
}

unsigned long CLookUpTable::LookUp_XY(const vector<string>& val_names_var, vector<su2double>& val_vars,
                                            su2double val_x, su2double val_y, unsigned long i_level) {
  vector<su2double*> look_up_data;

  for (long unsigned int i_var = 0; i_var < val_vars.size(); ++i_var) {
    look_up_data.push_back(&val_vars[i_var]);
  }

  unsigned long exit_code = LookUp_XY(val_names_var, look_up_data, val_x, val_y, i_level);

  return exit_code;
}

unsigned long CLookUpTable::LookUp_XY(const vector<string>& val_names_var, vector<su2double*>& val_vars,
                                            su2double val_x, su2double val_y, unsigned long i_level) {
  unsigned long exit_code = 0;
  unsigned long nearest_neighbor = 0;
  unsigned long id_triangle;
  std::array<su2double,3> interp_coeffs{0};
  std::array<unsigned long,3> triangle{0};

  /* check if x value is in table x-dimension range
   * and if y is in table y-dimension table range */
  if ((val_x >= limits_table_x[i_level].first && val_x <= limits_table_x[i_level].second) &&
      (val_y >= limits_table_y[i_level].first && val_y <= limits_table_y[i_level].second)){
    /* if so, try to find the triangle that holds the (prog, enth) point */
    id_triangle = trap_map_x_y[i_level].GetTriangle(val_x, val_y);

    /* check if point is inside a triangle (if table domain is non-rectangular,
     * the previous range check might be true but the point could still be outside of the domain) */
    if (IsInTriangle(val_x, val_y, id_triangle, i_level)) {
      /* if so, get interpolation coefficients for point in  the triangle */
      GetInterpCoeffs(val_x, val_y, interp_mat_inv_x_y[i_level][id_triangle], interp_coeffs);

      /* exit_code 0 means point was in triangle */
      exit_code = 0;

    } else {
      /* if point is not inside a triangle (outside table domain) search nearest neighbor */
      nearest_neighbor = FindNearestNeighborOnHull(val_x, val_y, i_level);
      exit_code = 1;
    }

  } else {
    /* if point is outside of table ranges, find nearest neighbor */
    nearest_neighbor = FindNearestNeighborOnHull(val_x, val_y, i_level);
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
          triangle[p] = triangles[i_level][id_triangle][p]; 
        *val_vars[i_var] = Interpolate(GetDataP(val_names_var[i_var]), triangle, interp_coeffs);
      }
      else
        *val_vars[i_var] = GetDataP(val_names_var[i_var], i_level)[nearest_neighbor];
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

unsigned long CLookUpTable::FindNearestNeighborOnHull(su2double val_x, su2double val_y, unsigned long i_level) {
  su2double min_distance = 1.e99;
  su2double next_distance = 1.e99;
  su2double next_x_norm;
  su2double next_y_norm;
  unsigned long neighbor_id = 0;

  const su2double* x_table = GetDataP(name_CV1, i_level);
  const su2double* y_table = GetDataP(name_CV2, i_level);

  su2double norm_coeff_x = 1. / (limits_table_x[i_level].second - limits_table_x[i_level].first);
  su2double norm_coeff_y = 1. / (limits_table_y[i_level].second - limits_table_y[i_level].first);
  su2double val_x_norm = val_x / (limits_table_x[i_level].second - limits_table_x[i_level].first);
  su2double val_y_norm = val_y / (limits_table_y[i_level].second - limits_table_y[i_level].first);

  for (unsigned long i_point = 0; i_point < n_hull_points[i_level]; ++i_point) {
    next_x_norm = x_table[hull[i_level][i_point]] * norm_coeff_x;
    next_y_norm = y_table[hull[i_level][i_point]] * norm_coeff_y;

    next_distance = sqrt(pow(val_x_norm - next_x_norm, 2) + pow(val_y_norm - next_y_norm, 2));

    if (next_distance < min_distance) {
      min_distance = next_distance;
      neighbor_id = hull[i_level][i_point];
    }
  }
  return neighbor_id;
}

bool CLookUpTable::IsInTriangle(su2double val_x, su2double val_y, unsigned long val_id_triangle, unsigned long i_level) {
  su2double tri_x_0 = GetDataP(name_CV1, i_level)[triangles[i_level][val_id_triangle][0]];
  su2double tri_y_0 = GetDataP(name_CV2, i_level)[triangles[i_level][val_id_triangle][0]];

  su2double tri_x_1 = GetDataP(name_CV1, i_level)[triangles[i_level][val_id_triangle][1]];
  su2double tri_y_1 = GetDataP(name_CV2, i_level)[triangles[i_level][val_id_triangle][1]];
  
  su2double tri_x_2 = GetDataP(name_CV1, i_level)[triangles[i_level][val_id_triangle][2]];
  su2double tri_y_2 = GetDataP(name_CV2, i_level)[triangles[i_level][val_id_triangle][2]];

  su2double area_tri = TriArea(tri_x_0, tri_y_0, tri_x_1, tri_y_1, tri_x_2, tri_y_2);

  su2double area_0 = TriArea(val_x, val_y, tri_x_1, tri_y_1, tri_x_2, tri_y_2);
  su2double area_1 = TriArea(tri_x_0, tri_y_0, val_x, val_y, tri_x_2, tri_y_2);
  su2double area_2 = TriArea(tri_x_0, tri_y_0, tri_x_1, tri_y_1, val_x, val_y);

  return (abs(area_tri - (area_0 + area_1 + area_2)) < area_tri * 1e-10);
}
