#include "../include/numerics/CLookUpTable.hpp"

CLookUpTable::CLookUpTable(string var_file_name_lut, string name_prog, string name_enth) {

  file_name_lut = var_file_name_lut;

  rank = SU2_MPI::GetRank();

  LoadTableRaw(var_file_name_lut);

  FindTableLimits(name_prog,name_enth);

  if (rank == MASTER_NODE)
    cout << "Detecting all unique edges and setting edge to triangle connectivity "
            "..." << endl;

  IdentifyUniqueEdges();

  if (rank == MASTER_NODE) cout << " done." << endl;

  PrintTableInfo();

  if (rank == MASTER_NODE)
    cout << "Building a trapezoidal map for the (progress variable, enthalpy) "
            "space ..." << endl;

  
  //cout << "expecting progress variable at table column " << GetIndexOfVar(name_prog)<<endl;
  //cout << "expecting enthalpy at table column " << GetIndexOfVar(name_enth)<<endl;

  trap_map_prog_enth = CTrapezoidalMap(GetData(name_prog), GetData(name_enth),
                                       GetEdges(), GetEdgeToTriangle());
  

  /* nijso todo writing the trapezoidal map*/
  //if (rank == MASTER_NODE)
  //  trap_map_prog_enth.writeVec();

  /*
  ofstream trapfile;
  string filename = "trapfile.dat";
  trapfile.open (filename.c_str(), ios::out|ios::binary|ios::app);
  */


  if (rank == MASTER_NODE) cout << " done." << endl;

  if (rank == MASTER_NODE) {
    cout << "Precomputing interpolation coefficients..." << endl;
  }

  ComputeInterpCoeffs(name_prog,name_enth);
  
  if (rank == MASTER_NODE) cout << " done." << endl;

  if (rank == MASTER_NODE) cout << "LUT fluid model ready for use" << endl;
  
}

void CLookUpTable::LoadTableRaw(string var_file_name_lut) {

  CFileReaderLUT file_reader;

  if (rank == MASTER_NODE)
  cout << "Loading look-up-table-file " << var_file_name_lut << " ..." << endl;

  file_reader.ReadRawDRG(var_file_name_lut);

  n_points       = file_reader.GetNPoints();
  n_triangles    = file_reader.GetNTriangles();
  n_variables    = file_reader.GetNVariables();
  n_hull_points  = file_reader.GetNHullPoints();
  type_lut       = file_reader.GetTypeLUT();
  version_lut    = file_reader.GetVersionLUT();
  version_reader = file_reader.GetVersionReader();

  names_var  = file_reader.GetNamesVar();
  table_data = file_reader.GetTableData();
  triangles  = file_reader.GetTriangles();
  hull       = file_reader.GetHull();

  if (rank == MASTER_NODE) cout << " done." << endl;
  
}

void CLookUpTable::FindTableLimits(string name_prog, string name_enth) {

  int ixEnth = GetIndexOfVar(name_enth);
  int ixProg = GetIndexOfVar(name_prog);
  //cout << "findtablelimits:: expecting progress variable " << name_prog<< " at table column " << ixProg<<endl;
  //cout << "findtablelimits:: expecting enthalpy "<<name_enth<<" at table column " << ixEnth<<endl;

  


  limits_table_enth[0] =
      *min_element(table_data.at(ixEnth).begin(), table_data.at(ixEnth).end());
  limits_table_enth[1] =
      *max_element(table_data.at(ixEnth).begin(), table_data.at(ixEnth).end());
  limits_table_prog[0] =
      *min_element(table_data.at(ixProg).begin(), table_data.at(ixProg).end());
  limits_table_prog[1] =
      *max_element(table_data.at(ixProg).begin(), table_data.at(ixProg).end());

}

void CLookUpTable::PrintTableInfo() {

  if (rank == MASTER_NODE) {
    cout << setfill(' ');
    cout << endl;
    cout << "+------------------------------------------------------------------+"      << endl;
    cout << "|                     Look-Up-Table (LUT) info                     |"      << endl;
    cout << "+------------------------------------------------------------------+"      << endl;
    cout << "| File name:"                 << setw(54) << right << file_name_lut        << " |" << endl;
    cout << "| Table type:"                << setw(53) << right << type_lut             << " |" << endl;
    cout << "| Table version:"             << setw(50) << right << version_lut          << " |" << endl;
    cout << "| Table reader version:"      << setw(43) << right << version_reader       << " |" << endl;
    cout << "| Number of variables:"       << setw(44) << right << n_variables          << " |" << endl;
    cout << "| Number of points:"          << setw(47) << right << n_points             << " |" << endl;
    cout << "| Number of triangles:"       << setw(44) << right << n_triangles          << " |" << endl;
    cout << "| Number of edges:"           << setw(48) << right << edges.size()         << " |" << endl;
    cout << "+------------------------------------------------------------------+"      << endl;
    cout << "| Minimum enthalpy:"          << setw(47) << right << limits_table_enth[0] << " |" << endl;
    cout << "| Maximum enthalpy:"          << setw(47) << right << limits_table_enth[1] << " |" << endl;
    cout << "| Minimum progress variable:" << setw(38) << right << limits_table_prog[0] << " |" << endl;
    cout << "| Maximum progress variable:" << setw(38) << right << limits_table_prog[1] << " |" << endl;
    cout << "+------------------------------------------------------------------+"      << endl;
    cout << "| Variable names:                                                  |"      << endl;
    cout << "|                                                                  |"      << endl;

    for (unsigned long i_var = 0; i_var < names_var.size(); i_var++)
      cout << "| " << right << setw(3) << i_var << ": " << left << setw(59) 
      << names_var.at(i_var) << " |" << endl;

    cout << "+------------------------------------------------------------------+" << endl;

    cout << endl;
  }
}

void CLookUpTable::IdentifyUniqueEdges() {
  
  /* we will fill these two data members with the point pair and adjacent
     elements, respectively, for each unique edge in the grid */
  edges.clear();
  edge_to_triangle.clear();
  
  /* loop through elements and store the element ID as
     a neighbor for each point in the element */
  vector<vector<unsigned long> > neighborElemsOfPoint;
  neighborElemsOfPoint.resize(n_points);
  for (unsigned long iElem = 0; iElem < n_triangles; iElem++) {
    
    /* loop over 3 points per triangle */
    for (unsigned long iPoint = 0; iPoint < N_POINTS_TRIANGLE; iPoint++) {
      
      /* get the global ID of the current point */
      const unsigned long GlobalIndex = triangles.at(iElem).at(iPoint);
      
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
        const unsigned long GlobalIndex = triangles.at(neighborElemsOfPoint[iPoint][iElem]).at(jPoint);
        
        /* add the current element ID to the neighbor list for this point */
        if (GlobalIndex != iPoint)
          neighborPointsOfPoint[iPoint].push_back(GlobalIndex);
        
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
        const unsigned long GlobalIndex = triangles.at(neighborElemsOfPoint[iPoint][iElem]).at(kPoint);
        
        /* Add the current element ID to the neighbor list for this point. */
        if (GlobalIndex == jPoint)
          edge_to_triangle[iEdge].push_back(neighborElemsOfPoint[iPoint][iElem]);
        
      }
    }
  }
  
}

void CLookUpTable::ComputeInterpCoeffs(string name_prog, string name_enth) {
  
  /* build KD tree for enthalpy, progress variable space */
  vector< unsigned long > points(n_points);
  vector< su2double > weights(2, 0);

  vector< su2double > prog_enth_pairs(2 * n_points);

  vector< unsigned long > next_triangle;

  const vector<su2double> &prog = GetData(name_prog);
  const vector<su2double> &enth = GetData(name_enth);

  vector< unsigned long > result_ids;
  vector< int > result_ranks;
  vector< su2double > best_dist;
  
  for (unsigned long i_point = 0; i_point < n_points; i_point++) {

    points.push_back(i_point);

    prog_enth_pairs.push_back(prog.at(i_point));
    prog_enth_pairs.push_back(enth.at(i_point));
  }

  /* calculate weights for each triangle (basically a distance function) and
   * build inverse interpolation matrices */
  for (unsigned long i_triangle = 0; i_triangle < n_triangles; i_triangle++) {

    next_triangle = triangles.at(i_triangle);

    /* the query point is the weighted average of the vertexes of the triangle */
    weights.at(0) = 0;
    weights.at(1) = 0;
    
    /* enthalpy */
    weights.at(0) += enth.at(next_triangle.at(0));
    weights.at(0) += enth.at(next_triangle.at(1));
    weights.at(0) += enth.at(next_triangle.at(2));
    weights.at(0) /= 3;

    /* progress variable */
    weights.at(1) += prog.at(next_triangle.at(0));
    weights.at(1) += prog.at(next_triangle.at(1));
    weights.at(1) += prog.at(next_triangle.at(2));
    weights.at(1) /= 3;

    interp_points.push_back(next_triangle);

    // Now use the nearest 16 points to construct an interpolation function
    // for each search pair option
    vector< vector< su2double > > prog_interp_mat_inv(3, vector< su2double >(3,0));
    GetInterpMatInv(prog, enth, next_triangle, prog_interp_mat_inv);
    interp_mat_inv_prog_enth.push_back(prog_interp_mat_inv);

  }
}

void CLookUpTable::GetInterpMatInv(const vector<su2double>       &vec_x,
                                   const vector<su2double>       &vec_y,
                                   vector<unsigned long>         &point_ids,
                                   vector< vector< su2double > > &interp_mat_inv) {

  vector<vector<su2double> > interp_mat(3, vector<su2double>(3, 0));

  /* setup LHM matrix for the interpolation */
  for (int i_point = 0; i_point < 3; i_point++) {

    su2double x = vec_x.at(point_ids.at(i_point));
    su2double y = vec_y.at(point_ids.at(i_point));

    interp_mat.at(i_point).at(0) = 1;
    interp_mat.at(i_point).at(1) = x;
    interp_mat.at(i_point).at(2) = y;
  }

  /* invert the Interpolation matrix using Gaussian elimination with pivoting */
  GaussianInverse(interp_mat, interp_mat_inv);

  /* transpose the inverse */
  su2double swap_helper;
  for (int i = 0; i < 2; i++) {
    for (int j = i + 1; j < 3; j++) {
      swap_helper                = interp_mat_inv.at(i).at(j);
      interp_mat_inv.at(i).at(j) = interp_mat_inv.at(j).at(i);
      interp_mat_inv.at(j).at(i) = swap_helper;
    }
  }
  
}

void CLookUpTable::GaussianInverse(vector< vector< su2double > > &mat,
                                   vector< vector< su2double > > &mat_inv) {

  /* temp provides memory to invert mat */
  vector< vector< su2double > > temp;

  /* number of dimensions of mat */
  int n_dim = (int)mat.size();

  temp.resize(n_dim, vector< su2double >(2 * n_dim, 0));

  /* copy original matrix into inverse */
  for (int i = 0; i < n_dim; i++) {
    for (int j = 0; j < n_dim; j++) {
      temp.at(i).at(j)         = mat.at(i).at(j);
      temp.at(i).at(n_dim + j) = 0;
    }
    temp.at(i).at(n_dim + i) = 1;
  }

  su2double max_val;
  int max_idx;
  /* pivot each column such that the largest number possible divides the other
   * rows The goal is to avoid zeros or small numbers in division. */
  for (int k = 0; k < n_dim - 1; k++) {
    max_idx = k;
    max_val = abs(temp.at(k).at(k));

    /* find largest value (pivot) in column */
    for (int j = k; j < n_dim; j++) {
      if (abs(temp.at(j).at(k)) > max_val) {
        max_idx = j;
        max_val = abs(temp.at(j).at(k));
      }
    }

    /* move row with the highest value up */
    for (int j = 0; j < (n_dim * 2); j++) {
      su2double d = temp.at(k).at(j);
      temp.at(k).at(j)       = temp.at(max_idx).at(j);
      temp.at(max_idx).at(j) = d;
    }

    /* subtract the moved row from all other rows */
    for (int i = k + 1; i < n_dim; i++) {
      su2double c = temp.at(i).at(k) / temp.at(k).at(k);
      for (int j = 0; j < (n_dim * 2); j++) {
        temp.at(i).at(j) = temp.at(i).at(j) - temp.at(k).at(j) * c;
      }
    }
  }

  /* perform back-substitution */
  for (int k = n_dim - 1; k > 0; k--) {
    if (temp.at(k).at(k) != 0) {
      for (int i = k - 1; i > -1; i--) {
        su2double c = temp.at(i).at(k) / temp.at(k).at(k);
        for (int j = 0; j < (n_dim * 2); j++) {
          temp.at(i).at(j) = temp.at(i).at(j) - temp.at(k).at(j) * c;
        }
      }
    }
  }

  /* normalize the inverse */
  for (int i = 0; i < n_dim; i++) {
    su2double c = temp.at(i).at(i);
    for (int j = 0; j < n_dim; j++) {
      temp.at(i).at(j + n_dim) = temp.at(i).at(j + n_dim) / c;
    }
  }

  /* copy inverse part into mat_inv */
  for (int i = 0; i < n_dim; i++) {
    for (int j = 0; j < n_dim; j++) {
      mat_inv.at(i).at(j) = temp.at(i).at(j + n_dim);
    }
  }
}

unsigned long CLookUpTable::LookUp_ProgEnth(string    val_name_var,
                                            su2double *val_var,
                                            su2double val_prog,
                                            su2double val_enth, string name_prog, string name_enth){

  unsigned long exit_code=0;
  
  if (val_name_var.compare("NULL") == 0) {
    //cout << "variable is null, returning nothing" << endl;
    *val_var  = 0.0;
    exit_code = 0;
    return exit_code;
  }

  /* check if progress variable and enthalpy value is in table range */
  if ( val_prog >= limits_table_prog[0] && val_prog <= limits_table_prog[1]){
       //&&
       //val_enth >= limits_table_enth[0] && val_enth <= limits_table_enth[1] ){

    /* find the triangle that holds the (prog, enth) point */
    unsigned long id_triangle = trap_map_prog_enth.GetTriangle(val_prog, val_enth);

    if (IsInTriangle(val_prog, val_enth,id_triangle, name_prog, name_enth)) {

      /* get interpolation coefficients for point on triangle */
      vector<su2double> interp_coeffs(3);
      GetInterpCoeffs(val_prog, val_enth, interp_mat_inv_prog_enth.at(id_triangle), interp_coeffs);

      /* DEBUG: */
      vector<su2double> corner_prog(3, 0);
      vector<su2double> corner_enth(3, 0);
      for (int iPoint = 0; iPoint < N_POINTS_TRIANGLE; ++iPoint) {
        corner_prog.at(iPoint) =
            GetData(name_prog).at(triangles.at(id_triangle).at(iPoint));
        corner_enth.at(iPoint) =
            GetData(name_enth).at(triangles.at(id_triangle).at(iPoint));
      }

      *val_var = Interpolate(GetData(val_name_var), (triangles.at(id_triangle)),
                           interp_coeffs);
      exit_code = 0;
    } else {
      //cout << "lookupprogenth: in bounding box but outside of table!, (c,h)="<<val_prog<<", "<<val_enth << endl;
      unsigned long nearest_neighbor = FindNearestNeighborOnHull(val_prog, val_enth,name_prog, name_enth);
      *val_var = GetData(val_name_var).at(nearest_neighbor);
      exit_code = 1;
    }
  } else {
    //if (rank == MASTER_NODE) cout << "WARNING: LookUp_ProgEnth: lookup is outside of table bounding box, c,h = "<< val_prog<< " "<<val_enth<<endl;
    unsigned long nearest_neighbor = FindNearestNeighborOnHull(val_prog, val_enth,name_prog,name_enth);
    *val_var = GetData(val_name_var).at(nearest_neighbor);
    exit_code = 1;
  }
  return exit_code;
}


unsigned long CLookUpTable::LookUp_ProgEnth(vector<string>     &val_names_var,
                                            vector<su2double>  &val_vars,
                                            su2double           val_prog,
                                            su2double           val_enth, string name_prog,string name_enth) {

  vector<su2double*> look_up_data;

  for (long unsigned int i_var=0; i_var < val_vars.size(); ++i_var) {
    look_up_data.push_back(&val_vars[i_var]);
  }

  unsigned long exit_code = LookUp_ProgEnth(val_names_var, look_up_data, val_prog, val_enth, name_prog, name_enth);

  return exit_code;
}

unsigned long CLookUpTable::LookUp_ProgEnth(vector<string>     &val_names_var,
                                            vector<su2double*> &val_vars,
                                            su2double           val_prog,
                                            su2double           val_enth, string name_prog,string name_enth) {
  unsigned long exit_code=0;
  unsigned long id_triangle;
  unsigned long nearest_neighbor;
  vector<su2double> interp_coeffs(3);

  /* check if progress variable value is in progress variable table range 
   * and if enthalpy is in enthalpy table range */
  if ( val_prog >= limits_table_prog[0] && val_prog <= limits_table_prog[1] &&
       val_enth >= limits_table_enth[0] && val_enth <= limits_table_enth[1] ){

    /* if so, try to find the triangle that holds the (prog, enth) point */
    id_triangle = trap_map_prog_enth.GetTriangle(val_prog, val_enth);

    /* check if point is inside a triangle (if table domain is non-rectangular, 
     * the previous range check might be true but the point is still outside of the domain) */
    if (IsInTriangle(val_prog, val_enth, id_triangle, name_prog, name_enth)) {

      /* if so, get interpolation coefficients for point in  the triangle */
      GetInterpCoeffs(val_prog, val_enth, interp_mat_inv_prog_enth.at(id_triangle), interp_coeffs);

      /* exit_code 0 means point was in triangle */
      exit_code = 0;

    } else {
      //cout << "lookup_progenth: outside table range, c,h = "<< val_prog<< " "<<val_enth<<endl;
      /* if point is not inside a triangle (outside table domain) search nearest neighbor */
      nearest_neighbor = FindNearestNeighborOnHull(val_prog, val_enth, name_prog,name_enth);
      exit_code = 1;

    }

  } else {
    //if (rank == MASTER_NODE) cout << "WARNING: LookUp_ProgEnth: lookup is outside of table bounding box, c,h = "<< val_prog<< " "<<val_enth<<endl;

    /* if point is outside of table range, search nearest neighbor */
    nearest_neighbor = FindNearestNeighborOnHull(val_prog, val_enth,name_prog,name_enth);
    exit_code = 1;

  }

  /* loop over variable names and interpolate / get values */
  for (long unsigned int i_var = 0; i_var < val_names_var.size(); ++i_var) {

    if (val_names_var.at(i_var).compare("NULL")==0) {

      *val_vars.at(i_var) = 0.0; 
    } else {

    if (exit_code == 0)
      *val_vars.at(i_var) = Interpolate(GetData(val_names_var.at(i_var)), (triangles.at(id_triangle)), interp_coeffs);
    else
      *val_vars.at(i_var) = GetData(val_names_var.at(i_var)).at(nearest_neighbor);

    }
  }
  return exit_code;
  
}




  void CLookUpTable::GetInterpCoeffs(su2double                     val_x,
                                     su2double                     val_y,
                                     vector< vector< su2double > > &interp_mat_inv,
                                     vector<su2double>             &interp_coeffs) {

    vector<su2double> query_vector;
    query_vector.push_back(1);
    query_vector.push_back(val_x);
    query_vector.push_back(val_y);

    su2double d;
    for (int i = 0; i < 3; i++) {
      d = 0;
      for (int j = 0; j < 3; j++) {
        d = d + interp_mat_inv.at(i).at(j) * query_vector.at(j);
      }
      interp_coeffs.at(i) = d;
    }

}

su2double CLookUpTable::Interpolate(const vector<su2double> &val_samples,
                                    vector<unsigned long>   &val_triangle,
                                    vector<su2double>       &val_interp_coeffs) {

  su2double result         = 0;
  su2double z              = 0;

  for (int i_point = 0; i_point < N_POINTS_TRIANGLE; i_point++) {
    z = val_samples.at(val_triangle.at(i_point));
    result += val_interp_coeffs.at(i_point) * z;
  }

  return result;

}

unsigned long CLookUpTable::FindNearestNeighborOnHull(su2double val_prog,
                                                      su2double val_enth, string name_prog, string name_enth){

  su2double min_distance  = 1.e99;
  su2double next_distance = 1.e99;
  su2double next_prog_norm;
  su2double next_enth_norm;
  unsigned long neighbor_id = 0;

  const vector<su2double> &prog_table = GetData(name_prog);
  const vector<su2double> &enth_table = GetData(name_enth);

  su2double norm_coeff_prog = 1.       / (limits_table_prog[1] - limits_table_prog[0]);
  su2double norm_coeff_enth = 1.       / (limits_table_enth[1] - limits_table_enth[0]);
  su2double val_prog_norm   = val_prog / (limits_table_prog[1] - limits_table_prog[0]);
  su2double val_enth_norm   = val_enth / (limits_table_enth[1] - limits_table_enth[0]);

  for (unsigned long i_point = 0; i_point < n_hull_points; ++i_point){

    next_prog_norm = prog_table.at(hull.at(i_point)) * norm_coeff_prog;
    next_enth_norm = enth_table.at(hull.at(i_point)) * norm_coeff_enth;
    
    next_distance  = sqrt( pow(val_prog_norm - next_prog_norm, 2) + 
                          pow(val_enth_norm - next_enth_norm, 2));

    if (next_distance < min_distance){
      min_distance = next_distance;
      neighbor_id  = hull.at(i_point);
    }
  }
  return neighbor_id;
}

bool CLookUpTable::IsInTriangle(su2double val_prog, su2double val_enth, unsigned long val_id_triangle, string name_prog, string name_enth){

  su2double tri_prog_0 = GetData(name_prog).at(triangles.at(val_id_triangle).at(0));
  su2double tri_enth_0 = GetData(name_enth).at(triangles.at(val_id_triangle).at(0));
  su2double tri_prog_1 = GetData(name_prog).at(triangles.at(val_id_triangle).at(1));
  su2double tri_enth_1 = GetData(name_enth).at(triangles.at(val_id_triangle).at(1));
  su2double tri_prog_2 = GetData(name_prog).at(triangles.at(val_id_triangle).at(2));
  su2double tri_enth_2 = GetData(name_enth).at(triangles.at(val_id_triangle).at(2));

  su2double area_tri = TriArea(tri_prog_0, tri_enth_0, tri_prog_1, tri_enth_1, tri_prog_2, tri_enth_2);
  su2double area_0   = TriArea(val_prog,   val_enth,   tri_prog_1, tri_enth_1, tri_prog_2, tri_enth_2);
  su2double area_1   = TriArea(tri_prog_0, tri_enth_0, val_prog,   val_enth,   tri_prog_2, tri_enth_2);
  su2double area_2   = TriArea(tri_prog_0, tri_enth_0, tri_prog_1, tri_enth_1, val_prog,   val_enth  );

  //if ( abs(area_tri - (area_0 + area_1 + area_2)) >= area_tri * 1e-10 ){
  //  cout << "id triangle = "<<val_id_triangle<<endl;
  //  cout << "p="<<val_prog<<", "<<val_enth << endl;
  //  cout << "p1="<<tri_prog_0<<", "<<tri_enth_0 << endl;
  //  cout << "p2="<<tri_prog_1<<", "<<tri_enth_1 << endl;
  //  cout << "p3="<<tri_prog_2<<", "<<tri_enth_2 << endl;
  //}

  return ( abs(area_tri - (area_0 + area_1 + area_2)) < area_tri * 1e-10 );

}
