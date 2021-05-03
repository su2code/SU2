
#include "../../Common/include/option_structure.hpp"

#include "../include/numerics/CTrapezoidalMap.hpp"

using namespace std;

  CTrapezoidalMap::CTrapezoidalMap(vector< su2double > const &samples_x,
                                   vector< su2double > const &samples_y,
                                   vector< vector< unsigned long > > const &edges,
                                   vector< vector< unsigned long > > const &val_edge_to_triangle){

  int rank = SU2_MPI::GetRank();
  //bool loadmap=false; /* load the trapezoidal map from file */
  clock_t build_start = clock();


  //ifstream edgetrianglefile("edge_to_triangle.bin");
  //ifstream edgelimitsxfile("edge_limits_x.bin");
  //ifstream edgelimitsyfile("edge_limits_x.bin");
  //ifstream uniquebandsfile("unique_bands.bin");
  //if (edgetrianglefile){
  //  cout << "edge_to_triangle.bin file exists"<<endl;
 // }
  //if (uniquebandsfile){
  //  cout << "unique_bands.bin file exists"<<endl;
  //}
  //if (edgelimitsxfile){
  //  cout << "edge_limits_x.bin file exists"<<endl;
  //}
  //if (edgelimitsyfile){
  //  cout << "edge_limits_y.bin file exists"<<endl;
  //}


  //edge_to_triangle = &val_edge_to_triangle;
  edge_to_triangle = vector< vector<unsigned long> > (val_edge_to_triangle);

  unique_bands_x = vector< su2double >(samples_x);

  /* sort x_bands and make them unique */
  sort(unique_bands_x.begin(), unique_bands_x.end());

  vector< su2double >::iterator iter;
  iter = unique(unique_bands_x.begin(), unique_bands_x.end());

  unique_bands_x.resize(distance(unique_bands_x.begin(), iter));

  edge_limits_x.resize(edges.size(), vector< su2double >(2, 0));
  edge_limits_y.resize(edges.size(), vector< su2double >(2, 0));

  /* store x and y values of each edge in a vector for a slight speed up
   * as it prevents some uncoalesced accesses */
  for (unsigned long j = 0; j < edges.size(); j++) {
    edge_limits_x[j][0] = samples_x[edges[j][0]];
    edge_limits_x[j][1] = samples_x[edges[j][1]];
    edge_limits_y[j][0] = samples_y[edges[j][0]];
    edge_limits_y[j][1] = samples_y[edges[j][1]];
  }

  /* number of bands */
  unsigned long n_bands_x = unique_bands_x.size() - 1;
  /* band index */
  unsigned long i_band = 0;
  /* number of edges */
  unsigned long n_edges = edges.size();
  /* edge index */
  unsigned long i_edge = 0;
  /* counter for edges intersects */
  unsigned long n_intersects = 0;
  /* lower and upper x value of each band */
  su2double band_lower_x = 0;
  su2double band_upper_x = 0;

  su2double x_0;
  su2double y_0;
  su2double dy_edge;
  su2double dx_edge;
  su2double x_band_mid;

  //cout << "y_edge_at_band_mid start, number of x-bands = "<< n_bands_x << endl;
  
  /* y values of all intersecting edges for every band */
  y_edge_at_band_mid.resize(unique_bands_x.size() - 1);

  /* loop over bands */
  while (i_band < n_bands_x) {
    band_lower_x = unique_bands_x[i_band];
    band_upper_x = unique_bands_x[i_band + 1];
    i_edge       = 0;
    n_intersects = 0;

    /* loop over edges and determine which edges appear in current band */
    while (i_edge < n_edges) {

      /* check if edge intersects the band 
       * (vertical edges are automatically discarded) */
      if (((edge_limits_x[i_edge][0] <= band_lower_x) and
           (edge_limits_x[i_edge][1] >= band_upper_x)) or
          ((edge_limits_x[i_edge][1] <= band_lower_x) and
           (edge_limits_x[i_edge][0] >= band_upper_x))) {

        y_edge_at_band_mid[i_band].push_back(
            make_pair(0.0, 0));

        x_0 = edge_limits_x[i_edge][0];
        y_0 = edge_limits_y[i_edge][0];

        dy_edge = edge_limits_y[i_edge][1] - edge_limits_y[i_edge][0];
        dx_edge = edge_limits_x[i_edge][1] - edge_limits_x[i_edge][0];
        x_band_mid = (band_lower_x + band_upper_x) / 2.0;

        y_edge_at_band_mid[i_band][n_intersects].first =
            y_0 + dy_edge / dx_edge * (x_band_mid - x_0);

        /* save edge index so it can later be recalled when searching */
        y_edge_at_band_mid[i_band][n_intersects].second =
            i_edge;                

        n_intersects++;
      }
      i_edge++;
    }

    /* sort edges by their y values.
     * note that these y values are unique (i.e. edges cannot
     * intersect in a band) */
    sort(y_edge_at_band_mid[i_band].begin(),
         y_edge_at_band_mid[i_band].end());

    i_band++;
  }

  su2double duration = ((su2double)clock() - (su2double)build_start) /
                       ((su2double)CLOCKS_PER_SEC);

  if (rank == MASTER_NODE)
    cout << "Construction of trapezoidal map took " << duration << " seconds\n" << endl;

  }

  CTrapezoidalMap::CTrapezoidalMap() {
  }

  CTrapezoidalMap::~CTrapezoidalMap(void) {}

/*
void CTrapezoidalMap::writeVec(void)
{
    const string trapfilename="trapmap.bin"; 

    int i,j;
    unsigned int size1;
    su2double x;
    unsigned long y;
    cout<<"writing trapezoidal map to file"<<endl;
    if(trapfilename != "")
    {
    ofstream myfile(trapfilename,std::ios::trunc | std::ios::binary);
    //myfile.open(trapfilename);
        if(!myfile.is_open()){
          cout <<"could not open trapezoidalmap file for writing" << endl;
        } else {
          size1=y_edge_at_band_mid.size();
          cout << "size of points = "<<size1<<endl;
          myfile.write(reinterpret_cast<char*>(&size1),sizeof(size1));

          for(i=0; i <y_edge_at_band_mid.size();i++ )
          {
            size1=y_edge_at_band_mid[i].size();
            myfile.write(reinterpret_cast<char*>(&size1),sizeof(size1));
            for(j=0; j <y_edge_at_band_mid[i].size();j++ ) {
              x = y_edge_at_band_mid[i][j].first;           
              myfile.write(reinterpret_cast<char*>(&x),sizeof(x));
              y = y_edge_at_band_mid[i][j].second;          
              myfile.write(reinterpret_cast<char*>(&y),sizeof(y));
            }  
          }
        }
        myfile.close();
    }
}
*/

/*
void CTrapezoidalMap::readVec(void){
  //readVec1file("unique_bands_x.bin",unique_bands_x);
  //readVec2file("edge_limits_x.bin",edge_limits_x);
  //readVec2file("edge_limits_y.bin",edge_limits_y);
  //readVec3file("edge_to_triangle.bin",edge_to_triangle);
}
*/

/*
void CTrapezoidalMap::writeVec2file(string path, vector<vector<su2double> >& myVector)
{
    ofstream FILE(path, std::ios::out | std::ofstream::binary);

    // Store size of the outer vector
    int s1 = myVector.size();
    FILE.write(reinterpret_cast<const char *>(&s1), sizeof(s1));    

    // Now write each vector one by one
    for (auto& v : myVector) {         
        // Store its size
        int size = v.size();
        FILE.write(reinterpret_cast<const char *>(&size), sizeof(size));

        // Store its contents
        FILE.write(reinterpret_cast<const char *>(&v[0]), v.size()*sizeof(su2double));
    }
    FILE.close();   
}
*/

/*
void CTrapezoidalMap::writeVec3file(string path, vector<vector<unsigned long> >& myVector)
{
    ofstream FILE(path, std::ios::out | std::ofstream::binary);

    int s1 = myVector.size();
    FILE.write(reinterpret_cast<const char *>(&s1), sizeof(s1));    
    for (auto& v : myVector) {         
        int size = v.size();
        FILE.write(reinterpret_cast<const char *>(&size), sizeof(size));
        FILE.write(reinterpret_cast<const char *>(&v[0]), v.size()*sizeof(unsigned long));
    }
    FILE.close();   
}
*/

/*
void CTrapezoidalMap::writeVec1file(string path, vector<su2double>& myVector)
{
    ofstream FILE(path, std::ios::out | std::ofstream::binary);

    int s1 = myVector.size();
    FILE.write(reinterpret_cast<const char *>(&s1), sizeof(s1));    
    for (auto& v : myVector) {         
        // Store its contents
        FILE.write(reinterpret_cast<const char *>(&v), sizeof(su2double));
    }
    FILE.close();   
}
*/

/*
void CTrapezoidalMap::readVec2file(string path,  vector<vector<double> >& myVector)
{
    ifstream FILE(path, std::ios::in | std::ifstream::binary);

    int size = 0;
    FILE.read(reinterpret_cast<char *>(&size), sizeof(size));
    myVector.resize(size);
    for (int n = 0; n < size; ++n) {
        int size2 = 0;
        FILE.read(reinterpret_cast<char *>(&size2), sizeof(size2));
        //su2double f;        
        float f;        
        for ( int k = 0; k < size2; ++k ) {
            FILE.read(reinterpret_cast<char *>(&f), sizeof(f));
            myVector[n].push_back(f);   
        }
    }
}
*/

/*
void CTrapezoidalMap::readVec3file(string path,  vector<vector<unsigned long> >& myVector)
{
    ifstream FILE(path, std::ios::in | std::ifstream::binary);

    int size = 0;
    FILE.read(reinterpret_cast<char *>(&size), sizeof(size));
    myVector.resize(size);
    for (int n = 0; n < size; ++n) {
        int size2 = 0;
        FILE.read(reinterpret_cast<char *>(&size2), sizeof(size2));
        //su2double f;        
        unsigned long f;        
        for ( int k = 0; k < size2; ++k ) {
            FILE.read(reinterpret_cast<char *>(&f), sizeof(f));
            myVector[n].push_back(f);   
        }
    }
}
*/

/*
void CTrapezoidalMap::readVec4file(string path,  vector<vector<pair<su2double,unsigned long> > >& myVector)
{
    
}
*/


/*
//void CTrapezoidalMap::writeVec4file(string path,  vector<vector<pair<su2double,unsigned long> > >& myVector)
//{
void writeVec4(std::string filename, std::vector<std::pair<std::string, std::vector<int>>> dataset){
    // Make a CSV file with one or more columns of integer values
    // Each column of data is represented by the pair <column name, column data>
    //   as std::pair<std::string, std::vector<int>>
    // The dataset is represented as a vector of these columns
    // Note that all columns should be the same size
    
    // Create an output filestream object
    std::ofstream myFile(filename);
    
    // Send column names to the stream
    for(int j = 0; j < dataset.size(); ++j)
    {
        myFile << dataset.at(j).first;
        if(j != dataset.size() - 1) myFile << ","; // No comma at end of line
    }
    myFile << "\n";
    
    // Send data to the stream
    for(int i = 0; i < dataset.at(0).second.size(); ++i)
    {
        for(int j = 0; j < dataset.size(); ++j)
        {
            myFile << dataset.at(j).second.at(i);
            if(j != dataset.size() - 1) myFile << ","; // No comma at end of line
        }
        myFile << "\n";
    }
    
    // Close the file
    myFile.close();
}    
*/

/*
void CTrapezoidalMap::readVec1file(string path,  vector<su2double>& myVector)
{
    ifstream FILE(path, std::ios::in | std::ifstream::binary);

    int size = 0;
    FILE.read(reinterpret_cast<char *>(&size), sizeof(size));
    //myVector.resize(size);
    for (int n = 0; n < size; ++n) {
        //int size2 = 0;
        //FILE.read(reinterpret_cast<char *>(&size2), sizeof(size2));
        su2double f;        
        //for ( int k = 0; k < size2; ++k ) {
        FILE.read(reinterpret_cast<char *>(&f), sizeof(f));
        myVector.push_back(f);   
    }
}
*/

  unsigned long CTrapezoidalMap::GetTriangle(su2double val_x, su2double val_y) {

    /* find x band in which val_x sits */
    pair<unsigned long, unsigned long> band = GetBand(val_x);

    /* within that band, find edges which enclose the (val_x, val_y) point */
    pair<unsigned long, unsigned long> edges = GetEdges(band, val_x, val_y);

    /* identify the triangle using the two edges */
    //vector<unsigned long> triangles_edge_low = edge_to_triangle->at(edges.first);
    //vector<unsigned long> triangles_edge_up  = edge_to_triangle->at(edges.second);
    vector<unsigned long> triangles_edge_low = edge_to_triangle.at(edges.first);
    vector<unsigned long> triangles_edge_up  = edge_to_triangle.at(edges.second);

    sort(triangles_edge_low.begin(), triangles_edge_low.end());

    sort(triangles_edge_up.begin(), triangles_edge_up.end());

    // The intersection of the faces to which upper or lower belongs is
    // the face that both belong to.
    vector<unsigned long> triangle(1);
    set_intersection(triangles_edge_up.begin(), triangles_edge_up.end(),
                     triangles_edge_low.begin(), triangles_edge_low.end(),
                     triangle.begin());

    return triangle.at(0);
  }

  pair<unsigned long, unsigned long> CTrapezoidalMap::GetBand(su2double val_x) {
    su2double x_low;
    su2double x_mid;
    su2double x_up;

    unsigned long i_low = 0;
    unsigned long i_mid = 0;
    unsigned long i_up = 0;

    /* start search at table limits */
    i_up = unique_bands_x.size() - 1;
    i_low = 0;

    /* check if val_x is in bounds of the table */
    // if (val_x < unique_bands_x.front() or val_x > unique_bands_x.back())
    //   SU2_MPI::Error("Table-look-up is out of bounds.", CURRENT_FUNCTION);

  if ( val_x < unique_bands_x.front() ) val_x = unique_bands_x.front();
  if ( val_x > unique_bands_x.back()  ) val_x = unique_bands_x.back();

  /* the next loop implements a binary search and computes i_low and i_up 
   * which are the band indices that include val_x */
  do {

    i_mid = (i_up + i_low) / 2;

    x_mid = unique_bands_x.at(i_mid);
    x_low = unique_bands_x.at(i_low);
    x_up  = unique_bands_x.at(i_up);

    /* check and restart the search on the low end */
    if ( (val_x < x_low) and (i_low > 0) ) {

      i_up  = i_low;
      i_low = i_low / 2;

    /* check and restart the search on the upper end */
    } else if ( (val_x > x_up) and (i_up < (unique_bands_x.size() - 1)) ) {
      
      i_low = i_up;
      i_up  = (i_up + (unique_bands_x.size() - 1)) / 2;

    /*  continue with regular binary search */
    } else if (val_x < x_mid) {
      i_up = i_mid;

    } else if (val_x > x_mid) {
      i_low = i_mid;

    } else if (x_mid == val_x) {
      i_low = i_mid;
      i_up  = i_low + 1;
      break;
    }

  } while (i_up - i_low > 1);

  return make_pair(i_low, i_up);
}

pair< unsigned long, unsigned long > CTrapezoidalMap::GetEdges(pair< unsigned long, unsigned long > val_band, su2double val_x, su2double val_y) {

  su2double next_y;
  su2double y_edge_low;
  su2double y_edge_up;
  su2double x_edge_low;
  su2double x_edge_up;

  unsigned long i_band_low = val_band.first;

  unsigned long next_edge;

  unsigned long j_low = 0;
  unsigned long j_mid = 0;
  unsigned long j_up  = 0;

  j_up  = y_edge_at_band_mid.at(i_band_low).size() - 1;
  j_low = 0;

  while (j_up - j_low > 1) {
    
    j_mid = (j_up + j_low) / 2;

    // Select the edge associated with the x band (i_band_low)
    // Search for the RunEdge in the y direction (second value is index of
    // edge)
    next_edge = y_edge_at_band_mid.at(i_band_low).at(j_mid).second;

    y_edge_low = edge_limits_y.at(next_edge).at(0);
    y_edge_up  = edge_limits_y.at(next_edge).at(1);
    x_edge_low = edge_limits_x.at(next_edge).at(0);
    x_edge_up  = edge_limits_x.at(next_edge).at(1);

    // The search variable in j should be interpolated in i as well
    next_y = y_edge_low + (y_edge_up - y_edge_low) / (x_edge_up - x_edge_low) * (val_x - x_edge_low);

    if (next_y > val_y) {
      j_up = j_mid;

    } else if (next_y < val_y) {
      j_low = j_mid;

    } else if (next_y == val_y) {
      j_low = j_mid;
      j_up = j_low + 1;
      break;
    }
  }

  unsigned long edge_low = y_edge_at_band_mid.at(i_band_low).at(j_low).second;
  unsigned long edge_up  = y_edge_at_band_mid.at(i_band_low).at(j_up).second;

  return make_pair(edge_low, edge_up); 
}
