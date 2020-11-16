#pragma once

#include<vector>
#include<iostream>
#include<fstream>
#include<iterator>

using namespace std;

class CTrapezoidalMap {
protected:

  /* The unique values of x which exist in the data */
  vector< su2double > unique_bands_x;

  vector< vector< su2double > > edge_limits_x;
  vector< vector< su2double > > edge_limits_y;


  vector< vector< unsigned long > > edge_to_triangle;

  /* The value that each edge which intersects the band takes within that
   * same band. Used to sort the edges */
  vector< vector< pair< su2double, unsigned long > > > y_edge_at_band_mid;

  
public:
  CTrapezoidalMap();

  CTrapezoidalMap(const vector< su2double >                &samples_x,
                  const vector< su2double >                &samples_y,
                  const vector< vector< unsigned long > >  &edges,
                  const vector< vector< unsigned long > >  &edge_to_triangle);

  ~CTrapezoidalMap(void);

  /* nijso: todo testing of writing the trapezoidal map to file */
  //void writeVec1file(string path, vector<su2double>& myVector);
  //void readVec1file(string path, vector<su2double>& myVector);
  //void writeVec2file(string path, vector<vector<su2double> >& myVector);
  //void readVec2file(string path, vector<vector<double> >& myVector);
  //void writeVec3file(string path, vector<vector<unsigned long> >& myVector);
  //void readVec3file(string path, vector<vector<unsigned long> > & myVector);
  //void writeVec4file(string path, vector<vector<pair<su2double,unsigned long> > >& myVector);
  //void readVec4file(string path, vector<vector<pair<su2double,unsigned long> > >& myVector);
  //void writeVec(void);
  //void readVec(void);

  void Search_Band_For_Edge(su2double val_x, su2double val_y);

  unsigned long GetTriangle(su2double val_x, su2double val_y);

  pair<unsigned long, unsigned long> GetBand(su2double val_x);
  pair<unsigned long, unsigned long> GetEdges(pair<unsigned long, unsigned long> val_band,
                                              su2double val_x, su2double val_y);

  inline bool IsInsideHullX(su2double val_x) {
    return ( val_x >= unique_bands_x.front() ) && ( val_x <= unique_bands_x.back() );
  }
};
