#pragma once

#include <string>
#include <vector>
#include <fstream>
#include "../../Common/include/parallelization/mpi_structure.hpp"
/*#include "../../Common/include/datatypes/primitive_structure.hpp"*/

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
  vector< string > names_var;
  
  /* !brief
   * Holds all data stored in the table. First index addresses the variable
   * while second index addresses the point.
   */
  vector< vector< su2double > > table_data;
  
  vector< vector< unsigned long > > triangles;

  vector< unsigned long > hull;

  string SkipToFlag(ifstream *file_stream, string flag);

  inline void SetTypeLUT(string value)    { type_lut      = value; }
  inline void SetVersionLUT(string value) { version_lut   = value; }
  inline void SetNPoints(unsigned long  value)       { n_points      = value; }
  inline void SetNTriangles(unsigned long  value)    { n_triangles   = value; }
  inline void SetNHullPoints(unsigned long  value)   { n_hull_points = value; }
  inline void SetNVariables(unsigned long  value)    { n_variables   = value; }
  
  inline void PushNameVar(string value) { names_var.push_back(value); }
  inline void PopNameVar() { names_var.pop_back(); }

  inline void AllocMemData() {
    table_data.resize(GetNVariables(), vector< su2double >(GetNPoints()));
  }
  
  inline void AllocMemTriangles() {
    triangles.resize(GetNTriangles(), vector< unsigned long >(3));
  }
  
  inline void AllocMemHull() {
    hull.resize(GetNHullPoints());
  }

 public:

  CFileReaderLUT();

  inline string        GetTypeLUT()       { return type_lut; }
  inline string        GetVersionLUT()    { return version_lut; }
  inline string        GetVersionReader() { return version_reader; }
  inline unsigned long GetNPoints()       { return n_points; }
  inline unsigned long GetNTriangles()    { return n_triangles; }
  inline unsigned long GetNHullPoints()   { return n_hull_points; }
  inline unsigned long GetNVariables()    { return n_variables; }

  inline const vector< string >                  &GetNamesVar()  const { return names_var; }

  inline const vector< vector< su2double > >     &GetTableData() const { return table_data; }

  inline const vector< vector< unsigned long > > &GetTriangles() const { return triangles; };

  inline const vector< unsigned long >           &GetHull()      const { return hull; };

  void ReadRawDRG(string file_name);

};
