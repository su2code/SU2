/*!
 * \file TrapezoidalMap.hpp
 * \brief Implementation of the trapezoidal map for tabulation and lookup of fluid properties
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
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

class CTrapezoidalMap {
 protected:
  /* The unique values of x which exist in the data */
  vector<su2double> unique_bands_x;

  vector<vector<su2double> > edge_limits_x;
  vector<vector<su2double> > edge_limits_y;

  vector<vector<unsigned long> > edge_to_triangle;

  /* The value that each edge which intersects the band takes within that
   * same band. Used to sort the edges */
  vector<vector<pair<su2double, unsigned long> > > y_edge_at_band_mid;

 public:
  CTrapezoidalMap();

  CTrapezoidalMap(const vector<su2double>& samples_x, const vector<su2double>& samples_y,
                  const vector<vector<unsigned long> >& edges, const vector<vector<unsigned long> >& edge_to_triangle);

  ~CTrapezoidalMap(void);

  void Search_Band_For_Edge(su2double val_x, su2double val_y);

  unsigned long GetTriangle(su2double val_x, su2double val_y);

  pair<unsigned long, unsigned long> GetBand(su2double val_x);
  pair<unsigned long, unsigned long> GetEdges(pair<unsigned long, unsigned long> val_band, su2double val_x,
                                              su2double val_y);

  inline bool IsInsideHullX(su2double val_x) {
    return (val_x >= unique_bands_x.front()) && (val_x <= unique_bands_x.back());
  }
};