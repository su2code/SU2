/*!
 * \file CTrapezoidalMap.cpp
 * \brief Implementation of the trapezoidal map for tabulation and lookup of fluid properties
 * \author D. Mayer, T. Economon, N. Beishuizen
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

#include <array>
#include <iomanip>

#include "../../Common/include/option_structure.hpp"
#include "../../Common/include/containers/CTrapezoidalMap.hpp"

using namespace std;

/* Trapezoidal map implementation. Reference:
 * M. de Berg, O. Cheong M. van Kreveld, M. Overmars,
 * Computational Geometry, Algorithms and Applications pp. 121-146 (2008)
 * NOTE: the current implementation is actually the simpler 'slab' approach.
 */
CTrapezoidalMap::CTrapezoidalMap(const su2double* samples_x, const su2double* samples_y, const unsigned long size,
                                 vector<std::array<unsigned long, 2> > const& edges,
                                 su2vector<vector<unsigned long> > const& val_edge_to_triangle, bool display) {
  int rank = SU2_MPI::GetRank();
  su2double startTime = SU2_MPI::Wtime();

  edge_to_triangle = su2vector<vector<unsigned long> >(val_edge_to_triangle);

  unique_bands_x.assign(samples_x, samples_x + size);

  /* sort x_bands and make them unique */
  sort(unique_bands_x.begin(), unique_bands_x.end());

  auto iter = unique(unique_bands_x.begin(), unique_bands_x.end());

  unique_bands_x.resize(distance(unique_bands_x.begin(), iter));

  edge_limits_x.resize(edges.size(), 2);
  edge_limits_y.resize(edges.size(), 2);

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
  unsigned long j_edge = 0;
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

  /* y values of all intersecting edges for every band */
  y_edge_at_band_mid.resize(unique_bands_x.size() - 1);

  /* loop over bands */
  while (i_band < n_bands_x) {
    band_lower_x = unique_bands_x[i_band];
    band_upper_x = unique_bands_x[i_band + 1];
    i_edge = 0;
    n_intersects = 0;

    /* loop over edges and determine which edges appear in current band */
    while (i_edge < n_edges) {
      /* check if edge intersects the band
       * (vertical edges are automatically discarded) */
      if (((edge_limits_x[i_edge][0] <= band_lower_x) and (edge_limits_x[i_edge][1] >= band_upper_x)) or
          ((edge_limits_x[i_edge][1] <= band_lower_x) and (edge_limits_x[i_edge][0] >= band_upper_x))) {
        y_edge_at_band_mid[i_band].emplace_back(0.0, 0);

        x_0 = edge_limits_x[i_edge][0];
        y_0 = edge_limits_y[i_edge][0];

        dy_edge = edge_limits_y[i_edge][1] - edge_limits_y[i_edge][0];
        dx_edge = edge_limits_x[i_edge][1] - edge_limits_x[i_edge][0];
        x_band_mid = (band_lower_x + band_upper_x) / 2.0;

        y_edge_at_band_mid[i_band][n_intersects].first = y_0 + dy_edge / dx_edge * (x_band_mid - x_0);

        /* save edge index so it can later be recalled when searching */
        y_edge_at_band_mid[i_band][n_intersects].second = i_edge;

        n_intersects++;
      }
      i_edge++;
    }

    /* sort edges by their y values.
     * note that these y values are unique (i.e. edges cannot
     * intersect in a band) */
    sort(y_edge_at_band_mid[i_band].begin(), y_edge_at_band_mid[i_band].end());

    i_band++;
  }

  su2double stopTime = SU2_MPI::Wtime();

  /* calculate size of trapezoidal map components */
  double size_unique_bands = sizeof(su2double) * unique_bands_x.size() / 1e6;
  double size_edge_limits_x = sizeof(su2double) * edge_limits_x.size() * 2 / 1e6;
  double size_edge_limits_y = sizeof(su2double) * edge_limits_y.size() * 2 / 1e6;

  double size_edge_to_triangle = 0;
  for (i_edge = 0; i_edge < edge_to_triangle.size(); i_edge++)
    for (j_edge = 0; j_edge < edge_to_triangle[i_edge].size(); j_edge++)
      size_edge_to_triangle += sizeof(unsigned long) / 1e6;

  double size_y_edge_at_band_mid = 0;
  for (unsigned long i_y = 0; i_y < y_edge_at_band_mid.size(); i_y++)
    for (unsigned long j_y = 0; j_y < y_edge_at_band_mid[i_y].size(); j_y++)
      size_y_edge_at_band_mid += sizeof(su2double) / 1e6 + sizeof(unsigned long) / 1e6;

  memory_footprint =
      size_unique_bands + size_edge_limits_x + size_edge_limits_y + size_edge_to_triangle + size_y_edge_at_band_mid;

  /* print size of trapezoidal map components to screen */
  if ((rank == MASTER_NODE) && display) {
    cout << setfill(' ');
    cout << "\n" << endl;
    cout << "+------------------------------------------------------------------+\n";
    cout << "|                       Trapezoidal map info                       |\n";
    cout << "+------------------------------------------------------------------+" << endl;

    cout << "| Time to construct trapezoidal map:    " << setw(22) << right << stopTime - startTime << " sec"
         << " |" << endl;
    cout << "| Size of unique_bands in memory:       " << setw(22) << size_unique_bands << " MB "
         << " |" << endl;
    cout << "| Size of edge_limits_x in memory:      " << setw(22) << size_edge_limits_x << " MB "
         << " |" << endl;
    cout << "| Size of edge_limits_y in memory:      " << setw(22) << size_edge_limits_y << " MB "
         << " |" << endl;
    cout << "| Size of edge_to_triangle in memory:   " << setw(22) << size_edge_to_triangle << " MB "
         << " |" << endl;
    cout << "| Size of y_edge_at_band_mid in memory: " << setw(22) << size_y_edge_at_band_mid << " MB "
         << " |" << endl;
    cout << "| Total:                                " << setw(22) << memory_footprint << " MB "
         << " |" << endl;
    cout << "+------------------------------------------------------------------+" << endl;
    cout << "\n" << endl;
  }
}

unsigned long CTrapezoidalMap::GetTriangle(su2double val_x, su2double val_y) {
  /* find x band in which val_x sits */
  pair<unsigned long, unsigned long> band = GetBand(val_x);

  /* within that band, find edges which enclose the (val_x, val_y) point */
  pair<unsigned long, unsigned long> edges = GetEdges(band, val_x, val_y);

  /* identify the adjacent triangles using the two edges */
  std::array<unsigned long, 2> triangles_edge_low;
  for (unsigned long i = 0; i < edge_to_triangle[edges.first].size(); i++)
    triangles_edge_low[i] = edge_to_triangle[edges.first][i];

  std::array<unsigned long, 2> triangles_edge_up;
  for (unsigned long i = 0; i < edge_to_triangle[edges.second].size(); i++)
    triangles_edge_up[i] = edge_to_triangle[edges.second][i];

  sort(triangles_edge_low.begin(), triangles_edge_low.end());
  sort(triangles_edge_up.begin(), triangles_edge_up.end());

  /* The intersection of the faces to which upper or lower belongs is the face that both belong to. */
  vector<unsigned long> triangle;
  set_intersection(triangles_edge_up.begin(), triangles_edge_up.end(), triangles_edge_low.begin(),
                   triangles_edge_low.end(), std::back_inserter(triangle));

  /*--- We failed to find an intersection, so take the lower triangle inside the band enclosing the point---*/
  if (triangle.size() < 1) {
    triangle.resize(1, triangles_edge_low[0]);
  }

  return triangle[0];
}

pair<unsigned long, unsigned long> CTrapezoidalMap::GetBand(su2double val_x) {
  unsigned long i_low = 0;
  unsigned long i_up = 0;

  /* check if val_x is in x-bounds of the table, if not then project val_x to either x-min or x-max */
  if (val_x < unique_bands_x.front()) val_x = unique_bands_x.front();
  if (val_x > unique_bands_x.back()) val_x = unique_bands_x.back();

  std::pair<std::vector<su2double>::iterator, std::vector<su2double>::iterator> bounds;
  bounds = std::equal_range(unique_bands_x.begin(), unique_bands_x.end(), val_x);

  /*--- if upper bound = 0, then use the range [0,1] ---*/
  i_up = max<unsigned long>(1, bounds.first - unique_bands_x.begin());
  i_low = i_up - 1;

  return make_pair(i_low, i_up);
}

pair<unsigned long, unsigned long> CTrapezoidalMap::GetEdges(pair<unsigned long, unsigned long> val_band,
                                                             su2double val_x, su2double val_y) const {
  su2double next_y;
  su2double y_edge_low;
  su2double y_edge_up;
  su2double x_edge_low;
  su2double x_edge_up;

  unsigned long i_band_low = val_band.first;

  unsigned long next_edge;

  unsigned long j_low = 0;
  unsigned long j_mid = 0;
  unsigned long j_up = 0;

  j_up = y_edge_at_band_mid[i_band_low].size() - 1;
  j_low = 0;

  while (j_up - j_low > 1) {
    j_mid = (j_up + j_low) / 2;

    // Select the edge associated with the x band (i_band_low)
    // Search for the RunEdge in the y direction (second value is index of
    // edge)
    next_edge = y_edge_at_band_mid[i_band_low][j_mid].second;

    y_edge_low = edge_limits_y[next_edge][0];
    y_edge_up = edge_limits_y[next_edge][1];
    x_edge_low = edge_limits_x[next_edge][0];
    x_edge_up = edge_limits_x[next_edge][1];

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

  unsigned long edge_low = y_edge_at_band_mid[i_band_low][j_low].second;
  unsigned long edge_up = y_edge_at_band_mid[i_band_low][j_up].second;

  return make_pair(edge_low, edge_up);
}
