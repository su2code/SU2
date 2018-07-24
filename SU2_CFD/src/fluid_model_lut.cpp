/*!
 * fluid_model_lut.cpp
 * \brief Source of the look-up table model.
 * \author S. Vitale, A. Rubino
 * \version 4.1.2 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/fluid_model_lut.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <string>
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "stdlib.h"
#include "stdio.h"
#include <iomanip>
#include <vector>

using namespace std;

CTrapezoidalMap::CTrapezoidalMap() {

	UpperEdge = 0;  LowerEdge = 0;
	LowerI    = 0;  UpperI    = 0;  MiddleI   = 0;
	LowerJ    = 0;  UpperJ    = 0; 	MiddleJ  = 0;
	rank = MASTER_NODE;

}
CTrapezoidalMap::~CTrapezoidalMap(void) {}
CTrapezoidalMap::CTrapezoidalMap(vector<su2double> const &x_samples,
		vector<su2double> const &y_samples,
		vector<vector<unsigned long> > const &unique_edges,
		vector<vector<unsigned long> > const &edge_to_face_connectivity) {

	UpperEdge = 0;  LowerEdge = 0;
	LowerI    = 0;  UpperI    = 0;  MiddleI   = 0;
	LowerJ    = 0;  UpperJ    = 0; 	MiddleJ  = 0;
	rank = MASTER_NODE;

#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	clock_t build_start = clock();
	CurrentFace.resize(1, 0);
	Edge_To_Face_Connectivity = edge_to_face_connectivity;
	Unique_X_Bands = x_samples; //copy the x_values in
	//Sort the x_bands and make them unique
	sort(Unique_X_Bands.begin(), Unique_X_Bands.end());
	vector<su2double>::iterator iter;
	iter = unique(Unique_X_Bands.begin(), Unique_X_Bands.end());
	Unique_X_Bands.resize(distance(Unique_X_Bands.begin(), iter));
	X_Limits_of_Edges.resize(unique_edges.size(), vector<su2double>(2, 0));
	Y_Limits_of_Edges.resize(unique_edges.size(), vector<su2double>(2, 0));

	//Store the x and y values of each edge into a vector for a slight speed up as it
	//prevents some uncoalesced accesses
	for (unsigned int j = 0; j < unique_edges.size(); j++) {
		X_Limits_of_Edges[j][0] = x_samples[unique_edges[j][0]];
		X_Limits_of_Edges[j][1] = x_samples[unique_edges[j][1]];
		Y_Limits_of_Edges[j][0] = y_samples[unique_edges[j][0]];
		Y_Limits_of_Edges[j][1] = y_samples[unique_edges[j][1]];
	}

	//How many bands to search?
	unsigned int b_max = Unique_X_Bands.size() - 1;
	//Start with band 0, obviously
	unsigned int b = 0;
	//How many edges to check for intersection with the band?
	unsigned int e_max = unique_edges.size();
	//Start with edge indexes as 0.
	unsigned int i = 0;
	//Count the how many edges intersect a band
	unsigned int k = 0;
	//The high and low x value of each band
	su2double x_low = 0;
	su2double x_hi = 0;

	//Store the y_values of edges as required for searching
	Y_Values_of_Edge_Within_Band_And_Index.resize(Unique_X_Bands.size() - 1);

	//Check which edges intersect the band
	while (b < (b_max)) {
		x_low = Unique_X_Bands[b];
		x_hi = Unique_X_Bands[b + 1];
		i = 0;
		k = 0;
		//This while loop determined which edges appear in a paritcular band
		//The index of the edge being tested is 'i'
		while (i < e_max) {
			//Check if edge intersects the band (vertical edges are automatically discared)
			if (((X_Limits_of_Edges[i][0] <= x_low)
					and (X_Limits_of_Edges[i][1] >= x_hi))
					or ((X_Limits_of_Edges[i][1] <= x_low)
							and (X_Limits_of_Edges[i][0] >= x_hi))) {
				Y_Values_of_Edge_Within_Band_And_Index[b].push_back(make_pair(0.0, 0));
				//Save the edge index so it can latter be recalled (when searching)
				Y_Values_of_Edge_Within_Band_And_Index[b][k].second = i;
				//Determine what y value the edge takes in the middle of the band
				Y_Values_of_Edge_Within_Band_And_Index[b][k].first =
						Y_Limits_of_Edges[i][0]
																 + (Y_Limits_of_Edges[i][1] - Y_Limits_of_Edges[i][0])
																 / (X_Limits_of_Edges[i][1] - X_Limits_of_Edges[i][0])
																 * ((x_low + x_hi) / 2.0 - X_Limits_of_Edges[i][0]);
				//k counts the number of edges which have been found to intersect with the band
				k++;
			}
			//increment i, which  moves the algorithm along to the next edge
			i++;
		}
		//Sort the edges in the band depending on the y values they were found to have
		//It is worth noting that these y values are unique (i.e. edges cannot intersect in a band)
		sort(Y_Values_of_Edge_Within_Band_And_Index[b].begin(),
				Y_Values_of_Edge_Within_Band_And_Index[b].end());
		//Move on to the next band of x values
		b++;
	}
	//Initialize the search to table limits
	UpperI = Unique_X_Bands.size() - 1;
	LowerI = 0;

	su2double duration = ((su2double) clock() - (su2double) build_start)
					/ ((su2double) CLOCKS_PER_SEC);
	if (rank == MASTER_NODE)
		cout << duration << " seconds\n";
}

void CTrapezoidalMap::Search_Simplexes(su2double x, su2double y) {
	//Find the x band in which the current x value sits
	Search_Bands_For(x);
	//Within that band find edges between which the points rests
	//these two edges uniquely identify the containing polygon
	Search_Band_For_Edge(x, y);
	//Now identify the simplex from the edges (cannot be ambiguous
	//if all faces have the same number of edges). Cases where
	//ambiguity might occur are not expected to occur in thermotables
	vector<unsigned long> upper_edge_belongs_to_faces;
	vector<unsigned long> lower_edge_belongs_to_faces;
	upper_edge_belongs_to_faces = Edge_To_Face_Connectivity[UpperEdge];
	sort(upper_edge_belongs_to_faces.begin(), upper_edge_belongs_to_faces.end());
	lower_edge_belongs_to_faces = Edge_To_Face_Connectivity[LowerEdge];
	sort(lower_edge_belongs_to_faces.begin(), lower_edge_belongs_to_faces.end());
	//The intersection of the faces to which upper or lower belongs is
	//the face that both belong to.
	set_intersection(upper_edge_belongs_to_faces.begin(),
			upper_edge_belongs_to_faces.end(), lower_edge_belongs_to_faces.begin(),
			lower_edge_belongs_to_faces.end(), CurrentFace.begin());
}

void CTrapezoidalMap::Search_Bands_For(su2double x) {
	su2double x_middle, x_lower, x_upper;

	//TODO Add error when point is outside and break! The loop is not robust.
	do {
		MiddleI = (UpperI + LowerI) / 2;
		x_middle = Unique_X_Bands[MiddleI];
		x_lower = Unique_X_Bands[LowerI];
		x_upper = Unique_X_Bands[UpperI];
		//Step used for restarting the search on the low end
		if (x < x_lower and (LowerI > 0)) {
			UpperI = LowerI;
			LowerI = LowerI / 2;
			//Step used for restarting the search on the upper end
		} else if (x > x_upper and (UpperI < (Unique_X_Bands.size() - 1))) {
			LowerI = UpperI;
			UpperI = (UpperI + (Unique_X_Bands.size() - 1)) / 2;
			//After the restart is cleared, do the normal binary search
		} else if (x < x_middle) {
			UpperI = MiddleI;
		} else if (x > x_middle) {
			LowerI = MiddleI;
		} else if (x_middle == x) {
			LowerI = MiddleI;
			UpperI = LowerI + 1;
			break;
		}

	} while (UpperI - LowerI > 1);
}

void CTrapezoidalMap::Search_Band_For_Edge(su2double x, su2double y) {

	su2double RunVal, y00, y10, x00, x10;
	unsigned int RunEdge;
	UpperJ = Y_Values_of_Edge_Within_Band_And_Index[LowerI].size() - 1;
	LowerJ = 0;

	while (UpperJ - LowerJ > 1) {
		MiddleJ= (UpperJ + LowerJ) / 2;
		//Select the edge associated with the current x band (LowerI)
		//Search for the RunEdge in the MiddleJdirection (second value is index of edge)
		RunEdge = Y_Values_of_Edge_Within_Band_And_Index[LowerI][MiddleJ].second;
		y00 = Y_Limits_of_Edges[RunEdge][0];
		y10 = Y_Limits_of_Edges[RunEdge][1];
		x00 = X_Limits_of_Edges[RunEdge][0];
		x10 = X_Limits_of_Edges[RunEdge][1];
		//The search variable in j should be interpolated in i as well
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (x - x00);
		if (RunVal > y) {
			UpperJ = MiddleJ;
		} else if (RunVal < y) {
			LowerJ = MiddleJ;
		} else if (RunVal == y) {
			LowerJ = MiddleJ;
			UpperJ = LowerJ + 1;
			break;
		}

	}
	UpperEdge = Y_Values_of_Edge_Within_Band_And_Index[LowerI][UpperJ].second;
	LowerEdge = Y_Values_of_Edge_Within_Band_And_Index[LowerI][LowerJ].second;
}

CLookUpTable::CLookUpTable(CConfig *config, bool dimensional) :
				CFluidModel() {
	LUT_Debug_Mode = false;
	rank = MASTER_NODE;
	//TODO this has to be generalize for multi-zone
	unsigned int SinglePhaseZone = 0;
	CurrentZone= SinglePhaseZone;
	nInterpPoints = 3;
	CurrentPoints.resize(nInterpPoints, 0);
	LUT_Debug_Mode = config->GetLUT_Debug_Mode();
	Interpolation_Matrix.resize(nInterpPoints,
			vector<su2double>(nInterpPoints, 0));
	Interpolation_Matrix_Inverse.resize(nInterpPoints,
			vector<su2double>(nInterpPoints, 0));

#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	if (dimensional) {
		Pressure_Reference_Value = 1;
		Temperature_Reference_Value = 1;
		Density_Reference_Value = 1;
		Velocity_Reference_Value = 1;
		Energy_Reference_Value = 1;
	} else {
		Pressure_Reference_Value = config->GetPressure_Ref();
		Temperature_Reference_Value = config->GetTemperature_Ref();
		Density_Reference_Value = config->GetDensity_Ref();
		Velocity_Reference_Value = config->GetVelocity_Ref();
		Energy_Reference_Value = config->GetEnergy_Ref();
	}

	if ((config->GetLUTFileName()).find(".tec") != string::npos) {
		if (rank == MASTER_NODE) {
			cout << ".tec type LUT found" << endl;
		}
		LookUpTable_Load_TEC(config->GetLUTFileName());
	} else {
		if (rank == MASTER_NODE) {
			cout << "No recognized LUT format found, exiting!" << endl;
		}
		exit(EXIT_FAILURE);
	}

	if (rank == MASTER_NODE) {
		// Give the user some information on the size of the table
		cout << "Number of stations  in zone 0: " << nTable_Zone_Stations[0]
																																			<< endl;
		cout << "Number of triangles in zone 0: " << nTable_Zone_Triangles[0]
																																			 << endl;
		cout << "Number of stations  in zone 1: " << nTable_Zone_Stations[1]
																																			<< endl;
		cout << "Number of triangles in zone 1: " << nTable_Zone_Triangles[1]
																																			 << endl;
		cout
		<< "Detecting all unique edges and setting edge to face connectivity..."
		<< endl;
	}
	Get_Unique_Edges();
	if (rank == MASTER_NODE) {
		cout << "Number of edges in zone 0: " << Table_Zone_Edges[0].size() << endl;
		cout << "Number of edges in zone 1: " << Table_Zone_Edges[1].size() << endl;

	}

	if (rank == MASTER_NODE) {
		// Building an KD_tree for the HS thermopair
		cout << "Building trapezoidal map for rhoe..." << endl;
	}
	//Buld a map for all search pairs
	//Currently only zone 1 is actually in use so one could
	//also skip zone 0
	rhoe_map[0] = CTrapezoidalMap(ThermoTables_Density[0],
			ThermoTables_StaticEnergy[0], Table_Zone_Edges[0],
			Table_Edge_To_Face_Connectivity[0]);
	rhoe_map[1] = CTrapezoidalMap(ThermoTables_Density[1],
			ThermoTables_StaticEnergy[1], Table_Zone_Edges[1],
			Table_Edge_To_Face_Connectivity[1]);

	if (rank == MASTER_NODE) {
		cout << "Building trapezoidal map for Prho..." << endl;
	}
	Prho_map[0] = CTrapezoidalMap(ThermoTables_Pressure[0],
			ThermoTables_Density[0], Table_Zone_Edges[0],
			Table_Edge_To_Face_Connectivity[0]);
	Prho_map[1] = CTrapezoidalMap(ThermoTables_Pressure[1],
			ThermoTables_Density[1], Table_Zone_Edges[1],
			Table_Edge_To_Face_Connectivity[1]);

	if (rank == MASTER_NODE) {
		cout << "Building trapezoidal map for hs..." << endl;
	}
	hs_map[0] = CTrapezoidalMap(ThermoTables_Enthalpy[0], ThermoTables_Entropy[0],
			Table_Zone_Edges[0], Table_Edge_To_Face_Connectivity[0]);
	hs_map[1] = CTrapezoidalMap(ThermoTables_Enthalpy[1], ThermoTables_Entropy[1],
			Table_Zone_Edges[1], Table_Edge_To_Face_Connectivity[1]);

	if (rank == MASTER_NODE) {
		cout << "Building trapezoidal map for Ps..." << endl;
	}
	Ps_map[0] = CTrapezoidalMap(ThermoTables_Pressure[0], ThermoTables_Entropy[0],
			Table_Zone_Edges[0], Table_Edge_To_Face_Connectivity[0]);
	Ps_map[1] = CTrapezoidalMap(ThermoTables_Pressure[1], ThermoTables_Entropy[1],
			Table_Zone_Edges[1], Table_Edge_To_Face_Connectivity[1]);

	if (rank == MASTER_NODE) {
		cout << "Building trapezoidal map for rhoT..." << endl;
	}
	rhoT_map[0] = CTrapezoidalMap(ThermoTables_Density[0],
			ThermoTables_Temperature[0], Table_Zone_Edges[0],
			Table_Edge_To_Face_Connectivity[0]);
	;
	rhoT_map[1] = CTrapezoidalMap(ThermoTables_Density[1],
			ThermoTables_Temperature[1], Table_Zone_Edges[1],
			Table_Edge_To_Face_Connectivity[1]);
	;

	if (rank == MASTER_NODE) {
		cout << "Building trapezoidal map for PT (in vapor region only)..." << endl;
	}
	PT_map[0] = CTrapezoidalMap(ThermoTables_Pressure[0],
			ThermoTables_Temperature[0], Table_Zone_Edges[0],
			Table_Edge_To_Face_Connectivity[0]);
	;
	PT_map[1] = PT_map[0];

	if (rank == MASTER_NODE) {
		cout << "Print LUT errors? (LUT_Debug_Mode):  " << LUT_Debug_Mode << endl;
	}

	if (rank == MASTER_NODE) {
		cout << "Precomputing interpolation coefficients..." << endl;
	}
	Compute_Interpolation_Coefficients();
	if (rank == MASTER_NODE) {
		cout << "LuT fluid model ready for use" << endl;
	}

}

CLookUpTable::~CLookUpTable(void) {
	// Using vectors so no need to deallocate
	delete KD_tree;
}

void CLookUpTable::Get_Unique_Edges() {
	//Import all potential edges into a vector (assumes only triangles are used)
	//Run through both zones of the lutmesh.tec (2 zones assumed)
	for (unsigned int j = 0; j < 2; j++) {
		Table_Zone_Edges[j].resize(3 * nTable_Zone_Triangles[j], vector<unsigned long>(3, 0));
		//Fill with edges (based on triangulation
		//For each zone, go through all the triangles
		for (unsigned int i = 0; i < nTable_Zone_Triangles[j]; i++) {
			unsigned int smaller_point, larger_point;
			//Each triangle has 3 edges, add all of them to the list of edges for that zone
			//By using smaller_point and larger_point ensures that edges which occur in 2
			//triangles are recorded in the same way both times when they are found in the loop
			//Duplicates can be easily filtered out afterwards to yield unique edges only
			//The index i of the face the edge is associated with is stored as well
			//so that it can later be used for setting the edge-to-face connectivity
			smaller_point = Table_Zone_Triangles[j][i][0];
			larger_point = Table_Zone_Triangles[j][i][1];
			Table_Zone_Edges[j][3 * i + 0][0] = min(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 0][1] = max(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 0][2] = i;
			smaller_point = Table_Zone_Triangles[j][i][1];
			larger_point = Table_Zone_Triangles[j][i][2];
			Table_Zone_Edges[j][3 * i + 1][0] = min(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 1][1] = max(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 1][2] = i;
			smaller_point = Table_Zone_Triangles[j][i][2];
			larger_point = Table_Zone_Triangles[j][i][0];
			Table_Zone_Edges[j][3 * i + 2][0] = min(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 2][1] = max(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 2][2] = i;
		}
		//Sort the edges to enable selecting unique entries only
		stable_sort(Table_Zone_Edges[j].begin(), Table_Zone_Edges[j].end());
		//Set connectivities of the edges to the faces of the triangulation.
		//This is necessary for the trapezoidal map search i.e. the search
		//will identify 2 edges between which the point is located. Those two
		//edges should then be associated with the correct triangle so that
		//the interpolation coefficients for that particular face can be used

		//The index the edge will have in the final edge list
		unsigned int k_final = 0;
		//The index of the edge in the current edge list (non-unique but sorted)
		unsigned int k_temp = 0;
		//Traverse the current edge list
		while (k_temp < Table_Zone_Edges[j].size() - 1) {
			//Each unique edge is connected to at least 1 face so push_back the connectivity arr.
			Table_Edge_To_Face_Connectivity[j].push_back(vector<unsigned long>(1, -1));
			//Set the connectivty of the edge with a final index of k_final to the index of the
			//face that the temporary edge with index k_temp is associated to
			Table_Edge_To_Face_Connectivity[j][k_final][0] =
					Table_Zone_Edges[j][k_temp][2];
			//If the next edge in the temporary list is the same as the current edge
			//Then skip the next edge and increment k_temp by 2 ...
			if ((Table_Zone_Edges[j][k_temp][0] == Table_Zone_Edges[j][k_temp + 1][0])
					and (Table_Zone_Edges[j][k_temp][1]
																					 == Table_Zone_Edges[j][k_temp + 1][1])) {
				//...and add the face the k_temp+1 edge is associated with to
				//the list of faces that the k_final edge is connected to
				//(Only edges on the periphery of the thermotable will
				//be connected to only a single edge)
				Table_Edge_To_Face_Connectivity[j][k_final].push_back(
						Table_Zone_Edges[j][k_temp + 1][2]);
				k_temp++;//sic!
			}
			//Move on to the next temporary edge
			k_temp++;//sic!
			//Advance to the next final unique edge
			k_final++;
		}
		//The triangle index (entry 2 in vector) is no longer required as connectivities have
		//been set already. Removing the last entry enables unique edges to be found
		//using the "unique" algorithm for vectors
		for (unsigned int i = 0; i < nTable_Zone_Triangles[j]; i++) {
			Table_Zone_Edges[j][3 * i + 0].erase(
					Table_Zone_Edges[j][3 * i + 0].begin() + 2);
			Table_Zone_Edges[j][3 * i + 1].erase(
					Table_Zone_Edges[j][3 * i + 1].begin() + 2);
			Table_Zone_Edges[j][3 * i + 2].erase(
					Table_Zone_Edges[j][3 * i + 2].begin() + 2);
		}
		//Make edges unique
		vector<vector<unsigned long> >::iterator iter;
		iter = unique(Table_Zone_Edges[j].begin(), Table_Zone_Edges[j].end());
		Table_Zone_Edges[j].resize(distance(Table_Zone_Edges[j].begin(), iter));
	}
}
void CLookUpTable::Compute_Interpolation_Coefficients() {

	//First build a KD tree for the current zone in Prho
	PointIDs.resize(ThermoTables_Density[CurrentZone].size(), 0);
	coors.resize(2 * ThermoTables_Density[CurrentZone].size(), 0);

	for (unsigned long i = 0; i < ThermoTables_Density[CurrentZone].size(); i++) {
		PointIDs[i] = i;
		coors[2 * i] = ThermoTables_Density[CurrentZone][i];
		coors[2 * i + 1] = ThermoTables_Pressure[CurrentZone][i];
	}
	KD_tree = new su2_adtPointsOnlyClass(2, PointIDs.size(), coors.data(),
			PointIDs.data());
	query.resize(2, 0);
	best_dist.resize(nInterpPoints, 0);
	result_IDs.resize(nInterpPoints, 0);
	result_ranks.resize(nInterpPoints, 0);
	//Allocate the space for all the interpolation coefficients to be stored
	Rhoe_Interpolation_Matrix_Inverse[CurrentZone].resize(
			Table_Zone_Triangles[CurrentZone].size(),
			vector<vector<su2double> >(nInterpPoints,
					vector<su2double>(nInterpPoints, 0)));
	PT_Interpolation_Matrix_Inverse[CurrentZone].resize(
			Table_Zone_Triangles[CurrentZone].size(),
			vector<vector<su2double> >(nInterpPoints,
					vector<su2double>(nInterpPoints, 0)));
	Prho_Interpolation_Matrix_Inverse[CurrentZone].resize(
			Table_Zone_Triangles[CurrentZone].size(),
			vector<vector<su2double> >(nInterpPoints,
					vector<su2double>(nInterpPoints, 0)));
	rhoT_Interpolation_Matrix_Inverse[CurrentZone].resize(
			Table_Zone_Triangles[CurrentZone].size(),
			vector<vector<su2double> >(nInterpPoints,
					vector<su2double>(nInterpPoints, 0)));
	hs_Interpolation_Matrix_Inverse[CurrentZone].resize(
			Table_Zone_Triangles[CurrentZone].size(),
			vector<vector<su2double> >(nInterpPoints,
					vector<su2double>(nInterpPoints, 0)));
	Ps_Interpolation_Matrix_Inverse[CurrentZone].resize(
			Table_Zone_Triangles[CurrentZone].size(),
			vector<vector<su2double> >(nInterpPoints,
					vector<su2double>(nInterpPoints, 0)));
	//Also store the indexes of the points on which the coefficients are based
	//as these directly yueld funciton values
	Interpolation_Points[CurrentZone].resize(
			Table_Zone_Triangles[CurrentZone].size(),
			vector<unsigned long>(nInterpPoints, 0));

	//Now for each triangle in the zone calculate the e.g. 16 nearest points
	for (unsigned int i = 0; i < Table_Zone_Triangles[CurrentZone].size(); i++) {
		vector<unsigned long> Points_in_Triangle = Table_Zone_Triangles[CurrentZone][i];
		//The query point is to be the weighted average of the vertexes of the
		//triangle
		query[0] = 0;
		query[1] = 0;
		query[0] += ThermoTables_Density[CurrentZone][Points_in_Triangle[0]];
		query[0] += ThermoTables_Density[CurrentZone][Points_in_Triangle[1]];
		query[0] += ThermoTables_Density[CurrentZone][Points_in_Triangle[2]];
		query[0] /= 3;
		query[1] += ThermoTables_Pressure[CurrentZone][Points_in_Triangle[0]];
		query[1] += ThermoTables_Pressure[CurrentZone][Points_in_Triangle[1]];
		query[1] += ThermoTables_Pressure[CurrentZone][Points_in_Triangle[2]];
		query[1] /= 3;
		if (nInterpPoints>3){
			//If more than 3 points are required, than search the tree for the closet points
			KD_tree->Determine_N_NearestNodes(nInterpPoints, query.data(),
					best_dist.data(), result_IDs.data(), result_ranks.data());
			//Set the found points as the current points
			CurrentPoints = result_IDs;
			Interpolation_Points[CurrentZone][i] = result_IDs;
		}
		else if (nInterpPoints==3)
		{
			result_IDs = Table_Zone_Triangles[CurrentZone][i];
			CurrentPoints = result_IDs;
			Interpolation_Points[CurrentZone][i] = result_IDs;
		}

		//Now use the nearest 16 points to construct an interpolation function
		//for each search pair option
		Rhoe_Interpolation_Matrix_Inverse[CurrentZone][i] =
				Interpolation_Matrix_Prepare_And_Invert(ThermoTables_Density,
						ThermoTables_StaticEnergy);
		PT_Interpolation_Matrix_Inverse[CurrentZone][i] =
				Interpolation_Matrix_Prepare_And_Invert(ThermoTables_Pressure,
						ThermoTables_Temperature);
		Prho_Interpolation_Matrix_Inverse[CurrentZone][i] =
				Interpolation_Matrix_Prepare_And_Invert(ThermoTables_Pressure,
						ThermoTables_Density);
		rhoT_Interpolation_Matrix_Inverse[CurrentZone][i] =
				Interpolation_Matrix_Prepare_And_Invert(ThermoTables_Density,
						ThermoTables_Temperature);
		hs_Interpolation_Matrix_Inverse[CurrentZone][i] =
				Interpolation_Matrix_Prepare_And_Invert(ThermoTables_Enthalpy,
						ThermoTables_Entropy);
		Ps_Interpolation_Matrix_Inverse[CurrentZone][i] =
				Interpolation_Matrix_Prepare_And_Invert(ThermoTables_Pressure,
						ThermoTables_Entropy);
	}
}

void CLookUpTable::Get_Bounding_Simplex_From_TrapezoidalMap(
		CTrapezoidalMap *t_map, su2double x, su2double y) {

	t_map[CurrentZone].Search_Simplexes(x, y);
	CurrentFace = t_map[CurrentZone].getCurrentFace();
	CurrentPoints = Interpolation_Points[CurrentZone][CurrentFace];

}

void CLookUpTable::SetTDState_rhoe(su2double rho, su2double e) {

	Get_Bounding_Simplex_From_TrapezoidalMap(rhoe_map, rho, e);
	Interpolation_Matrix_Inverse =
			Rhoe_Interpolation_Matrix_Inverse[CurrentZone][CurrentFace];
	Calculate_Query_Specific_Coefficients(rho, e);

	vector<su2double> Weights;
	Weights = CalculateWeight(ThermoTables_Density, ThermoTables_StaticEnergy, rho, e);

	//Interpolate the fluid properties
	StaticEnergy = e;
	Density = rho;
	Entropy = Interpolate_Function2D(ThermoTables_Entropy, Weights);
	Pressure = Interpolate_Function2D(ThermoTables_Pressure, Weights);
	SoundSpeed2 = Interpolate_Function2D(ThermoTables_SoundSpeed2, Weights);
	Temperature = Interpolate_Function2D(ThermoTables_Temperature, Weights);
	dPdrho_e = Interpolate_Function2D(ThermoTables_dPdrho_e, Weights);
	dPde_rho = Interpolate_Function2D(ThermoTables_dPde_rho, Weights);
	dTdrho_e = Interpolate_Function2D(ThermoTables_dTdrho_e, Weights);
	dTde_rho = Interpolate_Function2D(ThermoTables_dTde_rho, Weights);
	Cp = Interpolate_Function2D(ThermoTables_Cp, Weights);

}

void CLookUpTable::SetTDState_PT(su2double P, su2double T) {

	Get_Bounding_Simplex_From_TrapezoidalMap(PT_map, P, T);
	Interpolation_Matrix_Inverse =
			PT_Interpolation_Matrix_Inverse[CurrentZone][CurrentFace];
	Calculate_Query_Specific_Coefficients(P, T);

	vector<su2double> Weights;
	Weights = CalculateWeight(ThermoTables_Pressure, ThermoTables_Temperature, P, T);

	//Interpolate the fluid properties
	Pressure = P;
	Density = Interpolate_Function2D(ThermoTables_Density, Weights);
	StaticEnergy = Interpolate_Function2D(ThermoTables_StaticEnergy, Weights);
	Entropy = Interpolate_Function2D(ThermoTables_Entropy, Weights);
	SoundSpeed2 = Interpolate_Function2D(ThermoTables_SoundSpeed2, Weights);
	Temperature = T;
	dPdrho_e = Interpolate_Function2D(ThermoTables_dPdrho_e, Weights);
	dPde_rho = Interpolate_Function2D(ThermoTables_dPde_rho, Weights);
	dTdrho_e = Interpolate_Function2D(ThermoTables_dTdrho_e, Weights);
	dTde_rho = Interpolate_Function2D(ThermoTables_dTde_rho, Weights);
	Cp = Interpolate_Function2D(ThermoTables_Cp, Weights);

}

void CLookUpTable::SetTDState_Prho(su2double P, su2double rho) {

	Get_Bounding_Simplex_From_TrapezoidalMap(Prho_map, P, rho);
	Interpolation_Matrix_Inverse =
			Prho_Interpolation_Matrix_Inverse[CurrentZone][CurrentFace];
	Calculate_Query_Specific_Coefficients(P, rho);

	vector<su2double> Weights;
	Weights = CalculateWeight(ThermoTables_Pressure, ThermoTables_Density, P, rho);

	//Interpolate the fluid properties
	Pressure = P;
	Density = rho;
	StaticEnergy = Interpolate_Function2D(ThermoTables_StaticEnergy, Weights);
	Entropy = Interpolate_Function2D(ThermoTables_Entropy, Weights);
	SoundSpeed2 = Interpolate_Function2D(ThermoTables_SoundSpeed2, Weights);
	Temperature = Interpolate_Function2D(ThermoTables_Temperature, Weights);
	dPdrho_e = Interpolate_Function2D(ThermoTables_dPdrho_e, Weights);
	dPde_rho = Interpolate_Function2D(ThermoTables_dPde_rho, Weights);
	dTdrho_e = Interpolate_Function2D(ThermoTables_dTdrho_e, Weights);
	dTde_rho = Interpolate_Function2D(ThermoTables_dTde_rho, Weights);
	Cp = Interpolate_Function2D(ThermoTables_Cp, Weights);
}

void CLookUpTable::SetEnergy_Prho(su2double P, su2double rho) {

	Get_Bounding_Simplex_From_TrapezoidalMap(Prho_map, P, rho);
	Interpolation_Matrix_Inverse =
			Prho_Interpolation_Matrix_Inverse[CurrentZone][CurrentFace];
	Calculate_Query_Specific_Coefficients(P, rho);

	vector<su2double> Weights;
	Weights = CalculateWeight(ThermoTables_Pressure, ThermoTables_Density, P, rho);

	StaticEnergy = Interpolate_Function2D(ThermoTables_StaticEnergy, Weights);
	Pressure = P;
	Density = rho;

}

void CLookUpTable::SetTDState_hs(su2double h, su2double s) {

	Get_Bounding_Simplex_From_TrapezoidalMap(hs_map, h, s);
	Interpolation_Matrix_Inverse =
			hs_Interpolation_Matrix_Inverse[CurrentZone][CurrentFace];
	Calculate_Query_Specific_Coefficients(h, s);

	vector<su2double> Weights;
	Weights = CalculateWeight(ThermoTables_Enthalpy, ThermoTables_Entropy, h, s);

	//Interpolate the fluid properties
	Entropy = s;
	StaticEnergy = Interpolate_Function2D(ThermoTables_StaticEnergy, Weights);
	Pressure = Interpolate_Function2D(ThermoTables_Pressure, Weights);
	Density = Interpolate_Function2D(ThermoTables_Density, Weights);
	SoundSpeed2 = Interpolate_Function2D(ThermoTables_SoundSpeed2, Weights);
	Temperature = Interpolate_Function2D(ThermoTables_Temperature, Weights);
	dPdrho_e = Interpolate_Function2D(ThermoTables_dPdrho_e, Weights);
	dPde_rho = Interpolate_Function2D(ThermoTables_dPde_rho, Weights);
	dTdrho_e = Interpolate_Function2D(ThermoTables_dTdrho_e, Weights);
	dTde_rho = Interpolate_Function2D(ThermoTables_dTde_rho, Weights);
	Cp = Interpolate_Function2D(ThermoTables_Cp, Weights);

}

void CLookUpTable::SetTDState_Ps(su2double P, su2double s) {

	Get_Bounding_Simplex_From_TrapezoidalMap(Ps_map, P, s);
	Interpolation_Matrix_Inverse =
			Ps_Interpolation_Matrix_Inverse[CurrentZone][CurrentFace];
	Calculate_Query_Specific_Coefficients(P, s);

	vector<su2double> Weights;
	Weights = CalculateWeight(ThermoTables_Pressure, ThermoTables_Entropy, P, s);

	//Interpolate the fluid properties
	Entropy = s;
	Pressure = P;
	StaticEnergy = Interpolate_Function2D(ThermoTables_StaticEnergy, Weights);
	Density = Interpolate_Function2D(ThermoTables_Density, Weights);
	SoundSpeed2 = Interpolate_Function2D(ThermoTables_SoundSpeed2, Weights);
	Temperature = Interpolate_Function2D(ThermoTables_Temperature, Weights);
	dPdrho_e = Interpolate_Function2D(ThermoTables_dPdrho_e, Weights);
	dPde_rho = Interpolate_Function2D(ThermoTables_dPde_rho, Weights);
	dTdrho_e = Interpolate_Function2D(ThermoTables_dTdrho_e, Weights);
	dTde_rho = Interpolate_Function2D(ThermoTables_dTde_rho, Weights);
	Cp = Interpolate_Function2D(ThermoTables_Cp, Weights);

}

void CLookUpTable::SetTDState_rhoT(su2double rho, su2double T) {

	Get_Bounding_Simplex_From_TrapezoidalMap(rhoT_map, rho, T);
	Interpolation_Matrix_Inverse =
			rhoT_Interpolation_Matrix_Inverse[CurrentZone][CurrentFace];
	Calculate_Query_Specific_Coefficients(rho, T);

	vector<su2double> Weights;
	Weights = CalculateWeight(ThermoTables_Density, ThermoTables_Temperature, rho, T);

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	StaticEnergy = Interpolate_Function2D(ThermoTables_StaticEnergy, Weights);
	Entropy = Interpolate_Function2D(ThermoTables_Entropy, Weights);
	Pressure = Interpolate_Function2D(ThermoTables_Pressure, Weights);
	SoundSpeed2 = Interpolate_Function2D(ThermoTables_SoundSpeed2, Weights);
	dPdrho_e = Interpolate_Function2D(ThermoTables_dPdrho_e, Weights);
	dPde_rho = Interpolate_Function2D(ThermoTables_dPde_rho, Weights);
	dTdrho_e = Interpolate_Function2D(ThermoTables_dTdrho_e, Weights);
	dTde_rho = Interpolate_Function2D(ThermoTables_dTde_rho, Weights);
	Cp = Interpolate_Function2D(ThermoTables_Cp, Weights);

}

inline void CLookUpTable::Gaussian_Inverse(unsigned int nDim) {
	//A temporary matrix to hold the inverse
	vector<vector<su2double> > temp;
	temp.resize(nDim, vector<su2double>(2 * nDim, 0));

	//Copy the desired matrix into the temporary matrix
	for (int i = 0; i < nDim; i++) {
		for (int j = 0; j < nDim; j++) {
			temp[i][j] = Interpolation_Matrix[i][j];
			temp[i][nDim + j] = 0;
		}
		temp[i][nDim + i] = 1;
	}

	su2double max_val;
	int max_idx;
	//Pivot each column such that the largest number possible divides the oter rows
	//The goal is to avoid zeros or small numbers in division.
	for (int k = 0; k < nDim - 1; k++) {
		max_idx = k;
		max_val = abs(temp[k][k]);
		//Find the largest value (pivot) in the column
		for (int j = k; j < nDim; j++) {
			if (abs(temp[j][k]) > max_val) {
				max_idx = j;
				max_val = abs(temp[j][k]);
			}
		}
		//Move the row with the highest value up
		for (int j = 0; j < (nDim * 2); j++) {
			su2double d = temp[k][j];
			temp[k][j] = temp[max_idx][j];
			temp[max_idx][j] = d;
		}
		//Subtract the moved row from all other rows
		for (int i = k + 1; i < nDim; i++) {
			su2double c = temp[i][k] / temp[k][k];
			for (int j = 0; j < (nDim * 2); j++) {
				temp[i][j] = temp[i][j] - temp[k][j] * c;
			}
		}
	}

	//Back-substitution
	for (int k = nDim - 1; k > 0; k--) {
		if (temp[k][k] != 0) {
			for (int i = k - 1; i > -1; i--) {
				su2double c = temp[i][k] / temp[k][k];
				for (int j = 0; j < (nDim * 2); j++) {
					temp[i][j] = temp[i][j] - temp[k][j] * c;
				}
			}
		}
	}
	//Normalize the inverse
	for (int i = 0; i < nDim; i++) {
		su2double c = temp[i][i];
		for (int j = 0; j < nDim; j++) {
			temp[i][j + nDim] = temp[i][j + nDim] / c;
		}
	}
	//Copy the inverse back to the main program flow
	for (int i = 0; i < nDim; i++) {
		for (int j = 0; j < nDim; j++) {
			Interpolation_Matrix_Inverse[i][j] = temp[i][j + nDim];
		}
	}
	return;
}

vector<su2double> CLookUpTable::Evaluate_Interpolation_Vector(su2double x,
		su2double y) {
	vector<su2double> interpolation_vector;
	interpolation_vector.resize(nInterpPoints, 0);
	interpolation_vector[0] = 1;
	interpolation_vector[1] = x;
	interpolation_vector[2] = y;
//	interpolation_vector[3] = x * y;
//	interpolation_vector[4] = x * y * y;
//	interpolation_vector[5] = x * y * y * y;
//	interpolation_vector[6] = x * x;
//	interpolation_vector[7] = x * x * y;
//	interpolation_vector[8] = x * x * y * y;
//	interpolation_vector[9] = x * x * y * y * y;
//	interpolation_vector[10] = y * y;
//	interpolation_vector[11] = x * x * x;
//	interpolation_vector[12] = x * x * x * y;
//	interpolation_vector[13] = x * x * x * y * y;
//	interpolation_vector[14] = x * x * x * y * y * y;
//	interpolation_vector[15] = y * y * y;

//  interpolation_vector[4] = x * x;
//	interpolation_vector[5] = y * y;
//	interpolation_vector[4] = log(y);
//	interpolation_vector[5] = log(x);
//	interpolation_vector[8] = log(x + y);

//	interpolation_vector[8] = exp(x);
//	interpolation_vector[9] = exp(y);
//	interpolation_vector[8] = log(y);
//	interpolation_vector[10] = exp(x);
//	interpolation_vector[11] = exp(y);
//	interpolation_vector[12] = exp(x * y);
//	interpolation_vector[13] = exp(-x);
//	interpolation_vector[14] = exp(-y);
//	interpolation_vector[15] = exp(x * y);

	return interpolation_vector;
}

vector<vector<su2double> > CLookUpTable::Interpolation_Matrix_Prepare_And_Invert(
		vector<su2double> *ThermoTables_X, vector<su2double> *ThermoTables_Y) {

	//Setup the LHM matrix for the interpolation
	for (int i = 0; i < nInterpPoints; i++) {
		su2double x = ThermoTables_X[CurrentZone][CurrentPoints[i]];
		su2double y = ThermoTables_Y[CurrentZone][CurrentPoints[i]];
		Interpolation_Matrix[i] = Evaluate_Interpolation_Vector(x, y);
	}

	//Invert the Interpolation matrix using Gaussian elimination with pivoting
	Gaussian_Inverse(nInterpPoints);
	su2double d;

	//Transpose the inverse
	for (int i = 0; i < (nInterpPoints - 1); i++) {
		for (int j = i + 1; j < nInterpPoints; j++) {
			d = Interpolation_Matrix_Inverse[i][j];
			Interpolation_Matrix_Inverse[i][j] = Interpolation_Matrix_Inverse[j][i];
			Interpolation_Matrix_Inverse[j][i] = d;
		}
	}
	return Interpolation_Matrix_Inverse;
}

void CLookUpTable::Calculate_Query_Specific_Coefficients(su2double x,
		su2double y) {
	vector<su2double> query_vector = Evaluate_Interpolation_Vector(x, y);
	Query_Specific_Interpolation_Coefficients.resize(nInterpPoints, 0);
	su2double d;
	for (int i = 0; i < nInterpPoints; i++) {
		d = 0;
		for (int j = 0; j < nInterpPoints; j++) {
			d = d + Interpolation_Matrix_Inverse[i][j] * query_vector[j];
		}
		Query_Specific_Interpolation_Coefficients[i] = d;
	}

}

vector<su2double> CLookUpTable::CalculateWeight(
		vector<su2double> *ThermoTables_X, vector<su2double> *ThermoTables_Y, su2double x, su2double y) {

  // Shepard interpolation exponent
  su2double p = 1;
  su2double x2, y2;
	vector<su2double>weight; weight.reserve(nInterpPoints);

	for (int i = 0; i < nInterpPoints; ++i) {
		x2 = pow((x - ThermoTables_X[CurrentZone][CurrentPoints[i]]),2);
		x2 /= x*x;
		y2 = pow((y - ThermoTables_Y[CurrentZone][CurrentPoints[i]]),2);
		y2 /= y*y;
		weight.push_back(1/pow(sqrt(x2+y2),p));
	}
	return weight;
}

su2double CLookUpTable::Interpolate_Function2D(
		vector<su2double> *ThermoTables_Z,vector<su2double> Weights ) {

	su2double result = 0.;
	su2double result_shepard = 0.;
  su2double weight_sum = 0.;
  su2double z = 0.;

	for (int i = 0; i < nInterpPoints; i++) {
		z = ThermoTables_Z[CurrentZone][CurrentPoints[i]];
		result += Query_Specific_Interpolation_Coefficients[i] * z;
	}

		return result;
}

void CLookUpTable::RecordState(char* file) {
	//Record the state of the fluid model to a file for
	//verificaiton purposes
	fstream fs;
	fs.open(file, fstream::app);
	fs.precision(17);
	assert(fs.is_open());
	fs << Temperature << ", ";
	fs << Density << ", ";
	fs << StaticEnergy << ", ";
	fs << Entropy << ", ";
	fs << Pressure << ", ";
	fs << SoundSpeed2 << ", ";
	fs << dPdrho_e << ", ";
	fs << dPde_rho << ", ";
	fs << dTdrho_e << ", ";
	fs << dTde_rho << ", ";
	fs << Cp << ", ";
	fs << Mu << ", ";
	fs << Kt << " ";
	fs << "\n";
	fs.close();
}

void CLookUpTable::LookUpTable_Print_To_File(char* filename) {
	//Print the current table zone a file such that the mesh can be plotted
	//externally (for verification purposes)
	int i = CurrentZone;
	for (int j = 0; j < nTable_Zone_Stations[i]; j++) {
		Temperature = ThermoTables_Temperature[i][j];
		Density = ThermoTables_Density[i][j];
		StaticEnergy = ThermoTables_StaticEnergy[i][j];
		Entropy = ThermoTables_Entropy[i][j];
		Pressure = ThermoTables_Pressure[i][j];
		SoundSpeed2 = ThermoTables_SoundSpeed2[i][j];
		dPdrho_e = ThermoTables_dPdrho_e[i][j];
		dPde_rho = ThermoTables_dPde_rho[i][j];
		dTdrho_e = ThermoTables_dTdrho_e[i][j];
		dTde_rho = ThermoTables_dTde_rho[i][j];
		Cp = ThermoTables_Cp[i][j];
		Kt = ThermoTables_Kt[i][j];
		Mu = ThermoTables_Mu[i][j];
		RecordState(filename);
	}
}

void CLookUpTable::LookUpTable_Load_TEC(std::string filename) {
	string line;
	string value;
	int found;
	unsigned int zone_scanned;

	ifstream table(filename.c_str());
	if (!table.is_open()) {
		if (rank == MASTER_NODE) {
			cout << "The LUT file appears to be missing!! " << filename << endl;
		}
		exit(EXIT_FAILURE);
	}
	zone_scanned = 0;
	//Go through all lines in the table file.
	getline(table, line);	//Skip the header
	while (getline(table, line)) {
		found = line.find("ZONE");
		if (found != -1) {
			if (rank == MASTER_NODE and LUT_Debug_Mode) {
				cout << line << endl;
			}
			istringstream in(line);
			//Note down the dimensions of the table
			int nPoints_in_Zone, nTriangles_in_Zone;
			string c1, c2, c3, c4;
			in >> c1 >> c2 >> nPoints_in_Zone >> c3 >> c4 >> nTriangles_in_Zone;
			//Create the actual LUT of CThermoLists which is used in the FluidModel
			nTable_Zone_Stations[zone_scanned] = nPoints_in_Zone;
			nTable_Zone_Triangles[zone_scanned] = nTriangles_in_Zone;
			//Allocate the memory for the table
			LookUpTable_Malloc(zone_scanned);

			//Load the values of the themordynamic properties at each table station
			for (int j = 0; j < nTable_Zone_Stations[zone_scanned]; j++) {
				getline(table, line);
				istringstream in(line);
				in >> ThermoTables_Density[zone_scanned][j];
				in >> ThermoTables_Pressure[zone_scanned][j];
				in >> ThermoTables_SoundSpeed2[zone_scanned][j];
				in >> ThermoTables_Cp[zone_scanned][j];
				in >> ThermoTables_Entropy[zone_scanned][j];
				in >> ThermoTables_Mu[zone_scanned][j];
				in >> ThermoTables_Kt[zone_scanned][j];
				in >> ThermoTables_dPdrho_e[zone_scanned][j];
				in >> ThermoTables_dPde_rho[zone_scanned][j];
				in >> ThermoTables_dTdrho_e[zone_scanned][j];
				in >> ThermoTables_dTde_rho[zone_scanned][j];
				in >> ThermoTables_Temperature[zone_scanned][j];
				in >> ThermoTables_StaticEnergy[zone_scanned][j];
				in >> ThermoTables_Enthalpy[zone_scanned][j];
			}
			//Skip empty line
			getline(table, line);
			//Load the triangles i.e. how the data point in each zone are connected
			for (int j = 0; j < nTable_Zone_Triangles[zone_scanned]; j++) {
				getline(table, line);
				istringstream in(line);
				in >> Table_Zone_Triangles[zone_scanned][j][0]
																										>> Table_Zone_Triangles[zone_scanned][j][1]
																																														 >> Table_Zone_Triangles[zone_scanned][j][2];
				//Triangles in .tec file are indexed from 1
				//In cpp it is more convenient to start with 0.
				Table_Zone_Triangles[zone_scanned][j][0]--;
				Table_Zone_Triangles[zone_scanned][j][1]--;
				Table_Zone_Triangles[zone_scanned][j][2]--;
			}

			zone_scanned++;
		}
	}

	table.close();
	//NonDimensionalise
	NonDimensionalise_Table_Values();
}

void CLookUpTable::LookUpTable_Malloc(unsigned int Index_of_Zone) {
	ThermoTables_StaticEnergy[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Entropy[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Enthalpy[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Density[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Pressure[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_SoundSpeed2[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Temperature[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_dPdrho_e[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_dPde_rho[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_dTdrho_e[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_dTde_rho[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Cp[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Mu[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Kt[Index_of_Zone] = vector<su2double>(
			nTable_Zone_Stations[Index_of_Zone], 0);
	Table_Zone_Triangles[Index_of_Zone] = vector<vector<unsigned long> >(
			nTable_Zone_Triangles[Index_of_Zone]);
	for (int j = 0; j < nTable_Zone_Triangles[Index_of_Zone]; j++) {
		Table_Zone_Triangles[Index_of_Zone][j] = vector<unsigned long>(3, 0);
	}
}

void CLookUpTable::NonDimensionalise_Table_Values() {
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nTable_Zone_Stations[i]; j++) {
			ThermoTables_Density[i][j] /= Density_Reference_Value;
			ThermoTables_Pressure[i][j] /= Pressure_Reference_Value;
			ThermoTables_SoundSpeed2[i][j] = pow(ThermoTables_SoundSpeed2[i][j], 2);
			ThermoTables_SoundSpeed2[i][j] /= pow(Velocity_Reference_Value, 2);
			ThermoTables_Cp[i][j] *= (Temperature_Reference_Value
					/ Energy_Reference_Value);
			ThermoTables_Entropy[i][j] *= (Temperature_Reference_Value
					/ Energy_Reference_Value);
			ThermoTables_dPdrho_e[i][j] *= (Density_Reference_Value
					/ Pressure_Reference_Value);
			ThermoTables_dPde_rho[i][j] *= (Energy_Reference_Value
					/ Pressure_Reference_Value);
			ThermoTables_dTdrho_e[i][j] *= (Density_Reference_Value
					/ Temperature_Reference_Value);
			ThermoTables_dTde_rho[i][j] *= (Energy_Reference_Value
					/ Temperature_Reference_Value);
			ThermoTables_Temperature[i][j] /= Temperature_Reference_Value;
			ThermoTables_StaticEnergy[i][j] /= Energy_Reference_Value;
			ThermoTables_Enthalpy[i][j] /= Energy_Reference_Value;
		}
	}
}
