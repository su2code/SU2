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

CThermoList::CThermoList() {

	StaticEnergy = 0.0;
	Entropy = 0.0;
	Enthalpy = 0.0;
	Density = 0.0;
	Pressure = 0.0;
	SoundSpeed2 = 0.0;
	Temperature = 0.0;
	dPdrho_e = 0.0;
	dPde_rho = 0.0;
	dTdrho_e = 0.0;
	dTde_rho = 0.0;
	Cp = 0.0;
	Mu = 0.0;
	dmudrho_T = 0.0;
	dmudT_rho = 0.0;
	Kt = 0.0;
	dktdrho_T = 0.0;
	dktdT_rho = 0.0;

}

CThermoList::~CThermoList() {

}

//void CThermoList::CTLprint()
//{
//	cout<<"StaticEnergy:"<<StaticEnergy<<endl;
//	cout<<"Enthalpy    :"<<Enthalpy<<endl;
//	cout<<"Entropy     :"<<Entropy<<endl;
//	cout<<"Density     :"<<Density<<endl;
//	cout<<"Pressure    :"<<Pressure<<endl;
//	cout<<"SoundSpeed2 :"<<SoundSpeed2<<endl;
//	cout<<"Temperature :"<<Temperature<<endl;
//	cout<<"dPdrho_e    :"<<dPdrho_e<<endl;
//	cout<<"dPde_rho    :"<<dPde_rho<<endl;
//	cout<<"dTdrho_e    :"<<dTdrho_e<<endl;
//	cout<<"dTde_rho    :"<<dTde_rho<<endl;
//	cout<<"Cp          :"<<Cp<<endl;
//	cout<<"Mu          :"<<Mu<<endl;
//	cout<<"dmudrho_T   :"<<dmudrho_T<<endl;
//	cout<<"dmudT_rho   :"<<dmudT_rho<<endl;
//	cout<<"Kt          :"<<Kt<<endl;
//	cout<<"dktdrho_T   :"<<dktdrho_T<<endl;
//	cout<<"dktdT_rho   :"<<dktdT_rho<<endl;
//}

CLookUpTable::CLookUpTable() :
				CFluidModel() {

	ThermoTables = NULL;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Interpolation_Coeff[i][j] = 0.0;
		}
	}
	HS_tree = NULL;
	iIndex = -1;
	jIndex = -1;
	Table_Pressure_Stations = 0;
	Table_Density_Stations = 0;
}

CLookUpTable::CLookUpTable(CConfig *config, bool dimensional) :
				CFluidModel() {
	ThermoTables = NULL;
	if (dimensional) {
		Pressure_Reference_Value = 1;
		Temperature_Reference_Value = 1;
		Density_Reference_Value = 1;
		Velocity_Reference_Value = 1;
	} else {
		Pressure_Reference_Value = config->GetPressure_Ref();
		Temperature_Reference_Value = config->GetTemperature_Ref();
		Density_Reference_Value = config->GetDensity_Ref();
		Velocity_Reference_Value = config->GetVelocity_Ref();
	}

	LookUpTable_Load_CFX(config->GetLUTFileName(), true);
	Remove_Two_Phase_Region_CFX_Table(true);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Interpolation_Coeff[i][j] = -1.0;
		}
	}
	// Initialize to a negative value to indicate the index is new (no restart)
	iIndex = -1;
	jIndex = -1;

	// Give the user some information on the size of the table
	cout << "Table_Pressure_Stations  : " << Table_Pressure_Stations << endl;
	cout << "Table_Density_Stations: " << Table_Density_Stations << endl;

	// Building an KD_tree for the HS thermopair
	cout << "Building HS_tree" << endl;
	su2double* xtemp = new su2double[Table_Density_Stations
																	 * Table_Pressure_Stations];
	su2double* ytemp = new su2double[Table_Density_Stations
																	 * Table_Pressure_Stations];
	int* itemp = new int[Table_Density_Stations * Table_Pressure_Stations];
	// Deep copy the x,y, and index values for the KD_tree
	for (int i = 0; i < Table_Density_Stations; i++) {
		for (int j = 0; j < Table_Pressure_Stations; j++) {
			xtemp[Table_Pressure_Stations * i + j] = ThermoTables[i][j].Enthalpy;
			ytemp[Table_Pressure_Stations * i + j] = ThermoTables[i][j].Entropy;
			itemp[Table_Pressure_Stations * i + j] = Table_Pressure_Stations * i + j;
		}
	}
	// Recursive KD_tree build starts here
	HS_tree = KD_Tree(xtemp, ytemp, itemp,
			Table_Pressure_Stations * Table_Density_Stations, 0);
	cout << "HS_tree built" << endl;
}

CLookUpTable::~CLookUpTable(void) {
	// Delete the lut
	for (int i = 0; i < Table_Density_Stations; i++) {
		delete[] ThermoTables[i];
	}
	delete ThermoTables;
	// Recursively delete all KD_node structs in the KD tree
	free_KD_tree(HS_tree);
}

void CLookUpTable::free_KD_tree(KD_node* root) {
	if (root->Branch_Dimension > 1) {
		// Descend into upper and lower branches of the root
		free_KD_tree(root->upper);
		free_KD_tree(root->lower);
	}
	delete[] root->x_values;
	delete[] root->y_values;
	delete[] root->Flattened_Point_Index;
	delete root;
}

struct KD_node* CLookUpTable::KD_Tree(su2double* x_values, su2double* y_values,
		int* Flattened_Point_Index, int dim, int depth) {

	struct KD_node *kdn = new KD_node;
	kdn->x_values = x_values;
	kdn->y_values = y_values;
	kdn->Flattened_Point_Index = Flattened_Point_Index;
	// The depth is used to define the splitting direction of the KD_tree
	kdn->Branch_Splitting_Direction = depth;
	// The dimension of the KD_tree branch
	kdn->Branch_Dimension = dim;
	// If  branch dimension is larger than 1, the branch get's split
	if (dim > 1) {
		/*!
		 * \brief In order to build the KD_tree the values must first be sorted.
		 * Bubblesort is used here because it's the fastest to implement.
		 * The custom implementation is necessary because the x_values, y_values, and i_values
		 * are to be sorted simultaneously.
		 */
		su2double temp;/*!< \brief Variable in which to temporarily store the variable value during a swap. */
		int itemp = 0;/*!< \brief Variable in which to temporarily store the index value during a swap. */
		int swaps = 0; /*!< \brief How many swaps have been performed during the sort. */
		int number_passes = 0;/*!< \brief 50% speed up bubblesort by realizing, the the n-th pass sorts the n-th element */
		bool sorted = false; /*!< \brief Triggers when the array is sorted. */
		// If the depth of current branch is even, the sort along x_values
		if (depth % 2 == 0) {
			while (not sorted) {
				swaps = 0;
				for (int i = 0; i < dim - 1 - number_passes; i++) {
					if (x_values[i] > x_values[i + 1]) {
						// The x_values determine the sorting of all three arrays.
						temp = x_values[i];
						x_values[i] = x_values[i + 1];
						x_values[i + 1] = temp;

						temp = y_values[i];
						y_values[i] = y_values[i + 1];
						y_values[i + 1] = temp;

						itemp = Flattened_Point_Index[i];
						Flattened_Point_Index[i] = Flattened_Point_Index[i + 1];
						Flattened_Point_Index[i + 1] = itemp;
						//Keep a record of the number of swaps performed
						swaps++;
					}
				}
				number_passes = number_passes + 1;
				// If no elements have been swapped in the Bubblesort, the sorting is done
				if (swaps == 0)
					sorted = true;
			}
			// If the depth of the branch is odd, then sort along the y_values--*/
		} else if (depth % 2 == 1) {
			while (not sorted) {
				swaps = 0;

				for (int i = 0; i < dim - 1 - number_passes; i++) {
					if (y_values[i] > y_values[i + 1]) {
						// The y_values determine the sorting of all three arrays.
						temp = y_values[i];
						y_values[i] = y_values[i + 1];
						y_values[i + 1] = temp;

						temp = x_values[i];
						x_values[i] = x_values[i + 1];
						x_values[i + 1] = temp;

						itemp = Flattened_Point_Index[i];
						Flattened_Point_Index[i] = Flattened_Point_Index[i + 1];
						Flattened_Point_Index[i + 1] = itemp;
						//Keep a record of the number of swaps performed
						swaps++;
					}
				}
				number_passes = number_passes + 1;
				// If no elements have been swapped in the Bubblesort, the sorting is done
				if (swaps == 0)
					sorted = true;
			}
		}
		/*!
		 * \brief Now that the arrays have been sorted they get split.
		 * The values lower than the median will go into the lower branch of the
		 * current tree, and the ones above into the upper branch. The arrays are
		 * split in half, so identical values are not accounted for in the median.
		 * Dynamic allocation is used during the recursive insertion of points into the tree.
		 */
		su2double* upperx = new su2double[dim / 2];
		su2double* uppery = new su2double[dim / 2];
		int* upperi = new int[dim / 2];
		su2double* lowerx = new su2double[dim - dim / 2];
		su2double* lowery = new su2double[dim - dim / 2];
		int* loweri = new int[dim - dim / 2];
		for (int i = dim / 2; i < dim; i++) {
			upperx[i - dim / 2] = x_values[i];
			uppery[i - dim / 2] = y_values[i];
			upperi[i - dim / 2] = Flattened_Point_Index[i];
		}
		for (int i = 0; i < dim / 2; i++) {
			lowerx[i] = x_values[i];
			lowery[i] = y_values[i];
			loweri[i] = Flattened_Point_Index[i];
		}
		/*!
		 * \brief Trigger the recursion into the upper and lower branches.
		 * The depth increment allows the tree to track whether to split the
		 * branch along x, or y.
		 */
		kdn->upper = KD_Tree(upperx, uppery, upperi, dim / 2, depth + 1);
		kdn->lower = KD_Tree(lowerx, lowery, loweri, dim - dim / 2, depth + 1);
	}
	return kdn;
}

su2double CLookUpTable::Dist2_KD_Tree(su2double x, su2double y,
		KD_node *branch) {
	su2double dist;
	/*!
	 * The distance between the branch and the search point is characterized
	 * by the distance between branch median point of the branch to the search.
	 */
	dist = pow((branch->x_values[branch->Branch_Dimension / 2] - x) / x, 2)\

			+ pow((branch->y_values[branch->Branch_Dimension / 2] - y) / y, 2);
	return dist;
}

void CLookUpTable::N_Nearest_Neighbours_KD_Tree(int N, su2double thermo1,
		su2double thermo2, KD_node *root, su2double *best_dist) {

	su2double dist = Dist2_KD_Tree(thermo1, thermo2, root);/*!< \brief First compute the Euclidean branch distance to the search item using  */
	/*!
	 * This algorithm is kept general such that it works for N nearest neighbors.
	 * The following loop look at the current point and checks whether it is closer
	 * than any of the N current closest points.
	 */
	int i = 0;
	while (i < N) {
		if (dist == best_dist[i])
			i = N + 1;
		if (dist < best_dist[i]) {
			for (int j = N - 1; j > i; j--) {
				best_dist[j] = best_dist[j - 1];
				Nearest_Neighbour_iIndex[j] = Nearest_Neighbour_iIndex[j - 1];
				Nearest_Neighbour_jIndex[j] = Nearest_Neighbour_jIndex[j - 1];
			}
			best_dist[i] = dist;
			Nearest_Neighbour_iIndex[i] =
					root->Flattened_Point_Index[root->Branch_Dimension / 2]
																			/ Table_Pressure_Stations;
			Nearest_Neighbour_jIndex[i] =
					root->Flattened_Point_Index[root->Branch_Dimension / 2]
																			% Table_Pressure_Stations;
			i = N + 1;
		}
		i++;

	}
	/*!
	 * Propagate the search further down the tree based on whether the search value
	 * is above or below the median value. If the branch splitting direction is
	 * even, the x_values are compared, if it is odd y_values get compared. This
	 * corresponds to the sorting order.
	 */
	if ((root->Branch_Dimension > 1)) {
		if (root->Branch_Splitting_Direction % 2 == 0) {
			if (root->x_values[root->Branch_Dimension / 2] <= thermo1) {
				/*!
				 * Propagate into the upper branch according to x_values
				 */
				N_Nearest_Neighbours_KD_Tree(N, thermo1, thermo2, root->upper,
						best_dist);
				if (dist < best_dist[N - 1]) {
					/*!
					 * Unwinding the search back up the tree ensures that the closest points
					 * are found even when their x and y values lie below the search term (i.e.
					 * this would mean the search would not cover them on its downward pass)
					 */
					N_Nearest_Neighbours_KD_Tree(N, thermo1, thermo2, root->lower,
							best_dist);
				}
			} else if (root->x_values[root->Branch_Dimension / 2] > thermo1) {
				/*!
				 * Propagate the search into the lower branch according to x_values
				 */
				N_Nearest_Neighbours_KD_Tree(N, thermo1, thermo2, root->lower,
						best_dist);
				if (dist < best_dist[N - 1]) {
					/*!
					 * Unwinding the search; see above.
					 */
					N_Nearest_Neighbours_KD_Tree(N, thermo1, thermo2, root->upper,
							best_dist);
				}
			}
			/*!
			 * If depth is odd, split the search in the y direction.
			 */
		} else if (root->Branch_Splitting_Direction % 2 == 1) {
			if (root->y_values[root->Branch_Dimension / 2] <= thermo2) {
				/*!
				 * Propagate the search into the upper branch according to y_values
				 */
				N_Nearest_Neighbours_KD_Tree(N, thermo1, thermo2, root->upper,
						best_dist);
				if (dist < best_dist[N - 1]) {
					/*!
					 * Unwinding the search; see above.
					 */
					N_Nearest_Neighbours_KD_Tree(N, thermo1, thermo2, root->lower,
							best_dist);
				}
			} else if (root->y_values[root->Branch_Dimension / 2] > thermo2) {
				/*!
				 * Propagate the search into the lower branch according to y_values
				 */
				N_Nearest_Neighbours_KD_Tree(N, thermo1, thermo2, root->lower,
						best_dist);
				/*!
				 * Unwinding the search; see above.
				 */
				if (dist < best_dist[N - 1]) {
					N_Nearest_Neighbours_KD_Tree(N, thermo1, thermo2, root->upper,
							best_dist);
				}
			}
		}
	}
}

void CLookUpTable::SetTDState_rhoe(su2double rho, su2double e) {
	// Check if inputs are in total range (necessary but not sufficient condition)
	if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
		//cerr << "RHOE Input Density out of bounds\n";
	}
	if ((e > StaticEnergy_Table_Limits[1])
			or (e < StaticEnergy_Table_Limits[0])) {
		//cerr << "RHOE Input StaticEnergy out of bounds\n";
	}
	// Linear interpolation requires 4 neighbors to be selected from the LUT
	Nearest_Neighbour_iIndex = new int[4];
	Nearest_Neighbour_jIndex = new int[4];
	su2double RunVal;
	unsigned int LowerI, UpperI, LowerJ, UpperJ, middleI, middleJ;
	su2double grad, x00, y00, y10, x10;

	// Starting values for the search
	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	// Bilinear search for the I index (density), not assuming rho is equispaced
	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;/*!< \brief Splitting index for the search */
		x00 = ThermoTables[middleI][LowerJ].Density;
		grad = ThermoTables[middleI + 1][LowerJ].Density - x00; /*!< \brief Density gradient with increasing i */
		if (x00 * grad > rho * grad) {
			UpperI = middleI;
		} else if (x00 < rho) {
			LowerI = middleI;
		} else if (x00 == rho) {
			LowerI = middleI;
			UpperI = LowerI + 1;
			break;
		}
	}

	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2; /*!< \brief Splitting index for the search */
		/*
		 * The variable names is composed of a (i,j) pair, which is used to denote which on is
		 * incremented.
		 */
		y00 = ThermoTables[LowerI][middleJ].StaticEnergy;
		y10 = ThermoTables[UpperI][middleJ].StaticEnergy;
		x00 = ThermoTables[LowerI][middleJ].Density;
		x10 = ThermoTables[UpperI][middleJ].Density;
		/*
		 * As StaticEnergy also depends on the i_index (not just search j),
		 * the StaticEnergy should be interpolated between y00 and y10 to determine
		 * whether to search the upper range of jIndexes or lower.
		 */
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (rho - x00);
		grad = ThermoTables[LowerI][middleJ + 1].StaticEnergy - y00;
		if (RunVal * grad > e * grad) {
			UpperJ = middleJ;
		} else if (RunVal * grad < e * grad) {
			LowerJ = middleJ;
		} else if (RunVal == e) {
			LowerJ = middleJ;
			UpperJ = LowerJ + 1;
			break;
		}
	}

	iIndex = LowerI;
	jIndex = LowerJ;
	/*
	 *Now use the quadrilateral which contains the point to interpolate
	 */
	su2double x, y;
	x = rho;
	y = e;
	//Set the nearest neigbours to the adjacent i and j vertexes
	Nearest_Neighbour_iIndex[0] = iIndex;
	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[0] = jIndex;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;
	//Determine the interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(x, y, "RHOE");

	//Interpolate the fluid properties
	StaticEnergy = e;
	Density = rho;
	Entropy = Interpolate_2D_Bilinear("Entropy");
	Pressure = Interpolate_2D_Bilinear("Pressure");
	//Enthalpy = Interpolate_2D_Bilinear("Enthalpy");
	SoundSpeed2 = Interpolate_2D_Bilinear("SoundSpeed2");
	Temperature = Interpolate_2D_Bilinear("Temperature");
	dPdrho_e = Interpolate_2D_Bilinear("dPdrho_e");
	dPde_rho = Interpolate_2D_Bilinear("dPde_rho");
	dTdrho_e = Interpolate_2D_Bilinear("dTdrho_e");
	dTde_rho = Interpolate_2D_Bilinear("dTde_rho");
	Cp = Interpolate_2D_Bilinear("Cp");
	Mu = Interpolate_2D_Bilinear("Mu");
	dmudrho_T = Interpolate_2D_Bilinear("dmudrho_T");
	dmudT_rho = Interpolate_2D_Bilinear("dmudT_rho");
	Kt = Interpolate_2D_Bilinear("Kt");
	dktdrho_T = Interpolate_2D_Bilinear("dktdrho_T");
	dktdT_rho = Interpolate_2D_Bilinear("dktdT_rho");

	//Check that the interpolated density and pressure are within LUT limits
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		//cerr << "RHOE Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		//cerr << "RHOE Interpolated Pressure out of bounds\n";
	}
	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;
}

void CLookUpTable::SetTDState_PT(su2double P, su2double T) {
	// Check if inputs are in total range (necessary but not sufficient condition)
	if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
		//cerr << "PT Input Pressure out of bounds\n";
	}
	if ((T > Temperature_Table_Limits[1]) or (T < Temperature_Table_Limits[0])) {
		//cerr << "PT Input Temperature out of bounds\n";
	}
	// Linear interpolation requires 4 neighbors to be selected from the LUT
	Nearest_Neighbour_iIndex = new int[4];
	Nearest_Neighbour_jIndex = new int[4];
	unsigned int LowerI, UpperI, LowerJ, UpperJ, middleI, middleJ;

	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	su2double grad, x00, y00, y01, x01, RunVal;
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	//Determine the J index using a binary search, and not assuming P is equispaced
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		x00 = ThermoTables[LowerI][middleJ].Pressure;
		grad = ThermoTables[LowerI][middleJ + 1].Pressure - x00;
		if (x00 * grad > P * grad) {
			UpperJ = middleJ;
		} else if (x00 < P) {
			LowerJ = middleJ;
		} else if (x00 == P) {
			LowerJ = middleJ;
			UpperJ = LowerJ + 1;
			break;
		}
	}

	//Determine the I index (for T)
	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		//Use interpolated T as the running variable for the search (RunVal)
		y00 = ThermoTables[middleI][LowerJ].Pressure;
		y01 = ThermoTables[middleI][UpperJ].Pressure;
		x00 = ThermoTables[middleI][LowerJ].Temperature;
		x01 = ThermoTables[middleI][UpperJ].Temperature;
		grad = ThermoTables[UpperI][LowerJ].Temperature - x00;
		RunVal = x00 + (x01 - x00) / (y01 - y00) * (P - y00);
		if (RunVal * grad > T * grad) {
			UpperI = middleI;
		} else if (RunVal * grad < T * grad) {
			LowerI = middleI;
		} else if (RunVal == T) {
			LowerI = middleI;
			UpperI = LowerI + 1;
			break;
		}
	}

	iIndex = LowerI;
	jIndex = LowerJ;

	su2double x, y;
	x = T;
	y = P;
	//Set the nearest neigbours to the adjacent i and j vertexes
	Nearest_Neighbour_iIndex[0] = iIndex;
	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[0] = jIndex;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;
	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(x, y, "PT");

	//Interpolate the fluid properties
	Temperature = T;
	Pressure = P;
	StaticEnergy = Interpolate_2D_Bilinear("StaticEnergy");
	//Enthalpy = Interpolate_2D_Bilinear("Enthalpy");
	Entropy = Interpolate_2D_Bilinear("Entropy");
	Density = Interpolate_2D_Bilinear("Density");
	SoundSpeed2 = Interpolate_2D_Bilinear("SoundSpeed2");
	dPdrho_e = Interpolate_2D_Bilinear("dPdrho_e");
	dPde_rho = Interpolate_2D_Bilinear("dPde_rho");
	dTdrho_e = Interpolate_2D_Bilinear("dTdrho_e");
	dTde_rho = Interpolate_2D_Bilinear("dTde_rho");
	Cp = Interpolate_2D_Bilinear("Cp");
	Mu = Interpolate_2D_Bilinear("Mu");
	dmudrho_T = Interpolate_2D_Bilinear("dmudrho_T");
	dmudT_rho = Interpolate_2D_Bilinear("dmudT_rho");
	Kt = Interpolate_2D_Bilinear("Kt");
	dktdrho_T = Interpolate_2D_Bilinear("dktdrho_T");
	dktdT_rho = Interpolate_2D_Bilinear("dktdT_rho");

	//Check that the interpolated density and pressure are within LUT limits
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		//cerr << "PT Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		//cerr << "PT Interpolated Pressure out of bounds\n";
	}
	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;
}

void CLookUpTable::SetTDState_Prho(su2double P, su2double rho) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
		//cerr << "PRHO Input Pressure out of bounds\n";
	}
	if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
		//cerr << "PRHO Input Density out of bounds\n";
	}
	// Linear interpolation requires 4 neighbors to be selected from the LUT
	Nearest_Neighbour_iIndex = new int[4];
	Nearest_Neighbour_jIndex = new int[4];
	unsigned int LowerI, UpperI, LowerJ, UpperJ, middleI, middleJ;

	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	su2double grad, x00, y00;
	//  Determine the I index with binary search (rho is not assumed equispaced)
	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		x00 = ThermoTables[middleI][LowerJ].Density;
		grad = ThermoTables[middleI + 1][LowerJ].Density - x00;
		if (x00 * grad > rho * grad) {
			UpperI = middleI;
		} else if (x00 < rho) {
			LowerI = middleI;
		} else if (x00 == rho) {
			LowerI = middleI;
			UpperI = LowerI + 1;
			break;
		}
	}

	//Binary search in the J index
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		y00 = ThermoTables[LowerI][middleJ].Pressure;
		grad = ThermoTables[LowerI][middleJ + 1].Pressure - y00;
		if (y00 * grad > P * grad) {
			UpperJ = middleJ;
		} else if (y00 < P) {
			LowerJ = middleJ;
		} else if (y00 == P) {
			LowerJ = middleJ;
			UpperJ = LowerJ + 1;
			break;
		}
	}

	iIndex = LowerI;
	jIndex = LowerJ;

	//Set the nearest neigbours to the adjacent i and j vertexes
	//for the bilinear interpolation
	su2double x, y;
	x = rho;
	y = P;
	Nearest_Neighbour_iIndex[0] = iIndex;
	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[0] = jIndex;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;
	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(x, y, "PRHO");

	//Interpolate the fluid properties
	Pressure = P;
	Density = rho;
	StaticEnergy = Interpolate_2D_Bilinear("StaticEnergy");
	//Enthalpy = Interpolate_2D_Bilinear("Enthalpy");
	Entropy = Interpolate_2D_Bilinear("Entropy");
	SoundSpeed2 = Interpolate_2D_Bilinear("SoundSpeed2");
	Temperature = Interpolate_2D_Bilinear("Temperature");
	dPdrho_e = Interpolate_2D_Bilinear("dPdrho_e");
	dPde_rho = Interpolate_2D_Bilinear("dPde_rho");
	dTdrho_e = Interpolate_2D_Bilinear("dTdrho_e");
	dTde_rho = Interpolate_2D_Bilinear("dTde_rho");
	Cp = Interpolate_2D_Bilinear("Cp");
	Mu = Interpolate_2D_Bilinear("Mu");
	dmudrho_T = Interpolate_2D_Bilinear("dmudrho_T");
	dmudT_rho = Interpolate_2D_Bilinear("dmudT_rho");
	Kt = Interpolate_2D_Bilinear("Kt");
	dktdrho_T = Interpolate_2D_Bilinear("dktdrho_T");
	dktdT_rho = Interpolate_2D_Bilinear("dktdT_rho");

	//Check that the interpolated density and pressure are within LUT limits
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		//cerr << "PRHO Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		//cerr << "PRHO Interpolated Pressure out of bounds\n";
	}
	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;

}

void CLookUpTable::SetEnergy_Prho(su2double P, su2double rho) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
		//cerr << "PRHO Input Pressure out of bounds\n";
	}
	if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
		//cerr << "PRHO Input Density out of bounds\n";
	}
	// Linear interpolation requires 4 neighbors to be selected from the LUT
	Nearest_Neighbour_iIndex = new int[4];
	Nearest_Neighbour_jIndex = new int[4];
	unsigned int LowerI, UpperI, LowerJ, UpperJ, middleI, middleJ;
	su2double grad, x00, y00;

	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	//Binary search for I index, not assuming it is equispaced.
	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		x00 = ThermoTables[middleI][LowerJ].Density;
		grad = ThermoTables[middleI + 1][LowerJ].Density - x00;
		if (x00 * grad > rho * grad) {
			UpperI = middleI;
		} else if (x00 < rho) {
			LowerI = middleI;
		} else if (x00 == rho) {
			LowerI = middleI;
			UpperI = LowerI + 1;
			break;
		}
	}

	//Determine the J index with pure binary search
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		//Check current value
		y00 = ThermoTables[LowerI][middleJ].Pressure;
		grad = ThermoTables[LowerI][middleJ + 1].Pressure - y00;
		if (y00 * grad > P * grad) {
			UpperJ = middleJ;
		} else if (y00 < P) {
			LowerJ = middleJ;
		} else if (y00 == P) {
			LowerJ = middleJ;
			UpperJ = LowerJ + 1;
			break;
		}
	}

	iIndex = LowerI;
	jIndex = LowerJ;

	su2double x, y;
	x = rho;
	y = P;
	//Set the nearest neigbours to the adjacent i and j vertexes
	Nearest_Neighbour_iIndex[0] = iIndex;
	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[0] = jIndex;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;
	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(x, y, "PRHO");

	Pressure = P;
	Density = rho;
	StaticEnergy = Interpolate_2D_Bilinear("StaticEnergy");

	//Check that the interpolated density and pressure are within LUT limits
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		//cerr << "PRHO Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		//cerr << "PRHO Interpolated Pressure out of bounds\n";
	}
	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;

}

void CLookUpTable::SetTDState_hs(su2double h, su2double s) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if ((h > Enthalpy_Table_Limits[1]) or (h < Enthalpy_Table_Limits[0])) {
		//cerr << "HS Input Enthalpy out of bounds\n";
	}
	if ((s > Entropy_Table_Limits[1]) or (s < Entropy_Table_Limits[0])) {
		//cerr << "HS Input Entropy out of bounds\n";
	}

	iIndex = HS_tree->Flattened_Point_Index[HS_tree->Branch_Dimension / 2]
																					/ Table_Pressure_Stations;
	jIndex = HS_tree->Flattened_Point_Index[HS_tree->Branch_Dimension / 2]
																					% Table_Pressure_Stations;
	int N = 4;
	// Linear interpolation requires 4 neighbors to be selected from the LUT
	// More points may be used for the inverse distance interpolation
	Nearest_Neighbour_iIndex = new int[N];
	Nearest_Neighbour_jIndex = new int[N];
	su2double *best_dist = new su2double[N];
	for (int i = 0; i < N; i++) {
		Nearest_Neighbour_iIndex[i] = -1;
		Nearest_Neighbour_jIndex[i] = -1;
	}

	//Preset the distance variables to something large, so they can be subsituted
	//by any point in the table.
	for (int i = 0; i < N; i++) {
		best_dist[i] = 1E10;
	}

	//Search the HS_tree for the thermo-pair values
	N_Nearest_Neighbours_KD_Tree(1, h, s, HS_tree, best_dist);

	//If an upper or right edge point is found, decrement it
	if (Nearest_Neighbour_iIndex[0] == (Table_Density_Stations - 1)) {
		Nearest_Neighbour_iIndex[0]--;
	}
	if (Nearest_Neighbour_jIndex[0] == (Table_Pressure_Stations - 1)) {
		Nearest_Neighbour_jIndex[0]--;
	}

	//Set the nearest neigbors to the adjacent i and j vertexes
	iIndex = Nearest_Neighbour_iIndex[0];
	jIndex = Nearest_Neighbour_jIndex[0];

	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;

	//Using the closest element found in the KD_tree as a starting point, now find the closest
	//quadrilateral containing the point using a simple zigzag search method
	su2double dx, dy, x00, y00, dx10, dx01, dx11, dy10, dy01, dy11;
	bool BOTTOM, TOP, LEFT, RIGHT, found = false; //check if bellow BOTTOM, above, TOP, left of LEFT, and right of RIGHT
	for (int k = 0; k < 20 and not found; k++) //20 is arbitrary and used primarily to avoid a while loop which could get stuck
	{
		Nearest_Neighbour_iIndex[0] = iIndex;
		Nearest_Neighbour_jIndex[0] = jIndex;
		Nearest_Neighbour_iIndex[1] = iIndex + 1;
		Nearest_Neighbour_iIndex[2] = iIndex;
		Nearest_Neighbour_iIndex[3] = iIndex + 1;
		Nearest_Neighbour_jIndex[1] = jIndex;
		Nearest_Neighbour_jIndex[2] = jIndex + 1;
		Nearest_Neighbour_jIndex[3] = jIndex + 1;

		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Enthalpy;
		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Entropy;
		dx = h - x00;
		dy = s - y00;
		dx01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Enthalpy
				- x00;
		dy01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Entropy
				- y00;
		dx10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Enthalpy
				- x00;
		dy10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Entropy
				- y00;
		dx11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Enthalpy
				- x00;
		dy11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Entropy
				- y00;

		BOTTOM = (dy * dx10) < (dx * dy10);
		TOP = ((dy - dy01) * (dx11 - dx01)) > ((dy11 - dy01) * (dx - dx01));
		RIGHT = ((dx - dx10) * (dy11 - dy10)) > ((dx11 - dx10) * (dy - dy10));
		LEFT = (dx * dy01) < (dx01 * dy);
		//Check BOTTOM quad boundary
		if (BOTTOM and !TOP) {
			if (jIndex != 0)
				jIndex--;
		}
		//Check RIGHT quad boundary
		else if (RIGHT and !LEFT) {
			if (iIndex != (Table_Density_Stations - 1))
				iIndex++;
		}
		//Check TOP quad boundary
		else if (TOP and !BOTTOM) {
			if (jIndex != (Table_Pressure_Stations - 1))
				jIndex++;
		}
		//Check LEFT quad boundary
		else if (LEFT and !RIGHT) {
			if (iIndex != 0)
				iIndex--;
		} else {
			found = true;
		}
	}
	//Determine interpolation coefficients
	su2double x = h;
	su2double y = s;
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(x, y, "HS");

	//Interpolate the fluid properties
	Entropy = s;
	//Enthalpy = h;
	StaticEnergy = Interpolate_2D_Bilinear("StaticEnergy");
	Pressure = Interpolate_2D_Bilinear("Pressure");
	Density = Interpolate_2D_Bilinear("Density");
	SoundSpeed2 = Interpolate_2D_Bilinear("SoundSpeed2");
	Temperature = Interpolate_2D_Bilinear("Temperature");
	dPdrho_e = Interpolate_2D_Bilinear("dPdrho_e");
	dPde_rho = Interpolate_2D_Bilinear("dPde_rho");
	dTdrho_e = Interpolate_2D_Bilinear("dTdrho_e");
	dTde_rho = Interpolate_2D_Bilinear("dTde_rho");
	Cp = Interpolate_2D_Bilinear("Cp");
	Mu = Interpolate_2D_Bilinear("Mu");
	dmudrho_T = Interpolate_2D_Bilinear("dmudrho_T");
	dmudT_rho = Interpolate_2D_Bilinear("dmudT_rho");
	Kt = Interpolate_2D_Bilinear("Kt");
	dktdrho_T = Interpolate_2D_Bilinear("dktdrho_T");
	dktdT_rho = Interpolate_2D_Bilinear("dktdT_rho");

	//Check that the interpolated density and pressure are within LUT limits
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		//cerr << "HS Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		//cerr << "HS Interpolated Pressure out of bounds\n";
	}
	delete[] best_dist;
	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;

}

void CLookUpTable::SetTDState_Ps(su2double P, su2double s) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
		//cerr << "PS Input Pressure out of bounds\n";
	}
	if ((s > Entropy_Table_Limits[1]) or (s < Entropy_Table_Limits[0])) {
		//cerr << "PS Input Entropy  out of bounds\n";
	}

	unsigned int LowerI, UpperI, LowerJ, UpperJ, middleI, middleJ;
	// Linear interpolation requires 4 neighbors to be selected from the LUT
	Nearest_Neighbour_iIndex = new int[4];
	Nearest_Neighbour_jIndex = new int[4];
	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	su2double grad, x00, y00, y01, x01, RunVal;

	//Determine the I index: rho is not equispaced
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		//Check current value
		y00 = ThermoTables[LowerI][middleJ].Pressure;
		grad = ThermoTables[LowerI][middleJ + 1].Pressure - y00;
		if (y00 * grad > P * grad) {
			UpperJ = middleJ;
		} else if (y00 < P) {
			LowerJ = middleJ;
		} else if (y00 == P) {
			LowerJ = middleJ;
			UpperJ = LowerJ + 1;
			break;
		}
	}

	//Determine the I index (for s)
	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		//Check current value
		y00 = ThermoTables[middleI][LowerJ].Pressure;
		y01 = ThermoTables[middleI][UpperJ].Pressure;
		x00 = ThermoTables[middleI][LowerJ].Entropy;
		x01 = ThermoTables[middleI][UpperJ].Entropy;
		grad = ThermoTables[UpperI][LowerJ].Entropy - x00;
		RunVal = x00 + (x01 - x00) / (y01 - y00) * (P - y00);
		grad = ThermoTables[middleI + 1][LowerJ].Entropy - x00;
		if (RunVal * grad > s * grad) {
			UpperI = middleI;
		} else if (RunVal * grad < s * grad) {
			LowerI = middleI;
		} else if (RunVal == s) {
			LowerI = middleI;
			UpperI = LowerI + 1;
			break;
		}
	}

	iIndex = LowerI;
	jIndex = LowerJ;

	su2double x, y;
	y = P;
	x = s;
	//Set the nearest neigbours to the adjacent i and j vertexes
	Nearest_Neighbour_iIndex[0] = iIndex;
	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[0] = jIndex;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;
	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(x, y, "PS");

	//Interpolate the fluid properties
	Entropy = s;
	Pressure = P;
	StaticEnergy = Interpolate_2D_Bilinear("StaticEnergy");
	//Enthalpy = Interpolate_2D_Bilinear("Enthalpy");
	Density = Interpolate_2D_Bilinear("Density");
	SoundSpeed2 = Interpolate_2D_Bilinear("SoundSpeed2");
	Temperature = Interpolate_2D_Bilinear("Temperature");
	dPdrho_e = Interpolate_2D_Bilinear("dPdrho_e");
	dPde_rho = Interpolate_2D_Bilinear("dPde_rho");
	dTdrho_e = Interpolate_2D_Bilinear("dTdrho_e");
	dTde_rho = Interpolate_2D_Bilinear("dTde_rho");
	Cp = Interpolate_2D_Bilinear("Cp");
	Mu = Interpolate_2D_Bilinear("Mu");
	dmudrho_T = Interpolate_2D_Bilinear("dmudrho_T");
	dmudT_rho = Interpolate_2D_Bilinear("dmudT_rho");
	Kt = Interpolate_2D_Bilinear("Kt");
	dktdrho_T = Interpolate_2D_Bilinear("dktdrho_T");
	dktdT_rho = Interpolate_2D_Bilinear("dktdT_rho");

	//Check that the interpolated density and pressure are within LUT limits
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		//cerr << "PS Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		//cerr << "PS Interpolated Pressure out of bounds\n";
	}
	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;
}

void CLookUpTable::SetTDState_rhoT(su2double rho, su2double T) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
		//cerr << "RHOT Input Density out of bounds\n";
	}
	if ((T > Temperature_Table_Limits[1]) or (T < Temperature_Table_Limits[0])) {
		//cerr << "RHOT Input Temperature out of bounds\n";
	}

	unsigned int LowerI, UpperI, LowerJ, UpperJ, middleI, middleJ;
	// Linear interpolation requires 4 neighbors to be selected from the LUT
	Nearest_Neighbour_iIndex = new int[4];
	Nearest_Neighbour_jIndex = new int[4];
	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	su2double grad, x00, y00, y10, x10, RunVal;

	//Determine the I index: rho is not equispaced
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		//Check current value
		x00 = ThermoTables[middleI][LowerJ].Density;
		grad = ThermoTables[middleI + 1][LowerJ].Density - x00;
		if (x00 * grad > rho * grad) {
			UpperI = middleI;
		} else if (x00 < rho) {
			LowerI = middleI;
		} else if (x00 == rho) {
			LowerI = middleI;
			UpperI = LowerI + 1;
			break;
		}
	}

	//Determine the J index (for T)
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		//Check current value
		y00 = ThermoTables[LowerI][middleJ].Temperature;
		y10 = ThermoTables[UpperI][middleJ].Temperature;
		x00 = ThermoTables[LowerI][middleJ].Density;
		x10 = ThermoTables[UpperI][middleJ].Density;
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (rho - x00);
		grad = ThermoTables[LowerI][UpperJ].Temperature - y00;
		if (RunVal * grad > T * grad) {
			UpperJ = middleJ;
		} else if (RunVal * grad < T * grad) {
			LowerJ = middleJ;
		} else if (RunVal == T) {
			LowerJ = middleJ;
			UpperJ = LowerJ + 1;
			break;
		}
	}

	iIndex = LowerI;
	jIndex = LowerJ;

	su2double x, y;
	x = rho;
	y = T;
	//Set the nearest neigbours to the adjacent i and j vertexes
	Nearest_Neighbour_iIndex[0] = iIndex;
	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[0] = jIndex;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;
	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(x, y, "RHOT");

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	StaticEnergy = Interpolate_2D_Bilinear("StaticEnergy");
	//Enthalpy = Interpolate_2D_Bilinear("Enthalpy");
	Entropy = Interpolate_2D_Bilinear("Entropy");
	Pressure = Interpolate_2D_Bilinear("Pressure");
	SoundSpeed2 = Interpolate_2D_Bilinear("SoundSpeed2");
	dPdrho_e = Interpolate_2D_Bilinear("dPdrho_e");
	dPde_rho = Interpolate_2D_Bilinear("dPde_rho");
	dTdrho_e = Interpolate_2D_Bilinear("dTdrho_e");
	dTde_rho = Interpolate_2D_Bilinear("dTde_rho");
	Cp = Interpolate_2D_Bilinear("Cp");
	Mu = Interpolate_2D_Bilinear("Mu");
	dmudrho_T = Interpolate_2D_Bilinear("dmudrho_T");
	dmudT_rho = Interpolate_2D_Bilinear("dmudT_rho");
	Kt = Interpolate_2D_Bilinear("Kt");
	dktdrho_T = Interpolate_2D_Bilinear("dktdrho_T");
	dktdT_rho = Interpolate_2D_Bilinear("dktdT_rho");

	//Check that the interpolated density and pressure are within LUT limits
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		//cerr << "RHOT Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		//cerr << "RHOT Interpolated Pressure out of bounds\n";
	}
	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;
}

su2double CLookUpTable::Interp2D_Inv_Dist(int N, std::string interpolant_var,
		su2double* dist) {
	su2double interp_result = 0;
	//The function values to interpolate from
	su2double *Interpolation_RHS = new su2double[N];
	//For each thermopair combination the values are filled differently
	if (interpolant_var == "StaticEnergy") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].StaticEnergy;
		}

	} else if (interpolant_var == "Entropy") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Entropy;
		}
	} else if (interpolant_var == "Density") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Density;
		}
	} else if (interpolant_var == "Pressure") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Pressure;
		}
	} else if (interpolant_var == "SoundSpeed2") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].SoundSpeed2;
		}
	} else if (interpolant_var == "Temperature") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Temperature;
		}
	} else if (interpolant_var == "dPdrho_e") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dPdrho_e;
		}

	} else if (interpolant_var == "dPde_rho") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dPde_rho;
		}
	} else if (interpolant_var == "dTdrho_e") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dTdrho_e;
		}
	} else if (interpolant_var == "dTde_rho") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dTde_rho;
		}
	} else if (interpolant_var == "Cp") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Cp;
		}
	} else if (interpolant_var == "Mu") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Mu;
		}

	} else if (interpolant_var == "dmudrho_T") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dmudrho_T;
		}

	} else if (interpolant_var == "dmudT_rho") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dmudT_rho;
		}

	} else if (interpolant_var == "Kt") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Kt;
		}

	} else if (interpolant_var == "dktdrho_T") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dktdrho_T;
		}

	} else if (interpolant_var == "dktdT_rho") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dktdT_rho;
		}
	} else if (interpolant_var == "Enthalpy") {
		for (int i = 0; i < N; i++) {
			Interpolation_RHS[i] =
					ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Enthalpy;
		}

	}

	//Sum of the weights times scores. and weights alone
	su2double dist_sum = 0;
	for (int i = 0; i < N; i++) {
		interp_result += (1 / dist[i]) * Interpolation_RHS[i];
		dist_sum += 1 / dist[i];
	}

	interp_result = interp_result / dist_sum;
	delete[] Interpolation_RHS;
	return interp_result;
}
inline void CLookUpTable::Gaussian_Inverse(int nDim) {
	//A temporary matrix to hold the inverse, dynamically allocated
	su2double **temp = new su2double*[nDim];
	for (int i = 0; i < nDim; i++) {
		temp[i] = new su2double[2 * nDim];
	}

	//Copy the desired matrix into the temportary matrix
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
		max_val = temp[k][k];
		//Find the largest value in the column
		for (int j = k; j < nDim; j++) {
			if (abs(temp[j][k]) > max_val) {
				max_idx = j;
				max_val = temp[j][k];
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
			Interpolation_Coeff[i][j] = temp[i][j + nDim];
		}
	}
	//Delete dynamic template
	delete[] temp;
	return;
}
void CLookUpTable::Interpolate_2D_Bilinear_Read_Coordinates(su2double *coords,
		std::string grid_var) {
	su2double x00, y00, x10, x01, x11, y10, y01, y11;

	//Load in the coordinates of the qudrilateral
	if (grid_var == "RHOE") {
		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Density;

		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].StaticEnergy;

		x01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Density;

		y01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].StaticEnergy;

		x10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Density;

		y10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].StaticEnergy;

		x11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Density;

		y11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].StaticEnergy;

	} else if (grid_var == "PT") {
		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Pressure;

		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Temperature;

		y01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Pressure;

		x01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Temperature;

		y10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Pressure;

		x10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Temperature;

		y11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Pressure;

		x11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Temperature;

	} else if (grid_var == "PRHO") {
		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Pressure;

		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Density;

		y01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Pressure;

		x01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Density;

		y10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Pressure;

		x10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Density;

		y11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Pressure;

		x11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Density;

	} else if (grid_var == "RHOT") {
		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Density;

		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Temperature;

		x01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Density;

		y01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Temperature;

		x10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Density;

		y10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Temperature;

		x11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Density;

		y11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Temperature;

	} else if (grid_var == "PS") {
		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Pressure;

		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Entropy;

		y01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Pressure;

		x01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Entropy;

		y10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Pressure;

		x10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Entropy;

		y11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Pressure;

		x11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Entropy;

	} else if (grid_var == "HS") {
		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Enthalpy;

		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Entropy;

		x01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Enthalpy;

		y01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Entropy;

		x10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Enthalpy;

		y10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Entropy;

		x11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Enthalpy;

		y11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Entropy;

	}
	coords[0] = x00;
	coords[1] = y00;
	coords[2] = x10;
	coords[3] = y10;
	coords[4] = x01;
	coords[5] = y01;
	coords[6] = x11;
	coords[7] = y11;
	return;

}

void CLookUpTable::Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(su2double x,
		su2double y, std::string grid_var) {
	//The x,y coordinates of the quadrilateral
	su2double x00, y00, x10, x01, x11, y10, y01, y11;
	su2double coords[8];

	//Load in the coordinates of the qudrilateral
	Interpolate_2D_Bilinear_Read_Coordinates(coords, grid_var);
	x00 = coords[0];
	y00 = coords[1];
	x10 = coords[2];
	y10 = coords[3];
	x01 = coords[4];
	y01 = coords[5];
	x11 = coords[6];
	y11 = coords[7];

	//Check if x, y is indeed in the quad
	//The (true and not false) type of logic is needed as the both monotonically
	//increasing and monotonically decreasing functions need to pass the same test
	bool BOTTOM, TOP, LEFT, RIGHT, OUT_OF_BOUNDS;
	su2double dy, dx, dx10, dy10, dx01, dy01, dx11, dy11;
	dx = x - x00;
	dy = y - y00;
	dx10 = x10 - x00;
	dy10 = y10 - y00;
	dx01 = x01 - x00;
	dy01 = y01 - y00;
	dx11 = x11 - x00;
	dy11 = y11 - y00;
	BOTTOM = (dy * dx10) < (dx * dy10);
	TOP = ((dy - dy01) * (dx11 - dx01)) > ((dy11 - dy01) * (dx - dx01));
	RIGHT = ((dx - dx10) * (dy11 - dy10)) > ((dx11 - dx10) * (dy - dy10));
	LEFT = (dx * dy01) < (dx01 * dy);
	OUT_OF_BOUNDS = false;
	//Check BOTTOM quad boundary
	if (BOTTOM and !TOP) {
		OUT_OF_BOUNDS = true;
		//Check if the point is also outside the LUT
		if (Nearest_Neighbour_jIndex[0] == 0) {
			//Change the the 1 by 1 quad into a larger one for less extreme extrapolation
			Nearest_Neighbour_jIndex[1] = Table_Pressure_Stations;
			Nearest_Neighbour_jIndex[3] = Table_Pressure_Stations;
//			Nearest_Neighbour_jIndex[1] = Nearest_Neighbour_jIndex[1]+1 ;
//			Nearest_Neighbour_jIndex[3] = Nearest_Neighbour_jIndex[3]+1	 ;
			//cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
			//<< Nearest_Neighbour_jIndex[0]
			//<< " interpolation point lies below the LUT\n";
		} else {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies below bottom boundary of selected quad\n";
		}
	}
	//Check RIGHT quad boundary
	if (RIGHT and !LEFT) {
		OUT_OF_BOUNDS = true;
		//Check if the point is also outside the LUT
		if (Nearest_Neighbour_iIndex[0] == (Table_Density_Stations - 2)) {
			//Change the the 1 by 1 quad into a larger one for less extreme extrapolation
			Nearest_Neighbour_iIndex[0] = 0;// Table_Density_Stations;
			Nearest_Neighbour_iIndex[2] = 0;//Table_Density_Stations;
//			Nearest_Neighbour_iIndex[0] = Nearest_Neighbour_iIndex[0]-1;
//			Nearest_Neighbour_iIndex[2] = Nearest_Neighbour_iIndex[0]-1;
			//cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
			//<< Nearest_Neighbour_jIndex[0]
			//<< " interpolation point lies right of the LUT\n";
		} else {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies to the right of the boundary of selected quad\n";
		}
	}
	//Check TOP quad boundary
	if (TOP and !BOTTOM) {
		OUT_OF_BOUNDS = true;
		//Check if the point is also outside the LUT
		if (Nearest_Neighbour_jIndex[0] == (Table_Pressure_Stations - 2)) {
			//Change the the 1 by 1 quad into a larger one for less extreme extrapolation
			Nearest_Neighbour_jIndex[0] = 0;//Table_Pressure_Stations;
			Nearest_Neighbour_jIndex[2] = 0;//Table_Pressure_Stations;
//			Nearest_Neighbour_jIndex[0] = Nearest_Neighbour_jIndex[0]-1 ;
//			Nearest_Neighbour_jIndex[2] = Nearest_Neighbour_jIndex[2]-1 ;
			//cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
			//	<< Nearest_Neighbour_jIndex[0]
			//<< " interpolation point lies above the LUT\n";
		} else {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< +" interpolation point lies above the boundary of selected quad\n";
		}
	}
	//Check LEFT quad boundary
	if (LEFT and !RIGHT) {
		OUT_OF_BOUNDS = true;
		//Check if the point is also outside the LUT
		if (Nearest_Neighbour_iIndex[0] == 0) {
			//Change the the 1 by 1 quad into a larger one for less extreme extrapolation
			Nearest_Neighbour_iIndex[1] = Table_Density_Stations;
			Nearest_Neighbour_iIndex[3] = Table_Density_Stations;
//			Nearest_Neighbour_iIndex[1] = Nearest_Neighbour_iIndex[1]+1 ;
//			Nearest_Neighbour_iIndex[3] = Nearest_Neighbour_iIndex[3]+1 ;
			//cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
			//<< Nearest_Neighbour_jIndex[0]
			//<< " interpolation point lies left of the LUT\n";
		} else {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< +" interpolation point lies to the left of the boundary of selected quad\n";
		}
	}
	if (OUT_OF_BOUNDS) {
		Interpolate_2D_Bilinear_Read_Coordinates(coords, grid_var);
		x00 = coords[0];
		y00 = coords[1];
		x10 = coords[2];
		y10 = coords[3];
		x01 = coords[4];
		y01 = coords[5];
		x11 = coords[6];
		y11 = coords[7];
	}

	//Setup the LHM matrix for the interpolation (Vandermonde)
	Interpolation_Matrix[0][0] = 1;
	Interpolation_Matrix[0][1] = 0;
	Interpolation_Matrix[0][2] = 0;
	Interpolation_Matrix[0][3] = 0;

	Interpolation_Matrix[1][0] = 1;
	Interpolation_Matrix[1][1] = x10 - x00;
	Interpolation_Matrix[1][2] = y10 - y00;
	Interpolation_Matrix[1][3] = (x10 - x00) * (y10 - y00);

	Interpolation_Matrix[2][0] = 1;
	Interpolation_Matrix[2][1] = x01 - x00;
	Interpolation_Matrix[2][2] = y01 - y00;
	Interpolation_Matrix[2][3] = (x01 - x00) * (y01 - y00);

	Interpolation_Matrix[3][0] = 1;
	Interpolation_Matrix[3][1] = x11 - x00;
	Interpolation_Matrix[3][2] = y11 - y00;
	Interpolation_Matrix[3][3] = (x11 - x00) * (y11 - y00);

	//Invert the Interpolation matrix using Gaussian elimination with pivoting
	Gaussian_Inverse(4);
	su2double d;

	//Transpose the inverse
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 4; j++) {
			d = Interpolation_Coeff[i][j];
			Interpolation_Coeff[i][j] = Interpolation_Coeff[j][i];
			Interpolation_Coeff[j][i] = d;
		}
	}
	//The transpose allows the same coefficients to be used
	// for all Thermo variables (need only 4 coefficients)
	for (int i = 0; i < 4; i++) {
		d = 0;
		d = d + Interpolation_Coeff[i][0] * 1;
		d = d + Interpolation_Coeff[i][1] * (x - x00);
		d = d + Interpolation_Coeff[i][2] * (y - y00);
		d = d + Interpolation_Coeff[i][3] * (x - x00) * (y - y00);
		Interpolation_Coeff[i][0] = d;
	}

	return;
}

su2double CLookUpTable::Interpolate_2D_Bilinear(string interpolant_var) {
	//The function values at the 4 corners of the quad
	su2double func_value_at_i0j0, func_value_at_i1j0, func_value_at_i0j1,
	func_value_at_i1j1;
	if (interpolant_var == "StaticEnergy") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].StaticEnergy;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].StaticEnergy;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].StaticEnergy;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].StaticEnergy;
	} else if (interpolant_var == "Entropy") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Entropy;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Entropy;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Entropy;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Entropy;
	} else if (interpolant_var == "Density") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Density;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Density;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Density;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Density;
	} else if (interpolant_var == "Pressure") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Pressure;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Pressure;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Pressure;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Pressure;
	} else if (interpolant_var == "SoundSpeed2") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].SoundSpeed2;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].SoundSpeed2;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].SoundSpeed2;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].SoundSpeed2;
	} else if (interpolant_var == "Temperature") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Temperature;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Temperature;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Temperature;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Temperature;
	} else if (interpolant_var == "dPdrho_e") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dPdrho_e;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dPdrho_e;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dPdrho_e;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dPdrho_e;
	} else if (interpolant_var == "dPde_rho") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dPde_rho;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dPde_rho;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dPde_rho;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dPde_rho;
	} else if (interpolant_var == "dTdrho_e") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dTdrho_e;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dTdrho_e;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dTdrho_e;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dTdrho_e;
	} else if (interpolant_var == "dTde_rho") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dTde_rho;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dTde_rho;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dTde_rho;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dTde_rho;
	} else if (interpolant_var == "Cp") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Cp;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Cp;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Cp;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Cp;
	} else if (interpolant_var == "Mu") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Mu;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Mu;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Mu;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Mu;
	} else if (interpolant_var == "dmudrho_T") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dmudrho_T;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dmudrho_T;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dmudrho_T;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dmudrho_T;
	} else if (interpolant_var == "dmudT_rho") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dmudT_rho;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dmudT_rho;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dmudT_rho;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dmudT_rho;
	} else if (interpolant_var == "Kt") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Kt;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Kt;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Kt;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Kt;
	} else if (interpolant_var == "dktdrho_T") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dktdrho_T;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dktdrho_T;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dktdrho_T;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dktdrho_T;
	} else if (interpolant_var == "dktdT_rho") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dktdT_rho;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dktdT_rho;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dktdT_rho;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dktdT_rho;
	} else if (interpolant_var == "Enthalpy") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Enthalpy;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Enthalpy;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Enthalpy;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Enthalpy;
	}

	//The Interpolation_Coeff values depend on location alone
	//and are the same regardless of function values
	su2double result = 0;
	result = result + Interpolation_Coeff[0][0] * func_value_at_i0j0;
	result = result + Interpolation_Coeff[1][0] * func_value_at_i1j0;
	result = result + Interpolation_Coeff[2][0] * func_value_at_i0j1;
	result = result + Interpolation_Coeff[3][0] * func_value_at_i1j1;

	return result;
}

void CLookUpTable::RecordState(char* file) {
	//Record the state of the fluid model to a file for
	//verificaiton purposes
	//	fstream fs;
	//	fs.open(file, fstream::app);
	//	fs.precision(17);
	//	assert(fs.is_open());
	//	fs << Temperature << ", ";
	//	fs << Density << ", ";
	//	fs << Enthalpy << ", ";
	//	fs << StaticEnergy << ", ";
	//	fs << Entropy << ", ";
	//	fs << Pressure << ", ";
	//	fs << SoundSpeed2 << ", ";
	//	fs << dPdrho_e << ", ";
	//	fs << dPde_rho << ", ";
	//	fs << dTdrho_e << ", ";
	//	fs << dTde_rho << ", ";
	//	fs << Cp << ", ";
	//	fs << Mu << ", ";
	//	fs << dmudrho_T << ", ";
	//	fs << dmudT_rho << ", ";
	//	fs << Kt << ", ";
	//	fs << dktdrho_T << ", ";
	//	fs << dktdT_rho << ", ";
	//	fs << "\n";
	//	fs.close();
}

void CLookUpTable::LookUpTable_Print_To_File(char* filename) {
	//Print the entire table to a file such that the mesh can be plotted
	//externally (for verification purposes)
	//	for (int i = 0; i < Table_Density_Stations; i++) {
	//		for (int j = 0; j < Table_Pressure_Stations; j++) {
	//			iIndex = i;
	//			jIndex = j;
	//			Temperature = ThermoTables[iIndex][jIndex].Temperature;
	//			Density = ThermoTables[iIndex][jIndex].Density;
	//			Enthalpy = ThermoTables[iIndex][jIndex].Enthalpy;
	//			StaticEnergy = ThermoTables[iIndex][jIndex].StaticEnergy;
	//			Entropy = ThermoTables[iIndex][jIndex].Entropy;
	//			Pressure = ThermoTables[iIndex][jIndex].Pressure;
	//			SoundSpeed2 = ThermoTables[iIndex][jIndex].SoundSpeed2;
	//			dPdrho_e = ThermoTables[iIndex][jIndex].dPdrho_e;
	//			dPde_rho = ThermoTables[iIndex][jIndex].dPde_rho;
	//			dTdrho_e = ThermoTables[iIndex][jIndex].dTdrho_e;
	//			dTde_rho = ThermoTables[iIndex][jIndex].dTde_rho;
	//			Cp = ThermoTables[iIndex][jIndex].Cp;
	//			Mu = ThermoTables[iIndex][jIndex].Mu;
	//			dmudrho_T = ThermoTables[iIndex][jIndex].dmudrho_T;
	//			dmudT_rho = ThermoTables[iIndex][jIndex].dmudT_rho;
	//			Kt = ThermoTables[iIndex][jIndex].Kt;
	//			dktdrho_T = ThermoTables[iIndex][jIndex].dktdrho_T;
	//			dktdT_rho = ThermoTables[iIndex][jIndex].dktdT_rho;
	//			RecordState(filename);
	//		}
	//	}
	//
}

void CLookUpTable::Remove_Two_Phase_Region_CFX_Table(bool is_not_two_phase) {
	int** Indexes_of_two_phase = new int*[Table_Density_Stations];

	for (int i = 0; i < Table_Density_Stations; i++) {
		Indexes_of_two_phase[i] = new int[Table_Pressure_Stations];
	}
	for (int i = 0; i < Table_Density_Stations; i++) {
		for (int j = 0; j < Table_Pressure_Stations; j++) {
			Indexes_of_two_phase[i][j] = 0;
		}
	}
	//	//Edge detection going down
	//	for (int j = 0; j < Table_Pressure_Stations; j++) {
	//		for (int i = 0; i < Table_Density_Stations - 1; i++) {
	//			if (abs(
	//					ThermoTables[i + 1][j].Enthalpy - ThermoTables[i][j].Enthalpy)
	//					> 0.1 * ThermoTables[i + 1][j].Enthalpy) {
	//				Indexes_of_two_phase[i+1][j] = -10;
	//			}
	//		}
	//	}
	//	//Edge detection going up
	//	for (int j = 0; j < Table_Pressure_Stations; j++) {
	//		for (int i = Table_Density_Stations - 1; i > 0; i--) {
	//			if ((ThermoTables[i][j].Enthalpy - ThermoTables[i - 1][j].Enthalpy)
	//					> 1.1 * ThermoTables[i - 1][j].Enthalpy) {
	//				Indexes_of_two_phase[i][j] = -10;
	//			}
	//		}
	//	}
	//	for (int i =0;i<Table_Density_Stations; i++) {
	//			for (int j = 0; j < Table_Pressure_Stations; j++) {
	//				cout<<Indexes_of_two_phase[i][j]<<", ";
	//		}
	//		cout<<endl;
	//	}
	//

	delete[] SaturationTables;

	for (int i = 0; i < Table_Density_Stations; i++) {
		delete[] Indexes_of_two_phase[i];
	}
	delete[] Indexes_of_two_phase;
}

void CLookUpTable::LookUpTable_Load_CFX(string filename,
		bool read_saturation_properties) {
	//Load the table from a CFX type format file. However, the temperature
	//and the StaticEnergy have been added to the format as they are needed
	//directly in the table.
	int N_PARAM = 0;
	int set_x = 0;
	int set_y = 0;
	int var_steps = 0;
	int var_scanned = 0;

	string line;
	string value;

	ifstream table(filename.c_str());
	assert(table.is_open());
	//Go through all lines in the table file.
	while (getline(table, line)) {
		unsigned int found;
		found = line.find("$$PARAM");
		if (found < 10) {
			getline(table, line);
			istringstream in(line);
			in >> N_PARAM;
			N_PARAM++;
		}
		for (int var = var_scanned; var < N_PARAM + 1; var++) {
			string svar =
					static_cast<ostringstream*>(&(ostringstream() << var))->str();
			found = line.find("$TABLE_" + svar);
			if (found < 10) {
				var_scanned = var;
				getline(table, line);
				istringstream in(line);
				int x, y;
				in >> x >> y;
				//Create the actual LUT of CThermoLists which is used in the FluidModel
				if (var == 1) {
					ThermoTables = new CThermoList*[x];
					for (int i = 0; i < x; i++) {
						ThermoTables[i] = new CThermoList[y];
					}
					//If table is to be later split up into 2phase region and superheated vapor
					//the saturation properties should also be captured.
					if (read_saturation_properties) {
						SaturationTables = new CThermoList[y];
					}

					//Note down the dimensions of the table
					set_x = x;
					set_y = y;

					//The first axis is known to be the density values of the table
					//The table limits are stored for later checks of values becoming inconsistent
					Density_Table_Limits[0] = 10E15;
					Density_Table_Limits[1] = 0;
					var_steps = 10;
					su2double* vD = new su2double[set_x];

					//Go through all lines which are known to contain density
					//values (10 on each line)
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++) {
						getline(table, line);
						istringstream inD(line);
						if ((set_x - k * 10) < 10) {
							var_steps = (set_x - k * 10);
						}
						for (int i = 0; i < var_steps; i++) {
							inD >> vD[10 * k + i];

							vD[10 * k + i] = vD[10 * k + i] / Density_Reference_Value;

							if (vD[10 * k + i] > Density_Table_Limits[1]) {
								Density_Table_Limits[1] = vD[10 * k + i];
							}
							if (vD[10 * k + i] < Density_Table_Limits[0]) {
								Density_Table_Limits[0] = vD[10 * k + i];
							}
						}
					}
					for (int i = 0; i < set_x; i++) {
						for (int j = 0; j < set_y; j++) {
							ThermoTables[i][j].Density = vD[i];
						}
					}
					delete[] vD;

					//Fill in the pressures in the same way as the densities
					su2double* vP = new su2double[set_y];
					var_steps = 10;
					Pressure_Table_Limits[0] = 10E15;	//lower limit
					Pressure_Table_Limits[1] = 0;	//upper limit
					//Each line contains at most 10 pressure values
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++) {

						getline(table, line);
						istringstream inP(line);
						//Check if line contains less than 10 values
						if ((set_y - k * 10) < 10) {
							var_steps = (set_y - k * 10);
						}
						for (int j = 0; j < var_steps; j++) {
							inP >> vP[10 * k + j];

							vP[10 * k + j] = vP[10 * k + j] / Pressure_Reference_Value;

							if (vP[10 * k + j] > Pressure_Table_Limits[1]) {
								Pressure_Table_Limits[1] = vP[10 * k + j];
							}
							if (vP[10 * k + j] < Pressure_Table_Limits[0]) {
								Pressure_Table_Limits[0] = vP[10 * k + j];
							}
						}
					}
					for (int i = 0; i < set_x; i++) {
						for (int j = 0; j < set_y; j++) {
							ThermoTables[i][j].Pressure = vP[j];
						}
					}
					delete[] vP;
				}
				// Check that all encountered tables adhere to the same x,y dimensions
				// otherwise throw an error
				else if (x != set_x && y != set_y) {
					cerr
					<< "The encountered dimensions of the CFX table are not the same throughout. They should be; for this to work.\n";

				}
				//Go through each one of the variables of interest
				//ENTHALPY TABLE
				if (var == 1) {
					//The pressure and density lines have already been skipped for Var 1
					Enthalpy_Table_Limits[0] = 10E20;					//lower limit
					Enthalpy_Table_Limits[1] = 0;					//upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] / pow(Velocity_Reference_Value, 2);
								}
							}
							ThermoTables[i][j].Enthalpy = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Enthalpy_Table_Limits[1]) {
								Enthalpy_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Enthalpy_Table_Limits[0]) {
								Enthalpy_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}
					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//First skip the saturation temperature values
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
							}
						}
						//Now load the saturation property
						su2double inp[10];
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] / pow(Velocity_Reference_Value, 2);
								}
							}
							SaturationTables[j].Enthalpy = inp[j % 10];
						}
					}
				}

				//SOUNDS SPEED (SQUARED) TABLE
				if (var == 2) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density (already imported)
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure (already imported)

					SoundSpeed2_Table_Limits[0] = 10E20; //lower limit
					SoundSpeed2_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = pow(inp[z], 2);
									inp[z] = inp[z] / pow(Velocity_Reference_Value, 2);
								}
							}
							ThermoTables[i][j].SoundSpeed2 = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > SoundSpeed2_Table_Limits[1]) {
								SoundSpeed2_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < SoundSpeed2_Table_Limits[0]) {
								SoundSpeed2_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}
					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//First skip the saturation temperature values
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
							}
						}
						//Now load the saturation property
						su2double inp[10];
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = pow(inp[z], 2);
									inp[z] = inp[z] / pow(Velocity_Reference_Value, 2);
								}
							}
							SaturationTables[j].SoundSpeed2 = inp[j % 10];
						}
					}
				}
				//CP TABLE
				if (var == 5) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density (already imported)
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure (already imported)

					Cp_Table_Limits[0] = 10E20; //lower limit
					Cp_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * Temperature_Reference_Value
											/ pow(Velocity_Reference_Value, 2);
								}
							}
							ThermoTables[i][j].Cp = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Cp_Table_Limits[1]) {
								Cp_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Cp_Table_Limits[0]) {
								Cp_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}
					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//First skip the saturation temperature values
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
							}
						}
						//Now load the saturation property
						su2double inp[10];
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * Temperature_Reference_Value
											/ pow(Velocity_Reference_Value, 2);
								}
							}
							SaturationTables[j].Cp = inp[j % 10];
						}
					}
				}
				//ENTROPY TABLE
				if (var == 7) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density (already imported)
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure (already imported)

					Entropy_Table_Limits[0] = 10E20; //lower limit
					Entropy_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * Temperature_Reference_Value
											/ pow(Velocity_Reference_Value, 2);
								}
							}
							ThermoTables[i][j].Entropy = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Entropy_Table_Limits[1]) {
								Entropy_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Entropy_Table_Limits[0]) {
								Entropy_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}
					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//First skip the saturation temperature values
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
							}
						}
						//Now load the saturation property
						su2double inp[10];
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * Temperature_Reference_Value
											/ pow(Velocity_Reference_Value, 2);
								}
							}
							SaturationTables[j].Entropy = inp[j % 10];
						}
					}
				}
				//MU TABLE
				if (var == 8) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density (already imported)
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure (already imported)

					Mu_Table_Limits[0] = 10E20; //lower limit
					Mu_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].Mu = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Mu_Table_Limits[1]) {
								Mu_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Mu_Table_Limits[0]) {
								Mu_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}
					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//First skip the saturation temperature values
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
							}
						}
						//Now load the saturation property
						su2double inp[10];
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							SaturationTables[j].Mu = inp[j % 10];
						}
					}
				}
				//KT TABLE
				if (var == 9) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density (already imported)
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure (already imported)

					Kt_Table_Limits[0] = 10E20; //lower limit
					Kt_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].Kt = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Kt_Table_Limits[1]) {
								Kt_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Kt_Table_Limits[0]) {
								Kt_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}
					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//First skip the saturation temperature values
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
							}
						}
						//Now load the saturation property
						su2double inp[10];
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							SaturationTables[j].Kt = inp[j % 10];
						}
					}
				}
				//dPdrho_e TABLE
				if (var == 10) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density (already imported)
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure (already imported)

					dPdrho_e_Table_Limits[0] = 10E20; //lower limit
					dPdrho_e_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * Density_Reference_Value
											/ Pressure_Reference_Value;
								}
							}
							ThermoTables[i][j].dPdrho_e = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > dPdrho_e_Table_Limits[1]) {
								dPdrho_e_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < dPdrho_e_Table_Limits[0]) {
								dPdrho_e_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}
					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//First skip the saturation temperature values
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
							}
						}
						//Now load the saturation property
						su2double inp[10];
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * Density_Reference_Value
											/ Pressure_Reference_Value;
								}
							}
							SaturationTables[j].dPdrho_e = inp[j % 10];
						}
					}
				}
				//dPde_rho TABLE
				if (var == 11) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density (already imported)
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure (already imported)

					dPde_rho_Table_Limits[0] = 10E20; //lower limit
					dPde_rho_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * pow(Velocity_Reference_Value, 2)
													/ Pressure_Reference_Value;
								}
							}
							ThermoTables[i][j].dPde_rho = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > dPde_rho_Table_Limits[1]) {
								dPde_rho_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < dPde_rho_Table_Limits[0]) {
								dPde_rho_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}
					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//First skip the saturation temperature values
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
							}
						}
						//Now load the saturation property
						su2double inp[10];
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * pow(Velocity_Reference_Value, 2)
													/ Pressure_Reference_Value;
								}
							}
							SaturationTables[j].dPde_rho = inp[j % 10];
						}
					}
				}
				//dTdrho_e TABLE
				if (var == 12) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density (already imported)
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure (already imported)

					dTdrho_e_Table_Limits[0] = 10E20; //lower limit
					dTdrho_e_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * Density_Reference_Value
											/ Temperature_Reference_Value;
								}
							}
							ThermoTables[i][j].dTdrho_e = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > dTdrho_e_Table_Limits[1]) {
								dTdrho_e_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < dTdrho_e_Table_Limits[0]) {
								dTdrho_e_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}
					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//First skip the saturation temperature values
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
							}
						}
						//Now load the saturation property
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * Density_Reference_Value
											/ Temperature_Reference_Value;
								}
							}
							SaturationTables[j].dTdrho_e = inp[j % 10];
						}
					}
				}
				//dTde_rho TABLE
				if (var == 13) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density (already imported)
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure (already imported)

					dTde_rho_Table_Limits[0] = 10E20; //lower limit
					dTde_rho_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * pow(Velocity_Reference_Value, 2)
													/ Temperature_Reference_Value;
								}
							}
							ThermoTables[i][j].dTde_rho = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > dTde_rho_Table_Limits[1]) {
								dTde_rho_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < dTde_rho_Table_Limits[0]) {
								dTde_rho_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}
					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//First skip the saturation temperature values
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
							}
						}
						//Now load the saturation property
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] * pow(Velocity_Reference_Value, 2)
													/ Temperature_Reference_Value;
								}
							}
							SaturationTables[j].dTde_rho = inp[j % 10];
						}
					}
				}
				//TEMPERATURE TABLE
				if (var == 15) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density (already imported)
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure (already imported)

					Temperature_Table_Limits[0] = 10E20; //lower limit
					Temperature_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] / Temperature_Reference_Value;
								}
							}
							ThermoTables[i][j].Temperature = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Temperature_Table_Limits[1]) {
								Temperature_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Temperature_Table_Limits[0]) {
								Temperature_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}
					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//Load saturation temperature
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] / Temperature_Reference_Value;
								}
							}
							SaturationTables[j].Temperature = inp[j % 10];
						}
					}
				}
				//STATIC ENERGY TABLE
				if (var == 16) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density (already imported)
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure (already imported)

					StaticEnergy_Table_Limits[0] = 10E20; //lower limit
					StaticEnergy_Table_Limits[1] = 0; //upper limit

					su2double inp[10];
					//Gp through all lines of data
					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i));
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] / pow(Velocity_Reference_Value, 2);
								}
							}
							ThermoTables[i][j].StaticEnergy = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > StaticEnergy_Table_Limits[1]) {
								StaticEnergy_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < StaticEnergy_Table_Limits[0]) {
								StaticEnergy_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					//Also load the saturation properties if desired.
					if (read_saturation_properties) {
						//First skip the saturation temperature values
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
							}
						}
						//Now load the saturation property
						su2double inp[10];
						for (int j = 0; j < set_y; j++) {
							if ((j) % 10 == 0) {
								getline(table, line);
								istringstream in(line);
								var_steps = 10;
								if ((set_y - j) < 10)
									var_steps = (set_y - j);
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = inp[z] / pow(Velocity_Reference_Value, 2);
								}
							}
							SaturationTables[j].StaticEnergy = inp[j % 10];
						}
					}
				}
			}
		}

	}
	Table_Density_Stations = set_x;
	Table_Pressure_Stations = set_y;
	table.close();
}

