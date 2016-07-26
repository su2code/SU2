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

CLookUpTable::CLookUpTable(CConfig *config, bool dimensional) :
		CFluidModel() {
	LUT_Debug_Mode = config->GetLUT_Debug_Mode();
	rank = MASTER_NODE;
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
	skewed_linear_table = false;
	//Detect cfx filetype
	if ((config->GetLUTFileName()).find(".rgp") != string::npos) {
		if (rank == MASTER_NODE) {
			cout << "CFX type LUT found" << endl;
		}
		LookUpTable_Load_CFX(config->GetLUTFileName());
	}
	//Detect dat file type
	else if ((config->GetLUTFileName()).find(".dat") != string::npos) {
		if (rank == MASTER_NODE) {
			cout << "DAT type LUT found" << endl;
		}
		LookUpTable_Load_DAT(config->GetLUTFileName());
	} else {
		if (rank == MASTER_NODE) {
			cout << "No recognized LUT format found, exiting!" << endl;
		}
		exit(EXIT_FAILURE);
	}
	if (ThermoTables_Pressure[0][0]
			!= ThermoTables_Pressure[Table_Density_Stations - 1][0]) {
		skewed_linear_table = true;
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Interpolation_Coeff[i][j] = -1.0;
		}
	}
	if (rank == MASTER_NODE) {
		// Give the user some information on the size of the table
		cout << "Table_Pressure_Stations  : " << Table_Pressure_Stations << endl;
		cout << "Table_Density_Stations: " << Table_Density_Stations << endl;
		cout<<"Print LUT errors? (LUT_Debug_Mode):  "<<LUT_Debug_Mode<<endl;
		// Building an KD_tree for the HS thermopair
		cout << "Building HS_tree" << endl;
	}

	su2double* xtemp = new su2double[Table_Density_Stations
			* Table_Pressure_Stations];
	su2double* ytemp = new su2double[Table_Density_Stations
			* Table_Pressure_Stations];
	int* itemp = new int[Table_Density_Stations * Table_Pressure_Stations];
	// Deep copy the x,y, and index values for the KD_tree
	for (int i = 0; i < Table_Density_Stations; i++) {
		for (int j = 0; j < Table_Pressure_Stations; j++) {
			xtemp[Table_Pressure_Stations * i + j] = ThermoTables_Enthalpy[i][j];
			ytemp[Table_Pressure_Stations * i + j] = ThermoTables_Entropy[i][j];
			itemp[Table_Pressure_Stations * i + j] = Table_Pressure_Stations * i + j;
		}
	}
	// Recursive KD_tree build starts here
	HS_tree = KD_Tree(xtemp, ytemp, itemp,
			Table_Pressure_Stations * Table_Density_Stations, 0);
	if (rank == MASTER_NODE)
		cout << "HS_tree built" << endl;
}

CLookUpTable::~CLookUpTable(void) {
	// Delete the lut
	for (int i = 0; i < Table_Density_Stations; i++) {
		delete[] ThermoTables_StaticEnergy[i];
		delete[] ThermoTables_Entropy[i];
		delete[] ThermoTables_Enthalpy[i];
		delete[] ThermoTables_Density[i];
		delete[] ThermoTables_Pressure[i];
		delete[] ThermoTables_SoundSpeed2[i];
		delete[] ThermoTables_Temperature[i];
		delete[] ThermoTables_dPdrho_e[i];
		delete[] ThermoTables_dPde_rho[i];
		delete[] ThermoTables_dTdrho_e[i];
		delete[] ThermoTables_dTde_rho[i];
		delete[] ThermoTables_Cp[i];
		delete[] ThermoTables_Mu[i];
		delete[] ThermoTables_dmudrho_T[i];
		delete[] ThermoTables_dmudT_rho[i];
		delete[] ThermoTables_Kt[i];
		delete[] ThermoTables_dktdrho_T[i];
		delete[] ThermoTables_dktdT_rho[i];
	}
	delete[] ThermoTables_StaticEnergy;
	delete[] ThermoTables_Entropy;
	delete[] ThermoTables_Enthalpy;
	delete[] ThermoTables_Density;
	delete[] ThermoTables_Pressure;
	delete[] ThermoTables_SoundSpeed2;
	delete[] ThermoTables_Temperature;
	delete[] ThermoTables_dPdrho_e;
	delete[] ThermoTables_dPde_rho;
	delete[] ThermoTables_dTdrho_e;
	delete[] ThermoTables_dTde_rho;
	delete[] ThermoTables_Cp;
	delete[] ThermoTables_Mu;
	delete[] ThermoTables_dmudrho_T;
	delete[] ThermoTables_dmudT_rho;
	delete[] ThermoTables_Kt;
	delete[] ThermoTables_dktdrho_T;
	delete[] ThermoTables_dktdT_rho;

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
			i = i + 1;
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

void CLookUpTable::Search_NonEquispaced_Rho_Index(su2double rho) {
	{
		su2double grad, x00, y00;
		//  Determine the I index with binary search (rho is not assumed equispaced)
		while (UpperI - LowerI > 1) {
			middleI = (UpperI + LowerI) / 2;
			x00 = ThermoTables_Density[middleI][LowerJ];
			grad = ThermoTables_Density[middleI + 1][LowerJ] - x00;
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
	}
}
void CLookUpTable::Search_NonEquispaced_P_Index(su2double P) {
	su2double grad, x00, y00, y01, x01, RunVal;
	//Determine the J index using a binary search, and not assuming P is equispaced
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		x00 = ThermoTables_Pressure[LowerI][middleJ];
		grad = ThermoTables_Pressure[LowerI][middleJ + 1] - x00;
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
}

void CLookUpTable::Zig_Zag_Search(su2double x, su2double y,
		su2double **ThermoTables_X, su2double **ThermoTables_Y) {
	su2double dx, dy, x00, y00, dx10, dx01, dx11, dy10, dy01, dy11;
	bool BOTTOM, TOP, LEFT, RIGHT, found = false; //check if bellow BOTTOM, above, TOP, left of LEFT, and right of RIGHT
	for (int k = 0; k < 20 and not found; k++) //20 is arbitrary and used primarily to avoid a while loop which could get stuck
			{
		UpperI = LowerI + 1;
		UpperJ = LowerJ + 1;
		x00 = ThermoTables_X[LowerI][LowerJ];
		y00 = ThermoTables_Y[LowerI][LowerJ];
		dx = x - x00;
		dy = y - y00;
		dx01 = ThermoTables_X[LowerI][UpperJ] - x00;
		dy01 = ThermoTables_Y[LowerI][UpperJ] - y00;
		dx10 = ThermoTables_X[UpperI][LowerJ] - x00;
		dy10 = ThermoTables_Y[UpperI][LowerJ] - y00;
		dx11 = ThermoTables_X[UpperI][UpperJ] - x00;
		dy11 = ThermoTables_Y[UpperI][UpperJ] - y00;

		BOTTOM = (dy * dx10) < (dx * dy10);
		TOP = ((dy - dy01) * (dx11 - dx01)) > ((dy11 - dy01) * (dx - dx01));
		RIGHT = ((dx - dx10) * (dy11 - dy10)) > ((dx11 - dx10) * (dy - dy10));
		LEFT = (dx * dy01) < (dx01 * dy);
		//Check BOTTOM quad boundary
		if (BOTTOM and !TOP) {
			if (LowerJ != 0)
				LowerJ--;
		}
		//Check RIGHT quad boundary
		else if (RIGHT and !LEFT) {
			if (LowerI != (Table_Density_Stations - 2))
				LowerI++;
		}
		//Check TOP quad boundary
		else if (TOP and !BOTTOM) {
			if (LowerJ != (Table_Pressure_Stations - 2))
				LowerJ++;
		}
		//Check LEFT quad boundary
		else if (LEFT and !RIGHT) {
			if (LowerI != 0)
				LowerI--;
		} else {
			found = true;
		}
	}
}

void CLookUpTable::Search_i_for_X_given_j(su2double x, su2double y,
		su2double **ThermoTables_X, su2double **ThermoTables_Y) {
	su2double grad, x00, y00, y01, x01, RunVal;
	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		//Use interpolated T as the running variable for the search (RunVal)
		y00 = ThermoTables_Y[middleI][LowerJ];
		y01 = ThermoTables_Y[middleI][UpperJ];
		x00 = ThermoTables_X[middleI][LowerJ];
		x01 = ThermoTables_X[middleI][UpperJ];
		grad = ThermoTables_X[UpperI][LowerJ] - x00;
		RunVal = x00 + (x01 - x00) / (y01 - y00) * (y - y00);
		if (RunVal * grad > x * grad) {
			UpperI = middleI;
		} else if (RunVal * grad < x * grad) {
			LowerI = middleI;
		} else if (RunVal == x) {
			LowerI = middleI;
			UpperI = LowerI + 1;
			break;
		}
	}
}

void CLookUpTable::Search_j_for_Y_given_i(su2double x, su2double y,
		su2double **ThermoTables_X, su2double **ThermoTables_Y) {
	su2double RunVal;
	su2double grad, x00, y00, y10, x10;
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2; /*!< \brief Splitting index for the search */

		// The variable names is composed of a (i,j) pair
		y00 = ThermoTables_Y[LowerI][middleJ];
		y10 = ThermoTables_Y[UpperI][middleJ];
		x00 = ThermoTables_X[LowerI][middleJ];
		x10 = ThermoTables_X[UpperI][middleJ];
		//The search variable in j should be interpolated in i as well
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (x - x00);
		grad = ThermoTables_Y[LowerI][middleJ + 1] - y00;
		if (RunVal * grad > y * grad) {
			UpperJ = middleJ;
		} else if (RunVal * grad < y * grad) {
			LowerJ = middleJ;
		} else if (RunVal == y) {
			LowerJ = middleJ;
			UpperJ = LowerJ + 1;
			break;
		}
	}
}

void CLookUpTable::Search_Linear_Skewed_Table(su2double x, su2double P,
		su2double **ThermoTables_X) {
	su2double RunVal, rho;
	su2double grad, x00, y00, y10, x10;
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2; /*!< \brief Splitting index for the search */

		// The variable names is composed of a (i,j) pair
		y00 = ThermoTables_Pressure[0][middleJ];
		y10 = ThermoTables_Pressure[Table_Density_Stations - 1][middleJ];
		x00 = ThermoTables_Density[0][middleJ];
		x10 = ThermoTables_Density[Table_Density_Stations - 1][middleJ];
		//Using the input pressure and a given middleJ value,
		//the corresponding density and iIndex may be determined
		rho = (P - y00) * (x10 - x00) / (y10 - y00) + x00;
		//If density is greater than the limits of the table,
		//then search the upper part, if it is lower than the limits,
		//search the lower part of the table
		if (rho > Density_Table_Limits[1]) {
			if (y10 > y00) {
				LowerJ = middleJ;
			} else if (y10 < y00) {
				UpperJ = middleJ;
			}
		} else if (rho < Density_Table_Limits[0]) {
			if (y10 > y00) {
				UpperJ = middleJ;
			} else if (y10 < y00) {
				LowerJ = middleJ;
			}
		} else {
			//If the density is within limits, calculate the iIndex that
			//it corresponds to
			UpperI = Table_Density_Stations - 1;
			LowerI = 0;
			Search_NonEquispaced_Rho_Index(rho);
			//Now an i,j index pair has been found, and can be used to
			//guide the recursion further
			y00 = ThermoTables_X[LowerI][middleJ];
			y10 = ThermoTables_X[UpperI][middleJ];
			x00 = ThermoTables_Density[LowerI][middleJ];
			x10 = ThermoTables_Density[UpperI][middleJ];

			RunVal = y00 + (y10 - y00) / (x10 - x00) * (rho - x00);
			grad = ThermoTables_X[LowerI][middleJ + 1] - y00;
			if (RunVal * grad > x * grad) {
				UpperJ = middleJ;
			} else if (RunVal * grad < x * grad) {
				LowerJ = middleJ;
			} else if (RunVal == x) {
				LowerJ = middleJ;
				UpperJ = LowerJ + 1;
				break;
			}
		}
	}
}

void CLookUpTable::SetTDState_rhoe(su2double rho, su2double e) {
	// Check if inputs are in total range (necessary but not sufficient condition)
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "RHOE Input Density out of bounds\n";
		}
		if ((e > StaticEnergy_Table_Limits[1])
				or (e < StaticEnergy_Table_Limits[0])) {
			cerr << "RHOE Input StaticEnergy out of bounds\n";
		}
	}

	// Starting values for the search
	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	//Searches involving density are independent of the skewness of the table
	//Fix the i index first
	Search_NonEquispaced_Rho_Index(rho);
	//Having found the i index, search in j
	Search_j_for_Y_given_i(rho, e, ThermoTables_Density,
			ThermoTables_StaticEnergy);

	//Now use the quadrilateral which contains the point to interpolate
	//Determine the interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(rho, e, ThermoTables_Density,
			ThermoTables_StaticEnergy, "RHOE");

	//Interpolate the fluid properties
	StaticEnergy = e;
	Density = rho;
	Entropy = Interpolate_2D_Bilinear(ThermoTables_Entropy);
	Pressure = Interpolate_2D_Bilinear(ThermoTables_Pressure);
	//Enthalpy = Interpolate_2D_Bilinear(ThermoTables_Enthalpy);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	Temperature = Interpolate_2D_Bilinear(ThermoTables_Temperature);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);

	//Check that the interpolated density and pressure are within LUT limits
	Check_Interpolated_PRHO_Limits("RHOE");
}

void CLookUpTable::SetTDState_PT(su2double P, su2double T) {
	// Check if inputs are in total range (necessary but not sufficient condition)
	if (rank == MASTER_NODE) {
		if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
			cerr << "PT Input Pressure out of bounds\n";
		}
		if ((T > Temperature_Table_Limits[1])
				or (T < Temperature_Table_Limits[0])) {
			cerr << "PT Input Temperature out of bounds\n";
		}
	}

	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;
	if (not skewed_linear_table) {
		//Determine the j index first
		Search_NonEquispaced_P_Index(P);
		//Having fixes the j index, find the i index
		Search_i_for_X_given_j(T, P, ThermoTables_Temperature,
				ThermoTables_Pressure);
	} else if (skewed_linear_table) {
		Search_Linear_Skewed_Table(T, P, ThermoTables_Temperature);
		//Finish of the search to check that the correct quad has indeed been selected
		//Although brief inspection shows correct element to be always selected before zig zag.
		Zig_Zag_Search(T, P, ThermoTables_Temperature, ThermoTables_Pressure);
	}

	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(T, P, ThermoTables_Temperature,
			ThermoTables_Pressure, "PT");

	//Interpolate the fluid properties
	Pressure = P;
	Temperature = T;
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	//Enthalpy = Interpolate_2D_Bilinear(ThermoTables_Enthalpy);
	Entropy = Interpolate_2D_Bilinear(ThermoTables_Entropy);
	Density = Interpolate_2D_Bilinear(ThermoTables_Density);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);

	//Check that the interpolated density and pressure are within LUT limits
	Check_Interpolated_PRHO_Limits("PT");
}

void CLookUpTable::SetTDState_Prho(su2double P, su2double rho) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
			cerr << "PRHO Input Pressure out of bounds\n";
		}
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "PRHO Input Density out of bounds\n";
		}
	}

	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	if (not skewed_linear_table) {
		Search_NonEquispaced_Rho_Index(rho);
		Search_NonEquispaced_P_Index(P);
	} else if (skewed_linear_table) {
		Search_NonEquispaced_Rho_Index(rho);
		Search_j_for_Y_given_i(rho, P, ThermoTables_Density, ThermoTables_Pressure);
	}

	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(rho, P, ThermoTables_Density,
			ThermoTables_Pressure, "PRHO");
	//Interpolate the fluid properties
	Pressure = P;
	Density = rho;
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	//Enthalpy = Interpolate_2D_Bilinear(ThermoTables_Enthalpy);
	Entropy = Interpolate_2D_Bilinear(ThermoTables_Entropy);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	Temperature = Interpolate_2D_Bilinear(ThermoTables_Temperature);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);

	//Check that the interpolated density and pressure are within LUT limits
	Check_Interpolated_PRHO_Limits("PRHO");
}

void CLookUpTable::SetEnergy_Prho(su2double P, su2double rho) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
			cerr << "PRHO Input Pressure out of bounds\n";
		}
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "PRHO Input Density out of bounds\n";
		}
	}

	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	if (not skewed_linear_table) {
		Search_NonEquispaced_Rho_Index(rho);
		Search_NonEquispaced_P_Index(P);
	} else if (skewed_linear_table) {
		Search_NonEquispaced_Rho_Index(rho);
		Search_j_for_Y_given_i(rho, P, ThermoTables_Density, ThermoTables_Pressure);
	}

	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(rho, P, ThermoTables_Density,
			ThermoTables_Pressure, "PRHO");
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	Pressure = P;
	Density = rho;

	Check_Interpolated_PRHO_Limits("PRHO (energy)");
}

void CLookUpTable::SetTDState_hs(su2double h, su2double s) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		if ((h > Enthalpy_Table_Limits[1]) or (h < Enthalpy_Table_Limits[0])) {
			cerr << "HS Input Enthalpy out of bounds\n";
		}
		if ((s > Entropy_Table_Limits[1]) or (s < Entropy_Table_Limits[0])) {
			cerr << "HS Input Entropy out of bounds\n";
		}
	}

	LowerI = HS_tree->Flattened_Point_Index[HS_tree->Branch_Dimension / 2]
			/ Table_Pressure_Stations;
	LowerJ = HS_tree->Flattened_Point_Index[HS_tree->Branch_Dimension / 2]
			% Table_Pressure_Stations;

	//Simply find the closest point with the KD_tree and use it for the intrpolaiton
	int N = 1;
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
		best_dist[i] = HUGE_VAL;
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
	LowerI = Nearest_Neighbour_iIndex[0];
	LowerJ = Nearest_Neighbour_jIndex[0];

	//Using the closest element found in the KD_tree as a starting point, now find the closest
	//quadrilateral containing the point using a simple zigzag search method
	Zig_Zag_Search(h, s, ThermoTables_Enthalpy, ThermoTables_Entropy);

	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(h, s, ThermoTables_Enthalpy,
			ThermoTables_Entropy, "HS");

	//Interpolate the fluid properties
	//Enthalpy = h;
	Entropy = s;
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	Pressure = Interpolate_2D_Bilinear(ThermoTables_Pressure);
	Density = Interpolate_2D_Bilinear(ThermoTables_Density);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	Temperature = Interpolate_2D_Bilinear(ThermoTables_Temperature);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);

	Check_Interpolated_PRHO_Limits("HS");
	delete[] best_dist;
	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;

}

void CLookUpTable::SetTDState_Ps(su2double P, su2double s) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
			cerr << "PS Input Pressure out of bounds\n";
		}
		if ((s > Entropy_Table_Limits[1]) or (s < Entropy_Table_Limits[0])) {
			cerr << "PS Input Entropy  out of bounds\n";
		}
	}

	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	if (not skewed_linear_table) {
		//First fix the j index
		Search_NonEquispaced_P_Index(P);
		//Then find the i index, given a fixed j
		Search_i_for_X_given_j(s, P, ThermoTables_Entropy, ThermoTables_Pressure);
	} else if (skewed_linear_table) {
		Search_Linear_Skewed_Table(s, P, ThermoTables_Entropy);
		//Finish of the search to check that the correct quad has indeed been selected
		Zig_Zag_Search(s, P, ThermoTables_Entropy, ThermoTables_Pressure);
	}

	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(s, P, ThermoTables_Entropy,
			ThermoTables_Pressure, "PS");

	//Interpolate the fluid properties
	Entropy = s;
	Pressure = P;
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	//Enthalpy = Interpolate_2D_Bilinear(ThermoTables_Enthalpy);
	Density = Interpolate_2D_Bilinear(ThermoTables_Density);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	Temperature = Interpolate_2D_Bilinear(ThermoTables_Temperature);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);

	//Check that the interpolated density and pressure are within LUT limits
	Check_Interpolated_PRHO_Limits("PS");
}

void CLookUpTable::SetTDState_rhoT(su2double rho, su2double T) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "RHOT Input Density out of bounds\n";
		}
		if ((T > Temperature_Table_Limits[1])
				or (T < Temperature_Table_Limits[0])) {
			cerr << "RHOT Input Temperature out of bounds\n";
		}
	}
	// Linear interpolation requires 4 neighbors to be selected from the LUT

	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	if (not skewed_linear_table) {
		Search_NonEquispaced_Rho_Index(rho);
		Search_j_for_Y_given_i(rho, T, ThermoTables_Density,
				ThermoTables_Temperature);
	}

	//Determine the interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(rho, T, ThermoTables_Density,
			ThermoTables_Temperature, "RHOT");

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	//Enthalpy = Interpolate_2D_Bilinear(ThermoTables_Enthalpy");
	Entropy = Interpolate_2D_Bilinear(ThermoTables_Entropy);
	Pressure = Interpolate_2D_Bilinear(ThermoTables_Pressure);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);

	//Check that the interpolated density and pressure are within LUT limits
	Check_Interpolated_PRHO_Limits("RHOT");
}

void CLookUpTable::Check_Interpolated_PRHO_Limits(string interpolation_case) {
	//Check that the interpolated density and pressure are within LUT limits
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		if ((Density > Density_Table_Limits[1])
				or (Density < Density_Table_Limits[0])) {
			cerr << interpolation_case << " Interpolated Density out of bounds\n";
		}
		if ((Pressure > Pressure_Table_Limits[1])
				or (Pressure < Pressure_Table_Limits[0])) {
			cerr << interpolation_case << " Interpolated Pressure out of bounds\n";
		}
	}
}

inline void CLookUpTable::Gaussian_Inverse(int nDim) {
	//A temporary matrix to hold the inverse, dynamically allocated
	su2double **temp = new su2double*[nDim];
	for (int i = 0; i < nDim; i++) {
		temp[i] = new su2double[2 * nDim];
	}

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
			Interpolation_Coeff[i][j] = temp[i][j + nDim];
		}
	}
	//Delete dynamic template
	for (int i = 0; i < nDim; i++) {
		delete[] temp[i];
	}
	delete[] temp;
	return;
}

void CLookUpTable::Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(su2double x,
		su2double y, su2double **ThermoTables_X, su2double **ThermoTables_Y,
		std::string grid_var) {
	//The x,y coordinates of the quadrilateral
	su2double x00, y00, x10, x01, x11, y10, y01, y11;

	x00 = ThermoTables_X[LowerI][LowerJ];
	y00 = ThermoTables_Y[LowerI][LowerJ];
	x01 = ThermoTables_X[LowerI][UpperJ];
	y01 = ThermoTables_Y[LowerI][UpperJ];
	x10 = ThermoTables_X[UpperI][LowerJ];
	y10 = ThermoTables_Y[UpperI][LowerJ];
	x11 = ThermoTables_X[UpperI][UpperJ];
	y11 = ThermoTables_Y[UpperI][UpperJ];

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
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		BOTTOM = (dy * dx10) < (dx * dy10);
		TOP = ((dy - dy01) * (dx11 - dx01)) > ((dy11 - dy01) * (dx - dx01));
		RIGHT = ((dx - dx10) * (dy11 - dy10)) > ((dx11 - dx10) * (dy - dy10));
		LEFT = (dx * dy01) < (dx01 * dy);
		OUT_OF_BOUNDS = false;
		//Check BOTTOM quad boundary
		if (BOTTOM and !TOP) {
			OUT_OF_BOUNDS = true;
			//Check if the point is also outside the LUT
			if (LowerJ == 0) {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies below the LUT\n";
			} else {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies below bottom boundary of selected quad\n";
			}
		}
		//Check RIGHT quad boundary
		if (RIGHT and !LEFT) {
			OUT_OF_BOUNDS = true;
			//Check if the point is also outside the LUT
			if (LowerI == (Table_Density_Stations - 2)) {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies right of the LUT\n";
			} else {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies to the right of the boundary of selected quad\n";
			}
		}
		//Check TOP quad boundary
		if (TOP and !BOTTOM) {
			OUT_OF_BOUNDS = true;
			//Check if the point is also outside the LUT
			if (LowerJ == (Table_Pressure_Stations - 2)) {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies above the LUT\n";
			} else {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< +" interpolation point lies above the boundary of selected quad\n";
			}
		}
		//Check LEFT quad boundary
		if (LEFT and !RIGHT) {
			OUT_OF_BOUNDS = true;
			//Check if the point is also outside the LUT
			if (LowerI == 0) {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies left of the LUT\n";
			} else {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< +" interpolation point lies to the left of the boundary of selected quad\n";
			}
		}
	}

	//Setup the LHM matrix for the interpolation (Vandermonde)
//	Interpolation_Matrix[0][0] = 1;
//	Interpolation_Matrix[0][1] = 0;
//	Interpolation_Matrix[0][2] = 0;
//	Interpolation_Matrix[0][3] = 0;
//
//	Interpolation_Matrix[1][0] = 1;
//	Interpolation_Matrix[1][1] = (x10 - x00);
//	Interpolation_Matrix[1][2] = (y10 - y00);
//	Interpolation_Matrix[1][3] = (x10 - x00) * (y10 - y00);
//
//	Interpolation_Matrix[2][0] = 1;
//	Interpolation_Matrix[2][1] = x01 - x00;
//	Interpolation_Matrix[2][2] = y01 - y00;
//	Interpolation_Matrix[2][3] = (x01 - x00) * (y01 - y00);
//
//	Interpolation_Matrix[3][0] = 1;
//	Interpolation_Matrix[3][1] = x11 - x00;
//	Interpolation_Matrix[3][2] = y11 - y00;
//	Interpolation_Matrix[3][3] = (x11 - x00) * (y11 - y00);
//
//	Interpolation_Matrix[0][0] = 1;
//	Interpolation_Matrix[0][1] = 0;
//	Interpolation_Matrix[0][2] = 0;
//	Interpolation_Matrix[0][3] = 0;

//Normalization
	Interpolation_Matrix[0][0] = 1;
	Interpolation_Matrix[0][1] = 1;
	Interpolation_Matrix[0][2] = 1;
	Interpolation_Matrix[0][3] = 1;

	Interpolation_Matrix[1][0] = 1;
	Interpolation_Matrix[1][1] = (x10/x00);
	Interpolation_Matrix[1][2] = (y10/y00);
	Interpolation_Matrix[1][3] = (y10/y00)*(x10/x00);

	Interpolation_Matrix[2][0] = 1;
	Interpolation_Matrix[2][1] = (x01/x00);
	Interpolation_Matrix[2][2] = (y01/y00);
	Interpolation_Matrix[2][3] = (y01/y00)/(x01/x00);

	Interpolation_Matrix[3][0] = 1;
	Interpolation_Matrix[3][1] = (x11/x00);
	Interpolation_Matrix[3][2] = (y11/y00);
	Interpolation_Matrix[3][3] = (x11/x00)*(y11/y00);


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
		d = d + Interpolation_Coeff[i][1] * (x/x00);
		d = d + Interpolation_Coeff[i][2] * (y/y00);
		d = d + Interpolation_Coeff[i][3] * (x/x00)*(y/y00);
		Interpolation_Coeff[i][0] = d;
	}

//	//Fast inverse (Sympy CAS generated)
//	su2double a,b,c,d;
//	a  = (y10-y00)/(y01-y00);
//	b  = (x01-x00)/(x10-x00);
//	c  = (x11-x00)/(x10-x00);
//	d  = (y11-y00)/(y01-y00);
//
//	Interpolation_Coeff[0][0]= 1;
//	Interpolation_Coeff[0][1]= (a*b - a*c*d + a*d - a - b*d + c*d)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//	Interpolation_Coeff[0][2]= (a*b - a*c - b*c*d + b*c - b + c*d)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//	Interpolation_Coeff[0][3]= (-a*b + a*c + b*d - c - d + 1)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//	Interpolation_Coeff[1][0]= 0;
//	Interpolation_Coeff[1][1]= d*(b - c)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//	Interpolation_Coeff[1][2]= b*c*(d - 1)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//	Interpolation_Coeff[1][3]= (-b*d + c)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//	Interpolation_Coeff[2][0]= 0;
//	Interpolation_Coeff[2][1]= a*d*(c - 1)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//	Interpolation_Coeff[2][2]= c*(a - d)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//	Interpolation_Coeff[2][3]= (-a*c + d)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//	Interpolation_Coeff[3][0]= 0;
//	Interpolation_Coeff[3][1]= -a*(b - 1)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//	Interpolation_Coeff[3][2]= -b*(a - 1)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//	Interpolation_Coeff[3][3]= (a*b - 1)/(a*b*c*d - a*b*c - a*b*d + a*c + b*d - c*d);
//
//	//The transpose allows the same coefficients to be used
//	// for all Thermo variables (need only 4 coefficients)
//	su2double dt;
//	for (int i = 0; i < 4; i++) {
//		dt = 0;
//		dt = dt + Interpolation_Coeff[i][0] * 1;
//		dt = dt + Interpolation_Coeff[i][1] * (x - x00)/(x10-x00);
//		dt = dt + Interpolation_Coeff[i][2] * (y - y00)/(y01-y00);
//		dt = dt + Interpolation_Coeff[i][3] * (x - x00)/(x10-x00) * (y - y00)/(y01-y00);
//		Interpolation_Coeff[i][0] = dt;
//	}

	return;
}

su2double CLookUpTable::Interpolate_2D_Bilinear(su2double * *ThermoTables_Z) {
	//The function values at the 4 corners of the quad
	su2double func_value_at_i0j0, func_value_at_i1j0, func_value_at_i0j1,
			func_value_at_i1j1;

	func_value_at_i0j0 = ThermoTables_Z[LowerI][LowerJ];
	func_value_at_i1j0 = ThermoTables_Z[UpperI][LowerJ];
	func_value_at_i0j1 = ThermoTables_Z[LowerI][UpperJ];
	func_value_at_i1j1 = ThermoTables_Z[UpperI][UpperJ];
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
	//
	//
	//
	//			RecordState(filename);
	//		}
	//	}
	//
}

void CLookUpTable::LookUpTable_Load_CFX(string filename) {
	//Load the table from a CFX type format file. However, the temperature
	//and the StaticEnergy have been added to the format as they are needed
	//directly in the table.
	int N_PARAM = 0;
	int var_steps = 0;
	int var_scanned = 0;

	string line;
	string value;

	ifstream table(filename.c_str());

	if (!table.is_open()) {
		if (rank == MASTER_NODE) {
			cout << "The LUT file appears to be missing!! " << filename << endl;
		}
		exit(EXIT_FAILURE);
	}
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

				//Note down the dimensions of the table
				int x, y;
				in >> x >> y;
				//Create the actual LUT of CThermoLists which is used in the FluidModel
				if (var == 1) {
					Table_Density_Stations = x;
					Table_Pressure_Stations = y;
					//Allocate the memory for the table
					LookUpTable_Malloc();
					//The first axis is known to be the density values of the table
					var_steps = 10;
					su2double* vD = new su2double[Table_Density_Stations];

					//Go through all lines which are known to contain density
					//values (10 on each line)
					for (int k = 0; k < ceil(float(Table_Density_Stations) / 10.0); k++) {
						getline(table, line);
						istringstream inD(line);
						if ((Table_Density_Stations - k * 10) < 10) {
							var_steps = (Table_Density_Stations - k * 10);
						}
						for (int i = 0; i < var_steps; i++) {
							inD >> vD[10 * k + i];
						}
					}
					//Fill the loaded values into the LUT
					for (int i = 0; i < Table_Density_Stations; i++) {
						for (int j = 0; j < Table_Pressure_Stations; j++) {
							ThermoTables_Density[i][j] = vD[i];
						}
					}
					delete[] vD;
					//Fill in the pressures in the same way as the densities
					su2double* vP = new su2double[Table_Pressure_Stations];
					var_steps = 10;

					//Each line contains at most 10 pressure values
					for (int k = 0; k < ceil(float(Table_Pressure_Stations) / 10.0);
							k++) {

						getline(table, line);
						istringstream inP(line);
						//Check if line contains less than 10 values
						if ((Table_Pressure_Stations - k * 10) < 10) {
							var_steps = (Table_Pressure_Stations - k * 10);
						}
						for (int j = 0; j < var_steps; j++) {
							inP >> vP[10 * k + j];
						}
					}
					//Save the pressure values into the LUT
					for (int i = 0; i < Table_Density_Stations; i++) {
						for (int j = 0; j < Table_Pressure_Stations; j++) {
							ThermoTables_Pressure[i][j] = vP[j];
						}
					}
					delete[] vP;
					//Finally load the Enthalpy values that var 1 is associated with
					//ENTHALPY TABLE
					CFX_Import_Table_By_Number(&table, ThermoTables_Enthalpy, false);
				}
				// Check that all encountered tables adhere to the same x,y dimensions
				// otherwise throw an error
				else if (x != Table_Density_Stations && y != Table_Pressure_Stations) {
					if (rank == MASTER_NODE) {
						cerr
								<< "The encountered dimensions of the CFX table are not the same throughout. "
										"They should be; for this to work.\n";
					}

				}
				//Go through each one of the variables of interest
				//SOUNDS SPEED (SQUARED) TABLE
				if (var == 2) {
					CFX_Import_Table_By_Number(&table, ThermoTables_SoundSpeed2, true);
				}
				//CP TABLE
				if (var == 5) {
					CFX_Import_Table_By_Number(&table, ThermoTables_Cp, true);
				}
				//ENTROPY TABLE
				if (var == 7) {
					CFX_Import_Table_By_Number(&table, ThermoTables_Entropy, true);
				}
				//dPdrho_e TABLE
				if (var == 10) {
					CFX_Import_Table_By_Number(&table, ThermoTables_dPdrho_e, true);
				}
				//dPde_rho TABLE
				if (var == 11) {
					CFX_Import_Table_By_Number(&table, ThermoTables_dPde_rho, true);
				}
				//dTdrho_e TABLE
				if (var == 12) {
					CFX_Import_Table_By_Number(&table, ThermoTables_dTdrho_e, true);
				}
				//dTde_rho TABLE
				if (var == 13) {
					CFX_Import_Table_By_Number(&table, ThermoTables_dTde_rho, true);
				}
				//TEMPERATURE TABLE
				if (var == 15) {
					CFX_Import_Table_By_Number(&table, ThermoTables_Temperature, true);
				}
				//STATIC ENERGY TABLE
				if (var == 16) {
					CFX_Import_Table_By_Number(&table, ThermoTables_StaticEnergy, true);
				}
			}
		}
	}
	table.close();
	//Non dimensionalise the table, and find the value limits
	NonDimensionalise_Table_Values();
	Find_Table_Limits();
}

void CLookUpTable::CFX_Import_Table_By_Number(ifstream *tab,
		su2double **ThermoTables_X, bool skip_prho) {
	su2double inp[10];
	int var_steps = 10;
	string line;
	if (skip_prho) {
		for (int k = 0; k < ceil(float(Table_Density_Stations) / 10.0); k++)
			getline(*tab, line); //skip density (already imported)
		for (int k = 0; k < ceil(float(Table_Pressure_Stations) / 10.0); k++)
			getline(*tab, line); //skip pressure (already imported)
	}
	for (int j = 0; j < Table_Pressure_Stations; j++) {
		for (int i = 0; i < Table_Density_Stations; i++) {
			if ((j * Table_Density_Stations + i) % 10 == 0) {
				getline(*tab, line);
				istringstream in(line);
				var_steps = 10;
				if (((Table_Density_Stations * Table_Pressure_Stations)
						- (j * Table_Density_Stations + i)) < 10)
					var_steps = ((Table_Density_Stations * Table_Pressure_Stations)
							- (j * Table_Density_Stations + i));
				for (int z = 0; z < var_steps; z++) {
					in >> inp[z];
				}
			}
			ThermoTables_X[i][j] = inp[(j * Table_Density_Stations + i) % 10];
		}
	}
}

void CLookUpTable::LookUpTable_Load_DAT(std::string filename) {
	string line;
	string value;

	ifstream table(filename.c_str());
	if (!table.is_open()) {
		if (rank == MASTER_NODE) {
			cout << "The LUT file appears to be missing!! " << filename << endl;
		}
		exit(EXIT_FAILURE);
	}
	//Go through all lines in the table file.
	getline(table, line);	//Skip the name header
	getline(table, line);	//Line with the table dimensions
	istringstream in(line);
	in >> Table_Density_Stations >> Table_Pressure_Stations;
	//Allocate the memory for the table
	LookUpTable_Malloc();

	for (int i = 0; i < Table_Density_Stations; i++) {
		for (int j = 0; j < Table_Pressure_Stations; j++) {
			getline(table, line);
			istringstream in(line);
			in >> ThermoTables_Density[i][j];
			in >> ThermoTables_Pressure[i][j];
			in >> ThermoTables_SoundSpeed2[i][j];
			in >> ThermoTables_Cp[i][j];
			in >> ThermoTables_Entropy[i][j];
			in >> ThermoTables_Mu[i][j];
			in >> ThermoTables_Kt[i][j];
			in >> ThermoTables_dPdrho_e[i][j];
			in >> ThermoTables_dPde_rho[i][j];
			in >> ThermoTables_dTdrho_e[i][j];
			in >> ThermoTables_dTde_rho[i][j];
			in >> ThermoTables_Temperature[i][j];
			in >> ThermoTables_StaticEnergy[i][j];
			in >> ThermoTables_Enthalpy[i][j];
		}
	}

	table.close();
	//NonDimensionalise and find limits
	NonDimensionalise_Table_Values();
	Find_Table_Limits();
}

void CLookUpTable::Find_Table_Limits() {
	Density_Table_Limits[0] = HUGE_VAL;
	Density_Table_Limits[1] = -HUGE_VAL;
	Pressure_Table_Limits[0] = HUGE_VAL;	//lower limit
	Pressure_Table_Limits[1] = -HUGE_VAL;	//upper limit
	Enthalpy_Table_Limits[0] = HUGE_VAL;	//lower limit
	Enthalpy_Table_Limits[1] = -HUGE_VAL;	//upper limit
	SoundSpeed2_Table_Limits[0] = HUGE_VAL;	//lower limit
	SoundSpeed2_Table_Limits[1] = -HUGE_VAL;	//upper limit
	Cp_Table_Limits[0] = HUGE_VAL;	//lower limit
	Cp_Table_Limits[1] = -HUGE_VAL;	//upper limit
	Entropy_Table_Limits[0] = HUGE_VAL;	//lower limit
	Entropy_Table_Limits[1] = -HUGE_VAL;	//upper limit
	dPdrho_e_Table_Limits[0] = HUGE_VAL;	//lower limit
	dPdrho_e_Table_Limits[1] = -HUGE_VAL;	//upper limit
	dPde_rho_Table_Limits[0] = HUGE_VAL;	//lower limit
	dPde_rho_Table_Limits[1] = -HUGE_VAL;	//upper limit
	dTdrho_e_Table_Limits[0] = HUGE_VAL;	//lower limit
	dTdrho_e_Table_Limits[1] = -HUGE_VAL;	//upper limit
	dTde_rho_Table_Limits[0] = HUGE_VAL;	//lower limit
	dTde_rho_Table_Limits[1] = -HUGE_VAL;	//upper limit
	Temperature_Table_Limits[0] = HUGE_VAL;	//lower limit
	Temperature_Table_Limits[1] = -HUGE_VAL;	//upper limit
	StaticEnergy_Table_Limits[0] = HUGE_VAL;	//lower limit
	StaticEnergy_Table_Limits[1] = -HUGE_VAL;	//upper limit
	for (int i = 0; i < Table_Density_Stations; i++) {
		for (int j = 0; j < Table_Pressure_Stations; j++) {
			//The table limits are stored for later checks of values becoming inconsistent
			if (ThermoTables_Density[i][j] > Density_Table_Limits[1]) {
				Density_Table_Limits[1] = ThermoTables_Density[i][j];
			}
			if (ThermoTables_Density[i][j] < Density_Table_Limits[0]) {
				Density_Table_Limits[0] = ThermoTables_Density[i][j];
			}

			if (ThermoTables_Pressure[i][j] > Pressure_Table_Limits[1]) {
				Pressure_Table_Limits[1] = ThermoTables_Pressure[i][j];
			}
			if (ThermoTables_Pressure[i][j] < Pressure_Table_Limits[0]) {
				Pressure_Table_Limits[0] = ThermoTables_Pressure[i][j];
			}
			if (ThermoTables_Enthalpy[i][j] > Enthalpy_Table_Limits[1]) {
				Enthalpy_Table_Limits[1] = ThermoTables_Enthalpy[i][j];
			}
			if (ThermoTables_Enthalpy[i][j] < Enthalpy_Table_Limits[0]) {
				Enthalpy_Table_Limits[0] = ThermoTables_Enthalpy[i][j];
			}
			if (ThermoTables_SoundSpeed2[i][j] > SoundSpeed2_Table_Limits[1]) {
				SoundSpeed2_Table_Limits[1] = ThermoTables_SoundSpeed2[i][j];
			}
			if (ThermoTables_SoundSpeed2[i][j] < SoundSpeed2_Table_Limits[0]) {
				SoundSpeed2_Table_Limits[0] = ThermoTables_SoundSpeed2[i][j];
			}
			if (ThermoTables_Cp[i][j] > Cp_Table_Limits[1]) {
				Cp_Table_Limits[1] = ThermoTables_Cp[i][j];
			}
			if (ThermoTables_Cp[i][j] < Cp_Table_Limits[0]) {
				Cp_Table_Limits[0] = ThermoTables_Cp[i][j];
			}
			if (ThermoTables_Entropy[i][j] > Entropy_Table_Limits[1]) {
				Entropy_Table_Limits[1] = ThermoTables_Entropy[i][j];
			}
			if (ThermoTables_Entropy[i][j] < Entropy_Table_Limits[0]) {
				Entropy_Table_Limits[0] = ThermoTables_Entropy[i][j];
			}
			if (ThermoTables_dPdrho_e[i][j] > dPdrho_e_Table_Limits[1]) {
				dPdrho_e_Table_Limits[1] = ThermoTables_dPdrho_e[i][j];
			}
			if (ThermoTables_dPdrho_e[i][j] < dPdrho_e_Table_Limits[0]) {
				dPdrho_e_Table_Limits[0] = ThermoTables_dPdrho_e[i][j];
			}
			if (ThermoTables_dPde_rho[i][j] > dPde_rho_Table_Limits[1]) {
				dPde_rho_Table_Limits[1] = ThermoTables_dPde_rho[i][j];
			}
			if (ThermoTables_dPde_rho[i][j] < dPde_rho_Table_Limits[0]) {
				dPde_rho_Table_Limits[0] = ThermoTables_dPde_rho[i][j];
			}
			if (ThermoTables_dTdrho_e[i][j] > dTdrho_e_Table_Limits[1]) {
				dTdrho_e_Table_Limits[1] = ThermoTables_dTdrho_e[i][j];
			}
			if (ThermoTables_dTdrho_e[i][j] < dTdrho_e_Table_Limits[0]) {
				dTdrho_e_Table_Limits[0] = ThermoTables_dTdrho_e[i][j];
			}
			if (ThermoTables_dTde_rho[i][j] > dTde_rho_Table_Limits[1]) {
				dTde_rho_Table_Limits[1] = ThermoTables_dTde_rho[i][j];
			}
			if (ThermoTables_dTde_rho[i][j] < dTde_rho_Table_Limits[0]) {
				dTde_rho_Table_Limits[0] = ThermoTables_dTde_rho[i][j];
			}
			if (ThermoTables_Temperature[i][j] > Temperature_Table_Limits[1]) {
				Temperature_Table_Limits[1] = ThermoTables_Temperature[i][j];
			}
			if (ThermoTables_Temperature[i][j] < Temperature_Table_Limits[0]) {
				Temperature_Table_Limits[0] = ThermoTables_Temperature[i][j];
			}
			if (ThermoTables_StaticEnergy[i][j] > StaticEnergy_Table_Limits[1]) {
				StaticEnergy_Table_Limits[1] = ThermoTables_StaticEnergy[i][j];
			}
			if (ThermoTables_StaticEnergy[i][j] < StaticEnergy_Table_Limits[0]) {
				StaticEnergy_Table_Limits[0] = ThermoTables_StaticEnergy[i][j];
			}
		}
	}
}
void CLookUpTable::NonDimensionalise_Table_Values() {
	for (int i = 0; i < Table_Density_Stations; i++) {
		for (int j = 0; j < Table_Pressure_Stations; j++) {
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

void CLookUpTable::LookUpTable_Malloc() {
	ThermoTables_StaticEnergy = new su2double*[Table_Density_Stations];
	ThermoTables_Entropy = new su2double*[Table_Density_Stations];
	ThermoTables_Enthalpy = new su2double*[Table_Density_Stations];
	ThermoTables_Density = new su2double*[Table_Density_Stations];
	ThermoTables_Pressure = new su2double*[Table_Density_Stations];
	ThermoTables_SoundSpeed2 = new su2double*[Table_Density_Stations];
	ThermoTables_Temperature = new su2double*[Table_Density_Stations];
	ThermoTables_dPdrho_e = new su2double*[Table_Density_Stations];
	ThermoTables_dPde_rho = new su2double*[Table_Density_Stations];
	ThermoTables_dTdrho_e = new su2double*[Table_Density_Stations];
	ThermoTables_dTde_rho = new su2double*[Table_Density_Stations];
	ThermoTables_Cp = new su2double*[Table_Density_Stations];
	ThermoTables_Mu = new su2double*[Table_Density_Stations];
	ThermoTables_dmudrho_T = new su2double*[Table_Density_Stations];
	ThermoTables_dmudT_rho = new su2double*[Table_Density_Stations];
	ThermoTables_Kt = new su2double*[Table_Density_Stations];
	ThermoTables_dktdrho_T = new su2double*[Table_Density_Stations];
	ThermoTables_dktdT_rho = new su2double*[Table_Density_Stations];
	for (int i = 0; i < Table_Density_Stations; i++) {
		ThermoTables_StaticEnergy[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Entropy[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Enthalpy[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Density[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Pressure[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_SoundSpeed2[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Temperature[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dPdrho_e[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dPde_rho[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dTdrho_e[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dTde_rho[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Cp[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Mu[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dmudrho_T[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dmudT_rho[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Kt[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dktdrho_T[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dktdT_rho[i] = new su2double[Table_Pressure_Stations];
	}

}

