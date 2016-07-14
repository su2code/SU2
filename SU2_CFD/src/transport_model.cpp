/*!
 * transport_model.cpp
 * \brief Source of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 4.2.0 "Cardinal"
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

#include "../include/transport_model.hpp"

/*-------------------------------------------------*/
/*----------- Dynamic Viscosity Models ------------*/
/*-------------------------------------------------*/

CViscosityModel::CViscosityModel(void) {

	/*--- Attributes initialization ---*/

	Mu = 0.0;
	dmudrho_T = 0.0;
	dmudT_rho = 0.0;

}

CViscosityModel::~CViscosityModel(void) {
}

CConstantViscosity::CConstantViscosity(void) :
				CViscosityModel() {
}

CConstantViscosity::CConstantViscosity(su2double mu_const) :
				CViscosityModel() {

	/*--- Attributes initialization ---*/

	Mu = mu_const;
	dmudrho_T = 0.0;
	dmudT_rho = 0.0;

}

CConstantViscosity::~CConstantViscosity(void) {
}

CSutherland::CSutherland(void) :
				CViscosityModel() {
	Mu_ref = 0.0;
	T_ref = 0.0;
	S = 0.0;

}

CSutherland::CSutherland(su2double mu_ref, su2double t_ref, su2double s) :
				CViscosityModel() {

	Mu_ref = mu_ref;
	T_ref = t_ref;
	S = s;
}

CSutherland::~CSutherland(void) {
}

void CSutherland::SetViscosity(su2double T, su2double rho) {

	Mu = Mu_ref * pow((T / T_ref), (3.0 / 2.0)) * ((T_ref + S) / (T + S));

}

void CSutherland::SetDerViscosity(su2double T, su2double rho) {

	dmudrho_T = 0.0;
	dmudT_rho = Mu_ref
			* ((3.0 / 2.0) * pow((T / T_ref), (1.0 / 2.0)) * ((T_ref + S) / (T + S))
					- pow((T / T_ref), (3.0 / 2.0)) * (T_ref + S) / (T + S) / (T + S));

}
//The structured look-up table approach

CViscosityList::CViscosityList() {

	Density = 0.0;
	Temperature = 0.0;
	Mu = 0.0;
	dmudrho_T = 0.0;
	dmudT_rho = 0.0;
}

CViscosityList::~CViscosityList() {

}

CLookUpTable_Viscosity::CLookUpTable_Viscosity(CConfig *config, bool dimensional) : CViscosityModel() {
	ViscosityTables = NULL;
	if (dimensional) {
		Temperature_Reference_Value = 1;
		Density_Reference_Value = 1;
		Viscosity_Reference_Value = 1;
	} else {
		Temperature_Reference_Value = config->GetTemperature_Ref();
		Density_Reference_Value = config->GetDensity_Ref();
		Viscosity_Reference_Value = config->GetViscosity_Ref();
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

}

CLookUpTable_Viscosity::~CLookUpTable_Viscosity(void) {
	for (int i = 0; i < Table_Density_Stations; i++) {
		delete[] ViscosityTables[i];
	}
	delete[] ViscosityTables;
}

void CLookUpTable_Viscosity::SetViscosity(su2double T, su2double rho) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
		cerr << "MU_RHOT Input Density out of bounds\n";
	}
	if ((T > Temperature_Table_Limits[1]) or (T < Temperature_Table_Limits[0])) {
		cerr << "MU_RHOT Input Temperature out of bounds\n";
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
		x00 = ViscosityTables[middleI][LowerJ].Density;
		grad = ViscosityTables[middleI + 1][LowerJ].Density - x00;
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
		y00 = ViscosityTables[LowerI][middleJ].Temperature;
		y10 = ViscosityTables[UpperI][middleJ].Temperature;
		x00 = ViscosityTables[LowerI][middleJ].Density;
		x10 = ViscosityTables[UpperI][middleJ].Density;
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (rho - x00);
		grad = ViscosityTables[LowerI][UpperJ].Temperature - y00;
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
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(x, y, "MU_RHOT");

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	Mu = Interpolate_2D_Bilinear("Mu");

	//Check that the interpolated density and pressure are within LUT limits
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		cerr << "MU_RHOT Interpolated Density out of bounds\n";
	}

	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;
}

void CLookUpTable_Viscosity::SetDerViscosity(su2double T, su2double rho) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
		cerr << "dMU_RHOT Input Density out of bounds\n";
	}
	if ((T > Temperature_Table_Limits[1]) or (T < Temperature_Table_Limits[0])) {
		cerr << "dMU_RHOT Input Temperature out of bounds\n";
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
		x00 = ViscosityTables[middleI][LowerJ].Density;
		grad = ViscosityTables[middleI + 1][LowerJ].Density - x00;
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
		y00 = ViscosityTables[LowerI][middleJ].Temperature;
		y10 = ViscosityTables[UpperI][middleJ].Temperature;
		x00 = ViscosityTables[LowerI][middleJ].Density;
		x10 = ViscosityTables[UpperI][middleJ].Density;
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (rho - x00);
		grad = ViscosityTables[LowerI][UpperJ].Temperature - y00;
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
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(x, y, "dMU_RHOT");

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	dmudrho_T = Interpolate_2D_Bilinear("dmudrho_T");
	dmudT_rho = Interpolate_2D_Bilinear("dmudT_rho");

	//Check that the interpolated density and pressure are within LUT limits
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		cerr << "RHOT Interpolated Density out of bounds\n";
	}

	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;
}
/*!
 * \brief Calculate the inverse of a square matrix (e.g. the Vandermonde matrix) with pivoting Gaussian elimination
 * \param[in] nDim - the dimension of the square block to invert
 */

inline void CLookUpTable_Viscosity::Gaussian_Inverse(int nDim) {
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
		max_val = abs(temp[k][k]);	//fixed bug for dimensionless normalized coords!
		//Find the largest value in the column
		for (int j = k; j < nDim; j++) {
			if (abs(temp[j][k]) > max_val) {
				max_idx = j;
				max_val = abs(temp[j][k]);//fixed bug for dimensionless normalized coords!
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

void CLookUpTable_Viscosity::Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(
		su2double x, su2double y, std::string grid_var) {
	//The x,y coordinates of the quadrilateral
	su2double x00, y00, x10, x01, x11, y10, y01, y11;
	su2double coords[8];

	//Load in the coordinates of the qudrilateral
	x00 =
			ViscosityTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Density;

	y00 =
			ViscosityTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Temperature;

	x01 =
			ViscosityTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Density;

	y01 =
			ViscosityTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Temperature;

	x10 =
			ViscosityTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Density;

	y10 =
			ViscosityTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Temperature;

	x11 =
			ViscosityTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Density;

	y11 =
			ViscosityTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Temperature;

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
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies below the LUT\n";
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
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies right of the LUT\n";
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
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					//	<< Nearest_Neighbour_jIndex[0]
					<< " interpolation point lies above the LUT\n";
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
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies left of the LUT\n";
		} else {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< +" interpolation point lies to the left of the boundary of selected quad\n";
		}
	}
	if (OUT_OF_BOUNDS) {
	} //extrapolate

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

su2double CLookUpTable_Viscosity::Interpolate_2D_Bilinear(
		string interpolant_var) {
	//The function values at the 4 corners of the quad
	su2double func_value_at_i0j0, func_value_at_i1j0, func_value_at_i0j1,
	func_value_at_i1j1;

	if (interpolant_var == "Mu") {
		func_value_at_i0j0 =
				ViscosityTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Mu;
		func_value_at_i1j0 =
				ViscosityTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Mu;
		func_value_at_i0j1 =
				ViscosityTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Mu;
		func_value_at_i1j1 =
				ViscosityTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Mu;
	} else if (interpolant_var == "dmudrho_T") {
		func_value_at_i0j0 =
				ViscosityTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dmudrho_T;
		func_value_at_i1j0 =
				ViscosityTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dmudrho_T;
		func_value_at_i0j1 =
				ViscosityTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dmudrho_T;
		func_value_at_i1j1 =
				ViscosityTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dmudrho_T;
	} else if (interpolant_var == "dmudT_rho") {
		func_value_at_i0j0 =
				ViscosityTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dmudT_rho;
		func_value_at_i1j0 =
				ViscosityTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dmudT_rho;
		func_value_at_i0j1 =
				ViscosityTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dmudT_rho;
		func_value_at_i1j1 =
				ViscosityTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dmudT_rho;
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

void CLookUpTable_Viscosity::Remove_Two_Phase_Region_CFX_Table(
		bool is_not_two_phase) {
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
	//					ViscosityTables[i + 1][j].Enthalpy - ViscosityTables[i][j].Enthalpy)
	//					> 0.1 * ViscosityTables[i + 1][j].Enthalpy) {
	//				Indexes_of_two_phase[i+1][j] = -10;
	//			}
	//		}
	//	}
	//	//Edge detection going up
	//	for (int j = 0; j < Table_Pressure_Stations; j++) {
	//		for (int i = Table_Density_Stations - 1; i > 0; i--) {
	//			if ((ViscosityTables[i][j].Enthalpy - ViscosityTables[i - 1][j].Enthalpy)
	//					> 1.1 * ViscosityTables[i - 1][j].Enthalpy) {
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

void CLookUpTable_Viscosity::LookUpTable_Load_CFX(string filename,
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
				//Create the actual LUT of CConductivityLists which is used in the FluidModel
				if (var == 1) {
					ViscosityTables = new CViscosityList*[x];
					for (int i = 0; i < x; i++) {
						ViscosityTables[i] = new CViscosityList[y];
					}
					//If table is to be later split up into 2phase region and superheated vapor
					//the saturation properties should also be captured.
					if (read_saturation_properties) {
						SaturationTables = new CViscosityList[y];
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
							ViscosityTables[i][j].Density = vD[i];
						}
					}
					delete[] vD;

				}
				// Check that all encountered tables adhere to the same x,y dimensions
				// otherwise throw an error
				else if (x != set_x && y != set_y) {
					cerr
					<< "The encountered dimensions of the CFX table are not the same throughout. They should be; for this to work.\n";

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
							ViscosityTables[i][j].Mu = inp[(j * set_x + i) % 10];
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
									inp[z] = inp[z] / (Viscosity_Reference_Value);
								}
							}
							SaturationTables[j].Mu = inp[j % 10];
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
							ViscosityTables[i][j].Temperature = inp[(j * set_x + i) % 10];
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
			}
		}
	}
	Table_Density_Stations = set_x;
	Table_Pressure_Stations = set_y;
	table.close();
}

/*-------------------------------------------------*/
/*---------- Thermal Conductivity Models ----------*/
/*-------------------------------------------------*/

CConductivityModel::CConductivityModel(void) {

	/*--- Attributes initialization ---*/

	Kt = 0.0;
	dktdrho_T = 0.0;
	dktdT_rho = 0.0;

}

CConductivityModel::~CConductivityModel(void) {
}

CConstantConductivity::CConstantConductivity(void) :
				CConductivityModel() {
}

CConstantConductivity::CConstantConductivity(su2double kt_const) :
				CConductivityModel() {

	/*--- Attributes initialization ---*/

	Kt = kt_const;
	dktdrho_T = 0.0;
	dktdT_rho = 0.0;

}

CConstantConductivity::~CConstantConductivity(void) {
}

CConstantPrandtl::CConstantPrandtl(void) :
				CConductivityModel() {
}

CConstantPrandtl::CConstantPrandtl(su2double pr_const) :
				CConductivityModel() {

	/*--- Attributes initialization ---*/

	Pr_const = pr_const;

}

void CConstantPrandtl::SetConductivity(su2double T, su2double rho, su2double mu,
		su2double cp) {

	Kt = mu * cp / Pr_const;

}

void CConstantPrandtl::SetDerConductivity(su2double T, su2double rho,
		su2double dmudrho_T, su2double dmudT_rho, su2double cp) {

	dktdrho_T = dmudrho_T * cp / Pr_const;
	dktdT_rho = dmudT_rho * cp / Pr_const;

}

CConstantPrandtl::~CConstantPrandtl(void) {
}

//The look up table approach
CConductivityList::CConductivityList() {

	Density = 0.0;
	Temperature = 0.0;
	Kt = 0.0;
	dktdrho_T = 0.0;
	dktdT_rho = 0.0;

}

CConductivityList::~CConductivityList() {

}

CLookUpTable_Conductivity::CLookUpTable_Conductivity(CConfig *config) :
		CConductivityModel() {
	ConductivityTables = NULL;

	Temperature_Reference_Value = config->GetTemperature_Ref();
	Density_Reference_Value = config->GetDensity_Ref();
	Conductivity_Reference_Value = config->GetConductivity_Ref();


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

}

CLookUpTable_Conductivity::~CLookUpTable_Conductivity(void) {
	// Delete the lut
	for (int i = 0; i < Table_Density_Stations; i++) {
		delete[] ConductivityTables[i];
	}
	delete[] ConductivityTables;

}

void CLookUpTable_Conductivity::SetDerConductivity(su2double T, su2double rho) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
		cerr << "dKt_RHOT Input Density out of bounds\n";
	}
	if ((T > Temperature_Table_Limits[1]) or (T < Temperature_Table_Limits[0])) {
		cerr << "dKt_RHOT Input Temperature out of bounds\n";
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
		x00 = ConductivityTables[middleI][LowerJ].Density;
		grad = ConductivityTables[middleI + 1][LowerJ].Density - x00;
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
		y00 = ConductivityTables[LowerI][middleJ].Temperature;
		y10 = ConductivityTables[UpperI][middleJ].Temperature;
		x00 = ConductivityTables[LowerI][middleJ].Density;
		x10 = ConductivityTables[UpperI][middleJ].Density;
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (rho - x00);
		grad = ConductivityTables[LowerI][UpperJ].Temperature - y00;
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
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(x, y, "dKt_RHOT");

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	dktdrho_T = Interpolate_2D_Bilinear("dktdrho_T");
	dktdT_rho = Interpolate_2D_Bilinear("dktdT_rho");

	//Check that the interpolated density and pressure are within LUT limits
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		cerr << "dKt_RHOT Interpolated Density out of bounds\n";
	}

	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;
}

void CLookUpTable_Conductivity::SetConductivity(su2double T, su2double rho) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
		cerr << "Kt_RHOT Input Density out of bounds\n";
	}
	if ((T > Temperature_Table_Limits[1]) or (T < Temperature_Table_Limits[0])) {
		cerr << "Kt_RHOT Input Temperature out of bounds\n";
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
		x00 = ConductivityTables[middleI][LowerJ].Density;
		grad = ConductivityTables[middleI + 1][LowerJ].Density - x00;
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
		y00 = ConductivityTables[LowerI][middleJ].Temperature;
		y10 = ConductivityTables[UpperI][middleJ].Temperature;
		x00 = ConductivityTables[LowerI][middleJ].Density;
		x10 = ConductivityTables[UpperI][middleJ].Density;
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (rho - x00);
		grad = ConductivityTables[LowerI][UpperJ].Temperature - y00;
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
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(x, y, "Kt_RHOT");

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	Kt = Interpolate_2D_Bilinear("Kt");

	//Check that the interpolated density and pressure are within LUT limits
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		cerr << "Kt_RHOT Interpolated Density out of bounds\n";
	}

	delete[] Nearest_Neighbour_iIndex;
	delete[] Nearest_Neighbour_jIndex;
}

inline void CLookUpTable_Conductivity::Gaussian_Inverse(int nDim) {
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
		max_val = abs(temp[k][k]);	//fixed bug for dimensionless normalized coords!
		//Find the largest value in the column
		for (int j = k; j < nDim; j++) {
			if (abs(temp[j][k]) > max_val) {
				max_idx = j;
				max_val = abs(temp[j][k]);//fixed bug for dimensionless normalized coords!
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

void CLookUpTable_Conductivity::Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(
		su2double x, su2double y, std::string grid_var) {
	//The x,y coordinates of the quadrilateral
	su2double x00, y00, x10, x01, x11, y10, y01, y11;

	//Load in the coordinates of the qudrilateral
	x00 =
			ConductivityTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Density;

	y00 =
			ConductivityTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Temperature;

	x01 =
			ConductivityTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Density;

	y01 =
			ConductivityTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Temperature;

	x10 =
			ConductivityTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Density;

	y10 =
			ConductivityTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Temperature;

	x11 =
			ConductivityTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Density;

	y11 =
			ConductivityTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Temperature;
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
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies below the LUT\n";
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
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies right of the LUT\n";
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
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					//	<< Nearest_Neighbour_jIndex[0]
					<< " interpolation point lies above the LUT\n";
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
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies left of the LUT\n";
		} else {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< +" interpolation point lies to the left of the boundary of selected quad\n";
		}
	}
	if (OUT_OF_BOUNDS) {}

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

su2double CLookUpTable_Conductivity::Interpolate_2D_Bilinear(
		string interpolant_var) {
	//The function values at the 4 corners of the quad
	su2double func_value_at_i0j0, func_value_at_i1j0, func_value_at_i0j1,
	func_value_at_i1j1;

	if (interpolant_var == "Kt") {
		func_value_at_i0j0 =
				ConductivityTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Kt;
		func_value_at_i1j0 =
				ConductivityTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Kt;
		func_value_at_i0j1 =
				ConductivityTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Kt;
		func_value_at_i1j1 =
				ConductivityTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Kt;
	} else if (interpolant_var == "dktdrho_T") {
		func_value_at_i0j0 =
				ConductivityTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dktdrho_T;
		func_value_at_i1j0 =
				ConductivityTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dktdrho_T;
		func_value_at_i0j1 =
				ConductivityTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dktdrho_T;
		func_value_at_i1j1 =
				ConductivityTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dktdrho_T;
	} else if (interpolant_var == "dktdT_rho") {
		func_value_at_i0j0 =
				ConductivityTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dktdT_rho;
		func_value_at_i1j0 =
				ConductivityTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dktdT_rho;
		func_value_at_i0j1 =
				ConductivityTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dktdT_rho;
		func_value_at_i1j1 =
				ConductivityTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dktdT_rho;
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

void CLookUpTable_Conductivity::Remove_Two_Phase_Region_CFX_Table(
		bool is_not_two_phase) {
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
	//					ConductivityTables[i + 1][j].Enthalpy - ConductivityTables[i][j].Enthalpy)
	//					> 0.1 * ConductivityTables[i + 1][j].Enthalpy) {
	//				Indexes_of_two_phase[i+1][j] = -10;
	//			}
	//		}
	//	}
	//	//Edge detection going up
	//	for (int j = 0; j < Table_Pressure_Stations; j++) {
	//		for (int i = Table_Density_Stations - 1; i > 0; i--) {
	//			if ((ConductivityTables[i][j].Enthalpy - ConductivityTables[i - 1][j].Enthalpy)
	//					> 1.1 * ConductivityTables[i - 1][j].Enthalpy) {
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

void CLookUpTable_Conductivity::LookUpTable_Load_CFX(string filename,
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
				//Create the actual LUT of CConductivityLists which is used in the FluidModel
				if (var == 1) {
					ConductivityTables = new CConductivityList*[x];
					for (int i = 0; i < x; i++) {
						ConductivityTables[i] = new CConductivityList[y];
					}
					//If table is to be later split up into 2phase region and superheated vapor
					//the saturation properties should also be captured.
					if (read_saturation_properties) {
						SaturationTables = new CConductivityList[y];
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
							ConductivityTables[i][j].Density = vD[i];
						}
					}
					delete[] vD;
				}
				// Check that all encountered tables adhere to the same x,y dimensions
				// otherwise throw an error
				else if (x != set_x && y != set_y) {
					cerr
					<< "The encountered dimensions of the CFX table are not the same throughout. They should be; for this to work.\n";

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
									 inp[z] /=Conductivity_Reference_Value;
								}
							}
							ConductivityTables[i][j].Kt = inp[(j * set_x + i) % 10];
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
			}
		}
	}
	Table_Density_Stations = set_x;
	Table_Pressure_Stations = set_y;
	table.close();
}

