/*!
 * transport_model.cpp
 * \brief Source of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

CViscosityModel::CViscosityModel(void) {

  /*--- Attributes initialization ---*/

  Mu = 0.0;
  dmudrho_T = 0.0;
  dmudT_rho = 0.0;

}

CViscosityModel::~CViscosityModel(void) { }


CConstantViscosity::CConstantViscosity(void) : CViscosityModel() { }

CConstantViscosity::CConstantViscosity(su2double mu_const) : CViscosityModel() {

  /*--- Attributes initialization ---*/

  Mu = mu_const;
  dmudrho_T = 0.0;
  dmudT_rho = 0.0;

}

CConstantViscosity::~CConstantViscosity(void) { }




CSutherland::CSutherland(void) : CViscosityModel() {
  Mu_ref = 0.0;
  T_ref = 0.0;
  S = 0.0;

}

CSutherland::CSutherland(su2double mu_ref, su2double t_ref, su2double s) : CViscosityModel() {

  Mu_ref = mu_ref;
  T_ref = t_ref;
  S = s;
}

CSutherland::~CSutherland(void) { }


void CSutherland::SetViscosity(su2double T, su2double rho) {

  Mu = Mu_ref*pow((T/T_ref),(3.0/2.0))*((T_ref + S)/(T + S));

}

void CSutherland::SetDerViscosity(su2double T, su2double rho) {

  dmudrho_T = 0.0;
  dmudT_rho = Mu_ref*( (3.0/2.0)*pow( (T/T_ref),(1.0/2.0) )*( (T_ref + S)/(T + S) )
          -pow( (T/T_ref),(3.0/2.0) )*(T_ref + S)/(T + S)/(T + S) );

}



CLookUpTable_Viscosity::CLookUpTable_Viscosity(CConfig *config, bool dimensional) : CViscosityModel() {
	LUT_Debug_Mode = config->GetLUT_Debug_Mode();
	rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	if (dimensional) {
		Pressure_Reference_Value = 1;
		Temperature_Reference_Value = 1;
		Density_Reference_Value = 1;
		Viscosity_Reference_Value=1;

	} else {
		Pressure_Reference_Value = config->GetPressure_Ref();
		Temperature_Reference_Value = config->GetTemperature_Ref();
		Density_Reference_Value = config->GetDensity_Ref();
		Viscosity_Reference_Value=1;

	}
	skewed_linear_table = false;
	//Detect cfx filetype
	if ((config->GetLUTFileName()).find(".rgp") != string::npos) {
		if (rank == MASTER_NODE) {
			cout << "CFX type LUT Viscosity found" << endl;
		}
		LookUpTable_Load_CFX(config->GetLUTFileName());
	}
	//Detect dat file type
	else if ((config->GetLUTFileName()).find(".dat") != string::npos) {
		if (rank == MASTER_NODE) {
			cout << "DAT type LUT Viscosity found" << endl;
		}
		LookUpTable_Load_DAT(config->GetLUTFileName());
	} else {
		if (rank == MASTER_NODE) {
			cout << "No recognized LUT Viscosity format found, exiting!" << endl;
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
		cout<<"Print LUT Viscosity errors? (LUT_Debug_Mode):  "<<LUT_Debug_Mode<<endl;
	}


}

CLookUpTable_Viscosity::~CLookUpTable_Viscosity(void) {
	// Delete the lut
	for (int i = 0; i < Table_Density_Stations; i++) {

		delete[] ThermoTables_Density[i];
		delete[] ThermoTables_Pressure[i];
		delete[] ThermoTables_Temperature[i];
		delete[] ThermoTables_Mu[i];
		delete[] ThermoTables_dmudrho_T[i];
		delete[] ThermoTables_dmudT_rho[i];
	}

	delete[] ThermoTables_Density;
	delete[] ThermoTables_Pressure;
	delete[] ThermoTables_Temperature;
	delete[] ThermoTables_Mu;
	delete[] ThermoTables_dmudrho_T;
	delete[] ThermoTables_dmudT_rho;


}

void CLookUpTable_Viscosity::Search_NonEquispaced_Rho_Index(su2double rho) {
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


void CLookUpTable_Viscosity::Search_j_for_Y_given_i(su2double x, su2double y,
		su2double **ThermoTables_X, su2double **ThermoTables_Y) {
	su2double RunVal;
	su2double grad, x00, y00, y10, x10;
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;

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


void CLookUpTable_Viscosity::SetViscosity(su2double T, su2double rho){
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "Mu RHOT Input Density out of bounds\n";
		}
		if ((T > Temperature_Table_Limits[1])
				or (T < Temperature_Table_Limits[0])) {
			cerr << "Mu RHOT Input Temperature out of bounds\n";
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
			ThermoTables_Temperature, "Mu RHOT");

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	Mu = Interpolate_2D_Bilinear(ThermoTables_Mu);

	//Check that the interpolated density and pressure are within LUT limits
	if (LUT_Debug_Mode)
	{
		Pressure = Interpolate_2D_Bilinear(ThermoTables_Pressure);
		Check_Interpolated_PRHO_Limits("dMu RHOT");
	}
}


void CLookUpTable_Viscosity::SetDerViscosity(su2double T, su2double rho){
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "Mu RHOT Input Density out of bounds\n";
		}
		if ((T > Temperature_Table_Limits[1])
				or (T < Temperature_Table_Limits[0])) {
			cerr << "Mu RHOT Input Temperature out of bounds\n";
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
			ThermoTables_Temperature, "Mu RHOT");

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	dmudrho_T = Interpolate_2D_Bilinear(ThermoTables_dmudrho_T);
	dmudT_rho= Interpolate_2D_Bilinear(ThermoTables_dmudT_rho);

	//Check that the interpolated density and pressure are within LUT limits
	if (LUT_Debug_Mode)
	{
		Pressure = Interpolate_2D_Bilinear(ThermoTables_Pressure);
		Check_Interpolated_PRHO_Limits("dMu RHOT");
	}

}

void CLookUpTable_Viscosity::Check_Interpolated_PRHO_Limits(string interpolation_case) {
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

inline void CLookUpTable_Viscosity::Gaussian_Inverse(int nDim) {
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

void CLookUpTable_Viscosity::Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(su2double x,
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

su2double CLookUpTable_Viscosity::Interpolate_2D_Bilinear(su2double * *ThermoTables_Z) {
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


void CLookUpTable_Viscosity::LookUpTable_Load_CFX(string filename) {
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
							ThermoTables_dmudrho_T[i][j] = 0; // (UNUSED)
							ThermoTables_dmudT_rho[i][j] = 0; // (UNUSED)
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
				//TEMPERATURE TABLE
				if (var == 15) {
					CFX_Import_Table_By_Number(&table, ThermoTables_Temperature, true);
				}

				//MU
				if (var == 8) {
					CFX_Import_Table_By_Number(&table, ThermoTables_Mu, true);
				}
			}
		}
	}
	table.close();
	//Non dimensionalise the table, and find the value limits
	NonDimensionalise_Table_Values();
	Find_Table_Limits();
}

void CLookUpTable_Viscosity::CFX_Import_Table_By_Number(ifstream *tab,
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

void CLookUpTable_Viscosity::LookUpTable_Load_DAT(std::string filename) {
	string line;
	string value;
	su2double dummy; //discarded values stored here

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
			in >> dummy;//ThermoTables_SoundSpeed2[i][j];
			in >> dummy;//ThermoTables_Cp[i][j];
			in >> dummy;//ThermoTables_Entropy[i][j];
			in >> ThermoTables_Mu[i][j];
			in >> dummy;//ThermoTables_Kt[i][j];
			in >> dummy;//ThermoTables_dPdrho_e[i][j];
			in >> dummy;//ThermoTables_dPde_rho[i][j];
			in >> dummy;//ThermoTables_dTdrho_e[i][j];
			in >> dummy;//ThermoTables_dTde_rho[i][j];
			in >> ThermoTables_Temperature[i][j];
			in >> dummy;//ThermoTables_StaticEnergy[i][j];
			in >> dummy;//ThermoTables_Enthalpy[i][j];
			ThermoTables_dmudrho_T[i][j] = 0; // (UNUSED)
			ThermoTables_dmudT_rho[i][j] = 0; // (UNUSED)

		}
	}

	table.close();
	//NonDimensionalise and find limits
	NonDimensionalise_Table_Values();
	Find_Table_Limits();
}

void CLookUpTable_Viscosity::Find_Table_Limits() {
	Density_Table_Limits[0] = HUGE_VAL;
	Density_Table_Limits[1] = -HUGE_VAL;
	Pressure_Table_Limits[0] = HUGE_VAL;	//lower limit
	Pressure_Table_Limits[1] = -HUGE_VAL;	//upper limit
	Temperature_Table_Limits[0] = HUGE_VAL;	//lower limit
	Temperature_Table_Limits[1] = -HUGE_VAL;	//upper limit;

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

			if (ThermoTables_Temperature[i][j] > Temperature_Table_Limits[1]) {
				Temperature_Table_Limits[1] = ThermoTables_Temperature[i][j];
			}
			if (ThermoTables_Temperature[i][j] < Temperature_Table_Limits[0]) {
				Temperature_Table_Limits[0] = ThermoTables_Temperature[i][j];
			}


		}
	}
}
void CLookUpTable_Viscosity::NonDimensionalise_Table_Values() {
	for (int i = 0; i < Table_Density_Stations; i++) {
		for (int j = 0; j < Table_Pressure_Stations; j++) {
			ThermoTables_Density[i][j] /= Density_Reference_Value;
			ThermoTables_Pressure[i][j] /= Pressure_Reference_Value;
			ThermoTables_Temperature[i][j] /= Temperature_Reference_Value;
			ThermoTables_Mu[i][j] /= Viscosity_Reference_Value;
			ThermoTables_dmudrho_T[i][j] = 0; // (UNUSED)
			ThermoTables_dmudT_rho[i][j] = 0; // (UNUSED)


		}
	}
}

void CLookUpTable_Viscosity::LookUpTable_Malloc() {
	ThermoTables_Density = new su2double*[Table_Density_Stations];
	ThermoTables_Pressure = new su2double*[Table_Density_Stations];
	ThermoTables_Temperature = new su2double*[Table_Density_Stations];
	ThermoTables_Mu = new su2double*[Table_Density_Stations];
	ThermoTables_dmudrho_T = new su2double*[Table_Density_Stations];
	ThermoTables_dmudT_rho = new su2double*[Table_Density_Stations];
	//ThermoTables_Kt = new su2double*[Table_Density_Stations];
	//ThermoTables_dktdrho_T = new su2double*[Table_Density_Stations];
	//ThermoTables_dktdT_rho = new su2double*[Table_Density_Stations];
	for (int i = 0; i < Table_Density_Stations; i++) {
		ThermoTables_Density[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Pressure[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Temperature[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Mu[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dmudrho_T[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dmudT_rho[i] = new su2double[Table_Pressure_Stations];
		//ThermoTables_Kt[i] = new su2double[Table_Pressure_Stations];
		//ThermoTables_dktdrho_T[i] = new su2double[Table_Pressure_Stations];
		//ThermoTables_dktdT_rho[i] = new su2double[Table_Pressure_Stations];
	}

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


CLookUpTable_Conductivity::CLookUpTable_Conductivity(CConfig *config) : CConductivityModel() {
	LUT_Debug_Mode = config->GetLUT_Debug_Mode();
	rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	Pressure_Reference_Value = config->GetPressure_Ref();
	Temperature_Reference_Value = config->GetTemperature_Ref();
	Density_Reference_Value = config->GetDensity_Ref();
	Conductivity_Reference_Value=config->GetConductivity_Ref();

	skewed_linear_table = false;
	//Detect cfx filetype
	if ((config->GetLUTFileName()).find(".rgp") != string::npos) {
		if (rank == MASTER_NODE) {
			cout << "CFX type LUT Conductivity found" << endl;
		}
		LookUpTable_Load_CFX(config->GetLUTFileName());
	}
	//Detect dat file type
	else if ((config->GetLUTFileName()).find(".dat") != string::npos) {
		if (rank == MASTER_NODE) {
			cout << "DAT type LUT Conductivity found" << endl;
		}
		LookUpTable_Load_DAT(config->GetLUTFileName());
	} else {
		if (rank == MASTER_NODE) {
			cout << "No recognized LUT Conductivity format found, exiting!" << endl;
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
		cout<<"Print LUT Conductivity errors? (LUT_Debug_Mode):  "<<LUT_Debug_Mode<<endl;
	}
}

CLookUpTable_Conductivity::~CLookUpTable_Conductivity(void) {
	// Delete the lut
	for (int i = 0; i < Table_Density_Stations; i++) {

		delete[] ThermoTables_Density[i];
		delete[] ThermoTables_Pressure[i];
		delete[] ThermoTables_Temperature[i];
		delete[] ThermoTables_Kt[i];
		delete[] ThermoTables_dktdrho_T[i];
		delete[] ThermoTables_dktdT_rho[i];
	}

	delete[] ThermoTables_Density;
	delete[] ThermoTables_Pressure;
	delete[] ThermoTables_Temperature;
	delete[] ThermoTables_Kt;
	delete[] ThermoTables_dktdrho_T;
	delete[] ThermoTables_dktdT_rho;


}

void CLookUpTable_Conductivity::Search_NonEquispaced_Rho_Index(su2double rho) {
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


void CLookUpTable_Conductivity::Search_j_for_Y_given_i(su2double x, su2double y,
		su2double **ThermoTables_X, su2double **ThermoTables_Y) {
	su2double RunVal;
	su2double grad, x00, y00, y10, x10;
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;

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


void CLookUpTable_Conductivity::SetConductivity(su2double T, su2double rho, su2double mu,
		su2double cp){
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "Kt RHOT Input Density out of bounds\n";
		}
		if ((T > Temperature_Table_Limits[1])
				or (T < Temperature_Table_Limits[0])) {
			cerr << "Kt RHOT Input Temperature out of bounds\n";
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
			ThermoTables_Temperature, "Kt RHOT");

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	Kt = Interpolate_2D_Bilinear(ThermoTables_Kt);

	//Check that the interpolated density and pressure are within LUT limits
	{
		Pressure = Interpolate_2D_Bilinear(ThermoTables_Pressure);
		Check_Interpolated_PRHO_Limits("Kt RHOT");
	}

}



void CLookUpTable_Conductivity::SetDerConductivity(su2double T, su2double rho, su2double mu,
		su2double cp){
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == MASTER_NODE and LUT_Debug_Mode) {
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "Kt RHOT Input Density out of bounds\n";
		}
		if ((T > Temperature_Table_Limits[1])
				or (T < Temperature_Table_Limits[0])) {
			cerr << "Kt RHOT Input Temperature out of bounds\n";
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
			ThermoTables_Temperature, "Kt RHOT");

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	dktdrho_T = Interpolate_2D_Bilinear(ThermoTables_dktdrho_T);
	dktdT_rho= Interpolate_2D_Bilinear(ThermoTables_dktdT_rho);

	//Check that the interpolated density and pressure are within LUT limits
	{
		Pressure = Interpolate_2D_Bilinear(ThermoTables_Pressure);
		Check_Interpolated_PRHO_Limits("dKt RHOT");
	}
}

void CLookUpTable_Conductivity::Check_Interpolated_PRHO_Limits(string interpolation_case) {
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

inline void CLookUpTable_Conductivity::Gaussian_Inverse(int nDim) {
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

void CLookUpTable_Conductivity::Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(su2double x,
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

su2double CLookUpTable_Conductivity::Interpolate_2D_Bilinear(su2double * *ThermoTables_Z) {
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


void CLookUpTable_Conductivity::LookUpTable_Load_CFX(string filename) {
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
							ThermoTables_dktdrho_T[i][j] = 0; // (UNUSED)
							ThermoTables_dktdT_rho[i][j] = 0; // (UNUSED)
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
				//TEMPERATURE TABLE
				if (var == 15) {
					CFX_Import_Table_By_Number(&table, ThermoTables_Temperature, true);
				}

				//Kt
				if (var == 9) {
					CFX_Import_Table_By_Number(&table, ThermoTables_Kt, true);
				}
			}
		}
	}
	table.close();
	//Non dimensionalise the table, and find the value limits
	NonDimensionalise_Table_Values();
	Find_Table_Limits();
}

void CLookUpTable_Conductivity::CFX_Import_Table_By_Number(ifstream *tab,
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

void CLookUpTable_Conductivity::LookUpTable_Load_DAT(std::string filename) {
	string line;
	string value;
	su2double dummy; //discarded values stored here

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
			in >> dummy;//ThermoTables_SoundSpeed2[i][j];
			in >> dummy;//ThermoTables_Cp[i][j];
			in >> dummy;//ThermoTables_Entropy[i][j];
			in >> dummy;//ThermoTables_Mu[i][j];
			in >> ThermoTables_Kt[i][j];
			in >> dummy;//ThermoTables_dPdrho_e[i][j];
			in >> dummy;//ThermoTables_dPde_rho[i][j];
			in >> dummy;//ThermoTables_dTdrho_e[i][j];
			in >> dummy;//ThermoTables_dTde_rho[i][j];
			in >> ThermoTables_Temperature[i][j];
			in >> dummy;//ThermoTables_StaticEnergy[i][j];
			in >> dummy;//ThermoTables_Enthalpy[i][j];
			ThermoTables_dktdrho_T[i][j] = 0; // (UNUSED)
			ThermoTables_dktdT_rho[i][j] = 0; // (UNUSED)

		}
	}

	table.close();
	//NonDimensionalise and find limits
	NonDimensionalise_Table_Values();
	Find_Table_Limits();
}

void CLookUpTable_Conductivity::Find_Table_Limits() {
	Density_Table_Limits[0] = HUGE_VAL;
	Density_Table_Limits[1] = -HUGE_VAL;
	Pressure_Table_Limits[0] = HUGE_VAL;	//lower limit
	Pressure_Table_Limits[1] = -HUGE_VAL;	//upper limit
	Temperature_Table_Limits[0] = HUGE_VAL;	//lower limit
	Temperature_Table_Limits[1] = -HUGE_VAL;	//upper limit;

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

			if (ThermoTables_Temperature[i][j] > Temperature_Table_Limits[1]) {
				Temperature_Table_Limits[1] = ThermoTables_Temperature[i][j];
			}
			if (ThermoTables_Temperature[i][j] < Temperature_Table_Limits[0]) {
				Temperature_Table_Limits[0] = ThermoTables_Temperature[i][j];
			}


		}
	}
}
void CLookUpTable_Conductivity::NonDimensionalise_Table_Values() {
	for (int i = 0; i < Table_Density_Stations; i++) {
		for (int j = 0; j < Table_Pressure_Stations; j++) {
			ThermoTables_Density[i][j] /= Density_Reference_Value;
			ThermoTables_Pressure[i][j] /= Pressure_Reference_Value;
			ThermoTables_Temperature[i][j] /= Temperature_Reference_Value;
			ThermoTables_Kt[i][j] /= Conductivity_Reference_Value;
			ThermoTables_dktdrho_T[i][j] = 0; // (UNUSED)
			ThermoTables_dktdT_rho[i][j] = 0; // (UNUSED)
		}
	}
}

void CLookUpTable_Conductivity::LookUpTable_Malloc() {
	ThermoTables_Density = new su2double*[Table_Density_Stations];
	ThermoTables_Pressure = new su2double*[Table_Density_Stations];
	ThermoTables_Temperature = new su2double*[Table_Density_Stations];
	ThermoTables_Kt = new su2double*[Table_Density_Stations];
	ThermoTables_dktdrho_T = new su2double*[Table_Density_Stations];
	ThermoTables_dktdT_rho = new su2double*[Table_Density_Stations];
	for (int i = 0; i < Table_Density_Stations; i++) {
		ThermoTables_Density[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Pressure[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Temperature[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_Kt[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dktdrho_T[i] = new su2double[Table_Pressure_Stations];
		ThermoTables_dktdT_rho[i] = new su2double[Table_Pressure_Stations];

	}

}

