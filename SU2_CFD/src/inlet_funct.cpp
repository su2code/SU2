/*!
 * \file inlet_funct.cpp
 * \brief Functions to specify the inlet profile.
 * \author D. Manosalvas
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

#include "../include/inlet_funct.hpp"

/*!
 * \brief Transform the real location value to a number between -1 and 1.
 * \return A value between -1 and 1.
 */

su2double ScaleCoordinate(su2double y_max, su2double y_min, su2double y){
	su2double z = (y-y_min)/(y_max-y_min)*2 - 1;
	return z;
}

/*!
 * \brief Evaluates a polynomial thst represent the inlet velocity profile.
 * \return Normalized velocity profile from 0 to 1.
 */
su2double poly2D(su2double C1, su2double C2, su2double C3, su2double C4, su2double C5, su2double y){
	su2double Vel = C1*y*y*y*y + C2*y*y*y + C3*y*y + C4*y + C5;
	return Vel;
}


/*!
 * \brief Evaluates a piecewise velocity profile and includes the amplitude.
 * \return Velocity profile with amplitude A.
 */
su2double polydisc(su2double A , su2double y_max, su2double y_min, su2double y){
	su2double rho = 1.217; // Average jet density obtained from 2D @ T=290 P=101325
	su2double mu = 1.79820992909e-05; // Viscosity @ T=290
	su2double W = 0.168; // Truck width
	su2double Vel = 0;
	
	su2double Re = rho*A*W/mu; // Calculates the Reynolds Number to show jet development
	su2double d = 0.382*W/pow(Re,0.2); //Turbulent BL thickness for a flow that travel W
	
	su2double y_mid = y_min+(y_max-y_min)/2; // Calculates the middle of the real range
	
	// The real coordinates are scaled to meet the needs of the function to be between -1 and 1
	su2double z = ScaleCoordinate(y_max, y_min, y);
	su2double dz = ScaleCoordinate(y_max, y_min, y_mid + d);
	su2double hz = 2*dz;
	
	// Piecewise function
	if (z <= -1 + dz){
		Vel = A*( 1 - 4/pow(hz,2)*pow((z + 1 - dz),2));
	} else if (z >= 1 - dz){
		Vel = A*( 1 - 4/pow(hz,2)*pow((z - 1 + dz),2));
	} else {
		Vel = A;
	}
	return Vel;
}
