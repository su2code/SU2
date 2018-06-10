/*!
 * \file wall_model.cpp
 * \brief File, which contains the implementation for the wall model functions
 *        for large eddy simulations.
 * \author E. van der Weide, T. Economon, P. Urbanczyk
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
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

#include "../include/wall_model.hpp"

void CWallModel1DEQ::Initialize(const unsigned short *intInfo,
                                const su2double      *doubleInfo){

  /* Copy the data from the arguments into the member variables. */
  numPoints      = intInfo[0];
  thickness      = doubleInfo[0];
  expansionRatio = doubleInfo[1];

  /* Allocate the memory for the coordinates of the grid points used
     in the 1D equilibrium wall model. */
  coorGridPoints.resize(numPoints);

  /* Determine the scaled version of the normal coordinates, where the
     first normal coordinate is simply 1.0. */
  su2double currentHeight = 1.0;
  coorGridPoints[0] = 0.0;

  for(unsigned short i=1; i<numPoints; ++i) {
    coorGridPoints[i] = coorGridPoints[i-1] + currentHeight;
    currentHeight    *= expansionRatio;
  }

  /* Determine the scaling factor of the normal coordinates and
     apply the scaling to obtain the correct coordinates. */
  const su2double scaleFact = thickness/coorGridPoints[numPoints-1];

  for(unsigned short i=0; i<numPoints; ++i)
    coorGridPoints[i] *= scaleFact;
}

void CWallModel1DEQ::WallShearStressAndHeatFlux(const su2double rhoExchange,
                                                const su2double velExchange,
                                                const su2double pExchange,
                                                const su2double Wall_HeatFlux,
                                                const bool      HeatFlux_Prescribed,
                                                const su2double Wall_Temperature,
                                                const bool      Temperature_Prescribed,
                                                      su2double &tauWall,
                                                      su2double &qWall,
                                                      su2double &ViscosityWall,
                                                      su2double &kOverCvWall) {

  cout << "CWallModel1DEQ::WallShearStressAndHeatFlux: Not implemented yet" << endl;
  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}
