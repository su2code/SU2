/*!
 * \file data_manufactured_solutions.hpp
 * \brief Parameters for the manufactured solutions.
 * \author E. van der Weide
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../Common/include/mpi_structure.hpp"
#include <iostream>

#ifdef MANUFACTURED_SOLUTION

/*!
 * \brief Function, which determines the prescribed solution for the given
          input data.
 * \param[in]  nDim - Number of spatial dimensions.
 * \param[in]  Gam  - Specific heat ratio
 * \param[in]  RGas - Gas constant.
 * \param[in]  coor - Coordinates where the solution must be determined.
 * \param[out] sol  - Solution to be computed.
 */
void DetermineManufacturedSolution(const unsigned short nDim,
                                   const su2double      Gam,
                                   const su2double      RGas,
                                   const su2double      *coor,
                                         su2double      *sol);

/*!
 * \brief Function, which determines the source term of the manufactured
          solution for the given input data.
 * \param[in]  nDim   - Number of spatial dimensions.
 * \param[in]  Gam    - Specific heat ratio
 * \param[in]  RGas   - Gas constant.
 * \param[in]  coor   - Coordinates where the source term must be computed.
 * \param[out] source - Source term to be computed.
 */
void SourceTermManufacturedSolution(const unsigned short nDim,
                                    const su2double      Gam,
                                    const su2double      RGas,
                                    const su2double      mu,
                                    const su2double      k,
                                    const su2double      *coor,
                                          su2double      *source);

#ifdef MANUFACTURED_VISCOUS_UNIT_QUAD

/*--- Constants, which describe this manufactured solution. ---*/
/*--- This is a viscous solution on the unit quad, where    ---*/
/*--- the primitive variables vary as a combination of      ---*/
/*--- sine and cosine functions.                            ---*/
const su2double L       =      1.0;
const su2double a_Px    =      1.0;
const su2double a_Pxy   =      0.75;
const su2double a_Py    =      1.25;
const su2double a_rhox  =      0.75;
const su2double a_rhoxy =      1.25;
const su2double a_rhoy  =      1.0;
const su2double a_ux    =      1.6666666667;
const su2double a_uxy   =      0.6;
const su2double a_uy    =      1.5;
const su2double a_vx    =      1.5;
const su2double a_vxy   =      0.9;
const su2double a_vy    =      1.0;
const su2double P_0     = 100000.0;
const su2double P_x     = -30000.0;
const su2double P_xy    = -25000.0;
const su2double P_y     =  20000.0;
const su2double rho_0   =      1.0;
const su2double rho_x   =      0.1;
const su2double rho_xy  =      0.08;
const su2double rho_y   =      0.15;
const su2double u_0     =     70.0;
const su2double u_x     =      4.0;
const su2double u_xy    =      7.0;
const su2double u_y     =    -12.0;
const su2double v_0     =     90.0;
const su2double v_x     =    -20.0;
const su2double v_xy    =    -11.0;
const su2double v_y     =      4.0;

#endif

#endif  // MANUFACTURED_SOLUTION
