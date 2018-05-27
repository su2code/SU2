/*!
 * \file data_manufactured_solutions.cpp
 * \brief Functions to compute the solution and source terms for manufactured solutions.
 * \author E. van der Weide
 * \version 6.0.1 "Falcon"
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

#include "../include/data_manufactured_solutions.hpp"
using namespace std;

#ifdef MANUFACTURED_SOLUTION

#ifdef MANUFACTURED_VISCOUS_UNIT_QUAD

/*--- Functions for the manufactured solution for the
      viscous unit quad. ---*/
void DetermineManufacturedSolution(const unsigned short nDim,
                                   const su2double      Gam,
                                   const su2double      RGas,
                                   const su2double      *coor,
                                         su2double      *sol) {

  /* Easier storage of the x- and y-coordinates. */
  const su2double x = coor[0];
  const su2double y = coor[1];

  /* Determine the solution for the density, velocity
     components and pressure. */
  const su2double LInv    = 1.0/L;
  const su2double PiLInv  = 3.1415926535897931*LInv;
  const su2double PiL2inv = PiLInv*LInv;

  const su2double rho = rho_0 + rho_x *sin(a_rhox *PiLInv*x)
                      +         rho_y *cos(a_rhoy *PiLInv*y)
                      +         rho_xy*cos(a_rhoxy*PiL2Inv*x*y);

  const su2double u = u_0 + u_x *sin(a_ux *PiLInv*x)
                    +       u_y *cos(a_uy *PiLInv*y)
                    +       u_xy*cos(a_uxy*PiL2Inv*x*y);

  const su2double v = v_0 + v_x *cos(a_vx *PiLInv*x)
                    +       v_y *sin(a_vy *PiLInv*y)
                    +       v_xy*cos(a_vxy*PiL2Inv*x*y);

  const su2double p = P_0 + P_x *cos(a_Px *PiLInv*x)
                    +       P_y *sin(a_Py *PiLInv*y)
                    +       P_xy*sin(a_Pxy*PiL2Inv*x*y);

  /* Compute the conservative variables from the primitive ones.
     Note that the implementation below is valid for both 2D and 3D. */
  sol[0]      = rho;
  sol[1]      = rho*u;
  sol[2]      = rho*v;
  sol[3]      = 0.0;
  sol[nDim+1] = p/(Gam-1.0) + 0.5*rho*(u*u + v*v);
}

void SourceTermManufacturedSolution(const unsigned short nDim,
                                    const su2double      Gam,
                                    const su2double      RGas,
                                    const su2double      mu,
                                    const su2double      k,
                                    const su2double      *coor,
                                          su2double      *source) {

  SU2_MPI::Error("Not implemented yet.", CURRENT_FUNCTION);
}

#else

/*--- No ifdef specified for the actual manufactured solution.
      Print an error message and exit. ---*/
void DetermineManufacturedSolution(const unsigned short nDim,
                                   const su2double      Gam,
                                   const su2double      RGas,
                                   const su2double      *coor,
                                         su2double      *sol) {

  SU2_MPI::Error("No or wrong compiler directive specified for the actual manufactured solution.",
                 CURRENT_FUNCTION);
}

void SourceTermManufacturedSolution(const unsigned short nDim,
                                    const su2double      Gam,
                                    const su2double      RGas,
                                    const su2double      mu,
                                    const su2double      k,
                                    const su2double      *coor,
                                          su2double      *source) {

  SU2_MPI::Error("No or wrong compiler directive specified for the actual manufactured solution.",
                 CURRENT_FUNCTION);
}

#endif

#endif // MANUFACTURED_SOLUTION
