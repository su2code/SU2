/*!
 * \file CUserDefinedSolution.cpp
 * \brief Implementations of the member functions of CUserDefinedSolution.
 * \author T. Economon, E. van der Weide
 * \version 7.5.1 "Blackbird"
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

#include "../../../include/toolboxes/MMS/CUserDefinedSolution.hpp"

CUserDefinedSolution::CUserDefinedSolution() : CVerificationSolution() {}

CUserDefinedSolution::CUserDefinedSolution(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_iMesh,
                                           CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  /*--- Write a message that the solution is initialized for a
   user-defined verification case. ---*/

  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for a user-defined verification case!!!" << endl;
    cout << endl << flush;
  }

  gamma = config->GetGamma();

//  const su2double val = 0.5;  /** for uniform mesh 1x1x1 **/
  const su2double val = 27.5;   /** for stretched mesh 1x1x0.01 **/

  r0 = 1.0;
  r1 = 0.35;
  r2 = 0.525;
  r3 = 0.55;
  r4 = val;
  
  u0 = 0.1;
  u1 = 0.10;
  u2 = 0.125;
  u3 = 0.15;
  u4 = val;
  
  v0 = 0.2;
  v1 = 0.2;
  v2 = 0.225;
  v3 = 0.25;
  v4 = val;
  
  w0 = 0.3;
  w1 = 0.3;
  w2 = 0.325;
  w3 = 0.35;
  w4 = val;
  
  p0 = 0.4;
  p1 = 0.4;
  p2 = 0.425;
  p3 = 0.45;
  p4 = val;

  if (nDim == 2) {

    r4 = u4 = v4 = p4 = 0.0;

    w0 = 0.0;
    w1 = 0.0;
    w2 = 0.0;
    w3 = 0.0;
    w4 = 0.0;

  }

}

CUserDefinedSolution::~CUserDefinedSolution() = default;

void CUserDefinedSolution::GetBCState(const su2double* val_coords, const su2double val_t,
                                      su2double* val_solution) const {
  GetSolution(val_coords,val_t,val_solution);
}

void CUserDefinedSolution::GetSolution(const su2double* val_coords, const su2double val_t,
                                       su2double* val_solution) const {
  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  const su2double z = (nDim == 3) ? val_coords[2] : 0.0;

  /* Determine the solution for the density, velocity
     components and pressure. */

  const su2double rho = r1*exp(r2*x + r3*y + r4*z) + r0;

  const su2double u = u1*exp(u2*x + u3*y + u4*z) + u0;

  const su2double v = v1*exp(v2*x + v3*y + v4*z) + v0;

  const su2double w = (nDim == 3) ? w1*exp(w2*x + w3*y + w4*z) + w0 : 0.0;

  const su2double p = p1*exp(p2*x + p3*y + p4*z) + p0;


  /* Compute the conservative variables from the primitive ones.
     Note that the implementation below is valid for both 2D and 3D. */
  val_solution[0] = rho;
  val_solution[1] = rho * u;
  val_solution[2] = rho * v;
  val_solution[3] = rho * w;
  val_solution[nDim+1] = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v + w * w);
}

void CUserDefinedSolution::GetMMSSourceTerm(const su2double* val_coords, const su2double val_t,
                                            su2double* val_source) const {

  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  const su2double z = (nDim == 3) ? val_coords[2] : 0.0;

  val_source[0] = exp(u2*x + u3*y + u4*z)*r0*u1*u2 + 
                  exp(r2*x + u2*x + r3*y + u3*y + r4*z + u4*z)*r1*u1*(r2 + u2) +
                  exp(v2*x + v3*y + v4*z)*r0*v1*v3 + 
                  exp(r2*x + v2*x + r3*y + v3*y + r4*z + v4*z)*r1*v1*(r3 + v3) +
                  exp(r2*x + r3*y + r4*z)*r1*(r2*u0 + r3*v0 + r4*w0) +
                  exp(w2*x + w3*y + w4*z)*r0*w1*w4 + 
                  exp(r2*x + w2*x + r3*y + w3*y + r4*z + w4*z)*r1*w1*(r4 + w4);

  
  val_source[1] = exp(p2*x + p3*y + p4*z)*p1*p2 +
                  exp(r2*x + r3*y + r4*z)*r1*r2*pow(u0 + exp(u2*x + u3*y + u4*z)*u1,2) +
                  2*exp(u2*x + u3*y + u4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*u1*(u0 + exp(u2*x + u3*y + u4*z)*u1)*u2 +
                  exp(r2*x + r3*y + r4*z)*r1*r3*(u0 + exp(u2*x + u3*y + u4*z)*u1)*(v0 + exp(v2*x + v3*y + v4*z)*v1) +
                  exp(u2*x + u3*y + u4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*u1*u3*(v0 + exp(v2*x + v3*y + v4*z)*v1) +
                  exp(v2*x + v3*y + v4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*(u0 + exp(u2*x + u3*y + u4*z)*u1)*v1*v3 +
                  exp(r2*x + r3*y + r4*z)*r1*r4*(u0 + exp(u2*x + u3*y + u4*z)*u1)*(w0 + exp(w2*x + w3*y + w4*z)*w1) +
                  exp(u2*x + u3*y + u4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*u1*u4*(w0 + exp(w2*x + w3*y + w4*z)*w1) +
                  exp(w2*x + w3*y + w4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*(u0 + exp(u2*x + u3*y + u4*z)*u1)*w1*w4;
  
  val_source[2] = exp(p2*x + p3*y + p4*z)*p1*p3 + 
                  exp(r2*x + r3*y + r4*z)*r1*r2*(u0 + exp(u2*x + u3*y + u4*z)*u1)*(v0 + exp(v2*x + v3*y + v4*z)*v1) +
                  exp(u2*x + u3*y + u4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*u1*u2*(v0 + exp(v2*x + v3*y + v4*z)*v1) +
                  exp(r2*x + r3*y + r4*z)*r1*r3*pow(v0 + exp(v2*x + v3*y + v4*z)*v1,2) +
                  exp(v2*x + v3*y + v4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*(u0 + exp(u2*x + u3*y + u4*z)*u1)*v1*v2 +
                  2*exp(v2*x + v3*y + v4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*v1*(v0 + exp(v2*x + v3*y + v4*z)*v1)*v3 +
                  exp(r2*x + r3*y + r4*z)*r1*r4*(v0 + exp(v2*x + v3*y + v4*z)*v1)*(w0 + exp(w2*x + w3*y + w4*z)*w1) +
                  exp(v2*x + v3*y + v4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*v1*v4*(w0 + exp(w2*x + w3*y + w4*z)*w1) +
                  exp(w2*x + w3*y + w4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*(v0 + exp(v2*x + v3*y + v4*z)*v1)*w1*w4;


  val_source[3] = exp(p2*x + p3*y + p4*z)*p1*p4 + 
                  exp(r2*x + r3*y + r4*z)*r1*r2*(u0 + exp(u2*x + u3*y + u4*z)*u1)*(w0 + exp(w2*x + w3*y + w4*z)*w1) +
                  exp(u2*x + u3*y + u4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*u1*u2*(w0 + exp(w2*x + w3*y + w4*z)*w1) +
                  exp(r2*x + r3*y + r4*z)*r1*r3*(v0 + exp(v2*x + v3*y + v4*z)*v1)*(w0 + exp(w2*x + w3*y + w4*z)*w1) +
                  exp(v2*x + v3*y + v4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*v1*v3*(w0 + exp(w2*x + w3*y + w4*z)*w1) +
                  exp(r2*x + r3*y + r4*z)*r1*r4*pow(w0 + exp(w2*x + w3*y + w4*z)*w1,2) +
                  exp(w2*x + w3*y + w4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*(u0 + exp(u2*x + u3*y + u4*z)*u1)*w1*w2 +
                  exp(w2*x + w3*y + w4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*(v0 + exp(v2*x + v3*y + v4*z)*v1)*w1*w3 +
                  2*exp(w2*x + w3*y + w4*z)*(r0 + exp(r2*x + r3*y + r4*z)*r1)*w1*(w0 + exp(w2*x + w3*y + w4*z)*w1)*w4;


  val_source[nDim+1] = exp(u2*x + u3*y + u4*z)*u1*u2*
                             ((gamma*(p0 + exp(p2*x + p3*y + p4*z)*p1))/
                                  (-1 + gamma) + 0.5*
                                  (r0 + exp(r2*x + r3*y + r4*z)*r1)*
                                  (pow(u0 + exp(u2*x + u3*y + u4*z)*u1,
                                         2) + pow(v0 + 
                                             exp(v2*x + v3*y + v4*z)*v1,2) + 
                                   pow(w0 + exp(w2*x + w3*y + w4*z)*w1,2))
                                  ) + exp(v2*x + v3*y + v4*z)*v1*v3*
                             ((gamma*(p0 + exp(p2*x + p3*y + p4*z)*p1))/
                                  (-1 + gamma) + 0.5*
                                  (r0 + exp(r2*x + r3*y + r4*z)*r1)*
                                  (pow(u0 + exp(u2*x + u3*y + u4*z)*u1,
                                         2) + pow(v0 + 
                                             exp(v2*x + v3*y + v4*z)*v1,2) + 
                                   pow(w0 + exp(w2*x + w3*y + w4*z)*w1,2))
                                  ) + (u0 + exp(u2*x + u3*y + u4*z)*u1)*
                             ((exp(p2*x + p3*y + p4*z)*gamma*p1*p2)/(-1 + gamma) + 
                              0.5*exp(r2*x + r3*y + r4*z)*r1*r2*
                                  (pow(u0 + exp(u2*x + u3*y + u4*z)*u1,
                                         2) + pow(v0 + 
                                             exp(v2*x + v3*y + v4*z)*v1,2) + 
                                   pow(w0 + exp(w2*x + w3*y + w4*z)*w1,2))
                              + 1.*(r0 + exp(r2*x + r3*y + r4*z)*r1)*
                                    (exp(u2*x + u3*y + u4*z)*u1*
                                         (u0 + exp(u2*x + u3*y + u4*z)*u1)*u2 + 
                                     exp(v2*x + v3*y + v4*z)*v1*
                                         (v0 + exp(v2*x + v3*y + v4*z)*v1)*v2 + 
                                     exp(w2*x + w3*y + w4*z)*w1*
                                         (w0 + exp(w2*x + w3*y + w4*z)*w1)*w2))\
                         + (v0 + exp(v2*x + v3*y + v4*z)*v1)*
                               ((exp(p2*x + p3*y + p4*z)*gamma*p1*p3)/(-1 + gamma) + 
                                0.5*exp(r2*x + r3*y + r4*z)*r1*r3*
                                    (pow(u0 + exp(u2*x + u3*y + u4*z)*u1,
                                           2) + pow(v0 + 
                                               exp(v2*x + v3*y + v4*z)*v1,2) + 
                                     pow(w0 + exp(w2*x + w3*y + w4*z)*w1,2))
                                + 1.*(r0 + exp(r2*x + r3*y + r4*z)*r1)*
                                      (exp(u2*x + u3*y + u4*z)*u1*
                                           (u0 + exp(u2*x + u3*y + u4*z)*u1)*u3 + 
                                       exp(v2*x + v3*y + v4*z)*v1*
                                           (v0 + exp(v2*x + v3*y + v4*z)*v1)*v3 + 
                                       exp(w2*x + w3*y + w4*z)*w1*
                                           (w0 + exp(w2*x + w3*y + w4*z)*w1)*w3))\
                         + exp(w2*x + w3*y + w4*z)*w1*
                               ((gamma*(p0 + exp(p2*x + p3*y + p4*z)*p1))/
                                    (-1 + gamma) + 0.5*
                                    (r0 + exp(r2*x + r3*y + r4*z)*r1)*
                                    (pow(u0 + exp(u2*x + u3*y + u4*z)*u1,
                                           2) + pow(v0 + 
                                               exp(v2*x + v3*y + v4*z)*v1,2) + 
                                     pow(w0 + exp(w2*x + w3*y + w4*z)*w1,2))
                                    )*w4 + (w0 + exp(w2*x + w3*y + w4*z)*w1)*
                             ((exp(p2*x + p3*y + p4*z)*gamma*p1*p4)/(-1 + gamma) + 
                              0.5*exp(r2*x + r3*y + r4*z)*r1*r4*
                                  (pow(u0 + exp(u2*x + u3*y + u4*z)*u1,
                                         2) + pow(v0 + 
                                             exp(v2*x + v3*y + v4*z)*v1,2) + 
                                   pow(w0 + exp(w2*x + w3*y + w4*z)*w1,2))
                              + 1.*(r0 + exp(r2*x + r3*y + r4*z)*r1)*
                                    (exp(u2*x + u3*y + u4*z)*u1*
                                         (u0 + exp(u2*x + u3*y + u4*z)*u1)*u4 + 
                                     exp(v2*x + v3*y + v4*z)*v1*
                                         (v0 + exp(v2*x + v3*y + v4*z)*v1)*v4 + 
                                     exp(w2*x + w3*y + w4*z)*w1*
                                         (w0 + exp(w2*x + w3*y + w4*z)*w1)*w4));

}

bool CUserDefinedSolution::IsManufacturedSolution() const {
  return true;
}
