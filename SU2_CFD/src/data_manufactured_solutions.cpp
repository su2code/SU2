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
  const su2double PiL2Inv = PiLInv*LInv;

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

  /*--- The source code for the source terms is generated in Maple. ---*/
  const su2double Pi = 3.1415926535897931;
  const su2double fourThird = 4.0/3.0;
  const su2double x = coor[0];
  const su2double y = coor[1];
  
  const su2double t1 = rho_x * a_rhox;
  const su2double t2 = 1.0 / L;
  const su2double t3 = t2 * Pi;
  const su2double t4 = a_rhox * Pi;
  const su2double t5 = x * t2;
  const su2double t6 = t5 * t4;
  const su2double t7 = cos(t6);
  const su2double t10 = a_rhoxy * rho_xy;
  const su2double t11 = Pi * t10;
  const su2double t12 = L * L;
  const su2double t13 = 1.0 / t12;
  const su2double t14 = y * t13;
  const su2double t16 = x * t13;
  const su2double t17 = y * t16;
  const su2double t18 = t17 * a_rhoxy * Pi;
  const su2double t19 = sin(t18);
  const su2double t22 = t7 * t3 * t1 - t19 * t14 * t11;
  const su2double t24 = t5 * a_ux * Pi;
  const su2double t25 = sin(t24);
  const su2double t26 = t25 * u_x;
  const su2double t28 = y * t2;
  const su2double t29 = t28 * a_uy * Pi;
  const su2double t30 = cos(t29);
  const su2double t31 = t30 * u_y;
  const su2double t33 = t17 * a_uxy * Pi;
  const su2double t34 = cos(t33);
  const su2double t36 = t34 * u_xy + t26 + t31 + u_0;
  const su2double t37 = t36 * t22;
  const su2double t38 = sin(t6);
  const su2double t39 = t38 * rho_x;
  const su2double t40 = a_rhoy * Pi;
  const su2double t41 = t28 * t40;
  const su2double t42 = cos(t41);
  const su2double t43 = t42 * rho_y;
  const su2double t44 = cos(t18);
  const su2double t46 = t44 * rho_xy + rho_0 + t39 + t43;
  const su2double t47 = u_x * a_ux;
  const su2double t48 = cos(t24);
  const su2double t52 = u_xy * a_uxy * Pi;
  const su2double t53 = sin(t33);
  const su2double t56 = -t53 * t14 * t52 + t48 * t3 * t47;
  const su2double t57 = t56 * t46;
  const su2double t58 = rho_y * a_rhoy;
  const su2double t59 = sin(t41);
  const su2double t60 = t59 * t3;
  const su2double t64 = -t19 * t16 * t11 - t60 * t58;
  const su2double t66 = t5 * a_vx * Pi;
  const su2double t67 = cos(t66);
  const su2double t68 = t67 * v_x;
  const su2double t70 = t28 * a_vy * Pi;
  const su2double t71 = sin(t70);
  const su2double t72 = t71 * v_y;
  const su2double t74 = t17 * a_vxy * Pi;
  const su2double t75 = cos(t74);
  const su2double t77 = t75 * v_xy + t68 + t72 + v_0;
  const su2double t79 = v_y * a_vy;
  const su2double t80 = cos(t70);
  const su2double t84 = v_xy * a_vxy * Pi;
  const su2double t85 = sin(t74);
  const su2double t88 = -t85 * t16 * t84 + t80 * t3 * t79;
  const su2double t91 = t36 * t36;
  const su2double t93 = t36 * t46;
  const su2double t98 = t5 * a_Px * Pi;
  const su2double t99 = sin(t98);
  const su2double t101 = t99 * t3 * P_x * a_Px;
  const su2double t103 = P_xy * a_Pxy * Pi;
  const su2double t104 = a_Pxy * Pi;
  const su2double t105 = t17 * t104;
  const su2double t106 = cos(t105);
  const su2double t107 = t106 * t14;
  const su2double t108 = t107 * t103;
  const su2double t109 = mu * Pi;
  const su2double t110 = a_ux * a_ux;
  const su2double t113 = a_uxy * a_uxy;
  const su2double t115 = t13 * Pi * t113;
  const su2double t116 = y * y;
  const su2double t121 = t13 * (-u_xy * t34 * t116 * t115 - t26 * Pi * t110);
  const su2double t126 = u_y * a_uy;
  const su2double t127 = sin(t29);
  const su2double t132 = -t127 * t3 * t126 - t53 * t16 * t52;
  const su2double t136 = a_uy * a_uy;
  const su2double t139 = x * x;
  const su2double t146 = a_vxy * a_vxy;
  const su2double t148 = t13 * Pi * t146;
  const su2double t150 = v_xy * y;
  const su2double t153 = a_vxy * t85;
  const su2double t159 = (-t13 * (u_xy * t34 * t139 * t115 + t31 * Pi * t136) * Pi - t13 * (t150 * t75 * x * t148 + v_xy * t153) * Pi) * mu;
  const su2double t163 = v_x * a_vx;
  const su2double t164 = sin(t66);
  const su2double t169 = -t85 * t14 * t84 - t164 * t3 * t163;
  const su2double t172 = u_xy * x;
  const su2double t175 = a_uxy * t53;
  const su2double t180 = a_vx * a_vx;
  const su2double t190 = (-t13 * (t172 * t34 * y * t115 + u_xy * t175) * Pi - t13 * (v_xy * t75 * t116 * t148 + t68 * Pi * t180) * Pi) * mu;
  const su2double t191 = t77 * t77;
  const su2double t198 = t28 * a_Py * Pi;
  const su2double t199 = cos(t198);
  const su2double t201 = t199 * t3 * P_y * a_Py;
  const su2double t202 = t106 * t16;
  const su2double t203 = t202 * t103;
  const su2double t204 = a_vy * a_vy;
  const su2double t211 = t13 * (-v_xy * t75 * t139 * t148 - t72 * Pi * t204);
  const su2double t215 = -t101 + t108;
  const su2double t217 = 1.0/(Gam - 1.0);
  const su2double t219 = t91 + t191;
  const su2double t229 = cos(t98);
  const su2double t230 = t229 * P_x;
  const su2double t231 = sin(t198);
  const su2double t232 = t231 * P_y;
  const su2double t233 = sin(t105);
  const su2double t234 = t233 * P_xy;
  const su2double t235 = P_0 + t230 + t232 + t234;
  const su2double t239 = t217 * t235 + t219 * t46 / 2.0 + P_0 + t230 + t232 + t234;
  const su2double t241 = a_Px * a_Px;
  const su2double t242 = Pi * t241;
  const su2double t244 = a_Pxy * a_Pxy;
  const su2double t246 = P_xy * t116 * t244;
  const su2double t247 = t13 * Pi;
  const su2double t248 = t233 * t247;
  const su2double t254 = P_x * t99;
  const su2double t256 = a_Pxy * y;
  const su2double t257 = t106 * P_xy;
  const su2double t260 = (t254 * L * a_Px - t257 * t256) * rho_xy;
  const su2double t263 = t19 * y * t247;
  const su2double t268 = a_rhoxy * a_rhoxy;
  const su2double t269 = t268 * rho_xy;
  const su2double t273 = t44 * t13 * Pi * t235;
  const su2double t281 = t39 + t43 + rho_0;
  const su2double t284 = t233 * t13 * Pi * t281;
  const su2double t286 = a_rhox * a_rhox;
  const su2double t287 = t286 * rho_x;
  const su2double t292 = P_xy * t7;
  const su2double t302 = t230 + t232 + P_0;
  const su2double t304 = t38 * t3;
  const su2double t307 = t39 + rho_0;
  const su2double t315 = t44 * t44;
  const su2double t316 = rho_xy * rho_xy;
  const su2double t318 = t281 * rho_xy;
  const su2double t321 = t42 * t42;
  const su2double t322 = rho_y * rho_y;
  const su2double t324 = t307 * rho_y;
  const su2double t327 = t7 * t7;
  const su2double t328 = rho_x * rho_x;
  const su2double t333 = rho_0 * rho_0;
  const su2double t334 = 2.0 * rho_x * rho_0 * t38 + t316 * t315 + 2.0 * t44 * t318 + t322 * t321 + 2.0 * t42 * t324 - t328 * t327 + t328 + t333;
  const su2double t336 = 1.0 / RGas;
  const su2double t338 = k * t336 / t334;
  const su2double t345 = t106 * t281 * P_xy;
  const su2double t362 = t334 * t334;
  const su2double t364 = t336 / t362;
  const su2double t366 = a_rhoxy * t316 * t44;
  const su2double t373 = a_rhoxy * t318;
  const su2double t417 = (-t13 * (t126 * t127 * L + t172 * t175) * Pi - t13 * (t163 * t164 * L + t150 * t153) * Pi) * mu;
  const su2double t419 = t201 + t203;
  const su2double t431 = a_Py * a_Py;
  const su2double t435 = P_xy * t139 * t244;
  const su2double t441 = P_y * t199;
  const su2double t443 = a_Pxy * x;
  const su2double t446 = (t441 * L * a_Py + t257 * t443) * rho_xy;
  const su2double t449 = t19 * x * t247;
  const su2double t463 = a_rhoy * a_rhoy;
  const su2double t464 = t463 * rho_y;
  const su2double t469 = P_xy * t59;
  const su2double t475 = t231 * t2;
  const su2double t550 = t36 * (t217 * t215 + t219 * t22 / 2.0 + (t169 * t77 + t56 * t36) * t46 - t101 + t108)
                       + t56 * t239 + t338 * t13 * Pi * (t44 * (t230 * t242 + t248 * t246) * rho_xy - t263 * a_rhoxy * t260
                       - t19 * t215 * y * t10 - t273 * t116 * t269 - t106 * t7 * t2 * t4 * P_xy * rho_x * t256 + t284 * t246
                       + (t42 * P_x * t229 * t2 * Pi * rho_y * t241 - t234 * t38 * t2 * Pi * t287 + t107 * t104 * t292 * t1
                       + t307 * t230 * t2 * t242 - t304 * t302 * t287) * L) - 2.0 * (t44 * t7 * t3 * rho_xy * rho_x * a_rhox
                       + t42 * t7 * t3 * rho_y * rho_x * a_rhox + rho_x * rho_0 * t7 * t2 * t4 + t304 * a_rhox * t328 * t7
                       - t263 * t366 - t263 * t373) * k * t364 * t13 * Pi * (t44 * t260 - t19 * t235 * y * t10 - t345 * t256
                       + (t307 * P_x * t99 * a_Px + t42 * t254 * a_Px * rho_y + t233 * t292 * t1 + t7 * t302 * t1) * L)
                       - fourThird * t36 * t121 * t109 - fourThird * t56 * t13 * (t47 * t48 * L - u_xy * y * t175) * t109
                       - t77 * t190 - t169 * t417 + t77 * (t217 * t419 + t219 * t64 / 2.0 + (t132 * t36 + t88 * t77) * t46 + t201 + t203)
                       + t88 * t239 - t338 * t13 * Pi * (t44 * (-t232 * Pi * t431 - t248 * t435) * rho_xy - t449 * a_rhoxy * t446
                       + t19 * t419 * x * t10 + t273 * t139 * t269 - t106 * t59 * t2 * t40 * P_xy * rho_y * t443 - t284 * t435
                       + (-P_y * t38 * t475 * Pi * rho_x * t431 - t42 * P_y * t475 * Pi * rho_y * t431 - P_y * t475 * Pi * rho_0 * t431
                       + t234 * t42 * t2 * Pi * t464 + t202 * t104 * t469 * t58 + t42 * t3 * t302 * t464) * L)
                       + 2.0 * (-t44 * t59 * t3 * rho_xy * rho_y * a_rhoy - t60 * a_rhoy * t322 * t42 - t60 * a_rhoy * t324 - t449 * t366
                       - t449 * t373) * k * t364 * t13 * Pi * (t44 * t446 + t19 * t235 * x * t10 + t345 * t443
                       + (P_y * t38 * t199 * a_Py * rho_x + t42 * t441 * a_Py * rho_y + t441 * a_Py * rho_0 + t233 * t469 * t58
                       + t59 * t302 * t58) * L) - t36 * t159 - t132 * t417 - fourThird * t77 * t211 * t109
                       - fourThird * t88 * t13 * (t79 * t80 * L - v_xy * x * t153) * t109;

  source[0] = t88 * t46 + t77 * t64 + t37 + t57;
  source[1] = t91 * t22 + 2.0 * t56 * t93 - t101 + t108 - fourThird * t121 * t109 + t77 * t36 * t64 + t77 * t132 * t46 + t88 * t93 - t159;
  source[2] = t77 * t37 + t77 * t57 + t169 * t93 - t190 + t191 * t64 + 2.0 * t88 * t77 * t46 + t201 + t203 - fourThird * t211 * t109;
  source[3] = 0.0;
  source[nDim+1] = t550;
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
