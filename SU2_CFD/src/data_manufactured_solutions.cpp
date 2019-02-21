/*!
 * \file data_manufactured_solutions.cpp
 * \brief Functions to compute the solution and source terms for manufactured solutions.
 * \author E. van der Weide
 * \version 6.2.0 "Falcon"
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
  const su2double t2 = 0.1e1 / L;
  const su2double t3 = t2 * Pi;
  const su2double t4 = a_rhox * Pi;
  const su2double t5 = x * t2;
  const su2double t6 = t5 * t4;
  const su2double t7 = cos(t6);
  const su2double t10 = a_rhoxy * rho_xy;
  const su2double t11 = Pi * t10;
  const su2double t12 = L * L;
  const su2double t13 = 0.1e1 / t12;
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
  const su2double t48 = cos(t24);
  const su2double t52 = u_xy * a_uxy * Pi;
  const su2double t53 = sin(t33);
  const su2double t56 = t48 * t3 * u_x * a_ux - t53 * t14 * t52;
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
  const su2double t80 = cos(t70);
  const su2double t84 = v_xy * a_vxy * Pi;
  const su2double t85 = sin(t74);
  const su2double t88 = t80 * t3 * v_y * a_vy - t85 * t16 * t84;
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
  const su2double t109 = a_uxy * a_uxy;
  const su2double t111 = t13 * Pi * t109;
  const su2double t112 = y * y;
  const su2double t116 = a_vxy * a_vxy;
  const su2double t118 = t13 * Pi * t116;
  const su2double t120 = v_xy * x;
  const su2double t122 = t120 * t75 * y * t118;
  const su2double t124 = a_vxy * t85;
  const su2double t125 = v_xy * t124;
  const su2double t127 = a_ux * a_ux;
  const su2double t131 = (-u_xy * t34 * t112 * t111 + t122 / 0.2e1 + t125 / 0.2e1 - t26 * Pi * t127) * Pi;
  const su2double t132 = t13 * mu;
  const su2double t137 = u_y * a_uy;
  const su2double t138 = sin(t29);
  const su2double t143 = -t138 * t3 * t137 - t53 * t16 * t52;
  const su2double t147 = a_uy * a_uy;
  const su2double t150 = x * x;
  const su2double t161 = (-t13 * (u_xy * t34 * t150 * t111 + t31 * Pi * t147) * Pi - t13 * (t122 + t125) * Pi) * mu;
  const su2double t165 = v_x * a_vx;
  const su2double t166 = sin(t66);
  const su2double t171 = -t85 * t14 * t84 - t166 * t3 * t165;
  const su2double t174 = u_xy * x;
  const su2double t176 = t174 * t34 * y * t111;
  const su2double t177 = a_uxy * t53;
  const su2double t178 = u_xy * t177;
  const su2double t182 = a_vx * a_vx;
  const su2double t192 = (-t13 * (t176 + t178) * Pi - t13 * (v_xy * t75 * t112 * t118 + t68 * Pi * t182) * Pi) * mu;
  const su2double t193 = t77 * t77;
  const su2double t200 = t28 * a_Py * Pi;
  const su2double t201 = cos(t200);
  const su2double t203 = t201 * t3 * P_y * a_Py;
  const su2double t204 = t106 * t16;
  const su2double t205 = t204 * t103;
  const su2double t210 = a_vy * a_vy;
  const su2double t215 = (0.2e1 * v_xy * t75 * t150 * t118 + 0.2e1 * t72 * Pi * t210 - t176 - t178) * Pi;
  const su2double t219 = -t101 + t108;
  const su2double t221 = 1 / (Gam - 1);
  const su2double t223 = t91 + t193;
  const su2double t233 = cos(t98);
  const su2double t234 = t233 * P_x;
  const su2double t235 = sin(t200);
  const su2double t236 = t235 * P_y;
  const su2double t237 = sin(t105);
  const su2double t238 = t237 * P_xy;
  const su2double t239 = P_0 + t234 + t236 + t238;
  const su2double t243 =  t221 * t239 + t223 * t46 / 0.2e1 + P_0 + t234 + t236 + t238;
  const su2double t245 = a_Px * a_Px;
  const su2double t246 = Pi * t245;
  const su2double t248 = a_Pxy * a_Pxy;
  const su2double t250 = P_xy * t112 * t248;
  const su2double t251 = t13 * Pi;
  const su2double t252 = t237 * t251;
  const su2double t258 = P_x * t99;
  const su2double t260 = a_Pxy * y;
  const su2double t261 = t106 * P_xy;
  const su2double t264 = (t258 * L * a_Px - t261 * t260) * rho_xy;
  const su2double t267 = t19 * y * t251;
  const su2double t272 = a_rhoxy * a_rhoxy;
  const su2double t273 = t272 * rho_xy;
  const su2double t277 = t44 * t13 * Pi * t239;
  const su2double t285 = t39 + t43 + rho_0;
  const su2double t288 = t237 * t13 * Pi * t285;
  const su2double t290 = a_rhox * a_rhox;
  const su2double t291 = t290 * rho_x;
  const su2double t296 = P_xy * t7;
  const su2double t306 = t234 + t236 + P_0;
  const su2double t308 = t38 * t3;
  const su2double t311 = t39 + rho_0;
  const su2double t318 = t44 * t44;
  const su2double t319 = rho_xy * rho_xy;
  const su2double t321 = t285 * rho_xy;
  const su2double t324 = t42 * t42;
  const su2double t325 = rho_y * rho_y;
  const su2double t327 = t311 * rho_y;
  const su2double t330 = t7 * t7;
  const su2double t331 = rho_x * rho_x;
  const su2double t336 = rho_0 * rho_0;
  const su2double t337 = 0.2e1 * rho_x * rho_0 * t38 + t319 * t318 + 0.2e1 * t44 * t321 + t325 * t324 + 0.2e1 * t42 * t327 - t331 * t330 + t331 + t336;
  const su2double t338 = 0.1e1 / t337;
  const su2double t340 = 0.1e1 / RGas;
  const su2double t341 = t340 * t13;
  const su2double t349 = t106 * t285 * P_xy;
  const su2double t365 = t337 * t337;
  const su2double t366 = 0.1e1 / t365;
  const su2double t369 = a_rhoxy * t319 * t44;
  const su2double t376 = a_rhoxy * t321;
  const su2double t398 = u_xy * y * t177;
  const su2double t399 = t120 * t124;
  const su2double t402 = a_ux * u_x * t48;
  const su2double t404 = a_vy * v_y * t80;
  const su2double t428 = (-t13 * (t137 * t138 * L + t174 * t177) * Pi - t13 * (t165 * t166 * L + v_xy * y * t124) * Pi) * mu;
  const su2double t430 = t203 + t205;
  const su2double t442 = a_Py * a_Py;
  const su2double t446 = P_xy * t150 * t248;
  const su2double t452 = P_y * t201;
  const su2double t454 = a_Pxy * x;
  const su2double t457 = (t452 * L * a_Py + t261 * t454) * rho_xy;
  const su2double t460 = t19 * x * t251;
  const su2double t474 = a_rhoy * a_rhoy;
  const su2double t475 = t474 * rho_y;
  const su2double t480 = P_xy * t59;
  const su2double t486 = t235 * t2;
  const su2double t564 = t36 * ( t221 * t219 + t223 * t22 / 0.2e1 + (t171 * t77 + t56 * t36) * t46 - t101 + t108) + t56 * t243
                       + k * t341 * t338 * (t44 * (t234 * t246 + t252 * t250) * rho_xy - t267 * a_rhoxy * t264
                       - t19 * t219 * y * t10 - t277 * t112 * t273 - t106 * t7 * t2 * t4 * P_xy * rho_x * t260 + t288 * t250 
                       + (t42 * P_x * t233 * t2 * Pi * rho_y * t245 - t238 * t38 * t2 * Pi * t291 
                       + t107 * t104 * t296 * t1 + t311 * t234 * t2 * t246 - t308 * t306 * t291) * L) * Pi 
                       - 0.2e1 * (t44 * t7 * t3 * rho_xy * rho_x * a_rhox + t42 * t7 * t3 * rho_y * rho_x * a_rhox 
                       + rho_x * rho_0 * t7 * t2 * t4 + t308 * a_rhox * t331 * t7 - t267 * t369 - t267 * t376) * k * t341 * t366 
                       * (t44 * t264 - t19 * t239 * y * t10 - t349 * t260 + (t311 * P_x * t99 * a_Px + t42 * t258 * a_Px * rho_y 
                       + t237 * t296 * t1 + t7 * t306 * t1) * L) * Pi - fourThird * t36 * t132 * t131 
                       - fourThird * t56 * t132 * (-t398 + t399 / 0.2e1 + (t402 - t404 / 0.2e1) * L) * Pi 
                       - t77 * t192 - t171 * t428 + t77 * ( t221 * t430 + t223 * t64 / 0.2e1 + (t143 * t36 + t88 * t77) * t46 + t203 + t205)
                       + t88 * t243 - k * t340 * t338 * t13 * (t44 * (-t236 * Pi * t442 - t252 * t446) * rho_xy 
                       - t460 * a_rhoxy * t457 + t19 * t430 * x * t10 + t277 * t150 * t273 - t106 * t59 * t2 * t40 * P_xy * rho_y * t454 
                       - t288 * t446 + (-P_y * t38 * t486 * Pi * rho_x * t442 - t42 * P_y * t486 * Pi * rho_y * t442
                       - P_y * t486 * Pi * rho_0 * t442 + t238 * t42 * t2 * Pi * t475 + t204 * t104 * t480 * t58 
                       + t42 * t3 * t306 * t475) * L) * Pi + 0.2e1 * (-t44 * t59 * t3 * rho_xy * rho_y * a_rhoy 
                       - t60 * a_rhoy * t325 * t42 - t60 * a_rhoy * t327 - t460 * t369 - t460 * t376) * k * t340 * t366 * t13 
                       * (t44 * t457 + t19 * t239 * x * t10 + t349 * t454 + (P_y * t38 * t201 * a_Py * rho_x + t42 * t452 * a_Py * rho_y 
                       + t452 * a_Py * rho_0 + t237 * t480 * t58 + t59 * t306 * t58) * L) * Pi - t36 * t161 
                       - t143 * t428 + 0.2e1 / 0.3e1 * t77 * t132 * t215 + 0.2e1 / 0.3e1 * t88 * t132 * (-t398 + 0.2e1 * t399 
                       + (t402 - 0.2e1 * t404) * L) * Pi;

  source[0]      = t88 * t46 + t77 * t64 + t37 + t57;
  source[1]      = t91 * t22 + 0.2e1 * t56 * t93 - t101 + t108 - fourThird * t132 * t131 + t77 * t36 * t64 + t77 * t143 * t46 + t88 * t93 - t161;
  source[2]      = t77 * t37 + t77 * t57 + t171 * t93 - t192 + t193 * t64 + 0.2e1 * t88 * t77 * t46 + t203 + t205 + 0.2e1 / 0.3e1 * t132 * t215;
  source[3]      = 0.0;
  source[nDim+1] = t564;


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
