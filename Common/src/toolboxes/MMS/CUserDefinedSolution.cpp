/*!
 * \file CUserDefinedSolution.cpp
 * \brief Implementations of the member functions of CUserDefinedSolution.
 * \author T. Economon, E. van der Weide
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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
                                           CConfig *config)
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
    Pr = config->GetPrandtl_Lam();
    R = config->GetGas_ConstantND();
    C1 = config->GetMu_RefND() * (config->GetMu_Temperature_RefND() + config->GetMu_SND()) / pow(config->GetMu_Temperature_RefND(), 1.5);
    S = config->GetMu_SND();


////   const su2double val = 0.475;  /** for uniform mesh 1x1x1 **/
////    const su2double val = 5;   /** for stretched mesh 1x1x0.01 **/
//    const su2double val = 27.5;   /** for stretched mesh 1x1x0.01 **/
//
//    r0 = 1.0;
//    r1 = 0.35;
//    r2 = 0.525;
//    r3 = 0.55;
//    r4 = val;
//
//    u0 = 0.1;
//    u1 = 0.10;
//    u2 = 0.125;
//    u3 = 0.15;
//    u4 = val;
//
//    v0 = 0.2;
//    v1 = 0.2;
//    v2 = 0.225;
//    v3 = 0.25;
//    v4 = val;
//
//    w0 = 0.3;
//    w1 = 0.3;
//    w2 = 0.325;
//    w3 = 0.35;
//    w4 = val;
//
//    p0 = 0.4;
//    p1 = 0.35;
//    p2 = 0.425;
//    p3 = 0.45;
//    p4 = val;
//
//    if (nDim == 2) {
//
//        r4 = u4 = v4 = p4 = 0.0;
//
//        w0 = 0.0;
//        w1 = 0.0;
//        w2 = 0.0;
//        w3 = 0.0;
//        w4 = 0.0;
//
//    }

//    const su2double stretch = 1.0;
    const su2double stretch = 100.0;

    apx = 1.0;
    apy = 1.25;
    apz = 1.5 * stretch;
    apxyz = 0.75 * stretch;

    arx = 0.75;
    ary = 1.0;
    arz = 1.25 * stretch;
    arxyz = 1.25 * stretch;

    aux = 1.6666666667;
    auy = 1.5;
    auz = 1.4 * stretch;
    auxyz = 0.6 * stretch;

    avx = 1.5;
    avy = 1.0;
    avz = 0.5 * stretch;
    avxyz = 0.9 * stretch;

    awx = 1.4;
    awy = 1.1;
    awz = 0.7 * stretch;
    awxyz = 0.8 * stretch;

    p0 = 100000.0 / config->GetPressure_Ref();
    px = -30000.0 / config->GetPressure_Ref();
    py = 20000.0 / config->GetPressure_Ref();
    pz = 10000.0 / config->GetPressure_Ref();
    pxyz = -25000.0 / config->GetPressure_Ref();

    r0 = 1.0 / config->GetDensity_Ref();
    rx = 0.1 / config->GetDensity_Ref();
    ry = 0.15 / config->GetDensity_Ref();
    rz = 0.16 / config->GetDensity_Ref();
    rxyz = 0.08 / config->GetDensity_Ref();


    u0 = 140.0 / config->GetVelocity_Ref();
    ux = 4.0 / config->GetVelocity_Ref();
    uy = -12.0 / config->GetVelocity_Ref();
    uz = -6.0 / config->GetVelocity_Ref();
    uxyz = 7.0 / config->GetVelocity_Ref();

    v0 = 180.0 / config->GetVelocity_Ref();
    vx = -20.0 / config->GetVelocity_Ref();
    vy = 4.0 / config->GetVelocity_Ref();
    vz = -8.0 / config->GetVelocity_Ref();
    vxyz = -11.0 / config->GetVelocity_Ref();

    w0 = 160.0 / config->GetVelocity_Ref();
    wx = 20.0 / config->GetVelocity_Ref();
    wy = -10.0 / config->GetVelocity_Ref();
    wz = 15.0 / config->GetVelocity_Ref();
    wxyz = 11.0 / config->GetVelocity_Ref();

    if (nDim == 2) {
        rz = uz = vz = wz = pz = 0;
        arz = auz = avz = awz = apz = 0;
    }


}

CUserDefinedSolution::~CUserDefinedSolution() = default;

void CUserDefinedSolution::GetBCState(const su2double *val_coords, const su2double val_t,
                                      su2double *val_solution) const {
    GetSolution(val_coords, val_t, val_solution);
}

void CUserDefinedSolution::GetSolution(const su2double *val_coords, const su2double val_t,
                                       su2double *val_solution) const {
    /* Easier storage of the x- and y-coordinates. */
    const su2double x = val_coords[0];
    const su2double y = val_coords[1];
    const su2double z = (nDim == 3) ? val_coords[2] : 1.0;

    const su2double Pi = M_PI;

    /* Determine the solution for the density, velocity
     components and pressure. */

//    const su2double rho = r1 * exp(r2 * x + r3 * y + r4 * z) + r0;
//
//    const su2double u = u1 * exp(u2 * x + u3 * y + u4 * z) + u0;
//
//    const su2double v = v1 * exp(v2 * x + v3 * y + v4 * z) + v0;
//
//    const su2double w = (nDim == 3) ? w1 * exp(w2 * x + w3 * y + w4 * z) + w0 : 0.0;
//
//    const su2double p = p1 * exp(p2 * x + p3 * y + p4 * z) + p0;


    const su2double L = 1;

    const su2double rho = r0 + rx*sin(arx*Pi*x/L) + ry*cos(ary*Pi*y/L) + rz*cos(arz*Pi*z/L) + rxyz*cos(arxyz*Pi*x*y*z/pow(L,3));

    const su2double u = u0 + ux*sin(aux*Pi*x/L) + uy*cos(auy*Pi*y/L) + uz*sin(auz*Pi*z/L) + uxyz*cos(auxyz*Pi*x*y*z/pow(L,3));

    const su2double v = v0 + vx*cos(avx*Pi*x/L) + vy*sin(avy*Pi*y/L) + vz*sin(avz*Pi*z/L) + vxyz*cos(avxyz*Pi*x*y*z/pow(L,3));

    const su2double w = (nDim == 3) ? w0 + wx*sin(awx*Pi*x/L) + wy*cos(awy*Pi*y/L) + wz*cos(awz*Pi*z/L) + wxyz*sin(awxyz*Pi*x*y*z/pow(L,3)) : 0.0;

    const su2double p = p0 + px*cos(apx*Pi*x/L) + py*sin(apy*Pi*y/L) + pz*cos(apz*Pi*z/L) + pxyz*sin(apxyz*Pi*x*y*z/pow(L,3));


    /* Compute the conservative variables from the primitive ones.
     Note that the implementation below is valid for both 2D and 3D. */
    val_solution[0] = rho;
    val_solution[1] = rho * u;
    val_solution[2] = rho * v;
    val_solution[3] = rho * w;
    val_solution[nDim + 1] = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v + w * w);
}

void CUserDefinedSolution::GetMMSSourceTerm(const su2double *val_coords, const su2double val_t,
                                            su2double *val_source) const {

    /* Easier storage of the x- and y-coordinates. */
    const su2double x = val_coords[0];
    const su2double y = val_coords[1];
    const su2double z = (nDim == 3) ? val_coords[2] : 1.0;

    const su2double g = gamma;

    const su2double L = 1;
    const su2double Pi = M_PI;

    val_source[0] = (r0 + ry*cos((ary*Pi*y)/L) + rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    ((awxyz*Pi*wxyz*x*y*
                      cos((awxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3) -
                     (awz*Pi*wz*sin((awz*Pi*z)/L))/L) +
                    (v0 + vx*cos((avx*Pi*x)/L) +
                     vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                     vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L))*
                    (-((ary*Pi*ry*sin((ary*Pi*y)/L))/L) -
                     (arxyz*Pi*rxyz*x*z*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    (u0 + uy*cos((auy*Pi*y)/L) +
                     uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                     ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L))*
                    ((arx*Pi*rx*cos((arx*Pi*x)/L))/L -
                     (arxyz*Pi*rxyz*y*z*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    ((aux*Pi*ux*cos((aux*Pi*x)/L))/L -
                     (auxyz*Pi*uxyz*y*z*
                      sin((auxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    ((avy*Pi*vy*cos((avy*Pi*y)/L))/L -
                     (avxyz*Pi*vxyz*x*z*
                      sin((avxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    (-((arz*Pi*rz*sin((arz*Pi*z)/L))/L) -
                     (arxyz*Pi*rxyz*x*y*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                    (w0 + wy*cos((awy*Pi*y)/L) +
                     wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                     wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)));

    val_source[1] = (apxyz*Pi*pxyz*y*z*cos((apxyz*Pi*x*y*z)/pow(L,3)))/
                    pow(L,3) - (apx*Pi*px*sin((apx*Pi*x)/L))/L +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    (u0 + uy*cos((auy*Pi*y)/L) +
                     uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                     ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L))*
                    ((awxyz*Pi*wxyz*x*y*
                      cos((awxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3) -
                     (awz*Pi*wz*sin((awz*Pi*z)/L))/L) +
                    (u0 + uy*cos((auy*Pi*y)/L) +
                     uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                     ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L))*
                    (v0 + vx*cos((avx*Pi*x)/L) +
                     vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                     vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L))*
                    (-((ary*Pi*ry*sin((ary*Pi*y)/L))/L) -
                     (arxyz*Pi*rxyz*x*z*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    pow(u0 + uy*cos((auy*Pi*y)/L) +
                          uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                          ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L),2)*
                    ((arx*Pi*rx*cos((arx*Pi*x)/L))/L -
                     (arxyz*Pi*rxyz*y*z*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    (v0 + vx*cos((avx*Pi*x)/L) +
                     vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                     vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L))*
                    (-((auy*Pi*uy*sin((auy*Pi*y)/L))/L) -
                     (auxyz*Pi*uxyz*x*z*
                      sin((auxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    2*(r0 + ry*cos((ary*Pi*y)/L) +
                       rz*cos((arz*Pi*z)/L) +
                       rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                       rx*sin((arx*Pi*x)/L))*
                    (u0 + uy*cos((auy*Pi*y)/L) +
                     uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                     ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L))*
                    ((aux*Pi*ux*cos((aux*Pi*x)/L))/L -
                     (auxyz*Pi*uxyz*y*z*
                      sin((auxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    (u0 + uy*cos((auy*Pi*y)/L) +
                     uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                     ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L))*
                    ((avy*Pi*vy*cos((avy*Pi*y)/L))/L -
                     (avxyz*Pi*vxyz*x*z*
                      sin((avxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    (u0 + uy*cos((auy*Pi*y)/L) +
                     uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                     ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L))*
                    (-((arz*Pi*rz*sin((arz*Pi*z)/L))/L) -
                     (arxyz*Pi*rxyz*x*y*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                    (w0 + wy*cos((awy*Pi*y)/L) +
                     wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                     wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3))) +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    ((auz*Pi*uz*cos((auz*Pi*z)/L))/L -
                     (auxyz*Pi*uxyz*x*y*
                      sin((auxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                    (w0 + wy*cos((awy*Pi*y)/L) +
                     wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                     wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)));

    val_source[2] = (apy*Pi*py*cos((apy*Pi*y)/L))/L +
                    (apxyz*Pi*pxyz*x*z*cos((apxyz*Pi*x*y*z)/pow(L,3)))/
                    pow(L,3) + (r0 + ry*cos((ary*Pi*y)/L) +
                                rz*cos((arz*Pi*z)/L) +
                                rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                                rx*sin((arx*Pi*x)/L))*
                               (v0 + vx*cos((avx*Pi*x)/L) +
                                vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                                vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L))*
                               ((awxyz*Pi*wxyz*x*y*
                                 cos((awxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3) -
                                (awz*Pi*wz*sin((awz*Pi*z)/L))/L) +
                    pow(v0 + vx*cos((avx*Pi*x)/L) +
                          vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                          vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L),2)*
                    (-((ary*Pi*ry*sin((ary*Pi*y)/L))/L) -
                     (arxyz*Pi*rxyz*x*z*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    (u0 + uy*cos((auy*Pi*y)/L) +
                     uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                     ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L))*
                    (v0 + vx*cos((avx*Pi*x)/L) +
                     vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                     vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L))*
                    ((arx*Pi*rx*cos((arx*Pi*x)/L))/L -
                     (arxyz*Pi*rxyz*y*z*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    (v0 + vx*cos((avx*Pi*x)/L) +
                     vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                     vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L))*
                    ((aux*Pi*ux*cos((aux*Pi*x)/L))/L -
                     (auxyz*Pi*uxyz*y*z*
                      sin((auxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    2*(r0 + ry*cos((ary*Pi*y)/L) +
                       rz*cos((arz*Pi*z)/L) +
                       rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                       rx*sin((arx*Pi*x)/L))*
                    (v0 + vx*cos((avx*Pi*x)/L) +
                     vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                     vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L))*
                    ((avy*Pi*vy*cos((avy*Pi*y)/L))/L -
                     (avxyz*Pi*vxyz*x*z*
                      sin((avxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    (u0 + uy*cos((auy*Pi*y)/L) +
                     uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                     ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L))*
                    (-((avx*Pi*vx*sin((avx*Pi*x)/L))/L) -
                     (avxyz*Pi*vxyz*y*z*
                      sin((avxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)) +
                    (v0 + vx*cos((avx*Pi*x)/L) +
                     vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                     vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L))*
                    (-((arz*Pi*rz*sin((arz*Pi*z)/L))/L) -
                     (arxyz*Pi*rxyz*x*y*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                    (w0 + wy*cos((awy*Pi*y)/L) +
                     wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                     wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3))) +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    ((avz*Pi*vz*cos((avz*Pi*z)/L))/L -
                     (avxyz*Pi*vxyz*x*y*
                      sin((avxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                    (w0 + wy*cos((awy*Pi*y)/L) +
                     wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                     wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)));

    val_source[3] = (apxyz*Pi*pxyz*x*y*cos((apxyz*Pi*x*y*z)/pow(L,3)))/
                    pow(L,3) - (apz*Pi*pz*sin((apz*Pi*z)/L))/L +
                    ((awx*Pi*wx*cos((awx*Pi*x)/L))/L +
                     (awxyz*Pi*wxyz*y*z*
                      cos((awxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    (u0 + uy*cos((auy*Pi*y)/L) +
                     uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                     ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L)) +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    ((awxyz*Pi*wxyz*x*z*
                      cos((awxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3) -
                     (awy*Pi*wy*sin((awy*Pi*y)/L))/L)*
                    (v0 + vx*cos((avx*Pi*x)/L) +
                     vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                     vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L)) +
                    2*(r0 + ry*cos((ary*Pi*y)/L) +
                       rz*cos((arz*Pi*z)/L) +
                       rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                       rx*sin((arx*Pi*x)/L))*
                    ((awxyz*Pi*wxyz*x*y*
                      cos((awxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3) -
                     (awz*Pi*wz*sin((awz*Pi*z)/L))/L)*
                    (w0 + wy*cos((awy*Pi*y)/L) +
                     wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                     wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3))) +
                    (v0 + vx*cos((avx*Pi*x)/L) +
                     vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                     vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L))*
                    (-((ary*Pi*ry*sin((ary*Pi*y)/L))/L) -
                     (arxyz*Pi*rxyz*x*z*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                    (w0 + wy*cos((awy*Pi*y)/L) +
                     wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                     wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3))) +
                    (u0 + uy*cos((auy*Pi*y)/L) +
                     uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                     ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L))*
                    ((arx*Pi*rx*cos((arx*Pi*x)/L))/L -
                     (arxyz*Pi*rxyz*y*z*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                    (w0 + wy*cos((awy*Pi*y)/L) +
                     wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                     wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3))) +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    ((aux*Pi*ux*cos((aux*Pi*x)/L))/L -
                     (auxyz*Pi*uxyz*y*z*
                      sin((auxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                    (w0 + wy*cos((awy*Pi*y)/L) +
                     wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                     wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3))) +
                    (r0 + ry*cos((ary*Pi*y)/L) +
                     rz*cos((arz*Pi*z)/L) +
                     rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                     rx*sin((arx*Pi*x)/L))*
                    ((avy*Pi*vy*cos((avy*Pi*y)/L))/L -
                     (avxyz*Pi*vxyz*x*z*
                      sin((avxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                    (w0 + wy*cos((awy*Pi*y)/L) +
                     wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                     wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3))) +
                    (-((arz*Pi*rz*sin((arz*Pi*z)/L))/L) -
                     (arxyz*Pi*rxyz*x*y*
                      sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                    pow(w0 + wy*cos((awy*Pi*y)/L) +
                          wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                          wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)),2);

    val_source[nDim+1] = ((awxyz*Pi*wxyz*x*y*cos((awxyz*Pi*x*y*z)/pow(L,3)))/
                          pow(L,3) - (awz*Pi*wz*sin((awz*Pi*z)/L))/L)*
                         ((g*(p0 + px*cos((apx*Pi*x)/L) +
                              pz*cos((apz*Pi*z)/L) +
                              py*sin((apy*Pi*y)/L) +
                              pxyz*sin((apxyz*Pi*x*y*z)/pow(L,3))))/
                          (-1 + g) + 0.5*
                                     (r0 + ry*cos((ary*Pi*y)/L) +
                                      rz*cos((arz*Pi*z)/L) +
                                      rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                                      rx*sin((arx*Pi*x)/L))*
                                     (pow(u0 + uy*cos((auy*Pi*y)/L) +
                                            uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                                            ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L)
                                             ,2) + pow(v0 + vx*cos((avx*Pi*x)/L) +
                                                         vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                                                         vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L)
                                             ,2) + pow(w0 + wy*cos((awy*Pi*y)/L) +
                                                         wz*cos((awz*Pi*z)/L) +
                                                         wx*sin((awx*Pi*x)/L) +
                                                         wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)),2))) +
                         ((aux*Pi*ux*cos((aux*Pi*x)/L))/L -
                          (auxyz*Pi*uxyz*y*z*
                           sin((auxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                         ((g*(p0 + px*cos((apx*Pi*x)/L) +
                              pz*cos((apz*Pi*z)/L) +
                              py*sin((apy*Pi*y)/L) +
                              pxyz*sin((apxyz*Pi*x*y*z)/pow(L,3))))/
                          (-1 + g) + 0.5*
                                     (r0 + ry*cos((ary*Pi*y)/L) +
                                      rz*cos((arz*Pi*z)/L) +
                                      rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                                      rx*sin((arx*Pi*x)/L))*
                                     (pow(u0 + uy*cos((auy*Pi*y)/L) +
                                            uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                                            ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L)
                                             ,2) + pow(v0 + vx*cos((avx*Pi*x)/L) +
                                                         vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                                                         vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L)
                                             ,2) + pow(w0 + wy*cos((awy*Pi*y)/L) +
                                                         wz*cos((awz*Pi*z)/L) +
                                                         wx*sin((awx*Pi*x)/L) +
                                                         wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)),2))) +
                         ((avy*Pi*vy*cos((avy*Pi*y)/L))/L -
                          (avxyz*Pi*vxyz*x*z*
                           sin((avxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                         ((g*(p0 + px*cos((apx*Pi*x)/L) +
                              pz*cos((apz*Pi*z)/L) +
                              py*sin((apy*Pi*y)/L) +
                              pxyz*sin((apxyz*Pi*x*y*z)/pow(L,3))))/
                          (-1 + g) + 0.5*
                                     (r0 + ry*cos((ary*Pi*y)/L) +
                                      rz*cos((arz*Pi*z)/L) +
                                      rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                                      rx*sin((arx*Pi*x)/L))*
                                     (pow(u0 + uy*cos((auy*Pi*y)/L) +
                                            uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                                            ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L)
                                             ,2) + pow(v0 + vx*cos((avx*Pi*x)/L) +
                                                         vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                                                         vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L)
                                             ,2) + pow(w0 + wy*cos((awy*Pi*y)/L) +
                                                         wz*cos((awz*Pi*z)/L) +
                                                         wx*sin((awx*Pi*x)/L) +
                                                         wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)),2))) +
                         (w0 + wy*cos((awy*Pi*y)/L) +
                          wz*cos((awz*Pi*z)/L) + wx*sin((awx*Pi*x)/L) +
                          wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)))*
                         ((g*((apxyz*Pi*pxyz*x*y*
                               cos((apxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)
                              - (apz*Pi*pz*sin((apz*Pi*z)/L))/L))/
                          (-1 + g) + 0.5*
                                     (r0 + ry*cos((ary*Pi*y)/L) +
                                      rz*cos((arz*Pi*z)/L) +
                                      rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                                      rx*sin((arx*Pi*x)/L))*
                                     (2*(u0 + uy*cos((auy*Pi*y)/L) +
                                         uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                                         ux*sin((aux*Pi*x)/L) +
                                         uz*sin((auz*Pi*z)/L))*
                                      ((auz*Pi*uz*cos((auz*Pi*z)/L))/L -
                                       (auxyz*Pi*uxyz*x*y*
                                        sin((auxyz*Pi*x*y*z)/pow(L,3)))/
                                       pow(L,3)) +
                                      2*(v0 + vx*cos((avx*Pi*x)/L) +
                                         vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                                         vy*sin((avy*Pi*y)/L) +
                                         vz*sin((avz*Pi*z)/L))*
                                      ((avz*Pi*vz*cos((avz*Pi*z)/L))/L -
                                       (avxyz*Pi*vxyz*x*y*
                                        sin((avxyz*Pi*x*y*z)/pow(L,3)))/
                                       pow(L,3)) +
                                      2*((awxyz*Pi*wxyz*x*y*
                                          cos((awxyz*Pi*x*y*z)/pow(L,3)))/
                                         pow(L,3) -
                                         (awz*Pi*wz*sin((awz*Pi*z)/L))/L)*
                                      (w0 + wy*cos((awy*Pi*y)/L) +
                                       wz*cos((awz*Pi*z)/L) +
                                       wx*sin((awx*Pi*x)/L) +
                                       wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)))) +
                          0.5*(-((arz*Pi*rz*sin((arz*Pi*z)/L))/L) -
                               (arxyz*Pi*rxyz*x*y*
                                sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                          (pow(u0 + uy*cos((auy*Pi*y)/L) +
                                 uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                                 ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L)
                                  ,2) + pow(v0 + vx*cos((avx*Pi*x)/L) +
                                              vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                                              vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L)
                                  ,2) + pow(w0 + wy*cos((awy*Pi*y)/L) +
                                              wz*cos((awz*Pi*z)/L) +
                                              wx*sin((awx*Pi*x)/L) +
                                              wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)),2))) +
                         (v0 + vx*cos((avx*Pi*x)/L) +
                          vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                          vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L))*
                         ((g*((apy*Pi*py*cos((apy*Pi*y)/L))/L +
                              (apxyz*Pi*pxyz*x*z*
                               cos((apxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)
                         ))/(-1 + g) +
                          0.5*(r0 + ry*cos((ary*Pi*y)/L) +
                               rz*cos((arz*Pi*z)/L) +
                               rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                               rx*sin((arx*Pi*x)/L))*
                          (2*(u0 + uy*cos((auy*Pi*y)/L) +
                              uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                              ux*sin((aux*Pi*x)/L) +
                              uz*sin((auz*Pi*z)/L))*
                           (-((auy*Pi*uy*sin((auy*Pi*y)/L))/L) -
                            (auxyz*Pi*uxyz*x*z*
                             sin((auxyz*Pi*x*y*z)/pow(L,3)))/
                            pow(L,3)) +
                           2*(v0 + vx*cos((avx*Pi*x)/L) +
                              vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                              vy*sin((avy*Pi*y)/L) +
                              vz*sin((avz*Pi*z)/L))*
                           ((avy*Pi*vy*cos((avy*Pi*y)/L))/L -
                            (avxyz*Pi*vxyz*x*z*
                             sin((avxyz*Pi*x*y*z)/pow(L,3)))/
                            pow(L,3)) +
                           2*((awxyz*Pi*wxyz*x*z*
                               cos((awxyz*Pi*x*y*z)/pow(L,3)))/
                              pow(L,3) -
                              (awy*Pi*wy*sin((awy*Pi*y)/L))/L)*
                           (w0 + wy*cos((awy*Pi*y)/L) +
                            wz*cos((awz*Pi*z)/L) +
                            wx*sin((awx*Pi*x)/L) +
                            wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)))) +
                          0.5*(-((ary*Pi*ry*sin((ary*Pi*y)/L))/L) -
                               (arxyz*Pi*rxyz*x*z*
                                sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                          (pow(u0 + uy*cos((auy*Pi*y)/L) +
                                 uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                                 ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L)
                                  ,2) + pow(v0 + vx*cos((avx*Pi*x)/L) +
                                              vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                                              vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L)
                                  ,2) + pow(w0 + wy*cos((awy*Pi*y)/L) +
                                              wz*cos((awz*Pi*z)/L) +
                                              wx*sin((awx*Pi*x)/L) +
                                              wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)),2))) +
                         (u0 + uy*cos((auy*Pi*y)/L) +
                          uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                          ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L))*
                         ((g*((apxyz*Pi*pxyz*y*z*
                               cos((apxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3)
                              - (apx*Pi*px*sin((apx*Pi*x)/L))/L))/
                          (-1 + g) + 0.5*
                                     (r0 + ry*cos((ary*Pi*y)/L) +
                                      rz*cos((arz*Pi*z)/L) +
                                      rxyz*cos((arxyz*Pi*x*y*z)/pow(L,3)) +
                                      rx*sin((arx*Pi*x)/L))*
                                     (2*(u0 + uy*cos((auy*Pi*y)/L) +
                                         uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                                         ux*sin((aux*Pi*x)/L) +
                                         uz*sin((auz*Pi*z)/L))*
                                      ((aux*Pi*ux*cos((aux*Pi*x)/L))/L -
                                       (auxyz*Pi*uxyz*y*z*
                                        sin((auxyz*Pi*x*y*z)/pow(L,3)))/
                                       pow(L,3)) +
                                      2*(v0 + vx*cos((avx*Pi*x)/L) +
                                         vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                                         vy*sin((avy*Pi*y)/L) +
                                         vz*sin((avz*Pi*z)/L))*
                                      (-((avx*Pi*vx*sin((avx*Pi*x)/L))/L) -
                                       (avxyz*Pi*vxyz*y*z*
                                        sin((avxyz*Pi*x*y*z)/pow(L,3)))/
                                       pow(L,3)) +
                                      2*((awx*Pi*wx*cos((awx*Pi*x)/L))/L +
                                         (awxyz*Pi*wxyz*y*z*
                                          cos((awxyz*Pi*x*y*z)/pow(L,3)))/
                                         pow(L,3))*
                                      (w0 + wy*cos((awy*Pi*y)/L) +
                                       wz*cos((awz*Pi*z)/L) +
                                       wx*sin((awx*Pi*x)/L) +
                                       wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)))) +
                          0.5*((arx*Pi*rx*cos((arx*Pi*x)/L))/L -
                               (arxyz*Pi*rxyz*y*z*
                                sin((arxyz*Pi*x*y*z)/pow(L,3)))/pow(L,3))*
                          (pow(u0 + uy*cos((auy*Pi*y)/L) +
                                 uxyz*cos((auxyz*Pi*x*y*z)/pow(L,3)) +
                                 ux*sin((aux*Pi*x)/L) + uz*sin((auz*Pi*z)/L)
                                  ,2) + pow(v0 + vx*cos((avx*Pi*x)/L) +
                                              vxyz*cos((avxyz*Pi*x*y*z)/pow(L,3)) +
                                              vy*sin((avy*Pi*y)/L) + vz*sin((avz*Pi*z)/L)
                                  ,2) + pow(w0 + wy*cos((awy*Pi*y)/L) +
                                              wz*cos((awz*Pi*z)/L) +
                                              wx*sin((awx*Pi*x)/L) +
                                              wxyz*sin((awxyz*Pi*x*y*z)/pow(L,3)),2)));


}

bool CUserDefinedSolution::IsManufacturedSolution() const {
    return true;
}
