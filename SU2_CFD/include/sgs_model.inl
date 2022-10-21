/*!
 * \file sgs_model.inl
 * \brief In-Line subroutines of the <i>sgs_model.hpp</i> file.
 * \author E. van der Weide, T. Economon, P. Urbanczyk
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

inline CSGSModel::CSGSModel(void){}
inline CSGSModel::~CSGSModel(void){}

inline su2double CSGSModel::ComputeEddyViscosity_2D(const su2double rho,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double lenScale,
                                                    const su2double distToWall) {
  return 0.0;
}

inline su2double CSGSModel::ComputeEddyViscosity_3D(const su2double rho,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dudz,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double dvdz,
                                                    const su2double dwdx,
                                                    const su2double dwdy,
                                                    const su2double dwdz,
                                                    const su2double lenScale,
                                                    const su2double distToWall) {
  return 0.0;
}

inline void CSGSModel::ComputeGradEddyViscosity_2D(const su2double rho,
                                                   const su2double drhodx,
                                                   const su2double drhody,
                                                   const su2double dudx,
                                                   const su2double dudy,
                                                   const su2double dvdx,
                                                   const su2double dvdy,
                                                   const su2double d2udx2,
                                                   const su2double d2udy2,
                                                   const su2double d2udxdy,
                                                   const su2double d2vdx2,
                                                   const su2double d2vdy2,
                                                   const su2double d2vdxdy,
                                                   const su2double lenScale,
                                                   const su2double distToWall,
                                                         su2double &dMuTdx,
                                                         su2double &dMuTdy) {
  dMuTdx = dMuTdy = 0.0;
}

inline void CSGSModel::ComputeGradEddyViscosity_3D(const su2double rho,
                                                   const su2double drhodx,
                                                   const su2double drhody,
                                                   const su2double drhodz,
                                                   const su2double dudx,
                                                   const su2double dudy,
                                                   const su2double dudz,
                                                   const su2double dvdx,
                                                   const su2double dvdy,
                                                   const su2double dvdz,
                                                   const su2double dwdx,
                                                   const su2double dwdy,
                                                   const su2double dwdz,
                                                   const su2double d2udx2,
                                                   const su2double d2udy2,
                                                   const su2double d2udz2,
                                                   const su2double d2udxdy,
                                                   const su2double d2udxdz,
                                                   const su2double d2udydz,
                                                   const su2double d2vdx2,
                                                   const su2double d2vdy2,
                                                   const su2double d2vdz2,
                                                   const su2double d2vdxdy,
                                                   const su2double d2vdxdz,
                                                   const su2double d2vdydz,
                                                   const su2double d2wdx2,
                                                   const su2double d2wdy2,
                                                   const su2double d2wdz2,
                                                   const su2double d2wdxdy,
                                                   const su2double d2wdxdz,
                                                   const su2double d2wdydz,
                                                   const su2double lenScale,
                                                   const su2double distToWall,
                                                         su2double &dMuTdx,
                                                         su2double &dMuTdy,
                                                         su2double &dMuTdz) {
  dMuTdx = dMuTdy = dMuTdz = 0.0;
}

inline CSmagorinskyModel::CSmagorinskyModel(void) : CSGSModel() {
  const_smag  = 0.1;
  filter_mult = 2.0;
}

inline CSmagorinskyModel::~CSmagorinskyModel(void){}

inline su2double CSmagorinskyModel::ComputeEddyViscosity_2D(const su2double rho,
                                                            const su2double dudx,
                                                            const su2double dudy,
                                                            const su2double dvdx,
                                                            const su2double dvdy,
                                                            const su2double lenScale,
                                                            const su2double distToWall) {
  /* Constant coefficient Smagorinsky SGS is calculated:
   * ( C_s * L_c )^2 * |S(x,t)|
   * C_s = Smagorinsky constant
   * L_c = Filter width
   * S(x,t) = Rate of Strain Tensor ( 1/2 [ du_i/dx_j + du_j/dx_i] )
   */
  const su2double C_s_filter_width = const_smag*filter_mult*lenScale;

  const su2double S12          = 0.5*(dudy + dvdx);
  const su2double strain_rate2 = 2.0*(dudx*dudx + dvdy*dvdy + 2.0*S12*S12);

  /* Return the SGS dynamic viscosity. */
  return rho*C_s_filter_width*C_s_filter_width*sqrt(strain_rate2);
}

inline su2double CSmagorinskyModel::ComputeEddyViscosity_3D(const su2double rho,
                                                            const su2double dudx,
                                                            const su2double dudy,
                                                            const su2double dudz,
                                                            const su2double dvdx,
                                                            const su2double dvdy,
                                                            const su2double dvdz,
                                                            const su2double dwdx,
                                                            const su2double dwdy,
                                                            const su2double dwdz,
                                                            const su2double lenScale,
                                                            const su2double distToWall) {
  /* Constant coefficient Smagorinsky SGS is calculated:
   * ( C_s * L_c )^2 * |S(x,t)|
   * C_s = Smagorinsky constant
   * L_c = Filter width
   * S(x,t) = Rate of Strain Tensor ( 1/2 [ du_i/dx_j + du_j/dx_i] )
   */
  const su2double C_s_filter_width = const_smag*filter_mult*lenScale;

  const su2double S12 = 0.5*(dudy + dvdx);
  const su2double S13 = 0.5*(dudz + dwdx);
  const su2double S23 = 0.5*(dvdz + dwdy);

  const su2double strain_rate2 = 2.0*(dudx*dudx + dvdy*dvdy + dwdz*dwdz
                               +      2.0*(S12*S12 + S13*S13 + S23*S23));

  /* Return the SGS dynamic viscosity. */
  return rho*C_s_filter_width*C_s_filter_width*sqrt(strain_rate2);
}

inline void CSmagorinskyModel::ComputeGradEddyViscosity_2D(const su2double rho,
                                                           const su2double drhodx,
                                                           const su2double drhody,
                                                           const su2double dudx,
                                                           const su2double dudy,
                                                           const su2double dvdx,
                                                           const su2double dvdy,
                                                           const su2double d2udx2,
                                                           const su2double d2udy2,
                                                           const su2double d2udxdy,
                                                           const su2double d2vdx2,
                                                           const su2double d2vdy2,
                                                           const su2double d2vdxdy,
                                                           const su2double lenScale,
                                                           const su2double distToWall,
                                                                 su2double &dMuTdx,
                                                                 su2double &dMuTdy) {
  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

inline void CSmagorinskyModel::ComputeGradEddyViscosity_3D(const su2double rho,
                                                           const su2double drhodx,
                                                           const su2double drhody,
                                                           const su2double drhodz,
                                                           const su2double dudx,
                                                           const su2double dudy,
                                                           const su2double dudz,
                                                           const su2double dvdx,
                                                           const su2double dvdy,
                                                           const su2double dvdz,
                                                           const su2double dwdx,
                                                           const su2double dwdy,
                                                           const su2double dwdz,
                                                           const su2double d2udx2,
                                                           const su2double d2udy2,
                                                           const su2double d2udz2,
                                                           const su2double d2udxdy,
                                                           const su2double d2udxdz,
                                                           const su2double d2udydz,
                                                           const su2double d2vdx2,
                                                           const su2double d2vdy2,
                                                           const su2double d2vdz2,
                                                           const su2double d2vdxdy,
                                                           const su2double d2vdxdz,
                                                           const su2double d2vdydz,
                                                           const su2double d2wdx2,
                                                           const su2double d2wdy2,
                                                           const su2double d2wdz2,
                                                           const su2double d2wdxdy,
                                                           const su2double d2wdxdz,
                                                           const su2double d2wdydz,
                                                           const su2double lenScale,
                                                           const su2double distToWall,
                                                                 su2double &dMuTdx,
                                                                 su2double &dMuTdy,
                                                                 su2double &dMuTdz) {
  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

inline CWALEModel::CWALEModel(void) : CSGSModel() {
  const_WALE = 0.325;
}

inline CWALEModel::~CWALEModel(void){}

inline su2double CWALEModel::ComputeEddyViscosity_2D(const su2double rho,
                                                     const su2double dudx,
                                                     const su2double dudy,
                                                     const su2double dvdx,
                                                     const su2double dvdy,
                                                     const su2double lenScale,
                                                     const su2double distToWall) {
  /* Compute the length scale in WALE. */
  const su2double lenScaleWale = const_WALE*lenScale;

  /* Compute the strain rate tensor, which is symmetric. */
  const su2double S11 = dudx, S22 = dvdy;
  const su2double S12 = 0.5*(dudy + dvdx);

  /* Compute the values of the Sd tensor. First without the trace
     correction of the diagonal terms. */
  su2double Sd11 = dudx*dudx + dudy*dvdx;
  su2double Sd22 = dvdx*dudy + dvdy*dvdy;

  const su2double Sd12 = 0.5*(dudx*dudy + dudy*dvdy + dvdx*dudx + dvdy*dvdx);

  /* Correct the diagonal elements, such that the trace of the Sd tensor is zero
     Note that this comes from the 3D formulation. */
  const su2double thirdTrace = (Sd11 + Sd22)/3.0;

  Sd11 -= thirdTrace;
  Sd22 -= thirdTrace;

  /* Compute the summation of both tensors. */
  const su2double sumS  = S11 *S11  + S22 *S22  + 2.0*S12 *S12;
  const su2double sumSd = Sd11*Sd11 + Sd22*Sd22 + 2.0*Sd12*Sd12;

  /* Compute the kinematic eddy viscosity. */
  const su2double sumSdPow3_2 = sumSd*sqrt(sumSd);
  const su2double sumSdPow5_4 = sqrt(sumSdPow3_2*sumSd);
  const su2double sumSPow5_2  = sumS*sumS*sqrt(sumS);
  const su2double denom       = sumSPow5_2 + sumSdPow5_4;

  const su2double nuEddy = lenScaleWale*lenScaleWale*sumSdPow3_2
                         / max(denom, 1.e-20);

  /* Return the SGS dynamic viscosity. */
  return rho*nuEddy;
}

inline su2double CWALEModel::ComputeEddyViscosity_3D(const su2double rho,
                                                     const su2double dudx,
                                                     const su2double dudy,
                                                     const su2double dudz,
                                                     const su2double dvdx,
                                                     const su2double dvdy,
                                                     const su2double dvdz,
                                                     const su2double dwdx,
                                                     const su2double dwdy,
                                                     const su2double dwdz,
                                                     const su2double lenScale,
                                                     const su2double distToWall) {
  /* Compute the length scale in WALE. */
  const su2double lenScaleWale = const_WALE*lenScale;

  /* Compute the strain rate tensor, which is symmetric. */
  const su2double S11 = dudx, S22 = dvdy, S33 = dwdz;
  const su2double S12 = 0.5*(dudy + dvdx);
  const su2double S13 = 0.5*(dudz + dwdx);
  const su2double S23 = 0.5*(dvdz + dwdy);

  /* Compute the values of the Sd tensor. First without the trace
     correction of the diagonal terms. */
  su2double Sd11 = dudx*dudx + dudy*dvdx + dudz*dwdx;
  su2double Sd22 = dvdx*dudy + dvdy*dvdy + dvdz*dwdy;
  su2double Sd33 = dwdx*dudz + dwdy*dvdz + dwdz*dwdz;

  const su2double Sd12 = 0.5*(dudx*dudy + dudy*dvdy + dudz*dwdy
                       +      dvdx*dudx + dvdy*dvdx + dvdz*dwdx);
  const su2double Sd13 = 0.5*(dudx*dudz + dudy*dvdz + dudz*dwdz
                       +      dwdx*dudx + dwdy*dvdx + dwdz*dwdx);
  const su2double Sd23 = 0.5*(dvdx*dudz + dvdy*dvdz + dvdz*dwdz
                       +      dwdx*dudy + dwdy*dvdy + dwdz*dwdy);

  /* Correct the diagonal elements, such that the trace of the Sd tensor is zero. */
  const su2double thirdTrace = (Sd11 + Sd22 + Sd33)/3.0;

  Sd11 -= thirdTrace;
  Sd22 -= thirdTrace;
  Sd33 -= thirdTrace;

  /* Compute the summation of both tensors. */
  const su2double sumS  = S11*S11 + S22*S22 + S33*S33
                        + 2.0*(S12*S12 + S13*S13 + S23*S23);
  const su2double sumSd = Sd11*Sd11 + Sd22*Sd22 + Sd33*Sd33
                        + 2.0*(Sd12*Sd12 + Sd13*Sd13 + Sd23*Sd23);

  /* Compute the kinematic eddy viscosity. */
  const su2double sumSdPow3_2 = sumSd*sqrt(sumSd);
  const su2double sumSdPow5_4 = sqrt(sumSdPow3_2*sumSd);
  const su2double sumSPow5_2  = sumS*sumS*sqrt(sumS);
  const su2double denom       = sumSPow5_2 + sumSdPow5_4;

  const su2double nuEddy = lenScaleWale*lenScaleWale*sumSdPow3_2
                         / max(denom, 1.e-20);

  /* Return the SGS dynamic viscosity. */
  return rho*nuEddy;
}

inline void CWALEModel::ComputeGradEddyViscosity_2D(const su2double rho,
                                                    const su2double drhodx,
                                                    const su2double drhody,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double d2udx2,
                                                    const su2double d2udy2,
                                                    const su2double d2udxdy,
                                                    const su2double d2vdx2,
                                                    const su2double d2vdy2,
                                                    const su2double d2vdxdy,
                                                    const su2double lenScale,
                                                    const su2double distToWall,
                                                          su2double &dMuTdx,
                                                          su2double &dMuTdy) {
  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

inline void CWALEModel::ComputeGradEddyViscosity_3D(const su2double rho,
                                                    const su2double drhodx,
                                                    const su2double drhody,
                                                    const su2double drhodz,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dudz,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double dvdz,
                                                    const su2double dwdx,
                                                    const su2double dwdy,
                                                    const su2double dwdz,
                                                    const su2double d2udx2,
                                                    const su2double d2udy2,
                                                    const su2double d2udz2,
                                                    const su2double d2udxdy,
                                                    const su2double d2udxdz,
                                                    const su2double d2udydz,
                                                    const su2double d2vdx2,
                                                    const su2double d2vdy2,
                                                    const su2double d2vdz2,
                                                    const su2double d2vdxdy,
                                                    const su2double d2vdxdz,
                                                    const su2double d2vdydz,
                                                    const su2double d2wdx2,
                                                    const su2double d2wdy2,
                                                    const su2double d2wdz2,
                                                    const su2double d2wdxdy,
                                                    const su2double d2wdxdz,
                                                    const su2double d2wdydz,
                                                    const su2double lenScale,
                                                    const su2double distToWall,
                                                          su2double &dMuTdx,
                                                          su2double &dMuTdy,
                                                          su2double &dMuTdz) {
  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

inline CVremanModel::CVremanModel(void) : CSGSModel() {

  /* const_Vreman = 2.5*Cs*Cs where Cs is the Smagorinsky constant */
  const_Vreman = 0.07;
}

inline CVremanModel::~CVremanModel(void){}

inline su2double CVremanModel::ComputeEddyViscosity_2D(const su2double rho,
                                                     const su2double dudx,
                                                     const su2double dudy,
                                                     const su2double dvdx,
                                                     const su2double dvdy,
                                                     const su2double lenScale,
                                                     const su2double distToWall) {
  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
  return 0;
}

inline su2double CVremanModel::ComputeEddyViscosity_3D(const su2double rho,
                                                     const su2double dudx,
                                                     const su2double dudy,
                                                     const su2double dudz,
                                                     const su2double dvdx,
                                                     const su2double dvdy,
                                                     const su2double dvdz,
                                                     const su2double dwdx,
                                                     const su2double dwdy,
                                                     const su2double dwdz,
                                                     const su2double lenScale,
                                                     const su2double distToWall) {

  su2double alpha11 = dudx;
  su2double alpha22 = dvdy;
  su2double alpha33 = dwdz;

  //Check if it is necessary to remove the trace.
  const su2double tmp = (alpha11 + alpha22 + alpha33)/3.0;
  alpha11 -= tmp;
  alpha22 -= tmp;
  alpha33 -= tmp;

  const su2double lenScale2 = lenScale * lenScale;
  const su2double alpha12 = dudy;
  const su2double alpha13 = dudz;
  const su2double alpha23 = dvdz;

  const su2double alpha21 = dvdx;
  const su2double alpha31 = dwdx;
  const su2double alpha32 = dwdy;

  const su2double beta11  = lenScale2*alpha11*alpha11 + lenScale2*alpha12*alpha12 + lenScale2*alpha13*alpha13 ;
  const su2double beta12  = lenScale2*alpha11*alpha21 + lenScale2*alpha12*alpha22 + lenScale2*alpha13*alpha23 ;
  const su2double beta13  = lenScale2*alpha11*alpha31 + lenScale2*alpha12*alpha32 + lenScale2*alpha13*alpha33 ;
  const su2double beta22  = lenScale2*alpha21*alpha21 + lenScale2*alpha22*alpha22 + lenScale2*alpha23*alpha23 ;
  const su2double beta23  = lenScale2*alpha21*alpha31 + lenScale2*alpha22*alpha32 + lenScale2*alpha23*alpha33 ;
  const su2double beta33  = lenScale2*alpha31*alpha31 + lenScale2*alpha32*alpha32 + lenScale2*alpha33*alpha33 ;

  su2double B = beta11*beta22-beta12*beta12+beta11*beta33-beta13*beta13+beta22*beta33-beta23*beta23;
            B = (B + fabs(B))*0.5;
  const su2double denon    = alpha11*alpha11+alpha22*alpha22+alpha33*alpha33 +
                            alpha12*alpha12+alpha13*alpha13+alpha23*alpha23 +
                            alpha21*alpha21+alpha31*alpha31+alpha32*alpha32;

  const su2double nuEddy_Vreman = sqrt(B/(denon+1.0E-20));

 /* Return the SGS dynamic viscosity. */
 return rho*const_Vreman*nuEddy_Vreman;

}

inline void CVremanModel::ComputeGradEddyViscosity_2D(const su2double rho,
                                                    const su2double drhodx,
                                                    const su2double drhody,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double d2udx2,
                                                    const su2double d2udy2,
                                                    const su2double d2udxdy,
                                                    const su2double d2vdx2,
                                                    const su2double d2vdy2,
                                                    const su2double d2vdxdy,
                                                    const su2double lenScale,
                                                    const su2double distToWall,
                                                    su2double &dMuTdx,
                                                    su2double &dMuTdy) {
  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

inline void CVremanModel::ComputeGradEddyViscosity_3D(const su2double rho,
                                                    const su2double drhodx,
                                                    const su2double drhody,
                                                    const su2double drhodz,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dudz,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double dvdz,
                                                    const su2double dwdx,
                                                    const su2double dwdy,
                                                    const su2double dwdz,
                                                    const su2double d2udx2,
                                                    const su2double d2udy2,
                                                    const su2double d2udz2,
                                                    const su2double d2udxdy,
                                                    const su2double d2udxdz,
                                                    const su2double d2udydz,
                                                    const su2double d2vdx2,
                                                    const su2double d2vdy2,
                                                    const su2double d2vdz2,
                                                    const su2double d2vdxdy,
                                                    const su2double d2vdxdz,
                                                    const su2double d2vdydz,
                                                    const su2double d2wdx2,
                                                    const su2double d2wdy2,
                                                    const su2double d2wdz2,
                                                    const su2double d2wdxdy,
                                                    const su2double d2wdxdz,
                                                    const su2double d2wdydz,
                                                    const su2double lenScale,
                                                    const su2double distToWall,
                                                    su2double &dMuTdx,
                                                    su2double &dMuTdy,
                                                    su2double &dMuTdz) {
  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

inline CSIGMAModel::CSIGMAModel(void) : CSGSModel() {
  const_SIGMA = 1.35;
}

inline CSIGMAModel::~CSIGMAModel(void){}

inline su2double CSIGMAModel::ComputeEddyViscosity_2D(const su2double rho,
                                                     const su2double dudx,
                                                     const su2double dudy,
                                                     const su2double dvdx,
                                                     const su2double dvdy,
                                                     const su2double lenScale,
                                                     const su2double distToWall) {
  cout << "CSIGMAModel::ComputeEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);

}

inline su2double CSIGMAModel::ComputeEddyViscosity_3D(const su2double rho,
                                                     const su2double dudx,
                                                     const su2double dudy,
                                                     const su2double dudz,
                                                     const su2double dvdx,
                                                     const su2double dvdy,
                                                     const su2double dvdz,
                                                     const su2double dwdx,
                                                     const su2double dwdy,
                                                     const su2double dwdz,
                                                     const su2double lenScale,
                                                     const su2double distToWall) {
  /* Compute the length scale in SIGMA. */
  const su2double lenScaleSigma = const_SIGMA*lenScale;

  /* Velocity gradients (gij) */
  const su2double gij[3][3] = {{dudx, dudy, dudz},
                               {dvdx, dvdy, dvdz},
                               {dwdx, dwdy, dwdz}};

  /* Transpose of gji  */
  su2double gji[3][3];
  for (auto iDim = 0; iDim < 3; iDim++)
    for (auto jDim = 0; jDim < 3; jDim++)
      gji[jDim][iDim] = gij[iDim][jDim];

  /* Matrix Multiplication (gji,gij) */
  su2double G[3][3];
  for (auto iDim = 0; iDim < 3; iDim++){
    for (auto jDim = 0; jDim < 3; jDim++){
      G[iDim][jDim] = 0.;
      for (auto kDim = 0; kDim < 3; kDim++){
        G[iDim][jDim] += gji[iDim][kDim] * gij[kDim][jDim];
      }
    }
  }

  /* Matrix Multiplication (G,G) */
  su2double G2[3][3];
  for (auto iDim = 0; iDim < 3; iDim++){
    for (auto jDim = 0; jDim < 3; jDim++){
      G2[iDim][jDim] = 0.;
      for (auto kDim = 0; kDim < 3; kDim++){
        G2[iDim][jDim] += G[iDim][kDim] * G[kDim][jDim];
      }
    }
  }

  /* Compute Invariants:
     I1 = tr(G)
     I2 = 0.5 (tr(G)^2 - tr(G^2))
     I3 = det(G)*/
  su2double const invariant1 = G[0][0] + G[1][1] + G[2][2];
  su2double const invariant2 = 0.5*( pow(invariant1,2) - ( G2[0][0] + G2[1][1] + G2[2][2] ) );
  su2double const invariant3 = G[0][0]*(G[1][1]*G[2][2] - G[1][2]*G[2][1]) - G[0][1]*(G[1][0]*G[2][2] - G[1][2]*G[2][0]) + G[0][2]*(G[1][0]*G[2][1] - G[1][1]*G[2][0]);

  /* Compute angles alphai */
  su2double alpha1 = (pow(invariant1,2.) / 9.)  - (invariant2 / 3.);
  su2double alpha2 = (pow(invariant1,3.) / 27.) - (invariant1*invariant2/6.) + (0.5*invariant3);

  if ( alpha1 < 0.0 ) alpha1 = abs(alpha1); // avoid sqrt of negative value

  /* Compute alpha_3 */
  su2double alpha3 = 0.0;
  su2double pre_alpha3 = alpha2 / pow(alpha1,3./2);
  if ( pre_alpha3>=1. || std::isnan(pre_alpha3) ) alpha3 = 0.0;
  else if ( pre_alpha3<=-1. ) alpha3 = PI_NUMBER / 3.;
  else alpha3 = acos( pre_alpha3 )/3.;

  /* Compute singular values sigma_i */
  su2double sigma1 = sqrt( abs( invariant1/3. + 2.*sqrt(alpha1)*cos(alpha3) ) );
  su2double sigma2 = sqrt( abs( invariant1/3. - 2.*sqrt(alpha1)*cos(PI_NUMBER/3.+alpha3) ) );
  su2double sigma3 = sqrt( abs( invariant1/3. - 2.*sqrt(alpha1)*cos(PI_NUMBER/3.-alpha3) ) );

  /* Compute D_sigma */
  su2double D_sigma = sigma3*( sigma1 - sigma2 )*( sigma2 - sigma3 );

  const su2double nuEddy = lenScaleSigma*lenScaleSigma*D_sigma;

  /* Return the SGS dynamic viscosity. */
  return rho*nuEddy;
}

inline void CSIGMAModel::ComputeGradEddyViscosity_2D(const su2double rho,
                                                    const su2double drhodx,
                                                    const su2double drhody,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double d2udx2,
                                                    const su2double d2udy2,
                                                    const su2double d2udxdy,
                                                    const su2double d2vdx2,
                                                    const su2double d2vdy2,
                                                    const su2double d2vdxdy,
                                                    const su2double lenScale,
                                                    const su2double distToWall,
                                                          su2double &dMuTdx,
                                                          su2double &dMuTdy) {
  cout << "CSIGMAModel::ComputeGradEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
}

inline void CSIGMAModel::ComputeGradEddyViscosity_3D(const su2double rho,
                                                    const su2double drhodx,
                                                    const su2double drhody,
                                                    const su2double drhodz,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dudz,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double dvdz,
                                                    const su2double dwdx,
                                                    const su2double dwdy,
                                                    const su2double dwdz,
                                                    const su2double d2udx2,
                                                    const su2double d2udy2,
                                                    const su2double d2udz2,
                                                    const su2double d2udxdy,
                                                    const su2double d2udxdz,
                                                    const su2double d2udydz,
                                                    const su2double d2vdx2,
                                                    const su2double d2vdy2,
                                                    const su2double d2vdz2,
                                                    const su2double d2vdxdy,
                                                    const su2double d2vdxdz,
                                                    const su2double d2vdydz,
                                                    const su2double d2wdx2,
                                                    const su2double d2wdy2,
                                                    const su2double d2wdz2,
                                                    const su2double d2wdxdy,
                                                    const su2double d2wdxdz,
                                                    const su2double d2wdydz,
                                                    const su2double lenScale,
                                                    const su2double distToWall,
                                                          su2double &dMuTdx,
                                                          su2double &dMuTdy,
                                                          su2double &dMuTdz) {
  cout << "CSIGMAModel::ComputeGradEddyViscosity_3D: Not implemented yet" << endl;
  exit(1);
}

inline CAMDModel::CAMDModel(void) : CSGSModel() {

  const_AMD = 0.3;
}

inline CAMDModel::~CAMDModel(void){}

inline su2double CAMDModel::ComputeEddyViscosity_2D(const su2double rho,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double lenScale,
                                                    const su2double distToWall) {
  cout << "CAMDModel::ComputeEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
}

inline su2double CAMDModel::ComputeEddyViscosity_3D(const su2double rho,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dudz,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double dvdz,
                                                    const su2double dwdx,
                                                    const su2double dwdy,
                                                    const su2double dwdz,
                                                    const su2double lenScale,
                                                    const su2double distToWall) {

  su2double S11 = dudx;
  su2double S12 = 0.5*(dudy+dvdx);
  su2double S13 = 0.5*(dudz+dwdx);
  su2double S22 = dvdy;
  su2double S23 = 0.5*(dvdz+dwdy);
  su2double S33 = dwdy;

  const su2double beta11 =     (dudx*dudx + dudy*dudy + dudz*dudz)*S11;
  const su2double beta12 = 2.0*(dudx*dvdx + dudy*dvdy + dudz*dvdz)*S12;
  const su2double beta13 = 2.0*(dudx*dwdx + dudy*dwdy + dudz*dwdz)*S13;
  const su2double beta22 =     (dvdx*dvdx + dvdy*dvdy + dvdz*dvdz)*S22;
  const su2double beta23 = 2.0*(dvdx*dwdx + dvdy*dwdy + dvdz*dwdz)*S23;
  const su2double beta33 =     (dwdx*dwdx + dwdy*dwdy + dwdz*dwdz)*S33;

  const su2double alpha1 = dudx*dudx + dudy*dudy + dudz*dudz;
  const su2double alpha2 = dvdx*dvdx + dvdy*dvdy + dvdz*dvdz;
  const su2double alpha3 = dwdx*dwdx + dwdy*dwdy + dwdz*dwdz;

  const su2double lenScale2 = lenScale * lenScale;
  su2double B = -lenScale2*(beta11 + beta12 + beta13 + beta22 + beta23 + beta33);
      B = max(B,0.0);

  const su2double denon = alpha1 + alpha2 + alpha3;

  const su2double nuEddy_AMD = B/denon;

 /* Return the SGS dynamic viscosity. */
 return rho*const_AMD*nuEddy_AMD;

}

inline void CAMDModel::ComputeGradEddyViscosity_2D(const su2double rho,
                                                   const su2double drhodx,
                                                   const su2double drhody,
                                                   const su2double dudx,
                                                   const su2double dudy,
                                                   const su2double dvdx,
                                                   const su2double dvdy,
                                                   const su2double d2udx2,
                                                   const su2double d2udy2,
                                                   const su2double d2udxdy,
                                                   const su2double d2vdx2,
                                                   const su2double d2vdy2,
                                                   const su2double d2vdxdy,
                                                   const su2double lenScale,
                                                   const su2double distToWall,
                                                   su2double &dMuTdx,
                                                   su2double &dMuTdy) {
  cout << "CAMDModel::ComputeGradEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
}

inline void CAMDModel::ComputeGradEddyViscosity_3D(const su2double rho,
                                                   const su2double drhodx,
                                                   const su2double drhody,
                                                   const su2double drhodz,
                                                   const su2double dudx,
                                                   const su2double dudy,
                                                   const su2double dudz,
                                                   const su2double dvdx,
                                                   const su2double dvdy,
                                                   const su2double dvdz,
                                                   const su2double dwdx,
                                                   const su2double dwdy,
                                                   const su2double dwdz,
                                                   const su2double d2udx2,
                                                   const su2double d2udy2,
                                                   const su2double d2udz2,
                                                   const su2double d2udxdy,
                                                   const su2double d2udxdz,
                                                   const su2double d2udydz,
                                                   const su2double d2vdx2,
                                                   const su2double d2vdy2,
                                                   const su2double d2vdz2,
                                                   const su2double d2vdxdy,
                                                   const su2double d2vdxdz,
                                                   const su2double d2vdydz,
                                                   const su2double d2wdx2,
                                                   const su2double d2wdy2,
                                                   const su2double d2wdz2,
                                                   const su2double d2wdxdy,
                                                   const su2double d2wdxdz,
                                                   const su2double d2wdydz,
                                                   const su2double lenScale,
                                                   const su2double distToWall,
                                                   su2double &dMuTdx,
                                                   su2double &dMuTdy,
                                                   su2double &dMuTdz) {
  cout << "CAMDModel::ComputeGradEddyViscosity_3D: Not implemented yet" << endl;
  exit(1);
}
