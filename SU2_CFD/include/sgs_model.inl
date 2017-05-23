/*!
 * \file sgs_model.inl
 * \brief In-Line subroutines of the <i>sgs_model.hpp</i> file.
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
                                                         su2double &dMuTdy) {}

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
                                                         su2double &dMuTdz) {}

inline CSmagorinskyModel::CSmagorinskyModel(void) : CSGSModel() {}
inline CSmagorinskyModel::~CSmagorinskyModel(void){}

inline su2double CSmagorinskyModel::ComputeEddyViscosity_2D(const su2double rho,
                                                            const su2double dudx,
                                                            const su2double dudy,
                                                            const su2double dvdx,
                                                            const su2double dvdy,
                                                            const su2double lenScale,
                                                            const su2double distToWall) {
  cout << "CSmagorinskyModel::ComputeEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
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
  cout << "CSmagorinskyModel::ComputeEddyViscosity_3D: Not implemented yet" << endl;
  exit(1);
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
  cout << "CSmagorinskyModel::ComputeGradEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
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
  cout << "CSmagorinskyModel::ComputeGradEddyViscosity_3D: Not implemented yet" << endl;
  exit(1);
}

inline CWALEModel::CWALEModel(void) : CSGSModel() {}
inline CWALEModel::~CWALEModel(void){}

inline su2double CWALEModel::ComputeEddyViscosity_2D(const su2double rho,
                                                     const su2double dudx,
                                                     const su2double dudy,
                                                     const su2double dvdx,
                                                     const su2double dvdy,
                                                     const su2double lenScale,
                                                     const su2double distToWall) {
  cout << "CWALEModel::ComputeEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
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
  cout << "CWALEModel::ComputeEddyViscosity_3D: Not implemented yet" << endl;
  exit(1);
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
  cout << "CWALEModel::ComputeGradEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
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
  cout << "CWALEModel::ComputeGradEddyViscosity_3D: Not implemented yet" << endl;
  exit(1);
}
