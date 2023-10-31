/*!
 * \file sgs_model.hpp
 * \brief Headers of the LES subgrid scale models of the SU2 solvers.
 * \author E. van der Weide, T. Economon, P. Urbanczyk
 * \version 8.0.0 "Harrier"
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

#pragma once

#include "../../Common/include/parallelization/mpi_structure.hpp"

#include <iostream>
#include <cmath>

using namespace std;

/*!
 * \class CSGSModel
 * \brief Base class for defining the LES subgrid scale model.
 * \author: E. van der Weide, T. Economon, P. Urbanczyk
 * \version 8.0.0 "Harrier"
 */
class CSGSModel {

public:
  /*!
   * \brief Constructor of the class.
   */
  CSGSModel(void);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CSGSModel(void);

  /*!
   * \brief Virtual function to determine the eddy viscosity
            for the given function arguments for a 2D simulation.
   * \param[in] rho        - Density
   * \param[in] dudx       - x-derivative of the u-velocity.
   * \param[in] dudy       - y-derivative of the u-velocity.
   * \param[in] dvdx       - x-derivative of the v-velocity.
   * \param[in] dvdy       - y-derivative of the v-velocity.
   * \param[in] lenScale   - Length scale of the corresponding element.
   * \param[in] distToWall - Distance to the nearest wall.
   * \return Value of the dynamic eddy viscosity. For the base class 0 is returned.
   */
  virtual su2double ComputeEddyViscosity_2D(const su2double rho,
                                            const su2double dudx,
                                            const su2double dudy,
                                            const su2double dvdx,
                                            const su2double dvdy,
                                            const su2double lenScale,
                                            const su2double distToWall);

  /*!
   * \brief Virtual function to determine the eddy viscosity
            for the given function arguments for a 3D simulation.
   * \param[in] rho        - Density
   * \param[in] dudx       - x-derivative of the u-velocity.
   * \param[in] dudy       - y-derivative of the u-velocity.
   * \param[in] dudz       - z-derivative of the u-velocity.
   * \param[in] dvdx       - x-derivative of the v-velocity.
   * \param[in] dvdy       - y-derivative of the v-velocity.
   * \param[in] dvdz       - z-derivative of the v-velocity.
   * \param[in] dwdx       - x-derivative of the w-velocity.
   * \param[in] dwdy       - y-derivative of the w-velocity.
   * \param[in] dwdz       - z-derivative of the w-velocity.
   * \param[in] lenScale   - Length scale of the corresponding element.
   * \param[in] distToWall - Distance to the nearest wall.
   * \return Value of the dynamic eddy viscosity. For the base class 0 is returned.
   */
  virtual su2double ComputeEddyViscosity_3D(const su2double rho,
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
                                            const su2double distToWall);

  /*!
   * \brief Virtual function to determine the gradients of the eddy viscosity
            for the given function arguments for a 2D simulation.
   * \param[in]  rho        - Density.
   * \param[in]  drhodx     - x-derivative of the density.
   * \param[in]  drhody     - y-derivative of the density.
   * \param[in]  dudx       - x-derivative of the u-velocity.
   * \param[in]  dudy       - y-derivative of the u-velocity.
   * \param[in]  dvdx       - x-derivative of the v-velocity.
   * \param[in]  dvdy       - y-derivative of the v-velocity.
   * \param[in]  d2udx2     - 2nd x-derivative of the u-velocity.
   * \param[in]  d2udy2     - 2nd y-derivative of the u-velocity.
   * \param[in]  d2udxdy    - x-y cross-derivative of the u-velocity.
   * \param[in]  d2vdx2     - 2nd x-derivative of the v-velocity.
   * \param[in]  d2vdy2     - 2nd y-derivative of the v-velocity.
   * \param[in]  d2vdxdy    - x-y cross-derivative of the v-velocity.
   * \param[in]  lenScale   - Length scale of the corresponding element.
   * \param[in]  distToWall - Distance to the nearest wall.
   * \param[out] dMuTdx     - x-derivative of the turbulent viscosity.
   * \param[out] dMuTdy     - y-derivative of the turbulent viscosity.
   */
  virtual void ComputeGradEddyViscosity_2D(const su2double rho,
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
                                                 su2double &dMuTdy);

  /*!
   * \brief Virtual function to determine the gradients of the eddy viscosity
            for the given function arguments for a 3D simulation.
   * \param[in]  rho        - Density.
   * \param[in]  drhodx     - x-derivative of the density.
   * \param[in]  drhody     - y-derivative of the density.
   * \param[in]  drhodz     - z-derivative of the density.
   * \param[in]  dudx       - x-derivative of the u-velocity.
   * \param[in]  dudy       - y-derivative of the u-velocity.
   * \param[in]  dudz       - z-derivative of the u-velocity.
   * \param[in]  dvdx       - x-derivative of the v-velocity.
   * \param[in]  dvdy       - y-derivative of the v-velocity.
   * \param[in]  dvdz       - z-derivative of the v-velocity.
   * \param[in]  dwdx       - x-derivative of the w-velocity.
   * \param[in]  dwdy       - y-derivative of the w-velocity.
   * \param[in]  dwdz       - z-derivative of the w-velocity.
   * \param[in]  d2udx2     - 2nd x-derivative of the u-velocity.
   * \param[in]  d2udy2     - 2nd y-derivative of the u-velocity.
   * \param[in]  d2udz2     - 2nd z-derivative of the u-velocity.
   * \param[in]  d2udxdy    - x-y cross-derivative of the u-velocity.
   * \param[in]  d2udxdz    - x-z cross-derivative of the u-velocity.
   * \param[in]  d2udydz    - y-z cross-derivative of the u-velocity.
   * \param[in]  d2vdx2     - 2nd x-derivative of the v-velocity.
   * \param[in]  d2vdy2     - 2nd y-derivative of the v-velocity.
   * \param[in]  d2vdz2     - 2nd z-derivative of the v-velocity.
   * \param[in]  d2vdxdy    - x-y cross-derivative of the v-velocity.
   * \param[in]  d2vdxdz    - x-z cross-derivative of the v-velocity.
   * \param[in]  d2vdydz    - y-z cross-derivative of the v-velocity.
   * \param[in]  d2wdx2     - 2nd x-derivative of the w-velocity.
   * \param[in]  d2wdy2     - 2nd y-derivative of the w-velocity.
   * \param[in]  d2wdz2     - 2nd z-derivative of the w-velocity.
   * \param[in]  d2wdxdy    - x-y cross-derivative of the w-velocity.
   * \param[in]  d2wdxdz    - x-z cross-derivative of the w-velocity.
   * \param[in]  d2wdydz    - y-z cross-derivative of the w-velocity.
   * \param[in]  lenScale   - Length scale of the corresponding element.
   * \param[in]  distToWall - Distance to the nearest wall.
   * \param[out] dMuTdx     - x-derivative of the turbulent viscosity.
   * \param[out] dMuTdy     - y-derivative of the turbulent viscosity.
   * \param[out] dMuTdz     - z-derivative of the turbulent viscosity.
   */
  virtual void ComputeGradEddyViscosity_3D(const su2double rho,
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
                                                 su2double &dMuTdz);
};

/*!
 * \class CSmagorinskyModel
 * \brief Derived class for defining the Smagorinsky SGS model.
 * \author: E. van der Weide, T. Economon, P. Urbanczyk
 * \version 8.0.0 "Harrier"
 */
class CSmagorinskyModel : public CSGSModel {

public:

  su2double const_smag;  /*!< \brief Smagorinsky Constant C_s.  */
  su2double filter_mult; /*!< \brief Multiplier to get filter width from grid length scale. */
  /*!
   * \brief Constructor of the class.
   */
  CSmagorinskyModel(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CSmagorinskyModel(void) override;

  /*!
   * \brief Function to determine the eddy viscosity for
            the given function arguments for a 2D simulation.
   * \param[in] rho        - Density
   * \param[in] dudx       - x-derivative of the u-velocity.
   * \param[in] dudy       - y-derivative of the u-velocity.
   * \param[in] dvdx       - x-derivative of the v-velocity.
   * \param[in] dvdy       - y-derivative of the v-velocity.
   * \param[in] lenScale   - Length scale of the corresponding element.
   * \param[in] distToWall - Distance to the nearest wall.
   * \return Value of the dynamic eddy viscosity for the Smagorinsky model.
   */
  su2double ComputeEddyViscosity_2D(const su2double rho,
                                    const su2double dudx,
                                    const su2double dudy,
                                    const su2double dvdx,
                                    const su2double dvdy,
                                    const su2double lenScale,
                                    const su2double distToWall) override;

/*!
   * \brief Function to determine the eddy viscosity for
            the given function arguments for a 3D simulation.
   * \param[in] rho        - Density
   * \param[in] dudx       - x-derivative of the u-velocity.
   * \param[in] dudy       - y-derivative of the u-velocity.
   * \param[in] dudz       - z-derivative of the u-velocity.
   * \param[in] dvdx       - x-derivative of the v-velocity.
   * \param[in] dvdy       - y-derivative of the v-velocity.
   * \param[in] dvdz       - z-derivative of the v-velocity.
   * \param[in] dwdx       - x-derivative of the w-velocity.
   * \param[in] dwdy       - y-derivative of the w-velocity.
   * \param[in] dwdz       - z-derivative of the w-velocity.
   * \param[in] lenScale   - Length scale of the corresponding element.
   * \param[in] distToWall - Distance to the nearest wall.
   * \return Value of the dynamic eddy viscosity for the Smagorinsky model.
   */
  su2double ComputeEddyViscosity_3D(const su2double rho,
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
                                    const su2double distToWall) override;

  /*!
   * \brief Function to determine the gradients of the eddy viscosity
            for the given function arguments for a 2D simulation.
   * \param[in]  rho        - Density.
   * \param[in]  drhodx     - x-derivative of the density.
   * \param[in]  drhody     - y-derivative of the density.
   * \param[in]  dudx       - x-derivative of the u-velocity.
   * \param[in]  dudy       - y-derivative of the u-velocity.
   * \param[in]  dvdx       - x-derivative of the v-velocity.
   * \param[in]  dvdy       - y-derivative of the v-velocity.
   * \param[in]  d2udx2     - 2nd x-derivative of the u-velocity.
   * \param[in]  d2udy2     - 2nd y-derivative of the u-velocity.
   * \param[in]  d2udxdy    - x-y cross-derivative of the u-velocity.
   * \param[in]  d2vdx2     - 2nd x-derivative of the v-velocity.
   * \param[in]  d2vdy2     - 2nd y-derivative of the v-velocity.
   * \param[in]  d2vdxdy    - x-y cross-derivative of the v-velocity.
   * \param[in]  lenScale   - Length scale of the corresponding element.
   * \param[in]  distToWall - Distance to the nearest wall.
   * \param[out] dMuTdx     - x-derivative of the turbulent viscosity.
   * \param[out] dMuTdy     - y-derivative of the turbulent viscosity.
   */
  void ComputeGradEddyViscosity_2D(const su2double rho,
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
                                         su2double &dMuTdy) override;

  /*!
   * \brief function to determine the gradients of the eddy viscosity
            for the given function arguments for a 3D simulation.
   * \param[in]  rho        - Density.
   * \param[in]  drhodx     - x-derivative of the density.
   * \param[in]  drhody     - y-derivative of the density.
   * \param[in]  drhodz     - z-derivative of the density.
   * \param[in]  dudx       - x-derivative of the u-velocity.
   * \param[in]  dudy       - y-derivative of the u-velocity.
   * \param[in]  dudz       - z-derivative of the u-velocity.
   * \param[in]  dvdx       - x-derivative of the v-velocity.
   * \param[in]  dvdy       - y-derivative of the v-velocity.
   * \param[in]  dvdz       - z-derivative of the v-velocity.
   * \param[in]  dwdx       - x-derivative of the w-velocity.
   * \param[in]  dwdy       - y-derivative of the w-velocity.
   * \param[in]  dwdz       - z-derivative of the w-velocity.
   * \param[in]  d2udx2     - 2nd x-derivative of the u-velocity.
   * \param[in]  d2udy2     - 2nd y-derivative of the u-velocity.
   * \param[in]  d2udz2     - 2nd z-derivative of the u-velocity.
   * \param[in]  d2udxdy    - x-y cross-derivative of the u-velocity.
   * \param[in]  d2udxdz    - x-z cross-derivative of the u-velocity.
   * \param[in]  d2udydz    - y-z cross-derivative of the u-velocity.
   * \param[in]  d2vdx2     - 2nd x-derivative of the v-velocity.
   * \param[in]  d2vdy2     - 2nd y-derivative of the v-velocity.
   * \param[in]  d2vdz2     - 2nd z-derivative of the v-velocity.
   * \param[in]  d2vdxdy    - x-y cross-derivative of the v-velocity.
   * \param[in]  d2vdxdz    - x-z cross-derivative of the v-velocity.
   * \param[in]  d2vdydz    - y-z cross-derivative of the v-velocity.
   * \param[in]  d2wdx2     - 2nd x-derivative of the w-velocity.
   * \param[in]  d2wdy2     - 2nd y-derivative of the w-velocity.
   * \param[in]  d2wdz2     - 2nd z-derivative of the w-velocity.
   * \param[in]  d2wdxdy    - x-y cross-derivative of the w-velocity.
   * \param[in]  d2wdxdz    - x-z cross-derivative of the w-velocity.
   * \param[in]  d2wdydz    - y-z cross-derivative of the w-velocity.
   * \param[in]  lenScale   - Length scale of the corresponding element.
   * \param[in]  distToWall - Distance to the nearest wall.
   * \param[out] dMuTdx     - x-derivative of the turbulent viscosity.
   * \param[out] dMuTdy     - y-derivative of the turbulent viscosity.
   * \param[out] dMuTdz     - z-derivative of the turbulent viscosity.
   */
  void ComputeGradEddyViscosity_3D(const su2double rho,
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
                                         su2double &dMuTdz) override;
};

/*!
 * \class CWALEModel
 * \brief Derived class for defining the WALE SGS model.
 * \author: E. van der Weide, T. Economon, P. Urbanczyk
 * \version 8.0.0 "Harrier"
 */
class CWALEModel : public CSGSModel {

public:
  su2double const_WALE; /*!< \brief WALE Constant Cw.  */

  /*!
   * \brief Constructor of the class.
   */
  CWALEModel(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CWALEModel(void) override;

  /*!
   * \brief Function to determine the eddy viscosity for
            the given function arguments for a 2D simulation.
   * \param[in] rho        - Density
   * \param[in] dudx       - x-derivative of the u-velocity.
   * \param[in] dudy       - y-derivative of the u-velocity.
   * \param[in] dvdx       - x-derivative of the v-velocity.
   * \param[in] dvdy       - y-derivative of the v-velocity.
   * \param[in] lenScale   - Length scale of the corresponding element.
   * \param[in] distToWall - Distance to the nearest wall.
   * \return Value of the dynamic eddy viscosity for the WALE model.
   */
  su2double ComputeEddyViscosity_2D(const su2double rho,
                                    const su2double dudx,
                                    const su2double dudy,
                                    const su2double dvdx,
                                    const su2double dvdy,
                                    const su2double lenScale,
                                    const su2double distToWall) override;

/*!
   * \brief Function to determine the eddy viscosity for
            the given function arguments for a 3D simulation.
   * \param[in] rho        - Density
   * \param[in] dudx       - x-derivative of the u-velocity.
   * \param[in] dudy       - y-derivative of the u-velocity.
   * \param[in] dudz       - z-derivative of the u-velocity.
   * \param[in] dvdx       - x-derivative of the v-velocity.
   * \param[in] dvdy       - y-derivative of the v-velocity.
   * \param[in] dvdz       - z-derivative of the v-velocity.
   * \param[in] dwdx       - x-derivative of the w-velocity.
   * \param[in] dwdy       - y-derivative of the w-velocity.
   * \param[in] dwdz       - z-derivative of the w-velocity.
   * \param[in] lenScale   - Length scale of the corresponding element.
   * \param[in] distToWall - Distance to the nearest wall.
   * \return Value of the dynamic eddy viscosity for the WALE model.
   */
  su2double ComputeEddyViscosity_3D(const su2double rho,
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
                                    const su2double distToWall) override;

  /*!
   * \brief Function to determine the gradients of the eddy viscosity
            for the given function arguments for a 2D simulation.
   * \param[in]  rho        - Density.
   * \param[in]  drhodx     - x-derivative of the density.
   * \param[in]  drhody     - y-derivative of the density.
   * \param[in]  dudx       - x-derivative of the u-velocity.
   * \param[in]  dudy       - y-derivative of the u-velocity.
   * \param[in]  dvdx       - x-derivative of the v-velocity.
   * \param[in]  dvdy       - y-derivative of the v-velocity.
   * \param[in]  d2udx2     - 2nd x-derivative of the u-velocity.
   * \param[in]  d2udy2     - 2nd y-derivative of the u-velocity.
   * \param[in]  d2udxdy    - x-y cross-derivative of the u-velocity.
   * \param[in]  d2vdx2     - 2nd x-derivative of the v-velocity.
   * \param[in]  d2vdy2     - 2nd y-derivative of the v-velocity.
   * \param[in]  d2vdxdy    - x-y cross-derivative of the v-velocity.
   * \param[in]  lenScale   - Length scale of the corresponding element.
   * \param[in]  distToWall - Distance to the nearest wall.
   * \param[out] dMuTdx     - x-derivative of the turbulent viscosity.
   * \param[out] dMuTdy     - y-derivative of the turbulent viscosity.
   */
  void ComputeGradEddyViscosity_2D(const su2double rho,
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
                                         su2double &dMuTdy) override;

  /*!
   * \brief function to determine the gradients of the eddy viscosity
            for the given function arguments for a 3D simulation.
   * \param[in]  rho        - Density.
   * \param[in]  drhodx     - x-derivative of the density.
   * \param[in]  drhody     - y-derivative of the density.
   * \param[in]  drhodz     - z-derivative of the density.
   * \param[in]  dudx       - x-derivative of the u-velocity.
   * \param[in]  dudy       - y-derivative of the u-velocity.
   * \param[in]  dudz       - z-derivative of the u-velocity.
   * \param[in]  dvdx       - x-derivative of the v-velocity.
   * \param[in]  dvdy       - y-derivative of the v-velocity.
   * \param[in]  dvdz       - z-derivative of the v-velocity.
   * \param[in]  dwdx       - x-derivative of the w-velocity.
   * \param[in]  dwdy       - y-derivative of the w-velocity.
   * \param[in]  dwdz       - z-derivative of the w-velocity.
   * \param[in]  d2udx2     - 2nd x-derivative of the u-velocity.
   * \param[in]  d2udy2     - 2nd y-derivative of the u-velocity.
   * \param[in]  d2udz2     - 2nd z-derivative of the u-velocity.
   * \param[in]  d2udxdy    - x-y cross-derivative of the u-velocity.
   * \param[in]  d2udxdz    - x-z cross-derivative of the u-velocity.
   * \param[in]  d2udydz    - y-z cross-derivative of the u-velocity.
   * \param[in]  d2vdx2     - 2nd x-derivative of the v-velocity.
   * \param[in]  d2vdy2     - 2nd y-derivative of the v-velocity.
   * \param[in]  d2vdz2     - 2nd z-derivative of the v-velocity.
   * \param[in]  d2vdxdy    - x-y cross-derivative of the v-velocity.
   * \param[in]  d2vdxdz    - x-z cross-derivative of the v-velocity.
   * \param[in]  d2vdydz    - y-z cross-derivative of the v-velocity.
   * \param[in]  d2wdx2     - 2nd x-derivative of the w-velocity.
   * \param[in]  d2wdy2     - 2nd y-derivative of the w-velocity.
   * \param[in]  d2wdz2     - 2nd z-derivative of the w-velocity.
   * \param[in]  d2wdxdy    - x-y cross-derivative of the w-velocity.
   * \param[in]  d2wdxdz    - x-z cross-derivative of the w-velocity.
   * \param[in]  d2wdydz    - y-z cross-derivative of the w-velocity.
   * \param[in]  lenScale   - Length scale of the corresponding element.
   * \param[in]  distToWall - Distance to the nearest wall.
   * \param[out] dMuTdx     - x-derivative of the turbulent viscosity.
   * \param[out] dMuTdy     - y-derivative of the turbulent viscosity.
   * \param[out] dMuTdz     - z-derivative of the turbulent viscosity.
   */
  void ComputeGradEddyViscosity_3D(const su2double rho,
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
                                         su2double &dMuTdz) override;
};

/*!
 * \class CVremanModel
 * \brief Derived class for defining the WALE SGS model.
 * \author: E. van der Weide, T. Economon, P. Urbanczyk, E. Molina
 * \version 8.0.0 "Harrier"
 */
class CVremanModel : public CSGSModel {

public:
  su2double const_Vreman; /*!< \brief Vreman Constant c=2.5*Cs*Cs.  */

  /*!
   * \brief Constructor of the class.
   */
  CVremanModel(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CVremanModel(void) override;

  /*!
   * \brief Function to determine the eddy viscosity for
   the given function arguments for a 2D simulation.
   * \param[in] rho        - Density
   * \param[in] dudx       - x-derivative of the u-velocity.
   * \param[in] dudy       - y-derivative of the u-velocity.
   * \param[in] dvdx       - x-derivative of the v-velocity.
   * \param[in] dvdy       - y-derivative of the v-velocity.
   * \param[in] lenScale   - Length scale of the corresponding element.
   * \param[in] distToWall - Distance to the nearest wall.
   * \return Value of the dynamic eddy viscosity for the WALE model.
   */
  su2double ComputeEddyViscosity_2D(const su2double rho,
                                    const su2double dudx,
                                    const su2double dudy,
                                    const su2double dvdx,
                                    const su2double dvdy,
                                    const su2double lenScale,
                                    const su2double distToWall) override;

  /*!
   * \brief Function to determine the eddy viscosity for
   the given function arguments for a 3D simulation.
   * \param[in] rho        - Density
   * \param[in] dudx       - x-derivative of the u-velocity.
   * \param[in] dudy       - y-derivative of the u-velocity.
   * \param[in] dudz       - z-derivative of the u-velocity.
   * \param[in] dvdx       - x-derivative of the v-velocity.
   * \param[in] dvdy       - y-derivative of the v-velocity.
   * \param[in] dvdz       - z-derivative of the v-velocity.
   * \param[in] dwdx       - x-derivative of the w-velocity.
   * \param[in] dwdy       - y-derivative of the w-velocity.
   * \param[in] dwdz       - z-derivative of the w-velocity.
   * \param[in] lenScale   - Length scale of the corresponding element.
   * \param[in] distToWall - Distance to the nearest wall.
   * \return Value of the dynamic eddy viscosity for the WALE model.
   */
  su2double ComputeEddyViscosity_3D(const su2double rho,
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
                                    const su2double distToWall) override;

  /*!
   * \brief Function to determine the gradients of the eddy viscosity
   for the given function arguments for a 2D simulation.
   * \param[in]  rho        - Density.
   * \param[in]  drhodx     - x-derivative of the density.
   * \param[in]  drhody     - y-derivative of the density.
   * \param[in]  dudx       - x-derivative of the u-velocity.
   * \param[in]  dudy       - y-derivative of the u-velocity.
   * \param[in]  dvdx       - x-derivative of the v-velocity.
   * \param[in]  dvdy       - y-derivative of the v-velocity.
   * \param[in]  d2udx2     - 2nd x-derivative of the u-velocity.
   * \param[in]  d2udy2     - 2nd y-derivative of the u-velocity.
   * \param[in]  d2udxdy    - x-y cross-derivative of the u-velocity.
   * \param[in]  d2vdx2     - 2nd x-derivative of the v-velocity.
   * \param[in]  d2vdy2     - 2nd y-derivative of the v-velocity.
   * \param[in]  d2vdxdy    - x-y cross-derivative of the v-velocity.
   * \param[in]  lenScale   - Length scale of the corresponding element.
   * \param[in]  distToWall - Distance to the nearest wall.
   * \param[out] dMuTdx     - x-derivative of the turbulent viscosity.
   * \param[out] dMuTdy     - y-derivative of the turbulent viscosity.
   */
  void ComputeGradEddyViscosity_2D(const su2double rho,
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
                                   su2double &dMuTdy) override;

  /*!
   * \brief function to determine the gradients of the eddy viscosity
   for the given function arguments for a 3D simulation.
   * \param[in]  rho        - Density.
   * \param[in]  drhodx     - x-derivative of the density.
   * \param[in]  drhody     - y-derivative of the density.
   * \param[in]  drhodz     - z-derivative of the density.
   * \param[in]  dudx       - x-derivative of the u-velocity.
   * \param[in]  dudy       - y-derivative of the u-velocity.
   * \param[in]  dudz       - z-derivative of the u-velocity.
   * \param[in]  dvdx       - x-derivative of the v-velocity.
   * \param[in]  dvdy       - y-derivative of the v-velocity.
   * \param[in]  dvdz       - z-derivative of the v-velocity.
   * \param[in]  dwdx       - x-derivative of the w-velocity.
   * \param[in]  dwdy       - y-derivative of the w-velocity.
   * \param[in]  dwdz       - z-derivative of the w-velocity.
   * \param[in]  d2udx2     - 2nd x-derivative of the u-velocity.
   * \param[in]  d2udy2     - 2nd y-derivative of the u-velocity.
   * \param[in]  d2udz2     - 2nd z-derivative of the u-velocity.
   * \param[in]  d2udxdy    - x-y cross-derivative of the u-velocity.
   * \param[in]  d2udxdz    - x-z cross-derivative of the u-velocity.
   * \param[in]  d2udydz    - y-z cross-derivative of the u-velocity.
   * \param[in]  d2vdx2     - 2nd x-derivative of the v-velocity.
   * \param[in]  d2vdy2     - 2nd y-derivative of the v-velocity.
   * \param[in]  d2vdz2     - 2nd z-derivative of the v-velocity.
   * \param[in]  d2vdxdy    - x-y cross-derivative of the v-velocity.
   * \param[in]  d2vdxdz    - x-z cross-derivative of the v-velocity.
   * \param[in]  d2vdydz    - y-z cross-derivative of the v-velocity.
   * \param[in]  d2wdx2     - 2nd x-derivative of the w-velocity.
   * \param[in]  d2wdy2     - 2nd y-derivative of the w-velocity.
   * \param[in]  d2wdz2     - 2nd z-derivative of the w-velocity.
   * \param[in]  d2wdxdy    - x-y cross-derivative of the w-velocity.
   * \param[in]  d2wdxdz    - x-z cross-derivative of the w-velocity.
   * \param[in]  d2wdydz    - y-z cross-derivative of the w-velocity.
   * \param[in]  lenScale   - Length scale of the corresponding element.
   * \param[in]  distToWall - Distance to the nearest wall.
   * \param[out] dMuTdx     - x-derivative of the turbulent viscosity.
   * \param[out] dMuTdy     - y-derivative of the turbulent viscosity.
   * \param[out] dMuTdz     - z-derivative of the turbulent viscosity.
   */
  void ComputeGradEddyViscosity_3D(const su2double rho,
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
                                   su2double &dMuTdz) override;
};
#include "sgs_model.inl"
