/*!
 * \file sgs_model.hpp
 * \brief Headers of the LES subgrid scale models of the SU2 solvers.
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

#include "../../Common/include/mpi_structure.hpp"

#include <iostream>
#include <cmath>

using namespace std;

/*!
 * \class CSGSModel
 * \brief Base class for defining the LES subgrid scale model.
 * \author: E. van der Weide, T. Economon, P. Urbanczyk
 * \version 5.0.0 "Raven"
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
            for the given function arguments.
   * \param[in] nDim       - Number of dimensions of the problem.
   * \param[in] rho        - Density  
   * \param[in] velGrad    - Cartesian gradients of the velocity components.
                             valGrad[i][j] = du_i/dx_j.
   * \param[in] lenScale   - Length scale of the corresponding element.
   * \param[in] distToWall - Distance to the nearest wall.
   * \return Value of the dynamic eddy viscosity. For the base class 0 is returned.
   */
  virtual su2double ComputeEddyViscosity(const unsigned short nDim,
                                         const su2double      rho,
                                         const su2double      velGrad[3][3],
                                         const su2double      lenScale,
                                         const su2double      distToWall);

  /*!
   * \brief Virtual function to determine the gradients of the eddy viscosity
            for the given function arguments.
   * \param[in]  nDim              - Number of dimensions of the problem.
   * \param[in]  rho               - Density.
   * \param[in]  rhoGrad           - Cartesian gradients of the density.
   * \param[in]  velGrad           - Cartesian gradients of the velocity components.
                                     valGrad[i][j] = du_i/dx_j.
   * \param[in]  velHess           - Cartesian Hessian of the velocity components.
                                     valHess[i][j][k] = d2u_i/(dx_j dx_k).
   * \param[in]  lenScale          - Length scale of the corresponding element.
   * \param[in]  distToWall        - Distance to the nearest wall.
   * \param[out] ViscosityTurbGrad - Cartesian gradient of the turbulent viscosity.
   */
  virtual void ComputeGradEddyViscosity(const unsigned short nDim,
                                        const su2double      rho,
                                        const su2double      rhoGrad[3],
                                        const su2double      velGrad[3][3],
                                        const su2double      velHess[3][3][3],
                                        const su2double      lenScale,
                                        const su2double      distToWall,
                                              su2double      ViscosityTurbGrad[3]);
};

/*!
 * \class CSmagorinskyModel
 * \brief Derived class for defining the Smagorinsky SGS model.
 * \author: E. van der Weide, T. Economon, P. Urbanczyk
 * \version 5.0.0 "Raven"
 */
class CSmagorinskyModel : public CSGSModel {

public:

  su2double const_smag; /*!< \brief Smagorinsky Constant C_s.  */
  su2double filter_mult; /*!< \brief Multiplier to get filter width from grid length scale. */
  /*!
   * \brief Constructor of the class.
   */
  CSmagorinskyModel(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CSmagorinskyModel(void);

  /*!
   * \brief Function to determine the eddy viscosity for
            the given function arguments.
   * \param[in] nDim       - Number of dimensions of the problem.
   * \param[in] rho        - Density  
   * \param[in] velGrad    - Cartesian gradients of the velocity components.
                             valGrad[i][j] = du_i/dx_j.
   * \param[in] lenScale   - Length scale of the corresponding element.
   * \param[in] distToWall - Distance to the nearest wall.
   * \return Value of the dynamic eddy viscosity for the Smagorinsky model.
   */
  su2double ComputeEddyViscosity(const unsigned short nDim,
                                 const su2double      rho,
                                 const su2double      velGrad[3][3],
                                 const su2double      lenScale,
                                 const su2double      distToWall);

  /*!
   * \brief Function to determine the gradients of the eddy viscosity
            for the given function arguments.
   * \param[in]  nDim              - Number of dimensions of the problem.
   * \param[in]  rho               - Density.
   * \param[in]  rhoGrad           - Cartesian gradients of the density.
   * \param[in]  velGrad           - Cartesian gradients of the velocity components.
                                     valGrad[i][j] = du_i/dx_j.
   * \param[in]  velHess           - Cartesian Hessian of the velocity components.
                                     valHess[i][j][k] = d2u_i/(dx_j dx_k).
   * \param[in]  lenScale          - Length scale of the corresponding element.
   * \param[in]  distToWall        - Distance to the nearest wall.
   * \param[out] ViscosityTurbGrad - Cartesian gradient of the turbulent viscosity.
   */
  void ComputeGradEddyViscosity(const unsigned short nDim,
                                const su2double      rho,
                                const su2double      rhoGrad[3],
                                const su2double      velGrad[3][3],
                                const su2double      velHess[3][3][3],
                                const su2double      lenScale,
                                const su2double      distToWall,
                                      su2double      ViscosityTurbGrad[3]);
};

/*!
 * \class CWALEModel
 * \brief Derived class for defining the WALE SGS model.
 * \author: E. van der Weide, T. Economon, P. Urbanczyk
 * \version 5.0.0 "Raven"
 */
class CWALEModel : public CSGSModel {

public:
  /*!
   * \brief Constructor of the class.
   */
  CWALEModel(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CWALEModel(void);

  /*!
   * \brief Function to determine the eddy viscosity for
            the given function arguments.
   * \param[in] nDim       - Number of dimensions of the problem.
   * \param[in] rho        - Density  
   * \param[in] velGrad    - Cartesian gradients of the velocity components.
                             valGrad[i][j] = du_i/dx_j.
   * \param[in] lenScale   - Length scale of the corresponding element.
   * \param[in] distToWall - Distance to the nearest wall.
   * \return Value of the dynamic eddy viscosity for the WALE model.
   */
  su2double ComputeEddyViscosity(const unsigned short nDim,
                                 const su2double      rho,
                                 const su2double      velGrad[3][3],
                                 const su2double      lenScale,
                                 const su2double      distToWall);

  /*!
   * \brief Function to determine the gradients of the eddy viscosity
            for the given function arguments.
   * \param[in]  nDim              - Number of dimensions of the problem.
   * \param[in]  rho               - Density.
   * \param[in]  rhoGrad           - Cartesian gradients of the density.
   * \param[in]  velGrad           - Cartesian gradients of the velocity components.
                                     valGrad[i][j] = du_i/dx_j.
   * \param[in]  velHess           - Cartesian Hessian of the velocity components.
                                     valHess[i][j][k] = d2u_i/(dx_j dx_k).
   * \param[in]  lenScale          - Length scale of the corresponding element.
   * \param[in]  distToWall        - Distance to the nearest wall.
   * \param[out] ViscosityTurbGrad - Cartesian gradient of the turbulent viscosity.
   */
  void ComputeGradEddyViscosity(const unsigned short nDim,
                                const su2double      rho,
                                const su2double      rhoGrad[3],
                                const su2double      velGrad[3][3],
                                const su2double      velHess[3][3][3],
                                const su2double      lenScale,
                                const su2double      distToWall,
                                      su2double      ViscosityTurbGrad[3]);
};

#include "sgs_model.inl"
