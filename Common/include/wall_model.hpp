/*!
 * \file wall_model.hpp
 * \brief Headers for the wall model functions for large eddy simulations.
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

#include "./parallelization/mpi_structure.hpp"
#include "./CConfig.hpp"

#include <iostream>
#include <cmath>

using namespace std;

/* Forward declaration of the class CFluidModel to avoid any problems. */
class CFluidModel;

/*!
 * \class CWallModel
 * \brief Base class for defining the LES wall model.
 * \author: E. van der Weide, T. Economon, P. Urbanczyk
 * \version 8.0.0 "Harrier"
 */
class CWallModel {
 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CWallModel(CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CWallModel(void) = default;

  /*!
   * \brief Virtual function, which computes the wall shear stress and heat flux
            from the data at the exchange location.
   * \param[in]  tExchange              - Temperature at the exchange location.
   * \param[in]  velExchange            - Velocity at the exchange location.
   * \param[in]  muExchange             - Laminar viscosity at the exchange location.
   * \param[in]  pExchange              - Pressure at the exchange location.
   * \param[in]  Wall_HeatFlux          - Value of the wall heat flux, if prescribed.
   * \param[in]  HeatFlux_Prescribed    - Whether or not the wall heat flux is prescribed.
   * \param[in]  Wall_Temperature       - Value of the wall temperature, if prescribed.
   * \param[in]  Temperature_Prescribed - Wheter or not the wall temperature is prescribed.
   * \param[in]  FluidModel             - Fluid model used in the solver.
   * \param[out] tauWall                - Wall shear stress, to be computed.
   * \param[out] qWall                  - Wall hear flux, to be computed (if not prescribed).
   * \param[out] ViscosityWall          - Laminar viscosity at the wall, to be computed.
   * \param[out] OverCvWall             - Thermal conductivity divided by Cv at the wall,
                                          to be computed.
   */
  virtual void WallShearStressAndHeatFlux(const su2double tExchange, const su2double velExchange,
                                          const su2double muExchange, const su2double pExchange,
                                          const su2double Wall_HeatFlux, const bool HeatFlux_Prescribed,
                                          const su2double TWall, const bool Temperature_Prescribed,
                                          CFluidModel* FluidModel, su2double& tauWall, su2double& qWall,
                                          su2double& ViscosityWall, su2double& kOverCvWall);

 protected:
  su2double h_wm;    /*!< \brief The thickness of the wall model. This is also basically the exchange location */
  su2double Pr_lam;  /*!< \brief Laminar Prandtl number. */
  su2double Pr_turb; /*!< \brief Turbulent Prandtl number. */
  su2double karman;  /*!< \brief von Karman constant. */

 private:
  /*!
   * \brief Default constructor of the class, disabled.
   */
  CWallModel(void);
};

class CWallModel1DEQ : public CWallModel {
 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] config     - Definition of the particular problem.
   * \param[in] Marker_Tag - String, which identifies the boundary marker for
                             which the wall model is used.
   */
  CWallModel1DEQ(CConfig* config, const string& Marker_Tag);

  /*!
   * \brief Function, which computes the wall shear stress and heat flux
            from the data at the exchange location.
   * \param[in]  tExchange              - Temperature at the exchange location.
   * \param[in]  velExchange            - Velocity at the exchange location.
   * \param[in]  muExchange             - Laminar viscosity at the exchange location.
   * \param[in]  pExchange              - Pressure at the exchange location.
   * \param[in]  Wall_HeatFlux          - Value of the wall heat flux, if prescribed.
   * \param[in]  HeatFlux_Prescribed    - Whether or not the wall heat flux is prescribed.
   * \param[in]  Wall_Temperature       - Value of the wall temperature, if prescribed.
   * \param[in]  Temperature_Prescribed - Wheter or not the wall temperature is prescribed.
   * \param[in]  FluidModel             - Fluid model used in the solver.
   * \param[out] tauWall                - Wall shear stress, to be computed.
   * \param[out] qWall                  - Wall hear flux, to be computed (if not prescribed).
   * \param[out] ViscosityWall          - Laminar viscosity at the wall, to be computed.
   * \param[out] kOverCvWall             - Thermal conductivity divided by Cv at the wall,
                                          to be computed.
   */
  void WallShearStressAndHeatFlux(const su2double tExchange, const su2double velExchange, const su2double muExchange,
                                  const su2double pExchange, const su2double Wall_HeatFlux,
                                  const bool HeatFlux_Prescribed, const su2double Wall_Temperature,
                                  const bool Temperature_Prescribed, CFluidModel* FluidModel, su2double& tauWall,
                                  su2double& qWall, su2double& ViscosityWall, su2double& kOverCvWall) override;

 private:
  su2double expansionRatio; /*!< \brief  Stretching factor used for the wall model grid. */
  int numPoints;            /*!< \brief  Number of points used in the wall model grid. */

  vector<su2double> y_cv; /*!< \brief  The coordinates in normal direction of the wall model grid (control volumes). */
  vector<su2double> y_fa; /*!< \brief  The coordinates in normal direction of the wall model grid (faces of CV). */

  /*!
   * \brief Default constructor of the class, disabled.
   */
  CWallModel1DEQ(void);
};

class CWallModelLogLaw : public CWallModel {
 public:
  /*!
   * \brief Constructor of the class, which initializes the object.
   * \param[in] config     - Definition of the particular problem.
   * \param[in] Marker_Tag - String, which identifies the boundary marker for
                             which the wall model is used.
   */
  CWallModelLogLaw(CConfig* config, const string& Marker_Tag);

  /*!
   * \brief Function, which computes the wall shear stress and heat flux
   from the data at the exchange location.
   * \param[in]  tExchange              - Temperature at the exchange location.
   * \param[in]  velExchange            - Velocity at the exchange location.
   * \param[in]  muExchange             - Laminar Viscosity at the exchange location.
   * \param[in]  pExchange              - Pressure at the exchange location.
   * \param[in]  Wall_HeatFlux          - Value of the wall heat flux, if prescribed.
   * \param[in]  HeatFlux_Prescribed    - Whether or not the wall heat flux is prescribed.
   * \param[in]  Wall_Temperature       - Value of the wall temperature, if prescribed.
   * \param[in]  Temperature_Prescribed - Wheter or not the wall temperature is prescribed.
   * \param[in]  FluidModel             - Fluid model used in the solver.
   * \param[out] tauWall                - Wall shear stress, to be computed.
   * \param[out] qWall                  - Wall hear flux, to be computed (if not prescribed).
   * \param[out] ViscosityWall          - Laminar viscosity at the wall, to be computed.
   * \param[out] kOverCvWall             - Thermal conductivity divided by Cv at the wall,
   to be computed.
   */
  void WallShearStressAndHeatFlux(const su2double tExchange, const su2double velExchange, const su2double muExchange,
                                  const su2double pExchange, const su2double Wall_HeatFlux,
                                  const bool HeatFlux_Prescribed, const su2double Wall_Temperature,
                                  const bool Temperature_Prescribed, CFluidModel* FluidModel, su2double& tauWall,
                                  su2double& qWall, su2double& ViscosityWall, su2double& kOverCvWall) override;

 private:
  su2double C; /*!< \brief Constant to match the Reichardt BL profile. */

  /*!
   * \brief Default constructor of the class, disabled.
   */
  CWallModelLogLaw(void);
};
