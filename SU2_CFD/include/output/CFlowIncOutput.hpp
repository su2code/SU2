/*!
 * \file CFlowIncOutput.hpp
 * \brief  Headers of the incompressible flow output.
 * \author T. Albring, R. Sanchez
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

#include "CFlowOutput.hpp"

class CVariable;

/*! \class CFlowIncOutput
 *  \brief Output class for incompressible flow problems.
 *  \author R. Sanchez, T. Albring.
 *  \date May 30, 2018.
 */
class CFlowIncOutput final: public CFlowOutput {
private:
  TURB_MODEL turb_model;     /*!< \brief The kind of turbulence model*/
  bool heat;                 /*!< \brief Boolean indicating whether have a heat problem*/
  bool weakly_coupled_heat;  /*!< \brief Boolean indicating whether have a weakly coupled heat equation*/
  bool flamelet;  /*!< \brief Boolean indicating whether we solve the flamelet equations */
  unsigned short streamwisePeriodic;   /*!< \brief Boolean indicating whether it is a streamwise periodic simulation. */
  bool streamwisePeriodic_temperature; /*!< \brief Boolean indicating streamwise periodic temperature is used. */

public:
  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CFlowIncOutput(CConfig *config, unsigned short nDim);

  /*!
   * \brief Load the history output field values
   * \param[in] config - Definition of the particular problem.
   */
  void LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) override;

  /*!
   * \brief Set the available volume output fields
   * \param[in] config - Definition of the particular problem.
   */
  void SetVolumeOutputFields(CConfig *config) override;

  /*!
   * \brief Set the values of the volume output fields for a point.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The container holding all solution data.
   * \param[in] iPoint - Index of the point.
   */
  void LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint) override;

  /*!
   * \brief Set the available history output fields
   * \param[in] config - Definition of the particular problem.
   */
  void SetHistoryOutputFields(CConfig *config) override;

  /*!
   * \brief Check whether the base values for relative residuals should be initialized
   * \param[in] config - Definition of the particular problem.
   * \return <TRUE> if the residuals should be initialized.
   */
  bool SetInitResiduals(const CConfig *config) override ;

};
