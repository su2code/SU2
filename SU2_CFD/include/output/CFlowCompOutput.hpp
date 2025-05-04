/*!
 * \file CFlowCompOutput.hpp
 * \brief  Headers of the compressible flow output.
 * \author R. Sanchez, T. Albring.
 * \version 8.2.0 "Harrier"
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

#pragma once

#include "CFlowOutput.hpp"

class CVariable;

/*! \class CFlowCompOutput
 *  \brief Output class for compressible flow problems.
 *  \author R. Sanchez, T. Albring.
 *  \date May 30, 2018.
 */
class CFlowCompOutput final: public CFlowOutput {
private:
  TURB_MODEL turb_model; //!< Kind of turbulence model

public:
  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CFlowCompOutput(const CConfig *config, unsigned short nDim);

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

  /*!
   * \brief Write any additional output defined for the current solver.
   * \param[in] config - Definition of the particular problem per zone.
   */
  void SetAdditionalScreenOutput(const CConfig *config) override;

  /*!
   * \brief Determines if the history file output.
   * \param[in] config - Definition of the particular problem.
   */
  bool WriteHistoryFileOutput(const CConfig *config) override ;

  /*!
   * \brief Sets the turboperformance screen output
   * \param[in] TurboPerf - Turboperformance class 
   * \param[in] config - Definition of the particular problem
   * \param[in] TimeIter - Index of the current time-step
   * \param[in] OuterIter - Index of current outer iteration
   * \param[in] InnerIter - Index of current inner iteration
   */
  void SetTurboPerformance_Output(std::shared_ptr<CTurboOutput> TurboPerf, CConfig *config, unsigned long TimeIter, unsigned long OuterIter, unsigned long InnerIter) override;

  /*!
   * \brief Sets the multizone turboperformacne screen output
   * \param[in] TurboStagePerf - Stage turboperformance class
   * \param[in] TurboPerf - Turboperformance class
   * \param[in] config - Definition of the particular problem
   */
  void SetTurboMultiZonePerformance_Output(std::shared_ptr<CTurbomachineryStagePerformance> TurboStagePerf, std::shared_ptr<CTurboOutput> TurboPerf, CConfig *config) override;
  
  /*!
   * \brief Loads the turboperformacne history data
   * \param[in] TurboStagePerf - Stage turboperformance class
   * \param[in] TurboPerf - Turboperformance class
   * \param[in] config - Definition of the particular problem
   */
  void LoadTurboHistoryData(std::shared_ptr<CTurbomachineryStagePerformance> TurboStagePerf, std::shared_ptr<CTurboOutput> TurboPerf, CConfig *config) override;

  /*!
   * \brief Write the kinematic and thermodynamic variables at each spanwise division
   * \param[in] TurboPerf - Turboperformance class
   * \param[in] geometry - Geometrical definiton of the problem
   * \param[in] config - Descripiton of the particular problem
   * \param[in] val_iZone - Idientifier of current zone
  */
  void WriteTurboSpanwisePerformance(std::shared_ptr<CTurboOutput> TurboPerf, CGeometry *geometry, CConfig **config,
                                       unsigned short val_iZone) override;
};
