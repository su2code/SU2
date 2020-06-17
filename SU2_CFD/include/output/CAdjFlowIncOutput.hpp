/*!
 * \file CAdjFlowIncOutput.hpp
 * \brief Headers of the adjoint incompressible flow output.
 * \author T. Albring
 * \version 7.0.5 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "COutput.hpp"
#include "modules/CCommonModule.hpp"
#include "modules/CConvergenceModule.hpp"
#include "modules/CResidualModule.hpp"
#include "modules/CFVMBaseModule.hpp"

class CAdjFlowIncOutputModule : public CSolverOutputModule {
  bool heat;                 /*!< \brief Boolean indicating whether have a heat problem*/
  bool weakly_coupled_heat;  /*!< \brief Boolean indicating whether have a weakly coupled heat equation*/
  unsigned short rad_model;  /*!< \brief The kind of radiation model */

public:

  CAdjFlowIncOutputModule(CConfig *config, int nDim) : CSolverOutputModule(nDim),
    heat(config->GetEnergy_Equation()), weakly_coupled_heat(config->GetWeakly_Coupled_Heat()),
  rad_model(config->GetKind_RadiationModel()){}

  void LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                       const IterationInfo& iterationInfo) override;

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields) override;

  void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) override;

  void LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                      const IterationInfo& iterationInfo, const PointInfo& pointInfo) override;

  void LoadSurfaceData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                       const IterationInfo& iterationInfo, const PointInfo& pointInfo) override;
};

/*! \class CAdjFlowIncOutput
 *  \brief Output class for incompressible flow discrete adjoint problems.
 *  \author R. Sanchez, T. Albring.
 *  \date June 5, 2018.
 */
class CAdjFlowIncOutput final: public COutput {

  using Modules = ModuleList<CCommonModule,
                             CFVMBaseModule,
                             CAdjFlowIncOutputModule>;

  using Modifiers = ModuleList<CResidualModule,
                               CConvergenceModule>;

private:


public:


  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CAdjFlowIncOutput(CConfig *config, unsigned short nDim);

  /*!
   * \brief Destructor of the class.
   */
  ~CAdjFlowIncOutput(void) override;

  /*!
   * \brief Check whether the base values for relative residuals should be initialized
   * \param[in] config - Definition of the particular problem.
   * \return <TRUE> if the residuals should be initialized.
   */
  bool SetInit_Residuals(CConfig *config) override;

  /*!
   * \brief Check whether the averaged values should be updated
   * \param[in] config - Definition of the particular problem.
   * \return <TRUE> averages should be updated.
   */
  bool SetUpdate_Averages(CConfig *config) override;

};
