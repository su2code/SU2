/*!
 * \file CAdjFlowCompOutput.hpp
 * \brief Headers of the adjoint compressible flow output.
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

class CAdjFlowCompOutputModule : public CSolverOutputModule {

public:

  CAdjFlowCompOutputModule(CConfig *config, int nDim) : CSolverOutputModule(nDim){}

  void LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                       const IterationInfo& iterationInfo) override;

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields) override;

  void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) override;

  void LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                      const IterationInfo& iterationInfo, const PointInfo& pointInfo) override;

  void LoadSurfaceData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                       const IterationInfo& iterationInfo, const PointInfo& pointInfo) override;
};

/*! \class CAdjFlowCompOutput
 *  \brief Output class for compressible flow adjoint problems.
 *  \author R. Sanchez, T. Albring.
 *  \date June 5, 2018.
 */
class CAdjFlowCompOutput final: public COutput {

  using Modules = ModuleList<CCommonModule,
                             CFVMBaseModule,
                             CAdjFlowCompOutputModule>;

  using Modifiers = ModuleList<CResidualModule,
                               CConvergenceModule>;
private:

public:

  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CAdjFlowCompOutput(CConfig *config, unsigned short nDim);

  /*!
   * \brief Destructor of the class.
   */
  ~CAdjFlowCompOutput(void) override;

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
