/*!
 * \file NEMO_sources.hpp
 * \brief Declarations of numerics classes for source-term integration.
 * \author C. Garbacz, W. Maier, S. Copeland.
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

#include "CNEMONumerics.hpp"

/*!
 * \class CSource_NEMO
 * \brief Class for two-temperature model source terms.
 * \ingroup SourceDiscr
 * \author C. Garbacz, W. Maier, S. Copeland.
 * \version 8.0.0 "Harrier"
 */
class CSource_NEMO : public CNEMONumerics {
private:

  su2double *Y, **dYdr;                  // Mass fraction

  su2double*  residual = nullptr;        /*!< \brief The source residual. */
  su2double** jacobian = nullptr;
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSource_NEMO(unsigned short val_nDim,
               unsigned short val_nVar,
               unsigned short val_nPrimVar,
               unsigned short val_nPrimVarGrad,
               CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSource_NEMO(void);

  /*!
   * \brief Source residual of the chemistry.
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeChemistry(const CConfig* config) final;

  /*!
  * \brief Calculates constants used for Keq correlation.
  * \param[out] A - Pointer to coefficient array.
  * \param[in] val_reaction - Reaction number indicator.
  * \param[in] config - Definition of the particular problem.
  */
  void ComputeKeqConstants(su2double *A, unsigned short val_reaction, CConfig *config);

  /*!
   * \brief Residual of the translational to vibrational energy.
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeVibRelaxation(const CConfig* config) final;

  /*!
   * \brief Residual of axissymetric source term.
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeAxisymmetric(const CConfig* config) final;
};

