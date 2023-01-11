/*!
 * \file ausmplusm.hpp
 * \brief Declaration of numerics classes for the AUSM family of schemes in NEMO - AUSM+M.
 * \author F. Morgado
 * \version 7.5.0 "Blackbird"
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

#include "../CNEMONumerics.hpp"

/*!
 * \class CUpwAUSMPLUSM_NEMO
 * \brief Class for solving an approximate Riemann AUSM+ M, Two-Temperature Model.
 * https://doi.org/10.1016/j.apm.2019.09.005 \ingroup ConvDiscr \author F. Morgado
 */
class CUpwAUSMPLUSM_NEMO : public CNEMONumerics {
 private:
  su2double* FcL = nullptr;
  su2double* FcR = nullptr;
  su2double* dmLP = nullptr;
  su2double* dmRM = nullptr;
  su2double* dpLP = nullptr;
  su2double* dpRM = nullptr;
  su2double* daL = nullptr;
  su2double* daR = nullptr;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive variables of the problem.
   * \param[in] val_nPrimVarGrad - Number of primitive gradient variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSMPLUSM_NEMO(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar,
                     unsigned short val_nPrimVarGrad, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwAUSMPLUSM_NEMO(void);

  /*!
   * \brief Compute the AUSM+M flux between two nodes i and j.
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final;
};
