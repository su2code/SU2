/*!
 * \file CPoissonVariable.hpp
 * \brief Class for defining the variables of the finite-volume heat equation solver.
 * \author F. Palacios, T. Economon
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "CScalarVariable.hpp"

/*!
 * \class CPoissonVariable
 * \brief Class for defining the variables of the finite-volume poisson equation solver.
 * \author O. Burghardt
 * \version 8.5.0 "Harrier"
 */
class CPoissonVariable final : public CScalarVariable {
protected:
  VectorType MomCoeff;            /*!< \brief Momentum coefficients vol/A_p used as the diffusion coefficients in the poisson solver.  */
  MatrixType MomentumCorrection;  /*!< \brief (rho*u)' in the context of: (rho*u)** = (rho*u)* + (rho*u)'. */
  MatrixType HbyACorrection;      /*!< \brief H(rhou')/A = (sum_nb A_nb (rhou)'_nb) / A; used by the second pressure correction in the PISO algorithm. */
public:
  static constexpr size_t MAXNVAR = 1; /*!< \brief Max number of variables, for static arrays. */

  /*!
   * \brief Constructor of the class.
   * \param[in] value - Values of the poisson solution (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPoissonVariable(su2double value, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Get the momentum coefficient of the point.
   * \return Value of the momentum coefficient of the point.
   */
  inline su2double GetMomCoeff(unsigned long iPoint) final { return MomCoeff(iPoint);}
    
  /*!
   * \brief Set the momentum coefficient of the point.
   */
  inline void SetMomCoeff(unsigned long iPoint, su2double val_Mom_Coeff) final { MomCoeff(iPoint) = val_Mom_Coeff; }

  /*!
   * \brief Set H(u')/A for the point
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \param[in] val_HbyA - HbyA correction.
   */
  inline void SetHbyACorrection(unsigned long iPoint, unsigned short iDim, su2double val_HbyA) final { HbyACorrection(iPoint, iDim) = val_HbyA; }
  
  /*!
   * \brief Get H(u')/A for the point
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \return The H(u')/A for the point.
   */
  inline su2double GetHbyACorrection(unsigned long iPoint, unsigned short iDim) final { return HbyACorrection(iPoint, iDim); }

  /*!
   * \brief Set (rho*u)' for the point
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \param[in] val_mom - Momentum correction (rho*u)' value.
   */
  inline void SetVelocityCorrection(unsigned long iPoint, unsigned short iDim, su2double val_mom) final { MomentumCorrection(iPoint, iDim) = val_mom; }
  
  /*!
   * \brief Get (rho*u)' for the point
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \return The (rho*u)' for the point.
   */
  inline su2double GetVelocityCorrection(unsigned long iPoint, unsigned short iDim) final { return MomentumCorrection(iPoint, iDim); }

};
