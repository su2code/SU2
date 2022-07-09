/*!
 * \file ausmplusup2.hpp
 * \brief Declaration of numerics classes for the AUSM family of schemes in NEMO - AUSM+UP2.
 * \author W. Maier, A. Sachedeva, C. Garbacz.
 * \version 7.3.1 "Blackbird"
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
 * \class CUpwAUSMPLUSUP2_NEMO
 * \brief Class for solving an approximate Riemann AUSM+ -up2, Two-Temperature Model. https://doi.org/10.1016/j.jcp.2013.02.046
 * \ingroup ConvDiscr
 * \author W. Maier, A. Sachedeva, C. Garbacz
 */
class CUpwAUSMPLUSUP2_NEMO : public CNEMONumerics {
private:
  su2double *FcL, *FcR, *FcLR;
  su2double *dmLP, *dmRM, *dpLP, *dpRM;
  su2double *daL, *daR;
  su2double ProjVel_i, ProjVel_j;
  su2double mL, mR, mLP, mRM, mF, pLP, pRM, pFi, pF, Phi;
  su2double CstarL, CstarR, ChatL, ChatR, aF, rhoF, MFsq, Mrefsq, Mp, fa;
  su2double Kp, sigma, alpha, beta, param1, mfP, mfM;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSMPLUSUP2_NEMO(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,  CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwAUSMPLUSUP2_NEMO(void);

  /*!
   * \brief Compute the AUSM+ -up flux between two nodes i and j.
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final;
};
