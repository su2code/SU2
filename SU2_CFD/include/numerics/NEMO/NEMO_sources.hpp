/*!
 * \file NEMO_sources.hpp
 * \brief Delarations of numerics classes for source-term integration. The implementation is in NEMO_sources.cpp.
 * \author F. Palacios, T. Economon
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

#include "../CNumerics.hpp"
#include "../../variables/CNEMOEulerVariable.hpp"
/*!
 * \class CSource_NEMO
 * \brief Class for two-temperature model source terms.
 * \ingroup SourceDiscr
 * \author S. Copeland, W. Maier
 * \version 6.1.0 "Falcon"
 */
class CSource_NEMO : public CNumerics {
private:
  bool   implicit, ionization;
  unsigned short nSpecies, nVar, nPrimVar, nPrimVarGrad;
  int    *alphak, *betak;
  su2double *X; // Mole fraction
  su2double *Y, **dYdr; // Mass fraction
  su2double **RxnConstantTable;
  su2double *estar, *evib;
  su2double **tau_sr, *tauP, *tauMW, *taus;
  su2double *dkf, *dkb, *dRfok, *dRbok, *A;
  su2double *Cvvs, *Cvvsst;
  su2double *Cves;
  CNEMOEulerVariable *variable;

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
  * \param[out] val_residual - Pointer to the total residual.
  * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
  * \param[in] config - Definition of the particular problem.
  */
 void ComputeChemistry(su2double *val_residual, su2double *val_source, su2double **val_Jacobian_i, CConfig *config);

 /*!
  * \brief Calculates constants used for Keq correlation.
  * \param[out] A - Pointer to coefficient array.
  * \param[in] val_reaction - Reaction number indicator.
  * \param[in] config - Definition of the particular problem.
  */
 void GetKeqConstants(su2double *A, unsigned short val_reaction, CConfig *config);

 /*!
  * \brief Residual of the rotational frame source term.
  * \param[out] val_residual - Pointer to the total residual.
  * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
  * \param[in] config - Definition of the particular problem.
  */
 void ComputeVibRelaxation(su2double *val_residual,su2double *val_source, su2double **val_Jacobian_i, CConfig *config);

 void ComputeAxisymmetric(su2double *val_residual, su2double *val_source,su2double **val_Jacobian, CConfig *config);
};

