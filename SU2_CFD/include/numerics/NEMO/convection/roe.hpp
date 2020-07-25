/*!
 * \file roe.hpp
 * \brief Delarations of numerics classes for Roe-type schemes in NEMO,
 *        implemented in roe.cpp.
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

#include "../../CNumerics.hpp"

/*!
 * \class CUpwRoe_NEMO
 * \brief Class for evaluating the Riemann problem using Roe's scheme for a two-temperature model.
 * \ingroup ConvDiscr
 * \author S. Copeland, W. Maier
 * \version 6.1.0 "Falcon"
 */
class CUpwRoe_NEMO : public CNumerics {
private:
    bool implicit, ionization;
    su2double *Diff_U;
    su2double *RoeU, *RoeV;
    vector<su2double> roe_eves;
    su2double *ProjFlux_i, *ProjFlux_j;
    su2double *Lambda, *Epsilon;
    su2double **P_Tensor, **invP_Tensor;
    su2double RoeSoundSpeed;
    su2double ProjVelocity, ProjVelocity_i, ProjVelocity_j;
    su2double Proj_ModJac_Tensor_ij, R;
    su2double *RoedPdU;
    unsigned short nPrimVar, nPrimVarGrad;

    su2double* Flux = nullptr;        /*!< \brief The flux / residual across the edge. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwRoe_NEMO(unsigned short val_nDim, unsigned short val_nVar,
             unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
             CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoe_NEMO(void);

  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final;

};
