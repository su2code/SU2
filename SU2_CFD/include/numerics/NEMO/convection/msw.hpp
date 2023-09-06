/*!
 * \file msw.hpp
 * \brief Declaration of numerics classes for modified Steger-Warming scheme.
 * \author ADL Stanford, S.R. Copeland, W. Maier, C. Garbacz
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

#include "../CNEMONumerics.hpp"

/*!
 * \class CUpwMSW_NEMO
 * \brief Class for solving a flux-vector splitting method by Steger & Warming, modified version.
 * \ingroup ConvDiscr
 * \author ADL Stanford, S.R. Copeland, W. Maier, C. Garbacz
 * \version 8.0.0 "Harrier"
 */
class CUpwMSW_NEMO : public CNEMONumerics {
private:
    su2double *Diff_U;
    su2double *ust_i, *ust_j;
    su2double *Fc_i, *Fc_j;
    su2double *Lambda_i, *Lambda_j;
    su2double *rhosst_i, *rhosst_j;
    su2double ProjVel_i, ProjVel_j;
    su2double *Ust_i, *Ust_j, *Vst_i, *Vst_j;
    su2double *dPdUst_i, *dPdUst_j;
    su2double **P_Tensor, **invP_Tensor;

public:

    /*!
     * \brief Constructor of the class.
     * \param[in] val_nDim - Number of dimensions of the problem.
     * \param[in] val_nVar - Number of variables of the problem.
     * \param[in] val_nPrimVar - Number of primitive variables of the problem.
     * \param[in] val_nPrimVarGrad - Number primitive grad. variables of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    CUpwMSW_NEMO(unsigned short val_nDim, unsigned short val_nVar,
               unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
               CConfig *config);

    /*!
     * \brief Destructor of the class.
     */
    ~CUpwMSW_NEMO(void);

    /*!
     * \brief Compute the Roe's flux between two nodes i and j.
     * \param[in] config - Definition of the particular problem.
     */
    ResidualType<> ComputeResidual(const CConfig* config) final;

};
