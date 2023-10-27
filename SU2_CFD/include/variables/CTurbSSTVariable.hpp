/*!
 * \file CTurbSSTVariable.hpp
 * \brief Declaration of the variables of the SST turbulence model.
 * \author F. Palacios, T. Economon
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

#include "CTurbVariable.hpp"

/*!
 * \class CTurbSSTVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */

class CTurbSSTVariable final : public CTurbVariable {
protected:
  su2double sigma_om2;
  su2double beta_star;
  su2double prod_lim_const;
  VectorType F1;
  VectorType F2;    /*!< \brief Menter blending function for blending of k-w and k-eps. */
  VectorType CDkw;  /*!< \brief Cross-diffusion. */
  SST_ParsedOptions sstParsedOptions;

  VectorType ftilda_d;
  VectorType l_RANS;
  VectorType l_LES;
  VectorType r_dl;
  VectorType r_dt;
  VectorType r_d;
  VectorType FTrans;
  MatrixType VelocityLaplacian;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] kine - Turbulence kinetic energy (k) (initialization value).
   * \param[in] omega - Turbulent variable value (initialization value).
   * \param[in] mut - Eddy viscosity (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] constants - sst model constants
   * \param[in] config - Definition of the particular problem.
   */
  CTurbSSTVariable(su2double kine, su2double omega, su2double mut, unsigned long npoint,
                   unsigned long ndim, unsigned long nvar, const su2double* constants, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbSSTVariable() override = default;

  /*!
   * \brief Set the blending function for the blending of k-w and k-eps.
   * \param[in] val_viscosity - Value of the vicosity.
   * \param[in] val_dist - Value of the distance to the wall.
   * \param[in] val_density - Value of the density.
   */
  void SetBlendingFunc(unsigned long iPoint, su2double val_viscosity, su2double val_dist, su2double val_density, TURB_TRANS_MODEL trans_model) override;

  /*!
   * \brief Get the first blending function.
   */
  inline su2double GetF1blending(unsigned long iPoint) const override { return F1(iPoint); }

  /*!
   * \brief Get the second blending function.
   */
  inline su2double GetF2blending(unsigned long iPoint) const override { return F2(iPoint); }

  /*!
   * \brief Get the value of the cross diffusion of tke and omega.
   */
  inline su2double GetCrossDiff(unsigned long iPoint) const override { return CDkw(iPoint); }

  /*!
   * \brief Set the DES Length Scale.
   * \param[in] iPoint - Point index.
   */
  inline void SetDebug_Quantities(CConfig *config, unsigned long iPoint, su2double val_ftilda_d, su2double val_l_RANS, su2double val_l_LES, su2double val_r_d) override { 
    ftilda_d(iPoint) = val_ftilda_d;
    l_RANS(iPoint) = val_l_RANS;
    l_LES(iPoint) = val_l_LES;

    if ( config->GetKind_HybridRANSLES() == SST_DDES)
      r_d(iPoint) = val_r_d;
    else
      r_dt(iPoint) = val_r_d;
  }

  inline void SetDebug_Quantities(CConfig *config, unsigned long iPoint, su2double val_ftilda_d, su2double val_l_RANS, su2double val_l_LES, su2double val_r_dl, su2double val_r_dt) override { 
    ftilda_d(iPoint) = val_ftilda_d;
    l_RANS(iPoint) = val_l_RANS;
    l_LES(iPoint) = val_l_LES;
    r_dt(iPoint) = val_r_dt;
    r_dl(iPoint) = val_r_dl;
  }

  inline su2double Get_L_RANS(unsigned long iPoint) const override { return l_RANS(iPoint); }
  inline su2double Get_L_LES(unsigned long iPoint) const override { return l_LES(iPoint); }
  inline su2double Get_ftilda_d(unsigned long iPoint) const override { return ftilda_d(iPoint); }
  inline su2double Get_r_dt(unsigned long iPoint) const override { return r_dt(iPoint); }
  inline su2double Get_r_dl(unsigned long iPoint) const override { return r_dl(iPoint); }
  inline su2double Get_r_d(unsigned long iPoint) const override { return r_d(iPoint); }


  /*!
   * \brief Get the value of the FTrans.
   * \param[in] iPoint - Point index.
   */
  inline su2double GetFTrans(unsigned long iPoint) const override { return FTrans(iPoint); }

  /*!
   * \brief Set the value of the FTrans.
   * \param[in] iPoint - Point index.
   * \param[in] val_FTrans - Value of the FTrans variable.
   */
  inline void SetFTrans(unsigned long iPoint, su2double val_FTrans) override { 
    FTrans(iPoint) = val_FTrans;
  };

  /*!
   * \brief Get the value of the velocity laplacian.
   * \param[in] iPoint - Point index.
   */
  inline su2double* GetVelLapl(unsigned long iPoint) final { return VelocityLaplacian[iPoint]; }

  /*!
   * \brief Get the value of the velocity laplacian.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   */
  inline su2double GetVelLapl(unsigned long iPoint, unsigned short iDim) const final { return VelocityLaplacian[iPoint][iDim]; }

  /*!
   * \brief Incrementally add the velocity laplacian vector.
   * \param[in] iPoint - Point index.
   * \param[in] val_VelLapl_X - X-Component of the velocity laplacian.
   * \param[in] val_VelLapl_Y - Y-Component of the velocity laplacian.
   * \param[in] val_VelLapl_Z - Z-Component of the velocity laplacian.
   */
  inline void AddVelLapl(unsigned long iPoint, su2double val_VelLapl_X, su2double val_VelLapl_Y, su2double val_VelLapl_Z) override { 
    VelocityLaplacian(iPoint, 0) += val_VelLapl_X;
    VelocityLaplacian(iPoint, 1) += val_VelLapl_Y;
    if(nDim == 3) VelocityLaplacian(iPoint, 2) += val_VelLapl_Z;
  };

  /*!
   * \brief Set the value of the velocity laplacian.
   * \param[in] iPoint - Point index.
   * \param[in] val_VelLapl - Vector of velocity laplacian.
   */
  inline void SetVelLapl(unsigned long iPoint, su2double* val_VelLapl) override { 
    VelocityLaplacian(iPoint, 0) = val_VelLapl[0];
    VelocityLaplacian(iPoint, 1) = val_VelLapl[1];
    if(nDim == 3) VelocityLaplacian(iPoint, 2) = val_VelLapl[2];
  };

  /*!
   * \brief Set the value of the velocity laplacian.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \param[in] val_VelLapl - value of the velocity laplacian.
   */
  inline void SetVelLapl(unsigned long iPoint, unsigned short iDim, su2double val_VelLapl) override { 
    VelocityLaplacian(iPoint, iDim) = val_VelLapl;
  };


};
