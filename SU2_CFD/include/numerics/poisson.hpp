/*!
 * \file heat.hpp
 * \brief Delarations of numerics classes for heat transfer problems.
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
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

#include "CNumerics.hpp"

/*!
 * \class CAvgGrad_Poisson
 * \brief Class for computing viscous term using average of gradients without correction (Poisson equation).
 * \ingroup ViscDiscr
 */
class CAvgGrad_Poisson : public CNumerics {
private:

  su2double *Edge_Vector;
  bool implicit,direct;
  su2double **Mean_GradPoissonVar;
  su2double *Mom_Coeff_i,*Mom_Coeff_j;
  su2double *Proj_Mean_GradPoissonVar_Normal, *Proj_Mean_GradPoissonVar_Corrected;
  su2double dist_ij_2, proj_vector_ij, Poisson_Coeff_Mean ;
  unsigned short iVar, iDim;
  
  su2double *Flux = nullptr;
  su2double **Jacobian_i = nullptr;
  su2double **Jacobian_j = nullptr;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Poisson(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Poisson(void);

  /*!
   * \brief Compute the viscous heat residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;
  
  
  inline void SetInvMomCoeff(su2double *val_Mom_Coeff_i, su2double *val_Mom_Coeff_j) { 
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Mom_Coeff_i[iDim] = val_Mom_Coeff_i[iDim];
      Mom_Coeff_j[iDim] = val_Mom_Coeff_j[iDim];
    }
  }
};

/*!
 * \class CAvgGradCorrected_Poisson
 * \brief Class for computing viscous term using average of gradients with correction (Poisson equation).
 * \ingroup ViscDiscr
 */

class CAvgGradCorrected_Poisson : public CNumerics {
private:

  su2double *Edge_Vector;
  bool implicit,direct;
  su2double **Mean_GradPoissonVar;
  su2double *Mom_Coeff_i,*Mom_Coeff_j;
  su2double *Proj_Mean_GradPoissonVar_Kappa, *Proj_Mean_GradPoissonVar_Edge, *Proj_Mean_GradPoissonVar_Corrected;
  su2double dist_ij_2, proj_vector_ij, Poisson_Coeff_Mean;
  unsigned short iVar, iDim;
  
  su2double *Flux = nullptr;
  su2double **Jacobian_i = nullptr;
  su2double **Jacobian_j = nullptr;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_Poisson(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_Poisson(void);

  /*!
   * \brief Compute the viscous residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;
  
  
  inline void SetInvMomCoeff(su2double *val_Mom_Coeff_i, su2double *val_Mom_Coeff_j) { 
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Mom_Coeff_i[iDim] = val_Mom_Coeff_i[iDim];
      Mom_Coeff_j[iDim] = val_Mom_Coeff_j[iDim];
    }
  }
};

/*!
 * \class CAvgGrad_Poisson
 * \brief Class for computing viscous term using average of gradients without correction (Poisson equation).
 * \ingroup ViscDiscr
 */
class CPressure_Poisson : public CNumerics {
private:

  su2double *Edge_Vector;
  bool implicit,direct;
  su2double **Mean_GradPoissonVar;
  su2double *Proj_Mean_GradPoissonVar_Kappa, *Proj_Mean_GradPoissonVar_Edge, *Proj_Mean_GradPoissonVar_Corrected;
  su2double dist_ij_2, proj_vector_ij, Poisson_Coeff_Mean;
  su2double *Mom_Coeff_i,*Mom_Coeff_j;
  unsigned short iVar, iDim;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPressure_Poisson(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CPressure_Poisson(void);

  /*!
   * \brief Compute the viscous heat residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
  
  
  inline void SetInvMomCoeff(su2double *val_Mom_Coeff_i, su2double *val_Mom_Coeff_j)  { 
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Mom_Coeff_i[iDim] = val_Mom_Coeff_i[iDim];
      Mom_Coeff_j[iDim] = val_Mom_Coeff_j[iDim];
    }
  }
};

/*!
 * \class CSource_Poisson
 * \brief Class for source term of the Poisson equation.
 * \ingroup SourceDiscr
 */
class CSource_PoissonFVM : public CNumerics {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config -  Name of the input config file
   *
   */
  CSource_PoissonFVM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);


  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSource_PoissonFVM(void);
};

