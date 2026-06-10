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
  MatrixType MomCoeff; /*!< \brief Eddy viscosity. */
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
   * \brief Get the temperature of the point.
   * \return Value of the temperature of the point.
   */
  // inline su2double GetTemperature(unsigned long iPoint) const final { return Solution(iPoint, 0); }

  //TODO: place comments
  inline su2double GetMomCoeff(unsigned long iPoint, unsigned short val_Var) final { return MomCoeff(iPoint,val_Var);}

  inline void SetMomCoeff(unsigned long iPoint, const su2double *val_Mom_Coeff) { 
	  for (unsigned short iDim = 0; iDim < nDim; iDim++)
	        MomCoeff(iPoint,iDim) = val_Mom_Coeff[iDim]; 
  }
    
  inline void SetMomCoeff(unsigned long iPoint, unsigned short val_Var, su2double val_Mom_Coeff) final { MomCoeff(iPoint,val_Var) = val_Mom_Coeff; }

  inline void AddMomCoeff(unsigned long iPoint, su2double val_coeff, unsigned short val_Var) { MomCoeff(iPoint,val_Var) += val_coeff;}
    
  // inline su2double Get_Mom_Coeff_nb(unsigned long iPoint, unsigned short val_Var) { return Mom_Coeff_nb(iPoint,val_Var);}

  // inline void Set_Mom_Coeff_nb(unsigned long iPoint, su2double *val_Mom_Coeff) { 
	//   for (unsigned short iDim = 0; iDim < nDim; iDim++)
	//         Mom_Coeff_nb(iPoint,iDim) = val_Mom_Coeff[iDim]; 
  // }
  
  // inline void Set_Mom_Coeff_nb(unsigned long iPoint, unsigned short val_Var, su2double val_Mom_Coeff) { Mom_Coeff_nb(iPoint,val_Var) = val_Mom_Coeff; }
  
  // inline void Set_Mom_Coeff_nbZero(unsigned long iPoint) {
	//   for (unsigned short iDim = 0; iDim < nDim; iDim++)
	//         Mom_Coeff_nb(iPoint,iDim) = 0.0;
  // }
  
  // inline void Add_Mom_Coeff_nb(unsigned long iPoint, su2double val_coeff_nb, unsigned short val_Var) { Mom_Coeff_nb(iPoint,val_Var) += val_coeff_nb;}

  

};
