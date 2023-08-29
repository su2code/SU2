/*!
 * \file template.cpp
 * \brief Empty implementation of numerics templates, see .hpp file.
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

#include "../../include/numerics/template.hpp"

CConvective_Template::CConvective_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                      CNumerics(val_nDim, val_nVar, config) {

}

CConvective_Template::~CConvective_Template() = default;

void CConvective_Template::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                           su2double **val_Jacobian_j, CConfig *config) {

}

CViscous_Template::CViscous_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                   CNumerics(val_nDim, val_nVar, config) {

}

CViscous_Template::~CViscous_Template() = default;

void CViscous_Template::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                        su2double **val_Jacobian_j, CConfig *config) {

}


CSource_Template::CSource_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                  CNumerics(val_nDim, val_nVar, config) {

}

CSource_Template::~CSource_Template() = default;

void CSource_Template::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {

}
