/*!
 * \file species_sources.cpp
 * \brief Implementation of numerics classes for integration of
 *        species transport source-terms.
 * \author T. Kattmann
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

#include "../../../include/numerics/species/species_sources.hpp"

#include "../../../include/numerics/CNumerics.hpp"

#include "../../../include/variables/CEulerVariable.hpp"
#include "../../../include/variables/CIncEulerVariable.hpp"
#include "../../../include/variables/CNEMOEulerVariable.hpp"

CSourceBase_Species::CSourceBase_Species(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
    : CNumerics(val_nDim, val_nVar, config) {
  residual = new su2double[nVar]();
  jacobian = new su2double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    jacobian[iVar] = new su2double[nVar]();
  }
  bounded_scalar = config->GetBounded_Species();
}

CSourceBase_Species::~CSourceBase_Species() {
  delete[] residual;
  if (jacobian) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete[] jacobian[iVar];
    }
    delete[] jacobian;
  }
}

template <class T>
CSourceAxisymmetric_Species<T>::CSourceAxisymmetric_Species(unsigned short val_nDim, unsigned short val_nVar,
                                                         const CConfig* config)
    : CSourceBase_Species(val_nDim, val_nVar, config),
      idx(val_nDim, config->GetnSpecies()),
    implicit(config->GetKind_TimeIntScheme_Species() == EULER_IMPLICIT),
    viscous(config->GetViscous()),
    turbulence(config->GetKind_Turb_Model() != TURB_MODEL::NONE),
    incompressible(config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE),
    Sc_t(config->GetSchmidt_Number_Turbulent()) {
}

template <class T>
CNumerics::ResidualType<> CSourceAxisymmetric_Species<T>::ComputeResidual(const CConfig* config) {
  /*--- Preaccumulation ---*/
  AD::StartPreacc();
  AD::SetPreaccIn(ScalarVar_i, nVar);
  AD::SetPreaccIn(Volume);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+6);
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);
  }

  /*--- Initialization. ---*/
  for (auto iVar = 0u; iVar < nVar; iVar++) {
    residual[iVar] = 0.0;
    for (auto jVar = 0; jVar < nVar; jVar++) {
      jacobian[iVar][jVar] = 0.0;
    }
  }

  /*--- Contribution due to 2D axisymmetric formulation ---*/
  if (Coord_i[1] > EPS) {

    AD::SetPreaccIn(Coord_i[1]);
    AD::SetPreaccIn(Diffusion_Coeff_i, nVar);
    AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);

    const su2double yinv = 1.0 / Coord_i[1];

    const su2double Density_i = V_i[idx.Density()];

    su2double Velocity_i[3];
    for (auto iDim = 0u; iDim < nDim; iDim++)
      Velocity_i[iDim] = V_i[idx.Velocity() + iDim];

    /*--- Inviscid component of the source term. When div(v)=0, this term does not occur ---*/

    if (!bounded_scalar) {
      for (auto iVar = 0u; iVar < nVar; iVar++)
        residual[iVar] -= yinv * Volume * Density_i * ScalarVar_i[iVar] * Velocity_i[1];

      if (implicit) {
        for (auto iVar = 0u; iVar < nVar; iVar++) {
          jacobian[iVar][iVar] -= yinv * Volume * Velocity_i[1];
        }
      }
    }

    /*--- Add the viscous terms if necessary. ---*/

    if (config->GetViscous()) {
      su2double Mass_Diffusivity_Tur = 0.0;
      if (turbulence)
        Mass_Diffusivity_Tur = V_i[idx.EddyViscosity()] / Sc_t;

      for (auto iVar=0u; iVar < nVar; iVar++){
        residual[iVar] += yinv * Volume * (Density_i * Diffusion_Coeff_i[iVar] + Mass_Diffusivity_Tur) * ScalarVar_Grad_i[iVar][1];
      }
    }

  }

  AD::SetPreaccOut(residual, nVar);
  AD::EndPreacc();

  return ResidualType<>(residual, jacobian, nullptr);
}

/*--- Explicit instantiations until we don't move this to the hpp. ---*/
template class CSourceAxisymmetric_Species<CEulerVariable::CIndices<unsigned short> >;
template class CSourceAxisymmetric_Species<CIncEulerVariable::CIndices<unsigned short> >;
template class CSourceAxisymmetric_Species<CNEMOEulerVariable::CIndices<unsigned short> >;


CSourcePieceWise_transportedScalar_general::CSourcePieceWise_transportedScalar_general(unsigned short val_nDim,
                                                   unsigned short val_nVar,
                                                   const CConfig* config) :
                          CNumerics(val_nDim, val_nVar, config) {

  incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
  axisymmetric = config->GetAxisymmetric();
  viscous = config->GetViscous();
  implicit = (config->GetKind_TimeIntScheme_Species() == EULER_IMPLICIT);
  flamelet  = (config->GetKind_Species_Model() == SPECIES_MODEL::FLAMELET);
  inc_rans = (config->GetKind_Solver() == MAIN_SOLVER::INC_RANS) || (config->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_INC_RANS);

  Sc_t = config->GetSchmidt_Number_Turbulent();

  Residual       = new su2double [nVar];
  scalar_sources = new su2double [nVar];
  Jacobian_i     = new su2double* [nVar];

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar] ();
  }
}

CSourcePieceWise_transportedScalar_general::~CSourcePieceWise_transportedScalar_general(void){
  delete [] Residual;
  delete [] scalar_sources;
   if (Jacobian_i != nullptr) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
    }
    delete [] Jacobian_i;
  }
}

CNumerics::ResidualType<> CSourcePieceWise_transportedScalar_general::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  AD::SetPreaccIn(ScalarVar_i, nVar);
  AD::SetPreaccIn(scalar_sources, nVar);
  AD::SetPreaccIn(Volume); 

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+6);

  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);

  }

  /*--- Implicit part for production term (to do). ---*/
  for (auto i_var = 0; i_var < nVar; i_var++) {
    Residual[i_var] = scalar_sources[i_var] * Volume;
    for (auto j_var = 0; j_var < nVar; j_var++) {
      Jacobian_i[i_var][j_var] = 0.0;
    }
  }

   /*--- Contribution due to 2D axisymmetric formulation ---*/
   
   if (axisymmetric) ResidualAxisymmetric();

   /*--- Implicit part ---*/
   
  AD::SetPreaccOut(Residual, nVar);
  AD::EndPreacc();

  return ResidualType<>(Residual, Jacobian_i, nullptr);

}
