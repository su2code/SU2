/*!
 * \file CNumerics.cpp
 * \brief Implementation of the base for all numerics classes.
 *        Contains methods for common tasks, e.g. compute flux
 *        Jacobians.
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


#include "../../include/numerics/CNumerics.hpp"
#include "../../include/fluid/CFluidModel.hpp"

CNumerics::CNumerics() {

  Proj_Flux_Tensor  = nullptr;

  tau = nullptr;

  nemo = false;

}

CNumerics::CNumerics(unsigned short val_nDim, unsigned short val_nVar,
                     const CConfig* config) {

  unsigned short iDim;

  nDim = val_nDim;
  nVar = val_nVar;

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  Prandtl_Lam = config->GetPrandtl_Lam();
  Prandtl_Turb = config->GetPrandtl_Turb();
  Gas_Constant = config->GetGas_ConstantND();

  tau = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    tau[iDim] = new su2double [nDim] ();

  Proj_Flux_Tensor = new su2double [nVar] ();

  turb_ke_i = 0.0;
  turb_ke_j = 0.0;

  Dissipation_ij = 1.0;

  /* --- Initializing variables for the UQ methodology --- */
  sstParsedOptions = config->GetSSTParsedOptions();
  Eig_Val_Comp = config->GetEig_Val_Comp();
  uq_delta_b = config->GetUQ_Delta_B();
  uq_urlx = config->GetUQ_URLX();
  uq_permute = config->GetUQ_Permute();

}

CNumerics::~CNumerics() {

  // visc
  delete [] Proj_Flux_Tensor;

  if (tau) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      delete [] tau[iDim];
    delete [] tau;
  }
}

void CNumerics::GetInviscidProjFlux(const su2double *val_density,
                                    const su2double *val_velocity,
                                    const su2double *val_pressure,
                                    const su2double *val_enthalpy,
                                    const su2double *val_normal,
                                    su2double *val_Proj_Flux) const {

  su2double rhou, rhov, rhow;

  if (nDim == 2) {

    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];

    val_Proj_Flux[0] = rhou*val_normal[0];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
    val_Proj_Flux[3] = rhou*(*val_enthalpy)*val_normal[0];

    val_Proj_Flux[0] += rhov*val_normal[1];
    val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
    val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] += rhov*(*val_enthalpy)*val_normal[1];

  }
  else {

    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];
    rhow = (*val_density)*val_velocity[2];

    val_Proj_Flux[0] = rhou*val_normal[0];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
    val_Proj_Flux[3] = rhou*val_velocity[2]*val_normal[0];
    val_Proj_Flux[4] = rhou*(*val_enthalpy)*val_normal[0];

    val_Proj_Flux[0] += rhov*val_normal[1];
    val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
    val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] += rhov*val_velocity[2]*val_normal[1];
    val_Proj_Flux[4] += rhov*(*val_enthalpy)*val_normal[1];

    val_Proj_Flux[0] += rhow*val_normal[2];
    val_Proj_Flux[1] += rhow*val_velocity[0]*val_normal[2];
    val_Proj_Flux[2] += rhow*val_velocity[1]*val_normal[2];
    val_Proj_Flux[3] += (rhow*val_velocity[2]+(*val_pressure))*val_normal[2];
    val_Proj_Flux[4] += rhow*(*val_enthalpy)*val_normal[2];

  }

}

void CNumerics::GetInviscidIncProjFlux(const su2double *val_density,
                                       const su2double *val_velocity,
                                       const su2double *val_pressure,
                                       const su2double *val_betainc2,
                                       const su2double *val_enthalpy,
                                       const su2double *val_normal,
                                       su2double *val_Proj_Flux) const {
  su2double rhou, rhov, rhow;

  if (nDim == 2) {
    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];

    val_Proj_Flux[0] = rhou*val_normal[0] + rhov*val_normal[1];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0] + rhou*val_velocity[1]*val_normal[1];
    val_Proj_Flux[2] = rhov*val_velocity[0]*val_normal[0] + (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] = (rhou*val_normal[0] + rhov*val_normal[1])*(*val_enthalpy);
  }
  else {
    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];
    rhow = (*val_density)*val_velocity[2];

    val_Proj_Flux[0] = rhou*val_normal[0] + rhov*val_normal[1] + rhow*val_normal[2];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0] + rhou*val_velocity[1]*val_normal[1] + rhou*val_velocity[2]*val_normal[2];
    val_Proj_Flux[2] = rhov*val_velocity[0]*val_normal[0] + (rhov*val_velocity[1]+(*val_pressure))*val_normal[1] + rhov*val_velocity[2]*val_normal[2];
    val_Proj_Flux[3] = rhow*val_velocity[0]*val_normal[0] + rhow*val_velocity[1]*val_normal[1] + (rhow*val_velocity[2]+(*val_pressure))*val_normal[2];
    val_Proj_Flux[4] = (rhou*val_normal[0] + rhov*val_normal[1] + rhow*val_normal[2])*(*val_enthalpy);
  }

}

void CNumerics::GetInviscidProjJac(const su2double *val_velocity, const su2double *val_energy,
                                   const su2double *val_normal, su2double val_scale,
                                   su2double **val_Proj_Jac_Tensor) const {
  const bool wasActive = AD::BeginPassive();
  unsigned short iDim, jDim;
  su2double sqvel, proj_vel, phi, a1, a2;

  sqvel = 0.0; proj_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    sqvel    += val_velocity[iDim]*val_velocity[iDim];
    proj_vel += val_velocity[iDim]*val_normal[iDim];
  }

  phi = 0.5*Gamma_Minus_One*sqvel;
  a1 = Gamma*(*val_energy)-phi;
  a2 = Gamma-1.0;

  val_Proj_Jac_Tensor[0][0] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[0][iDim+1] = val_scale*val_normal[iDim];
  val_Proj_Jac_Tensor[0][nDim+1] = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    val_Proj_Jac_Tensor[iDim+1][0] = val_scale*(val_normal[iDim]*phi - val_velocity[iDim]*proj_vel);
    for (jDim = 0; jDim < nDim; jDim++)
      val_Proj_Jac_Tensor[iDim+1][jDim+1] = val_scale*(val_normal[jDim]*val_velocity[iDim]-a2*val_normal[iDim]*val_velocity[jDim]);
    val_Proj_Jac_Tensor[iDim+1][iDim+1] += val_scale*proj_vel;
    val_Proj_Jac_Tensor[iDim+1][nDim+1] = val_scale*a2*val_normal[iDim];
  }

  val_Proj_Jac_Tensor[nDim+1][0] = val_scale*proj_vel*(phi-a1);
  for (iDim = 0; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[nDim+1][iDim+1] = val_scale*(val_normal[iDim]*a1-a2*val_velocity[iDim]*proj_vel);
  val_Proj_Jac_Tensor[nDim+1][nDim+1] = val_scale*Gamma*proj_vel;
  AD::EndPassive(wasActive);
}


void CNumerics::GetInviscidProjJac(const su2double *val_velocity, const su2double *val_enthalpy,
                                   const su2double *val_chi, const su2double *val_kappa,
                                   const su2double *val_normal, su2double val_scale,
                                   su2double **val_Proj_Jac_Tensor) const {
  const bool wasActive = AD::BeginPassive();
  unsigned short iDim, jDim;
  su2double sqvel, proj_vel, phi, a1, a2;

  sqvel = 0.0; proj_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    sqvel += val_velocity[iDim]*val_velocity[iDim];
    proj_vel += val_velocity[iDim]*val_normal[iDim];
  }

  phi = *val_chi + 0.5*sqvel*(*val_kappa);
  a1 = *val_enthalpy;
  a2 = *val_kappa;

  val_Proj_Jac_Tensor[0][0] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[0][iDim+1] = val_scale*val_normal[iDim];
  val_Proj_Jac_Tensor[0][nDim+1] = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    val_Proj_Jac_Tensor[iDim+1][0] = val_scale*(val_normal[iDim]*phi - val_velocity[iDim]*proj_vel);
    for (jDim = 0; jDim < nDim; jDim++)
      val_Proj_Jac_Tensor[iDim+1][jDim+1] = val_scale*(val_normal[jDim]*val_velocity[iDim]-a2*val_normal[iDim]*val_velocity[jDim]);
    val_Proj_Jac_Tensor[iDim+1][iDim+1] += val_scale*proj_vel;
    val_Proj_Jac_Tensor[iDim+1][nDim+1] = val_scale*a2*val_normal[iDim];
  }

  val_Proj_Jac_Tensor[nDim+1][0] = val_scale*proj_vel*(phi-a1);
  for (iDim = 0; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[nDim+1][iDim+1] = val_scale*(val_normal[iDim]*a1-a2*val_velocity[iDim]*proj_vel);
  val_Proj_Jac_Tensor[nDim+1][nDim+1] = val_scale*(a2+1)*proj_vel;
  AD::EndPassive(wasActive);
}

void CNumerics::GetInviscidIncProjJac(const su2double *val_density, const su2double *val_velocity,
                                      const su2double *val_betainc2, const su2double *val_cp,
                                      const su2double *val_temperature, const su2double *val_dRhodT,
                                      const su2double *val_normal, su2double val_scale,
                                      su2double **val_Proj_Jac_Tensor) const {
  const bool wasActive = AD::BeginPassive();
  unsigned short iDim;
  su2double proj_vel;

  proj_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    proj_vel += val_velocity[iDim]*val_normal[iDim];

  if (nDim == 2) {

    val_Proj_Jac_Tensor[0][0] = val_scale*(proj_vel/(*val_betainc2));
    val_Proj_Jac_Tensor[0][1] = val_scale*(val_normal[0]*(*val_density));
    val_Proj_Jac_Tensor[0][2] = val_scale*(val_normal[1]*(*val_density));
    val_Proj_Jac_Tensor[0][3] = val_scale*((*val_dRhodT)*proj_vel);

    val_Proj_Jac_Tensor[1][0] = val_scale*(val_normal[0] + val_velocity[0]*proj_vel/(*val_betainc2));
    val_Proj_Jac_Tensor[1][1] = val_scale*((*val_density)*(val_normal[0]*val_velocity[0] + proj_vel));
    val_Proj_Jac_Tensor[1][2] = val_scale*(val_normal[1]*(*val_density)*val_velocity[0]);
    val_Proj_Jac_Tensor[1][3] = val_scale*((*val_dRhodT)*val_velocity[0]*proj_vel);

    val_Proj_Jac_Tensor[2][0] = val_scale*(val_normal[1] + val_velocity[1]*proj_vel/(*val_betainc2));
    val_Proj_Jac_Tensor[2][1] = val_scale*(val_normal[0]*(*val_density)*val_velocity[1]);
    val_Proj_Jac_Tensor[2][2] = val_scale*((*val_density)*(proj_vel + val_normal[1]*val_velocity[1]));
    val_Proj_Jac_Tensor[2][3] = val_scale*((*val_dRhodT)*val_velocity[1]*proj_vel);

    val_Proj_Jac_Tensor[3][0] = val_scale*((*val_cp)*(*val_temperature)*proj_vel/(*val_betainc2));
    val_Proj_Jac_Tensor[3][1] = val_scale*((*val_cp)*(*val_temperature)*val_normal[0]*(*val_density));
    val_Proj_Jac_Tensor[3][2] = val_scale*((*val_cp)*(*val_temperature)*val_normal[1]*(*val_density));
    val_Proj_Jac_Tensor[3][3] = val_scale*((*val_cp)*((*val_temperature)*(*val_dRhodT) + (*val_density))*proj_vel);

  } else {

    val_Proj_Jac_Tensor[0][0] = val_scale*(proj_vel/(*val_betainc2));
    val_Proj_Jac_Tensor[0][1] = val_scale*(val_normal[0]*(*val_density));
    val_Proj_Jac_Tensor[0][2] = val_scale*(val_normal[1]*(*val_density));
    val_Proj_Jac_Tensor[0][3] = val_scale*(val_normal[2]*(*val_density));
    val_Proj_Jac_Tensor[0][4] = val_scale*((*val_dRhodT)*proj_vel);

    val_Proj_Jac_Tensor[1][0] = val_scale*(val_normal[0] + val_velocity[0]*proj_vel/(*val_betainc2));
    val_Proj_Jac_Tensor[1][1] = val_scale*((*val_density)*(val_normal[0]*val_velocity[0] + proj_vel));
    val_Proj_Jac_Tensor[1][2] = val_scale*(val_normal[1]*(*val_density)*val_velocity[0]);
    val_Proj_Jac_Tensor[1][3] = val_scale*(val_normal[2]*(*val_density)*val_velocity[0]);
    val_Proj_Jac_Tensor[1][4] = val_scale*((*val_dRhodT)*val_velocity[0]*proj_vel);

    val_Proj_Jac_Tensor[2][0] = val_scale*(val_normal[1] + val_velocity[1]*proj_vel/(*val_betainc2));
    val_Proj_Jac_Tensor[2][1] = val_scale*(val_normal[0]*(*val_density)*val_velocity[1]);
    val_Proj_Jac_Tensor[2][2] = val_scale*((*val_density)*(proj_vel + val_normal[1]*val_velocity[1]));
    val_Proj_Jac_Tensor[2][3] = val_scale*(val_normal[2]*(*val_density)*val_velocity[1]);
    val_Proj_Jac_Tensor[2][4] = val_scale*((*val_dRhodT)*val_velocity[1]*proj_vel);

    val_Proj_Jac_Tensor[3][0] = val_scale*(val_normal[2] + val_velocity[2]*proj_vel/(*val_betainc2));
    val_Proj_Jac_Tensor[3][1] = val_scale*(val_normal[0]*(*val_density)*val_velocity[2]);
    val_Proj_Jac_Tensor[3][2] = val_scale*(val_normal[1]*(*val_density)*val_velocity[2]);
    val_Proj_Jac_Tensor[3][3] = val_scale*((*val_density)*(proj_vel + val_normal[2]*val_velocity[2]));
    val_Proj_Jac_Tensor[3][4] = val_scale*((*val_dRhodT)*val_velocity[2]*proj_vel);

    val_Proj_Jac_Tensor[4][0] = val_scale*((*val_cp)*(*val_temperature)*proj_vel/(*val_betainc2));
    val_Proj_Jac_Tensor[4][1] = val_scale*((*val_cp)*(*val_temperature)*val_normal[0]*(*val_density));
    val_Proj_Jac_Tensor[4][2] = val_scale*((*val_cp)*(*val_temperature)*val_normal[1]*(*val_density));
    val_Proj_Jac_Tensor[4][3] = val_scale*((*val_cp)*(*val_temperature)*val_normal[2]*(*val_density));
    val_Proj_Jac_Tensor[4][4] = val_scale*((*val_cp)*((*val_temperature)*(*val_dRhodT) + (*val_density))*proj_vel);

  }
  AD::EndPassive(wasActive);
}

void CNumerics::GetPreconditioner(const su2double *val_density, const su2double *val_velocity,
                                  const su2double *val_betainc2, const su2double *val_cp,
                                  const su2double *val_temperature, const su2double *val_drhodt,
                                  su2double **val_Precon) const {
  unsigned short iDim, jDim;

  val_Precon[0][0] = 1.0/(*val_betainc2);
  for (iDim = 0; iDim < nDim; iDim++)
    val_Precon[iDim+1][0] = val_velocity[iDim]/(*val_betainc2);
  val_Precon[nDim+1][0] = (*val_cp)*(*val_temperature)/(*val_betainc2);

  for (jDim = 0; jDim < nDim; jDim++) {
    val_Precon[0][jDim+1] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      if (iDim == jDim) val_Precon[iDim+1][jDim+1] = (*val_density);
      else val_Precon[iDim+1][jDim+1] = 0.0;
    }
    val_Precon[nDim+1][jDim+1] = 0.0;
  }

  val_Precon[0][nDim+1] = (*val_drhodt);
  for (iDim = 0; iDim < nDim; iDim++)
    val_Precon[iDim+1][nDim+1] = val_velocity[iDim]*(*val_drhodt);
  val_Precon[nDim+1][nDim+1] = (*val_cp)*((*val_drhodt)*(*val_temperature) + (*val_density));

}

void CNumerics::GetPreconditionedProjJac(const su2double *val_density, const su2double *val_lambda,
                                         const su2double *val_betainc2, const su2double *val_normal,
                                         su2double **val_invPrecon_A) const {
  unsigned short iDim, jDim, kDim;

  val_invPrecon_A[0][0] = val_lambda[nDim]/2.0 + val_lambda[nDim+1]/2.0;
  for (iDim = 0; iDim < nDim; iDim++)
    val_invPrecon_A[iDim+1][0] = val_normal[iDim]*(-val_lambda[nDim] + val_lambda[nDim+1])/(2.0*sqrt((*val_betainc2))*(*val_density));
  val_invPrecon_A[nDim+1][0] = 0.0;

  for (jDim = 0; jDim < nDim; jDim++) {
    val_invPrecon_A[0][jDim+1] = sqrt((*val_betainc2))*val_normal[jDim]*(*val_density)*(-val_lambda[nDim] + val_lambda[nDim+1])/(2.0);
    for (iDim = 0; iDim < nDim; iDim++) {
      if (iDim == jDim) {
        val_invPrecon_A[iDim+1][jDim+1] = (val_lambda[nDim] + val_lambda[nDim+1])*val_normal[iDim]*val_normal[iDim]/2.0;
        for (kDim = 0; kDim < nDim; kDim++) {
          if (kDim != iDim) val_invPrecon_A[iDim+1][jDim+1] += 2.0*val_lambda[0]*val_normal[kDim]*val_normal[kDim];
        }
      }
      else {
        val_invPrecon_A[iDim+1][jDim+1] = val_normal[iDim]*val_normal[jDim]*(-2.0*val_lambda[0] + val_lambda[nDim] + val_lambda[nDim+1])/(2.0);
      }
    }
    val_invPrecon_A[nDim+1][jDim+1] = 0.0;
  }

  val_invPrecon_A[0][nDim+1] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    val_invPrecon_A[iDim+1][nDim+1] = 0.0;
  val_invPrecon_A[nDim+1][nDim+1] = val_lambda[nDim-1];

}

void CNumerics::GetPMatrix(const su2double *val_density, const su2double *val_velocity,
                           const su2double *val_soundspeed, const su2double *val_normal,
                           su2double **val_p_tensor) const {

  su2double sqvel, rhooc, rhoxc;
  //su2double c2;

  rhooc = *val_density / *val_soundspeed;
  rhoxc = *val_density * *val_soundspeed;
  //c2 = *val_soundspeed * *val_soundspeed;

  if (nDim == 2) {

    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];

    val_p_tensor[0][0]=1.0;
    val_p_tensor[0][1]=0.0;
    val_p_tensor[0][2]=0.5*rhooc;
    val_p_tensor[0][3]=0.5*rhooc;

    val_p_tensor[1][0]=val_velocity[0];
    val_p_tensor[1][1]=*val_density*val_normal[1];
    val_p_tensor[1][2]=0.5*(val_velocity[0]*rhooc+val_normal[0]**val_density);
    val_p_tensor[1][3]=0.5*(val_velocity[0]*rhooc-val_normal[0]**val_density);

    val_p_tensor[2][0]=val_velocity[1];
    val_p_tensor[2][1]=-*val_density*val_normal[0];
    val_p_tensor[2][2]=0.5*(val_velocity[1]*rhooc+val_normal[1]**val_density);
    val_p_tensor[2][3]=0.5*(val_velocity[1]*rhooc-val_normal[1]**val_density);

    val_p_tensor[3][0]=0.5*sqvel;
    val_p_tensor[3][1]=*val_density*val_velocity[0]*val_normal[1]-*val_density*val_velocity[1]*val_normal[0];
    val_p_tensor[3][2]=0.5*(0.5*sqvel*rhooc+*val_density*val_velocity[0]*val_normal[0]+*val_density*val_velocity[1]*val_normal[1]+rhoxc/Gamma_Minus_One);
    val_p_tensor[3][3]=0.5*(0.5*sqvel*rhooc-*val_density*val_velocity[0]*val_normal[0]-*val_density*val_velocity[1]*val_normal[1]+rhoxc/Gamma_Minus_One);

  }
  else {

    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];

    val_p_tensor[0][0]=val_normal[0];
    val_p_tensor[0][1]=val_normal[1];
    val_p_tensor[0][2]=val_normal[2];
    val_p_tensor[0][3]=0.5*rhooc;
    val_p_tensor[0][4]=0.5*rhooc;

    val_p_tensor[1][0]=val_velocity[0]*val_normal[0];
    val_p_tensor[1][1]=val_velocity[0]*val_normal[1]-*val_density*val_normal[2];
    val_p_tensor[1][2]=val_velocity[0]*val_normal[2]+*val_density*val_normal[1];
    val_p_tensor[1][3]=0.5*(val_velocity[0]*rhooc+*val_density*val_normal[0]);
    val_p_tensor[1][4]=0.5*(val_velocity[0]*rhooc-*val_density*val_normal[0]);

    val_p_tensor[2][0]=val_velocity[1]*val_normal[0]+*val_density*val_normal[2];
    val_p_tensor[2][1]=val_velocity[1]*val_normal[1];
    val_p_tensor[2][2]=val_velocity[1]*val_normal[2]-*val_density*val_normal[0];
    val_p_tensor[2][3]=0.5*(val_velocity[1]*rhooc+*val_density*val_normal[1]);
    val_p_tensor[2][4]=0.5*(val_velocity[1]*rhooc-*val_density*val_normal[1]);

    val_p_tensor[3][0]=val_velocity[2]*val_normal[0]-*val_density*val_normal[1];
    val_p_tensor[3][1]=val_velocity[2]*val_normal[1]+*val_density*val_normal[0];
    val_p_tensor[3][2]=val_velocity[2]*val_normal[2];
    val_p_tensor[3][3]=0.5*(val_velocity[2]*rhooc+*val_density*val_normal[2]);
    val_p_tensor[3][4]=0.5*(val_velocity[2]*rhooc-*val_density*val_normal[2]);

    val_p_tensor[4][0]=0.5*sqvel*val_normal[0]+*val_density*val_velocity[1]*val_normal[2]-*val_density*val_velocity[2]*val_normal[1];
    val_p_tensor[4][1]=0.5*sqvel*val_normal[1]-*val_density*val_velocity[0]*val_normal[2]+*val_density*val_velocity[2]*val_normal[0];
    val_p_tensor[4][2]=0.5*sqvel*val_normal[2]+*val_density*val_velocity[0]*val_normal[1]-*val_density*val_velocity[1]*val_normal[0];
    val_p_tensor[4][3]=0.5*(0.5*sqvel*rhooc+*val_density*(val_velocity[0]*val_normal[0]+val_velocity[1]*val_normal[1]+val_velocity[2]*val_normal[2])+rhoxc/Gamma_Minus_One);
    val_p_tensor[4][4]=0.5*(0.5*sqvel*rhooc-*val_density*(val_velocity[0]*val_normal[0]+val_velocity[1]*val_normal[1]+val_velocity[2]*val_normal[2])+rhoxc/Gamma_Minus_One);

  }

}

void CNumerics::GetPMatrix(const su2double *val_density, const su2double *val_velocity,
                           const su2double *val_soundspeed, const su2double *val_enthalpy,
                           const su2double *val_chi, const su2double *val_kappa,
                           const su2double *val_normal, su2double **val_p_tensor) const {

  su2double sqvel, rhooc, zeta;
  //su2double rhoxc, c2;

  rhooc = *val_density / *val_soundspeed;
  //rhoxc = *val_density * *val_soundspeed;
  //c2 = *val_soundspeed * *val_soundspeed;

  if (nDim == 2) {
    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];
    zeta = sqvel - (*val_kappa*0.5*sqvel + *val_chi)/(*val_kappa);

    val_p_tensor[0][0] = 1.0;
    val_p_tensor[0][1]=0.0;
    val_p_tensor[0][2]=0.5*rhooc;
    val_p_tensor[0][3]=0.5*rhooc;

    val_p_tensor[1][0]=val_velocity[0];
    val_p_tensor[1][1]=*val_density*val_normal[1];
    val_p_tensor[1][2]=0.5*(val_velocity[0]*rhooc+val_normal[0]**val_density);
    val_p_tensor[1][3]=0.5*(val_velocity[0]*rhooc-val_normal[0]**val_density);

    val_p_tensor[2][0]=val_velocity[1];
    val_p_tensor[2][1]=-*val_density*val_normal[0];
    val_p_tensor[2][2]=0.5*(val_velocity[1]*rhooc+val_normal[1]**val_density);
    val_p_tensor[2][3]=0.5*(val_velocity[1]*rhooc-val_normal[1]**val_density);

    val_p_tensor[3][0]= zeta;
    val_p_tensor[3][1]=*val_density*val_velocity[0]*val_normal[1]-*val_density*val_velocity[1]*val_normal[0];
    val_p_tensor[3][2]=0.5*(*val_enthalpy*rhooc+*val_density*val_velocity[0]*val_normal[0]+*val_density*val_velocity[1]*val_normal[1]);
    val_p_tensor[3][3]=0.5*(*val_enthalpy*rhooc-*val_density*val_velocity[0]*val_normal[0]-*val_density*val_velocity[1]*val_normal[1]);
  }
  else {
    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];
    zeta = sqvel - (*val_kappa*0.5*sqvel + *val_chi)/(*val_kappa);

    val_p_tensor[0][0]=val_normal[0];
    val_p_tensor[0][1]=val_normal[1];
    val_p_tensor[0][2]=val_normal[2];
    val_p_tensor[0][3]=0.5*rhooc;
    val_p_tensor[0][4]=0.5*rhooc;

    val_p_tensor[1][0]=val_velocity[0]*val_normal[0];
    val_p_tensor[1][1]=val_velocity[0]*val_normal[1]-*val_density*val_normal[2];
    val_p_tensor[1][2]=val_velocity[0]*val_normal[2]+*val_density*val_normal[1];
    val_p_tensor[1][3]=0.5*(val_velocity[0]*rhooc+*val_density*val_normal[0]);
    val_p_tensor[1][4]=0.5*(val_velocity[0]*rhooc-*val_density*val_normal[0]);

    val_p_tensor[2][0]=val_velocity[1]*val_normal[0]+*val_density*val_normal[2];
    val_p_tensor[2][1]=val_velocity[1]*val_normal[1];
    val_p_tensor[2][2]=val_velocity[1]*val_normal[2]-*val_density*val_normal[0];
    val_p_tensor[2][3]=0.5*(val_velocity[1]*rhooc+*val_density*val_normal[1]);
    val_p_tensor[2][4]=0.5*(val_velocity[1]*rhooc-*val_density*val_normal[1]);

    val_p_tensor[3][0]=val_velocity[2]*val_normal[0]-*val_density*val_normal[1];
    val_p_tensor[3][1]=val_velocity[2]*val_normal[1]+*val_density*val_normal[0];
    val_p_tensor[3][2]=val_velocity[2]*val_normal[2];
    val_p_tensor[3][3]=0.5*(val_velocity[2]*rhooc+*val_density*val_normal[2]);
    val_p_tensor[3][4]=0.5*(val_velocity[2]*rhooc-*val_density*val_normal[2]);

    val_p_tensor[4][0]=zeta*val_normal[0]+*val_density*val_velocity[1]*val_normal[2]-*val_density*val_velocity[2]*val_normal[1];
    val_p_tensor[4][1]=zeta*val_normal[1]-*val_density*val_velocity[0]*val_normal[2]+*val_density*val_velocity[2]*val_normal[0];
    val_p_tensor[4][2]=zeta*val_normal[2]+*val_density*val_velocity[0]*val_normal[1]-*val_density*val_velocity[1]*val_normal[0];
    val_p_tensor[4][3]=0.5*(*val_enthalpy*rhooc+*val_density*(val_velocity[0]*val_normal[0]+val_velocity[1]*val_normal[1]+val_velocity[2]*val_normal[2]));
    val_p_tensor[4][4]=0.5*(*val_enthalpy*rhooc-*val_density*(val_velocity[0]*val_normal[0]+val_velocity[1]*val_normal[1]+val_velocity[2]*val_normal[2]));
  }

}

void CNumerics::GetPMatrix_inv(const su2double *val_density, const su2double *val_velocity,
                               const su2double *val_soundspeed, const su2double *val_normal,
                               su2double **val_invp_tensor) const {

  su2double rhoxc, c2, gm1, k0orho, k1orho, gm1_o_c2, gm1_o_rhoxc, sqvel;

  rhoxc = *val_density * *val_soundspeed;
  c2 = *val_soundspeed * *val_soundspeed;
  gm1 = Gamma_Minus_One;
  k0orho = val_normal[0] / *val_density;
  k1orho = val_normal[1] / *val_density;
  gm1_o_c2 = gm1/c2;
  gm1_o_rhoxc = gm1/rhoxc;

  if (nDim == 3) {

    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];

    val_invp_tensor[0][0]=val_normal[0]-val_normal[2]*val_velocity[1] / *val_density+val_normal[1]*val_velocity[2] / *val_density-val_normal[0]*0.5*gm1*sqvel/c2;
    val_invp_tensor[0][1]=val_normal[0]*gm1*val_velocity[0]/c2;
    val_invp_tensor[0][2]=val_normal[2] / *val_density+val_normal[0]*gm1*val_velocity[1]/c2;
    val_invp_tensor[0][3]=-val_normal[1] / *val_density+val_normal[0]*gm1*val_velocity[2]/c2;
    val_invp_tensor[0][4]=-val_normal[0]*gm1/c2;

    val_invp_tensor[1][0]=val_normal[1]+val_normal[2]*val_velocity[0] / *val_density-val_normal[0]*val_velocity[2] / *val_density-val_normal[1]*0.5*gm1*sqvel/c2;
    val_invp_tensor[1][1]=-val_normal[2] / *val_density+val_normal[1]*gm1*val_velocity[0]/c2;
    val_invp_tensor[1][2]=val_normal[1]*gm1*val_velocity[1]/c2;
    val_invp_tensor[1][3]=val_normal[0] / *val_density+val_normal[1]*gm1*val_velocity[2]/c2;
    val_invp_tensor[1][4]=-val_normal[1]*gm1/c2;

    val_invp_tensor[2][0]=val_normal[2]-val_normal[1]*val_velocity[0] / *val_density+val_normal[0]*val_velocity[1] / *val_density-val_normal[2]*0.5*gm1*sqvel/c2;
    val_invp_tensor[2][1]=val_normal[1] / *val_density+val_normal[2]*gm1*val_velocity[0]/c2;
    val_invp_tensor[2][2]=-val_normal[0] / *val_density+val_normal[2]*gm1*val_velocity[1]/c2;
    val_invp_tensor[2][3]=val_normal[2]*gm1*val_velocity[2]/c2;
    val_invp_tensor[2][4]=-val_normal[2]*gm1/c2;

    val_invp_tensor[3][0]=-(val_normal[0]*val_velocity[0]+val_normal[1]*val_velocity[1]+val_normal[2]*val_velocity[2]) / *val_density+0.5*gm1*sqvel/rhoxc;
    val_invp_tensor[3][1]=val_normal[0] / *val_density-gm1*val_velocity[0]/rhoxc;
    val_invp_tensor[3][2]=val_normal[1] / *val_density-gm1*val_velocity[1]/rhoxc;
    val_invp_tensor[3][3]=val_normal[2] / *val_density-gm1*val_velocity[2]/rhoxc;
    val_invp_tensor[3][4]=Gamma_Minus_One/rhoxc;

    val_invp_tensor[4][0]=(val_normal[0]*val_velocity[0]+val_normal[1]*val_velocity[1]+val_normal[2]*val_velocity[2]) / *val_density+0.5*gm1*sqvel/rhoxc;
    val_invp_tensor[4][1]=-val_normal[0] / *val_density-gm1*val_velocity[0]/rhoxc;
    val_invp_tensor[4][2]=-val_normal[1] / *val_density-gm1*val_velocity[1]/rhoxc;
    val_invp_tensor[4][3]=-val_normal[2] / *val_density-gm1*val_velocity[2]/rhoxc;
    val_invp_tensor[4][4]=Gamma_Minus_One/rhoxc;

  }
  else {

    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];

    val_invp_tensor[0][0] = 1.0-0.5*gm1_o_c2*sqvel;
    val_invp_tensor[0][1]=gm1_o_c2*val_velocity[0];
    val_invp_tensor[0][2]=gm1_o_c2*val_velocity[1];
    val_invp_tensor[0][3]=-gm1_o_c2;

    val_invp_tensor[1][0]=-k1orho*val_velocity[0]+k0orho*val_velocity[1];
    val_invp_tensor[1][1]=k1orho;
    val_invp_tensor[1][2]=-k0orho;
    val_invp_tensor[1][3]=0.0;

    val_invp_tensor[2][0]=-k0orho*val_velocity[0]-k1orho*val_velocity[1]+0.5*gm1_o_rhoxc*sqvel;
    val_invp_tensor[2][1]=k0orho-gm1_o_rhoxc*val_velocity[0];
    val_invp_tensor[2][2]=k1orho-gm1_o_rhoxc*val_velocity[1];
    val_invp_tensor[2][3]=gm1_o_rhoxc;

    val_invp_tensor[3][0]=k0orho*val_velocity[0]+k1orho*val_velocity[1]+0.5*gm1_o_rhoxc*sqvel;
    val_invp_tensor[3][1]=-k0orho-gm1_o_rhoxc*val_velocity[0];
    val_invp_tensor[3][2]=-k1orho-gm1_o_rhoxc*val_velocity[1];
    val_invp_tensor[3][3]=gm1_o_rhoxc;

  }
}

void CNumerics::GetPMatrix_inv(su2double **val_invp_tensor, const su2double *val_density,
                               const su2double *val_velocity, const su2double *val_soundspeed,
                               const su2double *val_chi, const su2double *val_kappa, const su2double *val_normal) const {

  su2double rhoxc, c2, k0orho, k1orho, sqvel, k_o_c2, k_o_rhoxc, dp_drho;

  rhoxc = *val_density * *val_soundspeed;
  c2 = *val_soundspeed * *val_soundspeed;
  k0orho = val_normal[0] / *val_density;
  k1orho = val_normal[1] / *val_density;
  k_o_c2 = (*val_kappa)/c2;
  k_o_rhoxc = (*val_kappa)/rhoxc;

  if (nDim == 3) {
    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];
    dp_drho = *val_chi + 0.5*sqvel*(*val_kappa);

    val_invp_tensor[0][0]=val_normal[0]-val_normal[2]*val_velocity[1] / *val_density + val_normal[1]*val_velocity[2] / *val_density - val_normal[0]*dp_drho/c2;
    val_invp_tensor[0][1]=val_normal[0]*val_velocity[0]*k_o_c2;
    val_invp_tensor[0][2]=val_normal[2] / *val_density + val_normal[0]*val_velocity[1]*k_o_c2;
    val_invp_tensor[0][3]=-val_normal[1] / *val_density + val_normal[0]*val_velocity[2]*k_o_c2;
    val_invp_tensor[0][4]=-val_normal[0]*k_o_c2;

    val_invp_tensor[1][0]=val_normal[1]+val_normal[2]*val_velocity[0] / *val_density - val_normal[0]*val_velocity[2] / *val_density - val_normal[1]*dp_drho/c2;
    val_invp_tensor[1][1]=-val_normal[2] / *val_density + val_normal[1]*val_velocity[0]*k_o_c2;
    val_invp_tensor[1][2]=val_normal[1]*val_velocity[1]*k_o_c2;
    val_invp_tensor[1][3]=val_normal[0] / *val_density + val_normal[1]*val_velocity[2]*k_o_c2;
    val_invp_tensor[1][4]=-val_normal[1]*k_o_c2;

    val_invp_tensor[2][0]=val_normal[2]-val_normal[1]*val_velocity[0] / *val_density + val_normal[0]*val_velocity[1] / *val_density - val_normal[2]*dp_drho/c2;
    val_invp_tensor[2][1]=val_normal[1] / *val_density + val_normal[2]*val_velocity[0]*k_o_c2;
    val_invp_tensor[2][2]=-val_normal[0] / *val_density + val_normal[2]*val_velocity[1]*k_o_c2;
    val_invp_tensor[2][3]=val_normal[2]*val_velocity[2]*k_o_c2;
    val_invp_tensor[2][4]=-val_normal[2]*k_o_c2;

    val_invp_tensor[3][0]=-(val_normal[0]*val_velocity[0]+val_normal[1]*val_velocity[1]+val_normal[2]*val_velocity[2]) / *val_density+ dp_drho/rhoxc;
    val_invp_tensor[3][1]=val_normal[0] / *val_density - val_velocity[0]*k_o_rhoxc;
    val_invp_tensor[3][2]=val_normal[1] / *val_density- val_velocity[1]*k_o_rhoxc;
    val_invp_tensor[3][3]=val_normal[2] / *val_density- val_velocity[2]*k_o_rhoxc;
    val_invp_tensor[3][4]= k_o_rhoxc;

    val_invp_tensor[4][0]=(val_normal[0]*val_velocity[0]+val_normal[1]*val_velocity[1]+val_normal[2]*val_velocity[2]) / *val_density+ dp_drho/rhoxc;
    val_invp_tensor[4][1]=-val_normal[0] / *val_density- val_velocity[0]*k_o_rhoxc;
    val_invp_tensor[4][2]=-val_normal[1] / *val_density- val_velocity[1]*k_o_rhoxc;
    val_invp_tensor[4][3]=-val_normal[2] / *val_density- val_velocity[2]*k_o_rhoxc;
    val_invp_tensor[4][4]= k_o_rhoxc;
  }
  else {
    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];
    dp_drho = *val_chi + 0.5*sqvel*(*val_kappa);

    val_invp_tensor[0][0] = 1.0 - dp_drho/c2;
    val_invp_tensor[0][1]= k_o_c2*val_velocity[0];
    val_invp_tensor[0][2]= k_o_c2*val_velocity[1];
    val_invp_tensor[0][3]=-k_o_c2;

    val_invp_tensor[1][0]=-k1orho*val_velocity[0]+k0orho*val_velocity[1];
    val_invp_tensor[1][1]=k1orho;
    val_invp_tensor[1][2]=-k0orho;
    val_invp_tensor[1][3]=0.0;

    val_invp_tensor[2][0]=-k0orho*val_velocity[0]-k1orho*val_velocity[1] + dp_drho/rhoxc;
    val_invp_tensor[2][1]=k0orho - k_o_rhoxc*val_velocity[0];
    val_invp_tensor[2][2]=k1orho - k_o_rhoxc*val_velocity[1];
    val_invp_tensor[2][3]=k_o_rhoxc;

    val_invp_tensor[3][0]=k0orho*val_velocity[0]+k1orho*val_velocity[1] + dp_drho/rhoxc;
    val_invp_tensor[3][1]=-k0orho - k_o_rhoxc*val_velocity[0];
    val_invp_tensor[3][2]=-k1orho - k_o_rhoxc*val_velocity[1];
    val_invp_tensor[3][3]= k_o_rhoxc;
  }
}

void CNumerics::GetinvRinvPe(su2double Beta2, su2double val_enthalpy,
                             su2double val_soundspeed, su2double val_density,
                             const su2double* val_velocity, su2double **invRinvPe) const {

  su2double sqvel;
  su2double factor = 1.0/(val_soundspeed*val_soundspeed*Beta2);

  if (nDim == 2) {

    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];

    invRinvPe[0][0] = factor;
    invRinvPe[0][1] = 0.0;
    invRinvPe[0][2] = 0.0;
    invRinvPe[0][3] = -val_density/Gamma;

    invRinvPe[1][0] = val_velocity[0]*factor;
    invRinvPe[1][1] = val_density;
    invRinvPe[1][2] = 0.0;
    invRinvPe[1][3] = -val_density*val_velocity[0]/Gamma;

    invRinvPe[2][0] = val_velocity[1]*factor;
    invRinvPe[2][1] = 0;
    invRinvPe[2][2] = val_density;
    invRinvPe[2][3] = -val_density*val_velocity[1]/Gamma;

    invRinvPe[3][0] = val_enthalpy*factor;
    invRinvPe[3][1] = val_density*val_velocity[0];
    invRinvPe[3][2] = val_density*val_velocity[1];
    invRinvPe[3][3] = -val_density*sqvel/(2.0*Gamma);
  }
  else {

    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];

    invRinvPe[0][0] =  factor;
    invRinvPe[0][1] = 0.0;
    invRinvPe[0][2] = 0.0;
    invRinvPe[0][3] = 0.0;
    invRinvPe[0][4] = -val_density/Gamma;

    invRinvPe[1][0] = val_velocity[0]*factor;
    invRinvPe[1][1] = val_density;
    invRinvPe[1][2] = 0.0;
    invRinvPe[1][3] = 0.0;
    invRinvPe[1][4] = -val_density*val_velocity[0]/Gamma;

    invRinvPe[2][0] = val_velocity[1]*factor;
    invRinvPe[2][1] = 0;
    invRinvPe[2][2] = val_density;
    invRinvPe[2][3] = 0.0;
    invRinvPe[2][4] = -val_density*val_velocity[1]/Gamma;


    invRinvPe[3][0] = val_velocity[2]*factor;
    invRinvPe[3][1] = 0;
    invRinvPe[3][2] = 0;
    invRinvPe[3][3] = val_density;
    invRinvPe[3][4] = -val_density*val_velocity[2]/Gamma;

    invRinvPe[4][0] = val_enthalpy*factor;
    invRinvPe[4][1] = val_density*val_velocity[0];
    invRinvPe[4][2] = val_density*val_velocity[1];
    invRinvPe[4][3] = val_density*val_velocity[2];
    invRinvPe[4][4] = -val_density*sqvel/(2.0*Gamma);

  }

}

void CNumerics::GetRMatrix(su2double val_pressure, su2double val_soundspeed, su2double val_density,
                           const su2double* val_velocity, su2double **R_Matrix) const {

  su2double sqvel;
  //su2double factor = 1.0/(val_soundspeed*val_soundspeed*1);
  su2double gm1 = Gamma - 1.0;

  if (nDim == 2) {

    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1];

    R_Matrix[0][0] =  0.5*gm1*sqvel;
    R_Matrix[0][1] = -val_velocity[0]*gm1;
    R_Matrix[0][2] = -val_velocity[1]*gm1;
    R_Matrix[0][3] = gm1;

    R_Matrix[1][0] = -val_velocity[0]/val_density;
    R_Matrix[1][1] = 1.0/val_density;
    R_Matrix[1][2] = 0.0;
    R_Matrix[1][3] = 0.0;

    R_Matrix[2][0] = -val_velocity[1]/val_density;
    R_Matrix[2][1] = 0.0;
    R_Matrix[2][2] = 1.0/val_density;
    R_Matrix[2][3] = 0.0;

    R_Matrix[3][0] = 0.5*gm1*sqvel/val_pressure - Gamma/val_density;
    R_Matrix[3][1] = -gm1*val_velocity[0]/val_pressure;
    R_Matrix[3][2] = -gm1*val_velocity[1]/val_pressure;
    R_Matrix[3][3] = gm1/val_pressure;
  }
  else {

    sqvel = val_velocity[0]*val_velocity[0]+val_velocity[1]*val_velocity[1]+val_velocity[2]*val_velocity[2];

    R_Matrix[0][0] =  0.5*gm1*sqvel;
    R_Matrix[0][1] = -val_velocity[0]*gm1;
    R_Matrix[0][2] = -val_velocity[1]*gm1;
    R_Matrix[0][3] = -val_velocity[2]*gm1;
    R_Matrix[0][4] = gm1;

    R_Matrix[1][0] = -val_velocity[0]/val_density;
    R_Matrix[1][1] = 1.0/val_density;
    R_Matrix[1][2] = 0.0;
    R_Matrix[1][3] = 0.0;
    R_Matrix[1][4] = 0.0;

    R_Matrix[2][0] = -val_velocity[1]/val_density;
    R_Matrix[2][1] = 0.0;
    R_Matrix[2][2] = 1.0/val_density;
    R_Matrix[2][3] = 0.0;
    R_Matrix[2][4] = 0.0;

    R_Matrix[3][0] = -val_velocity[2]/val_density;
    R_Matrix[3][1] = 0.0;
    R_Matrix[3][2] = 0.0;
    R_Matrix[3][3] = 1.0/val_density;
    R_Matrix[3][4] = 0.0;

    R_Matrix[4][0] = 0.5*gm1*sqvel/val_pressure - Gamma/val_density;
    R_Matrix[4][1] = -gm1*val_velocity[0]/val_pressure;
    R_Matrix[4][2] = -gm1*val_velocity[1]/val_pressure;
    R_Matrix[4][3] = -gm1*val_velocity[2]/val_pressure;
    R_Matrix[4][4] = gm1/val_pressure;

  }

}

void CNumerics::GetRMatrix(su2double val_soundspeed, su2double val_density, su2double **R_Matrix) const {

  su2double cc, rhoc;
  cc = val_soundspeed*val_soundspeed;
  rhoc = val_density*val_soundspeed;

  if (nDim == 2) {
    R_Matrix[0][0] = -1.0/cc;
    R_Matrix[0][1] = 0.0;
    R_Matrix[0][2] = 0.5/cc;
    R_Matrix[0][3] = 0.5/cc;

    R_Matrix[1][0] = 0.0;
    R_Matrix[1][1] = 0.0;
    R_Matrix[1][2] = 0.5/rhoc;
    R_Matrix[1][3] = -0.5/rhoc;

    R_Matrix[2][0] = 0.0;
    R_Matrix[2][1] = 1.0/rhoc;
    R_Matrix[2][2] = 0.0;
    R_Matrix[2][3] = 0.0;

    R_Matrix[3][0] = 0.0;
    R_Matrix[3][1] = 0.0;
    R_Matrix[3][2] = 0.5;
    R_Matrix[3][3] = 0.5;

  }
  else {

    R_Matrix[0][0] = -1.0/cc;
    R_Matrix[0][1] = 0.0;
    R_Matrix[0][2] = 0.0;
    R_Matrix[0][3] = 0.5/cc;
    R_Matrix[0][4] = 0.5/cc;

    R_Matrix[1][0] = 0.0;
    R_Matrix[1][1] = 0.0;
    R_Matrix[1][2] = 0.0;
    R_Matrix[1][3] = 0.5/rhoc;
    R_Matrix[1][4] = -0.5/rhoc;

    R_Matrix[2][0] = 0.0;
    R_Matrix[2][1] = 1.0/rhoc;
    R_Matrix[2][2] = 0.0;
    R_Matrix[2][3] = 0.0;
    R_Matrix[2][4] = 0.0;

    R_Matrix[3][0] = 0.0;
    R_Matrix[3][1] = 0.0;
    R_Matrix[3][2] = 1.0/rhoc;
    R_Matrix[3][3] = 0.0;
    R_Matrix[3][4] = 0.0;

    R_Matrix[4][0] = 0.0;
    R_Matrix[4][1] = 0.0;
    R_Matrix[4][2] = 0.0;
    R_Matrix[4][3] = 0.5;
    R_Matrix[4][4] = 0.5;

  }

}

void CNumerics::GetLMatrix(su2double val_soundspeed, su2double val_density, su2double **L_Matrix) const {

  su2double cc, rhoc;
  cc = val_soundspeed*val_soundspeed;
  rhoc = val_density*val_soundspeed;
  if (nDim == 2) {

    L_Matrix[0][0] = -cc;
    L_Matrix[0][1] = 0.0;
    L_Matrix[0][2] = 0.0;
    L_Matrix[0][3] = 1.0;

    L_Matrix[1][0] = 0.0;
    L_Matrix[1][1] = 0.0;
    L_Matrix[1][2] = rhoc;
    L_Matrix[1][3] = 0.0;

    L_Matrix[2][0] = 0.0;
    L_Matrix[2][1] = rhoc;
    L_Matrix[2][2] = 0.0;
    L_Matrix[2][3] = 1.0;

    L_Matrix[3][0] = 0.0;
    L_Matrix[3][1] = -rhoc;
    L_Matrix[3][2] = 0.0;
    L_Matrix[3][3] = 1.0;
  }
  else {

    L_Matrix[0][0] = -cc;
    L_Matrix[0][1] = 0.0;
    L_Matrix[0][2] = 0.0;
    L_Matrix[0][3] = 0.0;
    L_Matrix[0][4] = 1.0;

    L_Matrix[1][0] = 0.0;
    L_Matrix[1][1] = 0.0;
    L_Matrix[1][2] = rhoc;
    L_Matrix[1][3] = 0.0;
    L_Matrix[1][4] = 0.0;

    L_Matrix[2][0] = 0.0;
    L_Matrix[2][1] = 0.0;
    L_Matrix[2][2] = 0.0;
    L_Matrix[2][3] = rhoc;
    L_Matrix[2][4] = 0.0;

    L_Matrix[3][0] = 0.0;
    L_Matrix[3][1] = rhoc;
    L_Matrix[3][2] = 0.0;
    L_Matrix[3][3] = 0.0;
    L_Matrix[3][4] = 1.0;

    L_Matrix[4][0] = 0.0;
    L_Matrix[4][1] = -rhoc;
    L_Matrix[4][2] = 0.0;
    L_Matrix[4][3] = 0.0;
    L_Matrix[4][4] = 1.0;

  }

}

void CNumerics::ComputeResJacobianGiles(CFluidModel *FluidModel, su2double pressure, su2double density, const su2double *turboVel, su2double alphaInBC, su2double gammaInBC,  su2double **R_c, su2double **R_c_inv){
  su2double rhoc, cc;
  su2double dhdrho_P, dhdP_rho, dsdrho_P,dsdP_rho;

  FluidModel->ComputeDerivativeNRBC_Prho(pressure, density);
  cc   = FluidModel->GetSoundSpeed2();
  rhoc = density*sqrt(cc);


  dhdrho_P  = FluidModel->Getdhdrho_P();
  dhdP_rho  = FluidModel->GetdhdP_rho();
  dsdrho_P  = FluidModel->Getdsdrho_P();
  dsdP_rho  = FluidModel->GetdsdP_rho();

  if (nDim == 2) {

    R_c[0][0] = -1/cc*dsdrho_P;                   //a11
    R_c[0][1] = 0.0;                              //a12
    R_c[0][2] = 0.5/cc*dsdrho_P + 0.5*dsdP_rho;   //a13

    R_c[1][0] = 0.0;                              //a21
    R_c[1][1] = 1/rhoc;                           //a22
    R_c[1][2] = -0.5/rhoc*tan(alphaInBC);         //a23

    R_c[2][0] = -1/cc*dhdrho_P;                   //a31
    R_c[2][1] = turboVel[1]/rhoc;                 //a32
    R_c[2][2] = 0.5/cc*dhdrho_P + 0.5*turboVel[0]/rhoc + 0.5*dhdP_rho;  //a33

    InvMatrix3D(R_c, R_c_inv);
  }
  else {
    R_c[0][0] = -1/cc*dsdrho_P;                     //a11
    R_c[0][1] = 0.0;                                //a12
    R_c[0][2] = 0.0;                                //a13
    R_c[0][3] = 0.5/cc*dsdrho_P + 0.5*dsdP_rho;     //a14

    R_c[1][0] = 0.0;                                //a21
    R_c[1][1] = 1/rhoc;                             //a22
    R_c[1][2] = 0.0;                                //a23
    R_c[1][3] = -0.5/rhoc*tan(alphaInBC);           //a24

    R_c[2][0] = 0.0;                                //a31
    R_c[2][1] = 0.0;                                //a32
    R_c[2][2] = 1/rhoc;                             //a33
    R_c[2][3] = -0.5/rhoc*tan(gammaInBC);           //a34

    R_c[3][0] = -1/cc*dhdrho_P;                                          //a41
    R_c[3][1] = turboVel[1]/rhoc;                                        //a42
    R_c[3][2] = turboVel[2]/rhoc;                                        //a43
    R_c[3][3] = 0.5/cc*dhdrho_P + 0.5*turboVel[0]/rhoc + 0.5*dhdP_rho;   //a44

    InvMatrix4D(R_c, R_c_inv);
  }
}

void CNumerics::InvMatrix3D(su2double **matrix, su2double **invMatrix){

  su2double invDet;

  invDet = 1 /
      (- matrix[0][2]*matrix[1][1]*matrix[2][0] + matrix[0][1]*matrix[1][2]*matrix[2][0] + matrix[0][2]*matrix[1][0]*matrix[2][1] -
         matrix[0][0]*matrix[1][2]*matrix[2][1] - matrix[0][1]*matrix[1][0]*matrix[2][2] + matrix[0][0]*matrix[1][1]*matrix[2][2]);

  invMatrix[0][0] = invDet*( - matrix[1][2]*matrix[2][1] + matrix[1][1]*matrix[2][2] );
  invMatrix[0][1] = invDet*( + matrix[0][2]*matrix[2][1] - matrix[0][1]*matrix[2][2] );
  invMatrix[0][2] = invDet*( - matrix[0][2]*matrix[1][1] + matrix[0][1]*matrix[1][2] );

  invMatrix[1][0] = invDet*( + matrix[1][2]*matrix[2][0] - matrix[1][0]*matrix[2][2] );
  invMatrix[1][1] = invDet*( - matrix[0][2]*matrix[2][0] + matrix[0][0]*matrix[2][2] );
  invMatrix[1][2] = invDet*( + matrix[0][2]*matrix[1][0] - matrix[0][0]*matrix[1][2] );

  invMatrix[2][0] = invDet*( - matrix[1][1]*matrix[2][0] + matrix[1][0]*matrix[2][1] );
  invMatrix[2][1] = invDet*( + matrix[0][1]*matrix[2][0] - matrix[0][0]*matrix[2][1] );
  invMatrix[2][2] = invDet*( - matrix[0][1]*matrix[1][0] + matrix[0][0]*matrix[1][1] );

}

void CNumerics::InvMatrix4D(su2double **matrix, su2double **invMatrix){
  su2double invDet;

  invDet = 1 /
      (matrix[0][3]*matrix[1][2]*matrix[2][1]*matrix[3][0] - matrix[0][2]*matrix[1][3]*matrix[2][1]*matrix[3][0] - matrix[0][3]*matrix[1][1]*matrix[2][2]*matrix[3][0] +
          matrix[0][1]*matrix[1][3]*matrix[2][2]*matrix[3][0] + matrix[0][2]*matrix[1][1]*matrix[2][3]*matrix[3][0] - matrix[0][1]*matrix[1][2]*matrix[2][3]*matrix[3][0] -
          matrix[0][3]*matrix[1][2]*matrix[2][0]*matrix[3][1] + matrix[0][2]*matrix[1][3]*matrix[2][0]*matrix[3][1] + matrix[0][3]*matrix[1][0]*matrix[2][2]*matrix[3][1] -
          matrix[0][0]*matrix[1][3]*matrix[2][2]*matrix[3][1] - matrix[0][2]*matrix[1][0]*matrix[2][3]*matrix[3][1] + matrix[0][0]*matrix[1][2]*matrix[2][3]*matrix[3][1] +
          matrix[0][3]*matrix[1][1]*matrix[2][0]*matrix[3][2] - matrix[0][1]*matrix[1][3]*matrix[2][0]*matrix[3][2] - matrix[0][3]*matrix[1][0]*matrix[2][1]*matrix[3][2] +
          matrix[0][0]*matrix[1][3]*matrix[2][1]*matrix[3][2] + matrix[0][1]*matrix[1][0]*matrix[2][3]*matrix[3][2] - matrix[0][0]*matrix[1][1]*matrix[2][3]*matrix[3][2] -
          matrix[0][2]*matrix[1][1]*matrix[2][0]*matrix[3][3] + matrix[0][1]*matrix[1][2]*matrix[2][0]*matrix[3][3] + matrix[0][2]*matrix[1][0]*matrix[2][1]*matrix[3][3] -
          matrix[0][0]*matrix[1][2]*matrix[2][1]*matrix[3][3] - matrix[0][1]*matrix[1][0]*matrix[2][2]*matrix[3][3] + matrix[0][0]*matrix[1][1]*matrix[2][2]*matrix[3][3]);

  invMatrix[0][0] = invDet*(- matrix[1][3]*matrix[2][2]*matrix[3][1] + matrix[1][2]*matrix[2][3]*matrix[3][1] + matrix[1][3]*matrix[2][1]*matrix[3][2] - matrix[1][1]*matrix[2][3]*matrix[3][2] - matrix[1][2]*matrix[2][1]*matrix[3][3] + matrix[1][1]*matrix[2][2]*matrix[3][3]) ;
  invMatrix[0][1] = invDet*(  matrix[0][3]*matrix[2][2]*matrix[3][1] - matrix[0][2]*matrix[2][3]*matrix[3][1] - matrix[0][3]*matrix[2][1]*matrix[3][2] + matrix[0][1]*matrix[2][3]*matrix[3][2] + matrix[0][2]*matrix[2][1]*matrix[3][3] - matrix[0][1]*matrix[2][2]*matrix[3][3]) ;
  invMatrix[0][2] = invDet*(- matrix[0][3]*matrix[1][2]*matrix[3][1] + matrix[0][2]*matrix[1][3]*matrix[3][1] + matrix[0][3]*matrix[1][1]*matrix[3][2] - matrix[0][1]*matrix[1][3]*matrix[3][2] - matrix[0][2]*matrix[1][1]*matrix[3][3] + matrix[0][1]*matrix[1][2]*matrix[3][3]) ;
  invMatrix[0][3] = invDet*(  matrix[0][3]*matrix[1][2]*matrix[2][1] - matrix[0][2]*matrix[1][3]*matrix[2][1] - matrix[0][3]*matrix[1][1]*matrix[2][2] + matrix[0][1]*matrix[1][3]*matrix[2][2] + matrix[0][2]*matrix[1][1]*matrix[2][3] - matrix[0][1]*matrix[1][2]*matrix[2][3]) ;

  invMatrix[1][0] = invDet*(  matrix[1][3]*matrix[2][2]*matrix[3][0] - matrix[1][2]*matrix[2][3]*matrix[3][0] - matrix[1][3]*matrix[2][0]*matrix[3][2] + matrix[1][0]*matrix[2][3]*matrix[3][2] + matrix[1][2]*matrix[2][0]*matrix[3][3] - matrix[1][0]*matrix[2][2]*matrix[3][3]) ;
  invMatrix[1][1] = invDet*(- matrix[0][3]*matrix[2][2]*matrix[3][0] + matrix[0][2]*matrix[2][3]*matrix[3][0] + matrix[0][3]*matrix[2][0]*matrix[3][2] - matrix[0][0]*matrix[2][3]*matrix[3][2] - matrix[0][2]*matrix[2][0]*matrix[3][3] + matrix[0][0]*matrix[2][2]*matrix[3][3]) ;
  invMatrix[1][2] = invDet*(  matrix[0][3]*matrix[1][2]*matrix[3][0] - matrix[0][2]*matrix[1][3]*matrix[3][0] - matrix[0][3]*matrix[1][0]*matrix[3][2] + matrix[0][0]*matrix[1][3]*matrix[3][2] + matrix[0][2]*matrix[1][0]*matrix[3][3] - matrix[0][0]*matrix[1][2]*matrix[3][3]) ;
  invMatrix[1][3] = invDet*(- matrix[0][3]*matrix[1][2]*matrix[2][0] + matrix[0][2]*matrix[1][3]*matrix[2][0] + matrix[0][3]*matrix[1][0]*matrix[2][2] - matrix[0][0]*matrix[1][3]*matrix[2][2] - matrix[0][2]*matrix[1][0]*matrix[2][3] + matrix[0][0]*matrix[1][2]*matrix[2][3]) ;

  invMatrix[2][0] = invDet*(- matrix[1][3]*matrix[2][1]*matrix[3][0] + matrix[1][1]*matrix[2][3]*matrix[3][0] + matrix[1][3]*matrix[2][0]*matrix[3][1] - matrix[1][0]*matrix[2][3]*matrix[3][1] - matrix[1][1]*matrix[2][0]*matrix[3][3] + matrix[1][0]*matrix[2][1]*matrix[3][3]) ;
  invMatrix[2][1] = invDet*(  matrix[0][3]*matrix[2][1]*matrix[3][0] - matrix[0][1]*matrix[2][3]*matrix[3][0] - matrix[0][3]*matrix[2][0]*matrix[3][1] + matrix[0][0]*matrix[2][3]*matrix[3][1] + matrix[0][1]*matrix[2][0]*matrix[3][3] - matrix[0][0]*matrix[2][1]*matrix[3][3]) ;
  invMatrix[2][2] = invDet*(- matrix[0][3]*matrix[1][1]*matrix[3][0] + matrix[0][1]*matrix[1][3]*matrix[3][0] + matrix[0][3]*matrix[1][0]*matrix[3][1] - matrix[0][0]*matrix[1][3]*matrix[3][1] - matrix[0][1]*matrix[1][0]*matrix[3][3] + matrix[0][0]*matrix[1][1]*matrix[3][3]) ;
  invMatrix[2][3] = invDet*(  matrix[0][3]*matrix[1][1]*matrix[2][0] - matrix[0][1]*matrix[1][3]*matrix[2][0] - matrix[0][3]*matrix[1][0]*matrix[2][1] + matrix[0][0]*matrix[1][3]*matrix[2][1] + matrix[0][1]*matrix[1][0]*matrix[2][3] - matrix[0][0]*matrix[1][1]*matrix[2][3]) ;

  invMatrix[3][0] = invDet*(  matrix[1][2]*matrix[2][1]*matrix[3][0] - matrix[1][1]*matrix[2][2]*matrix[3][0] - matrix[1][2]*matrix[2][0]*matrix[3][1] + matrix[1][0]*matrix[2][2]*matrix[3][1] + matrix[1][1]*matrix[2][0]*matrix[3][2] - matrix[1][0]*matrix[2][1]*matrix[3][2]) ;
  invMatrix[3][1] = invDet*(- matrix[0][2]*matrix[2][1]*matrix[3][0] + matrix[0][1]*matrix[2][2]*matrix[3][0] + matrix[0][2]*matrix[2][0]*matrix[3][1] - matrix[0][0]*matrix[2][2]*matrix[3][1] - matrix[0][1]*matrix[2][0]*matrix[3][2] + matrix[0][0]*matrix[2][1]*matrix[3][2]) ;
  invMatrix[3][2] = invDet*(  matrix[0][2]*matrix[1][1]*matrix[3][0] - matrix[0][1]*matrix[1][2]*matrix[3][0] - matrix[0][2]*matrix[1][0]*matrix[3][1] + matrix[0][0]*matrix[1][2]*matrix[3][1] + matrix[0][1]*matrix[1][0]*matrix[3][2] - matrix[0][0]*matrix[1][1]*matrix[3][2]) ;
  invMatrix[3][3] = invDet*(- matrix[0][2]*matrix[1][1]*matrix[2][0] + matrix[0][1]*matrix[1][2]*matrix[2][0] + matrix[0][2]*matrix[1][0]*matrix[2][1] - matrix[0][0]*matrix[1][2]*matrix[2][1] - matrix[0][1]*matrix[1][0]*matrix[2][2] + matrix[0][0]*matrix[1][1]*matrix[2][2]) ;


}

void CNumerics::GetCharJump(su2double val_soundspeed, su2double val_density, const su2double *delta_prim, su2double *delta_char) const{

  su2double cc, rhoc;
  cc = val_soundspeed*val_soundspeed;
  rhoc = val_density*val_soundspeed;
  if (nDim == 2) {
    delta_char[0] = -cc*delta_prim[0] + delta_prim[3];
    delta_char[1] = rhoc*delta_prim[2];
    delta_char[2] = rhoc*delta_prim[1] + delta_prim[3];                                 ;
    delta_char[3] = -rhoc*delta_prim[1] + delta_prim[3];
  }else {
    delta_char[0] = -cc*delta_prim[0] + delta_prim[4];
    delta_char[1] = rhoc*delta_prim[2];
    delta_char[2] = rhoc*delta_prim[3];
    delta_char[3] = rhoc*delta_prim[1] + delta_prim[4];
    delta_char[4] = -rhoc*delta_prim[1] + delta_prim[4];
  }
}

void CNumerics::GetPrecondJacobian(su2double Beta2, su2double r_hat, su2double s_hat, su2double t_hat,
                                   su2double rB2a2, const su2double* Lambda, const su2double *val_normal, su2double **val_absPeJac) const {

  su2double lam1, lam2, lam3, lam4;
  lam1 = Lambda[0]; lam2 = Lambda[1]; lam3 = Lambda[2]; lam4 = Lambda[3];

  if (nDim == 2) {

    val_absPeJac[0][0] =  lam3*s_hat/(2.0*t_hat) - lam4*r_hat/(2.0*t_hat);
    val_absPeJac[0][1] = -lam3*rB2a2*val_normal[0]/(2.0*t_hat) + lam4*rB2a2*val_normal[0]/(2.0*t_hat);
    val_absPeJac[0][2] = -lam3*rB2a2*val_normal[1]/(2.0*t_hat) + lam4*rB2a2*val_normal[1]/(2.0*t_hat);
    val_absPeJac[0][3] =  0.0;

    val_absPeJac[1][0] = r_hat*val_normal[0]*lam3*s_hat/(2.0*t_hat*rB2a2) + s_hat*val_normal[0]*lam4*(-r_hat)/(2.0*t_hat*rB2a2);
    val_absPeJac[1][1] = lam2*(val_normal[1]*val_normal[1]) - lam3*r_hat*val_normal[0]*val_normal[0]/(2.0*t_hat) + lam4*s_hat*val_normal[0]*val_normal[0]/(2.0*t_hat);
    val_absPeJac[1][2] = -lam2*val_normal[0]*val_normal[1] - lam3*r_hat*val_normal[0]*val_normal[1]/(2.0*t_hat) + lam4*s_hat*val_normal[0]*val_normal[1]/(2.0*t_hat);
    val_absPeJac[1][3] = 0.0;

    val_absPeJac[2][0] = lam3*r_hat*val_normal[1]*s_hat/(2.0*t_hat*rB2a2) - s_hat*val_normal[1]*lam4*r_hat/(2.0*t_hat*rB2a2);
    val_absPeJac[2][1] = -val_normal[0]*val_normal[1]*lam2 - r_hat*val_normal[1]*val_normal[0]*lam3/(2.0*t_hat) + s_hat*val_normal[0]*val_normal[1]*lam4/(2.0*t_hat);
    val_absPeJac[2][2] = val_normal[0]*val_normal[0]*lam2 -r_hat*val_normal[1]*val_normal[1]*lam3/(2.0*t_hat) + s_hat*val_normal[1]*val_normal[1]*lam4/(2.0*t_hat);
    val_absPeJac[2][3] = 0.0;

    val_absPeJac[3][0] = 0.0;
    val_absPeJac[3][1] = 0.0;
    val_absPeJac[3][2] = 0.0;
    val_absPeJac[3][3] = lam1;

  }
  else {

    su2double lam5 = Lambda[4];

    val_absPeJac[0][0] =  lam4*s_hat/(2.0*t_hat) - lam5*r_hat/(2.0*t_hat);
    val_absPeJac[0][1] = -lam4*rB2a2*val_normal[0]/(2.0*t_hat) + lam5*rB2a2*val_normal[0]/(2.0*t_hat);
    val_absPeJac[0][2] = -lam4*rB2a2*val_normal[1]/(2.0*t_hat) + lam5*rB2a2*val_normal[1]/(2.0*t_hat);
    val_absPeJac[0][3] = -lam4*rB2a2*val_normal[2]/(2.0*t_hat) + lam5*rB2a2*val_normal[2]/(2.0*t_hat);
    val_absPeJac[0][4] =  0.0;

    val_absPeJac[1][0] = r_hat*val_normal[0]*lam4*s_hat/(2.0*t_hat*rB2a2) + s_hat*val_normal[0]*lam5*(-r_hat)/(2.0*t_hat*rB2a2);
    val_absPeJac[1][1] = lam2*(val_normal[2]*val_normal[2] + val_normal[1]*val_normal[1]) - lam4*r_hat*val_normal[0]*val_normal[0]/(2.0*t_hat) + lam5*s_hat*val_normal[0]*val_normal[0]/(2.0*t_hat);
    val_absPeJac[1][2] = -lam2*val_normal[0]*val_normal[1] - lam4*r_hat*val_normal[0]*val_normal[1]/(2.0*t_hat) + lam5*s_hat*val_normal[0]*val_normal[1]/(2.0*t_hat);
    val_absPeJac[1][3] = -lam2*val_normal[0]*val_normal[2] - lam4*r_hat*val_normal[0]*val_normal[2]/(2.0*t_hat) + lam5*s_hat*val_normal[0]*val_normal[2]/(2.0*t_hat);
    val_absPeJac[1][4] = 0.0;

    val_absPeJac[2][0] = lam4*r_hat*val_normal[1]*s_hat/(2.0*t_hat*rB2a2) - s_hat*val_normal[1]*lam5*r_hat/(2.0*t_hat*rB2a2);
    val_absPeJac[2][1] = -val_normal[0]*val_normal[1]*lam2 - r_hat*val_normal[1]*val_normal[0]*lam4/(2.0*t_hat) + s_hat*val_normal[0]*val_normal[1]*lam5/(2.0*t_hat);
    val_absPeJac[2][2] = val_normal[0]*val_normal[0]*lam2 + val_normal[2]*val_normal[2]*lam3 -r_hat*val_normal[1]*val_normal[1]*lam4/(2.0*t_hat) + s_hat*val_normal[1]*val_normal[1]*lam5/(2.0*t_hat);
    val_absPeJac[2][3] = -val_normal[2]*val_normal[1]*lam2 - r_hat*val_normal[2]*val_normal[1]*lam4/(2.0*t_hat) + s_hat*lam5*val_normal[1]*val_normal[2]/(2.0*t_hat);
    val_absPeJac[2][4] = 0.0;

    val_absPeJac[3][0] = r_hat*s_hat*val_normal[2]*lam4/(2.0*t_hat*rB2a2) - r_hat*s_hat*val_normal[2]*lam5/(2.0*t_hat*rB2a2);
    val_absPeJac[3][1] = -val_normal[0]*val_normal[2]*lam3 - lam4*val_normal[0]*val_normal[2]*r_hat/(2.0*t_hat) + lam5*val_normal[0]*val_normal[2]*s_hat/(2.0*t_hat);
    val_absPeJac[3][2] = -val_normal[1]*val_normal[2]*lam3 - lam4*val_normal[1]*val_normal[2]*r_hat/(2.0*t_hat) + lam5*val_normal[1]*val_normal[2]*s_hat/(2.0*t_hat);
    val_absPeJac[3][3] = (val_normal[1]*val_normal[1] + val_normal[0]*val_normal[0])*lam3 - lam4*val_normal[2]*val_normal[2]*r_hat/(2.0*t_hat) + lam5*val_normal[2]*val_normal[2]*s_hat/(2.0*t_hat);
    val_absPeJac[3][4] = 0.0;

    val_absPeJac[4][0] = 0.0;
    val_absPeJac[4][1] = 0.0;
    val_absPeJac[4][2] = 0.0;
    val_absPeJac[4][3] = 0.0;
    val_absPeJac[4][4] = lam1;

  }

}

void CNumerics::GetJacInviscidLambda_fabs(const su2double *val_velocity, su2double val_soundspeed,
                                          const su2double *val_normal, su2double *val_Lambda_Vector) const {
  su2double ProjVelocity = 0;

  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    ProjVelocity += val_velocity[iDim]*val_normal[iDim];

  if (nDim == 3) {
    val_Lambda_Vector[0] = fabs(ProjVelocity);
    val_Lambda_Vector[1] = fabs(ProjVelocity);
    val_Lambda_Vector[2] = fabs(ProjVelocity);
    val_Lambda_Vector[3] = fabs(ProjVelocity + val_soundspeed);
    val_Lambda_Vector[4] = fabs(ProjVelocity - val_soundspeed);
  }
  else {
    val_Lambda_Vector[0] = fabs(ProjVelocity);
    val_Lambda_Vector[1] = fabs(ProjVelocity);
    val_Lambda_Vector[2] = fabs(ProjVelocity + val_soundspeed);
    val_Lambda_Vector[3] = fabs(ProjVelocity - val_soundspeed);
  }
}

void CNumerics::GetAdjViscousFlux_Jac(su2double Pressure_i, su2double Pressure_j, su2double Density_i, su2double Density_j,
                                      su2double ViscDens_i, su2double ViscDens_j, const su2double *Velocity_i, const su2double *Velocity_j,
                                      su2double sq_vel_i, su2double sq_vel_j,
                                      su2double XiDens_i, su2double XiDens_j, su2double **Mean_GradPhi, const su2double *Mean_GradPsiE,
                                      su2double dPhiE_dn, const su2double *Normal, const su2double *Edge_Vector, su2double dist_ij_2, su2double *val_residual_i, su2double *val_residual_j,
                                      su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji,
                                      su2double **val_Jacobian_jj, bool implicit) const {

  su2double Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
  Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5,
  Sigma_yz5, Sigma_5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;
  su2double dSigmaxx_phi1, dSigmayy_phi1, dSigmazz_phi1, dSigmaxy_phi1, dSigmaxz_phi1, dSigmayz_phi1;
  su2double dSigmaxx_phi2, dSigmayy_phi2, dSigmazz_phi2, dSigmaxy_phi2, dSigmaxz_phi2, dSigmayz_phi2;
  su2double dSigmaxx_phi3, dSigmayy_phi3, dSigmazz_phi3, dSigmaxy_phi3, dSigmaxz_phi3, dSigmayz_phi3;
  su2double dSigma5_psi5;
  unsigned short iVar, jVar;

  if (nDim == 3) {

    /*--- Residual at iPoint ---*/

    Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
    Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
    Sigma_zz = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
    Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
    Sigma_xz = ViscDens_i * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
    Sigma_yz = ViscDens_i * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
    Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
    Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
    Sigma_zz5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_i[2] * Mean_GradPsiE[2]);
    Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
    Sigma_xz5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[0]);
    Sigma_yz5 = ViscDens_i * (Velocity_i[1] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[1]);
    Sigma_5   = XiDens_i * dPhiE_dn;
    eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
    eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;

    val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy + Velocity_i[2] * Normal[2] * eta_zz
                           + (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
                           + (Velocity_i[0] * Normal[2] + Velocity_i[2] * Normal[0]) * eta_xz
                           + (Velocity_i[2] * Normal[1] + Velocity_i[1] * Normal[2]) * eta_yz
                           - (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);

    val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_i[0] * Sigma_5);
    val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_i[1] * Sigma_5);
    val_residual_i[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_i[2] * Sigma_5);
    val_residual_i[4] = (Sigma_5);

    /*--- Computation of the Jacobians at Point i---*/

    if (implicit) {

      dSigmaxx_phi1 = -FOUR3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmaxx_phi2 =   TWO3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmaxx_phi3 =   TWO3 * ViscDens_i * Edge_Vector[2]/dist_ij_2;
      dSigmayy_phi1 =   TWO3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmayy_phi2 = -FOUR3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmayy_phi3 =   TWO3 * ViscDens_i * Edge_Vector[2]/dist_ij_2;
      dSigmazz_phi1 =   TWO3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmazz_phi2 =   TWO3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmazz_phi3 = -FOUR3 * ViscDens_i * Edge_Vector[2]/dist_ij_2;
      dSigmaxy_phi1 = -ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi2 = -ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmaxy_phi3 = 0;
      dSigmaxz_phi1 = -ViscDens_i * Edge_Vector[2]/dist_ij_2;
      dSigmaxz_phi2 = 0;
      dSigmaxz_phi3 = -ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmayz_phi1 = 0;
      dSigmayz_phi2 = -ViscDens_i * Edge_Vector[2]/dist_ij_2;
      dSigmayz_phi3 = -ViscDens_i * Edge_Vector[1]/dist_ij_2;

      //      dSigmaxx5_psi5 = -ViscDens_i * ( FOUR3*Velocity_i[0]*Edge_Vector[0] -  TWO3*Velocity_i[1]*Edge_Vector[1] -  TWO3*Velocity_i[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmayy5_psi5 = -ViscDens_i * (- TWO3*Velocity_i[0]*Edge_Vector[0] + FOUR3*Velocity_i[1]*Edge_Vector[1] -  TWO3*Velocity_i[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmazz5_psi5 = -ViscDens_i * (- TWO3*Velocity_i[0]*Edge_Vector[0] -  TWO3*Velocity_i[1]*Edge_Vector[1] + FOUR3*Velocity_i[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmaxy5_psi5 = -ViscDens_i * ( Velocity_i[0]*Edge_Vector[1] + Velocity_i[1]*Edge_Vector[0] )/dist_ij_2;
      //      dSigmaxz5_psi5 = -ViscDens_i * ( Velocity_i[0]*Edge_Vector[2] + Velocity_i[2]*Edge_Vector[0] )/dist_ij_2;
      //      dSigmayz5_psi5 = -ViscDens_i * ( Velocity_i[1]*Edge_Vector[2] + Velocity_i[2]*Edge_Vector[1] )/dist_ij_2;
      dSigma5_psi5   = -XiDens_i * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] + Edge_Vector[2]*Normal[2] )/dist_ij_2;

      val_Jacobian_ii[0][0] = 0;
      val_Jacobian_ii[0][1] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi1 + Velocity_i[1]*Normal[1]*dSigmayy_phi1 + Velocity_i[2]*Normal[2]*dSigmazz_phi1
                                + (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi1
                                + (Velocity_i[0]*Normal[2] + Velocity_i[2]*Normal[0])*dSigmaxz_phi1
                                + (Velocity_i[2]*Normal[1] + Velocity_i[1]*Normal[2])*dSigmayz_phi1 );
      val_Jacobian_ii[0][2] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi2 + Velocity_i[1]*Normal[1]*dSigmayy_phi2 + Velocity_i[2]*Normal[2]*dSigmazz_phi2
                                + (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi2
                                + (Velocity_i[0]*Normal[2] + Velocity_i[2]*Normal[0])*dSigmaxz_phi2
                                + (Velocity_i[2]*Normal[1] + Velocity_i[1]*Normal[2])*dSigmayz_phi2 );
      val_Jacobian_ii[0][3] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi3 + Velocity_i[1]*Normal[1]*dSigmayy_phi3 + Velocity_i[2]*Normal[2]*dSigmazz_phi3
                                + (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi3
                                + (Velocity_i[0]*Normal[2] + Velocity_i[2]*Normal[0])*dSigmaxz_phi3
                                + (Velocity_i[2]*Normal[1] + Velocity_i[1]*Normal[2])*dSigmayz_phi3 );
      val_Jacobian_ii[0][4] = (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * dSigma5_psi5;

      val_Jacobian_ii[1][0] = 0;
      val_Jacobian_ii[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1 + Normal[2]*dSigmaxz_phi1;
      val_Jacobian_ii[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2 + Normal[2]*dSigmaxz_phi2;
      val_Jacobian_ii[1][3] = Normal[0]*dSigmaxx_phi3 + Normal[1]*dSigmaxy_phi3 + Normal[2]*dSigmaxz_phi3;
      val_Jacobian_ii[1][4] = -Velocity_i[0]*dSigma5_psi5;

      val_Jacobian_ii[2][0] = 0;
      val_Jacobian_ii[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1 + Normal[2]*dSigmayz_phi1;
      val_Jacobian_ii[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2 + Normal[2]*dSigmayz_phi2;
      val_Jacobian_ii[2][3] = Normal[0]*dSigmaxy_phi3 + Normal[1]*dSigmayy_phi3 + Normal[2]*dSigmayz_phi3;
      val_Jacobian_ii[2][4] = -Velocity_i[1]*dSigma5_psi5;

      val_Jacobian_ii[3][0] = 0;
      val_Jacobian_ii[3][1] = Normal[0]*dSigmaxz_phi1 + Normal[1]*dSigmayz_phi1 + Normal[2]*dSigmazz_phi1;
      val_Jacobian_ii[3][2] = Normal[0]*dSigmaxz_phi2 + Normal[1]*dSigmayz_phi2 + Normal[2]*dSigmazz_phi2;
      val_Jacobian_ii[3][3] = Normal[0]*dSigmaxz_phi3 + Normal[1]*dSigmayz_phi3 + Normal[2]*dSigmazz_phi3;
      val_Jacobian_ii[3][4] = -Velocity_i[2]*dSigma5_psi5;

      val_Jacobian_ii[4][0] = 0;
      val_Jacobian_ii[4][1] = 0;
      val_Jacobian_ii[4][2] = 0;
      val_Jacobian_ii[4][3] = 0;
      val_Jacobian_ii[4][4] = dSigma5_psi5;

      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_ij[iVar][jVar] = -val_Jacobian_ii[iVar][jVar];
    }

    /*--- Residual at jPoint ---*/

    Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
    Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
    Sigma_zz = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
    Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
    Sigma_xz = ViscDens_j * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
    Sigma_yz = ViscDens_j * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
    Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
    Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
    Sigma_zz5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_j[2] * Mean_GradPsiE[2]);
    Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
    Sigma_xz5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[0]);
    Sigma_yz5 = ViscDens_j * (Velocity_j[1] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[1]);
    Sigma_5   = XiDens_j * dPhiE_dn;
    eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
    eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;

    val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy + Velocity_j[2] * Normal[2] * eta_zz
                           + (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
                           + (Velocity_j[0] * Normal[2] + Velocity_j[2] * Normal[0]) * eta_xz
                           + (Velocity_j[2] * Normal[1] + Velocity_j[1] * Normal[2]) * eta_yz
                           - (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
    val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_j[0] * Sigma_5);
    val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_j[1] * Sigma_5);
    val_residual_j[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_j[2] * Sigma_5);
    val_residual_j[4] = (Sigma_5);

    /*--- Computation of the Jacobians at Point j---*/

    if (implicit) {

      dSigmaxx_phi1 = FOUR3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmaxx_phi2 = -TWO3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmaxx_phi3 = -TWO3 * ViscDens_j * Edge_Vector[2]/dist_ij_2;
      dSigmayy_phi1 = -TWO3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmayy_phi2 = FOUR3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmayy_phi3 = -TWO3 * ViscDens_j * Edge_Vector[2]/dist_ij_2;
      dSigmazz_phi1 = -TWO3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmazz_phi2 = -TWO3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmazz_phi3 = FOUR3 * ViscDens_j * Edge_Vector[2]/dist_ij_2;
      dSigmaxy_phi1 = ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi2 = ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmaxy_phi3 = 0;
      dSigmaxz_phi1 = ViscDens_j * Edge_Vector[2]/dist_ij_2;
      dSigmaxz_phi2 = 0;
      dSigmaxz_phi3 = ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmayz_phi1 = 0;
      dSigmayz_phi2 = ViscDens_j * Edge_Vector[2]/dist_ij_2;
      dSigmayz_phi3 = ViscDens_j * Edge_Vector[1]/dist_ij_2;

      //      dSigmaxx5_psi5 = ViscDens_j * ( FOUR3*Velocity_j[0]*Edge_Vector[0] -  TWO3*Velocity_j[1]*Edge_Vector[1] -  TWO3*Velocity_j[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmayy5_psi5 = ViscDens_j * (- TWO3*Velocity_j[0]*Edge_Vector[0] + FOUR3*Velocity_j[1]*Edge_Vector[1] -  TWO3*Velocity_j[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmazz5_psi5 = ViscDens_j * (- TWO3*Velocity_j[0]*Edge_Vector[0] -  TWO3*Velocity_j[1]*Edge_Vector[1] + FOUR3*Velocity_j[2]*Edge_Vector[2])/dist_ij_2;
      //      dSigmaxy5_psi5 = ViscDens_j * ( Velocity_j[0]*Edge_Vector[1] + Velocity_j[1]*Edge_Vector[0] )/dist_ij_2;
      //      dSigmaxz5_psi5 = ViscDens_j * ( Velocity_j[0]*Edge_Vector[2] + Velocity_j[2]*Edge_Vector[0] )/dist_ij_2;
      //      dSigmayz5_psi5 = ViscDens_j * ( Velocity_j[1]*Edge_Vector[2] + Velocity_j[2]*Edge_Vector[1] )/dist_ij_2;
      dSigma5_psi5   = XiDens_j * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] + Edge_Vector[2]*Normal[2] )/dist_ij_2;

      val_Jacobian_jj[0][0] = 0;
      val_Jacobian_jj[0][1] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi1 + Velocity_j[1]*Normal[1]*dSigmayy_phi1 + Velocity_j[2]*Normal[2]*dSigmazz_phi1
                                + (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi1
                                + (Velocity_j[0]*Normal[2] + Velocity_j[2]*Normal[0])*dSigmaxz_phi1
                                + (Velocity_j[2]*Normal[1] + Velocity_j[1]*Normal[2])*dSigmayz_phi1 );
      val_Jacobian_jj[0][2] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi2 + Velocity_j[1]*Normal[1]*dSigmayy_phi2 + Velocity_j[2]*Normal[2]*dSigmazz_phi2
                                + (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi2
                                + (Velocity_j[0]*Normal[2] + Velocity_j[2]*Normal[0])*dSigmaxz_phi2
                                + (Velocity_j[2]*Normal[1] + Velocity_j[1]*Normal[2])*dSigmayz_phi2 );
      val_Jacobian_jj[0][3] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi3 + Velocity_j[1]*Normal[1]*dSigmayy_phi3 + Velocity_j[2]*Normal[2]*dSigmazz_phi3
                                + (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi3
                                + (Velocity_j[0]*Normal[2] + Velocity_j[2]*Normal[0])*dSigmaxz_phi3
                                + (Velocity_j[2]*Normal[1] + Velocity_j[1]*Normal[2])*dSigmayz_phi3 );
      val_Jacobian_jj[0][4] = (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * dSigma5_psi5;

      val_Jacobian_jj[1][0] = 0;
      val_Jacobian_jj[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1 + Normal[2]*dSigmaxz_phi1;
      val_Jacobian_jj[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2 + Normal[2]*dSigmaxz_phi2;
      val_Jacobian_jj[1][3] = Normal[0]*dSigmaxx_phi3 + Normal[1]*dSigmaxy_phi3 + Normal[2]*dSigmaxz_phi3;
      val_Jacobian_jj[1][4] = -Velocity_j[0]*dSigma5_psi5;

      val_Jacobian_jj[2][0] = 0;
      val_Jacobian_jj[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1 + Normal[2]*dSigmayz_phi1;
      val_Jacobian_jj[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2 + Normal[2]*dSigmayz_phi2;
      val_Jacobian_jj[2][3] = Normal[0]*dSigmaxy_phi3 + Normal[1]*dSigmayy_phi3 + Normal[2]*dSigmayz_phi3;
      val_Jacobian_jj[2][4] = -Velocity_j[1]*dSigma5_psi5;

      val_Jacobian_jj[3][0] = 0;
      val_Jacobian_jj[3][1] = Normal[0]*dSigmaxz_phi1 + Normal[1]*dSigmayz_phi1 + Normal[2]*dSigmazz_phi1;
      val_Jacobian_jj[3][2] = Normal[0]*dSigmaxz_phi2 + Normal[1]*dSigmayz_phi2 + Normal[2]*dSigmazz_phi2;
      val_Jacobian_jj[3][3] = Normal[0]*dSigmaxz_phi3 + Normal[1]*dSigmayz_phi3 + Normal[2]*dSigmazz_phi3;
      val_Jacobian_jj[3][4] = -Velocity_j[2]*dSigma5_psi5;

      val_Jacobian_jj[4][0] = 0;
      val_Jacobian_jj[4][1] = 0;
      val_Jacobian_jj[4][2] = 0;
      val_Jacobian_jj[4][3] = 0;
      val_Jacobian_jj[4][4] = dSigma5_psi5;

      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_ji[iVar][jVar] = -val_Jacobian_jj[iVar][jVar];
    }

  }
  else if (nDim == 2) {

    /*--- Residual at iPoint ---*/

    Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
    Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
    Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
    Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1]);
    Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1]);
    Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
    Sigma_5   = XiDens_i * dPhiE_dn;
    eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;

    val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy
                           + (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
                           - (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);
    val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] - Velocity_i[0] * Sigma_5);
    val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] - Velocity_i[1] * Sigma_5);
    val_residual_i[3] = (Sigma_5);

    /*--- Computation of the Jacobians at Point i---*/

    if (implicit) {

      dSigmaxx_phi1 = -FOUR3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmaxx_phi2 =   TWO3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmayy_phi1 =   TWO3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
      dSigmayy_phi2 = -FOUR3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi1 = -ViscDens_i * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi2 = -ViscDens_i * Edge_Vector[0]/dist_ij_2;

      //      dSigmaxx5_psi5 = -ViscDens_i * ( FOUR3*Velocity_i[0]*Edge_Vector[0] -  TWO3*Velocity_i[1]*Edge_Vector[1] )/dist_ij_2;
      //      dSigmayy5_psi5 = -ViscDens_i * (- TWO3*Velocity_i[0]*Edge_Vector[0] + FOUR3*Velocity_i[1]*Edge_Vector[1] )/dist_ij_2;
      //      dSigmaxy5_psi5 = -ViscDens_i * ( Velocity_i[0]*Edge_Vector[1] + Velocity_i[1]*Edge_Vector[0] )/dist_ij_2;
      dSigma5_psi5   = -XiDens_i * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] )/dist_ij_2;

      val_Jacobian_ii[0][0] = 0;

      val_Jacobian_ii[0][1] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi1 + Velocity_i[1]*Normal[1]*dSigmayy_phi1
                                + (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi1 );
      val_Jacobian_ii[0][2] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi2 + Velocity_i[1]*Normal[1]*dSigmayy_phi2
                                + (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi2 );
      val_Jacobian_ii[0][3] = (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * dSigma5_psi5;

      val_Jacobian_ii[1][0] = 0;
      val_Jacobian_ii[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1;
      val_Jacobian_ii[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2;
      val_Jacobian_ii[1][3] = -Velocity_i[0]*dSigma5_psi5;

      val_Jacobian_ii[2][0] = 0;
      val_Jacobian_ii[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1;
      val_Jacobian_ii[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2;
      val_Jacobian_ii[2][3] = -Velocity_i[1]*dSigma5_psi5;

      val_Jacobian_ii[3][0] = 0;
      val_Jacobian_ii[3][1] = 0;
      val_Jacobian_ii[3][2] = 0;
      val_Jacobian_ii[3][3] = dSigma5_psi5;

      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_ij[iVar][jVar] = -val_Jacobian_ii[iVar][jVar];
    }

    /*--- Residual at jPoint ---*/
    Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
    Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
    Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
    Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1]);
    Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1]);
    Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
    Sigma_5   = XiDens_j * dPhiE_dn;
    eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;

    val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy
                           + (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
                           - (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
    val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1]  - Velocity_j[0] * Sigma_5);
    val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1]  - Velocity_j[1] * Sigma_5);
    val_residual_j[3] = (Sigma_5);

    /*--- Computation of the Jacobians at Point j---*/
    if (implicit) {
      dSigmaxx_phi1 = FOUR3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmaxx_phi2 = -TWO3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmayy_phi1 = -TWO3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
      dSigmayy_phi2 = FOUR3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi1 = ViscDens_j * Edge_Vector[1]/dist_ij_2;
      dSigmaxy_phi2 = ViscDens_j * Edge_Vector[0]/dist_ij_2;

      //      dSigmaxx5_psi5 = ViscDens_j * ( FOUR3*Velocity_j[0]*Edge_Vector[0] -  TWO3*Velocity_j[1]*Edge_Vector[1] )/dist_ij_2;
      //      dSigmayy5_psi5 = ViscDens_j * (- TWO3*Velocity_j[0]*Edge_Vector[0] + FOUR3*Velocity_j[1]*Edge_Vector[1] )/dist_ij_2;
      //      dSigmaxy5_psi5 = ViscDens_j * ( Velocity_j[0]*Edge_Vector[1] + Velocity_j[1]*Edge_Vector[0] )/dist_ij_2;
      dSigma5_psi5   = XiDens_j * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] )/dist_ij_2;

      val_Jacobian_jj[0][0] = 0;
      val_Jacobian_jj[0][1] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi1 + Velocity_j[1]*Normal[1]*dSigmayy_phi1
                                + (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi1 );
      val_Jacobian_jj[0][2] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi2 + Velocity_j[1]*Normal[1]*dSigmayy_phi2
                                + (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi2 );
      val_Jacobian_jj[0][3] = (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * dSigma5_psi5;

      val_Jacobian_jj[1][0] = 0;
      val_Jacobian_jj[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1;
      val_Jacobian_jj[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2;
      val_Jacobian_jj[1][3] = -Velocity_j[0]*dSigma5_psi5;

      val_Jacobian_jj[2][0] = 0;
      val_Jacobian_jj[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1;
      val_Jacobian_jj[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2;
      val_Jacobian_jj[2][3] = -Velocity_j[1]*dSigma5_psi5;

      val_Jacobian_jj[3][0] = 0;
      val_Jacobian_jj[3][1] = 0;
      val_Jacobian_jj[3][2] = 0;
      val_Jacobian_jj[3][3] = dSigma5_psi5;

      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_ji[iVar][jVar] = -val_Jacobian_jj[iVar][jVar];
    }
  }

}

void CNumerics::GetPrimitive2Conservative (const su2double *val_Mean_PrimVar,
                                           const su2double *val_Mean_SecVar,
                                           su2double **val_Jac_PC) const {
  unsigned short iVar, jVar, iDim;

  // order of primitives: T, vx, vy, vz, P, rho, h, c, MuLam, MuEddy, kt, Cp
  // order of secondary: dPdrho_e, dPde_rho, dTdrho_e, dTde_rho, dmudrho_T, dmudT_rho, dktdrho_T, dktdT_rho

  su2double vx = val_Mean_PrimVar[1];
  su2double vy = val_Mean_PrimVar[2];
  su2double vz = val_Mean_PrimVar[3];
  su2double rho = val_Mean_PrimVar[nDim+2];
  su2double P = val_Mean_PrimVar[nDim+1];
  su2double e = val_Mean_PrimVar[nDim+3] - P/rho;
  su2double dTdrho_e = val_Mean_SecVar[2];
  su2double dTde_rho = val_Mean_SecVar[3];

  su2double sqvel = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
  }

  /*--- Initialize the Jacobian matrix ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_Jac_PC[iVar][jVar] = 0.0;
    }
  }

  /*--- Primitives to conservatives Jacobian matrix : (T, vx, vy, vz, rho) --> (u1, u2, u3, u4, u5) ---*/
  if (nDim == 2) {

  val_Jac_PC[0][0] = dTdrho_e - e/rho*dTde_rho + 0.5*dTde_rho*sqvel/rho;
  val_Jac_PC[0][1] = -1/rho*dTde_rho*vx;
  val_Jac_PC[0][2] = -1/rho*dTde_rho*vy;
  val_Jac_PC[0][3] = 1/rho*dTde_rho;

  val_Jac_PC[1][0] = -vx/rho;
  val_Jac_PC[1][1] = 1/rho;
  val_Jac_PC[1][2] = 0.0;
  val_Jac_PC[1][3] = 0.0;

  val_Jac_PC[2][0] = -vy/rho;
  val_Jac_PC[2][1] = 0.0;
  val_Jac_PC[2][2] = 1/rho;
  val_Jac_PC[2][3] = 0.0;

  val_Jac_PC[3][0] = 1.0;
  val_Jac_PC[3][1] = 0.0;
  val_Jac_PC[3][2] = 0.0;
  val_Jac_PC[3][3] = 0.0;

  }
  else {

  val_Jac_PC[0][0] = dTdrho_e - e/rho*dTde_rho + 0.5*dTde_rho*sqvel/rho;
  val_Jac_PC[0][1] = -1/rho*dTde_rho*vx;
  val_Jac_PC[0][2] = -1/rho*dTde_rho*vy;
  val_Jac_PC[0][3] = -1/rho*dTde_rho*vz;
  val_Jac_PC[0][4] = 1/rho*dTde_rho;

  val_Jac_PC[1][0] = -vx/rho;
  val_Jac_PC[1][1] = 1/rho;
  val_Jac_PC[1][2] = 0.0;
  val_Jac_PC[1][3] = 0.0;
  val_Jac_PC[1][4] = 0.0;

  val_Jac_PC[2][0] = -vy/rho;
  val_Jac_PC[2][1] = 0.0;
  val_Jac_PC[2][2] = 1/rho;
  val_Jac_PC[2][3] = 0.0;
  val_Jac_PC[2][4] = 0.0;

  val_Jac_PC[3][0] = -vz/rho;
  val_Jac_PC[3][1] = 0.0;
  val_Jac_PC[3][2] = 0.0;
  val_Jac_PC[3][3] = 1/rho;
  val_Jac_PC[3][4] = 0.0;

  val_Jac_PC[4][0] = 1.0;
  val_Jac_PC[4][1] = 0.0;
  val_Jac_PC[4][2] = 0.0;
  val_Jac_PC[4][3] = 0.0;
  val_Jac_PC[4][4] = 0.0;

  }
}

void CNumerics::CreateBasis(const su2double *val_Normal, su2double* l, su2double* m) {

  unsigned short iDim;
  su2double modm, modl;

  if (nDim ==2){

    /*--- Multiply Normal by [0 -1; 1 0] rotation matrix ---*/
    l[0] = -val_Normal[1];
    l[1] = val_Normal[0];

    /*--- Set m matrix to zero ---*/
    m[0] = 0.0;
    m[1] = 0.0;

    /*--- Normalize ---*/
    modl = 0.0;
    for (iDim = 0; iDim <nDim; iDim ++)
      modl += l[iDim]*l[iDim];
    modl = sqrt(modl);
    for (iDim =0; iDim<nDim; iDim++)
      l[iDim] = l[iDim]/modl;

  } else {

    /*--- Define l as a vector in the plane normal to the supplied vector ---*/
    l[0] = 0.0;
    l[1] = -val_Normal[2];
    l[2] = val_Normal[1];

    /*--- Check for the zero vector and re-assign if needed ---*/
    if (l[0] == 0.0 && l[1] == 0.0 && l[2] == 0.0) {
      l[0] = -val_Normal[2];
      l[1] = 0.0;
      l[2] = val_Normal[0];
    }

    /*--- Take vector product of n * l to make m ---*/
    m[0] = val_Normal[1]*l[2] - val_Normal[2]*l[1];
    m[1] = val_Normal[2]*l[0] - val_Normal[0]*l[2];
    m[2] = val_Normal[0]*l[1] - val_Normal[1]*l[0];

    /*--- Normalize ---*/
    modm =0 ; modl = 0;
    for (iDim =0 ; iDim < nDim; iDim++) {
      modm += m[iDim]*m[iDim];
      modl += l[iDim]*l[iDim];
    }
    modm = sqrt(modm);
    modl = sqrt(modl);
    for (iDim =0 ; iDim < nDim; iDim++) {
     l[iDim] = l[iDim]/modl;
     m[iDim] = m[iDim]/modm;
    }

  }
}

su2double CNumerics::GetRoe_Dissipation(const su2double Dissipation_i,
                                        const su2double Dissipation_j,
                                        const su2double Sensor_i,
                                        const su2double Sensor_j,
                                        const CConfig* config) const {

  /*--- Check for valid input ---*/

  unsigned short roe_low_diss = config->GetKind_RoeLowDiss();

  su2double Dissipation_ij = 0.0;

  assert((Dissipation_i >= 0) && (Dissipation_i <= 1));
  assert((Dissipation_j >= 0) && (Dissipation_j <= 1));
  if (roe_low_diss == FD_DUCROS || roe_low_diss == NTS_DUCROS) {
    assert((Sensor_i >= 0) && (Sensor_i <= 1));
    assert((Sensor_j >= 0) && (Sensor_j <= 1));
  }

  /*--- A minimum level of upwinding is used to enhance stability ---*/

  const su2double Min_Dissipation = 0.05;

  const su2double Mean_Dissipation = 0.5*(Dissipation_i + Dissipation_j);
  const su2double Mean_Sensor = 0.5*(Sensor_i + Sensor_j);

  if (roe_low_diss == FD || roe_low_diss == FD_DUCROS){

    Dissipation_ij = max(0.05,1.0 - (0.5 * (Dissipation_i + Dissipation_j)));

    if (roe_low_diss == FD_DUCROS){

      /*--- See Jonhsen et al. JCP 229 (2010) pag. 1234 ---*/

      su2double Ducros_ij;

      if (0.5*(Sensor_i + Sensor_j) > 0.65)
        Ducros_ij = 1.0;
      else
        Ducros_ij = 0.05;

      Dissipation_ij = max(Ducros_ij, Dissipation_ij);
    }

  } else if (roe_low_diss == NTS) {

    Dissipation_ij = max(Min_Dissipation, Mean_Dissipation);

  } else if (roe_low_diss == NTS_DUCROS) {

    /*--- See Xiao et al. INT J HEAT FLUID FL 51 (2015) pag. 141
     * https://doi.org/10.1016/j.ijheatfluidflow.2014.10.007 ---*/

    const su2double phi1 = Mean_Sensor;
    const su2double phi2 = Mean_Dissipation;

    Dissipation_ij = max(Min_Dissipation, phi1 + phi2 - (phi1*phi2));

  } else {

    SU2_MPI::Error("Unrecognized upwind/central blending scheme!",
                   CURRENT_FUNCTION);

  }
  return Dissipation_ij;
}
