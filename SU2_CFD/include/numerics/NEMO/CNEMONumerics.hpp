/*!
 * \file CNEMONumerics.hpp
 * \brief Base class template NEMO numerics.
 * \author C. Garbacz, W. Maier, S. R. Copeland
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

#include "../CNumerics.hpp"
#include "../../fluid/CNEMOGas.hpp"
#include "../../fluid/CMutationTCLib.hpp"
#include "../../fluid/CSU2TCLib.hpp"
#include "../../../../Common/include/toolboxes/geometry_toolbox.hpp"

/*!
 * \class CNEMONumerics
 * \brief Base class template NEMO numerics.
 * \author C. Garbacz., W. Maier
 */
class CNEMONumerics : public CNumerics {
public:
  bool implicit, ionization;
  su2double *rhos_i, *rhos_j;
  su2double Velocity_i[MAXNDIM] = {0.0}, Velocity_j[MAXNDIM] = {0.0};
  su2double e_ve_i, e_ve_j;
  su2double rhoCvtr_i, rhoCvtr_j;
  su2double rhoCvve_i, rhoCvve_j;
  su2double ProjVelocity_i, ProjVelocity_j;
  unsigned short nPrimVar, nPrimVarGrad;

  su2double* Flux = nullptr;            /*!< \brief The flux / residual across the edge. */
  su2double** Jacobian_i = nullptr;
  su2double** Jacobian_j = nullptr;

  unsigned short nSpecies, nHeavy, nEl; /*!< \brief Number of species present in plasma */

  /*--- Graidents w.r.t. conservative variables. ---*/
  su2double *dPdU_i, *dPdU_j;
  su2double *dTdU_i, *dTdU_j;
  su2double *dTvedU_i, *dTvedU_j;
  su2double Gamma_i, Gamma_j;

  su2double *eve_i, *eve_j, *Cvve_i, *Cvve_j;

  unsigned short RHOS_INDEX, T_INDEX, TVE_INDEX, VEL_INDEX, P_INDEX,
  RHO_INDEX, H_INDEX, A_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX,
  LAM_VISC_INDEX, EDDY_VISC_INDEX;

  CNEMOGas *fluidmodel;

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive variables of the problem.
   * \param[in] val_nPrimVarGrad - Number of primitive grad. variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CNEMONumerics(unsigned short val_nDim, unsigned short val_nVar,
                unsigned short val_nPrimVar,
                unsigned short val_nPrimVarGrad,
                const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CNEMONumerics(void);

  /*!
   * \Overload
   * \brief Compute the projected inviscid flux vector.
   * \param[in] val_U - Pointer to the conserved variables.
   * \param[in] val_V - Pointer to the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Flux - Pointer to the projected flux.
   */
  void GetInviscidProjFlux(const su2double *val_U, const su2double *val_V,
                           const su2double *val_normal, su2double *val_Proj_Flux);

  /*!
   * \overload
   * \brief Compute the projection of the inviscid Jacobian matrices for the two-temperature model.
   * \param[in] val_U - Vector conserved variables.
   * \param[in] val_V - Vector of primitive variables.
   * \param[in] val_dPdU - Vector of partial derivatives of pressure w.r.t. conserved vars.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidProjJac(const su2double *val_U, const su2double *val_V, const su2double *val_dPdU,
                          const su2double *val_normal, const su2double val_scale,
                          su2double **val_Proj_Jac_Tensor);

  /*!
   * \brief Compute the projection of the viscous fluxes into a direction.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of Primitive Variables.
   * \param[in] val_eve - Virbational-Electronical Energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_diffusioncoeff - Disffusion Coefficient.
   * \param[in] val_lam_viscosity - Laminar Viscosity
   * \param[in] val_eddy_viscosity - Eddy Viscosity
   * \param[in] val_thermal_conductivity - Thermal conductivity.
   * \param[in] val_thermal_conductivity_ve - Thermal conductivity of Vibe-Elec modes.
   * \param[in] config - Definition of the particular problem.
   */
  void GetViscousProjFlux(const su2double *val_primvar,
                          const su2double* const* val_gradprimvar,
                          su2double *val_eve,
                          const su2double *val_normal,
                          const su2double *val_diffusioncoeff,
                          su2double val_lam_viscosity,
                          su2double val_eddy_viscosity,
                          su2double val_therm_conductivity,
                          su2double val_therm_conductivity_ve,
                          const CConfig *config);

  /*!
   * \brief Staging function to compute viscous Jacobians.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_Mean_SecVar - Mean value of the secondary variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_thermal_conductivity - Value of the thermal conductivity.
   * \param[in] val_heat_capacity_cp - Value of the specific heat at constant pressure.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousProjJacs(const su2double *val_Mean_PrimVar,
                          su2double *val_Mean_Eve, const su2double *val_Mean_Cvve,
                          const su2double *val_diffusion_coeff, su2double val_laminar_viscosity,
                          su2double val_eddy_viscosity, su2double val_thermal_conductivity,
                          su2double val_thermal_conductivity_ve,
                          su2double val_dist_ij,
                          const su2double *val_normal,
                          su2double val_dS, const su2double *val_Fv,
                          su2double **val_Jac_i, su2double **val_Jac_j,
                          const CConfig *config);

  /*!
   * \brief TSL-Approximation of Viscous NS Jacobians for arbitrary equations of state.
   */
  template <unsigned short NVAR, unsigned short NSPECIES>
  void ComputeViscousJacs_impl(const su2double *val_Mean_PrimVar,
                               su2double *val_Mean_Eve,
                               const su2double *val_Mean_Cvve,
                               const su2double *val_diffusion_coeff,
                               su2double val_laminar_viscosity,
                               su2double val_eddy_viscosity,
                               su2double val_thermal_conductivity,
                               su2double val_thermal_conductivity_ve,
                               su2double val_dist_ij, const su2double *val_normal,
                               su2double val_dS, const su2double *val_Fv,
                               su2double **val_Jac_i, su2double **val_Jac_j,
                               const CConfig *config){

    /*--- Initialize matrices ---*/
    Matrix<NVAR, NVAR> dFdVi, dFdVj, dVdUi, dVdUj;

    /*--- Need to avoid size of 1 otherwise the type becomes a vector (without [i][j]). ---*/
    constexpr unsigned short nSpeciesDummy = (NSPECIES == 1) ? 2 : NSPECIES;
    Matrix<nSpeciesDummy, nSpeciesDummy> dJdr_i, dJdr_j;
    Matrix<1,NSPECIES> Ys, Ys_i, Ys_j;

    /*--- Play tricks on the compiler, in static mode use NVAR from the template, in dynamic mode
          use nVar from the class or from the arguments. ---*/
    const auto nVar = (NVAR != DynamicSize) ? NVAR : this->nVar;
    const auto nSpecies = (NSPECIES != DynamicSize) ? NSPECIES : this->nSpecies;

    /*--- Deduce nDim from the other variables to allow the compiler to also unroll the loops over nDim. ---*/
    const auto nDim = nVar - nSpecies - 2;

    // Allocate and initialize, for the static case the compiler optimizes this away.
    dFdVi.resize(nVar, nVar) = su2double(0.0);
    dFdVj.resize(nVar, nVar) = su2double(0.0);
    dVdUi.resize(nVar, nVar) = su2double(0.0);
    dVdUj.resize(nVar, nVar) = su2double(0.0);
    if (NSPECIES == DynamicSize) {
      dJdr_i.resize(nSpecies,nSpecies) = su2double(0.0);
      dJdr_j.resize(nSpecies,nSpecies) = su2double(0.0);
      Ys.resize(1,nSpecies) = su2double(0.0);
      Ys_i.resize(1,nSpecies) = su2double(0.0);
      Ys_j.resize(1,nSpecies) = su2double(0.0);
    }

    /*--- Calculate preliminary geometric quantities ---*/
    su2double dij = val_dist_ij;
    su2double theta = GeometryToolbox::SquaredNorm(nDim, val_normal);

    /*--- Rename for convenience ---*/
    su2double rho_i = V_i[RHO_INDEX];
    su2double rho_j = V_j[RHO_INDEX];
    su2double T   = val_Mean_PrimVar[T_INDEX];
    su2double Tve = val_Mean_PrimVar[TVE_INDEX];
    su2double mu  = val_laminar_viscosity;
    su2double ktr = val_thermal_conductivity;
    su2double kve = val_thermal_conductivity_ve;
    su2double RuSI= UNIVERSAL_GAS_CONSTANT;
    su2double Ru  = 1000.0*RuSI;
    const auto& Ds  = val_diffusion_coeff;
    const auto& hs = fluidmodel->ComputeSpeciesEnthalpy(T, Tve, val_Mean_Eve);
    const auto& Cvtr = fluidmodel->GetSpeciesCvTraRot();
    const auto& Ms = fluidmodel->GetSpeciesMolarMass();

    /*--- Calculate useful diffusion parameters ---*/
    // Mass fraction and summation of diffusion fluxes
    su2double sumY = 0.0;
    su2double sumY_i = 0.0;
    su2double sumY_j = 0.0;
    for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Ys[iSpecies]   = val_Mean_PrimVar[RHOS_INDEX+iSpecies];
      Ys_i[iSpecies] = V_i[RHOS_INDEX+iSpecies]/V_i[RHO_INDEX];
      Ys_j[iSpecies] = V_j[RHOS_INDEX+iSpecies]/V_j[RHO_INDEX];

      sumY_i += Ds[iSpecies]*theta/dij*Ys_i[iSpecies];
      sumY_j += Ds[iSpecies]*theta/dij*Ys_j[iSpecies];
      sumY   += Ds[iSpecies]*theta/dij*(Ys_j[iSpecies]-Ys_i[iSpecies]);
    }

    su2double vel[MAXNDIM];
    for (auto iDim = 0; iDim < nDim; iDim++)
      vel[iDim] = val_Mean_PrimVar[VEL_INDEX+iDim];

    for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (auto jSpecies  = 0; jSpecies < nSpecies; jSpecies++) {

        // first term
        dJdr_j[iSpecies][jSpecies] +=  0.5*(Ds[iSpecies]*theta/dij *
                                           (Ys_j[iSpecies]*rho_i/rho_j +
                                            Ys_i[iSpecies]));
        dJdr_i[iSpecies][jSpecies] += -0.5*(Ds[iSpecies]*theta/dij *
                                           (Ys_j[iSpecies] +
                                            Ys_i[iSpecies]*rho_j/rho_i));

        // second term
        dJdr_j[iSpecies][jSpecies] +=
              0.25*(Ys_i[iSpecies] - rho_i/rho_j*Ys_j[iSpecies])*sumY
            + 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*Ds[jSpecies]*theta/(dij*rho_j)
            - 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*sumY_j/rho_j;

        dJdr_i[iSpecies][jSpecies] +=
              0.25*(-rho_j/rho_i*Ys_i[iSpecies]+Ys_j[iSpecies])*sumY
            - 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*Ds[jSpecies]*theta/(dij*rho_i)
            + 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*sumY_i/rho_i;
      }

      // first term
      dJdr_j[iSpecies][iSpecies] += -0.5*Ds[iSpecies]*theta/dij*(1+rho_i/rho_j);
      dJdr_i[iSpecies][iSpecies] +=  0.5*Ds[iSpecies]*theta/dij*(1+rho_j/rho_i);

      // second term
      dJdr_j[iSpecies][iSpecies] += 0.25*(1.0+rho_i/rho_j)*sumY;
      dJdr_i[iSpecies][iSpecies] += 0.25*(1.0+rho_j/rho_i)*sumY;
    }

    /*--- Calculate transformation matrix ---*/
    for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      dVdUi[iSpecies][iSpecies] = 1.0;
      dVdUj[iSpecies][iSpecies] = 1.0;
    }
    for (auto iDim = 0; iDim < nDim; iDim++) {
      for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        dVdUi[nSpecies+iDim][iSpecies] = -V_i[VEL_INDEX+iDim]/V_i[RHO_INDEX];
        dVdUj[nSpecies+iDim][iSpecies] = -V_j[VEL_INDEX+iDim]/V_j[RHO_INDEX];
      }
      dVdUi[nSpecies+iDim][nSpecies+iDim] = 1.0/V_i[RHO_INDEX];
      dVdUj[nSpecies+iDim][nSpecies+iDim] = 1.0/V_j[RHO_INDEX];
    }
    for (auto iVar = 0; iVar < nVar; iVar++) {
      dVdUi[nSpecies+nDim][iVar]   = dTdU_i[iVar];
      dVdUj[nSpecies+nDim][iVar]   = dTdU_j[iVar];
      dVdUi[nSpecies+nDim+1][iVar] = dTvedU_i[iVar];
      dVdUj[nSpecies+nDim+1][iVar] = dTvedU_j[iVar];
    }

    if (nDim == 2) {

      /*--- Geometry parameters ---*/
      su2double thetax = theta + val_normal[0]*val_normal[0]/3.0;
      su2double thetay = theta + val_normal[1]*val_normal[1]/3.0;
      su2double etaz   = val_normal[0]*val_normal[1]/3.0;
      su2double pix    = mu/dij * (thetax*vel[0] + etaz*vel[1]  );
      su2double piy    = mu/dij * (etaz*vel[0]   + thetay*vel[1]);

      /*--- Populate primitive Jacobian ---*/

      // X-momentum
      dFdVj[nSpecies][nSpecies]     = mu*thetax/dij*val_dS;
      dFdVj[nSpecies][nSpecies+1]   = mu*etaz/dij*val_dS;

      // Y-momentum
      dFdVj[nSpecies+1][nSpecies]   = mu*etaz/dij*val_dS;
      dFdVj[nSpecies+1][nSpecies+1] = mu*thetay/dij*val_dS;

      // Energy
      dFdVj[nSpecies+2][nSpecies]   = pix*val_dS;
      dFdVj[nSpecies+2][nSpecies+1] = piy*val_dS;
      dFdVj[nSpecies+2][nSpecies+2] = ktr*theta/dij*val_dS;
      dFdVj[nSpecies+2][nSpecies+3] = kve*theta/dij*val_dS;

      // Vib-el Energy
      dFdVj[nSpecies+3][nSpecies+3] = kve*theta/dij*val_dS;

      for (auto iVar = 0; iVar < nVar; iVar++)
        for (auto jVar = 0; jVar < nVar; jVar++)
          dFdVi[iVar][jVar] = -dFdVj[iVar][jVar];

      // Common terms
      dFdVi[nSpecies+2][nSpecies]   += 0.5*val_Fv[nSpecies];
      dFdVj[nSpecies+2][nSpecies]   += 0.5*val_Fv[nSpecies];
      dFdVi[nSpecies+2][nSpecies+1] += 0.5*val_Fv[nSpecies+1];
      dFdVj[nSpecies+2][nSpecies+1] += 0.5*val_Fv[nSpecies+1];
      for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        dFdVi[nSpecies+2][nSpecies+2] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                               Cvtr[iSpecies]   );
        dFdVj[nSpecies+2][nSpecies+2] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                               Cvtr[iSpecies]   );
        dFdVi[nSpecies+2][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
        dFdVj[nSpecies+2][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
        dFdVi[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
        dFdVj[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      }

      // Unique terms
      for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        for (auto jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
          dFdVj[iSpecies][jSpecies]   += -dJdr_j[iSpecies][jSpecies]*val_dS;
          dFdVi[iSpecies][jSpecies]   += -dJdr_i[iSpecies][jSpecies]*val_dS;
          dFdVj[nSpecies+2][iSpecies] += -dJdr_j[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
          dFdVi[nSpecies+2][iSpecies] += -dJdr_i[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
          dFdVj[nSpecies+3][iSpecies] += -dJdr_j[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
          dFdVi[nSpecies+3][iSpecies] += -dJdr_i[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
        }
      }

    } //nDim == 2
    else {

      /*--- Geometry parameters ---*/
      su2double thetax = theta + val_normal[0]*val_normal[0]/3.0;
      su2double thetay = theta + val_normal[1]*val_normal[1]/3.0;
      su2double thetaz = theta + val_normal[2]*val_normal[2]/3.0;
      su2double etax   = val_normal[1]*val_normal[2]/3.0;
      su2double etay   = val_normal[0]*val_normal[2]/3.0;
      su2double etaz   = val_normal[0]*val_normal[1]/3.0;
      su2double pix    = mu/dij * (thetax*vel[0] + etaz*vel[1]   + etay*vel[2]  );
      su2double piy    = mu/dij * (etaz*vel[0]   + thetay*vel[1] + etax*vel[2]  );
      su2double piz    = mu/dij * (etay*vel[0]   + etax*vel[1]   + thetaz*vel[2]);

      /*--- Populate primitive Jacobian ---*/

      // X-momentum
      dFdVj[nSpecies][nSpecies]     = mu*thetax/dij*val_dS;
      dFdVj[nSpecies][nSpecies+1]   = mu*etaz/dij*val_dS;
      dFdVj[nSpecies][nSpecies+2]   = mu*etay/dij*val_dS;

      // Y-momentum
      dFdVj[nSpecies+1][nSpecies]   = mu*etaz/dij*val_dS;
      dFdVj[nSpecies+1][nSpecies+1] = mu*thetay/dij*val_dS;
      dFdVj[nSpecies+1][nSpecies+2] = mu*etax/dij*val_dS;

      // Z-momentum
      dFdVj[nSpecies+2][nSpecies]   = mu*etay/dij*val_dS;
      dFdVj[nSpecies+2][nSpecies+1] = mu*etax/dij*val_dS;
      dFdVj[nSpecies+2][nSpecies+2] = mu*thetaz/dij*val_dS;

      // Energy
      dFdVj[nSpecies+3][nSpecies]   = pix*val_dS;
      dFdVj[nSpecies+3][nSpecies+1] = piy*val_dS;
      dFdVj[nSpecies+3][nSpecies+2] = piz*val_dS;
      dFdVj[nSpecies+3][nSpecies+3] = ktr*theta/dij*val_dS;
      dFdVj[nSpecies+3][nSpecies+4] = kve*theta/dij*val_dS;

      // Vib.-el energy
      dFdVj[nSpecies+4][nSpecies+4] = kve*theta/dij*val_dS;

      for (auto iVar = 0; iVar < nVar; iVar++)
        for (auto jVar = 0; jVar < nVar; jVar++)
          dFdVi[iVar][jVar] = -dFdVj[iVar][jVar];

      // Common terms
      for (auto iDim = 0; iDim < nDim; iDim++) {
        dFdVi[nSpecies+3][nSpecies+iDim]   += 0.5*val_Fv[nSpecies+iDim];
        dFdVj[nSpecies+3][nSpecies+iDim]   += 0.5*val_Fv[nSpecies+iDim];
      }
      for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        dFdVi[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                               Cvtr[iSpecies]   );
        dFdVj[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                               Cvtr[iSpecies]   );
        dFdVi[nSpecies+3][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
        dFdVj[nSpecies+3][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
        dFdVi[nSpecies+4][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
        dFdVj[nSpecies+4][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      }

      // Unique terms
      for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        for (auto jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
          dFdVj[iSpecies][jSpecies]   += -dJdr_j[iSpecies][jSpecies]*val_dS;
          dFdVi[iSpecies][jSpecies]   += -dJdr_i[iSpecies][jSpecies]*val_dS;
          dFdVj[nSpecies+3][iSpecies] += -dJdr_j[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
          dFdVi[nSpecies+3][iSpecies] += -dJdr_i[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
          dFdVj[nSpecies+4][iSpecies] += -dJdr_j[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
          dFdVi[nSpecies+4][iSpecies] += -dJdr_i[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
        }
      }

    } // nDim == 3

    /*--- dFv/dUij = dFv/dVij * dVij/dUij ---*/
    for (auto iVar = 0; iVar < nVar; iVar++)
      for (auto jVar = 0; jVar < nVar; jVar++)
        for (auto kVar = 0; kVar < nVar; kVar++) {
          val_Jac_i[iVar][jVar] += dFdVi[iVar][kVar]*dVdUi[kVar][jVar];
          val_Jac_j[iVar][jVar] += dFdVj[iVar][kVar]*dVdUj[kVar][jVar];
        }
  }
  /*--- Template matrix creation for code optimization ---*/
  template <int nRows, int nCols>
  using Matrix = C2DContainer<unsigned short, su2double, StorageType::RowMajor, 8, nRows, nCols>;

  /*!
   * \overload
   * \brief Computation of the matrix P, this matrix diagonalizes the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] U - Vector of conserved variables (really only need rhoEve)
   * \param[in] V - Vector of primitive variables
   * \param[in] val_dPdU - Vector of derivatives of pressure w.r.t. conserved vars.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] l - Tangential vector to face.
   * \param[in] m - Tangential vector to face (mutually orthogonal to val_normal & l).
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix(const su2double *U, const su2double *V, const su2double *val_dPdU,
                  const su2double *val_normal, const su2double *l, const su2double *m,
                  su2double **val_p_tensor) const;

  /*!
   * \overload
   * \brief Computation of the matrix P^{-1}, this matrix diagonalizes the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] U - Vector of conserved variables.
   * \param[in] V - Vector of primitive variables.
   * \param[in] val_dPdU - Vector of derivatives of pressure w.r.t. conserved variables
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] l - Tangential vector to face.
   * \param[in] m - Tangential vector to face (mutually orthogonal to val_normal & l).
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix_inv(const su2double *U, const su2double *V, const su2double *val_dPdU,
                      const su2double *val_normal, const su2double *l, const su2double *m,
                      su2double **val_invp_tensor) const;


  /*!
   * \brief Set the pressure derivatives.
   * \param[in] val_dPdU_i - pressure derivatives at i.
   * \param[in] val_dPdU_j - pressure derivatives at j.
   */
  inline void SetdPdU(su2double *val_dPdU_i, su2double *val_dPdU_j)       final { dPdU_i = val_dPdU_i; dPdU_j = val_dPdU_j; }

  /*!
   * \brief Set the temperature derivatives.
   * \param[in] val_dTdU_i - temperature derivatives at i.
   * \param[in] val_dTdU_j - temperature derivatives at j.
   */
  inline void SetdTdU(su2double *val_dTdU_i, su2double *val_dTdU_j)       final { dTdU_i = val_dTdU_i; dTdU_j = val_dTdU_j; }

  /*!
   * \brief Set the vib-el temperature derivatives.
   * \param[in] val_dTvedU_i - t_ve derivatives at i.
   * \param[in] val_dTvedU_j - t_ve derivatives at j.
   */
  inline void SetdTvedU(su2double *val_dTvedU_i, su2double *val_dTvedU_j) final { dTvedU_i = val_dTvedU_i; dTvedU_j = val_dTvedU_j; }

  /*!
   * \brief Set the vib-el energy.
   * \param[in] val_Eve_i - vib-el energy at i.
   * \param[in] val_Eve_j - vib-el energy at j.
   */
  inline void SetEve(su2double *val_Eve_i, su2double *val_Eve_j)          final {eve_i = val_Eve_i; eve_j = val_Eve_j; }

  /*!
   * \brief Set the Cvve.
   * \param[in] val_Cvve_i - cvve at i.
   * \param[in] val_Cvve_j - cvve at j.
   */
  inline void SetCvve(su2double *val_Cvve_i, su2double *val_Cvve_j)       final {Cvve_i = val_Cvve_i; Cvve_j = val_Cvve_j; }

  /*!
   * \brief Set the ratio of specific heats.
   * \param[in] val_Gamma_i - Gamma at i.
   * \param[in] val_Gamma_j - Gamma at j.
   */
  inline void SetGamma(su2double val_Gamma_i, su2double val_Gamma_j)      final {Gamma_i = val_Gamma_i; Gamma_j = val_Gamma_j; }

};
