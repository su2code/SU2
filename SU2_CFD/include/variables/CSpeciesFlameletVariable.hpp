/*!
 * \file CSpeciesFlameletVariable.hpp
 * \brief Base class for defining the variables of the flamelet transport model.
 * \author D. Mayer, T. Economon, N. Beishuizen
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

#include "CSpeciesVariable.hpp"

/*!
 * \class CSpeciesFlameletVariable
 * \brief Base class for defining the variables of the flamelet model.
 */
class CSpeciesFlameletVariable final : public CSpeciesVariable {
 protected:
  MatrixType source_scalar; /*!< \brief Vector of the source terms from the lookup table for each scalar equation */
  MatrixType lookup_scalar; /*!< \brief Vector of the source terms from the lookup table for each scalar equation */
  MatrixType source_pd;     /*!< \brief PD closure source terms S_{phi_k} (Eq. 16) per control variable, for visualization. */
  MatrixType pd_flux_coeff; /*!< \brief Eq. (14) preferential diffusion flux coefficients (SOURCE_TERM method): per
                                 control variable k, columns [k*pd_terms_per_cv + s] hold the molecular coefficients
                                 D_{phi_k,s} for major species s, and column [k*pd_terms_per_cv + n_major_species]
                                 holds the thermal (Soret) coefficient D^T_{phi_k}. */
  unsigned short pd_terms_per_cv = 0; /*!< \brief Stride of pd_flux_coeff: one molecular coefficient per configured
                                           major species plus the single thermal coefficient. */
  su2vector<unsigned short> table_misses; /*!< \brief Vector of lookup table misses. */
  MatrixType source_cons_jac; /*!< \brief Consumption-rate Jacobian dS_aux_i/dY_aux_i = source_cons_i, one column per user scalar. */
  su2activevector hull_miss_dcv1_; /*!< \brief Signed CV1 deviation (query minus nearest hull node) at worst-miss Z level. */
  su2activevector hull_miss_dcv2_; /*!< \brief Signed CV2 deviation (query minus nearest hull node) at worst-miss Z level. */
  su2activevector z_level_dist_;   /*!< \brief Distance to nearest table Z level in physical Z units. */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] species_inf - species variable values (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSpeciesFlameletVariable(const su2double* species_inf, unsigned long npoint, unsigned long ndim, unsigned long nvar,
                           const CConfig* config);

  /*!
   * \brief Set the value of the transported scalar source term.
   * \param[in] iPoint - the location where the value has to be set.
   * \param[in] val_lookup_scalar - the value of the scalar to set.
   * \param[in] val_ivar - Eqn. index to the transport equation.
   */
  inline void SetLookupScalar(unsigned long iPoint, su2double val_lookup_scalar, unsigned short val_ivar) override {
    lookup_scalar(iPoint, val_ivar) = val_lookup_scalar;
  }

  /*!

   * \brief Set a source term for the specie transport equation.
   * \param[in] iPoint - Node index.
   * \param[in] val_ivar - Species index.
   * \param[in] val_source - Source term value.
  */
  inline void SetScalarSource(unsigned long iPoint, unsigned short val_ivar, su2double val_source) override {
    source_scalar(iPoint, val_ivar) = val_source;
  }

  /*!
   * \brief Get the value of the transported scalars source term.
   * \return Pointer to the transported scalars source term.
   */
  inline const su2double* GetScalarSources(unsigned long iPoint) const override { return source_scalar[iPoint]; }

  /*!
   * \brief Get the value of the looked up table based on the transported scalar.
   * \return Pointer to the transported scalars source term.
   */
  inline const su2double* GetScalarLookups(unsigned long iPoint) const override { return lookup_scalar[iPoint]; }

  /*!
   * \brief Store the PD closure source term S_{phi_k} for control variable iCV (for visualization).
   */
  inline void SetScalarSourcePD(unsigned long iPoint, unsigned short iCV, su2double val) {
    source_pd(iPoint, iCV) = val;
  }

  /*!
   * \brief Get the PD closure source terms S_{phi_k} for all control variables at iPoint.
   */
  inline const su2double* GetScalarSourcesPD(unsigned long iPoint) const override { return source_pd[iPoint]; }

  /*!
   * \brief Store an Eq. (14) preferential diffusion flux coefficient (SOURCE_TERM method).
   * \param[in] iPoint - Node index.
   * \param[in] iCV - Control variable index.
   * \param[in] iTerm - Major species index (molecular), or n_major_species for the thermal coefficient.
   * \param[in] val - Coefficient value from the manifold.
   */
  inline void SetPDFluxCoeff(unsigned long iPoint, unsigned short iCV, unsigned short iTerm, su2double val) {
    pd_flux_coeff(iPoint, iCV * pd_terms_per_cv + iTerm) = val;
  }

  /*!
   * \brief Get an Eq. (14) preferential diffusion flux coefficient (SOURCE_TERM method).
   * \param[in] iPoint - Node index.
   * \param[in] iCV - Control variable index.
   * \param[in] iTerm - Major species index (molecular), or n_major_species for the thermal coefficient.
   */
  inline su2double GetPDFluxCoeff(unsigned long iPoint, unsigned short iCV, unsigned short iTerm) const override {
    return pd_flux_coeff(iPoint, iCV * pd_terms_per_cv + iTerm);
  }

  /*!
   * \brief The whole Eq. (14) coefficient matrix, for the edge flux to gather from. Column
   *        iCV * GetPDTermsPerCV() + iTerm, with iTerm == n_major_species the thermal coefficient.
   */
  inline const MatrixType& GetPDFluxCoeffs() const { return pd_flux_coeff; }

  /*!
   * \brief Column stride of GetPDFluxCoeffs(): one coefficient per major species plus the thermal one.
   */
  inline unsigned short GetPDTermsPerCV() const { return pd_terms_per_cv; }

  inline void SetTableMisses(unsigned long iPoint, unsigned short misses) override { table_misses[iPoint] = misses; }

  inline unsigned short GetTableMisses(unsigned long iPoint) const override { return table_misses[iPoint]; }

  inline void SetHullMissDevCV1(unsigned long iPoint, su2double val) override { hull_miss_dcv1_(iPoint) = val; }
  inline su2double GetHullMissDevCV1(unsigned long iPoint) const override { return hull_miss_dcv1_(iPoint); }
  inline void SetHullMissDevCV2(unsigned long iPoint, su2double val) override { hull_miss_dcv2_(iPoint) = val; }
  inline su2double GetHullMissDevCV2(unsigned long iPoint) const override { return hull_miss_dcv2_(iPoint); }
  inline void SetZLevelDist(unsigned long iPoint, su2double val) override { z_level_dist_(iPoint) = val; }
  inline su2double GetZLevelDist(unsigned long iPoint) const override { return z_level_dist_(iPoint); }

  /*!
   * \brief Store the consumption-rate Jacobian dS_aux_i/dY_aux_i = source_cons_i for user scalar i_aux.
   */
  inline void SetAuxSourceCons(unsigned long iPoint, unsigned long i_aux, su2double val) {
    source_cons_jac(iPoint, i_aux) = val;
  }

  /*!
   * \brief Get the consumption-rate Jacobian entry for user scalar i_aux at iPoint.
   */
  inline su2double GetAuxSourceCons(unsigned long iPoint, unsigned long i_aux) const {
    return source_cons_jac(iPoint, i_aux);
  }
};
