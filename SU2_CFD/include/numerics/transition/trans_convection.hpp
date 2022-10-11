/*!
 * \file trams_convection.hpp
 * \brief Declarations of numerics classes for discretization of
 *        convective fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.4.0 "Blackbird"
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

#pragma once

#include "../CNumerics.hpp"
#include "../scalar/scalar_convection.hpp"

/*!
 * \class CUpwScalar
 * \brief Template class for scalar upwind fluxes between nodes i and j.
 * \details This class serves as a template for the scalar upwinding residual
 *   classes.  The general structure of a scalar upwinding calculation is the
 *   same for many different  models, which leads to a lot of repeated code.
 *   By using the template design pattern, these sections of repeated code are
 *   moved to this shared base class, and the specifics of each model
 *   are implemented by derived classes.  In order to add a new residual
 *   calculation for a convection residual, extend this class and implement
 *   the pure virtual functions with model-specific behavior.
 * \ingroup ConvDiscr
 * \author C. Pederson, A. Bueno., and A. Campos.
 */


/*!
 * \class my_CUpwSca_TransLM
 * \brief Class for doing a scalar upwind solver for the Langtry-Menter gamma-Re_theta transition model.
 * \ingroup ConvDiscr
 * \author A. Campos.
 */
template <class FlowIndices>
class CUpwSca_TransLM final : public CUpwScalar<FlowIndices> {
private:

    using Base = CUpwScalar<FlowIndices>;
    using Base::nDim;
    using Base::V_i;
    using Base::V_j;
    using Base::a0;
    using Base::a1;
    using Base::Flux;
    using Base::Jacobian_i;
    using Base::Jacobian_j;
    using Base::ScalarVar_i;
    using Base::ScalarVar_j;
    using Base::implicit;
    using Base::idx;

    // salatate le parti di AD

    /*!
   * \brief Adds any extra variables to AD
     */
    void ExtraADPreaccIn() override {}


    /*!
     * \brief LM specific steps in the ComputeResidual method
     * \param[in] config - Definition of the particular problem.
     */
    void FinishResidualCalc(const CConfig* config) override{

//      if(config->dummyVar == 0)
//        cout << "CUpwSca_TransLM::FinishResidualCalc" << endl;

      Flux[0] = a0*V_i[idx.Density()]*ScalarVar_i[0]+a1*V_j[idx.Density()]*ScalarVar_j[0];
      Flux[1] = a0*V_i[idx.Density()]*ScalarVar_i[1]+a1*V_j[idx.Density()]*ScalarVar_j[1];

//      cout << "Flux[0] = " << Flux[0] << endl;
//      cout << "Flux[1] = " << Flux[1] << endl;

      if (implicit) {
        Jacobian_i[0][0] = a0;    Jacobian_i[0][1] = 0.0;
        Jacobian_i[1][0] = 0.0;   Jacobian_i[1][1] = a0;

        Jacobian_j[0][0] = a1;    Jacobian_j[0][1] = 0.0;
        Jacobian_j[1][0] = 0.0;   Jacobian_j[1][1] = a1;
      }

      // è molto copia e incolla rispetto al turb, se è corretto come approccio allora bisognerà ripensare anche il sources!

    }

public:
    /*!
     * \brief Constructor of the class.
     * \param[in] val_nDim - Number of dimensions of the problem.
     * \param[in] val_nVar - Number of variables of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    CUpwSca_TransLM(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
     : CUpwScalar<FlowIndices>(val_nDim, val_nVar, config) {}


};