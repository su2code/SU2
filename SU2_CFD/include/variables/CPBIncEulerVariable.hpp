/*!
 * \file CPBIncEulerVariable.hpp
 * \brief Class for defining the variables of the pressure based incompressible Euler solver.
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

#include <limits>
#include "CFlowVariable.hpp"

/*!
 * \class CPBIncEulerVariable
 * \brief Class for defining the variables of the pressure based incompressible Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CPBIncEulerVariable : public CFlowVariable {
public:
  static constexpr size_t MAXNVAR = 13;

  template <class IndexType>
  struct CIndices {
    const IndexType nDim;
    CIndices(IndexType ndim, IndexType) : nDim(ndim) {}
    inline IndexType NDim() const { return nDim; }
    inline IndexType NSpecies() const { return 0; }
    inline IndexType Pressure() const { return 0; }
    inline IndexType Velocity() const { return 1; }
    inline IndexType Temperature() const { return nDim+1; }
    inline IndexType Density() const { return nDim+2; }
    inline IndexType Enthalpy() const { return nDim+3; }
    inline IndexType Beta() const { return nDim+4; }
    inline IndexType SoundSpeed() const { return Beta(); }
    inline IndexType LaminarViscosity() const { return nDim+5; }
    inline IndexType EddyViscosity() const { return nDim+6; }
    inline IndexType ThermalConductivity() const { return nDim+7; }
    inline IndexType CpTotal() const { return nDim+8; }
    inline IndexType CvTotal() const { return nDim+9; }

    /*--- For compatible interface with NEMO. ---*/
    inline IndexType SpeciesDensities() const { return std::numeric_limits<IndexType>::max(); }
    inline IndexType Temperature_ve() const { return std::numeric_limits<IndexType>::max(); }
  };

protected:
  const CIndices<unsigned long> indices;

public:

  CPBIncEulerVariable(su2double pressure, const su2double *velocity, su2double enthalpy,
                    unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config);

  inline CMatrixView<const su2double> GetVelocityGradient(unsigned long iPoint) const final {
    return Gradient_Primitive(iPoint, indices.Velocity());
  }


};
