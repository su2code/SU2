/*!
 * \file CVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../include/variables/CVariable.hpp"


CVariable::CVariable(Idx_t npoint, Idx_t nvar, CConfig *config) {

  /*--- Initialize the number of solution variables. This version
   of the constructor will be used primarily for converting the
   restart files into solution files (SU2_SOL). ---*/
  nPoint = npoint;
  nVar = nvar;

  /*--- Allocate the solution array - here it is also possible
   to allocate some extra flow variables that do not participate
   in the simulation ---*/
  Solution.resize(nPoint,nVar) = 0.0;

}

CVariable::CVariable(Idx_t npoint, Idx_t ndim, Idx_t nvar, CConfig *config) {

  /*--- Initializate the number of dimension and number of variables ---*/
  nPoint = npoint;
  nDim = ndim;
  nVar = nvar;

  /*--- Allocate solution, solution old, residual and gradient
   which is common for all the problems, here it is also possible
   to allocate some extra flow variables that do not participate
   in the simulation ---*/
  Solution.resize(nPoint,nVar) = 0.0;

  Solution_Old.resize(nPoint,nVar);

  Gradient.resize(nPoint,nVar,nDim,0.0);

  if (config->GetUnsteady_Simulation() != NO) {
    Solution_time_n.resize(nPoint,nVar);
    Solution_time_n1.resize(nPoint,nVar);
  }
  else if (config->GetDynamic_Analysis() == DYNAMIC) {
    Solution_time_n.resize(nPoint,nVar) = 0.0;
  }

	if (config->GetFSI_Simulation() && config->GetDiscrete_Adjoint()) {
	  Solution_Adj_Old.resize(nPoint,nVar);
	}

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    Gradient.resize(nPoint,nDim,nDim,0.0);
  }

}
