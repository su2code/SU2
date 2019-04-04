/*!
 * \file linear_solvers_structure_b.cpp
 * \brief Routines for the linear solver used in the reverse sweep of AD.
 * \author T. Albring
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

#include "../include/linear_solvers_structure_b.hpp"
#include "../include/linear_solvers_structure.hpp"
#include "../include/vector_structure.hpp"
#include "../include/matrix_structure.hpp"

#ifdef CODI_REVERSE_TYPE
void CSysSolve_b::Solve_b(AD::Tape* tape, AD::CheckpointHandler* data) {
  
  /*--- Extract data from the checkpoint handler ---*/

  su2double::GradientData *LinSysRes_Indices;
  su2double::GradientData *LinSysSol_Indices;
#if CODI_PRIMAL_INDEX_TAPE
  su2double::Real *oldValues;
#endif

  data->getData(LinSysRes_Indices);
  data->getData(LinSysSol_Indices);
#if CODI_PRIMAL_INDEX_TAPE
  data->getData(oldValues);
#endif

  unsigned long nBlk = 0, nVar = 0, nBlkDomain = 0, size = 0, i = 0;

  data->getData(size);
  data->getData(nBlk);
  data->getData(nVar);
  data->getData(nBlkDomain);

  CSysMatrix* Jacobian = NULL;
  data->getData(Jacobian);

  CGeometry* geometry  = NULL;
  data->getData(geometry);

  CConfig* config      = NULL;
  data->getData(config);

  CSysSolve* solver;
  data->getData(solver);

  /*--- Initialize the right-hand side with the gradient of the solution of the primal linear system ---*/

  CSysVector LinSysRes_b(nBlk, nBlkDomain, nVar, 0.0);
  CSysVector LinSysSol_b(nBlk, nBlkDomain, nVar, 0.0);

  for (i = 0; i < size; i ++) {
    su2double::GradientData& index = LinSysSol_Indices[i];
    LinSysRes_b[i] = AD::globalTape.getGradient(index);
    LinSysSol_b[i] = 0.0;
    AD::globalTape.gradient(index) = 0.0;
  }

  /*--- Solve the system ---*/

  solver->Solve_b(*Jacobian, LinSysRes_b, LinSysSol_b, geometry, config);

  /*--- Update the gradients of the right-hand side of the primal linear system ---*/

  for (i = 0; i < size; i ++) {
    su2double::GradientData& index = LinSysRes_Indices[i];
    AD::globalTape.gradient(index) += SU2_TYPE::GetValue(LinSysSol_b[i]);
  }

#if CODI_PRIMAL_INDEX_TAPE
  /*--- Set the old values that have been overwritten ---*/
  for (i = 0; i < size; i ++) {
    AD::globalTape.setExternalValueChange(LinSysSol_Indices[i], oldValues[i]);
  }
#endif
}

void CSysSolve_b::Delete_b(AD::Tape* tape, AD::CheckpointHandler* data) {

  su2double::GradientData *LinSysRes_Indices = NULL;
  su2double::GradientData *LinSysSol_Indices = NULL;
#if CODI_PRIMAL_INDEX_TAPE
  su2double::Real *oldValues;
#endif

  data->getData(LinSysRes_Indices);
  data->getData(LinSysSol_Indices);
#if CODI_PRIMAL_INDEX_TAPE
  data->getData(oldValues);
#endif

  delete [] LinSysRes_Indices;
  delete [] LinSysSol_Indices;
#if CODI_PRIMAL_INDEX_TAPE
  delete [] oldValues;
#endif

  unsigned long nBlk, nVar, nBlkDomain, size;

  data->getData(size);
  data->getData(nBlk);
  data->getData(nVar);
  data->getData(nBlkDomain);

  CSysMatrix* Jacobian;
  data->getData(Jacobian);

  CGeometry* geometry;
  data->getData(geometry);

  CConfig* config;
  data->getData(config);

  CSysSolve* solver;
  data->getData(solver);
}
#endif
