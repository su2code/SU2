/*!
 * \file dataype_structure.cpp
 * \brief Main subroutines for the datatype structures.
 * \author T. Albring
 * \version 4.0.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (francisco.palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/datatype_structure.hpp"

#if defined CODI_REVERSE_TYPE
namespace AD {
  int adjointVectorPosition = 0;

  std::vector<unsigned int> inputValues;

  codi::ChunkTape<double, int>& globalTape = codi::RealReverse::getGlobalTape();
}
#elif defined ADOLC_REVERSE_TYPE
namespace AD{
  /* --- Stores a copy of the input variables (since they might be overwritten) --- */

  std::vector<double> inputVariables;

  /* --- Stores the seeding vector for the adjoint computation (adjoints of output variables) --- */

  std::vector<double> seedVector;

  /* --- The current position in the adjoint vector --- */

  int adjointVector_Position = 0;

  /* --- Holds the adjoint values of the input variables after the adjoint computation --- */

  double* adjointVector = NULL;

  void ClearAdjoints(){
    if (adjointVector != NULL){
      delete [] adjointVector;
    }
    adjointVector = NULL;
    seedVector.clear();
  }

  void ComputeAdjoint(){
    size_t tape_stats[STAT_SIZE];

    /* --- Get information about the current number of inputs/outputs --- */

    tapestats(1, tape_stats);

    /* --- Create temporary arrays to hold the adjoints of the input/output variables --- */

    double* gradient = new double[tape_stats[0]];
    double* adjoint  = new double[tape_stats[1]];

    /* --- Initialize the adjoint values --- */

    for (int i = 0; i < tape_stats[0]; ++i) {
      gradient[i] = 0.0;
    }

    for(int i = 0; i < tape_stats[1]; ++i) {
      adjoint[i] = seedVector[i];
    }

    /* --- Reverse interpretation of the computational graph --- */

    fos_reverse(1,tape_stats[1], tape_stats[0],adjoint, gradient);

    adjointVector          = gradient;
    adjointVector_Position = 0;

    inputVariables.clear();
    seedVector.clear();

    delete [] adjoint;
  }
}

#endif
