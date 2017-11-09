/*!
 * \file load_hybrid_constants_test.cpp
 * \brief 
 * \author C. Pederson
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../include/hybrid_RANS_LES_model.hpp"

const std::string filename = "test_constants";
const std::string extension = ".dat";



int main() {

  /**-------------------------------------------------------------------------
   * SETUP
   *
   *--------------------------------------------------------------------------
   */
#ifdef HAVE_MPI
  MPI_Init(NULL,NULL);
#endif

  int return_flag=0;

  CConfig* test_config = new CConfig();

  const int nConstants = 15;

  /*--- Create constants files ---*/
  std::ofstream testfile;
  for (int i=0; i<3; i++) {
    std::stringstream temp_stream;
    temp_stream << filename << i << extension;
    testfile.open(temp_stream.str().c_str());
    for (int j=0; j<nConstants; j++) {
      testfile << j << std::endl;
    }
    testfile.close();
  }

  /**-------------------------------------------------------------------------
   * TEST
   *
   *--------------------------------------------------------------------------
   */

  CHybrid_Mediator mediator(3, test_config, filename);
  for (int i=0; i<3; i++) {
    for (int j=0; j<nConstants; j++) {
      if (mediator.GetConstants().at(i).at(j) != j) {
        cout << "ERROR: Constant at [" << i << "][" << j << "]";
        cout << "was wrong!" << endl;
        cout << "    Expected: " << j << endl;
        cout << "    Found;    " << mediator.GetConstants().at(i).at(j) << endl;
        return_flag = 1;
        break;
      }
    }
  }

  /**-------------------------------------------------------------------------
   * TEARDOWN
   *
   *--------------------------------------------------------------------------
   */
   delete test_config;
    


#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return return_flag;
}
