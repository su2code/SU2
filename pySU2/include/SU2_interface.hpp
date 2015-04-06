/*!
 * \file SU2_CFD.hpp
 * \brief Headers of the main subroutines of the code SU2_CFD.
 *        The subroutines and functions are in the <i>SU2_CFD.cpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.0 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

#ifdef HAVE_MPI
  #include "mpi.h"
#endif
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "../../SU2_CFD/include/solver_structure.hpp"
#include "../../SU2_CFD/include/integration_structure.hpp"
#include "../../SU2_CFD/include/output_structure.hpp"
#include "../../SU2_CFD/include/numerics_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../SU2_CFD/include/definition_structure.hpp"
#include "../../SU2_CFD/include/iteration_structure.hpp"



using namespace std;

/*!
 * \class SU2_interface
 * \brief Main class for defining the PDE solution, it requires
 * a child class for each particular solver (Euler, Navier-Stokes, etc.)
 * \author F. Palacios.
 * \version 3.2.0 "eagle"
 */



class SU2_interface {

    unsigned long ExtIter,ExtIter_f;
    unsigned short iMesh, iZone, iSol, nZone, nDim;
    ofstream ConvHist_file;
    bool StopCalc;
    double StartTime, StopTime, UsedTime;
    int rank,size;
    char runtime_file_name[MAX_STRING_SIZE];
    
    //MPI_Comm SU2_comm_local;
    
    public:
    
    //data members
    COutput *output;
    CIntegration ***integration_container;
    CGeometry ***geometry_container;
    CSolver ****solver_container;
    CNumerics *****numerics_container;
    CConfig **config_container;
    CSurfaceMovement **surface_movement;
    CVolumetricMovement **grid_movement;
    CFreeFormDefBox*** FFDBox;
    
    CGeometry *geometry_aux;
    
    
    //char config_file_name[200];
    
    CConfig *config;
    
    
    // Functions
    
    
    //SU2_interface(char case_filename[200],MPI_Comm SU2_comm_interface);
    SU2_interface(char case_filename[200]);
    
    void Run_Steady_Cfd(void);
    
    void Write_Output(void);
    
    void Run_Steady_Iteration(unsigned long no_of_iterations);
    
    void Deform_Mesh(void);
    
    void Write_Surface_Mesh(unsigned long iterations);
    
    void Write_Final_Output(void);
    
    //void Partition_Mesh(void);
    
};
