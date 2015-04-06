/*!
 * \file SU2_FSI.hpp
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

#include "../include/SU2_FSI.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    
    bool StopCalc = false;
    double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
    unsigned long ExtIter = 0;
    unsigned short iMesh, iZone, iSol, nZone, nDim;
    ofstream ConvHist_file;
    int rank = MASTER_NODE;
    int size = SINGLE_NODE;
    
    
    MPI_Comm commun = MPI_COMM_WORLD;
    MPI_Group orig_group, new_group;
    
#ifdef HAVE_MPI
    /*--- MPI initialization, and buffer setting ---*/
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    
    //cout<<"\n Rank : "<<rank<<endl;
    /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
     file is specified, default.cfg is used) ---*/
    
    char config_file_name[200];
    if (argc == 2){ strcpy(config_file_name,argv[1]); }
    else{ strcpy(config_file_name, "default.cfg"); }
    
    
    /*--- Read the name and format of the input mesh file ---*/
    
    SU2_interface *SU2 = NULL;
    SU2 = new SU2_interface(config_file_name);
    
    SU2->Run_Steady_Cfd();
    
    //SU2->Run_Steady_Iteration();
    
    //SU2->Deform_Mesh();
    
    SU2->Run_Steady_Cfd();
    
    
    
    
   
}

