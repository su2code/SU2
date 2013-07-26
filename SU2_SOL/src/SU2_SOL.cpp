/*!
 * \file SU2_SOL.cpp
 * \brief Main file for the solution export/conversion code (SU2_SOL).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/SU2_SOL.hpp"

using namespace std;

int main(int argc, char *argv[]) {
	/*--- Variable definitions ---*/
	unsigned short iZone, nZone;
  unsigned long iExtIter, nExtIter;
	ofstream ConvHist_file;
	char file_name[200];
	int rank = MASTER_NODE;
	
#ifndef NO_MPI
	/*--- MPI initialization, and buffer setting ---*/
	void *buffer, *old_buffer;
	int size, bufsize;
	bufsize = MAX_MPI_BUFFER;
	buffer = new char[bufsize];
	MPI::Init(argc, argv);
	MPI::Attach_buffer(buffer, bufsize);
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();
#endif
  
	/*--- Pointer to different structures that will be used throughout the entire code ---*/
	COutput *output = NULL;
	CGeometry **geometry = NULL;
	CSolver **solver = NULL;
	CConfig **config = NULL;
	
	/*--- Definition of the containers per zones ---*/
	solver = new CSolver*[MAX_ZONES];
	config = new CConfig*[MAX_ZONES];
	geometry = new CGeometry *[MAX_ZONES];
	
	/*--- Only one zone is allowed ---*/
	nZone = 1;
	
	for (iZone = 0; iZone < nZone; iZone++) {
		
		/*--- Definition of the configuration class per zones ---*/
		if (argc == 2) config[iZone] = new CConfig(argv[1], SU2_SOL, iZone, nZone, VERB_HIGH);
		else { strcpy (file_name, "default.cfg"); config[iZone] = new CConfig(file_name, SU2_SOL,
                                                                          iZone, nZone, VERB_HIGH); }
		
#ifndef NO_MPI
		/*--- Change the name of the input-output files for a parallel computation ---*/
		config[iZone]->SetFileNameDomain(rank+1);
#endif
    
		/*--- Definition of the geometry class and open the mesh file ---*/
		geometry[iZone] = new CPhysicalGeometry(config[iZone], config[iZone]->GetMesh_FileName(),
                                            config[iZone]->GetMesh_FileFormat(), iZone+1, nZone);
    
    /*--- Create the vertex structure (required for MPI) ---*/
    if (rank == MASTER_NODE) cout << "Identify vertices." <<endl;
    geometry[iZone]->SetVertex(config[iZone]);
    
  }
  
#ifndef NO_MPI
  /*--- Synchronization point after the geometrical definition subroutine ---*/
  MPI::COMM_WORLD.Barrier();
#endif
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solution Postprocessing -----------------------" << endl;
  
#ifndef NO_MPI
  /*--- Synchronization point after the solution subroutine ---*/
  MPI::COMM_WORLD.Barrier();
#endif
  
	/*--- Definition of the output class (one for all the zones) ---*/
	output = new COutput();
  
  /*---  Check whether this is an unsteady simulation, and call the
   solution merging routines accordingly.---*/
  
  if (config[ZONE_0]->GetUnsteady_Simulation() && config[ZONE_0]->GetWrt_Unsteady()) {
    
    /*--- Unsteady simulation: merge all unsteady time steps. First,
     find the frequency and total number of files to write. ---*/
    
    double Physical_dt, Physical_t;
    unsigned long iExtIter = 0;
    bool StopCalc = false;
    
    while (iExtIter < config[ZONE_0]->GetnExtIter()) {
      
      /*--- Check several conditions in order to merge the correct time step files. ---*/
      Physical_dt = config[ZONE_0]->GetDelta_UnstTime();
      Physical_t  = (iExtIter+1)*Physical_dt;
      if (Physical_t >=  config[ZONE_0]->GetTotal_UnstTime())
        StopCalc = true;
        
      if ((iExtIter+1 == config[ZONE_0]->GetnExtIter()) ||
          ((iExtIter % config[ZONE_0]->GetWrt_Sol_Freq() == 0) && (iExtIter != 0) &&
           !((config[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
             (config[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND))) ||
          (StopCalc) ||
          (((config[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
            (config[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) &&
           ((iExtIter == 0) || (iExtIter % config[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) {
          
        /*--- Set the current iteration number in the config class. ---*/
        config[ZONE_0]->SetExtIter(iExtIter);
        
        /*--- Read in the restart file for this time step ---*/
        for (iZone = 0; iZone < nZone; iZone++) {
          
          /*--- Either instantiate the solution class or load a restart file. ---*/
          if (iExtIter == 0)
            solver[iZone] = new CBaselineSolver(geometry[iZone], config[iZone], MESH_0);
          else
            solver[iZone]->GetRestart(geometry[iZone], config[iZone], MESH_0);
        }
        
        if (rank == MASTER_NODE)
          cout << "Writing the volume solution for time step " << iExtIter << "." << endl;
        output->SetBaselineResult_Files(solver, geometry, config, iExtIter, nZone);
      }
      
      iExtIter++;
      if (StopCalc) break;
    }
    
  } else {
    
    /*--- Steady simulation: merge the single solution file. ---*/
    
    for (iZone = 0; iZone < nZone; iZone++) {
      /*--- Definition of the solution class ---*/
      solver[iZone] = new CBaselineSolver(geometry[iZone], config[iZone], MESH_0);      
    }
    
    output->SetBaselineResult_Files(solver, geometry, config, 0, nZone);
    
  }
  
  /*--- Deallocate the solution and output objects. ---*/
  //  for (iZone = 0; iZone < nZone; iZone++) {
  //    delete [] solver[iZone];
  //  }
  //  delete solver;
  //  delete output;
  
  
#ifndef NO_MPI
  /*--- Finalize MPI parallelization ---*/
  old_buffer = buffer;
  MPI::Detach_buffer(old_buffer);
  //	delete [] buffer;
  MPI::Finalize();
#endif
  
  /*--- End solver ---*/
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_SOL) ------------------------" << endl << endl;
  
  return EXIT_SUCCESS;
}
