/*!
 * \file SU2_SOL.cpp
 * \brief ___________________________.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.2
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
	CSolution **solution = NULL;
	CConfig **config = NULL;
	
	/*--- Definition of the containers per zones ---*/
	solution = new CSolution*[MAX_ZONES];
	config = new CConfig*[MAX_ZONES];
	geometry = new CGeometry *[MAX_ZONES];
	
	/*--- Only one zone is allowed ---*/
	nZone = 1;
	
	for (iZone = 0; iZone < nZone; iZone++) {
		
		/*--- Definition of the configuration class per zones ---*/
		if (argc == 2) config[iZone] = new CConfig(argv[1], SU2_SOL, iZone, nZone, VERB_HIGH);
		else { strcpy (file_name, "default.cfg"); config[iZone] = new CConfig(file_name, SU2_CFD, iZone, nZone, VERB_HIGH); }
		
#ifndef NO_MPI
		/*--- Change the name of the input-output files for a parallel computation ---*/
		config[iZone]->SetFileNameDomain(rank+1);
#endif
				
		/*--- Definition of the geometry class and open the mesh file ---*/
		geometry[iZone] = new CPhysicalGeometry(config[iZone], config[iZone]->GetMesh_FileName(), config[iZone]->GetMesh_FileFormat(), iZone+1, nZone);
    }
    
#ifndef NO_MPI
    /*--- Synchronization point after the geometrical definition subroutine ---*/
    MPI::COMM_WORLD.Barrier();
#endif
    
    if (rank == MASTER_NODE)
        cout << endl <<"------------------------- Solution Postprocessing -----------------------" << endl;

  	for (iZone = 0; iZone < nZone; iZone++) {
        
		/*--- Definition of the solution class ---*/
      
    if (rank == MASTER_NODE) cout << "Read and store the solution." << endl;

		solution[iZone] = new CBaselineSolution(geometry[iZone], config[iZone], 0);
        
  	}
    
#ifndef NO_MPI
    /*--- Synchronization point after the solution subroutine ---*/
    MPI::COMM_WORLD.Barrier();
#endif
	
  if (rank == MASTER_NODE) cout << "Write the output file." << endl;

	/*--- Definition of the output class (one for all the zones) ---*/
	output = new COutput();
    output->SetBaselineResult_Files(solution, geometry, config, 0, nZone);
    

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
