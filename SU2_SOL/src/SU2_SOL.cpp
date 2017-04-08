/*!
 * \file SU2_SOL.cpp
 * \brief Main file for the solution export/conversion code (SU2_SOL).
 * \author F. Palacios, T. Economon
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
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/SU2_SOL.hpp"
#include "../../SU2_CFD/include/postprocessing_structure.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	unsigned short iZone, nZone = SINGLE_ZONE;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
	ofstream ConvHist_file;
	char config_file_name[MAX_STRING_SIZE];
	int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  int nProcessor = 0;
  su2double Objective_Function;

//  ofstream CFD_pressure_file ;
//  CFD_pressure_file.open("p_CFD.dat");

  /*--- MPI initialization ---*/

#ifdef HAVE_MPI
  SU2_MPI::Init(&argc,&argv);
  SU2_Comm MPICommunicator(MPI_COMM_WORLD);
  MPI_Comm_rank(MPICommunicator,&rank);
  MPI_Comm_size(MPICommunicator,&size);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#else
  SU2_Comm MPICommunicator(0);
#endif

	/*--- Pointer to different structures that will be used throughout the entire code ---*/

	COutput *output                = NULL;
	CGeometry **geometry_container = NULL;
	CSolver **solver_container     = NULL;
	CConfig **config_container     = NULL;
////////	FWHSolver **FWH_container      = NULL;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/

  if (argc == 2 || argc == 3) { strcpy(config_file_name,argv[1]); }
  else { strcpy(config_file_name, "default.cfg"); }

  CConfig *config = NULL;
  config = new CConfig(config_file_name, SU2_SOL);

  nZone = CConfig::GetnZone(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);

	/*--- Definition of the containers per zones ---*/

	solver_container = new CSolver*[nZone];
	config_container = new CConfig*[nZone];
	geometry_container = new CGeometry*[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone]       = NULL;
    config_container[iZone]       = NULL;
    geometry_container[iZone]     = NULL;
  }

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/

    config_container[iZone] = new CConfig(config_file_name, SU2_SOL, iZone, nZone, 0, VERB_HIGH);
//    config_container[iZone]->SetMPICommunicator(MPICommunicator);

    /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/

    CGeometry *geometry_aux = NULL;

    /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

    geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);

    /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

    geometry_aux->SetColorGrid_Parallel(config_container[iZone]);

    /*--- Allocate the memory of the current domain, and
     divide the grid between the nodes ---*/

    geometry_container[iZone] = new CPhysicalGeometry(geometry_aux, config_container[iZone]);

    /*--- Deallocate the memory of geometry_aux ---*/

    delete geometry_aux;

    /*--- Add the Send/Receive boundaries ---*/

    geometry_container[iZone]->SetSendReceive(config_container[iZone]);

    /*--- Add the Send/Receive boundaries ---*/

    geometry_container[iZone]->SetBoundaries(config_container[iZone]);

    /*--- Create the vertex structure (required for MPI) ---*/

    if (rank == MASTER_NODE) cout << "Identify vertices." <<endl;
    geometry_container[iZone]->SetVertex(config_container[iZone]);

    /*--- Store the global to local mapping after preprocessing. ---*/

    if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
    geometry_container[iZone]->SetGlobal_to_Local_Point();

  }

  /*--- Determine whether the simulation is a FSI simulation ---*/

  bool fsi = config_container[ZONE_0]->GetFSI_Simulation();

  /*--- Set up a timer for performance benchmarking (preprocessing time is included) ---*/

#ifdef HAVE_MPI
  StartTime = MPI_Wtime();
#else
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#endif




////////  if (rank == MASTER_NODE)
////////          cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;

/*--- Compute elements surrounding points, points surrounding points ---*/

////////if (rank == MASTER_NODE) cout << "Setting local point connectivity." <<endl;
////////geometry_container[ZONE_0]->SetPoint_Connectivity();

////////if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
////////geometry_container[ZONE_0]->SetRCM_Ordering(config_container[ZONE_0]);

/*--- recompute elements surrounding points, points surrounding points ---*/

////////if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
////////geometry_container[ZONE_0]->SetPoint_Connectivity();

/*--- Compute elements surrounding elements ---*/

////////if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
////////geometry_container[ZONE_0]->SetElement_Connectivity();

/*--- Check the orientation before computing geometrical quantities ---*/

////////if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation of the interior elements." <<endl;
////////geometry_container[ZONE_0]->SetBoundVolume();
////////geometry_container[ZONE_0]->Check_IntElem_Orientation(config_container[ZONE_0]);
////////geometry_container[ZONE_0]->Check_BoundElem_Orientation(config_container[ZONE_0]);
/*--- Create the edge structure ---*/

////////if (rank == MASTER_NODE) cout << "Identify edges and vertices." <<endl;
////////geometry_container[ZONE_0]->SetEdges(); geometry_container[ZONE_0]->SetVertex(config_container[ZONE_0]);

/*--- Compute center of gravity ---*/

////////if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
////////geometry_container[ZONE_0]->SetCoord_CG();

/*--- Create the dual control volume structures ---*/

////////if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
////////geometry_container[ZONE_0]->SetBoundControlVolume(config_container[ZONE_0], ALLOCATE);

////////  geometry_container[ZONE_0 ]->MatchNearField(config_container[ZONE_0 ]);

////////      FWH_container = new FWHSolver* [nZone];
////////for (iZone = 0; iZone < nZone; iZone++) {

////////   FWH_container[iZone]  = new FWHSolver(config_container[iZone],geometry_container[iZone]);

////////}





  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solution Postprocessing -----------------------" << endl;

	/*--- Definition of the output class (one for all the zones) ---*/
	output = new COutput();

  /*---  Check whether this is an FSI, fluid unsteady, harmonic balance or structural dynamic simulation and call the
   solution merging routines accordingly.---*/

	if (fsi){

		if (nZone < 2){
			cout << "For multizone computations, please add the number of zones as a second argument for SU2_SOL. " << endl;
			exit(EXIT_FAILURE);
		}

		su2double Physical_dt, Physical_t;
		unsigned long iExtIter = 0, iExtIterFlow = 0, iExtIterFEM = 0;
		bool StopCalc = false;
		bool SolutionInstantiatedFlow = false, SolutionInstantiatedFEM = false;

		/*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
		if (config_container[ZONE_0]->GetRestart()){
			iExtIterFlow = config_container[ZONE_0]->GetUnst_RestartIter();
			iExtIterFEM = config_container[ZONE_1]->GetDyn_RestartIter();
			if (iExtIterFlow != iExtIterFEM) {
				if (rank == MASTER_NODE) cout << "The restart iteration is different from Flow to Structure!!" << endl;
				exit(EXIT_FAILURE);
			}
			else {
				iExtIter = iExtIterFlow;
			}
		}


		while (iExtIter < config_container[ZONE_0]->GetnExtIter()) {

		  /*--- Check several conditions in order to merge the correct time step files. ---*/

		  Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
		  Physical_t  = (iExtIter+1)*Physical_dt;
		  if (Physical_t >=  config_container[ZONE_0]->GetTotal_UnstTime())
			  StopCalc = true;

		  if (
			  ((iExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
			  ((iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0) && (iExtIter != 0) &&
			   !((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				 (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND))) ||
			  (StopCalc) ||
			  (((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) &&
			   ((iExtIter == 0) || (iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0))))

			 &&

			 ((iExtIter+1 == config_container[ZONE_1]->GetnExtIter()) ||
			  (StopCalc) ||
			  ((config_container[ZONE_1]->GetDynamic_Analysis() == DYNAMIC) &&
			  ((iExtIter == 0) || (iExtIter % config_container[ZONE_1]->GetWrt_Sol_Freq_DualTime() == 0))))

			 ){

				/*--- Set the current iteration number in the config class. ---*/
				config_container[ZONE_0]->SetExtIter(iExtIter);
				config_container[ZONE_1]->SetExtIter(iExtIter);

				/*--- Read in the restart file for this time step ---*/

				/*--- For the fluid zone (ZONE_0) ---*/
				/*--- Either instantiate the solution class or load a restart file. ---*/
        if (SolutionInstantiatedFlow == false &&
            (iExtIter == 0 || ((config_container[ZONE_0]->GetRestart() && (SU2_TYPE::Int(iExtIter) == config_container[ZONE_0]->GetUnst_RestartIter())) ||
														iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0 ||
														iExtIter+1 == config_container[ZONE_0]->GetnExtIter()))) {
          solver_container[ZONE_0] = new CBaselineSolver(geometry_container[ZONE_0], config_container[ZONE_0]);
					SolutionInstantiatedFlow = true;
				}
					solver_container[ZONE_0]->LoadRestart_FSI(geometry_container[ZONE_0], &solver_container, config_container[ZONE_0], SU2_TYPE::Int(MESH_0));


				/*--- For the structural zone (ZONE_1) ---*/
				/*--- Either instantiate the solution class or load a restart file. ---*/
				/*--- Either instantiate the solution class or load a restart file. ---*/
        if (SolutionInstantiatedFEM == false &&
            (iExtIter == 0 || ((config_container[ZONE_1]->GetRestart() && (SU2_TYPE::Int(iExtIter) == config_container[ZONE_1]->GetDyn_RestartIter())) ||
																iExtIter % config_container[ZONE_1]->GetWrt_Sol_Freq_DualTime() == 0 ||
																iExtIter+1 == config_container[ZONE_1]->GetnExtIter()))) {
          solver_container[ZONE_1] = new CBaselineSolver(geometry_container[ZONE_1], config_container[ZONE_1]);
					SolutionInstantiatedFEM = true;
				}
					solver_container[ZONE_1]->LoadRestart_FSI(geometry_container[ZONE_1], &solver_container, config_container[ZONE_1], SU2_TYPE::Int(MESH_0));

				if (rank == MASTER_NODE) cout << "Writing the volume solution for time step " << iExtIter << "." << endl;
				output->SetBaselineResult_Files(solver_container, geometry_container, config_container, iExtIter, nZone);
		  }

		  iExtIter++;
		  if (StopCalc) break;
		}

	}
	else {

		if (config_container[ZONE_0]->GetWrt_Unsteady()) {

			/*--- Unsteady simulation: merge all unsteady time steps. First,
     	 	 	  find the frequency and total number of files to write. ---*/

			su2double Physical_dt, Physical_t;
			unsigned long iExtIter = 0;
			bool StopCalc = false;
			bool *SolutionInstantiated = new bool[nZone];

			for (iZone = 0; iZone < nZone; iZone++)
				SolutionInstantiated[iZone] = false;


			if (config_container[ZONE_0]->GetAD_Mode()) AD::StartRecording();


			/*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
			if (config_container[ZONE_0]->GetWrt_Unsteady() && config_container[ZONE_0]->GetRestart())
				iExtIter = config_container[ZONE_0]->GetUnst_RestartIter();

			while (iExtIter < config_container[ZONE_0]->GetnExtIter()) {

				/*--- Check several conditions in order to merge the correct time step files. ---*/
				Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
				Physical_t  = (iExtIter+1)*Physical_dt;
				if (Physical_t >=  config_container[ZONE_0]->GetTotal_UnstTime())
					StopCalc = true;

				if ((iExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
						((iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0) && (iExtIter != 0) &&
								!((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
										(config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND))) ||
										(StopCalc) ||
										(((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
												(config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) &&
												((iExtIter == 0) || (iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) {



					/*--- Read in the restart file for this time step ---*/
					for (iZone = 0; iZone < nZone; iZone++) {

						/*--- Set the current iteration number in the config class. ---*/
						config_container[iZone]->SetExtIter(iExtIter);

						/*--- Either instantiate the solution class or load a restart file. ---*/
                if (SolutionInstantiated[iZone] == false &&
                    (iExtIter == 0 || (config_container[ZONE_0]->GetRestart() && ((long)iExtIter == config_container[ZONE_0]->GetUnst_RestartIter() ||
										iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0 ||
										iExtIter+1 == config_container[ZONE_0]->GetnExtIter())))) {
                  solver_container[iZone] = new CBaselineSolver(geometry_container[iZone], config_container[iZone]);
							SolutionInstantiated[iZone] = true;
						}
                  solver_container[iZone]->LoadRestart(geometry_container, &solver_container, config_container[iZone], SU2_TYPE::Int(MESH_0), true);
					}
					//Extracting flow data on the FWH surface as well as static pressure at the observer locations (for validation only)
////////					 FWH_container[ZONE_0]->SetAeroacoustic_Analysis(solver_container[ZONE_0],config_container[ZONE_0],geometry_container[ZONE_0],CFD_pressure_file);


					if (rank == MASTER_NODE  && config_container[ZONE_0]->GetKind_ObjFunc() != NOISE  ){  // && config_container[ZONE_0]->GetKind_ObjFunc() == NOISE
						cout << "Writing the volume solution for time step " << iExtIter << "." << endl;
					output->SetBaselineResult_Files(solver_container, geometry_container, config_container, iExtIter, nZone);
				}
				}

				iExtIter++;
				if (StopCalc) break;
			}

			if (config_container[ZONE_0]->GetAD_Mode()) AD::StopRecording();

		}

		else if (config_container[ZONE_0]->GetUnsteady_Simulation() == HARMONIC_BALANCE) {

			/*--- Read in the restart file for this time step ---*/
			for (iZone = 0; iZone < nZone; iZone++) {

				/*--- Either instantiate the solution class or load a restart file. ---*/
        solver_container[iZone] = new CBaselineSolver(geometry_container[iZone], config_container[iZone]);
        solver_container[iZone]->LoadRestart(geometry_container, &solver_container, config_container[iZone], SU2_TYPE::Int(MESH_0), true);

				/*--- Print progress in solution writing to the screen. ---*/
				if (rank == MASTER_NODE) {
          cout << "Storing the volume solution for time instance " << iZone << "." << endl;
				}

			}

			output->SetBaselineResult_Files(solver_container, geometry_container, config_container, iZone, nZone);
		}

        else if (config_container[ZONE_0]->GetWrt_Dynamic()){

			/*--- Dynamic simulation: merge all unsteady time steps. First,
			 find the frequency and total number of files to write. ---*/

			su2double Physical_dt, Physical_t;
			unsigned long iExtIter = 0;
			bool StopCalc = false;
			bool SolutionInstantiated = false;



			/*--- Check for an dynamic restart (structural analysis). Update ExtIter if necessary. ---*/
			if (config_container[ZONE_0]->GetKind_Solver() == FEM_ELASTICITY &&
				  config_container[ZONE_0]->GetWrt_Dynamic() && config_container[ZONE_0]->GetRestart())
				  iExtIter = config_container[ZONE_0]->GetDyn_RestartIter();

			while (iExtIter < config_container[ZONE_0]->GetnExtIter()) {

			  /*--- Check several conditions in order to merge the correct time step files. ---*/
			  /*--- If the solver is structural, the total and delta_t are obtained from different functions. ---*/

				Physical_dt = config_container[ZONE_0]->GetDelta_DynTime();
				Physical_t  = (iExtIter+1)*Physical_dt;
				if (Physical_t >=  config_container[ZONE_0]->GetTotal_DynTime())
				  StopCalc = true;

			  if ((iExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
				  (StopCalc) ||
				  ((config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC) &&
				  ((iExtIter == 0) || (iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) {

				/*--- Set the current iteration number in the config class. ---*/
				config_container[ZONE_0]->SetExtIter(iExtIter);

				/*--- Read in the restart file for this time step ---*/
				for (iZone = 0; iZone < nZone; iZone++) {

				  /*--- Either instantiate the solution class or load a restart file. ---*/
                if (SolutionInstantiated == false &&
                    (iExtIter == 0 || ((config_container[ZONE_0]->GetRestart() && (SU2_TYPE::Int(iExtIter) == config_container[ZONE_0]->GetDyn_RestartIter())) ||
																  iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0 ||
																  iExtIter+1 == config_container[ZONE_0]->GetnExtIter()))) {
                  solver_container[iZone] = new CBaselineSolver(geometry_container[iZone], config_container[iZone]);
					SolutionInstantiated = true;
				  }
                  solver_container[iZone]->LoadRestart(geometry_container, &solver_container, config_container[iZone], SU2_TYPE::Int(MESH_0), true);
				}

					if (rank == MASTER_NODE)
				  cout << "Writing the volume solution for time step " << iExtIter << "." << endl;
				output->SetBaselineResult_Files(solver_container, geometry_container, config_container, iExtIter, nZone);
			  }

			  iExtIter++;
			  if (StopCalc) break;
			}

		  }

    else {

			  /*--- Steady simulation: merge the single solution file. ---*/

			  for (iZone = 0; iZone < nZone; iZone++) {
				  /*--- Definition of the solution class ---*/
        solver_container[iZone] = new CBaselineSolver(geometry_container[iZone], config_container[iZone]);
        solver_container[iZone]->LoadRestart(geometry_container, &solver_container, config_container[iZone], SU2_TYPE::Int(MESH_0), true);
			  }

			  output->SetBaselineResult_Files(solver_container, geometry_container, config_container, 0, nZone);

		  }

	}



        if (rank == MASTER_NODE)
          cout<<"Type= "<<   config_container[ZONE_0]->GetDiscrete_Adjoint() <<endl;  //  <--- returns 1! (cont adj)

         if (config_container[ZONE_0]->GetAD_Mode()){
/*             if (rank == MASTER_NODE)
               cout << endl <<"------------------------- Computing Far Field Noise (Primal+Adjoint) -----------------------" << endl;

             AD::StartRecording();
             FWH_container[ZONE_0]->Compute_FarfieldNoise(solver_container[ZONE_0],config_container[ZONE_0],geometry_container[ZONE_0]);
             Objective_Function = FWH_container[ZONE_0]-> SPL;

             //cout<<Objective_Function.getGradientData()<<endl;
             if (rank==MASTER_NODE){
             SU2_TYPE::SetDerivative(Objective_Function,1.0);
               }else{
                 SU2_TYPE::SetDerivative(Objective_Function,0.0);
               }
             AD::StopRecording();
             AD::ComputeAdjoint();

             cout<<"Finished Computing FWH Adjoint"<<endl;

             su2double extracted_derivative;

             for (int iSample=0; iSample<FWH_container[ZONE_0]->nSample; iSample++){
                for (int iPanel=0; iPanel<FWH_container[ZONE_0]->nPanel; iPanel++){
                  for (int i =0; i< FWH_container[ZONE_0]->nDim+3; i++){
                         FWH_container[ZONE_0]-> dJdU[i][iPanel][iSample]=SU2_TYPE::GetDerivative(extracted_derivative);
//                         cout<<"dJdU= "<< FWH_container[ZONE_0]-> dJdU[i][iPanel][iSample]<<endl;
                    }
               }
             }


             cout<<"Finished Extracting"<<endl;

            FWH_container[ZONE_0]->Write_Sensitivities(solver_container[ZONE_0],config_container[ZONE_0],geometry_container[ZONE_0]);

*/

             if (rank == MASTER_NODE)
               cout << endl <<"------------------------- Computing Far Field Noise (Primal+Adjoint) -----------------------" << endl;

             AD::StartRecording();
			 SUBoom boom(solver_container[ZONE_0], config_container[ZONE_0], geometry_container[ZONE_0]);
             if (rank == MASTER_NODE){
               cout << "SUBoom initialized." << endl;
             boom.ConditionAtmosphericData();
               cout << "Condition atmospheric data complete." << endl;
             boom.ScaleFactors();
               cout << "Scale factors complete." << endl;
             boom.InitialWaveNormals();
               cout << "Initial wave normals complete." << endl;
             boom.RayTracer();
               cout << "Ray tracer complete." << endl;
             boom.RayTubeArea();
               cout << "Ray tube area complete." << endl;
             boom.FindInitialRayTime();
               cout << "Find initial ray time complete." << endl;
             boom.ODETerms();
               cout << "ODE terms complete." << endl;
             boom.DistanceToTime();
               cout << "Distance to time complete." << endl;
             boom.CreateSignature();
               cout << "Create signature complete." << endl;
             boom.PropagateSignal();
             }
#ifdef HAVE_MPI
  SU2_MPI::Bcast(boom.p_int2, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#endif

             Objective_Function = boom.p_int2;
             if (rank==MASTER_NODE){
             SU2_TYPE::SetDerivative(Objective_Function,1.0);
               }else{
                 SU2_TYPE::SetDerivative(Objective_Function,0.0);
               }
             AD::StopRecording();
             AD::ComputeAdjoint();

             if (rank==MASTER_NODE)
               cout<<"Finished computing boom adjoint."<<endl;

            su2double extracted_derivative;

             for (int iSig=0; iSig<boom.nSig; iSig++){
                for (int i =0; i< boom.nDim+3; i++){
                    boom.dJdU[i][iSig]=SU2_TYPE::GetDerivative(extracted_derivative);
                }
             }

             if(rank==MASTER_NODE){
               cout<<"Finished extracting."<<endl;
             boom.WriteSensitivities(solver_container[ZONE_0],config_container[ZONE_0],geometry_container[ZONE_0]);
             }

           }
		   else{
             if (rank == MASTER_NODE)
               cout << endl <<"------------------------- Computing Far Field Noise (Primal Only) -----------------------" << endl;
////////             FWH_container[ZONE_0]->Compute_FarfieldNoise(solver_container[ZONE_0],config_container[ZONE_0],geometry_container[ZONE_0]);

             SUBoom boom(solver_container[ZONE_0], config_container[ZONE_0], geometry_container[ZONE_0]);
             if (rank == MASTER_NODE)
               cout << "SUBoom initialized." << endl;
             boom.ConditionAtmosphericData();
             if (rank == MASTER_NODE)
               cout << "Condition atmospheric data complete." << endl;
             boom.ScaleFactors();
             if (rank == MASTER_NODE)
               cout << "Scale factors complete." << endl;
             boom.InitialWaveNormals();
             if (rank == MASTER_NODE)
               cout << "Initial wave normals complete." << endl;
             boom.RayTracer();
             if (rank == MASTER_NODE)
               cout << "Ray tracer complete." << endl;
             boom.RayTubeArea();
             if (rank == MASTER_NODE)
               cout << "Ray tube area complete." << endl;
             boom.FindInitialRayTime();
             if (rank == MASTER_NODE)
               cout << "Find initial ray time complete." << endl;
             boom.ODETerms();
             if (rank == MASTER_NODE)
               cout << "ODE terms complete." << endl;
             boom.DistanceToTime();
             if (rank == MASTER_NODE)
               cout << "Distance to time complete." << endl;
             boom.CreateSignature();
             if (rank == MASTER_NODE)
               cout << "Create signature complete." << endl;
             boom.PropagateSignal();
           }


  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

#ifdef HAVE_MPI
  StopTime = MPI_Wtime();
#else
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#endif

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime-StartTime;
  if (rank == MASTER_NODE) {
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on "<< size;
    if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
  }

  /*--- Exit the solver cleanly ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_SOL) ------------------------" << endl << endl;

  /*--- Finalize MPI parallelization ---*/

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
