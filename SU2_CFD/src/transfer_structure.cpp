/*!
 * \file transfer_structure.cpp
 * \brief Main subroutines for physical transfer of information
 * \author R. Sanchez
 * \version 4.0.1 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
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

#include "../include/transfer_structure.hpp"

CTransfer::CTransfer(void) {

	Physical_Constants = NULL;
	Donor_Variable     = NULL;
	Target_Variable    = NULL;

	nVar = 0;

}

CTransfer::CTransfer(unsigned short val_nVar, unsigned short val_nConst, CConfig *config){

	unsigned short iVar;

	Physical_Constants = new su2double[val_nConst];
	Donor_Variable     = new su2double[val_nVar];
	Target_Variable    = new su2double[val_nVar];

	nVar = val_nVar;

	for (iVar = 0; iVar < nVar; iVar++){
		Donor_Variable[iVar]  = 0.0;
		Target_Variable[iVar] = 0.0;
	}

	for (iVar = 0; iVar < val_nConst; iVar++){
		Physical_Constants[iVar] = 0.0;
	}

}

CTransfer::~CTransfer(void) {

	if (Physical_Constants   != NULL) delete [] Physical_Constants;
	if (Donor_Variable       != NULL) delete [] Donor_Variable;
	if (Target_Variable      != NULL) delete [] Target_Variable;

}

void CTransfer::Scatter_InterfaceData_Matching(CSolver *donor_solution, CSolver *target_solution,
		   	   	   	   	   	   	   	   	   	   CGeometry *donor_geometry, CGeometry *target_geometry,
											   CConfig *donor_config, CConfig *target_config){

	unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;		// Number of markers on the interface, donor and target side
	unsigned short iMarkerInt, iMarkerDonor, iMarkerTarget;		// Variables for iteration over markers
	int Marker_Donor = -1, Marker_Target = -1;

	unsigned long nVertexDonor, nVertexTarget;					// Number of vertices on Donor and Target side
	unsigned long iVertex, iPoint;								// Variables for iteration over vertices and nodes
	unsigned long jVertex, jPoint;								// Variables for iteration over vertices and nodes

	unsigned short iVar, jDim;



	GetPhysical_Constants(donor_solution, target_solution, donor_geometry, target_geometry,
						  donor_config, target_config);

	unsigned long Check_Point_Global;
	unsigned long Point_Donor, Point_Target;
	su2double *Normal_Donor, *Normal_Target;

    int rank = MASTER_NODE;
    int size = SINGLE_NODE;

#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    unsigned long nLocalVertexDonor = 0, nLocalVertexTarget= 0;
    unsigned long iVertexDonor = 0, iVertexTarget;
    unsigned long nPoint_Total = 0;

	unsigned long MaxLocalVertexDonor = 0, MaxLocalVertexTarget = 0;

	unsigned long nBuffer_DonorVariables = 0, nBuffer_TargetVariables = 0;
	unsigned long nBuffer_DonorIndices = 0, nBuffer_TargetIndices = 0;

	unsigned long Processor_Donor, Processor_Target;

	int iProcessor, nProcessor = 0;
	unsigned long iVariable;

	/*--- Number of markers on the FSI interface ---*/

	nMarkerInt     = (donor_config->GetMarker_n_FSIinterface())/2;
	nMarkerTarget  = target_geometry->GetnMarker();
	nMarkerDonor   = donor_geometry->GetnMarker();

	nProcessor = size;

	/*--- Outer loop over the markers on the FSI interface: compute one by one ---*/
	/*--- The tags are always an integer greater than 1: loop from 1 to nMarkerFSI ---*/

	for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++){

		Marker_Donor = -1;
		Marker_Target = -1;

		/*--- Initialize pointer buffers inside the loop, so we can delete for each marker. ---*/
		unsigned long Buffer_Send_nVertexDonor[1], *Buffer_Recv_nVertexDonor = NULL;
		unsigned long Buffer_Send_nVertexTarget[1], *Buffer_Recv_nVertexTarget = NULL;

		/*--- The donor and target markers are tagged with the same index.
		 *--- This is independent of the MPI domain decomposition.
		 *--- We need to loop over all markers on both sides and get the number of nodes
		 *--- that belong to each FSI marker for each processor ---*/

		/*--- On the donor side ---*/

		for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++){
			/*--- If the tag GetMarker_All_FSIinterface(iMarkerDonor) equals the index we are looping at ---*/
			if ( donor_config->GetMarker_All_FSIinterface(iMarkerDonor) == iMarkerInt ){
				/*--- We have identified the local index of the Donor marker ---*/
				/*--- Store the number of local points that belong to Marker_Donor on each processor ---*/
				/*--- This are the number of points that will be sent from this particular processor ---*/
				/*--- This includes the halo nodes ---*/
				nLocalVertexDonor = donor_geometry->GetnVertex(iMarkerDonor);
				/*--- Store the identifier for the structural marker ---*/
				Marker_Donor = iMarkerDonor;
				/*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
				break;
			}
			else {
				/*--- If the tag hasn't matched any tag within the donor markers ---*/
				nLocalVertexDonor = 0;
				Marker_Donor = -1;
			}
		}

		/*--- On the target side ---*/

		for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++){
			/*--- If the tag GetMarker_All_FSIinterface(iMarkerFlow) equals the index we are looping at ---*/
			if ( target_config->GetMarker_All_FSIinterface(iMarkerTarget) == iMarkerInt ){
				/*--- We have identified the local index of the Flow marker ---*/
				/*--- Store the number of local points that belong to Marker_Flow on each processor ---*/
				/*--- This includes the halo nodes ---*/
				nLocalVertexTarget = target_geometry->GetnVertex(iMarkerTarget);
				/*--- Store the identifier for the fluid marker ---*/
				Marker_Target = iMarkerTarget;
				/*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
				break;
			}
			else {
				/*--- If the tag hasn't matched any tag within the Flow markers ---*/
				nLocalVertexTarget = 0;
				Marker_Target = -1;
			}
		}

		Buffer_Send_nVertexDonor[0] = nLocalVertexDonor;							   // Retrieve total number of vertices on Donor marker
		Buffer_Send_nVertexTarget[0] = nLocalVertexTarget;							   // Retrieve total number of vertices on Target marker
		if (rank == MASTER_NODE) Buffer_Recv_nVertexDonor = new unsigned long[size];   // Allocate memory to receive how many vertices are on each rank on the structural side
		if (rank == MASTER_NODE) Buffer_Recv_nVertexTarget = new unsigned long[size];  // Allocate memory to receive how many vertices are on each rank on the fluid side
#ifdef HAVE_MPI
		/*--- We receive MaxLocalVertexFEA as the maximum number of vertices in one single processor on the structural side---*/
		SU2_MPI::Allreduce(&nLocalVertexDonor, &MaxLocalVertexDonor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
		/*--- We receive MaxLocalVertexFlow as the maximum number of vertices in one single processor on the fluid side ---*/
		SU2_MPI::Allreduce(&nLocalVertexTarget, &MaxLocalVertexTarget, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

		/*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the structural side ---*/
		SU2_MPI::Gather(&Buffer_Send_nVertexDonor, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexDonor, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
		/*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the fluid side ---*/
		SU2_MPI::Gather(&Buffer_Send_nVertexTarget, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexTarget, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
		MaxLocalVertexDonor  = nLocalVertexDonor;
		MaxLocalVertexTarget = nLocalVertexTarget;

		Buffer_Recv_nVertexDonor[0] = nLocalVertexDonor;
		Buffer_Recv_nVertexTarget[0] = nLocalVertexTarget;

#endif

		/*--- We will be gathering the structural coordinates into the master node ---*/
		/*--- Then we will distribute them using a scatter operation into the appropriate fluid processor ---*/
		nBuffer_DonorVariables = MaxLocalVertexDonor * nVar;
		nBuffer_TargetVariables = MaxLocalVertexTarget * nVar;

		/*--- We will be gathering donor index and donor processor (for flow -> donor = structure) ---*/
		/*--- Then we will pass on to the structural side the index (fea point) to the appropriate processor ---*/
		nBuffer_DonorIndices = 2 * MaxLocalVertexDonor;
		nBuffer_TargetIndices = MaxLocalVertexTarget;

		/*--- Send and Recv buffers ---*/

		/*--- Buffers to send and receive the variables in the donor mesh ---*/
		su2double *Buffer_Send_DonorVariables = new su2double[nBuffer_DonorVariables];
		su2double *Buffer_Recv_DonorVariables = NULL;

		/*--- Buffers to send and receive the indices in the donor mesh ---*/
		long *Buffer_Send_DonorIndices = new long[nBuffer_DonorIndices];
		long *Buffer_Recv_DonorIndices = NULL;

		/*--- Buffers to send and receive the variables in the target mesh---*/
		su2double *Buffer_Send_TargetVariables = NULL;
		su2double *Buffer_Recv_TargetVariables = new su2double[nBuffer_TargetVariables];

		/*--- Buffers to send and receive the target indices ---*/
		long *Buffer_Send_TargetIndices = NULL;
		long *Buffer_Recv_TargetIndices = new long[nBuffer_TargetIndices];

		/*--- Prepare the receive buffers (1st step) and send buffers (2nd step) on the master node only. ---*/

		if (rank == MASTER_NODE) {
			Buffer_Recv_DonorVariables  = new su2double[size*nBuffer_DonorVariables];
			Buffer_Recv_DonorIndices    = new long[size*nBuffer_DonorIndices];
			Buffer_Send_TargetVariables = new su2double[size*nBuffer_TargetVariables];
			Buffer_Send_TargetIndices   = new long[size*nBuffer_TargetIndices];
		}

		/*--- On the fluid side ---*/

		/*--- If this processor owns the marker we are looping at on the structural side ---*/

		/*--- First we initialize all of the indices and processors to -1 ---*/
		/*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
		for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++)
			Buffer_Send_DonorIndices[iVertex] = -1;

		if (Marker_Donor >= 0){

			/*--- We have identified the local index of the FEA marker ---*/
			/*--- We loop over all the vertices in that marker and in that particular processor ---*/

			for (iVertex = 0; iVertex < nLocalVertexDonor; iVertex++){

		        Point_Donor = donor_geometry->vertex[Marker_Donor][iVertex]->GetNode();

		        Point_Target = donor_geometry->vertex[Marker_Donor][iVertex]->GetDonorPoint();

		        Processor_Target = donor_geometry->vertex[Marker_Donor][iVertex]->GetDonorProcessor();

				GetDonor_Variable(donor_solution, donor_geometry, donor_config, Marker_Donor, iVertex, Point_Donor);

				for (iVar = 0; iVar < nVar; iVar++){
					Buffer_Send_DonorVariables[iVertex*nVar+iVar] = Donor_Variable[iVar];
				}
				/*--- If this processor owns the node ---*/
				if (donor_geometry->node[Point_Donor]->GetDomain()){
					Buffer_Send_DonorIndices[2*iVertex]     = Point_Target;
					Buffer_Send_DonorIndices[2*iVertex + 1] = Processor_Target;
				}
				else{
					/*--- We set the values to be -1 to be able to identify them later as halo nodes ---*/
					Buffer_Send_DonorIndices[2*iVertex]     = -1;
					Buffer_Send_DonorIndices[2*iVertex + 1] = -1;
				}

			}
		}

#ifdef HAVE_MPI
		/*--- Once all the messages have been sent, we gather them all into the MASTER_NODE ---*/
		SU2_MPI::Gather(Buffer_Send_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, Buffer_Recv_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
		SU2_MPI::Gather(Buffer_Send_DonorIndices, nBuffer_DonorIndices, MPI_LONG, Buffer_Recv_DonorIndices, nBuffer_DonorIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);

#else
		for (iVariable = 0; iVariable < nBuffer_DonorVariables; iVariable++)
			Buffer_Recv_DonorVariables[iVariable] = Buffer_Send_DonorVariables[iVariable];
		for (iVariable = 0; iVariable < nBuffer_DonorIndices; iVariable++)
			Buffer_Recv_DonorIndices[iVariable] = Buffer_Send_DonorIndices[iVariable];
#endif

		/*--- Counter to determine where in the array we have to set the information ---*/
		long *Counter_Processor_Target = NULL;
		long iProcessor_Donor = 0, iIndex_Donor = 0;
		long iProcessor_Target = 0, iPoint_Target = 0, iIndex_Target = 0;

		/*--- Now we pack the information to send it over to the different processors ---*/

		if (rank == MASTER_NODE){

			/*--- We set the counter to 0 ---*/
			Counter_Processor_Target = new long[nProcessor];
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
				Counter_Processor_Target[iProcessor] = 0;
			}

			/*--- First we initialize the index vector to -1 ---*/
			/*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
			for (iVertex = 0; iVertex < nProcessor*nBuffer_TargetIndices; iVertex++)
				Buffer_Send_TargetIndices[iVertex] = -2;

			/*--- As of now we do the loop over the flow points ---*/
			/*--- The number of points for flow and structure does not necessarily have to match ---*/
			/*--- In fact, it's possible that a processor asks for nStruct nodes and there are only ---*/
			/*--- nFlow < nStruct available; this is due to halo nodes ---*/

			/*--- For every processor from which we have received information ---*/
			/*--- (This is, for every processor on the structural side) ---*/
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){

				/*--- This is the initial index on the coordinates buffer for that particular processor on the structural side ---*/
				iProcessor_Donor = iProcessor*nBuffer_DonorVariables;
				/*--- This is the initial index on the donor index/processor buffer for that particular processor on the structural side ---*/
				iIndex_Donor = iProcessor*nBuffer_DonorIndices;

				/*--- For every vertex in the information retreived from iProcessor ---*/
				for (iVertex = 0; iVertex < Buffer_Recv_nVertexDonor[iProcessor]; iVertex++) {

					/*--- The processor and index for the flow are: ---*/
					Processor_Target = Buffer_Recv_DonorIndices[iIndex_Donor+iVertex*2+1];
					Point_Target     = Buffer_Recv_DonorIndices[iIndex_Donor+iVertex*2];

					/*--- Load the buffer at the appropriate position ---*/
					/*--- This is determined on the fluid side by:
					 *--- Processor_Target*nBuffer_StructTraction -> Initial position of the processor array (fluid side)
					 *--- +
					 *--- Counter_Processor_Struct*nVar -> Initial position of the nVar array for the particular point on the fluid side
					 *--- +
					 *--- iVar -> Position within the nVar array that corresponds to a point
					 *---
					 *--- While on the structural side is:
					 *--- iProcessor*nBuffer_FlowTraction -> Initial position on the processor array (structural side)
					 *--- +
					 *--- iVertex*nVar -> Initial position of the nVar array for the particular point on the structural side
					 */

					/*--- We check that we are not setting the value for a halo node ---*/
					if (Point_Target != -1){
						iProcessor_Target = Processor_Target*nBuffer_TargetVariables;
						iIndex_Target = Processor_Target*nBuffer_TargetIndices;
						iPoint_Target = Counter_Processor_Target[Processor_Target]*nVar;

						for (iVar = 0; iVar < nVar; iVar++)
							Buffer_Send_TargetVariables[iProcessor_Target + iPoint_Target + iVar] = Buffer_Recv_DonorVariables[iProcessor_Donor + iVertex*nVar + iVar];

						/*--- We set the fluid index at an appropriate position matching the coordinates ---*/
						Buffer_Send_TargetIndices[iIndex_Target + Counter_Processor_Target[Processor_Target]] = Point_Target;

						Counter_Processor_Target[Processor_Target]++;
					}

				}

			}

		}

#ifdef HAVE_MPI
		/*--- Once all the messages have been prepared, we scatter them all from the MASTER_NODE ---*/
		SU2_MPI::Scatter(Buffer_Send_TargetVariables, nBuffer_TargetVariables, MPI_DOUBLE, Buffer_Recv_TargetVariables, nBuffer_TargetVariables, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
		SU2_MPI::Scatter(Buffer_Send_TargetIndices, nBuffer_TargetIndices, MPI_LONG, Buffer_Recv_TargetIndices, nBuffer_TargetIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
		for (iVariable = 0; iVariable < nBuffer_TargetVariables; iVariable++)
			Buffer_Recv_TargetVariables[iVariable] = Buffer_Send_TargetVariables[iVariable];
		for (iVariable = 0; iVariable < nBuffer_TargetIndices; iVariable++)
			Buffer_Recv_TargetIndices[iVariable] = Buffer_Send_TargetIndices[iVariable];
#endif

		long indexPoint_iVertex, Point_Target_Check;

		/*--- For the target marker we are studying ---*/
		if (Marker_Target >= 0){

			/*--- We have identified the local index of the Structural marker ---*/
			/*--- We loop over all the vertices in that marker and in that particular processor ---*/

			for (iVertex = 0; iVertex < nLocalVertexTarget; iVertex++){

				Point_Target = target_geometry->vertex[Marker_Target][iVertex]->GetNode();

				if (target_geometry->node[Point_Target]->GetDomain()){
					/*--- Find the index of the point Point_Struct in the buffer Buffer_Recv_SetIndex ---*/
					indexPoint_iVertex = std::distance(Buffer_Recv_TargetIndices, std::find(Buffer_Recv_TargetIndices, Buffer_Recv_TargetIndices + MaxLocalVertexTarget, Point_Target));

					Point_Target_Check = Buffer_Recv_TargetIndices[indexPoint_iVertex];

					if (Point_Target_Check < 0) {
						cout << "WARNING: A nonphysical point is being considered for traction transfer." << endl;
						exit(EXIT_FAILURE);
					}

					for (iVar = 0; iVar < nVar; iVar++)
						Target_Variable[iVar] = Buffer_Recv_TargetVariables[indexPoint_iVertex*nVar+iVar];

					SetTarget_Variable(target_solution, target_geometry, target_config, Marker_Target, iVertex, Point_Target);

				}

			}

		}

		delete [] Buffer_Send_DonorVariables;
		delete [] Buffer_Send_DonorIndices;
		delete [] Buffer_Recv_TargetVariables;
		delete [] Buffer_Recv_TargetIndices;

		if (rank == MASTER_NODE) {
			delete [] Buffer_Recv_nVertexDonor;
			delete [] Buffer_Recv_nVertexTarget;
			delete [] Buffer_Recv_DonorVariables;
			delete [] Buffer_Recv_DonorIndices;
			delete [] Buffer_Send_TargetVariables;
			delete [] Buffer_Send_TargetIndices;
			delete [] Counter_Processor_Target;
		}

	}


}

void CTransfer::Broadcast_InterfaceData_Matching(CSolver *donor_solution, CSolver *target_solution,
		   	   	   	   	   	   	   	   	   	   	 CGeometry *donor_geometry, CGeometry *target_geometry,
												 CConfig *donor_config, CConfig *target_config){



	unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;		// Number of markers on the interface, donor and target side
	unsigned short iMarkerInt, iMarkerDonor, iMarkerTarget;		// Variables for iteration over markers
	int Marker_Donor = -1, Marker_Target = -1;

	unsigned long nVertexDonor, nVertexTarget;					// Number of vertices on Donor and Target side
	unsigned long iVertex, iPoint;								// Variables for iteration over vertices and nodes
	unsigned long jVertex, jPoint;								// Variables for iteration over vertices and nodes

	unsigned short iVar, jDim;

	GetPhysical_Constants(donor_solution, target_solution, donor_geometry, target_geometry,
						  donor_config, target_config);

	unsigned long Point_Donor_Global, Donor_Global_Index;
	unsigned long Point_Donor, Point_Target;
	su2double *Normal_Donor, *Normal_Target;

    int rank = MASTER_NODE;
    int size = SINGLE_NODE;

#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    unsigned long iLocalVertex = 0;
    unsigned long nLocalVertexDonor = 0, nLocalVertexDonorOwned = 0;
    unsigned long iVertexDonor = 0, iVertexTarget = 0;
    unsigned long nPoint_Total = 0;

	unsigned long MaxLocalVertexDonor = 0, MaxLocalVertexTarget = 0;
	unsigned long TotalVertexDonor = 0;

	unsigned long nBuffer_DonorVariables = 0, nBuffer_TargetVariables = 0;
	unsigned long nBuffer_DonorIndices = 0, nBuffer_TargetIndices = 0;

	unsigned long nBuffer_BcastVariables = 0, nBuffer_BcastIndices = 0;

	unsigned long Processor_Donor, Processor_Target;

	int iProcessor, nProcessor = 0;
	unsigned long iVariable;

	/*--- Number of markers on the FSI interface ---*/

	nMarkerInt     = (donor_config->GetMarker_n_FSIinterface())/2;
	nMarkerTarget  = target_geometry->GetnMarker();
	nMarkerDonor   = donor_geometry->GetnMarker();

	nProcessor = size;

	/*--- Outer loop over the markers on the FSI interface: compute one by one ---*/
	/*--- The tags are always an integer greater than 1: loop from 1 to nMarkerFSI ---*/

	for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++){

		Marker_Donor = -1;
		Marker_Target = -1;

		/*--- Initialize pointer buffers inside the loop, so we can delete for each marker. ---*/
		unsigned long Buffer_Send_nVertexDonor[1], *Buffer_Recv_nVertexDonor = NULL;
		unsigned long Buffer_Send_nVertexDonorOwned[1], *Buffer_Recv_nVertexDonorOwned = NULL;
		unsigned long Buffer_Send_nVertexTarget[1], *Buffer_Recv_nVertexTarget = NULL;

		/*--- The donor and target markers are tagged with the same index.
		 *--- This is independent of the MPI domain decomposition.
		 *--- We need to loop over all markers on both sides and get the number of nodes
		 *--- that belong to each FSI marker for each processor ---*/

		/*--- On the donor side ---*/

		for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++){
			/*--- If the tag GetMarker_All_FSIinterface(iMarkerDonor) equals the index we are looping at ---*/
			if ( donor_config->GetMarker_All_FSIinterface(iMarkerDonor) == iMarkerInt ){
				/*--- We have identified the local index of the Donor marker ---*/
				/*--- Now we are going to store the number of local points that belong to Marker_Donor on each processor ---*/
				/*--- This are the number of points that will be sent from this particular processor ---*/
				/*--- nLocalVertexDonorOwned WILL NOT include halo nodes ---*/
				/*--- nLocalVertexDonor WILL include halo nodes ---*/
				nLocalVertexDonorOwned = 0;
				for (iVertex = 0; iVertex < donor_geometry->GetnVertex(iMarkerDonor); iVertex++){
					Point_Donor = donor_geometry->vertex[iMarkerDonor][iVertex]->GetNode();
					if (donor_geometry->node[Point_Donor]->GetDomain()){
						nLocalVertexDonorOwned++;
					}
				}
				nLocalVertexDonor = donor_geometry->GetnVertex(iMarkerDonor);
				/*--- Store the identifier for the structural marker ---*/
				Marker_Donor = iMarkerDonor;
				/*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
				break;
			}
			else {
				/*--- If the tag hasn't matched any tag within the donor markers ---*/
				nLocalVertexDonor = 0;
				nLocalVertexDonorOwned = 0;
				Marker_Donor = -1;
			}
		}

		/*--- On the target side we only have to identify the marker; then we'll loop over it and retrieve from the fluid points ---*/

		for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++){
			/*--- If the tag GetMarker_All_FSIinterface(iMarkerFlow) equals the index we are looping at ---*/
			if ( target_config->GetMarker_All_FSIinterface(iMarkerTarget) == iMarkerInt ){
				/*--- Store the identifier for the fluid marker ---*/
				Marker_Target = iMarkerTarget;
				/*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
				break;
			}
			else {
				/*--- If the tag hasn't matched any tag within the Flow markers ---*/
				Marker_Target = -1;
			}
		}

		Buffer_Send_nVertexDonor[0] = nLocalVertexDonor;							   // Retrieve total number of vertices on Donor marker
		if (rank == MASTER_NODE) Buffer_Recv_nVertexDonor = new unsigned long[size];   // Allocate memory to receive how many vertices are on each rank on the structural side

#ifdef HAVE_MPI
		/*--- We receive MaxLocalVertexDonor as the maximum number of vertices in one single processor on the donor side---*/
		SU2_MPI::Allreduce(&nLocalVertexDonor, &MaxLocalVertexDonor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
		/*--- We receive TotalVertexDonorOwned as the total (real) number of vertices in one single interface marker on the donor side ---*/
		SU2_MPI::Allreduce(&nLocalVertexDonorOwned, &TotalVertexDonor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		/*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the structural side ---*/
		SU2_MPI::Gather(&Buffer_Send_nVertexDonor, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexDonor, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
		MaxLocalVertexDonor         = nLocalVertexDonor;
		TotalVertexDonor            = nLocalVertexDonorOwned;
		Buffer_Recv_nVertexDonor[0] = nLocalVertexDonor;
#endif

		/*--- We will be gathering the donor information into the master node ---*/
		nBuffer_DonorVariables = MaxLocalVertexDonor * nVar;
		nBuffer_DonorIndices = MaxLocalVertexDonor;

		/*--- Then we will broadcasting it to all the processors so they can retrieve the info they need ---*/
		/*--- We only broadcast those nodes that we need ---*/
		nBuffer_BcastVariables = TotalVertexDonor * nVar;
		nBuffer_BcastIndices = TotalVertexDonor;

		/*--- Send and Recv buffers ---*/

		/*--- Buffers to send and receive the variables in the donor mesh ---*/
		su2double *Buffer_Send_DonorVariables = new su2double[nBuffer_DonorVariables];
		su2double *Buffer_Recv_DonorVariables = NULL;

		/*--- Buffers to send and receive the indices in the donor mesh ---*/
		long *Buffer_Send_DonorIndices = new long[nBuffer_DonorIndices];
		long *Buffer_Recv_DonorIndices = NULL;

		/*--- Buffers to broadcast the variables and the indices ---*/
		su2double *Buffer_Bcast_Variables = new su2double[nBuffer_BcastVariables];
		long *Buffer_Bcast_Indices        = new long[nBuffer_BcastIndices];

		/*--- Prepare the receive buffers (1st step) and send buffers (2nd step) on the master node only. ---*/

		if (rank == MASTER_NODE) {
			Buffer_Recv_DonorVariables  = new su2double[size*nBuffer_DonorVariables];
			Buffer_Recv_DonorIndices    = new long[size*nBuffer_DonorIndices];
		}

		/*--- On the donor side ---*/
		/*--- First we initialize all of the indices and processors to -1 ---*/
		/*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
		for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++)
			Buffer_Send_DonorIndices[iVertex] = -1;

		if (Marker_Donor >= 0){

			for (iVertex = 0; iVertex < nLocalVertexDonor; iVertex++){

				Point_Donor = donor_geometry->vertex[Marker_Donor][iVertex]->GetNode();

				GetDonor_Variable(donor_solution, donor_geometry, donor_config, Marker_Donor, iVertex, Point_Donor);

				for (iVar = 0; iVar < nVar; iVar++){
					Buffer_Send_DonorVariables[iVertex*nVar+iVar] = Donor_Variable[iVar];
				}

				/*--- If this processor owns the node ---*/
				if (donor_geometry->node[Point_Donor]->GetDomain()){
					Point_Donor_Global = donor_geometry->node[Point_Donor]->GetGlobalIndex();
					Buffer_Send_DonorIndices[iVertex]     = Point_Donor_Global;
				}
				else{
					/*--- We set the values to be -1 to be able to identify them later as halo nodes ---*/
					Buffer_Send_DonorIndices[iVertex]     = -1;
				}

			}

		}

#ifdef HAVE_MPI
		/*--- Once all the messages have been prepared, we gather them all into the MASTER_NODE ---*/
		SU2_MPI::Gather(Buffer_Send_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, Buffer_Recv_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
		SU2_MPI::Gather(Buffer_Send_DonorIndices, nBuffer_DonorIndices, MPI_LONG, Buffer_Recv_DonorIndices, nBuffer_DonorIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);

#else
		for (iVariable = 0; iVariable < nBuffer_DonorVariables; iVariable++)
			Buffer_Recv_DonorVariables[iVariable] = Buffer_Send_DonorVariables[iVariable];
		for (iVariable = 0; iVariable < nBuffer_DonorIndices; iVariable++)
			Buffer_Recv_DonorIndices[iVariable] = Buffer_Send_DonorIndices[iVariable];
#endif

		/*--- Now we pack the information to send it over to the different processors ---*/

		if (rank == MASTER_NODE){

			/*--- For all the data we have received ---*/
			/*--- We initialize a counter to determine the position in the broadcast vector ---*/
			iLocalVertex = 0;

			for (iVertex = 0; iVertex < nProcessor*nBuffer_DonorIndices; iVertex++){

				/*--- If the donor index is not -1 (this is, if the node is not originally a halo node) ---*/
				if (Buffer_Recv_DonorIndices[iVertex] != -1){

					/*--- We set the donor index ---*/
					Buffer_Bcast_Indices[iLocalVertex] = Buffer_Recv_DonorIndices[iVertex];

					for (iVar = 0; iVar < nVar; iVar++){
						Buffer_Bcast_Variables[iLocalVertex*nVar+iVar] = Buffer_Recv_DonorVariables[iVertex*nVar + iVar];
					}

					iLocalVertex++;

				}

				if (iLocalVertex == TotalVertexDonor) break;

			}

		}

#ifdef HAVE_MPI
	SU2_MPI::Bcast(Buffer_Bcast_Variables, nBuffer_BcastVariables, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
	SU2_MPI::Bcast(Buffer_Bcast_Indices, nBuffer_BcastIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
#endif

		long indexPoint_iVertex, Point_Target_Check;

		/*--- For the target marker we are studying ---*/
		if (Marker_Target >= 0){

			/*--- We have identified the local index of the Structural marker ---*/
			/*--- We loop over all the vertices in that marker and in that particular processor ---*/

			for (iVertex = 0; iVertex < target_geometry->GetnVertex(Marker_Target); iVertex++){

				Point_Target = target_geometry->vertex[Marker_Target][iVertex]->GetNode();

				if (target_geometry->node[Point_Target]->GetDomain()){

					/*--- Find the global index of the donor point for Point_Target ---*/
					Donor_Global_Index = target_geometry->vertex[Marker_Target][iVertex]->GetGlobalDonorPoint();

					/*--- Find the index of the global donor point in the buffer Buffer_Bcast_Indices ---*/
					indexPoint_iVertex = std::distance(Buffer_Bcast_Indices, std::find(Buffer_Bcast_Indices, Buffer_Bcast_Indices + nBuffer_BcastIndices, Donor_Global_Index));

					Point_Target_Check = Buffer_Bcast_Indices[indexPoint_iVertex];

					if (Point_Target_Check < 0) {
						cout << "WARNING: A nonphysical point is being considered for traction transfer." << endl;
						exit(EXIT_FAILURE);
					}

					for (iVar = 0; iVar < nVar; iVar++)
						Target_Variable[iVar] = Buffer_Bcast_Variables[indexPoint_iVertex*nVar+iVar];

					SetTarget_Variable(target_solution, target_geometry, target_config, Marker_Target, iVertex, Point_Target);

				}

			}

		}

		delete [] Buffer_Send_DonorVariables;
		delete [] Buffer_Send_DonorIndices;
		delete [] Buffer_Bcast_Variables;
		delete [] Buffer_Bcast_Indices;

		if (rank == MASTER_NODE) {
			delete [] Buffer_Recv_nVertexDonor;
			delete [] Buffer_Recv_nVertexTarget;
			delete [] Buffer_Recv_nVertexDonorOwned;
			delete [] Buffer_Recv_DonorVariables;
			delete [] Buffer_Recv_DonorIndices;
		}


	}

}


void CTransfer::Allgather_InterfaceData_Matching(CSolver *donor_solution, CSolver *target_solution,
		   	   	   	   	   	   	                 CGeometry *donor_geometry, CGeometry *target_geometry,
												 CConfig *donor_config, CConfig *target_config){


}

void CTransfer::Scatter_InterfaceData_Interpolate(CSolver *donor_solution, CSolver *target_solution,
		   	   	   	   	   	   	   	   	   	   	  CGeometry *donor_geometry, CGeometry *target_geometry,
												  CConfig *donor_config, CConfig *target_config){

}

void CTransfer::Broadcast_InterfaceData_Interpolate(CSolver *donor_solution, CSolver *target_solution,
		   	   	   	   	   	   	   	   	   	   	    CGeometry *donor_geometry, CGeometry *target_geometry,
												    CConfig *donor_config, CConfig *target_config){


}

void CTransfer::Allgather_InterfaceData_Interpolate(CSolver *donor_solution, CSolver *target_solution,
		   	   	   	   	   	   	                    CGeometry *donor_geometry, CGeometry *target_geometry,
												    CConfig *donor_config, CConfig *target_config){


}


CTransfer_FlowTraction::CTransfer_FlowTraction(void) : CTransfer() {

}

CTransfer_FlowTraction::CTransfer_FlowTraction(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_FlowTraction::~CTransfer_FlowTraction(void) {

}

void CTransfer_FlowTraction::GetPhysical_Constants(CSolver *flow_solution, CSolver *struct_solution,
		   	   	   	   	   	   	   	   	   	   	   CGeometry *flow_geometry, CGeometry *struct_geometry,
												   CConfig *flow_config, CConfig *struct_config){

	unsigned short iVar;

	/*--- We have to clear the traction before applying it, because we are "adding" to node and not "setting" ---*/

	for (unsigned long iPoint = 0; iPoint < struct_geometry->GetnPoint(); iPoint++){
		struct_solution->node[iPoint]->Clear_FlowTraction();
	}

  	/*--- Redimensionalize the pressure ---*/

	su2double *Velocity_ND, *Velocity_Real;
	su2double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;
	su2double factorForces;

    Velocity_Real = flow_config->GetVelocity_FreeStream();
    Density_Real = flow_config->GetDensity_FreeStream();

    Velocity_ND = flow_config->GetVelocity_FreeStreamND();
    Density_ND = flow_config->GetDensity_FreeStreamND();

	Velocity2_Real = 0.0;
	Velocity2_ND = 0.0;
    for (iVar = 0; iVar < nVar; iVar++){
    	Velocity2_Real += Velocity_Real[iVar]*Velocity_Real[iVar];
    	Velocity2_ND += Velocity_ND[iVar]*Velocity_ND[iVar];
    }

    Physical_Constants[0] = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);

	/*--- Apply a ramp to the transfer of the fluid loads ---*/

	su2double ModAmpl;
	su2double CurrentTime = struct_config->GetCurrent_DynTime();
	su2double Static_Time = struct_config->GetStatic_Time();

	bool Ramp_Load = struct_config->GetRamp_Load();
	su2double Ramp_Time = struct_config->GetRamp_Time();

	if (CurrentTime <= Static_Time){ ModAmpl=0.0; }
	else if((CurrentTime > Static_Time) &&
			(CurrentTime <= (Static_Time + Ramp_Time)) &&
			(Ramp_Load)){
		ModAmpl = (CurrentTime-Static_Time) / Ramp_Time;
		ModAmpl = max(ModAmpl,0.0);
		ModAmpl = min(ModAmpl,1.0);
		Physical_Constants[1] = ModAmpl;
	}
	else{ Physical_Constants[1] = 1.0; }

}

void CTransfer_FlowTraction::GetDonor_Variable(CSolver *flow_solution, CGeometry *flow_geometry, CConfig *flow_config,
					   	   	   	   	   	   	   unsigned long Marker_Flow, unsigned long Vertex_Flow, unsigned long Point_Struct){


	unsigned short iVar, jVar;
	unsigned long Point_Flow;
	su2double *Normal_Flow;

	// Check the kind of fluid problem
	bool compressible       = (flow_config->GetKind_Regime() == COMPRESSIBLE);
	bool incompressible     = (flow_config->GetKind_Regime() == INCOMPRESSIBLE);
	bool viscous_flow       = ((flow_config->GetKind_Solver() == NAVIER_STOKES) ||
							   (flow_config->GetKind_Solver() == RANS) );

	// Parameters for the calculations
	// Pn: Pressure
	// Pinf: Pressure_infinite
	// div_vel: Velocity divergence
	// Dij: Dirac delta
	su2double Pn = 0.0, div_vel = 0.0, Dij = 0.0;
	su2double Viscosity = 0.0, Density = 0.0;
	su2double **Grad_PrimVar;
	su2double Tau[3][3] = { {0.0, 0.0, 0.0} ,
							{0.0, 0.0, 0.0} ,
							{0.0, 0.0, 0.0} } ;

	su2double Pinf = flow_solution->GetPressure_Inf();

	Point_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNode();
		// Get the normal at the vertex: this normal goes inside the fluid domain.
	Normal_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNormal();

	// Retrieve the values of pressure, viscosity and density
	if (incompressible){

		Pn = flow_solution->node[Point_Flow]->GetPressureInc();

		if (viscous_flow){

			Grad_PrimVar = flow_solution->node[Point_Flow]->GetGradient_Primitive();
			Viscosity = flow_solution->node[Point_Flow]->GetLaminarViscosityInc();
			Density = flow_solution->node[Point_Flow]->GetDensityInc();

		}
	}
	else if (compressible){

		Pn = flow_solution->node[Point_Flow]->GetPressure();

		if (viscous_flow){

			Grad_PrimVar = flow_solution->node[Point_Flow]->GetGradient_Primitive();
			Viscosity = flow_solution->node[Point_Flow]->GetLaminarViscosity();
			Density = flow_solution->node[Point_Flow]->GetDensity();

		}
	}

	// Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
	for (iVar = 0; iVar < nVar; iVar++) {
		Donor_Variable[iVar] = -(Pn-Pinf)*Normal_Flow[iVar];
	}

	// Calculate tn in the fluid nodes for the viscous term

	if (viscous_flow){

		// Divergence of the velocity
		div_vel = 0.0; for (iVar = 0; iVar < nVar; iVar++) div_vel += Grad_PrimVar[iVar+1][iVar];
		if (incompressible) div_vel = 0.0;

		for (iVar = 0; iVar < nVar; iVar++) {

			for (jVar = 0 ; jVar < nVar; jVar++) {
				// Dirac delta
				Dij = 0.0; if (iVar == jVar) Dij = 1.0;

				// Viscous stress
				Tau[iVar][jVar] = Viscosity*(Grad_PrimVar[jVar+1][iVar] + Grad_PrimVar[iVar+1][jVar]) -
						TWO3*Viscosity*div_vel*Dij;

				// Viscous component in the tn vector --> Units of force (non-dimensional).
				Donor_Variable[iVar] += Tau[iVar][jVar]*Normal_Flow[jVar];
			}
		}
	}

	// Redimensionalize and take into account ramp transfer of the loads
	for (iVar = 0; iVar < nVar; iVar++){
		Donor_Variable[iVar] = Donor_Variable[iVar] * Physical_Constants[0] * Physical_Constants[1];
	}

}

void CTransfer_FlowTraction::SetTarget_Variable(CSolver *fea_solution, CGeometry *fea_geometry,
												CConfig *fea_config, unsigned long Marker_Struct,
												unsigned long Vertex_Struct, unsigned long Point_Struct){

	/*--- Add to the Flow traction ---*/
	fea_solution->node[Point_Struct]->Add_FlowTraction(Target_Variable);

}

CTransfer_StructuralDisplacements::CTransfer_StructuralDisplacements(void) : CTransfer() {

}

CTransfer_StructuralDisplacements::CTransfer_StructuralDisplacements(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_StructuralDisplacements::~CTransfer_StructuralDisplacements(void) {

}


void CTransfer_StructuralDisplacements::GetPhysical_Constants(CSolver *struct_solution, CSolver *flow_solution,
		   	   	   	   	   	   	   	   	   	   	   CGeometry *struct_geometry, CGeometry *flow_geometry,
												   CConfig *struct_config, CConfig *flow_config){

}

void CTransfer_StructuralDisplacements::GetDonor_Variable(CSolver *struct_solution, CGeometry *struct_geometry, CConfig *struct_config,
					   	   	   	   	   	   	   	          unsigned long Marker_Struct, unsigned long Vertex_Struct, unsigned long Point_Struct){


	su2double *Coord_Struct, *Displacement_Struct;
	unsigned short iVar;

    Coord_Struct = struct_geometry->node[Point_Struct]->GetCoord();

    /*--- The displacements come from the predicted solution ---*/
    Displacement_Struct = struct_solution->node[Point_Struct]->GetSolution_Pred();

//    cout << "For point " << Point_Struct << " we have coordinates " << Coord_Struct[0] << " " << Coord_Struct[1] << endl;
//    cout << "and displacements " << Displacement_Struct[0] << " " << Displacement_Struct[1] << endl;


	for (iVar = 0; iVar < nVar; iVar++){
		Donor_Variable[iVar] = Coord_Struct[iVar] + Displacement_Struct[iVar];
	}

}

void CTransfer_StructuralDisplacements::SetTarget_Variable(CSolver *flow_solution, CGeometry *flow_geometry,
														   CConfig *flow_config, unsigned long Marker_Flow,
														   unsigned long Vertex_Flow, unsigned long Point_Flow){

	su2double *Coord, VarCoord[3] = {0.0, 0.0, 0.0};
	unsigned short iVar;

	Coord = flow_geometry->node[Point_Flow]->GetCoord();

	for (iVar = 0; iVar < nVar; iVar++)
		VarCoord[iVar] = Target_Variable[iVar]-Coord[iVar];

	flow_geometry->vertex[Marker_Flow][Vertex_Flow]->SetVarCoord(VarCoord);

}
