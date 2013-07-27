/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
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

#include "../include/output_structure.hpp"


COutput::COutput(void) {

	/*--- Initialize point and connectivity counters to zero. ---*/
	nGlobal_Poin      = 0;
  nSurf_Poin        = 0;
	nGlobal_Elem      = 0;
	nSurf_Elem        = 0;
	nGlobal_Tria      = 0;
	nGlobal_Quad      = 0;
	nGlobal_Tetr      = 0;
	nGlobal_Hexa      = 0;
	nGlobal_Wedg      = 0;
	nGlobal_Pyra      = 0;
	nGlobal_Line      = 0;
	nGlobal_BoundTria = 0;
	nGlobal_BoundQuad = 0;

	/*--- Initialize CGNS write flag ---*/
	wrote_base_file = false;

	/*--- Initialize CGNS write flag ---*/
	wrote_CGNS_base = false;

	/*--- Initialize Tecplot write flag ---*/
	wrote_Tecplot_base = false;

}

COutput::~COutput(void) { }

void COutput::SetSurfaceCSV_Flow(CConfig *config, CGeometry *geometry, CSolver *FlowSolver, unsigned long iExtIter) {

#ifdef NO_MPI
	unsigned long iPoint, iVertex, Global_Index;
	unsigned short iMarker;
	double PressCoeff = 0.0, SkinFrictionCoeff, xCoord, yCoord, zCoord, Mach;
	char cstr[200], buffer [50];
	ofstream SurfFlow_file;
	unsigned short solver = config->GetKind_Solver();

	/*--- Write file name with extension if unsteady ---*/
	strcpy (cstr, config->GetSurfFlowCoeff_FileName().c_str());
	if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.csv", int(iExtIter));
		if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.csv", int(iExtIter));
		if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.csv", int(iExtIter));
		if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv", int(iExtIter));
		if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
	}
	else
		sprintf (buffer, ".csv");

	strcat (cstr, buffer);
	SurfFlow_file.precision(15);
	SurfFlow_file.open(cstr, ios::out);

	if (geometry->GetnDim() == 2) {
		switch (solver) {
		case EULER : case FREE_SURFACE_EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"Global_Index\"" << endl; break;
		case NAVIER_STOKES: case RANS: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\",\"Global_Index\"" << endl; break;
		}
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);
					PressCoeff = FlowSolver->GetCPressure(iMarker,iVertex);
					Global_Index = geometry->node[iPoint]->GetGlobalIndex();
					switch (solver) {
					case EULER : case FREE_SURFACE_EULER :
						Mach = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << Mach << "," << yCoord << ", " << Global_Index << endl;
						break;
					case NAVIER_STOKES: case RANS: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
						SkinFrictionCoeff = FlowSolver->GetCSkinFriction(iMarker,iVertex);
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << SkinFrictionCoeff << "," << yCoord << ", " << Global_Index << endl;
						break;
					}
				}
	}

	if (geometry->GetnDim() == 3) {
		switch (solver) {
		case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"z_coord\",\"Global_Index\"" << endl; break;
		case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\",\"z_coord\",\"Global_Index\"" << endl; break;
		}
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);
					zCoord = geometry->node[iPoint]->GetCoord(2);
					PressCoeff = FlowSolver->GetCPressure(iMarker,iVertex);
					Global_Index = geometry->node[iPoint]->GetGlobalIndex();
					switch (solver) {
					case EULER : case FREE_SURFACE_EULER :
						Mach = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << Mach <<"," << yCoord << "," << zCoord << ", " << Global_Index << endl;
						break;
					case NAVIER_STOKES: case RANS: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
						SkinFrictionCoeff = FlowSolver->GetCSkinFriction(iMarker,iVertex);
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << SkinFrictionCoeff << "," << yCoord << "," << zCoord << ", " << Global_Index << endl;
						break;
					}
				}
	}

	SurfFlow_file.close();

#else

	int rank = MPI::COMM_WORLD.Get_rank(), iProcessor, nProcessor = MPI::COMM_WORLD.Get_size();
	double PressCoeff = 0.0, SkinFrictionCoeff, xCoord, yCoord, zCoord, Mach;
	char cstr[200];
	unsigned short iMarker;
	unsigned long Buffer_Send_nVertex[1], iVertex, iPoint, nVertex_Surface = 0, nLocalVertex_Surface = 0,
			MaxLocalVertex_Surface = 0, nBuffer_Scalar, *Buffer_Receive_nVertex = NULL, position, Global_Index;
	ofstream SurfFlow_file;

	/*--- Write the surface .csv file, the information of all the vertices is
     sent to the MASTER_NODE for writing ---*/
	nLocalVertex_Surface = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Plotting(iMarker) == YES)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface++;
			}

	if (rank == MASTER_NODE)
		Buffer_Receive_nVertex = new unsigned long [nProcessor];

	Buffer_Send_nVertex[0] = nLocalVertex_Surface;

	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG, MASTER_NODE);

	double *Buffer_Send_Coord_x = new double [MaxLocalVertex_Surface];
	double *Buffer_Send_Coord_y = new double [MaxLocalVertex_Surface];
	double *Buffer_Send_Coord_z = new double [MaxLocalVertex_Surface];
	double *Buffer_Send_Press = new double [MaxLocalVertex_Surface];
	double *Buffer_Send_Mach = new double [MaxLocalVertex_Surface];
	double *Buffer_Send_SkinFriction = new double [MaxLocalVertex_Surface];
	unsigned long *Buffer_Send_GlobalIndex = new unsigned long [MaxLocalVertex_Surface];

	nVertex_Surface = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Plotting(iMarker) == YES)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					Buffer_Send_Press[nVertex_Surface] = FlowSolver->GetCPressure(iMarker,iVertex);
					Buffer_Send_Coord_x[nVertex_Surface] = geometry->node[iPoint]->GetCoord(0);
					Buffer_Send_Coord_y[nVertex_Surface] = geometry->node[iPoint]->GetCoord(1);
					Buffer_Send_GlobalIndex[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
					if (geometry->GetnDim() == 3) Buffer_Send_Coord_z[nVertex_Surface] = geometry->node[iPoint]->GetCoord(2);
					if (config->GetKind_Solver() == EULER)
						Buffer_Send_Mach[nVertex_Surface] = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
					if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS))
						Buffer_Send_SkinFriction[nVertex_Surface] = FlowSolver->GetCSkinFriction(iMarker,iVertex);
					nVertex_Surface++;
				}
			}

	double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Press = NULL,
			*Buffer_Receive_Mach = NULL, *Buffer_Receive_SkinFriction = NULL;
	unsigned long *Buffer_Receive_GlobalIndex = NULL;

	if (rank == MASTER_NODE) {
		Buffer_Receive_Coord_x = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Coord_y = new double [nProcessor*MaxLocalVertex_Surface];
		if (geometry->GetnDim() == 3)
			Buffer_Receive_Coord_z = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Press = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Mach = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_SkinFriction = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_GlobalIndex = new  unsigned long [nProcessor*MaxLocalVertex_Surface];
	}

	nBuffer_Scalar = MaxLocalVertex_Surface;

	/*--- Send the information to the Master node ---*/
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Gather(Buffer_Send_Coord_x, nBuffer_Scalar, MPI::DOUBLE,
			Buffer_Receive_Coord_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Coord_y, nBuffer_Scalar, MPI::DOUBLE,
			Buffer_Receive_Coord_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (geometry->GetnDim() == 3)
		MPI::COMM_WORLD.Gather(Buffer_Send_Coord_z, nBuffer_Scalar, MPI::DOUBLE,
				Buffer_Receive_Coord_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Press, nBuffer_Scalar, MPI::DOUBLE,
			Buffer_Receive_Press, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if ((config->GetKind_Solver() == EULER) || (config->GetKind_Solver() == FREE_SURFACE_EULER))
		MPI::COMM_WORLD.Gather(Buffer_Send_Mach, nBuffer_Scalar, MPI::DOUBLE,
				Buffer_Receive_Mach, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS) ||
      (config->GetKind_Solver() == FREE_SURFACE_NAVIER_STOKES) || (config->GetKind_Solver() == FREE_SURFACE_RANS))
		MPI::COMM_WORLD.Gather(Buffer_Send_SkinFriction, nBuffer_Scalar, MPI::DOUBLE,
				Buffer_Receive_SkinFriction, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);

	MPI::COMM_WORLD.Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
			Buffer_Receive_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG, MASTER_NODE);

	/*--- The master node writes the surface files ---*/
	if (rank == MASTER_NODE) {

		/*--- Write file name with extension if unsteady ---*/
		char buffer[50];
		string filename = config->GetSurfFlowCoeff_FileName();

		/*--- Remove the domain number from the surface csv filename ---*/
		if (nProcessor > 1) filename.erase (filename.end()-2, filename.end());

		/*--- Write file name with extension if unsteady ---*/
		strcpy (cstr, filename.c_str());
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
			if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv", int(iExtIter));
			if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
		}
		else
			sprintf (buffer, ".csv");

		strcat (cstr, buffer);
		SurfFlow_file.precision(15);
		SurfFlow_file.open(cstr, ios::out);

		/*--- Write the 2D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 2) {
			switch (config->GetKind_Solver()) {
			case EULER : case FREE_SURFACE_EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"Global_Index\"" << endl; break;
			case NAVIER_STOKES: case RANS: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\",\"Global_Index\"" << endl; break;
			}
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					xCoord = Buffer_Receive_Coord_x[position];
					yCoord = Buffer_Receive_Coord_y[position];
					PressCoeff = Buffer_Receive_Press[position];
					Global_Index = Buffer_Receive_GlobalIndex[position];
					switch (config->GetKind_Solver()) {
					case EULER : case FREE_SURFACE_EULER :
						Mach = Buffer_Receive_Mach[position];
						SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << Mach << ", " << yCoord << ", " << Global_Index << endl;
						break;
					case NAVIER_STOKES: case RANS: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
						SkinFrictionCoeff = Buffer_Receive_SkinFriction[position];
						SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << SkinFrictionCoeff << ", " << yCoord << ", " << Global_Index << endl;
						break;
					}
				}
		}

		/*--- Write the 3D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 3) {
			switch (config->GetKind_Solver()) {
			case EULER : case FREE_SURFACE_EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"z_coord\",\"Global_Index\"" << endl; break;
			case NAVIER_STOKES: case RANS: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\",\"z_coord\",\"Global_Index\"" << endl; break;
			}
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					xCoord = Buffer_Receive_Coord_x[position];
					yCoord = Buffer_Receive_Coord_y[position];
					zCoord = Buffer_Receive_Coord_z[position];
					PressCoeff = Buffer_Receive_Press[position];
					Global_Index = Buffer_Receive_GlobalIndex[position];
					switch (config->GetKind_Solver()) {
					case EULER : case FREE_SURFACE_EULER :
						Mach = Buffer_Receive_Mach[position];
						SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << Mach << ", " << yCoord << ", " << zCoord << ", " << Global_Index << endl;
						break;
					case NAVIER_STOKES: case RANS: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
						SkinFrictionCoeff = Buffer_Receive_SkinFriction[position];
						SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << SkinFrictionCoeff << ", " << yCoord << ", " << zCoord << ", " << Global_Index << endl;
						break;
					}
				}
		}
	}

	if (rank == MASTER_NODE) {
		delete [] Buffer_Receive_Coord_x;
		delete [] Buffer_Receive_Coord_y;
		if (geometry->GetnDim() == 3)
			delete [] Buffer_Receive_Coord_z;
		delete [] Buffer_Receive_Press;
		delete [] Buffer_Receive_Mach;
		delete [] Buffer_Receive_SkinFriction;
		delete [] Buffer_Receive_GlobalIndex;
	}

	delete [] Buffer_Send_Coord_x;
	delete [] Buffer_Send_Coord_y;
	delete [] Buffer_Send_Coord_z;
	delete [] Buffer_Send_Press;
	delete [] Buffer_Send_Mach;
	delete [] Buffer_Send_SkinFriction;
	delete [] Buffer_Send_GlobalIndex;

	SurfFlow_file.close();

#endif
}

void COutput::SetSurfaceCSV_Adjoint(CConfig *config, CGeometry *geometry, CSolver *AdjSolver, CSolver *FlowSolution, unsigned long iExtIter, unsigned short val_iZone) {

#ifdef NO_MPI

	unsigned long iPoint, iVertex;
	double *Solution, xCoord, yCoord, zCoord, *IntBoundary_Jump;
	unsigned short iMarker;
	char cstr[200], buffer[50];
	ofstream SurfAdj_file;

	/*--- Write file name with extension if unsteady ---*/
	strcpy (cstr, config->GetSurfAdjCoeff_FileName().c_str());

	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
		if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.csv", int(val_iZone));
		if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.csv", int(val_iZone));
		if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
		if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));

	} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
		if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
		if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
		if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
		if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
	}
	else
		sprintf (buffer, ".csv");

	strcat(cstr, buffer);
	SurfAdj_file.precision(15);
	SurfAdj_file.open(cstr, ios::out);

	if (geometry->GetnDim() == 2) {
		SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					Solution = AdjSolver->node[iPoint]->GetSolution();
					IntBoundary_Jump = AdjSolver->node[iPoint]->GetIntBoundary_Jump();
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);
					SurfAdj_file << scientific << iPoint << ", " << AdjSolver->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
							<< Solution[1] << ", " << Solution[2] << ", " << Solution[3] <<", " << xCoord <<", "<< yCoord << endl;
				}
		}
	}

	if (geometry->GetnDim() == 3) {
		SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					Solution = AdjSolver->node[iPoint]->GetSolution();
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);
					zCoord = geometry->node[iPoint]->GetCoord(2);

					SurfAdj_file << scientific << iPoint << ", " << AdjSolver->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
							<< Solution[1] << ", " << Solution[2] << ", " << Solution[3] << ", " << Solution[4] << ", "<< xCoord <<", "<< yCoord <<", "<< zCoord << endl;
				}
		}
	}

	SurfAdj_file.close();

#else

	int rank = MPI::COMM_WORLD.Get_rank(), iProcessor, nProcessor = MPI::COMM_WORLD.Get_size();
	unsigned short nDim = geometry->GetnDim(), iMarker;
	double *Solution, *Normal, *d, *Coord;
	unsigned long Buffer_Send_nVertex[1], iVertex, iPoint, nVertex_Surface = 0, nLocalVertex_Surface = 0,
			MaxLocalVertex_Surface = 0, nBuffer_Scalar;
	unsigned long *Buffer_Receive_nVertex = NULL;
	ofstream SurfAdj_file;

	/*--- Write the surface .csv file ---*/
	nLocalVertex_Surface = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Plotting(iMarker) == YES)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface ++;
			}

	if (rank == MASTER_NODE)
		Buffer_Receive_nVertex = new unsigned long [nProcessor];

	Buffer_Send_nVertex[0] = nLocalVertex_Surface;

	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG, MASTER_NODE);

	double *Buffer_Send_Coord_x = new double[MaxLocalVertex_Surface];
	double *Buffer_Send_Coord_y= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_Coord_z= new double[MaxLocalVertex_Surface];
	unsigned long *Buffer_Send_GlobalPoint= new unsigned long[MaxLocalVertex_Surface];
	double *Buffer_Send_Sensitivity= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_PsiRho= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_Phi_x= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_Phi_y= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_Phi_z= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_PsiE= new double[MaxLocalVertex_Surface];

	nVertex_Surface = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Plotting(iMarker) == YES)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					Solution = AdjSolver->node[iPoint]->GetSolution();
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Coord = geometry->node[iPoint]->GetCoord();
					d = AdjSolver->node[iPoint]->GetForceProj_Vector();
					Buffer_Send_GlobalPoint[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
					Buffer_Send_Coord_x[nVertex_Surface] = Coord[0];
					Buffer_Send_Coord_y[nVertex_Surface] = Coord[1];
					Buffer_Send_Sensitivity[nVertex_Surface] =  AdjSolver->GetCSensitivity(iMarker,iVertex);
					Buffer_Send_PsiRho[nVertex_Surface] = Solution[0];
					Buffer_Send_Phi_x[nVertex_Surface] = Solution[1];
					Buffer_Send_Phi_y[nVertex_Surface] = Solution[2];
					if (nDim == 2) Buffer_Send_PsiE[nVertex_Surface] = Solution[3];
					if (nDim == 3) {
						Buffer_Send_Coord_z[nVertex_Surface] = Coord[2];
						Buffer_Send_Phi_z[nVertex_Surface] = Solution[3];
						Buffer_Send_PsiE[nVertex_Surface] = Solution[4];
					}
					nVertex_Surface++;
				}
			}

	double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Sensitivity = NULL,
			*Buffer_Receive_PsiRho = NULL, *Buffer_Receive_Phi_x = NULL, *Buffer_Receive_Phi_y = NULL, *Buffer_Receive_Phi_z = NULL,
			*Buffer_Receive_PsiE = NULL;
	unsigned long *Buffer_Receive_GlobalPoint = NULL;

	if (rank == MASTER_NODE) {
		Buffer_Receive_Coord_x = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Coord_y = new double [nProcessor*MaxLocalVertex_Surface];
		if (nDim == 3) Buffer_Receive_Coord_z = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Sensitivity = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_PsiRho = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Phi_x = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Phi_y = new double [nProcessor*MaxLocalVertex_Surface];
		if (nDim == 3) Buffer_Receive_Phi_z = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_PsiE = new double [nProcessor*MaxLocalVertex_Surface];
	}

	nBuffer_Scalar = MaxLocalVertex_Surface;

	/*--- Send the information to the Master node ---*/
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Gather(Buffer_Send_Coord_x, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Coord_y, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (nDim == 3) MPI::COMM_WORLD.Gather(Buffer_Send_Coord_z, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_GlobalPoint, nBuffer_Scalar, MPI::UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Scalar, MPI::UNSIGNED_LONG, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Sensitivity, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Sensitivity, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_PsiRho, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_PsiRho, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Phi_x, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Phi_y, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (nDim == 3) MPI::COMM_WORLD.Gather(Buffer_Send_Phi_z, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_PsiE, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_PsiE, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);

	/*--- The master node is the one who writes the surface files ---*/
	if (rank == MASTER_NODE) {
		unsigned long iVertex, GlobalPoint, position;
		char cstr[200], buffer[50];
		ofstream SurfAdj_file;
		string filename = config->GetSurfAdjCoeff_FileName();

		/*--- Remove the domain number from the surface csv filename ---*/
		if (nProcessor > 1) filename.erase (filename.end()-2, filename.end());

		/*--- Write file name with extension if unsteady ---*/
		strcpy (cstr, filename.c_str());

		if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
			if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
			if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.csv", int(val_iZone));
			if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.csv", int(val_iZone));
			if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
			if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));

		} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
			if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv", int(iExtIter));
			if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
		}
		else
			sprintf (buffer, ".csv");

		strcat (cstr, buffer);
		SurfAdj_file.open(cstr, ios::out);
		SurfAdj_file.precision(15);

		/*--- Write the 2D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 2) {

			SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;

			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {

					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					GlobalPoint = Buffer_Receive_GlobalPoint[position];

					SurfAdj_file << scientific << GlobalPoint <<
							", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
							", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] <<
							", " << Buffer_Receive_PsiE[position] << ", " << Buffer_Receive_Coord_x[position] <<
							", "<< Buffer_Receive_Coord_y[position]  << endl;
				}
		}

		/*--- Write the 3D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 3) {

			SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;

			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					GlobalPoint = Buffer_Receive_GlobalPoint[position];

					SurfAdj_file << scientific << GlobalPoint <<
							", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
							", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] << ", " << Buffer_Receive_Phi_z[position] <<
							", " << Buffer_Receive_PsiE[position] <<", "<< Buffer_Receive_Coord_x[position] <<
							", "<< Buffer_Receive_Coord_y[position] <<", "<< Buffer_Receive_Coord_z[position] << endl;
				}
		}

	}

	if (rank == MASTER_NODE) {
		delete [] Buffer_Receive_nVertex;
		delete [] Buffer_Receive_Coord_x;
		delete [] Buffer_Receive_Coord_y;
		if (nDim == 3) delete [] Buffer_Receive_Coord_z;
		delete [] Buffer_Receive_Sensitivity;
		delete [] Buffer_Receive_PsiRho;
		delete [] Buffer_Receive_Phi_x;
		delete [] Buffer_Receive_Phi_y;
		if (nDim == 3) delete [] Buffer_Receive_Phi_z;
		delete [] Buffer_Receive_PsiE;
		delete [] Buffer_Receive_GlobalPoint;
	}

	delete [] Buffer_Send_Coord_x;
	delete [] Buffer_Send_Coord_y;
	delete [] Buffer_Send_Coord_z;
	delete [] Buffer_Send_GlobalPoint;
	delete [] Buffer_Send_Sensitivity;
	delete [] Buffer_Send_PsiRho;
	delete [] Buffer_Send_Phi_x;
	delete [] Buffer_Send_Phi_y;
	delete [] Buffer_Send_Phi_z;
	delete [] Buffer_Send_PsiE;

	SurfAdj_file.close();

#endif
}

void COutput::SetSurfaceCSV_Linearized(CConfig *config, CGeometry *geometry, CSolver *LinSolution, string val_filename, unsigned long iExtIter) { }

void COutput::MergeGeometry(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {

	int rank = MASTER_NODE;
  int size = 1;

#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
#endif

	/*--- Get the file output format ---*/
  
	unsigned short FileFormat = config->GetOutput_FileFormat();

	bool dynamic_mesh = (config->GetUnsteady_Simulation() &&
			config->GetGrid_Movement());

	/*--- Merge connectivity for each type of element (excluding halos). Note
     that we only need to merge the connectivity once, as it does not change
     during computation. Check whether the base file has been written. ---*/

	if (!wrote_base_file) {

    /*--- Merge volumetric grid. ---*/

		MergeVolumetricConnectivity(config, geometry, TRIANGLE    );
    if ((rank == MASTER_NODE) && (size != 1) && (nGlobal_Tria != 0)) cout <<"Merging volumetric triangle grid connectivity." << endl;

		MergeVolumetricConnectivity(config, geometry, RECTANGLE   );
    if ((rank == MASTER_NODE) && (size != 1) && (nGlobal_Quad != 0)) cout <<"Merging volumetric rectangle grid connectivity." << endl;

		MergeVolumetricConnectivity(config, geometry, TETRAHEDRON );
    if ((rank == MASTER_NODE) && (size != 1) && (nGlobal_Tetr != 0)) cout <<"Merging volumetric tetrahedron grid connectivity." << endl;

		MergeVolumetricConnectivity(config, geometry, HEXAHEDRON  );
    if ((rank == MASTER_NODE) && (size != 1) && (nGlobal_Hexa != 0)) cout <<"Merging volumetric hexahedron grid connectivity." << endl;

		MergeVolumetricConnectivity(config, geometry, WEDGE       );
    if ((rank == MASTER_NODE) && (size != 1) && (nGlobal_Wedg != 0)) cout <<"Merging volumetric wedge grid connectivity." << endl;

		MergeVolumetricConnectivity(config, geometry, PYRAMID     );
    if ((rank == MASTER_NODE) && (size != 1) && (nGlobal_Pyra != 0)) cout <<"Merging volumetric pyramid grid connectivity." << endl;

    /*--- Merge surface grid. ---*/
    
    MergeSurfaceConnectivity(config, geometry, LINE      );
    if ((rank == MASTER_NODE) && (size != 1) && (nGlobal_Line != 0)) cout <<"Merging surface line grid connectivity." << endl;

		MergeSurfaceConnectivity(config, geometry, TRIANGLE  );
    if ((rank == MASTER_NODE) && (size != 1) && (nGlobal_BoundTria != 0)) cout <<"Merging surface triangle grid connectivity." << endl;

		MergeSurfaceConnectivity(config, geometry, RECTANGLE );
    if ((rank == MASTER_NODE) && (size != 1) && (nGlobal_BoundQuad != 0)) cout <<"Merging surface rectangle grid connectivity." << endl;


		/*--- Update total number of volumetric elements after merge. ---*/

		nGlobal_Elem = nGlobal_Tria + nGlobal_Quad + nGlobal_Tetr +
				nGlobal_Hexa + nGlobal_Pyra + nGlobal_Wedg;
    
    /*--- Update total number of surface elements after merge. ---*/
    
    nSurf_Elem = nGlobal_Line + nGlobal_BoundTria + nGlobal_BoundQuad;

		/*--- Write the connectivity to the base binary output file, then
         clear the memory immediately for the rest of the computation. ---*/

		if (rank == MASTER_NODE && FileFormat == CGNS_SOL) {
			SetCGNS_Connectivity(config, geometry, val_iZone);
			DeallocateConnectivity(config, geometry, false);
		}

	}

	/*--- Merge coordinates of all grid nodes (excluding ghost points).
     Note that this also only needs to be done once, unless it is an
     unsteady simulation with grid motion. ---*/

	if (!wrote_base_file || dynamic_mesh) {
    if ((rank == MASTER_NODE) && (size != 1)) cout <<"Merging grid coordinates." << endl;
		MergeCoordinates(config, geometry);
  }

	if (rank == MASTER_NODE) {

		if (FileFormat == CGNS_SOL) {
			SetCGNS_Coordinates(config, geometry, val_iZone);
			if (!wrote_base_file || dynamic_mesh) DeallocateCoordinates(config, geometry);
		} else if (FileFormat == TECPLOT_BINARY) {
			SetTecplot_Mesh(config, geometry, val_iZone);
			if (!wrote_base_file) DeallocateConnectivity(config, geometry, false);
		}
	}

}

void COutput::MergeCoordinates(CConfig *config, CGeometry *geometry) {

	/*--- Local variables needed on all processors ---*/

	unsigned short iDim, nDim = geometry->GetnDim();
	unsigned long iPoint, jPoint;
  
#ifdef NO_MPI

	/*--- In serial, the single process has access to all geometry, so simply
     load the coordinates into the data structure. ---*/

	/*--- Total number of points in the mesh (excluding halos). ---*/
  nGlobal_Poin = geometry->GetnPointDomain();
	nGlobal_Doma = geometry->GetnPointDomain();

	/*--- Allocate the coordinates data structure. ---*/

	Coords = new double*[nDim];
	for (iDim = 0; iDim < nDim; iDim++) {
		Coords[iDim] = new double[nGlobal_Poin];
	}

	/*--- Loop over the mesh to collect the coords of the local points. ---*/

	jPoint = 0;
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Retrieve the current coordinates at this node. ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        Coords[iDim][jPoint] = geometry->node[iPoint]->GetCoord(iDim);
      }
      
      /*--- Increment a counter since we may be skipping over
       some halo nodes during this loop. ---*/
      jPoint++;
		}
	}

#else

	/*--- MPI preprocessing ---*/
	int iProcessor;
	int nProcessor = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();

  bool Wrt_Halo = config->GetWrt_Halo();

	/*--- Local variables needed for merging the geometry with MPI. ---*/

	unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
	unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
	unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;

	if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];

	/*--- Sum total number of nodes that belong to the domain ---*/

	Buffer_Send_nPoin[0] = geometry->GetnPointDomain();
	MPI::COMM_WORLD.Gather(&Buffer_Send_nPoin, 1, MPI::UNSIGNED_LONG,
			Buffer_Recv_nPoin, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
	if (rank == MASTER_NODE) {
		nGlobal_Doma = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
			nGlobal_Doma += Buffer_Recv_nPoin[iProcessor];
		}
	}

	/*--- Each processor sends its local number of nodes to the master. ---*/

  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else
  nLocalPoint = geometry->GetnPointDomain();
	Buffer_Send_nPoin[0] = nLocalPoint;
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalPoint, &MaxLocalPoint,
			1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nPoin, 1, MPI::UNSIGNED_LONG,
			Buffer_Recv_nPoin, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
	nBuffer_Scalar = MaxLocalPoint;

	/*--- Send and Recv buffers. ---*/

	double *Buffer_Send_X = new double[MaxLocalPoint];
	double *Buffer_Recv_X = NULL;

	double *Buffer_Send_Y = new double[MaxLocalPoint];
	double *Buffer_Recv_Y = NULL;

	double *Buffer_Send_Z, *Buffer_Recv_Z = NULL;
	if (nDim == 3) Buffer_Send_Z = new double[MaxLocalPoint];

	unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
	unsigned long *Buffer_Recv_GlobalIndex = NULL;

	/*--- Prepare the receive buffers in the master node only. ---*/

	if (rank == MASTER_NODE) {

		Buffer_Recv_X = new double[nProcessor*MaxLocalPoint];
		Buffer_Recv_Y = new double[nProcessor*MaxLocalPoint];
		if (nDim == 3) Buffer_Recv_Z = new double[nProcessor*MaxLocalPoint];
		Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];

		/*--- Sum total number of nodes to be written and allocate arrays ---*/
		nGlobal_Poin = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
			nGlobal_Poin += Buffer_Recv_nPoin[iProcessor];
		}
		Coords = new double*[nDim];
		for (iDim = 0; iDim < nDim; iDim++) {
			Coords[iDim] = new double[nGlobal_Poin];
		}
	}

	/*--- Main communication routine. Loop over each coordinate and perform
     the MPI comm. Temporary 1-D buffers are used to send the coordinates at
     all nodes on each partition to the master node. These are then unpacked
     by the master and sorted by global index in one large n-dim. array. ---*/

	/*--- Loop over this partition to collect the coords of the local points.
     Note that we are NOT including the halo nodes here. ---*/
	double *Coords_Local; jPoint = 0;
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check for halos and write only if requested ---*/
		if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
      
      /*--- Retrieve local coordinates at this node. ---*/
      Coords_Local = geometry->node[iPoint]->GetCoord();
      
      /*--- Load local coords into the temporary send buffer. ---*/
      Buffer_Send_X[jPoint] = Coords_Local[0];
      Buffer_Send_Y[jPoint] = Coords_Local[1];
      if (nDim == 3) Buffer_Send_Z[jPoint] = Coords_Local[2];
      
      /*--- Store the global index for this local node. ---*/
      Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Increment jPoint as the counter. We need this because iPoint
       may include halo nodes that we skip over during this loop. ---*/
      jPoint++;
		}
	}

	/*--- Gather the coordinate data on the master node using MPI. ---*/

	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Gather(Buffer_Send_X, nBuffer_Scalar, MPI::DOUBLE,
			Buffer_Recv_X, nBuffer_Scalar, MPI::DOUBLE,
			MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Y, nBuffer_Scalar, MPI::DOUBLE,
			Buffer_Recv_Y, nBuffer_Scalar, MPI::DOUBLE,
			MASTER_NODE);
	if (nDim == 3) {
		MPI::COMM_WORLD.Gather(Buffer_Send_Z, nBuffer_Scalar, MPI::DOUBLE,
				Buffer_Recv_Z, nBuffer_Scalar, MPI::DOUBLE,
				MASTER_NODE);
	}
	MPI::COMM_WORLD.Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
			Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
			MASTER_NODE);

	/*--- The master node unpacks and sorts this variable by global index ---*/

	if (rank == MASTER_NODE) {
		jPoint = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
			for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {

				/*--- Get global index, then loop over each variable and store ---*/
				iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
				Coords[0][iGlobal_Index] = Buffer_Recv_X[jPoint];
				Coords[1][iGlobal_Index] = Buffer_Recv_Y[jPoint];
				if (nDim == 3) Coords[2][iGlobal_Index] = Buffer_Recv_Z[jPoint];
				jPoint++;
			}
			/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
			jPoint = (iProcessor+1)*nBuffer_Scalar;
		}
	}

	/*--- Immediately release the temporary data buffers. ---*/

	delete [] Buffer_Send_X;
	delete [] Buffer_Send_Y;
	if (nDim == 3) delete [] Buffer_Send_Z;
	delete [] Buffer_Send_GlobalIndex;
	if (rank == MASTER_NODE) {
		delete [] Buffer_Recv_X;
		delete [] Buffer_Recv_Y;
		if (nDim == 3)  delete [] Buffer_Recv_Z;
		delete [] Buffer_Recv_GlobalIndex;
		delete [] Buffer_Recv_nPoin;
	}

#endif

}

void COutput::MergeVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	/*--- Local variables needed on all processors ---*/

	unsigned short NODES_PER_ELEMENT;

	unsigned long iPoint, iNode, jNode;
	unsigned long iElem = 0, jElem = 0;
	unsigned long nLocalElem = 0, nElem_Total = 0;

	int *Conn_Elem;

	/*--- Store the local number of this element type and the number of nodes
     per this element type. In serial, this will be the total number of this
     element type in the entire mesh. In parallel, it is the number on only
     the current partition. ---*/

	switch (Elem_Type) {
	case TRIANGLE:
		nLocalElem = geometry->GetnElemTria();
		NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
		break;
	case RECTANGLE:
		nLocalElem = geometry->GetnElemQuad();
		NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
		break;
	case TETRAHEDRON:
		nLocalElem = geometry->GetnElemTetr();
		NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
		break;
	case HEXAHEDRON:
		nLocalElem = geometry->GetnElemHexa();
		NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
		break;
	case WEDGE:
		nLocalElem = geometry->GetnElemWedg();
		NODES_PER_ELEMENT = N_POINTS_WEDGE;
		break;
	case PYRAMID:
		nLocalElem = geometry->GetnElemPyra();
		NODES_PER_ELEMENT = N_POINTS_PYRAMID;
		break;
	default:
		cout << "Error: Unrecognized element type \n";
		exit(0); break;
	}

	/*--- Merge the connectivity in serial or parallel. ---*/

#ifdef NO_MPI

	/*--- In serial, the single process has access to all connectivity,
     so simply load it into the data structure. ---*/

	/*--- Allocate a temporary array for the connectivity ---*/
	Conn_Elem = new int[nLocalElem*NODES_PER_ELEMENT];

	/*--- Load all elements of the current type into the buffer
     to be sent to the master node. ---*/
	jNode = 0; jElem = 0; nElem_Total = 0; bool isHalo;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		if(geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {

			/*--- Check if this is a halo node. ---*/
			isHalo = false;
			for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
				iPoint = geometry->elem[iElem]->GetNode(iNode);
				if (!geometry->node[iPoint]->GetDomain())
					isHalo = true;
			}

			/*--- Loop over all nodes in this element and load the
       connectivity into the temporary array. Do not merge any
       halo cells (periodic BC). Note that we are adding one to
       the index value because CGNS/Tecplot use 1-based indexing. ---*/
      
			if (!isHalo) {
				nElem_Total++;
				for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
					Conn_Elem[jNode] = (int)geometry->elem[iElem]->GetNode(iNode) + 1;

					/*--- Increment jNode as the counter. ---*/
					jNode++;
				}
			}
		}
	}

#else

	/*--- MPI preprocessing ---*/

	int iProcessor, jProcessor;
	int nProcessor = MPI::COMM_WORLD.Get_size();

	/*--- Local variables needed for merging the geometry with MPI. ---*/

  unsigned long iVertex, iMarker;
  
  int SendRecv, RecvFrom;
  
	unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
	unsigned long nBuffer_Scalar = 0;
	unsigned long kNode = 0, kElem = 0, pElem = 0;
	unsigned long MaxLocalElem = 0;

  bool Wrt_Halo = config->GetWrt_Halo();
	bool *Write_Elem;

	/*--- Find the max number of this element type among all
     partitions and set up buffers. ---*/

	Buffer_Send_nElem[0] = nLocalElem;
	if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[nProcessor];

	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalElem, &MaxLocalElem,
			1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nElem, 1, MPI::UNSIGNED_LONG,
			Buffer_Recv_nElem, 1, MPI::UNSIGNED_LONG, MASTER_NODE);

	nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;

	/*--- Send and Recv buffers ---*/

	int *Buffer_Send_Elem = new int[nBuffer_Scalar];
	int *Buffer_Recv_Elem = NULL;

	int *Buffer_Send_Halo = new int[MaxLocalElem];
	int *Buffer_Recv_Halo = NULL;

  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = false;
  
	/*--- Prepare the receive buffers on the master node only. ---*/

	if (rank == MASTER_NODE) {
		Buffer_Recv_Elem = new int[nProcessor*nBuffer_Scalar];
		Buffer_Recv_Halo = new int[nProcessor*MaxLocalElem];
		Conn_Elem = new int[nProcessor*MaxLocalElem*NODES_PER_ELEMENT];
	}

  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and 
   choose to keep only the halo cells from the lower rank processor. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			RecvFrom = abs(SendRecv)-1;
      if (SendRecv < 0 && RecvFrom < rank) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Local_Halo[iPoint] = true;
        }
      }
    }
  }
  
	/*--- Loop over all elements in this partition and load the
     elements of the current type into the buffer to be sent to
     the master node. ---*/
  
	jNode = 0; jElem = 0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		if(geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {

			/*--- Loop over all nodes in this element and load the
             connectivity into the send buffer. ---*/
      
			Buffer_Send_Halo[jElem] = false;
			for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {

				/*--- Store the global index values directly. ---*/
        
				iPoint = geometry->elem[iElem]->GetNode(iNode);
				Buffer_Send_Elem[jNode] = (int)geometry->node[iPoint]->GetGlobalIndex();

				/*--- Check if this is a halo node. If so, flag this element
         as a halo cell. We will use this later to sort and remove
         any duplicates from the connectivity list. ---*/
        
        if (Local_Halo[iPoint])
					Buffer_Send_Halo[jElem] = true;

				/*--- Increment jNode as the counter. We need this because iElem
         may include other elements that we skip over during this loop. ---*/
        
				jNode++;
			}
			jElem++;
		}
	}

	/*--- Gather the element connectivity information. ---*/

	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI::INT,
			Buffer_Recv_Elem, nBuffer_Scalar, MPI::INT,
			MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Halo, MaxLocalElem, MPI::INT,
			Buffer_Recv_Halo, MaxLocalElem, MPI::INT,
			MASTER_NODE);

	/*--- The master node unpacks and sorts the connectivity. ---*/

	if (rank == MASTER_NODE) {

		/*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/
    
		Write_Elem = new bool[nProcessor*MaxLocalElem];
		for (iElem = 0; iElem < nProcessor*MaxLocalElem; iElem++) {
			Write_Elem[iElem] = true;
		}

    /*--- Remove the rind layer from the solution only if requested ---*/
    
    if (!Wrt_Halo) {
      
      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/
      
      kElem = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
          
          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;
            
          }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }

		/*--- Store the unique connectivity list for this element type. ---*/

		jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
			for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

				/*--- Only write the elements that were flagged for it. ---*/
				if (Write_Elem[jElem+iElem]) {

					/*--- Increment total count for this element type ---*/
					nElem_Total++;

					/*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/
          
					for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
						Conn_Elem[kNode] = Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
						kNode++;
					}
				}
			}
			/*--- Adjust jNode to index of next proc's data in the buffers. ---*/
			jElem = (iProcessor+1)*MaxLocalElem;
			jNode = (iProcessor+1)*nBuffer_Scalar;
		}
	}

	/*--- Immediately release the temporary buffers. ---*/
	delete [] Buffer_Send_Elem;
	delete [] Buffer_Send_Halo;
  delete [] Local_Halo;
	if (rank == MASTER_NODE) {
		delete [] Buffer_Recv_nElem;
		delete [] Buffer_Recv_Elem;
		delete [] Buffer_Recv_Halo;
		delete [] Write_Elem;
	}

#endif

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

	if (rank == MASTER_NODE) {
		switch (Elem_Type) {
		case TRIANGLE:
			nGlobal_Tria = nElem_Total;
			if (nGlobal_Tria > 0) Conn_Tria = Conn_Elem;
			break;
		case RECTANGLE:
			nGlobal_Quad = nElem_Total;
			if (nGlobal_Quad > 0) Conn_Quad = Conn_Elem;
			break;
		case TETRAHEDRON:
			nGlobal_Tetr = nElem_Total;
			if (nGlobal_Tetr > 0) Conn_Tetr = Conn_Elem;
			break;
		case HEXAHEDRON:
			nGlobal_Hexa = nElem_Total;
			if (nGlobal_Hexa > 0) Conn_Hexa = Conn_Elem;
			break;
		case WEDGE:
			nGlobal_Wedg = nElem_Total;
			if (nGlobal_Wedg > 0) Conn_Wedg = Conn_Elem;
			break;
		case PYRAMID:
			nGlobal_Pyra = nElem_Total;
			if (nGlobal_Pyra > 0) Conn_Pyra = Conn_Elem;
			break;
		default:
			cout << "Error: Unrecognized element type \n";
			exit(0); break;
		}
	}

}

void COutput::MergeSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
	/*--- Local variables needed on all processors ---*/
  
	unsigned short NODES_PER_ELEMENT;
  
  unsigned short iMarker;
	unsigned long iPoint, iNode, jNode;
	unsigned long iElem = 0, jElem = 0;
	unsigned long nLocalElem = 0, nElem_Total = 0;
  
	int *Conn_Elem;
  
	/*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  nLocalElem = 0;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Plotting(iMarker) == YES) {
			for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          nLocalElem++;
        }
      }
    }
  }
  
  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case RECTANGLE:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      cout << "Error: Unrecognized element type \n";
      exit(0); break;
  }
  
	/*--- Merge the connectivity in serial or parallel. ---*/
  
#ifdef NO_MPI
  
	/*--- In serial, the single process has access to all connectivity,
   so simply load it into the data structure. ---*/
  
	/*--- Allocate a temporary array for the connectivity ---*/
	Conn_Elem = new int[nLocalElem*NODES_PER_ELEMENT];
  
	/*--- Load all elements of the current type into the buffer
   to be sent to the master node. ---*/
	jNode = 0; jElem = 0; nElem_Total = 0; bool isHalo;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Plotting(iMarker) == YES)
			for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        
        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          
          /*--- Check if this is a halo node. ---*/
          isHalo = false;
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
            if (!geometry->node[iPoint]->GetDomain())
              isHalo = true;
          }
          
          /*--- Loop over all nodes in this element and load the
           connectivity into the temporary array. Do not merge any
           halo cells (periodic BC). Note that we are adding one to
           the index value because CGNS/Tecplot use 1-based indexing. ---*/
          if (!isHalo) {
            nElem_Total++;
            for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
              Conn_Elem[jNode] = (int)geometry->bound[iMarker][iElem]->GetNode(iNode) + 1;
              
              /*--- Increment jNode as the counter. ---*/
              jNode++;
            }
          }
        }
      }
  
#else
  
	/*--- MPI preprocessing ---*/
  
	int iProcessor, jProcessor;
	int nProcessor = MPI::COMM_WORLD.Get_size();
  
	/*--- Local variables needed for merging the geometry with MPI. ---*/
  
  unsigned long iVertex;
  
  int SendRecv, RecvFrom;
  
	unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
	unsigned long nBuffer_Scalar = 0;
	unsigned long kNode = 0, kElem = 0, pElem = 0;
	unsigned long MaxLocalElem = 0;
  
  bool Wrt_Halo = config->GetWrt_Halo();
	bool *Write_Elem;
  
	/*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/
  
	Buffer_Send_nElem[0] = nLocalElem;
	if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[nProcessor];

	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalElem, &MaxLocalElem,
                            1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nElem, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nElem, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
  
	nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;
  
	/*--- Send and Recv buffers ---*/
  
	int *Buffer_Send_Elem = new int[nBuffer_Scalar];
	int *Buffer_Recv_Elem = NULL;
  
	int *Buffer_Send_Halo = new int[MaxLocalElem];
	int *Buffer_Recv_Halo = NULL;

  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = false;
  
	/*--- Prepare the receive buffers on the master node only. ---*/
  
	if (rank == MASTER_NODE) {
		Buffer_Recv_Elem = new int[nProcessor*nBuffer_Scalar];
		Buffer_Recv_Halo = new int[nProcessor*MaxLocalElem];
		Conn_Elem = new int[nProcessor*MaxLocalElem*NODES_PER_ELEMENT];
	}
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the lower rank processor. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			RecvFrom = abs(SendRecv)-1;
      if (SendRecv < 0 && RecvFrom < rank) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Local_Halo[iPoint] = true;
        }
      }
    }
  }
  
	/*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/
	jNode = 0; jElem = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Plotting(iMarker) == YES)
			for(iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        
        if(geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          
          /*--- Loop over all nodes in this element and load the
           connectivity into the send buffer. ---*/
          
          Buffer_Send_Halo[jElem] = false;
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            
            /*--- Store the global index values directly. ---*/
            
            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
            Buffer_Send_Elem[jNode] = (int)geometry->node[iPoint]->GetGlobalIndex();
            
            /*--- Check if this is a halo node. If so, flag this element
             as a halo cell. We will use this later to sort and remove
             any duplicates from the connectivity list. ---*/
            
            if (Local_Halo[iPoint])
              Buffer_Send_Halo[jElem] = true;
            
            /*--- Increment jNode as the counter. We need this because iElem
             may include other elements that we skip over during this loop. ---*/
            
            jNode++;
          }
          jElem++;
        }
      }

	/*--- Gather the element connectivity information. ---*/
  
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI::INT,
                         Buffer_Recv_Elem, nBuffer_Scalar, MPI::INT,
                         MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Halo, MaxLocalElem, MPI::INT,
                         Buffer_Recv_Halo, MaxLocalElem, MPI::INT,
                         MASTER_NODE);
  
	/*--- The master node unpacks and sorts the connectivity. ---*/
  
	if (rank == MASTER_NODE) {
    
    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/
    
    Write_Elem = new bool[nProcessor*MaxLocalElem];
    for (iElem = 0; iElem < nProcessor*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }
    
    /*--- Remove the rind layer from the solution only if requested ---*/
    
    if (!Wrt_Halo) {
      
      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/
      
      kElem = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
          
          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;
          
        }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }
    
		/*--- Store the unique connectivity list for this element type. ---*/
    
		jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
			for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
				/*--- Only write the elements that were flagged for it. ---*/
				if (Write_Elem[jElem+iElem]) {
          
					/*--- Increment total count for this element type ---*/
					nElem_Total++;
          
					/*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/
          
					for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
						Conn_Elem[kNode] = Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
						kNode++;
					}
				}
			}
			/*--- Adjust jNode to index of next proc's data in the buffers. ---*/
			jElem = (iProcessor+1)*MaxLocalElem;
			jNode = (iProcessor+1)*nBuffer_Scalar;
		}
	}
  
	/*--- Immediately release the temporary buffers. ---*/
	delete [] Buffer_Send_Elem;
	delete [] Buffer_Send_Halo;
  delete [] Local_Halo;
	if (rank == MASTER_NODE) {
		delete [] Buffer_Recv_nElem;
		delete [] Buffer_Recv_Elem;
		delete [] Buffer_Recv_Halo;
		delete [] Write_Elem;
	}
  
#endif
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
	if (rank == MASTER_NODE) {
		switch (Elem_Type) {
      case LINE:
        nGlobal_Line = nElem_Total;
        if (nGlobal_Line > 0) Conn_Line = Conn_Elem;
        break;
      case TRIANGLE:
        nGlobal_BoundTria = nElem_Total;
        if (nGlobal_BoundTria > 0) Conn_BoundTria = Conn_Elem;
        break;
      case RECTANGLE:
        nGlobal_BoundQuad = nElem_Total;
        if (nGlobal_BoundQuad > 0) Conn_BoundQuad = Conn_Elem;
        break;
      default:
        cout << "Error: Unrecognized element type \n";
        exit(0); break;
		}
	}
  
}

void COutput::MergeSolution(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
	/*--- Local variables needed on all processors ---*/
	unsigned short Kind_Solver  = config->GetKind_Solver();
	unsigned short iVar, jVar, iSpecies, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
	unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0, iVar_Eddy = 0;
	unsigned short iVar_GridVel = 0, iVar_PressMach = 0, iVar_Density = 0, iVar_TempLam = 0,
  iVar_Tempv = 0,iVar_MagF = 0, iVar_EF =0, iVar_Temp = 0, iVar_Lam =0, iVar_Mach = 0, iVar_Press = 0,
  iVar_ViscCoeffs = 0, iVar_Sens = 0;
  
	unsigned long iPoint = 0, jPoint = 0, iVertex = 0, iMarker = 0;
  
  double *Aux_Press, *Aux_Frict, *Aux_Heat, *Aux_yPlus, *Aux_Sens;
  
	bool grid_movement = ((config->GetUnsteady_Simulation() &&
                         config->GetWrt_Unsteady() &&
                         config->GetGrid_Movement()) ||
                        ((config->GetUnsteady_Simulation() == TIME_SPECTRAL) &&
                         config->GetGrid_Movement()));
	bool incompressible = config->GetIncompressible();
  bool transition = (config->GetKind_Trans_Model()==LM);
  
	if (Kind_Solver == AEROACOUSTIC_EULER) {
		if (val_iZone == ZONE_0) Kind_Solver = EULER;
		if (val_iZone == ZONE_1) Kind_Solver = WAVE_EQUATION;
	}
	if (Kind_Solver == PLASMA_EULER) {
		if (val_iZone == ZONE_0) Kind_Solver = PLASMA_EULER;
		if (val_iZone == ZONE_1) Kind_Solver = ELECTRIC_POTENTIAL;
	}
	if (Kind_Solver == PLASMA_NAVIER_STOKES) {
		if (val_iZone == ZONE_0) Kind_Solver = PLASMA_NAVIER_STOKES;
		if (val_iZone == ZONE_1) Kind_Solver = ELECTRIC_POTENTIAL;
	}
  
	/*--- Prepare send buffers for the conservative variables. Need to
   find the total number of conservative variables and also the
   index for their particular solution container. ---*/
	switch (Kind_Solver) {
    case EULER : case NAVIER_STOKES:
      FirstIndex = FLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
      FirstIndex = PLASMA_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case RANS :
      FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL;
      if (transition) ThirdIndex=TRANS_SOL;
      else ThirdIndex = NONE;
      break;
    case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES:
      FirstIndex = FLOW_SOL; SecondIndex = LEVELSET_SOL; ThirdIndex = NONE;
      break;
    case FREE_SURFACE_RANS:
      FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; ThirdIndex = LEVELSET_SOL;
      break;
    case ELECTRIC_POTENTIAL:
      FirstIndex = ELEC_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case WAVE_EQUATION:
      FirstIndex = WAVE_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case LINEAR_ELASTICITY:
      FirstIndex = FEA_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_EULER : case ADJ_NAVIER_STOKES :
      FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_PLASMA_EULER : case ADJ_PLASMA_NAVIER_STOKES :
      FirstIndex = ADJPLASMA_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_FREE_SURFACE_EULER : case ADJ_FREE_SURFACE_NAVIER_STOKES:  case ADJ_FREE_SURFACE_RANS:
      FirstIndex = ADJFLOW_SOL; SecondIndex = ADJLEVELSET_SOL; ThirdIndex = NONE;
      break;
    case ADJ_RANS :
      FirstIndex = ADJFLOW_SOL;
      if ((config->GetFrozen_Visc()) && (config->GetKind_Adjoint() != HYBRID))
        SecondIndex = NONE;
      else
        SecondIndex = ADJTURB_SOL;
      ThirdIndex = NONE;
      break;
    case LIN_EULER : case LIN_NAVIER_STOKES : ThirdIndex = NONE;
      FirstIndex = LINFLOW_SOL; SecondIndex = NONE;
      break;
    default: SecondIndex = NONE; ThirdIndex = NONE;
      break;
	}
	nVar_First = solver[FirstIndex]->GetnVar();
	if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
	if (ThirdIndex != NONE) nVar_Third = solver[ThirdIndex]->GetnVar();
	nVar_Consv = nVar_First + nVar_Second + nVar_Third;
  if (config->GetWrt_Residuals()) {
    nVar_Total = 2*nVar_Consv;
  } else {
    nVar_Total = nVar_Consv;
  }
  
	/*--- Add the grid velocity to the restart file for the unsteady adjoint ---*/
	if (grid_movement) {
		iVar_GridVel = nVar_Total;
		if (geometry->GetnDim() == 2) nVar_Total += 2;
		else if (geometry->GetnDim() == 3) nVar_Total += 3;
	}
  
  if ((Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) ||
			(Kind_Solver == FREE_SURFACE_RANS)) {
		/*--- Density ---*/
		iVar_Density = nVar_Total;
		nVar_Total += 1;
	}
  
	if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
			(Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) ||
			(Kind_Solver == FREE_SURFACE_RANS)) {
		/*--- Pressure, Cp, and Mach ---*/
		iVar_PressMach = nVar_Total;
		nVar_Total += 3;
	}
  
	if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
			(Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
		/*--- Temperature, Laminar Viscosity ---*/
		iVar_TempLam = nVar_Total;
		nVar_Total += 2;
    /*--- Skin Friction, Heat Flux, & yPlus ---*/
    iVar_ViscCoeffs = nVar_Total;
		nVar_Total += 3;
	}
  
	if ((Kind_Solver == RANS) || (Kind_Solver == FREE_SURFACE_RANS)) {
		/*--- Eddy Viscosity ---*/
		iVar_Eddy = nVar_Total;
		nVar_Total += 1;
	}
  
	if (Kind_Solver == ELECTRIC_POTENTIAL) {
		iVar_EF = geometry->GetnDim();
		nVar_Total += geometry->GetnDim();
	}
  
	if ((Kind_Solver == PLASMA_EULER) || (Kind_Solver == PLASMA_NAVIER_STOKES)) {
		iVar_Press  = nVar_Total;
		nVar_Total += config->GetnSpecies();
		iVar_Temp   = nVar_Total;
		nVar_Total += config->GetnSpecies();
    iVar_Tempv  = nVar_Total;
    nVar_Total += config->GetnDiatomics();
    iVar_Mach   = nVar_Total;
    nVar_Total += config->GetnSpecies();
	}
  
	if (Kind_Solver == PLASMA_NAVIER_STOKES) {
		iVar_Lam = nVar_Total;
		nVar_Total  += config->GetnSpecies();
		if((config->GetMagnetic_Force() == YES) && (geometry->GetnDim() ==3)) {
			iVar_MagF   = nVar_Total;
			nVar_Total  += 3;
		}
	}
  
  if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) ||
			(Kind_Solver == ADJ_FREE_SURFACE_EULER) || (Kind_Solver == ADJ_FREE_SURFACE_NAVIER_STOKES) ||
			(Kind_Solver == ADJ_FREE_SURFACE_RANS) || (Kind_Solver == ADJ_PLASMA_EULER) || (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES)) {
    /*--- Surface sensitivity coefficient ---*/
    iVar_Sens   = nVar_Total;
    nVar_Total += 1;
  }
  
	/*--- Merge the solution either in serial or parallel. ---*/
  
#ifdef NO_MPI
  
  /*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/
	nGlobal_Poin = geometry->GetnPointDomain();
	Data = new double*[nVar_Total];
	for (iVar = 0; iVar < nVar_Total; iVar++) {
		Data[iVar] = new double[nGlobal_Poin];
	}
  
	/*--- In case there is grid movement ---*/
	double *Grid_Vel;
  
  /*--- First, loop through the mesh in order to find and store the
   value of the coefficient of pressure at any surface nodes. They
   will be placed in an auxiliary vector and then communicated like
   all other volumetric variables. ---*/
  
  Aux_Press = new double [geometry->GetnPointDomain()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Press[iPoint] = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Aux_Press[iPoint] = solver[FLOW_SOL]->GetCPressure(iMarker,iVertex);
      }
  
	if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
    
    Aux_Frict = new double [geometry->GetnPointDomain()];
		Aux_Heat  = new double [geometry->GetnPointDomain()];
		Aux_yPlus = new double [geometry->GetnPointDomain()];
    
    for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
      Aux_Frict[iPoint] = 0.0;
      Aux_Heat[iPoint]  = 0.0;
      Aux_yPlus[iPoint] = 0.0;
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Frict[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker,iVertex);
          Aux_Heat[iPoint]  = solver[FLOW_SOL]->GetHeatTransferCoeff(iMarker,iVertex);
          Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker,iVertex);
        }
      }
  }
  
  if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) ||
			(Kind_Solver == ADJ_FREE_SURFACE_EULER) || (Kind_Solver == ADJ_FREE_SURFACE_NAVIER_STOKES) ||
			(Kind_Solver == ADJ_FREE_SURFACE_RANS) || (Kind_Solver == ADJ_PLASMA_EULER) || (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES)) {
    
    Aux_Sens = new double [geometry->GetnPointDomain()];
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Sens[iPoint] = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Sens[iPoint] = solver[ADJFLOW_SOL]->GetCSensitivity(iMarker,iVertex);
        }
      }
    
  }
  
	/*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (ignore periodic/sliding halo nodes). ---*/
	jPoint = 0;
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check for halo nodes & only write if requested ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Solution (first, second and third system of equations) ---*/
      jVar = 0;
      for (iVar = 0; iVar < nVar_First; iVar++) {
        Data[jVar][jPoint] = solver[FirstIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      for (iVar = 0; iVar < nVar_Second; iVar++) {
        Data[jVar][jPoint] = solver[SecondIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      for (iVar = 0; iVar < nVar_Third; iVar++) {
        Data[jVar][jPoint] = solver[ThirdIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      /*--- Residual (first, second and third system of equations) ---*/
      if (config->GetWrt_Residuals()) {
        for (iVar = 0; iVar < nVar_First; iVar++) {
          Data[jVar][jPoint] = solver[FirstIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
        
        for (iVar = 0; iVar < nVar_Second; iVar++) {
          Data[jVar][jPoint] = solver[SecondIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
        
        for (iVar = 0; iVar < nVar_Third; iVar++) {
          Data[jVar][jPoint] = solver[ThirdIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
      }
      
      /*--- For unsteady problems with grid movement, write the mesh velocities
       also, in case we need them for the unsteady adjoint. ---*/
      if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady() && grid_movement) {
        Grid_Vel = geometry->node[iPoint]->GetGridVel();
        for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
          Data[jVar][jPoint] = Grid_Vel[iDim];
          jVar++;
        }
      }
      
      /*--- Any extra data for output files ---*/
      switch (Kind_Solver) {
        case EULER: case FREE_SURFACE_EULER:
          /*--- Load buffers with the pressure, Cp, and mach variables. ---*/
          if (incompressible) {
            if (Kind_Solver == FREE_SURFACE_EULER) {
              Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc(); jVar++;
            }
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(incompressible); jVar++;
            Data[jVar][jPoint] = Aux_Press[iPoint]; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
          } else {
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(incompressible); jVar++;
            Data[jVar][jPoint] = Aux_Press[iPoint]; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
            solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
          }
          break;
          /*--- Write pressure, Cp, mach, temperature, laminar viscosity, skin friction, heat transfer, yplus ---*/
        case NAVIER_STOKES: case FREE_SURFACE_NAVIER_STOKES:
          if (incompressible) {
            if (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) {
              Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc(); jVar++;
            }
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(incompressible); jVar++;
            Data[jVar][jPoint] = Aux_Press[iPoint]; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
            Data[jVar][jPoint] = 0.0; jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
          } else {
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(incompressible); jVar++;
            Data[jVar][jPoint] = Aux_Press[iPoint]; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
            solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
          }
          break;
          /*--- Write pressure, Cp, mach, temperature, laminar viscosity, skin friction, heat transfer, yplus, eddy viscosity ---*/
        case RANS: case FREE_SURFACE_RANS:
          if (incompressible) {
            if (Kind_Solver == FREE_SURFACE_RANS) {
              Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc(); jVar++;
            }
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(incompressible); jVar++;
            Data[jVar][jPoint] = Aux_Press[iPoint]; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
            Data[jVar][jPoint] = 0.0; jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
          } else {
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(incompressible); jVar++;
            Data[jVar][jPoint] = Aux_Press[iPoint]; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
            solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
          }
          Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity(); jVar++;
          break;
          /*--- Write electric field. ---*/
        case ELECTRIC_POTENTIAL:
          for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
            Data[jVar][jPoint] = -1.0*solver[ELEC_SOL]->node[iPoint]->GetGradient(0,iDim);
            jVar++;
          }
          break;
          
        case PLASMA_EULER:
          /*--- Write partial pressures ---*/
          for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
            Data[jVar][jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies);
            jVar++;
          }
          /*--- Write translational-rotational temperature ---*/
          for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
            Data[jVar][jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetTemperature_tr(iSpecies);
            jVar++;
          }
          /*--- Write vibrational temperature ---*/
          for (iSpecies = 0; iSpecies < config->GetnDiatomics(); iSpecies++) {
            Data[jVar][jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetTemperature_vib(iSpecies);
            jVar++;
          }
          /*--- Write Mach number ---*/
          for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
            Data[jVar][jPoint] =  sqrt(solver[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies))
            / solver[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies);
            jVar++;
          }
          break;
          
        case PLASMA_NAVIER_STOKES:
          /*--- Write partial pressures ---*/
          for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
            Data[jVar][jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies);
            jVar++;
          }
          /*--- Write translational-rotational temperature ---*/
          for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
            Data[jVar][jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetTemperature_tr(iSpecies);
            jVar++;
          }
          /*--- Write vibrational temperature ---*/
          for (iSpecies = 0; iSpecies < config->GetnDiatomics(); iSpecies++) {
            Data[jVar][jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetTemperature_vib(iSpecies);
            jVar++;
          }
          /*--- Write Mach number ---*/
          for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
            Data[jVar][jPoint] =  sqrt(solver[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies))
            / solver[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies);
            jVar++;
          }
          /*--- Write laminar viscosity ---*/
          for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
            Data[jVar][jPoint] =  solver[PLASMA_SOL]->node[iPoint]->GetLaminarViscosity(iSpecies);
            jVar++;
          }
          /*--- Write magnetic force ---*/
          if((config->GetMagnetic_Force() == YES) && (geometry->GetnDim() == 3)) {
            for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
              Data[jVar][jPoint] =  solver[PLASMA_SOL]->node[iPoint]->GetMagneticField()[iDim];
              jVar++;
            }
          }
          break;
          
        case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
        case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
        case ADJ_PLASMA_EULER: case ADJ_PLASMA_NAVIER_STOKES:
          
          Data[jVar][jPoint] = Aux_Sens[iPoint]; jVar++;
          break;
          
      }
    }
    
    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/
    jPoint++;
    
	}
  
#else
  
	/*--- MPI preprocessing ---*/
  
	int rank = MPI::COMM_WORLD.Get_rank();
	int nProcessor = MPI::COMM_WORLD.Get_size();
	int iProcessor;
  
	/*--- Local variables needed for merging with MPI ---*/
	unsigned short CurrentIndex;
  
	unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
	unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
	unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  
	/*--- Each processor sends its local number of nodes to the master. ---*/

  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else
    nLocalPoint = geometry->GetnPointDomain();
	Buffer_Send_nPoint[0] = nLocalPoint;
	if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nPoint, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nPoint, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
	nBuffer_Scalar = MaxLocalPoint;
  
	/*--- Send and Recv buffers. ---*/
  
	double *Buffer_Send_Var = new double[MaxLocalPoint];
	double *Buffer_Recv_Var = NULL;
  
	double *Buffer_Send_Res = new double[MaxLocalPoint];
	double *Buffer_Recv_Res = NULL;
  
	double *Buffer_Send_Vol = new double[MaxLocalPoint];
	double *Buffer_Recv_Vol = NULL;
  
	unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
	unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Auxiliary vectors for surface coefficients ---*/
  
  Aux_Press = new double[geometry->GetnPoint()];
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
    Aux_Frict = new double[geometry->GetnPoint()];
    Aux_Heat  = new double[geometry->GetnPoint()];
    Aux_yPlus = new double[geometry->GetnPoint()];
  }
  
  if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) ||
			(Kind_Solver == ADJ_FREE_SURFACE_EULER) || (Kind_Solver == ADJ_FREE_SURFACE_NAVIER_STOKES) ||
			(Kind_Solver == ADJ_FREE_SURFACE_RANS) || (Kind_Solver == ADJ_PLASMA_EULER) || (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES)) {
    Aux_Sens = new double[geometry->GetnPoint()];
    
  }
  
	/*--- Prepare the receive buffers in the master node only. ---*/
  
	if (rank == MASTER_NODE) {
    
		Buffer_Recv_Var = new double[nProcessor*MaxLocalPoint];
		Buffer_Recv_Res = new double[nProcessor*MaxLocalPoint];
		Buffer_Recv_Vol = new double[nProcessor*MaxLocalPoint];
		Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
		/*--- Sum total number of nodes to be written and allocate arrays ---*/
		nGlobal_Poin = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
			nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
		}
		Data = new double*[nVar_Total];
		for (iVar = 0; iVar < nVar_Total; iVar++) {
			Data[iVar] = new double[nGlobal_Poin];
		}
	}
  
	/*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
  
	for (iVar = 0; iVar < nVar_Consv; iVar++) {
    
		/*--- Logic for which solution class to draw from. ---*/
    
		jVar = iVar;
		CurrentIndex = FirstIndex;
		if ((SecondIndex != NONE) && (iVar > nVar_First-1)) {
			jVar = iVar - nVar_First;
			CurrentIndex = SecondIndex;
		}
    if ((SecondIndex != NONE) && (ThirdIndex != NONE) && (iVar > (nVar_First + nVar_Second-1))) {
			jVar = iVar - nVar_First - nVar_Second;
			CurrentIndex = ThirdIndex;
		}
    
		/*--- Loop over this partition to collect the current variable ---*/
    
		jPoint = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Get this variable into the temporary send buffer. ---*/
        Buffer_Send_Var[jPoint] = solver[CurrentIndex]->node[iPoint]->GetSolution(jVar);
        if (config->GetWrt_Residuals()) {
          Buffer_Send_Res[jPoint] = solver[CurrentIndex]->LinSysRes.GetBlock(iPoint, jVar);
        }
        
        /*--- Only send/recv the volumes & global indices during the first loop ---*/
        if (iVar == 0) {
          Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
        }
        jPoint++;
      }
		}
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    if (config->GetWrt_Residuals()) {
      MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
    }
		if (iVar == 0) {
			MPI::COMM_WORLD.Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                             Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                             MASTER_NODE);
		}
    
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
			jPoint = 0;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          if (config->GetWrt_Residuals()) {
            Data[iVar+nVar_Consv][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          }
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
	}
  
	/*--- Additional communication routine for the grid velocity. Note that
   we are reusing the same temporary buffers from above for efficiency.
   Also, in the future more routines like this could be used to write
   an arbitrary number of additional variables to the file. ---*/
  
	if (grid_movement) {
    
		/*--- Loop over this partition to collect the current variable ---*/
		jPoint = 0; double *Grid_Vel;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the three grid velocity components. ---*/
        Grid_Vel = geometry->node[iPoint]->GetGridVel();
        Buffer_Send_Var[jPoint] = Grid_Vel[0];
        Buffer_Send_Res[jPoint] = Grid_Vel[1];
        if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Grid_Vel[2];
        jPoint++;
      }
    }
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		if (geometry->GetnDim() == 3) {
			MPI::COMM_WORLD.Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Vol, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
		}
    
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
			jPoint = 0; iVar = iVar_GridVel;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
					Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
					if (geometry->GetnDim() == 3)
						Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
	}
  
  /*--- Communicate the Density in Free-surface problems ---*/
	if ((Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
    
		/*--- Loop over this partition to collect the current variable ---*/
		jPoint = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the pressure and mach variables. ---*/
        if (incompressible) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc();
        }
        else {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc();
        }
        jPoint++;
      }
		}
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
			jPoint = 0; iVar = iVar_Density;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
    
	}
  
	/*--- Communicate Pressure, Cp, and Mach ---*/
  
	if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
    
    /*--- First, loop through the mesh in order to find and store the
     value of the coefficient of pressure at any surface nodes. They 
     will be placed in an auxiliary vector and then communicated like
     all other volumetric variables. ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Press[iPoint] = 0.0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					Aux_Press[iPoint] = solver[FLOW_SOL]->GetCPressure(iMarker,iVertex);
				}
    
		/*--- Loop over this partition to collect the current variable ---*/
		jPoint = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the pressure, Cp, and mach variables. ---*/
        if (incompressible) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(incompressible);
          Buffer_Send_Res[jPoint] = Aux_Press[iPoint];
          Buffer_Send_Vol[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/
          sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref()));
        }
        else {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(incompressible);
          Buffer_Send_Res[jPoint] = Aux_Press[iPoint];
          Buffer_Send_Vol[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
          solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
        }
        jPoint++;
      }
    }
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    MPI::COMM_WORLD.Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Vol, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
			jPoint = 0; iVar = iVar_PressMach;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
					Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
	}
  
	/*--- Communicate Temperature & Laminar Viscosity ---*/
	if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
    
		/*--- Loop over this partition to collect the current variable ---*/
		jPoint = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
        if (incompressible) {
          Buffer_Send_Var[jPoint] = 0.0;
          Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc();
        }
        else {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature();
          Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
        }
        jPoint++;
      }
		}
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
			jPoint = 0; iVar = iVar_TempLam;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
					Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
	}
  
  /*--- Communicate skin friction, heat transfer, y+ ---*/
	if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
    
    /*--- First, loop through the mesh in order to find and store the
     value of the viscous coefficients at any surface nodes. They
     will be placed in an auxiliary vector and then communicated like
     all other volumetric variables. ---*/
    
    for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Aux_Frict[iPoint] = 0.0;
      Aux_Heat[iPoint]  = 0.0;
      Aux_yPlus[iPoint] = 0.0;
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Frict[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker,iVertex);
          Aux_Heat[iPoint]  = solver[FLOW_SOL]->GetHeatTransferCoeff(iMarker,iVertex);
          Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker,iVertex);
        }
      }
    
		/*--- Loop over this partition to collect the current variable ---*/
		jPoint = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
        if (incompressible) {
          Buffer_Send_Var[jPoint] = Aux_Frict[iPoint];
          Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
          Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
        }
        else {
          Buffer_Send_Var[jPoint] = Aux_Frict[iPoint];
          Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
          Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
        }
        jPoint++;
      }
    }
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    MPI::COMM_WORLD.Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Vol, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
			jPoint = 0; iVar = iVar_ViscCoeffs;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
					Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
	}
  
	/*--- Communicate the Eddy Viscosity ---*/
	if ((Kind_Solver == RANS) || (Kind_Solver == FREE_SURFACE_RANS)) {
    
		/*--- Loop over this partition to collect the current variable ---*/
		jPoint = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the pressure and mach variables. ---*/
        if (incompressible) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity();
        }
        else {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity();
        }
        jPoint++;
      }
		}
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
			jPoint = 0; iVar = iVar_Eddy;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
    
	}
  
	/*--- Communicate additional variables for the plasma solvers ---*/
	if ((Kind_Solver == PLASMA_EULER) || (Kind_Solver == PLASMA_NAVIER_STOKES)) {
    
    /*--- Pressure ---*/
		for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies ++) {
			/*--- Loop over this partition to collect the current variable ---*/
			jPoint = 0;
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
          
          /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
          Buffer_Send_Var[jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies);
          jPoint++;
        }
        
			}
      
			/*--- Gather the data on the master node. ---*/
			MPI::COMM_WORLD.Barrier();
			MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
      
			/*--- The master node unpacks and sorts this variable by global index ---*/
			if (rank == MASTER_NODE) {
				jPoint = 0;
        iVar = iVar_Press + iSpecies;
				for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
					for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
						/*--- Get global index, then loop over each variable and store ---*/
						iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
						Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
						jPoint++;
					}
					/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
					jPoint = (iProcessor+1)*nBuffer_Scalar;
				}
			}
		}
    
    /*--- Translational-rotational temperature ---*/
		for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
			/*--- Loop over this partition to collect the current variable ---*/
			jPoint = 0;
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        
        /*--- Check for halos & write only if requested ---*/
        
        if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
          
          /*--- Load buffer with the temperature. ---*/
          Buffer_Send_Var[jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetTemperature_tr(iSpecies);
          jPoint++;
        }
			}
      
			/*--- Gather the data on the master node. ---*/
			MPI::COMM_WORLD.Barrier();
			MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
      
			/*--- The master node unpacks and sorts this variable by global index ---*/
			if (rank == MASTER_NODE) {
				jPoint = 0;
        iVar = iVar_Temp + iSpecies;
				for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
					for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
						/*--- Get global index, then loop over each variable and store ---*/
						iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
						Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
						jPoint++;
					}
					/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
					jPoint = (iProcessor+1)*nBuffer_Scalar;
				}
			}
		}
    
    /*--- Vibrational temperature ---*/
		for (iSpecies = 0; iSpecies < config->GetnDiatomics(); iSpecies++) {
			/*--- Loop over this partition to collect the current variable ---*/
			jPoint = 0;
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        
        /*--- Check for halos & write only if requested ---*/
        
        if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
          
          /*--- Load buffer with the vibrational temperature. ---*/
          Buffer_Send_Var[jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetTemperature_vib(iSpecies);
          jPoint++;
        }
			}
      
			/*--- Gather the data on the master node. ---*/
			MPI::COMM_WORLD.Barrier();
			MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
      
			/*--- The master node unpacks and sorts this variable by global index ---*/
			if (rank == MASTER_NODE) {
				jPoint = 0;
        iVar = iVar_Tempv + iSpecies;
				for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
					for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
						/*--- Get global index, then loop over each variable and store ---*/
						iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
						Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
						jPoint++;
					}
					/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
					jPoint = (iProcessor+1)*nBuffer_Scalar;
				}
			}
		}
    
    /*--- Mach number ---*/
		for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
			/*--- Loop over this partition to collect the current variable ---*/
			jPoint = 0;
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        
        /*--- Check for halos & write only if requested ---*/
        
        if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
          
          /*--- Load buffer with the vibrational temperature. ---*/
          Buffer_Send_Var[jPoint] = sqrt(solver[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies))
          / solver[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies);
          jPoint++;
        }
			}
      
			/*--- Gather the data on the master node. ---*/
			MPI::COMM_WORLD.Barrier();
			MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
      
			/*--- The master node unpacks and sorts this variable by global index ---*/
			if (rank == MASTER_NODE) {
				jPoint = 0;
        iVar = iVar_Mach + iSpecies;
				for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
					for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
						/*--- Get global index, then loop over each variable and store ---*/
						iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
						Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
						jPoint++;
					}
					/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
					jPoint = (iProcessor+1)*nBuffer_Scalar;
				}
			}
		}
	}
  
	if ( Kind_Solver == PLASMA_NAVIER_STOKES) {
    
		/*--- Laminar viscosity ---*/
		for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
			/*--- Loop over this partition to collect the current variable ---*/
			jPoint = 0;
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        
        /*--- Check for halos & write only if requested ---*/
        
        if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
          /*--- Load buffer with the laminar viscosity. ---*/
          Buffer_Send_Var[jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetLaminarViscosity(iSpecies);
          jPoint++;
          
        }
			}
      
			/*--- Gather the data on the master node. ---*/
			MPI::COMM_WORLD.Barrier();
			MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
      
			/*--- The master node unpacks and sorts this variable by global index ---*/
			if (rank == MASTER_NODE) {
				jPoint = 0;
        iVar = iVar_Lam + iSpecies;
				for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
					for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
						/*--- Get global index, then loop over each variable and store ---*/
						iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
						Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
						jPoint++;
					}
					/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
					jPoint = (iProcessor+1)*nBuffer_Scalar;
				}
			}
		}
	}
  
	/*--- Communicate the Eddy Viscosity ---*/
	if ( Kind_Solver == PLASMA_NAVIER_STOKES  && (config->GetMagnetic_Force() == YES) && (geometry->GetnDim() == 3)) {
    
		/*--- Loop over this partition to collect the current variable ---*/
		jPoint = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
        Buffer_Send_Var[jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetMagneticField()[0];
        Buffer_Send_Res[jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetMagneticField()[1];
        Buffer_Send_Vol[jPoint] = solver[PLASMA_SOL]->node[iPoint]->GetMagneticField()[2];
        jPoint++;
      }
		}
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		MPI::COMM_WORLD.Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Vol, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
			jPoint = 0; iVar = iVar_MagF;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
					Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
					Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
	}
  
  /*--- Communicate the surface sensitivity ---*/
  
  if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) ||
			(Kind_Solver == ADJ_FREE_SURFACE_EULER) || (Kind_Solver == ADJ_FREE_SURFACE_NAVIER_STOKES) ||
			(Kind_Solver == ADJ_FREE_SURFACE_RANS) || (Kind_Solver == ADJ_PLASMA_EULER) || (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES)) {
    
    /*--- First, loop through the mesh in order to find and store the
     value of the surface sensitivity at any surface nodes. They
     will be placed in an auxiliary vector and then communicated like
     all other volumetric variables. ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Sens[iPoint] = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Sens[iPoint] = solver[ADJFLOW_SOL]->GetCSensitivity(iMarker,iVertex);
        }
      }
    
		/*--- Loop over this partition to collect the current variable ---*/
		jPoint = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {
        
        /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
        Buffer_Send_Var[jPoint] = Aux_Sens[iPoint];
        
        jPoint++;
      }
    }
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
			jPoint = 0; iVar = iVar_Sens;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
	}
  
	/*--- Immediately release the temporary buffers. ---*/
  
	delete [] Buffer_Send_Var;
	delete [] Buffer_Send_Res;
	delete [] Buffer_Send_Vol;
	delete [] Buffer_Send_GlobalIndex;
	if (rank == MASTER_NODE) {
		delete [] Buffer_Recv_Var;
		delete [] Buffer_Recv_Res;
		delete [] Buffer_Recv_Vol;
		delete [] Buffer_Recv_GlobalIndex;
	}
  
#endif
  
  /*--- Release memory needed for surface coefficients ---*/
  
  delete [] Aux_Press;
	if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
    delete [] Aux_Frict; delete [] Aux_Heat; delete [] Aux_yPlus;
  }
  if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) ||
			(Kind_Solver == ADJ_FREE_SURFACE_EULER) || (Kind_Solver == ADJ_FREE_SURFACE_NAVIER_STOKES) ||
			(Kind_Solver == ADJ_FREE_SURFACE_RANS) || (Kind_Solver == ADJ_PLASMA_EULER) || (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES)) {
    delete [] Aux_Sens;
  }
  
}

void COutput::MergeBaselineSolution(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone) {
  
  /*--- Local variables needed on all processors ---*/
	unsigned short iVar;
	unsigned long iPoint = 0, jPoint = 0;
  
  nVar_Total = config->fields.size() - 1;
  
  /*--- Merge the solution either in serial or parallel. ---*/
#ifdef NO_MPI
  
	/*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/
  nGlobal_Poin = geometry->GetnPointDomain();
	Data = new double*[nVar_Total];
	for (iVar = 0; iVar < nVar_Total; iVar++) {
		Data[iVar] = new double[nGlobal_Poin];
	}
  
	/*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (ignore periodic/sliding halo nodes). ---*/
  jPoint = 0;
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		if (geometry->node[iPoint]->GetDomain()) {
			/*--- Solution (first, and second system of equations) ---*/
			unsigned short jVar = 0;
			for (iVar = 0; iVar < nVar_Total; iVar++) {
				Data[jVar][jPoint] = solver->node[iPoint]->GetSolution(iVar);
        jVar++;
			}
		}
    
    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/
    jPoint++;
    
	}
  
#else
  
	/*--- MPI preprocessing ---*/
  
  int rank = MPI::COMM_WORLD.Get_rank();
	int nProcessor = MPI::COMM_WORLD.Get_size();
	int iProcessor;
  
	/*--- Local variables needed for merging with MPI ---*/
	unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
	unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
	unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  
  /*--- Each processor sends its local number of nodes to the master. ---*/
	//!
	//! TO DO: MPI I/O for writing the solution files.
	//!
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else
    nLocalPoint = geometry->GetnPointDomain();
	Buffer_Send_nPoint[0] = nLocalPoint;
	if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nPoint, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nPoint, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
	nBuffer_Scalar = MaxLocalPoint;
  
	/*--- Send and Recv buffers. ---*/
  
	double *Buffer_Send_Var = new double[MaxLocalPoint];
	double *Buffer_Recv_Var = NULL;
  
	unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
	unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
	/*--- Prepare the receive buffers in the master node only. ---*/
	if (rank == MASTER_NODE) {
    
		Buffer_Recv_Var = new double[nProcessor*MaxLocalPoint];
		Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
		/*--- Sum total number of nodes to be written and allocate arrays ---*/
		nGlobal_Poin = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
			nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
		}
		Data = new double*[nVar_Total];
		for (iVar = 0; iVar < nVar_Total; iVar++) {
			Data[iVar] = new double[nGlobal_Poin];
		}
    
  }
  
	/*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
  
	for (iVar = 0; iVar < nVar_Total; iVar++) {
    
		/*--- Loop over this partition to collect the current variable ---*/
		jPoint = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			
      /*--- Check for halos and write only if requested ---*/
      if (geometry->node[iPoint]->GetDomain() || Wrt_Halo) {

      /*--- Get this variable into the temporary send buffer. ---*/
				Buffer_Send_Var[jPoint] = solver->node[iPoint]->GetSolution(iVar);

      /*--- Only send/recv the volumes & global indices during the first loop ---*/
				if (iVar == 0) {
					Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
				}
				jPoint++;
			}
		}
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		if (iVar == 0) {
			MPI::COMM_WORLD.Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                             Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                             MASTER_NODE);
		}
    
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
      jPoint = 0;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
	}
  
	/*--- Immediately release the temporary buffers. ---*/
  
	delete [] Buffer_Send_Var;
	delete [] Buffer_Send_GlobalIndex;
	if (rank == MASTER_NODE) {
		delete [] Buffer_Recv_Var;
		delete [] Buffer_Recv_GlobalIndex;
	}
  
#endif
  
}

void COutput::SetRestart(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {

	/*--- Local variables ---*/
	unsigned short iVar;
	unsigned long iPoint;
  unsigned short iDim, nDim = geometry->GetnDim();
	ofstream restart_file;
	string filename, AdjExt, UnstExt;
	unsigned long iExtIter = config->GetExtIter();
	char buffer[50];
	unsigned short Kind_ObjFunc = config->GetKind_ObjFunc();
  unsigned short Kind_Solver  = config->GetKind_Solver();
  bool grid_movement = ((config->GetUnsteady_Simulation() &&
                         config->GetWrt_Unsteady() &&
                         config->GetGrid_Movement()) ||
                        ((config->GetUnsteady_Simulation() == TIME_SPECTRAL) &&
                         config->GetGrid_Movement()));
  
	/*--- Retrieve filename from config ---*/
	if (config->GetAdjoint())
		filename = config->GetRestart_AdjFileName();
	else
		filename = config->GetRestart_FlowFileName();

	/*--- Remove restart filename extension ---*/
	filename.erase(filename.end()-4, filename.end());

	/*--- The adjoint problem requires a particular extension. ---*/
	if (config->GetAdjoint()) {
		switch (Kind_ObjFunc) {
		case DRAG_COEFFICIENT:      AdjExt = "_cd";   break;
		case LIFT_COEFFICIENT:      AdjExt = "_cl";   break;
		case SIDEFORCE_COEFFICIENT: AdjExt = "_csf";  break;
		case PRESSURE_COEFFICIENT:  AdjExt = "_cp";   break;
		case MOMENT_X_COEFFICIENT:  AdjExt = "_cmx";  break;
		case MOMENT_Y_COEFFICIENT:  AdjExt = "_cmy";  break;
		case MOMENT_Z_COEFFICIENT:  AdjExt = "_cmz";  break;
		case EFFICIENCY:            AdjExt = "_eff";  break;
		case EQUIVALENT_AREA:       AdjExt = "_ea";   break;
		case NEARFIELD_PRESSURE:    AdjExt = "_nfp";  break;
		case FORCE_X_COEFFICIENT:   AdjExt = "_cfx";  break;
		case FORCE_Y_COEFFICIENT:   AdjExt = "_cfy";  break;
		case FORCE_Z_COEFFICIENT:   AdjExt = "_cfz";  break;
		case THRUST_COEFFICIENT:    AdjExt = "_ct";   break;
		case TORQUE_COEFFICIENT:    AdjExt = "_cq";   break;
    case HEAT_LOAD:             AdjExt = "_Q";    break;
    case MAX_HEAT_FLUX:         AdjExt = "_qmax"; break;
		case FIGURE_OF_MERIT:       AdjExt = "_merit";break;
		case FREE_SURFACE:          AdjExt = "_fs";   break;
		case NOISE:                 AdjExt = "_fwh";  break;
		}
		filename.append(AdjExt);
	}

	/*--- Unsteady problems require the physical timestep to be appended. ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {

		if (int(val_iZone) < 10) sprintf (buffer, "_0000%d", int(val_iZone));
		if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d", int(val_iZone));
		if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d", int(val_iZone));
		if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d", int(val_iZone));
		if (int(val_iZone) >= 10000) sprintf (buffer, "_%d", int(val_iZone));
		UnstExt = string(buffer);
		filename.append(UnstExt);
	} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		if ((int(iExtIter) >= 0) && (int(iExtIter) < 10))
			sprintf (buffer, "_0000%d", int(iExtIter));
		if ((int(iExtIter) >= 10) && (int(iExtIter) < 100))
			sprintf (buffer, "_000%d", int(iExtIter));
		if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000))
			sprintf (buffer, "_00%d", int(iExtIter));
		if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000))
			sprintf (buffer, "_0%d", int(iExtIter));
		if (int(iExtIter) >= 10000)
			sprintf (buffer, "_%d", int(iExtIter));
		UnstExt = string(buffer);
		filename.append(UnstExt);
	}

	/*--- Lastly, add the .dat extension ---*/
	filename.append(".dat");

	/*--- Open the restart file and write the solution. ---*/
	restart_file.open(filename.c_str(), ios::out);
	restart_file.precision(15);


	/*--- Write the basic header based on the particular solution ----*/
	restart_file << "\"PointID\"";
	for (iVar = 0; iVar < nVar_Consv; iVar++) {
		restart_file << ", \"Conservative_" << iVar+1<<"\"";
	}
    if (config->GetWrt_Residuals()) {
        for (iVar = 0; iVar < nVar_Consv; iVar++) {
            restart_file << ", \"Residual_" << iVar+1<<"\"";
        }
    }

  /*--- Add names for any extra variables (this will need to be adjusted). ---*/
	if (grid_movement) {
    if (nDim == 2) {
      restart_file << ", \"Grid_Velx\", \"Grid_Vely\"";
    } else {
      restart_file << ", \"Grid_Velx\", \"Grid_Vely\", \"Grid_Velz\"";
    }
	}

  if ((Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
    restart_file << ", \"Density\"";
  }
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
    restart_file << ", \"Pressure\", \"Pressure_Coefficient\", \"Mach\"";
  }

  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
    restart_file << ", \"Temperature\", \"Laminar_Viscosity\", \"Skin_Friction_Coefficient\", \"Heat_Transfer\", \"Y_Plus\"";
  }

  if ((Kind_Solver == RANS) || (Kind_Solver == FREE_SURFACE_RANS)) {
    restart_file << ", \"Eddy_Viscosity\"";
  }

  if ((Kind_Solver == PLASMA_EULER) || (Kind_Solver == PLASMA_NAVIER_STOKES)) {
    unsigned short iSpecies;
    for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
      restart_file << ", \"Pressure_" << iSpecies << "\"";
    for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
      restart_file << ", \"Temperature_" << iSpecies << "\"";
    for (iSpecies = 0; iSpecies < config->GetnDiatomics(); iSpecies++)
      restart_file << ", \"TemperatureVib_" << iSpecies << "\"";
    for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
      restart_file << ", \"Mach_" << iSpecies << "\"";
  }
  
  if (Kind_Solver == PLASMA_NAVIER_STOKES) {
    unsigned short iSpecies;
    for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
      restart_file << ", \"LaminaryViscosity_" << iSpecies << "\"";
    
    if ( Kind_Solver == PLASMA_NAVIER_STOKES  && (config->GetMagnetic_Force() == YES) && (geometry->GetnDim() == 3)) {
      for (iDim = 0; iDim < nDim; iDim++)
        restart_file << ", \"Magnet_Field" << iDim << "\"";
    }
  }
  
  if (Kind_Solver == ELECTRIC_POTENTIAL) {
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
      restart_file << ", \"ElectricField_" << iDim+1 << "\"";
  }
  
  if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) ||
			(Kind_Solver == ADJ_FREE_SURFACE_EULER) || (Kind_Solver == ADJ_FREE_SURFACE_NAVIER_STOKES) ||
			(Kind_Solver == ADJ_FREE_SURFACE_RANS) || (Kind_Solver == ADJ_PLASMA_EULER) || (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES)) {
    restart_file << ", \"Surface_Sensitivity\"";
  }
  
  restart_file << endl;
  
  
	//  for (iPoint = 0; iPoint < nGlobal_Poin; iPoint++) { (F.P.)
	for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {

		/*--- Index of the point ---*/
		restart_file << iPoint << "\t";

		/*--- Loop over the variables and write the values to file ---*/
		for (iVar = 0; iVar < nVar_Total; iVar++) {
			restart_file << scientific << Data[iVar][iPoint] << "\t";
		}

		restart_file << endl;
	}

  restart_file.close();
  
	/*--- Binary version - eventually ---*/
	//  restart_file.close();
	//  restart_file.open("restart_file.bin", ios::binary);
	//  restart_file.precision(15);
	//  restart_file.write( (char *)&Data[0][0], sizeof(double)*(nGlobal_Poin*nVar_Consv*2) );
	//  restart_file.write( (char *)&Volume[0], sizeof(double)*(nGlobal_Poin) );
	//  restart_file.close();

}

void COutput::DeallocateCoordinates(CConfig *config, CGeometry *geometry) {

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	/*--- Local variables and initialization ---*/

	unsigned short iDim, nDim = geometry->GetnDim();

	/*--- The master node alone owns all data found in this routine. ---*/
	if (rank == MASTER_NODE) {

		/*--- Deallocate memory for coordinate data ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			delete [] Coords[iDim];
		}
		delete [] Coords;

	}
}

void COutput::DeallocateConnectivity(CConfig *config, CGeometry *geometry, bool surf_sol) {

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	/*--- The master node alone owns all data found in this routine. ---*/
	if (rank == MASTER_NODE) {

		/*--- Deallocate memory for connectivity data ---*/
    if (surf_sol) {
      if (nGlobal_Line > 0) delete [] Conn_Line;
      if (nGlobal_BoundTria > 0) delete [] Conn_BoundTria;
      if (nGlobal_BoundQuad > 0) delete [] Conn_BoundQuad;
    }
    else {
      if (nGlobal_Tria > 0) delete [] Conn_Tria;
      if (nGlobal_Quad > 0) delete [] Conn_Quad;
      if (nGlobal_Tetr > 0) delete [] Conn_Tetr;
      if (nGlobal_Hexa > 0) delete [] Conn_Hexa;
      if (nGlobal_Wedg > 0) delete [] Conn_Wedg;
      if (nGlobal_Pyra > 0) delete [] Conn_Pyra;
    }

	}
}

void COutput::DeallocateSolution(CConfig *config, CGeometry *geometry) {

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	/*--- Local variables and initilization ---*/

	unsigned short iVar;

	/*--- The master node alone owns all data found in this routine. ---*/
	if (rank == MASTER_NODE) {

		/*--- Deallocate memory for solution data ---*/
		for (iVar = 0; iVar < nVar_Total; iVar++) {
			delete [] Data[iVar];
		}
		delete [] Data;

		/*--- Deallocate memory for volume data (needed in restart) ---*/
		//delete [] Volume;
	}

}

void COutput::SetHistory_Header(ofstream *ConvHist_file, CConfig *config) {
	char cstr[200], buffer[50];

	bool rotating_frame = config->GetRotating_Frame();
	bool equiv_area = config->GetEquivArea();
	bool turbulent = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS));
  
  bool isothermal = false;
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) isothermal = true;
  
	unsigned short iSpecies;

	/*--- Write file name with extension ---*/
	strcpy (cstr, config->GetConv_FileName().data());
	if (config->GetOutput_FileFormat() == TECPLOT)  sprintf (buffer, ".plt");
	if (config->GetOutput_FileFormat() == TECPLOT_BINARY)  sprintf (buffer, ".plt");
	if (config->GetOutput_FileFormat() == CGNS_SOL)  sprintf (buffer, ".csv");
	strcat(cstr,buffer);

	ConvHist_file->open(cstr, ios::out);
	ConvHist_file->precision(15);

	char begin[]= "\"Iteration\"";

	char flow_coeff[]= ",\"CLift\",\"CDrag\",\"CSideForce\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\"";
	char heat_coeff[]= ",\"CHeat_Load\",\"CHeat_Max\"";
	char equivalent_area_coeff[]= ",\"CEquivArea\",\"CNearFieldOF\"";
	char rotating_frame_coeff[]= ",\"CMerit\",\"CT\",\"CQ\"";
	char free_surface_coeff[]= ",\"CFreeSurface\"";
	char plasma_coeff[]= ",\"CLift\",\"CDrag\",\"CSideForce\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\",\"Q\",\"PressDrag\",\"ViscDrag\",\"MagnetDrag\"";
	char wave_coeff[]= ",\"CWave\"";
	char fea_coeff[]= ",\"CFEA\"";
  char adj_coeff[]= ",\"Sens_Geo\",\"Sens_Mach\",\"Sens_AoA\",\"Sens_Press\",\"Sens_Temp\",\"Sens_AoS\"";
  char adj_plasma_coeff[]= ",\"Sens_Geo\",\"Sens_Mach\",\"Sens_AoA\",\"Sens_Press\",\"Sens_Temp\",\"Sens_AoS\"";
  
	char flow_resid[]= ",\"Res_Flow[0]\",\"Res_Flow[1]\",\"Res_Flow[2]\",\"Res_Flow[3]\",\"Res_Flow[4]\"";
	char adj_flow_resid[]= ",\"Res_AdjFlow[0]\",\"Res_AdjFlow[1]\",\"Res_AdjFlow[2]\",\"Res_AdjFlow[3]\",\"Res_AdjFlow[4]\"";
  
	char turb_resid[1000];
	switch (config->GetKind_Turb_Model()){
    case SA:	sprintf (turb_resid, ",\"Res_Turb[0]\""); break;
    case SST:	sprintf (turb_resid, ",\"Res_Turb[0]\",\"Res_Turb[1]\""); break;
	}

	char adj_turb_resid[]= ",\"Res_AdjTurb[0]\"";

	char levelset_resid[]= ",\"Res_LevelSet\"";
	char adj_levelset_resid[]= ",\"Res_AdjLevelSet\"";

	char wave_resid[]= ",\"Res_Wave[0]\",\"Res_Wave[1]\"";

	char fea_resid[]= ",\"Res_FEA\"";

	char end[]= ",\"Linear_Solver_Iterations\",\"Time(min)\"\n";

	if (config->GetOutput_FileFormat() == TECPLOT) {
		char title[]= "TITLE = \"SU2 Simulation\"";
		char variables[]= "VARIABLES = ";
		ConvHist_file[0] << title << endl;
		ConvHist_file[0] << variables;
	}

	switch (config->GetKind_Solver()) {
	case EULER : case NAVIER_STOKES: case RANS :
		ConvHist_file[0] << begin << flow_coeff;
    if (isothermal) ConvHist_file[0] << heat_coeff;
		if (equiv_area) ConvHist_file[0] << equivalent_area_coeff;
		if (rotating_frame) ConvHist_file[0] << rotating_frame_coeff;
		ConvHist_file[0] << flow_resid;
		if (turbulent) ConvHist_file[0] << turb_resid;
		ConvHist_file[0] << end;
		break;

	case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
		ConvHist_file[0] << begin << flow_coeff << free_surface_coeff;
		ConvHist_file[0] << flow_resid << levelset_resid << end;
		break;

	case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
		ConvHist_file[0] << begin << plasma_coeff;
		for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
			ConvHist_file[0] << ",\"Res_Density[" << iSpecies << "]\",\"Res_Energy[" << iSpecies << "]\"";
		}
		ConvHist_file[0] << end;
		break;

	case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS:
		ConvHist_file[0] << begin << adj_coeff << adj_flow_resid;
        if (turbulent)
            if (((config->GetKind_Adjoint() != HYBRID) && (!config->GetFrozen_Visc())) || (config->GetKind_Adjoint() == HYBRID) )
                ConvHist_file[0] << adj_turb_resid;
		ConvHist_file[0] << end;
		break;

	case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
		ConvHist_file[0] << begin << adj_coeff << adj_flow_resid << adj_levelset_resid << end;
		break;

	case ADJ_PLASMA_EULER : case ADJ_PLASMA_NAVIER_STOKES:
		ConvHist_file[0] << begin << adj_plasma_coeff;
		for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
			ConvHist_file[0] << ",\"Res_PsiDensity[" << iSpecies << "]\",\"Res_PsiEnergy[" << iSpecies << "]\"";
		}
		ConvHist_file[0] << end;
		break;

	case WAVE_EQUATION:
		ConvHist_file[0] << begin << wave_coeff;
		ConvHist_file[0] << wave_resid << end;
		break;

	case LINEAR_ELASTICITY:
		ConvHist_file[0] << begin << fea_coeff;
		ConvHist_file[0] << fea_resid << end;
		break;

	case AEROACOUSTIC_EULER:
		ConvHist_file[0] << begin << flow_coeff;
		ConvHist_file[0] << wave_coeff;
		ConvHist_file[0] << flow_resid;
		ConvHist_file[0] << wave_resid << end;
		break;

	case ADJ_AEROACOUSTIC_EULER:
		ConvHist_file[0] << begin << adj_coeff;
		ConvHist_file[0] << adj_flow_resid;
		ConvHist_file[0] << wave_coeff;
		ConvHist_file[0] << wave_resid << end;
		break;
	}

	if (config->GetOutput_FileFormat() == TECPLOT || config->GetOutput_FileFormat() == TECPLOT_BINARY) {
		char zone[]= "ZONE T= \"Convergence history\"";
		ConvHist_file[0] << zone << endl;
	}

}


void COutput::SetConvergence_History(ofstream *ConvHist_file, CGeometry ***geometry, CSolver ****solver_container, CConfig **config, CIntegration ***integration, bool DualTime_Iteration, unsigned long timeused, unsigned short val_iZone) {
  
	unsigned long iIntIter = config[val_iZone]->GetIntIter();
	unsigned long iExtIter = config[val_iZone]->GetExtIter();
  
#ifndef NO_MPI
  
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
  
#else
  
	int rank = MASTER_NODE;
  
#endif
  
	/*--- WARNING: These buffers have hard-coded lengths. Note that you
   may have to adjust them to be larger if adding more entries. ---*/
	char begin[1000], direct_coeff[1000], adjoint_coeff[1000], flow_resid[1000], adj_flow_resid[1000],
	turb_resid[1000], trans_resid[1000], adj_turb_resid[1000], plasma_resid[1000], adj_plasma_resid[1000], resid_aux[1000],
	levelset_resid[1000], adj_levelset_resid[1000], wave_coeff[1000], fea_coeff[1000], wave_resid[1000],
	fea_resid[1000], end[1000];
	double dummy = 0.0;
	unsigned short iVar;
  
	unsigned long LinSolvIter = 0;
	double timeiter = double(timeused)/double(iExtIter+1);
  
	unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
	unsigned short nDim = geometry[val_iZone][FinestMesh]->GetnDim();
	unsigned short nSpecies = config[val_iZone]->GetnSpecies();
  
	bool incompressible = config[val_iZone]->GetIncompressible();
	bool compressible = !config[val_iZone]->GetIncompressible();
	bool rotating_frame = config[val_iZone]->GetRotating_Frame();
	bool equiv_area = config[val_iZone]->GetEquivArea();
	bool transition = (config[val_iZone]->GetKind_Trans_Model() == LM);
  bool isothermal = false;
  for (unsigned short iMarker = 0; iMarker < config[val_iZone]->GetnMarker_All(); iMarker++)
    if (config[val_iZone]->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) isothermal = true;
  
  	bool turbulent = ((config[val_iZone]->GetKind_Solver() == RANS) || (config[val_iZone]->GetKind_Solver() == ADJ_RANS) ||
                    (config[val_iZone]->GetKind_Solver() == FREE_SURFACE_RANS) || (config[val_iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_RANS) ||
                    (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_RANS));
	bool adjoint = config[val_iZone]->GetAdjoint();
	bool fluid_structure = ((config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_EULER) || (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_NAVIER_STOKES) ||
                          (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_RANS));
    bool freesurface = config[val_iZone]->GetFreeSurface();
	bool aeroacoustic = ((config[val_iZone]->GetKind_Solver() == AEROACOUSTIC_EULER) || (config[val_iZone]->GetKind_Solver() == AEROACOUSTIC_NAVIER_STOKES) ||
                       (config[val_iZone]->GetKind_Solver() == AEROACOUSTIC_RANS));
	bool wave = (config[val_iZone]->GetKind_Solver() == WAVE_EQUATION);
	bool fea = (config[val_iZone]->GetKind_Solver() == LINEAR_ELASTICITY);
	bool plasma = ((config[val_iZone]->GetKind_Solver() == PLASMA_EULER) || (config[val_iZone]->GetKind_Solver() == PLASMA_NAVIER_STOKES) ||
                 (config[val_iZone]->GetKind_Solver() == ADJ_PLASMA_EULER) || (config[val_iZone]->GetKind_Solver() == ADJ_PLASMA_NAVIER_STOKES));
  
	/*--- Initialize variables to store information from all domains (direct solution) ---*/
  double Total_CLift = 0.0, Total_CDrag = 0.0, Total_CSideForce = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0,
  Total_CEquivArea = 0.0, Total_CNearFieldOF = 0.0, Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0, Total_CMerit = 0.0,
  Total_CT = 0.0, Total_CQ = 0.0, Total_CFreeSurface = 0.0, Total_CWave = 0.0, Total_CFEA = 0.0, PressureDrag = 0.0, ViscDrag = 0.0, MagDrag = 0.0, Total_Q = 0.0, Total_MaxQ = 0.0;
  
  /*--- Initialize variables to store information from all domains (adjoint solution) ---*/
  double Total_Sens_Geo = 0.0, Total_Sens_Mach = 0.0, Total_Sens_AoA = 0.0;
  double Total_Sens_Press = 0.0, Total_Sens_Temp = 0.0;
  
  /*--- Residual arrays ---*/
  double *residual_flow = NULL, *residual_turbulent = NULL, *residual_transition = NULL, *residual_levelset = NULL, *residual_plasma = NULL;
  double *residual_adjflow = NULL, *residual_adjturbulent = NULL, *residual_adjlevelset = NULL, *residual_adjplasma = NULL;
  double *residual_wave = NULL; double *residual_fea = NULL;
  
  /*--- Forces arrays ---*/
  double *sbuf_force = NULL,  *rbuf_force = NULL;
  unsigned long *sbuf_time = NULL, *rbuf_time = NULL;
  
  /*--- Initialize number of variables ---*/
  unsigned short nVar_Flow = 0, nVar_LevelSet = 0, nVar_Turb = 0, nVar_Trans = 0, nVar_Wave = 0, nVar_FEA = 0, nVar_Plasma = 0,
  nVar_AdjFlow = 0, nVar_AdjPlasma = 0, nVar_AdjLevelSet = 0, nVar_AdjTurb = 0, nVar_Force = 0, nVar_Time = 0;
  
  /*--- Direct problem variables ---*/
  if (compressible) nVar_Flow = nDim+2; else nVar_Flow = nDim+1;
  if (turbulent)
    switch (config[val_iZone]->GetKind_Turb_Model()){
      case SA:	nVar_Turb = 1; break;
      case SST: nVar_Turb = 2; break;
    }
  if (transition) nVar_Trans = 2;
  if (wave) nVar_Wave = 2;
  if (fea) nVar_FEA = nDim;
  if (plasma) nVar_Plasma = config[val_iZone]->GetnMonatomics()*(nDim+2) + config[val_iZone]->GetnDiatomics()*(nDim+3);
  if (freesurface) nVar_LevelSet = 1;
  
  /*--- Adjoint problem variables ---*/
  if (compressible) nVar_AdjFlow = nDim+2; else nVar_AdjFlow = nDim+1;
  if (turbulent)
    switch (config[val_iZone]->GetKind_Turb_Model()){
      case SA:	nVar_AdjTurb = 1; break;
      case SST: nVar_AdjTurb = 2; break;
    }
  if (plasma) nVar_AdjPlasma = config[val_iZone]->GetnMonatomics()*(nDim+2) + config[val_iZone]->GetnDiatomics()*(nDim+3);
  if (freesurface) nVar_AdjLevelSet = 1;
  
  /*--- Other vectors ---*/
  nVar_Force = 18; nVar_Time = 3;
  
  /*--- Allocate memory for send buffer ---*/
  sbuf_force = new double[nVar_Force];
  sbuf_time = new unsigned long[nVar_Time];
  
  for (iVar = 0; iVar < nVar_Force; iVar++) { sbuf_force[iVar] = 0.0; }
  for (iVar = 0; iVar < nVar_Time; iVar++) { sbuf_time[iVar] = 0; }
  
  /*--- Allocate memory for the residual ---*/
  residual_flow = new double[nVar_Flow];
  residual_turbulent = new double[nVar_Turb];
  residual_transition = new double[nVar_Trans];
  residual_plasma = new double[nVar_Plasma];
  residual_levelset = new double[nVar_LevelSet];
  residual_wave = new double[nVar_Wave];
  residual_fea = new double[nVar_FEA];
  
  residual_adjflow = new double[nVar_AdjFlow];
  residual_adjturbulent = new double[nVar_AdjTurb];
  residual_adjplasma = new double[nVar_AdjPlasma];
  residual_adjlevelset = new double[nVar_AdjLevelSet];
  
  /*--- Allocate a vector for the forces, and the time ---*/
  if (rank == MASTER_NODE) {
    
    rbuf_force = new double[nVar_Force];
    rbuf_time = new unsigned long[nVar_Time];
    
    for (iVar = 0; iVar < nVar_Force; iVar++) rbuf_force[iVar] = 0.0;
    for (iVar = 0; iVar < nVar_Time; iVar++) rbuf_time[iVar] = 0;
    
  }
  
  /*--- Write information from nodes ---*/
  switch (config[val_iZone]->GetKind_Solver()) {
      
    case EULER:                   case NAVIER_STOKES:                   case RANS:
    case FREE_SURFACE_EULER:      case FREE_SURFACE_NAVIER_STOKES:      case FREE_SURFACE_RANS:
    case FLUID_STRUCTURE_EULER:   case FLUID_STRUCTURE_NAVIER_STOKES:   case FLUID_STRUCTURE_RANS:
    case AEROACOUSTIC_EULER:      case AEROACOUSTIC_NAVIER_STOKES:      case AEROACOUSTIC_RANS:
    case ADJ_EULER:               case ADJ_NAVIER_STOKES:               case ADJ_RANS:
    case ADJ_FREE_SURFACE_EULER:  case ADJ_FREE_SURFACE_NAVIER_STOKES:  case ADJ_FREE_SURFACE_RANS:
    case ADJ_AEROACOUSTIC_EULER:
      
      /*--- Flow solution coefficients (serial) ---*/
      Total_CLift += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CLift();
      Total_CDrag += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CDrag();
      Total_CSideForce += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CSideForce();
      Total_CMx += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMx();
      Total_CMy += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMy();
      Total_CMz += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMz();
      Total_CFx += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFx();
      Total_CFy += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFy();
      Total_CFz += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFz();
      Total_CEff += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEff();
      
      if (freesurface) {
        Total_CFreeSurface += solver_container[val_iZone][FinestMesh][LEVELSET_SOL]->GetTotal_CFreeSurface();
      }
      
      if (isothermal) {
        Total_Q += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_Q();
        Total_MaxQ += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_MaxQ();
      }
      
      if (equiv_area) {
        Total_CEquivArea += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea();
        Total_CNearFieldOF += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();
      }
      
      if (rotating_frame) {
        Total_CMerit += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMerit();
        Total_CT += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CT();
        Total_CQ += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CQ();
      }
      
      if (aeroacoustic) {
        Total_CWave += solver_container[ZONE_1][FinestMesh][WAVE_SOL]->GetTotal_CWave();
      }
      
      if (fluid_structure) {
        Total_CFEA  += solver_container[ZONE_1][FinestMesh][FEA_SOL]->GetTotal_CFEA();
      }
      
      /*--- Flow solution coefficients (parallel) ---*/
      sbuf_force[0] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CLift();
      sbuf_force[1] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CDrag();
      sbuf_force[2] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CSideForce();
      sbuf_force[3] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMx();
      sbuf_force[4] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMy();
      sbuf_force[5] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMz();
      sbuf_force[6] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFx();
      sbuf_force[7] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFy();
      sbuf_force[8] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFz();
      
      if (freesurface) {
        sbuf_force[9] += solver_container[val_iZone][FinestMesh][LEVELSET_SOL]->GetTotal_CFreeSurface();
      }
      
      if (fluid_structure) {
        sbuf_force[9] += solver_container[ZONE_1][FinestMesh][FEA_SOL]->GetTotal_CFEA();
      }
      
      if (aeroacoustic) {
        sbuf_force[9]  += solver_container[ZONE_1][FinestMesh][WAVE_SOL]->GetTotal_CWave();
      }
      
      if (equiv_area) {
        sbuf_force[9] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea(),
        sbuf_force[10] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();
      }
      
      if (rotating_frame) {
        sbuf_force[9]  += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMerit();
        sbuf_force[10] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CT();
        sbuf_force[11] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CQ();
      }
      
      if (isothermal) {
        sbuf_force[12] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_Q();
        sbuf_force[13] += solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_MaxQ();
      }
      
      /*--- Flow Residuals ---*/
      for (iVar = 0; iVar < nVar_Flow; iVar++)
        residual_flow[iVar] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_RMS(iVar);
      
      /*--- Turbulent residual ---*/
      if (turbulent) {
        for (iVar = 0; iVar < nVar_Turb; iVar++)
          residual_turbulent[iVar] = solver_container[val_iZone][FinestMesh][TURB_SOL]->GetRes_RMS(iVar);
      }
      
      /*--- Transition residual ---*/
      if (transition) {
        for (iVar = 0; iVar < nVar_Trans; iVar++)
          residual_transition[iVar] = solver_container[val_iZone][FinestMesh][TRANS_SOL]->GetRes_RMS(iVar);
      }
      
      /*--- Free Surface residual ---*/
      if (freesurface) {
        for (iVar = 0; iVar < nVar_LevelSet; iVar++)
          residual_levelset[iVar] = solver_container[val_iZone][FinestMesh][LEVELSET_SOL]->GetRes_RMS(iVar);
      }
      
      /*--- FEA residual ---*/
      if (fluid_structure) {
        for (iVar = 0; iVar < nVar_FEA; iVar++)
          residual_fea[iVar] = solver_container[ZONE_1][FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
      }
      
      /*--- Aeroacoustic residual ---*/
      if (aeroacoustic) {
        for (iVar = 0; iVar < nVar_Wave; iVar++)
          residual_wave[iVar] = solver_container[ZONE_1][FinestMesh][WAVE_SOL]->GetRes_RMS(iVar);
      }
      
      /*--- Iterations of the linear solver ---*/
      LinSolvIter = (unsigned long) solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetIterLinSolver();
      
      /*--- Adjoint solver ---*/
      if (adjoint) {
        
        /*--- Adjoint solution coefficients (serial) ---*/
        Total_Sens_Geo  = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Geo();
        Total_Sens_Mach = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Mach();
        Total_Sens_AoA  = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_AoA();
        Total_Sens_Press = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Press();
        Total_Sens_Temp  = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Temp();
        
        /*--- Adjoint solution coefficients (parallel) ---*/
        sbuf_force[0] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Geo();
        sbuf_force[1] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Mach();
        sbuf_force[2] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_AoA();
        sbuf_force[3] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Press();
        sbuf_force[4] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Temp();
        
        /*--- Adjoint flow residuals ---*/
        for (iVar = 0; iVar < nVar_AdjFlow; iVar++) {
          residual_adjflow[iVar] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(iVar);
        }
        
        /*--- Adjoint turbulent residuals ---*/
        if (turbulent)
            if (((config[val_iZone]->GetKind_Adjoint() != HYBRID) && (!config[val_iZone]->GetFrozen_Visc())) || (config[val_iZone]->GetKind_Adjoint() == HYBRID)) {
                for (iVar = 0; iVar < nVar_AdjTurb; iVar++)
                    residual_adjturbulent[iVar] = solver_container[val_iZone][FinestMesh][ADJTURB_SOL]->GetRes_RMS(iVar);
            }
        
        /*--- Adjoint level set residuals ---*/
        if (freesurface) {
          for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++)
            residual_adjlevelset[iVar] = solver_container[val_iZone][FinestMesh][ADJLEVELSET_SOL]->GetRes_RMS(iVar);
        }
        
      }
      
      break;
      
    case PLASMA_EULER:      case PLASMA_NAVIER_STOKES:
    case ADJ_PLASMA_EULER:  case ADJ_PLASMA_NAVIER_STOKES:
      
      /*--- Plasma coefficients (serial) ---*/
      Total_CLift += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CLift();
      Total_CDrag += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CDrag();
      Total_CSideForce += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CSideForce();
      Total_CMx += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMx();
      Total_CMy += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMy();
      Total_CMz += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMz();
      Total_CFx += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFx();
      Total_CFy += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFy();
      Total_CFz += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFz();
      Total_CEff += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CEff();
      Total_Q += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_Q();
      PressureDrag += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_PressureDrag();
      ViscDrag += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_ViscDrag();
      MagDrag  += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_MagnetDrag();
      
      /*--- Plasma coefficients (parallel) ---*/
      sbuf_force[0] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CLift();
      sbuf_force[1] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CDrag();
      sbuf_force[2] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CSideForce();
      sbuf_force[3] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMx();
      sbuf_force[4] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMy();
      sbuf_force[5] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMz();
      sbuf_force[6] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFx();
      sbuf_force[7] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFy();
      sbuf_force[8] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFz();
      sbuf_force[9] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CQ();
      sbuf_force[10] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_PressureDrag();
      sbuf_force[11] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_ViscDrag();
      sbuf_force[12] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_MagnetDrag();
      sbuf_force[13] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_Q();
      
      /*--- Plasma Residuals ---*/
      for (iVar = 0; iVar < nVar_Plasma; iVar++)
        residual_plasma[iVar] = solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetRes_RMS(iVar);
      
      if (adjoint) {
        
        /*--- Adjoint solution coefficients (serial) ---*/
        Total_Sens_Geo  = 0.0; //solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Geo();
        Total_Sens_Mach = 0.0; //solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Mach();
        Total_Sens_AoA  = 0.0; //solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_AoA();
        Total_Sens_Press = 0.0; //solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Press();
        Total_Sens_Temp  = 0.0; //solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Temp();
        
        /*--- Adjoint solution coefficients (parallel) ---*/
        sbuf_force[0] = 0.0; //solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Geo();
        sbuf_force[1] = 0.0; //solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Mach();
        sbuf_force[2] = 0.0; //solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_AoA();
        sbuf_force[3] = 0.0; //solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Press();
        sbuf_force[4] = 0.0; //solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Temp();
        
        /*--- Adjoint plasma residuals ---*/
        for (iVar = 0; iVar < nVar_AdjPlasma; iVar++)
          residual_adjplasma[iVar] = solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetRes_RMS(iVar);
      }
      
      break;
      
    case WAVE_EQUATION:
      
      /*--- Wave coefficients (serial) ---*/
      Total_CWave += solver_container[val_iZone][FinestMesh][WAVE_SOL]->GetTotal_CWave();
      
      /*--- Wave coefficients (parallel) ---*/
      sbuf_force[0] += solver_container[val_iZone][FinestMesh][WAVE_SOL]->GetTotal_CWave();
      
      /*--- Wave Residuals ---*/
      for (iVar = 0; iVar < nVar_Wave; iVar++) {
        residual_wave[iVar] = solver_container[val_iZone][FinestMesh][WAVE_SOL]->GetRes_RMS(iVar);
      }
      
      break;
      
    case LINEAR_ELASTICITY:
      
      /*--- FEA coefficients (serial) ---*/
      Total_CFEA += solver_container[val_iZone][FinestMesh][FEA_SOL]->GetTotal_CFEA();
      
      /*--- FEA coefficients (parallel) ---*/
      sbuf_force[0] += solver_container[val_iZone][FinestMesh][FEA_SOL]->GetTotal_CFEA();
      
      /*--- Plasma Residuals ---*/
      for (iVar = 0; iVar < nVar_FEA; iVar++) {
        residual_fea[iVar] = solver_container[val_iZone][FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
      }
      
      break;
      
  }
  
  sbuf_time[0] = (unsigned long) timeiter;
  sbuf_time[1] = (unsigned long) timeused;
  sbuf_time[2] = LinSolvIter;
  
#ifndef NO_MPI
  
  /*--- Send/Receive information ---*/
  MPI::COMM_WORLD.Reduce(sbuf_force, rbuf_force, nVar_Force, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
  MPI::COMM_WORLD.Reduce(sbuf_time, rbuf_time, nVar_Time, MPI::UNSIGNED_LONG, MPI::MAX, MASTER_NODE);
  MPI::COMM_WORLD.Barrier();
  
  /*--- Set the value of the funcionals, and residuals ---*/
  if (rank == MASTER_NODE) {
    
    Total_CLift = rbuf_force[0];  Total_CDrag = rbuf_force[1];  Total_CSideForce = rbuf_force[2];
    Total_CMx = rbuf_force[3];    Total_CMy = rbuf_force[4];    Total_CMz = rbuf_force[5];
    Total_CFx = rbuf_force[6];    Total_CFy = rbuf_force[7];    Total_CFz = rbuf_force[8];
    Total_CEff = rbuf_force[0]/(rbuf_force[1]+EPS);
    
    if (freesurface) Total_CFreeSurface = rbuf_force[9];
    if (fluid_structure) Total_CFEA = rbuf_force[9];
    if (aeroacoustic) Total_CWave = rbuf_force[9];
    if (wave) Total_CWave = rbuf_force[9];
    if (fea) Total_CFEA = rbuf_force[9];
    if (isothermal) {
      Total_Q = rbuf_force[12];
      Total_MaxQ = rbuf_force[13];
    }
    
    /*--- Note that the nearfield based functionals should be redefined using the Cd weight ---*/
    if (equiv_area) {
      Total_CEquivArea  = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0-config[val_iZone]->GetWeightCd())*(rbuf_force[9]/double(size));
      Total_CNearFieldOF  = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0-config[val_iZone]->GetWeightCd())*rbuf_force[10];
    }
    
    if (rotating_frame) {
      Total_CMerit = rbuf_force[9];
      Total_CT = rbuf_force[10];
      Total_CQ = rbuf_force[11];
    }
    
    if (plasma) {
      Total_Q = rbuf_force[13];
      PressureDrag = rbuf_force[10];
      ViscDrag = rbuf_force[11];
      MagDrag = rbuf_force[12];
    }
    
    if (adjoint) {
      Total_Sens_Geo  = rbuf_force[0];
      Total_Sens_Mach = rbuf_force[1];
      Total_Sens_AoA  = rbuf_force[2];
      Total_Sens_Press = rbuf_force[3];
      Total_Sens_Temp  = rbuf_force[4];
    }
    
    timeused = rbuf_time[1];
    LinSolvIter = rbuf_time[2];
    
  }
  
#else
  
  /*--- Note that the nearfield based functionals should be redefined using the Cd weight ---*/
  if (equiv_area) {
    Total_CEquivArea  = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0-config[val_iZone]->GetWeightCd())*Total_CEquivArea;
    Total_CNearFieldOF  = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0-config[val_iZone]->GetWeightCd())*Total_CNearFieldOF;
  }
  
#endif
  
  /*--- Header frecuency ---*/
  bool Unsteady = ((config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                   (config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool In_NoDualTime = (!DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
  bool In_DualTime_0 = (DualTime_Iteration && (iIntIter % config[val_iZone]->GetWrt_Con_Freq_DualTime() == 0));
  bool In_DualTime_1 = (!DualTime_Iteration && Unsteady);
  bool In_DualTime_2 = (Unsteady && DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
  bool In_DualTime_3 = (Unsteady && !DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
  
  bool write_heads;
  if (Unsteady) write_heads = (iIntIter == 0);
  else write_heads = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*20)) == 0));
  
  if ((In_NoDualTime || In_DualTime_0 || In_DualTime_1) && (In_NoDualTime || In_DualTime_2 || In_DualTime_3)) {
    
    /*--- Output using all the processors ---*/
#ifndef NO_MPI
    MPI::COMM_WORLD.Barrier();
#endif
     
    /*--- Output using only the master node ---*/
    if (rank == MASTER_NODE) {
      
      /*--- Prepare the history file output, note that the dual
       time output don't write to the history file ---*/
      if (!DualTime_Iteration) {
        
        /*--- Write the begining of the history file ---*/
        sprintf (begin, "%12d", int(iExtIter));
        
        /*--- Write the end of the history file ---*/
        sprintf (end, ", %12.10f, %12.10f\n", double(LinSolvIter), double(timeused)/(CLOCKS_PER_SEC*60.0));
        
        /*--- Write the solution and residual of the history file ---*/
        switch (config[val_iZone]->GetKind_Solver()) {
            
          case EULER : case NAVIER_STOKES: case RANS:
          case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
          case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES: case FLUID_STRUCTURE_RANS:
          case AEROACOUSTIC_EULER: case AEROACOUSTIC_NAVIER_STOKES: case AEROACOUSTIC_RANS:
          case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
          case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
          case ADJ_AEROACOUSTIC_EULER:
            
            /*--- Direct coefficients ---*/
            sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                     Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                     Total_CFz, Total_CEff);
            if (isothermal)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy,
                       Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_Q, Total_MaxQ);
            if (equiv_area)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy,
                       Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CEquivArea, Total_CNearFieldOF);
            if (rotating_frame)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                       Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CMerit, Total_CT, Total_CQ);
            if (freesurface) {
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                       Total_CFz, Total_CEff, Total_CFreeSurface);
            }
            if (fluid_structure)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz,
                       Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CFEA);
            if (aeroacoustic)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz,
                       Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CWave);
            
            /*--- Flow residual ---*/
            if (nDim == 2) {
              if (incompressible) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), dummy, dummy );
              else sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), dummy );
            }
            else {
              if (incompressible) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), dummy );
              else sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), log10 (residual_flow[4]) );
            }
            
            /*--- Turbulent residual ---*/
            if (turbulent){
              switch(nVar_Turb) {
                case 1: sprintf (turb_resid, ", %12.10f", log10 (residual_turbulent[0])); break;
                case 2: sprintf (turb_resid, ", %12.10f, %12.10f", log10(residual_turbulent[0]), log10(residual_turbulent[1])); break;
              }
            }
            
            /*--- Transition residual ---*/
            if (transition){
                sprintf (trans_resid, ", %12.10f, %12.10f", log10(residual_transition[0]), log10(residual_transition[1]));
            }
            
            /*--- Free surface residual ---*/
            if (freesurface) {
              sprintf (levelset_resid, ", %12.10f", log10 (residual_levelset[0]));
            }
            
            /*--- Fluid structure residual ---*/
            if (fluid_structure) {
              if (nDim == 2) sprintf (levelset_resid, ", %12.10f, %12.10f, 0.0", log10 (residual_fea[0]), log10 (residual_fea[1]));
              else sprintf (levelset_resid, ", %12.10f, %12.10f, %12.10f", log10 (residual_fea[0]), log10 (residual_fea[1]), log10 (residual_fea[2]));
            }
            
            /*--- Aeroacoustics residual ---*/
            if (aeroacoustic) {
              sprintf (levelset_resid, ", %12.10f, %12.10f", log10 (residual_wave[0]), log10 (residual_wave[1]));
            }
            
            if (adjoint) {
              
              /*--- Adjoint coefficients ---*/
              sprintf (adjoint_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, 0.0", Total_Sens_Geo, Total_Sens_Mach, Total_Sens_AoA, Total_Sens_Press, Total_Sens_Temp);
              
              /*--- Adjoint flow residuals ---*/
              if (nDim == 2) {
                if (incompressible) sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, 0.0, 0.0", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]) );
                else sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]),log10 (residual_adjflow[3]) );
                
              }
              else {
                if (incompressible) sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]),log10 (residual_adjflow[3]) );
                else sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]),log10 (residual_adjflow[3]), log10 (residual_adjflow[4]) );
              }
              
              /*--- Adjoint turbulent residuals ---*/
              if (turbulent)
                  if (((config[val_iZone]->GetKind_Adjoint() != HYBRID) && (!config[val_iZone]->GetFrozen_Visc())) || (config[val_iZone]->GetKind_Adjoint() == HYBRID))
                      sprintf (adj_turb_resid, ", %12.10f", log10 (residual_adjturbulent[0]));
              
              /*--- Adjoint free surface residuals ---*/
              if (freesurface) sprintf (adj_levelset_resid, ", %12.10f", log10 (residual_adjlevelset[0]));
            }
            
            break;
            
          case PLASMA_EULER : case ADJ_PLASMA_EULER : case PLASMA_NAVIER_STOKES: case ADJ_PLASMA_NAVIER_STOKES:
            unsigned short iSpecies, loc;
            
            /*--- Direct problem coefficients ---*/
            sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                     Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                     Total_CFz, Total_CEff, Total_Q, PressureDrag, ViscDrag, MagDrag);
            //                        PressureDrag += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_PressureDrag();
            //                        ViscDrag += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_ViscDrag();
            //                        MagDrag += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_MagnetDrag();
            
            /*--- Direct problem residual ---*/
            for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
              if ( iSpecies < config[val_iZone]->GetnDiatomics() ) loc = (nDim+3)*iSpecies;
              else loc = (nDim+3)*config[val_iZone]->GetnDiatomics() + (nDim+2)*(iSpecies-config[val_iZone]->GetnDiatomics());
              sprintf (resid_aux, ", %12.10f, %12.10f", log10 (residual_plasma[loc+0]), log10 (residual_plasma[loc+nDim+1]));
              if (iSpecies == 0) strcpy(plasma_resid, resid_aux);
              else strcat(plasma_resid, resid_aux);
            }
            //                        sbuf_force[10] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_PressureDrag();
            //                        sbuf_force[11] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_ViscDrag();
            //                        sbuf_force[12] += solver_container[val_iZone][FinestMesh][PLASMA_SOL]->Get_MagnetDrag();
            
            /*--- Adjoint problem coefficients ---*/
            if (adjoint) {
              sprintf (adjoint_coeff, ", 0.0, 0.0, 0.0, 0.0");
              for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
                if ( iSpecies < config[val_iZone]->GetnDiatomics() ) loc = (nDim+3)*iSpecies;
                else loc = (nDim+3)*config[val_iZone]->GetnDiatomics() + (nDim+2)*(iSpecies-config[val_iZone]->GetnDiatomics());
                sprintf (resid_aux, ", %12.10f, %12.10f", log10 (residual_adjplasma[loc+0]),log10 (residual_adjplasma[loc+nDim+1]));
                if (iSpecies == 0) strcpy(adj_plasma_resid, resid_aux);
                else strcat(adj_plasma_resid, resid_aux);
              }
            }
            
            break;
            
          case WAVE_EQUATION:
            
            sprintf (direct_coeff, ", %12.10f", Total_CWave);
            sprintf (wave_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_wave[0]), log10 (residual_wave[1]), dummy, dummy, dummy );
            
            break;
            
          case LINEAR_ELASTICITY:
            
            sprintf (direct_coeff, ", %12.10f", Total_CFEA);
            sprintf (fea_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_fea[0]), dummy, dummy, dummy, dummy );
            
            break;
            
        }
      }
      
      /*--- Write the screen header---*/
      if ((write_heads) && !(!DualTime_Iteration && Unsteady)) {
        
        if (!Unsteady) {
          switch (config[val_iZone]->GetKind_Solver()) {
            case EULER :                  case NAVIER_STOKES:
            case FLUID_STRUCTURE_EULER :  case FLUID_STRUCTURE_NAVIER_STOKES:
            case AEROACOUSTIC_EULER :     case AEROACOUSTIC_NAVIER_STOKES:
              cout << endl << " Min Delta Time: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<<
              ". Max Delta Time: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() << ".";
              break;
              
            case PLASMA_EULER: case PLASMA_NAVIER_STOKES:
            case ADJ_PLASMA_EULER: case ADJ_PLASMA_NAVIER_STOKES:
              
              for (unsigned short iSpecies = 0; iSpecies < config[val_iZone]->GetnSpecies(); iSpecies++) {
                cout << endl << " Min Delta Time (" << iSpecies << "): " << solver_container[val_iZone][MESH_0][PLASMA_SOL]->GetMin_Delta_Time(iSpecies)<< ". Max Delta Time (" << iSpecies << "): " << solver_container[val_iZone][MESH_0][PLASMA_SOL]->GetMax_Delta_Time(iSpecies) << ".";
              }
              break;
          }
        }
        else {
          cout << endl << " Min DT: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<<
          ". Max DT: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() <<
          ". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
        }
                
        switch (config[val_iZone]->GetKind_Solver()) {
          case EULER :                  case NAVIER_STOKES:
          case FLUID_STRUCTURE_EULER :  case FLUID_STRUCTURE_NAVIER_STOKES:
          case AEROACOUSTIC_EULER :     case AEROACOUSTIC_NAVIER_STOKES:
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max(0) << "." << endl;

            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
                        
            if (!fluid_structure && !aeroacoustic) {
              if (incompressible) cout << "   Res[Press]" << "     Res[Velx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
              else if (rotating_frame && nDim == 3) cout << "     Res[Rho]" << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
              else if (equiv_area) cout << "     Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << "    CPress(N-F)" << endl;
              else cout << "     Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
            }
            else if (fluid_structure) cout << "     Res[Rho]" << "   Res[Displx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
            else if (aeroacoustic) cout << "     Res[Rho]" << "   Res[Wave]" << "   CLift(Total)" << "   CDrag(Total)" << "   CWave(Total)" << endl;
            break;
            
          case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max(0) << "." << endl;

            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            
            cout << "   Res[Press]" << "     Res[Dist]" << "   CLift(Total)" << "     CLevelSet" << endl;
            
            break;
            
          case RANS : case FLUID_STRUCTURE_RANS: case AEROACOUSTIC_RANS:
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max(0) << "." << endl;

            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            if (incompressible) cout << "   Res[Press]";
            else cout << "      Res[Rho]";

            switch (config[val_iZone]->GetKind_Turb_Model()){
              case SA:	cout << "       Res[nu]"; break;
              case SST:	cout << "     Res[kine]" << "     Res[omega]"; break;
            }
            
            if (transition) { cout << "      Res[Int]" << "       Res[Re]"; }
            if (rotating_frame && nDim == 3 ) cout << "   CThrust(Total)" << "   CTorque(Total)" << endl;
            else cout << "   CLift(Total)"   << "   CDrag(Total)"   << endl;
            break;
            
          case PLASMA_EULER : case PLASMA_NAVIER_STOKES :
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][PLASMA_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            if (config[val_iZone]->GetKind_GasModel() == ARGON)
              cout << "      Res[r1]" << "       Res[r2]" << "   Res(r3)"<<  endl;
            if (config[val_iZone]->GetKind_GasModel() == AIR21)
              cout << "      Res[r1]" << "       Res[r2]" << "   Res(r3)"<<  endl;
            if ((config[val_iZone]->GetKind_GasModel() == ARGON_SID) || (config[val_iZone]->GetKind_GasModel() == AIR7) || (config[val_iZone]->GetKind_GasModel() == AIR5) || (config[val_iZone]->GetKind_GasModel() == N2) || (config[val_iZone]->GetKind_GasModel() == O2))
              cout << "      Res[Rho0]" << "      Res[E0]" << "	  Q(Total)" << "	 CDrag(Total)" << endl;
            break;
            
          case WAVE_EQUATION :
            cout << "      Res[Wave]" << "   CWave(Total)"<<  endl;
            break;
            
          case LINEAR_ELASTICITY :
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            if (nDim == 2) cout << "    Res[Displx]" << "    Res[Disply]" << "   CFEA(Total)"<<  endl;
            if (nDim == 3) cout << "    Res[Displx]" << "    Res[Disply]" << "    Res[Displz]" << "   CFEA(Total)"<<  endl;
            break;
            
          case ADJ_PLASMA_EULER : case ADJ_PLASMA_NAVIER_STOKES :
            if ((config[val_iZone]->GetKind_GasModel() == ARGON_SID) || (config[val_iZone]->GetKind_GasModel() == AIR7) || (config[val_iZone]->GetKind_GasModel() == AIR5) || (config[val_iZone]->GetKind_GasModel() == N2) || (config[val_iZone]->GetKind_GasModel() == O2))
              
            /*--- Visualize the maximum residual ---*/
              cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetRes_Max(0))
              <<", located at point "<< solver_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetPoint_Max(0) << "." << endl;
            
              if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
              else cout << endl << " IntIter" << "  ExtIter";
            
              cout << "        Res[Psi_Rho0]" << "       Res[Psi_E0]" << "	      Sens_Geo" << endl;
            break;
            
          case ADJ_EULER :              case ADJ_NAVIER_STOKES :
          case ADJ_AEROACOUSTIC_EULER :
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            if (incompressible) cout << "   Res[Psi_Press]" << "   Res[Psi_Velx]";
            else cout << "   Res[Psi_Rho]" << "     Res[Psi_E]";
            cout << "     Sens_Geo" << "    Sens_Mach" << endl;
            break;
            
          case ADJ_FREE_SURFACE_EULER : case ADJ_FREE_SURFACE_NAVIER_STOKES :
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            
            cout << "   Res[Psi_Press]" << "   Res[Psi_Dist]" << "    Sens_Geo" << "   Sens_Mach" << endl;
            break;
            
          case ADJ_RANS :
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            if (incompressible) cout << "     Res[Psi_Press]";
            else cout << "     Res[Psi_Rho]";
                
            if (((config[val_iZone]->GetKind_Adjoint() != HYBRID) && (!config[val_iZone]->GetFrozen_Visc())) || (config[val_iZone]->GetKind_Adjoint() == HYBRID)) {
              cout << "      Res[Psi_nu]";
            }
            else {
              if (incompressible) cout << "   Res[Psi_Velx]";
              else cout << "     Res[Psi_E]";
            }
            cout << "     Sens_Geo" << "    Sens_Mach" << endl;
            break;
            
            
          case ADJ_FREE_SURFACE_RANS :
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            
            cout << "   Res[Psi_Press]" << "   Res[Psi_Dist]" << "    Sens_Geo" << "   Sens_Mach" << endl;
            break;
            
        }
        
      }
      
      /*--- Write the solution on the screen and history file ---*/
      cout.precision(6);
      cout.setf(ios::fixed,ios::floatfield);
      
      if (!Unsteady) {
        cout.width(5); cout << iExtIter;
        cout.width(11); cout << double(timeiter)/CLOCKS_PER_SEC;
        
      } else {
        cout.width(8); cout << iIntIter;
        cout.width(8); cout << iExtIter;
      }
      
      switch (config[val_iZone]->GetKind_Solver()) {
        case EULER : case NAVIER_STOKES:
        case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES:
        case AEROACOUSTIC_EULER: case AEROACOUSTIC_NAVIER_STOKES:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid;
            if (fluid_structure) ConvHist_file[0] << fea_resid;
            if (aeroacoustic) ConvHist_file[0] << levelset_resid;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(13); cout << log10(residual_flow[0]);
          if (!fluid_structure && !aeroacoustic && !equiv_area) {
            if (!incompressible ) {
              if (nDim == 2 ) { cout.width(14); cout << log10(residual_flow[3]); }
              else { cout.width(14); cout << log10(residual_flow[4]); }
            }
            else { cout.width(14); cout << log10(residual_flow[1]); }
          }
          else if (fluid_structure) { cout.width(14); cout << log10(residual_fea[0]); }
          else if (aeroacoustic) { cout.width(14); cout << log10(residual_wave[0]); }
          
          if (rotating_frame && nDim == 3 ) { cout.width(15); cout << Total_CT; cout.width(15); cout << Total_CQ; }
          else if (aeroacoustic) { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CDrag; cout.width(15); cout << Total_CWave; }
          else if (equiv_area) { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CDrag; cout.width(15);
            cout.precision(4);
            cout.setf(ios::scientific,ios::floatfield);
            cout << Total_CNearFieldOF; }
          else { cout.width(15); cout << min(1000.0,max(-1000.0, Total_CLift)); cout.width(15); cout << min(1000.0,max(-1000.0, Total_CDrag)); }
          cout << endl;
          break;
          
        case RANS :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid << turb_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          
          if (incompressible) cout.width(13);
          else  cout.width(14);
          cout << log10(residual_flow[0]);
          
          switch(nVar_Turb) {
            case 1: cout.width(14); cout << log10(residual_turbulent[0]); break;
            case 2: cout.width(14); cout << log10(residual_turbulent[0]);
              cout.width(14); cout << log10(residual_turbulent[1]); break;
          }
          
          if (transition) { cout.width(14); cout << log10(residual_transition[0]); cout.width(14); cout << log10(residual_transition[1]); }
      
          if (rotating_frame  && nDim == 3 ) { cout.width(15); cout << Total_CT; cout.width(15); cout << Total_CQ; }
          else { cout.width(15); cout << min(1000.0,max(-1000.0, Total_CLift)); cout.width(15); cout << min(1000.0,max(-1000.0, Total_CDrag)); }
          cout << endl;
          break;
          
        case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid << levelset_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(13); cout << log10(residual_flow[0]);
          cout.width(14); cout << log10(residual_levelset[0]);
          cout.width(15); cout << Total_CLift;
          cout.width(14); cout << Total_CFreeSurface;
          
          cout << endl;
          break;
          
        case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << plasma_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          if (config[val_iZone]->GetKind_GasModel() == ARGON || config[val_iZone]->GetKind_GasModel() == AIR21) {
            cout.width(14); cout << log10(residual_plasma[0]);
            cout.width(14); cout << log10(residual_plasma[nDim+2]);
            cout.width(14); cout << log10(residual_plasma[2*(nDim+2)]);
          }
          if ((config[val_iZone]->GetKind_GasModel() == ARGON_SID) || config[val_iZone]->GetKind_GasModel() == AIR7 || config[val_iZone]->GetKind_GasModel() == O2 || config[val_iZone]->GetKind_GasModel() == N2 || config[val_iZone]->GetKind_GasModel() == AIR5) {
            cout.width(14); cout << log10(residual_plasma[0]);
            cout.width(14); cout << log10(residual_plasma[nDim+1]);
            cout.width(14); cout << Total_Q;
            cout.width(14); cout << Total_CDrag;
          }
          cout << endl;
          break;
          
        case WAVE_EQUATION:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << wave_coeff << wave_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(14); cout << log10(residual_wave[0]);
          cout.width(14); cout << Total_CWave;
          cout << endl;
          break;
          
        case LINEAR_ELASTICITY:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << fea_coeff << fea_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(15); cout << log10(residual_fea[0]);
          cout.width(15); cout << log10(residual_fea[1]);
          if (nDim == 3) { cout.width(15); cout << log10(residual_fea[2]); }
          cout.width(14); cout << Total_CFEA;
          cout << endl;
          break;
          
        case ADJ_EULER :              case ADJ_NAVIER_STOKES :
        case ADJ_AEROACOUSTIC_EULER :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          if (!incompressible) {
            cout.width(15); cout << log10(residual_adjflow[0]);
            cout.width(15); cout << log10(residual_adjflow[nDim+1]);
          }
          else {
            cout.width(17); cout << log10(residual_adjflow[0]);
            cout.width(16); cout << log10(residual_adjflow[1]);
          }
          cout.precision(4);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(14); cout << Total_Sens_Geo;
          cout.width(14); cout << Total_Sens_Mach;
          cout << endl;
          break;
          
        case ADJ_FREE_SURFACE_EULER : case ADJ_FREE_SURFACE_NAVIER_STOKES :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << adj_levelset_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(17); cout << log10(residual_adjflow[0]);
          cout.width(16); cout << log10(residual_adjlevelset[0]);
          cout.precision(3);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(12); cout << Total_Sens_Geo;
          cout.width(12); cout << Total_Sens_Mach;
          cout << endl;
          break;
                    
        case ADJ_RANS :
          
          if (!DualTime_Iteration) {
              ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid;
              if (((config[val_iZone]->GetKind_Adjoint() != HYBRID) && (!config[val_iZone]->GetFrozen_Visc())) || (config[val_iZone]->GetKind_Adjoint() == HYBRID))
                  ConvHist_file[0] << adj_turb_resid;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(17); cout << log10(residual_adjflow[0]);
          if (((config[val_iZone]->GetKind_Adjoint() != HYBRID) && (!config[val_iZone]->GetFrozen_Visc())) || (config[val_iZone]->GetKind_Adjoint() == HYBRID)) {
            cout.width(17); cout << log10(residual_adjturbulent[0]);
          }
          else {
            if (!incompressible) {
              if (geometry[val_iZone][FinestMesh]->GetnDim() == 2 ) { cout.width(15); cout << log10(residual_adjflow[3]); }
              else { cout.width(15); cout << log10(residual_adjflow[4]); }
            }
            else {
              cout.width(15); cout << log10(residual_adjflow[1]);
            }
          }
          cout.precision(4);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(14); cout << Total_Sens_Geo;
          cout.width(14); cout << Total_Sens_Mach;
          cout << endl;
          break;
          
        case ADJ_FREE_SURFACE_RANS :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << adj_levelset_resid;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(17); cout << log10(residual_adjflow[0]);
          cout.width(16); cout << log10(residual_adjlevelset[0]);
          
          cout.precision(4);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(12); cout << Total_Sens_Geo;
          cout.width(12); cout << Total_Sens_Mach;
          cout << endl;
          break;
          
        case ADJ_PLASMA_EULER : case ADJ_PLASMA_NAVIER_STOKES:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_plasma_resid << end;
            ConvHist_file[0].flush();
          }
          
          if ((config[val_iZone]->GetKind_GasModel() == ARGON_SID) || config[val_iZone]->GetKind_GasModel() == AIR7 || config[val_iZone]->GetKind_GasModel() == O2 || config[val_iZone]->GetKind_GasModel() == N2 || config[val_iZone]->GetKind_GasModel() == AIR5) {
            cout.width(19); cout << log10(residual_adjplasma[0]);
            cout.width(19); cout << log10(residual_adjplasma[nDim+1]);
            cout.width(19); cout << Total_Sens_Geo;
          }
          cout << endl;
          break;
          
      }
      cout.unsetf(ios::fixed);
    }
    
    delete [] residual_flow;
    delete [] residual_levelset;
    delete [] residual_plasma;
    delete [] residual_turbulent;
    delete [] residual_transition;
    delete [] residual_wave;
    delete [] residual_fea;
    
    delete [] residual_adjflow;
    delete [] residual_adjlevelset;
    delete [] residual_adjplasma;
    delete [] residual_adjturbulent;
    
    delete [] sbuf_force;
    delete [] sbuf_time;
    
    if (rank == MASTER_NODE) {
      delete [] rbuf_force;
      delete [] rbuf_time;
    }
    
  }
}

void COutput::SetResult_Files(CSolver ****solver_container, CGeometry ***geometry, CConfig **config,
		unsigned long iExtIter, unsigned short val_nZone) {

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	unsigned short iZone;

	/*--- Allocate space for the local data arrays ---*/
	//nOutput_Vars = new unsigned short*[MAX_ZONES];
	//data_container = new double***[MAX_ZONES];

	for (iZone = 0; iZone < val_nZone; iZone++) {

		//nOutput_Vars[iZone] = new unsigned short[MAX_SOLS];
		//data_container[iZone] = new double**[MAX_SOLS];

		/*--- Flags identifying the types of files to be written. ---*/
		bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
    bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
    
#ifndef NO_MPI
    /*--- Do not merge the volume solutions if we are running in parallel.
     Force the use of SU2_SOL to merge the volume sols in this case. ---*/
    int size = MPI::COMM_WORLD.Get_size();
    if (size > 1) {
      Wrt_Vol = false;
      Wrt_Srf = false;
    }
#endif
    
		bool Wrt_Csv = config[iZone]->GetWrt_Csv_Sol();
		bool Wrt_Rst = config[iZone]->GetWrt_Restart();

		switch (config[iZone]->GetKind_Solver()) {

		case EULER : case NAVIER_STOKES : case RANS :
		case FREE_SURFACE_EULER : case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
		case FLUID_STRUCTURE_EULER : case FLUID_STRUCTURE_NAVIER_STOKES : case FLUID_STRUCTURE_RANS:

			if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][FLOW_SOL], iExtIter);
			break;

		case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS : case ADJ_FREE_SURFACE_EULER : case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
			if (Wrt_Csv) SetSurfaceCSV_Adjoint(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][ADJFLOW_SOL], solver_container[iZone][MESH_0][FLOW_SOL], iExtIter, iZone);
			break;

		case LIN_EULER : case LIN_NAVIER_STOKES :
			if (Wrt_Csv) SetSurfaceCSV_Linearized(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][LINFLOW_SOL], config[iZone]->GetSurfLinCoeff_FileName(), iExtIter);
			break;

		case AEROACOUSTIC_EULER : case AEROACOUSTIC_NAVIER_STOKES : case AEROACOUSTIC_RANS:
			if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][FLOW_SOL], iExtIter);
			break;
		case ADJ_AEROACOUSTIC_EULER:
			if (iZone == ZONE_0) {
				if (Wrt_Csv) SetSurfaceCSV_Adjoint(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][ADJFLOW_SOL], solver_container[iZone][MESH_0][FLOW_SOL], iExtIter, iZone);
			} else if (iZone == ZONE_1) {
				if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][FLOW_SOL], iExtIter);
			}
			break;
		}

		/*--- Get the file output format ---*/

		unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();

		bool dynamic_mesh = (config[iZone]->GetUnsteady_Simulation() &&
				config[iZone]->GetGrid_Movement());

		/*--- Merge the node coordinates and connectivity if necessary. This
         is only performed if a volume solution file is requested, and it
         is active by default. ---*/

    if (Wrt_Vol || Wrt_Srf)
			MergeGeometry(config[iZone], geometry[iZone][MESH_0], iZone);
  
		/*--- Merge the solution data needed for volume solutions and restarts ---*/

		if (Wrt_Vol || Wrt_Rst)
			MergeSolution(config[iZone], geometry[iZone][MESH_0],
					solver_container[iZone][MESH_0], iZone);
  
		/*--- Write restart, CGNS, or Tecplot files using the merged data.
         This data lives only on the master, and these routines are currently
         executed by the master proc alone (as if in serial). ---*/

		if (rank == MASTER_NODE) {
 
			/*--- Write a native restart file ---*/
			if (Wrt_Rst)
				SetRestart(config[iZone], geometry[iZone][MESH_0], iZone);
  
			if (Wrt_Vol) {

				switch (FileFormat) {
            
				case TECPLOT:

					/*--- Write a Tecplot ASCII file ---*/
					SetTecplot_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, false);
					DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
					break;

				case TECPLOT_BINARY:

					/*--- Write a Tecplot binary solution file ---*/
					SetTecplot_Solution(config[iZone], geometry[iZone][MESH_0], iZone);
					if (dynamic_mesh) DeallocateCoordinates(config[iZone], geometry[iZone][MESH_0]);
					break;

				case CGNS_SOL:

					/*--- Write a CGNS solution file ---*/
					SetCGNS_Solution(config[iZone], geometry[iZone][MESH_0], iZone);
					break;

				default:
					break;
				}
        
			}
      
      if (Wrt_Srf) {
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
            break;
            
          default:
            break;
        }
        
      }

			/*--- Release memory needed for merging the solution data. ---*/
      if (((Wrt_Vol) || (Wrt_Srf)) && (FileFormat == TECPLOT))
        DeallocateCoordinates(config[iZone], geometry[iZone][MESH_0]);
      
      if (Wrt_Vol || Wrt_Rst)
        DeallocateSolution(config[iZone], geometry[iZone][MESH_0]);
      
		}

		/*--- Final broadcast (informing other procs that the base output
         file was written) & barrier to sync up after master node writes
         output files. ---*/

#ifndef NO_MPI
		MPI::COMM_WORLD.Bcast(&wrote_base_file, 1, MPI::INT, MASTER_NODE);
		MPI::COMM_WORLD.Barrier();
#endif

	}
}

void COutput::SetBaselineResult_Files(CSolver **solver, CGeometry **geometry, CConfig **config,
                                      unsigned long iExtIter, unsigned short val_nZone) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
	unsigned short iZone;
  
	for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Flags identifying the types of files to be written. ---*/
    
		bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
		bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
    bool Wrt_Rst = config[iZone]->GetWrt_Restart();
    
    /*--- Get the file output format ---*/
    
    unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();
    
    bool dynamic_mesh = (config[iZone]->GetUnsteady_Simulation() &&
                         config[iZone]->GetGrid_Movement());
    
    /*--- Merge the node coordinates and connectivity if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/
    
    if (Wrt_Vol || Wrt_Srf) {
      if (rank == MASTER_NODE) cout <<"Merging geometry." << endl;
      MergeGeometry(config[iZone], geometry[iZone], iZone);
    }
    
    /*--- Merge the solution data needed for volume solutions and restarts ---*/
    
    if (Wrt_Vol || Wrt_Rst) {
      if (rank == MASTER_NODE) cout <<"Merging solution." << endl;
      MergeBaselineSolution(config[iZone], geometry[iZone], solver[iZone], iZone);
    }
    
    /*--- Write restart, CGNS, or Tecplot files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
    
    if (rank == MASTER_NODE) {
      
      if (Wrt_Vol) {
        
        if (rank == MASTER_NODE) cout <<"Writing volume solution." << endl;

        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config[iZone], geometry[iZone], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone], false);
            break;
            
          case TECPLOT_BINARY:
            
            /*--- Write a Tecplot binary solution file ---*/
            SetTecplot_Solution(config[iZone], geometry[iZone], iZone);
            if (dynamic_mesh) DeallocateCoordinates(config[iZone], geometry[iZone]);
            break;
            
          case CGNS_SOL:
            
            /*--- Write a CGNS solution file ---*/
            SetCGNS_Solution(config[iZone], geometry[iZone], iZone);
            break;
            
          default:
            break;
        }
        
      }
      
      if (Wrt_Srf) {
        
        if (rank == MASTER_NODE) cout <<"Writing surface solution." << endl;
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config[iZone], geometry[iZone], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone], true);
            break;
            
          default:
            break;
        }
        
      }
      
			/*--- Release memory needed for merging the solution data. ---*/
      if (((Wrt_Vol) || (Wrt_Srf)) && (FileFormat == TECPLOT))
        DeallocateCoordinates(config[iZone], geometry[iZone]);
      
      if (Wrt_Vol || Wrt_Rst)
        DeallocateSolution(config[iZone], geometry[iZone]);
    }
    
    /*--- Final broadcast (informing other procs that the base output
     file was written) & barrier to sync up after master node writes
     output files. ---*/
    
#ifndef NO_MPI
    MPI::COMM_WORLD.Bcast(&wrote_base_file, 1, MPI::INT, MASTER_NODE);
    MPI::COMM_WORLD.Barrier();
#endif
    
  }
}

void COutput::SetFlowRate(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {

	ifstream index_file;
	int integration_node[100];
	int nPointFlowRate = 0;

	/*--- Read list of nodes ---*/
	index_file.open("flowrate_nodes.dat", ios::in);
	while (index_file.good()) {
		index_file >> integration_node[nPointFlowRate];
		nPointFlowRate++;
	}
	nPointFlowRate--;
	index_file.close();


	/*--- Perform trapezoid integration ---*/
	double y1, y2, q1, q2;
	double integral = 0.0;
	for (int j=0; j<nPointFlowRate-1; j++) {
		y1 = geometry->node[integration_node[j]]->GetCoord(1);
		y2 = geometry->node[integration_node[j+1]]->GetCoord(1);
		q1 = solver_container->node[integration_node[j]]->GetVelocity(0, config->GetIncompressible());
		q2 = solver_container->node[integration_node[j+1]]->GetVelocity(0, config->GetIncompressible());

		integral = integral + 0.5*(q1+q2)*(y2-y1);
	}

	/*--- Store integral in solver_container ---*/
	solver_container->SetTotal_CEquivArea(integral);  // integral shows up as EquivArea in history
}

void COutput::SetEquivalentArea(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {

	ofstream EquivArea_file, FuncGrad_file;
	unsigned short iMarker = 0, iDim;
	short *AzimuthalAngle = NULL;
	double Gamma, auxXCoord, auxYCoord, auxZCoord, InverseDesign, DeltaX, Coord_i, Coord_j, jp1Coord, *Coord = NULL, MeanFuntion,
			*Face_Normal = NULL, auxArea, auxPress, Mach, Beta, R_Plane, Pressure_Inf, Density_Inf,
			RefAreaCoeff, ModVelocity_Inf, Velocity_Inf[3], factor, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL,
			*Pressure = NULL, *FaceArea = NULL, *EquivArea = NULL, *TargetArea = NULL, *NearFieldWeight = NULL,
			*Weight = NULL, jFunction, jp1Function;
	unsigned long jVertex, iVertex, iPoint, nVertex_NearField = 0, auxPoint,
			*IdPoint = NULL, *IdDomain = NULL, auxDomain;
	unsigned short iPhiAngle;
	ofstream NearFieldEA_file; ifstream TargetEA_file;

	bool incompressible = config->GetIncompressible();
	double XCoordBegin_OF = config->GetEA_IntLimit(0);
	double XCoordEnd_OF = config->GetEA_IntLimit(1);

	unsigned short nDim = geometry->GetnDim();
	double AoA = -(config->GetAoA()*PI_NUMBER/180.0);

	int rank = MESH_0;

	Mach  = config->GetMach_FreeStreamND();
	Gamma = config->GetGamma();
	Beta = sqrt(Mach*Mach-1.0);
	R_Plane = fabs(config->GetEA_IntLimit(2));
	Pressure_Inf = config->GetPressure_FreeStreamND();
	Density_Inf = config->GetDensity_FreeStreamND();
	RefAreaCoeff = config->GetRefAreaCoeff();
	Velocity_Inf[0] = config->GetVelocity_FreeStreamND()[0];
	Velocity_Inf[1] = config->GetVelocity_FreeStreamND()[1];
	Velocity_Inf[2] = config->GetVelocity_FreeStreamND()[2];
	ModVelocity_Inf = 0;
	for (iDim = 0; iDim < 3; iDim++)
		ModVelocity_Inf += Velocity_Inf[iDim] * Velocity_Inf[iDim];
	ModVelocity_Inf = sqrt(ModVelocity_Inf);

	factor = 4.0*sqrt(2.0*Beta*R_Plane) / (Gamma*Pressure_Inf*Mach*Mach);

#ifdef NO_MPI

/*--- Compute the total number of points on the near-field ---*/
	nVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();
				/*--- Using Face_Normal(z), and Coord(z) we identify only a surface,
				 note that there are 2 NEARFIELD_BOUNDARY surfaces ---*/
				if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) nVertex_NearField ++;
			}

	/*--- Create an array with all the coordinates, points, pressures, face area,
	 equivalent area, and nearfield weight ---*/
	Xcoord = new double[nVertex_NearField];
	Ycoord = new double[nVertex_NearField];
	Zcoord = new double[nVertex_NearField];
	AzimuthalAngle = new short[nVertex_NearField];
	IdPoint = new unsigned long[nVertex_NearField];
	IdDomain = new unsigned long[nVertex_NearField];
	Pressure = new double[nVertex_NearField];
	FaceArea = new double[nVertex_NearField];
	EquivArea = new double[nVertex_NearField];
	TargetArea = new double[nVertex_NearField];
	NearFieldWeight = new double[nVertex_NearField];
	Weight = new double[nVertex_NearField];

	/*--- Copy the boundary information to an array ---*/
	nVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) {

					IdPoint[nVertex_NearField] = iPoint;
					Xcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
					Ycoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(1);

					if (nDim ==2) {
						AzimuthalAngle[nVertex_NearField] = 0;
					}

					if (nDim == 3) {
						Zcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(2);

						/*--- Rotate the nearfield cylinder (AoA) only 3D ---*/
								double YcoordRot = Ycoord[nVertex_NearField];
								double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA);

								/* Compute the Azimuthal angle (resolution of degress in the Azimuthal angle)---*/
								double AngleDouble; short AngleInt;
								AngleDouble = atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER;
								AngleInt = (short) floor(AngleDouble + 0.5);
								if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
								else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
					}

					if (AzimuthalAngle[nVertex_NearField] <= 60) {
						Pressure[nVertex_NearField] = solver_container->node[iPoint]->GetPressure(incompressible);
						FaceArea[nVertex_NearField] = fabs(Face_Normal[nDim-1]);
						nVertex_NearField ++;
					}

				}
			}

#else

	int nProcessor = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();

	unsigned long nLocalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
	int iProcessor;

	unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
	unsigned long *Buffer_Send_nVertex = new unsigned long [1];

	/*--- Compute the total number of points of the near-field ghost nodes ---*/
	nLocalVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				if (geometry->node[iPoint]->GetDomain())
					if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0))
						nLocalVertex_NearField ++;
			}

	Buffer_Send_nVertex[0] = nLocalVertex_NearField;

	/*--- Send Near-Field vertex information --*/
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &nVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Allgather(Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG);

	double *Buffer_Send_Xcoord = new double[MaxLocalVertex_NearField];
	double *Buffer_Send_Ycoord = new double[MaxLocalVertex_NearField];
	double *Buffer_Send_Zcoord = new double[MaxLocalVertex_NearField];
	unsigned long *Buffer_Send_IdPoint = new unsigned long [MaxLocalVertex_NearField];
	double *Buffer_Send_Pressure = new double [MaxLocalVertex_NearField];
	double *Buffer_Send_FaceArea = new double[MaxLocalVertex_NearField];

	double *Buffer_Receive_Xcoord = new double[nProcessor*MaxLocalVertex_NearField];
	double *Buffer_Receive_Ycoord = new double[nProcessor*MaxLocalVertex_NearField];
	double *Buffer_Receive_Zcoord = new double[nProcessor*MaxLocalVertex_NearField];
	unsigned long *Buffer_Receive_IdPoint = new unsigned long[nProcessor*MaxLocalVertex_NearField];
	double *Buffer_Receive_Pressure = new double[nProcessor*MaxLocalVertex_NearField];
	double *Buffer_Receive_FaceArea = new double[nProcessor*MaxLocalVertex_NearField];

	unsigned long nBuffer_Xcoord = MaxLocalVertex_NearField;
	unsigned long nBuffer_Ycoord = MaxLocalVertex_NearField;
	unsigned long nBuffer_Zcoord = MaxLocalVertex_NearField;
	unsigned long nBuffer_IdPoint = MaxLocalVertex_NearField;
	unsigned long nBuffer_Pressure = MaxLocalVertex_NearField;
	unsigned long nBuffer_FaceArea = MaxLocalVertex_NearField;

	for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
		Buffer_Send_IdPoint[iVertex] = 0; Buffer_Send_Pressure[iVertex] = 0.0;
		Buffer_Send_FaceArea[iVertex] = 0.0; Buffer_Send_Xcoord[iVertex] = 0.0;
		Buffer_Send_Ycoord[iVertex] = 0.0; Buffer_Send_Zcoord[iVertex] = 0.0;
	}

	/*--- Copy coordinates, index points, and pressures to the auxiliar vector --*/
	nLocalVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				if (geometry->node[iPoint]->GetDomain())
					if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) {
						Buffer_Send_IdPoint[nLocalVertex_NearField] = iPoint;
						Buffer_Send_Xcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
						Buffer_Send_Ycoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(1);
						Buffer_Send_Zcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(2);
						Buffer_Send_Pressure[nLocalVertex_NearField] = solver_container->node[iPoint]->GetPressure(incompressible);
						Buffer_Send_FaceArea[nLocalVertex_NearField] = fabs(Face_Normal[nDim-1]);
						nLocalVertex_NearField++;
					}
			}

	/*--- Send all the information --*/
	MPI::COMM_WORLD.Gather(Buffer_Send_Xcoord, nBuffer_Xcoord, MPI::DOUBLE, Buffer_Receive_Xcoord, nBuffer_Xcoord, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Ycoord, nBuffer_Ycoord, MPI::DOUBLE, Buffer_Receive_Ycoord, nBuffer_Ycoord, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Zcoord, nBuffer_Zcoord, MPI::DOUBLE, Buffer_Receive_Zcoord, nBuffer_Zcoord, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_IdPoint, nBuffer_IdPoint, MPI::UNSIGNED_LONG, Buffer_Receive_IdPoint, nBuffer_IdPoint, MPI::UNSIGNED_LONG, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Pressure, nBuffer_Pressure, MPI::DOUBLE, Buffer_Receive_Pressure, nBuffer_Pressure, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_FaceArea, nBuffer_FaceArea, MPI::DOUBLE, Buffer_Receive_FaceArea, nBuffer_FaceArea, MPI::DOUBLE, MASTER_NODE);

	if (rank == MASTER_NODE) {

		Xcoord = new double[nVertex_NearField];
		Ycoord = new double[nVertex_NearField];
		Zcoord = new double[nVertex_NearField];
		AzimuthalAngle = new short[nVertex_NearField];
		IdPoint = new unsigned long[nVertex_NearField];
		IdDomain = new unsigned long[nVertex_NearField];
		Pressure = new double[nVertex_NearField];
		FaceArea = new double[nVertex_NearField];
		EquivArea = new double[nVertex_NearField];
		TargetArea = new double[nVertex_NearField];
		NearFieldWeight = new double[nVertex_NearField];
		Weight = new double[nVertex_NearField];

		nVertex_NearField = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
			for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
				Xcoord[nVertex_NearField] = Buffer_Receive_Xcoord[iProcessor*MaxLocalVertex_NearField+iVertex];
				Ycoord[nVertex_NearField] = Buffer_Receive_Ycoord[iProcessor*MaxLocalVertex_NearField+iVertex];

				if (nDim == 2) {
					AzimuthalAngle[nVertex_NearField] = 0;
				}

				if (nDim == 3) {
					Zcoord[nVertex_NearField] = Buffer_Receive_Zcoord[iProcessor*MaxLocalVertex_NearField+iVertex];

					/*--- Rotate the nearfield cylinder  ---*/
					double YcoordRot = Ycoord[nVertex_NearField];
					double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA);

					/*--- Compute the Azimuthal angle ---*/
					double AngleDouble; short AngleInt;
					AngleDouble = atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER;
					AngleInt = (short) floor(AngleDouble + 0.5);
					if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
					else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
				}

				if (AzimuthalAngle[nVertex_NearField] <= 60) {
					IdPoint[nVertex_NearField] = Buffer_Receive_IdPoint[iProcessor*MaxLocalVertex_NearField+iVertex];
					Pressure[nVertex_NearField] = Buffer_Receive_Pressure[iProcessor*MaxLocalVertex_NearField+iVertex];
					FaceArea[nVertex_NearField] = Buffer_Receive_FaceArea[iProcessor*MaxLocalVertex_NearField+iVertex];
					IdDomain[nVertex_NearField] = iProcessor;
					nVertex_NearField++;
				}

			}
	}

	delete [] Buffer_Receive_nVertex;
	delete [] Buffer_Send_nVertex;

	delete [] Buffer_Send_Xcoord;
	delete [] Buffer_Send_Ycoord;
	delete [] Buffer_Send_Zcoord;
	delete [] Buffer_Send_IdPoint;
	delete [] Buffer_Send_Pressure;
	delete [] Buffer_Send_FaceArea;

	delete [] Buffer_Receive_Xcoord;
	delete [] Buffer_Receive_IdPoint;
	delete [] Buffer_Receive_Pressure;
	delete [] Buffer_Receive_FaceArea;

#endif

	if (rank == MASTER_NODE) {

		vector<short> PhiAngleList;
		vector<short>::iterator IterPhiAngleList;

		for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
			PhiAngleList.push_back(AzimuthalAngle[iVertex]);

		sort( PhiAngleList.begin(), PhiAngleList.end());
		IterPhiAngleList = unique( PhiAngleList.begin(), PhiAngleList.end());
		PhiAngleList.resize( IterPhiAngleList - PhiAngleList.begin() );

		/*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
		vector<vector<double> > Xcoord_PhiAngle; Xcoord_PhiAngle.resize(PhiAngleList.size());
		vector<vector<double> > Ycoord_PhiAngle; Ycoord_PhiAngle.resize(PhiAngleList.size());
		vector<vector<double> > Zcoord_PhiAngle; Zcoord_PhiAngle.resize(PhiAngleList.size());
		vector<vector<unsigned long> > IdPoint_PhiAngle; IdPoint_PhiAngle.resize(PhiAngleList.size());
		vector<vector<unsigned long> > IdDomain_PhiAngle; IdDomain_PhiAngle.resize(PhiAngleList.size());
		vector<vector<double> > Pressure_PhiAngle; Pressure_PhiAngle.resize(PhiAngleList.size());
		vector<vector<double> > FaceArea_PhiAngle; FaceArea_PhiAngle.resize(PhiAngleList.size());
		vector<vector<double> > EquivArea_PhiAngle; EquivArea_PhiAngle.resize(PhiAngleList.size());
		vector<vector<double> > TargetArea_PhiAngle; TargetArea_PhiAngle.resize(PhiAngleList.size());
		vector<vector<double> > NearFieldWeight_PhiAngle; NearFieldWeight_PhiAngle.resize(PhiAngleList.size());
		vector<vector<double> > Weight_PhiAngle; Weight_PhiAngle.resize(PhiAngleList.size());

		/*--- Distribute the values among the different PhiAngles ---*/
		for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
			for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
				if (AzimuthalAngle[iVertex] == PhiAngleList[iPhiAngle]) {
					Xcoord_PhiAngle[iPhiAngle].push_back(Xcoord[iVertex]);
					Ycoord_PhiAngle[iPhiAngle].push_back(Ycoord[iVertex]);
					Zcoord_PhiAngle[iPhiAngle].push_back(Zcoord[iVertex]);
					IdPoint_PhiAngle[iPhiAngle].push_back(IdPoint[iVertex]);
					IdDomain_PhiAngle[iPhiAngle].push_back(IdDomain[iVertex]);
					Pressure_PhiAngle[iPhiAngle].push_back(Pressure[iVertex]);
					FaceArea_PhiAngle[iPhiAngle].push_back(FaceArea[iVertex]);
					EquivArea_PhiAngle[iPhiAngle].push_back(EquivArea[iVertex]);
					TargetArea_PhiAngle[iPhiAngle].push_back(TargetArea[iVertex]);
					NearFieldWeight_PhiAngle[iPhiAngle].push_back(NearFieldWeight[iVertex]);
					Weight_PhiAngle[iPhiAngle].push_back(Weight[iVertex]);
				}

		/*--- Order the arrays (x Coordinate, Pressure, Point, and Domain) ---*/
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
			for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++)
				for (jVertex = 0; jVertex < Xcoord_PhiAngle[iPhiAngle].size() - 1 - iVertex; jVertex++)
					if (Xcoord_PhiAngle[iPhiAngle][jVertex] > Xcoord_PhiAngle[iPhiAngle][jVertex+1]) {
						auxXCoord = Xcoord_PhiAngle[iPhiAngle][jVertex]; Xcoord_PhiAngle[iPhiAngle][jVertex] = Xcoord_PhiAngle[iPhiAngle][jVertex+1]; Xcoord_PhiAngle[iPhiAngle][jVertex+1] = auxXCoord;
						auxYCoord = Ycoord_PhiAngle[iPhiAngle][jVertex]; Ycoord_PhiAngle[iPhiAngle][jVertex] = Ycoord_PhiAngle[iPhiAngle][jVertex+1]; Ycoord_PhiAngle[iPhiAngle][jVertex+1] = auxYCoord;
						auxZCoord = Zcoord_PhiAngle[iPhiAngle][jVertex]; Zcoord_PhiAngle[iPhiAngle][jVertex] = Zcoord_PhiAngle[iPhiAngle][jVertex+1]; Zcoord_PhiAngle[iPhiAngle][jVertex+1] = auxZCoord;
						auxPress = Pressure_PhiAngle[iPhiAngle][jVertex]; Pressure_PhiAngle[iPhiAngle][jVertex] = Pressure_PhiAngle[iPhiAngle][jVertex+1]; Pressure_PhiAngle[iPhiAngle][jVertex+1] = auxPress;
						auxArea = FaceArea_PhiAngle[iPhiAngle][jVertex]; FaceArea_PhiAngle[iPhiAngle][jVertex] = FaceArea_PhiAngle[iPhiAngle][jVertex+1]; FaceArea_PhiAngle[iPhiAngle][jVertex+1] = auxArea;
						auxPoint = IdPoint_PhiAngle[iPhiAngle][jVertex]; IdPoint_PhiAngle[iPhiAngle][jVertex] = IdPoint_PhiAngle[iPhiAngle][jVertex+1]; IdPoint_PhiAngle[iPhiAngle][jVertex+1] = auxPoint;
						auxDomain = IdDomain_PhiAngle[iPhiAngle][jVertex]; IdDomain_PhiAngle[iPhiAngle][jVertex] = IdDomain_PhiAngle[iPhiAngle][jVertex+1]; IdDomain_PhiAngle[iPhiAngle][jVertex+1] = auxDomain;
					}


		/*--- Check that all the azimuth lists have the same size ---*/
		unsigned short nVertex = Xcoord_PhiAngle[0].size();
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
			unsigned short nVertex_aux = Xcoord_PhiAngle[iPhiAngle].size();
			if (nVertex_aux != nVertex) cout <<"Be careful!!! one azimuth list is shorter than the other"<< endl;
			nVertex = min(nVertex, nVertex_aux);
		}

		/*--- Compute equivalent area distribution at each azimuth angle ---*/
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
			EquivArea_PhiAngle[iPhiAngle][0] = 0.0;
			for (iVertex = 1; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
				EquivArea_PhiAngle[iPhiAngle][iVertex] = 0.0;

				Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][iVertex]*sin(AoA);

				for (jVertex = 0; jVertex < iVertex-1; jVertex++) {

					Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex]*sin(AoA);
					jp1Coord = Xcoord_PhiAngle[iPhiAngle][jVertex+1]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex+1]*sin(AoA);

					jFunction = factor*(Pressure_PhiAngle[iPhiAngle][jVertex] - Pressure_Inf)*sqrt(Coord_i-Coord_j);
					jp1Function = factor*(Pressure_PhiAngle[iPhiAngle][jVertex+1] - Pressure_Inf)*sqrt(Coord_i-jp1Coord);

					DeltaX = (jp1Coord-Coord_j);
					MeanFuntion = 0.5*(jp1Function + jFunction);
					EquivArea_PhiAngle[iPhiAngle][iVertex] += DeltaX * MeanFuntion;
				}
			}
		}

		/*--- Create a file with the equivalent area distribution at each azimuthal angle ---*/
		NearFieldEA_file.precision(15);
		NearFieldEA_file.open("NearFieldEA.plt", ios::out);
		NearFieldEA_file << "TITLE = \"Nearfield Equivalent Area at each azimuthal angle \"" << endl;
		NearFieldEA_file << "VARIABLES = \"Coord (local to the near-field cylinder)\"";

		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
			NearFieldEA_file << ", \"Equivalent Area, Phi= " << PhiAngleList[iPhiAngle] << " deg.\"";
		}

		NearFieldEA_file << endl;
		for (iVertex = 0; iVertex < EquivArea_PhiAngle[0].size(); iVertex++) {
			double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
			double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);
			NearFieldEA_file << scientific << XcoordRot - XcoordRot_init;
			for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
				NearFieldEA_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex];
			}
			NearFieldEA_file << endl;
		}
		NearFieldEA_file.close();

		/*--- Read target equivalent area from the configuration file,
		 this first implementation requires a complete table (same as the original
		 EA table). so... no interpolation. ---*/

		vector<vector<double> > TargetArea_PhiAngle_Trans;
		TargetEA_file.open("TargetEA.dat", ios::in);

		if (TargetEA_file.fail()) {
			if (iExtIter == 0) { cout << "There is no Target Equivalent Area file (TargetEA.dat)!!"<< endl;
			cout << "Using default parameters (Target Equiv Area = 0.0)" << endl;
			}
			/*--- Set the table to 0 ---*/
			for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
				for (iVertex = 0; iVertex < TargetArea_PhiAngle[iPhiAngle].size(); iVertex++)
					TargetArea_PhiAngle[iPhiAngle][iVertex] = 0.0;
		}
		else {

			/*--- skip header lines ---*/
			string line;
			getline(TargetEA_file, line);
			getline(TargetEA_file, line);

			while (TargetEA_file) {

				string line;
				getline(TargetEA_file, line);
				istringstream is(line);
				vector<double> row;
				unsigned short iter = 0;

				while (is.good()) {
					string token;
					getline(is,token,',');

					istringstream js(token);

					double data;
					js >> data;

					/*--- The first element in the table is the coordinate ---*/
					if (iter != 0) row.push_back(data);
					iter++;
				}
				TargetArea_PhiAngle_Trans.push_back(row);
			}

			for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
				for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++)
					TargetArea_PhiAngle[iPhiAngle][iVertex] = TargetArea_PhiAngle_Trans[iVertex][iPhiAngle];
		}

		/*--- Divide by the number of Phi angles in the nearfield ---*/
		double PhiFactor = 1.0/double(PhiAngleList.size());

		/*--- Evaluate the objective function ---*/
		InverseDesign = 0;
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
			for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
				Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;
				Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];

				double Difference = EquivArea_PhiAngle[iPhiAngle][iVertex]-TargetArea_PhiAngle[iPhiAngle][iVertex];
				if ((Coord_i < XCoordBegin_OF) || (Coord_i > XCoordEnd_OF)) Difference = 0.0;

				InverseDesign += PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*Difference*Difference;

			}

		/*--- Evaluate the weight of the nearfield pressure (adjoint input) ---*/
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
			for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
				Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];
				NearFieldWeight_PhiAngle[iPhiAngle][iVertex] = 0.0;
				for (jVertex = iVertex; jVertex < EquivArea_PhiAngle[iPhiAngle].size(); jVertex++) {
					Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex];
					Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;

					double Difference = EquivArea_PhiAngle[iPhiAngle][jVertex]-TargetArea_PhiAngle[iPhiAngle][jVertex];
					if ((Coord_j < XCoordBegin_OF) || (Coord_j > XCoordEnd_OF)) Difference = 0.0;

					NearFieldWeight_PhiAngle[iPhiAngle][iVertex] += PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*2.0*Difference*factor*sqrt(Coord_j-Coord_i);
				}
			}		

		/*--- Write the Nearfield pressure at each Azimuthal PhiAngle ---*/
		EquivArea_file.precision(15);
		EquivArea_file.open("nearfield_flow.plt", ios::out);
		EquivArea_file << "TITLE = \"SU2 Equivalent area computation at each azimuthal angle \"" << endl;
		EquivArea_file << "VARIABLES = \"Coord (local to the near-field cylinder)\",\"Equivalent area\",\"Target equivalent area\",\"NearField weight\",\"Pressure coefficient\"" << endl;

		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
			EquivArea_file << fixed << "ZONE T= \"Azimuthal angle " << PhiAngleList[iPhiAngle] << " deg.\"" << endl;
			for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++) {

				double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
				double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);

				EquivArea_file << scientific << XcoordRot - XcoordRot_init << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex]
				                                                                                                    << ", " << TargetArea_PhiAngle[iPhiAngle][iVertex] << ", " << NearFieldWeight_PhiAngle[iPhiAngle][iVertex] << ", " <<
				                                                                                                    (Pressure_PhiAngle[iPhiAngle][iVertex]-Pressure_Inf)/Pressure_Inf << endl;
			}
		}

		EquivArea_file.close();

		/*--- Write Weight file for adjoint computation ---*/
		FuncGrad_file.precision(15);
		FuncGrad_file.open("WeightNF.dat", ios::out);

		FuncGrad_file << scientific << "-1.0";
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
			FuncGrad_file << scientific << "\t" << PhiAngleList[iPhiAngle];
		FuncGrad_file << endl;

		for (iVertex = 0; iVertex < NearFieldWeight_PhiAngle[0].size(); iVertex++) {
			double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
			FuncGrad_file << scientific << XcoordRot;
			for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
				FuncGrad_file << scientific << "\t" << NearFieldWeight_PhiAngle[iPhiAngle][iVertex];
			FuncGrad_file << endl;
		}
		FuncGrad_file.close();		

		/*--- Delete structures ---*/
		delete [] Xcoord; delete [] Ycoord; delete [] Zcoord; 
		delete [] AzimuthalAngle; delete [] IdPoint; delete [] IdDomain;
		delete [] Pressure; delete [] FaceArea;
		delete [] EquivArea; delete [] TargetArea;
		delete [] NearFieldWeight; delete [] Weight;

	}

#ifdef NO_MPI

	/*--- Store the value of the NearField coefficient ---*/
	solver_container->SetTotal_CEquivArea(InverseDesign);

#else

	/*--- Send the value of the NearField coefficient to all the processors ---*/
	MPI::COMM_WORLD.Bcast (&InverseDesign, 1, MPI::DOUBLE, MASTER_NODE);

	/*--- Store the value of the NearField coefficient ---*/
	solver_container->SetTotal_CEquivArea(InverseDesign);

#endif

}
