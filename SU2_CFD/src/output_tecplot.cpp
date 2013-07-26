/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2013 Aerospace Design Laboratory
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
string AssembleVariableNames(bool GridMovement, bool Incompressible, unsigned short KindSolver, unsigned short nVar_Consv, unsigned short dims, unsigned short *NVar);

void COutput::SetTecplot_ASCII(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol) {
    
	/*--- Local variables and initialization ---*/
	unsigned short iDim, iVar, nDim = geometry->GetnDim();
	unsigned short Kind_Solver = config->GetKind_Solver();
    
	unsigned long iPoint, iElem, iNode;
	unsigned long iExtIter = config->GetExtIter();
  unsigned long *LocalIndex = NULL;
  bool *SurfacePoint = NULL;
  
	bool grid_movement  = config->GetGrid_Movement();
	bool adjoint = config->GetAdjoint();
    
	char cstr[200], buffer[50];
	string filename;
    
	/*--- Write file name with extension ---*/
  if (surf_sol) {
    if (adjoint)
      filename = config->GetSurfAdjCoeff_FileName();
    else
      filename = config->GetSurfFlowCoeff_FileName();
  }
  else {
    if (adjoint)
      filename = config->GetAdj_FileName();
    else
      filename = config->GetFlow_FileName();
  }
  
	if (Kind_Solver == LINEAR_ELASTICITY)
		filename = config->GetStructure_FileName().c_str();
	if (Kind_Solver == WAVE_EQUATION)
		filename = config->GetWave_FileName().c_str();
	if ((Kind_Solver == WAVE_EQUATION) && (Kind_Solver == ADJ_AEROACOUSTIC_EULER))
		filename = config->GetAdjWave_FileName().c_str();
	if (Kind_Solver == ELECTRIC_POTENTIAL)
		filename = config->GetStructure_FileName().c_str();
	if (Kind_Solver == PLASMA_EULER) {
		if (val_iZone == 0) Kind_Solver = PLASMA_EULER;
		if (val_iZone == 1) Kind_Solver = ELECTRIC_POTENTIAL;
	}
	if (Kind_Solver == PLASMA_NAVIER_STOKES) {
		if (val_iZone == 0) Kind_Solver = PLASMA_NAVIER_STOKES;
		if (val_iZone == 1) Kind_Solver = ELECTRIC_POTENTIAL;
	}
    
#ifndef NO_MPI
	/*--- Remove the domain number from the surface csv filename ---*/
	int nProcessor = MPI::COMM_WORLD.Get_size();
	if (nProcessor > 1) filename.erase (filename.end()-2, filename.end());
#endif
    
	strcpy (cstr, filename.c_str());
	if (Kind_Solver == ELECTRIC_POTENTIAL) strcpy (cstr, config->GetStructure_FileName().c_str());
    
	/*--- Special cases where a number needs to be appended to the file name. ---*/
	if ((Kind_Solver == EULER || Kind_Solver == NAVIER_STOKES || Kind_Solver == RANS) &&
        (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL)) {
		sprintf (buffer, "_%d", int(val_iZone));
		strcat(cstr,buffer);
	}
    
	/*--- Special cases where a number needs to be appended to the file name. ---*/
	if (((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) &&
        (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL)) {
		sprintf (buffer, "_%d", int(val_iZone));
		strcat(cstr,buffer);
	}
    
    /*--- Special cases where a number needs to be appended to the file name. ---*/
	if (((Kind_Solver == ELECTRIC_POTENTIAL) &&( (Kind_Solver == PLASMA_EULER) || (Kind_Solver == PLASMA_NAVIER_STOKES)) )
        && config->GetUnsteady_Simulation()) {
		sprintf (buffer, "_%d", int(iExtIter));
		strcat(cstr,buffer);
	}
    
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.dat", int(val_iZone));
		if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.dat", int(val_iZone));
        
	} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		if (int(iExtIter) < 10) sprintf (buffer, "_0000%d.dat", int(iExtIter));
		if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.dat", int(iExtIter));
		if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.dat", int(iExtIter));
		if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.dat", int(iExtIter));
		if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.dat", int(iExtIter));
	} else {
		sprintf (buffer, ".dat");
	}
    
	strcat(cstr,buffer);
    
	/*--- Open Tecplot ASCII file and write the header. ---*/
	ofstream Tecplot_File;
	Tecplot_File.open(cstr, ios::out);
  Tecplot_File.precision(6);
  if (surf_sol) Tecplot_File << "TITLE = \"Visualization of the surface solution\"" << endl;
  else Tecplot_File << "TITLE = \"Visualization of the volumetric solution\"" << endl;

	/*--- Prepare the variables lists. ---*/
	if (nDim == 2) {
		Tecplot_File << "VARIABLES = \"x\",\"y\"";
	} else {
		Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\"";
	}
  
  /*--- Write the list of the fields in the restart file.
   Without including the PointID---*/
  if (config->GetKind_SU2() == SU2_SOL) {
    
    /*--- If SU2_SOL called this routine, we already have a set of output
     variables with the appropriate string tags stored in the config class. ---*/
    nVar_Total = config->fields.size() - 1;
    for (unsigned short iField = 1; iField < config->fields.size(); iField++) {
      Tecplot_File << config->fields[iField];
    }
    Tecplot_File << endl;
    
  } else {
    
    /*--- Add names for conservative and residual variables ---*/
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
      Tecplot_File << ",\"Conservative_" << iVar+1 << "\"";
    }
    if (config->GetWrt_Residuals()) {
      for (iVar = 0; iVar < nVar_Consv; iVar++) {
        Tecplot_File << ",\"Residual_" << iVar+1 << "\"";
      }
    }
    
    /*--- Add names for any extra variables (this will need to be adjusted). ---*/
    if (grid_movement) {
      if (nDim == 2) {
        Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\"";
      } else {
        Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\",\"Grid_Velz\"";
      }
    }
    
    if ((Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
      Tecplot_File << ",\"Density\"";
    }
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
        (Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
      Tecplot_File << ",\"Pressure\",\"Pressure_Coefficient\",\"Mach\"";
    }
    
    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
        (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
      Tecplot_File << ", \"Temperature\", \"Laminar_Viscosity\", \"Skin_Friction_Coefficient\", \"Heat_Transfer\", \"Y_Plus\"";
    }
    
    if ((Kind_Solver == RANS) || (Kind_Solver == FREE_SURFACE_RANS)) {
      Tecplot_File << ",\"Eddy_Viscosity\"";
    }
    
    if ((Kind_Solver == PLASMA_EULER) || (Kind_Solver == PLASMA_NAVIER_STOKES)) {
      unsigned short iSpecies;
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
        Tecplot_File << ",\"Pressure_" << iSpecies << "\"";
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
        Tecplot_File << ",\"Temperature_" << iSpecies << "\"";
      for (iSpecies = 0; iSpecies < config->GetnDiatomics(); iSpecies++)
        Tecplot_File << ",\"TemperatureVib_" << iSpecies << "\"";
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
        Tecplot_File << ",\"Mach_" << iSpecies << "\"";
    }
    
    if (Kind_Solver == PLASMA_NAVIER_STOKES) {
      unsigned short iSpecies;
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
        Tecplot_File << ",\"LaminaryViscosity_" << iSpecies << "\"";
      
      if ( Kind_Solver == PLASMA_NAVIER_STOKES  && (config->GetMagnetic_Force() == YES) && (geometry->GetnDim() == 3)) {
        for (iDim = 0; iDim < nDim; iDim++)
          Tecplot_File << ",\"Magnet_Field" << iDim << "\"";
      }
    }
    
    if (Kind_Solver == ELECTRIC_POTENTIAL) {
      unsigned short iDim;
      for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
        Tecplot_File << ",\"ElectricField_" << iDim+1 << "\"";
    }
    
    if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) ||
        (Kind_Solver == ADJ_FREE_SURFACE_EULER) || (Kind_Solver == ADJ_FREE_SURFACE_NAVIER_STOKES) ||
        (Kind_Solver == ADJ_FREE_SURFACE_RANS) || (Kind_Solver == ADJ_PLASMA_EULER) || (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES)) {
      Tecplot_File << ", \"Surface_Sensitivity\"";
    }
    
    Tecplot_File << endl;
    
  }
  
  /*--- If it's a surface output, print only the points 
   that are in the element list, change the numbering ---*/
  
  if (surf_sol) {
        
    LocalIndex = new unsigned long [nGlobal_Poin+1];
    SurfacePoint = new bool [nGlobal_Poin+1];

    for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) SurfacePoint[iPoint] = false;

    for(iElem = 0; iElem < nGlobal_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      SurfacePoint[Conn_Line[iNode+0]] = true;
      SurfacePoint[Conn_Line[iNode+1]] = true;
    }
    for(iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      SurfacePoint[Conn_BoundTria[iNode+0]] = true;
      SurfacePoint[Conn_BoundTria[iNode+1]] = true;
      SurfacePoint[Conn_BoundTria[iNode+2]] = true;
    }
    for(iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      SurfacePoint[Conn_BoundQuad[iNode+0]] = true;
      SurfacePoint[Conn_BoundQuad[iNode+1]] = true;
      SurfacePoint[Conn_BoundQuad[iNode+2]] = true;
      SurfacePoint[Conn_BoundQuad[iNode+3]] = true;
    }
    
    nSurf_Poin = 0;
    for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) {
      LocalIndex[iPoint] = 0;
      if (SurfacePoint[iPoint]) { nSurf_Poin++; LocalIndex[iPoint] = nSurf_Poin; }
    }
    
  }
  
  /*--- Write the header ---*/
  Tecplot_File << "ZONE ";
	if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady())
		Tecplot_File << "STRANDID="<<int(iExtIter+1)<<", SOLUTIONTIME="<<config->GetDelta_UnstTime()*iExtIter<<", ";
  
	if (nDim == 2) {
    if (surf_sol) Tecplot_File << "NODES= "<< nSurf_Poin <<", ELEMENTS= "<< nSurf_Elem <<", DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
    else Tecplot_File << "NODES= "<< nGlobal_Poin <<", ELEMENTS= "<< nGlobal_Elem <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
	} else {
    if (surf_sol) Tecplot_File << "NODES= "<< nSurf_Poin<<", ELEMENTS= "<< nSurf_Elem <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
    else Tecplot_File << "NODES= "<< nGlobal_Poin <<", ELEMENTS= "<< nGlobal_Elem <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
	}
  
	/*--- Write surface and volumetric solution data. ---*/
  
  for (iPoint = 0; iPoint < nGlobal_Poin; iPoint++) {
    
    if (surf_sol) {
      
      if (LocalIndex[iPoint+1] != 0) {
        
        /*--- Write the node coordinates ---*/
        for(iDim = 0; iDim < nDim; iDim++)
          Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";
        
        /*--- Loop over the vars/residuals and write the values to file ---*/
        for (iVar = 0; iVar < nVar_Total; iVar++)
          Tecplot_File << scientific << Data[iVar][iPoint] << "\t";
        
        Tecplot_File << endl;

      }
      
    } else {
      
      /*--- Write the node coordinates ---*/
      for(iDim = 0; iDim < nDim; iDim++)
        Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";
      
      /*--- Loop over the vars/residuals and write the values to file ---*/
      for (iVar = 0; iVar < nVar_Total; iVar++)
        Tecplot_File << scientific << Data[iVar][iPoint] << "\t";
      
      Tecplot_File << endl;
      
    }

  }
  

	/*--- Write connectivity data. ---*/
  
  if (surf_sol) {
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      Tecplot_File << LocalIndex[Conn_Line[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_Line[iNode+1]] << "\n";
      
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+1]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+2]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+2]] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+1]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+2]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+3]] << "\n";
    }
    
  }
  else {
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Tria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      Tecplot_File << Conn_Tria[iNode+0] << "\t";
      Tecplot_File << Conn_Tria[iNode+1] << "\t";
      Tecplot_File << Conn_Tria[iNode+2] << "\t";
      Tecplot_File << Conn_Tria[iNode+2] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Quad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      Tecplot_File << Conn_Quad[iNode+0] << "\t";
      Tecplot_File << Conn_Quad[iNode+1] << "\t";
      Tecplot_File << Conn_Quad[iNode+2] << "\t";
      Tecplot_File << Conn_Quad[iNode+3] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Tetr; iElem++) {
      iNode = iElem*N_POINTS_TETRAHEDRON;
      Tecplot_File << Conn_Tetr[iNode+0] << "\t" << Conn_Tetr[iNode+1] << "\t";
      Tecplot_File << Conn_Tetr[iNode+2] << "\t" << Conn_Tetr[iNode+2] << "\t";
      Tecplot_File << Conn_Tetr[iNode+3] << "\t" << Conn_Tetr[iNode+3] << "\t";
      Tecplot_File << Conn_Tetr[iNode+3] << "\t" << Conn_Tetr[iNode+3] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Hexa; iElem++) {
      iNode = iElem*N_POINTS_HEXAHEDRON;
      Tecplot_File << Conn_Hexa[iNode+0] << "\t" << Conn_Hexa[iNode+1] << "\t";
      Tecplot_File << Conn_Hexa[iNode+2] << "\t" << Conn_Hexa[iNode+3] << "\t";
      Tecplot_File << Conn_Hexa[iNode+4] << "\t" << Conn_Hexa[iNode+5] << "\t";
      Tecplot_File << Conn_Hexa[iNode+6] << "\t" << Conn_Hexa[iNode+7] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Wedg; iElem++) {
      iNode = iElem*N_POINTS_WEDGE;
      Tecplot_File << Conn_Wedg[iNode+0] << "\t" << Conn_Wedg[iNode+1] << "\t";
      Tecplot_File << Conn_Wedg[iNode+1] << "\t" << Conn_Wedg[iNode+2] << "\t";
      Tecplot_File << Conn_Wedg[iNode+3] << "\t" << Conn_Wedg[iNode+4] << "\t";
      Tecplot_File << Conn_Wedg[iNode+4] << "\t" << Conn_Wedg[iNode+5] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Pyra; iElem++) {
      iNode = iElem*N_POINTS_PYRAMID;
      Tecplot_File << Conn_Pyra[iNode+0] << "\t" << Conn_Pyra[iNode+1] << "\t";
      Tecplot_File << Conn_Pyra[iNode+2] << "\t" << Conn_Pyra[iNode+3] << "\t";
      Tecplot_File << Conn_Pyra[iNode+4] << "\t" << Conn_Pyra[iNode+4] << "\t";
      Tecplot_File << Conn_Pyra[iNode+4] << "\t" << Conn_Pyra[iNode+4] << "\n";
    }
  }
    
	Tecplot_File.close();
  
  if (surf_sol) delete [] LocalIndex;

}

void COutput::SetTecplot_Mesh(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
    
#ifndef NO_TECIO
    
	double   t;
	INTEGER4 i, N, err, Debug, NPts, NElm, IsDouble, KMax;
	INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, StrandID, ParentZn, FileType;
	INTEGER4 *ShareFromZone, IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
	string buffer, variables;
	stringstream file;
	bool first_zone = true;
	unsigned short dims = geometry->GetnDim();
	enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };
	enum	 ZoneType { ORDERED=0, FELINESEG=1, FETRIANGLE=2, FEQUADRILATERAL=3, FETETRAHEDRON=4, FEBRICK=5, FEPOLYGON=6, FEPOLYHEDRON=7 };
    
	/*--- Consistent data for Tecplot zones ---*/
	Debug						= 0;
	IsDouble					= 1;
	NPts						= (INTEGER4)nGlobal_Poin;
	t							= 0.0;//iExtIter*config->GetDelta_UnstTimeND();
	KMax						= 0;
	ICellMax					= 0;
	JCellMax					= 0;
	KCellMax					= 0;
	StrandID					= 0;//(INTEGER4)iExtIter;
	ParentZn					= 0;
	IsBlock						= 1;
	NumFaceConnections			= 0;
	FaceNeighborMode			= 0;
	ShareConnectivityFromZone	= 0;
    
	/*--- Write Tecplot solution file ---*/
	if (!wrote_base_file) {
        
		file.str(string());
		buffer = config->GetFlow_FileName();
        
#ifndef NO_MPI
		/*--- Remove the domain number from the filename ---*/
		int nProcessor = MPI::COMM_WORLD.Get_size();
		if (nProcessor > 1) buffer.erase(buffer.end()-2, buffer.end());
#endif
        
		file << buffer << ".mesh.plt";
		FileType = GRID;
        
		if (dims == 2) variables = "x y";
		else if (dims == 3) variables = "x y z";
		else cout << "Error: wrong number of dimentsions: " << dims << endl;
        
		/*--- Open Tecplot file ---*/
		err = TECINI112((char *)config->GetFlow_FileName().c_str(),
                        (char *)variables.c_str(),
                        (char *)file.str().c_str(),
                        (char *)".",
                        &FileType,
                        &Debug,
                        &IsDouble);
		if (err) cout << "Error in opening Tecplot file" << endl;
        
		first_zone = true;
		ShareFromZone = new INTEGER4[dims];
		for (i = 0; i < dims; i++) ShareFromZone[i] = 0;
        
		if (nGlobal_Tria > 0) {
            
			/*--- Write the zone header information ---*/
			ZoneType = FETRIANGLE; NElm = (INTEGER4)nGlobal_Tria; N = NElm*N_POINTS_TRIANGLE;
            
			err = TECZNE112((char*)"Triangle Elements",
                            &ZoneType,
                            &NPts,
                            &NElm,
                            &KMax,
                            &ICellMax,
                            &JCellMax,
                            &KCellMax,
                            &t,
                            &StrandID,
                            &ParentZn,
                            &IsBlock,
                            &NumFaceConnections,
                            &FaceNeighborMode,
                            0,         /* TotalNumFaceNodes */
                            0,         /* NumConnectedBoundaryFaces */
                            0,         /* TotalNumBoundaryConnections */
                            NULL,      /* PassiveVarList */
                            NULL,      /* ValueLocation */
                            ShareFromZone,      /* ShareVarFromZone */
                            &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;
            
			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
                
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
				if (geometry->GetnDim() == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					ShareFromZone[2] = 1;
				}
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				first_zone = false;
			}
            
			err = TECNOD112(Conn_Tria);
			if (err) cout << "Error writing connectivity to Tecplot file" << endl;
            
		}
		if (nGlobal_Quad > 0) {
            
			/*--- Write the zone header information ---*/
			ZoneType = FEQUADRILATERAL; NElm = (INTEGER4)nGlobal_Quad; N = NElm*N_POINTS_QUADRILATERAL;
            
			err = TECZNE112((char*)"Quadrilateral Elements",
                            &ZoneType,
                            &NPts,
                            &NElm,
                            &KMax,
                            &ICellMax,
                            &JCellMax,
                            &KCellMax,
                            &t,
                            &StrandID,
                            &ParentZn,
                            &IsBlock,
                            &NumFaceConnections,
                            &FaceNeighborMode,
                            0,         /* TotalNumFaceNodes */
                            0,         /* NumConnectedBoundaryFaces */
                            0,         /* TotalNumBoundaryConnections */
                            NULL,      /* PassiveVarList */
                            NULL,      /* ValueLocation */
                            ShareFromZone,      /* ShareVarFromZone */
                            &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;
            
			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
                
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
				if (geometry->GetnDim() == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					ShareFromZone[2] = 1;
				}
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				first_zone = false;
			}
            
			err = TECNOD112(Conn_Quad);
			if (err) cout << "Error writing connectivity to Tecplot file" << endl;
            
		}
		if (nGlobal_Tetr > 0) {
            
			/*--- Write the zone header information ---*/
			ZoneType = FETETRAHEDRON; NElm = (INTEGER4)nGlobal_Tetr; N = NElm*N_POINTS_TETRAHEDRON;
            
			err = TECZNE112((char*)"Tetrahedral Elements",
                            &ZoneType,
                            &NPts,
                            &NElm,
                            &KMax,
                            &ICellMax,
                            &JCellMax,
                            &KCellMax,
                            &t,
                            &StrandID,
                            &ParentZn,
                            &IsBlock,
                            &NumFaceConnections,
                            &FaceNeighborMode,
                            0,         /* TotalNumFaceNodes */
                            0,         /* NumConnectedBoundaryFaces */
                            0,         /* TotalNumBoundaryConnections */
                            NULL,      /* PassiveVarList */
                            NULL,      /* ValueLocation */
                            NULL,      /* ShareVarFromZone */
                            &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;
            
			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
                
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
				if (geometry->GetnDim() == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					ShareFromZone[2] = 1;
				}
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				first_zone = false;
			}
            
			err = TECNOD112(Conn_Tetr);
			if (err) cout << "Error writing connectivity to Tecplot file" << endl;
            
		}
		if (nGlobal_Hexa > 0) {
            
			/*--- Write the zone header information ---*/
			ZoneType = FEBRICK; NElm = (INTEGER4)nGlobal_Hexa; N = NElm*N_POINTS_HEXAHEDRON;
            
			err = TECZNE112((char*)"Hexahedral Elements",
                            &ZoneType,
                            &NPts,
                            &NElm,
                            &KMax,
                            &ICellMax,
                            &JCellMax,
                            &KCellMax,
                            &t,
                            &StrandID,
                            &ParentZn,
                            &IsBlock,
                            &NumFaceConnections,
                            &FaceNeighborMode,
                            0,         /* TotalNumFaceNodes */
                            0,         /* NumConnectedBoundaryFaces */
                            0,         /* TotalNumBoundaryConnections */
                            NULL,      /* PassiveVarList */
                            NULL,      /* ValueLocation */
                            NULL,      /* ShareVarFromZone */
                            &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;
            
			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
                
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
				if (geometry->GetnDim() == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					ShareFromZone[2] = 1;
				}
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				first_zone = false;
			}
            
			err = TECNOD112(Conn_Hexa);
			if (err) cout << "Error writing connectivity to Tecplot file" << endl;
            
		}
		if (nGlobal_Pyra > 0) {
			cout << "Pyramid element type not yet supported; no zone written." << endl;
		}
		if (nGlobal_Wedg > 0) {
			cout << "Wedge element type not yet supported; no zone written." << endl;
		}
		if (nGlobal_Line > 0) {
            
			/*--- Write the zone header information ---*/
			ZoneType = FELINESEG; NElm = (INTEGER4)nGlobal_Line; N = NElm*N_POINTS_LINE;
            
			err = TECZNE112((char*)"Line Elements",
                            &ZoneType,
                            &NPts,
                            &NElm,
                            &KMax,
                            &ICellMax,
                            &JCellMax,
                            &KCellMax,
                            &t,
                            &StrandID,
                            &ParentZn,
                            &IsBlock,
                            &NumFaceConnections,
                            &FaceNeighborMode,
                            0,         /* TotalNumFaceNodes */
                            0,         /* NumConnectedBoundaryFaces */
                            0,         /* TotalNumBoundaryConnections */
                            NULL,      /* PassiveVarList */
                            NULL,      /* ValueLocation */
                            NULL,      /* ShareVarFromZone */
                            &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;
            
			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
                
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
				if (geometry->GetnDim() == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					ShareFromZone[2] = 1;
				}
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				first_zone = false; 
			}
            
			err = TECNOD112(Conn_Line);
			if (err) cout << "Error writing connectivity to Tecplot file" << endl;
            
		}
        
		delete [] ShareFromZone;
		wrote_base_file = true;
        
		err = TECEND112();
		if (err) cout << "Error in closing Tecplot file" << endl;
        
	}
    
    
#else // Not built with Tecplot binary support
    
	cout << "Tecplot binary file requested but SU^2 was built without Tecio support. No file written" << "\n";
    
#endif
    
}

void COutput::SetTecplot_Solution(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
    
#ifndef NO_TECIO
    
	double   t;
	INTEGER4 i, N, iVar, err, Debug, NPts, NElm, IsDouble, KMax;
	INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, StrandID, ParentZn, FileType;
	INTEGER4 *ShareFromZone, IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
	string buffer, variables;
	stringstream file;
	bool first_zone = true, unsteady = config->GetUnsteady_Simulation(), GridMovement = config->GetGrid_Movement();
	unsigned long iExtIter = config->GetExtIter();
	unsigned short NVar, dims = geometry->GetnDim();
	enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };
	enum	 ZoneType { ORDERED=0, FELINESEG=1, FETRIANGLE=2, FEQUADRILATERAL=3, FETETRAHEDRON=4, FEBRICK=5, FEPOLYGON=6, FEPOLYHEDRON=7 };
    
	/*--- Consistent data for Tecplot zones ---*/
	Debug						= 0;
	IsDouble					= 1;
	NPts						= (INTEGER4)nGlobal_Poin;
	t							= iExtIter*config->GetDelta_UnstTimeND();
	KMax						= 0;
	ICellMax					= 0;
	JCellMax					= 0;
	KCellMax					= 0;
	StrandID					= (INTEGER4)iExtIter+1;
	ParentZn					= 0;
	IsBlock						= 1;
	NumFaceConnections			= 0;
	FaceNeighborMode			= 0;
	ShareConnectivityFromZone	= 0;
    
    
	file.str(string());
	buffer = config->GetFlow_FileName();
    
#ifndef NO_MPI
	/*--- Remove the domain number from the filename ---*/
	int nProcessor = MPI::COMM_WORLD.Get_size();
	if (nProcessor > 1) buffer.erase(buffer.end()-2, buffer.end());
#endif
    
	file << buffer;
    
	if (unsteady) {
		if (((int)iExtIter >= 0) && ((int)iExtIter < 10))			file << "_0000" << iExtIter;
		if (((int)iExtIter >= 10) && ((int)iExtIter < 100))		file << "_000" << iExtIter;
		if (((int)iExtIter >= 100) && ((int)iExtIter < 1000))		file << "_00" << iExtIter;
		if (((int)iExtIter >= 1000) && ((int)iExtIter < 10000))	file << "_0" << iExtIter;
		if ((int)iExtIter >= 10000)							file << iExtIter;
	}
	file << ".sol.plt";
	FileType = SOLUTION;
	variables = AssembleVariableNames(config->GetGrid_Movement(), config->GetIncompressible(), config->GetKind_Solver(), nVar_Consv, dims, &NVar);
    
	/*--- Open Tecplot file ---*/
	err = TECINI112((char *)config->GetFlow_FileName().c_str(),
                    (char *)variables.c_str(),
                    (char *)file.str().c_str(),
                    (char *)".",
                    &FileType,
                    &Debug,
                    &IsDouble);
	if (err) cout << "Error in opening Tecplot file" << endl;
    
	first_zone = true;
	ShareFromZone = new INTEGER4[NVar];
	for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;
    
	if (nGlobal_Tria > 0) {
        
		/*--- Write the zone header information ---*/
		ZoneType = FETRIANGLE; NElm = (INTEGER4)nGlobal_Tria; N = NElm*N_POINTS_TRIANGLE;
        
		err = TECZNE112((char*)"Triangle Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        ShareFromZone,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
		if (err) cout << "Error writing Tecplot zone data" << endl;
        
		/*--- write node coordinates and data if not done already---*/
		if (first_zone) {
            
			i = 0;
			if (GridMovement) {
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				if (dims == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					ShareFromZone[i++] = 1;
				}
			}
			for (iVar = 0; iVar < nVar_Total; iVar++) {
				err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing data to Tecplot file" << endl;
			}
            
			first_zone = false;
		}
        
	}
	if (nGlobal_Quad > 0) {
        
		/*--- Write the zone header information ---*/
		ZoneType = FEQUADRILATERAL; NElm = (INTEGER4)nGlobal_Quad; N = NElm*N_POINTS_QUADRILATERAL;
        
		err = TECZNE112((char*)"Quadrilateral Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        ShareFromZone,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
		if (err) cout << "Error writing Tecplot zone data" << endl;
        
		/*--- write node coordinates and data if not done already---*/
		if (first_zone) {
            
			i = 0;
			if (GridMovement) {
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				if (dims == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					ShareFromZone[i++] = 1;
				}
			}
			for (iVar = 0; iVar < nVar_Total; iVar++) {
				err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing data to Tecplot file" << endl;
			}
            
			first_zone = false;
		}
        
	}
	if (nGlobal_Tetr > 0) {
        
		/*--- Write the zone header information ---*/
		ZoneType = FETETRAHEDRON; NElm = (INTEGER4)nGlobal_Tetr; N = NElm*N_POINTS_TETRAHEDRON;
        
		err = TECZNE112((char*)"Tetrahedral Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        NULL,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
		if (err) cout << "Error writing Tecplot zone data" << endl;
        
		/*--- write node coordinates and data if not done already---*/
		if (first_zone) {
            
			i = 0;
			if (GridMovement) {
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				if (dims == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					ShareFromZone[i++] = 1;
				}
			}
			for (iVar = 0; iVar < nVar_Total; iVar++) {
				err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing data to Tecplot file" << endl;
			}
            
			first_zone = false;
		}
        
	}
	if (nGlobal_Hexa > 0) {
        
		/*--- Write the zone header information ---*/
		ZoneType = FEBRICK; NElm = (INTEGER4)nGlobal_Hexa; N = NElm*N_POINTS_HEXAHEDRON;
        
		err = TECZNE112((char*)"Hexahedral Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        NULL,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
		if (err) cout << "Error writing Tecplot zone data" << endl;
        
		/*--- write node coordinates and data if not done already---*/
		if (first_zone) {
            
			i = 0;
			if (GridMovement) {
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				if (dims == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					ShareFromZone[i++] = 1;
				}
			}
			for (iVar = 0; iVar < nVar_Total; iVar++) {
				err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing data to Tecplot file" << endl;
			}
            
			first_zone = false;
		}
        
	}
	if (nGlobal_Pyra > 0) {
		cout << "Pyramid element type not yet supported; no zone written." << endl;
	}
	if (nGlobal_Wedg > 0) {
		cout << "Wedge element type not yet supported; no zone written." << endl;
	}
	if (nGlobal_Line > 0) {
        
		/*--- Write the zone header information ---*/
		ZoneType = FELINESEG; NElm = (INTEGER4)nGlobal_Line; N = NElm*N_POINTS_LINE;
        
		err = TECZNE112((char*)"Line Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        NULL,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
		if (err) cout << "Error writing Tecplot zone data" << endl;
        
		/*--- write node coordinates and data if not done already---*/
		if (first_zone) {
            
			i = 0;
			if (GridMovement) {
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				if (dims == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					ShareFromZone[i++] = 1;
				}
			}
			for (iVar = 0; iVar < nVar_Total; iVar++) {
				err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
				if (err) cout << "Error writing data to Tecplot file" << endl;
			}
            
			first_zone = false;
		}
        
	}
    
	delete [] ShareFromZone;
    
	err = TECEND112();
	if (err) cout << "Error in closing Tecplot file" << endl;
    
#else // Not built with Tecplot binary support
    
	cout << "Tecplot binary file requested but SU^2 was built without Tecio support. No file written" << "\n"; 
    
#endif
    
}

string AssembleVariableNames(bool GridMovement, bool Incompressible, unsigned short KindSolver, unsigned short nVar_Consv, unsigned short dims, unsigned short *NVar) {
    
	stringstream variables;
	unsigned short iVar;
	*NVar = 0;
    
	variables.str(string());
	if (GridMovement) {
		variables << "x y "; *NVar += 2;
		if (dims == 3) {
			variables << "z "; *NVar += 1;
		}
	}
	for (iVar = 0; iVar < nVar_Consv; iVar++){
		variables << "Conservative_Variable_" << iVar+1 << " ";
		*NVar += 1;
	}
//	for (iVar = nVar_Consv; iVar < 2*nVar_Consv; iVar++) {
//		variables << "Conservative_Residual_" << iVar-nVar_Consv+1 << " ";
//		*NVar += 1;
//	}
	if (GridMovement) {
		variables << "Grid_Velocity_X " << "Grid_Velocity_Y "; *NVar += 2;
		if (dims == 3) {
			variables << "Grid_Velocity_Z "; *NVar += 1;
		}
	}
    
	variables << "Pressure " << "Mach "; *NVar += 2;
	if (!Incompressible) {
		switch (KindSolver) {
                
                /*--- Include temperature and laminar viscosity, if applicable ---*/
            case NAVIER_STOKES: variables << "Temperature " << "Viscosity "; *NVar += 2; break;
                
                /*--- Include eddy viscosity, if applicable ---*/
            case RANS: variables << "Temperature " << "Viscosity " << "Eddy_Viscosity "; *NVar += 3; break;
		}
	}
    
	return variables.str(); 
    
}
