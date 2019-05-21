/*!
 * \file output_inria.cpp
 * \brief Main subroutines for output in Inria format.
 * \author V. Menier, B. Munguía
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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


#include "../include/output_structure.hpp"

enum BCVAR  { bcMach, bcTemp, bcPres, bcDens, bcGoal };

void COutput::SetInriaRestart(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
	
  unsigned short nZone = geometry->GetnZone();
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned short nVar_Buf = nVar_Par-nDim;
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool grid_movement = config->GetGrid_Movement();
  bool dynamic_fem = (config->GetDynamic_Analysis() == DYNAMIC);
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);
  string filename;
  
  unsigned long OutSol,i, npoin = geometry->GetGlobal_nPointDomain();
  unsigned long myPoint, offset, Global_Index;
  int VarTyp[GmfMaxTyp];
  passivedouble bufDbl[GmfMaxTyp];
  char OutNam[1024], BasNam[1024];
  char *ptr=NULL;
	
  /*--- Retrieve filename from config ---*/
  
  if ((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem){
    filename = config->GetRestart_FEMFileName();
  } else {
    filename = config->GetRestart_FlowFileName();
  }
	
	/*--- Get output name *.solb ---*/
	
  strcpy(BasNam, filename.c_str());
  ptr = strstr(BasNam,".dat");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
  ptr = strstr(BasNam,".solb");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
  SPRINTF (OutNam, "%s.solb", BasNam);

  /*--- Open the restart file and write the solution. ---*/
	
	OutSol = GmfOpenMesh(OutNam,GmfWrite,GmfDouble,nDim);
	
	if ( !OutSol ) {
	  printf("\n\n   !!! Error !!!\n" );
      printf("Unable to open %s", OutNam);
      printf("Now exiting...\n\n");
      exit(EXIT_FAILURE);
	}
	
  /*--- Write the restart file ---*/

	for (iVar = 0; iVar < nVar_Buf; iVar++) {
		VarTyp[iVar] = GmfSca;
	}
	
	if ( !GmfSetKwd(OutSol, GmfSolAtVertices, npoin, nVar_Buf, VarTyp) ) {
	  printf("\n\n   !!! Error !!!\n" );
      printf("Unable to write %s", OutNam);
      printf("Now exiting...\n\n");
      exit(EXIT_FAILURE);
	}

  myPoint = 0;
  offset = 0;
  for (unsigned short iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
        
        /*--- Global Index of the current point. (note outer loop over procs) ---*/
        
        Global_Index = iPoint + offset;
        
        /*--- Only write original domain points, i.e., exclude any periodic
         or halo nodes, even if they are output in the viz. files. ---*/
        
        if (Global_Index < nPoint_Restart) {
                    
          myPoint++;
          
          /*--- Loop over the variables and write the values to file,
           excluding mesh coordinates [0 - (nDim-1)] ---*/
          
          for (iVar = 0; iVar < nVar_Buf; iVar++){
            bufDbl[iVar] = SU2_TYPE::GetValue(Parallel_Data[iVar+nDim][iPoint]);
          }
          GmfSetLin(OutSol, GmfSolAtVertices, bufDbl);
        }
      }
    }
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
    
  }
  
	if ( !GmfCloseMesh(OutSol) ) {
	  printf("\n\n   !!! Error !!!\n" );
      printf("Cannot close solution file %s.", OutNam);
      printf("Now exiting...\n\n");
      exit(EXIT_FAILURE);
	}  
	
}

/*

	Write solutions of interest : mach, dens, pres etc.

*/
void COutput::WriteInriaOutputs(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
	
  unsigned short nZone = geometry->GetnZone();
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned short nVar_First, nVar_Second, nVar_Consv_Par;
  unsigned long iPoint, iExtIter = config->GetExtIter();
  unsigned long FirstIndex = NONE, SecondIndex = NONE;
  bool grid_movement = config->GetGrid_Movement();
  bool dynamic_fem = (config->GetDynamic_Analysis() == DYNAMIC);
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);
  string filename;
  
  unsigned long OutMach, OutPres, OutMetr;
  unsigned long i, npoin = geometry->GetGlobal_nPointDomain();
  unsigned long myPoint, offset, Global_Index;
  int VarTyp[GmfMaxTyp];
  passivedouble bufDbl[GmfMaxTyp];
  char OutNam[1024], BasNam[1024];
  char *ptr=NULL;
	
  int NbrVar, idxVar;

  /*--- Use a switch statement to decide how many solver containers we have
   in this zone for output. ---*/
  
  switch (config->GetKind_Solver()) {
    case EULER : case NAVIER_STOKES: FirstIndex = FLOW_SOL; SecondIndex = NONE; break;
    case RANS : FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; break;
    default: SecondIndex = NONE; break;
  }
  
  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  nVar_Consv_Par = nVar_First + nVar_Second;
	
  /* Get indices of mach, pres, etc. in the solution array */
  unsigned short *TagBc;
  TagBc = new unsigned short [100];
	
  idxVar=0;
  idxVar += nDim + nVar_Consv_Par; // Add coordinates and conservative variables
  
  if (!config->GetLow_MemoryOutput()) {
    
    if (config->GetWrt_Limiters()) {
      idxVar += nVar_Consv_Par; // Add limiters
    }
    if (config->GetWrt_Residuals()) {
      idxVar += nVar_Consv_Par; // Add residuals
    }
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      TagBc[bcPres] = idxVar;
      TagBc[bcTemp] = idxVar+1;
      TagBc[bcMach] = idxVar+2;
      idxVar += 4; // Add pressure, temperature, Mach, Cp
    }
		
  }

  /* Get index of adaptation parameter if performing error estimation */

  if(config->GetError_Estimate() && config->GetKind_SU2() == SU2_MET){
    
    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      idxVar += nDim + 2; // Add laminar viscosity, skin friction, heat flux
      if(config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR) idxVar += 1; // Add buffet sensor
    }
        
    if (Kind_Solver == RANS) idxVar += 2; // Add y-plus, eddy viscosity
        
    if (config->GetWrt_SharpEdges()) idxVar += 1; // Add sharp edges
        
    if (config->GetKind_Trans_Model() == BC) idxVar += 1; // Add the intermittency for the BC trans. model
    
    if (config->GetKind_HybridRANSLES()!=NO_HYBRIDRANSLES) idxVar += 2; // Add DES length scale and wall distance
    
    if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS) idxVar += 1; // Add Roe dissipation

    TagBc[bcGoal] = idxVar;
    idxVar += 1; // Add adaptation parameter
  }
	
  /*--- Retrieve filename from config ---*/
  
  if ((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem){
    filename = config->GetRestart_FEMFileName();
  } else {
    filename = config->GetRestart_FlowFileName();
  }

  /*--- Append the zone number if multizone problems ---*/
  if (nZone > 1)
    filename= config->GetMultizone_FileName(filename, val_iZone);

  /*--- Number of variables (1 for sensor files). ---*/
  NbrVar = 1;
  VarTyp[0]  = GmfSca;
	
  /*--- Open the restart file and write the solution. ---*/

  /*--- Write MACH ---*/

  SPRINTF (OutNam, "mach.solb");
  OutMach = GmfOpenMesh(OutNam,GmfWrite,GmfDouble,nDim);
	
  if ( !OutMach ) {
    printf("\n\n   !!! Error !!!\n" );
    printf("Unable to open %s", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
  }	
	
  if ( !GmfSetKwd(OutMach, GmfSolAtVertices, npoin, NbrVar, VarTyp) ) {
    printf("\n\n   !!! Error !!!\n" );
    printf("Unable to write Mach");
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
  }

  myPoint = 0;
  offset = 0;
  for (unsigned short iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
        
        /*--- Global Index of the current point. (note outer loop over procs) ---*/
        
        Global_Index = iPoint + offset;
        
        /*--- Only write original domain points, i.e., exclude any periodic
         or halo nodes, even if they are output in the viz. files. ---*/
        
        if (Global_Index < nPoint_Restart) {
                    
          myPoint++;
          
          /*--- Loop over the variables and write the values to file ---*/
          
          iVar = TagBc[bcMach];
          bufDbl[0] = SU2_TYPE::GetValue(Parallel_Data[iVar][iPoint]);
          GmfSetLin(OutMach, GmfSolAtVertices, bufDbl);
        }
      }
    }
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
    
  }
	
  if ( !GmfCloseMesh(OutMach) ) {
    printf("\n\n   !!! Error !!!\n" );
    printf("Cannot close solution file");
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
  }
	
  /*--- Write PRES ---*/

  SPRINTF (OutNam, "pres.solb");
  OutPres = GmfOpenMesh(OutNam,GmfWrite,GmfDouble,nDim);
	
  if ( !OutPres ) {
    printf("\n\n   !!! Error !!!\n" );
    printf("Unable to open %s", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
  }
	
  if ( !GmfSetKwd(OutPres, GmfSolAtVertices, npoin, NbrVar, VarTyp) ) {
    printf("\n\n   !!! Error !!!\n" );
    printf("Unable to write pressure");
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
  }

  myPoint = 0;
  offset = 0;
  for (unsigned short iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
        
        /*--- Global Index of the current point. (note outer loop over procs) ---*/
        
        Global_Index = iPoint + offset;
        
        /*--- Only write original domain points, i.e., exclude any periodic
         or halo nodes, even if they are output in the viz. files. ---*/
        
        if (Global_Index < nPoint_Restart) {
                    
          myPoint++;
          
          /*--- Loop over the variables and write the values to file ---*/
          
          iVar = TagBc[bcPres];
          bufDbl[0] = SU2_TYPE::GetValue(Parallel_Data[iVar][iPoint]);
          GmfSetLin(OutPres, GmfSolAtVertices, bufDbl);
        }
      }
    }
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
    
  }
		
	/*--- Close files ---*/
  	
  if ( !GmfCloseMesh(OutPres) ) {
    printf("\n\n   !!! Error !!!\n" );
    printf("Cannot close solution file");
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
  }

  /*--- Write metric ---*/

  if(config->GetError_Estimate() && config->GetKind_SU2() == SU2_MET){

    /*--- Write metric tensor ---*/

    SPRINTF (OutNam, "metr.solb");
    OutMetr = GmfOpenMesh(OutNam,GmfWrite,GmfDouble,nDim);
  
    if ( !OutMetr ) {
      printf("\n\n   !!! Error !!!\n" );
      printf("Unable to open %s", OutNam);
      printf("Now exiting...\n\n");
      exit(EXIT_FAILURE);
    }

    unsigned short iMetr, nMetr;
    if(nDim == 2) nMetr = 3;
    else          nMetr = 6;

    VarTyp[0] = GmfSymMat;
  
    if ( !GmfSetKwd(OutMetr, GmfSolAtVertices, npoin, NbrVar, VarTyp) ) {
      printf("\n\n   !!! Error !!!\n" );
      printf("Unable to write metric");
      printf("Now exiting...\n\n");
      exit(EXIT_FAILURE);
    }

    myPoint = 0;
    offset = 0;
    for (unsigned short iProcessor = 0; iProcessor < size; iProcessor++) {
      if (rank == iProcessor) {
        for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
          
          /*--- Global Index of the current point. (note outer loop over procs) ---*/
          
          Global_Index = iPoint + offset;
          
          /*--- Only write original domain points, i.e., exclude any periodic
           or halo nodes, even if they are output in the viz. files. ---*/
          
          if (Global_Index < nPoint_Restart) {
                      
            myPoint++;
            
            /*--- Loop over the variables and write the values to file ---*/
            
            iVar = TagBc[bcGoal];
            for (iMetr = 0; iMetr < nMetr; iMetr++)
              bufDbl[iMetr] = SU2_TYPE::GetValue(Parallel_Data[iVar+iMetr][iPoint]);
            GmfSetLin(OutMetr, GmfSolAtVertices, bufDbl);
          }
        }
      }
#ifdef HAVE_MPI
      SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
      SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
      
    }
    
    /*--- Close files ---*/
    
    if ( !GmfCloseMesh(OutMetr) ) {
      printf("\n\n   !!! Error !!!\n" );
      printf("Cannot close solution file");
      printf("Now exiting...\n\n");
      exit(EXIT_FAILURE);
    }

  }
	
  delete [] TagBc;
	
}




void COutput::SetInriaMesh(CConfig *config, CGeometry *geometry) {
  
  char cstr[MAX_STRING_SIZE], out_file[MAX_STRING_SIZE];
  unsigned long iElem, iPoint, iElem_Bound, nElem_Bound_, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], iNode, nElem;
  unsigned short iMarker, iDim, nDim = geometry->GetnDim(), iChar, iPeriodic, nPeriodic = 0, VTK_Type, nMarker_;
  su2double *center, *angles, *transl;
  ofstream output_file;
  ifstream input_file;
  string Grid_Marker, text_line, Marker_Tag, str;
  string::size_type position;
	
  unsigned short nMarker = config->GetnMarker_All();
  unsigned long cptElem = 0, nTri=0, nLin=0, nQua=0;
  unsigned long myPoint, offset, Global_Index;

  unsigned long OutMsh,i;
  int iVer,iTri,iEfr,iTet;
  passivedouble bufDbl[8];
  char OutNam[2014];
  int bufInt[8];
	
  unsigned long *PointSurface=NULL;
  unsigned long nPointSurface=0;
	
  CPrimalGrid* bnd = NULL;

  /*--- Read the name of the output and input file ---*/
  str = config->GetMesh_Out_FileName();

  unsigned short lastindex = str.find_last_of(".");
  str = str.substr(0, lastindex);

  SPRINTF (OutNam, "%s.meshb", str.c_str());
	
  OutMsh = GmfOpenMesh(OutNam,GmfWrite,GmfDouble,nDim);
  if ( !OutMsh ) {
    printf("  ## ERROR: Cannot open mesh file %s ! \n",OutNam);
	return;
  }
  
  /*--- Write vertices ---*/
	
  GmfSetKwd(OutMsh, GmfVertices, nPoint_Restart);

  myPoint = 0;
  offset = 0;
  for (unsigned short iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
        
        /*--- Global Index of the current point. (note outer loop over procs) ---*/
        
        Global_Index = iPoint + offset;
        
        /*--- Only write original domain points, i.e., exclude any periodic
         or halo nodes, even if they are output in the viz. files. ---*/
        
        if (Global_Index < nPoint_Restart) {
                    
          myPoint++;
          
          /*--- Loop over the variables and write the values to file ---*/
          
          for (iDim = 0; iDim < nDim; iDim++)
            bufDbl[iDim] = SU2_TYPE::GetValue(Parallel_Data[iDim][iPoint]);

          if ( nDim == 2 ) {
            GmfSetLin(OutMsh, GmfVertices,bufDbl[0],bufDbl[1],0); 
          }
          else {
            GmfSetLin(OutMsh, GmfVertices,bufDbl[0],bufDbl[1],bufDbl[2],0); 
          }
        }
      }
    }
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
    
  }

  /*--- Write 2D elements ---
    Note: in 3D, triangles/quads are boundary markers
  */
	
  if ( nDim == 2 ){
		
    /*--- Write triangles ---*/
		
	GmfSetKwd(OutMsh, GmfTriangles, nParallel_Tria);
  	for (iElem = 0; iElem < nParallel_Tria; iElem++) {
  	  iNode = iElem*N_POINTS_TRIANGLE;
	  GmfSetLin(OutMsh, GmfTriangles,Conn_Tria_Par[iNode+0],Conn_Tria_Par[iNode+1],Conn_Tria_Par[iNode+2], 1);  
  	}	

	/*--- Write quadrilaterals ---*/
		
	if ( nParallel_Quad > 0  ) {
	  GmfSetKwd(OutMsh, GmfQuadrilaterals, nParallel_Quad);
	  for (iElem = 0; iElem < nParallel_Quad; iElem++) {
  		  iNode = iElem*N_POINTS_QUADRILATERAL;
				GmfSetLin(OutMsh, GmfQuadrilaterals,Conn_Quad_Par[iNode+0],Conn_Quad_Par[iNode+1],Conn_Quad_Par[iNode+2], Conn_Quad_Par[iNode+3], 1);  
  		}
	  }
	
	}
	
	/*--- Write tetrahedra ---*/
	
	if ( nParallel_Tetr > 0  ) {
		GmfSetKwd(OutMsh, GmfTetrahedra, nParallel_Tetr);
		for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
	    iNode = iElem*N_POINTS_TETRAHEDRON;
			GmfSetLin(OutMsh, GmfTetrahedra,Conn_Tetr_Par[iNode+0],Conn_Tetr_Par[iNode+1],Conn_Tetr_Par[iNode+2], Conn_Tetr_Par[iNode+3], 1); 
	  }
	}
	
	/*--- Write hexahedra ---*/
	
	if ( nParallel_Hexa > 0 ) {
		GmfSetKwd(OutMsh, GmfHexahedra, nParallel_Hexa);
		for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
	    iNode = iElem*N_POINTS_HEXAHEDRON;
			GmfSetLin(OutMsh, GmfHexahedra,Conn_Hexa_Par[iNode+0],Conn_Hexa_Par[iNode+1], Conn_Hexa_Par[iNode+2], Conn_Hexa_Par[iNode+3], Conn_Hexa_Par[iNode+4],Conn_Hexa_Par[iNode+5],Conn_Hexa_Par[iNode+6], Conn_Hexa_Par[iNode+7],  1); 
	  }
	}
	
	/*--- Write prisms ---*/
	
	if ( nParallel_Pris > 0 ) {
		GmfSetKwd(OutMsh, GmfPrisms, nParallel_Pris);
		for (iElem = 0; iElem < nParallel_Pris; iElem++) {
	    iNode = iElem*N_POINTS_PRISM;
			GmfSetLin(OutMsh, GmfPrisms,Conn_Pris_Par[iNode+0],Conn_Pris_Par[iNode+1], Conn_Pris_Par[iNode+2], Conn_Pris_Par[iNode+3], Conn_Pris_Par[iNode+4],Conn_Pris_Par[iNode+5],  1); 
	  }
	}
	
	/*--- Write pyramids ---*/
	
	if ( nParallel_Pyra > 0 ) {
		GmfSetKwd(OutMsh, GmfPyramids, nParallel_Pyra);
		for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
	    iNode = iElem*N_POINTS_PYRAMID;
	  	GmfSetLin(OutMsh, GmfPyramids,Conn_Pyra_Par[iNode+0],Conn_Pyra_Par[iNode+1], Conn_Pyra_Par[iNode+2], Conn_Pyra_Par[iNode+3], Conn_Pyra_Par[iNode+4],1); 
		}
	}
	
	
	/* --- Boundary elements ---*/
	
	/*--- Get surface points ---*/
	
	nPointSurface = 0;
	PointSurface = new unsigned long[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    if (geometry->node[iPoint]->GetBoundary()) {
      PointSurface[nPointSurface] = iPoint;
      nPointSurface++;
    }
	}
	
	/*--- Count elements ---*/
	
	nLin = nTri = nQua = 0;
	
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
			bnd = geometry->bound[iMarker][iElem];
			switch ( bnd->GetVTK_Type() ) {
				case LINE:          nLin++; break;
				case TRIANGLE:      nTri++; break;
				case QUADRILATERAL: nQua++; break;
			}
		}
	}
	
	
	/*--- Write edges ---*/
	
	if ( nLin > 0 ) {
		
		GmfSetKwd(OutMsh, GmfEdges, nLin);
		
		cptElem = 0;
		
		for (iMarker = 0; iMarker < nMarker; iMarker++) {
			for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
				
				bnd = geometry->bound[iMarker][iElem];
				
				if ( bnd->GetVTK_Type() != LINE ) 
					continue;
				
				cptElem++;

				
				GmfSetLin(OutMsh, GmfEdges, bnd->GetNode(0)+1, bnd->GetNode(1)+1, iMarker+2); 	
			}
		}
	
		if ( cptElem != nLin ) {
			cout << "  !! Error Inria output:  Inconsistent number of edges\n" << endl;
			exit(EXIT_FAILURE);
		}
		
	}
	
	/*--- Write triangles ---*/
	
	if ( nTri > 0 ) {
		
		GmfSetKwd(OutMsh, GmfTriangles, nTri);
		
		cptElem = 0;
		
		for (iMarker = 0; iMarker < nMarker; iMarker++) {
			for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
				
				bnd = geometry->bound[iMarker][iElem];
				
				if ( bnd->GetVTK_Type() != TRIANGLE ) 
					continue;
				
				
				if ( cptElem < 10 ) 
					cout << bnd->GetNode(0)+1 << " " << bnd->GetNode(1)+1 << " " << bnd->GetNode(2)+1 << endl;
				
				cptElem++;
				
				GmfSetLin(OutMsh, GmfTriangles, bnd->GetNode(0)+1, bnd->GetNode(1)+1, 
                                        bnd->GetNode(2)+1, iMarker+2); 	
			}
		}
	
		if ( cptElem != nTri ) {
			cout << "  !! Error Inria output:  Inconsistent number of triangles\n" << endl;
			exit(EXIT_FAILURE);
		}
		
	}
	
	/*--- Write quadrilaterals ---*/
	
	if ( nQua > 0 ) {
		
		GmfSetKwd(OutMsh, GmfQuadrilaterals, nQua);
		
		cptElem = 0;
		
		for (iMarker = 0; iMarker < nMarker; iMarker++) {
			for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
				
				bnd = geometry->bound[iMarker][iElem];
				
				if ( bnd->GetVTK_Type() != QUADRILATERAL ) 
					continue;
				
				cptElem++;
				
				GmfSetLin(OutMsh, GmfQuadrilaterals, PointSurface[bnd->GetNode(0)]+1, PointSurface[bnd->GetNode(1)]+1, 
                                             PointSurface[bnd->GetNode(2)]+1, PointSurface[bnd->GetNode(3)]+1, iMarker+2); 	
			}
		}
	
		if ( cptElem != nQua ) {
			cout << "  !! Error Inria output:  Inconsistent number of quadrilaterals\n" << endl;
			exit(EXIT_FAILURE);
		}
		
	}
	
	if ( PointSurface )
		delete [] PointSurface;
	
	GmfCloseMesh(OutMsh);
  
}

