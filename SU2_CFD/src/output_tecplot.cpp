/*!
 * \file output_tecplot.cpp
 * \brief Main subroutines for output solver information.
 * \author F. Palacios, T. Economon, M. Colonno
 * \version 6.2.0 "Falcon"
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

void COutput::SetTecplotASCII(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol) {
  
  unsigned short iDim, iVar, nDim = geometry->GetnDim();
  unsigned short Kind_Solver = config->GetKind_Solver();
  
  unsigned long iPoint, iElem, iNode;
  unsigned long iExtIter = config->GetExtIter();
  unsigned long *LocalIndex = NULL;
  bool *SurfacePoint = NULL;
  
  bool grid_movement  = config->GetGrid_Movement();
  bool adjoint = config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint();

  char cstr[200], buffer[50];
  string filename;
  
  /*--- Write file name with extension ---*/
  
  if (surf_sol) {
    if (adjoint) filename = config->GetSurfAdjCoeff_FileName();
    else filename = config->GetSurfFlowCoeff_FileName();
  }
  else {
    if (adjoint)
    filename = config->GetAdj_FileName();
    else filename = config->GetFlow_FileName();
  }
  
  if (Kind_Solver == FEM_ELASTICITY) {
    if (surf_sol) filename = config->GetSurfStructure_FileName().c_str();
    else filename = config->GetStructure_FileName().c_str();
  }
  
  if (Kind_Solver == HEAT_EQUATION_FVM) {
    if (surf_sol) filename = config->GetSurfHeat_FileName().c_str();
    else filename = config->GetHeat_FileName().c_str();
  }
  
  if (config->GetKind_SU2() == SU2_DOT) {
    if (surf_sol) filename = config->GetSurfSens_FileName();
    else filename = config->GetVolSens_FileName();
  }
  strcpy (cstr, filename.c_str());
  
  /*--- Special cases where a number needs to be appended to the file name. ---*/
  
  if ((Kind_Solver == EULER || Kind_Solver == NAVIER_STOKES || Kind_Solver == RANS ||
       Kind_Solver == ADJ_EULER || Kind_Solver == ADJ_NAVIER_STOKES || Kind_Solver == ADJ_RANS ||
       Kind_Solver == DISC_ADJ_EULER || Kind_Solver == DISC_ADJ_NAVIER_STOKES || Kind_Solver == DISC_ADJ_RANS ||
       Kind_Solver == HEAT_EQUATION_FVM) &&
      (val_nZone > 1) ) {
    SPRINTF (buffer, "_%d", SU2_TYPE::Int(val_iZone));
    strcat(cstr, buffer);
  }
  
 if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady() && config->GetUnsteady_Simulation() != HARMONIC_BALANCE) {
    if (SU2_TYPE::Int(iExtIter) < 10) SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10) && (SU2_TYPE::Int(iExtIter) < 100)) SPRINTF (buffer, "_000%d.dat", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100) && (SU2_TYPE::Int(iExtIter) < 1000)) SPRINTF (buffer, "_00%d.dat", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.dat", SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(iExtIter));
  }
  else { SPRINTF (buffer, ".dat"); }
  
  strcat(cstr, buffer);
  
  /*--- Open Tecplot ASCII file and write the header. ---*/
  ofstream Tecplot_File;
  Tecplot_File.open(cstr, ios::out);
  Tecplot_File.precision(6);
  if (surf_sol) Tecplot_File << "TITLE = \"Visualization of the surface solution\"" << endl;
  else Tecplot_File << "TITLE = \"Visualization of the volumetric solution\"" << endl;
  
  /*--- Prepare the variable lists. ---*/
  
  /*--- Write the list of the fields in the restart file.
   Without including the PointID---*/
  if ((config->GetKind_SU2() == SU2_SOL) || (config->GetKind_SU2() == SU2_DOT)) {
    
    /*--- If SU2_SOL called this routine, we already have a set of output
     variables with the appropriate string tags stored in the config class. ---*/
    Tecplot_File << "VARIABLES = ";
    nVar_Total = config->fields.size() - 1;
    for (unsigned short iField = 1; iField < config->fields.size(); iField++) {
      Tecplot_File << config->fields[iField] << " ";
    }
    Tecplot_File << endl;
    
  } else {
    
    if (nDim == 2) {
      Tecplot_File << "VARIABLES = \"x\",\"y\"";
    } else {
      Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\"";
    }
    
    /*--- Add names for conservative and residual variables ---*/
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
      Tecplot_File << ",\"Conservative_" << iVar+1 << "\"";
    }
    
    if (!config->GetLow_MemoryOutput()) {
      
      if (config->GetWrt_Limiters()) {
        for (iVar = 0; iVar < nVar_Consv; iVar++) {
          Tecplot_File << ",\"Limiter_" << iVar+1 << "\"";
        }
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
      
      if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
        Tecplot_File << ",\"Pressure\",\"Temperature\",\"C<sub>p</sub>\",\"Mach\"";
      }
      
      if (((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))) {
        if (nDim == 2) Tecplot_File << ", \"<greek>m</greek>\", \"C<sub>f</sub>_x\", \"C<sub>f</sub>_y\", \"h\", \"y<sup>+</sup>\"";
        else Tecplot_File << ", \"<greek>m</greek>\", \"C<sub>f</sub>_x\", \"C<sub>f</sub>_y\", \"C<sub>f</sub>_z\", \"h\", \"y<sup>+</sup>\"";
      }
      
      if (Kind_Solver == RANS) {
        Tecplot_File << ", \"<greek>m</greek><sub>t</sub>\"";
      }
      
      if (config->GetWrt_SharpEdges()) {
        if (((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)))  {
          Tecplot_File << ", \"Sharp_Edge_Dist\"";
        }
      }
      
      if (( Kind_Solver == ADJ_EULER              ) ||
          ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
          ( Kind_Solver == ADJ_RANS               )   ) {
        Tecplot_File << ", \"Surface_Sensitivity\", \"Solution_Sensor\"";
      }

      if (( Kind_Solver == DISC_ADJ_EULER              ) ||
          ( Kind_Solver == DISC_ADJ_NAVIER_STOKES      ) ||
          ( Kind_Solver == DISC_ADJ_RANS               )) {
        Tecplot_File << ", \"Surface_Sensitivity\", \"Sensitivity_x\", \"Sensitivity_y\"";
        if (geometry->GetnDim() == 3) {
          Tecplot_File << ",\"Sensitivity_z\"";
        }
      }
      
      if (Kind_Solver == FEM_ELASTICITY) {
        Tecplot_File << ", \"Von_Mises_Stress\"";
      }
      
      if (config->GetExtraOutput()) {
        string *headings = NULL;
        //if (Kind_Solver == RANS) {
        headings = solver[TURB_SOL]->OutputHeadingNames;
        //}
        for (iVar = 0; iVar < nVar_Extra; iVar++) {
          //Tecplot_File << ", \"ExtraOutput_" << iVar+1<<"\"";
          if (headings == NULL) {
            Tecplot_File << ", \"ExtraOutput_" << iVar+1<<"\"";
          } else {
            Tecplot_File << ", \""<< headings[iVar] <<"\"";
          }
        }
      }
    }
    
    Tecplot_File << endl;
    
  }
  
  /*--- If it's a surface output, print only the points
   that are in the element list, change the numbering ---*/
  
  if (surf_sol) {
    
    LocalIndex = new unsigned long [nGlobal_Poin+1];
    SurfacePoint = new bool [nGlobal_Poin+1];
    
    for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) SurfacePoint[iPoint] = false;
    
    for (iElem = 0; iElem < nGlobal_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      SurfacePoint[Conn_Line[iNode+0]] = true;
      SurfacePoint[Conn_Line[iNode+1]] = true;
    }
    for (iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      SurfacePoint[Conn_BoundTria[iNode+0]] = true;
      SurfacePoint[Conn_BoundTria[iNode+1]] = true;
      SurfacePoint[Conn_BoundTria[iNode+2]] = true;
    }
    for (iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
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
  if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    Tecplot_File << "STRANDID="<<SU2_TYPE::Int(iExtIter+1)<<", SOLUTIONTIME="<<config->GetDelta_UnstTime()*iExtIter<<", ";
  } else if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
    su2double period = config->GetHarmonicBalance_Period();
    su2double deltaT = period/(su2double)(config->GetnTimeInstances());
    Tecplot_File << "STRANDID="<<SU2_TYPE::Int(val_iZone+1)<<", SOLUTIONTIME="<<deltaT*val_iZone<<", ";
  }
  
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
        if ((config->GetKind_SU2() != SU2_SOL) && (config->GetKind_SU2() != SU2_DOT)) {
          for (iDim = 0; iDim < nDim; iDim++)
          Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";
        }
        
        /*--- Loop over the vars/residuals and write the values to file ---*/
        for (iVar = 0; iVar < nVar_Total; iVar++)
        Tecplot_File << scientific << Data[iVar][iPoint] << "\t";
        
        Tecplot_File << endl;
        
      }
      
    } else {
      
      /*--- Write the node coordinates ---*/
      if ((config->GetKind_SU2() != SU2_SOL) && (config->GetKind_SU2() != SU2_DOT)) {
        for (iDim = 0; iDim < nDim; iDim++)
        Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";
      }
      
      /*--- Loop over the vars/residuals and write the values to file ---*/
      for (iVar = 0; iVar < nVar_Total; iVar++)
      Tecplot_File << scientific << Data[iVar][iPoint] << "\t";
      
      Tecplot_File << endl;
      
    }
    
  }
  
  
  /*--- Write connectivity data. ---*/
  if (surf_sol) {
    
    for (iElem = 0; iElem < nGlobal_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      Tecplot_File << LocalIndex[Conn_Line[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_Line[iNode+1]] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+1]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+2]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+2]] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+1]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+2]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+3]] << "\n";
    }
    
  } else {
    
    for (iElem = 0; iElem < nGlobal_Tria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      Tecplot_File << Conn_Tria[iNode+0] << "\t";
      Tecplot_File << Conn_Tria[iNode+1] << "\t";
      Tecplot_File << Conn_Tria[iNode+2] << "\t";
      Tecplot_File << Conn_Tria[iNode+2] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_Quad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      Tecplot_File << Conn_Quad[iNode+0] << "\t";
      Tecplot_File << Conn_Quad[iNode+1] << "\t";
      Tecplot_File << Conn_Quad[iNode+2] << "\t";
      Tecplot_File << Conn_Quad[iNode+3] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_Tetr; iElem++) {
      iNode = iElem*N_POINTS_TETRAHEDRON;
      Tecplot_File << Conn_Tetr[iNode+0] << "\t" << Conn_Tetr[iNode+1] << "\t";
      Tecplot_File << Conn_Tetr[iNode+2] << "\t" << Conn_Tetr[iNode+2] << "\t";
      Tecplot_File << Conn_Tetr[iNode+3] << "\t" << Conn_Tetr[iNode+3] << "\t";
      Tecplot_File << Conn_Tetr[iNode+3] << "\t" << Conn_Tetr[iNode+3] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_Hexa; iElem++) {
      iNode = iElem*N_POINTS_HEXAHEDRON;
      Tecplot_File << Conn_Hexa[iNode+0] << "\t" << Conn_Hexa[iNode+1] << "\t";
      Tecplot_File << Conn_Hexa[iNode+2] << "\t" << Conn_Hexa[iNode+3] << "\t";
      Tecplot_File << Conn_Hexa[iNode+4] << "\t" << Conn_Hexa[iNode+5] << "\t";
      Tecplot_File << Conn_Hexa[iNode+6] << "\t" << Conn_Hexa[iNode+7] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_Pris; iElem++) {
      iNode = iElem*N_POINTS_PRISM;
      Tecplot_File << Conn_Pris[iNode+0] << "\t" << Conn_Pris[iNode+1] << "\t";
      Tecplot_File << Conn_Pris[iNode+1] << "\t" << Conn_Pris[iNode+2] << "\t";
      Tecplot_File << Conn_Pris[iNode+3] << "\t" << Conn_Pris[iNode+4] << "\t";
      Tecplot_File << Conn_Pris[iNode+4] << "\t" << Conn_Pris[iNode+5] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_Pyra; iElem++) {
      iNode = iElem*N_POINTS_PYRAMID;
      Tecplot_File << Conn_Pyra[iNode+0] << "\t" << Conn_Pyra[iNode+1] << "\t";
      Tecplot_File << Conn_Pyra[iNode+2] << "\t" << Conn_Pyra[iNode+3] << "\t";
      Tecplot_File << Conn_Pyra[iNode+4] << "\t" << Conn_Pyra[iNode+4] << "\t";
      Tecplot_File << Conn_Pyra[iNode+4] << "\t" << Conn_Pyra[iNode+4] << "\n";
    }
  }
  
  Tecplot_File.close();
  
  if (surf_sol) {
    delete [] LocalIndex;
    delete[] SurfacePoint;
  }
  
}

void COutput::SetTecplotASCII_Mesh(CConfig *config, CGeometry *geometry, unsigned short val_iZone, bool surf_sol, bool new_file) {
  
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iElem, iNode;
  unsigned long *LocalIndex = NULL;
  bool *SurfacePoint = NULL;
  char cstr[200];
  ofstream Tecplot_File;

  if (surf_sol) strcpy(cstr, "surface_grid");
  else strcpy(cstr, "volumetric_grid");
  
  if (config->GetnZone() > 1){
    char appstr[200];
    SPRINTF(appstr, "_%u", val_iZone);
    strcat(cstr, appstr);
  }

  strcat(cstr,".dat");

  /*--- Open Tecplot ASCII file and write the header. ---*/
  
  if (new_file) {
    Tecplot_File.open(cstr, ios::out);
    Tecplot_File.precision(6);
    if (surf_sol) Tecplot_File << "TITLE = \"Visualization of the surface solution\"" << endl;
    else Tecplot_File << "TITLE = \"Visualization of the volumetric solution\"" << endl;
    
    if (nDim == 2) Tecplot_File << "VARIABLES = \"x\",\"y\"";
    else Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\"";
  }
  else Tecplot_File.open(cstr, ios::out | ios::app);
  Tecplot_File << endl;
  
  /*--- If it's a surface output, print only the points
   that are in the element list, change the numbering ---*/
  
  if (surf_sol) {
    
    LocalIndex = new unsigned long [nGlobal_Poin+1];
    SurfacePoint = new bool [nGlobal_Poin+1];
    
    for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) SurfacePoint[iPoint] = false;
    
    for (iElem = 0; iElem < nGlobal_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      SurfacePoint[Conn_Line[iNode+0]] = true;
      SurfacePoint[Conn_Line[iNode+1]] = true;
    }
    for (iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      SurfacePoint[Conn_BoundTria[iNode+0]] = true;
      SurfacePoint[Conn_BoundTria[iNode+1]] = true;
      SurfacePoint[Conn_BoundTria[iNode+2]] = true;
    }
    for (iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
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
  
  Tecplot_File << "ZONE T= ";
  if (new_file) Tecplot_File << "\"Original grid\", C=BLACK, ";
  else Tecplot_File << "\"Deformed grid\", C=RED, ";
  
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
          for (iDim = 0; iDim < nDim; iDim++)
            Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";
        
        Tecplot_File << endl;
        
      }
      
    } else {
      
      /*--- Write the node coordinates ---*/

      for (iDim = 0; iDim < nDim; iDim++)
          Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";
      
      
      Tecplot_File << endl;
      
    }
    
  }
  
  
  /*--- Write connectivity data. ---*/
  
  if (surf_sol) {
    
    for (iElem = 0; iElem < nGlobal_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      Tecplot_File << LocalIndex[Conn_Line[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_Line[iNode+1]] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+1]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+2]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+2]] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+1]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+2]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+3]] << "\n";
    }
    
  } else {
    
    for (iElem = 0; iElem < nGlobal_Tria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      Tecplot_File << Conn_Tria[iNode+0] << "\t";
      Tecplot_File << Conn_Tria[iNode+1] << "\t";
      Tecplot_File << Conn_Tria[iNode+2] << "\t";
      Tecplot_File << Conn_Tria[iNode+2] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_Quad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      Tecplot_File << Conn_Quad[iNode+0] << "\t";
      Tecplot_File << Conn_Quad[iNode+1] << "\t";
      Tecplot_File << Conn_Quad[iNode+2] << "\t";
      Tecplot_File << Conn_Quad[iNode+3] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_Tetr; iElem++) {
      iNode = iElem*N_POINTS_TETRAHEDRON;
      Tecplot_File << Conn_Tetr[iNode+0] << "\t" << Conn_Tetr[iNode+1] << "\t";
      Tecplot_File << Conn_Tetr[iNode+2] << "\t" << Conn_Tetr[iNode+2] << "\t";
      Tecplot_File << Conn_Tetr[iNode+3] << "\t" << Conn_Tetr[iNode+3] << "\t";
      Tecplot_File << Conn_Tetr[iNode+3] << "\t" << Conn_Tetr[iNode+3] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_Hexa; iElem++) {
      iNode = iElem*N_POINTS_HEXAHEDRON;
      Tecplot_File << Conn_Hexa[iNode+0] << "\t" << Conn_Hexa[iNode+1] << "\t";
      Tecplot_File << Conn_Hexa[iNode+2] << "\t" << Conn_Hexa[iNode+3] << "\t";
      Tecplot_File << Conn_Hexa[iNode+4] << "\t" << Conn_Hexa[iNode+5] << "\t";
      Tecplot_File << Conn_Hexa[iNode+6] << "\t" << Conn_Hexa[iNode+7] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_Pris; iElem++) {
      iNode = iElem*N_POINTS_PRISM;
      Tecplot_File << Conn_Pris[iNode+0] << "\t" << Conn_Pris[iNode+1] << "\t";
      Tecplot_File << Conn_Pris[iNode+1] << "\t" << Conn_Pris[iNode+2] << "\t";
      Tecplot_File << Conn_Pris[iNode+3] << "\t" << Conn_Pris[iNode+4] << "\t";
      Tecplot_File << Conn_Pris[iNode+4] << "\t" << Conn_Pris[iNode+5] << "\n";
    }
    
    for (iElem = 0; iElem < nGlobal_Pyra; iElem++) {
      iNode = iElem*N_POINTS_PYRAMID;
      Tecplot_File << Conn_Pyra[iNode+0] << "\t" << Conn_Pyra[iNode+1] << "\t";
      Tecplot_File << Conn_Pyra[iNode+2] << "\t" << Conn_Pyra[iNode+3] << "\t";
      Tecplot_File << Conn_Pyra[iNode+4] << "\t" << Conn_Pyra[iNode+4] << "\t";
      Tecplot_File << Conn_Pyra[iNode+4] << "\t" << Conn_Pyra[iNode+4] << "\n";
    }
  }
  
  Tecplot_File.close();
  
  if (surf_sol) {
  	delete [] LocalIndex;
  	delete [] SurfacePoint;
  }
  
}

void COutput::SetSTL_MeshASCII(CConfig *config, CGeometry *geometry) {
  
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iElem, iNode;
  char cstr[200];
  ofstream STL_File;
  su2double p[3] = {0.0,0.0,0.0}, u[3] = {0.0,0.0,0.0}, v[3] = {0.0,0.0,0.0}, n[3] = {0.0,0.0,0.0}, a;
  unsigned long Point_0, Point_1, Point_2;
  
  /*---	STL format:
   solid NAME
   ...
   facet normal 0.00 0.00 1.00
   outer loop
   vertex  2.00  2.00  0.00
   vertex -1.00  1.00  0.00
   vertex  0.00 -1.00  0.00
   endloop
   endfacet
   ...
   end solid
   ---*/
  
  if (nDim == 3) {
    
    strcpy(cstr, "surface_grid.stl");
    
    /*--- Open STL ASCII file and write the header. ---*/
    
    STL_File.open(cstr, ios::out);
    STL_File.precision(6);
    STL_File << "solid surface_mesh" << endl;
    
    for (iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
      
      iNode = iElem*N_POINTS_TRIANGLE;
      
      /*--- Compute Normal vectors ---*/
      
      Point_0 = Conn_BoundTria[iNode+0]-1;
      Point_1 = Conn_BoundTria[iNode+1]-1;
      Point_2 = Conn_BoundTria[iNode+2]-1;
      
      for (iDim = 0; iDim < nDim; iDim++) {
        p[0] = Coords[iDim][Point_0];
        p[1] = Coords[iDim][Point_1];
        p[2] = Coords[iDim][Point_2];
        u[iDim] = p[1]-p[0];
        v[iDim] = p[2]-p[0];
      }
      
      n[0] = u[1]*v[2]-u[2]*v[1];
      n[1] = u[2]*v[0]-u[0]*v[2];
      n[2] = u[0]*v[1]-u[1]*v[0];
      a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
      
      /*--- Print normal vector ---*/
      
      STL_File << "  facet normal ";
      for (iDim = 0; iDim < nDim; iDim++) {
        STL_File << n[iDim]/a << " ";
      }
      STL_File << endl;
      
      /*--- Print nodes for facet ---*/
      STL_File << "    outer loop" << endl;
      
      STL_File << "      vertex ";
      for (iDim = 0; iDim < nDim; iDim++) STL_File << Coords[iDim][Point_0] << " ";
      STL_File <<  endl;
      
      STL_File << "      vertex ";
      for (iDim = 0; iDim < nDim; iDim++) STL_File << Coords[iDim][Point_1] << " ";
      STL_File <<  endl;
      
      STL_File << "      vertex ";
      for (iDim = 0; iDim < nDim; iDim++) STL_File << Coords[iDim][Point_2] << " ";
      STL_File <<  endl;
      
      STL_File << "    endloop" << endl;
      STL_File << "  endfacet" << endl;
      
    }
    
    //  for (iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
    //      iNode = iElem*N_POINTS_QUADRILATERAL;
    //      STL_File << LocalIndex[Conn_BoundQuad[iNode+0]] << "\t";
    //      STL_File << LocalIndex[Conn_BoundQuad[iNode+1]] << "\t";
    //      STL_File << LocalIndex[Conn_BoundQuad[iNode+2]] << "\t";
    //      STL_File << LocalIndex[Conn_BoundQuad[iNode+3]] << "\n";
    //    }
    
    /*--- Done with Surface Mesh ---*/
    
    STL_File << "endsolid" << endl;
    
    STL_File.close();
    
    
  }
  
}

void COutput::SetCSV_MeshASCII(CConfig *config, CGeometry *geometry) {

	short iStation, nStation;
	unsigned short nDim = geometry->GetnDim();
	unsigned long iVertex;
	su2double *Plane_P0, *Plane_Normal;
	vector<su2double> Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil, Variable_Airfoil;
	ofstream csv_File;

	if (nDim == 3) {

		Plane_P0 = new su2double[3];
		Plane_Normal = new su2double[3];

		if (geometry->GetnDim() == 3) {

			nStation = config->GetnLocationStations();

			for (iStation = 0; iStation < nStation; iStation++) {

				/*--- Read the values from the config file ---*/

				Plane_Normal[0] = 0.0; Plane_P0[0] = 0.0;
				Plane_Normal[1] = 0.0; Plane_P0[1] = 0.0;
				Plane_Normal[2] = 0.0; Plane_P0[2] = 0.0;

     if (config->GetGeo_Description() == FUSELAGE) {
       Plane_Normal[0] = 1.0;
       Plane_P0[0] = config->GetLocationStations(iStation);
     }
     
     if (config->GetGeo_Description() == NACELLE) {
       Plane_Normal[0] = 0.0;
       Plane_Normal[1] = -sin(config->GetLocationStations(iStation)*PI_NUMBER/180.0);
       Plane_Normal[2] = cos(config->GetLocationStations(iStation)*PI_NUMBER/180.0);
     }
     
     if (config->GetGeo_Description() == WING) {
       Plane_Normal[1] = 1.0;
       Plane_P0[1] = config->GetLocationStations(iStation);
     }

				/*--- Compute the airfoil Stations (note that we feed in the Cp) ---*/

				geometry->ComputeAirfoil_Section(Plane_P0, Plane_Normal, -1E6, 1E6, -1E6, 1E6, -1E6, 1E6,
                                     NULL, Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil,
                                     Variable_Airfoil, true, config);

				if ((rank == MASTER_NODE) && (Xcoord_Airfoil.size() == 0)) {
					cout << "Please check the config file, the station (" << Plane_P0[0] << ", " << Plane_P0[1] << ", " << Plane_P0[2] << ") has not been detected." << endl;
				}

				/*--- Write Cp at each Station (csv format) ---*/

				if ((rank == MASTER_NODE) && (Xcoord_Airfoil.size() != 0)) {

					if (iStation == 0) csv_File.open("surface_grid.csv", ios::out);
					else csv_File.open("surface_grid.csv", ios::app);

					/*--- Coordinates value ---*/

					for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
						csv_File << Xcoord_Airfoil[iVertex] << " ," << Ycoord_Airfoil[iVertex] << " ," << Zcoord_Airfoil[iVertex];
						if (iVertex == 0) { if (iStation == 0) csv_File << ", 2"; else csv_File << ", 1"; }
						else csv_File << ", 0";
						csv_File << endl;
					}

					csv_File.close();

				}

			}

		}

		/*--- Delete dynamically allocated memory ---*/

		delete[] Plane_P0;
		delete[] Plane_Normal;

	}

}

namespace
{

  std::string GetTecplotFilename(CConfig *config, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol, const char *extension) {

    const bool adjoint = config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint();
    const unsigned short Kind_Solver = config->GetKind_Solver();
    string filename;

    if (config->GetKind_SU2() == SU2_DOT) {
      if (surf_sol) filename = config->GetSurfSens_FileName();
      else filename = config->GetVolSens_FileName();
    }
    else if (Kind_Solver == FEM_ELASTICITY) {
      if (surf_sol) filename = config->GetSurfStructure_FileName().c_str();
      else filename = config->GetStructure_FileName().c_str();
    }
    else if (surf_sol) {
      if (adjoint) filename = config->GetSurfAdjCoeff_FileName();
      else filename = config->GetSurfFlowCoeff_FileName();
    }
    else {
      if (adjoint)
        filename = config->GetAdj_FileName();
      else filename = config->GetFlow_FileName();
    }
    
    ostringstream string_stream;
    string_stream << filename;
    
    /*--- Special cases where a number needs to be appended to the file name. ---*/
    if ((Kind_Solver == EULER || Kind_Solver == NAVIER_STOKES || Kind_Solver == RANS ||
         Kind_Solver == ADJ_EULER || Kind_Solver == ADJ_NAVIER_STOKES || Kind_Solver == ADJ_RANS ||
         Kind_Solver == DISC_ADJ_EULER || Kind_Solver == DISC_ADJ_NAVIER_STOKES || Kind_Solver == DISC_ADJ_RANS) &&
        (val_nZone > 1) ) {
      string_stream << '_' << val_iZone;
    }
    
    const unsigned long iExtIter = config->GetExtIter();
    if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady() && config->GetUnsteady_Simulation() != HARMONIC_BALANCE) {
      string_stream << '_' << setfill('0') << setw(5) << iExtIter;
    }

    string_stream << extension;
    return string_stream.str();
  }

} /* namespace */

void COutput::WriteTecplotASCII_Parallel(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_nZone, unsigned short val_iInst, unsigned short val_nInst, bool surf_sol) {
  
  unsigned short iVar, nDim = geometry->GetnDim();
  
  unsigned long iPoint, iElem, iNode;
  unsigned long iExtIter = config->GetExtIter();
  
  int iProcessor;

  string filename = GetTecplotFilename(config, val_iZone, val_nZone, surf_sol, ".dat");
  ofstream Tecplot_File;
  
  /*--- Open Tecplot ASCII file and write the header. ---*/
  
  if (rank == MASTER_NODE) {
    Tecplot_File.open(filename.c_str(), ios::out);
    Tecplot_File.precision(6);
    if (surf_sol) Tecplot_File << "TITLE = \"Visualization of the surface solution\"" << endl;
    else Tecplot_File << "TITLE = \"Visualization of the volumetric solution\"" << endl;
    
    Tecplot_File << "VARIABLES = ";
    for (iVar = 0; iVar < Variable_Names.size()-1; iVar++) {
      Tecplot_File << "\"" << Variable_Names[iVar] << "\",";
    }
    Tecplot_File << "\"" << Variable_Names[Variable_Names.size()-1] << "\"" << endl;
    
    /*--- Write the header ---*/
    
    Tecplot_File << "ZONE ";
    if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
      Tecplot_File << "STRANDID="<<SU2_TYPE::Int(iExtIter+1)<<", SOLUTIONTIME="<<config->GetDelta_UnstTime()*iExtIter<<", ";
    } else if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
      /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
      su2double period = config->GetHarmonicBalance_Period();
      su2double deltaT = period/(su2double)(config->GetnTimeInstances());
      Tecplot_File << "STRANDID="<<SU2_TYPE::Int(val_iZone+1)<<", SOLUTIONTIME="<<deltaT*val_iZone<<", ";
    }
    if (nDim == 2) {
      if (surf_sol) Tecplot_File << "NODES= "<< nGlobal_Surf_Poin <<", ELEMENTS= "<< nSurf_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
      else Tecplot_File << "NODES= "<< nGlobal_Poin_Par <<", ELEMENTS= "<< nGlobal_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
    } else {
      if (surf_sol) Tecplot_File << "NODES= "<< nGlobal_Surf_Poin <<", ELEMENTS= "<< nSurf_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
      else Tecplot_File << "NODES= "<< nGlobal_Poin_Par <<", ELEMENTS= "<< nGlobal_Elem_Par <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
    }

    Tecplot_File.close();
    
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- Each processor opens the file. ---*/
  
  Tecplot_File.open(filename.c_str(), ios::out | ios::app);
  
  /*--- Write surface and volumetric solution data. ---*/
  
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      
        /*--- Write the node data from this proc ---*/
        
        if (surf_sol) {
            for (iPoint = 0; iPoint < nSurf_Poin_Par; iPoint++) {
        for (iVar = 0; iVar < nVar_Par; iVar++)
          Tecplot_File << scientific << Parallel_Surf_Data[iVar][iPoint] << "\t";
        Tecplot_File << endl;
            
          }
        } else {
          
           for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
          for (iVar = 0; iVar < nVar_Par; iVar++)
            Tecplot_File << scientific << Parallel_Data[iVar][iPoint] << "\t";
          Tecplot_File << endl;
        }
      }
    }
    Tecplot_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }
  
  /*--- Write connectivity data. ---*/
  
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      
      if (surf_sol) {
        
        for (iElem = 0; iElem < nParallel_Line; iElem++) {
          iNode = iElem*N_POINTS_LINE;
          Tecplot_File << Conn_BoundLine_Par[iNode+0] << "\t";
          Tecplot_File << Conn_BoundLine_Par[iNode+1] << "\n";
        }
        
        for (iElem = 0; iElem < nParallel_BoundTria; iElem++) {
          iNode = iElem*N_POINTS_TRIANGLE;
          Tecplot_File << Conn_BoundTria_Par[iNode+0] << "\t";
          Tecplot_File << Conn_BoundTria_Par[iNode+1] << "\t";
          Tecplot_File << Conn_BoundTria_Par[iNode+2] << "\t";
          Tecplot_File << Conn_BoundTria_Par[iNode+2] << "\n";
        }
        
        for (iElem = 0; iElem < nParallel_BoundQuad; iElem++) {
          iNode = iElem*N_POINTS_QUADRILATERAL;
          Tecplot_File << Conn_BoundQuad_Par[iNode+0] << "\t";
          Tecplot_File << Conn_BoundQuad_Par[iNode+1] << "\t";
          Tecplot_File << Conn_BoundQuad_Par[iNode+2] << "\t";
          Tecplot_File << Conn_BoundQuad_Par[iNode+3] << "\n";
        }
        
      } else {
        
      for (iElem = 0; iElem < nParallel_Tria; iElem++) {
        iNode = iElem*N_POINTS_TRIANGLE;
        Tecplot_File << Conn_Tria_Par[iNode+0] << "\t";
        Tecplot_File << Conn_Tria_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Tria_Par[iNode+2] << "\t";
        Tecplot_File << Conn_Tria_Par[iNode+2] << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Quad; iElem++) {
        iNode = iElem*N_POINTS_QUADRILATERAL;
        Tecplot_File << Conn_Quad_Par[iNode+0] << "\t";
        Tecplot_File << Conn_Quad_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Quad_Par[iNode+2] << "\t";
        Tecplot_File << Conn_Quad_Par[iNode+3] << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
        iNode = iElem*N_POINTS_TETRAHEDRON;
        Tecplot_File << Conn_Tetr_Par[iNode+0] << "\t" << Conn_Tetr_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Tetr_Par[iNode+2] << "\t" << Conn_Tetr_Par[iNode+2] << "\t";
        Tecplot_File << Conn_Tetr_Par[iNode+3] << "\t" << Conn_Tetr_Par[iNode+3] << "\t";
        Tecplot_File << Conn_Tetr_Par[iNode+3] << "\t" << Conn_Tetr_Par[iNode+3] << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
        iNode = iElem*N_POINTS_HEXAHEDRON;
        Tecplot_File << Conn_Hexa_Par[iNode+0] << "\t" << Conn_Hexa_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Hexa_Par[iNode+2] << "\t" << Conn_Hexa_Par[iNode+3] << "\t";
        Tecplot_File << Conn_Hexa_Par[iNode+4] << "\t" << Conn_Hexa_Par[iNode+5] << "\t";
        Tecplot_File << Conn_Hexa_Par[iNode+6] << "\t" << Conn_Hexa_Par[iNode+7] << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Pris; iElem++) {
        iNode = iElem*N_POINTS_PRISM;
        Tecplot_File << Conn_Pris_Par[iNode+0] << "\t" << Conn_Pris_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Pris_Par[iNode+1] << "\t" << Conn_Pris_Par[iNode+2] << "\t";
        Tecplot_File << Conn_Pris_Par[iNode+3] << "\t" << Conn_Pris_Par[iNode+4] << "\t";
        Tecplot_File << Conn_Pris_Par[iNode+4] << "\t" << Conn_Pris_Par[iNode+5] << "\n";
      }
      
      for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
        iNode = iElem*N_POINTS_PYRAMID;
        Tecplot_File << Conn_Pyra_Par[iNode+0] << "\t" << Conn_Pyra_Par[iNode+1] << "\t";
        Tecplot_File << Conn_Pyra_Par[iNode+2] << "\t" << Conn_Pyra_Par[iNode+3] << "\t";
        Tecplot_File << Conn_Pyra_Par[iNode+4] << "\t" << Conn_Pyra_Par[iNode+4] << "\t";
        Tecplot_File << Conn_Pyra_Par[iNode+4] << "\t" << Conn_Pyra_Par[iNode+4] << "\n";
      }
      }
      
    }
    Tecplot_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }
  
  Tecplot_File.close();
  
}

#ifdef HAVE_MPI

namespace
{

/*!
 * \brief Calculate the partitioning of nodes to determine:
 * (a) For a given global node number, to which partition does it belong and what is its local node number; and
 * (b) How many nodes are held by each partition.
 */
class NodePartitioner {
public:
  /*!
   * \param[in] global_num_nodes - The total number of nodes being output
   * \param[in] num_ranks - The number of MPI ranks involved in the output
   */
  NodePartitioner(unsigned long global_num_nodes, int num_ranks)
    : m_num_ranks(num_ranks) {
    /* rank i has (1-based) global nodes m_node_range[i] + 1 through m_node_range[i + 1] */
    unsigned long nodes_per_rank = global_num_nodes / num_ranks;
    unsigned long num_extra_nodes = global_num_nodes - nodes_per_rank * num_ranks;
    m_node_range.resize(num_ranks + 1);
    m_node_range[0] = 0;
    for (int ii = 1; ii <= num_ranks; ii++) {
      m_node_range[ii] = m_node_range[ii - 1] + nodes_per_rank;
      if (num_extra_nodes > 0) {
        ++m_node_range[ii];
        --num_extra_nodes;
      }
    }
    assert(m_node_range[num_ranks] == global_num_nodes);
  }

  /*!
   * \brief Determine the MPI rank that owns a global node number and its corresponding local node number.
   * \param global_node_number[in] - The global node number; global node numbers are sequential across all MPI ranks.
   * \param owning_rank[out] - The MPI rank that owns (will output) the global node
   * \param node_number[out] - The rank-local node number for the given global node number
   */
  void GetOwningRankAndNodeNumber(unsigned long global_node_number, int &owning_rank, unsigned long &node_number)
  {
    owning_rank = static_cast<int>(global_node_number / m_node_range[1]);
    if (owning_rank >= m_num_ranks)
      owning_rank = m_num_ranks - 1;
    while(global_node_number > m_node_range[owning_rank + 1])
      ++owning_rank;
    while(global_node_number <= m_node_range[owning_rank])
      --owning_rank;
    node_number = global_node_number - m_node_range[owning_rank];
  }

  /*!
   * \brief Determine the number of nodes to be output by a particular rank
   * \param which_rank[in] - The MPI rank
   * \ret - The number of nodes that will be output by the give MPI rank.
   */
  int64_t GetRankNumNodes(int which_rank)
  {
    return static_cast<int64_t>(m_node_range[which_rank + 1] - m_node_range[which_rank]);
  }

private:
  int m_num_ranks;
  vector<unsigned long> m_node_range;
};

int64_t GetHaloNodeNumber(unsigned long global_node_number, unsigned long last_local_node, vector<unsigned long> const &halo_node_list)
{
  vector<unsigned long>::const_iterator it = lower_bound(halo_node_list.begin(), halo_node_list.end(), global_node_number);
  assert(it != halo_node_list.end());
  assert(*it == global_node_number);
  /* When C++11 is universally available, replace the following mouthful with "auto" */
  iterator_traits<vector<unsigned long>::const_iterator>::difference_type offset = distance(halo_node_list.begin(), it);
  assert(offset >= 0);
  return (int64_t)(last_local_node + offset + 1);
}

} /* namespace */

#endif /* HAVE_MPI */

void COutput::WriteTecplotBinary_Parallel(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol) {

#ifdef HAVE_TECIO
  
  /*--- Open Tecplot binary file. ---*/
  
  string filename = GetTecplotFilename(config, val_iZone, val_nZone, surf_sol, ".szplt");
  
  string data_set_title = surf_sol
    ? "Visualization of the surface solution"
    : "Visualization of the volumetric solution";

  ostringstream tecplot_variable_names;
  for (size_t iVar = 0; iVar < Variable_Names.size()-1; ++iVar) {
    tecplot_variable_names << Variable_Names[iVar] << ",";
  }
  tecplot_variable_names << Variable_Names[Variable_Names.size()-1];

  void* file_handle = NULL;
  int32_t err = tecFileWriterOpen(filename.c_str(), data_set_title.c_str(), tecplot_variable_names.str().c_str(),
    FILEFORMAT_SZL, FILETYPE_FULL, (int32_t)FieldDataType_Double, NULL, &file_handle);
  if (err) cout << "Error opening Tecplot file '" << filename << "'" << endl;

#ifdef HAVE_MPI
  err = tecMPIInitialize(file_handle, MPI_COMM_WORLD, MASTER_NODE);
  if (err) cout << "Error initializing Tecplot parallel output." << endl;
#endif
  
  /*--- Define the zone(s). For 2D, and for 3D surfaces, each rank outputs a separate zone. ---*/

  int64_t num_nodes;
  int64_t num_cells;
  int32_t zone_type;
  if (surf_sol) {
    num_nodes = static_cast<int64_t>(nGlobal_Surf_Poin);
    num_cells = static_cast<int64_t>(nSurf_Elem_Par);
    if (geometry->GetnDim() == 2) 
      zone_type = ZONETYPE_FELINESEG;
    else
      zone_type = ZONETYPE_FEQUADRILATERAL;
  } else {
    num_nodes = static_cast<int64_t>(nGlobal_Poin_Par);
    num_cells = static_cast<int64_t>(nGlobal_Elem_Par);
    if (geometry->GetnDim() == 2)
      zone_type = ZONETYPE_FEQUADRILATERAL;
    else
      zone_type = ZONETYPE_FEBRICK;
  }

  bool is_unsteady = false;
  passivedouble solution_time = 0.0;
  if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    is_unsteady = true;
    solution_time = SU2_TYPE::GetValue(config->GetDelta_UnstTime()*config->GetExtIter());
  } else if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    is_unsteady = true;
    /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
    passivedouble period = SU2_TYPE::GetValue(config->GetHarmonicBalance_Period());
    passivedouble deltaT = period/SU2_TYPE::GetValue(config->GetnTimeInstances());
    solution_time = deltaT*val_iZone;
  }

  int32_t zone;
  vector<int32_t> value_locations(nVar_Par, 1); /* Nodal variables. */
  err = tecZoneCreateFE(file_handle, "Zone", zone_type, num_nodes, num_cells, NULL, NULL, &value_locations[0], NULL, 0, 0, 0, &zone);
  if (err) cout << rank << ": Error creating Tecplot zone." << endl;
  if (is_unsteady) {
    err = tecZoneSetUnsteadyOptions(file_handle, zone, solution_time, config->GetExtIter() + 1);
    if (err) cout << rank << ": Error setting Tecplot zone unsteady options." << std::endl;
  }

#ifdef HAVE_MPI

  unsigned short iVar;
  NodePartitioner node_partitioner(num_nodes, size);
  set<unsigned long> halo_nodes;
  vector<unsigned long> sorted_halo_nodes;
  vector<passivedouble> halo_var_data;
  vector<int> num_nodes_to_receive(size, 0);
  vector<int> values_to_receive_displacements(size);

  if (zone_type == ZONETYPE_FEBRICK) {

    /* We output a single, partitioned zone where each rank outputs one partition. */
    vector<int32_t> partition_owners;
    partition_owners.reserve(size);
    for (int32_t iRank = 0; iRank < size; ++iRank)
      partition_owners.push_back(iRank);
    err = tecZoneMapPartitionsToMPIRanks(file_handle, zone, size, &partition_owners[0]);
    if (err) cout << rank << ": Error assigning MPI ranks for Tecplot zone partitions." << endl;
  
    /* Gather a list of nodes we refer to but are not outputting. */

    for (unsigned long i = 0; i < nParallel_Tria * N_POINTS_TRIANGLE; ++i)
      if ((unsigned long)Conn_Tria_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Tria_Par[i])
        halo_nodes.insert(Conn_Tria_Par[i]);
  
    for (unsigned long i = 0; i < nParallel_Quad * N_POINTS_QUADRILATERAL; ++i)
      if ((unsigned long)Conn_Quad_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Quad_Par[i])
        halo_nodes.insert(Conn_Quad_Par[i]);
 
    for (unsigned long i = 0; i < nParallel_Tetr * N_POINTS_TETRAHEDRON; ++i)
      if ((unsigned long)Conn_Tetr_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Tetr_Par[i])
        halo_nodes.insert(Conn_Tetr_Par[i]);

    for (unsigned long i = 0; i < nParallel_Hexa * N_POINTS_HEXAHEDRON; ++i)
      if ((unsigned long)Conn_Hexa_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Hexa_Par[i])
        halo_nodes.insert(Conn_Hexa_Par[i]);
      
    for (unsigned long i = 0; i < nParallel_Pris * N_POINTS_PRISM; ++i)
      if ((unsigned long)Conn_Pris_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Pris_Par[i])
        halo_nodes.insert(Conn_Pris_Par[i]);
    
    for (unsigned long i = 0; i < nParallel_Pyra * N_POINTS_PYRAMID; ++i)
      if ((unsigned long)Conn_Pyra_Par[i] <= beg_node[rank] || end_node[rank] < (unsigned long)Conn_Pyra_Par[i])
        halo_nodes.insert(Conn_Pyra_Par[i]);

    /* Sorted list of halo nodes for this MPI rank. */
    sorted_halo_nodes.assign(halo_nodes.begin(), halo_nodes.end());
        
    /* Have to include all nodes our cells refer to or TecIO will barf, so add the halo node count to the number of local nodes. */
    int64_t partition_num_nodes = end_node[rank] - beg_node[rank] + static_cast<int64_t>(halo_nodes.size());
    int64_t partition_num_cells = nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;

    /*--- We effectively tack the halo nodes onto the end of the node list for this partition.
      TecIO will later replace them with references to nodes in neighboring partitions. */
    size_t num_halo_nodes = sorted_halo_nodes.size();
    vector<int64_t> halo_node_local_numbers(max((size_t)1, num_halo_nodes)); /* Min size 1 to avoid crashes when we access these vectors below. */
    vector<int32_t> neighbor_partitions(max((size_t)1, num_halo_nodes));
    vector<int64_t> neighbor_nodes(max((size_t)1, num_halo_nodes));
    for(int64_t i = 0; i < static_cast<int64_t>(num_halo_nodes); ++i) {
      halo_node_local_numbers[i] = end_node[rank] - beg_node[rank] + i + 1;
      int owning_rank;
      unsigned long node_number;
      node_partitioner.GetOwningRankAndNodeNumber(sorted_halo_nodes[i], owning_rank, node_number);
      neighbor_partitions[i] = owning_rank + 1; /* Partition numbers are 1-based. */
      neighbor_nodes[i] = static_cast<int64_t>(node_number);
    }
    err = tecFEPartitionCreate64(file_handle, zone, rank + 1, partition_num_nodes, partition_num_cells,
      static_cast<int64_t>(num_halo_nodes), &halo_node_local_numbers[0], &neighbor_partitions[0], &neighbor_nodes[0], 0, NULL);
    if (err) cout << rank << ": Error creating Tecplot zone partition." << endl;

    /* Gather halo node data. First, tell each rank how many nodes' worth of data we need from them. */
    for (size_t i = 0; i < num_halo_nodes; ++i)
      ++num_nodes_to_receive[neighbor_partitions[i] - 1];
    vector<int> num_nodes_to_send(size);
    SU2_MPI::Alltoall(&num_nodes_to_receive[0], 1, MPI_INT, &num_nodes_to_send[0], 1, MPI_INT, MPI_COMM_WORLD);

    /* Now send the global node numbers whose data we need,
       and receive the same from all other ranks.
       Each rank has globally consecutive node numbers,
       so we can just parcel out sorted_halo_nodes for send. */
    vector<int> nodes_to_send_displacements(size);
    vector<int> nodes_to_receive_displacements(size);
    nodes_to_send_displacements[0] = 0;
    nodes_to_receive_displacements[0] = 0;
    for(int iRank = 1; iRank < size; ++iRank) {
      nodes_to_send_displacements[iRank] = nodes_to_send_displacements[iRank - 1] + num_nodes_to_send[iRank - 1];
      nodes_to_receive_displacements[iRank] = nodes_to_receive_displacements[iRank - 1] + num_nodes_to_receive[iRank - 1];
    }
    int total_num_nodes_to_send = nodes_to_send_displacements[size - 1] + num_nodes_to_send[size - 1];
    vector<unsigned long> nodes_to_send(max(1, total_num_nodes_to_send));

    /* The terminology gets a bit confusing here. We're sending the node numbers
       (sorted_halo_nodes) whose data we need to receive, and receiving
       lists of nodes whose data we need to send. */
    if (sorted_halo_nodes.empty()) sorted_halo_nodes.resize(1); /* Avoid crash. */
    SU2_MPI::Alltoallv(&sorted_halo_nodes[0], &num_nodes_to_receive[0], &nodes_to_receive_displacements[0], MPI_UNSIGNED_LONG,
                       &nodes_to_send[0],     &num_nodes_to_send[0],    &nodes_to_send_displacements[0],    MPI_UNSIGNED_LONG,
                       MPI_COMM_WORLD);
    
    /* Now actually send and receive the data */
    vector<passivedouble> data_to_send(max(1, total_num_nodes_to_send * nVar_Par));
    halo_var_data.resize(max((size_t)1, nVar_Par * num_halo_nodes));
    vector<int> num_values_to_send(size);
    vector<int> values_to_send_displacements(size);
    vector<int> num_values_to_receive(size);
    size_t index = 0;
    for(int iRank = 0; iRank < size; ++iRank) {
      /* We send and receive nVar_Par values per node. */
      num_values_to_send[iRank]              = num_nodes_to_send[iRank] * nVar_Par;
      values_to_send_displacements[iRank]    = nodes_to_send_displacements[iRank] * nVar_Par;
      num_values_to_receive[iRank]           = num_nodes_to_receive[iRank] * nVar_Par;
      values_to_receive_displacements[iRank] = nodes_to_receive_displacements[iRank] * nVar_Par;
      for(iVar = 0; iVar < nVar_Par; ++iVar)
        for(int iNode = 0; iNode < num_nodes_to_send[iRank]; ++iNode) {
          unsigned long node_offset = nodes_to_send[nodes_to_send_displacements[iRank] + iNode] - beg_node[rank] - 1;
          data_to_send[index++] = SU2_TYPE::GetValue(Parallel_Data[iVar][node_offset]);
        }
    }
    SU2_MPI::Alltoallv(&data_to_send[0],  &num_values_to_send[0],    &values_to_send_displacements[0],    MPI_DOUBLE,
                       &halo_var_data[0], &num_values_to_receive[0], &values_to_receive_displacements[0], MPI_DOUBLE,
                       MPI_COMM_WORLD);
  }
  else {
    /* Zone will be gathered to and output by MASTER_NODE */
    int32_t partition_owner = MASTER_NODE;
    err = tecZoneMapPartitionsToMPIRanks(file_handle, zone, 1, &partition_owner);
  }

  /*--- Write surface and volumetric solution data. ---*/
  
  if (zone_type == ZONETYPE_FEBRICK) {
    std::vector<passivedouble> values_to_write(nParallel_Poin);
    for (iVar = 0; err == 0 && iVar < nVar_Par; iVar++) {
      for(unsigned long i = 0; i < nParallel_Poin; ++i)
        values_to_write[i] = SU2_TYPE::GetValue(Parallel_Data[iVar][i]);
      err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, rank + 1, nParallel_Poin, &values_to_write[0]);
      if (err) cout << rank << ": Error outputting Tecplot variable values." << endl;
      for (int iRank = 0; err == 0 && iRank < size; ++iRank) {
        if (num_nodes_to_receive[iRank] > 0) {
          int var_data_offset = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank] * iVar;
          err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, rank + 1, static_cast<int64_t>(num_nodes_to_receive[iRank]), &halo_var_data[var_data_offset]);
          if (err) cout << rank << ": Error outputting Tecplot halo values." << endl;
        }
      }
    }
  } else {
    if (rank == MASTER_NODE) {
      vector<passivedouble> var_data;
      vector<unsigned long> num_surface_points(size);
      if (surf_sol)
        SU2_MPI::Gather(&nSurf_Poin_Par, 1, MPI_UNSIGNED_LONG, &num_surface_points[0], 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      for(int iRank = 0; iRank < size; ++iRank) {
        int64_t rank_num_points;
        if (surf_sol)
          rank_num_points = num_surface_points[iRank];
        else
          rank_num_points = node_partitioner.GetRankNumNodes(iRank);
        if (rank_num_points > 0) {
          if (iRank == rank) { /* Output local data. */
            std::vector<passivedouble> values_to_write;
            for (iVar = 0; err == 0 && iVar < nVar_Par; iVar++) {
              if (surf_sol) {
                values_to_write.resize(nSurf_Poin_Par);
                for(unsigned long i = 0; i < nSurf_Poin_Par; ++i)
                  values_to_write[i] = SU2_TYPE::GetValue(Parallel_Surf_Data[iVar][i]);
                err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, nSurf_Poin_Par, &values_to_write[0]);
              }
              else {
                values_to_write.resize(rank_num_points);
                for(unsigned long i = 0; i < (unsigned long)rank_num_points; ++i)
                  values_to_write[i] = SU2_TYPE::GetValue(Parallel_Data[iVar][i]);
                err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, rank_num_points, &values_to_write[0]);
              }
              if (err) cout << rank << ": Error outputting Tecplot variable values." << endl;
            }
          }
          else { /* Receive data from other rank. */
            var_data.resize(max((int64_t)1, nVar_Par * rank_num_points));
            SU2_MPI::Recv(&var_data[0], nVar_Par * rank_num_points, MPI_DOUBLE, iRank, iRank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (iVar = 0; err == 0 && iVar < nVar_Par; iVar++) {
              err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, rank_num_points, &var_data[iVar * rank_num_points]);
              if (err) cout << rank << ": Error outputting Tecplot surface variable values." << endl;
            }
          }
        }
      }
    }
    else { /* Send data to MASTER_NODE */
      if (surf_sol)
        SU2_MPI::Gather(&nSurf_Poin_Par, 1, MPI_UNSIGNED_LONG, NULL, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      vector<passivedouble> var_data;
      size_t var_data_size = nVar_Par * (surf_sol ? nSurf_Poin_Par : nParallel_Poin);
      var_data.reserve(var_data_size);
      for (iVar = 0; err == 0 && iVar < nVar_Par; iVar++)
        if (surf_sol)
          for(unsigned long i = 0; i < nSurf_Poin_Par; ++i)
            var_data.push_back(SU2_TYPE::GetValue(Parallel_Surf_Data[iVar][i]));
        else
          for(unsigned long i = 0; i < nParallel_Poin; ++i)
            var_data.push_back(SU2_TYPE::GetValue(Parallel_Data[iVar][i]));
      if (var_data.size() > 0)
        SU2_MPI::Send(&var_data[0], static_cast<int>(var_data.size()), MPI_DOUBLE, MASTER_NODE, rank, MPI_COMM_WORLD);
    }
  }

#else

  unsigned short iVar;

  vector<passivedouble> var_data;
  size_t var_data_size = nVar_Par * (surf_sol ? nSurf_Poin_Par : nParallel_Poin);
  var_data.reserve(var_data_size);
  
  if (surf_sol) {
    for (iVar = 0; err == 0 && iVar < nVar_Par; iVar++) {
      for(unsigned long i = 0; i < nSurf_Poin_Par; ++i)
        var_data.push_back(SU2_TYPE::GetValue(Parallel_Surf_Data[iVar][i]));
      err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, nSurf_Poin_Par, &var_data[iVar * nSurf_Poin_Par]);
      if (err) cout << rank << ": Error outputting Tecplot variable value." << endl;
    }
  } else {
    for (iVar = 0; err == 0 && iVar < nVar_Par; iVar++) {
      for(unsigned long i = 0; i < nParallel_Poin; ++i)
        var_data.push_back(SU2_TYPE::GetValue(Parallel_Data[iVar][i]));
      err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, nParallel_Poin, &var_data[iVar * nParallel_Poin]);
      if (err) cout << rank << ": Error outputting Tecplot variable value." << endl;
    }
  }

#endif /* HAVE_MPI */
  
  /*--- Write connectivity data. ---*/

  unsigned long iElem, iNode;
  
#ifdef HAVE_MPI
  if (zone_type == ZONETYPE_FEBRICK) {

    int64_t nodes[8];

    /**
     *  Each rank writes node numbers relative to the partition it is outputting (starting with node number 1).
     *  Ghost (halo) nodes identified above are numbered sequentially just beyond the end of the actual, local nodes.
     *  Note that beg_node and end_node refer to zero-based node numbering, but Conn_* contain one-based node numbers.
     */
#define MAKE_LOCAL(n) beg_node[rank] < (unsigned long)n && (unsigned long)n <= end_node[rank] \
  ? (int64_t)((unsigned long)n - beg_node[rank]) \
  : GetHaloNodeNumber(n, end_node[rank] - beg_node[rank], sorted_halo_nodes)

    for (iElem = 0; err == 0 && iElem < nParallel_Tetr; iElem++) {
      iNode = iElem*N_POINTS_TETRAHEDRON;
      nodes[0] = MAKE_LOCAL(Conn_Tetr_Par[iNode+0]);
      nodes[1] = MAKE_LOCAL(Conn_Tetr_Par[iNode+1]);
      nodes[2] = MAKE_LOCAL(Conn_Tetr_Par[iNode+2]);
      nodes[3] = nodes[2];
      nodes[4] = MAKE_LOCAL(Conn_Tetr_Par[iNode+3]);
      nodes[5] = nodes[4];
      nodes[6] = nodes[4];
      nodes[7] = nodes[4];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Hexa; iElem++) {
      iNode = iElem*N_POINTS_HEXAHEDRON;
      nodes[0] = MAKE_LOCAL(Conn_Hexa_Par[iNode+0]);
      nodes[1] = MAKE_LOCAL(Conn_Hexa_Par[iNode+1]);
      nodes[2] = MAKE_LOCAL(Conn_Hexa_Par[iNode+2]);
      nodes[3] = MAKE_LOCAL(Conn_Hexa_Par[iNode+3]);
      nodes[4] = MAKE_LOCAL(Conn_Hexa_Par[iNode+4]);
      nodes[5] = MAKE_LOCAL(Conn_Hexa_Par[iNode+5]);
      nodes[6] = MAKE_LOCAL(Conn_Hexa_Par[iNode+6]);
      nodes[7] = MAKE_LOCAL(Conn_Hexa_Par[iNode+7]);
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
      
    for (iElem = 0; err == 0 && iElem < nParallel_Pris; iElem++) {
      iNode = iElem*N_POINTS_PRISM;
      nodes[0] = MAKE_LOCAL(Conn_Pris_Par[iNode+0]);
      nodes[1] = MAKE_LOCAL(Conn_Pris_Par[iNode+1]);
      nodes[2] = nodes[1];
      nodes[3] = MAKE_LOCAL(Conn_Pris_Par[iNode+2]);
      nodes[4] = MAKE_LOCAL(Conn_Pris_Par[iNode+3]);
      nodes[5] = MAKE_LOCAL(Conn_Pris_Par[iNode+4]);
      nodes[6] = nodes[5];
      nodes[7] = MAKE_LOCAL(Conn_Pris_Par[iNode+5]);
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
    
    for (iElem = 0; err == 0 && iElem < nParallel_Pyra; iElem++) {
      iNode = iElem*N_POINTS_PYRAMID;
      nodes[0] = MAKE_LOCAL(Conn_Pyra_Par[iNode+0]);
      nodes[1] = MAKE_LOCAL(Conn_Pyra_Par[iNode+1]);
      nodes[2] = MAKE_LOCAL(Conn_Pyra_Par[iNode+2]);
      nodes[3] = MAKE_LOCAL(Conn_Pyra_Par[iNode+3]);
      nodes[4] = MAKE_LOCAL(Conn_Pyra_Par[iNode+4]);
      nodes[5] = nodes[4];
      nodes[6] = nodes[4];
      nodes[7] = nodes[4];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
  } else {
    if (rank == MASTER_NODE) {

      /* Non-hexahedral output by the master node. Output local data directly, and gather other data from the other ranks. */

      int64_t nodes[4];

      vector<unsigned long> connectivity_sizes(size);
      unsigned long unused = 0;
      SU2_MPI::Gather(&unused, 1, MPI_UNSIGNED_LONG, &connectivity_sizes[0], 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      vector<int64_t> connectivity;
      for(int iRank = 0; iRank < size; ++iRank) {
        if (iRank == rank) {
          if (surf_sol) {
            for (iElem = 0; err == 0 && iElem < nParallel_Line; iElem++) {
              iNode = iElem*N_POINTS_LINE;
              nodes[0] = Conn_BoundLine_Par[iNode+0];
              nodes[1] = Conn_BoundLine_Par[iNode+1];
              err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 2, nodes);
              if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
            }
          
            for (iElem = 0; err == 0 && iElem < nParallel_BoundTria; iElem++) {
              iNode = iElem*N_POINTS_TRIANGLE;
              nodes[0] = Conn_BoundTria_Par[iNode+0];
              nodes[1] = Conn_BoundTria_Par[iNode+1];
              nodes[2] = Conn_BoundTria_Par[iNode+2];
              nodes[3] = Conn_BoundTria_Par[iNode+2];
              err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 4, nodes);
              if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
            }
          
            for (iElem = 0; err == 0 && iElem < nParallel_BoundQuad; iElem++) {
              iNode = iElem*N_POINTS_QUADRILATERAL;
              nodes[0] = Conn_BoundQuad_Par[iNode+0];
              nodes[1] = Conn_BoundQuad_Par[iNode+1];
              nodes[2] = Conn_BoundQuad_Par[iNode+2];
              nodes[3] = Conn_BoundQuad_Par[iNode+3];
              err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 4, nodes);
              if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
            }
          } else {
            for (iElem = 0; err == 0 && iElem < nParallel_Tria; iElem++) {
              iNode = iElem*N_POINTS_TRIANGLE;
              nodes[0] = Conn_Tria_Par[iNode+0];
              nodes[1] = Conn_Tria_Par[iNode+1];
              nodes[2] = Conn_Tria_Par[iNode+2];
              nodes[3] = Conn_Tria_Par[iNode+2];
              err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 4, nodes);
              if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
            }
        
            for (iElem = 0; err == 0 && iElem < nParallel_Quad; iElem++) {
              iNode = iElem*N_POINTS_QUADRILATERAL;
              nodes[0] = Conn_Quad_Par[iNode+0];
              nodes[1] = Conn_Quad_Par[iNode+1];
              nodes[2] = Conn_Quad_Par[iNode+2];
              nodes[3] = Conn_Quad_Par[iNode+3];
              err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 4, nodes);
              if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
            }
          }
        } else { /* Receive node map and write out. */
          connectivity.resize(max((unsigned long)1, connectivity_sizes[iRank]));
          SU2_MPI::Recv(&connectivity[0], connectivity_sizes[iRank], MPI_UNSIGNED_LONG, iRank, iRank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, connectivity_sizes[iRank], &connectivity[0]);
          if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
        }
      }
    } else {

      /* Non-hexahedral output by non-master node. Send what we've got to the master node. */

      unsigned long connectivity_size;
      if (surf_sol)
        connectivity_size = 2 * nParallel_Line + 4 * nParallel_BoundTria + 4 * nParallel_BoundQuad;
      else
        connectivity_size = 4 * (nParallel_Tria + nParallel_Quad);
      SU2_MPI::Gather(&connectivity_size, 1, MPI_UNSIGNED_LONG, NULL, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      vector<int64_t> connectivity;
      connectivity.reserve(connectivity_size);
      if (surf_sol) {
        for (iElem = 0; err == 0 && iElem < nParallel_Line; iElem++) {
          iNode = iElem*N_POINTS_LINE;
          connectivity.push_back(Conn_BoundLine_Par[iNode+0]);
          connectivity.push_back(Conn_BoundLine_Par[iNode+1]);
        }
      
        for (iElem = 0; err == 0 && iElem < nParallel_BoundTria; iElem++) {
          iNode = iElem*N_POINTS_TRIANGLE;
          connectivity.push_back(Conn_BoundTria_Par[iNode+0]);
          connectivity.push_back(Conn_BoundTria_Par[iNode+1]);
          connectivity.push_back(Conn_BoundTria_Par[iNode+2]);
          connectivity.push_back(Conn_BoundTria_Par[iNode+2]);
        }
      
        for (iElem = 0; err == 0 && iElem < nParallel_BoundQuad; iElem++) {
          iNode = iElem*N_POINTS_QUADRILATERAL;
          connectivity.push_back(Conn_BoundQuad_Par[iNode+0]);
          connectivity.push_back(Conn_BoundQuad_Par[iNode+1]);
          connectivity.push_back(Conn_BoundQuad_Par[iNode+2]);
          connectivity.push_back(Conn_BoundQuad_Par[iNode+3]);
        }
      } else {
        for (iElem = 0; err == 0 && iElem < nParallel_Tria; iElem++) {
          iNode = iElem*N_POINTS_TRIANGLE;
          connectivity.push_back(Conn_Tria_Par[iNode+0]);
          connectivity.push_back(Conn_Tria_Par[iNode+1]);
          connectivity.push_back(Conn_Tria_Par[iNode+2]);
          connectivity.push_back(Conn_Tria_Par[iNode+2]);
        }
    
        for (iElem = 0; err == 0 && iElem < nParallel_Quad; iElem++) {
          iNode = iElem*N_POINTS_QUADRILATERAL;
          connectivity.push_back(Conn_Quad_Par[iNode+0]);
          connectivity.push_back(Conn_Quad_Par[iNode+1]);
          connectivity.push_back(Conn_Quad_Par[iNode+2]);
          connectivity.push_back(Conn_Quad_Par[iNode+3]);
        }
      }
      if (connectivity.empty()) connectivity.resize(1); /* Avoid crash */
      SU2_MPI::Send(&connectivity[0], connectivity_size, MPI_UNSIGNED_LONG, MASTER_NODE, rank, MPI_COMM_WORLD);
    }
  }
#else
  if (surf_sol) {

    int64_t nodes[4];

    for (iElem = 0; err == 0 && iElem < nParallel_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      nodes[0] = Conn_BoundLine_Par[iNode+0];
      nodes[1] = Conn_BoundLine_Par[iNode+1];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 2, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
        
    for (iElem = 0; err == 0 && iElem < nParallel_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      nodes[0] = Conn_BoundTria_Par[iNode+0];
      nodes[1] = Conn_BoundTria_Par[iNode+1];
      nodes[2] = Conn_BoundTria_Par[iNode+2];
      nodes[3] = Conn_BoundTria_Par[iNode+2];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 4, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
        
    for (iElem = 0; err == 0 && iElem < nParallel_BoundQuad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      nodes[0] = Conn_BoundQuad_Par[iNode+0];
      nodes[1] = Conn_BoundQuad_Par[iNode+1];
      nodes[2] = Conn_BoundQuad_Par[iNode+2];
      nodes[3] = Conn_BoundQuad_Par[iNode+3];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 4, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }

  } else {

    int64_t nodes[8];

    for (iElem = 0; err == 0 && iElem < nParallel_Tria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      nodes[0] = Conn_Tria_Par[iNode+0];
      nodes[1] = Conn_Tria_Par[iNode+1];
      nodes[2] = Conn_Tria_Par[iNode+2];
      nodes[3] = Conn_Tria_Par[iNode+2];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 4, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Quad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      nodes[0] = Conn_Quad_Par[iNode+0];
      nodes[1] = Conn_Quad_Par[iNode+1];
      nodes[2] = Conn_Quad_Par[iNode+2];
      nodes[3] = Conn_Quad_Par[iNode+3];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 4, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Tetr; iElem++) {
      iNode = iElem*N_POINTS_TETRAHEDRON;
      nodes[0] = Conn_Tetr_Par[iNode+0];
      nodes[1] = Conn_Tetr_Par[iNode+1];
      nodes[2] = Conn_Tetr_Par[iNode+2];
      nodes[3] = Conn_Tetr_Par[iNode+2];
      nodes[4] = Conn_Tetr_Par[iNode+3];
      nodes[5] = Conn_Tetr_Par[iNode+3];
      nodes[6] = Conn_Tetr_Par[iNode+3];
      nodes[7] = Conn_Tetr_Par[iNode+3];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Hexa; iElem++) {
      iNode = iElem*N_POINTS_HEXAHEDRON;
      nodes[0] = Conn_Hexa_Par[iNode+0];
      nodes[1] = Conn_Hexa_Par[iNode+1];
      nodes[2] = Conn_Hexa_Par[iNode+2];
      nodes[3] = Conn_Hexa_Par[iNode+3];
      nodes[4] = Conn_Hexa_Par[iNode+4];
      nodes[5] = Conn_Hexa_Par[iNode+5];
      nodes[6] = Conn_Hexa_Par[iNode+6];
      nodes[7] = Conn_Hexa_Par[iNode+7];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
      
    for (iElem = 0; err == 0 && iElem < nParallel_Pris; iElem++) {
      iNode = iElem*N_POINTS_PRISM;
      nodes[0] = Conn_Pris_Par[iNode+0];
      nodes[1] = Conn_Pris_Par[iNode+1];
      nodes[2] = Conn_Pris_Par[iNode+1];
      nodes[3] = Conn_Pris_Par[iNode+2];
      nodes[4] = Conn_Pris_Par[iNode+3];
      nodes[5] = Conn_Pris_Par[iNode+4];
      nodes[6] = Conn_Pris_Par[iNode+4];
      nodes[7] = Conn_Pris_Par[iNode+5];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
    
    for (iElem = 0; err == 0 && iElem < nParallel_Pyra; iElem++) {
      iNode = iElem*N_POINTS_PYRAMID;
      nodes[0] = Conn_Pyra_Par[iNode+0];
      nodes[1] = Conn_Pyra_Par[iNode+1];
      nodes[2] = Conn_Pyra_Par[iNode+2];
      nodes[3] = Conn_Pyra_Par[iNode+3];
      nodes[4] = Conn_Pyra_Par[iNode+4];
      nodes[5] = Conn_Pyra_Par[iNode+4];
      nodes[6] = Conn_Pyra_Par[iNode+4];
      nodes[7] = Conn_Pyra_Par[iNode+4];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
      
  }

#endif
  
  err = tecFileWriterClose(&file_handle);
  if (err) cout << rank << ": Error finishing Tecplot file output." << endl;
  
#endif /* HAVE_TECIO */

}

void COutput::SetTecplotBinary_DomainMesh(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
#ifdef HAVE_TECIO
  
  passivedouble   t;
  INTEGER4 i, err, Debug, NPts, NElm, N2DElm, NVolElm, IsDouble, KMax;
  INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, StrandID, ParentZn, FileFormat, FileType;
  INTEGER4 *ShareFromZone = NULL, IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
  string buffer, variables;
  stringstream file;
  bool first_zone = true;
  bool adjoint = config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint();
  unsigned short dims = geometry->GetnDim();
  enum     FileFormat { PLT = 0, SZPLT = 1 };
  enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };
  enum   ZoneType { ORDERED=0, FELINESEG=1, FETRIANGLE=2, FEQUADRILATERAL=3, FETETRAHEDRON=4, FEBRICK=5, FEPOLYGON=6, FEPOLYHEDRON=7 };
  
  /*--- Consistent data for Tecplot zones ---*/
  
  Debug            = 0;
  IsDouble          = 1;
  NPts            = (INTEGER4)nGlobal_Poin;
  t              = 0.0;//iExtIter*config->GetDelta_UnstTimeND();
  KMax            = 0;
  ICellMax          = 0;
  JCellMax          = 0;
  KCellMax          = 0;
  StrandID          = 0;//(INTEGER4)iExtIter;
  ParentZn          = 0;
  IsBlock            = 1;
  NumFaceConnections      = 0;
  FaceNeighborMode      = 0;
  ShareConnectivityFromZone  = 0;
  
  /*--- Write Tecplot solution file ---*/
  
  if (!wrote_base_file) {
    
    file.str(string());

    if (adjoint)
      buffer = config->GetAdj_FileName();
    else buffer = config->GetFlow_FileName();

    if (config->GetKind_SU2() == SU2_DOT) {
      buffer = config->GetVolSens_FileName();
    }
    
    file << buffer << ".mesh.szplt";
    FileFormat = SZPLT;
    FileType = GRID;
    
    if (dims == 2) variables = "x y";
    else if (dims == 3) variables = "x y z";
    else cout << "Error: wrong number of dimensions: " << dims << endl;
    
    /*--- Open Tecplot file ---*/
    err = TECINI142((char *)config->GetFlow_FileName().c_str(),
                    (char *)variables.c_str(),
                    (char *)file.str().c_str(),
                    (char *)".",
                    &FileFormat,
                    &FileType,
                    &Debug,
                    &IsDouble);
    if (err) cout << "Error in opening Tecplot file" << endl;
    
    first_zone = true;
    
    N2DElm = (INTEGER4)(nGlobal_Tria + nGlobal_Quad);
    if (N2DElm > 0) {
      
      /*--- Write the zone header information ---*/
      
      if ((INTEGER4)nGlobal_Tria < N2DElm) {   /* Create a Quad zone with a mixed element types */

      
        ZoneType = FEQUADRILATERAL; NElm = N2DElm;
      
        err = TECZNE142((char*)"Mixed Elements",
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


      } else {   /* Create a Tria zone */

        ZoneType = FETRIANGLE; NElm = (INTEGER4)nGlobal_Tria;
      
        err = TECZNE142((char*)"Triangle Elements",
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

      }
      
      /*--- write node coordinates and data if not done already---*/
      
      if (first_zone) {
        
        ShareFromZone = new INTEGER4[dims];
        for (i = 0; i < dims; i++) ShareFromZone[i] = 0;
        
        if (config->GetKind_SU2() == SU2_SOL) {
          err = TECDAT142(&NPts, Data[0], &IsDouble); ShareFromZone[0] = 1;
          err = TECDAT142(&NPts, Data[1], &IsDouble); ShareFromZone[1] = 1;
          if (geometry->GetnDim() == 3) {
            err = TECDAT142(&NPts, Data[2], &IsDouble);
            ShareFromZone[2] = 1;
          }
          
        } else if (config->GetKind_SU2() == SU2_DOT) {
          
          passivedouble* PassiveData = new passivedouble[NPts];
          
          for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[0][i]);
          err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[0] = 1;
          
          for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[1][i]);
          err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[1] = 1;
          
          if (geometry->GetnDim() == 3) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[2][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble);
            ShareFromZone[2] = 1;
          }
          
          delete [] PassiveData;
          
        } else {
          err = TECDAT142(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
          err = TECDAT142(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
          if (geometry->GetnDim() == 3) {
            err = TECDAT142(&NPts, Coords[2], &IsDouble);
            ShareFromZone[2] = 1;
          }
        }
        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
        first_zone = false;
      }
      

      if (nGlobal_Tria > 0) {
        if ((INTEGER4)nGlobal_Tria < N2DElm) {   /* Write Tria connecivity as collapsed Quad */
    
          /*--- Convert the triangle connectivity from 3 nodes to 4 nodes for FEQUADRALATERL */
          int *Conn_Tria_Mod = new int[nGlobal_Tria*N_POINTS_QUADRILATERAL];
          unsigned long iNode_Tria, iNode_Quad;
          for (unsigned long iElem = 0; iElem < nGlobal_Tria; iElem++) {
            iNode_Tria = iElem*N_POINTS_TRIANGLE;
            iNode_Quad = iElem*N_POINTS_QUADRILATERAL;
            Conn_Tria_Mod[iNode_Quad+0] = Conn_Tria[iNode_Tria+0];
            Conn_Tria_Mod[iNode_Quad+1] = Conn_Tria[iNode_Tria+1];
            Conn_Tria_Mod[iNode_Quad+2] = Conn_Tria[iNode_Tria+2];
            Conn_Tria_Mod[iNode_Quad+3] = Conn_Tria[iNode_Tria+2];
          }
          NElm = (INTEGER4)(nGlobal_Tria*N_POINTS_QUADRILATERAL);
          err = TECNODE142(&NElm, Conn_Tria_Mod);
          if (err) cout << "Error writing triangle connectivity to Tecplot file" << endl;
          delete [] Conn_Tria_Mod;
  
        } else {   /* Write Tria connectivity */

          err = TECNOD142(Conn_Tria);
          if (err) cout << "Error writing connectivity to Tecplot file" << endl;
        }
      }

      if (nGlobal_Quad > 0) {
      
        NElm = (INTEGER4)(nGlobal_Quad*N_POINTS_QUADRILATERAL);
        err = TECNODE142(&NElm, Conn_Quad);
        if (err) cout << "Error writing connectivity to Tecplot file" << endl;

      }
      
    }

    /*--- Create 3D Volume Zone ---*/
    NVolElm = (INTEGER4)(nGlobal_Tetr + nGlobal_Pyra + nGlobal_Pris + nGlobal_Hexa);
    if (NVolElm > 0) {
      
      /*--- Write the zone header information ---*/
      
      if ((INTEGER4)nGlobal_Tetr < NVolElm) {   /* Create a Hexa zone with a mixed element types */
      
        /*--- Write the mixed-element zone header information ---*/
      
        ZoneType = FEBRICK; NElm = NVolElm;
      
        err = TECZNE142((char*)"Mixed Elements",
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
      
      } else {   /* Create a Tetra zone */

        /*--- Write the tetrahedral zone header information ---*/

        ZoneType = FETETRAHEDRON; NElm = (INTEGER4)nGlobal_Tetr;
        
        err = TECZNE142((char*)"Tetrahedral Elements",
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
      }
      
      /*--- write node coordinates and data if not done already---*/
      
      if (first_zone) {
        
        ShareFromZone = new INTEGER4[dims];
        for (i = 0; i < dims; i++) ShareFromZone[i] = 0;
        
        if (config->GetKind_SU2() == SU2_SOL) {
          err = TECDAT142(&NPts, Data[0], &IsDouble); ShareFromZone[0] = 1;
          err = TECDAT142(&NPts, Data[1], &IsDouble); ShareFromZone[1] = 1;
          if (geometry->GetnDim() == 3) {
            err = TECDAT142(&NPts, Data[2], &IsDouble);
            ShareFromZone[2] = 1;
          }
          
        } else if (config->GetKind_SU2() == SU2_DOT) {
          
          passivedouble* PassiveData = new passivedouble[NPts];
          
          for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[0][i]);
          err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[0] = 1;
          
          for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[1][i]);
          err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[1] = 1;
          
          if (geometry->GetnDim() == 3) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[2][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble);
            ShareFromZone[2] = 1;
          }
          
          delete [] PassiveData;
          
        }  else {
          err = TECDAT142(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
          err = TECDAT142(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
          if (geometry->GetnDim() == 3) {
            err = TECDAT142(&NPts, Coords[2], &IsDouble);
            ShareFromZone[2] = 1;
          }
        }
        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
        first_zone = false;
      }

    }
      
    if (nGlobal_Tetr > 0) {
      if ((INTEGER4)nGlobal_Tetr < NVolElm) {   /* Write tetra connecivity as collapsed hexa */
    
        /*--- Convert the tetrahedral connectivity from 4 nodes to 8 nodes for FEBRICK ---*/
        int *Conn_Tetr_Mod = new int[nGlobal_Tetr*N_POINTS_HEXAHEDRON];
        unsigned long iNode_Tetr, iNode_Hexa;
        for (unsigned long iElem = 0; iElem < nGlobal_Tetr; iElem++) {
          iNode_Tetr = iElem*N_POINTS_TETRAHEDRON;
          iNode_Hexa = iElem*N_POINTS_HEXAHEDRON;
          Conn_Tetr_Mod[iNode_Hexa+0] = Conn_Tetr[iNode_Tetr+0];
          Conn_Tetr_Mod[iNode_Hexa+1] = Conn_Tetr[iNode_Tetr+1];
          Conn_Tetr_Mod[iNode_Hexa+2] = Conn_Tetr[iNode_Tetr+2];
          Conn_Tetr_Mod[iNode_Hexa+3] = Conn_Tetr[iNode_Tetr+2];
          Conn_Tetr_Mod[iNode_Hexa+4] = Conn_Tetr[iNode_Tetr+3];
          Conn_Tetr_Mod[iNode_Hexa+5] = Conn_Tetr[iNode_Tetr+3];
          Conn_Tetr_Mod[iNode_Hexa+6] = Conn_Tetr[iNode_Tetr+3];
          Conn_Tetr_Mod[iNode_Hexa+7] = Conn_Tetr[iNode_Tetr+3];
        }
        NElm = (INTEGER4)(nGlobal_Tetr*N_POINTS_HEXAHEDRON);
        err = TECNODE142(&NElm, Conn_Tetr_Mod);
        if (err) cout << "Error writing tetrahedral connectivity to Tecplot file" << endl;
        delete [] Conn_Tetr_Mod;
  
      } else {   /* Write Tetra connectivity */

        err = TECNOD142(Conn_Tetr);
        if (err) cout << "Error writing connectivity to Tecplot file" << endl;
      }

    }

    if (nGlobal_Hexa > 0) {
      
      NElm = (INTEGER4)(nGlobal_Hexa*N_POINTS_HEXAHEDRON);
      err = TECNODE142(&NElm, Conn_Hexa);
      if (err) cout << "Error writing connectivity to Tecplot file" << endl;
      
    }
    
    if (nGlobal_Pyra > 0) {
      
      /*--- Convert the pyramid connectivity from 5 nodes to 8 nodes for FEBRICK ---*/
      int *Conn_Pyra_Mod = new int[nGlobal_Pyra*N_POINTS_HEXAHEDRON];
      unsigned long iNode_Pyra, iNode_Hexa;
      for (unsigned long iElem = 0; iElem < nGlobal_Pyra; iElem++) {
        iNode_Pyra = iElem*N_POINTS_PYRAMID;
        iNode_Hexa = iElem*N_POINTS_HEXAHEDRON;
        Conn_Pyra_Mod[iNode_Hexa+0] = Conn_Pyra[iNode_Pyra+4];
        Conn_Pyra_Mod[iNode_Hexa+1] = Conn_Pyra[iNode_Pyra+4];
        Conn_Pyra_Mod[iNode_Hexa+2] = Conn_Pyra[iNode_Pyra+4];
        Conn_Pyra_Mod[iNode_Hexa+3] = Conn_Pyra[iNode_Pyra+4];
        Conn_Pyra_Mod[iNode_Hexa+4] = Conn_Pyra[iNode_Pyra+0];
        Conn_Pyra_Mod[iNode_Hexa+5] = Conn_Pyra[iNode_Pyra+1];
        Conn_Pyra_Mod[iNode_Hexa+6] = Conn_Pyra[iNode_Pyra+2];
        Conn_Pyra_Mod[iNode_Hexa+7] = Conn_Pyra[iNode_Pyra+3];
      }
      NElm = (INTEGER4)(nGlobal_Pyra*N_POINTS_HEXAHEDRON);
      err = TECNODE142(&NElm, Conn_Pyra_Mod);
      if (err) cout << "Error writing pyramid connectivity to Tecplot file" << endl;
      delete [] Conn_Pyra_Mod;
      
    }
    
    if (nGlobal_Pris > 0) {
      
      /*--- Convert the prism connectivity from 6 nodes to 8 nodes for FEBRICK ---*/
      int *Conn_Pris_Mod = new int[nGlobal_Pris*N_POINTS_HEXAHEDRON];
      unsigned long iNode_Pris, iNode_Hexa;
      for (unsigned long iElem = 0; iElem < nGlobal_Pris; iElem++) {
        iNode_Pris = iElem*N_POINTS_PRISM;
        iNode_Hexa = iElem*N_POINTS_HEXAHEDRON;
        Conn_Pris_Mod[iNode_Hexa+0] = Conn_Pris[iNode_Pris+0];
        Conn_Pris_Mod[iNode_Hexa+1] = Conn_Pris[iNode_Pris+0];
        Conn_Pris_Mod[iNode_Hexa+2] = Conn_Pris[iNode_Pris+1];
        Conn_Pris_Mod[iNode_Hexa+3] = Conn_Pris[iNode_Pris+2];
        Conn_Pris_Mod[iNode_Hexa+4] = Conn_Pris[iNode_Pris+3];
        Conn_Pris_Mod[iNode_Hexa+5] = Conn_Pris[iNode_Pris+3];
        Conn_Pris_Mod[iNode_Hexa+6] = Conn_Pris[iNode_Pris+4];
        Conn_Pris_Mod[iNode_Hexa+7] = Conn_Pris[iNode_Pris+5];
      }
      NElm = (INTEGER4)(nGlobal_Pris*N_POINTS_HEXAHEDRON);
      err = TECNODE142(&NElm, Conn_Pris_Mod);
      if (err) cout << "Error writing prism connectivity to Tecplot file" << endl;
      delete [] Conn_Pris_Mod;
      
    }
    
    delete [] ShareFromZone;
    wrote_base_file = true;
    
    err = TECEND142();
    if (err) cout << "Error in closing Tecplot file" << endl;
    
  }
  
#endif
  
}

void COutput::SetTecplotBinary_DomainSolution(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
#ifdef HAVE_TECIO
  
  passivedouble   t;
  INTEGER4 i, iVar, err, Debug, NPts, NElm, N2DElm, NVolElm, IsDouble, KMax;
  INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, StrandID, ParentZn, FileFormat, FileType;
  INTEGER4 *ShareFromZone = NULL, IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
  string buffer, variables;
  stringstream file;
  bool first_zone = true, unsteady = config->GetUnsteady_Simulation(), GridMovement = config->GetGrid_Movement();
  bool Wrt_Unsteady = config->GetWrt_Unsteady();
  unsigned long iExtIter = config->GetExtIter();
  unsigned short NVar, dims = geometry->GetnDim();
  enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };
  enum     FileFormat { PLT = 0, SZPLT = 1 };
  enum   ZoneType { ORDERED=0, FELINESEG=1, FETRIANGLE=2, FEQUADRILATERAL=3, FETETRAHEDRON=4, FEBRICK=5, FEPOLYGON=6, FEPOLYHEDRON=7 };
  
  bool adjoint = config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint();
  unsigned short Kind_Solver = config->GetKind_Solver();

  /*--- Consistent data for Tecplot zones ---*/
  Debug            = 0;
  IsDouble          = 1;
  NPts            = (INTEGER4)nGlobal_Poin;
  t              = SU2_TYPE::GetValue(iExtIter*config->GetDelta_UnstTime());
  KMax            = 0;
  ICellMax          = 0;
  JCellMax          = 0;
  KCellMax          = 0;
  StrandID          = (INTEGER4)iExtIter+1;
  ParentZn          = 0;
  IsBlock            = 1;
  NumFaceConnections      = 0;
  FaceNeighborMode      = 0;
  ShareConnectivityFromZone  = 0;
  
  file.str(string());

  /*--- Write file name with extension ---*/
  
  if (adjoint)
    buffer = config->GetAdj_FileName();
  else buffer = config->GetFlow_FileName();
  
  if (Kind_Solver == FEM_ELASTICITY) {
    buffer = config->GetStructure_FileName().c_str();
  }
  
  if (Kind_Solver == HEAT_EQUATION_FVM) {
    buffer = config->GetHeat_FileName().c_str();
  }
  
  if (config->GetKind_SU2() == SU2_DOT) {
    buffer = config->GetVolSens_FileName();
  }
  file << buffer;
  
  if (unsteady) {
    if (((int)iExtIter >= 0) && ((int)iExtIter < 10))      file << "_0000" << iExtIter;
    if (((int)iExtIter >= 10) && ((int)iExtIter < 100))    file << "_000" << iExtIter;
    if (((int)iExtIter >= 100) && ((int)iExtIter < 1000))    file << "_00" << iExtIter;
    if (((int)iExtIter >= 1000) && ((int)iExtIter < 10000))  file << "_0" << iExtIter;
    if ((int)iExtIter >= 10000)              file << iExtIter;
  }

  file << ".sol.szplt";
  FileFormat = SZPLT;
  FileType = SOLUTION;
  variables = AssembleVariableNames(geometry, config, nVar_Consv, &NVar);
  if ((config->GetKind_SU2() == SU2_SOL) || (config->GetKind_SU2() == SU2_DOT)) {
    if (Wrt_Unsteady && GridMovement) nVar_Total = NVar;
    else nVar_Total = NVar+dims;
  }

  /*--- Open Tecplot file ---*/
  err = TECINI142((char *)config->GetFlow_FileName().c_str(),
                  (char *)variables.c_str(),
                  (char *)file.str().c_str(),
                  (char *)".",
                  &FileFormat,
                  &FileType,
                  &Debug,
                  &IsDouble);
  if (err) cout << "Error in opening Tecplot file" << endl;
  
  first_zone = true;
  ShareFromZone = new INTEGER4[NVar];
  for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;
  
  N2DElm = (INTEGER4)(nGlobal_Tria + nGlobal_Quad);
  if (N2DElm > 0) {
      
    /*--- Write the zone header information ---*/
      
    if ((INTEGER4)nGlobal_Tria < N2DElm) {   /* Create a Quad zone with a mixed element types */
      
      ZoneType = FEQUADRILATERAL; NElm = N2DElm;
      
      err = TECZNE142((char*)"Mixed Elements",
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
      
      
    } else {   /* Create a Tria zone */

      ZoneType = FETRIANGLE; NElm = (INTEGER4)nGlobal_Tria;
    
      err = TECZNE142((char*)"Triangle Elements",
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
    }
    
    /*--- write node coordinates and data if not done already---*/
    if (first_zone) {
      
      ShareFromZone = new INTEGER4[NVar];
      for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;

      i = 0;
      if (config->GetKind_SU2() == SU2_SOL) {
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
      } else if (config->GetKind_SU2() == SU2_DOT) {
        
        passivedouble* PassiveData = new passivedouble[NPts];
        
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
      
        delete [] PassiveData;
        
      } else {

        if (Wrt_Unsteady && GridMovement) {

          err = TECDAT142(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          err = TECDAT142(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          if (dims == 3) {
            err = TECDAT142(&NPts, Coords[2], &IsDouble);
            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
            ShareFromZone[i++] = 1;
          }
        }

        for (iVar = 0; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      }
      first_zone = false;
    }
    
  }
  if (nGlobal_Quad > 0) {
    
    /*--- write node coordinates and data if not done already---*/
    if (first_zone) {
      
      ShareFromZone = new INTEGER4[NVar];
      for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;
      
      i = 0;
      if (config->GetKind_SU2() == SU2_SOL) {
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
      } else if (config->GetKind_SU2() == SU2_DOT) {
        
        passivedouble* PassiveData = new passivedouble[NPts];
        
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
        
        delete [] PassiveData;
        
      } else {
        if (Wrt_Unsteady && GridMovement) {
          err = TECDAT142(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          err = TECDAT142(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          if (dims == 3) {
            err = TECDAT142(&NPts, Coords[2], &IsDouble);
            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
            ShareFromZone[i++] = 1;
          }
        }
        for (iVar = 0; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      }
      
      first_zone = false;
    }
    
  }

  /*--- Create 3D Volume Zone ---*/
  NVolElm = (INTEGER4)(nGlobal_Tetr + nGlobal_Pyra + nGlobal_Pris + nGlobal_Hexa);
  if (NVolElm > 0) {
      
    /*--- Write the zone header information ---*/
      
    if ((INTEGER4)nGlobal_Tetr < NVolElm) {   /* Create a Hexa zone with a mixed element types */
      
      /*--- Write the mixed-element zone header information ---*/
      
      ZoneType = FEBRICK; NElm = NVolElm;
      
      err = TECZNE142((char*)"Mixed Elements",
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
      
    } else {   /* Create a Tetra zone */

      /*--- Write the tetrahedral zone header information ---*/

      ZoneType = FETETRAHEDRON; NElm = (INTEGER4)nGlobal_Tetr;
    
      err = TECZNE142((char*)"Tetrahedral Elements",
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
    }
    
    /*--- write node coordinates and data if not done already---*/
    if (first_zone) {
      
      ShareFromZone = new INTEGER4[NVar];
      for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;
      
      i = 0;
      if (config->GetKind_SU2() == SU2_SOL) {
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
      } else if (config->GetKind_SU2() == SU2_DOT) {
        
        passivedouble* PassiveData = new passivedouble[NPts];
        
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
        
        delete [] PassiveData;
        
      } else {
        if (Wrt_Unsteady && GridMovement) {
          err = TECDAT142(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          err = TECDAT142(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          if (dims == 3) {
            err = TECDAT142(&NPts, Coords[2], &IsDouble);
            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
            ShareFromZone[i++] = 1;
          }
        }
        for (iVar = 0; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      }
      
      first_zone = false;
    }

  }

  if (nGlobal_Hexa > 0) {
    
    /*--- write node coordinates and data if not done already---*/
    if (first_zone) {
      
      ShareFromZone = new INTEGER4[NVar];
      for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;
      
      i = 0;
      if (config->GetKind_SU2() == SU2_SOL) {
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
      } else if (config->GetKind_SU2() == SU2_DOT) {
        
        passivedouble* PassiveData = new passivedouble[NPts];
        
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
        
        delete [] PassiveData;
        
      } else {
        if (Wrt_Unsteady && GridMovement) {
          err = TECDAT142(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          err = TECDAT142(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          if (dims == 3) {
            err = TECDAT142(&NPts, Coords[2], &IsDouble);
            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
            ShareFromZone[i++] = 1;
          }
        }
        for (iVar = 0; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      }
      
      first_zone = false;
    }
    
  }
  if (nGlobal_Pyra > 0) {
    
    /*--- write node coordinates and data if not done already---*/
    if (first_zone) {
      
      ShareFromZone = new INTEGER4[NVar];
      for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;
      
      i = 0;
      if (config->GetKind_SU2() == SU2_SOL) {
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
      } else if (config->GetKind_SU2() == SU2_DOT) {
        
        passivedouble* PassiveData = new passivedouble[NPts];
        
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
        
        delete [] PassiveData;
        
      } else {
        if (Wrt_Unsteady && GridMovement) {
          err = TECDAT142(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          err = TECDAT142(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          if (dims == 3) {
            err = TECDAT142(&NPts, Coords[2], &IsDouble);
            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
            ShareFromZone[i++] = 1;
          }
        }
        for (iVar = 0; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      }
      
      first_zone = false;
    }
    
  }
  if (nGlobal_Pris > 0) {

    /*--- write node coordinates and data if not done already---*/
    if (first_zone) {
      
      ShareFromZone = new INTEGER4[NVar];
      for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;
      
      i = 0;
      if (config->GetKind_SU2() == SU2_SOL) {
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
      } else if (config->GetKind_SU2() == SU2_DOT) {
        
        passivedouble* PassiveData = new passivedouble[NPts];
        
        if (Wrt_Unsteady && GridMovement) {
          for (iVar = 0; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        } else {
          for (iVar = dims; iVar < nVar_Total; iVar++) {
            for (i = 0; i < NPts; i++) PassiveData[i] = SU2_TYPE::GetValue(Data[iVar][i]);
            err = TECDAT142(&NPts, PassiveData, &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
        
        delete [] PassiveData;
        
      } else {
        if (Wrt_Unsteady && GridMovement) {
          err = TECDAT142(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          err = TECDAT142(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing coordinates to Tecplot file" << endl;
          if (dims == 3) {
            err = TECDAT142(&NPts, Coords[2], &IsDouble);
            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
            ShareFromZone[i++] = 1;
          }
        }
        for (iVar = 0; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      }
      
      first_zone = false;
    }
  }
  
  delete [] ShareFromZone;
  
  err = TECEND142();
  if (err) cout << "Error in closing Tecplot file" << endl;
  
#endif
  
}

void COutput::SetTecplotBinary_SurfaceMesh(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
#ifdef HAVE_TECIO
  
  passivedouble   t;
  INTEGER4 i, err, Debug, NPts, NElm, IsDouble, KMax;
  INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, StrandID, ParentZn, FileFormat, FileType;
  INTEGER4 *ShareFromZone, IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
  string buffer, variables;
  stringstream file;
  bool first_zone = true;
  unsigned short iDim, dims = geometry->GetnDim();
  unsigned long iPoint, iElem, iNode;
  enum     FileFormat { PLT = 0, SZPLT = 1 };
  enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };
  enum   ZoneType { ORDERED=0, FELINESEG=1, FETRIANGLE=2, FEQUADRILATERAL=3, FETETRAHEDRON=4, FEBRICK=5, FEPOLYGON=6, FEPOLYHEDRON=7 };
  
  /*--- Write Tecplot solution file ---*/
  if (!wrote_surf_file) {
    
    file.str(string());
    buffer = config->GetSurfFlowCoeff_FileName();
    if (config->GetKind_SU2() == SU2_DOT) {
      buffer = config->GetSurfSens_FileName();
    }
    
    FileFormat = PLT;

    file << buffer << ".mesh.szplt";
    FileType = GRID;
    
    if (dims == 2) variables = "x y";
    else if (dims == 3) variables = "x y z";
    else cout << "Error: wrong number of dimensions: " << dims << endl;
    
    first_zone = true;
    ShareFromZone = new INTEGER4[dims];
    for (i = 0; i < dims; i++) ShareFromZone[i] = 0;
    
    /*--- Perform a renumbering for the surface points/elements ---*/
    unsigned long *LocalIndex = new unsigned long [nGlobal_Poin+1];
    bool *SurfacePoint = new bool [nGlobal_Poin+1];
    
    for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) SurfacePoint[iPoint] = false;
    
    for (iElem = 0; iElem < nGlobal_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      SurfacePoint[Conn_Line[iNode+0]] = true;
      SurfacePoint[Conn_Line[iNode+1]] = true;
    }
    for (iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      SurfacePoint[Conn_BoundTria[iNode+0]] = true;
      SurfacePoint[Conn_BoundTria[iNode+1]] = true;
      SurfacePoint[Conn_BoundTria[iNode+2]] = true;
    }
    for (iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      SurfacePoint[Conn_BoundQuad[iNode+0]] = true;
      SurfacePoint[Conn_BoundQuad[iNode+1]] = true;
      SurfacePoint[Conn_BoundQuad[iNode+2]] = true;
      SurfacePoint[Conn_BoundQuad[iNode+3]] = true;
    }
    
    unsigned long nSurf_Poin = 0;
    for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) {
      LocalIndex[iPoint] = 0;
      if (SurfacePoint[iPoint]) { nSurf_Poin++; LocalIndex[iPoint] = nSurf_Poin; }
    }
    
    /*--- Collect surface coordinates into one array as well ---*/
    /*--- Note the -1 in the Coords/Data array in order to undo the 1-based indexing ---*/
    passivedouble **Surf_Coords = new passivedouble*[dims];
    for (iDim = 0; iDim < dims; iDim++)
      Surf_Coords[iDim] = new passivedouble[nSurf_Poin];
    
    unsigned long iSurf_Poin = 0;
    for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) {
      if (SurfacePoint[iPoint]) {
        for (iDim = 0; iDim < dims; iDim++) {
          if ((config->GetKind_SU2() == SU2_SOL) || (config->GetKind_SU2() == SU2_DOT))
            Surf_Coords[iDim][iSurf_Poin] = SU2_TYPE::GetValue(Data[iDim][iPoint-1]);
          else
            Surf_Coords[iDim][iSurf_Poin] = SU2_TYPE::GetValue(Coords[iDim][iPoint-1]);
        }
        iSurf_Poin++;
      }
    }
    
    /*--- Consistent data for Tecplot zones ---*/
    Debug            = 0;
    IsDouble          = 1;
    NPts            = (INTEGER4)nSurf_Poin;
    t              = 0.0;//iExtIter*config->GetDelta_UnstTimeND();
    KMax            = 0;
    ICellMax          = 0;
    JCellMax          = 0;
    KCellMax          = 0;
    StrandID          = 0;//(INTEGER4)iExtIter;
    ParentZn          = 0;
    IsBlock            = 1;
    NumFaceConnections      = 0;
    FaceNeighborMode      = 0;
    ShareConnectivityFromZone  = 0;
    
    /*--- Open Tecplot file ---*/
    err = TECINI142((char *)config->GetSurfFlowCoeff_FileName().c_str(),
                    (char *)variables.c_str(),
                    (char *)file.str().c_str(),
                    (char *)".",
                    &FileFormat,
                    &FileType,
                    &Debug,
                    &IsDouble);
    if (err) cout << "Error in opening Tecplot file" << endl;
    
    
    if (nGlobal_Line > 0) {
      
      /*--- Put the connectivity into a single array for writing ---*/
      int *Conn_Line_New = new int[nGlobal_Line*N_POINTS_LINE];
      iNode = 0;
      for (iElem = 0; iElem < nGlobal_Line; iElem++) {
        iNode = iElem*N_POINTS_LINE;
        Conn_Line_New[iNode+0] = LocalIndex[Conn_Line[iNode+0]];
        Conn_Line_New[iNode+1] = LocalIndex[Conn_Line[iNode+1]];
      }
      
      /*--- Write the zone header information ---*/
      ZoneType = FELINESEG; NElm = (INTEGER4)nGlobal_Line;
      
      err = TECZNE142((char*)"Line Elements",
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
        
        err = TECDAT142(&NPts, Surf_Coords[0], &IsDouble); ShareFromZone[0] = 1;
        err = TECDAT142(&NPts, Surf_Coords[1], &IsDouble); ShareFromZone[1] = 1;
        if (geometry->GetnDim() == 3) {
          err = TECDAT142(&NPts, Surf_Coords[2], &IsDouble);
          ShareFromZone[2] = 1;
        }
        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
        first_zone = false;
      }
      
      err = TECNOD142(Conn_Line_New);
      if (err) cout << "Error writing connectivity to Tecplot file" << endl;
      
      delete [] Conn_Line_New;
    }
    
    if (nGlobal_BoundTria > 0) {
      
      /*--- Put the connectivity into a single array for writing ---*/
      int *Conn_BoundTria_New = new int[nGlobal_BoundTria*N_POINTS_TRIANGLE];
      
      iNode = 0;
      for (iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
        iNode = iElem*N_POINTS_TRIANGLE;
        Conn_BoundTria_New[iNode+0] = LocalIndex[Conn_BoundTria[iNode+0]];
        Conn_BoundTria_New[iNode+1] = LocalIndex[Conn_BoundTria[iNode+1]];
        Conn_BoundTria_New[iNode+2] = LocalIndex[Conn_BoundTria[iNode+2]];
      }
      
      /*--- Write the zone header information ---*/
      ZoneType = FETRIANGLE; NElm = (INTEGER4)nGlobal_BoundTria;
      
      err = TECZNE142((char*)"Triangle Elements",
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
        
        err = TECDAT142(&NPts, Surf_Coords[0], &IsDouble); ShareFromZone[0] = 1;
        err = TECDAT142(&NPts, Surf_Coords[1], &IsDouble); ShareFromZone[1] = 1;
        if (geometry->GetnDim() == 3) {
          err = TECDAT142(&NPts, Surf_Coords[2], &IsDouble);
          ShareFromZone[2] = 1;
        }
        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
        first_zone = false;
      }
      
      err = TECNOD142(Conn_BoundTria_New);
      if (err) cout << "Error writing connectivity to Tecplot file" << endl;
      
      delete [] Conn_BoundTria_New;
    }
    
    if (nGlobal_BoundQuad > 0) {
      
      
      /*--- Put the connectivity into a single array for writing ---*/
      int *Conn_BoundQuad_New = new int[nGlobal_BoundQuad*N_POINTS_QUADRILATERAL];
      iNode = 0;
      for (iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
        iNode = iElem*N_POINTS_QUADRILATERAL;
        Conn_BoundQuad_New[iNode+0] = LocalIndex[Conn_BoundQuad[iNode+0]];
        Conn_BoundQuad_New[iNode+1] = LocalIndex[Conn_BoundQuad[iNode+1]];
        Conn_BoundQuad_New[iNode+2] = LocalIndex[Conn_BoundQuad[iNode+2]];
        Conn_BoundQuad_New[iNode+3] = LocalIndex[Conn_BoundQuad[iNode+3]];
      }
      
      /*--- Write the zone header information ---*/
      ZoneType = FEQUADRILATERAL; NElm = (INTEGER4)nGlobal_BoundQuad;
      
      err = TECZNE142((char*)"Quadrilateral Elements",
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
        
        err = TECDAT142(&NPts, Surf_Coords[0], &IsDouble); ShareFromZone[0] = 1;
        err = TECDAT142(&NPts, Surf_Coords[1], &IsDouble); ShareFromZone[1] = 1;
        if (geometry->GetnDim() == 3) {
          err = TECDAT142(&NPts, Surf_Coords[2], &IsDouble);
          ShareFromZone[2] = 1;
        }
        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
        first_zone = false;
      }
      
      err = TECNOD142(Conn_BoundQuad_New);
      if (err) cout << "Error writing connectivity to Tecplot file" << endl;
      
      delete [] Conn_BoundQuad_New;
    }
    
    for (iDim = 0; iDim < dims; iDim++)
      delete [] Surf_Coords[iDim];
    delete [] Surf_Coords;
    delete [] ShareFromZone;
    delete [] LocalIndex;
    delete [] SurfacePoint;
    
    wrote_surf_file = true;
    
    err = TECEND142();
    if (err) cout << "Error in closing Tecplot file" << endl;
    
  }
  
#endif
  
}

void COutput::SetTecplotBinary_SurfaceSolution(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
#ifdef HAVE_TECIO
  
  passivedouble   t;
  INTEGER4 i, iVar, err, Debug, NPts, NElm, IsDouble, KMax;
  INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, StrandID, ParentZn, FileFormat, FileType;
  INTEGER4 *ShareFromZone, IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
  string buffer, variables;
  stringstream file;
  bool first_zone = true, unsteady = config->GetUnsteady_Simulation(), GridMovement = config->GetGrid_Movement();
  bool Wrt_Unsteady = config->GetWrt_Unsteady();
  unsigned long iPoint, iElem, iNode, iSurf_Poin, iExtIter = config->GetExtIter();
  unsigned short iDim, NVar, dims = geometry->GetnDim();
  enum     FileFormat { PLT = 0, SZPLT = 1 };
  enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };
  enum   ZoneType { ORDERED=0, FELINESEG=1, FETRIANGLE=2, FEQUADRILATERAL=3, FETETRAHEDRON=4, FEBRICK=5, FEPOLYGON=6, FEPOLYHEDRON=7 };
  
  bool adjoint = config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint();
  unsigned short Kind_Solver = config->GetKind_Solver();
  
  file.str(string());
  
  /*--- Write file name with extension ---*/
  
  if (adjoint) buffer = config->GetSurfAdjCoeff_FileName();
  else buffer = config->GetSurfFlowCoeff_FileName();
  
  if (Kind_Solver == FEM_ELASTICITY) {
    buffer = config->GetSurfStructure_FileName().c_str();
  }
  
  if (Kind_Solver == HEAT_EQUATION_FVM) {
    buffer = config->GetSurfHeat_FileName().c_str();
  }
  
  if (config->GetKind_SU2() == SU2_DOT) {
    buffer = config->GetSurfSens_FileName();
  }
  
  file << buffer;
  
  if (unsteady) {
    if (((int)iExtIter >= 0) && ((int)iExtIter < 10))      file << "_0000" << iExtIter;
    if (((int)iExtIter >= 10) && ((int)iExtIter < 100))    file << "_000" << iExtIter;
    if (((int)iExtIter >= 100) && ((int)iExtIter < 1000))    file << "_00" << iExtIter;
    if (((int)iExtIter >= 1000) && ((int)iExtIter < 10000))  file << "_0" << iExtIter;
    if ((int)iExtIter >= 10000)              file << iExtIter;
  }
  file << ".sol.szplt";
  FileFormat = PLT;
  FileType = SOLUTION;
  variables = AssembleVariableNames(geometry, config, nVar_Consv, &NVar);
  if ((config->GetKind_SU2() == SU2_SOL) || (config->GetKind_SU2() == SU2_DOT)) {
    if (Wrt_Unsteady && GridMovement) nVar_Total = NVar;
    else nVar_Total = NVar+dims;
  }
  
  first_zone = true;
  ShareFromZone = new INTEGER4[NVar];
  for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;
  
  
  /*--- Perform a renumbering for the surface points/elements ---*/
  unsigned long *LocalIndex = new unsigned long [nGlobal_Poin+1];
  bool *SurfacePoint = new bool [nGlobal_Poin+1];
  
  for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) SurfacePoint[iPoint] = false;
  
  for (iElem = 0; iElem < nGlobal_Line; iElem++) {
    iNode = iElem*N_POINTS_LINE;
    SurfacePoint[Conn_Line[iNode+0]] = true;
    SurfacePoint[Conn_Line[iNode+1]] = true;
  }
  for (iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
    iNode = iElem*N_POINTS_TRIANGLE;
    SurfacePoint[Conn_BoundTria[iNode+0]] = true;
    SurfacePoint[Conn_BoundTria[iNode+1]] = true;
    SurfacePoint[Conn_BoundTria[iNode+2]] = true;
  }
  for (iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
    iNode = iElem*N_POINTS_QUADRILATERAL;
    SurfacePoint[Conn_BoundQuad[iNode+0]] = true;
    SurfacePoint[Conn_BoundQuad[iNode+1]] = true;
    SurfacePoint[Conn_BoundQuad[iNode+2]] = true;
    SurfacePoint[Conn_BoundQuad[iNode+3]] = true;
  }
  
  unsigned long nSurf_Poin = 0;
  for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) {
    LocalIndex[iPoint] = 0;
    if (SurfacePoint[iPoint]) { nSurf_Poin++; LocalIndex[iPoint] = nSurf_Poin; }
  }
  
  /*--- Collect surface coordinates into one array as well ---*/
  /*--- Note the -1 in the Coords/Data array in order to undo the 1-based indexing ---*/
  passivedouble **Surf_Coords = NULL;
  if (Wrt_Unsteady && GridMovement) {
    Surf_Coords = new passivedouble*[dims];
    for (iDim = 0; iDim < dims; iDim++)
    Surf_Coords[iDim] = new passivedouble[nSurf_Poin];
    
    iSurf_Poin = 0;
    for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) {
      if (SurfacePoint[iPoint]) {
        for (iDim = 0; iDim < dims; iDim++) {
          if ((config->GetKind_SU2() == SU2_SOL) || (config->GetKind_SU2() == SU2_DOT))
            Surf_Coords[iDim][iSurf_Poin] = SU2_TYPE::GetValue(Data[iDim][iPoint-1]);
          else
            Surf_Coords[iDim][iSurf_Poin] = SU2_TYPE::GetValue(Coords[iDim][iPoint-1]);
        }
        iSurf_Poin++;
      }
    }
  }
  
  /*--- Collect surface data into one array for the surface as well ---*/
  /*--- Note the -1 in the Coords/Data array in order to undo the 1-based indexing ---*/
  passivedouble **Surf_Data = new passivedouble*[nVar_Total];
  for (iVar = 0; iVar < nVar_Total; iVar++)
  Surf_Data[iVar] = new passivedouble[nSurf_Poin];
  
  iSurf_Poin = 0;
  for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) {
    if (SurfacePoint[iPoint]) {
      for (iVar = 0; iVar < nVar_Total; iVar++) {
        if ((config->GetKind_SU2() == SU2_SOL) || (config->GetKind_SU2() == SU2_DOT)) {
          if (Wrt_Unsteady && GridMovement)
            Surf_Data[iVar][iSurf_Poin] = SU2_TYPE::GetValue(Data[iVar][iPoint-1]);
          else
          Surf_Data[iVar][iSurf_Poin] = SU2_TYPE::GetValue(Data[iVar][iPoint-1]);
        } else
        Surf_Data[iVar][iSurf_Poin] = SU2_TYPE::GetValue(Data[iVar][iPoint-1]);
      }
      iSurf_Poin++;
    }
  }
  
  /*--- Consistent data for Tecplot zones ---*/
  Debug            = 0;
  IsDouble          = 1;
  NPts            = (INTEGER4)nSurf_Poin;
  t              = SU2_TYPE::GetValue(iExtIter*config->GetDelta_UnstTime());
  KMax            = 0;
  ICellMax          = 0;
  JCellMax          = 0;
  KCellMax          = 0;
  StrandID          = (INTEGER4)iExtIter+1;
  ParentZn          = 0;
  IsBlock            = 1;
  NumFaceConnections      = 0;
  FaceNeighborMode      = 0;
  ShareConnectivityFromZone  = 0;
  
  
  /*--- Open Tecplot file ---*/
  err = TECINI142((char *)config->GetFlow_FileName().c_str(),
                  (char *)variables.c_str(),
                  (char *)file.str().c_str(),
                  (char *)".",
                  &FileFormat,
                  &FileType,
                  &Debug,
                  &IsDouble);
  if (err) cout << "Error in opening Tecplot file" << endl;
  
  
  if (nGlobal_Line > 0) {
    
    /*--- Write the zone header information ---*/
    ZoneType = FELINESEG; NElm = (INTEGER4)nGlobal_Line;
    
    err = TECZNE142((char*)"Line Elements",
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
      if ((config->GetKind_SU2() == SU2_SOL) || (config->GetKind_SU2() == SU2_DOT)) {
        if (Wrt_Unsteady && GridMovement) {
          for (iDim = 0; iDim < dims; iDim++) {
            err = TECDAT142(&NPts, Surf_Data[iDim], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
        for (iVar = dims; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      } else {
        if (Wrt_Unsteady && GridMovement) {
          for (iDim = 0; iDim < dims; iDim++) {
            err = TECDAT142(&NPts, Surf_Coords[iDim], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
        for (iVar = 0; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      }
      first_zone = false;
    }
    
  }
  
  if (nGlobal_BoundTria > 0) {
    
    /*--- Write the zone header information ---*/
    ZoneType = FETRIANGLE; NElm = (INTEGER4)nGlobal_BoundTria;
    
    err = TECZNE142((char*)"Triangle Elements",
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
      if ((config->GetKind_SU2() == SU2_SOL) || (config->GetKind_SU2() == SU2_DOT)) {
        if (Wrt_Unsteady && GridMovement) {
          for (iDim = 0; iDim < dims; iDim++) {
            err = TECDAT142(&NPts, Surf_Data[iDim], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
        for (iVar = dims; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      } else {
        if (Wrt_Unsteady && GridMovement) {
          for (iDim = 0; iDim < dims; iDim++) {
            err = TECDAT142(&NPts, Surf_Coords[iDim], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
        for (iVar = 0; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      }
      first_zone = false;
    }
    
  }
  
  if (nGlobal_BoundQuad > 0) {
    
    /*--- Write the zone header information ---*/
    ZoneType = FEQUADRILATERAL; NElm = (INTEGER4)nGlobal_BoundQuad;
    
    err = TECZNE142((char*)"Quadrilateral Elements",
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
      if ((config->GetKind_SU2() == SU2_SOL) || (config->GetKind_SU2() == SU2_DOT)) {
        if (Wrt_Unsteady && GridMovement) {
          for (iDim = 0; iDim < dims; iDim++) {
            err = TECDAT142(&NPts, Surf_Data[iDim], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
        for (iVar = dims; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      } else {
        if (Wrt_Unsteady && GridMovement) {
          for (iDim = 0; iDim < dims; iDim++) {
            err = TECDAT142(&NPts, Surf_Coords[iDim], &IsDouble); ShareFromZone[i++] = 1;
            if (err) cout << "Error writing data to Tecplot file" << endl;
          }
        }
        for (iVar = 0; iVar < nVar_Total; iVar++) {
          err = TECDAT142(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
          if (err) cout << "Error writing data to Tecplot file" << endl;
        }
      }
      first_zone = false;
    }
    
  }
  
  for (iVar = 0; iVar < nVar_Total; iVar++)
  delete [] Surf_Data[iVar];
  delete [] Surf_Data;
  
  if (Surf_Coords != NULL) {
    for (iDim = 0; iDim < dims; iDim++) delete [] Surf_Coords[iDim];
    delete [] Surf_Coords;
  }
  delete [] LocalIndex;
  delete [] SurfacePoint;
  delete [] ShareFromZone;
  
  err = TECEND142();
  if (err) cout << "Error in closing Tecplot file" << endl;
  
#endif
  
}

string COutput::AssembleVariableNames(CGeometry *geometry, CConfig *config, unsigned short nVar_Consv, unsigned short *NVar) {
  
  /*--- Local variables ---*/
  stringstream variables; variables.str(string());
  unsigned short iVar;
  *NVar = 0;
  unsigned short nDim = geometry->GetnDim();
  unsigned short Kind_Solver  = config->GetKind_Solver();
  bool grid_movement = config->GetGrid_Movement();
  bool Wrt_Unsteady = config->GetWrt_Unsteady();
  
  
  /*--- Write the basic variable header based on the particular solution ----*/
  
  /*--- Write the list of the fields in the restart file.
   Without including the PointID---*/
  if ((config->GetKind_SU2() == SU2_SOL) || (config->GetKind_SU2() == SU2_DOT)) {
    
    /*--- If SU2_SOL called this routine, we already have a set of output
     variables with the appropriate string tags stored in the config class.
     We simply read in and remove the quotation marks from the var names. ---*/
    
    /*--- Set the number of variables to be written. Subtract off an index for
     the PointID as well as each coordinate (x, y, z). ---*/
    string varname;
    
    if (Wrt_Unsteady && grid_movement) {
      
      *NVar = config->fields.size()-1;
      for (unsigned short iField = 1; iField < config->fields.size(); iField++) {
        varname = config->fields[iField];
        varname.erase (varname.begin(), varname.begin()+1);
        varname.erase (varname.end()-1, varname.end());
        variables << varname << " ";
      }
    } else {
      
      *NVar = config->fields.size()-1-nDim;
      for (unsigned short iField = 1+nDim; iField < config->fields.size(); iField++) {
        varname = config->fields[iField];
        varname.erase (varname.begin(), varname.begin()+1);
        varname.erase (varname.end()-1, varname.end());
        variables << varname << " ";
      }
    }
    
  } else {
    
    if (Wrt_Unsteady && grid_movement) {
      if (nDim == 2) {
        variables << "x y "; *NVar += 2;
      } else {
        variables << "x y z "; *NVar += 3;
      }
    }
    
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
      variables << "Conservative_" << iVar+1<<" "; *NVar += 1;
    }
    if (config->GetWrt_Limiters()) {
      for (iVar = 0; iVar < nVar_Consv; iVar++) {
        variables << "Limiter_" << iVar+1<<" "; *NVar += 1;
      }
    }
    if (config->GetWrt_Residuals()) {
      for (iVar = 0; iVar < nVar_Consv; iVar++) {
        variables << "Residual_" << iVar+1<<" "; *NVar += 1;
      }
    }
    
    /*--- Add names for any extra variables (this will need to be adjusted). ---*/
    if (grid_movement) {
      if (nDim == 2) {
        variables << "Grid_Velx Grid_Vely "; *NVar += 2;
      } else {
        variables << "Grid_Velx Grid_Vely Grid_Velz "; *NVar += 3;
      }
    }
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      variables << "Pressure Temperature Pressure_Coefficient Mach ";
      *NVar += 4;
    }
    
    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      if (nDim == 2) {
        variables << "Laminar_Viscosity Skin_Friction_Coefficient_x Skin_Friction_Coefficient_y Heat_Flux Y_Plus ";
        *NVar += 5;
      }
      else {
        variables << "Laminar_Viscosity Skin_Friction_Coefficient_x Skin_Friction_Coefficient_y Skin_Friction_Coefficient_z Heat_Flux Y_Plus ";
        *NVar += 6;
      }
    }
    
    if (Kind_Solver == RANS) {
      variables << "Eddy_Viscosity ";
      *NVar += 1;
    }
    
    if (config->GetWrt_SharpEdges()) {
      if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
        variables << "Sharp_Edge_Dist ";
        *NVar += 1;
      }
    }
    

    
    if (( Kind_Solver == ADJ_EULER              ) ||
        ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
        ( Kind_Solver == ADJ_RANS               )   ) {
      variables << "Surface_Sensitivity Solution_Sensor ";
      *NVar += 2;
    }
  }
  
  return variables.str();
  
}
