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

void COutput::WriteTecplotASCII_Parallel(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_nZone, unsigned short val_iInst, unsigned short val_nInst, bool surf_sol) {
  
  unsigned short iVar, nDim = geometry->GetnDim();
  unsigned short Kind_Solver = config->GetKind_Solver();
  
  unsigned long iPoint, iElem, iNode;
  unsigned long iExtIter = config->GetExtIter();
  
  bool adjoint = config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint();
  
  int iProcessor;

  char cstr[200], buffer[50];
  string filename;
  ofstream Tecplot_File;

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
  
  if (config->GetKind_SU2() == SU2_DOT) {
    if (surf_sol) filename = config->GetSurfSens_FileName();
    else filename = config->GetVolSens_FileName();
  }
  
  strcpy (cstr, filename.c_str());
  
  /*--- Special cases where a number needs to be appended to the file name. ---*/
  
  if ((Kind_Solver == EULER || Kind_Solver == NAVIER_STOKES || Kind_Solver == RANS ||
       Kind_Solver == ADJ_EULER || Kind_Solver == ADJ_NAVIER_STOKES || Kind_Solver == ADJ_RANS ||
       Kind_Solver == DISC_ADJ_EULER || Kind_Solver == DISC_ADJ_NAVIER_STOKES || Kind_Solver == DISC_ADJ_RANS) &&
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
  
  if (rank == MASTER_NODE) {
    Tecplot_File.open(cstr, ios::out);
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
  
  Tecplot_File.open(cstr, ios::out | ios::app);
  
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

    file << buffer << ".mesh.plt";
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
  file << ".sol.plt";
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
