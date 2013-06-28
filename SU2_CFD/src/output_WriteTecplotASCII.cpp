/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.2
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

void COutput::WriteTecplotASCII(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone) {
    
    /*--- Local variables and initialization ---*/
    unsigned short iDim, iVar, nDim = geometry->GetnDim();
    unsigned short Kind_Solver = config->GetKind_Solver();
    
    unsigned long iPoint, iElem, iNode;
    unsigned long iExtIter = config->GetExtIter();
    
    bool grid_movement  = config->GetGrid_Movement();
    bool isAdjoint = config->IsAdjoint();
    
    char cstr[200], buffer[50];
    string filename;
    
    /*--- Write file name with extension ---*/
    if (isAdjoint)
        filename = config->GetAdj_FileName();
    else
        filename = config->GetFlow_FileName();
    
    if (Kind_Solver == LINEAR_ELASTICITY)
        filename = config->GetStructure_FileName().c_str();
    if (Kind_Solver == WAVE_EQUATION)
        filename = config->GetWave_FileName().c_str();
    if ((Kind_Solver == WAVE_EQUATION) && (Kind_Solver == ADJ_AEROACOUSTIC_EULER))
        filename = config->GetAdjWave_FileName().c_str();
    if (Kind_Solver == ELECTRIC_POTENTIAL)
        filename = config->GetStructure_FileName().c_str();
    
#ifndef NO_MPI
    /*--- Remove the domain number from the surface csv filename ---*/
    int nProcessor = MPI::COMM_WORLD.Get_size();
    if (nProcessor > 1) filename.erase (filename.end()-2, filename.end());
#endif
    
    strcpy (cstr, filename.c_str());
    
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
    Tecplot_File << "TITLE = \"Visualization of the volumetric solution\"" << endl;
    
    /*--- Prepare the variables lists. ---*/
    if (nDim == 2) {
        Tecplot_File << "VARIABLES = \"x\",\"y\"";
    } else {
        Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\"";
    }
    
    /*--- Add temporary names for conservative and residual variables ---*/
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
        Tecplot_File << ",\"Conservative_" << iVar+1 << "\"";
    }
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
        Tecplot_File << ",\"Residual_" << iVar+1 << "\"";
    }
    
    /*--- Add names for any extra variables (this will need to be adjusted). ---*/
	if (grid_movement) {
        if (nDim == 2) {
            Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\"";
        } else {
            Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\",\"Grid_Velz\"";
        }
	}
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
        (Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
        Tecplot_File << ",\"Pressure\",\"Mach\"";
    }
    
    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
        (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
        Tecplot_File << ",\"Temperature\",\"Laminar Viscosity\"";
    }
    
    if ((Kind_Solver == RANS) || (Kind_Solver == FREE_SURFACE_RANS)) {
        Tecplot_File << ",\"Eddy Viscosity\"";
    }
    
    Tecplot_File << endl;
    
    Tecplot_File << "ZONE ";
    if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady())
        Tecplot_File << "STRANDID="<<int(iExtIter+1)<<", SOLUTIONTIME="<<config->GetDelta_UnstTime()*iExtIter<<", ";
    
    if (nDim == 2) {
        Tecplot_File << "NODES= "<< nGlobal_Poin <<", ELEMENTS= "<< nGlobal_Elem <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
    } else {
        Tecplot_File << "NODES= "<< nGlobal_Poin <<", ELEMENTS= "<< nGlobal_Elem <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
    }
    
    /*--- Write solution data. ---*/
    for (iPoint = 0; iPoint < nGlobal_Poin; iPoint++) {
        
        /*--- Write the node coordinates ---*/
        for(iDim = 0; iDim < nDim; iDim++)
            Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";
        
        /*--- Loop over the vars/residuals and write the values to file ---*/
        for (iVar = 0; iVar < 2*nVar_Consv; iVar++) {
            Tecplot_File << scientific << Data[iVar][iPoint] << "\t";
        }
        
        /*--- Loop over any remaining variables ---*/
        for (iVar = 2*nVar_Consv; iVar < nVar_Total; iVar++) {
            Tecplot_File << scientific << Data[iVar][iPoint] << "\t";
        }
        Tecplot_File << endl;
    }
    
    /*--- Write connectivity data. ---*/
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Tria; iElem++) {
        iNode = iElem*N_POINTS_TRIANGLE;
        Tecplot_File << Conn_Tria[iNode+0]+1 << "\t";
        Tecplot_File << Conn_Tria[iNode+1]+1 << "\t";
        Tecplot_File << Conn_Tria[iNode+2]+1 << "\t";
        Tecplot_File << Conn_Tria[iNode+2]+1 << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Quad; iElem++) {
        iNode = iElem*N_POINTS_QUADRILATERAL;
        Tecplot_File << Conn_Quad[iNode+0]+1 << "\t";
        Tecplot_File << Conn_Quad[iNode+1]+1 << "\t";
        Tecplot_File << Conn_Quad[iNode+2]+1 << "\t";
        Tecplot_File << Conn_Quad[iNode+3]+1 << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Tetr; iElem++) {
        iNode = iElem*N_POINTS_TETRAHEDRON;
        Tecplot_File << Conn_Tetr[iNode+0]+1 << "\t" << Conn_Tetr[iNode+1]+1 << "\t";
        Tecplot_File << Conn_Tetr[iNode+2]+1 << "\t" << Conn_Tetr[iNode+2]+1 << "\t";
        Tecplot_File << Conn_Tetr[iNode+3]+1 << "\t" << Conn_Tetr[iNode+3]+1 << "\t";
        Tecplot_File << Conn_Tetr[iNode+3]+1 << "\t" << Conn_Tetr[iNode+3]+1 << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Hexa; iElem++) {
        iNode = iElem*N_POINTS_HEXAHEDRON;
        Tecplot_File << Conn_Hexa[iNode+0]+1 << "\t" << Conn_Hexa[iNode+1]+1 << "\t";
        Tecplot_File << Conn_Hexa[iNode+2]+1 << "\t" << Conn_Hexa[iNode+3]+1 << "\t";
        Tecplot_File << Conn_Hexa[iNode+4]+1 << "\t" << Conn_Hexa[iNode+5]+1 << "\t";
        Tecplot_File << Conn_Hexa[iNode+6]+1 << "\t" << Conn_Hexa[iNode+7]+1 << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Wedg; iElem++) {
        iNode = iElem*N_POINTS_WEDGE;
        Tecplot_File << Conn_Wedg[iNode+0]+1 << "\t" << Conn_Wedg[iNode+1]+1 << "\t";
        Tecplot_File << Conn_Wedg[iNode+1]+1 << "\t" << Conn_Wedg[iNode+2]+1 << "\t";
        Tecplot_File << Conn_Wedg[iNode+3]+1 << "\t" << Conn_Wedg[iNode+4]+1 << "\t";
        Tecplot_File << Conn_Wedg[iNode+4]+1 << "\t" << Conn_Wedg[iNode+5]+1 << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Pyra; iElem++) {
        iNode = iElem*N_POINTS_PYRAMID;
        Tecplot_File << Conn_Pyra[iNode+0]+1 << "\t" << Conn_Pyra[iNode+1]+1 << "\t";
        Tecplot_File << Conn_Pyra[iNode+2]+1 << "\t" << Conn_Pyra[iNode+3]+1 << "\t";
        Tecplot_File << Conn_Pyra[iNode+4]+1 << "\t" << Conn_Pyra[iNode+4]+1 << "\t";
        Tecplot_File << Conn_Pyra[iNode+4]+1 << "\t" << Conn_Pyra[iNode+4]+1 << "\n";
    }
    
    Tecplot_File.close();
    
}