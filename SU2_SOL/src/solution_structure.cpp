/*!
 * \file solution_direct_mean.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

#include "../include/solution_structure.hpp"

CSolution::CSolution(void) {}

CSolution::~CSolution(void) {}

unsigned short CSolution::GetnVar(void) { return nVar; }

CBaselineSolution::CBaselineSolution(void) : CSolution() { }

CBaselineSolution::CBaselineSolution(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolution() {
	unsigned long iPoint, index, flowIter, iPoint_Global;
  long iPoint_Local;
  double Solution[10];
  char buffer[50];
  unsigned short iField;
  string Tag, text_line;
  
	/*--- Define geometry constants in the solver structure ---*/
	nDim = geometry->GetnDim();
  
  /*--- The number of variables should be readed from the config file ---*/
	nVar = nDim + 2;
  
	/*--- Allocate the node variables ---*/
	node = new CVariable*[geometry->GetnPoint()];
  
  /*--- Restart the solution from file information ---*/
  ifstream restart_file;
  string filename = config->GetSolution_FlowFileName();
  
  /*--- Append time step for unsteady restart ---*/
  if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    flowIter = config->GetnExtIter() - 1;
    filename.erase (filename.end()-4, filename.end());
    if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
    if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
    if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
    if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
    if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
    string UnstExt = string(buffer);
    filename.append(UnstExt);
  }
  restart_file.open(filename.data(), ios::in);
  
  /*--- In case there is no restart file ---*/
  if (restart_file.fail()) {
    cout << "There is no flow restart file!!" << endl;
    cout << "Press any key to exit..." << endl;
    cin.get(); exit(1);
  }
  
  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
  long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];
  
  /*--- First, set all indices to a negative value by default ---*/
  for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
    Global2Local[iPoint] = -1;
  
  /*--- Now fill array with the transform values only for local points ---*/
  for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
  
  
  /*--- Identify the number of fields (and names) in the restart file ---*/
  getline (restart_file, text_line);
  stringstream ss(text_line);
  while (ss >> Tag) {
    config->fields.push_back(Tag);
    if (ss.peek() == ',') ss.ignore();
  }

  /*--- Set the number of variables, one per field in the 
   restart file (without including the PointID) ---*/
  nVar = config->fields.size() - 1;
  
  /*--- Read all lines in the restart file ---*/
  iPoint_Global = 0;
  while (getline (restart_file, text_line)) {
    istringstream point_line(text_line);

    /*--- Retrieve local index. If this node from the restart file lives
     on a different processor, the value of iPoint_Local will be -1.
     Otherwise, the local index for this node on the current processor
     will be returned and used to instantiate the vars. ---*/
    iPoint_Local = Global2Local[iPoint_Global];
    if (iPoint_Local >= 0) {
      
      /*--- The PointID is not stored --*/
      point_line >> index;
      
      /*--- Store the solution --*/
      for (iField = 0; iField < nVar; iField++)
        point_line >> Solution[iField];
      
      node[iPoint_Local] = new CBaselineVariable(Solution, nDim, nVar, config);
    }
    iPoint_Global++;
  }
  
  /*--- Instantiate the variable class with an arbitrary solution
   at any halo/periodic nodes. The initial solution can be arbitrary,
   because a send/recv is performed immediately in the solver. ---*/
  for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++)
    node[iPoint] = new CBaselineVariable(Solution, nDim, nVar, config);
  
  /*--- Close the restart file ---*/
  restart_file.close();
  
  /*--- Free memory needed for the transformation ---*/
  delete [] Global2Local;
  
}

CBaselineSolution::~CBaselineSolution(void) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < nPoint; iPoint++)
		delete node[iPoint];
	delete [] node;

}