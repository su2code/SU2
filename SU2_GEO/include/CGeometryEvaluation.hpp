/*!
 * \file CGeometryEvaluation.hpp
 * \brief Headers of the main subroutines of the CGeometryEvaluation class
 *        The subroutines and functions are in the <i>CGeometryEvaluation.cpp</i> file.
 * \author J. Mukhopadhaya
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#pragma once

#include "../../Common/include/mpi_structure.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/CConfig.hpp"
#include "../../Common/include/grid_movement_structure.hpp"

using namespace std;

class CGeometryEvaluation {
private:
	CGeometry* geometry;
	CConfig* config;
	unsigned short iZone, nZone = SINGLE_ZONE;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  unsigned short iDV, iFFDBox, iPlane, nPlane, iVar, iLoc, iFile;
  su2double *ObjectiveFunc, *ObjectiveFunc_New, *Gradient, delta_eps,
  **Plane_P0, **Plane_Normal,
  
  Fuselage_Volume = 0.0, Fuselage_WettedArea = 0.0, Fuselage_MinWidth = 0.0, Fuselage_MaxWidth = 0.0, Fuselage_MinWaterLineWidth = 0.0, Fuselage_MaxWaterLineWidth = 0.0, Fuselage_MinHeight = 0.0, Fuselage_MaxHeight = 0.0, Fuselage_MaxCurvature = 0.0,
  Fuselage_Volume_New = 0.0, Fuselage_WettedArea_New = 0.0, Fuselage_MinWidth_New = 0.0, Fuselage_MaxWidth_New = 0.0, Fuselage_MinWaterLineWidth_New = 0.0, Fuselage_MaxWaterLineWidth_New = 0.0, Fuselage_MinHeight_New = 0.0, Fuselage_MaxHeight_New = 0.0, Fuselage_MaxCurvature_New = 0.0,
  Fuselage_Volume_Grad = 0.0, Fuselage_WettedArea_Grad = 0.0, Fuselage_MinWidth_Grad = 0.0, Fuselage_MaxWidth_Grad = 0.0, Fuselage_MinWaterLineWidth_Grad = 0.0, Fuselage_MaxWaterLineWidth_Grad = 0.0, Fuselage_MinHeight_Grad = 0.0, Fuselage_MaxHeight_Grad = 0.0, Fuselage_MaxCurvature_Grad = 0.0,
  
  Wing_Volume = 0.0, Wing_MinThickness = 0.0, Wing_MaxThickness = 0.0, Wing_MinChord = 0.0, Wing_MaxChord = 0.0, Wing_MinLERadius = 0.0, Wing_MaxLERadius = 0.0, Wing_MinToC = 0.0, Wing_MaxToC = 0.0, Wing_ObjFun_MinToC = 0.0, Wing_MaxTwist = 0.0, Wing_MaxCurvature = 0.0, Wing_MaxDihedral = 0.0,
  Wing_Volume_New = 0.0, Wing_MinThickness_New = 0.0, Wing_MaxThickness_New = 0.0, Wing_MinChord_New = 0.0, Wing_MaxChord_New = 0.0, Wing_MinLERadius_New = 0.0, Wing_MaxLERadius_New = 0.0, Wing_MinToC_New = 0.0, Wing_MaxToC_New = 0.0, Wing_ObjFun_MinToC_New = 0.0, Wing_MaxTwist_New = 0.0, Wing_MaxCurvature_New = 0.0, Wing_MaxDihedral_New = 0.0,
  Wing_Volume_Grad = 0.0, Wing_MinThickness_Grad = 0.0, Wing_MaxThickness_Grad = 0.0, Wing_MinChord_Grad = 0.0, Wing_MaxChord_Grad = 0.0, Wing_MinLERadius_Grad = 0.0, Wing_MaxLERadius_Grad = 0.0, Wing_MinToC_Grad = 0.0, Wing_MaxToC_Grad = 0.0, Wing_ObjFun_MinToC_Grad = 0.0, Wing_MaxTwist_Grad = 0.0, Wing_MaxCurvature_Grad = 0.0, Wing_MaxDihedral_Grad = 0.0,
  
  Nacelle_Volume = 0.0, Nacelle_MinThickness = 0.0, Nacelle_MaxThickness = 0.0, Nacelle_MinChord = 0.0, Nacelle_MaxChord = 0.0, Nacelle_MinLERadius = 0.0, Nacelle_MaxLERadius = 0.0, Nacelle_MinToC = 0.0, Nacelle_MaxToC = 0.0, Nacelle_ObjFun_MinToC = 0.0, Nacelle_MaxTwist = 0.0,
  Nacelle_Volume_New = 0.0, Nacelle_MinThickness_New = 0.0, Nacelle_MaxThickness_New = 0.0, Nacelle_MinChord_New = 0.0, Nacelle_MaxChord_New = 0.0, Nacelle_MinLERadius_New = 0.0, Nacelle_MaxLERadius_New = 0.0, Nacelle_MinToC_New = 0.0, Nacelle_MaxToC_New = 0.0, Nacelle_ObjFun_MinToC_New = 0.0, Nacelle_MaxTwist_New = 0.0,
  Nacelle_Volume_Grad = 0.0, Nacelle_MinThickness_Grad = 0.0, Nacelle_MaxThickness_Grad = 0.0, Nacelle_MinChord_Grad = 0.0, Nacelle_MaxChord_Grad = 0.0, Nacelle_MinLERadius_Grad = 0.0, Nacelle_MaxLERadius_Grad = 0.0, Nacelle_MinToC_Grad = 0.0, Nacelle_MaxToC_Grad = 0.0, Nacelle_ObjFun_MinToC_Grad = 0.0, Nacelle_MaxTwist_Grad = 0.0;

  vector<su2double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;
  vector<su2double> Xcoord_Fan, Ycoord_Fan, Zcoord_Fan;
  char config_file_name[MAX_STRING_SIZE];
  bool Local_MoveSurface, MoveSurface = false;
  ofstream Gradient_file, ObjFunc_file;
  int rank, size;
  bool fuselage = false, wing = false, nacelle = false, thickness = false, airfoil = false, tabTecplot = false;

  CSurfaceMovement *surface_movement  = NULL;
  CFreeFormDefBox** FFDBox            = NULL;

  vector<vector<su2double>> thicc, thicc_New, thicc_Grad;

  SU2_Comm MPICommunicator;


public:

	CGeometryEvaluation(CGeometry* geo_val, CConfig* config_val, SU2_Comm val_MPICommunicator);

  ~CGeometryEvaluation(void);

	// void SetPlaneStructure(void);

	void ComputeGeometry(void);

  void EvaluateObjectiveFunction(void);

  void OutputFunctionFile(void);

  void EverythingGradient(void);

};
