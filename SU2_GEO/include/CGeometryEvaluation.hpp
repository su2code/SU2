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
#include <map>

#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/CConfig.hpp"
#include "../../Common/include/grid_movement_structure.hpp"

using namespace std;

class CGeometryEvaluation {
private:
  CConfig **config_container          = nullptr;
  CGeometry **geometry_container      = nullptr;
  CGeometry* geometry                 = nullptr;
  CConfig* config                     = nullptr;

	unsigned short iZone, nZone = SINGLE_ZONE;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  unsigned short iDV, iFFDBox, iPlane, nPlane, iVar, iLoc, iFile;
  su2double *ObjectiveFunc, *ObjectiveFunc_New, *Gradient, delta_eps,
  **Plane_P0, **Plane_Normal;
  
  vector<su2double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;
  vector<su2double> Xcoord_Fan, Ycoord_Fan, Zcoord_Fan;
  char* config_file_name;                       /*!< \brief Configuration file name of the problem.*/  
  bool Local_MoveSurface, MoveSurface = false;
  ofstream Gradient_file, ObjFunc_file;
  int rank, size;
  bool fuselage = false, wing = false, nacelle = false, thickness = false, airfoil = false, tabTecplot = false;

  CSurfaceMovement *surface_movement  = nullptr;
  CFreeFormDefBox** FFDBox            = nullptr;

  vector<vector<su2double>> thicc, thicc_New, thicc_Grad;

  SU2_Comm MPICommunicator;

  map<string, su2double> Wing_ObjectiveFuncs = 
  {
    {"WING_VOLUME", 0.0},
    {"WING_MIN_THICKNESS", 0.0},
    {"WING_MAX_THICKNESS", 0.0},
    {"WING_MIN_CHORD", 0.0},
    {"WING_MAX_CHORD", 0.0},
    {"WING_MIN_LE_RADIUS", 0.0},
    {"WING_MAX_LE_RADIUS", 0.0},
    {"WING_MIN_TOC", 0.0},
    {"WING_MAX_TOC", 0.0},
    {"WING_OBJFUN_MIN_TOC", 0.0},
    {"WING_MAX_TWIST", 0.0},
    {"WING_MAX_CURVATURE", 0.0},
    {"WING_MAX_DIHEDRAL", 0.0}
  };
  map<string, su2double> Wing_ObjectiveFuncs_FD;
  map<string, su2double> Wing_ObjectiveFuncs_Grad;

  set<string> WingSection_Funcs = 
  {
    "AREA",
    "THICKNESS",
    "CHORD",
    "LE_RADIUS",
    "THICKNESS_OVER_CHORD",
    "AOA"
  };

  map<string, su2double> Fuselage_ObjectiveFuncs = 
  {
    {"FUSELAGE_VOLUME", 0.0},
    {"FUSELAGE_WETTED_AREA", 0.0},
    {"FUSELAGE_MIN_WIDTH", 0.0},
    {"FUSELAGE_MAX_WIDTH", 0.0},
    {"FUSELAGE_MIN_WATERLINE_WIDTH", 0.0},
    {"FUSELAGE_MAX_WATERLINE_WIDTH", 0.0},
    {"FUSELAGE_MIN_HEIGHT", 0.0},
    {"FUSELAGE_MAX_HEIGHT", 0.0},
    {"FUSELAGE_MAX_CURVATURE", 0.0},
  };
  map<string, su2double> Fuselage_ObjectiveFuncs_FD;
  map<string, su2double> Fuselage_ObjectiveFuncs_Grad;

  set<string> FuselageSection_Funcs = 
  {
    "AREA",
    "LENGTH",
    "WIDTH",
    "WATERLINE_WIDTH",
    "HEIGHT"
  };

  map<string, su2double> Nacelle_ObjectiveFuncs =
  {
    {"NACELLE_VOLUME", 0.0},
    {"NACELLE_MIN_THICKNESS", 0.0},
    {"NACELLE_MAX_THICKNESS", 0.0},
    {"NACELLE_MIN_CHORD", 0.0},
    {"NACELLE_MAX_CHORD", 0.0},
    {"NACELLE_MIN_LE_RADIUS", 0.0},
    {"NACELLE_MAX_LE_RADIUS", 0.0},
    {"NACELLE_MIN_TOC", 0.0},
    {"NACELLE_MAX_TOC", 0.0},
    {"NACELLE_OBJFUN_MIN_TOC", 0.0},
    {"NACELLE_MAX_TWIST", 0.0},
  };
  map<string, su2double> Nacelle_ObjectiveFuncs_FD;
  map<string, su2double> Nacelle_ObjectiveFuncs_Grad;

  set<string> NacelleSection_Funcs = 
  {
    "AREA",
    "THICKNESS",
    "CHORD",
    "LE_RADIUS",
    "THICKNESS_OVER_CHORD",
    "AOA"
  };


public:

	CGeometryEvaluation(char* confFile, unsigned short val_nZone, SU2_Comm val_MPICommunicator);

  ~CGeometryEvaluation(void);

  /*!
   * \brief Read in the config files.
   */
  void Config_Preprocessing();

  /*!
   * \brief Construction of the edge-based data structure
   */
  void Geometry_Preprocessing();

  void SetMapToZero(map<string,su2double> &Function_Map);
	
  void SetSectioningVariables(unsigned short nPlane_val);

	void ComputeGeometry(void);

  void EvaluateObjectiveFunction(void);

  void OutputFunctionFile(void);

  void EverythingGradient(void);

};
