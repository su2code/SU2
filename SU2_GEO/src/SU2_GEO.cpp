/*!
 * \file SU2_GEO.cpp
 * \brief Main file of the Geometry Definition Code (SU2_GEO).
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../include/SU2_GEO.hpp"
using namespace std;

int main(int argc, char* argv[]) {
  unsigned short iZone, nZone = SINGLE_ZONE;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  unsigned short iDV, iFFDBox, iPlane, nPlane, iVar;
  su2double *ObjectiveFunc, *ObjectiveFunc_New, *Gradient, delta_eps, **Plane_P0, **Plane_Normal,

      Fuselage_Volume = 0.0, Fuselage_WettedArea = 0.0, Fuselage_MinWidth = 0.0, Fuselage_MaxWidth = 0.0,
      Fuselage_MinWaterLineWidth = 0.0, Fuselage_MaxWaterLineWidth = 0.0, Fuselage_MinHeight = 0.0,
      Fuselage_MaxHeight = 0.0, Fuselage_MaxCurvature = 0.0, Fuselage_Volume_New = 0.0, Fuselage_WettedArea_New = 0.0,
      Fuselage_MinWidth_New = 0.0, Fuselage_MaxWidth_New = 0.0, Fuselage_MinWaterLineWidth_New = 0.0,
      Fuselage_MaxWaterLineWidth_New = 0.0, Fuselage_MinHeight_New = 0.0, Fuselage_MaxHeight_New = 0.0,
      Fuselage_MaxCurvature_New = 0.0, Fuselage_Volume_Grad = 0.0, Fuselage_WettedArea_Grad = 0.0,
      Fuselage_MinWidth_Grad = 0.0, Fuselage_MaxWidth_Grad = 0.0, Fuselage_MinWaterLineWidth_Grad = 0.0,
      Fuselage_MaxWaterLineWidth_Grad = 0.0, Fuselage_MinHeight_Grad = 0.0, Fuselage_MaxHeight_Grad = 0.0,
      Fuselage_MaxCurvature_Grad = 0.0,

      Wing_Volume = 0.0, Wing_MinThickness = 0.0, Wing_MaxThickness = 0.0, Wing_MinChord = 0.0, Wing_MaxChord = 0.0,
      Wing_MinLERadius = 0.0, Wing_MaxLERadius = 0.0, Wing_MinToC = 0.0, Wing_MaxToC = 0.0, Wing_ObjFun_MinToC = 0.0,
      Wing_MaxTwist = 0.0, Wing_MaxCurvature = 0.0, Wing_MaxDihedral = 0.0, Wing_Volume_New = 0.0,
      Wing_MinThickness_New = 0.0, Wing_MaxThickness_New = 0.0, Wing_MinChord_New = 0.0, Wing_MaxChord_New = 0.0,
      Wing_MinLERadius_New = 0.0, Wing_MaxLERadius_New = 0.0, Wing_MinToC_New = 0.0, Wing_MaxToC_New = 0.0,
      Wing_ObjFun_MinToC_New = 0.0, Wing_MaxTwist_New = 0.0, Wing_MaxCurvature_New = 0.0, Wing_MaxDihedral_New = 0.0,
      Wing_Volume_Grad = 0.0, Wing_MinThickness_Grad = 0.0, Wing_MaxThickness_Grad = 0.0, Wing_MinChord_Grad = 0.0,
      Wing_MaxChord_Grad = 0.0, Wing_MinLERadius_Grad = 0.0, Wing_MaxLERadius_Grad = 0.0, Wing_MinToC_Grad = 0.0,
      Wing_MaxToC_Grad = 0.0, Wing_ObjFun_MinToC_Grad = 0.0, Wing_MaxTwist_Grad = 0.0, Wing_MaxCurvature_Grad = 0.0,
      Wing_MaxDihedral_Grad = 0.0,

      Nacelle_Volume = 0.0, Nacelle_MinThickness = 0.0, Nacelle_MaxThickness = 0.0, Nacelle_MinChord = 0.0,
      Nacelle_MaxChord = 0.0, Nacelle_MinLERadius = 0.0, Nacelle_MaxLERadius = 0.0, Nacelle_MinToC = 0.0,
      Nacelle_MaxToC = 0.0, Nacelle_ObjFun_MinToC = 0.0, Nacelle_MaxTwist = 0.0, Nacelle_Volume_New = 0.0,
      Nacelle_MinThickness_New = 0.0, Nacelle_MaxThickness_New = 0.0, Nacelle_MinChord_New = 0.0,
      Nacelle_MaxChord_New = 0.0, Nacelle_MinLERadius_New = 0.0, Nacelle_MaxLERadius_New = 0.0,
      Nacelle_MinToC_New = 0.0, Nacelle_MaxToC_New = 0.0, Nacelle_ObjFun_MinToC_New = 0.0, Nacelle_MaxTwist_New = 0.0,
      Nacelle_Volume_Grad = 0.0, Nacelle_MinThickness_Grad = 0.0, Nacelle_MaxThickness_Grad = 0.0,
      Nacelle_MinChord_Grad = 0.0, Nacelle_MaxChord_Grad = 0.0, Nacelle_MinLERadius_Grad = 0.0,
      Nacelle_MaxLERadius_Grad = 0.0, Nacelle_MinToC_Grad = 0.0, Nacelle_MaxToC_Grad = 0.0,
      Nacelle_ObjFun_MinToC_Grad = 0.0, Nacelle_MaxTwist_Grad = 0.0;

  vector<su2double>*Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;
  vector<su2double> Xcoord_Fan, Ycoord_Fan, Zcoord_Fan;
  char config_file_name[MAX_STRING_SIZE];
  bool Local_MoveSurface, MoveSurface = false;
  ofstream Gradient_file, ObjFunc_file;
  int rank, size;

  /*--- MPI initialization ---*/

  SU2_MPI::Init(&argc, &argv);
  SU2_MPI::Comm MPICommunicator = SU2_MPI::GetComm();

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  /*--- Pointer to different structures that will be used throughout the entire code ---*/

  CConfig** config_container = nullptr;
  CGeometry** geometry_container = nullptr;
  CSurfaceMovement* surface_movement = nullptr;
  CFreeFormDefBox** FFDBox = nullptr;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/

  if (argc == 2) {
    strcpy(config_file_name, argv[1]);
  } else {
    strcpy(config_file_name, "default.cfg");
  }

  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation)  ---*/

  CConfig* config = nullptr;
  config = new CConfig(config_file_name, SU2_COMPONENT::SU2_GEO);

  nZone = config->GetnZone();

  /*--- Definition of the containers per zones ---*/

  config_container = new CConfig*[nZone];
  geometry_container = new CGeometry*[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone] = nullptr;
    geometry_container[iZone] = nullptr;
  }

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/

    config_container[iZone] = new CConfig(config_file_name, SU2_COMPONENT::SU2_GEO, true);
    config_container[iZone]->SetMPICommunicator(MPICommunicator);

    /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/

    CGeometry* geometry_aux = nullptr;

    /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

    geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);

    /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

    geometry_aux->SetColorGrid_Parallel(config_container[iZone]);

    /*--- Build the grid data structures using the ParMETIS coloring. ---*/

    geometry_container[iZone] = new CPhysicalGeometry(geometry_aux, config_container[iZone]);

    /*--- Deallocate the memory of geometry_aux ---*/

    delete geometry_aux;

    /*--- Add the Send/Receive boundaries ---*/

    geometry_container[iZone]->SetSendReceive(config_container[iZone]);

    /*--- Add the Send/Receive boundaries ---*/

    geometry_container[iZone]->SetBoundaries(config_container[iZone]);
  }

  bool tabTecplot = config_container[ZONE_0]->GetTabular_FileFormat() == TAB_OUTPUT::TAB_TECPLOT;

  /*--- Set up a timer for performance benchmarking (preprocessing time is included) ---*/

  StartTime = SU2_MPI::Wtime();

  /*--- Evaluation of the objective function ---*/

  if (rank == MASTER_NODE)
    cout << endl << "----------------------- Preprocessing computations ----------------------" << endl;

  /*--- Set the number of sections, and allocate the memory ---*/

  if (geometry_container[ZONE_0]->GetnDim() == 2)
    nPlane = 1;
  else
    nPlane = config_container[ZONE_0]->GetnLocationStations();

  Xcoord_Airfoil = new vector<su2double>[nPlane];
  Ycoord_Airfoil = new vector<su2double>[nPlane];
  Zcoord_Airfoil = new vector<su2double>[nPlane];
  Variable_Airfoil = new vector<su2double>[nPlane];

  Plane_P0 = new su2double*[nPlane];
  Plane_Normal = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    Plane_P0[iPlane] = new su2double[3];
    Plane_Normal[iPlane] = new su2double[3];
  }

  ObjectiveFunc = new su2double[nPlane * 20];
  ObjectiveFunc_New = new su2double[nPlane * 20];
  Gradient = new su2double[nPlane * 20];

  for (iVar = 0; iVar < nPlane * 20; iVar++) {
    ObjectiveFunc[iVar] = 0.0;
    ObjectiveFunc_New[iVar] = 0.0;
    Gradient[iVar] = 0.0;
  }

  /*--- Compute elements surrounding points, points surrounding points ---*/

  if (rank == MASTER_NODE) cout << "Setting local point connectivity." << endl;
  geometry_container[ZONE_0]->SetPoint_Connectivity();

  /*--- Check the orientation before computing geometrical quantities ---*/

  if (config_container[ZONE_0]->GetReorientElements()) {
    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation of the interior elements." << endl;
    geometry_container[ZONE_0]->Check_IntElem_Orientation(config_container[ZONE_0]);
  }

  /*--- Create the edge structure ---*/

  if (rank == MASTER_NODE) cout << "Identify edges and vertices." << endl;
  geometry_container[ZONE_0]->SetEdges();
  geometry_container[ZONE_0]->SetVertex(config_container[ZONE_0]);

  /*--- Create the dual control volume structures ---*/

  if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
  geometry_container[ZONE_0]->SetBoundControlVolume(config_container[ZONE_0], ALLOCATE);

  /*--- Compute the surface curvature ---*/

  if (rank == MASTER_NODE) cout << "Compute the surface curvature." << endl;
  geometry_container[ZONE_0]->ComputeSurf_Curvature(config_container[ZONE_0]);

  /*--- Computation of positive surface area in the z-plane ---*/

  if (rank == MASTER_NODE) cout << "Setting reference area and span." << endl;
  geometry_container[ZONE_0]->SetPositive_ZArea(config_container[ZONE_0]);

  /*--- Create the point-to-point MPI communication structures. ---*/
  geometry_container[ZONE_0]->PreprocessP2PComms(geometry_container[ZONE_0], config_container[ZONE_0]);

  /*--- Create plane structure ---*/

  if (rank == MASTER_NODE) cout << "Set plane structure." << endl;

  if (geometry_container[ZONE_0]->GetnDim() == 2) {
    Plane_Normal[0][0] = 0.0;
    Plane_P0[0][0] = 0.0;
    Plane_Normal[0][1] = 1.0;
    Plane_P0[0][1] = 0.0;
    Plane_Normal[0][2] = 0.0;
    Plane_P0[0][2] = 0.0;
  } else if (geometry_container[ZONE_0]->GetnDim() == 3) {
    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      Plane_Normal[iPlane][0] = 0.0;
      Plane_P0[iPlane][0] = 0.0;
      Plane_Normal[iPlane][1] = 0.0;
      Plane_P0[iPlane][1] = 0.0;
      Plane_Normal[iPlane][2] = 0.0;
      Plane_P0[iPlane][2] = 0.0;

      if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
        Plane_Normal[iPlane][0] = 1.0;
        Plane_P0[iPlane][0] = config_container[ZONE_0]->GetLocationStations(iPlane);
      } else if (config_container[ZONE_0]->GetGeo_Description() == NACELLE) {
        Plane_Normal[iPlane][0] = 0.0;
        Plane_Normal[iPlane][1] = -sin(config_container[ZONE_0]->GetLocationStations(iPlane) * PI_NUMBER / 180.0);
        Plane_Normal[iPlane][2] = cos(config_container[ZONE_0]->GetLocationStations(iPlane) * PI_NUMBER / 180.0);

        /*--- Apply tilt angle to the plane ---*/

        su2double Tilt_Angle = config_container[ZONE_0]->GetNacelleLocation(3) * PI_NUMBER / 180;
        su2double Plane_NormalX_Tilt =
            Plane_Normal[iPlane][0] * cos(Tilt_Angle) + Plane_Normal[iPlane][2] * sin(Tilt_Angle);
        su2double Plane_NormalY_Tilt = Plane_Normal[iPlane][1];
        su2double Plane_NormalZ_Tilt =
            Plane_Normal[iPlane][2] * cos(Tilt_Angle) - Plane_Normal[iPlane][0] * sin(Tilt_Angle);

        /*--- Apply toe angle to the plane ---*/

        su2double Toe_Angle = config_container[ZONE_0]->GetNacelleLocation(4) * PI_NUMBER / 180;
        su2double Plane_NormalX_Tilt_Toe = Plane_NormalX_Tilt * cos(Toe_Angle) - Plane_NormalY_Tilt * sin(Toe_Angle);
        su2double Plane_NormalY_Tilt_Toe = Plane_NormalX_Tilt * sin(Toe_Angle) + Plane_NormalY_Tilt * cos(Toe_Angle);
        su2double Plane_NormalZ_Tilt_Toe = Plane_NormalZ_Tilt;

        /*--- Update normal vector ---*/

        Plane_Normal[iPlane][0] = Plane_NormalX_Tilt_Toe;
        Plane_Normal[iPlane][1] = Plane_NormalY_Tilt_Toe;
        Plane_Normal[iPlane][2] = Plane_NormalZ_Tilt_Toe;

        /*--- Point in the plane ---*/

        Plane_P0[iPlane][0] = config_container[ZONE_0]->GetNacelleLocation(0);
        Plane_P0[iPlane][1] = config_container[ZONE_0]->GetNacelleLocation(1);
        Plane_P0[iPlane][2] = config_container[ZONE_0]->GetNacelleLocation(2);
      } else {
        Plane_Normal[iPlane][1] = 1.0;
        Plane_P0[iPlane][1] = config_container[ZONE_0]->GetLocationStations(iPlane);
      }
    }
  }

  /*--- Compute the wing and fan description (only 3D). ---*/

  if (geometry_container[ZONE_0]->GetnDim() == 3) {
    if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
      if (rank == MASTER_NODE) {
        cout << "Computing the fuselage continuous description." << endl << endl;
      }

      geometry_container[ZONE_0]->Compute_Fuselage(config_container[ZONE_0], true, Fuselage_Volume, Fuselage_WettedArea,
                                                   Fuselage_MinWidth, Fuselage_MaxWidth, Fuselage_MinWaterLineWidth,
                                                   Fuselage_MaxWaterLineWidth, Fuselage_MinHeight, Fuselage_MaxHeight,
                                                   Fuselage_MaxCurvature);

      /*--- Screen output for the wing definition ---*/

      if (rank == MASTER_NODE) {
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Fuselage volume: " << Fuselage_Volume << " in^3. ";
        else
          cout << "Fuselage volume: " << Fuselage_Volume << " m^3. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Fuselage wetted area: " << Fuselage_WettedArea << " in^2. " << endl;
        else
          cout << "Fuselage wetted area: " << Fuselage_WettedArea << " m^2. " << endl;
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Fuselage min. width: " << Fuselage_MinWidth << " in. ";
        else
          cout << "Fuselage min. width: " << Fuselage_MinWidth << " m. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Fuselage max. width: " << Fuselage_MaxWidth << " in. " << endl;
        else
          cout << "Fuselage max. width: " << Fuselage_MaxWidth << " m. " << endl;
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Fuselage min. waterline width: " << Fuselage_MinWaterLineWidth << " in. ";
        else
          cout << "Fuselage min. waterline width: " << Fuselage_MinWaterLineWidth << " m. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Fuselage max. waterline width: " << Fuselage_MaxWaterLineWidth << " in. " << endl;
        else
          cout << "Fuselage max. waterline width: " << Fuselage_MaxWaterLineWidth << " m. " << endl;
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Fuselage min. height: " << Fuselage_MinHeight << " in. ";
        else
          cout << "Fuselage min. height: " << Fuselage_MinHeight << " m. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Fuselage max. height: " << Fuselage_MaxHeight << " in. " << endl;
        else
          cout << "Fuselage max. height: " << Fuselage_MaxHeight << " m. " << endl;
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Fuselage max. curvature: " << Fuselage_MaxCurvature << " 1/in. " << endl;
        else
          cout << "Fuselage max. curvature: " << Fuselage_MaxCurvature << " 1/m. " << endl;
      }

    }

    else if (config_container[ZONE_0]->GetGeo_Description() == NACELLE) {
      if (rank == MASTER_NODE) {
        cout << "Computing the nacelle continuous description." << endl << endl;
      }

      geometry_container[ZONE_0]->Compute_Nacelle(config_container[ZONE_0], true, Nacelle_Volume, Nacelle_MinThickness,
                                                  Nacelle_MaxThickness, Nacelle_MinChord, Nacelle_MaxChord,
                                                  Nacelle_MinLERadius, Nacelle_MaxLERadius, Nacelle_MinToC,
                                                  Nacelle_MaxToC, Nacelle_ObjFun_MinToC, Nacelle_MaxTwist);

      /*--- Screen output for the wing definition ---*/

      if (rank == MASTER_NODE) {
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Nacelle volume: " << Nacelle_Volume << " in^3. ";
        else
          cout << "Nacelle volume: " << Nacelle_Volume << " m^3. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Nacelle min. thickness: " << Nacelle_MinThickness << " in. ";
        else
          cout << "Nacelle min. thickness: " << Nacelle_MinThickness << " m. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Nacelle max. thickness: " << Nacelle_MaxThickness << " in. " << endl;
        else
          cout << "Nacelle max. thickness: " << Nacelle_MaxThickness << " m. " << endl;
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Nacelle min. chord: " << Nacelle_MinChord << " in. ";
        else
          cout << "Nacelle min. chord: " << Nacelle_MinChord << " m. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Nacelle max. chord: " << Nacelle_MaxChord << " in. ";
        else
          cout << "Nacelle max. chord: " << Nacelle_MaxChord << " m. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Nacelle min. LE radius: " << Nacelle_MinLERadius << " 1/in. ";
        else
          cout << "Nacelle min. LE radius: " << Nacelle_MinLERadius << " 1/m. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Nacelle max. LE radius: " << Nacelle_MaxLERadius << " 1/in. " << endl;
        else
          cout << "Nacelle max. LE radius: " << Nacelle_MaxLERadius << " 1/m. " << endl;
        cout << "Nacelle min. ToC: " << Nacelle_MinToC << ". ";
        cout << "Nacelle max. ToC: " << Nacelle_MaxToC << ". ";
        cout << "Nacelle delta ToC: " << Nacelle_ObjFun_MinToC << ". ";
        cout << "Nacelle max. twist: " << Nacelle_MaxTwist << " deg. " << endl;
      }

    }

    else {
      if (rank == MASTER_NODE) {
        cout << "Computing the wing continuous description." << endl << endl;
      }

      geometry_container[ZONE_0]->Compute_Wing(config_container[ZONE_0], true, Wing_Volume, Wing_MinThickness,
                                               Wing_MaxThickness, Wing_MinChord, Wing_MaxChord, Wing_MinLERadius,
                                               Wing_MaxLERadius, Wing_MinToC, Wing_MaxToC, Wing_ObjFun_MinToC,
                                               Wing_MaxTwist, Wing_MaxCurvature, Wing_MaxDihedral);

      /*--- Screen output for the wing definition ---*/

      if (rank == MASTER_NODE) {
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Wing volume: " << Wing_Volume << " in^3. ";
        else
          cout << "Wing volume: " << Wing_Volume << " m^3. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Wing min. thickness: " << Wing_MinThickness << " in. ";
        else
          cout << "Wing min. thickness: " << Wing_MinThickness << " m. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Wing max. thickness: " << Wing_MaxThickness << " in. " << endl;
        else
          cout << "Wing max. thickness: " << Wing_MaxThickness << " m. " << endl;
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Wing min. chord: " << Wing_MinChord << " in. ";
        else
          cout << "Wing min. chord: " << Wing_MinChord << " m. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Wing max. chord: " << Wing_MaxChord << " in. ";
        else
          cout << "Wing max. chord: " << Wing_MaxChord << " m. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Wing min. LE radius: " << Wing_MinLERadius << " 1/in. ";
        else
          cout << "Wing min. LE radius: " << Wing_MinLERadius << " 1/m. ";
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Wing max. LE radius: " << Wing_MaxLERadius << " 1/in. " << endl;
        else
          cout << "Wing max. LE radius: " << Wing_MaxLERadius << " 1/m. " << endl;
        cout << "Wing min. ToC: " << Wing_MinToC << ". ";
        cout << "Wing max. ToC: " << Wing_MaxToC << ". ";
        cout << "Wing delta ToC: " << Wing_ObjFun_MinToC << ". ";
        cout << "Wing max. twist: " << Wing_MaxTwist << " deg. " << endl;
        if (config_container[ZONE_0]->GetSystemMeasurements() == US)
          cout << "Wing max. curvature: " << Wing_MaxCurvature << " 1/in. ";
        else
          cout << "Wing max. curvature: " << Wing_MaxCurvature << " 1/m. ";
        cout << "Wing max. dihedral: " << Wing_MaxDihedral << " deg." << endl;
      }
    }
  }

  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    geometry_container[ZONE_0]->ComputeAirfoil_Section(
        Plane_P0[iPlane], Plane_Normal[iPlane], -1E6, 1E6, -1E6, 1E6, -1E6, 1E6, nullptr, Xcoord_Airfoil[iPlane],
        Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], Variable_Airfoil[iPlane], true, config_container[ZONE_0]);
  }

  if (rank == MASTER_NODE)
    cout << endl << "-------------------- Objective function evaluation ----------------------" << endl;

  if (rank == MASTER_NODE) {
    /*--- Evaluate objective function ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      if (Xcoord_Airfoil[iPlane].size() > 1) {
        cout << "\nStation " << (iPlane + 1);
        if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << ". XCoord: " << Plane_P0[iPlane][0] << " in, ";
          else
            cout << ". XCoord: " << Plane_P0[iPlane][0] << " m, ";
        }
        if (config_container[ZONE_0]->GetGeo_Description() == WING) {
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << ". YCoord: " << Plane_P0[iPlane][1] << " in, ";
          else
            cout << ". YCoord: " << Plane_P0[iPlane][1] << " m, ";
        }
        if (config_container[ZONE_0]->GetGeo_Description() == TWOD_AIRFOIL) {
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << ". ZCoord: " << Plane_P0[iPlane][2] << " in, ";
          else
            cout << ". ZCoord: " << Plane_P0[iPlane][2] << " m, ";
        }
        if (config_container[ZONE_0]->GetGeo_Description() == NACELLE)
          cout << ". Theta: " << atan2(Plane_Normal[iPlane][1], -Plane_Normal[iPlane][2]) / PI_NUMBER * 180 + 180
               << " deg, ";

        if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
          ObjectiveFunc[0 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Area(
              Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
              Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[1 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Length(
              Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
              Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[2 * nPlane + iPlane] =
              geometry_container[ZONE_0]->Compute_Width(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane],
                                                        Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[3 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_WaterLineWidth(
              Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
              Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[4 * nPlane + iPlane] =
              geometry_container[ZONE_0]->Compute_Height(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane],
                                                         Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "Area: " << ObjectiveFunc[0 * nPlane + iPlane] << " in^2, ";
          else
            cout << "Area: " << ObjectiveFunc[0 * nPlane + iPlane] << " m^2, ";
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "Length: " << ObjectiveFunc[1 * nPlane + iPlane] << " in, ";
          else
            cout << "Length: " << ObjectiveFunc[1 * nPlane + iPlane] << " m, ";
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "Width: " << ObjectiveFunc[2 * nPlane + iPlane] << " in, ";
          else
            cout << "Width: " << ObjectiveFunc[2 * nPlane + iPlane] << " m, ";
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "Waterline width: " << ObjectiveFunc[3 * nPlane + iPlane] << " in, ";
          else
            cout << "Waterline width: " << ObjectiveFunc[3 * nPlane + iPlane] << " m, ";
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "Height: " << ObjectiveFunc[4 * nPlane + iPlane] << " in.";
          else
            cout << "Height: " << ObjectiveFunc[4 * nPlane + iPlane] << " m.";
        } else if (config_container[ZONE_0]->GetGeo_Description() == NACELLE) {
          ObjectiveFunc[0 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Area(
              Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
              Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[1 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_MaxThickness(
              Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
              Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[2 * nPlane + iPlane] =
              geometry_container[ZONE_0]->Compute_Chord(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane],
                                                        Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[3 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_LERadius(
              Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane],
              Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[4 * nPlane + iPlane] = ObjectiveFunc[1 * nPlane + iPlane] / ObjectiveFunc[2 * nPlane + iPlane];
          ObjectiveFunc[5 * nPlane + iPlane] =
              geometry_container[ZONE_0]->Compute_Twist(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane],
                                                        Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "Area: " << ObjectiveFunc[0 * nPlane + iPlane] << " in^2, ";
          else
            cout << "Area: " << ObjectiveFunc[0 * nPlane + iPlane] << " m^2, ";
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "Thickness: " << ObjectiveFunc[1 * nPlane + iPlane] << " in, " << endl;
          else
            cout << "Thickness: " << ObjectiveFunc[1 * nPlane + iPlane] << " m, " << endl;
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "Chord: " << ObjectiveFunc[2 * nPlane + iPlane] << " in, ";
          else
            cout << "Chord: " << ObjectiveFunc[2 * nPlane + iPlane] << " m, ";
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "LE radius: " << ObjectiveFunc[3 * nPlane + iPlane] << " 1/in, ";
          else
            cout << "LE radius: " << ObjectiveFunc[3 * nPlane + iPlane] << " 1/m, ";
          cout << "ToC: " << ObjectiveFunc[4 * nPlane + iPlane] << ", ";
          if (geometry_container[ZONE_0]->GetnDim() == 2)
            cout << "Alpha: " << ObjectiveFunc[5 * nPlane + iPlane] << " deg.";
          else if (geometry_container[ZONE_0]->GetnDim() == 3)
            cout << "Twist angle: " << ObjectiveFunc[5 * nPlane + iPlane] << " deg.";
        } else {
          ObjectiveFunc[0 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Area(
              Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
              Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[1 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_MaxThickness(
              Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
              Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[2 * nPlane + iPlane] =
              geometry_container[ZONE_0]->Compute_Chord(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane],
                                                        Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[3 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_LERadius(
              Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane],
              Zcoord_Airfoil[iPlane]);
          ObjectiveFunc[4 * nPlane + iPlane] = ObjectiveFunc[1 * nPlane + iPlane] / ObjectiveFunc[2 * nPlane + iPlane];
          ObjectiveFunc[5 * nPlane + iPlane] =
              geometry_container[ZONE_0]->Compute_Twist(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane],
                                                        Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "Area: " << ObjectiveFunc[0 * nPlane + iPlane] << " in^2, ";
          else
            cout << "Area: " << ObjectiveFunc[0 * nPlane + iPlane] << " m^2, ";
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "Thickness: " << ObjectiveFunc[1 * nPlane + iPlane] << " in, " << endl;
          else
            cout << "Thickness: " << ObjectiveFunc[1 * nPlane + iPlane] << " m, " << endl;
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "Chord: " << ObjectiveFunc[2 * nPlane + iPlane] << " in, ";
          else
            cout << "Chord: " << ObjectiveFunc[2 * nPlane + iPlane] << " m, ";
          if (config_container[ZONE_0]->GetSystemMeasurements() == US)
            cout << "LE radius: " << ObjectiveFunc[3 * nPlane + iPlane] << " 1/in, ";
          else
            cout << "LE radius: " << ObjectiveFunc[3 * nPlane + iPlane] << " 1/m, ";
          cout << "ToC: " << ObjectiveFunc[4 * nPlane + iPlane] << ", ";
          if (geometry_container[ZONE_0]->GetnDim() == 2)
            cout << "Alpha: " << ObjectiveFunc[5 * nPlane + iPlane] << " deg.";
          else if (geometry_container[ZONE_0]->GetnDim() == 3)
            cout << "Twist angle: " << ObjectiveFunc[5 * nPlane + iPlane] << " deg.";
        }
      }
    }

    /*--- Write the objective function in a external file ---*/
    string filename = config_container[ZONE_0]->GetObjFunc_Value_FileName();
    unsigned short lastindex = filename.find_last_of('.');
    filename = filename.substr(0, lastindex);
    if (tabTecplot)
      filename += ".dat";
    else
      filename += ".csv";
    ObjFunc_file.open(filename.c_str(), ios::out);
    if (tabTecplot) ObjFunc_file << "TITLE = \"SU2_GEO Evaluation\"" << endl;

    if (geometry_container[ZONE_0]->GetnDim() == 2) {
      if (tabTecplot) ObjFunc_file << "VARIABLES =//" << endl;
      ObjFunc_file << "\"AIRFOIL_AREA\",\"AIRFOIL_THICKNESS\",\"AIRFOIL_CHORD\",\"AIRFOIL_LE_RADIUS\",\"AIRFOIL_TOC\","
                      "\"AIRFOIL_ALPHA\"";
    } else if (geometry_container[ZONE_0]->GetnDim() == 3) {
      if (tabTecplot) ObjFunc_file << "VARIABLES = //" << endl;

      if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
        ObjFunc_file << "\"FUSELAGE_VOLUME\",\"FUSELAGE_WETTED_AREA\",\"FUSELAGE_MIN_WIDTH\",\"FUSELAGE_MAX_WIDTH\","
                        "\"FUSELAGE_MIN_WATERLINE_WIDTH\",\"FUSELAGE_MAX_WATERLINE_WIDTH\",\"FUSELAGE_MIN_HEIGHT\","
                        "\"FUSELAGE_MAX_HEIGHT\",\"FUSELAGE_MAX_CURVATURE\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_AREA\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_LENGTH\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_WIDTH\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++)
          ObjFunc_file << "\"STATION" << (iPlane + 1) << "_WATERLINE_WIDTH\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) {
          ObjFunc_file << "\"STATION" << (iPlane + 1) << "_HEIGHT\"";
          if (iPlane != nPlane - 1) ObjFunc_file << ",";
        }
      } else if (config_container[ZONE_0]->GetGeo_Description() == NACELLE) {
        ObjFunc_file << "\"NACELLE_VOLUME\",\"NACELLE_MIN_THICKNESS\",\"NACELLE_MAX_THICKNESS\",\"NACELLE_MIN_CHORD\","
                        "\"NACELLE_MAX_CHORD\",\"NACELLE_MIN_LE_RADIUS\",\"NACELLE_MAX_LE_RADIUS\",\"NACELLE_MIN_TOC\","
                        "\"NACELLE_MAX_TOC\",\"NACELLE_OBJFUN_MIN_TOC\",\"NACELLE_MAX_TWIST\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_AREA\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_THICKNESS\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_CHORD\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_LE_RADIUS\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_TOC\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) {
          ObjFunc_file << "\"STATION" << (iPlane + 1) << "_TWIST\"";
          if (iPlane != nPlane - 1) ObjFunc_file << ",";
        }
      } else {
        ObjFunc_file << "\"WING_VOLUME\",\"WING_MIN_THICKNESS\",\"WING_MAX_THICKNESS\",\"WING_MIN_CHORD\",\"WING_MAX_"
                        "CHORD\",\"WING_MIN_LE_RADIUS\",\"WING_MAX_LE_RADIUS\",\"WING_MIN_TOC\",\"WING_MAX_TOC\","
                        "\"WING_OBJFUN_MIN_TOC\",\"WING_MAX_TWIST\",\"WING_MAX_CURVATURE\",\"WING_MAX_DIHEDRAL\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_AREA\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_THICKNESS\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_CHORD\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_LE_RADIUS\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"STATION" << (iPlane + 1) << "_TOC\",";
        for (iPlane = 0; iPlane < nPlane; iPlane++) {
          ObjFunc_file << "\"STATION" << (iPlane + 1) << "_TWIST\"";
          if (iPlane != nPlane - 1) ObjFunc_file << ",";
        }
      }
    }

    if (tabTecplot)
      ObjFunc_file << "\nZONE T= \"Geometrical variables (value)\"" << endl;
    else
      ObjFunc_file << endl;

    if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
      if (geometry_container[ZONE_0]->GetnDim() == 3) {
        ObjFunc_file << Fuselage_Volume << ", " << Fuselage_WettedArea << ", " << Fuselage_MinWidth << ", "
                     << Fuselage_MaxWidth << ", " << Fuselage_MinWaterLineWidth << ", " << Fuselage_MaxWaterLineWidth
                     << ", " << Fuselage_MinHeight << ", " << Fuselage_MaxHeight << ", " << Fuselage_MaxCurvature
                     << ", ";
      }
      for (iPlane = 0; iPlane < nPlane * 5; iPlane++) {
        ObjFunc_file << ObjectiveFunc[iPlane];
        if (iPlane != (nPlane * 5) - 1) ObjFunc_file << ", ";
      }
    } else if (config_container[ZONE_0]->GetGeo_Description() == NACELLE) {
      if (geometry_container[ZONE_0]->GetnDim() == 3) {
        ObjFunc_file << Nacelle_Volume << ", " << Nacelle_MinThickness << ", " << Nacelle_MaxThickness << ", "
                     << Nacelle_MinChord << ", " << Nacelle_MaxChord << ", " << Nacelle_MinLERadius << ", "
                     << Nacelle_MaxLERadius << ", " << Nacelle_MinToC << ", " << Nacelle_MaxToC << ", "
                     << Nacelle_ObjFun_MinToC << ", " << Nacelle_MaxTwist << ", ";
      }
      for (iPlane = 0; iPlane < nPlane * 6; iPlane++) {
        ObjFunc_file << ObjectiveFunc[iPlane];
        if (iPlane != (nPlane * 6) - 1) ObjFunc_file << ", ";
      }
    } else {
      if (geometry_container[ZONE_0]->GetnDim() == 3) {
        ObjFunc_file << Wing_Volume << ", " << Wing_MinThickness << ", " << Wing_MaxThickness << ", " << Wing_MinChord
                     << ", " << Wing_MaxChord << ", " << Wing_MinLERadius << ", " << Wing_MaxLERadius << ", "
                     << Wing_MinToC << ", " << Wing_MaxToC << ", " << Wing_ObjFun_MinToC << ", " << Wing_MaxTwist
                     << ", " << Wing_MaxCurvature << ", " << Wing_MaxDihedral << ", ";
      }
      for (iPlane = 0; iPlane < nPlane * 6; iPlane++) {
        ObjFunc_file << ObjectiveFunc[iPlane];
        if (iPlane != (nPlane * 6) - 1) ObjFunc_file << ", ";
      }
    }

    ObjFunc_file.close();
  }

  if (config_container[ZONE_0]->GetGeometryMode() == GRADIENT) {
    /*--- Definition of the Class for surface deformation ---*/
    surface_movement = new CSurfaceMovement();

    /*--- Copy coordinates to the surface structure ---*/
    surface_movement->CopyBoundary(geometry_container[ZONE_0], config_container[ZONE_0]);

    /*--- Definition of the FFD deformation class ---*/
    FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];
    for (iFFDBox = 0; iFFDBox < MAX_NUMBER_FFD; iFFDBox++) FFDBox[iFFDBox] = nullptr;

    if (rank == MASTER_NODE)
      cout << endl << endl << "------------- Gradient evaluation using finite differences --------------" << endl;

    /*--- Write the gradient in a external file ---*/
    if (rank == MASTER_NODE) {
      string filename = config_container[ZONE_0]->GetObjFunc_Grad_FileName();
      unsigned short lastindex = filename.find_last_of('.');
      filename = filename.substr(0, lastindex);
      if (tabTecplot)
        filename += ".dat";
      else
        filename += ".csv";
      Gradient_file.open(filename.c_str(), ios::out);
    }

    for (iDV = 0; iDV < config_container[ZONE_0]->GetnDV(); iDV++) {
      /*--- Free Form deformation based ---*/

      if ((config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_CONTROL_POINT_2D) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_CAMBER_2D) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_THICKNESS_2D) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_CONTROL_POINT) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_NACELLE) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_GULL) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_TWIST) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_ROTATION) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_CAMBER) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_THICKNESS)) {
        /*--- Read the FFD information in the first iteration ---*/

        if (iDV == 0) {
          if (rank == MASTER_NODE) cout << "Read the FFD information from mesh file." << endl;

          /*--- Read the FFD information from the grid file ---*/

          surface_movement->ReadFFDInfo(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox,
                                        config_container[ZONE_0]->GetMesh_FileName());

          /*--- Modify the control points for polar based computations ---*/

          if (config_container[ZONE_0]->GetFFD_CoordSystem() == CYLINDRICAL) {
            for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetCart2Cyl_ControlPoints(config_container[ZONE_0]);
            }
          } else if (config_container[ZONE_0]->GetFFD_CoordSystem() == SPHERICAL) {
            for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetCart2Sphe_ControlPoints(config_container[ZONE_0]);
            }
          } else if (config_container[ZONE_0]->GetFFD_CoordSystem() == POLAR) {
            for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetCart2Sphe_ControlPoints(config_container[ZONE_0]);
            }
          }

          /*--- If the FFDBox was not defined in the input file ---*/

          if (!surface_movement->GetFFDBoxDefinition()) {
            SU2_MPI::Error("The input grid doesn't have the entire FFD information!", CURRENT_FUNCTION);
          }

          for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
            if (rank == MASTER_NODE) cout << "Checking FFD box dimension." << endl;
            surface_movement->CheckFFDDimension(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox],
                                                iFFDBox);

            if (rank == MASTER_NODE) cout << "Check the FFD box intersections with the solid surfaces." << endl;
            surface_movement->CheckFFDIntersections(geometry_container[ZONE_0], config_container[ZONE_0],
                                                    FFDBox[iFFDBox], iFFDBox);
          }

          if (rank == MASTER_NODE)
            cout << "-------------------------------------------------------------------------" << endl;
        }

        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number " << iDV << "." << endl;
          cout << "Perform 3D deformation of the surface." << endl;
        }

        /*--- Apply the control point change ---*/

        MoveSurface = false;

        for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
          Local_MoveSurface = false;

          switch (config_container[ZONE_0]->GetDesign_Variable(iDV)) {
            case FFD_CONTROL_POINT_2D:
              Local_MoveSurface = surface_movement->SetFFDCPChange_2D(
                  geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true);
              break;
            case FFD_CAMBER_2D:
              Local_MoveSurface = surface_movement->SetFFDCamber_2D(
                  geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true);
              break;
            case FFD_THICKNESS_2D:
              Local_MoveSurface = surface_movement->SetFFDThickness_2D(
                  geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true);
              break;
            case FFD_CONTROL_POINT:
              Local_MoveSurface = surface_movement->SetFFDCPChange(geometry_container[ZONE_0], config_container[ZONE_0],
                                                                   FFDBox[iFFDBox], FFDBox, iDV, true);
              break;
            case FFD_NACELLE:
              Local_MoveSurface = surface_movement->SetFFDNacelle(geometry_container[ZONE_0], config_container[ZONE_0],
                                                                  FFDBox[iFFDBox], FFDBox, iDV, true);
              break;
            case FFD_GULL:
              Local_MoveSurface = surface_movement->SetFFDGull(geometry_container[ZONE_0], config_container[ZONE_0],
                                                               FFDBox[iFFDBox], FFDBox, iDV, true);
              break;
            case FFD_TWIST:
              Local_MoveSurface = surface_movement->SetFFDTwist(geometry_container[ZONE_0], config_container[ZONE_0],
                                                                FFDBox[iFFDBox], FFDBox, iDV, true);
              break;
            case FFD_ROTATION:
              Local_MoveSurface = surface_movement->SetFFDRotation(geometry_container[ZONE_0], config_container[ZONE_0],
                                                                   FFDBox[iFFDBox], FFDBox, iDV, true);
              break;
            case FFD_CAMBER:
              Local_MoveSurface = surface_movement->SetFFDCamber(geometry_container[ZONE_0], config_container[ZONE_0],
                                                                 FFDBox[iFFDBox], FFDBox, iDV, true);
              break;
            case FFD_THICKNESS:
              Local_MoveSurface = surface_movement->SetFFDThickness(
                  geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true);
              break;
            case FFD_CONTROL_SURFACE:
              Local_MoveSurface = surface_movement->SetFFDControl_Surface(
                  geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true);
              break;
          }

          /*--- Recompute cartesian coordinates using the new control points position ---*/

          if (Local_MoveSurface) {
            MoveSurface = true;
            surface_movement->SetCartesianCoord(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox],
                                                iFFDBox, true);
          }
        }

      }

      /*--- Hicks Henne design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == HICKS_HENNE) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number " << iDV << "." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetHicksHenne(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- Surface bump design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == SURFACE_BUMP) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number " << iDV << "." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetSurface_Bump(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- CST design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == CST) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number " << iDV << "." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetCST(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- Translation design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == TRANSLATION) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number " << iDV << "." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetTranslation(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- Scale design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == SCALE) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number " << iDV << "." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetScale(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- Rotation design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == ROTATION) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number " << iDV << "." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetRotation(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- NACA_4Digits design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == NACA_4DIGITS) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number " << iDV << "." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetNACA_4Digits(geometry_container[ZONE_0], config_container[ZONE_0]);
      }

      /*--- Parabolic design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == PARABOLIC) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number " << iDV << "." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetParabolic(geometry_container[ZONE_0], config_container[ZONE_0]);
      }

      /*--- Design variable not implement ---*/

      else {
        if (rank == MASTER_NODE) cout << "Design Variable not implemented yet" << endl;
      }

      if (MoveSurface) {
        /*--- Compute the gradient for the volume. In 2D this is just
         the gradient of the area. ---*/

        if (geometry_container[ZONE_0]->GetnDim() == 3) {
          if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
            geometry_container[ZONE_0]->Compute_Fuselage(
                config_container[ZONE_0], false, Fuselage_Volume_New, Fuselage_WettedArea_New, Fuselage_MinWidth_New,
                Fuselage_MaxWidth_New, Fuselage_MinWaterLineWidth_New, Fuselage_MaxWaterLineWidth_New,
                Fuselage_MinHeight_New, Fuselage_MaxHeight_New, Fuselage_MaxCurvature_New);
          } else if (config_container[ZONE_0]->GetGeo_Description() == NACELLE) {
            geometry_container[ZONE_0]->Compute_Nacelle(
                config_container[ZONE_0], false, Nacelle_Volume_New, Nacelle_MinThickness_New, Nacelle_MaxThickness_New,
                Nacelle_MinChord_New, Nacelle_MaxChord_New, Nacelle_MinLERadius_New, Nacelle_MaxLERadius_New,
                Nacelle_MinToC_New, Nacelle_MaxToC_New, Nacelle_ObjFun_MinToC_New, Nacelle_MaxTwist_New);
          } else {
            geometry_container[ZONE_0]->Compute_Wing(config_container[ZONE_0], false, Wing_Volume_New,
                                                     Wing_MinThickness_New, Wing_MaxThickness_New, Wing_MinChord_New,
                                                     Wing_MaxChord_New, Wing_MinLERadius_New, Wing_MaxLERadius_New,
                                                     Wing_MinToC_New, Wing_MaxToC_New, Wing_ObjFun_MinToC_New,
                                                     Wing_MaxTwist_New, Wing_MaxCurvature_New, Wing_MaxDihedral_New);
          }
        }

        /*--- Create airfoil structure ---*/

        for (iPlane = 0; iPlane < nPlane; iPlane++) {
          geometry_container[ZONE_0]->ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane], -1E6, 1E6, -1E6,
                                                             1E6, -1E6, 1E6, nullptr, Xcoord_Airfoil[iPlane],
                                                             Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane],
                                                             Variable_Airfoil[iPlane], false, config_container[ZONE_0]);
        }
      }

      /*--- Compute gradient ---*/

      if (rank == MASTER_NODE) {
        delta_eps = config_container[ZONE_0]->GetDV_Value(iDV);

        if (delta_eps == 0) {
          SU2_MPI::Error("The finite difference steps is zero!!", CURRENT_FUNCTION);
        }

        if (MoveSurface) {
          if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
            Fuselage_Volume_Grad = (Fuselage_Volume_New - Fuselage_Volume) / delta_eps;
            Fuselage_WettedArea_Grad = (Fuselage_WettedArea_New - Fuselage_WettedArea) / delta_eps;
            Fuselage_MinWidth_Grad = (Fuselage_MinWidth_New - Fuselage_MinWidth) / delta_eps;
            Fuselage_MaxWidth_Grad = (Fuselage_MaxWidth_New - Fuselage_MaxWidth) / delta_eps;
            Fuselage_MinWaterLineWidth_Grad = (Fuselage_MinWaterLineWidth_New - Fuselage_MinWaterLineWidth) / delta_eps;
            Fuselage_MaxWaterLineWidth_Grad = (Fuselage_MaxWaterLineWidth_New - Fuselage_MaxWaterLineWidth) / delta_eps;
            Fuselage_MinHeight_Grad = (Fuselage_MinHeight_New - Fuselage_MinHeight) / delta_eps;
            Fuselage_MaxHeight_Grad = (Fuselage_MaxHeight_New - Fuselage_MaxHeight) / delta_eps;
            Fuselage_MaxCurvature_Grad = (Fuselage_MaxCurvature_New - Fuselage_MaxCurvature) / delta_eps;

          } else if (config_container[ZONE_0]->GetGeo_Description() == NACELLE) {
            Nacelle_Volume_Grad = (Nacelle_Volume_New - Nacelle_Volume) / delta_eps;
            Nacelle_MinThickness_Grad = (Nacelle_MinThickness_New - Nacelle_MinThickness) / delta_eps;
            Nacelle_MaxThickness_Grad = (Nacelle_MaxThickness_New - Nacelle_MaxThickness) / delta_eps;
            Nacelle_MinChord_Grad = (Nacelle_MinChord_New - Nacelle_MinChord) / delta_eps;
            Nacelle_MaxChord_Grad = (Nacelle_MaxChord_New - Nacelle_MaxChord) / delta_eps;
            Nacelle_MinLERadius_Grad = (Nacelle_MinLERadius_New - Nacelle_MinLERadius) / delta_eps;
            Nacelle_MaxLERadius_Grad = (Nacelle_MaxLERadius_New - Nacelle_MaxLERadius) / delta_eps;
            Nacelle_MinToC_Grad = (Nacelle_MinToC_New - Nacelle_MinToC) / delta_eps;
            Nacelle_MaxToC_Grad = (Nacelle_MaxToC_New - Nacelle_MaxToC) / delta_eps;
            Nacelle_ObjFun_MinToC_Grad = (Nacelle_ObjFun_MinToC_New - Nacelle_ObjFun_MinToC) / delta_eps;
            Nacelle_MaxTwist_Grad = (Nacelle_MaxTwist_New - Nacelle_MaxTwist) / delta_eps;
          } else {
            Wing_Volume_Grad = (Wing_Volume_New - Wing_Volume) / delta_eps;
            Wing_MinThickness_Grad = (Wing_MinThickness_New - Wing_MinThickness) / delta_eps;
            Wing_MaxThickness_Grad = (Wing_MaxThickness_New - Wing_MaxThickness) / delta_eps;
            Wing_MinChord_Grad = (Wing_MinChord_New - Wing_MinChord) / delta_eps;
            Wing_MaxChord_Grad = (Wing_MaxChord_New - Wing_MaxChord) / delta_eps;
            Wing_MinLERadius_Grad = (Wing_MinLERadius_New - Wing_MinLERadius) / delta_eps;
            Wing_MaxLERadius_Grad = (Wing_MaxLERadius_New - Wing_MaxLERadius) / delta_eps;
            Wing_MinToC_Grad = (Wing_MinToC_New - Wing_MinToC) / delta_eps;
            Wing_MaxToC_Grad = (Wing_MaxToC_New - Wing_MaxToC) / delta_eps;
            Wing_ObjFun_MinToC_Grad = (Wing_ObjFun_MinToC_New - Wing_ObjFun_MinToC) / delta_eps;
            Wing_MaxTwist_Grad = (Wing_MaxTwist_New - Wing_MaxTwist) / delta_eps;
            Wing_MaxCurvature_Grad = (Wing_MaxCurvature_New - Wing_MaxCurvature) / delta_eps;
            Wing_MaxDihedral_Grad = (Wing_MaxDihedral_New - Wing_MaxDihedral) / delta_eps;
          }

          for (iPlane = 0; iPlane < nPlane; iPlane++) {
            if (Xcoord_Airfoil[iPlane].size() > 1) {
              if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
                ObjectiveFunc_New[0 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Area(
                    Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
                    Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
                Gradient[0 * nPlane + iPlane] =
                    (ObjectiveFunc_New[0 * nPlane + iPlane] - ObjectiveFunc[0 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[1 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Length(
                    Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
                    Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
                Gradient[1 * nPlane + iPlane] =
                    (ObjectiveFunc_New[1 * nPlane + iPlane] - ObjectiveFunc[1 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[2 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Width(
                    Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane],
                    Zcoord_Airfoil[iPlane]);
                Gradient[2 * nPlane + iPlane] =
                    (ObjectiveFunc_New[2 * nPlane + iPlane] - ObjectiveFunc[2 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[3 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_WaterLineWidth(
                    Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
                    Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
                Gradient[3 * nPlane + iPlane] =
                    (ObjectiveFunc_New[3 * nPlane + iPlane] - ObjectiveFunc[3 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[4 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Height(
                    Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane],
                    Zcoord_Airfoil[iPlane]);
                Gradient[4 * nPlane + iPlane] =
                    (ObjectiveFunc_New[4 * nPlane + iPlane] - ObjectiveFunc[4 * nPlane + iPlane]) / delta_eps;

              }

              else if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
                ObjectiveFunc_New[0 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Area(
                    Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
                    Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
                Gradient[0 * nPlane + iPlane] =
                    (ObjectiveFunc_New[0 * nPlane + iPlane] - ObjectiveFunc[0 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[1 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_MaxThickness(
                    Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
                    Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
                Gradient[1 * nPlane + iPlane] =
                    (ObjectiveFunc_New[1 * nPlane + iPlane] - ObjectiveFunc[1 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[2 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Chord(
                    Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane],
                    Zcoord_Airfoil[iPlane]);
                Gradient[2 * nPlane + iPlane] =
                    (ObjectiveFunc_New[2 * nPlane + iPlane] - ObjectiveFunc[2 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[3 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_LERadius(
                    Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane],
                    Zcoord_Airfoil[iPlane]);
                Gradient[3 * nPlane + iPlane] =
                    (ObjectiveFunc_New[3 * nPlane + iPlane] - ObjectiveFunc[3 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[4 * nPlane + iPlane] =
                    ObjectiveFunc_New[1 * nPlane + iPlane] / ObjectiveFunc_New[2 * nPlane + iPlane];
                Gradient[4 * nPlane + iPlane] =
                    (ObjectiveFunc_New[4 * nPlane + iPlane] - ObjectiveFunc[4 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[5 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Twist(
                    Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane],
                    Zcoord_Airfoil[iPlane]);
                Gradient[5 * nPlane + iPlane] =
                    (ObjectiveFunc_New[5 * nPlane + iPlane] - ObjectiveFunc[5 * nPlane + iPlane]) / delta_eps;

              }

              else {
                ObjectiveFunc_New[0 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Area(
                    Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
                    Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
                Gradient[0 * nPlane + iPlane] =
                    (ObjectiveFunc_New[0 * nPlane + iPlane] - ObjectiveFunc[0 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[1 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_MaxThickness(
                    Plane_P0[iPlane], Plane_Normal[iPlane], config_container[ZONE_0], Xcoord_Airfoil[iPlane],
                    Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
                Gradient[1 * nPlane + iPlane] =
                    (ObjectiveFunc_New[1 * nPlane + iPlane] - ObjectiveFunc[1 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[2 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Chord(
                    Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane],
                    Zcoord_Airfoil[iPlane]);
                Gradient[2 * nPlane + iPlane] =
                    (ObjectiveFunc_New[2 * nPlane + iPlane] - ObjectiveFunc[2 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[3 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_LERadius(
                    Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane],
                    Zcoord_Airfoil[iPlane]);
                Gradient[3 * nPlane + iPlane] =
                    (ObjectiveFunc_New[3 * nPlane + iPlane] - ObjectiveFunc[3 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[4 * nPlane + iPlane] =
                    ObjectiveFunc_New[1 * nPlane + iPlane] / ObjectiveFunc_New[2 * nPlane + iPlane];
                Gradient[4 * nPlane + iPlane] =
                    (ObjectiveFunc_New[4 * nPlane + iPlane] - ObjectiveFunc[4 * nPlane + iPlane]) / delta_eps;

                ObjectiveFunc_New[5 * nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Twist(
                    Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane],
                    Zcoord_Airfoil[iPlane]);
                Gradient[5 * nPlane + iPlane] =
                    (ObjectiveFunc_New[5 * nPlane + iPlane] - ObjectiveFunc[5 * nPlane + iPlane]) / delta_eps;
              }
            }
          }

        }

        else {
          if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
            Fuselage_Volume_Grad = 0.0;
            Fuselage_WettedArea_Grad = 0.0;
            Fuselage_MinWidth_Grad = 0.0;
            Fuselage_MaxWidth_Grad = 0.0;
            Fuselage_MinWaterLineWidth_Grad = 0.0;
            Fuselage_MaxWaterLineWidth_Grad = 0.0;
            Fuselage_MinHeight_Grad = 0.0;
            Fuselage_MaxHeight_Grad = 0.0;
            Fuselage_MaxCurvature_Grad = 0.0;

            for (iPlane = 0; iPlane < nPlane; iPlane++) {
              Gradient[0 * nPlane + iPlane] = 0.0;
              Gradient[1 * nPlane + iPlane] = 0.0;
              Gradient[2 * nPlane + iPlane] = 0.0;
              Gradient[3 * nPlane + iPlane] = 0.0;
              Gradient[4 * nPlane + iPlane] = 0.0;
            }

          } else if (config_container[ZONE_0]->GetGeo_Description() == NACELLE) {
            Nacelle_Volume_Grad = 0.0;
            Nacelle_MinThickness_Grad = 0.0;
            Nacelle_MaxThickness_Grad = 0.0;
            Nacelle_MinChord_Grad = 0.0;
            Nacelle_MaxChord_Grad = 0.0;
            Nacelle_MinLERadius_Grad = 0.0;
            Nacelle_MaxLERadius_Grad = 0.0;
            Nacelle_MinToC_Grad = 0.0;
            Nacelle_MaxToC_Grad = 0.0;
            Nacelle_ObjFun_MinToC_Grad = 0.0;
            Nacelle_MaxTwist_Grad = 0.0;

            for (iPlane = 0; iPlane < nPlane; iPlane++) {
              Gradient[0 * nPlane + iPlane] = 0.0;
              Gradient[1 * nPlane + iPlane] = 0.0;
              Gradient[2 * nPlane + iPlane] = 0.0;
              Gradient[3 * nPlane + iPlane] = 0.0;
              Gradient[4 * nPlane + iPlane] = 0.0;
              Gradient[5 * nPlane + iPlane] = 0.0;
            }
          } else {
            Wing_Volume_Grad = 0.0;
            Wing_MinThickness_Grad = 0.0;
            Wing_MaxThickness_Grad = 0.0;
            Wing_MinChord_Grad = 0.0;
            Wing_MaxChord_Grad = 0.0;
            Wing_MinLERadius_Grad = 0.0;
            Wing_MaxLERadius_Grad = 0.0;
            Wing_MinToC_Grad = 0.0;
            Wing_MaxToC_Grad = 0.0;
            Wing_ObjFun_MinToC_Grad = 0.0;
            Wing_MaxTwist_Grad = 0.0;
            Wing_MaxCurvature_Grad = 0.0;
            Wing_MaxDihedral_Grad = 0.0;

            for (iPlane = 0; iPlane < nPlane; iPlane++) {
              Gradient[0 * nPlane + iPlane] = 0.0;
              Gradient[1 * nPlane + iPlane] = 0.0;
              Gradient[2 * nPlane + iPlane] = 0.0;
              Gradient[3 * nPlane + iPlane] = 0.0;
              Gradient[4 * nPlane + iPlane] = 0.0;
              Gradient[5 * nPlane + iPlane] = 0.0;
            }
          }
        }

        /*--- Screen output ---*/

        if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
          if (geometry_container[ZONE_0]->GetnDim() == 3) {
            cout << "\nFuselage volume grad.: " << Fuselage_Volume_Grad << ". ";
            cout << "Fuselage wetted area grad.: " << Fuselage_WettedArea_Grad << ". ";
            cout << "Fuselage min. width grad.: " << Fuselage_MinWidth_Grad << ". ";
            cout << "Fuselage max. width grad.: " << Fuselage_MaxWidth_Grad << "." << endl;
            cout << "Fuselage min. waterline width grad.: " << Fuselage_MinWaterLineWidth_Grad << ". ";
            cout << "Fuselage max. waterline width grad.: " << Fuselage_MaxWaterLineWidth_Grad << "." << endl;
            cout << "Fuselage min. height grad.: " << Fuselage_MinHeight_Grad << ". ";
            cout << "Fuselage max. height grad.: " << Fuselage_MaxHeight_Grad << ". ";
            cout << "Fuselage max. curv. grad.: " << Fuselage_MaxCurvature_Grad << ".";
          }

          for (iPlane = 0; iPlane < nPlane; iPlane++) {
            if (Xcoord_Airfoil[iPlane].size() > 1) {
              cout << "\nStation " << (iPlane + 1) << ". XCoord: " << Plane_P0[iPlane][0] << ". ";
              cout << "Area grad.: " << Gradient[0 * nPlane + iPlane] << ". ";
              cout << "Length grad.: " << Gradient[1 * nPlane + iPlane] << ". ";
              cout << "Width grad.: " << Gradient[2 * nPlane + iPlane] << ". ";
              cout << "Waterline width grad.: " << Gradient[3 * nPlane + iPlane] << ". ";
              cout << "Height grad.: " << Gradient[4 * nPlane + iPlane] << ". ";
            }
          }
        } else if (config_container[ZONE_0]->GetGeo_Description() == NACELLE) {
          if (geometry_container[ZONE_0]->GetnDim() == 3) {
            cout << "\nNacelle volume grad.: " << Nacelle_Volume_Grad << ". ";
            cout << "Nacelle min. thickness grad.: " << Nacelle_MinThickness_Grad << ". ";
            cout << "Nacelle max. thickness grad.: " << Nacelle_MaxThickness_Grad << ". ";
            cout << "Nacelle min. chord grad.: " << Nacelle_MinChord_Grad << ". ";
            cout << "Nacelle max. chord grad.: " << Nacelle_MaxChord_Grad << "." << endl;
            cout << "Nacelle min. LE radius grad.: " << Nacelle_MinChord_Grad << ". ";
            cout << "Nacelle max. LE radius grad.: " << Nacelle_MaxChord_Grad << ". ";
            cout << "Nacelle min. ToC grad.: " << Nacelle_MinToC_Grad << ". ";
            cout << "Nacelle max. ToC grad.: " << Nacelle_MaxToC_Grad << ". ";
            cout << "Nacelle delta ToC grad.: " << Nacelle_ObjFun_MinToC_Grad << "." << endl;
            cout << "Nacelle max. twist grad.: " << Nacelle_MaxTwist_Grad << ". ";
          }

          for (iPlane = 0; iPlane < nPlane; iPlane++) {
            if (Xcoord_Airfoil[iPlane].size() > 1) {
              cout << "\nStation " << (iPlane + 1) << ". YCoord: " << Plane_P0[iPlane][1] << ". ";
              cout << "Area grad.: " << Gradient[0 * nPlane + iPlane] << ". ";
              cout << "Thickness grad.: " << Gradient[1 * nPlane + iPlane] << ". ";
              cout << "Chord grad.: " << Gradient[2 * nPlane + iPlane] << ". ";
              cout << "LE radius grad.: " << Gradient[3 * nPlane + iPlane] << ". ";
              cout << "ToC grad.: " << Gradient[4 * nPlane + iPlane] << ". ";
              cout << "Twist angle grad.: " << Gradient[5 * nPlane + iPlane] << ". ";
            }
          }
        } else {
          if (geometry_container[ZONE_0]->GetnDim() == 3) {
            cout << "\nWing volume grad.: " << Wing_Volume_Grad << ". ";
            cout << "Wing min. thickness grad.: " << Wing_MinThickness_Grad << ". ";
            cout << "Wing max. thickness grad.: " << Wing_MaxThickness_Grad << ". ";
            cout << "Wing min. chord grad.: " << Wing_MinChord_Grad << ". ";
            cout << "Wing max. chord grad.: " << Wing_MaxChord_Grad << "." << endl;
            cout << "Wing min. LE radius grad.: " << Wing_MinChord_Grad << ". ";
            cout << "Wing max. LE radius grad.: " << Wing_MaxChord_Grad << ". ";
            cout << "Wing min. ToC grad.: " << Wing_MinToC_Grad << ". ";
            cout << "Wing max. ToC grad.: " << Wing_MaxToC_Grad << ". ";
            cout << "Wing delta ToC grad.: " << Wing_ObjFun_MinToC_Grad << "." << endl;
            cout << "Wing max. twist grad.: " << Wing_MaxTwist_Grad << ". ";
            cout << "Wing max. curv. grad.: " << Wing_MaxCurvature_Grad << ". ";
            cout << "Wing max. dihedral grad.: " << Wing_MaxDihedral_Grad << "." << endl;
          }

          for (iPlane = 0; iPlane < nPlane; iPlane++) {
            if (Xcoord_Airfoil[iPlane].size() > 1) {
              cout << "\nStation " << (iPlane + 1) << ". YCoord: " << Plane_P0[iPlane][1] << ". ";
              cout << "Area grad.: " << Gradient[0 * nPlane + iPlane] << ". ";
              cout << "Thickness grad.: " << Gradient[1 * nPlane + iPlane] << ". ";
              cout << "Chord grad.: " << Gradient[2 * nPlane + iPlane] << ". ";
              cout << "LE radius grad.: " << Gradient[3 * nPlane + iPlane] << ". ";
              cout << "ToC grad.: " << Gradient[4 * nPlane + iPlane] << ". ";
              cout << "Twist angle grad.: " << Gradient[5 * nPlane + iPlane] << ". ";
            }
          }
        }

        cout << endl;

        if (iDV == 0) {
          if (tabTecplot) Gradient_file << "TITLE = \"SU2_GEO Gradient\"" << endl;
          if (tabTecplot) Gradient_file << "VARIABLES = //" << endl;

          if (geometry_container[ZONE_0]->GetnDim() == 2) {
            Gradient_file << "\"DESIGN_VARIABLE\",\"AIRFOIL_AREA\",\"AIRFOIL_THICKNESS\",\"AIRFOIL_CHORD\",\"AIRFOIL_"
                             "LE_RADIUS\",\"AIRFOIL_TOC\",\"AIRFOIL_ALPHA\"";
          } else if (geometry_container[ZONE_0]->GetnDim() == 3) {
            if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
              Gradient_file << "\"DESIGN_VARIABLE\",";
              Gradient_file << "\"FUSELAGE_VOLUME\",\"FUSELAGE_WETTED_AREA\",\"FUSELAGE_MIN_WIDTH\",\"FUSELAGE_MAX_"
                               "WIDTH\",\"FUSELAGE_MIN_WATERLINE_WIDTH\",\"FUSELAGE_MAX_WATERLINE_WIDTH\",\"FUSELAGE_"
                               "MIN_HEIGHT\",\"FUSELAGE_MAX_HEIGHT\",\"FUSELAGE_MAX_CURVATURE\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"STATION" << (iPlane + 1) << "_AREA\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"STATION" << (iPlane + 1) << "_LENGTH\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"STATION" << (iPlane + 1) << "_WIDTH\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++)
                Gradient_file << "\"STATION" << (iPlane + 1) << "_WATERLINE_WIDTH\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) {
                Gradient_file << "\"STATION" << (iPlane + 1) << "_HEIGHT\"";
                if (iPlane != nPlane - 1) Gradient_file << ",";
              }
            } else if (config_container[ZONE_0]->GetGeo_Description() == NACELLE) {
              Gradient_file << "\"DESIGN_VARIABLE\",";
              Gradient_file
                  << "\"NACELLE_VOLUME\",\"NACELLE_MIN_THICKNESS\",\"NACELLE_MAX_THICKNESS\",\"NACELLE_MIN_CHORD\","
                     "\"NACELLE_MAX_CHORD\",\"NACELLE_MIN_LE_RADIUS\",\"NACELLE_MAX_LE_RADIUS\",\"NACELLE_MIN_TOC\","
                     "\"NACELLE_MAX_TOC\",\"NACELLE_OBJFUN_MIN_TOC\",\"NACELLE_MAX_TWIST\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"STATION" << (iPlane + 1) << "_AREA\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++)
                Gradient_file << "\"STATION" << (iPlane + 1) << "_THICKNESS\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"STATION" << (iPlane + 1) << "_CHORD\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++)
                Gradient_file << "\"STATION" << (iPlane + 1) << "_LE_RADIUS\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"STATION" << (iPlane + 1) << "_TOC\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) {
                Gradient_file << "\"STATION" << (iPlane + 1) << "_TWIST\"";
                if (iPlane != nPlane - 1) Gradient_file << ",";
              }
            } else {
              Gradient_file << "\"DESIGN_VARIABLE\",";
              Gradient_file
                  << "\"WING_VOLUME\",\"WING_MIN_THICKNESS\",\"WING_MAX_THICKNESS\",\"WING_MIN_CHORD\",\"WING_MAX_"
                     "CHORD\",\"WING_MIN_LE_RADIUS\",\"WING_MAX_LE_RADIUS\",\"WING_MIN_TOC\",\"WING_MAX_TOC\",\"WING_"
                     "OBJFUN_MIN_TOC\",\"WING_MAX_TWIST\",\"WING_MAX_CURVATURE\",\"WING_MAX_DIHEDRAL\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"STATION" << (iPlane + 1) << "_AREA\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++)
                Gradient_file << "\"STATION" << (iPlane + 1) << "_THICKNESS\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"STATION" << (iPlane + 1) << "_CHORD\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++)
                Gradient_file << "\"STATION" << (iPlane + 1) << "_LE_RADIUS\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"STATION" << (iPlane + 1) << "_TOC\",";
              for (iPlane = 0; iPlane < nPlane; iPlane++) {
                Gradient_file << "\"STATION" << (iPlane + 1) << "_TWIST\"";
                if (iPlane != nPlane - 1) Gradient_file << ",";
              }
            }
          }

          if (tabTecplot)
            Gradient_file << "\nZONE T= \"Geometrical variables (gradient)\"" << endl;
          else
            Gradient_file << endl;
        }

        Gradient_file << (iDV) << ",";

        if (config_container[ZONE_0]->GetGeo_Description() == FUSELAGE) {
          if (geometry_container[ZONE_0]->GetnDim() == 3) {
            Gradient_file << Fuselage_Volume_Grad << "," << Fuselage_WettedArea_Grad << "," << Fuselage_MinWidth_Grad
                          << "," << Fuselage_MaxWidth_Grad << "," << Fuselage_MinWaterLineWidth_Grad << ","
                          << Fuselage_MaxWaterLineWidth_Grad << "," << Fuselage_MinHeight_Grad << ","
                          << Fuselage_MaxHeight_Grad << "," << Fuselage_MaxCurvature_Grad << ",";
          }
          for (iPlane = 0; iPlane < nPlane * 5; iPlane++) {
            Gradient_file << Gradient[iPlane];
            if (iPlane != (nPlane * 5) - 1) Gradient_file << ",";
          }
        } else if (config_container[ZONE_0]->GetGeo_Description() == NACELLE) {
          if (geometry_container[ZONE_0]->GetnDim() == 3) {
            Gradient_file << Nacelle_Volume_Grad << "," << Nacelle_MinThickness_Grad << "," << Nacelle_MaxThickness_Grad
                          << "," << Nacelle_MinChord_Grad << "," << Nacelle_MaxChord_Grad << ","
                          << Nacelle_MinLERadius_Grad << "," << Nacelle_MaxLERadius_Grad << "," << Nacelle_MinToC_Grad
                          << "," << Nacelle_MaxToC_Grad << "," << Nacelle_ObjFun_MinToC_Grad << ","
                          << Nacelle_MaxTwist_Grad << ",";
          }
          for (iPlane = 0; iPlane < nPlane * 6; iPlane++) {
            Gradient_file << Gradient[iPlane];
            if (iPlane != (nPlane * 6) - 1) Gradient_file << ",";
          }
        } else {
          if (geometry_container[ZONE_0]->GetnDim() == 3) {
            Gradient_file << Wing_Volume_Grad << "," << Wing_MinThickness_Grad << "," << Wing_MaxThickness_Grad << ","
                          << Wing_MinChord_Grad << "," << Wing_MaxChord_Grad << "," << Wing_MinLERadius_Grad << ","
                          << Wing_MaxLERadius_Grad << "," << Wing_MinToC_Grad << "," << Wing_MaxToC_Grad << ","
                          << Wing_ObjFun_MinToC_Grad << "," << Wing_MaxTwist_Grad << "," << Wing_MaxCurvature_Grad
                          << "," << Wing_MaxDihedral_Grad << ",";
          }
          for (iPlane = 0; iPlane < nPlane * 6; iPlane++) {
            Gradient_file << Gradient[iPlane];
            if (iPlane != (nPlane * 6) - 1) Gradient_file << ",";
          }
        }

        Gradient_file << endl;

        if (iDV != (config_container[ZONE_0]->GetnDV() - 1))
          cout << "-------------------------------------------------------------------------" << endl;
      }
    }

    if (rank == MASTER_NODE) Gradient_file.close();
  }

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Finalize Solver -------------------------" << endl;

  delete[] Xcoord_Airfoil;
  delete[] Ycoord_Airfoil;
  delete[] Zcoord_Airfoil;

  delete[] ObjectiveFunc;
  delete[] ObjectiveFunc_New;
  delete[] Gradient;

  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    delete Plane_P0[iPlane];
    delete Plane_Normal[iPlane];
  }
  delete[] Plane_P0;
  delete[] Plane_Normal;

  delete config;
  config = nullptr;

  if (rank == MASTER_NODE) cout << "Deleted main variables." << endl;

  if (geometry_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (geometry_container[iZone] != nullptr) {
        delete geometry_container[iZone];
      }
    }
    delete[] geometry_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;

  delete surface_movement;
  if (rank == MASTER_NODE) cout << "Deleted CSurfaceMovement class." << endl;

  if (FFDBox != nullptr) {
    for (iFFDBox = 0; iFFDBox < MAX_NUMBER_FFD; iFFDBox++) {
      if (FFDBox[iFFDBox] != nullptr) {
        delete FFDBox[iFFDBox];
      }
    }
    delete[] FFDBox;
  }
  if (rank == MASTER_NODE) cout << "Deleted CFreeFormDefBox class." << endl;

  if (config_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone] != nullptr) {
        delete config_container[iZone];
      }
    }
    delete[] config_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

  StopTime = SU2_MPI::Wtime();

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime - StartTime;
  if (rank == MASTER_NODE) {
    cout << "\n\nCompleted in " << fixed << UsedTime << " seconds on " << size;
    if (size == 1)
      cout << " core." << endl;
    else
      cout << " cores." << endl;
  }

  /*--- Exit the solver cleanly ---*/

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Exit Success (SU2_GEO) ------------------------" << endl << endl;

  /*--- Finalize MPI parallelization ---*/
  SU2_MPI::Finalize();

  return EXIT_SUCCESS;
}
