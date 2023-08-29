/*!
 * \file CSurfaceMovement.cpp
 * \brief Subroutines for moving mesh surface elements
 * \author F. Palacios, T. Economon, S. Padron
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

#include "../../include/grid_movement/CSurfaceMovement.hpp"
#include "../../include/toolboxes/C1DInterpolation.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"

CSurfaceMovement::CSurfaceMovement() : CGridMovement() {
  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  nFFDBox = 0;
  nLevel = 0;
  FFDBoxDefinition = false;
}

CSurfaceMovement::~CSurfaceMovement() = default;

vector<vector<su2double> > CSurfaceMovement::SetSurface_Deformation(CGeometry* geometry, CConfig* config) {
  unsigned short iFFDBox, iDV, iLevel, iChild, iParent, jFFDBox, iMarker;
  unsigned short Degree_Unitary[] = {1, 1, 1}, BSpline_Unitary[] = {2, 2, 2};
  su2double MaxDiff, Current_Scale, Ratio, New_Scale;
  string FFDBoxTag;
  bool allmoving;

  const bool cylindrical = (config->GetFFD_CoordSystem() == CYLINDRICAL);
  const bool spherical = (config->GetFFD_CoordSystem() == SPHERICAL);
  const bool polar = (config->GetFFD_CoordSystem() == POLAR);
  const bool cartesian = (config->GetFFD_CoordSystem() == CARTESIAN);
  const su2double BoundLimit = config->GetOpt_LineSearch_Bound();

  vector<vector<su2double> > totaldeformation;

  /*--- Setting the Free Form Deformation ---*/

  if (config->GetDesign_Variable(0) == FFD_SETTING) {
    /*--- Definition of the FFD deformation class ---*/

    FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];

    /*--- Read the FFD information from the config file ---*/

    ReadFFDInfo(geometry, config, FFDBox);

    /*--- If there is a FFDBox in the input file ---*/

    if (nFFDBox != 0) {
      /*--- if polar coordinates, trnasform the corner to polar ---*/

      if (cylindrical) {
        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
          FFDBox[iFFDBox]->SetCart2Cyl_CornerPoints(config);
        }
      } else if (spherical || polar) {
        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
          FFDBox[iFFDBox]->SetCart2Sphe_CornerPoints(config);
        }
      }

      /*--- If the FFDBox was not defined in the input file ---*/

      if ((rank == MASTER_NODE) && (GetnFFDBox() != 0)) {
        if (cartesian)
          cout << endl << "----------------- FFD technique (cartesian -> parametric) ---------------" << endl;
        else if (cylindrical)
          cout << endl << "----------------- FFD technique (cylinder -> parametric) ---------------" << endl;
        else if (spherical)
          cout << endl << "----------------- FFD technique (spherical -> parametric) ---------------" << endl;
        else if (polar)
          cout << endl << "----------------- FFD technique (polar -> parametric) ---------------" << endl;
      }

      /*--- Create a unitary FFDBox as baseline for other FFDBoxes shapes ---*/

      CFreeFormDefBox FFDBox_unitary(Degree_Unitary, BSpline_Unitary, BEZIER);
      FFDBox_unitary.SetUnitCornerPoints();

      /*--- Compute the control points of the unitary box, in this case the degree is 1 and the order is 2 ---*/

      FFDBox_unitary.SetControlPoints_Parallelepiped();

      for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
        /*--- Compute the support control points for the final FFD using the unitary box ---*/

        FFDBox_unitary.SetSupportCP(FFDBox[iFFDBox]);

        /*--- Compute control points in the support box ---*/

        FFDBox_unitary.SetSupportCPChange(FFDBox[iFFDBox]);

        /*--- Compute the parametric coordinates, it also find the points in
         the FFDBox using the parametrics coordinates ---*/

        SetParametricCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);

        /*--- If polar coordinates, transform the corners and control points to cartesians ---*/

        if (cylindrical) {
          FFDBox[iFFDBox]->SetCyl2Cart_CornerPoints(config);
          FFDBox[iFFDBox]->SetCyl2Cart_ControlPoints(config);
        } else if (spherical || polar) {
          FFDBox[iFFDBox]->SetSphe2Cart_CornerPoints(config);
          FFDBox[iFFDBox]->SetSphe2Cart_ControlPoints(config);
        }
      }
      /*--- Output original FFD FFDBox ---*/

      if (rank == MASTER_NODE) {
        for (unsigned short iFile = 0; iFile < config->GetnVolumeOutputFiles(); iFile++) {
          auto FileFormat = config->GetVolumeOutputFiles();
          if (isParaview(FileFormat[iFile])) {
            cout << "Writing a Paraview file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetParaview(geometry, iFFDBox, true);
            }
          } else if (isTecplot(FileFormat[iFile])) {
            cout << "Writing a Tecplot file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetTecplot(geometry, iFFDBox, true);
            }
          } else if (FileFormat[iFile] == OUTPUT_TYPE::CGNS) {
            cout << "Writing a CGNS file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetCGNS(geometry, iFFDBox, true);
            }
          }
        }
      }
    }

    else {
      SU2_MPI::Error("There are no FFD boxes in the mesh file!!", CURRENT_FUNCTION);
    }
  }

  /*--- Free Form deformation based ---*/

  if ((config->GetDesign_Variable(0) == FFD_CONTROL_POINT_2D) || (config->GetDesign_Variable(0) == FFD_CAMBER_2D) ||
      (config->GetDesign_Variable(0) == FFD_THICKNESS_2D) || (config->GetDesign_Variable(0) == FFD_CONTROL_POINT) ||
      (config->GetDesign_Variable(0) == FFD_NACELLE) || (config->GetDesign_Variable(0) == FFD_GULL) ||
      (config->GetDesign_Variable(0) == FFD_TWIST) || (config->GetDesign_Variable(0) == FFD_ROTATION) ||
      (config->GetDesign_Variable(0) == FFD_CONTROL_SURFACE) || (config->GetDesign_Variable(0) == FFD_CAMBER) ||
      (config->GetDesign_Variable(0) == FFD_THICKNESS) || (config->GetDesign_Variable(0) == FFD_ANGLE_OF_ATTACK)) {
    /*--- Definition of the FFD deformation class ---*/

    FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];

    /*--- Read the FFD information from the grid file ---*/

    ReadFFDInfo(geometry, config, FFDBox, config->GetMesh_FileName());

    /*--- If there is a FFDBox in the input file ---*/

    if (nFFDBox != 0) {
      /*--- If the FFDBox was not defined in the input file ---*/

      if (!GetFFDBoxDefinition()) {
        SU2_MPI::Error(
            string("There is not FFD box definition in the mesh file,\n") + string("run DV_KIND=FFD_SETTING first !!"),
            CURRENT_FUNCTION);
      }

      /* --- Check if the FFD boxes referenced in the design variable definition can be found --- */

      for (iDV = 0; iDV < config->GetnDV(); iDV++) {
        if (!CheckFFDBoxDefinition(config, iDV)) {
          SU2_MPI::Error(string("There is no FFD box with tag \"") + config->GetFFDTag(iDV) +
                             string("\" defined in the mesh file.\n") +
                             string("Check the definition of the design variables and/or the FFD settings !!"),
                         CURRENT_FUNCTION);
        }
      }

      /*--- Check that the user has specified a non-zero number of surfaces to move with DV_MARKER. ---*/

      if (config->GetnMarker_DV() == 0) {
        SU2_MPI::Error(string("No markers are specified in DV_MARKER, so no deformation will occur.\n") +
                           string("List markers to be deformed in DV_MARKER."),
                       CURRENT_FUNCTION);
      }

      /*--- Output original FFD FFDBox ---*/

      if ((rank == MASTER_NODE) && (config->GetKind_SU2() != SU2_COMPONENT::SU2_DOT)) {
        for (unsigned short iFile = 0; iFile < config->GetnVolumeOutputFiles(); iFile++) {
          auto FileFormat = config->GetVolumeOutputFiles();

          if (isParaview(FileFormat[iFile])) {
            cout << "Writing a Paraview file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetParaview(geometry, iFFDBox, true);
            }
          } else if (isTecplot(FileFormat[iFile])) {
            cout << "Writing a Tecplot file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetTecplot(geometry, iFFDBox, true);
            }
          } else if (FileFormat[iFile] == OUTPUT_TYPE::CGNS) {
            cout << "Writing a CGNS file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetCGNS(geometry, iFFDBox, true);
            }
          }
        }
      }

      /*--- If polar FFD, change the coordinates system ---*/

      if (cylindrical) {
        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
          FFDBox[iFFDBox]->SetCart2Cyl_CornerPoints(config);
          FFDBox[iFFDBox]->SetCart2Cyl_ControlPoints(config);
        }
      } else if (spherical || polar) {
        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
          FFDBox[iFFDBox]->SetCart2Sphe_CornerPoints(config);
          FFDBox[iFFDBox]->SetCart2Sphe_ControlPoints(config);
        }
      }

      /*--- Apply the deformation to the orifinal FFD box ---*/

      if ((rank == MASTER_NODE) && (GetnFFDBox() != 0))
        cout << endl << "----------------- FFD technique (parametric -> cartesian) ---------------" << endl;

      /*--- Loop over all the FFD boxes levels ---*/

      for (iLevel = 0; iLevel < GetnLevel(); iLevel++) {
        /*--- Loop over all FFD FFDBoxes ---*/

        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
          /*--- Check the level of the FFD box ---*/

          if (FFDBox[iFFDBox]->GetLevel() == iLevel) {
            /*--- Check the dimension of the FFD compared with the design variables ---*/

            if (rank == MASTER_NODE) cout << "Checking FFD box dimension." << endl;
            CheckFFDDimension(geometry, config, FFDBox[iFFDBox], iFFDBox);

            /*--- Compute intersections of the FFD box with the surface to eliminate design
             variables and satisfy surface continuity ---*/

            if (rank == MASTER_NODE) cout << "Checking FFD box intersections with the solid surfaces." << endl;
            CheckFFDIntersections(geometry, config, FFDBox[iFFDBox], iFFDBox);

            /*--- Compute the parametric coordinates of the child box
             control points (using the parent FFDBox)  ---*/

            for (iChild = 0; iChild < FFDBox[iFFDBox]->GetnChildFFDBox(); iChild++) {
              FFDBoxTag = FFDBox[iFFDBox]->GetChildFFDBoxTag(iChild);
              for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
                if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
              SetParametricCoordCP(geometry, config, FFDBox[iFFDBox], FFDBox[jFFDBox]);
            }

            /*--- Update the parametric coordinates if it is a child FFDBox ---*/

            if (iLevel > 0) UpdateParametricCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);

            /*--- Apply the design variables to the control point position ---*/
            ApplyDesignVariables(geometry, config, FFDBox, iFFDBox);

            /*--- Recompute cartesian coordinates using the new control point location ---*/

            MaxDiff = SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, false);

            if ((MaxDiff > BoundLimit) && (config->GetKind_SU2() == SU2_COMPONENT::SU2_DEF)) {
              if (rank == MASTER_NODE)
                cout << "Out-of-bounds, re-adjusting scale factor to safisfy line search limit." << endl;

              Current_Scale = config->GetOpt_RelaxFactor();
              Ratio = (BoundLimit / MaxDiff);
              New_Scale = Current_Scale * (Ratio - 1.0);
              config->SetOpt_RelaxFactor(New_Scale);

              /*--- Apply the design variables to the control point position ---*/
              ApplyDesignVariables(geometry, config, FFDBox, iFFDBox);

              /*--- Recompute cartesian coordinates using the new control point location ---*/

              MaxDiff = SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, false);
            }

            /*--- Set total deformation values in config ---*/
            if (config->GetKind_SU2() == SU2_COMPONENT::SU2_DEF) {
              totaldeformation.resize(config->GetnDV());
              for (iDV = 0; iDV < config->GetnDV(); iDV++) {
                totaldeformation[iDV].resize(config->GetnDV_Value(iDV));
                for (auto iDV_Value = 0u; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++) {
                  totaldeformation[iDV][iDV_Value] = config->GetDV_Value(iDV, iDV_Value);
                }
              }

              /*--- Get parameters for FFD self-intersection prevention ---*/
              bool FFD_IntPrev;
              unsigned short FFD_IntPrev_MaxIter, FFD_IntPrev_MaxDepth;

              tie(FFD_IntPrev, FFD_IntPrev_MaxIter, FFD_IntPrev_MaxDepth) = config->GetFFD_IntPrev();

              /*--- Calculate Jacobian determinant for FFD-deformation to check for self-intersection in FFD box.
              Procedure from J.E. Gain & N.A. Dodgson, "Preventing Self-Intersection under Free Form Deformation".
              IEEE transactions on Visualization and Computer Graphics, vol. 7 no. 4. October-December 2001 ---*/
              unsigned long nNegativeDeterminants = 0, nNegativeDeterminants_previous;
              su2double DeformationDifference = 1.0, DeformationFactor = 1.0;

              if (FFD_IntPrev) {
                nNegativeDeterminants = calculateJacobianDeterminant(geometry, config, FFDBox[iFFDBox]);
              }

              /*--- If enabled: start recursive procedure to decrease deformation magnitude
              to remove self-intersections in FFD box ---*/
              if (nNegativeDeterminants > 0) {
                if (rank == MASTER_NODE) {
                  cout << "Self-intersections within FFD box present. ";
                  cout << "Performing iterative deformation reduction procedure." << endl;
                }

                /*--- Lower the deformation magnitude and add this to the total deformation value in config ---*/
                for (iDV = 0; iDV < config->GetnDV(); iDV++) {
                  for (auto iDV_Value = 0u; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++) {
                    auto dv_value = config->GetDV_Value(iDV, iDV_Value);
                    config->SetDV_Value(iDV, iDV_Value, -dv_value / 2);
                    totaldeformation[iDV][iDV_Value] -= dv_value / 2;
                  }
                }

                DeformationDifference /= 2.0;
                DeformationFactor -= DeformationDifference;
                nNegativeDeterminants_previous = nNegativeDeterminants;

                /*--- Recursively check for self-intersections. ---*/
                unsigned short FFD_IntPrev_Iter, FFD_IntPrev_Depth = 0;

                for (FFD_IntPrev_Iter = 1; FFD_IntPrev_Iter <= FFD_IntPrev_MaxIter; FFD_IntPrev_Iter++) {
                  if (rank == MASTER_NODE) cout << "Checking FFD box intersections with the solid surfaces." << endl;
                  CheckFFDIntersections(geometry, config, FFDBox[iFFDBox], iFFDBox);

                  /*--- Compute the parametric coordinates of the child box
                  control points (using the parent FFDBox)  ---*/
                  for (iChild = 0; iChild < FFDBox[iFFDBox]->GetnChildFFDBox(); iChild++) {
                    FFDBoxTag = FFDBox[iFFDBox]->GetChildFFDBoxTag(iChild);
                    for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
                      if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
                    SetParametricCoordCP(geometry, config, FFDBox[iFFDBox], FFDBox[jFFDBox]);
                  }

                  /*--- Update the parametric coordinates if it is a child FFDBox ---*/
                  iLevel = FFDBox[iFFDBox]->GetLevel();
                  if (iLevel > 0) UpdateParametricCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);

                  /*--- Apply the design variables to the control point position ---*/
                  ApplyDesignVariables(geometry, config, FFDBox, iFFDBox);

                  /*--- Recompute cartesian coordinates using the new control point location ---*/
                  MaxDiff = SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, false);

                  /*--- Check for self-intersections in FFD box ---*/
                  nNegativeDeterminants = calculateJacobianDeterminant(geometry, config, FFDBox[iFFDBox]);

                  if (rank == MASTER_NODE) {
                    cout << "Amount of points with negative Jacobian determinant for iteration ";
                    cout << FFD_IntPrev_Iter << ": " << nNegativeDeterminants << endl;
                    cout << "Remaining amount of original deformation: " << DeformationFactor * 100.0 << " percent."
                         << endl;
                  }

                  /*--- Recursively change deformation magnitude.
                  Increase if there are no points with negative determinants, decrease otherwise. ---*/
                  if (nNegativeDeterminants == 0) {
                    DeformationDifference = abs(DeformationDifference / 2.0);

                    /*--- Update recursion depth if there are no points with negative determinant.
                    Quit if maximum depth is reached. ---*/
                    FFD_IntPrev_Depth++;

                    if (FFD_IntPrev_Depth == FFD_IntPrev_MaxDepth) {
                      if (rank == MASTER_NODE) {
                        cout << "Maximum recursion depth reached." << endl;
                        cout << "Remaining amount of original deformation: " << endl;
                        cout << DeformationFactor * 100.0 << " percent." << endl;
                      }
                      break;
                    }
                  } else {
                    DeformationDifference = -abs(DeformationDifference / 2.0);
                  }

                  if (FFD_IntPrev_Iter < FFD_IntPrev_MaxIter) {
                    DeformationFactor += DeformationDifference;
                  }

                  /*--- Set DV values for next iteration. Always decrease absolute value,
                  since starting point of iteration is previously deformed value. ---*/
                  if (FFD_IntPrev_Iter < FFD_IntPrev_MaxIter) {
                    su2double sign = -1.0;
                    if ((nNegativeDeterminants_previous > 0 && nNegativeDeterminants > 0) ||
                        (nNegativeDeterminants_previous == 0 && nNegativeDeterminants == 0)) {
                      sign = 1.0;
                    }
                    for (iDV = 0; iDV < config->GetnDV(); iDV++) {
                      for (auto iDV_Value = 0u; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++) {
                        auto dv_value = sign * config->GetDV_Value(iDV, iDV_Value);
                        config->SetDV_Value(iDV, iDV_Value, dv_value / 2.0);
                        totaldeformation[iDV][iDV_Value] += dv_value / 2.0;
                      }
                    }
                  }
                  nNegativeDeterminants_previous = nNegativeDeterminants;
                }
              }

            }  // end SU2_DEF

            /*--- Reparametrization of the parent FFD box ---*/

            for (iParent = 0; iParent < FFDBox[iFFDBox]->GetnParentFFDBox(); iParent++) {
              FFDBoxTag = FFDBox[iFFDBox]->GetParentFFDBoxTag(iParent);
              for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
                if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
              UpdateParametricCoord(geometry, config, FFDBox[jFFDBox], jFFDBox);
            }

            /*--- Compute the new location of the control points of the child boxes
             (using the parent FFDBox) ---*/

            for (iChild = 0; iChild < FFDBox[iFFDBox]->GetnChildFFDBox(); iChild++) {
              FFDBoxTag = FFDBox[iFFDBox]->GetChildFFDBoxTag(iChild);
              for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
                if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
              GetCartesianCoordCP(geometry, config, FFDBox[iFFDBox], FFDBox[jFFDBox]);
            }
          }
        }

        /*--- If polar, compute the cartesians coordinates ---*/

        if (cylindrical) {
          for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
            FFDBox[iFFDBox]->SetCyl2Cart_CornerPoints(config);
            FFDBox[iFFDBox]->SetCyl2Cart_ControlPoints(config);
          }
        } else if (spherical || polar) {
          for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
            FFDBox[iFFDBox]->SetSphe2Cart_CornerPoints(config);
            FFDBox[iFFDBox]->SetSphe2Cart_ControlPoints(config);
          }
        }

        /*--- Output the deformed FFD Boxes ---*/

        if ((rank == MASTER_NODE) && (config->GetKind_SU2() != SU2_COMPONENT::SU2_DOT)) {
          for (unsigned short iFile = 0; iFile < config->GetnVolumeOutputFiles(); iFile++) {
            auto FileFormat = config->GetVolumeOutputFiles();

            if (isParaview(FileFormat[iFile])) {
              cout << "Writing a Paraview file of the FFD boxes." << endl;
              for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
                FFDBox[iFFDBox]->SetParaview(geometry, iFFDBox, false);
              }
            } else if (isTecplot(FileFormat[iFile])) {
              cout << "Writing a Tecplot file of the FFD boxes." << endl;
              for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
                FFDBox[iFFDBox]->SetTecplot(geometry, iFFDBox, false);
              }
            } else if (FileFormat[iFile] == OUTPUT_TYPE::CGNS) {
              cout << "Writing a CGNS file of the FFD boxes." << endl;
              for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
                FFDBox[iFFDBox]->SetCGNS(geometry, iFFDBox, false);
              }
            }
          }
        }
      }
    }

    else {
      SU2_MPI::Error("There are no FFD Boxes in the mesh file!!", CURRENT_FUNCTION);
    }

  }

  /*--- External surface file based ---*/

  else if (config->GetDesign_Variable(0) == SURFACE_FILE) {
    /*--- Check whether a surface file exists for input ---*/
    ofstream Surface_File;
    string filename = config->GetDV_Filename();
    Surface_File.open(filename.c_str(), ios::in);

    /*--- A surface file does not exist, so write a new one for the
     markers that are specified as part of the motion. ---*/
    if (Surface_File.fail()) {
      if (rank == MASTER_NODE && size == SINGLE_NODE) {
        cout << "No surface positions file found. Writing a template file: " << filename << "." << endl;

        Surface_File.open(filename.c_str(), ios::out);
        Surface_File.precision(15);
        unsigned long iMarker, jPoint, GlobalIndex, iVertex;
        su2double* Coords;
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
          if (config->GetMarker_All_DV(iMarker) == YES) {
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
              jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              GlobalIndex = geometry->nodes->GetGlobalIndex(jPoint);
              Coords = geometry->nodes->GetCoord(jPoint);
              Surface_File << GlobalIndex << "\t" << Coords[0] << "\t" << Coords[1];
              if (geometry->GetnDim() == 2)
                Surface_File << endl;
              else
                Surface_File << "\t" << Coords[2] << endl;
            }
          }
        }
        Surface_File.close();

      } else {
        SU2_MPI::Error(
            "No surface positions file found and template writing not yet supported in parallel.\n To generate a "
            "template surface positions file, run SU2_DEF again in serial.",
            CURRENT_FUNCTION);
      }
    }

    else {
      /*--- A surface file exists, so read in the coordinates ---*/
      Surface_File.close();
      if (rank == MASTER_NODE) cout << "Updating the surface coordinates from the input file." << endl;
      SetExternal_Deformation(geometry, config, ZONE_0, 0);
    }

  }

  else if ((config->GetDesign_Variable(0) == ROTATION) || (config->GetDesign_Variable(0) == TRANSLATION) ||
           (config->GetDesign_Variable(0) == SCALE) || (config->GetDesign_Variable(0) == HICKS_HENNE) ||
           (config->GetDesign_Variable(0) == SURFACE_BUMP) || (config->GetDesign_Variable(0) == ANGLE_OF_ATTACK)) {
    /*--- Apply rotation, displacement and stretching design variables (this
     should be done before the bump function design variables) ---*/

    for (iDV = 0; iDV < config->GetnDV(); iDV++) {
      switch (config->GetDesign_Variable(iDV)) {
        case SCALE:
          SetScale(geometry, config, iDV, false);
          break;
        case TRANSLATION:
          SetTranslation(geometry, config, iDV, false);
          break;
        case ROTATION:
          SetRotation(geometry, config, iDV, false);
          break;
      }
    }

    /*--- Apply the design variables to the control point position ---*/

    for (iDV = 0; iDV < config->GetnDV(); iDV++) {
      switch (config->GetDesign_Variable(iDV)) {
        case HICKS_HENNE:
          SetHicksHenne(geometry, config, iDV, false);
          break;
      }
    }

    /*--- Apply the design variables to the control point position ---*/

    for (iDV = 0; iDV < config->GetnDV(); iDV++) {
      switch (config->GetDesign_Variable(iDV)) {
        case SURFACE_BUMP:
          SetSurface_Bump(geometry, config, iDV, false);
          break;
      }
    }

    /*--- Apply the angle of attack design variable ---*/

    for (iDV = 0; iDV < config->GetnDV(); iDV++) {
      switch (config->GetDesign_Variable(iDV)) {
        case ANGLE_OF_ATTACK:
          SetAngleOfAttack(geometry, config, iDV, false);
          break;
      }
    }

  }

  /*--- NACA_4Digits design variable ---*/

  else if (config->GetDesign_Variable(0) == NACA_4DIGITS) {
    SetNACA_4Digits(geometry, config);
  }

  /*--- Parabolic airfoil design variable ---*/

  else if (config->GetDesign_Variable(0) == PARABOLIC) {
    SetParabolic(geometry, config);
  }

  /*--- Airfoil from file design variable ---*/

  else if (config->GetDesign_Variable(0) == AIRFOIL) {
    SetAirfoil(geometry, config);
  }

  /*--- FFD setting ---*/

  else if (config->GetDesign_Variable(0) == FFD_SETTING) {
    if (rank == MASTER_NODE) cout << "No surface deformation (setting FFD)." << endl;
  }

  /*--- Scale, Translate, and Rotate will be done with rigid mesh transforms. ---*/

  else if ((config->GetDesign_Variable(0) == ROTATION) || (config->GetDesign_Variable(0) == TRANSLATION) ||
           (config->GetDesign_Variable(0) == SCALE)) {
    /*--- If all markers are deforming, use volume method.
     If only some are deforming, use surface method ---*/

    /*--- iDV was uninitialized, so hard-coding to one. Check intended
     behavior (might want to loop over all iDV in case we have trans & rotate. ---*/
    iDV = 0;
    allmoving = true;

    /*--- Loop over markers ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_DV(iMarker) == NO) allmoving = false;
    }

    if (!allmoving) {
      /*---Only some markers are moving, use the surface method ---*/
      if (config->GetDesign_Variable(0) == ROTATION) SetRotation(geometry, config, iDV, false);
      if (config->GetDesign_Variable(0) == SCALE) SetScale(geometry, config, iDV, false);
      if (config->GetDesign_Variable(0) == TRANSLATION) SetTranslation(geometry, config, iDV, false);
    } else {
      if (rank == MASTER_NODE) cout << "No surface deformation (scaling, rotation, or translation)." << endl;
    }
  }

  /*--- Design variable not implement ---*/

  else {
    if (rank == MASTER_NODE) cout << "Design Variable not implemented yet" << endl;
  }

  return totaldeformation;
}

void CSurfaceMovement::SetSurface_Derivative(CGeometry* geometry, CConfig* config) {
  su2double DV_Value = 0.0;

  unsigned short iDV = 0, iDV_Value = 0;

  for (iDV = 0; iDV < config->GetnDV(); iDV++) {
    for (iDV_Value = 0; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++) {
      DV_Value = config->GetDV_Value(iDV, iDV_Value);

      /*--- If value of the design variable is not 0.0 we apply the differentation.
       *     Note if multiple variables are non-zero, we end up with the sum of all the derivatives. ---*/

      if (DV_Value != 0.0) {
        DV_Value = 0.0;

        SU2_TYPE::SetDerivative(DV_Value, 1.0);

        config->SetDV_Value(iDV, iDV_Value, DV_Value);
      }
    }
  }

  /*--- Run the surface deformation with DV_Value = 0.0 (no deformation at all) ---*/

  SetSurface_Deformation(geometry, config);
}

void CSurfaceMovement::CopyBoundary(CGeometry* geometry, CConfig* config) {
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  su2double* Coord;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Coord = geometry->nodes->GetCoord(iPoint);
      geometry->vertex[iMarker][iVertex]->SetCoord(Coord);
    }
  }
}

void CSurfaceMovement::SetParametricCoord(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                          unsigned short iFFDBox) {
  const auto nDim = geometry->GetnDim();
  const bool cartesian = (config->GetFFD_CoordSystem() == CARTESIAN);
  const bool cylindrical = (config->GetFFD_CoordSystem() == CYLINDRICAL);
  const bool spherical = (config->GetFFD_CoordSystem() == SPHERICAL);
  const bool polar = (config->GetFFD_CoordSystem() == POLAR);

  /*--- Change order and control points reduce the complexity of the point inversion
   (this only works with boxes, in case of Bezier curves, and we maintain an internal copy). ---*/

  const bool BoxFFD = true;
  if (BoxFFD && (config->GetFFD_Blending() == BEZIER)) {
    for (int iOrder = 0; iOrder < 2; iOrder++) {
      for (int jOrder = 0; jOrder < 2; jOrder++) {
        for (int kOrder = 0; kOrder < 2; kOrder++) {
          unsigned short lOrder = 0;
          unsigned short mOrder = 0;
          unsigned short nOrder = 0;
          if (iOrder == 1) {
            lOrder = FFDBox->GetlOrder() - 1;
          }
          if (jOrder == 1) {
            mOrder = FFDBox->GetmOrder() - 1;
          }
          if (kOrder == 1) {
            nOrder = FFDBox->GetnOrder() - 1;
          }
          const auto* Coord = FFDBox->GetCoordControlPoints(lOrder, mOrder, nOrder);

          FFDBox->SetCoordControlPoints(Coord, iOrder, jOrder, kOrder);
        }
      }
    }

    FFDBox->SetlOrder(2);
    FFDBox->SetmOrder(2);
    FFDBox->SetnOrder(2);
    FFDBox->SetnControlPoints();
    FFDBox->BlendingFunction[0]->SetOrder(2, 2);
    FFDBox->BlendingFunction[1]->SetOrder(2, 2);
    FFDBox->BlendingFunction[2]->SetOrder(2, 2);
  }

  /*--- Point inversion algorithm with a basic box ---*/

  su2double my_MaxDiff = 0.0;
  unsigned long TotalVertex = 0;
  unsigned long VisitedVertex = 0;
  unsigned long MappedVertex = 0;
  su2double ParamCoordGuess[3] = {0.5, 0.5, 0.5};

  /*--- Check that the box is defined correctly for the preliminary point containment check,
   * by checking that the midpoint of the box is considered to be inside it. ---*/

  su2double BoxMidPoint[3] = {};
  for (int iOrder = 0; iOrder < 2; iOrder++) {
    for (int jOrder = 0; jOrder < 2; jOrder++) {
      for (int kOrder = 0; kOrder < 2; kOrder++) {
        const auto* Coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
        BoxMidPoint[0] += 0.125 * Coord[0];
        BoxMidPoint[1] += 0.125 * Coord[1];
        BoxMidPoint[2] += 0.125 * Coord[2];
      }
    }
  }
  if (!FFDBox->CheckPointInsideFFD(BoxMidPoint)) {
    SU2_MPI::Error("The FFD box '" + FFDBox->GetTag() +
                       "' is not properly defined. The first 4 points must be listed counter\n"
                       "clockwise, such that applying the right-hand rule results in a vector into the box.\n"
                       "This is according to the VTK hexahedron ordering:\n"
                       "    7 +----+ 6 \n"
                       "     /|   /|   \n"
                       "  4 +----+5|   \n"
                       "    |3+--|-+ 2 \n"
                       "    |/   |/    \n"
                       "  0 +----+ 1   \n"
                       "The CCW convention also applies in 2D, where only the bottom face is specified.",
                   CURRENT_FUNCTION);
  }

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      TotalVertex += geometry->nVertex[iMarker];

      for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        /*--- Get the cartesian coordinates ---*/

        su2double CartCoord[3] = {};
        for (auto iDim = 0u; iDim < nDim; iDim++) CartCoord[iDim] = geometry->vertex[iMarker][iVertex]->GetCoord(iDim);

        /*--- Transform the cartesian into polar ---*/

        if (!cartesian) {
          const su2double X_0 = config->GetFFD_Axis(0);
          const su2double Y_0 = config->GetFFD_Axis(1);
          const su2double Z_0 = config->GetFFD_Axis(2);

          const su2double Xbar = CartCoord[0] - X_0;
          const su2double Ybar = CartCoord[1] - Y_0;
          const su2double Zbar = CartCoord[2] - Z_0;

          CartCoord[1] = atan2(Zbar, Ybar);
          if (CartCoord[1] > PI_NUMBER / 2.0) CartCoord[1] -= 2.0 * PI_NUMBER;

          if (cylindrical) {
            CartCoord[0] = sqrt(Ybar * Ybar + Zbar * Zbar);
            CartCoord[2] = Xbar;
          } else if (spherical || polar) {
            CartCoord[0] = sqrt(Xbar * Xbar + Ybar * Ybar + Zbar * Zbar);
            CartCoord[2] = acos(Xbar / CartCoord[0]);
          }
        }

        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- If the point is inside the FFD, compute the value of the parametric coordinate. ---*/

        if (FFDBox->CheckPointInsideFFD(CartCoord)) {
          /*--- Find the parametric coordinate ---*/

          ++VisitedVertex;
          auto* ParamCoord = FFDBox->GetParametricCoord_Iterative(iPoint, CartCoord, ParamCoordGuess, config);

          /*--- Compute the cartesian coordinates using the parametric coordinates
           to check that everything is correct ---*/

          const auto* CartCoordNew = FFDBox->EvalCartesianCoord(ParamCoord);

          /*--- Compute max difference between original value and the recomputed value ---*/

          const su2double Diff = GeometryToolbox::Distance(nDim, CartCoordNew, CartCoord);
          my_MaxDiff = max(my_MaxDiff, Diff);

          /*--- If the parametric coordinates are in (-tol, 1+tol) the point belongs to the FFDBox ---*/

          if (((ParamCoord[0] >= -config->GetFFD_Tol()) && (ParamCoord[0] <= 1.0 + config->GetFFD_Tol())) &&
              ((ParamCoord[1] >= -config->GetFFD_Tol()) && (ParamCoord[1] <= 1.0 + config->GetFFD_Tol())) &&
              ((ParamCoord[2] >= -config->GetFFD_Tol()) && (ParamCoord[2] <= 1.0 + config->GetFFD_Tol()))) {
            /*--- Rectification of the initial tolerance (we have detected situations
             where 0.0 and 1.0 do not work properly. ---*/

            const su2double lower_limit = config->GetFFD_Tol();
            const su2double upper_limit = 1.0 - config->GetFFD_Tol();

            ParamCoord[0] = fmin(fmax(lower_limit, ParamCoord[0]), upper_limit);
            ParamCoord[1] = fmin(fmax(lower_limit, ParamCoord[1]), upper_limit);
            ParamCoord[2] = fmin(fmax(lower_limit, ParamCoord[2]), upper_limit);

            /*--- Set the value of the parametric coordinate ---*/

            ++MappedVertex;
            FFDBox->Set_MarkerIndex(iMarker);
            FFDBox->Set_VertexIndex(iVertex);
            FFDBox->Set_PointIndex(iPoint);
            FFDBox->Set_ParametricCoord(ParamCoord);
            FFDBox->Set_CartesianCoord(CartCoord);

            ParamCoordGuess[0] = ParamCoord[0];
            ParamCoordGuess[1] = ParamCoord[1];
            ParamCoordGuess[2] = ParamCoord[2];
          }

          if (Diff >= config->GetFFD_Tol()) {
            cout << "Please check this point: Local (" << ParamCoord[0] << " " << ParamCoord[1] << " " << ParamCoord[2]
                 << ") <-> Global (" << CartCoord[0] << " " << CartCoord[1] << " " << CartCoord[2] << ") <-> Error "
                 << Diff << " vs " << config->GetFFD_Tol() << "." << endl;
          }
        }
      }
    }
  }

  su2double MaxDiff = 0.0;
  SU2_MPI::Reduce(&my_MaxDiff, &MaxDiff, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, SU2_MPI::GetComm());

  unsigned long GlobalVertex = 0, GlobalVisited = 0, GlobalMapped = 0;
  SU2_MPI::Reduce(&TotalVertex, &GlobalVertex, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
  SU2_MPI::Reduce(&VisitedVertex, &GlobalVisited, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
  SU2_MPI::Reduce(&MappedVertex, &GlobalMapped, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());

  if (rank == MASTER_NODE) {
    cout << "Computed parametric coord for FFD box '" << FFDBox->GetTag() << "'\n";
    cout << "  Number of vertices (Total, Inside FFD, Mapped to FFD): " << GlobalVertex;
    cout << ", " << GlobalVisited << ", " << GlobalMapped << "\n";
    cout << "  Max coord difference: " << MaxDiff << "\n";
  }

  /*--- After the point inversion, copy the original information back. ---*/

  if (BoxFFD) {
    FFDBox->SetOriginalControlPoints();
    if (config->GetFFD_Blending() == BEZIER) {
      FFDBox->BlendingFunction[0]->SetOrder(FFDBox->GetlOrder(), FFDBox->GetlOrder());
      FFDBox->BlendingFunction[1]->SetOrder(FFDBox->GetmOrder(), FFDBox->GetmOrder());
      FFDBox->BlendingFunction[2]->SetOrder(FFDBox->GetnOrder(), FFDBox->GetnOrder());
    }
  }
}

void CSurfaceMovement::SetParametricCoordCP(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBoxParent,
                                            CFreeFormDefBox* FFDBoxChild) {
  unsigned short iOrder, jOrder, kOrder;
  su2double *CartCoord, *ParamCoord, ParamCoordGuess[3];

  for (iOrder = 0; iOrder < FFDBoxChild->GetlOrder(); iOrder++)
    for (jOrder = 0; jOrder < FFDBoxChild->GetmOrder(); jOrder++)
      for (kOrder = 0; kOrder < FFDBoxChild->GetnOrder(); kOrder++) {
        CartCoord = FFDBoxChild->GetCoordControlPoints(iOrder, jOrder, kOrder);
        ParamCoord = FFDBoxParent->GetParametricCoord_Iterative(0, CartCoord, ParamCoordGuess, config);
        FFDBoxChild->SetParCoordControlPoints(ParamCoord, iOrder, jOrder, kOrder);
      }

  if (rank == MASTER_NODE)
    cout << "Compute parametric coord (CP) | FFD parent box: " << FFDBoxParent->GetTag()
         << ". FFD child box: " << FFDBoxChild->GetTag() << "." << endl;
}

void CSurfaceMovement::GetCartesianCoordCP(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBoxParent,
                                           CFreeFormDefBox* FFDBoxChild) {
  unsigned short iOrder, jOrder, kOrder, iDim;
  su2double *CartCoord, *ParamCoord;

  for (iOrder = 0; iOrder < FFDBoxChild->GetlOrder(); iOrder++)
    for (jOrder = 0; jOrder < FFDBoxChild->GetmOrder(); jOrder++)
      for (kOrder = 0; kOrder < FFDBoxChild->GetnOrder(); kOrder++) {
        ParamCoord = FFDBoxChild->GetParCoordControlPoints(iOrder, jOrder, kOrder);

        /*--- Clip the value of the parametric coordinates (just in case)  ---*/
        for (iDim = 0; iDim < 3; iDim++) {
          if (ParamCoord[iDim] >= 1.0) ParamCoord[iDim] = 1.0;
          if (ParamCoord[iDim] <= 0.0) ParamCoord[iDim] = 0.0;
        }

        CartCoord = FFDBoxParent->EvalCartesianCoord(ParamCoord);
        FFDBoxChild->SetCoordControlPoints(CartCoord, iOrder, jOrder, kOrder);
        FFDBoxChild->SetCoordControlPoints_Copy(CartCoord, iOrder, jOrder, kOrder);
      }

  if (rank == MASTER_NODE)
    cout << "Update cartesian coord (CP)   | FFD parent box: " << FFDBoxParent->GetTag()
         << ". FFD child box: " << FFDBoxChild->GetTag() << "." << endl;
}

void CSurfaceMovement::CheckFFDDimension(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                         unsigned short iFFDBox) {
  unsigned short iIndex, jIndex, kIndex, lDegree, mDegree, nDegree, iDV;
  bool OutOffLimits;
  bool polar = (config->GetFFD_CoordSystem() == POLAR);

  lDegree = FFDBox->GetlOrder() - 1;
  mDegree = FFDBox->GetmOrder() - 1;
  nDegree = FFDBox->GetnOrder() - 1;

  OutOffLimits = false;
  for (iDV = 0; iDV < config->GetnDV(); iDV++) {
    if (config->GetFFDTag(iDV) == FFDBox->GetTag()) {
      switch (config->GetDesign_Variable(iDV)) {
        case FFD_CONTROL_POINT_2D:
          if (polar) {
            iIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
            kIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 2)));
            if ((iIndex > lDegree) || (kIndex > nDegree)) OutOffLimits = true;
          } else {
            iIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
            jIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 2)));
            if ((iIndex > lDegree) || (jIndex > mDegree)) OutOffLimits = true;
          }
          break;
        case FFD_CAMBER:
        case FFD_THICKNESS:
          iIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
          jIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 2)));
          if ((iIndex > lDegree) || (jIndex > mDegree)) OutOffLimits = true;
          break;
        case FFD_CAMBER_2D:
        case FFD_THICKNESS_2D:
          iIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
          if (iIndex > lDegree) OutOffLimits = true;
          break;
        case FFD_CONTROL_POINT:
        case FFD_NACELLE:
          iIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
          jIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 2)));
          kIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 3)));
          if ((iIndex > lDegree) || (jIndex > mDegree) || (kIndex > nDegree)) OutOffLimits = true;
          break;
        case FFD_GULL:
        case FFD_TWIST:
          jIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
          if (jIndex > mDegree) OutOffLimits = true;
          break;
      }
    }
  }

  if (rank == MASTER_NODE) {
    if (OutOffLimits) {
      char buf1[100], buf2[100];
      SPRINTF(buf1, "Design variables out off FFD limits (%u, %u, %u).\n", lDegree, mDegree, nDegree);
      SPRINTF(buf2, "Please check the ijk indices of the design variables.");
      SU2_MPI::Error(string(buf1) + string(buf2), CURRENT_FUNCTION);
    }
  }

  /*--- This barrier is important to guaranty that we will stop the software in a clean way ---*/

#ifdef HAVE_MPI
  SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
}

void CSurfaceMovement::CheckFFDIntersections(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                             unsigned short iFFDBox) {
  su2double Coord_0[] = {0, 0, 0}, Coord_1[] = {0, 0, 0};
  unsigned short index, iMarker, iNode, jNode, lDegree, mDegree, nDegree, iDim;
  unsigned long iElem, iPoint, jPoint;
  bool IPlane_Intersect_A = false, IPlane_Intersect_B = false;
  bool JPlane_Intersect_A = false, JPlane_Intersect_B = false;
  bool KPlane_Intersect_A = false, KPlane_Intersect_B = false;
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;

  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();
  bool FFD_Symmetry_Plane = config->GetFFD_Symmetry_Plane();
  bool cylindrical = (config->GetFFD_CoordSystem() == CYLINDRICAL);
  bool spherical = (config->GetFFD_CoordSystem() == SPHERICAL);
  bool polar = (config->GetFFD_CoordSystem() == POLAR);
  bool cartesian = (config->GetFFD_CoordSystem() == CARTESIAN);

  lDegree = FFDBox->GetlOrder() - 1;
  mDegree = FFDBox->GetmOrder() - 1;
  nDegree = FFDBox->GetnOrder() - 1;

  if (config->GetFFD_Continuity() != USER_INPUT) {
    /*--- Check intersection with plane i=0 ---*/

    su2double* IPlane_Coord_0_A = FFDBox->GetCoordControlPoints(0, 0, 0);
    su2double* IPlane_Coord_1_A = FFDBox->GetCoordControlPoints(0, 0, nDegree);
    su2double* IPlane_Coord_2_A = FFDBox->GetCoordControlPoints(0, mDegree, 0);

    su2double* IPlane_Coord_0_A_ = FFDBox->GetCoordControlPoints(0, mDegree, nDegree);
    su2double* IPlane_Coord_1_A_ = FFDBox->GetCoordControlPoints(0, mDegree, 0);
    su2double* IPlane_Coord_2_A_ = FFDBox->GetCoordControlPoints(0, 0, nDegree);

    /*--- Check intersection with plane i=lDegree ---*/

    su2double* IPlane_Coord_0_B = FFDBox->GetCoordControlPoints(lDegree, 0, 0);
    su2double* IPlane_Coord_1_B = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);
    su2double* IPlane_Coord_2_B = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);

    su2double* IPlane_Coord_0_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, nDegree);
    su2double* IPlane_Coord_1_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);
    su2double* IPlane_Coord_2_B_ = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);

    /*--- Check intersection with plane j=0 ---*/

    su2double* JPlane_Coord_0_A = FFDBox->GetCoordControlPoints(0, 0, 0);
    su2double* JPlane_Coord_1_A = FFDBox->GetCoordControlPoints(0, 0, nDegree);
    su2double* JPlane_Coord_2_A = FFDBox->GetCoordControlPoints(lDegree, 0, 0);

    su2double* JPlane_Coord_0_A_ = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);
    su2double* JPlane_Coord_1_A_ = FFDBox->GetCoordControlPoints(lDegree, 0, 0);
    su2double* JPlane_Coord_2_A_ = FFDBox->GetCoordControlPoints(0, 0, nDegree);

    /*--- Check intersection with plane j=mDegree ---*/

    su2double* JPlane_Coord_0_B = FFDBox->GetCoordControlPoints(0, mDegree, 0);
    su2double* JPlane_Coord_1_B = FFDBox->GetCoordControlPoints(0, mDegree, nDegree);
    su2double* JPlane_Coord_2_B = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);

    su2double* JPlane_Coord_0_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, nDegree);
    su2double* JPlane_Coord_1_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);
    su2double* JPlane_Coord_2_B_ = FFDBox->GetCoordControlPoints(0, mDegree, nDegree);

    /*--- Check intersection with plane k=0 ---*/

    su2double* KPlane_Coord_0_A = FFDBox->GetCoordControlPoints(0, 0, 0);
    su2double* KPlane_Coord_1_A = FFDBox->GetCoordControlPoints(0, mDegree, 0);
    su2double* KPlane_Coord_2_A = FFDBox->GetCoordControlPoints(lDegree, 0, 0);

    su2double* KPlane_Coord_0_A_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);
    su2double* KPlane_Coord_1_A_ = FFDBox->GetCoordControlPoints(lDegree, 0, 0);
    su2double* KPlane_Coord_2_A_ = FFDBox->GetCoordControlPoints(0, mDegree, 0);

    /*--- Check intersection with plane k=nDegree ---*/

    su2double* KPlane_Coord_0_B = FFDBox->GetCoordControlPoints(0, 0, nDegree);
    su2double* KPlane_Coord_1_B = FFDBox->GetCoordControlPoints(0, mDegree, nDegree);
    su2double* KPlane_Coord_2_B = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);

    su2double* KPlane_Coord_0_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, nDegree);
    su2double* KPlane_Coord_1_B_ = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);
    su2double* KPlane_Coord_2_B_ = FFDBox->GetCoordControlPoints(0, mDegree, nDegree);

    /*--- Loop over all the grid triangles ---*/

    IPlane_Intersect_A = false;
    IPlane_Intersect_B = false;
    JPlane_Intersect_A = false;
    JPlane_Intersect_B = false;
    KPlane_Intersect_A = false;
    KPlane_Intersect_B = false;

    /*--- Only the markers in the moving list ---*/

    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DEF)) ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_GEO)) ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DOT)) ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (config->GetDirectDiff() == D_DESIGN))) {
        for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
          for (iNode = 0; iNode < geometry->bound[iMarker][iElem]->GetnNodes(); iNode++) {
            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);

            for (jNode = 0; jNode < geometry->bound[iMarker][iElem]->GetnNodes(); jNode++) {
              jPoint = geometry->bound[iMarker][iElem]->GetNode(jNode);

              if (jPoint > iPoint) {
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  Coord_0[iDim] = geometry->nodes->GetCoord(iPoint, iDim);
                  Coord_1[iDim] = geometry->nodes->GetCoord(jPoint, iDim);
                }

                /*--- Write the coordinates in the right parametric system ---*/

                if (cylindrical) {
                  X_0 = config->GetFFD_Axis(0);
                  Y_0 = config->GetFFD_Axis(1);
                  Z_0 = config->GetFFD_Axis(2);

                  Xbar = Coord_0[0] - X_0;
                  Ybar = Coord_0[1] - Y_0;
                  Zbar = Coord_0[2] - Z_0;

                  Coord_0[0] = sqrt(Ybar * Ybar + Zbar * Zbar);
                  Coord_0[1] = atan2(Zbar, Ybar);
                  if (Coord_0[1] > PI_NUMBER / 2.0) Coord_0[1] -= 2.0 * PI_NUMBER;
                  Coord_0[2] = Xbar;

                  Xbar = Coord_1[0] - X_0;
                  Ybar = Coord_1[1] - Y_0;
                  Zbar = Coord_1[2] - Z_0;

                  Coord_1[0] = sqrt(Ybar * Ybar + Zbar * Zbar);
                  Coord_1[1] = atan2(Zbar, Ybar);
                  if (Coord_1[1] > PI_NUMBER / 2.0) Coord_1[1] -= 2.0 * PI_NUMBER;
                  Coord_1[2] = Xbar;

                }

                else if (spherical || polar) {
                  X_0 = config->GetFFD_Axis(0);
                  Y_0 = config->GetFFD_Axis(1);
                  Z_0 = config->GetFFD_Axis(2);

                  Xbar = Coord_0[0] - X_0;
                  Ybar = Coord_0[1] - Y_0;
                  Zbar = Coord_0[2] - Z_0;

                  Coord_0[0] = sqrt(Xbar * Xbar + Ybar * Ybar + Zbar * Zbar);
                  Coord_0[1] = atan2(Zbar, Ybar);
                  if (Coord_0[1] > PI_NUMBER / 2.0) Coord_0[1] -= 2.0 * PI_NUMBER;
                  Coord_0[2] = acos(Xbar / Coord_0[0]);

                  Xbar = Coord_1[0] - X_0;
                  Ybar = Coord_1[1] - Y_0;
                  Zbar = Coord_1[2] - Z_0;

                  Coord_1[0] = sqrt(Xbar * Xbar + Ybar * Ybar + Zbar * Zbar);
                  Coord_1[1] = atan2(Zbar, Ybar);
                  if (Coord_1[1] > PI_NUMBER / 2.0) Coord_1[1] -= 2.0 * PI_NUMBER;
                  Coord_1[2] = acos(Xbar / Coord_1[0]);
                }

                if (geometry->GetnDim() == 3) {
                  if (!IPlane_Intersect_A) {
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_A, IPlane_Coord_1_A,
                                                            IPlane_Coord_2_A)) {
                      IPlane_Intersect_A = true;
                    }
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_A_, IPlane_Coord_1_A_,
                                                            IPlane_Coord_2_A_)) {
                      IPlane_Intersect_A = true;
                    }
                  }

                  if (!IPlane_Intersect_B) {
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_B, IPlane_Coord_1_B,
                                                            IPlane_Coord_2_B)) {
                      IPlane_Intersect_B = true;
                    }
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_B_, IPlane_Coord_1_B_,
                                                            IPlane_Coord_2_B_)) {
                      IPlane_Intersect_B = true;
                    }
                  }

                  if ((!JPlane_Intersect_A) && (!FFD_Symmetry_Plane)) {
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_A, JPlane_Coord_1_A,
                                                            JPlane_Coord_2_A)) {
                      JPlane_Intersect_A = true;
                    }
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_A_, JPlane_Coord_1_A_,
                                                            JPlane_Coord_2_A_)) {
                      JPlane_Intersect_A = true;
                    }
                  }

                  if (cartesian) {
                    if ((!JPlane_Intersect_B) && (!FFD_Symmetry_Plane)) {
                      if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_B, JPlane_Coord_1_B,
                                                              JPlane_Coord_2_B)) {
                        JPlane_Intersect_B = true;
                      }
                      if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_B_, JPlane_Coord_1_B_,
                                                              JPlane_Coord_2_B_)) {
                        JPlane_Intersect_B = true;
                      }
                    }
                  } else {
                    if (!JPlane_Intersect_B) {
                      if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_B, JPlane_Coord_1_B,
                                                              JPlane_Coord_2_B)) {
                        JPlane_Intersect_B = true;
                      }
                      if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_B_, JPlane_Coord_1_B_,
                                                              JPlane_Coord_2_B_)) {
                        JPlane_Intersect_B = true;
                      }
                    }
                  }

                  if (!KPlane_Intersect_A) {
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_A, KPlane_Coord_1_A,
                                                            KPlane_Coord_2_A)) {
                      KPlane_Intersect_A = true;
                    }
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_A_, KPlane_Coord_1_A_,
                                                            KPlane_Coord_2_A_)) {
                      KPlane_Intersect_A = true;
                    }
                  }

                  if (!KPlane_Intersect_B) {
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_B, KPlane_Coord_1_B,
                                                            KPlane_Coord_2_B)) {
                      KPlane_Intersect_B = true;
                    }
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_B_, KPlane_Coord_1_B_,
                                                            KPlane_Coord_2_B_)) {
                      KPlane_Intersect_B = true;
                    }
                  }

                } else {
                  if (!IPlane_Intersect_A) {
                    if (geometry->SegmentIntersectsLine(Coord_0, Coord_1, IPlane_Coord_0_A, IPlane_Coord_2_A)) {
                      IPlane_Intersect_A = true;
                    }
                  }
                  if (!IPlane_Intersect_B) {
                    if (geometry->SegmentIntersectsLine(Coord_0, Coord_1, IPlane_Coord_0_B, IPlane_Coord_2_B)) {
                      IPlane_Intersect_B = true;
                    }
                  }
                  if (!JPlane_Intersect_A) {
                    if (geometry->SegmentIntersectsLine(Coord_0, Coord_1, JPlane_Coord_0_A, JPlane_Coord_2_A)) {
                      JPlane_Intersect_A = true;
                    }
                  }
                  if (!JPlane_Intersect_B) {
                    if (geometry->SegmentIntersectsLine(Coord_0, Coord_1, JPlane_Coord_0_B, JPlane_Coord_2_B)) {
                      JPlane_Intersect_B = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    /*--- Comunicate the planes that interesect the surface ---*/

    unsigned short MyCode[6] = {0, 0, 0, 0, 0, 0}, Code[6] = {0, 0, 0, 0, 0, 0};

    if (IPlane_Intersect_A) MyCode[0] = 1;
    if (IPlane_Intersect_B) MyCode[1] = 1;
    if (JPlane_Intersect_A) MyCode[2] = 1;
    if (JPlane_Intersect_B) MyCode[3] = 1;
    if (KPlane_Intersect_A) MyCode[4] = 1;
    if (KPlane_Intersect_B) MyCode[5] = 1;

#ifdef HAVE_MPI

    /*--- Add SU2_MPI::Allreduce information using all the nodes ---*/

    SU2_MPI::Allreduce(&MyCode, &Code, 6, MPI_UNSIGNED_SHORT, MPI_SUM, SU2_MPI::GetComm());

#else

    Code[0] = MyCode[0];
    Code[1] = MyCode[1];
    Code[2] = MyCode[2];
    Code[3] = MyCode[3];
    Code[4] = MyCode[4];
    Code[5] = MyCode[5];

#endif

    IPlane_Intersect_A = Code[0] != 0;
    IPlane_Intersect_B = Code[1] != 0;
    JPlane_Intersect_A = Code[2] != 0;
    JPlane_Intersect_B = Code[3] != 0;
    KPlane_Intersect_A = Code[4] != 0;
    KPlane_Intersect_B = Code[5] != 0;

    /*--- Screen output ---*/

    if (rank == MASTER_NODE) {
      if (IPlane_Intersect_A || IPlane_Intersect_B || JPlane_Intersect_A || JPlane_Intersect_B || KPlane_Intersect_A ||
          KPlane_Intersect_B) {
        cout << "The FFD planes ";

        if (cartesian) {
          if (IPlane_Intersect_A) cout << "i=0, ";
          if (IPlane_Intersect_B) cout << "i=" << lDegree << ", ";
          if (JPlane_Intersect_A) cout << "j=0, ";
          if (JPlane_Intersect_B) cout << "j=" << mDegree << ", ";
          if (KPlane_Intersect_A) cout << "k=0, ";
          if (KPlane_Intersect_B) cout << "k=" << nDegree << ", ";
        } else if (cylindrical) {
          if (IPlane_Intersect_A) cout << "r=0, ";
          if (IPlane_Intersect_B) cout << "r=" << lDegree << ", ";
          if (JPlane_Intersect_A) cout << "theta=0, ";
          if (JPlane_Intersect_B) cout << "theta=" << mDegree << ", ";
          if (KPlane_Intersect_A) cout << "z=0, ";
          if (KPlane_Intersect_B) cout << "z=" << nDegree << ", ";
        } else if (spherical) {
          if (IPlane_Intersect_A) cout << "r=0, ";
          if (IPlane_Intersect_B) cout << "r=" << lDegree << ", ";
          if (JPlane_Intersect_A) cout << "theta=0, ";
          if (JPlane_Intersect_B) cout << "theta=" << mDegree << ", ";
          if (KPlane_Intersect_A) cout << "phi=0, ";
          if (KPlane_Intersect_B) cout << "phi=" << nDegree << ", ";
        } else if (polar) {
          if (IPlane_Intersect_A) cout << "r=0, ";
          if (IPlane_Intersect_B) cout << "r=" << lDegree << ", ";
          if (KPlane_Intersect_A) cout << "theta=0, ";
          if (KPlane_Intersect_B) cout << "theta=" << nDegree << ", ";
        }

        cout << "intersect solid surfaces." << endl;
      }
    }
  }

  /*--- Fix the FFD planes based on the intersections with solid surfaces,
   and the continuity level, check that we have enough degree for the continuity
   that we are looking for ---*/

  if (config->GetFFD_Continuity() == USER_INPUT) {
    if (rank == MASTER_NODE) cout << "SU2 is fixing user's input planes." << endl;

    for (index = 0; index < config->GetnFFD_Fix_IDir(); index++)
      if ((config->GetFFD_Fix_IDir(index) <= lDegree) && (config->GetFFD_Fix_IDir(index) >= 0))
        FFDBox->Set_Fix_IPlane(config->GetFFD_Fix_IDir(index));
    for (index = 0; index < config->GetnFFD_Fix_JDir(); index++)
      if ((config->GetFFD_Fix_JDir(index) <= mDegree) && (config->GetFFD_Fix_JDir(index) >= 0))
        FFDBox->Set_Fix_JPlane(config->GetFFD_Fix_JDir(index));
    for (index = 0; index < config->GetnFFD_Fix_KDir(); index++)
      if ((config->GetFFD_Fix_KDir(index) <= nDegree) && (config->GetFFD_Fix_KDir(index) >= 0))
        FFDBox->Set_Fix_KPlane(config->GetFFD_Fix_KDir(index));
  }

  if (config->GetFFD_Continuity() == DERIVATIVE_NONE) {
    if (rank == MASTER_NODE) cout << "SU2 is fixing the planes to maintain a continuous surface." << endl;

    if (IPlane_Intersect_A) {
      FFDBox->Set_Fix_IPlane(0);
    }
    if (IPlane_Intersect_B) {
      FFDBox->Set_Fix_IPlane(lDegree);
    }
    if (JPlane_Intersect_A) {
      FFDBox->Set_Fix_JPlane(0);
    }
    if (JPlane_Intersect_B) {
      FFDBox->Set_Fix_JPlane(mDegree);
    }
    if (KPlane_Intersect_A) {
      FFDBox->Set_Fix_KPlane(0);
    }
    if (KPlane_Intersect_B) {
      FFDBox->Set_Fix_KPlane(nDegree);
    }
  }

  if (config->GetFFD_Continuity() == DERIVATIVE_1ST) {
    if (rank == MASTER_NODE) cout << "SU2 is fixing the planes to maintain a continuous 1st order derivative." << endl;

    if (IPlane_Intersect_A) {
      FFDBox->Set_Fix_IPlane(0);
      FFDBox->Set_Fix_IPlane(1);
    }
    if (IPlane_Intersect_B) {
      FFDBox->Set_Fix_IPlane(lDegree);
      FFDBox->Set_Fix_IPlane(lDegree - 1);
    }
    if (JPlane_Intersect_A) {
      FFDBox->Set_Fix_JPlane(0);
      FFDBox->Set_Fix_JPlane(1);
    }
    if (JPlane_Intersect_B) {
      FFDBox->Set_Fix_JPlane(mDegree);
      FFDBox->Set_Fix_JPlane(mDegree - 1);
    }
    if (KPlane_Intersect_A) {
      FFDBox->Set_Fix_KPlane(0);
      FFDBox->Set_Fix_KPlane(1);
    }
    if (KPlane_Intersect_B) {
      FFDBox->Set_Fix_KPlane(nDegree);
      FFDBox->Set_Fix_KPlane(nDegree - 1);
    }
  }

  if (config->GetFFD_Continuity() == DERIVATIVE_2ND) {
    if (rank == MASTER_NODE) cout << "SU2 is fixing the planes to maintain a continuous 2nd order derivative." << endl;

    if ((IPlane_Intersect_A) && (lDegree > 1)) {
      FFDBox->Set_Fix_IPlane(0);
      FFDBox->Set_Fix_IPlane(1);
      FFDBox->Set_Fix_IPlane(2);
    }
    if ((IPlane_Intersect_B) && (lDegree > 1)) {
      FFDBox->Set_Fix_IPlane(lDegree);
      FFDBox->Set_Fix_IPlane(lDegree - 1);
      FFDBox->Set_Fix_IPlane(lDegree - 2);
    }
    if ((JPlane_Intersect_A) && (mDegree > 1)) {
      FFDBox->Set_Fix_JPlane(0);
      FFDBox->Set_Fix_JPlane(1);
      FFDBox->Set_Fix_JPlane(2);
    }
    if ((JPlane_Intersect_B) && (mDegree > 1)) {
      FFDBox->Set_Fix_JPlane(mDegree);
      FFDBox->Set_Fix_JPlane(mDegree - 1);
      FFDBox->Set_Fix_JPlane(mDegree - 2);
    }
    if ((KPlane_Intersect_A) && (nDegree > 1)) {
      FFDBox->Set_Fix_KPlane(0);
      FFDBox->Set_Fix_KPlane(1);
      FFDBox->Set_Fix_KPlane(2);
    }
    if ((KPlane_Intersect_B) && (nDegree > 1)) {
      FFDBox->Set_Fix_KPlane(nDegree);
      FFDBox->Set_Fix_KPlane(nDegree - 1);
      FFDBox->Set_Fix_KPlane(nDegree - 2);
    }
  }
}

void CSurfaceMovement::UpdateParametricCoord(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                             unsigned short iFFDBox) {
  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint, iSurfacePoints;
  su2double CartCoord[3] = {0.0, 0.0, 0.0}, *CartCoordNew, *CartCoordOld;
  su2double *ParamCoord, *var_coord, ParamCoordGuess[3] = {0.0, 0.0, 0.0};
  su2double MaxDiff, my_MaxDiff = 0.0, Diff;

  /*--- Recompute the parametric coordinates ---*/

  for (iSurfacePoints = 0; iSurfacePoints < FFDBox->GetnSurfacePoint(); iSurfacePoints++) {
    /*--- Get the marker of the surface point ---*/

    iMarker = FFDBox->Get_MarkerIndex(iSurfacePoints);

    if (config->GetMarker_All_DV(iMarker) == YES) {
      /*--- Get the vertex of the surface point ---*/

      iVertex = FFDBox->Get_VertexIndex(iSurfacePoints);
      iPoint = FFDBox->Get_PointIndex(iSurfacePoints);

      /*--- Get the parametric and cartesians coordinates of the
       surface point (they don't mach) ---*/

      ParamCoord = FFDBox->Get_ParametricCoord(iSurfacePoints);

      /*--- Compute and set the cartesian coord using the variation computed
       with the previous deformation ---*/

      var_coord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
      CartCoordOld = geometry->nodes->GetCoord(iPoint);
      for (iDim = 0; iDim < 3; iDim++) CartCoord[iDim] = CartCoordOld[iDim] + var_coord[iDim];
      FFDBox->Set_CartesianCoord(CartCoord, iSurfacePoints);

      /*--- Find the parametric coordinate using as ParamCoordGuess the previous value ---*/

      ParamCoordGuess[0] = ParamCoord[0];
      ParamCoordGuess[1] = ParamCoord[1];
      ParamCoordGuess[2] = ParamCoord[2];
      ParamCoord = FFDBox->GetParametricCoord_Iterative(iPoint, CartCoord, ParamCoordGuess, config);

      /*--- Set the new value of the parametric coordinates ---*/

      FFDBox->Set_ParametricCoord(ParamCoord, iSurfacePoints);

      /*--- Compute the cartesian coordinates using the parametric coordinates
       to check that everything is correct ---*/

      CartCoordNew = FFDBox->EvalCartesianCoord(ParamCoord);

      /*--- Compute max difference between original value and the recomputed value ---*/

      Diff = 0.0;
      for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
        Diff += (CartCoordNew[iDim] - CartCoord[iDim]) * (CartCoordNew[iDim] - CartCoord[iDim]);
      Diff = sqrt(Diff);
      my_MaxDiff = max(my_MaxDiff, Diff);
    }
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&my_MaxDiff, &MaxDiff, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
#else
  MaxDiff = my_MaxDiff;
#endif

  if (rank == MASTER_NODE)
    cout << "Update parametric coord       | FFD box: " << FFDBox->GetTag() << ". Max Diff: " << MaxDiff << "." << endl;
}

void CSurfaceMovement::ApplyDesignVariables(CGeometry* geometry, CConfig* config, CFreeFormDefBox** FFDBox,
                                            unsigned short iFFDBox) {
  unsigned short iDV;

  for (iDV = 0; iDV < config->GetnDV(); iDV++) {
    switch (config->GetDesign_Variable(iDV)) {
      case FFD_CONTROL_POINT_2D:
        SetFFDCPChange_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
      case FFD_CAMBER_2D:
        SetFFDCamber_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
      case FFD_THICKNESS_2D:
        SetFFDThickness_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
      case FFD_CONTROL_POINT:
        SetFFDCPChange(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
      case FFD_NACELLE:
        SetFFDNacelle(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
      case FFD_GULL:
        SetFFDGull(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
      case FFD_TWIST:
        SetFFDTwist(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
      case FFD_ROTATION:
        SetFFDRotation(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
      case FFD_CONTROL_SURFACE:
        SetFFDControl_Surface(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
      case FFD_CAMBER:
        SetFFDCamber(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
      case FFD_THICKNESS:
        SetFFDThickness(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
      case FFD_ANGLE_OF_ATTACK:
        SetFFDAngleOfAttack(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false);
        break;
    }
  }
}

su2double CSurfaceMovement::SetCartesianCoord(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                              unsigned short iFFDBox, bool ResetDef) {
  su2double *CartCoordNew, Diff, my_MaxDiff = 0.0, MaxDiff, *ParamCoord, VarCoord[3] = {0.0, 0.0, 0.0},
                                 CartCoordOld[3] = {0.0, 0.0, 0.0};
  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint, iSurfacePoints;

  bool cylindrical = (config->GetFFD_CoordSystem() == CYLINDRICAL);
  bool spherical = (config->GetFFD_CoordSystem() == SPHERICAL);
  bool polar = (config->GetFFD_CoordSystem() == POLAR);
  unsigned short nDim = geometry->GetnDim();

  /*--- Set to zero all the porints in VarCoord, this is important when we are dealing with different boxes
    because a loop over GetnSurfacePoint is no sufficient ---*/

  if (ResetDef) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
    }
  }

  /*--- Recompute the cartesians coordinates ---*/

  for (iSurfacePoints = 0; iSurfacePoints < FFDBox->GetnSurfacePoint(); iSurfacePoints++) {
    /*--- Get the marker of the surface point ---*/

    iMarker = FFDBox->Get_MarkerIndex(iSurfacePoints);

    if (config->GetMarker_All_DV(iMarker) == YES) {
      /*--- Get the vertex of the surface point ---*/

      iVertex = FFDBox->Get_VertexIndex(iSurfacePoints);
      iPoint = FFDBox->Get_PointIndex(iSurfacePoints);

      /*--- Set to zero the variation of the coordinates ---*/

      geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

      /*--- Get the parametric coordinate of the surface point ---*/

      ParamCoord = FFDBox->Get_ParametricCoord(iSurfacePoints);

      /*--- Compute the new cartesian coordinate, and set the value in
       the FFDBox structure ---*/

      CartCoordNew = FFDBox->EvalCartesianCoord(ParamCoord);

      /*--- If polar coordinates, compute the cartesians from the polar value ---*/

      if (cylindrical) {
        su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
        X_0 = config->GetFFD_Axis(0);
        Y_0 = config->GetFFD_Axis(1);
        Z_0 = config->GetFFD_Axis(2);

        Xbar = CartCoordNew[2];
        Ybar = CartCoordNew[0] * cos(CartCoordNew[1]);
        Zbar = CartCoordNew[0] * sin(CartCoordNew[1]);

        CartCoordNew[0] = Xbar + X_0;
        CartCoordNew[1] = Ybar + Y_0;
        CartCoordNew[2] = Zbar + Z_0;

      } else if (spherical || polar) {
        su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
        X_0 = config->GetFFD_Axis(0);
        Y_0 = config->GetFFD_Axis(1);
        Z_0 = config->GetFFD_Axis(2);

        Xbar = CartCoordNew[0] * cos(CartCoordNew[2]);
        Ybar = CartCoordNew[0] * cos(CartCoordNew[1]) * sin(CartCoordNew[2]);
        Zbar = CartCoordNew[0] * sin(CartCoordNew[1]) * sin(CartCoordNew[2]);

        CartCoordNew[0] = Xbar + X_0;
        CartCoordNew[1] = Ybar + Y_0;
        CartCoordNew[2] = Zbar + Z_0;
      }

      FFDBox->Set_CartesianCoord(CartCoordNew, iSurfacePoints);

      /*--- Get the original cartesian coordinates of the surface point ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        CartCoordOld[iDim] = geometry->nodes->GetCoord(iPoint, iDim);
      }

      /*--- Set the value of the variation of the coordinates ---*/

      Diff = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        VarCoord[iDim] = CartCoordNew[iDim] - CartCoordOld[iDim];
        if ((fabs(VarCoord[iDim]) <= EPS) && (config->GetDirectDiff() != D_DESIGN) && (!config->GetAD_Mode()))
          VarCoord[iDim] = 0.0;
        Diff += (VarCoord[iDim] * VarCoord[iDim]);
      }
      Diff = sqrt(Diff);

      my_MaxDiff = max(my_MaxDiff, Diff);

      /*--- Set the variation of the coordinates ---*/

      geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
    }
  }

  SU2_MPI::Allreduce(&my_MaxDiff, &MaxDiff, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());

  if (rank == MASTER_NODE)
    cout << "Update cartesian coord        | FFD box: " << FFDBox->GetTag() << ". Max Diff: " << MaxDiff << "." << endl;

  return MaxDiff;
}

bool CSurfaceMovement::SetFFDCPChange_2D(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                         CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const {
  su2double movement[3] = {0.0, 0.0, 0.0}, Ampl;
  unsigned short index[3], i, j, iFFDBox, iPlane;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();
  bool polar = (config->GetFFD_CoordSystem() == POLAR);

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/

  if (ResetDef) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }

  design_FFDBox = config->GetFFDTag(iDV);

  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    /*--- Compute deformation ---*/

    /*--- If we have only design value, than this value is the amplitude,
     * otherwise we have a general movement. ---*/

    if (config->GetnDV_Value(iDV) == 1) {
      Ampl = config->GetDV_Value(iDV) * Scale;

      if (polar) {
        movement[0] = config->GetParamDV(iDV, 3) * Ampl;
        movement[1] = 0.0;
        movement[2] = config->GetParamDV(iDV, 4) * Ampl;
      } else {
        movement[0] = config->GetParamDV(iDV, 3) * Ampl;
        movement[1] = config->GetParamDV(iDV, 4) * Ampl;
        movement[2] = 0.0;
      }

    } else {
      if (polar) {
        movement[0] = config->GetDV_Value(iDV, 0);
        movement[1] = 0.0;
        movement[2] = config->GetDV_Value(iDV, 1);
      } else {
        movement[0] = config->GetDV_Value(iDV, 0);
        movement[1] = config->GetDV_Value(iDV, 1);
        movement[2] = 0.0;
      }
    }

    if (polar) {
      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = 0;
      index[2] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
    } else {
      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
      index[2] = 0;
    }

    /*--- Check that it is possible to move the control point ---*/

    for (iPlane = 0; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
      if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) return false;
    }

    for (iPlane = 0; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
      if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
    }

    for (iPlane = 0; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
      if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) return false;
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        FFDBox->SetControlPoints(index, movement);
      }
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1)) {
      for (j = 0; j < FFDBox->GetmOrder(); j++) {
        index[1] = j;
        FFDBox->SetControlPoints(index, movement);
      }
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        for (j = 0; j < FFDBox->GetmOrder(); j++) {
          index[1] = j;
          FFDBox->SetControlPoints(index, movement);
        }
      }
    }
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1)) {
      FFDBox->SetControlPoints(index, movement);
    }

    /*--- Upper surface ---*/

    if (polar)
      index[1] = 1;
    else
      index[2] = 1;

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        FFDBox->SetControlPoints(index, movement);
      }
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1)) {
      for (j = 0; j < FFDBox->GetmOrder(); j++) {
        index[1] = j;
        FFDBox->SetControlPoints(index, movement);
      }
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        for (j = 0; j < FFDBox->GetmOrder(); j++) {
          index[1] = j;
          FFDBox->SetControlPoints(index, movement);
        }
      }
    }
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1)) {
      FFDBox->SetControlPoints(index, movement);
    }
  } else {
    return false;
  }

  return true;
}

bool CSurfaceMovement::SetFFDCPChange(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                      CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const {
  su2double movement[3] = {0.0, 0.0, 0.0}, Ampl;
  unsigned short index[3], i, j, k, iPlane, iFFDBox;
  bool CheckIndex;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/

  if (ResetDef) {
    FFDBox->SetOriginalControlPoints();
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }

  design_FFDBox = config->GetFFDTag(iDV);

  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    /*--- Compute deformation ---*/

    /*--- If we have only design value, than this value is the amplitude,
     * otherwise we have a general movement. ---*/

    if (config->GetnDV_Value(iDV) == 1) {
      Ampl = config->GetDV_Value(iDV) * Scale;

      movement[0] = config->GetParamDV(iDV, 4) * Ampl;
      movement[1] = config->GetParamDV(iDV, 5) * Ampl;
      movement[2] = config->GetParamDV(iDV, 6) * Ampl;

    } else {
      movement[0] = config->GetDV_Value(iDV, 0);
      movement[1] = config->GetDV_Value(iDV, 1);
      movement[2] = config->GetDV_Value(iDV, 2);
    }

    index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
    index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
    index[2] = SU2_TYPE::Int(config->GetParamDV(iDV, 3));

    /*--- Check that it is possible to move the control point ---*/

    for (iPlane = 0; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
      if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) return false;
    }

    for (iPlane = 0; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
      if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
    }

    for (iPlane = 0; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
      if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) return false;
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;

        CheckIndex = true;
        for (iPlane = 0; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
          if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) CheckIndex = false;
        }

        if (CheckIndex) FFDBox->SetControlPoints(index, movement);
      }
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
      for (j = 0; j < FFDBox->GetmOrder(); j++) {
        index[1] = j;

        CheckIndex = true;
        for (iPlane = 0; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
          if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) CheckIndex = false;
        }

        if (CheckIndex) FFDBox->SetControlPoints(index, movement);
      }
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
      for (k = 0; k < FFDBox->GetnOrder(); k++) {
        index[2] = k;

        CheckIndex = true;
        for (iPlane = 0; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
          if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) CheckIndex = false;
        }

        if (CheckIndex) FFDBox->SetControlPoints(index, movement);
      }
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        for (j = 0; j < FFDBox->GetmOrder(); j++) {
          index[1] = j;
          FFDBox->SetControlPoints(index, movement);
        }
      }
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
      for (j = 0; j < FFDBox->GetmOrder(); j++) {
        index[1] = j;
        for (k = 0; k < FFDBox->GetnOrder(); k++) {
          index[2] = k;
          FFDBox->SetControlPoints(index, movement);
        }
      }
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        for (k = 0; k < FFDBox->GetnOrder(); k++) {
          index[2] = k;
          FFDBox->SetControlPoints(index, movement);
        }
      }
    }

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
      FFDBox->SetControlPoints(index, movement);
    }

  } else {
    return false;
  }

  return true;
}

bool CSurfaceMovement::SetFFDGull(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                  CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const {
  su2double movement[3] = {0.0, 0.0, 0.0}, Ampl;
  unsigned short index[3], i, k, iPlane, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/

  if (ResetDef) {
    FFDBox->SetOriginalControlPoints();
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }

  design_FFDBox = config->GetFFDTag(iDV);

  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    /*--- Compute deformation ---*/

    Ampl = config->GetDV_Value(iDV) * Scale;

    movement[0] = 0.0;
    movement[1] = 0.0;
    movement[2] = Ampl;

    /*--- Change the control points ---*/

    index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));

    /*--- Check that it is possible to move the control point ---*/

    for (iPlane = 0; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
      if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
    }

    for (i = 0; i < FFDBox->GetlOrder(); i++) {
      index[0] = i;
      for (k = 0; k < FFDBox->GetnOrder(); k++) {
        index[2] = k;
        FFDBox->SetControlPoints(index, movement);
      }
    }

  } else {
    return false;
  }

  return true;
}

bool CSurfaceMovement::SetFFDNacelle(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                     CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const {
  su2double movement[3] = {0.0, 0.0, 0.0}, Ampl;
  unsigned short index[3], i, j, k, iPlane, iFFDBox, Theta, ThetaMax;
  string design_FFDBox;
  bool SameCP = false;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/

  if (ResetDef) {
    FFDBox->SetOriginalControlPoints();
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }

  design_FFDBox = config->GetFFDTag(iDV);

  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    /*--- Compute deformation ---*/

    Ampl = config->GetDV_Value(iDV) * Scale;

    movement[0] = config->GetParamDV(iDV, 4) * Ampl;
    movement[1] = 0.0;
    movement[2] = config->GetParamDV(iDV, 5) * Ampl;

    index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
    index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
    index[2] = SU2_TYPE::Int(config->GetParamDV(iDV, 3));
    if (index[1] == SU2_TYPE::Int(FFDBox->GetmOrder()) - index[1] - 1) SameCP = true;

    ThetaMax = 2;
    if (SameCP) ThetaMax = 1;

    for (Theta = 0; Theta < ThetaMax; Theta++) {
      if (Theta == 1) index[1] = SU2_TYPE::Int(FFDBox->GetmOrder()) - index[1] - 1;

      /*--- Check that it is possible to move the control point ---*/

      for (iPlane = 0; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
        if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) return false;
      }

      for (iPlane = 0; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
        if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
      }

      for (iPlane = 0; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
        if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) return false;
      }

      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
        for (i = 0; i < FFDBox->GetlOrder(); i++) {
          index[0] = i;
          FFDBox->SetControlPoints(index, movement);
        }
      }

      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
        for (j = 0; j < FFDBox->GetmOrder(); j++) {
          index[1] = j;
          FFDBox->SetControlPoints(index, movement);
        }
      }

      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
        for (k = 0; k < FFDBox->GetnOrder(); k++) {
          index[2] = k;
          FFDBox->SetControlPoints(index, movement);
        }
      }

      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
        for (i = 0; i < FFDBox->GetlOrder(); i++) {
          index[0] = i;
          for (j = 0; j < FFDBox->GetmOrder(); j++) {
            index[1] = j;
            FFDBox->SetControlPoints(index, movement);
          }
        }
      }

      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
        for (j = 0; j < FFDBox->GetmOrder(); j++) {
          index[1] = j;
          for (k = 0; k < FFDBox->GetnOrder(); k++) {
            index[2] = k;
            FFDBox->SetControlPoints(index, movement);
          }
        }
      }

      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
        for (i = 0; i < FFDBox->GetlOrder(); i++) {
          index[0] = i;
          for (k = 0; k < FFDBox->GetnOrder(); k++) {
            index[2] = k;
            FFDBox->SetControlPoints(index, movement);
          }
        }
      }

      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) && (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
        FFDBox->SetControlPoints(index, movement);
      }
    }

  } else {
    return false;
  }

  return true;
}

bool CSurfaceMovement::SetFFDCamber_2D(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                       CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const {
  su2double Ampl, movement[3] = {0.0, 0.0, 0.0};
  unsigned short index[3], kIndex, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/

  if (ResetDef) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }

  design_FFDBox = config->GetFFDTag(iDV);

  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    for (kIndex = 0; kIndex < 2; kIndex++) {
      Ampl = config->GetDV_Value(iDV) * Scale;

      movement[0] = 0.0;
      if (kIndex == 0)
        movement[1] = Ampl;
      else
        movement[1] = Ampl;
      movement[2] = 0.0;

      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = kIndex;
      index[2] = 0;
      FFDBox->SetControlPoints(index, movement);

      index[2] = 1;
      FFDBox->SetControlPoints(index, movement);
    }

  } else {
    return false;
  }

  return true;
}

bool CSurfaceMovement::SetFFDThickness_2D(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                          CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const {
  su2double Ampl, movement[3] = {0.0, 0.0, 0.0};
  unsigned short index[3], kIndex, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/

  if (ResetDef) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }

  design_FFDBox = config->GetFFDTag(iDV);

  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    for (kIndex = 0; kIndex < 2; kIndex++) {
      Ampl = config->GetDV_Value(iDV) * Scale;

      movement[0] = 0.0;
      if (kIndex == 0)
        movement[1] = -Ampl;
      else
        movement[1] = Ampl;
      movement[2] = 0.0;

      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = kIndex;
      index[2] = 0;
      FFDBox->SetControlPoints(index, movement);

      index[2] = 1;
      FFDBox->SetControlPoints(index, movement);
    }

  } else {
    return false;
  }

  return true;
}

bool CSurfaceMovement::SetFFDCamber(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                    CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const {
  su2double Ampl, movement[3] = {0.0, 0.0, 0.0};
  unsigned short index[3], kIndex, iPlane, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/

  if (ResetDef) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }

  design_FFDBox = config->GetFFDTag(iDV);

  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    /*--- Check that it is possible to move the control point ---*/

    for (kIndex = 0; kIndex < 2; kIndex++) {
      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
      index[2] = kIndex;

      for (iPlane = 0; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
        if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) return false;
      }

      for (iPlane = 0; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
        if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
      }

      for (iPlane = 0; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
        if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) return false;
      }
    }

    for (kIndex = 0; kIndex < 2; kIndex++) {
      Ampl = config->GetDV_Value(iDV) * Scale;

      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
      index[2] = kIndex;

      movement[0] = 0.0;
      movement[1] = 0.0;
      if (kIndex == 0)
        movement[2] = Ampl;
      else
        movement[2] = Ampl;

      FFDBox->SetControlPoints(index, movement);
    }

  } else {
    return false;
  }

  return true;
}

void CSurfaceMovement::SetFFDAngleOfAttack(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                           CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) {
  su2double Scale = config->GetOpt_RelaxFactor();

  su2double Ampl = config->GetDV_Value(iDV) * Scale;

  config->SetAoA_Offset(Ampl);
}

bool CSurfaceMovement::SetFFDThickness(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                       CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const {
  su2double Ampl, movement[3] = {0.0, 0.0, 0.0};
  unsigned short index[3], kIndex, iPlane, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/

  if (ResetDef) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }

  design_FFDBox = config->GetFFDTag(iDV);

  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    /*--- Check that it is possible to move the control point ---*/

    for (kIndex = 0; kIndex < 2; kIndex++) {
      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
      index[2] = kIndex;

      for (iPlane = 0; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
        if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) return false;
      }

      for (iPlane = 0; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
        if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
      }

      for (iPlane = 0; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
        if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) return false;
      }
    }

    for (kIndex = 0; kIndex < 2; kIndex++) {
      Ampl = config->GetDV_Value(iDV) * Scale;

      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
      index[2] = kIndex;

      movement[0] = 0.0;
      movement[1] = 0.0;
      if (kIndex == 0)
        movement[2] = -Ampl;
      else
        movement[2] = Ampl;

      FFDBox->SetControlPoints(index, movement);
    }

  } else {
    return false;
  }

  return true;
}

bool CSurfaceMovement::SetFFDTwist(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                   CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const {
  unsigned short iOrder, jOrder, kOrder;
  su2double x, y, z, movement[3], Segment_P0[3], Segment_P1[3], Plane_P0[3], Plane_Normal[3], Variable_P0, Variable_P1,
      Intersection[3], Variable_Interp;
  unsigned short index[3], iPlane, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/

  if (ResetDef) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }

  design_FFDBox = config->GetFFDTag(iDV);

  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    /*--- Check that it is possible to move the control point ---*/

    jOrder = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
    for (iPlane = 0; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
      if (jOrder == FFDBox->Get_Fix_JPlane(iPlane)) return false;
    }

    /*--- Line plane intersection to find the origin of rotation ---*/

    Segment_P0[0] = config->GetParamDV(iDV, 2);
    Segment_P0[1] = config->GetParamDV(iDV, 3);
    Segment_P0[2] = config->GetParamDV(iDV, 4);

    Segment_P1[0] = config->GetParamDV(iDV, 5);
    Segment_P1[1] = config->GetParamDV(iDV, 6);
    Segment_P1[2] = config->GetParamDV(iDV, 7);

    iOrder = 0;
    jOrder = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
    kOrder = 0;
    su2double* coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
    Plane_P0[0] = coord[0];
    Plane_P0[1] = coord[1];
    Plane_P0[2] = coord[2];
    Plane_Normal[0] = 0.0;
    Plane_Normal[1] = 1.0;
    Plane_Normal[2] = 0.0;

    Variable_P0 = 0.0;
    Variable_P1 = 0.0;

    Intersection[0] = 0.0;
    Intersection[1] = 0.0;
    Intersection[2] = 0.0;

    bool result = geometry->SegmentIntersectsPlane(Segment_P0, Segment_P1, Variable_P0, Variable_P1, Plane_P0,
                                                   Plane_Normal, Intersection, Variable_Interp);

    if (result) {
      /*--- xyz-coordinates of a point on the line of rotation. ---*/

      su2double a = Intersection[0];
      su2double b = Intersection[1];
      su2double c = Intersection[2];

      /*--- xyz-coordinate of the line's direction vector. ---*/

      su2double u = Plane_Normal[0];
      su2double v = Plane_Normal[1];
      su2double w = Plane_Normal[2];

      /*--- The angle of rotation is computed based on a characteristic length of the wing,
       otherwise it is difficult to compare with other length based design variables. ---*/

      su2double RefLength = config->GetRefLength();
      su2double theta = atan(config->GetDV_Value(iDV) * Scale / RefLength);

      /*--- An intermediate value used in computations. ---*/

      su2double u2 = u * u;
      su2double v2 = v * v;
      su2double w2 = w * w;
      su2double l2 = u2 + v2 + w2;
      su2double l = sqrt(l2);
      su2double cosT;
      su2double sinT;

      /*--- Change the value of the control point if move is true ---*/

      jOrder = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++)
        for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
          index[0] = iOrder;
          index[1] = jOrder;
          index[2] = kOrder;
          su2double* coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
          x = coord[0];
          y = coord[1];
          z = coord[2];

          cosT = cos(theta);
          sinT = sin(theta);

          movement[0] = a * (v2 + w2) + u * (-b * v - c * w + u * x + v * y + w * z) +
                        (-a * (v2 + w2) + u * (b * v + c * w - v * y - w * z) + (v2 + w2) * x) * cosT +
                        l * (-c * v + b * w - w * y + v * z) * sinT;
          movement[0] = movement[0] / l2 - x;

          movement[1] = b * (u2 + w2) + v * (-a * u - c * w + u * x + v * y + w * z) +
                        (-b * (u2 + w2) + v * (a * u + c * w - u * x - w * z) + (u2 + w2) * y) * cosT +
                        l * (c * u - a * w + w * x - u * z) * sinT;
          movement[1] = movement[1] / l2 - y;

          movement[2] = c * (u2 + v2) + w * (-a * u - b * v + u * x + v * y + w * z) +
                        (-c * (u2 + v2) + w * (a * u + b * v - u * x - v * y) + (u2 + v2) * z) * cosT +
                        l * (-b * u + a * v - v * x + u * y) * sinT;
          movement[2] = movement[2] / l2 - z;

          /*--- Check that it is possible to move the control point ---*/

          for (iPlane = 0; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
            if (iOrder == FFDBox->Get_Fix_IPlane(iPlane)) {
              movement[0] = 0.0;
              movement[1] = 0.0;
              movement[2] = 0.0;
            }
          }

          for (iPlane = 0; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
            if (kOrder == FFDBox->Get_Fix_KPlane(iPlane)) {
              movement[0] = 0.0;
              movement[1] = 0.0;
              movement[2] = 0.0;
            }
          }

          FFDBox->SetControlPoints(index, movement);
        }
    }

  } else {
    return false;
  }

  return true;
}

bool CSurfaceMovement::SetFFDRotation(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                      CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const {
  unsigned short iOrder, jOrder, kOrder;
  su2double movement[3] = {0.0, 0.0, 0.0}, x, y, z;
  unsigned short index[3], iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/

  if (ResetDef) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }

  design_FFDBox = config->GetFFDTag(iDV);

  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    /*--- xyz-coordinates of a point on the line of rotation. ---*/

    su2double a = config->GetParamDV(iDV, 1);
    su2double b = config->GetParamDV(iDV, 2);
    su2double c = config->GetParamDV(iDV, 3);

    /*--- xyz-coordinate of the line's direction vector. ---*/

    su2double u = config->GetParamDV(iDV, 4) - config->GetParamDV(iDV, 1);
    su2double v = config->GetParamDV(iDV, 5) - config->GetParamDV(iDV, 2);
    su2double w = config->GetParamDV(iDV, 6) - config->GetParamDV(iDV, 3);

    /*--- The angle of rotation. ---*/

    su2double theta = config->GetDV_Value(iDV) * Scale * PI_NUMBER / 180.0;

    /*--- An intermediate value used in computations. ---*/

    su2double u2 = u * u;
    su2double v2 = v * v;
    su2double w2 = w * w;
    su2double cosT = cos(theta);
    su2double sinT = sin(theta);
    su2double l2 = u2 + v2 + w2;
    su2double l = sqrt(l2);

    /*--- Change the value of the control point if move is true ---*/

    for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++)
      for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++)
        for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
          index[0] = iOrder;
          index[1] = jOrder;
          index[2] = kOrder;
          su2double* coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
          x = coord[0];
          y = coord[1];
          z = coord[2];
          movement[0] = a * (v2 + w2) + u * (-b * v - c * w + u * x + v * y + w * z) +
                        (-a * (v2 + w2) + u * (b * v + c * w - v * y - w * z) + (v2 + w2) * x) * cosT +
                        l * (-c * v + b * w - w * y + v * z) * sinT;
          movement[0] = movement[0] / l2 - x;

          movement[1] = b * (u2 + w2) + v * (-a * u - c * w + u * x + v * y + w * z) +
                        (-b * (u2 + w2) + v * (a * u + c * w - u * x - w * z) + (u2 + w2) * y) * cosT +
                        l * (c * u - a * w + w * x - u * z) * sinT;
          movement[1] = movement[1] / l2 - y;

          movement[2] = c * (u2 + v2) + w * (-a * u - b * v + u * x + v * y + w * z) +
                        (-c * (u2 + v2) + w * (a * u + b * v - u * x - v * y) + (u2 + v2) * z) * cosT +
                        l * (-b * u + a * v - v * x + u * y) * sinT;
          movement[2] = movement[2] / l2 - z;

          FFDBox->SetControlPoints(index, movement);
        }
  } else {
    return false;
  }

  return true;
}

bool CSurfaceMovement::SetFFDControl_Surface(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                                             CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const {
  unsigned short iOrder, jOrder, kOrder;
  su2double movement[3] = {0.0, 0.0, 0.0}, x, y, z;
  unsigned short index[3], iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/

  if (ResetDef) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }

  design_FFDBox = config->GetFFDTag(iDV);

  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    /*--- xyz-coordinates of a point on the line of rotation. ---*/

    su2double a = config->GetParamDV(iDV, 1);
    su2double b = config->GetParamDV(iDV, 2);
    su2double c = config->GetParamDV(iDV, 3);

    /*--- xyz-coordinate of the line's direction vector. ---*/

    su2double u = config->GetParamDV(iDV, 4) - config->GetParamDV(iDV, 1);
    su2double v = config->GetParamDV(iDV, 5) - config->GetParamDV(iDV, 2);
    su2double w = config->GetParamDV(iDV, 6) - config->GetParamDV(iDV, 3);

    /*--- The angle of rotation. ---*/

    su2double theta = -config->GetDV_Value(iDV) * Scale * PI_NUMBER / 180.0;

    /*--- An intermediate value used in computations. ---*/

    su2double u2 = u * u;
    su2double v2 = v * v;
    su2double w2 = w * w;
    su2double cosT = cos(theta);
    su2double sinT = sin(theta);
    su2double l2 = u2 + v2 + w2;
    su2double l = sqrt(l2);

    /*--- Change the value of the control point if move is true ---*/

    for (iOrder = 0; iOrder < FFDBox->GetlOrder() - 2; iOrder++)
      for (jOrder = 2; jOrder < FFDBox->GetmOrder() - 2; jOrder++)
        for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
          index[0] = iOrder;
          index[1] = jOrder;
          index[2] = kOrder;
          su2double* coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
          x = coord[0];
          y = coord[1];
          z = coord[2];
          movement[0] = a * (v2 + w2) + u * (-b * v - c * w + u * x + v * y + w * z) +
                        (-a * (v2 + w2) + u * (b * v + c * w - v * y - w * z) + (v2 + w2) * x) * cosT +
                        l * (-c * v + b * w - w * y + v * z) * sinT;
          movement[0] = movement[0] / l2 - x;

          movement[1] = b * (u2 + w2) + v * (-a * u - c * w + u * x + v * y + w * z) +
                        (-b * (u2 + w2) + v * (a * u + c * w - u * x - w * z) + (u2 + w2) * y) * cosT +
                        l * (c * u - a * w + w * x - u * z) * sinT;
          movement[1] = movement[1] / l2 - y;

          movement[2] = c * (u2 + v2) + w * (-a * u - b * v + u * x + v * y + w * z) +
                        (-c * (u2 + v2) + w * (a * u + b * v - u * x - v * y) + (u2 + v2) * z) * cosT +
                        l * (-b * u + a * v - v * x + u * y) * sinT;
          movement[2] = movement[2] / l2 - z;

          FFDBox->SetControlPoints(index, movement);
        }
  } else {
    return false;
  }

  return true;
}

void CSurfaceMovement::SetAngleOfAttack(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef) {
  su2double Scale = config->GetOpt_RelaxFactor();
  su2double Ampl = config->GetDV_Value(iDV) * Scale;
  config->SetAoA_Offset(Ampl);
}

void CSurfaceMovement::SetHicksHenne(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0, 0.0, 0.0}, VarCoord_[3] = {0.0, 0.0, 0.0}, *Coord_, *Normal_, ek, fk,
            Coord[3] = {0.0, 0.0, 0.0}, Normal[3] = {0.0, 0.0, 0.0}, TPCoord[2] = {0.0, 0.0}, LPCoord[2] = {0.0, 0.0},
            Distance, Chord, AoA, ValCos, ValSin;

  bool upper = true;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

  if ((iDV == 0) || (ResetDef)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0;
        VarCoord[1] = 0.0;
        VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }

  /*--- Compute the angle of attack to apply the deformation ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      Coord_ = boundary->vertex[iMarker][0]->GetCoord();
      TPCoord[0] = Coord_[0];
      TPCoord[1] = Coord_[1];
      for (iVertex = 1; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        if (Coord_[0] > TPCoord[0]) {
          TPCoord[0] = Coord_[0];
          TPCoord[1] = Coord_[1];
        }
      }
    }
  }

#ifdef HAVE_MPI

  int iProcessor, nProcessor = size;
  su2double *Buffer_Send_Coord, *Buffer_Receive_Coord;

  Buffer_Receive_Coord = new su2double[nProcessor * 2];
  Buffer_Send_Coord = new su2double[2];

  Buffer_Send_Coord[0] = TPCoord[0];
  Buffer_Send_Coord[1] = TPCoord[1];

  SU2_MPI::Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, SU2_MPI::GetComm());

  TPCoord[0] = Buffer_Receive_Coord[0];
  TPCoord[1] = Buffer_Receive_Coord[1];
  for (iProcessor = 1; iProcessor < nProcessor; iProcessor++) {
    Coord[0] = Buffer_Receive_Coord[iProcessor * 2 + 0];
    Coord[1] = Buffer_Receive_Coord[iProcessor * 2 + 1];
    if (Coord[0] > TPCoord[0]) {
      TPCoord[0] = Coord[0];
      TPCoord[1] = Coord[1];
    }
  }

  delete[] Buffer_Send_Coord;
  delete[] Buffer_Receive_Coord;

#endif

  Chord = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        Distance = sqrt(pow(Coord_[0] - TPCoord[0], 2.0) + pow(Coord_[1] - TPCoord[1], 2.0));
        if (Chord < Distance) {
          Chord = Distance;
          LPCoord[0] = Coord_[0];
          LPCoord[1] = Coord_[1];
        }
      }
    }
  }

#ifdef HAVE_MPI

  Buffer_Receive_Coord = new su2double[nProcessor * 2];
  Buffer_Send_Coord = new su2double[2];

  Buffer_Send_Coord[0] = LPCoord[0];
  Buffer_Send_Coord[1] = LPCoord[1];

  SU2_MPI::Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, SU2_MPI::GetComm());

  Chord = 0.0;
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    Coord[0] = Buffer_Receive_Coord[iProcessor * 2 + 0];
    Coord[1] = Buffer_Receive_Coord[iProcessor * 2 + 1];
    Distance = sqrt(pow(Coord[0] - TPCoord[0], 2.0) + pow(Coord[1] - TPCoord[1], 2.0));
    if (Chord < Distance) {
      Chord = Distance;
      LPCoord[0] = Coord[0];
      LPCoord[1] = Coord[1];
    }
  }

  delete[] Buffer_Send_Coord;
  delete[] Buffer_Receive_Coord;

#endif

  AoA = atan((LPCoord[1] - TPCoord[1]) / (TPCoord[0] - LPCoord[0])) * 180 / PI_NUMBER;

  /*--- WARNING: AoA currently overwritten to zero. ---*/
  AoA = 0.0;

  /*--- Perform multiple airfoil deformation ---*/

  su2double Ampl = config->GetDV_Value(iDV) * Scale;
  su2double xk = config->GetParamDV(iDV, 1);
  const su2double t2 = 3.0;

  if (config->GetParamDV(iDV, 0) == NO) {
    upper = false;
  }
  if (config->GetParamDV(iDV, 0) == YES) {
    upper = true;
  }

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0;
      VarCoord[1] = 0.0;
      VarCoord[2] = 0.0;

      if (config->GetMarker_All_DV(iMarker) == YES) {
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        Normal_ = boundary->vertex[iMarker][iVertex]->GetNormal();

        /*--- The Hicks Henne bump functions should be applied to a basic airfoil without AoA,
         and unitary chord, a tranformation is required ---*/

        ValCos = cos(AoA * PI_NUMBER / 180.0);
        ValSin = sin(AoA * PI_NUMBER / 180.0);

        Coord[0] = Coord_[0] * ValCos - Coord_[1] * ValSin;
        Coord[0] = max(0.0, Coord[0]);  // Coord x should be always positive
        Coord[1] = Coord_[1] * ValCos + Coord_[0] * ValSin;

        Normal[0] = Normal_[0] * ValCos - Normal_[1] * ValSin;
        Normal[1] = Normal_[1] * ValCos + Normal_[0] * ValSin;

        /*--- Bump computation ---*/

        ek = log10(0.5) / log10(xk);
        if (Coord[0] > 10 * EPS)
          fk = pow(sin(PI_NUMBER * pow(Coord[0], ek)), t2);
        else
          fk = 0.0;

        /*--- Upper and lower surface ---*/

        if ((upper) && (Normal[1] > 0)) {
          VarCoord[1] = Ampl * fk;
        }
        if ((!upper) && (Normal[1] < 0)) {
          VarCoord[1] = -Ampl * fk;
        }
      }

      /*--- Apply the transformation to the coordinate variation ---*/

      ValCos = cos(-AoA * PI_NUMBER / 180.0);
      ValSin = sin(-AoA * PI_NUMBER / 180.0);

      VarCoord_[0] = VarCoord[0] * ValCos - VarCoord[1] * ValSin;
      VarCoord_[1] = VarCoord[1] * ValCos + VarCoord[0] * ValSin;

      boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord_);
    }
  }
}

void CSurfaceMovement::SetSurface_Bump(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0, 0.0, 0.0}, ek, fk, *Coord, xCoord;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

  if ((iDV == 0) || (ResetDef)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0;
        VarCoord[1] = 0.0;
        VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }

  /*--- Perform multiple airfoil deformation ---*/

  su2double Ampl = config->GetDV_Value(iDV) * Scale;
  su2double x_start = config->GetParamDV(iDV, 0);
  su2double x_end = config->GetParamDV(iDV, 1);
  su2double BumpSize = x_end - x_start;
  su2double BumpLoc = x_start;
  su2double xk = config->GetParamDV(iDV, 2);
  const su2double t2 = 3.0;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0;
      VarCoord[1] = 0.0;
      VarCoord[2] = 0.0;

      if (config->GetMarker_All_DV(iMarker) == YES) {
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();

        xCoord = (Coord[0] - BumpLoc);
        ek = log10(0.5) / log10((xk - BumpLoc + EPS) / BumpSize);
        if (xCoord > 0.0)
          fk = pow(sin(PI_NUMBER * pow((xCoord + EPS) / BumpSize, ek)), t2);
        else
          fk = 0.0;

        if ((xCoord <= 0.0) || (xCoord >= BumpSize))
          VarCoord[1] = 0.0;
        else {
          VarCoord[1] = Ampl * fk;
        }
      }

      boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
    }
  }
}

void CSurfaceMovement::SetCST(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0, 0.0, 0.0}, VarCoord_[3] = {0.0, 0.0, 0.0}, *Coord_, *Normal_, fk,
            Coord[3] = {0.0, 0.0, 0.0}, Normal[3] = {0.0, 0.0, 0.0}, TPCoord[2] = {0.0, 0.0}, LPCoord[2] = {0.0, 0.0},
            Distance, Chord, AoA, ValCos, ValSin;

  bool upper = true;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

  if ((iDV == 0) || (ResetDef)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0;
        VarCoord[1] = 0.0;
        VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }

  /*--- Compute the angle of attack to apply the deformation ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      Coord_ = boundary->vertex[iMarker][0]->GetCoord();
      TPCoord[0] = Coord_[0];
      TPCoord[1] = Coord_[1];
      for (iVertex = 1; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        if (Coord_[0] > TPCoord[0]) {
          TPCoord[0] = Coord_[0];
          TPCoord[1] = Coord_[1];
        }
      }
    }
  }

#ifdef HAVE_MPI

  int iProcessor, nProcessor = size;
  su2double *Buffer_Send_Coord, *Buffer_Receive_Coord;

  Buffer_Receive_Coord = new su2double[nProcessor * 2];
  Buffer_Send_Coord = new su2double[2];

  Buffer_Send_Coord[0] = TPCoord[0];
  Buffer_Send_Coord[1] = TPCoord[1];

  SU2_MPI::Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, SU2_MPI::GetComm());

  TPCoord[0] = Buffer_Receive_Coord[0];
  TPCoord[1] = Buffer_Receive_Coord[1];
  for (iProcessor = 1; iProcessor < nProcessor; iProcessor++) {
    Coord[0] = Buffer_Receive_Coord[iProcessor * 2 + 0];
    Coord[1] = Buffer_Receive_Coord[iProcessor * 2 + 1];
    if (Coord[0] > TPCoord[0]) {
      TPCoord[0] = Coord[0];
      TPCoord[1] = Coord[1];
    }
  }

  delete[] Buffer_Send_Coord;
  delete[] Buffer_Receive_Coord;

#endif

  Chord = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        Distance = sqrt(pow(Coord_[0] - TPCoord[0], 2.0) + pow(Coord_[1] - TPCoord[1], 2.0));
        if (Chord < Distance) {
          Chord = Distance;
          LPCoord[0] = Coord_[0];
          LPCoord[1] = Coord_[1];
        }
      }
    }
  }

#ifdef HAVE_MPI

  Buffer_Receive_Coord = new su2double[nProcessor * 2];
  Buffer_Send_Coord = new su2double[2];

  Buffer_Send_Coord[0] = LPCoord[0];
  Buffer_Send_Coord[1] = LPCoord[1];

  SU2_MPI::Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, SU2_MPI::GetComm());

  Chord = 0.0;
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    Coord[0] = Buffer_Receive_Coord[iProcessor * 2 + 0];
    Coord[1] = Buffer_Receive_Coord[iProcessor * 2 + 1];
    Distance = sqrt(pow(Coord[0] - TPCoord[0], 2.0) + pow(Coord[1] - TPCoord[1], 2.0));
    if (Chord < Distance) {
      Chord = Distance;
      LPCoord[0] = Coord[0];
      LPCoord[1] = Coord[1];
    }
  }

  delete[] Buffer_Send_Coord;
  delete[] Buffer_Receive_Coord;

#endif

  AoA = atan((LPCoord[1] - TPCoord[1]) / (TPCoord[0] - LPCoord[0])) * 180 / PI_NUMBER;

  /*--- WARNING: AoA currently overwritten to zero. ---*/
  AoA = 0.0;

  /*--- Perform multiple airfoil deformation ---*/

  su2double Ampl = config->GetDV_Value(iDV) * Scale;
  su2double KulfanNum = config->GetParamDV(iDV, 1) - 1.0;
  su2double maxKulfanNum = config->GetParamDV(iDV, 2) - 1.0;
  if (KulfanNum < 0) {
    std::cout << "Warning: Kulfan number should be greater than 1." << std::endl;
  }
  if (KulfanNum > maxKulfanNum) {
    std::cout << "Warning: Kulfan number should be less than provided maximum." << std::endl;
  }

  if (config->GetParamDV(iDV, 0) == NO) {
    upper = false;
  }
  if (config->GetParamDV(iDV, 0) == YES) {
    upper = true;
  }

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0;
      VarCoord[1] = 0.0;
      VarCoord[2] = 0.0;

      if (config->GetMarker_All_DV(iMarker) == YES) {
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        Normal_ = boundary->vertex[iMarker][iVertex]->GetNormal();

        /*--- The CST functions should be applied to a basic airfoil without AoA,
         and unitary chord, a tranformation is required ---*/

        ValCos = cos(AoA * PI_NUMBER / 180.0);
        ValSin = sin(AoA * PI_NUMBER / 180.0);

        Coord[0] = Coord_[0] * ValCos - Coord_[1] * ValSin;
        Coord[0] = max(0.0, Coord[0]);  // Coord x should be always positive
        Coord[1] = Coord_[1] * ValCos + Coord_[0] * ValSin;

        Normal[0] = Normal_[0] * ValCos - Normal_[1] * ValSin;
        Normal[1] = Normal_[1] * ValCos + Normal_[0] * ValSin;

        /*--- CST computation ---*/
        su2double fact_n = 1;
        su2double fact_cst = 1;
        su2double fact_cst_n = 1;

        for (int i = 1; i <= maxKulfanNum; i++) {
          fact_n = fact_n * i;
        }
        for (int i = 1; i <= KulfanNum; i++) {
          fact_cst = fact_cst * i;
        }
        for (int i = 1; i <= maxKulfanNum - KulfanNum; i++) {
          fact_cst_n = fact_cst_n * i;
        }

        // CST method only for 2D NACA type airfoils
        su2double N1, N2;
        N1 = 0.5;
        N2 = 1.0;

        /*--- Upper and lower surface change in coordinates based on CST equations by Kulfan et. al
         * (www.brendakulfan.com/docs/CST3.pdf)  ---*/
        fk = pow(Coord[0], N1) * pow((1 - Coord[0]), N2) * fact_n / (fact_cst * (fact_cst_n)) *
             pow(Coord[0], KulfanNum) * pow((1 - Coord[0]), (maxKulfanNum - (KulfanNum)));

        if ((upper) && (Normal[1] > 0)) {
          VarCoord[1] = Ampl * fk;
        }

        if ((!upper) && (Normal[1] < 0)) {
          VarCoord[1] = Ampl * fk;
        }
      }

      /*--- Apply the transformation to the coordinate variation ---*/

      ValCos = cos(-AoA * PI_NUMBER / 180.0);
      ValSin = sin(-AoA * PI_NUMBER / 180.0);

      VarCoord_[0] = VarCoord[0] * ValCos - VarCoord[1] * ValSin;
      VarCoord_[1] = VarCoord[1] * ValCos + VarCoord[0] * ValSin;

      boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord_);
    }
  }
}

void CSurfaceMovement::SetRotation(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0, 0.0, 0.0}, *Coord;
  su2double movement[3] = {0.0, 0.0, 0.0}, x, y, z;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

  if ((iDV == 0) || (ResetDef)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0;
        VarCoord[1] = 0.0;
        VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }

  /*--- xyz-coordinates of a point on the line of rotation. */

  su2double a = config->GetParamDV(iDV, 0);
  su2double b = config->GetParamDV(iDV, 1);
  su2double c = 0.0;
  if (boundary->GetnDim() == 3) c = config->GetParamDV(0, 2);

  /*--- xyz-coordinate of the line's direction vector. ---*/

  su2double u = config->GetParamDV(iDV, 3) - config->GetParamDV(iDV, 0);
  su2double v = config->GetParamDV(iDV, 4) - config->GetParamDV(iDV, 1);
  su2double w = 1.0;
  if (boundary->GetnDim() == 3) w = config->GetParamDV(iDV, 5) - config->GetParamDV(iDV, 2);

  /*--- The angle of rotation. ---*/

  su2double theta = config->GetDV_Value(iDV) * Scale * PI_NUMBER / 180.0;

  /*--- An intermediate value used in computations. ---*/

  su2double u2 = u * u;
  su2double v2 = v * v;
  su2double w2 = w * w;
  su2double cosT = cos(theta);
  su2double sinT = sin(theta);
  su2double l2 = u2 + v2 + w2;
  su2double l = sqrt(l2);

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0;
      VarCoord[1] = 0.0;
      VarCoord[2] = 0.0;
      if (config->GetMarker_All_DV(iMarker) == YES) {
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
        x = Coord[0];
        y = Coord[1];
        z = Coord[2];

        movement[0] = a * (v2 + w2) + u * (-b * v - c * w + u * x + v * y + w * z) +
                      (-a * (v2 + w2) + u * (b * v + c * w - v * y - w * z) + (v2 + w2) * x) * cosT +
                      l * (-c * v + b * w - w * y + v * z) * sinT;
        movement[0] = movement[0] / l2 - x;

        movement[1] = b * (u2 + w2) + v * (-a * u - c * w + u * x + v * y + w * z) +
                      (-b * (u2 + w2) + v * (a * u + c * w - u * x - w * z) + (u2 + w2) * y) * cosT +
                      l * (c * u - a * w + w * x - u * z) * sinT;
        movement[1] = movement[1] / l2 - y;

        movement[2] = c * (u2 + v2) + w * (-a * u - b * v + u * x + v * y + w * z) +
                      (-c * (u2 + v2) + w * (a * u + b * v - u * x - v * y) + (u2 + v2) * z) * cosT +
                      l * (-b * u + a * v - v * x + u * y) * sinT;
        if (boundary->GetnDim() == 3)
          movement[2] = movement[2] / l2 - z;
        else
          movement[2] = 0.0;

        VarCoord[0] = movement[0];
        VarCoord[1] = movement[1];
        if (boundary->GetnDim() == 3) VarCoord[2] = movement[2];
      }
      boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
    }
}

void CSurfaceMovement::SetTranslation(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0, 0.0, 0.0};
  su2double Scale = config->GetOpt_RelaxFactor();
  su2double Ampl = config->GetDV_Value(iDV) * Scale;

  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

  if ((iDV == 0) || (ResetDef)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0;
        VarCoord[1] = 0.0;
        VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }

  su2double xDispl = config->GetParamDV(iDV, 0);
  su2double yDispl = config->GetParamDV(iDV, 1);
  su2double zDispl = 0;
  if (boundary->GetnDim() == 3) zDispl = config->GetParamDV(iDV, 2);

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0;
      VarCoord[1] = 0.0;
      VarCoord[2] = 0.0;
      if (config->GetMarker_All_DV(iMarker) == YES) {
        VarCoord[0] = Ampl * xDispl;
        VarCoord[1] = Ampl * yDispl;
        if (boundary->GetnDim() == 3) VarCoord[2] = Ampl * zDispl;
      }
      boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
    }
}

void CSurfaceMovement::SetScale(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0, 0.0, 0.0}, x, y, z, *Coord;
  su2double Scale = config->GetOpt_RelaxFactor();
  su2double Ampl = config->GetDV_Value(iDV) * Scale;

  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

  if ((iDV == 0) || (ResetDef)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0;
        VarCoord[1] = 0.0;
        VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0;
      VarCoord[1] = 0.0;
      VarCoord[2] = 0.0;
      if (config->GetMarker_All_DV(iMarker) == YES) {
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
        x = Coord[0];
        y = Coord[1];
        z = Coord[2];
        VarCoord[0] = (Ampl - 1.0) * x;
        VarCoord[1] = (Ampl - 1.0) * y;
        if (boundary->GetnDim() == 3) VarCoord[2] = (Ampl - 1.0) * z;
      }
      boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
    }
}

void CSurfaceMovement::AeroelasticDeform(CGeometry* geometry, CConfig* config, unsigned long TimeIter,
                                         unsigned short iMarker, unsigned short iMarker_Monitoring,
                                         vector<su2double>& displacements) {
  /* The sign conventions of these are those of the Typical Section Wing Model, below the signs are corrected */
  su2double dh = -displacements[0];      // relative plunge
  su2double dalpha = -displacements[1];  // relative pitch
  su2double dh_x, dh_y;
  su2double Center[2];
  unsigned short iDim;
  su2double Lref = config->GetLength_Ref();
  su2double* Coord;
  unsigned long iPoint, iVertex;
  su2double x_new, y_new;
  su2double VarCoord[3];
  string Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);

  /*--- Calculate the plunge displacement for the Typical Section Wing Model taking into account rotation ---*/
  if (config->GetKind_GridMovement() == AEROELASTIC_RIGID_MOTION) {
    su2double Omega, dt, psi;
    dt = config->GetDelta_UnstTimeND();
    Omega = config->GetRotation_Rate(2) / config->GetOmega_Ref();
    psi = Omega * (dt * TimeIter);

    /*--- Correct for the airfoil starting position (This is hardcoded in here) ---*/
    if (Monitoring_Tag == "Airfoil1") {
      psi = psi + 0.0;
    } else if (Monitoring_Tag == "Airfoil2") {
      psi = psi + 2.0 / 3.0 * PI_NUMBER;
    } else if (Monitoring_Tag == "Airfoil3") {
      psi = psi + 4.0 / 3.0 * PI_NUMBER;
    } else
      cout << "WARNING: There is a marker that we are monitoring that doesn't match the values hardcoded above!"
           << endl;

    dh_x = -dh * sin(psi);
    dh_y = dh * cos(psi);

  } else {
    dh_x = 0;
    dh_y = dh;
  }

  /*--- Pitching origin from config. ---*/

  Center[0] = config->GetRefOriginMoment_X(iMarker_Monitoring);
  Center[1] = config->GetRefOriginMoment_Y(iMarker_Monitoring);

  for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
    /*--- Coordinates of the current point ---*/
    Coord = geometry->nodes->GetCoord(iPoint);

    /*--- Calculate non-dim. position from rotation center ---*/
    su2double r[2] = {0, 0};
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++) r[iDim] = (Coord[iDim] - Center[iDim]) / Lref;

    /*--- Compute delta of transformed point coordinates ---*/
    // The deltas are needed for the FEA grid deformation Method.
    // rotation contribution - previous position + plunging contribution
    x_new = cos(dalpha) * r[0] - sin(dalpha) * r[1] - r[0] + dh_x;
    y_new = sin(dalpha) * r[0] + cos(dalpha) * r[1] - r[1] + dh_y;

    VarCoord[0] = x_new;
    VarCoord[1] = y_new;
    VarCoord[2] = 0.0;

    /*--- Store new delta node locations for the surface ---*/
    geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
  }
  /*--- Set the elastic axis to the new location after incrementing the position with the plunge ---*/
  config->SetRefOriginMoment_X(iMarker_Monitoring, Center[0] + dh_x);
  config->SetRefOriginMoment_Y(iMarker_Monitoring, Center[1] + dh_y);
}

void CSurfaceMovement::SetBoundary_Flutter3D(CGeometry* geometry, CConfig* config, CFreeFormDefBox** FFDBox,
                                             unsigned long iter, unsigned short iZone) {
  su2double omega, deltaT;
  su2double alpha, alpha_new, alpha_old;
  su2double time_new, time_old;
  su2double Omega[3], Ampl[3];
  su2double DEG2RAD = PI_NUMBER / 180.0;
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());
  unsigned short iDim = 0;

  /*--- Retrieve values from the config file ---*/

  deltaT = config->GetDelta_UnstTimeND();

  /*--- Pitching origin, frequency, and amplitude from config. ---*/

  for (iDim = 0; iDim < 3; iDim++) {
    Omega[iDim] = config->GetPitching_Omega(iDim) / config->GetOmega_Ref();
    Ampl[iDim] = config->GetPitching_Ampl(iDim) * DEG2RAD;
  }

  /*--- Compute delta time based on physical time step ---*/

  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/

    unsigned long nFlowIter = config->GetnTime_Iter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<su2double>(directIter) * deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<su2double>(directIter) + 1.0) * deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/

    time_new = static_cast<su2double>(iter) * deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<su2double>(iter) - 1.0) * deltaT;
  }

  /*--- Update the pitching angle at this time step. Flip sign for
   nose-up positive convention. ---*/

  omega = Omega[2];
  alpha_new = Ampl[2] * sin(omega * time_new);
  alpha_old = Ampl[2] * sin(omega * time_old);
  alpha = (1E-10 + (alpha_new - alpha_old)) * (-PI_NUMBER / 180.0);

  if (rank == MASTER_NODE) cout << "New dihedral angle (alpha): " << alpha_new / DEG2RAD << " degrees." << endl;

  unsigned short iOrder, jOrder, kOrder;
  short iFFDBox;
  su2double movement[3] = {0.0, 0.0, 0.0};
  bool* move = new bool[nFFDBox];
  auto* index = new unsigned short[3];

  move[0] = true;
  move[1] = true;
  move[2] = true;

  /*--- Change the value of the control point if move is true ---*/

  for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
    if (move[iFFDBox])
      for (iOrder = 0; iOrder < FFDBox[iFFDBox]->GetlOrder(); iOrder++)
        for (jOrder = 0; jOrder < FFDBox[iFFDBox]->GetmOrder(); jOrder++)
          for (kOrder = 0; kOrder < FFDBox[iFFDBox]->GetnOrder(); kOrder++) {
            index[0] = iOrder;
            index[1] = jOrder;
            index[2] = kOrder;
            su2double* coord = FFDBox[iFFDBox]->GetCoordControlPoints(iOrder, jOrder, kOrder);
            movement[0] = 0.0;
            movement[1] = 0.0;
            movement[2] = coord[1] * tan(alpha);
            FFDBox[iFFDBox]->SetControlPoints(index, movement);
          }

  /*--- Recompute cartesian coordinates using the new control points position ---*/

  for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, false);

  delete[] index;
  delete[] move;
}

void CSurfaceMovement::SetExternal_Deformation(CGeometry* geometry, CConfig* config, unsigned short iZone,
                                               unsigned long iter) {
  /*--- Local variables ---*/

  unsigned short iDim, nDim;
  unsigned long iPoint = 0, flowIter = 0;
  unsigned long jPoint, GlobalIndex;
  su2double VarCoord[3], *Coord_Old = nullptr, *Coord_New = nullptr, Center[3] = {0.0, 0.0, 0.0};
  su2double Lref = config->GetLength_Ref();
  su2double NewCoord[3] = {0.0, 0.0, 0.0}, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  su2double r[3] = {0.0, 0.0, 0.0}, rotCoord[3] = {0.0, 0.0, 0.0};
  unsigned long iVertex;
  unsigned short iMarker;
  char buffer[50];
  string DV_Filename, UnstExt, text_line;
  ifstream surface_positions;
  bool unsteady = config->GetTime_Marching() != TIME_MARCHING::STEADY;
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());

  /*--- Load stuff from config ---*/

  nDim = geometry->GetnDim();
  DV_Filename = config->GetDV_Filename();

  /*--- Set the extension for the correct unsteady mesh motion file ---*/

  if (unsteady) {
    if (adjoint) {
      /*--- For the unsteady adjoint, we integrate backwards through
       physical time, so perform mesh motion in reverse. ---*/
      unsigned long nFlowIter = config->GetnTime_Iter() - 1;
      flowIter = nFlowIter - iter;
      unsigned short lastindex = DV_Filename.find_last_of('.');
      DV_Filename = DV_Filename.substr(0, lastindex);
      if ((SU2_TYPE::Int(flowIter) >= 0) && (SU2_TYPE::Int(flowIter) < 10))
        SPRINTF(buffer, "_0000%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 10) && (SU2_TYPE::Int(flowIter) < 100))
        SPRINTF(buffer, "_000%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 100) && (SU2_TYPE::Int(flowIter) < 1000))
        SPRINTF(buffer, "_00%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 1000) && (SU2_TYPE::Int(flowIter) < 10000))
        SPRINTF(buffer, "_0%d.dat", SU2_TYPE::Int(flowIter));
      if (SU2_TYPE::Int(flowIter) >= 10000) SPRINTF(buffer, "_%d.dat", SU2_TYPE::Int(flowIter));
      UnstExt = string(buffer);
      DV_Filename.append(UnstExt);
    } else {
      /*--- Forward time for the direct problem ---*/
      flowIter = iter;
      unsigned short lastindex = DV_Filename.find_last_of('.');
      DV_Filename = DV_Filename.substr(0, lastindex);
      if ((SU2_TYPE::Int(flowIter) >= 0) && (SU2_TYPE::Int(flowIter) < 10))
        SPRINTF(buffer, "_0000%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 10) && (SU2_TYPE::Int(flowIter) < 100))
        SPRINTF(buffer, "_000%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 100) && (SU2_TYPE::Int(flowIter) < 1000))
        SPRINTF(buffer, "_00%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 1000) && (SU2_TYPE::Int(flowIter) < 10000))
        SPRINTF(buffer, "_0%d.dat", SU2_TYPE::Int(flowIter));
      if (SU2_TYPE::Int(flowIter) >= 10000) SPRINTF(buffer, "_%d.dat", SU2_TYPE::Int(flowIter));
      UnstExt = string(buffer);
      DV_Filename.append(UnstExt);
    }

    if (rank == MASTER_NODE)
      cout << "Reading in the arbitrary mesh motion from direct iteration " << flowIter << "." << endl;
  }

  /*--- Open the motion file ---*/

  surface_positions.open(DV_Filename.data(), ios::in);

  /*--- Throw error if there is no file ---*/

  if (surface_positions.fail()) {
    SU2_MPI::Error(string("There is no surface positions file ") + DV_Filename, CURRENT_FUNCTION);
  }

  /*--- Read in and store the new mesh node locations ---*/

  while (getline(surface_positions, text_line)) {
    istringstream point_line(text_line);
    if (nDim == 2) point_line >> iPoint >> NewCoord[0] >> NewCoord[1];
    if (nDim == 3) point_line >> iPoint >> NewCoord[0] >> NewCoord[1] >> NewCoord[2];
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_DV(iMarker) == YES && config->GetKind_SU2() == SU2_COMPONENT::SU2_DEF) ||
          (config->GetMarker_All_Moving(iMarker) == YES && config->GetKind_SU2() == SU2_COMPONENT::SU2_CFD)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex = geometry->nodes->GetGlobalIndex(jPoint);
          if (GlobalIndex == iPoint) {
            geometry->vertex[iMarker][iVertex]->SetVarCoord(NewCoord);
            break;
          }
        }
      }
    }
  }

  /*--- Close the surface positions file ---*/

  surface_positions.close();

  /*--- If rotating as well, prepare the rotation matrix ---*/

  if (config->GetKind_GridMovement() == EXTERNAL_ROTATION) {
    /*--- Variables needed only for rotation ---*/

    su2double Omega[3], dt;
    su2double dtheta, dphi, dpsi, cosTheta, sinTheta;
    su2double cosPhi, sinPhi, cosPsi, sinPsi;

    /*--- Center of rotation & angular velocity vector from config ---*/
    Center[0] = config->GetMotion_Origin(0);
    Center[1] = config->GetMotion_Origin(1);
    Center[2] = config->GetMotion_Origin(2);

    /*--- Angular velocity vector from config ---*/

    dt = static_cast<su2double>(iter) * config->GetDelta_UnstTimeND();
    Omega[0] = config->GetRotation_Rate(0);
    Omega[1] = config->GetRotation_Rate(1);
    Omega[2] = config->GetRotation_Rate(2);

    /*--- For the unsteady adjoint, use reverse time ---*/
    if (adjoint) {
      /*--- Set the first adjoint mesh position to the final direct one ---*/
      if (iter == 0) dt = ((su2double)config->GetnTime_Iter() - 1) * dt;
      /*--- Reverse the rotation direction for the adjoint ---*/
      else
        dt = -1.0 * dt;
    } else {
      /*--- No rotation at all for the first direct solution ---*/
      if (iter == 0) dt = 0;
    }

    /*--- Compute delta change in the angle about the x, y, & z axes. ---*/

    dtheta = Omega[0] * dt;
    dphi = Omega[1] * dt;
    dpsi = Omega[2] * dt;

    /*--- Store angles separately for clarity. Compute sines/cosines. ---*/

    cosTheta = cos(dtheta);
    cosPhi = cos(dphi);
    cosPsi = cos(dpsi);
    sinTheta = sin(dtheta);
    sinPhi = sin(dphi);
    sinPsi = sin(dpsi);

    /*--- Compute the rotation matrix. Note that the implicit
     ordering is rotation about the x-axis, y-axis, then z-axis. ---*/

    rotMatrix[0][0] = cosPhi * cosPsi;
    rotMatrix[1][0] = cosPhi * sinPsi;
    rotMatrix[2][0] = -sinPhi;

    rotMatrix[0][1] = sinTheta * sinPhi * cosPsi - cosTheta * sinPsi;
    rotMatrix[1][1] = sinTheta * sinPhi * sinPsi + cosTheta * cosPsi;
    rotMatrix[2][1] = sinTheta * cosPhi;

    rotMatrix[0][2] = cosTheta * sinPhi * cosPsi + sinTheta * sinPsi;
    rotMatrix[1][2] = cosTheta * sinPhi * sinPsi - sinTheta * cosPsi;
    rotMatrix[2][2] = cosTheta * cosPhi;
  }

  /*--- Loop through to find only moving surface markers ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_DV(iMarker) == YES && config->GetKind_SU2() == SU2_COMPONENT::SU2_DEF) ||
        (config->GetMarker_All_Moving(iMarker) == YES && config->GetKind_SU2() == SU2_COMPONENT::SU2_CFD)) {
      /*--- Loop over all surface points for this marker ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Get current and new coordinates from file ---*/

        Coord_Old = geometry->nodes->GetCoord(iPoint);
        Coord_New = geometry->vertex[iMarker][iVertex]->GetVarCoord();

        /*--- If we're also rotating, multiply each point by the
         rotation matrix. It is assumed that the coordinates in
         Coord_Old have already been rotated using SetRigid_Rotation(). ---*/

        if (config->GetKind_GridMovement() == EXTERNAL_ROTATION) {
          /*--- Calculate non-dim. position from rotation center ---*/

          for (iDim = 0; iDim < nDim; iDim++) r[iDim] = (Coord_New[iDim] - Center[iDim]) / Lref;
          if (nDim == 2) r[nDim] = 0.0;

          /*--- Compute transformed point coordinates ---*/

          rotCoord[0] = rotMatrix[0][0] * r[0] + rotMatrix[0][1] * r[1] + rotMatrix[0][2] * r[2] + Center[0];

          rotCoord[1] = rotMatrix[1][0] * r[0] + rotMatrix[1][1] * r[1] + rotMatrix[1][2] * r[2] + Center[1];

          rotCoord[2] = rotMatrix[2][0] * r[0] + rotMatrix[2][1] * r[1] + rotMatrix[2][2] * r[2] + Center[2];

          /*--- Copy rotated coords back to original array for consistency ---*/
          for (iDim = 0; iDim < nDim; iDim++) Coord_New[iDim] = rotCoord[iDim];
        }

        /*--- Calculate delta change in the x, y, & z directions ---*/
        for (iDim = 0; iDim < nDim; iDim++) VarCoord[iDim] = (Coord_New[iDim] - Coord_Old[iDim]) / Lref;
        if (nDim == 2) VarCoord[nDim] = 0.0;

        /*--- Set position changes to be applied by the spring analogy ---*/
        geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
    }
  }
}

void CSurfaceMovement::SetNACA_4Digits(CGeometry* boundary, CConfig* config) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3], *Coord, *Normal, Ycurv, Yesp;

  if (config->GetnDV() != 1) {
    cout << "This kind of design variable is not prepared for multiple deformations.";
    cin.get();
  }

  su2double Ya = config->GetParamDV(0, 0) / 100.0; /*--- Maximum camber as a fraction of the chord
           (100 m is the first of the four digits) ---*/
  su2double Xa = config->GetParamDV(0, 1) / 10.0;  /*--- Location of maximum camber as a fraction of
            the chord (10 p is the second digit in the NACA xxxx description) ---*/
  su2double t = config->GetParamDV(0, 2) / 100.0;  /*--- Maximum thickness as a fraction of the
              chord (so 100 t gives the last two digits in
              the NACA 4-digit denomination) ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0;
      VarCoord[1] = 0.0;
      VarCoord[2] = 0.0;
      if (config->GetMarker_All_DV(iMarker) == YES) {
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
        Normal = boundary->vertex[iMarker][iVertex]->GetNormal();

        if (Coord[0] < Xa)
          Ycurv = (2.0 * Xa * Coord[0] - pow(Coord[0], 2.0)) * (Ya / pow(Xa, 2.0));
        else
          Ycurv = ((1.0 - 2.0 * Xa) + 2.0 * Xa * Coord[0] - pow(Coord[0], 2.0)) * (Ya / pow((1.0 - Xa), 2.0));

        Yesp = t * (1.4845 * sqrt(Coord[0]) - 0.6300 * Coord[0] - 1.7580 * pow(Coord[0], 2.0) +
                    1.4215 * pow(Coord[0], 3.0) - 0.518 * pow(Coord[0], 4.0));

        if (Normal[1] > 0) VarCoord[1] = (Ycurv + Yesp) - Coord[1];
        if (Normal[1] < 0) VarCoord[1] = (Ycurv - Yesp) - Coord[1];
      }
      boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
    }
}

void CSurfaceMovement::SetParabolic(CGeometry* boundary, CConfig* config) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3], *Coord, *Normal;

  if (config->GetnDV() != 1) {
    cout << "This kind of design variable is not prepared for multiple deformations.";
    cin.get();
  }

  su2double c = config->GetParamDV(0, 0);         /*--- Center of the parabola ---*/
  su2double t = config->GetParamDV(0, 1) / 100.0; /*--- Thickness of the parabola ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0;
      VarCoord[1] = 0.0;
      VarCoord[2] = 0.0;
      if (config->GetMarker_All_DV(iMarker) == YES) {
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
        Normal = boundary->vertex[iMarker][iVertex]->GetNormal();

        if (Normal[1] > 0) {
          VarCoord[1] = t * (Coord[0] * Coord[0] - Coord[0]) / (2.0 * (c * c - c)) - Coord[1];
        }
        if (Normal[1] < 0) {
          VarCoord[1] = t * (Coord[0] - Coord[0] * Coord[0]) / (2.0 * (c * c - c)) - Coord[1];
        }
      }
      boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
    }
}

void CSurfaceMovement::SetAirfoil(CGeometry* boundary, CConfig* config) {
  unsigned long iVertex, n_Airfoil = 0;
  unsigned short iMarker, nUpper, nLower, iUpper, iLower, iVar, iDim;
  su2double *VarCoord, *Coord, NewYCoord, NewXCoord, *Coord_i, *Coord_ip1, yp1, ypn,
      Airfoil_Coord[2] = {0.0, 0.0}, factor, coeff = 10000, Upper, Lower, Arch = 0.0, TotalArch = 0.0, x_i, x_ip1, y_i,
      y_ip1;
  passivedouble AirfoilScale;
  vector<su2double> Svalue, Xcoord, Ycoord, Xcoord2, Ycoord2, Xcoord_Aux, Ycoord_Aux;
  bool AddBegin = true, AddEnd = true;
  char AirfoilFile[256], AirfoilFormat[15], MeshOrientation[15], AirfoilClose[15];
  ifstream airfoil_file;
  string text_line;
  int ierr = 0;

  unsigned short nDim = boundary->GetnDim();

  VarCoord = new su2double[nDim];
  for (iDim = 0; iDim < nDim; iDim++) VarCoord[iDim] = 0.0;

  /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
   deforming meshes (MARKER_MOVING), while SU2_DEF will use it for deforming
   meshes after imposing design variable surface deformations (DV_MARKER). ---*/

  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();

  /*--- Read the coordinates. Two main formats:
   - Selig are in an x, y format starting from trailing edge, along the upper surface to the leading
   edge and back around the lower surface to trailing edge.
   - Lednicer are upper surface points leading edge to trailing edge and then lower surface leading
   edge to trailing edge.
   ---*/

  /*--- Open the restart file, throw an error if this fails. ---*/

  cout << "Enter the name of file with the airfoil information: ";
  ierr = scanf("%255s", AirfoilFile);
  if (ierr == 0) {
    SU2_MPI::Error("No input read!!", CURRENT_FUNCTION);
  }
  airfoil_file.open(AirfoilFile, ios::in);
  if (airfoil_file.fail()) {
    SU2_MPI::Error(string("There is no airfoil file ") + string(AirfoilFile), CURRENT_FUNCTION);
  }
  cout << "Enter the format of the airfoil (Selig or Lednicer): ";
  ierr = scanf("%14s", AirfoilFormat);
  if (ierr == 0) {
    SU2_MPI::Error("No input read!!", CURRENT_FUNCTION);
  }

  cout << "Thickness scaling (1.0 means no scaling)?: ";
  ierr = scanf("%lf", &AirfoilScale);
  if (ierr == 0) {
    SU2_MPI::Error("No input read!!", CURRENT_FUNCTION);
  }

  cout << "Close the airfoil (Yes or No)?: ";
  ierr = scanf("%14s", AirfoilClose);
  if (ierr == 0) {
    SU2_MPI::Error("No input read!!", CURRENT_FUNCTION);
  }

  cout << "Surface mesh orientation (clockwise, or anticlockwise): ";
  ierr = scanf("%14s", MeshOrientation);
  if (ierr == 0) {
    SU2_MPI::Error("No input read!!", CURRENT_FUNCTION);
  }

  /*--- The first line is the header ---*/

  getline(airfoil_file, text_line);
  cout << "File info: " << text_line << endl;

  if (strcmp(AirfoilFormat, "Selig") == 0) {
    while (getline(airfoil_file, text_line)) {
      istringstream point_line(text_line);

      /*--- Read the x & y coordinates from this line of the file (anticlockwise) ---*/

      point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];

      /*--- Close the arifoil ---*/

      if (strcmp(AirfoilClose, "Yes") == 0)
        factor = -atan(coeff * (Airfoil_Coord[0] - 1.0)) * 2.0 / PI_NUMBER;
      else
        factor = 1.0;

      /*--- Store the coordinates in vectors ---*/

      Xcoord.push_back(Airfoil_Coord[0]);
      Ycoord.emplace_back(Airfoil_Coord[1] * factor * AirfoilScale);
    }
  }
  if (strcmp(AirfoilFormat, "Lednicer") == 0) {
    /*--- The second line is the number of points ---*/

    getline(airfoil_file, text_line);
    istringstream point_line(text_line);
    point_line >> Upper >> Lower;

    nUpper = SU2_TYPE::Int(Upper);
    nLower = SU2_TYPE::Int(Lower);

    Xcoord.resize(nUpper + nLower - 1);
    Ycoord.resize(nUpper + nLower - 1);

    /*--- White line ---*/

    getline(airfoil_file, text_line);

    for (iUpper = 0; iUpper < nUpper; iUpper++) {
      getline(airfoil_file, text_line);
      istringstream point_line(text_line);
      point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];
      Xcoord[nUpper - iUpper - 1] = Airfoil_Coord[0];

      if (strcmp(AirfoilClose, "Yes") == 0)
        factor = -atan(coeff * (Airfoil_Coord[0] - 1.0)) * 2.0 / PI_NUMBER;
      else
        factor = 1.0;

      Ycoord[nUpper - iUpper - 1] = Airfoil_Coord[1] * AirfoilScale * factor;
    }

    getline(airfoil_file, text_line);

    for (iLower = 0; iLower < nLower; iLower++) {
      getline(airfoil_file, text_line);
      istringstream point_line(text_line);
      point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];

      if (strcmp(AirfoilClose, "Yes") == 0)
        factor = -atan(coeff * (Airfoil_Coord[0] - 1.0)) * 2.0 / PI_NUMBER;
      else
        factor = 1.0;

      Xcoord[nUpper + iLower - 1] = Airfoil_Coord[0];
      Ycoord[nUpper + iLower - 1] = Airfoil_Coord[1] * AirfoilScale * factor;
    }
  }

  /*--- Check the coordinate (1,0) at the beginning and end of the file ---*/

  if (Xcoord[0] == 1.0) AddBegin = false;
  if (Xcoord[Xcoord.size() - 1] == 1.0) AddEnd = false;

  if (AddBegin) {
    Xcoord.insert(Xcoord.begin(), 1.0);
    Ycoord.insert(Ycoord.begin(), 0.0);
  }
  if (AddEnd) {
    Xcoord.emplace_back(1.0);
    Ycoord.emplace_back(0.0);
  }

  /*--- Change the orientation (depend on the input file, and the mesh file) ---*/

  if (strcmp(MeshOrientation, "clockwise") == 0) {
    for (iVar = 0; iVar < Xcoord.size(); iVar++) {
      Xcoord_Aux.push_back(Xcoord[iVar]);
      Ycoord_Aux.push_back(Ycoord[iVar]);
    }

    for (iVar = 0; iVar < Xcoord.size(); iVar++) {
      Xcoord[iVar] = Xcoord_Aux[Xcoord.size() - iVar - 1];
      Ycoord[iVar] = Ycoord_Aux[Xcoord.size() - iVar - 1];
    }
  }

  /*--- Compute the total arch length ---*/

  Arch = 0.0;
  Svalue.push_back(Arch);

  for (iVar = 0; iVar < Xcoord.size() - 1; iVar++) {
    x_i = Xcoord[iVar];
    x_ip1 = Xcoord[iVar + 1];
    y_i = Ycoord[iVar];
    y_ip1 = Ycoord[iVar + 1];
    Arch += sqrt((x_ip1 - x_i) * (x_ip1 - x_i) + (y_ip1 - y_i) * (y_ip1 - y_i));
    Svalue.push_back(Arch);
  }
  x_i = Xcoord[Xcoord.size() - 1];
  x_ip1 = Xcoord[0];
  y_i = Ycoord[Xcoord.size() - 1];
  y_ip1 = Ycoord[0];
  Arch += sqrt((x_ip1 - x_i) * (x_ip1 - x_i) + (y_ip1 - y_i) * (y_ip1 - y_i));

  /*--- Non dimensionalization ---*/

  for (iVar = 0; iVar < Svalue.size(); iVar++) {
    Svalue[iVar] /= Arch;
  }

  /*--- Close the restart file ---*/

  airfoil_file.close();

  /*--- Create a spline for X and Y coordiantes using the arch length ---*/

  n_Airfoil = Svalue.size();
  yp1 = (Xcoord[1] - Xcoord[0]) / (Svalue[1] - Svalue[0]);
  ypn = (Xcoord[n_Airfoil - 1] - Xcoord[n_Airfoil - 2]) / (Svalue[n_Airfoil - 1] - Svalue[n_Airfoil - 2]);

  CCubicSpline splineX(Svalue, Xcoord, CCubicSpline::FIRST, yp1, CCubicSpline::FIRST, ypn);

  n_Airfoil = Svalue.size();
  yp1 = (Ycoord[1] - Ycoord[0]) / (Svalue[1] - Svalue[0]);
  ypn = (Ycoord[n_Airfoil - 1] - Ycoord[n_Airfoil - 2]) / (Svalue[n_Airfoil - 1] - Svalue[n_Airfoil - 2]);

  CCubicSpline splineY(Svalue, Ycoord, CCubicSpline::FIRST, yp1, CCubicSpline::FIRST, ypn);

  TotalArch = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) ||
        ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DEF))) {
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker] - 1; iVertex++) {
        Coord_i = boundary->vertex[iMarker][iVertex]->GetCoord();
        Coord_ip1 = boundary->vertex[iMarker][iVertex + 1]->GetCoord();

        x_i = Coord_i[0];
        x_ip1 = Coord_ip1[0];
        y_i = Coord_i[1];
        y_ip1 = Coord_ip1[1];

        TotalArch += sqrt((x_ip1 - x_i) * (x_ip1 - x_i) + (y_ip1 - y_i) * (y_ip1 - y_i));
      }
      Coord_i = boundary->vertex[iMarker][boundary->nVertex[iMarker] - 1]->GetCoord();
      Coord_ip1 = boundary->vertex[iMarker][0]->GetCoord();
      x_i = Coord_i[0];
      x_ip1 = Coord_ip1[0];
      y_i = Coord_i[1];
      y_ip1 = Coord_ip1[1];
      TotalArch += sqrt((x_ip1 - x_i) * (x_ip1 - x_i) + (y_ip1 - y_i) * (y_ip1 - y_i));
    }
  }

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Arch = 0.0;
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0;
      VarCoord[1] = 0.0;
      VarCoord[2] = 0.0;
      if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DEF))) {
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();

        if (iVertex == 0)
          Arch = 0.0;
        else {
          Coord_i = boundary->vertex[iMarker][iVertex - 1]->GetCoord();
          Coord_ip1 = boundary->vertex[iMarker][iVertex]->GetCoord();
          x_i = Coord_i[0];
          x_ip1 = Coord_ip1[0];
          y_i = Coord_i[1];
          y_ip1 = Coord_ip1[1];
          Arch += sqrt((x_ip1 - x_i) * (x_ip1 - x_i) + (y_ip1 - y_i) * (y_ip1 - y_i)) / TotalArch;
        }

        NewXCoord = splineX(Arch);
        NewYCoord = splineY(Arch);

        /*--- Store the delta change in the x & y coordinates ---*/

        VarCoord[0] = NewXCoord - Coord[0];
        VarCoord[1] = NewYCoord - Coord[1];
      }

      boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
    }
  }

  delete[] VarCoord;
}

void CSurfaceMovement::ReadFFDInfo(CGeometry* geometry, CConfig* config, CFreeFormDefBox** FFDBox,
                                   const string& val_mesh_filename) {
  string text_line, iTag;
  ifstream mesh_file;
  su2double CPcoord[3], coord[] = {0, 0, 0};
  unsigned short degree[3], iFFDBox, iCornerPoints, iControlPoints, iMarker, iDegree, jDegree, kDegree, iChar,
      LevelFFDBox, nParentFFDBox, iParentFFDBox, nChildFFDBox, iChildFFDBox, nMarker, *nCornerPoints, *nControlPoints;
  unsigned long iSurfacePoints, iPoint, jPoint, iVertex, nVertex, nPoint, iElem = 0, nElem, my_nSurfPoints, nSurfPoints,
                                                                          *nSurfacePoints;
  su2double XCoord, YCoord;

  bool polar = (config->GetFFD_CoordSystem() == POLAR);
  unsigned short nDim = geometry->GetnDim(), iDim;
  unsigned short SplineOrder[3];
  unsigned short Blending = 0;

  mesh_file.open(val_mesh_filename);
  if (mesh_file.fail()) {
    SU2_MPI::Error("There is no geometry file (ReadFFDInfo)!!", CURRENT_FUNCTION);
  }

  while (getline(mesh_file, text_line)) {
    /*--- Read the inner elements ---*/

    string::size_type position = text_line.find("NELEM=", 0);
    if (position != string::npos) {
      text_line.erase(0, 6);
      nElem = atoi(text_line.c_str());
      for (iElem = 0; iElem < nElem; iElem++) {
        getline(mesh_file, text_line);
      }
    }

    /*--- Read the inner points ---*/

    position = text_line.find("NPOIN=", 0);
    if (position != string::npos) {
      text_line.erase(0, 6);
      nPoint = atoi(text_line.c_str());
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        getline(mesh_file, text_line);
      }
    }

    /*--- Read the boundaries  ---*/

    position = text_line.find("NMARK=", 0);
    if (position != string::npos) {
      text_line.erase(0, 6);
      nMarker = atoi(text_line.c_str());
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        getline(mesh_file, text_line);
        getline(mesh_file, text_line);
        text_line.erase(0, 13);
        nVertex = atoi(text_line.c_str());
        for (iVertex = 0; iVertex < nVertex; iVertex++) {
          getline(mesh_file, text_line);
        }
      }
    }

    /*--- Read the FFDBox information  ---*/

    position = text_line.find("FFD_NBOX=", 0);
    if (position != string::npos) {
      text_line.erase(0, 9);
      nFFDBox = atoi(text_line.c_str());

      if (rank == MASTER_NODE) cout << nFFDBox << " Free Form Deformation boxes." << endl;

      nCornerPoints = new unsigned short[nFFDBox];
      nControlPoints = new unsigned short[nFFDBox];
      nSurfacePoints = new unsigned long[nFFDBox];

      getline(mesh_file, text_line);
      text_line.erase(0, 11);
      nLevel = atoi(text_line.c_str());

      if (rank == MASTER_NODE) cout << nLevel << " Free Form Deformation nested levels." << endl;

      for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {
        /*--- Read the name of the FFD box ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 8);

        /*--- Remove extra data from the FFDBox name ---*/

        string::size_type position;
        for (iChar = 0; iChar < 20; iChar++) {
          position = text_line.find(' ', 0);
          if (position != string::npos) text_line.erase(position, 1);
          position = text_line.find('\r', 0);
          if (position != string::npos) text_line.erase(position, 1);
          position = text_line.find('\n', 0);
          if (position != string::npos) text_line.erase(position, 1);
        }

        string TagFFDBox = text_line;

        if (rank == MASTER_NODE) cout << "FFD box tag: " << TagFFDBox << ". ";

        /*--- Read the level of the FFD box ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 10);
        LevelFFDBox = atoi(text_line.c_str());

        if (rank == MASTER_NODE) cout << "FFD box level: " << LevelFFDBox << ". ";

        /*--- Read the degree of the FFD box ---*/

        if (nDim == 2) {
          if (polar) {
            getline(mesh_file, text_line);
            text_line.erase(0, 13);
            degree[0] = atoi(text_line.c_str());
            degree[1] = 1;
            getline(mesh_file, text_line);
            text_line.erase(0, 13);
            degree[2] = atoi(text_line.c_str());
          } else {
            getline(mesh_file, text_line);
            text_line.erase(0, 13);
            degree[0] = atoi(text_line.c_str());
            getline(mesh_file, text_line);
            text_line.erase(0, 13);
            degree[1] = atoi(text_line.c_str());
            degree[2] = 1;
          }
        } else {
          getline(mesh_file, text_line);
          text_line.erase(0, 13);
          degree[0] = atoi(text_line.c_str());
          getline(mesh_file, text_line);
          text_line.erase(0, 13);
          degree[1] = atoi(text_line.c_str());
          getline(mesh_file, text_line);
          text_line.erase(0, 13);
          degree[2] = atoi(text_line.c_str());
        }

        if (rank == MASTER_NODE) {
          if (nDim == 2) {
            if (polar)
              cout << "Degrees: " << degree[0] << ", " << degree[2] << "." << endl;
            else
              cout << "Degrees: " << degree[0] << ", " << degree[1] << "." << endl;
          } else
            cout << "Degrees: " << degree[0] << ", " << degree[1] << ", " << degree[2] << "." << endl;
        }

        getline(mesh_file, text_line);
        if (text_line.substr(0, 12) != "FFD_BLENDING") {
          SU2_MPI::Error(
              string("Deprecated FFD information found in mesh file.\n") +
                  string(
                      "FFD information generated with SU2 version <= 4.3 is incompatible with the current version.") +
                  string("Run SU2_DEF again with DV_KIND= FFD_SETTING."),
              CURRENT_FUNCTION);
        }
        text_line.erase(0, 14);
        if (text_line == "BEZIER") {
          Blending = BEZIER;
        }
        if (text_line == "BSPLINE_UNIFORM") {
          Blending = BSPLINE_UNIFORM;
        }

        if (Blending == BSPLINE_UNIFORM) {
          getline(mesh_file, text_line);
          text_line.erase(0, 17);
          SplineOrder[0] = atoi(text_line.c_str());
          getline(mesh_file, text_line);
          text_line.erase(0, 17);
          SplineOrder[1] = atoi(text_line.c_str());
          if (nDim == 3) {
            getline(mesh_file, text_line);
            text_line.erase(0, 17);
            SplineOrder[2] = atoi(text_line.c_str());
          } else {
            SplineOrder[2] = 2;
          }
        }
        if (rank == MASTER_NODE) {
          if (Blending == BSPLINE_UNIFORM) {
            cout << "FFD Blending using B-Splines. ";
            cout << "Order: " << SplineOrder[0] << ", " << SplineOrder[1];
            if (nDim == 3) cout << ", " << SplineOrder[2];
            cout << ". " << endl;
          }
          if (Blending == BEZIER) {
            cout << "FFD Blending using Bezier Curves." << endl;
          }
        }

        FFDBox[iFFDBox] = new CFreeFormDefBox(degree, SplineOrder, Blending);
        FFDBox[iFFDBox]->SetTag(TagFFDBox);
        FFDBox[iFFDBox]->SetLevel(LevelFFDBox);

        /*--- Read the number of parents boxes ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 12);
        nParentFFDBox = atoi(text_line.c_str());
        if (rank == MASTER_NODE) cout << "Number of parent boxes: " << nParentFFDBox << ". ";
        for (iParentFFDBox = 0; iParentFFDBox < nParentFFDBox; iParentFFDBox++) {
          getline(mesh_file, text_line);

          /*--- Remove extra data from the FFDBox name ---*/

          string::size_type position;
          for (iChar = 0; iChar < 20; iChar++) {
            position = text_line.find(' ', 0);
            if (position != string::npos) text_line.erase(position, 1);
            position = text_line.find('\r', 0);
            if (position != string::npos) text_line.erase(position, 1);
            position = text_line.find('\n', 0);
            if (position != string::npos) text_line.erase(position, 1);
          }

          string ParentFFDBox = text_line;
          FFDBox[iFFDBox]->SetParentFFDBox(ParentFFDBox);
        }

        /*--- Read the number of children boxes ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 13);
        nChildFFDBox = atoi(text_line.c_str());
        if (rank == MASTER_NODE) cout << "Number of child boxes: " << nChildFFDBox << "." << endl;

        for (iChildFFDBox = 0; iChildFFDBox < nChildFFDBox; iChildFFDBox++) {
          getline(mesh_file, text_line);

          /*--- Remove extra data from the FFDBox name ---*/

          string::size_type position;
          for (iChar = 0; iChar < 20; iChar++) {
            position = text_line.find(' ', 0);
            if (position != string::npos) text_line.erase(position, 1);
            position = text_line.find('\r', 0);
            if (position != string::npos) text_line.erase(position, 1);
            position = text_line.find('\n', 0);
            if (position != string::npos) text_line.erase(position, 1);
          }

          string ChildFFDBox = text_line;
          FFDBox[iFFDBox]->SetChildFFDBox(ChildFFDBox);
        }

        /*--- Read the number of the corner points ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 18);
        nCornerPoints[iFFDBox] = atoi(text_line.c_str());
        if (rank == MASTER_NODE) cout << "Corner points: " << nCornerPoints[iFFDBox] << ". ";
        if (nDim == 2) nCornerPoints[iFFDBox] = nCornerPoints[iFFDBox] * SU2_TYPE::Int(2);

        /*--- Read the coordinates of the corner points ---*/

        if (nDim == 2) {
          if (polar) {
            getline(mesh_file, text_line);
            istringstream FFDBox_line_1(text_line);
            FFDBox_line_1 >> XCoord;
            FFDBox_line_1 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = -sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 4);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 7);

            getline(mesh_file, text_line);
            istringstream FFDBox_line_2(text_line);
            FFDBox_line_2 >> XCoord;
            FFDBox_line_2 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = -sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 0);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 3);

            getline(mesh_file, text_line);
            istringstream FFDBox_line_3(text_line);
            FFDBox_line_3 >> XCoord;
            FFDBox_line_3 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = -sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 1);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 2);

            getline(mesh_file, text_line);
            istringstream FFDBox_line_4(text_line);
            FFDBox_line_4 >> XCoord;
            FFDBox_line_4 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = -sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 5);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 6);

          } else {
            for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iFFDBox]; iCornerPoints++) {
              if (iCornerPoints < nCornerPoints[iFFDBox] / SU2_TYPE::Int(2)) {
                getline(mesh_file, text_line);
                istringstream FFDBox_line(text_line);
                FFDBox_line >> CPcoord[0];
                FFDBox_line >> CPcoord[1];
                CPcoord[2] = -0.5;
              } else {
                CPcoord[0] =
                    FFDBox[iFFDBox]->GetCoordCornerPoints(0, iCornerPoints - nCornerPoints[iFFDBox] / SU2_TYPE::Int(2));
                CPcoord[1] =
                    FFDBox[iFFDBox]->GetCoordCornerPoints(1, iCornerPoints - nCornerPoints[iFFDBox] / SU2_TYPE::Int(2));
                CPcoord[2] = 0.5;
              }
              FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, iCornerPoints);
            }
          }

        } else {
          for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iFFDBox]; iCornerPoints++) {
            getline(mesh_file, text_line);
            istringstream FFDBox_line(text_line);
            FFDBox_line >> CPcoord[0];
            FFDBox_line >> CPcoord[1];
            FFDBox_line >> CPcoord[2];
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, iCornerPoints);
          }
        }

        /*--- Read the number of the control points ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 19);
        nControlPoints[iFFDBox] = atoi(text_line.c_str());

        if (rank == MASTER_NODE) cout << "Control points: " << nControlPoints[iFFDBox] << ". ";

        /*--- Method to identify if there is a FFDBox definition ---*/

        if (nControlPoints[iFFDBox] != 0) FFDBoxDefinition = true;

        /*--- Read the coordinates of the control points ---*/

        for (iControlPoints = 0; iControlPoints < nControlPoints[iFFDBox]; iControlPoints++) {
          getline(mesh_file, text_line);
          istringstream FFDBox_line(text_line);
          FFDBox_line >> iDegree;
          FFDBox_line >> jDegree;
          FFDBox_line >> kDegree;
          FFDBox_line >> CPcoord[0];
          FFDBox_line >> CPcoord[1];
          FFDBox_line >> CPcoord[2];
          FFDBox[iFFDBox]->SetCoordControlPoints(CPcoord, iDegree, jDegree, kDegree);
          FFDBox[iFFDBox]->SetCoordControlPoints_Copy(CPcoord, iDegree, jDegree, kDegree);
        }

        getline(mesh_file, text_line);
        text_line.erase(0, 19);
        nSurfacePoints[iFFDBox] = atoi(text_line.c_str());

        /*--- The surface points parametric coordinates, all the nodes read the FFD
         information but they only store their part ---*/

        my_nSurfPoints = 0;
        for (iSurfacePoints = 0; iSurfacePoints < nSurfacePoints[iFFDBox]; iSurfacePoints++) {
          getline(mesh_file, text_line);
          istringstream FFDBox_line(text_line);
          FFDBox_line >> iTag;
          FFDBox_line >> iPoint;

          if (config->GetMarker_All_TagBound(iTag) != -1) {
            iMarker = config->GetMarker_All_TagBound(iTag);
            FFDBox_line >> CPcoord[0];
            FFDBox_line >> CPcoord[1];
            FFDBox_line >> CPcoord[2];

            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
              jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              if (iPoint == geometry->nodes->GetGlobalIndex(jPoint)) {
                for (iDim = 0; iDim < nDim; iDim++) {
                  coord[iDim] = geometry->nodes->GetCoord(jPoint, iDim);
                }
                FFDBox[iFFDBox]->Set_MarkerIndex(iMarker);
                FFDBox[iFFDBox]->Set_VertexIndex(iVertex);
                FFDBox[iFFDBox]->Set_PointIndex(jPoint);
                FFDBox[iFFDBox]->Set_ParametricCoord(CPcoord);
                FFDBox[iFFDBox]->Set_CartesianCoord(coord);
                my_nSurfPoints++;
              }
            }
          }
        }

        nSurfacePoints[iFFDBox] = my_nSurfPoints;

#ifdef HAVE_MPI
        nSurfPoints = 0;
        SU2_MPI::Allreduce(&my_nSurfPoints, &nSurfPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
        if (rank == MASTER_NODE) cout << "Surface points: " << nSurfPoints << "." << endl;
#else
        nSurfPoints = my_nSurfPoints;
        if (rank == MASTER_NODE) cout << "Surface points: " << nSurfPoints << "." << endl;
#endif
      }

      delete[] nCornerPoints;
      delete[] nControlPoints;
      delete[] nSurfacePoints;
    }
  }
  mesh_file.close();

  if (nFFDBox == 0) {
    if (rank == MASTER_NODE) cout << "There is no FFD box definition. Just in case, check the .su2 file" << endl;
  }
}

void CSurfaceMovement::ReadFFDInfo(CGeometry* geometry, CConfig* config, CFreeFormDefBox** FFDBox) {
  string text_line, iTag;
  ifstream mesh_file;
  su2double coord[3];
  unsigned short degree[3], iFFDBox, iCornerPoints, LevelFFDBox, nParentFFDBox, iParentFFDBox, nChildFFDBox,
      iChildFFDBox, *nCornerPoints;

  bool polar = (config->GetFFD_CoordSystem() == POLAR);
  unsigned short nDim = geometry->GetnDim(), iDim;
  unsigned short SplineOrder[3] = {2, 2, 2};

  for (iDim = 0; iDim < 3; iDim++) {
    SplineOrder[iDim] = SU2_TYPE::Short(config->GetFFD_BSplineOrder()[iDim]);
  }

  /*--- Read the FFDBox information from the config file ---*/

  nFFDBox = config->GetnFFDBox();

  if (rank == MASTER_NODE) cout << nFFDBox << " Free Form Deformation boxes." << endl;

  nCornerPoints = new unsigned short[nFFDBox];

  nLevel = 1;  // Nested FFD is not active

  if (rank == MASTER_NODE) cout << nLevel << " Free Form Deformation nested levels." << endl;

  for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {
    /*--- Read the name of the FFD box ---*/

    string TagFFDBox = config->GetTagFFDBox(iFFDBox);

    if (rank == MASTER_NODE) cout << "FFD box tag: " << TagFFDBox << ". ";

    /*--- Read the level of the FFD box ---*/

    LevelFFDBox = 0;  // Nested FFD is not active

    if (rank == MASTER_NODE) cout << "FFD box level: " << LevelFFDBox << ". ";

    /*--- Read the degree of the FFD box ---*/

    if (nDim == 2) {
      if (polar) {
        degree[0] = config->GetDegreeFFDBox(iFFDBox, 0);
        degree[1] = 1;
        degree[2] = config->GetDegreeFFDBox(iFFDBox, 1);
      } else {
        degree[0] = config->GetDegreeFFDBox(iFFDBox, 0);
        degree[1] = config->GetDegreeFFDBox(iFFDBox, 1);
        degree[2] = 1;
      }
    } else {
      degree[0] = config->GetDegreeFFDBox(iFFDBox, 0);
      degree[1] = config->GetDegreeFFDBox(iFFDBox, 1);
      degree[2] = config->GetDegreeFFDBox(iFFDBox, 2);
    }

    if (rank == MASTER_NODE) {
      if (nDim == 2) {
        if (polar)
          cout << "Degrees: " << degree[0] << ", " << degree[2] << "." << endl;
        else
          cout << "Degrees: " << degree[0] << ", " << degree[1] << "." << endl;
      } else
        cout << "Degrees: " << degree[0] << ", " << degree[1] << ", " << degree[2] << "." << endl;
    }

    if (rank == MASTER_NODE) {
      if (config->GetFFD_Blending() == BSPLINE_UNIFORM) {
        cout << "FFD Blending using B-Splines. ";
        cout << "Order: " << SplineOrder[0] << ", " << SplineOrder[1];
        if (nDim == 3) cout << ", " << SplineOrder[2];
        cout << ". " << endl;
      }
      if (config->GetFFD_Blending() == BEZIER) {
        cout << "FFD Blending using Bezier Curves." << endl;
      }
    }

    FFDBox[iFFDBox] = new CFreeFormDefBox(degree, SplineOrder, config->GetFFD_Blending());
    FFDBox[iFFDBox]->SetTag(TagFFDBox);
    FFDBox[iFFDBox]->SetLevel(LevelFFDBox);

    /*--- Read the number of parents boxes ---*/

    nParentFFDBox = 0;  // Nested FFD is not active
    if (rank == MASTER_NODE) cout << "Number of parent boxes: " << nParentFFDBox << ". ";

    for (iParentFFDBox = 0; iParentFFDBox < nParentFFDBox; iParentFFDBox++) {
      string ParentFFDBox = "NONE";  // Nested FFD is not active
      FFDBox[iFFDBox]->SetParentFFDBox(ParentFFDBox);
    }

    /*--- Read the number of children boxes ---*/

    nChildFFDBox = 0;  // Nested FFD is not active
    if (rank == MASTER_NODE) cout << "Number of child boxes: " << nChildFFDBox << "." << endl;

    for (iChildFFDBox = 0; iChildFFDBox < nChildFFDBox; iChildFFDBox++) {
      string ChildFFDBox = "NONE";  // Nested FFD is not active
      FFDBox[iFFDBox]->SetChildFFDBox(ChildFFDBox);
    }

    /*--- Read the number of the corner points ---*/

    nCornerPoints[iFFDBox] = 8;

    /*--- Read the coordinates of the corner points ---*/

    for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iFFDBox]; iCornerPoints++) {
      if (nDim == 2) {
        if (polar) {
          coord[0] = config->GetCoordFFDBox(iFFDBox, 1 * 3);
          coord[1] = cos(0.1) * config->GetCoordFFDBox(iFFDBox, 1 * 3 + 1);
          coord[2] = -sin(0.1) * config->GetCoordFFDBox(iFFDBox, 1 * 3 + 1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 0);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 2 * 3);
          coord[1] = cos(0.1) * config->GetCoordFFDBox(iFFDBox, 2 * 3 + 1);
          coord[2] = -sin(0.1) * config->GetCoordFFDBox(iFFDBox, 2 * 3 + 1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 1);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 2 * 3);
          coord[1] = cos(0.1) * config->GetCoordFFDBox(iFFDBox, 2 * 3 + 1);
          coord[2] = sin(0.1) * config->GetCoordFFDBox(iFFDBox, 2 * 3 + 1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 2);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 1 * 3);
          coord[1] = cos(0.1) * config->GetCoordFFDBox(iFFDBox, 1 * 3 + 1);
          coord[2] = sin(0.1) * config->GetCoordFFDBox(iFFDBox, 1 * 3 + 1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 3);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 0 * 3);
          coord[1] = cos(0.1) * config->GetCoordFFDBox(iFFDBox, 0 * 3 + 1);
          coord[2] = -sin(0.1) * config->GetCoordFFDBox(iFFDBox, 0 * 3 + 1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 4);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 3 * 3);
          coord[1] = cos(0.1) * config->GetCoordFFDBox(iFFDBox, 3 * 3 + 1);
          coord[2] = -sin(0.1) * config->GetCoordFFDBox(iFFDBox, 3 * 3 + 1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 5);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 3 * 3);
          coord[1] = cos(0.1) * config->GetCoordFFDBox(iFFDBox, 3 * 3 + 1);
          coord[2] = sin(0.1) * config->GetCoordFFDBox(iFFDBox, 3 * 3 + 1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 6);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 0 * 3);
          coord[1] = cos(0.1) * config->GetCoordFFDBox(iFFDBox, 0 * 3 + 1);
          coord[2] = sin(0.1) * config->GetCoordFFDBox(iFFDBox, 0 * 3 + 1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 7);

        }

        else {
          if (iCornerPoints < nCornerPoints[iFFDBox] / SU2_TYPE::Int(2)) {
            coord[0] = config->GetCoordFFDBox(iFFDBox, iCornerPoints * 3);
            coord[1] = config->GetCoordFFDBox(iFFDBox, iCornerPoints * 3 + 1);
            coord[2] = -0.5;
          } else {
            coord[0] =
                FFDBox[iFFDBox]->GetCoordCornerPoints(0, iCornerPoints - nCornerPoints[iFFDBox] / SU2_TYPE::Int(2));
            coord[1] =
                FFDBox[iFFDBox]->GetCoordCornerPoints(1, iCornerPoints - nCornerPoints[iFFDBox] / SU2_TYPE::Int(2));
            coord[2] = 0.5;
          }
        }

      } else {
        coord[0] = config->GetCoordFFDBox(iFFDBox, iCornerPoints * 3);
        coord[1] = config->GetCoordFFDBox(iFFDBox, iCornerPoints * 3 + 1);
        coord[2] = config->GetCoordFFDBox(iFFDBox, iCornerPoints * 3 + 2);
      }

      FFDBox[iFFDBox]->SetCoordCornerPoints(coord, iCornerPoints);
    }

    /*--- Method to identify if there is a FFDBox definition ---*/

    FFDBoxDefinition = false;
  }

  delete[] nCornerPoints;

  if (nFFDBox == 0) {
    SU2_MPI::Error("There is no FFD box definition. Check the config file.", CURRENT_FUNCTION);
  }
}

void CSurfaceMovement::MergeFFDInfo(CGeometry* geometry, CConfig* config) {
  /*--- Local variables needed on all processors ---*/

  unsigned long iPoint;
  unsigned short iFFDBox;

#ifndef HAVE_MPI

  /*--- In serial, the single process has access to all geometry, so simply
   load the coordinates into the data structure. ---*/

  /*--- Total number of points in each FFD box. ---*/

  for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {
    /*--- Loop over the mesh to collect the coords of the local points. ---*/

    for (iPoint = 0; iPoint < FFDBox[iFFDBox]->GetnSurfacePoint(); iPoint++) {
      /*--- Retrieve the current parametric coordinates at this node. ---*/

      GlobalCoordX[iFFDBox].push_back(FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[0]);
      GlobalCoordY[iFFDBox].push_back(FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[1]);
      GlobalCoordZ[iFFDBox].push_back(FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[2]);
      GlobalPoint[iFFDBox].push_back(FFDBox[iFFDBox]->Get_PointIndex(iPoint));

      /*--- Marker of the boundary in the local domain. ---*/

      unsigned short MarkerIndex = FFDBox[iFFDBox]->Get_MarkerIndex(iPoint);
      string TagBound = config->GetMarker_All_TagBound(MarkerIndex);

      /*--- Find the Marker of the boundary in the config file. ---*/

      unsigned short MarkerIndex_CfgFile = config->GetMarker_CfgFile_TagBound(TagBound);
      string TagBound_CfgFile = config->GetMarker_CfgFile_TagBound(MarkerIndex_CfgFile);

      /*--- Set the value of the tag at this node. ---*/

      GlobalTag[iFFDBox].push_back(TagBound_CfgFile);
    }
  }

#else

  /*--- MPI preprocessing ---*/

  int iProcessor, nProcessor = size;

  /*--- Local variables needed for merging the geometry with MPI. ---*/

  unsigned long jPoint, iPointLocal;
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = nullptr;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long nBuffer_Scalar = 0;

  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];

  for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {
    nLocalPoint = 0;
    for (iPoint = 0; iPoint < FFDBox[iFFDBox]->GetnSurfacePoint(); iPoint++) {
      iPointLocal = FFDBox[iFFDBox]->Get_PointIndex(iPoint);

      if (iPointLocal < geometry->GetnPointDomain()) {
        nLocalPoint++;
      }
    }
    Buffer_Send_nPoint[0] = nLocalPoint;

    /*--- Communicate the total number of nodes on this domain. ---*/

    SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE,
                    SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());

    nBuffer_Scalar = MaxLocalPoint;

    /*--- Send and Recv buffers. ---*/

    auto* Buffer_Send_X = new su2double[MaxLocalPoint];
    su2double* Buffer_Recv_X = nullptr;

    auto* Buffer_Send_Y = new su2double[MaxLocalPoint];
    su2double* Buffer_Recv_Y = nullptr;

    auto* Buffer_Send_Z = new su2double[MaxLocalPoint];
    su2double* Buffer_Recv_Z = nullptr;

    auto* Buffer_Send_Point = new unsigned long[MaxLocalPoint];
    unsigned long* Buffer_Recv_Point = nullptr;

    auto* Buffer_Send_MarkerIndex_CfgFile = new unsigned short[MaxLocalPoint];
    unsigned short* Buffer_Recv_MarkerIndex_CfgFile = nullptr;

    /*--- Prepare the receive buffers in the master node only. ---*/

    if (rank == MASTER_NODE) {
      Buffer_Recv_X = new su2double[nProcessor * MaxLocalPoint];
      Buffer_Recv_Y = new su2double[nProcessor * MaxLocalPoint];
      Buffer_Recv_Z = new su2double[nProcessor * MaxLocalPoint];
      Buffer_Recv_Point = new unsigned long[nProcessor * MaxLocalPoint];
      Buffer_Recv_MarkerIndex_CfgFile = new unsigned short[nProcessor * MaxLocalPoint];
    }

    /*--- Main communication routine. Loop over each coordinate and perform
     the MPI comm. Temporary 1-D buffers are used to send the coordinates at
     all nodes on each partition to the master node. These are then unpacked
     by the master and sorted by global index in one large n-dim. array. ---*/

    /*--- Loop over this partition to collect the coords of the local points. ---*/

    jPoint = 0;
    for (iPoint = 0; iPoint < FFDBox[iFFDBox]->GetnSurfacePoint(); iPoint++) {
      iPointLocal = FFDBox[iFFDBox]->Get_PointIndex(iPoint);

      if (iPointLocal < geometry->GetnPointDomain()) {
        /*--- Load local coords into the temporary send buffer. ---*/

        Buffer_Send_X[jPoint] = FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[0];
        Buffer_Send_Y[jPoint] = FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[1];
        Buffer_Send_Z[jPoint] = FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[2];

        /*--- Store the global index for this local node. ---*/

        Buffer_Send_Point[jPoint] = geometry->nodes->GetGlobalIndex(FFDBox[iFFDBox]->Get_PointIndex(iPoint));

        /*--- Marker of the boundary in the local domain. ---*/

        unsigned short MarkerIndex = FFDBox[iFFDBox]->Get_MarkerIndex(iPoint);
        string TagBound = config->GetMarker_All_TagBound(MarkerIndex);

        /*--- Find the Marker of the boundary in the config file.---*/

        unsigned short MarkerIndex_CfgFile = config->GetMarker_CfgFile_TagBound(TagBound);
        Buffer_Send_MarkerIndex_CfgFile[jPoint] = MarkerIndex_CfgFile;

        jPoint++;
      }
    }

    /*--- Gather the coordinate data on the master node using MPI. ---*/

    SU2_MPI::Gather(Buffer_Send_X, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_X, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE,
                    SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_Y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE,
                    SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_Z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE,
                    SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_Point, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Point, nBuffer_Scalar,
                    MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_MarkerIndex_CfgFile, nBuffer_Scalar, MPI_UNSIGNED_SHORT,
                    Buffer_Recv_MarkerIndex_CfgFile, nBuffer_Scalar, MPI_UNSIGNED_SHORT, MASTER_NODE,
                    SU2_MPI::GetComm());

    /*--- The master node unpacks and sorts this variable by global index ---*/

    if (rank == MASTER_NODE) {
      jPoint = 0;

      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          /*--- Get global index, then loop over each variable and store ---*/

          GlobalCoordX[iFFDBox].push_back(Buffer_Recv_X[jPoint]);
          GlobalCoordY[iFFDBox].push_back(Buffer_Recv_Y[jPoint]);
          GlobalCoordZ[iFFDBox].push_back(Buffer_Recv_Z[jPoint]);
          GlobalPoint[iFFDBox].push_back(Buffer_Recv_Point[jPoint]);

          string TagBound_CfgFile = config->GetMarker_CfgFile_TagBound(Buffer_Recv_MarkerIndex_CfgFile[jPoint]);
          GlobalTag[iFFDBox].push_back(TagBound_CfgFile);
          jPoint++;
        }

        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

        jPoint = (iProcessor + 1) * nBuffer_Scalar;
      }
    }

    /*--- Immediately release the temporary data buffers. ---*/

    delete[] Buffer_Send_X;
    delete[] Buffer_Send_Y;
    delete[] Buffer_Send_Z;
    delete[] Buffer_Send_Point;
    delete[] Buffer_Send_MarkerIndex_CfgFile;

    if (rank == MASTER_NODE) {
      delete[] Buffer_Recv_X;
      delete[] Buffer_Recv_Y;
      delete[] Buffer_Recv_Z;
      delete[] Buffer_Recv_Point;
      delete[] Buffer_Recv_MarkerIndex_CfgFile;
    }
  }

  if (rank == MASTER_NODE) {
    delete[] Buffer_Recv_nPoint;
  }

#endif
}

void CSurfaceMovement::WriteFFDInfo(CSurfaceMovement** surface_movement, CGeometry**** geometry, CConfig** config) {
  unsigned short iOrder, jOrder, kOrder, iFFDBox, iCornerPoints, iParentFFDBox, iChildFFDBox, iZone;
  unsigned long iSurfacePoints;
  ofstream output_file;
  su2double* coord;
  string text_line;

  bool polar = (config[ZONE_0]->GetFFD_CoordSystem() == POLAR);

  unsigned short nDim = geometry[ZONE_0][INST_0][MESH_0]->GetnDim();

  for (iZone = 0; iZone < config[ZONE_0]->GetnZone(); iZone++) {
    /*--- Merge the parallel FFD info ---*/

    surface_movement[iZone]->MergeFFDInfo(geometry[iZone][INST_0][MESH_0], config[iZone]);

    if (iZone > 0) {
      /* --- Merge the per-zone FFD info from the other zones into ZONE_0 ---*/

      for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {
        surface_movement[ZONE_0]->GlobalCoordX[iFFDBox].insert(surface_movement[ZONE_0]->GlobalCoordX[iFFDBox].end(),
                                                               surface_movement[iZone]->GlobalCoordX[iFFDBox].begin(),
                                                               surface_movement[iZone]->GlobalCoordX[iFFDBox].end());
        surface_movement[ZONE_0]->GlobalCoordY[iFFDBox].insert(surface_movement[ZONE_0]->GlobalCoordY[iFFDBox].end(),
                                                               surface_movement[iZone]->GlobalCoordY[iFFDBox].begin(),
                                                               surface_movement[iZone]->GlobalCoordY[iFFDBox].end());
        surface_movement[ZONE_0]->GlobalCoordZ[iFFDBox].insert(surface_movement[ZONE_0]->GlobalCoordZ[iFFDBox].end(),
                                                               surface_movement[iZone]->GlobalCoordZ[iFFDBox].begin(),
                                                               surface_movement[iZone]->GlobalCoordZ[iFFDBox].end());
        surface_movement[ZONE_0]->GlobalTag[iFFDBox].insert(surface_movement[ZONE_0]->GlobalTag[iFFDBox].end(),
                                                            surface_movement[iZone]->GlobalTag[iFFDBox].begin(),
                                                            surface_movement[iZone]->GlobalTag[iFFDBox].end());
        surface_movement[ZONE_0]->GlobalPoint[iFFDBox].insert(surface_movement[ZONE_0]->GlobalPoint[iFFDBox].end(),
                                                              surface_movement[iZone]->GlobalPoint[iFFDBox].begin(),
                                                              surface_movement[iZone]->GlobalPoint[iFFDBox].end());
      }
    }
  }

  /*--- Attach to the mesh file the FFD information (all information is in ZONE_0) ---*/

  if (rank == MASTER_NODE) {
    /*--- Read the name of the output file ---*/

    auto str = config[ZONE_0]->GetMesh_Out_FileName();
    unsigned short lastindex = str.find_last_of('.');
    str = str.substr(0, lastindex) + ".su2";

    output_file.precision(15);
    output_file.open(str, ios::out | ios::app);

    if (nFFDBox != 0) {
      output_file << "FFD_NBOX= " << nFFDBox << endl;
      output_file << "FFD_NLEVEL= " << nLevel << endl;
    }

    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {
      output_file << "FFD_TAG= " << FFDBox[iFFDBox]->GetTag() << endl;
      output_file << "FFD_LEVEL= " << FFDBox[iFFDBox]->GetLevel() << endl;

      output_file << "FFD_DEGREE_I= " << FFDBox[iFFDBox]->GetlOrder() - 1 << endl;
      if (polar)
        output_file << "FFD_DEGREE_J= " << FFDBox[iFFDBox]->GetnOrder() - 1 << endl;
      else
        output_file << "FFD_DEGREE_J= " << FFDBox[iFFDBox]->GetmOrder() - 1 << endl;
      if (nDim == 3) output_file << "FFD_DEGREE_K= " << FFDBox[iFFDBox]->GetnOrder() - 1 << endl;
      if (config[ZONE_0]->GetFFD_Blending() == BSPLINE_UNIFORM) {
        output_file << "FFD_BLENDING= BSPLINE_UNIFORM" << endl;
        output_file << "BSPLINE_ORDER_I= " << FFDBox[iFFDBox]->BlendingFunction[0]->GetOrder() << endl;
        if (polar)
          output_file << "BSPLINE_ORDER_J= " << FFDBox[iFFDBox]->BlendingFunction[2]->GetOrder() << endl;
        else
          output_file << "BSPLINE_ORDER_J= " << FFDBox[iFFDBox]->BlendingFunction[1]->GetOrder() << endl;
        if (nDim == 3) output_file << "BSPLINE_ORDER_K= " << FFDBox[iFFDBox]->BlendingFunction[2]->GetOrder() << endl;
      }
      if (config[ZONE_0]->GetFFD_Blending() == BEZIER) {
        output_file << "FFD_BLENDING= BEZIER" << endl;
      }

      output_file << "FFD_PARENTS= " << FFDBox[iFFDBox]->GetnParentFFDBox() << endl;
      for (iParentFFDBox = 0; iParentFFDBox < FFDBox[iFFDBox]->GetnParentFFDBox(); iParentFFDBox++)
        output_file << FFDBox[iFFDBox]->GetParentFFDBoxTag(iParentFFDBox) << endl;
      output_file << "FFD_CHILDREN= " << FFDBox[iFFDBox]->GetnChildFFDBox() << endl;
      for (iChildFFDBox = 0; iChildFFDBox < FFDBox[iFFDBox]->GetnChildFFDBox(); iChildFFDBox++)
        output_file << FFDBox[iFFDBox]->GetChildFFDBoxTag(iChildFFDBox) << endl;

      if (nDim == 2) {
        output_file << "FFD_CORNER_POINTS= " << FFDBox[iFFDBox]->GetnCornerPoints() / SU2_TYPE::Int(2) << endl;
        if (polar) {
          coord = FFDBox[iFFDBox]->GetCoordCornerPoints(4);
          output_file << coord[0] << "\t" << sqrt(coord[1] * coord[1] + coord[2] * coord[2]) << endl;

          coord = FFDBox[iFFDBox]->GetCoordCornerPoints(0);
          output_file << coord[0] << "\t" << sqrt(coord[1] * coord[1] + coord[2] * coord[2]) << endl;

          coord = FFDBox[iFFDBox]->GetCoordCornerPoints(1);
          output_file << coord[0] << "\t" << sqrt(coord[1] * coord[1] + coord[2] * coord[2]) << endl;

          coord = FFDBox[iFFDBox]->GetCoordCornerPoints(5);
          output_file << coord[0] << "\t" << sqrt(coord[1] * coord[1] + coord[2] * coord[2]) << endl;
        } else {
          for (iCornerPoints = 0; iCornerPoints < FFDBox[iFFDBox]->GetnCornerPoints() / SU2_TYPE::Int(2);
               iCornerPoints++) {
            coord = FFDBox[iFFDBox]->GetCoordCornerPoints(iCornerPoints);
            output_file << coord[0] << "\t" << coord[1] << endl;
          }
        }
      } else {
        output_file << "FFD_CORNER_POINTS= " << FFDBox[iFFDBox]->GetnCornerPoints() << endl;
        for (iCornerPoints = 0; iCornerPoints < FFDBox[iFFDBox]->GetnCornerPoints(); iCornerPoints++) {
          coord = FFDBox[iFFDBox]->GetCoordCornerPoints(iCornerPoints);
          output_file << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
        }
      }

      /*--- Writing control points ---*/

      if (FFDBox[iFFDBox]->GetnControlPoints() == 0) {
        output_file << "FFD_CONTROL_POINTS= 0" << endl;
      } else {
        output_file << "FFD_CONTROL_POINTS= " << FFDBox[iFFDBox]->GetnControlPoints() << endl;
        for (iOrder = 0; iOrder < FFDBox[iFFDBox]->GetlOrder(); iOrder++)
          for (jOrder = 0; jOrder < FFDBox[iFFDBox]->GetmOrder(); jOrder++)
            for (kOrder = 0; kOrder < FFDBox[iFFDBox]->GetnOrder(); kOrder++) {
              coord = FFDBox[iFFDBox]->GetCoordControlPoints(iOrder, jOrder, kOrder);
              output_file << iOrder << "\t" << jOrder << "\t" << kOrder << "\t" << coord[0] << "\t" << coord[1] << "\t"
                          << coord[2] << endl;
            }
      }

      /*--- Writing surface points ---*/

      if (FFDBox[iFFDBox]->GetnControlPoints() == 0) {
        output_file << "FFD_SURFACE_POINTS= 0" << endl;
      } else {
        output_file << "FFD_SURFACE_POINTS= " << GlobalTag[iFFDBox].size() << endl;

        for (iSurfacePoints = 0; iSurfacePoints < GlobalTag[iFFDBox].size(); iSurfacePoints++) {
          output_file << scientific << GlobalTag[iFFDBox][iSurfacePoints] << "\t"
                      << GlobalPoint[iFFDBox][iSurfacePoints] << "\t" << GlobalCoordX[iFFDBox][iSurfacePoints] << "\t"
                      << GlobalCoordY[iFFDBox][iSurfacePoints] << "\t" << GlobalCoordZ[iFFDBox][iSurfacePoints] << endl;
        }
      }
    }

    output_file.close();
  }
}

unsigned long CSurfaceMovement::calculateJacobianDeterminant(CGeometry* geometry, CConfig* config,
                                                             CFreeFormDefBox* FFDBox) const {
  unsigned long iSurfacePoints;
  unsigned short iMarker;
  unsigned long negative_determinants = 0;

  /*--- Loop over the surface points ---*/
  for (iSurfacePoints = 0; iSurfacePoints < FFDBox->GetnSurfacePoint(); iSurfacePoints++) {
    /*--- Get the marker of the surface point ---*/
    iMarker = FFDBox->Get_MarkerIndex(iSurfacePoints);

    if (config->GetMarker_All_DV(iMarker) == YES) {
      const auto ParamCoord = FFDBox->Get_ParametricCoord(iSurfacePoints);

      /*--- Calculate partial derivatives ---*/
      unsigned short iDegree, jDegree, kDegree;
      su2double Ba, Bb, Bc, Ba_der, Bb_der, Bc_der;
      su2double determinant, d_du[3] = {0.0}, d_dv[3] = {0.0}, d_dw[3] = {0.0};

      for (iDegree = 0; iDegree <= FFDBox->lDegree; iDegree++) {
        Ba = FFDBox->BlendingFunction[0]->GetBasis(iDegree, ParamCoord[0]);
        Ba_der = FFDBox->BlendingFunction[0]->GetDerivative(iDegree, ParamCoord[0], 1);

        for (jDegree = 0; jDegree <= FFDBox->mDegree; jDegree++) {
          Bb = FFDBox->BlendingFunction[1]->GetBasis(jDegree, ParamCoord[1]);
          Bb_der = FFDBox->BlendingFunction[1]->GetDerivative(jDegree, ParamCoord[1], 1);

          for (kDegree = 0; kDegree <= FFDBox->nDegree; kDegree++) {
            Bc = FFDBox->BlendingFunction[2]->GetBasis(kDegree, ParamCoord[2]);
            Bc_der = FFDBox->BlendingFunction[2]->GetDerivative(kDegree, ParamCoord[2], 1);

            for (int i = 0; i < 3; ++i) {
              d_du[i] += Ba_der * Bb * Bc * FFDBox->Coord_Control_Points[iDegree][jDegree][kDegree][i];
              d_dv[i] += Ba * Bb_der * Bc * FFDBox->Coord_Control_Points[iDegree][jDegree][kDegree][i];
              d_dw[i] += Ba * Bb * Bc_der * FFDBox->Coord_Control_Points[iDegree][jDegree][kDegree][i];
            }
          }
        }
      }

      /*--- Calculate determinant ---*/
      determinant = d_du[0] * (d_dv[1] * d_dw[2] - d_dv[2] * d_dw[1]) -
                    d_dv[0] * (d_du[1] * d_dw[2] - d_du[2] * d_dw[1]) +
                    d_dw[0] * (d_du[1] * d_dv[2] - d_du[2] * d_dv[1]);

      if (determinant < 0) {
        negative_determinants++;
      }
    }
  }

  unsigned long tmp = negative_determinants;
  SU2_MPI::Allreduce(&tmp, &negative_determinants, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  return negative_determinants;
}
