/*!
 * \file CSurfaceMovement.hpp
 * \brief Headers of the CSurfaceMovement class.
 * \author F. Palacios, T. Economon.
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

#pragma once

#include "CGridMovement.hpp"
#include "CFreeFormDefBox.hpp"

/*!
 * \class CSurfaceMovement
 * \brief Class for moving the surface numerical grid.
 * \author F. Palacios, T. Economon.
 */
class CSurfaceMovement : public CGridMovement {
 protected:
  CFreeFormDefBox** FFDBox; /*!< \brief Definition of the Free Form Deformation Box. */
  unsigned short nFFDBox;   /*!< \brief Number of FFD FFDBoxes. */
  unsigned short nLevel;    /*!< \brief Level of the FFD FFDBoxes (parent/child). */
  bool FFDBoxDefinition;    /*!< \brief If the FFD FFDBox has been defined in the input file. */

 public:
  vector<su2double> GlobalCoordX[MAX_NUMBER_FFD];
  vector<su2double> GlobalCoordY[MAX_NUMBER_FFD];
  vector<su2double> GlobalCoordZ[MAX_NUMBER_FFD];
  vector<string> GlobalTag[MAX_NUMBER_FFD];
  vector<unsigned long> GlobalPoint[MAX_NUMBER_FFD];

  /*!
   * \brief Constructor of the class.
   */
  CSurfaceMovement(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CSurfaceMovement(void) override;

  /*!
   * \brief Set a Hicks-Henne deformation bump functions on an airfoil.
   * \param[in] boundary - Geometry of the boundary.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  void SetHicksHenne(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef);

  /*!
   * \brief Set a Hicks-Henne deformation bump functions on an airfoil.
   * \param[in] boundary - Geometry of the boundary.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  void SetSurface_Bump(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef);

  /*!
   * \brief Set a Hicks-Henne deformation bump functions on an airfoil.
   * \param[in] boundary - Geometry of the boundary.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  void SetAngleOfAttack(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef);

  /*!
   * \brief Set a deformation based on a change in the Kulfan parameters for an airfoil.
   * \param[in] boundary - Geometry of the boundary.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  void SetCST(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef);

  /*!
   * \brief Set a NACA 4 digits airfoil family for airfoil deformation.
   * \param[in] boundary - Geometry of the boundary.
   * \param[in] config - Definition of the particular problem.
   */
  void SetNACA_4Digits(CGeometry* boundary, CConfig* config);

  /*!
   * \brief Set a parabolic family for airfoil deformation.
   * \param[in] boundary - Geometry of the boundary.
   * \param[in] config - Definition of the particular problem.
   */
  void SetParabolic(CGeometry* boundary, CConfig* config);

  /*!
   * \brief Set a obstacle in a channel.
   * \param[in] boundary - Geometry of the boundary.
   * \param[in] config - Definition of the particular problem.
   */
  void SetAirfoil(CGeometry* boundary, CConfig* config);

  /*!
   * \brief Set a rotation for surface movement.
   * \param[in] boundary - Geometry of the boundary.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  void SetRotation(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef);

  /*!
   * \brief Computes the displacement of a rotating surface for a dynamic mesh simulation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iter - Current physical time iteration.
   * \param[in] iZone - Zone number in the mesh.
   */
  void HTP_Rotation(CGeometry* geometry, CConfig* config, unsigned long iter, unsigned short iZone);

  /*!
   * \brief Unsteady aeroelastic grid movement by deforming the mesh.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - Physical iteration number.
   * \param[in] iMarker - Marker to deform.
   * \param[in] iMarker_Monitoring - Marker we are monitoring.
   * \param[in] displacements - solution of typical section wing model.
   */
  void AeroelasticDeform(CGeometry* geometry, CConfig* config, unsigned long TimeIter, unsigned short iMarker,
                         unsigned short iMarker_Monitoring, vector<su2double>& displacements);

  /*!
   * \brief Deforms a 3-D flutter/pitching surface during an unsteady simulation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iter - Current physical time iteration.
   * \param[in] iZone - Zone number in the mesh.
   */
  void SetBoundary_Flutter3D(CGeometry* geometry, CConfig* config, CFreeFormDefBox** FFDBox, unsigned long iter,
                             unsigned short iZone);

  /*!
   * \brief Set the collective pitch for a blade surface movement.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCollective_Pitch(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Set any surface deformationsbased on an input file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Current physical time iteration.
   */
  void SetExternal_Deformation(CGeometry* geometry, CConfig* config, unsigned short iZone, unsigned long iter);

  /*!
   * \brief Set a displacement for surface movement.
   * \param[in] boundary - Geometry of the boundary.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  void SetTranslation(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef);

  /*!
   * \brief Set a displacement for surface movement.
   * \param[in] boundary - Geometry of the boundary.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  void SetScale(CGeometry* boundary, CConfig* config, unsigned short iDV, bool ResetDef);

  /*!
   * \brief Copy the boundary coordinates to each vertex.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void CopyBoundary(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Set the surface/boundary deformation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \return Total deformation applied, which may be less than target if intersection prevention is used.
   */
  vector<vector<su2double> > SetSurface_Deformation(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Compute the parametric coordinates of a grid point using a point inversion strategy
   *        in the free form FFDBox.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   */
  void SetParametricCoord(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, unsigned short iFFDBox);

  /*!
   * \brief Update the parametric coordinates of a grid point using a point inversion strategy
   *        in the free form FFDBox.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iFFDBox - Index of FFD box.
   */
  void UpdateParametricCoord(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, unsigned short iFFDBox);

  /*!
   * \brief Check the intersections of the FFD with the surface
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iFFDBox - Index of FFD box.
   */
  void CheckFFDIntersections(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, unsigned short iFFDBox);

  /*!
   * \brief Check the intersections of the FFD with the surface
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iFFDBox - Index of FFD box.
   */
  void CheckFFDDimension(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, unsigned short iFFDBox);

  /*!
   * \brief Set the Parametric coordinates.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBoxParent - Array with parent FFDBoxes of the computation.
   * \param[in] FFDBoxChild - Array with child FFDBoxes of the computation.
   */
  void SetParametricCoordCP(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBoxParent,
                            CFreeFormDefBox* FFDBoxChild);

  /*!
   * \brief Get the cartes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBoxParent - Array with parent FFDBoxes of the computation.
   * \param[in] FFDBoxChild - Array with child FFDBoxes of the computation.
   */
  void GetCartesianCoordCP(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBoxParent,
                           CFreeFormDefBox* FFDBoxChild);

  /*!
   * \brief Apply the design variables to the control point position
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iFFDBox - Index of FFD box.
   */
  void ApplyDesignVariables(CGeometry* geometry, CConfig* config, CFreeFormDefBox** FFDBox, unsigned short iFFDBox);

  /*!
   * \brief Recompute the cartesian coordinates using the control points position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iFFDBox - Index of FFD box.
   */
  su2double SetCartesianCoord(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, unsigned short iFFDBox,
                              bool ResetDef);

  /*!
   * \brief Set the deformation of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  bool SetFFDCPChange_2D(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, CFreeFormDefBox** ResetFFDBox,
                         unsigned short iDV, bool ResetDef) const;

  /*!
   * \brief Set the deformation of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  bool SetFFDCPChange(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, CFreeFormDefBox** ResetFFDBox,
                      unsigned short iDV, bool ResetDef) const;

  /*!
   * \brief Set the deformation of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  bool SetFFDGull(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, CFreeFormDefBox** ResetFFDBox,
                  unsigned short iDV, bool ResetDef) const;

  /*!
   * \brief Set the deformation of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  bool SetFFDNacelle(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, CFreeFormDefBox** ResetFFDBox,
                     unsigned short iDV, bool ResetDef) const;

  /*!
   * \brief Set a camber deformation of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  bool SetFFDCamber_2D(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, CFreeFormDefBox** ResetFFDBox,
                       unsigned short iDV, bool ResetDef) const;

  /*!
   * \brief Set a thickness deformation of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  bool SetFFDThickness_2D(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, CFreeFormDefBox** ResetFFDBox,
                          unsigned short iDV, bool ResetDef) const;

  /*!
   * \brief Set a camber deformation of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  bool SetFFDCamber(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, CFreeFormDefBox** ResetFFDBox,
                    unsigned short iDV, bool ResetDef) const;

  /*!
   * \brief Set a thickness deformation of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  bool SetFFDThickness(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, CFreeFormDefBox** ResetFFDBox,
                       unsigned short iDV, bool ResetDef) const;

  /*!
   * \brief Set a thickness deformation of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  void SetFFDAngleOfAttack(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, CFreeFormDefBox** ResetFFDBox,
                           unsigned short iDV, bool ResetDef);

  /*!
   * \brief Set a twist angle deformation of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  bool SetFFDTwist(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, CFreeFormDefBox** ResetFFDBox,
                   unsigned short iDV, bool ResetDef) const;

  /*!
   * \brief Set a rotation angle deformation of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  bool SetFFDRotation(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox, CFreeFormDefBox** ResetFFDBox,
                      unsigned short iDV, bool ResetDef) const;

  /*!
   * \brief Set a rotation angle deformation in a control surface of the Free From box using the control point position.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] iDV - Index of the design variable.
   * \param[in] ResetDef - Reset the deformation before starting a new one.
   */
  bool SetFFDControl_Surface(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox,
                             CFreeFormDefBox** ResetFFDBox, unsigned short iDV, bool ResetDef) const;

  /*!
   * \brief Read the free form information from the grid input file.
   * \note If there is no control point information, and no parametric
   *       coordinates information, the code will compute that information.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   * \param[in] val_mesh_filename - Name of the grid input file.
   */
  void ReadFFDInfo(CGeometry* geometry, CConfig* config, CFreeFormDefBox** FFDBox, const string& val_mesh_filename);

  /*!
   * \brief Read the free form information from the grid input file.
   * \note If there is no control point information, and no parametric
   *       coordinates information, the code will compute that information.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
   */
  void ReadFFDInfo(CGeometry* geometry, CConfig* config, CFreeFormDefBox** FFDBox);

  /*!
   * \brief Merge the Free Form information in the SU2 file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_mesh_filename - Name of the grid output file.
   */
  void MergeFFDInfo(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Write the Free Form information in the SU2 file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_mesh_filename - Name of the grid output file.
   */
  void WriteFFDInfo(CSurfaceMovement** surface_movement, CGeometry**** geometry, CConfig** config);

  /*!
   * \brief Get information about if there is a complete FFDBox definition, or it is necessary to
   *        compute the parametric coordinates.
   * \return <code>TRUE</code> if the input grid file has a complete information; otherwise <code>FALSE</code>.
   */
  inline bool GetFFDBoxDefinition(void) const { return FFDBoxDefinition; }

  /*!
   * \brief Check if the design variable definition matches the FFD box definition.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iDV - Index of the design variable.
   * \return <code>TRUE</code> if the FFD box name referenced with DV_PARAM can be found in the FFD box definition;
   * otherwise <code>FALSE</code>.
   */
  inline bool CheckFFDBoxDefinition(CConfig* config, unsigned short iDV) {
    for (unsigned short iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
      if (FFDBox[iFFDBox]->GetTag() == config->GetFFDTag(iDV)) {
        return true;
      }
    }
    return false;
  }

  /*!
   * \brief Obtain the number of FFDBoxes.
   * \return Number of FFD FFDBoxes.
   */
  inline unsigned short GetnFFDBox(void) const { return nFFDBox; }

  /*!
   * \brief Obtain the number of levels.
   * \return Number of FFD levels.
   */
  inline unsigned short GetnLevel(void) const { return nLevel; }

  /*!
   * \brief Set derivatives of the surface/boundary deformation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSurface_Derivative(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Calculate the determinant of the Jacobian matrix for the FFD problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] FFDBox - Free form deformation box.
   * \return Number of points with negative Jacobian determinant.
   */
  unsigned long calculateJacobianDeterminant(CGeometry* geometry, CConfig* config, CFreeFormDefBox* FFDBox) const;
};
