/*!
 * \file CFlowOutput.hpp
 * \brief  Headers of the flow output.
 * \author F. Palacios, T. Economon, M. Colonno
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "CFVMOutput.hpp"
#include "../variables/CVariable.hpp"

class CFlowOutput : public CFVMOutput{
protected:
  unsigned long lastInnerIter;

  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CFlowOutput(CConfig *config, unsigned short nDim, bool femOutput);

  /*!
   * \brief Add flow surface output fields
   * \param[in] config - Definition of the particular problem.
   */
  void AddAnalyzeSurfaceOutput(CConfig *config);

  /*!
   * \brief Set flow surface output field values
   * \param[in] solver - The container holding all solution data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] output - Boolean indicating whether information should be written to screen
   */
  void SetAnalyzeSurface(CSolver *solver, CGeometry *geometry, CConfig *config, bool output);

  /*!
   * \brief Add aerodynamic coefficients as output fields
   * \param[in] config - Definition of the particular problem.
   */
  void AddAerodynamicCoefficients(CConfig *config);

  /*!
   * \brief  Set the value of the aerodynamic coefficients
   * \param[in] config - Definition of the particular problem.
   * \param[in] flow_solver - The container holding all solution data.
   */
  void SetAerodynamicCoefficients(CConfig *config, CSolver *flow_solver);

  /*!
   * \brief  Set the value of the rotating frame coefficients (CT, CQ and CMerit).
   * \param[in] config - Definition of the particular problem.
   * \param[in] flow_solver - The container holding all solution data.
   */
  void SetRotatingFrameCoefficients(CConfig *config, CSolver *flow_solver);

  /*!
   * \brief Add CP inverse design output as history fields
   */
  void Add_CpInverseDesignOutput();

  /*!
   * \brief Set CP inverse design output field values (and also into the solver).
   * \param[in,out] solver - The container holding all solution data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_CpInverseDesign(CSolver *solver, const CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Compute value of the Q criteration for vortex idenfitication
   * \param[in] VelocityGradient - Velocity gradients
   * \return Value of the Q criteration at the node
   */
  template<class T>
  su2double GetQ_Criterion(const T& VelocityGradient) const {

    /*--- Make a 3D copy of the gradient so we do not have worry about nDim ---*/

    su2double Grad_Vel[3][3] = {{0.0}};

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        Grad_Vel[iDim][jDim] = VelocityGradient[iDim][jDim];

    /*--- Q Criterion Eq 1.2 of HALLER, G. (2005). An objective definition of a vortex.
     Journal of Fluid Mechanics, 525, 1-26. doi:10.1017/S0022112004002526 ---*/

    /*--- Components of the strain rate tensor (symmetric) ---*/
    su2double s11 = Grad_Vel[0][0];
    su2double s12 = 0.5 * (Grad_Vel[0][1] + Grad_Vel[1][0]);
    su2double s13 = 0.5 * (Grad_Vel[0][2] + Grad_Vel[2][0]);
    su2double s22 = Grad_Vel[1][1];
    su2double s23 = 0.5 * (Grad_Vel[1][2] + Grad_Vel[2][1]);
    su2double s33 = Grad_Vel[2][2];

    /*--- Components of the spin tensor (skew-symmetric) ---*/
    su2double omega12 = 0.5 * (Grad_Vel[0][1] - Grad_Vel[1][0]);
    su2double omega13 = 0.5 * (Grad_Vel[0][2] - Grad_Vel[2][0]);
    su2double omega23 = 0.5 * (Grad_Vel[1][2] - Grad_Vel[2][1]);

    /*--- Q = ||Omega|| - ||Strain|| ---*/
    su2double Q = 2*(pow(omega12,2) + pow(omega13,2) + pow(omega23,2)) -
      (pow(s11,2) + pow(s22,2) + pow(s33,2) + 2*(pow(s12,2) + pow(s13,2) + pow(s23,2)));

    return Q;
  }

  /*!
   * \brief Write information to meta data file
   * \param[in] config - Definition of the particular problem per zone.
   */
  void WriteMetaData(const CConfig *config);

  /*!
   * \brief Write any additional files defined for the current solver.
   * \param[in] config - Definition of the particular problem per zone.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - The container holding all solution data.
   */
  void WriteAdditionalFiles(CConfig *config, CGeometry *geometry, CSolver **solver_container) override;

  /*!
   * \brief Determines if the the volume output should be written.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Iter - Current iteration index.
   * \param[in] force_writing - boolean that forces writing of volume output
   */
  bool WriteVolume_Output(CConfig *config, unsigned long Iter, bool force_writing) override;

  /*!
   * \brief Write the forces breakdown file
   * \param[in] config - Definition of the particular problem per zone.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - The container holding all solution data.
   */
  void WriteForcesBreakdown(CConfig *config, CGeometry *geometry, CSolver **solver_container);

  /*!
   * \brief Set the time averaged output fields.
   */
  void SetTimeAveragedFields();

  /*!
   * \brief Load the time averaged output fields.
   * \param iPoint
   * \param node_flow
   */
  void LoadTimeAveragedData(unsigned long iPoint, CVariable *node_flow);

  /*!
   * \brief Write additional output for fixed CL mode.
   * \param[in] config - Definition of the particular problem per zone.
   */
  void SetFixedCLScreenOutput(const CConfig *config);

};
