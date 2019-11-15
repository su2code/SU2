/*!
 * \file CTurbChannelSolution.hpp
 * \brief Header file for the class CTurbChannelSolution.
 *        The implementations are in the <i>CTurbChannelSolution.cpp</i> file.
 * \author E.Molina
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

#pragma once

#include <iomanip>
#include <cmath>
#include "CVerificationSolution.hpp"

/*!
 * \class CTurbChannelSolution
 * \brief Class to define the required data for the Turbulent Channel case
 * \author E. van der Weide,  E.Molina
 */
class CTurbChannelSolution: public CVerificationSolution {
  
protected:
  
  /*--- Specific conditions. ---*/
  
  su2double ReynoldsFriction;    /*!< \brief Friction Reynolds Number. */
  su2double ReynoldsMeanVelocity; /*!< \brief Reynolds Number. */
  su2double ovGm1, RGas;
  su2double halfChan;
  su2double alpha, a0;
  su2double rhoRef, uRef;
  su2double TMiddle, TWall;
  su2double uTau, lTau, tauWall, fBodyX;
  
  // Definition of the parameters for the numerical solution.
  const su2double convergenceThreshold = 1.e-12;
  const int nGridPoints     = 101;
  const su2double yPlusWall = (su2double) 0.5;

  // Constants of the SST turbulence model.
  const su2double sigmaK1  = (su2double) 0.85;
  const su2double sigmaK2  = (su2double) 1.0;
  const su2double sigmaOm1 = (su2double) 0.5;
  const su2double sigmaOm2 = (su2double) 0.856;
  const su2double beta1    = (su2double) 0.075;
  const su2double beta2    = (su2double) 0.0828;
  const su2double betaStar = (su2double) 0.09;
  const su2double kappa    = (su2double) 0.41;
  const su2double a1       = (su2double) 0.31;
  
  const su2double gam1 = beta1/betaStar - sigmaOm1*kappa*kappa/sqrt(betaStar);
  const su2double gam2 = beta2/betaStar - sigmaOm2*kappa*kappa/sqrt(betaStar);
  
  // Define the vectors to store the fully developed RANS solution.
  vector<su2double> yRANS, rhoRANS, uRANS, kRANS, omegaRANS;
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CTurbChannelSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CTurbChannelSolution(unsigned short val_nDim,
               unsigned short val_nvar,
               unsigned short val_iMesh,
               CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CTurbChannelSolution(void);
  
  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const su2double *val_coords,
                   const su2double val_t,
                   su2double       *val_solution);
  
  /*!
   * \brief Whether or not the exact solution is known for this verification solution.
   * \return  - False, because the exact solution is not known for the TGV case.
   */
  bool ExactSolutionKnown(void);
  
  void ComputeFullyDevelopedRANS(CConfig*       config,
                                 vector<su2double> &yRANS,
                                 vector<su2double> &rhoRANS,
                                 vector<su2double> &uRANS,
                                 vector<su2double> &kRANS,
                                 vector<su2double> &omegaRANS);
  
  void PointDistributionVinokur(const su2double len,
                                const su2double spacingEnds,
                                vector<su2double> &yPoints);
  void BlockThomas(const int n,
                   vector<su2double> &a,
                   vector<su2double> &b,
                   vector<su2double> &c,
                   vector<su2double> &d);

};
