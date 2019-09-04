/*!
 * \file integration_structure.inl
 * \brief In-Line subroutines of the <i>integration_structure.hpp</i> file.
 * \author F. Palacios, T. Economon
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

inline su2double CIntegration::GetCauchy_Value(void) { return Cauchy_Value; }

inline bool CIntegration::GetConvergence(void) { return Convergence; }

inline bool CIntegration::GetConvergence_FSI(void) { return Convergence_FSI; }

inline bool CIntegration::GetConvergence_FullMG(void) { return Convergence_FullMG; }

inline void CIntegration::SetConvergence(bool value) { Convergence = value; }

inline void CIntegration::SetConvergence_FSI(bool valueFSI) { Convergence_FSI = valueFSI; }

inline void CIntegration::MultiGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container, 
                        CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) { }
  
inline void CIntegration::MultiGrid_Cycle(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container,
                 CConfig **config, unsigned short iMesh, unsigned short mu, unsigned short RunTime_EqSystem,
                 unsigned short iZone, unsigned short iInst) { }
                    
inline void CIntegration::NonDimensional_Parameters(CGeometry **geometry, CSolver ***solver_container, CNumerics ****numerics_container, 
                                                      CConfig *config, unsigned short FinestMesh, unsigned short RunTime_EqSystem, 
                                                      su2double *monitor) { }
  
inline void CIntegration::SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine, CConfig *config, unsigned short iMesh) { }

inline void CIntegration::SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                          CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) { }

inline void CIntegration::SetRestricted_Residual(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, 
                           CGeometry *geo_coarse, CConfig *config) { }

inline void CIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) { }

inline void CIntegration::SetRestricted_EddyVisc(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) { }

inline void CIntegration::SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, 
                         CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) { }
  
inline void CIntegration::SetResidual_Term(CGeometry *geometry, CSolver *flow) { }

inline void CIntegration::SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, 
                      CConfig *config, unsigned short iMesh) { }

inline void CIntegration::SingleGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container, 
                        CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) { }

inline void CIntegration::Structural_Iteration(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container, 
                        CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) { }


inline void CIntegration::SetPotential_Solver(CGeometry ****geometry, CSolver *****solver_container, CNumerics ******numerics_container, 
                                              CConfig **config, unsigned short RunTime_EqSystem, unsigned short iMesh, unsigned short iZone) { }
                                              
inline void CIntegration::Smooth_Solution(unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry, unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config) { }
