/*!
 * \file CFEM_DG_SolverBase.hpp
 * \brief Headers of the CFEM_DG_SolverBase class
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 7.1.0 "Blackbird"
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

#include "CSolver.hpp"

/*!
 * \class CFEM_DG_SolverBase
 * \brief Base class for defining a Discontinuous Galerkin finite element flow solver.
 * \ingroup Euler_Equations
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 7.1.0 "Blackbird"
 */
class CFEM_DG_SolverBase : public CSolver {
protected:
  static constexpr size_t MAXNDIM = 3; /*!< \brief Max number of space dimensions, used in some static arrays. */

  CFluidModel *FluidModel = nullptr; /*!< \brief fluid model used in the solver */

  su2double Mach_Inf        = 0.0;      /*!< \brief Mach number at the infinity. */
  su2double Density_Inf     = 0.0;      /*!< \brief Density at the infinity. */
  su2double Energy_Inf      = 0.0;      /*!< \brief Energy at the infinity. */
  su2double Temperature_Inf = 0.0;      /*!< \brief Energy at the infinity. */
  su2double Pressure_Inf    = 0.0;      /*!< \brief Pressure at the infinity. */
  su2double *Velocity_Inf   = nullptr;  /*!< \brief Flow Velocity vector at the infinity. */

  su2double Viscosity_Inf; /*!< \brief Viscosity at the infinity. */
  su2double Tke_Inf;       /*!< \brief Turbulent kinetic energy at the infinity. */

  su2double StrainMag_Max; /*!< \brief Maximum Strain Rate magnitude. */
  su2double Omega_Max;     /*!< \brief Maximum Omega. */

  vector<su2double> ConsVarFreeStream; /*!< \brief Vector, which contains the free stream
                                                   conservative variables. */

  AeroCoeffsArray InvCoeff;        /*!< \brief Inviscid pressure contributions for each boundary. */
  AeroCoeffsArray SurfaceInvCoeff; /*!< \brief Inviscid pressure contributions for each monitoring boundary. */
  AeroCoeffs AllBoundInvCoeff;     /*!< \brief Total pressure contribution for all the boundaries. */

  AeroCoeffsArray ViscCoeff;        /*!< \brief Viscous contributions for each boundary. */
  AeroCoeffsArray SurfaceViscCoeff; /*!< \brief Viscous contributions for each monitoring boundary. */
  AeroCoeffs AllBoundViscCoeff;     /*!< \brief Total pressure contribution for all the boundaries. */

  AeroCoeffsArray SurfaceCoeff; /*!< \brief Totals for each monitoring surface. */
  AeroCoeffs TotalCoeff;        /*!< \brief Totals for all boundaries. */

  unsigned long nVolElemTot   = 0;    /*!< \brief Total number of local volume elements, including halos. */
  unsigned long nVolElemOwned = 0;    /*!< \brief Number of owned local volume elements. */

  unsigned long nDOFsLocOwned = 0;    /*!< \brief Number of owned local DOFs. */
  unsigned long nDOFsLocTot   = 0;    /*!< \brief Total number of local DOFs, including halos. */
  unsigned long nDOFsGlobal   = 0;    /*!< \brief Number of global DOFs. */

  unsigned int sizeWorkArray = 0;     /*!< \brief The size of the work array needed. */

  CVolumeElementFEM_DG *volElem = nullptr;   /*!< \brief Array of the local volume elements, including halos. */

  const unsigned long *nVolElemOwnedPerTimeLevel = nullptr;    /*!< \brief Number of owned local volume elements
                                                                           per time level. Cumulative storage. */
  const unsigned long *nVolElemInternalPerTimeLevel = nullptr; /*!< \brief Number of internal local volume elements per
                                                                           time level. Internal means that the solution
                                                                 data does not need to be communicated. */
  const unsigned long *nVolElemHaloPerTimeLevel = nullptr;     /*!< \brief Number of halo volume elements
                                                                           per time level. Cumulative storage. */
 
  vector<vector<unsigned long> > ownedElemAdjLowTimeLevel; /*!< \brief List of owned elements per time level that are
                                                                       adjacent to elements of the lower time level. */
  vector<vector<unsigned long> > haloElemAdjLowTimeLevel;  /*!< \brief List of halo elements per time level that are
                                                                       adjacent to elements of the lower time level. */
  
  vector<vector<unsigned long> > adjMatchingFacesOfVolElem; /*!< \brief The ID's of the adjacent matching faces
                                                                        for all volume elements. */
  
  vector<vector<vector<unsigned long> > > adjSurfacesOfVolElem; /*!< \brief The ID's of the adjacent surfaces of the
                                                                            physical boundaries for all volume elements. */

  unsigned long nMeshPoints = 0;    /*!< \brief Number of mesh points in the local part of the grid. */
  CPointFEM *meshPoints = nullptr;  /*!< \brief Array of the points of the FEM mesh. */

  const unsigned long *nMatchingInternalFacesWithHaloElem = nullptr;  /*!< \brief Number of local matching internal faces per time level
                                                                                  between an owned and a halo element. Cumulative storage. */
  const unsigned long *nMatchingInternalFacesLocalElem    = nullptr;  /*!< \brief Number of local matching internal faces per time level
                                                                                  between local elements. Cumulative storage. */

  CInternalFaceFEM_DG *matchingInternalFaces = nullptr; /*!< \brief Array of the local matching internal faces. */
  CBoundaryFEM *boundaries = nullptr;                   /*!< \brief Array of the boundaries of the FEM mesh. */

  CSquareMatrixCM timeCoefADER_DG;                                      /*!< \brief The time coefficients in the iteration matrix
                                                                                    of the ADER-DG predictor step. */
  ColMajorMatrix<passivedouble> timeInterpolDOFToIntegrationADER_DG;    /*!< \brief The interpolation matrix between the time DOFs
                                                                                    and the time integration points for ADER-DG. */
  ColMajorMatrix<passivedouble> timeInterpolAdjDOFToIntegrationADER_DG; /*!< \brief The interpolation matrix between the time DOFs
                                                                                    of adjacent elements of a higher time level and
                                                                                    the time integration points for ADER-DG. */

  vector<su2double> TolSolADER;   /*!< \brief Vector, which stores the tolerances for the
                                              variables in the ADER predictor step. */

  bool symmetrizingTermsPresent = true;  /*!< \brief Whether or not symmetrizing terms are present in the
                                                     discretization. */

  vector<unsigned long> nDOFsPerRank;    /*!< \brief Number of DOFs per rank in cumulative storage format. */

  vector<vector<unsigned long> > nonZeroEntriesJacobian; /*!< \brief The ID's of the DOFs for the
                                                                     non-zero entries of the Jacobian. */

#ifdef HAVE_MPI
  vector<vector<SU2_MPI::Request> > commRequests;  /*!< \brief Communication requests in the communication of the solution for all
                                                               time levels. These are both sending and receiving requests. */

  vector<vector<vector<unsigned long> > > elementsRecvMPIComm;  /*!< \brief Triple vector, which contains the halo elements
                                                                            for MPI communication for all time levels. */
  vector<vector<vector<unsigned long> > > elementsSendMPIComm;  /*!< \brief Triple vector, which contains the donor elements
                                                                            for MPI communication for all time levels. */

  vector<vector<int> > ranksRecvMPI; /*!< \brief Double vector, which contains the ranks from which the halo elements
                                                 are received for all time levels. */
  vector<vector<int> > ranksSendMPI; /*!< \brief Double vector, which contains the ranks to which the donor elements
                                                 are sent for all time levels. */

  vector<vector<vector<su2double> > > commRecvBuf;  /*!< \brief Receive buffers used to receive the solution data
                                                                in the communication pattern for all time levels. */
  vector<vector<vector<su2double> > > commSendBuf;  /*!< \brief Send buffers used to send the solution data
                                                                in the communication pattern for all time levels. */
#endif

  vector<vector<unsigned long> > elementsRecvSelfComm;  /*!< \brief Double vector, which contains the halo elements
                                                                    for self communication for all time levels. */
  vector<vector<unsigned long> > elementsSendSelfComm;  /*!< \brief Double vector, which contains the donor elements
                                                                    for self communication for all time levels. */

  vector<su2double> rotationMatricesPeriodicity;    /*!< \brief Vector, which contains the rotation matrices
                                                                for the rotational periodic transformations. */
  vector<vector<vector<unsigned long> > > halosRotationalPeriodicity; /*!< \brief Triple vector, which contains the indices
                                                                                  of halo elements for which a periodic
                                                                                  transformation must be applied for all
                                                                                  time levels. */

private:

  CVariable* GetBaseClassPointerToNodes() final {return nullptr;}

public:

  /*!
   * \brief Constructor of the class.
   */
  CFEM_DG_SolverBase(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CFEM_DG_SolverBase(CGeometry      *geometry,
                     CConfig        *config,
                     unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_DG_SolverBase(void) override;

protected:

private:

  /*!
   * \brief Function, which computes the graph of the spatial discretization
            for the locally owned DOFs.
   * \param[in] DGGeometry - Geometrical definition of the DG problem.
   * \param[in] config     - Definition of the particular problem.
   */
  void DetermineGraphDOFs(const CMeshFEM_DG *DGGeometry,
                          CConfig           *config);
};
