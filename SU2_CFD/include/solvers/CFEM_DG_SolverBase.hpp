/*!
 * \file CFEM_DG_SolverBase.hpp
 * \brief Headers of the CFEM_DG_SolverBase class
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 7.1.1 "Blackbird"
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
 * \version 7.1.1 "Blackbird"
 */
class CFEM_DG_SolverBase : public CSolver {
protected:
  static constexpr size_t MAXNDIM = 3; /*!< \brief Max number of space dimensions, used in some static arrays. */

  vector<CFluidModel*> FluidModel;   /*!< \brief fluid model used in the solver. */

  su2double Mach_Inf        = 0.0;      /*!< \brief Mach number at the infinity. */
  su2double Density_Inf     = 0.0;      /*!< \brief Density at the infinity. */
  su2double Energy_Inf      = 0.0;      /*!< \brief Energy at the infinity. */
  su2double Temperature_Inf = 0.0;      /*!< \brief Energy at the infinity. */
  su2double Pressure_Inf    = 0.0;      /*!< \brief Pressure at the infinity. */
  su2double *Velocity_Inf   = nullptr;  /*!< \brief Flow Velocity vector at the infinity. */

  su2double Viscosity_Inf = 0; /*!< \brief Viscosity at the infinity. */
  su2double Tke_Inf = 0;       /*!< \brief Turbulent kinetic energy at the infinity. */
  su2double Prandtl_Lam = 0;   /*!< \brief Laminar Prandtl number. */
  su2double Prandtl_Turb = 0;  /*!< \brief Turbulent Prandtl number. */

  su2double StrainMag_Max = 0; /*!< \brief Maximum Strain Rate magnitude. */
  su2double Omega_Max = 0;     /*!< \brief Maximum Omega. */

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

  su2double AeroCoeffForceRef = 1.0;    /*!< \brief Reference force for aerodynamic coefficients. */
  su2double DynamicPressureRef = 1.0;   /*!< \brief Reference dynamic pressure. */

  unsigned long counter = 0;          /*!< \brief Shared variable to be able to carry out some
                                                  global counting in OpenMP. ---*/

  unsigned long nVolElemTot   = 0;    /*!< \brief Total number of local volume elements, including halos. */
  unsigned long nVolElemOwned = 0;    /*!< \brief Number of owned local volume elements. */

  unsigned long nDOFsLocOwned = 0;    /*!< \brief Number of owned local DOFs. */
  unsigned long nDOFsLocTot   = 0;    /*!< \brief Total number of local DOFs, including halos. */
  unsigned long nDOFsGlobal   = 0;    /*!< \brief Number of global DOFs. */

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

  bool symmetrizingTermsPresent = false; /*!< \brief Whether or not symmetrizing terms are present in the
                                                     discretization. */

  vector<unsigned long> nDOFsPerRank;    /*!< \brief Number of DOFs per rank in cumulative storage format. */

  vector<vector<unsigned long> > nonZeroEntriesJacobian; /*!< \brief The ID's of the DOFs for the
                                                                     non-zero entries of the Jacobian. */

  int nGlobalColors = 0;          /*!< \brief Number of global colors for the Jacobian computation. */

  vector<vector<unsigned long> > localDOFsPerColor;   /*!< \brief Double vector, which contains for every
                                                                  color the local DOFs. */
  vector<vector<int> > colorToIndEntriesJacobian;     /*!< \brief Double vector, which contains for every
                                                                  local DOF the mapping from the color to the
                                                                  entry in the Jacobian. A -1 indicates that
                                                                  the color does not contribute to the Jacobian
                                                                  of the DOF. */

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

  vector<CTaskDefinition> tasksList; /*!< \brief List of tasks to be carried out in the computationally
                                                 intensive part of the solver. */

  vector<bool> taskCompleted;        /*!< \brief Bool vector, that indicates whether or not the tasks
                                                 from the list tasksList have been completed. ---*/

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

  /*!
   * \brief Function, which returns the fluid model for this thread.
   * \return Pointer to the fluid model to be used.
   */
  inline CFluidModel* GetFluidModel(void) const final { return FluidModel[omp_get_thread_num()]; }

  /*!
   * \brief Function, which determines the communication pattern of the flow
   *        variables and the reverse pattern for the residuals.
   * \param[in] DGGeometry - Geometrical definition of the DG problem.
   * \param[in] config     - Definition of the particular problem.
   */
  void Prepare_MPI_Communication(const CMeshFEM_DG *DGGeometry,
                                 CConfig           *config);

  /*!
   * \brief Provide the non dimensional lift coefficient (inviscid contribution).
   * \param val_marker Surface where the coefficient is going to be computed.
   * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCL_Inv(unsigned short val_marker) const final { return InvCoeff.CL[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCD_Inv(unsigned short val_marker) const final { return InvCoeff.CD[val_marker]; }

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL(unsigned short val_marker) const final { return SurfaceCoeff.CL[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD(unsigned short val_marker) const final { return SurfaceCoeff.CD[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF(unsigned short val_marker) const final { return SurfaceCoeff.CSF[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff(unsigned short val_marker) const final { return SurfaceCoeff.CEff[val_marker]; }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx(unsigned short val_marker) const final { return SurfaceCoeff.CFx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy(unsigned short val_marker) const final { return SurfaceCoeff.CFy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz(unsigned short val_marker) const final { return SurfaceCoeff.CFz[val_marker]; }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx(unsigned short val_marker) const final { return SurfaceCoeff.CMx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy(unsigned short val_marker) const final { return SurfaceCoeff.CMy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz(unsigned short val_marker) const final { return SurfaceCoeff.CMz[val_marker]; }

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CL[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CD[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CSF[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff_Inv(unsigned short val_marker) const final {
    return SurfaceInvCoeff.CEff[val_marker];
  }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CFx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CFy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CFz[val_marker]; }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CMx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CMy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz_Inv(unsigned short val_marker) const final { return SurfaceInvCoeff.CMz[val_marker]; }

  /*!
   * \brief Provide the non dimensional sideforce coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCSF_Inv(unsigned short val_marker) const final { return InvCoeff.CSF[val_marker]; }

  /*!
   * \brief Provide the non dimensional efficiency coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the efficiency coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCEff_Inv(unsigned short val_marker) const final { return InvCoeff.CEff[val_marker]; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CSF() const final { return TotalCoeff.CSF; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CEff() const final { return TotalCoeff.CEff; }

  /*!
   * \brief Get the reference force used to compute CL, CD, etc.
   */
  inline su2double GetAeroCoeffsReferenceForce() const final { return AeroCoeffForceRef; }

  /*!
   * \brief Get the reference dynamic pressure, for Cp, Cf, etc.
   */
  inline su2double GetReferenceDynamicPressure() const final { return DynamicPressureRef; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
   * \return Value of the lift coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CL() const final { return TotalCoeff.CL; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CD() const final { return TotalCoeff.CD; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional x moment coefficient.
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMx() const final { return TotalCoeff.CMx; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y moment coefficient.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMy() const final { return TotalCoeff.CMy; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z moment coefficient.
   * \return Value of the moment z coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMz() const final { return TotalCoeff.CMz; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional x moment coefficient.
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CoPx() const final { return TotalCoeff.CoPx; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y moment coefficient.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CoPy() const final { return TotalCoeff.CoPy; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z moment coefficient.
   * \return Value of the moment z coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CoPz() const final { return TotalCoeff.CoPz; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional x force coefficient.
   * \return Value of the force x coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFx() const final { return TotalCoeff.CFx; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y force coefficient.
   * \return Value of the force y coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFy() const final { return TotalCoeff.CFy; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z force coefficient.
   * \return Value of the force z coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFz() const final { return TotalCoeff.CFz; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional thrust coefficient.
   * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CT() const final { return TotalCoeff.CT; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional thrust coefficient.
   * \param[in] val_Total_CT - Value of the total thrust coefficient.
   */
  inline void SetTotal_CT(su2double val_Total_CT) final { TotalCoeff.CT = val_Total_CT; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional torque coefficient.
   * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CQ() const final { return TotalCoeff.CQ; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional torque coefficient.
   * \param[in] val_Total_CQ - Value of the total torque coefficient.
   */
  inline void SetTotal_CQ(su2double val_Total_CQ) final { TotalCoeff.CQ = val_Total_CQ; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional rotor Figure of Merit.
   * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMerit() const final { return TotalCoeff.CMerit; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_CD(su2double val_Total_CD) final { TotalCoeff.CD = val_Total_CD; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional lift coefficient.
   * \param[in] val_Total_CL - Value of the total lift coefficient.
   */
  inline void SetTotal_CL(su2double val_Total_CL) final { TotalCoeff.CL = val_Total_CL; }

  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CL_Inv() const final { return AllBoundInvCoeff.CL; }

  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CD_Inv() const final { return AllBoundInvCoeff.CD; }

  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CSF_Inv() const final { return AllBoundInvCoeff.CSF; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CEff_Inv() const final { return AllBoundInvCoeff.CEff; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMx_Inv() const final { return AllBoundInvCoeff.CMx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMy_Inv() const final { return AllBoundInvCoeff.CMy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMz_Inv() const final { return AllBoundInvCoeff.CMz; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPx_Inv() const final { return AllBoundInvCoeff.CoPx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPy_Inv() const final { return AllBoundInvCoeff.CoPy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPz_Inv() const final { return AllBoundInvCoeff.CoPz; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFx_Inv() const final { return AllBoundInvCoeff.CFx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFy_Inv() const final { return AllBoundInvCoeff.CFy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFz_Inv() const final { return AllBoundInvCoeff.CFz; }

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CL[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD_Visc(unsigned short val_marker) const final { return SurfaceViscCoeff.CD[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF_Visc(unsigned short val_marker) const final {
    return SurfaceViscCoeff.CSF[val_marker];
  }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff_Visc(unsigned short val_marker) const final {
    return SurfaceViscCoeff.CEff[val_marker];
  }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx_Visc(unsigned short val_marker) const final {
    return SurfaceViscCoeff.CFx[val_marker];
  }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy_Visc(unsigned short val_marker) const final {
    return SurfaceViscCoeff.CFy[val_marker];
  }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz_Visc(unsigned short val_marker) const final {
    return SurfaceViscCoeff.CFz[val_marker];
  }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx_Visc(unsigned short val_marker) const final {
    return SurfaceViscCoeff.CMx[val_marker];
  }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy_Visc(unsigned short val_marker) const final {
    return SurfaceViscCoeff.CMy[val_marker];
  }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz_Visc(unsigned short val_marker) const final {
    return SurfaceViscCoeff.CMz[val_marker];
  }

  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CL_Visc() const final { return AllBoundViscCoeff.CL; }

  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CD_Visc() const final { return AllBoundViscCoeff.CD; }

  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CSF_Visc() const final { return AllBoundViscCoeff.CSF; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CEff_Visc() const final { return AllBoundViscCoeff.CEff; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMx_Visc() const final { return AllBoundViscCoeff.CMx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMy_Visc() const final { return AllBoundViscCoeff.CMy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMz_Visc() const final { return AllBoundViscCoeff.CMz; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPx_Visc() const final { return AllBoundViscCoeff.CoPx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPy_Visc() const final { return AllBoundViscCoeff.CoPy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPz_Visc() const final { return AllBoundViscCoeff.CoPz; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFx_Visc() const final { return AllBoundViscCoeff.CFx; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFy_Visc() const final { return AllBoundViscCoeff.CFy; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFz_Visc() const final { return AllBoundViscCoeff.CFz; }

  /*!
   * \brief Get the non dimensional lift coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCL_Visc(unsigned short val_marker) const final { return ViscCoeff.CL[val_marker]; }

  /*!
   * \brief Get the non dimensional sideforce coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the sideforce coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCSF_Visc(unsigned short val_marker) const final { return ViscCoeff.CSF[val_marker]; }

  /*!
   * \brief Get the non dimensional drag coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCD_Visc(unsigned short val_marker) const final { return ViscCoeff.CD[val_marker]; }

private:

  /*!
   * \brief Function, which computes the graph of the spatial discretization
            for the locally owned DOFs.
   * \param[in] DGGeometry - Geometrical definition of the DG problem.
   * \param[in] config     - Definition of the particular problem.
   */
  void DetermineGraphDOFs(const CMeshFEM_DG *DGGeometry,
                          CConfig           *config);

  /*!
   * \brief Function, which determines the meta data needed for the computation
   *        of the Jacobian of the spatial residual.
   * \param[in] DGGeometry     - Geometrical definition of the DG problem.
   * \param[in] colorLocalDOFs - Color of the locally stored DOFs.
   */
  void MetaDataJacobianComputation(const CMeshFEM_DG *DGGeometry,
                                   const vector<int> &colorLocalDOFs);
};
