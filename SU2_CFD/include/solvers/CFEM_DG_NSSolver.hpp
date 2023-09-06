/*!
 * \file CFEM_DG_NSSolver.hpp
 * \brief Headers of the CFEM_DG_NSSolver class
 * \author E. van der Weide, T. Economon, J. Alonso
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

#include "CFEM_DG_EulerSolver.hpp"

/*!
 * \class CFEM_DG_NSSolver
 * \brief Main class for defining the Navier-Stokes Discontinuous Galerkin finite element flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 8.0.0 "Harrier"
 */
class CFEM_DG_NSSolver final : public CFEM_DG_EulerSolver {
private:
  su2double Viscosity_Inf; /*!< \brief Viscosity at the infinity. */
  su2double Tke_Inf;       /*!< \brief Turbulent kinetic energy at the infinity. */
  su2double Prandtl_Lam,   /*!< \brief Laminar Prandtl number. */
  Prandtl_Turb;            /*!< \brief Turbulent Prandtl number. */

  CSGSModel *SGSModel;     /*!< \brief LES Subgrid Scale model. */
  bool SGSModelUsed;       /*!< \brief Whether or not an LES Subgrid Scale model is used. */

  su2double
  *CL_Visc,                 /*!< \brief Lift coefficient (viscous contribution) for each boundary. */
  *CD_Visc,                 /*!< \brief Drag coefficient (viscous contribution) for each boundary. */
  *CSF_Visc,                /*!< \brief Side force coefficient (viscous contribution) for each boundary. */
  *CMx_Visc,                /*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
  *CMy_Visc,                /*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
  *CMz_Visc,                /*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
  *CFx_Visc,                /*!< \brief Force x coefficient (viscous contribution) for each boundary. */
  *CFy_Visc,                /*!< \brief Force y coefficient (viscous contribution) for each boundary. */
  *CFz_Visc,                /*!< \brief Force z coefficient (viscous contribution) for each boundary. */
  *CEff_Visc,               /*!< \brief Efficiency (Cl/Cd) (Viscous contribution) for each boundary. */
  *Surface_CL_Visc,         /*!< \brief Lift coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CD_Visc,         /*!< \brief Drag coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CSF_Visc,        /*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CEff_Visc,       /*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFx_Visc,        /*!< \brief Force x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFy_Visc,        /*!< \brief Force y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFz_Visc,        /*!< \brief Force z coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMx_Visc,        /*!< \brief Moment x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMy_Visc,        /*!< \brief Moment y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMz_Visc,        /*!< \brief Moment z coefficient (viscous contribution) for each monitoring surface. */
  *Heat_Visc,               /*!< \brief Heat load (viscous contribution) for each boundary. */
  *MaxHeatFlux_Visc;        /*!< \brief Maximum heat flux (viscous contribution) for each boundary. */

  su2double
  AllBound_CD_Visc,         /*!< \brief Drag coefficient (viscous contribution) for all the boundaries. */
  AllBound_CL_Visc,         /*!< \brief Lift coefficient (viscous contribution) for all the boundaries. */
  AllBound_CSF_Visc,        /*!< \brief Sideforce coefficient (viscous contribution) for all the boundaries. */
  AllBound_CMx_Visc,        /*!< \brief Moment x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Visc,        /*!< \brief Moment y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Visc,        /*!< \brief Moment z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Visc,       /*!< \brief Efficient coefficient (Viscous contribution) for all the boundaries. */
  AllBound_CFx_Visc,        /*!< \brief Force x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Visc,        /*!< \brief Force y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Visc,        /*!< \brief Force z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_HeatFlux_Visc,     /*!< \brief Heat load (viscous contribution) for all the boundaries. */
  AllBound_MaxHeatFlux_Visc;  /*!< \brief Maximum heat flux (viscous contribution) for all boundaries. */
  su2double StrainMag_Max,
  Omega_Max;                  /*!< \brief Maximum Strain Rate magnitude and Omega. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CFEM_DG_NSSolver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_DG_NSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_DG_NSSolver(void) override;

  /*!
   * \brief Function to compute the time step for solving the Navier-Stokes equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Value of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry,
                    CSolver **solver_container,
                    CConfig *config,
                    unsigned short iMesh,
                    unsigned long Iteration) override;

  /*!
   * \brief Compute the artificial viscosity for shock capturing in DG.
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  void Shock_Capturing_DG(CConfig             *config,
                          const unsigned long elemBeg,
                          const unsigned long elemEnd,
                          su2double           *workArray) override;

  /*!
   * \brief Per-Olof Persson's method for capturing shock in DG
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  void Shock_Capturing_DG_Persson(const unsigned long elemBeg,
                                  const unsigned long elemEnd,
                                  su2double           *workArray);

  /*!
   * \brief Compute the volume contributions to the spatial residual.
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  void Volume_Residual(CConfig             *config,
                       const unsigned long elemBeg,
                       const unsigned long elemEnd,
                       su2double           *workArray) override;

  /*!
   * \brief Compute the spatial residual for the given range of faces.
   * \param[in]     config      - Definition of the particular problem.
   * \param[in]     indFaceBeg  - Starting index in the matching faces.
   * \param[in]     indFaceEnd  - End index in the matching faces.
   * \param[in,out] indResFaces - Index where to store the residuals in
                                  the vector of face residuals.
   * \param[in]     numerics    - Description of the numerical method.
   * \param[out]    workArray   - Work array.
   */
  void ResidualFaces(CConfig             *config,
                     const unsigned long indFaceBeg,
                     const unsigned long indFaceEnd,
                     unsigned long       &indResFaces,
                     CNumerics           *numerics,
                     su2double           *workArray) override;

  /*!
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray - Work array.
   */
  void BC_Euler_Wall(CConfig                  *config,
                     const unsigned long      surfElemBeg,
                     const unsigned long      surfElemEnd,
                     const CSurfaceElementFEM *surfElem,
                     su2double                *resFaces,
                     CNumerics                *conv_numerics,
                     su2double                *workArray) override;

  /*!
   * \brief Impose the far-field boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  void BC_Far_Field(CConfig                  *config,
                    const unsigned long      surfElemBeg,
                    const unsigned long      surfElemEnd,
                    const CSurfaceElementFEM *surfElem,
                    su2double                *resFaces,
                    CNumerics                *conv_numerics,
                    su2double                *workArray) override;

  /*!
   * \brief Impose the symmetry boundary condition using the residual.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  void BC_Sym_Plane(CConfig                  *config,
                    const unsigned long      surfElemBeg,
                    const unsigned long      surfElemEnd,
                    const CSurfaceElementFEM *surfElem,
                    su2double                *resFaces,
                    CNumerics                *conv_numerics,
                    su2double                *workArray) override;

 /*!
   * \brief Impose the supersonic outlet boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  void BC_Supersonic_Outlet(CConfig                  *config,
                            const unsigned long      surfElemBeg,
                            const unsigned long      surfElemEnd,
                            const CSurfaceElementFEM *surfElem,
                            su2double                *resFaces,
                            CNumerics                *conv_numerics,
                            su2double                *workArray) override;

  /*!
   * \brief Impose the subsonic inlet boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  void BC_Inlet(CConfig                  *config,
                const unsigned long      surfElemBeg,
                const unsigned long      surfElemEnd,
                const CSurfaceElementFEM *surfElem,
                su2double                *resFaces,
                CNumerics                *conv_numerics,
                unsigned short           val_marker,
                su2double                *workArray) override;

  /*!
   * \brief Impose the outlet boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  void BC_Outlet(CConfig                  *config,
                 const unsigned long      surfElemBeg,
                 const unsigned long      surfElemEnd,
                 const CSurfaceElementFEM *surfElem,
                 su2double                *resFaces,
                 CNumerics                *conv_numerics,
                 unsigned short           val_marker,
                 su2double                *workArray) override;

  /*!
   * \brief Impose a constant heat-flux condition at the wall.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  void BC_HeatFlux_Wall(CConfig                  *config,
                        const unsigned long      surfElemBeg,
                        const unsigned long      surfElemEnd,
                        const CSurfaceElementFEM *surfElem,
                        su2double                *resFaces,
                        CNumerics                *conv_numerics,
                        unsigned short           val_marker,
                        su2double                *workArray) override;

  /*!
   * \brief Impose an isothermal condition at the wall.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  void BC_Isothermal_Wall(CConfig                  *config,
                          const unsigned long      surfElemBeg,
                          const unsigned long      surfElemEnd,
                          const CSurfaceElementFEM *surfElem,
                          su2double                *resFaces,
                          CNumerics                *conv_numerics,
                          unsigned short           val_marker,
                          su2double                *workArray) override;

  /*!
   * \brief Impose the boundary condition using characteristic reconstruction.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  void BC_Riemann(CConfig                  *config,
                  const unsigned long      surfElemBeg,
                  const unsigned long      surfElemEnd,
                  const CSurfaceElementFEM *surfElem,
                  su2double                *resFaces,
                  CNumerics                *conv_numerics,
                  unsigned short           val_marker,
                  su2double                *workArray) override;

  /*!
   * \brief Impose the user customized boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  void BC_Custom(CConfig                  *config,
                 const unsigned long      surfElemBeg,
                 const unsigned long      surfElemEnd,
                 const CSurfaceElementFEM *surfElem,
                 su2double                *resFaces,
                 CNumerics                *conv_numerics,
                 su2double                *workArray) override;

  /*!
   * \brief Compute the viscosity at the infinity.
   * \return Value of the viscosity at the infinity.
   */
  inline su2double GetViscosity_Inf(void) const override { return Viscosity_Inf; }

  /*!
   * \brief Get the turbulent kinetic energy at the infinity.
   * \return Value of the turbulent kinetic energy at the infinity.
   */
  inline su2double GetTke_Inf(void) const override { return Tke_Inf; }

  /*!
   * \brief Compute the viscous forces and all the addimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Friction_Forces(const CGeometry* geometry, const CConfig* config) override;

  /*!
   * \brief Get the non dimensional lift coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCL_Visc(unsigned short val_marker) const override { return CL_Visc[val_marker]; }

  /*!
   * \brief Get the non dimensional sideforce coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the sideforce coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCSF_Visc(unsigned short val_marker) const override { return CSF_Visc[val_marker]; }

  /*!
   * \brief Get the non dimensional drag coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCD_Visc(unsigned short val_marker) const override { return CD_Visc[val_marker]; }

  /*!
   * \brief Get the total non dimensional lift coefficient (viscous contribution).
   * \return Value of the lift coefficient (viscous contribution).
   */
  inline su2double GetAllBound_CL_Visc() const override { return AllBound_CL_Visc; }

  /*!
   * \brief Get the total non dimensional sideforce coefficient (viscous contribution).
   * \return Value of the lift coefficient (viscous contribution).
   */
  inline su2double GetAllBound_CSF_Visc() const override { return AllBound_CSF_Visc; }

  /*!
   * \brief Get the total non dimensional drag coefficient (viscous contribution).
   * \return Value of the drag coefficient (viscous contribution).
   */
  inline su2double GetAllBound_CD_Visc() const override { return AllBound_CD_Visc; }

  /*!
   * \brief Get the max Omega.
   * \return Value of the max Omega.
   */
  inline su2double GetOmega_Max(void) const override { return Omega_Max; }

  /*!
   * \brief Get the max Strain rate magnitude.
   * \return Value of the max Strain rate magnitude.
   */
  inline su2double GetStrainMag_Max(void) const override { return StrainMag_Max; }

private:

  /*!
   * \brief Function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using an
            aliased discretization in 2D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  void ADER_DG_AliasedPredictorResidual_2D(CConfig              *config,
                                           CVolumeElementFEM    *elem,
                                           const su2double      *sol,
                                           const unsigned short nSimul,
                                           const unsigned short NPad,
                                           su2double            *res,
                                           su2double            *work) override;

/*!
   * \brief Function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using an
            aliased discretization in 3D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  void ADER_DG_AliasedPredictorResidual_3D(CConfig              *config,
                                           CVolumeElementFEM    *elem,
                                           const su2double      *sol,
                                           const unsigned short nSimul,
                                           const unsigned short NPad,
                                           su2double            *res,
                                           su2double            *work) override;
  /*!
   * \brief Function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using a
            non-aliased discretization in 2D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  void ADER_DG_NonAliasedPredictorResidual_2D(CConfig              *config,
                                              CVolumeElementFEM    *elem,
                                              const su2double      *sol,
                                              const unsigned short nSimul,
                                              const unsigned short NPad,
                                              su2double            *res,
                                              su2double            *work) override;

  /*!
   * \brief Function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using a
            non-aliased discretization in 3D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  void ADER_DG_NonAliasedPredictorResidual_3D(CConfig              *config,
                                              CVolumeElementFEM    *elem,
                                              const su2double      *sol,
                                              const unsigned short nSimul,
                                              const unsigned short NPad,
                                              su2double            *res,
                                              su2double            *work) override;
  /*!
   * \brief Function to compute the penalty terms in the integration
            points of a face.
   * \param[in]  indFaceChunk        - Index of the face in the chunk of fused faces.
   * \param[in]  nInt                - Number of integration points of the face.
   * \param[in]  NPad                - Value of the padding parameter to obtain optimal
                                       performance in the gemm computations.
   * \param[in]  solInt0             - Solution in the integration points of side 0.
   * \param[in]  solInt1             - Solution in the integration points of side 1.
   * \param[in]  viscosityInt0       - Viscosity in the integration points of side 0.
   * \param[in]  viscosityInt1       - Viscosity in the integration points of side 1.
   * \param[in]  kOverCvInt0         - Heat conductivity divided by Cv in the
                                       integration points of side 0.
   * \param[in]  kOverCvInt1         - Heat conductivity divided by Cv in the
                                       integration points of side 1.
   * \param[in]  ConstPenFace        - Penalty constant for this face.
   * \param[in]  lenScale0           - Length scale of the element of side 0.
   * \param[in]  lenScale1           - Length scale of the element of side 1.
   * \param[in]  metricNormalsFace   - Metric terms in the integration points, which
                                       contain the normals.
   * \param[out] penaltyFluxes       - Penalty fluxes in the integration points.
   */
  void PenaltyTermsFluxFace(const unsigned short indFaceChunk,
                            const unsigned short nInt,
                            const unsigned short NPad,
                            const su2double      *solInt0,
                            const su2double      *solInt1,
                            const su2double      *viscosityInt0,
                            const su2double      *viscosityInt1,
                            const su2double      *kOverCvInt0,
                            const su2double      *kOverCvInt1,
                            const su2double      ConstPenFace,
                            const su2double      lenScale0,
                            const su2double      lenScale1,
                            const su2double      *metricNormalsFace,
                                  su2double      *penaltyFluxes);

  /*!
   * \brief Function, which performs the treatment of the boundary faces for
            the Navier-Stokes equations for the most of the boundary conditions.
   * \param[in]     config                 - Definition of the particular problem.
   * \param[in]     conv_numerics          - Description of the numerical method.
   * \param[in]     nFaceSimul             - Number of faces that are treated simultaneously
                                             to improve performance.
   * \param[in]     NPad                   - Value of the padding parameter to obtain optimal
                                             performance in the gemm computations.
   * \param[in]     Wall_HeatFlux          - The value of the prescribed heat flux.
   * \param[in]     HeatFlux_Prescribed    - Whether or not the heat flux is prescribed by
                                             e.g. the boundary conditions.
   * \param[in]     Wall_Temperature       - The value of the prescribed wall temperature.
   * \param[in]     Temperature_Prescribed - Whether or not the temperature is precribed
                                             by e.g. the boundary conditions.
   * \param[in]     surfElem               - Surface boundary elements for which the
                                             residuals mut be computed.
   * \param[in]     solIntL                - Left states in the integration points of the face.
   * \param[in]     solIntR                - Right states in the integration points of the face.
   * \param[out]    workArray              - Storage for the local arrays.
   * \param[out]    resFaces               - Array to store the residuals of the face.
   * \param[in,out] indResFaces            - Index in resFaces, where the current residual
                                             should be stored. It is updated in the function
                                             for the next boundary element.
   * \param[in,out] wallModel              - Possible pointer to the wall model treatment.
                                             NULL pointer indicates no wall model treatment.
   */
  void ViscousBoundaryFacesBCTreatment(CConfig                  *config,
                                       CNumerics                *conv_numerics,
                                       const unsigned short     nFaceSimul,
                                       const unsigned short     NPad,
                                       const su2double          Wall_HeatFlux,
                                       const bool               HeatFlux_Prescribed,
                                       const su2double          Wall_Temperature,
                                       const bool               Temperature_Prescribed,
                                       const CSurfaceElementFEM *surfElem,
                                       const su2double          *solIntL,
                                       const su2double          *solIntR,
                                             su2double          *workArray,
                                             su2double          *resFaces,
                                             unsigned long      &indResFaces,
                                             CWallModel         *wallModel);

  /*!
   * \brief Function, which computes the viscous fluxes in the integration
            points for the boundary faces that must be treated simulaneously.
            This function uses the standard approach for computing the fluxes,
            i.e. no wall modeling.
   * \param[in]  config               - Definition of the particular problem.
   * \param[in]  nFaceSimul           - Number of faces that are treated simultaneously
                                        to improve performance.
   * \param[in]  NPad                 - Value of the padding parameter to obtain optimal
                                        performance in the gemm computations.
   * \param[in]  nInt                 - Number of integration points on the face.
   * \param[in]  nDOFsElem            - Number of DOFs of the adjacent element.
   * \param[in]  Wall_HeatFlux        - The value of the prescribed heat flux.
   * \param[in]  HeatFlux_Prescribed  - Whether or not the heat flux is prescribed by
                                        e.g. the boundary conditions.
   * \param[in]  derBasisElem         - Array, which contains the derivatives of the
                                        basis functions of the adjacent element
                                        in the integration points.
   * \param[in]  surfElem             - Surface boundary elements for which the
                                        viscous fluxes must be computed.
   * \param[in]  solIntL              - Left states in the integration points of the face.
   * \param[out] solElem              - Storage for the solution in the adjacent elements.
   * \param[out] gradSolInt           - Storage for the gradients of the solution in the
                                        integration points of the face.
   * \param[out] viscFluxes           - To be computed viscous fluxes in the
                                        integration points.
   * \param[out] viscosityInt         - To be computed viscosity in the integration points.
   * \param[out] kOverCvInt           - To be computed thermal conductivity in the
                                        integration points.
   */
  void ComputeViscousFluxesBoundaryFaces(CConfig                  *config,
                                         const unsigned short     nFaceSimul,
                                         const unsigned short     NPad,
                                         const unsigned short     nInt,
                                         const unsigned short     nDOFsElem,
                                         const su2double          Wall_HeatFlux,
                                         const bool               HeatFlux_Prescribed,
                                         const su2double          *derBasisElem,
                                         const CSurfaceElementFEM *surfElem,
                                         const su2double          *solIntL,
                                               su2double          *solElem,
                                               su2double          *gradSolInt,
                                               su2double          *viscFluxes,
                                               su2double          *viscosityInt,
                                               su2double          *kOverCvInt);

  /*!
   * \brief Function, which computes the viscous fluxes in the integration
            points for the boundary faces that must be treated simulaneously.
            The viscous fluxes are computed via a wall modeling approach.
   * \param[in]  config                 - Definition of the particular problem.
   * \param[in]  nFaceSimul             - Number of faces that are treated simultaneously
                                          to improve performance.
   * \param[in]  NPad                   - Value of the padding parameter to obtain optimal
                                          performance in the gemm computations.
   * \param[in]  nInt                   - Number of integration points on the face.
   * \param[in]  Wall_HeatFlux          - The value of the prescribed heat flux.
   * \param[in]  HeatFlux_Prescribed    - Whether or not the heat flux is prescribed by
                                          the boundary conditions.
   * \param[in]  Wall_Temperature       - The value of the prescribed wall temperature.
   * \param[in]  Temperature_Prescribed - Whether or not the temperature is precribed
                                          by  the boundary conditions
   * \param[in]  surfElem               - Surface boundary elements for which the
                                          viscous fluxes must be computed.
   * \param[in]  solIntL                - Left states in the integration points of the face.
   * \param[out] workArray              - Storage array
   * \param[out] viscFluxes             - To be computed viscous fluxes in the
                                          integration points.
   * \param[out] viscosityInt           - To be computed viscosity in the integration points.
   * \param[out] kOverCvInt             - To be computed thermal conductivity in the
                                          integration points.
   * \param[in,out] wallModel           - Pointer to the wall model treatment.
   */
  void WallTreatmentViscousFluxes(CConfig                  *config,
                                  const unsigned short     nFaceSimul,
                                  const unsigned short     NPad,
                                  const unsigned short     nInt,
                                  const su2double          Wall_HeatFlux,
                                  const bool               HeatFlux_Prescribed,
                                  const su2double          Wall_Temperature,
                                  const bool               Temperature_Prescribed,
                                  const CSurfaceElementFEM *surfElem,
                                  const su2double          *solIntL,
                                        su2double          *workArray,
                                        su2double          *viscFluxes,
                                        su2double          *viscosityInt,
                                        su2double          *kOverCvInt,
                                        CWallModel         *wallModel);

  /*!
   * \brief Function, which computes the residual contribution from a boundary
   face in a viscous computation when the boundary conditions have
   already been applied.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     nFaceSimul    - Number of fused faces, i.e. the number of faces
                                    that are treated simultaneously to improve performance.
   * \param[in]     NPad          - Value of the padding parameter to obtain optimal
                                    performance in the gemm computations.
   * \param[in]     surfElem      - Surface boundary elements for which the
                                    contribution to the residual must be computed.
   * \param[in]     solInt0       - Solution in the integration points of side 0.
   * \param[in]     solInt1       - Solution in the integration points of side 1.
   * \param[out]    paramFluxes   - Array used for temporary storage.
   * \param[out]    fluxes        - Temporary storage for the fluxes in the
                                    integration points.
   * \param[in,out] viscFluxes    - On input this array contains the viscous fluxes
                                    in the integration points. It is also used for
                                    temporary storage.
   * \param[in]     viscosityInt  - Temporary storage for the viscosity in the
                                    integration points.
   * \param[in]     kOverCvInt    - Temporary storage for the thermal conductivity
                                    over Cv in the integration points.
   * \param[out]    resFaces      - Array to store the residuals of the face.
   * \param[in,out] indResFaces   - Index in resFaces, where the current residual
                                    should be stored. It is updated in the function
                                    for the next boundary element.
   */
  void ResidualViscousBoundaryFace(CConfig                  *config,
                                   CNumerics                *conv_numerics,
                                   const unsigned short     nFaceSimul,
                                   const unsigned short     NPad,
                                   const CSurfaceElementFEM *surfElem,
                                   const su2double          *solInt0,
                                   const su2double          *solInt1,
                                   su2double                *paramFluxes,
                                   su2double                *fluxes,
                                   su2double                *viscFluxes,
                                   const su2double          *viscosityInt,
                                   const su2double          *kOverCvInt,
                                   su2double                *resFaces,
                                   unsigned long            &indResFaces);

  /*!
   * \brief Function to compute the symmetrizing terms in the integration
            points of a face.
   * \param[in]  indFaceChunk      - Index of the face in the chunk of fused faces.
   * \param[in]  nInt              - Number of integration points of the face.
   * \param[in]  NPad              - Value of the padding parameter to obtain optimal
                                     performance in the gemm computations.
   * \param[in]  solInt0           - Solution in the integration points of side 0.
   * \param[in]  solInt1           - Solution in the integration points of side 1.
   * \param[in]  viscosityInt0     - Viscosity in the integration points of side 0.
   * \param[in]  viscosityInt1     - Viscosity in the integration points of side 1.
   * \param[in]  kOverCvInt0       - Heat conductivity divided by Cv in the
                                     integration points of side 0.
   * \param[in]  kOverCvInt1       - Heat conductivity divided by Cv in the
                                     integration points of side 1.
   * \param[in]  metricNormalsFace - Metric terms in the integration points, which
                                     contain the normals.
   * \param[out] symmFluxes        - Symmetrizing fluxes in the integration points.
   */
  void SymmetrizingFluxesFace(const unsigned short indFaceChunk,
                              const unsigned short nInt,
                              const unsigned short NPad,
                              const su2double      *solInt0,
                              const su2double      *solInt1,
                              const su2double      *viscosityInt0,
                              const su2double      *viscosityInt1,
                              const su2double      *kOverCvInt0,
                              const su2double      *kOverCvInt1,
                              const su2double      *metricNormalsFace,
                                    su2double      *symmFluxes);

  /*!
   * \brief Function, which transforms the symmetrizing fluxes in the integration points
            such that they are suited to be multiplied by the parametric gradients of
            the basis functions.
   * \param[in]  indFaceChunk   - Index of the face in the chunk of fused faces.
   * \param[in]  nInt           - Number of integration points of the face.
   * \param[in]  NPad           - Value of the padding parameter to obtain optimal
                                  performance in the gemm computations.
   * \param[in]  halfTheta      - Half times the theta parameter in the symmetrizing terms.
   * \param[in]  symmFluxes     - Symmetrizing fluxes to be multiplied by the Cartesian
                                  gradients of the basis functions.
   * \param[in]  weights        - Integration weights of the integration points.
   * \param[in]  metricCoorFace - Derivatives of the parametric coordinates w.r.t. the
                                  Cartesian coordinates in the integration points of
                                  the face.
   * \param[out] paramFluxes    - Parametric fluxes in the integration points.
   */
  void TransformSymmetrizingFluxes(const unsigned short indFaceChunk,
                                   const unsigned short nInt,
                                   const unsigned short NPad,
                                   const su2double      halfTheta,
                                   const su2double      *symmFluxes,
                                   const su2double      *weights,
                                   const su2double      *metricCoorFace,
                                         su2double      *paramFluxes);

  /*!
   * \brief Function to compute the viscous normal fluxes in the integration points of a face.
   * \param[in]   adjVolElem          - Pointer to the adjacent volume.
   * \param[in]   indFaceChunk        - Index of the face in the chunk of fused faces.
   * \param[in]   nInt                - Number of integration points of the face.
   * \param[in]   NPad                - Value of the padding parameter to obtain optimal
                                        performance in the gemm computations.
   * \param[in]   Wall_HeatFlux       - The value of the prescribed heat flux.
   * \param[in]   HeatFlux_Prescribed - Whether or not the heat flux is prescribed by
                                        e.g. the boundary conditions.
   * \param[in]   solInt              - Solution in the integration points.
   * \param[in]   gradSolInt          - Gradient of the solution in the integration points.
   * \param[in]   metricCoorDerivFace - Metric terms in the integration points, which
                                        contain the derivatives of the parametric
                                        coordinates w.r.t. the Cartesian coordinates.
                                        Needed to compute the Cartesian gradients.
   * \param[in]   metricNormalsFace   - Metric terms in the integration points, which
                                        contain the normals.
   * \param[in]   wallDistanceInt     - Wall distances in the integration points of the face.
   * \param[out]  viscNormFluxes      - Viscous normal fluxes in the integration points.
   * \param[out]  viscosityInt        - Viscosity in the integration points, which is
                                        needed for other terms in the discretization.
   * \param[out]  kOverCvInt          - Thermal conductivity over Cv in the integration points,
                                        which is needed for other terms in the discretization.
   */
  void ViscousNormalFluxFace(const CVolumeElementFEM *adjVolElem,
                             const unsigned short    indFaceChunk,
                             const unsigned short    nInt,
                             const unsigned short    NPad,
                             const su2double         Wall_HeatFlux,
                             const bool              HeatFlux_Prescribed,
                             const su2double         *solInt,
                             const su2double         *gradSolInt,
                             const su2double         *metricCoorDerivFace,
                             const su2double         *metricNormalsFace,
                             const su2double         *wallDistanceInt,
                                   su2double         *viscNormFluxes,
                                   su2double         *viscosityInt,
                                   su2double         *kOverCvInt);

  /*!
   * \brief Function to compute the viscous normal flux in one integration point for a
            2D simulation.
   * \param[in]  sol            - Conservative variables.
   * \param[in]  solGradCart   - Cartesian gradients of the conservative variables.
   * \param[in]  normal        - Normal vector
   * \param[in]  HeatFlux      - Value of the prescribed heat flux. If not
                                 prescribed, this value should be zero.
   * \param[in]  factHeatFlux  - Multiplication factor for the heat flux. It is zero
                                 when the heat flux is prescribed and one when it has
                                 to be computed.
   * \param[in]  wallDist      - Distance to the nearest viscous wall, if appropriate.
   * \param[in   lenScale_LES  - LES length scale, if appropriate.
   * \param[out] Viscosity     - Total viscosity, to be computed.
   * \param[out] kOverCv       - Total thermal conductivity over Cv, to be computed.
   * \param[out] normalFlux    - Viscous normal flux, to be computed.
   */
  void ViscousNormalFluxIntegrationPoint_2D(const su2double *sol,
                                            const su2double solGradCart[4][2],
                                            const su2double *normal,
                                            const su2double HeatFlux,
                                            const su2double factHeatFlux,
                                            const su2double wallDist,
                                            const su2double lenScale_LES,
                                                  su2double &Viscosity,
                                                  su2double &kOverCv,
                                                  su2double *normalFlux);

  /*!
   * \brief Function to compute the viscous normal flux in one integration point for a
            3D simulation.
   * \param[in]  sol           - Conservative variables.
   * \param[in]  solGradCart   - Cartesian gradients of the conservative variables.
   * \param[in]  normal        - Normal vector
   * \param[in]  HeatFlux      - Value of the prescribed heat flux. If not
                                 prescribed, this value should be zero.
   * \param[in]  factHeatFlux  - Multiplication factor for the heat flux. It is zero
                                 when the heat flux is prescribed and one when it has
                                 to be computed.
   * \param[in]  wallDist      - Distance to the nearest viscous wall, if appropriate.
   * \param[in   lenScale_LES  - LES length scale, if appropriate.
   * \param[out] Viscosity     - Total viscosity, to be computed.
   * \param[out] kOverCv       - Total thermal conductivity over Cv, to be computed.
   * \param[out] normalFlux    - Viscous normal flux, to be computed.
   */
  void ViscousNormalFluxIntegrationPoint_3D(const su2double *sol,
                                            const su2double solGradCart[5][3],
                                            const su2double *normal,
                                            const su2double HeatFlux,
                                            const su2double factHeatFlux,
                                            const su2double wallDist,
                                            const su2double lenScale_LES,
                                                  su2double &Viscosity,
                                                  su2double &kOverCv,
                                                  su2double *normalFlux);
};
