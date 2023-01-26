/*!
 * \file CFEM_DG_NSSolver.hpp
 * \brief Headers of the CFEM_DG_NSSolver class
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 7.5.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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
 * \version 7.5.0 "Blackbird"
 */
class CFEM_DG_NSSolver final : public CFEM_DG_EulerSolver {
private:
  CSGSModel *SGSModel = nullptr; /*!< \brief LES Subgrid Scale model. */
  bool SGSModelUsed = false;     /*!< \brief Whether or not an LES Subgrid Scale model is used. */
public:

  /*!
   * \brief Constructor of the class.
   */
  CFEM_DG_NSSolver(void) = default;

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
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
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elemBeg - Begin index of the element range to be computed.
   * \param[in]  elemEnd - End index (not included) of the element range to be computed.
   */
  void Shock_Capturing_DG(CConfig             *config,
                          const unsigned long elemBeg,
                          const unsigned long elemEnd) override;

  /*!
   * \brief Per-Olof Persson's method for capturing shock in DG
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   */
  void Shock_Capturing_DG_Persson(const unsigned long elemBeg,
                                  const unsigned long elemEnd);

  /*!
   * \brief Compute the volume contributions to the spatial residual.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elemBeg - Begin index of the element range to be computed.
   * \param[in]  elemEnd - End index (not included) of the element range to be computed.
   */
  void Volume_Residual(CConfig             *config,
                       const unsigned long elemBeg,
                       const unsigned long elemEnd) override;

  /*!
   * \brief Compute the spatial residual for the given range of faces.
   * \param[in] config      - Definition of the particular problem.
   * \param[in] indFaceBeg  - Starting index in the matching faces.
   * \param[in] indFaceEnd  - End index in the matching faces.
   * \param[in] numerics    - Description of the numerical method.
   */
  void ResidualFaces(CConfig             *config,
                     const unsigned long indFaceBeg,
                     const unsigned long indFaceEnd,
                     CNumerics           *numerics) override;

  /*!
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   */
  void BC_Euler_Wall(CConfig             *config,
                     const unsigned long surfElemBeg,
                     const unsigned long surfElemEnd,
                     CSurfaceElementFEM  *surfElem,
                     CNumerics           *conv_numerics) override;

  /*!
   * \brief Impose the far-field boundary condition.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   */
  void BC_Far_Field(CConfig             *config,
                    const unsigned long surfElemBeg,
                    const unsigned long surfElemEnd,
                    CSurfaceElementFEM  *surfElem,
                    CNumerics           *conv_numerics) override;

  /*!
   * \brief Impose the symmetry boundary condition using the residual.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   */
  void BC_Sym_Plane(CConfig             *config,
                    const unsigned long surfElemBeg,
                    const unsigned long surfElemEnd,
                    CSurfaceElementFEM  *surfElem,
                    CNumerics           *conv_numerics) override;

 /*!
   * \brief Impose the supersonic outlet boundary condition.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   */
  void BC_Supersonic_Outlet(CConfig             *config,
                            const unsigned long surfElemBeg,
                            const unsigned long surfElemEnd,
                            CSurfaceElementFEM  *surfElem,
                            CNumerics           *conv_numerics) override;

  /*!
   * \brief Impose the subsonic inlet boundary condition.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     val_marker    - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet(CConfig             *config,
                const unsigned long surfElemBeg,
                const unsigned long surfElemEnd,
                CSurfaceElementFEM  *surfElem,
                CNumerics           *conv_numerics,
                unsigned short      val_marker) override;

  /*!
   * \brief Impose the outlet boundary condition.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     val_marker    - Surface marker where the boundary condition is applied.
   */
  void BC_Outlet(CConfig             *config,
                 const unsigned long surfElemBeg,
                 const unsigned long surfElemEnd,
                 CSurfaceElementFEM  *surfElem,
                 CNumerics           *conv_numerics,
                 unsigned short      val_marker) override;

  /*!
   * \brief Impose a constant heat-flux condition at the wall.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     val_marker    - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CConfig             *config,
                        const unsigned long surfElemBeg,
                        const unsigned long surfElemEnd,
                        CSurfaceElementFEM  *surfElem,
                        CNumerics           *conv_numerics,
                        unsigned short      val_marker) override;

  /*!
   * \brief Impose an isothermal condition at the wall.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     val_marker    - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CConfig             *config,
                          const unsigned long surfElemBeg,
                          const unsigned long surfElemEnd,
                          CSurfaceElementFEM  *surfElem,
                          CNumerics           *conv_numerics,
                          unsigned short      val_marker) override;

  /*!
   * \brief Impose the boundary condition using characteristic reconstruction.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     val_marker    - Surface marker where the boundary condition is applied.
   */
  void BC_Riemann(CConfig             *config,
                  const unsigned long surfElemBeg,
                  const unsigned long surfElemEnd,
                  CSurfaceElementFEM  *surfElem,
                  CNumerics           *conv_numerics,
                  unsigned short      val_marker) override;

  /*!
   * \brief Impose the user customized boundary condition.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     surfElemBeg   - Start index in the list of surface elements.
   * \param[in]     surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in,out] surfElem      - Array of surface elements for which the boundary
                                    conditions must be imposed.
   * \param[in]     conv_numerics - Description of the numerical method.
   */
  void BC_Custom(CConfig             *config,
                 const unsigned long surfElemBeg,
                 const unsigned long surfElemEnd,
                 CSurfaceElementFEM  *surfElem,
                 CNumerics           *conv_numerics) override;

  /*!
   * \brief Compute the viscous forces and all the addimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Friction_Forces(const CGeometry* geometry, const CConfig* config) override;

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
                                           CVolumeElementFEM_DG *elem,
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
                                           CVolumeElementFEM_DG *elem,
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
                                              CVolumeElementFEM_DG *elem,
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
                                              CVolumeElementFEM_DG *elem,
                                              const su2double      *sol,
                                              const unsigned short nSimul,
                                              const unsigned short NPad,
                                              su2double            *res,
                                              su2double            *work) override;

  /*!
   * \brief Function, which computes the viscous (including penalty) fluxes
            and the symmetrizing fluxes in face points.
   * \param[in]     config             - Definition of the particular problem.
   * \param[in]     heatFluxPrescribed - Whether or not the heat flux is prescribed.
   * \param[in]     heatFlux           - Value of the heat flux, if prescribed.
   * \param[in]     lenScale           - Length scale for the penalty terms.
   * \param[in]     lenScaleLES        - Length scale for LES.
   * \param[in]     solLeft            - Left solution in the points.
   * \param[in]     solRight           - Right solution in the points.
   * \param[in]     gradSol            - The gradients of the entropy variables in the points.
   * \param[in]     JacobiansFace      - The jacobians of the face.
   * \param[in]     normalsFace        - The normals of the face.
   * \param[in,out] fluxes             - The viscous and penalty fluxes are added to
                                         the already stored inviscid fluxes.
   * \param[out]    symFluxes          - The symmetrizing fluxes.
   */
  void ComputeViscousFluxesFace(CConfig                            *config,
                                bool                               heatFluxPrescribed,
                                su2double                          heatFlux,
                                su2double                          lenScale,
                                su2double                          lenScaleLES,
                                ColMajorMatrix<su2double>          &solLeft,
                                ColMajorMatrix<su2double>          &solRight,
                                vector<ColMajorMatrix<su2double> > &gradSol,
                                su2activevector                    &JacobiansFace,
                                ColMajorMatrix<su2double>          &normalsFace,
                                ColMajorMatrix<su2double>          &fluxes,
                                vector<ColMajorMatrix<su2double> > &symFluxes);
};
