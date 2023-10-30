/*!
 * \file CBaslineSolver_FEM.hpp
 * \brief Headers of the CBaselineSolver class
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

#include "CSolver.hpp"

/*!
 * \class CBaselineSolver_FEM
 * \brief Main class for defining a baseline solution from a restart file for the DG-FEM solver output.
 * \author T. Economon.
 * \version 8.0.0 "Harrier"
 */
class CBaselineSolver_FEM final : public CSolver {
protected:

  unsigned long nDOFsLocTot;    /*!< \brief Total number of local DOFs, including halos. */
  unsigned long nDOFsLocOwned;  /*!< \brief Number of owned local DOFs. */
  unsigned long nDOFsGlobal;    /*!< \brief Number of global DOFs. */

  unsigned long nVolElemTot;    /*!< \brief Total number of local volume elements, including halos. */
  unsigned long nVolElemOwned;  /*!< \brief Number of owned local volume elements. */
  CVolumeElementFEM *volElem;   /*!< \brief Array of the local volume elements, including halos. */

  vector<su2double> VecSolDOFs;    /*!< \brief Vector, which stores the solution variables in all the DOFs. */

  CVariable* GetBaseClassPointerToNodes() override {return nullptr;}

public:

  /*!
   * \brief Constructor of the class.
   */
  CBaselineSolver_FEM(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CBaselineSolver_FEM(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CBaselineSolver_FEM(void) override;

  /*!
   * \brief Set the number of variables and string names from the restart file.
   * \param[in] config - Definition of the particular problem.
   */
  void SetOutputVariables(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry,
                   CSolver ***solver,
                   CConfig *config,
                   int val_iter,
                   bool val_update_geo) override;

  /*!
   * \brief Get a pointer to the vector of the solution degrees of freedom.
   * \return Pointer to the vector of the solution degrees of freedom.
   */
  inline su2double* GetVecSolDOFs(void) override { return VecSolDOFs.data(); }

};
