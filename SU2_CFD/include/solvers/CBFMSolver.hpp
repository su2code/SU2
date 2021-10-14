/*!
 * \file CBFMSolver.hpp
 * \brief Headers of the CBFMSolver class
 * \author E.C. Bunschoten
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
#include "../variables/CBFMVariable.hpp"
#include "../numerics/ReadBFMInput.hpp"
#include "../numerics/BFMInterpolator.hpp"

/*!
 * \class CBFMSolver
 * \brief Main class for defining the Body-Force Model solver.
 * \ingroup BFM equations
 * \author E.C. Bunschoten
 */
class CBFMSolver final : public CSolver {
private:
  
  CBFMVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */
  ReadBFMInput* BFM_File_Reader;  /*! \brief Pointer to the file reader class which contains the information in the BFM input file. */
  BFMInterpolator* Interpolator;  /*! \brief Pointer to the interpolator class which interpolates the blade geometry parameters to the nodes. */
  vector<string> BFM_Parameter_Names{};  /*! \brief Vector containing the names of the blade geometric paremeters. */
  vector<unsigned short> BFM_Parameter_Indices{};
  su2double Omega{}; /*! \brief Rotation rate of the rotor blades in the BFM problem. */
  bool constant_viscosity{false};
  su2double mu_constant;
  su2double *Body_Force_Cart;

  unsigned short BFM_formulation{HALL}; /*! \brief BFM formulation used for analysis. Set to Hall as default. */
  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

  /*!
   * \brief Computes the relative flow velocity with respect to the local blade at iPoint.
   * \param[in] solver_container - pointer to solver container.
   * \param[in] iPoint - node index.
   * \param[in] W_cyl - Vector containing relative velocity component pointers.
  */
  void ComputeRelativeVelocity(CSolver **solver_container, unsigned long iPoint, vector<su2double*>&W_cyl);


  /*!
   * \brief Computes the cylindrical projections of cartesian node coordinates.
   * \param[in] geometry - pointer to geometry class.
   * \param[in] config - pointer to config class. 
  */
  void ComputeCylProjections(const CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Computes the flow source terms according to Halls body-force model at node iPoint
   * \param[in] solver_container - pointer to solver container.
   * \param[in] iPoint - node index.
   * \param[in] BFM_sources - vector containing BFM source terms.
   * \param[in] W_cyl - Vector containing relative velocity component pointers.
  */
  void ComputeBFMSources_Hall(CSolver **solver_container, unsigned long iPoint, vector<su2double>&BFM_sources, vector<su2double*>&W_cyl);

  su2double ComputeNormalForce_Thollet(CSolver **solver_container, unsigned long iPoint, su2double * W_cyl);

  su2double ComputeNormalForce_Hall(CSolver **solver_container, unsigned long iPoint, su2double * W_cyl);

  su2double ComputeParallelForce_Thollet(CSolver **solver_container, unsigned long iPoint, su2double * W_cyl);

  void ComputeBFM_Sources(CSolver **solver_container, unsigned long iPoint, vector<su2double>&BFM_sources, vector<su2double*>&W_cyl);
  /*!
   * \brief Computes the flow source terms according to Thollets body-force model at node iPoint
   * \param[in] solver_container - pointer to solver container.
   * \param[in] iPoint - node index.
   * \param[in] BFM_sources - vector containing BFM source terms.
   * \param[in] W_cyl - Vector containing relative velocity component pointers.
  */
  void ComputeBFMSources_Thollet(CSolver **solver_container, unsigned long iPoint, vector<su2double>&BFM_sources, vector<su2double*>&W_cyl);

  /*!
   * \brief Computes the flow source terms due to metal blockage
   * \param[in] solver_container - pointer to solver container.
   * \param[in] iPoint - node index.
   * \param[in] BFM_sources - vector containing BFM source terms.
  */
  void ComputeBlockageSources(CSolver **solver_container, unsigned long iPoint, vector<su2double>&BFM_sources);

  /*!
   * \brief Computes the compressibility factor in Thollets BFM
   * \param[in] solver_container - pointer to solver container.
   * \param[in] iPoint - node index.
   * \param[in] W_cyl - Vector containing relative velocity component pointers.
  */
  su2double ComputeKMach(CSolver **solver_container, unsigned long iPoint,  vector<su2double*>W_cyl);

  su2double ComputeKMach(CSolver **solver_container, unsigned long iPoint,  su2double * W_cyl);
public:

  /*!
   * \brief Constructor of the class.
   */
  CBFMSolver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CBFMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CBFMSolver(void) override;

  /*!
  * \brief Compute the source terms according to the BFM on the point iPoint
  * \param[in] solver_container - Container vector with all the solutions.
  * \param[in] iPoint - Point index at which the BFM source terms are to be calcluated.
  * \param[in] BFM_sources - Pointer to array where the BFM source terms are stored.
  */
  inline void ComputeBFMSources(CSolver **solver_container, unsigned long iPoint, vector<su2double>&BFM_sources) override;

  inline su2double GetBody_Force(unsigned short iDim){
    return Body_Force_Cart[iDim];
  }
};
