/*!
 * \file CRadialBasisFunctionInterpolation.hpp
 * \brief Headers of the CRadialBasisFunctionInterpolation class.
 * \author F. van Steen
 * \version 8.0.1 "Harrier"
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
#include "CVolumetricMovement.hpp"
#include "CRadialBasisFunctionNode.hpp"
#include "../../include/toolboxes/CSymmetricMatrix.hpp"

/*!
 * \class CLinearElasticity
 * \brief Class for moving the volumetric numerical grid using Radial Basis Function interpolation.
 * \author F. van Steen
 */

class CRadialBasisFunctionInterpolation : public CVolumetricMovement {
protected:

  vector<CRadialBasisFunctionNode*> boundaryNodes;  /*!< \brief Vector with boundary nodes.*/
  vector<unsigned long> internalNodes;              /*!< \brief Vector with internal nodes.*/
  unsigned long nBoundaryNodes,                     /*!< \brief Number of boundary nodes*/
  nInternalNodes;                                   /*!< \brief Number of internal nodes*/

  vector<CRadialBasisFunctionNode*>* controlNodes;  /*!< \brief Vector with control nodes*/
  
  vector<passivedouble> deformationVector;  /*!< \brief Deformation vector.*/

  vector<passivedouble> coefficients;       /*!< \brief Control node interpolation coefficients.*/
  CSymmetricMatrix interpMat;               /*!< \brief Interpolation matrix.*/  
  // su2activematrix interpMat;

  RADIAL_BASIS kindRBF; /*!< \brief Type of Radial Basis Function.*/
  su2double radius;     /*!< \brief Support radius of compact Radial Basis Function.*/


  /*--- data reduction parameters ---*/
  vector<CRadialBasisFunctionNode*> greedyNodes;  /*!< \brief Vector with selected control nodes in greedy algorithm. */
  bool dataReduction;                             /*!< \brief Determines whether data reduction is used. */
  unsigned long MaxErrorNode;
  su2double MaxError;


  unsigned long Global_nControlNodes{0};
  /*--- mpi related*/
  
  unsigned long Local_nControlNodes;
  vector<unsigned long> Local_nControlNodesVec;

  vector<su2double> LocalCoords;
  vector<su2double> GlobalCoords;



  
public:

  /*!
  * \brief Constructor of the class.
  */
  CRadialBasisFunctionInterpolation(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CRadialBasisFunctionInterpolation(void) override;

  /*!
   * \brief Grid deformation using the spring analogy method.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] UpdateGeo - Update geometry.
   * \param[in] Derivative - Compute the derivative (disabled by default). Does not actually deform the grid if enabled.
   */
  void SetVolume_Deformation(CGeometry* geometry, CConfig* config, bool UpdateGeo, bool Derivative,
                                                bool ForwardProjectionDerivative);
  
  /*! 
  * \brief Obtaining the interpolation coefficients based on displacement of the control nodes.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] iNonlinear_Iter - Surface deformation step.
  */
  void GetInterpolationCoefficients(CGeometry* geometry, CConfig* config, unsigned long iNonlinear_Iter);

  /*!
  * \brief Compute the interpolation matrix with Radial Basis Function evaluations.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  */
  void SetInterpolationMatrix(CGeometry* geometry, CConfig* config);

  /*!
  * \brief Build the deformation vector with surface displacements of the control nodes.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  */
  void SetDeformationVector(CGeometry* geometry, CConfig* config);

  /*!
  * \brief Selecting unique set of control nodes based on marker information.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  */
  void SetControlNodes(CGeometry* geometry, CConfig* config);

  /*!
  * \brief Selecting internal nodes for the volumetric deformation.
  * \param[in] geometry - Geometrical definition of the problem.
  */
  void SetInternalNodes(CGeometry* geometry, CConfig* config);

  /*!
  * \brief Solving of the Radial Basis Function interpolation system, yielding the interpolation coefficients
  */
  void SolveRBF_System(void);

  
  /*!
  * \brief Updating the grid coordinates.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  */
  void UpdateGridCoord(CGeometry* geometry, CConfig* config);

  /*!
  * \brief Custom comparison function, for sorting the CRadialBasisFunctionNode objects based on their index.
  * \param[in] a - First considered Radial Basis Function Node.
  * \param[in] b - second considered Radial Basis Function Node.
  */
  inline static bool Compare(CRadialBasisFunctionNode* a, CRadialBasisFunctionNode* b){
    return a->GetIndex() < b->GetIndex();
  }

  /*!
  * \brief Custom equality function, for obtaining a unique set of CRadialBasisFunctionNode objects.
  * \param[in] a - First considered Radial Basis Function Node.
  * \param[in] b - second considered Radial Basis Function Node.
  */
  inline static bool Equal(CRadialBasisFunctionNode* a, CRadialBasisFunctionNode* b){
    return a->GetIndex() == b->GetIndex();
  }

  inline static bool Equal2(unsigned long a, unsigned long b){
      return a == b;    
  };


  /*--- Data reduction functions ---*/
  void GreedyIteration(CGeometry* geometry, CConfig* config);

  void GetInitMaxErrorNode(CGeometry* geometry, CConfig* config);
};