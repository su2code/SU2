/*!
 * \file element_structure.hpp
 * \brief Headers of the Finite Element structure (elements)
 *        The subroutines and functions are in the <i>element_structure.cpp</i>
 *        and <i>element_linear.cpp</i> files.
 * \author R. Sanchez
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

#include "mpi_structure.hpp"

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "gauss_structure.hpp"

using namespace std;

/*!
 * \class CElement
 * \brief Main class for defining the element structure.
 * \author R. Sanchez
 */
class CElement {
protected:
  unsigned short nGaussPoints;    /*!< \brief Number of gaussian points. */
  unsigned short nNodes;          /*!< \brief Number of gaussian points. */
  static unsigned short nDim;     /*!< \brief Number of dimension of the problem. */
  su2double **CurrentCoord,       /*!< \brief Coordinates in the current frame. */
  **RefCoord;                     /*!< \brief Coordinates in the reference frame. */
  su2double **GaussCoord,         /*!< \brief Parent coordinates of the Gaussian Points. */
  *GaussWeight;                   /*!< \brief Weight of the Gaussian Points for the integration. */
  su2double **NodalExtrap;        /*!< \brief Coordinates of the nodal points for Gaussian extrapolation. */
  su2double **NodalStress;        /*!< \brief Stress at the nodes. */
  CGaussVariable **GaussPoint;    /*!< \brief Structure for the Gaussian Points. */
  su2double **Mab;                /*!< \brief Structure for the nodal components of the mass matrix. */
  su2double ***Kab;               /*!< \brief Structure for the constitutive component of the tangent matrix. */
  su2double **Ks_ab;              /*!< \brief Structure for the stress component of the tangent matrix. */
  su2double **Kt_a;               /*!< \brief Structure for the nodal stress term for the residual computation. */
  su2double **FDL_a;              /*!< \brief Structure for the dead loads for the residual computation. */
  su2double el_Pressure;          /*!< \brief Pressure in the element. */
  su2double ***dNiXj;             /*!< \brief Shape function derivatives. */
  unsigned short iDe;             /*!< \brief ID of the dielectric elastomer. */
  unsigned long iDV;              /*!< \brief ID of the Design Variable (if it is element based). */
  unsigned long iProp;            /*!< \brief ID of the Element Property. */
public:
  enum FrameType {REFERENCE=1, CURRENT=2}; /*!< \brief Type of nodal coordinates. */

protected:
  /*!
   * \brief Allocate element matrices and vectors, to be called by constructors of children classes.
   * \param[in] body_forces - If we need dead loads.
   */
  void AllocateStructures(const bool body_forces);

  /*!
   * \brief Compute gradients for 2D elements.
   * \param[in] mode - Type of coordinates to consider.
   */
  void ComputeGrad_2D(const FrameType mode);

  /*!
   * \brief Compute gradients for 3D elements.
   * \param[in] mode - Type of coordinates to consider.
   */
  void ComputeGrad_3D(const FrameType mode);

public:
  /*!
   * \brief Default constructor of the class.
   */
  CElement(void);

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CElement(unsigned short val_nDim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CElement(void);

  /*!
   * \brief Retrieve the number of nodes of the element.
   * \param[out] nNodes - Number of nodes of the element.
   */
  inline unsigned short GetnNodes(void) const {return nNodes;}

  /*!
   * \brief Retrieve the number of nodes of the element.
   * \param[out] nGaussPoints - Number of Gaussian Points of the element.
   */
  inline unsigned short GetnGaussPoints(void) const {return nGaussPoints;}

  /*!
   * \brief Set the value of the coordinate of the nodes in the reference configuration.
   * \param[in] val_CoordRef - Value of the coordinate.
   * \param[in] iNode - Number of node.
   * \param[in] iDim - Dimension
   */
  inline void SetRef_Coord(su2double val_CoordRef, unsigned short iNode, unsigned short iDim) {
    RefCoord[iNode][iDim] = val_CoordRef;
  }

  /*!
   * \brief Set the value of the coordinate of the nodes in the current configuration.
   * \param[in] val_CoordRef - Value of the coordinate.
   * \param[in] iNode - Number of node.
   * \param[in] iDim - Dimension
   */
  inline void SetCurr_Coord(su2double val_CoordCurr, unsigned short iNode, unsigned short iDim) {
    CurrentCoord[iNode][iDim] = val_CoordCurr;
  }

  /*!
   * \brief Set the value of the coordinate of the nodes in the reference configuration.
   * \param[in] val_CoordRef - Value of the coordinate.
   * \param[in] iNode - Number of node.
   * \param[in] iDim - Dimension
   * \param[out] Coordinate
   */
  inline su2double GetRef_Coord(unsigned short iNode, unsigned short iDim) const {
    return RefCoord[iNode][iDim];
  }

  /*!
   * \brief Get the value of the coordinate of the nodes in the current configuration.
   * \param[in] val_CoordRef - Value of the coordinate.
   * \param[in] iNode - Number of node.
   * \param[in] iDim - Dimension
   * \param[out] Coordinate
   */
  inline su2double GetCurr_Coord(unsigned short iNode, unsigned short iDim) const {
    return CurrentCoord[iNode][iDim];
  }

  /*!
   * \brief Get the weight of the corresponding Gaussian Point.
   * \param[in] iGauss - index of the Gaussian point.
   * \param[out] Weight.
   */
  inline su2double GetWeight(unsigned short iGauss) const {
    return GaussWeight[iGauss];
  }

  /*!
   * \brief Get the jacobian respect to the reference configuration for the Gaussian Point iGauss.
   * \param[in] iGauss - index of the Gaussian point.
   * \param[out] Weight.
   */
  inline su2double GetJ_X(unsigned short iGauss) const {
    return GaussPoint[iGauss]->GetJ_X();
  }

  /*!
   * \brief Get the jacobian respect to the current configuration for the Gaussian Point iGauss.
   * \param[in] iGauss - index of the Gaussian point.
   * \param[out] Weight.
   */
  inline su2double GetJ_x(unsigned short iGauss) const {
    return GaussPoint[iGauss]->GetJ_x();
  }

  /*!
   * \brief Retrieve the value of the pressure in the element for incompressible materials.
   * \param[out] Value of the pressure.
   */
  inline su2double GetElement_Pressure(void) const {
    return el_Pressure;
  }

  /*!
   * \brief Add the value of the diagonal term for the mass matrix.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \param[in] val_Ks_ab - value of the term that will constitute the diagonal of the stress contribution.
   */
  inline void Add_Mab(su2double val_Mab, unsigned short nodeA, unsigned short nodeB) {
    Mab[nodeA][nodeB] += val_Mab;
  }

  /*!
   * \brief Add the value of a submatrix K relating nodes a and b, for the constitutive term.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \param[in] val_Kab - value of the matrix K.
   */
  void Add_Kab(su2double **val_Kab, unsigned short nodeA, unsigned short nodeB);

  /*!
   * \brief Add the value of a submatrix K relating nodes a and b, for the constitutive term (symmetric terms need transpose)
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \param[in] val_Kab - value of the matrix K.
   */
  void Add_Kab_T(su2double **val_Kab, unsigned short nodeA, unsigned short nodeB);

  /*!
   * \brief Add the value of the diagonal term for the stress contribution to the stiffness of the system.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \param[in] val_Ks_ab - value of the term that will constitute the diagonal of the stress contribution.
   */
  inline void Add_Ks_ab(su2double val_Ks_ab, unsigned short nodeA, unsigned short nodeB) {
    Ks_ab[nodeA][nodeB] += val_Ks_ab;
  }

  /*!
   * \brief Add the value of the nodal stress term for the computation of the residual.
   * \param[in] nodeA - index of Node a.
   * \param[in] val_Kt_a - value of the term that will constitute the diagonal of the stress contribution.
   */
  void Add_Kt_a(su2double *val_Kt_a, unsigned short nodeA);

  /*!
   * \brief Add the value of the dead load for the computation of the residual.
   * \param[in] nodeA - index of Node a.
   * \param[in] val_FDL_a - value of the term that will constitute the diagonal of the stress contribution.
   */
  void Add_FDL_a(su2double *val_FDL_a, unsigned short nodeA);

  /*!
   * \brief Restarts the values in the element.
   */
  void clearElement(void);

  /*!
   * \brief Restarts the values of stress in the element.
   */
  void clearStress(void);

  /*!
   * \brief Return the value of the diagonal term for the mass matrix, relating nodes a and b.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \param[out] val_Mab - value of the diagonal term of Mab.
   */
  inline su2double Get_Mab(unsigned short nodeA, unsigned short nodeB) const {
    return Mab[nodeA][nodeB];
  }

  /*!
   * \brief Return the value of the submatrix K relating nodes a and b.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \param[out] val_Kab - value of the matrix K.
   */
  inline const su2double *Get_Kab(unsigned short nodeA, unsigned short nodeB) const {
    return Kab[nodeA][nodeB];
  }

  /*!
   * \brief Return the value of the diagonal term for the stress contribution, relating nodes a and b.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \param[out] val_Kab - value of the matrix K.
   */
  inline su2double Get_Ks_ab(unsigned short nodeA, unsigned short nodeB) const {
    return Ks_ab[nodeA][nodeB];
  }

  /*!
   * \brief Return the value of the nodal stress component of the residual for node a.
   * \param[in] nodeA - index of Node a.
   * \param[out] val_Kt_a - value of the stress term.
   */
  inline const su2double *Get_Kt_a(unsigned short nodeA) const {
    return Kt_a[nodeA];
  }

  /*!
   * \brief Return the value of the dead load component of the residual for node a.
   * \param[in] nodeA - index of Node a.
   * \param[out] val_Kt_a - value of the stress term.
   */
  inline const su2double *Get_FDL_a(unsigned short nodeA) const {
    return FDL_a[nodeA];
  }

  /*!
   * \brief Retrieve the value of the shape functions.
   * \param[in] iNode - Index of the node.
   * \param[in] iNode - Index of the Gaussian Point.
   * \param[out] GradNi_X - Gradient of the shape function related to node iNode and evaluated at Gaussian Point iGauss
   */
  inline su2double GetNi(unsigned short iNode, unsigned short iGauss) const {
    return GaussPoint[iGauss]->GetNi(iNode);
  }

  /*!
   * \brief Retrieve the value of the gradient of the shape functions respect to the reference configuration.
   * \param[in] iNode - Index of the node.
   * \param[in] iNode - Index of the Gaussian Point.
   * \param[out] GradNi_X - Gradient of the shape function related to node iNode and evaluated at Gaussian Point iGauss
   */
  inline su2double GetGradNi_X(unsigned short iNode, unsigned short iGauss, unsigned short iDim) const {
    return GaussPoint[iGauss]->GetGradNi_Xj(iNode,iDim);
  }

  /*!
   * \brief Retrieve the value of the gradient of the shape functions respect to the current configuration.
   * \param[in] iNode - Index of the node.
   * \param[in] iNode - Index of the Gaussian Point.
   * \param[out] GradNi_X - Gradient of the shape function related to node iNode and evaluated at Gaussian Point iGauss
   */
  inline su2double GetGradNi_x(unsigned short iNode, unsigned short iGauss, unsigned short iDim) const {
    return GaussPoint[iGauss]->GetGradNi_xj(iNode,iDim);
  }

  /*!
   * \brief Retrieve the value of the gradient of the shape functions respect to the reference configuration.
   * \param[in] iNode - Index of the node.
   * \param[in] iGauss - Index of the Gaussian Point.
   * \param[out] val_Ni_Ext - Value of the shape function at the nodes for extrapolation purposes
   */
  inline su2double GetNi_Extrap(unsigned short iNode, unsigned short iGauss) const {
    return NodalExtrap[iNode][iGauss];
  }

  /*!
   * \brief Add a value to the nodal stress for an element.
   * \param[in] iNode - Index of the node.
   * \param[in] iGauss - Index of the variable.
   * \param[in] val_Stress - Value of the stress added.
   */
  inline void Add_NodalStress(su2double val_Stress, unsigned short iNode, unsigned short iVar) {
    NodalStress[iNode][iVar] += val_Stress;
  }

  /*!
   * \brief Retrieve the value of the nodal stress for an element.
   * \param[in] iNode - Index of the node.
   * \param[in] iGauss - Index of the variable.
   * \param[in] val_Stress - Value of the stress added.
   */
  inline su2double Get_NodalStress(unsigned short iNode, unsigned short iVar) const {
    return NodalStress[iNode][iVar];
  }

  /*!
   * \brief Store the value of the identifier for the Dielectric Elastomers.
   * \param[in] val_iDe - identifier of the DE property.
   */
  void Set_ElProperties(const CProperty *element_property);

  /*!
   * \brief Store the value of the identifier for the Dielectric Elastomers.
   * \param[in] val_iDe - identifier of the DE property.
   */
  inline void Set_iDe(unsigned short val_iDe) {iDe = val_iDe;}

  /*!
   * \brief Return the value of the identifier for the Dielectric Elastomers.
   * \param[out] val_iDe - identifier of the DE property.
   */
  inline unsigned short Get_iDe(void) const {return iDe;}

  /*!
   * \brief Return the value of the identifier for the Design Variable.
   * \param[out] val_iDV - identifier of the DV.
   */
  inline unsigned long Get_iDV(void) const {return iDV;}

  /*!
   * \brief Return the value of the identifier for the Element Property.
   * \param[out] val_iProp - identifier of the property.
   */
  inline unsigned long Get_iProp(void) const {return iProp;}

  /*!
   * \brief Compute the value of the area of the element
   * \param[in] mode - Type of coordinates to consider in the computation
   * \param[out] val_Area - Area of the element
   */
  inline virtual su2double ComputeArea(const FrameType mode = REFERENCE) const {return 0.0;}

  /*!
   * \brief Compute the value of the volume of the element
   * \param[in] mode - Type of coordinates to consider in the computation
   * \param[out] val_Volume - Volume of the element
   */
  inline virtual su2double ComputeVolume(const FrameType mode = REFERENCE) const {return 0.0;}

  /*!
   * \brief Compute the value of the area of the element in current coordinates (wrapper to ComputeArea(CURRENT)).
   * \param[out] val_Area - Area of the element
   */
  inline su2double ComputeCurrentArea(void) const {return ComputeArea(CURRENT);}

  /*!
   * \brief Compute the value of the volume of the element in current coordinates (wrapper to ComputeVolume(CURRENT)).
   * \param[out] val_Volume - Volume of the element
   */
  inline su2double ComputeCurrentVolume(void) const {return ComputeVolume(CURRENT);}

  /*!
   * \brief Set the value of the gradient of the shape functions respect to the reference configuration.
   * \param[in] val_solution - Solution of the problem.
   * \param[out] J_X - Jacobian of the element evaluated at the current Gauss Point respect to the reference configuration
   */
  void ComputeGrad_Linear(void);

  /*!
   * \brief Set the value of the gradient of the shape functions respect to the current configuration.
   * \param[in] val_solution - Solution of the problem.
   * \param[out] J_x - Jacobian of the element evaluated at the current Gauss Point respect to the current configuration
   */
  void ComputeGrad_NonLinear(void);

  /*!
   * \brief Register the current and reference coordinates of the element as pre-accumulation inputs
   * the latter are needed for compatibility with shape derivatives, there is no problem registering
   * because inactive variables are ignored.
   */
  inline void SetPreaccIn_Coords(void) {
    AD::SetPreaccIn(RefCoord,nNodes,nDim);
    AD::SetPreaccIn(CurrentCoord,nNodes,nDim);
  }

  /*!
   * \brief Register the stress residual as a pre-accumulation output. When computing the element
   * stiffness matrix this is the only term that sees its way into the RHS of the system.
   */
  inline void SetPreaccOut_Kt_a(void) {
    AD::SetPreaccOut(Kt_a,nNodes,nDim);
  }

};

/*!
 * \class CTRIA1
 * \brief Tria element with 1 Gauss Points
 * \author R. Sanchez
 */
class CTRIA1 final : public CElement {

public:

  /*!
   * \brief Constructor of the class.
   */
  CTRIA1(void);

  /*!
   * \overload
   * \param[in] val_fea - Values of the fea solution (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTRIA1(unsigned short val_nDim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTRIA1(void);

  /*!
   * \brief Compute the value of the area of the element
   * \param[in] mode - Type of coordinates to consider in the computation
   * \param[out] val_Area - Area of the element
   */
  su2double ComputeArea(const FrameType mode = REFERENCE) const override;

};


/*!
 * \class CQUAD4
 * \brief Quadrilateral element with 4 Gauss Points
 * \author R. Sanchez
 */
class CQUAD4 final : public CElement {

public:

  /*!
   * \brief Constructor of the class.
   */
  CQUAD4(void);

  /*!
   * \overload
   * \param[in] val_fea - Values of the fea solution (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CQUAD4(unsigned short val_nDim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CQUAD4(void);

  /*!
   * \brief Compute the value of the area of the element
   * \param[in] mode - Type of coordinates to consider in the computation
   * \param[out] val_Area - Area of the element
   */
  su2double ComputeArea(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CTETRA1
 * \brief Tetrahedral element with 1 Gauss Point
 * \author R. Sanchez
 */
class CTETRA1 final : public CElement {

public:

  /*!
   * \brief Constructor of the class.
   */
  CTETRA1(void);

  /*!
   * \overload
   * \param[in] val_fea - Values of the fea solution (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTETRA1(unsigned short val_nDim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTETRA1(void);

  /*!
   * \brief Compute the value of the volume of the element
   * \param[out] val_Volume - Volume of the element
   */
  su2double ComputeVolume(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CHEXA8
 * \brief Hexahedral element with 8 Gauss Points
 * \author R. Sanchez
 */
class CHEXA8 final : public CElement {

public:

  /*!
   * \brief Constructor of the class.
   */
  CHEXA8(void);

  /*!
   * \overload
   * \param[in] val_fea - Values of the fea solution (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CHEXA8(unsigned short val_nDim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CHEXA8(void);

  /*!
   * \brief Compute the value of the volume of the element
   * \param[in] mode - Type of coordinates to consider in the computation
   * \param[out] val_Volume - Volume of the element
   */
  su2double ComputeVolume(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CPYRAM5
 * \brief Pyramid element with 5 Gauss Points
 * \author R. Sanchez, F. Palacios, A. Bueno, T. Economon, S. Padron.
 */
class CPYRAM5 final : public CElement {

public:

  /*!
   * \brief Constructor of the class.
   */
  CPYRAM5(void);

  /*!
   * \overload
   * \param[in] val_fea - Values of the fea solution (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPYRAM5(unsigned short val_nDim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CPYRAM5(void);

  /*!
   * \brief Compute the value of the volume of the element
   * \param[in] mode - Type of coordinates to consider in the computation
   * \param[out] val_Volume - Volume of the element
   */
  su2double ComputeVolume(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CPRISM6
 * \brief Prism element with 6 Gauss Points
 * \author R. Sanchez, F. Palacios, A. Bueno, T. Economon, S. Padron.
 * \version 6.2.0 "Falcon"
 */
class CPRISM6 final : public CElement {

public:

  /*!
   * \brief Constructor of the class.
   */
  CPRISM6(void);

  /*!
   * \overload
   * \param[in] val_fea - Values of the fea solution (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPRISM6(unsigned short val_nDim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CPRISM6(void);

  /*!
   * \brief Compute the value of the volume of the element
   * \param[in] mode - Type of coordinates to consider in the computation
   * \param[out] val_Volume - Volume of the element
   */
  su2double ComputeVolume(const FrameType mode = REFERENCE) const override;

};
