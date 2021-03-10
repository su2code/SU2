﻿/*!
 * \file CElement.hpp
 * \brief Main header of the Finite Element structure declaring the abstract
 *        interface and the available finite element types.
 * \author R. Sanchez
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

#include "CGaussVariable.hpp"
#include "CElementProperty.hpp"
#include <vector>

/*!
 * \class CElement
 * \brief Abstract class for defining finite elements.
 * \note Usage: Element instances are used to compute gradients (in reference or current
 *       coordinates), element matrices (stiffness, mass, etc.), and nodal residuals.
 *       The use-case is to set the coordinates, compute the gradients, and then
 *       incrementally build the matrices using the various "Add_XXX" methods.
 * \note Implementation: This base class provides a generic interface to an element,
 *       so that client numerics do not need to care about specific details (e.g. shape
 *       functions). However, more efficient code can be generated by caring about them,
 *       we try to strike a balance by tailoring the gradient computation to the specific
 *       sizes (number of integration points, nodes, and dimensions) of derived elements.
 *       This is done via the intermediate template class CElementWithKnownSizes from
 *       which derived elements can inherit the implementation of "ComputeGrad_XXX".
 *       To do so their constructors should consist of:
 *        - Call CElementWithKnownSizes's ctor to allocate parent class members.
 *        - Populate the "GaussWeights", and set the shape function values (Ni) in "GaussPoint".
 *        - Populate the "GaussCoords" and "dNiXj" fields of CElementWithKnownSizes.
 *        - Define the coordinates for nodal extrapolation of quantities at the Gauss points.
 *       Derived elements can also implement cheap methods of computing their area/volume.
 * \author R. Sanchez
 */
class CElement {
protected:
  std::vector<CGaussVariable> GaussPoint;  /*!< \brief Vector of Gaussian Points. */

  su2activematrix CurrentCoord;       /*!< \brief Coordinates in the current frame. */
  su2activematrix RefCoord;           /*!< \brief Coordinates in the reference frame. */
  su2activevector GaussWeight;        /*!< \brief Weight of the Gaussian Points for the integration. */
  su2activematrix NodalExtrap;        /*!< \brief Coordinates of the nodal points for Gaussian extrapolation. */
  su2activematrix NodalStress;        /*!< \brief Stress at the nodes. */

  /*--- Stiffness and load matrices. ---*/
  std::vector<su2activematrix> Kab;   /*!< \brief Structure for the constitutive component of the tangent matrix. */
  su2activematrix Mab;                /*!< \brief Structure for the nodal components of the mass matrix. */
  su2activematrix Ks_ab;              /*!< \brief Structure for the stress component of the tangent matrix. */
  su2activematrix Kt_a;               /*!< \brief Matrix of nodal stress terms for residual computation. */
  su2activematrix FDL_a;              /*!< \brief Matrix of dead loads for residual computation. */

  su2double el_Pressure = 0.0;        /*!< \brief Pressure in the element. */

  unsigned long iProp = 0;            /*!< \brief ID of the Element Property. */
  unsigned long iDV = 0;              /*!< \brief ID of the Design Variable (if it is element based). */
  unsigned long iDe = 0;              /*!< \brief ID of the dielectric elastomer. */

  unsigned short nGaussPoints;        /*!< \brief Number of gaussian points. */
  unsigned short nNodes;              /*!< \brief Number of geometric points. */
  unsigned short nDim;                /*!< \brief Number of dimension of the problem. */

  su2activematrix HiHj = 0.0;                     /*!< \brief Scalar product of 2 ansatz functions. */
  std::vector<std::vector<su2activematrix>> DHiDHj;             /*!< \brief Scalar product of the gradients of 2 ansatz functions. */

public:
  enum FrameType {REFERENCE=1, CURRENT=2}; /*!< \brief Type of nodal coordinates. */

  /*!
   * \brief Default constructor of the class, deleted to make sure derived
   *        classes always call the allocating constructor.
   */
  CElement(void) = delete;

  /*!
   * \brief Constructor of the class.
   * \param[in] ngauss - Number of Gaussian integration points.
   * \param[in] nnodes - Number of nodes of the element.
   * \param[in] ndim - Number of dimensions of the problem (2D, 3D).
   */
  CElement(unsigned short ngauss, unsigned short nnodes, unsigned short ndim);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CElement(void) = default;

  /*!
   * \brief Set the value of the gradient of the shape functions wrt the reference configuration.
   * \note The reference coordinates need to be set before calling this method.
   */
  virtual void ComputeGrad_Linear(void) = 0;

  /*!
   * \brief Set the value of the gradient of the shape functions wrt the current configuration.
   * \note The current coordinates need to be set before calling this method.
   */
  virtual void ComputeGrad_NonLinear(void) = 0;

  /*!
   * \brief Overload needed for deformed 2D elements on a surface in 3D.
   * \note The Coord are handed down from the numerics.
   */
  virtual void ComputeGrad_Linear(su2activematrix& Coord) = 0;

  /*!
   * \brief Sets matrices to 0.
   */
  void ClearElement(void);

  /*!
   * \brief Retrieve the number of nodes of the element.
   * \return Number of nodes of the element.
   */
  inline unsigned short GetnNodes(void) const {return nNodes;}

  /*!
   * \brief Retrieve the number of nodes of the element.
   * \return Number of Gaussian Points of the element.
   */
  inline unsigned short GetnGaussPoints(void) const {return nGaussPoints;}

  /*!
   * \brief Set the value of the coordinate of the nodes in the reference configuration.
   * \param[in] iNode - Index of node.
   * \param[in] iDim - Dimension.
   * \param[in] val_CoordRef - Value of the coordinate.
   */
  inline void SetRef_Coord(unsigned short iNode, unsigned short iDim, su2double val_CoordRef) {
    RefCoord(iNode,iDim) = val_CoordRef;
  }

  /*!
   * \brief Set the value of the coordinate of the nodes in the current configuration.
   * \param[in] iNode - Index of node.
   * \param[in] iDim - Dimension.
   * \param[in] val_CoordRef - Value of the coordinate.
   */
  inline void SetCurr_Coord(unsigned short iNode, unsigned short iDim, su2double val_CoordCurr) {
    CurrentCoord(iNode,iDim) = val_CoordCurr;
  }

  /*!
   * \brief Get the value of the coordinate of the nodes in the reference configuration.
   * \param[in] iNode - Index of node.
   * \param[in] iDim - Dimension.
   * \return Reference coordinate.
   */
  inline su2double GetRef_Coord(unsigned short iNode, unsigned short iDim) const {
    return RefCoord(iNode,iDim);
  }

  /*!
   * \brief Get the value of the coordinate of the nodes in the current configuration.
   * \param[in] iNode - Index of node.
   * \param[in] iDim - Dimension.
   * \return Current coordinate.
   */
  inline su2double GetCurr_Coord(unsigned short iNode, unsigned short iDim) const {
    return CurrentCoord(iNode,iDim);
  }

  /*!
   * \brief Get the weight of the corresponding Gaussian Point.
   * \param[in] iGauss - index of the Gaussian point.
   * \return Weight.
   */
  inline su2double GetWeight(unsigned short iGauss) const {
    return GaussWeight(iGauss);
  }

  /*!
   * \brief Get the Jacobian respect to the reference configuration for the Gaussian Point iGauss.
   * \param[in] iGauss - index of the Gaussian point.
   * \return Jacobian.
   */
  inline su2double GetJ_X(unsigned short iGauss) const {
    return GaussPoint[iGauss].GetJ_X();
  }

  /*!
   * \brief Get the jacobian respect to the current configuration for the Gaussian Point iGauss.
   * \param[in] iGauss - index of the Gaussian point.
   * \return Jacobian.
   */
  inline su2double GetJ_x(unsigned short iGauss) const {
    return GaussPoint[iGauss].GetJ_x();
  }

  /*!
   * \brief Retrieve the value of the pressure in the element for incompressible materials.
   * \return Value of the pressure.
   */
  inline su2double GetElement_Pressure(void) const { return el_Pressure; }

  /*!
   * \brief Add the value of the diagonal term for the mass matrix.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \param[in] val_Ks_ab - value of the term that will constitute the diagonal of the stress contribution.
   */
  inline void Add_Mab(unsigned short nodeA, unsigned short nodeB, su2double val_Mab) {
    Mab(nodeA,nodeB) += val_Mab;
  }

  /*!
   * \brief Add the value of a submatrix K relating nodes a and b, for the constitutive term.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \param[in] val_Kab - value of the matrix K.
   */
  inline void Add_Kab(unsigned short nodeA, unsigned short nodeB, su2double **val_Kab) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0; jDim < nDim; jDim++)
        Kab[nodeA](nodeB, iDim*nDim+jDim) += val_Kab[iDim][jDim];
  }

  /*!
   * \brief Add the value of a submatrix K relating nodes a and b, for the constitutive term (symmetric terms need transpose)
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \param[in] val_Kab - value of the matrix K.
   */
  inline void Add_Kab_T(unsigned short nodeA, unsigned short nodeB, su2double **val_Kab) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0; jDim < nDim; jDim++)
        Kab[nodeA](nodeB, iDim*nDim+jDim) += val_Kab[jDim][iDim];
  }

  /*!
   * \brief Add the value of the diagonal term for the stress contribution to the stiffness of the system.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \param[in] val_Ks_ab - value of the term that will constitute the diagonal of the stress contribution.
   */
  inline void Add_Ks_ab(unsigned short nodeA, unsigned short nodeB, su2double val_Ks_ab) {
    Ks_ab(nodeA,nodeB) += val_Ks_ab;
  }

  /*!
   * \brief Add the value of the nodal stress term for the computation of the residual.
   * \param[in] nodeA - index of Node a.
   * \param[in] val_Kt_a - value of the term that will constitute the diagonal of the stress contribution.
   */
  inline void Add_Kt_a(unsigned short nodeA, const su2double *val_Kt_a) {
    for(unsigned short iDim = 0; iDim < nDim; iDim++)
      Kt_a(nodeA,iDim) += val_Kt_a[iDim];
  }

  /*!
   * \brief Add the value of the dead load for the computation of the residual.
   * \param[in] nodeA - index of Node a.
   * \param[in] val_FDL_a - value of the term that will constitute the diagonal of the stress contribution.
   */
  inline void Add_FDL_a(unsigned short nodeA, const su2double *val_FDL_a) {
    for(unsigned short iDim = 0; iDim < nDim; iDim++)
      FDL_a(nodeA,iDim) += val_FDL_a[iDim];
  }

  /*!
   * \brief Restarts the values of stress in the element.
   */
  inline void ClearStress(void) {NodalStress.setConstant(0.0);}

  /*!
   * \brief Return the value of the diagonal term for the mass matrix, relating nodes a and b.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \return Value of the diagonal term of Mab.
   */
  inline su2double Get_Mab(unsigned short nodeA, unsigned short nodeB) const {
    return Mab(nodeA,nodeB);
  }

  /*!
   * \brief Return the value of the submatrix K relating nodes a and b.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \return Values of the matrix K.
   */
  inline const su2double *Get_Kab(unsigned short nodeA, unsigned short nodeB) const {
    return Kab[nodeA][nodeB];
  }

  /*!
   * \brief Return the value of the diagonal term for the stress contribution, relating nodes a and b.
   * \param[in] nodeA - index of Node a.
   * \param[in] nodeB - index of Node b.
   * \return Value of the matrix Ks.
   */
  inline su2double Get_Ks_ab(unsigned short nodeA, unsigned short nodeB) const {
    return Ks_ab(nodeA,nodeB);
  }

  /*!
   * \brief Return the values of the nodal stress components of the residual for node a.
   * \param[in] nodeA - index of Node a.
   * \return Values of the stress term.
   */
  inline const su2double *Get_Kt_a(unsigned short nodeA) const {
    return Kt_a[nodeA];
  }

  /*!
   * \brief Return the values of the dead load components of the residual for node a.
   * \param[in] nodeA - index of Node a.
   * \return Value of the dead loads.
   */
  inline const su2double *Get_FDL_a(unsigned short nodeA) const {
    return FDL_a[nodeA];
  }

  /*!
   * \brief Retrieve the value of the shape functions.
   * \param[in] iNode - Index of the node.
   * \param[in] iGauss - Index of the Gaussian Point.
   * \return Gradient of the shape function related to node iNode and evaluated at Gaussian Point iGauss
   */
  inline su2double GetNi(unsigned short iNode, unsigned short iGauss) const {
    return GaussPoint[iGauss].GetNi(iNode);
  }

  /*!
   * \brief Retrieve the value of the gradient of the shape functions respect to the reference configuration.
   * \param[in] iNode - Index of the node.
   * \param[in] iGauss - Index of the Gaussian Point.
   * \param[in] iDim - Coordinate index.
   * \return Gradient of the shape function related to node iNode and evaluated at Gaussian Point iGauss
   */
  inline su2double GetGradNi_X(unsigned short iNode, unsigned short iGauss, unsigned short iDim) const {
    return GaussPoint[iGauss].GetGradNi_Xj(iNode,iDim);
  }

  /*!
   * \brief Retrieve the value of the gradient of the shape functions respect to the current configuration.
   * \param[in] iNode - Index of the node.
   * \param[in] iGauss - Index of the Gaussian Point.
   * \param[in] iDim - Coordinate index.
   * \return Gradient of the shape function related to node iNode and evaluated at Gaussian Point iGauss
   */
  inline su2double GetGradNi_x(unsigned short iNode, unsigned short iGauss, unsigned short iDim) const {
    return GaussPoint[iGauss].GetGradNi_xj(iNode,iDim);
  }

  /*!
   * \brief Retrieve the value of the gradient of the shape functions respect to the reference configuration.
   * \param[in] iNode - Index of the node.
   * \param[in] iGauss - Index of the Gaussian Point.
   * \return Value of the shape function at the nodes for extrapolation purposes
   */
  inline su2double GetNi_Extrap(unsigned short iNode, unsigned short iGauss) const {
    return NodalExtrap(iNode,iGauss);
  }

  /*!
   * \brief Add a value to the nodal stress for an element.
   * \param[in] iNode - Index of the node.
   * \param[in] iVar - Variable index.
   * \param[in] val_Stress - Value of the stress added.
   */
  inline void Add_NodalStress(unsigned short iNode, unsigned short iVar, su2double val_Stress) {
    NodalStress(iNode,iVar) += val_Stress;
  }

  /*!
   * \brief Retrieve the value of the nodal stress for an element.
   * \param[in] iNode - Index of the node.
   * \param[in] iVar - Variable index.
   * \return Value of the stress.
   */
  inline su2double Get_NodalStress(unsigned short iNode, unsigned short iVar) const {
    return NodalStress(iNode,iVar);
  }

  /*!
   * \brief Store the values of the identifiers for element properties.
   * \param[in] element_property - element properties container.
   */
  inline void Set_ElProperties(const CProperty *element_property) {
    iDV   = element_property->GetDV();
    iProp = element_property->GetMat_Prop();
    iDe   = element_property->GetElectric_Prop();
  }

  /*!
   * \brief Store the value of the identifier for the Dielectric Elastomers.
   * \param[in] val_iDe - identifier of the DE property.
   */
  inline void Set_iDe(unsigned long val_iDe) {iDe = val_iDe;}

  /*!
   * \brief Return the value of the identifier for the Dielectric Elastomers.
   * \return Identifier of the DE property.
   */
  inline unsigned long Get_iDe(void) const {return iDe;}

  /*!
   * \brief Return the value of the identifier for the Design Variable.
   * \return Identifier of the DV.
   */
  inline unsigned long Get_iDV(void) const {return iDV;}

  /*!
   * \brief Return the value of the identifier for the Element Property.
   * \return Identifier of the property.
   */
  inline unsigned long Get_iProp(void) const {return iProp;}

  /*!
   * \brief Compute the value of the length of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Length of the (1D) element.
   */
  inline virtual su2double ComputeLength(const FrameType mode = REFERENCE) const {return 0.0;}

  /*!
   * \brief Compute the value of the area of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Area of the (2D) element.
   */
  inline virtual su2double ComputeArea(const FrameType mode = REFERENCE) const {return 0.0;}

  /*!
   * \brief Compute the value of the volume of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Volume of the (3D) element.
   */
  inline virtual su2double ComputeVolume(const FrameType mode = REFERENCE) const {return 0.0;}

  /*!
   * \brief Compute the value of the area of the element in current coordinates (wrapper to ComputeArea(CURRENT)).
   * \return Current area of the (2D) element.
   */
  inline su2double ComputeCurrentArea(void) const {return ComputeArea(CURRENT);}

  /*!
   * \brief Compute the value of the volume of the element in current coordinates (wrapper to ComputeVolume(CURRENT)).
   * \return Current volume of the (3D) element.
   */
  inline su2double ComputeCurrentVolume(void) const {return ComputeVolume(CURRENT);}

  /*!
   * \brief Register the current and reference coordinates of the element as pre-accumulation inputs
   * the latter are needed for compatibility with shape derivatives, there is no problem registering
   * because inactive variables are ignored.
   */
  inline void SetPreaccIn_Coords(bool nonlinear = true) {
    AD::SetPreaccIn(RefCoord.data(), nNodes*nDim);
    if (nonlinear)
      AD::SetPreaccIn(CurrentCoord.data(), nNodes*nDim);
  }

  /*!
   * \brief Register the stress residual as a pre-accumulation output. When computing the element
   * stiffness matrix this is the only term that sees its way into the RHS of the system.
   */
  inline void SetPreaccOut_Kt_a(void) {
    AD::SetPreaccOut(Kt_a.data(), nNodes*nDim);
  }

  /*!
   * \brief Register the mass matrix as a pre-accumulation output.
   */
  inline void SetPreaccOut_Mab(void) {
    AD::SetPreaccOut(Mab.data(), nNodes*nNodes);
  }

  /*!
   * \brief Register the dead load as a pre-accumulation output.
   */
  inline void SetPreaccOut_FDL_a(void) {
    AD::SetPreaccOut(FDL_a.data(), nNodes*nDim);
  }

 /*!
  * \brief Add the scalar product of the shape functions to the tangent matrix.
  * \param[in] nodeA - index of Node a.
  * \param[in] nodeB - index of Node b.
  * \param[in] val - value of the scalar product of ansatz function.
  */
 inline void Add_HiHj(su2double val, unsigned short nodeA, unsigned short nodeB) { HiHj[nodeA][nodeB] += val; }

 /*!
  * \brief Add the scalar product of the gradients of shape functions to the tangent matrix.
  * \param[in] nodeA - index of Node a.
  * \param[in] nodeB - index of Node b.
  * \param[in] val - value of the scalar product of gradients of ansatz function.
  */
 inline void Add_DHiDHj(su2double **val, unsigned short nodeA, unsigned short nodeB){
   unsigned short iDim, jDim;
   for(iDim = 0; iDim < nDim; iDim++) {
     for (jDim = 0; jDim < nDim; jDim++) {
       DHiDHj[nodeA][nodeB][iDim][jDim] += val[iDim][jDim];
     }
   }
 }

 /*!
  * \brief Add the transposed scalar product of the gradients of shape functions to the tangent matrix.
  * \param[in] nodeA - index of Node a.
  * \param[in] nodeB - index of Node b.
  * \param[in] val - value of the term that will contribute.
  */
 inline void Add_DHiDHj_T(su2double **val, unsigned short nodeA, unsigned short nodeB)  {
   unsigned short iDim, jDim;
   for(iDim = 0; iDim < nDim; iDim++) {
     for (jDim = 0; jDim < nDim; jDim++) {
       DHiDHj[nodeA][nodeB][iDim][jDim] += val[jDim][iDim];
     }
   }
 }

 /*!
  * \brief Get the scalar product of the shape functions to the tangent matrix.
  * \param[in] nodeA - index of Node a.
  * \param[in] nodeB - index of Node b.
  * \param[out] val - value of the scalar product of ansatz function.
  */
 inline su2double Get_HiHj(unsigned short nodeA, unsigned short nodeB)  { return HiHj[nodeA][nodeB]; }

 /*!
  * \brief Get the scalar product of the gradients of shape functions to the tangent matrix.
  * \param[in] nodeA - index of Node a.
  * \param[in] nodeB - index of Node b.
  * \return val - value of the scalar product of gradients of ansatz function.
  */
 inline su2activematrix& Get_DHiDHj(unsigned short nodeA, unsigned short nodeB) { return DHiDHj[nodeA][nodeB];}

};

/*!
 * \class CElementWithKnownSizes
 * \brief Templated class to implement the computation of gradients for specific element sizes.
 * \author P. Gomes, R. Sanchez
 */
template<unsigned short NGAUSS, unsigned short NNODE, unsigned short NDIM>
class CElementWithKnownSizes : public CElement {
private:

  FORCEINLINE static su2double JacobianAdjoint(const su2double Jacobian[][1], su2double[][1]) {
    /*--- Adjoint to Jacobian, we put 1.0 here so that ad/detJac is the inverse later ---*/
    ad[0][0] = 1.0;
    /*--- Determinant of Jacobian ---*/
    detJac = Jacobian[0][0];
  }


  FORCEINLINE static su2double JacobianAdjoint(const su2double Jacobian[][2], su2double ad[][2]) {
    ad[0][0] =  Jacobian[1][1];  ad[0][1] = -Jacobian[0][1];
    ad[1][0] = -Jacobian[1][0];  ad[1][1] =  Jacobian[0][0];
    /*--- Determinant of Jacobian ---*/
    return ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
  }

  FORCEINLINE static su2double JacobianAdjoint(const su2double Jacobian[][3], su2double ad[][3]) {
    ad[0][0] = Jacobian[1][1]*Jacobian[2][2]-Jacobian[1][2]*Jacobian[2][1];
    ad[0][1] = Jacobian[0][2]*Jacobian[2][1]-Jacobian[0][1]*Jacobian[2][2];
    ad[0][2] = Jacobian[0][1]*Jacobian[1][2]-Jacobian[0][2]*Jacobian[1][1];
    ad[1][0] = Jacobian[1][2]*Jacobian[2][0]-Jacobian[1][0]*Jacobian[2][2];
    ad[1][1] = Jacobian[0][0]*Jacobian[2][2]-Jacobian[0][2]*Jacobian[2][0];
    ad[1][2] = Jacobian[0][2]*Jacobian[1][0]-Jacobian[0][0]*Jacobian[1][2];
    ad[2][0] = Jacobian[1][0]*Jacobian[2][1]-Jacobian[1][1]*Jacobian[2][0];
    ad[2][1] = Jacobian[0][1]*Jacobian[2][0]-Jacobian[0][0]*Jacobian[2][1];
    ad[2][2] = Jacobian[0][0]*Jacobian[1][1]-Jacobian[0][1]*Jacobian[1][0];
    /*--- Determinant of Jacobian ---*/
    return Jacobian[0][0]*ad[0][0]+Jacobian[0][1]*ad[1][0]+Jacobian[0][2]*ad[2][0];
  }

protected:
  static_assert(NDIM==1 || NDIM==2 || NDIM==3, "ComputeGrad_impl expects 1D, 2D or 3D");

  su2double GaussCoord[NGAUSS][NDIM];   /*!< \brief Coordinates of the integration points. */
  su2double dNiXj[NGAUSS][NNODE][NDIM]; /*!< \brief Shape function derivatives evaluated at the Gauss points. */

  /*!
   * \brief Constructor of the class.
   * \note This constructor is protected to prevent this class from being instantiated.
   *       Its only role is to ensure the parent constructor is called with the right arguments.
   */
  CElementWithKnownSizes() : CElement(NGAUSS, NNODE, NDIM) {}

  /*!
   * \brief Implementation of gradient computation leveraging the static sizes.
   * \param[in] FRAME - template, REFERENCE or CURRENT coordinates.
   */
  template<FrameType FRAME>
  void ComputeGrad_impl(void) {

    su2double Jacobian[NDIM][NDIM], ad[NDIM][NDIM];
    unsigned short iNode, iDim, jDim, iGauss;

    /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
          for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed) ---*/
    const su2activematrix& Coord = (FRAME==REFERENCE) ? RefCoord : CurrentCoord;

    for (iGauss = 0; iGauss < NGAUSS; iGauss++) {

      /*--- Jacobian transformation ---*/
      /*--- This does dX/dXi transpose ---*/

      for (iDim = 0; iDim < NDIM; iDim++)
        for (jDim = 0; jDim < NDIM; jDim++)
          Jacobian[iDim][jDim] = 0.0;

      for (iNode = 0; iNode < NNODE; iNode++)
        for (iDim = 0; iDim < NDIM; iDim++)
          for (jDim = 0; jDim < NDIM; jDim++)
            Jacobian[iDim][jDim] += Coord(iNode,jDim) * dNiXj[iGauss][iNode][iDim];

      if (NDIM==1) {
        /*--- Obviously the Jacobian is the slope of the line ---*/
        Jacobian[0][0] = (Coord[1][0] - Coord[0][0]);
      }

      /*--- Adjoint to the Jacobian and determinant ---*/

      auto detJac = JacobianAdjoint(Jacobian, ad);

      if (FRAME==REFERENCE)
        GaussPoint[iGauss].SetJ_X(detJac);
      else
        GaussPoint[iGauss].SetJ_x(detJac);

      /*--- Jacobian inverse (it was already computed as transpose) ---*/

      for (iDim = 0; iDim < NDIM; iDim++)
        for (jDim = 0; jDim < NDIM; jDim++)
          Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;

      /*--- Derivatives with respect to global coordinates ---*/

      for (iNode = 0; iNode < NNODE; iNode++) {
        for (iDim = 0; iDim < NDIM; iDim++) {
          su2double GradNi_Xj = 0.0;
          for (jDim = 0; jDim < NDIM; jDim++)
            GradNi_Xj += Jacobian[iDim][jDim] * dNiXj[iGauss][iNode][jDim];

          if (FRAME==REFERENCE)
            GaussPoint[iGauss].SetGradNi_Xj(GradNi_Xj, iDim, iNode);
          else
            GaussPoint[iGauss].SetGradNi_xj(GradNi_Xj, iDim, iNode);
        }
      }

    }

  }

  template<FrameType FRAME>
  void ComputeGrad_impl_embedded(su2activematrix& Coord) {

    unsigned short iNode, iDim, iGauss;

    if (NDIM==1) {

      su2double Jacobian[2];
      su2double val_grad;

      Jacobian[0] = (Coord[1][0] - Coord[0][0]);
      Jacobian[1] = (Coord[1][1] - Coord[0][1]);
      const su2double JTJ = Jacobian[0]*Jacobian[0]+Jacobian[1]*Jacobian[1];
      const su2double volJacobian = sqrt(JTJ);

      for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

        GaussPoint[iGauss].SetJ_X(volJacobian);

        for (iNode = 0; iNode < nNodes; iNode++) {
          for (iDim=0; iDim<2; iDim++) {
            val_grad = Jacobian[iDim]*dNiXj[iGauss][iNode][0]/JTJ;
            GaussPoint[iGauss].SetGradNi_Xj( val_grad, iDim, iNode);
          }
        }
      }

    } else if (NDIM==2) {

      su2double Jacobian[3][2], JTJ[2][2];
      su2double volJacobian, GradN_Xi;
      unsigned short iShape, iDim, jDim, jVertex, kDim, iGauss;

      for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

        /*--- Jacobian transformation ---*/
        /*--- This does dX/dXi transpose ---*/

        for (iDim = 0; iDim < 3; iDim++) {
          for (jVertex = 0; jVertex < 2; jVertex++) {
              Jacobian[iDim][jVertex] = Coord[jVertex+1][iDim] - Coord[0][iDim];
          }
        }

        for (iDim = 0; iDim < 2; iDim++) {
          for (jDim = 0; jDim < 2; jDim++) {
            JTJ[iDim][jDim] = 0.0;
            for (kDim = 0; kDim < 3; kDim++) {
              JTJ[iDim][jDim] += Jacobian[kDim][iDim]*Jacobian[kDim][jDim];
            }
          }
        }

        /*--- matrix volume of Jacobian ---*/
        const su2double detJacobian=JTJ[0][0]*JTJ[1][1]-JTJ[0][1]*JTJ[1][0];

        volJacobian = sqrt(detJacobian);

        GaussPoint[iGauss].SetJ_X(volJacobian);

        /*--- Derivatives with respect to global coordinates ---*/

        su2double JTJinv[2][2];
        JTJinv[0][0]=JTJ[1][1]/detJacobian;
        JTJinv[1][1]=JTJ[0][0]/detJacobian;
        JTJinv[0][1]=-JTJ[0][1]/detJacobian;
        JTJinv[1][0]=-JTJ[1][0]/detJacobian;

        su2double Jdagger[2][3];

        for (iDim = 0; iDim < 2; iDim++) {
          for (jDim = 0; jDim < 3; jDim++) {
            Jdagger[iDim][jDim] = 0.0;
            for (kDim = 0; kDim < 2; kDim++) {
              Jdagger[iDim][jDim] += JTJinv[iDim][kDim]*Jacobian[jDim][kDim];
            }
          }
        }

        su2double GradOut[3];
        for (iShape = 0; iShape < nNodes; iShape++) {
          for (iDim = 0; iDim < 3; iDim++) {
            GradOut[iDim] = 0.0;
            for (jDim = 0; jDim < 2; jDim++) {
              GradOut[iDim] += Jdagger[jDim][iDim]*dNiXj[iGauss][iShape][jDim];
            }
            GaussPoint[iGauss].SetGradNi_Xj(GradOut[iDim], iDim, iShape);
          }
        }

      }

    }

  }


public:
  /*!
   * \brief Set the value of the gradient of the shape functions wrt the reference configuration.
   */
  void ComputeGrad_Linear(void) final { ComputeGrad_impl<REFERENCE>(); }

  /*!
   * \brief Set the value of the gradient of the shape functions wrt the current configuration.
   */
  void ComputeGrad_NonLinear(void) final { ComputeGrad_impl<CURRENT>(); }

  /*!
   * \brief Overload needed for deformed 2D elements on a surface in 3D.
   */
  void ComputeGrad_Linear(su2activematrix& Coord) final { ComputeGrad_impl_embedded<REFERENCE>(Coord); }


};

/*!
 * \class CTRIA1
 * \brief Tria element with 1 Gauss Points
 * \author R. Sanchez
 */
class CTRIA1 final : public CElementWithKnownSizes<1,3,2> {
private:
  enum : unsigned short {NGAUSS = 1};
  enum : unsigned short {NNODE = 3};
  enum : unsigned short {NDIM = 2};

public:
  /*!
   * \brief Constructor of the class.
   */
  CTRIA1();

  /*!
   * \brief Compute the value of the area of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Area of the element.
   */
  su2double ComputeArea(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CQUAD4
 * \brief Quadrilateral element with 4 Gauss Points
 * \author R. Sanchez
 */
class CQUAD4 final : public CElementWithKnownSizes<4,4,2> {
private:
  enum : unsigned short {NGAUSS = 4};
  enum : unsigned short {NNODE = 4};
  enum : unsigned short {NDIM = 2};

public:
  /*!
   * \brief Constructor of the class.
   */
  CQUAD4();

  /*!
   * \brief Shape functions (Ni) evaluated at point Xi,Eta.
   */
  inline static void ShapeFunctions(su2double Xi, su2double Eta, su2double* Ni) {
    Ni[0] = 0.25*(1.0-Xi)*(1.0-Eta);
    Ni[1] = 0.25*(1.0+Xi)*(1.0-Eta);
    Ni[2] = 0.25*(1.0+Xi)*(1.0+Eta);
    Ni[3] = 0.25*(1.0-Xi)*(1.0+Eta);
  }

  /*!
   * \brief Shape function Jacobian (dNi) evaluated at point Xi,Eta.
   */
  inline static void ShapeFunctionJacobian(su2double Xi, su2double Eta, su2double dNi[][2]) {
    dNi[0][0] = -0.25*(1.0-Eta);  dNi[0][1] = -0.25*(1.0-Xi);
    dNi[1][0] =  0.25*(1.0-Eta);  dNi[1][1] = -0.25*(1.0+Xi);
    dNi[2][0] =  0.25*(1.0+Eta);  dNi[2][1] =  0.25*(1.0+Xi);
    dNi[3][0] = -0.25*(1.0+Eta);  dNi[3][1] =  0.25*(1.0-Xi);
  }

  /*!
   * \brief Compute the value of the area of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Area of the element.
   */
  su2double ComputeArea(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CTETRA1
 * \brief Tetrahedral element with 1 Gauss Point
 * \author R. Sanchez
 */
class CTETRA1 final : public CElementWithKnownSizes<1,4,3> {
private:
  enum : unsigned short {NGAUSS = 1};
  enum : unsigned short {NNODE = 4};
  enum : unsigned short {NDIM = 3};

public:
  /*!
   * \brief Constructor of the class.
   */
  CTETRA1();

  /*!
   * \brief Compute the value of the volume of the element.
   * \return Volume of the element.
   */
  su2double ComputeVolume(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CHEXA8
 * \brief Hexahedral element with 8 Gauss Points
 * \author R. Sanchez
 */
class CHEXA8 final : public CElementWithKnownSizes<8,8,3> {
private:
  enum : unsigned short {NGAUSS = 8};
  enum : unsigned short {NNODE = 8};
  enum : unsigned short {NDIM = 3};

public:
  /*!
   * \brief Constructor of the class.
   */
  CHEXA8();

  /*!
   * \brief Compute the value of the volume of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Volume of the element.
   */
  su2double ComputeVolume(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CPYRAM5
 * \brief Pyramid element with 5 Gauss Points
 * \author R. Sanchez, F. Palacios, A. Bueno, T. Economon, S. Padron.
 */
class CPYRAM5 final : public CElementWithKnownSizes<5,5,3> {
private:
  enum : unsigned short {NGAUSS = 5};
  enum : unsigned short {NNODE = 5};
  enum : unsigned short {NDIM = 3};

public:
  /*!
   * \brief Constructor of the class.
   */
  CPYRAM5();

  /*!
   * \brief Compute the value of the volume of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Volume of the element.
   */
  su2double ComputeVolume(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CPRISM6
 * \brief Prism element with 6 Gauss Points
 * \author R. Sanchez, F. Palacios, A. Bueno, T. Economon, S. Padron.
 * \version 7.1.1 "Blackbird"
 */
class CPRISM6 final : public CElementWithKnownSizes<6,6,3> {
private:
  enum : unsigned short {NGAUSS = 6};
  enum : unsigned short {NNODE = 6};
  enum : unsigned short {NDIM = 3};

public:
  /*!
   * \brief Constructor of the class.
   */
  CPRISM6();

  /*!
   * \brief Compute the value of the volume of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Volume of the element.
   */
  su2double ComputeVolume(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CTRIA3
 * \brief Tria element with 3 Gauss Points
 * \author T.Dick
 */
class CTRIA3 final : public CElementWithKnownSizes<3,3,2> {
private:
  enum : unsigned short {NGAUSS = 3};
  enum : unsigned short {NNODE = 3};
  enum : unsigned short {NDIM = 2};

public:
  /*!
   * \brief Constructor of the class.
   */
  CTRIA3();

  /*!
   * \brief Compute the value of the area of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Area of the element.
   */
  su2double ComputeArea(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CTETRA4
 * \brief Tetrahedral element with 4 Gauss Points
 * \author T.Dick
 */
class CTETRA4 final : public CElementWithKnownSizes<4,4,3> {
private:
  enum : unsigned short {NGAUSS = 4};
  enum : unsigned short {NNODE = 4};
  enum : unsigned short {NDIM = 3};

public:
  /*!
   * \brief Constructor of the class.
   */
  CTETRA4();

  /*!
   * \brief Compute the value of the area of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Area of the element.
   */
  su2double ComputeVolume(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CPYRAM6
 * \brief Pyramid element with 6 Gauss Points
 * \author T.Dick
 */
class CPYRAM6 final : public CElementWithKnownSizes<6,5,3> {
private:
  enum : unsigned short {NGAUSS = 6};
  enum : unsigned short {NNODE = 5};
  enum : unsigned short {NDIM = 3};

public:
  /*!
   * \brief Constructor of the class.
   */
  CPYRAM6();

  /*!
   * \brief Compute the value of the volume of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Volume of the element.
   */
  su2double ComputeVolume(const FrameType mode = REFERENCE) const override;

};

/*!
 * \class CLINE
 * \brief line element with 2 Gauss Points
 * \author T.Dick
 */
class CLINE final : public CElementWithKnownSizes<2,2,1> {
private:
  enum : unsigned short {NGAUSS = 2};
  enum : unsigned short {NNODE = 2};
  enum : unsigned short {NDIM = 1};

public:
  /*!
   * \brief Constructor of the class.
   */
  CLINE();

  /*!
   * \brief Compute the value of the area of the element.
   * \param[in] mode - Type of coordinates to consider in the computation.
   * \return Area of the element.
   */
  su2double ComputeLength(const FrameType mode = REFERENCE) const override;

};
