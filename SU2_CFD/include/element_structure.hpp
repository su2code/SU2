/*!
 * \file element_structure.hpp
 * \brief Headers of the finite element structure (elements)
 *        The subroutines and functions are in the <i>element_structure.cpp</i> file.
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#ifdef HAVE_MPI
  #include "mpi.h"
#endif
#include <cmath>
#include <iostream>
#include <cstdlib>

#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "gauss_structure.hpp"

using namespace std;

/*!
 * \class CElement
 * \brief Main class for defining the element structure.
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
 */

class CElement {
protected:
	unsigned short nGaussPoints;		/*!< \brief Number of gaussian points. */
	unsigned short nGaussPointsP;		/*!< \brief Number of gaussian points for the pressure term. */
	unsigned short nNodes;				/*!< \brief Number of gaussian points. */
	static unsigned short nDim;		/*!< \brief Number of dimension of the problem. */
	double **CurrentCoord,				/*!< \brief Coordinates in the current frame. */
	**RefCoord;							/*!< \brief Coordinates in the reference frame. */
	double **GaussCoord,				/*!< \brief Parent coordinates of the Gaussian Points. */
	*GaussWeight;						/*!< \brief Weight of the Gaussian Points for the integration. */
	double	**GaussCoordP,				/*!< \brief Parent coordinates of the Gaussian Points for the pressure subintegration.. */
	*GaussWeightP;						/*!< \brief Weight of the Gaussian Points for the pressure subintegration. */
	CGaussVariable **GaussPoint;		/*!< \brief Structure for the Gaussian Points. */
	CGaussVariable **GaussPointP;		/*!< \brief Structure for the Gaussian Points for the pressure subintegration. */
	double ***Kab;						/*!< \brief Structure for the constitutive component of the tangent matrix. */
	double **Ks_ab;						/*!< \brief Structure for the stress component of the tangent matrix. */
	double ***Kk_ab;					/*!< \brief Structure for the pressure component of the tangent matrix. */
	double **Kt_a;						/*!< \brief Structure for the nodal stress term for the residual computation. */

public:
	/*!
	 * \brief Constructor of the class.
	 */
	CElement(void);

	/*!
	 * \overload
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
	unsigned short GetnNodes(void);

	/*!
	 * \brief Retrieve the number of nodes of the element.
	 * \param[out] nGaussPoints - Number of Gaussian Points of the element.
	 */
	unsigned short GetnGaussPoints(void);

	/*!
	 * \brief Retrieve the number of nodes of the element.
	 * \param[out] nGaussPointsP - Number of Gaussian Points for pressure underintegration.
	 */
	unsigned short GetnGaussPointsP(void);

	/*!
	 * \brief Set the value of the coordinate of the nodes in the reference configuration.
	 * \param[in] val_CoordRef - Value of the coordinate.
	 * \param[in] iNode - Number of node.
	 * \param[in] iDim - Dimension
	 */
	void SetRef_Coord(double val_CoordRef, unsigned short iNode, unsigned short iDim);

	/*!
	 * \brief Set the value of the coordinate of the nodes in the current configuration.
	 * \param[in] val_CoordRef - Value of the coordinate.
	 * \param[in] iNode - Number of node.
	 * \param[in] iDim - Dimension
	 */
	void SetCurr_Coord(double val_CoordCurr, unsigned short iNode, unsigned short iDim);

	/*!
	 * \brief Set the value of the coordinate of the nodes in the reference configuration.
	 * \param[in] val_CoordRef - Value of the coordinate.
	 * \param[in] iNode - Number of node.
	 * \param[in] iDim - Dimension
	 * \param[out] Coordinate
	 */
	double GetRef_Coord(unsigned short iNode, unsigned short iDim);

	/*!
	 * \brief Get the value of the coordinate of the nodes in the current configuration.
	 * \param[in] val_CoordRef - Value of the coordinate.
	 * \param[in] iNode - Number of node.
	 * \param[in] iDim - Dimension
	 * \param[out] Coordinate
	 */
	double GetCurr_Coord(unsigned short iNode, unsigned short iDim);

	/*!
	 * \brief Get the weight of the corresponding Gaussian Point.
	 * \param[in] iGauss - index of the Gaussian point.
	 * \param[out] Weight.
	 */
	double GetWeight(unsigned short iGauss);

	/*!
	 * \brief Get the weight of the corresponding Gaussian Point for pressure subintegration.
	 * \param[in] iGaussP - index of the Gaussian point.
	 * \param[out] Weight.
	 */
	double GetWeight_P(unsigned short iGaussP);

	/*!
	 * \brief Get the jacobian respect to the reference configuration for the Gaussian Point iGauss.
	 * \param[in] iGauss - index of the Gaussian point.
	 * \param[out] Weight.
	 */
	double GetJ_X(unsigned short iGauss);

	/*!
	 * \brief Get the jacobian respect to the current configuration for the Gaussian Point iGauss.
	 * \param[in] iGauss - index of the Gaussian point.
	 * \param[out] Weight.
	 */
	double GetJ_x(unsigned short iGauss);

	/*!
	 * \brief Get the jacobian respect to the reference configuration for the Gaussian Point iGauss and the pressure term.
	 * \param[in] iGauss - index of the Gaussian point.
	 * \param[out] Weight.
	 */
	double GetJ_X_P(unsigned short iGauss);

	/*!
	 * \brief Get the jacobian respect to the current configuration for the Gaussian Point iGauss and the pressure term.
	 * \param[in] iGauss - index of the Gaussian point.
	 * \param[out] Weight.
	 */
	double GetJ_x_P(unsigned short iGauss);

	/*!
	 * \brief Add the value of a submatrix K relating nodes a and b, for the constitutive term.
	 * \param[in] nodeA - index of Node a.
	 * \param[in] nodeB - index of Node b.
	 * \param[in] val_Kab - value of the matrix K.
	 */
	void Add_Kab(double **val_Kab, unsigned short nodeA, unsigned short nodeB);

	/*!
	 * \brief Add the value of a submatrix K relating nodes a and b, for the constitutive term (symmetric terms need transpose)
	 * \param[in] nodeA - index of Node a.
	 * \param[in] nodeB - index of Node b.
	 * \param[in] val_Kab - value of the matrix K.
	 */
	void Add_Kab_T(double **val_Kab, unsigned short nodeA, unsigned short nodeB);


	/*!
	 * \brief Add the value of the diagonal term for the stress contribution to the stiffness of the system.
	 * \param[in] nodeA - index of Node a.
	 * \param[in] nodeB - index of Node b.
	 * \param[in] val_Ks_ab - value of the term that will constitute the diagonal of the stress contribution.
	 */
	void Add_Ks_ab(double val_Ks_ab, unsigned short nodeA, unsigned short nodeB);

	/*!
	 * \brief Add the value of the nodal stress term for the computation of the residual.
	 * \param[in] nodeA - index of Node a.
	 * \param[in] val_Kt_a - value of the term that will constitute the diagonal of the stress contribution.
	 */
	void Add_Kt_a(double *val_Kt_a, unsigned short nodeA);

	/*!
	 * \brief Set the value of a submatrix K relating nodes a and b, for the pressure term (this term is subintegrated).
	 * \param[in] nodeA - index of Node a.
	 * \param[in] nodeB - index of Node b.
	 * \param[in] val_Kab - value of the matrix K.
	 */
	void Set_Kk_ab(double **val_Kk_ab, unsigned short nodeA, unsigned short nodeB);

	/*!
	 * \brief Restarts the values in the element.
	 */
	void clearElement(void);

	/*!
	 * \brief Return the value of the submatrix K relating nodes a and b.
	 * \param[in] nodeA - index of Node a.
	 * \param[in] nodeB - index of Node b.
	 * \param[out] val_Kab - value of the matrix K.
	 */
	double *Get_Kab(unsigned short nodeA, unsigned short nodeB);

	/*!
	 * \brief Return the value of the diagonal term for the stress contribution, relating nodes a and b.
	 * \param[in] nodeA - index of Node a.
	 * \param[in] nodeB - index of Node b.
	 * \param[out] val_Kab - value of the matrix K.
	 */
	double Get_Ks_ab(unsigned short nodeA, unsigned short nodeB);

	/*!
	 * \brief Return the value of a submatrix K relating nodes a and b, for the pressure term (this term is subintegrated).
	 * \param[in] nodeA - index of Node a.
	 * \param[in] nodeB - index of Node b.
	 * \param[in] val_Kab - value of the matrix K.
	 */
	double *Get_Kk_ab(unsigned short nodeA, unsigned short nodeB);

	/*!
	 * \brief Return the value of the nodal stress component of the residual for node a.
	 * \param[in] nodeA - index of Node a.
	 * \param[out] val_Kt_a - value of the stress term.
	 */
	double *Get_Kt_a(unsigned short nodeA);

	/*!
	 * \brief Retrieve the value of the gradient of the shape functions respect to the reference configuration.
	 * \param[in] iNode - Index of the node.
	 * \param[in] iNode - Index of the Gaussian Point.
	 * \param[out] GradNi_X - Gradient of the shape function related to node iNode and evaluated at Gaussian Point iGauss
	 */
	double GetGradNi_X(unsigned short iNode, unsigned short iGauss, unsigned short iDim);

	/*!
	 * \brief Retrieve the value of the gradient of the shape functions respect to the reference configuration.
	 * \param[in] iNode - Index of the node.
	 * \param[in] iNode - Index of the Gaussian Point.
	 * \param[out] GradNi_X - Gradient of the shape function related to node iNode and evaluated at Gaussian Point iGauss
	 */
	double GetGradNi_x(unsigned short iNode, unsigned short iGauss, unsigned short iDim);

	/*!
	 * \brief Retrieve the value of the gradient of the shape functions respect to the reference configuration.
	 * \param[in] iNode - Index of the node.
	 * \param[in] iNode - Index of the Gaussian Point.
	 * \param[out] GradNi_x - Gradient of the shape function related to node iNode and evaluated at Gaussian Point iGauss
	 */
	double GetGradNi_x_P(unsigned short iNode, unsigned short iGaussP, unsigned short iDim);

	/*!
	 * \brief Set the value of the gradient of the shape functions respect to the reference configuration.
	 * \param[in] val_solution - Solution of the problem.
	 * \param[out] J_X - Jacobian of the element evaluated at the current Gauss Point respect to the reference configuration
	 */
	virtual void ComputeGrad_Linear(void);

	/*!
	 * \brief Set the value of the gradient of the shape functions respect to the current configuration.
	 * \param[in] val_solution - Solution of the problem.
	 * \param[out] J_x - Jacobian of the element evaluated at the current Gauss Point respect to the current configuration
	 */
	virtual void ComputeGrad_NonLinear(void);

	/*!
	 * \brief Virtual member
	 */
	virtual void ComputeGrad_Pressure(void);


};

/*!
 * \class CTRIA1
 * \brief Tria element with 1 Gauss Points
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
 */

class CTRIA1 : public CElement {

protected:

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

};


/*!
 * \class CQUAD4
 * \brief Quadrilateral element with 4 Gauss Points
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
 */

class CQUAD4 : public CElement {

protected:

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
	 * \brief Virtual member.
	 */
	virtual void ComputeGrad_Pressure(void);


};

/*!
 * \class CQUAD4P1
 * \brief Quadrilateral element with 4 Gauss Points and 1 Gauss Point for pressure subintegration
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
 */

class CQUAD4P1 : public CQUAD4 {

protected:

public:

	/*!
	 * \brief Constructor of the class.
	 */
	CQUAD4P1(void);

	/*!
	 * \overload
	 * \param[in] val_fea - Values of the fea solution (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CQUAD4P1(unsigned short val_nDim, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CQUAD4P1(void);

	/*!
	 * \brief Set the value of the gradient of the shape functions respect to the current configuration on 1 Gauss Point.
	 * \param[in] val_solution - Solution of the problem.
	 * \param[out] J_X - Jacobian of the element evaluated at the current Gauss Point respect to the reference configuration
	 */
	void ComputeGrad_Pressure(void);


};

/*!
 * \class CTETRA1
 * \brief Tetrahedral element with 1 Gauss Point
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
 */

class CTETRA1 : public CElement {

protected:

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

};

/*!
 * \class CHEXA8
 * \brief Hexahedral element with 8 Gauss Points
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
 */

class CHEXA8 : public CElement {

protected:

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
	 * \brief Virtual member.
	 */
	virtual void ComputeGrad_Pressure(void);


};

/*!
 * \class CHEXA8P1
 * \brief Hexahedral element with 8 Gauss Points and 1 Gauss Point for pressure subintegration
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
 */

class CHEXA8P1 : public CHEXA8 {

protected:

public:

	/*!
	 * \brief Constructor of the class.
	 */
	CHEXA8P1(void);

	/*!
	 * \overload
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CHEXA8P1(unsigned short val_nDim, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CHEXA8P1(void);

	/*!
	 * \brief Set the value of the gradient of the shape functions respect to the current configuration on 1 Gauss Point.
	 * \param[in] val_solution - Solution of the problem.
	 * \param[out] J_X - Jacobian of the element evaluated at the current Gauss Point respect to the reference configuration
	 */
	void ComputeGrad_Pressure(void);


};

#include "element_structure.inl"
