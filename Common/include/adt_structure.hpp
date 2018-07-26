/*!
 * \file adt_structure.hpp
 * \brief Headers of the subroutines for carrying out geometrical searches using an
 *        alternating digital tree (ADT).
 *        The subroutines and functions are in the <i>adt_structure.cpp</i> file.
 * \author E. van der Weide
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>

#include "./mpi_structure.hpp"
#include "./option_structure.hpp"

using namespace std;

/*! 
 * \class su2_adtComparePointClass
 * \brief  Functor, used for the sorting of the points when building an ADT.
 * \author E. van der Weide
 * \version 4.1.3 "Cardinal"
 */
class su2_adtComparePointClass {
private:
  const su2double      *pointCoor;      /*!< \brief Pointer to the coordinates of the points. */
  const unsigned short splitDirection;  /*!< \brief Split direction used in the sorting. */
  const unsigned short nDim;            /*!< \brief Number of spatial dimensions stored in the coordinates. */

public:
  /*!
   * \brief Constructor of the class. The member variables are initialized.
   * \param[in] coor      Pointer to the coordinates of the points.
   * \param[in] splitDir  Direction that must be used to sort the coordinates.
   * \param[in] nDimADT   Number of spatial dimensions of the ADT and coordinates.
   */
  su2_adtComparePointClass(const su2double      *coor,
                           const unsigned short splitDir,
                           const unsigned short nDimADT);
  /*!
   * \brief Destructor, nothing to be done.
   */
  ~su2_adtComparePointClass();

  /*!
   * \brief Operator used for the sorting of the points.
   * \param[in] p0  Index of the first point to be compared.
   * \param[in] p1  Index of the second point to be compared.
   */
  bool operator()(const unsigned long p0,
                  const unsigned long p1) const;
private:
  /*!
   * \brief Default constructor of the class, disabled.
   */
  su2_adtComparePointClass();
};

/*! 
 * \class su2_adtNodeClass
 * \brief  Class for storing the information needed in a node of an ADT.
 * \author E. van der Weide
 * \version 4.1.3 "Cardinal"
 */
class su2_adtNodeClass {
public:
  bool          childrenAreTerminal[2];  /*!< \brief Whether or not the child leaves are terminal. */
  unsigned long children[2];             /*!< \brief Child leaves. If childrenAreTerminal is true the children
                                                     contain the point ID's or bounding box ID's. Note that it
                                                     is allowed that one child is termimal and the other is not. */
  unsigned long centralNodeID;           /*!< \brief ID of a node, which is near the center of the leaf. */

  su2double *xMin;  /*!< \brief The minimum coordinates of this leaf. It points to a position in the large
                                vector, which contains the coordinates of all leaves. */
  su2double *xMax;  /*!< \brief The maximum coordinates of this leaf. It points to a position in the large
                                vector, which contains the coordinates of all leaves. */

  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  su2_adtNodeClass();

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~su2_adtNodeClass();

  /*!
   * \brief Copy constructor of the class.
   * \param[in] other  Object from which the data must be copied.
   */
  su2_adtNodeClass(const su2_adtNodeClass &other);

  /*!
   * \brief Assignment operator.
   * \param[in] other  Object from which the data must be copied.
   */
  su2_adtNodeClass& operator=(const su2_adtNodeClass &other);

private:

  /*!
   * \brief Copy function, which copies the data from the given object.
   * \param[in] other  Object from which the data must be copied.
   */
  void Copy(const su2_adtNodeClass &other);
};

/*! 
 * \class su2_adtBaseClass
 * \brief  Base class for storing an ADT in an arbitrary number of dimensions.
 * \author E. van der Weide
 * \version 4.1.3 "Cardinal"
 */
class su2_adtBaseClass {
protected:
  unsigned long nLeaves;    /*!< \brief Number of leaves in the ADT. */
  unsigned short nDimADT;   /*!< \brief Number of dimensions of the ADT. */
  bool           isEmpty;   /*!< \brief Whether or not the ADT is empty. */

  vector<su2_adtNodeClass> leaves; /*!< \brief Vector, which contains all the leaves of the ADT. */

private:
  vector<su2double> coorMinLeaves; /*!< \brief Vector, which contains all the minimum coordinates
                                               of the leaves. */
  vector<su2double> coorMaxLeaves; /*!< \brief Vector, which contains all the maximum coordinates
                                               of the leaves. */
protected:
  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  su2_adtBaseClass();

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  virtual ~su2_adtBaseClass();  

  /*!
   * \brief Function, which builds the ADT of the given coordinates.
   * \param[in] nDim    Number of dimensions of the ADT.
   * \param[in] nPoints Number of points present in the ADT.
   * \param[in] coor    Coordinates of the points.
   */
  void BuildADT(unsigned short  nDim,
                unsigned long   nPoints,
                const su2double *coor);
public:
  /*!
   * \brief Function, which returns whether or not the ADT is empty.
   * \return  Whether or not the ADT is empty.
   *
   */
  bool IsEmpty(void) const;

private:
  /*!
   * \brief Copy constructor of the class, disabled.
   */
  su2_adtBaseClass(const su2_adtBaseClass &);

  /*!
   * \brief Assignment operator, disabled.
   */
  su2_adtBaseClass& operator=(const su2_adtBaseClass &);
};

/*! 
 * \class su2_adtPointsOnlyClass
 * \brief  Class for storing an ADT of only points in an arbitrary number of dimensions.
 * \author E. van der Weide
 * \version 4.1.3 "Cardinal"
 */
class su2_adtPointsOnlyClass : public su2_adtBaseClass {
private:
  vector<unsigned long> frontLeaves;    /*!< \brief Vector used in the tree traversal. */
  vector<unsigned long> frontLeavesNew; /*!< \brief Vector used in the tree traversal. */

  vector<su2double>     coorPoints;    /*!< \brief Vector, which contains the coordinates
                                                   of the points in the ADT. */
  vector<unsigned long> localPointIDs; /*!< \brief Vector, which contains the local point ID's
                                                   of the points in the ADT. */
  vector<int>           ranksOfPoints; /*!< \brief Vector, which contains the ranks
                                                   of the points in the ADT. */
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] nDim    Number of spatial dimensions of the problem.
   * \param[in] nPoints Number of local points to be stored in the ADT.
   * \param[in] coor    Coordinates of the local points.
   * \param[in] pointID Local point IDs of the local points. 
   */
  su2_adtPointsOnlyClass(unsigned short      nDim,
                         unsigned long       nPoints,
                         const su2double     *coor,
                         const unsigned long *pointID);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~su2_adtPointsOnlyClass();

  /*!
   * \brief Function, which determines the nearest node in the ADT for the
            given coordinate.
   * \param[in]  coor    Coordinate for which the nearest node in the ADT must be determined.
   * \param[out] dist    Distance to the nearest node in the ADT.
   * \param[out] pointID Local point ID of the nearest node in the ADT.
   * \param[out] rankID  Rank on which the nearest node in the ADT is stored.
   */
  void DetermineNearestNode(const su2double *coor,
                            su2double       &dist,
                            unsigned long   &pointID,
                            int             &rankID);

  void Determine_N_NearestNodes(int N, const su2double *coor,
                              su2double *       dist,
                              unsigned long*  pointID,
                              int*             rankID);


private:
  /*!
   * \brief Default constructor of the class, disabled.
   */
  su2_adtPointsOnlyClass();

  /*!
   * \brief Copy constructor of the class, disabled.
   */
  su2_adtPointsOnlyClass(const su2_adtPointsOnlyClass &);

  /*!
   * \brief Assignment operator, disabled.
   */
  su2_adtPointsOnlyClass& operator=(const su2_adtPointsOnlyClass &);
};

#include "adt_structure.inl"
