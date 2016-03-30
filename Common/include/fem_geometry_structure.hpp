/*!
 * \file fem_geometry_structure.hpp
 * \brief Headers of the main subroutines for creating the geometrical structure for the FEM solver.
 *        The subroutines and functions are in the <i>fem_geometry_structure.cpp</i> file.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
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

#include "./geometry_structure.hpp"

using namespace std;

/*!
 * \class CVolumeElementFEM
 * \brief Class to store a volume element for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CVolumeElementFEM {
public:
  bool elemIsOwned;             /*!< \brief Whether or not this is an owned element. */
  bool JacIsConsideredConstant; /*!< \brief Whether or not the Jacobian of the transformation
                                     to the standard element is considered constant. */

  int rankOriginal;            /*!< \brief The rank where the original volume is stored. For
                                    the original volume, this is simply the current rank. */

  short periodIndexToDonor;    /*!< \brief The index of the periodic transformation to the donor
                                    element. Only for halo elements. A -1 indicates no periodic
                                    transformation. */

  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;    /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nPolySol;     /*!< \brief Polynomial degree for the solution of the element. */
  unsigned short nDOFsGrid;    /*!< \brief Number of DOFs for the geometry of the element. */
  unsigned short nDOFsSol;     /*!< \brief Number of DOFs for the solution of the element. */
  unsigned short nFaces;       /*!< \brief Number of faces of the element. */

  unsigned short indStandardElement; /*!< \brief Index in the vector of standard elements. */

  unsigned long elemIDGlobal;        /*!< \brief Global element ID of this element. */
  unsigned long offsetDOFsSolGlobal; /*!< \brief Global offset of the solution DOFs of this element. */
  unsigned long offsetDOFsSolLocal;  /*!< \brief Local offset of the solution DOFs of this element. */

  vector<long> neighborElemIDMatchingFaces;  /*!< \brief Vector with the neighbor element IDs of the faces of this
                                                   element. A negative value indicates no matching neighbor. */
  vector<short> periodIndexNeighbors;        /*!< \brief Vector with the periodic indices of the neighbors to the
                                                  faces. A -1 indicates no periodic transformation to the neighbor. */
  vector<bool> JacFacesIsConsideredConstant; /*!< \brief Vector with the booleans whether the Jacobian of the
                                                  transformation to the standard element is constant for the faces. */

  vector<unsigned long> nodeIDsGrid; /*!< \brief Vector with the node IDs of the grid for this element. */

  /*!
   * \brief Constructor of the class. Nothing to be done
   */
  CVolumeElementFEM(void);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CVolumeElementFEM(void);
};

/*!
 * \class CPointFEM
 * \brief Class to a point for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CPointFEM {
public:
  unsigned long globalID;    /*!< \brief The global ID of this point in the grid. */
  short periodIndexToDonor;  /*!< \brief The index of the periodic transformation to the donor
                                  element. Only for halo elements. A -1 indicates no periodic
                                  transformation. */
  su2double coor[3];         /*!< \brief Array with the coordinates of the node. */

  /*!
   * \brief Default constructor of the class. Initialize the coordinates to zero
            to avoid a valgrind warning in two space dimensions.
   */
  CPointFEM(void);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CPointFEM(void);

  /*!
   * \brief Copy constructor of the class.
   */
  CPointFEM(const CPointFEM &other);

  /*!
   * \brief Assignment operator of the class.
   */
  CPointFEM& operator=(const CPointFEM &other);

  /*!
   * \brief Less than operator of the class. Needed for the sorting.
   */
  bool operator<(const CPointFEM &other) const;

  /*!
   * \brief Equal operator of the class. Needed for the removal of double entities.
   */
  bool operator==(const CPointFEM &other) const;

private:
  /*!
   * \brief Copy function. Needed for the copy constructor and assignment operator.
   */
  void Copy(const CPointFEM &other);
};

/*!
 * \class CSurfaceElementFEM
 * \brief Class to store a surface element for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CSurfaceElementFEM {
public:
  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;    /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nDOFsGrid;    /*!< \brief Number of DOFs for the geometry of the element. */

  unsigned short indStandardElement; /*!< \brief Index in the vector of standard elements. */

  unsigned long volElemID;         /*!< \brief ID of the corresponding volume element. */
  unsigned long boundElemIDGlobal; /*!< \brief Global ID of this surface element inside
                                        the boundary to which it belongs. */

  vector<unsigned long> nodeIDsGrid; /*!< \brief Vector with the node IDs of the grid for this element. */

  /*!
   * \brief Constructor of the class. Initialize some variables to avoid a valgrid warning.
   */
  CSurfaceElementFEM(void);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CSurfaceElementFEM(void);

  /*!
   * \brief Copy constructor of the class.
   */
  CSurfaceElementFEM(const CSurfaceElementFEM &other);

  /*!
   * \brief Assignment operator of the class.
   */
  CSurfaceElementFEM& operator=(const CSurfaceElementFEM &other);

  /*!
   * \brief Less than operator of the class. Needed for the sorting.
   */
  bool operator<(const CSurfaceElementFEM &other) const;

private:
  /*!
   * \brief Copy function. Needed for the copy constructor and assignment operator.
   */
  void Copy(const CSurfaceElementFEM &other);
};

/*!
 * \class CBoundaryFEM
 * \brief Class to store a boundary for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CBoundaryFEM {
public:
  string markerTag;  /*!< \brief Marker tag of this boundary. */

  vector<CSurfaceElementFEM> surfElem; /*!< \brief Vector of the local surface elements. */

  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  CBoundaryFEM(void);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CBoundaryFEM(void);
};

/*!
 * \class CMeshFEM
 * \brief Base class for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CMeshFEM: public CGeometry {
protected:
  unsigned long nVolElemOwned;  /*!< \brief Number of local volume elements owned by this rank. */
  unsigned long nVolElemTot;    /*!< \brief Total number of local volume elements, including halos. */

  vector<CVolumeElementFEM> volElem; /*!< \brief Vector of the local volume elements, including halos. */

  vector<CPointFEM> meshPoints;      /*!< \brief Vector of the points of the FEM mesh. */

  vector<CBoundaryFEM> boundaries;   /*!< \brief Vector of the boundaries of the FEM mesh. */
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CMeshFEM(void);

  /*!
	 * \overload
	 * \brief Redistributes the grid over the ranks and creates the halo layer.
  * \param[in] geometry - The linear distributed grid that must be redistributed.
	 * \param[in] config   - Definition of the particular problem.
	 */
  CMeshFEM(CGeometry *geometry, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CMeshFEM(void);
};

/*!
 * \class CMeshFEM_DG
 * \brief Class which contains all the variables for the DG FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CMeshFEM_DG: public CMeshFEM {

public:

 /*!
  * \brief Constructor of the class.
  */
 CMeshFEM_DG(void);

  /*!
  * \overload
  * \brief Redistributes the grid over the ranks and creates the halo layer.
  * \param[in] geometry - The linear distributed grid that must be redistributed.
  * \param[in] config   - Definition of the particular problem.
  */
  CMeshFEM_DG(CGeometry *geometry, CConfig *config);
 
 /*!
  * \brief Destructor of the class.
  */
 ~CMeshFEM_DG(void);
};

#include "fem_geometry_structure.inl"
