/*!
 * \file transfer_structure.hpp
 * \brief Headers of the transfer structure
 *        The subroutines and functions are in the <i>transfer_structure.cpp</i> and <i>transfer_physics.cpp</i> files.
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

#include "../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/vector_structure.hpp"
#include "solver_structure.hpp"

using namespace std;

/*!
 * \class CTransfer
 * \brief Main class for defining the physical transfer of information.
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 */

class CTransfer {
protected:

  int rank, 	/*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */

  su2double *Physical_Constants;
  su2double *Donor_Variable;
  su2double *Target_Variable;
  bool valAggregated;

  /*--- Mixing Plane interface variable ---*/
  su2double 	  *SpanValueCoeffTarget;
  unsigned short *SpanLevelDonor;
  unsigned short nSpanMaxAllZones;




  unsigned short nVar;

public:
  /*!
   * \brief Constructor of the class.
   */
  CTransfer(void);

  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] val_nConst - Number of physical constants that need to be taken into account.
   * \param[in] config - Definition of the particular problem.
   */
  CTransfer(unsigned short val_nVar, unsigned short val_nConst, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CTransfer(void);

  /*!
   * \brief Interpolate data and broadcast it into all processors, for nonmatching meshes.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void Broadcast_InterfaceData(CSolver *donor_solution, CSolver *target_solution,
                               CGeometry *donor_geometry, CGeometry *target_geometry,
                               CConfig *donor_config, CConfig *target_config);
  /*!
   * \brief A virtual member.
   */

  virtual void GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                                       CGeometry *donor_geometry, CGeometry *target_geometry,
                     CConfig *donor_config, CConfig *target_config);
  /*!
   * \brief A virtual member.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] Marker_Donor - Index of the donor marker.
   * \param[in] Vertex_Donor - Index of the donor vertex.
   */
  virtual void GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry,
                   CConfig *donor_config, unsigned long Marker_Donor,
                   unsigned long Vertex_Donor, unsigned long Point_Donor);

  /*!
   * \brief Initializes the target variable.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] nDonorPoints - Number of donor points.
   */
  virtual void InitializeTarget_Variable(CSolver *target_solution, unsigned long Marker_Target,
                                         unsigned long Vertex_Target, unsigned short nDonorPoints);

  /*!
   * \brief Recovers the target variable from the buffer of su2doubles that was broadcasted.
   * \param[in] indexPoint_iVertex - index of the vertex in the buffer array.
   * \param[in] Buffer_Bcast_Variables - full broadcasted buffer array of doubles.
   * \param[in] donorCoeff - value of the donor coefficient.
   */
  virtual void RecoverTarget_Variable(long indexPoint_iVertex, su2double *Buffer_Bcast_Variables,
                                      su2double donorCoeff);

  /*!
   * \brief A virtual member.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] Point_Target - Index of the target point.
   */
  virtual void SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                  CConfig *target_config, unsigned long Marker_Target,
                  unsigned long Vertex_Target, unsigned long Point_Target);

  /*!
   * \brief A virtual member.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_zone      - Index of the donorZone.
   */

  virtual void SetAverageValues(CSolver *donor_solution, CSolver *target_solution,  unsigned short donorZone);

  /*!
   * \brief A virtual member.
   * \param[in] donor_geometry - Geometry of the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_zone      - Index of the donorZone.
   */
  virtual void SetAverageTurboGeoValues(CGeometry *donor_geometry, CGeometry *target_geometry, unsigned short donorZone);

  /*!
   * \brief A virtual member.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  virtual void SetSpanWiseLevels(CConfig *donor_config, CConfig *target_config);


  /*!
   * \brief Transfer pre-processing for the mixing plane inteface.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void Preprocessing_InterfaceAverage(CGeometry *donor_geometry, CGeometry *target_geometry,
                                      CConfig *donor_config, CConfig *target_config, unsigned short iMarkerInt);


  /*!
   * \brief Interpolate data and scatter it into different processors, for matching meshes.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void Allgather_InterfaceAverage(CSolver *donor_solution, CSolver *target_solution,
                                  CGeometry *donor_geometry, CGeometry *target_geometry,
                                  CConfig *donor_config, CConfig *target_config, unsigned short iMarkerInt);

  /*!
   * \brief Interpolate data and scatter it into different processors, for matching meshes.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GatherAverageValues(CSolver *donor_solution, CSolver *target_solution, unsigned short donorZone);

  /*!
   * \brief Exchange Average geometrical value beteween zones .
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GatherAverageTurboGeoValues(CGeometry *donor_geometry, CGeometry *target_geometry, unsigned short donorZone);


};

/*!
 * \class CTransfer_FlowTraction
 * \brief Transfer flow tractions from a fluid zone into a structural zone
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 */

class CTransfer_FlowTraction : public CTransfer {

protected:
  bool consistent_interpolation;
  
  /*!
   * \brief Sets the dimensional factor for pressure and the consistent_interpolation flag
   * \param[in] flow_config - Definition of the fluid (donor) problem.
   */
  void Preprocess(CConfig *flow_config);

public:

  /*!
   * \brief Constructor of the class.
   */
  CTransfer_FlowTraction(void);

  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] config - Definition of the particular problem.
   */
  CTransfer_FlowTraction(unsigned short val_nVar, unsigned short val_nConst, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CTransfer_FlowTraction(void);

  /*!
   * \brief Retrieve some constants needed for the calculations.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                 CGeometry *donor_geometry, CGeometry *target_geometry,
                 CConfig *donor_config, CConfig *target_config);

  /*!
   * \brief Retrieve the variable that will be sent from donor mesh to target mesh.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] Marker_Donor - Index of the donor marker.
   * \param[in] Vertex_Donor - Index of the donor vertex.
   */
  void GetDonor_Variable(CSolver *flow_solution, CGeometry *flow_geometry, CConfig *flow_config,
               unsigned long Marker_Flow, unsigned long Vertex_Flow, unsigned long Point_Flow);

  /*!
   * \brief Set the variable that has been received from the target mesh into the target mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] Point_Target - Index of the target point.
   */
  void SetTarget_Variable(CSolver *fea_solution, CGeometry *fea_geometry,
              CConfig *fea_config, unsigned long Marker_Struct,
              unsigned long Vertex_Struct, unsigned long Point_Struct);

};

/*!
 * \class CTransfer_StructuralDisplacements
 * \brief Transfer structural displacements from a structural zone into a fluid zone
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 */

class CTransfer_StructuralDisplacements : public CTransfer {

protected:

public:

  /*!
   * \brief Constructor of the class.
   */
  CTransfer_StructuralDisplacements(void);

  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] config - Definition of the particular problem.
   */
  CTransfer_StructuralDisplacements(unsigned short val_nVar, unsigned short val_nConst, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CTransfer_StructuralDisplacements(void);

  /*!
   * \brief Retrieve some constants needed for the calculations.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                 CGeometry *donor_geometry, CGeometry *target_geometry,
                 CConfig *donor_config, CConfig *target_config);

  /*!
   * \brief Retrieve the variable that will be sent from donor mesh to target mesh.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] Marker_Donor - Index of the donor marker.
   * \param[in] Vertex_Donor - Index of the donor vertex.
   */
  void GetDonor_Variable(CSolver *struct_solution, CGeometry *struct_geometry, CConfig *struct_config,
               unsigned long Marker_Struct, unsigned long Vertex_Struct, unsigned long Point_Struct);

  /*!
   * \brief Set the variable that has been received from the target mesh into the target mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] Point_Target - Index of the target point.
   */
  void SetTarget_Variable(CSolver *flow_solution, CGeometry *flow_geometry,
              CConfig *flow_config, unsigned long Marker_Flow,
              unsigned long Vertex_Flow, unsigned long Point_Flow);

};

/*!
 * \class CTransfer_FlowTraction_DiscAdj
 * \brief Transfer flow tractions from a fluid zone into a structural zone in a discrete adjoint simulation
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 */

class CTransfer_FlowTraction_DiscAdj : public CTransfer_FlowTraction {

protected:

public:

  /*!
   * \brief Constructor of the class.
   */
  CTransfer_FlowTraction_DiscAdj(void);

  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] config - Definition of the particular problem.
   */
  CTransfer_FlowTraction_DiscAdj(unsigned short val_nVar, unsigned short val_nConst, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CTransfer_FlowTraction_DiscAdj(void);

  /*!
   * \brief Retrieve some constants needed for the calculations.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                 CGeometry *donor_geometry, CGeometry *target_geometry,
                 CConfig *donor_config, CConfig *target_config);

};

/*!
 * \class CTransfer_StructuralDisplacements_DiscAdj
 * \brief Transfer structural displacements from a structural zone into a fluid zone in a discrete adjoint simulation
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 */

class CTransfer_StructuralDisplacements_DiscAdj : public CTransfer {

protected:

public:

  /*!
   * \brief Constructor of the class.
   */
  CTransfer_StructuralDisplacements_DiscAdj(void);

  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] config - Definition of the particular problem.
   */
  CTransfer_StructuralDisplacements_DiscAdj(unsigned short val_nVar, unsigned short val_nConst, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CTransfer_StructuralDisplacements_DiscAdj(void);

  /*!
   * \brief Retrieve some constants needed for the calculations.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                 CGeometry *donor_geometry, CGeometry *target_geometry,
                 CConfig *donor_config, CConfig *target_config);

  /*!
   * \brief Retrieve the variable that will be sent from donor mesh to target mesh.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] Marker_Donor - Index of the donor marker.
   * \param[in] Vertex_Donor - Index of the donor vertex.
   */
  void GetDonor_Variable(CSolver *struct_solution, CGeometry *struct_geometry, CConfig *struct_config,
               unsigned long Marker_Struct, unsigned long Vertex_Struct, unsigned long Point_Struct);

  /*!
   * \brief Set the variable that has been received from the target mesh into the target mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] Point_Target - Index of the target point.
   */
  void SetTarget_Variable(CSolver *flow_solution, CGeometry *flow_geometry,
              CConfig *flow_config, unsigned long Marker_Flow,
              unsigned long Vertex_Flow, unsigned long Point_Flow);

};

/*!
 * \class CTransfer_ConservativeVars
 * \brief Transfer conservative variables from a generic zone into another
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 */

class CTransfer_ConservativeVars : public CTransfer {

protected:

public:

  /*!
   * \brief Constructor of the class.
   */
  CTransfer_ConservativeVars(void);

  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] config - Definition of the particular problem.
   */
  CTransfer_ConservativeVars(unsigned short val_nVar, unsigned short val_nConst, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CTransfer_ConservativeVars(void);

  /*!
   * \brief Retrieve some constants needed for the calculations.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                 CGeometry *donor_geometry, CGeometry *target_geometry,
                 CConfig *donor_config, CConfig *target_config);

  /*!
   * \brief Retrieve the variable that will be sent from donor mesh to target mesh.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] Marker_Donor - Index of the donor marker.
   * \param[in] Vertex_Donor - Index of the donor vertex.
   * \param[in] Point_Donor - Index of the donor point.
   */
  void GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry, CConfig *donor_config,
               unsigned long Marker_Donor, unsigned long Vertex_Donor, unsigned long Point_Donor);

  /*!
   * \brief Set the variable that has been received from the target mesh into the target mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] Point_Target - Index of the target point.
   */
  void SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry, CConfig *target_config,
              unsigned long Marker_Target, unsigned long Vertex_Target, unsigned long Point_Target);


};

/*!
 * \class CTransfer_MixingPlaneInterface
 * \brief Transfer average variables needed for MixingPlane computation from a generic zone into another one
 * \author S. Vitale
 * \version 6.2.0 "Falcon"
 */


class CTransfer_MixingPlaneInterface : public CTransfer {

public:

	/*!
	 * \brief Constructor of the class.
	 */
	CTransfer_MixingPlaneInterface(void);

	/*!
	 * \overload
	 * \param[in] val_nVar - Number of variables that need to be transferred.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTransfer_MixingPlaneInterface(unsigned short val_nVar, unsigned short val_nConst, CConfig *donor_config, CConfig *target_config);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CTransfer_MixingPlaneInterface(void);


/*!
 * \brief Initialize quantities for spanwise sections for interpolation.
 * \param[in] donor_config - Definition of the problem at the donor mesh.
 * \param[in] target_config - Definition of the problem at the target mesh.
 */
void SetSpanWiseLevels(CConfig *donor_config, CConfig *target_config);

	/*!
	 * \brief Retrieve the variable that will be sent from donor mesh to target mesh.
	 * \param[in] donor_solution - Solution from the donor mesh.
	 * \param[in] donor_geometry - Geometry of the donor mesh.
	 * \param[in] donor_config - Definition of the problem at the donor mesh.
	 * \param[in] Marker_Donor - Index of the donor marker.
	 * \param[in] Vertex_Donor - Index of the donor vertex.
	 * \param[in] Point_Donor - Index of the donor point.
	 */
	void GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry, CConfig *donor_config,
						   unsigned long Marker_Donor, unsigned long val_Span, unsigned long Point_Donor);

	/*!
	 * \brief Set the variable that has been received from the target mesh into the target mesh.
	 * \param[in] target_solution - Solution from the target mesh.
	 * \param[in] target_geometry - Geometry of the target mesh.
	 * \param[in] target_config - Definition of the problem at the target mesh.
	 * \param[in] Marker_Target - Index of the target marker.
	 * \param[in] Vertex_Target - Index of the target vertex.
	 * \param[in] Point_Target - Index of the target point.
	 */
	void SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry, CConfig *target_config,
							unsigned long Marker_Target, unsigned long val_Span, unsigned long Point_Target);

	/*!
	 * \brief Store all the turboperformance in the solver in ZONE_0.
	 * \param[in] donor_solution  - Solution from the donor mesh.
	 * \param[in] target_solution - Solution from the target mesh.
	 * \param[in] donorZone       - counter of the donor solution
	 */
	void SetAverageValues(CSolver *donor_solution, CSolver *target_solution, unsigned short donorZone);

	/*!
	 * \brief Store all the turboperformance in the solver in ZONE_0.
	 * \param[in] donor_geometry  - Solution from the donor mesh.
	 * \param[in] target_geometry - Solution from the target mesh.
	 * \param[in] donorZone       - counter of the donor solution
	 */
	void SetAverageTurboGeoValues(CGeometry *donor_geometry, CGeometry *target_geometry, unsigned short donorZone);



};


/*!
 * \class CTransfer_SlidingInterface
 * \brief Transfer conservative variables from a generic zone into another
 * \author G. Gori Politecnico di Milano
 * \version 6.2.0 "Falcon"
 */

class CTransfer_SlidingInterface : public CTransfer {

protected:

public:

  /*!
   * \brief Constructor of the class.
   */
  CTransfer_SlidingInterface(void);

  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] config - Definition of the particular problem.
   */
  CTransfer_SlidingInterface(unsigned short val_nVar, unsigned short val_nConst, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CTransfer_SlidingInterface(void);

  /*!
   * \brief Retrieve some constants needed for the calculations.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                 CGeometry *donor_geometry, CGeometry *target_geometry,
                 CConfig *donor_config, CConfig *target_config);

  /*!
   * \brief Retrieve the variable that will be sent from donor mesh to target mesh.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] Marker_Donor - Index of the donor marker.
   * \param[in] Vertex_Donor - Index of the donor vertex.
   * \param[in] Point_Donor - Index of the donor point.
   */
  void GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry, CConfig *donor_config,
               unsigned long Marker_Donor, unsigned long Vertex_Donor, unsigned long Point_Donor);

  /*!
   * \brief A virtual member, initializes the target variable for sliding mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] nDonorPoints - Number of donor points.
   */
  void InitializeTarget_Variable(CSolver *target_solution, unsigned long Marker_Target,
                                 unsigned long Vertex_Target, unsigned short nDonorPoints);

  /*!
   * \brief Recovers the target variable from the buffer of su2doubles that was broadcasted.
   * \param[in] indexPoint_iVertex - index of the vertex in the buffer array.
   * \param[in] Buffer_Bcast_Variables - full broadcasted buffer array of doubles.
   * \param[in] donorCoeff - value of the donor coefficient.
   */
  void RecoverTarget_Variable(long indexPoint_iVertex, su2double *Buffer_Bcast_Variables,
                              su2double donorCoeff);

  /*!
   * \brief Set the variable that has been received from the target mesh into the target mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] Point_Target - Index of the target point.
   */
  void SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry, CConfig *target_config,
              unsigned long Marker_Target, unsigned long Vertex_Target, unsigned long Point_Target);


};

/*!
 * \class CTransfer_ConjugateHeatVars
 * \brief Transfer temperature and heatflux density for conjugate heat interfaces between structure and fluid zones.
 * \author O. Burghardt
 * \version 6.2.0 "Falcon"
 */
class CTransfer_ConjugateHeatVars : public CTransfer {

protected:

public:

  /*!
   * \brief Constructor of the class.
   */
  CTransfer_ConjugateHeatVars(void);

  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] config - Definition of the particular problem.
   */
  CTransfer_ConjugateHeatVars(unsigned short val_nVar, unsigned short val_nConst, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CTransfer_ConjugateHeatVars(void);

  /*!
   * \brief Retrieve the variable that will be sent from donor mesh to target mesh.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] Marker_Donor - Index of the donor marker.
   * \param[in] Vertex_Donor - Index of the donor vertex.
   * \param[in] Point_Donor - Index of the donor point.
   */
  void GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry, CConfig *donor_config,
               unsigned long Marker_Donor, unsigned long Vertex_Donor, unsigned long Point_Donor);

  /*!
   * \brief Set the variable that has been received from the target mesh into the target mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   * \param[in] Marker_Target - Index of the target marker.
   * \param[in] Vertex_Target - Index of the target vertex.
   * \param[in] Point_Target - Index of the target point.
   */
  void SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry, CConfig *target_config,
              unsigned long Marker_Target, unsigned long Vertex_Target, unsigned long Point_Target);
};

#include "transfer_structure.inl"
