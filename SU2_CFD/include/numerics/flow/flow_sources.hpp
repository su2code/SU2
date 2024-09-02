/*!
 * \file flow_sources.hpp
 * \brief Declarations of numerics classes for source-term integration.
 * \author F. Palacios, T. Economon
 * \version 8.0.1 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

#include "../CNumerics.hpp"

/*!
 * \class CSourceBase_Flow
 * \brief Intermediate source term class to allocate the internally
 *        stored residual and Jacobian. Not for stand alone use,
 *        just a helper to build more complicated classes.
 * \ingroup SourceDiscr
 */
class CSourceBase_Flow : public CNumerics {
protected:
  su2double* residual = nullptr;
  su2double** jacobian = nullptr;
  struct StreamwisePeriodicValues SPvals;

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceBase_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

public:
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceBase_Flow() override;

  /*!
   * \brief Set massflow, heatflow & inlet temperature for streamwise periodic flow.
   * \param[in] SolverSPvals - Struct holding the values.
   */
  void SetStreamwisePeriodicValues(const StreamwisePeriodicValues SolverSPvals) final { SPvals = SolverSPvals; }

};

/*!
 * \class CSourceAxisymmetric_Flow
 * \brief Class for source term for solving axisymmetric problems.
 * \ingroup SourceDiscr
 * \author F. Palacios
 */
class CSourceAxisymmetric_Flow : public CSourceBase_Flow {
protected:
    bool implicit, viscous, rans;
    su2double yinv{0.0};

  /*!
  * \brief Diffusion residual of the axisymmetric source term.
  */
  void ResidualDiffusion();

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Residual of the axisymmetric source term.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CSourceGeneralAxisymmetric_Flow
 * \brief Class for source term for solving axisymmetric problems for a general (non ideal) fluid.
 * \ingroup SourceDiscr
 * \author F. Dittmann
 */
class CSourceGeneralAxisymmetric_Flow final : public CSourceAxisymmetric_Flow {
public:

  using CSourceAxisymmetric_Flow::CSourceAxisymmetric_Flow;
  /*!
   * \brief Residual of the general axisymmetric source term.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CSourceIncAxisymmetric_Flow
 * \brief Class for source term for solving incompressible axisymmetric problems.
 * \ingroup SourceDiscr
 * \author T. Economon
 */
class CSourceIncAxisymmetric_Flow final : public CSourceBase_Flow {
  bool implicit, /*!< \brief Implicit calculation. */
  viscous,       /*!< \brief Viscous incompressible flows. */
  energy;        /*!< \brief computation with the energy equation. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceIncAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Residual of the rotational frame source term.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CSourceBodyForce
 * \brief Class for the source term integration of a body force.
 * \ingroup SourceDiscr
 * \author T. Economon
 */
class CSourceBodyForce final : public CSourceBase_Flow {
  su2double Body_Force_Vector[3];

public:
  /*!
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceBodyForce(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Source term integration for a body force.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CSourceIncBodyForce
 * \brief Class for the source term integration of a body force in the incompressible solver.
 * \ingroup SourceDiscr
 * \author T. Economon
 * \version 8.0.1 "Harrier"
 */
class CSourceIncBodyForce final : public CSourceBase_Flow {
  su2double Body_Force_Vector[3];

public:
  /*!
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceIncBodyForce(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Source term integration for a body force.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CSourceBoussinesq
 * \brief Class for the source term integration of the Boussinesq approximation for incompressible flow.
 * \ingroup SourceDiscr
 * \author T. Economon
 * \version 8.0.1 "Harrier"
 */
class CSourceBoussinesq final : public CSourceBase_Flow {
  su2double Gravity_Vector[3];

public:
  /*!
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceBoussinesq(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Source term integration for the Boussinesq approximation.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CSourceGravity
 * \brief Class for the source term integration of the gravity force.
 * \ingroup SourceDiscr
 * \author F. Palacios
 */
class CSourceGravity final : public CSourceBase_Flow {
  su2double Force_Ref;

public:
  /*!
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceGravity(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Source term integration for the poissonal potential.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CSourceRotatingFrame_Flow
 * \brief Class for a rotating frame source term.
 * \ingroup SourceDiscr
 * \author F. Palacios, T. Economon.
 */
class CSourceRotatingFrame_Flow final : public CSourceBase_Flow {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceRotatingFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Residual of the rotational frame source term.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CSourceIncRotatingFrame_Flow
 * \brief Class for a rotating frame source term.
 * \ingroup SourceDiscr
 */
class CSourceIncRotatingFrame_Flow final : public CSourceBase_Flow {

private:
  su2double Omega[3];  /*!< \brief Angular velocity */
  bool implicit; /*!< \brief Implicit calculation. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceIncRotatingFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Residual of the rotational frame source term.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CSourceVorticityConfinement
 * \brief Class for a source term due to vorticity confinement.
 * \ingroup SourceDiscr
 * \brief Vorticity Confinement (VC) technique to counter the numerical
 * diffusion offered by the numerical scheme
 * \author: Y Chandukrishna, Josy P Pullockara, T N Venkatesh.
 * Computational and Theoretical Fluid Dynamics division,
 * CSIR-National Aerospace Laboratories (NAL), Bangalore.
 * Academy of Scientific and Innovative Research (AcSIR), Ghaziabad.
 * First release date : 23 December 2022
 * modified on:
 *
 * VC technique introduces an additional term to the N-S equations that
 * counters the numerical difusion offerd by the numerical schemes. The
 * additional term is introduced as a source term. VC technique requires an
 * input confinement parameter that controls the magnitude of the source
 * term. A suitable Confinement parameter is problem, scheme and grid
 * dependant.
 *
 * Required changes in config file -
 * VORTICITY_CONFINEMENT = YES
 * CONFINEMENT_PARAMETER = <0.0(default)>
 * Currently have support for Viscous, Inviscid cases for both 2D and 3D
 * problems. Can be used along with other source terms also.
 *
 * The current implementation follows,
 * R. Loehner and C. Yang, Vorticity confinement on unstructured grids, 2002.
 * N. Butsuntorn and A. Jameson, Time Spectral Method for Rotorcraft Flow, 2008.
 */
class CSourceVorticityConfinement final : public CSourceBase_Flow {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceVorticityConfinement(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Residual of the rotational frame source term.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;
};

/*!
 * \class CSourceWindGust
 * \brief Class for a source term due to a wind gust.
 * \ingroup SourceDiscr
 * \author S. Padrón
 */
class CSourceWindGust final : public CSourceBase_Flow {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceWindGust(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Residual of the wind gust source term.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CSourceIncStreamwise_Periodic
 * \brief Class for the source term integration of a streamwise periodic body force in the incompressible solver.
 * \ingroup SourceDiscr
 * \author T. Kattmann
 */
class CSourceIncStreamwise_Periodic final : public CSourceBase_Flow {
private:

  bool turbulent; /*!< \brief Turbulence model used. */
  bool energy;    /*!< \brief Energy equation on. */
  bool streamwisePeriodic_temperature; /*!< \brief Periodicity in energy equation */
  su2double Streamwise_Coord_Vector[MAXNDIM] = {0.0}; /*!< \brief Translation vector between streamwise periodic surfaces. */

  su2double norm2_translation, /*!< \brief Square of distance between the 2 periodic surfaces. */
            dot_product, /*!< \brief Container for various dot-products. */
            scalar_factor; /*!< \brief Holds scalar factors to simplify final equations. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceIncStreamwise_Periodic(unsigned short val_nDim,
                                unsigned short val_nVar,
                                CConfig        *config);

  /*!
   * \brief Source term integration for a body force.
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig *config) override;

};

/*!
 * \class CSourceIncStreamwisePeriodic_Outlet
 * \brief Class for the outlet heat sink. Acts like a heatflux boundary on the outlet and not as a volume source.
 * \ingroup SourceDiscr
 * \author T. Kattmann
 */
class CSourceIncStreamwisePeriodic_Outlet : public CSourceBase_Flow {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceIncStreamwisePeriodic_Outlet(unsigned short val_nDim,
                                      unsigned short val_nVar,
                                      CConfig        *config);

  /*!
   * \brief Source term integration for boundary heat sink.
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig *config) override;

};


/*!
 * \class CSourceRadiation
 * \brief Class for a source term due to radiation.
 * \ingroup SourceDiscr
 * \author Ruben Sanchez
 */
class CSourceRadiation : public CSourceBase_Flow {
private:
  bool implicit;

public:

  CSourceRadiation(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Source term integration for a radiation source.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CSourceBAYModel
 * \brief Source term for BAY model for Vortex Generators
 * \ingroup SourceDiscr
 * \author M. Lagorio
 */
class CSourceBAYModel : public CSourceBase_Flow {
 private:
  bool implicit;                 /*!< \brief Implicit Calculation */
  bool first_computation{true}; /*!< \brief Flag for MPI comunication*/
  su2double calibrationConstant; /*!< \bried Calibration constant for the BAY model force calculation*/
  unsigned short nVgs{0}; /*!< \brief Number of Vortex Generators */

  /*!
   *\class Edge_info_VGModel
   *\brief Data structure to store the edge information for the BAY model
   */
  struct Edge_info_VGModel {
    unsigned long iPoint, jPoint;
    su2double iDistance, jDistance, cells_volume;
};
  /*!
   *\class Vortex_Generator
   *\brief Data structure to store the Vortex Generator information for the BAY model
   */
  public:
  class Vortex_Generator {
  private:

    public:
    su2double** coords_vg = nullptr; /*!< \brief Data structure to store the points defining the VG */
    su2double* t=nullptr; /*!< \brief Vectors defining the VG directions */
    su2double* b=nullptr;
    su2double* n=nullptr;

    su2double Vtot{0.0};        /*!< \brief Total Volume of the cells defining the Vortex Generator Volume*/
    unsigned short nPoints;


    map<unsigned long, Edge_info_VGModel*>
        EdgesBay; /*!< \brief Data structure to store edge informations defining the BAY model */

    map<unsigned long, unsigned long> pointsBAY; /*!< \brief Structure to map the mesh nodes to edges */

    su2double Svg; /*!< \brief Platform area of the VG */

    /*!
     * \brief VG constructor
     * \param[in] config - Definition of the particular problem
     * \param[in] iVG - ndex of the vortex generator 
     */
    Vortex_Generator(const CConfig* config, const unsigned short iVG);

    /*!
     * \brief Class deconstructor
     */
     ~Vortex_Generator();

    /*!
     * \brief Add cell volume to total volume of the VG
     * \param[in] Vcell - Cell volume
     */
    void addVGcellVolume(const su2double Vcell) { Vtot += Vcell; }

    /*!
     * \brief Checks if the point is alredy selected by anothe edge
     * \param[in] Point - Point to check
     * \param[in] jEdge - Variable to store the previous edge
     * \param[in] distanceOld - Variable to store the distance between the Vg and the previous point
     */
    bool Check_edge_map(const unsigned long Point, unsigned long &jEdge,su2double &distanceOld);

    /*!
     * \brief Checks if the edge intersects the vortex generator
     * \param[in] Coord_i - Coordinates of the first point
     * \param[in] Coord_j - Coodinates of the second point
     * \param[in] Normal - Normal vector
     * \param[out] distanceToVg - Distance between the edge and the VG
     */
    bool EdgeIntersectsVG(const su2double *Coord_i,const su2double *Coord_j, const su2double *Normal,su2double &distanceToVg);

    /*!
     * \brief Method returns if a point is inside the VG
     * \param[in] vgModel - Type of VG model
     * \param[in] Point_i - Index of the point
     * \param[in] Edge - Index of the edge
     * \return True if the point is inside the VG
     */
    const bool PointInVG(const ENUM_VG_MODEL vgModel,const unsigned long& Point_i,const unsigned long& Edge);

    /*!
     * \brief Method returns interpolation coefficient for obtaining velocity at VG plane
     * \param[in] vgModel - Type of VG model
     * \param[in] Edge - Index of the edge
     * \param[in] Point_i - Index of the point
     * \return Interpolation coefficient
     */
    const su2double Get_interpolationCoefficient(ENUM_VG_MODEL vgModel,const unsigned long Edge, const unsigned long Point_i);

    /*!
     * \brief Method returns jBAY redistribution constant
     * \param[in] vgModel - Type of VG model
     * \param[in] Edge - Index of the edge
     * \param[in] Point_i - Index of the point
     * \param[in] Volume - Volume of the cells
     * \return Interpolation coefficient
     */
    const su2double Get_redistributionConstant(const ENUM_VG_MODEL vgModel,const unsigned long Edge, const unsigned long Point_i, const su2double Volume);
  };

  vector<Vortex_Generator*> VGs{nullptr}; /*!< \brief Vector storing all VGs data structures */


public:
  CSourceBAYModel(unsigned short val_ndim, unsigned short val_nVar, const CConfig* config);

  ~CSourceBAYModel() {
    for (auto* iVg : VGs) {
      delete iVg;
    }
    VGs.clear();
  }

  /*!
   * \brief Source term integration for a VG Bay model source.
   * \param[in] config - Definition of the particular problem.
   * \return Lightweight const-view of residual and Jacobian.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

  void ReduceVariables(void);

  /*!
   * \brief Source term initialization
   */
  void UpdateSource(const CConfig* config) override;

};
