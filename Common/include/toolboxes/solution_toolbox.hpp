/*!
 * \file solution_toolbox.hpp
 * \brief Header file for the verification solution classes.
 *        The implementations are in the <i>solution_toolbox.cpp</i> file.
 * \author T. Economon, E. van der Weide
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

#include <cmath>
#include "../config_structure.hpp"

/*!
 * \class CVerificationSolution
 * \brief Class for holding verification PDE solutions, e.g., phi = phi(x,y,z,t),
 *        used for initial conditions, analytic solutions, manufactured solutions.
 * \author T. Economon, E. van der Weide
 */
class CVerificationSolution {
  
protected:
  
  int rank;  /*!< \brief MPI Rank. */
  int size;  /*!< \brief MPI Size. */
  
  unsigned short nDim;  /*!< \brief Number of dimension of the problem. */
  unsigned short nVar;  /*!< \brief Number of variables of the problem  */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CVerificationSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config   - Definition of the particular problem.
   */
  CVerificationSolution(unsigned short val_nDim,
                        unsigned short val_nvar,
                        CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CVerificationSolution(void);

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  virtual void GetSolution(const unsigned short val_nParams,
                           const su2double      *val_params,
                           const su2double      *val_coords,
                           const su2double      val_t,
                           su2double            *val_solution);
  
  /*!
   * \brief Get the exact solution at the current position and t = 0.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetInitialCondition(const unsigned short val_nParams,
                           const su2double      *val_params,
                           const su2double      *val_coords,
                           su2double            *val_solution);
  
  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  virtual void GetBCState(const unsigned short val_nParams,
                          const su2double      *val_params,
                          const su2double      *val_coords,
                          const su2double      val_t,
                          su2double            *val_solution);
  
  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  virtual void GetMMSSourceTerm(const unsigned short val_nParams,
                                const su2double      *val_params,
                                const su2double      *val_coords,
                                const su2double      val_t,
                                su2double            *val_source);

  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - False as default value. Overwrite this function for a
                manufactured solution.
   */
  virtual bool IsManufacturedSolution(void);

  /*!
   * \brief Whether or not the exact solution is known for this verification solution.
   * \return  - True as default value. Overwrite this function if the exacti
                solution is not known.
   */
  virtual bool ExactSolutionKnown(void);
  
  /*!
   * \brief Get the local error defined as the local solution minus the verification solution.
   * \param[in]  val_nParams  - Number of additional input parameters.
   * \param[in]  val_params   - Array of additional input parameters.
   * \param[in]  val_coords   - Cartesian coordinates of the current position.
   * \param[in]  val_solution - Array where the exact solution is stored.
   * \param[out] val_error    - Array where the local error is stored.
   */
  void GetLocalError(const unsigned short val_nParams,
                     const su2double      *val_params,
                     const su2double      *val_coords,
                     const su2double      val_t,
                     const su2double      *GetLocalErrorval_solution,
                     su2double            *val_error);
};

/*!
 * \class CInviscidVortexSolution
 * \brief Class to define the required data for the Inviscid Vortex.
 * \author E. van der Weide, T. Economon
 */
class CInviscidVortexSolution: public CVerificationSolution {

protected:

  /*--- Specific conditions for the inviscid vortex. ---*/
  su2double MachVortex;     /*!< \brief Mach number of the undisturbed flow. */
  su2double x0Vortex;       /*!< \brief Initial x-coordinate of the vortex center. */
  su2double y0Vortex;       /*!< \brief Initial y-coordinate of the vortex center. */
  su2double RVortex;        /*!< \brief Radius of the vortex. */
  su2double epsVortex;      /*!< \brief Strength of the vortex. */
  su2double thetaVortex;    /*!< \brief Advection angle (in degrees) of the vortex. */

  /*--- Variables involving gamma. */
  su2double Gamma;        /*!< \brief Gamma */
  su2double Gm1;          /*!< \brief Gamma minus 1 */
  su2double ovGm1;        /*!< \brief 1 over Gamma minus 1 */
  su2double gamOvGm1;     /*!< \brief Gamma over Gamma minus 1 */

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CInviscidVortexSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config   - Configuration of the particular problem.
   */
  CInviscidVortexSolution(unsigned short val_nDim,
                          unsigned short val_nvar,
                          CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CInviscidVortexSolution(void);

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const unsigned short val_nParams,
                   const su2double      *val_params,
                   const su2double      *val_coords,
                   const su2double      val_t,
                   su2double            *val_solution);

  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const unsigned short val_nParams,
                  const su2double      *val_params,
                  const su2double      *val_coords,
                  const su2double      val_t,
                  su2double            *val_solution);
};

/*!
 * \class CRinglebSolution
 * \brief Class to define the required data for the Ringleb flow.
 * \author E. van der Weide, T. Economon
 */
class CRinglebSolution: public CVerificationSolution {

protected:

  /*--- Variables involving gamma. ---*/
  su2double Gamma;        /*!< \brief Gamma */
  su2double Gm1;          /*!< \brief Gamma minus 1 */
  su2double tovGm1;       /*!< \brief 2 over Gamma minus 1 */
  su2double tGamOvGm1;    /*!< \brief 2 Gamma over Gamma minus 1 */

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CRinglebSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config   - Configuration of the particular problem.
   */
  CRinglebSolution(unsigned short val_nDim,
                   unsigned short val_nvar,
                   CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CRinglebSolution(void);

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const unsigned short val_nParams,
                   const su2double      *val_params,
                   const su2double      *val_coords,
                   const su2double      val_t,
                   su2double            *val_solution);

  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const unsigned short val_nParams,
                  const su2double      *val_params,
                  const su2double      *val_coords,
                  const su2double      val_t,
                  su2double            *val_solution);
};

/*!
 * \class CNSUnitQuadSolution
 * \brief Class to define the required data for the Navier-Stokes solution
          on a unit quad, heat conduction is neglected.
 * \author E. van der Weide, T. Economon
 */
class CNSUnitQuadSolution: public CVerificationSolution {

protected:

  /*--- Variables that define the soltion. ---*/
  su2double Gm1;          /*!< \brief Gamma minus 1 */
  su2double flowAngle;    /*!< \brief Angle of the velocity vector in radians. */
  su2double Viscosity;    /*!< \brief Viscosity, must be constant. */

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CNSUnitQuadSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config   - Configuration of the particular problem.
   */
  CNSUnitQuadSolution(unsigned short val_nDim,
                      unsigned short val_nvar,
                      CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CNSUnitQuadSolution(void);

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const unsigned short val_nParams,
                   const su2double      *val_params,
                   const su2double      *val_coords,
                   const su2double      val_t,
                   su2double            *val_solution);
  
  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const unsigned short val_nParams,
                  const su2double      *val_params,
                  const su2double      *val_coords,
                  const su2double      val_t,
                  su2double            *val_solution);
};

/*!
 * \class CTGVSolution
 * \brief Class to define the required data for the Taylor Green Vortex.
 * \author E. van der Weide, T. Economon
 */
class CTGVSolution: public CVerificationSolution {
  
protected:
  
  /*--- TGV specific conditions. ---*/
  
  su2double tgvLength;    /*!< \brief Taylor-Green length scale. */
  su2double tgvVelocity;  /*!< \brief Taylor-Green velocity. */
  su2double tgvDensity;   /*!< \brief Taylor-Green density. */
  su2double tgvPressure;  /*!< \brief Taylor-Green pressure. */
  su2double ovGm1;        /*!< \brief 1 over Gamma minus 1 */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CTGVSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config   - Configuration of the particular problem.
   */
  CTGVSolution(unsigned short val_nDim,
               unsigned short val_nvar,
               CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CTGVSolution(void);
  
  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const unsigned short val_nParams,
                   const su2double      *val_params,
                   const su2double      *val_coords,
                   const su2double      val_t,
                   su2double            *val_solution);
  
  /*!
   * \brief Whether or not the exact solution is known for this verification solution.
   * \return  - False, because the exact solution is not known for the TGV case.
   */
  bool ExactSolutionKnown(void);
};

/*!
 * \class CIncTGVSolution
 * \brief Class to define the required data for the incompressible Taylor Green Vortex.
 * \author T. Economon, E. van der Weide
 */
class CIncTGVSolution: public CVerificationSolution {
  
protected:
  
  /*--- TGV specific conditions. ---*/
  
  su2double tgvLength;    /*!< \brief Taylor-Green length scale. */
  su2double tgvVelocity;  /*!< \brief Taylor-Green velocity. */
  su2double tgvDensity;   /*!< \brief Taylor-Green density. */
  su2double tgvViscosity; /*!< \brief Taylor-Green viscosity. */
  
  su2double Temperature;  /*!< \brief Temperature, just to be safe. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CIncTGVSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config   - Configuration of the particular problem.
   */
  CIncTGVSolution(unsigned short val_nDim,
                  unsigned short val_nvar,
                  CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CIncTGVSolution(void);
  
  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const unsigned short val_nParams,
                   const su2double      *val_params,
                   const su2double      *val_coords,
                   const su2double      val_t,
                   su2double            *val_solution);
  
  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const unsigned short val_nParams,
                  const su2double      *val_params,
                  const su2double      *val_coords,
                  const su2double      val_t,
                  su2double            *val_solution);
};

/*!
 * \class CMMSNSUnitQuadSolution
 * \brief Class to define the required data for the manufactured solution of the
          laminar Navier-Stokes equations on a unit quad.
 * \author E. van der Weide, T. Economon
 */
class CMMSNSUnitQuadSolution: public CVerificationSolution {

protected:

  /*--- Variables that define the solution and MMS source term. ---*/
  su2double Gamma;        /*!< \brief Specific heat ratio. */
  su2double RGas;         /*!< \brief Gas constant. */
  su2double Viscosity;    /*!< \brief Viscosity, must be constant. */
  su2double Conductivity; /*!< \brief Thermal conductivity, must be constant. */

  /*--- Constants, which describe this manufactured solution. This is a viscous
        solution on the unit quad, where the primitive variables vary as a
        combination of sine and cosine functions. The unit quad is probably not
        necessary, and an arbitrary domain should work as well. ---*/

  su2double L;        /*!< \brief Length scale. */
  su2double a_Px;     /*!< \brief Parameter for the pressure solution. */
  su2double a_Pxy;    /*!< \brief Parameter for the pressure solution. */
  su2double a_Py;     /*!< \brief Parameter for the pressure solution. */
  su2double a_rhox;   /*!< \brief Parameter for the density solution. */
  su2double a_rhoxy;  /*!< \brief Parameter for the density solution. */
  su2double a_rhoy;   /*!< \brief Parameter for the density solution. */
  su2double a_ux;     /*!< \brief Parameter for the x-velocity solution. */
  su2double a_uxy;    /*!< \brief Parameter for the x-velocity solution. */
  su2double a_uy;     /*!< \brief Parameter for the x-velocity solution. */
  su2double a_vx;     /*!< \brief Parameter for the y-velocity solution. */
  su2double a_vxy;    /*!< \brief Parameter for the y-velocity solution. */
  su2double a_vy;     /*!< \brief Parameter for the y-velocity solution. */
  su2double P_0;      /*!< \brief Parameter for the pressure solution. */
  su2double P_x;      /*!< \brief Parameter for the pressure solution. */
  su2double P_xy;     /*!< \brief Parameter for the pressure solution. */
  su2double P_y;      /*!< \brief Parameter for the pressure solution. */ 
  su2double rho_0;    /*!< \brief Parameter for the density solution. */
  su2double rho_x;    /*!< \brief Parameter for the density solution. */
  su2double rho_xy;   /*!< \brief Parameter for the density solution. */
  su2double rho_y;    /*!< \brief Parameter for the density solution. */
  su2double u_0;      /*!< \brief Parameter for the x-velocity solution. */
  su2double u_x;      /*!< \brief Parameter for the x-velocity solution. */
  su2double u_xy;     /*!< \brief Parameter for the x-velocity solution. */
  su2double u_y;      /*!< \brief Parameter for the x-velocity solution. */
  su2double v_0;      /*!< \brief Parameter for the y-velocity solution. */
  su2double v_x;      /*!< \brief Parameter for the y-velocity solution. */
  su2double v_xy;     /*!< \brief Parameter for the y-velocity solution. */
  su2double v_y;      /*!< \brief Parameter for the y-velocity solution. */

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CMMSNSUnitQuadSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config   - Configuration of the particular problem.
   */
  CMMSNSUnitQuadSolution(unsigned short val_nDim,
                         unsigned short val_nvar,
                         CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CMMSNSUnitQuadSolution(void);

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const unsigned short val_nParams,
                   const su2double      *val_params,
                   const su2double      *val_coords,
                   const su2double      val_t,
                   su2double            *val_solution);

  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const unsigned short val_nParams,
                  const su2double      *val_params,
                  const su2double      *val_coords,
                  const su2double      val_t,
                  su2double            *val_solution);

  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetMMSSourceTerm(const unsigned short val_nParams,
                        const su2double      *val_params,
                        const su2double      *val_coords,
                        const su2double      val_t,
                        su2double            *val_source);

  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - True, because this is a manufactured solution.
   */
  bool IsManufacturedSolution(void);
};

/*!
 * \class CMMSIncEulerSolution
 * \brief Class to define the required data for the manufactured solution of the
 *        incompressible Euler equations.
 * \author T. Economon, E. van der Weide
 */
class CMMSIncEulerSolution: public CVerificationSolution {
  
protected:
  
  /*--- Variables that define the solution and MMS source term. ---*/
  su2double Density;      /*!< \brief Density, must be constant. */
  su2double Temperature;  /*!< \brief Temperature, just to be safe. */
  
  /*--- Constants, which describe this manufactured solution. This is a
   solution where the primitive variables vary as a combination
   of sine and cosine functions. The solution is from Salari K, and
   Knupp P, "Code verification by the method of manufactured solutions,"
   SAND 2000-1444, Sandia National Laboratories, Albuquerque, NM, 2000. ---*/
  
  su2double P_0;      /*!< \brief Parameter for the pressure solution. */
  su2double u_0;      /*!< \brief Parameter for the x-velocity solution. */
  su2double v_0;      /*!< \brief Parameter for the y-velocity solution. */
  su2double epsilon;  /*!< \brief Parameter for the velocity solutions. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CMMSIncEulerSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config   - Configuration of the particular problem.
   */
  CMMSIncEulerSolution(unsigned short val_nDim,
                       unsigned short val_nvar,
                       CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CMMSIncEulerSolution(void);
  
  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const unsigned short val_nParams,
                   const su2double      *val_params,
                   const su2double      *val_coords,
                   const su2double      val_t,
                   su2double            *val_solution);
  
  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const unsigned short val_nParams,
                  const su2double      *val_params,
                  const su2double      *val_coords,
                  const su2double      val_t,
                  su2double            *val_solution);
  
  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetMMSSourceTerm(const unsigned short val_nParams,
                        const su2double      *val_params,
                        const su2double      *val_coords,
                        const su2double      val_t,
                        su2double            *val_source);
  
  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - True, because this is a manufactured solution.
   */
  bool IsManufacturedSolution(void);
};

/*!
 * \class CMMSIncNSSolution
 * \brief Class to define the required data for the manufactured solution of the
 *        laminar incompressible Navier-Stokes equations.
 * \author T. Economon, E. van der Weide
 */
class CMMSIncNSSolution: public CVerificationSolution {
  
protected:
  
  /*--- Variables that define the solution and MMS source term. ---*/
  su2double Viscosity;    /*!< \brief Viscosity, must be constant. */
  su2double Density;      /*!< \brief Density, must be constant. */
  su2double Temperature;  /*!< \brief Temperature, just to be safe. */
  
  /*--- Constants, which describe this manufactured solution. This is a
   viscous solution where the primitive variables vary as a combination
   of sine and cosine functions. The solution is from Salari K, and
   Knupp P, "Code verification by the method of manufactured solutions,"
   SAND 2000-1444, Sandia National Laboratories, Albuquerque, NM, 2000. ---*/
  
  su2double P_0;      /*!< \brief Parameter for the pressure solution. */
  su2double u_0;      /*!< \brief Parameter for the x-velocity solution. */
  su2double v_0;      /*!< \brief Parameter for the y-velocity solution. */
  su2double epsilon;  /*!< \brief Parameter for the velocity solutions. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CMMSIncNSSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config   - Configuration of the particular problem.
   */
  CMMSIncNSSolution(unsigned short val_nDim,
                    unsigned short val_nvar,
                    CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CMMSIncNSSolution(void);
  
  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const unsigned short val_nParams,
                   const su2double      *val_params,
                   const su2double      *val_coords,
                   const su2double      val_t,
                   su2double            *val_solution);
  
  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const unsigned short val_nParams,
                  const su2double      *val_params,
                  const su2double      *val_coords,
                  const su2double      val_t,
                  su2double            *val_solution);
  
  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetMMSSourceTerm(const unsigned short val_nParams,
                        const su2double      *val_params,
                        const su2double      *val_coords,
                        const su2double      val_t,
                        su2double            *val_source);
  
  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - True, because this is a manufactured solution.
   */
  bool IsManufacturedSolution(void);
};

/*!
 * \class CUserDefinedSolution
 * \brief Class to define the required data for a user defined solution.
 * \author E. van der Weide, T. Economon
 */
class CUserDefinedSolution: public CVerificationSolution {

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CUserDefinedSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config   - Configuration of the particular problem.
   */
  CUserDefinedSolution(unsigned short val_nDim,
                       unsigned short val_nvar,
                       CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUserDefinedSolution(void);

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const unsigned short val_nParams,
                   const su2double      *val_params,
                   const su2double      *val_coords,
                   const su2double      val_t,
                   su2double            *val_solution);

  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const unsigned short val_nParams,
                  const su2double      *val_params,
                  const su2double      *val_coords,
                  const su2double      val_t,
                  su2double            *val_solution);

  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_nParams  - Number of additional input parameters.
   * \param[in] val_params   - Array of additional input parameters.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetMMSSourceTerm(const unsigned short val_nParams,
                        const su2double      *val_params,
                        const su2double      *val_coords,
                        const su2double      val_t,
                        su2double            *val_source);

  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - True if this is a manufactured solution and false otherwise.
   */
  bool IsManufacturedSolution(void);
};
