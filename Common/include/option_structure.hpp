/*!
 * \file option_structure.hpp
 * \brief Defines classes for referencing options for easy input in CConfig
 * \author J. Hicken, B. Tracey
 * \version 7.0.4 "Blackbird"
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

#include "./mpi_structure.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>

using namespace std;

/*!
 * \class CEmptyMap
 * \brief We use this dummy class instead of std::map when
 * we only need the enum definition and not the string to
 * enum maps, this makes compilation much faster.
 */
template <typename T, typename U>
struct CEmptyMap {
  CEmptyMap(initializer_list<pair<const T, U> >) {}
};

#ifdef ENABLE_MAPS
template<class T, class U>
using MapType = map<T,U>;
#define MakePair(a,b) {a,b},
#else
template<class T, class U>
using MapType = CEmptyMap<T,U>;
#define MakePair(a,b)
#endif

/*!
 * \brief Different software components of SU2
 */
enum SU2_COMPONENT {
  SU2_CFD = 1,	/*!< \brief Running the SU2_CFD software. */
  SU2_DEF = 2,	/*!< \brief Running the SU2_DEF software. */
  SU2_DOT = 3,	/*!< \brief Running the SU2_DOT software. */
  SU2_MSH = 4,	/*!< \brief Running the SU2_MSH software. */
  SU2_GEO = 5,	/*!< \brief Running the SU2_GEO software. */
  SU2_SOL = 6 	/*!< \brief Running the SU2_SOL software. */
};

const unsigned int EXIT_DIVERGENCE = 2;   /*!< \brief Exit code (divergence). */

const unsigned int BUFSIZE = 3000000;         /*!< \brief MPI buffer. */
const unsigned int MAX_PARAMETERS = 10;       /*!< \brief Maximum number of parameters for a design variable definition. */
const unsigned int MAX_NUMBER_PERIODIC = 10;  /*!< \brief Maximum number of periodic boundary conditions. */
const unsigned int MAX_STRING_SIZE = 200;     /*!< \brief Maximum number of domains. */
const unsigned int MAX_NUMBER_FFD = 15;       /*!< \brief Maximum number of FFDBoxes for the FFD. */
const unsigned int MAX_SOLS = 12;             /*!< \brief Maximum number of solutions at the same time (dimension of solution container array). */
const unsigned int MAX_TERMS = 6;             /*!< \brief Maximum number of terms in the numerical equations (dimension of solver container array). */
const unsigned int MAX_ZONES = 3;             /*!< \brief Maximum number of zones. */
const unsigned int MAX_FE_KINDS = 4;          /*!< \brief Maximum number of Finite Elements. */
const unsigned int NO_RK_ITER = 0;            /*!< \brief No Runge-Kutta iteration. */

const unsigned int OVERHEAD = 4;    /*!< \brief Overhead space above nMarker when allocating space for boundary elems (MPI + periodic). */

const unsigned int MESH_0 = 0;  /*!< \brief Definition of the finest grid level. */
const unsigned int MESH_1 = 1;  /*!< \brief Definition of the finest grid level. */
const unsigned int ZONE_0 = 0;  /*!< \brief Definition of the first grid domain. */
const unsigned int ZONE_1 = 1;  /*!< \brief Definition of the second grid domain. */
const unsigned int INST_0 = 0;  /*!< \brief Definition of the first instance per grid level. */

const su2double STANDARD_GRAVITY = 9.80665;           /*!< \brief Acceleration due to gravity at surface of earth. */
const su2double UNIVERSAL_GAS_CONSTANT = 8.3144598;   /*!< \brief Universal gas constant in J/(mol*K) */

const su2double EPS = 1.0E-16;        /*!< \brief Error scale. */
const su2double TURB_EPS = 1.0E-16;   /*!< \brief Turbulent Error scale. */

const su2double ONE2 = 0.5;         /*!< \brief One divided by two. */
const su2double TWO3 = 2.0 / 3.0;   /*!< \brief Two divided by three. */
const su2double FOUR3 = 4.0 / 3.0;  /*!< \brief Four divided by three. */

const su2double PI_NUMBER = 4.0 * atan(1.0);  /*!< \brief Pi number. */

const su2double STEFAN_BOLTZMANN = 5.670367E-08;  /*!< \brief Stefan-Boltzmann constant in W/(m^2*K^4). */

const int MASTER_NODE = 0;			/*!< \brief Master node for MPI parallelization. */
const int SINGLE_NODE = 1;			/*!< \brief There is only a node in the MPI parallelization. */
const int SINGLE_ZONE = 1;			/*!< \brief There is only a zone. */

const unsigned short COMM_TYPE_UNSIGNED_LONG  = 1;  /*!< \brief Communication type for unsigned long. */
const unsigned short COMM_TYPE_LONG           = 2;  /*!< \brief Communication type for long. */
const unsigned short COMM_TYPE_UNSIGNED_SHORT = 3;  /*!< \brief Communication type for unsigned short. */
const unsigned short COMM_TYPE_DOUBLE         = 4;  /*!< \brief Communication type for double. */
const unsigned short COMM_TYPE_CHAR           = 5;  /*!< \brief Communication type for char. */
const unsigned short COMM_TYPE_SHORT          = 6;  /*!< \brief Communication type for short. */
const unsigned short COMM_TYPE_INT            = 7;  /*!< \brief Communication type for int. */

const unsigned short N_ELEM_TYPES = 7;           /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_LINE = 2;          /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_TRIANGLE = 3;      /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_QUADRILATERAL = 4; /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_TETRAHEDRON = 4;   /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_HEXAHEDRON = 8;    /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_PYRAMID = 5;       /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_PRISM = 6;         /*!< \brief General output & CGNS defines. */
enum: unsigned short{N_POINTS_MAXIMUM = 8};      /*!< \brief Max. out of the above, used for static arrays, keep it up to date. */

const int CGNS_STRING_SIZE = 33; /*!< \brief Length of strings used in the CGNS format. */
const int SU2_CONN_SIZE   = 10;  /*!< \brief Size of the connectivity array that is allocated for each element
                                             that we read from a mesh file in the format [[globalID vtkType n0 n1 n2 n3 n4 n5 n6 n7 n8]. */
const int SU2_CONN_SKIP   = 2;   /*!< \brief Offset to skip the globalID and VTK type at the start of the element connectivity list for each CGNS element. */

const su2double COLORING_EFF_THRESH = 0.875;  /*!< \brief Below this value fallback strategies are used instead. */

/*!
 * \brief Boolean answers
 */
enum ANSWER {
  NONE = 0,
  NO = 0,   /*!< \brief Boolean definition of no. */
  YES = 1   /*!< \brief Boolean definition of yes. */
};

/*!
 * \brief Average method for marker analyze
 */
enum AVERAGE_TYPE {
  AVERAGE_AREA = 1,     /*!< \brief Area-weighted average. */
  AVERAGE_MASSFLUX = 2  /*!< \brief Mass-flux weighted average. */
};
static const MapType<string, AVERAGE_TYPE> Average_Map = {
  MakePair("AREA", AVERAGE_AREA)
  MakePair("MASSFLUX", AVERAGE_MASSFLUX)
};

/*!
 * \brief different solver types for the CFD component
 */
enum ENUM_MAIN_SOLVER {
  NO_SOLVER = 0,                    /*!< \brief Definition of no solver. */
  EULER = 1,                        /*!< \brief Definition of the Euler's solver. */
  NAVIER_STOKES = 2,                /*!< \brief Definition of the Navier-Stokes' solver. */
  RANS = 3,                         /*!< \brief Definition of the Reynolds-averaged Navier-Stokes' (RANS) solver. */
  INC_EULER = 4,                    /*!< \brief Definition of the incompressible Euler's solver. */
  INC_NAVIER_STOKES =5,             /*!< \brief Definition of the incompressible Navier-Stokes' solver. */
  INC_RANS = 6,                     /*!< \brief Definition of the incompressible Reynolds-averaged Navier-Stokes' (RANS) solver. */
  HEAT_EQUATION = 7,                /*!< \brief Definition of the finite volume heat solver. */
  FLUID_STRUCTURE_INTERACTION = 8,  /*!< \brief Definition of a FSI solver. */
  FEM_ELASTICITY = 9,               /*!< \brief Definition of a FEM solver. */
  ADJ_EULER = 10,                   /*!< \brief Definition of the continuous adjoint Euler's solver. */
  ADJ_NAVIER_STOKES = 11,           /*!< \brief Definition of the continuous adjoint Navier-Stokes' solver. */
  ADJ_RANS = 12,                    /*!< \brief Definition of the continuous adjoint Reynolds-averaged Navier-Stokes' (RANS) solver. */
  TEMPLATE_SOLVER = 13,             /*!< \brief Definition of template solver. */
  DISC_ADJ_EULER = 15,              /*!< \brief Definition of the discrete adjoint Euler solver. */
  DISC_ADJ_RANS = 16,               /*!< \brief Definition of the discrete adjoint Reynolds-averaged Navier-Stokes' (RANS) solver. */
  DISC_ADJ_NAVIER_STOKES = 17,      /*!< \brief Definition of the discrete adjoint Navier-Stokes' solver. */
  DISC_ADJ_INC_EULER = 18,          /*!< \brief Definition of the discrete adjoint incompressible Euler solver. */
  DISC_ADJ_INC_RANS = 19,           /*!< \brief Definition of the discrete adjoint imcompressible Reynolds-averaged Navier-Stokes' (RANS) solver. */
  DISC_ADJ_INC_NAVIER_STOKES = 20,  /*!< \brief Definition of the discrete adjoint imcompressible Navier-Stokes'. */
  DISC_ADJ_HEAT = 21,               /*!< \brief Definition of the discrete adjoint heat solver. */
  DISC_ADJ_FEM_EULER = 22,          /*!< \brief Definition of the discrete adjoint FEM Euler solver. */
  DISC_ADJ_FEM_RANS = 23,           /*!< \brief Definition of the discrete adjoint FEM Reynolds-averaged Navier-Stokes' (RANS) solver. */
  DISC_ADJ_FEM_NS = 24,             /*!< \brief Definition of the discrete adjoint FEM Navier-Stokes' solver. */
  DISC_ADJ_FEM = 25,                /*!< \brief Definition of the discrete adjoint FEM solver. */
  FEM_EULER = 26,                   /*!< \brief Definition of the finite element Euler's solver. */
  FEM_NAVIER_STOKES = 27,           /*!< \brief Definition of the finite element Navier-Stokes' solver. */
  FEM_RANS = 28,                    /*!< \brief Definition of the finite element Reynolds-averaged Navier-Stokes' (RANS) solver. */
  FEM_LES = 29,                     /*!< \brief Definition of the finite element Large Eddy Simulation Navier-Stokes' (LES) solver. */
  MULTIPHYSICS = 30
};
static const MapType<string, ENUM_MAIN_SOLVER> Solver_Map = {
  MakePair("NONE", NO_SOLVER)
  MakePair("EULER", EULER)
  MakePair("NAVIER_STOKES", NAVIER_STOKES)
  MakePair("RANS", RANS)
  MakePair("INC_EULER", INC_EULER)
  MakePair("INC_NAVIER_STOKES", INC_NAVIER_STOKES)
  MakePair("INC_RANS", INC_RANS)
  MakePair("FEM_EULER", FEM_EULER)
  MakePair("FEM_NAVIER_STOKES", FEM_NAVIER_STOKES)
  MakePair("FEM_RANS", FEM_RANS)
  MakePair("FEM_LES", FEM_LES)
  MakePair("ADJ_EULER", ADJ_EULER)
  MakePair("ADJ_NAVIER_STOKES", ADJ_NAVIER_STOKES)
  MakePair("ADJ_RANS", ADJ_RANS )
  MakePair("HEAT_EQUATION", HEAT_EQUATION)
  MakePair("ELASTICITY", FEM_ELASTICITY)
  MakePair("DISC_ADJ_EULER", DISC_ADJ_EULER)
  MakePair("DISC_ADJ_RANS", DISC_ADJ_RANS)
  MakePair("DISC_ADJ_NAVIERSTOKES", DISC_ADJ_NAVIER_STOKES)
  MakePair("DISC_ADJ_INC_EULER", DISC_ADJ_INC_EULER)
  MakePair("DISC_ADJ_INC_RANS", DISC_ADJ_INC_RANS)
  MakePair("DISC_ADJ_INC_NAVIERSTOKES", DISC_ADJ_INC_NAVIER_STOKES)
  MakePair("DISC_ADJ_HEAT_EQUATION", DISC_ADJ_HEAT)
  MakePair("DISC_ADJ_FEM_EULER", DISC_ADJ_FEM_EULER)
  MakePair("DISC_ADJ_FEM_RANS", DISC_ADJ_FEM_RANS)
  MakePair("DISC_ADJ_FEM_NS", DISC_ADJ_FEM_NS)
  MakePair("DISC_ADJ_FEM", DISC_ADJ_FEM)
  MakePair("FLUID_STRUCTURE_INTERACTION", FLUID_STRUCTURE_INTERACTION)
  MakePair("TEMPLATE_SOLVER", TEMPLATE_SOLVER)
  MakePair("MULTIPHYSICS", MULTIPHYSICS)
};

/*!
 * \brief different solver types for the multizone environment component
 */
enum ENUM_MULTIZONE {
  MZ_BLOCK_GAUSS_SEIDEL = 0,   /*!< \brief Definition of a Block-Gauss-Seidel multizone solver. */
  MZ_BLOCK_JACOBI = 1          /*!< \brief Definition of a Block-Jacobi solver. */
};
static const MapType<string, ENUM_MULTIZONE> Multizone_Map = {
  MakePair("BLOCK_GAUSS_SEIDEL", MZ_BLOCK_GAUSS_SEIDEL)
  MakePair("BLOCK_JACOBI", MZ_BLOCK_JACOBI)
};

/*!
 * \brief Types of fluid solvers
 */
enum ENUM_FSI_FLUID_PROBLEM {
  NO_SOLVER_FFSI = 0,      /*!< \brief Definition of no solver. */
  EULER_FFSI = 1,          /*!< \brief Euler equations for the FSI problem */
  NAVIER_STOKES_FFSI = 2,  /*!< \brief NS equations for the FSI problem */
  RANS_FFSI = 3            /*!< \brief RANS equations for the FSI problem */
};
static const MapType<string, ENUM_FSI_FLUID_PROBLEM> FSI_Fluid_Solver_Map = {
  MakePair("NONE", NO_SOLVER_FFSI)
  MakePair("EULER", EULER_FFSI)
  MakePair("NAVIER_STOKES", NAVIER_STOKES_FFSI)
  MakePair("RANS", RANS_FFSI)
};

/*!
 * \brief Types of structural solvers
 */
enum ENUM_FSI_STRUC_PROBLEM {
  NO_SOLVER_SFSI = 0,           /*!< \brief Definition of no solver. */
  FEM_ELASTICITY_SFSI = 9,      /*!< \brief Nonlinear elasticity equations for the FSI problem */
};
static const MapType<string, ENUM_FSI_STRUC_PROBLEM> FSI_Struc_Solver_Map = {
  MakePair("NONE", NO_SOLVER_SFSI)
  MakePair("ELASTICITY", FEM_ELASTICITY_SFSI)
};

/*!
 * \brief Material geometric conditions
 */
enum ENUM_STRUCT_SOLVER {
  SMALL_DEFORMATIONS = 0,       /*!< \brief Definition of linear elastic material. */
  LARGE_DEFORMATIONS = 1,       /*!< \brief Definition of Neo-Hookean material. */
};
static const MapType<string, ENUM_STRUCT_SOLVER> Struct_Map = {
  MakePair("SMALL_DEFORMATIONS", SMALL_DEFORMATIONS)
  MakePair("LARGE_DEFORMATIONS", LARGE_DEFORMATIONS)
};

/*!
 * \brief Material model
 */
enum ENUM_MATERIAL_MODEL {
  LINEAR_ELASTIC = 0,   /*!< \brief Definition of linear elastic material. */
  NEO_HOOKEAN = 1,      /*!< \brief Definition of Neo-Hookean material. */
  KNOWLES = 2,          /*!< \brief Definition of Knowles stored-energy potential */
  IDEAL_DE = 3          /*!< \brief Definition of ideal Dielectric Elastomer */
};
static const MapType<string, ENUM_MATERIAL_MODEL> Material_Map = {
  MakePair("LINEAR_ELASTIC", LINEAR_ELASTIC)
  MakePair("NEO_HOOKEAN", NEO_HOOKEAN)
  MakePair("KNOWLES", KNOWLES)
  MakePair("IDEAL_DE", IDEAL_DE)
};

/*!
 * \brief Material compressibility
 */
enum ENUM_MAT_COMPRESS {
  COMPRESSIBLE_MAT = 0,           /*!< \brief Definition of compressible material. */
  NEARLY_INCOMPRESSIBLE_MAT = 1,  /*!< \brief Definition of nearly incompressible material. */
};
static const MapType<string, ENUM_MAT_COMPRESS> MatComp_Map = {
  MakePair("COMPRESSIBLE", COMPRESSIBLE_MAT)
  MakePair("NEARLY_INCOMPRESSIBLE", NEARLY_INCOMPRESSIBLE_MAT)
};

/*!
 * \brief Types of interpolators
 */
enum ENUM_INTERPOLATOR {
  NEAREST_NEIGHBOR = 0,      /*!< \brief Nearest Neigbhor interpolation */
  ISOPARAMETRIC = 1,         /*!< \brief Isoparametric interpolation, use CONSERVATIVE_INTERPOLATION=YES for conservative interpolation (S.A. Brown 1997).*/
  WEIGHTED_AVERAGE = 3,      /*!< \brief Sliding Mesh Approach E. Rinaldi 2015 */
  RADIAL_BASIS_FUNCTION = 4, /*!< \brief Radial basis function interpolation. */
};
static const MapType<string, ENUM_INTERPOLATOR> Interpolator_Map = {
  MakePair("NEAREST_NEIGHBOR", NEAREST_NEIGHBOR)
  MakePair("ISOPARAMETRIC",    ISOPARAMETRIC)
  MakePair("WEIGHTED_AVERAGE", WEIGHTED_AVERAGE)
  MakePair("RADIAL_BASIS_FUNCTION", RADIAL_BASIS_FUNCTION)
};

/*!
 * \brief Types of radial basis functions
 */
enum ENUM_RADIALBASIS {
  WENDLAND_C2 = 0,        /*!< \brief Wendland C2 radial basis function. */
  INV_MULTI_QUADRIC = 1,  /*!< \brief Inversed multi quartic biharmonic spline. */
  GAUSSIAN = 2,           /*!< \brief Gaussian basis function. */
  THIN_PLATE_SPLINE = 3,  /*!< \brief Thin plate spline. */
  MULTI_QUADRIC = 4,      /*!< \brief Multi quartic biharmonic spline. */
};
static const MapType<string, ENUM_RADIALBASIS> RadialBasisFunction_Map = {
  MakePair("WENDLAND_C2", WENDLAND_C2)
  MakePair("INV_MULTI_QUADRIC", INV_MULTI_QUADRIC)
  MakePair("GAUSSIAN", GAUSSIAN)
  MakePair("THIN_PLATE_SPLINE", THIN_PLATE_SPLINE)
  MakePair("MULTI_QUADRIC", MULTI_QUADRIC)
};

/*!
 * \brief type of radial spanwise interpolation function for the inlet face
 */
enum ENUM_INLET_SPANWISEINTERPOLATION {
  NO_INTERPOLATION = 0,
  LINEAR_1D = 1,
  AKIMA_1D = 2,
};
static const map<string, ENUM_INLET_SPANWISEINTERPOLATION> Inlet_SpanwiseInterpolation_Map = {
  MakePair("NONE", NO_INTERPOLATION)
  MakePair("LINEAR_1D",LINEAR_1D)
  MakePair("AKIMA_1D",AKIMA_1D)
};

/*!
 * \brief type of radial spanwise interpolation data type for the inlet face
 */
enum ENUM_INLET_INTERPOLATIONTYPE {
  VR_VTHETA = 0,
  ALPHA_PHI = 1,
};
static const map<string, ENUM_INLET_INTERPOLATIONTYPE> Inlet_SpanwiseInterpolationType_Map = {
  MakePair("VR_VTHETA",VR_VTHETA)
  MakePair("ALPHA_PHI",ALPHA_PHI)
};

/*!
 * \brief types of (coupling) transfers between distinct physical zones
 */
enum ENUM_TRANSFER {
  ZONES_ARE_EQUAL                   = 0,    /*!< \brief Zones are equal - no transfer. */
  NO_COMMON_INTERFACE               = 1,    /*!< \brief No common interface between the zones (geometrical). */
  NO_TRANSFER                       = 2,    /*!< \brief Zones may share a boundary, but still no coupling desired. */
  FLOW_TRACTION                     = 10,   /*!< \brief Flow traction coupling (between fluids and solids). */
  BOUNDARY_DISPLACEMENTS            = 21,   /*!< \brief Boundary displacements (between fluids and solids) */
  SLIDING_INTERFACE                 = 13,   /*!< \brief Sliding interface (between fluids). */
  CONSERVATIVE_VARIABLES            = 14,   /*!< \brief General coupling that simply transfers the conservative variables (between same solvers). */
  MIXING_PLANE                      = 15,   /*!< \brief Mixing plane between fluids. */
  CONJUGATE_HEAT_FS                 = 16,   /*!< \brief Conjugate heat transfer (between compressible fluids and solids). */
  CONJUGATE_HEAT_WEAKLY_FS          = 17,   /*!< \brief Conjugate heat transfer (between incompressible fluids and solids). */
  CONJUGATE_HEAT_SF                 = 18,   /*!< \brief Conjugate heat transfer (between solids and compressible fluids). */
  CONJUGATE_HEAT_WEAKLY_SF          = 19,   /*!< \brief Conjugate heat transfer (between solids and incompressible fluids). */
};

/*!
 * \brief different regime modes
 */
enum ENUM_REGIME {
  COMPRESSIBLE = 0,		/*!< \brief Definition of compressible solver. */
  INCOMPRESSIBLE = 1,	/*!< \brief Definition of incompressible solver. */
  NO_FLOW = 2
};

/*!
 * \brief different non-dimensional modes
 */
enum ENUM_KIND_NONDIM {
  DIMENSIONAL = 0,              /*!< \brief Dimensional simulation (compressible or incompressible). */
  FREESTREAM_PRESS_EQ_ONE = 1,  /*!< \brief Non-dimensional compressible simulation with freestream pressure equal to 1.0. */
  FREESTREAM_VEL_EQ_MACH = 2,   /*!< \brief Non-dimensional compressible simulation with freestream velocity equal to Mach number. */
  FREESTREAM_VEL_EQ_ONE = 3,    /*!< \brief Non-dimensional compressible simulation with freestream pressure equal to 1.0. */
  INITIAL_VALUES   = 4,         /*!< \brief Non-dimensional incompressible simulation based on intial values for external flow. */
  REFERENCE_VALUES = 5          /*!< \brief Non-dimensional incompressible simulation based on custom reference values. */
};
static const MapType<string, ENUM_KIND_NONDIM> NonDim_Map = {
  MakePair("DIMENSIONAL", DIMENSIONAL)
  MakePair("FREESTREAM_PRESS_EQ_ONE", FREESTREAM_PRESS_EQ_ONE)
  MakePair("FREESTREAM_VEL_EQ_MACH",  FREESTREAM_VEL_EQ_MACH)
  MakePair("FREESTREAM_VEL_EQ_ONE",   FREESTREAM_VEL_EQ_ONE)
  MakePair("INITIAL_VALUES",   INITIAL_VALUES)
  MakePair("REFERENCE_VALUES", REFERENCE_VALUES)
};

/*!
 * \brief different system of measurements
 */
enum ENUM_MEASUREMENTS {
  SI = 0,			/*!< \brief Definition of compressible solver. */
  US = 1			/*!< \brief Definition of incompressible solver. */
};
static const MapType<string, ENUM_MEASUREMENTS> Measurements_Map = {
  MakePair("SI", SI)
  MakePair("US", US)
};

/*!
 * \brief different types of systems
 */
enum RUNTIME_TYPE {
  RUNTIME_FLOW_SYS = 2,       /*!< \brief One-physics case, the code is solving the flow equations(Euler and Navier-Stokes). */
  RUNTIME_TURB_SYS = 3,       /*!< \brief One-physics case, the code is solving the turbulence model. */
  RUNTIME_ADJPOT_SYS = 5,     /*!< \brief One-physics case, the code is solving the adjoint potential flow equation. */
  RUNTIME_ADJFLOW_SYS = 6,    /*!< \brief One-physics case, the code is solving the adjoint equations is being solved (Euler and Navier-Stokes). */
  RUNTIME_ADJTURB_SYS = 7,    /*!< \brief One-physics case, the code is solving the adjoint turbulence model. */
  RUNTIME_MULTIGRID_SYS = 14, /*!< \brief Full Approximation Storage Multigrid system of equations. */
  RUNTIME_FEA_SYS = 20,       /*!< \brief One-physics case, the code is solving the FEA equation. */
  RUNTIME_ADJFEA_SYS = 30,    /*!< \brief One-physics case, the code is solving the adjoint FEA equation. */
  RUNTIME_HEAT_SYS = 21,      /*!< \brief One-physics case, the code is solving the heat equation. */
  RUNTIME_ADJHEAT_SYS = 31,   /*!< \brief One-physics case, the code is solving the adjoint heat equation. */
  RUNTIME_TRANS_SYS = 22,     /*!< \brief One-physics case, the code is solving the turbulence model. */
  RUNTIME_RADIATION_SYS = 23, /*!< \brief One-physics case, the code is solving the radiation model. */
  RUNTIME_ADJRAD_SYS = 24,    /*!< \brief One-physics case, the code is solving the adjoint radiation model. */
};

const int FLOW_SOL = 0;     /*!< \brief Position of the mean flow solution in the solver container array. */
const int ADJFLOW_SOL = 1;  /*!< \brief Position of the continuous adjoint flow solution in the solver container array. */

const int TURB_SOL = 2;     /*!< \brief Position of the turbulence model solution in the solver container array. */
const int ADJTURB_SOL = 3;  /*!< \brief Position of the continuous adjoint turbulence solution in the solver container array. */

const int TRANS_SOL = 4;    /*!< \brief Position of the transition model solution in the solver container array. */
const int HEAT_SOL = 5;     /*!< \brief Position of the heat equation in the solution solver array. */
const int ADJHEAT_SOL = 6;  /*!< \brief Position of the adjoint heat equation in the solution solver array. */
const int RAD_SOL = 7;      /*!< \brief Position of the radiation equation in the solution solver array. */
const int ADJRAD_SOL = 8;   /*!< \brief Position of the continuous adjoint turbulence solution in the solver container array. */

const int MESH_SOL = 9;      /*!< \brief Position of the mesh solver. */
const int ADJMESH_SOL = 10;   /*!< \brief Position of the adjoint of the mesh solver. */

const int FEA_SOL = 0;      /*!< \brief Position of the FEA equation in the solution solver array. */
const int ADJFEA_SOL = 1;   /*!< \brief Position of the FEA adjoint equation in the solution solver array. */

const int TEMPLATE_SOL = 0; /*!< \brief Position of the template solution. */

const int CONV_TERM = 0;           /*!< \brief Position of the convective terms in the numerics container array. */
const int VISC_TERM = 1;           /*!< \brief Position of the viscous terms in the numerics container array. */
const int SOURCE_FIRST_TERM = 2;   /*!< \brief Position of the first source term in the numerics container array. */
const int SOURCE_SECOND_TERM = 3;  /*!< \brief Position of the second source term in the numerics container array. */
const int CONV_BOUND_TERM = 4;     /*!< \brief Position of the convective boundary terms in the numerics container array. */
const int VISC_BOUND_TERM = 5;     /*!< \brief Position of the viscous boundary terms in the numerics container array. */

const int FEA_TERM = 0;      /*!< \brief Position of the finite element analysis terms in the numerics container array. */
const int DE_TERM = 1;       /*!< \brief Position of the dielectric terms in the numerics container array. */

const int MAT_NHCOMP  = 2;   /*!< \brief Position of the Neo-Hookean compressible material model. */
const int MAT_IDEALDE = 3;   /*!< \brief Position of the Ideal-DE material model. */
const int MAT_KNOWLES = 4;   /*!< \brief Position of the Knowles material model. */

/*!
 * \brief Types of finite elements (in 2D or 3D)
 */
const int EL_TRIA = 0;    /*!< \brief Elements of three nodes (2D). */
const int EL_QUAD = 1;    /*!< \brief Elements of four nodes (2D). */

const int EL_TETRA = 0;   /*!< \brief Elements of four nodes (3D). */
const int EL_HEXA  = 1;   /*!< \brief Elements of eight nodes (3D). */
const int EL_PYRAM = 2;   /*!< \brief Elements of five nodes (3D). */
const int EL_PRISM = 3;   /*!< \brief Elements of six nodes (3D). */


/*!
 * \brief Types of mathematical problem to solve
 */
enum ENUM_MATH_PROBLEM {
  DIRECT = 0,               /*!< \brief Direct problem */
  CONTINUOUS_ADJOINT = 1,   /*!< \brief Continuous adjoint problem */
  DISCRETE_ADJOINT = 2      /*!< \brief AD-based discrete adjoint problem. */
};
static const MapType<string, ENUM_MATH_PROBLEM> Math_Problem_Map = {
  MakePair("DIRECT", DIRECT)
  MakePair("CONTINUOUS_ADJOINT", CONTINUOUS_ADJOINT)
  MakePair("DISCRETE_ADJOINT", DISCRETE_ADJOINT)
};

/*!
 * \brief Types of spatial discretizations
 */
enum ENUM_SPACE {
  NO_CONVECTIVE = 0,   /*!< \brief No convective scheme is used. */
  SPACE_CENTERED = 1,  /*!< \brief Space centered convective numerical method. */
  SPACE_UPWIND = 2,    /*!< \brief Upwind convective numerical method. */
  FINITE_ELEMENT = 3   /*!< \brief Finite element convective numerical method. */
};
static const MapType<string, ENUM_SPACE> Space_Map = {
  MakePair("NONE", NO_CONVECTIVE)
  MakePair("SPACE_CENTERED", SPACE_CENTERED)
  MakePair("SPACE_UPWIND", SPACE_UPWIND)
  MakePair("FINITE_ELEMENT", FINITE_ELEMENT)
};

/*!
 * \brief Types of fluid model
 */
enum ENUM_FLUIDMODEL {
  STANDARD_AIR = 0,       /*!< \brief Standard air gas model. */
  IDEAL_GAS = 1,          /*!< \brief Ideal gas model. */
  VW_GAS = 2,             /*!< \brief Van Der Waals gas model. */
  PR_GAS = 3,             /*!< \brief Perfect Real gas model. */
  CONSTANT_DENSITY = 4,   /*!< \brief Constant density gas model. */
  INC_IDEAL_GAS = 5,      /*!< \brief Incompressible ideal gas model. */
  INC_IDEAL_GAS_POLY = 6  /*!< \brief Inc. ideal gas, polynomial gas model. */
};
static const MapType<string, ENUM_FLUIDMODEL> FluidModel_Map = {
  MakePair("STANDARD_AIR", STANDARD_AIR)
  MakePair("IDEAL_GAS", IDEAL_GAS)
  MakePair("VW_GAS", VW_GAS)
  MakePair("PR_GAS", PR_GAS)
  MakePair("CONSTANT_DENSITY", CONSTANT_DENSITY)
  MakePair("INC_IDEAL_GAS", INC_IDEAL_GAS)
  MakePair("INC_IDEAL_GAS_POLY", INC_IDEAL_GAS_POLY)
};

/*!
 * \brief Types of density models
 */
enum ENUM_DENSITYMODEL {
  CONSTANT = 0,
  BOUSSINESQ = 1,  /*!< \brief BoussinesQ density model. */
  VARIABLE = 2     /*!< \brief Variable density model. */
};
static const MapType<string, ENUM_DENSITYMODEL> DensityModel_Map = {
  MakePair("CONSTANT", CONSTANT)
  MakePair("BOUSSINESQ", BOUSSINESQ)
  MakePair("VARIABLE", VARIABLE)
};

/*!
 * \brief Types of initialization option
 */
enum ENUM_INIT_OPTION {
  REYNOLDS = 0,      /*!< \brief Reynold's number initalization. */
  TD_CONDITIONS = 1  /*!< \brief Total conditions initalization. */
};
static const MapType<string, ENUM_INIT_OPTION> InitOption_Map = {
  MakePair("REYNOLDS", REYNOLDS)
  MakePair("TD_CONDITIONS", TD_CONDITIONS)
};

/*!
 * \brief Types of initialization option
 */
enum ENUM_FREESTREAM_OPTION {
  TEMPERATURE_FS = 0,  /*!< \brief Temperature initialization. */
  DENSITY_FS = 1       /*!< \brief Density initalization. */
};
static const MapType<string, ENUM_FREESTREAM_OPTION> FreeStreamOption_Map = {
  MakePair("TEMPERATURE_FS", TEMPERATURE_FS)
  MakePair("DENSITY_FS", DENSITY_FS)
};

/*!
 * \brief Types of viscosity model
 */
enum ENUM_VISCOSITYMODEL {
  CONSTANT_VISCOSITY = 0,   /*!< \brief Constant viscosity. */
  SUTHERLAND = 1,           /*!< \brief Sutherlands Law viscosity. */
  POLYNOMIAL_VISCOSITY = 2  /*!< \brief Polynomial viscosity. */
};
static const MapType<string, ENUM_VISCOSITYMODEL> ViscosityModel_Map = {
  MakePair("CONSTANT_VISCOSITY", CONSTANT_VISCOSITY)
  MakePair("SUTHERLAND", SUTHERLAND)
  MakePair("POLYNOMIAL_VISCOSITY", POLYNOMIAL_VISCOSITY)
};

/*!
 * \brief Types of thermal conductivity model
 */
enum ENUM_CONDUCTIVITYMODEL {
  CONSTANT_CONDUCTIVITY = 0,   /*!< \brief Constant thermal conductivity. */
  CONSTANT_PRANDTL = 1,        /*!< \brief Constant Prandtl number. */
  POLYNOMIAL_CONDUCTIVITY = 2  /*!< \brief Polynomial thermal conductivity. */
};
static const MapType<string, ENUM_CONDUCTIVITYMODEL> ConductivityModel_Map = {
  MakePair("CONSTANT_CONDUCTIVITY", CONSTANT_CONDUCTIVITY)
  MakePair("CONSTANT_PRANDTL", CONSTANT_PRANDTL)
  MakePair("POLYNOMIAL_CONDUCTIVITY", POLYNOMIAL_CONDUCTIVITY)
};

/*!
 * \brief Types of turbulent thermal conductivity model
 */
enum ENUM_CONDUCTIVITYMODEL_TURB {
  NO_CONDUCTIVITY_TURB  = 0,  /*!< \brief No turbulent contribution to the effective thermal conductivity for RANS. */
  CONSTANT_PRANDTL_TURB = 1   /*!< \brief Include contribution to effective conductivity using constant turbulent Prandtl number for RANS. */
};
static const MapType<string, ENUM_CONDUCTIVITYMODEL_TURB> TurbConductivityModel_Map = {
  MakePair("NONE", NO_CONDUCTIVITY_TURB)
  MakePair("CONSTANT_PRANDTL_TURB", CONSTANT_PRANDTL_TURB)
};

/*!
 * \brief Types of unsteady mesh motion
 */
enum ENUM_GRIDMOVEMENT {
  NO_MOVEMENT = 0,          /*!< \brief Simulation on a static mesh. */
  RIGID_MOTION = 2,         /*!< \brief Simulation with rigid mesh motion (plunging/pitching/rotation). */
  ROTATING_FRAME = 8,       /*!< \brief Simulation in a rotating frame. */
  ELASTICITY = 9,           /*!< \brief Linear Elasticity. */
  STEADY_TRANSLATION = 11,  /*!< \brief Simulation in a steadily translating frame. */
  GUST = 12,                /*!< \brief Simulation on a static mesh with a gust. */
  MOVING_HTP = 13,          /*!< \brief Simulation with moving HTP (rotation). */
};
static const MapType<string, ENUM_GRIDMOVEMENT> GridMovement_Map = {
  MakePair("NONE", NO_MOVEMENT)
  MakePair("RIGID_MOTION", RIGID_MOTION)
  MakePair("ROTATING_FRAME", ROTATING_FRAME)
  MakePair("ELASTICITY", ELASTICITY)
  MakePair("MOVING_HTP", MOVING_HTP)
  MakePair("STEADY_TRANSLATION", STEADY_TRANSLATION)
  MakePair("GUST", GUST)
};

enum ENUM_SURFACEMOVEMENT {
  DEFORMING = 1,                 /*!< \brief Simulation with deformation. */
  MOVING_WALL = 2,               /*!< \brief Simulation with moving wall. */
  AEROELASTIC = 3,               /*!< \brief Simulation with aeroelastic motion. */
  AEROELASTIC_RIGID_MOTION = 4,  /*!< \brief Simulation with rotation and aeroelastic motion. */
  FLUID_STRUCTURE = 5,           /*!< \brief Fluid structure deformation. */
  EXTERNAL = 6,                  /*!< \brief Simulation with external motion. */
  EXTERNAL_ROTATION = 7,         /*!< \brief Simulation with external rotation motion. */
};
static const MapType<string, ENUM_SURFACEMOVEMENT> SurfaceMovement_Map = {
  MakePair("DEFORMING", DEFORMING)
  MakePair("MOVING_WALL", MOVING_WALL)
  MakePair("AEROELASTIC_RIGID_MOTION", AEROELASTIC_RIGID_MOTION)
  MakePair("AEROELASTIC", AEROELASTIC)
  MakePair("FLUID_STRUCTURE", FLUID_STRUCTURE)
  MakePair("EXTERNAL", EXTERNAL)
  MakePair("EXTERNAL_ROTATION", EXTERNAL_ROTATION)
};

/*!
 * \brief Type of wind gusts
 */
enum ENUM_GUST_TYPE {
  NO_GUST = 0,      /*!< \brief No gust. */
  TOP_HAT = 1,      /*!< \brief Top-hat function shaped gust  */
  SINE = 2,         /*!< \brief Sine shaped gust */
  ONE_M_COSINE = 3, /*!< \brief 1-cosine shaped gust */
  VORTEX = 4,       /*!< \brief A gust made from vortices */
  EOG = 5           /*!< \brief An extreme operating gust */
};
static const MapType<string, ENUM_GUST_TYPE> Gust_Type_Map = {
  MakePair("NONE", NO_GUST)
  MakePair("TOP_HAT", TOP_HAT)
  MakePair("SINE", SINE)
  MakePair("ONE_M_COSINE", ONE_M_COSINE)
  MakePair("VORTEX", VORTEX)
  MakePair("EOG", EOG)
};

/*!
 * \brief Type of wind direction
 */
enum ENUM_GUST_DIR {
  X_DIR = 0,  /*!< \brief Gust direction-X. */
  Y_DIR = 1   /*!< \brief Gust direction-Y. */
};
static const MapType<string, ENUM_GUST_DIR> Gust_Dir_Map = {
  MakePair("X_DIR", X_DIR)
  MakePair("Y_DIR", Y_DIR)
};

// If you add to ENUM_CENTERED, you must also add the option to ENUM_CONVECTIVE
/*!
 * \brief Types of centered spatial discretizations
 */
enum ENUM_CENTERED {
  NO_CENTERED = 0,    /*!< \brief No centered scheme is used. */
  JST = 1,            /*!< \brief Jameson-Smith-Turkel centered numerical method. */
  LAX = 2,            /*!< \brief Lax-Friedrich centered numerical method. */
  JST_KE = 4          /*!< \brief Kinetic Energy preserving Jameson-Smith-Turkel centered numerical method. */
};
static const MapType<string, ENUM_CENTERED> Centered_Map = {
  MakePair("NONE", NO_CENTERED)
  MakePair("JST", JST)
  MakePair("JST_KE", JST_KE)
  MakePair("LAX-FRIEDRICH", LAX)
};


// If you add to ENUM_UPWIND, you must also add the option to ENUM_CONVECTIVE
/*!
 * \brief Types of upwind spatial discretizations
 */
enum ENUM_UPWIND {
  NO_UPWIND = 0,              /*!< \brief No upwind scheme is used. */
  ROE = 1,                    /*!< \brief Roe's upwind numerical method. */
  SCALAR_UPWIND = 2,          /*!< \brief Scalar upwind numerical method. */
  AUSM = 3,                   /*!< \brief AUSM numerical method. */
  HLLC = 4,                   /*!< \brief HLLC numerical method. */
  SW = 5,                     /*!< \brief Steger-Warming method. */
  MSW = 6,                    /*!< \brief Modified Steger-Warming method. */
  TURKEL = 7,                 /*!< \brief Roe-Turkel's upwind numerical method. */
  SLAU = 8,                   /*!< \brief Simple Low-Dissipation AUSM numerical method. */
  CUSP = 9,                   /*!< \brief Convective upwind and split pressure numerical method. */
  CONVECTIVE_TEMPLATE = 10,   /*!< \brief Template for new numerical method . */
  L2ROE = 11,                 /*!< \brief L2ROE numerical method . */
  LMROE = 12,                 /*!< \brief Rieper's Low Mach ROE numerical method . */
  SLAU2 = 13,                 /*!< \brief Simple Low-Dissipation AUSM 2 numerical method. */
  FDS = 14,                   /*!< \brief Flux difference splitting upwind method (incompressible flows). */
  LAX_FRIEDRICH = 15,         /*!< \brief Lax-Friedrich numerical method. */
  AUSMPLUSUP = 16,            /*!< \brief AUSM+ -up numerical method (All Speed) */
  AUSMPLUSUP2 = 17            /*!< \brief AUSM+ -up2 numerical method (All Speed) */
};
static const MapType<string, ENUM_UPWIND> Upwind_Map = {
  MakePair("NONE", NO_UPWIND)
  MakePair("ROE", ROE)
  MakePair("TURKEL_PREC", TURKEL)
  MakePair("AUSM", AUSM)
  MakePair("AUSMPLUSUP", AUSMPLUSUP)
  MakePair("AUSMPLUSUP2", AUSMPLUSUP2)
  MakePair("SLAU", SLAU)
  MakePair("HLLC", HLLC)
  MakePair("SW", SW)
  MakePair("MSW", MSW)
  MakePair("CUSP", CUSP)
  MakePair("SCALAR_UPWIND", SCALAR_UPWIND)
  MakePair("CONVECTIVE_TEMPLATE", CONVECTIVE_TEMPLATE)
  MakePair("L2ROE", L2ROE)
  MakePair("LMROE", LMROE)
  MakePair("SLAU2", SLAU2)
  MakePair("FDS", FDS)
  MakePair("LAX-FRIEDRICH", LAX_FRIEDRICH)
};

/*!
 * \brief Types of FEM spatial discretizations
 */
enum ENUM_FEM {
  NO_FEM = 0,  /*!< \brief No finite element scheme is used. */
  DG = 1       /*!< \brief Discontinuous Galerkin numerical method. */
};
static const MapType<string, ENUM_FEM> FEM_Map = {
  MakePair("NONE", NO_FEM)
  MakePair("DG", DG)
};

/*!
 * \brief Types of shock capturing method in Discontinuous Galerkin numerical method.
 */
enum ENUM_SHOCK_CAPTURING_DG {
  NO_SHOCK_CAPTURING = 0,     /*!< \brief Shock capturing is not used. */
  PERSSON = 1                 /*!< \brief Per-Olof Persson's sub-cell shock capturing method. */
};
static const MapType<string, ENUM_SHOCK_CAPTURING_DG> ShockCapturingDG_Map = {
  MakePair("NONE", NO_SHOCK_CAPTURING)
  MakePair("PERSSON", PERSSON)
};

/*!
 * \brief Types of matrix coloring to compute a sparse Jacobian matrix.
 */
enum ENUM_MATRIX_COLORING {
  GREEDY_COLORING = 0,            /*!< \brief Greedy type of algorithm for the coloring. */
  NATURAL_COLORING = 1            /*!< \brief One color for every DOF, very slow. Only to be used for debugging. */
};
static const MapType<string, ENUM_MATRIX_COLORING> MatrixColoring_Map = {
  MakePair("GREEDY_COLORING", GREEDY_COLORING)
  MakePair("NATURAL_COLORING", NATURAL_COLORING)
};

/*!
 * \brief Types of slope limiters
 */
enum ENUM_LIMITER {
  NO_LIMITER           = 0, /*!< \brief No limiter. */
  VENKATAKRISHNAN      = 1, /*!< \brief Slope limiter using Venkatakrisnan method (stencil formulation). */
  VENKATAKRISHNAN_WANG = 2, /*!< \brief Slope limiter using Venkatakrisnan method, eps based on solution (stencil formulation). */
  BARTH_JESPERSEN      = 3, /*!< \brief Slope limiter using Barth-Jespersen method (stencil formulation). */
  VAN_ALBADA_EDGE      = 4, /*!< \brief Slope limiter using Van Albada method (edge formulation). */
  SHARP_EDGES          = 5, /*!< \brief Slope limiter using sharp edges. */
  WALL_DISTANCE        = 6  /*!< \brief Slope limiter using wall distance. */
};
static const MapType<string, ENUM_LIMITER> Limiter_Map = {
  MakePair("NONE", NO_LIMITER)
  MakePair("VENKATAKRISHNAN", VENKATAKRISHNAN)
  MakePair("VENKATAKRISHNAN_WANG", VENKATAKRISHNAN_WANG)
  MakePair("BARTH_JESPERSEN", BARTH_JESPERSEN)
  MakePair("VAN_ALBADA_EDGE", VAN_ALBADA_EDGE)
  MakePair("SHARP_EDGES", SHARP_EDGES)
  MakePair("WALL_DISTANCE", WALL_DISTANCE)
};

/*!
 * \brief Types of turbulent models
 */
enum ENUM_TURB_MODEL {
  NO_TURB_MODEL = 0, /*!< \brief No turbulence model. */
  SA        = 1,     /*!< \brief Kind of Turbulent model (Spalart-Allmaras). */
  SA_NEG    = 2,     /*!< \brief Kind of Turbulent model (Spalart-Allmaras). */
  SA_E      = 3,     /*!< \brief Kind of Turbulent model (Spalart-Allmaras Edwards). */
  SA_COMP   = 4,     /*!< \brief Kind of Turbulent model (Spalart-Allmaras Compressibility Correction). */
  SA_E_COMP = 5,     /*!< \brief Kind of Turbulent model (Spalart-Allmaras Edwards with Compressibility Correction). */
  SST       = 6,     /*!< \brief Kind of Turbulence model (Menter SST). */
  SST_SUST  = 7      /*!< \brief Kind of Turbulence model (Menter SST with sustaining terms for free-stream preservation). */
};
static const MapType<string, ENUM_TURB_MODEL> Turb_Model_Map = {
  MakePair("NONE", NO_TURB_MODEL)
  MakePair("SA", SA)
  MakePair("SA_NEG", SA_NEG)
  MakePair("SA_E", SA_E)
  MakePair("SA_COMP", SA_COMP)
  MakePair("SA_E_COMP", SA_E_COMP)
  MakePair("SST", SST)
  MakePair("SST_SUST", SST_SUST)
};

/*!
 * \brief Types of transition models
 */
enum ENUM_TRANS_MODEL {
  NO_TRANS_MODEL = 0,  /*!< \brief No transition model. */
  LM = 1,              /*!< \brief Kind of transition model (Langtry-Menter (LM) for SST and Spalart-Allmaras). */
  BC = 2               /*!< \brief Kind of transition model (BAS-CAKMAKCIOGLU (BC) for Spalart-Allmaras). */
};
static const MapType<string, ENUM_TRANS_MODEL> Trans_Model_Map = {
  MakePair("NONE", NO_TRANS_MODEL)
  MakePair("LM", LM)
  MakePair("BC", BC)
};

/*!
 * \brief Types of subgrid scale models
 */
enum ENUM_SGS_MODEL {
  NO_SGS_MODEL = 0, /*!< \brief No subgrid scale model. */
  IMPLICIT_LES = 1, /*!< \brief Implicit LES, i.e. no explicit SGS model. */
  SMAGORINSKY  = 2, /*!< \brief Smagorinsky SGS model. */
  WALE         = 3, /*!< \brief Wall-Adapting Local Eddy-viscosity SGS model. */
  VREMAN       = 4  /*!< \brief Vreman SGS model. */
};
static const MapType<string, ENUM_SGS_MODEL> SGS_Model_Map = {
  MakePair("NONE",         NO_SGS_MODEL)
  MakePair("IMPLICIT_LES", IMPLICIT_LES)
  MakePair("SMAGORINSKY",  SMAGORINSKY)
  MakePair("WALE",         WALE)
  MakePair("VREMAN",       VREMAN)
};


/*!
 * \brief Types of window (weight) functions for cost functional
 */
enum WINDOW_FUNCTION {
  SQUARE = 0,        /*!< \brief No weight function  (order 1)*/
  HANN = 1,          /*!< \brief Hann-type weight function (order 3) */
  HANN_SQUARE = 2,   /*!< \brief Hann-squared type weight function (order 5)*/
  BUMP = 3,          /*!< \brief bump type weight function (exponential order of convergence) */
};
static const MapType<string, WINDOW_FUNCTION> Window_Map = {
  MakePair("SQUARE", SQUARE)
  MakePair("HANN", HANN)
  MakePair("HANN_SQUARE", HANN_SQUARE)
  MakePair("BUMP", BUMP)
};

/*!
 * \brief Types of hybrid RANS/LES models
 */
enum ENUM_HYBRIDRANSLES {
  NO_HYBRIDRANSLES = 0,  /*!< \brief No turbulence model. */
  SA_DES   = 1,          /*!< \brief Kind of Hybrid RANS/LES (SA - Detached Eddy Simulation (DES)). */
  SA_DDES  = 2,          /*!< \brief Kind of Hybrid RANS/LES (SA - Delayed DES (DDES) with Delta_max SGS ). */
  SA_ZDES  = 3,          /*!< \brief Kind of Hybrid RANS/LES (SA - Delayed DES (DDES) with Vorticity based SGS like Zonal DES). */
  SA_EDDES = 4           /*!< \brief Kind of Hybrid RANS/LES (SA - Delayed DES (DDES) with Shear Layer Adapted SGS: Enhanced DDES). */
};
static const MapType<string, ENUM_HYBRIDRANSLES> HybridRANSLES_Map = {
  MakePair("NONE", NO_HYBRIDRANSLES)
  MakePair("SA_DES", SA_DES)
  MakePair("SA_DDES", SA_DDES)
  MakePair("SA_ZDES", SA_ZDES)
  MakePair("SA_EDDES", SA_EDDES)
};

/*!
 * \brief Types of Roe Low Dissipation Schemes
 */
enum ENUM_ROELOWDISS {
    NO_ROELOWDISS = 0, /*!< \brief No Roe Low Dissipation model. */
    FD            = 1, /*!< \brief Numerical Blending based on DDES's F_d function */
    NTS           = 2, /*!< \brief Numerical Blending of Travin and Shur. */
    NTS_DUCROS    = 3, /*!< \brief Numerical Blending of Travin and Shur + Ducros' Shock Sensor. */
    FD_DUCROS     = 4  /*!< \brief Numerical Blending based on DDES's F_d function + Ducros' Shock Sensor */
};
static const MapType<string, ENUM_ROELOWDISS> RoeLowDiss_Map = {
  MakePair("NONE", NO_ROELOWDISS)
  MakePair("FD", FD)
  MakePair("NTS", NTS)
  MakePair("NTS_DUCROS", NTS_DUCROS)
  MakePair("FD_DUCROS", FD_DUCROS)
};

/*!
 * \brief Types of wall functions.
 */
enum ENUM_WALL_FUNCTIONS {
  NO_WALL_FUNCTION          = 0,   /*!< \brief No wall function treatment, integration to the wall. Default behavior. */
  STANDARD_WALL_FUNCTION    = 1,   /*!< \brief Standard wall function. */
  ADAPTIVE_WALL_FUNCTION    = 2,   /*!< \brief Adaptive wall function. Formulation depends on y+. */
  SCALABLE_WALL_FUNCTION    = 3,   /*!< \brief Scalable wall function. */
  EQUILIBRIUM_WALL_MODEL    = 4,   /*!< \brief Equilibrium wall model for LES. */
  NONEQUILIBRIUM_WALL_MODEL = 5,   /*!< \brief Non-equilibrium wall model for LES. */
  LOGARITHMIC_WALL_MODEL    = 6    /*!< \brief Logarithmic law-of-the-wall model for LES. */
};
static const MapType<string, ENUM_WALL_FUNCTIONS> Wall_Functions_Map = {
  MakePair("NO_WALL_FUNCTION",          NO_WALL_FUNCTION)
  MakePair("STANDARD_WALL_FUNCTION",    STANDARD_WALL_FUNCTION)
  MakePair("ADAPTIVE_WALL_FUNCTION",    ADAPTIVE_WALL_FUNCTION)
  MakePair("SCALABLE_WALL_FUNCTION",    SCALABLE_WALL_FUNCTION)
  MakePair("EQUILIBRIUM_WALL_MODEL",    EQUILIBRIUM_WALL_MODEL)
  MakePair("NONEQUILIBRIUM_WALL_MODEL", NONEQUILIBRIUM_WALL_MODEL)
  MakePair("LOGARITHMIC_WALL_MODEL", LOGARITHMIC_WALL_MODEL)
};

/*!
 * \brief Type of time integration schemes
 */
enum ENUM_TIME_INT {
  RUNGE_KUTTA_EXPLICIT = 1,   /*!< \brief Explicit Runge-Kutta time integration definition. */
  EULER_EXPLICIT = 2,         /*!< \brief Explicit Euler time integration definition. */
  EULER_IMPLICIT = 3,         /*!< \brief Implicit Euler time integration definition. */
  CLASSICAL_RK4_EXPLICIT = 4, /*!< \brief Classical RK4 time integration definition. */
  ADER_DG = 5                 /*!< \brief ADER-DG time integration definition. */
};
static const MapType<string, ENUM_TIME_INT> Time_Int_Map = {
  MakePair("RUNGE-KUTTA_EXPLICIT", RUNGE_KUTTA_EXPLICIT)
  MakePair("EULER_EXPLICIT", EULER_EXPLICIT)
  MakePair("EULER_IMPLICIT", EULER_IMPLICIT)
  MakePair("CLASSICAL_RK4_EXPLICIT", CLASSICAL_RK4_EXPLICIT)
  MakePair("ADER_DG", ADER_DG)
};

/*!
 * \brief Type of predictor for the ADER-DG time integration scheme.
 */
enum ENUM_ADER_PREDICTOR {
  ADER_ALIASED_PREDICTOR     = 1, /*!< \brief Aliased predictor, easiest to do. */
  ADER_NON_ALIASED_PREDICTOR = 2  /*!< \brief Non-aliased predictor. Consistent, but more difficult. */
};
static const MapType<string, ENUM_ADER_PREDICTOR> Ader_Predictor_Map = {
  MakePair("ADER_ALIASED_PREDICTOR", ADER_ALIASED_PREDICTOR)
  MakePair("ADER_NON_ALIASED_PREDICTOR", ADER_NON_ALIASED_PREDICTOR)
};

/*!
 * \brief Type of heat timestep calculation
 */
enum ENUM_HEAT_TIMESTEP {
  MINIMUM = 1,     /*!< \brief Local time stepping based on minimum lambda.*/
  CONVECTIVE = 2,  /*!< \brief Local time stepping based on convective spectral radius.*/
  VISCOUS = 3,     /*!< \brief Local time stepping based on viscous spectral radius.*/
  BYFLOW = 4,      /*!< \brief Unsing the mean solvers time step. */
};
static const MapType<string, ENUM_HEAT_TIMESTEP> Heat_TimeStep_Map = {
  MakePair("LOCAL", MINIMUM)
  MakePair("CONVECTIVE", CONVECTIVE)
  MakePair("VISCOUS", VISCOUS)
  MakePair("BYFLOW", BYFLOW)
};

/*!
 * \brief Type of time integration schemes
 */
enum ENUM_TIME_INT_FEA {
  CD_EXPLICIT = 1,       /*!< \brief Support for implementing an explicit method. */
  NEWMARK_IMPLICIT = 2,  /*!< \brief Implicit Newmark integration definition. */
  GENERALIZED_ALPHA = 3  /*!< \brief Support for implementing another implicit method. */
};
static const MapType<string, ENUM_TIME_INT_FEA> Time_Int_Map_FEA = {
  MakePair("CD_EXPLICIT", CD_EXPLICIT)
  MakePair("NEWMARK_IMPLICIT", NEWMARK_IMPLICIT)
  MakePair("GENERALIZED_ALPHA", GENERALIZED_ALPHA)
};

/*!
 * \brief Type of time integration schemes
 */
enum ENUM_SPACE_ITE_FEA {
  NEWTON_RAPHSON = 1,           /*!< \brief Full Newton-Rapshon method. */
  MODIFIED_NEWTON_RAPHSON = 2   /*!< \brief Modified Newton-Raphson method. */
};
static const MapType<string, ENUM_SPACE_ITE_FEA> Space_Ite_Map_FEA = {
  MakePair("NEWTON_RAPHSON", NEWTON_RAPHSON)
  MakePair("MODIFIED_NEWTON_RAPHSON", MODIFIED_NEWTON_RAPHSON)
};

/*!
 * \brief Types of schemes to compute the flow gradient
 */
enum ENUM_FLOW_GRADIENT {
  NO_GRADIENT            = 0,   /*!< \brief No gradient method. Only possible for reconstruction gradient, in which case, the option chosen for NUM_METHOD_GRAD is used. */
  GREEN_GAUSS            = 1,   /*!< \brief Gradient computation using Green-Gauss theorem. */
  LEAST_SQUARES          = 2,   /*!< \brief Gradient computation using unweighted least squares. */
  WEIGHTED_LEAST_SQUARES = 3    /*!< \brief Gradients computation using inverse-distance weighted least squares. */
};
static const MapType<string, ENUM_FLOW_GRADIENT> Gradient_Map = {
  MakePair("NONE", NO_GRADIENT)
  MakePair("GREEN_GAUSS", GREEN_GAUSS)
  MakePair("LEAST_SQUARES", LEAST_SQUARES)
  MakePair("WEIGHTED_LEAST_SQUARES", WEIGHTED_LEAST_SQUARES)
};

/*!
 * \brief Types of action to take on a geometry structure
 */
enum GEOMETRY_ACTION {
  ALLOCATE = 0,     /*!< \brief Allocate geometry structure. */
  UPDATE = 1        /*!< \brief Update geometry structure (grid moving, adaptation, etc.). */
};

/*!
 * \brief Types of action to perform when doing the geometry evaluation
 */
enum GEOMETRY_MODE {
  FUNCTION = 0,     /*!< \brief Geometrical analysis. */
  GRADIENT = 1      /*!< \brief Geometrical analysis and gradient using finite differences. */
};
static const MapType<string, GEOMETRY_MODE> GeometryMode_Map = {
  MakePair("FUNCTION", FUNCTION)
  MakePair("GRADIENT", GRADIENT)
};

/*!
 * \brief Types of boundary conditions
 */
enum BC_TYPE {
  EULER_WALL = 1,             /*!< \brief Boundary Euler wall definition. */
  FAR_FIELD = 2,              /*!< \brief Boundary far-field definition. */
  SYMMETRY_PLANE = 3,         /*!< \brief Boundary symmetry plane definition. */
  INLET_FLOW = 4,             /*!< \brief Boundary inlet flow definition. */
  OUTLET_FLOW = 5,            /*!< \brief Boundary outlet flow definition. */
  PERIODIC_BOUNDARY = 6,      /*!< \brief Periodic boundary definition. */
  NEARFIELD_BOUNDARY = 7,     /*!< \brief Near-Field boundary definition. */
  ELECTRODE_BOUNDARY = 8,     /*!< \brief Electrode boundary definition. */
  DIELEC_BOUNDARY = 9,        /*!< \brief Dipoisson boundary definition. */
  CUSTOM_BOUNDARY = 10,       /*!< \brief custom boundary definition. */
  INTERFACE_BOUNDARY = 11,    /*!< \brief Domain interface boundary definition. */
  DIRICHLET = 12,             /*!< \brief Boundary Euler wall definition. */
  NEUMANN = 13,               /*!< \brief Boundary Neumann definition. */
  DISPLACEMENT_BOUNDARY = 14, /*!< \brief Boundary displacement definition. */
  LOAD_BOUNDARY = 15,         /*!< \brief Boundary Load definition. */
  FLOWLOAD_BOUNDARY = 16,     /*!< \brief Boundary Load definition. */
  SUPERSONIC_INLET = 19,      /*!< \brief Boundary supersonic inlet definition. */
  SUPERSONIC_OUTLET = 20,     /*!< \brief Boundary supersonic inlet definition. */
  ENGINE_INFLOW = 21,         /*!< \brief Boundary nacelle inflow. */
  ENGINE_EXHAUST = 22,        /*!< \brief Boundary nacelle exhaust. */
  RIEMANN_BOUNDARY= 24,       /*!< \brief Riemann Boundary definition. */
  ISOTHERMAL = 25,            /*!< \brief No slip isothermal wall boundary condition. */
  HEAT_FLUX  = 26,            /*!< \brief No slip constant heat flux wall boundary condition. */
  ACTDISK_INLET = 32,         /*!< \brief Actuator disk inlet boundary definition. */
  ACTDISK_OUTLET = 33,        /*!< \brief Actuator disk outlet boundary definition. */
  CLAMPED_BOUNDARY = 34,      /*!< \brief Clamped Boundary definition. */
  LOAD_DIR_BOUNDARY = 35,     /*!< \brief Boundary Load definition. */
  LOAD_SINE_BOUNDARY = 36,    /*!< \brief Sine-waveBoundary Load definition. */
  GILES_BOUNDARY= 37,         /*!< \brief Giles Boundary definition. */
  INTERNAL_BOUNDARY= 38,      /*!< \brief Internal Boundary definition. */
  FLUID_INTERFACE = 39,       /*!< \brief Domain interface definition. */
  DISP_DIR_BOUNDARY = 40,     /*!< \brief Boundary displacement definition. */
  DAMPER_BOUNDARY = 41,       /*!< \brief Damper. */
  CHT_WALL_INTERFACE = 50,    /*!< \brief Domain interface definition. */
  SEND_RECEIVE = 99,          /*!< \brief Boundary send-receive definition. */
};

/*!
 * \brief Different regime modes
 */
enum ENUM_2DFORM {
  PLANE_STRESS = 0,     /*!< \brief Definition of plane stress solver. */
  PLANE_STRAIN = 1      /*!< \brief Definition of plane strain solver. */
};
static const MapType<string, ENUM_2DFORM> ElasForm_2D = {
  MakePair("PLANE_STRESS", PLANE_STRESS)
  MakePair("PLANE_STRAIN", PLANE_STRAIN)
};

/*!
 * \brief Kinds of relaxation for FSI problem
 */
enum ENUM_AITKEN {
  NO_RELAXATION = 0,        /*!< \brief No relaxation in the strongly coupled approach. */
  FIXED_PARAMETER = 1,      /*!< \brief Relaxation with a fixed parameter. */
  AITKEN_DYNAMIC = 2        /*!< \brief Relaxation using Aitken's dynamic parameter. */
};
static const MapType<string, ENUM_AITKEN> AitkenForm_Map = {
  MakePair("NONE", NO_RELAXATION)
  MakePair("FIXED_PARAMETER", FIXED_PARAMETER)
  MakePair("AITKEN_DYNAMIC", AITKEN_DYNAMIC)
};

/*!
 * \brief Types of dynamic transfer methods
 */
enum ENUM_DYN_TRANSFER_METHOD {
  INSTANTANEOUS = 1,   /*!< \brief No ramp, load is transfer instantaneously. */
  POL_ORDER_1 = 2,     /*!< \brief The load is transferred using a ramp. */
  POL_ORDER_3 = 3,     /*!< \brief The load is transferred using an order 3 polynomial function */
  POL_ORDER_5 = 4,     /*!< \brief The load is transferred using an order 5 polynomial function */
  SIGMOID_10 = 5,      /*!< \brief The load is transferred using a sigmoid with parameter 10 */
  SIGMOID_20 = 6       /*!< \brief The load is transferred using a sigmoid with parameter 20 */
};
static const MapType<string, ENUM_DYN_TRANSFER_METHOD> Dyn_Transfer_Method_Map = {
  MakePair("INSTANTANEOUS", INSTANTANEOUS)
  MakePair("RAMP", POL_ORDER_1)
  MakePair("CUBIC", POL_ORDER_3)
  MakePair("QUINTIC", POL_ORDER_5)
  MakePair("SIGMOID_10", SIGMOID_10)
  MakePair("SIGMOID_20", SIGMOID_20)
};

/*!
 * \brief Kinds of Design Variables for FEA problems
 */
enum ENUM_DVFEA {
  NODV_FEA = 0,         /*!< \brief No design variable for FEA problems. */
  YOUNG_MODULUS = 1,    /*!< \brief Young modulus (E) as design variable. */
  POISSON_RATIO = 2,    /*!< \brief Poisson ratio (Nu) as design variable. */
  DENSITY_VAL = 3,      /*!< \brief Density (Rho) as design variable. */
  DEAD_WEIGHT = 4,      /*!< \brief Dead Weight (Rho_DL) as design variable. */
  ELECTRIC_FIELD = 5    /*!< \brief Electric field (E) as design variable. */
};
static const MapType<string, ENUM_DVFEA> DVFEA_Map = {
  MakePair("NONE", NODV_FEA)
  MakePair("YOUNG_MODULUS", YOUNG_MODULUS)
  MakePair("POISSON_RATIO", POISSON_RATIO)
  MakePair("DENSITY", DENSITY_VAL)
  MakePair("DEAD_WEIGHT", DEAD_WEIGHT)
  MakePair("ELECTRIC_FIELD", ELECTRIC_FIELD)
};

/*!
 * \brief Kinds of radiation models
 */
enum ENUM_RADIATION {
  NO_RADIATION = 0,      /*!< \brief No radiation model */
  P1_MODEL = 1           /*!< \brief P1 Radiation model. */
};
static const MapType<string, ENUM_RADIATION> Radiation_Map = {
  MakePair("NONE", NO_RADIATION)
  MakePair("P1", P1_MODEL)
};

/*!
 * \brief Kinds of P1 initialization
 */
enum ENUM_P1_INIT {
  P1_INIT_ZERO = 0,      /*!< \brief Initialize the P1 model from zero values */
  P1_INIT_TEMP = 1       /*!< \brief Initialize the P1 model from blackbody energy computed from the initial temperature. */
};
static const MapType<string, ENUM_P1_INIT> P1_Init_Map = {
  MakePair("ZERO", P1_INIT_ZERO)
  MakePair("TEMPERATURE_INIT", P1_INIT_TEMP)
};

/*!
 * \brief Kinds of coupling methods at CHT interfaces.
 * The first (temperature) part determines the BC method on the fluid side, the second (heatflux) part determines
 * the BC method on the solid side of the CHT interface.
 */
enum ENUM_CHT_COUPLING {
  DIRECT_TEMPERATURE_NEUMANN_HEATFLUX = 0,
  AVERAGED_TEMPERATURE_NEUMANN_HEATFLUX = 1,
  DIRECT_TEMPERATURE_ROBIN_HEATFLUX = 2,
  AVERAGED_TEMPERATURE_ROBIN_HEATFLUX = 3
};
static const MapType<string, ENUM_CHT_COUPLING> CHT_Coupling_Map = {
  MakePair("DIRECT_TEMPERATURE_NEUMANN_HEATFLUX", DIRECT_TEMPERATURE_NEUMANN_HEATFLUX)
  MakePair("AVERAGED_TEMPERATURE_NEUMANN_HEATFLUX", AVERAGED_TEMPERATURE_NEUMANN_HEATFLUX)
  MakePair("DIRECT_TEMPERATURE_ROBIN_HEATFLUX", DIRECT_TEMPERATURE_ROBIN_HEATFLUX)
  MakePair("AVERAGED_TEMPERATURE_ROBIN_HEATFLUX", AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)
};

/*!
 * \brief Types Riemann boundary treatments
 */
enum RIEMANN_TYPE {
  TOTAL_CONDITIONS_PT = 1,          /*!< \brief User specifies total pressure, total temperature, and flow direction. */
  DENSITY_VELOCITY = 2,             /*!< \brief User specifies density and velocity, and flow direction. */
  STATIC_PRESSURE = 3,              /*!< \brief User specifies static pressure. */
  TOTAL_SUPERSONIC_INFLOW = 4,	    /*!< \brief User specifies total pressure, total temperature and Velocity components. */
  STATIC_SUPERSONIC_INFLOW_PT = 5,  /*!< \brief User specifies static pressure, static temperature, and Mach components. */
  STATIC_SUPERSONIC_INFLOW_PD = 6,  /*!< \brief User specifies static pressure, static temperature, and Mach components. */
  MIXING_IN = 7,                    /*!< \brief User does not specify anything; information is retrieved from the other domain */
  MIXING_OUT = 8,                   /*!< \brief User does not specify anything; information is retrieved from the other domain */
  SUPERSONIC_OUTFLOW = 9,
  RADIAL_EQUILIBRIUM = 10,
  TOTAL_CONDITIONS_PT_1D = 11,
  STATIC_PRESSURE_1D = 12,
  MIXING_IN_1D = 13,
  MIXING_OUT_1D =14
};
static const MapType<string, RIEMANN_TYPE> Riemann_Map = {
  MakePair("TOTAL_CONDITIONS_PT", TOTAL_CONDITIONS_PT)
  MakePair("DENSITY_VELOCITY", DENSITY_VELOCITY)
  MakePair("STATIC_PRESSURE", STATIC_PRESSURE)
  MakePair("TOTAL_SUPERSONIC_INFLOW", TOTAL_SUPERSONIC_INFLOW)
  MakePair("STATIC_SUPERSONIC_INFLOW_PT", STATIC_SUPERSONIC_INFLOW_PT)
  MakePair("STATIC_SUPERSONIC_INFLOW_PD", STATIC_SUPERSONIC_INFLOW_PD)
  MakePair("MIXING_IN", MIXING_IN)
  MakePair("MIXING_OUT", MIXING_OUT)
  MakePair("MIXING_IN_1D", MIXING_IN_1D)
  MakePair("MIXING_OUT_1D", MIXING_OUT_1D)
  MakePair("SUPERSONIC_OUTFLOW", SUPERSONIC_OUTFLOW)
  MakePair("RADIAL_EQUILIBRIUM", RADIAL_EQUILIBRIUM)
  MakePair("TOTAL_CONDITIONS_PT_1D", TOTAL_CONDITIONS_PT_1D)
  MakePair("STATIC_PRESSURE_1D", STATIC_PRESSURE_1D)
};

static const MapType<string, RIEMANN_TYPE> Giles_Map = {
  MakePair("TOTAL_CONDITIONS_PT", TOTAL_CONDITIONS_PT)
  MakePair("DENSITY_VELOCITY", DENSITY_VELOCITY)
  MakePair("STATIC_PRESSURE", STATIC_PRESSURE)
  MakePair("TOTAL_SUPERSONIC_INFLOW", TOTAL_SUPERSONIC_INFLOW)
  MakePair("STATIC_SUPERSONIC_INFLOW_PT", STATIC_SUPERSONIC_INFLOW_PT)
  MakePair("STATIC_SUPERSONIC_INFLOW_PD", STATIC_SUPERSONIC_INFLOW_PD)
  MakePair("MIXING_IN", MIXING_IN)
  MakePair("MIXING_OUT", MIXING_OUT)
  MakePair("MIXING_IN_1D", MIXING_IN_1D)
  MakePair("MIXING_OUT_1D", MIXING_OUT_1D)
  MakePair("SUPERSONIC_OUTFLOW", SUPERSONIC_OUTFLOW)
  MakePair("RADIAL_EQUILIBRIUM", RADIAL_EQUILIBRIUM)
  MakePair("TOTAL_CONDITIONS_PT_1D", TOTAL_CONDITIONS_PT_1D)
  MakePair("STATIC_PRESSURE_1D", STATIC_PRESSURE_1D)
};

/*!
 * \brief Types of mixing process for averaging quantities at the boundaries.
 */
enum AVERAGEPROCESS_TYPE {
  ALGEBRAIC = 1,  /*!< \brief an algebraic average is computed at the boundary of interest. */
  AREA = 2,       /*!< \brief an area average is computed at the boundary of interest. */
  MIXEDOUT = 3,	  /*!< \brief an mixed-out average is computed at the boundary of interest. */
  MASSFLUX = 4    /*!< \brief a mass flow average is computed at the boundary of interest. */
};
static const MapType<string, AVERAGEPROCESS_TYPE> AverageProcess_Map = {
  MakePair("ALGEBRAIC", ALGEBRAIC)
  MakePair("AREA", AREA)
  MakePair("MIXEDOUT", MIXEDOUT)
  MakePair("MASSFLUX", MASSFLUX)
};

/*!
 * \brief Types of mixing process for averaging quantities at the boundaries.
 */
enum MIXINGPLANE_INTERFACE_TYPE {
  MATCHING = 1,             /*!< \brief an algebraic average is computed at the boundary of interest. */
  NEAREST_SPAN = 2,         /*!< \brief an area average is computed at the boundary of interest. */
  LINEAR_INTERPOLATION = 3  /*!< \brief an mixed-out average is computed at the boundary of interest. */
};
static const MapType<string, MIXINGPLANE_INTERFACE_TYPE> MixingPlaneInterface_Map = {
  MakePair("MATCHING", MATCHING)
  MakePair("NEAREST_SPAN",  NEAREST_SPAN)
  MakePair("LINEAR_INTERPOLATION", LINEAR_INTERPOLATION)
};

/*!
 * \brief this option allow to compute the span-wise section in different ways.
 */
enum SPANWISE_TYPE {
  AUTOMATIC = 1,      /*!< \brief number of span-wise section are computed automatically */
  EQUISPACED = 2      /*!< \brief number of span-wise section are specified from the user */
};
static const MapType<string, SPANWISE_TYPE> SpanWise_Map = {
  MakePair("AUTOMATIC", AUTOMATIC)
  MakePair("EQUISPACED", EQUISPACED)
};

/*!
 * \brief Types of mixing process for averaging quantities at the boundaries.
 */
enum TURBOMACHINERY_TYPE {
  AXIAL = 1,              /*!< \brief axial turbomachinery. */
  CENTRIFUGAL = 2,        /*!< \brief centrifugal turbomachinery. */
  CENTRIPETAL = 3,        /*!< \brief centripetal turbomachinery. */
  CENTRIPETAL_AXIAL = 4,  /*!< \brief mixed flow turbine. */
  AXIAL_CENTRIFUGAL = 5   /*!< \brief mixed flow turbine. */
};
static const MapType<string, TURBOMACHINERY_TYPE> TurboMachinery_Map = {
  MakePair("AXIAL", AXIAL)
  MakePair("CENTRIFUGAL", CENTRIFUGAL)
  MakePair("CENTRIPETAL",  CENTRIPETAL)
  MakePair("CENTRIPETAL_AXIAL",  CENTRIPETAL_AXIAL)
  MakePair("AXIAL_CENTRIFUGAL",  AXIAL_CENTRIFUGAL)
};

/*!
 * \brief Types of Turbomachinery performance flag.
 */
enum TURBO_MARKER_TYPE{
  INFLOW   = 1,	  /*!< \brief flag for inflow marker for compute turboperformance. */
  OUTFLOW = 2     /*!< \brief flag for outflow marker for compute turboperformance. */
};

/*!
 * \brief Types inlet boundary treatments
 */
enum INLET_TYPE {
  TOTAL_CONDITIONS = 1,	  /*!< \brief User specifies total pressure, total temperature, and flow direction. */
  MASS_FLOW = 2,          /*!< \brief User specifies density and velocity (mass flow). */
  INPUT_FILE = 3,         /*!< \brief User specifies an input file. */
  VELOCITY_INLET = 4,     /*!< \brief Velocity inlet for an incompressible flow. */
  PRESSURE_INLET = 5      /*!< \brief Total pressure inlet for an incompressible flow. */
};
static const MapType<string, INLET_TYPE> Inlet_Map = {
  MakePair("TOTAL_CONDITIONS", TOTAL_CONDITIONS)
  MakePair("MASS_FLOW", MASS_FLOW)
  MakePair("INPUT_FILE", INPUT_FILE)
  MakePair("VELOCITY_INLET", VELOCITY_INLET)
  MakePair("PRESSURE_INLET", PRESSURE_INLET)
};

/*!
 * \brief Types outlet boundary treatments
 */
enum OUTLET_TYPE {
  PRESSURE_OUTLET = 1,    /*!< \brief Gauge pressure outlet for incompressible flow */
  MASS_FLOW_OUTLET = 2,   /*!< \brief Mass flow outlet for incompressible flow. */
};
static const MapType<string, OUTLET_TYPE> Outlet_Map = {
  MakePair("PRESSURE_OUTLET", PRESSURE_OUTLET)
  MakePair("MASS_FLOW_OUTLET", MASS_FLOW_OUTLET)
};

/*!
 * \brief Types engine inflow boundary treatments
 */
enum ENGINE_INFLOW_TYPE {
  FAN_FACE_MACH = 1,          /*!< \brief User specifies fan face mach number. */
  FAN_FACE_MDOT = 2,          /*!< \brief User specifies Static pressure. */
  FAN_FACE_PRESSURE = 3       /*!< \brief User specifies Static pressure. */
};
static const MapType<string, ENGINE_INFLOW_TYPE> Engine_Inflow_Map = {
  MakePair("FAN_FACE_MACH", FAN_FACE_MACH)
  MakePair("FAN_FACE_MDOT", FAN_FACE_MDOT)
  MakePair("FAN_FACE_PRESSURE", FAN_FACE_PRESSURE)
};

/*!
 * \brief Types actuator disk boundary treatments
 */
enum ACTDISK_TYPE {
  VARIABLES_JUMP = 1,     /*!< \brief User specifies the variables jump. */
  BC_THRUST = 2,          /*!< \brief User specifies the BC thrust. */
  NET_THRUST = 3,         /*!< \brief User specifies the Net thrust. */
  DRAG_MINUS_THRUST = 4,  /*!< \brief User specifies the D-T. */
  MASSFLOW = 5,           /*!< \brief User specifies the massflow. */
  POWER = 6               /*!< \brief User specifies the power. */
};
static const MapType<string, ACTDISK_TYPE> ActDisk_Map = {
  MakePair("VARIABLES_JUMP", VARIABLES_JUMP)
  MakePair("BC_THRUST", BC_THRUST)
  MakePair("NET_THRUST", NET_THRUST)
  MakePair("DRAG_MINUS_THRUST", DRAG_MINUS_THRUST)
  MakePair("MASSFLOW", MASSFLOW)
  MakePair("POWER", POWER)
};

/*!
 * \brief Types of geometric entities based on VTK nomenclature
 */
enum GEO_TYPE {
  VERTEX = 1,         /*!< \brief VTK nomenclature for defining a vertex element. */
  LINE = 3,           /*!< \brief VTK nomenclature for defining a line element. */
  TRIANGLE = 5,       /*!< \brief VTK nomenclature for defining a triangle element. */
  QUADRILATERAL = 9,  /*!< \brief VTK nomenclature for defining a quadrilateral element. */
  TETRAHEDRON = 10,   /*!< \brief VTK nomenclature for defining a tetrahedron element. */
  HEXAHEDRON = 12,    /*!< \brief VTK nomenclature for defining a hexahedron element. */
  PRISM = 13,         /*!< \brief VTK nomenclature for defining a prism element. */
  PYRAMID = 14        /*!< \brief VTK nomenclature for defining a pyramid element. */
};

/*!
 * \brief Types of objective functions
 */
enum ENUM_OBJECTIVE {
  DRAG_COEFFICIENT = 1,         /*!< \brief Drag objective function definition. */
  LIFT_COEFFICIENT = 2,         /*!< \brief Lift objective function definition. */
  SIDEFORCE_COEFFICIENT = 3,    /*!< \brief Side force objective function definition. */
  EFFICIENCY = 4,               /*!< \brief Efficiency objective function definition. */
  INVERSE_DESIGN_PRESSURE = 5,	/*!< \brief Pressure objective function definition (inverse design). */
  INVERSE_DESIGN_HEATFLUX = 6,  /*!< \brief Heat flux objective function definition (inverse design). */
  TOTAL_HEATFLUX = 7,           /*!< \brief Total heat flux. */
  MAXIMUM_HEATFLUX = 8,         /*!< \brief Maximum heat flux. */
  TOTAL_AVG_TEMPERATURE = 70,   /*!< \brief Total averaged temperature. */
  MOMENT_X_COEFFICIENT = 9,     /*!< \brief Pitching moment objective function definition. */
  MOMENT_Y_COEFFICIENT = 10,    /*!< \brief Rolling moment objective function definition. */
  MOMENT_Z_COEFFICIENT = 11,    /*!< \brief Yawing objective function definition. */
  EQUIVALENT_AREA = 12,         /*!< \brief Equivalent area objective function definition. */
  NEARFIELD_PRESSURE = 13,      /*!< \brief NearField Pressure objective function definition. */
  FORCE_X_COEFFICIENT = 14,     /*!< \brief X-direction force objective function definition. */
  FORCE_Y_COEFFICIENT = 15,     /*!< \brief Y-direction force objective function definition. */
  FORCE_Z_COEFFICIENT = 16,     /*!< \brief Z-direction force objective function definition. */
  THRUST_COEFFICIENT = 17,      /*!< \brief Thrust objective function definition. */
  TORQUE_COEFFICIENT = 18,      /*!< \brief Torque objective function definition. */
  FIGURE_OF_MERIT = 19,         /*!< \brief Rotor Figure of Merit objective function definition. */
  BUFFET_SENSOR = 20,           /*!< \brief Sensor for detecting separation. */
  SURFACE_TOTAL_PRESSURE = 28,  /*!< \brief Total Pressure objective function definition. */
  SURFACE_STATIC_PRESSURE = 29, /*!< \brief Static Pressure objective function definition. */
  SURFACE_MASSFLOW = 30,        /*!< \brief Mass Flow Rate objective function definition. */
  SURFACE_MACH = 51,            /*!< \brief Mach number objective function definition. */
  SURFACE_UNIFORMITY = 52,      /*!< \brief Flow uniformity objective function definition. */
  SURFACE_SECONDARY = 53,       /*!< \brief Secondary flow strength objective function definition. */
  SURFACE_MOM_DISTORTION = 54,  /*!< \brief Momentum distortion objective function definition. */
  SURFACE_SECOND_OVER_UNIFORM = 55, /*!< \brief Secondary over uniformity (relative secondary strength) objective function definition. */
  SURFACE_PRESSURE_DROP = 56, 	/*!< \brief Pressure drop objective function definition. */
  CUSTOM_OBJFUNC = 31, 	        /*!< \brief Custom objective function definition. */
  TOTAL_PRESSURE_LOSS = 39,
  KINETIC_ENERGY_LOSS = 40,
  TOTAL_EFFICIENCY = 41,
  TOTAL_STATIC_EFFICIENCY = 42,
  EULERIAN_WORK = 43,
  TOTAL_ENTHALPY_IN = 44,
  FLOW_ANGLE_IN = 45,
  FLOW_ANGLE_OUT = 46,
  MASS_FLOW_IN = 47,
  MASS_FLOW_OUT = 48,
  PRESSURE_RATIO = 49,
  ENTROPY_GENERATION = 50,
  REFERENCE_GEOMETRY = 60,      /*!< \brief Norm of displacements with respect to target geometry. */
  REFERENCE_NODE = 61,          /*!< \brief Objective function defined as the difference of a particular node respect to a reference position. */
  VOLUME_FRACTION = 62,         /*!< \brief Volume average physical density, for material-based topology optimization applications. */
  TOPOL_DISCRETENESS = 63,      /*!< \brief Measure of the discreteness of the current topology. */
  TOPOL_COMPLIANCE = 64         /*!< \brief Measure of the discreteness of the current topology. */
};
static const MapType<string, ENUM_OBJECTIVE> Objective_Map = {
  MakePair("DRAG", DRAG_COEFFICIENT)
  MakePair("LIFT", LIFT_COEFFICIENT)
  MakePair("SIDEFORCE", SIDEFORCE_COEFFICIENT)
  MakePair("EFFICIENCY", EFFICIENCY)
  MakePair("INVERSE_DESIGN_PRESSURE", INVERSE_DESIGN_PRESSURE)
  MakePair("INVERSE_DESIGN_HEATFLUX", INVERSE_DESIGN_HEATFLUX)
  MakePair("MOMENT_X", MOMENT_X_COEFFICIENT)
  MakePair("MOMENT_Y", MOMENT_Y_COEFFICIENT)
  MakePair("MOMENT_Z", MOMENT_Z_COEFFICIENT)
  MakePair("EQUIVALENT_AREA", EQUIVALENT_AREA)
  MakePair("NEARFIELD_PRESSURE", NEARFIELD_PRESSURE)
  MakePair("FORCE_X", FORCE_X_COEFFICIENT)
  MakePair("FORCE_Y", FORCE_Y_COEFFICIENT)
  MakePair("FORCE_Z", FORCE_Z_COEFFICIENT)
  MakePair("THRUST", THRUST_COEFFICIENT)
  MakePair("TORQUE", TORQUE_COEFFICIENT)
  MakePair("TOTAL_HEATFLUX", TOTAL_HEATFLUX)
  MakePair("MAXIMUM_HEATFLUX", MAXIMUM_HEATFLUX)
  MakePair("TOTAL_AVG_TEMPERATURE", TOTAL_AVG_TEMPERATURE)
  MakePair("FIGURE_OF_MERIT", FIGURE_OF_MERIT)
  MakePair("BUFFET", BUFFET_SENSOR)
  MakePair("SURFACE_TOTAL_PRESSURE", SURFACE_TOTAL_PRESSURE)
  MakePair("SURFACE_STATIC_PRESSURE", SURFACE_STATIC_PRESSURE)
  MakePair("SURFACE_MASSFLOW", SURFACE_MASSFLOW)
  MakePair("SURFACE_MACH", SURFACE_MACH)
  MakePair("SURFACE_UNIFORMITY", SURFACE_UNIFORMITY)
  MakePair("SURFACE_SECONDARY", SURFACE_SECONDARY)
  MakePair("SURFACE_MOM_DISTORTION", SURFACE_MOM_DISTORTION)
  MakePair("SURFACE_SECOND_OVER_UNIFORM", SURFACE_SECOND_OVER_UNIFORM)
  MakePair("SURFACE_PRESSURE_DROP", SURFACE_PRESSURE_DROP)
  MakePair("CUSTOM_OBJFUNC", CUSTOM_OBJFUNC)
  MakePair("TOTAL_EFFICIENCY", TOTAL_EFFICIENCY)
  MakePair("TOTAL_STATIC_EFFICIENCY", TOTAL_STATIC_EFFICIENCY)
  MakePair("TOTAL_PRESSURE_LOSS", TOTAL_PRESSURE_LOSS)
  MakePair("EULERIAN_WORK", EULERIAN_WORK)
  MakePair("TOTAL_ENTHALPY_IN", TOTAL_ENTHALPY_IN)
  MakePair("FLOW_ANGLE_IN", FLOW_ANGLE_IN)
  MakePair("FLOW_ANGLE_OUT", FLOW_ANGLE_OUT)
  MakePair("MASS_FLOW_IN", MASS_FLOW_IN)
  MakePair("MASS_FLOW_OUT", MASS_FLOW_OUT)
  MakePair("PRESSURE_RATIO",  PRESSURE_RATIO)
  MakePair("ENTROPY_GENERATION",  ENTROPY_GENERATION)
  MakePair("KINETIC_ENERGY_LOSS", KINETIC_ENERGY_LOSS)
  MakePair("REFERENCE_GEOMETRY", REFERENCE_GEOMETRY)
  MakePair("REFERENCE_NODE", REFERENCE_NODE)
  MakePair("VOLUME_FRACTION", VOLUME_FRACTION)
  MakePair("TOPOL_DISCRETENESS", TOPOL_DISCRETENESS)
  MakePair("TOPOL_COMPLIANCE", TOPOL_COMPLIANCE)
};

/*!
 * \brief Types of residual criteria equations
 */
enum ENUM_RESIDUAL {
    RHO_RESIDUAL = 1, 	     /*!< \brief Rho equation residual criteria equation. */
    RHO_ENERGY_RESIDUAL = 2  /*!< \brief RhoE equation residual criteria equation. */
};
static const MapType<string, ENUM_RESIDUAL> Residual_Map = {
  MakePair("RHO", RHO_RESIDUAL)
  MakePair("RHO_ENERGY", RHO_ENERGY_RESIDUAL)
};

/*!
 * \brief Types of residual criteria for structural problems
 */
enum ENUM_RESFEM {
  RESFEM_RELATIVE = 1,         /*!< \brief Relative criteria: Res/Res0. */
  RESFEM_ABSOLUTE = 2          /*!< \brief Absolute criteria: abs(Res). */
};
static const MapType<string, ENUM_RESFEM> ResFem_Map = {
  MakePair("RELATIVE", RESFEM_RELATIVE)
  MakePair("ABSOLUTE", RESFEM_ABSOLUTE)
};

/*!
 * \brief Types of sensitivities to compute
 */
enum ENUM_SENS {
  SENS_GEOMETRY = 1,    /*!< \brief Geometrical sensitivity. */
  SENS_MACH = 2,        /*!< \brief Mach number sensitivity. */
  SENS_AOA = 3,         /*!< \brief Angle of attack sensitivity. */
  SENS_AOS = 4          /*!< \brief Angle of Sideslip sensitivity. */
};
static const MapType<string, ENUM_SENS> Sens_Map = {
  MakePair("SENS_GEOMETRY", SENS_GEOMETRY)
  MakePair("SENS_MACH", SENS_MACH)
  MakePair("SENS_AOA", SENS_AOA)
  MakePair("SENS_AOS", SENS_AOS)
};

/*!
 * \brief Types of grid adaptation/refinement
 */
enum ENUM_ADAPT {
  NO_ADAPT = 0,           /*!< \brief No grid adaptation. */
  FULL = 1,               /*!< \brief Do a complete grid refinement of all the computational grids. */
  FULL_FLOW = 2,          /*!< \brief Do a complete grid refinement of the flow grid. */
  FULL_ADJOINT = 3,       /*!< \brief Do a complete grid refinement of the adjoint grid. */
  GRAD_FLOW = 5,          /*!< \brief Do a gradient based grid adaptation of the flow grid. */
  GRAD_ADJOINT = 6,       /*!< \brief Do a gradient based grid adaptation of the adjoint grid. */
  GRAD_FLOW_ADJ = 7,      /*!< \brief Do a gradient based grid adaptation of the flow and adjoint grid. */
  COMPUTABLE = 9,         /*!< \brief Apply a computable error grid adaptation. */
  REMAINING = 10,         /*!< \brief Apply a remaining error grid adaptation. */
  WAKE = 12,              /*!< \brief Do a grid refinement on the wake. */
  SMOOTHING = 14,         /*!< \brief Do a grid smoothing of the geometry. */
  SUPERSONIC_SHOCK = 15,  /*!< \brief Do a grid smoothing. */
  PERIODIC = 17           /*!< \brief Add the periodic halo cells. */
};
static const MapType<string, ENUM_ADAPT> Adapt_Map = {
  MakePair("NONE", NO_ADAPT)
  MakePair("FULL", FULL)
  MakePair("FULL_FLOW", FULL_FLOW)
  MakePair("FULL_ADJOINT", FULL_ADJOINT)
  MakePair("GRAD_FLOW", GRAD_FLOW)
  MakePair("GRAD_ADJOINT", GRAD_ADJOINT)
  MakePair("GRAD_FLOW_ADJ", GRAD_FLOW_ADJ)
  MakePair("COMPUTABLE", COMPUTABLE)
  MakePair("REMAINING", REMAINING)
  MakePair("WAKE", WAKE)
  MakePair("SMOOTHING", SMOOTHING)
  MakePair("SUPERSONIC_SHOCK", SUPERSONIC_SHOCK)
  MakePair("PERIODIC", PERIODIC)
};

/*!
 * \brief Types of input file formats
 */
enum ENUM_INPUT {
  SU2       = 1,  /*!< \brief SU2 input format. */
  CGNS_GRID = 2,  /*!< \brief CGNS input format for the computational grid. */
  RECTANGLE = 3,  /*!< \brief 2D rectangular mesh with N x M points of size Lx x Ly. */
  BOX       = 4   /*!< \brief 3D box mesh with N x M x L points of size Lx x Ly x Lz. */
};
static const MapType<string, ENUM_INPUT> Input_Map = {
  MakePair("SU2", SU2)
  MakePair("CGNS", CGNS_GRID)
  MakePair("RECTANGLE", RECTANGLE)
  MakePair("BOX", BOX)
};

/*!
 * \brief Type of solution output file formats
 */
enum ENUM_OUTPUT {
  TECPLOT                 = 1,  /*!< \brief Tecplot format for the solution output. */
  TECPLOT_BINARY          = 2,  /*!< \brief Tecplot binary format for the solution output. */
  SURFACE_TECPLOT         = 3,  /*!< \brief Tecplot format for the solution output. */
  SURFACE_TECPLOT_BINARY  = 4,  /*!< \brief Tecplot binary format for the solution output. */
  CSV                     = 5,  /*!< \brief Comma-separated values format for the solution output. */
  SURFACE_CSV             = 6,  /*!< \brief Comma-separated values format for the solution output. */
  PARAVIEW                = 7,  /*!< \brief Paraview ASCII format for the solution output. */
  PARAVIEW_BINARY         = 8,  /*!< \brief Paraview binary format for the solution output. */
  SURFACE_PARAVIEW        = 9,  /*!< \brief Paraview ASCII format for the solution output. */
  SURFACE_PARAVIEW_BINARY = 10, /*!< \brief Paraview binary format for the solution output. */
  MESH                    = 11, /*!< \brief SU2 mesh format. */
  RESTART_BINARY          = 12, /*!< \brief SU2 binary restart format. */
  RESTART_ASCII           = 13, /*!< \brief SU2 ASCII restart format. */
  CGNS                    = 14, /*!< \brief CGNS format. */
  STL                     = 15, /*!< \brief STL ASCII format for surface solution output. */
  STL_BINARY              = 16, /*!< \brief STL binary format for surface solution output. Not implemented yet. */
  PARAVIEW_XML            = 17, /*!< \brief Paraview XML with binary data format */
  SURFACE_PARAVIEW_XML    = 18, /*!< \brief Surface Paraview XML with binary data format */
  PARAVIEW_MULTIBLOCK     = 19  /*!< \brief Paraview XML Multiblock */
};
static const MapType<string, ENUM_OUTPUT> Output_Map = {
  MakePair("TECPLOT_ASCII", TECPLOT)
  MakePair("TECPLOT", TECPLOT_BINARY)
  MakePair("SURFACE_TECPLOT_ASCII", SURFACE_TECPLOT)
  MakePair("SURFACE_TECPLOT", SURFACE_TECPLOT_BINARY)
  MakePair("CSV", CSV)
  MakePair("SURFACE_CSV", SURFACE_CSV)
  MakePair("PARAVIEW_ASCII", PARAVIEW)
  MakePair("PARAVIEW_LEGACY", PARAVIEW_BINARY)
  MakePair("SURFACE_PARAVIEW_ASCII", SURFACE_PARAVIEW)
  MakePair("SURFACE_PARAVIEW_LEGACY", SURFACE_PARAVIEW_BINARY)
  MakePair("PARAVIEW", PARAVIEW_XML)
  MakePair("SURFACE_PARAVIEW", SURFACE_PARAVIEW_XML)
  MakePair("PARAVIEW_MULTIBLOCK", PARAVIEW_MULTIBLOCK)
  MakePair("RESTART_ASCII", RESTART_ASCII)
  MakePair("RESTART", RESTART_BINARY)
  MakePair("CGNS", CGNS)
  MakePair("STL", STL)
  MakePair("STL_BINARY", STL_BINARY)
};

/*!
 * \brief Type of solution output file formats
 */
enum ENUM_TAB_OUTPUT {
  TAB_CSV = 1,            /*!< \brief Comma-separated values format for the solution output. */
  TAB_TECPLOT = 2         /*!< \brief Tecplot format for the solution output. */
};
static const MapType<string, ENUM_TAB_OUTPUT> TabOutput_Map = {
  MakePair("CSV", TAB_CSV)
  MakePair("TECPLOT", TAB_TECPLOT)
};

/*!
 * \brief Type of volume sensitivity file formats (inout to SU2_DOT)
 */
enum ENUM_SENSITIVITY {
  SU2_NATIVE = 1,       /*!< \brief SU2 native binary format for the volume sensitivity input. */
  UNORDERED_ASCII = 2   /*!< \brief Unordered ASCII list (x,y,z,dJ/dx,dJ/dy/dJ/dz) format for the volume sensitivity input. */
};
static const MapType<string, ENUM_SENSITIVITY> Sensitivity_Map = {
  MakePair("SU2_NATIVE", SU2_NATIVE)
  MakePair("UNORDERED_ASCII", UNORDERED_ASCII)
};

/*!
 * \brief Type of jump definition
 */
enum JUMP_DEFINITION {
  DIFFERENCE = 1,     /*!< \brief Jump given by a difference in values. */
  RATIO = 2           /*!< \brief Jump given by a ratio. */
};
static const MapType<string, JUMP_DEFINITION> Jump_Map = {
  MakePair("DIFFERENCE", DIFFERENCE)
  MakePair("RATIO", RATIO)
};

/*!
 * \brief Type of multigrid cycle
 */
enum MG_CYCLE {
  V_CYCLE = 0,        /*!< \brief V cycle. */
  W_CYCLE = 1,        /*!< \brief W cycle. */
  FULLMG_CYCLE = 2    /*!< \brief FullMG cycle. */
};
static const MapType<string, MG_CYCLE> MG_Cycle_Map = {
  MakePair("V_CYCLE", V_CYCLE)
  MakePair("W_CYCLE", W_CYCLE)
  MakePair("FULLMG_CYCLE", FULLMG_CYCLE)
};

/*!
 * \brief Type of solution output variables
 */
enum ENUM_OUTPUT_VARS {
  DENSITY = 1,      /*!< \brief Density. */
  VEL_X = 2,        /*!< \brief X-component of velocity. */
  VEL_Y = 3,        /*!< \brief Y-component of velocity. */
  VEL_Z = 4,        /*!< \brief Z-component of velocity. */
  PRESSURE = 5,     /*!< \brief Static pressure. */
  MACH = 6,         /*!< \brief Mach number. */
  TEMPERATURE = 7,  /*!< \brief Temperature. */
  LAM_VISC = 8,     /*!< \brief Laminar viscosity. */
  EDDY_VISC = 9     /*!< \brief Eddy viscosity. */
};
static const MapType<string, ENUM_OUTPUT_VARS> Output_Vars_Map = {
  MakePair("DENSITY", DENSITY)
  MakePair("VEL_X", VEL_X)
  MakePair("VEL_Y", VEL_Y)
  MakePair("VEL_Z", VEL_Z)
  MakePair("PRESSURE", PRESSURE)
  MakePair("MACH", MACH)
  MakePair("TEMPERATURE", TEMPERATURE)
  MakePair("LAM_VISC", LAM_VISC)
  MakePair("EDDY_VISC", EDDY_VISC)
};

/*!
 * \brief Types of design parameterizations
 */
enum ENUM_PARAM {
  NO_DEFORMATION = 0,         /*!< \brief No deformation. */
  TRANSLATION = 1,            /*!< \brief Surface movement as design variable. */
  ROTATION = 2,               /*!< \brief Surface rotation as design variable. */
  SCALE = 3,                  /*!< \brief Surface rotation as design variable. */
  FFD_SETTING = 10,           /*!< \brief No surface deformation. */
  FFD_CONTROL_POINT = 11,     /*!< \brief Free form deformation for 3D design (change a control point). */
  FFD_NACELLE = 12,           /*!< \brief Free form deformation for 3D design (change a control point). */
  FFD_GULL = 13,              /*!< \brief Free form deformation for 3D design (change a control point). */
  FFD_CAMBER = 14,            /*!< \brief Free form deformation for 3D design (camber change). */
  FFD_TWIST = 15,             /*!< \brief Free form deformation for 3D design (change the twist angle of a section). */
  FFD_THICKNESS = 16,         /*!< \brief Free form deformation for 3D design (thickness change). */
  FFD_ROTATION = 18,          /*!< \brief Free form deformation for 3D design (rotation around a line). */
  FFD_CONTROL_POINT_2D = 19,  /*!< \brief Free form deformation for 2D design (change a control point). */
  FFD_CAMBER_2D = 20,         /*!< \brief Free form deformation for 3D design (camber change). */
  FFD_THICKNESS_2D = 21,      /*!< \brief Free form deformation for 3D design (thickness change). */
  FFD_TWIST_2D = 22,          /*!< \brief Free form deformation for 3D design (camber change). */
  FFD_CONTROL_SURFACE = 23,   /*!< \brief Free form deformation for 3D design (control surface). */
  FFD_ANGLE_OF_ATTACK = 24,   /*!< \brief Angle of attack for FFD problem. */
  HICKS_HENNE = 30,           /*!< \brief Hicks-Henne bump function for airfoil deformation. */
  PARABOLIC = 31,             /*!< \brief Parabolic airfoil definition as design variables. */
  NACA_4DIGITS = 32,          /*!< \brief The four digits NACA airfoil family as design variables. */
  AIRFOIL = 33,               /*!< \brief Airfoil definition as design variables. */
  CST = 34,                   /*!< \brief CST method with Kulfan parameters for airfoil deformation. */
  SURFACE_BUMP = 35,          /*!< \brief Surfacebump function for flat surfaces deformation. */
  SURFACE_FILE = 36,          /*!< \brief Nodal coordinates for surface set using a file (external parameterization). */
  DV_EFIELD = 40,             /*!< \brief Electric field in deformable membranes. */
  DV_YOUNG = 41,
  DV_POISSON = 42,
  DV_RHO = 43,
  DV_RHO_DL = 44,
  TRANSLATE_GRID = 50,        /*!< \brief Translate the volume grid. */
  ROTATE_GRID = 51,           /*!< \brief Rotate the volume grid */
  SCALE_GRID = 52,            /*!< \brief Scale the volume grid. */
  ANGLE_OF_ATTACK = 101       /*!< \brief Angle of attack for airfoils. */
};
static const MapType<string, ENUM_PARAM> Param_Map = {
  MakePair("FFD_SETTING", FFD_SETTING)
  MakePair("FFD_CONTROL_POINT_2D", FFD_CONTROL_POINT_2D)
  MakePair("FFD_TWIST_2D", FFD_TWIST_2D)
  MakePair("FFD_ANGLE_OF_ATTACK", FFD_ANGLE_OF_ATTACK)
  MakePair("FFD_CAMBER_2D", FFD_CAMBER_2D)
  MakePair("FFD_THICKNESS_2D", FFD_THICKNESS_2D)
  MakePair("HICKS_HENNE", HICKS_HENNE)
  MakePair("SURFACE_BUMP", SURFACE_BUMP)
  MakePair("ANGLE_OF_ATTACK", ANGLE_OF_ATTACK)
  MakePair("NACA_4DIGITS", NACA_4DIGITS)
  MakePair("TRANSLATION", TRANSLATION)
  MakePair("ROTATION", ROTATION)
  MakePair("SCALE", SCALE)
  MakePair("FFD_CONTROL_POINT", FFD_CONTROL_POINT)
  MakePair("FFD_ROTATION", FFD_ROTATION)
  MakePair("FFD_CONTROL_SURFACE", FFD_CONTROL_SURFACE)
  MakePair("FFD_NACELLE", FFD_NACELLE)
  MakePair("FFD_GULL", FFD_GULL)
  MakePair("FFD_TWIST", FFD_TWIST)
  MakePair("FFD_CAMBER", FFD_CAMBER)
  MakePair("FFD_THICKNESS", FFD_THICKNESS)
  MakePair("PARABOLIC", PARABOLIC)
  MakePair("AIRFOIL", AIRFOIL)
  MakePair("SURFACE_FILE", SURFACE_FILE)
  MakePair("NO_DEFORMATION", NO_DEFORMATION)
  MakePair("CST", CST)
  MakePair("ELECTRIC_FIELD", DV_EFIELD)
  MakePair("YOUNG_MODULUS", DV_YOUNG)
  MakePair("POISSON_RATIO", DV_POISSON)
  MakePair("STRUCTURAL_DENSITY", DV_RHO)
  MakePair("DEAD_WEIGHT", DV_RHO_DL)
  MakePair("TRANSLATE_GRID", TRANSLATE_GRID)
  MakePair("ROTATE_GRID", ROTATE_GRID)
  MakePair("SCALE_GRID", SCALE_GRID)
};

/*!
 * \brief Types of FFD Blending function
 */
enum ENUM_FFD_BLENDING{
  BSPLINE_UNIFORM = 0,  /*!< \brief BSpline blending */
  BEZIER = 1,           /*!< \brief Bezier blending */
};
static const MapType<string, ENUM_FFD_BLENDING> Blending_Map = {
  MakePair("BSPLINE_UNIFORM", BSPLINE_UNIFORM)
  MakePair("BEZIER", BEZIER)
};

/*!
 * \brief Types of solvers for solving linear systems
 */
enum ENUM_LINEAR_SOLVER {
  STEEPEST_DESCENT = 1,     /*!< \brief Steepest descent method for point inversion algoritm (Free-Form). */
  NEWTON = 2,               /*!< \brief Newton method for point inversion algorithm (Free-Form). */
  QUASI_NEWTON = 3,         /*!< \brief Quasi Newton method for point inversion algorithm (Free-Form). */
  CONJUGATE_GRADIENT = 4,   /*!< \brief Preconditionated conjugate gradient method for grid deformation. */
  FGMRES = 5,               /*!< \brief Flexible Generalized Minimal Residual method. */
  BCGSTAB = 6,              /*!< \brief BCGSTAB - Biconjugate Gradient Stabilized Method (main solver). */
  RESTARTED_FGMRES = 7,     /*!< \brief Flexible Generalized Minimal Residual method with restart. */
  SMOOTHER = 8,             /*!< \brief Iterative smoother. */
  PASTIX_LDLT = 9,          /*!< \brief PaStiX LDLT (complete) factorization. */
  PASTIX_LU = 10,           /*!< \brief PaStiX LU (complete) factorization. */
};
static const MapType<string, ENUM_LINEAR_SOLVER> Linear_Solver_Map = {
  MakePair("STEEPEST_DESCENT", STEEPEST_DESCENT)
  MakePair("NEWTON", NEWTON)
  MakePair("QUASI_NEWTON", QUASI_NEWTON)
  MakePair("CONJUGATE_GRADIENT", CONJUGATE_GRADIENT)
  MakePair("BCGSTAB", BCGSTAB)
  MakePair("FGMRES", FGMRES)
  MakePair("RESTARTED_FGMRES", RESTARTED_FGMRES)
  MakePair("SMOOTHER", SMOOTHER)
  MakePair("PASTIX_LDLT", PASTIX_LDLT)
  MakePair("PASTIX_LU", PASTIX_LU)
};

/*!
 * \brief Types surface continuity at the intersection with the FFD
 */
enum ENUM_FFD_CONTINUITY {
  DERIVATIVE_NONE = 0,    /*!< \brief No derivative continuity. */
  DERIVATIVE_1ST = 1,     /*!< \brief First derivative continuity. */
  DERIVATIVE_2ND = 2,     /*!< \brief Second derivative continuity. */
  USER_INPUT = 3          /*!< \brief User input. */
};
static const MapType<string, ENUM_FFD_CONTINUITY> Continuity_Map = {
  MakePair("NO_DERIVATIVE", DERIVATIVE_NONE)
  MakePair("1ST_DERIVATIVE", DERIVATIVE_1ST)
  MakePair("2ND_DERIVATIVE", DERIVATIVE_2ND)
  MakePair("USER_INPUT", USER_INPUT)
};

/*!
 * \brief Types of coordinates systems for the FFD
 */
enum ENUM_FFD_COORD_SYSTEM {
  CARTESIAN = 0,    /*!< \brief Cartesian coordinate system. */
  CYLINDRICAL = 1,  /*!< \brief Cylindrical coordinate system. */
  SPHERICAL = 2,    /*!< \brief Spherical coordinate system. */
  POLAR = 3         /*!< \brief Polar coordinate system. */
};
static const MapType<string, ENUM_FFD_COORD_SYSTEM> CoordSystem_Map = {
  MakePair("CARTESIAN", CARTESIAN)
  MakePair("CYLINDRICAL", CYLINDRICAL)
  MakePair("SPHERICAL", SPHERICAL)
  MakePair("POLAR", POLAR)
};

/*!
 * \brief Types of sensitivity smoothing
 */
enum ENUM_SENS_SMOOTHING {
  NO_SMOOTH = 0,  /*!< \brief No smoothing. */
  SOBOLEV = 1,    /*!< \brief Sobolev gradient smoothing. */
  BIGRID = 2      /*!< \brief Bi-grid technique smoothing. */
};
static const MapType<string, ENUM_SENS_SMOOTHING> Sens_Smoothing_Map = {
  MakePair("NONE", NO_SMOOTH)
  MakePair("SOBOLEV", SOBOLEV)
  MakePair("BIGRID", BIGRID)
};

/*!
 * \brief Types of preconditioners for the linear solver
 */
enum ENUM_LINEAR_SOLVER_PREC {
  JACOBI = 1,        /*!< \brief Jacobi preconditioner. */
  LU_SGS = 2,        /*!< \brief LU SGS preconditioner. */
  LINELET = 3,       /*!< \brief Line implicit preconditioner. */
  ILU = 4,           /*!< \brief ILU(k) preconditioner. */
  PASTIX_ILU= 5,     /*!< \brief PaStiX ILU(k) preconditioner. */
  PASTIX_LU_P= 6,    /*!< \brief PaStiX LU as preconditioner. */
  PASTIX_LDLT_P= 7,  /*!< \brief PaStiX LDLT as preconditioner. */
};
static const MapType<string, ENUM_LINEAR_SOLVER_PREC> Linear_Solver_Prec_Map = {
  MakePair("JACOBI", JACOBI)
  MakePair("LU_SGS", LU_SGS)
  MakePair("LINELET", LINELET)
  MakePair("ILU", ILU)
  MakePair("PASTIX_ILU", PASTIX_ILU)
  MakePair("PASTIX_LU", PASTIX_LU_P)
  MakePair("PASTIX_LDLT", PASTIX_LDLT_P)
};

/*!
 * \brief Types of analytic definitions for various geometries
 */
enum ENUM_GEO_ANALYTIC {
  NO_GEO_ANALYTIC = 0,   /*!< \brief No analytic definition of the geometry. */
  NACA0012_AIRFOIL = 1,  /*!< \brief Use the analytical definition of the NACA0012 for doing the grid adaptation. */
  NACA4412_AIRFOIL = 2,  /*!< \brief Use the analytical definition of the NACA4412 for doing the grid adaptation. */
  CYLINDER = 3,          /*!< \brief Use the analytical definition of a cylinder for doing the grid adaptation. */
  BIPARABOLIC = 4        /*!< \brief Use the analytical definition of a biparabolic airfoil for doing the grid adaptation. */
};
static const MapType<string, ENUM_GEO_ANALYTIC> Geo_Analytic_Map = {
  MakePair("NONE", NO_GEO_ANALYTIC)
  MakePair("NACA0012_AIRFOIL", NACA0012_AIRFOIL)
  MakePair("NACA4412_AIRFOIL", NACA4412_AIRFOIL)
  MakePair("CYLINDER", CYLINDER)
  MakePair("BIPARABOLIC", BIPARABOLIC)
};

/*!
 * \brief Types of axis orientation
 */
enum ENUM_GEO_DESCRIPTION {
  TWOD_AIRFOIL = 0, /*!< \brief Airfoil analysis. */
  WING = 1, 	    /*!< \brief Wing analysis. */
  FUSELAGE = 2,     /*!< \brief Fuselage analysis. */
  NACELLE = 3       /*!< \brief Nacelle analysis. */
};
static const MapType<string, ENUM_GEO_DESCRIPTION> Geo_Description_Map = {
  MakePair("AIRFOIL", TWOD_AIRFOIL)
  MakePair("WING", WING)
  MakePair("FUSELAGE", FUSELAGE)
  MakePair("NACELLE", NACELLE)
};

/*!
 * \brief Types of schemes for unsteady computations
 */
enum ENUM_UNSTEADY {
  STEADY = 0,            /*!< \brief A steady computation. */
  TIME_STEPPING = 1,     /*!< \brief Use a time stepping strategy for unsteady computations. */
  DT_STEPPING_1ST = 2,	 /*!< \brief Use a dual time stepping strategy for unsteady computations (1st order). */
  DT_STEPPING_2ND = 3,	 /*!< \brief Use a dual time stepping strategy for unsteady computations (2nd order). */
  ROTATIONAL_FRAME = 4,  /*!< \brief Use a rotational source term. */
  HARMONIC_BALANCE = 5   /*!< \brief Use a harmonic balance source term. */
};
static const MapType<string, ENUM_UNSTEADY> TimeMarching_Map = {
  MakePair("NO", STEADY)
  MakePair("TIME_STEPPING", TIME_STEPPING)
  MakePair("DUAL_TIME_STEPPING-1ST_ORDER", DT_STEPPING_1ST)
  MakePair("DUAL_TIME_STEPPING-2ND_ORDER", DT_STEPPING_2ND)
  MakePair("HARMONIC_BALANCE", HARMONIC_BALANCE)
  MakePair("ROTATIONAL_FRAME", ROTATIONAL_FRAME)
};

/*!
 * \brief Types of criteria to determine when the solution is converged
 */
enum ENUM_CONVERGE_CRIT {
  CAUCHY = 1,       /*!< \brief Cauchy criteria to establish the convergence of the code. */
  RESIDUAL = 2      /*!< \brief Residual criteria to establish the convergence of the code. */
};
static const MapType<string, ENUM_CONVERGE_CRIT> Converge_Crit_Map = {
  MakePair("CAUCHY", CAUCHY)
  MakePair("RESIDUAL", RESIDUAL)
};

/*!
 * \brief Types of element stiffnesses imposed for FEA mesh deformation
 */
enum ENUM_DEFORM_STIFFNESS {
  CONSTANT_STIFFNESS = 0,     /*!< \brief Impose a constant stiffness for each element (steel). */
  INVERSE_VOLUME = 1,         /*!< \brief Impose a stiffness for each element that is inversely proportional to cell volume. */
  SOLID_WALL_DISTANCE = 2     /*!< \brief Impose a stiffness for each element that is proportional to the distance from the solid surface. */
};
static const MapType<string, ENUM_DEFORM_STIFFNESS> Deform_Stiffness_Map = {
  MakePair("CONSTANT_STIFFNESS", CONSTANT_STIFFNESS)
  MakePair("INVERSE_VOLUME", INVERSE_VOLUME)
  MakePair("WALL_DISTANCE", SOLID_WALL_DISTANCE)
};

/*!
 * \brief The direct differentation variables.
 */
enum ENUM_DIRECTDIFF_VAR {
  NO_DERIVATIVE = 0,
  D_MACH = 1,         /*!< \brief Derivative w.r.t. the Mach number */
  D_AOA = 2,          /*!< \brief Derivative w.r.t. the angle of attack */
  D_PRESSURE = 3,     /*!< \brief Derivative w.r.t. the freestream pressure */
  D_TEMPERATURE = 4,  /*!< \brief Derivative w.r.t. the freestream temperature */
  D_DENSITY = 5,      /*!< \brief Derivative w.r.t. the freestream density */
  D_TURB2LAM = 6,     /*!< \brief Derivative w.r.t. the turb2lam */
  D_SIDESLIP = 7,     /*!< \brief Derivative w.r.t. the sideslip angle */
  D_VISCOSITY = 8,    /*!< \brief Derivative w.r.t. the viscosity */
  D_REYNOLDS = 9,     /*!< \brief Derivative w.r.t. the reynolds number */
  D_DESIGN = 10,      /*!< \brief Derivative w.r.t. the design?? */
  D_YOUNG = 11,       /*!< \brief Derivative w.r.t. the Young's modulus */
  D_POISSON = 12,     /*!< \brief Derivative w.r.t. the Poisson's ratio */
  D_RHO = 13,         /*!< \brief Derivative w.r.t. the solid density (inertial) */
  D_RHO_DL = 14,      /*!< \brief Derivative w.r.t. the density for dead loads */
  D_EFIELD = 15       /*!< \brief Derivative w.r.t. the electric field */
};
static const MapType<string, ENUM_DIRECTDIFF_VAR> DirectDiff_Var_Map = {
  MakePair("NONE", NO_DERIVATIVE)
  MakePair("MACH", D_MACH)
  MakePair("AOA", D_AOA)
  MakePair("PRESSURE", D_PRESSURE)
  MakePair("TEMPERATURE", D_TEMPERATURE)
  MakePair("DENSITY", D_DENSITY)
  MakePair("TURB2LAM", D_TURB2LAM)
  MakePair("SIDESLIP", D_SIDESLIP)
  MakePair("VISCOSITY", D_VISCOSITY)
  MakePair("REYNOLDS", D_REYNOLDS)
  MakePair("DESIGN_VARIABLES", D_DESIGN)
  MakePair("YOUNG_MODULUS", D_YOUNG)
  MakePair("POISSON_RATIO", D_POISSON)
  MakePair("STRUCTURAL_DENSITY", D_RHO)
  MakePair("STRUCTURAL_DEAD_LOAD", D_RHO_DL)
  MakePair("ELECTRIC_FIELD", D_EFIELD)
};


enum ENUM_RECORDING {
  SOLUTION_VARIABLES = 1,
  MESH_COORDS = 2,
  MESH_DEFORM = 3,
  SOLUTION_AND_MESH = 4
};

/*!
 * \brief Types of schemes for dynamic structural computations
 */
enum ENUM_DYNAMIC {
  STATIC = 0,     /*!< \brief A static structural computation. */
  DYNAMIC = 1     /*!< \brief Use a time stepping strategy for dynamic computations. */
};
static const MapType<string, ENUM_DYNAMIC> Dynamic_Map = {
  MakePair("NO", STATIC)
  MakePair("YES", DYNAMIC)
};

/*!
 * \brief Types of input file formats
 */
enum ENUM_INPUT_REF {
  SU2_REF = 1,              /*!< \brief SU2 input format (from a restart). */
  CUSTOM_REF = 2            /*!< \brief CGNS input format for the computational grid. */
};
static const MapType<string, ENUM_INPUT_REF> Input_Ref_Map = {
  MakePair("SU2", SU2_REF)
  MakePair("CUSTOM", CUSTOM_REF)
};

/*!
 * \brief Vertex-based quantities exchanged during periodic marker communications.
 */
enum PERIODIC_QUANTITIES {
  PERIODIC_NONE       = 99,  /*!< \brief No periodic communication required. */
  PERIODIC_VOLUME     =  1,  /*!< \brief Volume communication for summing total CV (periodic only). */
  PERIODIC_NEIGHBORS  =  2,  /*!< \brief Communication of the number of neighbors for centered schemes (periodic only). */
  PERIODIC_RESIDUAL   =  3,  /*!< \brief Residual and Jacobian communication (periodic only). */
  PERIODIC_LAPLACIAN  =  4,  /*!< \brief Undivided Laplacian communication for JST (periodic only). */
  PERIODIC_MAX_EIG    =  5,  /*!< \brief Maximum eigenvalue communication (periodic only). */
  PERIODIC_SENSOR     =  6,  /*!< \brief Dissipation sensor communication (periodic only). */
  PERIODIC_SOL_GG     =  7,  /*!< \brief Solution gradient communication for Green-Gauss (periodic only). */
  PERIODIC_PRIM_GG    =  8,  /*!< \brief Primitive gradient communication for Green-Gauss (periodic only). */
  PERIODIC_SOL_LS     =  9,  /*!< \brief Solution gradient communication for weighted Least Squares (periodic only). */
  PERIODIC_PRIM_LS    = 10,  /*!< \brief Primitive gradient communication for weighted Least Squares (periodic only). */
  PERIODIC_LIM_SOL_1  = 11,  /*!< \brief Solution limiter communication phase 1 of 2 (periodic only). */
  PERIODIC_LIM_SOL_2  = 12,  /*!< \brief Solution limiter communication phase 2 of 2 (periodic only). */
  PERIODIC_LIM_PRIM_1 = 13,  /*!< \brief Primitive limiter communication phase 1 of 2 (periodic only). */
  PERIODIC_LIM_PRIM_2 = 14,  /*!< \brief Primitive limiter communication phase 2 of 2 (periodic only). */
  PERIODIC_IMPLICIT   = 15,  /*!< \brief Implicit update communication to ensure consistency across periodic boundaries. */
  PERIODIC_SOL_ULS    = 16,  /*!< \brief Solution gradient communication for unwieghted Least Squares (periodic only). */
  PERIODIC_PRIM_ULS   = 17   /*!< \brief Primitive gradient communication for unweighted Least Squares (periodic only). */
};

/*!
 * \brief Vertex-based quantities exchanged in MPI point-to-point communications.
 */
enum MPI_QUANTITIES {
  SOLUTION             =  0,  /*!< \brief Conservative solution communication. */
  SOLUTION_OLD         =  1,  /*!< \brief Conservative solution old communication. */
  SOLUTION_GRADIENT    =  2,  /*!< \brief Conservative solution gradient communication. */
  SOLUTION_LIMITER     =  3,  /*!< \brief Conservative solution limiter communication. */
  SOLUTION_PRED        =  5,  /*!< \brief Solution predicted communication. */
  SOLUTION_PRED_OLD    =  6,  /*!< \brief Solution predicted old communication. */
  SOLUTION_GEOMETRY    =  7,  /*!< \brief Geometry solution communication. */
  PRIMITIVE_GRADIENT   =  8,  /*!< \brief Primitive gradient communication. */
  PRIMITIVE_LIMITER    =  9,  /*!< \brief Primitive limiter communication. */
  UNDIVIDED_LAPLACIAN  = 10,  /*!< \brief Undivided Laplacian communication. */
  MAX_EIGENVALUE       = 11,  /*!< \brief Maximum eigenvalue communication. */
  SENSOR               = 12,  /*!< \brief Dissipation sensor communication. */
  AUXVAR_GRADIENT      = 13,  /*!< \brief Auxiliary variable gradient communication. */
  COORDINATES          = 14,  /*!< \brief Vertex coordinates communication. */
  COORDINATES_OLD      = 15,  /*!< \brief Old vertex coordinates communication. */
  MAX_LENGTH           = 16,  /*!< \brief Maximum length communication. */
  GRID_VELOCITY        = 17,  /*!< \brief Grid velocity communication. */
  CROSS_TERM           = 18,  /*!< \brief Cross term communication. */
  CROSS_TERM_GEOMETRY  = 19,  /*!< \brief Geometric cross term communication. */
  REF_GEOMETRY         = 20,  /*!< \brief Reference geometry communication. */
  SOLUTION_EDDY        = 21,  /*!< \brief Turbulent solution plus eddy viscosity communication. */
  SOLUTION_MATRIX      = 22,  /*!< \brief Matrix solution communication. */
  SOLUTION_MATRIXTRANS = 23,  /*!< \brief Matrix transposed solution communication. */
  NEIGHBORS            = 24,  /*!< \brief Neighbor point count communication (for JST). */
  SOLUTION_FEA         = 25,  /*!< \brief FEA solution communication. */
  SOLUTION_FEA_OLD     = 26,  /*!< \brief FEA solution old communication. */
  MESH_DISPLACEMENTS   = 27,  /*!< \brief Mesh displacements at the interface. */
  SOLUTION_TIME_N      = 28,  /*!< \brief Solution at time n. */
  SOLUTION_TIME_N1     = 29   /*!< \brief Solution at time n-1. */
};

/*!
 * \brief MPI communication level
 */
enum COMM_LEVEL {
  COMM_NONE    = 0,   /*!< \brief Disable all MPI comms. Purely for testing, as results are incorrect. */
  COMM_MINIMAL = 1,   /*!< \brief Perform only the minimal set of MPI communications for correctness. Disables many console and output comms. */
  COMM_FULL    = 2    /*!< \brief Perform all MPI communications. */
};
static const MapType<string, COMM_LEVEL> Comm_Map = {
  MakePair("NONE",    COMM_NONE)
  MakePair("MINIMAL", COMM_MINIMAL)
  MakePair("FULL",    COMM_FULL)
};

/*
 * \brief Types of filter kernels, initially intended for structural topology optimization applications
 */
enum ENUM_FILTER_KERNEL {
  CONSTANT_WEIGHT_FILTER = 0,      /*!< \brief Uniform weight. */
  CONICAL_WEIGHT_FILTER  = 1,      /*!< \brief Linear decay with distance from center point [Bruns and Tortorelli, 2001]. */
  GAUSSIAN_WEIGHT_FILTER = 2,      /*!< \brief Bell shape around center point [Bruns and Tortorelli, 2003]. */
  DILATE_MORPH_FILTER    = 3,      /*!< \brief Continuous version of the dilate morphology operator [Sigmund 2007]. */
  ERODE_MORPH_FILTER     = 4,      /*!< \brief Continuous version of the erode morphology operator [Sigmund 2007].*/
};
static const MapType<string, ENUM_FILTER_KERNEL> Filter_Kernel_Map = {
  MakePair("CONSTANT", CONSTANT_WEIGHT_FILTER)
  MakePair("CONICAL" , CONICAL_WEIGHT_FILTER)
  MakePair("GAUSSIAN", GAUSSIAN_WEIGHT_FILTER)
  MakePair("DILATE"  , DILATE_MORPH_FILTER)
  MakePair("ERODE"   , ERODE_MORPH_FILTER)
};

/*!
 * \brief Types of projection function, initially intended for structural topology optimization applications
 */
enum ENUM_PROJECTION_FUNCTION {
  NO_PROJECTION  = 0,      /*!< \brief No projection. */
  HEAVISIDE_UP   = 1,      /*!< \brief Project values towards 1. */
  HEAVISIDE_DOWN = 2,      /*!< \brief Project values towards 0. */
};
static const MapType<string, ENUM_PROJECTION_FUNCTION> Projection_Function_Map = {
  MakePair("NO_PROJECTION" , NO_PROJECTION)
  MakePair("HEAVISIDE_UP"  , HEAVISIDE_UP)
  MakePair("HEAVISIDE_DOWN", HEAVISIDE_DOWN)
};

/*!
 * \brief the different validation solution
 */
enum ENUM_VERIFICATION_SOLUTIONS {
  NO_VERIFICATION_SOLUTION =  0,       /*!< \brief No verification solution, standard solver mode. */
  INVISCID_VORTEX          =  1,       /*!< \brief Inviscid vortex. Exact solution of the unsteady Euler equations. */
  RINGLEB                  =  2,       /*!< \brief Ringleb flow. Exact solution of the steady Euler equations. */
  NS_UNIT_QUAD             = 31,       /*!< \brief Exact solution of the laminar Navier Stokes equations without heat conduction. */
  TAYLOR_GREEN_VORTEX      = 32,       /*!< \brief Taylor Green Vortex. */
  INC_TAYLOR_GREEN_VORTEX  = 33,       /*!< \brief Incompressible Taylor Green Vortex (2D). */
  MMS_NS_UNIT_QUAD         = 61,       /*!< \brief Manufactured solution of the laminar Navier Stokes equations on a unit quad. */
  MMS_NS_UNIT_QUAD_WALL_BC = 62,       /*!< \brief Manufactured solution of the laminar Navier Stokes equations on a unit quad with wall BC's. */
  MMS_NS_TWO_HALF_CIRCLES  = 63,       /*!< \brief Manufactured solution of the laminar Navier Stokes equations between two half circles. */
  MMS_NS_TWO_HALF_SPHERES  = 64,       /*!< \brief Manufactured solution of the laminar Navier Stokes equations between two half spheres. */
  MMS_INC_EULER            = 65,       /*!< \brief Manufactured solution of the incompressible Euler equations. */
  MMS_INC_NS               = 66,       /*!< \brief Manufactured solution of the laminar incompressible Navier Stokes equations. */
  USER_DEFINED_SOLUTION    = 99,       /*!< \brief User defined solution. */
};
static const MapType<string, ENUM_VERIFICATION_SOLUTIONS> Verification_Solution_Map = {
  MakePair("NO_VERIFICATION_SOLUTION", NO_VERIFICATION_SOLUTION)
  MakePair("INVISCID_VORTEX",          INVISCID_VORTEX)
  MakePair("RINGLEB",                  RINGLEB)
  MakePair("NS_UNIT_QUAD",             NS_UNIT_QUAD)
  MakePair("TAYLOR_GREEN_VORTEX",      TAYLOR_GREEN_VORTEX)
  MakePair("INC_TAYLOR_GREEN_VORTEX",  INC_TAYLOR_GREEN_VORTEX)
  MakePair("MMS_NS_UNIT_QUAD",         MMS_NS_UNIT_QUAD)
  MakePair("MMS_NS_UNIT_QUAD_WALL_BC", MMS_NS_UNIT_QUAD_WALL_BC)
  MakePair("MMS_NS_TWO_HALF_CIRCLES",  MMS_NS_TWO_HALF_CIRCLES)
  MakePair("MMS_NS_TWO_HALF_SPHERES",  MMS_NS_TWO_HALF_SPHERES)
  MakePair("MMS_INC_EULER",            MMS_INC_EULER)
  MakePair("MMS_INC_NS",               MMS_INC_NS)
  MakePair("USER_DEFINED_SOLUTION",    USER_DEFINED_SOLUTION)
};

#undef MakePair
/* END_CONFIG_ENUMS */

class COptionBase {
private:
  vector<string> value;
public:
  COptionBase() {};
  virtual  ~COptionBase() = 0;

  virtual string SetValue(vector<string> value){this->value = value; return "";}
  vector<string> GetValue() {return value;}
  virtual void SetDefault() = 0;

  string optionCheckMultipleValues(vector<string> & option_value, string type_id, string option_name) {
    if (option_value.size() != 1) {
      string newString;
      newString.append(option_name);
      newString.append(": multiple values for type ");
      newString.append(type_id);
      return newString;
    }
    return "";
  }

  string badValue(vector<string> & option_value, string type_id, string option_name) {
    string newString;
    newString.append(option_name);
    newString.append(": improper option value for type ");
    newString.append(type_id);
    return newString;
  }
};

inline COptionBase::~COptionBase() {}

#ifdef ENABLE_MAPS
#include "option_structure.inl"
#endif
