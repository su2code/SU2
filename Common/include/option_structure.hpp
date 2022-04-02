/*!
 * \file option_structure.hpp
 * \brief Defines classes for referencing options for easy input in CConfig
 * \author J. Hicken, B. Tracey
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "./parallelization/mpi_structure.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <cassert>

/*!
 * \class CEmptyMap
 * \brief We use this dummy class instead of std::map when
 * we only need the enum definition and not the string to
 * enum maps, this makes compilation much faster.
 */
template <typename T, typename U>
struct CEmptyMap {
  CEmptyMap(std::initializer_list<std::pair<const T, U> >) {}
};

#ifdef ENABLE_MAPS
template<class T, class U>
using MapType = std::map<T,U>;
#define MakePair(a,b) {a,b},
#else
template<class T, class U>
using MapType = CEmptyMap<T,U>;
#define MakePair(a,b)
#endif

/*!
 * \brief Different software components of SU2
 */
enum class SU2_COMPONENT {
  SU2_CFD, /*!< \brief Running the SU2_CFD software. */
  SU2_DEF, /*!< \brief Running the SU2_DEF software. */
  SU2_DOT, /*!< \brief Running the SU2_DOT software. */
  SU2_GEO, /*!< \brief Running the SU2_GEO software. */
  SU2_SOL  /*!< \brief Running the SU2_SOL software. */
};

const unsigned int EXIT_DIVERGENCE = 2;   /*!< \brief Exit code (divergence). */

const unsigned int MAX_PARAMETERS = 10;       /*!< \brief Maximum number of parameters for a design variable definition. */
const unsigned int MAX_NUMBER_PERIODIC = 10;  /*!< \brief Maximum number of periodic boundary conditions. */
const unsigned int MAX_STRING_SIZE = 200;     /*!< \brief Maximum number of domains. */
const unsigned int MAX_NUMBER_FFD = 15;       /*!< \brief Maximum number of FFDBoxes for the FFD. */
enum: unsigned int{MAX_SOLS = 13};            /*!< \brief Maximum number of solutions at the same time (dimension of solution container array). */
const unsigned int MAX_TERMS = 7;             /*!< \brief Maximum number of terms in the numerical equations (dimension of solver container array). */
const unsigned int MAX_ZONES = 3;             /*!< \brief Maximum number of zones. */
const unsigned int MAX_FE_KINDS = 7;          /*!< \brief Maximum number of Finite Elements. */
const unsigned int NO_RK_ITER = 0;            /*!< \brief No Runge-Kutta iteration. */

const unsigned int OVERHEAD = 4;    /*!< \brief Overhead space above nMarker when allocating space for boundary elems (MPI + periodic). */

const unsigned int MESH_0 = 0;  /*!< \brief Definition of the finest grid level. */
const unsigned int MESH_1 = 1;  /*!< \brief Definition of the finest grid level. */
const unsigned int ZONE_0 = 0;  /*!< \brief Definition of the first grid domain. */
const unsigned int ZONE_1 = 1;  /*!< \brief Definition of the second grid domain. */
const unsigned int INST_0 = 0;  /*!< \brief Definition of the first instance per grid level. */

const su2double STANDARD_GRAVITY = 9.80665;           /*!< \brief Acceleration due to gravity at surface of earth. */
const su2double UNIVERSAL_GAS_CONSTANT = 8.3144598;   /*!< \brief Universal gas constant in J/(mol*K) */
const su2double BOLTZMANN_CONSTANT = 1.3806503E-23;   /*! \brief Boltzmann's constant [J K^-1] */
const su2double AVOGAD_CONSTANT = 6.0221415E26; /*!< \brief Avogardro's constant, number of particles in one kmole. */

const su2double EPS = 1.0E-16;        /*!< \brief Error scale. */
const su2double TURB_EPS = 1.0E-16;   /*!< \brief Turbulent Error scale. */

const su2double ONE2 = 0.5;         /*!< \brief One divided by two. */
const su2double ONE3 = 1.0 / 3.0;   /*!< \brief One divided by three. */
const su2double TWO3 = 2.0 / 3.0;   /*!< \brief Two divided by three. */
const su2double FOUR3 = 4.0 / 3.0;  /*!< \brief Four divided by three. */

const su2double PI_NUMBER = 4.0 * atan(1.0);  /*!< \brief Pi number. */

const su2double STEFAN_BOLTZMANN = 5.670367E-08;  /*!< \brief Stefan-Boltzmann constant in W/(m^2*K^4). */

const int MASTER_NODE = 0;      /*!< \brief Master node for MPI parallelization. */
const int SINGLE_NODE = 1;      /*!< \brief There is only a node in the MPI parallelization. */
const int SINGLE_ZONE = 1;      /*!< \brief There is only a zone. */

const unsigned short COMM_TYPE_UNSIGNED_LONG  = 1;  /*!< \brief Communication type for unsigned long. */
const unsigned short COMM_TYPE_LONG           = 2;  /*!< \brief Communication type for long. */
const unsigned short COMM_TYPE_UNSIGNED_SHORT = 3;  /*!< \brief Communication type for unsigned short. */
const unsigned short COMM_TYPE_DOUBLE         = 4;  /*!< \brief Communication type for double. */
const unsigned short COMM_TYPE_CHAR           = 5;  /*!< \brief Communication type for char. */
const unsigned short COMM_TYPE_SHORT          = 6;  /*!< \brief Communication type for short. */
const unsigned short COMM_TYPE_INT            = 7;  /*!< \brief Communication type for int. */

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
constexpr unsigned short N_ELEM_TYPES = 7;           /*!< \brief General output & CGNS defines. */

constexpr unsigned short N_POINTS_LINE = 2;          /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_TRIANGLE = 3;      /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_QUADRILATERAL = 4; /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_TETRAHEDRON = 4;   /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_HEXAHEDRON = 8;    /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_PYRAMID = 5;       /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_PRISM = 6;         /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_MAXIMUM = 8;       /*!< \brief Max. out of the above, used for static arrays, keep it up to date. */

constexpr unsigned short N_FACES_LINE = 1;           /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_TRIANGLE = 3;       /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_QUADRILATERAL = 4;  /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_TETRAHEDRON = 4;    /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_PYRAMID = 5;        /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_PRISM = 5;          /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_HEXAHEDRON = 6;     /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_MAXIMUM = 6;        /*!< \brief Max. out of the above, used for static arrays, keep it up to date. */

/*!
 * \brief Get the number of faces of the element.
 * \param[in] elementType - element type
 * \return number of faces
 */
inline unsigned short nFacesOfElementType(unsigned short elementType) {
  switch (elementType) {
    case LINE: return N_FACES_LINE;
    case TRIANGLE: return N_FACES_TRIANGLE;
    case QUADRILATERAL: return N_FACES_QUADRILATERAL;
    case TETRAHEDRON: return N_FACES_TETRAHEDRON;
    case HEXAHEDRON: return N_FACES_HEXAHEDRON;
    case PYRAMID: return N_FACES_PYRAMID;
    case PRISM: return N_FACES_PRISM;
    default: assert(false && "Invalid element type."); return 0;
  }
}

/*!
 * \brief Get the number of points of the element.
 * \param[in] elementType - element type
 * \return number of points
 */
inline unsigned short nPointsOfElementType(unsigned short elementType) {
  switch (elementType) {
    case LINE: return N_POINTS_LINE;
    case TRIANGLE: return N_POINTS_TRIANGLE;
    case QUADRILATERAL: return N_POINTS_QUADRILATERAL;
    case TETRAHEDRON: return N_POINTS_TETRAHEDRON;
    case HEXAHEDRON: return N_POINTS_HEXAHEDRON;
    case PYRAMID: return N_POINTS_PYRAMID;
    case PRISM: return N_POINTS_PRISM;
    default: assert(false && "Invalid element type."); return 0;
  }
}

const int CGNS_STRING_SIZE = 33; /*!< \brief Length of strings used in the CGNS format. */
const int SU2_CONN_SIZE   = 10;  /*!< \brief Size of the connectivity array that is allocated for each element
                                             that we read from a mesh file in the format [[globalID vtkType n0 n1 n2 n3 n4 n5 n6 n7 n8]. */
const int SU2_CONN_SKIP   = 2;   /*!< \brief Offset to skip the globalID and VTK type at the start of the element connectivity list for each CGNS element. */

const su2double COLORING_EFF_THRESH = 0.875;  /*!< \brief Below this value fallback strategies are used instead. */

/*--- All temperature polynomial fits for the fluid models currently
   assume a quartic form (5 coefficients). For example,
   Cp(T) = b0 + b1*T + b2*T^2 + b3*T^3 + b4*T^4. By default, all coeffs
   are set to zero and will be properly non-dim. in the solver. ---*/
constexpr int N_POLY_COEFFS = 5; /*!< \brief Number of coefficients in temperature polynomial fits for fluid models. */

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
static const MapType<std::string, AVERAGE_TYPE> Average_Map = {
  MakePair("AREA", AVERAGE_AREA)
  MakePair("MASSFLUX", AVERAGE_MASSFLUX)
};

/*!
 * \brief different solver types for the CFD component
 */
enum class MAIN_SOLVER {
  NONE,                        /*!< \brief Definition of no solver. */
  EULER,                       /*!< \brief Definition of the Euler's solver. */
  NAVIER_STOKES,               /*!< \brief Definition of the Navier-Stokes' solver. */
  RANS,                        /*!< \brief Definition of the Reynolds-averaged Navier-Stokes' (RANS) solver. */
  INC_EULER,                   /*!< \brief Definition of the incompressible Euler's solver. */
  INC_NAVIER_STOKES,           /*!< \brief Definition of the incompressible Navier-Stokes' solver. */
  INC_RANS,                    /*!< \brief Definition of the incompressible Reynolds-averaged Navier-Stokes' (RANS) solver. */
  HEAT_EQUATION,               /*!< \brief Definition of the finite volume heat solver. */
  FEM_ELASTICITY,              /*!< \brief Definition of a FEM solver. */
  ADJ_EULER,                   /*!< \brief Definition of the continuous adjoint Euler's solver. */
  ADJ_NAVIER_STOKES,           /*!< \brief Definition of the continuous adjoint Navier-Stokes' solver. */
  ADJ_RANS,                    /*!< \brief Definition of the continuous adjoint Reynolds-averaged Navier-Stokes' (RANS) solver. */
  TEMPLATE_SOLVER,             /*!< \brief Definition of template solver. */
  DISC_ADJ_EULER,              /*!< \brief Definition of the discrete adjoint Euler solver. */
  DISC_ADJ_RANS,               /*!< \brief Definition of the discrete adjoint Reynolds-averaged Navier-Stokes' (RANS) solver. */
  DISC_ADJ_NAVIER_STOKES,      /*!< \brief Definition of the discrete adjoint Navier-Stokes' solver. */
  DISC_ADJ_INC_EULER,          /*!< \brief Definition of the discrete adjoint incompressible Euler solver. */
  DISC_ADJ_INC_RANS,           /*!< \brief Definition of the discrete adjoint incompressible Reynolds-averaged Navier-Stokes' (RANS) solver. */
  DISC_ADJ_INC_NAVIER_STOKES,  /*!< \brief Definition of the discrete adjoint incompressible Navier-Stokes'. */
  DISC_ADJ_HEAT,               /*!< \brief Definition of the discrete adjoint heat solver. */
  DISC_ADJ_FEM_EULER,          /*!< \brief Definition of the discrete adjoint FEM Euler solver. */
  DISC_ADJ_FEM_RANS,           /*!< \brief Definition of the discrete adjoint FEM Reynolds-averaged Navier-Stokes' (RANS) solver. */
  DISC_ADJ_FEM_NS,             /*!< \brief Definition of the discrete adjoint FEM Navier-Stokes' solver. */
  DISC_ADJ_FEM,                /*!< \brief Definition of the discrete adjoint FEM solver. */
  FEM_EULER,                   /*!< \brief Definition of the finite element Euler's solver. */
  FEM_NAVIER_STOKES,           /*!< \brief Definition of the finite element Navier-Stokes' solver. */
  FEM_RANS,                    /*!< \brief Definition of the finite element Reynolds-averaged Navier-Stokes' (RANS) solver. */
  FEM_LES,                     /*!< \brief Definition of the finite element Large Eddy Simulation Navier-Stokes' (LES) solver. */
  MULTIPHYSICS,
  NEMO_EULER,                  /*!< \brief Definition of the NEMO Euler solver. */
  NEMO_NAVIER_STOKES,          /*!< \brief Definition of the NEMO NS solver. */
};
static const MapType<std::string, MAIN_SOLVER> Solver_Map = {
  MakePair("NONE", MAIN_SOLVER::NONE)
  MakePair("EULER", MAIN_SOLVER::EULER)
  MakePair("NAVIER_STOKES", MAIN_SOLVER::NAVIER_STOKES)
  MakePair("RANS", MAIN_SOLVER::RANS)
  MakePair("INC_EULER", MAIN_SOLVER::INC_EULER)
  MakePair("INC_NAVIER_STOKES", MAIN_SOLVER::INC_NAVIER_STOKES)
  MakePair("INC_RANS", MAIN_SOLVER::INC_RANS)
  MakePair("FEM_EULER", MAIN_SOLVER::FEM_EULER)
  MakePair("FEM_NAVIER_STOKES", MAIN_SOLVER::FEM_NAVIER_STOKES)
  MakePair("FEM_RANS", MAIN_SOLVER::FEM_RANS)
  MakePair("FEM_LES", MAIN_SOLVER::FEM_LES)
  MakePair("NEMO_EULER",MAIN_SOLVER::NEMO_EULER)
  MakePair("NEMO_NAVIER_STOKES",MAIN_SOLVER::NEMO_NAVIER_STOKES)
  MakePair("ADJ_EULER", MAIN_SOLVER::ADJ_EULER)
  MakePair("ADJ_NAVIER_STOKES", MAIN_SOLVER::ADJ_NAVIER_STOKES)
  MakePair("ADJ_RANS", MAIN_SOLVER::ADJ_RANS )
  MakePair("HEAT_EQUATION", MAIN_SOLVER::HEAT_EQUATION)
  MakePair("ELASTICITY", MAIN_SOLVER::FEM_ELASTICITY)
  MakePair("DISC_ADJ_EULER", MAIN_SOLVER::DISC_ADJ_EULER)
  MakePair("DISC_ADJ_RANS", MAIN_SOLVER::DISC_ADJ_RANS)
  MakePair("DISC_ADJ_NAVIERSTOKES", MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES)
  MakePair("DISC_ADJ_INC_EULER", MAIN_SOLVER::DISC_ADJ_INC_EULER)
  MakePair("DISC_ADJ_INC_RANS", MAIN_SOLVER::DISC_ADJ_INC_RANS)
  MakePair("DISC_ADJ_INC_NAVIERSTOKES", MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES)
  MakePair("DISC_ADJ_HEAT_EQUATION", MAIN_SOLVER::DISC_ADJ_HEAT)
  MakePair("DISC_ADJ_FEM_EULER", MAIN_SOLVER::DISC_ADJ_FEM_EULER)
  MakePair("DISC_ADJ_FEM_RANS", MAIN_SOLVER::DISC_ADJ_FEM_RANS)
  MakePair("DISC_ADJ_FEM_NS", MAIN_SOLVER::DISC_ADJ_FEM_NS)
  MakePair("DISC_ADJ_FEM", MAIN_SOLVER::DISC_ADJ_FEM)
  MakePair("TEMPLATE_SOLVER", MAIN_SOLVER::TEMPLATE_SOLVER)
  MakePair("MULTIPHYSICS", MAIN_SOLVER::MULTIPHYSICS)
};

/*!
 * \brief Different solver types for multizone problems
 */
enum class ENUM_MULTIZONE {
  MZ_BLOCK_GAUSS_SEIDEL, /*!< \brief Definition of a Block-Gauss-Seidel multizone solver. */
  MZ_BLOCK_JACOBI,       /*!< \brief Definition of a Block-Jacobi solver. */
};
static const MapType<std::string, ENUM_MULTIZONE> Multizone_Map = {
  MakePair("BLOCK_GAUSS_SEIDEL", ENUM_MULTIZONE::MZ_BLOCK_GAUSS_SEIDEL)
  MakePair("BLOCK_JACOBI", ENUM_MULTIZONE::MZ_BLOCK_JACOBI)
};

/*!
 * \brief Material geometric conditions
 */
enum class STRUCT_DEFORMATION {
  SMALL,       /*!< \brief Definition of linear elastic material. */
  LARGE,       /*!< \brief Definition of Neo-Hookean material. */
};
static const MapType<std::string, STRUCT_DEFORMATION> Struct_Map = {
  MakePair("SMALL_DEFORMATIONS", STRUCT_DEFORMATION::SMALL)
  MakePair("LARGE_DEFORMATIONS", STRUCT_DEFORMATION::LARGE)
};

/*!
 * \brief Material model
 */
enum class STRUCT_MODEL {
  LINEAR_ELASTIC,   /*!< \brief Definition of linear elastic material. */
  NEO_HOOKEAN,      /*!< \brief Definition of Neo-Hookean material. */
  KNOWLES,          /*!< \brief Definition of Knowles stored-energy potential */
  IDEAL_DE,         /*!< \brief Definition of ideal Dielectric Elastomer */
};
static const MapType<std::string, STRUCT_MODEL> Material_Map = {
  MakePair("LINEAR_ELASTIC", STRUCT_MODEL::LINEAR_ELASTIC)
  MakePair("NEO_HOOKEAN", STRUCT_MODEL::NEO_HOOKEAN)
  MakePair("KNOWLES", STRUCT_MODEL::KNOWLES)
  MakePair("IDEAL_DE", STRUCT_MODEL::IDEAL_DE)
};

/*!
 * \brief Material compressibility
 */
enum class STRUCT_COMPRESS {
  COMPRESSIBLE,     /*!< \brief Definition of compressible material. */
  NEARLY_INCOMP,    /*!< \brief Definition of nearly incompressible material. */
};
static const MapType<std::string, STRUCT_COMPRESS> MatComp_Map = {
  MakePair("COMPRESSIBLE", STRUCT_COMPRESS::COMPRESSIBLE)
  MakePair("NEARLY_INCOMPRESSIBLE", STRUCT_COMPRESS::NEARLY_INCOMP)
};

/*!
 * \brief Types of interpolators
 */
enum class INTERFACE_INTERPOLATOR {
  NEAREST_NEIGHBOR,      /*!< \brief Nearest Neigbhor interpolation */
  ISOPARAMETRIC,         /*!< \brief Isoparametric interpolation, use CONSERVATIVE_INTERPOLATION=YES for conservative interpolation (S.A. Brown 1997).*/
  WEIGHTED_AVERAGE,      /*!< \brief Sliding Mesh Approach E. Rinaldi 2015 */
  RADIAL_BASIS_FUNCTION, /*!< \brief Radial basis function interpolation. */
};
static const MapType<std::string, INTERFACE_INTERPOLATOR> Interpolator_Map = {
  MakePair("NEAREST_NEIGHBOR", INTERFACE_INTERPOLATOR::NEAREST_NEIGHBOR)
  MakePair("ISOPARAMETRIC",    INTERFACE_INTERPOLATOR::ISOPARAMETRIC)
  MakePair("WEIGHTED_AVERAGE", INTERFACE_INTERPOLATOR::WEIGHTED_AVERAGE)
  MakePair("RADIAL_BASIS_FUNCTION", INTERFACE_INTERPOLATOR::RADIAL_BASIS_FUNCTION)
};

/*!
 * \brief Types of radial basis functions
 */
enum class RADIAL_BASIS {
  WENDLAND_C2,        /*!< \brief Wendland C2 radial basis function. */
  INV_MULTI_QUADRIC,  /*!< \brief Inversed multi quartic biharmonic spline. */
  GAUSSIAN,           /*!< \brief Gaussian basis function. */
  THIN_PLATE_SPLINE,  /*!< \brief Thin plate spline. */
  MULTI_QUADRIC,      /*!< \brief Multi quartic biharmonic spline. */
};
static const MapType<std::string, RADIAL_BASIS> RadialBasisFunction_Map = {
  MakePair("WENDLAND_C2", RADIAL_BASIS::WENDLAND_C2)
  MakePair("INV_MULTI_QUADRIC", RADIAL_BASIS::INV_MULTI_QUADRIC)
  MakePair("GAUSSIAN", RADIAL_BASIS::GAUSSIAN)
  MakePair("THIN_PLATE_SPLINE", RADIAL_BASIS::THIN_PLATE_SPLINE)
  MakePair("MULTI_QUADRIC", RADIAL_BASIS::MULTI_QUADRIC)
};

/*!
 * \brief type of radial spanwise interpolation function for the inlet face
 */
enum class INLET_SPANWISE_INTERP {
  NONE,
  LINEAR_1D,
  AKIMA_1D,
  CUBIC_1D,
};
static const MapType<std::string, INLET_SPANWISE_INTERP> Inlet_SpanwiseInterpolation_Map = {
  MakePair("NONE", INLET_SPANWISE_INTERP::NONE)
  MakePair("LINEAR_1D", INLET_SPANWISE_INTERP::LINEAR_1D)
  MakePair("AKIMA_1D", INLET_SPANWISE_INTERP::AKIMA_1D)
  MakePair("CUBIC_1D", INLET_SPANWISE_INTERP::CUBIC_1D)
};

/*!
 * \brief type of radial spanwise interpolation data type for the inlet face
 */
enum class INLET_INTERP_TYPE {
  VR_VTHETA,
  ALPHA_PHI,
};
static const MapType<std::string, INLET_INTERP_TYPE> Inlet_SpanwiseInterpolationType_Map = {
  MakePair("VR_VTHETA", INLET_INTERP_TYPE::VR_VTHETA)
  MakePair("ALPHA_PHI", INLET_INTERP_TYPE::ALPHA_PHI)
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
enum class ENUM_REGIME {
  COMPRESSIBLE = 0,   /*!< \brief Definition of compressible solver. */
  INCOMPRESSIBLE = 1, /*!< \brief Definition of incompressible solver. */
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
static const MapType<std::string, ENUM_KIND_NONDIM> NonDim_Map = {
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
  SI = 0,     /*!< \brief Definition of compressible solver. */
  US = 1      /*!< \brief Definition of incompressible solver. */
};
static const MapType<std::string, ENUM_MEASUREMENTS> Measurements_Map = {
  MakePair("SI", SI)
  MakePair("US", US)
};

/*!
 * \brief different types of systems
 */
enum RUNTIME_TYPE {
  RUNTIME_FLOW_SYS = 2,       /*!< \brief One-physics case, the code is solving the flow equations(Euler and Navier-Stokes). */
  RUNTIME_TURB_SYS = 3,       /*!< \brief One-physics case, the code is solving the turbulence model. */
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
  RUNTIME_SPECIES_SYS = 25,   /*!< \brief One-physics case, the code is solving the species model. */
  RUNTIME_ADJSPECIES_SYS = 26,/*!< \brief One-physics case, the code is solving the adjoint species model. */
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

const int SPECIES_SOL = 11;    /*!< \brief Position of the species solver. */
const int ADJSPECIES_SOL = 12; /*!< \brief Position of the adjoint of the species solver. */

const int FEA_SOL = 0;      /*!< \brief Position of the FEA equation in the solution solver array. */
const int ADJFEA_SOL = 1;   /*!< \brief Position of the FEA adjoint equation in the solution solver array. */

const int TEMPLATE_SOL = 0; /*!< \brief Position of the template solution. */

const int CONV_TERM = 0;           /*!< \brief Position of the convective terms in the numerics container array. */
const int VISC_TERM = 1;           /*!< \brief Position of the viscous terms in the numerics container array. */
const int SOURCE_FIRST_TERM = 2;   /*!< \brief Position of the first source term in the numerics container array. */
const int SOURCE_SECOND_TERM = 3;  /*!< \brief Position of the second source term in the numerics container array. */
const int CONV_BOUND_TERM = 4;     /*!< \brief Position of the convective boundary terms in the numerics container array. */
const int VISC_BOUND_TERM = 5;     /*!< \brief Position of the viscous boundary terms in the numerics container array. */
const int GRAD_TERM = 6;           /*!< \brief Position of the gradient smoothing terms in the numerics container array. */

const int FEA_TERM = 0;      /*!< \brief Position of the finite element analysis terms in the numerics container array. */
const int DE_TERM = 1;       /*!< \brief Position of the dielectric terms in the numerics container array. */

const int MAT_NHCOMP  = 2;   /*!< \brief Position of the Neo-Hookean compressible material model. */
const int MAT_IDEALDE = 3;   /*!< \brief Position of the Ideal-DE material model. */
const int MAT_KNOWLES = 4;   /*!< \brief Position of the Knowles material model. */

/*!
 * \brief Types of finite elements (in 1D or 2D or 3D)
 */
const int EL_LINE = 6;    /*!< \brief Elements of two nodes, with second order gauss quadrature (1D). */

const int EL_TRIA = 0;    /*!< \brief Elements of three nodes (2D). */
const int EL_QUAD = 1;    /*!< \brief Elements of four nodes (2D). */
const int EL_TRIA2 = 2;   /*!< \brief Elements of three nodes (2D), with second order gauss quadrature. */

const int EL_TETRA = 0;   /*!< \brief Elements of four nodes (3D). */
const int EL_HEXA  = 1;   /*!< \brief Elements of eight nodes (3D). */
const int EL_PYRAM = 2;   /*!< \brief Elements of five nodes (3D). */
const int EL_PRISM = 3;   /*!< \brief Elements of six nodes (3D). */
const int EL_TETRA2 = 4;  /*!< \brief Elements of four nodes, with second order gauss quadrature (3D). */
const int EL_PYRAM2 = 5;  /*!< \brief Elements of five nodes, with third order gauss quadrature (3D). */

/*!
 * \brief Types of spatial discretizations
 */
enum ENUM_SPACE {
  NO_CONVECTIVE = 0,   /*!< \brief No convective scheme is used. */
  SPACE_CENTERED = 1,  /*!< \brief Space centered convective numerical method. */
  SPACE_UPWIND = 2,    /*!< \brief Upwind convective numerical method. */
  FINITE_ELEMENT = 3   /*!< \brief Finite element convective numerical method. */
};
static const MapType<std::string, ENUM_SPACE> Space_Map = {
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
  INC_IDEAL_GAS_POLY = 6, /*!< \brief Inc. ideal gas, polynomial gas model. */
  MUTATIONPP = 7,         /*!< \brief Mutation++ gas model for nonequilibrium flow. */
  SU2_NONEQ = 8           /*!< \brief User defined gas model for nonequilibrium flow. */
};
static const MapType<std::string, ENUM_FLUIDMODEL> FluidModel_Map = {
  MakePair("STANDARD_AIR", STANDARD_AIR)
  MakePair("IDEAL_GAS", IDEAL_GAS)
  MakePair("VW_GAS", VW_GAS)
  MakePair("PR_GAS", PR_GAS)
  MakePair("CONSTANT_DENSITY", CONSTANT_DENSITY)
  MakePair("INC_IDEAL_GAS", INC_IDEAL_GAS)
  MakePair("INC_IDEAL_GAS_POLY", INC_IDEAL_GAS_POLY)
  MakePair("MUTATIONPP", MUTATIONPP)
  MakePair("SU2_NONEQ", SU2_NONEQ)
};

/*!
 * \brief types of gas models
 */
enum ENUM_GASMODEL {
   NO_MODEL   = 0,
   ARGON      = 1,
   AIR7       = 2,
   AIR21      = 3,
   O2         = 4,
   N2         = 5,
   AIR5       = 6,
   ARGON_SID  = 7,
   ONESPECIES = 8
};
static const MapType<std::string, ENUM_GASMODEL> GasModel_Map = {
MakePair("NONE", NO_MODEL)
MakePair("ARGON", ARGON)
MakePair("AIR-7", AIR7)
MakePair("AIR-21", AIR21)
MakePair("O2", O2)
MakePair("N2", N2)
MakePair("AIR-5", AIR5)
MakePair("ARGON-SID",ARGON_SID)
MakePair("ONESPECIES", ONESPECIES)
};

/*!
 * \brief types of coefficient transport model
 */
enum class TRANSCOEFFMODEL {
  WILKE,
  GUPTAYOS,
  CHAPMANN_ENSKOG
};
static const MapType<std::string, TRANSCOEFFMODEL> TransCoeffModel_Map = {
MakePair("WILKE", TRANSCOEFFMODEL::WILKE)
MakePair("GUPTA-YOS", TRANSCOEFFMODEL::GUPTAYOS)
MakePair("CHAPMANN-ENSKOG", TRANSCOEFFMODEL::CHAPMANN_ENSKOG)
};

/*!
 * \brief Types of density models
 */
enum class INC_DENSITYMODEL {
  CONSTANT,   /*!< \brief Constant density. */
  BOUSSINESQ, /*!< \brief Boussinesq density model. */
  VARIABLE,   /*!< \brief Variable density model. */
};
static const MapType<std::string, INC_DENSITYMODEL> DensityModel_Map = {
  MakePair("CONSTANT", INC_DENSITYMODEL::CONSTANT)
  MakePair("BOUSSINESQ", INC_DENSITYMODEL::BOUSSINESQ)
  MakePair("VARIABLE", INC_DENSITYMODEL::VARIABLE)
};

/*!
 * \brief Types of initialization option
 */
enum ENUM_INIT_OPTION {
  REYNOLDS = 0,      /*!< \brief Reynold's number initalization. */
  TD_CONDITIONS = 1  /*!< \brief Total conditions initalization. */
};
static const MapType<std::string, ENUM_INIT_OPTION> InitOption_Map = {
  MakePair("REYNOLDS", REYNOLDS)
  MakePair("TD_CONDITIONS", TD_CONDITIONS)
};

/*!
 * \brief Types of initialization option
 */
enum class FREESTREAM_OPTION {
  TEMPERATURE_FS, /*!< \brief Temperature initialization. */
  DENSITY_FS, /*!< \brief Density initalization. */
};
static const MapType<std::string, FREESTREAM_OPTION> FreeStreamOption_Map = {
  MakePair("TEMPERATURE_FS", FREESTREAM_OPTION::TEMPERATURE_FS)
  MakePair("DENSITY_FS", FREESTREAM_OPTION::DENSITY_FS)
};

/*!
 * \brief Types of viscosity model
 */
enum class VISCOSITYMODEL {
  CONSTANT, /*!< \brief Constant viscosity. */
  SUTHERLAND, /*!< \brief Sutherlands Law viscosity. */
  POLYNOMIAL, /*!< \brief Polynomial viscosity. */
};
static const MapType<std::string, VISCOSITYMODEL> ViscosityModel_Map = {
  MakePair("CONSTANT_VISCOSITY", VISCOSITYMODEL::CONSTANT)
  MakePair("SUTHERLAND", VISCOSITYMODEL::SUTHERLAND)
  MakePair("POLYNOMIAL_VISCOSITY", VISCOSITYMODEL::POLYNOMIAL)
};

/*!
 * \brief Types of thermal conductivity model
 */
enum class CONDUCTIVITYMODEL {
  CONSTANT, /*!< \brief Constant thermal conductivity. */
  CONSTANT_PRANDTL, /*!< \brief Constant Prandtl number. */
  POLYNOMIAL, /*!< \brief Polynomial thermal conductivity. */
};
static const MapType<std::string, CONDUCTIVITYMODEL> ConductivityModel_Map = {
  MakePair("CONSTANT_CONDUCTIVITY", CONDUCTIVITYMODEL::CONSTANT)
  MakePair("CONSTANT_PRANDTL", CONDUCTIVITYMODEL::CONSTANT_PRANDTL)
  MakePair("POLYNOMIAL_CONDUCTIVITY", CONDUCTIVITYMODEL::POLYNOMIAL)
};

/*!
 * \brief Types of turbulent thermal conductivity model
 */
enum class CONDUCTIVITYMODEL_TURB {
  NONE, /*!< \brief No turbulent contribution to the effective thermal conductivity for RANS. */
  CONSTANT_PRANDTL, /*!< \brief Include contribution to effective conductivity using constant turbulent Prandtl number for RANS. */
};
static const MapType<std::string, CONDUCTIVITYMODEL_TURB> TurbConductivityModel_Map = {
  MakePair("NONE", CONDUCTIVITYMODEL_TURB::NONE)
  MakePair("CONSTANT_PRANDTL_TURB", CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL)
};

/*!
 * \brief types of mass diffusivity models
 */
enum class DIFFUSIVITYMODEL {
  CONSTANT_DIFFUSIVITY, /*!< \brief Constant mass diffusivity for scalar transport. */
  CONSTANT_SCHMIDT,     /*!< \brief Constant Schmidt number for mass diffusion in scalar transport. */
  UNITY_LEWIS,          /*!< \brief Unity Lewis model */
};

static const MapType<std::string, DIFFUSIVITYMODEL> Diffusivity_Model_Map = {
  MakePair("CONSTANT_DIFFUSIVITY", DIFFUSIVITYMODEL::CONSTANT_DIFFUSIVITY)
  MakePair("CONSTANT_SCHMIDT", DIFFUSIVITYMODEL::CONSTANT_SCHMIDT)
  MakePair("UNITY_LEWIS", DIFFUSIVITYMODEL::UNITY_LEWIS)
};

/*!
 * \brief Types of unsteady mesh motion
 */
enum ENUM_GRIDMOVEMENT {
  NO_MOVEMENT = 0,          /*!< \brief Simulation on a static mesh. */
  RIGID_MOTION = 2,         /*!< \brief Simulation with rigid mesh motion (plunging/pitching/rotation). */
  ROTATING_FRAME = 8,       /*!< \brief Simulation in a rotating frame. */
  STEADY_TRANSLATION = 11,  /*!< \brief Simulation in a steadily translating frame. */
  GUST = 12,                /*!< \brief Simulation on a static mesh with a gust. */
};
static const MapType<std::string, ENUM_GRIDMOVEMENT> GridMovement_Map = {
  MakePair("NONE", NO_MOVEMENT)
  MakePair("RIGID_MOTION", RIGID_MOTION)
  MakePair("ROTATING_FRAME", ROTATING_FRAME)
  MakePair("STEADY_TRANSLATION", STEADY_TRANSLATION)
  MakePair("GUST", GUST)
};

enum ENUM_SURFACEMOVEMENT {
  DEFORMING = 1,                 /*!< \brief Simulation with deformation. */
  MOVING_WALL = 2,               /*!< \brief Simulation with moving wall. */
  AEROELASTIC = 3,               /*!< \brief Simulation with aeroelastic motion. */
  AEROELASTIC_RIGID_MOTION = 4,  /*!< \brief Simulation with rotation and aeroelastic motion. */
  EXTERNAL = 6,                  /*!< \brief Simulation with external motion. */
  EXTERNAL_ROTATION = 7,         /*!< \brief Simulation with external rotation motion. */
};
static const MapType<std::string, ENUM_SURFACEMOVEMENT> SurfaceMovement_Map = {
  MakePair("DEFORMING", DEFORMING)
  MakePair("MOVING_WALL", MOVING_WALL)
  MakePair("AEROELASTIC_RIGID_MOTION", AEROELASTIC_RIGID_MOTION)
  MakePair("AEROELASTIC", AEROELASTIC)
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
static const MapType<std::string, ENUM_GUST_TYPE> Gust_Type_Map = {
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
static const MapType<std::string, ENUM_GUST_DIR> Gust_Dir_Map = {
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
  JST_MAT = 3,        /*!< \brief JST with matrix dissipation. */
  JST_KE = 4          /*!< \brief Kinetic Energy preserving Jameson-Smith-Turkel centered numerical method. */
};
static const MapType<std::string, ENUM_CENTERED> Centered_Map = {
  MakePair("NONE", NO_CENTERED)
  MakePair("JST", JST)
  MakePair("JST_KE", JST_KE)
  MakePair("JST_MAT", JST_MAT)
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
  AUSMPLUSUP2 = 17,            /*!< \brief AUSM+ -up2 numerical method (All Speed) */
  AUSMPWPLUS = 18            /*!< \brief AUSMplus numerical method. (MAYBE for TNE2 ONLY)*/
};
static const MapType<std::string, ENUM_UPWIND> Upwind_Map = {
  MakePair("NONE", NO_UPWIND)
  MakePair("ROE", ROE)
  MakePair("TURKEL_PREC", TURKEL)
  MakePair("AUSM", AUSM)
  MakePair("AUSMPLUSUP", AUSMPLUSUP)
  MakePair("AUSMPLUSUP2", AUSMPLUSUP2)
  MakePair("AUSMPWPLUS", AUSMPWPLUS)
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
static const MapType<std::string, ENUM_FEM> FEM_Map = {
  MakePair("NONE", NO_FEM)
  MakePair("DG", DG)
};

/*!
 * \brief Types of shock capturing method in Discontinuous Galerkin numerical method.
 */
enum class FEM_SHOCK_CAPTURING_DG {
  NONE,     /*!< \brief Shock capturing is not used. */
  PERSSON   /*!< \brief Per-Olof Persson's sub-cell shock capturing method. */
};
static const MapType<std::string, FEM_SHOCK_CAPTURING_DG> ShockCapturingDG_Map = {
  MakePair("NONE", FEM_SHOCK_CAPTURING_DG::NONE)
  MakePair("PERSSON", FEM_SHOCK_CAPTURING_DG::PERSSON)
};

/*!
 * \brief Types of matrix coloring to compute a sparse Jacobian matrix.
 */
enum ENUM_MATRIX_COLORING {
  GREEDY_COLORING = 0,            /*!< \brief Greedy type of algorithm for the coloring. */
  NATURAL_COLORING = 1            /*!< \brief One color for every DOF, very slow. Only to be used for debugging. */
};
static const MapType<std::string, ENUM_MATRIX_COLORING> MatrixColoring_Map = {
  MakePair("GREEDY_COLORING", GREEDY_COLORING)
  MakePair("NATURAL_COLORING", NATURAL_COLORING)
};

/*!
 * \brief Types of slope limiters
 */
enum class LIMITER {
  NONE                 , /*!< \brief No limiter. */
  VENKATAKRISHNAN      , /*!< \brief Slope limiter using Venkatakrisnan method (stencil formulation). */
  VENKATAKRISHNAN_WANG , /*!< \brief Slope limiter using Venkatakrisnan method, eps based on solution (stencil formulation). */
  BARTH_JESPERSEN      , /*!< \brief Slope limiter using Barth-Jespersen method (stencil formulation). */
  VAN_ALBADA_EDGE      , /*!< \brief Slope limiter using Van Albada method (edge formulation). */
  SHARP_EDGES          , /*!< \brief Slope limiter using sharp edges. */
  WALL_DISTANCE          /*!< \brief Slope limiter using wall distance. */
};
static const MapType<std::string, LIMITER> Limiter_Map = {
  MakePair("NONE", LIMITER::NONE)
  MakePair("VENKATAKRISHNAN", LIMITER::VENKATAKRISHNAN)
  MakePair("VENKATAKRISHNAN_WANG", LIMITER::VENKATAKRISHNAN_WANG)
  MakePair("BARTH_JESPERSEN", LIMITER::BARTH_JESPERSEN)
  MakePair("VAN_ALBADA_EDGE", LIMITER::VAN_ALBADA_EDGE)
  MakePair("SHARP_EDGES", LIMITER::SHARP_EDGES)
  MakePair("WALL_DISTANCE", LIMITER::WALL_DISTANCE)
};

/*!
 * \brief Types of turbulent models
 */
enum class TURB_MODEL {
  NONE,      /*!< \brief No turbulence model. */
  SA,        /*!< \brief Kind of Turbulent model (Spalart-Allmaras). */
  SA_NEG,    /*!< \brief Kind of Turbulent model (Spalart-Allmaras). */
  SA_E,      /*!< \brief Kind of Turbulent model (Spalart-Allmaras Edwards). */
  SA_COMP,   /*!< \brief Kind of Turbulent model (Spalart-Allmaras Compressibility Correction). */
  SA_E_COMP, /*!< \brief Kind of Turbulent model (Spalart-Allmaras Edwards with Compressibility Correction). */
  SST,       /*!< \brief Kind of Turbulence model (Menter SST). */
  SST_SUST   /*!< \brief Kind of Turbulence model (Menter SST with sustaining terms for free-stream preservation). */
};
static const MapType<std::string, TURB_MODEL> Turb_Model_Map = {
  MakePair("NONE", TURB_MODEL::NONE)
  MakePair("SA", TURB_MODEL::SA)
  MakePair("SA_NEG", TURB_MODEL::SA_NEG)
  MakePair("SA_E", TURB_MODEL::SA_E)
  MakePair("SA_COMP", TURB_MODEL::SA_COMP)
  MakePair("SA_E_COMP", TURB_MODEL::SA_E_COMP)
  MakePair("SST", TURB_MODEL::SST)
  MakePair("SST_SUST", TURB_MODEL::SST_SUST)
};

/*!
 * \brief Families of turbulence models
 */
enum class TURB_FAMILY {
  NONE,   /*!< \brief No turbulence model. */
  SA,     /*!< \brief Spalart-Allmaras variants. */
  KW,     /*!< \brief k-w models. */
};
/*!
 * \brief Associate turb models with their family
 */
inline TURB_FAMILY TurbModelFamily(TURB_MODEL model) {
  switch (model) {
    case TURB_MODEL::NONE:
      return TURB_FAMILY::NONE;
    case TURB_MODEL::SA:
    case TURB_MODEL::SA_NEG:
    case TURB_MODEL::SA_E:
    case TURB_MODEL::SA_COMP:
    case TURB_MODEL::SA_E_COMP:
      return TURB_FAMILY::SA;
    case TURB_MODEL::SST:
    case TURB_MODEL::SST_SUST:
      return TURB_FAMILY::KW;
  }
  return TURB_FAMILY::NONE;
}

/*!
 * \brief Types of transition models
 */
enum class TURB_TRANS_MODEL {
  NONE,  /*!< \brief No transition model. */
  LM,    /*!< \brief Kind of transition model (Langtry-Menter (LM) for SST and Spalart-Allmaras). */
  BC    /*!< \brief Kind of transition model (BAS-CAKMAKCIOGLU (BC) for Spalart-Allmaras). */
};
static const MapType<std::string, TURB_TRANS_MODEL> Trans_Model_Map = {
  MakePair("NONE", TURB_TRANS_MODEL::NONE)
  MakePair("LM", TURB_TRANS_MODEL::LM)
  MakePair("BC", TURB_TRANS_MODEL::BC)
};

/*!
 * \brief types of species transport models
 */
enum class SPECIES_MODEL {
  NONE,              /*!< \brief No scalar transport model. */
  PASSIVE_SCALAR,    /*!< \brief Passive scalar transport model. */
};
static const MapType<std::string, SPECIES_MODEL> Species_Model_Map = {
  MakePair("NONE", SPECIES_MODEL::NONE)
  MakePair("PASSIVE_SCALAR", SPECIES_MODEL::PASSIVE_SCALAR)
};

/*!
 * \brief Types of subgrid scale models
 */
enum class TURB_SGS_MODEL {
  NONE        , /*!< \brief No subgrid scale model. */
  IMPLICIT_LES, /*!< \brief Implicit LES, i.e. no explicit SGS model. */
  SMAGORINSKY , /*!< \brief Smagorinsky SGS model. */
  WALE        , /*!< \brief Wall-Adapting Local Eddy-viscosity SGS model. */
  VREMAN        /*!< \brief Vreman SGS model. */
};
static const MapType<std::string, TURB_SGS_MODEL> SGS_Model_Map = {
  MakePair("NONE",         TURB_SGS_MODEL::NONE)
  MakePair("IMPLICIT_LES", TURB_SGS_MODEL::IMPLICIT_LES)
  MakePair("SMAGORINSKY",  TURB_SGS_MODEL::SMAGORINSKY)
  MakePair("WALE",         TURB_SGS_MODEL::WALE)
  MakePair("VREMAN",       TURB_SGS_MODEL::VREMAN)
};


/*!
 * \brief Types of window (weight) functions for cost functional
 */
enum class WINDOW_FUNCTION {
  SQUARE,        /*!< \brief No weight function  (order 1)*/
  HANN,          /*!< \brief Hann-type weight function (order 3) */
  HANN_SQUARE,   /*!< \brief Hann-squared type weight function (order 5)*/
  BUMP,          /*!< \brief bump type weight function (exponential order of convergence) */
};
static const MapType<std::string, WINDOW_FUNCTION> Window_Map = {
  MakePair("SQUARE", WINDOW_FUNCTION::SQUARE)
  MakePair("HANN", WINDOW_FUNCTION::HANN)
  MakePair("HANN_SQUARE", WINDOW_FUNCTION::HANN_SQUARE)
  MakePair("BUMP", WINDOW_FUNCTION::BUMP)
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
static const MapType<std::string, ENUM_HYBRIDRANSLES> HybridRANSLES_Map = {
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
static const MapType<std::string, ENUM_ROELOWDISS> RoeLowDiss_Map = {
  MakePair("NONE", NO_ROELOWDISS)
  MakePair("FD", FD)
  MakePair("NTS", NTS)
  MakePair("NTS_DUCROS", NTS_DUCROS)
  MakePair("FD_DUCROS", FD_DUCROS)
};

/*!
 * \brief Types of wall functions.
 */
enum class WALL_FUNCTIONS {
  NONE                 ,   /*!< \brief No wall function treatment, integration to the wall. Default behavior. */
  STANDARD_FUNCTION    ,   /*!< \brief Standard wall function. */
  ADAPTIVE_FUNCTION    ,   /*!< \brief Adaptive wall function. Formulation depends on y+. */
  SCALABLE_FUNCTION    ,   /*!< \brief Scalable wall function. */
  EQUILIBRIUM_MODEL    ,   /*!< \brief Equilibrium wall model for LES. */
  NONEQUILIBRIUM_MODEL ,   /*!< \brief Non-equilibrium wall model for LES. */
  LOGARITHMIC_MODEL        /*!< \brief Logarithmic law-of-the-wall model for LES. */
};
static const MapType<std::string, WALL_FUNCTIONS> Wall_Functions_Map = {
  MakePair("NO_WALL_FUNCTION",          WALL_FUNCTIONS::NONE)
  MakePair("STANDARD_WALL_FUNCTION",    WALL_FUNCTIONS::STANDARD_FUNCTION)
  MakePair("ADAPTIVE_WALL_FUNCTION",    WALL_FUNCTIONS::ADAPTIVE_FUNCTION)
  MakePair("SCALABLE_WALL_FUNCTION",    WALL_FUNCTIONS::SCALABLE_FUNCTION)
  MakePair("EQUILIBRIUM_WALL_MODEL",    WALL_FUNCTIONS::EQUILIBRIUM_MODEL)
  MakePair("NONEQUILIBRIUM_WALL_MODEL", WALL_FUNCTIONS::NONEQUILIBRIUM_MODEL)
  MakePair("LOGARITHMIC_WALL_MODEL",    WALL_FUNCTIONS::LOGARITHMIC_MODEL)
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
static const MapType<std::string, ENUM_TIME_INT> Time_Int_Map = {
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
static const MapType<std::string, ENUM_ADER_PREDICTOR> Ader_Predictor_Map = {
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
static const MapType<std::string, ENUM_HEAT_TIMESTEP> Heat_TimeStep_Map = {
  MakePair("LOCAL", MINIMUM)
  MakePair("CONVECTIVE", CONVECTIVE)
  MakePair("VISCOUS", VISCOUS)
  MakePair("BYFLOW", BYFLOW)
};

/*!
 * \brief Type of time integration schemes
 */
enum class STRUCT_TIME_INT {
  CD_EXPLICIT,       /*!< \brief Support for implementing an explicit method. */
  NEWMARK_IMPLICIT,  /*!< \brief Implicit Newmark integration definition. */
  GENERALIZED_ALPHA, /*!< \brief Support for implementing another implicit method. */
};
static const MapType<std::string, STRUCT_TIME_INT> Time_Int_Map_FEA = {
  MakePair("CD_EXPLICIT", STRUCT_TIME_INT::CD_EXPLICIT)
  MakePair("NEWMARK_IMPLICIT", STRUCT_TIME_INT::NEWMARK_IMPLICIT)
  MakePair("GENERALIZED_ALPHA", STRUCT_TIME_INT::GENERALIZED_ALPHA)
};

/*!
 * \brief Type of time integration schemes
 */
enum class STRUCT_SPACE_ITE {
  NEWTON,       /*!< \brief Full Newton-Rapshon method. */
  MOD_NEWTON,   /*!< \brief Modified Newton-Raphson method. */
};
static const MapType<std::string, STRUCT_SPACE_ITE> Space_Ite_Map_FEA = {
  MakePair("NEWTON_RAPHSON", STRUCT_SPACE_ITE::NEWTON)
  MakePair("MODIFIED_NEWTON_RAPHSON", STRUCT_SPACE_ITE::MOD_NEWTON)
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
static const MapType<std::string, ENUM_FLOW_GRADIENT> Gradient_Map = {
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
static const MapType<std::string, GEOMETRY_MODE> GeometryMode_Map = {
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
  CUSTOM_BOUNDARY = 10,       /*!< \brief custom boundary definition. */
  DISPLACEMENT_BOUNDARY = 14, /*!< \brief Boundary displacement definition. */
  LOAD_BOUNDARY = 15,         /*!< \brief Boundary Load definition. */
  FLOWLOAD_BOUNDARY = 16,     /*!< \brief Boundary Load definition. */
  SUPERSONIC_INLET = 19,      /*!< \brief Boundary supersonic inlet definition. */
  SUPERSONIC_OUTLET = 20,     /*!< \brief Boundary supersonic inlet definition. */
  ENGINE_INFLOW = 21,         /*!< \brief Boundary nacelle inflow. */
  ENGINE_EXHAUST = 22,        /*!< \brief Boundary nacelle exhaust. */
  RIEMANN_BOUNDARY= 24,       /*!< \brief Riemann Boundary definition. */
  ISOTHERMAL = 25,            /*!< \brief No slip isothermal wall boundary condition. */
  HEAT_FLUX = 26,             /*!< \brief No slip constant heat flux wall boundary condition. */
  HEAT_TRANSFER = 27,         /*!< \brief No slip heat transfer boundary condition. */
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
  SMOLUCHOWSKI_MAXWELL = 55,  /*!< \brief Smoluchoski/Maxwell wall boundary condition. */
  SEND_RECEIVE = 99,          /*!< \brief Boundary send-receive definition. */
};

/*!
 * \brief 2D Formulation for structural problems
 */
enum class STRUCT_2DFORM {
  PLANE_STRESS,     /*!< \brief Definition of plane stress solver. */
  PLANE_STRAIN      /*!< \brief Definition of plane strain solver. */
};
static const MapType<std::string, STRUCT_2DFORM> ElasForm_2D = {
  MakePair("PLANE_STRESS", STRUCT_2DFORM::PLANE_STRESS)
  MakePair("PLANE_STRAIN", STRUCT_2DFORM::PLANE_STRAIN)
};

/*!
 * \brief Kinds of relaxation for multizone problems
 */
enum class BGS_RELAXATION {
  NONE,       /*!< \brief No relaxation in the strongly coupled approach. */
  FIXED,      /*!< \brief Relaxation with a fixed parameter. */
  AITKEN,     /*!< \brief Relaxation using Aitken's dynamic parameter. */
};
static const MapType<std::string, BGS_RELAXATION> AitkenForm_Map = {
  MakePair("NONE", BGS_RELAXATION::NONE)
  MakePair("FIXED_PARAMETER", BGS_RELAXATION::FIXED)
  MakePair("AITKEN_DYNAMIC", BGS_RELAXATION::AITKEN)
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
static const MapType<std::string, ENUM_DYN_TRANSFER_METHOD> Dyn_Transfer_Method_Map = {
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
static const MapType<std::string, ENUM_DVFEA> DVFEA_Map = {
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
enum class RADIATION_MODEL {
  NONE,   /*!< \brief No radiation model */
  P1,     /*!< \brief P1 Radiation model. */
};
static const MapType<std::string, RADIATION_MODEL> Radiation_Map = {
  MakePair("NONE", RADIATION_MODEL::NONE)
  MakePair("P1", RADIATION_MODEL::P1)
};

/*!
 * \brief Kinds of P1 initialization
 */
enum class P1_INIT {
  ZERO,         /*!< \brief Initialize the P1 model from zero values */
  TEMPERATURE,  /*!< \brief Initialize the P1 model from blackbody energy computed from the initial temperature. */
};
static const MapType<std::string, P1_INIT> P1_Init_Map = {
  MakePair("ZERO", P1_INIT::ZERO)
  MakePair("TEMPERATURE_INIT", P1_INIT::TEMPERATURE)
};

/*!
 * \brief Kinds of coupling methods at CHT interfaces.
 * The first (temperature) part determines the BC method on the fluid side, the second (heatflux) part determines
 * the BC method on the solid side of the CHT interface.
 */
enum CHT_COUPLING {
  DIRECT_TEMPERATURE_NEUMANN_HEATFLUX,
  AVERAGED_TEMPERATURE_NEUMANN_HEATFLUX,
  DIRECT_TEMPERATURE_ROBIN_HEATFLUX,
  AVERAGED_TEMPERATURE_ROBIN_HEATFLUX,
};
static const MapType<std::string, CHT_COUPLING> CHT_Coupling_Map = {
  MakePair("DIRECT_TEMPERATURE_NEUMANN_HEATFLUX", CHT_COUPLING::DIRECT_TEMPERATURE_NEUMANN_HEATFLUX)
  MakePair("AVERAGED_TEMPERATURE_NEUMANN_HEATFLUX", CHT_COUPLING::AVERAGED_TEMPERATURE_NEUMANN_HEATFLUX)
  MakePair("DIRECT_TEMPERATURE_ROBIN_HEATFLUX", CHT_COUPLING::DIRECT_TEMPERATURE_ROBIN_HEATFLUX)
  MakePair("AVERAGED_TEMPERATURE_ROBIN_HEATFLUX", CHT_COUPLING::AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)
};

/*!
 * \brief Types Riemann boundary treatments
 */
enum RIEMANN_TYPE {
  TOTAL_CONDITIONS_PT = 1,          /*!< \brief User specifies total pressure, total temperature, and flow direction. */
  DENSITY_VELOCITY = 2,             /*!< \brief User specifies density and velocity, and flow direction. */
  STATIC_PRESSURE = 3,              /*!< \brief User specifies static pressure. */
  TOTAL_SUPERSONIC_INFLOW = 4,      /*!< \brief User specifies total pressure, total temperature and Velocity components. */
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
static const MapType<std::string, RIEMANN_TYPE> Riemann_Map = {
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

static const MapType<std::string, RIEMANN_TYPE> Giles_Map = {
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
  MIXEDOUT = 3,   /*!< \brief an mixed-out average is computed at the boundary of interest. */
  MASSFLUX = 4    /*!< \brief a mass flow average is computed at the boundary of interest. */
};
static const MapType<std::string, AVERAGEPROCESS_TYPE> AverageProcess_Map = {
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
static const MapType<std::string, MIXINGPLANE_INTERFACE_TYPE> MixingPlaneInterface_Map = {
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
static const MapType<std::string, SPANWISE_TYPE> SpanWise_Map = {
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
static const MapType<std::string, TURBOMACHINERY_TYPE> TurboMachinery_Map = {
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
  INFLOW  = 1,    /*!< \brief flag for inflow marker for compute turboperformance. */
  OUTFLOW = 2     /*!< \brief flag for outflow marker for compute turboperformance. */
};

/*!
 * \brief Types inlet boundary treatments
 */
enum class INLET_TYPE {
  TOTAL_CONDITIONS, /*!< \brief User specifies total pressure, total temperature, and flow direction. */
  MASS_FLOW,        /*!< \brief User specifies density and velocity (mass flow). */
  INPUT_FILE,       /*!< \brief User specifies an input file. */
  VELOCITY_INLET,   /*!< \brief Velocity inlet for an incompressible flow. */
  PRESSURE_INLET,   /*!< \brief Total pressure inlet for an incompressible flow. */
};
static const MapType<std::string, INLET_TYPE> Inlet_Map = {
  MakePair("TOTAL_CONDITIONS", INLET_TYPE::TOTAL_CONDITIONS)
  MakePair("MASS_FLOW", INLET_TYPE::MASS_FLOW)
  MakePair("INPUT_FILE", INLET_TYPE::INPUT_FILE)
  MakePair("VELOCITY_INLET", INLET_TYPE::VELOCITY_INLET)
  MakePair("PRESSURE_INLET", INLET_TYPE::PRESSURE_INLET)
};

/*!
 * \brief Types outlet boundary treatments
 */
enum class INC_OUTLET_TYPE {
  PRESSURE_OUTLET,    /*!< \brief Gauge pressure outlet for incompressible flow */
  MASS_FLOW_OUTLET,   /*!< \brief Mass flow outlet for incompressible flow. */
};
static const MapType<std::string, INC_OUTLET_TYPE> Inc_Outlet_Map = {
  MakePair("PRESSURE_OUTLET",  INC_OUTLET_TYPE::PRESSURE_OUTLET)
  MakePair("MASS_FLOW_OUTLET", INC_OUTLET_TYPE::MASS_FLOW_OUTLET)
};

/*!
 * \brief Types engine inflow boundary treatments
 */
enum ENGINE_INFLOW_TYPE {
  FAN_FACE_MACH = 1,          /*!< \brief User specifies fan face mach number. */
  FAN_FACE_MDOT = 2,          /*!< \brief User specifies Static pressure. */
  FAN_FACE_PRESSURE = 3       /*!< \brief User specifies Static pressure. */
};
static const MapType<std::string, ENGINE_INFLOW_TYPE> Engine_Inflow_Map = {
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
  POWER = 6,              /*!< \brief User specifies the power. */
  VARIABLE_LOAD = 7       /*!< \brief User specifies the load distribution. */
};
static const MapType<std::string, ACTDISK_TYPE> ActDisk_Map = {
  MakePair("VARIABLES_JUMP", VARIABLES_JUMP)
  MakePair("BC_THRUST", BC_THRUST)
  MakePair("NET_THRUST", NET_THRUST)
  MakePair("DRAG_MINUS_THRUST", DRAG_MINUS_THRUST)
  MakePair("MASSFLOW", MASSFLOW)
  MakePair("POWER", POWER)
  MakePair("VARIABLE_LOAD", VARIABLE_LOAD)
};

/*!
 * \brief types of wall boundary condition - smooth or rough
 */
enum class WALL_TYPE {
  SMOOTH,  /*!< \brief Smooth wall */
  ROUGH,   /*!< \brief Rough wall */
};
static const MapType<std::string, WALL_TYPE> WallType_Map = {
  MakePair("SMOOTH", WALL_TYPE::SMOOTH)
  MakePair("ROUGH", WALL_TYPE::ROUGH)
};

/*!
 * \brief Types of objective functions
 */
enum ENUM_OBJECTIVE {
  DRAG_COEFFICIENT = 1,         /*!< \brief Drag objective function definition. */
  LIFT_COEFFICIENT = 2,         /*!< \brief Lift objective function definition. */
  SIDEFORCE_COEFFICIENT = 3,    /*!< \brief Side force objective function definition. */
  EFFICIENCY = 4,               /*!< \brief Efficiency objective function definition. */
  INVERSE_DESIGN_PRESSURE = 5,  /*!< \brief Pressure objective function definition (inverse design). */
  INVERSE_DESIGN_HEATFLUX = 6,  /*!< \brief Heat flux objective function definition (inverse design). */
  TOTAL_HEATFLUX = 7,           /*!< \brief Total heat flux. */
  MAXIMUM_HEATFLUX = 8,         /*!< \brief Maximum heat flux. */
  AVG_TEMPERATURE = 70,         /*!< \brief Total averaged temperature. */
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
  SURFACE_STATIC_TEMPERATURE = 57, /*!< \brief Static Temperature objective function definition. */
  SURFACE_MASSFLOW = 30,        /*!< \brief Mass Flow Rate objective function definition. */
  SURFACE_MACH = 51,            /*!< \brief Mach number objective function definition. */
  SURFACE_UNIFORMITY = 52,      /*!< \brief Flow uniformity objective function definition. */
  SURFACE_SECONDARY = 53,       /*!< \brief Secondary flow strength objective function definition. */
  SURFACE_MOM_DISTORTION = 54,  /*!< \brief Momentum distortion objective function definition. */
  SURFACE_SECOND_OVER_UNIFORM = 55, /*!< \brief Secondary over uniformity (relative secondary strength) objective function definition. */
  SURFACE_PRESSURE_DROP = 56,   /*!< \brief Pressure drop objective function definition. */
  SURFACE_SPECIES_0 = 58,       /*!< \brief Surface Avg. Species_0 objective function definition. */
  SURFACE_SPECIES_VARIANCE = 59,/*!< \brief Species Variance objective function definition. */
  CUSTOM_OBJFUNC = 31,          /*!< \brief Custom objective function definition. */
  FLOW_ANGLE_OUT = 46,
  MASS_FLOW_IN = 47,
  ENTROPY_GENERATION = 50,
  REFERENCE_GEOMETRY = 60,      /*!< \brief Norm of displacements with respect to target geometry. */
  REFERENCE_NODE = 61,          /*!< \brief Objective function defined as the difference of a particular node respect to a reference position. */
  VOLUME_FRACTION = 62,         /*!< \brief Volume average physical density, for material-based topology optimization applications. */
  TOPOL_DISCRETENESS = 63,      /*!< \brief Measure of the discreteness of the current topology. */
  TOPOL_COMPLIANCE = 64,        /*!< \brief Measure of the discreteness of the current topology. */
  STRESS_PENALTY = 65,          /*!< \brief Penalty function of VM stresses above a maximum value. */
};
static const MapType<std::string, ENUM_OBJECTIVE> Objective_Map = {
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
  MakePair("AVG_TEMPERATURE", AVG_TEMPERATURE)
  MakePair("FIGURE_OF_MERIT", FIGURE_OF_MERIT)
  MakePair("BUFFET", BUFFET_SENSOR)
  MakePair("SURFACE_TOTAL_PRESSURE", SURFACE_TOTAL_PRESSURE)
  MakePair("SURFACE_STATIC_PRESSURE", SURFACE_STATIC_PRESSURE)
  MakePair("SURFACE_STATIC_TEMPERATURE", SURFACE_STATIC_TEMPERATURE)
  MakePair("SURFACE_MASSFLOW", SURFACE_MASSFLOW)
  MakePair("SURFACE_MACH", SURFACE_MACH)
  MakePair("SURFACE_UNIFORMITY", SURFACE_UNIFORMITY)
  MakePair("SURFACE_SECONDARY", SURFACE_SECONDARY)
  MakePair("SURFACE_MOM_DISTORTION", SURFACE_MOM_DISTORTION)
  MakePair("SURFACE_SECOND_OVER_UNIFORM", SURFACE_SECOND_OVER_UNIFORM)
  MakePair("SURFACE_PRESSURE_DROP", SURFACE_PRESSURE_DROP)
  MakePair("SURFACE_SPECIES_0", SURFACE_SPECIES_0)
  MakePair("SURFACE_SPECIES_VARIANCE", SURFACE_SPECIES_VARIANCE)
  MakePair("CUSTOM_OBJFUNC", CUSTOM_OBJFUNC)
  MakePair("FLOW_ANGLE_OUT", FLOW_ANGLE_OUT)
  MakePair("MASS_FLOW_IN", MASS_FLOW_IN)
  MakePair("ENTROPY_GENERATION",  ENTROPY_GENERATION)
  MakePair("REFERENCE_GEOMETRY", REFERENCE_GEOMETRY)
  MakePair("REFERENCE_NODE", REFERENCE_NODE)
  MakePair("VOLUME_FRACTION", VOLUME_FRACTION)
  MakePair("TOPOL_DISCRETENESS", TOPOL_DISCRETENESS)
  MakePair("TOPOL_COMPLIANCE", TOPOL_COMPLIANCE)
  MakePair("STRESS_PENALTY", STRESS_PENALTY)
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
static const MapType<std::string, ENUM_INPUT> Input_Map = {
  MakePair("SU2", SU2)
  MakePair("CGNS", CGNS_GRID)
  MakePair("RECTANGLE", RECTANGLE)
  MakePair("BOX", BOX)
};


/*!
 * \brief Type of solution output file formats
 */
enum class OUTPUT_TYPE {
  TECPLOT_ASCII,           /*!< \brief Tecplot format for the solution output. */
  TECPLOT_BINARY,          /*!< \brief Tecplot binary format for the solution output. */
  SURFACE_TECPLOT_ASCII,   /*!< \brief Tecplot format for the solution output. */
  SURFACE_TECPLOT_BINARY,  /*!< \brief Tecplot binary format for the solution output. */
  CSV,                     /*!< \brief Comma-separated values format for the solution output. */
  SURFACE_CSV,             /*!< \brief Comma-separated values format for the solution output. */
  PARAVIEW_ASCII,          /*!< \brief Paraview ASCII format for the solution output. */
  PARAVIEW_LEGACY_BINARY,  /*!< \brief Paraview binary format for the solution output. */
  SURFACE_PARAVIEW_ASCII,  /*!< \brief Paraview ASCII format for the solution output. */
  SURFACE_PARAVIEW_LEGACY_BINARY, /*!< \brief Paraview binary format for the solution output. */
  MESH,                    /*!< \brief SU2 mesh format. */
  RESTART_BINARY,          /*!< \brief SU2 binary restart format. */
  RESTART_ASCII,           /*!< \brief SU2 ASCII restart format. */
  PARAVIEW_XML,            /*!< \brief Paraview XML with binary data format */
  SURFACE_PARAVIEW_XML,    /*!< \brief Surface Paraview XML with binary data format */
  PARAVIEW_MULTIBLOCK,     /*!< \brief Paraview XML Multiblock */
  CGNS,                    /*!< \brief CGNS format. */
  SURFACE_CGNS,            /*!< \brief CGNS format. */
  STL_ASCII,               /*!< \brief STL ASCII format for surface solution output. */
  STL_BINARY,              /*!< \brief STL binary format for surface solution output. Not implemented yet. */
};
static const MapType<std::string, OUTPUT_TYPE> Output_Map = {
  MakePair("TECPLOT_ASCII", OUTPUT_TYPE::TECPLOT_ASCII)
  MakePair("TECPLOT", OUTPUT_TYPE::TECPLOT_BINARY)
  MakePair("SURFACE_TECPLOT_ASCII", OUTPUT_TYPE::SURFACE_TECPLOT_ASCII)
  MakePair("SURFACE_TECPLOT", OUTPUT_TYPE::SURFACE_TECPLOT_BINARY)
  MakePair("CSV", OUTPUT_TYPE::CSV)
  MakePair("SURFACE_CSV", OUTPUT_TYPE::SURFACE_CSV)
  MakePair("PARAVIEW_ASCII", OUTPUT_TYPE::PARAVIEW_ASCII)
  MakePair("PARAVIEW_LEGACY", OUTPUT_TYPE::PARAVIEW_LEGACY_BINARY)
  MakePair("SURFACE_PARAVIEW_ASCII", OUTPUT_TYPE::SURFACE_PARAVIEW_ASCII)
  MakePair("SURFACE_PARAVIEW_LEGACY", OUTPUT_TYPE::SURFACE_PARAVIEW_LEGACY_BINARY)
  MakePair("PARAVIEW", OUTPUT_TYPE::PARAVIEW_XML)
  MakePair("SURFACE_PARAVIEW", OUTPUT_TYPE::SURFACE_PARAVIEW_XML)
  MakePair("PARAVIEW_MULTIBLOCK", OUTPUT_TYPE::PARAVIEW_MULTIBLOCK)
  MakePair("MESH", OUTPUT_TYPE::MESH)
  MakePair("RESTART_ASCII", OUTPUT_TYPE::RESTART_ASCII)
  MakePair("RESTART", OUTPUT_TYPE::RESTART_BINARY)
  MakePair("CGNS", OUTPUT_TYPE::CGNS)
  MakePair("SURFACE_CGNS", OUTPUT_TYPE::SURFACE_CGNS)
  MakePair("STL_ASCII", OUTPUT_TYPE::STL_ASCII)
  MakePair("STL_BINARY", OUTPUT_TYPE::STL_BINARY)
};

/*!
 * \brief Return true if format is one of the Paraview options.
 */
inline bool isParaview(OUTPUT_TYPE format) {
  switch(format) {
    case OUTPUT_TYPE::PARAVIEW_ASCII:
    case OUTPUT_TYPE::PARAVIEW_LEGACY_BINARY:
    case OUTPUT_TYPE::SURFACE_PARAVIEW_ASCII:
    case OUTPUT_TYPE::SURFACE_PARAVIEW_LEGACY_BINARY:
    case OUTPUT_TYPE::PARAVIEW_XML:
    case OUTPUT_TYPE::SURFACE_PARAVIEW_XML:
    case OUTPUT_TYPE::PARAVIEW_MULTIBLOCK:
      return true;
    default:
      return false;
  }
}

/*!
 * \brief Return true if format is one of the Tecplot options.
 */
inline bool isTecplot(OUTPUT_TYPE format) {
  switch(format) {
    case OUTPUT_TYPE::TECPLOT_ASCII:
    case OUTPUT_TYPE::TECPLOT_BINARY:
    case OUTPUT_TYPE::SURFACE_TECPLOT_ASCII:
    case OUTPUT_TYPE::SURFACE_TECPLOT_BINARY:
      return true;
    default:
      return false;
  }
}

/*!
 * \brief Type of solution output file formats
 */
enum class TAB_OUTPUT {
  TAB_CSV,            /*!< \brief Comma-separated values format for the solution output. */
  TAB_TECPLOT         /*!< \brief Tecplot format for the solution output. */
};
static const MapType<std::string, TAB_OUTPUT> TabOutput_Map = {
  MakePair("CSV", TAB_OUTPUT::TAB_CSV)
  MakePair("TECPLOT", TAB_OUTPUT::TAB_TECPLOT)
};

/*!
 * \brief Type of volume sensitivity file formats (inout to SU2_DOT)
 */
enum ENUM_SENSITIVITY {
  SU2_NATIVE = 1,       /*!< \brief SU2 native binary format for the volume sensitivity input. */
  UNORDERED_ASCII = 2   /*!< \brief Unordered ASCII list (x,y,z,dJ/dx,dJ/dy/dJ/dz) format for the volume sensitivity input. */
};
static const MapType<std::string, ENUM_SENSITIVITY> Sensitivity_Map = {
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
static const MapType<std::string, JUMP_DEFINITION> Jump_Map = {
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
static const MapType<std::string, MG_CYCLE> MG_Cycle_Map = {
  MakePair("V_CYCLE", V_CYCLE)
  MakePair("W_CYCLE", W_CYCLE)
  MakePair("FULLMG_CYCLE", FULLMG_CYCLE)
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
static const MapType<std::string, ENUM_PARAM> Param_Map = {
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
static const MapType<std::string, ENUM_FFD_BLENDING> Blending_Map = {
  MakePair("BSPLINE_UNIFORM", BSPLINE_UNIFORM)
  MakePair("BEZIER", BEZIER)
};

/*!
 * \brief Types of solvers for solving linear systems
 */
enum ENUM_LINEAR_SOLVER {
  CONJUGATE_GRADIENT,   /*!< \brief Preconditionated conjugate gradient method for grid deformation. */
  FGMRES,               /*!< \brief Flexible Generalized Minimal Residual method. */
  BCGSTAB,              /*!< \brief BCGSTAB - Biconjugate Gradient Stabilized Method (main solver). */
  RESTARTED_FGMRES,     /*!< \brief Flexible Generalized Minimal Residual method with restart. */
  SMOOTHER,             /*!< \brief Iterative smoother. */
  PASTIX_LDLT,          /*!< \brief PaStiX LDLT (complete) factorization. */
  PASTIX_LU,            /*!< \brief PaStiX LU (complete) factorization. */
};
static const MapType<std::string, ENUM_LINEAR_SOLVER> Linear_Solver_Map = {
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
static const MapType<std::string, ENUM_FFD_CONTINUITY> Continuity_Map = {
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
static const MapType<std::string, ENUM_FFD_COORD_SYSTEM> CoordSystem_Map = {
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
static const MapType<std::string, ENUM_SENS_SMOOTHING> Sens_Smoothing_Map = {
  MakePair("NONE", NO_SMOOTH)
  MakePair("SOBOLEV", SOBOLEV)
  MakePair("BIGRID", BIGRID)
};

/*!
 * \brief Types of preconditioners for the linear solver
 */
enum ENUM_LINEAR_SOLVER_PREC {
  JACOBI,         /*!< \brief Jacobi preconditioner. */
  LU_SGS,         /*!< \brief LU SGS preconditioner. */
  LINELET,        /*!< \brief Line implicit preconditioner. */
  ILU,            /*!< \brief ILU(k) preconditioner. */
  PASTIX_ILU=10,  /*!< \brief PaStiX ILU(k) preconditioner. */
  PASTIX_LU_P,    /*!< \brief PaStiX LU as preconditioner. */
  PASTIX_LDLT_P,  /*!< \brief PaStiX LDLT as preconditioner. */
};
static const MapType<std::string, ENUM_LINEAR_SOLVER_PREC> Linear_Solver_Prec_Map = {
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
static const MapType<std::string, ENUM_GEO_ANALYTIC> Geo_Analytic_Map = {
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
  WING = 1,         /*!< \brief Wing analysis. */
  FUSELAGE = 2,     /*!< \brief Fuselage analysis. */
  NACELLE = 3       /*!< \brief Nacelle analysis. */
};
static const MapType<std::string, ENUM_GEO_DESCRIPTION> Geo_Description_Map = {
  MakePair("AIRFOIL", TWOD_AIRFOIL)
  MakePair("WING", WING)
  MakePair("FUSELAGE", FUSELAGE)
  MakePair("NACELLE", NACELLE)
};

/*!
 * \brief Types of schemes for unsteady computations
 */
enum class TIME_MARCHING {
  STEADY,           /*!< \brief A steady computation. */
  TIME_STEPPING,    /*!< \brief Use a time stepping strategy for unsteady computations. */
  DT_STEPPING_1ST,  /*!< \brief Use a dual time stepping strategy for unsteady computations (1st order). */
  DT_STEPPING_2ND,  /*!< \brief Use a dual time stepping strategy for unsteady computations (2nd order). */
  ROTATIONAL_FRAME, /*!< \brief Use a rotational source term. */
  HARMONIC_BALANCE, /*!< \brief Use a harmonic balance source term. */
};
static const MapType<std::string, TIME_MARCHING> TimeMarching_Map = {
  MakePair("NO", TIME_MARCHING::STEADY)
  MakePair("TIME_STEPPING", TIME_MARCHING::TIME_STEPPING)
  MakePair("DUAL_TIME_STEPPING-1ST_ORDER", TIME_MARCHING::DT_STEPPING_1ST)
  MakePair("DUAL_TIME_STEPPING-2ND_ORDER", TIME_MARCHING::DT_STEPPING_2ND)
  MakePair("HARMONIC_BALANCE", TIME_MARCHING::HARMONIC_BALANCE)
  MakePair("ROTATIONAL_FRAME", TIME_MARCHING::ROTATIONAL_FRAME)
};

/*!
 * \brief Types of element stiffnesses imposed for FEA mesh deformation
 */
enum ENUM_DEFORM_STIFFNESS {
  CONSTANT_STIFFNESS = 0,     /*!< \brief Impose a constant stiffness for each element (steel). */
  INVERSE_VOLUME = 1,         /*!< \brief Impose a stiffness for each element that is inversely proportional to cell volume. */
  SOLID_WALL_DISTANCE = 2     /*!< \brief Impose a stiffness for each element that is proportional to the distance from the solid surface. */
};
static const MapType<std::string, ENUM_DEFORM_STIFFNESS> Deform_Stiffness_Map = {
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
static const MapType<std::string, ENUM_DIRECTDIFF_VAR> DirectDiff_Var_Map = {
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


enum class RECORDING {
  CLEAR_INDICES,
  SOLUTION_VARIABLES,
  MESH_COORDS,
  MESH_DEFORM,
  SOLUTION_AND_MESH,
};

/*!
 * \brief Types of schemes for dynamic structural computations
 */
enum ENUM_DYNAMIC {
  STATIC = 0,     /*!< \brief A static structural computation. */
  DYNAMIC = 1     /*!< \brief Use a time stepping strategy for dynamic computations. */
};
static const MapType<std::string, ENUM_DYNAMIC> Dynamic_Map = {
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
static const MapType<std::string, ENUM_INPUT_REF> Input_Ref_Map = {
  MakePair("SU2", SU2_REF)
  MakePair("CUSTOM", CUSTOM_REF)
};

/*!
 * \brief Vertex-based quantities exchanged during periodic marker communications.
 */
enum PERIODIC_QUANTITIES {
  PERIODIC_NONE       ,  /*!< \brief No periodic communication required. */
  PERIODIC_VOLUME     ,  /*!< \brief Volume communication for summing total CV (periodic only). */
  PERIODIC_NEIGHBORS  ,  /*!< \brief Communication of the number of neighbors for centered schemes (periodic only). */
  PERIODIC_RESIDUAL   ,  /*!< \brief Residual and Jacobian communication (periodic only). */
  PERIODIC_LAPLACIAN  ,  /*!< \brief Undivided Laplacian communication for JST (periodic only). */
  PERIODIC_MAX_EIG    ,  /*!< \brief Maximum eigenvalue communication (periodic only). */
  PERIODIC_SENSOR     ,  /*!< \brief Dissipation sensor communication (periodic only). */
  PERIODIC_SOL_GG     ,  /*!< \brief Solution gradient communication for Green-Gauss (periodic only). */
  PERIODIC_PRIM_GG    ,  /*!< \brief Primitive gradient communication for Green-Gauss (periodic only). */
  PERIODIC_SOL_LS     ,  /*!< \brief Solution gradient communication for weighted Least Squares (periodic only). */
  PERIODIC_PRIM_LS    ,  /*!< \brief Primitive gradient communication for weighted Least Squares (periodic only). */
  PERIODIC_SOL_ULS    ,  /*!< \brief Solution gradient communication for unwieghted Least Squares (periodic only). */
  PERIODIC_PRIM_ULS   ,  /*!< \brief Primitive gradient communication for unweighted Least Squares (periodic only). */
  PERIODIC_SOL_GG_R   ,  /*!< \brief Same but reconstruction. */
  PERIODIC_PRIM_GG_R  ,  /*!< \brief Same but reconstruction. */
  PERIODIC_SOL_LS_R   ,  /*!< \brief Same but reconstruction. */
  PERIODIC_PRIM_LS_R  ,  /*!< \brief Same but reconstruction. */
  PERIODIC_SOL_ULS_R  ,  /*!< \brief Same but reconstruction. */
  PERIODIC_PRIM_ULS_R ,  /*!< \brief Same but reconstruction. */
  PERIODIC_LIM_SOL_1  ,  /*!< \brief Solution limiter communication phase 1 of 2 (periodic only). */
  PERIODIC_LIM_SOL_2  ,  /*!< \brief Solution limiter communication phase 2 of 2 (periodic only). */
  PERIODIC_LIM_PRIM_1 ,  /*!< \brief Primitive limiter communication phase 1 of 2 (periodic only). */
  PERIODIC_LIM_PRIM_2 ,  /*!< \brief Primitive limiter communication phase 2 of 2 (periodic only). */
  PERIODIC_IMPLICIT   ,  /*!< \brief Implicit update communication to ensure consistency across periodic boundaries. */
};

/*!
 * \brief Vertex-based quantities exchanged in MPI point-to-point communications.
 */
enum MPI_QUANTITIES {
  SOLUTION             ,  /*!< \brief Conservative solution communication. */
  SOLUTION_OLD         ,  /*!< \brief Conservative solution old communication. */
  SOLUTION_GRADIENT    ,  /*!< \brief Conservative solution gradient communication. */
  SOLUTION_GRAD_REC    ,  /*!< \brief Conservative solution reconstruction gradient communication. */
  SOLUTION_LIMITER     ,  /*!< \brief Conservative solution limiter communication. */
  SOLUTION_GEOMETRY    ,  /*!< \brief Geometry solution communication. */
  PRIMITIVE_GRADIENT   ,  /*!< \brief Primitive gradient communication. */
  PRIMITIVE_GRAD_REC   ,  /*!< \brief Primitive reconstruction gradient communication. */
  PRIMITIVE_LIMITER    ,  /*!< \brief Primitive limiter communication. */
  UNDIVIDED_LAPLACIAN  ,  /*!< \brief Undivided Laplacian communication. */
  MAX_EIGENVALUE       ,  /*!< \brief Maximum eigenvalue communication. */
  SENSOR               ,  /*!< \brief Dissipation sensor communication. */
  AUXVAR_GRADIENT      ,  /*!< \brief Auxiliary variable gradient communication. */
  COORDINATES          ,  /*!< \brief Vertex coordinates communication. */
  COORDINATES_OLD      ,  /*!< \brief Old vertex coordinates communication. */
  MAX_LENGTH           ,  /*!< \brief Maximum length communication. */
  GRID_VELOCITY        ,  /*!< \brief Grid velocity communication. */
  SOLUTION_EDDY        ,  /*!< \brief Turbulent solution plus eddy viscosity communication. */
  SOLUTION_MATRIX      ,  /*!< \brief Matrix solution communication. */
  SOLUTION_MATRIXTRANS ,  /*!< \brief Matrix transposed solution communication. */
  NEIGHBORS            ,  /*!< \brief Neighbor point count communication (for JST). */
  SOLUTION_FEA         ,  /*!< \brief FEA solution communication. */
  MESH_DISPLACEMENTS   ,  /*!< \brief Mesh displacements at the interface. */
  SOLUTION_TIME_N      ,  /*!< \brief Solution at time n. */
  SOLUTION_TIME_N1     ,  /*!< \brief Solution at time n-1. */
};

/*!
 * \brief MPI communication level
 */
enum COMM_LEVEL {
  COMM_NONE    = 0,   /*!< \brief Disable all MPI comms. Purely for testing, as results are incorrect. */
  COMM_MINIMAL = 1,   /*!< \brief Perform only the minimal set of MPI communications for correctness. Disables many console and output comms. */
  COMM_FULL    = 2    /*!< \brief Perform all MPI communications. */
};
static const MapType<std::string, COMM_LEVEL> Comm_Map = {
  MakePair("NONE",    COMM_NONE)
  MakePair("MINIMAL", COMM_MINIMAL)
  MakePair("FULL",    COMM_FULL)
};

/*!
 * \brief Types of filter kernels, initially intended for structural topology optimization applications
 */
enum class ENUM_FILTER_KERNEL {
  CONSTANT_WEIGHT,  /*!< \brief Uniform weight. */
  CONICAL_WEIGHT,   /*!< \brief Linear decay with distance from center point [Bruns and Tortorelli, 2001]. */
  GAUSSIAN_WEIGHT,  /*!< \brief Bell shape around center point [Bruns and Tortorelli, 2003]. */
  DILATE_MORPH,     /*!< \brief Continuous version of the dilate morphology operator [Sigmund 2007]. */
  ERODE_MORPH,      /*!< \brief Continuous version of the erode morphology operator [Sigmund 2007].*/
};
static const MapType<std::string, ENUM_FILTER_KERNEL> Filter_Kernel_Map = {
  MakePair("CONSTANT", ENUM_FILTER_KERNEL::CONSTANT_WEIGHT)
  MakePair("CONICAL", ENUM_FILTER_KERNEL::CONICAL_WEIGHT)
  MakePair("GAUSSIAN", ENUM_FILTER_KERNEL::GAUSSIAN_WEIGHT)
  MakePair("DILATE", ENUM_FILTER_KERNEL::DILATE_MORPH)
  MakePair("ERODE", ENUM_FILTER_KERNEL::ERODE_MORPH)
};

/*!
 * \brief Types of projection function, initially intended for structural topology optimization applications
 */
enum class ENUM_PROJECTION_FUNCTION {
  NONE,           /*!< \brief No projection. */
  HEAVISIDE_UP,   /*!< \brief Project values towards 1. */
  HEAVISIDE_DOWN, /*!< \brief Project values towards 0. */
};
static const MapType<std::string, ENUM_PROJECTION_FUNCTION> Projection_Function_Map = {
  MakePair("NO_PROJECTION", ENUM_PROJECTION_FUNCTION::NONE)
  MakePair("HEAVISIDE_UP", ENUM_PROJECTION_FUNCTION::HEAVISIDE_UP)
  MakePair("HEAVISIDE_DOWN", ENUM_PROJECTION_FUNCTION::HEAVISIDE_DOWN)
};

/*!
 * \brief the different validation solution
 */
enum class VERIFICATION_SOLUTION {
  NONE,                     /*!< \brief No verification solution, standard solver mode. */
  INVISCID_VORTEX,          /*!< \brief Inviscid vortex. Exact solution of the unsteady Euler equations. */
  RINGLEB,                  /*!< \brief Ringleb flow. Exact solution of the steady Euler equations. */
  NS_UNIT_QUAD,             /*!< \brief Exact solution of the laminar Navier Stokes equations without heat conduction. */
  TAYLOR_GREEN_VORTEX,      /*!< \brief Taylor Green Vortex. */
  INC_TAYLOR_GREEN_VORTEX,  /*!< \brief Incompressible Taylor Green Vortex (2D). */
  MMS_NS_UNIT_QUAD,         /*!< \brief Manufactured solution of the laminar Navier Stokes equations on a unit quad. */
  MMS_NS_UNIT_QUAD_WALL_BC, /*!< \brief Manufactured solution of the laminar Navier Stokes equations on a unit quad with wall BC's. */
  MMS_NS_TWO_HALF_CIRCLES,  /*!< \brief Manufactured solution of the laminar Navier Stokes equations between two half circles. */
  MMS_NS_TWO_HALF_SPHERES,  /*!< \brief Manufactured solution of the laminar Navier Stokes equations between two half spheres. */
  MMS_INC_EULER,            /*!< \brief Manufactured solution of the incompressible Euler equations. */
  MMS_INC_NS,               /*!< \brief Manufactured solution of the laminar incompressible Navier Stokes equations. */
  USER_DEFINED_SOLUTION,    /*!< \brief User defined solution. */
};
static const MapType<std::string, VERIFICATION_SOLUTION> Verification_Solution_Map = {
  MakePair("NO_VERIFICATION_SOLUTION", VERIFICATION_SOLUTION::NONE)
  MakePair("INVISCID_VORTEX",          VERIFICATION_SOLUTION::INVISCID_VORTEX)
  MakePair("RINGLEB",                  VERIFICATION_SOLUTION::RINGLEB)
  MakePair("NS_UNIT_QUAD",             VERIFICATION_SOLUTION::NS_UNIT_QUAD)
  MakePair("TAYLOR_GREEN_VORTEX",      VERIFICATION_SOLUTION::TAYLOR_GREEN_VORTEX)
  MakePair("INC_TAYLOR_GREEN_VORTEX",  VERIFICATION_SOLUTION::INC_TAYLOR_GREEN_VORTEX)
  MakePair("MMS_NS_UNIT_QUAD",         VERIFICATION_SOLUTION::MMS_NS_UNIT_QUAD)
  MakePair("MMS_NS_UNIT_QUAD_WALL_BC", VERIFICATION_SOLUTION::MMS_NS_UNIT_QUAD_WALL_BC)
  MakePair("MMS_NS_TWO_HALF_CIRCLES",  VERIFICATION_SOLUTION::MMS_NS_TWO_HALF_CIRCLES)
  MakePair("MMS_NS_TWO_HALF_SPHERES",  VERIFICATION_SOLUTION::MMS_NS_TWO_HALF_SPHERES)
  MakePair("MMS_INC_EULER",            VERIFICATION_SOLUTION::MMS_INC_EULER)
  MakePair("MMS_INC_NS",               VERIFICATION_SOLUTION::MMS_INC_NS)
  MakePair("USER_DEFINED_SOLUTION",    VERIFICATION_SOLUTION::USER_DEFINED_SOLUTION)
};

/*!
 * \brief Types of streamwise periodicity.
 */
enum class ENUM_STREAMWISE_PERIODIC {
  NONE,          /*!< \brief No streamwise periodic flow. */
  PRESSURE_DROP, /*!< \brief Prescribed pressure drop. */
  MASSFLOW,      /*!< \brief Prescribed massflow. */
};
static const MapType<std::string, ENUM_STREAMWISE_PERIODIC> Streamwise_Periodic_Map = {
  MakePair("NONE",          ENUM_STREAMWISE_PERIODIC::NONE)
  MakePair("PRESSURE_DROP", ENUM_STREAMWISE_PERIODIC::PRESSURE_DROP)
  MakePair("MASSFLOW",      ENUM_STREAMWISE_PERIODIC::MASSFLOW)
};

/*!
 * \brief Container to hold Variables for streamwise Periodic flow as they are often used together in places.
 */
struct StreamwisePeriodicValues {
  su2double Streamwise_Periodic_PressureDrop;       /*!< \brief Value of prescribed pressure drop [Pa] which results in an artificial body force vector. */
  su2double Streamwise_Periodic_MassFlow;           /*!< \brief Value of current massflow [kg/s] which results in a delta p and therefore an artificial body force vector. */
  su2double Streamwise_Periodic_IntegratedHeatFlow; /*!< \brief Value of of the net sum of heatflow [W] into the domain. */
  su2double Streamwise_Periodic_InletTemperature;   /*!< \brief Area avg static Temp [K] at the periodic inlet. Used for adaptive outlet heatsink. */
  su2double Streamwise_Periodic_BoundaryArea;       /*!< \brief Global Surface area of the streamwise periodic interface. */
  su2double Streamwise_Periodic_AvgDensity;         /*!< \brief Area avg density on the periodic interface. */
};

/*!
 * \brief Type of POD basis generation (for use with libROM)
 */
enum class POD_KIND {
  STATIC,            /*!< \brief Use static SVD for POD basis generation. */
  INCREMENTAL,       /*!< \brief Use incremental SVD for POD basis generation. */
};
static const MapType<std::string, POD_KIND> POD_Map = {
  MakePair("STATIC_POD",      POD_KIND::STATIC)
  MakePair("INCREMENTAL_POD", POD_KIND::INCREMENTAL)
};

/*!
 * \brief Type of operation for the linear system solver, changes the source of solver options.
 */
enum class LINEAR_SOLVER_MODE {
  STANDARD,        /*!< \brief Operate in standard mode. */
  MESH_DEFORM,     /*!< \brief Operate in mesh deformation mode. */
  GRADIENT_MODE,   /*!< \brief Operate in gradient smoothing mode. */
};

/*!
 * \brief mode of operation for the sobolev smoothing solver.
 */
enum class ENUM_SOBOLEV_MODUS {
  NONE,                 /*!< \brief Default option if none is choosen. */
  PARAM_LEVEL_COMPLETE, /*!< \brief Operate on parameter level. */
  MESH_LEVEL,           /*!< \brief Operate on mesh level. */
  ONLY_GRAD,            /*!< \brief Flag to only compute the original gradient. */
};
static const MapType<std::string, ENUM_SOBOLEV_MODUS> Sobolev_Modus_Map = {
  MakePair("NONE",                 ENUM_SOBOLEV_MODUS::NONE)
  MakePair("PARAM_LEVEL_COMPLETE", ENUM_SOBOLEV_MODUS::PARAM_LEVEL_COMPLETE)
  MakePair("MESH_LEVEL",           ENUM_SOBOLEV_MODUS::MESH_LEVEL)
  MakePair("ONLY_GRADIENT",        ENUM_SOBOLEV_MODUS::ONLY_GRAD)
};

#undef MakePair
/* END_CONFIG_ENUMS */

class COptionBase {
private:
  std::vector<std::string> value;
public:
  virtual ~COptionBase() = default;

  const std::vector<std::string>& GetValue() const {return value;}

  virtual std::string SetValue(const std::vector<std::string>& val) {
    value = val;
    return "";
  }
  virtual void SetDefault() = 0;

  std::string optionCheckMultipleValues(const std::vector<std::string>& option_value,
                                        std::string type_id, const std::string& option_name) {
    if (option_value.size() != 1) {
      std::string newString(option_name);
      newString.append(": multiple values for type ");
      newString.append(type_id);
      return newString;
    }
    return "";
  }

  std::string badValue(std::string type_id, const std::string& option_name) {
    std::string newString(option_name);
    newString.append(": improper option value for type ");
    newString.append(type_id);
    return newString;
  }
};

#ifdef ENABLE_MAPS
#include "option_structure.inl"
#endif
