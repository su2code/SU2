/*!
 * \file option_structure.hpp
 * \brief Defines classes for referencing options for easy input in CConfig
 * \author J. Hicken, B. Tracey
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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
 * \class CCreateMap
 * \brief creates a map from a list by overloading operator()
 * \tparam T - the key type in the map
 * \tparam U - the mapped value type in the map
 * \author Boost.Assign and anonymous person on stackoverflow
 *
 * We need this to create static const maps that map strings to enum
 * types.  The implementation is based on the Boost.Assign library.  This
 * particular version is taken from
 * http://stackoverflow.com/questions/138600/initializing-a-static-stdmapint-int-in-c
 */
template <typename T, typename U>
class CCreateMap {
private:
  std::map<T, U> m_map;
public:
  CCreateMap(const T& key, const U& val) {
    m_map[key] = val;
  }
  CCreateMap<T, U>& operator()(const T& key, const U& val) {
    m_map[key] = val;
    return *this;
  }
  operator std::map<T, U>() {
    return m_map;
  }
};

/*!
 * \brief utility function for converting strings to uppercase
 * \param[in, out] str - string we want to convert
 */
inline void StringToUpperCase(string & str) {
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

/*!
 * \brief utility function for converting strings to uppercase
 * \param[in] str - string we want a copy of converted to uppercase
 * \returns a copy of str in uppercase
 */
inline string StringToUpperCase(const string & str) {
  string upp_str(str);
  std::transform(upp_str.begin(), upp_str.end(), upp_str.begin(), ::toupper);
  return upp_str;
}

/*!
 * \brief different software components of SU2
 */
enum SU2_COMPONENT {
  SU2_CFD = 1,	/*!< \brief Running the SU2_CFD software. */
  SU2_DEF = 2,	/*!< \brief Running the SU2_DEF software. */
  SU2_DOT = 3,	/*!< \brief Running the SU2_DOT software. */
  SU2_MSH = 4,	/*!< \brief Running the SU2_MSH software. */
  SU2_GEO = 5,	/*!< \brief Running the SU2_GEO software. */
  SU2_SOL = 6 	/*!< \brief Running the SU2_SOL software. */
};

const unsigned int EXIT_DIVERGENCE = 2;       /*!< \brief Exit code (divergence). */

const unsigned int BUFSIZE = 3000000;		  /*!< \brief MPI buffer. */
const unsigned int MAX_PARAMETERS = 10;		  /*!< \brief Maximum number of parameters for a design variable definition. */
const unsigned int MAX_NUMBER_PERIODIC = 10;  /*!< \brief Maximum number of periodic boundary conditions. */
const unsigned int MAX_STRING_SIZE = 200;     /*!< \brief Maximum number of domains. */
const unsigned int MAX_NUMBER_FFD = 15;	      /*!< \brief Maximum number of FFDBoxes for the FFD. */
const unsigned int MAX_SOLS = 10;		      /*!< \brief Maximum number of solutions at the same time (dimension of solution container array). */
const unsigned int MAX_TERMS = 6;		      /*!< \brief Maximum number of terms in the numerical equations (dimension of solver container array). */
const unsigned int MAX_TERMS_FEA = 10;        /*!< \brief Maximum number of terms in the numerical equations (dimension of solver container array). */
const unsigned int MAX_ZONES = 3;             /*!< \brief Maximum number of zones. */
const unsigned int MAX_FE_KINDS = 4;          /*!< \brief Maximum number of Finite Elements. */
const unsigned int NO_RK_ITER = 0;		      /*!< \brief No Runge-Kutta iteration. */

const unsigned int OVERHEAD = 4;   /*!< \brief Overhead space above nMarker when allocating space for boundary elems (MPI + periodic). */

const unsigned int MESH_0 = 0;  /*!< \brief Definition of the finest grid level. */
const unsigned int MESH_1 = 1;  /*!< \brief Definition of the finest grid level. */
const unsigned int ZONE_0 = 0;  /*!< \brief Definition of the first grid domain. */
const unsigned int ZONE_1 = 1;  /*!< \brief Definition of the second grid domain. */
const unsigned int INST_0 = 0;  /*!< \brief Definition of the first instance per grid level. */

const su2double STANDARD_GRAVITY = 9.80665;           /*!< \brief Acceleration due to gravity at surface of earth. */
const su2double UNIVERSAL_GAS_CONSTANT = 8.3144598;  /*!< \brief Universal gas constant in J/(mol*K) */

const su2double EPS = 1.0E-16;		 /*!< \brief Error scale. */
const su2double TURB_EPS = 1.0E-16;  /*!< \brief Turbulent Error scale. */

const su2double ONE2 = 0.5;			/*!< \brief One divided by two. */
const su2double TWO3 = 2.0 / 3.0;	/*!< \brief Two divided by three. */
const su2double FOUR3 = 4.0 / 3.0;  /*!< \brief Four divided by three. */

const su2double PI_NUMBER = 4.0 * atan(1.0);  /*!< \brief Pi number. */

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

const int CGNS_STRING_SIZE = 33; /*!< \brief Length of strings used in the CGNS format. */
const int SU2_CONN_SIZE   = 10;  /*!< \brief Size of the connectivity array that is allocated for each element
                                             that we read from a mesh file in the format [[globalID vtkType n0 n1 n2 n3 n4 n5 n6 n7 n8]. */
const int SU2_CONN_SKIP   = 2;   /*!< \brief Offset to skip the globalID and VTK type at the start of the element connectivity list for each CGNS element. */

/*!
 * \brief Boolean answers
 */
enum ANSWER {
  NONE = 0,
  NO = 0,   /*!< \brief Boolean definition of no. */
  YES = 1	/*!< \brief Boolean definition of yes. */
};

/*!
 * \brief Average method for marker analyze
 */
enum AVERAGE_TYPE {
  AVERAGE_AREA = 1,     /*!< \brief Area-weighted average. */
  AVERAGE_MASSFLUX = 2  /*!< \brief Mass-flux weighted average. */
};
static const map<string, AVERAGE_TYPE> Average_Map = CCreateMap<string, AVERAGE_TYPE>
("AREA", AVERAGE_AREA)
("MASSFLUX", AVERAGE_MASSFLUX);

/*!
 * \brief different solver types for the CFD component
 */
enum ENUM_SOLVER {
  NO_SOLVER = 0,					/*!< \brief Definition of no solver. */
  EULER = 1,						/*!< \brief Definition of the Euler's solver. */
  NAVIER_STOKES = 2,				/*!< \brief Definition of the Navier-Stokes' solver. */
  RANS = 3,							/*!< \brief Definition of the Reynolds-averaged Navier-Stokes' (RANS) solver. */
  INC_EULER = 4,					/*!< \brief Definition of the incompressible Euler's solver. */
  INC_NAVIER_STOKES =5,				/*!< \brief Definition of the incompressible Navier-Stokes' solver. */
  INC_RANS = 6,						/*!< \brief Definition of the incompressible Reynolds-averaged Navier-Stokes' (RANS) solver. */
  HEAT_EQUATION_FVM = 7,            /*!< \brief Definition of the finite volume heat solver. */
  FLUID_STRUCTURE_INTERACTION = 8,  /*!< \brief Definition of a FSI solver. */
  FEM_ELASTICITY = 9,				/*!< \brief Definition of a FEM solver. */
  ADJ_EULER = 10,					/*!< \brief Definition of the continuous adjoint Euler's solver. */
  ADJ_NAVIER_STOKES = 11,			/*!< \brief Definition of the continuous adjoint Navier-Stokes' solver. */
  ADJ_RANS = 12,					/*!< \brief Definition of the continuous adjoint Reynolds-averaged Navier-Stokes' (RANS) solver. */
  TEMPLATE_SOLVER = 13,             /*!< \brief Definition of template solver. */
  DISC_ADJ_EULER = 15,              /*!< \brief Definition of the discrete adjoint Euler solver. */
  DISC_ADJ_RANS = 16,               /*!< \brief Definition of the discrete adjoint Reynolds-averaged Navier-Stokes' (RANS) solver. */
  DISC_ADJ_NAVIER_STOKES = 17,      /*!< \brief Definition of the discrete adjoint Navier-Stokes' solver. */
  DISC_ADJ_INC_EULER = 18,          /*!< \brief Definition of the discrete adjoint incompressible Euler solver. */
  DISC_ADJ_INC_RANS = 19,           /*!< \brief Definition of the discrete adjoint imcompressible Reynolds-averaged Navier-Stokes' (RANS) solver. */
  DISC_ADJ_INC_NAVIER_STOKES = 20,  /*!< \brief Definition of the doscrete adjoint imcompressible Navier-Stokes'. */
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
/* BEGIN_CONFIG_ENUMS */
static const map<string, ENUM_SOLVER> Solver_Map = CCreateMap<string, ENUM_SOLVER>
("NONE", NO_SOLVER)
("EULER", EULER)
("NAVIER_STOKES", NAVIER_STOKES)
("RANS", RANS)
("INC_EULER", INC_EULER)
("INC_NAVIER_STOKES", INC_NAVIER_STOKES)
("INC_RANS", INC_RANS)
("FEM_EULER", FEM_EULER)
("FEM_NAVIER_STOKES", FEM_NAVIER_STOKES)
("FEM_RANS", FEM_RANS)
("FEM_LES", FEM_LES)
("ADJ_EULER", ADJ_EULER)
("ADJ_NAVIER_STOKES", ADJ_NAVIER_STOKES)
("ADJ_RANS", ADJ_RANS )
("HEAT_EQUATION_FVM", HEAT_EQUATION_FVM)
("ELASTICITY", FEM_ELASTICITY)
("DISC_ADJ_EULER", DISC_ADJ_EULER)
("DISC_ADJ_RANS", DISC_ADJ_RANS)
("DISC_ADJ_NAVIERSTOKES", DISC_ADJ_NAVIER_STOKES)
("DISC_ADJ_INC_EULER", DISC_ADJ_INC_EULER)
("DISC_ADJ_INC_RANS", DISC_ADJ_INC_RANS)
("DISC_ADJ_INC_NAVIERSTOKES", DISC_ADJ_INC_NAVIER_STOKES)
("DISC_ADJ_HEAT_EQUATION_FVM", DISC_ADJ_HEAT)
("DISC_ADJ_FEM_EULER", DISC_ADJ_FEM_EULER)
("DISC_ADJ_FEM_RANS", DISC_ADJ_FEM_RANS)
("DISC_ADJ_FEM_NS", DISC_ADJ_FEM_NS)
("DISC_ADJ_FEM", DISC_ADJ_FEM)
("FLUID_STRUCTURE_INTERACTION", FLUID_STRUCTURE_INTERACTION)
("TEMPLATE_SOLVER", TEMPLATE_SOLVER)
("MULTIPHYSICS", MULTIPHYSICS);

/*!
 * \brief different solver types for the multizone environment component
 */
enum ENUM_MULTIZONE {
  MZ_BLOCK_GAUSS_SEIDEL = 0,   /*!< \brief Definition of a Block-Gauss-Seidel multizone solver. */
  MZ_BLOCK_JACOBI = 1          /*!< \brief Definition of a Block-Jacobi solver. */
};
/* BEGIN_CONFIG_ENUMS */
static const map<string, ENUM_MULTIZONE> Multizone_Map = CCreateMap<string, ENUM_MULTIZONE>
("BLOCK_GAUSS_SEIDEL", MZ_BLOCK_GAUSS_SEIDEL)
("BLOCK_JACOBI", MZ_BLOCK_JACOBI);

/*!
 * \brief types of fluid solvers
 */
enum ENUM_FSI_FLUID_PROBLEM {
      NO_SOLVER_FFSI = 0,      /*!< \brief Definition of no solver. */
      EULER_FFSI = 1,          /*!< \brief Euler equations for the FSI problem */
      NAVIER_STOKES_FFSI = 2,  /*!< \brief NS equations for the FSI problem */
      RANS_FFSI = 3            /*!< \brief RANS equations for the FSI problem */
};
static const map<string, ENUM_FSI_FLUID_PROBLEM> FSI_Fluid_Solver_Map = CCreateMap<string, ENUM_FSI_FLUID_PROBLEM>
("NONE", NO_SOLVER_FFSI)
("EULER", EULER_FFSI)
("NAVIER_STOKES", NAVIER_STOKES_FFSI)
("RANS", RANS_FFSI);

/*!
 * \brief types of structural solvers
 */
enum ENUM_FSI_STRUC_PROBLEM {
  NO_SOLVER_SFSI = 0,		/*!< \brief Definition of no solver. */
  FEM_ELASTICITY_SFSI = 9,	/*!< \brief Nonlinear elasticity equations for the FSI problem */
};
static const map<string, ENUM_FSI_STRUC_PROBLEM> FSI_Struc_Solver_Map = CCreateMap<string, ENUM_FSI_STRUC_PROBLEM>
("NONE", NO_SOLVER_SFSI)
("ELASTICITY", FEM_ELASTICITY_SFSI);

/*!
 * \brief Material geometric conditions
 */
enum ENUM_STRUCT_SOLVER {
    SMALL_DEFORMATIONS = 0,	 /*!< \brief Definition of linear elastic material. */
    LARGE_DEFORMATIONS = 1,	 /*!< \brief Definition of Neo-Hookean material. */
};
static const map<string, ENUM_STRUCT_SOLVER> Struct_Map = CCreateMap<string, ENUM_STRUCT_SOLVER>
("SMALL_DEFORMATIONS", SMALL_DEFORMATIONS)
("LARGE_DEFORMATIONS", LARGE_DEFORMATIONS);

/*!
 * \brief Material model
 */
enum ENUM_MATERIAL_MODEL {
	LINEAR_ELASTIC = 0,			/*!< \brief Definition of linear elastic material. */
	NEO_HOOKEAN = 1,			/*!< \brief Definition of Neo-Hookean material. */
	KNOWLES = 2,				/*!< \brief Definition of Knowles stored-energy potential */
	IDEAL_DE = 3				/*!< \brief Definition of ideal Dielectric Elastomer */
};
static const map<string, ENUM_MATERIAL_MODEL> Material_Map = CCreateMap<string, ENUM_MATERIAL_MODEL>
("LINEAR_ELASTIC", LINEAR_ELASTIC)
("NEO_HOOKEAN", NEO_HOOKEAN)
("KNOWLES", KNOWLES)
("IDEAL_DE", IDEAL_DE);

/*!
 * \brief Material compressibility
 */
enum ENUM_MAT_COMPRESS {
  COMPRESSIBLE_MAT = 0,		      /*!< \brief Definition of compressible material. */
  NEARLY_INCOMPRESSIBLE_MAT = 1,  /*!< \brief Definition of nearly incompressible material. */
};
static const map<string, ENUM_MAT_COMPRESS> MatComp_Map = CCreateMap<string, ENUM_MAT_COMPRESS>
("COMPRESSIBLE", COMPRESSIBLE_MAT)
("NEARLY_INCOMPRESSIBLE", NEARLY_INCOMPRESSIBLE_MAT);

/*!
 * \brief types of interpolators
 */
enum ENUM_INTERPOLATOR {
  NEAREST_NEIGHBOR 	= 0,   	/*!< \brief Nearest Neigbhor interpolation */
  ISOPARAMETRIC 	= 1,	/*!< \brief Isoparametric interpolation, use CONSERVATIVE_INTERPOLATION=YES for conservative interpolation (S.A. Brown 1997).*/
  WEIGHTED_AVERAGE  = 3, 	/*!< \brief Sliding Mesh Approach E. Rinaldi 2015 */
  RADIAL_BASIS_FUNCTION= 4, /*!< \brief Radial basis function interpolation. */
};
static const map<string, ENUM_INTERPOLATOR> Interpolator_Map = CCreateMap<string, ENUM_INTERPOLATOR>
("NEAREST_NEIGHBOR", NEAREST_NEIGHBOR)
("ISOPARAMETRIC",    ISOPARAMETRIC)
("WEIGHTED_AVERAGE", WEIGHTED_AVERAGE)
("RADIAL_BASIS_FUNCTION", RADIAL_BASIS_FUNCTION);

/*!
 * \brief types of radial basis functions
 */
enum ENUM_RADIALBASIS {
  WENDLAND_C2 = 0,        /*!< \brief Wendland C2 radial basis function. */
  INV_MULTI_QUADRIC = 1,  /*!< \brief Inversed multi quartic biharmonic spline. */
  GAUSSIAN = 2,           /*!< \brief Gaussian basis function. */
  THIN_PLATE_SPLINE = 3,  /*!< \brief Thin plate spline. */
  MULTI_QUADRIC = 4,      /*!< \brief Multi quartic biharmonic spline. */
};
static const map<string, ENUM_RADIALBASIS> RadialBasisFunction_Map = CCreateMap<string, ENUM_RADIALBASIS>
("WENDLAND_C2", WENDLAND_C2)
("INV_MULTI_QUADRIC", INV_MULTI_QUADRIC)
("GAUSSIAN", GAUSSIAN)
("THIN_PLATE_SPLINE", THIN_PLATE_SPLINE)
("MULTI_QUADRIC", MULTI_QUADRIC);

/*!
 * \brief types of (coupling) transfers between distinct physical zones
 */
enum ENUM_TRANSFER {
  ZONES_ARE_EQUAL                   = 0,    /*!< \brief Zones are equal - no transfer. */
  NO_COMMON_INTERFACE               = 1,    /*!< \brief No common interface between the zones (geometrical). */
  NO_TRANSFER                       = 2,    /*!< \brief Zones may share a boundary, but still no coupling desired. */
  FLOW_TRACTION                     = 10,   /*!< \brief Flow traction coupling (between fluids and solids). */
  STRUCTURAL_DISPLACEMENTS_LEGACY   = 11,   /*!< \brief Structural displacements (between fluids and solids) - legacy version (to be removed). */
  BOUNDARY_DISPLACEMENTS            = 21,   /*!< \brief Boundary displacements (between fluids and solids) */
  STRUCTURAL_DISPLACEMENTS_DISC_ADJ = 12,   /*!< \brief Adjoints of structural displacements (between fluids and solids). */
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
  DIMENSIONAL = 0,			    /*!< \brief Dimensional simulation (compressible or incompressible). */
  FREESTREAM_PRESS_EQ_ONE = 1,  /*!< \brief Non-dimensional compressible simulation with freestream pressure equal to 1.0. */
  FREESTREAM_VEL_EQ_MACH = 2,   /*!< \brief Non-dimensional compressible simulation with freestream velocity equal to Mach number. */
  FREESTREAM_VEL_EQ_ONE = 3,    /*!< \brief Non-dimensional compressible simulation with freestream pressure equal to 1.0. */
  INITIAL_VALUES   = 4,         /*!< \brief Non-dimensional incompressible simulation based on intial values for external flow. */
  REFERENCE_VALUES = 5          /*!< \brief Non-dimensional incompressible simulation based on custom reference values. */
};
static const map<string, ENUM_KIND_NONDIM> NonDim_Map = CCreateMap<string, ENUM_KIND_NONDIM>
("DIMENSIONAL", DIMENSIONAL)
("FREESTREAM_PRESS_EQ_ONE", FREESTREAM_PRESS_EQ_ONE)
("FREESTREAM_VEL_EQ_MACH",  FREESTREAM_VEL_EQ_MACH)
("FREESTREAM_VEL_EQ_ONE",   FREESTREAM_VEL_EQ_ONE)
("INITIAL_VALUES",   INITIAL_VALUES)
("REFERENCE_VALUES", REFERENCE_VALUES);

/*!
 * \brief different system of measurements
 */
enum ENUM_MEASUREMENTS {
  SI = 0,			/*!< \brief Definition of compressible solver. */
  US = 1			/*!< \brief Definition of incompressible solver. */
};
static const map<string, ENUM_MEASUREMENTS> Measurements_Map = CCreateMap<string, ENUM_MEASUREMENTS>
("SI", SI)
("US", US);

/*!
 * \brief different types of systems
 */
enum RUNTIME_TYPE {
  RUNTIME_FLOW_SYS = 2,			/*!< \brief One-physics case, the code is solving the flow equations(Euler and Navier-Stokes). */
  RUNTIME_TURB_SYS = 3,			/*!< \brief One-physics case, the code is solving the turbulence model. */
  RUNTIME_ADJPOT_SYS = 5,		/*!< \brief One-physics case, the code is solving the adjoint potential flow equation. */
  RUNTIME_ADJFLOW_SYS = 6,		/*!< \brief One-physics case, the code is solving the adjoint equations is being solved (Euler and Navier-Stokes). */
  RUNTIME_ADJTURB_SYS = 7,		/*!< \brief One-physics case, the code is solving the adjoint turbulence model. */
  RUNTIME_MULTIGRID_SYS = 14,   /*!< \brief Full Approximation Storage Multigrid system of equations. */
  RUNTIME_FEA_SYS = 20,		    /*!< \brief One-physics case, the code is solving the FEA equation. */
  RUNTIME_ADJFEA_SYS = 30,		/*!< \brief One-physics case, the code is solving the adjoint FEA equation. */
  RUNTIME_HEAT_SYS = 21,		/*!< \brief One-physics case, the code is solving the heat equation. */
  RUNTIME_ADJHEAT_SYS = 31,     /*!< \brief One-physics case, the code is solving the adjoint heat equation. */
  RUNTIME_TRANS_SYS = 22,		/*!< \brief One-physics case, the code is solving the turbulence model. */
};

const int FLOW_SOL = 0;		/*!< \brief Position of the mean flow solution in the solver container array. */
const int ADJFLOW_SOL = 1;	/*!< \brief Position of the continuous adjoint flow solution in the solver container array. */

const int TURB_SOL = 2;		/*!< \brief Position of the turbulence model solution in the solver container array. */
const int ADJTURB_SOL = 3;	/*!< \brief Position of the continuous adjoint turbulence solution in the solver container array. */

const int TRANS_SOL = 4;	/*!< \brief Position of the transition model solution in the solver container array. */
const int HEAT_SOL = 5;		/*!< \brief Position of the heat equation in the solution solver array. */
const int ADJHEAT_SOL = 6;  /*!< \brief Position of the adjoint heat equation in the solution solver array. */

const int FEA_SOL = 0;		/*!< \brief Position of the FEA equation in the solution solver array. */
const int ADJFEA_SOL = 1;	/*!< \brief Position of the FEA adjoint equation in the solution solver array. */

const int TEMPLATE_SOL = 0;  /*!< \brief Position of the template solution. */

const int CONV_TERM = 0;	       /*!< \brief Position of the convective terms in the numerics container array. */
const int VISC_TERM = 1;           /*!< \brief Position of the viscous terms in the numerics container array. */
const int SOURCE_FIRST_TERM = 2;   /*!< \brief Position of the first source term in the numerics container array. */
const int SOURCE_SECOND_TERM = 3;  /*!< \brief Position of the second source term in the numerics container array. */
const int CONV_BOUND_TERM = 4;     /*!< \brief Position of the convective boundary terms in the numerics container array. */
const int VISC_BOUND_TERM = 5;     /*!< \brief Position of the viscous boundary terms in the numerics container array. */

const int FEA_TERM = 0;			/*!< \brief Position of the finite element analysis terms in the numerics container array. */
const int DE_TERM = 1;			/*!< \brief Position of the dielectric terms in the numerics container array. */

const int MAT_NHCOMP  = 2;   /*!< \brief Position of the Neo-Hookean compressible material model. */
const int MAT_IDEALDE = 3;   /*!< \brief Position of the Ideal-DE material model. */
const int MAT_KNOWLES = 4;   /*!< \brief Position of the Knowles material model. */

const int MESH_SOL = 8;      /*!< \brief Position of the mesh solver. */
const int ADJMESH_SOL = 9;   /*!< \brief Position of the adjoint of the mesh solver. */


/*!
 * \brief types of finite elements (in 2D or 3D)
 */
const int EL_TRIA = 0;		/*!< \brief Elements of three nodes (2D). */
const int EL_QUAD = 1;		/*!< \brief Elements of four nodes (2D). */

const int EL_TETRA = 0;		/*!< \brief Elements of four nodes (3D). */
const int EL_HEXA  = 1;		/*!< \brief Elements of eight nodes (3D). */
const int EL_PYRAM = 2;     /*!< \brief Elements of five nodes (3D). */
const int EL_PRISM = 3;     /*!< \brief Elements of six nodes (3D). */


/*!
 * \brief types of mathematical problem to solve
 */
enum ENUM_MATH_PROBLEM {
  DIRECT = 0,		       /*!< \brief Direct problem */
  CONTINUOUS_ADJOINT = 1,  /*!< \brief Continuous adjoint problem */
  DISCRETE_ADJOINT = 2     /*!< \brief AD-based discrete adjoint problem. */
};
static const map<string, ENUM_MATH_PROBLEM> Math_Problem_Map = CCreateMap<string, ENUM_MATH_PROBLEM>
("DIRECT", DIRECT)
("CONTINUOUS_ADJOINT", CONTINUOUS_ADJOINT)
("DISCRETE_ADJOINT", DISCRETE_ADJOINT);

/*!
 * \brief types of spatial discretizations
 */
enum ENUM_SPACE {
  NO_CONVECTIVE = 0,   /*!< \brief No convective scheme is used. */
  SPACE_CENTERED = 1,  /*!< \brief Space centered convective numerical method. */
  SPACE_UPWIND = 2,	   /*!< \brief Upwind convective numerical method. */
  FINITE_ELEMENT = 3   /*!< \brief Finite element convective numerical method. */
};
static const map<string, ENUM_SPACE> Space_Map = CCreateMap<string, ENUM_SPACE>
("NONE", NO_CONVECTIVE)
("SPACE_CENTERED", SPACE_CENTERED)
("SPACE_UPWIND", SPACE_UPWIND)
("FINITE_ELEMENT", FINITE_ELEMENT);

/*!
 * \brief types of fluid model
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
static const map<string, ENUM_FLUIDMODEL> FluidModel_Map = CCreateMap<string, ENUM_FLUIDMODEL>
("STANDARD_AIR", STANDARD_AIR)
("IDEAL_GAS", IDEAL_GAS)
("VW_GAS", VW_GAS)
("PR_GAS", PR_GAS)
("CONSTANT_DENSITY", CONSTANT_DENSITY)
("INC_IDEAL_GAS", INC_IDEAL_GAS)
("INC_IDEAL_GAS_POLY", INC_IDEAL_GAS_POLY);

/*!
 * \brief types of density models
 */
enum ENUM_DENSITYMODEL {
  CONSTANT = 0,
  BOUSSINESQ = 1,  /*!< \brief BoussinesQ density model. */
  VARIABLE = 2     /*!< \brief Variable density model. */
};
static const map<string, ENUM_DENSITYMODEL> DensityModel_Map = CCreateMap<string, ENUM_DENSITYMODEL>
("CONSTANT", CONSTANT)
("BOUSSINESQ", BOUSSINESQ)
("VARIABLE", VARIABLE);

/*!
 * \brief types of initialization option
 */
enum ENUM_INIT_OPTION {
  REYNOLDS = 0,      /*!< \brief Reynold's number initalization. */
  TD_CONDITIONS = 1  /*!< \brief Total conditions initalization. */
};
static const map<string, ENUM_INIT_OPTION> InitOption_Map = CCreateMap<string, ENUM_INIT_OPTION>
("REYNOLDS", REYNOLDS)
("TD_CONDITIONS", TD_CONDITIONS);

/*!
 * \brief types of initialization option
 */
enum ENUM_FREESTREAM_OPTION {
  TEMPERATURE_FS = 0,  /*!< \brief Temperature initialization. */
  DENSITY_FS = 1       /*!< \brief Density initalization. */
};
static const map<string, ENUM_FREESTREAM_OPTION> FreeStreamOption_Map = CCreateMap<string, ENUM_FREESTREAM_OPTION>
("TEMPERATURE_FS", TEMPERATURE_FS)
("DENSITY_FS", DENSITY_FS);

/*!
 * \brief types of viscosity model
 */
enum ENUM_VISCOSITYMODEL {
  CONSTANT_VISCOSITY = 0,   /*!< \brief Constant viscosity. */
  SUTHERLAND = 1,           /*!< \brief Sutherlands Law viscosity. */
  POLYNOMIAL_VISCOSITY = 2  /*!< \brief Polynomial viscosity. */
};
static const map<string, ENUM_VISCOSITYMODEL> ViscosityModel_Map = CCreateMap<string, ENUM_VISCOSITYMODEL>
("CONSTANT_VISCOSITY", CONSTANT_VISCOSITY)
("SUTHERLAND", SUTHERLAND)
("POLYNOMIAL_VISCOSITY", POLYNOMIAL_VISCOSITY);

/*!
 * \brief types of thermal conductivity model
 */
enum ENUM_CONDUCTIVITYMODEL {
  CONSTANT_CONDUCTIVITY = 0,   /*!< \brief Constant thermal conductivity. */
  CONSTANT_PRANDTL = 1,        /*!< \brief Constant Prandtl number. */
  POLYNOMIAL_CONDUCTIVITY = 2  /*!< \brief Polynomial thermal conductivity. */
};
static const map<string, ENUM_CONDUCTIVITYMODEL> ConductivityModel_Map = CCreateMap<string, ENUM_CONDUCTIVITYMODEL>
("CONSTANT_CONDUCTIVITY", CONSTANT_CONDUCTIVITY)
("CONSTANT_PRANDTL", CONSTANT_PRANDTL)
("POLYNOMIAL_CONDUCTIVITY", POLYNOMIAL_CONDUCTIVITY);

/*!
 * \brief types of turbulent thermal conductivity model
 */
enum ENUM_CONDUCTIVITYMODEL_TURB {
  NO_CONDUCTIVITY_TURB  = 0,  /*!< \brief No turbulent contribution to the effective thermal conductivity for RANS. */
  CONSTANT_PRANDTL_TURB = 1   /*!< \brief Include contribution to effective conductivity using constant turbulent Prandtl number for RANS. */
};
static const map<string, ENUM_CONDUCTIVITYMODEL_TURB> TurbConductivityModel_Map = CCreateMap<string, ENUM_CONDUCTIVITYMODEL_TURB>
("NONE", NO_CONDUCTIVITY_TURB)
("CONSTANT_PRANDTL_TURB", CONSTANT_PRANDTL_TURB);

/*!
 * \brief types of unsteady mesh motion
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
static const map<string, ENUM_GRIDMOVEMENT> GridMovement_Map = CCreateMap<string, ENUM_GRIDMOVEMENT>
("NONE", NO_MOVEMENT)
("RIGID_MOTION", RIGID_MOTION)
("ROTATING_FRAME", ROTATING_FRAME)
("ELASTICITY", ELASTICITY)
("MOVING_HTP", MOVING_HTP)
("STEADY_TRANSLATION", STEADY_TRANSLATION)
("GUST", GUST);

enum ENUM_SURFACEMOVEMENT {
  DEFORMING = 1,                 /*!< \brief Simulation with deformation. */
  MOVING_WALL = 2,               /*!< \brief Simulation with moving wall. */
  AEROELASTIC = 3,               /*!< \brief Simulation with aeroelastic motion. */
  AEROELASTIC_RIGID_MOTION = 4,  /*!< \brief Simulation with rotation and aeroelastic motion. */
  FLUID_STRUCTURE = 5,		     /*!< \brief Fluid structure deformation. */
  EXTERNAL = 6,                  /*!< \brief Simulation with external motion. */
  EXTERNAL_ROTATION = 7,         /*!< \brief Simulation with external rotation motion. */
  FLUID_STRUCTURE_STATIC = 8     /*!< \brief Fluid structure deformation with no grid velocity. */
};
static const map<string, ENUM_SURFACEMOVEMENT> SurfaceMovement_Map = CCreateMap<string, ENUM_SURFACEMOVEMENT>
("DEFORMING", DEFORMING)
("MOVING_WALL", MOVING_WALL)
("AEROELASTIC_RIGID_MOTION", AEROELASTIC_RIGID_MOTION)
("AEROELASTIC", AEROELASTIC)
("FLUID_STRUCTURE_STATIC", FLUID_STRUCTURE_STATIC)
("FLUID_STRUCTURE", FLUID_STRUCTURE)
("EXTERNAL", EXTERNAL)
("EXTERNAL_ROTATION", EXTERNAL_ROTATION);

/*!
 * \brief type of wind gusts
 */
enum ENUM_GUST_TYPE {
  NO_GUST = 0,      /*!< \brief No gust. */
  TOP_HAT = 1,      /*!< \brief Top-hat function shaped gust  */
  SINE = 2,         /*!< \brief Sine shaped gust */
  ONE_M_COSINE = 3, /*!< \brief 1-cosine shaped gust */
  VORTEX = 4,       /*!< \brief A gust made from vortices */
  EOG = 5           /*!< \brief An extreme operating gust */
};
static const map<string, ENUM_GUST_TYPE> Gust_Type_Map = CCreateMap<string, ENUM_GUST_TYPE>
("NONE", NO_GUST)
("TOP_HAT", TOP_HAT)
("SINE", SINE)
("ONE_M_COSINE", ONE_M_COSINE)
("VORTEX", VORTEX)
("EOG", EOG);

/*!
 * \brief type of wind direction
 */
enum ENUM_GUST_DIR {
  X_DIR = 0,  /*!< \brief Gust direction-X. */
  Y_DIR = 1   /*!< \brief Gust direction-Y. */
};
static const map<string, ENUM_GUST_DIR> Gust_Dir_Map = CCreateMap<string, ENUM_GUST_DIR>
("X_DIR", X_DIR)
("Y_DIR", Y_DIR);

// If you add to ENUM_CENTERED, you must also add the option to ENUM_CONVECTIVE
/*!
 * \brief types of centered spatial discretizations
 */
enum ENUM_CENTERED {
  NO_CENTERED = 0,    /*!< \brief No centered scheme is used. */
  JST = 1,            /*!< \brief Jameson-Smith-Turkel centered numerical method. */
  LAX = 2,            /*!< \brief Lax-Friedrich centered numerical method. */
  JST_KE = 4          /*!< \brief Kinetic Energy preserving Jameson-Smith-Turkel centered numerical method. */
};
static const map<string, ENUM_CENTERED> Centered_Map = CCreateMap<string, ENUM_CENTERED>
("NONE", NO_CENTERED)
("JST", JST)
("JST_KE", JST_KE)
("LAX-FRIEDRICH", LAX);


// If you add to ENUM_UPWIND, you must also add the option to ENUM_CONVECTIVE
/*!
 * \brief types of upwind spatial discretizations
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
static const map<string, ENUM_UPWIND> Upwind_Map = CCreateMap<string, ENUM_UPWIND>
("NONE", NO_UPWIND)
("ROE", ROE)
("TURKEL_PREC", TURKEL)
("AUSM", AUSM)
("AUSMPLUSUP", AUSMPLUSUP)
("AUSMPLUSUP2", AUSMPLUSUP2)
("SLAU", SLAU)
("HLLC", HLLC)
("SW", SW)
("MSW", MSW)
("CUSP", CUSP)
("SCALAR_UPWIND", SCALAR_UPWIND)
("CONVECTIVE_TEMPLATE", CONVECTIVE_TEMPLATE)
("L2ROE", L2ROE)
("LMROE", LMROE)
("SLAU2", SLAU2)
("FDS", FDS)
("LAX-FRIEDRICH", LAX_FRIEDRICH);

/*!
 * \brief types of FEM spatial discretizations
 */
enum ENUM_FEM {
  NO_FEM = 0,  /*!< \brief No finite element scheme is used. */
  DG = 1       /*!< \brief Discontinuous Galerkin numerical method. */
};
static const map<string, ENUM_FEM> FEM_Map = CCreateMap<string, ENUM_FEM>
("NONE", NO_FEM)
("DG", DG);

/*!
 * \brief types of shock capturing method in Discontinuous Galerkin numerical method.
 */
enum ENUM_SHOCK_CAPTURING_DG {
  NO_SHOCK_CAPTURING = 0,            /*!< \brief Shock capturing is not used. */
  PERSSON = 1                        /*!< \brief Per-Olof Persson's sub-cell shock capturing method. */
};
static const map<string, ENUM_SHOCK_CAPTURING_DG> ShockCapturingDG_Map = CCreateMap<string, ENUM_SHOCK_CAPTURING_DG>
("NONE", NO_SHOCK_CAPTURING)
("PERSSON", PERSSON);

/*!
 * \brief types of matrix coloring to compute a sparse Jacobian matrix.
 */
enum ENUM_MATRIX_COLORING {
  GREEDY_COLORING = 0,            /*!< \brief Greedy type of algorithm for the coloring. */
  NATURAL_COLORING = 1            /*!< \brief One color for every DOF, very slow. Only to be used for debugging. */
};
static const map<string, ENUM_MATRIX_COLORING> MatrixColoring_Map = CCreateMap<string, ENUM_MATRIX_COLORING>
("GREEDY_COLORING", GREEDY_COLORING)
("NATURAL_COLORING", NATURAL_COLORING);

/*!
 * \brief types of slope limiters
 */
enum ENUM_LIMITER {
  NO_LIMITER           = 0, /*!< \brief No limiter. */
  VENKATAKRISHNAN      = 1,	/*!< \brief Slope limiter using Venkatakrisnan method (stencil formulation). */
  VENKATAKRISHNAN_WANG = 2,	/*!< \brief Slope limiter using Venkatakrisnan method, eps based on solution (stencil formulation). */
  BARTH_JESPERSEN      = 3, /*!< \brief Slope limiter using Barth-Jespersen method (stencil formulation). */
  VAN_ALBADA_EDGE      = 4, /*!< \brief Slope limiter using Van Albada method (edge formulation). */
  SHARP_EDGES          = 5, /*!< \brief Slope limiter using sharp edges. */
  WALL_DISTANCE        = 6  /*!< \brief Slope limiter using wall distance. */
};
static const map<string, ENUM_LIMITER> Limiter_Map = CCreateMap<string, ENUM_LIMITER>
("NONE", NO_LIMITER)
("VENKATAKRISHNAN", VENKATAKRISHNAN)
("VENKATAKRISHNAN_WANG", VENKATAKRISHNAN_WANG)
("BARTH_JESPERSEN", BARTH_JESPERSEN)
("VAN_ALBADA_EDGE", VAN_ALBADA_EDGE)
("SHARP_EDGES", SHARP_EDGES)
("WALL_DISTANCE", WALL_DISTANCE);

/*!
 * \brief types of turbulent models
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
static const map<string, ENUM_TURB_MODEL> Turb_Model_Map = CCreateMap<string, ENUM_TURB_MODEL>
("NONE", NO_TURB_MODEL)
("SA", SA)
("SA_NEG", SA_NEG)
("SA_E", SA_E)
("SA_COMP", SA_COMP)
("SA_E_COMP", SA_E_COMP)
("SST", SST)
("SST_SUST", SST_SUST);

/*!
 * \brief types of transition models
 */
enum ENUM_TRANS_MODEL {
  NO_TRANS_MODEL = 0,  /*!< \brief No transition model. */
  LM = 1,	           /*!< \brief Kind of transition model (Langtry-Menter (LM) for SST and Spalart-Allmaras). */
  BC = 2	           /*!< \brief Kind of transition model (BAS-CAKMAKCIOGLU (BC) for Spalart-Allmaras). */
};
static const map<string, ENUM_TRANS_MODEL> Trans_Model_Map = CCreateMap<string, ENUM_TRANS_MODEL>
("NONE", NO_TRANS_MODEL)
("LM", LM)
("BC", BC); //BAS-CAKMAKCIOGLU

/*!
 * \brief types of subgrid scale models
 */
enum ENUM_SGS_MODEL {
  NO_SGS_MODEL = 0, /*!< \brief No subgrid scale model. */
  IMPLICIT_LES = 1, /*!< \brief Implicit LES, i.e. no explicit SGS model. */
  SMAGORINSKY  = 2, /*!< \brief Smagorinsky SGS model. */
  WALE         = 3, /*!< \brief Wall-Adapting Local Eddy-viscosity SGS model. */
  VREMAN       = 4  /*!< \brief Vreman SGS model. */
};
static const map<string, ENUM_SGS_MODEL> SGS_Model_Map = CCreateMap<string, ENUM_SGS_MODEL>
("NONE",         NO_SGS_MODEL)
("IMPLICIT_LES", IMPLICIT_LES)
("SMAGORINSKY",  SMAGORINSKY)
("WALE",         WALE)
("VREMAN",       VREMAN);


/*!
 * \brief types of window (weight) functions for cost functional
 */
enum WINDOW_FUNCTION {
  SQUARE = 0,          /*!< \brief No weight function  (order 1)*/
  HANN   = 1,           /*!< \brief Hann-type weight function (order 3) */
  HANN_SQUARE  = 2,    /*!< \brief Hann-squared type weight function (order 5)*/
  BUMP  = 3,            /*!< \brief bump type weight function (exponential order of convergence) */
};

static const map<string, WINDOW_FUNCTION> Window_Map = CCreateMap<string, WINDOW_FUNCTION>
("SQUARE", SQUARE)
("HANN", HANN)
("HANN_SQUARE", HANN_SQUARE)
("BUMP", BUMP);

/*!
 * \brief types of hybrid RANS/LES models
 */
enum ENUM_HYBRIDRANSLES {
  NO_HYBRIDRANSLES = 0,  /*!< \brief No turbulence model. */
  SA_DES   = 1,          /*!< \brief Kind of Hybrid RANS/LES (SA - Detached Eddy Simulation (DES)). */
  SA_DDES  = 2,          /*!< \brief Kind of Hybrid RANS/LES (SA - Delayed DES (DDES) with Delta_max SGS ). */
  SA_ZDES  = 3,          /*!< \brief Kind of Hybrid RANS/LES (SA - Delayed DES (DDES) with Vorticity based SGS like Zonal DES). */
  SA_EDDES = 4           /*!< \brief Kind of Hybrid RANS/LES (SA - Delayed DES (DDES) with Shear Layer Adapted SGS: Enhanced DDES). */
};
static const map<string, ENUM_HYBRIDRANSLES> HybridRANSLES_Map = CCreateMap<string, ENUM_HYBRIDRANSLES>
("NONE", NO_HYBRIDRANSLES)
("SA_DES", SA_DES)
("SA_DDES", SA_DDES)
("SA_ZDES", SA_ZDES)
("SA_EDDES", SA_EDDES);

/*!
 * \brief types of Roe Low Dissipation Schemes
 */
enum ENUM_ROELOWDISS {
    NO_ROELOWDISS = 0, /*!< \brief No Roe Low Dissipation model. */
    FD            = 1, /*!< \brief Numerical Blending based on DDES's F_d function */
    NTS           = 2, /*!< \brief Numerical Blending of Travin and Shur. */
    NTS_DUCROS    = 3, /*!< \brief Numerical Blending of Travin and Shur + Ducros' Shock Sensor. */
    FD_DUCROS     = 4  /*!< \brief Numerical Blending based on DDES's F_d function + Ducros' Shock Sensor */
};
static const map<string, ENUM_ROELOWDISS> RoeLowDiss_Map = CCreateMap<string, ENUM_ROELOWDISS>
("NONE", NO_ROELOWDISS)
("FD", FD)
("NTS", NTS)
("NTS_DUCROS", NTS_DUCROS)
("FD_DUCROS", FD_DUCROS);

/*!
 * \brief types of wall functions.
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
static const map<string, ENUM_WALL_FUNCTIONS> Wall_Functions_Map = CCreateMap<string, ENUM_WALL_FUNCTIONS>
("NO_WALL_FUNCTION",          NO_WALL_FUNCTION)
("STANDARD_WALL_FUNCTION",    STANDARD_WALL_FUNCTION)
("ADAPTIVE_WALL_FUNCTION",    ADAPTIVE_WALL_FUNCTION)
("SCALABLE_WALL_FUNCTION",    SCALABLE_WALL_FUNCTION)
("EQUILIBRIUM_WALL_MODEL",    EQUILIBRIUM_WALL_MODEL)
("NONEQUILIBRIUM_WALL_MODEL", NONEQUILIBRIUM_WALL_MODEL)
("LOGARITHMIC_WALL_MODEL", LOGARITHMIC_WALL_MODEL);

/*!
 * \brief type of time integration schemes
 */
enum ENUM_TIME_INT {
  RUNGE_KUTTA_EXPLICIT = 1,	   /*!< \brief Explicit Runge-Kutta time integration definition. */
  EULER_EXPLICIT = 2,   	   /*!< \brief Explicit Euler time integration definition. */
  EULER_IMPLICIT = 3,   	   /*!< \brief Implicit Euler time integration definition. */
  CLASSICAL_RK4_EXPLICIT = 4,  /*!< \brief Classical RK4 time integration definition. */
  ADER_DG = 5                  /*!< \brief ADER-DG time integration definition. */
};
static const map<string, ENUM_TIME_INT> Time_Int_Map = CCreateMap<string, ENUM_TIME_INT>
("RUNGE-KUTTA_EXPLICIT", RUNGE_KUTTA_EXPLICIT)
("EULER_EXPLICIT", EULER_EXPLICIT)
("EULER_IMPLICIT", EULER_IMPLICIT)
("CLASSICAL_RK4_EXPLICIT", CLASSICAL_RK4_EXPLICIT)
("ADER_DG", ADER_DG);

/*!
 * \brief type of predictor for the ADER-DG time integration scheme.
 */
enum ENUM_ADER_PREDICTOR {
  ADER_ALIASED_PREDICTOR     = 1, /*!< \brief Aliased predictor, easiest to do. */
  ADER_NON_ALIASED_PREDICTOR = 2  /*!< \brief Non-aliased predictor. Consistent, but more difficult. */
};
static const map<string, ENUM_ADER_PREDICTOR> Ader_Predictor_Map = CCreateMap<string, ENUM_ADER_PREDICTOR>
("ADER_ALIASED_PREDICTOR", ADER_ALIASED_PREDICTOR)
("ADER_NON_ALIASED_PREDICTOR", ADER_NON_ALIASED_PREDICTOR);

/*!
 * \brief type of heat timestep calculation
 */
enum ENUM_HEAT_TIMESTEP {
  MINIMUM = 1,     /*!< \brief Local time stepping based on minimum lambda.*/
  CONVECTIVE = 2,  /*!< \brief Local time stepping based on convective spectral radius.*/
  VISCOUS = 3,     /*!< \brief Local time stepping based on viscous spectral radius.*/
  BYFLOW = 4,      /*!< \brief Unsing the mean solvers time step. */
};
static const map<string, ENUM_HEAT_TIMESTEP> Heat_TimeStep_Map = CCreateMap<string, ENUM_HEAT_TIMESTEP>
("LOCAL", MINIMUM)
("CONVECTIVE", CONVECTIVE)
("VISCOUS", VISCOUS)
("BYFLOW", BYFLOW);

/*!
 * \brief type of time integration schemes
 */
enum ENUM_TIME_INT_FEA {
  CD_EXPLICIT = 1,       /*!< \brief Support for implementing an explicit method. */
  NEWMARK_IMPLICIT = 2,  /*!< \brief Implicit Newmark integration definition. */
  GENERALIZED_ALPHA = 3  /*!< \brief Support for implementing another implicit method. */
};
static const map<string, ENUM_TIME_INT_FEA> Time_Int_Map_FEA = CCreateMap<string, ENUM_TIME_INT_FEA>
("CD_EXPLICIT", CD_EXPLICIT)
("NEWMARK_IMPLICIT", NEWMARK_IMPLICIT)
("GENERALIZED_ALPHA", GENERALIZED_ALPHA);

/*!
 * \brief type of time integration schemes
 */
enum ENUM_SPACE_ITE_FEA {
  NEWTON_RAPHSON = 1,			/*!< \brief Full Newton-Rapshon method. */
  MODIFIED_NEWTON_RAPHSON = 2   /*!< \brief Modified Newton-Raphson method. */
};
static const map<string, ENUM_SPACE_ITE_FEA> Space_Ite_Map_FEA = CCreateMap<string, ENUM_SPACE_ITE_FEA>
("NEWTON_RAPHSON", NEWTON_RAPHSON)
("MODIFIED_NEWTON_RAPHSON", MODIFIED_NEWTON_RAPHSON);

/*!
 * \brief types of schemes to compute the flow gradient
 */
enum ENUM_FLOW_GRADIENT {
  NO_GRADIENT            = 0, /*!< \brief No gradient method. Only possible for reconstruction gradient, in which case, the option chosen for NUM_METHOD_GRAD is used. */
  GREEN_GAUSS            = 1,	/*!< \brief Gradient computation using Green-Gauss theorem. */
  LEAST_SQUARES          = 2, /*!< \brief Gradient computation using unweighted least squares. */
  WEIGHTED_LEAST_SQUARES = 3	/*!< \brief Gradients computation using inverse-distance weighted least squares. */
};
static const map<string, ENUM_FLOW_GRADIENT> Gradient_Map = CCreateMap<string, ENUM_FLOW_GRADIENT>
("NONE", NO_GRADIENT)
("GREEN_GAUSS", GREEN_GAUSS)
("LEAST_SQUARES", LEAST_SQUARES)
("WEIGHTED_LEAST_SQUARES", WEIGHTED_LEAST_SQUARES);

/*!
 * \brief types of action to take on a geometry structure
 */
enum GEOMETRY_ACTION {
  ALLOCATE = 0,     /*!<  \brief Allocate geometry structure. */
  UPDATE = 1        /*!<  \brief Update geometry structure (grid moving, adaptation, etc.). */
};

/*!
 * \brief types of action to perform when doing the geometry evaluation
 */
enum GEOMETRY_MODE {
  FUNCTION = 0,     /*!<  \brief Geometrical analysis. */
  GRADIENT = 1      /*!<  \brief Geometrical analysis and gradient using finite differences. */
};
static const map<string, GEOMETRY_MODE> GeometryMode_Map = CCreateMap<string, GEOMETRY_MODE>
("FUNCTION", FUNCTION)
("GRADIENT", GRADIENT);

/*!
 * \brief types of boundary conditions
 */
enum BC_TYPE {
  EULER_WALL = 1,		       /*!< \brief Boundary Euler wall definition. */
  FAR_FIELD = 2,		       /*!< \brief Boundary far-field definition. */
  SYMMETRY_PLANE = 3,   	   /*!< \brief Boundary symmetry plane definition. */
  INLET_FLOW = 4,		       /*!< \brief Boundary inlet flow definition. */
  OUTLET_FLOW = 5,		       /*!< \brief Boundary outlet flow definition. */
  PERIODIC_BOUNDARY = 6,	   /*!< \brief Periodic boundary definition. */
  NEARFIELD_BOUNDARY = 7,	   /*!< \brief Near-Field boundary definition. */
  ELECTRODE_BOUNDARY = 8,	   /*!< \brief Electrode boundary definition. */
  DIELEC_BOUNDARY = 9,	       /*!< \brief Dipoisson boundary definition. */
  CUSTOM_BOUNDARY = 10,        /*!< \brief custom boundary definition. */
  INTERFACE_BOUNDARY = 11,	   /*!< \brief Domain interface boundary definition. */
  DIRICHLET = 12,		       /*!< \brief Boundary Euler wall definition. */
  NEUMANN = 13,		           /*!< \brief Boundary Neumann definition. */
  DISPLACEMENT_BOUNDARY = 14,  /*!< \brief Boundary displacement definition. */
  LOAD_BOUNDARY = 15,		   /*!< \brief Boundary Load definition. */
  FLOWLOAD_BOUNDARY = 16,	   /*!< \brief Boundary Load definition. */
  SUPERSONIC_INLET = 19,	   /*!< \brief Boundary supersonic inlet definition. */
  SUPERSONIC_OUTLET = 20,	   /*!< \brief Boundary supersonic inlet definition. */
  ENGINE_INFLOW = 21,		   /*!< \brief Boundary nacelle inflow. */
  ENGINE_EXHAUST = 22,		   /*!< \brief Boundary nacelle exhaust. */
  RIEMANN_BOUNDARY= 24,        /*!< \brief Riemann Boundary definition. */
  ISOTHERMAL = 25,             /*!< \brief No slip isothermal wall boundary condition. */
  HEAT_FLUX  = 26,             /*!< \brief No slip constant heat flux wall boundary condition. */
  ACTDISK_INLET = 32,	       /*!< \brief Actuator disk inlet boundary definition. */
  ACTDISK_OUTLET = 33,	       /*!< \brief Actuator disk outlet boundary definition. */
  CLAMPED_BOUNDARY = 34,	   /*!< \brief Clamped Boundary definition. */
  LOAD_DIR_BOUNDARY = 35,	   /*!< \brief Boundary Load definition. */
  LOAD_SINE_BOUNDARY = 36,	   /*!< \brief Sine-waveBoundary Load definition. */
  GILES_BOUNDARY= 37,          /*!< \brief Giles Boundary definition. */
  INTERNAL_BOUNDARY= 38,       /*!< \brief Internal Boundary definition. */
  FLUID_INTERFACE = 39,	       /*!< \brief Domain interface definition. */
  DISP_DIR_BOUNDARY = 40,      /*!< \brief Boundary displacement definition. */
  DAMPER_BOUNDARY = 41,        /*!< \brief Damper. */
  CHT_WALL_INTERFACE = 50,     /*!< \brief Domain interface definition. */
  SEND_RECEIVE = 99,		   /*!< \brief Boundary send-receive definition. */
};


/*!
 * \brief different regime modes
 */
enum ENUM_2DFORM {
  PLANE_STRESS = 0,			/*!< \brief Definition of plane stress solver. */
  PLANE_STRAIN = 1			/*!< \brief Definition of plane strain solver. */
};
static const map<string, ENUM_2DFORM> ElasForm_2D = CCreateMap<string, ENUM_2DFORM>
("PLANE_STRESS", PLANE_STRESS)
("PLANE_STRAIN", PLANE_STRAIN);

/*!
 * \brief Kinds of relaxation for FSI problem
 */
enum ENUM_AITKEN {
  NO_RELAXATION = 0,			/*!< \brief No relaxation in the strongly coupled approach. */
  FIXED_PARAMETER = 1,			/*!< \brief Relaxation with a fixed parameter. */
  AITKEN_DYNAMIC = 2			/*!< \brief Relaxation using Aitken's dynamic parameter. */
};
static const map<string, ENUM_AITKEN> AitkenForm_Map = CCreateMap<string, ENUM_AITKEN>
("NONE", NO_RELAXATION)
("FIXED_PARAMETER", FIXED_PARAMETER)
("AITKEN_DYNAMIC", AITKEN_DYNAMIC);

/*!
 * \brief types of dynamic transfer methods
 */
enum ENUM_DYN_TRANSFER_METHOD {
  INSTANTANEOUS = 1,   /*!< \brief No ramp, load is transfer instantaneously. */
  POL_ORDER_1 = 2,     /*!< \brief The load is transferred using a ramp. */
  POL_ORDER_3 = 3,     /*!< \brief The load is transferred using an order 3 polynomial function */
  POL_ORDER_5 = 4,     /*!< \brief The load is transferred using an order 5 polynomial function */
  SIGMOID_10 = 5,      /*!< \brief The load is transferred using a sigmoid with parameter 10 */
  SIGMOID_20 = 6       /*!< \brief The load is transferred using a sigmoid with parameter 20 */
};
static const map<string, ENUM_DYN_TRANSFER_METHOD> Dyn_Transfer_Method_Map = CCreateMap<string, ENUM_DYN_TRANSFER_METHOD>
("INSTANTANEOUS", INSTANTANEOUS)
("RAMP", POL_ORDER_1)
("CUBIC", POL_ORDER_3)
("QUINTIC", POL_ORDER_5)
("SIGMOID_10", SIGMOID_10)
("SIGMOID_20", SIGMOID_20);

/*!
 * \brief Kinds of Design Variables for FEA problems 
 */
enum ENUM_DVFEA {
  NODV_FEA = 0,         /*!< \brief No design variable for FEA problems. */
  YOUNG_MODULUS = 1,	/*!< \brief Young modulus (E) as design variable. */
  POISSON_RATIO = 2,  	/*!< \brief Poisson ratio (Nu) as design variable. */
  DENSITY_VAL = 3,      /*!< \brief Density (Rho) as design variable. */
  DEAD_WEIGHT = 4,      /*!< \brief Dead Weight (Rho_DL) as design variable. */
  ELECTRIC_FIELD = 5    /*!< \brief Electric field (E) as design variable. */
};
static const map<string, ENUM_DVFEA> DVFEA_Map = CCreateMap<string, ENUM_DVFEA>
("NONE", NODV_FEA)
("YOUNG_MODULUS", YOUNG_MODULUS)
("POISSON_RATIO", POISSON_RATIO)
("DENSITY", DENSITY_VAL)
("DEAD_WEIGHT", DEAD_WEIGHT)
("ELECTRIC_FIELD", ELECTRIC_FIELD);

/*!
 * \brief types Riemann boundary treatments
 */
enum RIEMANN_TYPE {
  TOTAL_CONDITIONS_PT = 1,		    /*!< \brief User specifies total pressure, total temperature, and flow direction. */
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

static const map<string, RIEMANN_TYPE> Riemann_Map = CCreateMap<string, RIEMANN_TYPE>
("TOTAL_CONDITIONS_PT", TOTAL_CONDITIONS_PT)
("DENSITY_VELOCITY", DENSITY_VELOCITY)
("STATIC_PRESSURE", STATIC_PRESSURE)
("TOTAL_SUPERSONIC_INFLOW", TOTAL_SUPERSONIC_INFLOW)
("STATIC_SUPERSONIC_INFLOW_PT", STATIC_SUPERSONIC_INFLOW_PT)
("STATIC_SUPERSONIC_INFLOW_PD", STATIC_SUPERSONIC_INFLOW_PD)
("MIXING_IN", MIXING_IN)
("MIXING_OUT", MIXING_OUT)
("MIXING_IN_1D", MIXING_IN_1D)
("MIXING_OUT_1D", MIXING_OUT_1D)
("SUPERSONIC_OUTFLOW", SUPERSONIC_OUTFLOW)
("RADIAL_EQUILIBRIUM", RADIAL_EQUILIBRIUM)
("TOTAL_CONDITIONS_PT_1D", TOTAL_CONDITIONS_PT_1D)
("STATIC_PRESSURE_1D", STATIC_PRESSURE_1D);


static const map<string, RIEMANN_TYPE> Giles_Map = CCreateMap<string, RIEMANN_TYPE>
("TOTAL_CONDITIONS_PT", TOTAL_CONDITIONS_PT)
("DENSITY_VELOCITY", DENSITY_VELOCITY)
("STATIC_PRESSURE", STATIC_PRESSURE)
("TOTAL_SUPERSONIC_INFLOW", TOTAL_SUPERSONIC_INFLOW)
("STATIC_SUPERSONIC_INFLOW_PT", STATIC_SUPERSONIC_INFLOW_PT)
("STATIC_SUPERSONIC_INFLOW_PD", STATIC_SUPERSONIC_INFLOW_PD)
("MIXING_IN", MIXING_IN)
("MIXING_OUT", MIXING_OUT)
("MIXING_IN_1D", MIXING_IN_1D)
("MIXING_OUT_1D", MIXING_OUT_1D)
("SUPERSONIC_OUTFLOW", SUPERSONIC_OUTFLOW)
("RADIAL_EQUILIBRIUM", RADIAL_EQUILIBRIUM)
("TOTAL_CONDITIONS_PT_1D", TOTAL_CONDITIONS_PT_1D)
("STATIC_PRESSURE_1D", STATIC_PRESSURE_1D);

/*!
 * \brief types of mixing process for averaging quantities at the boundaries.
 */
enum AVERAGEPROCESS_TYPE {
  ALGEBRAIC = 1,  /*!< \brief an algebraic average is computed at the boundary of interest. */
  AREA = 2,       /*!< \brief an area average is computed at the boundary of interest. */
  MIXEDOUT = 3,	  /*!< \brief an mixed-out average is computed at the boundary of interest. */
  MASSFLUX = 4    /*!< \brief a mass flow average is computed at the boundary of interest. */
};
static const map<string, AVERAGEPROCESS_TYPE> AverageProcess_Map = CCreateMap<string, AVERAGEPROCESS_TYPE>
("ALGEBRAIC", ALGEBRAIC)
("AREA", AREA)
("MIXEDOUT",  MIXEDOUT)
("MASSFLUX", MASSFLUX);

/*!
 * \brief types of mixing process for averaging quantities at the boundaries.
 */
enum MIXINGPLANE_INTERFACE_TYPE {
  MATCHING = 1,		        /*!< \brief an algebraic average is computed at the boundary of interest. */
  NEAREST_SPAN = 2,         /*!< \brief an area average is computed at the boundary of interest. */
  LINEAR_INTERPOLATION = 3	/*!< \brief an mixed-out average is computed at the boundary of interest. */
};
static const map<string, MIXINGPLANE_INTERFACE_TYPE> MixingPlaneInterface_Map = CCreateMap<string, MIXINGPLANE_INTERFACE_TYPE>
("MATCHING", MATCHING)
("NEAREST_SPAN",  NEAREST_SPAN)
("LINEAR_INTERPOLATION", LINEAR_INTERPOLATION);

/*!
 * \brief this option allow to compute the span-wise section in different ways.
 */
enum SPANWISE_TYPE {
  AUTOMATIC = 1,		/*!< \brief number of span-wise section are computed automatically */
  EQUISPACED = 2        /*!< \brief number of span-wise section are specified from the user */
};
static const map<string, SPANWISE_TYPE> SpanWise_Map = CCreateMap<string, SPANWISE_TYPE>
("AUTOMATIC", AUTOMATIC)
("EQUISPACED", EQUISPACED);

/*!
 * \brief types of mixing process for averaging quantities at the boundaries.
 */
enum TURBOMACHINERY_TYPE {
  AXIAL       = 1,		  /*!< \brief axial turbomachinery. */
  CENTRIFUGAL = 2,        /*!< \brief centrifugal turbomachinery. */
  CENTRIPETAL = 3,		  /*!< \brief centripetal turbomachinery. */
  CENTRIPETAL_AXIAL = 4,  /*!< \brief mixed flow turbine. */
  AXIAL_CENTRIFUGAL = 5	  /*!< \brief mixed flow turbine. */
};
static const map<string, TURBOMACHINERY_TYPE> TurboMachinery_Map = CCreateMap<string, TURBOMACHINERY_TYPE>
("AXIAL", AXIAL)
("CENTRIFUGAL", CENTRIFUGAL)
("CENTRIPETAL",  CENTRIPETAL)
("CENTRIPETAL_AXIAL",  CENTRIPETAL_AXIAL)
("AXIAL_CENTRIFUGAL",  AXIAL_CENTRIFUGAL);

/*!
 * \brief types of Turbomachinery performance flag.
 */
enum TURBO_MARKER_TYPE{
  INFLOW   = 1,	  /*!< \brief flag for inflow marker for compute turboperformance. */
  OUTFLOW = 2     /*!< \brief flag for outflow marker for compute turboperformance. */
};

/*!
 * \brief types inlet boundary treatments
 */
enum INLET_TYPE {
  TOTAL_CONDITIONS = 1,	  /*!< \brief User specifies total pressure, total temperature, and flow direction. */
  MASS_FLOW = 2,          /*!< \brief User specifies density and velocity (mass flow). */
  INPUT_FILE = 3,         /*!< \brief User specifies an input file. */
  VELOCITY_INLET = 4,     /*!< \brief Velocity inlet for an incompressible flow. */
  PRESSURE_INLET = 5      /*!< \brief Total pressure inlet for an incompressible flow. */
};
static const map<string, INLET_TYPE> Inlet_Map = CCreateMap<string, INLET_TYPE>
("TOTAL_CONDITIONS", TOTAL_CONDITIONS)
("MASS_FLOW", MASS_FLOW)
("INPUT_FILE", INPUT_FILE)
("VELOCITY_INLET", VELOCITY_INLET)
("PRESSURE_INLET", PRESSURE_INLET);

/*!
 * \brief types outlet boundary treatments
 */
enum OUTLET_TYPE {
  PRESSURE_OUTLET = 1,    /*!< \brief Gauge pressure outlet for incompressible flow */
  MASS_FLOW_OUTLET = 2,   /*!< \brief Mass flow outlet for incompressible flow. */
};
static const map<string, OUTLET_TYPE> Outlet_Map = CCreateMap<string, OUTLET_TYPE>
("PRESSURE_OUTLET", PRESSURE_OUTLET)
("MASS_FLOW_OUTLET", MASS_FLOW_OUTLET);

/*!
 * \brief types engine inflow boundary treatments
 */
enum ENGINE_INFLOW_TYPE {
  FAN_FACE_MACH = 1,	       /*!< \brief User specifies fan face mach number. */
  FAN_FACE_MDOT = 2,           /*!< \brief User specifies Static pressure. */
  FAN_FACE_PRESSURE = 3        /*!< \brief User specifies Static pressure. */
};
static const map<string, ENGINE_INFLOW_TYPE> Engine_Inflow_Map = CCreateMap<string, ENGINE_INFLOW_TYPE>
("FAN_FACE_MACH", FAN_FACE_MACH)
("FAN_FACE_MDOT", FAN_FACE_MDOT)
("FAN_FACE_PRESSURE", FAN_FACE_PRESSURE);

/*!
 * \brief types actuator disk boundary treatments
 */
enum ACTDISK_TYPE {
  VARIABLES_JUMP = 1,	  /*!< \brief User specifies the variables jump. */
  BC_THRUST = 2,          /*!< \brief User specifies the BC thrust. */
  NET_THRUST = 3,         /*!< \brief User specifies the Net thrust. */
  DRAG_MINUS_THRUST = 4,  /*!< \brief User specifies the D-T. */
  MASSFLOW = 5,           /*!< \brief User specifies the massflow. */
  POWER = 6               /*!< \brief User specifies the power. */
};
static const map<string, ACTDISK_TYPE> ActDisk_Map = CCreateMap<string, ACTDISK_TYPE>
("VARIABLES_JUMP", VARIABLES_JUMP)
("BC_THRUST", BC_THRUST)
("NET_THRUST", NET_THRUST)
("DRAG_MINUS_THRUST", DRAG_MINUS_THRUST)
("MASSFLOW", MASSFLOW)
("POWER", POWER);

/*!
 * \brief types of geometric entities based on VTK nomenclature
 */
enum GEO_TYPE {
  VERTEX = 1,   	  /*!< \brief VTK nomenclature for defining a vertex element. */
  LINE = 3,			  /*!< \brief VTK nomenclature for defining a line element. */
  TRIANGLE = 5, 	  /*!< \brief VTK nomenclature for defining a triangle element. */
  QUADRILATERAL = 9,  /*!< \brief VTK nomenclature for defining a quadrilateral element. */
  TETRAHEDRON = 10,   /*!< \brief VTK nomenclature for defining a tetrahedron element. */
  HEXAHEDRON = 12,    /*!< \brief VTK nomenclature for defining a hexahedron element. */
  PRISM = 13,     	  /*!< \brief VTK nomenclature for defining a prism element. */
  PYRAMID = 14  	  /*!< \brief VTK nomenclature for defining a pyramid element. */
};

/*!
 * \brief types of objective functions
 */
enum ENUM_OBJECTIVE {
  DRAG_COEFFICIENT = 1, 	    /*!< \brief Drag objective function definition. */
  LIFT_COEFFICIENT = 2, 	    /*!< \brief Lift objective function definition. */
  SIDEFORCE_COEFFICIENT = 3,    /*!< \brief Side force objective function definition. */
  EFFICIENCY = 4,		        /*!< \brief Efficiency objective function definition. */
  INVERSE_DESIGN_PRESSURE = 5,	/*!< \brief Pressure objective function definition (inverse design). */
  INVERSE_DESIGN_HEATFLUX = 6,  /*!< \brief Heat flux objective function definition (inverse design). */
  TOTAL_HEATFLUX = 7,           /*!< \brief Total heat flux. */
  MAXIMUM_HEATFLUX = 8,         /*!< \brief Maximum heat flux. */
  TOTAL_AVG_TEMPERATURE = 70,   /*!< \brief Total averaged temperature. */
  MOMENT_X_COEFFICIENT = 9,	    /*!< \brief Pitching moment objective function definition. */
  MOMENT_Y_COEFFICIENT = 10,    /*!< \brief Rolling moment objective function definition. */
  MOMENT_Z_COEFFICIENT = 11,    /*!< \brief Yawing objective function definition. */
  EQUIVALENT_AREA = 12,		    /*!< \brief Equivalent area objective function definition. */
  NEARFIELD_PRESSURE = 13,	    /*!< \brief NearField Pressure objective function definition. */
  FORCE_X_COEFFICIENT = 14,	    /*!< \brief X-direction force objective function definition. */
  FORCE_Y_COEFFICIENT = 15,	    /*!< \brief Y-direction force objective function definition. */
  FORCE_Z_COEFFICIENT = 16,	    /*!< \brief Z-direction force objective function definition. */
  THRUST_COEFFICIENT = 17,		/*!< \brief Thrust objective function definition. */
  TORQUE_COEFFICIENT = 18,		/*!< \brief Torque objective function definition. */
  FIGURE_OF_MERIT = 19,		    /*!< \brief Rotor Figure of Merit objective function definition. */
  BUFFET_SENSOR = 20,           /*!< \brief Sensor for detecting separation. */
  SURFACE_TOTAL_PRESSURE = 28,  /*!< \brief Total Pressure objective function definition. */
  SURFACE_STATIC_PRESSURE = 29, /*!< \brief Static Pressure objective function definition. */
  SURFACE_MASSFLOW = 30,        /*!< \brief Mass Flow Rate objective function definition. */
  SURFACE_MACH = 51,            /*!< \brief Mach number objective function definition. */
  SURFACE_UNIFORMITY = 52,      /*!< \brief Flow uniformity objective function definition. */
  SURFACE_SECONDARY = 53,       /*!< \brief Secondary flow strength objective function definition. */
  SURFACE_MOM_DISTORTION = 54,  /*!< \brief Momentum distortion objective function definition. */
  SURFACE_SECOND_OVER_UNIFORM = 55, /*!< \brief Secondary over uniformity (relative secondary strength) objective function definition. */
  SURFACE_PRESSURE_DROP = 56, 	    /*!< \brief Pressure drop objective function definition. */
  CUSTOM_OBJFUNC = 31, 	            /*!< \brief Custom objective function definition. */
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
  REFERENCE_GEOMETRY=60,          /*!< \brief Norm of displacements with respect to target geometry. */
  REFERENCE_NODE=61,              /*!< \brief Objective function defined as the difference of a particular node respect to a reference position. */
  VOLUME_FRACTION=62,             /*!< \brief Volume average physical density, for material-based topology optimization applications. */
  TOPOL_DISCRETENESS=63,          /*!< \brief Measure of the discreteness of the current topology. */
  TOPOL_COMPLIANCE=64             /*!< \brief Measure of the discreteness of the current topology. */
};
static const map<string, ENUM_OBJECTIVE> Objective_Map = CCreateMap<string, ENUM_OBJECTIVE>
("DRAG", DRAG_COEFFICIENT)
("LIFT", LIFT_COEFFICIENT)
("SIDEFORCE", SIDEFORCE_COEFFICIENT)
("EFFICIENCY", EFFICIENCY)
("INVERSE_DESIGN_PRESSURE", INVERSE_DESIGN_PRESSURE)
("INVERSE_DESIGN_HEATFLUX", INVERSE_DESIGN_HEATFLUX)
("MOMENT_X", MOMENT_X_COEFFICIENT)
("MOMENT_Y", MOMENT_Y_COEFFICIENT)
("MOMENT_Z", MOMENT_Z_COEFFICIENT)
("EQUIVALENT_AREA", EQUIVALENT_AREA)
("NEARFIELD_PRESSURE", NEARFIELD_PRESSURE)
("FORCE_X", FORCE_X_COEFFICIENT)
("FORCE_Y", FORCE_Y_COEFFICIENT)
("FORCE_Z", FORCE_Z_COEFFICIENT)
("THRUST", THRUST_COEFFICIENT)
("TORQUE", TORQUE_COEFFICIENT)
("TOTAL_HEATFLUX", TOTAL_HEATFLUX)
("MAXIMUM_HEATFLUX", MAXIMUM_HEATFLUX)
("TOTAL_AVG_TEMPERATURE", TOTAL_AVG_TEMPERATURE)
("FIGURE_OF_MERIT", FIGURE_OF_MERIT)
("BUFFET", BUFFET_SENSOR)
("SURFACE_TOTAL_PRESSURE", SURFACE_TOTAL_PRESSURE)
("SURFACE_STATIC_PRESSURE", SURFACE_STATIC_PRESSURE)
("SURFACE_MASSFLOW", SURFACE_MASSFLOW)
("SURFACE_MACH", SURFACE_MACH)
("SURFACE_UNIFORMITY", SURFACE_UNIFORMITY)
("SURFACE_SECONDARY", SURFACE_SECONDARY)
("SURFACE_MOM_DISTORTION", SURFACE_MOM_DISTORTION)
("SURFACE_SECOND_OVER_UNIFORM", SURFACE_SECOND_OVER_UNIFORM)
("SURFACE_PRESSURE_DROP", SURFACE_PRESSURE_DROP)
("CUSTOM_OBJFUNC", CUSTOM_OBJFUNC)
("TOTAL_EFFICIENCY", TOTAL_EFFICIENCY)
("TOTAL_STATIC_EFFICIENCY", TOTAL_STATIC_EFFICIENCY)
("TOTAL_PRESSURE_LOSS", TOTAL_PRESSURE_LOSS)
("EULERIAN_WORK", EULERIAN_WORK)
("TOTAL_ENTHALPY_IN", TOTAL_ENTHALPY_IN)
("FLOW_ANGLE_IN", FLOW_ANGLE_IN)
("FLOW_ANGLE_OUT", FLOW_ANGLE_OUT)
("MASS_FLOW_IN", MASS_FLOW_IN)
("MASS_FLOW_OUT", MASS_FLOW_OUT)
("PRESSURE_RATIO",  PRESSURE_RATIO)
("ENTROPY_GENERATION",  ENTROPY_GENERATION)
("KINETIC_ENERGY_LOSS", KINETIC_ENERGY_LOSS)
("REFERENCE_GEOMETRY", REFERENCE_GEOMETRY)
("REFERENCE_NODE", REFERENCE_NODE)
("VOLUME_FRACTION", VOLUME_FRACTION)
("TOPOL_DISCRETENESS", TOPOL_DISCRETENESS)
("TOPOL_COMPLIANCE", TOPOL_COMPLIANCE);

/*!
 * \brief types of residual criteria equations
 */
enum ENUM_RESIDUAL {
    RHO_RESIDUAL = 1, 	     /*!< \brief Rho equation residual criteria equation. */
    RHO_ENERGY_RESIDUAL = 2  /*!< \brief RhoE equation residual criteria equation. */
};
static const map<string, ENUM_RESIDUAL> Residual_Map = CCreateMap<string, ENUM_RESIDUAL>
("RHO", RHO_RESIDUAL)
("RHO_ENERGY", RHO_ENERGY_RESIDUAL);

/*!
 * \brief types of residual criteria for structural problems
 */
enum ENUM_RESFEM {
  RESFEM_RELATIVE = 1,         /*!< \brief Relative criteria: Res/Res0. */
  RESFEM_ABSOLUTE = 2          /*!< \brief Absolute criteria: abs(Res). */
};
static const map<string, ENUM_RESFEM> ResFem_Map = CCreateMap<string, ENUM_RESFEM>
("RELATIVE", RESFEM_RELATIVE)
("ABSOLUTE", RESFEM_ABSOLUTE);

/*!
 * \brief types of sensitivities to compute
 */
enum ENUM_SENS {
  SENS_GEOMETRY = 1,   	/*!< \brief Geometrical sensitivity. */
  SENS_MACH = 2,		/*!< \brief Mach number sensitivity. */
  SENS_AOA = 3, 		/*!< \brief Angle of attack sensitivity. */
  SENS_AOS = 4  		/*!< \brief Angle of Sideslip sensitivity. */
};
static const map<string, ENUM_SENS> Sens_Map = CCreateMap<string, ENUM_SENS>
("SENS_GEOMETRY", SENS_GEOMETRY)
("SENS_MACH", SENS_MACH)
("SENS_AOA", SENS_AOA)
("SENS_AOS", SENS_AOS);

/*!
 * \brief types of grid adaptation/refinement
 */
enum ENUM_ADAPT {
  NO_ADAPT = 0,           /*!< \brief No grid adaptation. */
  FULL = 1,		          /*!< \brief Do a complete grid refinement of all the computational grids. */
  FULL_FLOW = 2,		  /*!< \brief Do a complete grid refinement of the flow grid. */
  FULL_ADJOINT = 3,		  /*!< \brief Do a complete grid refinement of the adjoint grid. */
  GRAD_FLOW = 5,		  /*!< \brief Do a gradient based grid adaptation of the flow grid. */
  GRAD_ADJOINT = 6,		  /*!< \brief Do a gradient based grid adaptation of the adjoint grid. */
  GRAD_FLOW_ADJ = 7,	  /*!< \brief Do a gradient based grid adaptation of the flow and adjoint grid. */
  COMPUTABLE = 9,		  /*!< \brief Apply a computable error grid adaptation. */
  REMAINING = 10,		  /*!< \brief Apply a remaining error grid adaptation. */
  WAKE = 12,			  /*!< \brief Do a grid refinement on the wake. */
  SMOOTHING = 14,		  /*!< \brief Do a grid smoothing of the geometry. */
  SUPERSONIC_SHOCK = 15,  /*!< \brief Do a grid smoothing. */
  PERIODIC = 17			  /*!< \brief Add the periodic halo cells. */
};
static const map<string, ENUM_ADAPT> Adapt_Map = CCreateMap<string, ENUM_ADAPT>
("NONE", NO_ADAPT)
("FULL", FULL)
("FULL_FLOW", FULL_FLOW)
("FULL_ADJOINT", FULL_ADJOINT)
("GRAD_FLOW", GRAD_FLOW)
("GRAD_ADJOINT", GRAD_ADJOINT)
("GRAD_FLOW_ADJ", GRAD_FLOW_ADJ)
("COMPUTABLE", COMPUTABLE)
("REMAINING", REMAINING)
("WAKE", WAKE)
("SMOOTHING", SMOOTHING)
("SUPERSONIC_SHOCK", SUPERSONIC_SHOCK)
("PERIODIC", PERIODIC);

/*!
 * \brief types of input file formats
 */
enum ENUM_INPUT {
  SU2       = 1,  /*!< \brief SU2 input format. */
  CGNS_GRID = 2,  /*!< \brief CGNS input format for the computational grid. */
  RECTANGLE = 3,  /*!< \brief 2D rectangular mesh with N x M points of size Lx x Ly. */
  BOX       = 4   /*!< \brief 3D box mesh with N x M x L points of size Lx x Ly x Lz. */
};
static const map<string, ENUM_INPUT> Input_Map = CCreateMap<string, ENUM_INPUT>
("SU2", SU2)
("CGNS", CGNS_GRID)
("RECTANGLE", RECTANGLE)
("BOX", BOX);

/*!
 * \brief type of solution output file formats
 */
enum ENUM_OUTPUT {
  TECPLOT = 1,  		         /*!< \brief Tecplot format for the solution output. */
  TECPLOT_BINARY = 2,            /*!< \brief Tecplot binary format for the solution output. */
  SURFACE_TECPLOT = 3,  	     /*!< \brief Tecplot format for the solution output. */
  SURFACE_TECPLOT_BINARY = 4,    /*!< \brief Tecplot binary format for the solution output. */
  CSV = 5,			             /*!< \brief Comma-separated values format for the solution output. */
  SURFACE_CSV = 6,			     /*!< \brief Comma-separated values format for the solution output. */
  PARAVIEW = 7,  		         /*!< \brief Paraview ASCII format for the solution output. */
  PARAVIEW_BINARY = 8,           /*!< \brief Paraview binary format for the solution output. */
  SURFACE_PARAVIEW = 9,  	     /*!< \brief Paraview ASCII format for the solution output. */
  SURFACE_PARAVIEW_BINARY = 10,  /*!< \brief Paraview binary format for the solution output. */
  MESH      = 11,                /*!< \brief SU2 mesh format. */
  RESTART_BINARY = 12,           /*!< \brief SU2 binary restart format. */
  RESTART_ASCII = 13,            /*!< \brief SU2 ASCII restart format. */
  CGNS = 14                      /*!< \brief CGNS format. */
};
static const map<string, ENUM_OUTPUT> Output_Map = CCreateMap<string, ENUM_OUTPUT>
("TECPLOT_ASCII", TECPLOT)
("TECPLOT", TECPLOT_BINARY)
("SURFACE_TECPLOT_ASCII", SURFACE_TECPLOT)
("SURFACE_TECPLOT", SURFACE_TECPLOT_BINARY)
("CSV", CSV)
("SURFACE_CSV", SURFACE_CSV)
("PARAVIEW_ASCII", PARAVIEW)
("PARAVIEW", PARAVIEW_BINARY)
("SURFACE_PARAVIEW_ASCII", SURFACE_PARAVIEW)
("SURFACE_PARAVIEW", SURFACE_PARAVIEW_BINARY)
("RESTART_ASCII", RESTART_ASCII)
("RESTART", RESTART_BINARY)
("CGNS", CGNS);

/*!
 * \brief type of solution output file formats
 */
enum ENUM_TAB_OUTPUT {
  TAB_CSV = 1,			    /*!< \brief Comma-separated values format for the solution output. */
  TAB_TECPLOT = 2           /*!< \brief Tecplot format for the solution output. */
};
static const map<string, ENUM_TAB_OUTPUT> TabOutput_Map = CCreateMap<string, ENUM_TAB_OUTPUT>
("CSV", TAB_CSV)
("TECPLOT", TAB_TECPLOT);

/*!
 * \brief type of volume sensitivity file formats (inout to SU2_DOT)
 */
enum ENUM_SENSITIVITY {
  SU2_NATIVE = 1,       /*!< \brief SU2 native binary format for the volume sensitivity input. */
  UNORDERED_ASCII = 2   /*!< \brief Unordered ASCII list (x,y,z,dJ/dx,dJ/dy/dJ/dz) format for the volume sensitivity input. */
};
static const map<string, ENUM_SENSITIVITY> Sensitivity_Map = CCreateMap<string, ENUM_SENSITIVITY>
("SU2_NATIVE", SU2_NATIVE)
("UNORDERED_ASCII", UNORDERED_ASCII);

/*!
 * \brief type of jump definition
 */
enum JUMP_DEFINITION {
  DIFFERENCE = 1,     /*!< \brief Jump given by a difference in values. */
  RATIO = 2           /*!< \brief Jump given by a ratio. */
};
static const map<string, JUMP_DEFINITION> Jump_Map = CCreateMap<string, JUMP_DEFINITION>
("DIFFERENCE", DIFFERENCE)
("RATIO", RATIO);

/*!
 * \brief type of multigrid cycle
 */
enum MG_CYCLE {
  V_CYCLE = 0,  		/*!< \brief V cycle. */
  W_CYCLE = 1,			/*!< \brief W cycle. */
  FULLMG_CYCLE = 2      /*!< \brief FullMG cycle. */
};
static const map<string, MG_CYCLE> MG_Cycle_Map = CCreateMap<string, MG_CYCLE>
("V_CYCLE", V_CYCLE)
("W_CYCLE", W_CYCLE)
("FULLMG_CYCLE", FULLMG_CYCLE);

/*!
 * \brief type of solution output variables
 */
enum ENUM_OUTPUT_VARS {
  DENSITY = 1,      /*!< \brief Density. */
  VEL_X = 2,        /*!< \brief X-component of velocity. */
  VEL_Y = 3,        /*!< \brief Y-component of velocity. */
  VEL_Z = 4,        /*!< \brief Z-component of velocity. */
  PRESSURE = 5, 	/*!< \brief Static pressure. */
  MACH = 6,         /*!< \brief Mach number. */
  TEMPERATURE = 7,  /*!< \brief Temperature. */
  LAM_VISC = 8,     /*!< \brief Laminar viscosity. */
  EDDY_VISC = 9     /*!< \brief Eddy viscosity. */
};
static const map<string, ENUM_OUTPUT_VARS> Output_Vars_Map = CCreateMap<string, ENUM_OUTPUT_VARS>
("DENSITY", DENSITY)
("VEL_X", VEL_X)
("VEL_Y", VEL_Y)
("VEL_Z", VEL_Z)
("PRESSURE", PRESSURE)
("MACH", MACH)
("TEMPERATURE", TEMPERATURE)
("LAM_VISC", LAM_VISC)
("EDDY_VISC", EDDY_VISC);

/*!
 * \brief types of design parameterizations
 */
enum ENUM_PARAM {
  NO_DEFORMATION = 0,        /*!< \brief No deformation. */

  TRANSLATION = 1,		     /*!< \brief Surface movement as design variable. */
  ROTATION = 2,			     /*!< \brief Surface rotation as design variable. */
  SCALE = 3,			     /*!< \brief Surface rotation as design variable. */
  
  FFD_SETTING = 10,		     /*!< \brief No surface deformation. */
  FFD_CONTROL_POINT = 11,	 /*!< \brief Free form deformation for 3D design (change a control point). */
  FFD_NACELLE = 12,	         /*!< \brief Free form deformation for 3D design (change a control point). */
  FFD_GULL = 13,	         /*!< \brief Free form deformation for 3D design (change a control point). */
  FFD_CAMBER = 14,		     /*!< \brief Free form deformation for 3D design (camber change). */
  FFD_TWIST = 15,		     /*!< \brief Free form deformation for 3D design (change the twist angle of a section). */
  FFD_THICKNESS = 16,		 /*!< \brief Free form deformation for 3D design (thickness change). */
  FFD_ROTATION = 18,		 /*!< \brief Free form deformation for 3D design (rotation around a line). */
  FFD_CONTROL_POINT_2D = 19, /*!< \brief Free form deformation for 2D design (change a control point). */
  FFD_CAMBER_2D = 20,		 /*!< \brief Free form deformation for 3D design (camber change). */
  FFD_THICKNESS_2D = 21,	 /*!< \brief Free form deformation for 3D design (thickness change). */
  FFD_TWIST_2D = 22,		 /*!< \brief Free form deformation for 3D design (camber change). */
  FFD_CONTROL_SURFACE = 23,	 /*!< \brief Free form deformation for 3D design (control surface). */
  FFD_ANGLE_OF_ATTACK = 24,  /*!< \brief Angle of attack for FFD problem. */

  HICKS_HENNE = 30,	         /*!< \brief Hicks-Henne bump function for airfoil deformation. */
  PARABOLIC = 31,		     /*!< \brief Parabolic airfoil definition as design variables. */
  NACA_4DIGITS = 32,	     /*!< \brief The four digits NACA airfoil family as design variables. */
  AIRFOIL = 33,		         /*!< \brief Airfoil definition as design variables. */
  CST = 34,                  /*!< \brief CST method with Kulfan parameters for airfoil deformation. */
  SURFACE_BUMP = 35,	     /*!< \brief Surfacebump function for flat surfaces deformation. */
  SURFACE_FILE = 36,	     /*!< \brief Nodal coordinates for surface set using a file (external parameterization). */
  
  DV_EFIELD = 40,            /*!< \brief Electric field in deformable membranes. */
  DV_YOUNG = 41,
  DV_POISSON = 42,
  DV_RHO = 43,
  DV_RHO_DL = 44,
  
  TRANSLATE_GRID = 50,       /*!< \brief Translate the volume grid. */
  ROTATE_GRID = 51,          /*!< \brief Rotate the volume grid */
  SCALE_GRID = 52,           /*!< \brief Scale the volume grid. */
  
  ANGLE_OF_ATTACK = 101	     /*!< \brief Angle of attack for airfoils. */
};
static const map<string, ENUM_PARAM> Param_Map = CCreateMap<string, ENUM_PARAM>
("FFD_SETTING", FFD_SETTING)
("FFD_CONTROL_POINT_2D", FFD_CONTROL_POINT_2D)
("FFD_TWIST_2D", FFD_TWIST_2D)
("FFD_ANGLE_OF_ATTACK", FFD_ANGLE_OF_ATTACK)
("FFD_CAMBER_2D", FFD_CAMBER_2D)
("FFD_THICKNESS_2D", FFD_THICKNESS_2D)
("HICKS_HENNE", HICKS_HENNE)
("SURFACE_BUMP", SURFACE_BUMP)
("ANGLE_OF_ATTACK", ANGLE_OF_ATTACK)
("NACA_4DIGITS", NACA_4DIGITS)
("TRANSLATION", TRANSLATION)
("ROTATION", ROTATION)
("SCALE", SCALE)
("FFD_CONTROL_POINT", FFD_CONTROL_POINT)
("FFD_ROTATION", FFD_ROTATION)
("FFD_CONTROL_SURFACE", FFD_CONTROL_SURFACE)
("FFD_NACELLE", FFD_NACELLE)
("FFD_GULL", FFD_GULL)
("FFD_TWIST", FFD_TWIST)
("FFD_CAMBER", FFD_CAMBER)
("FFD_THICKNESS", FFD_THICKNESS)
("PARABOLIC", PARABOLIC)
("AIRFOIL", AIRFOIL)
("SURFACE_FILE", SURFACE_FILE)
("NO_DEFORMATION", NO_DEFORMATION)
("CST", CST)
("ELECTRIC_FIELD", DV_EFIELD)
("YOUNG_MODULUS", DV_YOUNG)
("POISSON_RATIO", DV_POISSON)
("STRUCTURAL_DENSITY", DV_RHO)
("DEAD_WEIGHT", DV_RHO_DL)
("TRANSLATE_GRID", TRANSLATE_GRID)
("ROTATE_GRID", ROTATE_GRID)
("SCALE_GRID", SCALE_GRID)
;

/*!
 * \brief types of FFD Blending function
 */
enum ENUM_FFD_BLENDING{
  BSPLINE_UNIFORM = 0,  /*!< \brief BSpline blending */
  BEZIER = 1,           /*!< \brief Bezier blending */
};
static const map<string, ENUM_FFD_BLENDING> Blending_Map = CCreateMap<string, ENUM_FFD_BLENDING>
("BSPLINE_UNIFORM", BSPLINE_UNIFORM)
("BEZIER", BEZIER);

/*!
 * \brief types of solvers for solving linear systems
 */
enum ENUM_LINEAR_SOLVER {
  STEEPEST_DESCENT = 1,	   /*!< \brief Steepest descent method for point inversion algoritm (Free-Form). */
  NEWTON = 2,			   /*!< \brief Newton method for point inversion algorithm (Free-Form). */
  QUASI_NEWTON = 3,		   /*!< \brief Quasi Newton method for point inversion algorithm (Free-Form). */
  CONJUGATE_GRADIENT = 4,  /*!< \brief Preconditionated conjugate gradient method for grid deformation. */
  FGMRES = 5,    	       /*!< \brief Flexible Generalized Minimal Residual method. */
  BCGSTAB = 6,	           /*!< \brief BCGSTAB - Biconjugate Gradient Stabilized Method (main solver). */
  RESTARTED_FGMRES = 7,    /*!< \brief Flexible Generalized Minimal Residual method with restart. */
  SMOOTHER = 8,            /*!< \brief Iterative smoother. */
  PASTIX_LDLT = 9,         /*!< \brief PaStiX LDLT (complete) factorization. */
  PASTIX_LU = 10,          /*!< \brief PaStiX LU (complete) factorization. */
};
static const map<string, ENUM_LINEAR_SOLVER> Linear_Solver_Map = CCreateMap<string, ENUM_LINEAR_SOLVER>
("STEEPEST_DESCENT", STEEPEST_DESCENT)
("NEWTON", NEWTON)
("QUASI_NEWTON", QUASI_NEWTON)
("CONJUGATE_GRADIENT", CONJUGATE_GRADIENT)
("BCGSTAB", BCGSTAB)
("FGMRES", FGMRES)
("RESTARTED_FGMRES", RESTARTED_FGMRES)
("SMOOTHER", SMOOTHER)
("PASTIX_LDLT", PASTIX_LDLT)
("PASTIX_LU", PASTIX_LU);

/*!
 * \brief types surface continuity at the intersection with the FFD
 */
enum ENUM_FFD_CONTINUITY {
  DERIVATIVE_NONE = 0,	/*!< \brief No derivative continuity. */
  DERIVATIVE_1ST = 1,	/*!< \brief First derivative continuity. */
  DERIVATIVE_2ND = 2,	/*!< \brief Second derivative continuity. */
  USER_INPUT = 3		/*!< \brief User input. */
};
static const map<string, ENUM_FFD_CONTINUITY> Continuity_Map = CCreateMap<string, ENUM_FFD_CONTINUITY>
("NO_DERIVATIVE", DERIVATIVE_NONE)
("1ST_DERIVATIVE", DERIVATIVE_1ST)
("2ND_DERIVATIVE", DERIVATIVE_2ND)
("USER_INPUT", USER_INPUT);

/*!
 * \brief types of coordinates systems for the FFD
 */
enum ENUM_FFD_COORD_SYSTEM {
  CARTESIAN = 0,    /*!< \brief Cartesian coordinate system. */
  CYLINDRICAL = 1,  /*!< \brief Cylindrical coordinate system. */
  SPHERICAL = 2,    /*!< \brief Spherical coordinate system. */
  POLAR = 3         /*!< \brief Polar coordinate system. */
};
static const map<string, ENUM_FFD_COORD_SYSTEM> CoordSystem_Map = CCreateMap<string, ENUM_FFD_COORD_SYSTEM>
("CARTESIAN", CARTESIAN)
("CYLINDRICAL", CYLINDRICAL)
("SPHERICAL", SPHERICAL)
("POLAR", POLAR);

/*!
 * \brief types of sensitivity smoothing
 */
enum ENUM_SENS_SMOOTHING {
  NO_SMOOTH = 0,  /*!< \brief No smoothing. */
  SOBOLEV = 1,	  /*!< \brief Sobolev gradient smoothing. */
  BIGRID = 2	  /*!< \brief Bi-grid technique smoothing. */
};
static const map<string, ENUM_SENS_SMOOTHING> Sens_Smoothing_Map = CCreateMap<string, ENUM_SENS_SMOOTHING>
("NONE", NO_SMOOTH)
("SOBOLEV", SOBOLEV)
("BIGRID", BIGRID);

/*!
 * \brief types of preconditioners for the linear solver
 */
enum ENUM_LINEAR_SOLVER_PREC {
  JACOBI = 1,	     /*!< \brief Jacobi preconditioner. */
  LU_SGS = 2,	     /*!< \brief LU SGS preconditioner. */
  LINELET = 3,       /*!< \brief Line implicit preconditioner. */
  ILU = 4,           /*!< \brief ILU(k) preconditioner. */
  PASTIX_ILU= 5,     /*!< \brief PaStiX ILU(k) preconditioner. */
  PASTIX_LU_P= 6,    /*!< \brief PaStiX LU as preconditioner. */
  PASTIX_LDLT_P= 7,  /*!< \brief PaStiX LDLT as preconditioner. */
};
static const map<string, ENUM_LINEAR_SOLVER_PREC> Linear_Solver_Prec_Map = CCreateMap<string, ENUM_LINEAR_SOLVER_PREC>
("JACOBI", JACOBI)
("LU_SGS", LU_SGS)
("LINELET", LINELET)
("ILU", ILU)
("PASTIX_ILU", PASTIX_ILU)
("PASTIX_LU", PASTIX_LU_P)
("PASTIX_LDLT", PASTIX_LDLT_P);

/*!
 * \brief types of analytic definitions for various geometries
 */
enum ENUM_GEO_ANALYTIC {
  NO_GEO_ANALYTIC = 0,   /*!< \brief No analytic definition of the geometry. */
  NACA0012_AIRFOIL = 1,  /*!< \brief Use the analytical definition of the NACA0012 for doing the grid adaptation. */
  NACA4412_AIRFOIL = 2,  /*!< \brief Use the analytical definition of the NACA4412 for doing the grid adaptation. */
  CYLINDER = 3, 	     /*!< \brief Use the analytical definition of a cylinder for doing the grid adaptation. */
  BIPARABOLIC = 4        /*!< \brief Use the analytical definition of a biparabolic airfoil for doing the grid adaptation. */
};
static const map<string, ENUM_GEO_ANALYTIC> Geo_Analytic_Map = CCreateMap<string, ENUM_GEO_ANALYTIC>
("NONE", NO_GEO_ANALYTIC)
("NACA0012_AIRFOIL", NACA0012_AIRFOIL)
("NACA4412_AIRFOIL", NACA4412_AIRFOIL)
("CYLINDER", CYLINDER)
("BIPARABOLIC", BIPARABOLIC);

/*!
 * \brief types of axis orientation
 */
enum ENUM_GEO_DESCRIPTION {
  TWOD_AIRFOIL = 0, /*!< \brief Airfoil analysis. */
  WING = 1, 	    /*!< \brief Wing analysis. */
  FUSELAGE = 2,     /*!< \brief Fuselage analysis. */
  NACELLE = 3       /*!< \brief Nacelle analysis. */
};
static const map<string, ENUM_GEO_DESCRIPTION> Geo_Description_Map = CCreateMap<string, ENUM_GEO_DESCRIPTION>
("AIRFOIL", TWOD_AIRFOIL)
("WING", WING)
("FUSELAGE", FUSELAGE)
("NACELLE", NACELLE);

/*!
 * \brief types of schemes for unsteady computations
 */
enum ENUM_UNSTEADY {
  STEADY = 0,            /*!< \brief A steady computation. */
  TIME_STEPPING = 1,	 /*!< \brief Use a time stepping strategy for unsteady computations. */
  DT_STEPPING_1ST = 2,	 /*!< \brief Use a dual time stepping strategy for unsteady computations (1st order). */
  DT_STEPPING_2ND = 3,	 /*!< \brief Use a dual time stepping strategy for unsteady computations (2nd order). */
  ROTATIONAL_FRAME = 4,  /*!< \brief Use a rotational source term. */
  HARMONIC_BALANCE = 5   /*!< \brief Use a harmonic balance source term. */
};
static const map<string, ENUM_UNSTEADY> TimeMarching_Map = CCreateMap<string, ENUM_UNSTEADY>
("NO", STEADY)
("TIME_STEPPING", TIME_STEPPING)
("DUAL_TIME_STEPPING-1ST_ORDER", DT_STEPPING_1ST)
("DUAL_TIME_STEPPING-2ND_ORDER", DT_STEPPING_2ND)
("HARMONIC_BALANCE", HARMONIC_BALANCE)
("ROTATIONAL_FRAME", ROTATIONAL_FRAME);

/*!
 * \brief types of criteria to determine when the solution is converged
 */
enum ENUM_CONVERGE_CRIT {
  CAUCHY = 1,			/*!< \brief Cauchy criteria to establish the convergence of the code. */
  RESIDUAL = 2			/*!< \brief Residual criteria to establish the convergence of the code. */
};
static const map<string, ENUM_CONVERGE_CRIT> Converge_Crit_Map = CCreateMap<string, ENUM_CONVERGE_CRIT>
("CAUCHY", CAUCHY)
("RESIDUAL", RESIDUAL);

/*!
 * \brief types of element stiffnesses imposed for FEA mesh deformation
 */
enum ENUM_DEFORM_STIFFNESS {
  CONSTANT_STIFFNESS = 0,      /*!< \brief Impose a constant stiffness for each element (steel). */
  INVERSE_VOLUME = 1,	       /*!< \brief Impose a stiffness for each element that is inversely proportional to cell volume. */
  SOLID_WALL_DISTANCE = 2      /*!< \brief Impose a stiffness for each element that is proportional to the distance from the solid surface. */
};
static const map<string, ENUM_DEFORM_STIFFNESS> Deform_Stiffness_Map = CCreateMap<string, ENUM_DEFORM_STIFFNESS>
("CONSTANT_STIFFNESS", CONSTANT_STIFFNESS)
("INVERSE_VOLUME", INVERSE_VOLUME)
("WALL_DISTANCE", SOLID_WALL_DISTANCE);

/*!
 * \brief The direct differentation variables.
 */
enum ENUM_DIRECTDIFF_VAR {
  NO_DERIVATIVE = 0,
  D_MACH = 1,         /*!< \brief Derivative with respect to the mach number */
  D_AOA = 2,		  /*!< \brief Derivative with respect to the angle of attack */
  D_PRESSURE = 3,     /*!< \brief Derivative with respect to the freestream pressure */
  D_TEMPERATURE = 4,  /*!< \brief Derivative with respect to the freestream temperature */
  D_DENSITY = 5,      /*!< \brief Derivative with respect to the freestream density */
  D_TURB2LAM = 6,     /*!< \brief Derivative with respect to the turb2lam */
  D_SIDESLIP = 7,     /*!< \brief Derivative with respect to the sideslip angle */
  D_VISCOSITY = 8,    /*!< \brief Derivative with respect to the viscosity */
  D_REYNOLDS = 9,     /*!< \brief Derivative with respect to the reynolds number */
  D_DESIGN = 10,      /*!< \brief Derivative with respect to the design?? */
  D_YOUNG = 11,       /*!< \brief Derivative with respect to the Young's Modulus */
  D_POISSON = 12,
  D_RHO = 13,
  D_RHO_DL = 14,
  D_EFIELD = 15
};
static const map<string, ENUM_DIRECTDIFF_VAR> DirectDiff_Var_Map = CCreateMap<string, ENUM_DIRECTDIFF_VAR>
("NONE", NO_DERIVATIVE)
("MACH", D_MACH)
("AOA", D_AOA)
("PRESSURE", D_PRESSURE)
("TEMPERATURE", D_TEMPERATURE)
("DENSITY", D_DENSITY)
("TURB2LAM", D_TURB2LAM)
("SIDESLIP", D_SIDESLIP)
("VISCOSITY", D_VISCOSITY)
("REYNOLDS", D_REYNOLDS)
("DESIGN_VARIABLES", D_DESIGN)
("YOUNG_MODULUS", D_YOUNG)
("POISSON_RATIO", D_POISSON)
("STRUCTURAL_DENSITY", D_RHO)
("STRUCTURAL_DEAD_LOAD", D_RHO_DL)
("ELECTRIC_FIELD", D_EFIELD);


enum ENUM_RECORDING {
  FLOW_CONS_VARS   = 1,
  MESH_COORDS = 2,
  COMBINED    = 3,
  FEA_DISP_VARS = 4,
  FLOW_CROSS_TERM = 5,
  FEM_CROSS_TERM_GEOMETRY = 6,
  GEOMETRY_CROSS_TERM = 7,
  ALL_VARIABLES = 8,
  MESH_DEFORM = 9
};

/*!
 * \brief types of schemes for dynamic structural computations
 */
enum ENUM_DYNAMIC {
  STATIC = 0,             /*!< \brief A static structural computation. */
  DYNAMIC = 1		      /*!< \brief Use a time stepping strategy for dynamic computations. */
};
static const map<string, ENUM_DYNAMIC> Dynamic_Map = CCreateMap<string, ENUM_DYNAMIC>
("NO", STATIC)
("YES", DYNAMIC);

/*!
 * \brief types of input file formats
 */
enum ENUM_INPUT_REF {
  SU2_REF = 1,                     /*!< \brief SU2 input format (from a restart). */
  CUSTOM_REF = 2                   /*!< \brief CGNS input format for the computational grid. */
};
static const map<string, ENUM_INPUT_REF> Input_Ref_Map = CCreateMap<string, ENUM_INPUT_REF>
("SU2", SU2_REF)
("CUSTOM", CUSTOM_REF);

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
  SOLUTION_DISPONLY    =  4,  /*!< \brief Solution displacement only communication. */
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
static const map<string, COMM_LEVEL> Comm_Map = CCreateMap<string, COMM_LEVEL>
("NONE",    COMM_NONE)
("MINIMAL", COMM_MINIMAL)
("FULL",    COMM_FULL);

/*
 * \brief types of filter kernels, initially intended for structural topology optimization applications
 */
enum ENUM_FILTER_KERNEL {
  CONSTANT_WEIGHT_FILTER = 0,      /*!< \brief Uniform weight. */
  CONICAL_WEIGHT_FILTER  = 1,      /*!< \brief Linear decay with distance from center point [Bruns and Tortorelli, 2001]. */
  GAUSSIAN_WEIGHT_FILTER = 2,      /*!< \brief Bell shape around center point [Bruns and Tortorelli, 2003]. */
  DILATE_MORPH_FILTER    = 3,      /*!< \brief Continuous version of the dilate morphology operator [Sigmund 2007]. */
  ERODE_MORPH_FILTER     = 4,      /*!< \brief Continuous version of the erode morphology operator [Sigmund 2007].*/
};
static const map<string, ENUM_FILTER_KERNEL> Filter_Kernel_Map = CCreateMap<string, ENUM_FILTER_KERNEL>
("CONSTANT", CONSTANT_WEIGHT_FILTER)
("CONICAL" , CONICAL_WEIGHT_FILTER)
("GAUSSIAN", GAUSSIAN_WEIGHT_FILTER)
("DILATE"  , DILATE_MORPH_FILTER)
("ERODE"   , ERODE_MORPH_FILTER);

/*!
 * \brief types of projection function, initially intended for structural topology optimization applications
 */
enum ENUM_PROJECTION_FUNCTION {
  NO_PROJECTION  = 0,      /*!< \brief No projection. */
  HEAVISIDE_UP   = 1,      /*!< \brief Project values towards 1. */
  HEAVISIDE_DOWN = 2,      /*!< \brief Project values towards 0. */
};
static const map<string, ENUM_PROJECTION_FUNCTION> Projection_Function_Map = CCreateMap<string, ENUM_PROJECTION_FUNCTION>
("NO_PROJECTION" , NO_PROJECTION)
("HEAVISIDE_UP"  , HEAVISIDE_UP)
("HEAVISIDE_DOWN", HEAVISIDE_DOWN);

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
static const map<string, ENUM_VERIFICATION_SOLUTIONS> Verification_Solution_Map = CCreateMap<string, ENUM_VERIFICATION_SOLUTIONS>
("NO_VERIFICATION_SOLUTION", NO_VERIFICATION_SOLUTION)
("INVISCID_VORTEX",          INVISCID_VORTEX)
("RINGLEB",                  RINGLEB)
("NS_UNIT_QUAD",             NS_UNIT_QUAD)
("TAYLOR_GREEN_VORTEX",      TAYLOR_GREEN_VORTEX)
("INC_TAYLOR_GREEN_VORTEX",  INC_TAYLOR_GREEN_VORTEX)
("MMS_NS_UNIT_QUAD",         MMS_NS_UNIT_QUAD)
("MMS_NS_UNIT_QUAD_WALL_BC", MMS_NS_UNIT_QUAD_WALL_BC)
("MMS_NS_TWO_HALF_CIRCLES",  MMS_NS_TWO_HALF_CIRCLES)
("MMS_NS_TWO_HALF_SPHERES",  MMS_NS_TWO_HALF_SPHERES)
("MMS_INC_EULER",            MMS_INC_EULER)
("MMS_INC_NS",               MMS_INC_NS)
("USER_DEFINED_SOLUTION",    USER_DEFINED_SOLUTION);

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

template <class Tenum>
class COptionEnum : public COptionBase {

  map<string, Tenum> m;
  unsigned short & field; // Reference to the feildname
  Tenum def; // Default value
  string name; // identifier for the option

public:
  COptionEnum(string option_field_name, const map<string, Tenum> m, unsigned short & option_field, Tenum default_value) : field(option_field) {
    this->m = m;
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionEnum() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // Check if there is more than one string
    string out = optionCheckMultipleValues(option_value, "enum", this->name);
    if (out.compare("") != 0) {
      return out;
    }

    // Check to see if the enum value is in the map
    if (this->m.find(option_value[0]) == m.end()) {
      string str;
      str.append(this->name);
      str.append(": invalid option value ");
      str.append(option_value[0]);
      str.append(". Check current SU2 options in config_template.cfg.");
      return str;
    }
    // If it is there, set the option value
    Tenum val = this->m[option_value[0]];
    this->field = val;
    return "";
  }

  void SetDefault() {
    this->field = this->def;
  }
};

class COptionDouble : public COptionBase {
  su2double & field; // Reference to the fieldname
  su2double def; // Default value
  string name; // identifier for the option

public:
  COptionDouble(string option_field_name, su2double & option_field, su2double default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionDouble() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // check if there is more than one value
    string out = optionCheckMultipleValues(option_value, "su2double", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    su2double val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "su2double", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionString : public COptionBase {
  string & field; // Reference to the fieldname
  string def; // Default value
  string name; // identifier for the option

public:
  COptionString(string option_field_name, string & option_field, string default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionString() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // check if there is more than one value
    string out = optionCheckMultipleValues(option_value, "su2double", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    this->field.assign(option_value[0]);
    return "";
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionInt : public COptionBase {
  int & field; // Reference to the feildname
  int def; // Default value
  string name; // identifier for the option

public:
  COptionInt(string option_field_name, int & option_field, int default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionInt() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    string out = optionCheckMultipleValues(option_value, "int", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    int val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "int", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionULong : public COptionBase {
  unsigned long & field; // Reference to the feildname
  unsigned long def; // Default value
  string name; // identifier for the option

public:
  COptionULong(string option_field_name, unsigned long & option_field, unsigned long default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionULong() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    string out = optionCheckMultipleValues(option_value, "unsigned long", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    unsigned long val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "unsigned long", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionUShort : public COptionBase {
  unsigned short & field; // Reference to the feildname
  unsigned short def; // Default value
  string name; // identifier for the option

public:
  COptionUShort(string option_field_name, unsigned short & option_field, unsigned short default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionUShort() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    unsigned short val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "unsigned short", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionLong : public COptionBase {
  long & field; // Reference to the feildname
  long def; // Default value
  string name; // identifier for the option

public:
  COptionLong(string option_field_name, long & option_field, long default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionLong() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    string out = optionCheckMultipleValues(option_value, "long", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    long val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "long", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionBool : public COptionBase {
  bool & field; // Reference to the feildname
  bool def; // Default value
  string name; // identifier for the option

public:
  COptionBool(string option_field_name, bool & option_field, bool default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionBool() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // check if there is more than one value
    string out = optionCheckMultipleValues(option_value, "bool", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    if (option_value[0].compare("YES") == 0) {
      this->field = true;
      return "";
    }
    if (option_value[0].compare("NO") == 0) {
      this->field = false;
      return "";
    }
    return badValue(option_value, "bool", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

template <class Tenum>
class COptionEnumList : public COptionBase {

  map<string, Tenum> m;
  unsigned short * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionEnumList(string option_field_name, const map<string, Tenum> m, unsigned short * & option_field, unsigned short & list_size) : field(option_field) , size(list_size) {
    this->m = m;
    this->name = option_field_name;
  }

  ~COptionEnumList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    if (option_value.size() == 1 && option_value[0].compare("NONE") == 0) {
      this->size = 0;
      return "";
    }
    // size is the length of the option list
    this->size = option_value.size();
    unsigned short * enums = new unsigned short[size];
    for (int i  = 0; i < this->size; i++) {
      // Check to see if the enum value is in the map
      if (this->m.find(option_value[i]) == m.end()) {
        string str;
        str.append(this->name);
        str.append(": invalid option value ");
        str.append(option_value[i]);
        str.append(". Check current SU2 options in config_template.cfg.");
        return str;
      }
      // If it is there, set the option value
      enums[i] = this->m[option_value[i]];
    }
    this->field = enums;
    return "";
  }

  void SetDefault() {
    // No default to set
    size = 0;
  }
};

class COptionDoubleArray : public COptionBase {
  su2double * & field; // Reference to the feildname
  string name; // identifier for the option
  const int size;
  su2double * def;
  su2double * vals;
  su2double * default_value;

public:
  COptionDoubleArray(string option_field_name, const int list_size, su2double * & option_field, su2double * default_value) : field(option_field), size(list_size) {
    this->name = option_field_name;
    this->default_value = default_value;
    def  = NULL;
    vals = NULL;
  }

  ~COptionDoubleArray() {
     if(def  != NULL) delete [] def; 
     if(vals != NULL) delete [] vals; 
  };
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // Check that the size is correct
    if (option_value.size() != (unsigned long)this->size) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": wrong number of arguments: ");
      stringstream ss;
      ss << this->size;
      newstring.append(ss.str());
      newstring.append(" expected, ");
      stringstream ss2;
      ss2 << option_value.size();
      newstring.append(ss2.str());
      newstring.append(" found");
      return newstring;
    }
    vals = new su2double[this->size];
    for (int i  = 0; i < this->size; i++) {
      istringstream is(option_value[i]);
      su2double val;
      if (!(is >> val)) {
        delete [] vals;
        return badValue(option_value, "su2double array", this->name);
      }
      vals[i] = val;
    }
    this->field = vals;
    return "";
  }

  void SetDefault() {
    def = new su2double [size];
    for (int i = 0; i < size; i++) {
      def[i] = default_value[i];
    }
    this->field = def;
  }
};

class COptionDoubleList : public COptionBase {
  su2double * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionDoubleList(string option_field_name, unsigned short & list_size, su2double * & option_field) : field(option_field), size(list_size) {
    this->name = option_field_name;
  }

  ~COptionDoubleList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // The size is the length of option_value
    unsigned short option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE") == 0) {
      // No options
      this->size = 0;
      return "";
    }

    this->size = option_size;

    // Parse all of the options
    su2double * vals = new su2double[option_size];
    for (unsigned long i  = 0; i < option_size; i++) {
      istringstream is(option_value[i]);
      su2double val;
      if (!(is >> val)) {
        delete [] vals;
        return badValue(option_value, "su2double list", this->name);
      }
      vals[i] = val;
    }
    this->field = vals;
    return "";
  }

  void SetDefault() {
    this->size = 0; // There is no default value for list
  }
};

class COptionShortList : public COptionBase {
  short * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;
  
public:
  COptionShortList(string option_field_name, unsigned short & list_size,  short * & option_field) : field(option_field), size(list_size) {
    this->name = option_field_name;
  }
  
  ~COptionShortList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // The size is the length of option_value
    unsigned short option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE") == 0) {
      // No options
      this->size = 0;
      return "";
    }
    this->size = option_size;
    
    // Parse all of the options
    short * vals = new  short[option_size];
    for (unsigned long i  = 0; i < option_size; i++) {
      istringstream is(option_value[i]);
      unsigned short val;
      if (!(is >> val)) {
        delete [] vals;
        return badValue(option_value, "short", this->name);
      }
      vals[i] = val;
    }
    this->field = vals;
    return "";
  }
  
  void SetDefault() {
    this->size = 0; // There is no default value for list
  }
};

class COptionUShortList : public COptionBase {
  unsigned short * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionUShortList(string option_field_name, unsigned short & list_size, unsigned short * & option_field) : field(option_field), size(list_size) {
    this->name = option_field_name;
  }

  ~COptionUShortList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // The size is the length of option_value
    unsigned short option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE") == 0) {
      // No options
      this->size = 0;
      return "";
    }
    this->size = option_size;

    // Parse all of the options
    unsigned short * vals = new unsigned short[option_size];
    for (unsigned long i  = 0; i < option_size; i++) {
      istringstream is(option_value[i]);
      unsigned short val;
      if (!(is >> val)) {
        delete [] vals;
        return badValue(option_value, "unsigned short", this->name);
      }
      vals[i] = val;
    }
    this->field = vals;
    return "";
  }

  void SetDefault() {
    this->size = 0; // There is no default value for list
  }
};

class COptionStringList : public COptionBase {
  string * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionStringList(string option_field_name, unsigned short & list_size, string * & option_field) : field(option_field), size(list_size) {
    this->name = option_field_name;
  }

  ~COptionStringList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // The size is the length of option_value
    unsigned short option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE") == 0) {
      this->size = 0;
      return "";
    }
    this->size = option_size;

    // Parse all of the options
    string * vals = new string[option_size];
    for (unsigned long i  = 0; i < option_size; i++) {
      vals[i].assign(option_value[i]);
    }
    this->field = vals;
    return "";
  }

  void SetDefault() {
    this->size = 0; // There is no default value for list
  }
};

class COptionConvect : public COptionBase {
  string name; // identifier for the option
  unsigned short & space;
  unsigned short & centered;
  unsigned short & upwind;

public:
  COptionConvect(string option_field_name, unsigned short & space_field, unsigned short & centered_field, unsigned short & upwind_field) : space(space_field), centered(centered_field), upwind(upwind_field) {
    this->name = option_field_name;
  }

  ~COptionConvect() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);

    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0) {
      return out;
    }

    if (Centered_Map.count(option_value[0])) {
      this->space = Space_Map.find("SPACE_CENTERED")->second;
      this->centered = Centered_Map.find(option_value[0])->second;
      this->upwind = NO_UPWIND;
      return "";
    }
    if (Upwind_Map.count(option_value[0])) {
      this->space = Space_Map.find("SPACE_UPWIND")->second;
      this->upwind = Upwind_Map.find(option_value[0])->second;
      this->centered = NO_CENTERED;
      return "";
    }
    // Make them defined in case something weird happens
    this->centered = NO_CENTERED;
    this->upwind = NO_UPWIND;
    this->space = SPACE_CENTERED;
    return badValue(option_value, "convect", this->name);

  }

  void SetDefault() {
    this->centered = NO_CENTERED;
    this->upwind = NO_UPWIND;
    this->space = SPACE_CENTERED;
  }
};

class COptionFEMConvect : public COptionBase{
  string name; // identifier for the option
  unsigned short & space;
  unsigned short & fem;

public:
  COptionFEMConvect(string option_field_name, unsigned short & space_field, unsigned short & fem_field) : space(space_field), fem(fem_field) {
    this->name = option_field_name;
  }

  ~COptionFEMConvect() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);

    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0) {
      return out;
    }

    if (FEM_Map.count(option_value[0])) {
      this->space = Space_Map.find("FINITE_ELEMENT")->second;
      this->fem = FEM_Map.find(option_value[0])->second;
      return "";
    }

    // Make them defined in case something weird happens
    this->fem = NO_FEM;
    return badValue(option_value, "convect", this->name);

  }

  void SetDefault() {
    this->fem = NO_FEM;
  }
};

class COptionMathProblem : public COptionBase {
  string name; // identifier for the option
  bool & cont_adjoint;
  bool cont_adjoint_def;
  bool & disc_adjoint;
  bool disc_adjoint_def;
  bool & restart;
  bool restart_def;

public:
  COptionMathProblem(string option_field_name, bool & cont_adjoint_field, bool cont_adjoint_default, bool & disc_adjoint_field, bool disc_adjoint_default, bool & restart_field, bool restart_default) : cont_adjoint(cont_adjoint_field), disc_adjoint(disc_adjoint_field), restart(restart_field) {
    this->name = option_field_name;
    this->cont_adjoint_def = cont_adjoint_default;
    this->disc_adjoint_def = disc_adjoint_default;
    this->restart_def = restart_default;
  }

  ~COptionMathProblem() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    if (option_value[0] == "ADJOINT") {
      return badValue(option_value, "math problem (try CONTINUOUS_ADJOINT)", this->name);
    }
    if (Math_Problem_Map.find(option_value[0]) == Math_Problem_Map.end()) {
      return badValue(option_value, "math problem", this->name);
    }
    if (option_value[0] == "DIRECT") {
      this->cont_adjoint = false;
      this->disc_adjoint = false;
      this->restart = false;
      return "";
    }
    if (option_value[0] == "CONTINUOUS_ADJOINT") {
      this->cont_adjoint= true;
      this->disc_adjoint = false;
      this->restart= true;
      return "";
    }
    if (option_value[0] == "DISCRETE_ADJOINT") {
      this->disc_adjoint = true;
      this->cont_adjoint= false;
      this->restart = true;
      return "";
    }
    return "option in math problem map not considered in constructor";
  }

  void SetDefault() {
    this->cont_adjoint = this->cont_adjoint_def;
    this->disc_adjoint = this->disc_adjoint_def;
    this->restart = this->restart_def;
  }
  
};

class COptionDVParam : public COptionBase {
  string name; // identifier for the option
  unsigned short & nDV;
  su2double ** & paramDV;
  string * & FFDTag;
  unsigned short* & design_variable;

public:
  COptionDVParam(string option_field_name, unsigned short & nDV_field, su2double** & paramDV_field, string* & FFDTag_field, unsigned short * & design_variable_field) : nDV(nDV_field), paramDV(paramDV_field), FFDTag(FFDTag_field), design_variable(design_variable_field) {
    this->name = option_field_name;
  }

  ~COptionDVParam() {};
  
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    if ((option_value.size() == 1) && (option_value[0].compare("NONE") == 0)) {
      this->nDV = 0;
      return "";
    }

    // Cannot have ; at the beginning or the end
    if (option_value[0].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have beginning semicolon");
      return newstring;
    }
    if (option_value[option_value.size()-1].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have ending semicolon");
      return newstring;
    }


    // use the ";" token to determine the number of design variables
    // This works because semicolon is not one of the delimiters in tokenize string
    this->nDV = 0;
    //unsigned int num_semi = 0;
    for (unsigned int i = 0; i < static_cast<unsigned int>(option_value.size()); i++) {
      if (option_value[i].compare(";") == 0) {
        this->nDV++;
        //        num_semi++;
      }
    }

    // One more design variable than semicolon
    this->nDV++;

    if ( (this->nDV > 0) && (this->design_variable == NULL) ) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": Design_Variable array has not been allocated. Check that DV_KIND appears before DV_PARAM in configuration file.");
      return newstring;
    }

    this->paramDV = new su2double*[this->nDV];
    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) {
      this->paramDV[iDV] = new su2double[MAX_PARAMETERS];
    }

    this->FFDTag = new string[this->nDV];

   vector<unsigned short> nParamDV(nDV, 0);
   unsigned short totalnParamDV = 0;
   stringstream ss;
   unsigned int i = 0;
    
    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) {
      switch (this->design_variable[iDV]) {
        case NO_DEFORMATION:       nParamDV[iDV] = 0; break;
        case FFD_SETTING:          nParamDV[iDV] = 0; break;
        case FFD_CONTROL_POINT_2D: nParamDV[iDV] = 5; break;
        case FFD_CAMBER_2D:        nParamDV[iDV] = 2; break;
        case FFD_THICKNESS_2D:     nParamDV[iDV] = 2; break;
        case FFD_TWIST_2D:         nParamDV[iDV] = 3; break;
        case HICKS_HENNE:          nParamDV[iDV] = 2; break;
        case SURFACE_BUMP:         nParamDV[iDV] = 3; break;
        case CST:                  nParamDV[iDV] = 3; break;
        case ANGLE_OF_ATTACK:      nParamDV[iDV] = 1; break;
        case SCALE:                nParamDV[iDV] = 0; break;
        case TRANSLATION:          nParamDV[iDV] = 3; break;
        case ROTATION:             nParamDV[iDV] = 6; break;
        case NACA_4DIGITS:         nParamDV[iDV] = 3; break;
        case PARABOLIC:            nParamDV[iDV] = 2; break;
        case AIRFOIL:              nParamDV[iDV] = 2; break;
        case FFD_CONTROL_POINT:    nParamDV[iDV] = 7; break;
        case FFD_NACELLE:          nParamDV[iDV] = 6; break;
        case FFD_GULL:             nParamDV[iDV] = 2; break;
        case FFD_TWIST:            nParamDV[iDV] = 8; break;
        case FFD_ROTATION:         nParamDV[iDV] = 7; break;
        case FFD_CONTROL_SURFACE:  nParamDV[iDV] = 7; break;
        case FFD_CAMBER:           nParamDV[iDV] = 3; break;
        case FFD_THICKNESS:        nParamDV[iDV] = 3; break;
        case FFD_ANGLE_OF_ATTACK:  nParamDV[iDV] = 2; break;
        case SURFACE_FILE:         nParamDV[iDV] = 0; break;
        case DV_EFIELD:            nParamDV[iDV] = 2; break;
        case DV_YOUNG:             nParamDV[iDV] = 0; break;
        case DV_POISSON:           nParamDV[iDV] = 0; break;
        case DV_RHO:               nParamDV[iDV] = 0; break;
        case DV_RHO_DL:            nParamDV[iDV] = 0; break;
        case SCALE_GRID:           nParamDV[iDV] = 0; break;
        case TRANSLATE_GRID:       nParamDV[iDV] = 3; break;
        case ROTATE_GRID:          nParamDV[iDV] = 6; break;
        default : {
          string newstring;
          newstring.append(this->name);
          newstring.append(": undefined design variable type found in configuration file.");
          return newstring;
        }
      }
      totalnParamDV += nParamDV[iDV];
    }
    
    if (totalnParamDV > option_value.size()){
      SU2_MPI::Error("Wrong number of arguments for DV_PARAM!", CURRENT_FUNCTION);
    }
    
    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) { 
      for (unsigned short iParamDV = 0; iParamDV < nParamDV[iDV]; iParamDV++) {

        ss << option_value[i] << " ";

        if ((iParamDV == 0) &&
            ((this->design_variable[iDV] == NO_DEFORMATION) ||
             (this->design_variable[iDV] == FFD_SETTING) ||
             (this->design_variable[iDV] == FFD_ANGLE_OF_ATTACK)||
             (this->design_variable[iDV] == FFD_CONTROL_POINT_2D) ||
             (this->design_variable[iDV] == FFD_CAMBER_2D) ||
             (this->design_variable[iDV] == FFD_TWIST_2D) ||
             (this->design_variable[iDV] == FFD_THICKNESS_2D) ||
             (this->design_variable[iDV] == FFD_CONTROL_POINT) ||
             (this->design_variable[iDV] == FFD_NACELLE) ||
             (this->design_variable[iDV] == FFD_GULL) ||
             (this->design_variable[iDV] == FFD_TWIST) ||
             (this->design_variable[iDV] == FFD_ROTATION) ||
             (this->design_variable[iDV] == FFD_CONTROL_SURFACE) ||
             (this->design_variable[iDV] == FFD_CAMBER) ||
             (this->design_variable[iDV] == FFD_THICKNESS))) {
              ss >> this->FFDTag[iDV];
              this->paramDV[iDV][iParamDV] = 0;
            }
        else
          ss >> this->paramDV[iDV][iParamDV];

        i++;
      }
      if (iDV < (this->nDV-1)) {
        if (option_value[i].compare(";") != 0) {
          string newstring;
          newstring.append(this->name);
          newstring.append(": a design variable in the configuration file has the wrong number of parameters");
          return newstring;
        }
        i++;
      }
    }

    // Need to return something...
    return "";
  }

  void SetDefault() {
    this->nDV = 0;
    this->paramDV = NULL;
    this->FFDTag = NULL;
    // Don't mess with the Design_Variable because it's an input, not modified
  }
};

class COptionDVValue : public COptionBase {
  string name; // identifier for the option
  unsigned short* & nDV_Value;
  su2double ** & valueDV;
  unsigned short & nDV;
  su2double ** & paramDV;
  unsigned short* & design_variable;

public:
  COptionDVValue(string option_field_name, unsigned short* & nDVValue_field, su2double** & valueDV_field, unsigned short & nDV_field,  su2double** & paramDV_field, unsigned short * & design_variable_field) : nDV_Value(nDVValue_field), valueDV(valueDV_field), nDV(nDV_field), paramDV(paramDV_field), design_variable(design_variable_field) {
    this->name = option_field_name;
  }

  ~COptionDVValue() {};

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    if ((option_value.size() == 1) && (option_value[0].compare("NONE") == 0)) {
      this->nDV_Value = NULL;
      return "";
    }

    if ( (this->nDV > 0) && (this->design_variable == NULL) ) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": Design_Variable array has not been allocated. Check that DV_KIND appears before DV_VALUE in configuration file.");
      return newstring;
    }
    if ( (this->nDV > 0) && (this->paramDV == NULL) ) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": Design_Parameter array has not been allocated. Check that DV_PARAM appears before DV_VALUE in configuration file.");
      return newstring;
    }

    this->valueDV = new su2double*[this->nDV];
    this->nDV_Value = new unsigned short[this->nDV];

    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) {
      this->valueDV[iDV] = new su2double[3];
    }

    unsigned short nValueDV = 0;
    unsigned short totalnValueDV = 0;
    stringstream ss;
    unsigned int i = 0;
    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) {
      switch (this->design_variable[iDV]) {
        case FFD_CONTROL_POINT:
          if((this->paramDV[iDV][4] == 0) &&
             (this->paramDV[iDV][5] == 0) &&
             (this->paramDV[iDV][6] == 0)) {
            nValueDV = 3;
          } else {
            nValueDV = 1;
          }
          break;
        case FFD_CONTROL_POINT_2D:
          if((this->paramDV[iDV][3] == 0) &&
             (this->paramDV[iDV][4] == 0)) {
            nValueDV = 2;
          } else {
            nValueDV = 1;
          }
          break;
        default :
          nValueDV = 1;
      }

      this->nDV_Value[iDV] = nValueDV;

      totalnValueDV += nValueDV;

      for (unsigned short iValueDV = 0; iValueDV < nValueDV; iValueDV++) {

        if (i >= option_value.size()) {
          string newstring;
          newstring.append(this->name);
          newstring.append(": DV_VALUE does not contain enough entries to match DV_KIND or DV_PARAM.");
          return newstring;
        }
        
        ss << option_value[i] << " ";

        ss >> this->valueDV[iDV][iValueDV];

        i++;
      }
    }

    if (i != totalnValueDV) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": a design variable in the configuration file has the wrong number of values");
      return newstring;
    }

    // Need to return something...
    return "";
  }

  void SetDefault() {
    this->nDV_Value = 0;
    this->valueDV = NULL;
    // Don't mess with the Design_Variable because it's an input, not modified
  }
};

class COptionFFDDef : public COptionBase {
  string name;
  unsigned short & nFFD;
  su2double ** & CoordFFD;
  string * & FFDTag;
  
public:
  COptionFFDDef(string option_field_name, unsigned short & nFFD_field, su2double** & coordFFD_field, string* & FFDTag_field) : nFFD(nFFD_field), CoordFFD(coordFFD_field), FFDTag(FFDTag_field) {
    this->name = option_field_name;
  }
  
  ~COptionFFDDef() {};
  
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    if ((option_value.size() == 1) && (option_value[0].compare("NONE") == 0)) {
      this->nFFD = 0;
      return "";
    }
    
    // Cannot have ; at the beginning or the end
    if (option_value[0].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have beginning semicolon");
      return newstring;
    }
    if (option_value[option_value.size()-1].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have ending semicolon");
      return newstring;
    }
    
    
    // use the ";" token to determine the number of design variables
    // This works because semicolon is not one of the delimiters in tokenize string
    this->nFFD = 0;
    for (unsigned int i = 0; i < static_cast<unsigned int>(option_value.size()); i++) {
      if (option_value[i].compare(";") == 0) {
        this->nFFD++;
      }
    }
    
    // One more design variable than semicolon
    this->nFFD++;
    
    this->CoordFFD = new su2double*[this->nFFD];
    for (unsigned short iFFD = 0; iFFD < this->nFFD; iFFD++) {
      this->CoordFFD[iFFD] = new su2double[25];
    }
    
    this->FFDTag = new string[this->nFFD];
    
    unsigned short nCoordFFD = 0;
    stringstream ss;
    unsigned int i = 0;
    
    for (unsigned short iFFD = 0; iFFD < this->nFFD; iFFD++) {
      
      nCoordFFD = 25;
      
      for (unsigned short iCoordFFD = 0; iCoordFFD < nCoordFFD; iCoordFFD++) {
        
        ss << option_value[i] << " ";
        
        if (iCoordFFD == 0) ss >> this->FFDTag[iFFD];
        else ss >> this->CoordFFD[iFFD][iCoordFFD-1];
        
        i++;
      }
      
      if (iFFD < (this->nFFD-1)) {
        if (option_value[i].compare(";") != 0) {
          string newstring;
          newstring.append(this->name);
          newstring.append(": a FFD box in the configuration file has the wrong number of parameters");
          return newstring;
        }
        i++;
      }
      
    }
    
    // Need to return something...
    return "";
  }
  
  void SetDefault() {
    this->nFFD = 0;
    this->CoordFFD = NULL;
    this->FFDTag = NULL;
  }
  
};

class COptionFFDDegree : public COptionBase {
  string name;
  unsigned short & nFFD;
  unsigned short ** & DegreeFFD;
  
public:
  COptionFFDDegree(string option_field_name, unsigned short & nFFD_field, unsigned short** & degreeFFD_field) : nFFD(nFFD_field), DegreeFFD(degreeFFD_field) {
    this->name = option_field_name;
  }
  
  ~COptionFFDDegree() {};
  
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    if ((option_value.size() == 1) && (option_value[0].compare("NONE") == 0)) {
      this->nFFD = 0;
      return "";
    }
    
    // Cannot have ; at the beginning or the end
    if (option_value[0].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have beginning semicolon");
      return newstring;
    }
    if (option_value[option_value.size()-1].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have ending semicolon");
      return newstring;
    }
    
    
    // use the ";" token to determine the number of design variables
    // This works because semicolon is not one of the delimiters in tokenize string
    this->nFFD = 0;
    for (unsigned int i = 0; i < static_cast<unsigned int>(option_value.size()); i++) {
      if (option_value[i].compare(";") == 0) {
        this->nFFD++;
      }
    }
    
    // One more design variable than semicolon
    this->nFFD++;
    
    this->DegreeFFD = new unsigned short*[this->nFFD];
    for (unsigned short iFFD = 0; iFFD < this->nFFD; iFFD++) {
      this->DegreeFFD[iFFD] = new unsigned short[3];
    }
    
    unsigned short nDegreeFFD = 0;
    stringstream ss;
    unsigned int i = 0;
    
    for (unsigned short iFFD = 0; iFFD < this->nFFD; iFFD++) {
      
      nDegreeFFD = 3;
      
      for (unsigned short iDegreeFFD = 0; iDegreeFFD < nDegreeFFD; iDegreeFFD++) {
        ss << option_value[i] << " ";
        ss >> this->DegreeFFD[iFFD][iDegreeFFD];
        i++;
      }
      
      if (iFFD < (this->nFFD-1)) {
        if (option_value[i].compare(";") != 0) {
          string newstring;
          newstring.append(this->name);
          newstring.append(": a FFD degree in the configuration file has the wrong number of parameters");
          return newstring;
        }
        i++;
      }
      
    }
    
    // Need to return something...
    return "";
  }
  
  void SetDefault() {
    this->nFFD = 0;
    this->DegreeFFD = NULL;
  }
  
};

// Class where the option is represented by (String, su2double, string, su2double, ...)
class COptionStringDoubleList : public COptionBase {
  string name; // identifier for the option
  unsigned short & size; // how many strings are there (same as number of su2doubles)

  string * & s_f; // Reference to the string fields
  su2double* & d_f; // reference to the su2double fields

public:
  COptionStringDoubleList(string option_field_name, unsigned short & list_size, string * & string_field, su2double* & double_field) : size(list_size), s_f(string_field), d_f(double_field) {
    this->name = option_field_name;
  }

  ~COptionStringDoubleList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // There must be an even number of entries (same number of strings and doubles
    unsigned short totalVals = option_value.size();
    if ((totalVals % 2) != 0) {
      if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
        // It's okay to say its NONE
        this->size = 0;
        return "";
      }
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have an even number of entries");
      return newstring;
    }
    unsigned short nVals = totalVals / 2;
    this->size = nVals;
    this->s_f = new string[nVals];
    this->d_f = new su2double[nVals];

    for (unsigned long i = 0; i < nVals; i++) {
      this->s_f[i].assign(option_value[2*i]); // 2 because have su2double and string
      istringstream is(option_value[2*i + 1]);
      su2double val;
      if (!(is >> val)) {
        return badValue(option_value, "string su2double", this->name);
      }
      this->d_f[i] = val;
    }
    // Need to return something...
    return "";
  }

  void SetDefault() {
    this->size = 0; // There is no default value for list
  }
};

class COptionInlet : public COptionBase {
  string name; // identifier for the option
  unsigned short & size;
  string * & marker;
  su2double * & ttotal;
  su2double * & ptotal;
  su2double ** & flowdir;

public:
  COptionInlet(string option_field_name, unsigned short & nMarker_Inlet, string* & Marker_Inlet, su2double* & Ttotal, su2double* & Ptotal, su2double** & FlowDir) : size(nMarker_Inlet), marker(Marker_Inlet), ttotal(Ttotal), ptotal(Ptotal), flowdir(FlowDir) {
    this->name = option_field_name;
  }

  ~COptionInlet() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    unsigned short totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker = NULL;
      this->ttotal = NULL;
      this->ptotal = NULL;
      this->flowdir = NULL;
      return "";
    }

    if (totalVals % 6 != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 6");
      this->size = 0;
      this->marker = NULL;
      this->ttotal = NULL;
      this->ptotal = NULL;
      this->flowdir = NULL;
      return newstring;
    }

    unsigned short nVals = totalVals / 6;
    this->size = nVals;
    this->marker = new string[nVals];
    this->ttotal = new su2double[nVals];
    this->ptotal = new su2double[nVals];
    this->flowdir = new su2double*[nVals];
    for (unsigned long i = 0; i < nVals; i++) {
      this->flowdir[i] = new su2double[3];
    }

    for (unsigned long i = 0; i < nVals; i++) {
      this->marker[i].assign(option_value[6*i]);
      istringstream ss_1st(option_value[6*i + 1]);
      if (!(ss_1st >> this->ttotal[i])) {
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_2nd(option_value[6*i + 2]);
      if (!(ss_2nd >> this->ptotal[i])) {
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_3rd(option_value[6*i + 3]);
      if (!(ss_3rd >> this->flowdir[i][0])) {
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_4th(option_value[6*i + 4]);
      if (!(ss_4th >> this->flowdir[i][1])) {
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_5th(option_value[6*i + 5]);
      if (!(ss_5th >> this->flowdir[i][2])) {
        return badValue(option_value, "inlet", this->name);
      }
    }

    return "";
  }

  void SetDefault() {
    this->marker = NULL;
    this->ttotal = NULL;
    this->ptotal = NULL;
    this->flowdir = NULL;
    this->size = 0; // There is no default value for list
  }
};

template <class Tenum>
class COptionRiemann : public COptionBase {

protected:
  map<string, Tenum> m;
  string name; // identifier for the option
  unsigned short & size;
  string * & marker;
  unsigned short* & field; // Reference to the field name
  su2double * & var1;
  su2double * & var2;
  su2double ** & flowdir;

public:
  COptionRiemann(string option_field_name, unsigned short & nMarker_Riemann, string* & Marker_Riemann, unsigned short* & option_field, const map<string, Tenum> m, su2double* & var1, su2double* & var2, su2double** & FlowDir) : size(nMarker_Riemann),
  	  	  	  	  marker(Marker_Riemann), field(option_field), var1(var1), var2(var2), flowdir(FlowDir) {
    this->name = option_field_name;
    this->m = m;
  }
  ~COptionRiemann() {};

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    unsigned short totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker = NULL;
      this->field = 0;
      this->var1 = NULL;
      this->var2 = NULL;
      this->flowdir = NULL;
      return "";
    }

    if (totalVals % 7 != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 7");
      this->size = 0;
      this->marker = NULL;
      this->var1 = NULL;
      this->var2 = NULL;
      this->flowdir = NULL;
      this->field = NULL;
      return newstring;
    }

    unsigned short nVals = totalVals / 7;
    this->size = nVals;
    this->marker = new string[nVals];
    this->var1 = new su2double[nVals];
    this->var2 = new su2double[nVals];
    this->flowdir = new su2double*[nVals];
    this->field = new unsigned short[nVals];

    for (unsigned long i = 0; i < nVals; i++) {
      this->flowdir[i] = new su2double[3];
    }

    for (unsigned long i = 0; i < nVals; i++) {
      this->marker[i].assign(option_value[7*i]);
        // Check to see if the enum value is in the map
    if (this->m.find(option_value[7*i + 1]) == m.end()) {
      string str;
      str.append(this->name);
      str.append(": invalid option value ");
      str.append(option_value[0]);
      str.append(". Check current SU2 options in config_template.cfg.");
      return str;
    }
      Tenum val = this->m[option_value[7*i + 1]];
      this->field[i] = val;

      istringstream ss_1st(option_value[7*i + 2]);
      if (!(ss_1st >> this->var1[i])) {
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_2nd(option_value[7*i + 3]);
      if (!(ss_2nd >> this->var2[i])) {
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_3rd(option_value[7*i + 4]);
      if (!(ss_3rd >> this->flowdir[i][0])) {
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_4th(option_value[7*i + 5]);
      if (!(ss_4th >> this->flowdir[i][1])) {
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_5th(option_value[7*i + 6]);
      if (!(ss_5th >> this->flowdir[i][2])) {
        return badValue(option_value, "Riemann", this->name);
      }
    }

    return "";
  }

  void SetDefault() {
    this->marker = NULL;
    this->var1 = NULL;
    this->var2 = NULL;
    this->flowdir = NULL;
    this->size = 0; // There is no default value for list
  }
};

template <class Tenum>
class COptionGiles : public COptionBase{

  map<string, Tenum> m;
  unsigned short & size;
  string * & marker;
  unsigned short* & field; // Reference to the fieldname
  string name; // identifier for the option
  su2double * & var1;
  su2double * & var2;
  su2double ** & flowdir;
  su2double * & relfac1;
  su2double * & relfac2;

public:
  COptionGiles(string option_field_name, unsigned short & nMarker_Giles, string* & Marker_Giles, unsigned short* & option_field, const map<string, Tenum> m, su2double* & var1, su2double* & var2, su2double** & FlowDir, su2double* & relfac1, su2double* & relfac2) : size(nMarker_Giles),
  	  	  	  	  marker(Marker_Giles), field(option_field), var1(var1), var2(var2), flowdir(FlowDir), relfac1(relfac1), relfac2(relfac2) {
    this->name = option_field_name;
    this->m = m;
  }
  ~COptionGiles() {};

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    unsigned long totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker = NULL;
      this->field = 0;
      this->var1 = NULL;
      this->var2 = NULL;
      this->flowdir = NULL;
      this->relfac1 = NULL;
      this->relfac2 = NULL;
      return "";
    }

    if (totalVals % 9 != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 9");
      this->size = 0;
      this->marker = NULL;
      this->var1 = NULL;
      this->var2 = NULL;
      this->flowdir = NULL;
      this->field = NULL;
      this->relfac1 = NULL;
      this->relfac2 = NULL;
      return newstring;
    }

    unsigned long nVals = totalVals / 9;
    this->size = nVals;
    this->marker = new string[nVals];
    this->var1 = new su2double[nVals];
    this->var2 = new su2double[nVals];
    this->flowdir = new su2double*[nVals];
    this->field = new unsigned short[nVals];
    this->relfac1 = new su2double[nVals];
    this->relfac2 = new su2double[nVals];

    for (unsigned int i = 0; i < nVals; i++) {
      this->flowdir[i] = new su2double[3];
    }

    for (unsigned int i = 0; i < nVals; i++) {
      this->marker[i].assign(option_value[9*i]);
        // Check to see if the enum value is in the map
    if (this->m.find(option_value[9*i + 1]) == m.end()) {
      string str;
      str.append(this->name);
      str.append(": invalid option value ");
      str.append(option_value[0]);
      str.append(". Check current SU2 options in config_template.cfg.");
      return str;
    }
      Tenum val = this->m[option_value[9*i + 1]];
      this->field[i] = val;

      istringstream ss_1st(option_value[9*i + 2]);
      if (!(ss_1st >> this->var1[i])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_2nd(option_value[9*i + 3]);
      if (!(ss_2nd >> this->var2[i])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_3rd(option_value[9*i + 4]);
      if (!(ss_3rd >> this->flowdir[i][0])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_4th(option_value[9*i + 5]);
      if (!(ss_4th >> this->flowdir[i][1])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_5th(option_value[9*i + 6]);
      if (!(ss_5th >> this->flowdir[i][2])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_6th(option_value[9*i + 7]);
      if (!(ss_6th >> this->relfac1[i])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_7th(option_value[9*i + 8]);
      if (!(ss_7th >> this->relfac2[i])) {
        return badValue(option_value, "Giles BC", this->name);
      }
    }

    return "";
  }

  void SetDefault() {
    this->marker = NULL;
    this->var1 = NULL;
    this->var2 = NULL;
    this->relfac1 = NULL;
    this->relfac2 = NULL;
    this->flowdir = NULL;
    this->size = 0; // There is no default value for list
  }
};

//Inlet condition where the input direction is assumed
class COptionExhaust : public COptionBase {
  string name; // identifier for the option
  unsigned short & size;
  string * & marker;
  su2double * & ttotal;
  su2double * & ptotal;

public:
  COptionExhaust(string option_field_name, unsigned short & nMarker_Exhaust, string* & Marker_Exhaust, su2double* & Ttotal, su2double* & Ptotal) : size(nMarker_Exhaust), marker(Marker_Exhaust), ttotal(Ttotal), ptotal(Ptotal) {
    this->name = option_field_name;
  }

  ~COptionExhaust() {};
  
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    unsigned short totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker = NULL;
      this->ttotal = NULL;
      this->ptotal = NULL;
      return "";
    }

    if (totalVals % 3 != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 3");
      this->size = 0;
      this->marker = NULL;
      this->ttotal = NULL;
      this->ptotal = NULL;
      return newstring;
    }

    unsigned short nVals = totalVals / 3;
    this->size = nVals;
    this->marker = new string[nVals];
    this->ttotal = new su2double[nVals];
    this->ptotal = new su2double[nVals];

    for (unsigned long i = 0; i < nVals; i++) {
      this->marker[i].assign(option_value[3*i]);
      istringstream ss_1st(option_value[3*i + 1]);
      if (!(ss_1st >> this->ttotal[i]))
        return badValue(option_value, "exhaust fixed", this->name);
      istringstream ss_2nd(option_value[3*i + 2]);
      if (!(ss_2nd >> this->ptotal[i]))
        return badValue(option_value, "exhaust fixed", this->name);
    }
    
    return "";
  }

  void SetDefault() {
    this->marker = NULL;
    this->ttotal = NULL;
    this->ptotal = NULL;
    this->size = 0; // There is no default value for list
  }
  
};

class COptionPeriodic : public COptionBase {
  string name; // identifier for the option
  unsigned short & size;
  string * & marker_bound;
  string * & marker_donor;
  su2double ** & rot_center;
  su2double ** & rot_angles;
  su2double ** & translation;

public:
  COptionPeriodic(const string option_field_name, unsigned short & nMarker_PerBound,
                  string* & Marker_PerBound, string* & Marker_PerDonor,
                  su2double** & RotCenter, su2double** & RotAngles, su2double** & Translation) : size(nMarker_PerBound), marker_bound(Marker_PerBound), marker_donor(Marker_PerDonor), rot_center(RotCenter), rot_angles(RotAngles), translation(Translation) {
    this->name = option_field_name;
  }

  ~COptionPeriodic() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    const int mod_num = 11;

    unsigned short totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker_bound = NULL;
      this->marker_donor = NULL;
      this->rot_center = NULL;
      this->rot_angles = NULL;
      this->translation = NULL;
      return "";
    }

    if (totalVals % mod_num != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 11");
      this->size = 0;
      this->marker_bound = NULL;
      this->marker_donor = NULL;
      this->rot_center = NULL;
      this->rot_angles = NULL;
      this->translation = NULL;
      return newstring;
    }

    unsigned short nVals = 2 * (totalVals / mod_num); // To account for periodic and donor
    this->size = nVals;
    this->marker_bound = new string[nVals];
    this->marker_donor = new string[nVals];
    this->rot_center = new su2double*[nVals];
    this->rot_angles = new su2double*[nVals];
    this->translation = new su2double*[nVals];
    for (unsigned long i = 0; i < nVals; i++) {
      this->rot_center[i] = new su2double[3];
      this->rot_angles[i] = new su2double[3];
      this->translation[i] = new su2double[3];
    }

    su2double deg2rad = PI_NUMBER/180.0;

    for (unsigned long i = 0; i < (nVals/2); i++) {
      this->marker_bound[i].assign(option_value[mod_num*i]);
      this->marker_donor[i].assign(option_value[mod_num*i+1]);
      istringstream ss_1st(option_value[mod_num*i + 2]);
      if (!(ss_1st >> this->rot_center[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_2nd(option_value[mod_num*i + 3]);
      if (!(ss_2nd >> this->rot_center[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_3rd(option_value[mod_num*i + 4]);
      if (!(ss_3rd >> this->rot_center[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_4th(option_value[mod_num*i + 5]);
      if (!(ss_4th >> this->rot_angles[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_5th(option_value[mod_num*i + 6]);
      if (!(ss_5th >> this->rot_angles[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_6th(option_value[mod_num*i + 7]);
      if (!(ss_6th >> this->rot_angles[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_7th(option_value[mod_num*i + 8]);
      if (!(ss_7th >> this->translation[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_8th(option_value[mod_num*i + 9]);
      if (!(ss_8th >> this->translation[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_9th(option_value[mod_num*i + 10]);
      if (!(ss_9th >> this->translation[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      this->rot_angles[i][0] *= deg2rad;
      this->rot_angles[i][1] *= deg2rad;
      this->rot_angles[i][2] *= deg2rad;
    }

    for (unsigned long i = (nVals/2); i < nVals; i++) {
      this->marker_bound[i].assign(option_value[mod_num*(i-nVals/2)+1]);
      this->marker_donor[i].assign(option_value[mod_num*(i-nVals/2)]);
      istringstream ss_1st(option_value[mod_num*(i-nVals/2) + 2]);
      if (!(ss_1st >> this->rot_center[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_2nd(option_value[mod_num*(i-nVals/2) + 3]);
      if (!(ss_2nd >> this->rot_center[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_3rd(option_value[mod_num*(i-nVals/2) + 4]);
      if (!(ss_3rd >> this->rot_center[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_4th(option_value[mod_num*(i-nVals/2) + 5]);
      if (!(ss_4th >> this->rot_angles[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_5th(option_value[mod_num*(i-nVals/2) + 6]);
      if (!(ss_5th >> this->rot_angles[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_6th(option_value[mod_num*(i-nVals/2) + 7]);
      if (!(ss_6th >> this->rot_angles[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_7th(option_value[mod_num*(i-nVals/2) + 8]);
      if (!(ss_7th >> this->translation[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_8th(option_value[mod_num*(i-nVals/2) + 9]);
      if (!(ss_8th >> this->translation[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_9th(option_value[mod_num*(i-nVals/2) + 10]);
      if (!(ss_9th >> this->translation[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      /*--- Mirror the rotational angles and translation vector (rotational
       center does not need to move) ---*/
      this->rot_center[i][0] *= 1.0;
      this->rot_center[i][1] *= 1.0;
      this->rot_center[i][2] *= 1.0;
      this->rot_angles[i][0] *= -deg2rad;
      this->rot_angles[i][1] *= -deg2rad;
      this->rot_angles[i][2] *= -deg2rad;
      this->translation[i][0] *= -1.0;
      this->translation[i][1] *= -1.0;
      this->translation[i][2] *= -1.0;
    }

    return "";
  }

  void SetDefault() {
    this->size = 0;
    this->marker_bound = NULL;
    this->marker_donor = NULL;
    this->rot_center = NULL;
    this->rot_angles = NULL;
    this->translation = NULL;
  }
};

class COptionTurboPerformance : public COptionBase {
  string name; // identifier for the option
  unsigned short & size;
  string * & marker_turboIn;
  string * & marker_turboOut;

public:
  COptionTurboPerformance(const string option_field_name, unsigned short & nMarker_TurboPerf,
                          string* & Marker_TurboBoundIn, string* & Marker_TurboBoundOut) : size(nMarker_TurboPerf), marker_turboIn(Marker_TurboBoundIn), marker_turboOut(Marker_TurboBoundOut){
    this->name = option_field_name;
  }

  ~COptionTurboPerformance() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    const int mod_num = 2;

    unsigned long totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker_turboIn= NULL;
      this->marker_turboOut = NULL;
      return "";
    }

    if (totalVals % mod_num != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 2");
      this->size = 0;
      this->marker_turboIn= NULL;
      this->marker_turboOut = NULL;;
      return newstring;
    }

    unsigned long nVals = totalVals / mod_num;
    this->size = nVals;
    this->marker_turboIn = new string[nVals];
    this->marker_turboOut = new string[nVals];
    for (unsigned long i = 0; i < nVals; i++) {
      this->marker_turboIn[i].assign(option_value[mod_num*i]);
      this->marker_turboOut[i].assign(option_value[mod_num*i+1]);
     }


    return "";
  }

  void SetDefault() {
    this->size = 0;
    this->marker_turboIn= NULL;
    this->marker_turboOut = NULL;
  }
};

class COptionPython : public COptionBase {
  string name;
public:
  COptionPython(const string name) {
    this->name = name;
  }
  ~COptionPython() {};
  // No checking happens with python options
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    return "";
  }
  // No defaults with python options
  void SetDefault() {
    return;
  };
};

class COptionActDisk : public COptionBase {
  string name; // identifier for the option
  unsigned short & inlet_size;
  unsigned short & outlet_size;
  string * & marker_inlet;
  string * & marker_outlet;
  su2double ** & press_jump;
  su2double ** & temp_jump;
  su2double ** & omega;
  
public:
  COptionActDisk(const string name,
                 unsigned short & nMarker_ActDiskInlet, unsigned short & nMarker_ActDiskOutlet, string * & Marker_ActDiskInlet, string * & Marker_ActDiskOutlet,
                 su2double ** & ActDisk_PressJump, su2double ** & ActDisk_TempJump, su2double ** & ActDisk_Omega) :
  inlet_size(nMarker_ActDiskInlet), outlet_size(nMarker_ActDiskOutlet), marker_inlet(Marker_ActDiskInlet), marker_outlet(Marker_ActDiskOutlet),
  press_jump(ActDisk_PressJump), temp_jump(ActDisk_TempJump), omega(ActDisk_Omega) {
    this->name = name;
  }
  
  ~COptionActDisk() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    const int mod_num = 8;
    unsigned short totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->SetDefault();
      return "";
    }
    
    if (totalVals % mod_num != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 8");
      this->SetDefault();
      return newstring;
    }
    
    unsigned short nVals = totalVals / mod_num;
    this->inlet_size = nVals;
    this->outlet_size = nVals;
    this->marker_inlet = new string[this->inlet_size];
    this->marker_outlet = new string[this->outlet_size];
    
    this->press_jump = new su2double*[this->inlet_size];
    this->temp_jump = new su2double*[this->inlet_size];
    this->omega = new su2double*[this->inlet_size];
    for (int i = 0; i < this->inlet_size; i++) {
      this->press_jump[i] = new su2double[2];
      this->temp_jump[i] = new su2double[2];
      this->omega[i] = new su2double[2];
    }
    
    string tname = "actuator disk";
    
    for (int i = 0; i < this->inlet_size; i++) {
      this->marker_inlet[i].assign(option_value[mod_num*i]);
      this->marker_outlet[i].assign(option_value[mod_num*i+1]);
      istringstream ss_1st(option_value[mod_num*i + 2]);
      if (!(ss_1st >> this->press_jump[i][0])) {
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_2nd(option_value[mod_num*i + 3]);
      if (!(ss_2nd >> this->temp_jump[i][0])) {
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_3rd(option_value[mod_num*i + 4]);
      if (!(ss_3rd >> this->omega[i][0])) {
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_4th(option_value[mod_num*i + 5]);
      if (!(ss_4th >> this->press_jump[i][1])) {
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_5th(option_value[mod_num*i + 6]);
      if (!(ss_5th >> this->temp_jump[i][1])) {
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_6th(option_value[mod_num*i + 7]);
      if (!(ss_6th >> this->omega[i][1])) {
        return badValue(option_value, tname, this->name);
      }
    }
    return "";
  }
  void SetDefault() {
    this->inlet_size = 0;
    this->outlet_size = 0;
    this->marker_inlet = NULL;
    this->marker_outlet = NULL;
    this->press_jump = NULL;
    this->temp_jump = NULL;
    this->omega = NULL;
  }
};

class COptionWallFunction : public COptionBase {
  string name; // identifier for the option
  unsigned short &nMarkers;
  string* &markers;
  unsigned short*  &walltype;
  unsigned short** &intInfo;
  su2double**      &doubleInfo;

public:
  COptionWallFunction(const string name, unsigned short &nMarker_WF, 
                      string* &Marker_WF, unsigned short* &type_WF,
                      unsigned short** &intInfo_WF, su2double** &doubleInfo_WF) :
  nMarkers(nMarker_WF), markers(Marker_WF), walltype(type_WF),
  intInfo(intInfo_WF), doubleInfo(doubleInfo_WF) {
    this->name = name;
  }

  ~COptionWallFunction(){}

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    /*--- First check if NONE is specified. ---*/
    unsigned short totalSize = option_value.size();
    if ((totalSize == 1) && (option_value[0].compare("NONE") == 0)) {
      this->SetDefault();
      return "";
    }

    /*--- Determine the number of markers, for which a wall
          function treatment has been specified. ---*/
    unsigned short counter = 0, nVals = 0;
    while (counter < totalSize ) {

      /* Update the counter for the number of markers specified
         and store the current index for possible error messages. */
      ++nVals;
      const unsigned short indMarker = counter;

      /* Check if a wall function type has been specified for this marker.
         If not, create an error message and return. */
      ++counter;
      const unsigned short indWallType = counter;
      unsigned short typeWF = NO_WALL_FUNCTION;
      bool validWF = true;
      if (counter == totalSize) validWF = false;
      else {
        map<string, ENUM_WALL_FUNCTIONS>::const_iterator it;
        it = Wall_Functions_Map.find(option_value[counter]);
        if(it == Wall_Functions_Map.end()) validWF = false;
        else                               typeWF  = it->second;
      }

      if (!validWF ) {
        string newstring;
        newstring.append(this->name);
        newstring.append(": Invalid wall function type, ");
        newstring.append(option_value[counter]);
        newstring.append(", encountered for marker ");
        newstring.append(option_value[indMarker]);
        return newstring;
      }

      /* Update the counter, as the wall function type is valid. */
      ++counter;

      /*--- For some wall function types some additional info
            must be specified. Hence the counter must be updated
            accordingly. ---*/
      switch( typeWF ) {
        case EQUILIBRIUM_WALL_MODEL:    counter += 3; break;
        case NONEQUILIBRIUM_WALL_MODEL: counter += 2; break;
        case LOGARITHMIC_WALL_MODEL: counter += 3; break;
        default: break;
      }

      /* In case the counter is larger than totalSize, the data for
         this wall function type has not been specified correctly. */
      if (counter > totalSize) {
        string newstring;
        newstring.append(this->name);
        newstring.append(", marker ");
        newstring.append(option_value[indMarker]);
        newstring.append(", wall function type ");
        newstring.append(option_value[indWallType]);
        newstring.append(": Additional information is missing.");
        return newstring;
      }
    }

    /* Allocate the memory to store the data for the wall function markers. */
    this->nMarkers   = nVals;
    this->markers    = new string[nVals];
    this->walltype   = new unsigned short[nVals];
    this->intInfo    = new unsigned short*[nVals];
    this->doubleInfo = new su2double*[nVals];

    for (unsigned short i=0; i<nVals; i++) {
      this->intInfo[i]    = NULL;
      this->doubleInfo[i] = NULL;
    }

    /*--- Loop over the wall markers and store the info in the
          appropriate arrays. ---*/
    counter = 0;
    for (unsigned short i=0; i<nVals; i++) {

      /* Set the name of the wall function marker. */
      this->markers[i].assign(option_value[counter++]);

      /* Determine the wall function type. As their validaties have
         already been tested, there is no need to do so again. */
      map<string, ENUM_WALL_FUNCTIONS>::const_iterator it;
      it = Wall_Functions_Map.find(option_value[counter++]);

      this->walltype[i] = it->second;

      /*--- For some wall function types, some additional info
            is needed, which is extracted from option_value. ---*/
      switch( this->walltype[i] ) {

        case EQUILIBRIUM_WALL_MODEL: {

          /* LES equilibrium wall model. The exchange distance, stretching
             factor and number of points in the wall model must be specified. */
          this->intInfo[i]    = new unsigned short[1];
          this->doubleInfo[i] = new su2double[2];

          istringstream ss_1st(option_value[counter++]);
          if (!(ss_1st >> this->doubleInfo[i][0])) {
            return badValue(option_value, "su2double", this->name);
          }

          istringstream ss_2nd(option_value[counter++]);
          if (!(ss_2nd >> this->doubleInfo[i][1])) {
            return badValue(option_value, "su2double", this->name);
          }

          istringstream ss_3rd(option_value[counter++]);
          if (!(ss_3rd >> this->intInfo[i][0])) {
            return badValue(option_value, "unsigned short", this->name);
          }

          break;
        }

        case NONEQUILIBRIUM_WALL_MODEL: {

          /* LES non-equilibrium model. The RANS turbulence model and
             the exchange distance need to be specified. */
          this->intInfo[i]    = new unsigned short[1];
          this->doubleInfo[i] = new su2double[1];

          /* Check for a valid RANS turbulence model. */
          map<string, ENUM_TURB_MODEL>::const_iterator iit;
          iit = Turb_Model_Map.find(option_value[counter++]);
          if(iit == Turb_Model_Map.end()) {
            string newstring;
            newstring.append(this->name);
            newstring.append(", marker ");
            newstring.append(this->markers[i]);
            newstring.append(", wall function type ");
            newstring.append(option_value[counter-2]);
            newstring.append(": Invalid RANS turbulence model, ");
            newstring.append(option_value[counter-1]);
            newstring.append(", specified");
            return newstring;
          }

          this->intInfo[i][0] = iit->second;

          /* Extract the exchange distance. */
          istringstream ss_1st(option_value[counter++]);
          if (!(ss_1st >> this->doubleInfo[i][0])) {
            return badValue(option_value, "su2double", this->name);
          }

          break;
        }
        case LOGARITHMIC_WALL_MODEL: {
          
          /* LES Logarithmic law-of-the-wall model. The exchange distance, stretching
           factor and number of points in the wall model must be specified. */
          this->intInfo[i]    = new unsigned short[1];
          this->doubleInfo[i] = new su2double[2];
          
          istringstream ss_1st(option_value[counter++]);
          if (!(ss_1st >> this->doubleInfo[i][0])) {
            return badValue(option_value, "su2double", this->name);
          }
          
          istringstream ss_2nd(option_value[counter++]);
          if (!(ss_2nd >> this->doubleInfo[i][1])) {
            return badValue(option_value, "su2double", this->name);
          }
          
          istringstream ss_3rd(option_value[counter++]);
          if (!(ss_3rd >> this->intInfo[i][0])) {
            return badValue(option_value, "unsigned short", this->name);
          }
          
          break;
        }

        default: // Just to avoid a compiler warning.
          break;
      }
    }

    // Need to return something...
    return "";
  }

  void SetDefault() {
    this->nMarkers   = 0;
    this->markers    = NULL;
    this->walltype   = NULL;
    this->intInfo    = NULL;
    this->doubleInfo = NULL;
  }
};
