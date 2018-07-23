/*!
 * \file option_structure.hpp
 * \brief Defines classes for referencing options for easy input in CConfig
 * \author J. Hicken, B. Tracey
 * \version 5.0.0 "Raven"
 *
 * Many of the classes in this file are templated, and therefore must
 * be declared and defined here; to keep all elements together, there
 * is no corresponding .cpp file at this time.
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

const unsigned int EXIT_DIVERGENCE = 2; /*!< \brief Exit code (divergence). */

const unsigned int BUFSIZE = 3000000;		     /*!< \brief MPI buffer. */
const unsigned int MAX_PARAMETERS = 10;		   /*!< \brief Maximum number of parameters for a design variable definition. */
const unsigned int MAX_NUMBER_PERIODIC = 10; /*!< \brief Maximum number of periodic boundary conditions. */
const unsigned int MAX_STRING_SIZE = 200;    /*!< \brief Maximum number of domains. */
const unsigned int MAX_NUMBER_FFD = 10;	     /*!< \brief Maximum number of FFDBoxes for the FFD. */
const unsigned int MAX_SOLS = 6;		         /*!< \brief Maximum number of solutions at the same time (dimension of solution container array). */
const unsigned int MAX_TERMS = 6;		         /*!< \brief Maximum number of terms in the numerical equations (dimension of solver container array). */
const unsigned int MAX_ZONES = 3;            /*!< \brief Maximum number of zones. */
const unsigned int MAX_FE_KINDS = 4;            	/*!< \brief Maximum number of Finite Elements. */
const unsigned int NO_RK_ITER = 0;		       /*!< \brief No Runge-Kutta iteration. */

const unsigned int OVERHEAD = 4; /*!< \brief Overhead space above nMarker when allocating space for boundary elems (MPI + periodic). */

const unsigned int MESH_0 = 0; /*!< \brief Definition of the finest grid level. */
const unsigned int MESH_1 = 1; /*!< \brief Definition of the finest grid level. */
const unsigned int ZONE_0 = 0; /*!< \brief Definition of the first grid domain. */
const unsigned int ZONE_1 = 1; /*!< \brief Definition of the first grid domain. */

const su2double STANDART_GRAVITY = 9.80665;           /*!< \brief Acceleration due to gravity at surface of earth. */

const su2double EPS = 1.0E-16;		   /*!< \brief Error scale. */
const su2double TURB_EPS = 1.0E-16; /*!< \brief Turbulent Error scale. */

const su2double ONE2 = 0.5;			   /*!< \brief One divided by two. */
const su2double TWO3 = 2.0 / 3.0;	 /*!< \brief Two divided by three. */
const su2double FOUR3 = 4.0 / 3.0;  /*!< \brief Four divided by three. */

const su2double PI_NUMBER = 4.0 * atan(1.0);	/*!< \brief Pi number. */

const int MASTER_NODE = 0;			/*!< \brief Master node for MPI parallelization. */
const int SINGLE_NODE = 1;			/*!< \brief There is only a node in the MPI parallelization. */
const int SINGLE_ZONE = 1;			/*!< \brief There is only a zone. */

const unsigned short N_ELEM_TYPES = 7;           /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_LINE = 2;          /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_TRIANGLE = 3;      /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_QUADRILATERAL = 4; /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_TETRAHEDRON = 4;   /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_HEXAHEDRON = 8;    /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_PYRAMID = 5;       /*!< \brief General output & CGNS defines. */
const unsigned short N_POINTS_PRISM = 6;         /*!< \brief General output & CGNS defines. */

/*!
 * \brief Boolean answers
 */
enum ANSWER {
  NONE = 0,
  NO = 0,    /*!< \brief Boolean definition of no. */
  YES = 1	/*!< \brief Boolean definition of yes. */
};

/*!
 * \brief Verbosity level
 */
enum VERB_LEVEL {
  VERB_NONE = 0,   /*!< \brief No verbosity. */
  VERB_MEDIUM = 1,   /*!< \brief Medium level of verbosity. */
  VERB_HIGH = 2			/*!< \brief High level of verbosity. */
};
static const map<string, VERB_LEVEL> Verb_Map = CCreateMap<string, VERB_LEVEL>
("NONE", VERB_NONE)
("MEDIUM", VERB_MEDIUM)
("HIGH", VERB_HIGH);

/*!
 * \brief Average method for marker analyze
 */
enum AVERAGE_TYPE {
  AVERAGE_AREA = 1, /*!< \brief Area-weighted average. */
  AVERAGE_MASSFLUX = 2 /*!< \brief Mass-flux weighted average. */
};
static const map<string, AVERAGE_TYPE> Average_Map = CCreateMap<string, AVERAGE_TYPE>
("AREA", AVERAGE_AREA)
("MASSFLUX", AVERAGE_MASSFLUX);

/*!
 * \brief different solver types for the CFD component
 */
enum ENUM_SOLVER {
  NO_SOLVER = 0,						/*!< \brief Definition of no solver. */
  EULER = 1,							/*!< \brief Definition of the Euler's solver. */
  NAVIER_STOKES = 2,					/*!< \brief Definition of the Navier-Stokes' solver. */
  RANS = 3,								/*!< \brief Definition of the Reynolds-averaged Navier-Stokes' (RANS) solver. */
  POISSON_EQUATION = 4,       			/*!< \brief Definition of the poisson potential solver. */
  WAVE_EQUATION = 10,					/*!< \brief Definition of the wave solver. */
  HEAT_EQUATION = 29,					/*!< \brief Definition of the heat solver. */
  FLUID_STRUCTURE_INTERACTION = 12,		/*!< \brief Definition of a FSI solver. */
  FEM_ELASTICITY = 13,					/*!< \brief Definition of a FEM solver. */
  ADJ_EULER = 18,						/*!< \brief Definition of the continuous adjoint Euler's solver. */
  ADJ_NAVIER_STOKES = 19,				/*!< \brief Definition of the continuous adjoint Navier-Stokes' solver. */
  ADJ_RANS = 20,						/*!< \brief Definition of the continuous adjoint Reynolds-averaged Navier-Stokes' (RANS) solver. */
  TEMPLATE_SOLVER = 30,                 /*!< \brief Definition of template solver. */
  DISC_ADJ_EULER = 35,
  DISC_ADJ_RANS = 36,
  DISC_ADJ_NAVIER_STOKES = 37
};
/* BEGIN_CONFIG_ENUMS */
static const map<string, ENUM_SOLVER> Solver_Map = CCreateMap<string, ENUM_SOLVER>
("NONE", NO_SOLVER)
("EULER", EULER)
("NAVIER_STOKES", NAVIER_STOKES)
("RANS", RANS)
("POISSON_EQUATION", POISSON_EQUATION)
("ADJ_EULER", ADJ_EULER)
("ADJ_NAVIER_STOKES", ADJ_NAVIER_STOKES)
("ADJ_RANS", ADJ_RANS )
("WAVE_EQUATION", WAVE_EQUATION)
("HEAT_EQUATION", HEAT_EQUATION)
("FEM_ELASTICITY", FEM_ELASTICITY)
("DISC_ADJ_EULER", DISC_ADJ_EULER)
("DISC_ADJ_RANS", DISC_ADJ_RANS)
("DISC_ADJ_NAVIERSTOKES", DISC_ADJ_EULER)
("FLUID_STRUCTURE_INTERACTION", FLUID_STRUCTURE_INTERACTION)

("TEMPLATE_SOLVER", TEMPLATE_SOLVER);


/*!
 * \brief types of fluid solvers
 */
enum ENUM_FSI_FLUID_PROBLEM {
	  NO_SOLVER_FFSI = 0,			/*!< \brief Definition of no solver. */
	  EULER_FFSI = 1,				/*!< \brief Euler equations for the FSI problem */
	  NAVIER_STOKES_FFSI = 2,		/*!< \brief NS equations for the FSI problem */
	  RANS_FFSI = 3 				/*!< \brief RANS equations for the FSI problem */
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
  NO_SOLVER_SFSI = 0,				/*!< \brief Definition of no solver. */
  FEM_ELASTICITY_SFSI = 13,		/*!< \brief Nonlinear elasticity equations for the FSI problem */
};
static const map<string, ENUM_FSI_STRUC_PROBLEM> FSI_Struc_Solver_Map = CCreateMap<string, ENUM_FSI_STRUC_PROBLEM>
("NONE", NO_SOLVER_SFSI)
("FEM_ELASTICITY", FEM_ELASTICITY_SFSI);

/*!
 * \brief Material geometric conditions
 */
enum ENUM_STRUCT_SOLVER {
	SMALL_DEFORMATIONS = 0,			/*!< \brief Definition of linear elastic material. */
	LARGE_DEFORMATIONS = 1,			/*!< \brief Definition of Neo-Hookean material. */
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
};
static const map<string, ENUM_MATERIAL_MODEL> Material_Map = CCreateMap<string, ENUM_MATERIAL_MODEL>
("LINEAR_ELASTIC", LINEAR_ELASTIC)
("NEO_HOOKEAN", NEO_HOOKEAN);

/*!
 * \brief Material compressibility
 */
enum ENUM_MAT_COMPRESS {
  COMPRESSIBLE_MAT = 0,			/*!< \brief Definition of compressible material. */
  INCOMPRESSIBLE_MAT = 1,		/*!< \brief Definition of incompressible material. */
};
static const map<string, ENUM_MAT_COMPRESS> MatComp_Map = CCreateMap<string, ENUM_MAT_COMPRESS>
("COMPRESSIBLE", COMPRESSIBLE_MAT)
("INCOMPRESSIBLE", INCOMPRESSIBLE_MAT);



/*!
 * \brief types of interpolators
 */
enum ENUM_INTERPOLATOR {
  NEAREST_NEIGHBOR 	= 0,   	/*!< \brief Nearest Neigbhor interpolation */
  ISOPARAMETRIC 	= 1,	/*!< \brief Isoparametric interpolation */
  CONSISTCONSERVE 	= 2,	/*!< \brief Consistent & Conservative interpolation (S.A. Brown 1997). Utilizes Isoparametric interpolation. */
  WEIGHTED_AVERAGE  = 3, 	/*!< \brief Sliding Mesh Approach E. Rinaldi 2015 */
};

static const map<string, ENUM_INTERPOLATOR> Interpolator_Map = CCreateMap<string, ENUM_INTERPOLATOR>
("NEAREST_NEIGHBOR", NEAREST_NEIGHBOR)
("ISOPARAMETRIC",    ISOPARAMETRIC)
("CONSISTCONSERVE",  CONSISTCONSERVE)
("WEIGHTED_AVERAGE", WEIGHTED_AVERAGE);

/*!
 * \brief different regime modes
 */
enum ENUM_REGIME {
  COMPRESSIBLE = 0,			/*!< \brief Definition of compressible solver. */
  INCOMPRESSIBLE = 1,				/*!< \brief Definition of incompressible solver. */
};
static const map<string, ENUM_REGIME> Regime_Map = CCreateMap<string, ENUM_REGIME>
("COMPRESSIBLE", COMPRESSIBLE)
("INCOMPRESSIBLE", INCOMPRESSIBLE);
/*!
 * \brief different non-dimensional modes
 */
enum ENUM_KIND_NONDIM {
  DIMENSIONAL = 0,			    /*!< \brief Dimensional simulation. */
  FREESTREAM_PRESS_EQ_ONE = 1, /*!< \brief Non-dimensional simulation. */
  FREESTREAM_VEL_EQ_MACH = 2, /*!< \brief Non-dimensional simulation. */
  FREESTREAM_VEL_EQ_ONE = 3 /*!< \brief Non-dimensional simulation. */
};
static const map<string, ENUM_KIND_NONDIM> NonDim_Map = CCreateMap<string, ENUM_KIND_NONDIM>
("DIMENSIONAL", DIMENSIONAL)
("FREESTREAM_PRESS_EQ_ONE", FREESTREAM_PRESS_EQ_ONE)
("FREESTREAM_VEL_EQ_MACH", FREESTREAM_VEL_EQ_MACH)
("FREESTREAM_VEL_EQ_ONE", FREESTREAM_VEL_EQ_ONE);

/*!
 * \brief different system of measurements
 */
enum ENUM_MEASUREMENTS {
  SI = 0,			/*!< \brief Definition of compressible solver. */
  US = 1				/*!< \brief Definition of incompressible solver. */
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
  RUNTIME_POISSON_SYS = 4,			/*!< \brief One-physics case, the code is solving the poissonal potential equation. */
  RUNTIME_ADJPOT_SYS = 5,		/*!< \brief One-physics case, the code is solving the adjoint potential flow equation. */
  RUNTIME_ADJFLOW_SYS = 6,		/*!< \brief One-physics case, the code is solving the adjoint equations is being solved (Euler and Navier-Stokes). */
  RUNTIME_ADJTURB_SYS = 7,		/*!< \brief One-physics case, the code is solving the adjoint turbulence model. */
  RUNTIME_WAVE_SYS = 8,		/*!< \brief One-physics case, the code is solving the wave equation. */
  RUNTIME_MULTIGRID_SYS = 14,   	/*!< \brief Full Approximation Storage Multigrid system of equations. */
  RUNTIME_FEA_SYS = 20,		/*!< \brief One-physics case, the code is solving the FEA equation. */
  RUNTIME_HEAT_SYS = 21,		/*!< \brief One-physics case, the code is solving the heat equation. */
  RUNTIME_TRANS_SYS = 22,			/*!< \brief One-physics case, the code is solving the turbulence model. */
};

const int FLOW_SOL = 0;		/*!< \brief Position of the mean flow solution in the solver container array. */
const int ADJFLOW_SOL = 1;	/*!< \brief Position of the continuous adjoint flow solution in the solver container array. */

const int TURB_SOL = 2;		/*!< \brief Position of the turbulence model solution in the solver container array. */
const int ADJTURB_SOL = 3;	/*!< \brief Position of the continuous adjoint turbulence solution in the solver container array. */

const int TRANS_SOL = 4;	/*!< \brief Position of the transition model solution in the solver container array. */
const int POISSON_SOL = 2;		/*!< \brief Position of the electronic potential solution in the solver container array. */
const int WAVE_SOL = 1;		/*!< \brief Position of the wave equation in the solution solver array. */
const int HEAT_SOL = 2;		/*!< \brief Position of the heat equation in the solution solver array. */
const int FEA_SOL = 1;		/*!< \brief Position of the FEA equation in the solution solver array. */

const int TEMPLATE_SOL = 0;     /*!< \brief Position of the template solution. */

const int CONV_TERM = 0;	/*!< \brief Position of the convective terms in the numerics container array. */
const int VISC_TERM = 1;        /*!< \brief Position of the viscous terms in the numerics container array. */
const int SOURCE_FIRST_TERM = 2;        /*!< \brief Position of the first source term in the numerics container array. */
const int SOURCE_SECOND_TERM = 3;   /*!< \brief Position of the second source term in the numerics container array. */
const int CONV_BOUND_TERM = 4;       /*!< \brief Position of the convective boundary terms in the numerics container array. */
const int VISC_BOUND_TERM = 5;       /*!< \brief Position of the viscous boundary terms in the numerics container array. */

const int FEA_TERM = 0;			/*!< \brief Position of the finite element analysis terms in the numerics container array. */


/*!
 * \brief types of finite elements (in 2D or 3D)
 */

const int EL_TRIA = 0;		/*!< \brief Elements of three nodes (2D). */
const int EL_QUAD = 1;		/*!< \brief Elements of four nodes (2D). */

const int EL_TETRA = 0;		/*!< \brief Elements of four nodes (3D). */
const int EL_HEXA = 1;		/*!< \brief Elements of eight nodes (3D). */


/*!
 * \brief types of mathematical problem to solve
 */
enum ENUM_MATH_PROBLEM {
  DIRECT = 0,		/*!< \brief Direct problem */
  CONTINUOUS_ADJOINT = 1,		/*!< \brief Continuous adjoint problem */
  DISCRETE_ADJOINT = 2 /*< \brief AD-based discrete adjoint problem. */
};
static const map<string, ENUM_MATH_PROBLEM> Math_Problem_Map = CCreateMap<string, ENUM_MATH_PROBLEM>
("DIRECT", DIRECT)
("CONTINUOUS_ADJOINT", CONTINUOUS_ADJOINT)
("DISCRETE_ADJOINT", DISCRETE_ADJOINT);

/*!
 * \brief types of spatial discretizations
 */
enum ENUM_SPACE {
  NO_CONVECTIVE = 0, /*!< \brief No convective scheme is used. */
  SPACE_CENTERED = 1,		/*!< \brief Space centered convective numerical method. */
  SPACE_UPWIND = 2		/*!< \brief Upwind convective numerical method. */
};
static const map<string, ENUM_SPACE> Space_Map = CCreateMap<string, ENUM_SPACE>
("NONE", NO_CONVECTIVE)
("SPACE_CENTERED", SPACE_CENTERED)
("SPACE_UPWIND", SPACE_UPWIND);

/*!
 * \brief types of fluid model
 */
enum ENUM_FLUIDMODEL {
	STANDARD_AIR = 0,
	IDEAL_GAS = 1, /*!< \brief _____. */
	VW_GAS = 2,
	PR_GAS = 3,
	LUT =4
};

static const map<string, ENUM_FLUIDMODEL> FluidModel_Map = CCreateMap<string, ENUM_FLUIDMODEL>
("STANDARD_AIR", STANDARD_AIR)
("IDEAL_GAS", IDEAL_GAS)
("VW_GAS", VW_GAS)
("PR_GAS", PR_GAS)
("LUT", LUT);

/*!
 * \brief types of initialization option
 */

enum ENUM_INIT_OPTION {
	REYNOLDS = 0, /*!< \brief _____. */
	TD_CONDITIONS = 1

};

static const map<string, ENUM_INIT_OPTION> InitOption_Map = CCreateMap<string, ENUM_INIT_OPTION>
("REYNOLDS", REYNOLDS)
("TD_CONDITIONS", TD_CONDITIONS);

/*!
 * \brief types of initialization option
 */

enum ENUM_FREESTREAM_OPTION {
	TEMPERATURE_FS = 0, /*!< \brief _____. */
	DENSITY_FS = 1

};

static const map<string, ENUM_FREESTREAM_OPTION> FreeStreamOption_Map = CCreateMap<string, ENUM_FREESTREAM_OPTION>
("TEMPERATURE_FS", TEMPERATURE_FS)
("DENSITY_FS", DENSITY_FS);

/*!
 * \brief types of viscosity model
 */
enum ENUM_VISCOSITYMODEL {
	CONSTANT_VISCOSITY = 0, /*!< \brief _____. */
	SUTHERLAND = 1,
	LUT_VISCOSITY
};

static const map<string, ENUM_VISCOSITYMODEL> ViscosityModel_Map = CCreateMap<string, ENUM_VISCOSITYMODEL>
("CONSTANT_VISCOSITY", CONSTANT_VISCOSITY)
("SUTHERLAND", SUTHERLAND)
("LUT_VISCOSITY", LUT_VISCOSITY);

/*!
 * \brief types of thermal conductivity model
 */
enum ENUM_CONDUCTIVITYMODEL {
	CONSTANT_CONDUCTIVITY = 0, /*!< \brief _____. */
	CONSTANT_PRANDTL = 1,
	LUT_CONDUCTIVITY
};

static const map<string, ENUM_CONDUCTIVITYMODEL> ConductivityModel_Map = CCreateMap<string, ENUM_CONDUCTIVITYMODEL>
("CONSTANT_CONDUCTIVITY", CONSTANT_CONDUCTIVITY)
("CONSTANT_PRANDTL", CONSTANT_PRANDTL)
("LUT_CONDUCTIVITY", LUT_CONDUCTIVITY);

/*!
 * \brief types of unsteady mesh motion
 */
enum ENUM_GRIDMOVEMENT {
  NO_MOVEMENT = 0, /*!< \brief Simulation on a static mesh. */
  DEFORMING = 1,		/*!< \brief Simulation with dynamically deforming meshes (plunging/pitching/rotation). */
  RIGID_MOTION = 2,		/*!< \brief Simulation with rigid mesh motion (plunging/pitching/rotation). */
  FLUID_STRUCTURE = 3,		/*!< \brief Fluid structure defromation. */
  EXTERNAL = 4,  /*!< \brief Arbitrary grid motion specified by external files at each time step. */
  EXTERNAL_ROTATION = 5,  /*!< \brief Arbitrary grid motion specified by external files at each time step with rigid rotation. */
  AEROELASTIC = 6,    /*!< \brief Simulation with aeroelastic motion. */
  MOVING_WALL = 7,    /*!< \brief Simulation with moving walls (translation/rotation). */
  ROTATING_FRAME = 8,    /*!< \brief Simulation in a rotating frame. */
  ELASTICITY = 9,    /*!< \brief Linear Elasticity. */
  AEROELASTIC_RIGID_MOTION = 10, /*!< \brief Simulation with rotation and aeroelastic motion. */
  STEADY_TRANSLATION = 11,    /*!< \brief Simulation in a steadily translating frame. */
  GUST = 12, /*!< \brief Simulation on a static mesh with a gust. */
  MOVING_HTP = 13    /*!< \brief Simulation with moving HTP (rotation). */

};

static const map<string, ENUM_GRIDMOVEMENT> GridMovement_Map = CCreateMap<string, ENUM_GRIDMOVEMENT>
("NONE", NO_MOVEMENT)
("DEFORMING", DEFORMING)
("RIGID_MOTION", RIGID_MOTION)
("FLUID_STRUCTURE", FLUID_STRUCTURE)
("EXTERNAL", EXTERNAL)
("EXTERNAL_ROTATION", EXTERNAL_ROTATION)
("AEROELASTIC", AEROELASTIC)
("ROTATING_FRAME", ROTATING_FRAME)
("ELASTICITY", ELASTICITY)
("MOVING_WALL", MOVING_WALL)
("MOVING_HTP", MOVING_HTP)
("AEROELASTIC_RIGID_MOTION", AEROELASTIC_RIGID_MOTION)
("STEADY_TRANSLATION", STEADY_TRANSLATION)
("GUST", GUST);

/*!
 * \brief type of wind gusts
 */
enum ENUM_GUST_TYPE {
  NO_GUST = 0,      /*!< \brief _______. */
  TOP_HAT = 1,      /*!< \brief Top-hat function shaped gust  */
  SINE = 2,         /*!< \brief  Sine shaped gust */
  ONE_M_COSINE = 3, /*!< \brief  1-cosine shaped gust */
  VORTEX = 4,       /*!< \brief  A gust made from vortices */
  EOG = 5           /*!< \brief  An extreme operating gust */
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
  X_DIR = 0,        /*!< \brief _______. */
  Y_DIR = 1 		 /*!< \brief _______. */
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
  AUSMPWPLUS = 8,             /*!< \brief AUSMPW+ numerical method. */
  CUSP = 9,                   /*!< \brief Convective upwind and split pressure numerical method. */
  CONVECTIVE_TEMPLATE = 10    /*!< \brief Template for new numerical method . */
};
static const map<string, ENUM_UPWIND> Upwind_Map = CCreateMap<string, ENUM_UPWIND>
("NONE", NO_UPWIND)
("ROE", ROE)
("TURKEL_PREC", TURKEL)
("AUSM", AUSM)
("AUSMPW+", AUSMPWPLUS)
("HLLC", HLLC)
("SW", SW)
("MSW", MSW)
("CUSP", CUSP)
("SCALAR_UPWIND", SCALAR_UPWIND)
("CONVECTIVE_TEMPLATE", CONVECTIVE_TEMPLATE);

/*!
 * \brief types of slope limiters
 */
enum ENUM_LIMITER {
  NO_LIMITER           = 0, /*!< \brief No limiter. */
  VENKATAKRISHNAN      = 1,	/*!< \brief Slope limiter using Venkatakrisnan method (stencil formulation). */
  VENKATAKRISHNAN_WANG = 2,	/*!< \brief Slope limiter using Venkatakrisnan method, eps based on solution (stencil formulation). */
  BARTH_JESPERSEN      = 3, /*!< \brief Slope limiter using Barth-Jespersen method (stencil formulation). */
  VAN_ALBADA_EDGE      = 4, /*!< \brief Slope limiter using Van Albada method (edge formulation). */
  MINMOD_EDGE          = 5, /*!< \brief Slope limiter using Van Albada method (edge formulation). */
  SHARP_EDGES          = 6, /*!< \brief Slope limiter using sharp edges. */
  WALL_DISTANCE        = 7  /*!< \brief Slope limiter using wall distance. */
};
static const map<string, ENUM_LIMITER> Limiter_Map = CCreateMap<string, ENUM_LIMITER>
("NONE", NO_LIMITER)
("VENKATAKRISHNAN", VENKATAKRISHNAN)
("VENKATAKRISHNAN_WANG", VENKATAKRISHNAN_WANG)
("BARTH_JESPERSEN", BARTH_JESPERSEN)
("VAN_ALBADA_EDGE", VAN_ALBADA_EDGE)
("MINMOD_EDGE", MINMOD_EDGE)
("SHARP_EDGES", SHARP_EDGES)
("WALL_DISTANCE", WALL_DISTANCE);

/*!
 * \brief types of turbulent models
 */
enum ENUM_TURB_MODEL {
  NO_TURB_MODEL = 0, /*!< \brief No turbulence model. */
  SA      = 1, /*!< \brief Kind of Turbulent model (Spalart-Allmaras). */
  SA_NEG  = 2, /*!< \brief Kind of Turbulent model (Spalart-Allmaras). */
  SST     = 3, /*!< \brief Kind of Turbulence model (Menter SST). */
};
static const map<string, ENUM_TURB_MODEL> Turb_Model_Map = CCreateMap<string, ENUM_TURB_MODEL>
("NONE", NO_TURB_MODEL)
("SA", SA)
("SA_NEG", SA_NEG)
("SST", SST);

/*!
 * \brief types of transition models
 */
enum ENUM_TRANS_MODEL {
  NO_TRANS_MODEL = 0,            /*!< \brief No transition model. */
  LM = 1,	/*!< \brief Kind of transition model (LM for Spalart-Allmaras). */
  BC = 2	/*!< \brief Kind of transition model (BAS-CAKMAKCIOGLU (BC) for Spalart-Allmaras). */
};
static const map<string, ENUM_TRANS_MODEL> Trans_Model_Map = CCreateMap<string, ENUM_TRANS_MODEL>
("NONE", NO_TRANS_MODEL)
("LM", LM)
("BC", BC); //BAS-CAKMAKCIOGLU

/*!
 * \brief types of wall functions.
 */
enum ENUM_WALL_FUNCTIONS {
  NO_WALL_FUNCTION          = 0,   /*!< \brief No wall function treatment, integration to the wall. Default behavior. */
  STANDARD_WALL_FUNCTION    = 1,   /*!< \brief Standard wall function. */
  ADAPTIVE_WALL_FUNCTION    = 2,   /*!< \brief Adaptive wall function. Formulation depends on y+. */
  SCALABLE_WALL_FUNCTION    = 3,   /*!< \brief Scalable wall function. */
  EQUILIBRIUM_WALL_MODEL    = 4,   /*!< \brief Equilibrium wall model for LES. */
  NONEQUILIBRIUM_WALL_MODEL = 5    /*!< \brief Non-equilibrium wall model for LES. */
};
static const map<string, ENUM_WALL_FUNCTIONS> Wall_Functions_Map = CCreateMap<string, ENUM_WALL_FUNCTIONS>
("NO_WALL_FUNCTION",          NO_WALL_FUNCTION)
("STANDARD_WALL_FUNCTION",    STANDARD_WALL_FUNCTION)
("ADAPTIVE_WALL_FUNCTION",    ADAPTIVE_WALL_FUNCTION)
("SCALABLE_WALL_FUNCTION",    SCALABLE_WALL_FUNCTION)
("EQUILIBRIUM_WALL_MODEL",    EQUILIBRIUM_WALL_MODEL)
("NONEQUILIBRIUM_WALL_MODEL", NONEQUILIBRIUM_WALL_MODEL);

/*!
 * \brief type of time integration schemes
 */
enum ENUM_TIME_INT {
  RUNGE_KUTTA_EXPLICIT = 1,	/*!< \brief Explicit Runge-Kutta time integration definition. */
  EULER_EXPLICIT = 2,   	/*!< \brief Explicit Euler time integration definition. */
  EULER_IMPLICIT = 3,   	/*!< \brief Implicit Euler time integration definition. */
  CLASSICAL_RK4_EXPLICIT = 4,   	/*!< \brief Calssical RK4 time integration definition. */
};
static const map<string, ENUM_TIME_INT> Time_Int_Map = CCreateMap<string, ENUM_TIME_INT>
("RUNGE-KUTTA_EXPLICIT", RUNGE_KUTTA_EXPLICIT)
("EULER_EXPLICIT", EULER_EXPLICIT)
("EULER_IMPLICIT", EULER_IMPLICIT)
("CLASSICAL_RK4_EXPLICIT", CLASSICAL_RK4_EXPLICIT);

/*!
 * \brief type of time integration schemes
 */
enum ENUM_TIME_INT_FEA {
  CD_EXPLICIT = 1,			/*!< \brief Support for implementing an explicit method. */
  NEWMARK_IMPLICIT = 2,   	/*!< \brief Implicit Newmark integration definition. */
  GENERALIZED_ALPHA = 3   		/*!< \brief Support for implementing another implicit method. */
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
 * \brief types of transfer methods
 */
enum ENUM_TRANSFER_METHOD {
  BROADCAST_DATA = 1,	/*!< \brief Gather data on one processor and broadcast it into all of them, relating to global nodes. */
  SCATTER_DATA = 2,   	/*!< \brief Gather data on one processor and scatter it into the one that needs it. */
  ALLGATHER_DATA = 3,   /*!< \brief All processors gather data (this will be useful for operations over a group of data - averaging) */
  LEGACY_METHOD = 4		/*!< \brief Original transfer method, maintained to check . */
};
static const map<string, ENUM_TRANSFER_METHOD> Transfer_Method_Map = CCreateMap<string, ENUM_TRANSFER_METHOD>
("BROADCAST_DATA", BROADCAST_DATA)
("SCATTER_DATA", SCATTER_DATA)
("ALLGATHER_DATA", ALLGATHER_DATA)
("LEGACY_METHOD", LEGACY_METHOD);

/*!
 * \brief types of schemes to compute the flow gradient
 */
enum ENUM_FLOW_GRADIENT {
  GREEN_GAUSS = 1,		/*!< \brief Gradients computation using Green Gauss theorem. */
  WEIGHTED_LEAST_SQUARES = 2	/*!< \brief Gradients computation using Weighted Least Squares. */
};
static const map<string, ENUM_FLOW_GRADIENT> Gradient_Map = CCreateMap<string, ENUM_FLOW_GRADIENT>
("GREEN_GAUSS", GREEN_GAUSS)
("WEIGHTED_LEAST_SQUARES", WEIGHTED_LEAST_SQUARES);

/*!
 * \brief types of action to take on a geometry structure
 */
enum GEOMETRY_ACTION {
  ALLOCATE = 0,     /*!<  \brief Allocate geometry structure. */
  UPDATE = 1       /*!<  \brief Update geometry structure (grid moving, adaptation, etc.). */
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
  EULER_WALL = 1,		/*!< \brief Boundary Euler wall definition. */
  FAR_FIELD = 2,		/*!< \brief Boundary far-field definition. */
  SYMMETRY_PLANE = 3,   	/*!< \brief Boundary symmetry plane definition. */
  INLET_FLOW = 4,		/*!< \brief Boundary inlet flow definition. */
  OUTLET_FLOW = 5,		/*!< \brief Boundary outlet flow definition. */
  PERIODIC_BOUNDARY = 6,	/*!< \brief Periodic boundary definition. */
  NEARFIELD_BOUNDARY = 7,	/*!< \brief Near-Field boundary definition. */
  ELECTRODE_BOUNDARY = 8,	/*!< \brief Electrode boundary definition. */
  DIELEC_BOUNDARY = 9,	/*!< \brief Dipoisson boundary definition. */
  CUSTOM_BOUNDARY = 10,         /*!< \brief custom boundary definition. */
  INTERFACE_BOUNDARY = 11,	/*!< \brief Domain interface boundary definition. */
  DIRICHLET = 12,		/*!< \brief Boundary Euler wall definition. */
  NEUMANN = 13,		/*!< \brief Boundary Neumann definition. */
  DISPLACEMENT_BOUNDARY = 14,		/*!< \brief Boundary displacement definition. */
  LOAD_BOUNDARY = 15,		/*!< \brief Boundary Load definition. */
  FLOWLOAD_BOUNDARY = 16,		/*!< \brief Boundary Load definition. */
  SUPERSONIC_INLET = 19,		/*!< \brief Boundary supersonic inlet definition. */
  SUPERSONIC_OUTLET = 20,		/*!< \brief Boundary supersonic inlet definition. */
  ENGINE_INFLOW = 21,		/*!< \brief Boundary nacelle inflow. */
  ENGINE_EXHAUST = 22,		/*!< \brief Boundary nacelle exhaust. */
  RIEMANN_BOUNDARY= 24,   /*!< \brief Riemann Boundary definition. */
  ISOTHERMAL = 25,      /*!< \brief No slip isothermal wall boundary condition. */
  HEAT_FLUX  = 26,      /*!< \brief No slip constant heat flux wall boundary condition. */
  PRESSURE_BOUNDARY = 27,   	/*!< \brief Pressure boundary condition. */
  ACTDISK_INLET = 32,	/*!< \brief Actuator disk inlet boundary definition. */
  ACTDISK_OUTLET = 33,	/*!< \brief Actuator disk outlet boundary definition. */
  CLAMPED_BOUNDARY = 34,		/*!< \brief Clamped Boundary definition. */
  LOAD_DIR_BOUNDARY = 35,		/*!< \brief Boundary Load definition. */
  LOAD_SINE_BOUNDARY = 36,		/*!< \brief Sine-waveBoundary Load definition. */
  GILES_BOUNDARY= 37,   /*!< \brief Giles Boundary definition. */
  INTERNAL_BOUNDARY= 38,   /*!< \brief Internal Boundary definition. */
  FLUID_INTERFACE = 39,	/*!< \brief Domain interface definition. */
  SEND_RECEIVE = 99,		/*!< \brief Boundary send-receive definition. */
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
 * \brief different regime modes
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
 * \brief types Riemann boundary treatments
 */
enum RIEMANN_TYPE {
  TOTAL_CONDITIONS_PT = 1,		/*!< \brief User specifies total pressure, total temperature, and flow direction. */
  DENSITY_VELOCITY = 2,         /*!< \brief User specifies density and velocity, and flow direction. */
  STATIC_PRESSURE = 3,           /*!< \brief User specifies static pressure. */
  TOTAL_SUPERSONIC_INFLOW = 4,	/*!< \brief User specifies total pressure, total temperature and Velocity components. */
  STATIC_SUPERSONIC_INFLOW_PT = 5, /*!< \brief User specifies static pressure, static temperature, and Mach components. */
  STATIC_SUPERSONIC_INFLOW_PD = 6, /*!< \brief User specifies static pressure, static temperature, and Mach components. */
  MIXING_IN = 7, /*!< \brief User does not specify anything information are retrieved from the other domain */
  MIXING_OUT = 8, /*!< \brief User does not specify anything information are retrieved from the other domain */
  SUPERSONIC_OUTFLOW = 9,
  RADIAL_EQUILIBRIUM = 10,
  TOTAL_CONDITIONS_PT_1D = 11,
  STATIC_PRESSURE_1D = 12,
  MIXING_IN_1D = 13,
  MIXING_OUT_1D =14,
  SPANWISE_TOTAL_CONDITIONS_PT=15,
  SPANWISE_STATIC_PRESSURE=16
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
("STATIC_PRESSURE_1D", STATIC_PRESSURE_1D)
("SPANWISE_TOTAL_CONDITIONS_PT", SPANWISE_TOTAL_CONDITIONS_PT)
("SPANWISE_STATIC_PRESSURE", SPANWISE_STATIC_PRESSURE);



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
("STATIC_PRESSURE_1D", STATIC_PRESSURE_1D)
("SPANWISE_TOTAL_CONDITIONS_PT", SPANWISE_TOTAL_CONDITIONS_PT)
("SPANWISE_STATIC_PRESSURE", SPANWISE_STATIC_PRESSURE);

/*!
 * \brief types of mixing process for averaging quantities at the boundaries.
 */
enum AVERAGEPROCESS_TYPE {
  ALGEBRAIC = 1,		/*!< \brief an algebraic average is computed at the boundary of interest. */
  AREA = 2,           /*!< \brief an area average is computed at the boundary of interest. */
  MIXEDOUT = 3,		 /*!< \brief an mixed-out average is computed at the boundary of interest. */
  MASSFLUX = 4           /*!< \brief a mass flow average is computed at the boundary of interest. */

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
  MATCHING = 1,		/*!< \brief an algebraic average is computed at the boundary of interest. */
  NEAREST_SPAN = 2,           /*!< \brief an area average is computed at the boundary of interest. */
  LINEAR_INTERPOLATION = 3		 /*!< \brief an mixed-out average is computed at the boundary of interest. */
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
  EQUISPACED = 2           /*!< \brief number of span-wise section are specified from the user */
};

static const map<string, SPANWISE_TYPE> SpanWise_Map = CCreateMap<string, SPANWISE_TYPE>
("AUTOMATIC", AUTOMATIC)
("EQUISPACED", EQUISPACED);

/*!
 * \brief types of mixing process for averaging quantities at the boundaries.
 */
enum TURBOMACHINERY_TYPE {
  AXIAL       = 1,		/*!< \brief axial turbomachinery. */
  CENTRIFUGAL = 2,    /*!< \brief centrifugal turbomachinery. */
  CENTRIPETAL = 3,		 /*!< \brief centripetal turbomachinery. */
  CENTRIPETAL_AXIAL = 4,		 /*!< \brief mixed flow turbine. */
  AXIAL_CENTRIFUGAL = 5		 /*!< \brief mixed flow turbine. */
};

static const map<string, TURBOMACHINERY_TYPE> TurboMachinery_Map = CCreateMap<string, TURBOMACHINERY_TYPE>
("AXIAL", AXIAL)
("CENTRIFUGAL", CENTRIFUGAL)
("CENTRIPETAL",  CENTRIPETAL)
("CENTRIPETAL_AXIAL",  CENTRIPETAL_AXIAL)
("AXIAL_CENTRIFUGAL",  AXIAL_CENTRIFUGAL);

///*!
// * \brief types of Turbomachinery performance indicators.
// */
//enum TURBO_PERFORMANCE_TYPE {
//  BLADE   = 1,		/*!< \brief Turbomachinery blade performances. */
//  STAGE = 2,      /*!< \brief Turbomachinery blade stage performances. */
//  TURBINE              = 3		/*!< \brief Turbomachinery turbine performances. */
//};
//
//static const map<string, TURBO_PERFORMANCE_TYPE> TurboPerformance_Map = CCreateMap<string, TURBO_PERFORMANCE_TYPE>
//("BLADE", BLADE)
//("STAGE", STAGE)
//("TURBINE", TURBINE);

/*!
 * \brief types of Turbomachinery performance flag.
 */
enum TURBO_MARKER_TYPE{
  INFLOW   = 1,		/*!< \brief flag for inflow marker for compute turboperformance. */
  OUTFLOW = 2     /*!< \brief flag for outflow marker for compute turboperformance. */
};

/*!
 * \brief types inlet boundary treatments
 */
enum INLET_TYPE {
  TOTAL_CONDITIONS = 1,		/*!< \brief User specifies total pressure, total temperature, and flow direction. */
  MASS_FLOW = 2,           /*!< \brief User specifies density and velocity (mass flow). */
  INPUT_FILE = 3           /*!< \brief User specifies an input file. */
};
static const map<string, INLET_TYPE> Inlet_Map = CCreateMap<string, INLET_TYPE>
("TOTAL_CONDITIONS", TOTAL_CONDITIONS)
("MASS_FLOW", MASS_FLOW)
("INPUT_FILE", INPUT_FILE);

/*!
 * \brief types engine inflow boundary treatments
 */
enum ENGINE_INFLOW_TYPE {
  FAN_FACE_MACH = 1,	         /*!< \brief User specifies fan face mach number. */
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
  VARIABLES_JUMP = 1,		/*!< \brief User specifies the variables jump. */
  BC_THRUST = 2,     /*!< \brief User specifies the BC thrust. */
  NET_THRUST = 3,     /*!< \brief User specifies the Net thrust. */
  DRAG_MINUS_THRUST = 4,     /*!< \brief User specifies the D-T. */
  MASSFLOW = 5,     /*!< \brief User specifies the massflow. */
  POWER = 6     /*!< \brief User specifies the power. */
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
  VERTEX = 1,   		/*!< \brief VTK nomenclature for defining a vertex element. */
  LINE = 3,			/*!< \brief VTK nomenclature for defining a line element. */
  TRIANGLE = 5, 		/*!< \brief VTK nomenclature for defining a triangle element. */
  QUADRILATERAL = 9,		/*!< \brief VTK nomenclature for defining a quadrilateral element. */
  TETRAHEDRON = 10,     	/*!< \brief VTK nomenclature for defining a tetrahedron element. */
  HEXAHEDRON = 12,      	/*!< \brief VTK nomenclature for defining a hexahedron element. */
  PRISM = 13,     		/*!< \brief VTK nomenclature for defining a prism element. */
  PYRAMID = 14  		/*!< \brief VTK nomenclature for defining a pyramid element. */
};

/*!
 * \brief types of objective functions
 */
enum ENUM_OBJECTIVE {
  DRAG_COEFFICIENT = 1, 	      /*!< \brief Drag objective function definition. */
  LIFT_COEFFICIENT = 2, 	      /*!< \brief Lift objective function definition. */
  SIDEFORCE_COEFFICIENT = 3,	  /*!< \brief Side force objective function definition. */
  EFFICIENCY = 4,		            /*!< \brief Efficiency objective function definition. */
  INVERSE_DESIGN_PRESSURE = 5,	/*!< \brief Pressure objective function definition (inverse design). */
  INVERSE_DESIGN_HEATFLUX = 6,  /*!< \brief Heat flux objective function definition (inverse design). */
  TOTAL_HEATFLUX = 7,           /*!< \brief Total heat flux. */
  MAXIMUM_HEATFLUX = 8,         /*!< \brief Maximum heat flux. */
  MOMENT_X_COEFFICIENT = 9,	    /*!< \brief Pitching moment objective function definition. */
  MOMENT_Y_COEFFICIENT = 10,    /*!< \brief Rolling moment objective function definition. */
  MOMENT_Z_COEFFICIENT = 11,    /*!< \brief Yawing objective function definition. */
  EQUIVALENT_AREA = 12,		      /*!< \brief Equivalent area objective function definition. */
  NEARFIELD_PRESSURE = 13,	    /*!< \brief NearField Pressure objective function definition. */
  FORCE_X_COEFFICIENT = 14,	    /*!< \brief X-direction force objective function definition. */
  FORCE_Y_COEFFICIENT = 15,	    /*!< \brief Y-direction force objective function definition. */
  FORCE_Z_COEFFICIENT = 16,	    /*!< \brief Z-direction force objective function definition. */
  THRUST_COEFFICIENT = 17,		  /*!< \brief Thrust objective function definition. */
  TORQUE_COEFFICIENT = 18,		  /*!< \brief Torque objective function definition. */
  FIGURE_OF_MERIT = 19,		      /*!< \brief Rotor Figure of Merit objective function definition. */
  SURFACE_TOTAL_PRESSURE = 28, 	    /*!< \brief Total Pressure objective function definition. */
  SURFACE_STATIC_PRESSURE = 29,      /*!< \brief Static Pressure objective function definition. */
  SURFACE_MASSFLOW = 30,           /*!< \brief Mass Flow Rate objective function definition. */
  SURFACE_MACH = 51,           /*!< \brief Mach number objective function definition. */
  CUSTOM_OBJFUNC = 31, 	           /*!< \brief Custom objective function definition. */
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
  ENTROPY_GENERATION = 50
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
("FIGURE_OF_MERIT", FIGURE_OF_MERIT)
("SURFACE_TOTAL_PRESSURE", SURFACE_TOTAL_PRESSURE)
("SURFACE_STATIC_PRESSURE", SURFACE_STATIC_PRESSURE)
("SURFACE_MASSFLOW", SURFACE_MASSFLOW)
("SURFACE_MACH", SURFACE_MACH)
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
("KINETIC_ENERGY_LOSS", KINETIC_ENERGY_LOSS);

/*!
 * \brief types of residual criteria equations
 */

enum ENUM_RESIDUAL {
	RHO_RESIDUAL = 1, 	      /*!< \brief Rho equation residual criteria equation. */
	RHO_ENERGY_RESIDUAL = 2 	      /*!< \brief RhoE equation residual criteria equation. */
};

static const map<string, ENUM_RESIDUAL> Residual_Map = CCreateMap<string, ENUM_RESIDUAL>
("RHO", RHO_RESIDUAL)
("RHO_ENERGY", RHO_ENERGY_RESIDUAL);

/*!
 * \brief types of sensitivities to compute
 */
enum ENUM_SENS {
  SENS_GEOMETRY = 1,    	/*!< \brief Geometrical sensitivity. */
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
  NO_ADAPT = 0,                 /*!< \brief No grid adaptation. */
  FULL = 1,			/*!< \brief Do a complete grid refinement of all the computational grids. */
  FULL_FLOW = 2,		/*!< \brief Do a complete grid refinement of the flow grid. */
  FULL_ADJOINT = 3,		/*!< \brief Do a complete grid refinement of the adjoint grid. */
  GRAD_FLOW = 5,		/*!< \brief Do a gradient based grid adaptation of the flow grid. */
  GRAD_ADJOINT = 6,		/*!< \brief Do a gradient based grid adaptation of the adjoint grid. */
  GRAD_FLOW_ADJ = 7,		/*!< \brief Do a gradient based grid adaptation of the flow and adjoint grid. */
  COMPUTABLE = 9,		/*!< \brief Apply a computable error grid adaptation. */
  REMAINING = 10,		/*!< \brief Apply a remaining error grid adaptation. */
  WAKE = 12,			/*!< \brief Do a grid refinement on the wake. */
  SMOOTHING = 14,		/*!< \brief Do a grid smoothing of the geometry. */
  SUPERSONIC_SHOCK = 15,	/*!< \brief Do a grid smoothing. */
  PERIODIC = 17			/*!< \brief Add the periodic halo cells. */
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
  SU2 = 1,                       /*!< \brief SU2 input format. */
  CGNS = 2                     /*!< \brief CGNS input format for the computational grid. */
};
static const map<string, ENUM_INPUT> Input_Map = CCreateMap<string, ENUM_INPUT>
("SU2", SU2)
("CGNS", CGNS);

const int CGNS_STRING_SIZE = 33;/*!< \brief Length of strings used in the CGNS format. */

/*!
 * \brief type of solution output file formats
 */
enum ENUM_OUTPUT {
  TECPLOT = 1,  		     /*!< \brief Tecplot format for the solution output. */
  TECPLOT_BINARY = 2,    /*!< \brief Tecplot binary format for the solution output. */
  FIELDVIEW = 3,  		   /*!< \brief FieldView format for the solution output. */
  FIELDVIEW_BINARY = 4,  /*!< \brief FieldView binary format for the solution output. */
  CSV = 5,			         /*!< \brief Comma-separated values format for the solution output. */
  CGNS_SOL = 6,  	     	 /*!< \brief CGNS format for the solution output. */
  PARAVIEW = 7  		     /*!< \brief Paraview format for the solution output. */
};
static const map<string, ENUM_OUTPUT> Output_Map = CCreateMap<string, ENUM_OUTPUT>
("TECPLOT", TECPLOT)
("TECPLOT_BINARY", TECPLOT_BINARY)
("FIELDVIEW", FIELDVIEW)
("FIELDVIEW_BINARY", FIELDVIEW_BINARY)
("CSV", CSV)
("CGNS", CGNS_SOL)
("PARAVIEW", PARAVIEW);

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
  FULLMG_CYCLE = 2  /*!< \brief FullMG cycle. */
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
  PRESSURE = 5, 		/*!< \brief Static pressure. */
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
  TRANSLATION = 0,		       /*!< \brief Surface movement as design variable. */
  ROTATION = 1,			         /*!< \brief Surface rotation as design variable. */
  SCALE = 2,			           /*!< \brief Surface rotation as design variable. */
  FFD_SETTING = 3,		       /*!< \brief No surface deformation. */
  FFD_CONTROL_POINT = 4,	   /*!< \brief Free form deformation for 3D design (change a control point). */
  FFD_NACELLE = 5,	         /*!< \brief Free form deformation for 3D design (change a control point). */
  FFD_GULL = 6,	             /*!< \brief Free form deformation for 3D design (change a control point). */
  FFD_CAMBER = 7,		         /*!< \brief Free form deformation for 3D design (camber change). */
  FFD_TWIST = 8,		         /*!< \brief Free form deformation for 3D design (change the twist angle of a section). */
  FFD_THICKNESS = 9,		     /*!< \brief Free form deformation for 3D design (thickness change). */
  FFD_DIHEDRAL_ANGLE = 10,	 /*!< \brief Free form deformation for 3D design (change the dihedral angle). */
  FFD_ROTATION = 11,		     /*!< \brief Free form deformation for 3D design (rotation around a line). */
  FFD_CONTROL_POINT_2D = 12, /*!< \brief Free form deformation for 2D design (change a control point). */
  FFD_CAMBER_2D = 13,		     /*!< \brief Free form deformation for 3D design (camber change). */
  FFD_THICKNESS_2D = 14,		 /*!< \brief Free form deformation for 3D design (thickness change). */
  FFD_TWIST_2D = 15,		     /*!< \brief Free form deformation for 3D design (camber change). */
  FFD_CONTROL_SURFACE = 16,	 /*!< \brief Free form deformation for 3D design (control surface). */
  HICKS_HENNE = 17,	         /*!< \brief Hicks-Henne bump function for airfoil deformation. */
  PARABOLIC = 18,		         /*!< \brief Parabolic airfoil definition as design variables. */
  NACA_4DIGITS = 19,	       /*!< \brief The four digits NACA airfoil family as design variables. */
  AIRFOIL = 20,		           /*!< \brief Airfoil definition as design variables. */
  CST = 21,                  /*!< \brief CST method with Kulfan parameters for airfoil deformation. */
  SURFACE_BUMP = 22,	       /*!< \brief Surfacebump function for flat surfaces deformation. */
  SURFACE_FILE = 23,		     /*!< Nodal coordinates set using a surface file. */
  NO_DEFORMATION = 24,		   /*!< \brief No Deformation. */
  ANGLE_OF_ATTACK = 101,	   /*!< \brief Angle of attack for airfoils. */
  FFD_ANGLE_OF_ATTACK = 102	 /*!< \brief Angle of attack for FFD problem. */
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
("FFD_DIHEDRAL_ANGLE", FFD_DIHEDRAL_ANGLE)
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
("CST", CST);


/*!
 * \brief types of FFD Blending function
 */
enum ENUM_FFD_BLENDING{
  BSPLINE_UNIFORM = 0,
  BEZIER = 1,
};
static const map<string, ENUM_FFD_BLENDING> Blending_Map = CCreateMap<string, ENUM_FFD_BLENDING>
("BSPLINE_UNIFORM", BSPLINE_UNIFORM)
("BEZIER", BEZIER);

/*!
 * \brief types of solvers for solving linear systems
 */
enum ENUM_LINEAR_SOLVER {
  STEEPEST_DESCENT = 1,		/*!< \brief Steepest descent method for point inversion algoritm (Free-Form). */
  NEWTON = 2,			/*!< \brief Newton method for point inversion algorithm (Free-Form). */
  QUASI_NEWTON = 3,		/*!< \brief Quasi Newton method for point inversion algorithm (Free-Form). */
  CONJUGATE_GRADIENT = 4,	/*!< \brief Preconditionated conjugate gradient method for grid deformation. */
  FGMRES = 5,    	/*!< \brief Flexible Generalized Minimal Residual method. */
  BCGSTAB = 6,	/*!< \brief BCGSTAB - Biconjugate Gradient Stabilized Method (main solver). */
  RESTARTED_FGMRES = 7,  /*!< \brief Flexible Generalized Minimal Residual method with restart. */
  SMOOTHER_LUSGS = 8,  /*!< \brief LU_SGS smoother. */
  SMOOTHER_JACOBI = 9,  /*!< \brief Jacobi smoother. */
  SMOOTHER_ILU = 10,  /*!< \brief ILU smoother. */
  SMOOTHER_LINELET = 11  /*!< \brief Linelet smoother. */
};
static const map<string, ENUM_LINEAR_SOLVER> Linear_Solver_Map = CCreateMap<string, ENUM_LINEAR_SOLVER>
("STEEPEST_DESCENT", STEEPEST_DESCENT)
("NEWTON", NEWTON)
("QUASI_NEWTON", QUASI_NEWTON)
("CONJUGATE_GRADIENT", CONJUGATE_GRADIENT)
("BCGSTAB", BCGSTAB)
("FGMRES", FGMRES)
("RESTARTED_FGMRES", RESTARTED_FGMRES)
("SMOOTHER_LUSGS", SMOOTHER_LUSGS)
("SMOOTHER_JACOBI", SMOOTHER_JACOBI)
("SMOOTHER_LINELET", SMOOTHER_LINELET)
("SMOOTHER_ILU", SMOOTHER_ILU);

/*!
 * \brief types surface continuity at the intersection with the FFD
 */
enum ENUM_FFD_CONTINUITY {
  DERIVATIVE_NONE = 0,		/*!< \brief No derivative continuity. */
  DERIVATIVE_1ST = 1,		/*!< \brief First derivative continuity. */
  DERIVATIVE_2ND = 2,	/*!< \brief Second derivative continuity. */
  USER_INPUT = 3			      /*!< \brief User input. */
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
  CARTESIAN = 0,
  CYLINDRICAL = 1,
  SPHERICAL = 2,
  POLAR = 3
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
  NO_SMOOTH = 0,		/*!< \brief No smoothing. */
  SOBOLEV = 1,		/*!< \brief Sobolev gradient smoothing. */
  BIGRID = 2	/*!< \brief Bi-grid technique smoothing. */
};
static const map<string, ENUM_SENS_SMOOTHING> Sens_Smoothing_Map = CCreateMap<string, ENUM_SENS_SMOOTHING>
("NONE", NO_SMOOTH)
("SOBOLEV", SOBOLEV)
("BIGRID", BIGRID);

/*!
 * \brief types of preconditioners for the linear solver
 */
enum ENUM_LINEAR_SOLVER_PREC {
  JACOBI = 1,		/*!< \brief Jacobi preconditioner. */
  LU_SGS = 2,		/*!< \brief LU SGS preconditioner. */
  LINELET = 3,  /*!< \brief Line implicit preconditioner. */
  ILU = 4       /*!< \brief ILU(0) preconditioner. */
};
static const map<string, ENUM_LINEAR_SOLVER_PREC> Linear_Solver_Prec_Map = CCreateMap<string, ENUM_LINEAR_SOLVER_PREC>
("JACOBI", JACOBI)
("LU_SGS", LU_SGS)
("LINELET", LINELET)
("ILU", ILU);

/*!
 * \brief types of analytic definitions for various geometries
 */
enum ENUM_GEO_ANALYTIC {
  NO_GEO_ANALYTIC = 0,          /*!< \brief No analytic definition of the geometry. */
  NACA0012_AIRFOIL = 1, 	/*!< \brief Use the analytical definition of the NACA0012 for doing the grid adaptation. */
  NACA4412_AIRFOIL = 2, 	/*!< \brief Use the analytical definition of the NACA4412 for doing the grid adaptation. */
  CYLINDER = 3, 		/*!< \brief Use the analytical definition of a cylinder for doing the grid adaptation. */
  BIPARABOLIC = 4       	/*!< \brief Use the analytical definition of a biparabolic airfoil for doing the grid adaptation. */
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
	TWOD_AIRFOIL = 0,   /*!< \brief Airfoil analysis. */
  WING = 1, 	/*!< \brief Wing analysis. */
  FUSELAGE = 2    /*!< \brief Fuselage analysis. */
};
static const map<string, ENUM_GEO_DESCRIPTION> Geo_Description_Map = CCreateMap<string, ENUM_GEO_DESCRIPTION>
("AIRFOIL", TWOD_AIRFOIL)
("WING", WING)
("FUSELAGE", FUSELAGE);

/*!
 * \brief types of schemes for unsteady computations
 */
enum ENUM_UNSTEADY {
  STEADY = 0,            /*!< \brief A steady computation. */
  TIME_STEPPING = 1,		 /*!< \brief Use a time stepping strategy for unsteady computations. */
  DT_STEPPING_1ST = 2,	 /*!< \brief Use a dual time stepping strategy for unsteady computations (1st order). */
  DT_STEPPING_2ND = 3,	 /*!< \brief Use a dual time stepping strategy for unsteady computations (2nd order). */
  ROTATIONAL_FRAME = 4,  /*!< \brief Use a rotational source term. */
  HARMONIC_BALANCE = 5    /*!< \brief Use a harmonic balance source term. */

};
static const map<string, ENUM_UNSTEADY> Unsteady_Map = CCreateMap<string, ENUM_UNSTEADY>
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
  CONSTANT_STIFFNESS = 0,               /*!< \brief Impose a constant stiffness for each element (steel). */
  INVERSE_VOLUME = 1,			/*!< \brief Impose a stiffness for each element that is inversely proportional to cell volume. */
  SOLID_WALL_DISTANCE = 2			/*!< \brief Impose a stiffness for each element that is proportional to the distance from the solid surface. */
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
  D_MACH = 1,   /*!< \brief Derivative with respect to the mach number */
  D_AOA = 2,		 /*!< \brief Derivative with respect to the angle of attack */
  D_PRESSURE = 3, /*!< \brief Derivative with respect to the freestream pressure */
  D_TEMPERATURE = 4,/*!< \brief Derivative with respect to the freestream temperature */
  D_DENSITY = 5,
  D_TURB2LAM = 6,
  D_SIDESLIP = 7,
  D_VISCOSITY = 8,
  D_REYNOLDS = 9,
  D_DESIGN = 10
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
("DESIGN_VARIABLES", D_DESIGN);


enum ENUM_RECORDING {
  CONS_VARS   = 1,
  MESH_COORDS = 2,
  COMBINED    = 3
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

/* END_CONFIG_ENUMS */

class COptionBase {
private:
public:
  COptionBase() {};
  virtual  ~COptionBase() = 0;
  //  virtual string SetValue(string) {SU2MPI::PrintAndFinalize("shouldn't be here"); return "";};
  virtual string SetValue(vector<string>) = 0;
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

    unsigned short nParamDV = 0;
    stringstream ss;
    unsigned int i = 0;
    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) {
      switch (this->design_variable[iDV]) {
        case NO_DEFORMATION:       nParamDV = 0; break;
        case FFD_SETTING:          nParamDV = 0; break;
        case FFD_CONTROL_POINT_2D: nParamDV = 5; break;
        case FFD_CAMBER_2D:        nParamDV = 2; break;
        case FFD_THICKNESS_2D:     nParamDV = 2; break;
        case FFD_TWIST_2D:         nParamDV = 3; break;
        case HICKS_HENNE:          nParamDV = 2; break;
        case SURFACE_BUMP:         nParamDV = 3; break;
        case CST:                  nParamDV = 3; break;
        case ANGLE_OF_ATTACK:      nParamDV = 1; break;
        case SCALE:                nParamDV = 0; break;
        case TRANSLATION:          nParamDV = 3; break;
        case ROTATION:             nParamDV = 6; break;
        case NACA_4DIGITS:         nParamDV = 3; break;
        case PARABOLIC:            nParamDV = 2; break;
        case AIRFOIL:              nParamDV = 2; break;
        case FFD_CONTROL_POINT:    nParamDV = 7; break;
        case FFD_NACELLE:          nParamDV = 6; break;
        case FFD_GULL:             nParamDV = 2; break;
        case FFD_TWIST:            nParamDV = 8; break;
        case FFD_ROTATION:         nParamDV = 7; break;
        case FFD_CONTROL_SURFACE:  nParamDV = 7; break;
        case FFD_CAMBER:           nParamDV = 3; break;
        case FFD_THICKNESS:        nParamDV = 3; break;
        case FFD_ANGLE_OF_ATTACK:  nParamDV = 2; break;
        case SURFACE_FILE:         nParamDV = 0; break;
        default : {
          string newstring;
          newstring.append(this->name);
          newstring.append(": undefined design variable type found in configuration file.");
          return newstring;
        }
      }

      for (unsigned short iParamDV = 0; iParamDV < nParamDV; iParamDV++) {

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
  string SetValue(vector<string>) {
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
