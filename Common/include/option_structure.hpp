/*!
 * \file option_structure.hpp
 * \brief Defines classes for referencing options for easy input in CConfig
 * \author J. Hicken, B. Tracey
 * \version 3.2.7.3 "eagle"
 *
 * Many of the classes in this file are templated, and therefore must
 * be declared and defined here; to keep all elements together, there
 * is no corresponding .cpp file at this time.
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (fpalacios@stanford.edu).
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
 * \param[in,out] str - string we want to convert
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
const unsigned int NO_RK_ITER = 0;		       /*!< \brief No Runge-Kutta iteration. */

const unsigned int MESH_0 = 0; /*!< \brief Definition of the finest grid level. */
const unsigned int MESH_1 = 1; /*!< \brief Definition of the finest grid level. */
const unsigned int ZONE_0 = 0; /*!< \brief Definition of the first grid domain. */
const unsigned int ZONE_1 = 1; /*!< \brief Definition of the first grid domain. */

const double AVOGAD_CONSTANT = 6.0221415E26;	     /*!< \brief Avogardro's constant, number of particles in one kmole. */
const double BOLTZMANN_CONSTANT = 1.3806503E-23;   /*! \brief Boltzmann's constant [J K^-1] */
const double UNIVERSAL_GAS_CONSTANT = 8314.462175; /*! \brief Universal gas constant [J kmol^-1 K^-1] */
const double ELECTRON_CHARGE = 1.60217646E-19;	   /*!< \brief Electronic charge constant. */
const double ELECTRON_MASS = 9.10938188E-31;	     /*!< \brief Mass of an electron. */
const double FREE_PERMITTIVITY = 8.8541878176E-12; /*!< \brief Premittivity of free space. */
const double MAGNETIC_CONSTANT = 1.25663706E-6;    /*!< \brief magnetic permeability of free space. */
const double STANDART_GRAVITY = 9.80665;           /*!< \brief Acceleration due to gravity at surface of earth. */

const double EPS = 1.0E-16;		   /*!< \brief Error scale. */
const double TURB_EPS = 1.0E-16; /*!< \brief Turbulent Error scale. */

const double ONE2 = 0.5;			   /*!< \brief One divided by two. */
const double TWO3 = 2.0 / 3.0;	 /*!< \brief Two divided by three. */
const double FOUR3 = 4.0 / 3.0;  /*!< \brief Four divided by three. */

const double PI_NUMBER = 4.0 * atan(1.0);	/*!< \brief Pi number. */

const int MASTER_NODE = 0;			/*!< \brief Master node for MPI parallelization. */
const int SINGLE_NODE = 1;			/*!< \brief There is only a node in the MPI parallelization. */
const int SINGLE_ZONE = 1;			/*!< \brief There is only a zone. */

const unsigned int N_ELEM_TYPES = 7;           /*!< \brief General output & CGNS defines. */
const unsigned int N_POINTS_LINE = 2;          /*!< \brief General output & CGNS defines. */
const unsigned int N_POINTS_TRIANGLE = 3;      /*!< \brief General output & CGNS defines. */
const unsigned int N_POINTS_QUADRILATERAL = 4; /*!< \brief General output & CGNS defines. */
const unsigned int N_POINTS_TETRAHEDRON = 4;   /*!< \brief General output & CGNS defines. */
const unsigned int N_POINTS_HEXAHEDRON = 8;    /*!< \brief General output & CGNS defines. */
const unsigned int N_POINTS_PYRAMID = 5;       /*!< \brief General output & CGNS defines. */
const unsigned int N_POINTS_WEDGE = 6;         /*!< \brief General output & CGNS defines. */

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
 * \brief different solver types for the CFD component
 */
enum ENUM_SOLVER {
  NO_SOLVER = 0,			/*!< \brief Definition of no solver. */
  EULER = 1,				/*!< \brief Definition of the Euler's solver. */
  NAVIER_STOKES = 2,			/*!< \brief Definition of the Navier-Stokes' solver. */
  RANS = 3,				/*!< \brief Definition of the Reynolds-averaged Navier-Stokes' (RANS) solver. */
  POISSON_EQUATION = 4,       	/*!< \brief Definition of the poisson potential solver. */
  WAVE_EQUATION = 10,	/*!< \brief Definition of the wave solver. */
  HEAT_EQUATION = 29,								/*!< \brief Definition of the heat solver. */
  LINEAR_ELASTICITY = 11,	/*!< \brief Definition of the FEA solver. */
  FLUID_STRUCTURE_EULER = 12,	/*!< \brief Definition of the FEA solver. */
  FLUID_STRUCTURE_NAVIER_STOKES = 13,	/*!< \brief Definition of the FEA solver. */
  FLUID_STRUCTURE_RANS = 14,	/*!< \brief Definition of the FEA solver. */
  ADJ_EULER = 18,			/*!< \brief Definition of the continuous adjoint Euler's solver. */
  ADJ_NAVIER_STOKES = 19,		/*!< \brief Definition of the continuous adjoint Navier-Stokes' solver. */
  ADJ_RANS = 20,				/*!< \brief Definition of the continuous adjoint Reynolds-averaged Navier-Stokes' (RANS) solver. */
  LIN_EULER = 21,			/*!< \brief Definition of the linear Euler's solver. */
  LIN_NAVIER_STOKES = 22,		/*!< \brief Definition of the linear Navier-Stokes' solver. */
  TEMPLATE_SOLVER = 30,                  /*!< \brief Definition of template solver. */
  TNE2_EULER = 31,
  TNE2_NAVIER_STOKES = 32,
  ADJ_TNE2_EULER = 33,
  ADJ_TNE2_NAVIER_STOKES = 34
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
("LIN_EULER", LIN_EULER)
("LIN_NAVIER_STOKES", LIN_NAVIER_STOKES)
("TNE2_EULER", TNE2_EULER)
("TNE2_NAVIER_STOKES", TNE2_NAVIER_STOKES)
("ADJ_TNE2_EULER", ADJ_TNE2_EULER)
("ADJ_TNE2_NAVIER_STOKES", ADJ_TNE2_NAVIER_STOKES)
("WAVE_EQUATION", WAVE_EQUATION)
("HEAT_EQUATION", HEAT_EQUATION)
("LINEAR_ELASTICITY", LINEAR_ELASTICITY)
("FLUID_STRUCTURE_EULER", FLUID_STRUCTURE_EULER)
("FLUID_STRUCTURE_NAVIER_STOKES", FLUID_STRUCTURE_NAVIER_STOKES)
("FLUID_STRUCTURE_RANS", FLUID_STRUCTURE_RANS)
("TEMPLATE_SOLVER", TEMPLATE_SOLVER);

/*!
 * \brief different regime modes
 */
enum ENUM_REGIME {
  COMPRESSIBLE = 0,			/*!< \brief Definition of compressible solver. */
  INCOMPRESSIBLE = 1,				/*!< \brief Definition of incompressible solver. */
  FREESURFACE = 2			/*!< \brief Definition of freesurface solver (incompressible). */
};
static const map<string, ENUM_REGIME> Regime_Map = CCreateMap<string, ENUM_REGIME>
("COMPRESSIBLE", COMPRESSIBLE)
("INCOMPRESSIBLE", INCOMPRESSIBLE)
("FREESURFACE", FREESURFACE);

/*!
 * \brief different system of measurements
 */
enum ENUM_MEASUREMENTS {
  SI = 0,			/*!< \brief Definition of compressible solver. */
  US = 1,				/*!< \brief Definition of incompressible solver. */
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
  RUNTIME_LINPOT_SYS = 9,		/*!< \brief One-physics case, the code is solving the linear potential flow equations. */
  RUNTIME_LINFLOW_SYS = 10,		/*!< \brief One-physics case, the code is solving the linear equations is being solved (Euler and Navier-Stokes). */
  RUNTIME_MULTIGRID_SYS = 14,   	/*!< \brief Full Approximation Storage Multigrid system of equations. */
  RUNTIME_FEA_SYS = 20,		/*!< \brief One-physics case, the code is solving the FEA equation. */
  RUNTIME_HEAT_SYS = 21,		/*!< \brief One-physics case, the code is solving the heat equation. */
  RUNTIME_TRANS_SYS = 22,			/*!< \brief One-physics case, the code is solving the turbulence model. */
  RUNTIME_TNE2_SYS = 23,  /*!< \brief One-physics case, the code is solving the two-temperature model. */
  RUNTIME_ADJTNE2_SYS = 24  /*!< \brief One-physics case, the code is solving the two-temperature model. */
};

const int FLOW_SOL = 0;		/*!< \brief Position of the mean flow solution in the solver container array. */
const int ADJFLOW_SOL = 1;	/*!< \brief Position of the continuous adjoint flow solution in the solver container array. */
const int LINFLOW_SOL = 1;	/*!< \brief Position of the linearized flow solution in the solution solver array. */

const int TURB_SOL = 2;		/*!< \brief Position of the turbulence model solution in the solver container array. */
const int ADJTURB_SOL = 3;	/*!< \brief Position of the continuous adjoint turbulence solution in the solver container array. */
const int LINTURB_SOL = 3;	/*!< \brief Position of the linearized turbulence model in the solver container array. */

const int TNE2_SOL = 0;		/*!< \brief Position of the mean flow solution in the solution container array. */
const int ADJTNE2_SOL = 1;	/*!< \brief Position of the continuous adjoint flow solution in the solution container array. */
const int LINTNE2_SOL = 1;	/*!< \brief Position of the linearized flow solution in the solution container array. */

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

/*!
 * \brief types of mathematical problem to solve
 */
enum ENUM_MATH_PROBLEM {
  DIRECT_PROBLEM = 0,		/*!< \brief Direct problem */
  ADJOINT_PROBLEM = 1,		/*!< \brief Adjoint problem */
  LINEARIZED_PROBLEM = 2 /*< \brief Linearized numerical method */
};
static const map<string, ENUM_MATH_PROBLEM> Math_Problem_Map = CCreateMap<string, ENUM_MATH_PROBLEM>
("DIRECT", DIRECT_PROBLEM)
("ADJOINT", ADJOINT_PROBLEM)
("LINEARIZED", LINEARIZED_PROBLEM);


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
	PR_GAS = 3
};

static const map<string, ENUM_FLUIDMODEL> FluidModel_Map = CCreateMap<string, ENUM_FLUIDMODEL>
("STANDARD_AIR", STANDARD_AIR)
("IDEAL_GAS", IDEAL_GAS)
("VW_GAS", VW_GAS)
("PR_GAS", PR_GAS);


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
	SUTHERLAND = 1
};

static const map<string, ENUM_VISCOSITYMODEL> ViscosityModel_Map = CCreateMap<string, ENUM_VISCOSITYMODEL>
("CONSTANT_VISCOSITY", CONSTANT_VISCOSITY)
("SUTHERLAND", SUTHERLAND);

/*!
 * \brief types of thermal conductivity model
 */
enum ENUM_CONDUCTIVITYMODEL {
	CONSTANT_CONDUCTIVITY = 0, /*!< \brief _____. */
	CONSTANT_PRANDTL = 1
};

static const map<string, ENUM_CONDUCTIVITYMODEL> ConductivityModel_Map = CCreateMap<string, ENUM_CONDUCTIVITYMODEL>
("CONSTANT_CONDUCTIVITY", CONSTANT_CONDUCTIVITY)
("CONSTANT_PRANDTL", CONSTANT_PRANDTL);

/*!
 * \brief types of spatial discretizations
 */
enum ENUM_GASMODEL {
  NO_MODEL = 0, /*!< \brief _____. */
  ARGON = 1,		/*!< \brief _____. */
  AIR7 = 2,		/*!< \brief _______. */
  AIR21 = 3,		/*!< \brief _______. */
  O2 = 4,
  N2 = 5,
  AIR5 = 6,
  ARGON_SID = 7,
  ONESPECIES = 8

};
static const map<string, ENUM_GASMODEL> GasModel_Map = CCreateMap<string, ENUM_GASMODEL>
("NONE", NO_MODEL)
("ARGON", ARGON)
("AIR-7", AIR7)
("AIR-21", AIR21)
("O2", O2)
("N2", N2)
("AIR-5", AIR5)
("ARGON-SID",ARGON_SID)
("ONESPECIES", ONESPECIES);

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
  AEROELASTIC_RIGID_MOTION = 10 /*!< \brief Simulation with rotation and aeroelastic motion. */

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
("AEROELASTIC_RIGID_MOTION", AEROELASTIC_RIGID_MOTION);

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
 * \brief Spatial numerical order integration
 */
enum ENUM_SPATIAL_ORDER {
  FIRST_ORDER = 0,        /*!< \brief First order */
  SECOND_ORDER = 1,        /*!< \brief Second order. */
  SECOND_ORDER_LIMITER = 2 /*!< \brief Second order with limiter. */
};
static const map<string, ENUM_SPATIAL_ORDER> SpatialOrder_Map = CCreateMap<string, ENUM_SPATIAL_ORDER>
("1ST_ORDER", FIRST_ORDER)
("2ND_ORDER", SECOND_ORDER)
("2ND_ORDER_LIMITER", SECOND_ORDER_LIMITER);

/*!
 * \brief types of slope limiters
 */
enum ENUM_LIMITER {
  VENKATAKRISHNAN = 0,	/*!< \brief Slope limiter using Venkatakrisnan method. */
  BARTH_JESPERSEN = 1,  /*!< \brief Slope limiter using Barth-Jespersen method. */
  SHARP_EDGES = 2       /*!< \brief Slope limiter using sharp edges. */
};
static const map<string, ENUM_LIMITER> Limiter_Map = CCreateMap<string, ENUM_LIMITER>
("VENKATAKRISHNAN", VENKATAKRISHNAN)
("BARTH_JESPERSEN", BARTH_JESPERSEN)
("SHARP_EDGES", SHARP_EDGES);

/*!
 * \brief types of methods used to calculate source term Jacobians
 */
enum ENUM_SOURCEJAC {
  NO_JACOBIAN = 0,             /*!< \brief No source term Jacobian. */
  FINITE_DIFF = 1,	/*!< \brief Numerical method for source term in flow equations. */
  AUTO_DIFF = 2		/*!< \brief Numerical method for source term in charge distribution. */
};
static const map<string, ENUM_SOURCEJAC> SourceJac_Map = CCreateMap<string, ENUM_SOURCEJAC>
("NO_JACOBIAN", NO_JACOBIAN)
("FINITE_DIFF", FINITE_DIFF)
("AUTO_DIFF", AUTO_DIFF);

/*!
 * \brief types of turbulent models
 */
enum ENUM_TURB_MODEL {
  NO_TURB_MODEL = 0,            /*!< \brief No turbulence model. */
  SA = 1,                       /*!< \brief Kind of Turbulent model (Spalart-Allmaras). */
  SST = 2,       		/*!< \brief Kind of Turbulence model (Menter SST). */
  ML  = 3       		/*!< \brief Kind of Turbulence model (Machine Learning). */
};
static const map<string, ENUM_TURB_MODEL> Turb_Model_Map = CCreateMap<string, ENUM_TURB_MODEL>
("NONE", NO_TURB_MODEL)
("SA", SA)
("SST", SST)
("ML", ML);

/*!
 * \brief types of transition models
 */
enum ENUM_TRANS_MODEL {
  NO_TRANS_MODEL = 0,            /*!< \brief No transition model. */
  LM = 1												/*!< \brief Kind of transition model (LM for Spalart-Allmaras). */
};
static const map<string, ENUM_TRANS_MODEL> Trans_Model_Map = CCreateMap<string, ENUM_TRANS_MODEL>
("NONE", NO_TRANS_MODEL)
("LM", LM);

/*!
 * \brief type of time integration schemes
 */
enum ENUM_TIME_INT {
  RUNGE_KUTTA_EXPLICIT = 1,	/*!< \brief Explicit Runge-Kutta time integration definition. */
  EULER_EXPLICIT = 2,   	/*!< \brief Explicit Euler time integration definition. */
  EULER_IMPLICIT = 3   	/*!< \brief Implicit Euler time integration definition. */
};
static const map<string, ENUM_TIME_INT> Time_Int_Map = CCreateMap<string, ENUM_TIME_INT>
("RUNGE-KUTTA_EXPLICIT", RUNGE_KUTTA_EXPLICIT)
("EULER_EXPLICIT", EULER_EXPLICIT)
("EULER_IMPLICIT", EULER_IMPLICIT);

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
  ELEC_DIELEC_BOUNDARY = 17,	/*!< \brief Dipoisson boundary definition for the poissonal potential. */
  ELEC_NEUMANN = 18,		/*!< \brief Boundary Neumann definition. */
  SUPERSONIC_INLET = 19,		/*!< \brief Boundary supersonic inlet definition. */
  SUPERSONIC_OUTLET = 20,		/*!< \brief Boundary supersonic inlet definition. */
  ENGINE_INFLOW = 21,		/*!< \brief Boundary nacelle inflow. */
  ENGINE_EXHAUST = 22,		/*!< \brief Boundary nacelle exhaust. */
  ENGINE_BLEED = 23,		/*!< \brief Boundary engine bleed. */
  RIEMANN_BOUNDARY= 24,   /*!< \brief Riemann Boundary definition. */
  ISOTHERMAL = 25,      /*!< \brief No slip isothermal wall boundary condition. */
  HEAT_FLUX  = 26,      /*!< \brief No slip constant heat flux wall boundary condition. */
  PRESSURE_BOUNDARY = 27,   	/*!< \brief Pressure boundary condition. */
  HEAT_FLUX_NONCATALYTIC = 28, /*!< \brief No-slip, constant heat flux, noncatalytic bc. */
  HEAT_FLUX_CATALYTIC= 29, /*!< \brief No-slip, constant heat flux, catalytic bc. */
  ISOTHERMAL_NONCATALYTIC = 30, /*!< \brief No-slip, constant temperature, noncatalytic bc. */
  ISOTHERMAL_CATALYTIC = 31, /*!< \brief No-slip, constant temperature, catalytic bc. */
  ACTDISK_INLET = 32,	/*!< \brief Actuator disk inlet boundary definition. */
  ACTDISK_OUTLET = 33,	/*!< \brief Actuator disk outlet boundary definition. */
  SEND_RECEIVE = 99,		/*!< \brief Boundary send-receive definition. */
};

/*!
 * \brief types inlet boundary treatments
 */
enum RIEMANN_TYPE {
  TOTAL_CONDITIONS_PT = 1,		/*!< \brief User specifies total pressure, total temperature, and flow direction. */
  DENSITY_VELOCITY = 2,         /*!< \brief User specifies density and velocity, and flow direction. */
  STATIC_PRESSURE = 3,           /*!< \brief User specifies static pressure. */
  TOTAL_SUPERSONIC_INFLOW = 4,	/*!< \brief User specifies total pressure, total temperature and Velocity components. */
  STATIC_SUPERSONIC_INFLOW_PT = 5, /*!< \brief User specifies static pressure, static temperature, and Mach components. */
  STATIC_SUPERSONIC_INFLOW_PD = 6 /*!< \brief User specifies static pressure, static temperature, and Mach components. */
};

static const map<string, RIEMANN_TYPE> Riemann_Map = CCreateMap<string, RIEMANN_TYPE>
("TOTAL_CONDITIONS_PT", TOTAL_CONDITIONS_PT)
("DENSITY_VELOCITY", DENSITY_VELOCITY)
("STATIC_PRESSURE", STATIC_PRESSURE)
("TOTAL_SUPERSONIC_INFLOW", TOTAL_SUPERSONIC_INFLOW)
("STATIC_SUPERSONIC_INFLOW_PT", STATIC_SUPERSONIC_INFLOW_PT)
("STATIC_SUPERSONIC_INFLOW_PD", STATIC_SUPERSONIC_INFLOW_PD);

/*!
 * \brief types inlet boundary treatments
 */
enum INLET_TYPE {
  TOTAL_CONDITIONS = 1,		/*!< \brief User specifies total pressure, total temperature, and flow direction. */
  MASS_FLOW = 2           /*!< \brief User specifies density and velocity (mass flow). */
};
static const map<string, INLET_TYPE> Inlet_Map = CCreateMap<string, INLET_TYPE>
("TOTAL_CONDITIONS", TOTAL_CONDITIONS)
("MASS_FLOW", MASS_FLOW);

/*!
 * \brief types of geometric entities based on VTK nomenclature
 */
enum GEO_TYPE {
  VERTEX = 1,   		/*!< \brief VTK nomenclature for defining a vertex element. */
  LINE = 3,			/*!< \brief VTK nomenclature for defining a line element. */
  TRIANGLE = 5, 		/*!< \brief VTK nomenclature for defining a triangle element. */
  RECTANGLE = 9,		/*!< \brief VTK nomenclature for defining a rectangle element. */
  TETRAHEDRON = 10,     	/*!< \brief VTK nomenclature for defining a tetrahedron element. */
  HEXAHEDRON = 12,      	/*!< \brief VTK nomenclature for defining a hexahedron element. */
  WEDGE = 13,     		/*!< \brief VTK nomenclature for defining a wedge element. */
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
  FREE_SURFACE = 20,				    /*!< \brief Free Surface objective function definition. */
  MAX_THICKNESS = 21,           /*!< \brief Maximum thickness. */
  MIN_THICKNESS = 22,           /*!< \brief Minimum thickness. */
  MAX_THICK_SEC1 = 23,          /*!< \brief Maximum thickness in section 1. */
  MAX_THICK_SEC2 = 24,          /*!< \brief Maximum thickness in section 2. */
  MAX_THICK_SEC3 = 25,          /*!< \brief Maximum thickness in section 3. */
  MAX_THICK_SEC4 = 26,          /*!< \brief Maximum thickness in section 4. */
  MAX_THICK_SEC5 = 27,           /*!< \brief Maximum thickness in section 5. */
  AVG_TOTAL_PRESSURE = 28, 	    /*!< \brief Total Pressure objective function definition. */
  AVG_OUTLET_PRESSURE = 29,      /*!< \brief Static Pressure objective function definition. */
  MASS_FLOW_RATE = 30           /*!< \brief Mass Flow Rate objective function definition. */
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
("FREE_SURFACE", FREE_SURFACE)
("MAX_THICKNESS", MAX_THICKNESS)
("MIN_THICKNESS", MIN_THICKNESS)
("MAX_THICK_SEC1", MAX_THICK_SEC1)
("MAX_THICK_SEC2", MAX_THICK_SEC2)
("MAX_THICK_SEC3", MAX_THICK_SEC3)
("MAX_THICK_SEC4", MAX_THICK_SEC4)
("MAX_THICK_SEC5", MAX_THICK_SEC5)
("AVG_TOTAL_PRESSURE", AVG_TOTAL_PRESSURE)
("AVG_OUTLET_PRESSURE",AVG_OUTLET_PRESSURE)
("MASS_FLOW_RATE", MASS_FLOW_RATE);

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
 * \brief types of Continuous equations
 */
enum ENUM_CONTINUOUS_EQNS {
  EULER_EQNS= 1, 	/*!< \brief Euler equations. */
  NAVIER_STOKES_EQNS = 2 	/*!< \brief Navier Stokes equations. */
};
static const map<string, ENUM_CONTINUOUS_EQNS> ContinuousEqns_Map = CCreateMap<string, ENUM_CONTINUOUS_EQNS>
("EULER", EULER_EQNS)
("NAVIER_STOKES", NAVIER_STOKES_EQNS);

/*!
 * \brief types of Discrete equations
 */
enum ENUM_DISCRETE_EQNS {
  NONE_EQNS= 1, 	/*!< \brief No equations. */
  SA_EQNS = 2, 	/*!< \brief Spallart-Almaras equations. */
  SST_EQNS = 2 	/*!< \brief SST equations. */
};
static const map<string, ENUM_DISCRETE_EQNS> DiscreteEqns_Map = CCreateMap<string, ENUM_DISCRETE_EQNS>
("NONE", NONE_EQNS)
("SA", SA_EQNS)
("SST", SST_EQNS);

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
 * \brief types of linearized objective functions
 */
enum ENUM_LINEAR_OBJ {
  DELTA_DRAG_COEFFICIENT = 1,	/*!< \brief Linearized Drag objective function definition. */
  DELTA_LIFT_COEFFICIENT = 2	/*!< \brief Linearized Lift objective function definition. */
};
static const map<string, ENUM_LINEAR_OBJ> Linear_Obj_Map = CCreateMap<string, ENUM_LINEAR_OBJ>
("DELTA_DRAG", DELTA_DRAG_COEFFICIENT)
("DELTA_LIFT", DELTA_LIFT_COEFFICIENT);

/*!
 * \brief types of grid adaptation/refinement
 */
enum ENUM_ADAPT {
  NO_ADAPT = 0,                 /*!< \brief No grid adaptation. */
  FULL = 1,			/*!< \brief Do a complete grid refinement of all the computational grids. */
  FULL_FLOW = 2,		/*!< \brief Do a complete grid refinement of the flow grid. */
  FULL_ADJOINT = 3,		/*!< \brief Do a complete grid refinement of the adjoint grid. */
  FULL_LINEAR = 4,		/*!< \brief Do a complete grid refinement of the linear grid. */
  GRAD_FLOW = 5,		/*!< \brief Do a gradient based grid adaptation of the flow grid. */
  GRAD_ADJOINT = 6,		/*!< \brief Do a gradient based grid adaptation of the adjoint grid. */
  GRAD_FLOW_ADJ = 7,		/*!< \brief Do a gradient based grid adaptation of the flow and adjoint grid. */
  ROBUST = 8,			/*!< \brief Apply a robust grid adaptation (residual based). */
  COMPUTABLE = 9,		/*!< \brief Apply a computable error grid adaptation. */
  REMAINING = 10,		/*!< \brief Apply a remaining error grid adaptation. */
  COMPUTABLE_ROBUST = 11,	/*!< \brief Apply a computable robust grid adaptation (using linearized problem). */
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
("FULL_LINEAR", FULL_LINEAR)
("GRAD_FLOW", GRAD_FLOW)
("GRAD_ADJOINT", GRAD_ADJOINT)
("GRAD_FLOW_ADJ", GRAD_FLOW_ADJ)
("ROBUST", ROBUST)
("COMPUTABLE", COMPUTABLE)
("REMAINING", REMAINING)
("COMPUTABLE_ROBUST", COMPUTABLE_ROBUST)
("WAKE", WAKE)
("SMOOTHING", SMOOTHING)
("SUPERSONIC_SHOCK", SUPERSONIC_SHOCK)
("PERIODIC", PERIODIC);

/*!
 * \brief types of input file formats
 */
enum ENUM_INPUT {
  SU2 = 1,                       /*!< \brief SU2 input format. */
  CGNS = 2,                     /*!< \brief CGNS input format for the computational grid. */
  NETCDF_ASCII = 3      	/*!< \brief ASCII NETCDF input format for the computational grid. */
};
static const map<string, ENUM_INPUT> Input_Map = CCreateMap<string, ENUM_INPUT>
("SU2", SU2)
("CGNS", CGNS)
("NETCDF_ASCII", NETCDF_ASCII);

const int CGNS_STRING_SIZE = 33;/*!< \brief Length of strings used in the CGNS format. */

/*!
 * \brief type of solution output file formats
 */
enum ENUM_OUTPUT {
  TECPLOT = 1,  		/*!< \brief Tecplot format for the solution output. */
  EXCEL = 2,			/*!< \brief Excel format for the solution output. */
  CSV = 3,			/*!< \brief Comma-separated values format for the solution output. */
  TECPLOT_BINARY = 4,  		/*!< \brief Tecplot binary format for the solution output. */
  CGNS_SOL = 5,  		/*!< \brief CGNS format for the solution output. */
  PARAVIEW = 6  		/*!< \brief Paraview format for the solution output. */
};
static const map<string, ENUM_OUTPUT> Output_Map = CCreateMap<string, ENUM_OUTPUT>
("TECPLOT", TECPLOT)
("EXCEL", EXCEL)
("CSV", CSV)
("TECPLOT_BINARY", TECPLOT_BINARY)
("CGNS", CGNS_SOL)
("PARAVIEW", PARAVIEW);

/*!
 * \brief type of multigrid cycle
 */
enum MG_CYCLE {
  V_CYCLE = 0,  		/*!< \brief V cycle. */
  W_CYCLE = 1,			/*!< \brief W cycle. */
  FULLMG_CYCLE = 2,			/*!< \brief FullMG cycle. */
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
  FFD_SETTING = 0,		/*!< \brief No surface deformation. */
  HICKS_HENNE = 1,		/*!< \brief Hicks-Henne bump function for airfoil deformation. */
  NACA_4DIGITS = 6,		/*!< \brief The four digits NACA airfoil family as design variables. */
  DISPLACEMENT = 8,		/*!< \brief Surface movement as design variable. */
  ROTATION = 9,			/*!< \brief Surface rotation as design variable. */
  FFD_CONTROL_POINT = 10,	/*!< \brief Free form deformation for 3D design (change a control point). */
  FFD_DIHEDRAL_ANGLE = 11,	/*!< \brief Free form deformation for 3D design (change the dihedral angle). */
  FFD_TWIST_ANGLE = 12,		/*!< \brief Free form deformation for 3D design (change the twist angle). */
  FFD_ROTATION = 13,		/*!< \brief Free form deformation for 3D design (rotation around a line). */
  FFD_CAMBER = 14,		/*!< \brief Free form deformation for 3D design (camber change). */
  FFD_THICKNESS = 15,		/*!< \brief Free form deformation for 3D design (thickness change). */
  PARABOLIC = 17,		/*!< \brief Parabolic airfoil definition as design variables. */
  OBSTACLE = 18,		        /*!< \brief Obstacle for free surface optimization. */
  STRETCH = 19,		        /*!< \brief Stretch one side of a channel. */
  SURFACE_FILE = 20,		   /*!< Nodal coordinates set using a surface file. */
  COSINE_BUMP = 21,		/*!< \brief Gauss bump function for airfoil deformation. */
  FOURIER = 22,		/*!< \brief Fourier function for airfoil deformation. */
  SPHERICAL = 23,		/*!< \brief Spherical geometry parameterization with spline-based radial profile. */
  AIRFOIL = 24,		/*!< \brief Airfoil definition as design variables. */
  FFD_CONTROL_POINT_2D = 25,	/*!< \brief Free form deformation for 2D design (change a control point). */
  FFD_CAMBER_2D = 26,		/*!< \brief Free form deformation for 3D design (camber change). */
  FFD_THICKNESS_2D = 27,		/*!< \brief Free form deformation for 3D design (thickness change). */
  FFD_CONTROL_SURFACE = 28,		/*!< \brief Free form deformation for 3D design (control surface). */
  FFD_RADIUS_2D = 29		/*!< \brief Free form deformation for 2D design (radious change). */
};
static const map<string, ENUM_PARAM> Param_Map = CCreateMap<string, ENUM_PARAM>
("FFD_SETTING", FFD_SETTING)
("FFD_CONTROL_POINT_2D", FFD_CONTROL_POINT_2D)
("FFD_RADIUS_2D", FFD_RADIUS_2D)
("FFD_CAMBER_2D", FFD_CAMBER_2D)
("FFD_THICKNESS_2D", FFD_THICKNESS_2D)
("HICKS_HENNE", HICKS_HENNE)
("SPHERICAL", SPHERICAL)
("NACA_4DIGITS", NACA_4DIGITS)
("DISPLACEMENT", DISPLACEMENT)
("ROTATION", ROTATION)
("FFD_CONTROL_POINT", FFD_CONTROL_POINT)
("FFD_DIHEDRAL_ANGLE", FFD_DIHEDRAL_ANGLE)
("FFD_TWIST_ANGLE", FFD_TWIST_ANGLE)
("FFD_ROTATION", FFD_ROTATION)
("FFD_CONTROL_SURFACE", FFD_CONTROL_SURFACE)
("FFD_CAMBER", FFD_CAMBER)
("FFD_THICKNESS", FFD_THICKNESS)
("PARABOLIC", PARABOLIC)
("OBSTACLE", OBSTACLE)
("STRETCH", STRETCH)
("COSINE_BUMP", COSINE_BUMP)
("FOURIER", FOURIER)
("AIRFOIL", AIRFOIL)
("SURFACE_FILE", SURFACE_FILE);

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
  RFGMRES = 7,  /*!< \brief Flexible Generalized Minimal Residual method with restart. */
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
("RFGMRES", RFGMRES)
("SMOOTHER_LUSGS", SMOOTHER_LUSGS)
("SMOOTHER_JACOBI", SMOOTHER_JACOBI)
("SMOOTHER_LINELET", SMOOTHER_LINELET)
("SMOOTHER_ILU0", SMOOTHER_ILU);

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
("ILU0", ILU);

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
enum ENUM_AXIS_ORIENTATION {
  X_AXIS = 0,   /*!< \brief X axis orientation. */
  Y_AXIS = 1, 	/*!< \brief Y axis orientation. */
  Z_AXIS = 2    /*!< \brief Z axis orientation. */
};
static const map<string, ENUM_AXIS_ORIENTATION> Axis_Orientation_Map = CCreateMap<string, ENUM_AXIS_ORIENTATION>
("X_AXIS", X_AXIS)
("Y_AXIS", Y_AXIS)
("Z_AXIS", Z_AXIS);

/*!
 * \brief types of schemes for unsteady computations
 */
enum ENUM_UNSTEADY {
  STEADY = 0,             /*!< \brief A steady computation. */
  TIME_STEPPING = 1,		/*!< \brief Use a time stepping strategy for unsteady computations. */
  DT_STEPPING_1ST = 2,	/*!< \brief Use a dual time stepping strategy for unsteady computations (1st order). */
  DT_STEPPING_2ND = 3,	/*!< \brief Use a dual time stepping strategy for unsteady computations (2nd order). */
  ROTATIONAL_FRAME = 4,   /*!< \brief Use a rotational source term. */
  TIME_SPECTRAL = 5       	/*!< \brief Use a time spectral source term. */

};
static const map<string, ENUM_UNSTEADY> Unsteady_Map = CCreateMap<string, ENUM_UNSTEADY>
("NO", STEADY)
("TIME_STEPPING", TIME_STEPPING)
("DUAL_TIME_STEPPING-1ST_ORDER", DT_STEPPING_1ST)
("DUAL_TIME_STEPPING-2ND_ORDER", DT_STEPPING_2ND)
("TIME_SPECTRAL", TIME_SPECTRAL)
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
  WALL_DISTANCE = 2			/*!< \brief Impose a stiffness for each element that is proportional to the distance from the deforming surface. */
};
static const map<string, ENUM_DEFORM_STIFFNESS> Deform_Stiffness_Map = CCreateMap<string, ENUM_DEFORM_STIFFNESS>
("CONSTANT_STIFFNESS", CONSTANT_STIFFNESS)
("INVERSE_VOLUME", INVERSE_VOLUME)
("WALL_DISTANCE", WALL_DISTANCE);

/* END_CONFIG_ENUMS */












class COptionBase{
private:
public:
  COptionBase(){};
  virtual  ~COptionBase() = 0;
  //  virtual string SetValue(string){SU2MPI::PrintAndFinalize("shouldn't be here"); return "";};
  virtual string SetValue(vector<string>) = 0;
  virtual void SetDefault() = 0;

  string optionCheckMultipleValues(vector<string> & option_value, string type_id, string option_name){
    if (option_value.size() != 1){
      string newString;
      newString.append(option_name);
      newString.append(": multiple values for type ");
      newString.append(type_id);
      return newString;
    }
    return "";
  }

  string badValue(vector<string> & option_value, string type_id, string option_name){
    string newString;
    newString.append(option_name);
    newString.append(": improper option value for type ");
    newString.append(type_id);
    return newString;
  }
};

inline COptionBase::~COptionBase(){}


template <class Tenum>
class COptionEnum : public COptionBase{

  map<string, Tenum> m;
  unsigned short & field; // Reference to the feildname
  Tenum def; // Default value
  string name; // identifier for the option

public:
  COptionEnum(string option_field_name, const map<string, Tenum> m, unsigned short & option_field, Tenum default_value) : field(option_field){
    this->m = m;
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionEnum(){};
  string SetValue(vector<string> option_value){
    // Check if there is more than one string
    string out = optionCheckMultipleValues(option_value, "enum", this->name);
    if (out.compare("") != 0){
      return out;
    }

    // Check to see if the enum value is in the map
    if (this->m.find(option_value[0]) == m.end()){
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

  void SetDefault(){
    this->field = this->def;
  }
};

class COptionDouble : public COptionBase{
  double & field; // Reference to the fieldname
  double def; // Default value
  string name; // identifier for the option

public:
  COptionDouble(string option_field_name, double & option_field, double default_value) : field(option_field){
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionDouble(){};
  string SetValue(vector<string> option_value){
    // check if there is more than one value
    string out = optionCheckMultipleValues(option_value, "double", this->name);
    if (out.compare("") != 0){
      return out;
    }
    istringstream is(option_value[0]);
    double val;
    if (is >> val){
      this->field = val;
      return "";
    }
    return badValue(option_value, "double", this->name);
  }
  void SetDefault(){
    this->field = this->def;
  }
};

class COptionString : public COptionBase{
  string & field; // Reference to the fieldname
  string def; // Default value
  string name; // identifier for the option

public:
  COptionString(string option_field_name, string & option_field, string default_value) : field(option_field){
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionString(){};
  string SetValue(vector<string> option_value){
    // check if there is more than one value
    string out = optionCheckMultipleValues(option_value, "double", this->name);
    if (out.compare("") != 0){
      return out;
    }
    this->field.assign(option_value[0]);
    return "";
  }
  void SetDefault(){
    this->field = this->def;
  }
};

class COptionInt : public COptionBase{
  int & field; // Reference to the feildname
  int def; // Default value
  string name; // identifier for the option

public:
  COptionInt(string option_field_name, int & option_field, int default_value) : field(option_field){
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionInt(){};
  string SetValue(vector<string> option_value){
    string out = optionCheckMultipleValues(option_value, "int", this->name);
    if (out.compare("") != 0){
      return out;
    }
    istringstream is(option_value[0]);
    int val;
    if (is >> val){
      this->field = val;
      return "";
    }
    return badValue(option_value, "int", this->name);
  }
  void SetDefault(){
    this->field = this->def;
  }
};

class COptionULong : public COptionBase{
  unsigned long & field; // Reference to the feildname
  unsigned long def; // Default value
  string name; // identifier for the option

public:
  COptionULong(string option_field_name, unsigned long & option_field, unsigned long default_value) : field(option_field){
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionULong(){};
  string SetValue(vector<string> option_value){
    string out = optionCheckMultipleValues(option_value, "unsigned long", this->name);
    if (out.compare("") != 0){
      return out;
    }
    istringstream is(option_value[0]);
    unsigned long val;
    if (is >> val){
      this->field = val;
      return "";
    }
    return badValue(option_value, "unsigned long", this->name);
  }
  void SetDefault(){
    this->field = this->def;
  }
};

class COptionUShort : public COptionBase{
  unsigned short & field; // Reference to the feildname
  unsigned short def; // Default value
  string name; // identifier for the option

public:
  COptionUShort(string option_field_name, unsigned short & option_field, unsigned short default_value) : field(option_field){
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionUShort(){};
  string SetValue(vector<string> option_value){
    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0){
      return out;
    }
    istringstream is(option_value[0]);
    unsigned short val;
    if (is >> val){
      this->field = val;
      return "";
    }
    return badValue(option_value, "unsigned short", this->name);
  }
  void SetDefault(){
    this->field = this->def;
  }
};

class COptionLong : public COptionBase{
  long & field; // Reference to the feildname
  long def; // Default value
  string name; // identifier for the option

public:
  COptionLong(string option_field_name, long & option_field, long default_value) : field(option_field){
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionLong(){};
  string SetValue(vector<string> option_value){
    string out = optionCheckMultipleValues(option_value, "long", this->name);
    if (out.compare("") != 0){
      return out;
    }
    istringstream is(option_value[0]);
    long val;
    if (is >> val){
      this->field = val;
      return "";
    }
    return badValue(option_value, "long", this->name);
  }
  void SetDefault(){
    this->field = this->def;
  }
};


class COptionBool : public COptionBase{
  bool & field; // Reference to the feildname
  bool def; // Default value
  string name; // identifier for the option

public:
  COptionBool(string option_field_name, bool & option_field, bool default_value) : field(option_field){
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionBool(){};
  string SetValue(vector<string> option_value){
    // check if there is more than one value
    string out = optionCheckMultipleValues(option_value, "bool", this->name);
    if (out.compare("") != 0){
      return out;
    }
    if (option_value[0].compare("YES") == 0){
      this->field = true;
      return "";
    }
    if (option_value[0].compare("NO") == 0){
      this->field = false;
      return "";
    }
    return badValue(option_value, "bool", this->name);
  }
  void SetDefault(){
    this->field = this->def;
  }
};

template <class Tenum>
class COptionEnumList : public COptionBase{

  map<string, Tenum> m;
  unsigned short * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionEnumList(string option_field_name, const map<string, Tenum> m, unsigned short * & option_field, unsigned short & list_size) : field(option_field) , size(list_size){
    this->m = m;
    this->name = option_field_name;
  }

  ~COptionEnumList(){};
  string SetValue(vector<string> option_value){
    if (option_value.size() == 1 && option_value[0].compare("NONE")==0){
      this->size = 0;
      return "";
    }
    // size is the length of the option list
    this->size = option_value.size();
    unsigned short * enums = new unsigned short[size];
    for(int i  = 0; i < this->size; i++){
      // Check to see if the enum value is in the map
      if (this->m.find(option_value[i]) == m.end()){
        string str;
        str.append(this->name);
        str.append(": invalid option value ");
        str.append(option_value[0]);
        str.append(". Check current SU2 options in config_template.cfg.");
        return str;
      }
      // If it is there, set the option value
      enums[i] = this->m[option_value[0]];
    }
    this->field = enums;
    return "";
  }

  void SetDefault(){
    // No default to set
    size = 0;
  }
};

class COptionDoubleArray : public COptionBase{
  double * & field; // Reference to the feildname
  string name; // identifier for the option
  const int size;
  double * default_value;

public:
  COptionDoubleArray(string option_field_name, const int list_size, double * & option_field, double * default_value) : field(option_field), size(list_size){
    this->name = option_field_name;
    this->default_value = default_value;
  }

  ~COptionDoubleArray(){};
  string SetValue(vector<string> option_value){
    // Check that the size is correct
    if (option_value.size() != this->size){
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
    double * vals = new double[this->size];
    for(int i  = 0; i < this->size; i++){
      istringstream is(option_value[i]);
      double val;
      if (!(is >> val)){
        delete [] vals;
        return badValue(option_value, "double array", this->name);
      }
      vals[i] = val;
    }
    this->field = vals;
    return "";
  }

  void SetDefault(){
    this->field = this->default_value;
  }
};

class COptionDoubleList : public COptionBase{
  double * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionDoubleList(string option_field_name, unsigned short & list_size, double * & option_field) : field(option_field), size(list_size){
    this->name = option_field_name;
  }

  ~COptionDoubleList(){};
  string SetValue(vector<string> option_value){
    // The size is the length of option_value
    unsigned long option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE")==0){
      // No options
      this->size = 0;
      return "";
    }

    this->size = option_size;

    // Parse all of the options
    double * vals = new double[option_size];
    for(int i  = 0; i < option_size; i++){
      istringstream is(option_value[i]);
      double val;
      if (!(is >> val)){
        delete [] vals;
        return badValue(option_value, "double list", this->name);
      }
      vals[i] = val;
    }
    this->field = vals;
    return "";
  }

  void SetDefault(){
    this->size = 0; // There is no default value for list
  }
};

class COptionUShortList : public COptionBase{
  unsigned short * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionUShortList(string option_field_name, unsigned short & list_size, unsigned short * & option_field) : field(option_field), size(list_size){
    this->name = option_field_name;
  }

  ~COptionUShortList(){};
  string SetValue(vector<string> option_value){
    // The size is the length of option_value
    unsigned long option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE")==0){
      // No options
      this->size = 0;
      return "";
    }
    this->size = option_size;

    // Parse all of the options
    unsigned short * vals = new unsigned short[option_size];
    for(int i  = 0; i < option_size; i++){
      istringstream is(option_value[i]);
      unsigned short val;
      if (!(is >> val)){
        delete [] vals;
        return badValue(option_value, "unsigned short", this->name);
      }
      vals[i] = val;
    }
    this->field = vals;
    return "";
  }

  void SetDefault(){
    this->size = 0; // There is no default value for list
  }
};

class COptionStringList : public COptionBase{
  string * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionStringList(string option_field_name, unsigned short & list_size, string * & option_field) : field(option_field), size(list_size){
    this->name = option_field_name;
  }

  ~COptionStringList(){};
  string SetValue(vector<string> option_value){
    // The size is the length of option_value
    unsigned long option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE")==0){
      this->size = 0;
      return "";
    }
    this->size = option_size;

    // Parse all of the options
    string * vals = new string[option_size];
    for(int i  = 0; i < option_size; i++){
      vals[i].assign(option_value[i]);
    }
    this->field = vals;
    return "";
  }

  void SetDefault(){
    this->size = 0; // There is no default value for list
  }
};

class COptionConvect : public COptionBase{
  string name; // identifier for the option
  unsigned short & space;
  unsigned short & centered;
  unsigned short & upwind;

public:
  COptionConvect(string option_field_name, unsigned short & space_field, unsigned short & centered_field, unsigned short & upwind_field) : space(space_field), centered(centered_field), upwind(upwind_field){
    this->name = option_field_name;
  }

  ~COptionConvect(){};
  string SetValue(vector<string> option_value){

    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0){
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

  void SetDefault(){
    this->centered = NO_CENTERED;
    this->upwind = NO_UPWIND;
    this->space = SPACE_CENTERED;
  }
};

class COptionMathProblem : public COptionBase{
  string name; // identifier for the option
  bool & adjoint;
  bool & linearized;
  bool & restart;
  bool adjoint_def;
  bool linearized_def;
  bool restart_def;

public:
  COptionMathProblem(string option_field_name, bool & adjoint_field, bool adjoint_default, bool & linearized_field, bool linearized_default, bool & restart_field, bool restart_default) : adjoint(adjoint_field), linearized(linearized_field), restart(restart_field) {
    this->name = option_field_name;
    this->adjoint_def = adjoint_default;
    this->linearized_def = linearized_default;
    this->restart_def = restart_default;
  }

  ~COptionMathProblem(){};
  string SetValue(vector<string> option_value){
    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0){
      return out;
    }
    if (Math_Problem_Map.find(option_value[0]) == Math_Problem_Map.end()){
      return badValue(option_value, "math problem", this->name);
    }
    if (option_value[0] == "DIRECT") {
      this->adjoint = false;
      this->linearized = false;
      this->restart = false;
      return "";
    }
    if (option_value[0] == "ADJOINT") {
      this->adjoint= true;
      this->restart= true;
      this->linearized = false;
      return "";
    }
    if (option_value[0] == "LINEARIZED") {
      this->linearized = true;
      this->restart = true;
      this->adjoint= false;
      return "";
    }
    return "option in math problem map not considered in constructor";
  }

  void SetDefault(){
    this->adjoint = this->adjoint_def;
    this->linearized = this->linearized_def;
    this->restart = this->restart_def;
  }
};

class COptionDVParam : public COptionBase{
  string name; // identifier for the option
  unsigned short & nDV;
  double ** & paramDV;
  string * & FFDTag;
  unsigned short* & design_variable;

public:
  COptionDVParam(string option_field_name, unsigned short & nDV_field, double** & paramDV_field, string* & FFDTag_field, unsigned short * & design_variable_field) : nDV(nDV_field), paramDV(paramDV_field), FFDTag(FFDTag_field), design_variable(design_variable_field){
    this->name = option_field_name;
  }

  ~COptionDVParam(){};
  string SetValue(vector<string> option_value){
    if ((option_value.size() == 1) && (option_value[0].compare("NONE") == 0)){
      this->nDV = 0;
      return "";
    }

    // Cannot have ; at the beginning or the end
    if (option_value[0].compare(";") == 0){
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have beginning semicolon");
      return newstring;
    }
    if (option_value[option_value.size()-1].compare(";") == 0){
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

    this->paramDV = new double*[this->nDV];
    for (unsigned short iDV = 0; iDV < this->nDV; iDV++){
      this->paramDV[iDV] = new double[MAX_PARAMETERS];
    }

    this->FFDTag = new string[this->nDV];

    unsigned short nParamDV = 0;
    stringstream ss;
    unsigned int i = 0;
    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) {
      switch (this->design_variable[iDV]) {
        case FFD_SETTING: nParamDV = 0; break;
        case FFD_CONTROL_POINT_2D: nParamDV = 5; break;
        case FFD_RADIUS_2D: nParamDV = 1; break;
        case FFD_CAMBER_2D: nParamDV = 2; break;
        case FFD_THICKNESS_2D: nParamDV = 2; break;
        case HICKS_HENNE: nParamDV = 2; break;
        case SPHERICAL: nParamDV = 3; break;
        case COSINE_BUMP: nParamDV = 3; break;
        case FOURIER: nParamDV = 3; break;
        case DISPLACEMENT: nParamDV = 3; break;
        case ROTATION: nParamDV = 6; break;
        case NACA_4DIGITS: nParamDV = 3; break;
        case PARABOLIC: nParamDV = 2; break;
        case OBSTACLE: nParamDV = 2; break;
        case AIRFOIL: nParamDV = 2; break;
        case STRETCH: nParamDV = 2; break;
        case FFD_CONTROL_POINT: nParamDV = 7; break;
        case FFD_DIHEDRAL_ANGLE: nParamDV = 7; break;
        case FFD_TWIST_ANGLE: nParamDV = 7; break;
        case FFD_ROTATION: nParamDV = 7; break;
        case FFD_CONTROL_SURFACE: nParamDV = 7; break;
        case FFD_CAMBER: nParamDV = 3; break;
        case FFD_THICKNESS: nParamDV = 3; break;
        case SURFACE_FILE: nParamDV = 0; break;
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
            ((this->design_variable[iDV] == FFD_SETTING) ||
             (this->design_variable[iDV] == FFD_CONTROL_POINT_2D) ||
             (this->design_variable[iDV] == FFD_CAMBER_2D) ||
             (this->design_variable[iDV] == FFD_THICKNESS_2D) ||
             (this->design_variable[iDV] == FFD_RADIUS_2D) ||
             (this->design_variable[iDV] == FFD_CONTROL_POINT) ||
             (this->design_variable[iDV] == FFD_DIHEDRAL_ANGLE) ||
             (this->design_variable[iDV] == FFD_TWIST_ANGLE) ||
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

  void SetDefault(){
    this->nDV = 0;
    this->paramDV = NULL;
    this->FFDTag = NULL;
    // Don't mess with the Design_Variable because it's an input, not modified
  }
};


// Class where the option is represented by (String, double, string, double, ...)
class COptionStringDoubleList : public COptionBase{
  string name; // identifier for the option
  unsigned short & size; // how many strings are there (same as number of doubles)

  string * & s_f; // Reference to the string fields
  double* & d_f; // reference to the double fields

public:
  COptionStringDoubleList(string option_field_name, unsigned short & list_size, string * & string_field, double* & double_field) : size(list_size), s_f(string_field), d_f(double_field){
    this->name = option_field_name;
  }

  ~COptionStringDoubleList(){};
  string SetValue(vector<string> option_value){
    // There must be an even number of entries (same number of strings and doubles
    unsigned long totalVals = option_value.size();
    if ((totalVals % 2) != 0){
      if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)){
        // It's okay to say its NONE
        this->size = 0;
        return "";
      }
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have an even number of entries");
      return newstring;
    }
    unsigned long nVals = totalVals / 2;
    this->size = nVals;
    this->s_f = new string[nVals];
    this->d_f = new double[nVals];

    for (int i = 0; i < nVals; i++){
      this->s_f[i].assign(option_value[2*i]); // 2 because have double and string
      istringstream is(option_value[2*i + 1]);
      double val;
      if (!(is >> val)){
        return badValue(option_value, "string double", this->name);
      }
      this->d_f[i] = val;
    }
    // Need to return something...
    return "";
  }

  void SetDefault(){
    this->size = 0; // There is no default value for list
  }
};

class COptionInlet : public COptionBase{
  string name; // identifier for the option
  unsigned short & size;
  string * & marker;
  double * & ttotal;
  double * & ptotal;
  double ** & flowdir;

public:
  COptionInlet(string option_field_name, unsigned short & nMarker_Inlet, string* & Marker_Inlet, double* & Ttotal, double* & Ptotal, double** & FlowDir) : size(nMarker_Inlet), marker(Marker_Inlet), ttotal(Ttotal), ptotal(Ptotal), flowdir(FlowDir){
    this->name = option_field_name;
  }

  ~COptionInlet(){};
  string SetValue(vector<string> option_value){

    unsigned long totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)){
      this->size = 0;
      this->marker = NULL;
      this->ttotal = NULL;
      this->ptotal = NULL;
      this->flowdir = NULL;
      return "";
    }

    if (totalVals % 6 != 0){
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

    unsigned long nVals = totalVals / 6;
    this->size = nVals;
    this->marker = new string[nVals];
    this->ttotal = new double[nVals];
    this->ptotal = new double[nVals];
    this->flowdir = new double*[nVals];
    for (int i = 0; i < nVals; i++){
      this->flowdir[i] = new double[3];
    }

    for (int i = 0; i < nVals; i++){
      this->marker[i].assign(option_value[6*i]);
      istringstream ss_1st(option_value[6*i + 1]);
      if(!(ss_1st >> this->ttotal[i])){
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_2nd(option_value[6*i + 2]);
      if(!(ss_2nd >> this->ptotal[i])){
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_3rd(option_value[6*i + 3]);
      if (!(ss_3rd >> this->flowdir[i][0])){
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_4th(option_value[6*i + 4]);
      if (!(ss_4th >> this->flowdir[i][1])){
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_5th(option_value[6*i + 5]);
      if (!(ss_5th >> this->flowdir[i][2])){
        return badValue(option_value, "inlet", this->name);
      }
    }

    return "";
  }

  void SetDefault(){
    this->marker = NULL;
    this->ttotal = NULL;
    this->ptotal = NULL;
    this->flowdir = NULL;
    this->size = 0; // There is no default value for list
  }
};

template <class Tenum>
class COptionRiemann : public COptionBase{

  map<string, Tenum> m;
  unsigned short* & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;
  string * & marker;
  double * & var1;
  double * & var2;
  double ** & flowdir;

public:
  COptionRiemann(string option_field_name, unsigned short & nMarker_Riemann, string* & Marker_Riemann, unsigned short* & option_field, const map<string, Tenum> m, double* & var1, double* & var2, double** & FlowDir) : size(nMarker_Riemann),
  	  	  	  	  marker(Marker_Riemann), field(option_field), var1(var1), var2(var2), flowdir(FlowDir){
    this->name = option_field_name;
    this->m = m;
  }
  ~COptionRiemann(){};

  string SetValue(vector<string> option_value){

    unsigned long totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)){
      this->size = 0;
      this->marker = NULL;
      this->field = 0;
      this->var1 = NULL;
      this->var2 = NULL;
      this->flowdir = NULL;
      return "";
    }

    if (totalVals % 7 != 0){
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

    unsigned long nVals = totalVals / 7;
    this->size = nVals;
    this->marker = new string[nVals];
    this->var1 = new double[nVals];
    this->var2 = new double[nVals];
    this->flowdir = new double*[nVals];
    this->field = new unsigned short[nVals];

    for (int i = 0; i < nVals; i++){
      this->flowdir[i] = new double[3];
    }

    for (int i = 0; i < nVals; i++){
      this->marker[i].assign(option_value[7*i]);
        // Check to see if the enum value is in the map
    if (this->m.find(option_value[7*i + 1]) == m.end()){
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
      if(!(ss_1st >> this->var1[i])){
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_2nd(option_value[7*i + 3]);
      if(!(ss_2nd >> this->var2[i])){
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_3rd(option_value[7*i + 4]);
      if (!(ss_3rd >> this->flowdir[i][0])){
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_4th(option_value[7*i + 5]);
      if (!(ss_4th >> this->flowdir[i][1])){
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_5th(option_value[7*i + 6]);
      if (!(ss_5th >> this->flowdir[i][2])){
        return badValue(option_value, "Riemann", this->name);
      }
    }

    return "";
  }

  void SetDefault(){
    this->marker = NULL;
    this->var1 = NULL;
    this->var2 = NULL;
    this->flowdir = NULL;
    this->size = 0; // There is no default value for list
  }
};

//Inlet condition where the input direction is assumed
class COptionExhaust : public COptionBase{
  string name; // identifier for the option
  unsigned short & size;
  string * & marker;
  double * & ttotal;
  double * & ptotal;

public:
  COptionExhaust(string option_field_name, unsigned short & nMarker_Exhaust, string* & Marker_Exhaust, double* & Ttotal, double* & Ptotal) : size(nMarker_Exhaust), marker(Marker_Exhaust), ttotal(Ttotal), ptotal(Ptotal){
    this->name = option_field_name;
  }

  ~COptionExhaust(){};
  
  string SetValue(vector<string> option_value){

    unsigned long totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)){
      this->size = 0;
      this->marker = NULL;
      this->ttotal = NULL;
      this->ptotal = NULL;
      return "";
    }

    if (totalVals % 3 != 0){
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 3");
      this->size = 0;
      this->marker = NULL;
      this->ttotal = NULL;
      this->ptotal = NULL;
      return newstring;
    }

    unsigned long nVals = totalVals / 3;
    this->size = nVals;
    this->marker = new string[nVals];
    this->ttotal = new double[nVals];
    this->ptotal = new double[nVals];

    for (int i = 0; i < nVals; i++){
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

  void SetDefault(){
    this->marker = NULL;
    this->ttotal = NULL;
    this->ptotal = NULL;
    this->size = 0; // There is no default value for list
  }
  
};

//Inlet condition where the input direction is assumed
class COptionBleed : public COptionBase{
  string name; // identifier for the option
  unsigned short & size;
  string * & marker;
  double * & massflow_target;
  double * & temp_target;
  
public:
  COptionBleed(string option_field_name, unsigned short & nMarker_Bleed, string* & Marker_Bleed, double* & MassFlow_Target, double* & Temp_Target) : size(nMarker_Bleed), marker(Marker_Bleed), massflow_target(MassFlow_Target), temp_target(Temp_Target){
    this->name = option_field_name;
  }
  
  ~COptionBleed(){};
  
  string SetValue(vector<string> option_value){
    
    unsigned long totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)){
      this->size = 0;
      this->marker = NULL;
      this->massflow_target = NULL;
      this->temp_target = NULL;
      return "";
    }
    
    if (totalVals % 3 != 0){
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 3");
      this->size = 0;
      this->marker = NULL;
      this->massflow_target = NULL;
      this->temp_target = NULL;
      return newstring;
    }
    
    unsigned long nVals = totalVals / 3;
    this->size = nVals;
    this->marker = new string[nVals];
    this->massflow_target = new double[nVals];
    this->temp_target = new double[nVals];
    
    for (int i = 0; i < nVals; i++){
      this->marker[i].assign(option_value[3*i]);
      istringstream ss_1st(option_value[3*i + 1]);
      if (!(ss_1st >> this->massflow_target[i]))
        return badValue(option_value, "bleed fixed", this->name);
      istringstream ss_2nd(option_value[3*i + 2]);
      if (!(ss_2nd >> this->temp_target[i]))
        return badValue(option_value, "bleed fixed", this->name);
    }
    
    return "";
  }
  
  void SetDefault(){
    this->marker = NULL;
    this->massflow_target = NULL;
    this->temp_target = NULL;
    this->size = 0; // There is no default value for list
  }
  
};

class COptionPeriodic : public COptionBase{
  string name; // identifier for the option
  unsigned short & size;
  string * & marker_bound;
  string * & marker_donor;
  double ** & rot_center;
  double ** & rot_angles;
  double ** & translation;

public:
  COptionPeriodic(const string option_field_name, unsigned short & nMarker_PerBound,
                  string* & Marker_PerBound, string* & Marker_PerDonor,
                  double** & RotCenter, double** & RotAngles, double** & Translation) : size(nMarker_PerBound), marker_bound(Marker_PerBound), marker_donor(Marker_PerDonor), rot_center(RotCenter), rot_angles(RotAngles), translation(Translation){
    this->name = option_field_name;
  }

  ~COptionPeriodic(){};
  string SetValue(vector<string> option_value){

    const int mod_num = 11;

    unsigned long totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)){
      this->size = 0;
      this->marker_bound = NULL;
      this->marker_donor = NULL;
      this->rot_center = NULL;
      this->rot_angles = NULL;
      this->translation = NULL;
      return "";
    }

    if (totalVals % mod_num != 0){
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

    unsigned long nVals = 2 * (totalVals / mod_num); // To account for periodic and donor
    this->size = nVals;
    this->marker_bound = new string[nVals];
    this->marker_donor = new string[nVals];
    this->rot_center = new double*[nVals];
    this->rot_angles = new double*[nVals];
    this->translation = new double*[nVals];
    for (int i = 0; i < nVals; i++){
      this->rot_center[i] = new double[3];
      this->rot_angles[i] = new double[3];
      this->translation[i] = new double[3];
    }

    double deg2rad = PI_NUMBER/180.0;

    for (int i = 0; i < (nVals/2); i++){
      this->marker_bound[i].assign(option_value[mod_num*i]);
      this->marker_donor[i].assign(option_value[mod_num*i+1]);
      istringstream ss_1st(option_value[mod_num*i + 2]);
      if (!(ss_1st >> this->rot_center[i][0])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_2nd(option_value[mod_num*i + 3]);
      if (!(ss_2nd >> this->rot_center[i][1])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_3rd(option_value[mod_num*i + 4]);
      if (!(ss_3rd >> this->rot_center[i][2])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_4th(option_value[mod_num*i + 5]);
      if (!(ss_4th >> this->rot_angles[i][0])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_5th(option_value[mod_num*i + 6]);
      if (!(ss_5th >> this->rot_angles[i][1])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_6th(option_value[mod_num*i + 7]);
      if (!(ss_6th >> this->rot_angles[i][2])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_7th(option_value[mod_num*i + 8]);
      if (!(ss_7th >> this->translation[i][0])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_8th(option_value[mod_num*i + 9]);
      if (!(ss_8th >> this->translation[i][1])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_9th(option_value[mod_num*i + 10]);
      if (!(ss_9th >> this->translation[i][2])){
        return badValue(option_value, "periodic", this->name);
      }
      this->rot_angles[i][0] *= deg2rad;
      this->rot_angles[i][1] *= deg2rad;
      this->rot_angles[i][2] *= deg2rad;
    }

    for (unsigned long i = (nVals/2); i < nVals; i++){
      this->marker_bound[i].assign(option_value[mod_num*(i-nVals/2)+1]);
      this->marker_donor[i].assign(option_value[mod_num*(i-nVals/2)]);
      istringstream ss_1st(option_value[mod_num*(i-nVals/2) + 2]);
      if (!(ss_1st >> this->rot_center[i][0])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_2nd(option_value[mod_num*(i-nVals/2) + 3]);
      if (!(ss_2nd >> this->rot_center[i][1])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_3rd(option_value[mod_num*(i-nVals/2) + 4]);
      if (!(ss_3rd >> this->rot_center[i][2])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_4th(option_value[mod_num*(i-nVals/2) + 5]);
      if (!(ss_4th >> this->rot_angles[i][0])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_5th(option_value[mod_num*(i-nVals/2) + 6]);
      if (!(ss_5th >> this->rot_angles[i][1])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_6th(option_value[mod_num*(i-nVals/2) + 7]);
      if (!(ss_6th >> this->rot_angles[i][2])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_7th(option_value[mod_num*(i-nVals/2) + 8]);
      if (!(ss_7th >> this->translation[i][0])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_8th(option_value[mod_num*(i-nVals/2) + 9]);
      if (!(ss_8th >> this->translation[i][1])){
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_9th(option_value[mod_num*(i-nVals/2) + 10]);
      if (!(ss_9th >> this->translation[i][2])){
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

  void SetDefault(){
    this->size = 0;
    this->marker_bound = NULL;
    this->marker_donor = NULL;
    this->rot_center = NULL;
    this->rot_angles = NULL;
    this->translation = NULL;
  }
};

class COptionPython : public COptionBase{
  string name;
public:
  COptionPython(const string name){
    this->name = name;
  }
  ~COptionPython(){};
  // No checking happens with python options
  string SetValue(vector<string> option_value){
    return "";
  }
  // No defaults with python options
  void SetDefault(){
    return;
  };
};

class COptionActuatorDisk : public COptionBase{
  string name; // identifier for the option
  unsigned short & inlet_size;
  unsigned short & outlet_size;
  string * & marker_inlet;
  string * & marker_outlet;
  double ** & origin;
  double * & root_radius;
  double * & tip_radius;
  double * & press_jump;
  double * & temp_jump;
  double * & omega;

public:
  COptionActuatorDisk(const string name, unsigned short & nMarker_ActDisk_Inlet, unsigned short & nMarker_ActDisk_Outlet, string * & Marker_ActDisk_Inlet, string * & Marker_ActDisk_Outlet, double ** & ActDisk_Origin, double * & ActDisk_RootRadius, double * & ActDisk_TipRadius, double * & ActDisk_PressJump, double * & ActDisk_TempJump, double * & ActDisk_Omega) : inlet_size(nMarker_ActDisk_Inlet),outlet_size(nMarker_ActDisk_Outlet), marker_inlet(Marker_ActDisk_Inlet), marker_outlet(Marker_ActDisk_Outlet), origin(ActDisk_Origin), root_radius(ActDisk_RootRadius), tip_radius(ActDisk_TipRadius), press_jump(ActDisk_PressJump), temp_jump(ActDisk_TempJump), omega(ActDisk_Omega) {
    this->name = name;
  }

  ~COptionActuatorDisk(){};
  string SetValue(vector<string> option_value){
    const int mod_num = 10;
    unsigned long totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)){
      this->SetDefault();
      return "";
    }

    if (totalVals % mod_num != 0){
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 10");
      this->SetDefault();
      return newstring;
    }

    unsigned long nVals = totalVals / mod_num;
    this->inlet_size = nVals;
    this->outlet_size = nVals;
    this->marker_inlet = new string[this->inlet_size];
    this->marker_outlet = new string[this->outlet_size];
    this->root_radius = new double[this->inlet_size];
    this->tip_radius = new double[this->inlet_size];
    this->press_jump = new double[this->outlet_size];
    this->temp_jump = new double[this->outlet_size];
    this->omega = new double[this->inlet_size];

    this->origin = new double*[this->inlet_size];
    for (int i = 0; i < this->inlet_size; i++){
      this->origin[i] = new double[3];
    }

    string tname = "actuator disk";

    for (int i = 0; i < this->inlet_size; i++){
      this->marker_inlet[i].assign(option_value[mod_num*i]);
      this->marker_outlet[i].assign(option_value[mod_num*i+1]);
      istringstream ss_1st(option_value[mod_num*i + 2]);
      if (!(ss_1st >> this->origin[i][0])){
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_2nd(option_value[mod_num*i + 3]);
      if (!(ss_2nd >> this->origin[i][1])){
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_3rd(option_value[mod_num*i + 4]);
      if (!(ss_3rd >> this->origin[i][2])){
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_4th(option_value[mod_num*i + 5]);
      if (!(ss_4th >> this->root_radius[i])){
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_5th(option_value[mod_num*i + 6]);
      if (!(ss_5th >> this->tip_radius[i])){
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_6th(option_value[mod_num*i + 7]);
      if (!(ss_6th >> this->press_jump[i])){
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_7th(option_value[mod_num*i + 8]);
      if (!(ss_7th >> this->temp_jump[i])){
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_8th(option_value[mod_num*i + 9]);
      if (!(ss_8th >> this->omega[i])){
        return badValue(option_value, tname, this->name);
      }
    }
    return "";
  }
  void SetDefault(){
    this->inlet_size = 0;
    this->outlet_size = 0;
    this->marker_inlet = NULL;
    this->marker_outlet = NULL;
    this->origin = NULL;
    this->root_radius = NULL;
    this->tip_radius = NULL;
    this->press_jump = NULL;
    this->temp_jump = NULL;
    this->omega = NULL;
  }
};

