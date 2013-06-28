/*!
 * \file option_structure.hpp
 * \brief Defines classes for referencing options for easy input in CConfig
 * \author J. Hicken
 * \version 1.1.
 *
 * Many of the classes in this file are templated, and therefore must
 * be declared and defined here; to keep all elements together, there
 * is no corresponding .cpp file at this time.
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
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
	SU2_MDC = 2,	/*!< \brief Running the SU2_MDC software. */
	SU2_GPC = 3,	/*!< \brief Running the SU2_GPC software. */
	SU2_DDC = 4,	/*!< \brief Running the SU2_DDC software. */
	SU2_MAC = 5,	/*!< \brief Running the SU2_MAC software. */
	SU2_GDC = 6,	/*!< \brief Running the SU2_GDC software. */
	SU2_PBC = 7	  /*!< \brief Running the SU2_PBC software. */
};

const unsigned int MAX_PROCESSORS = 1000;	/*!< \brief Maximum number of processors. */
const unsigned int MAX_PARAMETERS = 10;		/*!< \brief Maximum number of parameters for a design variable definition. */
const unsigned int MAX_INDEX_VALUE = 100;	/*!< \brief Maximum value for a marker index. */
const unsigned int MAX_NUMBER_MARKER = 200;	/*!< \brief Maximum number of domains. */
const unsigned int MAX_NUMBER_CHUNCK = 10;	/*!< \brief Maximum number of chunks for the FFD. */
const unsigned int MAX_NEIGHBORS = 150;		/*!< \brief Maximum number of neighbors. */
const unsigned int MAX_SOLS = 6;		/*!< \brief Maximum number of solutions at the same time (dimension of solution container array). */
const unsigned int MAX_TERMS = 5;		/*!< \brief Maximum number of terms in the numerical equations (dimension of solver container array). */
const unsigned int NO_RK_ITER = 0;		/*!< \brief No Runge-Kutta iteration. */
const unsigned int MESH_0 = 0;			/*!< \brief Definition of the finest grid level. */
const unsigned int DOMAIN_0 = 0;			/*!< \brief Definition of the first grid domain. */
const unsigned int MAX_MPI_BUFFER = 3000000; /*!< \brief Buffer size for parallel simulations. */

const double PRANDTL = 0.72;	        	/*!< \brief Fluid's Prandtl constant (air). */
const double PRANDTL_TURB = 0.90;	/*!< \brief Fluid's turbulent Prandtl constant (air). */
const double AVOGAD_CONSTANT = 6.0221415E26;	/*!< \brief Avogardro's constant, number of particles in one kmole. */
const double BOLTZMANN_CONSTANT = 1.3806503E-23; /*! \brief Boltzmann's constant [J K^-1] */
const double UNIVERSAL_GAS_CONSTANT = 8314.462175; /*! \brief Universal gas constant [J kmol^-1 K^-1] */
const double ELECTRON_CHARGE = 1.60217646E-19;	/*!< \brief Electronic charge constant. */
const double ELECTRON_MASS = 9.10938188E-31;	/*!< \brief Mass of an electron. */
const double FREE_PERMITTIVITY = 8.8541878176E-12; /*!< \brief Premittivity of free space. */
const double MAGNETIC_CONSTANT = 1.25663706E-6;  /*!< \brief magnetic permeability of free space. */
const double STANDART_GRAVITY = 9.80665;        /*!< \brief Acceleration due to gravity at surface of earth. */
const double EPS = 1.0E-16;			/*!< \brief Error scale. */
const double TURB_EPS = 1.0E-16;		/*!< \brief Turbulent Error scale. */
const double ONE2 = 0.5;			/*!< \brief One divided by two. */
const double TWO3 = 2.0 / 3.0;			/*!< \brief Two divided by three. */
const double FOUR3 = 4.0 / 3.0;			/*!< \brief Four divided by three. */
const double PI_NUMBER = 4.0 * atan(1.0);	/*!< \brief Pi number. */

/*!
 * \brief Boolean answers
 */
enum ANSWER {
	NONE = 0,
	NO = 0,    /*!< \brief Boolean definition of no. */
	YES = 1	/*!< \brief Boolean definition of yes. */
};

const unsigned int MAX_NUMBER_DOMAIN = 100;	/*!< \brief Maximum number of domains. */
const unsigned int MAX_NUMBER_PERIODIC = 10;	/*!< \brief Maximum number of periodic boundary conditions. */
const int MASTER_NODE = 0;			/*!< \brief Master node for MPI parallelization. */
const int AUX_NODE = 1;			/*!< \brief Computational node that is used for IO stuff. */

/*!
 * \brief types of MPI communications
 */
enum COMM_TYPE {
	SEND = 1,					/*!< \brief Boolean definition of send (parallelization). */
	RECEIVE = 2				/*!< \brief Boolean definition of receive (parallelization). */
};

/*!
 * \brief different solver types for the CFD component
 */
enum ENUM_SOLVER {
	NO_SOLVER = 0,			/*!< \brief Definition of no solver. */
	EULER = 2,				/*!< \brief Definition of the Euler's solver. */
	NAVIER_STOKES = 3,			/*!< \brief Definition of the Navier-Stokes' solver. */
	RANS = 4,				/*!< \brief Definition of the Reynolds-averaged Navier-Stokes' (RANS) solver. */
	ELECTRIC_POTENTIAL = 5,       	/*!< \brief Definition of the electric potential solver. */
	ADJ_EULER = 7,			/*!< \brief Definition of the continuous adjoint Euler's solver. */
	ADJ_NAVIER_STOKES = 8,		/*!< \brief Definition of the continuous adjoint Navier-Stokes' solver. */
	ADJ_RANS = 9,				/*!< \brief Definition of the continuous adjoint Reynolds-averaged Navier-Stokes' (RANS) solver. */
	LIN_EULER = 12,			/*!< \brief Definition of the linear Euler's solver. */
	LIN_NAVIER_STOKES = 13,		/*!< \brief Definition of the linear Navier-Stokes' solver. */
	FREE_SURF_EULER = 18,			/*!< \brief Definition of the Free Surface Euler solver. */
	FREE_SURF_NAVIER_STOKES = 19,			/*!< \brief Definition of the Free Surface Navier-Stokes solver. */
	FREE_SURF_RANS = 20,			/*!< \brief Definition of the Free Surface RANS solver. */
	NS_PLASMA = 21,	/*!< \brief Definition of the plasma solver. */
	ADJ_FREE_SURF_EULER = 23,			/*!< \brief Definition of the adjoint Free Surface Euler solver. */
	ADJ_FREE_SURF_NAVIER_STOKES = 24,			/*!< \brief Definition of the adjoint Free Surface Navier-Stokes solver. */
	ADJ_FREE_SURF_RANS = 25,			/*!< \brief Definition of the adjoint Free Surface RANS solver. */
	ADJ_NS_PLASMA = 26,	/*!< \brief Definition of the adjoint plasma solver. */
	WAVE = 27,	/*!< \brief Definition of the wave solver. */
	TEMPLATE_SOLVER = 28                  /*!< \brief Definition of template solver. */
	
};
static const map<string, ENUM_SOLVER> Solver_Map = CCreateMap<string, ENUM_SOLVER>
("NONE", NO_SOLVER)
("EULER", EULER)
("NAVIER_STOKES", NAVIER_STOKES)
("RANS", RANS)
("ELECTRIC_POTENTIAL", ELECTRIC_POTENTIAL)
("ADJ_EULER", ADJ_EULER)
("ADJ_NAVIER_STOKES", ADJ_NAVIER_STOKES)
("ADJ_RANS", ADJ_RANS )
("LIN_EULER", LIN_EULER)
("LIN_NAVIER_STOKES", LIN_NAVIER_STOKES)
("NS_PLASMA", NS_PLASMA)
("FREE_SURF_EULER", FREE_SURF_EULER)
("FREE_SURF_NAVIER_STOKES", FREE_SURF_NAVIER_STOKES)
("FREE_SURF_RANS", FREE_SURF_RANS)
("WAVE", WAVE)
("TEMPLATE_SOLVER", TEMPLATE_SOLVER);

/*!
 * \brief different types of systems
 */
enum RUNTIME_TYPE {
	RUNTIME_POT_SYS = 1,			/*!< \brief One-physics case, the code is solving the potential equation. */
	RUNTIME_FLOW_SYS = 2,			/*!< \brief One-physics case, the code is solving the flow equations(Euler and Navier-Stokes). */
	RUNTIME_TURB_SYS = 3,			/*!< \brief One-physics case, the code is solving the turbulence model. */
	RUNTIME_ELEC_SYS = 4,			/*!< \brief One-physics case, the code is solving the electrical potential equation. */
	RUNTIME_PLASMA_SYS = 15,		/*!< \brief One-physics case, the code is solving the plasma equations. */
	RUNTIME_LEVELSET_SYS = 16,		/*!< \brief One-physics case, the code is solving the level set equations. */
	RUNTIME_WAVE_SYS = 8,		/*!< \brief One-physics case, the code is solving the wave equation. */
	RUNTIME_ADJPOT_SYS = 5,		/*!< \brief One-physics case, the code is solving the adjoint potential flow equation. */
	RUNTIME_ADJFLOW_SYS = 6,		/*!< \brief One-physics case, the code is solving the adjoint equations is being solved (Euler and Navier-Stokes). */
	RUNTIME_ADJTURB_SYS = 7,		/*!< \brief One-physics case, the code is solving the adjoint turbulence model. */
	RUNTIME_ADJLEVELSET_SYS = 18,		/*!< \brief One-physics case, the code is solving the adjoint evel set equations. */
	RUNTIME_LINPOT_SYS = 9,		/*!< \brief One-physics case, the code is solving the linear potential flow equations. */
	RUNTIME_LINFLOW_SYS = 10,		/*!< \brief One-physics case, the code is solving the linear equations is being solved (Euler and Navier-Stokes). */
	RUNTIME_FLOW_ADJFLOW_SYS = 12,	/*!< \brief Multi-physics case, we are solving together the flow and adjoint equation. */
	RUNTIME_FLOW_LINFLOW_SYS = 13,	/*!< \brief Multi-physics case, we are solving together the flow and linearized equation. */
	RUNTIME_MULTIGRID_SYS = 14,   	/*!< \brief Full Approximation Storage Multigrid system of equations. */
	RUNTIME_ADJPLASMA_SYS = 19		/*!< \brief One-physics case, the code is solving the plasma equations. */
};

const int FLOW_SOL = 0;		/*!< \brief Position of the flow solution (potential, Euler and Navier-Stokes) in the solution container array. */
const int TURB_SOL = 2;		/*!< \brief Position of the turbulence solution in the solution container array. */
const int ELEC_SOL = 2;		/*!< \brief Position of the electronic potential solution in the solution container array. */
const int WAVE_SOL = 1;		/*!< \brief Position of the wave equation in the solution container array. */
const int ADJFLOW_SOL = 1;	/*!< \brief Position of the continuous adjoint flow solution in the solution container array. */
const int ADJTURB_SOL = 3;	/*!< \brief Position of the continuous adjoint turbulence solution in the solution container array. */
const int ADJLEVELSET_SOL = 3;	/*!< \brief Position of the continuous adjoint level set solution in the solution container array. */
const int LINFLOW_SOL = 1;	/*!< \brief Position of the linear flow solution in the solution container array. */
const int PLASMA_SOL = 1;	/*!< \brief Position of the plasma solution in the solution container array. */
const int ADJPLASMA_SOL = 2;	/*!< \brief Position of the adjoint plasma solution in the solution container array. */
const int LEVELSET_SOL = 2;	/*!< \brief Position of the level set solution in the solution container array. */
const int TEMPLATE_SOL = 0;     /*!< \brief Position of the template solution. */
const int CONV_TERM = 0;	/*!< \brief Position of the convective terms in the solver container array. */
const int VISC_TERM = 1;        /*!< \brief Position of the viscous terms in the solver container array. */
const int SOUR_TERM = 2;        /*!< \brief Position of the source terms in the solver container array. */
const int CONS_SOUR_TERM = 3;   /*!< \brief Position of the conservative source terms in the solver container array. */
const int BOUND_TERM = 4;       /*!< \brief Position of the boundary terms in the solver container array. */

/*!
 * \brief types of spatial discretizations
 */
enum ENUM_SPACE {
	NO_CONVECTIVE = 0, /*!< \brief No convective scheme is used. */
	SPACE_CENTRED = 1,		/*!< \brief Space centred convective numerical method. */
	SPACE_UPWIND = 2		/*!< \brief Upwind convective numerical method. */
};
static const map<string, ENUM_SPACE> Space_Map = CCreateMap<string, ENUM_SPACE>
("NONE", NO_CONVECTIVE)
("SPACE_CENTRED", SPACE_CENTRED)
("SPACE_UPWIND", SPACE_UPWIND);

/*!
 * \brief types of spatial discretizations
 */
enum ENUM_GASMODEL {
	NO_MODEL = 0, /*!< \brief _____. */
	ARGON = 1,		/*!< \brief _____. */
	AIR7 = 2,		/*!< \brief _______. */
	AIR21 = 3		/*!< \brief _______. */
};
static const map<string, ENUM_GASMODEL> GasModel_Map = CCreateMap<string, ENUM_GASMODEL>
("NONE", NO_MODEL)       
("ARGON", ARGON)
("AIR-7", AIR7)
("AIR-21", AIR21);

/*!
 * \brief types of centered spatial discretizations
 */
enum ENUM_CENTRED {
	NO_CENTRED = 0,               /*!< \brief No centred scheme is used. */
	JST = 1,			/*!< \brief Jameson-Smith-Turkel centered numerical method. */
	LAX = 2			/*!< \brief Lax-Friedrich centered numerical method. */
};
static const map<string, ENUM_CENTRED> Centred_Map = CCreateMap<string, ENUM_CENTRED>
("NONE", NO_CENTRED)
("JST", JST)
("LAX-FRIEDRICH", LAX);

/*!
 * \brief types of upwind spatial discretizations
 */
enum ENUM_UPWIND {
	NO_UPWIND = 0,                /*!< \brief No upwind scheme is used. */
	ROE_1ST = 1,			/*!< \brief First order Roe's upwind numerical method. */
	ROE_2ND = 2,			/*!< \brief Second order Roe's upwind numerical method. */
	SCALAR_UPWIND_1ST = 3,	/*!< \brief First order scalar upwind numerical method. */
	SCALAR_UPWIND_2ND = 4,	/*!< \brief Second order scalar upwind numerical method. */
	CONVECTIVE_TEMPLATE = 5,       /*!< \brief Template for new numerical method . */
	AUSM_1ST = 6,			/*!< \brief First order AUSM numerical method. */
	AUSM_2ND = 7,			/*!< \brief Second order AUSM numerical method. */
	HLLC_1ST = 8,			/*!< \brief First order HLLC numerical method. */
	HLLC_2ND = 9,			/*!< \brief Second order HLLC numerical method. */
};
static const map<string, ENUM_UPWIND> Upwind_Map = CCreateMap<string, ENUM_UPWIND>
("NONE", NO_UPWIND)
("ROE-1ST_ORDER", ROE_1ST)
("ROE-2ND_ORDER", ROE_2ND)
("AUSM-1ST_ORDER", AUSM_1ST)
("AUSM-2ND_ORDER", AUSM_2ND)
("HLLC-1ST_ORDER", HLLC_1ST)
("HLLC-2ND_ORDER", HLLC_2ND)
("SCALAR_UPWIND-1ST_ORDER", SCALAR_UPWIND_1ST)
("SCALAR_UPWIND-2ND_ORDER", SCALAR_UPWIND_2ND)
("CONVECTIVE_TEMPLATE", CONVECTIVE_TEMPLATE);

/*!
 * \brief types of slope limiters
 */
enum ENUM_LIMITER {
	NO_LIMITER = 0,               /*!< \brief No slope limiter */
	VENKATAKRISHNAN = 1,		/*!< \brief Slope limiter using Venkatakrisnan method. */
	BARTH = 2							/*!< \brief Slope limiter using Barth and Jespersen method. */
};
static const map<string, ENUM_LIMITER> Limiter_Map = CCreateMap<string, ENUM_LIMITER>
("NONE", NO_LIMITER)
("VENKATAKRISHNAN", VENKATAKRISHNAN)
("BARTH", BARTH);

/*!
 * \brief types of viscous term discretizations
 */
enum ENUM_VISCOUS {
	NO_VISCOUS = 0,               /*!< \brief No viscous term computation. */
	AVG_GRAD = 1,			/*!< \brief Average of gradients method for viscous term computation. */
	AVG_GRAD_CORRECTED = 2,	/*!< \brief Average of gradients with correction for viscous term computation. */
	GALERKIN = 3			/*!< \brief Galerkin method for viscous term computation. */
};
static const map<string, ENUM_VISCOUS> Viscous_Map = CCreateMap<string, ENUM_VISCOUS>
("NONE", NO_VISCOUS)
("AVG_GRAD", AVG_GRAD)
("AVG_GRAD_CORRECTED", AVG_GRAD_CORRECTED)
("GALERKIN", GALERKIN);

/*!
 * \brief types of source term methods to use
 */
enum ENUM_SOURCE {
	NO_SOURCE = 0,                /*!< \brief No source term. */
	PIECEWISE_CONSTANT = 1,	/*!< \brief Numerical method for source term in flow equations. */
	CHARGE_DIST = 2,		/*!< \brief Numerical method for source term in charge distribution. */
	SOURCE_TEMPLATE = 4           /*!< \brief Template for New numerical method for source term. */
};
static const map<string, ENUM_SOURCE> Source_Map = CCreateMap<string, ENUM_SOURCE>
("NONE", NO_SOURCE)
("PIECEWISE_CONSTANT", PIECEWISE_CONSTANT)
("CHARGE_DIST", CHARGE_DIST)
("TEMPLATE_SOURCE_METHOD", SOURCE_TEMPLATE);

/*!
 * \brief types of turbulent models
 */
enum ENUM_TURB_MODEL {
	NO_TURB_MODEL = 0,            /*!< \brief No turbulence model. */
	SA = 1,                       /*!< \brief Kind of Turbulent model (Spalart-Allmaras). */
	SA_COMP = 2,          	/*!< \brief Kind of Turbulent model (Spalart-Allmaras) for compresibles flows. */
	SST = 3       		/*!< \brief Kind of Turbulence model (Menter SST). */
};
static const map<string, ENUM_TURB_MODEL> Turb_Model_Map = CCreateMap<string, ENUM_TURB_MODEL>
("NONE", NO_TURB_MODEL)
("SA", SA)
("SA_COMP", SA_COMP)
("SST", SST);

/*!
 * \brief type of time integration schemes
 */
enum ENUM_TIME_INT {
	RUNGE_KUTTA_EXPLICIT = 1,	/*!< \brief Explicit Runge-Kutta time integration definition. */
	EULER_EXPLICIT = 2,   	/*!< \brief Explicit Euler time integration definition. */
	EULER_IMPLICIT = 3,   	/*!< \brief Implicit Euler time integration definition. */
};
static const map<string, ENUM_TIME_INT> Time_Int_Map = CCreateMap<string, ENUM_TIME_INT>
("RUNGE-KUTTA_EXPLICIT", RUNGE_KUTTA_EXPLICIT)
("EULER_EXPLICIT", EULER_EXPLICIT)
("EULER_IMPLICIT", EULER_IMPLICIT);

/*!
 * \brief types of schemes to compute the gradient
 */
enum ENUM_GRADIENT {
	GREEN_GAUSS = 1,		/*!< \brief Gradients computation using Green Gauss theorem. */
	LEAST_SQUARES = 2,		/*!< \brief Gradients computation using Least Squares. */
	WEIGHTED_LEAST_SQUARES = 3	/*!< \brief Gradients computation using Weighted Least Squares. */
};
static const map<string, ENUM_GRADIENT> Gradient_Map = CCreateMap<string, ENUM_GRADIENT>
("GREEN_GAUSS", GREEN_GAUSS)
("LEAST_SQUARES", LEAST_SQUARES)
("WEIGHTED_LEAST_SQUARES", WEIGHTED_LEAST_SQUARES);

/*!
 * \brief types of action to take on a geometry structure
 */
enum GEOMETRY_ACTION {
	ALLOCATE = 0,         	/*!<  \brief Allocate geometry structure. */
	UPDATE = 1            	/*!<  \brief Update geometry structure (grid moving, adaptation, etc.). */
};

/*!
 * \brief types of boundary conditions
 */
enum BC_TYPE {
	EULER_WALL = 1,		/*!< \brief Boundary Euler wall definition. */
	NO_SLIP_WALL = 4,		/*!< \brief Boundary no slip wall definition. */
	FAR_FIELD = 2,		/*!< \brief Boundary far-field definition. */
	SYMMETRY_PLANE = 3,   	/*!< \brief Boundary symmetry plane definition. */
	INLET_FLOW = 5,		/*!< \brief Boundary inlet flow definition. */
	OUTLET_FLOW = 6,		/*!< \brief Boundary outlet flow definition. */
	PERIODIC_BOUNDARY = 7,	/*!< \brief Periodic boundary definition. */
	NEARFIELD_BOUNDARY = 8,	/*!< \brief Near-Field boundary definition. */
	ELECTRODE_BOUNDARY = 9,	/*!< \brief Electrode boundary definition. */
	DIELECTRIC_BOUNDARY = 10,	/*!< \brief Dielectric boundary definition. */
	CUSTOM_BOUNDARY = 11,         /*!< \brief custom boundary definition. */
	INTERFACE_BOUNDARY = 12,	/*!< \brief Domain interface boundary definition. */
	DIRICHLET = 13,		/*!< \brief Boundary Euler wall definition. */
	NEUMANN = 14,		/*!< \brief Boundary Neumann definition. */
	SEND_RECEIVE = 99		/*!< \brief Boundary send-receive definition. */
};

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
	DRAG_COEFFICIENT = 1, 	/*!< \brief Drag objective function definition. */
	LIFT_COEFFICIENT = 2, 	/*!< \brief Lift objective function definition. */
	SIDEFORCE_COEFFICIENT = 3,	/*!< \brief Side force objective function definition. */
	EFFICIENCY = 4,		/*!< \brief Efficiency objective function definition. */
	PRESSURE_COEFFICIENT = 5,	/*!< \brief Pressure objective function definition. */
	MOMENT_X_COEFFICIENT = 6,	/*!< \brief Pitching moment objective function definition. */
	MOMENT_Y_COEFFICIENT = 7,	/*!< \brief Rolling moment objective function definition. */
	MOMENT_Z_COEFFICIENT = 8,	/*!< \brief Yawing objective function definition. */
	EQUIVALENT_AREA = 9,		/*!< \brief Equivalent area objective function definition. */
	NEARFIELD_PRESSURE = 10,	/*!< \brief NearField Pressure objective function definition. */
	ELECTRIC_CHARGE = 11,		/*!< \brief Charge distribution objective function definition. */
  FORCE_X_COEFFICIENT = 12,	/*!< \brief X-direction force objective function definition. */
  FORCE_Y_COEFFICIENT = 13,	/*!< \brief Y-direction force objective function definition. */
  FORCE_Z_COEFFICIENT = 14,	/*!< \brief Z-direction force objective function definition. */
  THRUST_COEFFICIENT = 15,		/*!< \brief Thrust objective function definition. */
  TORQUE_COEFFICIENT = 16,		/*!< \brief Torque objective function definition. */
  FIGURE_OF_MERIT = 17,		/*!< \brief Rotor Figure of Merit objective function definition. */
  FREESURFACE = 18,		/*!< \brief Free Surface objective function definition. */

};
static const map<string, ENUM_OBJECTIVE> Objective_Map = CCreateMap<string, ENUM_OBJECTIVE>
("DRAG", DRAG_COEFFICIENT)
("LIFT", LIFT_COEFFICIENT)
("SIDEFORCE", SIDEFORCE_COEFFICIENT)
("EFFICIENCY", EFFICIENCY)
("PRESSURE", PRESSURE_COEFFICIENT)
("MOMENT_X", MOMENT_X_COEFFICIENT)
("MOMENT_Y", MOMENT_Y_COEFFICIENT)
("MOMENT_Z", MOMENT_Z_COEFFICIENT)
("EQUIVALENT_AREA", EQUIVALENT_AREA)
("NEARFIELD_PRESSURE", NEARFIELD_PRESSURE)
("ELECTRIC_CHARGE", ELECTRIC_CHARGE)
("FORCE_X", FORCE_X_COEFFICIENT)
("FORCE_Y", FORCE_Y_COEFFICIENT)
("FORCE_Z", FORCE_Z_COEFFICIENT)
("THRUST", THRUST_COEFFICIENT)
("TORQUE", TORQUE_COEFFICIENT)
("FIGURE_OF_MERIT", FIGURE_OF_MERIT)
("FREESURFACE", FREESURFACE);

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
	HORIZONTAL_PLANE = 13,	/*!< \brief Do a grid refinement on a near field limit (line (y coordinate), or a plane (z coordinate). */
	SMOOTHING = 14,		/*!< \brief Do a grid smoothing of the geometry. */
	SUPERSONIC_SHOCK = 15,	/*!< \brief Do a grid smoothing. */
	TWOPHASE = 16,			/*!< \brief Do a grid refinement on the free surface interphase. */
	CURVE_SURFACE = 17,	/*!< \brief Create a curve surface using the current grid (hexa). */
	DOUBLE_SURFACE = 18	/*!< \brief Create another surface using a previous defined surface. */
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
("HORIZONTAL_PLANE", HORIZONTAL_PLANE)
("SMOOTHING", SMOOTHING)
("SUPERSONIC_SHOCK", SUPERSONIC_SHOCK)
("TWOPHASE", TWOPHASE)
("DOUBLE_SURFACE", DOUBLE_SURFACE)
("CURVE_SURFACE", CURVE_SURFACE);

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
	PARAVIEW = 1, 		/*!< \brief Paraview format for the solution output. */
	TECPLOT = 2,  		/*!< \brief Tecplot format for the solution output. */
	EXCEL = 3,			/*!< \brief Excel format for the solution output. */
	CSV = 4,			/*!< \brief Comma-separated values format for the solution output. */
	STL = 5				/*!< \brief STL CAD format for the solution output. */
};
static const map<string, ENUM_OUTPUT> Output_Map = CCreateMap<string, ENUM_OUTPUT>
("PARAVIEW", PARAVIEW)
("TECPLOT", TECPLOT)
("EXCEL", EXCEL)
("CSV", CSV)
("STL", STL);

/*!
 * \brief types of design parameterizations
 */
enum ENUM_PARAM {
  NO_DEFORMATION = 0,		/*!< \brief No surface deformation. */
  HICKS_HENNE = 1,		/*!< \brief Hicks-Henne bump function for airfoil deformation. */
  HICKS_HENNE_SHOCK = 2,	/*!< \brief Special deformation at the shocks as design variable. */
  HICKS_HENNE_NORMAL = 4,	/*!< \brief Normal Hicks-Henne bump function for airfoil deformation. */
  MACH_NUMBER = 5,		/*!< \brief Mach number as design variable. */
  NACA_4DIGITS = 6,		/*!< \brief The four digits NACA airfoil family as design variables. */
  DISPLACEMENT = 8,		/*!< \brief Surface movement as design variable. */
  ROTATION = 9,			/*!< \brief Surface rotation as design variable. */
  FFD_CONTROL_POINT = 10,	/*!< \brief Free form deformation for 3D design (change a control point). */
  FFD_DIHEDRAL_ANGLE = 11,	/*!< \brief Free form deformation for 3D design (change the dihedral angle). */
  FFD_TWIST_ANGLE = 12,		/*!< \brief Free form deformation for 3D design (change the twist angle). */
  FFD_ROTATION = 13,		/*!< \brief Free form deformation for 3D design (rotation around a line). */
  FFD_CAMBER = 14,		/*!< \brief Free form deformation for 3D design (camber change). */
  FFD_THICKNESS = 15,		/*!< \brief Free form deformation for 3D design (thickness change). */
  FFD_VOLUME = 16,		/*!< \brief Free form deformation for 3D design (volume change). */
  PARABOLIC = 17,		/*!< \brief Parabolic airfoil definition as design variables. */
  OBSTACLE = 18,		        /*!< \brief Obstacle for free surface optimization. */
	STRETCH = 19		        /*!< \brief Stretch one side of a channel. */
};
static const map<string, ENUM_PARAM> Param_Map = CCreateMap<string, ENUM_PARAM>
    ("NO_DEFORMATION", NO_DEFORMATION)
    ("HICKS_HENNE", HICKS_HENNE)
    ("HICKS_HENNE_SHOCK", HICKS_HENNE_SHOCK)
    ("HICKS_HENNE_NORMAL", HICKS_HENNE_NORMAL)
    ("MACH_NUMBER", MACH_NUMBER)
    ("NACA_4DIGITS", NACA_4DIGITS)
    ("DISPLACEMENT", DISPLACEMENT)
    ("ROTATION", ROTATION)
    ("FFD_CONTROL_POINT", FFD_CONTROL_POINT)
    ("FFD_DIHEDRAL_ANGLE", FFD_DIHEDRAL_ANGLE)
    ("FFD_TWIST_ANGLE", FFD_TWIST_ANGLE)
    ("FFD_ROTATION", FFD_ROTATION)
    ("FFD_CAMBER", FFD_CAMBER)
    ("FFD_THICKNESS", FFD_THICKNESS)
    ("FFD_VOLUME", FFD_VOLUME)
    ("PARABOLIC", PARABOLIC)
    ("OBSTACLE", OBSTACLE)
		("STRETCH", STRETCH);
    
/*!
 * \brief types of solvers for the grid deformation and free-form point inversion
 */
enum ENUM_LINEAR_SOLVER {
	STEEPEST_DESCENT = 1,		/*!< \brief Steepest descent method for point inversion algoritm (Free-Form). */
	NEWTON = 2,			/*!< \brief Newton method for point inversion algorithm (Free-Form). */
	QUASI_NEWTON = 3,		/*!< \brief Quasi Newton method for point inversion algorithm (Free-Form). */
	SYM_GAUSS_SEIDEL = 4,		/*!< \brief Symmetric Gauss-Seidel method for grid deformation. */
	CONJUGATE_GRADIENT = 5,	/*!< \brief Preconditionated conjugate gradient method for grid deformation. */
	FLEXIBLE_GMRES = 6,    	/*!< \brief Flexible Generalized Minimal RESidual method. */
	LU_SGS = 7,		/*!< \brief LU - Symmetric Gauss-Seidel method (main solver). */
	BCGSTAB = 8	/*!< \brief BCGSTAB - Biconjugate Gradient Stabilized Method (main solver). */
};
static const map<string, ENUM_LINEAR_SOLVER> Linear_Solver_Map = CCreateMap<string, ENUM_LINEAR_SOLVER>
("STEEPEST_DESCENT", STEEPEST_DESCENT)
("NEWTON", NEWTON)
("QUASI_NEWTON", QUASI_NEWTON)
("SYM_GAUSS_SEIDEL", SYM_GAUSS_SEIDEL)
("LU_SGS", LU_SGS)
("CONJUGATE_GRADIENT", CONJUGATE_GRADIENT)
("BCGSTAB", BCGSTAB)
("FLEXIBLE_GMRES", FLEXIBLE_GMRES);

/*!
 * \brief types of grid deformation techniques
 */
enum ENUM_DEFORM {
	SPRING = 1,  	         	/*!< \brief Classical spring analogy as the grid deformation technique. */
	TORSIONAL_SPRING = 2, 	/*!< \brief Torsional spring analogy as the grid deformation technique. */
	ALGEBRAIC = 3                 /*!< \brief Movement of the grid using an algebraic based method. */
};
static const map<string, ENUM_DEFORM> Deform_Map = CCreateMap<string, ENUM_DEFORM>
("SPRING", SPRING)
("TORSIONAL_SPRING", TORSIONAL_SPRING)
("ALGEBRAIC", ALGEBRAIC);

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
 * \brief types of schemes for unsteady computations
 */
enum ENUM_UNSTEADY {
	STEADY = 0,                   /*!< \brief A steady computation. */
	TIME_STEPPING = 1,		/*!< \brief Use a time stepping strategy for unsteady computations. */
	DT_STEPPING_1ST = 2,	/*!< \brief Use a dual time stepping strategy for unsteady computations (1st order). */
	DT_STEPPING_2ND = 3,	/*!< \brief Use a dual time stepping strategy for unsteady computations (2nd order). */
	ROTATIONAL_FRAME = 4    	/*!< \brief Use a rotational source term. */
};
static const map<string, ENUM_UNSTEADY> Unsteady_Map = CCreateMap<string, ENUM_UNSTEADY>
("NO", STEADY)
("TIME_STEPPING", TIME_STEPPING)
("DUAL_TIME_STEPPING-1ST_ORDER", DT_STEPPING_1ST)
("DUAL_TIME_STEPPING-2ND_ORDER", DT_STEPPING_2ND)
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
 * \class CAnyOptionRef
 * \brief provides a means of referencing variables of any type
 * \author J. Hicken
 *
 * In order to build a map that associates option names (strings) with
 * options, we need a way of referencing artbitrary types; otherwise,
 * we would need a separate map for ints, doubles, etc.  This class is
 * an abstract base class designed to accommodate this requirement.
 */
class CAnyOptionRef {
public:
	virtual ~CAnyOptionRef() = 0;
	virtual void WriteValue() = 0;
	virtual void SetValue(const vector<string> & value) = 0;
};
inline CAnyOptionRef::~CAnyOptionRef() {};

/*!
 * \class COptionRef
 * \brief a typed version of the base class for standard types
 * \tparam T - an arbitary standard type (short, int, double, etc)
 * \author J. Hicken
 *
 * This class can accommodate scalars and arrays of constant length.
 * If your option requies a variable length, use class CListOptionRef.
 * Why not use CListOptionRef for arrays?  You could, but COptionRef
 * lets you set default values for arrays (because it does not do any
 * memory management).
 */
template <class T>
class COptionRef : public CAnyOptionRef {
private:
	unsigned short ndim_; /*!< \brief number of array dimensions */
	unsigned int* dim_; /*!< \brief length of each array dimension */
	T* ref_; /*!< \brief pointer to the option */

	/*!< \brief function pointer used (if necessary) to set the value of the option */
	void (*user_set_value_)(T* ref, const vector<string> & value);

public:

	/*!
	 * \brief constructor for scalar options
	 * \param[in] value - variable to create a reference to
	 */
	COptionRef(T & value) {
		ndim_ = 0;
		dim_ = new unsigned int[ndim_];
		ref_ = &value;
		user_set_value_ = NULL;
	}

	/*!
	 * \brief constructor for 1D fixed length array options
	 * \param[in] value_ptr - pointer to variable to create a reference to
	 * \param[in] size - length of the array option
	 */
	COptionRef(T * & value_ptr, const int & size) {
		ndim_ = 1;
		dim_ = new unsigned int[ndim_];
		if (size <= 0) {
			cerr << "COptionRef::COptionRef(T&, const int &): "
					<< "invalid input: size = " << size << endl;
			throw(-1);
		}
		dim_[0] = size;
		ref_ = value_ptr;
		user_set_value_ = NULL;
	}

	/*!
	 * \brief constructor for scalar options that require special parsing
	 * \param[in] value - variable to create a reference to
	 * \param[in] set_value - function that sets options based on vector of strings
	 *
	 * If the option requires a special parsing of the configuration
	 * file, this constructor can be used to set a pointer to a function
	 * that performs the necessary parsing
	 */
	COptionRef(T & value, void (*set_value)(T*, const vector<string>&)) {
		ndim_ = 0;
		dim_ = new unsigned int[ndim_];
		ref_ = &value;
		user_set_value_ = set_value;
	}

	/*!
	 * \brief class destructor
	 */
	~COptionRef() {
		if (ndim_ > 0) {
			delete [] dim_;
		}
		user_set_value_ = NULL;
	}

	/*!
	 * \brief sets the value of the referenced option using vector of strings
	 * \param[in] value - a set of strings used to define the option
	 */
	void SetValue(const vector<string> & value) {
		stringstream ss;
		if (user_set_value_ == NULL) {
			if (ndim_ == 0) {
				ss << value[0];
				ss >> (*ref_);
			} else if (ndim_ == 1) {
				if (value.size() != dim_[0]) {
					cerr << "COptionRef::SetValue(const vector<string>&): "
							<< "number of input values does not match size of option array."
							<< endl;
					throw(-1);
				}
				for (unsigned int i = 0; i < dim_[0]; i++) {
					ss << value[i] << " ";
					ss >> (ref_[i]);
				}
			} else {
				cerr << "COptionRef::SetValue(const vector<string>&): "
						<< "cannot handle 2d arrays yet" << endl;
				throw(-1);
			}
		} else {
			// this option requires a special function for parsing
			(*user_set_value_)(ref_, value);
		}
	}

	/*!
	 * \brief write the value of the option to std out (mostly for debugging)
	 */
	void WriteValue() {
		if (ndim_ == 0) {
			cout << *ref_ << endl;
		} else if (ndim_ == 1) {
			for (unsigned int i = 0; i < dim_[0]; i++)
				cout << ref_[i] << ", ";
			cout << endl;
		}
	}
};

/*!
 * \class CListOptionRef
 * \brief for options of variable array length
 * \author J. Hicken
 */
template <class T>
class CListOptionRef : public CAnyOptionRef {
private:
	T** ref_; /*!< \brief pointer to the memory holding the list option */
	unsigned short* ref_size_; /*!< \brief number of items in list */

public:

	/*!
	 * \brief constructor for list-type options
	 * \param[in] value - variable we want to create a reference to
	 * \param[in] size - number of elements the list WILL have
	 */
	CListOptionRef(unsigned short & size, T* & value) {
		ref_ = &value;
		*ref_ = NULL;
		size = 0;
		ref_size_ = &size;
	}

	/*!
	 * \brief sets the value of the referenced convection option using vector of strings
	 * \param[in] value - a set of strings used to define the option
	 */
	void SetValue(const vector<string> & value) {
		stringstream ss;
		*ref_size_ = static_cast<unsigned short>(value.size());
		if (*ref_ != NULL) {
			cerr << "Error in CListOptionRef(SetValue): "
					<< "list option has already been allocated."
					<< endl;
			throw(-1);
		}
		(*ref_) = new T[*ref_size_];
		for (unsigned short i = 0; i < *ref_size_; i++) {
			ss << value[i] << " ";
			ss >> (*ref_)[i];
		}
	}

	/*!
	 * \brief write the value of the option to std out (mostly for debugging)
	 */
	void WriteValue() {
		for (unsigned short i = 0; i < *ref_size_; i++)
			cout << (*ref_)[i] << ", ";
		cout << endl;
	}
};

/*!
 * \class CEnumOptionRef
 * \brief a typed version of the base class for options associated with enumerations
 * \tparam T - the type of option (usually unsigned short)
 * \tparam Tenum - an enumeration assocatied with T
 * \author J. Hicken
 */
template <class T, class Tenum>
class CEnumOptionRef : public CAnyOptionRef {
private:
	T** ref_; /*!< \brief pointer to memory for the option */
	unsigned short* ref_dim_; /*!< \brief number of elemets in ref_ for array enums*/
	const map<string, Tenum> * Tmap_; /*!< \brief map between strings and enums */

public:

	/*!
	 * \brief constructor for scalar enum options
	 * \param[in] value - variable to create a reference to
	 * \param[in] Tmap - map between strings and enums
	 */
	CEnumOptionRef(T & value, const map<string, Tenum> & Tmap) {
		ref_ = new T*;
		*ref_ = &value;
		ref_dim_ = NULL;
		Tmap_ = &Tmap;
	}

	/*!
	 * \brief constructor for list enum options
	 * \param[in] value - variable to create a reference to
	 * \param[in] Tmap - map between strings and enums
	 * \param[in] size - length of the array option
	 */
	CEnumOptionRef(unsigned short & size, T* & value, const map<string, Tenum> & Tmap) {
		ref_ = &value;
		ref_dim_ = &size;
		Tmap_ = &Tmap;
	}

	/*!
	 * \brief write the value of the option to std out (mostly for debugging)
	 */
	void WriteValue() {
		if (ref_dim_ == NULL) {
			cout << *(*ref_) << endl;
		} else {
			for (unsigned int i = 0; i < *ref_dim_; i++)
				cout << (*ref_)[i] << ", ";
			cout << endl;
		}
	}

	/*!
	 * \brief sets the value of the referenced enum option using vector of strings
	 * \param[in] value - a set of strings used to define the option
	 */
	void SetValue(const vector<string> & value) {
		typename map<string,Tenum>::const_iterator it;
		if (ref_dim_ == NULL) {
			// this is a scalar enum option
			it = Tmap_->find(StringToUpperCase(value[0]));
			if (it == Tmap_->end()) {
				cerr << "Error in CEnumOptionRef::SetValue(const vector<string> &): "
						<< "cannot find value " << value[0] << " in given map."
						<< endl;
				throw(-1);
			}
			*(*ref_) = it->second;
		} else {
			// this is an array enum option
			(*ref_dim_) = static_cast<unsigned short>(value.size());
			(*ref_) = new T[*ref_dim_];
			for (unsigned short i = 0; i < *ref_dim_; i++) {
				it = Tmap_->find(StringToUpperCase(value[i]));
				if (it == Tmap_->end()) {
					cerr << "Error in CEnumOptionRef::SetValue(const vector<string> &): "
							<< "cannot find value " << value[i] << " in given map."
							<< endl;
					throw(-1);
				}
				(*ref_)[i] = it->second;
			}
		}
	}

};

/*!
 * \class CMarkerOptionRef
 * \brief a typed version of the base class for marker options
 * \author J. Hicken
 */
class CMarkerOptionRef : public CAnyOptionRef {
private:
	string** marker_ref_; /*!< \brief pointer to the memory for the marker option */
	unsigned short* num_marker_; /*!< \brief number of markers */

public:

	/*!
	 * \brief constructor for marker options
	 * \param[in] value - variable to create a reference to
	 * \param[in] size - variable refering to the number of markers
	 */
	CMarkerOptionRef(string* & value, unsigned short & size) {
		marker_ref_ = &value;
		*marker_ref_ = NULL;
		num_marker_ = &size;
	}

	/*!
	 * \brief sets the value of the referenced marker option using vector of strings
	 * \param[in] value - a set of strings used to define the option
	 */
	void SetValue(const vector<string> & value) {
		if (value.size() == 0) {
			cerr << "Error in CMarkerOptionRef::SetValue(): "
					<< "marker option in config file with no value" << endl;
			cerr << "Use NONE for no markers." << endl;
			throw(-1);
		}
		if (marker_ref_ != NULL) delete [] (*marker_ref_);
		(*marker_ref_) = new string[value.size()];
		for (unsigned int i = 0; i < value.size(); i++)
			(*marker_ref_)[i] = value[i];
		if ( (value.size() == 1) && ((*marker_ref_)[0] == "NONE") ) {
			*num_marker_ = 0;
		} else {
			*num_marker_ = value.size();
		}
	}

	/*!
	 * \brief write the value of the option to std out (mostly for debugging)
	 */
	void WriteValue() {
		for (unsigned int i = 0; i < *num_marker_; i++)
			cout << (*marker_ref_)[i] << ", ";
		cout << endl;
	}
};

/*!
 * \class CMarkerPeriodicRef
 * \brief Specialized option for periodic boundary markers
 * \author J. Hicken
 */
class CMarkerPeriodicRef : public CAnyOptionRef {
private:
	unsigned short* nMarker_PerBound_; /*!< \brief Number of periodic boundary markers. */
	string** Marker_PerBound_; /*!< \brief Periodic boundaries markers. */
	string** Marker_PerDonor_; /*!< \brief Rotationally periodic boundary donor markers. */
	double*** Periodic_RotCenter_; /*!< \brief Rotational center for each periodic boundary. */
	double*** Periodic_RotAngles_; /*!< \brief Rotation angles for each periodic boundary. */
	double*** Periodic_Translation_; /*!< \brief Translation vector for each periodic boundary. */

public:

	/*!
	 * \brief constructor for periodic marker option
	 * \param[in] nMarker_PerBound - number of periodic boundary markers
	 * \param[in] Marker_PerBound - periodic boundary markers
	 * \param[in] Marker_PerDonor - boundary marker of rotationally periodic donor
	 * \param[in] Periodic_RotCenter - rotational center for each periodic boundary
	 * \param[in] Periodic_RotAngles - rotational angles for each periodic boundary
	 * \param[in] Periodic_Translation - translation vector for each periodic boundary
	 */
	CMarkerPeriodicRef(unsigned short & nMarker_PerBound, string* & Marker_PerBound,
			string* & Marker_PerDonor, double** & Periodic_RotCenter,
			double** & Periodic_RotAngles, double** & Periodic_Translation) {
		nMarker_PerBound_ = &nMarker_PerBound;
		Marker_PerBound_ = &Marker_PerBound;
		*Marker_PerBound_ = NULL;
		Marker_PerDonor_ = &Marker_PerDonor;
		*Marker_PerDonor_ = NULL;
		Periodic_RotCenter_ = &Periodic_RotCenter;
		*Periodic_RotCenter_ = NULL;
		Periodic_RotAngles_ = &Periodic_RotAngles;
		*Periodic_RotAngles_ = NULL;
		Periodic_Translation_ = &Periodic_Translation;
		*Periodic_Translation_ = NULL;
	}

	/*!
	 * \brief sets the value of the periodic boundary parameters given the vector of strings
	 * \param[in] value - a set of strings used to define the option
	 */
	void SetValue(const vector<string> & value) {
		if ( (*Marker_PerBound_ != NULL) || (*Marker_PerDonor_ != NULL) ||
				(*Periodic_RotCenter_ != NULL) || (*Periodic_RotAngles_ != NULL) ||
				(*Periodic_Translation_ != NULL) ) {
			cerr << "Error in CMarkerPeriodicRef::SetValue(): "
					<< "one or more periodic-marker option arrays have already been allocated."
					<< endl;
			throw(-1);
		}
		if (static_cast<int>(value.size()) % 11 != 0) {
			if (value[0].compare("NONE") == 0) {
				*nMarker_PerBound_ = 0;
				return;
			}
			cerr << "Error in CMarkerPeriodicRef::SetValue(): "
					<< "incorrect number of MARKER_PERIODIC parameters in the configuration file."
					<< endl;
			throw(-1);
		}
		*nMarker_PerBound_ = static_cast<unsigned short>(value.size())/11;
		(*Marker_PerBound_)      = new string[*nMarker_PerBound_];
		(*Marker_PerDonor_)      = new string[*nMarker_PerBound_];
		(*Periodic_RotCenter_)   = new double*[*nMarker_PerBound_];
		(*Periodic_RotAngles_)   = new double*[*nMarker_PerBound_];
		(*Periodic_Translation_) = new double*[*nMarker_PerBound_];
		for (unsigned short iMarker_PerBound = 0; iMarker_PerBound < *nMarker_PerBound_; iMarker_PerBound++) {
			(*Periodic_RotCenter_)[iMarker_PerBound]   = new double[3];
			(*Periodic_RotAngles_)[iMarker_PerBound]   = new double[3];
			(*Periodic_Translation_)[iMarker_PerBound] = new double[3];
		}

		stringstream ss;
		unsigned short i = 0;
		double deg2rad = PI_NUMBER/180.0;
		for (unsigned short iMarker_PerBound = 0; iMarker_PerBound < *nMarker_PerBound_; iMarker_PerBound++) {
			ss << value[i++] << " ";
			ss >> (*Marker_PerBound_)[iMarker_PerBound];
			ss << value[i++] << " ";
			ss >> (*Marker_PerDonor_)[iMarker_PerBound];
			ss << value[i++] << " ";
			ss >> (*Periodic_RotCenter_)[iMarker_PerBound][0];
			ss << value[i++] << " ";
			ss >> (*Periodic_RotCenter_)[iMarker_PerBound][1];
			ss << value[i++] << " ";
			ss >> (*Periodic_RotCenter_)[iMarker_PerBound][2];
			ss << value[i++] << " ";
			ss >> (*Periodic_RotAngles_)[iMarker_PerBound][0];
			ss << value[i++] << " ";
			ss >> (*Periodic_RotAngles_)[iMarker_PerBound][1];
			ss << value[i++] << " ";
			ss >> (*Periodic_RotAngles_)[iMarker_PerBound][2];
			ss << value[i++] << " ";
			ss >> (*Periodic_Translation_)[iMarker_PerBound][0];
			ss << value[i++] << " ";
			ss >> (*Periodic_Translation_)[iMarker_PerBound][1];
			ss << value[i++] << " ";
			ss >> (*Periodic_Translation_)[iMarker_PerBound][2];

			/*--- Convert the angles from degrees to radians ---*/
			(*Periodic_RotAngles_)[iMarker_PerBound][0] *= deg2rad;
			(*Periodic_RotAngles_)[iMarker_PerBound][1] *= deg2rad;
			(*Periodic_RotAngles_)[iMarker_PerBound][2] *= deg2rad;
		}
	}

	/*!
	 * \brief write the value of the option to std out (mostly for debugging)
	 */
	void WriteValue() {
		cout << "CMarkerPeriodicRef::WriteValue(): not implemented yet" << endl;
	}

};

/*!
 * \class CMarkerInletRef
 * \brief Specialized option for inlet boundary markers
 * \author J. Hicken
 */
class CMarkerInletRef : public CAnyOptionRef {
private:
	unsigned short* nMarker_Inlet_; /*!< \brief the number of inlet marker boundaries */
	string** Marker_Inlet_;         /*!< \brief string names of inlet boundaries */
	double** Ttotal_;               /*!< \brief specified total temperatures for inlet boundaries */
	double** Ptotal_;               /*!< \brief specified total pressures for inlet boundaries */
	double*** FlowDir_;             /*!< \brief specified flow direction vector (unit vector) for inlet boundaries */

public:

	/*!
	 * \brief constructor for inlet marker option
	 * \param[in] nMarker_Inlet - the number of inlet marker boundaries
	 * \param[in] Marker_Inlet - string names of inlet boundaries
	 * \param[in] Ttotal - specified total temperatures for inlet boundaries
	 * \param[in] Ptotal - specified total pressures for inlet boundaries
	 * \param[in] FlowDir - specified flow direction vector (unit vector) for inlet boundaries
	 */
	CMarkerInletRef(unsigned short & nMarker_Inlet, string* & Marker_Inlet, double* & Ttotal,
			double* & Ptotal, double** & FlowDir) {
		nMarker_Inlet_ = &nMarker_Inlet;
		Marker_Inlet_ = &Marker_Inlet;
		*Marker_Inlet_ = NULL;
		Ttotal_ = &Ttotal;
		*Ttotal_ = NULL;
		Ptotal_ = &Ptotal;
		*Ptotal_ = NULL;
		FlowDir_ = &FlowDir;
		*FlowDir_ = NULL;
	}

	/*!
	 * \brief sets the value of the inlet parameters given the vector of strings
	 * \param[in] value - a set of strings used to define the option
	 */
	void SetValue(const vector<string> & value) {
		if ( (*Marker_Inlet_ != NULL) || (*Ttotal_ != NULL) || (*Ptotal_ != NULL) ||
				(*FlowDir_ != NULL) ) {
			cerr << "Error in CMarkerInletRef::SetValue(): "
					<< "one or more inlet-marker option arrays have already been allocated."
					<< endl;
			throw(-1);
		}
		if (static_cast<int>(value.size()) % 6 != 0) {
			if (value[0].compare("NONE") == 0) {
				*nMarker_Inlet_ = 0;
				return;
			}
			cerr << "Error in CMarkerInletRef::SetValue(): "
					<< "incorrect number of MARKER_INLET parameters in the configuration file."
					<< endl;
			throw(-1);
		}
		*nMarker_Inlet_ = static_cast<unsigned short>(value.size())/6;
		(*Marker_Inlet_) = new string[*nMarker_Inlet_];
		(*Ttotal_)       = new double[*nMarker_Inlet_];
		(*Ptotal_)       = new double[*nMarker_Inlet_];
		(*FlowDir_)      = new double*[*nMarker_Inlet_];
		for (unsigned short iMarker = 0; iMarker < *nMarker_Inlet_; iMarker++)
			(*FlowDir_)[iMarker]   = new double[3];

		stringstream ss;
		unsigned short i = 0;
		for (unsigned short iMarker = 0; iMarker < *nMarker_Inlet_; iMarker++) {
			ss << value[i++] << " ";
			ss >> (*Marker_Inlet_)[iMarker];
			ss << value[i++] << " ";
			ss >> (*Ttotal_)[iMarker];
			ss << value[i++] << " ";
			ss >> (*Ptotal_)[iMarker];
			ss << value[i++] << " ";
			ss >> (*FlowDir_)[iMarker][0];
			ss << value[i++] << " ";
			ss >> (*FlowDir_)[iMarker][1];
			ss << value[i++] << " ";
			ss >> (*FlowDir_)[iMarker][2];
		}
	}

	/*!
	 * \brief write the value of the option to std out (mostly for debugging)
	 */
	void WriteValue() {
		cout << "CMarkerInletRef::WriteValue(): not implemented yet" << endl;
	}

};

/*!
 * \class CMarkerOutletRef
 * \brief Specialized option for outlet boundary markers
 * \author J. Hicken
 */
class CMarkerOutletRef : public CAnyOptionRef {
private:
	unsigned short* nMarker_Outlet_; /*!< \brief the number of outlet marker boundaries */
	string** Marker_Outlet_;         /*!< \brief string names of outlet boundaries */
	double** Pressure_;              /*!< \brief specified back pressures (static) for outlet boundaries */

public:

	/*!
	 * \brief constructor for outlet marker option
	 * \param[in] nMarker_Outlet - the number of outlet marker boundaries
	 * \param[in] Marker_Outlet - string names of outlet boundaries
	 * \param[in] Pressure - specified back pressures (static) for outlet boundaries
	 */
	CMarkerOutletRef(unsigned short & nMarker_Outlet, string* & Marker_Outlet, double* & Pressure) {
		nMarker_Outlet_ = &nMarker_Outlet;
		Marker_Outlet_ = &Marker_Outlet;
		*Marker_Outlet_ = NULL;
		Pressure_ = &Pressure;
		*Pressure_ = NULL;
	}

	/*!
	 * \brief sets the value of the outlet parameters given the vector of strings
	 * \param[in] value - a set of strings used to define the option
	 */
	void SetValue(const vector<string> & value) {
		if ( (*Marker_Outlet_ != NULL) || (*Pressure_ != NULL) ) {
			cerr << "Error in CMarkerOutletRef::SetValue(): "
					<< "one or more outlet-marker option arrays have already been allocated."
					<< endl;
			throw(-1);
		}
		if (static_cast<int>(value.size()) % 2 != 0) {
			if (value[0].compare("NONE") == 0) {
				*nMarker_Outlet_ = 0;
				return;
			}
			cerr << "Error in CMarkerOutletRef::SetValue(): "
					<< "incorrect number of MARKER_OUTLET parameters in the configuration file."
					<< endl;
			throw(-1);
		}
		*nMarker_Outlet_  = static_cast<unsigned short>(value.size())/2;
		(*Marker_Outlet_) = new string[*nMarker_Outlet_];
		(*Pressure_)      = new double[*nMarker_Outlet_];
		stringstream ss;
		unsigned short i = 0;
		for (unsigned short iMarker = 0; iMarker < *nMarker_Outlet_; iMarker++) {
			ss << value[i++] << " ";
			ss >> (*Marker_Outlet_)[iMarker];
			ss << value[i++] << " ";
			ss >> (*Pressure_)[iMarker];
		}
	}

	/*!
	 * \brief write the value of the option to std out (mostly for debugging)
	 */
	void WriteValue() {
		cout << "Outlet markers (" << (*nMarker_Outlet_) << ")" << endl;
		for (unsigned short iMarker = 0; iMarker < *nMarker_Outlet_; iMarker++)
			cout << "name = " << (*Marker_Outlet_)[iMarker]
			                                       << ": back press. = " << (*Pressure_)[iMarker] << endl;
	}

};

/*!
 * \class CConvOptionRef
 * \brief a typed version of the base class for convection discretization options
 * \author J. Hicken
 */
class CConvOptionRef : public CAnyOptionRef {
private:
	unsigned short* space_ref_; /*!< \brief pointer to the space discretization type */
	unsigned short* centred_ref_; /*!< \brief pointer to the centred discretization type */
	unsigned short* upwind_ref_; /*!< \brief pointer to the upwind discretization type */

public:

	/*!
	 * \brief constructor for convection options
	 * \param[in] space - space discretization variable to create a reference to
	 * \param[in] centred - centred discretization variable to create a reference to
	 * \param[in] upwind - upwind discretization variable to create a reference to
	 */
	CConvOptionRef(unsigned short & space, unsigned short & centred,
			unsigned short & upwind) {
		space_ref_ = &space;
		centred_ref_ = &centred;
		upwind_ref_ = &upwind;
	}

	/*!
	 * \brief sets the value of the referenced convection option using vector of strings
	 * \param[in] value - a set of strings used to define the option
	 */
	void SetValue(const vector<string> & value) {
		if (Centred_Map.count(value[0])) {
			*space_ref_ = Space_Map.find("SPACE_CENTRED")->second;
			*centred_ref_ = Centred_Map.find(value[0])->second;
			*upwind_ref_ = NO_UPWIND;
		} else if (Upwind_Map.count(value[0])) {
			*space_ref_ = Space_Map.find("SPACE_UPWIND")->second;
			*upwind_ref_ = Upwind_Map.find(value[0])->second;
			*centred_ref_ = NO_CENTRED;
		} else {
			cerr << "Error in CConvOptionRef::SetValue(): "
					<< value[0] << " is an invalid space discretization" << endl;
			throw(-1);
		}
	}

	/*!
	 * \brief write the value of the option to std out (mostly for debugging)
	 */
	void WriteValue() {
		cout << "space discretization: " << *space_ref_ << endl;
		cout << "centred discretization: " << *centred_ref_ << endl;
		cout << "upwind discretization: " << *upwind_ref_ << endl;
	}
};

/*!
 * \class CMathProblemRef
 * \brief this very specialized option for MATH_PROBLEM set several variables
 * \author J. Hicken
 */
class CMathProblemRef : public CAnyOptionRef {
private:
	bool* ContAdj_; /*!< \brief pointer to the ContAdj variable */
	bool* OneShot_; /*!< \brief pointer to the OneShot variable */
	bool* Linearized_; /*!< \brief pointer to the Linerized variable */
	bool* Restart_Flow_; /*!< \brief pointer to the Restart_Flow variable */

public:

	/*!
	 * \brief constructor for math problem option
	 * \param[in] ContAdj - is the continuous adjoint being solved
	 * \param[in] OneShot - is a one shot problem being solved
	 * \param[in] Linearized - is a linearized problem being solved
	 * \param[in] Restart_Flow - restart flow solution for adjoint and linearized problems
	 */
	CMathProblemRef(bool & ContAdj, bool & OneShot, bool & Linearized,
			bool & Restart_Flow) {
		ContAdj_ = &ContAdj;
		OneShot_ = &OneShot;
		Linearized_ = &Linearized;
		Restart_Flow_ = &Restart_Flow;
	}

	/*!
	 * \brief sets the value of the math problem given the vector of strings
	 * \param[in] value - a set of strings used to define the option
	 */
	void SetValue(const vector<string> & value) {
		if (value[0] == "DIRECT") {
			*ContAdj_ = false; *OneShot_ = false; *Linearized_ = false; *Restart_Flow_ = false;
		} else if (value[0] == "ADJOINT") {
			*ContAdj_ = true; *Restart_Flow_ = true;
		}
		if (value[0] == "LINEARIZED") {
			*Linearized_ = true; *Restart_Flow_ = true;
		}
	}

	/*!
	 * \brief write the value of the option to std out (mostly for debugging)
	 */
	void WriteValue() {
		cout << "ContAdj_ = " << *ContAdj_;
		cout << ": OneShot_ = " << *OneShot_;
		cout << ": Linearized_ = " << *Linearized_;
		cout << ": Restart_Flow_ = " << *Restart_Flow_ << endl;
	}
};

/*!
 * \class CDVParamOptionRef
 * \brief Specialized option for grid deformation parameters
 * \author J. Hicken
 */
class CDVParamOptionRef : public CAnyOptionRef {
private:
	unsigned short* nDV_; /*!< \brief Number of variables. */
	double*** ParamDV_; /*!< \brief Parameters of the variables. */
	unsigned short** Design_Variable_; /*!< \brief design variable types; should already be allocated */

public:

	/*!
	 * \brief constructor for DVParam option
	 * \param[in] nDV - number of design variables
	 * \param[in] ParamDV - the paramters of each design variable
	 * \param[in] Design_Variable - DV types
	 */
	CDVParamOptionRef(unsigned short & nDV, double** & ParamDV,
			unsigned short* & Design_Variable) {
		nDV_ = &nDV;
		ParamDV_ = &ParamDV;
		Design_Variable_ = &Design_Variable;
	}

	/*!
	 * \brief sets the value of the design variables given the vector of strings
	 * \param[in] value - a set of strings used to define the option
	 */
	void SetValue(const vector<string> & value) {

		// use the ";" token to determine the number of design variables
		*nDV_ = 0;
		unsigned int num_semi = 0;
		for (unsigned int i = 0; i < static_cast<unsigned int>(value.size()); i++) {
			if (value[i].compare(";") == 0) {
				(*nDV_)++;
				num_semi++;
			}
		}
		// if ";" are at both ends, we over-counted
		if ( (value[0].compare(";") == 0) && (value[value.size()-1].compare(";") == 0) )
			(*nDV_)--;
		// if no ";" at either end, we under-counted
		if ( (value[0].compare(";") != 0) && (value[value.size()-1].compare(";") != 0) )
			(*nDV_)++;

		if ( (*nDV_ > 0) && (*Design_Variable_ == NULL) ) {
			cerr << "Error in CDVParamOptionRef::SetValue(): "
					<< "Design_Variable array has not been allocated." << endl;
			cerr << "Check that DV_KIND appears before DV_PARAM in configuration file." << endl;
			throw(-1);
		}

#if 0
		cout << "Found " << (*nDV_) << " DV parameters" << endl;
		cout << "DV param value = ";
		for (unsigned int i = 0; i < value.size(); i++)
			cout << value[i] << ", ";
		cout << endl;
#endif
    
		
		
		(*ParamDV_) = new double*[*nDV_];
		for (unsigned short iDV = 0; iDV < *nDV_; iDV++)
			(*ParamDV_)[iDV] = new double[MAX_PARAMETERS];
		
		unsigned short nParamDV = 0;
		stringstream ss;
		unsigned int i = 0;
		for (unsigned short iDV = 0; iDV < *nDV_; iDV++) {
			switch ((*Design_Variable_)[iDV]) {
				case NO_DEFORMATION: nParamDV = 0; break;
				case HICKS_HENNE: nParamDV = 2; break;
				case DISPLACEMENT: nParamDV = 3; break;
				case ROTATION: nParamDV = 6; break;
				case NACA_4DIGITS: nParamDV = 3; break;
				case PARABOLIC: nParamDV = 2; break;
				case OBSTACLE: nParamDV = 2; break;
				case STRETCH: nParamDV = 2; break;
				case FFD_CONTROL_POINT: nParamDV = 7; break;
				case FFD_DIHEDRAL_ANGLE: nParamDV = 7; break;
				case FFD_TWIST_ANGLE: nParamDV = 7; break;
				case FFD_ROTATION: nParamDV = 7; break;
				case FFD_CAMBER: nParamDV = 3; break;
				case FFD_THICKNESS: nParamDV = 3; break;
				case FFD_VOLUME: nParamDV = 3; break;
				default : {
					cerr << "Error in CDVParamOptionRef::SetValue(): "
					<< "undefined design variable type found in configuration file." << endl;
				}
			}
			for (unsigned short iParamDV = 0; iParamDV < nParamDV; iParamDV++) {
				ss << value[i] << " ";
				ss >> (*ParamDV_)[iDV][iParamDV];
				i++;
			}
			if (iDV < (*nDV_)-1) {
				if (value[i].compare(";") != 0) {
					cerr << "Error in CDVParamOptionRef::SetValue(): "
					<< "a design variable in the configuration file "
					<< "has the wrong number of parameters" << endl;
					throw(-1);
				}
				i++;
			}
		}
	}
	
	
	/*!
	 * \brief write the value of the option to std out (mostly for debugging)
	 */
	void WriteValue() {
		//cout << "CDVParamOptionRef::WriteValue(): not implemented yet" << endl;
		for (unsigned short iDV = 0; iDV < *nDV_; iDV++) {
			unsigned short nParamDV = 0;
			switch ((*Design_Variable_)[iDV]) {
				case NO_DEFORMATION: nParamDV = 0; break;
				case HICKS_HENNE: nParamDV = 2; break;
				case DISPLACEMENT: nParamDV = 3; break;
				case ROTATION: nParamDV = 6; break;
				case NACA_4DIGITS: nParamDV = 3; break;
				case PARABOLIC: nParamDV = 2; break;
				case FFD_CONTROL_POINT: nParamDV = 7; break;
				case FFD_DIHEDRAL_ANGLE: nParamDV = 7; break;
				case FFD_TWIST_ANGLE: nParamDV = 7; break;
				case FFD_ROTATION: nParamDV = 7; break;
				case FFD_CAMBER: nParamDV = 3; break;
				case FFD_THICKNESS: nParamDV = 3; break;
				case FFD_VOLUME: nParamDV = 3; break;
				default : {
					cerr << "Error in CDVParamOptionRef::SetValue(): "
					<< "undefined design variable type found in configuration file." << endl;
				}
			}
			cout << "DV param type: " << (*Design_Variable_)[iDV] << ": values = ";
			for (unsigned short iParamDV = 0; iParamDV < nParamDV; iParamDV++)
				cout << (*ParamDV_)[iDV][iParamDV] << ", ";
			cout << endl;
		}
	}
};
