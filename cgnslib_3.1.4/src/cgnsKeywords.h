/*
   @#@#@ $Id: cgnsKeywords.h,v 1.3 2008/07/11 00:25:35 brucewedan Exp $
   @#@#@ CHECKSUM: 296864779 40788
   @#@#@ (cat cgnsKeywords.h | grep -v '@#@#@' | cksum)
*/

/*-------------------------------------------------------------------------
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from
the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not
   be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution.
-------------------------------------------------------------------------*/

#ifndef CGNSLIB_KEYWORDS_H
#define CGNSLIB_KEYWORDS_H

/*
   THIS HEADER CAN BE READ BY ANY C LANGAGE PRE-PROCESSOR

   ENUMERATES DEFINITIONS ARE DECLARED, IN THE CASE YOU DO NOT WANT
   THESE DEFINITIONS (e.g. you compile Fortran) THEN YOU HAVE TO UNSET
   THE __CGNS_ENUMS__

#undef __CGNS_ENUMS__

   BY DEFAULT YOU WANT TO DECLARE THESE ENUMERATES

#define __CGNS_ENUMS__

*/

#define __CGNS_ENUMS__

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      VERSION NUMBER                                                   *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define CGNS_VERSION 3130
#define CGNS_DOTVERS 3.13
#define CGNS_COMPATVERSION 2540
#define CGNS_COMPATDOTVERS 2.54

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      modes for cgns file                                              *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      modes for cgns file                                              *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define CG_MODE_READ	0
#define CG_MODE_WRITE	1
#define CG_MODE_MODIFY	2
#define CG_MODE_CLOSED	3

/* function return codes */

#define CG_OK		  0
#define CG_ERROR	  1
#define CG_NODE_NOT_FOUND 2
#define CG_INCORRECT_PATH 3
#define CG_NO_INDEX_DIM   4

/* Null and UserDefined enums */

#define CG_Null        0
#define CG_UserDefined 1

/* max goto depth */

#define CG_MAX_GOTO_DEPTH 20

/* configuration options */

#define CG_CONFIG_ERROR     1
#define CG_CONFIG_COMPRESS  2
#define CG_CONFIG_SET_PATH  3
#define CG_CONFIG_ADD_PATH  4
#define CG_CONFIG_FILE_TYPE 5

/* The Null and the UserDefined string constants are always the same,
   whatever enumerate you consider. */
#define CG_Null_s                      "Null"
#define CG_UserDefined_s               "Userdefined"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Dimensional Units                                                *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Mass ---*/
#ifdef __CGNS_ENUMS__
typedef enum {
	MassUnitsNull, MassUnitsUserDefined,
	Kilogram, Gram, Slug, PoundMass
} MassUnits_t;
#endif
#define NofValidMassUnits              6

#define Kilogram_s                     "Kilogram"
#define Gram_s                         "Gram"
#define Slug_s                         "Slug"
#define PoundMass_s                    "PoundMass"
#define Meter_s                        "Meter"

/* Length ---*/
#ifdef __CGNS_ENUMS__
typedef enum {
	LengthUnitsNull, LengthUnitsUserDefined,
	Meter, Centimeter, Millimeter, Foot, Inch
} LengthUnits_t;
#endif
#define NofValidLengthUnits            7

#define Centimeter_s                   "Centimeter"
#define Millimeter_s                   "Millimeter"
#define Foot_s                         "Foot"
#define Inch_s                         "Inch"

/* Mass ---*/
#ifdef __CGNS_ENUMS__
typedef enum {
	TimeUnitsNull, TimeUnitsUserDefined, Second
} TimeUnits_t;
#endif
#define NofValidTimeUnits              3

#define Second_s                       "Second"

/* Temperature ---*/
#ifdef __CGNS_ENUMS__
typedef enum {
	TemperatureUnitsNull, TemperatureUnitsUserDefined,
	Kelvin, Celsius, Rankine, Fahrenheit
} TemperatureUnits_t;
#endif
#define NofValidTemperatureUnits       6

#define Kelvin_s                       "Kelvin"
#define Celsius_s                      "Celsius"
#define Rankine_s                      "Rankine"
#define Fahrenheit_s                   "Fahrenheit"

/* Angle ---*/
#ifdef __CGNS_ENUMS__
typedef enum {
	AngleUnitsNull, AngleUnitsUserDefined, Degree, Radian
} AngleUnits_t;
#endif
#define NofValidAngleUnits             4

#define Degree_s                       "Degree"
#define Radian_s                       "Radian"

/* ElectricCurrent ---*/
#ifdef __CGNS_ENUMS__
typedef enum {
	ElectricCurrentUnitsNull, ElectricCurrentUnitsUserDefined,
	Ampere, Abampere, Statampere, Edison, auCurrent
} ElectricCurrentUnits_t;
#endif
#define NofValidElectricCurrentUnits   7

#define Ampere_s                       "Ampere"
#define Abampere_s                     "Abampere"
#define Statampere_s                   "Statampere"
#define Edison_s                       "Edison"
#define auCurrent_s                    "auCurrent"

/* SubstanceAmount ---*/
#ifdef __CGNS_ENUMS__
typedef enum {
	SubstanceAmountUnitsNull, SubstanceAmountUnitsUserDefined,
	Mole, Entities, StandardCubicFoot, StandardCubicMeter
} SubstanceAmountUnits_t;
#endif
#define NofValidSubstanceAmountUnits   6

#define Mole_s                         "Mole"
#define Entities_s                     "Entities"
#define StandardCubicFoot_s            "StandardCubicFoot"
#define StandardCubicMeter_s           "StandardCubicMeter"

/* LuminousIntensity ---*/
#ifdef __CGNS_ENUMS__
typedef enum {
	LuminousIntensityUnitsNull, LuminousIntensityUnitsUserDefined,
	Candela, Candle, Carcel, Hefner, Violle
} LuminousIntensityUnits_t;
#endif
#define NofValidLuminousIntensityUnits 7

#define Candela_s                      "Candela"
#define Candle_s                       "Candle"
#define Carcel_s                       "Carcel"
#define Hefner_s                       "Hefner"
#define Violle_s                       "Violle"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Data Class                                                       *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	DataClassNull, DataClassUserDefined,
	Dimensional, NormalizedByDimensional,
	NormalizedByUnknownDimensional,
	NondimensionalParameter, DimensionlessConstant
} DataClass_t;
#endif
#define NofValidDataClass 7

#define Dimensional_s                    "Dimensional"
#define NormalizedByDimensional_s        "NormalizedByDimensional"
#define NormalizedByUnknownDimensional_s "NormalizedByUnknownDimensional"
#define NondimensionalParameter_s        "NondimensionalParameter"
#define DimensionlessConstant_s          "DimensionlessConstant"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *	Grid Location
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	GridLocationNull, GridLocationUserDefined,
        Vertex, CellCenter, FaceCenter,
        IFaceCenter, JFaceCenter, KFaceCenter, EdgeCenter
} GridLocation_t;
#endif
#define NofValidGridLocation 9

#define Vertex_s                       "Vertex"
#define CellCenter_s                   "CellCenter"
#define FaceCenter_s                   "FaceCenter"
#define IFaceCenter_s                  "IFaceCenter"
#define JFaceCenter_s                  "JFaceCenter"
#define KFaceCenter_s                  "KFaceCenter"
#define EdgeCenter_s                   "EdgeCenter"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      BCData Types                                                     *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	BCDataTypeNull, BCDataTypeUserDefined,
	Dirichlet, Neumann
} BCDataType_t;
#endif
#define NofValidBCDataTypes 4

#define Dirichlet_s                    "Dirichlet"
#define Neumann_s                      "Neumann"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *	Grid Connectivity Types 					 *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	GridConnectivityTypeNull, GridConnectivityTypeUserDefined,
	Overset, Abutting, Abutting1to1
} GridConnectivityType_t;
#endif
#define NofValidGridConnectivityTypes 5

#define Overset_s                      "Overset"
#define Abutting_s                     "Abutting"
#define Abutting1to1_s                 "Abutting1to1"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *	Point Set Types							 *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	PointSetTypeNull, PointSetTypeUserDefined,
        PointList,  PointListDonor,
        PointRange, PointRangeDonor,
	ElementRange, ElementList, CellListDonor
} PointSetType_t;

#endif
#define NofValidPointSetTypes 9

#define PointList_s                    "PointList"
#define PointListDonor_s               "PointListDonor"
#define PointRange_s                   "PointRange"
#define PointRangeDonor_s              "PointRangeDonor"
#define ElementRange_s                 "ElementRange"
#define ElementList_s                  "ElementList"
#define CellListDonor_s                "CellListDonor"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Governing Equations and Physical Models Types                    *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	GoverningEquationsNull, GoverningEquationsUserDefined,
	FullPotential, Euler, NSLaminar, NSTurbulent,
	NSLaminarIncompressible, NSTurbulentIncompressible
} GoverningEquationsType_t;
#endif
#define NofValidGoverningEquationsTypes 8

#define FullPotential_s                "FullPotential"
#define Euler_s                        "Euler"
#define NSLaminar_s                    "NSLaminar"
#define NSTurbulent_s                  "NSTurbulent"
#define NSLaminarIncompressible_s      "NSLaminarIncompressible"
#define NSTurbulentIncompressible_s    "NSTurbulentIncompressible"

/* Any model type will accept both ModelTypeNull and ModelTypeUserDefined.
** The following models will accept these values as vaild...
**
** GasModel_t: Ideal, VanderWaals, CaloricallyPerfect, ThermallyPerfect,
**    ConstantDensity, RedlichKwong
**
** ViscosityModel_t: Constant, PowerLaw, SutherlandLaw
**
** ThermalConductivityModel_t: PowerLaw, SutherlandLaw, ConstantPrandtl
**
** TurbulenceModel_t: Algebraic_BaldwinLomax, Algebraic_CebeciSmith,
**    HalfEquation_JohnsonKing, OneEquation_BaldwinBarth,
**    OneEquation_SpalartAllmaras, TwoEquation_JonesLaunder,
**    TwoEquation_MenterSST,TwoEquation_Wilcox
**
** TurbulenceClosure_t: EddyViscosity, ReynoldsStress, ReynoldsStressAlgebraic
**
** ThermalRelaxationModel_t: Frozen, ThermalEquilib, ThermalNonequilib
**
** ChemicalKineticsModel_t: Frozen, ChemicalEquilibCurveFit,
**    ChemicalEquilibMinimization, ChemicalNonequilib
**
** EMElectricFieldModel_t: Voltage, Interpolated, Constant, Frozen
**
** EMMagneticFieldModel_t: Interpolated, Constant, Frozen
**
** EMConductivityModel_t: Constant, Frozen, Equilibrium_LinRessler,
**				Chemistry_LinRessler
*/

#ifdef __CGNS_ENUMS__
typedef enum {
	ModelTypeNull, ModelTypeUserDefined,
	Ideal, VanderWaals,
	Constant,
	PowerLaw, SutherlandLaw,
	ConstantPrandtl,
	EddyViscosity, ReynoldsStress, ReynoldsStressAlgebraic,
	Algebraic_BaldwinLomax, Algebraic_CebeciSmith,
	HalfEquation_JohnsonKing, OneEquation_BaldwinBarth,
	OneEquation_SpalartAllmaras, TwoEquation_JonesLaunder,
	TwoEquation_MenterSST, TwoEquation_Wilcox,
	CaloricallyPerfect, ThermallyPerfect,
	ConstantDensity, RedlichKwong,
	Frozen, ThermalEquilib, ThermalNonequilib,
	ChemicalEquilibCurveFit, ChemicalEquilibMinimization,
	ChemicalNonequilib,
	EMElectricField, EMMagneticField, EMConductivity,
	Voltage, Interpolated, Equilibrium_LinRessler, Chemistry_LinRessler
} ModelType_t;
#endif
#define NofValidModelTypes 36

#define Ideal_s                        "Ideal"
#define VanderWaals_s                  "VanderWaals"
#define Constant_s                     "Constant"
#define PowerLaw_s                     "PowerLaw"    
#define SutherlandLaw_s                "SutherlandLaw"
#define ConstantPrandtl_s              "ConstantPrandtl"
#define EddyViscosity_s                "EddyViscosity"
#define ReynoldsStress_s               "ReynoldsStress"
#define Algebraic_s                    "Algebraic"
#define BaldwinLomax_s                 "BaldwinLomax"
#define ReynoldsStressAlgebraic_s      "ReynoldsStressAlgebraic"
#define Algebraic_CebeciSmith_s	       "Algebraic_CebeciSmith"
#define HalfEquation_JohnsonKing_s     "HalfEquation_JohnsonKing"
#define OneEquation_BaldwinBarth_s     "OneEquation_BaldwinBarth"
#define OneEquation_SpalartAllmaras_s  "OneEquation_SpalartAllmaras"
#define TwoEquation_JonesLaunder_s     "TwoEquation_JonesLaunder"
#define TwoEquation_MenterSST_s        "TwoEquation_MenterSST"
#define TwoEquation_Wilcox_s           "TwoEquation_Wilcox"
#define	CaloricallyPerfect_s           "CaloricallyPerfect"
#define ThermallyPerfect_s             "ThermallyPerfect"
#define ConstantDensity_s              "ConstantDensity"
#define RedlichKwong_s                 "RedlichKwong"
#define Frozen_s                       "Frozen"
#define ThermalEquilib_s               "ThermalEquilib"
#define ThermalNonequilib_s            "ThermalNonequilib"
#define ChemicalEquilibCurveFit_s      "ChemicalEquilibCurveFit"
#define ChemicalEquilibMinimization_s  "ChemicalEquilibMinimization"
#define ChemicalNonequilib_s           "ChemicalNonequilib"
#define EMElectricField_s              "EMElectricField"
#define EMMagneticField_s              "EMMagneticField"
#define EMConductivity_s               "EMConductivity"
#define Voltage_s                      "Voltage"
#define Interpolated_s                 "Interpolated"
#define Equilibrium_LinRessler_s       "Equilibrium_LinRessler"
#define Chemistry_LinRessler_s         "Chemistry_LinRessler"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 * 	Boundary Condition Types					 *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	BCTypeNull, BCTypeUserDefined,
	BCAxisymmetricWedge, BCDegenerateLine, BCDegeneratePoint,
	BCDirichlet, BCExtrapolate, BCFarfield, BCGeneral, BCInflow,
	BCInflowSubsonic,  BCInflowSupersonic, BCNeumann, BCOutflow,
	BCOutflowSubsonic, BCOutflowSupersonic, BCSymmetryPlane,
	BCSymmetryPolar, BCTunnelInflow, BCTunnelOutflow, BCWall,
	BCWallInviscid, BCWallViscous, BCWallViscousHeatFlux,
	BCWallViscousIsothermal, FamilySpecified
} BCType_t;
#endif
#define NofValidBCTypes 26

#define BCAxisymmetricWedge_s          "BCAxisymmetricWedge"
#define BCDegenerateLine_s             "BCDegenerateLine"
#define BCDegeneratePoint_s            "BCDegeneratePoint"
#define BCDirichlet_s                  "BCDirichlet"
#define BCExtrapolate_s                "BCExtrapolate"
#define BCFarfield_s                   "BCFarfield"
#define BCGeneral_s                    "BCGeneral"
#define BCInflow_s                     "BCInflow"
#define BCInflowSubsonic_s             "BCInflowSubsonic"
#define BCInflowSupersonic_s           "BCInflowSupersonic"
#define BCNeumann_s                    "BCNeumann"
#define BCOutflow_s                    "BCOutflow"
#define BCOutflowSubsonic_s            "BCOutflowSubsonic"
#define BCOutflowSupersonic_s          "BCOutflowSupersonic"
#define BCSymmetryPlane_s              "BCSymmetryPlane"
#define BCSymmetryPolar_s              "BCSymmetryPolar"
#define BCTunnelInflow_s               "BCTunnelInflow"
#define BCTunnelOutflow_s              "BCTunnelOutflow"
#define BCWall_s                       "BCWall"
#define BCWallInviscid_s               "BCWallInviscid"
#define BCWallViscous_s                "BCWallViscous"
#define BCWallViscousHeatFlux_s        "BCWallViscousHeatFlux"
#define BCWallViscousIsothermal_s      "BCWallViscousIsothermal"
#define FamilySpecified_s              "FamilySpecified_s"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Data types:  Can not add data types and stay forward compatible  *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	DataTypeNull, DataTypeUserDefined, Integer, RealSingle,
	RealDouble, Character, LongInteger
} DataType_t;
#endif
#define NofValidDataTypes 7

#define Integer_s                      "Integer"
#define RealSingle_s                   "RealSingle"
#define RealDouble_s                   "RealDouble"
#define Character_s                    "Character"
#define LongInteger_s                  "LongInteger"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Element types                                                    *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	ElementTypeNull, ElementTypeUserDefined,	/* 0, 1,	*/
	NODE, BAR_2, BAR_3, 				/* 2, 3, 4, 	*/
	TRI_3, TRI_6,					/* 5, 6,	*/
	QUAD_4, QUAD_8, QUAD_9,				/* 7, 8, 9,	*/
	TETRA_4, TETRA_10, 				/* 10, 11,	*/
	PYRA_5, PYRA_14, 				/* 12, 13,	*/
	PENTA_6, PENTA_15, PENTA_18,			/* 14, 15, 16,	*/
	HEXA_8, HEXA_20, HEXA_27, 			/* 17, 18, 19,	*/
	MIXED, PYRA_13, NGON_n, NFACE_n			/* 20, 21, 22, 23*/
} ElementType_t;
#endif
#define NofValidElementTypes 24

#define NODE_s                         "NODE"
#define BAR_2_s                        "BAR_2"
#define BAR_3_s                        "BAR_3"
#define TRI_3_s                        "TRI_3"
#define TRI_6_s                        "TRI_6"
#define QUAD_4_s                       "QUAD_4"
#define QUAD_8_s                       "QUAD_8"
#define QUAD_9_s                       "QUAD_9"
#define TETRA_4_s                      "TETRA_4"
#define TETRA_10_s                     "TETRA_10"
#define PYRA_5_s                       "PYRA_5"
#define PYRA_14_s                      "PYRA_14"
#define PENTA_6_s                      "PENTA_6"
#define PENTA_15_s                     "PENTA_15"
#define PENTA_18_s                     "PENTA_18"
#define HEXA_8_s                       "HEXA_8"
#define HEXA_20_s                      "HEXA_20"
#define HEXA_27_s                      "HEXA_27"
#define MIXED_s                        "MIXED"
#define PYRA_13_s                      "PYRA_13"
#define NGON_n_s                       "NGON_n"
#define NFACE_n_s                      "NFACE_n"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Zone types                                                       *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	ZoneTypeNull, ZoneTypeUserDefined,
	Structured, Unstructured
} ZoneType_t;
#endif
#define NofValidZoneTypes 4

#define Structured_s                   "Structured"
#define Unstructured_s                 "Unstructured"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Rigid Grid Motion types						 *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	RigidGridMotionTypeNull, RigidGridMotionTypeUserDefined,
	ConstantRate, VariableRate
} RigidGridMotionType_t;
#endif
#define NofValidRigidGridMotionTypes 4

#define ConstantRate_s                 "ConstantRate"
#define VariableRate_s                 "VariableRate"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Arbitrary Grid Motion types                                      *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
        ArbitraryGridMotionTypeNull, ArbitraryGridMotionTypeUserDefined,
        NonDeformingGrid, DeformingGrid
} ArbitraryGridMotionType_t;
#endif
#define NofValidArbitraryGridMotionTypes 4

#define NonDeformingGrid_s             "NonDeformingGrid"
#define DeformingGrid_s                "DeformingGrid"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Simulation types					         *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	SimulationTypeNull, SimulationTypeUserDefined,
	TimeAccurate, NonTimeAccurate
} SimulationType_t;
#endif
#define NofValidSimulationTypes 4

#define TimeAccurate_s                 "TimeAccurate"
#define NonTimeAccurate_s              "NonTimeAccurate"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *	BC Property types						 *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	WallFunctionTypeNull, WallFunctionTypeUserDefined,
	Generic
} WallFunctionType_t;
#endif
#define NofValidWallFunctionTypes 3

#define Generic_s                      "Generic"

#ifdef __CGNS_ENUMS__
typedef enum {
	AreaTypeNull, AreaTypeUserDefined,
	BleedArea, CaptureArea
} AreaType_t;
#endif
#define NofValidAreaTypes 4

#define BleedArea_s                    "BleedArea"
#define CaptureArea_s                  "CaptureArea"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Grid Connectivity Property types				 *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef __CGNS_ENUMS__
typedef enum {
	AverageInterfaceTypeNull, AverageInterfaceTypeUserDefined,
	AverageAll, AverageCircumferential, AverageRadial, AverageI,
	AverageJ, AverageK
} AverageInterfaceType_t;
#endif
#define NofValidAverageInterfaceTypes 8

#define AverageAll_s                   "AverageAll"
#define AverageCircumferential_s       "AverageCircumferential"
#define AverageRadial_s                "AverageRadial"
#define AverageI_s                     "AverageI"
#define AverageJ_s                     "AverageJ"
#define AverageK_s                     "AverageK"



/* The strings defined below are node names or node name patterns */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Coordinate system                                                *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define CoordinateX_s                  "CoordinateX"
#define CoordinateY_s                  "CoordinateY"
#define CoordinateZ_s                  "CoordinateZ"
#define CoordinateR_s                  "CoordinateR"
#define CoordinateTheta_s              "CoordinateTheta"
#define CoordinatePhi_s                "CoordinatePhi"
#define CoordinateNormal_s             "CoordinateNormal"
#define CoordinateTangential_s         "CoordinateTangential"
#define CoordinateXi_s                 "CoordinateXi"
#define CoordinateEta_s                "CoordinateEta"
#define CoordinateZeta_s               "CoordinateZeta"
#define CoordinateTransform_s          "CoordinateTransform"
#define InterpolantsDonor_s            "InterpolantsDonor"
#define ElementConnectivity_s          "ElementConnectivity"
#define ParentData_s                   "ParentData"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      FlowSolution Quantities                                          *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Patterns --- */
#define VectorX_ps                     "%sX"
#define VectorY_ps                     "%sY"
#define VectorZ_ps                     "%sZ"
#define VectorTheta_ps                 "%sTheta"
#define VectorPhi_ps                   "%sPhi"
#define VectorMagnitude_ps             "%sMagnitude"
#define VectorNormal_ps                "%sNormal"
#define VectorTangential_ps            "%sTangential"

#define Potential_s                    "Potential"
#define StreamFunction_s               "StreamFunction"
#define Density_s                      "Density"
#define Pressure_s                     "Pressure"
#define Temperature_s                  "Temperature"
#define EnergyInternal_s               "EnergyInternal"
#define Enthalpy_s                     "Enthalpy"
#define Entropy_s                      "Entropy"
#define EntropyApprox_s                "EntropyApprox"
#define DensityStagnation_s            "DensityStagnation"
#define PressureStagnation_s           "PressureStagnation"
#define TemperatureStagnation_s        "TemperatureStagnation"
#define EnergyStagnation_s             "EnergyStagnation"
#define EnthalpyStagnation_s           "EnthalpyStagnation"
#define EnergyStagnationDensity_s      "EnergyStagnationDensity"
#define VelocityX_s                    "VelocityX"
#define VelocityY_s                    "VelocityY"
#define VelocityZ_s                    "VelocityZ"
#define VelocityR_s                    "VelocityR"
#define VelocityTheta_s                "VelocityTheta"
#define VelocityPhi_s                  "VelocityPhi"
#define VelocityMagnitude_s            "VelocityMagnitude"
#define VelocityNormal_s               "VelocityNormal"
#define VelocityTangential_s           "VelocityTangential"
#define VelocitySound_s                "VelocitySound"
#define VelocitySoundStagnation_s      "VelocitySoundStagnation"
#define MomentumX_s                    "MomentumX"
#define MomentumY_s                    "MomentumY"
#define MomentumZ_s                    "MomentumZ"
#define MomentumMagnitude_s            "MomentumMagnitude"
#define RotatingVelocityX_s            "RotatingVelocityX"
#define RotatingVelocityY_s            "RotatingVelocityY"
#define RotatingVelocityZ_s            "RotatingVelocityZ"
#define RotatingMomentumX_s            "RotatingMomentumX"
#define RotatingMomentumY_s            "RotatingMomentumY"
#define RotatingMomentumZ_s            "RotatingMomentumZ"
#define RotatingVelocityMagnitude_s    "RotatingVelocityMagnitude"
#define RotatingPressureStagnation_s   "RotatingPressureStagnation"
#define RotatingEnergyStagnation_s     "RotatingEnergyStagnation"
#define RotatingEnergyStagnationDensity_s  "RotatingEnergyStagnationDensity"
#define RotatingEnthalpyStagnation_s     "RotatingEnthalpyStagnation"
#define EnergyKinetic_s                "EnergyKinetic"
#define PressureDynamic_s              "PressureDynamic"
#define SoundIntensityDB_s             "SoundIntensityDB"
#define SoundIntensity_s               "SoundIntensity"

#define VorticityX_s                   "VorticityX"
#define VorticityY_s                   "VorticityY"
#define VorticityZ_s                   "VorticityZ"
#define VorticityMagnitude_s           "VorticityMagnitude"
#define SkinFrictionX_s                "SkinFrictionX"
#define SkinFrictionY_s                "SkinFrictionY"
#define SkinFrictionZ_s                "SkinFrictionZ"
#define SkinFrictionMagnitude_s        "SkinFrictionMagnitude"
#define VelocityAngleX_s               "VelocityAngleX"
#define VelocityAngleY_s               "VelocityAngleY"
#define VelocityAngleZ_s               "VelocityAngleZ"
#define VelocityUnitVectorX_s          "VelocityUnitVectorX"
#define VelocityUnitVectorY_s          "VelocityUnitVectorY"
#define VelocityUnitVectorZ_s          "VelocityUnitVectorZ"
#define MassFlow_s                     "MassFlow"
#define ViscosityKinematic_s           "ViscosityKinematic"
#define ViscosityMolecular_s           "ViscosityMolecular"
#define ViscosityEddyDynamic_s         "ViscosityEddyDynamic"
#define ViscosityEddy_s                "ViscosityEddy"
#define ThermalConductivity_s          "ThermalConductivity"
#define PowerLawExponent_s             "PowerLawExponent"
#define SutherlandLawConstant_s        "SutherlandLawConstant"
#define TemperatureReference_s         "TemperatureReference"
#define ViscosityMolecularReference_s  "ViscosityMolecularReference"
#define ThermalConductivityReference_s "ThermalConductivityReference"
#define IdealGasConstant_s             "IdealGasConstant"
#define SpecificHeatPressure_s         "SpecificHeatPressure"
#define SpecificHeatVolume_s           "SpecificHeatVolume"
#define ReynoldsStressXX_s             "ReynoldsStressXX"
#define ReynoldsStressXY_s             "ReynoldsStressXY"
#define ReynoldsStressXZ_s             "ReynoldsStressXZ"
#define ReynoldsStressYY_s             "ReynoldsStressYY"
#define ReynoldsStressYZ_s             "ReynoldsStressYZ"
#define ReynoldsStressZZ_s             "ReynoldsStressZZ"
#define LengthReference_s              "LengthReference"

#define MolecularWeight_s              "MolecularWeight"
#define MolecularWeight_ps             "MolecularWeight%s"
#define HeatOfFormation_s              "HeatOfFormation"
#define HeatOfFormation_ps             "HeatOfFormation%s"
#define FuelAirRatio_s                 "FuelAirRatio"
#define ReferenceTemperatureHOF_s      "ReferenceTemperatureHOF"
#define MassFraction_s                 "MassFraction"
#define MassFraction_ps                "MassFraction%s"
#define LaminarViscosity_s             "LaminarViscosity"
#define LaminarViscosity_ps            "LaminarViscosity%s"
#define ThermalConductivity_ps         "ThermalConductivity%s"
#define EnthalpyEnergyRatio_s          "EnthalpyEnergyRatio"
#define CompressibilityFactor_s        "CompressibilityFactor"
#define VibrationalElectronEnergy_s    "VibrationalElectronEnergy"
#define VibrationalElectronTemperature_s  "VibrationalElectronTemperature"
#define SpeciesDensity_s               "SpeciesDensity"
#define SpeciesDensity_ps              "SpeciesDensity%s"
#define MoleFraction_s                 "MoleFraction"
#define MoleFraction_ps                "MoleFraction%s"

#define ElectricFieldX_s               "ElectricFieldX"
#define ElectricFieldY_s               "ElectricFieldY"
#define ElectricFieldZ_s               "ElectricFieldZ"
#define MagneticFieldX_s               "MagneticFieldX"
#define MagneticFieldY_s               "MagneticFieldY"
#define MagneticFieldZ_s               "MagneticFieldZ"
#define CurrentDensityX_s              "CurrentDensityX"
#define CurrentDensityY_s              "CurrentDensityY"
#define CurrentDensityZ_s              "CurrentDensityZ"
#define LorentzForceX_s                "LorentzForceX"
#define LorentzForceY_s                "LorentzForceY"
#define LorentzForceZ_s                "LorentzForceZ"
#define ElectricConductivity_s         "ElectricConductivity"
#define JouleHeating_s                 "JouleHeating"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Typical Turbulence Models                                        *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define TurbulentDistance_s            "TurbulentDistance"
#define TurbulentEnergyKinetic_s       "TurbulentEnergyKinetic"
#define TurbulentDissipation_s         "TurbulentDissipation"
#define TurbulentDissipationRate_s     "TurbulentDissipationRate"
#define TurbulentBBReynolds_s          "TurbulentBBReynolds"
#define TurbulentSANuTilde_s           "TurbulentSANuTilde"


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Nondimensional Parameters                                        *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define Mach_s                         "Mach"
#define Mach_Velocity_s                "Mach_Velocity"
#define Mach_VelocitySound_s           "Mach_VelocitySound"
#define Reynolds_s                     "Reynolds"
#define Reynolds_Velocity_s            "Reynolds_Velocity"
#define Reynolds_Length_s              "Reynolds_Length"
#define Reynolds_ViscosityKinematic_s  "Reynolds_ViscosityKinematic"
#define Prandtl_s                      "Prandtl"
#define Prandtl_ThermalConductivity_s  "Prandtl_ThermalConductivity"
#define Prandtl_ViscosityMolecular_s   "Prandtl_ViscosityMolecular"
#define Prandtl_SpecificHeatPressure_s "Prandtl_SpecificHeatPressure"
#define PrandtlTurbulent_s             "PrandtlTurbulent"
#define SpecificHeatRatio_s            "SpecificHeatRatio"
#define SpecificHeatRatio_Pressure_s   "SpecificHeatRatio_Pressure"
#define SpecificHeatRatio_Volume_s     "SpecificHeatRatio_Volume"
#define CoefPressure_s                 "CoefPressure"
#define CoefSkinFrictionX_s            "CoefSkinFrictionX"
#define CoefSkinFrictionY_s            "CoefSkinFrictionY"
#define CoefSkinFrictionZ_s            "CoefSkinFrictionZ"
#define Coef_PressureDynamic_s         "Coef_PressureDynamic"
#define Coef_PressureReference_s       "Coef_PressureReference"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Characteristics and Riemann invariant                            *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define Vorticity_s                    "Vorticity"
#define Acoustic_s                     "Acoustic"

#define RiemannInvariantPlus_s         "RiemannInvariantPlus"
#define RiemannInvariantMinus_s        "RiemannInvariantMinus"
#define CharacteristicEntropy_s        "CharacteristicEntropy"
#define CharacteristicVorticity1_s     "CharacteristicVorticity1"
#define CharacteristicVorticity2_s     "CharacteristicVorticity2"
#define CharacteristicAcousticPlus_s   "CharacteristicAcousticPlus"
#define CharacteristicAcousticMinus_s  "CharacteristicAcousticMinus"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Forces and Moments                                               *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define ForceX_s                       "ForceX"
#define ForceY_s                       "ForceY"
#define ForceZ_s                       "ForceZ"
#define ForceR_s                       "ForceR"
#define ForceTheta_s                   "ForceTheta"
#define ForcePhi_s                     "ForcePhi"
#define Lift_s                         "Lift"
#define Drag_s                         "Drag"
#define MomentX_s                      "MomentX"
#define MomentY_s                      "MomentY"
#define MomentZ_s                      "MomentZ"
#define MomentR_s                      "MomentR"
#define MomentTheta_s                  "MomentTheta"
#define MomentPhi_s                    "MomentPhi"
#define MomentXi_s                     "MomentXi"
#define MomentEta_s                    "MomentEta"
#define MomentZeta_s                   "MomentZeta"
#define Moment_CenterX_s               "Moment_CenterX"
#define Moment_CenterY_s               "Moment_CenterY"
#define Moment_CenterZ_s               "Moment_CenterZ"
#define CoefLift_s                     "CoefLift"
#define CoefDrag_s                     "CoefDrag"
#define CoefMomentX_s                  "CoefMomentX"
#define CoefMomentY_s                  "CoefMomentY"
#define CoefMomentZ_s                  "CoefMomentZ"
#define CoefMomentR_s                  "CoefMomentR"
#define CoefMomentTheta_s              "CoefMomentTheta"
#define CoefMomentPhi_s                "CoefMomentPhi"
#define CoefMomentXi_s                 "CoefMomentXi"
#define CoefMomentEta_s                "CoefMomentEta"
#define CoefMomentZeta_s               "CoefMomentZeta"
#define Coef_PressureDynamic_s         "Coef_PressureDynamic"
#define Coef_Area_s                    "Coef_Area"
#define Coef_Length_s                  "Coef_Length"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *       Time dependent flow                                             *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define TimeValues_s                   "TimeValues"
#define IterationValues_s              "IterationValues"
#define NumberOfZones_s                "NumberOfZones"
#define NumberOfFamilies_s             "NumberOfFamilies"
#define ZonePointers_s                 "ZonePointers"
#define FamilyPointers_s               "FamilyPointers"
#define RigidGridMotionPointers_s      "RigidGridMotionPointers"
#define ArbitraryGridMotionPointers_s  "ArbitraryGridMotionPointers"
#define GridCoordinatesPointers_s      "GridCoordinatesPointers"
#define FlowSolutionPointers_s         "FlowSolutionPointers"
#define OriginLocation_s               "OriginLocation"
#define RigidRotationAngle_s           "RigidRotationAngle"
#define RigidVelocity_s                "RigidVelocity"
#define RigidRotationRate_s            "RigidRotationRate"
#define GridVelocityX_s                "GridVelocityX"
#define GridVelocityY_s                "GridVelocityY"
#define GridVelocityZ_s                "GridVelocityZ"
#define GridVelocityR_s                "GridVelocityR"
#define GridVelocityTheta_s            "GridVelocityTheta"
#define GridVelocityPhi_s              "GridVelocityPhi"
#define GridVelocityXi_s               "GridVelocityXi"
#define GridVelocityEta_s              "GridVelocityEta"
#define GridVelocityZeta_s             "GridVelocityZeta"


/* The strings defined below are type names used for node labels */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *       Types as strings                                                *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#define ArbitraryGridMotion_ts         "ArbitraryGridMotion_t"
#define Area_ts                        "Area_t"
#define AverageInterface_ts            "AverageInterface_t"
#define Axisymmetry_ts                 "Axisymmetry_t"
#define BCDataSet_ts                   "BCDataSet_t"
#define BCData_ts                      "BCData_t"
#define BCProperty_ts                  "BCProperty_t"
#define BC_ts                          "BC_t"
#define BaseIterativeData_ts           "BaseIterativeData_t"
#define CGNSBase_ts                    "CGNSBase_t"
#define CGNSLibraryVersion_ts          "CGNSLibraryVersion_t"
#define ChemicalKineticsModel_ts       "ChemicalKineticsModel_t"
#define ConvergenceHistory_ts          "ConvergenceHistory_t"
#define DataArray_ts                   "DataArray_t"
#define DataClass_ts                   "DataClass_t"
#define DataConversion_ts              "DataConversion_t"
#define Descriptor_ts                  "Descriptor_t"
#define DimensionalExponents_ts        "DimensionalExponents_t"
#define DimensionalUnits_ts            "DimensionalUnits_t"   
#define DiscreteData_ts                "DiscreteData_t"
#define Elements_ts                    "Elements_t"
#define FamilyBC_ts                    "FamilyBC_t"
#define FamilyName_ts                  "FamilyName_t"
#define Family_ts                      "Family_t"
#define FlowEquationSet_ts             "FlowEquationSet_t"
#define FlowSolution_ts                "FlowSolution_t"
#define GasModel_ts                    "GasModel_t"
#define GeometryEntity_ts              "GeometryEntity_t"
#define GeometryFile_ts                "GeometryFile_t"
#define GeometryFormat_ts              "GeometryFormat_t"
#define GeometryReference_ts           "GeometryReference_t"
#define GoverningEquations_ts          "GoverningEquations_t"
#define Gravity_ts                     "Gravity_t"
#define GridConnectivity1to1_ts        "GridConnectivity1to1_t"
#define GridConnectivityProperty_ts    "GridConnectivityProperty_t"
#define GridConnectivityType_ts        "GridConnectivityType_t"
#define GridConnectivity_ts            "GridConnectivity_t"
#define GridCoordinates_ts             "GridCoordinates_t"
#define GridLocation_ts                "GridLocation_t"
#define IndexArray_ts                  "IndexArray_t"
#define IndexRange_ts                  "IndexRange_t"   
#define IntegralData_ts                "IntegralData_t"
#define InwardNormalList_ts            "InwardNormalList_t"
#define Ordinal_ts                     "Ordinal_t"
#define OversetHoles_ts                "OversetHoles_t"
#define Periodic_ts                    "Periodic_t"
#define ReferenceState_ts              "ReferenceState_t"
#define RigidGridMotion_ts             "RigidGridMotion_t"
#define Rind_ts                        "Rind_t"   
#define RotatingCoordinates_ts         "RotatingCoordinates_t"
#define SimulationType_ts              "SimulationType_t"
#define ThermalConductivityModel_ts    "ThermalConductivityModel_t"
#define ThermalRelaxationModel_ts      "ThermalRelaxationModel_t"
#define TurbulenceClosure_ts           "TurbulenceClosure_t"
#define TurbulenceModel_ts             "TurbulenceModel_t"
#define UserDefinedData_ts             "UserDefinedData_t"
#define ViscosityModel_ts              "ViscosityModel_t"
#define WallFunction_ts                "WallFunction_t"
#define ZoneBC_ts                      "ZoneBC_t"
#define ZoneGridConnectivity_ts        "ZoneGridConnectivity_t"
#define ZoneIterativeData_ts           "ZoneIterativeData_t"
#define ZoneType_ts                    "ZoneType_t"
#define Zone_ts                        "Zone_t"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *       No line after this comment -                                    *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#endif
