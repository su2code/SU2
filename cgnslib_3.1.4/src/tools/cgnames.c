#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cgnslib.h"

/*=======================================================================*/

#if CGNS_VERSION < 2500

static char const *InvalidName = "<invalid>";
static char const *NotImplemented = "<not implemented>";

char const *cg_MassUnitsName (int type)
{
    if (type < 0 || type >= NofValidMassUnits)
        return InvalidName;
    return MassUnitsName[type];
}

char const *cg_LengthUnitsName (int type)
{
    if (type < 0 || type >= NofValidLengthUnits)
        return InvalidName;
    return LengthUnitsName[type];
}

char const *cg_TimeUnitsName (int type)
{
    if (type < 0 || type >= NofValidTimeUnits)
        return InvalidName;
    return TimeUnitsName[type];
}

char const *cg_TemperatureUnitsName (int type)
{
    if (type < 0 || type >= NofValidTemperatureUnits)
        return InvalidName;
    return TemperatureUnitsName[type];
}

char const *cg_AngleUnitsName (int type)
{
    if (type < 0 || type >= NofValidAngleUnits)
        return InvalidName;
    return AngleUnitsName[type];
}

char const *cg_ElectricCurrentUnitsName (int type)
{
#ifdef NofValidElectricCurrentUnits
    if (type < 0 || type >= NofValidElectricCurrentUnits)
        return InvalidName;
    return ElectricCurrentUnitsName[type];
#else
    return NotImplemented;
#endif
}

char const *cg_SubstanceAmountUnitsName (int type)
{
#ifdef NofValidSubstanceAmountUnits
    if (type < 0 || type >= NofValidSubstanceAmountUnits)
        return InvalidName;
    return SubstanceAmountUnitsName[type];
#else
    return NotImplemented;
#endif
}

char const *cg_LuminousIntensityUnitsName (int type)
{
#ifdef NofValidLuminousIntensityUnits
    if (type < 0 || type >= NofValidLuminousIntensityUnits)
        return InvalidName;
    return LuminousIntensityUnitsName[type];
#else
    return NotImplemented;
#endif
}

char const *cg_DataClassName (int type)
{
    if (type < 0 || type >= NofValidDataClass)
        return InvalidName;
    return DataClassName[type];
}

char const *cg_GridLocationName (int type)
{
    if (type < 0 || type >= NofValidGridLocation)
        return InvalidName;
    return GridLocationName[type];
}

char const *cg_BCDataTypeName (int type)
{
    if (type < 0 || type >= NofValidBCDataTypes)
        return InvalidName;
    return BCDataTypeName[type];
}

char const *cg_GridConnectivityTypeName (int type)
{
    if (type < 0 || type >= NofValidGridConnectivityTypes)
        return InvalidName;
    return GridConnectivityTypeName[type];
}

char const *cg_PointSetTypeName (int type)
{
    if (type < 0 || type >= NofValidPointSetTypes)
        return InvalidName;
    return PointSetTypeName[type];
}

char const *cg_GoverningEquationsTypeName (int type)
{
    if (type < 0 || type >= NofValidGoverningEquationsTypes)
        return InvalidName;
    return GoverningEquationsTypeName[type];
}

char const *cg_ModelTypeName (int type)
{
    if (type < 0 || type >= NofValidModelTypes)
        return InvalidName;
    return ModelTypeName[type];
}

char const *cg_BCTypeName (int type)
{
    if (type < 0 || type >= NofValidBCTypes)
        return InvalidName;
    return BCTypeName[type];
}

char const *cg_DataTypeName (int type)
{
    if (type < 0 || type >= NofValidDataTypes)
        return InvalidName;
    return DataTypeName[type];
}

char const *cg_ElementTypeName (int type)
{
    if (type < 0 || type >= NofValidElementTypes)
        return InvalidName;
    return ElementTypeName[type];
}

char const *cg_ZoneTypeName (int type)
{
    if (type < 0 || type >= NofValidZoneTypes)
        return InvalidName;
    return ZoneTypeName[type];
}

char const *cg_RigidGridMotionTypeName (int type)
{
#ifdef NofValidRigidGridMotionTypes
    if (type < 0 || type >= NofValidRigidGridMotionTypes)
        return InvalidName;
    return RigidGridMotionTypeName[type];
#else
    return NotImplemented;
#endif
}

char const *cg_ArbitraryGridMotionTypeName (int type)
{
#ifdef NofValidArbitraryGridMotionTypes
    if (type < 0 || type >= NofValidArbitraryGridMotionTypes)
        return InvalidName;
    return ArbitraryGridMotionTypeName[type];
#else
    return NotImplemented;
#endif
}

char const *cg_SimulationTypeName (int type)
{
#ifdef NofValidSimulationTypes
    if (type < 0 || type >= NofValidSimulationTypes)
        return InvalidName;
    return SimulationTypeName[type];
#else
    return NotImplemented;
#endif
}

char const *cg_WallFunctionTypeName (int type)
{
#ifdef NofValidWallFunctionTypes
    if (type < 0 || type >= NofValidWallFunctionTypes)
        return InvalidName;
    return WallFunctionTypeName[type];
#else
    return NotImplemented;
#endif
}

char const *cg_AreaTypeName (int type)
{
#ifdef NofValidAreaTypes
    if (type < 0 || type >= NofValidAreaTypes)
        return InvalidName;
    return AreaTypeName[type];
#else
    return NotImplemented;
#endif
}

char const *cg_AverageInterfaceTypeName (int type)
{
#ifdef NofValidAverageInterfaceTypes
    if (type < 0 || type >= NofValidAverageInterfaceTypes)
        return InvalidName;
    return AverageInterfaceTypeName[type];
#else
    return NotImplemented;
#endif
}

#endif

/*======================================================================*/

typedef struct {
    char *name;
    int flags;
    int nexps;
    int exps[8];
} IDENTIFIER;

/* exponents are: Mass, Length, Time, Temperature, Angle,
                  ElectricCurrent, SubstanceAmount, LuminousIntensity */

static IDENTIFIER Identifier[] = {
{"AxisymmetryAngle",                 0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"AxisymmetryAxisVector",            0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"AxisymmetryReferencePoint",        0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"CharacteristicAcousticMinus",      0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"CharacteristicAcousticPlus",       0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"CharacteristicEntropy",            0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"CharacteristicVorticity1",         0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"CharacteristicVorticity2",         0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"CoefDrag",                         0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefLift",                         0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentEta",                    0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentMagnitude",              0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentNormal",                 0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentPhi",                    0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentR",                      0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentTangential",             0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentTheta",                  0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentX",                      0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentXi",                     0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentY",                      0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentZ",                      0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefMomentZeta",                   0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefPressure",                     0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefSkinFrictionMagnitude",        0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefSkinFrictionNormal",           0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefSkinFrictionPhi",              0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefSkinFrictionR",                0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefSkinFrictionTangential",       0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefSkinFrictionTheta",            0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefSkinFrictionX",                0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefSkinFrictionY",                0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoefSkinFrictionZ",                0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"Coef_Area",                        0, 5, { 0,  2,  0,  0,  0,  0,  0,  0}},
{"Coef_Length",                      0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"Coef_PressureDynamic",             0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"Coef_PressureReference",           0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"CompressibilityFactor",            0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"CoordinateEta",                    0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"CoordinateNormal",                 0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"CoordinatePhi",                    0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"CoordinateR",                      0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"CoordinateTangential",             0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"CoordinateTheta",                  0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"CoordinateX",                      0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"CoordinateXi",                     0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"CoordinateY",                      0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"CoordinateZ",                      0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"CoordinateZeta",                   0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"CurrentDensityX",                  0, 8, { 0, -2,  0,  0,  0,  1,  0,  0}},
{"CurrentDensityY",                  0, 8, { 0, -2,  0,  0,  0,  1,  0,  0}},
{"CurrentDensityZ",                  0, 8, { 0, -2,  0,  0,  0,  1,  0,  0}},
{"Density",                          0, 5, { 1, -3,  0,  0,  0,  0,  0,  0}},
{"DensityStagnation",                0, 5, { 1, -3,  0,  0,  0,  0,  0,  0}},
{"Drag",                             0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"ElectricConductivity",             0, 8, { 1,  1, -3,  0,  0, -2,  0,  0}},
/* electric field (volts/meter) */
{"ElectricFieldX",                   0, 8, { 1,  1, -1,  0,  0, -1,  0,  0}},
{"ElectricFieldY",                   0, 8, { 1,  1, -1,  0,  0, -1,  0,  0}},
{"ElectricFieldZ",                   0, 8, { 1,  1, -1,  0,  0, -1,  0,  0}},
{"EnergyInternal",                   0, 5, { 0,  2, -2,  0,  0,  0,  0,  0}},
{"EnergyKinetic",                    0, 5, { 0,  2, -2,  0,  0,  0,  0,  0}},
{"EnergyStagnation",                 0, 5, { 0,  2, -2,  0,  0,  0,  0,  0}},
{"EnergyStagnationDensity",          0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"Enthalpy",                         0, 5, { 0,  2, -2,  0,  0,  0,  0,  0}},
{"EnthalpyEnergyRatio",              0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"EnthalpyStagnation",               0, 5, { 0,  2, -2,  0,  0,  0,  0,  0}},
{"Entropy",                          0, 5, { 1,  2, -2, -1,  0,  0,  0,  0}},
 /* the length and mass exponents depend on gamma */
{"EntropyApprox",                    0,-5, {-1,  1, -2,  0,  0,  0,  0,  0}},
{"ForceMagnitude",                   0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"ForceNormal",                      0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"ForcePhi",                         0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"ForceR",                           0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"ForceTangential",                  0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"ForceTheta",                       0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"ForceX",                           0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"ForceY",                           0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"ForceZ",                           0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"FuelAirRatio",                     0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"GravityVector",                    0, 5, { 0,  1, -2,  0,  0,  0,  0,  0}},
{"GridVelocityEta",                  0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"GridVelocityMagnitude",            0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"GridVelocityNormal",               0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"GridVelocityPhi",                  0, 5, { 0,  0, -1,  0,  1,  0,  0,  0}},
{"GridVelocityR",                    0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"GridVelocityTangential",           0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"GridVelocityTheta",                0, 5, { 0,  0, -1,  0,  1,  0,  0,  0}},
{"GridVelocityX",                    0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"GridVelocityXi",                   0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"GridVelocityY",                    0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"GridVelocityZ",                    0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"GridVelocityZeta",                 0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"HeatOfFormation",                  1, 5, { 0,  2, -2,  0,  0,  0,  0,  0}},
{"IdealGasConstant",                 0, 5, { 0,  2, -2, -1,  0,  0,  0,  0}},
{"JouleHeating",                     0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"LaminarViscosity",                 1, 5, { 1, -1, -1,  0,  0,  0,  0,  0}},
{"LengthReference",                  0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"Lift",                             0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"LorenzForceX",                     0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"LorenzForceY",                     0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"LorenzForceZ",                     0, 5, { 1,  1, -2,  0,  0,  0,  0,  0}},
{"Mach",                             0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"Mach_Velocity",                    0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"Mach_VelocitySound",               0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
/* magnetic field strength (tesla)
{"MagneticFieldX",                   0, 8, { 1,  0, -2,  0,  0, -1,  0,  0}},
{"MagneticFieldY",                   0, 8, { 1,  0, -2,  0,  0, -1,  0,  0}},
{"MagneticFieldZ",                   0, 8, { 1,  0, -2,  0,  0, -1,  0,  0}},
   magnetic field (Amperes/meter) */
{"MagneticFieldX",                   0, 8, { 0, -1,  0,  0,  0,  1,  0,  0}},
{"MagneticFieldY",                   0, 8, { 0, -1,  0,  0,  0,  1,  0,  0}},
{"MagneticFieldZ",                   0, 8, { 0, -1,  0,  0,  0,  1,  0,  0}},
{"MassFlow",                         0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"MassFraction",                     1, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"MoleFraction",                     1, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"MolecularWeight",                  1, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"MomentEta",                        0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"MomentMagnitude",                  0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"MomentNormal",                     0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"MomentPhi",                        0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"MomentR",                          0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"MomentTangential",                 0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"MomentTheta",                      0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"MomentX",                          0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"MomentXi",                         0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"MomentY",                          0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"MomentZ",                          0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"MomentZeta",                       0, 5, { 1,  2, -2,  0,  0,  0,  0,  0}},
{"Moment_CenterNormal",              0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"Moment_CenterPhi",                 0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"Moment_CenterR",                   0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"Moment_CenterTangential",          0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"Moment_CenterTheta",               0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"Moment_CenterX",                   0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"Moment_CenterY",                   0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"Moment_CenterZ",                   0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"MomentumMagnitude",                0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"MomentumNormal",                   0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"MomentumPhi",                      0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"MomentumR",                        0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"MomentumTangential",               0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"MomentumTheta",                    0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"MomentumX",                        0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"MomentumY",                        0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"MomentumZ",                        0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"OriginLocation",                   0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"Potential",                        0, 5, { 0,  2, -1,  0,  0,  0,  0,  0}},
{"PowerLawExponent",                 0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"Prandtl",                          0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"PrandtlTurbulent",                 0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"Prandtl_SpecificHeatPressure",     0, 5, { 0,  2, -2, -1,  0,  0,  0,  0}},
{"Prandtl_ThermalConductivity",      0, 5, { 1,  1, -3, -1,  0,  0,  0,  0}},
{"Prandtl_ViscosityMolecular",       0, 5, { 1, -1, -1,  0,  0,  0,  0,  0}},
{"Pressure",                         0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"PressureDynamic",                  0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"PressureStagnation",               0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"ReferenceTemperatureHOF",          0, 5, { 0,  0,  0,  1,  0,  0,  0,  0}},
{"Reynolds",                         0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"ReynoldsStressXX",                 0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"ReynoldsStressXY",                 0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"ReynoldsStressXZ",                 0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"ReynoldsStressYY",                 0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"ReynoldsStressYZ",                 0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"ReynoldsStressZZ",                 0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"Reynolds_Length",                  0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"Reynolds_Velocity",                0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"Reynolds_ViscosityKinematic",      0, 5, { 0,  2, -1,  0,  0,  0,  0,  0}},
{"RiemannInvariantMinus",            0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RiemannInvariantPlus",             0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RigidRotationAngle",               0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"RigidRotationRate",                0, 5, { 0,  0, -1,  0,  1,  0,  0,  0}},
{"RigidVelocity",                    0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RotatingEnergyStagnation",         0, 5, { 0,  2, -2,  0,  0,  0,  0,  0}},
{"RotatingEnergyStagnationDensity",  0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"RotatingEnthalpyStagnation",       0, 5, { 0,  2, -2,  0,  0,  0,  0,  0}},
{"RotatingMach",                     0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"RotatingMomentumMagnitude",        0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"RotatingMomentumNormal",           0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"RotatingMomentumPhi",              0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"RotatingMomentumR",                0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"RotatingMomentumTangential",       0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"RotatingMomentumTheta",            0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"RotatingMomentumX",                0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"RotatingMomentumY",                0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"RotatingMomentumZ",                0, 5, { 1, -2, -1,  0,  0,  0,  0,  0}},
{"RotatingPressureStagnation",       0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"RotatingVelocityMagnitude",        0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RotatingVelocityNormal",           0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RotatingVelocityPhi",              0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RotatingVelocityR",                0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RotatingVelocityTangential",       0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RotatingVelocityTheta",            0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RotatingVelocityX",                0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RotatingVelocityY",                0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RotatingVelocityZ",                0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"RotationAngle",                    0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"RotationCenter",                   0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"RotationRateVector",               0, 5, { 0,  0, -1,  0,  1,  0,  0,  0}},
/* sigma - surface charge density (coulombs/square meter) ?
{"Sigma",                            0, 8, { 0, -2,  1,  0,  0,  1,  0,  0}},
*/
{"SkinFrictionMagnitude",            0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"SkinFrictionNormal",               0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"SkinFrictionPhi",                  0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"SkinFrictionR",                    0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"SkinFrictionTangential",           0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"SkinFrictionTheta",                0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"SkinFrictionX",                    0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"SkinFrictionY",                    0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"SkinFrictionZ",                    0, 5, { 1, -1, -2,  0,  0,  0,  0,  0}},
{"SoundIntensity",                   0, 5, { 1,  0, -3,  0,  0,  0,  0,  0}},
{"SoundIntensityDB",                 0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"SpeciesDensity",                   1, 5, { 1, -3,  0,  0,  0,  0,  0,  0}},
{"SpecificHeatPressure",             0, 5, { 0,  2, -2, -1,  0,  0,  0,  0}},
{"SpecificHeatRatio",                0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"SpecificHeatRatio_Pressure",       0, 5, { 0,  2, -2, -1,  0,  0,  0,  0}},
{"SpecificHeatRatio_Volume",         0, 5, { 0,  2, -2, -1,  0,  0,  0,  0}},
{"SpecificHeatVolume",               0, 5, { 0,  2, -2, -1,  0,  0,  0,  0}},
{"StreamFunction",                   0, 5, { 0,  2, -1,  0,  0,  0,  0,  0}},
{"SurfaceArea",                      0, 5, { 0,  2,  0,  0,  0,  0,  0,  0}},
{"SutherlandLawConstant",            0, 5, { 0,  0,  0,  1,  0,  0,  0,  0}},
{"Temperature",                      0, 5, { 0,  0,  0,  1,  0,  0,  0,  0}},
{"TemperatureReference",             0, 5, { 0,  0,  0,  1,  0,  0,  0,  0}},
{"TemperatureStagnation",            0, 5, { 0,  0,  0,  1,  0,  0,  0,  0}},
{"ThermalConductivity",              1, 5, { 1,  1, -3, -1,  0,  0,  0,  0}},
{"ThermalConductivityReference",     0, 5, { 1,  1, -3, -1,  0,  0,  0,  0}},
{"Translation",                      0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"TurbulentBBReynolds",              0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"TurbulentDissipation",             0, 5, { 0,  2, -3,  0,  0,  0,  0,  0}},
{"TurbulentDissipationRate",         0, 5, { 0,  0, -1,  0,  0,  0,  0,  0}},
{"TurbulentDistance",                0, 5, { 0,  1,  0,  0,  0,  0,  0,  0}},
{"TurbulentEnergyKinetic",           0, 5, { 0,  2, -2,  0,  0,  0,  0,  0}},
{"TurbulentSANuTilde",               0, 5, { 0,  2, -1,  0,  0,  0,  0,  0}},
{"VelocityAngleMagnitude",           0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"VelocityAngleNormal",              0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"VelocityAnglePhi",                 0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"VelocityAngleR",                   0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"VelocityAngleTangential",          0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"VelocityAngleTheta",               0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"VelocityAngleX",                   0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"VelocityAngleY",                   0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"VelocityAngleZ",                   0, 5, { 0,  0,  0,  0,  1,  0,  0,  0}},
{"VelocityMagnitude",                0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"VelocityNormal",                   0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"VelocityPhi",                      0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"VelocityR",                        0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"VelocitySound",                    0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"VelocitySoundStagnation",          0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"VelocityTangential",               0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"VelocityTheta",                    0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"VelocityUnitVectorMagnitude",      0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"VelocityUnitVectorNormal",         0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"VelocityUnitVectorPhi",            0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"VelocityUnitVectorR",              0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"VelocityUnitVectorTangential",     0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"VelocityUnitVectorTheta",          0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"VelocityUnitVectorX",              0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"VelocityUnitVectorY",              0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"VelocityUnitVectorZ",              0, 0, { 0,  0,  0,  0,  0,  0,  0,  0}},
{"VelocityX",                        0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"VelocityY",                        0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"VelocityZ",                        0, 5, { 0,  1, -1,  0,  0,  0,  0,  0}},
{"VibrationalElectronEnergy",        0, 5, { 0,  2, -2,  0,  0,  0,  0,  0}},
{"VibrationalElectronTemperature",   0, 5, { 0,  0,  0,  1,  0,  0,  0,  0}},
{"ViscosityEddy",                    0, 5, { 1, -1, -1,  0,  0,  0,  0,  0}},
{"ViscosityEddyKinematic",           0, 5, { 0,  2, -1,  0,  0,  0,  0,  0}},
{"ViscosityKinematic",               0, 5, { 0,  2, -1,  0,  0,  0,  0,  0}},
{"ViscosityMolecular",               0, 5, { 1, -1, -1,  0,  0,  0,  0,  0}},
{"ViscosityMolecularReference",      0, 5, { 1, -1, -1,  0,  0,  0,  0,  0}},
{"Voltage",                          0, 8, { 1,  2, -1,  0,  0, -1,  0,  0}},
{"VorticityMagnitude",               0, 5, { 0,  0, -1,  0,  0,  0,  0,  0}},
{"VorticityNormal",                  0, 5, { 0,  0, -1,  0,  0,  0,  0,  0}},
{"VorticityPhi",                     0, 5, { 0,  0, -1,  0,  0,  0,  0,  0}},
{"VorticityR",                       0, 5, { 0,  0, -1,  0,  0,  0,  0,  0}},
{"VorticityTangential",              0, 5, { 0,  0, -1,  0,  0,  0,  0,  0}},
{"VorticityTheta",                   0, 5, { 0,  0, -1,  0,  0,  0,  0,  0}},
{"VorticityX",                       0, 5, { 0,  0, -1,  0,  0,  0,  0,  0}},
{"VorticityY",                       0, 5, { 0,  0, -1,  0,  0,  0,  0,  0}},
{"VorticityZ",                       0, 5, { 0,  0, -1,  0,  0,  0,  0,  0}}
};

#define NUM_IDENTIFIER (sizeof(Identifier)/sizeof(IDENTIFIER))

/*-----------------------------------------------------------------------*/

int cg_get_identifier (const char *name, int *nexps, float *exps)
{
    int n, cmp, lo = 0, hi = NUM_IDENTIFIER - 1, mid;
    IDENTIFIER *ident = NULL;

    if (NULL == name || !*name) return 1;
    if (0 == strcmp (Identifier[lo].name, name))
        ident = &Identifier[lo];
    else if (0 == strcmp (Identifier[hi].name, name))
        ident = &Identifier[hi];
    else {
        while (lo <= hi) {
            mid = (lo + hi) >> 1;
            if ((Identifier[mid].flags & 1) == 1)
                cmp = strncmp (Identifier[mid].name, name,
                               strlen(Identifier[mid].name));
            else
                cmp = strcmp (Identifier[mid].name, name);
            if (0 == cmp) {
                ident = &Identifier[mid];
                break;
            }
            if (cmp > 0)
                hi = mid - 1;
            else
                lo = mid + 1;
        }
    }

    if (ident == NULL) return 1;

    *nexps = ident->nexps;
    if (ident->nexps && exps != NULL) {
        cmp = abs (ident->nexps);
        for (n = 0; n < cmp; n++)
            exps[n] = (float)ident->exps[n];
    }
    return 0;
}

/*---------------------------------------------------------------------*/

static int matches (char *p, char *s)
{
    char *cmp;
    int rev, strt, n;

    while (*p) {
        if (*s == '\\')
            s++;
        switch (*p) {

            /* match any single character */

            case '?':
                if (!*s++)
                    return 0;
                break;

            /* match 0 or more characters */

            case '*':
                if (!*++p)
                    return 1;
                while (*s) {
                    if ((n = matches (p, s)) != 0)
                        break;
                    s++;
                }
                return (*s ? n : 0);

            /* match list of characters */

            case '[':
                if (*++p == '^') {
                    rev = 1;
                    p++;
                }
                else
                    rev = 0;
                for (cmp = p; *cmp != ']'; cmp++) {
                    if (!*cmp)
                        return -1;
                    if (*cmp == '-') {
                        strt = *(cmp-1);
                        if (cmp == p || !*++cmp || (*cmp == '\\' && !*++cmp))
                            return -1;
                        if (((strt <= *s && *s <= *cmp ? 1 : 0) ^ rev) != 0)
                            break;
                    }
                    if (*cmp == '\\' && !*++cmp)
                        return -1;
                    if (((*s == *cmp ? 1 : 0) ^ rev) != 0)
                        break;
                }
                if (*cmp == ']')
                    return 0;
                while (*cmp && (*cmp != ']' || *(cmp-1) == '\\'))
                    cmp++;
                if (!*cmp)
                    return 0;
                p = cmp;
                s++;
                break;

            /* single character */

            case '\\':
                if (!*++p)
                    return -1;
            default:
                if (*s++ != *p)
                    return 0;
                break;
        }
        p++;
    }
    return (*s ? 0 : 1);
}

/*-----------------------------------------------------------------------*/

int cg_find_identifier (const char *pattern, int *nnames, char ***names)
{
    int n, cnt = 0;
    char *p, **pp;

    for (n = 0; n < NUM_IDENTIFIER; n++) {
        if (matches ((char *)pattern, Identifier[n].name) > 0) {
            Identifier[n].flags |= 8;
            cnt++;
        }
        else
            Identifier[n].flags &= 7;
    }
    *nnames = cnt;
    if (!cnt) return 0;

    if (names != NULL) {
        pp = (char **) malloc (cnt * (33 + sizeof(char *)));
        if (pp == NULL) {
            fprintf (stderr, "malloc failed for indentifier names\n");
            return 1;
        }
        p = (char *)(pp + cnt);
        for (cnt = 0, n = 0; n < NUM_IDENTIFIER; n++) {
            if ((Identifier[n].flags & 8) == 8) {
                if ((Identifier[n].flags & 1) == 1)
                    sprintf (p, "%s#", Identifier[n].name);
                else
                    strcpy (p, Identifier[n].name);
                pp[cnt++] = p;
                p += 33;
            }
        }
    }

    return 0;
}

/*-----------------------------------------------------------------------*/

int cg_enum_identifier (int (*callback)(char *name,
    int nexps, float *exps, void *user), void *user)
{
    int n, i, ierr;
    float exps[8];
    char name[33];

    for (n = 0; n < NUM_IDENTIFIER; n++) {
        for (i = 0; i < 8; i++)
            exps[i] = (float)Identifier[n].exps[i];
        if ((Identifier[n].flags & 1) == 1)
            sprintf (name, "%s#", Identifier[n].name);
        else
            strcpy (name, Identifier[n].name);
        ierr = (*callback) (name, Identifier[n].nexps, exps, user);
        if (ierr) return ierr;
    }
    return 0;
}

/*-----------------------------------------------------------------------*/

#ifdef STANDALONE

/* Removed by KSH 2009.05.15 */

#endif

