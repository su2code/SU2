#if CGNS_VERSION < 2500
char const *cg_MassUnitsName (int);
char const *cg_LengthUnitsName (int);
char const *cg_TimeUnitsName (int);
char const *cg_TemperatureUnitsName (int);
char const *cg_AngleUnitsName (int);
char const *cg_ElectricCurrentUnitsName (int);
char const *cg_SubstanceAmountUnitsName (int);
char const *cg_LuminousIntensityUnitsName (int);
char const *cg_DataClassName (int);
char const *cg_GridLocationName (int);
char const *cg_BCDataTypeName (int);
char const *cg_GridConnectivityTypeName (int);
char const *cg_PointSetTypeName (int);
char const *cg_GoverningEquationsTypeName (int);
char const *cg_ModelTypeName (int);
char const *cg_BCTypeName (int);
char const *cg_DataTypeName (int);
char const *cg_ElementTypeName (int);
char const *cg_ZoneTypeName (int);
char const *cg_RigidGridMotionTypeName (int);
char const *cg_ArbitraryGridMotionTypeName (int);
char const *cg_SimulationTypeName (int);
char const *cg_WallFunctionTypeName (int);
char const *cg_AreaTypeName (int);
char const *cg_AverageInterfaceTypeName (int);
#endif

int cg_get_identifier (const char *name, int *nexps, float *exps);
int cg_find_identifier (const char *pattern, int maxnames, char **names);
int cg_enum_identifier (int (*callback)(char *name,
    int nexps, float *exps, void *user), void *user);

