inline void CSolver::SetIterLinSolver(unsigned short val_iterlinsolver) { IterLinSolver = val_iterlinsolver; }

inline void CSolver::SetResLinSolver(su2double val_reslinsolver) { ResLinSolver = val_reslinsolver; }

inline void CSolver::SetNondimensionalization(CConfig *config, unsigned short iMesh) { }

inline bool CSolver::GetAdjoint(void) { return adjoint; }

inline unsigned short CSolver::GetIterLinSolver(void) { return IterLinSolver; }

inline su2double CSolver::GetCSensitivity(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
                                     unsigned short iMesh, unsigned short RunTime_EqSystem) { }

inline void CSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) { }

inline void CSolver::ResetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) { }

inline void CSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) { }

inline void CSolver::LoadRestart_FSI(CGeometry *geometry, CConfig *config, int val_iter) { }

inline void CSolver::PredictStruct_Displacement(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) { }

inline void CSolver::ComputeAitken_Coefficient(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution, unsigned long iOuterIter) { }

inline void CSolver::SetAitken_Relaxation(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) { }

inline void CSolver::Update_StructSolution(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) { }

inline void CSolver::Compute_OFRefGeom(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Compute_OFRefNode(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Compute_OFVolFrac(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Compute_OFCompliance(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::SetForceCoeff(su2double val_forcecoeff_history) { }

inline void CSolver::SetFSI_Residual(su2double val_FSI_residual) { }

inline void CSolver::SetRelaxCoeff(su2double val_relaxecoeff_history) { }

inline su2double CSolver::GetRelaxCoeff(void) { return 0.0; }

inline su2double CSolver::GetForceCoeff(void) { return 0.0; }

inline su2double CSolver::GetFSI_Residual(void) { return 0.0; }

inline void CSolver::Stiffness_Penalty(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container, CConfig *config) { }

inline void CSolver::SetCSensitivity(unsigned short val_marker, unsigned long val_vertex, su2double val_sensitivity) { }

inline void CSolver::Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline void CSolver::Smooth_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline void CSolver::Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline su2double CSolver::GetPhi_Inf(unsigned short val_dim) { return 0; }

inline su2double CSolver::GetPsiRho_Inf(void) { return 0; }

inline su2double* CSolver::GetPsiRhos_Inf(void) { return NULL; }

inline su2double CSolver::GetPsiE_Inf(void) { return 0; }

inline void CSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config, bool reconstruction) { }

inline void CSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config, bool reconstruction) { }

inline void CSolver::SetPrimitive_Limiter_MPI(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPreconditioner(CConfig *config, unsigned long iPoint) { }

inline void CSolver::SetDistance(CGeometry *geometry, CConfig *config) { };

inline su2double CSolver::GetCD_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCL_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CL(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CD(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CSF(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CEff(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFx(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFy(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFz(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMx(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMy(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMz(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CL_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CD_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CSF_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CL_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CD_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CSF_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CEff_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFx_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFy_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFz_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMx_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMy_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMz_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_Buffet_Metric(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CL_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CD_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CSF_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CEff_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFx_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFy_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFz_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMx_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMy_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMz_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetInflow_MassFlow(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetExhaust_MassFlow(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetInflow_Pressure(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetInflow_Mach(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCSF_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCEff_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_HF_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_MaxHF_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCL_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCSF_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCD_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetAllBound_CL_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CD_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CSF_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CEff_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CMx_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CMy_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CMz_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CoPx_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CoPy_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CoPz_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CFx_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CFy_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CFz_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CL_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CD_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CSF_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CEff_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CMx_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CMy_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CMz_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CoPx_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CoPy_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CoPz_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CFx_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CFy_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CFz_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CL_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CD_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CSF_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CEff_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CMx_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CMy_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CMz_Visc() { return 0; }

inline su2double CSolver::GetTotal_Buffet_Metric() { return 0; }

inline su2double CSolver::GetAllBound_CoPx_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CoPy_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CoPz_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CFx_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CFy_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CFz_Visc() { return 0; }

inline void CSolver::SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::SetIntBoundary_Jump(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline su2double CSolver::GetTotal_CL() { return 0; }

inline su2double CSolver::GetTotal_CD() { return 0; }

inline su2double CSolver::GetTotal_NetThrust() { return 0; }

inline su2double CSolver::GetTotal_Power() { return 0; }

inline su2double CSolver::GetTotal_SolidCD() { return 0; }

inline su2double CSolver::GetTotal_ReverseFlow() { return 0; }

inline su2double CSolver::GetTotal_MFR() { return 0; }

inline su2double CSolver::GetTotal_Prop_Eff() { return 0; }

inline su2double CSolver::GetTotal_ByPassProp_Eff() { return 0; }

inline su2double CSolver::GetTotal_Adiab_Eff() { return 0; }

inline su2double CSolver::GetTotal_Poly_Eff() { return 0; }

inline su2double CSolver::GetTotal_IDC_Mach() { return 0; }

inline su2double CSolver::GetTotal_DC60() { return 0; }

inline su2double CSolver::GetTotal_Custom_ObjFunc() { return 0; }

inline su2double CSolver::GetTotal_CMx() { return 0; }

inline su2double CSolver::GetTotal_CMy() { return 0; }

inline su2double CSolver::GetTotal_CMz() { return 0; }

inline su2double CSolver::GetTotal_CoPx() { return 0; }

inline su2double CSolver::GetTotal_CoPy() { return 0; }

inline su2double CSolver::GetTotal_CoPz() { return 0; }

inline su2double CSolver::GetTotal_CFx() { return 0; }

inline su2double CSolver::GetTotal_CFy() { return 0; }

inline su2double CSolver::GetTotal_CFz() { return 0; }

inline su2double CSolver::GetTotal_CSF() { return 0; }

inline su2double CSolver::GetTotal_CEff() { return 0; }

inline su2double CSolver::GetTotal_CT() { return 0; }

inline void CSolver::SetTotal_CT(su2double val_Total_CT) { }

inline su2double CSolver::GetTotal_CQ() { return 0; }

inline su2double CSolver::GetTotal_HeatFlux() { return 0; }

inline su2double CSolver::GetTotal_AvgTemperature() { return 0; }

inline su2double CSolver::GetTotal_MaxHeatFlux() { return 0; }

inline su2double CSolver::Get_PressureDrag() { return 0; }

inline su2double CSolver::Get_ViscDrag() { return 0; }

inline void CSolver::SetTotal_CQ(su2double val_Total_CQ) { }

inline void CSolver::SetTotal_HeatFlux(su2double val_Total_Heat) { }

inline void CSolver::SetTotal_MaxHeatFlux(su2double val_Total_Heat) { }

inline su2double CSolver::GetTotal_CMerit() { return 0; }

inline su2double CSolver::GetTotal_CEquivArea() { return 0; }

inline su2double CSolver::GetTotal_AeroCD() { return 0; }

inline su2double CSolver::GetTotal_IDR() { return 0; }

inline su2double CSolver::GetTotal_IDC() { return 0; }

inline su2double CSolver::GetTotal_CpDiff() { return 0; }

inline su2double CSolver::GetTotal_HeatFluxDiff() { return 0; }

inline su2double CSolver::GetTotal_CFEA() { return 0; }

inline su2double CSolver::GetTotal_CNearFieldOF() { return 0; }

inline su2double CSolver::GetTotal_OFRefGeom() { return 0; }

inline su2double CSolver::GetTotal_OFRefNode() { return 0; }

inline su2double CSolver::GetTotal_OFVolFrac() { return 0; }

inline su2double CSolver::GetTotal_OFCompliance() { return 0; }

inline bool CSolver::IsElementBased(void){ return false; }

inline void CSolver::AddTotal_ComboObj(su2double val_obj) {}

inline void CSolver::SetTotal_CEquivArea(su2double val_cequivarea) { }

inline void CSolver::SetTotal_AeroCD(su2double val_aerocd) { }

inline void CSolver::SetTotal_CpDiff(su2double val_pressure) { }

inline void CSolver::SetTotal_HeatFluxDiff(su2double val_heat) { }

inline void CSolver::SetTotal_CFEA(su2double val_cfea) { }

inline void CSolver::SetTotal_OFRefGeom(su2double val_ofrefgeom) { }

inline void CSolver::SetTotal_OFRefNode(su2double val_ofrefnode) { }

inline su2double CSolver::GetWAitken_Dyn(void) { return 0; }

inline su2double CSolver::GetWAitken_Dyn_tn1(void) { return 0; }

inline void CSolver::SetWAitken_Dyn(su2double waitk) {  }

inline void CSolver::SetWAitken_Dyn_tn1(su2double waitk_tn1) {  }

inline void CSolver::SetLoad_Increment(su2double val_loadIncrement) {  }

inline su2double CSolver::GetLoad_Increment() { return 0; }

inline void CSolver::SetTotal_CNearFieldOF(su2double val_cnearfieldpress) { }

inline su2double CSolver::GetTotal_CWave() { return 0; }

inline su2double CSolver::GetTotal_CHeat() { return 0; }

inline void CSolver::SetTotal_CL(su2double val_Total_CL) { }

inline void CSolver::SetTotal_CD(su2double val_Total_CD) { }

inline void CSolver::SetTotal_NetThrust(su2double val_Total_NetThrust) { }

inline void CSolver::SetTotal_Power(su2double val_Total_Power) { }

inline void CSolver::SetTotal_SolidCD(su2double val_Total_SolidCD) { }

inline void CSolver::SetTotal_ReverseFlow(su2double val_Total_ReverseFlow) { }

inline void CSolver::SetTotal_MFR(su2double val_Total_MFR) { }

inline void CSolver::SetTotal_Prop_Eff(su2double val_Total_Prop_Eff) { }

inline void CSolver::SetTotal_ByPassProp_Eff(su2double val_Total_ByPassProp_Eff) { }

inline void CSolver::SetTotal_Adiab_Eff(su2double val_Total_Adiab_Eff) { }

inline void CSolver::SetTotal_Poly_Eff(su2double val_Total_Poly_Eff) { }

inline void CSolver::SetTotal_IDC(su2double val_Total_IDC) { }

inline void CSolver::SetTotal_IDC_Mach(su2double val_Total_IDC_Mach) { }

inline void CSolver::SetTotal_IDR(su2double val_Total_IDR) { }

inline void CSolver::SetTotal_DC60(su2double val_Total_DC60) { }

inline void CSolver::SetTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { }

inline void CSolver::AddTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { }

inline su2double CSolver::GetCPressure(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetCPressureTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_pressure) { }

inline void CSolver::SetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_heat) { }

inline su2double *CSolver::GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { }

inline su2double *CSolver::GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { }

inline void CSolver::SetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { }

inline su2double CSolver::GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var) { return 0; }

inline su2double *CSolver::GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var) { return 0; }

inline unsigned long CSolver::GetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex, unsigned long val_index) { }

inline su2double CSolver::GetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex, su2double val_deltap) { }

inline su2double CSolver::GetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex, su2double val_deltat) { }

inline su2double CSolver::GetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim) { return 0; }

inline void CSolver::SetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ttotal) { }

inline void CSolver::SetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ptotal) { }

inline void CSolver::SetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_flowdir) { }

inline void CSolver::SetInlet_TurbVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_turb_var) { }

inline void CSolver::SetUniformInlet(CConfig* config, unsigned short iMarker) {};

inline void CSolver::SetInletAtVertex(su2double *val_inlet, unsigned short iMarker, unsigned long iVertex) { };

inline su2double CSolver::GetInletAtVertex(su2double *val_inlet, unsigned long val_inlet_point, unsigned short val_kind_marker, string val_marker, CGeometry *geometry, CConfig *config) { return 0; }

inline void CSolver::UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config) { }

inline su2double CSolver::GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim) { return 0; }

inline su2double CSolver::GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetBuffetSensor(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetYPlus(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetStrainMag_Max(void) { return 0; }

inline su2double CSolver::GetOmega_Max(void) { return 0; }

inline void CSolver::SetStrainMag_Max(su2double val_strainmag_max) { }

inline void CSolver::SetOmega_Max(su2double val_omega_max) { }

inline void CSolver::Viscous_Residual(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CNumerics *numerics, CConfig
                                      *config, unsigned short iMesh,
                                      unsigned short iRKstep) { }

inline void CSolver::AddStiffMatrix(su2double ** StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3) { }

inline void CSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                          CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) { }

inline void CSolver::Source_Template(CGeometry *geometry, CSolver **solver_container,
                          CNumerics *numerics, CConfig *config, unsigned short iMesh) { }

inline su2double CSolver::GetTotal_Sens_Geo() { return 0; }

inline su2double CSolver::GetTotal_Sens_Mach() { return 0; }

inline su2double CSolver::GetTotal_Sens_AoA() { return 0; }

inline su2double CSolver::GetTotal_Sens_Press() { return 0; }

inline su2double CSolver::GetTotal_Sens_Temp() { return 0; }

inline su2double CSolver::GetTotal_Sens_BPress() { return 0; }

inline su2double CSolver::GetTotal_Sens_Density() { return 0; }

inline su2double CSolver::GetTotal_Sens_ModVel() { return 0; }

inline su2double CSolver::GetDensity_Inf(void) { return 0; }

inline su2double CSolver::GetDensity_Inf(unsigned short val_var) { return 0; }

inline su2double CSolver::GetModVelocity_Inf(void) { return 0; }

inline su2double CSolver::GetDensity_Energy_Inf(void) { return 0; }

inline su2double CSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return 0; }

inline su2double CSolver::GetDensity_Velocity_Inf(unsigned short val_dim, unsigned short val_var) { return 0; }

inline su2double CSolver::GetVelocity_Inf(unsigned short val_dim) { return 0; }

inline su2double* CSolver::GetVelocity_Inf(void) { return 0; }

inline su2double CSolver::GetPressure_Inf(void) { return 0; }

inline su2double CSolver::GetViscosity_Inf(void) { return 0; }

inline su2double CSolver::GetNuTilde_Inf(void) { return 0; }

inline su2double CSolver::GetTke_Inf(void) { return 0; }

inline su2double CSolver::GetOmega_Inf(void) { return 0; }

inline su2double CSolver::GetTotal_Sens_E(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetTotal_Sens_Nu(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetTotal_Sens_Rho(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetTotal_Sens_Rho_DL(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetTotal_Sens_EField(unsigned short iEField) { return 0.0; }

inline su2double CSolver::GetTotal_Sens_DVFEA(unsigned short iDVFEA) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_E(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_Nu(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_Rho(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_Rho_DL(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_EField(unsigned short iEField) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_DVFEA(unsigned short iDVFEA) { return 0.0; }

inline su2double CSolver::GetVal_Young(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetVal_Poisson(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetVal_Rho(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetVal_Rho_DL(unsigned short iVal) { return 0.0; }

inline unsigned short CSolver::GetnEField(void) { return 0; }

inline unsigned short CSolver::GetnDVFEA(void) { return 0; }

inline void CSolver::ReadDV(CConfig *config) { }

inline su2double CSolver::GetVal_EField(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetVal_DVFEA(unsigned short iVal) { return 0.0; }

inline su2double* CSolver::GetConstants() { return NULL;}

inline void CSolver::SetTotal_ComboObj(su2double ComboObj) {}

inline su2double CSolver::GetTotal_ComboObj(void) { return 0;}

inline void CSolver::Set_Heatflux_Areas(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Evaluate_ObjFunc(CConfig *config) {};

inline void CSolver::Solve_System(CGeometry *geometry, CConfig *config) { }

inline void CSolver::BC_Euler_Wall(CGeometry      *geometry,
                                   CSolver        **solver_container,
                                   CNumerics      *conv_numerics,
                                   CNumerics      *visc_numerics,
                                   CConfig        *config,
                                   unsigned short val_marker) { }

inline void CSolver::BC_Clamped(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_DispDir(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Clamped_Post(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Normal_Displacement(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Normal_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Dir_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Sine_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Damper(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Deforming(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                  unsigned short val_marker) { }

inline void CSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                         CConfig *config) { }

inline void CSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                           CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                           CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config) { }

inline void CSolver::BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                      CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                       CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker, bool val_inlet_surface) { }

inline void CSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                    CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                  CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_TurboRiemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
										 CConfig *config, unsigned short val_marker) { }

inline void CSolver::PreprocessBC_Giles(CGeometry *geometry, CConfig *config,
																								CNumerics *conv_numerics,unsigned short marker_flag){}

inline void CSolver::BC_Giles(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                      CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                         CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                      CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                      CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                      CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Dielec(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                   CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Electrode(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                  CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                      CConfig *config, unsigned short val_marker) { }

inline void CSolver::GetPower_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) { }

inline void CSolver::GetOutlet_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) { }

inline void CSolver::GetEllipticSpanLoad_Diff(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                     CConfig *config, unsigned short iMesh, bool Output) { }

inline bool CSolver::FixedCL_Convergence(CConfig *config, bool convergence) { return false; }

inline bool CSolver::GetStart_AoA_FD(void) { return false; }

inline bool CSolver::GetEnd_AoA_FD(void) { return false; }

inline unsigned long CSolver::GetIter_Update_AoA(void) { return 0; }

inline su2double CSolver::GetPrevious_AoA(void) { return 0.0; }

inline su2double CSolver::GetAoA_inc(void) { return 0.0; }

inline void CSolver::SetActDisk_BCThrust(CGeometry *geometry, CSolver **solver_container,
                                         CConfig *config, unsigned short iMesh, bool Output) { }

inline void CSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
							        unsigned short iMesh, unsigned long Iteration) { }

inline void CSolver::CheckTimeSynchronization(CConfig         *config,
                                              const su2double TimeSync,
                                              su2double       &timeEvolved,
                                              bool            &syncTimeReached) {}

inline void CSolver::ProcessTaskList_DG(CGeometry *geometry,  CSolver **solver_container,
                                        CNumerics **numerics, CConfig *config,
                                        unsigned short iMesh) {}

inline void CSolver::ADER_SpaceTimeIntegration(CGeometry *geometry,  CSolver **solver_container,
                                               CNumerics **numerics, CConfig *config,
                                               unsigned short iMesh, unsigned short RunTime_EqSystem) {}

inline void CSolver::ComputeSpatialJacobian(CGeometry *geometry,  CSolver **solver_container,
                                            CNumerics **numerics, CConfig *config,
                                            unsigned short iMesh, unsigned short RunTime_EqSystem) {}

inline void CSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                      unsigned short iMesh) { }

inline void CSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
                      unsigned short iMesh) { }

inline void CSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                    CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

inline void CSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                     CConfig *config, unsigned short iMesh) { }

inline void CSolver::Convective_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                       CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

inline void CSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) { }

inline void CSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) { }

inline void CSolver::SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetUpwind_Ducros_Sensor(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Nearfield(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Pressure_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Momentum_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Friction_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Buffet_Monitoring(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Heat_Fluxes(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Inviscid_DeltaForces(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Viscous_DeltaForces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Wave_Strength(CGeometry *geometry, CConfig *config) { }

inline void CSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                      CConfig *config, unsigned short iRKStep) { }

inline void CSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                            CConfig *config, unsigned short iRKStep) { }

inline void CSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::ComputeUnderRelaxationFactor(CSolver **solver_container, CConfig *config) { }

inline void CSolver::ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::ImplicitNewmark_Update(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::ImplicitNewmark_Relaxation(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::GeneralizedAlpha_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::GeneralizedAlpha_UpdateDisp(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::GeneralizedAlpha_UpdateSolution(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::GeneralizedAlpha_UpdateLoads(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Compute_Residual(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                    unsigned short iMesh) { }

inline void CSolver::SetRes_RMS(unsigned short val_var, su2double val_residual) { Residual_RMS[val_var] = val_residual; }

inline void CSolver::AddRes_RMS(unsigned short val_var, su2double val_residual) { Residual_RMS[val_var] += val_residual; }

inline su2double CSolver::GetRes_RMS(unsigned short val_var) { return Residual_RMS[val_var]; }

inline void CSolver::SetRes_Max(unsigned short val_var, su2double val_residual, unsigned long val_point) { Residual_Max[val_var] = val_residual; Point_Max[val_var] = val_point; }

inline void CSolver::AddRes_Max(unsigned short val_var, su2double val_residual, unsigned long val_point, su2double* val_coord) {
  if (val_residual > Residual_Max[val_var]) {
  Residual_Max[val_var] = val_residual;
  Point_Max[val_var] = val_point;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Point_Max_Coord[val_var][iDim] = val_coord[iDim];
  }
}

inline void CSolver::AddRes_Max(unsigned short val_var, su2double val_residual, unsigned long val_point, const su2double* val_coord) {
  if (val_residual > Residual_Max[val_var]) {
    Residual_Max[val_var] = val_residual;
    Point_Max[val_var] = val_point;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Point_Max_Coord[val_var][iDim] = val_coord[iDim];
  }
}

inline su2double CSolver::GetRes_Max(unsigned short val_var) { return Residual_Max[val_var]; }

inline void CSolver::SetRes_BGS(unsigned short val_var, su2double val_residual) { Residual_BGS[val_var] = val_residual; }

inline void CSolver::AddRes_BGS(unsigned short val_var, su2double val_residual) { Residual_BGS[val_var] += val_residual; }

inline su2double CSolver::GetRes_BGS(unsigned short val_var) { return Residual_BGS[val_var]; }

inline void CSolver::SetRes_Max_BGS(unsigned short val_var, su2double val_residual, unsigned long val_point) { Residual_Max_BGS[val_var] = val_residual; Point_Max_BGS[val_var] = val_point; }

inline void CSolver::AddRes_Max_BGS(unsigned short val_var, su2double val_residual, unsigned long val_point, su2double* val_coord) {
  if (val_residual > Residual_Max_BGS[val_var]) {
  Residual_Max_BGS[val_var] = val_residual;
  Point_Max_BGS[val_var] = val_point;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Point_Max_Coord_BGS[val_var][iDim] = val_coord[iDim];
  }
}

inline su2double CSolver::GetRes_Max_BGS(unsigned short val_var) { return Residual_Max_BGS[val_var]; }

inline su2double CSolver::GetRes_FEM(unsigned short val_var) { return 0.0; }

inline unsigned long CSolver::GetPoint_Max(unsigned short val_var) { return Point_Max[val_var]; }

inline su2double* CSolver::GetPoint_Max_Coord(unsigned short val_var) { return Point_Max_Coord[val_var]; }

inline unsigned long CSolver::GetPoint_Max_BGS(unsigned short val_var) { return Point_Max_BGS[val_var]; }

inline su2double* CSolver::GetPoint_Max_Coord_BGS(unsigned short val_var) { return Point_Max_Coord_BGS[val_var]; }

inline void CSolver::Set_OldSolution(CGeometry *geometry) { base_nodes->Set_OldSolution(); }

inline void CSolver::Set_NewSolution(CGeometry *geometry) { }

inline unsigned short CSolver::GetnVar(void) { return nVar; }

inline unsigned short CSolver::GetnOutputVariables(void) { return nOutputVariables; }

inline unsigned short CSolver::GetnPrimVar(void) { return nPrimVar; }

inline unsigned short CSolver::GetnPrimVarGrad(void) { return nPrimVarGrad; }

inline unsigned short CSolver::GetnSecondaryVar(void) { return nSecondaryVar; }

inline unsigned short CSolver::GetnSecondaryVarGrad(void) { return nSecondaryVarGrad; }

inline su2double CSolver::GetMax_Delta_Time(void) { return Max_Delta_Time; }

inline su2double CSolver::GetMin_Delta_Time(void) { return Min_Delta_Time; }

inline su2double CSolver::GetMax_Delta_Time(unsigned short val_Species) { return 0.0; }

inline su2double CSolver::GetMin_Delta_Time(unsigned short val_Species) { return 0.0; }

inline void CSolver::Copy_Zone_Solution(CSolver ***solver1_solution, CGeometry **solver1_geometry, CConfig *solver1_config,
                      CSolver ***solver2_solution, CGeometry **solver2_geometry, CConfig *solver2_config) {};

inline CFluidModel* CSolver::GetFluidModel(void) { return NULL;}

inline su2double* CSolver::GetVecSolDOFs(void) {return NULL;}

inline unsigned long CSolver::GetnDOFsGlobal(void) {return 0;}

inline void CSolver::Set_ReferenceGeometry(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_Prestretch(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_ElementProperties(CGeometry *geometry, CConfig *config) { }

inline su2double CSolver::Compute_LoadCoefficient(su2double CurrentTime, su2double RampTime, CConfig *config) { return 0.0; }

inline su2double CSolver::Get_ValCoord(CGeometry *geometry, unsigned long indexNode, unsigned short iDim) {return 0.0;}

inline void CSolver::Compute_StiffMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_MassMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_MassRes(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_NodalStressRes(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_NodalStress(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_DeadLoad(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Compute_IntegrationConstants(CConfig *config) { }

inline void CSolver::SetFSI_ConvValue(unsigned short val_index, su2double val_criteria) { };

inline su2double CSolver::GetFSI_ConvValue(unsigned short val_index) { return 0.0; }

inline void CSolver::RegisterSolution(CGeometry *geometry_container, CConfig *config){}

inline void CSolver::RegisterOutput(CGeometry *geometry_container, CConfig *config){}

inline void CSolver::SetAdjoint_Output(CGeometry *geometry, CConfig *config){}

inline void CSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){}

inline void CSolver::RegisterObj_Func(CConfig *config){}

inline void CSolver::SetSurface_Sensitivity(CGeometry *geometry, CConfig *config){}

inline void CSolver::SetSensitivity(CGeometry *geometry, CSolver **solver, CConfig *config){}

inline void CSolver::SetAdj_ObjFunc(CGeometry *geometry, CConfig *config){}

inline unsigned long CSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {return 0;}

inline void CSolver::SetRecording(CGeometry *geometry, CConfig *config){}

inline void CSolver::SetPressure_Inf(su2double p_inf){}

inline void CSolver::SetTemperature_Inf(su2double t_inf){}

inline void CSolver::SetDensity_Inf(su2double rho_inf){}

inline void CSolver::SetVelocity_Inf(unsigned short val_dim, su2double val_velocity) { }

inline void CSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset){}

inline void CSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config){}

inline void CSolver::SetFreeStream_Solution(CConfig *config){}

inline void CSolver::SetTauWall_WF(CGeometry *geometry, CSolver** solver_container, CConfig* config){}

inline void CSolver::SetNuTilde_WF(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                           CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {}

inline void CSolver::InitTurboContainers(CGeometry *geometry, CConfig *config){}

inline void CSolver::PreprocessAverage(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag){}

inline void CSolver::TurboAverageProcess(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag){}

inline void CSolver::GatherInOutAverageValues(CConfig *config, CGeometry *geometry){ }

inline su2double CSolver::GetAverageDensity(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetAveragePressure(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double* CSolver::GetAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan){return NULL;}

inline su2double CSolver::GetAverageNu(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetAverageKine(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetAverageOmega(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetExtAverageNu(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetExtAverageKine(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetExtAverageOmega(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline void CSolver::SetExtAverageDensity(unsigned short valMarker, unsigned short valSpan, su2double valDensity){ }

inline void CSolver::SetExtAveragePressure(unsigned short valMarker, unsigned short valSpan, su2double valPressure){ }

inline void CSolver::SetExtAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan, unsigned short valIndex, su2double valTurboVelocity){ }

inline void CSolver::SetExtAverageNu(unsigned short valMarker, unsigned short valSpan, su2double valNu){ }

inline void CSolver::SetExtAverageKine(unsigned short valMarker, unsigned short valSpan, su2double valKine){ }

inline void CSolver::SetExtAverageOmega(unsigned short valMarker, unsigned short valSpan, su2double valOmega){ }

inline su2double CSolver::GetDensityIn(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetPressureIn(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double* CSolver::GetTurboVelocityIn(unsigned short inMarkerTP, unsigned short valSpan){return NULL;}

inline su2double CSolver::GetDensityOut(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetPressureOut(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double* CSolver::GetTurboVelocityOut(unsigned short inMarkerTP, unsigned short valSpan){return NULL;}

inline su2double CSolver::GetKineIn(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetOmegaIn(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetNuIn(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetKineOut(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetOmegaOut(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetNuOut(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline void CSolver::SetDensityIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetPressureIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetTurboVelocityIn(su2double *value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetDensityOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetPressureOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetTurboVelocityOut(su2double *value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetKineIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetOmegaIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetNuIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetKineOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetOmegaOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetNuOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetFreeStream_TurboSolution(CConfig *config){ }

inline void CSolver::SetBeta_Parameter(CGeometry *geometry, CSolver **solver_container,
                                     CConfig *config, unsigned short iMesh) { }

inline void CSolver::SetRoe_Dissipation(CGeometry *geometry, CConfig *config) {}

inline void CSolver::SetDES_LengthScale(CSolver** solver, CGeometry *geometry, CConfig *config) { }

inline void CSolver::DeformMesh(CGeometry **geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::SetMesh_Stiffness(CGeometry **geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var, su2double relaxation_factor, su2double val_var) { }

inline su2double CSolver::GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var) { return 0.0; }

inline void CSolver::ComputeVerificationError(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetImplicitPeriodic(bool val_implicit_periodic) { implicit_periodic = val_implicit_periodic; }

inline void CSolver::SetRotatePeriodic(bool val_rotate_periodic) { rotate_periodic = val_rotate_periodic; }

inline string CSolver::GetSolverName(void) {return SolverName;}

inline void CSolver::RefGeom_Sensitivity(CGeometry *geometry, CSolver **solver_container, CConfig *config){ }

inline void CSolver::DE_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){ }

inline void CSolver::Stiffness_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){ }

inline unsigned short CSolver::Get_iElem_iDe(unsigned long iElem){ return 0; }

inline void CSolver::Set_DV_Val(su2double val_EField, unsigned short i_DV){ }

inline su2double CSolver::Get_DV_Val(unsigned short i_DV){ return 0.0; }

inline su2double CSolver::Get_val_I(void){ return 0.0; }

inline su2double CSolver::Get_MassMatrix(unsigned long iPoint, unsigned long jPoint, unsigned short iVar, unsigned short jVar){ return 0.0; }

inline void CSolver::SetAdjoint_OutputMesh(CGeometry *geometry, CConfig *config) {}

inline void CSolver::ExtractAdjoint_Geometry(CGeometry *geometry, CConfig *config) {}

inline void CSolver::ExtractAdjoint_CrossTerm(CGeometry *geometry, CConfig *config) {}

inline void CSolver::ExtractAdjoint_CrossTerm_Geometry(CGeometry *geometry, CConfig *config) {}

inline void CSolver::ExtractAdjoint_CrossTerm_Geometry_Flow(CGeometry *geometry, CConfig *config) {}

inline void CSolver::SetMesh_Recording(CGeometry **geometry, CVolumetricMovement *grid_movement, CConfig *config) {}

inline void CSolver::SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component){ }

inline su2double CSolver::GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index) { return 0; }

inline int CSolver::GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex){ return 0; }

inline void CSolver::SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value){}

inline void CSolver::SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex){}

inline void CSolver::SetDualTime_Mesh(void){ }

inline vector<string> CSolver::GetSolutionFields(){return fields;}

