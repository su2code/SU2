/*!
 * \file variable_structure.inl
 * \brief In-Line subroutines of the <i>variable_structure.hpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

inline bool CVariable::SetDensity(void) { return 0; }

inline bool CVariable::SetDensity(su2double val_density){ return 0; }

inline void CVariable::SetVelSolutionOldDVector(void) { }

inline void CVariable::SetVelSolutionDVector(void) { }

inline void CVariable::SetTraction(unsigned short iVar, unsigned short jVar, su2double val_traction) { }

inline void CVariable::AddTraction(unsigned short iVar, unsigned short jVar, su2double val_traction) { }

inline su2double **CVariable::GetTraction(void) { return NULL; }

inline void CVariable::SetStress_FEM(unsigned short iVar, su2double val_stress) { }

inline void CVariable::AddStress_FEM(unsigned short iVar, su2double val_stress) { }

inline su2double *CVariable::GetStress_FEM(void) { return NULL; }

inline void CVariable::SetVonMises_Stress(su2double val_stress) { }

inline su2double CVariable::GetVonMises_Stress(void) { return 0; }

inline void CVariable::Add_SurfaceLoad_Res(su2double *val_surfForce) { }

inline void CVariable::Set_SurfaceLoad_Res(unsigned short iVar, su2double val_surfForce) { }

inline su2double *CVariable::Get_SurfaceLoad_Res(void) { return NULL;}

inline su2double CVariable::Get_SurfaceLoad_Res(unsigned short iVar) { return 0.0;}

inline void CVariable::Clear_SurfaceLoad_Res(void) { }

inline void CVariable::Set_SurfaceLoad_Res_n(void) { }

inline su2double CVariable::Get_SurfaceLoad_Res_n(unsigned short iVar) { return 0.0;}

inline void CVariable::Add_BodyForces_Res(su2double *val_bodyForce) { }

inline su2double *CVariable::Get_BodyForces_Res(void) { return NULL;}

inline su2double CVariable::Get_BodyForces_Res(unsigned short iVar) { return 0.0;}

inline void CVariable::Clear_BodyForces_Res(void) { }

inline void CVariable::Set_FlowTraction(su2double *val_flowTraction) { }

inline void CVariable::Add_FlowTraction(su2double *val_flowTraction) { }

inline su2double *CVariable::Get_FlowTraction(void) { return NULL;}

inline su2double CVariable::Get_FlowTraction(unsigned short iVar) { return 0.0;}

inline void CVariable::Clear_FlowTraction(void) { }

inline void CVariable::Set_FlowTraction_n(void) { }

inline su2double CVariable::Get_FlowTraction_n(unsigned short iVar) { return 0.0; }

inline su2double CVariable::GetBetaInc2(void) { return 0; }

inline su2double CVariable::GetMassFraction(unsigned short val_Species) { return 0; }

inline void CVariable::SetNon_Physical(bool val_value) { Non_Physical = !val_value; }

inline su2double CVariable::GetNon_Physical(void) { return su2double(Non_Physical); }

inline void CVariable::SetSolution(unsigned short val_var, su2double val_solution) { Solution[val_var] = val_solution; }

inline void CVariable::Add_DeltaSolution(unsigned short val_var, su2double val_solution) { Solution[val_var] += val_solution; }

inline void CVariable::SetUndivided_Laplacian(unsigned short val_var, su2double val_undivided_laplacian) { Undivided_Laplacian[val_var] = val_undivided_laplacian; }

inline void CVariable::SetAuxVar(su2double val_auxvar) { AuxVar = val_auxvar; }

inline void CVariable::SetSolution_Old(unsigned short val_var, su2double val_solution_old) { Solution_Old[val_var] = val_solution_old; }

inline void CVariable::SetLimiter(unsigned short val_var, su2double val_limiter) { Limiter[val_var] = val_limiter; }

inline void CVariable::SetLimiterPrimitive(unsigned short val_species, unsigned short val_var, su2double val_limiter) { }

inline su2double CVariable::GetLimiterPrimitive(unsigned short val_species, unsigned short val_var) { return 0.0; }

inline void CVariable::SetSolution_Max(unsigned short val_var, su2double val_solution) { Solution_Max[val_var] = val_solution; }

inline void CVariable::SetSolution_Min(unsigned short val_var, su2double val_solution) { Solution_Min[val_var] = val_solution; }

inline void CVariable::SetAuxVarGradient(unsigned short iDim, su2double val_gradient) { Grad_AuxVar[iDim] = val_gradient; }

inline su2double *CVariable::GetSolution(void) { return Solution; }

inline su2double *CVariable::GetSolution_Old(void) { return Solution_Old; }

inline su2double *CVariable::GetSolution_time_n(void) { return Solution_time_n; }

inline su2double *CVariable::GetSolution_time_n1(void) { return Solution_time_n1; }

inline su2double CVariable::GetAuxVar(void) { return AuxVar; }

inline su2double *CVariable::GetUndivided_Laplacian(void) { return Undivided_Laplacian; }

inline su2double CVariable::GetUndivided_Laplacian(unsigned short val_var) { return Undivided_Laplacian[val_var]; }

inline su2double CVariable::GetSolution(unsigned short val_var) { return Solution[val_var]; }

inline su2double CVariable::GetSolution_Old(unsigned short val_var) { return Solution_Old[val_var]; }

inline su2double CVariable::GetSolution_Old_Adj(unsigned short val_var) { return Solution_Adj_Old[val_var]; }

inline su2double *CVariable::GetResidual_Sum(void) { return Residual_Sum; }

inline su2double *CVariable::GetResidual_Old(void) { return Residual_Old; }

inline void CVariable::SetGradient(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient[val_var][val_dim] = val_value; }

inline void CVariable::AddGradient(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient[val_var][val_dim] += val_value; }

inline void CVariable::SubtractGradient(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient[val_var][val_dim] -= val_value; }

inline void CVariable::AddAuxVarGradient(unsigned short val_dim, su2double val_value) { Grad_AuxVar[val_dim] += val_value; }

inline void CVariable::SubtractAuxVarGradient(unsigned short val_dim, su2double val_value) { Grad_AuxVar[val_dim] -= val_value; }

inline su2double CVariable::GetGradient(unsigned short val_var, unsigned short val_dim) { return Gradient[val_var][val_dim]; }

inline su2double CVariable::GetLimiter(unsigned short val_var) { return Limiter[val_var]; }

inline su2double CVariable::GetSolution_Max(unsigned short val_var) { return Solution_Max[val_var]; }

inline su2double CVariable::GetSolution_Min(unsigned short val_var) { return Solution_Min[val_var]; }

inline su2double CVariable::GetPreconditioner_Beta() { return 0; }

inline void CVariable::SetPreconditioner_Beta( su2double val_Beta) { }

inline su2double* CVariable::GetWindGust() { return 0; }

inline void CVariable::SetWindGust( su2double* val_WindGust) {}

inline su2double* CVariable::GetWindGustDer() { return 0; }

inline void CVariable::SetWindGustDer( su2double* val_WindGustDer) {}

inline su2double **CVariable::GetGradient(void) { return Gradient; }

inline su2double *CVariable::GetLimiter(void) { return Limiter; }

inline su2double *CVariable::GetAuxVarGradient(void) { return Grad_AuxVar; }

inline su2double CVariable::GetAuxVarGradient(unsigned short val_dim) { return Grad_AuxVar[val_dim]; }

inline su2double *CVariable::GetResTruncError(void) { return Res_TruncError; }

inline void CVariable::SetDelta_Time(su2double val_delta_time) { Delta_Time = val_delta_time; }

inline void CVariable::SetDelta_Time(su2double val_delta_time, unsigned short iSpecies) {  }

inline su2double CVariable::GetDelta_Time(void) { return Delta_Time; }

inline su2double CVariable::GetDelta_Time(unsigned short iSpecies) { return 0;}

inline void CVariable::SetMax_Lambda(su2double val_max_lambda) { Max_Lambda = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Inv(su2double val_max_lambda) { Max_Lambda_Inv = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Inv(su2double val_max_lambda, unsigned short val_species) { }

inline void CVariable::SetMax_Lambda_Visc(su2double val_max_lambda) { Max_Lambda_Visc = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Visc(su2double val_max_lambda, unsigned short val_species) { }

inline void CVariable::SetLambda(su2double val_lambda) { Lambda = val_lambda; }

inline void CVariable::SetLambda(su2double val_lambda, unsigned short iSpecies) {}

inline void CVariable::AddMax_Lambda(su2double val_max_lambda) { Max_Lambda += val_max_lambda; }

inline void CVariable::AddMax_Lambda_Inv(su2double val_max_lambda) { Max_Lambda_Inv += val_max_lambda; }

inline void CVariable::AddMax_Lambda_Visc(su2double val_max_lambda) { Max_Lambda_Visc += val_max_lambda; }

inline void CVariable::AddLambda(su2double val_lambda) { Lambda += val_lambda; }

inline void CVariable::AddLambda(su2double val_lambda, unsigned short iSpecies) {}

inline su2double CVariable::GetMax_Lambda(void) { return Max_Lambda; }

inline su2double CVariable::GetMax_Lambda_Inv(void) { return Max_Lambda_Inv; }

inline su2double CVariable::GetMax_Lambda_Visc(void) { return Max_Lambda_Visc; }

inline su2double CVariable::GetLambda(void) { return Lambda; }

inline su2double CVariable::GetLambda(unsigned short iSpecies) { return 0; }

inline su2double CVariable::GetSensor(void) { return Sensor; }

inline su2double CVariable::GetSensor(unsigned short iSpecies) { return 0;}

inline void CVariable::AddMax_Lambda_Inv(su2double val_max_lambda, unsigned short iSpecies) { }

inline void CVariable::AddMax_Lambda_Visc(su2double val_max_lambda, unsigned short iSpecies) { }

inline void CVariable::SetSensor(su2double val_sensor) { Sensor = val_sensor; }

inline void CVariable::SetSensor(su2double val_sensor, unsigned short val_iSpecies) {}

inline su2double CVariable::GetDensity(void) {  return 0; }

inline su2double CVariable::GetDensity_Old(void) {  return 0; }

inline su2double CVariable::GetDensity(unsigned short val_iSpecies) {  return 0; }

inline su2double CVariable::GetEnergy(void) { return 0; }

inline su2double *CVariable::GetForceProj_Vector(void) { return NULL; }

inline su2double *CVariable::GetObjFuncSource(void) { return NULL; }

inline su2double *CVariable::GetIntBoundary_Jump(void) { return NULL; }

inline su2double CVariable::GetEddyViscosity(void) { return 0; }

inline void CVariable::SetGammaEff(void) { }

inline void CVariable::SetGammaSep(su2double gamma_sep) { }

inline su2double CVariable::GetIntermittency(void) { return 0; }

inline su2double CVariable::GetEnthalpy(void) { return 0; }

inline su2double CVariable::GetPressure(void) { return 0; }

inline su2double CVariable::GetProjVel(su2double *val_vector) { return 0; }

inline su2double CVariable::GetProjVel(su2double *val_vector, unsigned short val_species) { return 0; }

inline su2double CVariable::GetSoundSpeed(void) { return 0; }

inline su2double CVariable::GetTemperature(void) { return 0; }

inline su2double CVariable::GetTemperature_ve(void) { return 0; }

inline su2double CVariable::GetRhoCv_tr(void) { return 0; }

inline su2double CVariable::GetRhoCv_ve(void) { return 0; }

inline su2double CVariable::GetVelocity(unsigned short val_dim) { return 0; }

inline su2double CVariable::GetVelocity2(void) { return 0; }

inline su2double CVariable::GetVelocity2(unsigned short val_species) { return 0;}

inline su2double CVariable::GetLaminarViscosity(void) { return 0; }

inline su2double CVariable::GetLaminarViscosity(unsigned short iSpecies) { return 0; }

inline su2double* CVariable::GetDiffusionCoeff(void) { return NULL; }

inline su2double CVariable::GetThermalConductivity(void) { return 0; }

inline su2double CVariable::GetSpecificHeatCp(void) { return 0; }

inline su2double CVariable::GetSpecificHeatCv(void) { return 0; }

inline su2double CVariable::GetThermalConductivity_ve(void) { return 0; }

inline su2double* CVariable::GetVorticity(void) { return 0; }

inline su2double CVariable::GetStrainMag(void) { return 0; }

inline void CVariable::SetForceProj_Vector(su2double *val_ForceProj_Vector) { }

inline void CVariable::SetObjFuncSource(su2double *val_ObjFuncSource) { }

inline void CVariable::SetIntBoundary_Jump(su2double *val_IntBoundary_Jump) { }

inline su2double CVariable::GetGammaBC(void) { return 0; }

inline void CVariable::SetGammaBC(su2double val_gamma) { }

inline void CVariable::SetEnthalpy(void) { }

inline bool CVariable::SetPrimVar(su2double SharpEdge_Distance, bool check, CConfig *config) { return true; }

inline bool CVariable::SetPrimVar(CConfig *config) { return true; }

inline bool CVariable::SetPrimVar(CFluidModel *FluidModel) { return true; }

inline void CVariable::SetSecondaryVar(CFluidModel *FluidModel) { }

inline bool CVariable::SetPrimVar(su2double eddy_visc, su2double turb_ke, CConfig *config) { return true; }

inline bool CVariable::SetPrimVar(su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) { return true; }

inline bool CVariable::SetPrimVar(su2double Density_Inf, CConfig *config) { return true; }

inline bool CVariable::SetPrimVar(su2double Density_Inf, su2double Viscosity_Inf, su2double eddy_visc, su2double turb_ke, CConfig *config) { return true; }

inline su2double CVariable::GetPrimitive(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetPrimitive(void) { return NULL; }

inline void CVariable::SetPrimitive(unsigned short val_var, su2double val_prim) { }

inline void CVariable::SetPrimitive(su2double *val_prim) { }

inline su2double CVariable::GetSecondary(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetSecondary(void) { return NULL; }

inline void CVariable::SetSecondary(unsigned short val_var, su2double val_secondary) { }

inline void CVariable::SetSecondary(su2double *val_prim) { }

inline bool CVariable::Cons2PrimVar(CConfig *config, su2double *U, su2double *V,
                                    su2double *val_dPdU, su2double *val_dTdU,
                                    su2double *val_dTvedU) { return false; }

inline void CVariable::Prim2ConsVar(CConfig *config, su2double *V, su2double *U) { return; }

inline void CVariable::SetBetaInc2(su2double val_betainc2) { }

inline void CVariable::SetPhi_Old(su2double *val_phi) { }

inline void CVariable::SetdPdrho_e(su2double dPdrho_e) { }

inline void CVariable::SetdPde_rho(su2double dPde_rho) { }

inline void CVariable::SetdTdrho_e(su2double dTdrho_e) { }

inline void CVariable::SetdTde_rho(su2double dTde_rho) { }

inline void CVariable::Setdmudrho_T(su2double dmudrho_T) { }

inline void CVariable::SetdmudT_rho(su2double dmudT_rho) { }

inline void CVariable::Setdktdrho_T(su2double dktdrho_T) { }

inline void CVariable::SetdktdT_rho(su2double dktdT_rho) { }

inline bool CVariable::SetPressure(su2double Gamma) { return false; }

inline bool CVariable::SetPressure(CConfig *config) { return false; }

inline bool CVariable::SetPressure(su2double Gamma, su2double turb_ke) { return false; }

inline void CVariable::SetPressure() { }

inline su2double *CVariable::GetdPdU() { return NULL; }

inline su2double *CVariable::GetdTdU() { return NULL; }

inline su2double *CVariable::GetdTvedU() { return NULL; }

inline su2double CVariable::CalcEve(su2double *V, CConfig *config, unsigned short val_Species) { return 0; }

inline su2double CVariable::CalcHs(su2double *V, CConfig *config, unsigned short val_Species) { return 0; }

inline su2double CVariable::CalcCvve(su2double val_Tve, CConfig *config, unsigned short val_Species) { return 0; }

inline void CVariable::CalcdPdU(su2double *V, CConfig *config, su2double *dPdU) { }

inline void CVariable::CalcdTdU(su2double *V, CConfig *config, su2double *dTdU) { }

inline void CVariable::CalcdTvedU(su2double *V, CConfig *config, su2double *dTvedU) { }

inline void CVariable::SetDeltaPressure(su2double *val_velocity, su2double Gamma) { }

inline bool CVariable::SetSoundSpeed(CConfig *config) { return false; }

inline bool CVariable::SetSoundSpeed() { return false; }

inline bool CVariable::SetSoundSpeed(su2double Gamma) { return false; }

inline bool CVariable::SetTemperature(su2double Gas_Constant) { return false; }

inline bool CVariable::SetTemperature_ve(su2double val_Tve) { return false; }

inline bool CVariable::SetTemperature(CConfig *config) { return false; }

inline void CVariable::SetPrimitive(CConfig *config) { }

inline void CVariable::SetPrimitive(CConfig *config, su2double *Coord) { }

inline void CVariable::SetWallTemperature(su2double Temperature_Wall) { }

inline void CVariable::SetWallTemperature(su2double* Temperature_Wall) { }

inline void CVariable::SetThermalCoeff(CConfig *config) { }

inline void CVariable::SetVelocity(void) { }

inline void CVariable::SetVelocity2(void) { }

inline void CVariable::SetVelocity_Old(su2double *val_velocity) { }

inline void CVariable::SetVel_ResTruncError_Zero(unsigned short iSpecies) { }

inline void CVariable::SetLaminarViscosity(su2double laminarViscosity) { }

inline void CVariable::SetLaminarViscosity(CConfig *config) { }

inline void CVariable::SetEddyViscosity(su2double eddy_visc) { }

inline void CVariable::SetThermalConductivity(su2double thermalConductivity) { }

inline void CVariable::SetThermalConductivity(CConfig *config) { }

inline void CVariable::SetSpecificHeatCp(su2double Cp) { }

inline void CVariable::SetSpecificHeatCv(su2double Cv) { }

inline bool CVariable::SetVorticity(void) { return false; }

inline bool CVariable::SetStrainMag(void) { return false; }

inline void CVariable::SetTauWall(su2double val_tau_wall) { }

inline su2double CVariable::GetTauWall(void) { return 0; }

inline void CVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) { }

inline void CVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline void CVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline su2double CVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return 0; }

inline su2double CVariable::GetLimiter_Primitive(unsigned short val_var) { return 0; }

inline void CVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline void CVariable::SetLimiter_Primitive(unsigned short val_var, su2double val_value) { }

inline su2double **CVariable::GetGradient_Primitive(void) { return NULL; }

inline su2double *CVariable::GetLimiter_Primitive(void) { return NULL; }

inline void CVariable::SetGradient_SecondaryZero(unsigned short val_secondaryvar) { }

inline void CVariable::AddGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline void CVariable::SubtractGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline su2double CVariable::GetGradient_Secondary(unsigned short val_var, unsigned short val_dim) { return 0; }

inline su2double CVariable::GetLimiter_Secondary(unsigned short val_var) { return 0; }

inline void CVariable::SetGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline void CVariable::SetLimiter_Secondary(unsigned short val_var, su2double val_value) { }

inline su2double **CVariable::GetGradient_Secondary(void) { return NULL; }

inline su2double *CVariable::GetLimiter_Secondary(void) { return NULL; }

inline void CVariable::SetBlendingFunc(su2double val_viscosity, su2double val_dist, su2double val_density) { }

inline su2double CVariable::GetF1blending(void) { return 0; }

inline su2double CVariable::GetF2blending(void) { return 0; }

inline su2double CVariable::GetmuT() { return 0;}

inline void CVariable::SetmuT(su2double val_muT) { }

inline su2double* CVariable::GetSolution_Direct() { return NULL; }

inline void CVariable::SetSolution_Direct(su2double *val_solution_direct) { }

inline su2double* CVariable::GetGeometry_Direct() { return NULL; }

inline su2double CVariable::GetGeometry_Direct(unsigned short val_dim) { return 0.0; }

inline void CVariable::SetGeometry_Direct(su2double *val_geometry_direct) { }

inline su2double CVariable::GetSolution_Geometry(unsigned short val_var) { return 0.0;}

inline void CVariable::SetSolution_Geometry(su2double *val_solution_geometry) { }

inline void CVariable::SetSolution_Geometry(unsigned short val_var, su2double val_solution_geometry) { }

inline su2double CVariable::GetGeometry_CrossTerm_Derivative(unsigned short val_var) { return 0.0;}

inline void CVariable::SetGeometry_CrossTerm_Derivative(unsigned short iDim, su2double der) { }

inline su2double CVariable::GetGeometry_CrossTerm_Derivative_Flow(unsigned short val_var) { return 0.0;}

inline void CVariable::SetGeometry_CrossTerm_Derivative_Flow(unsigned short iDim, su2double der) { }

inline void CVariable::Set_OldSolution_Geometry(void) { }

inline su2double CVariable::Get_OldSolution_Geometry(unsigned short iDim) { return 0.0;}

inline su2double CVariable::GetCross_Term_Derivative(unsigned short iVar) { return 0.0; }

inline void CVariable::SetCross_Term_Derivative(unsigned short iVar, su2double der) { }

inline void CVariable::Set_BGSSolution(unsigned short iDim, su2double val_solution) { }

inline void CVariable::Set_BGSSolution_k(void) { }

inline su2double CVariable::Get_BGSSolution(unsigned short iDim) { return 0.0;}

inline su2double CVariable::Get_BGSSolution_k(unsigned short iDim) { return 0.0;}

inline void CVariable::Set_BGSSolution_Geometry(void) { }

inline su2double CVariable::Get_BGSSolution_Geometry(unsigned short iDim) { return 0.0;}

inline su2double* CVariable::GetSolution_Vel_Direct() { return NULL; }

inline void CVariable::SetSolution_Vel_Direct(su2double *val_solution_direct) { }

inline su2double* CVariable::GetSolution_Accel_Direct() { return NULL; }

inline void CVariable::SetSolution_Accel_Direct(su2double *val_solution_direct) { }

inline void CVariable::SetHarmonicBalance_Source(unsigned short val_var, su2double val_source) { }

inline su2double CVariable::GetHarmonicBalance_Source(unsigned short val_var) { return 0; }

inline void CVariable::SetEddyViscSens(su2double *val_EddyViscSens, unsigned short numTotalVar) { }

inline su2double *CVariable::GetEddyViscSens(void) { return NULL; }

inline void CVariable::SetSolution_time_n(void) { }

inline void CVariable::SetSolution_time_n(unsigned short val_var, su2double val_solution_time_n) { }

inline void CVariable::SetSolution_Vel(su2double *val_solution_vel) { }

inline void CVariable::SetSolution_Vel(unsigned short val_var, su2double val_solution_vel) { }

inline void CVariable::SetSolution_Vel_time_n(su2double *val_solution_vel_time_n) { }

inline void CVariable::SetSolution_Vel_time_n(void) { }

inline void CVariable::SetSolution_Vel_time_n(unsigned short val_var, su2double val_solution_vel_time_n) { }

inline void CVariable::Set_OldSolution_Accel(void) { }

inline void CVariable::Set_OldSolution_Vel(void) { }

inline su2double CVariable::GetSolution_time_n(unsigned short val_var) { return Solution_time_n[val_var]; }

inline su2double CVariable::GetSolution_Vel(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetSolution_Vel(void) { return NULL; }

inline su2double CVariable::GetSolution_Vel_time_n(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetSolution_Vel_time_n(void) { return NULL; }

inline void CVariable::SetSolution_Accel(su2double *val_solution_accel) { }

inline void CVariable::SetSolution_Accel(unsigned short val_var, su2double val_solution_accel) { }

inline void CVariable::SetSolution_Accel_time_n(su2double *val_solution_accel_time_n) { }

inline void CVariable::SetSolution_Accel_time_n(void) { }

inline void CVariable::SetSolution_Accel_time_n(unsigned short val_var, su2double val_solution_accel_time_n) { }

inline su2double CVariable::GetSolution_Accel(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetSolution_Accel(void) { return NULL; }

inline su2double CVariable::GetSolution_Accel_time_n(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetSolution_Accel_time_n(void) { return NULL; }

inline void CVariable::SetSolution_Pred(unsigned short val_var, su2double val_solution_pred) {  }

inline void CVariable::SetSolution_Pred(su2double *val_solution_pred) {  }

inline void CVariable::SetSolution_Pred(void) { }

inline su2double CVariable::GetSolution_Pred(unsigned short val_var) { return 0.0; }

inline su2double *CVariable::GetSolution_Pred(void) { return NULL; }

inline void CVariable::SetSolution_Pred_Old(unsigned short val_var, su2double val_solution_pred_old) {  }

inline void CVariable::SetSolution_Pred_Old(su2double *val_solution_pred_Old) {  }

inline void CVariable::SetSolution_Pred_Old(void) { }

inline su2double CVariable::GetSolution_Pred_Old(unsigned short val_var) { return 0.0; }

inline su2double *CVariable::GetSolution_Pred_Old(void) { return NULL; }

inline void CVariable::SetReference_Geometry(unsigned short iVar, su2double ref_geometry){  }

inline su2double *CVariable::GetReference_Geometry(void){ return NULL; }

inline su2double CVariable::GetReference_Geometry(unsigned short iVar){ return 0.0; }

inline void CVariable::SetPrestretch(unsigned short iVar, su2double val_prestretch) {  }

inline su2double *CVariable::GetPrestretch(void) { return NULL; }

inline su2double CVariable::GetPrestretch(unsigned short iVar) { return 0.0; }

inline void CVariable::Register_femSolution_time_n() { }

inline void CVariable::RegisterSolution_Vel(bool input) { }

inline void CVariable::RegisterSolution_Vel_time_n() { }

inline void CVariable::RegisterSolution_Accel(bool input) { }

inline void CVariable::RegisterSolution_Accel_time_n() { }

inline void CVariable::SetAdjointSolution_Vel(su2double *adj_sol) { }

inline void CVariable::GetAdjointSolution_Vel(su2double *adj_sol) { }

inline void CVariable::SetAdjointSolution_Vel_time_n(su2double *adj_sol) { }

inline void CVariable::GetAdjointSolution_Vel_time_n(su2double *adj_sol) { }

inline void CVariable::SetAdjointSolution_Accel(su2double *adj_sol) { }

inline void CVariable::GetAdjointSolution_Accel(su2double *adj_sol) { }

inline void CVariable::SetAdjointSolution_Accel_time_n(su2double *adj_sol) { }

inline void CVariable::GetAdjointSolution_Accel_time_n(su2double *adj_sol) { }

inline su2double CVariable::GetSolution_New(unsigned short val_var) { return 0.0; }


inline su2double CVariable::GetRoe_Dissipation(void) { return 0.0; }

inline void CVariable::SetRoe_Dissipation_FD(su2double val_wall_dist) { }

inline void CVariable::SetRoe_Dissipation_NTS(su2double val_delta, su2double val_const_DES) { }

inline su2double CVariable::GetDES_LengthScale(void) { return 0.0; }

inline void CVariable::SetDES_LengthScale(su2double val_des_lengthscale) { }

inline void CVariable::SetSolution_New(void) { }

inline void CVariable::AddSolution_New(unsigned short val_var, su2double val_solution) { }

inline void CVariable::SetRoe_Dissipation(su2double val_dissipation) { }

inline void CVariable::SetVortex_Tilting(su2double **PrimGrad_Flow, su2double* Vorticity, su2double LaminarViscosity) { }

inline su2double CVariable::GetVortex_Tilting() { return 0.0; }

inline su2double CEulerVariable::GetSolution_New(unsigned short val_var) { return Solution_New[val_var]; }

inline su2double CNSVariable::GetRoe_Dissipation(void) { return Roe_Dissipation; }

inline su2double CNSVariable::GetDES_LengthScale(void) { return DES_LengthScale; }

inline void CNSVariable::SetDES_LengthScale(su2double val_des_lengthscale) { DES_LengthScale = val_des_lengthscale; }

inline void CIncNSVariable::SetDES_LengthScale(su2double val_des_lengthscale) { DES_LengthScale = val_des_lengthscale; }

inline su2double CIncNSVariable::GetDES_LengthScale(void) { return DES_LengthScale; }

inline void CEulerVariable::SetSolution_New(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_New[iVar] = Solution[iVar];
}

inline void CEulerVariable::AddSolution_New(unsigned short val_var, su2double val_solution) {
  Solution_New[val_var] += val_solution;
}

inline su2double CEulerVariable::GetDensity(void) { return Solution[0]; }

inline su2double CEulerVariable::GetEnergy(void) { return Solution[nVar-1]/Solution[0]; };

inline su2double CEulerVariable::GetEnthalpy(void) { return Primitive[nDim+3]; }

inline su2double CEulerVariable::GetPressure(void) { return Primitive[nDim+1]; }

inline su2double CEulerVariable::GetSoundSpeed(void) { return Primitive[nDim+4]; }

inline su2double CEulerVariable::GetTemperature(void) { return Primitive[0]; }

inline su2double CEulerVariable::GetVelocity(unsigned short val_dim) { return Primitive[val_dim+1]; }

inline su2double CEulerVariable::GetVelocity2(void) { return Velocity2; }

inline bool CEulerVariable::SetDensity(void) {
  Primitive[nDim+2] = Solution[0];
  if (Primitive[nDim+2] > 0.0) return false;
  else return true;
}

inline bool CEulerVariable::SetPressure(su2double pressure) {
  Primitive[nDim+1] = pressure;
  if (Primitive[nDim+1] > 0.0) return false;
  else return true;
}

inline void CEulerVariable::SetVelocity(void) {
  Velocity2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Primitive[iDim+1] = Solution[iDim+1] / Solution[0];
    Velocity2 += Primitive[iDim+1]*Primitive[iDim+1];
  }
}

inline void CEulerVariable::SetEnthalpy(void) { Primitive[nDim+3] = (Solution[nVar-1] + Primitive[nDim+1]) / Solution[0]; }

inline bool CEulerVariable::SetSoundSpeed(su2double soundspeed2) {
  su2double radical = soundspeed2;
  if (radical < 0.0) return true;
  else {
    Primitive[nDim+4] = sqrt(radical);
    return false;
  }
}

inline bool CEulerVariable::SetTemperature(su2double temperature) {
  Primitive[0] = temperature;
  if (Primitive[0] > 0.0) return false;
  else return true;
}

inline void CEulerVariable::SetdPdrho_e(su2double dPdrho_e) {
  Secondary[0] = dPdrho_e;
}

inline void CEulerVariable::SetdPde_rho(su2double dPde_rho) {
  Secondary[1] = dPde_rho;
}

inline su2double CEulerVariable::GetPrimitive(unsigned short val_var) { return Primitive[val_var]; }

inline void CEulerVariable::SetPrimitive(unsigned short val_var, su2double val_prim) { Primitive[val_var] = val_prim; }

inline void CEulerVariable::SetPrimitive(su2double *val_prim) {
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
    Primitive[iVar] = val_prim[iVar];
}

inline su2double *CEulerVariable::GetPrimitive(void) { return Primitive; }

inline su2double CEulerVariable::GetSecondary(unsigned short val_var) { return Secondary[val_var]; }

inline void CEulerVariable::SetSecondary(unsigned short val_var, su2double val_secondary) { Secondary[val_var] = val_secondary; }

inline void CEulerVariable::SetSecondary(su2double *val_secondary) {
  for (unsigned short iVar = 0; iVar < nSecondaryVar; iVar++)
    Secondary[iVar] = val_secondary[iVar];
}

inline su2double *CEulerVariable::GetSecondary(void) { return Secondary; }

inline void CEulerVariable::SetVelocity_Old(su2double *val_velocity) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Solution_Old[iDim+1] = val_velocity[iDim]*Solution[0];
}

inline void CEulerVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] += val_value; }

inline void CEulerVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] -= val_value; }

inline su2double CEulerVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return Gradient_Primitive[val_var][val_dim]; }

inline su2double CEulerVariable::GetLimiter_Primitive(unsigned short val_var) { return Limiter_Primitive[val_var]; }

inline void CEulerVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] = val_value; }

inline void CEulerVariable::SetLimiter_Primitive(unsigned short val_var, su2double val_value) { Limiter_Primitive[val_var] = val_value; }

inline su2double **CEulerVariable::GetGradient_Primitive(void) { return Gradient_Primitive; }

inline su2double *CEulerVariable::GetLimiter_Primitive(void) { return Limiter_Primitive; }

inline void CEulerVariable::AddGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Secondary[val_var][val_dim] += val_value; }

inline void CEulerVariable::SubtractGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Secondary[val_var][val_dim] -= val_value; }

inline su2double CEulerVariable::GetGradient_Secondary(unsigned short val_var, unsigned short val_dim) { return Gradient_Secondary[val_var][val_dim]; }

inline su2double CEulerVariable::GetLimiter_Secondary(unsigned short val_var) { return Limiter_Secondary[val_var]; }

inline void CEulerVariable::SetGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Secondary[val_var][val_dim] = val_value; }

inline void CEulerVariable::SetLimiter_Secondary(unsigned short val_var, su2double val_value) { Limiter_Secondary[val_var] = val_value; }

inline su2double **CEulerVariable::GetGradient_Secondary(void) { return Gradient_Secondary; }

inline su2double *CEulerVariable::GetLimiter_Secondary(void) { return Limiter_Secondary; }

inline void CEulerVariable::SetHarmonicBalance_Source(unsigned short val_var, su2double val_source) { HB_Source[val_var] = val_source; }

inline su2double CEulerVariable::GetHarmonicBalance_Source(unsigned short val_var) { return HB_Source[val_var]; }

inline su2double CEulerVariable::GetPreconditioner_Beta() { return Precond_Beta; }

inline void CEulerVariable::SetPreconditioner_Beta(su2double val_Beta) { Precond_Beta = val_Beta; }

inline void CEulerVariable::SetWindGust( su2double* val_WindGust) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    WindGust[iDim] = val_WindGust[iDim];}

inline su2double* CEulerVariable::GetWindGust() { return WindGust;}

inline void CEulerVariable::SetWindGustDer( su2double* val_WindGustDer) {
  for (unsigned short iDim = 0; iDim < nDim+1; iDim++)
    WindGustDer[iDim] = val_WindGustDer[iDim];}

inline su2double* CEulerVariable::GetWindGustDer() { return WindGustDer;}

inline su2double CEulerVariable::Get_BGSSolution_k(unsigned short iDim) { return Solution_BGS_k[iDim];}

inline void CEulerVariable::Set_BGSSolution_k(void) { 
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_BGS_k[iVar] = Solution[iVar];
}

inline su2double CNSVariable::GetEddyViscosity(void) { return Primitive[nDim+6]; }

inline su2double CNSVariable::GetLaminarViscosity(void) { return Primitive[nDim+5]; }

inline su2double CNSVariable::GetThermalConductivity(void) { return Primitive[nDim+7]; }

inline su2double CNSVariable::GetSpecificHeatCp(void) { return Primitive[nDim+8]; }

inline su2double* CNSVariable::GetVorticity(void) { return Vorticity; }

inline su2double CNSVariable::GetStrainMag(void) { return StrainMag; }

inline void CNSVariable::SetLaminarViscosity(su2double laminarViscosity) {
  Primitive[nDim+5] = laminarViscosity;
}

inline void CNSVariable::SetThermalConductivity(su2double thermalConductivity) {
  Primitive[nDim+7] = thermalConductivity;
}

inline void CNSVariable::SetSpecificHeatCp(su2double Cp) {
  Primitive[nDim+8] = Cp;
}

inline void CNSVariable::SetdTdrho_e(su2double dTdrho_e) {
  Secondary[2] = dTdrho_e;
}

inline void CNSVariable::SetdTde_rho(su2double dTde_rho) {
  Secondary[3] = dTde_rho;
}

inline void CNSVariable::Setdmudrho_T(su2double dmudrho_T) {
  Secondary[4] = dmudrho_T;
}

inline void CNSVariable::SetdmudT_rho(su2double dmudT_rho) {
  Secondary[5] = dmudT_rho;
}

inline void CNSVariable::Setdktdrho_T(su2double dktdrho_T) {
  Secondary[6] = dktdrho_T;
}

inline void CNSVariable::SetdktdT_rho(su2double dktdT_rho) {
  Secondary[7] = dktdT_rho;
}
inline void CNSVariable::SetTauWall(su2double val_tau_wall) { Tau_Wall = val_tau_wall; }

inline su2double CNSVariable::GetTauWall(void) { return Tau_Wall; }

inline void CNSVariable::SetEddyViscosity(su2double eddy_visc) { Primitive[nDim+6] = eddy_visc; }

inline void CNSVariable::SetWallTemperature(su2double Temperature_Wall ) { Primitive[0] = Temperature_Wall; }

inline void CNSVariable::SetRoe_Dissipation(su2double val_dissipation) { Roe_Dissipation = val_dissipation; }

inline su2double CTurbSAVariable::GetVortex_Tilting() { return Vortex_Tilting; }

inline su2double *CAdjEulerVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline su2double *CAdjEulerVariable::GetObjFuncSource(void) { return ObjFuncSource; }

inline su2double *CAdjEulerVariable::GetIntBoundary_Jump(void) { return IntBoundary_Jump; }

inline void CAdjEulerVariable::SetForceProj_Vector(su2double *val_ForceProj_Vector) { for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

inline void CAdjEulerVariable::SetObjFuncSource(su2double *val_ObjFuncSource) { for (unsigned short iVar = 0; iVar < nVar; iVar++) ObjFuncSource[iVar] = val_ObjFuncSource[iVar]; }

inline void CAdjEulerVariable::SetIntBoundary_Jump(su2double *val_IntBoundary_Jump) { for (unsigned short iVar = 0; iVar < nVar; iVar++) IntBoundary_Jump[iVar] = val_IntBoundary_Jump[iVar]; }

inline void CAdjEulerVariable::SetPhi_Old(su2double *val_phi) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1]=val_phi[iDim]; };

inline void CAdjEulerVariable::SetHarmonicBalance_Source(unsigned short val_var, su2double val_source) { HB_Source[val_var] = val_source; }

inline su2double CAdjEulerVariable::GetHarmonicBalance_Source(unsigned short val_var) { return HB_Source[val_var]; }

inline su2double *CAdjNSVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline void CAdjNSVariable::SetForceProj_Vector(su2double *val_ForceProj_Vector) {  for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

inline void CAdjNSVariable::SetPhi_Old(su2double *val_phi) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1] = val_phi[iDim]; };

inline void CAdjNSVariable::SetVelSolutionOldDVector(void) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1] = ForceProj_Vector[iDim]; };

inline void CAdjNSVariable::SetVelSolutionDVector(void) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution[iDim+1] = ForceProj_Vector[iDim]; };

inline su2double CIncEulerVariable::GetDensity(void) { return Primitive[nDim+2]; }

inline su2double CIncEulerVariable::GetDensity_Old(void) { return Density_Old; }

inline su2double CIncEulerVariable::GetBetaInc2(void) { return Primitive[nDim+3]; }

inline su2double CIncEulerVariable::GetPressure(void) { return Primitive[0]; }

inline su2double CIncEulerVariable::GetTemperature(void) { return Primitive[nDim+1]; }

inline su2double CIncEulerVariable::GetVelocity(unsigned short val_dim) { return Primitive[val_dim+1]; }

inline su2double CIncEulerVariable::GetVelocity2(void) { return Velocity2; }

inline bool CIncEulerVariable::SetDensity(su2double val_density) { 
  Primitive[nDim+2] = val_density; 
  if (Primitive[nDim+2] > 0.0) return false;
  else return true;
}

inline void CIncEulerVariable::SetPressure(void) { Primitive[0] = Solution[0]; }

inline bool CIncEulerVariable::SetTemperature(su2double val_temperature) { 
  Primitive[nDim+1] = val_temperature;
  if (Primitive[nDim+1] > 0.0) return false;
  else return true; 
}

inline void CIncEulerVariable::SetVelocity(void) {
  Velocity2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Primitive[iDim+1] = Solution[iDim+1];
    Velocity2 += Primitive[iDim+1]*Primitive[iDim+1];
  }
}

inline void CIncEulerVariable::SetBetaInc2(su2double val_betainc2) { Primitive[nDim+3] = val_betainc2; }

inline su2double CIncEulerVariable::GetPrimitive(unsigned short val_var) { return Primitive[val_var]; }

inline void CIncEulerVariable::SetPrimitive(unsigned short val_var, su2double val_prim) { Primitive[val_var] = val_prim; }

inline void CIncEulerVariable::SetPrimitive(su2double *val_prim) {
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
    Primitive[iVar] = val_prim[iVar];
}

inline su2double *CIncEulerVariable::GetPrimitive(void) { return Primitive; }

inline void CIncEulerVariable::SetVelocity_Old(su2double *val_velocity) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Solution_Old[iDim+1] = val_velocity[iDim];
}

inline void CIncEulerVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] += val_value; }

inline void CIncEulerVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] -= val_value; }

inline su2double CIncEulerVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return Gradient_Primitive[val_var][val_dim]; }

inline su2double CIncEulerVariable::GetLimiter_Primitive(unsigned short val_var) { return Limiter_Primitive[val_var]; }

inline void CIncEulerVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] = val_value; }

inline void CIncEulerVariable::SetLimiter_Primitive(unsigned short val_var, su2double val_value) { Limiter_Primitive[val_var] = val_value; }

inline su2double **CIncEulerVariable::GetGradient_Primitive(void) { return Gradient_Primitive; }

inline su2double *CIncEulerVariable::GetLimiter_Primitive(void) { return Limiter_Primitive; }

inline void CIncEulerVariable::SetSpecificHeatCp(su2double val_Cp) {
  Primitive[nDim+7] = val_Cp;
}

inline void CIncEulerVariable::SetSpecificHeatCv(su2double val_Cv) {
  Primitive[nDim+8] = val_Cv;
}

inline su2double CIncEulerVariable::GetSpecificHeatCp(void) { return Primitive[nDim+7]; }

inline su2double CIncEulerVariable::GetSpecificHeatCv(void) { return Primitive[nDim+8]; }

inline su2double CIncEulerVariable::Get_BGSSolution_k(unsigned short iDim) { return Solution_BGS_k[iDim];}

inline void CIncEulerVariable::Set_BGSSolution_k(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_BGS_k[iVar] = Solution[iVar];
}

inline su2double CIncNSVariable::GetEddyViscosity(void) { return Primitive[nDim+5]; }

inline su2double CIncNSVariable::GetLaminarViscosity(void) { return Primitive[nDim+4]; }

inline su2double CIncNSVariable::GetThermalConductivity(void) { return Primitive[nDim+6]; }

inline su2double* CIncNSVariable::GetVorticity(void) { return Vorticity; }

inline su2double CIncNSVariable::GetStrainMag(void) { return StrainMag; }

inline void CIncNSVariable::SetLaminarViscosity(su2double val_laminar_viscosity_inc) { Primitive[nDim+4] = val_laminar_viscosity_inc; }

inline void CIncNSVariable::SetEddyViscosity(su2double eddy_visc) { Primitive[nDim+5] = eddy_visc; }

inline void CIncNSVariable::SetThermalConductivity(su2double val_thermal_conductivity) {
  Primitive[nDim+6] = val_thermal_conductivity;
}

inline su2double CTransLMVariable::GetIntermittency() { return Solution[0]; }

inline void CTransLMVariable::SetGammaSep(su2double gamma_sep_in) {gamma_sep = gamma_sep_in;}

inline void CFEAVariable::SetStress_FEM(unsigned short iVar, su2double val_stress) { Stress[iVar] = val_stress; }

inline void CFEAVariable::AddStress_FEM(unsigned short iVar, su2double val_stress) { Stress[iVar] += val_stress; }

inline su2double *CFEAVariable::GetStress_FEM(void) { return Stress; }

inline void CFEAVariable::Add_SurfaceLoad_Res(su2double *val_surfForce) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Residual_Ext_Surf[iVar] += val_surfForce[iVar];
}

inline void CFEAVariable::Set_SurfaceLoad_Res(unsigned short iVar, su2double val_surfForce) { Residual_Ext_Surf[iVar] = val_surfForce;}

inline su2double *CFEAVariable::Get_SurfaceLoad_Res(void) { return Residual_Ext_Surf;}

inline su2double CFEAVariable::Get_SurfaceLoad_Res(unsigned short iVar) { return Residual_Ext_Surf[iVar];}

inline void CFEAVariable::Clear_SurfaceLoad_Res(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  Residual_Ext_Surf[iVar] = 0.0;
}

inline void CFEAVariable::Set_SurfaceLoad_Res_n(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  Residual_Ext_Surf_n[iVar] = Residual_Ext_Surf[iVar];
}

inline su2double CFEAVariable::Get_SurfaceLoad_Res_n(unsigned short iVar) { return Residual_Ext_Surf_n[iVar];}

inline void CFEAVariable::Add_BodyForces_Res(su2double *val_bodyForce) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Residual_Ext_Body[iVar] += val_bodyForce[iVar];
}

inline su2double *CFEAVariable::Get_BodyForces_Res(void) { return Residual_Ext_Body;}

inline su2double CFEAVariable::Get_BodyForces_Res(unsigned short iVar) { return Residual_Ext_Body[iVar];}

inline void CFEAVariable::Clear_BodyForces_Res(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  Residual_Ext_Body[iVar] = 0.0;
}

inline void CFEAVariable::Set_FlowTraction(su2double *val_flowTraction) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    FlowTraction[iVar] = val_flowTraction[iVar];
}

inline void CFEAVariable::Add_FlowTraction(su2double *val_flowTraction) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    FlowTraction[iVar] += val_flowTraction[iVar];
}

inline su2double *CFEAVariable::Get_FlowTraction(void) { return FlowTraction;}

inline su2double CFEAVariable::Get_FlowTraction(unsigned short iVar) { return FlowTraction[iVar];}

inline void CFEAVariable::Clear_FlowTraction(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  FlowTraction[iVar] = 0.0;
}

inline void CFEAVariable::Set_FlowTraction_n(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  FlowTraction_n[iVar] = FlowTraction[iVar];
}

inline su2double CFEAVariable::Get_FlowTraction_n(unsigned short iVar) { return FlowTraction_n[iVar];}

inline void CFEAVariable::SetSolution_time_n(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  Solution_time_n[iVar] = Solution[iVar];
}

inline void CFEAVariable::SetSolution_time_n(su2double *val_solution_time_n) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  Solution_time_n[iVar] = val_solution_time_n[iVar];
}

inline void CFEAVariable::SetSolution_time_n(unsigned short val_var, su2double val_solution_time_n) { Solution_time_n[val_var] = val_solution_time_n; }

inline void CFEAVariable::SetSolution_Vel(unsigned short val_var, su2double val_solution_vel) { Solution_Vel[val_var] = val_solution_vel; }

inline void CFEAVariable::SetSolution_Vel(su2double *val_solution_vel) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  Solution_Vel[iVar] = val_solution_vel[iVar];
}

inline void CFEAVariable::SetSolution_Vel_time_n(unsigned short val_var, su2double val_solution_vel_time_n) { Solution_Vel_time_n[val_var] = val_solution_vel_time_n; }

inline void CFEAVariable::SetSolution_Vel_time_n(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  Solution_Vel_time_n[iVar] = Solution_Vel[iVar];
}

inline void CFEAVariable::SetSolution_Vel_time_n(su2double *val_solution_vel_time_n) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  Solution_Vel_time_n[iVar] = val_solution_vel_time_n[iVar];
}

inline void CFEAVariable::SetSolution_Accel(unsigned short val_var, su2double val_solution_accel) { Solution_Accel[val_var] = val_solution_accel;  }

inline void CFEAVariable::SetSolution_Accel(su2double *val_solution_accel) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  Solution_Accel[iVar] = val_solution_accel[iVar];
}

inline void CFEAVariable::SetSolution_Accel_time_n(unsigned short val_var, su2double val_solution_accel_time_n) { Solution_Accel_time_n[val_var] = val_solution_accel_time_n; }

inline void CFEAVariable::SetSolution_Accel_time_n(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  Solution_Accel_time_n[iVar] = Solution_Accel[iVar];
}

inline void CFEAVariable::SetSolution_Accel_time_n(su2double *val_solution_accel_time_n) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)  Solution_Accel_time_n[iVar] = val_solution_accel_time_n[iVar];
}

inline void CFEAVariable::SetSolution_Pred(unsigned short val_var, su2double val_solution_pred) { Solution_Pred[val_var] = val_solution_pred;  }

inline void CFEAVariable::SetSolution_Pred(su2double *val_solution_pred) { Solution_Pred = val_solution_pred;  }

inline void CFEAVariable::SetSolution_Pred(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Pred[iVar] = Solution[iVar];
}

inline void CFEAVariable::SetSolution_Pred_Old(unsigned short val_var, su2double val_solution_pred_old) { Solution_Pred_Old[val_var] = val_solution_pred_old;  }

inline void CFEAVariable::SetSolution_Pred_Old(su2double *val_solution_pred_Old) { Solution_Pred_Old = val_solution_pred_Old;  }

inline void CFEAVariable::SetSolution_Pred_Old(void) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Pred_Old[iVar] = Solution_Pred[iVar];
}

inline su2double *CFEAVariable::GetSolution_Vel(void) { return Solution_Vel; }

inline su2double CFEAVariable::GetSolution_Vel(unsigned short val_var) { return Solution_Vel[val_var]; }

inline su2double *CFEAVariable::GetSolution_Vel_time_n(void) { return Solution_Vel_time_n; }

inline su2double CFEAVariable::GetSolution_Vel_time_n(unsigned short val_var) { return Solution_Vel_time_n[val_var]; }

inline su2double *CFEAVariable::GetSolution_Accel(void) { return Solution_Accel; }

inline su2double CFEAVariable::GetSolution_Accel(unsigned short val_var) { return Solution_Accel[val_var]; }

inline su2double *CFEAVariable::GetSolution_Accel_time_n(void) { return Solution_Accel_time_n; }

inline su2double CFEAVariable::GetSolution_Accel_time_n(unsigned short val_var) { return Solution_Accel_time_n[val_var]; }

inline su2double *CFEAVariable::GetSolution_Pred(void) { return Solution_Pred; }

inline su2double CFEAVariable::GetSolution_Pred(unsigned short val_var) { return Solution_Pred[val_var]; }

inline su2double *CFEAVariable::GetSolution_Pred_Old(void) { return Solution_Pred_Old; }

inline su2double CFEAVariable::GetSolution_Pred_Old(unsigned short val_var) { return Solution_Pred_Old[val_var]; }

inline void CFEAVariable::SetVonMises_Stress(su2double val_stress) { VonMises_Stress = val_stress; }

inline su2double CFEAVariable::GetVonMises_Stress(void) { return VonMises_Stress; }

inline void CFEAVariable::SetReference_Geometry(unsigned short iVar, su2double ref_geometry){ Reference_Geometry[iVar] = ref_geometry;}

inline su2double *CFEAVariable::GetReference_Geometry(void){ return Reference_Geometry; }

inline su2double CFEAVariable::GetReference_Geometry(unsigned short iVar){ return Reference_Geometry[iVar]; }

inline void CFEAVariable::Register_femSolution_time_n() {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
	    AD::RegisterInput(Solution_time_n[iVar]);
}

inline void CFEAVariable::RegisterSolution_Vel(bool input) {
	  if (input) {
	    for (unsigned short iVar = 0; iVar < nVar; iVar++)
	      AD::RegisterInput(Solution_Vel[iVar]);
	  }
	  else { for (unsigned short iVar = 0; iVar < nVar; iVar++)
	      AD::RegisterOutput(Solution_Vel[iVar]);}
}

inline void CFEAVariable::RegisterSolution_Vel_time_n() {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
	    AD::RegisterInput(Solution_Vel_time_n[iVar]);
}

inline void CFEAVariable::RegisterSolution_Accel(bool input) {
	  if (input) {
	    for (unsigned short iVar = 0; iVar < nVar; iVar++)
	      AD::RegisterInput(Solution_Accel[iVar]);
	  }
	  else { for (unsigned short iVar = 0; iVar < nVar; iVar++)
	      AD::RegisterOutput(Solution_Accel[iVar]);}
}

inline void CFEAVariable::RegisterSolution_Accel_time_n() {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
	    AD::RegisterInput(Solution_Accel_time_n[iVar]);
}

inline su2double CFEAVariable::Get_BGSSolution_k(unsigned short iDim) { return Solution_BGS_k[iDim];}

inline void CFEAVariable::Set_BGSSolution_k(void) { 
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_BGS_k[iVar] = Solution[iVar];
}

inline void CFEAVariable::SetPrestretch(unsigned short iVar, su2double val_prestretch) { Prestretch[iVar] = val_prestretch;}

inline su2double *CFEAVariable::GetPrestretch(void) { return Prestretch; }

inline su2double CFEAVariable::GetPrestretch(unsigned short iVar) { return Prestretch[iVar]; }

inline void CFEABoundVariable::SetTraction(unsigned short iVar, unsigned short jVar, su2double val_traction) { Traction[iVar][jVar] = val_traction; }

inline void CFEABoundVariable::AddTraction(unsigned short iVar, unsigned short jVar, su2double val_traction) { Traction[iVar][jVar] += val_traction; }

inline su2double **CFEABoundVariable::GetTraction(void) { return Traction; }

inline su2double* CWaveVariable::GetSolution_Direct() { return Solution_Direct;}

inline void CWaveVariable::SetSolution_Direct(su2double *val_solution_direct) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Direct[iVar] += val_solution_direct[iVar];}

inline su2double* CPotentialVariable::GetChargeDensity() { return Charge_Density;}

inline void CPotentialVariable::SetChargeDensity(su2double positive_charge, su2double negative_charge) {Charge_Density[0] = positive_charge; Charge_Density[1] = negative_charge;}

inline void CTurbSAVariable::SetHarmonicBalance_Source(unsigned short val_var, su2double val_source) { HB_Source[val_var] = val_source; }

inline su2double CTurbSAVariable::GetHarmonicBalance_Source(unsigned short val_var) { return HB_Source[val_var]; }

inline su2double CTurbSAVariable::GetGammaBC(void) { return gamma_BC; }

inline void CTurbSAVariable::SetGammaBC(su2double val_gamma) { gamma_BC = val_gamma; }

inline su2double CTurbSAVariable::GetDES_LengthScale(void) { return DES_LengthScale; }

inline void CTurbSAVariable::SetDES_LengthScale(su2double val_des_lengthscale) { DES_LengthScale = val_des_lengthscale; }

inline su2double CTurbSSTVariable::GetF1blending() { return F1; }

inline su2double CTurbSSTVariable::GetF2blending() { return F2; }

inline su2double CTurbSSTVariable::GetCrossDiff() { return CDkw; }

inline void CAdjTurbVariable::SetEddyViscSens(su2double *val_EddyViscSens, unsigned short numTotalVar) {
  for (unsigned short iVar = 0; iVar < numTotalVar; iVar++) {
    EddyViscSens[iVar] = val_EddyViscSens[iVar];}
}

inline su2double *CAdjTurbVariable::GetEddyViscSens(void) { return EddyViscSens; }

inline void CVariable::RegisterSolution(bool input) {
  if (input) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      AD::RegisterInput(Solution[iVar]);
  }
  else { for (unsigned short iVar = 0; iVar < nVar; iVar++)
      AD::RegisterOutput(Solution[iVar]);}
}

inline void CVariable::RegisterSolution_time_n() {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    AD::RegisterInput(Solution_time_n[iVar]);
}

inline void CVariable::RegisterSolution_time_n1() {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    AD::RegisterInput(Solution_time_n1[iVar]);
}

inline void CVariable::SetAdjointSolution(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
        SU2_TYPE::SetDerivative(Solution[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
}


inline void CVariable::GetAdjointSolution(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution[iVar]);
    }
}

inline void CVariable::SetAdjointSolution_time_n(su2double *adj_sol) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_time_n[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
}


inline void CVariable::GetAdjointSolution_time_n(su2double *adj_sol) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_time_n[iVar]);
  }
}

inline void CVariable::SetAdjointSolution_time_n1(su2double *adj_sol) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_time_n1[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
}


inline void CVariable::GetAdjointSolution_time_n1(su2double *adj_sol) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_time_n1[iVar]);
  }
}
inline void CVariable::SetDual_Time_Derivative(unsigned short iVar, su2double der) {}

inline void CDiscAdjVariable::SetDual_Time_Derivative(unsigned short iVar, su2double der) {DualTime_Derivative[iVar] = der;}

inline void CVariable::SetDual_Time_Derivative_n(unsigned short iVar, su2double der) {}

inline void CDiscAdjVariable::SetDual_Time_Derivative_n(unsigned short iVar, su2double der) {DualTime_Derivative_n[iVar] = der;}

inline su2double CVariable::GetDual_Time_Derivative(unsigned short iVar) { return 0.0;}

inline su2double CDiscAdjVariable::GetDual_Time_Derivative(unsigned short iVar) { return DualTime_Derivative[iVar];}

inline su2double CVariable::GetDual_Time_Derivative_n(unsigned short iVar) { return 0.0;}

inline su2double CDiscAdjVariable::GetDual_Time_Derivative_n(unsigned short iVar) { return DualTime_Derivative_n[iVar];}

inline void CVariable::SetSensitivity(unsigned short iDim, su2double val) {}

inline su2double CVariable::GetSensitivity(unsigned short iDim) { return 0.0; }

inline void CDiscAdjVariable::SetSensitivity(unsigned short iDim, su2double val) {Sensitivity[iDim] = val;}

inline su2double CDiscAdjVariable::GetSensitivity(unsigned short iDim) { return Sensitivity[iDim];}

inline su2double* CDiscAdjVariable::GetSolution_Direct() { return Solution_Direct; }

inline void CDiscAdjVariable::SetSolution_Direct(su2double *val_solution_direct) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Solution_Direct[iVar] = val_solution_direct[iVar];
  }
}

inline su2double* CDiscAdjVariable::GetGeometry_Direct() { return Geometry_Direct; }

inline su2double CDiscAdjVariable::GetGeometry_Direct(unsigned short val_dim) { return Geometry_Direct[val_dim]; }

inline void CDiscAdjVariable::SetGeometry_Direct(su2double *val_geometry_direct) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    Geometry_Direct[iDim] = val_geometry_direct[iDim];
  }
}

inline su2double CDiscAdjVariable::GetSolution_Geometry(unsigned short val_var) { return Solution_Geometry[val_var];}

inline void CDiscAdjVariable::SetSolution_Geometry(su2double *val_solution_geometry) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    Solution_Geometry[iDim] = val_solution_geometry[iDim];
  }
}

inline void CDiscAdjVariable::SetSolution_Geometry(unsigned short val_var, su2double val_solution_geometry) {
    Solution_Geometry[val_var] = val_solution_geometry;
}

inline void CDiscAdjVariable::Set_OldSolution_Geometry(void) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    Solution_Geometry_Old[iDim] = Solution_Geometry[iDim];
  }
}

inline su2double CDiscAdjVariable::GetGeometry_CrossTerm_Derivative(unsigned short val_var) { return Geometry_CrossTerm_Derivative[val_var];}

inline void CDiscAdjVariable::SetGeometry_CrossTerm_Derivative(unsigned short iDim, su2double der) { Geometry_CrossTerm_Derivative[iDim] = der;}

inline su2double CDiscAdjVariable::GetGeometry_CrossTerm_Derivative_Flow(unsigned short val_var) { return Geometry_CrossTerm_Derivative_Flow[val_var];}

inline void CDiscAdjVariable::SetGeometry_CrossTerm_Derivative_Flow(unsigned short iDim, su2double der) { Geometry_CrossTerm_Derivative_Flow[iDim] = der;}

inline su2double CDiscAdjFEAVariable::GetGeometry_CrossTerm_Derivative(unsigned short val_var) { return Geometry_CrossTerm_Derivative[val_var];}

inline void CDiscAdjFEAVariable::SetGeometry_CrossTerm_Derivative(unsigned short iDim, su2double der) { Geometry_CrossTerm_Derivative[iDim] = der;}

inline su2double CDiscAdjVariable::Get_OldSolution_Geometry(unsigned short iDim) { return Solution_Geometry_Old[iDim];}

inline void CVariable::SetDynamic_Derivative(unsigned short iVar, su2double der) { }

inline void CVariable::SetDynamic_Derivative_n(unsigned short iVar, su2double der) { }

inline su2double CVariable::GetDynamic_Derivative(unsigned short iVar) { return 0.0; }

inline su2double CVariable::GetDynamic_Derivative_n(unsigned short iVar) { return 0.0; }

inline void CVariable::SetDynamic_Derivative_Vel(unsigned short iVar, su2double der) { }

inline void CVariable::SetDynamic_Derivative_Vel_n(unsigned short iVar, su2double der) { }

inline su2double CVariable::GetDynamic_Derivative_Vel(unsigned short iVar) { return 0.0; }

inline su2double CVariable::GetDynamic_Derivative_Vel_n(unsigned short iVar) { return 0.0; }

inline void CVariable::SetDynamic_Derivative_Accel(unsigned short iVar, su2double der) { }

inline void CVariable::SetDynamic_Derivative_Accel_n(unsigned short iVar, su2double der) { }

inline su2double CVariable::GetDynamic_Derivative_Accel(unsigned short iVar) { return 0.0; }

inline su2double CVariable::GetDynamic_Derivative_Accel_n(unsigned short iVar) { return 0.0; }

inline su2double CVariable::GetSolution_Old_Vel(unsigned short iVar){ return 0.0; }

inline su2double CVariable::GetSolution_Old_Accel(unsigned short iVar){ return 0.0; }

inline void CDiscAdjFEAVariable::SetDynamic_Derivative(unsigned short iVar, su2double der) { Dynamic_Derivative[iVar] = der; }

inline void CDiscAdjFEAVariable::SetDynamic_Derivative_n(unsigned short iVar, su2double der) { Dynamic_Derivative_n[iVar] = der; }

inline su2double CDiscAdjFEAVariable::GetDynamic_Derivative(unsigned short iVar) { return Dynamic_Derivative[iVar]; }

inline su2double CDiscAdjFEAVariable::GetDynamic_Derivative_n(unsigned short iVar) { return Dynamic_Derivative_n[iVar]; }

inline void CDiscAdjFEAVariable::SetDynamic_Derivative_Vel(unsigned short iVar, su2double der) { Dynamic_Derivative_Vel[iVar] = der; }

inline void CDiscAdjFEAVariable::SetDynamic_Derivative_Vel_n(unsigned short iVar, su2double der) { Dynamic_Derivative_Vel_n[iVar] = der; }

inline su2double CDiscAdjFEAVariable::GetDynamic_Derivative_Vel(unsigned short iVar) { return Dynamic_Derivative_Vel[iVar]; }

inline su2double CDiscAdjFEAVariable::GetDynamic_Derivative_Vel_n(unsigned short iVar) { return Dynamic_Derivative_Vel_n[iVar]; }

inline void CDiscAdjFEAVariable::SetDynamic_Derivative_Accel(unsigned short iVar, su2double der) { Dynamic_Derivative_Accel[iVar] = der; }

inline void CDiscAdjFEAVariable::SetDynamic_Derivative_Accel_n(unsigned short iVar, su2double der) { Dynamic_Derivative_Accel_n[iVar] = der; }

inline su2double CDiscAdjFEAVariable::GetDynamic_Derivative_Accel(unsigned short iVar) { return Dynamic_Derivative_Accel[iVar]; }

inline su2double CDiscAdjFEAVariable::GetDynamic_Derivative_Accel_n(unsigned short iVar) { return Dynamic_Derivative_Accel_n[iVar]; }

inline su2double CDiscAdjFEAVariable::GetSolution_Old_Vel(unsigned short iVar){ return Solution_Old_Vel[iVar]; }

inline su2double CDiscAdjFEAVariable::GetSolution_Old_Accel(unsigned short iVar){ return Solution_Old_Accel[iVar]; }


inline void CDiscAdjFEAVariable::SetSolution_Accel(su2double *val_solution_accel) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)	Solution_Accel[iVar] = val_solution_accel[iVar];
}

inline void CDiscAdjFEAVariable::SetSolution_Vel(su2double *val_solution_vel) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)	Solution_Vel[iVar] = val_solution_vel[iVar];
}

inline void CDiscAdjFEAVariable::Set_OldSolution_Accel(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)  Solution_Old_Accel[iVar] = Solution_Accel[iVar];
 }

inline void CDiscAdjFEAVariable::Set_OldSolution_Vel(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)  Solution_Old_Vel[iVar] = Solution_Vel[iVar];
}

inline void CDiscAdjFEAVariable::SetSolution_Accel_time_n(su2double *val_solution_accel_time_n) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)	Solution_Accel_time_n[iVar] = val_solution_accel_time_n[iVar];
}

inline void CDiscAdjFEAVariable::SetSolution_Vel_time_n(su2double *val_solution_vel_time_n) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)	Solution_Vel_time_n[iVar] = val_solution_vel_time_n[iVar];
}

inline void CDiscAdjFEAVariable::SetSolution_Direct(su2double *val_solution_direct) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++){
    Solution_Direct[iVar] = val_solution_direct[iVar];
  }
}

inline void CDiscAdjFEAVariable::SetSensitivity(unsigned short iDim, su2double val){Sensitivity[iDim] = val;}

inline su2double CDiscAdjFEAVariable::GetSensitivity(unsigned short iDim){return Sensitivity[iDim];}

inline su2double* CDiscAdjFEAVariable::GetSolution_Direct() { return Solution_Direct; }

inline su2double* CDiscAdjFEAVariable::GetSolution_Vel_Direct() { return Solution_Direct_Vel; }

inline void CDiscAdjFEAVariable::SetSolution_Vel_Direct(su2double *val_solution_direct) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++){
	    Solution_Direct_Vel[iVar] = val_solution_direct[iVar];
	  }
 }

inline su2double* CDiscAdjFEAVariable::GetSolution_Accel_Direct() { return Solution_Direct_Accel; }

inline void CDiscAdjFEAVariable::SetSolution_Accel_Direct(su2double *val_solution_direct) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++){
	    Solution_Direct_Accel[iVar] = val_solution_direct[iVar];
	  }
}

inline void CDiscAdjFEAVariable::SetSolution_time_n(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)	Solution_time_n[iVar] = Solution[iVar];
}

inline su2double CDiscAdjFEAVariable::GetSolution_Accel(unsigned short val_var) { return Solution_Accel[val_var]; }
inline su2double CDiscAdjFEAVariable::GetSolution_Accel_time_n(unsigned short val_var) { return Solution_Accel_time_n[val_var]; }
inline su2double CDiscAdjFEAVariable::GetSolution_Vel_time_n(unsigned short val_var) { return Solution_Vel_time_n[val_var]; }
inline su2double CDiscAdjFEAVariable::GetSolution_Vel(unsigned short val_var) { return Solution_Vel[val_var]; }

inline void CFEAVariable::SetAdjointSolution_Vel(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
        SU2_TYPE::SetDerivative(Solution_Vel[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
}

inline void CFEAVariable::GetAdjointSolution_Vel(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++){
        adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Vel[iVar]);
    }
}

inline void CFEAVariable::SetAdjointSolution_Vel_time_n(su2double *adj_sol) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
	      SU2_TYPE::SetDerivative(Solution_Vel_time_n[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
}

inline void CFEAVariable::GetAdjointSolution_Vel_time_n(su2double *adj_sol) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++){
	      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Vel_time_n[iVar]);
	  }
}

inline void CFEAVariable::SetAdjointSolution_Accel(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
        SU2_TYPE::SetDerivative(Solution_Accel[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
}

inline void CFEAVariable::GetAdjointSolution_Accel(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++){
        adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Accel[iVar]);
    }
}

inline void CFEAVariable::SetAdjointSolution_Accel_time_n(su2double *adj_sol) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
	      SU2_TYPE::SetDerivative(Solution_Accel_time_n[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
}

inline void CFEAVariable::GetAdjointSolution_Accel_time_n(su2double *adj_sol) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++){
	      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Accel_time_n[iVar]);
	  }
}

inline su2double CDiscAdjVariable::GetCross_Term_Derivative(unsigned short iVar) { return Cross_Term_Derivative[iVar]; }

inline void CDiscAdjVariable::SetCross_Term_Derivative(unsigned short iVar, su2double der) { Cross_Term_Derivative[iVar] = der; }

inline su2double CDiscAdjFEAVariable::GetCross_Term_Derivative(unsigned short iVar) { return Cross_Term_Derivative[iVar]; }

inline void CDiscAdjFEAVariable::SetCross_Term_Derivative(unsigned short iVar, su2double der) { Cross_Term_Derivative[iVar] = der; }

inline void CDiscAdjVariable::Set_BGSSolution_k(void) { 
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_BGS_k[iVar] = Solution_BGS[iVar];
}

inline void CDiscAdjVariable::Set_BGSSolution(unsigned short iDim, su2double val_solution) { 
  Solution_BGS[iDim] = val_solution;
}

inline su2double CDiscAdjVariable::Get_BGSSolution(unsigned short iDim) { return Solution_BGS[iDim];}

inline su2double CDiscAdjVariable::Get_BGSSolution_k(unsigned short iDim) { return Solution_BGS_k[iDim];}

inline void CDiscAdjVariable::Set_BGSSolution_Geometry(void) { 
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Solution_Geometry_BGS_k[iDim] = Solution_Geometry[iDim];
}

inline su2double CDiscAdjVariable::Get_BGSSolution_Geometry(unsigned short iDim) { return Solution_Geometry_BGS_k[iDim];}

inline void CDiscAdjFEAVariable::Set_BGSSolution(unsigned short iDim, su2double val_solution) { 
  Solution_BGS[iDim] = val_solution;
}

inline void CDiscAdjFEAVariable::Set_BGSSolution_k(void) { 
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Solution_BGS_k[iDim] = Solution_BGS[iDim];
}

inline su2double CDiscAdjFEAVariable::Get_BGSSolution(unsigned short iDim) { return Solution_BGS[iDim];}

inline su2double CDiscAdjFEAVariable::Get_BGSSolution_k(unsigned short iDim) { return Solution_BGS_k[iDim];}

