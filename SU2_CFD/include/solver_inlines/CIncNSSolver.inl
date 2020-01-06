inline su2double CIncNSSolver::GetViscosity_Inf(void) { return Viscosity_Inf; }

inline su2double CIncNSSolver::GetTke_Inf(void) { return Tke_Inf; }

inline su2double CIncNSSolver::GetSurface_HF_Visc(unsigned short val_marker) { return Surface_HF_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_MaxHF_Visc(unsigned short val_marker) { return Surface_MaxHF_Visc[val_marker]; }

inline su2double CIncNSSolver::GetCL_Visc(unsigned short val_marker) { return CL_Visc[val_marker]; }

inline su2double CIncNSSolver::GetCSF_Visc(unsigned short val_marker) { return CSF_Visc[val_marker]; }

inline su2double CIncNSSolver::GetCD_Visc(unsigned short val_marker) { return CD_Visc[val_marker]; }

inline su2double CIncNSSolver::GetAllBound_CL_Visc() { return AllBound_CL_Visc; }

inline su2double CIncNSSolver::GetAllBound_CSF_Visc() { return AllBound_CSF_Visc; }

inline su2double CIncNSSolver::GetAllBound_CD_Visc() { return AllBound_CD_Visc; }

inline su2double CIncNSSolver::GetAllBound_CEff_Visc() { return AllBound_CEff_Visc; }

inline su2double CIncNSSolver::GetAllBound_CMx_Visc() { return AllBound_CMx_Visc; }

inline su2double CIncNSSolver::GetAllBound_CMy_Visc() { return AllBound_CMy_Visc; }

inline su2double CIncNSSolver::GetAllBound_CMz_Visc() { return AllBound_CMz_Visc; }

inline su2double CIncNSSolver::GetAllBound_CoPx_Visc() { return AllBound_CoPx_Visc; }

inline su2double CIncNSSolver::GetAllBound_CoPy_Visc() { return AllBound_CoPy_Visc; }

inline su2double CIncNSSolver::GetAllBound_CoPz_Visc() { return AllBound_CoPz_Visc; }

inline su2double CIncNSSolver::GetAllBound_CFx_Visc() { return AllBound_CFx_Visc; }

inline su2double CIncNSSolver::GetAllBound_CFy_Visc() { return AllBound_CFy_Visc; }

inline su2double CIncNSSolver::GetAllBound_CFz_Visc() { return AllBound_CFz_Visc; }

inline su2double CIncNSSolver::GetSurface_CL_Visc(unsigned short val_marker) { return Surface_CL_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CD_Visc(unsigned short val_marker) { return Surface_CD_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CSF_Visc(unsigned short val_marker) { return Surface_CSF_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CEff_Visc(unsigned short val_marker) { return Surface_CEff_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CFx_Visc(unsigned short val_marker) { return Surface_CFx_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CFy_Visc(unsigned short val_marker) { return Surface_CFy_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CFz_Visc(unsigned short val_marker) { return Surface_CFz_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CMx_Visc(unsigned short val_marker) { return Surface_CMx_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CMy_Visc(unsigned short val_marker) { return Surface_CMy_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CMz_Visc(unsigned short val_marker) { return Surface_CMz_Visc[val_marker]; }

inline su2double CIncNSSolver::GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim) { return CSkinFriction[val_marker][val_dim][val_vertex]; }

inline su2double CIncNSSolver::GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) { return HeatFlux[val_marker][val_vertex]; }

inline su2double CIncNSSolver::GetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex) { return HeatFluxTarget[val_marker][val_vertex]; }

inline void CIncNSSolver::SetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_heat) { HeatFluxTarget[val_marker][val_vertex] = val_heat; }

inline su2double CIncNSSolver::GetYPlus(unsigned short val_marker, unsigned long val_vertex) { return YPlus[val_marker][val_vertex]; }

inline su2double CIncNSSolver::GetStrainMag_Max(void) { return StrainMag_Max; }

inline su2double CIncNSSolver::GetOmega_Max(void) { return Omega_Max; }

inline void CIncNSSolver::SetStrainMag_Max(su2double val_strainmag_max) { StrainMag_Max = val_strainmag_max; }

inline void CIncNSSolver::SetOmega_Max(su2double val_omega_max) { Omega_Max = val_omega_max; }

inline su2double CIncNSSolver::GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var) { return HeatConjugateVar[val_marker][val_vertex][pos_var]; }

inline void CIncNSSolver::SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var, su2double relaxation_factor, su2double val_var) {
  HeatConjugateVar[val_marker][val_vertex][pos_var] = relaxation_factor*val_var + (1.0-relaxation_factor)*HeatConjugateVar[val_marker][val_vertex][pos_var]; }

