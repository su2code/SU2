inline su2double CFEM_DG_NSSolver::GetViscosity_Inf(void) { return Viscosity_Inf; }

inline su2double CFEM_DG_NSSolver::GetTke_Inf(void) { return Tke_Inf; }

inline su2double CFEM_DG_NSSolver::GetCL_Visc(unsigned short val_marker) { return CL_Visc[val_marker]; }

inline su2double CFEM_DG_NSSolver::GetCMz_Visc(unsigned short val_marker) { return CMz_Visc[val_marker]; }

inline su2double CFEM_DG_NSSolver::GetCSF_Visc(unsigned short val_marker) { return CSF_Visc[val_marker]; }

inline su2double CFEM_DG_NSSolver::GetCD_Visc(unsigned short val_marker) { return CD_Visc[val_marker]; }

inline su2double CFEM_DG_NSSolver::GetAllBound_CL_Visc() { return AllBound_CL_Visc; }

inline su2double CFEM_DG_NSSolver::GetAllBound_CSF_Visc() { return AllBound_CSF_Visc; }

inline su2double CFEM_DG_NSSolver::GetAllBound_CD_Visc() { return AllBound_CD_Visc; }

inline su2double CFEM_DG_NSSolver::GetStrainMag_Max(void) { return StrainMag_Max; }

inline su2double CFEM_DG_NSSolver::GetOmega_Max(void) { return Omega_Max; }

inline void CFEM_DG_NSSolver::SetStrainMag_Max(su2double val_strainmag_max) { StrainMag_Max = val_strainmag_max; }

inline void CFEM_DG_NSSolver::SetOmega_Max(su2double val_omega_max) { Omega_Max = val_omega_max; }

