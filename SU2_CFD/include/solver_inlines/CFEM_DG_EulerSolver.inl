inline CFluidModel* CFEM_DG_EulerSolver::GetFluidModel(void) { return FluidModel;}

inline su2double* CFEM_DG_EulerSolver::GetVecSolDOFs(void) {return VecSolDOFs.data();}

inline unsigned long CFEM_DG_EulerSolver::GetnDOFsGlobal(void) {return nDOFsGlobal;}

inline su2double CFEM_DG_EulerSolver::GetDensity_Inf(void) { return Density_Inf; }

inline su2double CFEM_DG_EulerSolver::GetModVelocity_Inf(void) {
  su2double Vel2 = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  return sqrt(Vel2);
}

inline su2double CFEM_DG_EulerSolver::GetDensity_Energy_Inf(void) { return Density_Inf*Energy_Inf; }

inline su2double CFEM_DG_EulerSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return Density_Inf*Velocity_Inf[val_dim]; }

inline su2double CFEM_DG_EulerSolver::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline su2double *CFEM_DG_EulerSolver::GetVelocity_Inf(void) { return Velocity_Inf; }

inline su2double CFEM_DG_EulerSolver::GetPressure_Inf(void) { return Pressure_Inf; }

inline su2double CFEM_DG_EulerSolver::GetCL_Inv(unsigned short val_marker) { return CL_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetCMz_Inv(unsigned short val_marker) { return CMz_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetCD_Inv(unsigned short val_marker) { return CD_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CL(unsigned short val_marker) { return Surface_CL[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CD(unsigned short val_marker) { return Surface_CD[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CSF(unsigned short val_marker) { return Surface_CSF[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CEff(unsigned short val_marker) { return Surface_CEff[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFx(unsigned short val_marker) { return Surface_CFx[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFy(unsigned short val_marker) { return Surface_CFy[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFz(unsigned short val_marker) { return Surface_CFz[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMx(unsigned short val_marker) { return Surface_CMx[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMy(unsigned short val_marker) { return Surface_CMy[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMz(unsigned short val_marker) { return Surface_CMz[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CL_Inv(unsigned short val_marker) { return Surface_CL_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CD_Inv(unsigned short val_marker) { return Surface_CD_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CSF_Inv(unsigned short val_marker) { return Surface_CSF_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return Surface_CEff_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return Surface_CFx_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return Surface_CFy_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return Surface_CFz_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return Surface_CMx_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return Surface_CMy_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return Surface_CMz_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetCSF_Inv(unsigned short val_marker) { return CSF_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetCEff_Inv(unsigned short val_marker) { return CEff_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CL() { return Total_CL; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CD() { return Total_CD; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CMx() { return Total_CMx; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CMy() { return Total_CMy; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CMz() { return Total_CMz; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CFx() { return Total_CFx; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CFy() { return Total_CFy; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CFz() { return Total_CFz; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CSF() { return Total_CSF; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CEff() { return Total_CEff; }

inline void CFEM_DG_EulerSolver::SetTotal_CL(su2double val_Total_CL) { Total_CL = val_Total_CL; }

inline void CFEM_DG_EulerSolver::SetTotal_CD(su2double val_Total_CD) { Total_CD = val_Total_CD; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CL_Inv() { return AllBound_CL_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CD_Inv() { return AllBound_CD_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CSF_Inv() { return AllBound_CSF_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CEff_Inv() { return AllBound_CEff_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CMx_Inv() { return AllBound_CMx_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CMy_Inv() { return AllBound_CMy_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CMz_Inv() { return AllBound_CMz_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CFx_Inv() { return AllBound_CFx_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CFy_Inv() { return AllBound_CFy_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CFz_Inv() { return AllBound_CFz_Inv; }

inline void CFEM_DG_EulerSolver::SetPressure_Inf(su2double p_inf){Pressure_Inf = p_inf;}

inline void CFEM_DG_EulerSolver::SetTemperature_Inf(su2double t_inf){Temperature_Inf = t_inf;}

inline void CFEM_DG_EulerSolver::BC_HeatFlux_Wall(CConfig                  *config,
                                                  const unsigned long      surfElemBeg,
                                                  const unsigned long      surfElemEnd,
                                                  const CSurfaceElementFEM *surfElem,
                                                  su2double                *resFaces,
                                                  CNumerics                *conv_numerics,
                                                  unsigned short           val_marker,
                                                  su2double                *workArray) {}

inline void CFEM_DG_EulerSolver::BC_Isothermal_Wall(CConfig                  *config,
                                                    const unsigned long      surfElemBeg,
                                                    const unsigned long      surfElemEnd,
                                                    const CSurfaceElementFEM *surfElem,
                                                    su2double                *resFaces,
                                                    CNumerics                *conv_numerics,
                                                    unsigned short           val_marker,
                                                    su2double                *workArray) {}

