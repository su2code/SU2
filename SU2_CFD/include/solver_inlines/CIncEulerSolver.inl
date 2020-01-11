inline su2double CIncEulerSolver::GetDensity_Inf(void) { return Density_Inf; }

inline su2double CIncEulerSolver::GetModVelocity_Inf(void) {
  su2double Vel2 = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  return sqrt(Vel2);
}

inline CFluidModel* CIncEulerSolver::GetFluidModel(void) { return FluidModel;}

inline su2double CIncEulerSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return Density_Inf*Velocity_Inf[val_dim]; }

inline su2double CIncEulerSolver::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline void CIncEulerSolver::SetVelocity_Inf(unsigned short val_dim, su2double val_velocity) { Velocity_Inf[val_dim] = val_velocity; }

inline su2double *CIncEulerSolver::GetVelocity_Inf(void) { return Velocity_Inf; }

inline su2double CIncEulerSolver::GetPressure_Inf(void) { return Pressure_Inf; }

inline su2double CIncEulerSolver::GetTemperature_Inf(void) { return Temperature_Inf; }

inline su2double CIncEulerSolver::GetCPressure(unsigned short val_marker, unsigned long val_vertex) { return CPressure[val_marker][val_vertex]; }

inline su2double CIncEulerSolver::GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex) { return CPressureTarget[val_marker][val_vertex]; }

inline void CIncEulerSolver::SetCPressureTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_pressure) { CPressureTarget[val_marker][val_vertex] = val_pressure; }

inline su2double *CIncEulerSolver::GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex) { return CharacPrimVar[val_marker][val_vertex]; }

inline su2double CIncEulerSolver::GetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex) { return Inlet_Ttotal[val_marker][val_vertex]; }

inline su2double CIncEulerSolver::GetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex) { return Inlet_Ptotal[val_marker][val_vertex]; }

inline su2double CIncEulerSolver::GetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim) { return Inlet_FlowDir[val_marker][val_vertex][val_dim]; }

inline su2double CIncEulerSolver::GetCD_Inv(unsigned short val_marker) { return CD_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CL(unsigned short val_marker) { return Surface_CL[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CD(unsigned short val_marker) { return Surface_CD[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CSF(unsigned short val_marker) { return Surface_CSF[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CEff(unsigned short val_marker) { return Surface_CEff[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFx(unsigned short val_marker) { return Surface_CFx[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFy(unsigned short val_marker) { return Surface_CFy[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFz(unsigned short val_marker) { return Surface_CFz[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMx(unsigned short val_marker) { return Surface_CMx[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMy(unsigned short val_marker) { return Surface_CMy[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMz(unsigned short val_marker) { return Surface_CMz[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CL_Inv(unsigned short val_marker) { return Surface_CL_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CD_Inv(unsigned short val_marker) { return Surface_CD_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CSF_Inv(unsigned short val_marker) { return Surface_CSF_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return Surface_CEff_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return Surface_CFx_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return Surface_CFy_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return Surface_CFz_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return Surface_CMx_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return Surface_CMy_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return Surface_CMz_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetCSF_Inv(unsigned short val_marker) { return CSF_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetCEff_Inv(unsigned short val_marker) { return CEff_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetTotal_CL() { return Total_CL; }

inline su2double CIncEulerSolver::GetTotal_CD() { return Total_CD; }

inline su2double CIncEulerSolver::GetTotal_CMx() { return Total_CMx; }

inline su2double CIncEulerSolver::GetTotal_CMy() { return Total_CMy; }

inline su2double CIncEulerSolver::GetTotal_CMz() { return Total_CMz; }

inline su2double CIncEulerSolver::GetTotal_CoPx() { return Total_CoPx; }

inline su2double CIncEulerSolver::GetTotal_CoPy() { return Total_CoPy; }

inline su2double CIncEulerSolver::GetTotal_CoPz() { return Total_CoPz; }

inline su2double CIncEulerSolver::GetTotal_CFx() { return Total_CFx; }

inline su2double CIncEulerSolver::GetTotal_CFy() { return Total_CFy; }

inline su2double CIncEulerSolver::GetTotal_CFz() { return Total_CFz; }

inline su2double CIncEulerSolver::GetTotal_CSF() { return Total_CSF; }

inline su2double CIncEulerSolver::GetTotal_CEff() { return Total_CEff; }

inline su2double CIncEulerSolver::GetTotal_CT() { return Total_CT; }

inline void CIncEulerSolver::SetTotal_CT(su2double val_Total_CT) { Total_CT = val_Total_CT; }

inline su2double CIncEulerSolver::GetTotal_CQ() { return Total_CQ; }

inline su2double CIncEulerSolver::GetTotal_HeatFlux() { return Total_Heat; }

inline su2double CIncEulerSolver::GetTotal_MaxHeatFlux() { return Total_MaxHeat; }

inline void CIncEulerSolver::SetTotal_CQ(su2double val_Total_CQ) { Total_CQ = val_Total_CQ; }

inline void CIncEulerSolver::SetTotal_HeatFlux(su2double val_Total_Heat) { Total_Heat = val_Total_Heat; }

inline void CIncEulerSolver::SetTotal_MaxHeatFlux(su2double val_Total_MaxHeat) { Total_MaxHeat = val_Total_MaxHeat; }

inline su2double CIncEulerSolver::GetTotal_CMerit() { return Total_CMerit; }

inline su2double CIncEulerSolver::GetTotal_CpDiff() { return Total_CpDiff; }

inline su2double CIncEulerSolver::GetTotal_HeatFluxDiff() { return Total_HeatFluxDiff; }

inline void CIncEulerSolver::SetTotal_CpDiff(su2double val_pressure) { Total_CpDiff = pressure; }

inline void CIncEulerSolver::SetTotal_HeatFluxDiff(su2double val_heat) { Total_HeatFluxDiff = heat; }

inline void CIncEulerSolver::SetTotal_CD(su2double val_Total_CD) { Total_CD = val_Total_CD; }

inline su2double CIncEulerSolver::GetTotal_Custom_ObjFunc() { return Total_Custom_ObjFunc; }

inline void CIncEulerSolver::SetTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { Total_Custom_ObjFunc = val_total_custom_objfunc*val_weight; }

inline void CIncEulerSolver::AddTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { Total_Custom_ObjFunc += val_total_custom_objfunc*val_weight; }

inline su2double CIncEulerSolver::GetAllBound_CL_Inv() { return AllBound_CL_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CD_Inv() { return AllBound_CD_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CSF_Inv() { return AllBound_CSF_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CEff_Inv() { return AllBound_CEff_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CMx_Inv() { return AllBound_CMx_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CMy_Inv() { return AllBound_CMy_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CMz_Inv() { return AllBound_CMz_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CoPx_Inv() { return AllBound_CoPx_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CoPy_Inv() { return AllBound_CoPy_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CoPz_Inv() { return AllBound_CoPz_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CFx_Inv() { return AllBound_CFx_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CFy_Inv() { return AllBound_CFy_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CFz_Inv() { return AllBound_CFz_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CL_Mnt() { return AllBound_CL_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CD_Mnt() { return AllBound_CD_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CSF_Mnt() { return AllBound_CSF_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CEff_Mnt() { return AllBound_CEff_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CMx_Mnt() { return AllBound_CMx_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CMy_Mnt() { return AllBound_CMy_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CMz_Mnt() { return AllBound_CMz_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CoPx_Mnt() { return AllBound_CoPx_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CoPy_Mnt() { return AllBound_CoPy_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CoPz_Mnt() { return AllBound_CoPz_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CFx_Mnt() { return AllBound_CFx_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CFy_Mnt() { return AllBound_CFy_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CFz_Mnt() { return AllBound_CFz_Mnt; }

inline su2double CIncEulerSolver::GetSurface_CL_Mnt(unsigned short val_marker) { return Surface_CL_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CD_Mnt(unsigned short val_marker) { return Surface_CD_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CSF_Mnt(unsigned short val_marker) { return Surface_CSF_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CEff_Mnt(unsigned short val_marker) { return Surface_CEff_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFx_Mnt(unsigned short val_marker) { return Surface_CFx_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFy_Mnt(unsigned short val_marker) { return Surface_CFy_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFz_Mnt(unsigned short val_marker) { return Surface_CFz_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMx_Mnt(unsigned short val_marker) { return Surface_CMx_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMy_Mnt(unsigned short val_marker) { return Surface_CMy_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMz_Mnt(unsigned short val_marker) { return Surface_CMz_Mnt[val_marker]; }

inline void CIncEulerSolver::SetPressure_Inf(su2double p_inf){Pressure_Inf = p_inf;}

inline void CIncEulerSolver::SetTemperature_Inf(su2double t_inf){Temperature_Inf = t_inf;}

inline void CIncEulerSolver::SetDensity_Inf(su2double rho_inf){Density_Inf = rho_inf;}

inline void CIncEulerSolver::SetTotal_ComboObj(su2double ComboObj) {Total_ComboObj = ComboObj; }

inline su2double CIncEulerSolver::GetTotal_ComboObj() { return Total_ComboObj; }

inline void CIncEulerSolver::SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component){
  SlidingState[val_marker][val_vertex][val_state][donor_index] = component;
}

inline su2double CIncEulerSolver::GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index) { return SlidingState[val_marker][val_vertex][val_state][donor_index]; }

inline int CIncEulerSolver::GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex){ return SlidingStateNodes[val_marker][val_vertex]; }

inline void CIncEulerSolver::SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value){ SlidingStateNodes[val_marker][val_vertex] = value; }

inline void CIncEulerSolver::SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex){
  int iVar;

  for( iVar = 0; iVar < nPrimVar+1; iVar++){
    if( SlidingState[val_marker][val_vertex][iVar] != NULL )
      delete [] SlidingState[val_marker][val_vertex][iVar];
  }

  for( iVar = 0; iVar < nPrimVar+1; iVar++)
    SlidingState[val_marker][val_vertex][iVar] = new su2double[ GetnSlidingStates(val_marker, val_vertex) ];
}



