inline void CEulerSolver::Set_NewSolution(CGeometry *geometry) { nodes->SetSolution_New(); }

inline su2double CEulerSolver::GetDensity_Inf(void) { return Density_Inf; }

inline su2double CEulerSolver::GetModVelocity_Inf(void) {
  su2double Vel2 = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  return sqrt(Vel2);
}

inline su2double CEulerSolver::GetDensity_Energy_Inf(void) { return Density_Inf*Energy_Inf; }

inline su2double CEulerSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return Density_Inf*Velocity_Inf[val_dim]; }

inline su2double CEulerSolver::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline su2double *CEulerSolver::GetVelocity_Inf(void) { return Velocity_Inf; }

inline su2double CEulerSolver::GetPressure_Inf(void) { return Pressure_Inf; }

inline su2double CEulerSolver::GetCPressure(unsigned short val_marker, unsigned long val_vertex) { return CPressure[val_marker][val_vertex]; }

inline su2double CEulerSolver::GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex) { return CPressureTarget[val_marker][val_vertex]; }

inline void CEulerSolver::SetCPressureTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_pressure) { CPressureTarget[val_marker][val_vertex] = val_pressure; }

inline su2double *CEulerSolver::GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex) { return CharacPrimVar[val_marker][val_vertex]; }

inline void CEulerSolver::SetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { CharacPrimVar[val_marker][val_vertex][val_var] = val_value; }

inline su2double *CEulerSolver::GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex) { return DonorPrimVar[val_marker][val_vertex]; }

inline void CEulerSolver::SetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { DonorPrimVar[val_marker][val_vertex][val_var] = val_value; }

inline su2double CEulerSolver::GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var) { return DonorPrimVar[val_marker][val_vertex][val_var]; }

inline unsigned long CEulerSolver::GetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex) { return DonorGlobalIndex[val_marker][val_vertex]; }

inline void CEulerSolver::SetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex, unsigned long val_index) { DonorGlobalIndex[val_marker][val_vertex] = val_index; }

inline su2double CEulerSolver::GetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex) { return ActDisk_DeltaP[val_marker][val_vertex]; }

inline void CEulerSolver::SetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex, su2double val_deltap) { ActDisk_DeltaP[val_marker][val_vertex] = val_deltap; }

inline su2double CEulerSolver::GetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex) { return ActDisk_DeltaT[val_marker][val_vertex]; }

inline void CEulerSolver::SetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex, su2double val_deltat) { ActDisk_DeltaT[val_marker][val_vertex] = val_deltat; }

inline su2double CEulerSolver::GetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex) { return Inlet_Ttotal[val_marker][val_vertex]; }

inline su2double CEulerSolver::GetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex) { return Inlet_Ptotal[val_marker][val_vertex]; }

inline su2double CEulerSolver::GetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim) { return Inlet_FlowDir[val_marker][val_vertex][val_dim]; }

inline void CEulerSolver::SetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ttotal) {
  /*--- Since this call can be accessed indirectly using python, do some error
   * checking to prevent segmentation faults ---*/
  if (val_marker >= nMarker)
    SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
  else if (Inlet_Ttotal == NULL || Inlet_Ttotal[val_marker] == NULL)
    SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
  else if (val_vertex >= nVertex[val_marker])
    SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
  else
    Inlet_Ttotal[val_marker][val_vertex] = val_ttotal;
}

inline void CEulerSolver::SetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ptotal) {
  /*--- Since this call can be accessed indirectly using python, do some error
   * checking to prevent segmentation faults ---*/
  if (val_marker >= nMarker)
    SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
  else if (Inlet_Ptotal == NULL || Inlet_Ptotal[val_marker] == NULL)
    SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
  else if (val_vertex >= nVertex[val_marker])
    SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
  else
    Inlet_Ptotal[val_marker][val_vertex] = val_ptotal;
}

inline void CEulerSolver::SetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_flowdir) {
  /*--- Since this call can be accessed indirectly using python, do some error
   * checking to prevent segmentation faults ---*/
  if (val_marker >= nMarker)
    SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
  else if (Inlet_FlowDir == NULL || Inlet_FlowDir[val_marker] == NULL)
      SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
  else if (val_vertex >= nVertex[val_marker])
    SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
  else
    Inlet_FlowDir[val_marker][val_vertex][val_dim] = val_flowdir;
}

inline su2double CEulerSolver::GetCL_Inv(unsigned short val_marker) { return CL_Inv[val_marker]; }

inline su2double CEulerSolver::GetCD_Inv(unsigned short val_marker) { return CD_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CL(unsigned short val_marker) { return Surface_CL[val_marker]; }

inline su2double CEulerSolver::GetSurface_CD(unsigned short val_marker) { return Surface_CD[val_marker]; }

inline su2double CEulerSolver::GetSurface_CSF(unsigned short val_marker) { return Surface_CSF[val_marker]; }

inline su2double CEulerSolver::GetSurface_CEff(unsigned short val_marker) { return Surface_CEff[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFx(unsigned short val_marker) { return Surface_CFx[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFy(unsigned short val_marker) { return Surface_CFy[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFz(unsigned short val_marker) { return Surface_CFz[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMx(unsigned short val_marker) { return Surface_CMx[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMy(unsigned short val_marker) { return Surface_CMy[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMz(unsigned short val_marker) { return Surface_CMz[val_marker]; }

inline su2double CEulerSolver::GetSurface_CL_Inv(unsigned short val_marker) { return Surface_CL_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CD_Inv(unsigned short val_marker) { return Surface_CD_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CSF_Inv(unsigned short val_marker) { return Surface_CSF_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return Surface_CEff_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return Surface_CFx_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return Surface_CFy_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return Surface_CFz_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return Surface_CMx_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return Surface_CMy_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return Surface_CMz_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CL_Mnt(unsigned short val_marker) { return Surface_CL_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CD_Mnt(unsigned short val_marker) { return Surface_CD_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CSF_Mnt(unsigned short val_marker) { return Surface_CSF_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CEff_Mnt(unsigned short val_marker) { return Surface_CEff_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFx_Mnt(unsigned short val_marker) { return Surface_CFx_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFy_Mnt(unsigned short val_marker) { return Surface_CFy_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFz_Mnt(unsigned short val_marker) { return Surface_CFz_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMx_Mnt(unsigned short val_marker) { return Surface_CMx_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMy_Mnt(unsigned short val_marker) { return Surface_CMy_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMz_Mnt(unsigned short val_marker) { return Surface_CMz_Mnt[val_marker]; }

inline su2double CEulerSolver::GetInflow_MassFlow(unsigned short val_marker) { return Inflow_MassFlow[val_marker]; }

inline su2double CEulerSolver::GetExhaust_MassFlow(unsigned short val_marker) { return Exhaust_MassFlow[val_marker]; }

inline su2double CEulerSolver::GetInflow_Pressure(unsigned short val_marker) { return Inflow_Pressure[val_marker]; }

inline su2double CEulerSolver::GetInflow_Mach(unsigned short val_marker) { return Inflow_Mach[val_marker]; }

inline su2double CEulerSolver::GetCSF_Inv(unsigned short val_marker) { return CSF_Inv[val_marker]; }

inline su2double CEulerSolver::GetCEff_Inv(unsigned short val_marker) { return CEff_Inv[val_marker]; }

inline su2double CEulerSolver::GetTotal_CL() { return Total_CL; }

inline void CEulerSolver::SetTotal_ComboObj(su2double ComboObj) {Total_ComboObj = ComboObj; }

inline su2double CEulerSolver::GetTotal_ComboObj() { return Total_ComboObj; }

inline su2double CEulerSolver::GetTotal_CD() { return Total_CD; }

inline su2double CEulerSolver::GetTotal_NetThrust() { return Total_NetThrust; }

inline su2double CEulerSolver::GetTotal_Power() { return Total_Power; }

inline su2double CEulerSolver::GetTotal_SolidCD() { return Total_SolidCD; }

inline su2double CEulerSolver::GetTotal_ReverseFlow() { return Total_ReverseFlow; }

inline su2double CEulerSolver::GetTotal_MFR() { return Total_MFR; }

inline su2double CEulerSolver::GetTotal_Prop_Eff() { return Total_Prop_Eff; }

inline su2double CEulerSolver::GetTotal_ByPassProp_Eff() { return Total_ByPassProp_Eff; }

inline su2double CEulerSolver::GetTotal_Adiab_Eff() { return Total_Adiab_Eff; }

inline su2double CEulerSolver::GetTotal_Poly_Eff() { return Total_Poly_Eff; }

inline su2double CEulerSolver::GetTotal_IDC_Mach() { return Total_IDC_Mach; }

inline su2double CEulerSolver::GetTotal_DC60() { return Total_DC60; }

inline su2double CEulerSolver::GetTotal_Custom_ObjFunc() { return Total_Custom_ObjFunc; }

inline su2double CEulerSolver::GetTotal_CMx() { return Total_CMx; }

inline su2double CEulerSolver::GetTotal_CMy() { return Total_CMy; }

inline su2double CEulerSolver::GetTotal_CMz() { return Total_CMz; }

inline su2double CEulerSolver::GetTotal_CoPx() { return Total_CoPx; }

inline su2double CEulerSolver::GetTotal_CoPy() { return Total_CoPy; }

inline su2double CEulerSolver::GetTotal_CoPz() { return Total_CoPz; }

inline su2double CEulerSolver::GetTotal_CFx() { return Total_CFx; }

inline su2double CEulerSolver::GetTotal_CFy() { return Total_CFy; }

inline su2double CEulerSolver::GetTotal_CFz() { return Total_CFz; }

inline su2double CEulerSolver::GetTotal_CSF() { return Total_CSF; }

inline su2double CEulerSolver::GetTotal_CEff() { return Total_CEff; }

inline su2double CEulerSolver::GetTotal_CT() { return Total_CT; }

inline void CEulerSolver::SetTotal_CT(su2double val_Total_CT) { Total_CT = val_Total_CT; }

inline su2double CEulerSolver::GetTotal_CQ() { return Total_CQ; }

inline su2double CEulerSolver::GetTotal_HeatFlux() { return Total_Heat; }

inline su2double CEulerSolver::GetTotal_MaxHeatFlux() { return Total_MaxHeat; }

inline void CEulerSolver::SetTotal_CQ(su2double val_Total_CQ) { Total_CQ = val_Total_CQ; }

inline void CEulerSolver::SetTotal_HeatFlux(su2double val_Total_Heat) { Total_Heat = val_Total_Heat; }

inline void CEulerSolver::SetTotal_MaxHeatFlux(su2double val_Total_MaxHeat) { Total_MaxHeat = val_Total_MaxHeat; }

inline su2double CEulerSolver::GetTotal_CMerit() { return Total_CMerit; }

inline su2double CEulerSolver::GetTotal_CEquivArea() { return Total_CEquivArea; }

inline su2double CEulerSolver::GetTotal_AeroCD() { return Total_AeroCD; }

inline su2double CEulerSolver::GetTotal_IDR() { return Total_IDR; }

inline su2double CEulerSolver::GetTotal_IDC() { return Total_IDC; }

inline su2double CEulerSolver::GetTotal_CpDiff() { return Total_CpDiff; }

inline su2double CEulerSolver::GetTotal_HeatFluxDiff() { return Total_HeatFluxDiff; }

inline su2double CEulerSolver::GetTotal_CNearFieldOF() { return Total_CNearFieldOF; }

inline void CEulerSolver::AddTotal_ComboObj(su2double val_obj) {Total_ComboObj +=val_obj;}

inline void CEulerSolver::SetTotal_CEquivArea(su2double val_cequivarea) { Total_CEquivArea = val_cequivarea; }

inline void CEulerSolver::SetTotal_AeroCD(su2double val_aerocd) { Total_AeroCD = val_aerocd; }

inline void CEulerSolver::SetTotal_CpDiff(su2double pressure) { Total_CpDiff = pressure; }

inline void CEulerSolver::SetTotal_HeatFluxDiff(su2double heat) { Total_HeatFluxDiff = heat; }

inline void CEulerSolver::SetTotal_CNearFieldOF(su2double cnearfieldpress) { Total_CNearFieldOF = cnearfieldpress; }

inline void CEulerSolver::SetTotal_CL(su2double val_Total_CL) { Total_CL = val_Total_CL; }

inline void CEulerSolver::SetTotal_CD(su2double val_Total_CD) { Total_CD = val_Total_CD; }

inline void CEulerSolver::SetTotal_NetThrust(su2double val_Total_NetThrust) { Total_NetThrust = val_Total_NetThrust; }

inline void CEulerSolver::SetTotal_Power(su2double val_Total_Power) { Total_Power = val_Total_Power; }

inline void CEulerSolver::SetTotal_SolidCD(su2double val_Total_SolidCD) { Total_SolidCD = val_Total_SolidCD; }

inline void CEulerSolver::SetTotal_ReverseFlow(su2double val_Total_ReverseFlow) { Total_ReverseFlow = val_Total_ReverseFlow; }

inline void CEulerSolver::SetTotal_MFR(su2double val_Total_MFR) { Total_MFR = val_Total_MFR; }

inline void CEulerSolver::SetTotal_Prop_Eff(su2double val_Total_Prop_Eff) { Total_Prop_Eff = val_Total_Prop_Eff; }

inline void CEulerSolver::SetTotal_ByPassProp_Eff(su2double val_Total_ByPassProp_Eff) { Total_ByPassProp_Eff = val_Total_ByPassProp_Eff; }

inline void CEulerSolver::SetTotal_Adiab_Eff(su2double val_Total_Adiab_Eff) { Total_Adiab_Eff = val_Total_Adiab_Eff; }

inline void CEulerSolver::SetTotal_Poly_Eff(su2double val_Total_Poly_Eff) { Total_Poly_Eff = val_Total_Poly_Eff; }

inline void CEulerSolver::SetTotal_IDC(su2double val_Total_IDC) { Total_IDC = val_Total_IDC; }

inline void CEulerSolver::SetTotal_IDC_Mach(su2double val_Total_IDC_Mach) { Total_IDC_Mach = val_Total_IDC_Mach; }

inline void CEulerSolver::SetTotal_IDR(su2double val_Total_IDR) { Total_IDR = val_Total_IDR; }

inline void CEulerSolver::SetTotal_DC60(su2double val_Total_DC60) { Total_DC60 = val_Total_DC60; }

inline void CEulerSolver::SetTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { Total_Custom_ObjFunc = val_total_custom_objfunc*val_weight; }

inline void CEulerSolver::AddTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { Total_Custom_ObjFunc += val_total_custom_objfunc*val_weight; }

inline su2double CEulerSolver::GetAllBound_CL_Inv() { return AllBound_CL_Inv; }

inline su2double CEulerSolver::GetAllBound_CD_Inv() { return AllBound_CD_Inv; }

inline su2double CEulerSolver::GetAllBound_CSF_Inv() { return AllBound_CSF_Inv; }

inline su2double CEulerSolver::GetAllBound_CEff_Inv() { return AllBound_CEff_Inv; }

inline su2double CEulerSolver::GetAllBound_CMx_Inv() { return AllBound_CMx_Inv; }

inline su2double CEulerSolver::GetAllBound_CMy_Inv() { return AllBound_CMy_Inv; }

inline su2double CEulerSolver::GetAllBound_CMz_Inv() { return AllBound_CMz_Inv; }

inline su2double CEulerSolver::GetAllBound_CoPx_Inv() { return AllBound_CoPx_Inv; }

inline su2double CEulerSolver::GetAllBound_CoPy_Inv() { return AllBound_CoPy_Inv; }

inline su2double CEulerSolver::GetAllBound_CoPz_Inv() { return AllBound_CoPz_Inv; }

inline su2double CEulerSolver::GetAllBound_CFx_Inv() { return AllBound_CFx_Inv; }

inline su2double CEulerSolver::GetAllBound_CFy_Inv() { return AllBound_CFy_Inv; }

inline su2double CEulerSolver::GetAllBound_CFz_Inv() { return AllBound_CFz_Inv; }

inline su2double CEulerSolver::GetAllBound_CL_Mnt() { return AllBound_CL_Mnt; }

inline su2double CEulerSolver::GetAllBound_CD_Mnt() { return AllBound_CD_Mnt; }

inline su2double CEulerSolver::GetAllBound_CSF_Mnt() { return AllBound_CSF_Mnt; }

inline su2double CEulerSolver::GetAllBound_CEff_Mnt() { return AllBound_CEff_Mnt; }

inline su2double CEulerSolver::GetAllBound_CMx_Mnt() { return AllBound_CMx_Mnt; }

inline su2double CEulerSolver::GetAllBound_CMy_Mnt() { return AllBound_CMy_Mnt; }

inline su2double CEulerSolver::GetAllBound_CMz_Mnt() { return AllBound_CMz_Mnt; }

inline su2double CEulerSolver::GetAllBound_CoPx_Mnt() { return AllBound_CoPx_Mnt; }

inline su2double CEulerSolver::GetAllBound_CoPy_Mnt() { return AllBound_CoPy_Mnt; }

inline su2double CEulerSolver::GetAllBound_CoPz_Mnt() { return AllBound_CoPz_Mnt; }

inline su2double CEulerSolver::GetAllBound_CFx_Mnt() { return AllBound_CFx_Mnt; }

inline su2double CEulerSolver::GetAllBound_CFy_Mnt() { return AllBound_CFy_Mnt; }

inline su2double CEulerSolver::GetAllBound_CFz_Mnt() { return AllBound_CFz_Mnt; }

inline su2double CEulerSolver::GetAverageDensity(unsigned short valMarker, unsigned short valSpan){return AverageDensity[valMarker][valSpan];}

inline su2double CEulerSolver::GetAveragePressure(unsigned short valMarker, unsigned short valSpan){return AveragePressure[valMarker][valSpan];}

inline su2double* CEulerSolver::GetAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan){return AverageTurboVelocity[valMarker][valSpan];}

inline su2double CEulerSolver::GetAverageNu(unsigned short valMarker, unsigned short valSpan){return AverageNu[valMarker][valSpan];}

inline su2double CEulerSolver::GetAverageKine(unsigned short valMarker, unsigned short valSpan){return AverageKine[valMarker][valSpan];}

inline su2double CEulerSolver::GetAverageOmega(unsigned short valMarker, unsigned short valSpan){return AverageOmega[valMarker][valSpan];}

inline su2double CEulerSolver::GetExtAverageNu(unsigned short valMarker, unsigned short valSpan){return ExtAverageNu[valMarker][valSpan];}

inline su2double CEulerSolver::GetExtAverageKine(unsigned short valMarker, unsigned short valSpan){return ExtAverageKine[valMarker][valSpan];}

inline su2double CEulerSolver::GetExtAverageOmega(unsigned short valMarker, unsigned short valSpan){return ExtAverageOmega[valMarker][valSpan];}

inline void CEulerSolver::SetExtAverageDensity(unsigned short valMarker, unsigned short valSpan, su2double valDensity){ExtAverageDensity[valMarker][valSpan] = valDensity;}

inline void CEulerSolver::SetExtAveragePressure(unsigned short valMarker, unsigned short valSpan, su2double valPressure){ExtAveragePressure[valMarker][valSpan] = valPressure;}

inline void CEulerSolver::SetExtAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan, unsigned short valIndex, su2double valTurboVelocity){ExtAverageTurboVelocity[valMarker][valSpan][valIndex] = valTurboVelocity;}

inline void CEulerSolver::SetExtAverageNu(unsigned short valMarker, unsigned short valSpan, su2double valNu){ExtAverageNu[valMarker][valSpan] = valNu;}

inline void CEulerSolver::SetExtAverageKine(unsigned short valMarker, unsigned short valSpan, su2double valKine){ExtAverageKine[valMarker][valSpan] = valKine;}

inline void CEulerSolver::SetExtAverageOmega(unsigned short valMarker, unsigned short valSpan, su2double valOmega){ExtAverageOmega[valMarker][valSpan] = valOmega;}

inline su2double  CEulerSolver::GetDensityIn(unsigned short inMarkerTP, unsigned short valSpan){return DensityIn[inMarkerTP][valSpan];}

inline su2double  CEulerSolver::GetPressureIn(unsigned short inMarkerTP, unsigned short valSpan){return PressureIn[inMarkerTP][valSpan];}

inline su2double* CEulerSolver::GetTurboVelocityIn(unsigned short inMarkerTP, unsigned short valSpan){return TurboVelocityIn[inMarkerTP][valSpan];}

inline su2double  CEulerSolver::GetDensityOut(unsigned short inMarkerTP, unsigned short valSpan){return DensityOut[inMarkerTP][valSpan];}

inline su2double  CEulerSolver::GetPressureOut(unsigned short inMarkerTP, unsigned short valSpan){return PressureOut[inMarkerTP][valSpan];}

inline su2double* CEulerSolver::GetTurboVelocityOut(unsigned short inMarkerTP, unsigned short valSpan){return TurboVelocityOut[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetKineIn(unsigned short inMarkerTP, unsigned short valSpan){return KineIn[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetOmegaIn(unsigned short inMarkerTP, unsigned short valSpan){return OmegaIn[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetNuIn(unsigned short inMarkerTP, unsigned short valSpan){return NuIn[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetKineOut(unsigned short inMarkerTP, unsigned short valSpan){return KineOut[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetOmegaOut(unsigned short inMarkerTP, unsigned short valSpan){return OmegaOut[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetNuOut(unsigned short inMarkerTP, unsigned short valSpan){return NuOut[inMarkerTP][valSpan];}

inline void CEulerSolver::SetDensityIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){DensityIn[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetPressureIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){PressureIn[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetTurboVelocityIn(su2double *value, unsigned short inMarkerTP, unsigned short valSpan){
  unsigned short iDim;

  for(iDim = 0; iDim < nDim; iDim++)
    TurboVelocityIn[inMarkerTP][valSpan][iDim] = value[iDim];
}

inline void CEulerSolver::SetDensityOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){DensityOut[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetPressureOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){PressureOut[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetTurboVelocityOut(su2double *value, unsigned short inMarkerTP, unsigned short valSpan){
  unsigned short iDim;

  for(iDim = 0; iDim < nDim; iDim++)
    TurboVelocityOut[inMarkerTP][valSpan][iDim] = value[iDim];
}

inline void CEulerSolver::SetKineIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){KineIn[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetOmegaIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){OmegaIn[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetNuIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){NuIn[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetKineOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){KineOut[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetOmegaOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){OmegaOut[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetNuOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){NuOut[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::ComputeTurboVelocity(const su2double *cartesianVelocity, const su2double *turboNormal, su2double *turboVelocity,
                                               unsigned short marker_flag, unsigned short kind_turb) {

  if ((kind_turb == AXIAL && nDim == 3) || (kind_turb == CENTRIPETAL_AXIAL && marker_flag == OUTFLOW) || (kind_turb == AXIAL_CENTRIFUGAL && marker_flag == INFLOW) ){
    turboVelocity[2] =  turboNormal[0]*cartesianVelocity[0] + cartesianVelocity[1]*turboNormal[1];
    turboVelocity[1] =  turboNormal[0]*cartesianVelocity[1] - turboNormal[1]*cartesianVelocity[0];
    turboVelocity[0] = cartesianVelocity[2];
  }
  else{
    turboVelocity[0] =  turboNormal[0]*cartesianVelocity[0] + cartesianVelocity[1]*turboNormal[1];
    turboVelocity[1] =  turboNormal[0]*cartesianVelocity[1] - turboNormal[1]*cartesianVelocity[0];
    if (marker_flag == INFLOW){
      turboVelocity[0] *= -1.0;
      turboVelocity[1] *= -1.0;
    }
    if(nDim == 3)
      turboVelocity[2] = cartesianVelocity[2];
  }
}

inline void CEulerSolver::ComputeBackVelocity(const su2double *turboVelocity, const su2double *turboNormal, su2double *cartesianVelocity,
                                              unsigned short marker_flag, unsigned short kind_turb){

  if ((kind_turb == AXIAL && nDim == 3) || (kind_turb == CENTRIPETAL_AXIAL && marker_flag == OUTFLOW) || (kind_turb == AXIAL_CENTRIFUGAL && marker_flag == INFLOW)){
    cartesianVelocity[0] = turboVelocity[2]*turboNormal[0] - turboVelocity[1]*turboNormal[1];
    cartesianVelocity[1] = turboVelocity[2]*turboNormal[1] + turboVelocity[1]*turboNormal[0];
    cartesianVelocity[2] = turboVelocity[0];
  }
  else{
    cartesianVelocity[0] =  turboVelocity[0]*turboNormal[0] - turboVelocity[1]*turboNormal[1];
    cartesianVelocity[1] =  turboVelocity[0]*turboNormal[1] + turboVelocity[1]*turboNormal[0];

    if (marker_flag == INFLOW){
      cartesianVelocity[0] *= -1.0;
      cartesianVelocity[1] *= -1.0;
    }

    if(nDim == 3)
      cartesianVelocity[2] = turboVelocity[2];
  }
}


inline CFluidModel* CEulerSolver::GetFluidModel(void) { return FluidModel;}

inline void CEulerSolver::SetPressure_Inf(su2double p_inf) {Pressure_Inf = p_inf;}

inline void CEulerSolver::SetTemperature_Inf(su2double t_inf) {Temperature_Inf = t_inf;}

inline bool CEulerSolver::GetStart_AoA_FD(void) { return Start_AoA_FD; }

inline bool CEulerSolver::GetEnd_AoA_FD(void) { return End_AoA_FD; }

inline unsigned long CEulerSolver::GetIter_Update_AoA(void) { return Iter_Update_AoA; }

inline su2double CEulerSolver::GetPrevious_AoA(void) { return AoA_Prev; }

inline su2double CEulerSolver::GetAoA_inc(void) { return AoA_inc; }

inline void CEulerSolver::SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component){
  SlidingState[val_marker][val_vertex][val_state][donor_index] = component;
}

inline su2double CEulerSolver::GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index) { return SlidingState[val_marker][val_vertex][val_state][donor_index]; }

inline int CEulerSolver::GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex){ return SlidingStateNodes[val_marker][val_vertex]; }

inline void CEulerSolver::SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value){ SlidingStateNodes[val_marker][val_vertex] = value; }

inline void CEulerSolver::SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex){
  int iVar;

  for( iVar = 0; iVar < nPrimVar+1; iVar++){
    if( SlidingState[val_marker][val_vertex][iVar] != NULL )
      delete [] SlidingState[val_marker][val_vertex][iVar];
  }

  for( iVar = 0; iVar < nPrimVar+1; iVar++)
    SlidingState[val_marker][val_vertex][iVar] = new su2double[ GetnSlidingStates(val_marker, val_vertex) ];
}


