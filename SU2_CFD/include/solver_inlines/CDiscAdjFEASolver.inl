inline su2double CDiscAdjFEASolver::GetTotal_Sens_E(unsigned short iVal) { return Total_Sens_E[iVal]; }

inline su2double CDiscAdjFEASolver::GetTotal_Sens_Nu(unsigned short iVal) { return Total_Sens_Nu[iVal]; }

inline su2double CDiscAdjFEASolver::GetTotal_Sens_Rho(unsigned short iVal) { return Total_Sens_Rho[iVal]; }

inline su2double CDiscAdjFEASolver::GetTotal_Sens_Rho_DL(unsigned short iVal) { return Total_Sens_Rho_DL[iVal]; }

inline su2double CDiscAdjFEASolver::GetTotal_Sens_EField(unsigned short iEField) { return Total_Sens_EField[iEField]; }

inline su2double CDiscAdjFEASolver::GetTotal_Sens_DVFEA(unsigned short iDVFEA) { return Total_Sens_DV[iDVFEA]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_E(unsigned short iVal) { return Global_Sens_E[iVal]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_Nu(unsigned short iVal) { return Global_Sens_Nu[iVal]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_Rho(unsigned short iVal) { return Global_Sens_Rho[iVal]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_Rho_DL(unsigned short iVal) { return Global_Sens_Rho_DL[iVal]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_EField(unsigned short iEField) { return Global_Sens_EField[iEField]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_DVFEA(unsigned short iDVFEA) { return Global_Sens_DV[iDVFEA]; }

inline su2double CDiscAdjFEASolver::GetVal_Young(unsigned short iVal) { return E_i[iVal]; }

inline su2double CDiscAdjFEASolver::GetVal_Poisson(unsigned short iVal) { return Nu_i[iVal]; }

inline su2double CDiscAdjFEASolver::GetVal_Rho(unsigned short iVal) { return Rho_i[iVal]; }

inline su2double CDiscAdjFEASolver::GetVal_Rho_DL(unsigned short iVal) { return Rho_DL_i[iVal]; }

inline unsigned short CDiscAdjFEASolver::GetnEField(void) { return nEField; }

inline unsigned short CDiscAdjFEASolver::GetnDVFEA(void) { return nDV; }

inline su2double CDiscAdjFEASolver::GetVal_EField(unsigned short iVal) { return EField[iVal]; }

inline su2double CDiscAdjFEASolver::GetVal_DVFEA(unsigned short iVal) { return DV_Val[iVal]; }

