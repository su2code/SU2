inline su2double CAdjEulerSolver::GetCSensitivity(unsigned short val_marker, unsigned long val_vertex) { return CSensitivity[val_marker][val_vertex]; }

inline void CAdjEulerSolver::SetCSensitivity(unsigned short val_marker, unsigned long val_vertex, su2double val_sensitivity) { CSensitivity[val_marker][val_vertex] = val_sensitivity; }

inline unsigned long CAdjEulerSolver::GetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex) { return DonorGlobalIndex[val_marker][val_vertex]; }

inline void CAdjEulerSolver::SetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex, unsigned long val_index) { DonorGlobalIndex[val_marker][val_vertex] = val_index; }

inline void CAdjEulerSolver::SetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { DonorAdjVar[val_marker][val_vertex][val_var] = val_value; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Geo() { return Total_Sens_Geo; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Mach() { return Total_Sens_Mach; }

inline su2double CAdjEulerSolver::GetTotal_Sens_AoA() { return Total_Sens_AoA; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Press() { return Total_Sens_Press; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Temp() { return Total_Sens_Temp; }

inline su2double CAdjEulerSolver::GetTotal_Sens_BPress() { return Total_Sens_BPress; }

inline su2double CAdjEulerSolver::GetPsiRho_Inf(void) { return PsiRho_Inf; }

inline su2double CAdjEulerSolver::GetPsiE_Inf(void) { return PsiE_Inf; }

inline su2double *CAdjEulerSolver::GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex) { return DonorAdjVar[val_marker][val_vertex]; }

inline su2double CAdjEulerSolver::GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var) { return DonorAdjVar[val_marker][val_vertex][val_var]; }

inline su2double CAdjEulerSolver::GetPhi_Inf(unsigned short val_dim) { return Phi_Inf[val_dim]; }

