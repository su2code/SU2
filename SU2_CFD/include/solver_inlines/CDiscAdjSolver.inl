inline su2double CDiscAdjSolver::GetTotal_Sens_Geo() { return Total_Sens_Geo; }

inline su2double CDiscAdjSolver::GetTotal_Sens_Mach() { return Total_Sens_Mach; }

inline su2double CDiscAdjSolver::GetTotal_Sens_AoA() { return Total_Sens_AoA; }

inline su2double CDiscAdjSolver::GetTotal_Sens_Press() { return Total_Sens_Press; }

inline su2double CDiscAdjSolver::GetTotal_Sens_Temp() { return Total_Sens_Temp; }

inline su2double CDiscAdjSolver::GetTotal_Sens_BPress() { return Total_Sens_BPress; }

inline su2double CDiscAdjSolver::GetTotal_Sens_Density() { return Total_Sens_Density; }

inline su2double CDiscAdjSolver::GetTotal_Sens_ModVel() { return Total_Sens_ModVel; }

inline su2double CDiscAdjSolver::GetCSensitivity(unsigned short val_marker, unsigned long val_vertex) { return CSensitivity[val_marker][val_vertex]; }

