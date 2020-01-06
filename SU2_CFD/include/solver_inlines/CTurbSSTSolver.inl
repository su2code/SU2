inline void CTurbSSTSolver::SetFreeStream_Solution(CConfig *config){
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
    nodes->SetSolution(iPoint, 0, kine_Inf);
    nodes->SetSolution(iPoint, 1, omega_Inf);
  }
}

inline su2double CTurbSSTSolver::GetTke_Inf(void) { return kine_Inf; }

inline su2double CTurbSSTSolver::GetOmega_Inf(void) { return omega_Inf; }

