inline void CTurbSASolver::SetFreeStream_Solution(CConfig *config) {
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) nodes->SetSolution(iPoint, 0, nu_tilde_Inf);
}

inline su2double CTurbSASolver::GetNuTilde_Inf(void) { return nu_tilde_Inf; }

