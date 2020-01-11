inline su2double CFEASolver::Get_MassMatrix(unsigned long iPoint, unsigned long jPoint, unsigned short iVar, unsigned short jVar){
  return MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar); }

inline unsigned short CFEASolver::Get_iElem_iDe(unsigned long iElem){ return iElem_iDe[iElem]; }

inline su2double CFEASolver::GetRes_FEM(unsigned short val_var) { return Conv_Check[val_var]; }

inline su2double CFEASolver::GetTotal_CFEA() { return Total_CFEA; }

inline void CFEASolver::SetTotal_CFEA(su2double cfea) { Total_CFEA = cfea; }

inline void CFEASolver::SetTotal_OFRefGeom(su2double val_ofrefgeom) { Total_OFRefGeom = val_ofrefgeom; }

inline void CFEASolver::SetTotal_OFRefNode(su2double val_ofrefnode) { Total_OFRefNode = val_ofrefnode; }

inline su2double CFEASolver::GetWAitken_Dyn(void) { return WAitken_Dyn; }

inline su2double CFEASolver::GetWAitken_Dyn_tn1(void) { return WAitken_Dyn_tn1; }

inline void CFEASolver::SetWAitken_Dyn(su2double waitk) { WAitken_Dyn = waitk; }

inline void CFEASolver::SetWAitken_Dyn_tn1(su2double waitk_tn1) { WAitken_Dyn_tn1 = waitk_tn1; }

inline void CFEASolver::SetLoad_Increment(su2double val_loadIncrement) { loadIncrement = val_loadIncrement; }

inline su2double CFEASolver::GetLoad_Increment(void) { return loadIncrement; }

inline void CFEASolver::SetFSI_ConvValue(unsigned short val_index, su2double val_criteria) { FSI_Conv[val_index] = val_criteria; }

inline su2double CFEASolver::GetFSI_ConvValue(unsigned short val_index){ return FSI_Conv[val_index]; }

inline su2double CFEASolver::GetTotal_OFRefGeom(void){ return Total_OFRefGeom; }

inline su2double CFEASolver::GetTotal_OFRefNode(void){ return Total_OFRefNode; }

inline su2double CFEASolver::GetTotal_OFVolFrac(void){ return Total_OFVolFrac; }

inline su2double CFEASolver::GetTotal_OFCompliance(void){ return Total_OFCompliance; }

inline bool CFEASolver::IsElementBased(void){ return element_based; }

inline void CFEASolver::SetForceCoeff(su2double val_forcecoeff_history) { ForceCoeff = val_forcecoeff_history; }

inline void CFEASolver::SetRelaxCoeff(su2double val_relaxcoeff_history) { RelaxCoeff = val_relaxcoeff_history; }

inline void CFEASolver::SetFSI_Residual(su2double val_FSI_residual) { FSI_Residual = val_FSI_residual; }

inline su2double CFEASolver::GetForceCoeff(void) { return ForceCoeff; }

inline su2double CFEASolver::GetRelaxCoeff(void) { return RelaxCoeff; }

inline su2double CFEASolver::GetFSI_Residual(void) { return FSI_Residual; }

inline su2double CFEASolver::Get_ValCoord(CGeometry *geometry, unsigned long indexNode, unsigned short iDim) {return geometry->node[indexNode]->GetCoord(iDim);}

