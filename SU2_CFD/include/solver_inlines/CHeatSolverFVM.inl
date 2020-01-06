inline su2double CHeatSolverFVM::GetTotal_HeatFlux() { return Total_HeatFlux; }

inline su2double CHeatSolverFVM::GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) { return HeatFlux[val_marker][val_vertex]; }

inline su2double CHeatSolverFVM::GetTotal_AvgTemperature() { return Total_AverageT; }

inline su2double CHeatSolverFVM::GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var) { return ConjugateVar[val_marker][val_vertex][pos_var]; }

inline void CHeatSolverFVM::SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var, su2double relaxation_factor, su2double val_var) {
  ConjugateVar[val_marker][val_vertex][pos_var] = relaxation_factor*val_var + (1.0-relaxation_factor)*ConjugateVar[val_marker][val_vertex][pos_var]; }

