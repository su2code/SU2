inline void CVerificationSolution::SetError_RMS(unsigned short val_var, su2double val_error) { Error_RMS[val_var] = val_error; }

inline void CVerificationSolution::AddError_RMS(unsigned short val_var, su2double val_error) { Error_RMS[val_var] += val_error; }

inline su2double CVerificationSolution::GetError_RMS(unsigned short val_var) { return Error_RMS[val_var]; }

inline void CVerificationSolution::SetError_Max(unsigned short val_var, su2double val_error, unsigned long val_point) {
  Error_Max[val_var]       = val_error;
  Error_Point_Max[val_var] = val_point;
}

inline void CVerificationSolution::AddError_Max(unsigned short val_var, su2double val_error, unsigned long val_point, su2double* val_coord) {
  if (val_error > Error_Max[val_var]) {
    Error_Max[val_var] = val_error;
    Error_Point_Max[val_var] = val_point;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Error_Point_Max_Coord[val_var][iDim] = val_coord[iDim];
  }
}

inline su2double CVerificationSolution::GetError_Max(unsigned short val_var) { return Error_Max[val_var]; }

inline unsigned long CVerificationSolution::GetError_Point_Max(unsigned short val_var) { return Error_Point_Max[val_var]; }

inline su2double* CVerificationSolution::GetError_Point_Max_Coord(unsigned short val_var) { return Error_Point_Max_Coord[val_var]; }
