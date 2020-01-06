inline void CTurbSolver::SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component){
  SlidingState[val_marker][val_vertex][val_state][donor_index] = component;
}

inline int CTurbSolver::GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex){ return SlidingStateNodes[val_marker][val_vertex]; }

inline void CTurbSolver::SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex){
  int iVar;

  for( iVar = 0; iVar < nVar+1; iVar++){
    if( SlidingState[val_marker][val_vertex][iVar] != NULL )
      delete [] SlidingState[val_marker][val_vertex][iVar];
  }

  for( iVar = 0; iVar < nVar+1; iVar++)
    SlidingState[val_marker][val_vertex][iVar] = new su2double[ GetnSlidingStates(val_marker, val_vertex) ];
}

inline void CTurbSolver::SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value){ SlidingStateNodes[val_marker][val_vertex] = value; }

inline su2double CTurbSolver::GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index) { return SlidingState[val_marker][val_vertex][val_state][donor_index]; }

inline void CTurbSolver::SetInlet_TurbVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_turb_var) {
  /*--- Since this call can be accessed indirectly using python, do some error
   * checking to prevent segmentation faults ---*/
  if (val_marker >= nMarker)
    SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
  else if (Inlet_TurbVars == NULL || Inlet_TurbVars[val_marker] == NULL)
    SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
  else if (val_vertex >= nVertex[val_marker])
    SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
  else if (val_dim >= nVar)
    SU2_MPI::Error("Out-of-bounds index used for inlet turbulence variable.", CURRENT_FUNCTION);
  else
    Inlet_TurbVars[val_marker][val_vertex][val_dim] = val_turb_var;
}

