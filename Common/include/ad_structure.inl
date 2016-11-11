#pragma once


namespace AD{
#if defined CODI_REVERSE_TYPE

  typedef codi::DataStore CheckpointHandler;

  /*--- Stores the indices of the input variables (they might be overwritten) ---*/

  extern std::vector<su2double::GradientData> inputValues;

  /*--- Current position inside the adjoint vector ---*/

  extern int adjointVectorPosition;

  /*--- Reference to the tape ---*/

  extern su2double::TapeType& globalTape;

  extern bool Status;

  extern bool PreaccActive;

  extern su2double::TapeType::Position StartPosition, EndPosition;

  extern std::vector<su2double::GradientData> localInputValues;

  extern std::vector<su2double*> localOutputValues;

  inline void RegisterInput(su2double &data){AD::globalTape.registerInput(data);
                                             inputValues.push_back(data.getGradientData());}

  inline void RegisterOutput(su2double& data){AD::globalTape.registerOutput(data);}

  inline void ResetInput(su2double &data){data.getGradientData() = su2double::GradientData();}

  inline void StartRecording(){AD::globalTape.setActive();}

  inline void StopRecording(){AD::globalTape.setPassive();}

  inline void ClearAdjoints(){AD::globalTape.clearAdjoints(); }

  inline void ComputeAdjoint(){AD::globalTape.evaluate();
                               adjointVectorPosition = 0;}

  inline void Reset(){
    if (inputValues.size() != 0){
      globalTape.reset();
      adjointVectorPosition = 0;
      inputValues.clear();
    }
  }

  inline void SetPreaccIn(const su2double &data){
    if (PreaccActive){
      if (data.isActive()){
        localInputValues.push_back(data.getGradientData());
      }
    }
  }

  inline void SetPreaccIn(su2double* data, const int size){
    if (PreaccActive){
      for (unsigned short i = 0; i < size; i++){
        if (data[i].isActive()){
          localInputValues.push_back(data[i].getGradientData());
        }
      }
    }
  }

  inline void SetPreaccIn(su2double** data, const int size_x, const int size_y){
    if (PreaccActive){
      for (unsigned short i = 0; i < size_x; i++){
        for (unsigned short j = 0; j < size_y; j++){
          if (data[i][j].isActive()){
            localInputValues.push_back(data[i][j].getGradientData());
          }
        }
      }
    }
  }

  inline void StartPreacc(){
    if (globalTape.isActive()){
      StartPosition = globalTape.getPosition();
      PreaccActive = true;
    }
  }

  inline void SetPreaccOut(su2double& data){
    if (PreaccActive){
      if (data.isActive()){
        localOutputValues.push_back(&data);
      }
    }
  }

  inline void SetPreaccOut(su2double* data, const int size){
    if (PreaccActive){
      for (unsigned short i = 0; i < size; i++){
        if (data[i].isActive()){
          localOutputValues.push_back(&data[i]);
        }
      }
    }
  }

  inline void SetPreaccOut(su2double** data, const int size_x, const int size_y){
    if (PreaccActive){
      for (unsigned short i = 0; i < size_x; i++){
        for (unsigned short j = 0; j < size_y; j++){
          if (data[i][j].isActive()){
            localOutputValues.push_back(&data[i][j]);
          }
        }
      }
    }
  }


  inline void delete_handler(void *handler){
    CheckpointHandler *checkpoint = static_cast<CheckpointHandler*>(handler);
    checkpoint->clear();
  }
#else

  /*--- Default implementation if reverse mode is disabled ---*/

  inline void RegisterInput(su2double &data){}

  inline void RegisterOutput(su2double& data){}

  inline void StartRecording(){}

  inline void StopRecording(){}

  inline void ClearAdjoints(){}

  inline void ComputeAdjoint(){}

  inline void Reset(){}

  inline void ResetInput(su2double &data){}

  inline void SetPreaccIn(const su2double &data){}

  inline void SetPreaccIn(su2double* data, const int size){}

  inline void SetPreaccIn(su2double** data, const int size_x, const int size_y){}

  inline void SetPreaccOut(su2double &data){}

  inline void SetPreaccOut(su2double* data, const int size){}

  inline void SetPreaccOut(su2double** data, const int size_x, const int size_y){}

  inline void StartPreacc(){}

  inline void EndPreacc(){}
#endif
}
