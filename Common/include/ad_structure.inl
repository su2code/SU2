#pragma once


namespace AD{
#ifdef CODI_REVERSE_TYPE

  typedef codi::DataStore CheckpointHandler;

  /*--- Stores the indices of the input variables (they might be overwritten) ---*/

  extern std::vector<unsigned int> inputValues;

  /*--- Current position inside the adjoint vector ---*/

  extern int adjointVectorPosition;

  /*--- Reference to the tape ---*/

  extern codi::ChunkTape<double, int>& globalTape;

  extern bool Status;

  extern bool PreaccActive;

  extern codi::ChunkTape<double, int>::Position StartPosition, EndPosition;

  extern std::vector<unsigned int> localInputValues;

  extern std::vector<su2double*> localOutputValues;

  inline void RegisterInput(su2double &data){AD::globalTape.registerInput(data);
                                             inputValues.push_back(data.getGradientData());}

  inline void RegisterOutput(su2double& data){AD::globalTape.registerOutput(data);}

  inline void ResetInput(su2double &data){data.getGradientData() = 0;}

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

  inline void SetPreaccInput(const su2double &data){
    if (data.getGradientData() != 0){
      localInputValues.push_back(data.getGradientData());
    }
  }

  inline void SetPreaccInput(const CArray1D& data){
    for (unsigned short i = 0; i < data.size; i++){
      if (data.vec[i].getGradientData() != 0){
        localInputValues.push_back(data.vec[i].getGradientData());
      }
    }
  }

  inline void SetPreaccInput(const CArray2D& data){
    for (unsigned short i = 0; i < data.size_x; i++){
      for (unsigned short j = 0; j < data.size_y; j++){
        if (data.mat[i][j].getGradientData() != 0){
          localInputValues.push_back(data.mat[i][j].getGradientData());
        }
      }
    }
  }

  inline void SetPreaccInput_Variadic(){}

  template <typename Arg1, typename ... Args>
  inline void SetPreaccInput_Variadic(const Arg1& arg1, Args& ... args){
    SetPreaccInput(arg1);
    SetPreaccInput_Variadic(args...);
  }

  template <typename ... Args>
  inline void StartPreacc(Args && ... args){
    if (globalTape.isActive()){
      SetPreaccInput_Variadic(args...);
      StartPosition = globalTape.getPosition();
      PreaccActive = true;
    }
  }

  inline void SetPreaccOutput(su2double& data){
    if (data.getGradientData() != 0){
      localOutputValues.push_back(&data);
    }
  }

  inline void SetPreaccOutput(CArray1D& data){
    for (unsigned short i = 0; i < data.size; i++){
      if (data.vec[i].getGradientData() != 0){
        localOutputValues.push_back(&data.vec[i]);
      }
    }
  }

  inline void SetPreaccOutput(CArray2D& data){
    for (unsigned short i = 0; i < data.size_x; i++){
      for (unsigned short j = 0; j < data.size_y; j++){
        if (data.mat[i][j].getGradientData() != 0){
          localOutputValues.push_back(&data.mat[i][j]);
        }
      }
    }
  }

  inline void SetPreaccOutput_Variadic(){}

  template <typename Arg1, typename ... Args>
  inline void SetPreaccOutput_Variadic(Arg1& arg1, Args& ... args){
    SetPreaccOutput(arg1);
    SetPreaccOutput_Variadic(args...);
  }


  template <typename ... Args>
  inline void EndPreacc(Args && ... args){
    if (PreaccActive){
      SetPreaccOutput_Variadic(args...);
      Preaccumulate();
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

  inline void Preaccumulate(){}

  inline void SetPreaccInput(const su2double &data){}

  inline void SetPreaccInput(const CArray1D &data){}

  inline void SetPreaccInput(const CArray2D &data){}

  inline void SetPreaccInput_Variadic(){}

  template <typename Arg1, typename ... Args>
  inline void SetPreaccInput_Variadic(const Arg1& arg1, Args& ... args){}

  template <typename ... Args>
  inline void StartPreacc(Args && ... args){}

  inline void SetPreaccOutput(su2double &data){}

  inline void SetPreaccOutput(CArray1D &data){}

  inline void SetPreaccOutput(CArray2D &data){}

  inline void SetPreaccOutput_Variadic(){}

  template <typename Arg1, typename ... Args>
  inline void SetPreaccOutput_Variadic(Arg1& arg1, Args& ... args){}

  template <typename ... Args>
  inline void EndPreacc(Args && ... args){}

#endif
}
