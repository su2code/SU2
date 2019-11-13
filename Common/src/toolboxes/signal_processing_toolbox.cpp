#include "../../include/toolboxes/signal_processing_toolbox.hpp"

//WindowingTools
/*! \brief Returns the value of a windowing function given by fctIdx at CurTimeIdx with given endTimeIdx (i.e. endTimeIdx=nTimeIter, if one starts  windowing at time t =0.) */
su2double WindowingTools::GetWndWeight(int fctIdx, unsigned long CurTimeIdx, unsigned long endTimeIdx){
    switch (fctIdx){
      case 1: return HannWindow(CurTimeIdx, endTimeIdx);
      case 2: return HannSquaredWindow(CurTimeIdx, endTimeIdx);
      case 3: return BumpWindow(CurTimeIdx, endTimeIdx);
      default:return 1.0;
    }
}
su2double WindowingTools::HannWindow(unsigned long i, unsigned long endTimeIdx){
  if(i==0) return 0; //catch div by zero error
  su2double currTime = static_cast<su2double>(i);
  su2double endTime = static_cast<su2double>(endTimeIdx);
  su2double tau = currTime/endTime;
  return 1.0-cos(2*PI_NUMBER*tau);
}

su2double WindowingTools::HannSquaredWindow(unsigned long i, unsigned long endTimeIdx){
  if(i==0) return 0; //catch div by zero error
  su2double currTime = static_cast<su2double>(i);
  su2double endTime = static_cast<su2double>(endTimeIdx);
  su2double tau = currTime/endTime;
  return 2.0/3.0*(1-cos(2*PI_NUMBER*tau))*(1-cos(2*PI_NUMBER*tau));
}

su2double WindowingTools::BumpWindow(unsigned long i, unsigned long endTimeIdx){
  if(i==0) return 0;
  if(i==endTimeIdx) return 0;
  su2double currTime = static_cast<su2double>(i);
  su2double endTime = static_cast<su2double>(endTimeIdx);
  su2double tau = currTime/endTime;
  return 1.0/0.00702986*(exp(-1/(tau-tau*tau)));
}


//WindowedAverage
WindowedAverage::WindowedAverage(){
  this->Reset();
}

su2double WindowedAverage::Update(su2double valIn){
      su2double scaling = 1. / (su2double)(count + 1);
      val = valIn * scaling + val * (1. - scaling);
      count++;
      return val;
}

su2double WindowedAverage::Get(){
  return val;
}

unsigned long WindowedAverage::Count(){
  return count;
}

void WindowedAverage::Reset(){

  val = 0.;
  count = 0;
}


/*
void WindowedAverage::addValue(su2double valIn, unsigned long currentIter,unsigned long startIter){
    if(currentIter >= startIter)values.push_back(valIn);
}

su2double WindowedAverage::WindowedUpdate(int fctIdx){ //Computes a windowed time average (integral)
  if(values.size()>0){
      switch (fctIdx){
        case 1: val    = HannWindowing();         return val;
        case 2: val    = HannSquaredWindowing();  return val;
        case 3: val    = BumpWindowing();         return val;
        default: val   = NoWindowing();           return val;
      }
  }
  return 0.0;
}

su2double WindowedAverage::NoWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned long i=0; i<values.size(); i++){
      wnd_timeAvg+=values[i];
    }
  return wnd_timeAvg/static_cast<su2double>(values.size());
}

su2double WindowedAverage::HannWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned long i=0; i<values.size(); i++){
      wnd_timeAvg+=values[i]*HannWindow(i,values.size()-1);
  }
  return wnd_timeAvg/static_cast<su2double>(values.size());
}

su2double WindowedAverage::HannSquaredWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned long i=0; i<values.size(); i++){
      wnd_timeAvg+=values[i]*HannSquaredWindow(i,values.size()-1);
    }
  return wnd_timeAvg/static_cast<su2double>(values.size());
}

su2double WindowedAverage::BumpWindowing(){
  su2double wnd_timeAvg = 0.0;
  for(unsigned long i=0; i<values.size(); i++){
      wnd_timeAvg+=values[i]*BumpWindow(i,values.size()-1);
    }
  return wnd_timeAvg/static_cast<su2double>(values.size());
}
*/

