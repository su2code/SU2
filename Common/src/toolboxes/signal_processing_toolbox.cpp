#include "../../include/toolboxes/signal_processing_toolbox.hpp"


su2double Signal_Processing::Average(const std::vector<su2double> &data){
  su2double avg = 0.0;
  for (const su2double& val : data){
    avg += val;
  }
  return avg/data.size();
}
