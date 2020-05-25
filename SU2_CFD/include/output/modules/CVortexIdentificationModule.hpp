#pragma once

#include "COutputModule.hpp"

class CVortexIdentificationModule : public CSolverOutputModule{

public:
  CVortexIdentificationModule(CConfig* config, int nDim) : CSolverOutputModule(nDim, config->GetViscous()){}

  void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) override;

  void LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                      const IterationInfo& iterationInfo, const PointInfo& pointInfo) override;

  /*!
   * \brief Compute value of the Q criteration for vortex idenfitication
   * \param[in] VelocityGradient - Velocity gradients
   * \return Value of the Q criteration at the node
   */
  su2double GetQ_Criterion(su2double** VelocityGradient) const;

};
