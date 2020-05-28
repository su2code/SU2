#include "../../../include/output/modules/CVortexIdentificationModule.hpp"

#include "../../../../Common/include/CConfig.hpp"
#include "../../../include/solvers/CSolver.hpp"

void CVortexIdentificationModule::DefineVolumeFields(CVolumeOutFieldManager& volumeFields) {
  if (nDim == 3) {
    volumeFields.AddField("VORTICITY_X", "Vorticity_x", "VORTEX_IDENTIFICATION", "x-component of the vorticity vector",
                          FieldType::DEFAULT);
    volumeFields.AddField("VORTICITY_Y", "Vorticity_y", "VORTEX_IDENTIFICATION", "y-component of the vorticity vector",
                          FieldType::DEFAULT);
    volumeFields.AddField("VORTICITY_Z", "Vorticity_z", "VORTEX_IDENTIFICATION", "z-component of the vorticity vector",
                          FieldType::DEFAULT);
  } else {
    volumeFields.AddField("VORTICITY", "Vorticity", "VORTEX_IDENTIFICATION", "Value of the vorticity",
                          FieldType::DEFAULT);
  }
  volumeFields.AddField("Q_CRITERION", "Q_Criterion", "VORTEX_IDENTIFICATION", "Value of the Q-Criterion",
                        FieldType::DEFAULT);
}

void CVortexIdentificationModule::LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                                 const IterationInfo&, const PointInfo& pointInfo) {

  auto* flow_solver = solverData.solver[FLOW_SOL];
  auto* Node_Flow = flow_solver->GetNodes();
  const auto iPoint = pointInfo.iPoint;

  if (nDim == 3) {
    volumeFields.SetFieldValue("VORTICITY_X", Node_Flow->GetVorticity(iPoint)[0]);
    volumeFields.SetFieldValue("VORTICITY_Y", Node_Flow->GetVorticity(iPoint)[1]);
    volumeFields.SetFieldValue("VORTICITY_Z", Node_Flow->GetVorticity(iPoint)[2]);
  } else {
    volumeFields.SetFieldValue("VORTICITY", Node_Flow->GetVorticity(iPoint)[2]);
  }
  volumeFields.SetFieldValue("Q_CRITERION", GetQ_Criterion(&(Node_Flow->GetGradient_Primitive(iPoint)[1])));
}

su2double CVortexIdentificationModule::GetQ_Criterion(su2double** VelocityGradient) const {
  /*--- Make a 3D copy of the gradient so we do not have worry about nDim ---*/

  su2double Grad_Vel[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    for (unsigned short jDim = 0; jDim < nDim; jDim++) Grad_Vel[iDim][jDim] = VelocityGradient[iDim][jDim];

  /*--- Q Criterion Eq 1.2 of HALLER, G. (2005). An objective definition of a vortex.
   Journal of Fluid Mechanics, 525, 1-26. doi:10.1017/S0022112004002526 ---*/

  /*--- Components of the strain rate tensor (symmetric) ---*/
  const su2double s11 = Grad_Vel[0][0];
  const su2double s12 = 0.5 * (Grad_Vel[0][1] + Grad_Vel[1][0]);
  const su2double s13 = 0.5 * (Grad_Vel[0][2] + Grad_Vel[2][0]);
  const su2double s22 = Grad_Vel[1][1];
  const su2double s23 = 0.5 * (Grad_Vel[1][2] + Grad_Vel[2][1]);
  const su2double s33 = Grad_Vel[2][2];

  /*--- Components of the spin tensor (skew-symmetric) ---*/
  const su2double omega12 = 0.5 * (Grad_Vel[0][1] - Grad_Vel[1][0]);
  const su2double omega13 = 0.5 * (Grad_Vel[0][2] - Grad_Vel[2][0]);
  const su2double omega23 = 0.5 * (Grad_Vel[1][2] - Grad_Vel[2][1]);

  /*--- Q = ||Omega|| - ||Strain|| ---*/
  const su2double Q = 2 * (pow(omega12, 2) + pow(omega13, 2) + pow(omega23, 2)) -
                      (pow(s11, 2) + pow(s22, 2) + pow(s33, 2) + 2 * (pow(s12, 2) + pow(s13, 2) + pow(s23, 2)));

  return Q;
}