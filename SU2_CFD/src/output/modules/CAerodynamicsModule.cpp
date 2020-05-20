#include "../../../include/output/modules/CAerodynamicsModule.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../../Common/include/fem_geometry_structure.hpp"

#include "../../../include/solvers/CSolver.hpp"

CAerodynamicsModule::CAerodynamicsModule(CConfig *config, int nDim): CSolverOutputModule(nDim),
  Alpha(config->GetAoA()*PI_NUMBER/180.0),
  Beta(config->GetAoS()*PI_NUMBER/180.0),
  RefLength(config->GetRefLength()),
  Gas_Constant(config->GetGas_Constant()),
  Gamma(config->GetGamma()),
  Prandtl_Lam(config->GetPrandtl_Lam()),
  axisymmetric(config->GetAxisymmetric()){

}

void CAerodynamicsModule::LoadHistoryDataPerSurface(CHistoryOutFieldManager& historyFields){

  const int nDim = solverData.geometry->GetnDim();
  const su2double Force[3] = {historyFields.GetFieldValue("FORCE_X"),
                              historyFields.GetFieldValue("FORCE_Y"),
                              historyFields.GetFieldValue("FORCE_Z")};

  const su2double alpha = solverData.config->GetAoA()*PI_NUMBER/180.0;
  if (nDim == 2){
    const su2double Drag = Force[0]*cos(alpha) + Force[1]*sin(alpha);
    const su2double Lift = -Force[0]*sin(alpha) + Force[1]*cos(alpha);
    historyFields.SetFieldValue("DRAG", Drag);
    historyFields.SetFieldValue("LIFT", Lift);
  } else {
    const su2double beta = solverData.config->GetAoS()*PI_NUMBER/180.0;
    const su2double Drag = Force[0]*cos(alpha)*cos(beta)  + Force[1]*sin(beta) + Force[2]*sin(alpha)*cos(beta);
    const su2double Lift = -Force[0]*sin(alpha)*cos(beta) + Force[2]*cos(alpha);
    const su2double Sideforce = Force[0]*sin(beta)*cos(alpha) + Force[1]*cos(beta) - Force[2]*sin(beta)*sin(alpha);
    historyFields.SetFieldValue("DRAG", Drag);
    historyFields.SetFieldValue("LIFT", Lift);
    historyFields.SetFieldValue("SIDEFORCE", Sideforce);
  }
}

void CAerodynamicsModule::DefineHistoryFields(CHistoryOutFieldManager& historyFields){

  const std::string aeroGroupName = "AERO_COEFF";

  historyFields.AddField("LIFT",  "CL", ScreenOutputFormat::FIXED, aeroGroupName,
                         "Lift Coefficient", FieldType::COEFFICIENT);
  historyFields.AddField("DRAG", "CD", ScreenOutputFormat::FIXED, aeroGroupName,
                         "Drag Coefficient", FieldType::COEFFICIENT);
  historyFields.AddField("SIDEFORCE", "CSF", ScreenOutputFormat::FIXED, aeroGroupName,
                         "Sideforce Coefficient", FieldType::COEFFICIENT);

}


void CAerodynamicsModule::DefineVolumeFields(CVolumeOutFieldManager &volumeFields){

  std::string aeroGroupName = "AERO_COEFF";

  volumeFields.AddField("MOMENT_X", "CMx", aeroGroupName, "Moment around the x-axis", FieldType::SURFACE_INTEGRATE);
  volumeFields.AddField("MOMENT_Y", "CMy", aeroGroupName, "Moment around the y-axis", FieldType::SURFACE_INTEGRATE);
  volumeFields.AddField("MOMENT_Z", "CMz", aeroGroupName, "Moment around the z-axis", FieldType::SURFACE_INTEGRATE);

  volumeFields.AddField("FORCE_X", "Fx", aeroGroupName, "Force in x-direction", FieldType::SURFACE_INTEGRATE);
  volumeFields.AddField("FORCE_Y", "Fy", aeroGroupName, "Force in y-direction", FieldType::SURFACE_INTEGRATE);
  volumeFields.AddField("FORCE_Z", "Fz", aeroGroupName, "Force in z-direction", FieldType::SURFACE_INTEGRATE);

  volumeFields.AddField("SKIN_FRICTION_X", "Skin_Friction_Coefficient_x",  aeroGroupName, "Skin friction coefficient x-component", FieldType::SURFACE_INTEGRATE);
  volumeFields.AddField("SKIN_FRICTION_Y", "Skin_Friction_Coefficient_y",  aeroGroupName, "Skin friction coefficient y-component", FieldType::SURFACE_INTEGRATE);
  volumeFields.AddField("SKIN_FRICTION_Z", "Skin_Friction_Coefficient_z",  aeroGroupName, "Skin friction coefficient z-component", FieldType::SURFACE_INTEGRATE);

  volumeFields.AddField("HEATFLUX", "Heatflux", aeroGroupName, "Heatflux", FieldType::SURFACE_INTEGRATE);
  volumeFields.AddField("Y_PLUS", "Y_Plus", aeroGroupName, "Non-dim. wall distance (Y-Plus)", FieldType::DEFAULT);
}

void CAerodynamicsModule::LoadSurfaceData(CVolumeOutFieldManager& volumeFields){

  std::fill(std::begin(Moment), std::end(Moment), 0.0);
  std::fill(std::begin(Force), std::end(Force), 0.0);
  std::fill(std::begin(SkinFriction), std::end(SkinFriction), 0.0);
  Heatflux = 0.0;

  const su2double RefDensity     = solverData.config->GetDensity_FreeStreamND();
  const su2double RefTemperature = solverData.config->GetTemperature_FreeStreamND();
  const su2double* RefVelocity   = solverData.config->GetVelocity_FreeStreamND();
  const su2double RefHeatFlux    = solverData.config->GetHeat_Flux_Ref();
  const su2double RefArea        = solverData.config->GetRefArea();

  su2double RefVel2;
  if (solverData.config->GetDynamic_Grid()) {
    const su2double Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemperature);
    const su2double Mach_Motion = solverData.config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (int iDim = 0; iDim < nDim; iDim++)
      RefVel2  += RefVelocity[iDim]*RefVelocity[iDim];
  }
  const su2double factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);
  const unsigned long iPoint = solverData.iPoint;
  const unsigned long iPointNormal = solverData.vertex->GetNormal_Neighbor();

  const su2double* Normal = solverData.vertex->GetNormal();
  const su2double* Coord = solverData.geometry->nodes->GetCoord(iPoint);
  const su2double* Coord_Normal = solverData.geometry->nodes->GetCoord(iPointNormal);
  const su2double Pressure_Inf = solverData.solver[FLOW_SOL]->GetPressure_Inf();

  // TODO: Get ref origin from config
  su2double MomentDist[3] = {0.0}, Origin[3] = {0.0};

//  Origin[0] = data.config->GetRefOriginMoment_X(0);
//  Origin[1] = data.config->GetRefOriginMoment_Y(0);
//  Origin[2] = data.config->GetRefOriginMoment_Z(0);

  for (int iDim = 0; iDim < nDim; iDim++) {
    MomentDist[iDim] = Coord[iDim] - Origin[iDim];
  }

  /*--- Axisymmetric simulations ---*/
  su2double AxiFactor = 0.0;
  if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*Coord[1];
  else AxiFactor = 1.0;

  const su2double Pressure = solverData.solver[FLOW_SOL]->GetNodes()->GetPressure(iPoint);

  for (int iDim = 0; iDim < nDim; iDim++) {
    Force[iDim] = -(Pressure - Pressure_Inf) * Normal[iDim] * factor * AxiFactor;
  }

  if (solverData.config->GetViscous()){

    const auto& Grad_Primitive = solverData.solver[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);
    su2double div_vel = 0.0; for (int iDim = 0; iDim < nDim; iDim++) div_vel += Grad_Primitive[iDim+1][iDim];
    const su2double Viscosity = solverData.solver[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    const su2double Density = solverData.solver[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    constexpr passivedouble delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
    su2double Area = 0.0; for (int iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    su2double UnitNormal[3];
    for (int iDim = 0; iDim < nDim; iDim++) {
      UnitNormal[iDim] = Normal[iDim]/Area;
    }
    su2double Tau[3][3] = {{0.0}};

    for (int iDim = 0; iDim < nDim; iDim++) {
      for (int jDim = 0; jDim < nDim; jDim++) {
        Tau[iDim][jDim] = Viscosity*(Grad_Primitive[jDim+1][iDim] + Grad_Primitive[iDim+1][jDim]) - TWO3*Viscosity*div_vel*delta[iDim][jDim];
      }
    }

    if (solverData.config->GetQCR()) {
      su2double den_aux, c_cr1=0.3, O_ik, O_jk;
      unsigned short kDim;

      /*--- Denominator Antisymmetric normalized rotation tensor ---*/

      den_aux = 0.0;
      for (int iDim = 0 ; iDim < nDim; iDim++)
        for (int jDim = 0 ; jDim < nDim; jDim++)
          den_aux += Grad_Primitive[iDim+1][jDim] * Grad_Primitive[iDim+1][jDim];
      den_aux = sqrt(max(den_aux,1E-10));

      /*--- Adding the QCR contribution ---*/

      for (int iDim = 0 ; iDim < nDim; iDim++){
        for (int jDim = 0 ; jDim < nDim; jDim++){
          for (kDim = 0 ; kDim < nDim; kDim++){
            O_ik = (Grad_Primitive[iDim+1][kDim] - Grad_Primitive[kDim+1][iDim])/ den_aux;
            O_jk = (Grad_Primitive[jDim+1][kDim] - Grad_Primitive[kDim+1][jDim])/ den_aux;
            Tau[iDim][jDim] -= c_cr1 * (O_ik * Tau[jDim][kDim] + O_jk * Tau[iDim][kDim]);
          }
        }
      }
    }

    su2double TauElem[3] = {0.0};
    for (int iDim = 0; iDim < nDim; iDim++) {
      for (int jDim = 0; jDim < nDim; jDim++) {
        TauElem[iDim] += Tau[iDim][jDim]*UnitNormal[jDim];
      }
    }
    su2double TauNormal = 0.0; for (int iDim = 0; iDim < nDim; iDim++) TauNormal += TauElem[iDim] * UnitNormal[iDim];
    su2double WallShearStress = 0.0;
    su2double TauTangent[3] = {0.0};
    for (int iDim = 0; iDim < nDim; iDim++) {
      TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
      SkinFriction[iDim] = TauTangent[iDim] / (0.5*RefDensity*RefVel2);
      WallShearStress += TauTangent[iDim] * TauTangent[iDim];
    }
    WallShearStress = sqrt(WallShearStress);

    su2double WallDist[3];
    for (int iDim = 0; iDim < nDim; iDim++) WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
    su2double WallDistMod = 0.0; for (int iDim = 0; iDim < nDim; iDim++) WallDistMod += WallDist[iDim]*WallDist[iDim];
    WallDistMod = sqrt(WallDistMod);

    /*--- Compute y+ and non-dimensional velocity ---*/

    const su2double FrictionVel = sqrt(fabs(WallShearStress)/Density);
    Y_Plus = WallDistMod*FrictionVel/(Viscosity/Density);

    su2double GradTemperature = 0.0;
    for (int iDim = 0; iDim < nDim; iDim++)
      GradTemperature -= Grad_Primitive[0][iDim]*UnitNormal[iDim];

    const su2double Cp = Gamma / (Gamma - 1) * Gas_Constant;
    const su2double thermal_conductivity = Cp * Viscosity/Prandtl_Lam;
    Heatflux = -thermal_conductivity*GradTemperature*RefHeatFlux;

    for(int iDim = 0; iDim < nDim; iDim++){
      Force[iDim] += TauElem[iDim]*Area*factor*AxiFactor;
    }
  }


  if (nDim == 3){
    Moment[0] += (Force[2]*MomentDist[1] - Force[1]*MomentDist[2])/RefLength;
    Moment[1] += (Force[0]*MomentDist[2] - Force[2]*MomentDist[0])/RefLength;
  }
  Moment[2] += (Force[1]*MomentDist[0] - Force[0]*MomentDist[1])/RefLength;

  volumeFields.SetFieldValue("MOMENT_X", Moment[0]);
  volumeFields.SetFieldValue("MOMENT_Y", Moment[1]);
  volumeFields.SetFieldValue("MOMENT_Z", Moment[2]);
  volumeFields.SetFieldValue("FORCE_X", Force[0]);
  volumeFields.SetFieldValue("FORCE_Y", Force[1]);
  volumeFields.SetFieldValue("FORCE_Z", Force[2]);
  volumeFields.SetFieldValue("HEATFLUX", Heatflux);
  volumeFields.SetFieldValue("SKIN_FRICTION_X", SkinFriction[0]);
  volumeFields.SetFieldValue("SKIN_FRICTION_Y", SkinFriction[1]);
  volumeFields.SetFieldValue("SKIN_FRICTION_Z", SkinFriction[2]);
  volumeFields.SetFieldValue("Y_PLUS", Y_Plus);
}


//void CAerodynamicsModuleFEM::LoadSurfaceData(COutFieldCollection &fieldCollection){

//  const int nDim = solverData.geometry->GetnDim();

//  const su2double RefDensity     = solverData.config->GetDensity_FreeStreamND();
//  const su2double RefTemperature = solverData.config->GetTemperature_FreeStreamND();
//  const su2double* RefVelocity   = solverData.config->GetVelocity_FreeStreamND();
//  const su2double RefHeatFlux    = solverData.config->GetHeat_Flux_Ref();
//  su2double RefVel2;
//  if (solverData.config->GetDynamic_Grid()) {
//    const su2double Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemperature);
//    const su2double Mach_Motion = solverData.config->GetMach_Motion();
//    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
//  }
//  else {
//    RefVel2 = 0.0;
//    for (int iDim = 0; iDim < nDim; iDim++)
//      RefVel2  += RefVelocity[iDim]*RefVelocity[iDim];
//  }


//  std::fill(std::begin(Moment), std::end(Moment), 0.0);
//  std::fill(std::begin(Force), std::end(Force), 0.0);
//  std::fill(std::begin(SkinFriction), std::end(SkinFriction), 0.0);
//  Heatflux = 0.0;

//  /* Easier storage of the solution, the normals and the coordinates
//     for this integration point. */
//  const su2double *sol     = solverData.sol;
//  const su2double *normals = solverData.surfElem->metricNormalsFace.data()
//                           + solverData.i*(nDim+1);
//  const su2double *Coord   = solverData.surfElem->coorIntegrationPoints.data()
//                           + solverData.i*nDim;

//  const auto weight = solverData.weight;
//  const su2double factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

//  /*--- Compute the pressure in this integration point. ---*/
//  const su2double DensityInv   = 1.0/sol[0];
//  const su2double u            = sol[1]*DensityInv;
//  const su2double v            = sol[2]*DensityInv;
//  const su2double StaticEnergy = sol[3]*DensityInv - 0.5*(u*u + v*v);

//  solverData.solver[FLOW_SOL]->GetFluidModel()->SetTDState_rhoe(sol[0], StaticEnergy);
//  const su2double Pressure =  solverData.solver[FLOW_SOL]->GetFluidModel()->GetPressure();

//  su2double Origin[3] = {0.0};
//  const su2double Pressure_Inf = solverData.solver[FLOW_SOL]->GetPressure_Inf();

//  /*-- Compute the vector from the reference point to the integration
//       point and update the inviscid force. Note that the normal points
//       into the geometry, hence no minus sign. ---*/
//  const su2double DistX = Coord[0] - Origin[0];
//  const su2double DistY = Coord[1] - Origin[1];

//  const su2double ForceMag = (Pressure - Pressure_Inf)*weight
//                           * normals[nDim]*factor;
//  Force[0]   = ForceMag*normals[0];
//  Force[1]   = ForceMag*normals[1];

//  /*--- Update the inviscid moment. ---*/
//  Moment[2] += (Force[1]*DistX - Force[0]*DistY)/RefLength;

//}