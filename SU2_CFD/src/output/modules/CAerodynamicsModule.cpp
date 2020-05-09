#include "../../../include/output/modules/CAerodynamicsModule.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../../Common/include/fem_geometry_structure.hpp"

#include "../../../include/solvers/CSolver.hpp"

CAerodynamicsModule::CAerodynamicsModule(CConfig *config):
  Alpha(config->GetAoA()*PI_NUMBER/180.0),
  Beta(config->GetAoS()*PI_NUMBER/180.0),
  RefArea(config->GetRefArea()),
  RefLength(config->GetRefLength()),
  Gas_Constant(config->GetGas_Constant()),
  Gamma(config->GetGamma()),
  Prandtl_Lam(config->GetPrandtl_Lam()),
  axisymmetric(config->GetAxisymmetric()){

}

void CAerodynamicsModule::LoadHistoryDataPerSurface(COutFieldCollection &fieldCollection){

  const int nDim = solverData.geometry->GetnDim();
  const su2double Force[3] = {fieldCollection.GetValueByKey("FORCE_X"),
                              fieldCollection.GetValueByKey("FORCE_Y"),
                              fieldCollection.GetValueByKey("FORCE_Z")};

  const su2double alpha = solverData.config->GetAoA()*PI_NUMBER/180.0;
  if (nDim == 2){
    const su2double Drag = Force[0]*cos(alpha) + Force[1]*sin(alpha);
    const su2double Lift = -Force[0]*sin(alpha) + Force[1]*cos(alpha);
    fieldCollection.SetValueByKey("DRAG", Drag);
    fieldCollection.SetValueByKey("LIFT", Lift);
  } else {
    const su2double beta = solverData.config->GetAoS()*PI_NUMBER/180.0;
    const su2double Drag = Force[0]*cos(alpha)*cos(beta)  + Force[1]*sin(beta) + Force[2]*sin(alpha)*cos(beta);
    const su2double Lift = -Force[0]*sin(alpha)*cos(beta) + Force[2]*cos(alpha);
    const su2double Sideforce = Force[0]*sin(beta)*cos(alpha) + Force[1]*cos(beta) - Force[2]*sin(beta)*sin(alpha);
    fieldCollection.SetValueByKey("DRAG", Drag);
    fieldCollection.SetValueByKey("LIFT", Lift);
    fieldCollection.SetValueByKey("SIDEFORCE", Sideforce);
  }
}

void CAerodynamicsModule::DefineHistoryFields(COutFieldCollection &fieldCollection){

  const std::string aeroGroupName = "AERO_COEFF";

  fieldCollection.AddItem("LIFT",      COutputField("CL", ScreenOutputFormat::FIXED, aeroGroupName,
                                                     FieldType::COEFFICIENT, "Lift Coefficient"));
  fieldCollection.AddItem("DRAG",      COutputField("CL", ScreenOutputFormat::FIXED, aeroGroupName,
                                                     FieldType::COEFFICIENT, "Drag Coefficient"));
  fieldCollection.AddItem("SIDEFORCE", COutputField("CL", ScreenOutputFormat::FIXED, aeroGroupName,
                                                     FieldType::COEFFICIENT, "Sideforce Coefficient"));

}


void CAerodynamicsModule::DefineVolumeFields(COutFieldCollection &fieldCollection){

  std::string aeroGroupName = "AERO_COEFF";

  fieldCollection.AddItem("MOMENT_X", COutputField("CMx", -1, aeroGroupName, "Moment around the x-axis", FieldType::SURFACE_INTEGRATE));
  fieldCollection.AddItem("MOMENT_Y", COutputField("CMy", -1, aeroGroupName, "Moment around the y-axis", FieldType::SURFACE_INTEGRATE));
  fieldCollection.AddItem("MOMENT_Z", COutputField("CMz", -1, aeroGroupName, "Moment around the z-axis", FieldType::SURFACE_INTEGRATE));

  fieldCollection.AddItem("FORCE_X", COutputField("Fx", -1, aeroGroupName, "Force in x-direction", FieldType::SURFACE_INTEGRATE));
  fieldCollection.AddItem("FORCE_Y", COutputField("Fy", -1, aeroGroupName, "Force in y-direction", FieldType::SURFACE_INTEGRATE));
  fieldCollection.AddItem("FORCE_Z", COutputField("Fz", -1, aeroGroupName, "Force in z-direction", FieldType::SURFACE_INTEGRATE));

  fieldCollection.AddItem("SKIN_FRICTION_X", COutputField("Skin_Friction_Coefficient_x", -1, aeroGroupName, "Skin friction coefficient x-component", FieldType::SURFACE_INTEGRATE));
  fieldCollection.AddItem("SKIN_FRICTION_Y", COutputField("Skin_Friction_Coefficient_y", -1, aeroGroupName, "Skin friction coefficient y-component", FieldType::SURFACE_INTEGRATE));
  fieldCollection.AddItem("SKIN_FRICTION_Z", COutputField("Skin_Friction_Coefficient_z", -1, aeroGroupName, "Skin friction coefficient z-component", FieldType::SURFACE_INTEGRATE));

  fieldCollection.AddItem("HEATFLUX", COutputField("Heatflux", -1, aeroGroupName, "Heatflux", FieldType::SURFACE_INTEGRATE));
  fieldCollection.AddItem("Y_PLUS", COutputField("Y_Plus", -1, aeroGroupName, "Non-dim. wall distance (Y-Plus)", FieldType::DEFAULT));
}

void CAerodynamicsModule::LoadSurfaceData(COutFieldCollection& fieldCollection){

  std::fill(std::begin(Moment), std::end(Moment), 0.0);
  std::fill(std::begin(Force), std::end(Force), 0.0);
  std::fill(std::begin(SkinFriction), std::end(SkinFriction), 0.0);
  Heatflux = 0.0;

  const su2double RefDensity     = solverData.config->GetDensity_FreeStreamND();
  const su2double RefTemperature = solverData.config->GetTemperature_FreeStreamND();
  const su2double* RefVelocity   = solverData.config->GetVelocity_FreeStreamND();
  const su2double RefHeatFlux    = solverData.config->GetHeat_Flux_Ref();
  const int nDim = solverData.geometry->GetnDim();
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
  const su2double* Coord = solverData.geometry->node[iPoint]->GetCoord();
  const su2double* Coord_Normal = solverData.geometry->node[iPointNormal]->GetCoord();
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

  fieldCollection.SetValueByKey("MOMENT_X", Moment[0]);
  fieldCollection.SetValueByKey("MOMENT_Y", Moment[1]);
  fieldCollection.SetValueByKey("MOMENT_Z", Moment[2]);
  fieldCollection.SetValueByKey("FORCE_X", Force[0]);
  fieldCollection.SetValueByKey("FORCE_Y", Force[1]);
  fieldCollection.SetValueByKey("FORCE_Z", Force[2]);
  fieldCollection.SetValueByKey("HEATFLUX", Heatflux);
  fieldCollection.SetValueByKey("SKIN_FRICTION_X", SkinFriction[0]);
  fieldCollection.SetValueByKey("SKIN_FRICTION_Y", SkinFriction[1]);
  fieldCollection.SetValueByKey("SKIN_FRICTION_Z", SkinFriction[2]);
  fieldCollection.SetValueByKey("Y_PLUS", Y_Plus);
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