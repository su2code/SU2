/*!
 * \file CConjugateHeatInterface.cpp
 * \brief Declaration and inlines of the class to transfer temperature and heatflux
 *        density for conjugate heat interfaces between structure and fluid zones.
 * \author O. Burghardt
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../../include/interfaces/cht/CConjugateHeatInterface.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../include/solvers/CSolver.hpp"

CConjugateHeatInterface::CConjugateHeatInterface(unsigned short val_nVar, unsigned short val_nConst) :
  CInterface(val_nVar, val_nConst) {
}

void CConjugateHeatInterface::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry,
                                                const CConfig *donor_config, unsigned long Marker_Donor,
                                                unsigned long Vertex_Donor, unsigned long Point_Donor) {

  /*--- Compute distance of donor point to PointNormal for T-gradient/heatflux computation ---*/
  const auto nDim = donor_geometry->GetnDim();
  const auto Coord = donor_geometry->nodes->GetCoord(Point_Donor);
  const auto PointNormal = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNormal_Neighbor();
  const auto Coord_Normal = donor_geometry->nodes->GetCoord(PointNormal);

  su2double Edge_Vector[MAXNDIM] = {0.0};
  GeometryToolbox::Distance(nDim, Coord_Normal, Coord, Edge_Vector);
  const su2double dist = GeometryToolbox::Norm(nDim, Edge_Vector);

  /*--- Retrieve temperature solution and its gradient ---*/

  const su2double Twall   = donor_solution->GetNodes()->GetTemperature(Point_Donor);
  const su2double Tnormal = donor_solution->GetNodes()->GetTemperature(PointNormal);

  const su2double dTdn = (Twall - Tnormal)/dist;
  // TODO: Check if these improve accuracy, if needed at all
  //    const auto Normal = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNormal();
  //    su2double Area = GeometryToolbox::Norm(nDim, Normal);
  //    for (iDim = 0; iDim < nDim; iDim++) {
  //      dTdn += (Twall - Tnormal)/dist * (Edge_Vector[iDim]/dist) * (Normal[iDim]/Area);
  //    }

  /*--- Calculate the heat flux density (temperature gradient times thermal conductivity) and
        thermal conductivity divided by distance. ---*/
  su2double thermal_conductivity = 0.0;
  su2double heat_flux_density = 0.0;
  su2double conductivity_over_dist = 0.0;

  const bool compressible_flow = (donor_config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);
  const bool incompressible_flow = (donor_config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE) &&
                                   (donor_config->GetEnergy_Equation() || (donor_config->GetKind_FluidModel() == ENUM_FLUIDMODEL::FLUID_FLAMELET));

  if (compressible_flow) {

    const su2double Gamma = donor_config->GetGamma();
    const su2double Gas_Constant = donor_config->GetGas_ConstantND();
    const su2double Cp = (Gamma / (Gamma - 1.0)) * Gas_Constant;

    const su2double Prandtl_Lam = donor_config->GetPrandtl_Lam();
    const su2double laminar_viscosity = donor_solution->GetNodes()->GetLaminarViscosity(Point_Donor); // TDE check for consistency

    const su2double thermal_conductivityND = Cp*(laminar_viscosity/Prandtl_Lam);
    heat_flux_density = thermal_conductivityND*dTdn;

    if ((donor_config->GetKind_CHT_Coupling() == CHT_COUPLING::DIRECT_TEMPERATURE_ROBIN_HEATFLUX) ||
        (donor_config->GetKind_CHT_Coupling() == CHT_COUPLING::AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)) {

      thermal_conductivity   = thermal_conductivityND*donor_config->GetViscosity_Ref();
      conductivity_over_dist = thermal_conductivity/dist;
    }
  }
  else if (incompressible_flow) {

    const auto iPoint = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNode();

    const su2double thermal_conductivityND  = donor_solution->GetNodes()->GetThermalConductivity(iPoint);
    heat_flux_density       = thermal_conductivityND*dTdn;

    if ((donor_config->GetKind_CHT_Coupling() == CHT_COUPLING::DIRECT_TEMPERATURE_ROBIN_HEATFLUX) ||
        (donor_config->GetKind_CHT_Coupling() == CHT_COUPLING::AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)) {

      switch (donor_config->GetKind_ConductivityModel()) {

        case CONDUCTIVITYMODEL::CONSTANT:
        case CONDUCTIVITYMODEL::FLAMELET:
        case CONDUCTIVITYMODEL::COOLPROP:
          thermal_conductivity = thermal_conductivityND*donor_config->GetThermal_Conductivity_Ref();
          break;
        case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
          thermal_conductivity = thermal_conductivityND*donor_config->GetGas_Constant_Ref()
                                 *donor_config->GetViscosity_Ref();
          break;
        case CONDUCTIVITYMODEL::POLYNOMIAL:
          SU2_MPI::Error("Polynomial Conductivity model not implemented for CHT interface.", CURRENT_FUNCTION);
          break;
      }

      conductivity_over_dist  = thermal_conductivity/dist;
    }
  }
  else if (donor_config->GetHeatProblem()) {

    /*--- Heat solver stand-alone case ---*/

    const su2double thermal_diffusivity = donor_config->GetThermalDiffusivity();
    heat_flux_density = thermal_diffusivity*dTdn;


    if ((donor_config->GetKind_CHT_Coupling() == CHT_COUPLING::DIRECT_TEMPERATURE_ROBIN_HEATFLUX) ||
        (donor_config->GetKind_CHT_Coupling() == CHT_COUPLING::AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)) {

      const su2double rho_cp_solid = donor_config->GetSpecific_Heat_Cp()*donor_config->GetMaterialDensity(0);
      conductivity_over_dist  = thermal_diffusivity*rho_cp_solid/dist;
    }
  }

  /*--- Set the conjugate heat variables that are transferred by default ---*/

  Donor_Variable[0] = Twall*donor_config->GetTemperature_Ref();
  Donor_Variable[1] = heat_flux_density*donor_config->GetHeat_Flux_Ref();

  /*--- We only need these for the Robin BC option ---*/

  if ((donor_config->GetKind_CHT_Coupling() == CHT_COUPLING::DIRECT_TEMPERATURE_ROBIN_HEATFLUX) ||
      (donor_config->GetKind_CHT_Coupling() == CHT_COUPLING::AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)) {

    Donor_Variable[2] = conductivity_over_dist;
    Donor_Variable[3] = Tnormal*donor_config->GetTemperature_Ref();
  }

}

void CConjugateHeatInterface::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                                                 const CConfig *target_config, unsigned long Marker_Target,
                                                 unsigned long Vertex_Target, unsigned long Point_Target) {

  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 0,
                                            target_config->GetRelaxation_Factor_CHT(), Target_Variable[0]);
  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 1,
                                            target_config->GetRelaxation_Factor_CHT(), Target_Variable[1]);

  if ((target_config->GetKind_CHT_Coupling() == CHT_COUPLING::DIRECT_TEMPERATURE_ROBIN_HEATFLUX) ||
      (target_config->GetKind_CHT_Coupling() == CHT_COUPLING::AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)) {

    target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 2,
                                              target_config->GetRelaxation_Factor_CHT(), Target_Variable[2]);
    target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 3,
                                              target_config->GetRelaxation_Factor_CHT(), Target_Variable[3]);
  }
}
