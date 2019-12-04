/*!
 * \file CConjugateHeatInterface.cpp
 * \brief Declaration and inlines of the class to transfer temperature and heatflux
 *        density for conjugate heat interfaces between structure and fluid zones.
 * \author O. Burghardt
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

CConjugateHeatInterface::CConjugateHeatInterface(void) : CInterface() { }

CConjugateHeatInterface::CConjugateHeatInterface(unsigned short val_nVar, unsigned short val_nConst,
                                                 CConfig *config) : CInterface(val_nVar, val_nConst, config) { }

CConjugateHeatInterface::~CConjugateHeatInterface(void) { }

void CConjugateHeatInterface::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry,
                                                CConfig *donor_config, unsigned long Marker_Donor,
                                                unsigned long Vertex_Donor, unsigned long Point_Donor) {

  unsigned short nDim, iDim;
  unsigned long iPoint, PointNormal;

  su2double *Coord, *Coord_Normal, *Normal, *Edge_Vector, dist, dist2, Area,
      Twall, Tnormal, dTdn, rho_cp_solid, Prandtl_Lam, laminar_viscosity,
      thermal_diffusivity, thermal_conductivity=0.0, thermal_conductivityND,
      heat_flux_density=0.0, conductivity_over_dist=0.0;

  nDim = donor_geometry->GetnDim();

  Edge_Vector = new su2double[nDim];

  /*--- Check whether the current zone is a solid zone or a fluid zone ---*/

  bool flow = ((donor_config->GetKind_Solver() == NAVIER_STOKES)
               || (donor_config->GetKind_Solver() == RANS)
               || (donor_config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (donor_config->GetKind_Solver() == DISC_ADJ_RANS)
               || (donor_config->GetKind_Solver() == INC_NAVIER_STOKES)
               || (donor_config->GetKind_Solver() == INC_RANS)
               || (donor_config->GetKind_Solver() == DISC_ADJ_INC_NAVIER_STOKES)
               || (donor_config->GetKind_Solver() == DISC_ADJ_INC_RANS));

  bool compressible_flow    = (donor_config->GetKind_Regime() == COMPRESSIBLE) && flow;
  bool incompressible_flow  = (donor_config->GetEnergy_Equation()) && flow;
  bool heat_equation        = (donor_config->GetKind_Solver() == HEAT_EQUATION_FVM
                               || donor_config->GetKind_Solver() == DISC_ADJ_HEAT);

  Coord         = donor_geometry->node[Point_Donor]->GetCoord();

  Normal        = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNormal();
  PointNormal   = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNormal_Neighbor();
  Coord_Normal  = donor_geometry->node[PointNormal]->GetCoord();

  Twall = 0.0; Tnormal = 0.0; dTdn = 0.0; dist2 = 0.0; Area = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_Normal[iDim] - Coord[iDim];
    dist2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  dist = sqrt(dist2);
  Area = sqrt(Area);

  /*--- Retrieve temperature solution and its gradient ---*/

  if (compressible_flow) {

    Twall   = donor_solution->GetNodes()->GetPrimitive(Point_Donor,0);
    Tnormal = donor_solution->GetNodes()->GetPrimitive(PointNormal,0);

    dTdn = (Twall - Tnormal)/dist;
  }
  else if (incompressible_flow) {

    Twall   = donor_solution->GetNodes()->GetTemperature(Point_Donor);
    Tnormal = donor_solution->GetNodes()->GetTemperature(PointNormal);

    dTdn = (Twall - Tnormal)/dist;
  }
  else if (heat_equation) {
    Twall   = donor_solution->GetNodes()->GetSolution(Point_Donor,0);
    Tnormal = donor_solution->GetNodes()->GetSolution(PointNormal,0);

    // TODO: Check if these improve accuracy, if needed at all
    //    for (iDim = 0; iDim < nDim; iDim++) {
    //      dTdn += (Twall - Tnormal)/dist * (Edge_Vector[iDim]/dist) * (Normal[iDim]/Area);
    //    }

    dTdn = (Twall - Tnormal)/dist;
  }
  else {

    SU2_MPI::Error("Transfer of conjugate heat variables failed (non-supported donor solver).", CURRENT_FUNCTION);
  }

  /*--- Calculate the heat flux density (temperature gradient times thermal conductivity) ---*/

  if (compressible_flow) {

    su2double Gamma         = donor_config->GetGamma();
    su2double Gas_Constant  = donor_config->GetGas_ConstantND();
    su2double Cp            = (Gamma / (Gamma - 1.0)) * Gas_Constant;

    Prandtl_Lam             = donor_config->GetPrandtl_Lam();
    laminar_viscosity       = donor_solution->GetNodes()->GetLaminarViscosity(Point_Donor); // TDE check for consistency
    Cp                      = (Gamma / (Gamma - 1.0)) * Gas_Constant;

    thermal_conductivityND  = Cp*(laminar_viscosity/Prandtl_Lam);
    heat_flux_density       = thermal_conductivityND*dTdn;

    if (donor_config->GetCHT_Robin()) {

      thermal_conductivity    = thermal_conductivityND*donor_config->GetViscosity_Ref();
      conductivity_over_dist  = thermal_conductivity/dist;
    }
  }
  else if (incompressible_flow) {

    iPoint = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNode();

    thermal_conductivityND  = donor_solution->GetNodes()->GetThermalConductivity(iPoint);
    heat_flux_density       = thermal_conductivityND*dTdn;

    if (donor_config->GetCHT_Robin()) {

      switch (donor_config->GetKind_ConductivityModel()) {

        case CONSTANT_CONDUCTIVITY:
          thermal_conductivity = thermal_conductivityND*donor_config->GetConductivity_Ref();
          break;

        case CONSTANT_PRANDTL:
          thermal_conductivity = thermal_conductivityND*donor_config->GetGas_Constant_Ref()
                                 *donor_config->GetViscosity_Ref();
          break;
      }

      conductivity_over_dist  = thermal_conductivity/dist;
    }
  }
  else if (heat_equation) {

    /*--- Heat solver stand-alone case ---*/

    thermal_diffusivity     = donor_config->GetThermalDiffusivity_Solid();
    heat_flux_density       = thermal_diffusivity*dTdn;

    if (donor_config->GetCHT_Robin()) {
      rho_cp_solid            = donor_config->GetSpecific_Heat_Cp()*donor_config->GetDensity_Solid();
      conductivity_over_dist  = thermal_diffusivity*rho_cp_solid/dist;
    }
  }

  /*--- Set the conjugate heat variables that are transferred by default ---*/

  Donor_Variable[0] = Twall*donor_config->GetTemperature_Ref();
  Donor_Variable[1] = heat_flux_density*donor_config->GetHeat_Flux_Ref();

  /*--- We only need these for the Robin BC option ---*/

  if (donor_config->GetCHT_Robin()) {

    Donor_Variable[2] = conductivity_over_dist;
    Donor_Variable[3] = Tnormal*donor_config->GetTemperature_Ref();
  }

  if (Edge_Vector != NULL) {

    delete [] Edge_Vector;
  }
}

void CConjugateHeatInterface::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                                                 CConfig *target_config, unsigned long Marker_Target,
                                                 unsigned long Vertex_Target, unsigned long Point_Target) {

  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 0,
                                            target_config->GetRelaxation_Factor_CHT(), Target_Variable[0]);
  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 1,
                                            target_config->GetRelaxation_Factor_CHT(), Target_Variable[1]);

  if (target_config->GetCHT_Robin()) {

    target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 2,
                                              target_config->GetRelaxation_Factor_CHT(), Target_Variable[2]);
    target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 3,
                                              target_config->GetRelaxation_Factor_CHT(), Target_Variable[3]);
  }
}
