/*!
 * \file CConjugateHeatInterface.cpp
 * \brief Declaration and inlines of the class to transfer temperature and heatflux
 *        density for conjugate heat interfaces between structure and fluid zones.
 * \author O. Burghardt
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

CConjugateHeatInterface::CConjugateHeatInterface(void) : CInterface() {

}

CConjugateHeatInterface::CConjugateHeatInterface(unsigned short val_nVar, unsigned short val_nConst,
                                                 CConfig *config) : CInterface(val_nVar, val_nConst, config) {

}

CConjugateHeatInterface::~CConjugateHeatInterface(void) {

}

void CConjugateHeatInterface::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry,
                                                CConfig *donor_config, unsigned long Marker_Donor,
                                                unsigned long Vertex_Donor, unsigned long Point_Donor) {

  unsigned long iPoint;
  unsigned long PointNormal;
  unsigned short nDim, iDim;

  nDim = donor_geometry->GetnDim();

  su2double *Coord, *Coord_Normal, *Normal, *Edge_Vector, dist, dist2, Area, Twall = 0.0, Tnormal = 0.0,
      dTdn, cp_fluid, rho_cp_solid, Prandtl_Lam, laminar_viscosity,
      thermal_diffusivity, thermal_conductivity, thermal_conductivityND,
      heat_flux_density, conductivity_over_dist, Temperature_Ref;
  su2double Gamma = donor_config->GetGamma();
  su2double Gas_Constant = donor_config->GetGas_ConstantND();
  su2double Cp = (Gamma / (Gamma - 1.0)) * Gas_Constant;
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
  bool compressible_flow  = (donor_config->GetKind_Regime() == COMPRESSIBLE) && flow;
  bool incompressible_flow = (donor_config->GetEnergy_Equation()) && flow;
  bool heat_equation      = donor_config->GetKind_Solver() == HEAT_EQUATION_FVM;

  Temperature_Ref   = donor_config->GetTemperature_Ref();
  Prandtl_Lam       = donor_config->GetPrandtl_Lam();
  laminar_viscosity = donor_config->GetMu_ConstantND(); // TDE check for consistency
  cp_fluid          = donor_config->GetSpecific_Heat_Cp();
  rho_cp_solid      = donor_config->GetSpecific_Heat_Cp_Solid()*donor_config->GetDensity_Solid();

  PointNormal   = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNormal_Neighbor();
  Coord         = donor_geometry->node[Point_Donor]->GetCoord();
  Coord_Normal  = donor_geometry->node[PointNormal]->GetCoord();
  Normal        = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNormal();

  dist2 = 0.0;
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_Normal[iDim] - Coord[iDim];
    dist2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  dist = sqrt(dist2);
  Area = sqrt(Area);

  /*--- Retrieve temperature solution (later set is as the first donor variable) and its gradient ---*/

  dTdn = 0.0;

  if (compressible_flow) {

    Twall   = donor_solution->node[Point_Donor]->GetPrimitive(0)*Temperature_Ref;
    Tnormal = donor_solution->node[PointNormal]->GetPrimitive(0)*Temperature_Ref;

    dTdn = (Twall - Tnormal)/dist;
  }
  else if (incompressible_flow) {

    Twall   = donor_solution->node[Point_Donor]->GetTemperature()*Temperature_Ref;
    Tnormal = donor_solution->node[PointNormal]->GetTemperature()*Temperature_Ref;

    dTdn = (Twall - Tnormal)/dist;
  }
  else if (flow || heat_equation) {
    Twall   = donor_solution->node[Point_Donor]->GetSolution(0)*Temperature_Ref;
    Tnormal = donor_solution->node[PointNormal]->GetSolution(0)*Temperature_Ref;

    //    for (iDim = 0; iDim < nDim; iDim++) {
    //      dTdn += (Twall - Tnormal)/dist * (Edge_Vector[iDim]/dist) * (Normal[iDim]/Area);
    //    }

    dTdn = (Twall - Tnormal)/dist;
  }
  else {
    cout << "WARNING: Transfer of conjugate heat variables is called with non-supported donor solver!" << endl;
  }

  /*--- Calculate the heat flux density (temperature gradient times thermal conductivity)
   *--- and set it as second donor variable ---*/
  if (compressible_flow) {

    iPoint = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNode();

    thermal_conductivityND  = Cp*(laminar_viscosity/Prandtl_Lam);
    thermal_conductivity    = thermal_conductivityND*donor_config->GetViscosity_Ref();

    heat_flux_density       = thermal_conductivity*dTdn;
    conductivity_over_dist  = thermal_conductivity/dist;
  }
  else if (incompressible_flow) {

    iPoint = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNode();

    thermal_conductivityND  = donor_solution->node[iPoint]->GetThermalConductivity();
    thermal_conductivity = thermal_conductivityND*donor_config->GetConductivity_Ref();

    switch (donor_config->GetKind_ConductivityModel()) {

      case CONSTANT_CONDUCTIVITY:
        thermal_conductivity = thermal_conductivityND*donor_config->GetConductivity_Ref();
        break;

      case CONSTANT_PRANDTL:
        thermal_conductivity = thermal_conductivityND*donor_config->GetGas_Constant_Ref()
                               *donor_config->GetViscosity_Ref();
        break;
    }

    heat_flux_density       = thermal_conductivity*dTdn;
    conductivity_over_dist  = thermal_conductivity/dist;
  }
  else if (flow) {

    iPoint = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNode();

    thermal_conductivityND  = laminar_viscosity/Prandtl_Lam;
    thermal_conductivity    = thermal_conductivityND*donor_config->GetViscosity_Ref()*cp_fluid;

    heat_flux_density       = thermal_conductivity*dTdn;
    conductivity_over_dist  = thermal_conductivity/dist;
  }
  else {

    thermal_diffusivity     = donor_config->GetThermalDiffusivity_Solid();
    heat_flux_density       = (thermal_diffusivity*dTdn)*rho_cp_solid;
    conductivity_over_dist  = thermal_diffusivity*rho_cp_solid/dist;
  }

  Donor_Variable[0] = Twall;
  Donor_Variable[1] = heat_flux_density;
  Donor_Variable[2] = conductivity_over_dist;
  Donor_Variable[3] = Tnormal;

  delete [] Edge_Vector;
}

void CConjugateHeatInterface::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                                                 CConfig *target_config, unsigned long Marker_Target,
                                                 unsigned long Vertex_Target, unsigned long Point_Target) {

  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 0,
                                            target_config->GetRelaxation_Factor_CHT(), Target_Variable[0]);
  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 1,
                                            target_config->GetRelaxation_Factor_CHT(), Target_Variable[1]);
  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 2,
                                            target_config->GetRelaxation_Factor_CHT(), Target_Variable[2]);
  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 3,
                                            target_config->GetRelaxation_Factor_CHT(), Target_Variable[3]);
}
