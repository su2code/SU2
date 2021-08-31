/*!
 * \file CTurbSolver.cpp
 * \brief Main subrotuines of CTurbSolver class
 * \author F. Palacios, A. Bueno
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/solvers/CTurbSolver.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"


CTurbSolver::CTurbSolver(bool conservative) : CScalarSolver(conservative) { }

CTurbSolver::CTurbSolver(CGeometry* geometry, CConfig *config, bool conservative) : CScalarSolver(geometry, config, conservative) { }

CTurbSolver::~CTurbSolver() {

  for (auto& mat : SlidingState) {
    for (auto ptr : mat) delete [] ptr;
  }

}

void CTurbSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Riemann(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT: case STATIC_SUPERSONIC_INFLOW_PT: case STATIC_SUPERSONIC_INFLOW_PD: case DENSITY_VELOCITY:
    BC_Inlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  case STATIC_PRESSURE:
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}

void CTurbSolver::BC_TurboRiemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Riemann(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT: case STATIC_SUPERSONIC_INFLOW_PT: case STATIC_SUPERSONIC_INFLOW_PD: case DENSITY_VELOCITY:
    BC_Inlet_Turbo(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  case STATIC_PRESSURE:
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}


void CTurbSolver::BC_Giles(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Giles(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT:case TOTAL_CONDITIONS_PT_1D: case DENSITY_VELOCITY:
    BC_Inlet_Turbo(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  case MIXING_IN:
    if (config->GetBoolTurbMixingPlane()){
      BC_Inlet_MixingPlane(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    }
    else{
      BC_Inlet_Turbo(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    }
    break;

  case STATIC_PRESSURE: case MIXING_OUT: case STATIC_PRESSURE_1D: case RADIAL_EQUILIBRIUM:
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}

void CTurbSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config) {

  const bool sst = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
  const auto nPrimVar = solver_container[FLOW_SOL]->GetnPrimVar();
  su2double *PrimVar_j = new su2double[nPrimVar];
  su2double solution_j[MAXNVAR] = {0.0};

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) != FLUID_INTERFACE) continue;

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->nVertex[iMarker]; iVertex++) {

      const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

      if (!geometry->nodes->GetDomain(iPoint)) continue;

      const auto Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
      const auto nDonorVertex = GetnSlidingStates(iMarker,iVertex);

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[iMarker][iVertex]->GetNormal()[iDim];

      su2double* PrimVar_i = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      auto Jacobian_i = Jacobian.GetBlock(iPoint,iPoint);

      /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/

      for (auto jVertex = 0; jVertex < nDonorVertex; jVertex++) {

        for (auto iVar = 0u; iVar < nPrimVar; iVar++)
          PrimVar_j[iVar] = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, iVar, jVertex);

        /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/

        const su2double weight = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

        /*--- Set primitive variables ---*/

        conv_numerics->SetPrimitive( PrimVar_i, PrimVar_j );

        /*--- Set the turbulent variable states ---*/

        for (auto iVar = 0u; iVar < nVar; ++iVar)
          solution_j[iVar] = GetSlidingState(iMarker, iVertex, iVar, jVertex);

        conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);

        /*--- Set the normal vector ---*/

        conv_numerics->SetNormal(Normal);

        if (dynamic_grid)
          conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

        auto residual = conv_numerics->ComputeResidual(config);

        /*--- Accumulate the residuals to compute the average ---*/

        for (auto iVar = 0u; iVar < nVar; iVar++) {
          LinSysRes(iPoint,iVar) += weight*residual[iVar];
          for (auto jVar = 0u; jVar < nVar; jVar++)
            Jacobian_i[iVar*nVar+jVar] += SU2_TYPE::GetValue(weight*residual.jacobian_i[iVar][jVar]);
        }
      }

      /*--- Set the normal vector and the coordinates ---*/

      visc_numerics->SetNormal(Normal);
      su2double Coord_Reflected[MAXNDIM];
      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);

      /*--- Primitive variables ---*/

      visc_numerics->SetPrimitive(PrimVar_i, PrimVar_j);

      /*--- Turbulent variables and their gradients ---*/

      visc_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);
      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/

      if(sst) visc_numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(iPoint));

      /*--- Compute and update residual ---*/

      auto residual = visc_numerics->ComputeResidual(config);

      LinSysRes.SubtractBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

  delete [] PrimVar_j;

}

void CTurbSolver::Impose_Fixed_Values(const CGeometry *geometry, const CConfig *config){

  /*--- Check whether turbulence quantities are fixed to far-field values on a half-plane. ---*/
  if(config->GetTurb_Fixed_Values()){

    const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

    /*--- Form normalized far-field velocity ---*/
    const su2double* velocity_inf = config->GetVelocity_FreeStreamND();
    su2double velmag_inf = GeometryToolbox::Norm(nDim, velocity_inf);
    if(velmag_inf==0)
      SU2_MPI::Error("Far-field velocity is zero, cannot fix turbulence quantities to inflow values.", CURRENT_FUNCTION);
    su2double unit_velocity_inf[MAXNDIM];
    for(unsigned short iDim=0; iDim<nDim; iDim++)
      unit_velocity_inf[iDim] = velocity_inf[iDim] / velmag_inf;

    SU2_OMP_FOR_DYN(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if( GeometryToolbox::DotProduct(nDim, geometry->nodes->GetCoord(iPoint), unit_velocity_inf)
        < config->GetTurb_Fixed_Values_MaxScalarProd() ) {
        /*--- Set the solution values and zero the residual ---*/
        nodes->SetSolution_Old(iPoint, Solution_Inf);
        nodes->SetSolution(iPoint, Solution_Inf);
        LinSysRes.SetBlock_Zero(iPoint);
        if (implicit) {
          /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
          for(unsigned long iVar=0; iVar<nVar; iVar++)
            Jacobian.DeleteValsRowi(iPoint*nVar+iVar);
        }
      }
    }
    END_SU2_OMP_FOR
  }

}
