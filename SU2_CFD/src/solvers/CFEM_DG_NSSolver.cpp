/*!
 * \file CFEM_DG_NSSolver.cpp
 * \brief Main subroutines for solving finite element Navier-Stokes flow problems
 * \author J. Alonso, E. van der Weide, T. Economon
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


#include "../../include/solvers/CFEM_DG_NSSolver.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"

enum {
SIZE_ARR_NORM = 8
};

CFEM_DG_NSSolver::CFEM_DG_NSSolver() : CFEM_DG_EulerSolver() {

  /*--- Basic array initialization ---*/
  CD_Visc  = nullptr; CL_Visc  = nullptr; CSF_Visc = nullptr; CEff_Visc = nullptr;
  CMx_Visc = nullptr; CMy_Visc = nullptr; CMz_Visc = nullptr;
  CFx_Visc = nullptr; CFy_Visc = nullptr; CFz_Visc = nullptr;

  /*--- Surface-based array initialization ---*/
  Surface_CL_Visc  = nullptr; Surface_CD_Visc  = nullptr; Surface_CSF_Visc = nullptr; Surface_CEff_Visc = nullptr;
  Surface_CFx_Visc = nullptr; Surface_CFy_Visc = nullptr; Surface_CFz_Visc = nullptr;
  Surface_CMx_Visc = nullptr; Surface_CMy_Visc = nullptr; Surface_CMz_Visc = nullptr;
  MaxHeatFlux_Visc = nullptr; Heat_Visc = nullptr;

  /*--- Set the SGS model to NULL and indicate that no SGS model is used. ---*/
  SGSModel     = nullptr;
  SGSModelUsed = false;
}

CFEM_DG_NSSolver::CFEM_DG_NSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
 : CFEM_DG_EulerSolver(geometry, config, iMesh) {

  /*--- Array initialization ---*/
  CD_Visc = nullptr;  CL_Visc = nullptr;  CSF_Visc = nullptr; CEff_Visc = nullptr;
  CMx_Visc = nullptr; CMy_Visc = nullptr; CMz_Visc = nullptr;
  CFx_Visc = nullptr; CFy_Visc = nullptr; CFz_Visc = nullptr;

  Surface_CL_Visc  = nullptr; Surface_CD_Visc = nullptr;  Surface_CSF_Visc = nullptr; Surface_CEff_Visc = nullptr;
  Surface_CFx_Visc = nullptr; Surface_CFy_Visc = nullptr; Surface_CFz_Visc = nullptr;
  Surface_CMx_Visc = nullptr; Surface_CMy_Visc = nullptr; Surface_CMz_Visc = nullptr;
  MaxHeatFlux_Visc = nullptr; Heat_Visc = nullptr;

  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/

  //LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  //LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Non dimensional coefficients ---*/
  CD_Visc       = new su2double[nMarker];
  CL_Visc       = new su2double[nMarker];
  CSF_Visc      = new su2double[nMarker];
  CMx_Visc      = new su2double[nMarker];
  CMy_Visc      = new su2double[nMarker];
  CMz_Visc      = new su2double[nMarker];
  CEff_Visc     = new su2double[nMarker];
  CFx_Visc      = new su2double[nMarker];
  CFy_Visc      = new su2double[nMarker];
  CFz_Visc      = new su2double[nMarker];

  Surface_CL_Visc   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Visc   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Visc = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Visc  = new su2double[config->GetnMarker_Monitoring()];

  Heat_Visc        = new su2double[nMarker];
  MaxHeatFlux_Visc = new su2double[nMarker];

  /*--- Init total coefficients ---*/

  Total_CD   = 0.0; Total_CL  = 0.0; Total_CSF = 0.0;
  Total_CMx  = 0.0; Total_CMy = 0.0; Total_CMz = 0.0;
  Total_CEff = 0.0;
  Total_CFx  = 0.0; Total_CFy = 0.0; Total_CFz = 0.0;

  /*--- Read farfield conditions from config ---*/

  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  Prandtl_Lam   = config->GetPrandtl_Lam();
  Prandtl_Turb  = config->GetPrandtl_Turb();
  Tke_Inf       = config->GetTke_FreeStreamND();

  /*--- Set the SGS model in case an LES simulation is carried out ---*/

  if(config->GetKind_Solver() == MAIN_SOLVER::FEM_LES) {

    /* Make a distinction between the SGS models used and set SGSModel and
       SGSModelUsed accordingly. */
    switch( config->GetKind_SGS_Model() ) {

      case TURB_SGS_MODEL::IMPLICIT_LES:
        SGSModel     = nullptr;
        SGSModelUsed = false;
        break;

      case TURB_SGS_MODEL::SMAGORINSKY:
        SGSModel     = new CSmagorinskyModel;
        SGSModelUsed = true;
        break;

      case TURB_SGS_MODEL::WALE:
        SGSModel     = new CWALEModel;
        SGSModelUsed = true;
        break;

      case TURB_SGS_MODEL::VREMAN:
        SGSModel     = new CVremanModel;
        SGSModelUsed = true;
        break;

      default:
        SU2_MPI::Error("Unknown SGS model encountered", CURRENT_FUNCTION);
    }
  }
  else {

    /* No LES, so no SGS model needed.
       Set the pointer to NULL and the boolean to false. */
    SGSModel     = nullptr;
    SGSModelUsed = false;
  }
}

CFEM_DG_NSSolver::~CFEM_DG_NSSolver() {

        delete [] CD_Visc;
        delete [] CL_Visc;
       delete [] CSF_Visc;
       delete [] CMx_Visc;
       delete [] CMy_Visc;
       delete [] CMz_Visc;
       delete [] CFx_Visc;
       delete [] CFy_Visc;
       delete [] CFz_Visc;
      delete [] CEff_Visc;

    delete [] Surface_CL_Visc;
    delete [] Surface_CD_Visc;
   delete [] Surface_CSF_Visc;
  delete [] Surface_CEff_Visc;
   delete [] Surface_CFx_Visc;
   delete [] Surface_CFy_Visc;
   delete [] Surface_CFz_Visc;
   delete [] Surface_CMx_Visc;
   delete [] Surface_CMy_Visc;
   delete [] Surface_CMz_Visc;

   delete [] Heat_Visc;
   delete [] MaxHeatFlux_Visc;

  delete SGSModel;
}

void CFEM_DG_NSSolver::Friction_Forces(const CGeometry* geometry, const CConfig* config) {

  /* Allocate the memory for the work array and initialize it to zero to avoid
     warnings in debug mode  about uninitialized memory when padding is applied. */
  vector<su2double> workArrayVec(sizeWorkArray, 0.0);
  su2double *workArray = workArrayVec.data();

  /*--------------------------------------------------------------------------*/
  /*--- The skin friction is computed using the laminar viscosity and      ---*/
  /*--- velocity gradients. This is correct when integration to the wall   ---*/
  /*--- is performed, but not when wall functions are used. Hence, this    ---*/
  /*--- function must be modified when wall functions are implemented.     ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /* Constant factor present in the heat flux vector. */
  const su2double factHeatFlux_Lam = Gamma/Prandtl_Lam;

  /*--- Get the information of the angle of attack, reference area, etc. ---*/
  const su2double Alpha        = config->GetAoA()*PI_NUMBER/180.0;
  const su2double Beta         = config->GetAoS()*PI_NUMBER/180.0;
  const su2double RefLength    = config->GetRefLength();
  auto Origin      = config->GetRefOriginMoment(0);

  /*--- Evaluate reference values for non-dimensionalization. ---*/
  const su2double RefHeatFlux = config->GetHeat_Flux_Ref();

  const su2double factor = 1.0 / AeroCoeffForceRef;

  /*--- Variables initialization ---*/
  AllBound_CD_Visc = 0.0;  AllBound_CL_Visc  = 0.0; AllBound_CSF_Visc = 0.0;
  AllBound_CMx_Visc = 0.0; AllBound_CMy_Visc = 0.0; AllBound_CMz_Visc = 0.0;
  AllBound_CFx_Visc = 0.0; AllBound_CFy_Visc = 0.0; AllBound_CFz_Visc = 0.0;

  AllBound_HeatFlux_Visc = 0.0; AllBound_MaxHeatFlux_Visc = 0.0; AllBound_CEff_Visc = 0.0;

  for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring(); ++iMarker_Monitoring) {
    Surface_CL_Visc[iMarker_Monitoring]  = 0.0; Surface_CD_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CSF_Visc[iMarker_Monitoring] = 0.0; Surface_CEff_Visc[iMarker_Monitoring] = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring] = 0.0; Surface_CFy_Visc[iMarker_Monitoring]  = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring] = 0.0; Surface_CMx_Visc[iMarker_Monitoring]  = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring] = 0.0; Surface_CMz_Visc[iMarker_Monitoring]  = 0.0;
  }

  /*--- Loop over the Navier-Stokes markers ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    /* Check if this boundary must be monitored. */
    const unsigned short Monitoring = config->GetMarker_All_Monitoring(iMarker);
    if(Monitoring == YES) {

      /* Easier storage of the boundary condition. */
      const unsigned short Boundary = config->GetMarker_All_KindBC(iMarker);

      /*--- Obtain the origin for the moment computation for a particular marker ---*/
      for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                       ++iMarker_Monitoring) {
        string Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        string Marker_Tag     = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }

      /* Check for a boundary for which the viscous forces must be computed. */
      if((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {

        /*--- Determine the prescribed heat flux or prescribed temperature. ---*/
        bool HeatFlux_Prescribed = false, Temperature_Prescribed = false;
        su2double Wall_HeatFlux = 0.0, Wall_Temperature = 0.0;

        const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if(Boundary == HEAT_FLUX) {
          HeatFlux_Prescribed = true;
          Wall_HeatFlux       = config->GetWall_HeatFlux(Marker_Tag);
        }
        else {
          Temperature_Prescribed = true;
          Wall_Temperature       = config->GetIsothermal_Temperature(Marker_Tag)
                                 / config->GetTemperature_Ref();
        }

        /*--- Forces initialization at each Marker ---*/
        CD_Visc[iMarker]  = 0.0; CL_Visc[iMarker]  = 0.0; CSF_Visc[iMarker] = 0.0;
        CMx_Visc[iMarker] = 0.0; CMy_Visc[iMarker] = 0.0; CMz_Visc[iMarker] = 0.0;
        CFx_Visc[iMarker] = 0.0; CFy_Visc[iMarker] = 0.0; CFz_Visc[iMarker] = 0.0;

        Heat_Visc[iMarker]  = 0.0; MaxHeatFlux_Visc[iMarker] = 0.0; CEff_Visc[iMarker] = 0.0;

        su2double ForceViscous[]  = {0.0, 0.0, 0.0};
        su2double MomentViscous[] = {0.0, 0.0, 0.0};

        /* Easier storage of the boundary faces for this boundary marker. */
        const unsigned long      nSurfElem = boundaries[iMarker].surfElem.size();
        const CSurfaceElementFEM *surfElem = boundaries[iMarker].surfElem.data();

        /* Check if a wall treatment is used. */
        if( boundaries[iMarker].wallModel ) {

          /*--- Wall treatment is used, so the wall shear stress and heat flux
                are computed using the wall model. As the interpolation of data
                of the exchange point is different for each element, it is not
                possible to treat multiple faces simultaneously. So here just
                a loop over the number of faces is carried out. ---*/
          for(unsigned long l=0; l<nSurfElem; ++l) {

            /* Get the required information from the corresponding standard face. */
            const unsigned short ind  = surfElem[l].indStandardElement;
            const su2double *weights  = standardBoundaryFacesSol[ind].GetWeightsIntegration();

            /* Loop over the donors for this boundary face. */
            for(unsigned long j=0; j<surfElem[l].donorsWallFunction.size(); ++j) {

              /* Easier storage of the element ID of the donor and set the pointer
                 where the solution of this element starts. Note that VecWorkSolDOFs
                 must be used and not VecSolDOFs, because it is possible that a donor
                 is a halo element, which are not stored in VecSolDOFs. */
              const unsigned long  donorID   = surfElem[l].donorsWallFunction[j];
              const unsigned short timeLevel = volElem[donorID].timeLevel;
              const unsigned short nDOFsElem = volElem[donorID].nDOFsSol;
              const su2double *solDOFsElem   = VecWorkSolDOFs[timeLevel].data()
                                             + nVar*volElem[donorID].offsetDOFsSolThisTimeLevel;

              /* Determine the number of integration points for this donor and
                 interpolate the solution for the corresponding exchange points. */
              const unsigned short nIntThisDonor = surfElem[l].nIntPerWallFunctionDonor[j+1]
                                                 - surfElem[l].nIntPerWallFunctionDonor[j];

              blasFunctions->gemm(nIntThisDonor, nVar, nDOFsElem, surfElem[l].matWallFunctionDonor[j].data(),
                                  solDOFsElem, workArray, config);

              /* Loop over the integration points for this donor element. */
              for(unsigned short i=surfElem[l].nIntPerWallFunctionDonor[j];
                                 i<surfElem[l].nIntPerWallFunctionDonor[j+1]; ++i) {

                /* Easier storage of the actual integration point. */
                const unsigned short ii = surfElem[l].intPerWallFunctionDonor[i];

                /* Determine the normal, wall velocity and coordinates
                   for this integration point. */
                const su2double *normals = surfElem[l].metricNormalsFace.data() + ii*(nDim+1);
                const su2double *gridVel = surfElem[l].gridVelocities.data() + ii*nDim;
                const su2double *Coord   = surfElem[l].coorIntegrationPoints.data() + ii*nDim;

                /* Determine the velocities and pressure in the exchange point. */
                const su2double *solInt = workArray
                                        + nVar*(i-surfElem[l].nIntPerWallFunctionDonor[j]);

                su2double rhoInv = 1.0/solInt[0];
                su2double vel[]  = {0.0, 0.0, 0.0};
                for(unsigned short k=0; k<nDim; ++k) vel[k] = rhoInv*solInt[k+1];

                su2double vel2Mag = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
                su2double eInt    = rhoInv*solInt[nVar-1] - 0.5*vel2Mag;

                FluidModel->SetTDState_rhoe(solInt[0], eInt);
                const su2double Pressure = FluidModel->GetPressure();
                const su2double Temperature = FluidModel->GetTemperature();
                const su2double LaminarViscosity= FluidModel->GetLaminarViscosity();

                /* Subtract the prescribed wall velocity, i.e. grid velocity
                   from the velocity in the exchange point. */
                for(unsigned short k=0; k<nDim; ++k) vel[k] -= gridVel[k];

                /* Determine the tangential velocity by subtracting the normal
                   velocity component. */
                su2double velNorm = 0.0;
                for(unsigned short k=0; k<nDim; ++k) velNorm += normals[k]*vel[k];
                for(unsigned short k=0; k<nDim; ++k) vel[k]  -= normals[k]*velNorm;

                /* Determine the magnitude of the tangential velocity as well
                   as its direction (unit vector). */
                su2double velTan = sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
                velTan = max(velTan,1.e-25);

                su2double dirTan[] = {0.0, 0.0, 0.0};
                for(unsigned short k=0; k<nDim; ++k) dirTan[k] = vel[k]/velTan;

                /* Compute the wall shear stress and heat flux vector using
                   the wall model. */
                su2double tauWall, qWall, ViscosityWall, kOverCvWall;

                boundaries[iMarker].wallModel->WallShearStressAndHeatFlux(Temperature, velTan,
                                                                          LaminarViscosity, Pressure,
                                                                          Wall_HeatFlux, HeatFlux_Prescribed,
                                                                          Wall_Temperature, Temperature_Prescribed,
                                                                          FluidModel, tauWall, qWall,
                                                                          ViscosityWall, kOverCvWall);

                /* Update the viscous forces and moments. Note that the force direction
                   is the direction of the tangential velocity. */
                const su2double dForceMag = tauWall*weights[ii]*normals[nDim]*factor;
                su2double dForces[] = {0.0, 0.0, 0.0}, dCoor[] = {0.0, 0.0, 0.0};
                for(unsigned short k=0; k<nDim; ++k) {
                  dCoor[k]         = Coord[k] - Origin[k];
                  dForces[k]       = dForceMag*dirTan[k];
                  ForceViscous[k] += dForces[k];
                }

                if(nDim == 2) {
                  MomentViscous[2] += (dForces[1]*dCoor[0] - dForces[0]*dCoor[1])/RefLength;
                }
                else {
                  MomentViscous[0] += (dForces[2]*dCoor[1] - dForces[1]*dCoor[2])/RefLength;
                  MomentViscous[1] += (dForces[0]*dCoor[2] - dForces[2]*dCoor[0])/RefLength;
                  MomentViscous[2] += (dForces[1]*dCoor[0] - dForces[0]*dCoor[1])/RefLength;
                }

                /* Update the heat flux and maximum heat flux for this marker. */
                Heat_Visc[iMarker]       += qWall*weights[i]*normals[nDim];
                MaxHeatFlux_Visc[iMarker] = max(MaxHeatFlux_Visc[iMarker], fabs(qWall));
              }
            }
          }

        } else {

          /*--- Integration to the wall is used so the wall data must be computed.
                Loop over the faces of this boundary. Multiple faces are treated
                simultaneously to improve the performance of the matrix
                multiplications. As a consequence, the update of the counter l
                happens at the end of this loop section. ---*/
          for(unsigned long l=0; l<nSurfElem;) {

            /* Determine the end index for this chunk of faces and the padded
               N value in the gemm computations. */
            unsigned long lEnd;
            unsigned short ind, llEnd, NPad;

            MetaDataChunkOfElem(surfElem, l, nSurfElem, nFaceSimul,
                                nPadMin, lEnd, ind, llEnd, NPad);

            /* Get the required information from the corresponding standard face. */
            const unsigned short nInt      = standardBoundaryFacesSol[ind].GetNIntegration();
            const unsigned short nDOFsFace = standardBoundaryFacesSol[ind].GetNDOFsFace();
            const unsigned short nDOFsElem = standardBoundaryFacesSol[ind].GetNDOFsElem();
            const su2double *basisFace     = standardBoundaryFacesSol[ind].GetBasisFaceIntegration();
            const su2double *derBasisElem  = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();
            const su2double *weights       = standardBoundaryFacesSol[ind].GetWeightsIntegration();

            /* Set the pointers for the local work arrays. */
            su2double *solInt     = workArray;
            su2double *gradSolInt = solInt     + NPad*nInt;
            su2double *solCopy    = gradSolInt + NPad*nInt*nDim;

            /* Loop over the faces that are treated simultaneously. */
            for(unsigned short ll=0; ll<llEnd; ++ll) {
              const unsigned short llNVar = ll*nVar;

              /* Easier storage of the DOFs of the face. */
              const unsigned long *DOFs = surfElem[l+ll].DOFsSolFace.data();

              /* Copy the solution of the DOFs of the face such that it is
                 contiguous in memory. */
              for(unsigned short i=0; i<nDOFsFace; ++i) {
                const su2double *solDOF = VecSolDOFs.data() + nVar*DOFs[i];
                su2double       *sol    = solCopy + NPad*i + llNVar;
                for(unsigned short mm=0; mm<nVar; ++mm)
                  sol[mm] = solDOF[mm];
              }
            }

            /* Call the general function to carry out the matrix product to determine
               the solution in the integration points. */
            blasFunctions->gemm(nInt, NPad, nDOFsFace, basisFace, solCopy, solInt, config);

            /*--- Store the solution of the DOFs of the adjacent elements in contiguous
                  memory such that the function blasFunctions->gemm can be used to compute
                  the gradients solution variables in the integration points of the face. ---*/
            for(unsigned short ll=0; ll<llEnd; ++ll) {
              const unsigned short llNVar = ll*nVar;
              const unsigned long  lll    = l + ll;

              for(unsigned short i=0; i<nDOFsElem; ++i) {
                const su2double *solDOF = VecSolDOFs.data() + nVar*surfElem[lll].DOFsSolElement[i];
                su2double       *sol    = solCopy + NPad*i + llNVar;
                for(unsigned short mm=0; mm<nVar; ++mm)
                  sol[mm] = solDOF[mm];
              }
            }

            /* Compute the gradients in the integration points. Call the general function to
               carry out the matrix product. */
            blasFunctions->gemm(nInt*nDim, NPad, nDOFsElem, derBasisElem, solCopy, gradSolInt, config);

            /* Determine the offset between r- and -s-derivatives, which is also the
               offset between s- and t-derivatives. */
            const unsigned short offDeriv = NPad*nInt;

            /* Make a distinction between two and three space dimensions
               in order to have the most efficient code. */
            switch( nDim ) {

              case 2: {

                /* Two dimensional simulation. Loop over the number of faces treated
                   simultaneously. */
                for(unsigned short ll=0; ll<llEnd; ++ll) {
                  const unsigned short llNVar = ll*nVar;
                  const unsigned long  lll    = l + ll;

                  /* Loop over the integration points of this surface element. */
                  for(unsigned short i=0; i<nInt; ++i) {

                    /* Easier storage of the solution, its gradients, the normals,
                       the metric terms and the coordinates of this integration point. */
                    const su2double *sol         = solInt     + i*NPad + llNVar;
                    const su2double *solDOFDr    = gradSolInt + i*NPad + llNVar;
                    const su2double *solDOFDs    = solDOFDr   + offDeriv;
                    const su2double *normals     = surfElem[lll].metricNormalsFace.data()
                                                 + i*(nDim+1);
                    const su2double *metricTerms = surfElem[lll].metricCoorDerivFace.data()
                                                 + i*nDim*nDim;
                    const su2double *Coord       = surfElem[lll].coorIntegrationPoints.data()
                                                 + i*nDim;

                    /* Easier storage of the metric terms. */
                    const su2double drdx = metricTerms[0];
                    const su2double drdy = metricTerms[1];
                    const su2double dsdx = metricTerms[2];
                    const su2double dsdy = metricTerms[3];

                    /* Compute the Cartesian gradients of the solution. */
                    const su2double drhodx = solDOFDr[0]*drdx + solDOFDs[0]*dsdx;
                    const su2double drudx  = solDOFDr[1]*drdx + solDOFDs[1]*dsdx;
                    const su2double drvdx  = solDOFDr[2]*drdx + solDOFDs[2]*dsdx;
                    const su2double drEdx  = solDOFDr[3]*drdx + solDOFDs[3]*dsdx;

                    const su2double drhody = solDOFDr[0]*drdy + solDOFDs[0]*dsdy;
                    const su2double drudy  = solDOFDr[1]*drdy + solDOFDs[1]*dsdy;
                    const su2double drvdy  = solDOFDr[2]*drdy + solDOFDs[2]*dsdy;
                    const su2double drEdy  = solDOFDr[3]*drdy + solDOFDs[3]*dsdy;

                    /* Compute the velocities and static energy in this
                       integration point. */
                    const su2double DensityInv   = 1.0/sol[0];
                    const su2double u            = DensityInv*sol[1];
                    const su2double v            = DensityInv*sol[2];
                    const su2double TotalEnergy  = DensityInv*sol[3];
                    const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v);

                    /* Compute the Cartesian gradients of the velocities and
                       static energy in this integration point and also the
                       divergence of the velocity. */
                    const su2double dudx = DensityInv*(drudx - u*drhodx);
                    const su2double dudy = DensityInv*(drudy - u*drhody);

                    const su2double dvdx = DensityInv*(drvdx - v*drhodx);
                    const su2double dvdy = DensityInv*(drvdy - v*drhody);

                    const su2double dStaticEnergydx = DensityInv*(drEdx - TotalEnergy*drhodx)
                                                    - u*dudx - v*dvdx;
                    const su2double dStaticEnergydy = DensityInv*(drEdy - TotalEnergy*drhody)
                                                    - u*dudy - v*dvdy;
                    const su2double divVel = dudx + dvdy;

                    /* Compute the laminar viscosity. */
                    FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
                    const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

                    /* Set the value of the second viscosity and compute the
                       divergence term in the viscous normal stresses. */
                    const su2double lambda     = -TWO3*ViscosityLam;
                    const su2double lamDivTerm =  lambda*divVel;

                    /* Compute the viscous stress tensor and the normal flux.
                       Note that there is a plus sign for the heat flux, because
                       the normal points into the geometry. */
                    const su2double tauxx = 2.0*ViscosityLam*dudx + lamDivTerm;
                    const su2double tauyy = 2.0*ViscosityLam*dvdy + lamDivTerm;
                    const su2double tauxy = ViscosityLam*(dudy + dvdx);

                    const su2double qHeatNorm = ViscosityLam*factHeatFlux_Lam
                                              * (dStaticEnergydx*normals[0]
                                              +  dStaticEnergydy*normals[1]);

                    /* Update the viscous force and moment. Note that the normal
                       points into the geometry, hence the minus sign for the stress. */
                    const su2double scaleFac =  weights[i]*normals[nDim]*factor;
                    const su2double Fx       = -scaleFac*(tauxx*normals[0] + tauxy*normals[1]);
                    const su2double Fy       = -scaleFac*(tauxy*normals[0] + tauyy*normals[1]);

                    ForceViscous[0] += Fx;
                    ForceViscous[1] += Fy;

                    const su2double dx = Coord[0] - Origin[0];
                    const su2double dy = Coord[1] - Origin[1];

                    MomentViscous[2] += (Fy*dx - Fx*dy)/RefLength;

                    /* Update the heat flux and maximum heat flux for this marker. */
                    Heat_Visc[iMarker] += qHeatNorm*weights[i]*normals[nDim]*RefHeatFlux;
                    MaxHeatFlux_Visc[iMarker] = max(MaxHeatFlux_Visc[iMarker], fabs(qHeatNorm));
                  }
                }

                break;
              }

              /*------------------------------------------------------------------*/

              case 3: {

                /* Three dimensional simulation. Loop over the number of faces treated
                   simultaneously. */
                for(unsigned short ll=0; ll<llEnd; ++ll) {
                  const unsigned short llNVar = ll*nVar;
                  const unsigned long  lll    = l + ll;

                  /* Loop over the integration points of this surface element. */
                  for(unsigned short i=0; i<nInt; ++i) {

                    /* Easier storage of the solution, its gradients, the normals,
                       the metric terms and the coordinates of this integration point. */
                    const su2double *sol         = solInt     + i*NPad + llNVar;
                    const su2double *solDOFDr    = gradSolInt + i*NPad + llNVar;
                    const su2double *solDOFDs    = solDOFDr   + offDeriv;
                    const su2double *solDOFDt    = solDOFDs   + offDeriv;
                    const su2double *normals     = surfElem[lll].metricNormalsFace.data()
                                                 + i*(nDim+1);
                    const su2double *metricTerms = surfElem[lll].metricCoorDerivFace.data()
                                                 + i*nDim*nDim;
                    const su2double *Coord       = surfElem[lll].coorIntegrationPoints.data()
                                                 + i*nDim;

                    /* Easier storage of the metric terms. */
                    const su2double drdx = metricTerms[0];
                    const su2double drdy = metricTerms[1];
                    const su2double drdz = metricTerms[2];

                    const su2double dsdx = metricTerms[3];
                    const su2double dsdy = metricTerms[4];
                    const su2double dsdz = metricTerms[5];

                    const su2double dtdx = metricTerms[6];
                    const su2double dtdy = metricTerms[7];
                    const su2double dtdz = metricTerms[8];

                    /* Compute the Cartesian gradients of the solution. */
                    const su2double drhodx = solDOFDr[0]*drdx + solDOFDs[0]*dsdx + solDOFDt[0]*dtdx;
                    const su2double drudx  = solDOFDr[1]*drdx + solDOFDs[1]*dsdx + solDOFDt[1]*dtdx;
                    const su2double drvdx  = solDOFDr[2]*drdx + solDOFDs[2]*dsdx + solDOFDt[2]*dtdx;
                    const su2double drwdx  = solDOFDr[3]*drdx + solDOFDs[3]*dsdx + solDOFDt[3]*dtdx;
                    const su2double drEdx  = solDOFDr[4]*drdx + solDOFDs[4]*dsdx + solDOFDt[4]*dtdx;

                    const su2double drhody = solDOFDr[0]*drdy + solDOFDs[0]*dsdy + solDOFDt[0]*dtdy;
                    const su2double drudy  = solDOFDr[1]*drdy + solDOFDs[1]*dsdy + solDOFDt[1]*dtdy;
                    const su2double drvdy  = solDOFDr[2]*drdy + solDOFDs[2]*dsdy + solDOFDt[2]*dtdy;
                    const su2double drwdy  = solDOFDr[3]*drdy + solDOFDs[3]*dsdy + solDOFDt[3]*dtdy;
                    const su2double drEdy  = solDOFDr[4]*drdy + solDOFDs[4]*dsdy + solDOFDt[4]*dtdy;

                    const su2double drhodz = solDOFDr[0]*drdz + solDOFDs[0]*dsdz + solDOFDt[0]*dtdz;
                    const su2double drudz  = solDOFDr[1]*drdz + solDOFDs[1]*dsdz + solDOFDt[1]*dtdz;
                    const su2double drvdz  = solDOFDr[2]*drdz + solDOFDs[2]*dsdz + solDOFDt[2]*dtdz;
                    const su2double drwdz  = solDOFDr[3]*drdz + solDOFDs[3]*dsdz + solDOFDt[3]*dtdz;
                    const su2double drEdz  = solDOFDr[4]*drdz + solDOFDs[4]*dsdz + solDOFDt[4]*dtdz;

                    /* Compute the velocities and static energy in this
                       integration point. */
                    const su2double DensityInv   = 1.0/sol[0];
                    const su2double u            = DensityInv*sol[1];
                    const su2double v            = DensityInv*sol[2];
                    const su2double w            = DensityInv*sol[3];
                    const su2double TotalEnergy  = DensityInv*sol[4];
                    const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v + w*w);

                    /* Compute the Cartesian gradients of the velocities and
                       static energy in this integration point and also the
                       divergence of the velocity. */
                    const su2double dudx = DensityInv*(drudx - u*drhodx);
                    const su2double dudy = DensityInv*(drudy - u*drhody);
                    const su2double dudz = DensityInv*(drudz - u*drhodz);

                    const su2double dvdx = DensityInv*(drvdx - v*drhodx);
                    const su2double dvdy = DensityInv*(drvdy - v*drhody);
                    const su2double dvdz = DensityInv*(drvdz - v*drhodz);

                    const su2double dwdx = DensityInv*(drwdx - w*drhodx);
                    const su2double dwdy = DensityInv*(drwdy - w*drhody);
                    const su2double dwdz = DensityInv*(drwdz - w*drhodz);

                    const su2double dStaticEnergydx = DensityInv*(drEdx - TotalEnergy*drhodx)
                                                    - u*dudx - v*dvdx - w*dwdx;
                    const su2double dStaticEnergydy = DensityInv*(drEdy - TotalEnergy*drhody)
                                                    - u*dudy - v*dvdy - w*dwdy;
                    const su2double dStaticEnergydz = DensityInv*(drEdz - TotalEnergy*drhodz)
                                                    - u*dudz - v*dvdz - w*dwdz;
                    const su2double divVel = dudx + dvdy + dwdz;

                    /* Compute the laminar viscosity. */
                    FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
                    const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

                    /* Set the value of the second viscosity and compute the
                       divergence term in the viscous normal stresses. */
                    const su2double lambda     = -TWO3*ViscosityLam;
                    const su2double lamDivTerm =  lambda*divVel;

                    /* Compute the viscous stress tensor and the normal flux.
                       Note that there is a plus sign for the heat flux, because
                       the normal points into the geometry. */
                    const su2double tauxx = 2.0*ViscosityLam*dudx + lamDivTerm;
                    const su2double tauyy = 2.0*ViscosityLam*dvdy + lamDivTerm;
                    const su2double tauzz = 2.0*ViscosityLam*dwdz + lamDivTerm;

                    const su2double tauxy = ViscosityLam*(dudy + dvdx);
                    const su2double tauxz = ViscosityLam*(dudz + dwdx);
                    const su2double tauyz = ViscosityLam*(dvdz + dwdy);

                    const su2double qHeatNorm = ViscosityLam*factHeatFlux_Lam
                                              * (dStaticEnergydx*normals[0]
                                              +  dStaticEnergydy*normals[1]
                                              +  dStaticEnergydz*normals[2]);

                    /* Update the viscous force and moment. Note that the normal
                       points into the geometry, hence the minus sign for the stress. */
                    const su2double scaleFac =  weights[i]*normals[nDim]*factor;

                    const su2double Fx = -scaleFac*(tauxx*normals[0] + tauxy*normals[1]
                                       +            tauxz*normals[2]);
                    const su2double Fy = -scaleFac*(tauxy*normals[0] + tauyy*normals[1]
                                       +            tauyz*normals[2]);
                    const su2double Fz = -scaleFac*(tauxz*normals[0] + tauyz*normals[1]
                                       +            tauzz*normals[2]);

                    ForceViscous[0] += Fx;
                    ForceViscous[1] += Fy;
                    ForceViscous[2] += Fz;

                    const su2double dx = Coord[0] - Origin[0];
                    const su2double dy = Coord[1] - Origin[1];
                    const su2double dz = Coord[2] - Origin[2];

                    MomentViscous[0] += (Fz*dy - Fy*dz)/RefLength;
                    MomentViscous[1] += (Fx*dz - Fz*dx)/RefLength;
                    MomentViscous[2] += (Fy*dx - Fx*dy)/RefLength;

                    /* Update the heat flux and maximum heat flux for this marker. */
                    Heat_Visc[iMarker] += qHeatNorm*weights[i]*normals[nDim]*RefHeatFlux;
                    MaxHeatFlux_Visc[iMarker] = max(MaxHeatFlux_Visc[iMarker], fabs(qHeatNorm));
                  }
                }

                break;
              }
            }

            /* Update the value of the counter l to the end index of the
               current chunk. */
            l = lEnd;
          }
        }

        /*--- Project forces and store the non-dimensional coefficients ---*/
        if (nDim == 2) {
          CD_Visc[iMarker]   =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
          CL_Visc[iMarker]   = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
          CEff_Visc[iMarker] = CL_Visc[iMarker] / (CD_Visc[iMarker]+EPS);
          CMz_Visc[iMarker]  = MomentViscous[2];
          CFx_Visc[iMarker]  = ForceViscous[0];
          CFy_Visc[iMarker]  = ForceViscous[1];
        }
        if (nDim == 3) {
          CD_Visc[iMarker]   =  ForceViscous[0]*cos(Alpha)*cos(Beta)
                             +  ForceViscous[1]*sin(Beta)
                             +  ForceViscous[2]*sin(Alpha)*cos(Beta);
          CL_Visc[iMarker]   = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
          CSF_Visc[iMarker]  = -ForceViscous[0]*sin(Beta)*cos(Alpha)
                             +  ForceViscous[1]*cos(Beta)
                             -  ForceViscous[2]*sin(Beta)*sin(Alpha);
          CEff_Visc[iMarker] = CL_Visc[iMarker]/(CD_Visc[iMarker] + EPS);
          CMx_Visc[iMarker]  = MomentViscous[0];
          CMy_Visc[iMarker]  = MomentViscous[1];
          CMz_Visc[iMarker]  = MomentViscous[2];
          CFx_Visc[iMarker]  = ForceViscous[0];
          CFy_Visc[iMarker]  = ForceViscous[1];
          CFz_Visc[iMarker]  = ForceViscous[2];
        }

        AllBound_CD_Visc   += CD_Visc[iMarker];
        AllBound_CL_Visc   += CL_Visc[iMarker];
        AllBound_CSF_Visc  += CSF_Visc[iMarker];
        AllBound_CMx_Visc  += CMx_Visc[iMarker];
        AllBound_CMy_Visc  += CMy_Visc[iMarker];
        AllBound_CMz_Visc  += CMz_Visc[iMarker];
        AllBound_CFx_Visc  += CFx_Visc[iMarker];
        AllBound_CFy_Visc  += CFy_Visc[iMarker];
        AllBound_CFz_Visc  += CFz_Visc[iMarker];

        AllBound_HeatFlux_Visc   += Heat_Visc[iMarker];
        AllBound_MaxHeatFlux_Visc = max(AllBound_MaxHeatFlux_Visc,
                                        MaxHeatFlux_Visc[iMarker]);

        /*--- Compute the coefficients per surface ---*/
        for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                         ++iMarker_Monitoring) {
          string Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CL_Visc[iMarker_Monitoring]   += CL_Visc[iMarker];
            Surface_CD_Visc[iMarker_Monitoring]   += CD_Visc[iMarker];
            Surface_CSF_Visc[iMarker_Monitoring]  += CSF_Visc[iMarker];
            Surface_CEff_Visc[iMarker_Monitoring] += CEff_Visc[iMarker];
            Surface_CFx_Visc[iMarker_Monitoring]  += CFx_Visc[iMarker];
            Surface_CFy_Visc[iMarker_Monitoring]  += CFy_Visc[iMarker];
            Surface_CFz_Visc[iMarker_Monitoring]  += CFz_Visc[iMarker];
            Surface_CMx_Visc[iMarker_Monitoring]  += CMx_Visc[iMarker];
            Surface_CMy_Visc[iMarker_Monitoring]  += CMy_Visc[iMarker];
            Surface_CMz_Visc[iMarker_Monitoring]  += CMz_Visc[iMarker];
          }
        }
      }
    }
  }

#ifdef HAVE_MPI

  /*--- Parallel mode. The data from all ranks must be gathered.
        Determine the size of the communication buffer. ---*/
  const unsigned long nCommSize = 9*config->GetnMarker_Monitoring() + 10;

  /*--- Define the communication buffers and store to local data in
        the local buffer. ---*/
  vector<su2double> locBuf(nCommSize), globBuf(nCommSize);

  unsigned long ii = 0;
  locBuf[ii++] = AllBound_CD_Visc;  locBuf[ii++] = AllBound_CL_Visc;
  locBuf[ii++] = AllBound_CSF_Visc; locBuf[ii++] = AllBound_CMx_Visc;
  locBuf[ii++] = AllBound_CMy_Visc; locBuf[ii++] = AllBound_CMz_Visc;
  locBuf[ii++] = AllBound_CFx_Visc; locBuf[ii++] = AllBound_CFy_Visc;
  locBuf[ii++] = AllBound_CFz_Visc; locBuf[ii++] = AllBound_HeatFlux_Visc;

  for(unsigned short i=0; i<config->GetnMarker_Monitoring(); ++i) {
    locBuf[ii++] = Surface_CL_Visc[i];  locBuf[ii++] = Surface_CD_Visc[i];
    locBuf[ii++] = Surface_CSF_Visc[i]; locBuf[ii++] = Surface_CFx_Visc[i];
    locBuf[ii++] = Surface_CFy_Visc[i]; locBuf[ii++] = Surface_CFz_Visc[i];
    locBuf[ii++] = Surface_CMx_Visc[i]; locBuf[ii++] = Surface_CMy_Visc[i];
    locBuf[ii++] = Surface_CMz_Visc[i];
  }

  /* Sum up all the data from all ranks. The result will be available on all ranks. */
  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Allreduce(locBuf.data(), globBuf.data(), nCommSize, MPI_DOUBLE,
                       MPI_SUM, SU2_MPI::GetComm());
  }

  /*--- Copy the data back from globBuf into the required variables. ---*/
  ii = 0;
  AllBound_CD_Visc  = globBuf[ii++]; AllBound_CL_Visc       = globBuf[ii++];
  AllBound_CSF_Visc = globBuf[ii++]; AllBound_CMx_Visc      = globBuf[ii++];
  AllBound_CMy_Visc = globBuf[ii++]; AllBound_CMz_Visc      = globBuf[ii++];
  AllBound_CFx_Visc = globBuf[ii++]; AllBound_CFy_Visc      = globBuf[ii++];
  AllBound_CFz_Visc = globBuf[ii++]; AllBound_HeatFlux_Visc = globBuf[ii++];

  AllBound_CEff_Visc = AllBound_CL_Visc/(AllBound_CD_Visc + EPS);

  for(unsigned short i=0; i<config->GetnMarker_Monitoring(); ++i) {
    Surface_CL_Visc[i]  = globBuf[ii++]; Surface_CD_Visc[i]  = globBuf[ii++];
    Surface_CSF_Visc[i] = globBuf[ii++]; Surface_CFx_Visc[i] = globBuf[ii++];
    Surface_CFy_Visc[i] = globBuf[ii++]; Surface_CFz_Visc[i] = globBuf[ii++];
    Surface_CMx_Visc[i] = globBuf[ii++]; Surface_CMy_Visc[i] = globBuf[ii++];
    Surface_CMz_Visc[i] = globBuf[ii++];

    Surface_CEff_Visc[i] = Surface_CL_Visc[i]/(Surface_CD_Visc[i] + EPS);
  }

  /* Determine the maximum heat flux over all ranks. */
  su2double localMax = AllBound_MaxHeatFlux_Visc;
  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Allreduce(&localMax, &AllBound_MaxHeatFlux_Visc, 1, MPI_DOUBLE,
                       MPI_MAX, SU2_MPI::GetComm());
  }
#endif

  /*--- Update the total coefficients (note that all the nodes have the same value)---*/
  Total_CD   += AllBound_CD_Visc;
  Total_CL   += AllBound_CL_Visc;
  Total_CSF  += AllBound_CSF_Visc;
  Total_CEff  = Total_CL / (Total_CD + EPS);
  Total_CMx  += AllBound_CMx_Visc;
  Total_CMy  += AllBound_CMy_Visc;
  Total_CMz  += AllBound_CMz_Visc;
  Total_CFx  += AllBound_CFx_Visc;
  Total_CFy  += AllBound_CFy_Visc;
  Total_CFz  += AllBound_CFz_Visc;

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  for (unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                    ++iMarker_Monitoring) {
    Surface_CL[iMarker_Monitoring]   += Surface_CL_Visc[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]   += Surface_CD_Visc[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring]  += Surface_CSF_Visc[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]  = Surface_CL[iMarker_Monitoring] / (Surface_CD[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]  += Surface_CFx_Visc[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]  += Surface_CFy_Visc[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]  += Surface_CFz_Visc[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]  += Surface_CMx_Visc[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]  += Surface_CMy_Visc[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]  += Surface_CMz_Visc[iMarker_Monitoring];
  }
}

void CFEM_DG_NSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  /* Check whether or not a time stepping scheme is used. */
  const bool time_stepping = config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING;

  /* Allocate the memory for the work array and initialize it to zero to avoid
     warnings in debug mode  about uninitialized memory when padding is applied. */
  vector<su2double> workArrayVec(sizeWorkArray, 0.0);
  su2double *workArray = workArrayVec.data();

  /* Constant factor present in the heat flux vector, namely the ratio of
     thermal conductivity and viscosity. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /* Constant ratio of the second viscosity and the viscosity itself. */
  const su2double lambdaOverMu = -TWO3;

  /* The eigenvalues of the viscous Jacobian, scaled by the kinematic viscosity,
     are 1.0, 2.0 + lambdaOverMu and kOverCv/Mu. The last is variable due to the
     possible presence of an eddy viscosity, but the first two are constant and
     the maximum can be determined. */
  const su2double radOverNuTerm = max(1.0, 2.0+lambdaOverMu);

  /* Store the number of metric points per DOF, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

  /* Determine the number of elements that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nElemSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /* Initialize the minimum and maximum time step. */
  Min_Delta_Time = 1.e25; Max_Delta_Time = 0.0;

  /* Easier storage of the CFL number. Note that if we are using explicit
   time stepping, the regular CFL condition has been overwritten with the
   unsteady CFL condition in the config post-processing (if non-zero). */

  const su2double CFL = config->GetCFL(iMesh);

  /*--- Explicit time stepping with imposed time step (eventually will
   allow for local time stepping with this value imposed as the time
   for syncing the cells). If the unsteady CFL is set to zero (default),
   it uses the defined unsteady time step, otherwise it computes the time
   step based on the provided unsteady CFL. Note that the regular CFL
   option in the config is always ignored with time stepping. ---*/
  if (time_stepping && (config->GetUnst_CFL() == 0.0)) {

    /*--- Loop over the owned volume elements and set the fixed dt. ---*/
    for(unsigned long l=0; l<nVolElemOwned; ++l)
      VecDeltaTime[l] = config->GetDelta_UnstTimeND();

  } else {

    /*--- Check for a compressible solver. ---*/
    if(config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) {

      /*--- Loop over the owned volume elements. Multiple elements are treated
            simultaneously to improve the performance of the matrix
            multiplications. As a consequence, the update of the counter l
            happens at the end of this loop section. ---*/
      for(unsigned long l=0; l<nVolElemOwned;) {

        /* Determine the end index for this chunk of elements and the padded
           N value in the gemm computations. */
        unsigned long lEnd;
        unsigned short ind, llEnd, NPad;

        MetaDataChunkOfElem(volElem, l, nVolElemOwned, nElemSimul, nPadMin,
                            lEnd, ind, llEnd, NPad);

        /* Get the required data from the corresponding standard element. */
        const unsigned short nDOFs          = volElem[l].nDOFsSol;
        const su2double *matDerBasisSolDOFs = standardElementsSol[ind].GetMatDerBasisFunctionsSolDOFs();

        unsigned short nPoly = standardElementsSol[ind].GetNPoly();
        if(nPoly == 0) nPoly = 1;

        /*--- Set the pointers for the local arrays. ---*/
        su2double *solDOFs     = workArray;
        su2double *gradSolDOFs = solDOFs + nDOFs*NPad;

        /* Determine the offset between the r-derivatives and s-derivatives, which is
           also the offset between s- and t-derivatives, of the solution in the DOFs. */
        const unsigned short offDerivSol = NPad*volElem[l].nDOFsSol;

        /*--- Compute the gradients of the conserved variables if a subgrid
              scale model for LES is used. ---*/
        if( SGSModelUsed ) {

          /* Copy the solution of the DOFs into solDOFs for this chunk of elements. */
          for(unsigned short ll=0; ll<llEnd; ++ll) {

            /* Easier storage of the solution of this element. */
            const unsigned long lInd = l + ll;
            const su2double *sol = VecSolDOFs.data() + nVar*volElem[lInd].offsetDOFsSolLocal;

            /* Loop over the DOFs and copy the data. */
            const unsigned short llNVar = ll*nVar;
            for(unsigned short i=0; i<nDOFs; ++i)
              for(unsigned short mm=0; mm<nVar; ++mm)
                solDOFs[i*NPad+llNVar+mm] =  sol[i*nVar+mm];
          }

          /* Call the general function to carry out the matrix product to determine
             the gradients in the DOFs of this chunk of elements. */
          blasFunctions->gemm(nDOFs*nDim, NPad, nDOFs, matDerBasisSolDOFs, solDOFs, gradSolDOFs, config);
        }

        /*--- Make a distinction between 2D and 3D for optimal performance. ---*/
        switch( nDim ) {

          case 2: {

            /*--- 2D simulation. Loop over the chunk of elements. ---*/
            for(unsigned short ll=0; ll<llEnd; ++ll) {
              const unsigned short llNVar = ll*nVar;
              const unsigned long  lInd   = l + ll;

              /* Compute the length scale of this element and initialize the
                 inviscid and viscous spectral radii to zero. */
              const su2double lenScaleInv = nPoly/volElem[lInd].lenScale;
              const su2double lenScale    = 1.0/lenScaleInv;

              su2double charVel2Max = 0.0, radViscMax = 0.0;

              /* Loop over the DOFs of the element to determine the time step. */
              for(unsigned short i=0; i<nDOFs; ++i) {
                const su2double *solDOF  = VecSolDOFs.data()
                                         + nVar*(volElem[lInd].offsetDOFsSolLocal + i);
                const su2double *gridVel = volElem[lInd].gridVelocitiesSolDOFs.data()
                                         + i*nDim;

                /* Compute the velocities and the internal energy per unit mass. */
                const su2double DensityInv   = 1.0/solDOF[0];
                const su2double u            = DensityInv*solDOF[1];
                const su2double v            = DensityInv*solDOF[2];
                const su2double StaticEnergy = DensityInv*solDOF[3] - 0.5*(u*u + v*v);

                /*--- Compute the maximum value of the wave speed. This is a rather
                      conservative estimate. ---*/
                FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
                const su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
                const su2double SoundSpeed  = sqrt(fabs(SoundSpeed2));

                const su2double radx     = fabs(u-gridVel[0]) + SoundSpeed;
                const su2double rady     = fabs(v-gridVel[1]) + SoundSpeed;
                const su2double charVel2 = radx*radx + rady*rady;

                charVel2Max = max(charVel2Max, charVel2);

                /* Compute the laminar kinematic viscosity and check if an eddy
                   viscosity must be determined. */
                const su2double muLam = FluidModel->GetLaminarViscosity();
                su2double muTurb      = 0.0;

                if( SGSModelUsed ) {

                  /* Set the pointers to the locations where the gradients
                     of this DOF start. */
                  const su2double *solDOFDr = gradSolDOFs + i*NPad + llNVar;
                  const su2double *solDOFDs = solDOFDr    + offDerivSol;

                  /* Compute the true value of the metric terms in this DOF. Note that in
                     metricTerms the metric terms scaled by the Jacobian are stored. */
                  const su2double *metricTerms = volElem[lInd].metricTermsSolDOFs.data()
                                               + i*nMetricPerPoint;
                  const su2double JacInv       = 1.0/metricTerms[0];

                  const su2double drdx = JacInv*metricTerms[1];
                  const su2double drdy = JacInv*metricTerms[2];

                  const su2double dsdx = JacInv*metricTerms[3];
                  const su2double dsdy = JacInv*metricTerms[4];

                  /*--- Compute the Cartesian gradients of the independent solution
                        variables from the gradients in parametric coordinates and the metric
                        terms in this DOF. ---*/
                  const su2double drhodx = solDOFDr[0]*drdx + solDOFDs[0]*dsdx;
                  const su2double drudx  = solDOFDr[1]*drdx + solDOFDs[1]*dsdx;
                  const su2double drvdx  = solDOFDr[2]*drdx + solDOFDs[2]*dsdx;

                  const su2double drhody = solDOFDr[0]*drdy + solDOFDs[0]*dsdy;
                  const su2double drudy  = solDOFDr[1]*drdy + solDOFDs[1]*dsdy;
                  const su2double drvdy  = solDOFDr[2]*drdy + solDOFDs[2]*dsdy;

                  /*--- Compute the Cartesian gradients of the velocities. ---*/
                  const su2double dudx = DensityInv*(drudx - u*drhodx);
                  const su2double dvdx = DensityInv*(drvdx - v*drhodx);
                  const su2double dudy = DensityInv*(drudy - u*drhody);
                  const su2double dvdy = DensityInv*(drvdy - v*drhody);

                  /* Compute the eddy viscosity. */
                  const su2double dist = volElem[lInd].wallDistanceSolDOFs[i];
                  muTurb = SGSModel->ComputeEddyViscosity_2D(solDOF[0], dudx, dudy,
                                                             dvdx, dvdy, lenScale,
                                                             dist);
                }

                /*--- Determine the viscous spectral radius. ---*/
                const su2double mu           = muLam + muTurb;
                const su2double kOverCv      = muLam*factHeatFlux_Lam
                                             + muTurb*factHeatFlux_Turb;
                const su2double factHeatFlux = kOverCv/mu;

                const su2double radVisc = DensityInv*mu*max(radOverNuTerm, factHeatFlux);

                /* Update the maximum value of the viscous spectral radius. */
                radViscMax = max(radViscMax, radVisc);
              }

              /*--- Compute the time step for the element and update the minimum and
                    maximum value. Take the factor for time accurate local time
                    stepping into account for the minimum and maximum. ---*/
              const su2double dtInv = lenScaleInv*(sqrt(charVel2Max) + radViscMax*lenScaleInv);

              VecDeltaTime[lInd] = CFL/dtInv;

              const su2double dtEff = volElem[lInd].factTimeLevel*VecDeltaTime[lInd];
              Min_Delta_Time = min(Min_Delta_Time, dtEff);
              Max_Delta_Time = max(Max_Delta_Time, dtEff);
            }

            break;
          }

          /*------------------------------------------------------------------*/

          case 3: {

            /*--- 3D simulation. Loop over the chunk of elements. ---*/
            for(unsigned short ll=0; ll<llEnd; ++ll) {
              const unsigned short llNVar = ll*nVar;
              const unsigned long  lInd   = l + ll;

              /* Compute the length scale of this element and initialize the
                 inviscid and viscous spectral radii to zero. */
              const su2double lenScaleInv = nPoly/volElem[lInd].lenScale;
              const su2double lenScale    = 1.0/lenScaleInv;

              su2double charVel2Max = 0.0, radViscMax = 0.0;

              /* Loop over the DOFs of the element to determine the time step. */
              for(unsigned short i=0; i<nDOFs; ++i) {
                const su2double *solDOF  = VecSolDOFs.data()
                                         + nVar*(volElem[lInd].offsetDOFsSolLocal + i);
                const su2double *gridVel = volElem[lInd].gridVelocitiesSolDOFs.data()
                                         + i*nDim;

                /* Compute the velocities and the internal energy per unit mass. */
                const su2double DensityInv   = 1.0/solDOF[0];
                const su2double u            = DensityInv*solDOF[1];
                const su2double v            = DensityInv*solDOF[2];
                const su2double w            = DensityInv*solDOF[3];
                const su2double StaticEnergy = DensityInv*solDOF[4] - 0.5*(u*u + v*v + w*w);

                /*--- Compute the maximum value of the wave speed. This is a rather
                      conservative estimate. ---*/
                FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
                const su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
                const su2double SoundSpeed  = sqrt(fabs(SoundSpeed2));

                const su2double radx     = fabs(u-gridVel[0]) + SoundSpeed;
                const su2double rady     = fabs(v-gridVel[1]) + SoundSpeed;
                const su2double radz     = fabs(w-gridVel[2]) + SoundSpeed;
                const su2double charVel2 = radx*radx + rady*rady + radz*radz;

                charVel2Max = max(charVel2Max, charVel2);

                /* Compute the laminar kinematic viscosity and check if an eddy
                   viscosity must be determined. */
                const su2double muLam = FluidModel->GetLaminarViscosity();
                su2double muTurb      = 0.0;

                if( SGSModelUsed ) {

                  /* Set the pointers to the locations where the gradients
                     of this DOF start. */
                  const su2double *solDOFDr = gradSolDOFs + i*NPad + llNVar;
                  const su2double *solDOFDs = solDOFDr    + offDerivSol;
                  const su2double *solDOFDt = solDOFDs    + offDerivSol;

                  /* Compute the true value of the metric terms in this DOF. Note that in
                     metricTerms the metric terms scaled by the Jacobian are stored. */
                  const su2double *metricTerms = volElem[lInd].metricTermsSolDOFs.data()
                                               + i*nMetricPerPoint;
                  const su2double JacInv       = 1.0/metricTerms[0];

                  const su2double drdx = JacInv*metricTerms[1];
                  const su2double drdy = JacInv*metricTerms[2];
                  const su2double drdz = JacInv*metricTerms[3];

                  const su2double dsdx = JacInv*metricTerms[4];
                  const su2double dsdy = JacInv*metricTerms[5];
                  const su2double dsdz = JacInv*metricTerms[6];

                  const su2double dtdx = JacInv*metricTerms[7];
                  const su2double dtdy = JacInv*metricTerms[8];
                  const su2double dtdz = JacInv*metricTerms[9];

                  /*--- Compute the Cartesian gradients of the independent solution
                        variables from the gradients in parametric coordinates and the metric
                        terms in this DOF. ---*/
                  const su2double drhodx = solDOFDr[0]*drdx + solDOFDs[0]*dsdx + solDOFDt[0]*dtdx;
                  const su2double drux   = solDOFDr[1]*drdx + solDOFDs[1]*dsdx + solDOFDt[1]*dtdx;
                  const su2double drvx   = solDOFDr[2]*drdx + solDOFDs[2]*dsdx + solDOFDt[2]*dtdx;
                  const su2double drwx   = solDOFDr[3]*drdx + solDOFDs[3]*dsdx + solDOFDt[3]*dtdx;

                  const su2double drhody = solDOFDr[0]*drdy + solDOFDs[0]*dsdy + solDOFDt[0]*dtdy;
                  const su2double druy   = solDOFDr[1]*drdy + solDOFDs[1]*dsdy + solDOFDt[1]*dtdy;
                  const su2double drvy   = solDOFDr[2]*drdy + solDOFDs[2]*dsdy + solDOFDt[2]*dtdy;
                  const su2double drwy   = solDOFDr[3]*drdy + solDOFDs[3]*dsdy + solDOFDt[3]*dtdy;

                  const su2double drhodz = solDOFDr[0]*drdz + solDOFDs[0]*dsdz + solDOFDt[0]*dtdz;
                  const su2double druz   = solDOFDr[1]*drdz + solDOFDs[1]*dsdz + solDOFDt[1]*dtdz;
                  const su2double drvz   = solDOFDr[2]*drdz + solDOFDs[2]*dsdz + solDOFDt[2]*dtdz;
                  const su2double drwz   = solDOFDr[3]*drdz + solDOFDs[3]*dsdz + solDOFDt[3]*dtdz;

                  /*--- Compute the Cartesian gradients of the velocities. ---*/
                  const su2double dudx = DensityInv*(drux - u*drhodx);
                  const su2double dudy = DensityInv*(druy - u*drhody);
                  const su2double dudz = DensityInv*(druz - u*drhodz);

                  const su2double dvdx = DensityInv*(drvx - v*drhodx);
                  const su2double dvdy = DensityInv*(drvy - v*drhody);
                  const su2double dvdz = DensityInv*(drvz - v*drhodz);

                  const su2double dwdx = DensityInv*(drwx - w*drhodx);
                  const su2double dwdy = DensityInv*(drwy - w*drhody);
                  const su2double dwdz = DensityInv*(drwz - w*drhodz);

                  /* Compute the eddy viscosity. */
                  const su2double dist = volElem[lInd].wallDistanceSolDOFs[i];
                  muTurb = SGSModel->ComputeEddyViscosity_3D(solDOF[0], dudx, dudy, dudz,
                                                             dvdx, dvdy, dvdz, dwdx, dwdy,
                                                             dwdz, lenScale, dist);
                }

                /*--- Determine the viscous spectral radius. ---*/
                const su2double mu           = muLam + muTurb;
                const su2double kOverCv      = muLam*factHeatFlux_Lam
                                             + muTurb*factHeatFlux_Turb;
                const su2double factHeatFlux = kOverCv/mu;

                const su2double radVisc = DensityInv*mu*max(radOverNuTerm, factHeatFlux);

                /* Update the maximum value of the viscous spectral radius. */
                radViscMax = max(radViscMax, radVisc);
              }

              /*--- Compute the time step for the element and update the minimum and
                    maximum value. Take the factor for time accurate local time
                    stepping into account for the minimum and maximum. ---*/
              const su2double dtInv = lenScaleInv*(sqrt(charVel2Max) + radViscMax*lenScaleInv);

              VecDeltaTime[lInd] = CFL/dtInv;

              const su2double dtEff = volElem[lInd].factTimeLevel*VecDeltaTime[lInd];
              Min_Delta_Time = min(Min_Delta_Time, dtEff);
              Max_Delta_Time = max(Max_Delta_Time, dtEff);
            }

            break;
          }
        }

        /* Update the value of the counter l to the end index of the
           current chunk. */
        l = lEnd;
      }
    }
    else {

      /*--- Incompressible solver. ---*/

      SU2_MPI::Error("Incompressible solver not implemented yet", CURRENT_FUNCTION);
    }

    /*--- Compute the max and the min dt (in parallel). Note that we only
     do this for steady calculations if the high verbosity is set, but we
     always perform the reduction for unsteady calculations where the CFL
     limit is used to set the global time step. ---*/
    if ((config->GetComm_Level() == COMM_FULL) || time_stepping) {
#ifdef HAVE_MPI
      su2double rbuf_time = Min_Delta_Time;
      SU2_MPI::Allreduce(&rbuf_time, &Min_Delta_Time, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());

      rbuf_time = Max_Delta_Time;
      SU2_MPI::Allreduce(&rbuf_time, &Max_Delta_Time, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
#endif
    }

    /*--- For explicit time stepping with an unsteady CFL imposed, use the
          minimum delta time of the entire mesh. As Min_Delta_Time is scaled to
          the time step of the largest time level, a correction must be used
          for the time level when time accurate local time stepping is used. ---*/
    if (time_stepping) {
      for(unsigned long l=0; l<nVolElemOwned; ++l)
        VecDeltaTime[l] = Min_Delta_Time/volElem[l].factTimeLevel;

      config->SetDelta_UnstTimeND(Min_Delta_Time);
    }
  }
}

void CFEM_DG_NSSolver::ADER_DG_AliasedPredictorResidual_2D(CConfig              *config,
                                                           CVolumeElementFEM    *elem,
                                                           const su2double      *sol,
                                                           const unsigned short nSimul,
                                                           const unsigned short NPad,
                                                           su2double            *res,
                                                           su2double            *work) {
  /* Constant factor present in the heat flux vector. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /* Get the necessary information from the standard element. */
  const unsigned short ind                = elem->indStandardElement;
  const unsigned short nInt               = standardElementsSol[ind].GetNIntegration();
  const unsigned short nDOFs              = elem->nDOFsSol;
  const su2double *matBasisInt            = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
  const su2double *matDerBasisInt         = matBasisInt + nDOFs*nInt;
  const su2double *matDerBasisSolDOFs     = standardElementsSol[ind].GetMatDerBasisFunctionsSolDOFs();
  const su2double *basisFunctionsIntTrans = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
  const su2double *weights                = standardElementsSol[ind].GetWeightsIntegration();

  unsigned short nPoly = standardElementsSol[ind].GetNPoly();
  if(nPoly == 0) nPoly = 1;

  /* Compute the length scale of the current element for the LES. */
  const su2double lenScale = elem->lenScale/nPoly;

  /* Set the pointers for fluxes in the DOFs, the gradient of the fluxes in
     the integration points, the gradient of the solution in the DOFs and the
     divergence of the fluxes in the integration points. Note that some pointers
     point to the same physical location. This is because this memory can be
     used for different purposes. */
  su2double *fluxXDOF     = work;
  su2double *fluxYDOF     = fluxXDOF + NPad*nDOFs;
  su2double *gradFluxXInt = fluxYDOF + NPad*nDOFs;
  su2double *gradFluxYInt = gradFluxXInt + nDim*NPad*nInt;
  su2double *gradSolDOFs  = gradFluxXInt;
  su2double *divFlux      = work;

  /* Determine the offset between the r-derivatives and s-derivatives of the
     fluxes in the integration points and the offset between the r-derivatives
     and s-derivatives of the solution in the DOFs. */
  const unsigned short offDerivSol    = NPad*nDOFs;
  const unsigned short offDerivFluxes = NPad*nInt;

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 5;  /* nDim*nDim + 1. */

  /*--------------------------------------------------------------------------*/
  /*---          Construct the Cartesian fluxes in the DOFs.               ---*/
  /*--------------------------------------------------------------------------*/

  /* Compute the derivatives of the solution variables w.r.t. the parametric
     coordinates in the DOFs. */
  blasFunctions->gemm(nDOFs*nDim, NPad, nDOFs, matDerBasisSolDOFs, sol, gradSolDOFs, config);

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /*--- Loop over the DOFs to compute the Cartesian fluxes in the DOFs. ---*/
    for(unsigned short i=0; i<nDOFs; ++i) {

      /* Set the pointers for the solution, its gradients and the fluxes
         for this DOF. */
      const unsigned short offDOF = i*NPad + simul*nVar;
      const su2double *solDOF     = sol         + offDOF;
      const su2double *solDOFDr   = gradSolDOFs + offDOF;
      const su2double *solDOFDs   = solDOFDr    + offDerivSol;
      su2double *fluxX            = fluxXDOF    + offDOF;
      su2double *fluxY            = fluxYDOF    + offDOF;

      /* Set the pointer to the grid velocities at the location of the
         solution DOFS. THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL
         MOTION IS SPECIFIED, THE DATA FOR THIS DOF FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *gridVel = elem->gridVelocitiesSolDOFs.data() + 2*i; /* nDim*i. */

      /* Compute the true value of the metric terms in this DOF. Note that in
         metricTerms the metric terms scaled by the Jacobian are stored. THIS
         IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED, THE
         DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST BE TAKEN. */
      const su2double *metricTerms = elem->metricTermsSolDOFs.data()
                                   + i*nMetricPerPoint;
      const su2double JacInv       = 1.0/metricTerms[0];

      const su2double drdx = JacInv*metricTerms[1];
      const su2double drdy = JacInv*metricTerms[2];

      const su2double dsdx = JacInv*metricTerms[3];
      const su2double dsdy = JacInv*metricTerms[4];

      /* Compute the Cartesian gradients of the independent solution variables
         from the gradients in parametric coordinates and the metric terms. */
      const su2double drhodx = solDOFDr[0]*drdx + solDOFDs[0]*dsdx;
      const su2double drudx  = solDOFDr[1]*drdx + solDOFDs[1]*dsdx;
      const su2double drvdx  = solDOFDr[2]*drdx + solDOFDs[2]*dsdx;
      const su2double drEdx  = solDOFDr[3]*drdx + solDOFDs[3]*dsdx;

      const su2double drhody = solDOFDr[0]*drdy + solDOFDs[0]*dsdy;
      const su2double drudy  = solDOFDr[1]*drdy + solDOFDs[1]*dsdy;
      const su2double drvdy  = solDOFDr[2]*drdy + solDOFDs[2]*dsdy;
      const su2double drEdy  = solDOFDr[3]*drdy + solDOFDs[3]*dsdy;

      /* Compute the velocities, pressure and laminar viscosity in this DOF. */
      const su2double DensityInv   = 1.0/solDOF[0];
      const su2double u            = DensityInv*solDOF[1];
      const su2double v            = DensityInv*solDOF[2];
      const su2double TotalEnergy  = DensityInv*solDOF[3];
      const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v);

      FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
      const su2double Pressure     = FluidModel->GetPressure();
      const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

      /* Compute the Cartesian gradients of the velocities and static energy. */
      const su2double dudx = DensityInv*(drudx - u*drhodx);
      const su2double dvdx = DensityInv*(drvdx - v*drhodx);
      const su2double dudy = DensityInv*(drudy - u*drhody);
      const su2double dvdy = DensityInv*(drvdy - v*drhody);

      const su2double dedx = DensityInv*(drEdx - TotalEnergy*drhodx) - u*dudx - v*dvdx;
      const su2double dedy = DensityInv*(drEdy - TotalEnergy*drhody) - u*dudy - v*dvdy;

      /* Compute the eddy viscosity, if needed, and the total viscosity. */
      su2double ViscosityTurb = 0.0;
      if( SGSModelUsed )
        ViscosityTurb = SGSModel->ComputeEddyViscosity_2D(solDOF[0], dudx, dudy,
                                                          dvdx, dvdy, lenScale,
                                                          elem->wallDistanceSolDOFs[i]);
      const su2double Viscosity = ViscosityLam + ViscosityTurb;

      /* Compute the total thermal conductivity divided by Cv. */
      const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                              + ViscosityTurb*factHeatFlux_Turb;

      /* Set the value of the second viscosity and compute the divergence
         term in the viscous normal stresses. */
      const su2double lambda     = -TWO3*Viscosity;
      const su2double lamDivTerm =  lambda*(dudx + dvdy);

      /* Compute the viscous stress tensor. */
      const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
      const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
      const su2double tauxy = Viscosity*(dudy + dvdx);

      /* The Cartesian fluxes in the x-direction. */
      const su2double uRel = u - gridVel[0];
      fluxX[0] = solDOF[0]*uRel;
      fluxX[1] = solDOF[1]*uRel + Pressure - tauxx;
      fluxX[2] = solDOF[2]*uRel - tauxy;
      fluxX[3] = solDOF[3]*uRel + Pressure*u - kOverCv*dedx - u*tauxx - v*tauxy;;

      /* The Cartesian fluxes in the y-direction. */
      const su2double vRel = v - gridVel[1];
      fluxY[0] = solDOF[0]*vRel;
      fluxY[1] = solDOF[1]*vRel - tauxy;
      fluxY[2] = solDOF[2]*vRel + Pressure - tauyy;
      fluxY[3] = solDOF[3]*vRel + Pressure*v - kOverCv*dedy - u*tauxy - v*tauyy;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the derivatives of the Cartesian fluxes w.r.t. the         ---*/
  /*--- parametric coordinates in the integration points.                  ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nInt*nDim, NPad, nDOFs, matDerBasisInt, fluxXDOF, gradFluxXInt, config);
  blasFunctions->gemm(nInt*nDim, NPad, nDOFs, matDerBasisInt, fluxYDOF, gradFluxYInt, config);

  /*--------------------------------------------------------------------------*/
  /*--- Compute the divergence of the fluxes in the integration points,    ---*/
  /*--- multiplied by the integration weight.                              ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /*--- Loop over the integration points of the element. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Set the pointers for this integration point. */
      const unsigned short offInt  = i*NPad + simul*nVar;
      const su2double *gradFluxXDr = gradFluxXInt + offInt;
      const su2double *gradFluxXDs = gradFluxXDr  + offDerivFluxes;
      const su2double *gradFluxYDr = gradFluxYInt + offInt;
      const su2double *gradFluxYDs = gradFluxYDr  + offDerivFluxes;
      su2double       *divFluxInt  = divFlux      + offInt;

      /* Easier storage of the metric terms in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
         BE TAKEN. */
      const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

      /* Compute the metric terms multiplied by the integration weight. Note that the
         first term in the metric terms is the Jacobian. */
      const su2double wDrdx = weights[i]*metricTerms[1];
      const su2double wDrdy = weights[i]*metricTerms[2];

      const su2double wDsdx = weights[i]*metricTerms[3];
      const su2double wDsdy = weights[i]*metricTerms[4];

      /* Compute the divergence of the fluxes, multiplied by the
         integration weight. */
      divFluxInt[0] = gradFluxXDr[0]*wDrdx + gradFluxXDs[0]*wDsdx
                    + gradFluxYDr[0]*wDrdy + gradFluxYDs[0]*wDsdy;
      divFluxInt[1] = gradFluxXDr[1]*wDrdx + gradFluxXDs[1]*wDsdx
                    + gradFluxYDr[1]*wDrdy + gradFluxYDs[1]*wDsdy;
      divFluxInt[2] = gradFluxXDr[2]*wDrdx + gradFluxXDs[2]*wDsdx
                    + gradFluxYDr[2]*wDrdy + gradFluxYDs[2]*wDsdy;
      divFluxInt[3] = gradFluxXDr[3]*wDrdx + gradFluxXDs[3]*wDsdx
                    + gradFluxYDr[3]*wDrdy + gradFluxYDs[3]*wDsdy;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the body force to the divergence of the fluxes, if such a      ---*/
  /*--- force is present.                                                  ---*/
  /*--------------------------------------------------------------------------*/

  if( config->GetBody_Force() ) {

    /* Easier storage of the body force. */
    const su2double *body_force_vector = config->GetBody_Force_Vector();

    /* Compute the solution in the integration points of the element.
       Use gradFluxYInt to store this solution. */
    su2double *solInt = gradFluxYInt;

    blasFunctions->gemm(nInt, NPad, nDOFs, matBasisInt, sol, solInt, config);

    /*--- Loop over the number of entities that are treated simultaneously. */
    for(unsigned short simul=0; simul<nSimul; ++simul) {

      /*--- Loop over the integration points of the element. ---*/
      for(unsigned short i=0; i<nInt; ++i) {

        /* Set the pointers for this integration point. */
        const unsigned short offInt = i*NPad + simul*nVar;
        const su2double *solThisInt = solInt  + offInt;
        su2double       *divFluxInt = divFlux + offInt;

        /* Easier storage of the metric terms in this integration point.
           THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
           THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
           BE TAKEN. */
        const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

        /* Compute the velocities. */
        const su2double rhoInv = 1.0/solThisInt[0];
        const su2double u      = solThisInt[1]*rhoInv;
        const su2double v      = solThisInt[2]*rhoInv;

        /* Add the body force to the flux divergence for the momentum and energy
           equation. Note that the source terms are multiplied with minus the
           integration weight in order to be consistent with the formulation of
           the residual. Also note that for the energy source term the absolute
           velocity must be taken and not the relative. */
        const su2double weightJac = weights[i]*metricTerms[0];

        divFluxInt[1] -= weightJac*body_force_vector[0];
        divFluxInt[2] -= weightJac*body_force_vector[1];
        divFluxInt[3] -= weightJac*(u*body_force_vector[0] + v*body_force_vector[1]);
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the source terms of the manufactured solution to the divergence---*/
  /*--- of the fluxes, if a manufactured solution is used.                 ---*/
  /*--------------------------------------------------------------------------*/

  if( VerificationSolution ) {
    if( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Loop over the number of entities that are treated simultaneously. */
      for(unsigned short simul=0; simul<nSimul; ++simul) {

        /*--- Loop over the integration points of the element. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /* Set the pointers for this integration point. */
          const unsigned short offInt  = i*NPad + simul*nVar;
          su2double       *divFluxInt = divFlux + offInt;

          /* Set the pointer to the coordinates in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *coor = elem->coorIntegrationPoints.data() + i*nDim;

          /* Easier storage of the metric terms in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
          const su2double weightJac    = weights[i]*metricTerms[0];

          /* Compute the source terms of the manufactured solution.
             THIS IS A TEMPORARY IMPLEMENTATION. FOR AN ACTUAL TIME ACCURATE
             SIMULATION THE CORRECT TIME MUST BE GIVEN TO THIS FUNCTION. */
          su2double sourceMan[4];
          VerificationSolution->GetMMSSourceTerm(coor, 0.0, sourceMan);

          /* Add the source terms to the flux divergence. Note that the source
             terms are multiplied with minus the integration weight in order
             to be consistent with the formulation of the residual. */
          divFluxInt[0] -= weightJac*sourceMan[0];
          divFluxInt[1] -= weightJac*sourceMan[1];
          divFluxInt[2] -= weightJac*sourceMan[2];
          divFluxInt[3] -= weightJac*sourceMan[3];
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the residual in the DOFs, which is the matrix product of   ---*/
  /*--- basisFunctionsIntTrans and divFlux.                                ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nDOFs, NPad, nInt, basisFunctionsIntTrans, divFlux, res, config);
}

void CFEM_DG_NSSolver::ADER_DG_AliasedPredictorResidual_3D(CConfig              *config,
                                                           CVolumeElementFEM    *elem,
                                                           const su2double      *sol,
                                                           const unsigned short nSimul,
                                                           const unsigned short NPad,
                                                           su2double            *res,
                                                           su2double            *work) {
  /* Constant factor present in the heat flux vector. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /*--- Get the necessary information from the standard element. ---*/
  const unsigned short ind                = elem->indStandardElement;
  const unsigned short nInt               = standardElementsSol[ind].GetNIntegration();
  const unsigned short nDOFs              = elem->nDOFsSol;
  const su2double *matBasisInt            = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
  const su2double *matDerBasisInt         = matBasisInt + nDOFs*nInt;
  const su2double *matDerBasisSolDOFs     = standardElementsSol[ind].GetMatDerBasisFunctionsSolDOFs();
  const su2double *basisFunctionsIntTrans = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
  const su2double *weights                = standardElementsSol[ind].GetWeightsIntegration();

  unsigned short nPoly = standardElementsSol[ind].GetNPoly();
  if(nPoly == 0) nPoly = 1;

  /* Compute the length scale of the current element for the LES. */
  const su2double lenScale = elem->lenScale/nPoly;

  /* Set the pointers for fluxes in the DOFs, the gradient of the fluxes in
     the integration points, the gradient of the solution in the DOFs and the
     divergence of the fluxes in the integration points. Note that some pointers
     point to the same physical location. This is because this memory can be
     used for different purposes. */
  su2double *fluxXDOF     = work;
  su2double *fluxYDOF     = fluxXDOF + NPad*nDOFs;
  su2double *fluxZDOF     = fluxYDOF + NPad*nDOFs;
  su2double *gradFluxXInt = fluxZDOF + NPad*nDOFs;
  su2double *gradFluxYInt = gradFluxXInt + nDim*NPad*nInt;
  su2double *gradFluxZInt = gradFluxYInt + nDim*NPad*nInt;
  su2double *gradSolDOFs  = gradFluxXInt;
  su2double *divFlux      = work;

  /* Determine the offset between the r-derivatives and s-derivatives of the
     fluxes in the integration points and the offset between the r-derivatives
     and s-derivatives of the solution in the DOFs. */
  const unsigned short offDerivSol    = NPad*nDOFs;
  const unsigned short offDerivFluxes = NPad*nInt;

  /* Store the number of metric points per integration point/DOF for readability. */
  const unsigned short nMetricPerPoint = 10;  /* nDim*nDim + 1. */

  /*--------------------------------------------------------------------------*/
  /*---          Construct the Cartesian fluxes in the DOFs.               ---*/
  /*--------------------------------------------------------------------------*/

  /* Compute the derivatives of the solution variables w.r.t. the parametric
     coordinates in the DOFs. */
  blasFunctions->gemm(nDOFs*nDim, NPad, nDOFs, matDerBasisSolDOFs, sol, gradSolDOFs, config);

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /*--- Loop over the DOFs to compute the Cartesian fluxes in the DOFs. ---*/
    for(unsigned short i=0; i<nDOFs; ++i) {

      /* Set the pointers for the solution, its gradients and the fluxes
         for this DOF. */
      const unsigned short offDOF = i*NPad + simul*nVar;
      const su2double *solDOF     = sol         + offDOF;
      const su2double *solDOFDr   = gradSolDOFs + offDOF;
      const su2double *solDOFDs   = solDOFDr    + offDerivSol;
      const su2double *solDOFDt   = solDOFDs    + offDerivSol;
      su2double *fluxX            = fluxXDOF    + offDOF;
      su2double *fluxY            = fluxYDOF    + offDOF;
      su2double *fluxZ            = fluxZDOF    + offDOF;

      /* Set the pointer to the grid velocities at the location of the
         solution DOFS. THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL
         MOTION IS SPECIFIED, THE DATA FOR THIS DOF FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *gridVel = elem->gridVelocitiesSolDOFs.data() + 3*i; /* nDim*i. */

      /* Compute the true value of the metric terms in this DOF. Note that in
         metricTerms the metric terms scaled by the Jacobian are stored. THIS
         IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED, THE
         DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST BE TAKEN. */
      const su2double *metricTerms = elem->metricTermsSolDOFs.data()
                                   + i*nMetricPerPoint;
      const su2double JacInv       = 1.0/metricTerms[0];

      const su2double drdx = JacInv*metricTerms[1];
      const su2double drdy = JacInv*metricTerms[2];
      const su2double drdz = JacInv*metricTerms[3];

      const su2double dsdx = JacInv*metricTerms[4];
      const su2double dsdy = JacInv*metricTerms[5];
      const su2double dsdz = JacInv*metricTerms[6];

      const su2double dtdx = JacInv*metricTerms[7];
      const su2double dtdy = JacInv*metricTerms[8];
      const su2double dtdz = JacInv*metricTerms[9];

      /* Compute the Cartesian gradients of the independent solution variables
         from the gradients in parametric coordinates and the metric terms. */
      const su2double drhodx = solDOFDr[0]*drdx + solDOFDs[0]*dsdx + solDOFDt[0]*dtdx;
      const su2double drudx  = solDOFDr[1]*drdx + solDOFDs[1]*dsdx + solDOFDt[1]*dtdx;
      const su2double drvdx  = solDOFDr[2]*drdx + solDOFDs[2]*dsdx + solDOFDt[2]*dtdx;
      const su2double drwdx  = solDOFDr[3]*drdx + solDOFDs[3]*dsdx + solDOFDt[3]*dtdx;
      const su2double drEdx  = solDOFDr[4]*drdx + solDOFDs[4]*dsdx + solDOFDt[4]*dtdx;

      const su2double drhody = solDOFDr[0]*drdy + solDOFDs[0]*dsdy + solDOFDt[0]*dtdy;
      const su2double drudy  = solDOFDr[1]*drdy + solDOFDs[1]*dsdy + solDOFDt[1]*dtdy;
      const su2double drvdy  = solDOFDr[2]*drdy + solDOFDs[2]*dsdy + solDOFDt[2]*dtdy;
      const su2double drwdy  = solDOFDr[3]*drdy + solDOFDs[3]*dsdy + solDOFDt[3]*dtdy;
      const su2double drEdy  = solDOFDr[4]*drdy + solDOFDs[4]*dsdy + solDOFDt[4]*dtdy;

      const su2double drhodz = solDOFDr[0]*drdz + solDOFDs[0]*dsdz + solDOFDt[0]*dtdz;
      const su2double drudz  = solDOFDr[1]*drdz + solDOFDs[1]*dsdz + solDOFDt[1]*dtdz;
      const su2double drvdz  = solDOFDr[2]*drdz + solDOFDs[2]*dsdz + solDOFDt[2]*dtdz;
      const su2double drwdz  = solDOFDr[3]*drdz + solDOFDs[3]*dsdz + solDOFDt[3]*dtdz;
      const su2double drEdz  = solDOFDr[4]*drdz + solDOFDs[4]*dsdz + solDOFDt[4]*dtdz;

      /* Compute the velocities, pressure and laminar viscosity in this DOF. */
      const su2double DensityInv   = 1.0/solDOF[0];
      const su2double u            = DensityInv*solDOF[1];
      const su2double v            = DensityInv*solDOF[2];
      const su2double w            = DensityInv*solDOF[3];
      const su2double TotalEnergy  = DensityInv*solDOF[4];
      const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v + w*w);

      FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
      const su2double Pressure     = FluidModel->GetPressure();
      const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

      /* Compute the Cartesian gradients of the velocities and static energy. */
      const su2double dudx = DensityInv*(drudx - u*drhodx);
      const su2double dudy = DensityInv*(drudy - u*drhody);
      const su2double dudz = DensityInv*(drudz - u*drhodz);

      const su2double dvdx = DensityInv*(drvdx - v*drhodx);
      const su2double dvdy = DensityInv*(drvdy - v*drhody);
      const su2double dvdz = DensityInv*(drvdz - v*drhodz);

      const su2double dwdx = DensityInv*(drwdx - w*drhodx);
      const su2double dwdy = DensityInv*(drwdy - w*drhody);
      const su2double dwdz = DensityInv*(drwdz - w*drhodz);

      const su2double dedx = DensityInv*(drEdx - TotalEnergy*drhodx) - u*dudx - v*dvdx - w*dwdx;
      const su2double dedy = DensityInv*(drEdy - TotalEnergy*drhody) - u*dudy - v*dvdy - w*dwdy;
      const su2double dedz = DensityInv*(drEdz - TotalEnergy*drhodz) - u*dudz - v*dvdz - w*dwdz;

      /* Compute the eddy viscosity, if needed, and the total viscosity. */
      su2double ViscosityTurb = 0.0;
      if( SGSModelUsed )
        ViscosityTurb = SGSModel->ComputeEddyViscosity_3D(solDOF[0], dudx, dudy, dudz,
                                                          dvdx, dvdy, dvdz, dwdx,
                                                          dwdy, dwdz, lenScale,
                                                          elem->wallDistanceSolDOFs[i]);
      const su2double Viscosity = ViscosityLam + ViscosityTurb;

      /* Compute the total thermal conductivity divided by Cv. */
      const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                              + ViscosityTurb*factHeatFlux_Turb;

      /* Set the value of the second viscosity and compute the divergence
         term in the viscous normal stresses. */
      const su2double lambda     = -TWO3*Viscosity;
      const su2double lamDivTerm =  lambda*(dudx + dvdy + dwdz);

      /* Compute the viscous stress tensor. */
      const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
      const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
      const su2double tauzz = 2.0*Viscosity*dwdz + lamDivTerm;

      const su2double tauxy = Viscosity*(dudy + dvdx);
      const su2double tauxz = Viscosity*(dudz + dwdx);
      const su2double tauyz = Viscosity*(dvdz + dwdy);

      /* The Cartesian fluxes in the x-direction. */
      const su2double uRel = u - gridVel[0];
      fluxX[0] = solDOF[0]*uRel;
      fluxX[1] = solDOF[1]*uRel + Pressure - tauxx;
      fluxX[2] = solDOF[2]*uRel - tauxy;
      fluxX[3] = solDOF[3]*uRel - tauxz;
      fluxX[4] = solDOF[4]*uRel + Pressure*u - kOverCv*dedx - u*tauxx - v*tauxy - w*tauxz;

      /* The Cartesian fluxes in the y-direction. */
      const su2double vRel = v - gridVel[1];
      fluxY[0] = solDOF[0]*vRel;
      fluxY[1] = solDOF[1]*vRel - tauxy;
      fluxY[2] = solDOF[2]*vRel + Pressure - tauyy;
      fluxY[3] = solDOF[3]*vRel - tauyz;
      fluxY[4] = solDOF[4]*vRel + Pressure*v - kOverCv*dedy - u*tauxy - v*tauyy - w*tauyz;

      /* The Cartesian fluxes in the z-direction. */
      const su2double wRel = w - gridVel[2];
      fluxZ[0] = solDOF[0]*wRel;
      fluxZ[1] = solDOF[1]*wRel - tauxz;
      fluxZ[2] = solDOF[2]*wRel - tauyz;
      fluxZ[3] = solDOF[3]*wRel + Pressure - tauzz;
      fluxZ[4] = solDOF[4]*wRel + Pressure*w - kOverCv*dedz - u*tauxz - v*tauyz - w*tauzz;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the derivatives of the Cartesian fluxes w.r.t. the         ---*/
  /*--- parametric coordinates in the integration points.                  ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nInt*nDim, NPad, nDOFs, matDerBasisInt, fluxXDOF, gradFluxXInt, config);
  blasFunctions->gemm(nInt*nDim, NPad, nDOFs, matDerBasisInt, fluxYDOF, gradFluxYInt, config);
  blasFunctions->gemm(nInt*nDim, NPad, nDOFs, matDerBasisInt, fluxZDOF, gradFluxZInt, config);

  /*--------------------------------------------------------------------------*/
  /*--- Compute the divergence of the fluxes in the integration points,    ---*/
  /*--- multiplied by the integration weight.                              ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /*--- Loop over the integration points of the element. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Set the pointers for this integration point. */
      const unsigned short offInt  = i*NPad + simul*nVar;
      const su2double *gradFluxXDr = gradFluxXInt + offInt;
      const su2double *gradFluxXDs = gradFluxXDr  + offDerivFluxes;
      const su2double *gradFluxXDt = gradFluxXDs  + offDerivFluxes;
      const su2double *gradFluxYDr = gradFluxYInt + offInt;
      const su2double *gradFluxYDs = gradFluxYDr  + offDerivFluxes;
      const su2double *gradFluxYDt = gradFluxYDs  + offDerivFluxes;
      const su2double *gradFluxZDr = gradFluxZInt + offInt;
      const su2double *gradFluxZDs = gradFluxZDr  + offDerivFluxes;
      const su2double *gradFluxZDt = gradFluxZDs  + offDerivFluxes;
      su2double       *divFluxInt  = divFlux      + offInt;

      /* Easier storage of the metric terms in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

      /* Compute the metric terms multiplied by the integration weight. Note that the
         first term in the metric terms is the Jacobian. */
      const su2double wDrdx = weights[i]*metricTerms[1];
      const su2double wDrdy = weights[i]*metricTerms[2];
      const su2double wDrdz = weights[i]*metricTerms[3];

      const su2double wDsdx = weights[i]*metricTerms[4];
      const su2double wDsdy = weights[i]*metricTerms[5];
      const su2double wDsdz = weights[i]*metricTerms[6];

      const su2double wDtdx = weights[i]*metricTerms[7];
      const su2double wDtdy = weights[i]*metricTerms[8];
      const su2double wDtdz = weights[i]*metricTerms[9];

      /* Compute the divergence of the fluxes, multiplied by the integration weight. */
      divFluxInt[0] = gradFluxXDr[0]*wDrdx + gradFluxXDs[0]*wDsdx + gradFluxXDt[0]*wDtdx
                    + gradFluxYDr[0]*wDrdy + gradFluxYDs[0]*wDsdy + gradFluxYDt[0]*wDtdy
                    + gradFluxZDr[0]*wDrdz + gradFluxZDs[0]*wDsdz + gradFluxZDt[0]*wDtdz;
      divFluxInt[1] = gradFluxXDr[1]*wDrdx + gradFluxXDs[1]*wDsdx + gradFluxXDt[1]*wDtdx
                    + gradFluxYDr[1]*wDrdy + gradFluxYDs[1]*wDsdy + gradFluxYDt[1]*wDtdy
                    + gradFluxZDr[1]*wDrdz + gradFluxZDs[1]*wDsdz + gradFluxZDt[1]*wDtdz;
      divFluxInt[2] = gradFluxXDr[2]*wDrdx + gradFluxXDs[2]*wDsdx + gradFluxXDt[2]*wDtdx
                    + gradFluxYDr[2]*wDrdy + gradFluxYDs[2]*wDsdy + gradFluxYDt[2]*wDtdy
                    + gradFluxZDr[2]*wDrdz + gradFluxZDs[2]*wDsdz + gradFluxZDt[2]*wDtdz;
      divFluxInt[3] = gradFluxXDr[3]*wDrdx + gradFluxXDs[3]*wDsdx + gradFluxXDt[3]*wDtdx
                    + gradFluxYDr[3]*wDrdy + gradFluxYDs[3]*wDsdy + gradFluxYDt[3]*wDtdy
                    + gradFluxZDr[3]*wDrdz + gradFluxZDs[3]*wDsdz + gradFluxZDt[3]*wDtdz;
      divFluxInt[4] = gradFluxXDr[4]*wDrdx + gradFluxXDs[4]*wDsdx + gradFluxXDt[4]*wDtdx
                    + gradFluxYDr[4]*wDrdy + gradFluxYDs[4]*wDsdy + gradFluxYDt[4]*wDtdy
                    + gradFluxZDr[4]*wDrdz + gradFluxZDs[4]*wDsdz + gradFluxZDt[4]*wDtdz;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the body force to the divergence of the fluxes, if such a      ---*/
  /*--- force is present.                                                  ---*/
  /*--------------------------------------------------------------------------*/

  if( config->GetBody_Force() ) {

    /* Easier storage of the body force. */
    const su2double *body_force_vector = config->GetBody_Force_Vector();

    /* Compute the solution in the integration points of the element.
       Use gradFluxYInt to store this solution. */
    su2double *solInt = gradFluxYInt;

    blasFunctions->gemm(nInt, NPad, nDOFs, matBasisInt, sol, solInt, config);

    /*--- Loop over the number of entities that are treated simultaneously. */
    for(unsigned short simul=0; simul<nSimul; ++simul) {

      /*--- Loop over the integration points of the element. ---*/
      for(unsigned short i=0; i<nInt; ++i) {

        /* Set the pointers for this integration point. */
        const unsigned short offInt = i*NPad + simul*nVar;
        const su2double *solThisInt = solInt  + offInt;
        su2double       *divFluxInt = divFlux + offInt;

        /* Easier storage of the metric terms in this integration point.
           THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
           THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
           BE TAKEN. */
        const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

        /* Compute the velocities. */
        const su2double rhoInv = 1.0/solThisInt[0];
        const su2double u      = solThisInt[1]*rhoInv;
        const su2double v      = solThisInt[2]*rhoInv;
        const su2double w      = solThisInt[3]*rhoInv;

        /* Add the body force to the flux divergence for the momentum and energy
           equation. Note that the source terms are multiplied with minus the
           integration weight in order to be consistent with the formulation of
           the residual. Also note that for the energy source term the absolute
           velocity must be taken and not the relative. */
        const su2double weightJac = weights[i]*metricTerms[0];

        divFluxInt[1] -= weightJac*body_force_vector[0];
        divFluxInt[2] -= weightJac*body_force_vector[1];
        divFluxInt[3] -= weightJac*body_force_vector[2];
        divFluxInt[4] -= weightJac*(u*body_force_vector[0] + v*body_force_vector[1]
                       +            w*body_force_vector[2]);
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the source terms of the manufactured solution to the divergence---*/
  /*--- of the fluxes, if a manufactured solution is used.                 ---*/
  /*--------------------------------------------------------------------------*/

  if( VerificationSolution ) {
    if( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Loop over the number of entities that are treated simultaneously. */
      for(unsigned short simul=0; simul<nSimul; ++simul) {

        /*--- Loop over the integration points of the element. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /* Set the pointers for this integration point. */
          const unsigned short offInt  = i*NPad + simul*nVar;
          su2double       *divFluxInt = divFlux + offInt;

          /* Set the pointer to the coordinates in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *coor = elem->coorIntegrationPoints.data() + i*nDim;

          /* Easier storage of the metric terms in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
          const su2double weightJac    = weights[i]*metricTerms[0];

          /* Compute the source terms of the manufactured solution.
             THIS IS A TEMPORARY IMPLEMENTATION. FOR AN ACTUAL TIME ACCURATE
             SIMULATION THE CORRECT TIME MUST BE GIVEN TO THIS FUNCTION. */
          su2double sourceMan[5];
          VerificationSolution->GetMMSSourceTerm(coor, 0.0, sourceMan);

          /* Add the source terms to the flux divergence. Note that the source
             terms are multiplied with minus the integration weight in order
             to be consistent with the formulation of the residual. */
          divFluxInt[0] -= weightJac*sourceMan[0];
          divFluxInt[1] -= weightJac*sourceMan[1];
          divFluxInt[2] -= weightJac*sourceMan[2];
          divFluxInt[3] -= weightJac*sourceMan[3];
          divFluxInt[4] -= weightJac*sourceMan[4];
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the residual in the DOFs, which is the matrix product of   ---*/
  /*--- basisFunctionsIntTrans and divFlux.                                ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nDOFs, NPad, nInt, basisFunctionsIntTrans, divFlux, res, config);
}

void CFEM_DG_NSSolver::ADER_DG_NonAliasedPredictorResidual_2D(CConfig              *config,
                                                              CVolumeElementFEM    *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {

  /* Constant factor present in the heat flux vector, the inverse of
     the specific heat at constant volume and ratio lambdaOverMu. */
  const su2double factHeatFlux_Lam  =  Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb =  Gamma/Prandtl_Turb;
  const su2double Gas_Constant      =  config->GetGas_ConstantND();
  const su2double CvInv             =  Gamma_Minus_One/Gas_Constant;
  const su2double lambdaOverMu      = -TWO3;

  /*--- Get the necessary information from the standard element. ---*/
  const unsigned short ind                = elem->indStandardElement;
  const unsigned short nInt               = standardElementsSol[ind].GetNIntegration();
  const unsigned short nDOFs              = elem->nDOFsSol;
  const su2double *matBasisInt            = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
  const su2double *mat2ndDerBasisInt      = standardElementsSol[ind].GetMat2ndDerBasisFunctionsInt();
  const su2double *basisFunctionsIntTrans = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
  const su2double *weights                = standardElementsSol[ind].GetWeightsIntegration();

  unsigned short nPoly = standardElementsSol[ind].GetNPoly();
  if(nPoly == 0) nPoly = 1;

  /* Check if a body force is present and set it accordingly. */
  su2double bodyForceX = 0.0, bodyForceY = 0.0;
  if( config->GetBody_Force() ) {
    const su2double *body_force_vector = config->GetBody_Force_Vector();
    bodyForceX = body_force_vector[0];
    bodyForceY = body_force_vector[1];
  }

  /* Compute the length scale of the current element for the LES. */
  const su2double lenScale = elem->lenScale/nPoly;

  /* Set the pointers for solAndGradInt and divFlux to work. The same array
     can be used for both help arrays. */
  su2double *solAndGradInt = work;
  su2double *divFlux       = work;

  /* Determine the offset between the solution variables and the r-derivatives,
     which is also the offset between the r- and s-derivatives in the
     integration points. */
  const unsigned short offDerivInt = NPad*nInt;

  /* Set the pointer for the second derivatives such that they are stored
     after the first derivatives. */
  su2double *secDerSol = solAndGradInt + 3*NPad*nInt;   /*(nDim+1)*NPad*nInt. */

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 5;  /* nDim*nDim + 1. */

  /* Store the number of additional metric points per integration point, which
     are needed to compute the second derivatives. These terms take the
     non-constant metric into account. */
  const unsigned short nMetric2ndDerPerPoint = 6; /*nDim*(nDim + nDim*(nDim-1)/2). */

  /*--------------------------------------------------------------------------*/
  /*--- Interpolate the solution variables to the integration points and   ---*/
  /*--- also determine the first and second derivatives of these variables ---*/
  /*--- in the integration points. All derivatives are w.r.t. the          ---*/
  /*--- parametric coordinates.                                            ---*/
  /*--------------------------------------------------------------------------*/

  /* Compute the solution and the derivatives w.r.t. the parametric coordinates
     in the integration points. The first argument is nInt*(nDim+1). */
  blasFunctions->gemm(nInt*3, NPad, nDOFs, matBasisInt, sol, solAndGradInt, config);

  /* Compute the second derivatives w.r.t. the parametric coordinates
     in the integration points. */
  blasFunctions->gemm(nInt*3, NPad, nDOFs, mat2ndDerBasisInt, sol, secDerSol, config);

  /*--------------------------------------------------------------------------*/
  /*--- Compute the divergence of viscous fluxes, multiplied by the        ---*/
  /*--- integration weight in the integration points of the element.       ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /*--- Loop over the integration points. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the locations where the solution and the gradient
         data of this integration point start. */
      const unsigned short offInt = NPad*i + simul*nVar;
      const su2double *sol   = solAndGradInt + offInt;
      const su2double *solDr = sol   + offDerivInt;
      const su2double *solDs = solDr + offDerivInt;

      /* Easier storage of the locations where the second derivatives are
         stored for this integration point. */
      const su2double *solDrDr = secDerSol + offInt;
      const su2double *solDrDs = solDrDr   + offDerivInt;
      const su2double *solDsDs = solDrDs   + offDerivInt;

      /* Compute the velocities, pressure and total enthalpy
         in this integration point. */
      const su2double rho = sol[0];
      const su2double ru  = sol[1];
      const su2double rv  = sol[2];
      const su2double rE  = sol[3];

      const su2double rhoInv       = 1.0/rho;
      const su2double u            = rhoInv*ru;
      const su2double v            = rhoInv*rv;
      const su2double kinEnergy    = 0.5*(u*u + v*v);
      const su2double TotalEnergy  = rhoInv*rE;
      const su2double StaticEnergy = TotalEnergy - kinEnergy;

      FluidModel->SetTDState_rhoe(rho, StaticEnergy);
      const su2double Pressure = FluidModel->GetPressure();
      const su2double Htot     = rhoInv*(rE + Pressure);

      /* Compute the laminar viscosity and its derivative w.r.t. temperature. */
      const su2double ViscosityLam = FluidModel->GetLaminarViscosity();
      const su2double dViscLamdT   = FluidModel->GetdmudT_rho();

      /* Set the pointer to the grid velocities in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *gridVel = elem->gridVelocities.data() + 2*i; /* nDim*i. */

      /* Easier storage of the metric terms in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

      /* Compute the true metric terms. Note in metricTerms the actual metric
         terms multiplied by the Jacobian are stored. */
      const su2double Jac    = metricTerms[0];
      const su2double JacInv = 1.0/Jac;

      const su2double drdx = JacInv*metricTerms[1];
      const su2double drdy = JacInv*metricTerms[2];
      const su2double dsdx = JacInv*metricTerms[3];
      const su2double dsdy = JacInv*metricTerms[4];

      /* Compute the Cartesian gradients of the independent solution
         variables from the gradients in parametric coordinates and the
         metric terms in this integration point. */
      const su2double drhodx = solDr[0]*drdx + solDs[0]*dsdx;
      const su2double drudx  = solDr[1]*drdx + solDs[1]*dsdx;
      const su2double drvdx  = solDr[2]*drdx + solDs[2]*dsdx;
      const su2double drEdx  = solDr[3]*drdx + solDs[3]*dsdx;

      const su2double drhody = solDr[0]*drdy + solDs[0]*dsdy;
      const su2double drudy  = solDr[1]*drdy + solDs[1]*dsdy;
      const su2double drvdy  = solDr[2]*drdy + solDs[2]*dsdy;
      const su2double drEdy  = solDr[3]*drdy + solDs[3]*dsdy;

      /* Pointer to the necessary additional metric terms needed to compute
         the Cartesian second derivatives for this integration point.
         HIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *metricTerms2ndDer = elem->metricTerms2ndDer.data()
                                         + i*nMetric2ndDerPerPoint;

      /* Compute the Cartesian second derivatives of the independent solution
         variables from the gradients and second derivatives in parametric
         coordinates and the metric terms and its derivatives w.r.t. the
         parametric coordinates. */
      const su2double d2rhodxdx = solDrDr[0]*drdx*drdx + solDsDs[0]*dsdx*dsdx
                                + 2.0*solDrDs[0]*drdx*dsdx
                                + solDr[0]*metricTerms2ndDer[0] + solDs[0]*metricTerms2ndDer[1];
      const su2double d2rudxdx  = solDrDr[1]*drdx*drdx + solDsDs[1]*dsdx*dsdx
                                + 2.0*solDrDs[1]*drdx*dsdx
                                + solDr[1]*metricTerms2ndDer[0] + solDs[1]*metricTerms2ndDer[1];
      const su2double d2rvdxdx  = solDrDr[2]*drdx*drdx + solDsDs[2]*dsdx*dsdx
                                + 2.0*solDrDs[2]*drdx*dsdx
                                + solDr[2]*metricTerms2ndDer[0] + solDs[2]*metricTerms2ndDer[1];
      const su2double d2rEdxdx  = solDrDr[3]*drdx*drdx + solDsDs[3]*dsdx*dsdx
                                + 2.0*solDrDs[3]*drdx*dsdx
                                + solDr[3]*metricTerms2ndDer[0] + solDs[3]*metricTerms2ndDer[1];

      const su2double d2rhodydy = solDrDr[0]*drdy*drdy + solDsDs[0]*dsdy*dsdy
                                + 2.0*solDrDs[0]*drdy*dsdy
                                + solDr[0]*metricTerms2ndDer[4] + solDs[0]*metricTerms2ndDer[5];
      const su2double d2rudydy  = solDrDr[1]*drdy*drdy + solDsDs[1]*dsdy*dsdy
                                + 2.0*solDrDs[1]*drdy*dsdy
                                + solDr[1]*metricTerms2ndDer[4] + solDs[1]*metricTerms2ndDer[5];
      const su2double d2rvdydy  = solDrDr[2]*drdy*drdy + solDsDs[2]*dsdy*dsdy
                                + 2.0*solDrDs[2]*drdy*dsdy
                                + solDr[2]*metricTerms2ndDer[4] + solDs[2]*metricTerms2ndDer[5];
      const su2double d2rEdydy  = solDrDr[3]*drdy*drdy + solDsDs[3]*dsdy*dsdy
                                + 2.0*solDrDs[3]*drdy*dsdy
                                + solDr[3]*metricTerms2ndDer[4] + solDs[3]*metricTerms2ndDer[5];

      const su2double d2rhodxdy = solDrDr[0]*drdx*drdy + solDsDs[0]*dsdx*dsdy
                                + solDrDs[0]*(drdx*dsdy + dsdx*drdy)
                                + solDr[0]*metricTerms2ndDer[2] + solDs[0]*metricTerms2ndDer[3];
      const su2double d2rudxdy  = solDrDr[1]*drdx*drdy + solDsDs[1]*dsdx*dsdy
                                + solDrDs[1]*(drdx*dsdy + dsdx*drdy)
                                + solDr[1]*metricTerms2ndDer[2] + solDs[1]*metricTerms2ndDer[3];
      const su2double d2rvdxdy  = solDrDr[2]*drdx*drdy + solDsDs[2]*dsdx*dsdy
                                + solDrDs[2]*(drdx*dsdy + dsdx*drdy)
                                + solDr[2]*metricTerms2ndDer[2] + solDs[2]*metricTerms2ndDer[3];

      /* Compute the Cartesian gradients of the pressure, velocity components,
         static energy and dynamic viscosity. */
      const su2double dpdx = Gamma_Minus_One*(drEdx + kinEnergy*drhodx
                           -                  u*drudx - v*drvdx);
      const su2double dpdy = Gamma_Minus_One*(drEdy + kinEnergy*drhody
                           -                  u*drudy - v*drvdy);

      const su2double dudx = rhoInv*(drudx - u*drhodx);
      const su2double dudy = rhoInv*(drudy - u*drhody);
      const su2double dvdx = rhoInv*(drvdx - v*drhodx);
      const su2double dvdy = rhoInv*(drvdy - v*drhody);

      const su2double dedx = rhoInv*(drEdx - TotalEnergy*drhodx) - u*dudx - v*dvdx;
      const su2double dedy = rhoInv*(drEdy - TotalEnergy*drhody) - u*dudy - v*dvdy;

      const su2double dViscLamdx = CvInv*dedx*dViscLamdT;
      const su2double dViscLamdy = CvInv*dedy*dViscLamdT;

      /* Compute the second derivatives of the velocity components. */
      const su2double d2udxdx = rhoInv*(d2rudxdx - u*d2rhodxdx
                              +         2.0*rhoInv*drhodx*(u*drhodx - drudx));
      const su2double d2udydy = rhoInv*(d2rudydy - u*d2rhodydy
                              +         2.0*rhoInv*drhody*(u*drhody - drudy));
      const su2double d2udxdy = rhoInv*(d2rudxdy - u*d2rhodxdy
                              +         rhoInv*(drhodx*(u*drhody - drudy)
                              +                 drhody*(u*drhodx - drudx)));

      const su2double d2vdxdx = rhoInv*(d2rvdxdx - v*d2rhodxdx
                              +         2.0*rhoInv*drhodx*(v*drhodx - drvdx));
      const su2double d2vdydy = rhoInv*(d2rvdydy - v*d2rhodydy
                              +         2.0*rhoInv*drhody*(v*drhody - drvdy));
      const su2double d2vdxdy = rhoInv*(d2rvdxdy - v*d2rhodxdy
                              +         rhoInv*(drhodx*(v*drhody - drvdy)
                              +                 drhody*(v*drhodx - drvdx)));

      /* Compute the second derivatives of the static energy. Note that this
         term appears in the heat flux and therefore only the pure second
         derivatives are needed. Hence, the cross-derivatives are omitted. */
      const su2double d2edxdx = rhoInv*(d2rEdxdx - TotalEnergy*d2rhodxdx
                              +         2.0*rhoInv*drhodx*(TotalEnergy*drhodx - drEdx))
                              -         u*d2udxdx - dudx*dudx - v*d2vdxdx - dvdx*dvdx;
      const su2double d2edydy = rhoInv*(d2rEdydy - TotalEnergy*d2rhodydy
                              +         2.0*rhoInv*drhody*(TotalEnergy*drhody - drEdy))
                              -         u*d2udydy - dudy*dudy - v*d2vdydy - dvdy*dvdy;

      /* If an SGS model is used the eddy viscosity and its spatial
         derivatives must be computed. */
      su2double ViscosityTurb = 0.0;
      su2double dViscTurbdx = 0.0, dViscTurbdy = 0.0;

      if( SGSModelUsed ) {
        const su2double dist = elem->wallDistance[i];
        ViscosityTurb = SGSModel->ComputeEddyViscosity_2D(rho, dudx, dudy, dvdx,
                                                          dvdy, lenScale, dist);

        SGSModel->ComputeGradEddyViscosity_2D(rho, drhodx, drhody, dudx, dudy,
                                              dvdx, dvdy, d2udxdx, d2udydy, d2udxdy,
                                              d2vdxdx, d2vdydy, d2vdxdy, lenScale,
                                              dist, dViscTurbdx, dViscTurbdy);
      }

      /* Compute the total viscosity, the total heat conductivity and their
         gradients. Note that the heat conductivity is divided by the Cv,
         because gradients of internal energy are computed and not temperature. */
      const su2double Viscosity = ViscosityLam + ViscosityTurb;
      const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                              + ViscosityTurb*factHeatFlux_Turb;

      const su2double dViscDx = dViscLamdx + dViscTurbdx;
      const su2double dViscDy = dViscLamdy + dViscTurbdy;

      const su2double dkOverCvdx = dViscLamdx *factHeatFlux_Lam
                                 + dViscTurbdx*factHeatFlux_Turb;
      const su2double dkOverCvdy = dViscLamdy *factHeatFlux_Lam
                                 + dViscTurbdy*factHeatFlux_Turb;

      /* Abbreviations, which make it easier to compute the divergence term. */
      const su2double abv1 = drudx + drvdy;
      const su2double abv2 = u*drhodx + v*drhody;
      const su2double abv3 = u*(drEdx + dpdx) + v*(drEdy + dpdy);
      const su2double abv4 = dudx + dvdy;

      /* Compute the divergence of the grid velocity.
         SET TO ZERO FOR NOW. THIS IS NOT CORRECT!!!!. */
      const su2double divGridVel = 0.0;

      /* Set the pointer to store the divergence terms for this integration
         point and compute these terms, multiplied by the integration weight
         and Jacobian. */
      const su2double weightJac = weights[i]*Jac;
      su2double *divFluxInt     = divFlux + offInt;

      divFluxInt[0] = weightJac*(abv1 - rho*divGridVel
                    -            gridVel[0]*drhodx - gridVel[1]*drhody);
      divFluxInt[1] = weightJac*(dpdx + u*(abv1-abv2) - lambdaOverMu*abv4*dViscDx
                    +            u*drudx + v*drudy
                    -            lambdaOverMu*Viscosity*(d2udxdx + d2vdxdy)
                    -            Viscosity*(2.0*d2udxdx + d2udydy + d2vdxdy)
                    -            2.0*dViscDx*dudx - dViscDy*(dudy+dvdx)
                    -            ru*divGridVel
                    -            gridVel[0]*drudx - gridVel[1]*drudy);
      divFluxInt[2] = weightJac*(dpdy + v*(abv1-abv2) - lambdaOverMu*abv4*dViscDy
                    +            u*drvdx + v*drvdy
                    -            lambdaOverMu*Viscosity*(d2udxdy + d2vdydy)
                    -            Viscosity*(2.0*d2vdydy + d2vdxdx + d2udxdy)
                    -            dViscDx*(dudy + dvdx) - 2.0*dViscDy*dvdy
                    -            rv*divGridVel
                    -            gridVel[0]*drvdx - gridVel[1]*drvdy);
      divFluxInt[3] = weightJac*(abv3 + Htot*(abv1 - abv2)
                    -            abv4*lambdaOverMu*(Viscosity*abv4 + u*dViscDx + v*dViscDy)
                    -            dkOverCvdx*dedx - dkOverCvdy*dedy - kOverCv*(d2edxdx + d2edydy)
                    -            (Viscosity*dudx + u*dViscDx)*2.0*dudx
                    -            (Viscosity*dvdy + v*dViscDy)*2.0*dvdy
                    -            (Viscosity*dudy + u*dViscDy + Viscosity*dvdx + v*dViscDx)*(dudy + dvdx)
                    -            Viscosity*u*(d2udxdx+d2udydy + (1.0+lambdaOverMu)*(d2udxdx+d2vdxdy))
                    -            Viscosity*v*(d2vdxdx+d2vdydy + (1.0+lambdaOverMu)*(d2udxdy+d2vdydy))
                    -            rE*divGridVel
                    -            gridVel[0]*drEdx - gridVel[1]*drEdy);

      /* Add the body force to the flux divergence for the momentum and energy
         equation. Note that the source terms are multiplied with minus the
         integration weight in order to be consistent with the formulation of
         the residual. Also note that for the energy source term the absolute
         velocity must be taken and not the relative. */
      divFluxInt[1] -= weightJac*bodyForceX;
      divFluxInt[2] -= weightJac*bodyForceY;
      divFluxInt[3] -= weightJac*(u*bodyForceX + v*bodyForceY);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the source terms of the manufactured solution to the divergence---*/
  /*--- of the fluxes, if a manufactured solution is used.                 ---*/
  /*--------------------------------------------------------------------------*/

  if( VerificationSolution ) {
    if( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Loop over the number of entities that are treated simultaneously. */
      for(unsigned short simul=0; simul<nSimul; ++simul) {

        /*--- Loop over the integration points of the element. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /* Set the pointers for this integration point. */
          const unsigned short offInt  = i*NPad + simul*nVar;
          su2double       *divFluxInt = divFlux + offInt;

          /* Set the pointer to the coordinates in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *coor = elem->coorIntegrationPoints.data() + i*nDim;

          /* Easier storage of the metric terms in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
          const su2double weightJac    = weights[i]*metricTerms[0];

          /* Compute the source terms of the manufactured solution.
             THIS IS A TEMPORARY IMPLEMENTATION. FOR AN ACTUAL TIME ACCURATE
             SIMULATION THE CORRECT TIME MUST BE GIVEN TO THIS FUNCTION. */
          su2double sourceMan[4];
          VerificationSolution->GetMMSSourceTerm(coor, 0.0, sourceMan);

          /* Add the source terms to the flux divergence. Note that the source
             terms are multiplied with minus the integration weight in order
             to be consistent with the formulation of the residual. */
          divFluxInt[0] -= weightJac*sourceMan[0];
          divFluxInt[1] -= weightJac*sourceMan[1];
          divFluxInt[2] -= weightJac*sourceMan[2];
          divFluxInt[3] -= weightJac*sourceMan[3];
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the residual in the DOFs, which is the matrix product of   ---*/
  /*--- basisFunctionsIntTrans and divFlux.                                ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nDOFs, NPad, nInt, basisFunctionsIntTrans, divFlux, res, config);
}

void CFEM_DG_NSSolver::ADER_DG_NonAliasedPredictorResidual_3D(CConfig              *config,
                                                              CVolumeElementFEM    *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {

  /* Constant factor present in the heat flux vector, the inverse of
     the specific heat at constant volume and ratio lambdaOverMu. */
  const su2double factHeatFlux_Lam  =  Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb =  Gamma/Prandtl_Turb;
  const su2double Gas_Constant      =  config->GetGas_ConstantND();
  const su2double CvInv             =  Gamma_Minus_One/Gas_Constant;
  const su2double lambdaOverMu      = -TWO3;

  /*--- Get the necessary information from the standard element. ---*/
  const unsigned short ind                = elem->indStandardElement;
  const unsigned short nInt               = standardElementsSol[ind].GetNIntegration();
  const unsigned short nDOFs              = elem->nDOFsSol;
  const su2double *matBasisInt            = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
  const su2double *mat2ndDerBasisInt      = standardElementsSol[ind].GetMat2ndDerBasisFunctionsInt();
  const su2double *basisFunctionsIntTrans = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
  const su2double *weights                = standardElementsSol[ind].GetWeightsIntegration();

  unsigned short nPoly = standardElementsSol[ind].GetNPoly();
  if(nPoly == 0) nPoly = 1;

  /* Check if a body force is present and set it accordingly. */
  su2double bodyForceX = 0.0, bodyForceY = 0.0, bodyForceZ = 0.0;
  if( config->GetBody_Force() ) {
    const su2double *body_force_vector = config->GetBody_Force_Vector();
    bodyForceX = body_force_vector[0];
    bodyForceY = body_force_vector[1];
    bodyForceZ = body_force_vector[2];
  }

  /* Compute the length scale of the current element for the LES. */
  const su2double lenScale = elem->lenScale/nPoly;

  /* Set the pointers for solAndGradInt and divFlux to work. The same array
     can be used for both help arrays. */
  su2double *solAndGradInt = work;
  su2double *divFlux       = work;

  /* Determine the offset between the solution variables and the r-derivatives,
     which is also the offset between the r- and s-derivatives in the
     integration points. */
  const unsigned short offDerivInt = NPad*nInt;

  /* Set the pointer for the second derivatives such that they are stored
     after the first derivatives. */
  su2double *secDerSol = solAndGradInt + 4*NPad*nInt;  /*(nDim+1)*NPad*nInt. */

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 10;  /* nDim*nDim + 1. */

  /* Store the number of additional metric points per integration point, which
     are needed to compute the second derivatives. These terms take the
     non-constant metric into account. */
  const unsigned short nMetric2ndDerPerPoint = 18; /*nDim*(nDim + nDim*(nDim-1)/2). */

  /*--------------------------------------------------------------------------*/
  /*--- Interpolate the solution variables to the integration points and   ---*/
  /*--- also determine the first and second derivatives of these variables ---*/
  /*--- in the integration points. All derivatives are w.r.t. the          ---*/
  /*--- parametric coordinates.                                            ---*/
  /*--------------------------------------------------------------------------*/

  /* Compute the solution and the derivatives w.r.t. the parametric coordinates
     in the integration points. The first argument is nInt*(nDim+1). */
  blasFunctions->gemm(nInt*4, NPad, nDOFs, matBasisInt, sol, solAndGradInt, config);

  /* Compute the second derivatives w.r.t. the parametric coordinates
     in the integration points. */
  blasFunctions->gemm(nInt*6, NPad, nDOFs, mat2ndDerBasisInt, sol, secDerSol, config);

  /*--------------------------------------------------------------------------*/
  /*--- Compute the divergence of viscous fluxes, multiplied by the        ---*/
  /*--- integration weight in the integration points of the element.       ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /*--- Loop over the integration points. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the locations where the solution and the gradient
         data of this integration point start. */
      const unsigned short offInt = NPad*i + simul*nVar;
      const su2double *sol   = solAndGradInt + offInt;
      const su2double *solDr = sol   + offDerivInt;
      const su2double *solDs = solDr + offDerivInt;
      const su2double *solDt = solDs + offDerivInt;

      /* Easier storage of the locations where the second derivatives are
         stored for this integration point. */
      const su2double *solDrDr = secDerSol + offInt;
      const su2double *solDrDs = solDrDr   + offDerivInt;
      const su2double *solDsDs = solDrDs   + offDerivInt;
      const su2double *solDrDt = solDsDs   + offDerivInt;
      const su2double *solDsDt = solDrDt   + offDerivInt;
      const su2double *solDtDt = solDsDt   + offDerivInt;

      /* Compute the velocities, pressure and total enthalpy
         in this integration point. */
      const su2double rho = sol[0];
      const su2double ru  = sol[1];
      const su2double rv  = sol[2];
      const su2double rw  = sol[3];
      const su2double rE  = sol[4];

      const su2double rhoInv       = 1.0/rho;
      const su2double u            = rhoInv*ru;
      const su2double v            = rhoInv*rv;
      const su2double w            = rhoInv*rw;
      const su2double kinEnergy    = 0.5*(u*u + v*v + w*w);
      const su2double TotalEnergy  = rhoInv*rE;
      const su2double StaticEnergy = TotalEnergy - kinEnergy;

      FluidModel->SetTDState_rhoe(rho, StaticEnergy);
      const su2double Pressure = FluidModel->GetPressure();
      const su2double Htot     = rhoInv*(rE + Pressure);

       /* Compute the laminar viscosity and its derivative w.r.t. temperature. */
      const su2double ViscosityLam = FluidModel->GetLaminarViscosity();
      const su2double dViscLamdT   = FluidModel->GetdmudT_rho();

      /* Set the pointer to the grid velocities in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *gridVel = elem->gridVelocities.data() + 3*i; /* nDim*i. */

      /* Easier storage of the metric terms in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

      /* Compute the true metric terms. Note in metricTerms the actual metric
         terms multiplied by the Jacobian are stored. */
      const su2double Jac    = metricTerms[0];
      const su2double JacInv = 1.0/Jac;

      const su2double drdx = JacInv*metricTerms[1];
      const su2double drdy = JacInv*metricTerms[2];
      const su2double drdz = JacInv*metricTerms[3];

      const su2double dsdx = JacInv*metricTerms[4];
      const su2double dsdy = JacInv*metricTerms[5];
      const su2double dsdz = JacInv*metricTerms[6];

      const su2double dtdx = JacInv*metricTerms[7];
      const su2double dtdy = JacInv*metricTerms[8];
      const su2double dtdz = JacInv*metricTerms[9];

      /* Compute the Cartesian gradients of the independent solution
         variables from the gradients in parametric coordinates and the
         metric terms in this integration point. */
      const su2double drhodx = solDr[0]*drdx + solDs[0]*dsdx + solDt[0]*dtdx;
      const su2double drudx  = solDr[1]*drdx + solDs[1]*dsdx + solDt[1]*dtdx;
      const su2double drvdx  = solDr[2]*drdx + solDs[2]*dsdx + solDt[2]*dtdx;
      const su2double drwdx  = solDr[3]*drdx + solDs[3]*dsdx + solDt[3]*dtdx;
      const su2double drEdx  = solDr[4]*drdx + solDs[4]*dsdx + solDt[4]*dtdx;

      const su2double drhody = solDr[0]*drdy + solDs[0]*dsdy + solDt[0]*dtdy;
      const su2double drudy  = solDr[1]*drdy + solDs[1]*dsdy + solDt[1]*dtdy;
      const su2double drvdy  = solDr[2]*drdy + solDs[2]*dsdy + solDt[2]*dtdy;
      const su2double drwdy  = solDr[3]*drdy + solDs[3]*dsdy + solDt[3]*dtdy;
      const su2double drEdy  = solDr[4]*drdy + solDs[4]*dsdy + solDt[4]*dtdy;

      const su2double drhodz = solDr[0]*drdz + solDs[0]*dsdz + solDt[0]*dtdz;
      const su2double drudz  = solDr[1]*drdz + solDs[1]*dsdz + solDt[1]*dtdz;
      const su2double drvdz  = solDr[2]*drdz + solDs[2]*dsdz + solDt[2]*dtdz;
      const su2double drwdz  = solDr[3]*drdz + solDs[3]*dsdz + solDt[3]*dtdz;
      const su2double drEdz  = solDr[4]*drdz + solDs[4]*dsdz + solDt[4]*dtdz;

      /* Pointer to the necessary additional metric terms needed to compute
         the Cartesian second derivatives for this integration point.
         HIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *metricTerms2ndDer = elem->metricTerms2ndDer.data()
                                         + i*nMetric2ndDerPerPoint;

      /* Compute the Cartesian second derivatives of the independent solution
         variables from the gradients and second derivatives in parametric
         coordinates and the metric terms and its derivatives w.r.t. the
         parametric coordinates. */
      const su2double d2rhodxdx = solDrDr[0]*drdx*drdx + solDsDs[0]*dsdx*dsdx + solDtDt[0]*dtdx*dtdx
                                + 2.0*(solDrDs[0]*drdx*dsdx + solDrDt[0]*drdx*dtdx + solDsDt[0]*dsdx*dtdx)
                                + solDr[0]*metricTerms2ndDer[0] + solDs[0]*metricTerms2ndDer[1]
                                + solDt[0]*metricTerms2ndDer[2];
      const su2double d2rudxdx  = solDrDr[1]*drdx*drdx + solDsDs[1]*dsdx*dsdx + solDtDt[1]*dtdx*dtdx
                                + 2.0*(solDrDs[1]*drdx*dsdx + solDrDt[1]*drdx*dtdx + solDsDt[1]*dsdx*dtdx)
                                + solDr[1]*metricTerms2ndDer[0] + solDs[1]*metricTerms2ndDer[1]
                                + solDt[1]*metricTerms2ndDer[2];
      const su2double d2rvdxdx  = solDrDr[2]*drdx*drdx + solDsDs[2]*dsdx*dsdx + solDtDt[2]*dtdx*dtdx
                                + 2.0*(solDrDs[2]*drdx*dsdx + solDrDt[2]*drdx*dtdx + solDsDt[2]*dsdx*dtdx)
                                + solDr[2]*metricTerms2ndDer[0] + solDs[2]*metricTerms2ndDer[1]
                                + solDt[2]*metricTerms2ndDer[2];
      const su2double d2rwdxdx  = solDrDr[3]*drdx*drdx + solDsDs[3]*dsdx*dsdx + solDtDt[3]*dtdx*dtdx
                                + 2.0*(solDrDs[3]*drdx*dsdx + solDrDt[3]*drdx*dtdx + solDsDt[3]*dsdx*dtdx)
                                + solDr[3]*metricTerms2ndDer[0] + solDs[3]*metricTerms2ndDer[1]
                                + solDt[3]*metricTerms2ndDer[2];
      const su2double d2rEdxdx  = solDrDr[4]*drdx*drdx + solDsDs[4]*dsdx*dsdx + solDtDt[4]*dtdx*dtdx
                                + 2.0*(solDrDs[4]*drdx*dsdx + solDrDt[4]*drdx*dtdx + solDsDt[4]*dsdx*dtdx)
                                + solDr[4]*metricTerms2ndDer[0] + solDs[4]*metricTerms2ndDer[1]
                                + solDt[4]*metricTerms2ndDer[2];

      const su2double d2rhodydy = solDrDr[0]*drdy*drdy + solDsDs[0]*dsdy*dsdy + solDtDt[0]*dtdy*dtdy
                                + 2.0*(solDrDs[0]*drdy*dsdy + solDrDt[0]*drdy*dtdy + solDsDt[0]*dsdy*dtdy)
                                + solDr[0]*metricTerms2ndDer[6] + solDs[0]*metricTerms2ndDer[7]
                                + solDt[0]*metricTerms2ndDer[8];
      const su2double d2rudydy  = solDrDr[1]*drdy*drdy + solDsDs[1]*dsdy*dsdy + solDtDt[1]*dtdy*dtdy
                                + 2.0*(solDrDs[1]*drdy*dsdy + solDrDt[1]*drdy*dtdy + solDsDt[1]*dsdy*dtdy)
                                + solDr[1]*metricTerms2ndDer[6] + solDs[1]*metricTerms2ndDer[7]
                                + solDt[1]*metricTerms2ndDer[8];
      const su2double d2rvdydy  = solDrDr[2]*drdy*drdy + solDsDs[2]*dsdy*dsdy + solDtDt[2]*dtdy*dtdy
                                + 2.0*(solDrDs[2]*drdy*dsdy + solDrDt[2]*drdy*dtdy + solDsDt[2]*dsdy*dtdy)
                                + solDr[2]*metricTerms2ndDer[6] + solDs[2]*metricTerms2ndDer[7]
                                + solDt[2]*metricTerms2ndDer[8];
      const su2double d2rwdydy  = solDrDr[3]*drdy*drdy + solDsDs[3]*dsdy*dsdy + solDtDt[3]*dtdy*dtdy
                                + 2.0*(solDrDs[3]*drdy*dsdy + solDrDt[3]*drdy*dtdy + solDsDt[3]*dsdy*dtdy)
                                + solDr[3]*metricTerms2ndDer[6] + solDs[3]*metricTerms2ndDer[7]
                                + solDt[3]*metricTerms2ndDer[8];
      const su2double d2rEdydy  = solDrDr[4]*drdy*drdy + solDsDs[4]*dsdy*dsdy + solDtDt[4]*dtdy*dtdy
                                + 2.0*(solDrDs[4]*drdy*dsdy + solDrDt[4]*drdy*dtdy + solDsDt[4]*dsdy*dtdy)
                                + solDr[4]*metricTerms2ndDer[6] + solDs[4]*metricTerms2ndDer[7]
                                + solDt[4]*metricTerms2ndDer[8];

      const su2double d2rhodzdz = solDrDr[0]*drdz*drdz + solDsDs[0]*dsdz*dsdz + solDtDt[0]*dtdz*dtdz
                                + 2.0*(solDrDs[0]*drdz*dsdz + solDrDt[0]*drdz*dtdz + solDsDt[0]*dsdz*dtdz)
                                + solDr[0]*metricTerms2ndDer[15] + solDs[0]*metricTerms2ndDer[16]
                                + solDt[0]*metricTerms2ndDer[17];
      const su2double d2rudzdz  = solDrDr[1]*drdz*drdz + solDsDs[1]*dsdz*dsdz + solDtDt[1]*dtdz*dtdz
                                + 2.0*(solDrDs[1]*drdz*dsdz + solDrDt[1]*drdz*dtdz + solDsDt[1]*dsdz*dtdz)
                                + solDr[1]*metricTerms2ndDer[15] + solDs[1]*metricTerms2ndDer[16]
                                + solDt[1]*metricTerms2ndDer[17];
      const su2double d2rvdzdz  = solDrDr[2]*drdz*drdz + solDsDs[2]*dsdz*dsdz + solDtDt[2]*dtdz*dtdz
                                + 2.0*(solDrDs[2]*drdz*dsdz + solDrDt[2]*drdz*dtdz + solDsDt[2]*dsdz*dtdz)
                                + solDr[2]*metricTerms2ndDer[15] + solDs[2]*metricTerms2ndDer[16]
                                + solDt[2]*metricTerms2ndDer[17];
      const su2double d2rwdzdz  = solDrDr[3]*drdz*drdz + solDsDs[3]*dsdz*dsdz + solDtDt[3]*dtdz*dtdz
                                + 2.0*(solDrDs[3]*drdz*dsdz + solDrDt[3]*drdz*dtdz + solDsDt[3]*dsdz*dtdz)
                                + solDr[3]*metricTerms2ndDer[15] + solDs[3]*metricTerms2ndDer[16]
                                + solDt[3]*metricTerms2ndDer[17];
      const su2double d2rEdzdz  = solDrDr[4]*drdz*drdz + solDsDs[4]*dsdz*dsdz + solDtDt[4]*dtdz*dtdz
                                + 2.0*(solDrDs[4]*drdz*dsdz + solDrDt[4]*drdz*dtdz + solDsDt[4]*dsdz*dtdz)
                                + solDr[4]*metricTerms2ndDer[15] + solDs[4]*metricTerms2ndDer[16]
                                + solDt[4]*metricTerms2ndDer[17];

      const su2double d2rhodxdy = solDrDr[0]*drdx*drdy + solDsDs[0]*dsdx*dsdy + solDtDt[0]*dtdx*dtdy
                                + solDrDs[0]*(drdx*dsdy + dsdx*drdy) + solDrDt[0]*(drdx*dtdy + dtdx*drdy)
                                + solDsDt[0]*(dsdx*dtdy + dtdx*dsdy) + solDr[0]*metricTerms2ndDer[3]
                                + solDs[0]*metricTerms2ndDer[4] + solDt[0]*metricTerms2ndDer[5];
      const su2double d2rudxdy  = solDrDr[1]*drdx*drdy + solDsDs[1]*dsdx*dsdy + solDtDt[1]*dtdx*dtdy
                                + solDrDs[1]*(drdx*dsdy + dsdx*drdy) + solDrDt[1]*(drdx*dtdy + dtdx*drdy)
                                + solDsDt[1]*(dsdx*dtdy + dtdx*dsdy) + solDr[1]*metricTerms2ndDer[3]
                                + solDs[1]*metricTerms2ndDer[4] + solDt[1]*metricTerms2ndDer[5];
      const su2double d2rvdxdy  = solDrDr[2]*drdx*drdy + solDsDs[2]*dsdx*dsdy + solDtDt[2]*dtdx*dtdy
                                + solDrDs[2]*(drdx*dsdy + dsdx*drdy) + solDrDt[2]*(drdx*dtdy + dtdx*drdy)
                                + solDsDt[2]*(dsdx*dtdy + dtdx*dsdy) + solDr[2]*metricTerms2ndDer[3]
                                + solDs[2]*metricTerms2ndDer[4] + solDt[2]*metricTerms2ndDer[5];
      const su2double d2rwdxdy  = solDrDr[3]*drdx*drdy + solDsDs[3]*dsdx*dsdy + solDtDt[3]*dtdx*dtdy
                                + solDrDs[3]*(drdx*dsdy + dsdx*drdy) + solDrDt[3]*(drdx*dtdy + dtdx*drdy)
                                + solDsDt[3]*(dsdx*dtdy + dtdx*dsdy) + solDr[3]*metricTerms2ndDer[3]
                                + solDs[3]*metricTerms2ndDer[4] + solDt[3]*metricTerms2ndDer[5];

      const su2double d2rhodxdz = solDrDr[0]*drdx*drdz + solDsDs[0]*dsdx*dsdz + solDtDt[0]*dtdx*dtdz
                                + solDrDs[0]*(drdx*dsdz + dsdx*drdz) + solDrDt[0]*(drdx*dtdz + dtdx*drdz)
                                + solDsDt[0]*(dsdx*dtdz + dtdx*dsdz) + solDr[0]*metricTerms2ndDer[9]
                                + solDs[0]*metricTerms2ndDer[10] + solDt[0]*metricTerms2ndDer[11];
      const su2double d2rudxdz  = solDrDr[1]*drdx*drdz + solDsDs[1]*dsdx*dsdz + solDtDt[1]*dtdx*dtdz
                                + solDrDs[1]*(drdx*dsdz + dsdx*drdz) + solDrDt[1]*(drdx*dtdz + dtdx*drdz)
                                + solDsDt[1]*(dsdx*dtdz + dtdx*dsdz) + solDr[1]*metricTerms2ndDer[9]
                                + solDs[1]*metricTerms2ndDer[10] + solDt[1]*metricTerms2ndDer[11];
      const su2double d2rvdxdz  = solDrDr[2]*drdx*drdz + solDsDs[2]*dsdx*dsdz + solDtDt[2]*dtdx*dtdz
                                + solDrDs[2]*(drdx*dsdz + dsdx*drdz) + solDrDt[2]*(drdx*dtdz + dtdx*drdz)
                                + solDsDt[2]*(dsdx*dtdz + dtdx*dsdz) + solDr[2]*metricTerms2ndDer[9]
                                + solDs[2]*metricTerms2ndDer[10] + solDt[2]*metricTerms2ndDer[11];
      const su2double d2rwdxdz  = solDrDr[3]*drdx*drdz + solDsDs[3]*dsdx*dsdz + solDtDt[3]*dtdx*dtdz
                                + solDrDs[3]*(drdx*dsdz + dsdx*drdz) + solDrDt[3]*(drdx*dtdz + dtdx*drdz)
                                + solDsDt[3]*(dsdx*dtdz + dtdx*dsdz) + solDr[3]*metricTerms2ndDer[9]
                                + solDs[3]*metricTerms2ndDer[10] + solDt[3]*metricTerms2ndDer[11];

      const su2double d2rhodydz = solDrDr[0]*drdy*drdz + solDsDs[0]*dsdy*dsdz + solDtDt[0]*dtdy*dtdz
                                + solDrDs[0]*(drdy*dsdz + dsdy*drdz) + solDrDt[0]*(drdy*dtdz + dtdy*drdz)
                                + solDsDt[0]*(dsdy*dtdz + dtdy*dsdz) + solDr[0]*metricTerms2ndDer[12]
                                + solDs[0]*metricTerms2ndDer[13] + solDt[0]*metricTerms2ndDer[14];
      const su2double d2rudydz  = solDrDr[1]*drdy*drdz + solDsDs[1]*dsdy*dsdz + solDtDt[1]*dtdy*dtdz
                                + solDrDs[1]*(drdy*dsdz + dsdy*drdz) + solDrDt[1]*(drdy*dtdz + dtdy*drdz)
                                + solDsDt[1]*(dsdy*dtdz + dtdy*dsdz) + solDr[1]*metricTerms2ndDer[12]
                                + solDs[1]*metricTerms2ndDer[13] + solDt[1]*metricTerms2ndDer[14];
      const su2double d2rvdydz  = solDrDr[2]*drdy*drdz + solDsDs[2]*dsdy*dsdz + solDtDt[2]*dtdy*dtdz
                                + solDrDs[2]*(drdy*dsdz + dsdy*drdz) + solDrDt[2]*(drdy*dtdz + dtdy*drdz)
                                + solDsDt[2]*(dsdy*dtdz + dtdy*dsdz) + solDr[2]*metricTerms2ndDer[12]
                                + solDs[2]*metricTerms2ndDer[13] + solDt[2]*metricTerms2ndDer[14];
      const su2double d2rwdydz  = solDrDr[3]*drdy*drdz + solDsDs[3]*dsdy*dsdz + solDtDt[3]*dtdy*dtdz
                                + solDrDs[3]*(drdy*dsdz + dsdy*drdz) + solDrDt[3]*(drdy*dtdz + dtdy*drdz)
                                + solDsDt[3]*(dsdy*dtdz + dtdy*dsdz) + solDr[3]*metricTerms2ndDer[12]
                                + solDs[3]*metricTerms2ndDer[13] + solDt[3]*metricTerms2ndDer[14];

      /* Compute the Cartesian gradients of the pressure, velocity components,
         static energy and dynamic viscosity. */
      const su2double dpdx = Gamma_Minus_One*(drEdx + kinEnergy*drhodx
                           -                  u*drudx - v*drvdx - w*drwdx);
      const su2double dpdy = Gamma_Minus_One*(drEdy + kinEnergy*drhody
                           -                  u*drudy - v*drvdy - w*drwdy);
      const su2double dpdz = Gamma_Minus_One*(drEdz + kinEnergy*drhodz
                           -                  u*drudz - v*drvdz - w*drwdz);

      const su2double dudx = rhoInv*(drudx - u*drhodx);
      const su2double dudy = rhoInv*(drudy - u*drhody);
      const su2double dudz = rhoInv*(drudz - u*drhodz);

      const su2double dvdx = rhoInv*(drvdx - v*drhodx);
      const su2double dvdy = rhoInv*(drvdy - v*drhody);
      const su2double dvdz = rhoInv*(drvdz - v*drhodz);

      const su2double dwdx = rhoInv*(drwdx - w*drhodx);
      const su2double dwdy = rhoInv*(drwdy - w*drhody);
      const su2double dwdz = rhoInv*(drwdz - w*drhodz);

      const su2double dedx = rhoInv*(drEdx - TotalEnergy*drhodx) - u*dudx - v*dvdx - w*dwdx;
      const su2double dedy = rhoInv*(drEdy - TotalEnergy*drhody) - u*dudy - v*dvdy - w*dwdy;
      const su2double dedz = rhoInv*(drEdz - TotalEnergy*drhodz) - u*dudz - v*dvdz - w*dwdz;

      const su2double dViscLamdx = CvInv*dedx*dViscLamdT;
      const su2double dViscLamdy = CvInv*dedy*dViscLamdT;
      const su2double dViscLamdz = CvInv*dedz*dViscLamdT;

      /*--- Compute the second derivatives of the velocity components. ---*/
      const su2double d2udxdx = rhoInv*(d2rudxdx - u*d2rhodxdx
                              +         2.0*rhoInv*drhodx*(u*drhodx - drudx));
      const su2double d2udydy = rhoInv*(d2rudydy - u*d2rhodydy
                              +         2.0*rhoInv*drhody*(u*drhody - drudy));
      const su2double d2udzdz = rhoInv*(d2rudzdz - u*d2rhodzdz
                              +         2.0*rhoInv*drhodz*(u*drhodz - drudz));
      const su2double d2udxdy = rhoInv*(d2rudxdy - u*d2rhodxdy
                              +         rhoInv*(drhodx*(u*drhody - drudy)
                              +                 drhody*(u*drhodx - drudx)));
      const su2double d2udxdz = rhoInv*(d2rudxdz - u*d2rhodxdz
                              +         rhoInv*(drhodx*(u*drhodz - drudz)
                              +                 drhodz*(u*drhodx - drudx)));
      const su2double d2udydz = rhoInv*(d2rudydz - u*d2rhodydz
                              +         rhoInv*(drhody*(u*drhodz - drudz)
                              +                 drhodz*(u*drhody - drudy)));

      const su2double d2vdxdx = rhoInv*(d2rvdxdx - v*d2rhodxdx
                              +         2.0*rhoInv*drhodx*(v*drhodx - drvdx));
      const su2double d2vdydy = rhoInv*(d2rvdydy - v*d2rhodydy
                              +         2.0*rhoInv*drhody*(v*drhody - drvdy));
      const su2double d2vdzdz = rhoInv*(d2rvdzdz - v*d2rhodzdz
                              +         2.0*rhoInv*drhodz*(v*drhodz - drvdz));
      const su2double d2vdxdy = rhoInv*(d2rvdxdy - v*d2rhodxdy
                              +         rhoInv*(drhodx*(v*drhody - drvdy)
                              +                 drhody*(v*drhodx - drvdx)));
      const su2double d2vdxdz = rhoInv*(d2rvdxdz - v*d2rhodxdz
                              +         rhoInv*(drhodx*(v*drhodz - drvdz)
                              +                 drhodz*(v*drhodx - drvdx)));
      const su2double d2vdydz = rhoInv*(d2rvdydz - v*d2rhodydz
                              +         rhoInv*(drhody*(v*drhodz - drvdz)
                              +                 drhodz*(v*drhody - drvdy)));

      const su2double d2wdxdx = rhoInv*(d2rwdxdx - w*d2rhodxdx
                              +         2.0*rhoInv*drhodx*(w*drhodx - drwdx));
      const su2double d2wdydy = rhoInv*(d2rwdydy - w*d2rhodydy
                              +         2.0*rhoInv*drhody*(w*drhody - drwdy));
      const su2double d2wdzdz = rhoInv*(d2rwdzdz - w*d2rhodzdz
                              +         2.0*rhoInv*drhodz*(w*drhodz - drwdz));
      const su2double d2wdxdy = rhoInv*(d2rwdxdy - w*d2rhodxdy
                              +         rhoInv*(drhodx*(w*drhody - drwdy)
                              +                 drhody*(w*drhodx - drwdx)));
      const su2double d2wdxdz = rhoInv*(d2rwdxdz - w*d2rhodxdz
                              +         rhoInv*(drhodx*(w*drhodz - drwdz)
                              +                 drhodz*(w*drhodx - drwdx)));
      const su2double d2wdydz = rhoInv*(d2rwdydz - w*d2rhodydz
                              +         rhoInv*(drhody*(w*drhodz - drwdz)
                              +                 drhodz*(w*drhody - drwdy)));

      /* Compute the second derivatives of the static energy. Note that this
         term appears in the heat flux and therefore only the pure second
         derivatives are needed. Hence, the cross-derivatives are omitted. */
      const su2double d2edxdx = rhoInv*(d2rEdxdx - TotalEnergy*d2rhodxdx
                              +         2.0*rhoInv*drhodx*(TotalEnergy*drhodx - drEdx))
                              -         u*d2udxdx - dudx*dudx - v*d2vdxdx - dvdx*dvdx
                              -         w*d2wdxdx - dwdx*dwdx;
      const su2double d2edydy = rhoInv*(d2rEdydy - TotalEnergy*d2rhodydy
                              +         2.0*rhoInv*drhody*(TotalEnergy*drhody - drEdy))
                              -         u*d2udydy - dudy*dudy - v*d2vdydy - dvdy*dvdy
                              -         w*d2wdydy - dwdy*dwdy;
      const su2double d2edzdz = rhoInv*(d2rEdzdz - TotalEnergy*d2rhodzdz
                              +         2.0*rhoInv*drhodz*(TotalEnergy*drhodz - drEdz))
                              -         u*d2udzdz - dudz*dudz - v*d2vdzdz - dvdz*dvdz
                              -         w*d2wdzdz - dwdz*dwdz;

      /*--- If an SGS model is used the eddy viscosity and its spatial
            derivatives must be computed. ---*/
      su2double ViscosityTurb = 0.0;
      su2double dViscTurbdx = 0.0, dViscTurbdy = 0.0, dViscTurbdz = 0.0;

      if( SGSModelUsed ) {
        const su2double dist = elem->wallDistance[i];
        ViscosityTurb = SGSModel->ComputeEddyViscosity_3D(rho, dudx, dudy, dudz,
                                                          dvdx, dvdy, dvdz, dwdx,
                                                          dwdy, dwdz, lenScale, dist);

        SGSModel->ComputeGradEddyViscosity_3D(rho, drhodx, drhody, drhodz, dudx, dudy,
                                              dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz,
                                              d2udxdx, d2udydy, d2udzdz, d2udxdy,
                                              d2udxdz, d2udydz, d2vdxdx, d2vdydy,
                                              d2vdzdz, d2vdxdy, d2vdxdz, d2vdydz,
                                              d2wdxdx, d2wdydy, d2wdzdz, d2wdxdy,
                                              d2wdxdz, d2wdydz, lenScale, dist,
                                              dViscTurbdx, dViscTurbdy, dViscTurbdz);
      }

      /*--- Compute the total viscosity, the total heat conductivity and their
            gradients. Note that the heat conductivity is divided by the Cv,
            because gradients of internal energy are computed and not temperature. ---*/
      const su2double Viscosity = ViscosityLam + ViscosityTurb;
      const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                              + ViscosityTurb*factHeatFlux_Turb;

      const su2double dViscDx = dViscLamdx + dViscTurbdx;
      const su2double dViscDy = dViscLamdy + dViscTurbdy;
      const su2double dViscDz = dViscLamdz + dViscTurbdz;

      const su2double dkOverCvdx = dViscLamdx *factHeatFlux_Lam
                                 + dViscTurbdx*factHeatFlux_Turb;
      const su2double dkOverCvdy = dViscLamdy *factHeatFlux_Lam
                                 + dViscTurbdy*factHeatFlux_Turb;
      const su2double dkOverCvdz = dViscLamdz *factHeatFlux_Lam
                                 + dViscTurbdz*factHeatFlux_Turb;

      /* Abbreviations, which make it easier to compute the divergence term. */
      const su2double abv1 = drudx + drvdy + drwdz;
      const su2double abv2 = u*drhodx + v*drhody + w*drhodz;
      const su2double abv3 = u*(drEdx + dpdx) + v*(drEdy + dpdy) + w*(drEdz + dpdz);
      const su2double abv4 = dudx + dvdy + dwdz;

      /*--- Compute the divergence of the grid velocity.
            SET TO ZERO FOR NOW. THIS IS NOT CORRECT!!!!. ---*/
      const su2double divGridVel = 0.0;

      /* Set the pointer to store the divergence terms for this integration
         point and compute these terms, multiplied by the integration weight
         and Jacobian. */
      const su2double weightJac = weights[i]*Jac;
      su2double *divFluxInt     = divFlux + offInt;

      divFluxInt[0] = weightJac*(abv1 - rho*divGridVel - gridVel[0]*drhodx
                    -            gridVel[1]*drhody - gridVel[2]*drhodz);
      divFluxInt[1] = weightJac*(dpdx + u*(abv1-abv2) - lambdaOverMu*abv4*dViscDx
                    +            u*drudx + v*drudy + w*drudz
                    -            lambdaOverMu*Viscosity*(d2udxdx + d2vdxdy + d2wdxdz)
                    -            Viscosity*(2.0*d2udxdx + d2udydy + d2vdxdy + d2udzdz + d2wdxdz)
                    -            2.0*dViscDx*dudx - dViscDy*(dudy+dvdx) - dViscDz*(dudz+dwdx)
                    -            ru*divGridVel
                    -            gridVel[0]*drudx - gridVel[1]*drudy - gridVel[2]*drudz);
      divFluxInt[2] = weightJac*(dpdy + v*(abv1-abv2) - lambdaOverMu*abv4*dViscDy
                    +            u*drvdx + v*drvdy + w*drvdz
                    -            lambdaOverMu*Viscosity*(d2udxdy + d2vdydy + d2wdydz)
                    -            Viscosity*(d2udxdy + d2vdxdx + 2.0*d2vdydy + d2vdzdz + d2wdydz)
                    -            dViscDx*(dudy + dvdx) - 2.0*dViscDy*dvdy - dViscDz*(dvdz+dwdy)
                    -            rv*divGridVel
                    -            gridVel[0]*drvdx - gridVel[1]*drvdy - gridVel[2]*drvdz);
      divFluxInt[3] = weightJac*(dpdz + w*(abv1-abv2) - lambdaOverMu*abv4*dViscDz
                    +            u*drwdx + v*drwdy + w*drwdz
                    -            lambdaOverMu*Viscosity*(d2udxdz + d2vdydz + d2wdzdz)
                    -            Viscosity*(d2udxdz + d2wdxdx + d2vdydz + d2wdydy + 2.0*d2wdzdz)
                    -            dViscDx*(dudz+dwdx) - dViscDy*(dvdz+dwdy) - 2.0*dViscDz*dwdz
                    -            rw*divGridVel
                    -            gridVel[0]*drwdx - gridVel[1]*drwdy - gridVel[2]*drwdz);
      divFluxInt[4] = weightJac*(abv3 + Htot*(abv1 - abv2)
                    -            abv4*lambdaOverMu*(Viscosity*abv4 + u*dViscDx + v*dViscDy + w*dViscDz)
                    -            dkOverCvdx*dedx - dkOverCvdy*dedy - dkOverCvdz*dedz
                    -            kOverCv*(d2edxdx + d2edydy + d2edzdz)
                    -            (Viscosity*dudx + u*dViscDx)*2.0*dudx
                    -            (Viscosity*dvdy + v*dViscDy)*2.0*dvdy
                    -            (Viscosity*dwdz + w*dViscDz)*2.0*dwdz
                    -            (Viscosity*dudy + u*dViscDy + Viscosity*dvdx + v*dViscDx)*(dudy + dvdx)
                    -            (Viscosity*dudz + u*dViscDz + Viscosity*dwdx + w*dViscDx)*(dudz + dwdx)
                    -            (Viscosity*dvdz + v*dViscDz + Viscosity*dwdy + w*dViscDy)*(dvdz + dwdy)
                    -            Viscosity*u*(d2udxdx+d2udydy+d2udzdz + (1.0+lambdaOverMu)*(d2udxdx+d2vdxdy+d2wdxdz))
                    -            Viscosity*v*(d2vdxdx+d2vdydy+d2vdzdz + (1.0+lambdaOverMu)*(d2udxdy+d2vdydy+d2wdydz))
                    -            Viscosity*w*(d2wdxdx+d2wdydy+d2wdzdz + (1.0+lambdaOverMu)*(d2udxdz+d2vdydz+d2wdzdz))
                    -            rE*divGridVel
                    -            gridVel[0]*drEdx - gridVel[1]*drEdy - gridVel[2]*drEdz);

      /* Add the body force to the flux divergence for the momentum and energy
         equation. Note that the source terms are multiplied with minus the
         integration weight in order to be consistent with the formulation of
         the residual. Also note that for the energy source term the absolute
         velocity must be taken and not the relative. */
      divFluxInt[1] -= weightJac*bodyForceX;
      divFluxInt[2] -= weightJac*bodyForceY;
      divFluxInt[3] -= weightJac*bodyForceZ;
      divFluxInt[4] -= weightJac*(u*bodyForceX + v*bodyForceY + w*bodyForceZ);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the source terms of the manufactured solution to the divergence---*/
  /*--- of the fluxes, if a manufactured solution is used.                 ---*/
  /*--------------------------------------------------------------------------*/

  if( VerificationSolution ) {
    if( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Loop over the number of entities that are treated simultaneously. */
      for(unsigned short simul=0; simul<nSimul; ++simul) {

        /*--- Loop over the integration points of the element. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /* Set the pointers for this integration point. */
          const unsigned short offInt  = i*NPad + simul*nVar;
          su2double       *divFluxInt = divFlux + offInt;

          /* Set the pointer to the coordinates in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *coor = elem->coorIntegrationPoints.data() + i*nDim;

          /* Easier storage of the metric terms in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
          const su2double weightJac    = weights[i]*metricTerms[0];

          /* Compute the source terms of the manufactured solution.
             THIS IS A TEMPORARY IMPLEMENTATION. FOR AN ACTUAL TIME ACCURATE
             SIMULATION THE CORRECT TIME MUST BE GIVEN TO THIS FUNCTION. */
          su2double sourceMan[5];
          VerificationSolution->GetMMSSourceTerm(coor, 0.0, sourceMan);

          /* Add the source terms to the flux divergence. Note that the source
             terms are multiplied with minus the integration weight in order
             to be consistent with the formulation of the residual. */
          divFluxInt[0] -= weightJac*sourceMan[0];
          divFluxInt[1] -= weightJac*sourceMan[1];
          divFluxInt[2] -= weightJac*sourceMan[2];
          divFluxInt[3] -= weightJac*sourceMan[3];
          divFluxInt[4] -= weightJac*sourceMan[4];
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the residual in the DOFs, which is the matrix product of   ---*/
  /*--- basisFunctionsIntTrans and divFlux.                                ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nDOFs, NPad, nInt, basisFunctionsIntTrans, divFlux, res, config);
}

void CFEM_DG_NSSolver::Shock_Capturing_DG(CConfig             *config,
                                          const unsigned long elemBeg,
                                          const unsigned long elemEnd,
                                          su2double           *workArray) {

  /*--- Run shock capturing algorithm ---*/
  switch( config->GetKind_FEM_DG_Shock() ) {
    case FEM_SHOCK_CAPTURING_DG::NONE:
      break;
    case FEM_SHOCK_CAPTURING_DG::PERSSON:
      Shock_Capturing_DG_Persson(elemBeg, elemEnd, workArray);
      break;
  }

}
void CFEM_DG_NSSolver::Shock_Capturing_DG_Persson(const unsigned long elemBeg,
                                                  const unsigned long elemEnd,
                                                  su2double           *workArray) {

  /*--- Dummy variable for storing shock sensor value temporarily ---*/
  su2double sensorVal, sensorLowerBound, machNorm, machMax;
  su2double DensityInv, Velocity2, StaticEnergy, SoundSpeed2, Velocity2Rel;

  bool shockExist;
  unsigned short nDOFsPm1;       // Number of DOFs up to polynomial degree p-1

  /*--- Loop over the given range of elements to sense the shock. If shock exists,
        add artificial viscosity for DG FEM formulation to the residual.  ---*/
  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /* Get the data from the corresponding standard element. */
    const unsigned short ind          = volElem[l].indStandardElement;
    const unsigned short nDOFs        = volElem[l].nDOFsSol;
    const unsigned short VTK_TypeElem = volElem[l].VTK_Type;
    const unsigned short nPoly        = standardElementsSol[ind].GetNPoly();
    const su2double *matVanderInv     = standardElementsSol[ind].GetMatVandermondeInv();

    /*----------------------------------------------------------------------------*/
    /*--- Step 1: Calculate the number of DOFs up to polynomial degree p-1.    ---*/
    /*----------------------------------------------------------------------------*/

    nDOFsPm1 = 0;
    switch( VTK_TypeElem ) {
      case TRIANGLE:
        nDOFsPm1 = nPoly*(nPoly+1)/2;
        break;
      case QUADRILATERAL:
        nDOFsPm1 = nPoly*nPoly;
        break;
      case TETRAHEDRON:
        nDOFsPm1 = nPoly*(nPoly+1)*(nPoly+2)/6;
        break;
      case PYRAMID:
        nDOFsPm1 = nPoly*(nPoly+1)*(2*nPoly+1)/6;
        break;
      case PRISM:
        nDOFsPm1 = nPoly*nPoly*(nPoly+1)/2;
        break;
      case HEXAHEDRON:
        nDOFsPm1 = nPoly*nPoly*nPoly;
        break;
    }

    /*---------------------------------------------------------------------*/
    /*--- Step 2: Calculate the shock sensor value for this element.    ---*/
    /*---------------------------------------------------------------------*/

    /* Initialize dummy variable for this volume element */
    sensorVal = 0;
    machMax = -1;
    shockExist = false;
    sensorLowerBound = 1.e15;

    /* Easier storage of the solution variables for this element. */
    const unsigned short timeLevel = volElem[l].timeLevel;
    const su2double *solDOFs = VecWorkSolDOFs[timeLevel].data()
                             + nVar*volElem[l].offsetDOFsSolThisTimeLevel;

    /* Temporary storage of mach number for DOFs in this element. */
    su2double *machSolDOFs = workArray;
    su2double *vecTemp     = machSolDOFs + nDOFs;

    /* Calculate primitive variables and mach number for DOFs in this element.
       Also, track the maximum mach number in this element. */
    for(unsigned short iInd=0; iInd<nDOFs; ++iInd) {

      const su2double *sol     = solDOFs + iInd*nVar;
      const su2double *gridVel = volElem[l].gridVelocitiesSolDOFs.data() + iInd*nDim;
      DensityInv = 1.0/sol[0];
      Velocity2 = 0.0;
      Velocity2Rel = 0.0;
      for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
        const su2double vel    = sol[iDim]*DensityInv;
        const su2double velRel = vel - gridVel[iDim-1];
        Velocity2    += vel*vel;
        Velocity2Rel += velRel*velRel;
      }

      StaticEnergy = sol[nDim+1]*DensityInv - 0.5*Velocity2;

      FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
      SoundSpeed2 = FluidModel->GetSoundSpeed2();
      machSolDOFs[iInd] = sqrt( Velocity2Rel/SoundSpeed2 );
      machMax = max(machSolDOFs[iInd],machMax);
    }

    /* Change the solution coefficients to modal form from nodal form */
    for(unsigned short i=0; i<nDOFs; ++i) {
      vecTemp[i] = 0.0;
      for (unsigned short j=0; j<nDOFs; ++j)
        vecTemp[i] += matVanderInv[i+j*nDOFs]*machSolDOFs[j];
    }

    /* Get the L2 norm of solution coefficients for the highest polynomial order. */
    for(unsigned short i=nDOFsPm1; i<nDOFs; ++i) {
        sensorVal += vecTemp[i]*vecTemp[i];
    }

    /* If the maximum mach number is greater than 1.0, try to calculate the shockSensorValue.
       Otherwise, assign default value. */
    if ( machMax > 1.0) {
      // !!!!!Threshold value for sensorVal should be further investigated
      if(sensorVal > 1.e-15) {
        machNorm = 0.0;

        /*--- Get L2 norm square of vecTemp ---*/
        for (unsigned short i=0; i<nDOFs; ++i) {
          machNorm += vecTemp[i]*vecTemp[i];
        }
        if (machNorm < 1.e-15) {
          // This should not happen
          volElem[l].shockSensorValue = 1000.0;
        }
        else {
          volElem[l].shockSensorValue = log(sensorVal/machNorm);
          shockExist = true;
        }
      }
      else {
        // There is no shock in this element
        volElem[l].shockSensorValue = -1000.0;
      }
    }
    else {
      volElem[l].shockSensorValue = -1000.0;
    }

    /*---------------------------------------------------------------------*/
    /*--- Step 3: Determine artificial viscosity for this element.      ---*/
    /*---------------------------------------------------------------------*/
    if (shockExist) {
      // Following if-else clause is purely empirical from NACA0012 case.
      // Need to develop thorough method for general problems.
      switch ( nPoly ) {
        case 1:  sensorLowerBound =  -6.0; break;
        case 2:  sensorLowerBound = -12.0; break;
        case 3:  sensorLowerBound = -12.0; break;
        case 4:  sensorLowerBound = -17.0; break;
        default: sensorLowerBound = -17.0; break;
      }

      // Assign artificial viscosity based on shockSensorValue
      if ( volElem[l].shockSensorValue > sensorLowerBound ) {
        // Following value is initial guess.
        volElem[l].shockArtificialViscosity = 1.e-10;
      }
      else {
        volElem[l].shockArtificialViscosity = 0.0;
      }
    }
    else {
      volElem[l].shockArtificialViscosity = 0.0;
    }
  }
}

void CFEM_DG_NSSolver::Volume_Residual(CConfig             *config,
                                       const unsigned long elemBeg,
                                       const unsigned long elemEnd,
                                       su2double           *workArray) {

  /*--- Determine whether a body force term is present. ---*/
  bool body_force = config->GetBody_Force();
  const su2double *body_force_vector = body_force ? config->GetBody_Force_Vector() : nullptr;

  /*--- Get the physical time if necessary. ---*/
  su2double time = 0.0;
  if (config->GetTime_Marching() != TIME_MARCHING::STEADY) time = config->GetPhysicalTime();

  /* Constant factor present in the heat flux vector. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /* Determine the number of elements that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nElemSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /* Store the number of metric points per integration point, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

  /*--- Loop over the given element range to compute the contribution of the
        volume integral in the DG FEM formulation to the residual. Multiple
        elements are treated simultaneously to improve the performance
        of the matrix multiplications. As a consequence, the update of the
        counter l happens at the end of this loop section. ---*/
  for(unsigned long l=elemBeg; l<elemEnd;) {

    /* Determine the end index for this chunk of elements and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(volElem, l, elemEnd, nElemSimul, nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the required data from the corresponding standard element. */
    const unsigned short nInt            = standardElementsSol[ind].GetNIntegration();
    const unsigned short nDOFs           = volElem[l].nDOFsSol;
    const su2double *matBasisInt         = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
    const su2double *matDerBasisIntTrans = standardElementsSol[ind].GetDerMatBasisFunctionsIntTrans();
    const su2double *matBasisIntTrans    = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
    const su2double *weights             = standardElementsSol[ind].GetWeightsIntegration();

    unsigned short nPoly = standardElementsSol[ind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solDOFs       = workArray;
    su2double *sources       = solDOFs       + nDOFs*NPad;
    su2double *solAndGradInt = sources       + nInt *NPad;
    su2double *fluxes        = solAndGradInt + nInt *NPad*(nDim+1);

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Determine the solution variables and their gradients     ---*/
    /*---         w.r.t. the parametric coordinates in the integration     ---*/
    /*---         points of the element.                                   ---*/
    /*------------------------------------------------------------------------*/

    /* Copy the solution of the DOFs into solDOFs for this chunk of elements. */
    for(unsigned short ll=0; ll<llEnd; ++ll) {

      /* Easier storage of the solution variables for this element. */
      const unsigned long lInd = l + ll;
      const unsigned short timeLevel = volElem[lInd].timeLevel;
      const su2double *solDOFsElem = VecWorkSolDOFs[timeLevel].data()
                                   + nVar*volElem[lInd].offsetDOFsSolThisTimeLevel;

      /* Loop over the DOFs and copy the data. */
      const unsigned short llNVar = ll*nVar;
      for(unsigned short i=0; i<nDOFs; ++i)
        for(unsigned short mm=0; mm<nVar; ++mm)
          solDOFs[i*NPad+llNVar+mm] = solDOFsElem[i*nVar+mm];
    }

    /* Call the general function to carry out the matrix product to determine
       the solution and gradients in the integration points of the chunk
       of elements. */
    blasFunctions->gemm(nInt*(nDim+1), NPad, nDOFs, matBasisInt, solDOFs, solAndGradInt, config);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the total fluxes (inviscid fluxes minus the      ---*/
    /*---         viscous fluxes), multiplied by minus the integration     ---*/
    /*---         weight, in the integration points.                       ---*/
    /*------------------------------------------------------------------------*/

    /* Determine the offset between the solution variables and the r-derivatives,
       which is also the offset between the r- and s-derivatives and the offset
       between s- and t-derivatives. */
    const unsigned short offDeriv = NPad*nInt;

    /* Make a distinction between two and three space dimensions
        in order to have the most efficient code. */
    switch( nDim ) {

      case 2: {

        /* 2D simulation. Loop over the chunk of elements and loop over the
           integration points of the elements to compute the fluxes. */
        for(unsigned short ll=0; ll<llEnd; ++ll) {
          const unsigned short llNVar = ll*nVar;
          const unsigned long  lInd   = l + ll;

          for(unsigned short i=0; i<nInt; ++i) {
            const unsigned short iNPad = i*NPad;

            /* Easier storage of the metric terms and grid velocities in this
               integration point and compute the inverse of the Jacobian. */
            const su2double *metricTerms = volElem[lInd].metricTerms.data()
                                         + i*nMetricPerPoint;
            const su2double *gridVel     = volElem[lInd].gridVelocities.data() + i*nDim;
            const su2double Jac          = metricTerms[0];
            const su2double JacInv       = 1.0/Jac;

            /* Compute the true metric terms in this integration point. */
            const su2double drdx = JacInv*metricTerms[1];
            const su2double drdy = JacInv*metricTerms[2];

            const su2double dsdx = JacInv*metricTerms[3];
            const su2double dsdy = JacInv*metricTerms[4];

            /* Compute the metric terms multiplied by minus the integration weight.
               The minus sign comes from the integration by parts in the weak
               formulation. */
            const su2double wDrdx = -weights[i]*metricTerms[1];
            const su2double wDrdy = -weights[i]*metricTerms[2];

            const su2double wDsdx = -weights[i]*metricTerms[3];
            const su2double wDsdy = -weights[i]*metricTerms[4];

            /* Easier storage of the location where the solution data of this
               integration point starts. */
            const su2double *sol    = solAndGradInt + iNPad + llNVar;
            const su2double *dSolDr = sol    + offDeriv;
            const su2double *dSolDs = dSolDr + offDeriv;

            /*--- Compute the Cartesian gradients of the independent solution
                  variables from the gradients in parametric coordinates and the
                  metric terms in this integration point. ---*/
            const su2double dRhoDx  = dSolDr[0]*drdx + dSolDs[0]*dsdx;
            const su2double dRhoUDx = dSolDr[1]*drdx + dSolDs[1]*dsdx;
            const su2double dRhoVDx = dSolDr[2]*drdx + dSolDs[2]*dsdx;
            const su2double dRhoEDx = dSolDr[3]*drdx + dSolDs[3]*dsdx;

            const su2double dRhoDy  = dSolDr[0]*drdy + dSolDs[0]*dsdy;
            const su2double dRhoUDy = dSolDr[1]*drdy + dSolDs[1]*dsdy;
            const su2double dRhoVDy = dSolDr[2]*drdy + dSolDs[2]*dsdy;
            const su2double dRhoEDy = dSolDr[3]*drdy + dSolDs[3]*dsdy;

            /*--- Compute the velocities and static energy in this integration point. ---*/
            const su2double rhoInv       = 1.0/sol[0];
            const su2double u            = sol[1]*rhoInv;
            const su2double v            = sol[2]*rhoInv;
            const su2double TotalEnergy  = sol[3]*rhoInv;
            const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v);

            /*--- Compute the Cartesian gradients of the velocities and static energy
                  in this integration point and also the divergence of the velocity. ---*/
            const su2double dudx = rhoInv*(dRhoUDx - u*dRhoDx);
            const su2double dudy = rhoInv*(dRhoUDy - u*dRhoDy);

            const su2double dvdx = rhoInv*(dRhoVDx - v*dRhoDx);
            const su2double dvdy = rhoInv*(dRhoVDy - v*dRhoDy);

            const su2double dStaticEnergyDx = rhoInv*(dRhoEDx - TotalEnergy*dRhoDx)
                                            - u*dudx - v*dvdx;
            const su2double dStaticEnergyDy = rhoInv*(dRhoEDy - TotalEnergy*dRhoDy)
                                            - u*dudy - v*dvdy;

            const su2double divVel = dudx + dvdy;

            /*--- Compute the pressure and the laminar viscosity. ---*/
            FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
            const su2double Pressure     = FluidModel->GetPressure();
            const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

            /*--- If an SGS model is used the eddy viscosity must be computed. ---*/
            su2double ViscosityTurb = 0.0;
            if( SGSModelUsed ) {
              const su2double lenScale = volElem[lInd].lenScale/nPoly;
              ViscosityTurb = SGSModel->ComputeEddyViscosity_2D(sol[0], dudx, dudy,
                                                                dvdx, dvdy, lenScale,
                                                                volElem[lInd].wallDistance[i]);
            }

            /* Compute the total viscosity and heat conductivity. Note that the heat
               conductivity is divided by the Cv, because gradients of internal energy
               are computed and not temperature. */
            const su2double Viscosity = ViscosityLam + ViscosityTurb;
            const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                                    + ViscosityTurb*factHeatFlux_Turb;

            /*--- Set the value of the second viscosity and compute the divergence
                  term in the viscous normal stresses. ---*/
            const su2double lambda     = -TWO3*Viscosity;
            const su2double lamDivTerm =  lambda*divVel;

            /*--- Compute the viscous stress tensor and minus the heatflux vector. ---*/
            const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
            const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
            const su2double tauxy = Viscosity*(dudy + dvdx);

            const su2double qx = kOverCv*dStaticEnergyDx;
            const su2double qy = kOverCv*dStaticEnergyDy;

            /* Compute the relative velocities w.r.t. the grid. */
            const su2double uRel = u - gridVel[0];
            const su2double vRel = v - gridVel[1];

            /* Compute the viscous normal stress minus the pressure. */
            const su2double tauxxMP = tauxx - Pressure;
            const su2double tauyyMP = tauyy - Pressure;

            /* Set the pointer for the fluxes in this integration point. */
            su2double *flux = fluxes + nDim*iNPad + llNVar;

            /*--- Fluxes in r-direction. */
            const su2double Ur = uRel*wDrdx + vRel*wDrdy;

            flux[0] = sol[0]*Ur;
            flux[1] = sol[1]*Ur - tauxxMP*wDrdx - tauxy*wDrdy;
            flux[2] = sol[2]*Ur - tauxy*wDrdx - tauyyMP*wDrdy;
            flux[3] = sol[3]*Ur - (u*tauxxMP + v*tauxy + qx)*wDrdx
                                - (u*tauxy + v*tauyyMP + qy)*wDrdy;

            /*--- Fluxes in s-direction. */
            flux = flux + NPad;
            const su2double Us = uRel*wDsdx + vRel*wDsdy;

            flux[0] = sol[0]*Us;
            flux[1] = sol[1]*Us - tauxxMP*wDsdx - tauxy*wDsdy;
            flux[2] = sol[2]*Us - tauxy*wDsdx - tauyyMP*wDsdy;
            flux[3] = sol[3]*Us - (u*tauxxMP + v*tauxy + qx)*wDsdx
                                - (u*tauxy + v*tauyyMP + qy)*wDsdy;

            /*--- If needed, compute the body forces in this integration point.
                  Note that the source terms are multiplied with minus the
                  integration weight in order to be consistent with the
                  formulation of the residual. Note that for the energy source
                  term the absolute velocity must be taken and not the
                  relative. ---*/
            if( body_force ) {
              su2double *source         = sources + iNPad + llNVar;
              const su2double weightJac = weights[i]*Jac;

              source[0] =  0.0;
              source[1] = -weightJac*body_force_vector[0];
              source[2] = -weightJac*body_force_vector[1];
              source[3] = -weightJac*(u*body_force_vector[0] + v*body_force_vector[1]);
            }
          }
        }

        break;
      }

      /*----------------------------------------------------------------------*/

      case 3: {

        /* 3D simulation. Loop over the chunk of elements and loop over the
           integration points of the elements to compute the fluxes. */
        for(unsigned short ll=0; ll<llEnd; ++ll) {
          const unsigned short llNVar = ll*nVar;
          const unsigned long  lInd   = l + ll;

          for(unsigned short i=0; i<nInt; ++i) {
            const unsigned short iNPad = i*NPad;

            /* Easier storage of the metric terms and grid velocities in this
               integration point and compute the inverse of the Jacobian. */
            const su2double *metricTerms = volElem[lInd].metricTerms.data()
                                         + i*nMetricPerPoint;
            const su2double *gridVel     = volElem[lInd].gridVelocities.data() + i*nDim;
            const su2double Jac          = metricTerms[0];
            const su2double JacInv       = 1.0/Jac;

            /* Compute the true metric terms in this integration point. */
            const su2double drdx = JacInv*metricTerms[1];
            const su2double drdy = JacInv*metricTerms[2];
            const su2double drdz = JacInv*metricTerms[3];

            const su2double dsdx = JacInv*metricTerms[4];
            const su2double dsdy = JacInv*metricTerms[5];
            const su2double dsdz = JacInv*metricTerms[6];

            const su2double dtdx = JacInv*metricTerms[7];
            const su2double dtdy = JacInv*metricTerms[8];
            const su2double dtdz = JacInv*metricTerms[9];

            /* Compute the metric terms multiplied by minus the integration weight.
               The minus sign comes from the integration by parts in the weak
               formulation. */
            const su2double wDrdx = -weights[i]*metricTerms[1];
            const su2double wDrdy = -weights[i]*metricTerms[2];
            const su2double wDrdz = -weights[i]*metricTerms[3];

            const su2double wDsdx = -weights[i]*metricTerms[4];
            const su2double wDsdy = -weights[i]*metricTerms[5];
            const su2double wDsdz = -weights[i]*metricTerms[6];

            const su2double wDtdx = -weights[i]*metricTerms[7];
            const su2double wDtdy = -weights[i]*metricTerms[8];
            const su2double wDtdz = -weights[i]*metricTerms[9];

            /* Easier storage of the location where the solution data of this
               integration point starts. */
            const su2double *sol    = solAndGradInt + iNPad + llNVar;
            const su2double *dSolDr = sol    + offDeriv;
            const su2double *dSolDs = dSolDr + offDeriv;
            const su2double *dSolDt = dSolDs + offDeriv;

            /*--- Compute the Cartesian gradients of the independent solution
                  variables from the gradients in parametric coordinates and the
                  metric terms in this integration point. ---*/
            const su2double dRhoDx  = dSolDr[0]*drdx + dSolDs[0]*dsdx + dSolDt[0]*dtdx;
            const su2double dRhoUDx = dSolDr[1]*drdx + dSolDs[1]*dsdx + dSolDt[1]*dtdx;
            const su2double dRhoVDx = dSolDr[2]*drdx + dSolDs[2]*dsdx + dSolDt[2]*dtdx;
            const su2double dRhoWDx = dSolDr[3]*drdx + dSolDs[3]*dsdx + dSolDt[3]*dtdx;
            const su2double dRhoEDx = dSolDr[4]*drdx + dSolDs[4]*dsdx + dSolDt[4]*dtdx;

            const su2double dRhoDy  = dSolDr[0]*drdy + dSolDs[0]*dsdy + dSolDt[0]*dtdy;
            const su2double dRhoUDy = dSolDr[1]*drdy + dSolDs[1]*dsdy + dSolDt[1]*dtdy;
            const su2double dRhoVDy = dSolDr[2]*drdy + dSolDs[2]*dsdy + dSolDt[2]*dtdy;
            const su2double dRhoWDy = dSolDr[3]*drdy + dSolDs[3]*dsdy + dSolDt[3]*dtdy;
            const su2double dRhoEDy = dSolDr[4]*drdy + dSolDs[4]*dsdy + dSolDt[4]*dtdy;

            const su2double dRhoDz  = dSolDr[0]*drdz + dSolDs[0]*dsdz + dSolDt[0]*dtdz;
            const su2double dRhoUDz = dSolDr[1]*drdz + dSolDs[1]*dsdz + dSolDt[1]*dtdz;
            const su2double dRhoVDz = dSolDr[2]*drdz + dSolDs[2]*dsdz + dSolDt[2]*dtdz;
            const su2double dRhoWDz = dSolDr[3]*drdz + dSolDs[3]*dsdz + dSolDt[3]*dtdz;
            const su2double dRhoEDz = dSolDr[4]*drdz + dSolDs[4]*dsdz + dSolDt[4]*dtdz;

            /*--- Compute the velocities and static energy in this integration point. ---*/
            const su2double rhoInv       = 1.0/sol[0];
            const su2double u            = sol[1]*rhoInv;
            const su2double v            = sol[2]*rhoInv;
            const su2double w            = sol[3]*rhoInv;
            const su2double TotalEnergy  = sol[4]*rhoInv;
            const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v + w*w);

            /*--- Compute the Cartesian gradients of the velocities and static energy
                  in this integration point and also the divergence of the velocity. ---*/
            const su2double dudx = rhoInv*(dRhoUDx - u*dRhoDx);
            const su2double dudy = rhoInv*(dRhoUDy - u*dRhoDy);
            const su2double dudz = rhoInv*(dRhoUDz - u*dRhoDz);

            const su2double dvdx = rhoInv*(dRhoVDx - v*dRhoDx);
            const su2double dvdy = rhoInv*(dRhoVDy - v*dRhoDy);
            const su2double dvdz = rhoInv*(dRhoVDz - v*dRhoDz);

            const su2double dwdx = rhoInv*(dRhoWDx - w*dRhoDx);
            const su2double dwdy = rhoInv*(dRhoWDy - w*dRhoDy);
            const su2double dwdz = rhoInv*(dRhoWDz - w*dRhoDz);

            const su2double dStaticEnergyDx = rhoInv*(dRhoEDx - TotalEnergy*dRhoDx)
                                            - u*dudx - v*dvdx - w*dwdx;
            const su2double dStaticEnergyDy = rhoInv*(dRhoEDy - TotalEnergy*dRhoDy)
                                            - u*dudy - v*dvdy - w*dwdy;
            const su2double dStaticEnergyDz = rhoInv*(dRhoEDz - TotalEnergy*dRhoDz)
                                            - u*dudz - v*dvdz - w*dwdz;

            const su2double divVel = dudx + dvdy + dwdz;

            /*--- Compute the pressure and the laminar viscosity. ---*/
            FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
            const su2double Pressure     = FluidModel->GetPressure();
            const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

            /*--- If an SGS model is used the eddy viscosity must be computed. ---*/
            su2double ViscosityTurb = 0.0;
            if( SGSModelUsed ) {
              const su2double lenScale = volElem[lInd].lenScale/nPoly;
              ViscosityTurb = SGSModel->ComputeEddyViscosity_3D(sol[0], dudx, dudy, dudz,
                                                                dvdx, dvdy, dvdz, dwdx,
                                                                dwdy, dwdz, lenScale,
                                                                volElem[lInd].wallDistance[i]);
            }

            /* Compute the total viscosity and heat conductivity. Note that the heat
               conductivity is divided by the Cv, because gradients of internal energy
               are computed and not temperature. */
            const su2double Viscosity = ViscosityLam + ViscosityTurb;
            const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                                    + ViscosityTurb*factHeatFlux_Turb;

            /*--- Set the value of the second viscosity and compute the divergence
                  term in the viscous normal stresses. ---*/
            const su2double lambda     = -TWO3*Viscosity;
            const su2double lamDivTerm =  lambda*divVel;

            /*--- Compute the viscous stress tensor and minus the heatflux vector. ---*/
            const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
            const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
            const su2double tauzz = 2.0*Viscosity*dwdz + lamDivTerm;

            const su2double tauxy = Viscosity*(dudy + dvdx);
            const su2double tauxz = Viscosity*(dudz + dwdx);
            const su2double tauyz = Viscosity*(dvdz + dwdy);

            const su2double qx = kOverCv*dStaticEnergyDx;
            const su2double qy = kOverCv*dStaticEnergyDy;
            const su2double qz = kOverCv*dStaticEnergyDz;

            /* Compute the relative velocities w.r.t. the grid. */
            const su2double uRel = u - gridVel[0];
            const su2double vRel = v - gridVel[1];
            const su2double wRel = w - gridVel[2];

            /* Compute the viscous normal stress minus the pressure. */
            const su2double tauxxMP = tauxx - Pressure;
            const su2double tauyyMP = tauyy - Pressure;
            const su2double tauzzMP = tauzz - Pressure;

            /* Set the pointer for the fluxes in this integration point. */
            su2double *flux = fluxes + nDim*iNPad + llNVar;

            /*--- Fluxes in r-direction. */
            const su2double Ur = uRel*wDrdx + vRel*wDrdy + wRel*wDrdz;

            flux[0] = sol[0]*Ur;
            flux[1] = sol[1]*Ur - tauxxMP*wDrdx - tauxy*wDrdy - tauxz*wDrdz;
            flux[2] = sol[2]*Ur - tauxy*wDrdx - tauyyMP*wDrdy - tauyz*wDrdz;
            flux[3] = sol[3]*Ur - tauxz*wDrdx - tauyz*wDrdy - tauzzMP*wDrdz;
            flux[4] = sol[4]*Ur - (u*tauxxMP + v*tauxy + w*tauxz + qx)*wDrdx
                                - (u*tauxy + v*tauyyMP + w*tauyz + qy)*wDrdy
                                - (u*tauxz + v*tauyz + w*tauzzMP + qz)*wDrdz;

            /*--- Fluxes in s-direction. */
            flux = flux + NPad;
            const su2double Us = uRel*wDsdx + vRel*wDsdy + wRel*wDsdz;

            flux[0] = sol[0]*Us;
            flux[1] = sol[1]*Us - tauxxMP*wDsdx - tauxy*wDsdy - tauxz*wDsdz;
            flux[2] = sol[2]*Us - tauxy*wDsdx - tauyyMP*wDsdy - tauyz*wDsdz;
            flux[3] = sol[3]*Us - tauxz*wDsdx - tauyz*wDsdy - tauzzMP*wDsdz;
            flux[4] = sol[4]*Us - (u*tauxxMP + v*tauxy + w*tauxz + qx)*wDsdx
                                - (u*tauxy + v*tauyyMP + w*tauyz + qy)*wDsdy
                                - (u*tauxz + v*tauyz + w*tauzzMP + qz)*wDsdz;

            /*--- Fluxes in t-direction. */
            flux = flux + NPad;
            const su2double Ut = uRel*wDtdx + vRel*wDtdy + wRel*wDtdz;

            flux[0] = sol[0]*Ut;
            flux[1] = sol[1]*Ut - tauxxMP*wDtdx - tauxy*wDtdy - tauxz*wDtdz;
            flux[2] = sol[2]*Ut - tauxy*wDtdx - tauyyMP*wDtdy - tauyz*wDtdz;
            flux[3] = sol[3]*Ut - tauxz*wDtdx - tauyz*wDtdy - tauzzMP*wDtdz;
            flux[4] = sol[4]*Ut - (u*tauxxMP + v*tauxy + w*tauxz + qx)*wDtdx
                                - (u*tauxy + v*tauyyMP + w*tauyz + qy)*wDtdy
                                - (u*tauxz + v*tauyz + w*tauzzMP + qz)*wDtdz;

            /*--- If needed, compute the body forces in this integration point.
                  Note that the source terms are multiplied with minus the
                  integration weight in order to be consistent with the
                  formulation of the residual. Note that for the energy source
                  term the absolute velocity must be taken and not the
                  relative. ---*/
            if( body_force ) {
              su2double *source         = sources + iNPad + llNVar;
              const su2double weightJac = weights[i]*Jac;

              source[0] =  0.0;
              source[1] = -weightJac*body_force_vector[0];
              source[2] = -weightJac*body_force_vector[1];
              source[3] = -weightJac*body_force_vector[2];
              source[4] = -weightJac*(u*body_force_vector[0] + v*body_force_vector[1]
                        +             w*body_force_vector[2]);
            }
          }
        }

        break;
      }
    }

    /* Initialize addSourceTerms to body_force. The value of addSourceTerms
       is set to true when a manufactured solution is computed. */
    bool addSourceTerms = body_force;

    /* Check whether or not a manufactured solution is used. */
    if( VerificationSolution ) {
      if( VerificationSolution->IsManufacturedSolution() ) {

        /*--- For the manufactured solutions a source term must be added. If a
              standard source term has not been specified, initialize the source
              terms to zero and set addSourceTerms to true. ---*/
        addSourceTerms = true;
        if( !body_force ) {
          for(unsigned short i=0; i<(nInt*NPad); ++i)
            sources[i] = 0.0;
        }

        /*--- Loop over the chunk of elements and its integration points. ---*/
        for(unsigned short ll=0; ll<llEnd; ++ll) {
          const unsigned short llNVar = ll*nVar;
          const unsigned long  lInd   = l + ll;
          for(unsigned short i=0; i<nInt; ++i) {
            const unsigned short iNPad = i*NPad;

            /* Determine the integration weight multiplied by the Jacobian. */
            const su2double *metricTerms = volElem[lInd].metricTerms.data()
                                         + i*nMetricPerPoint;
            const su2double weightJac    = weights[i]*metricTerms[0];

            /* Set the pointer to the coordinates in this integration point and
               call the function to compute the source terms for the manufactured
               solution. */
            const su2double *coor = volElem[lInd].coorIntegrationPoints.data() + i*nDim;

            su2double sourceMan[5];

            VerificationSolution->GetMMSSourceTerm(coor, time, sourceMan);

            /*--- Subtract the source term of the manufactured solution, multiplied
                  by the appropriate weight, from the possibly earlier computed
                  source term. It is subtracted in order to be consistent with
                  the definition of the residual used in this code. ---*/
            su2double *source = sources + iNPad + llNVar;
            for(unsigned short k=0; k<nVar; ++k)
              source[k] -= weightJac*sourceMan[k];
          }
        }
      }
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the contribution to the residuals from the       ---*/
    /*---         integration over the volume element.                     ---*/
    /*------------------------------------------------------------------------*/

    /* Call the general function to carry out the matrix product.
       Use solDOFs as a temporary storage for the matrix product. */
    blasFunctions->gemm(nDOFs, NPad, nInt*nDim, matDerBasisIntTrans, fluxes, solDOFs, config);

    /* Add the contribution from the source terms, if needed. Use solAndGradInt
       as temporary storage for the matrix product. */
    if( addSourceTerms ) {

      /* Call the general function to carry out the matrix product. */
      blasFunctions->gemm(nDOFs, NPad, nInt, matBasisIntTrans, sources, solAndGradInt, config);

      /* Add the residuals due to source terms to the volume residuals */
      for(unsigned short i=0; i<(nDOFs*NPad); ++i)
        solDOFs[i] += solAndGradInt[i];
    }

    /* Loop over the elements in this chunk to store the residuals
       in the appropriate locations. */
    for(unsigned short ll=0; ll<llEnd; ++ll) {
      const unsigned short llNVar = ll*nVar;
      const unsigned long  lInd   = l + ll;

      /* Easier storage of the residuals for this volume element. */
      su2double *res = VecResDOFs.data() + nVar*volElem[lInd].offsetDOFsSolLocal;

      /* Loop over the DOFs and copy the data. */
      for(unsigned short i=0; i<nDOFs; ++i)
        for(unsigned short mm=0; mm<nVar; ++mm)
           res[i*nVar+mm] = solDOFs[i*NPad+llNVar+mm];
    }

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::ResidualFaces(CConfig             *config,
                                     const unsigned long indFaceBeg,
                                     const unsigned long indFaceEnd,
                                     unsigned long       &indResFaces,
                                     CNumerics           *numerics,
                                     su2double           *workArray) {

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of matching faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=indFaceBeg; l<indFaceEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(matchingInternalFaces, l, indFaceEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the necessary data from the standard face. */
    const unsigned short nInt = standardMatchingFacesSol[ind].GetNIntegration();
    const su2double *weights  = standardMatchingFacesSol[ind].GetWeightsIntegration();

    const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
    const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();

    const unsigned short nDOFsElem0  = standardMatchingFacesSol[ind].GetNDOFsElemSide0();
    const unsigned short nDOFsElem1  = standardMatchingFacesSol[ind].GetNDOFsElemSide1();

    /*--- Set the pointers for the local arrays. ---*/
    unsigned int sizeFluxes = nInt*nDim;
    sizeFluxes = NPad*max(sizeFluxes, (unsigned int) max(nDOFsElem0, nDOFsElem1));

    const unsigned int sizeGradSolInt = nInt*nDim*NPad;

    su2double *solIntL       = workArray;
    su2double *solIntR       = solIntL       + NPad*max(nInt, nDOFsElem0);
    su2double *viscosityIntL = solIntR       + NPad*max(nInt, nDOFsElem1);
    su2double *kOverCvIntL   = viscosityIntL + llEnd*nInt;
    su2double *viscosityIntR = kOverCvIntL   + llEnd*nInt;
    su2double *kOverCvIntR   = viscosityIntR + llEnd*nInt;
    su2double *gradSolInt    = kOverCvIntR   + llEnd*nInt;
    su2double *fluxes        = gradSolInt    + sizeGradSolInt;
    su2double *viscFluxes    = fluxes        + sizeFluxes;

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Compute the inviscid fluxes in the integration points of ---*/
    /*---         this chunk of matching faces.                            ---*/
    /*------------------------------------------------------------------------*/

    InviscidFluxesInternalMatchingFace(config, l, lEnd, NPad,
                                       solIntL, solIntR, fluxes, numerics);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the viscous fluxes in the integration points of  ---*/
    /*---         this chunk of matching faces and subtract them from the  ---*/
    /*---         already computed inviscid fluxes.                        ---*/
    /*------------------------------------------------------------------------*/

    /*---------------------------*/
    /*--- Side 0 of the face. ---*/
    /*---------------------------*/

    /* Get the matrix, which contains the derivatives w.r.t. the
       parametric coordinates of the basis functions. */
    const su2double *derBasisElem = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationSide0();

    /* Set the pointer solElem to viscFluxes. This is just for readability,
       as the same memory can be used for the storage of the solution of the
       DOFs of the element and the viscous fluxes to be computed. */
    su2double *solElem = viscFluxes;

    /*--- Loop over the faces in this chunk to set the conserved variables
          in the DOFs of the faces. */
    for(unsigned long ll=l; ll<lEnd; ++ll) {
      const unsigned short llRel  = ll - l;
      const unsigned short llNVar = llRel*nVar;

      /* Determine the time level of the face, which is the minimum
         of the levels of the adjacent elements. */
      const unsigned long elemID0 = matchingInternalFaces[ll].elemID0;
      const unsigned long elemID1 = matchingInternalFaces[ll].elemID1;
      const unsigned short timeLevelFace = min(volElem[elemID0].timeLevel,
                                               volElem[elemID1].timeLevel);

      /* Determine the offset that must be applied to access the correct data
         for this elemID0 in the working vector of the solution. If the element
         has the same time level as the face offsetDOFsSolThisTimeLevel is used,
         otherwise offsetDOFsSolPrevTimeLevel is the correct value. */
      unsigned long offset;
      if(volElem[elemID0].timeLevel == timeLevelFace) {
        offset = volElem[elemID0].offsetDOFsSolLocal
               - volElem[elemID0].offsetDOFsSolThisTimeLevel;
      }
      else {
        offset = volElem[elemID0].offsetDOFsSolLocal
               - volElem[elemID0].offsetDOFsSolPrevTimeLevel;
      }

      /*--- Store the solution of the DOFs of the elemID0 in the correct
            sequence, such that a generalized approach is possible. ---*/
      const unsigned long *DOFsElem = matchingInternalFaces[ll].DOFsSolElementSide0.data();
      for(unsigned short i=0; i<nDOFsElem0; ++i) {
        const su2double *solDOF = VecWorkSolDOFs[timeLevelFace].data()
                                + nVar*(DOFsElem[i] - offset);
        su2double       *sol    = solElem + NPad*i + llNVar;
        for(unsigned short mm=0; mm<nVar; ++mm)
          sol[mm] = solDOF[mm];
      }
    }

    /* Compute the gradients w.r.t. the parametric coordinates in the integration
       points. Call the general function to carry out the matrix product. */
    blasFunctions->gemm(nInt*nDim, NPad, nDOFsElem0, derBasisElem, solElem, gradSolInt, config);

    /*--- Loop over the faces in this chunk to compute the viscous flux
          vector for side 0. */
    for(unsigned long ll=l; ll<lEnd; ++ll) {

      /* Determine the ID of the adjacent element. */
      const unsigned short llRel = ll - l;
      const unsigned long elemID0 = matchingInternalFaces[ll].elemID0;

      /* Call the general function to compute the viscous flux in normal
         direction for side 0. */
      ViscousNormalFluxFace(&volElem[elemID0], llRel, nInt, NPad,
                            0.0, false, solIntL, gradSolInt,
                            matchingInternalFaces[ll].metricCoorDerivFace0.data(),
                            matchingInternalFaces[ll].metricNormalsFace.data(),
                            matchingInternalFaces[ll].wallDistance.data(),
                            viscFluxes, viscosityIntL, kOverCvIntL);
    }

    /*--- Subtract half of the viscous fluxes from the inviscid fluxes. The
          factor 0.5 comes from the fact that the average of the viscous fluxes
          of side 0 and side 1 must be taken in the DG-FEM formulation. ---*/
    for(unsigned short j=0; j<(NPad*nInt); ++j) fluxes[j] -= 0.5*viscFluxes[j];

    /*---------------------------*/
    /*--- Side 1 of the face. ---*/
    /*---------------------------*/

    /* Get the matrix, which contains the derivatives w.r.t. the
       parametric coordinates of the basis functions. */
    derBasisElem = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationSide1();

    /*--- Loop over the faces in this chunk to set the conserved variables
          in the DOFs of the faces. */
    for(unsigned long ll=l; ll<lEnd; ++ll) {
      const unsigned short llRel  = ll - l;
      const unsigned short llNVar = llRel*nVar;

      /* Determine the time level of the face, which is the minimum
         of the levels of the adjacent elements. */
      const unsigned long elemID0 = matchingInternalFaces[ll].elemID0;
      const unsigned long elemID1 = matchingInternalFaces[ll].elemID1;
      const unsigned short timeLevelFace = min(volElem[elemID0].timeLevel,
                                               volElem[elemID1].timeLevel);

      /* Determine the offset that must be applied to access the correct data
         for this elemID0 in the working vector of the solution. If the element
         has the same time level as the face offsetDOFsSolThisTimeLevel is used,
         otherwise offsetDOFsSolPrevTimeLevel is the correct value. */
      unsigned long offset;
      if(volElem[elemID1].timeLevel == timeLevelFace) {
        offset = volElem[elemID1].offsetDOFsSolLocal
               - volElem[elemID1].offsetDOFsSolThisTimeLevel;
      }
      else {
        offset = volElem[elemID1].offsetDOFsSolLocal
               - volElem[elemID1].offsetDOFsSolPrevTimeLevel;
      }

      /*--- Store the solution of the DOFs of the elemID1 in the correct
            sequence, such that a generalized approach is possible. ---*/
      const unsigned long *DOFsElem = matchingInternalFaces[ll].DOFsSolElementSide1.data();
      for(unsigned short i=0; i<nDOFsElem1; ++i) {
        const su2double *solDOF = VecWorkSolDOFs[timeLevelFace].data()
                                + nVar*(DOFsElem[i] - offset);
        su2double       *sol    = solElem + NPad*i + llNVar;
        for(unsigned short mm=0; mm<nVar; ++mm)
          sol[mm] = solDOF[mm];
      }
    }

    /* Compute the gradients w.r.t. the parametric coordinates in the integration
       points. Call the general function to carry out the matrix product. */
    blasFunctions->gemm(nInt*nDim, NPad, nDOFsElem1, derBasisElem, solElem, gradSolInt, config);

    /*--- Loop over the faces in this chunk to compute the viscous flux
          vector for side 1. */
    for(unsigned long ll=l; ll<lEnd; ++ll) {

      /* Determine the ID of the adjacent element. */
      const unsigned short llRel = ll - l;
      const unsigned long elemID1 = matchingInternalFaces[ll].elemID1;

      /* Call the general function to compute the viscous flux in normal
         direction for side 1. */
      ViscousNormalFluxFace(&volElem[elemID1], llRel, nInt, NPad,
                            0.0, false, solIntR, gradSolInt,
                            matchingInternalFaces[ll].metricCoorDerivFace1.data(),
                            matchingInternalFaces[ll].metricNormalsFace.data(),
                            matchingInternalFaces[ll].wallDistance.data(),
                            viscFluxes, viscosityIntR, kOverCvIntR);
    }

    /*--- Subtract half of the viscous fluxes from the inviscid fluxes. ---*/
    for(unsigned short j=0; j<(NPad*nInt); ++j) fluxes[j] -= 0.5*viscFluxes[j];

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the penalty terms in the integration points of   ---*/
    /*---         this chunk of matching faces and add them to the already ---*/
    /*---         stored inviscid and viscous fluxes.                      ---*/
    /*------------------------------------------------------------------------*/

    /* Get the required constant needed for the penalty terms. */
    const su2double ConstPenFace = standardMatchingFacesSol[ind].GetPenaltyConstant();

    /*--- Loop over the faces in this chunk to compute the penalty fluxes. ---*/
    for(unsigned long ll=l; ll<lEnd; ++ll) {

      /* Get the length scales of the adjacent elements. */
      const unsigned short llRel = ll - l;
      const unsigned long elemID0 = matchingInternalFaces[ll].elemID0;
      const unsigned long elemID1 = matchingInternalFaces[ll].elemID1;

      const su2double lenScale0 = volElem[elemID0].lenScale;
      const su2double lenScale1 = volElem[elemID1].lenScale;

      /* Call the function PenaltyTermsFluxFace to compute the actual penalty
         terms. Use the array viscFluxes as storage. */
      PenaltyTermsFluxFace(llRel, nInt, NPad, solIntL, solIntR, viscosityIntL,
                           viscosityIntR, kOverCvIntL, kOverCvIntR,
                           ConstPenFace, lenScale0, lenScale1,
                           matchingInternalFaces[ll].metricNormalsFace.data(),
                           viscFluxes);
    }

    /* Add the penalty fluxes to the earlier computed fluxes. */
    for(unsigned short j=0; j<(NPad*nInt); ++j) fluxes[j] += viscFluxes[j];

    /* Multiply the fluxes with the integration weight of the corresponding
       integration point. */
    for(unsigned short i=0; i<nInt; ++i) {
      su2double *flux = fluxes + i*NPad;

      for(unsigned short j=0; j<NPad; ++j)
        flux[j] *= weights[i];
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 4: Compute the contribution to the residuals from the       ---*/
    /*---         integration over this internal matching face.            ---*/
    /*------------------------------------------------------------------------*/

    /* Set the value of the offset of the residual between each of the fused
       faces. This depends whether or not symmetrizing terms are stored.
       Also set the location where the symmetrizing terms of the first face
       must be stored. */
    unsigned long offsetRes = nDOFsFace0 + nDOFsFace1;
    if( symmetrizingTermsPresent ) offsetRes += nDOFsElem0 + nDOFsElem1;

    unsigned long indResSym = indResFaces + nDOFsFace0 + nDOFsFace1;

    /* Get the correct form of the basis functions needed for the matrix
       multiplication to compute the residual. Use gradSolInt as a temporary
       buffer to store this product. */
    const su2double *basisFaceTrans = standardMatchingFacesSol[ind].GetBasisFaceIntegrationTransposeSide0();
    su2double *resSide0 = gradSolInt;

    /* Call the general function to carry out the matrix product. */
    blasFunctions->gemm(nDOFsFace0, NPad, nInt, basisFaceTrans, fluxes, resSide0, config);

    /* Check if the number of DOFs on both sides of the face is different.
       In that case also the matrix product with the basis functions on side 1
       must be carried out. Use viscFluxes as a temporary buffer to store
       this product. */
    su2double *resSide1 = viscFluxes;
    if(nDOFsFace1 != nDOFsFace0) {
      basisFaceTrans = standardMatchingFacesSol[ind].GetBasisFaceIntegrationTransposeSide1();
      blasFunctions->gemm(nDOFsFace1, NPad, nInt, basisFaceTrans, fluxes, resSide1, config);
    }

    /* Loop over the number of faces in this chunk. */
    for(unsigned short ll=0; ll<llEnd; ++ll) {
      const unsigned short llNVar = ll*nVar;

      /* Set the pointer in the residual array where the residual of side 0
         and side 1 of this face must be stored. Update the counter. */
      su2double *resFace0 = VecResFaces.data() + indResFaces*nVar;
      su2double *resFace1 = VecResFaces.data() + (indResFaces+nDOFsFace0)*nVar;
      indResFaces        += offsetRes;

      /* Loop over the DOFs and copy the data for side 0. */
      for(unsigned short i=0; i<nDOFsFace0; ++i)
        for(unsigned short mm=0; mm<nVar; ++mm)
          resFace0[nVar*i+mm] = resSide0[NPad*i+llNVar+mm];

      /* If the number of DOFs on both sides is the same, then the residual
         of side 1 is obtained by simply negating the residual of side 0.
         Otherwise a copy must be carried out and the residual negated,
         because the normal is pointing into the adjacent element for side 1. */
      if(nDOFsFace1 == nDOFsFace0) {
        for(unsigned short i=0; i<(nVar*nDOFsFace1); ++i)
          resFace1[i] = -resFace0[i];
      }
      else {
        for(unsigned short i=0; i<nDOFsFace1; ++i)
          for(unsigned short mm=0; mm<nVar; ++mm)
            resFace1[nVar*i+mm] = -resSide1[NPad*i+llNVar+mm];
      }
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 5: Compute and distribute the symmetrizing terms, if        ---*/
    /*---         present. Note that these terms must be distributed to    ---*/
    /*---         DOFs of the adjacent element, not only of the face.      ---*/
    /*------------------------------------------------------------------------*/

    if( symmetrizingTermsPresent ) {

      /* Use gradSolInt as a buffer to store the parametric fluxes. */
      su2double *paramFluxes = gradSolInt;

      /* The symmetrizing fluxes will be multiplied by 0.5 times the theta
         parameter. Store this factor a bit easier. */
      const su2double halfTheta = 0.5*config->GetTheta_Interior_Penalty_DGFEM();

      /* Loop over the faces in this chunk to compute the symmetrizing fluxes. */
      for(unsigned long ll=l; ll<lEnd; ++ll) {

        /* Compute the symmetrizing fluxes in the nDim directions in the
           integration points of the face. */
        const unsigned short llRel = ll - l;
        SymmetrizingFluxesFace(llRel, nInt, NPad, solIntL, solIntR, viscosityIntL,
                               viscosityIntR, kOverCvIntL, kOverCvIntR,
                               matchingInternalFaces[ll].metricNormalsFace.data(),
                               fluxes);

        /* Transform the fluxes, such that they must be multiplied with the
           gradients w.r.t. the parametric coordinates rather than the
           Cartesian coordinates of the basis functions. Also a multiplication
           with the integration weight and the theta parameter is carried out.
           In this loop only side 0 of the face is treated. */
        TransformSymmetrizingFluxes(llRel, nInt, NPad, halfTheta, fluxes, weights,
                                    matchingInternalFaces[ll].metricCoorDerivFace0.data(),
                                    paramFluxes);
      }

      /* Call the general function to carry out the matrix product to compute
         the residual for side 0. Use solIntL as storage for the residual. */
      const su2double *derBasisElemTrans = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationTransposeSide0();
      blasFunctions->gemm(nDOFsElem0, NPad, nInt*nDim, derBasisElemTrans, paramFluxes, solIntL, config);

      /* Loop over the faces in this chunk to compute the transformed
         symmetrizing fluxes for side 1. */
      for(unsigned long ll=l; ll<lEnd; ++ll) {
        const unsigned short llRel = ll - l;
        TransformSymmetrizingFluxes(llRel, nInt, NPad, halfTheta, fluxes, weights,
                                    matchingInternalFaces[ll].metricCoorDerivFace1.data(),
                                    paramFluxes);
      }

      /* Call the general function to carry out the matrix product to compute
         the residual for side 1. Note that the symmetrizing residual should not
         be negated, because two minus signs enter the formulation for side 1,
         which cancel each other. Use solIntR as storage for the residual. */
      derBasisElemTrans = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationTransposeSide1();
      blasFunctions->gemm(nDOFsElem1, NPad, nInt*nDim, derBasisElemTrans, paramFluxes, solIntR, config);

      /* Loop over the number of faces in this chunk. */
      for(unsigned short ll=0; ll<llEnd; ++ll) {
        const unsigned short llNVar = ll*nVar;

        /* Set the pointer in the residual array where the residual of side 0
           and side 1 of this face must be stored. Update the counter. */
        su2double *resElem0 = VecResFaces.data() + indResSym*nVar;
        su2double *resElem1 = VecResFaces.data() + (indResSym+nDOFsElem0)*nVar;
        indResSym          += offsetRes;

        /* Loop over the DOFs of the element and copy the data for side 0. */
        for(unsigned short i=0; i<nDOFsElem0; ++i)
          for(unsigned short mm=0; mm<nVar; ++mm)
            resElem0[nVar*i+mm] = solIntL[NPad*i+llNVar+mm];

        /* Loop over the DOFs of the element and copy the data for side 1. */
        for(unsigned short i=0; i<nDOFsElem1; ++i)
          for(unsigned short mm=0; mm<nVar; ++mm)
            resElem1[nVar*i+mm] = solIntR[NPad*i+llNVar+mm];
      }
    }

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::ViscousNormalFluxFace(const CVolumeElementFEM *adjVolElem,
                                             const unsigned short    indFaceChunk,
                                             const unsigned short    nInt,
                                             const unsigned short    NPad,
                                             const su2double         Wall_HeatFlux,
                                             const bool              HeatFlux_Prescribed,
                                             const su2double         *solInt,
                                             const su2double         *gradSolInt,
                                             const su2double         *metricCoorDerivFace,
                                             const su2double         *metricNormalsFace,
                                             const su2double         *wallDistanceInt,
                                                   su2double         *viscNormFluxes,
                                                   su2double         *viscosityInt,
                                                   su2double         *kOverCvInt) {

  /* Multiplication factor for the heat flux. It is set to zero if the wall heat flux
     is prescribed, such that the computed heat flux is zero, and to one otherwise. */
  const su2double factHeatFlux = HeatFlux_Prescribed ? su2double(0.0): su2double(1.0);

  /* Set the value of the prescribed heat flux for the same reason. */
  const su2double HeatFlux = HeatFlux_Prescribed ? Wall_HeatFlux : su2double(0.0);

  /* Compute the length scale for the LES of the adjacent element. */
  const unsigned short iind  = adjVolElem->indStandardElement;
  unsigned short       nPoly = standardElementsSol[iind].GetNPoly();
  if(nPoly == 0) nPoly = 1;

  const su2double lenScale_LES = adjVolElem->lenScale/nPoly;

  /* Determine the offset between r- and -s-derivatives, which is also the
     offset between s- and t-derivatives. */
  const unsigned short offDeriv = NPad*nInt;

  /* Make a distinction between two and three space dimensions
     in order to have the most efficient code. */
  switch( nDim ) {

    case 2: {

      /* 2D simulation. Loop over the integration points to
         compute the viscous fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the metric terms needed to compute the Cartesian
           gradients in this integration point and the starting locations of
           the solution and the gradients, w.r.t. the parametric coordinates
           of this solution. */
        const unsigned short offPointer = NPad*i + nVar*indFaceChunk;
        const su2double *metricTerms = metricCoorDerivFace + i*nDim*nDim;
        const su2double *sol         = solInt     + offPointer;
        const su2double *dSolDr      = gradSolInt + offPointer;
        const su2double *dSolDs      = dSolDr     + offDeriv;

        /* Easier storage of the metric terms in this integration point. */
        const su2double drdx = metricTerms[0];
        const su2double drdy = metricTerms[1];

        const su2double dsdx = metricTerms[2];
        const su2double dsdy = metricTerms[3];

        /*--- Compute the Cartesian gradients of the solution. ---*/
        su2double solGradCart[4][2];

        solGradCart[0][0] = dSolDr[0]*drdx + dSolDs[0]*dsdx;
        solGradCart[1][0] = dSolDr[1]*drdx + dSolDs[1]*dsdx;
        solGradCart[2][0] = dSolDr[2]*drdx + dSolDs[2]*dsdx;
        solGradCart[3][0] = dSolDr[3]*drdx + dSolDs[3]*dsdx;

        solGradCart[0][1] = dSolDr[0]*drdy + dSolDs[0]*dsdy;
        solGradCart[1][1] = dSolDr[1]*drdy + dSolDs[1]*dsdy;
        solGradCart[2][1] = dSolDr[2]*drdy + dSolDs[2]*dsdy;
        solGradCart[3][1] = dSolDr[3]*drdy + dSolDs[3]*dsdy;

        /*--- Call the function ViscousNormalFluxIntegrationPoint to compute the
              actual normal viscous flux. The viscosity and thermal conductivity
              are stored for later use. ---*/
        const su2double *normal  = metricNormalsFace + i*(nDim+1);
        su2double *normalFlux    = viscNormFluxes + offPointer;
        const su2double wallDist = wallDistanceInt ? wallDistanceInt[i] : 0.0;

        su2double Viscosity, kOverCv;

        ViscousNormalFluxIntegrationPoint_2D(sol, solGradCart, normal, HeatFlux,
                                             factHeatFlux, wallDist, lenScale_LES,
                                             Viscosity, kOverCv, normalFlux);

        const unsigned short ind = indFaceChunk*nInt + i;
        viscosityInt[ind] = Viscosity;
        kOverCvInt[ind]   = kOverCv;
      }

      break;
    }

    /*------------------------------------------------------------------------*/

    case 3: {

      /* 3D simulation. Loop over the integration points to
         compute the viscous fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the metric terms needed to compute the Cartesian
           gradients in this integration point and the starting locations of
           the solution and the gradients, w.r.t. the parametric coordinates
           of this solution. */
        const unsigned short offPointer = NPad*i + nVar*indFaceChunk;
        const su2double *metricTerms = metricCoorDerivFace + i*nDim*nDim;
        const su2double *sol         = solInt     + offPointer;
        const su2double *dSolDr      = gradSolInt + offPointer;
        const su2double *dSolDs      = dSolDr     + offDeriv;
        const su2double *dSolDt      = dSolDs     + offDeriv;

        /* Easier storage of the metric terms in this integration point. */
        const su2double drdx = metricTerms[0];
        const su2double drdy = metricTerms[1];
        const su2double drdz = metricTerms[2];

        const su2double dsdx = metricTerms[3];
        const su2double dsdy = metricTerms[4];
        const su2double dsdz = metricTerms[5];

        const su2double dtdx = metricTerms[6];
        const su2double dtdy = metricTerms[7];
        const su2double dtdz = metricTerms[8];

        /*--- Compute the Cartesian gradients of the solution. ---*/
        su2double solGradCart[5][3];

        solGradCart[0][0] = dSolDr[0]*drdx + dSolDs[0]*dsdx + dSolDt[0]*dtdx;
        solGradCart[1][0] = dSolDr[1]*drdx + dSolDs[1]*dsdx + dSolDt[1]*dtdx;
        solGradCart[2][0] = dSolDr[2]*drdx + dSolDs[2]*dsdx + dSolDt[2]*dtdx;
        solGradCart[3][0] = dSolDr[3]*drdx + dSolDs[3]*dsdx + dSolDt[3]*dtdx;
        solGradCart[4][0] = dSolDr[4]*drdx + dSolDs[4]*dsdx + dSolDt[4]*dtdx;

        solGradCart[0][1] = dSolDr[0]*drdy + dSolDs[0]*dsdy + dSolDt[0]*dtdy;
        solGradCart[1][1] = dSolDr[1]*drdy + dSolDs[1]*dsdy + dSolDt[1]*dtdy;
        solGradCart[2][1] = dSolDr[2]*drdy + dSolDs[2]*dsdy + dSolDt[2]*dtdy;
        solGradCart[3][1] = dSolDr[3]*drdy + dSolDs[3]*dsdy + dSolDt[3]*dtdy;
        solGradCart[4][1] = dSolDr[4]*drdy + dSolDs[4]*dsdy + dSolDt[4]*dtdy;

        solGradCart[0][2] = dSolDr[0]*drdz + dSolDs[0]*dsdz + dSolDt[0]*dtdz;
        solGradCart[1][2] = dSolDr[1]*drdz + dSolDs[1]*dsdz + dSolDt[1]*dtdz;
        solGradCart[2][2] = dSolDr[2]*drdz + dSolDs[2]*dsdz + dSolDt[2]*dtdz;
        solGradCart[3][2] = dSolDr[3]*drdz + dSolDs[3]*dsdz + dSolDt[3]*dtdz;
        solGradCart[4][2] = dSolDr[4]*drdz + dSolDs[4]*dsdz + dSolDt[4]*dtdz;

        /*--- Call the function ViscousNormalFluxIntegrationPoint to compute the
              actual normal viscous flux. The viscosity and thermal conductivity
              are stored for later use. ---*/
        const su2double *normal  = metricNormalsFace + i*(nDim+1);
        su2double *normalFlux    = viscNormFluxes + offPointer;
        const su2double wallDist = wallDistanceInt ? wallDistanceInt[i] : 0.0;

        su2double Viscosity, kOverCv;

        ViscousNormalFluxIntegrationPoint_3D(sol, solGradCart, normal, HeatFlux,
                                             factHeatFlux, wallDist, lenScale_LES,
                                             Viscosity, kOverCv, normalFlux);

        const unsigned short ind = indFaceChunk*nInt + i;
        viscosityInt[ind] = Viscosity;
        kOverCvInt[ind]   = kOverCv;
      }

      break;
    }
  }
}

void CFEM_DG_NSSolver::ViscousNormalFluxIntegrationPoint_2D(const su2double *sol,
                                                            const su2double solGradCart[4][2],
                                                            const su2double *normal,
                                                            const su2double HeatFlux,
                                                            const su2double factHeatFlux,
                                                            const su2double wallDist,
                                                            const su2double lenScale_LES,
                                                                  su2double &Viscosity,
                                                                  su2double &kOverCv,
                                                                  su2double *normalFlux) {

  /* Constant factor present in the heat flux vector, namely the ratio of
     thermal conductivity and viscosity. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /*--- Compute the velocities and static energy in this integration point. ---*/
  const su2double rhoInv = 1.0/sol[0];
  const su2double u = rhoInv*sol[1];
  const su2double v = rhoInv*sol[2];

  const su2double TotalEnergy  = rhoInv*sol[3];
  const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v);

  /*--- Compute the Cartesian gradients of the velocities and static energy
        in this integration point and also the divergence of the velocity. ---*/
  const su2double dudx = rhoInv*(solGradCart[1][0] - u*solGradCart[0][0]);
  const su2double dudy = rhoInv*(solGradCart[1][1] - u*solGradCart[0][1]);

  const su2double dvdx = rhoInv*(solGradCart[2][0] - v*solGradCart[0][0]);
  const su2double dvdy = rhoInv*(solGradCart[2][1] - v*solGradCart[0][1]);

  const su2double dStaticEnergyDx = rhoInv*(solGradCart[3][0]
                                  -         TotalEnergy*solGradCart[0][0])
                                  - u*dudx - v*dvdx;
  const su2double dStaticEnergyDy = rhoInv*(solGradCart[3][1]
                                  -         TotalEnergy*solGradCart[0][1])
                                  - u*dudy - v*dvdy;

  const su2double divVel = dudx + dvdy;

  /*--- Compute the laminar viscosity. ---*/
  FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
  const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

  /*--- Compute the eddy viscosity, if needed. ---*/
  su2double ViscosityTurb = 0.0;
  if( SGSModelUsed )
    ViscosityTurb = SGSModel->ComputeEddyViscosity_2D(sol[0], dudx, dudy, dvdx,
                                                      dvdy, lenScale_LES, wallDist);

  /* Compute the total viscosity and heat conductivity. Note that the heat
     conductivity is divided by the Cv, because gradients of internal energy
     are computed and not temperature. */
  Viscosity = ViscosityLam + ViscosityTurb;
  kOverCv   = ViscosityLam*factHeatFlux_Lam + ViscosityTurb*factHeatFlux_Turb;

  /*--- Set the value of the second viscosity and compute the divergence
        term in the viscous normal stresses. ---*/
  const su2double lambda     = -TWO3*Viscosity;
  const su2double lamDivTerm =  lambda*divVel;

  /*--- Compute the viscous stress tensor and minus the heatflux vector.
        The heat flux vector is multiplied by factHeatFlux, such that the
        case of a prescribed heat flux is treated correctly. ---*/
  const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
  const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
  const su2double tauxy = Viscosity*(dudy + dvdx);

  const su2double qx = factHeatFlux*kOverCv*dStaticEnergyDx;
  const su2double qy = factHeatFlux*kOverCv*dStaticEnergyDy;

  /* Compute the unscaled normal vector. */
  const su2double nx = normal[0]*normal[2];
  const su2double ny = normal[1]*normal[2];

  /*--- Compute the viscous normal flux. Note that the energy flux get a
        contribution from both the prescribed and the computed heat flux.
        At least one of these terms is zero. ---*/
  normalFlux[0] = 0.0;
  normalFlux[1] = tauxx*nx + tauxy*ny;
  normalFlux[2] = tauxy*nx + tauyy*ny;
  normalFlux[3] = normal[2]*HeatFlux
                + (u*tauxx + v*tauxy + qx)*nx + (u*tauxy + v*tauyy + qy)*ny;
}

void CFEM_DG_NSSolver::ViscousNormalFluxIntegrationPoint_3D(const su2double *sol,
                                                            const su2double solGradCart[5][3],
                                                            const su2double *normal,
                                                            const su2double HeatFlux,
                                                            const su2double factHeatFlux,
                                                            const su2double wallDist,
                                                            const su2double lenScale_LES,
                                                                  su2double &Viscosity,
                                                                  su2double &kOverCv,
                                                                  su2double *normalFlux) {

  /* Constant factor present in the heat flux vector, namely the ratio of
     thermal conductivity and viscosity. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /*--- Compute the velocities and static energy in this integration point. ---*/
  const su2double rhoInv = 1.0/sol[0];
  const su2double u = rhoInv*sol[1];
  const su2double v = rhoInv*sol[2];
  const su2double w = rhoInv*sol[3];

  const su2double TotalEnergy  = rhoInv*sol[4];
  const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v + w*w);

  /*--- Compute the Cartesian gradients of the velocities and static energy
        in this integration point and also the divergence of the velocity. ---*/
  const su2double dudx = rhoInv*(solGradCart[1][0] - u*solGradCart[0][0]);
  const su2double dudy = rhoInv*(solGradCart[1][1] - u*solGradCart[0][1]);
  const su2double dudz = rhoInv*(solGradCart[1][2] - u*solGradCart[0][2]);

  const su2double dvdx = rhoInv*(solGradCart[2][0] - v*solGradCart[0][0]);
  const su2double dvdy = rhoInv*(solGradCart[2][1] - v*solGradCart[0][1]);
  const su2double dvdz = rhoInv*(solGradCart[2][2] - v*solGradCart[0][2]);

  const su2double dwdx = rhoInv*(solGradCart[3][0] - w*solGradCart[0][0]);
  const su2double dwdy = rhoInv*(solGradCart[3][1] - w*solGradCart[0][1]);
  const su2double dwdz = rhoInv*(solGradCart[3][2] - w*solGradCart[0][2]);

  const su2double dStaticEnergyDx = rhoInv*(solGradCart[4][0]
                                  -         TotalEnergy*solGradCart[0][0])
                                  - u*dudx - v*dvdx - w*dwdx;
  const su2double dStaticEnergyDy = rhoInv*(solGradCart[4][1]
                                  -         TotalEnergy*solGradCart[0][1])
                                  - u*dudy - v*dvdy - w*dwdy;
  const su2double dStaticEnergyDz = rhoInv*(solGradCart[4][2]
                                  -         TotalEnergy*solGradCart[0][2])
                                  - u*dudz - v*dvdz - w*dwdz;

  const su2double divVel = dudx + dvdy + dwdz;

  /*--- Compute the laminar viscosity. ---*/
  FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
  const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

  /*--- Compute the eddy viscosity, if needed. ---*/
  su2double ViscosityTurb = 0.0;
  if( SGSModelUsed )
    ViscosityTurb = SGSModel->ComputeEddyViscosity_3D(sol[0], dudx, dudy, dudz,
                                                      dvdx, dvdy, dvdz, dwdx,
                                                      dwdy, dwdz, lenScale_LES,
                                                      wallDist);

  /* Compute the total viscosity and heat conductivity. Note that the heat
     conductivity is divided by the Cv, because gradients of internal energy
     are computed and not temperature. */
  Viscosity = ViscosityLam + ViscosityTurb;
  kOverCv   = ViscosityLam*factHeatFlux_Lam + ViscosityTurb*factHeatFlux_Turb;

  /*--- Set the value of the second viscosity and compute the divergence
        term in the viscous normal stresses. ---*/
  const su2double lambda     = -TWO3*Viscosity;
  const su2double lamDivTerm =  lambda*divVel;

  /*--- Compute the viscous stress tensor and minus the heatflux vector.
        The heat flux vector is multiplied by factHeatFlux, such that the
        case of a prescribed heat flux is treated correctly. ---*/
  const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
  const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
  const su2double tauzz = 2.0*Viscosity*dwdz + lamDivTerm;

  const su2double tauxy = Viscosity*(dudy + dvdx);
  const su2double tauxz = Viscosity*(dudz + dwdx);
  const su2double tauyz = Viscosity*(dvdz + dwdy);

  const su2double qx = factHeatFlux*kOverCv*dStaticEnergyDx;
  const su2double qy = factHeatFlux*kOverCv*dStaticEnergyDy;
  const su2double qz = factHeatFlux*kOverCv*dStaticEnergyDz;

  /* Compute the unscaled normal vector. */
  const su2double nx = normal[0]*normal[3];
  const su2double ny = normal[1]*normal[3];
  const su2double nz = normal[2]*normal[3];

  /*--- Compute the viscous normal flux. Note that the energy flux get a
        contribution from both the prescribed and the computed heat flux.
        At least one of these terms is zero. ---*/
  normalFlux[0] = 0.0;
  normalFlux[1] = tauxx*nx + tauxy*ny + tauxz*nz;
  normalFlux[2] = tauxy*nx + tauyy*ny + tauyz*nz;
  normalFlux[3] = tauxz*nx + tauyz*ny + tauzz*nz;
  normalFlux[4] = normal[3]*HeatFlux
                + (u*tauxx + v*tauxy + w*tauxz + qx)*nx
                + (u*tauxy + v*tauyy + w*tauyz + qy)*ny
                + (u*tauxz + v*tauyz + w*tauzz + qz)*nz;
}

void CFEM_DG_NSSolver::PenaltyTermsFluxFace(const unsigned short indFaceChunk,
                                            const unsigned short nInt,
                                            const unsigned short NPad,
                                            const su2double      *solInt0,
                                            const su2double      *solInt1,
                                            const su2double      *viscosityInt0,
                                            const su2double      *viscosityInt1,
                                            const su2double      *kOverCvInt0,
                                            const su2double      *kOverCvInt1,
                                            const su2double      ConstPenFace,
                                            const su2double      lenScale0,
                                            const su2double      lenScale1,
                                            const su2double      *metricNormalsFace,
                                                  su2double      *penaltyFluxes) {

  /* Constant ratio of the second viscosity and the viscosity itself. */
  const su2double lambdaOverMu = -TWO3;

  /* The eigenvalues of the viscous Jacobian, scaled by the kinematic viscosity,
     are 1.0, 2.0 + lambdaOverMu and kOverCv/Mu. The last is variable due to the
     possible presence of an eddy viscosity, but the first two are constant and
     the maximum can be determined. */
  const su2double radOverNuTerm = max(1.0, 2.0+lambdaOverMu);

  /*--- Make a distinction between 2D and 3D for efficiency. ---*/
  switch ( nDim ) {
    case 2: {

      /* 2D simulation. Loop over the integration points to compute
         the penalty fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the variables for this integration point. */
        const unsigned short offPointer = NPad*i + nVar*indFaceChunk;
        const su2double *sol0   = solInt0 + offPointer;
        const su2double *sol1   = solInt1 + offPointer;
        const su2double *normal = metricNormalsFace + i*(nDim+1);
        su2double       *flux   = penaltyFluxes + offPointer;

        /* Determine the ratio of kOverCv and mu for both sides and compute the
           spectral radius of the viscous terms, scaled by kinematic viscosity. */
        const unsigned short ind = indFaceChunk*nInt + i;
        const su2double factHeatFlux0 = kOverCvInt0[ind]/viscosityInt0[ind];
        const su2double factHeatFlux1 = kOverCvInt1[ind]/viscosityInt1[ind];

        const su2double radOverNu0 = max(radOverNuTerm, factHeatFlux0);
        const su2double radOverNu1 = max(radOverNuTerm, factHeatFlux1);

        /* Compute the kinematic viscosities of both sides. Multiply it by
           ConstPenFace and divide by the length scale. */
        const su2double nu0 = ConstPenFace*viscosityInt0[ind]/(lenScale0*sol0[0]);
        const su2double nu1 = ConstPenFace*viscosityInt1[ind]/(lenScale1*sol1[0]);

        /* Compute the penalty parameter of this face as the maximum of the
           penalty parameter from both sides. Multiply by the area to obtain the
           correct expression. */
        const su2double pen0 = radOverNu0*nu0;
        const su2double pen1 = radOverNu1*nu1;

        const su2double penFace = normal[2]*max(pen0, pen1);

        /* Compute the penalty flux, where it is assumed that the normal points from
           side 0 to side 1. */
        flux[0] = penFace*(sol0[0] - sol1[0]);
        flux[1] = penFace*(sol0[1] - sol1[1]);
        flux[2] = penFace*(sol0[2] - sol1[2]);
        flux[3] = penFace*(sol0[3] - sol1[3]);
      }

      break;
    }

    /*----------------------------------------------------------------------*/

    case 3: {

      /* 3D simulation. Loop over the integration points to compute
         the penalty fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the variables for this integration point. */
        const unsigned short offPointer = NPad*i + nVar*indFaceChunk;
        const su2double *sol0   = solInt0 + offPointer;
        const su2double *sol1   = solInt1 + offPointer;
        const su2double *normal = metricNormalsFace + i*(nDim+1);
        su2double       *flux   = penaltyFluxes + offPointer;

        /* Determine the ratio of kOverCv and mu for both sides and compute the
           spectral radius of the viscous terms, scaled by kinematic viscosity. */
        const unsigned short ind = indFaceChunk*nInt + i;
        const su2double factHeatFlux0 = kOverCvInt0[ind]/viscosityInt0[ind];
        const su2double factHeatFlux1 = kOverCvInt1[ind]/viscosityInt1[ind];

        const su2double radOverNu0 = max(radOverNuTerm, factHeatFlux0);
        const su2double radOverNu1 = max(radOverNuTerm, factHeatFlux1);

        /* Compute the kinematic viscosities of both sides. Multiply it by
           ConstPenFace and divide by the length scale. */
        const su2double nu0 = ConstPenFace*viscosityInt0[ind]/(lenScale0*sol0[0]);
        const su2double nu1 = ConstPenFace*viscosityInt1[ind]/(lenScale1*sol1[0]);

        /* Compute the penalty parameter of this face as the maximum of the
           penalty parameter from both sides. Multiply by the area to obtain the
           correct expression. */
        const su2double pen0 = radOverNu0*nu0;
        const su2double pen1 = radOverNu1*nu1;

        const su2double penFace = normal[3]*max(pen0, pen1);

        /* Compute the penalty flux, where it is assumed that the normal points from
           side 0 to side 1. */
        flux[0] = penFace*(sol0[0] - sol1[0]);
        flux[1] = penFace*(sol0[1] - sol1[1]);
        flux[2] = penFace*(sol0[2] - sol1[2]);
        flux[3] = penFace*(sol0[3] - sol1[3]);
        flux[4] = penFace*(sol0[4] - sol1[4]);
      }

      break;
    }
  }
}

void CFEM_DG_NSSolver::SymmetrizingFluxesFace(const unsigned short indFaceChunk,
                                              const unsigned short nInt,
                                              const unsigned short NPad,
                                              const su2double      *solInt0,
                                              const su2double      *solInt1,
                                              const su2double      *viscosityInt0,
                                              const su2double      *viscosityInt1,
                                              const su2double      *kOverCvInt0,
                                              const su2double      *kOverCvInt1,
                                              const su2double      *metricNormalsFace,
                                                    su2double      *symmFluxes) {

  /* Constant ratio of the second viscosity and the viscosity itself. */
  const su2double lambdaOverMu = -TWO3;

  /*--- Set two factors such that either the original or the transposed diffusion
        tensor is taken in the symmetrizing fluxes. ---*/

  /* Use the following line for the original formulation. */
  const su2double alpha = lambdaOverMu;

  /* Use the following line for the transposed formulation. */
  //const su2double alpha = 1.0;

  /* Other constants, which appear in the symmetrizing fluxes. */
  const su2double beta     = lambdaOverMu + 1.0 - alpha;
  const su2double alphaP1  = alpha + 1.0;
  const su2double lambdaP1 = lambdaOverMu + 1.0;

  /* Make a distinction between two and three space dimensions
     in order to have the most efficient code. */
  switch( nDim ) {

    case 2: {

      /*--- 2D simulation. Loop over the number of integration points
            of the face. ---*/
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the variables for this integration point. */
        const unsigned short offPointer = NPad*i + nVar*indFaceChunk;
        const su2double *sol0   = solInt0 + offPointer;
        const su2double *sol1   = solInt1 + offPointer;
        const su2double *normal = metricNormalsFace + i*(nDim+1);

        su2double *flux = symmFluxes + i*nDim*NPad + nVar*indFaceChunk;

         /* Determine the difference in conservative variables. Multiply these
           differences by the length of the normal vector to obtain the correct
           dimensions for the symmetrizing fluxes. */
        const su2double dSol[] = {normal[2]*(sol0[0] - sol1[0]),
                                  normal[2]*(sol0[1] - sol1[1]),
                                  normal[2]*(sol0[2] - sol1[2]),
                                  normal[2]*(sol0[3] - sol1[3])};

        /*--- Compute the terms that occur in the symmetrizing fluxes
              for state 0 and state 1. ---*/
        const unsigned short ind = indFaceChunk*nInt + i;

        const su2double DensityInv0 = 1.0/sol0[0],             DensityInv1 = 1.0/sol1[0];
        const su2double Etot0 = DensityInv0*sol0[3],           Etot1 = DensityInv1*sol1[3];
        const su2double nu0 = DensityInv0*viscosityInt0[ind],  nu1 = DensityInv1*viscosityInt1[ind];
        const su2double kScal0 = DensityInv0*kOverCvInt0[ind], kScal1 = DensityInv1*kOverCvInt1[ind];

        const su2double vel0[] = {DensityInv0*sol0[1], DensityInv0*sol0[2]};
        const su2double vel1[] = {DensityInv1*sol1[1], DensityInv1*sol1[2]};

        const su2double velNorm0 = vel0[0]*normal[0] + vel0[1]*normal[1];
        const su2double velNorm1 = vel1[0]*normal[0] + vel1[1]*normal[1];

        const su2double velSquared0 = vel0[0]*vel0[0] + vel0[1]*vel0[1];
        const su2double velSquared1 = vel1[0]*vel1[0] + vel1[1]*vel1[1];

        /*--- Compute the average of the terms that occur in the symmetrizing
              fluxes. The average of the left and right terms is taken, rather
              than the terms evaluated at the average state, because the viscous
              fluxes are also computed as the average of the fluxes and not the
              fluxes of the averaged state. ---*/
        const su2double nuAvg    = 0.5*(nu0    + nu1);
        const su2double kScalAvg = 0.5*(kScal0 + kScal1);
        const su2double nuVelSquaredAvg     = 0.5*(nu0*velSquared0 + nu1*velSquared1);
        const su2double nuVelNormAve        = 0.5*(nu0*velNorm0    + nu1*velNorm1);
        const su2double kScalEminVelSquaredAve = 0.5*(kScal0*(Etot0-velSquared0)
                                               +      kScal1*(Etot1-velSquared1));

        const su2double nuVelAvg[] = {0.5*(nu0*vel0[0] + nu1*vel1[0]),
                                      0.5*(nu0*vel0[1] + nu1*vel1[1])};
        const su2double kScalVelAvg[] = {0.5*(kScal0*vel0[0] + kScal1*vel1[0]),
                                         0.5*(kScal0*vel0[1] + kScal1*vel1[1])};
        const su2double nuVelVelAvg[] = {0.5*(nu0*vel0[0]*velNorm0 + nu1*vel1[0]*velNorm1),
                                         0.5*(nu0*vel0[1]*velNorm0 + nu1*vel1[1]*velNorm1)};

        /*--- Abbreviations to make the flux computations a bit more efficient. ---*/
        const su2double abv1 = normal[0]  *dSol[1] + normal[1]  *dSol[2];
        const su2double abv2 = nuVelAvg[0]*dSol[1] + nuVelAvg[1]*dSol[2];

        const su2double abv2kScal = kScalVelAvg[0]*dSol[1] + kScalVelAvg[1]*dSol[2];

        const su2double abv3 = beta*(nuAvg*abv1 - nuVelNormAve*dSol[0]);
        const su2double abv4 = kScalAvg*dSol[3] - abv2kScal
                             - kScalEminVelSquaredAve*dSol[0] + abv2;

        /*--- Compute the symmetrizing fluxes in x-direction. ---*/
        flux[0] = 0.0;
        flux[1] = abv3 + alphaP1*normal[0]*(nuAvg*dSol[1] - nuVelAvg[0]*dSol[0]);
        flux[2] = nuAvg*(normal[0]*dSol[2] + alpha*normal[1]*dSol[1])
                - (normal[0]*nuVelAvg[1] + alpha*normal[1]*nuVelAvg[0])*dSol[0];
        flux[3] = normal[0]*abv4
                - (lambdaP1*nuVelVelAvg[0] + nuVelSquaredAvg*normal[0])*dSol[0]
                + alpha*nuVelNormAve*dSol[1] + beta*nuVelAvg[0]*abv1;

        /*--- Compute the symmetrizing fluxes in y-direction. */
        flux = flux + NPad;
        flux[0] = 0.0;
        flux[1] = nuAvg*(normal[1]*dSol[1] + alpha*normal[0]*dSol[2])
                - (normal[1]*nuVelAvg[0] + alpha*normal[0]*nuVelAvg[1])*dSol[0];
        flux[2] = abv3 + alphaP1*normal[1]*(nuAvg*dSol[2] - nuVelAvg[1]*dSol[0]);
        flux[3] = normal[1]*abv4
                - (lambdaP1*nuVelVelAvg[1] + nuVelSquaredAvg*normal[1])*dSol[0]
                + alpha*nuVelNormAve*dSol[2] + beta*nuVelAvg[1]*abv1;
      }

      break;
    }

    /*------------------------------------------------------------------------*/

    case 3: {

      /*--- 3D simulation. Loop over the number of integration points
            of the face. ---*/
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the variables for this integration point. */
        const unsigned short offPointer = NPad*i + nVar*indFaceChunk;
        const su2double *sol0   = solInt0 + offPointer;
        const su2double *sol1   = solInt1 + offPointer;
        const su2double *normal = metricNormalsFace + i*(nDim+1);

        su2double *flux = symmFluxes + i*nDim*NPad + nVar*indFaceChunk;

        /* Determine the difference in conservative variables. Multiply these
           differences by the length of the normal vector to obtain the correct
           dimensions for the symmetrizing fluxes. */
        const su2double dSol[] = {normal[3]*(sol0[0] - sol1[0]),
                                  normal[3]*(sol0[1] - sol1[1]),
                                  normal[3]*(sol0[2] - sol1[2]),
                                  normal[3]*(sol0[3] - sol1[3]),
                                  normal[3]*(sol0[4] - sol1[4])};

        /*--- Compute the terms that occur in the symmetrizing fluxes
              for state 0 and state 1. ---*/
        const unsigned short ind = indFaceChunk*nInt + i;

        const su2double DensityInv0 = 1.0/sol0[0],             DensityInv1 = 1.0/sol1[0];
        const su2double Etot0 = DensityInv0*sol0[4],           Etot1 = DensityInv1*sol1[4];
        const su2double nu0 = DensityInv0*viscosityInt0[ind],  nu1 = DensityInv1*viscosityInt1[ind];
        const su2double kScal0 = DensityInv0*kOverCvInt0[ind], kScal1 = DensityInv1*kOverCvInt1[ind];

        const su2double vel0[] = {DensityInv0*sol0[1], DensityInv0*sol0[2], DensityInv0*sol0[3]};
        const su2double vel1[] = {DensityInv1*sol1[1], DensityInv1*sol1[2], DensityInv1*sol1[3]};

        const su2double velNorm0 = vel0[0]*normal[0] + vel0[1]*normal[1] + vel0[2]*normal[2];
        const su2double velNorm1 = vel1[0]*normal[0] + vel1[1]*normal[1] + vel1[2]*normal[2];

        const su2double velSquared0 = vel0[0]*vel0[0] + vel0[1]*vel0[1] + vel0[2]*vel0[2];
        const su2double velSquared1 = vel1[0]*vel1[0] + vel1[1]*vel1[1] + vel1[2]*vel1[2];

        /*--- Compute the average of the terms that occur in the symmetrizing
              fluxes. The average of the left and right terms is taken, rather
              than the terms evaluated at the average state, because the viscous
              fluxes are also computed as the average of the fluxes and not the
              fluxes of the averaged state. ---*/
        const su2double nuAvg    = 0.5*(nu0    + nu1);
        const su2double kScalAvg = 0.5*(kScal0 + kScal1);
        const su2double nuVelSquaredAvg     = 0.5*(nu0*velSquared0 + nu1*velSquared1);
        const su2double nuVelNormAve        = 0.5*(nu0*velNorm0    + nu1*velNorm1);
        const su2double kScalEminVelSquaredAve = 0.5*(kScal0*(Etot0-velSquared0)
                                               +      kScal1*(Etot1-velSquared1));

        const su2double nuVelAvg[] = {0.5*(nu0*vel0[0] + nu1*vel1[0]),
                                      0.5*(nu0*vel0[1] + nu1*vel1[1]),
                                      0.5*(nu0*vel0[2] + nu1*vel1[2])};
        const su2double kScalVelAvg[] = {0.5*(kScal0*vel0[0] + kScal1*vel1[0]),
                                         0.5*(kScal0*vel0[1] + kScal1*vel1[1]),
                                         0.5*(kScal0*vel0[2] + kScal1*vel1[2])};
        const su2double nuVelVelAvg[] = {0.5*(nu0*vel0[0]*velNorm0 + nu1*vel1[0]*velNorm1),
                                         0.5*(nu0*vel0[1]*velNorm0 + nu1*vel1[1]*velNorm1),
                                         0.5*(nu0*vel0[2]*velNorm0 + nu1*vel1[2]*velNorm1)};

        /*--- Abbreviations to make the flux computations a bit more efficient. ---*/
        const su2double abv1 = normal[0]  *dSol[1] + normal[1]  *dSol[2] + normal[2]  *dSol[3];
        const su2double abv2 = nuVelAvg[0]*dSol[1] + nuVelAvg[1]*dSol[2] + nuVelAvg[2]*dSol[3];

        const su2double abv2kScal = kScalVelAvg[0]*dSol[1] + kScalVelAvg[1]*dSol[2]
                                  + kScalVelAvg[2]*dSol[3];

        const su2double abv3 = beta*(nuAvg*abv1 - nuVelNormAve*dSol[0]);
        const su2double abv4 = kScalAvg*dSol[4] - abv2kScal
                             - kScalEminVelSquaredAve*dSol[0] + abv2;

        /*--- Compute the symmetrizing fluxes in x-direction. ---*/
        flux[0] = 0.0;
        flux[1] = abv3 + alphaP1*normal[0]*(nuAvg*dSol[1] - nuVelAvg[0]*dSol[0]);
        flux[2] = nuAvg*(normal[0]*dSol[2] + alpha*normal[1]*dSol[1])
                - (normal[0]*nuVelAvg[1] + alpha*normal[1]*nuVelAvg[0])*dSol[0];
        flux[3] = nuAvg*(normal[0]*dSol[3] + alpha*normal[2]*dSol[1])
                - (normal[0]*nuVelAvg[2] + alpha*normal[2]*nuVelAvg[0])*dSol[0];
        flux[4] = normal[0]*abv4
                - (lambdaP1*nuVelVelAvg[0] + nuVelSquaredAvg*normal[0])*dSol[0]
                + alpha*nuVelNormAve*dSol[1] + beta*nuVelAvg[0]*abv1;

        /*--- Compute the symmetrizing fluxes in y-direction. ---*/
        flux    = flux + NPad;
        flux[0] = 0.0;
        flux[1] = nuAvg*(normal[1]*dSol[1] + alpha*normal[0]*dSol[2])
                - (normal[1]*nuVelAvg[0] + alpha*normal[0]*nuVelAvg[1])*dSol[0];
        flux[2] = abv3 + alphaP1*normal[1]*(nuAvg*dSol[2] - nuVelAvg[1]*dSol[0]);
        flux[3] = nuAvg*(normal[1]*dSol[3] + alpha*normal[2]*dSol[2])
                - (normal[1]*nuVelAvg[2] + alpha*normal[2]*nuVelAvg[1])*dSol[0];
        flux[4] = normal[1]*abv4
                - (lambdaP1*nuVelVelAvg[1] + nuVelSquaredAvg*normal[1])*dSol[0]
                + alpha*nuVelNormAve*dSol[2] + beta*nuVelAvg[1]*abv1;

        /*--- Compute the symmetrizing fluxes in z-direction. ---*/
        flux    = flux + NPad;
        flux[0] = 0.0;
        flux[1] = nuAvg*(normal[2]*dSol[1] + alpha*normal[0]*dSol[3])
                - (normal[2]*nuVelAvg[0] + alpha*normal[0]*nuVelAvg[2])*dSol[0];
        flux[2] = nuAvg*(normal[2]*dSol[2] + alpha*normal[1]*dSol[3])
                - (normal[2]*nuVelAvg[1] + alpha*normal[1]*nuVelAvg[2])*dSol[0];
        flux[3] = abv3 + alphaP1*normal[2]*(nuAvg*dSol[3] - nuVelAvg[2]*dSol[0]);
        flux[4] = normal[2]*abv4
                - (lambdaP1*nuVelVelAvg[2] + nuVelSquaredAvg*normal[2])*dSol[0]
                + alpha*nuVelNormAve*dSol[3] + beta*nuVelAvg[2]*abv1;
      }

      break;
    }
  }
}

void CFEM_DG_NSSolver::TransformSymmetrizingFluxes(const unsigned short indFaceChunk,
                                                   const unsigned short nInt,
                                                   const unsigned short NPad,
                                                   const su2double      halfTheta,
                                                   const su2double      *symmFluxes,
                                                   const su2double      *weights,
                                                   const su2double      *metricCoorFace,
                                                         su2double      *paramFluxes) {

  /*--- Transform the fluxes, such that they must be multiplied with the
        gradients w.r.t. the parametric coordinates rather than the
        Cartesian coordinates of the basis functions. This involves the
        multiplication with the metric terms. Also multiply the fluxes with
        their integration weights and -theta/2. The parameter theta is the
        parameter in the Interior Penalty formulation, the factor 1/2 comes in
        from the averaging and the minus sign is from the convention that the
        viscous fluxes come with a minus sign in this code. ---*/

  /* Make a distinction between two and three space dimensions
     in order to have the most efficient code. */
  switch( nDim ) {

    case 2: {
      /* Two-dimensional computation. Loop over the integration
         points to correct the fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the location where the the fluxes
           are stored for this integration point as well as an
           abbreviation for the integration weights time halfTheta. */
        const unsigned short offPointer = i*nDim*NPad + indFaceChunk*nVar;
        const su2double *fluxX =  symmFluxes + offPointer;
        const su2double *fluxY =  fluxX + NPad;
        const su2double wTheta = -halfTheta*weights[i];
        su2double *paramFlux   =  paramFluxes + offPointer;

        /* Compute the modified metric terms. */
        const su2double *metricTerms = metricCoorFace + 4*i;   // The 4 is nDim*nDim;
        const su2double drdx = wTheta*metricTerms[0];
        const su2double drdy = wTheta*metricTerms[1];
        const su2double dsdx = wTheta*metricTerms[2];
        const su2double dsdy = wTheta*metricTerms[3];

        /* Parametric fluxes in r-direction. */
        paramFlux[0] = fluxX[0]*drdx + fluxY[0]*drdy;
        paramFlux[1] = fluxX[1]*drdx + fluxY[1]*drdy;
        paramFlux[2] = fluxX[2]*drdx + fluxY[2]*drdy;
        paramFlux[3] = fluxX[3]*drdx + fluxY[3]*drdy;

        /* Parametric fluxes in s-direction. */
        paramFlux    = paramFlux + NPad;
        paramFlux[0] = fluxX[0]*dsdx + fluxY[0]*dsdy;
        paramFlux[1] = fluxX[1]*dsdx + fluxY[1]*dsdy;
        paramFlux[2] = fluxX[2]*dsdx + fluxY[2]*dsdy;
        paramFlux[3] = fluxX[3]*dsdx + fluxY[3]*dsdy;
      }

      break;
    }

    case 3: {
      /* Three-dimensional computation. Loop over the integration
         points to correct the fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the location where the the fluxes
           are stored for this integration point as well as an
           abbreviation for the integration weights time halfTheta. */
        const unsigned short offPointer = i*nDim*NPad + indFaceChunk*nVar;
        const su2double *fluxX =  symmFluxes + offPointer;
        const su2double *fluxY =  fluxX + NPad;
        const su2double *fluxZ =  fluxY + NPad;
        const su2double wTheta = -halfTheta*weights[i];
        su2double *paramFlux   =  paramFluxes + offPointer;

        /* Compute the modified metric terms. */
        const su2double *metricTerms = metricCoorFace + 9*i;   // The 9 is nDim*nDim;
        const su2double drdx = wTheta*metricTerms[0];
        const su2double drdy = wTheta*metricTerms[1];
        const su2double drdz = wTheta*metricTerms[2];

        const su2double dsdx = wTheta*metricTerms[3];
        const su2double dsdy = wTheta*metricTerms[4];
        const su2double dsdz = wTheta*metricTerms[5];

        const su2double dtdx = wTheta*metricTerms[6];
        const su2double dtdy = wTheta*metricTerms[7];
        const su2double dtdz = wTheta*metricTerms[8];

        /* Parametric fluxes in r-direction. */
        paramFlux[0] = fluxX[0]*drdx + fluxY[0]*drdy + fluxZ[0]*drdz;
        paramFlux[1] = fluxX[1]*drdx + fluxY[1]*drdy + fluxZ[1]*drdz;
        paramFlux[2] = fluxX[2]*drdx + fluxY[2]*drdy + fluxZ[2]*drdz;
        paramFlux[3] = fluxX[3]*drdx + fluxY[3]*drdy + fluxZ[3]*drdz;
        paramFlux[4] = fluxX[4]*drdx + fluxY[4]*drdy + fluxZ[4]*drdz;

        /* Parametric fluxes in s-direction. */
        paramFlux    = paramFlux + NPad;
        paramFlux[0] = fluxX[0]*dsdx + fluxY[0]*dsdy + fluxZ[0]*dsdz;
        paramFlux[1] = fluxX[1]*dsdx + fluxY[1]*dsdy + fluxZ[1]*dsdz;
        paramFlux[2] = fluxX[2]*dsdx + fluxY[2]*dsdy + fluxZ[2]*dsdz;
        paramFlux[3] = fluxX[3]*dsdx + fluxY[3]*dsdy + fluxZ[3]*dsdz;
        paramFlux[4] = fluxX[4]*dsdx + fluxY[4]*dsdy + fluxZ[4]*dsdz;

        /* Parametric fluxes in t-direction. */
        paramFlux    = paramFlux + NPad;
        paramFlux[0] = fluxX[0]*dtdx + fluxY[0]*dtdy + fluxZ[0]*dtdz;
        paramFlux[1] = fluxX[1]*dtdx + fluxY[1]*dtdy + fluxZ[1]*dtdz;
        paramFlux[2] = fluxX[2]*dtdx + fluxY[2]*dtdy + fluxZ[2]*dtdz;
        paramFlux[3] = fluxX[3]*dtdx + fluxY[3]*dtdy + fluxZ[3]*dtdz;
        paramFlux[4] = fluxX[4]*dtdx + fluxY[4]*dtdy + fluxZ[4]*dtdz;
      }

      break;
    }
  }
}

void CFEM_DG_NSSolver::BC_Euler_Wall(CConfig                  *config,
                                     const unsigned long      surfElemBeg,
                                     const unsigned long      surfElemEnd,
                                     const CSurfaceElementFEM *surfElem,
                                     su2double                *resFaces,
                                     CNumerics                *conv_numerics,
                                     su2double                *workArray){

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /*--- Get the information from the standard element, which is the same
          for all the faces in the chunks considered. ---*/
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array workarray is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Compute the right state by applying the inviscid wall BC's. */
    BoundaryStates_Euler_Wall(config, llEnd, NPad, &surfElem[l],
                              solIntL, solIntR);

    /* The remainder of the boundary treatment is the same for all
       boundary conditions (except the symmetry plane). */
    ViscousBoundaryFacesBCTreatment(config, conv_numerics, llEnd, NPad,
                                    0.0, false, 0.0, false, &surfElem[l],
                                    solIntL, solIntR, work, resFaces,
                                    indResFaces, nullptr);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::BC_Far_Field(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    su2double                *workArray){

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /*--- Get the information from the standard element, which is the same
          for all the faces in the chunks considered. ---*/
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array workarray is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Set the right state in the integration points to the free stream value. */
    for(unsigned short ll=0; ll<llEnd; ++ll) {
      const unsigned short llNVar = ll*nVar;
      for(unsigned short i=0; i<nInt; ++i) {
        su2double *UR = solIntR + i*NPad + llNVar;
        for(unsigned short j=0; j<nVar; ++j)
          UR[j] = ConsVarFreeStream[j];
      }
    }

    /* The remainder of the boundary treatment is the same for all
       boundary conditions (except the symmetry plane). */
    ViscousBoundaryFacesBCTreatment(config, conv_numerics, llEnd, NPad,
                                    0.0, false, 0.0, false, &surfElem[l],
                                    solIntL, solIntR, work, resFaces,
                                    indResFaces, nullptr);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::BC_Sym_Plane(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    su2double                *workArray){

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /*--------------------------------------------------------------------------*/
    /*--- Step 1: Initialization and computation of the left states.         ---*/
    /*--------------------------------------------------------------------------*/

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /*--- Get the information from the standard element, which is the same
          for all the faces in the chunks considered. ---*/
    const unsigned short nInt         = standardBoundaryFacesSol[ind].GetNIntegration();
    const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
    const su2double     *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL      = workArray;
    su2double *solIntR      = solIntL + NPad*nInt;
    su2double *viscosityInt = solIntR + NPad*nInt;
    su2double *kOverCvInt   = viscosityInt + llEnd*nInt;
    su2double *gradSolInt   = kOverCvInt   + llEnd*nInt;
    su2double *fluxes       = gradSolInt   + NPad*nInt*nDim;
    su2double *viscFluxes   = fluxes       + NPad*max(nInt*nDim, (int) nDOFsElem);

    /* Compute the left states in the integration points of this chunk of
       faces.  Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            fluxes, solIntL);

    /*-----------------------------------------------------------------------*/
    /*--- Step 2: Computation of the gradients of the conservative        ---*/
    /*---         variables w.r.t. the parametric coordinates.            ---*/
    /*-----------------------------------------------------------------------*/

    /* Set the pointer solElem to fluxes, just for readability. It is a temporary
       storage of the solution of the elements, such that the gradients can be
       computed efficiently. */
    su2double *solElem = fluxes;

    /* Loop over the faces of the chunk to set the conservative variables
       in the appropriate sequence. */
    for(unsigned long ll=l; ll<lEnd; ++ll) {
      const unsigned short llRel  = ll - l;
      const unsigned short llNVar = llRel*nVar;

      /* Determine the time level of the boundary face, which is equal
         to the time level of the adjacent element. */
      const unsigned long  elemID    = surfElem[ll].volElemID;
      const unsigned short timeLevel = volElem[elemID].timeLevel;

      /* Determine the offset that must be applied to access the correct data for
         the adjacent element in the working vector of the solution. This is a
         boundary face, which has the same time level as the adjacent element. */
      const unsigned long offset = volElem[elemID].offsetDOFsSolLocal
                                 - volElem[elemID].offsetDOFsSolThisTimeLevel;

      /* Loop over the DOFs of the element and copy the solution. */
      for(unsigned short i=0; i<nDOFsElem; ++i) {
        const su2double *solDOF = VecWorkSolDOFs[timeLevel].data()
                                + nVar*(surfElem[ll].DOFsSolElement[i] - offset);
        for(unsigned short mm=0; mm<nVar; ++mm)
          solElem[NPad*i+llNVar+mm] = solDOF[mm];
      }
    }

    /* Compute the left gradients in the integration points. Call the general
       function to carry out the matrix product. */
    blasFunctions->gemm(nInt*nDim, NPad, nDOFsElem, derBasisElem, solElem, gradSolInt, config);

    /*-----------------------------------------------------------------------*/
    /*--- Step 3: Computation of the viscous fluxes in the integration    ---*/
    /*---         points.                                                 ---*/
    /*-----------------------------------------------------------------------*/

    /* Determine the offset between r- and -s-derivatives, which is also the
       offset between s- and t-derivatives. */
    const unsigned short offDeriv = NPad*nInt;

    /* Loop over the faces of the chunk. */
    for(unsigned long ll=l; ll<lEnd; ++ll) {
      const unsigned short llRel  = ll - l;
      const unsigned short llNVar = llRel*nVar;

      /* Compute the length scale for the LES of the adjacent element. */
      const unsigned long  elemID = surfElem[ll].volElemID;
      const unsigned short iind   = volElem[elemID].indStandardElement;

      unsigned short nPoly = standardElementsSol[iind].GetNPoly();
      if(nPoly == 0) nPoly = 1;

      const su2double lenScale_LES = volElem[elemID].lenScale/nPoly;

      /* Easier storage of the wall distance array for this surface element. */
      const su2double *wallDistance = surfElem[ll].wallDistance.data();

      /* Make a distinction between two and three space dimensions
         in order to have the most efficient code. */
      switch( nDim ) {

        case 2: {

          /*--- 2D simulation. Loop over the integration points to apply the
                symmetry condition and to compute the viscous normal fluxes. This
                is done in the same loop, because the gradients in the right
                integration points are constructed using the symmetry boundary
                condition. ---*/
          for(unsigned short i=0; i<nInt; ++i) {

            /* Easier storage of the left and right solution and the normals
               for this integration point. */
            const unsigned short pointerOffset = NPad*i + llNVar;
            const su2double *UL      = solIntL + pointerOffset;
                  su2double *UR      = solIntR + pointerOffset;
            const su2double *normals = surfElem[ll].metricNormalsFace.data() + i*(nDim+1);

            /* Compute twice the normal component of the momentum variables. The
               factor 2 comes from the fact that the velocity must be mirrored. */
            const su2double rVn = 2.0*(UL[1]*normals[0] + UL[2]*normals[1]);

            /* Set the right state. Note that the total energy of the right state is
               identical to the left state, because the magnitude of the velocity
               remains the same. */
            UR[0] = UL[0];
            UR[1] = UL[1] - rVn*normals[0];
            UR[2] = UL[2] - rVn*normals[1];
            UR[3] = UL[3];

            /* Easier storage of the metric terms and left gradients of this
               integration point. */
            const su2double *metricTerms = surfElem[ll].metricCoorDerivFace.data()
                                         + i*nDim*nDim;
            const su2double *dULDr       = gradSolInt + pointerOffset;
            const su2double *dULDs       = dULDr + offDeriv;

            /* Easier storage of the metric terms in this integration point. */
            const su2double drdx = metricTerms[0];
            const su2double drdy = metricTerms[1];
            const su2double dsdx = metricTerms[2];
            const su2double dsdy = metricTerms[3];

            /* Compute the Cartesian gradients of the left solution. */
            su2double ULGradCart[4][2];
            ULGradCart[0][0] = dULDr[0]*drdx + dULDs[0]*dsdx;
            ULGradCart[1][0] = dULDr[1]*drdx + dULDs[1]*dsdx;
            ULGradCart[2][0] = dULDr[2]*drdx + dULDs[2]*dsdx;
            ULGradCart[3][0] = dULDr[3]*drdx + dULDs[3]*dsdx;

            ULGradCart[0][1] = dULDr[0]*drdy + dULDs[0]*dsdy;
            ULGradCart[1][1] = dULDr[1]*drdy + dULDs[1]*dsdy;
            ULGradCart[2][1] = dULDr[2]*drdy + dULDs[2]*dsdy;
            ULGradCart[3][1] = dULDr[3]*drdy + dULDs[3]*dsdy;

            /* Determine the normal gradients of all variables. */
            su2double ULGradNorm[4];
            ULGradNorm[0] = ULGradCart[0][0]*normals[0] + ULGradCart[0][1]*normals[1];
            ULGradNorm[1] = ULGradCart[1][0]*normals[0] + ULGradCart[1][1]*normals[1];
            ULGradNorm[2] = ULGradCart[2][0]*normals[0] + ULGradCart[2][1]*normals[1];
            ULGradNorm[3] = ULGradCart[3][0]*normals[0] + ULGradCart[3][1]*normals[1];

            /* For the construction of the gradients of the right solution, also
               the Cartesian gradients and the normal gradient of the normal
               momentum is needed. This is computed below. */
            su2double GradCartNormMomL[2];
            GradCartNormMomL[0] = ULGradCart[1][0]*normals[0] + ULGradCart[2][0]*normals[1];
            GradCartNormMomL[1] = ULGradCart[1][1]*normals[0] + ULGradCart[2][1]*normals[1];

            const su2double GradNormNormMomL = ULGradNorm[1]*normals[0] + ULGradNorm[2]*normals[1];

            /* Abbreviate twice the normal vector. */
            const su2double tnx = 2.0*normals[0], tny = 2.0*normals[1];

            /*--- Construct the gradients of the right solution. The tangential
                  gradients of the normal momentum and normal gradients of the other
                  variables, density, energy and tangential momentum, must be negated.
                  For the common situation that the symmetry plane coincides with an
                  x- or y-plane this boils down to a multiplication by -1 or +1,
                  depending on the variable and gradient. However, the implementation
                  below is valid for an arbitrary orientation of the symmetry plane. ---*/
            su2double URGradCart[4][2];
            URGradCart[0][0] = ULGradCart[0][0] - tnx*ULGradNorm[0];
            URGradCart[1][0] = ULGradCart[1][0] - tnx*ULGradNorm[1] - tnx*GradCartNormMomL[0] + tnx*tnx*GradNormNormMomL;
            URGradCart[2][0] = ULGradCart[2][0] - tnx*ULGradNorm[2] - tny*GradCartNormMomL[0] + tnx*tny*GradNormNormMomL;
            URGradCart[3][0] = ULGradCart[3][0] - tnx*ULGradNorm[3];

            URGradCart[0][1] = ULGradCart[0][1] - tny*ULGradNorm[0];
            URGradCart[1][1] = ULGradCart[1][1] - tny*ULGradNorm[1] - tnx*GradCartNormMomL[1] + tny*tnx*GradNormNormMomL;
            URGradCart[2][1] = ULGradCart[2][1] - tny*ULGradNorm[2] - tny*GradCartNormMomL[1] + tny*tny*GradNormNormMomL;
            URGradCart[3][1] = ULGradCart[3][1] - tny*ULGradNorm[3];

            /*--- Compute the viscous fluxes of the left and right state and average
                  them to compute the correct viscous flux for a symmetry boundary. ---*/
            const su2double wallDist = wallDistance ? wallDistance[i] : 0.0;
            su2double viscFluxL[4], viscFluxR[4];
            su2double Viscosity, kOverCv;

            ViscousNormalFluxIntegrationPoint_2D(UR, URGradCart, normals, 0.0, 1.0,
                                                 wallDist, lenScale_LES, Viscosity,
                                                 kOverCv, viscFluxR);
            ViscousNormalFluxIntegrationPoint_2D(UL, ULGradCart, normals, 0.0, 1.0,
                                                 wallDist, lenScale_LES, Viscosity,
                                                 kOverCv, viscFluxL);
            viscosityInt[nInt*llRel+i] = Viscosity;
            kOverCvInt[nInt*llRel+i]   = kOverCv;

            su2double *viscNormalFlux = viscFluxes + pointerOffset;
            viscNormalFlux[0] = 0.5*(viscFluxL[0] + viscFluxR[0]);
            viscNormalFlux[1] = 0.5*(viscFluxL[1] + viscFluxR[1]);
            viscNormalFlux[2] = 0.5*(viscFluxL[2] + viscFluxR[2]);
            viscNormalFlux[3] = 0.5*(viscFluxL[3] + viscFluxR[3]);
          }

          break;
        }

        /*--------------------------------------------------------------------*/

        case 3: {

          /*--- 3D simulation. Loop over the integration points to apply the
                symmetry condition and to compute the viscous normal fluxes. This
                is done in the same loop, because the gradients in the right
                integration points are constructed using the symmetry boundary
                condition. ---*/
          for(unsigned short i=0; i<nInt; ++i) {

            /* Easier storage of the left and right solution and the normals
               for this integration point. */
            const unsigned short pointerOffset = NPad*i + llNVar;
            const su2double *UL      = solIntL + pointerOffset;
                  su2double *UR      = solIntR + pointerOffset;
            const su2double *normals = surfElem[ll].metricNormalsFace.data() + i*(nDim+1);

            /* Compute twice the normal component of the momentum variables. The
               factor 2 comes from the fact that the velocity must be mirrored. */
            const su2double rVn = 2.0*(UL[1]*normals[0] + UL[2]*normals[1] + UL[3]*normals[2]);

            /* Set the right state. Note that the total energy of the right state is
               identical to the left state, because the magnitude of the velocity
               remains the same. */
            UR[0] = UL[0];
            UR[1] = UL[1] - rVn*normals[0];
            UR[2] = UL[2] - rVn*normals[1];
            UR[3] = UL[3] - rVn*normals[2];
            UR[4] = UL[4];

            /* Easier storage of the metric terms and left gradients of this
               integration point. */
            const su2double *metricTerms = surfElem[ll].metricCoorDerivFace.data()
                                         + i*nDim*nDim;
            const su2double *dULDr       = gradSolInt + pointerOffset;
            const su2double *dULDs       = dULDr + offDeriv;
            const su2double *dULDt       = dULDs + offDeriv;

            /* Easier storage of the metric terms in this integration point. */
            const su2double drdx = metricTerms[0];
            const su2double drdy = metricTerms[1];
            const su2double drdz = metricTerms[2];

            const su2double dsdx = metricTerms[3];
            const su2double dsdy = metricTerms[4];
            const su2double dsdz = metricTerms[5];

            const su2double dtdx = metricTerms[6];
            const su2double dtdy = metricTerms[7];
            const su2double dtdz = metricTerms[8];

            /* Compute the Cartesian gradients of the left solution. */
            su2double ULGradCart[5][3];
            ULGradCart[0][0] = dULDr[0]*drdx + dULDs[0]*dsdx + dULDt[0]*dtdx;
            ULGradCart[1][0] = dULDr[1]*drdx + dULDs[1]*dsdx + dULDt[1]*dtdx;
            ULGradCart[2][0] = dULDr[2]*drdx + dULDs[2]*dsdx + dULDt[2]*dtdx;
            ULGradCart[3][0] = dULDr[3]*drdx + dULDs[3]*dsdx + dULDt[3]*dtdx;
            ULGradCart[4][0] = dULDr[4]*drdx + dULDs[4]*dsdx + dULDt[4]*dtdx;

            ULGradCart[0][1] = dULDr[0]*drdy + dULDs[0]*dsdy + dULDt[0]*dtdy;
            ULGradCart[1][1] = dULDr[1]*drdy + dULDs[1]*dsdy + dULDt[1]*dtdy;
            ULGradCart[2][1] = dULDr[2]*drdy + dULDs[2]*dsdy + dULDt[2]*dtdy;
            ULGradCart[3][1] = dULDr[3]*drdy + dULDs[3]*dsdy + dULDt[3]*dtdy;
            ULGradCart[4][1] = dULDr[4]*drdy + dULDs[4]*dsdy + dULDt[4]*dtdy;

            ULGradCart[0][2] = dULDr[0]*drdz + dULDs[0]*dsdz + dULDt[0]*dtdz;
            ULGradCart[1][2] = dULDr[1]*drdz + dULDs[1]*dsdz + dULDt[1]*dtdz;
            ULGradCart[2][2] = dULDr[2]*drdz + dULDs[2]*dsdz + dULDt[2]*dtdz;
            ULGradCart[3][2] = dULDr[3]*drdz + dULDs[3]*dsdz + dULDt[3]*dtdz;
            ULGradCart[4][2] = dULDr[4]*drdz + dULDs[4]*dsdz + dULDt[4]*dtdz;

            /* Determine the normal gradients of all variables. */
            su2double ULGradNorm[5];
            ULGradNorm[0] = ULGradCart[0][0]*normals[0] + ULGradCart[0][1]*normals[1] + ULGradCart[0][2]*normals[2];
            ULGradNorm[1] = ULGradCart[1][0]*normals[0] + ULGradCart[1][1]*normals[1] + ULGradCart[1][2]*normals[2];
            ULGradNorm[2] = ULGradCart[2][0]*normals[0] + ULGradCart[2][1]*normals[1] + ULGradCart[2][2]*normals[2];
            ULGradNorm[3] = ULGradCart[3][0]*normals[0] + ULGradCart[3][1]*normals[1] + ULGradCart[3][2]*normals[2];
            ULGradNorm[4] = ULGradCart[4][0]*normals[0] + ULGradCart[4][1]*normals[1] + ULGradCart[4][2]*normals[2];

            /* For the construction of the gradients of the right solution, also
               the Cartesian gradients and the normal gradient of the normal
               momentum is needed. This is computed below. */
            su2double GradCartNormMomL[3];
            GradCartNormMomL[0] = ULGradCart[1][0]*normals[0] + ULGradCart[2][0]*normals[1] + ULGradCart[3][0]*normals[2];
            GradCartNormMomL[1] = ULGradCart[1][1]*normals[0] + ULGradCart[2][1]*normals[1] + ULGradCart[3][1]*normals[2];
            GradCartNormMomL[2] = ULGradCart[1][2]*normals[0] + ULGradCart[2][2]*normals[1] + ULGradCart[3][2]*normals[2];

            const su2double GradNormNormMomL = ULGradNorm[1]*normals[0]
                                             + ULGradNorm[2]*normals[1]
                                             + ULGradNorm[3]*normals[2];

            /* Abbreviate twice the normal vector. */
            const su2double tnx = 2.0*normals[0], tny = 2.0*normals[1], tnz = 2.0*normals[2];

            /*--- Construct the gradients of the right solution. The tangential
                  gradients of the normal momentum and normal gradients of the other
                  variables, density, energy and tangential momentum, must be negated.
                  For the common situation that the symmetry plane coincides with an
                  x-, y- or z-plane this boils down to a multiplication by -1 or +1,
                  depending on the variable and gradient. However, the implementation
                  below is valid for an arbitrary orientation of the symmetry plane. ---*/
            su2double URGradCart[5][3];
            URGradCart[0][0] = ULGradCart[0][0] - tnx*ULGradNorm[0];
            URGradCart[1][0] = ULGradCart[1][0] - tnx*ULGradNorm[1] - tnx*GradCartNormMomL[0] + tnx*tnx*GradNormNormMomL;
            URGradCart[2][0] = ULGradCart[2][0] - tnx*ULGradNorm[2] - tny*GradCartNormMomL[0] + tnx*tny*GradNormNormMomL;
            URGradCart[3][0] = ULGradCart[3][0] - tnx*ULGradNorm[3] - tnz*GradCartNormMomL[0] + tnx*tnz*GradNormNormMomL;
            URGradCart[4][0] = ULGradCart[4][0] - tnx*ULGradNorm[4];

            URGradCart[0][1] = ULGradCart[0][1] - tny*ULGradNorm[0];
            URGradCart[1][1] = ULGradCart[1][1] - tny*ULGradNorm[1] - tnx*GradCartNormMomL[1] + tny*tnx*GradNormNormMomL;
            URGradCart[2][1] = ULGradCart[2][1] - tny*ULGradNorm[2] - tny*GradCartNormMomL[1] + tny*tny*GradNormNormMomL;
            URGradCart[3][1] = ULGradCart[3][1] - tny*ULGradNorm[3] - tnz*GradCartNormMomL[1] + tny*tnz*GradNormNormMomL;
            URGradCart[4][1] = ULGradCart[4][1] - tny*ULGradNorm[4];

            URGradCart[0][2] = ULGradCart[0][2] - tnz*ULGradNorm[0];
            URGradCart[1][2] = ULGradCart[1][2] - tnz*ULGradNorm[1] - tnx*GradCartNormMomL[2] + tnz*tnx*GradNormNormMomL;
            URGradCart[2][2] = ULGradCart[2][2] - tnz*ULGradNorm[2] - tny*GradCartNormMomL[2] + tnz*tny*GradNormNormMomL;
            URGradCart[3][2] = ULGradCart[3][2] - tnz*ULGradNorm[3] - tnz*GradCartNormMomL[2] + tnz*tnz*GradNormNormMomL;
            URGradCart[4][2] = ULGradCart[4][2] - tnz*ULGradNorm[4];

            /*--- Compute the viscous fluxes of the left and right state and average
                  them to compute the correct viscous flux for a symmetry boundary. ---*/
            const su2double wallDist = wallDistance ? wallDistance[i] : 0.0;
            su2double viscFluxL[5], viscFluxR[5];
            su2double Viscosity, kOverCv;

            ViscousNormalFluxIntegrationPoint_3D(UR, URGradCart, normals, 0.0, 1.0,
                                                 wallDist, lenScale_LES, Viscosity,
                                                 kOverCv, viscFluxR);
            ViscousNormalFluxIntegrationPoint_3D(UL, ULGradCart, normals, 0.0, 1.0,
                                                 wallDist, lenScale_LES, Viscosity,
                                                 kOverCv, viscFluxL);
            viscosityInt[nInt*llRel+i] = Viscosity;
            kOverCvInt[nInt*llRel+i]   = kOverCv;

            su2double *viscNormalFlux = viscFluxes + pointerOffset;
            viscNormalFlux[0] = 0.5*(viscFluxL[0] + viscFluxR[0]);
            viscNormalFlux[1] = 0.5*(viscFluxL[1] + viscFluxR[1]);
            viscNormalFlux[2] = 0.5*(viscFluxL[2] + viscFluxR[2]);
            viscNormalFlux[3] = 0.5*(viscFluxL[3] + viscFluxR[3]);
            viscNormalFlux[4] = 0.5*(viscFluxL[4] + viscFluxR[4]);
          }

          break;
        }
      }
    }

    /* The remainder of the boundary condition treatment is the same for all
       types of boundary conditions, including the symmetry plane. The function
       ResidualViscousBoundaryFace will carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, llEnd, NPad, &surfElem[l],
                                solIntL, solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::BC_Supersonic_Outlet(CConfig                  *config,
                                            const unsigned long      surfElemBeg,
                                            const unsigned long      surfElemEnd,
                                            const CSurfaceElementFEM *surfElem,
                                            su2double                *resFaces,
                                            CNumerics                *conv_numerics,
                                            su2double                *workArray){

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /*--- Get the information from the standard element, which is the same
          for all the faces in the chunks considered. ---*/
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array workarray is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Set the right state in the integration points to the left state, i.e.
       no boundary condition is applied for a supersonic outlet. */
    for(unsigned short mm=0; mm<(NPad*nInt); ++mm)
      solIntR[mm] = solIntL[mm];

    /* The remainder of the boundary treatment is the same for all
       boundary conditions (except the symmetry plane). */
    ViscousBoundaryFacesBCTreatment(config, conv_numerics, llEnd, NPad,
                                    0.0, false, 0.0, false, &surfElem[l],
                                    solIntL, solIntR, work, resFaces,
                                    indResFaces, nullptr);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::BC_Inlet(CConfig                  *config,
                                const unsigned long      surfElemBeg,
                                const unsigned long      surfElemEnd,
                                const CSurfaceElementFEM *surfElem,
                                su2double                *resFaces,
                                CNumerics                *conv_numerics,
                                unsigned short           val_marker,
                                su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /*--- Get the information from the standard element, which is the same
          for all the faces in the chunks considered. ---*/
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array workarray is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Compute the right state by applying the subsonic inlet BC's. */
    BoundaryStates_Inlet(config, llEnd, NPad, &surfElem[l], val_marker, solIntL, solIntR);

    /* The remainder of the boundary treatment is the same for all
       boundary conditions (except the symmetry plane). */
    ViscousBoundaryFacesBCTreatment(config, conv_numerics, llEnd, NPad,
                                    0.0, false, 0.0, false, &surfElem[l],
                                    solIntL, solIntR, work, resFaces,
                                    indResFaces, nullptr);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::BC_Outlet(CConfig                  *config,
                                 const unsigned long      surfElemBeg,
                                 const unsigned long      surfElemEnd,
                                 const CSurfaceElementFEM *surfElem,
                                 su2double                *resFaces,
                                 CNumerics                *conv_numerics,
                                 unsigned short           val_marker,
                                 su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /*--- Get the information from the standard element, which is the same
          for all the faces in the chunks considered. ---*/
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array workarray is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Compute the right state by applying the subsonic outlet BC's. */
    BoundaryStates_Outlet(config, llEnd, NPad, &surfElem[l], val_marker, solIntL, solIntR);

    /* The remainder of the boundary treatment is the same for all
       boundary conditions (except the symmetry plane). */
    ViscousBoundaryFacesBCTreatment(config, conv_numerics, llEnd, NPad,
                                    0.0, false, 0.0, false, &surfElem[l],
                                    solIntL, solIntR, work, resFaces,
                                    indResFaces, nullptr);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::BC_HeatFlux_Wall(CConfig                  *config,
                                        const unsigned long      surfElemBeg,
                                        const unsigned long      surfElemEnd,
                                        const CSurfaceElementFEM *surfElem,
                                        su2double                *resFaces,
                                        CNumerics                *conv_numerics,
                                        unsigned short           val_marker,
                                        su2double                *workArray) {

  /* Set the factor for the wall velocity. For factWallVel = 0, the right state
     contains the wall velocity. For factWallVel = 1.0, the velocity of the
     right state is obtained by negating the interior velocity w.r.t. the
     velocity of the wall. */
  const su2double factWallVel = 0.0;
  // const su2double factWallVel = 1.0;

  /* Get the wall heat flux. */
  const string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
  const su2double Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /*--- Get the information from the standard element, which is the same
          for all the faces in the chunks considered. ---*/
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array workarray is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /*--- Compute the right states. A distinction must be made whether or not
          all functions are used. ---*/
    if( boundaries[val_marker].wallModel ) {

      /*--- Wall functions are used to model the lower part of the boundary
            layer. Consequently, a slip boundary condition is used for the
            velocities and the skin friction and heat flux from the wall
            should drive the velocities towards the no-slip and prescribed
            heat flux. Therefore the right states for the actual flow solver
            can be computed using BoundaryStates_Euler_Wall. ---*/
      BoundaryStates_Euler_Wall(config, llEnd, NPad, &surfElem[l],
                                solIntL, solIntR);
    }
    else {

      /*--- Integration to the wall is used, so the no-slip condition is enforced,
            albeit weakly. Loop over the number of faces in this chunk and the
            number of integration points and apply the heat flux wall boundary
            conditions to compute the right state. There are two options. Either
            the velocity is negated or it is set to zero. Some experiments are
            needed to see which formulation gives better results. ---*/
      for(unsigned short ll=0; ll<llEnd; ++ll) {
        const unsigned long  lll    = l + ll;
        const unsigned short llNVar = ll*nVar;

        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the grid velocity and the left and right solution
             for this integration point. */
          const su2double *gridVel = surfElem[lll].gridVelocities.data() + i*nDim;
          const su2double *UL      = solIntL + NPad*i + llNVar;
                su2double *UR      = solIntR + NPad*i + llNVar;

          /* Set the right state. The initial value of the total energy is the
             energy of the left state. Also compute the difference in kinetic
             energy between the left and right state. */
          UR[0]      = UL[0];
          UR[nDim+1] = UL[nDim+1];

          su2double DensityInv = 1.0/UL[0];
          su2double diffKin    = 0.0;
          for(unsigned short iDim=0; iDim<nDim; ++iDim) {
            const su2double velL = DensityInv*UL[iDim+1];
            const su2double dV   = factWallVel*(velL-gridVel[iDim]);
            const su2double velR = gridVel[iDim] - dV;

            UR[iDim+1] = UR[0]*velR;
            diffKin   += velL*velL - velR*velR;
          }

          /* As only the internal energy of UR is equal to UL, the difference
             in kinetic energy must be subtracted. */
          UR[nDim+1] -= 0.5*UR[0]*diffKin;
        }
      }
    }

    /* The remainder of the boundary treatment is the same for all
       boundary conditions (except the symmetry plane). */
    ViscousBoundaryFacesBCTreatment(config, conv_numerics, llEnd, NPad,
                                    Wall_HeatFlux, true, 0.0, false,
                                    &surfElem[l], solIntL, solIntR,
                                    work, resFaces, indResFaces,
                                    boundaries[val_marker].wallModel);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::BC_Isothermal_Wall(CConfig                  *config,
                                          const unsigned long      surfElemBeg,
                                          const unsigned long      surfElemEnd,
                                          const CSurfaceElementFEM *surfElem,
                                          su2double                *resFaces,
                                          CNumerics                *conv_numerics,
                                          unsigned short           val_marker,
                                          su2double                *workArray) {

  /* Set the factor for the wall velocity. For factWallVel = 0, the right state
     contains the wall velocity. For factWallVel = 1.0, the velocity of the
     right state is obtained by negating the interior velocity w.r.t. the
     velocity of the wall. */
  const su2double factWallVel = 0.0;
  // const su2double factWallVel = 1.0;

  /* Get the wall temperature. */
  const string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  const su2double TWall   = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

  /* Compute the prescribed value of the energy (per unit mass). */
  const su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double Cv           = Gas_Constant/Gamma_Minus_One;
  const su2double StaticEnergy = Cv*TWall;

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /*--- Get the information from the standard element, which is the same
          for all the faces in the chunks considered. ---*/
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array workarray is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /*--- Integration to the wall is used, so the no-slip condition is enforced,
          albeit weakly. Loop over the number of faces in this chunk and the
          number of integration points and apply the isothermal wall boundary
          conditions to compute the right state. There are two options. Either
          the velocity is negated or it is set to zero. Some experiments are
          needed to see which formulation gives better results. ---*/
    for(unsigned short ll=0; ll<llEnd; ++ll) {
      const unsigned long  lll    = l + ll;
      const unsigned short llNVar = ll*nVar;

      for(unsigned short i=0; i<nInt; ++i) {

         /* Easier storage of the grid velocity and the left and right solution
           for this integration point. */
        const su2double *gridVel = surfElem[lll].gridVelocities.data() + i*nDim;
        const su2double *UL      = solIntL + NPad*i + llNVar;
              su2double *UR      = solIntR + NPad*i + llNVar;

        /* Set the right state for the density and the momentum variables of the
           right state. Compute twice the possible kinetic energy. */
        UR[0] = UL[0];
        su2double DensityInv = 1.0/UL[0];
        su2double kinEner    = 0.0;
        for(unsigned short iDim=0; iDim<nDim; ++iDim) {
          const su2double velL = DensityInv*UL[iDim+1];
          const su2double dV   = factWallVel*(velL-gridVel[iDim]);
          const su2double velR = gridVel[iDim] - dV;

          UR[iDim+1] = UR[0]*velR;
          kinEner   += velR*velR;
        }

        /* Compute the total energy of the right state. */
        UR[nDim+1] = UR[0]*(StaticEnergy + 0.5*kinEner);
      }
    }

    /* The remainder of the boundary treatment is the same for all
       boundary conditions (except the symmetry plane). */
    ViscousBoundaryFacesBCTreatment(config, conv_numerics, llEnd, NPad,
                                    0.0, false, TWall, true,
                                    &surfElem[l], solIntL, solIntR,
                                    work, resFaces, indResFaces,
                                    boundaries[val_marker].wallModel);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::BC_Riemann(CConfig                  *config,
                                  const unsigned long      surfElemBeg,
                                  const unsigned long      surfElemEnd,
                                  const CSurfaceElementFEM *surfElem,
                                  su2double                *resFaces,
                                  CNumerics                *conv_numerics,
                                  unsigned short           val_marker,
                                  su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /*--- Get the information from the standard element, which is the same
          for all the faces in the chunks considered. ---*/
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array workarray is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Compute the right state by applying the Riemann BC's. */
    BoundaryStates_Riemann(config, llEnd, NPad, &surfElem[l], val_marker, solIntL, solIntR);

    /* The remainder of the boundary treatment is the same for all
       boundary conditions (except the symmetry plane). */
    ViscousBoundaryFacesBCTreatment(config, conv_numerics, llEnd, NPad,
                                    0.0, false, 0.0, false, &surfElem[l],
                                    solIntL, solIntR, work,
                                    resFaces, indResFaces, nullptr);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::BC_Custom(CConfig                  *config,
                                 const unsigned long      surfElemBeg,
                                 const unsigned long      surfElemEnd,
                                 const CSurfaceElementFEM *surfElem,
                                 su2double                *resFaces,
                                 CNumerics                *conv_numerics,
                                 su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Get the physical time if necessary. ---*/
  su2double time = 0.0;
  if (config->GetTime_Marching() != TIME_MARCHING::STEADY) time = config->GetPhysicalTime();

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /*--- Get the information from the standard element, which is the same
          for all the faces in the chunks considered. ---*/
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array workarray is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Check for a verification solution. */
    if( VerificationSolution ) {

      /*--- Loop over the number of simultaneously treated faces and integration points
            to compute the right state for the boundary conditions. ---*/
      for(unsigned short ll=0; ll<llEnd; ++ll) {
        for(unsigned short i=0; i<nInt; ++i) {

          /* Determine the pointer to the coordinates of this integration
             point and the pointer to the solution and call the function
             GetBCState to determine the actual boundary state. */
          const su2double *coor = surfElem[ll+l].coorIntegrationPoints.data() + i*nDim;
          su2double *UR   = solIntR + NPad*i + ll*nVar;

          VerificationSolution->GetBCState(coor, time, UR);
        }
      }
    }
    else {

      /* The user must specify the custom BC's here. */
      SU2_MPI::Error("Implement customized boundary conditions here.", CURRENT_FUNCTION);
    }


    /* The remainder of the boundary treatment is the same for all
       boundary conditions (except the symmetry plane). */
    ViscousBoundaryFacesBCTreatment(config, conv_numerics, llEnd, NPad,
                                    0.0, false, 0.0, false, &surfElem[l],
                                    solIntL, solIntR, work,
                                    resFaces, indResFaces, nullptr);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_NSSolver::ViscousBoundaryFacesBCTreatment(
                                       CConfig                  *config,
                                       CNumerics                *conv_numerics,
                                       const unsigned short     nFaceSimul,
                                       const unsigned short     NPad,
                                       const su2double          Wall_HeatFlux,
                                       const bool               HeatFlux_Prescribed,
                                       const su2double          Wall_Temperature,
                                       const bool               Temperature_Prescribed,
                                       const CSurfaceElementFEM *surfElem,
                                       const su2double          *solIntL,
                                       const su2double          *solIntR,
                                             su2double          *workArray,
                                             su2double          *resFaces,
                                             unsigned long      &indResFaces,
                                             CWallModel         *wallModel) {

  /*--- Get the information from the standard element, which is the same
        for all the faces in the chunks considered. ---*/
  const unsigned short ind       = surfElem[0].indStandardElement;
  const unsigned short nInt      = standardBoundaryFacesSol[ind].GetNIntegration();
  const unsigned short nDOFsElem = standardBoundaryFacesSol[ind].GetNDOFsElem();
  const su2double *derBasisElem  = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

  /*--- Set the pointers for the local arrays. ---*/
  su2double *viscosityInt = workArray;
  su2double *kOverCvInt   = viscosityInt + nFaceSimul*nInt;
  su2double *gradSolInt   = kOverCvInt   + nFaceSimul*nInt;
  su2double *fluxes       = gradSolInt   + NPad*nInt*nDim;
  su2double *viscFluxes   = fluxes       + NPad*max(nInt*nDim, (int) nDOFsElem);

  /* Compute the viscous fluxes in the integration points of the faces that
     are treated simulaneously. Make a distinction between a wall function
     treatment and a standard computation of the viscous fluxes. */
  if( wallModel ) {
    WallTreatmentViscousFluxes(config, nFaceSimul, NPad, nInt, Wall_HeatFlux,
                               HeatFlux_Prescribed, Wall_Temperature,
                               Temperature_Prescribed, surfElem, solIntL,
                               gradSolInt, viscFluxes, viscosityInt,
                               kOverCvInt, wallModel);
  }
  else {
    ComputeViscousFluxesBoundaryFaces(config, nFaceSimul, NPad, nInt, nDOFsElem,
                                      Wall_HeatFlux, HeatFlux_Prescribed,
                                      derBasisElem, surfElem, solIntL,
                                      fluxes, gradSolInt, viscFluxes,
                                      viscosityInt, kOverCvInt);
  }

  /* The remainder of the boundary condition treatment is the same for all
     types of boundary conditions, including the symmetry plane and the
     wall function treatment. The function ResidualViscousBoundaryFace will
     carry out this task. */

  ResidualViscousBoundaryFace(config, conv_numerics, nFaceSimul, NPad, surfElem,
                              solIntL, solIntR, gradSolInt, fluxes, viscFluxes,
                              viscosityInt, kOverCvInt, resFaces, indResFaces);
}

void CFEM_DG_NSSolver::ComputeViscousFluxesBoundaryFaces(
                                       CConfig                  *config,
                                       const unsigned short     nFaceSimul,
                                       const unsigned short     NPad,
                                       const unsigned short     nInt,
                                       const unsigned short     nDOFsElem,
                                       const su2double          Wall_HeatFlux,
                                       const bool               HeatFlux_Prescribed,
                                       const su2double          *derBasisElem,
                                       const CSurfaceElementFEM *surfElem,
                                       const su2double          *solIntL,
                                             su2double          *solElem,
                                             su2double          *gradSolInt,
                                             su2double          *viscFluxes,
                                             su2double          *viscosityInt,
                                             su2double          *kOverCvInt) {

  /*---------------------------------------------------------------------------*/
  /*--- Step 1: Compute the gradients of the conservative variables in the  ---*/
  /*---         integration points of the faces.                            ---*/
  /*---------------------------------------------------------------------------*/

  /* Loop over the simultaneously treated faces to set the solution of the elements. */
  for(unsigned short l=0; l<nFaceSimul; ++l) {
    const unsigned short llNVar = l*nVar;

    /* Determine the offset that must be applied to access the correct data for
       this element in the working vector of the solution. Per definition the
       boundary face has the same time level as the adjacent element, hence
       offsetDOFsSolThisTimeLevel is used. */
    const unsigned long elemID     = surfElem[l].volElemID;
    const unsigned short timeLevel = volElem[elemID].timeLevel;
    const unsigned long offset     = volElem[elemID].offsetDOFsSolLocal
                                   - volElem[elemID].offsetDOFsSolThisTimeLevel;

    /* Loop over the DOFs of the element and copy the solution to the appropriate
       position in solElem. */
    const unsigned long *DOFsElem = surfElem[l].DOFsSolElement.data();
    for(unsigned short i=0; i<nDOFsElem; ++i) {
      const su2double *solDOF = VecWorkSolDOFs[timeLevel].data()
                              + nVar*(DOFsElem[i] - offset);
      su2double       *sol    = solElem + NPad*i + llNVar;
      for(unsigned short mm=0; mm<nVar; ++mm)
        sol[mm] = solDOF[mm];
    }
  }

  /* Compute the gradients in the integration points. Call the general function to
     carry out the matrix product. */
  blasFunctions->gemm(nInt*nDim, NPad, nDOFsElem, derBasisElem, solElem, gradSolInt, config);

  /*---------------------------------------------------------------------------*/
  /*--- Step 2: Compute the viscous normal fluxes in the integration points ---*/
  /*---         of the faces.                                               ---*/
  /*---------------------------------------------------------------------------*/

  /* Loop over the simultaneously treated faces. */
  for(unsigned short l=0; l<nFaceSimul; ++l) {

    /* Determine the ID of the adjacent element. */
    const unsigned long elemID = surfElem[l].volElemID;

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    ViscousNormalFluxFace(&volElem[elemID], l, nInt, NPad, Wall_HeatFlux,
                          HeatFlux_Prescribed, solIntL, gradSolInt,
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(),
                          viscFluxes, viscosityInt, kOverCvInt);
  }
}

void CFEM_DG_NSSolver::WallTreatmentViscousFluxes(
                                  CConfig                  *config,
                                  const unsigned short     nFaceSimul,
                                  const unsigned short     NPad,
                                  const unsigned short     nInt,
                                  const su2double          Wall_HeatFlux,
                                  const bool               HeatFlux_Prescribed,
                                  const su2double          Wall_Temperature,
                                  const bool               Temperature_Prescribed,
                                  const CSurfaceElementFEM *surfElem,
                                  const su2double          *solIntL,
                                        su2double          *workArray,
                                        su2double          *viscFluxes,
                                        su2double          *viscosityInt,
                                        su2double          *kOverCvInt,
                                        CWallModel         *wallModel) {

  /* Loop over the simultaneously treated faces. */
  for(unsigned short l=0; l<nFaceSimul; ++l) {
    const unsigned short llNVar = l*nVar;

    /* Loop over the donors for this boundary face. */
    for(unsigned long j=0; j<surfElem[l].donorsWallFunction.size(); ++j) {

      /* Easier storage of the element ID of the donor and set the pointer
         where the solution of this element starts. Note that by construction
         the time level of the donor element is the same as the time level
         of the boundary face. */
      const unsigned long donorID    = surfElem[l].donorsWallFunction[j];
      const unsigned short timeLevel = volElem[donorID].timeLevel;
      const unsigned short nDOFsElem = volElem[donorID].nDOFsSol;
      const su2double *solDOFsElem   = VecWorkSolDOFs[timeLevel].data()
                                     + nVar*volElem[donorID].offsetDOFsSolThisTimeLevel;

      /* Determine the number of integration points for this donor and
         interpolate the solution for the corresponding exchange points. */
      const unsigned short nIntThisDonor = surfElem[l].nIntPerWallFunctionDonor[j+1]
                                         - surfElem[l].nIntPerWallFunctionDonor[j];

      blasFunctions->gemm(nIntThisDonor, nVar, nDOFsElem, surfElem[l].matWallFunctionDonor[j].data(),
                          solDOFsElem, workArray, config);

      /* Loop over the integration points for this donor element. */
      for(unsigned short i=surfElem[l].nIntPerWallFunctionDonor[j];
                         i<surfElem[l].nIntPerWallFunctionDonor[j+1]; ++i) {

        /* Easier storage of the actual integration point. */
        const unsigned short ii = surfElem[l].intPerWallFunctionDonor[i];

        /* Determine the normal and the wall velocity for this integration point. */
        const su2double *normals = surfElem[l].metricNormalsFace.data() + ii*(nDim+1);
        const su2double *gridVel = surfElem[l].gridVelocities.data() + ii*nDim;

        /* Determine the velocities and pressure in the exchange point. */
        const su2double *solInt = workArray
                                + nVar*(i-surfElem[l].nIntPerWallFunctionDonor[j]);

        su2double rhoInv = 1.0/solInt[0];
        su2double vel[]  = {0.0, 0.0, 0.0};
        for(unsigned short k=0; k<nDim; ++k) vel[k] = rhoInv*solInt[k+1];

        su2double vel2Mag = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
        su2double eInt    = rhoInv*solInt[nVar-1] - 0.5*vel2Mag;

        FluidModel->SetTDState_rhoe(solInt[0], eInt);
        const su2double Pressure = FluidModel->GetPressure();
        const su2double Temperature = FluidModel->GetTemperature();
        const su2double LaminarViscosity= FluidModel->GetLaminarViscosity();

        /* Subtract the prescribed wall velocity, i.e. grid velocity
           from the velocity in the exchange point. */
        for(unsigned short k=0; k<nDim; ++k) vel[k] -= gridVel[k];

        /* Determine the tangential velocity by subtracting the normal
           velocity component. */
        su2double velNorm = 0.0;
        for(unsigned short k=0; k<nDim; ++k) velNorm += normals[k]*vel[k];
        for(unsigned short k=0; k<nDim; ++k) vel[k]  -= normals[k]*velNorm;

        /* Determine the magnitude of the tangential velocity as well
           as its direction (unit vector). */
        su2double velTan = sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
        velTan = max(velTan,1.e-25);

        su2double dirTan[] = {0.0, 0.0, 0.0};
        for(unsigned short k=0; k<nDim; ++k) dirTan[k] = vel[k]/velTan;

        /* Compute the wall shear stress and heat flux vector using
           the wall model. */
        su2double tauWall, qWall, ViscosityWall, kOverCvWall;

        wallModel->WallShearStressAndHeatFlux(Temperature, velTan, LaminarViscosity, Pressure,
                                              Wall_HeatFlux, HeatFlux_Prescribed,
                                              Wall_Temperature, Temperature_Prescribed,
                                              FluidModel, tauWall, qWall, ViscosityWall,
                                              kOverCvWall);

        /* Compute the wall velocity in tangential direction. */
        const su2double *solWallInt = solIntL + NPad*ii + llNVar;
        su2double velWallTan = 0.0;
        for(unsigned short k=0; k<nDim; ++k)
          velWallTan += solWallInt[k+1]*dirTan[k];
        velWallTan /= solWallInt[0];

        /* Determine the position where the viscous fluxes, viscosity and
           thermal conductivity must be stored. */
        su2double *normalFlux = viscFluxes + NPad*ii + llNVar;

        const unsigned short ind = l*nInt + ii;
        viscosityInt[ind] = ViscosityWall;
        kOverCvInt[ind]   = kOverCvWall;

        /* Compute the viscous normal flux. Note that the unscaled normals
           must be used, hence the multiplication with normals[nDim]. */
        normalFlux[0] = 0.0;
        for(unsigned short k=0; k<nDim; ++k)
          normalFlux[k+1] = -normals[nDim]*tauWall*dirTan[k];
        normalFlux[nVar-1] = normals[nDim]*(qWall - tauWall*velWallTan);
      }
    }
  }
}

void CFEM_DG_NSSolver::ResidualViscousBoundaryFace(
                                      CConfig                  *config,
                                      CNumerics                *conv_numerics,
                                      const unsigned short     nFaceSimul,
                                      const unsigned short     NPad,
                                      const CSurfaceElementFEM *surfElem,
                                      const su2double          *solInt0,
                                      const su2double          *solInt1,
                                      su2double                *paramFluxes,
                                      su2double                *fluxes,
                                      su2double                *viscFluxes,
                                      const su2double          *viscosityInt,
                                      const su2double          *kOverCvInt,
                                      su2double                *resFaces,
                                      unsigned long            &indResFaces) {

  /*--- Get the required information from the standard element, which is the
        same for all faces considered. ---*/
  const unsigned short ind          = surfElem[0].indStandardElement;
  const unsigned short nInt         = standardBoundaryFacesSol[ind].GetNIntegration();
  const unsigned short nDOFs        = standardBoundaryFacesSol[ind].GetNDOFsFace();
  const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
  const su2double      ConstPenFace = standardBoundaryFacesSol[ind].GetPenaltyConstant();
  const su2double     *weights      = standardBoundaryFacesSol[ind].GetWeightsIntegration();

  /*------------------------------------------------------------------------*/
  /*--- Step 1: Compute the sum of the inviscid, viscous and penalty     ---*/
  /*---         fluxes in the integration points of the faces.           ---*/
  /*------------------------------------------------------------------------*/

  /* Store the normals and grid velocities of the faces in an array of
     pointers to be used in the call to ComputeInviscidFluxesFace. */
  const su2double* arrNorm[SIZE_ARR_NORM];
  const su2double* arrGridVel[SIZE_ARR_NORM];

  if(nFaceSimul > SIZE_ARR_NORM)
    SU2_MPI::Error("SIZE_ARR_NORM is too small. Increase it or decrease ALIGNED_BYTES_MATMUL",
                   CURRENT_FUNCTION);

  for(unsigned short l=0; l<nFaceSimul; ++l) {
    arrNorm[l]    = surfElem[l].metricNormalsFace.data();
    arrGridVel[l] = surfElem[l].gridVelocities.data();
  }

  /* General function to compute the inviscid fluxes, using an approximate
     Riemann solver in the integration points. */
  ComputeInviscidFluxesFace(config, nFaceSimul, NPad, nInt, arrNorm, arrGridVel,
                            solInt0, solInt1, fluxes, conv_numerics);

  /* Subtract the viscous fluxes from the inviscid fluxes. */
  for(unsigned short j=0; j<(NPad*nInt); ++j) fluxes[j] -= viscFluxes[j];

  /*--- Loop over the faces in this chunk to compute the penalty fluxes. ---*/
  for(unsigned short l=0; l<nFaceSimul; ++l) {

    /* Get the length scale for the adjacent element. */
    const su2double lenScale = volElem[surfElem[l].volElemID].lenScale;

    /* Call the function PenaltyTermsFluxFace to compute the actual penalty
       terms. Use the array viscFluxes as storage. */
    PenaltyTermsFluxFace(l, nInt, NPad, solInt0, solInt1, viscosityInt,
                         viscosityInt, kOverCvInt, kOverCvInt, ConstPenFace,
                         lenScale, lenScale, surfElem[l].metricNormalsFace.data(),
                         viscFluxes);
  }

  /* Add the penalty fluxes to the earlier computed fluxes. */
  for(unsigned short j=0; j<(NPad*nInt); ++j) fluxes[j] += viscFluxes[j];

  /* Multiply the fluxes with the integration weight of the corresponding
     integration point. */
  for(unsigned short i=0; i<nInt; ++i) {
    su2double *flux = fluxes + i*NPad;

    for(unsigned short j=0; j<NPad; ++j)
      flux[j] *= weights[i];
  }

  /*------------------------------------------------------------------------*/
  /*--- Step 2: Compute the contribution to the residuals from the       ---*/
  /*---         integration over the surface element of the invisid      ---*/
  /*---         fluxes, viscous fluxes and penalty terms.                ---*/
  /*------------------------------------------------------------------------*/

  /* Set the value of the offset of the residual between each of the
     faces. This depends whether or not symmetrizing terms are stored.
     Also set the location where the symmetrizing terms of the first face
     must be stored. */
  unsigned long offsetRes = nDOFs;
  if( symmetrizingTermsPresent ) offsetRes += nDOFsElem;

  unsigned long indResSym = indResFaces + nDOFs;

  /* Get the correct form of the basis functions needed for the matrix
     multiplication to compute the residual. */
  const su2double *basisFaceTrans = standardBoundaryFacesSol[ind].GetBasisFaceIntegrationTranspose();

  /* Call the general function to carry out the matrix product. Use viscFluxes
     as temporary storage for the result. */
  blasFunctions->gemm(nDOFs, NPad, nInt, basisFaceTrans, fluxes, viscFluxes, config);

  /* Loop over the number of faces in this chunk to store the residual in
     the correct locations in resFaces. */
  for(unsigned short l=0; l<nFaceSimul; ++l) {
    const unsigned short llNVar = l*nVar;

    /* Easier storage of the position in the residual array for this face
       and update the corresponding counter. */
    su2double *resFace = resFaces + indResFaces*nVar;
    indResFaces       += offsetRes;

    /* Loop over the DOFs and copy the data. */
    for(unsigned short i=0; i<nDOFs; ++i)
      for(unsigned short mm=0; mm<nVar; ++mm)
        resFace[nVar*i+mm] = viscFluxes[NPad*i+llNVar+mm];
  }

  /*------------------------------------------------------------------------*/
  /*--- Step 3: Compute and distribute the symmetrizing terms, if        ---*/
  /*---         present. Note that these terms must be distributed to    ---*/
  /*---         DOFs of the adjacent element, not only of the face.      ---*/
  /*------------------------------------------------------------------------*/

  if( symmetrizingTermsPresent ) {

    /* The symmetrizing fluxes will be multiplied by 0.5 times the theta
       parameter. Store this factor a bit easier. */
    const su2double halfTheta = 0.5*config->GetTheta_Interior_Penalty_DGFEM();

    /* Loop over the simultaneously treated faces. */
    for(unsigned short l=0; l<nFaceSimul; ++l) {

      /* Compute the symmetrizing fluxes in the nDim directions. */
      SymmetrizingFluxesFace(l, nInt, NPad, solInt0, solInt1, viscosityInt,
                             viscosityInt, kOverCvInt, kOverCvInt,
                             surfElem[l].metricNormalsFace.data(), fluxes);

      /* Transform the fluxes, such that they must be multiplied with the
         gradients w.r.t. the parametric coordinates rather than the
         Cartesian coordinates of the basis functions. Also a multiplication
         with the integration weight and the theta parameter is carried out. */
      TransformSymmetrizingFluxes(l, nInt, NPad, halfTheta, fluxes, weights,
                                  surfElem[l].metricCoorDerivFace.data(),
                                  paramFluxes);
    }

    /* Call the general function to carry out the matrix product to compute
       the residual. Use fluxes as temporary storage for the result. */
    const su2double *derBasisElemTrans = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegrationTranspose();

    blasFunctions->gemm(nDOFsElem, NPad, nInt*nDim, derBasisElemTrans, paramFluxes, fluxes, config);

    /* Loop over the faces of this chunk to store the residual in
       the correct locations in resFaces. */
    for(unsigned short l=0; l<nFaceSimul; ++l) {
      const unsigned short llNVar = l*nVar;

      /* Easier storage of the position in the residual array for this face
         and update the corresponding counter. */
      su2double *resElem = resFaces + indResSym*nVar;
      indResSym         += offsetRes;

      /* Loop over the DOFs of the element and copy the data. */
      for(unsigned short i=0; i<nDOFsElem; ++i)
        for(unsigned short mm=0; mm<nVar; ++mm)
          resElem[nVar*i+mm] = fluxes[NPad*i+llNVar+mm];
    }
  }
}
