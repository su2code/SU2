/*!
 * \file CFEM_DG_NSSolver.cpp
 * \brief Main subroutines for solving finite element Navier-Stokes flow problems
 * \author J. Alonso, E. van der Weide, T. Economon
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

CFEM_DG_NSSolver::CFEM_DG_NSSolver(CGeometry      *geometry,
                                   CConfig        *config,
                                   unsigned short iMesh)
 : CFEM_DG_EulerSolver(geometry, config, iMesh) {

  /*--- Check if the symmetrizing terms are present. ---*/
  if(fabs(config->GetTheta_Interior_Penalty_DGFEM()) > 1.e-8)
    symmetrizingTermsPresent = true;

  /*--- Get the viscous data at the farfield from config. ---*/
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  Prandtl_Lam   = config->GetPrandtl_Lam();
  Prandtl_Turb  = config->GetPrandtl_Turb();
  Tke_Inf       = config->GetTke_FreeStreamND();

  /*--- Set the SGS model in case an LES simulation is carried out ---*/
  if(config->GetKind_Solver() == FEM_LES) {

    /*--- Make a distinction between the SGS models used and set SGSModel and
          SGSModelUsed accordingly. ---*/
    switch( config->GetKind_SGS_Model() ) {

      case IMPLICIT_LES:
        SGSModel     = nullptr;
        SGSModelUsed = false;
        break;

      case SMAGORINSKY:
        SGSModel     = new CSmagorinskyModel;
        SGSModelUsed = true;
        break;

      case WALE:
        SGSModel     = new CWALEModel;
        SGSModelUsed = true;
        break;

      case VREMAN:
        SGSModel     = new CVremanModel;
        SGSModelUsed = true;
        break;

      default:
        SU2_MPI::Error(string("Unknown SGS model encountered"),
                       CURRENT_FUNCTION);
    }
  }
}

CFEM_DG_NSSolver::~CFEM_DG_NSSolver(void) {

  delete SGSModel;
}

void CFEM_DG_NSSolver::Friction_Forces(const CGeometry *geometry, const CConfig *config) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  /*--- Initialize the minimum and maximum time step. ---*/
  SU2_OMP_SINGLE
  {
    Min_Delta_Time = 1.e25;
    Max_Delta_Time = 0.0;
  }

  /*--- Determine the chunk size for the OMP loops, if supported. ---*/
#ifdef HAVE_OMP
  const size_t omp_chunk_size_elem = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif

  /*--- Check whether or not a time stepping scheme is used and store
        the CFL number a bit easier. Note that if we are using explicit
        time stepping, the regular CFL condition has been overwritten with the
        unsteady CFL condition in the config post-processing (if non-zero). ---*/
  const bool time_stepping = config->GetTime_Marching() == TIME_STEPPING;
  const su2double CFL = config->GetCFL(iMesh);

  /*--- Constant factor present in the heat flux vector, namely the ratio of
        thermal conductivity and viscosity. ---*/
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /*--- Constant ratio of the second viscosity and the viscosity itself. ---*/
  const su2double lambdaOverMu = -TWO3;

  /*--- The eigenvalues of the viscous Jacobian, scaled by the kinematic viscosity,
        are 1.0, 2.0 + lambdaOverMu and kOverCv/Mu. The last is variable due to the
        possible presence of an eddy viscosity, but the first two are constant and
        the maximum can be determined. */
  const su2double radOverNuTerm = max(1.0, 2.0+lambdaOverMu);

  /*--- Check for explicit time stepping with imposed time step. If the unsteady
        CFL is set to zero (default), it uses the defined unsteady time step,
        otherwise it computes the time step based on the provided unsteady CFL.
        Note that the regular CFL option in the config is always ignored with
        time stepping. ---*/
  if(time_stepping && (config->GetUnst_CFL() == 0.0)) {

    /*--- Loop over the owned volume elements and set the fixed dt. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size_elem)
    for(unsigned long i=0; i<nVolElemOwned; ++i)
      volElem[i].deltaTime = config->GetDelta_UnstTimeND();
  } else {

    /*--- Define the thread local variables for the minimum and maximum time step. ---*/
    su2double MaxDeltaT = 0.0, MinDeltaT = 1.e25;

    /*--- Loop over the owned volume elements. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size_elem)
    for(unsigned long l=0; l<nVolElemOwned; ++l) {

      /*--- Determine the solution in the integration points and
            convert it to primitive variables. ---*/
      ColMajorMatrix<su2double> &solInt = volElem[l].ComputeSolIntPoints();
      EntropyToPrimitiveVariables(solInt);

      /*--- Determine the gradients in the integration points if an SGS
            model is used. Otherwise set the reference to the appropriate
            memory. ---*/
      vector<ColMajorMatrix<su2double> > &gradSolInt = SGSModelUsed ?
        volElem[l].ComputeGradSolIntPoints() :
        volElem[l].standardElemFlow->workGradSolInt[omp_get_thread_num()];

      /*--- Abbreviate the number of integration points and its padded version. ---*/
      const unsigned short nInt    = volElem[l].standardElemFlow->GetNIntegration();
      const unsigned short nIntPad = volElem[l].standardElemFlow->GetNIntegrationPad();

      /*--- Check if an SGS model must be used. ---*/
      if( SGSModelUsed ) {

        /*--- Determine the length scale for the LES. This is the length scale
              of the element corrected with its polynomial degree. ---*/
        unsigned short nPoly = volElem[l].standardElemFlow->GetPolyDegree();
        if(nPoly == 0) nPoly = 1;
        const su2double lenScale = volElem[l].lenScale/nPoly;

        /*--- Make a distinction between two and three space dimensions
              in order to have the most efficient code. ---*/
        switch( nDim ) {

          case 2: {

            /*--- Two dimensional simulation. Easier storage of the metric terms. ---*/
            ColMajorMatrix<su2double> &dParDx = volElem[l].metricTermsInt[0];
            ColMajorMatrix<su2double> &dParDy = volElem[l].metricTermsInt[1];

            /*--- Loop over the integration points to compute the Cartesian
                  gradients of the velocity. ---*/
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nIntPad; ++i) {

              /*--- Easier storage of some of the primitive variables. ---*/
              const su2double u     =  solInt(i,1);
              const su2double v     =  solInt(i,2);
              const su2double V3Inv = -solInt(i,3)/solInt(i,0);

              /*--- Compute the gradients in computational space. ---*/
              const su2double dudr = -V3Inv*(gradSolInt[0](i,1) + u*gradSolInt[0](i,3));
              const su2double dvdr = -V3Inv*(gradSolInt[0](i,2) + v*gradSolInt[0](i,3));

              const su2double duds = -V3Inv*(gradSolInt[1](i,1) + u*gradSolInt[1](i,3));
              const su2double dvds = -V3Inv*(gradSolInt[1](i,2) + v*gradSolInt[1](i,3));

              /*--- Compute the true value of the metric terms in this integration point.
                    Note that the metric terms stored are scaled by the Jacobian. ---*/
              const su2double JacInv = 1.0/volElem[l].JacobiansInt(i);

              const su2double drdx = JacInv*dParDx(i,0);
              const su2double dsdx = JacInv*dParDx(i,1);

              const su2double drdy = JacInv*dParDy(i,0);
              const su2double dsdy = JacInv*dParDy(i,1);

              /*--- Compute the Cartesian gradients. ---*/
              gradSolInt[0](i,1) = dudr*drdx + duds*dsdx;
              gradSolInt[0](i,2) = dvdr*drdx + dvds*dsdx;

              gradSolInt[1](i,1) = dudr*drdy + duds*dsdy;
              gradSolInt[1](i,2) = dvdr*drdy + dvds*dsdy;
            }

            /*--- Loop over the integration points to compute the SGS
                  eddy viscosity for the 2D situation, which is stored
                  in gradSolInt[0](i,0). ---*/
            for(unsigned short i=0; i<nIntPad; ++i)
              gradSolInt[0](i,0) = SGSModel->ComputeEddyViscosity_2D(solInt(i,0),
                                                                     gradSolInt[0](i,1),
                                                                     gradSolInt[1](i,1),
                                                                     gradSolInt[0](i,2),
                                                                     gradSolInt[1](i,2),
                                                                     lenScale,
                                                                     volElem[l].wallDistanceInt(i));
            break;
          }

          case 3: {

            /*--- Three dimensional simulation. Easier storage of the metric terms. ---*/
            ColMajorMatrix<su2double> &dParDx = volElem[l].metricTermsInt[0];
            ColMajorMatrix<su2double> &dParDy = volElem[l].metricTermsInt[1];
            ColMajorMatrix<su2double> &dParDz = volElem[l].metricTermsInt[2];

            /*--- Loop over the integration points to compute the Cartesian
                  gradients of the velocity. ---*/
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nIntPad; ++i) {

              /*--- Easier storage of some of the primitive variables. ---*/
              const su2double u     =  solInt(i,1);
              const su2double v     =  solInt(i,2);
              const su2double w     =  solInt(i,3);
              const su2double V4Inv = -solInt(i,4)/solInt(i,0);

              /*--- Compute the gradients in computational space. ---*/
              const su2double dudr = -V4Inv*(gradSolInt[0](i,1) + u*gradSolInt[0](i,4));
              const su2double dvdr = -V4Inv*(gradSolInt[0](i,2) + v*gradSolInt[0](i,4));
              const su2double dwdr = -V4Inv*(gradSolInt[0](i,3) + w*gradSolInt[0](i,4));

              const su2double duds = -V4Inv*(gradSolInt[1](i,1) + u*gradSolInt[1](i,4));
              const su2double dvds = -V4Inv*(gradSolInt[1](i,2) + v*gradSolInt[1](i,4));
              const su2double dwds = -V4Inv*(gradSolInt[1](i,3) + w*gradSolInt[1](i,4));

              const su2double dudt = -V4Inv*(gradSolInt[2](i,1) + u*gradSolInt[2](i,4));
              const su2double dvdt = -V4Inv*(gradSolInt[2](i,2) + v*gradSolInt[2](i,4));
              const su2double dwdt = -V4Inv*(gradSolInt[2](i,3) + w*gradSolInt[2](i,4));

              /*--- Compute the true value of the metric terms in this integration point.
                    Note that the metric terms stored are scaled by the Jacobian. ---*/
              const su2double JacInv = 1.0/volElem[l].JacobiansInt(i);

              const su2double drdx = JacInv*dParDx(i,0);
              const su2double dsdx = JacInv*dParDx(i,1);
              const su2double dtdx = JacInv*dParDx(i,2);

              const su2double drdy = JacInv*dParDy(i,0);
              const su2double dsdy = JacInv*dParDy(i,1);
              const su2double dtdy = JacInv*dParDy(i,2);

              const su2double drdz = JacInv*dParDz(i,0);
              const su2double dsdz = JacInv*dParDz(i,1);
              const su2double dtdz = JacInv*dParDz(i,2);

              /*--- Compute the Cartesian gradients. ---*/
              gradSolInt[0](i,1) = dudr*drdx + duds*dsdx + dudt*dtdx;
              gradSolInt[0](i,2) = dvdr*drdx + dvds*dsdx + dvdt*dtdx;
              gradSolInt[0](i,3) = dwdr*drdx + dwds*dsdx + dwdt*dtdx;

              gradSolInt[1](i,1) = dudr*drdy + duds*dsdy + dudt*dtdy;
              gradSolInt[1](i,2) = dvdr*drdy + dvds*dsdy + dvdt*dtdy;
              gradSolInt[1](i,3) = dwdr*drdy + dwds*dsdy + dwdt*dtdy;

              gradSolInt[2](i,1) = dudr*drdz + duds*dsdz + dudt*dtdz;
              gradSolInt[2](i,2) = dvdr*drdz + dvds*dsdz + dvdt*dtdz;
              gradSolInt[2](i,3) = dwdr*drdz + dwds*dsdz + dwdt*dtdz;
            }

            /*--- Loop over the integration points to compute the SGS
                  eddy viscosity for the 3D situation, which is stored
                  in gradSolInt[0](i,0). ---*/
            for(unsigned short i=0; i<nIntPad; ++i)
              gradSolInt[0](i,0) = SGSModel->ComputeEddyViscosity_3D(solInt(i,0),
                                                                     gradSolInt[0](i,1),
                                                                     gradSolInt[1](i,1),
                                                                     gradSolInt[2](i,1),
                                                                     gradSolInt[0](i,2),
                                                                     gradSolInt[1](i,2),
                                                                     gradSolInt[2](i,2),
                                                                     gradSolInt[0](i,3),
                                                                     gradSolInt[1](i,3),
                                                                     gradSolInt[2](i,3),
                                                                     lenScale,
                                                                     volElem[l].wallDistanceInt(i));
            break;
          }
        }

        /*--- Set the eddy viscosity for the padded values to avoid problems. ---*/
        for(unsigned short i=nInt; i<nIntPad; ++i)
          gradSolInt[0](i,0) = gradSolInt[0](0,0);

        SU2_MPI::Error(string("SGS model not there yet in time step computation"), CURRENT_FUNCTION);
      }
      else {

        /*--- No subgrid scale model is used. Set the eddy viscosity to zero, which is
              stored in the first entry of gradSolInt. ---*/
        SU2_OMP_SIMD
        for(unsigned short i=0; i<nIntPad; ++i)
          gradSolInt[0](i,0) = 0.0;
      }

      /*--- Compute the laminar viscosities in the integration points,
            which is stored in gradSolInt[1](i,0). This loop is not
            vectorized, because of the call to the gas model. ---*/
      for(unsigned short i=0; i<nInt; ++i) {
        GetFluidModel()->SetTDState_Prho(solInt(i,nDim+1), solInt(i,0));
        gradSolInt[1](i,0) = GetFluidModel()->GetLaminarViscosity();
      }

      for(unsigned short i=nInt; i<nIntPad; ++i)
        gradSolInt[1](i,0) = gradSolInt[1](0,0);

      /*--- Make a distinction between two and three space dimensions
            in order to have the most efficient code. ---*/
      switch( nDim ) {

        case 2: {

          /*--- Two dimensional simulation. Loop over the padded integration
                points and compute the inviscid and viscous spectral radius. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nIntPad; ++i) {

            /*--- Compute the inviscid spectral radius squared, which is a
                  rather conservative estimate. ---*/
            const su2double rhoInv = 1.0/solInt(i,0);
            const su2double u      = solInt(i,1) - volElem[l].gridVelocitiesInt(i,0);
            const su2double v      = solInt(i,2) - volElem[l].gridVelocitiesInt(i,1);
            const su2double a      = sqrt(fabs(Gamma*solInt(i,3)*rhoInv));

            const su2double radx = fabs(u) + a;
            const su2double rady = fabs(v) + a;

            solInt(i,0) = radx*radx + rady*rady;

            /*--- Compute the viscous spectral radius. ---*/
            const su2double muLam        = gradSolInt[1](i,0);
            const su2double muTurb       = gradSolInt[0](i,0);
            const su2double mu           = muLam + muTurb;
            const su2double kOverCv      = muLam*factHeatFlux_Lam
                                         + muTurb*factHeatFlux_Turb;
            const su2double factHeatFlux = kOverCv/mu;

            solInt(i,1) = rhoInv*mu*max(radOverNuTerm, factHeatFlux);
          }

          break;
        }

        case 3: {

          /*--- Three dimensional simulation. Loop over the padded integration
                points and compute the inviscid and viscous spectral radius. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nIntPad; ++i) {

            /*--- Compute the inviscid spectral radius squared, which is a
                  rather conservative estimate. ---*/
            const su2double rhoInv = 1.0/solInt(i,0);
            const su2double u      = solInt(i,1) - volElem[l].gridVelocitiesInt(i,0);
            const su2double v      = solInt(i,2) - volElem[l].gridVelocitiesInt(i,1);
            const su2double w      = solInt(i,3) - volElem[l].gridVelocitiesInt(i,2);
            const su2double a      = sqrt(fabs(Gamma*solInt(i,4)*rhoInv));

            const su2double radx = fabs(u) + a;
            const su2double rady = fabs(v) + a;
            const su2double radz = fabs(w) + a;

            solInt(i,0) = radx*radx + rady*rady + radz*radz;

            /*--- Compute the viscous spectral radius. ---*/
            const su2double muLam        = gradSolInt[1](i,0);
            const su2double muTurb       = gradSolInt[0](i,0);
            const su2double mu           = muLam + muTurb;
            const su2double kOverCv      = muLam*factHeatFlux_Lam
                                         + muTurb*factHeatFlux_Turb;
            const su2double factHeatFlux = kOverCv/mu;

            solInt(i,1) = rhoInv*mu*max(radOverNuTerm, factHeatFlux);
          }

          break;
        }
      }

      /*--- Determine the maximum value of the inviscid spectral radius squared
            and the viscous spectra radisu for the integration points. ---*/
      su2double charVel2Max = 0.0, radViscMax = 0.0;
      for(unsigned short i=0; i<nInt; ++i) {
        charVel2Max = max(charVel2Max, solInt(i,0));
        radViscMax  = max(radViscMax,  solInt(i,1));
      }

      /*--- Compute the time step for the element. Note that for the spectral
            radii correction factors, which are a function of the polynomial degree
            and the element type, must be taken into account. ---*/
      const passivedouble factInv = volElem[l].standardElemFlow->GetFactorInviscidSpectralRadius();
      const passivedouble factVis = volElem[l].standardElemFlow->GetFactorViscousSpectralRadius();

      const su2double invLen = 1.0/volElem[l].lenScale;
      const su2double dtInv  = invLen*(factInv*sqrt(charVel2Max) + factVis*radViscMax*invLen);

      volElem[l].deltaTime = CFL/dtInv;

      /*--- Update the minimum and maximum value, for which the factor for
            time accurate local time stepping must be taken into account. ---*/
      const su2double dtEff = volElem[l].factTimeLevel*volElem[l].deltaTime;
      MinDeltaT = min(MinDeltaT, dtEff);
      MaxDeltaT = max(MaxDeltaT, dtEff);
    }

    /*--- Update the shared variables Min_Delta_Time and Max_Delta_Time. ---*/
    SU2_OMP_CRITICAL
    {
      Min_Delta_Time = min(Min_Delta_Time, MinDeltaT);
      Max_Delta_Time = max(Max_Delta_Time, MaxDeltaT);
    }

    /*--- Compute the max and the min dt (in parallel). Note that we only
          do this for steady calculations if the high verbosity is set, but we
          always perform the reduction for unsteady calculations where the CFL
          limit is used to set the global time step. ---*/
#ifdef HAVE_MPI
    SU2_OMP_SINGLE
    {
      if ((config->GetComm_Level() == COMM_FULL) || time_stepping) {
        su2double rbuf_time = Min_Delta_Time;
        SU2_MPI::Allreduce(&rbuf_time, &Min_Delta_Time, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());

        rbuf_time = Max_Delta_Time;
        SU2_MPI::Allreduce(&rbuf_time, &Max_Delta_Time, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      }
    }
#endif

    /*--- For explicit time stepping with an unsteady CFL imposed, use the
          minimum delta time of the entire mesh. As Min_Delta_Time is scaled to
          the time step of the largest time level, a correction must be used
          for the time level when time accurate local time stepping is used. ---*/
    if (time_stepping) {
      SU2_OMP_FOR_STAT(omp_chunk_size_elem)
      for(unsigned long l=0; l<nVolElemOwned; ++l)
        volElem[l].deltaTime = Min_Delta_Time/volElem[l].factTimeLevel;

      SU2_OMP_SINGLE
      config->SetDelta_UnstTimeND(Min_Delta_Time);
    }
  }
}

void CFEM_DG_NSSolver::ADER_DG_AliasedPredictorResidual_2D(CConfig              *config,
                                                           CVolumeElementFEM_DG *elem,
                                                           const su2double      *sol,
                                                           const unsigned short nSimul,
                                                           const unsigned short NPad,
                                                           su2double            *res,
                                                           su2double            *work) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::ADER_DG_AliasedPredictorResidual_3D(CConfig              *config,
                                                           CVolumeElementFEM_DG *elem,
                                                           const su2double      *sol,
                                                           const unsigned short nSimul,
                                                           const unsigned short NPad,
                                                           su2double            *res,
                                                           su2double            *work) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE  
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::ADER_DG_NonAliasedPredictorResidual_2D(CConfig              *config,
                                                              CVolumeElementFEM_DG *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::ADER_DG_NonAliasedPredictorResidual_3D(CConfig              *config,
                                                              CVolumeElementFEM_DG *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::Shock_Capturing_DG(CConfig             *config,
                                          const unsigned long elemBeg,
                                          const unsigned long elemEnd) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}
void CFEM_DG_NSSolver::Shock_Capturing_DG_Persson(const unsigned long elemBeg,
                                                  const unsigned long elemEnd,
                                                  su2double           *workArray) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::Volume_Residual(CConfig             *config,
                                       const unsigned long elemBeg,
                                       const unsigned long elemEnd) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::ResidualFaces(CConfig             *config,
                                     const unsigned long indFaceBeg,
                                     const unsigned long indFaceEnd,
                                     CNumerics           *numerics) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::ViscousNormalFluxFace(const CVolumeElementFEM_DG *adjVolElem,
                                             const unsigned short       indFaceChunk,
                                             const unsigned short       nInt,
                                             const unsigned short       NPad,
                                             const su2double            Wall_HeatFlux,
                                             const bool                 HeatFlux_Prescribed,
                                             const su2double            *solInt,
                                             const su2double            *gradSolInt,
                                             const su2double            *metricCoorDerivFace,
                                             const su2double            *metricNormalsFace,
                                             const su2double            *wallDistanceInt,
                                                   su2double            *viscNormFluxes,
                                                   su2double            *viscosityInt,
                                                   su2double            *kOverCvInt) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
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
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
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
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
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
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
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
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::TransformSymmetrizingFluxes(const unsigned short indFaceChunk,
                                                   const unsigned short nInt,
                                                   const unsigned short NPad,
                                                   const su2double      halfTheta,
                                                   const su2double      *symmFluxes,
                                                   const su2double      *weights,
                                                   const su2double      *metricCoorFace,
                                                         su2double      *paramFluxes) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::BC_Euler_Wall(CConfig                  *config,
                                     const unsigned long      surfElemBeg,
                                     const unsigned long      surfElemEnd,
                                     const CSurfaceElementFEM *surfElem,
                                     su2double                *resFaces,
                                     CNumerics                *conv_numerics,
                                     su2double                *workArray){
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::BC_Far_Field(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    su2double                *workArray){
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::BC_Sym_Plane(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    su2double                *workArray){
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::BC_Supersonic_Outlet(CConfig                  *config,
                                            const unsigned long      surfElemBeg,
                                            const unsigned long      surfElemEnd,
                                            const CSurfaceElementFEM *surfElem,
                                            su2double                *resFaces,
                                            CNumerics                *conv_numerics,
                                            su2double                *workArray){
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::BC_Inlet(CConfig                  *config,
                                const unsigned long      surfElemBeg,
                                const unsigned long      surfElemEnd,
                                const CSurfaceElementFEM *surfElem,
                                su2double                *resFaces,
                                CNumerics                *conv_numerics,
                                unsigned short           val_marker,
                                su2double                *workArray) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::BC_Outlet(CConfig                  *config,
                                 const unsigned long      surfElemBeg,
                                 const unsigned long      surfElemEnd,
                                 const CSurfaceElementFEM *surfElem,
                                 su2double                *resFaces,
                                 CNumerics                *conv_numerics,
                                 unsigned short           val_marker,
                                 su2double                *workArray) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::BC_HeatFlux_Wall(CConfig                  *config,
                                        const unsigned long      surfElemBeg,
                                        const unsigned long      surfElemEnd,
                                        const CSurfaceElementFEM *surfElem,
                                        su2double                *resFaces,
                                        CNumerics                *conv_numerics,
                                        unsigned short           val_marker,
                                        su2double                *workArray) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::BC_Isothermal_Wall(CConfig                  *config,
                                          const unsigned long      surfElemBeg,
                                          const unsigned long      surfElemEnd,
                                          const CSurfaceElementFEM *surfElem,
                                          su2double                *resFaces,
                                          CNumerics                *conv_numerics,
                                          unsigned short           val_marker,
                                          su2double                *workArray) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::BC_Riemann(CConfig                  *config,
                                  const unsigned long      surfElemBeg,
                                  const unsigned long      surfElemEnd,
                                  const CSurfaceElementFEM *surfElem,
                                  su2double                *resFaces,
                                  CNumerics                *conv_numerics,
                                  unsigned short           val_marker,
                                  su2double                *workArray) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_NSSolver::BC_Custom(CConfig                  *config,
                                 const unsigned long      surfElemBeg,
                                 const unsigned long      surfElemEnd,
                                 const CSurfaceElementFEM *surfElem,
                                 su2double                *resFaces,
                                 CNumerics                *conv_numerics,
                                 su2double                *workArray) {
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
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
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
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
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
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
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
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
  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}
