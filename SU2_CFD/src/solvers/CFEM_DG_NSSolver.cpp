/*!
 * \file CFEM_DG_NSSolver.cpp
 * \brief Main subroutines for solving finite element Navier-Stokes flow problems
 * \author J. Alonso, E. van der Weide, T. Economon
 * \version 7.1.1 "Blackbird"
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  /*--- Initialize the minimum and maximum time step. ---*/
  SU2_OMP_SINGLE
  {
    Min_Delta_Time = 1.e25;
    Max_Delta_Time = 0.0;
  }
  END_SU2_OMP_SINGLE

  /*--- Determine the chunk size for the OMP loops, if supported. ---*/
#ifdef HAVE_OMP
  const size_t omp_chunk_size_elem = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif

  /*--- Check whether or not a time stepping scheme is used and store
        the CFL number a bit easier. Note that if we are using explicit
        time stepping, the regular CFL condition has been overwritten with the
        unsteady CFL condition in the config post-processing (if non-zero). ---*/
  const bool time_stepping = config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING;
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
    END_SU2_OMP_FOR

  } else {

    /*--- Define the thread local variables for the minimum and maximum time step. ---*/
    su2double MaxDeltaT = 0.0, MinDeltaT = 1.e25;

    /*--- Loop over the owned volume elements. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size_elem)
    for(unsigned long l=0; l<nVolElemOwned; ++l) {

      /*--- Determine the solution in the integration points and
            convert it to primitive variables. ---*/
      ColMajorMatrix<su2double> &solInt = volElem[l].ComputeSolIntPoints(volElem[l].solDOFs);
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
            for(unsigned short i=0; i<nInt; ++i)
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
            for(unsigned short i=0; i<nInt; ++i)
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
    END_SU2_OMP_FOR

    /*--- Update the shared variables Min_Delta_Time and Max_Delta_Time. ---*/
    SU2_OMP_CRITICAL
    {
      Min_Delta_Time = min(Min_Delta_Time, MinDeltaT);
      Max_Delta_Time = max(Max_Delta_Time, MaxDeltaT);
    }
    END_SU2_OMP_CRITICAL

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
    END_SU2_OMP_SINGLE
#endif

    /*--- For explicit time stepping with an unsteady CFL imposed, use the
          minimum delta time of the entire mesh. As Min_Delta_Time is scaled to
          the time step of the largest time level, a correction must be used
          for the time level when time accurate local time stepping is used. ---*/
    if (time_stepping) {
      SU2_OMP_FOR_STAT(omp_chunk_size_elem)
      for(unsigned long l=0; l<nVolElemOwned; ++l)
        volElem[l].deltaTime = Min_Delta_Time/volElem[l].factTimeLevel;
      END_SU2_OMP_FOR

      SU2_OMP_SINGLE
      config->SetDelta_UnstTimeND(Min_Delta_Time);
      END_SU2_OMP_SINGLE
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE  
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::Shock_Capturing_DG(CConfig             *config,
                                          const unsigned long elemBeg,
                                          const unsigned long elemEnd) {

  /*--- Run shock capturing algorithm ---*/
  switch( config->GetKind_FEM_DG_Shock() ) {
    case NONE:
      break;
    case PERSSON:
      Shock_Capturing_DG_Persson(elemBeg, elemEnd);
      break;
  }
}

void CFEM_DG_NSSolver::Shock_Capturing_DG_Persson(const unsigned long elemBeg,
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::Volume_Residual(CConfig             *config,
                                       const unsigned long elemBeg,
                                       const unsigned long elemEnd) {

  /*--- Abbreviation of 1/(Gamma-1) and compute the constant
        factors present in the heat flux vector. ---*/
  const su2double ovgm1             = 1.0/Gamma_Minus_One;
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /*--- Determine the chunk size for the OpenMP parallelization. ---*/
#ifdef HAVE_OMP
    const unsigned long nElem = elemEnd - elemBeg;
    const size_t omp_chunk_size = computeStaticChunkSize(nElem, omp_get_num_threads(), 64);
#endif

  /*--- Loop over the given element range. ---*/
  SU2_OMP_FOR_DYN(omp_chunk_size)
  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /*--- Abbreviate the number of integration points, its padded version
          and the integration weights. ---*/
    const unsigned short nInt    = volElem[l].standardElemFlow->GetNIntegration();
    const unsigned short nIntPad = volElem[l].standardElemFlow->GetNIntegrationPad();
    const passivedouble *weights = volElem[l].standardElemFlow->GetIntegrationWeights();

    /*--- Determine the length scale for the LES. This is the length scale
          of the element corrected with its polynomial degree. ---*/
    unsigned short nPoly = volElem[l].standardElemFlow->GetPolyDegree();
    if(nPoly == 0) nPoly = 1;
    const su2double lenScale = volElem[l].lenScale/nPoly;

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Determine the primitive solution and the gradients of    ---*/
    /*---         the entropy variables in the integration points.         ---*/
    /*------------------------------------------------------------------------*/

    ColMajorMatrix<su2double>          &solInt     = volElem[l].ComputeSolIntPoints(volElem[l].solDOFsWork);
    vector<ColMajorMatrix<su2double> > &gradSolInt = volElem[l].ComputeGradSolIntPoints();

    EntropyToPrimitiveVariables(solInt);

    /*--- Compute the transformation matrix between conservative and
          entropy variables in the integration poins. ---*/
    VolumeTransformationMatrix(&volElem[l], solInt);

    /*--- Make a distinction between two and three space dimensions
          in order to have the most efficient code. ---*/
    switch( nDim ) {
      case 2: {

        /*--- Two dimensional simulation. Easier storage of the metric terms. ---*/
        ColMajorMatrix<su2double> &dParDx = volElem[l].metricTermsInt[0];
        ColMajorMatrix<su2double> &dParDy = volElem[l].metricTermsInt[1];

        /*--------------------------------------------------------------------*/
        /*--- Step 2: Convert the gradients of the entropy variables to    ---*/
        /*---         the gradients of velocity and internal energy.       ---*/
        /*--------------------------------------------------------------------*/

        /*--- Loop over the padded number of integration points. ---*/
        SU2_OMP_SIMD_IF_NOT_AD
        for(unsigned short i=0; i<nIntPad; ++i) {

          /*--- Easier storage of some of the primitive variables. ---*/
          const su2double u     =  solInt(i,1);
          const su2double v     =  solInt(i,2);
          const su2double V4Inv = -solInt(i,3)/solInt(i,0);

          /*--- Compute the velocity gradients in computational space. ---*/
          const su2double dudr = -V4Inv*(gradSolInt[0](i,1) + u*gradSolInt[0](i,3));
          const su2double dvdr = -V4Inv*(gradSolInt[0](i,2) + v*gradSolInt[0](i,3));

          const su2double duds = -V4Inv*(gradSolInt[1](i,1) + u*gradSolInt[1](i,3));
          const su2double dvds = -V4Inv*(gradSolInt[1](i,2) + v*gradSolInt[1](i,3));

          /*--- Compute the gradient of the internal energy in computational space. ---*/
          const su2double fact = V4Inv*V4Inv*ovgm1;
          const su2double dedr = fact*gradSolInt[0](i,3);
          const su2double deds = fact*gradSolInt[1](i,3);

          /*--- Compute the true value of the metric terms in this integration point.
                Note that the metric terms stored are scaled by the Jacobian. ---*/
          const su2double JacInv = 1.0/volElem[l].JacobiansInt(i);

          const su2double drdx = JacInv*dParDx(i,0);
          const su2double dsdx = JacInv*dParDx(i,1);

          const su2double drdy = JacInv*dParDy(i,0);
          const su2double dsdy = JacInv*dParDy(i,1);

          /*--- Compute the Cartesian gradients of the velocities and
                internal energy. ---*/
          gradSolInt[0](i,1) = dudr*drdx + duds*dsdx;
          gradSolInt[0](i,2) = dvdr*drdx + dvds*dsdx;
          gradSolInt[0](i,3) = dedr*drdx + deds*dsdx;

          gradSolInt[1](i,1) = dudr*drdy + duds*dsdy;
          gradSolInt[1](i,2) = dvdr*drdy + dvds*dsdy;
          gradSolInt[1](i,3) = dedr*drdy + deds*dsdy;
        }

        /*--------------------------------------------------------------------*/
        /*--- Step 3: Compute the laminar viscosity and, if needed, the    ---*/
        /*---         eddy viscosity of the subgrid scale model.           ---*/
        /*---         These are stored in gradSolInt[1](i,0) and           ---*/
        /*---         gradSolInt[0](i,0), respectively.                    ---*/
        /*--------------------------------------------------------------------*/

        /*--- Loop over the number of integration points. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /*--- Compute the laminar viscosity from the gas model. ---*/
          GetFluidModel()->SetTDState_Prho(solInt(i,3), solInt(i,0));
          gradSolInt[1](i,0) = GetFluidModel()->GetLaminarViscosity();

          /*--- Compute the eddy viscosity, if needed. 
                Otherwise set it to zero. ---*/
          if( SGSModelUsed )
            gradSolInt[0](i,0) = SGSModel->ComputeEddyViscosity_2D(solInt(i,0),
                                                                   gradSolInt[0](i,1),
                                                                   gradSolInt[1](i,1),
                                                                   gradSolInt[0](i,2),
                                                                   gradSolInt[1](i,2),
                                                                   lenScale,
                                                                   volElem[l].wallDistanceInt(i));
          else
            gradSolInt[0](i,0) = 0.0;
        }

        /*--- Set the padded values to avoid problems. ---*/
        for(unsigned short i=nInt; i<nIntPad; ++i) {
          gradSolInt[0](i,0) = gradSolInt[0](0,0);
          gradSolInt[1](i,0) = gradSolInt[1](0,0);
        }

        /*--------------------------------------------------------------------*/
        /*--- Step 4: Compute the total fluxes (inviscid fluxes minus the  ---*/
        /*---         viscous fluxes), multiplied by minus the integration ---*/
        /*---         weight, in the integration points.                   ---*/
        /*--------------------------------------------------------------------*/

        /*--- Loop over the padded number of integration points. ---*/
        SU2_OMP_SIMD_IF_NOT_AD
        for(unsigned short i=0; i<nIntPad; ++i) {

          /*--- Easier storage of the primitive variables. ---*/
          const su2double rho = solInt(i,0);
          const su2double u   = solInt(i,1);
          const su2double v   = solInt(i,2);
          const su2double p   = solInt(i,3);

          /*--- Easier storage of the gradients. ---*/
          const su2double dudx = gradSolInt[0](i,1);
          const su2double dvdx = gradSolInt[0](i,2);
          const su2double dedx = gradSolInt[0](i,3);

          const su2double dudy = gradSolInt[1](i,1);
          const su2double dvdy = gradSolInt[1](i,2);
          const su2double dedy = gradSolInt[1](i,3);

          /*--- Compute the total energy per unit volume
                and the relative velocities. ---*/
          const su2double rE   = ovgm1*p + 0.5*rho*(u*u + v*v);
          const su2double uRel = u - volElem[l].gridVelocitiesInt(i,0);
          const su2double vRel = v - volElem[l].gridVelocitiesInt(i,1);

          /*--- Compute the total viscosity and total factor for the heat flux. ---*/
          const su2double Viscosity = gradSolInt[1](i,0) + gradSolInt[0](i,0);
          const su2double kOverCv   = gradSolInt[1](i,0)*factHeatFlux_Lam
                                    + gradSolInt[0](i,0)*factHeatFlux_Turb;

          /*--- Set the value of the second viscosity and compute the divergence
                term in the viscous normal stresses. ---*/
          const su2double lambda     = -TWO3*Viscosity;
          const su2double lamDivTerm =  lambda*(dudx + dvdy);

          /*--- Compute the viscous stress tensor and minus the heatflux vector. ---*/
          const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
          const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
          const su2double tauxy = Viscosity*(dudy + dvdx);

          const su2double qx = kOverCv*dedx;
          const su2double qy = kOverCv*dedy;

          /*--- Compute the viscous normal stress minus the pressure. ---*/
          const su2double tauxxMP = tauxx - p;
          const su2double tauyyMP = tauyy - p;

          /*--- Compute the metric terms multiplied by minus the integration weight.
                The minus sign comes from the integration by parts in the weak
                formulation. ---*/
          const su2double drdx = -weights[i]*dParDx(i,0);
          const su2double dsdx = -weights[i]*dParDx(i,1);

          const su2double drdy = -weights[i]*dParDy(i,0);
          const su2double dsdy = -weights[i]*dParDy(i,1);

          /*--- Fluxes in r-direction, which are stored in gradSolInt[0]. */
          const su2double Ur = uRel*drdx + vRel*drdy;

          gradSolInt[0](i,0) = Ur*rho;
          gradSolInt[0](i,1) = Ur*rho*u - tauxxMP*drdx - tauxy*drdy;
          gradSolInt[0](i,2) = Ur*rho*v - tauxy*drdx - tauyyMP*drdy;
          gradSolInt[0](i,3) = Ur*rE - (u*tauxxMP + v*tauxy + qx)*drdx
                                     - (u*tauxy + v*tauyyMP + qy)*drdy;

          /*--- Fluxes in s-direction, which are stored in gradSolInt[1]. */
          const su2double Us = uRel*dsdx + vRel*dsdy;

          gradSolInt[1](i,0) = Us*rho;
          gradSolInt[1](i,1) = Us*rho*u - tauxxMP*dsdx - tauxy*dsdy;
          gradSolInt[1](i,2) = Us*rho*v - tauxy*dsdx - tauyyMP*dsdy;
          gradSolInt[1](i,3) = Us*rE - (u*tauxxMP + v*tauxy + qx)*dsdx
                                     - (u*tauxy + v*tauyyMP + qy)*dsdy;
        }

        break;
      }

      case 3: {

        /*--- Three dimensional simulation. Easier storage of the metric terms. ---*/
        ColMajorMatrix<su2double> &dParDx = volElem[l].metricTermsInt[0];
        ColMajorMatrix<su2double> &dParDy = volElem[l].metricTermsInt[1];
        ColMajorMatrix<su2double> &dParDz = volElem[l].metricTermsInt[2];

        /*--------------------------------------------------------------------*/
        /*--- Step 2: Convert the gradients of the entropy variables to    ---*/
        /*---         the gradients of velocity and internal energy.       ---*/
        /*--------------------------------------------------------------------*/

        /*--- Loop over the padded number of integration points. ---*/
        SU2_OMP_SIMD_IF_NOT_AD
        for(unsigned short i=0; i<nIntPad; ++i) {

          /*--- Easier storage of some of the primitive variables. ---*/
          const su2double u     =  solInt(i,1);
          const su2double v     =  solInt(i,2);
          const su2double w     =  solInt(i,3);
          const su2double V4Inv = -solInt(i,4)/solInt(i,0);

          /*--- Compute the velocity gradients in computational space. ---*/
          const su2double dudr = -V4Inv*(gradSolInt[0](i,1) + u*gradSolInt[0](i,4));
          const su2double dvdr = -V4Inv*(gradSolInt[0](i,2) + v*gradSolInt[0](i,4));
          const su2double dwdr = -V4Inv*(gradSolInt[0](i,3) + w*gradSolInt[0](i,4));

          const su2double duds = -V4Inv*(gradSolInt[1](i,1) + u*gradSolInt[1](i,4));
          const su2double dvds = -V4Inv*(gradSolInt[1](i,2) + v*gradSolInt[1](i,4));
          const su2double dwds = -V4Inv*(gradSolInt[1](i,3) + w*gradSolInt[1](i,4));

          const su2double dudt = -V4Inv*(gradSolInt[2](i,1) + u*gradSolInt[2](i,4));
          const su2double dvdt = -V4Inv*(gradSolInt[2](i,2) + v*gradSolInt[2](i,4));
          const su2double dwdt = -V4Inv*(gradSolInt[2](i,3) + w*gradSolInt[2](i,4));

          /*--- Compute the gradient of the internal energy in computational space. ---*/
          const su2double fact = V4Inv*V4Inv*ovgm1;
          const su2double dedr = fact*gradSolInt[0](i,4);
          const su2double deds = fact*gradSolInt[1](i,4);
          const su2double dedt = fact*gradSolInt[2](i,4);

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

          /*--- Compute the Cartesian gradients of the velocities and
                internal energy. ---*/
          gradSolInt[0](i,1) = dudr*drdx + duds*dsdx + dudt*dtdx;
          gradSolInt[0](i,2) = dvdr*drdx + dvds*dsdx + dvdt*dtdx;
          gradSolInt[0](i,3) = dwdr*drdx + dwds*dsdx + dwdt*dtdx;
          gradSolInt[0](i,4) = dedr*drdx + deds*dsdx + dedt*dtdx;

          gradSolInt[1](i,1) = dudr*drdy + duds*dsdy + dudt*dtdy;
          gradSolInt[1](i,2) = dvdr*drdy + dvds*dsdy + dvdt*dtdy;
          gradSolInt[1](i,3) = dwdr*drdy + dwds*dsdy + dwdt*dtdy;
          gradSolInt[1](i,4) = dedr*drdy + deds*dsdy + dedt*dtdy;

          gradSolInt[2](i,1) = dudr*drdz + duds*dsdz + dudt*dtdz;
          gradSolInt[2](i,2) = dvdr*drdz + dvds*dsdz + dvdt*dtdz;
          gradSolInt[2](i,3) = dwdr*drdz + dwds*dsdz + dwdt*dtdz;
          gradSolInt[2](i,4) = dedr*drdz + deds*dsdz + dedt*dtdz;
        }

        /*--------------------------------------------------------------------*/
        /*--- Step 3: Compute the laminar viscosity and, if needed, the    ---*/
        /*---         eddy viscosity of the subgrid scale model.           ---*/
        /*---         These are stored in gradSolInt[1](i,0) and           ---*/
        /*---         gradSolInt[0](i,0), respectively.                    ---*/
        /*--------------------------------------------------------------------*/

        /*--- Loop over the number of integration points. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /*--- Compute the laminar viscosity from the gas model. ---*/
          GetFluidModel()->SetTDState_Prho(solInt(i,4), solInt(i,0));
          gradSolInt[1](i,0) = GetFluidModel()->GetLaminarViscosity();

          /*--- Compute the eddy viscosity, if needed. 
                Otherwise set it to zero. ---*/
          if( SGSModelUsed )
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
          else
            gradSolInt[0](i,0) = 0.0;
        }

        /*--- Set the padded values to avoid problems. ---*/
        for(unsigned short i=nInt; i<nIntPad; ++i) {
          gradSolInt[0](i,0) = gradSolInt[0](0,0);
          gradSolInt[1](i,0) = gradSolInt[1](0,0);
        }

        /*--------------------------------------------------------------------*/
        /*--- Step 4: Compute the total fluxes (inviscid fluxes minus the  ---*/
        /*---         viscous fluxes), multiplied by minus the integration ---*/
        /*---         weight, in the integration points.                   ---*/
        /*--------------------------------------------------------------------*/

        /*--- Loop over the padded number of integration points. ---*/
        SU2_OMP_SIMD_IF_NOT_AD
        for(unsigned short i=0; i<nIntPad; ++i) {

          /*--- Easier storage of the primitive variables. ---*/
          const su2double rho = solInt(i,0);
          const su2double u   = solInt(i,1);
          const su2double v   = solInt(i,2);
          const su2double w   = solInt(i,3);
          const su2double p   = solInt(i,4);

          /*--- Easier storage of the gradients. ---*/
          const su2double dudx = gradSolInt[0](i,1);
          const su2double dvdx = gradSolInt[0](i,2);
          const su2double dwdx = gradSolInt[0](i,3);
          const su2double dedx = gradSolInt[0](i,4);

          const su2double dudy = gradSolInt[1](i,1);
          const su2double dvdy = gradSolInt[1](i,2);
          const su2double dwdy = gradSolInt[1](i,3);
          const su2double dedy = gradSolInt[1](i,4);

          const su2double dudz = gradSolInt[2](i,1);
          const su2double dvdz = gradSolInt[2](i,2);
          const su2double dwdz = gradSolInt[2](i,3);
          const su2double dedz = gradSolInt[2](i,4);

          /*--- Compute the total energy per unit volume
                and the relative velocities. ---*/
          const su2double rE   = ovgm1*p + 0.5*rho*(u*u + v*v + w*w);
          const su2double uRel = u - volElem[l].gridVelocitiesInt(i,0);
          const su2double vRel = v - volElem[l].gridVelocitiesInt(i,1);
          const su2double wRel = w - volElem[l].gridVelocitiesInt(i,2);

          /*--- Compute the total viscosity and total factor for the heat flux. ---*/
          const su2double Viscosity = gradSolInt[1](i,0) + gradSolInt[0](i,0);
          const su2double kOverCv   = gradSolInt[1](i,0)*factHeatFlux_Lam
                                    + gradSolInt[0](i,0)*factHeatFlux_Turb;

          /*--- Set the value of the second viscosity and compute the divergence
                term in the viscous normal stresses. ---*/
          const su2double lambda     = -TWO3*Viscosity;
          const su2double lamDivTerm =  lambda*(dudx + dvdy + dwdz);

          /*--- Compute the viscous stress tensor and minus the heatflux vector. ---*/
          const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
          const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
          const su2double tauzz = 2.0*Viscosity*dwdz + lamDivTerm;

          const su2double tauxy = Viscosity*(dudy + dvdx);
          const su2double tauxz = Viscosity*(dudz + dwdx);
          const su2double tauyz = Viscosity*(dvdz + dwdy);

          const su2double qx = kOverCv*dedx;
          const su2double qy = kOverCv*dedy;
          const su2double qz = kOverCv*dedz;

          /*--- Compute the viscous normal stress minus the pressure. ---*/
          const su2double tauxxMP = tauxx - p;
          const su2double tauyyMP = tauyy - p;
          const su2double tauzzMP = tauzz - p;

          /*--- Compute the metric terms multiplied by minus the integration weight.
                The minus sign comes from the integration by parts in the weak
                formulation. ---*/
          const su2double drdx = -weights[i]*dParDx(i,0);
          const su2double dsdx = -weights[i]*dParDx(i,1);
          const su2double dtdx = -weights[i]*dParDx(i,2);

          const su2double drdy = -weights[i]*dParDy(i,0);
          const su2double dsdy = -weights[i]*dParDy(i,1);
          const su2double dtdy = -weights[i]*dParDy(i,2);

          const su2double drdz = -weights[i]*dParDz(i,0);
          const su2double dsdz = -weights[i]*dParDz(i,1);
          const su2double dtdz = -weights[i]*dParDz(i,2);

          /*--- Fluxes in r-direction, which are stored in gradSolInt[0]. */
          const su2double Ur = uRel*drdx + vRel*drdy + wRel*drdz;

          gradSolInt[0](i,0) = Ur*rho;
          gradSolInt[0](i,1) = Ur*rho*u - tauxxMP*drdx - tauxy*drdy - tauxz*drdz;
          gradSolInt[0](i,2) = Ur*rho*v - tauxy*drdx - tauyyMP*drdy - tauyz*drdz;
          gradSolInt[0](i,3) = Ur*rho*w - tauxz*drdx - tauyz*drdy - tauzzMP*drdz;
          gradSolInt[0](i,4) = Ur*rE - (u*tauxxMP + v*tauxy + w*tauxz + qx)*drdx
                                     - (u*tauxy + v*tauyyMP + w*tauyz + qy)*drdy
                                     - (u*tauxz + v*tauyz + w*tauzzMP + qz)*drdz;

          /*--- Fluxes in s-direction, which are stored in gradSolInt[1]. */
          const su2double Us = uRel*dsdx + vRel*dsdy + wRel*dsdz;

          gradSolInt[1](i,0) = Us*rho;
          gradSolInt[1](i,1) = Us*rho*u - tauxxMP*dsdx - tauxy*dsdy - tauxz*dsdz;
          gradSolInt[1](i,2) = Us*rho*v - tauxy*dsdx - tauyyMP*dsdy - tauyz*dsdz;
          gradSolInt[1](i,3) = Us*rho*w - tauxz*dsdx - tauyz*dsdy - tauzzMP*dsdz;
          gradSolInt[1](i,4) = Us*rE - (u*tauxxMP + v*tauxy + w*tauxz + qx)*dsdx
                                     - (u*tauxy + v*tauyyMP + w*tauyz + qy)*dsdy
                                     - (u*tauxz + v*tauyz + w*tauzzMP + qz)*dsdz;

          /*--- Fluxes in t-direction, which are stored in gradSolInt[2]. */
          const su2double Ut = uRel*dtdx + vRel*dtdy + wRel*dtdz;

          gradSolInt[2](i,0) = Ut*rho;
          gradSolInt[2](i,1) = Ut*rho*u - tauxxMP*dtdx - tauxy*dtdy - tauxz*dtdz;
          gradSolInt[2](i,2) = Ut*rho*v - tauxy*dtdx - tauyyMP*dtdy - tauyz*dtdz;
          gradSolInt[2](i,3) = Ut*rho*w - tauxz*dtdx - tauyz*dtdy - tauzzMP*dtdz;
          gradSolInt[2](i,4) = Ut*rE - (u*tauxxMP + v*tauxy + w*tauxz + qx)*dtdx
                                     - (u*tauxy + v*tauyyMP + w*tauyz + qy)*dtdy
                                     - (u*tauxz + v*tauyz + w*tauzzMP + qz)*dtdz;
        }

        break;
      }
    }

    /*--- Compute the volume source terms, if needed. ---*/
    const bool addSourceTerms = VolumeSourceTerms(config, &volElem[l], solInt);

    /*------------------------------------------------------------------------*/
    /*--- Step 5: Compute the contribution to the residuals from the       ---*/
    /*---         integration over the volume element.                     ---*/
    /*------------------------------------------------------------------------*/

    /*--- Initialize the residual to zero. ---*/
    volElem[l].resDOFs.setConstant(0.0);

    /*--- Add the contribution from the fluxes,
          which are stored in gradSolInt. ---*/
    volElem[l].ResidualGradientBasisFunctions(gradSolInt);

    /*--- Add the contribution from the source terms, if needed.
          The source terms are stored in solInt. ---*/
    if( addSourceTerms ) volElem[l].ResidualBasisFunctions(solInt);
  }
  END_SU2_OMP_FOR
}

void CFEM_DG_NSSolver::ResidualFaces(CConfig             *config,
                                     const unsigned long indFaceBeg,
                                     const unsigned long indFaceEnd,
                                     CNumerics           *numerics) {
 
  /*--- Determine the chunk size for the OpenMP parallelization. ---*/
#ifdef HAVE_OMP
    const unsigned long nFaces = indFaceEnd - indFaceBeg;
    const size_t omp_chunk_size = computeStaticChunkSize(nFaces, omp_get_num_threads(), 64);
#endif

  /*--- Loop over the given face range. ---*/
  SU2_OMP_FOR_DYN(omp_chunk_size)
  for(unsigned long l=indFaceBeg; l<indFaceEnd; ++l) {

    /*--- Abbreviate the number of padded integration points
          and the integration weights. ---*/
    const unsigned short nIntPad = matchingInternalFaces[l].standardElemFlow->GetNIntegrationPad();
    const passivedouble *weights = matchingInternalFaces[l].standardElemFlow->GetIntegrationWeights();

    /*----------------------------------------------------------------------*/
    /*--- Step 1: Compute the sum of the inviscid, viscous and penalty   ---*/
    /*---         fluxes in the integration points of the faces.         ---*/
    /*----------------------------------------------------------------------*/

    /*--- Compute the primitive variables of the left and right state
          in the integration point of the face. ---*/
    ColMajorMatrix<su2double> &solIntLeft  = matchingInternalFaces[l].ComputeSolSide0IntPoints(volElem);
    ColMajorMatrix<su2double> &solIntRight = matchingInternalFaces[l].ComputeSolSide1IntPoints(volElem);

    EntropyToPrimitiveVariables(solIntLeft);
    EntropyToPrimitiveVariables(solIntRight);

    /*--- Compute the gradients of the entropy variables of the left
          and right state in the integration point of the face. ---*/
    vector<ColMajorMatrix<su2double> > &gradSolIntLeft  = matchingInternalFaces[l].ComputeGradSolSide0IntPoints(volElem);
    vector<ColMajorMatrix<su2double> > &gradSolIntRight = matchingInternalFaces[l].ComputeGradSolSide1IntPoints(volElem);

    /*--- The gradients computed above are the gradients w.r.t. the
          parametric coordinates. Convert them to gradients w.r.t.
          Cartesian coordinates and average them afterwards. Make
          a distinction between two and three space dimensions
          in order to have the most efficient code. ---*/
    switch( nDim ) {
      case 2: {

        /*--- Two dimensional simulation. Easier storage of the metric terms. ---*/
        ColMajorMatrix<su2double> &dParDx0 =  matchingInternalFaces[l].metricCoorDerivFace0[0];
        ColMajorMatrix<su2double> &dParDy0 =  matchingInternalFaces[l].metricCoorDerivFace0[1];

        ColMajorMatrix<su2double> &dParDx1 =  matchingInternalFaces[l].metricCoorDerivFace1[0];
        ColMajorMatrix<su2double> &dParDy1 =  matchingInternalFaces[l].metricCoorDerivFace1[1];

        /*--- Loop over the number of variables and padded number of integration points. ---*/
        for(unsigned short l=0; l<nVar; ++l) {
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nIntPad; ++i) {

            /*--- Left state. Store the gradients w.r.t. the parametric coordinates. ---*/
            su2double dvdr = gradSolIntLeft[0](i,l);
            su2double dvds = gradSolIntLeft[1](i,l);

            /*--- Compute the Cartesian gradients of the left state. ---*/
            gradSolIntLeft[0](i,l) = dvdr*dParDx0(i,0) + dvds*dParDx0(i,1);
            gradSolIntLeft[1](i,l) = dvdr*dParDy0(i,0) + dvds*dParDy0(i,1);

            /*--- Right state. Store the gradients w.r.t. the parametric coordinates. ---*/
            dvdr = gradSolIntRight[0](i,l);
            dvds = gradSolIntRight[1](i,l);

            /*--- Compute the Cartesian gradients of the right state. ---*/
            gradSolIntRight[0](i,l) = dvdr*dParDx1(i,0) + dvds*dParDx1(i,1);
            gradSolIntRight[1](i,l) = dvdr*dParDy1(i,0) + dvds*dParDy1(i,1);

            /*--- Store the averaged gradients in gradSolIntLeft. ---*/
            gradSolIntLeft[0](i,l) = 0.5*(gradSolIntLeft[0](i,l) + gradSolIntRight[0](i,l));
            gradSolIntLeft[1](i,l) = 0.5*(gradSolIntLeft[1](i,l) + gradSolIntRight[1](i,l));
          }
        }

        break;
      }

      case 3: {

        /*--- Three dimensional simulation. Easier storage of the metric terms. ---*/
        ColMajorMatrix<su2double> &dParDx0 =  matchingInternalFaces[l].metricCoorDerivFace0[0];
        ColMajorMatrix<su2double> &dParDy0 =  matchingInternalFaces[l].metricCoorDerivFace0[1];
        ColMajorMatrix<su2double> &dParDz0 =  matchingInternalFaces[l].metricCoorDerivFace0[2];

        ColMajorMatrix<su2double> &dParDx1 =  matchingInternalFaces[l].metricCoorDerivFace1[0];
        ColMajorMatrix<su2double> &dParDy1 =  matchingInternalFaces[l].metricCoorDerivFace1[1];
        ColMajorMatrix<su2double> &dParDz1 =  matchingInternalFaces[l].metricCoorDerivFace1[2];

        /*--- Loop over the number of variables and padded number of integration points. ---*/
        for(unsigned short l=0; l<nVar; ++l) {
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nIntPad; ++i) {

            /*--- Left state. Store the gradients w.r.t. the parametric coordinates. ---*/
            su2double dvdr = gradSolIntLeft[0](i,l);
            su2double dvds = gradSolIntLeft[1](i,l);
            su2double dvdt = gradSolIntLeft[2](i,l);

            /*--- Compute the Cartesian gradients of the left state. ---*/
            gradSolIntLeft[0](i,l) = dvdr*dParDx0(i,0) + dvds*dParDx0(i,1) + dvdt*dParDx0(i,2);
            gradSolIntLeft[1](i,l) = dvdr*dParDy0(i,0) + dvds*dParDy0(i,1) + dvdt*dParDy0(i,2);
            gradSolIntLeft[2](i,l) = dvdr*dParDz0(i,0) + dvds*dParDz0(i,1) + dvdt*dParDz0(i,2);

            /*--- Right state. Store the gradients w.r.t. the parametric coordinates. ---*/
            dvdr = gradSolIntRight[0](i,l);
            dvds = gradSolIntRight[1](i,l);
            dvdt = gradSolIntRight[2](i,l);

            /*--- Compute the Cartesian gradients of the right state. ---*/
            gradSolIntRight[0](i,l) = dvdr*dParDx1(i,0) + dvds*dParDx1(i,1) + dvdt*dParDx1(i,2);
            gradSolIntRight[1](i,l) = dvdr*dParDy1(i,0) + dvds*dParDy1(i,1) + dvdt*dParDy1(i,2);
            gradSolIntRight[2](i,l) = dvdr*dParDz1(i,0) + dvds*dParDz1(i,1) + dvdt*dParDz1(i,2);

            /*--- Store the averaged gradients in gradSolIntLeft. ---*/
            gradSolIntLeft[0](i,l) = 0.5*(gradSolIntLeft[0](i,l) + gradSolIntRight[0](i,l));
            gradSolIntLeft[1](i,l) = 0.5*(gradSolIntLeft[1](i,l) + gradSolIntRight[1](i,l));
            gradSolIntLeft[2](i,l) = 0.5*(gradSolIntLeft[2](i,l) + gradSolIntRight[2](i,l));
          }
        }

        break;
      }
    }
    
    /*--- Compute the invisid fluxes in the integration points of the face. ---*/
    const unsigned int indFlux = omp_get_num_threads() + omp_get_thread_num();
    ColMajorMatrix<su2double> &fluxes = matchingInternalFaces[l].standardElemFlow->elem0->workSolInt[indFlux];
    ComputeInviscidFluxesFace(config, solIntLeft, solIntRight, matchingInternalFaces[l].JacobiansFace,
                              matchingInternalFaces[l].metricNormalsFace,
                              matchingInternalFaces[l].gridVelocities, numerics, fluxes);

    /*--- Compute the viscous fluxes and the symmetrizing fluxes in the
          integration points of the face. ---*/
  

    /*----------------------------------------------------------------------*/
    /*--- Step 2: Create the final form of the symmetrizing fluxes and   ---*/
    /*---         multiply all fluxes by the integration weights.        ---*/
    /*----------------------------------------------------------------------*/

    /*--- The symmetrizing fluxes are multiplied by the gradients of the basis
          functions for the residual. These must be the Cartesian gradients,
          but for the standard element the gradients w.r.t. the parametric
          coordinates are stored. To account for this, the Cartesian symmetrizing
          fluxes are converted parametric fluxes. Note that this conversion is
          usually different for the left and the right element. Therefore the final
          symmetrizing fluxes for the left element are stored in gradSolIntLeft
          and for the right element in gradSolIntRight. ---*/
  }
  END_SU2_OMP_FOR

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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::BC_Euler_Wall(CConfig             *config,
                                     const unsigned long surfElemBeg,
                                     const unsigned long surfElemEnd,
                                     CSurfaceElementFEM  *surfElem,
                                     CNumerics           *conv_numerics){
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::BC_Far_Field(CConfig             *config,
                                    const unsigned long surfElemBeg,
                                    const unsigned long surfElemEnd,
                                    CSurfaceElementFEM  *surfElem,
                                    CNumerics           *conv_numerics){
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::BC_Sym_Plane(CConfig             *config,
                                    const unsigned long surfElemBeg,
                                    const unsigned long surfElemEnd,
                                    CSurfaceElementFEM  *surfElem,
                                    CNumerics           *conv_numerics){
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::BC_Supersonic_Outlet(CConfig             *config,
                                            const unsigned long surfElemBeg,
                                            const unsigned long surfElemEnd,
                                            CSurfaceElementFEM  *surfElem,
                                            CNumerics           *conv_numerics){
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::BC_Inlet(CConfig             *config,
                                const unsigned long surfElemBeg,
                                const unsigned long surfElemEnd,
                                CSurfaceElementFEM  *surfElem,
                                CNumerics           *conv_numerics,
                                unsigned short      val_marker) {
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::BC_Outlet(CConfig             *config,
                                 const unsigned long surfElemBeg,
                                 const unsigned long surfElemEnd,
                                 CSurfaceElementFEM  *surfElem,
                                 CNumerics           *conv_numerics,
                                 unsigned short      val_marker) {
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::BC_HeatFlux_Wall(CConfig             *config,
                                        const unsigned long surfElemBeg,
                                        const unsigned long surfElemEnd,
                                        CSurfaceElementFEM  *surfElem,
                                        CNumerics           *conv_numerics,
                                        unsigned short      val_marker) {
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::BC_Isothermal_Wall(CConfig             *config,
                                          const unsigned long surfElemBeg,
                                          const unsigned long surfElemEnd,
                                          CSurfaceElementFEM  *surfElem,
                                          CNumerics           *conv_numerics,
                                          unsigned short      val_marker) {
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::BC_Riemann(CConfig             *config,
                                  const unsigned long surfElemBeg,
                                  const unsigned long surfElemEnd,
                                  CSurfaceElementFEM  *surfElem,
                                  CNumerics           *conv_numerics,
                                  unsigned short      val_marker) {
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_NSSolver::BC_Custom(CConfig             *config,
                                 const unsigned long surfElemBeg,
                                 const unsigned long surfElemEnd,
                                 CSurfaceElementFEM  *surfElem,
                                 CNumerics           *conv_numerics) {
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
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}