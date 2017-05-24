/*!
 * \file postprocessing_structure.hpp
 * \brief Headers of the post processing structure.
 * \author T. Albring, Beckett Y. Zhou
 * \version 4.3.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#pragma once


#include "../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

//#include "fluid_model.hpp"
//#include "numerics_structure.hpp"
//#include "variable_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/matrix_structure.hpp"
#include "../../Common/include/vector_structure.hpp"
#include "../../Common/include/linear_solvers_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../SU2_CFD/include/solver_structure.hpp"
//#include "numerics_machine_learning.hpp"

using namespace std;


class FWHSolver {
  public:
  su2double  CFD_PressureFluctuation;
  su2double  CAA_PressureFluctuation;
  su2double U1, U2, U3, M, a_inf, AOA, beta_sq, FreeStreamDensity, FreeStreamPressure;
  complex <su2double> ***G, ***dGdy1, ***dGdy2, ***dGdy3,**fpp;
  su2double ***dJdU;
  su2double **surface_geo;
  complex <su2double>  **fF1, **fF2, **fQ  ;
  su2double   **F1, **F2, **Q, *F1_mean, *F2_mean, *Q_mean, *Hanning_W, Hanning_Scale;
  su2double ***F, **F_mean;
  complex <su2double>  ***fF;
  complex <su2double>  **pp;
  su2double  **fpp_r, **fpp_i, **fpp_r_root,**fpp_i_root, **pp_CFD, *pp_CFD_mean ;
  unsigned long nObserver, nPanel, nSample, idx_window_l, idx_window_r, nDim, SamplingFreq;
  unsigned long totFWH;
  unsigned long *PointID;
  su2double **Observer_Locations, **closet_coord_AllObs;
  su2double SPL;
  su2double ***Fr, **Fr_mean, ***pp_ret, ***t_Obs, **t_interp, ***pp_interp, **pp_TimeDomain , **pp_TimeDomain_root, **r_minmax;
  bool TimeDomain3D, UseAnalytic;
  su2double T_a, Freq_a, Amp_a;


	/*!
	 * \brief Constructor of the  class.
	 */
	FWHSolver(CConfig *config, CGeometry *geometry);

	/*!
	 * \brief Destructor of the class.
	 */
	~FWHSolver(void);



        void SetAeroacoustic_Analysis(CSolver *solver, CConfig *config, CGeometry *geometry, ofstream &CFD_pressure_file);
        void SetCFD_PressureFluctuation(CSolver *solver, CConfig *config, CGeometry *geometry, unsigned long iObserver);
        void Extract_NoiseSources(CSolver *solver, CConfig* config, CGeometry *geometry);
        void Compute_FarfieldNoise(CSolver *solver, CConfig* config, CGeometry *geometry);
        void Compute_TimeDomainPanelSignal(CConfig* config);
        void Compute_ObserverTime(CConfig* config);
        void Interpolate_PressureSignal(CGeometry *geometry);
        void Compute_GreensFunction2D (CConfig* config);
        void Compute_GreensFunction3D (CConfig* config);
        void Window_SourceTerms ();
        void FFT_SourceTermsR2 ();
        void Integrate_Sources (CConfig* config);
        void iFFT_SignalR2 ();
        void Compute_SPL ();
        void Write_Sensitivities(CSolver *solver, CConfig *config, CGeometry *geometry);
        su2double bessj0 (su2double x);
        su2double bessj1 (su2double x);
        su2double bessy0 (su2double x);
        su2double bessy1 (su2double x);
        su2double GetCFD_PressureFluctuation();

//        void iFFT_Signal (CSolver *solver, CConfig* config, CGeometry *geometry);
//        void FFT_SourceTerms (CSolver *solver, CConfig* config, CGeometry *geometry);
//        void SetCAA_PressureFluctuation(CSolver *solver, CConfig* config, CGeometry *geometry, ofstream &SRC_p_file, ofstream &SRC_ux_file, ofstream &SRC_uy_file, ofstream &SRC_rho_file, ofstream &time_file);
//        su2double GetCAA_PressureFluctuation();

};



#include "postprocessing_structure.inl"
