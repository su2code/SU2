/*!
 * \file COutputLegacy.hpp
 * \brief Headers of the main subroutines for generating the file outputs.
 *        The subroutines and functions are in the <i>output_structure_legacy.cpp</i> file.
 * \author F. Palacios, T. Economon, M. Colonno
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/parallelization/mpi_structure.hpp"

#ifdef HAVE_CGNS
  #include "cgnslib.h"
#endif
#ifdef HAVE_TECIO
  #include "TECIO.h"
#endif
#include <fstream>
#include <cmath>
#include <vector>

#include "../../../Common/include/option_structure.hpp"
class CGeometry;
class CConfig;
class CSolver;
class CIntegration;

using namespace std;

/*!
 * \class COutputLegacy
 * \brief Class for writing the flow, adjoint and linearized solver
 *        solution (including the history solution, and parallel stuff).
 * \author F. Palacios, T. Economon, M. Colonno.
 */
class COutputLegacy {

  su2double Sum_Total_RadialDistortion, Sum_Total_CircumferentialDistortion; // Add all the distortion to compute a run average.
  bool turbo;
  unsigned short nSpanWiseSections, nMarkerTurboPerf;

  su2double **TotalStaticEfficiency,
        **TotalTotalEfficiency,
        **KineticEnergyLoss,
        **TRadius,
        **TotalPressureLoss,
        **MassFlowIn,
        **MassFlowOut,
        **FlowAngleIn,
        **FlowAngleIn_BC,
        **FlowAngleOut,
        **EulerianWork,
        **TotalEnthalpyIn,
        **TotalEnthalpyIn_BC,
        **EntropyIn,
        **EntropyOut,
        **EntropyIn_BC,
        **PressureRatio,
        **TotalTemperatureIn,
        **EnthalpyOut,
        ***MachIn,
        ***MachOut,
        **VelocityOutIs,
        **DensityIn,
        **PressureIn,
        ***TurboVelocityIn,
        **DensityOut,
        **PressureOut,
        ***TurboVelocityOut,
        **EnthalpyOutIs,
        **EntropyGen,
        **AbsFlowAngleIn,
        **TotalEnthalpyOut,
        **RothalpyIn,
        **RothalpyOut,
        **TotalEnthalpyOutIs,
        **AbsFlowAngleOut,
        **PressureOut_BC,
        **TemperatureIn,
        **TemperatureOut,
        **TotalPressureIn,
        **TotalPressureOut,
        **TotalTemperatureOut,
        **EnthalpyIn,
        **TurbIntensityIn,
        **Turb2LamViscRatioIn,
        **TurbIntensityOut,
        **Turb2LamViscRatioOut,
        **NuFactorIn,
        **NuFactorOut;

protected:

  int rank,   /*!< \brief MPI Rank. */
  size;         /*!< \brief MPI Size. */

public:

  /*!
   * \brief Constructor of the class.
   */
  COutputLegacy(CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~COutputLegacy(void);

  /*!
   * \brief Writes inverse design.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iExtIter - Current external (time) iteration.
   */
  void SetCp_InverseDesign(CSolver *solver_container, CGeometry *geometry, CConfig *config,
                         unsigned long iExtIter);

  /*!
   * \brief Writes inverse design.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iExtIter - Current external (time) iteration.
   */
  void SetHeatFlux_InverseDesign(CSolver *solver_container, CGeometry *geometry, CConfig *config,
                        unsigned long iExtIter);

  /*!
   * \brief Writes forces at different sections.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] output - Create output files.
   */
  void SpecialOutput_SpanLoad(CSolver *solver, CGeometry *geometry, CConfig *config, bool output) const;

  /*!
   * \brief Writes one dimensional output.
   * \author H. Kline
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] output - Create output files.
   */
  void SpecialOutput_AnalyzeSurface(CSolver *solver, CGeometry *geometry, CConfig *config, bool output) const;

  /*!
   * \brief Create and write the file with the flow coefficient on the surface.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] FlowSolution - Flow solution.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Current zone number in the grid file.
   * \param[in] output - Create output files.
   */
  void SpecialOutput_Distortion(CSolver *solver, CGeometry *geometry, CConfig *config, bool output) const;

  /*!
   * \brief Write the header of the history file.
   * \param[in] ConvHist_file - Pointer to the convergence history file (which is defined in the main subroutine).
   * \param[in] config - Definition of the particular problem.
   */
  void SetConvHistory_Header(ofstream *ConvHist_file, CConfig *config, unsigned short val_iZone, unsigned short val_iInst);

  /*!
   * \brief Write the history file and the convergence on the screen for serial computations.
   * \param[in] ConvHist_file - Pointer to the convergence history file (which is defined in the main subroutine).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] integration - Generic subroutines for space integration, time integration, and monitoring.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] timeused - Current number of clock tick in the computation (related with total time).
   * \param[in] val_nZone - iZone index.
   */
  void SetConvHistory_Body(ofstream *ConvHist_file, CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
                              CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst);

  /*!
   * \brief Write the history file and the convergence on the screen for serial computations.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] output - Create output files.
   */
  void SpecialOutput_ForcesBreakdown(CSolver *****solver, CGeometry ****geometry, CConfig **config, unsigned short val_iZone, bool output) const;

  /*!
   * \brief Compute .
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iExtIter - Current external (time) iteration.
   */
  void ComputeTurboPerformance(CSolver *solver_container, CGeometry *geometry, CConfig *config);

  /*!
   * \brief Compute .
   * \param[in] config - Definition of the particular problem.
   */
  void WriteTurboPerfConvHistory(CConfig *config);

  /*!
   * \brief Write the output file for spanwise turboperformance.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - iZone index.
   * \param[in] output - Create output files.
   */
  void SpecialOutput_Turbo(CSolver *****solver_container, CGeometry ****geometry, CConfig **config, unsigned short val_iZone, bool output);

  /*!
   * \brief Give the Entropy Generation performance parameters for turbomachinery.
   * \param[in] iMarkerTP - Marker turbo-performance.
   * \param[in] iSpan - span section.
   */
  su2double GetEntropyGen(unsigned short iMarkerTP, unsigned short iSpan);

  /*!
   * \brief Give the Entropy Generation performance parameters for turbomachinery.
   * \param[in] iMarkerTP - Marker turbo-performance.
   * \param[in] iSpan - span section.
   */
  su2double GetFlowAngleOut(unsigned short iMarkerTP, unsigned short iSpan);

  /*!
   * \brief Give the Entropy Generation performance parameters for turbomachinery.
   * \param[in] iMarkerTP - Marker turbo-performance.
   * \param[in] iSpan - span section.
   */
  su2double GetMassFlowIn(unsigned short iMarkerTP, unsigned short iSpan);

  /*!
   * \brief Write the output file for harmonic balance for each time-instance.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Number of Zones.
   * \param[in] val_iZone - Zone index.
   * \param[in] output - Create output files.
   */
  void SpecialOutput_HarmonicBalance(CSolver *****solver, CGeometry ****geometry, CConfig **config, unsigned short iZone, unsigned short val_nZone, bool output) const;

  /*!
   * \brief Writes the special output files.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Total number of domains in the grid file.
   * \param[in] val_nZone - Total number of domains in the grid file.
   */
  void SetSpecial_Output(CSolver *****solver_container, CGeometry ****geometry, CConfig **config,
                         unsigned long iExtIter, unsigned short val_nZone);

};

inline su2double COutputLegacy::GetEntropyGen(unsigned short iMarkerTP, unsigned short iSpan) { return EntropyGen[iMarkerTP][iSpan]; }

inline su2double COutputLegacy::GetFlowAngleOut(unsigned short iMarkerTP, unsigned short iSpan) { return FlowAngleOut[iMarkerTP][iSpan]*180.0/PI_NUMBER; }

inline su2double COutputLegacy::GetMassFlowIn(unsigned short iMarkerTP, unsigned short iSpan) { return MassFlowIn[iMarkerTP][iSpan]; }
