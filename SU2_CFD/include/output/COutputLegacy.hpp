/*!
 * \file output.hpp
 * \brief Headers of the main subroutines for generating the file outputs.
 *        The subroutines and functions are in the <i>output_structure.cpp</i> file.
 * \author F. Palacios, T. Economon, M. Colonno
 * \version 6.1.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#pragma once

#include "../../../Common/include/mpi_structure.hpp"

#ifdef HAVE_CGNS
  #include "cgnslib.h"
#endif
#ifdef HAVE_TECIO
  #include "TECIO.h"
#endif
#include <fstream>
#include <cmath>
#include <time.h>
#include <fstream>
#include <vector>

//#include "../solver_structure.hpp"
//#include "../integration_structure.hpp"
//#include "../../../Common/include/geometry_structure.hpp"
//#include "../../../Common/include/fem_geometry_structure.hpp"
//#include "../../../Common/include/fem_standard_element.hpp"

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

  unsigned long nGlobal_Poin;   // Global number of nodes with halos
  unsigned long nSurf_Poin;   // Global number of nodes of the surface
  unsigned long nGlobal_Doma;   // Global number of nodes without halos
  unsigned long nGlobal_Elem;  // Global number of elems without halos
  unsigned long nSurf_Elem,  // Global number of surface elems without halos
  nGlobal_Line,
  nGlobal_BoundTria,
  nGlobal_BoundQuad,
  nGlobal_Tria,
  nGlobal_Quad,
  nGlobal_Tetr,
  nGlobal_Hexa,
  nGlobal_Pris,
  nGlobal_Pyra;
  su2double **Coords;              // node i (x, y, z) = (Coords[0][i], Coords[1][i], Coords[2][i])
  int *Conn_Line;
  int *Conn_BoundTria;
  int *Conn_BoundQuad;
  int *Conn_Tria;  // triangle 1 = Conn_Tria[0], Conn_Tria[1], Conn_Tria[3]
  int *Conn_Quad;
  int *Conn_Tetr;
  int *Conn_Hexa;
  int *Conn_Pris;
  int *Conn_Pyra;
  
  
  unsigned long nGlobal_Poin_Par;   // Global number of nodes with halos
  unsigned long nGlobal_Elem_Par;  // Global number of elems without halos
  unsigned long nGlobal_Surf_Poin;
  unsigned long nSurf_Elem_Par;
  unsigned long nSurf_Poin_Par;
  unsigned long nParallel_Poin;
  unsigned long nParallel_Line,
  nParallel_BoundTria,
  nParallel_BoundQuad,
  nParallel_Tria,
  nParallel_Quad,
  nParallel_Tetr,
  nParallel_Hexa,
  nParallel_Pris,
  nParallel_Pyra;
  int *Conn_BoundLine_Par;
  int *Conn_BoundTria_Par;
  int *Conn_BoundQuad_Par;
  int *Conn_Tria_Par;  // triangle 1 = Conn_Tria[0], Conn_Tria[1], Conn_Tria[3]
  int *Conn_Quad_Par;
  int *Conn_Tetr_Par;
  int *Conn_Hexa_Par;
  int *Conn_Pris_Par;
  int *Conn_Pyra_Par;

  unsigned long nGlobalPoint_Sort;
  unsigned long nLocalPoint_Sort;
  unsigned long nPoint_Restart;
  int *Local_Halo_Sort;

  unsigned long *beg_node;
  unsigned long *end_node;

  unsigned long *nPointLinear;
  unsigned long *nPointCumulative;
  
  unsigned short nVar_Par;
  su2double **Local_Data;
  su2double **Local_Data_Copy;      // Local data copy for cte. lift mode
  su2double **Parallel_Data;        // node i (x, y, z) = (Coords[0][i], Coords[1][i], Coords[2][i])
  su2double **Parallel_Surf_Data;   // node i (x, y, z) = (Coords[0][i], Coords[1][i], Coords[2][i])
  vector<string> Variable_Names;

  su2double **Data;
  unsigned short nVar_Consv, nVar_Total, nVar_Extra, nZones;
  bool wrote_surf_file, wrote_CGNS_base, wrote_Tecplot_base, wrote_Paraview_base;
  unsigned short wrote_base_file;
  su2double RhoRes_New, *RhoRes_Old;
  int cgns_base, cgns_zone, cgns_base_results, cgns_zone_results;
  
  su2double Sum_Total_RadialDistortion, Sum_Total_CircumferentialDistortion; // Add all the distortion to compute a run average.
  bool turbo;
  unsigned short   nSpanWiseSections,
		   nMarkerTurboPerf;

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

  unsigned long nMarker_InletFile;       /*!< \brief Counter for total number of inlet boundaries written to inlet profile file. */
  vector<string> Marker_Tags_InletFile;   /*!< \brief Marker tags for the strings of the markers in the inlet profile file. */
  unsigned long *nRow_InletFile;         /*!< \brief Counters for the number of points per marker in the inlet profile file. */
  unsigned long *nRowCum_InletFile;      /*!< \brief Counters for the number of points per marker in cumulative storage format in the inlet profile file. */
  su2double **InletCoords;  /*!< \brief Data structure for holding the merged inlet boundary coordinates from all ranks. */

protected:

  int rank, 	/*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */

public:

  /*! 
   * \brief Constructor of the class. 
   */
  COutputLegacy(CConfig *congig);

  /*! 
   * \brief Destructor of the class. 
   */
  ~COutputLegacy(void);

  /*! 
   * \brief Writes and organizes the all the output files, except the history one, for serial computations.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Total number of domains in the grid file.
   * \param[in] val_nZone - Total number of domains in the grid file.
   */
  void SetResult_Files(CSolver *****solver_container, CGeometry ****geometry, CConfig **config,
                       unsigned long iExtIter, unsigned short val_nZone);
  
  /*!
   * \brief Writes and organizes the all the output files, except the history one, for serial computations with the FEM solver.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Total number of domains in the grid file.
   * \param[in] val_nZone - Total number of domains in the grid file.
   */
  void SetResult_Files_FEM(CSolver ****solver_container, CGeometry ***geometry, CConfig **config,
                       unsigned long iExtIter, unsigned short val_nZone);
  
  /*!
   * \brief Writes and organizes the all the output files, except the history one, for serial computations.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Total number of domains in the grid file.
   * \param[in] val_nZone - Total number of domains in the grid file.
   */
  void SetBaselineResult_Files(CSolver ***solver, CGeometry ***geometry, CConfig **config,
                               unsigned long iExtIter, unsigned short val_nZone);

  /*!
   * \brief Writes and organizes the all the output files, except the history one, for DG-FEM simulations (SU2_SOL).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Total number of domains in the grid file.
   * \param[in] val_nZone - Total number of domains in the grid file.
   */
  void SetBaselineResult_Files_FEM(CSolver ***solver, CGeometry ***geometry, CConfig **config,
                                   unsigned long iExtIter, unsigned short val_nZone);

  /*!
   * \brief Writes and organizes the all the output files, except the history one, for serial computations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Total number of domains in the grid file.
   */
  void SetMesh_Files(CGeometry **geometry, CConfig **config, unsigned short val_nZone, bool new_file, bool su2_file);

  /*!
   * \brief Writes equivalent area.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] output - Create output files.
   */
  void SpecialOutput_SonicBoom(CSolver *solver, CGeometry *geometry, CConfig *config, bool output);
  
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
  void SpecialOutput_SpanLoad(CSolver *solver, CGeometry *geometry, CConfig *config, bool output);
  
  /*!
   * \brief Writes one dimensional output.
   * \author H. Kline
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] output - Create output files.
   */
  void SpecialOutput_AnalyzeSurface(CSolver *solver, CGeometry *geometry, CConfig *config, bool output);
  
  /*!
   * \brief Create and write the file with the flow coefficient on the surface.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] FlowSolution - Flow solution.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Current zone number in the grid file.
   * \param[in] output - Create output files.
   */
  void SpecialOutput_Distortion(CSolver *solver, CGeometry *geometry, CConfig *config, bool output);

  /*!
   * \brief Create and write the file with the FSI convergence history.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Solver for all physical problems.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Current zone number in the grid file.
   */
  void SpecialOutput_FSI(ofstream *FSIHist_file, CGeometry ****geometry, CSolver *****solver_container, CConfig **config, CIntegration ****integration,
                         unsigned long iExtIter, unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT, bool header);

  /*!
   * \brief Create and write the file with the FSI convergence history.
   * \param[in] iIter - Current iteration.
   * \param[in] iFreq - Frequency of output printing.
   */
  bool PrintOutput(unsigned long iIter, unsigned long iFreq);

  /*! 
   * \brief Create and write the file with the flow coefficient on the surface.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] FlowSolution - Flow solution.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Current zone number in the grid file.
   */
  void SetSurfaceCSV_Flow(CConfig *config, CGeometry *geometry, CSolver *FlowSolver, unsigned long iExtIter, unsigned short val_iZone, unsigned short val_iInst);

  /*! 
   * \brief Create and write the file with the adjoint coefficients on the surface for serial computations.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] AdjSolution - Adjoint solution.
   * \param[in] FlowSolution - Flow solution.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Current zone number in the grid file.
   */
  void SetSurfaceCSV_Adjoint(CConfig *config, CGeometry *geometry, CSolver *AdjSolver, CSolver *FlowSolution, unsigned long iExtIter, unsigned short val_iZone, unsigned short val_iInst);

  /*!
   * \brief Merge the geometry into a data structure used for output file writing.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_nZone - iZone index.
   */
  void MergeConnectivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone);
  
  /*!
   * \brief Merge the FEM geometry into a data structure used for output file writing.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_nZone - iZone index.
   */
  void MergeConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short val_iZone);
  
  /*!
   * \brief Merge the node coordinates from all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void MergeCoordinates(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Merge the node coordinates from all processors for the FEM solver.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void MergeCoordinates_FEM(CConfig *config, CGeometry *geometry);
  
  /*!
   * \brief Merge the connectivity for a single element type from all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void MergeVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);
  
  /*!
   * \brief Merge the connectivity for a single element type from all processors for the FEM solver.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void MergeVolumetricConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);
  
  /*!
   * \brief Merge the connectivity for a single element type from all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void MergeSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);
  
  /*!
   * \brief Merge the connectivity for a single element type from all processors for the FEM solver.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void MergeSurfaceConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);
  
  /*!
   * \brief Merge the solution into a data structure used for output file writing.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Flow, adjoint or linearized solution.
   * \param[in] val_nZone - iZone index.
   */
  void MergeSolution(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone);

  /*!
   * \brief Merge the FEM solution into a data structure used for output file writing.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Flow, adjoint or linearized solution.
   * \param[in] val_nZone - iZone index.
   */
  void MergeSolution_FEM(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone);
  
  /*!
   * \brief Merge the solution into a data structure used for output file writing.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Flow, adjoint or linearized solution.
   * \param[in] val_nZone - iZone index.
   */
  void MergeBaselineSolution(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone);

  /*!
   * \brief Merge the solution into a data structure used for output file writing (DG-FEM).
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Flow, adjoint or linearized solution.
   * \param[in] val_nZone - iZone index.
   */
  void MergeBaselineSolution_FEM(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone);

  /*!
   * \brief Write a native SU2 restart file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Flow, adjoint or linearized solution.
   * \param[in] val_iZone - iZone index.
   */
  void SetRestart(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone);

  /*!
   * \brief Write a native SU2 restart file (ASCII) in parallel.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Flow, adjoint or linearized solution.
   * \param[in] val_iZone - iZone index.
   */
  void WriteRestart_Parallel_ASCII(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_iInst);

  /*!
   * \brief Write a native SU2 restart file (binary) in parallel.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Flow, adjoint or linearized solution.
   * \param[in] val_iZone - iZone index.
   */
  void WriteRestart_Parallel_Binary(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_iInst);

  /*!
   * \brief Write the x, y, & z coordinates to a CGNS output file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetCGNS_Coordinates(CConfig *config, CGeometry *geometry, unsigned short val_iZone);
  
  /*!
   * \brief Write the element connectivity to a CGNS output file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetCGNS_Connectivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone);
  
  /*!
   * \brief Write solution data to a CGNS output file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetCGNS_Solution(CConfig *config, CGeometry *geometry, unsigned short val_iZone);
  
//  /*!
//   * \brief Write a Paraview ASCII solution file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - Current zone.
//   * \param[in] val_nZone - Total number of zones.
//   */
//  void SetParaview_ASCII(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol);

//  /*!
//   * \brief Write a Paraview ASCII solution file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - Current zone.
//   * \param[in] val_nZone - Total number of zones.
//   */
//  void SetParaview_MeshASCII(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol, bool new_file);

//  /*!
//   * \brief Write a Paraview ASCII solution file with parallel output.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - Current zone.
//   * \param[in] val_nZone - Total number of zones.
//   * \param[in] surf_sol - Flag controlling whether this is a volume or surface file.
//   */
//  void WriteParaViewASCII_Parallel(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_nZone, unsigned short val_iInst, unsigned short val_nInst, bool surf_sol);

//  /*!
//   * \brief Write a Paraview binary solution file with parallel output.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - Current zone.
//   * \param[in] val_nZone - Total number of zones.
//   * \param[in] surf_sol - Flag controlling whether this is a volume or surface file.
//   */
//  void WriteParaViewBinary_Parallel(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol);
  
//  /*!
//   * \brief Write a Tecplot ASCII solution file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - Current zone.
//   * \param[in] val_nZone - Total number of zones.
//   */
//  void SetTecplotASCII(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol);
  
//  /*!
//   * \brief Write the nodal coordinates and connectivity to a Tecplot ASCII mesh file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - iZone index.
//   */
//  void SetTecplot_MeshASCII(CConfig *config, CGeometry *geometry, bool surf_sol, bool new_file);

//  /*!
//   * \brief Write the nodal coordinates and connectivity to a stl ASCII mesh file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - iZone index.
//   */
//  void SetSTL_MeshASCII(CConfig *config, CGeometry *geometry);

//  /*!
//   * \brief Write the nodal coordinates and connectivity to a n3d ASCII mesh file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   */
//  void SetCSV_MeshASCII(CConfig *config, CGeometry *geometry);
  
//  /*!
//   * \brief Write the nodal coordinates and connectivity to a n3d ASCII mesh file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - iZone index.
//   */
//  void SetTecplotASCII_Mesh(CConfig *config, CGeometry *geometry, unsigned short val_iZone, bool surf_sol, bool new_file);

//  /*!
//   * \brief Write the solution data and connectivity to a Tecplot ASCII mesh file in parallel.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - Current zone.
//   * \param[in] val_nZone - Total number of zones.
//   * \param[in] surf_sol - Flag controlling whether this is a volume or surface file.
//   */
//  void WriteTecplotASCII_Parallel(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_nZone, unsigned short val_iInst, unsigned short val_nInst, bool surf_sol);
  
  /*!
   * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  string AssembleVariableNames(CGeometry *geometry, CConfig *config, unsigned short nVar_Consv, unsigned short *NVar);

  /*!
   * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetSU2_MeshASCII(CConfig *config, CGeometry *geometry, unsigned short val_iZone, ofstream &output_file);
  
  /*!
   * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetSU2_MeshBinary(CConfig *config, CGeometry *geometry);
  
  /*!
   * \brief Write the nodal coordinates to a binary file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void WriteCoordinates_Binary(CConfig *config, CGeometry *geometry, unsigned short val_iZone);

//  /*!
//   * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - iZone index.
//   */
//  void SetTecplotBinary_DomainMesh(CConfig *config, CGeometry *geometry, unsigned short val_iZone);
  
//  /*!
//   * \brief Write the coordinates and connectivity to a Tecplot binary surface mesh file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - iZone index.
//   */
//  void SetTecplotBinary_SurfaceMesh(CConfig *config, CGeometry *geometry, unsigned short val_iZone);
  
//  /*!
//   * \brief Write solution data to a Tecplot binary volume solution file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - iZone index.
//   */
//  void SetTecplotBinary_DomainSolution(CConfig *config, CGeometry *geometry, unsigned short val_iZone);

//  /*!
//   * \brief Write solution data to a Tecplot binary surface solution file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - iZone index.
//   */
//  void SetTecplotBinary_SurfaceSolution(CConfig *config, CGeometry *geometry, unsigned short val_iZone);
  
//  /*!
//   * \brief Write a Tecplot ASCII solution file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - Current zone.
//   * \param[in] val_nZone - Total number of zones.
//   */
//  void SetFieldViewASCII(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone);
  
//  /*!
//   * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - iZone index.
//   */
//  void SetFieldViewASCII_Mesh(CConfig *config, CGeometry *geometry);
  
//  /*!
//   * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - iZone index.
//   */
//  void SetFieldViewBinary_Mesh(CConfig *config, CGeometry *geometry);
  
//  /*!
//   * \brief Write solution data to a Tecplot binary volume solution file.
//   * \param[in] config - Definition of the particular problem.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] val_iZone - iZone index.
//   */
//  void SetFieldViewBinary(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone);
  
  /*!
   * \brief Deallocate temporary memory needed for merging and writing coordinates.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void DeallocateCoordinates(CConfig *config, CGeometry *geometry);
  
  /*!
   * \brief Deallocate temporary memory needed for merging and writing connectivity.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void DeallocateConnectivity(CConfig *config, CGeometry *geometry, bool surf_sol);
  
  /*!
   * \brief Deallocate temporary memory needed for merging and writing solution variables.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void DeallocateSolution(CConfig *config, CGeometry *geometry);
  
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
  void SpecialOutput_ForcesBreakdown(CSolver *****solver, CGeometry ****geometry, CConfig **config, unsigned short val_iZone, bool output);
  
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
  void SetCFL_Number(CSolver *****solver_container, CConfig **config, unsigned short val_iZone);
  
  /*!
   * \brief Write the sensitivity (including mesh sensitivity) computed with the discrete adjoint method
   *  on the surface and in the volume to a file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Number of Zones.
   */
  void SetSensitivity_Files(CGeometry ***geometry, CConfig **config, unsigned short val_nZone);

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
  void SpecialOutput_HarmonicBalance(CSolver *****solver, CGeometry ****geometry, CConfig **config, unsigned short iZone, unsigned short val_nZone, bool output);

  /*!
   * \brief Writes and organizes the all the output files, except the history one, for parallel computations.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Total number of domains in the grid file.
   * \param[in] val_nZone - Total number of domains in the grid file.
   */
  void SetResult_Files_Parallel(CSolver *****solver_container, CGeometry ****geometry, CConfig **config,
                                unsigned long iExtIter, unsigned short val_nZone);
  
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

  /*!
   * \brief Load the desired solution data into a structure used for parallel reordering and output file writing for flow problems.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Flow, adjoint or linearized solution.
   * \param[in] val_nZone - iZone index.
   */
  void LoadLocalData_Flow(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone);
  
  /*!
   * \brief Load the desired solution data into a structure used for parallel reordering and output file writing for incmopressible flow problems.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Flow, adjoint or linearized solution.
   * \param[in] val_iZone - iZone index.
   */
  void LoadLocalData_IncFlow(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone);

  /*!
   * \brief Load the desired solution data into a structure used for parallel reordering and output file writing for adjoint flow problems.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Flow, adjoint or linearized solution.
   * \param[in] val_nZone - iZone index.
   */
  void LoadLocalData_AdjFlow(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone);
  
  /*!
   * \brief Load the desired solution data into a structure used for parallel reordering and output file writing for elasticity problems.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Flow, adjoint or linearized solution.
   * \param[in] val_nZone - iZone index.
   */
  void LoadLocalData_Elasticity(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone);
  
  /*!
   * \brief Load the desired solution data into a structure used for parallel reordering and output file writing for generic problems.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Flow, adjoint or linearized solution.
   * \param[in] val_nZone - iZone index.
   */
  void LoadLocalData_Base(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone);

  /*!
   * \brief Load the desired solution data into a structure used for parallel reordering and output file writing for DG-FEM flow problems.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Flow, adjoint or linearized solution.
   * \param[in] val_nZone - iZone index.
   */
  void LoadLocalData_FEM(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone);

  /*!
   * \brief Prepare the number of points and offsets for linear partitioning that are needed for output.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void PrepareOffsets(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Sort the connectivities (volume and surface) into data structures used for output file writing.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_nZone - iZone index.
   */
  void SortConnectivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone);

  /*!
   * \brief Sort the connectivities (volume and surface) into data structures used for output file writing (DG-FEM solver).
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_nZone - iZone index.
   */
  void SortConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short val_iZone);

  /*!
   * \brief Sort the connectivity for a single volume element type into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void SortVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);

  /*!
   * \brief Sort the connectivity for a single volume element type into a linear partitioning across all processors (DG-FEM solver).
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void SortVolumetricConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);

  /*!
   * \brief Sort the connectivity for a single surface element type into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void SortSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);

  /*!
   * \brief Sort the connectivity for a single surface element type into a linear partitioning across all processors (DG-FEM solver).
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   */
  void SortSurfaceConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type);

  /*!
   * \brief Sort the output data for each grid node into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SortOutputData(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Sort the output data for each grid node into a linear partitioning across all processors (DG-FEM solver).
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SortOutputData_FEM(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Sort the surface output data for each grid node into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SortOutputData_Surface(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Sort the surface output data for each grid node into a linear partitioning across all processors (DG-FEM solver).
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SortOutputData_Surface_FEM(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Deallocate temporary memory needed for merging and writing connectivity in parallel.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void DeallocateConnectivity_Parallel(CConfig *config, CGeometry *geometry, bool surf_sol);
  
  /*!
   * \brief Deallocate temporary memory needed for merging and writing output data in parallel.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void DeallocateData_Parallel(CConfig *config, CGeometry *geometry);
  
  /*!
   * \brief Deallocate temporary memory needed for merging and writing output data in parallel.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void DeallocateSurfaceData_Parallel(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Merge the node coordinates of all inlet boundaries from all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void MergeInletCoordinates(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Write a template inlet profile file for all inlets for flow problems.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Solver container.
   */
  void Write_InletFile_Flow(CConfig *config, CGeometry *geometry, CSolver **solver);

  /*!
   * \brief Deallocate temporary memory needed for merging and writing inlet boundary coordinates.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void DeallocateInletCoordinates(CConfig *config, CGeometry *geometry);
  
  /*! 
   * \brief Create and write a CSV file with a slice of data.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] FlowSolution - Flow solution.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] val_iZone - Current zone number in the grid file.
   * \param[in] val_direction - Controls the slice direction (0 for constant x/vertical, 1 for constant y/horizontal.
   */
  void WriteCSV_Slice(CConfig *config, CGeometry *geometry, CSolver *FlowSolver, unsigned long iExtIter, unsigned short val_iZone, unsigned short val_direction);

};

#include "output_structure_legacy.inl"
