/*!
 * \file output.hpp
 * \brief Headers of the main subroutines for generating the file outputs.
 *        The subroutines and functions are in the <i>output_structure.cpp</i> file.
 * \author F. Palacios, T. Economon, M. Colonno
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
#include <sstream>
#include <iomanip>
#include <limits>

#include "../solver_structure.hpp"
#include "../integration_structure.hpp"
#include "../../../Common/include/geometry_structure.hpp"
#include "../../../Common/include/fem_geometry_structure.hpp"
#include "../../../Common/include/fem_standard_element.hpp"
#include "../../../Common/include/config_structure.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../../Common/include/toolboxes/signal_processing_toolbox.hpp"

#include "output_structure_legacy.hpp"

using namespace std;

/*! 
 * \class COutput
 * \brief Class for writing the flow, adjoint and linearized solver 
 *        solution (including the history solution, and parallel stuff).
 * \author F. Palacios, T. Economon, M. Colonno.
 */
class COutput {

protected:
  
  COutputLegacy *output_legacy;

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
  map<unsigned long,unsigned long> Global2Renumber, 
                                   Renumber2Global;
  
  unsigned long nGlobalPoint_Sort;
  unsigned long nLocalPoint_Sort;
  unsigned long nPoint_Restart;
  int *Local_Halo_Sort;

  unsigned long *beg_node;
  unsigned long *end_node;

  unsigned long *nPoint_Lin;
  unsigned long *nPoint_Cum;
  
  std::vector< std::vector<su2double> > Local_Data;
  su2double **Local_Data_Copy;      // Local data copy for cte. lift mode
  su2double **Parallel_Data;        // node i (x, y, z) = (Coords[0][i], Coords[1][i], Coords[2][i])
  su2double **Parallel_Surf_Data;   // node i (x, y, z) = (Coords[0][i], Coords[1][i], Coords[2][i])
  vector<string> Variable_Names;
  int* Local_Halo;

  su2double **Data;
  unsigned short nVar_Consv, nVar_Total, nVar_Extra, nZones;
  bool wrote_surf_file, wrote_CGNS_base, wrote_Tecplot_base, wrote_Paraview_base;
  unsigned short wrote_base_file;
  su2double RhoRes_New, *RhoRes_Old;
  int cgns_base, cgns_zone, cgns_base_results, cgns_zone_results;
  
  su2double Sum_Total_RadialDistortion, Sum_Total_CircumferentialDistortion; // Add all the distortion to compute a run average.

  unsigned long nMarker_InletFile;       /*!< \brief Counter for total number of inlet boundaries written to inlet profile file. */
  vector<string> Marker_Tags_InletFile;   /*!< \brief Marker tags for the strings of the markers in the inlet profile file. */
  unsigned long *nRow_InletFile;         /*!< \brief Counters for the number of points per marker in the inlet profile file. */
  unsigned long *nRowCum_InletFile;      /*!< \brief Counters for the number of points per marker in cumulative storage format in the inlet profile file. */
  su2double **InletCoords;  /*!< \brief Data structure for holding the merged inlet boundary coordinates from all ranks. */

protected:

  int rank, 	  /*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */

  unsigned short nDim;
  
  unsigned short GlobalField_Counter; /*!< \brief Number of fields in the volume output */ 

  unsigned short field_width;         /*!< \brief Width of each column for the screen output (hardcoded for now) */
  
  /** \brief Enum to identify the screen output format. */
  enum ScreenOutputFormat {           
    FORMAT_INTEGER,         /*!< \brief Integer format. Example: 34 */
    FORMAT_FIXED,           /*!< \brief Format with fixed precision for floating point values. Example: 344.54  */
    FORMAT_SCIENTIFIC       /*!< \brief Scientific format for floating point values. Example: 3.4454E02 */  
  };
  
  enum HistoryFieldType {           
    TYPE_RESIDUAL,         /*!< \brief Integer format. Example: 34 */
    TYPE_COEFFICIENT,           /*!< \brief Format with fixed precision for floating point values. Example: 344.54  */
    TYPE_DEFAULT       /*!< \brief Scientific format for floating point values. Example: 3.4454E02 */  
  };
  
  string HistorySep;                  /*!< \brief Character which separates values in the history file */
  
  /** \brief Structure to store information for a history output field.
   * 
   *  The stored information is printed to the history file and to screen. 
   * Each individual instance represents a single field (i.e. column) in the history file or on screen.   
   */
  struct HistoryOutputField {
    string              FieldName;    /*!< \brief The name of the field, i.e. the name that is printed in the screen or file header.*/
    su2double           Value;        /*!< \brief The value of the field. */
    unsigned short      ScreenFormat; /*!< \brief The format that is used to print this value to screen. */
    string              OutputGroup;  /*!< \brief The group this field belongs to. */
    unsigned short      FieldType;
    HistoryOutputField() {}           /*!< \brief Default constructor. */
    HistoryOutputField(string fieldname, unsigned short screenformat, string historyoutputgroup, unsigned short fieldtype):
      FieldName(fieldname), Value(0.0), ScreenFormat(screenformat), OutputGroup(historyoutputgroup), FieldType(fieldtype){}
  };
  
  /** \brief Structure to store information for a volume output field.
   * 
   *  The stored information is used to create the volume solution file.   
   */
  struct VolumeOutputField {
    string FieldName;          /*!< \brief The name of the field, i.e. the name that is printed in the file header.*/
    int    Offset;             /*!< \brief This value identifies the position of the values of this field at each node in the ::Local_Data array. */
    string OutputGroup;        /*!< \brief The group this field belongs to. */
    VolumeOutputField () {}    /*!< \brief Default constructor. */
    VolumeOutputField(string fieldname, int offset, string volumeoutputgroup):
      FieldName(fieldname), Offset(offset), OutputGroup(volumeoutputgroup){}
  };

  std::map<string, HistoryOutputField >         HistoryOutput_Map;    /*!< \brief Associative map to access data stored in the history output fields by a string identifier. */
  std::vector<string>                           HistoryOutput_List;   /*!< \brief Vector that contains the keys of the HistoryOutput_Map in the order of their insertion. */
  std::map<string, vector<HistoryOutputField> > HistoryOutputPerSurface_Map; /*!< \brief Associative map to access data stored in the history per surface output fields by a string identifier. */
  std::vector<string>                           HistoryOutputPerSurface_List;  /*!< \brief Vector that contains the keys of the HistoryOutputPerSurface_Map in the order of their insertion. */

  std::map<string, VolumeOutputField >          VolumeOutput_Map;
  std::vector<string>                           VolumeOutput_List;
  
  std::vector<string> RequestedHistoryFields;
  unsigned short nRequestedHistoryFields;
  std::vector<string> RequestedScreenFields;
  unsigned short nRequestedScreenFields;
  std::vector<string> RequestedVolumeFields;
  unsigned short nRequestedVolumeFields;
  char char_histfile[200];

  ofstream HistFile;
  
  PrintingToolbox::CTablePrinter* ConvergenceTable;
  PrintingToolbox::CTablePrinter* MultiZoneHeaderTable;
  
  std::string MultiZoneHeaderString;
  
  std::map<string, su2double> Init_Residuals;
  
  map<string, Signal_Processing::RunningAverage> RunningAverages;
  
  bool multizone, grid_movement, fem_output;
  
  string VolumeFilename, SurfaceFilename, RestartFilename;
  
  
  
  
public:
  
  /*----------------------------- Public member functions ----------------------------*/
  
  /*! 
   * \brief Constructor of the class. 
   */
  COutput(CConfig *config);
  
  /*! 
   * \brief Destructor of the class. 
   */
  virtual ~COutput(void);  

  /*!
   * \brief Preprocess the volume output by setting the requested volume output fields.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void PreprocessVolumeOutput(CConfig *config, CGeometry *geometry);   
  
  /*!
   * \brief Load the data from the solvers into the local data array and sort it for the linear partitioning.
   * 
   * After calling this method the data is distributed to all processors based on a linear partition 
   * and is ready to be written in parallel to file using the SetVolume_Output or SetSurface_Output routines.
   * 
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - The container holding all solution data.
   */
  void Load_Data(CGeometry *geometry, CConfig *config, CSolver **solver_container);
  
  /*!
   * \brief Write the linear partitioned volume data in parallel to file. ::Load_Data() has to be called before!
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] format - The data format of the output files.
   */
  void SetVolume_Output(CGeometry *geometry, CConfig *config, unsigned short format);
  
  /*!
   * \brief Write the linear partitioned surface data in parallel to file. ::Load_Data() has to be called before!
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] format - The data format of the output files.
   */
  void SetSurface_Output(CGeometry *geometry, CConfig *config, unsigned short format);
  
  /*!
   * \brief Deallocate temporary memory needed for merging and writing output data in parallel.
   */
  void DeallocateData_Parallel();
    
  /*!
   * \brief Preprocess the history output by setting the history fields and opening the history file.
   * \param[in] config - Definition of the particular problem.
   */
  void PreprocessHistoryOutput(CConfig *config);  
  
  /*! 
   * \brief Write the history file and the convergence on the screen for serial computations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] integration - Generic subroutines for space integration, time integration, and monitoring.
   * \param[in] iExtIter - Current external (time) iteration.
   * \param[in] timeused - Current number of clock tick in the computation (related with total time).
   * \param[in] val_nZone - iZone index.
   */
  void SetConvHistory_Body(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
                              CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst);  
  
  /*! 
   * \brief Returns a pointer to the legacy output class needed for some old driver implementations.
   * \return - Returns a pointer to the legacy output class.
   */
  inline COutputLegacy* GetLegacyOutput(){return output_legacy;}

  /*!
   * \brief Sets the volume output filename
   * \param[in] filename - the new filename
   */
  inline void SetVolume_Filename(string filename) {VolumeFilename = filename;}
  
  /*!
   * \brief Sets the surface output filename
   * \param[in] filename - the new filename
   */
  inline void SetSurface_Filename(string filename) {SurfaceFilename = filename;}
  
  /*!
   * \brief Returns the current volume filename
   * \return - The current volume filename
   */
  inline string GetVolume_Filename() {return VolumeFilename;}
  
  /*!
   * \brief Returns the current surface filename
   * \return - The current surface filename
   */
  inline string GetSurface_Filename() {return SurfaceFilename;}
  
  
  /*!
   * \brief Sets the restart filename
   * \param[in] filename - the new filename
   */
  inline void SetRestart_Filename(string filename) {RestartFilename = filename;}
    
  /*!
   * \brief Returns the current restart filename
   * \return - The current restart filename
   */
  inline string GetRestart_Filename() {return RestartFilename;}
  
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
  
protected:
  
  /*----------------------------- Protected member functions ----------------------------*/  

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
  void WriteRestart_Parallel_ASCII(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Write a native SU2 restart file (binary) in parallel.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Flow, adjoint or linearized solution.
   * \param[in] val_iZone - iZone index.
   */
  void WriteRestart_Parallel_Binary(CConfig *config, CGeometry *geometry);

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
  
  /*!
   * \brief Write a Paraview ASCII solution file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - Current zone.
   * \param[in] val_nZone - Total number of zones.
   */
  void SetParaview_ASCII(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol);

  /*!
   * \brief Write a Paraview ASCII solution file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - Current zone.
   * \param[in] val_nZone - Total number of zones.
   */
  void SetParaview_MeshASCII(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol, bool new_file);

  /*!
   * \brief Write a Paraview ASCII solution file with parallel output.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - Current zone.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] surf_sol - Flag controlling whether this is a volume or surface file.
   */
  void WriteParaViewASCII_Parallel(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol);

  /*!
   * \brief Write a Paraview binary solution file with parallel output.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - Current zone.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] surf_sol - Flag controlling whether this is a volume or surface file.
   */
  void WriteParaViewBinary_Parallel(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol);
  
  /*!
   * \brief Write a Tecplot ASCII solution file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - Current zone.
   * \param[in] val_nZone - Total number of zones.
   */
  void SetTecplotASCII(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol);
  
  /*!
   * \brief Write the nodal coordinates and connectivity to a Tecplot ASCII mesh file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetTecplot_MeshASCII(CConfig *config, CGeometry *geometry, bool surf_sol, bool new_file);

  /*!
   * \brief Write the nodal coordinates and connectivity to a stl ASCII mesh file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetSTL_MeshASCII(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Write the nodal coordinates and connectivity to a n3d ASCII mesh file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SetCSV_MeshASCII(CConfig *config, CGeometry *geometry);
  
  /*!
   * \brief Write the nodal coordinates and connectivity to a n3d ASCII mesh file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetTecplotASCII_Mesh(CConfig *config, CGeometry *geometry, unsigned short val_iZone, bool surf_sol, bool new_file);

  /*!
   * \brief Write the solution data and connectivity to a Tecplot ASCII mesh file in parallel.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - Current zone.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] surf_sol - Flag controlling whether this is a volume or surface file.
   */
  void WriteTecplotASCII_Parallel(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol);

  /*!
   * \brief Write a Tecplot binary solution file with parallel output.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - Current zone.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] surf_sol - Flag controlling whether this is a volume or surface file.
   */
  void WriteTecplotBinary_Parallel(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol);
  
  /*!
   * \brief Write the solution data and connectivity to a CSV file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void WriteSurface_CSV(CConfig *config, CGeometry *geometry);

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
  void SetSU2_MeshASCII(CConfig *config, CGeometry *geometry);
  
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

  /*!
   * \brief Write a file with the adjoint sensitivities projected onto each surface node.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - Current zone.
   * \param[in] val_nZone - Total number of zones.
   */
  void WriteProjectedSensitivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone);
  
  /*!
   * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetTecplotBinary_DomainMesh(CConfig *config, CGeometry *geometry, unsigned short val_iZone);
  
  /*!
   * \brief Write the coordinates and connectivity to a Tecplot binary surface mesh file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetTecplotBinary_SurfaceMesh(CConfig *config, CGeometry *geometry, unsigned short val_iZone);
  
  /*!
   * \brief Write solution data to a Tecplot binary volume solution file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetTecplotBinary_DomainSolution(CConfig *config, CGeometry *geometry, unsigned short val_iZone);

  /*!
   * \brief Write solution data to a Tecplot binary surface solution file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetTecplotBinary_SurfaceSolution(CConfig *config, CGeometry *geometry, unsigned short val_iZone);
  
  /*!
   * \brief Write a Tecplot ASCII solution file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - Current zone.
   * \param[in] val_nZone - Total number of zones.
   */
  void SetFieldViewASCII(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone);
  
  /*!
   * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetFieldViewASCII_Mesh(CConfig *config, CGeometry *geometry);
  
  /*!
   * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetFieldViewBinary_Mesh(CConfig *config, CGeometry *geometry);
  
  /*!
   * \brief Write solution data to a Tecplot binary volume solution file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_iZone - iZone index.
   */
  void SetFieldViewBinary(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone);
  
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
   * \brief Write the sensitivity (including mesh sensitivity) computed with the discrete adjoint method
   *  on the surface and in the volume to a file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Number of Zones.
   */
  void SetSensitivity_Files(CGeometry ***geometry, CConfig **config, unsigned short val_nZone);

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
   * \param[in] surf - boolean controlling whether surface <TRUE> or volume connectivity <FALSE> should be sorted.
   * \param[in] val_sort - boolean controlling whether the elements are sorted or simply loaded by their owning rank.
   */
  void SortConnectivity(CConfig *config, CGeometry *geometry, bool surf, bool val_sort);

  /*!
   * \brief Sort the connectivities (volume and surface) into data structures used for output file writing (DG-FEM solver).
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] surf - boolean controlling whether surface <TRUE> or volume connectivity <FALSE> should be sorted.
   */
  void SortConnectivity_FEM(CConfig *config, CGeometry *geometry, bool surf);

  /*!
   * \brief Sort the connectivity for a single volume element type into a linear partitioning across all processors.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Elem_Type - VTK index of the element type being merged.
   * \param[in] val_sort - boolean controlling whether the elements are sorted or simply loaded by their owning rank.
   */
  void SortVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type, bool val_sort);

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
   * \param[in] surf_sol - if <TRUE>, surface connectivity is deallocated, otherwise the volume connectivity.
   */
  void DeallocateConnectivity_Parallel(bool surf_sol);
  
  /*!
   * \brief Deallocate temporary memory needed for merging and writing output data in parallel.
   */
  void DeallocateSurfaceData_Parallel();

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
  
  /*!
   * \brief Set the history file header
   * \param[in] config - Definition of the particular problem.
   */
  void SetHistoryFile_Header(CConfig *config);

  /*!
   * \brief Write the history file output
   * \param[in] config - Definition of the particular problem.
   */
  void SetHistoryFile_Output(CConfig *config);

  /*!
   * \brief Write the screen header.
   * \param[in] config - Definition of the particular problem.
   */
  void SetScreen_Header(CConfig *config);


  /*!
   * \brief Write the screen output.
   * \param[in] config - Definition of the particular problem.
   */
  void SetScreen_Output(CConfig *config);

  /*!
   * \brief Add a new field to the history output.
   * \param[in] name - Name for referencing it (in the config file and in the code).
   * \param[in] field_name - Header that is printed on screen and in the history file.
   * \param[in] format - The screen printing format (::ScreenOutputFormat).
   * \param[in] groupname - The name of the group this field belongs to.
   * \param[in] field_type - The type of the field (::HistoryFieldType).
   */
  inline void AddHistoryOutput(string name, string field_name, unsigned short format, string groupname, unsigned short field_type = TYPE_DEFAULT ){
    HistoryOutput_Map[name] = HistoryOutputField(field_name, format, groupname, field_type);
    HistoryOutput_List.push_back(name);
  }
  
  /*!
   * \brief Set the value of a history output field
   * \param[in] name - Name of the field.
   * \param[in] value - The new value of this field.
   */
  inline void SetHistoryOutputValue(string name, su2double value){
    if (HistoryOutput_Map.count(name) > 0){
      HistoryOutput_Map[name].Value = value;
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
    }
  }
  
  /*!
   * \brief Add a new field per surface marker to the history output.
   * \param[in] name - Name for referencing it (in the config file and in the code).
   * \param[in] field_name - Header that is printed on screen and in the history file.
   * \param[in] format - The screen printing format (::ScreenOutputFormat).
   * \param[in] groupname - The name of the group this field belongs to.
   * \param[in] marker_names - A list of markers. For every marker in this list a new field is created with "field_name + _marker_names[i]".
   * \param[in] field_type - The type of the field (::HistoryFieldType).
   */
  inline void AddHistoryOutputPerSurface(string name, string field_name, unsigned short format, string groupname, vector<string> marker_names, unsigned short field_type = TYPE_DEFAULT){
    if (marker_names.size() != 0){
      HistoryOutputPerSurface_List.push_back(name);
      for (unsigned short i = 0; i < marker_names.size(); i++){
        HistoryOutputPerSurface_Map[name].push_back(HistoryOutputField(field_name+"("+marker_names[i]+")", format, groupname, field_type));
      }
    }
  }
  
  /*!
   * \brief Set the value of a history output field for a specific surface marker
   * \param[in] name - Name of the field.
   * \param[in] value - The new value of this field.
   * \param[in] iMarker - The index of the marker.
   */
  inline void SetHistoryOutputPerSurfaceValue(string name, su2double value, unsigned short iMarker){
    if (HistoryOutputPerSurface_Map.count(name) > 0){
      HistoryOutputPerSurface_Map[name][iMarker].Value = value;
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
    }
  }
  
  /*!
   * \brief Add a new field to the volume output.
   * \param[in] name - Name for referencing it (in the config file and in the code).
   * \param[in] field_name - Header that is printed in the output files.
   * \param[in] groupname - The name of the group this field belongs to.
   */
  inline void AddVolumeOutput(string name, string field_name, string groupname){
    VolumeOutput_Map[name] = VolumeOutputField(field_name, -1, groupname);
    VolumeOutput_List.push_back(name);
  }
  
  /*!
   * \brief Set the value of a volume output field
   * \param[in] name - Name of the field.
   * \param[in] value - The new value of this field.
   */
  inline void SetVolumeOutputValue(string name, unsigned long iPoint, su2double value){
    if (VolumeOutput_Map.count(name) > 0){
      if (VolumeOutput_Map[name].Offset != -1){
        Local_Data[iPoint][VolumeOutput_Map[name].Offset] = value;
      }
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);    
    }
  }
  

  /*!
   * \brief Load up the values of the requested volume fields into ::Local_Data array. 
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - The container holding all solution data.
   */
  void CollectVolumeData(CConfig* config, CGeometry* geometry, CSolver** solver);
      
  /*!
   * \brief Postprocess_HistoryData
   * \param[in] config - Definition of the particular problem.
   * \param[in] dualtime - TODO: REMOVE PARAMETER
   */
  void Postprocess_HistoryData(CConfig *config, bool dualtime);

  /*!
   * \brief Postprocess_HistoryFields
   * \param[in] config - Definition of the particular problem.
   */
  void Postprocess_HistoryFields(CConfig *config);
  
  /*!
   * \brief Create and write the file with the FSI convergence history.
   * \param[in] iIter - Current iteration.
   * \param[in] iFreq - Frequency of output printing.
   */
  inline bool PrintOutput(unsigned long iIter, unsigned long iFreq) {return (iIter % iFreq == 0);} 
  
  /*--------------------------------- Virtual functions ---------------------------------------- */
  
  /*!
   * \brief Determines if the history file output.
   * \param[in] config - Definition of the particular problem.
   */
  virtual bool WriteHistoryFile_Output(CConfig *config, bool write_dualtime);
  
  /*!
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  virtual bool WriteScreen_Header(CConfig *config);
  
  /*!
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  virtual bool WriteScreen_Output(CConfig *config, bool write_dualtime);

  /*!
   * \brief LoadVolumeData
   * \param config
   * \param geometry
   * \param solver
   * \param iPoint
   */
  inline virtual void LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){}
  
  /*!
   * \brief LoadVolumeDataFEM
   * \param config
   * \param geometry
   * \param solver
   * \param iElem
   * \param index
   * \param dof
   */
  inline virtual void LoadVolumeDataFEM(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iElem, unsigned long index, unsigned short dof){}  
  
  /*!
   * \brief SetInit_Residuals
   * \param config
   * \return 
   */
  inline virtual bool SetInit_Residuals(CConfig *config) {return false;}
  
  /*!
   * \brief SetUpdate_Averages
   * \param config
   * \param dualtime
   * \return 
   */
  inline virtual bool SetUpdate_Averages(CConfig *config, bool dualtime){return false;}
  
  /*!
   * \brief LoadSurfaceData
   * \param config
   * \param geometry
   * \param solver
   * \param iPoint
   * \param iMarker
   * \param iVertex
   */
  inline virtual void LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){}
  
  /*!
   * \brief SetVolumeOutputFields
   * \param config
   */
  inline virtual void SetVolumeOutputFields(CConfig *config){}
  
  
  /*!
   * \brief Load the output data to the containers in each subclass
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void LoadHistoryData(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
                                      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) {}

  /*!
   * \brief SetHistoryOutputFields
   * \param config
   */
  inline virtual void SetHistoryOutputFields(CConfig *config) {}
  
};

