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

#include <fstream>
#include <cmath>
#include <time.h>
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>

#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../../Common/include/toolboxes/signal_processing_toolbox.hpp"

class CGeometry;
class CSolver;
class CFileWriter;
class CParallelDataSorter;

#include "COutputLegacy.hpp"

using namespace std;

/*! 
 * \class COutput
 * \brief Class for writing the flow, adjoint and linearized solver 
 *        solution (including the history solution, and parallel stuff).
 * \author F. Palacios, T. Economon, M. Colonno.
 */
class COutput {

protected:
  
  CParallelDataSorter* data_sorter;
  
  COutputLegacy *output_legacy;

  std::vector< std::vector<su2double> > Local_Data;

  vector<string> Variable_Names;

  su2double RhoRes_New, *RhoRes_Old;

  bool no_writing;

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
    TYPE_AUTO_RESIDUAL,         /*!< \brief Integer format. Example: 34 */
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
  std::vector<short>                            Offset_Cache;
  unsigned short                                Offset_Cache_Index;
  bool                                          Offset_Cache_Checked,
                                                Build_Offset_Cache;
  
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
  
  string VolumeFilename, SurfaceFilename, RestartFilename, HistoryFilename;
  
  unsigned long curr_TimeIter, curr_OuterIter, curr_InnerIter;
    
  su2double Cauchy_Value,  /*!< \brief Summed value of the convergence indicator. */
  Cauchy_Func;      /*!< \brief Current value of the convergence indicator at one iteration. */
  unsigned short Cauchy_Counter;  /*!< \brief Number of elements of the Cauchy serial. */
  su2double *Cauchy_Serie;      /*!< \brief Complete Cauchy serial. */
  su2double Old_Func,  /*!< \brief Old value of the objective function (the function which is monitored). */
  New_Func;      /*!< \brief Current value of the objective function (the function which is monitored). */
  bool Convergence,    /*!< \brief To indicate if the flow solver (direct, adjoint, or linearized) has converged or not. */
  Convergence_FSI,    /*!< \brief To indicate if the FSI problem has converged or not. */
  Convergence_FullMG;    /*!< \brief To indicate if the Full Multigrid has converged and it is necessary to add a new level. */
  su2double InitResidual;  /*!< \brief Initial value of the residual to evaluate the convergence level. */
  string Conv_Field;
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
   * \param[in] wrt - If <TRUE> prepares history file for writing.
   */
  void PreprocessHistoryOutput(CConfig *config, bool wrt = true);  
  
  /*!
   * \brief Preprocess the history output by setting the history fields and opening the history file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] wrt - If <TRUE> prepares history file for writing.
   */
  void PreprocessMultizoneHistoryOutput(COutput **output, CConfig **config, bool wrt = true);  
  
  
  /*!
   * \brief SetHistory_Output
   * \param geometry
   * \param solver_container
   * \param config
   * \param TimeIter
   * \param OuterIter
   * \param InnerIter
   */
  void SetHistory_Output(CGeometry *geometry, CSolver **solver_container, CConfig *config,unsigned long TimeIter, unsigned long OuterIter, unsigned long InnerIter);  
  
  
  /*!
   * \brief SetHistory_Output
   * \param geometry
   * \param solver_container
   * \param config
   * \param TimeIter
   * \param OuterIter
   * \param InnerIter
   */
  void SetHistory_Output(CGeometry *geometry, CSolver **solver_container, CConfig *config);  
  
  
  void SetMultizoneHistory_Output(COutput** output, CConfig **config, CConfig *driver_config, unsigned long TimeIter, unsigned long OuterIter);
  
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
  
  /*!
   * \brief Set the current iteration indices
   * \param[in] TimeIter  - Timer iteration index
   * \param[in] OuterIter - Outer iteration index
   * \param[in] InnerIter - Inner iteration index
   */
  inline void SetIteration(unsigned long TimeIter, unsigned long OuterIter, unsigned long InnerIter){
    curr_TimeIter  = TimeIter;
    curr_OuterIter = OuterIter;
    curr_InnerIter = InnerIter;
  }
  
  su2double GetHistoryFieldValue(string field){
    return HistoryOutput_Map[field].Value;
  }
  
  vector<HistoryOutputField> GetHistoryGroup(string groupname){
    vector<HistoryOutputField> HistoryGroup;  
    for (unsigned short iField = 0; iField < HistoryOutput_Map.size(); iField++){   
      if (HistoryOutput_Map[HistoryOutput_List[iField]].OutputGroup == groupname){
        HistoryGroup.push_back((HistoryOutput_Map[HistoryOutput_List[iField]]));
      }
    }
    return HistoryGroup;
  }
  
  vector<string> GetHistoryOutput_List(){
    return HistoryOutput_List;
  }
  
  map<string, HistoryOutputField> GetHistoryFields(){
    return HistoryOutput_Map;
  }
  
  bool Convergence_Monitoring(CConfig *config, unsigned long Iteration);

  bool GetConvergence() {return Convergence;}

  void SetConvergence(bool conv) {Convergence = conv;}

protected:
  
  /*----------------------------- Protected member functions ----------------------------*/  

  /*!
   * \brief Allocates the appropriate file writer based on the chosen format.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] sorter - The parallel file sorter.
   * \param[out] filewriter - The allocated filewriter.
   * \param[in] format - The output format.
   */
  void SetFileWriter(CConfig *config, CGeometry *geomery, CParallelDataSorter *sorter, CFileWriter *&filewriter, unsigned short format);  
  
  
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
  void SetVolumeOutputValue(string name, unsigned long iPoint, su2double value);
  
  /*!
   * \brief CheckOffsetCache
   */
  void CheckOffsetCache();
  
  /*!
   * \brief CheckHistoryOutput
   */
  void CheckHistoryOutput();
  
  void PrepareHistoryFile(CConfig *config);

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
  void Postprocess_HistoryData(CConfig *config);

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
  inline bool PrintOutput(unsigned long iIter, unsigned long iFreq) {
    if (iFreq == 0.0){
      return false;
    }
    
    return (iIter % iFreq == 0);
  } 
  
  /*--------------------------------- Virtual functions ---------------------------------------- */
  
  /*!
   * \brief Determines if the history file output.
   * \param[in] config - Definition of the particular problem.
   */
  virtual bool WriteHistoryFile_Output(CConfig *config);
  
  /*!
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  virtual bool WriteScreen_Header(CConfig *config);
  
  /*!
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  virtual bool WriteScreen_Output(CConfig *config);

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
  inline virtual bool SetUpdate_Averages(CConfig *config){return false;}
  
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
  inline virtual void LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {}

  /*!
   * \brief Load the output data to the containers in each subclass
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void LoadMultizoneHistoryData(COutput **output, CConfig **config) {}  
  
  /*!
   * \brief SetHistoryOutputFields
   * \param config
   */
  inline virtual void SetHistoryOutputFields(CConfig *config) {}
  
  /*!
   * \brief SetMultizoneHistoryOutputFields
   * \param output
   * \param solver
   * \param config
   */
  inline virtual void SetMultizoneHistoryOutputFields(COutput **output, CConfig **config) {}

  inline void SetCommonHistoryFields(CConfig *config);
  
  inline void LoadCommonHistoryData(CConfig *config);
};

