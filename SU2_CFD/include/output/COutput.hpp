/*!
 * \file COutput.hpp
 * \brief Headers of the output class.
 * \author T.Albring
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include <fstream>
#include <cmath>
#include <map>
#include <time.h>
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>

#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "tools/CWindowingTools.hpp"
#include "../../../Common/include/option_structure.hpp"

class CGeometry;
class CSolver;
class CFileWriter;
class CParallelDataSorter;
class CConfig;

using namespace std;

/*!
 * \class COutput
 * \brief Class for writing the convergence history and to write solution data to file.
 * \author T.Albring
 */
class COutput {
protected:

  /*----------------------------- General ----------------------------*/

  int rank,     /*!< \brief MPI Rank. */
  size;         /*!< \brief MPI Size. */

  unsigned short nDim;   /*!< \brief Physical Dimension */

  bool multiZone,       /*!< \brief Boolean to store whether we are running a multizone problem */
  gridMovement,         /*!< \brief Boolean to store whether we have grid movement enabled */
  femOutput;            /*!< \brief Boolean to store whether we should use the FEM routines */

  /*----------------------------- Screen and history output ----------------------------*/

  string historySep;              /*!< \brief Character which separates values in the history file */
  unsigned short fieldWidth;      /*!< \brief Width of each column for the screen output (hardcoded for now) */
  bool noWriting;                 /*!< \brief Boolean indicating whether a screen/history output should be written */
  unsigned long curTimeIter,      /*!< \brief Current value of the time iteration index */
  curAbsTimeIter,                 /*!< \brief Current value of the time iteration index */
  curOuterIter,                   /*!< \brief Current value of the outer iteration index */
  curInnerIter;                   /*!< \brief Current value of the inner iteration index */

  string historyFilename;   /*!< \brief The history filename*/
  char char_histfile[200];  /*! \brief Temporary variable to store the history filename */
  ofstream histFile;        /*! \brief Output file stream for the history */

  /** \brief Enum to identify the screen output format. */
  enum class ScreenOutputFormat {
    INTEGER,         /*!< \brief Integer format. Example: 34 */
    FIXED,           /*!< \brief Format with fixed precision for floating point values. Example: 344.54  */
    SCIENTIFIC       /*!< \brief Scientific format for floating point values. Example: 3.4454E02 */
  };

  /** \brief Enum to identify the screen/history field type. */
  enum class HistoryFieldType {
    RESIDUAL,         /*!< \brief A user-defined residual field type*/
    AUTO_RESIDUAL,    /*!< \brief An automatically generated residual field type */
    COEFFICIENT,      /*!< \brief User defined coefficient field type  */
    AUTO_COEFFICIENT, /*!< \brief Automatically generated coefficient field type  */
    DEFAULT           /*!< \brief Default field type */
  };

  /** \brief Structure to store information for a history output field.
   *
   *  The stored information is printed to the history file and to screen.
   * Each individual instance represents a single field (i.e. column) in the history file or on screen.
   */
  struct HistoryOutputField {
    /*! \brief The name of the field, i.e. the name that is printed in the screen or file header.*/
    string              fieldName = "";
    /*! \brief The value of the field. */
    su2double           value = 0.0;
    /*! \brief The format that is used to print this value to screen. */
    ScreenOutputFormat  screenFormat = ScreenOutputFormat::FIXED;
    /*! \brief The group this field belongs to. */
    string              outputGroup  ="";
    /*! \brief The field type*/
    HistoryFieldType     fieldType = HistoryFieldType::DEFAULT;
    /*! \brief String containing the description of the field */
    string              description = "";
    /*! \brief Default constructor. */
    HistoryOutputField() {}
    /*! \brief Constructor to initialize all members. */
    HistoryOutputField(string fieldName_, ScreenOutputFormat screenFormat_, string OutputGroup_,
                       HistoryFieldType fieldType_, string description_):
      fieldName(std::move(fieldName_)), value(0.0), screenFormat(screenFormat_),
      outputGroup(std::move(OutputGroup_)), fieldType(fieldType_), description(std::move(description_)){}
  };

  /*! \brief Associative map to access data stored in the history output fields by a string identifier. */
  std::map<string, HistoryOutputField >         historyOutput_Map;
  /*! \brief Vector that contains the keys of the ::historyOutput_Map in the order of their insertion. */
  std::vector<string>                           historyOutput_List;
  /*! \brief Associative map to access data stored in the history per surface output fields by a string identifier. */
  std::map<string, vector<HistoryOutputField> > historyOutputPerSurface_Map;
  /*! \brief Vector that contains the keys of the ::historyOutputPerSurface_Map in the order of their insertion. */
  std::vector<string>                           historyOutputPerSurface_List;

  /*! \brief Requested history field names in the config file. */
  std::vector<string> requestedHistoryFields;
  /*! \brief Number of requested history field names in the config file. */
  unsigned short nRequestedHistoryFields;
  /*! \brief Requested screen field names in the config file. */
  std::vector<string> requestedScreenFields;
  /*! \brief Number of requested screen field names in the config file. */
  unsigned short nRequestedScreenFields;

  PrintingToolbox::CTablePrinter* convergenceTable;     //!< Convergence  output table structure
  PrintingToolbox::CTablePrinter* multiZoneHeaderTable; //!< Multizone header output structure
  PrintingToolbox::CTablePrinter* historyFileTable;     //!< Table structure for writing to history file
  PrintingToolbox::CTablePrinter* fileWritingTable;     //!< File writing header
  std::string multiZoneHeaderString;                    //!< Multizone header string
  bool headerNeeded;                                    //!< Boolean that stores whether a screen header is needed

  //! Structure to store the value of the running averages
  map<string, CWindowedAverage> windowedTimeAverages;

  //! Structure to store the value initial residuals for relative residual computation
  std::map<string, su2double> initialResiduals;

   /*----------------------------- Volume output ----------------------------*/

   CParallelDataSorter* volumeDataSorter;    //!< Volume data sorter
   CParallelDataSorter* surfaceDataSorter;   //!< Surface data sorter

   vector<string> volumeFieldNames;     //!< Vector containing the volume field names
   unsigned short nVolumeFields;        /*!< \brief Number of fields in the volume output */

   string volumeFilename,               //!< Volume output filename
   surfaceFilename,                     //!< Surface output filename
   restartFilename;                     //!< Restart output filename

  /** \brief Structure to store information for a volume output field.
   *
   *  The stored information is used to create the volume solution file.
   */
  struct VolumeOutputField {
    /*! \brief The name of the field, i.e. the name that is printed in the file header.*/
    string fieldName;
    /*! \brief This value identifies the position of the values of this field at each node in the ::Local_Data array. */
    short  offset;
    /*! \brief The group this field belongs to. */
    string outputGroup;
    /*! \brief String containing the description of the field */
    string description;
    /*! \brief Default constructor. */
    VolumeOutputField () {}
    /*! \brief Constructor to initialize all members. */
    VolumeOutputField(string fieldName_, int offset_, string volumeOutputGroup_, string description_):
      fieldName(std::move(fieldName_)), offset(std::move(offset_)),
      outputGroup(std::move(volumeOutputGroup_)), description(std::move(description_)){}
  };

  /*! \brief Associative map to access data stored in the volume output fields by a string identifier. */
  std::map<string, VolumeOutputField >          volumeOutput_Map;
  /*! \brief Vector that contains the keys of the ::volumeOutput_Map in the order of their insertion. */
  std::vector<string>                           volumeOutput_List;

  /*! \brief Vector to cache the positions of the field in the data array */
  std::vector<short>                            fieldIndexCache;
  /*! \brief Current value of the cache index */
  unsigned short                                cachePosition;
  /*! \brief Boolean to store whether the field index cache should be build. */
  bool                                          buildFieldIndexCache;
  /*! \brief Vector to cache the positions of the field in the data array */
  std::vector<short>                            fieldGetIndexCache;
  /*! \brief Current value of the cache index */
  unsigned short                                curGetFieldIndex;

  /*! \brief Requested volume field names in the config file. */
  std::vector<string> requestedVolumeFields;
  /*! \brief Number of requested volume field names in the config file. */
  unsigned short nRequestedVolumeFields;

  /*----------------------------- Convergence monitoring ----------------------------*/

  su2double cauchyValue,         /*!< \brief Summed value of the convergence indicator. */
  cauchyFunc;                    /*!< \brief Current value of the convergence indicator at one iteration. */
  unsigned short Cauchy_Counter; /*!< \brief Number of elements of the Cauchy serial. */
  vector<vector<su2double> > cauchySerie;        /*!< \brief Complete Cauchy serial. */
  unsigned long nCauchy_Elems;   /*!< \brief Total number of cauchy elems to monitor */
  su2double cauchyEps;           /*!< \brief Defines the threshold when to stop the solver. */
  su2double minLogResidual;      /*!< \brief Minimum value of the residual to reach */
  vector<su2double> oldFunc,     /*!< \brief Old value of the coefficient. */
  newFunc;                       /*!< \brief Current value of the coefficient. */
  bool convergence;              /*!< \brief To indicate if the solver has converged or not. */
  su2double initResidual;        /*!< \brief Initial value of the residual to evaluate the convergence level. */
  vector<string> convFields;     /*!< \brief Name of the field to be monitored for convergence. */

  /*----------------------------- Adaptive CFL ----------------------------*/     

  su2double rhoResNew,    /*!< New value of the residual for adaptive CFL routine. */
  rhoResOld;              /*!< Old value of the residual for adaptive CFL routine. */

  /*----------------------------- Time Convergence monitoring ----------------------------*/
  vector<string> wndConvFields;                /*!< \brief Name of the field to be monitored for convergence. */
  vector<vector<su2double> > WndCauchy_Serie;  /*!< \brief Complete Cauchy serial. */
  unsigned long nWndCauchy_Elems;              /*!< \brief Total number of cauchy elems to monitor */
  su2double wndCauchyEps;                      /*!< \brief Defines the threshold when to stop the solver. */

  vector<su2double> WndOld_Func;  /*!< \brief Old value of the objective function (the function which is monitored). */
  vector<su2double> WndNew_Func;  /*!< \brief Current value of the objective function (the function which is monitored). */
  su2double WndCauchy_Func;       /*!< \brief Current value of the convergence indicator at one iteration. */
  su2double WndCauchy_Value;      /*!< \brief Summed value of the convergence indicator. */
  bool TimeConvergence;   /*!< \brief To indicate, if the windowed time average of the time loop has converged*/

public:

  /*----------------------------- Public member functions ----------------------------*/

  /*!
   * \brief Constructor of the class.
   */
  COutput(CConfig *config, unsigned short nDim, bool femOutput);

  /*!
   * \brief Preprocess the volume output by setting the requested volume output fields.
   * \param[in] config - Definition of the particular problem.
   */
  void PreprocessVolumeOutput(CConfig *config);

  /*!
   * \brief Load the data from the solvers into the data sorters and sort it for the linear partitioning.
   *
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - The container holding all solution data.
   */
  void Load_Data(CGeometry *geometry, CConfig *config, CSolver **solver_container);

  /*!
   * \brief Preprocess the history output by setting the history fields and opening the history file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] wrt - If <TRUE> prepares history file for writing.
   */
  void PreprocessHistoryOutput(CConfig *config, bool wrt = true);

  /*!
   * \brief Preprocess the history output by setting the history fields and opening the history file.
   * \param[in] output - Container holding the output instances per zone.
   * \param[in] config - Definition of the particular problem per zone.
   * \param[in] wrt - If <TRUE> prepares history file for writing.
   */
  void PreprocessMultizoneHistoryOutput(COutput **output, CConfig **config, CConfig *driver_config, bool wrt = true);

  /*!
   * \brief Collects history data from the solvers, monitors the convergence and writes to screen and history file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] TimeIter - Value of the time iteration index
   * \param[in] OuterIter - Value of outer iteration index
   * \param[in] InnerIter - Value of the inner iteration index
   */
  void SetHistory_Output(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                         unsigned long TimeIter, unsigned long OuterIter, unsigned long InnerIter);

  /*!
   * \brief Collects history data from the solvers and monitors the convergence. Does not write to screen or file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetHistory_Output(CGeometry *geometry, CSolver **solver_container, CConfig *config);

  /*!
   *  Collects history data from the individual output per zone,
   *  monitors the convergence and writes to screen and history file.

   * \param[in] output - Container holding the output instances per zone.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem per zone.
   * \param[in] driver_config - Base definition of the particular problem.
   * \param[in] TimeIter - Value of the time iteration index
   * \param[in] OuterIter - Value of outer iteration index
   */
  void SetMultizoneHistory_Output(COutput** output, CConfig **config, CConfig *driver_config,
                                  unsigned long TimeIter, unsigned long OuterIter);

  /*!
   * \brief Sets the volume output filename
   * \param[in] filename - the new filename
   */
  inline void SetVolume_Filename(string filename) {volumeFilename = filename;}

  /*!
   * \brief Sets the surface output filename
   * \param[in] filename - the new filename
   */
  inline void SetSurface_Filename(string filename) {surfaceFilename = filename;}

  /*!
   * \brief Returns the current volume filename
   * \return - The current volume filename
   */
  inline string GetVolume_Filename() {return volumeFilename;}

  /*!
   * \brief Returns the current surface filename
   * \return - The current surface filename
   */
  inline string GetSurface_Filename() {return surfaceFilename;}

  /*!
   * \brief Sets the restart filename
   * \param[in] filename - the new filename
   */
  inline void SetRestart_Filename(string filename) {restartFilename = filename;}

  /*!
   * \brief Returns the current restart filename
   * \return - The current restart filename
   */
  inline string GetRestart_Filename() {return restartFilename;}

  /*!
   * \brief Set the current iteration indices
   * \param[in] TimeIter  - Timer iteration index
   * \param[in] OuterIter - Outer iteration index
   * \param[in] InnerIter - Inner iteration index
   */
  inline void SetIteration(unsigned long TimeIter, unsigned long OuterIter, unsigned long InnerIter){
    curTimeIter  = TimeIter;
    curOuterIter = OuterIter;
    curInnerIter = InnerIter;
  }

  /*!
   * \brief Get the value of particular history output field
   * \param[in] field - Name of the field
   * \return Value of the field
   */
  su2double GetHistoryFieldValue(string field){
    return historyOutput_Map[field].value;
  }

  su2double GetHistoryFieldValuePerSurface(string field, unsigned short iMarker){
    return historyOutputPerSurface_Map[field][iMarker].value;
  }

  /*!
   * \brief Get a vector with all output fields in a particular group
   * \param groupname - Name of the history group
   * \return Vector containing all output fields of a group
   */
  vector<HistoryOutputField> GetHistoryGroup(string groupname){
    vector<HistoryOutputField> HistoryGroup;
    for (unsigned short iField = 0; iField < historyOutput_Map.size(); iField++){
      if (historyOutput_Map[historyOutput_List[iField]].outputGroup == groupname){
        HistoryGroup.push_back((historyOutput_Map[historyOutput_List[iField]]));
      }
    }
    return HistoryGroup;
  }

  /*!
   * \brief Get the list of all output fields
   * \return Vector container all output fields
   */
  vector<string> GetHistoryOutput_List(){
    return historyOutput_List;
  }

  /*!
   * \brief Get the map containing all output fields
   * \return Map containing all output fields
   */
  map<string, HistoryOutputField> GetHistoryFields(){
    return historyOutput_Map;
  }

  /*!
   * \brief Monitor the convergence of an output field
   * \param[in] config - Definition of the particular problem.
   * \param[in] Iteration - Index of the current iteration.
   * \return Boolean indicating whether the problem is converged.
   */
  bool Convergence_Monitoring(CConfig *config, unsigned long Iteration);

  /*!
   * \brief Print a summary of the convergence to screen.
   */
  void PrintConvergenceSummary();

  /*!
   * \brief Get convergence of the problem.
   * \return Boolean indicating whether the problem is converged.
   */
  bool GetConvergence() {return convergence;}

  /*!
     * \brief  Monitor the time convergence of the specified windowed-time-averaged ouput
     * \param[in] config - Definition of the particular problem.
     * \param[in] Iteration - Index of the current iteration.
     * \return Boolean indicating whether the problem is converged.
     */
  bool MonitorTimeConvergence(CConfig *config, unsigned long Iteration);

  /*!
   * \brief Get convergence time convergence of the specified windowed-time-averaged ouput of the problem.
   * \return Boolean indicating whether the problem is converged.
   */
  bool GetTimeConvergence()const {return TimeConvergence;} /*! \brief Indicates, if the time loop is converged. COnvergence criterion: Windowed time average */


  /*!
   * \brief Set the value of the convergence flag.
   * \param[in] conv - New value of the convergence flag.
   */
  void SetConvergence(bool conv) {convergence = conv;}

  /*!
   * \brief Print a list of all history output fields to screen.
   */
  void PrintHistoryFields();

  /*!
   * \brief Print a list of all volume output fields to screen.
   */
  void PrintVolumeFields();

  /*!
   * \brief Loop through all requested output files and write the volume output data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] iter - The current time, outer or inner iteration index.
   * \param[in] force_writing - If <TRUE>, writing of output files is forced without checking the output frequency.
   * \return <TRUE> if output files have been written to disk.
   */
  bool SetResult_Files(CGeometry *geometry, CConfig *config, CSolver** solver_container,
                       unsigned long iter, bool force_writing = false);

  /*!
   * \brief Allocates the appropriate file writer based on the chosen format and writes sorted data to file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] format - The output format.
   * \param[in] fileName - The file name. If empty, the filenames are automatically determined.
   */
  void WriteToFile(CConfig *config, CGeometry *geomery, unsigned short format, string fileName = "");

protected:

  /*----------------------------- Protected member functions ----------------------------*/

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
   * \param[in] description - A description of the field.
   * \param[in] field_type - The type of the field (::HistoryFieldType).
   */
  inline void AddHistoryOutput(string name, string field_name, ScreenOutputFormat format,
                               string groupname, string description,
                               HistoryFieldType field_type = HistoryFieldType::DEFAULT ){
    historyOutput_Map[name] = HistoryOutputField(field_name, format, groupname, field_type, description);
    historyOutput_List.push_back(name);
  }

  /*!
   * \brief Set the value of a history output field
   * \param[in] name - Name of the field.
   * \param[in] value - The new value of this field.
   */
  inline void SetHistoryOutputValue(string name, su2double value){
    if (historyOutput_Map.count(name) > 0){
      historyOutput_Map[name].value = value;
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
  inline void AddHistoryOutputPerSurface(string name, string field_name, ScreenOutputFormat format,
                                         string groupname, vector<string> marker_names,
                                         HistoryFieldType field_type = HistoryFieldType::DEFAULT){
    if (marker_names.size() != 0){
      historyOutputPerSurface_List.push_back(name);
      for (unsigned short i = 0; i < marker_names.size(); i++){
        historyOutputPerSurface_Map[name].push_back(HistoryOutputField(field_name+"("+marker_names[i]+")", format, groupname, field_type, ""));
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
    if (historyOutputPerSurface_Map.count(name) > 0){
      historyOutputPerSurface_Map[name][iMarker].value = value;
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
    }
  }

  /*!
   * \brief Add a new field to the volume output.
   * \param[in] name - Name for referencing it (in the config file and in the code).
   * \param[in] field_name - Header that is printed in the output files.
   * \param[in] groupname - The name of the group this field belongs to.
   * \param[in] description - Description of the volume field.
   */
  inline void AddVolumeOutput(string name, string field_name, string groupname, string description){
    volumeOutput_Map[name] = VolumeOutputField(field_name, -1, groupname, description);
    volumeOutput_List.push_back(name);
  }


  /*!
   * \brief Set the value of a volume output field
   * \param[in] name - Name of the field.
   * \param[in] value - The new value of this field.
   */
  su2double GetVolumeOutputValue(string name, unsigned long iPoint);

  /*!
   * \brief Set the value of a volume output field
   * \param[in] name - Name of the field.
   * \param[in] value - The new value of this field.
   */
  void SetVolumeOutputValue(string name, unsigned long iPoint, su2double value);

  /*!
   * \brief Set the value of a volume output field
   * \param[in] name - Name of the field.
   * \param[in] value - The new value of this field.
   */
  void SetAvgVolumeOutputValue(string name, unsigned long iPoint, su2double value);

  /*!
   * \brief CheckHistoryOutput
   */
  void CheckHistoryOutput();

  /*!
   * \brief Open the history file and write the header.
   * \param[in] config - Definition of the particular problem.
   */
  void PrepareHistoryFile(CConfig *config);

  /*!
   * \brief Load up the values of the requested volume fields into ::Local_Data array.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - The container holding all solution data.
   */
  void LoadDataIntoSorter(CConfig* config, CGeometry* geometry, CSolver** solver);

  /*!
   * \brief Postprocess_HistoryData
   * \param[in] config - Definition of the particular problem.
   */
  void Postprocess_HistoryData(CConfig *config);

  /*!
   * \brief Postprocess_HistoryFields
   * \param[in] config - Definition of the particular problem.
   */
  void Postprocess_HistoryFields(CConfig *config);

  /*!
   * \brief Check whether we should print output.
   * \param[in] iIter - Current iteration.
   * \param[in] iFreq - Frequency of output printing.
   */
  inline bool PrintOutput(unsigned long iIter, unsigned long iFreq) {
    if (iFreq == 0){
      return false;
    }

    return (iIter % iFreq == 0);
  }

  /*!
   * \brief Set the history fields common for all solvers.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCommonHistoryFields(CConfig *config);

  /*!
   * \brief Load values of the history fields common for all solvers.
   * \param[in] config - Definition of the particular problem.
   */
  void LoadCommonHistoryData(CConfig *config);

  /*!
   * \brief Allocates the data sorters if necessary.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void AllocateDataSorters(CConfig *config, CGeometry *geometry);

  /*--------------------------------- Virtual functions ---------------------------------------- */
public:

  /*!
   * \brief Destructor of the class.
   */
  virtual ~COutput(void);

protected:

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
   * \brief Determines if the the volume output should be written.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Iter - Current iteration index.
   * \param[in] force_writing - boolean that forces writing of volume output
   */
  virtual bool WriteVolume_Output(CConfig *config, unsigned long Iter, bool force_writing);

  /*!
   * \brief Set the values of the volume output fields for a point.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The container holding all solution data.
   * \param[in] iPoint - Index of the point.
   */
  inline virtual void LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){}

  /*!
   * \brief Set the values of the volume output fields for a point.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The container holding all solution data.
   * \param[in] iElem - Index of the element.
   * \param[in] index - Index of the value.
   * \param[in] dof - Index of the local degree of freedom.
   */
  inline virtual void LoadVolumeDataFEM(CConfig *config, CGeometry *geometry, CSolver **solver,
                                        unsigned long iElem, unsigned long index, unsigned short dof){}

  /*!
   * \brief Check whether the base values for relative residuals should be initialized
   * \param[in] config - Definition of the particular problem.
   * \return <TRUE> if the residuals should be initialized.
   */
  inline virtual bool SetInit_Residuals(CConfig *config) {return false;}

  /*!
   * \brief Check whether the averaged values should be updated
   * \param[in] config - Definition of the particular problem.
   * \return <TRUE> averages should be updated.
   */
  inline virtual bool SetUpdate_Averages(CConfig *config){return false;}

  /*!
   * \brief Set the values of the volume output fields for a surface point.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The container holding all solution data.
   * \param[in] iPoint - Index of the point.
   * \param[in] iMarker - Index of the surface marker.
   * \param[in] iVertex - Index of the vertex on the marker.
   */
  inline virtual void LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver,
                                      unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){}

  /*!
   * \brief Set the available volume output fields
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetVolumeOutputFields(CConfig *config){}


  /*!
   * \brief Load the history output field values
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {}

  /*!
   * \brief Load the multizone history output field values
   * \param[in] output - Container holding the output instances per zone.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void LoadMultizoneHistoryData(COutput **output, CConfig **config) {}

  /*!
   * \brief Set the available history output fields
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetHistoryOutputFields(CConfig *config) {}

  /*!
   * \brief Set the available multizone history output fields
   * \param[in] output - Container holding the output instances per zone.
   * \param[in] config - Definition of the particular problem per zone.
   */
  inline virtual void SetMultizoneHistoryOutputFields(COutput **output, CConfig **config) {}

  /*!
   * \brief Write any additional files defined for the current solver.
   * \param[in] config - Definition of the particular problem per zone.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - The container holding all solution data.
   */
  inline virtual void WriteAdditionalFiles(CConfig *config, CGeometry* geometry, CSolver** solver_container){}

  /*!
   * \brief Write any additional output defined for the current solver.
   * \param[in] config - Definition of the particular problem per zone.
   */
  inline virtual void SetAdditionalScreenOutput(CConfig *config){}

};

