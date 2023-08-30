/*!
 * \file COutput.hpp
 * \brief Headers of the output class.
 * \author T.Albring
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

#pragma once

#include <fstream>
#include <cmath>
#include <map>
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>

#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "tools/CWindowingTools.hpp"
#include "../../../Common/include/option_structure.hpp"

/*--- AD workaround for a cmath function not defined in CoDi. ---*/
namespace mel {
namespace internal {
inline su2double hypot(const su2double& a, const su2double& b) {
  return sqrt(a*a + b*b);
}
}
}
#include "mel.hpp"

class CGeometry;
class CSolver;
class CFileWriter;
class CParallelDataSorter;
class CConfig;

using namespace std;

/*!
 * \brief Class for writing the convergence history and to write solution data to file.
 * \author T.Albring
 * \ingroup Output
 */
class COutput {
protected:

  /*----------------------------- General ----------------------------*/

  const int rank;     /*!< \brief MPI Rank. */
  const int size;     /*!< \brief MPI Size. */

  const unsigned short nDim;   /*!< \brief Physical Dimension */

  const bool multiZone;     /*!< \brief Boolean to store whether we are running a multizone problem */
  const bool gridMovement;  /*!< \brief Boolean to store whether we have grid movement enabled */
  const bool femOutput;     /*!< \brief Boolean to store whether we should use the FEM routines */
  const bool si_units;
  const bool us_units;

  /*----------------------------- Screen and history output ----------------------------*/

  const unsigned short fieldWidth = 12; /*!< \brief Width of each column for the screen output (hardcoded for now) */
  string historySep;              /*!< \brief Character which separates values in the history file */
  bool noWriting;                 /*!< \brief Boolean indicating whether a screen/history output should be written */
  unsigned long curTimeIter,      /*!< \brief Current value of the time iteration index */
  curAbsTimeIter,                 /*!< \brief Current value of the time iteration index */
  curOuterIter,                   /*!< \brief Current value of the outer iteration index */
  curInnerIter;                   /*!< \brief Current value of the inner iteration index */

  string historyFilename;   /*!< \brief The history filename*/
  ofstream histFile;        /*!< \brief Output file stream for the history */

  bool cauchyTimeConverged; /*! \brief: Flag indicating that solver is already converged. Needed for writing restart files. */

  /** \brief Enum to identify the screen output format. */
  enum class ScreenOutputFormat {
    INTEGER,         /*!< \brief Integer format. Example: 34 */
    FIXED,           /*!< \brief Format with fixed precision for floating point values. Example: 344.54  */
    SCIENTIFIC,      /*!< \brief Scientific format for floating point values. Example: 3.4454E02 */
    PERCENT          /*!< \brief Format with fixed precision for floating point values with a % signs. Example: 99.52% */
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

  /*! \brief Caches to avoid hashing the output maps to retrieve field values. */
  std::vector<const su2double*> requestedHistoryFieldCache;
  std::vector<const HistoryOutputField*> requestedScreenFieldCache;

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

  /*! \brief Struct to hold a parsed user-defined expression. */
  struct CustomHistoryOutput {
    mel::ExpressionTree<passivedouble> expression;
    /*--- Pointers to values in the history output maps, to avoid key lookup every time. ---*/
    std::vector<const su2double*> symbolValues;
    bool ready = false;

    su2double Eval() const {
      return mel::Eval<su2double>(expression, [&](int i) {return *symbolValues[i];});
    }
  };

  CustomHistoryOutput customObjFunc;  /*!< \brief User-defined expression for a custom objective. */

  /*! \brief Type of operation for custom outputs. */
  enum class OperationType { MACRO, FUNCTION, AREA_AVG, AREA_INT, MASSFLOW_AVG, MASSFLOW_INT, PROBE };

  /*! \brief Struct to hold a parsed custom output function. */
  struct CustomOutput {
    /*--- First level of parsing the syntax "name : type{func}[markers];". ---*/
    std::string name;
    OperationType type;
    std::string func;
    std::vector<std::string> markers;

    /*--- Second level, func into expression, and acceleration structures. ---*/
    mel::ExpressionTree<passivedouble> expression;
    std::vector<std::string> varSymbols;
    std::vector<unsigned short> markerIndices;
    static constexpr long PROBE_NOT_SETUP = -2;
    static constexpr long PROBE_NOT_OWNED = -1;
    long iPoint = PROBE_NOT_SETUP;

    /*--- The symbols (strings) are associated with an integer index for efficiency. For evaluation this index
     is passed to a functor that returns the value associated with the symbol. This functor is an input to "eval()"
     and needs to be generated on-the-fly for each point. The functor approach is more generic than a pointer, for
     example it allows wrapping the access to multiple solvers. The interpretation of these indices is dictated by
     the functor used in eval, for example indices may be established as 32 * solver_idx + variable_idx.
     The parts of the code that assign and interpret indices need to be in sync. ---*/
    std::vector<unsigned long> varIndices;

    /*--- Offset between varIndices of different solvers (see above). Power of 2 to make decoding faster. ---*/
    static constexpr unsigned long MAX_VARS_PER_SOLVER = 32;

    /*--- Arbitrary number to indicate that a string did not match a variable. ---*/
    static constexpr unsigned long NOT_A_VARIABLE = MAX_SOLS * MAX_VARS_PER_SOLVER;

    /*--- Other outputs can be referenced in expressions, e.g. to compute variance.
     We store pointers to the required outputs to speed-up access. ---*/
    std::vector<const su2double*> otherOutputs;

    /*--- For discrete adjoint we may need to skip some expressions because there is one output class
     for the primal solver and one for the discrete adjoint (each with different variables). ---*/
    bool skip = false;

    /*--- For evaluation, "vars" is a functor (i.e. has operator()) that returns the value of a variable at a given
     point. For example, it can be a wrapper to the primitives pointer, in which case varIndices needs to be setup
     with primitive indices. ---*/
    template <class Variables>
    su2double Eval(const Variables& vars) const {
      return mel::Eval<su2double>(expression, [&](int iSymbol) {return vars(varIndices[iSymbol]);});
    }
  };

  std::vector<CustomOutput> customOutputs;  /*!< \brief User-defined outputs. */

  /*----------------------------- Volume output ----------------------------*/

  CParallelDataSorter* volumeDataSorter;    //!< Volume data sorter
  CParallelDataSorter* surfaceDataSorter;   //!< Surface data sorter

  vector<string> volumeFieldNames;     //!< Vector containing the volume field names
  unsigned short nVolumeFields;        //!< Number of fields in the volume output

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
  COutput(const CConfig *config, unsigned short nDim, bool femOutput);

  /*!
   * \brief Preprocess the volume output by setting the requested volume output fields.
   * \param[in] config - Definition of the particular problem.
   */
  void PreprocessVolumeOutput(CConfig *config);

  /*!
   * \brief Load the data from the solvers into the data sorters and sort it for the linear partitioning.
   *
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver_container - The container holding all solution data.
   */
  void LoadData(CGeometry *geometry, CConfig *config, CSolver **solver_container);

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
   * \param[in] driver_config - Base definition of the particular problem.
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
  void SetHistoryOutput(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                         unsigned long TimeIter, unsigned long OuterIter, unsigned long InnerIter);

  /*!
   * \brief Collects history data from the solvers and monitors the convergence. Does not write to screen or file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetHistoryOutput(CGeometry *geometry, CSolver **solver_container, CConfig *config);

  /*!
   *  Collects history data from the individual output per zone,
   *  monitors the convergence and writes to screen and history file.

   * \param[in] output - Container holding the output instances per zone.
   * \param[in] config - Definition of the particular problem per zone.
   * \param[in] driver_config - Base definition of the particular problem.
   * \param[in] TimeIter - Value of the time iteration index
   * \param[in] OuterIter - Value of outer iteration index
   */
  void SetMultizoneHistoryOutput(COutput** output, CConfig **config, CConfig *driver_config,
                                  unsigned long TimeIter, unsigned long OuterIter);

  /*!
   * \brief Sets the volume output filename
   * \param[in] filename - the new filename
   */
  inline void SetVolumeFilename(string filename) {volumeFilename = filename;}

  /*!
   * \brief Sets the surface output filename
   * \param[in] filename - the new filename
   */
  inline void SetSurfaceFilename(string filename) {surfaceFilename = filename;}

  /*!
   * \brief Returns the current volume filename
   * \return - The current volume filename
   */
  inline string GetVolumeFilename() {return volumeFilename;}

  /*!
   * \brief Returns the current surface filename
   * \return - The current surface filename
   */
  inline string GetSurfaceFilename() {return surfaceFilename;}

  /*!
   * \brief Sets the restart filename
   * \param[in] filename - the new filename
   */
  inline void SetRestartFilename(string filename) {restartFilename = filename;}

  /*!
   * \brief Returns the current restart filename
   * \return - The current restart filename
   */
  inline string GetRestartFilename() {return restartFilename;}

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
  su2double GetHistoryFieldValue(const string& name) const {
    auto it = historyOutput_Map.find(name);
    if (it != historyOutput_Map.end()) return it->second.value;
    SU2_MPI::Error("Cannot find output field with name " + name, CURRENT_FUNCTION);
    return 0;
  }

 /*!
  * \brief Get the value of particular surface history output field
  * \param[in] field - Name of the field
  * \param[in] iMarker - Index of the surface marker
  * \return Value of the field
  */
  su2double GetHistoryFieldValuePerSurface(const string& name, unsigned short iMarker) const {
    auto it = historyOutputPerSurface_Map.find(name);
    if (it != historyOutputPerSurface_Map.end()) return it->second[iMarker].value;
    SU2_MPI::Error("Cannot find output field with name " + name, CURRENT_FUNCTION);
    return 0;
  }

  /*!
   * \brief Get a vector with all output fields in a particular group
   * \param groupname - Name of the history group
   * \return Vector containing all output fields of a group
   */
  vector<HistoryOutputField> GetHistoryGroup(string groupname){
    vector<HistoryOutputField> HistoryGroup;
    for (unsigned short iField = 0; iField < historyOutput_Map.size(); iField++){
      if (historyOutput_Map.at(historyOutput_List[iField]).outputGroup == groupname){
        HistoryGroup.push_back((historyOutput_Map[historyOutput_List[iField]]));
      }
    }
    return HistoryGroup;
  }

  /*!
   * \brief Get the list of all output fields
   * \return Vector container all output fields
   */
  const vector<string>& GetHistoryOutputList() const {
    return historyOutput_List;
  }

  /*!
   * \brief Get the list of all per-surface fields
   * \return Vector container all output per-surface fields
   */
  const vector<string>& GetHistoryOutputPerSurfaceList() const {
    return historyOutputPerSurface_List;
  }

  /*!
   * \brief Get the map containing all output fields
   * \return Map containing all output fields
   */
  const map<string, HistoryOutputField>& GetHistoryFields() const {
    return historyOutput_Map;
  }

  /*!
   * \brief Get the map containing all output per-surface fields
   * \return Map containing all output per-surface fields
   */
  const map<string, vector<HistoryOutputField>>& GetHistoryPerSurfaceFields() const {
    return historyOutputPerSurface_Map;
  }

  /*!
   * \brief Monitor the convergence of an output field
   * \param[in] config - Definition of the particular problem.
   * \param[in] Iteration - Index of the current iteration.
   * \return Boolean indicating whether the problem is converged.
   */
  bool ConvergenceMonitoring(CConfig *config, unsigned long Iteration);

  /*!
   * \brief Print a summary of the convergence to screen.
   */
  void PrintConvergenceSummary();

  /*!
   * \brief Get convergence of the problem.
   * \return Boolean indicating whether the problem is converged.
   */
  bool GetConvergence() const {return convergence;}

  /*!
   * \brief Set the value of the convergence flag.
   * \param[in] conv - New value of the convergence flag.
   */
  void SetConvergence(const bool conv) {convergence = conv;}

  /*!
   * \brief  Monitor the time convergence of the specified windowed-time-averaged ouput
   * \param[in] config - Definition of the particular problem.
   * \param[in] Iteration - Index of the current iteration.
   * \return Boolean indicating whether the problem is converged.
   */
  bool MonitorTimeConvergence(CConfig *config, unsigned long Iteration);

  /*!
   * \brief Print a list of all history output fields to screen.
   */
  void PrintHistoryFields() const;

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
  bool SetResultFiles(CGeometry *geometry, CConfig *config, CSolver** solver_container,
                       unsigned long iter, bool force_writing = false);

  /*!
   * \brief Get convergence time convergence of the specified windowed-time-averaged ouput of the problem.
   *        Delays solver stop, if Cauchy time convergence criterion is fullfilled, but 2nd order
   *        time marching is active, to ensure that enough restart files are written.
   * \param[in] config - Definition of the particular problem.
   * \return <TRUE> if Solver has converged and has run another iteration.
   */
  bool GetCauchyCorrectedTimeConvergence(const CConfig *config);

  /*!
   * \brief Allocates the appropriate file writer based on the chosen format and writes sorted data to file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] format - The output format.
   * \param[in] fileName - The file name. If empty, the filenames are automatically determined.
   */
  void WriteToFile(CConfig *config, CGeometry *geometry, OUTPUT_TYPE format, string fileName = "");

protected:

  /*----------------------------- Protected member functions ----------------------------*/

  /*!
   * \brief Set the history file header
   * \param[in] config - Definition of the particular problem.
   */
  void SetHistoryFileHeader(const CConfig *config);

  /*!
   * \brief Write the history file output
   * \param[in] config - Definition of the particular problem.
   */
  void SetHistoryFileOutput(const CConfig *config);

  /*!
   * \brief Write the screen header.
   * \param[in] config - Definition of the particular problem.
   */
  void SetScreenHeader(const CConfig *config);

  /*!
   * \brief Write the screen output.
   * \param[in] config - Definition of the particular problem.
   */
  void SetScreenOutput(const CConfig *config);

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
                               HistoryFieldType field_type = HistoryFieldType::DEFAULT) {
    historyOutput_Map[name] = HistoryOutputField(field_name, format, groupname, field_type, description);
    historyOutput_List.push_back(name);
  }

  /*!
   * \brief Set the value of a history output field
   * \param[in] name - Name of the field.
   * \param[in] value - The new value of this field.
   */
  inline void SetHistoryOutputValue(string name, su2double value){
    auto it = historyOutput_Map.find(name);
    if (it != historyOutput_Map.end()){
      it->second.value = value;
    } else {
      SU2_MPI::Error("Cannot find output field with name " + name, CURRENT_FUNCTION);
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
                                         string groupname, const vector<string>& marker_names,
                                         HistoryFieldType field_type = HistoryFieldType::DEFAULT) {
    if (!marker_names.empty()) {
      historyOutputPerSurface_List.push_back(name);
      vector<HistoryOutputField> fields;
      fields.reserve(marker_names.size());
      for (const auto& marker : marker_names) {
        fields.push_back(HistoryOutputField(field_name+"("+marker+")", format, groupname, field_type, ""));
      }
      historyOutputPerSurface_Map[name] = std::move(fields);
    }
  }

  /*!
   * \brief Set the value of a history output field for a specific surface marker
   * \param[in] name - Name of the field.
   * \param[in] value - The new value of this field.
   * \param[in] iMarker - The index of the marker.
   */
  inline void SetHistoryOutputPerSurfaceValue(string name, su2double value, unsigned short iMarker) {
    auto it = historyOutputPerSurface_Map.find(name);
    if (it != historyOutputPerSurface_Map.end()) {
      it->second[iMarker].value = value;
    } else {
      SU2_MPI::Error("Cannot find output field with name " + name, CURRENT_FUNCTION);
    }
  }

  /*!
   * \brief Returns a pointer to the value of an history output.
   * \note For per-surface outputs the marker index is specified as "name[index]".
   */
  inline const su2double* GetPtrToHistoryOutput(const string& name) const {
    /*--- Decide if it should be per surface. ---*/
    const auto pos = name.find('[');
    const su2double* ptr = nullptr;
    if (pos == std::string::npos) {
      const auto it = historyOutput_Map.find(name);
      if (it != historyOutput_Map.end()) {
        ptr = &(it->second.value);
      }
    } else {
      const auto idx = std::stoi(std::string(name.begin()+pos+1, name.end()-1));
      const auto it = historyOutputPerSurface_Map.find(std::string(name, 0, pos));
      if (it != historyOutputPerSurface_Map.end()) {
        ptr = &(it->second[idx].value);
      }
    }
    return ptr;
  }

  /*!
   * \brief Setup a custom history output object for a given expression.
   * \param[in] expression - Some user-defined math with the history field names as variables.
   * \param[out] output - Custom output ready to evaluate.
   */
  void SetupCustomHistoryOutput(const string& expression, CustomHistoryOutput& output) const;

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
   * \param[in] iPoint - The point location in the field.
   */
  su2double GetVolumeOutputValue(const string& name, unsigned long iPoint);

  /*!
   * \brief Set the value of a volume output field
   * \param[in] name - Name of the field.
   * \param[in] iPoint - The point location in the field.
   * \param[in] value - The new value of this field.
   */
  void SetVolumeOutputValue(const string& name, unsigned long iPoint, su2double value);

  /*!
   * \brief Set the value of a volume output field
   * \param[in] name - Name of the field.
   * \param[in] iPoint - The point location in the field.
   * \param[in] value - The new value of this field.
   */
  void SetAvgVolumeOutputValue(const string& name, unsigned long iPoint, su2double value);

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
  void PostprocessHistoryData(CConfig *config);

  /*!
   * \brief Postprocess_HistoryFields
   * \param[in] config - Definition of the particular problem.
   */
  void PostprocessHistoryFields(CConfig *config);

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
   * \brief Write screen and history output.
   * \param[in] config - Definition of the particular problem.
   */
  void OutputScreenAndHistory(CConfig *config);

  /*!
   * \brief Set the history fields common for all solvers.
   */
  void SetCommonHistoryFields();

  /*!
   * \brief Parses user-defined outputs.
   */
  void SetCustomOutputs(const CConfig *config);

  /*!
   * \brief Evaluates function-type custom outputs.
   * Derived classes can use this to compute simple expressions of other outputs if they
   * do not implement surface averages. This should be called just before evaluating the
   * custom objective function.
   */
  void ComputeSimpleCustomOutputs(const CConfig *config);

  /*!
   * \brief Load values of the history fields common for all solvers.
   * \param[in] config - Definition of the particular problem.
   */
  void LoadCommonHistoryData(const CConfig *config);

  /*!
   * \brief Allocates the data sorters if necessary.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void AllocateDataSorters(CConfig *config, CGeometry *geometry);

  /*!
   * \brief Computes the custom and combo objectives.
   * \note To be called after all other history outputs are set.
   * \param[in] idxSol - Index of the main solver.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - The container holding all solution data.
   */
  void SetCustomAndComboObjectives(int idxSol, const CConfig *config, CSolver **solver);

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
  virtual bool WriteHistoryFileOutput(const CConfig *config);

  /*!
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  virtual bool WriteScreenHeader(const CConfig *config);

  /*!
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  virtual bool WriteScreenOutput(const CConfig *config);

  /*!
   * \brief Determines if the the volume output should be written.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Iter - Current iteration index.
   * \param[in] force_writing - boolean that forces writing of volume output
   * \param[in] iFile - index to the file that we need to consider for volume output
   */
  virtual bool WriteVolumeOutput(CConfig *config, unsigned long Iter, bool force_writing, unsigned short iFile);

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
  inline virtual bool SetInitResiduals(const CConfig *config) {return false;}

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
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The container holding all solution data.
   */
  inline virtual void LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {}

  /*!
   * \brief Load the multizone history output field values
   * \param[in] output - Container holding the output instances per zone.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void LoadMultizoneHistoryData(const COutput* const* output, const CConfig* const* config) {}

  /*!
   * \brief Set the available history output fields
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetHistoryOutputFields(CConfig *config) {}

  /*!
   * \brief Set the available multizone history output fields
   * \param[in] output - Container holding the output instances per zone.
   */
  inline virtual void SetMultizoneHistoryOutputFields(const COutput* const* output, const CConfig* const* config) {}

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
  inline virtual void SetAdditionalScreenOutput(const CConfig *config){}

};
