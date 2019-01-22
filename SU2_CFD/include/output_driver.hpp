/*!
 * \file output_driver.hpp
 * \brief Headers of the main subroutines for screen and history output in multizone problems.
 * \author R. Sanchez, T. Albring
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

#include "../../Common/include/mpi_structure.hpp"

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

#include "output_structure.hpp"
#include "../../Common/include/config_structure.hpp"

using namespace std;

/*!
 * \class CDriverOutput
 * \brief Class for writing the multizone output.
 * \author R. Sanchez, T. Albring.
 */
class CDriverOutput {

protected:

  int rank, 	  /*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */

  unsigned short nZone;

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

  std::map<string, HistoryOutputField >         HistoryOutput_Map;    /*!< \brief Associative map to access data stored in the history output fields by a string identifier. */
  std::vector<string>                           HistoryOutput_List;   /*!< \brief Vector that contains the keys of the HistoryOutput_Map in the order of their insertion. */

  std::vector<string> RequestedHistoryFields;
  unsigned short nRequestedHistoryFields;
  std::vector<string> RequestedScreenFields;
  unsigned short nRequestedScreenFields;
  char char_histfile[200];

  ofstream HistFile;

  PrintingToolbox::CTablePrinter* OuterConvergenceTable;

  std::string MultiZoneHeaderString;

  std::map<string, su2double> Init_Residuals;

  map<string, Signal_Processing::RunningAverage> RunningAverages;

public:

  /*!
   * \brief Constructor of the class.
   */
  CDriverOutput(CConfig **config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CDriverOutput(void);

  /*!
   * \brief Write the header of the screen and history file.
   * \param[in] config - Definition of the particular problem.
   */

  void SetHeader(COutput **output, CSolver *****solver, CConfig **config);

  /*!
   * \brief Write the history file and the convergence on the screen for multizone computations.
   * \param[in] output - Output information of the particular multizone problem.
   * \param[in] config - Definition of the particular multizone problem.
   */
  void SetBody(COutput **output, CSolver *****solver, CConfig *driver_config, CConfig **config);

  /*!
   * \brief Load the history file fields for multizone computations.
   * \param[in] output - Output information of the particular multizone problem.
   * \param[in] config - Definition of the particular multizone problem.
   */
  void LoadHistoryData(COutput **output, CSolver *****solver, CConfig *config);

  /*!
   * \brief Write the history output fields for multizone computations.
   * \param[in] output - Output information of the particular multizone problem.
   * \param[in] config - Definition of the particular multizone problem.
   */
  void SetHistoryOutputFields(COutput **output, CSolver *****solver, CConfig **config);

  inline void AddHistoryOutput(string name, string field_name, unsigned short format, string groupname , unsigned short field_type = TYPE_DEFAULT){
    HistoryOutput_Map[name] = HistoryOutputField(field_name, format, groupname, field_type);
    HistoryOutput_List.push_back(name);
  }

  inline void SetHistoryOutputValue(string name, su2double value){
    if (HistoryOutput_Map.count(name) > 0){
      HistoryOutput_Map[name].Value = value;
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
    }
  }

  void SetScreen_Header(CConfig *driver_config, CConfig **config);

  void SetScreen_Output(COutput **output, CConfig *driver_config, CConfig **config);

};
